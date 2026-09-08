## ------------------------------- ##
## Script name: 01-PA_candidates
##
## Purpose of script: Prefilter PA pairs by edge-to-edge (patch boundary to
# patch boundary) distance, per species' dispersal capacity, using a
# negative exponential dispersal kernel to estimate the search radius
# capturing 95% of *possible* movement. Cheap to compute (via a spatial
# index, not a full pairwise distance matrix), and avoids running
# expensive cost-distance on pairs that could never be within range --
# cost-distance is always >= edge-to-edge distance when resistance >= 1.
#
# Saves one .csv per species to
# data/intermediate/{species}_PApairs_dkernel_cutoff.csv
#
## Author: Andrew Habrich
##
## Notes ------------------------- ##
source(file.path("scripts", "setup_script.R"))

## 1. Load the protected area model shapefile ##################################
pa_mod <- vect(here(interm_dir, "pa_model_5km2.shp"),
               crs = target_crs)

# remove duplicated polygons (if any) based on geometry (may occur during merging)
#pa_clean <- tidyterra::distinct(pa_mod, geometry, .keep_all = T)
# save it to file for later use
#writeVector(pa_clean, here(interm_dir, "pa_model_clean.shp"), overwrite = T)

## Extract the centroid of each PA polygon as the nodes for the minimum spanning tree
nodes <- centroids(pa_mod, inside = T) #inside ensures the node is inside the polygon
nodes$site_id <- seq_len(nrow(nodes))

## Also persist the PA polygons themselves (not just centroids), with the
## same site_id -- scripts 03/05 need the actual patch shapes for
## edge-to-edge (patch boundary to patch boundary) cost-distance and
## corridor tracing, rather than centroid-to-centroid.
pa_mod$site_id <- seq_len(nrow(pa_mod))
writeVector(pa_mod, here(output_dir, "pa_patches.gpkg"), overwrite = TRUE)

## Quick visual check of the centroids
#plot(pa_mod)
#plot(nodes, add = T, col = "red")
#hist(nodes$area)   # distribution of patch areas (this should be the same as above)

## Save the centroid nodes to a vector for later
writeVector(nodes, here(output_dir, "pa_ctroidnodes.shp"), overwrite = TRUE)

## 2. Eligible PA pairs using dispersal kernel cutoff ##########################
## Uses a negative exponential dispersal kernel to model probability of
## movement, and finds each species' candidate pairs via
## st_is_within_distance() (spatial-index based) rather than a full N x N
## distance matrix: most PA pairs are far beyond any species' search
## radius and would just be discarded, so pruning them via the index
## before computing any exact edge-to-edge distance is far cheaper than
## computing every pair exhaustively. Search radius differs per species
## (a few tens of metres up to hundreds of km), so this runs once per
## species using that species' own radius.
specieslist <- readxl::read_excel(here(data_dir, "specieslist.xlsx"))

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to quickly test/compare just those species instead.
species_subset <- c("MAAM", "DOOR")  # e.g. NULL for all species
if (!is.null(species_subset)) {
  specieslist <- specieslist |> filter(species %in% species_subset)
}

pa_sf <- st_as_sf(pa_mod)

kernel_params <- list()

## Persistent cluster for the exact-distance step below. sf geometries
## (unlike terra SpatVector/SpatRaster) serialize cleanly to PSOCK workers,
## and pa_sf doesn't change across species, so it's exported once here
## rather than per-species. st_is_within_distance() itself stays on the
## master -- it's already fast (spatial-index based), so only the exact
## st_distance() step that follows it needs to be spread across cores.
n_cores <- max(1, parallel::detectCores() - 10)  # tune down if you want to keep using your machine for other things
message("Setting up a ", n_cores, "-worker cluster for parallel distance calculation.")
cl <- parallel::makeCluster(n_cores)
parallel::clusterEvalQ(cl, { library(sf) })
parallel::clusterExport(cl, varlist = "pa_sf")

compute_chunk_dist <- function(rows) {
  units::drop_units(
    st_distance(pa_sf[rows$from_idx, ], pa_sf[rows$to_idx, ], by_element = TRUE)
  )
}

tryCatch({

for (i in seq_len(nrow(specieslist))) {
  sp <- specieslist$species[i]
  #d_ref is the dispersal_km of the species, derived from the median literature value
  #p_ref describes what proportion of dispersers the d_ref refers to (median = 50%)
  alpha <- dispersal_kernel_alpha(
    d_ref = specieslist$dispersal_km[i]*1000,
    p_ref = 0.5
  )
  # p_threshold describes the probability cutoff we are looking to describe as potential movement
  p_threshold <- 0.05

  # distance at which the kernel probability equals p_threshold -- the
  # search radius for this species' prefilter
  search_radius <- -log(p_threshold) / alpha

  message("Processing species ", i, "/", nrow(specieslist), ": ", sp,
          " (search radius = ", round(search_radius / 1000, 1), " km)")

  kernel_params[[sp]] <- tibble(
    species             = sp,
    kernel_ref_dist_m = specieslist$dispersal_km[i]*1000,
    kernel_ref_prob     = 0.5,
    alpha               = alpha,
    kernel_p_threshold  = p_threshold,
    search_radius       = search_radius
  )

  ## Sparse neighbor list within this species' search radius (spatial-
  ## index based -- fast), then exact edge-to-edge distance for just those
  ## pairs, not the full N x N matrix. NOTE: this prefilter only helps as
  ## much as the search radius actually excludes -- a species with a huge
  ## radius (e.g. LYCA, ~100 km dispersal) can still end up with a
  ## candidate set approaching the full N^2/2 pairs, in which case the
  ## exact-distance step below is genuinely a lot of work, not a hang.
  nbrs <- st_is_within_distance(pa_sf, dist = search_radius, sparse = TRUE)
  pairs <- bind_rows(lapply(seq_along(nbrs), function(idx) {
    js <- nbrs[[idx]][nbrs[[idx]] > idx]  # from < to only -- drop self-matches and duplicate pairs
    if (length(js) == 0) return(NULL)
    tibble(from_idx = idx, to_idx = js)
  }))

  if (nrow(pairs) == 0) {
    warning(sp, ": no candidate pairs at all, even under the kernel's ",
            "long-distance tail (p >= ", p_threshold, ") -- every node is ",
            "isolated for this species at this site spacing. Consider ",
            "lowering kernel_p_threshold if that's ecologically defensible.")
    next
  }

  n_pairs <- nrow(pairs)
  message("  ", sp, ": ", n_pairs, " candidate pair(s) after spatial-index ",
          "prefilter -- computing exact distances...")

  ## Split into chunks and dispatch them across the cluster (load-balanced),
  ## in batches so progress can still be reported as batches complete --
  ## more chunks than workers (n_cores * 4) keeps workers busy even when
  ## chunks finish at uneven speeds.
  n_chunks <- min(n_pairs, max(20, n_cores * 4))
  chunk_size <- ceiling(n_pairs / n_chunks)
  chunk_ids <- ceiling(seq_len(n_pairs) / chunk_size)
  chunks <- split(pairs, chunk_ids)

  exact_dist <- numeric(n_pairs)
  dispatch_batches <- split(seq_along(chunks), ceiling(seq_along(chunks) / (n_cores * 2)))
  n_done <- 0
  for (b in dispatch_batches) {
    batch_out <- parallel::parLapplyLB(cl, chunks[b], compute_chunk_dist)
    exact_dist[unlist(lapply(chunks[b], function(x) as.integer(rownames(x))))] <- unlist(batch_out)
    n_done <- n_done + sum(vapply(chunks[b], nrow, integer(1)))
    message("  ", sp, ": ", n_done, "/", n_pairs, " pair(s) processed")
  }

  candidates <- tibble(
    from = pa_mod$site_id[pairs$from_idx],
    to   = pa_mod$site_id[pairs$to_idx],
    euclidean_dist = exact_dist
  ) %>%
    filter(euclidean_dist <= search_radius) %>%  # defensive: guards any index/exact boundary mismatch
    mutate(euclidean_kernel_prob = dispersal_probability(euclidean_dist, alpha))

  out_path <- file.path(interm_dir, paste0(sp, "_PApairs_dkernel_cutoff.csv"))
  write_csv(candidates, out_path)

  message(sp, ": alpha = ", signif(alpha, 4), ", search radius = ",
          round(search_radius), " map units -> ", nrow(candidates),
          " candidate pair(s) -> ", out_path)
}

}, finally = {
  parallel::stopCluster(cl)
  message("Cluster stopped.")
})

kernel_params_tbl <- bind_rows(kernel_params)
write_csv(kernel_params_tbl, file.path(interm_dir, "dispersal_kernel_params.csv"))

