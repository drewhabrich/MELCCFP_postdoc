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
# Species with dispersal_km > `sparsify_dispersal_km` (search radius ~ the
# whole study area, so nothing gets pruned) are further sparsified to each
# PA's `sparsify_knn` nearest partners plus the Euclidean-MST edges, so the
# downstream cost-distance stage stays feasible and the graph stays connected.
# (Both settings live in setup_script.R; the pair search itself is
# find_candidate_pairs() in pipeline_functions.R.)
#
# Saves per species to Intermediate/01_PA_candidates/:
#   {species}_PApairs_dkernel_cutoff.csv   (from/to site_id + euclidean_dist/kernel_prob)
#   {species}_PApairs_edge_lines.gpkg      (same pairs as edge-to-edge line geometry, for GIS spot-checking)
# and once to Output/01_PA_patches/: pa_patches.gpkg, pa_ctroidnodes.shp
#
## Author: Andrew Habrich
##
## Notes ------------------------- ##
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))

## 1. Load the protected area model shapefile ##################################
pa_mod <- vect(here(pa_dir, "pa_model_5km2.shp"),
               crs = target_crs)

if (!dir.exists(pa_patches_dir)) dir.create(pa_patches_dir)
if (!dir.exists(pa_pairs_dir)) dir.create(pa_pairs_dir)

## Extract the centroid of each PA polygon as the nodes for the minimum spanning tree
nodes <- centroids(pa_mod, inside = T) #inside ensures the node is inside the polygon
nodes$site_id <- seq_len(nrow(nodes))
## Save the centroid nodes to a vector for later (this is mostly for plotting utility)
writeVector(nodes, here(pa_patches_dir, "pa_ctroidnodes.shp"), overwrite = TRUE)

## Also save the PA polygons themselves (not just centroids), with the
## same site_id -- scripts 03/05 need the actual patch shapes for
## edge-to-edge (patch boundary to patch boundary) cost-distance and
## corridor tracing, rather than centroid-to-centroid.
pa_mod$site_id <- seq_len(nrow(pa_mod))
writeVector(pa_mod, here(pa_patches_dir, "pa_patches.gpkg"), overwrite = TRUE)

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
specieslist <- readxl::read_excel(here(tabular_dir, "specieslist.xlsx"))

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to quickly test/compare just those species instead.
species_subset <- c("MAAM", "ASFL", "DOOR", "DRPI", "AQCH")  # e.g. NULL for all species
if (!is.null(species_subset)) {
  specieslist <- specieslist |> filter(species %in% species_subset)
}

pa_sf <- st_as_sf(pa_mod)

kernel_params <- list()

## Persistent cluster for the exact-distance step. sf geometries (unlike
## terra SpatVector/SpatRaster) serialize cleanly to PSOCK workers, and
## pa_sf doesn't change across species, so it's exported once here rather
## than per-species. st_is_within_distance() itself stays on the master --
## it's already fast (spatial-index based), so only the exact
## st_distance() step that follows it needs to be spread across cores.
n_cores <- max(1, parallel::detectCores() - 10)  # tune down if you want to keep using your machine for other things
message("Setting up a ", n_cores, "-worker cluster for parallel distance calculation.")
cl <- parallel::makeCluster(n_cores)
parallel::clusterEvalQ(cl, { library(sf) })
parallel::clusterExport(cl, varlist = "pa_sf")

tryCatch({

for (i in seq_len(nrow(specieslist))) {
  sp <- specieslist$species[i]
  #d_ref is the dispersal_km of the species, derived from the median literature value
  #p_ref describes what proportion of dispersers the d_ref refers to (median = 50%)
  # p_threshold (kernel_p_threshold, setup_script.R) is the probability
  # cutoff we are looking to describe as potential movement; the search
  # radius is the distance at which the kernel probability equals it
  kernel <- species_kernel(specieslist$dispersal_km[i])
  alpha <- kernel$alpha
  p_threshold <- kernel$p_threshold
  search_radius <- kernel$search_radius

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

  candidates <- find_candidate_pairs(pa_sf, pa_mod$site_id, specieslist$dispersal_km[i],
                                     kernel, cl, n_cores, sp)

  if (is.null(candidates)) {
    warning(sp, ": no candidate pairs at all, even under the kernel's ",
            "long-distance tail (p >= ", p_threshold, ") -- every node is ",
            "isolated for this species at this site spacing. Consider ",
            "lowering kernel_p_threshold if that's ecologically defensible.")
    next
  }

  out_path <- file.path(pa_pairs_dir, paste0(sp, "_PApairs_dkernel_cutoff.csv"))
  write_csv(candidates, out_path)

  message(sp, ": alpha = ", signif(alpha, 4), ", search radius = ",
          round(search_radius), " map units -> ", nrow(candidates),
          " candidate pair(s) -> ", out_path)

  ## Edge-to-edge line geometry for each surviving candidate pair (nearest
  ## boundary point to nearest boundary point), tagged with from/to site_id and the computed
  ## distance -- lets the reported euclidean_dist be spot-checked by
  ## manually measuring these lines in GIS.
  from_polys <- pa_sf[match(candidates$from, pa_mod$site_id), ]
  to_polys   <- pa_sf[match(candidates$to,   pa_mod$site_id), ]

  pair_lines <- st_nearest_points(from_polys, to_polys, pairwise = TRUE) %>%
    st_as_sf() %>%
    mutate(from = candidates$from,
           to = candidates$to,
           euclidean_dist = candidates$euclidean_dist,
           euclidean_kernel_prob = candidates$euclidean_kernel_prob)

  lines_out_path <- file.path(pa_pairs_dir, paste0(sp, "_PApairs_edge_lines.gpkg"))
  st_write(pair_lines, lines_out_path, delete_dsn = TRUE, quiet = TRUE)
  message(sp, ": saved edge-to-edge candidate-pair lines to ", lines_out_path)
}

}, finally = {
  parallel::stopCluster(cl)
  message("Cluster stopped.")
})

kernel_params_tbl <- bind_rows(kernel_params)
write_csv(kernel_params_tbl, file.path(table_dir, "01_dispersal_kernel_params.csv"))
