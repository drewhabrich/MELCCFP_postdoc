## ------------------------------- ##
## Script name: pipeline_functions
##
## Purpose of script: The per-species core of every pipeline stage
##   (setup_species_from_xlsx.R, 01-07) as plain functions, shared by the
##   stage scripts AND the sensitivity analysis (08/09) -- so a sensitivity
##   run executes exactly the production code, with only its input
##   parameters changed. Sourced at the end of setup_script.R; never source
##   it on its own (it relies on the packages, template_raster, kernel
##   helpers and shared parameters defined there).
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
# PARALLEL WORKERS: several functions below run ON PSOCK workers
# (compute_chunk_dist, crop_extent_for_origin, cost_distance_for_origin,
# attempt_edge_lcp, trace_edge_lcp, lcp_worker, edge_density_worker). They
# are deliberately defined at top level (global environment): they
# reference per-species objects (pa_sf, candidate_pairs, resistance_path,
# mst_edges, edge_coords_list, ...) as FREE variables, which clusterExport()
# places in each worker's global environment. Don't nest them inside
# another function -- their enclosing environment would then be serialized
# with them and the free-variable lookup would stop reaching the exported
# copies. The orchestrating functions (find_candidate_pairs(),
# compute_cost_distances(), trace_lcp_corridors(), corridor_density())
# export their per-species objects with clusterExport(..., envir =
# environment()) for the same reason.
# ============================================================

## 1. Species rank tables (setup_species_from_xlsx.R) ##########################
## Albert et al. 2017 doubling scheme -- see setup_species_from_xlsx.R's
## header for the full description of the lookup tables these read.

non_forest_categories <- c(wetland_open = "wetland_open", wetland_forested = "wetland_forested",
                           grass = "grass", shrub = "shrub", agriculture = "agriculture")
forest_type_cols <- c(conif = "conif", decid = "decid", mixed = "mixed")

# Water (LULC_code 4) isn't in non_forest_categories -- it uses a shared
# fixed_rank (landcover_class_lookup.csv) like Urban/roads, adjusted per
# species by dispersal_type: species that cross water cheaply (flying,
# aquatic) get a lower rank; fossorial species get a higher one.
dispersal_type_water_adjustment <- c(flying = -1, aquatic = -1, fossorial = 1, terrestrial = 0)

# Roads/railway (LULC_codes 100/200/300/400/500/600) also use a shared
# fixed_rank -- flying species can cross these cheaply (fly over), so they
# get a rank discount; other dispersal types are unaffected.
road_rail_codes <- c(100, 200, 300, 400, 500, 600)  # Major/Minor/Forestry A-C roads, Railway
dispersal_type_road_adjustment <- c(flying = -1, aquatic = 0, fossorial = 0, terrestrial = 0)

## One species' landcover / forest-type / forest-age-penalty rank tables.
## `sp_row` is that species' single row of the xlsx "top10" sheet.
build_rank_tables <- function(sp_row, landcover_lookup, forest_type_lookup, forest_age_lookup,
                              rank_scale) {
  sp <- sp_row$species[1]

  # ---- forest type rank table (applied only within forest LULC cells) ----
  # (computed before the land-cover table below, so its average can seed
  # the Forest fallback rank used there)
  sp_rank_by_type <- tibble(
    type_category    = names(forest_type_cols),
    suitability = unlist(sp_row[1, forest_type_cols], use.names = FALSE)
  ) %>%
    left_join(rank_scale %>% dplyr::select(suitability, rank), by = "suitability")

  # Forest (LULC_code 2) fallback -- used by build_resistance() only for
  # forest cells where forest_type/forest_age data is missing. This
  # species' average rank across its own conif/decid/mixed suitability,
  # i.e. assume typical forest suitability when the exact type/age can't
  # be determined.
  forest_fallback_rank <- round(mean(sp_rank_by_type$rank, na.rm = TRUE))

  forest_type_rank_tbl <- forest_type_lookup %>%
    left_join(sp_rank_by_type %>% dplyr::select(type_category, rank), by = "type_category") %>%
    dplyr::select(forest_type_class, rank)

  if (any(is.na(forest_type_rank_tbl$rank))) {
    warning(sp, ": some forest_type classes have no resolvable rank -- check ",
            "forest_type_class_lookup.csv's type_category values match conif/decid/mixed.")
  }

  # ---- non-forest land-cover rank table ----
  sp_rank_by_category <- tibble(
    habitat_category = names(non_forest_categories),
    suitability = unlist(sp_row[1, non_forest_categories], use.names = FALSE)
  ) %>%
    left_join(rank_scale %>% dplyr::select(suitability, rank), by = "suitability")

  if (any(is.na(sp_rank_by_category$rank))) {
    warning(sp, ": suitability tier not found in suitability_rank_scale.csv for some ",
            "habitat categories -- check spelling matches exactly (e.g. 'Highly suitable').")
  }

  water_adjustment <- dispersal_type_water_adjustment[sp_row$dispersal_type[1]]
  if (is.na(water_adjustment)) {
    warning(sp, ": dispersal_type '", sp_row$dispersal_type[1],
            "' not recognized -- no water rank adjustment applied.")
    water_adjustment <- 0
  }

  road_adjustment <- dispersal_type_road_adjustment[sp_row$dispersal_type[1]]
  if (is.na(road_adjustment)) {
    warning(sp, ": dispersal_type '", sp_row$dispersal_type[1],
            "' not recognized -- no road/rail rank adjustment applied.")
    road_adjustment <- 0
  }

  landcover_rank_tbl <- landcover_lookup %>%
    filter(LULC_code != 2) %>%
    left_join(sp_rank_by_category %>% dplyr::select(habitat_category, rank), by = "habitat_category") %>%
    mutate(rank = coalesce(fixed_rank, rank)) %>%
    mutate(rank = if_else(LULC_code == 4, pmax(rank + water_adjustment, 0), rank)) %>%
    mutate(rank = if_else(LULC_code %in% road_rail_codes, pmax(rank + road_adjustment, 0), rank)) %>%
    dplyr::select(cover_class = LULC_code, rank) %>%
    bind_rows(tibble(cover_class = 2, rank = forest_fallback_rank))

  if (any(is.na(landcover_rank_tbl$rank))) {
    warning(sp, ": ", sum(is.na(landcover_rank_tbl$rank)),
            " non-forest land-cover class(es) have no resolvable rank -- check ",
            "landcover_class_lookup.csv's habitat_category/fixed_rank values.")
  }

  # ---- forest age penalty table (applied only within forest LULC cells) ----
  # sp_row's forest_age column is a comma-separated list of this species'
  # preferred age classes, e.g. "mature, forêt ancienne". Penalty per
  # raster class is the minimum ordinal distance (forest_age_class_lookup.csv)
  # across all of a species' preferred categories -- best-case match wins.
  preferred_ages <- sp_row$forest_age[1] %>%
    str_split(",") %>%
    unlist() %>%
    str_trim()

  age_penalty_tbl <- forest_age_lookup %>%
    filter(age_category %in% preferred_ages) %>%
    group_by(forest_age_class) %>%
    summarize(penalty = min(distance), .groups = "drop")

  if (nrow(age_penalty_tbl) != 5) {
    warning(sp, ": forest age penalty table doesn't cover all 5 classes -- check ",
            "specieslist.xlsx's forest_age value(s) match forest_age_class_lookup.csv's ",
            "age_category values exactly (', ' is required as list separator).")
  }

  list(landcover = landcover_rank_tbl, forest_type = forest_type_rank_tbl,
       age_penalty = age_penalty_tbl)
}

## 2. Resistance raster (02) ###################################################

## The three land-cover-derived input rasters, aligned to the shared 30 m
## reference grid (template_raster) once, up front, if needed (nearest
## neighbour -- all categorical), plus the forest mask (LULC class 2).
load_resistance_inputs <- function() {
  align_categorical <- function(r, target, label) {
    if (!compareGeom(target, r, stopOnError = FALSE)) {
      message(label, " grid doesn't match the reference grid -- resampling (nearest neighbor).")
      r <- resample(r, target, method = "near")
    }
    r
  }
  landcover   <- align_categorical(rast(here(output_dir, "lulc_road_model.tif")), template_raster, "Land cover")
  forest_type <- align_categorical(rast(here(output_dir, "cover_type_model.tif")), template_raster, "Forest type")
  forest_age  <- align_categorical(rast(here(output_dir, "age_model.tif")), template_raster, "Forest age")
  list(landcover = landcover, forest_type = forest_type, forest_age = forest_age,
       forest_mask = landcover == 2)
}

## ALBERT ET AL. 2017 DOUBLING SCALE: resistance = 2^total_rank (rank 0 ->
## 1, each worse rank step doubles the cost). `rank_tbls` is
## list(landcover, forest_type, age_penalty), as returned by
## build_rank_tables() or read back from its CSVs. `forest_unsuitable_rank`
## is the rank of the "Unsuitable" tier; forest types at or above it get
## their age penalty capped at unsuitable_age_penalty_cap (NULL = no cap).
## Writes to `out_path` and returns this species' worst-case rank (from its
## tables, not hardcoded).
build_resistance <- function(rank_tbls, landcover, forest_type, forest_age, forest_mask,
                             out_path, forest_unsuitable_rank = NULL) {
  lulc_rank <- classify(landcover, as.matrix(rank_tbls$landcover))
  forest_type_rank <- classify(forest_type, as.matrix(rank_tbls$forest_type))
  age_penalty <- classify(forest_age, as.matrix(rank_tbls$age_penalty))

  # Forest types the species rates Unsuitable (type rank >= the Unsuitable
  # tier's rank) get at most unsuitable_age_penalty_cap
  # of age penalty -- unsuitable forest tops out level with urban. The cap is
  # per forest-type cell, so suitable types keep their full age penalty.
  if (!is.null(forest_unsuitable_rank)) {
    age_penalty <- ifel(forest_type_rank >= forest_unsuitable_rank,
                        min(age_penalty, unsuitable_age_penalty_cap), age_penalty)
  }

  forest_total_rank <- forest_type_rank + age_penalty
  # Forest cells where forest_type/forest_age data is missing (forest_total_rank
  # is NA there) fall back to lulc_rank's Forest row -- this species' average
  # conif/decid/mixed rank (see build_rank_tables()) -- instead of
  # propagating NA into the resistance raster.
  total_rank <- ifel(forest_mask,
                      ifel(is.na(forest_total_rank), lulc_rank, forest_total_rank),
                      lulc_rank)

  # worst-case forest rank: a type at/above the Unsuitable rank only takes the capped penalty
  type_ranks <- rank_tbls$forest_type$rank
  max_penalty <- max(rank_tbls$age_penalty$penalty, na.rm = TRUE)
  forest_max <- if (is.null(forest_unsuitable_rank)) {
    max(type_ranks, na.rm = TRUE) + max_penalty
  } else {
    max(type_ranks + ifelse(type_ranks >= forest_unsuitable_rank,
                            min(max_penalty, unsuitable_age_penalty_cap), max_penalty), na.rm = TRUE)
  }
  max_rank <- max(max(rank_tbls$landcover$rank, na.rm = TRUE), forest_max)

  # true doubling scale: rank 0 -> 1, each rank step doubles the cost
  species_resistance <- 2^total_rank

  writeRaster(species_resistance, out_path, overwrite = TRUE)
  invisible(max_rank)
}

## 3. Dispersal kernel + candidate PA pairs (01) ###############################

## A species' negative-exponential kernel: dispersal_km is the median
## (p_ref = 0.5) dispersal distance; search_radius is where the kernel
## drops to p_threshold (= dispersal_km * log2(1 / p_threshold)).
species_kernel <- function(dispersal_km, p_threshold = kernel_p_threshold) {
  alpha <- dispersal_kernel_alpha(d_ref = dispersal_km * 1000, p_ref = 0.5)
  list(alpha = alpha, p_threshold = p_threshold, search_radius = -log(p_threshold) / alpha)
}

## Worker: exact edge-to-edge distance for one chunk of pairs. `pa_sf` is a
## free variable -- clusterExport() it to the cluster first.
compute_chunk_dist <- function(rows) {
  units::drop_units(
    st_distance(pa_sf[rows$from_idx, ], pa_sf[rows$to_idx, ], by_element = TRUE)
  )
}

## Candidate PA pairs within `kernel$search_radius`, edge-to-edge: sparse
## neighbour list from st_is_within_distance() (spatial index -- fast),
## then exact st_distance() for just those pairs, spread over `cl` (which
## must already have `pa_sf` exported). Species with dispersal_km >
## sparsify_dispersal_km are sparsified to each PA's sparsify_knn nearest
## partners plus the Euclidean-MST edges (see 01's header). Returns a
## tibble(from, to, euclidean_dist, euclidean_kernel_prob), or NULL if no
## pair is within range.
find_candidate_pairs <- function(pa_sf, site_ids, dispersal_km, kernel, cl, n_workers, sp = "") {
  search_radius <- kernel$search_radius

  ## NOTE: this prefilter only helps as much as the search radius actually
  ## excludes -- a species with a huge radius can still end up with a
  ## candidate set approaching the full N^2/2 pairs, in which case the
  ## exact-distance step below is genuinely a lot of work, not a hang.
  nbrs <- st_is_within_distance(pa_sf, dist = search_radius, sparse = TRUE)
  pairs <- bind_rows(lapply(seq_along(nbrs), function(idx) {
    js <- nbrs[[idx]][nbrs[[idx]] > idx]  # from < to only -- drop self-matches and duplicate pairs
    if (length(js) == 0) return(NULL)
    tibble(from_idx = idx, to_idx = js)
  }))

  if (nrow(pairs) == 0) return(NULL)

  n_pairs <- nrow(pairs)
  message("  ", sp, ": ", n_pairs, " candidate pair(s) after spatial-index ",
          "prefilter -- computing exact distances...")

  ## Split into chunks and dispatch them across the cluster (load-balanced),
  ## in batches so progress can still be reported as batches complete --
  ## more chunks than workers (n_workers * 4) keeps workers busy even when
  ## chunks finish at uneven speeds.
  n_chunks <- min(n_pairs, max(20, n_workers * 4))
  chunk_size <- ceiling(n_pairs / n_chunks)
  chunk_ids <- ceiling(seq_len(n_pairs) / chunk_size)
  ## row_id tracks each pair's true position in `pairs` explicitly --
  ## tibbles (unlike base data.frame) reset rownames() to a fresh 1:n on
  ## every subset, so split()'s chunks can't be re-placed into exact_dist
  ## via rownames() (that silently scrambled/zeroed results before).
  pairs$row_id <- seq_len(n_pairs)
  chunks <- split(pairs, chunk_ids)

  exact_dist <- numeric(n_pairs)
  dispatch_batches <- split(seq_along(chunks), ceiling(seq_along(chunks) / (n_workers * 2)))
  n_done <- 0
  for (b in dispatch_batches) {
    batch_out <- parallel::parLapplyLB(cl, chunks[b], compute_chunk_dist)
    exact_dist[unlist(lapply(chunks[b], function(x) x$row_id))] <- unlist(batch_out)
    n_done <- n_done + sum(vapply(chunks[b], nrow, integer(1)))
    message("  ", sp, ": ", n_done, "/", n_pairs, " pair(s) processed")
  }

  candidates <- tibble(
    from = site_ids[pairs$from_idx],
    to   = site_ids[pairs$to_idx],
    euclidean_dist = exact_dist
  ) %>%
    filter(euclidean_dist <= search_radius) %>%  # defensive: guards any index/exact boundary mismatch
    mutate(euclidean_kernel_prob = dispersal_probability(euclidean_dist, kernel$alpha))

  ## A very large dispersal distance makes the search radius approach or
  ## exceed the study area, so essentially every PA pair survives the
  ## prefilter -- and cost-distance for pairs hundreds of km apart at 30 m
  ## is infeasible (each crop spans ~1e8+ cells). An MST only needs each
  ## node's short links plus enough bridges to stay connected, so keep each
  ## PA's `sparsify_knn` nearest partners (edge-to-edge) plus the Euclidean-
  ## MST edges (which guarantee the pruned graph stays connected).
  if (dispersal_km > sparsify_dispersal_km) {
    n_before <- nrow(candidates)
    candidates <- candidates %>% mutate(pair_id = row_number())

    knn_ids <- bind_rows(
      candidates %>% transmute(node = from, pair_id, euclidean_dist),
      candidates %>% transmute(node = to,   pair_id, euclidean_dist)
    ) %>%
      group_by(node) %>%
      slice_min(euclidean_dist, n = sparsify_knn, with_ties = FALSE) %>%
      pull(pair_id)

    euc_graph <- igraph::graph_from_data_frame(
      candidates %>% dplyr::select(from, to, euclidean_dist, pair_id), directed = FALSE
    )
    euc_mst <- igraph::mst(euc_graph, weights = igraph::E(euc_graph)$euclidean_dist)
    mst_ids <- igraph::E(euc_mst)$pair_id

    candidates <- candidates %>%
      filter(pair_id %in% union(knn_ids, mst_ids)) %>%
      dplyr::select(-pair_id)

    message("  ", sp, ": dispersal_km = ", dispersal_km, " > ", sparsify_dispersal_km,
            " -- sparsified to k-NN (k = ", sparsify_knn, ") + Euclidean-MST edges: ",
            n_before, " -> ", nrow(candidates), " pair(s)")
  }

  candidates
}

## 4. Cost-distance (03) #######################################################

## terra decides how much of a raster to keep resident in memory (vs.
## spilling to disk) by querying the NODE's total physical RAM -- on a
## shared HPC node that can be hundreds of GB, far more than this job's
## actual --mem-per-cpu x --cpus-per-task cgroup allocation, so it can keep
## a huge crop fully resident right up until the cgroup OOM-kills the
## process. This reads the job's real allocation from SLURM's own env vars
## (NA outside a SLURM job, e.g. interactively -- terra's default
## node-RAM-based heuristic is then reasonable again).
slurm_job_mem_gb <- function() {
  mem_per_cpu_mb <- suppressWarnings(as.numeric(Sys.getenv("SLURM_MEM_PER_CPU", NA)))
  cpus_per_task  <- suppressWarnings(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
  if (is.na(mem_per_cpu_mb) || is.na(cpus_per_task)) NA_real_ else (mem_per_cpu_mb / 1024) * cpus_per_task
}

## Bounding-box extent of one origin plus all its candidate destinations,
## padded by a buffer sized as a fraction of the farthest destination.
## Shared by cost_distance_for_origin() (which crops the real raster to it)
## and estimate_origin_crop_cells() (which uses it, before any raster is
## touched, to route unusually wide-spread origins to the low-concurrency
## cluster).
crop_extent_for_origin <- function(origin_poly, dest_polys, dest_pairs, buffer_mult = 1) {
  buffer_dist <- max(max(dest_pairs$euclidean_dist) * cd_buffer_frac, cd_buffer_min_m) * buffer_mult
  ext(vect(rbind(origin_poly, dest_polys))) + buffer_dist
}

## Worker: cost-distance from one origin to all its candidate destinations,
## edge-to-edge (patch boundary to patch boundary): the whole origin
## polygon is marked as a zero-cost source region (terra::costDist()
## already supports multi-cell sources), and each destination's cost is
## the MINIMUM accumulated cost over its own polygon's cells -- i.e. the
## cost to whichever point on that patch's boundary is cheapest to reach,
## not to a single centroid. Both on a crop of the resistance raster sized
## to just this origin's search area. References `candidate_pairs` /
## `resistance_path` / `pa_patches` as free variables (exported by
## setup_costdist_cluster() / compute_cost_distances()). Returns
## list(data = <the pair-cost tibble>, timing = <per-step seconds for this
## origin, for profiling>) -- timing is split out from data so the
## diagnostic columns never leak into {sp}_PA_costs.csv. If a destination
## comes back unreachable within the crop, retries once with a doubled
## buffer before accepting the NA.
cost_distance_for_origin <- function(origin_id, buffer_mult = 1) {
  t0 <- Sys.time()
  origin_poly <- pa_patches[as.character(pa_patches$site_id) == as.character(origin_id), ]

  dest_pairs <- candidate_pairs  |>  filter(from == origin_id)
  dest_ids <- dest_pairs |> pull(to)
  # match(), not %in% filtering, to guarantee dest_polys is in the same
  # order as dest_ids -- extract()'s output order follows the polygon
  # input order, and it needs to line up with dest_ids below.
  dest_polys <- pa_patches[match(as.character(dest_ids), as.character(pa_patches$site_id)), ]

  crop_ext <- crop_extent_for_origin(origin_poly, dest_polys, dest_pairs, buffer_mult)
  t_setup <- Sys.time()

  friction <- crop(rast(resistance_path), crop_ext)
  t_crop <- Sys.time()

  origin_mask <- rasterize(vect(origin_poly), friction)
  friction <- ifel(!is.na(origin_mask), sentinel_value, friction)
  t_mask <- Sys.time()

  cost_surface <- costDist(friction, target = sentinel_value, scale = 1)
  t_costdist <- Sys.time()

  cost_vals <- terra::extract(cost_surface, vect(dest_polys), fun = "min", na.rm = TRUE)[[2]]
  t_extract <- Sys.time()

  terra::tmpFiles(remove = TRUE)
  t_cleanup <- Sys.time()

  timing <- tibble(
    origin_id    = origin_id,
    buffer_mult  = buffer_mult,
    n_dest       = length(dest_ids),
    n_cells      = ncell(friction),
    setup_sec    = as.numeric(difftime(t_setup, t0, units = "secs")),
    crop_sec     = as.numeric(difftime(t_crop, t_setup, units = "secs")),
    mask_sec     = as.numeric(difftime(t_mask, t_crop, units = "secs")),
    costdist_sec = as.numeric(difftime(t_costdist, t_mask, units = "secs")),
    extract_sec  = as.numeric(difftime(t_extract, t_costdist, units = "secs")),
    cleanup_sec  = as.numeric(difftime(t_cleanup, t_extract, units = "secs")),
    total_sec    = as.numeric(difftime(t_cleanup, t0, units = "secs"))
  )

  if (any(is.na(cost_vals)) && buffer_mult == 1) {
    retry <- cost_distance_for_origin(origin_id, buffer_mult = 2)
    retry$timing <- bind_rows(timing, retry$timing)
    return(retry)
  }

  list(
    data = tibble(from = origin_id, to = dest_ids,
                  euclidean_dist = dest_pairs$euclidean_dist,
                  euclidean_kernel_prob = dest_pairs$euclidean_kernel_prob,
                  cost_dist = cost_vals),
    timing = timing
  )
}

## Builds a fresh cost-distance cluster with everything that doesn't change
## across species already exported (needs `pa_patches` in the global
## environment). Each worker gets an explicit share of the job's real
## memory budget (see slurm_job_mem_gb()) instead of letting terra assume
## it has the whole node, and spills to the job's local scratch rather than
## a default tempdir that may be slower or shared. 80% leaves headroom for
## R/OS overhead so terra spills before the cgroup limit is actually hit.
## `memmax_gb` (optional) sets that per-worker cap explicitly instead --
## e.g. on a desktop, where one MAAM-sized crop can commit ~14 GB and a few
## concurrent workers overcommit RAM (08 uses this locally).
setup_costdist_cluster <- function(ncores, memmax_gb = NA_real_) {
  new_cl <- parallel::makeCluster(ncores)
  parallel::clusterEvalQ(new_cl, {
    library(terra); library(sf); library(dplyr); library(tibble)
  })
  job_mem_gb <- slurm_job_mem_gb()
  worker_memmax_gb <- if (!is.na(memmax_gb)) memmax_gb else if (is.na(job_mem_gb)) NA_real_ else (job_mem_gb * 0.8) / ncores
  slurm_tmpdir <- Sys.getenv("SLURM_TMPDIR", "")
  parallel::clusterCall(new_cl, function(memmax_gb, tmp) {
    if (!is.na(memmax_gb)) terra::terraOptions(memmax = memmax_gb)
    if (nzchar(tmp)) terra::terraOptions(tempdir = tmp)
  }, worker_memmax_gb, slurm_tmpdir)
  parallel::clusterExport(new_cl, varlist = c(
    "pa_patches", "sentinel_value", "cd_buffer_frac", "cd_buffer_min_m",
    "crop_extent_for_origin", "cost_distance_for_origin"
  ))
  new_cl
}

## Cheap estimate of one origin's crop size in cells -- same bbox + buffer
## logic as cost_distance_for_origin(), but pure geometry, no raster read.
estimate_origin_crop_cells <- function(origin_id, candidate_pairs, pa_patches, res_m) {
  origin_poly <- pa_patches[as.character(pa_patches$site_id) == as.character(origin_id), ]
  dest_pairs  <- candidate_pairs |> filter(from == origin_id)
  dest_polys  <- pa_patches[match(as.character(dest_pairs$to), as.character(pa_patches$site_id)), ]
  e <- as.vector(crop_extent_for_origin(origin_poly, dest_polys, dest_pairs))
  ((e["xmax"] - e["xmin"]) / res_m) * ((e["ymax"] - e["ymin"]) / res_m)
}

## Cost-distance for every candidate pair of one species. Origins whose
## estimated crop exceeds `big_origin_cells` go to the low-concurrency
## `cl_big` cluster, the rest to `cl` (see 03's header). Returns
## list(costs, timing, cl, cl_big) -- a cluster whose worker died is
## rebuilt, so the caller must keep the returned clusters.
compute_cost_distances <- function(candidate_pairs, resistance_path, cl, cl_big,
                                   n_cores, big_n_cores, batch_size, big_origin_cells, sp = "",
                                   memmax_gb = NA_real_) {
  cd_env <- environment()
  ## `candidate_pairs` / `resistance_path` change every species -- export to both clusters.
  parallel::clusterExport(cl, varlist = c("candidate_pairs", "resistance_path"), envir = cd_env)
  parallel::clusterExport(cl_big, varlist = c("candidate_pairs", "resistance_path"), envir = cd_env)

  origins <- unique(candidate_pairs$from)
  n_origins <- length(origins)

  ## Cheaply estimate each origin's crop size so wide-destination-spread
  ## origins can be routed to the low-concurrency cluster before any worker
  ## actually crops a huge raster.
  res_m <- terra::res(terra::rast(resistance_path))[1]
  origin_cells <- vapply(origins, estimate_origin_crop_cells, numeric(1),
                         candidate_pairs = candidate_pairs, pa_patches = pa_patches, res_m = res_m)
  ## Largest crops first, so the slowest origins start early instead of
  ## trailing at the end of a tier while other workers sit idle.
  size_order   <- order(origin_cells, decreasing = TRUE)
  origins      <- origins[size_order]
  origin_cells <- origin_cells[size_order]
  big_origins   <- origins[origin_cells > big_origin_cells]
  small_origins <- origins[origin_cells <= big_origin_cells]

  message("Computing cost-distances for: ", sp, " (", nrow(candidate_pairs), " candidate pairs, ",
          n_origins, " origin nodes: ", length(small_origins), " normal on ", n_cores, " worker(s), ",
          length(big_origins), " big-crop outlier(s) on ", big_n_cores, " worker(s))")

  results <- list()
  timing_results <- list()

  ## Runs one tier's origins through parLapplyLB() in batches on the given
  ## cluster, appending each batch's results/timing. If a worker dies
  ## mid-batch (e.g. an unusually large search area exhausted its
  ## memory), parLapplyLB() throws rather than returning a result --
  ## rebuild that tier's cluster and mark this batch's origins as failed
  ## rather than letting that abort the whole run. Returns the rebuilt
  ## cluster (unchanged if no failure occurred).
  run_tier <- function(tier_cl, tier_origins, tier_batch_size, tier_ncores, tier_label) {
    if (length(tier_origins) == 0) return(tier_cl)
    batches <- split(seq_along(tier_origins), ceiling(seq_along(tier_origins) / tier_batch_size))
    n_done <- 0
    for (b in batches) {
      batch_results <- tryCatch(
        parallel::parLapplyLB(tier_cl, tier_origins[b], cost_distance_for_origin),
        error = function(e) {
          warning("  ", sp, " (", tier_label, "): cluster worker failed on this batch (",
                  conditionMessage(e), ") -- rebuilding the cluster and marking these ",
                  length(b), " origin(s) as failed, then continuing.")
          ## The dead worker breaks parLapplyLB() for the whole batch, but
          ## the OTHER workers in this cluster are still alive -- stop it
          ## explicitly before replacing it, or those survivors leak as
          ## orphan processes that keep consuming memory on top of the
          ## new cluster, compounding with every subsequent failure.
          try(parallel::stopCluster(tier_cl), silent = TRUE)
          tier_cl <<- setup_costdist_cluster(tier_ncores, memmax_gb)
          parallel::clusterExport(tier_cl, varlist = c("candidate_pairs", "resistance_path"), envir = cd_env)
          vector("list", length(b))
        }
      )
      results       <<- c(results, lapply(batch_results, function(x) if (is.null(x)) NULL else x$data))
      timing_results <<- c(timing_results, lapply(batch_results, function(x) if (is.null(x)) NULL else x$timing))
      n_done <- n_done + length(b)
      message("  ", sp, " (", tier_label, "): ", n_done, "/", length(tier_origins), " origin(s) processed")
      ## Persistent PSOCK workers are reused across every batch in this
      ## tier -- terra's C-level allocations aren't reliably reclaimed by
      ## R's GC between calls, so a long tier can creep in memory even
      ## without any single origin being oversized. Force a GC on each
      ## surviving worker between batches as a cheap safety net.
      try(parallel::clusterEvalQ(tier_cl, invisible(gc())), silent = TRUE)
    }
    tier_cl
  }

  cl     <- run_tier(cl, small_origins, batch_size, n_cores, "normal")
  cl_big <- run_tier(cl_big, big_origins, big_n_cores, big_n_cores, "big-crop")

  results <- results[!vapply(results, is.null, logical(1))]
  timing_results <- timing_results[!vapply(timing_results, is.null, logical(1))]

  list(costs = bind_rows(results), timing = bind_rows(timing_results), cl = cl, cl_big = cl_big)
}

## 5. MST (04) #################################################################

## Per-edge distribution stats for comparing species (written to
## mst_edge_stats.csv; mst_summary.csv keeps the network-level totals).
## Edge cost and detour factor are strongly right-skewed (SD > mean), so
## describe them with median [q25, q75] and min-max; use mean with its 95% CI
## (t-interval -- matches a bootstrap CI to ~1% at n ~ 550) only to compare
## means between species. The CI treats the MST edges as a sample, but edges
## are spatially non-independent, so it's optimistic -- a guide, not a test.
describe_values <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0) {
    return(tibble(n = 0L, mean = NA_real_, sd = NA_real_, ci95_low = NA_real_,
                  ci95_high = NA_real_, median = NA_real_, q25 = NA_real_,
                  q75 = NA_real_, min = NA_real_, max = NA_real_))
  }
  ci_half <- if (n > 1) qt(0.975, df = n - 1) * sd(x) / sqrt(n) else NA_real_
  tibble(
    n         = n,
    mean      = mean(x),
    sd        = if (n > 1) sd(x) else NA_real_,
    ci95_low  = mean(x) - ci_half,
    ci95_high = mean(x) + ci_half,
    median    = median(x),
    q25       = quantile(x, 0.25, names = FALSE),
    q75       = quantile(x, 0.75, names = FALSE),
    min       = min(x),
    max       = max(x)
  )
}

## PA areas in km^2, named by site_id. Uses the stored `area` column (the
## true area, computed in 00 before the 30 m polygon simplification), not
## st_area() of the simplified geometry -- they differ by up to ~7 km^2 per PA.
pa_area_lookup <- function(pa_patches) {
  setNames(pa_patches$area / 1e6, as.character(pa_patches$site_id))
}

## Protected area held by an MST's network (one-row tibble):
##   largest_component_area_* -- PA area in the biggest connected component
##     (and its share of all PA area); < 1 only when the MST is a forest
##   isolated_pa_area_km2 -- area of PAs with no MST edge at all
##   eca_* -- equivalent connected area (Saura et al. 2011):
##     ECA = sqrt(sum_ij a_i a_j p*_ij), p*_ij = product of the link
##     probabilities along the path from i to j (0 between components;
##     the i = j terms are a_i^2). On a tree there is exactly one path per
##     pair, so this is a lower bound on the candidate graph's ECA. Link
##     probability is the euclidean-distance kernel (same as link_prob), so
##     ECA responds to dispersal distance and PA sizes but barely to
##     resistance. Its floor (no links at all) is sqrt(sum a_i^2), ~0.19 of
##     total PA area for the 5 km^2 PA set -- read eca_frac against that,
##     not against 0.
## `mst_g` needs a `dispersal_probability` edge attribute and vertex names
## that are site_ids in `pa_area_km2` (from pa_area_lookup()).
mst_area_stats <- function(mst_g, pa_area_km2) {
  a <- unname(pa_area_km2[igraph::V(mst_g)$name])
  if (anyNA(a)) stop("mst_area_stats(): no PA area for some site_ids in the MST")
  total_area <- sum(a)
  comp_area <- tapply(a, igraph::components(mst_g)$membership, sum)
  link_p <- igraph::E(mst_g)$dispersal_probability
  ## -log turns the product of probabilities along a path into a sum, so
  ## shortest-path distances give p*_ij = exp(-D_ij)
  D <- igraph::distances(mst_g, weights = if (length(link_p) > 0) -log(link_p) else NA)
  eca <- sqrt(sum(outer(a, a) * exp(-D)))  # exp(-Inf) = 0 between components
  tibble(
    total_pa_area_km2           = total_area,
    largest_component_area_km2  = max(comp_area),
    largest_component_area_frac = max(comp_area) / total_area,
    isolated_pa_area_km2        = sum(a[igraph::degree(mst_g) == 0]),
    eca_km2                     = eca,
    eca_frac                    = eca / total_area
  )
}

## Re-apply the dispersal cutoff to the costed pairs, then build the MST
## (or minimum spanning FOREST if disconnected) with cost_dist as the edge
## weight (Prim's algorithm). `pair_costs` is a {sp}_PA_costs.csv table,
## `site_ids` every PA node, `pa_area_km2` their areas (pa_area_lookup()).
## Returns list(graph, mst, edges, summary, edge_stats) -- `summary` is this
## species' mst_summary.csv row and `edge_stats` its mst_edge_stats.csv rows.
build_mst <- function(pair_costs, site_ids, alpha, p_threshold, sp = "", pa_area_km2) {
  edges <- pair_costs %>%
    filter(!is.na(cost_dist)) %>%
    ## SHOULD WE USE THE CUTOFF FROM EUCL. DISTANCE? COST DIST DOESNT REALLY MAKE SENSE...CAN BE EITHER ONE
    mutate(dispersal_probability = dispersal_probability(euclidean_dist, alpha)) %>%
    filter(dispersal_probability >= p_threshold)

  g <- igraph::graph_from_data_frame(
    edges %>% dplyr::select(from, to, weight = cost_dist, euclidean_dist, dispersal_probability),
    directed = FALSE,
    vertices = site_ids
  )

  n_components <- igraph::components(g)$no
  if (n_components > 1) {
    message("  ", n_components, " disconnected components -- result is a ",
            "minimum spanning FOREST, not a single tree.")
  }

  ## Generate the MST using the graph of PA nodes and the weighted edges (of cost)
  # The minimum spanning tree will have the smallest sum of edge weights (lowest costs)
  mst_g <- igraph::mst(g, weights = igraph::E(g)$weight, algorithm = "prim")
  mst_edges <- igraph::as_data_frame(mst_g, what = "edges") %>% as_tibble()

  weakest_link_prob <- if (nrow(mst_edges) > 0) min(mst_edges$dispersal_probability) else NA_real_
  mean_link_prob <- if (nrow(mst_edges) > 0) mean(mst_edges$dispersal_probability) else NA_real_

  ## detour factor: how much farther the cheapest cost-path is than a
  ## straight line, per MST edge -- flags edges where resistance forces a
  ## long detour relative to geography alone. NA (not Inf) for edges whose
  ## patches already touch (euclidean_dist = 0), since the ratio is
  ## undefined there.
  detour_factor <- ifelse(mst_edges$euclidean_dist > 0,
                           mst_edges$weight / mst_edges$euclidean_dist, NA_real_)
  n_detour_valid <- sum(!is.na(detour_factor))
  mst_degree <- igraph::degree(mst_g)
  ## cost-weighted diameter: worst-case cumulative cost between the two
  ## most distant connected PAs (NA, not 0, when there are no MST edges).
  mst_diameter_cost <- if (nrow(mst_edges) > 0) igraph::diameter(mst_g, weights = igraph::E(mst_g)$weight) else NA_real_

  summary_row <- tibble(
    species                    = sp,
    kernel_p_threshold         = p_threshold,
    n_nodes                    = igraph::vcount(g),
    n_candidate_edges          = igraph::ecount(g),
    n_isolated_nodes           = sum(igraph::degree(g) == 0),  # nodes with no candidate pair within dispersal range
    n_mst_edges                = igraph::ecount(mst_g),
    n_components               = n_components,
    total_mst_cost             = sum(igraph::E(mst_g)$weight),
    mean_mst_edge_cost         = if (nrow(mst_edges) > 0) mean(mst_edges$weight) else NA_real_,
    median_mst_edge_cost       = if (nrow(mst_edges) > 0) median(mst_edges$weight) else NA_real_,
    max_mst_edge_cost          = if (nrow(mst_edges) > 0) max(mst_edges$weight) else NA_real_,
    total_mst_euclidean_dist_m = if (nrow(mst_edges) > 0) sum(mst_edges$euclidean_dist) else NA_real_,
    mean_detour_factor         = if (n_detour_valid > 0) mean(detour_factor, na.rm = TRUE) else NA_real_,
    max_detour_factor          = if (n_detour_valid > 0) max(detour_factor, na.rm = TRUE) else NA_real_,
    mst_diameter_cost          = mst_diameter_cost,
    max_node_degree            = if (length(mst_degree) > 0) max(mst_degree) else NA_real_,
    weakest_link_prob          = weakest_link_prob,  # lowest dispersal probability among the MST's own edges
    mean_link_prob             = mean_link_prob
  ) %>%
    bind_cols(mst_area_stats(mst_g, pa_area_km2))

  edge_stats <- bind_rows(
    describe_values(mst_edges$weight)                %>% mutate(metric = "edge_cost",           unit = "resistance x m"),
    describe_values(mst_edges$euclidean_dist)        %>% mutate(metric = "edge_euclidean_dist", unit = "m"),
    describe_values(detour_factor)                   %>% mutate(metric = "detour_factor",       unit = "cost / m"),
    describe_values(mst_edges$dispersal_probability) %>% mutate(metric = "link_prob",           unit = "probability")
  ) %>%
    mutate(species = sp) %>%
    relocate(species, metric, unit)

  list(graph = g, mst = mst_g, edges = mst_edges, summary = summary_row, edge_stats = edge_stats)
}

## MST edges as boundary-to-boundary (edge-to-edge) STRAIGHT lines: each
## runs between the nearest points on the two patches' actual boundaries
## (sf::st_nearest_points(), same endpoint method as trace_edge_lcp()).
## Not the routed least-cost corridor -- that's 05's job.
mst_edge_lines <- function(mst_edges, pa_patches) {
  edge_lines <- Map(function(from_id, to_id) {
    origin_poly <- pa_patches[as.character(pa_patches$site_id) == as.character(from_id), ]
    dest_poly   <- pa_patches[as.character(pa_patches$site_id) == as.character(to_id), ]
    st_nearest_points(origin_poly, dest_poly)
  }, mst_edges$from, mst_edges$to)

  mst_lines <- mst_edges
  st_geometry(mst_lines) <- do.call(c, edge_lines)
  st_crs(mst_lines) <- st_crs(pa_patches)
  mst_lines
}

## 6. Least-cost-path corridors (05) ###########################################
## See 05's header for why the transition surface is built per edge on a
## buffered, corridor-masked crop, and for the resistance -> conductance
## convention.

## Attempt to trace one edge's LCP on a crop of the raster at `resistance_path`,
## buffered by `buffer_mult` x the edge's normal buffer. Returns NULL (never
## throws) if the crop, transition build, or path tracing fails -- safe to
## call from a parallel worker.
attempt_edge_lcp <- function(from_xy, to_xy, resistance_path, buffer_mult) {
  edge_dist <- sqrt(sum((from_xy - to_xy)^2))
  buffer_dist <- min(max(edge_dist * lcp_buffer_frac, lcp_buffer_min_m), lcp_buffer_max_m) * buffer_mult

  edge_ext <- ext(
    min(from_xy[, 1], to_xy[, 1]) - buffer_dist,
    max(from_xy[, 1], to_xy[, 1]) + buffer_dist,
    min(from_xy[, 2], to_xy[, 2]) - buffer_dist,
    max(from_xy[, 2], to_xy[, 2]) + buffer_dist
  )

  tryCatch({
    resistance <- rast(resistance_path)
    resistance_crop <- crop(resistance, edge_ext)

    ## Mask down to a corridor around the straight line between the two
    ## nodes, not the full rectangle -- transition() excludes NA cells from
    ## the graph it builds, so this keeps graph size proportional to edge
    ## length rather than edge length squared.
    edge_line <- vect(rbind(from_xy, to_xy), type = "lines", crs = crs(resistance))
    edge_corridor <- buffer(edge_line, width = buffer_dist)
    resistance_crop <- mask(resistance_crop, edge_corridor)

    resistance_crop_r <- raster::raster(resistance_crop)

    tr <- transition(resistance_crop_r, transitionFunction = function(x) 1 / mean(x), directions = 8)
    tr <- geoCorrection(tr, type = "c", multpl = FALSE)

    path <- gdistance::shortestPath(tr, from_xy, to_xy, output = "SpatialLines")
    path_sf <- st_as_sf(path)
    st_crs(path_sf) <- crs(resistance)

    list(path = path_sf, resistance_crop = resistance_crop)
  }, error = function(e) NULL)
}

## Trace one MST edge's LCP, retrying once with a doubled buffer on failure.
## Edge-to-edge (patch boundary to patch boundary): from_xy/to_xy are the
## geometrically nearest points between the two PA polygons' boundaries
## (sf::st_nearest_points()), not centroids. `pa_patches` is a free variable.
trace_edge_lcp <- function(from_id, to_id, resistance_path) {
  origin_poly <- pa_patches[as.character(pa_patches$site_id) == from_id, ]
  dest_poly   <- pa_patches[as.character(pa_patches$site_id) == to_id, ]

  nearest_line <- sf::st_nearest_points(origin_poly, dest_poly)
  nearest_coords <- sf::st_coordinates(nearest_line)

  from_xy <- nearest_coords[1, c("X", "Y"), drop = FALSE]
  to_xy   <- nearest_coords[2, c("X", "Y"), drop = FALSE]

  result <- attempt_edge_lcp(from_xy, to_xy, resistance_path, buffer_mult = 1)
  if (is.null(result)) {
    result <- attempt_edge_lcp(from_xy, to_xy, resistance_path, buffer_mult = 2)
  }
  result
}

## Worker: one MST edge's full task -- trace it, then attach its edge
## metadata. References `mst_edges` / `resistance_path` as free variables
## (exported by trace_lcp_corridors()).
lcp_worker <- function(j) {
  from_id <- as.character(mst_edges$from[j])
  to_id   <- as.character(mst_edges$to[j])

  result <- trace_edge_lcp(from_id, to_id, resistance_path)
  if (is.null(result)) return(NULL)

  lcp_sf <- result$path
  lcp_sf$from <- from_id
  lcp_sf$to <- to_id
  lcp_sf$cost_dist <- mst_edges$weight[j]
  lcp_sf$dispersal_probability <- mst_edges$dispersal_probability[j]
  lcp_sf
}

## Builds a fresh LCP cluster with everything that doesn't change across
## species already exported (needs `pa_patches` in the global environment).
setup_lcp_cluster <- function(n_workers) {
  new_cl <- parallel::makeCluster(n_workers)
  parallel::clusterEvalQ(new_cl, {
    library(terra); library(sf); library(gdistance); library(raster)
  })
  parallel::clusterExport(new_cl, varlist = c(
    "pa_patches", "lcp_buffer_frac", "lcp_buffer_min_m", "lcp_buffer_max_m",
    "attempt_edge_lcp", "trace_edge_lcp", "lcp_worker"
  ))
  new_cl
}

## Trace the LCPs of MST edges `edge_ids` (rows of `mst_edges`, which needs
## from/to/weight/dispersal_probability) across `cl`, in batches with
## progress messages. If a worker dies mid-batch (e.g. an outlier-long edge
## exhausted its memory), parLapplyLB() throws rather than returning a
## result -- the cluster is rebuilt and that batch marked failed rather
## than aborting the run. Returns list(paths, cl): `paths` is one entry per
## row of mst_edges (NULL where failed or not in edge_ids); keep the
## returned cluster.
trace_lcp_corridors <- function(mst_edges, resistance_path, cl, n_workers,
                                edge_ids = seq_len(nrow(mst_edges)), batch_size = 100, sp = "") {
  lcp_env <- environment()
  ## `mst_edges` / `resistance_path` change every species -- re-export.
  parallel::clusterExport(cl, varlist = c("mst_edges", "resistance_path"), envir = lcp_env)

  n_edges <- nrow(mst_edges)
  edge_paths <- vector("list", n_edges)

  batches <- split(edge_ids, ceiling(seq_along(edge_ids) / batch_size))
  n_done <- n_edges - length(edge_ids)

  for (b in batches) {
    batch_results <- tryCatch(
      parallel::parLapplyLB(cl, b, lcp_worker),
      error = function(e) {
        warning("  ", sp, ": cluster worker failed on this batch (", conditionMessage(e),
                ") -- rebuilding the cluster and marking these ", length(b),
                " edge(s) as failed, then continuing.")
        cl <<- setup_lcp_cluster(n_workers)
        ## setup_lcp_cluster() only exports the species-invariant objects --
        ## this species' mst_edges/resistance_path need re-exporting to
        ## the new cluster too, or the next batch fails the same way.
        parallel::clusterExport(cl, varlist = c("mst_edges", "resistance_path"), envir = lcp_env)
        vector("list", length(b))
      }
    )
    edge_paths[b] <- batch_results
    n_done <- n_done + length(b)
    message("  ", sp, ": ", n_done, "/", n_edges, " edge(s) traced")
  }

  list(paths = edge_paths, cl = cl)
}

## 7. Corridor density (06) ####################################################
## See 06's header for the SIGMA / SCALE / PARALLELIZATION design notes.

## Worker: one MST edge's density contribution -- rebuild its geometry and
## this species' coarse grid from plain (serializable) values, rasterize +
## smooth on a small crop around the edge, write the result to a temp .tif,
## and return the path (a SpatRaster can't be returned from a worker).
## References `edge_coords_list`, `sp_grid_ext_vec`, `sp_grid_res`,
## `crop_buffer_m`, `gauss_1d`, `cell_area_km2`, `template_crs_wkt` as free
## variables (exported by setup_density_cluster() / corridor_density()).
edge_density_worker <- function(e) {
  edge_v <- vect(edge_coords_list[[e]], type = "lines", crs = template_crs_wkt)
  grid_template <- rast(ext(sp_grid_ext_vec), resolution = sp_grid_res, crs = template_crs_wkt)

  crop_ext <- ext(edge_v) + crop_buffer_m
  grid_crop <- crop(grid_template, crop_ext)

  length_per_cell_km <- rasterizeGeom(edge_v, grid_crop, fun = "length", unit = "km")
  density_raw <- length_per_cell_km / cell_area_km2  # km corridor per km^2, unsmoothed

  ## Separable Gaussian smoothing: one pass along rows, one along columns.
  density_smoothed <- focal(density_raw, w = matrix(gauss_1d, nrow = 1), fun = "sum", na.rm = TRUE)
  density_smoothed <- focal(density_smoothed, w = matrix(gauss_1d, ncol = 1), fun = "sum", na.rm = TRUE)

  out_path <- tempfile(fileext = ".tif")
  writeRaster(density_smoothed, out_path, overwrite = TRUE)
  out_path
}

## Builds a fresh density cluster: workers get the species-invariant pieces
## (CRS, the worker function); everything that varies per species is
## exported by corridor_density().
setup_density_cluster <- function(n_workers) {
  template_crs_wkt <- crs(template_raster)
  new_cl <- parallel::makeCluster(n_workers)
  parallel::clusterExport(new_cl, varlist = c("template_crs_wkt", "edge_density_worker"),
                          envir = environment())
  parallel::clusterEvalQ(new_cl, library(terra))
  new_cl
}

## Kernel-smoothed corridor density for one species' traced LCPs
## (`corridors`, sf), written to `out_path`. Sigma = corridor_sigma(dispersal_km),
## on the species' own grid (~sigma / cells_per_sigma, snapped to a multiple
## of 30 m, but never coarser than common_res_m = the SDMs' 1 km, so the 07
## overlay keeps the SDM's detail). Each edge is rasterized + smoothed on its own crop in
## parallel, then combined with mosaic(fun = "sum") -- valid because
## Gaussian smoothing is linear. Returns list(path, cl): path is NULL if no
## edge could be rasterized; keep the returned cluster.
corridor_density <- function(corridors, dispersal_km, cl, n_workers, out_path,
                             batch_size = 50, sp = "") {
  dens_env <- environment()
  template_crs_wkt <- crs(template_raster)

  corridors_v <- vect(corridors)
  n_edges <- nrow(corridors_v)

  ## Plain coordinate matrices, one per edge -- what actually gets
  ## exported to the cluster (SpatVector itself can't be).
  edge_coords_list <- lapply(seq_len(n_edges), function(e) {
    terra::geom(corridors_v[e, ])[, c("x", "y")]
  })

  sigma <- corridor_sigma(dispersal_km)
  ## capped at common_res_m: a capped species' grid is then exactly the
  ## common 1 km grid (same extent/resolution as to_common_grid())
  sp_grid_res_val <- min(common_res_m, max(30, round(sigma / cells_per_sigma / 30) * 30))
  sp_grid_res <- c(sp_grid_res_val, sp_grid_res_val)
  sp_grid_ext_vec <- as.vector(ext(template_raster))
  crop_buffer_m <- 3 * sigma
  cell_area_km2 <- (sp_grid_res_val * sp_grid_res_val) / 1e6

  ## 1D Gaussian profile, for separable smoothing. Extracted from terra's own
  ## 2D matrix (its center row, renormalized to sum to 1) rather than
  ## re-derived by hand, so it reproduces exactly what the 2D kernel would
  ## have -- a separable Gaussian's 2D matrix is the outer product of this
  ## profile with itself.
  sp_grid_template <- rast(ext(sp_grid_ext_vec), resolution = sp_grid_res, crs = template_crs_wkt)
  gauss_2d <- focalMat(sp_grid_template, sigma, "Gauss")
  gauss_1d <- gauss_2d[(nrow(gauss_2d) + 1) / 2, ]
  gauss_1d <- gauss_1d / sum(gauss_1d)

  sp_export_vars <- c("edge_coords_list", "sp_grid_ext_vec", "sp_grid_res",
                      "crop_buffer_m", "gauss_1d", "cell_area_km2")
  parallel::clusterExport(cl, varlist = sp_export_vars, envir = dens_env)

  message("Computing corridor density for: ", sp, " (", n_edges, " edge(s), ", n_workers, " worker(s)) -- ",
          "dispersal_km = ", dispersal_km, ", sigma = ", round(sigma), " m, ",
          "grid res = ", sp_grid_res_val, " m, crop buffer = ", round(crop_buffer_m), " m")

  edge_tif_paths <- vector("list", n_edges)
  batches <- split(seq_len(n_edges), ceiling(seq_len(n_edges) / batch_size))
  n_done <- 0

  for (b in batches) {
    batch_results <- tryCatch(
      parallel::parLapplyLB(cl, b, edge_density_worker),
      error = function(e) {
        warning("  ", sp, ": cluster worker failed on this batch (", conditionMessage(e),
                ") -- rebuilding the cluster and marking these ", length(b),
                " edge(s) as failed, then continuing.")
        cl <<- setup_density_cluster(n_workers)
        parallel::clusterExport(cl, varlist = sp_export_vars, envir = dens_env)
        vector("list", length(b))
      }
    )
    edge_tif_paths[b] <- batch_results
    n_done <- n_done + length(b)
    message("  ", sp, ": ", n_done, "/", n_edges, " edge(s) processed")
  }

  n_failed <- sum(vapply(edge_tif_paths, is.null, logical(1)))
  edge_tif_paths <- edge_tif_paths[!vapply(edge_tif_paths, is.null, logical(1))]

  if (length(edge_tif_paths) == 0) {
    warning(sp, ": no edges could be rasterized -- nothing saved.")
    return(list(path = NULL, cl = cl, n_failed = n_failed))
  }

  edge_density_rasters <- lapply(edge_tif_paths, rast)
  density_smoothed <- mosaic(sprc(edge_density_rasters), fun = "sum")

  writeRaster(density_smoothed, out_path, overwrite = TRUE)
  unlink(unlist(edge_tif_paths))  # done reading these back -- clean up the temp per-edge rasters
  message("  -> saved ", out_path, " (", length(edge_density_rasters), " of ", n_edges, " edge(s) combined",
          if (n_failed > 0) paste0(", ", n_failed, " skipped") else "", ")")

  list(path = out_path, cl = cl, n_failed = n_failed)
}

## Resample a raster onto the shared `common_res_m` grid over the full
## modelling extent (1 km = the SDMs' native resolution) -- averaging when
## downsampling (raster finer than the grid), bilinear when upsampling.
## Used for the multi-species map (06) and for comparing sensitivity runs.
to_common_grid <- function(r, res_m = common_res_m) {
  common_grid <- rast(ext(template_raster), resolution = res_m, crs = crs(template_raster))
  resample_method <- if (xres(r) <= res_m) "average" else "bilinear"
  resample(r, common_grid, method = resample_method)
}

## 8. SDM overlap (07) ##########################################################

rescale_raster <- function(r) {
  mm <- minmax(r)
  (r - mm[1, 1]) / (mm[2, 1] - mm[1, 1])
}

## Overlay a corridor density raster with a (pre-cropped, script 00) SDM:
##   1. drop the bottom `quantile_cutoff` of corridor density values as
##      noise (set to NA), then rescale the remainder to [0, 1]
##   2. reproject/crop/mask the SDM onto the corridor raster's grid
##   3. multiply corridor x SDM suitability, rescale the product to [0, 1]
sdm_overlap <- function(corridor, sdm_path, quantile_cutoff) {
  thresh <- as.numeric(global(corridor, fun = quantile, probs = quantile_cutoff, na.rm = TRUE))
  corridor_q <- ifel(corridor <= thresh, NA, corridor)
  corridor_q <- rescale_raster(corridor_q)

  # bring the SDM onto the corridor raster's (already small) grid
  sdm <- rast(sdm_path)
  sdm <- project(sdm, corridor_q)
  sdm <- crop(sdm, corridor_q)
  sdm <- mask(sdm, corridor_q)

  ### Combine corridor and species suitability by multiplication
  ## This assumes that connectivity is modulated by habitat suitability.
  #-#-#-# NOTE THIS CAN BE CHANGED TO WHATEVER IS APPROPRIATE #-#-#-#
  overlap <- corridor_q * sdm #straight multiplication used currently
  rescale_raster(overlap)
}

## 9. Sensitivity analysis (08/09) #############################################

## Comparisons of a sensitivity scenario with production. All return
## similarities in [0, 1] (1 = identical), so 1 - value is the change.

## Undirected edge keys ("small-large" site_id pairs) for an MST edge table.
mst_edge_keys <- function(edges) {
  a <- as.numeric(edges$from); b <- as.numeric(edges$to)
  paste(pmin(a, b), pmax(a, b), sep = "-")
}

## Topology agreement between two MSTs over the same node set:
##   edge_jaccard -- |shared edges| / |union of edges|
##   ref_edge_retention -- share of `ref`'s edges also in `x`
##   degree_spearman / betweenness_spearman -- rank agreement of per-node
##     degree / betweenness (hub and stepping-stone stability); NA when a
##     vector is constant (e.g. no edges)
compare_msts <- function(x, ref, site_ids) {
  kx <- mst_edge_keys(x); kr <- mst_edge_keys(ref)
  n_union <- length(union(kx, kr))
  n_shared <- length(intersect(kx, kr))
  node_stats <- function(edges) {
    g <- igraph::graph_from_data_frame(edges[, c("from", "to")], directed = FALSE,
                                       vertices = data.frame(name = site_ids))
    list(degree = igraph::degree(g), betweenness = igraph::betweenness(g, weights = NA))
  }
  safe_spearman <- function(a, b) {
    if (sd(a) == 0 || sd(b) == 0) return(NA_real_)
    cor(a, b, method = "spearman")
  }
  sx <- node_stats(x); sr <- node_stats(ref)
  tibble(
    edge_jaccard         = if (n_union == 0) 1 else n_shared / n_union,
    ref_edge_retention   = if (length(kr) == 0) 1 else n_shared / length(kr),
    degree_spearman      = safe_spearman(sx$degree, sr$degree),
    betweenness_spearman = safe_spearman(sx$betweenness, sr$betweenness)
  )
}

## A run's LCP network prepared for compare_corridors(): the dissolved
## lines (so corridor stretches shared by several MST edges count once),
## their `buffer_m` tolerance zone, and total dissolved length. Built once
## per run -- the buffer is the expensive part, so it's built from lines
## simplified to buffer_m / 10 (LCPs have a vertex every 30 m cell; that
## detail is irrelevant to a zone buffer_m wide, and the simplified line
## stays within buffer_m / 10 of the original).
corridor_set <- function(lcp, buffer_m = 1000) {
  lines <- st_union(st_geometry(lcp))
  list(lines = lines, zone = st_buffer(st_simplify(lines, dTolerance = buffer_m / 10), buffer_m),
       length = sum(as.numeric(st_length(lines))))
}

## Corridor-location agreement (Beier et al. 2009-style): share of one LCP
## network's (dissolved) length lying within the buffer of the other's,
## averaged over both directions (symmetric). Takes corridor_set()s; NULL
## = a run with no corridors.
compare_corridors <- function(a, b) {
  if (is.null(a) || is.null(b)) {
    return(tibble(corridor_agreement = if (is.null(a) && is.null(b)) 1 else 0))
  }
  share_within <- function(x, y) sum(as.numeric(st_length(st_intersection(x$lines, y$zone)))) / x$length
  tibble(corridor_agreement = mean(c(share_within(a, b), share_within(b, a))))
}

## Agreement between two surfaces on the SAME grid (to_common_grid() output),
## NA treated as 0 (no corridor), over cells positive in either:
##   spearman -- cell-wise rank correlation (descriptive only: cells are
##               spatially autocorrelated, so no p-value)
##   schoener_d -- 1 - 0.5 * sum|p - q| with each surface normalized to sum 1
##               (Warren, Glor & Turelli 2008)
##   top_jaccard -- Jaccard overlap of each surface's top `top_frac` of
##               positive cells (hotspot stability)
compare_surfaces <- function(r_a, r_b, top_frac = 0.1) {
  if (!compareGeom(r_a, r_b, stopOnError = FALSE)) stop("compare_surfaces(): rasters must share one grid -- use to_common_grid() first.")
  a <- values(r_a, mat = FALSE); b <- values(r_b, mat = FALSE)
  a[is.na(a)] <- 0; b[is.na(b)] <- 0
  keep <- a > 0 | b > 0
  a <- a[keep]; b <- b[keep]
  if (length(a) == 0 || sum(a) == 0 || sum(b) == 0) {
    return(tibble(spearman = NA_real_, schoener_d = NA_real_, top_jaccard = NA_real_))
  }
  top_set <- function(v) {
    pos <- which(v > 0)
    pos[v[pos] >= quantile(v[pos], 1 - top_frac, names = FALSE)]
  }
  ta <- top_set(a); tb <- top_set(b)
  tibble(
    spearman    = if (sd(a) == 0 || sd(b) == 0) NA_real_ else cor(a, b, method = "spearman"),
    schoener_d  = 1 - 0.5 * sum(abs(a / sum(a) - b / sum(b))),
    top_jaccard = length(intersect(ta, tb)) / length(union(ta, tb))
  )
}
