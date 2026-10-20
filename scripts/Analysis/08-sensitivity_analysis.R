## ------------------------------- ##
## Script name: 08-sensitivity_runs
##
## Purpose of script: One-at-a-time sensitivity analysis of dispersal
# distance and resistance contrast. Each SCENARIO changes exactly one of
# them from the production setup, reruns only the stages it affects
# (through the same functions the stage scripts use, in
# pipeline_functions.R), and is compared with the PRODUCTION outputs of
# scripts 04/05/07. 09 turns the results into tables and figures.
#
# SCENARIOS (per species):
#   dispersal_mult   dispersal_km x 0.5, 0.75, 1.33, 2. Moves the pair
#                    cutoff (01/04) and the corridor smoothing sigma (06).
#                    (p_threshold needs no scenarios of its own: the cutoff
#                    radius is dispersal_km * log2(1/p), so p = 0.01 is the
#                    same cutoff as dispersal x1.54, p = 0.10 as x0.77.)
#   contrast_mult    resistance = 2^(m * rank), m = 0.5, 0.75, 1.25, 1.5 --
#                    how steep the doubling curve is (m < 1 flattens the
#                    contrast between suitable and unsuitable cover, m > 1
#                    sharpens it).
# Literature basis: baseline-centred one-at-a-time perturbation of
# resistance values and dispersal parameters, as in Beier et al. 2009 (Ecol
# Appl 19:2067) and Rayfield et al. 2010 (Landscape Ecol 25:519).
#
# WHAT EACH SCENARIO RERUNS (everything else is read from production):
#   dispersal scenarios -- the resistance surface doesn't change, so the
#     production resistance raster (02), cost-distances (03) and corridors
#     (05) are reused. A shorter cutoff just drops pairs from the production
#     cost table; a longer one finds the extra pairs and computes
#     cost-distance for those only. Only MST edges production didn't trace
#     get a new LCP. Then corridor density (with this dispersal's sigma)
#     and SDM overlap.
#   contrast scenarios -- same candidate pairs as production (01). Since
#     2^(m * rank) = (2^rank)^m, the scenario's resistance raster is simply
#     the production raster raised to m (exact: ranks are integers). Then
#     cost-distance on all pairs, the MST, every corridor, density and
#     overlap are recomputed.
#
# OUTPUTS -- three files per scenario, in Intermediate/08_sensitivity/{species}/:
#   {scenario}_result.csv      one row: the input changed, the MST's shape,
#                              and its similarity to production (written
#                              LAST -- its presence marks the scenario done)
#   {scenario}_mst_edges.csv   the scenario's MST edges
#   {scenario}_overlap_1km.tif SDM overlap on the common 1 km grid
#                              ({scenario}_density_1km.tif for species with
#                              no SDM -- GLIN, ASFL)
# Everything else (resistance raster, cost table, corridors, full-resolution
# density/overlap) is written to a temp folder and deleted.
#
# Production outputs (01-07) must be current for the species first.
#
# RUNNING ON FIR: Scripts/Analysis/fir_sensitivity.sh runs one scenario per SLURM
# array task (the numbered task list printed below). Copy Scripts/ and the
# species' production outputs there first, and rerun
# setup_species_from_xlsx.R on Fir so species_metadata.csv's SDM paths are
# Fir paths (this script stops if they point to a missing file). Worker
# counts, terra memory and temp files adapt to the job automatically.
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))

## ---- Toggles ----
species_subset <- c("MAAM", "ASFL", "DOOR", "DRPI", "AQCH")   # NULL = all species
scenario_subset <- NULL       # e.g. c("disp_x0.5", "disp_x0.75") for a quick test; NULL = all
overwrite <- FALSE            # TRUE reruns scenarios that already have a _result.csv
## Test only: add an "identity_check" scenario (contrast m = 1) that changes
## nothing but runs the full 03 -> 07 path. It should reproduce production
## (edge_jaccard = 1, corridor_agreement = 1, surface_schoener_d = 1).
include_identity_check <- FALSE

## ---- Scenario levels ----
disp_mults <- c(0.5, 0.75, 4 / 3, 2)
contrast_mults <- c(0.5, 0.75, 1.25, 1.5)

## ---- Compute knobs ----
## On Fir (any SLURM job) worker counts follow the job's --cpus-per-task,
## each cost-distance worker's terra memory is capped to its share of the
## job's allocation (setup_costdist_cluster()), and temp files go to the
## node's fast local disk ($SLURM_TMPDIR). On this desktop: few cost-
## distance workers -- they commit up to ~14 GB each on MAAM's biggest
## crops, and the machine has 32 GB.
on_slurm <- nzchar(Sys.getenv("SLURM_JOB_ID"))
n_cpus <- if (on_slurm) as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")) else parallel::detectCores()
cd_big_n_cores <- if (on_slurm) 2 else 1
cd_n_cores <- if (on_slurm) max(1, n_cpus - cd_big_n_cores - 1) else 2  # e.g. 13 CPUs -> 10 + 2 big + master
pairs_n_cores <- if (on_slurm) max(1, n_cpus - 1) else 4
lcp_n_cores <- if (on_slurm) max(1, n_cpus - 1) else 8
dens_n_cores <- lcp_n_cores
cd_batch_size <- 25
cd_big_origin_cells <- 1e8
lcp_batch_size <- 100
edge_batch_size <- 50

## Where each scenario's temp folder (resistance raster, costs, density)
## and terra's own temp files go
scratch_root <- if (on_slurm && nzchar(Sys.getenv("SLURM_TMPDIR"))) Sys.getenv("SLURM_TMPDIR") else tempdir()
terraOptions(tempdir = scratch_root)
if (on_slurm) {
  ## the master raises the full-extent resistance raster to m before any
  ## worker exists -- cap terra to the job's own allocation, not the node's RAM
  job_mem_gb <- slurm_job_mem_gb()
  if (!is.na(job_mem_gb)) terraOptions(memmax = 0.5 * job_mem_gb)
}

## ---- Inputs ----
pa_patches <- read_sf(here(pa_patches_dir, "pa_patches.gpkg"))
pa_sf <- pa_patches  # find_candidate_pairs()' workers read it as `pa_sf`
site_ids <- pa_patches$site_id
pa_area_km2 <- pa_area_lookup(pa_patches)
kernel_params <- read_csv(here(table_dir, "01_dispersal_kernel_params.csv"), show_col_types = FALSE)
species_meta <- read_csv(here(tabular_dir, "species_metadata.csv"), show_col_types = FALSE)

## ---- Scenario list (same for every species) ----
make_scenarios <- function() {
  bind_rows(
    tibble(parameter = "dispersal_mult", value = disp_mults),
    tibble(parameter = "contrast_mult", value = contrast_mults)
  ) %>%
    mutate(scenario = if_else(parameter == "dispersal_mult",
                              paste0("disp_x", signif(value, 3)),
                              paste0("contrast_", value)), .before = 1)
}

## ---- Production outputs a species' scenarios are compared with / reuse ----
load_production <- function(sp) {
  sdm_path <- species_meta$sdm_raster_path[species_meta$species == sp]
  if (length(sdm_path) == 1 && !is.na(sdm_path) && !file.exists(sdm_path)) {
    ## e.g. a Windows path in species_metadata.csv read on Fir -- don't
    ## silently fall back to comparing corridor density instead of overlap
    stop(sp, ": species_metadata.csv points to an SDM that isn't here (", sdm_path,
         ") -- rerun setup_species_from_xlsx.R on this machine.")
  }
  has_sdm <- length(sdm_path) == 1 && !is.na(sdm_path)
  files <- c(
    pairs     = here(pa_pairs_dir, paste0(sp, "_PApairs_dkernel_cutoff.csv")),
    res       = here(resistance_dir, paste0(sp, "_res.tif")),
    costs     = here(costdist_dir, paste0(sp, "_PA_costs.csv")),
    mst_edges = here(mst_dir, paste0(sp, "_mst_edges.csv")),
    lcps      = here(corridor_dir, paste0(sp, "_lcp_corridors.gpkg")),
    surface   = if (has_sdm) here(sdm_overlap_dir, paste0(sp, "_habitat_overlap.tif"))
                else here(density_dir, paste0(sp, "_corridor_density.tif"))
  )
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop(sp, ": production output(s) missing -- run 01-07 for this species first:\n  ",
         paste(missing, collapse = "\n  "))
  }
  dispersal_km <- kernel_params$kernel_ref_dist_m[kernel_params$species == sp] / 1000
  lcps <- st_read(files[["lcps"]], quiet = TRUE)
  list(
    dispersal_km = dispersal_km,
    kernel       = species_kernel(dispersal_km),
    res_path     = files[["res"]],
    pairs        = read.csv(files[["pairs"]]),
    pair_costs   = read_csv(files[["costs"]], show_col_types = FALSE),
    edges        = read_csv(files[["mst_edges"]], show_col_types = FALSE),
    lcp_by_key   = split(lcps, mst_edge_keys(lcps)),
    corr         = corridor_set(lcps),
    has_sdm      = has_sdm,
    sdm_path     = sdm_path,
    surface_type = if (has_sdm) "overlap" else "density",
    ## in memory: each scenario clears terra's temp files when it finishes
    surface_1km  = terra::toMemory(to_common_grid(rast(files[["surface"]])))
  )
}

## ---- Stage helpers (each builds its cluster, runs, stops it) ----
scenario_cost_distances <- function(pairs, res_path, sp) {
  cl <- setup_costdist_cluster(cd_n_cores)
  cl_big <- setup_costdist_cluster(cd_big_n_cores)
  cd <- NULL
  tryCatch({
    cd <- compute_cost_distances(pairs, res_path, cl, cl_big, cd_n_cores, cd_big_n_cores,
                                 cd_batch_size, cd_big_origin_cells, sp)
  }, finally = {
    try(parallel::stopCluster(if (is.null(cd)) cl else cd$cl), silent = TRUE)
    try(parallel::stopCluster(if (is.null(cd)) cl_big else cd$cl_big), silent = TRUE)
  })
  cd$costs
}

scenario_candidate_pairs <- function(dispersal_km, kern, sp) {
  cl <- parallel::makeCluster(pairs_n_cores)
  on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, library(sf))
  parallel::clusterExport(cl, "pa_sf")
  find_candidate_pairs(pa_sf, site_ids, dispersal_km, kern, cl, pairs_n_cores, sp)
}

scenario_lcps <- function(edges, res_path, edge_ids, sp) {
  cl <- setup_lcp_cluster(lcp_n_cores)
  lcp <- NULL
  tryCatch({
    lcp <- trace_lcp_corridors(edges, res_path, cl, lcp_n_cores, edge_ids = edge_ids,
                               batch_size = lcp_batch_size, sp = sp)
  }, finally = try(parallel::stopCluster(if (is.null(lcp)) cl else lcp$cl), silent = TRUE))
  lcp$paths
}

scenario_density <- function(corridors, dispersal_km, out_path, sp) {
  cl <- setup_density_cluster(dens_n_cores)
  dens <- NULL
  tryCatch({
    dens <- corridor_density(corridors, dispersal_km, cl, dens_n_cores, out_path,
                             batch_size = edge_batch_size, sp = sp)
  }, finally = try(parallel::stopCluster(if (is.null(dens)) cl else dens$cl), silent = TRUE))
  dens$path
}

## ---- One scenario: rerun the affected stages, compare, write 3 files ----
run_scenario <- function(sp, sc, prod, out_dir) {
  t0 <- Sys.time()
  tmp <- file.path(scratch_root, paste0("sa_", sp, "_", sc$scenario))
  dir.create(tmp, showWarnings = FALSE)
  on.exit({ gc(); unlink(tmp, recursive = TRUE); terra::tmpFiles(remove = TRUE) }, add = TRUE)

  is_disp <- sc$parameter == "dispersal_mult"
  dispersal_km <- prod$dispersal_km * if (is_disp) sc$value else 1
  kern <- species_kernel(dispersal_km)

  ## 02/03: resistance raster and pair costs
  if (is_disp) {
    res_path <- prod$res_path
    pair_costs <- prod$pair_costs
    n_pairs_costed <- 0
    if (kern$search_radius > prod$kernel$search_radius) {
      pairs <- scenario_candidate_pairs(dispersal_km, kern, sp)
      new_pairs <- if (is.null(pairs)) pairs[0, ] else anti_join(pairs, pair_costs, by = c("from", "to"))
      n_pairs_costed <- NROW(new_pairs)
      if (n_pairs_costed > 0) pair_costs <- bind_rows(pair_costs, scenario_cost_distances(new_pairs, res_path, sp))
    }
  } else {
    ## contrast: resistance = 2^(m * rank) = (production resistance)^m
    res_path <- file.path(tmp, "resistance.tif")
    writeRaster(rast(prod$res_path)^sc$value, res_path, overwrite = TRUE)
    pair_costs <- scenario_cost_distances(prod$pairs, res_path, sp)
    n_pairs_costed <- nrow(prod$pairs)
  }

  ## 04: MST at this scenario's cutoff
  mst <- build_mst(pair_costs, site_ids, kern$alpha, kern$p_threshold, sp, pa_area_km2)
  edges <- mst$edges
  write_csv(edges %>% dplyr::select(from, to, cost_dist = weight, euclidean_dist, dispersal_probability),
            file.path(out_dir, paste0(sc$scenario, "_mst_edges.csv")))

  ## 05: corridors (dispersal scenarios reuse production's LCP of any edge it traced)
  keys <- mst_edge_keys(edges)
  paths <- if (is_disp) unname(prod$lcp_by_key[keys]) else vector("list", nrow(edges))
  to_trace <- which(vapply(paths, is.null, logical(1)))
  if (length(to_trace) > 0) paths[to_trace] <- scenario_lcps(edges, res_path, to_trace, sp)[to_trace]
  traced <- !vapply(paths, is.null, logical(1))

  ## 06/07: density -> SDM overlap (or density itself if no SDM), on the 1 km grid
  corridors <- NULL
  surface_1km <- NULL
  if (any(traced)) {
    geoms <- do.call(c, lapply(paths[traced], function(p) st_set_crs(st_geometry(p)[1], NA)))
    corridors <- st_sf(from = edges$from[traced], to = edges$to[traced],
                       geometry = st_set_crs(geoms, st_crs(pa_patches)))
    density_path <- scenario_density(corridors, dispersal_km, file.path(tmp, "density.tif"), sp)
    if (!is.null(density_path)) {
      density_r <- rast(density_path)
      surface <- if (prod$has_sdm) sdm_overlap(density_r, prod$sdm_path, quantile_cutoff) else density_r
      surface_1km <- to_common_grid(surface)
      writeRaster(surface_1km, file.path(out_dir, paste0(sc$scenario, "_", prod$surface_type, "_1km.tif")),
                  overwrite = TRUE)
    }
  }

  ## compare with production, write the result row last
  surface_cmp <- if (is.null(surface_1km)) tibble(spearman = NA_real_, schoener_d = NA_real_, top_jaccard = NA_real_)
                 else compare_surfaces(surface_1km, prod$surface_1km)
  row <- bind_cols(
    tibble(species = sp, scenario = sc$scenario, parameter = sc$parameter, value = sc$value,
           dispersal_km = dispersal_km, surface = prod$surface_type),
    mst$summary %>% dplyr::select(n_mst_edges, n_components, total_mst_euclidean_dist_m,
                                  eca_frac, largest_component_area_frac),
    compare_msts(edges, prod$edges, site_ids) %>%
      dplyr::select(edge_jaccard, baseline_edge_retention = ref_edge_retention, degree_spearman),
    compare_corridors(if (is.null(corridors)) NULL else corridor_set(corridors), prod$corr),
    surface_cmp %>% rename_with(~ paste0("surface_", .x)),
    tibble(n_pairs_costed = n_pairs_costed, n_edges_traced = length(to_trace),
           n_edges_untraced = sum(!traced),
           minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  )
  write_csv(row, file.path(out_dir, paste0(sc$scenario, "_result.csv")))
  row
}

## ============================================================
## Run
## ============================================================
sa_species <- if (is.null(species_subset)) kernel_params$species else intersect(kernel_params$species, species_subset)

## Every (species, scenario) pair to run, numbered in a fixed order. Locally
## they all run in sequence. In a SLURM array job (Scripts/Analysis/fir_sensitivity.sh)
## each task runs only the row whose number = SLURM_ARRAY_TASK_ID, so submit
## with --array=1-<number of rows printed below> (same toggles as here).
tasks <- bind_rows(lapply(sa_species, function(sp) {
  sc <- make_scenarios()
  if (include_identity_check) {
    sc <- bind_rows(tibble(scenario = "identity_check", parameter = "contrast_mult", value = 1), sc)
  }
  sc %>% mutate(species = sp, .before = 1)
}))
if (!is.null(scenario_subset)) tasks <- tasks %>% filter(scenario %in% scenario_subset)
tasks <- tasks %>% mutate(task = row_number(), .before = 1)
message(nrow(tasks), " scenario(s) -- as a Fir array job: --array=1-", nrow(tasks))
print(tasks, n = Inf)

task_id <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA)))
if (!is.na(task_id)) {
  if (task_id < 1 || task_id > nrow(tasks)) {
    stop("SLURM_ARRAY_TASK_ID = ", task_id, " is outside 1-", nrow(tasks), " (the task list above).")
  }
  tasks <- tasks[task_id, ]
}

prod <- NULL
for (j in seq_len(nrow(tasks))) {
  sc <- tasks[j, ]
  sp <- sc$species
  out_dir <- file.path(sensitivity_dir, sp)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  if (file.exists(file.path(out_dir, paste0(sc$scenario, "_result.csv"))) && !overwrite) {
    message(sp, " ", sc$scenario, ": already done -- skipping (set overwrite <- TRUE to redo).")
    next
  }
  if (is.null(prod) || prod$species != sp) prod <- c(list(species = sp), load_production(sp))

  message("==== task ", sc$task, ": ", sp, " ", sc$scenario)
  row <- run_scenario(sp, sc, prod, out_dir)
  message(sprintf("==== %s %s done in %.1f min -- MST edges shared with production: %.0f%%, %s similarity (Schoener's D): %.2f",
                  sp, sc$scenario, row$minutes, 100 * row$edge_jaccard, row$surface, row$surface_schoener_d))
}
