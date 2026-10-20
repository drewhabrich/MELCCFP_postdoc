## ------------------------------- ##
## Script name: 03-costdistance_calculation
##
## Purpose of script: Cost-distance for the surviving candidate pairs, using that
##  species' resistance raster from script 2.
#
# For each origin node, crops the resistance raster to a buffered bounding
# box around just that origin and its own candidate destinations (from
# script 01's dispersal-kernel prefilter) -- not the whole province -- then
# marks the origin with the shared sentinel value and runs terra::costDist()
# once to get accumulated cost to every one of its destinations in a single
# pass.
#
# WHY CROP: costDist() needs the whole extent it's given, and at the full
# modelling extent (~1.1 billion cells) each call produces a several-GB
# raster. terra spills rasters too big for memory to temp files on disk,
# and relies on R's garbage collector to clean them up -- which doesn't run
# immediately after each loop iteration, so hundreds of origins in a row
# can pile up hundreds of GB of temp files before anything gets reclaimed
# (this is what was exhausting disk space). Cropping to each origin's own
# search area (typically tens of km, not the whole province) cuts that
# footprint by roughly two orders of magnitude, and also makes it safe to
# parallelize across origins (below) since N workers' combined footprint
# is now N x (a few hundred MB), not N x (several GB).
# `terra::tmpFiles(remove = TRUE)` after each origin is a cheap extra
# safety net on top of that.
#
# This "tens of km per origin" assumption only holds if the candidate pairs
# are local. Species with a very large dispersal distance (dispersal_km >
# `sparsify_dispersal_km` in script 01) are sparsified there to k-NN +
# Euclidean-MST edges; without that, every origin's crop spans the whole
# study area (~3e8 cells, ~1.5 h each at 30 m).
#
# If a destination comes back with cost_dist = NA (unreachable within the
# crop -- the buffer was too tight to route around some barrier), that
# origin is retried once with a doubled buffer before accepting the NA.
#
# The per-origin worker, the crop-size estimate, the two-tier dispatch and
# the job-aware terra memory cap are compute_cost_distances() and friends
# in pipeline_functions.R; the crop buffer (cd_buffer_frac/cd_buffer_min_m)
# is set in setup_script.R.
#
## Author: Andrew Habrich
##
## Notes ------------------------- ##
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))
library(dplyr)

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to quickly test/compare just those species instead (e.g. before
## committing to a full multi-hour run across all 10).
species_subset <- c("MAAM", "ASFL", "DOOR", "DRPI", "AQCH")  # e.g. NULL for all species

overwrite <- FALSE  # set TRUE to force recomputation of species already done

## How many origins to hand to the cluster per parLapplyLB() call.
cd_batch_size <- 25

#n_cores <- max(1, parallel::detectCores() - 12)  # tune down if you want to keep using your machine for other things
# SERVER HAS 64 CORES
n_cores <- 12

## Some origins' candidate destinations are spread across an unusually wide
## area (e.g. a species with a large dispersal kernel producing
## near-fully-connected PA pairs) -- their crop extent approaches the full
## modelling extent instead of the "tens of km" this cropping scheme
## assumes, which can be several GB for a single origin and OOM-kill a
## normal worker even though the vast majority of origins crop to a few
## hundred MB. Origins whose estimated crop exceeds this many cells are
## routed to a separate, low-concurrency cluster instead, so the memory-
## hungry outliers don't need every normal worker over-provisioned for them.
## Keep this high enough that only true outliers land there: that tier has
## few workers, so routing many origins to it makes it the bottleneck (with
## sparsified ASFL pairs at k = 20, crops are ~3e7 cells median, ~1.5e8 max).
cd_big_origin_cells <- 1e8
cd_big_n_cores <- 2

## 1. Load required files ---------------------------------------
# pa_patches: sf POLYGON object of protected areas, must have a unique
# `site_id` column -- loaded as sf (not terra::vect) specifically because
# sf objects serialize cleanly to parallel workers (see
# setup_costdist_cluster()), unlike a SpatVector.
pa_patches <- read_sf(here(pa_patches_dir, "pa_patches.gpkg"))

specieslist <- readxl::read_excel(here(tabular_dir, "specieslist.xlsx"))
if (!is.null(species_subset)) {
  specieslist <- specieslist |> filter(species %in% species_subset)
}

message("Setting up a ", n_cores, "-worker cluster for parallel cost-distance calculation.")
cl <- setup_costdist_cluster(n_cores)
message("Setting up a ", cd_big_n_cores, "-worker cluster for big-crop outlier origins.")
cl_big <- setup_costdist_cluster(cd_big_n_cores)

## 2. Loop over every species and compute cost-distances for its candidate pairs ####
if (!dir.exists(costdist_dir)) dir.create(costdist_dir)

tryCatch({

  for (i in seq_len(nrow(specieslist))) {
    sp <- specieslist$species[i]
    out_path <- here(costdist_dir, paste0(sp, "_PA_costs.csv"))

    if (file.exists(out_path) && !overwrite) {
      message(sp, ": ", out_path, " already exists -- skipping (set overwrite <- TRUE to redo).")
      next
    }

    ## Read the list of candidate node pairs with their Euclidean distances
    pairs_path <- here(pa_pairs_dir, paste0(sp, "_PApairs_dkernel_cutoff.csv"))
    if (!file.exists(pairs_path)) {
      warning(sp, ": missing ", pairs_path, " -- run script 01 first. Skipping.")
      next
    }
    candidate_pairs <- read.csv(pairs_path)

    ## This species' resistance raster (from script 02)
    resistance_path <- here(resistance_dir, paste0(sp, "_res.tif"))
    if (!file.exists(resistance_path)) {
      warning(sp, ": missing ", resistance_path, " -- run script 02 first. Skipping.")
      next
    }

    ## Normal-tier origins on `cl`, big-crop outliers on `cl_big` (either
    ## may come back rebuilt if a worker died, so keep the returned ones).
    cd <- compute_cost_distances(candidate_pairs, resistance_path, cl, cl_big,
                                 n_cores, cd_big_n_cores, cd_batch_size, cd_big_origin_cells, sp)
    cl <- cd$cl
    cl_big <- cd$cl_big
    sp_costs <- cd$costs

    write_csv(sp_costs, out_path)
    message("  -> ", nrow(sp_costs), " pair(s) saved to ", out_path)

    ## Profiling: which step of cost_distance_for_origin() actually costs
    ## the time -- reading/cropping the raster from disk, rasterizing the
    ## origin mask, the costDist() flood-fill itself, or the zonal
    ## extract() -- broken down per origin (one row per attempt, so a
    ## retried origin contributes two rows: buffer_mult 1 and 2).
    sp_timing <- cd$timing
    timing_out_path <- here(costdist_dir, paste0(sp, "_costdist_profiling.csv"))
    write_csv(sp_timing, timing_out_path)
    message("  -> profiling (mean seconds/origin): crop = ", round(mean(sp_timing$crop_sec), 3),
            ", mask = ", round(mean(sp_timing$mask_sec), 3),
            ", costDist = ", round(mean(sp_timing$costdist_sec), 3),
            ", extract = ", round(mean(sp_timing$extract_sec), 3),
            " (", nrow(sp_timing), " attempt(s) -> ", timing_out_path, ")")
  }

}, finally = {
  parallel::stopCluster(cl)
  parallel::stopCluster(cl_big)
  message("Clusters stopped.")
})
