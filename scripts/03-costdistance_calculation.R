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
# If a destination comes back with cost_dist = NA (unreachable within the
# crop -- the buffer was too tight to route around some barrier), that
# origin is retried once with a doubled buffer before accepting the NA.
#
## Author: Andrew Habrich
##
## Notes ------------------------- ##
source(file.path("scripts", "setup_script.R"))
library(dplyr)

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to quickly test/compare just those species instead (e.g. before
## committing to a full multi-hour run across all 10).
species_subset <- c("MAAM", "DOOR")  # e.g. NULL for all species

overwrite <- FALSE  # set TRUE to force recomputation of species already done

## How far beyond an origin's candidate-destination bounding box to crop
## the resistance raster, as a fraction of that origin's farthest candidate
## distance (plus a minimum, so origins with few/close destinations still
## get a sane amount of room). Doubled once automatically if any
## destination comes back unreachable within the crop.
cd_buffer_frac <- 0.25
cd_buffer_min_m <- 1000

## How many origins to hand to the cluster per parLapplyLB() call.
cd_batch_size <- 50

n_cores <- max(1, parallel::detectCores() - 12)  # tune down if you want to keep using your machine for other things

## 1. Load required files ---------------------------------------
# nodes: sf POINT object (or polygon centroids) of protected areas,
#        must have a unique `site_id` column
nodes <- read_sf(here(output_dir, "pa_ctroidnodes.shp"))
node_coords <- st_coordinates(nodes)
rownames(node_coords) <- nodes$site_id

specieslist <- readxl::read_excel(here(data_dir, "specieslist.xlsx"))
if (!is.null(species_subset)) {
  specieslist <- specieslist |> filter(species %in% species_subset)
}

## Compute cost-distance from one origin to all its candidate destinations,
## on a crop of the resistance raster sized to just that origin's search
## area. References `candidate_pairs` / `resistance_path` as free variables
## -- exported to the cluster (setup_cluster(), below) each time they
## change (i.e. each species).
cost_distance_for_origin <- function(origin_id, buffer_mult = 1) {
  origin_xy <- node_coords[as.character(origin_id), , drop = FALSE]

  dest_pairs <- candidate_pairs  |>  filter(from == origin_id)
  dest_ids <- dest_pairs |> pull(to)
  dest_xy <- node_coords[as.character(dest_ids), , drop = FALSE]

  all_xy <- rbind(origin_xy, dest_xy)
  buffer_dist <- max(max(dest_pairs$euclidean_dist) * cd_buffer_frac, cd_buffer_min_m) * buffer_mult

  crop_ext <- ext(
    min(all_xy[, 1]) - buffer_dist, max(all_xy[, 1]) + buffer_dist,
    min(all_xy[, 2]) - buffer_dist, max(all_xy[, 2]) + buffer_dist
  )

  friction <- crop(rast(resistance_path), crop_ext)
  origin_cell <- cellFromXY(friction, origin_xy)
  friction[origin_cell] <- sentinel_value

  cost_surface <- costDist(friction, target = sentinel_value, scale = 1)
  cost_vals <- terra::extract(cost_surface, dest_xy)[[1]]

  terra::tmpFiles(remove = TRUE)

  if (any(is.na(cost_vals)) && buffer_mult == 1) {
    return(cost_distance_for_origin(origin_id, buffer_mult = 2))
  }

  tibble(from = origin_id, to = dest_ids,
         euclidean_dist = dest_pairs$euclidean_dist,
         euclidean_kernel_prob = dest_pairs$euclidean_kernel_prob,
         cost_dist = cost_vals)
}

## Builds a fresh cluster with everything that doesn't change across
## species already exported. Used for the initial cluster and again if a
## worker dies mid-run.
setup_cluster <- function() {
  new_cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(new_cl, {
    library(terra); library(sf); library(dplyr); library(tibble)
  })
  parallel::clusterExport(new_cl, varlist = c(
    "node_coords", "sentinel_value", "cd_buffer_frac", "cd_buffer_min_m",
    "cost_distance_for_origin"
  ))
  new_cl
}

message("Setting up a ", n_cores, "-worker cluster for parallel cost-distance calculation.")
cl <- setup_cluster()

## 2. Loop over every species and compute cost-distances for its candidate pairs ####
costoutput_dir <- file.path(interm_dir, "cost_distances")
if (!dir.exists(costoutput_dir)) dir.create(costoutput_dir)

tryCatch({

  for (i in seq_len(nrow(specieslist))) {
    sp <- specieslist$species[i]
    out_path <- here(costoutput_dir, paste0(sp, "_PA_costs.csv"))

    if (file.exists(out_path) && !overwrite) {
      message(sp, ": ", out_path, " already exists -- skipping (set overwrite <- TRUE to redo).")
      next
    }

    ## Read the list of candidate node pairs with their Euclidean distances
    pairs_path <- here(interm_dir, paste0(sp, "_PApairs_dkernel_cutoff.csv"))
    if (!file.exists(pairs_path)) {
      warning(sp, ": missing ", pairs_path, " -- run script 01 first. Skipping.")
      next
    }
    candidate_pairs <- read.csv(pairs_path)

    ## This species' resistance raster (from script 02)
    resistance_path <- here(output_dir, "resistance", paste0(sp, "_res.tif"))
    if (!file.exists(resistance_path)) {
      warning(sp, ": missing ", resistance_path, " -- run script 02 first. Skipping.")
      next
    }

    ## `candidate_pairs` / `resistance_path` change every species -- re-export.
    parallel::clusterExport(cl, varlist = c("candidate_pairs", "resistance_path"))

    origins <- unique(candidate_pairs$from)
    n_origins <- length(origins)
    message("Computing cost-distances for: ", sp, " (", nrow(candidate_pairs), " candidate pairs, ",
            n_origins, " origin nodes, ", n_cores, " worker(s))")

    tictoc::tic(paste(sp, "cost-distance calculation"))

    results <- vector("list", n_origins)
    batches <- split(seq_len(n_origins), ceiling(seq_len(n_origins) / cd_batch_size))
    n_done <- 0

    for (b in batches) {
      ## If a worker dies mid-batch (e.g. an unusually large search area
      ## exhausted its memory), parLapplyLB() throws rather than returning
      ## a result. Rebuild the cluster and mark this batch's origins as
      ## failed rather than letting that abort the whole run.
      batch_results <- tryCatch(
        parallel::parLapplyLB(cl, origins[b], cost_distance_for_origin),
        error = function(e) {
          warning("  ", sp, ": cluster worker failed on this batch (", conditionMessage(e),
                  ") -- rebuilding the cluster and marking these ", length(b),
                  " origin(s) as failed, then continuing.")
          cl <<- setup_cluster()
          parallel::clusterExport(cl, varlist = c("candidate_pairs", "resistance_path"))
          vector("list", length(b))
        }
      )
      results[b] <- batch_results
      n_done <- n_done + length(b)
      message("  ", sp, ": ", n_done, "/", n_origins, " origin(s) processed")
    }

    tictoc::toc()

    results <- results[!vapply(results, is.null, logical(1))]
    sp_costs <- bind_rows(results)

    write_csv(sp_costs, out_path)
    message("  -> ", nrow(sp_costs), " pair(s) saved to ", out_path)
  }

}, finally = {
  parallel::stopCluster(cl)
  message("Cluster stopped.")
})

