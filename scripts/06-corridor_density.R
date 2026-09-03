## ------------------------------- ##
## Script name: 06-corridor_density
##
## Purpose of script: Build a kernel-smoothed corridor density surface per
# species -- each pixel's value is (after smoothing) roughly "km of
# corridor per km^2", the same concept as ArcGIS's "Line Density" tool --
# plus one combined multi-species "hotspot" map (each species normalized
# to [0,1] by its own max before summing, so a species with a denser/more
# numerous corridor network doesn't dominate the total).
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
# Steps per edge (an MST edge's traced LCP corridor from script 05):
#   1. terra::rasterizeGeom(fun="length") -- total corridor length
#      falling within each cell, in km
#   2. divide by cell area (km^2) to get a RAW, unsmoothed density
#   3. terra::focal() with a Gaussian weight matrix (terra::focalMat,
#      sigma = corridor_density_sigma) -- a normalized (sum-to-1) kernel,
#      so this is a proper density-preserving smoothing convolution, not
#      just an arbitrary blur
# Every edge's result is then combined into one {species}_corridor_density.tif.
#
# SCALE: rasterizing/smoothing the FULL modelling-extent grid (~1.1
# billion cells) for a corridor network that only occupies a small corner
# of it wastes the same order of time/disk that script 03's cost-distance
# and script 05's LCP tracing already had to be fixed for. Cropping to a
# whole SPECIES' corridor network at once isn't enough, though -- a
# species' MST edges can be scattered across most of the study area even
# though each individual edge is short (it connects two geographically
# close PAs), so the union bounding box of all of a species' edges is
# still close to the full extent. So, like script 05, this crops per
# EDGE: each edge's own bounding box, padded by 3 x corridor_density_sigma
# (the padding matters specifically for the Gaussian focal smoothing,
# which needs real neighboring cells out to a few sigma to avoid
# artificially low smoothed values right at the crop edge -- 3 sigma is
# the standard "effectively zero beyond this" cutoff for a Gaussian).
# Each edge's small raster is rasterized AND smoothed on its own crop,
# then all of a species' edges are combined with terra::mosaic(fun="sum")
# -- valid because Gaussian smoothing is linear (smoothing each edge's
# contribution separately and summing gives the same result as summing
# first and smoothing once), and because the only operation touching the
# full per-species extent is that final mosaic, a cheap overlay/combine
# rather than a per-cell convolution.
#
# The Gaussian kernel itself is also applied as two 1D passes (row, then
# column) rather than one 2D matrix -- mathematically identical for a
# separable kernel like this one, at a fraction of the cost (~2k vs ~k^2
# per cell, where k is the kernel's cell-radius).
#
# PARALLELIZATION: edges are traced in parallel across CPU cores, the same
# PSOCK-cluster pattern used in scripts 03 and 05 (base R's `parallel`
# package -- ships with R, no install needed), since each edge's
# rasterize+smooth is independent of every other. Two terra-specific
# wrinkles this raises, beyond what 03/05 already had to handle:
#   - `template_raster` (a SpatRaster) can't be exported to a worker any
#     more than a loaded resistance raster could in 03/05 -- it wraps a
#     C++ pointer, not serializable data. Its plain extent/resolution/CRS
#     values (ordinary numbers and a string) ARE exportable, so each
#     worker rebuilds an identical `template_raster` from those once, at
#     cluster setup.
#   - A `SpatRaster` (or `SpatVector`) can't be *returned* from a worker
#     either, for the same reason -- unlike script 05's edges, which
#     return `sf` objects (sf stores geometry as plain serializable data,
#     not a live C++ pointer, so that direction was never a problem there).
#     So each worker writes its edge's small smoothed raster to a temp
#     .tif and returns the file path (a plain string); the master reads
#     each path back with `rast()` and deletes the temp file once read.
#     Similarly, an edge's LINESTRING geometry is passed to its worker as
#     a plain coordinate matrix (`terra::geom(v)[, c("x","y")]`), not a
#     SpatVector, and rebuilt into one with `vect(..., type = "lines")` on
#     the worker side.
#
# MULTI-SPECIES SUM: species are cropped to different extents (each to its
# own corridors), so combining them uses terra::mosaic(..., fun = "sum")
# rather than a plain `+`/`sum()` across a raster stack, which requires
# identical extents. Each species' density is normalized to [0,1] (divided
# by its own max) first, so the combined map reflects relative corridor
# importance per species, not raw density scale (which varies a lot with
# how many MST edges/PAs a species has).
# ============================================================
source(file.path("scripts", "setup_script.R"))

dir_corridors <- here(interm_dir, "corridors")
if (!dir.exists(dir_corridors)) dir.create(dir_corridors)

specieslist <- readxl::read_excel(here(data_dir, "specieslist.xlsx"))
if (!is.null(species_subset)) {
  specieslist <- specieslist %>% filter(species %in% species_subset)
}

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to test/compare just those species instead.
species_subset <- c("MAAM", "DOOR")  # e.g. NULL for all species

overwrite <- TRUE  # set TRUE to force recomputation of species already done

## Gaussian smoothing bandwidth, map units (metres). No principled default
## -- a visualization choice trading off how tightly the density surface
## hugs the traced corridors vs. how much nearby corridors blur together.
## Inspect the output and tune for this project's PA spacing.

#-#-#-# SHOULD THIS BE A SPECIES SPECIFIC PARAMETER? #-#-#-#
corridor_density_sigma <- 5000

## How many edges to hand to the cluster per parLapplyLB() call.
edge_batch_size <- 50

n_cores <- max(1, parallel::detectCores() - 10)  # tune down if you want to keep using your machine for other things

## Depend only on cell resolution/sigma, not on any per-species/per-edge
## crop -- computed once rather than per species like the original design did.
cell_area_km2 <- (xres(template_raster) * yres(template_raster)) / 1e6
crop_buffer_m <- 3 * corridor_density_sigma

## 1D Gaussian profile, for separable smoothing (two 1D passes instead of
## one 2D matrix -- see SCALE note above). Extracted from terra's own 2D
## matrix (its center row, renormalized to sum to 1) rather than
## re-derived by hand, so it reproduces exactly what the 2D kernel would
## have -- a separable Gaussian's 2D matrix is the outer product of this
## profile with itself.
gauss_2d <- focalMat(template_raster, corridor_density_sigma, "Gauss")
gauss_1d <- gauss_2d[(nrow(gauss_2d) + 1) / 2, ]
gauss_1d <- gauss_1d / sum(gauss_1d)

## Plain (serializable) description of template_raster's grid, for workers
## to rebuild it from -- see PARALLELIZATION note above.
template_ext_vec <- as.vector(ext(template_raster))
template_res <- res(template_raster)
template_crs_wkt <- crs(template_raster)

## One MST edge's density contribution: rebuild its geometry from a plain
## coordinate matrix, rasterize + smooth on a small crop around it, write
## the result to a temp .tif, and return the path -- a SpatRaster can't be
## returned directly from a parallel worker (see PARALLELIZATION note).
## References `edge_coords_list` as a free variable -- exported to the
## cluster each time it changes (i.e. each species).
edge_density_worker <- function(e) {
  edge_v <- vect(edge_coords_list[[e]], type = "lines", crs = template_crs_wkt)

  crop_ext <- ext(edge_v) + crop_buffer_m
  grid_crop <- crop(template_raster, crop_ext)

  length_per_cell_km <- rasterizeGeom(edge_v, grid_crop, fun = "length", unit = "km")
  density_raw <- length_per_cell_km / cell_area_km2  # km corridor per km^2, unsmoothed

  ## Separable Gaussian smoothing: one pass along rows, one along columns.
  density_smoothed <- focal(density_raw, w = matrix(gauss_1d, nrow = 1), fun = "sum", na.rm = TRUE)
  density_smoothed <- focal(density_smoothed, w = matrix(gauss_1d, ncol = 1), fun = "sum", na.rm = TRUE)

  out_path <- tempfile(fileext = ".tif")
  writeRaster(density_smoothed, out_path, overwrite = TRUE)
  out_path
}

## Builds a fresh cluster: each worker rebuilds template_raster from plain
## values (see PARALLELIZATION note), and gets everything species-invariant
## exported once. Used for the initial cluster and again if a worker dies
## mid-run.
setup_cluster <- function() {
  new_cl <- parallel::makeCluster(n_cores)
  parallel::clusterExport(new_cl, varlist = c(
    "template_ext_vec", "template_res", "template_crs_wkt",
    "cell_area_km2", "crop_buffer_m", "gauss_1d", "edge_density_worker"
  ))
  parallel::clusterEvalQ(new_cl, {
    library(terra)
    template_raster <- rast(ext(template_ext_vec), resolution = template_res, crs = template_crs_wkt)
  })
  new_cl
}

message("Setting up a ", n_cores, "-worker cluster for parallel corridor density calculation.")
cl <- setup_cluster()

## Collects each species' smoothed density raster (freshly computed, or
## re-read from disk if skipped as already-done) for the combined
## multi-species map built after the loop.
density_rasters <- list()

tryCatch({

  for (i in seq_len(nrow(specieslist))) {
    sp <- specieslist$species[i]
    out_path <- file.path(dir_corridors, paste0(sp, "_corridor_density.tif"))

    if (file.exists(out_path) && !overwrite) {
      message(sp, ": ", out_path, " already exists -- skipping (set overwrite <- TRUE to redo).")
      density_rasters[[sp]] <- rast(out_path)
      next
    }

    corridors_path <- file.path(dir_corridors, paste0(sp, "_lcp_corridors.gpkg"))
    if (!file.exists(corridors_path)) {
      warning(sp, ": missing ", corridors_path, " -- run script 05 first. Skipping.")
      next
    }

    corridors <- st_read(corridors_path, quiet = TRUE)
    corridors_v <- vect(corridors)
    n_edges <- nrow(corridors_v)

    ## Plain coordinate matrices, one per edge -- what actually gets
    ## exported to the cluster (SpatVector itself can't be).
    edge_coords_list <- lapply(seq_len(n_edges), function(e) {
      terra::geom(corridors_v[e, ])[, c("x", "y")]
    })
    parallel::clusterExport(cl, varlist = "edge_coords_list")

    message("Computing corridor density for: ", sp, " (", n_edges, " edge(s), ", n_cores, " worker(s))")

    edge_tif_paths <- vector("list", n_edges)
    batches <- split(seq_len(n_edges), ceiling(seq_len(n_edges) / edge_batch_size))
    n_done <- 0

    for (b in batches) {
      ## If a worker dies mid-batch, parLapplyLB() throws rather than
      ## returning a result. Rebuild the cluster and mark this batch's
      ## edges as failed rather than letting that abort the whole run.
      batch_results <- tryCatch(
        parallel::parLapplyLB(cl, b, edge_density_worker),
        error = function(e) {
          warning("  ", sp, ": cluster worker failed on this batch (", conditionMessage(e),
                  ") -- rebuilding the cluster and marking these ", length(b),
                  " edge(s) as failed, then continuing.")
          cl <<- setup_cluster()
          parallel::clusterExport(cl, varlist = "edge_coords_list")
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
      next
    }

    ## Combine every edge's smoothed contribution into one species-level
    ## density surface. Valid because Gaussian smoothing is linear: smoothing
    ## each edge separately and summing equals summing first and smoothing
    ## once (to the same small truncation tolerance already accepted by the
    ## 3-sigma crop padding).
    edge_density_rasters <- lapply(edge_tif_paths, rast)
    density_smoothed <- mosaic(sprc(edge_density_rasters), fun = "sum")

    writeRaster(density_smoothed, out_path, overwrite = TRUE)
    unlink(unlist(edge_tif_paths))  # done reading these back -- clean up the temp per-edge rasters
    message("  -> saved ", out_path, " (", length(edge_density_rasters), " of ", n_edges, " edge(s) combined",
            if (n_failed > 0) paste0(", ", n_failed, " skipped") else "", ")")

    ## Re-read from out_path rather than keeping the mosaic() result --
    ## terra rasters can stay lazily linked to their source files, and the
    ## temp edge rasters just got deleted above, so hanging onto the
    ## original object risks a dangling reference later (the cross-species
    ## step below).
    density_rasters[[sp]] <- rast(out_path)
  }

}, finally = {
  parallel::stopCluster(cl)
  message("Cluster stopped.")
})

## Combined multi-species map: normalize each species to [0,1] by its own
## max (so no single species' raw density scale dominates), then sum via
## mosaic() since species are cropped to different extents.
density_rasters <- density_rasters[!vapply(density_rasters, is.null, logical(1))]

if (length(density_rasters) > 0) {
  normalized <- lapply(density_rasters, function(r) {
    r_max <- as.numeric(global(r, "max", na.rm = TRUE)[1, 1])
    if (is.na(r_max) || r_max == 0) return(r)  # nothing to normalize against (all-zero raster)
    r / r_max
  })

  multi_species_density <- mosaic(sprc(normalized), fun = "sum")

  multi_out_path <- here(dir_corridors, "all_species_corridor_density.tif")
  writeRaster(multi_species_density, multi_out_path, overwrite = TRUE)
  message("-> saved normalized multi-species corridor density (species: ",
          paste(names(density_rasters), collapse = ", "), ") to ", multi_out_path)
} else {
  warning("No species' corridor density available -- multi-species map not created.")
}
