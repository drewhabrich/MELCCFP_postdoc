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
# SIGMA IS SPECIES-SPECIFIC: sigma = corridor_sigma(dispersal_km), defined
# in setup_script.R next to the dispersal-kernel helpers. It's the MEAN of
# the species' negative-exponential dispersal kernel (median_dispersal /
# ln(2)) rather than the median itself, so the smoothed surface reflects a
# species' broader potential-movement footprint around its traced least-
# cost paths, not just the paths themselves. Ranges from a 60 m floor
# (near-sessile species) to ~144 km (species with dispersal_km >= 100,
# which are all treated as 100 km, matching `sparsify_dispersal_km` in
# script 01). See corridor_sigma() for the exact formula.
#
# Steps per edge (an MST edge's traced LCP corridor from script 05):
#   1. terra::rasterizeGeom(fun="length") -- total corridor length
#      falling within each cell, in km
#   2. divide by cell area (km^2) to get a RAW, unsmoothed density
#   3. terra::focal() with a Gaussian weight matrix (terra::focalMat,
#      sigma = this species' sigma) -- a normalized (sum-to-1) kernel, so
#      this is a proper density-preserving smoothing convolution, not
#      just an arbitrary blur
# Every edge's result is then combined into one {species}_corridor_density.tif.
# (All of this is corridor_density() / edge_density_worker() in
# pipeline_functions.R; cells_per_sigma is set in setup_script.R.)
#
# SCALE: rasterizing/smoothing the FULL modelling-extent grid at its native
# 30 m (~1.1 billion cells) is infeasible even per edge once sigma runs
# into the tens of km (a species like MAAM, sigma ~17 km, would need a
# ~3,460-cell-wide 1D kernel over a ~14M-cell crop, per edge, x ~500+
# edges). So each species is smoothed on its OWN coarser grid, scaled to
# its sigma (roughly sigma / 10, snapped to a multiple of 30 m) -- fine
# detail a species could never resolve at that sigma anyway is dropped
# before the convolution rather than computed and thrown away. The grid
# is never coarser than common_res_m (1 km, the SDMs' resolution), so
# wide-ranging species (sigma > ~10 km) are all smoothed on the common
# 1 km grid and script 07's overlay keeps the SDM's full detail -- e.g.
# ASFL's sigma ~144 km would otherwise give ~14 km cells. Measured cost
# at 1 km: ~1.4 s for ASFL's largest edge crop. A corridor network that only
# occupies a small corner of the full extent still wastes time/disk if
# processed at the whole-species scale, though -- a species' MST edges can
# be scattered across most of the study area even though each individual
# edge is short, so the union bounding box of all of a species' edges is
# still close to the full extent. So, like script 05, this crops per EDGE:
# each edge's own bounding box on the species' coarse grid, padded by
# 3 x sigma (the padding matters specifically for the Gaussian focal
# smoothing, which needs real neighboring cells out to a few sigma to
# avoid artificially low smoothed values right at the crop edge -- 3 sigma
# is the standard "effectively zero beyond this" cutoff for a Gaussian).
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
#   - A SpatRaster can't be exported to a worker any more than a loaded
#     resistance raster could in 03/05 -- it wraps a C++ pointer, not
#     serializable data. Each species' coarse grid is described instead by
#     plain extent/resolution/CRS values (ordinary numbers and a string),
#     exported fresh whenever the grid resolution changes (i.e. every
#     species, since sigma -- and so grid resolution -- is species-
#     specific), and each worker rebuilds the SpatRaster from those inside
#     edge_density_worker() itself.
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
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to test/compare just those species instead.
species_subset <- c("MAAM", "ASFL", "DOOR", "DRPI", "AQCH")  # e.g. NULL for all species

overwrite <- TRUE  # set TRUE to force recomputation of species already done

## reads script 05's corridor_dir, writes density_dir (both set in setup_script.R)
if (!dir.exists(density_dir)) dir.create(density_dir)

specieslist <- readxl::read_excel(here(tabular_dir, "specieslist.xlsx"))
if (!is.null(species_subset)) {
  specieslist <- specieslist %>% filter(species %in% species_subset)
}

## How many edges to hand to the cluster per parLapplyLB() call.
edge_batch_size <- 50

n_cores <- max(1, parallel::detectCores() - 2)  # tune down if you want to keep using your machine for other things

message("Setting up a ", n_cores, "-worker cluster for parallel corridor density calculation.")
cl <- setup_density_cluster(n_cores)

## Collects each species' smoothed density raster (freshly computed, or
## re-read from disk if skipped as already-done) for the combined
## multi-species map built after the loop.
density_rasters <- list()

tryCatch({

  for (i in seq_len(nrow(specieslist))) {
    sp <- specieslist$species[i]
    out_path <- file.path(density_dir, paste0(sp, "_corridor_density.tif"))

    if (file.exists(out_path) && !overwrite) {
      message(sp, ": ", out_path, " already exists -- skipping (set overwrite <- TRUE to redo).")
      density_rasters[[sp]] <- rast(out_path)
      next
    }

    corridors_path <- file.path(corridor_dir, paste0(sp, "_lcp_corridors.gpkg"))
    if (!file.exists(corridors_path)) {
      warning(sp, ": missing ", corridors_path, " -- run script 05 first. Skipping.")
      next
    }

    corridors <- st_read(corridors_path, quiet = TRUE)

    ## Rasterize + smooth every edge in parallel, mosaic-sum, write out_path
    ## (keep the returned cluster -- it's rebuilt if a worker died).
    dens <- corridor_density(corridors, specieslist$dispersal_km[i], cl, n_cores, out_path,
                             batch_size = edge_batch_size, sp = sp)
    cl <- dens$cl
    if (is.null(dens$path)) next

    ## Re-read from out_path rather than keeping the mosaic() result --
    ## terra rasters can stay lazily linked to their source files, and the
    ## temp edge rasters were deleted after the mosaic, so hanging onto the
    ## original object risks a dangling reference later (the cross-species
    ## step below).
    density_rasters[[sp]] <- rast(out_path)
  }

}, finally = {
  parallel::stopCluster(cl)
  message("Cluster stopped.")
})

## Combined multi-species map: each species was smoothed on its OWN
## dispersal-scaled grid (see SIGMA/SCALE notes above), so resolutions
## differ across species and have to be reconciled onto one common grid
## before combining -- to_common_grid() (1 km = common_res_m, matching the
## SDM's native resolution; see script 07/00) averages species whose own
## grid is finer than 1 km (species capped at 1 km are already on it).
## Each species is then normalized to
## [0,1] by its own max (so no single species' raw density scale
## dominates) before summing via mosaic(), since species still differ in
## extent (each cropped to just its own corridors).
density_rasters <- density_rasters[!vapply(density_rasters, is.null, logical(1))]

if (length(density_rasters) > 0) {
  normalized <- lapply(density_rasters, function(r) {
    r_common <- to_common_grid(r)

    r_max <- as.numeric(global(r_common, "max", na.rm = TRUE)[1, 1])
    if (is.na(r_max) || r_max == 0) return(r_common)  # nothing to normalize against (all-zero raster)
    r_common / r_max
  })

  multi_species_density <- mosaic(sprc(normalized), fun = "sum")

  multi_out_path <- here(density_dir, "all_species_corridor_density.tif")
  writeRaster(multi_species_density, multi_out_path, overwrite = TRUE)
  message("-> saved normalized multi-species corridor density (species: ",
          paste(names(density_rasters), collapse = ", "), ") to ", multi_out_path)
} else {
  warning("No species' corridor density available -- multi-species map not created.")
}
