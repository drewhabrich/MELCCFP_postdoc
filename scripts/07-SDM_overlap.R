## ------------------------------- ##
## Script name: 07-SDM_overlap
##
## Purpose of script: Overlay the kernel-smoothed corridor density surface
# (script 06) with an externally-produced species distribution model (SDM)
# raster, to get a species-specific connectivity map -- how well does the
# traced corridor network line up with independently-modeled suitable
# habitat? Only runs for species with an SDM raster available
# (species_metadata.csv's sdm_raster_path, set by setup_species_from_xlsx.R
# by matching each species' scientific name against data/drive_data/SDMs/
# -- species without a match there are skipped, since SDMs aren't expected
# for every species).
#
# Adapted from a worked example script (Suitable_habitat_overlap.R):
#   1. drop the bottom `quantile_cutoff` of corridor density values as
#      noise (set to NA), then rescale the remainder to [0, 1]
#   2. reproject/crop/mask the SDM raster onto the corridor raster's grid
#   3. multiply corridor x SDM suitability, rescale the product to [0, 1]
#
## Author: Andrew Habrich
##
## Notes ------------------------- ##
# SCALE: the original design (from the worked example this was adapted
# from) started by extend()-ing the corridor raster out to the FULL study
# area, filling the new cells with the corridor raster's own minimum. At
# this project's scale that reintroduces the exact full-modelling-extent
# cost scripts 03/05/06 already had to be fixed for --
# {species}_corridor_density.tif is already cropped to that species' own
# corridor network (script 06's per-edge design), and blowing it back up
# to province scale just to pad mostly-empty space with a near-constant
# background value isn't worth the cost, especially since nothing
# downstream of this script needs a full-extent raster. So this version
# skips extend()/NA-fill entirely and works on the already-cropped extent
# throughout. The "how does this sit in the bigger picture" visual context
# isn't lost: the QA PNG still draws the full study area boundary as a
# vector layer, and ggplot's default panel extent expands to fit every
# layer -- so the saved PNG shows the small raster patch positioned within
# the full study area outline regardless of the raster's own extent. Only
# the raster DATA stays small; the figure's context doesn't change.
# ============================================================
source(file.path("scripts", "setup_script.R"))
library(ggplot2)
library(tidyterra)

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to test/compare just those species instead.
species_subset <- c("MAAM", "DOOR")  # e.g. NULL for all species

overwrite <- FALSE  # set TRUE to force recomputation of species already done

## Fraction of the corridor density raster's lowest values to treat as
## noise (set to NA before rescaling), e.g. 0.1 = drop the bottom 10%.
## Adapted from the worked example this was based on, not a validated
## threshold -- check whether it's dropping meaningful low-density
## corridor or genuinely negligible noise for your data before trusting it.
quantile_cutoff <- 0.1

dir_corridors <- here(interm_dir, "corridors")
dir_habitat_overlap <- here(output_dir, "mst_SDM_overlap")
if (!dir.exists(dir_habitat_overlap)) dir.create(dir_habitat_overlap)

rescale_raster <- function(r) {
  mm <- minmax(r)
  (r - mm[1, 1]) / (mm[2, 1] - mm[1, 1])
}

species_meta <- read_csv(here(data_dir, "species_metadata.csv"), show_col_types = FALSE)
if (!is.null(species_subset)) {
  species_meta <- species_meta %>% filter(species %in% species_subset)
}

# Load the study extent and PA nodes for the loop
study_area <- vect(here(data_dir, "study_area.shp"))
nodes_v <- vect(here(output_dir, "pa_ctroidnodes.shp"))

# Loop over each species to get the overlap with the LC corridor and the SDM model
for (i in seq_len(nrow(species_meta))) {
  sp <- species_meta$species[i]
  out_path <- file.path(dir_habitat_overlap, paste0(sp, "_habitat_overlap.tif"))

  if (file.exists(out_path) && !overwrite) {
    message(sp, ": ", out_path, " already exists -- skipping (set overwrite <- TRUE to redo).")
    next
  }

  sdm_path <- species_meta$sdm_raster_path[i]
  if (is.na(sdm_path) || !file.exists(sdm_path)) {
    message(sp, ": no SDM raster available -- skipping.")
    next
  }

  corridor_path <- file.path(dir_corridors, paste0(sp, "_corridor_density.tif"))
  if (!file.exists(corridor_path)) {
    warning(sp, ": missing ", corridor_path, " -- run script 06 first. Skipping.")
    next
  }

  message("Computing habitat overlap for: ", sp)
  corridor <- rast(corridor_path)

  # drop the bottom quantile_cutoff of values as noise, then rescale (remove bottom 10%)
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
  overlap <- corridor_q * sdm
  overlap <- rescale_raster(overlap)

  writeRaster(overlap, out_path, overwrite = TRUE)
  message("  -> saved ", out_path)
}

# # quick QA figure -- full study area boundary + PA nodes for context,
# # even though the raster data itself only covers this species' corridor extent
grad <- hypso.colors(15, "dem_screen")
overlap <- rast(here(dir_habitat_overlap, "MAAM_habitat_overlap.tif"))

p <- ggplot() +
  geom_spatraster(data = overlap, maxcell = 500000) +
  geom_spatvector(data = study_area, fill = NA, colour = "gray20", linewidth = 0.8) +
  geom_spatvector(data = nodes_v, fill = "lightblue", colour = "grey20", alpha = 0.3) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours = grad, na.value = NA, name = "Connectivity") +
  theme_bw()

plot(p)

#save to file
png_path <- file.path(fig_dir, paste0(sp, "_habitat_overlap.png"))
ggsave(plot = p, width = 10, height = 8, filename = png_path)
message("  -> saved ", png_path)