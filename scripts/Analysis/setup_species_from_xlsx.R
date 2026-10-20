## ------------------------------- ##
## Script name: setup_species_from_xlsx
##
## Purpose of script: Create the species metadata file that the other scripts refer to.
##                    Extracts the suitabilty from the sheet to parameterize the appropriate resistance scoring.
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
# ============================================================
# Stage 0: Convert Data/Tabular/specieslist.xlsx ("top10" sheet) into this
# pipeline's species_metadata.csv and per-species lookup tables for the
# "lookup" method, following the general scheme in Albert et al. 2017
# (doi: 10.1111/cobi.12943): preferred habitat gets a resistance value
# of 1, and resistance doubles with each step down in suitability
# (1, 2, 4, 8, 16, 32, ...). We extend this to layer forest age on top
# of land-cover suitability, since specieslist.xlsx has both.
#
#
# Also NOT implemented here (need data we don't have):
#     we substitute a fixed high rank ("Other" in suitability_rank_scale.csv)
#     for any land-cover class outside decid/mixed/conif/wetland_open/
#     wetland_forested/grass/shrub
#
# Run this ONCE (or whenever specieslist.xlsx changes), BEFORE stage 01.
# -- review its output before running the rest of the pipeline on it.
#
# Requires three small mapping tables you fill in once
#
#   landcover_class_lookup.csv
#     Maps YOUR landcover.tif class codes to decid/mixed/conif/
#     wetland_open/wetland_forested/grass/shrub/other. decid/mixed/conif are treated as "forest" --
#     only those cells get an age penalty applied.
#
#   forest_age_class_lookup.csv
#     Long-format ordinal-distance table: one row per (forest_age.tif class
#     code x xlsx forest_age category), giving the rank-step "distance"
#     between them (0 = exact match, +1 per ordinal step away along
#     jeune -> medium1 -> medium2 -> forêt ancienne). Class 5 (mixed/
#     inégale -- uneven-aged, containing both young and old trees) is
#     special-cased: distance 0 against jeune or forêt ancienne, 1 against
#     mature. A species' penalty for each raster class is the minimum
#     distance across its (possibly multiple) preferred categories.
#
#   suitability_rank_scale.csv
#     Maps Unsuitable/Suitable/Highly suitable/Other to an ordinal rank
#     (0 = best). See 02_build_resistance_rasters.R for how rank becomes
#     resistance.
#
# The per-species table-building logic itself (tiers -> ranks, forest
# fallback, water/road adjustments by dispersal_type, age penalties) lives
# in build_rank_tables() in pipeline_functions.R, so the sensitivity
# analysis (08) builds its perturbed tables with exactly the same code.
# ============================================================

source(file.path("Scripts", "Analysis", "setup_script.R"))

# Load in the raw data
species_raw <- readxl::read_excel(here(tabular_dir, "specieslist.xlsx"), sheet = "top10")

landcover_lookup <- read_csv(here(tabular_dir, "Lookup_tables", "landcover_class_lookup.csv"), show_col_types = FALSE)
forest_age_lookup <- read_csv(here(tabular_dir, "Lookup_tables", "forest_age_class_lookup.csv"), show_col_types = FALSE)
forest_type_lookup <- read_csv(here(tabular_dir, "Lookup_tables", "forest_type_class_lookup.csv"), show_col_types = FALSE)
rank_scale <- read_csv(here(tabular_dir, "Lookup_tables", "suitability_rank_scale.csv"), show_col_types = FALSE)

metadata_rows <- list()

## Run the loop to generate
for (i in seq_len(nrow(species_raw))) {
  sp <- species_raw$species[i]
  message("Processing: ", sp, " (", species_raw$en_name[i], ")")

  ## forest type ranks (forest cells only), non-forest land-cover ranks
  ## (incl. the Forest fallback row), and forest age penalties
  rank_tbls <- build_rank_tables(species_raw[i, ], landcover_lookup, forest_type_lookup,
                                 forest_age_lookup, rank_scale)

  forest_type_rank_path <- file.path(tabular_dir, "Suitability", paste0(sp, "_forest_type_rank.csv"))
  write_csv(rank_tbls$forest_type, forest_type_rank_path)

  landcover_rank_path <- file.path(tabular_dir, "Suitability", paste0(sp, "_landcover_rank.csv"))
  write_csv(rank_tbls$landcover, landcover_rank_path)

  age_penalty_path <- file.path(tabular_dir, "Suitability", paste0(sp, "_age_penalty.csv"))
  write_csv(rank_tbls$age_penalty, age_penalty_path)

  # ---- SDM raster path (used only by script 07) ----
  # Intermediate/00_SDM_cropped/{Genus species}.tif -- the pre-cropped/
  # reprojected copy 00-PA_preparation.R produces from
  # Data/Spatial/SDMs/, named by scientific name (not species code), so
  # derive it from sci_name (genus + species, dropping any subspecies
  # epithet). Species without a matching file (SDMs aren't expected for
  # every species) get NA -- script 07 skips those.
  sdm_genus_species <- stringr::word(species_raw$sci_name[i], 1, 2)
  sdm_path_candidate <- file.path(sdm_cropped_dir, paste0(sdm_genus_species, ".tif"))
  sdm_raster_path <- if (file.exists(sdm_path_candidate)) sdm_path_candidate else NA_character_
  if (is.na(sdm_raster_path)) {
    warning(sp, ": no SDM raster found at ", sdm_path_candidate,
            " -- run 00-PA_preparation.R's SDM cropping step first if this ",
            "is unexpected, otherwise script 07 will skip this species.")
  }

  metadata_rows[[sp]] <- tibble(
    species              = sp,
    landcover_rank_csv   = landcover_rank_path,
    forest_type_rank_csv = forest_type_rank_path,
    age_penalty_csv      = age_penalty_path,
    sdm_raster_path      = sdm_raster_path,
    kernel_ref_mode      = "median",
    kernel_ref_distance  = species_raw$dispersal_km[i] * 1000,  # km -> m
    kernel_ref_prob      = NA_real_,
    kernel_p_threshold   = kernel_p_threshold
  )
}

species_metadata <- bind_rows(metadata_rows)
write_csv(species_metadata, paste0(here(tabular_dir), "/species_metadata.csv"))
