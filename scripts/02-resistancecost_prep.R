## ------------------------------- ##
## Script name: 02-resistancecost_prep
##
## Purpose of script: Generate the species specific trasition matrix that reflects
##  the movement cost across different land cover types, forest types, and forest ages
#
# Saves:
#   data/output_data/resistance/{species}_res.tif        (resistance raster, one per species)
#   output/tables/landcover_resistance_by_species.csv     (one combined table: species x landcover class -> resistance, for QA)
#
## Author: Andrew Habrich
##
## Notes ------------------------- ##
# ALBERT ET AL. 2017 DOUBLING SCALE: resistance = 2^total_rank, where
# total_rank comes from this species' landcover/forest-type/forest-age
# rank tables (rank 0 = most suitable -> resistance 1; each worse rank
# step doubles the cost). No external base resistance layer is blended
# in -- this project previously combined this with a generic Pither et
# al. 2023 layer via a weighted geometric mean, but that's been dropped
# in favour of scoring resistance purely via Albert et al.'s method.
# pither_model.tif / cost_rast_30m.tif are left on disk as a reference/
# test layer but are no longer read by this script.
#
# Adding a new species only ever means adding a row to
# species_metadata.csv -- this script doesn't change either way.
# ============================================================
source(file.path("scripts", "setup_script.R"))

# 1. Load required files ---------------------------------------
landcover <- rast(here(data_dir, "lulc_road_model.tif"))
fcl_rast <- rast(here(data_dir, "cover_type_model.tif"))
age_rast <- rast(here(data_dir, "age_model.tif"))
species_meta <- read_csv(here(data_dir, "species_metadata.csv"), show_col_types = FALSE)
# species_table:
specieslist <- readxl::read_excel(here(data_dir, "specieslist.xlsx"))

# landcover / forest-type class names, for labeling the per-species rank/resistance table below
landcover_lookup <- read_csv(here("data", "lookup_tables", "landcover_class_lookup.csv"), show_col_types = FALSE)
forest_type_lookup <- read_csv(here("data", "lookup_tables", "forest_type_class_lookup.csv"), show_col_types = FALSE)

# align the three land-cover-derived rasters to the shared 30 m reference
# grid (template_raster, from setup_script.R) once, up front, if needed
# (nearest neighbor -- all categorical)
align_categorical <- function(r, target, label) {
  if (!compareGeom(target, r, stopOnError = FALSE)) {
    message(label, " grid doesn't match the reference grid -- resampling (nearest neighbor).")
    r <- resample(r, target, method = "near")
  }
  r
}
landcover <- align_categorical(landcover, template_raster, "Land cover")
forest_type <- align_categorical(fcl_rast, template_raster, "Forest type")
forest_age  <- align_categorical(age_rast, template_raster, "Forest age")
forest_mask <- landcover == 2

## Define the output directories
resistance_dir <- file.path(output_dir, "resistance")
if (!dir.exists(resistance_dir)) dir.create(resistance_dir)
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

## Collects each species' landcover_resistance_tbl (below); combined into
## one species x landcover table after the loop.
species_landcover_resistance <- list()

## TEST SUBSET WITH MAAM AND DOOR
species_subset <- c("MAAM", "DOOR")  # e.g. NULL for all species
if (!is.null(species_subset)) {
  species_meta <- species_meta %>% filter(species %in% species_subset)
}

for (i in seq_len(nrow(species_meta))) {
  sp <- species_meta$species[i]
  message("Building resistance raster for: ", sp)
  
  landcover_rank_tbl <- read_csv(species_meta$landcover_rank_csv[i], show_col_types = FALSE)
  forest_type_rank_tbl <- read_csv(species_meta$forest_type_rank_csv[i], show_col_types = FALSE)
  age_penalty_tbl <- read_csv(species_meta$age_penalty_csv[i], show_col_types = FALSE)

  ## Rank -> resistance (2^rank) for this species' landcover classes AND
  ## forest types (conif/decid/mixed -- each type's own resistance, not
  ## combined with an age penalty), labeled with class names -- stashed for
  ## the combined table written after the loop.
  landcover_long <- landcover_rank_tbl %>%
    left_join(landcover_lookup, by = c("cover_class" = "LULC_code")) %>%
    mutate(resistance = 2^rank) %>%
    dplyr::select(class = class_en, resistance)

  forest_type_long <- forest_type_rank_tbl %>%
    left_join(forest_type_lookup, by = "forest_type_class") %>%
    mutate(resistance = 2^rank) %>%
    dplyr::select(class = type_category, resistance)

  species_landcover_resistance[[sp]] <- bind_rows(landcover_long, forest_type_long) %>%
    mutate(species = sp) %>%
    dplyr::select(species, class, resistance)

  lulc_rank <- classify(landcover, as.matrix(landcover_rank_tbl))
  forest_type_rank <- classify(forest_type, as.matrix(forest_type_rank_tbl))
  age_penalty <- classify(forest_age, as.matrix(age_penalty_tbl))
  
  forest_total_rank <- forest_type_rank + age_penalty
  # Forest cells where forest_type/forest_age data is missing (forest_total_rank
  # is NA there) fall back to lulc_rank's Forest row -- this species' average
  # conif/decid/mixed rank (see setup_species_from_xlsx.R) -- instead of
  # propagating NA into the resistance raster.
  total_rank <- ifel(forest_mask,
                      ifel(is.na(forest_total_rank), lulc_rank, forest_total_rank),
                      lulc_rank)
  
  # this species' own worst-case rank, from its actual tables -- NOT
  # hardcoded, so it stays correct if the tables change
  max_rank <- max(
    max(landcover_rank_tbl$rank, na.rm = TRUE),
    max(forest_type_rank_tbl$rank, na.rm = TRUE) + max(age_penalty_tbl$penalty, na.rm = TRUE)
  )
  
  # true doubling scale: rank 0 -> 1, each rank step doubles the cost
  species_resistance <- 2^total_rank

  out_path <- here(output_dir, "resistance", paste0(sp, "_res.tif"))
  writeRaster(species_resistance, out_path, overwrite = TRUE)

  message("  -> saved ", out_path, " (max_rank = ", max_rank, ")")
}

## Combine every species' landcover_resistance into one table: species in
## rows, landcover class in columns, resistance (2^rank) as the value.
landcover_resistance_table <- bind_rows(species_landcover_resistance) %>%
  tidyr::pivot_wider(names_from = class, values_from = resistance)

landcover_resistance_path <- here(table_dir, "landcover_resistance_by_species.csv")
write_csv(landcover_resistance_table, landcover_resistance_path)
message("-> saved combined species x landcover resistance table to ", landcover_resistance_path)
