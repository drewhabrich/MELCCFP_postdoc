## ------------------------------- ##
## Script name: 02-resistancecost_prep
##
## Purpose of script: Generate the species specific trasition matrix that reflects
##  the movement cost across different land cover types, forest types, and forest ages
#
# Saves:
#   Output/02_resistance/{species}_res.tif     (resistance raster, one per species)
#   Output/Tables/02_landcover_resistance_by_species.csv  (one combined table: species x landcover class -> resistance, for QA)
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
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))

# 1. Load required files ---------------------------------------
# land cover / forest type / forest age, aligned to the shared 30 m
# reference grid (template_raster, from setup_script.R) once, up front, if
# needed (nearest neighbor -- all categorical) -- see load_resistance_inputs()
resistance_inputs <- load_resistance_inputs()
landcover   <- resistance_inputs$landcover
forest_type <- resistance_inputs$forest_type
forest_age  <- resistance_inputs$forest_age
forest_mask <- resistance_inputs$forest_mask
species_meta <- read_csv(here(tabular_dir, "species_metadata.csv"), show_col_types = FALSE)
# species_table:
specieslist <- readxl::read_excel(here(tabular_dir, "specieslist.xlsx"))

# landcover / forest-type class names, for labeling the per-species rank/resistance table below
landcover_lookup <- read_csv(here(tabular_dir, "Lookup_tables", "landcover_class_lookup.csv"), show_col_types = FALSE)
forest_type_lookup <- read_csv(here(tabular_dir, "Lookup_tables", "forest_type_class_lookup.csv"), show_col_types = FALSE)

## Create the output directories (resistance_dir is set in setup_script.R)
if (!dir.exists(resistance_dir)) dir.create(resistance_dir)
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

## Collects each species' landcover_resistance_tbl (below); combined into
## one species x landcover table after the loop.
species_landcover_resistance <- list()

## TEST SUBSET WITH MAAM AND DOOR
species_subset <- c("MAAM", "ASFL", "DOOR", "DRPI", "AQCH")  # e.g. NULL for all species
if (!is.null(species_subset)) {
  species_meta <- species_meta %>% filter(species %in% species_subset)
}

overwrite <- TRUE  # set FALSE to skip species that already have a {species}_res.tif

## rank of the "Unsuitable" tier -- forest types at this rank get a capped age penalty
unsuitable_rank <- read_csv(here(tabular_dir, "Lookup_tables", "suitability_rank_scale.csv"),
                            show_col_types = FALSE) %>%
  filter(suitability == "Unsuitable") %>% pull(rank)

for (i in seq_len(nrow(species_meta))) {
  sp <- species_meta$species[i]
  out_path <- here(resistance_dir, paste0(sp, "_res.tif"))

  landcover_rank_tbl <- read_csv(species_meta$landcover_rank_csv[i], show_col_types = FALSE)
  forest_type_rank_tbl <- read_csv(species_meta$forest_type_rank_csv[i], show_col_types = FALSE)
  age_penalty_tbl <- read_csv(species_meta$age_penalty_csv[i], show_col_types = FALSE)

  ## Rank -> resistance (2^rank) for this species' landcover classes AND
  ## forest types (conif/decid/mixed -- each type's own resistance, not
  ## combined with an age penalty), labeled with class names -- stashed for
  ## the combined table written after the loop. Kept outside the
  ## overwrite check below (cheap, CSV-only) so the combined table stays
  ## complete even for species whose raster is being skipped.
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

  if (file.exists(out_path) && !overwrite) {
    message(sp, ": ", out_path, " already exists -- skipping (set overwrite <- TRUE to redo).")
    next
  }

  message("Building resistance raster for: ", sp)

  ## resistance = 2^total_rank, forest cells = forest-type rank + age
  ## penalty (falling back to the Forest landcover row where type/age data
  ## is missing) -- see build_resistance() in pipeline_functions.R. Returns
  ## this species' own worst-case rank, from its actual tables -- NOT
  ## hardcoded, so it stays correct if the tables change.
  max_rank <- build_resistance(
    list(landcover = landcover_rank_tbl, forest_type = forest_type_rank_tbl, age_penalty = age_penalty_tbl),
    landcover, forest_type, forest_age, forest_mask, out_path,
    forest_unsuitable_rank = unsuitable_rank
  )

  message("  -> saved ", out_path, " (max_rank = ", max_rank, ")")
}

## Combine every species' landcover_resistance into one table: species in
## rows, landcover class in columns, resistance (2^rank) as the value.
landcover_resistance_table <- bind_rows(species_landcover_resistance) %>%
  tidyr::pivot_wider(names_from = class, values_from = resistance)

landcover_resistance_path <- here(table_dir, "02_landcover_resistance_by_species.csv")
write_csv(landcover_resistance_table, landcover_resistance_path)
message("-> saved combined species x landcover resistance table to ", landcover_resistance_path)
