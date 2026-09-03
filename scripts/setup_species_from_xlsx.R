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
# Stage 0: Convert data/raw/specieslist.xlsx ("top10" sheet) into this
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
#     for any land-cover class outside decid/mixed/conif/wetland/grass/shrub
#
# Run this ONCE (or whenever specieslist.xlsx changes), BEFORE stage 01.
# -- review its output before running the rest of the pipeline on it.
#
# Requires three small mapping tables you fill in once
#
#   landcover_class_lookup.csv
#     Maps YOUR landcover.tif class codes to decid/mixed/conif/wetland/
#     grass/shrub/other. decid/mixed/conif are treated as "forest" --
#     only those cells get an age penalty applied.
#
#   forest_age_class_lookup.csv
#     Maps YOUR forest_age.tif class codes to jeune/mature/forêt ancienne
#     (matching the spreadsheet's forest_age column).
#
#   suitability_rank_scale.csv
#     Maps Unsuitable/Suitable/Highly suitable/Other to an ordinal rank
#     (0 = best). See 02_build_resistance_rasters.R for how rank becomes
#     resistance.
# ============================================================

source(file.path("scripts", "setup_script.R"))

# Load in the raw data
species_raw <- readxl::read_excel(here(data_dir, "specieslist.xlsx"), sheet = "top10")

landcover_lookup <- read_csv(here("data", "lookup_tables", "landcover_class_lookup.csv"), show_col_types = FALSE)
forest_age_lookup <- read_csv(here("data", "lookup_tables", "forest_age_class_lookup.csv"), show_col_types = FALSE)
forest_type_lookup <- read_csv(here("data", "lookup_tables", "forest_type_class_lookup.csv"), show_col_types = FALSE)
rank_scale <- read_csv(here("data", "lookup_tables", "suitability_rank_scale.csv"), show_col_types = FALSE)

non_forest_categories <- c(wetland = "wetland", grass = "grass", shrub = "shrub",
                           agriculture = "agriculture")
forest_type_cols <- c(conif = "conif", decid = "decid", mixed = "mixed")

# Water (LULC_code 4) isn't in non_forest_categories -- it uses a shared
# fixed_rank (landcover_class_lookup.csv) like Urban/roads, adjusted per
# species by dispersal_type: species that cross water cheaply (flying,
# aquatic) get a lower rank; fossorial species get a higher one.
dispersal_type_water_adjustment <- c(flying = -1, aquatic = -1, fossorial = 1, terrestrial = 0)

# default penalty (in rank steps, i.e. doublings) applied to forest cells
# whose age class isn't in a species' preferred set -- edit the resulting
# {species}_age_penalty.csv files directly if you want this to vary by
# species rather than using one global default
age_mismatch_penalty <- 1

metadata_rows <- list()

## Run the loop to generate
for (i in seq_len(nrow(species_raw))) {
  sp <- species_raw$species[i]
  message("Processing: ", sp, " (", species_raw$en_name[i], ")")
  
  # ---- forest type rank table (applied only within forest LULC cells) ----
  # (computed before the land-cover table below, so its average can seed
  # the Forest fallback rank used there)
  sp_rank_by_type <- tibble(
    type_category    = names(forest_type_cols),
    suitability = unlist(species_raw[i, forest_type_cols], use.names = FALSE)
  ) %>%
    left_join(rank_scale %>% dplyr::select(suitability, rank), by = "suitability")

  forest_type_rank_tbl <- forest_type_lookup %>%
    left_join(sp_rank_by_type %>% dplyr::select(type_category, rank), by = "type_category") %>%
    dplyr::select(forest_type_class, rank)

  if (any(is.na(forest_type_rank_tbl$rank))) {
    warning(sp, ": some forest_type classes have no resolvable rank -- check ",
            "forest_type_class_lookup.csv's type_category values match conif/decid/mixed.")
  }

  forest_type_rank_path <- file.path(data_dir, "suitability", paste0(sp, "_forest_type_rank.csv"))
  write_csv(forest_type_rank_tbl, forest_type_rank_path)

  # ---- non-forest land-cover rank table ----
  sp_rank_by_category <- tibble(
    habitat_category = names(non_forest_categories),
    suitability = unlist(species_raw[i, non_forest_categories], use.names = FALSE)
  ) %>%
    left_join(rank_scale %>% dplyr::select(suitability, rank), by = "suitability")

  if (any(is.na(sp_rank_by_category$rank))) {
    warning(sp, ": suitability tier not found in suitability_rank_scale.csv for some ",
            "habitat categories -- check spelling matches exactly (e.g. 'Highly suitable').")
  }

  # Forest (LULC_code 2) fallback -- used by 02-resistancecost_prep.R only
  # for forest cells where forest_type/forest_age data is missing. This
  # species' average rank across its own conif/decid/mixed suitability,
  # i.e. assume typical forest suitability when the exact type/age can't
  # be determined.
  forest_fallback_rank <- mean(sp_rank_by_type$rank, na.rm = TRUE)

  water_adjustment <- dispersal_type_water_adjustment[species_raw$dispersal_type[i]]
  if (is.na(water_adjustment)) {
    warning(sp, ": dispersal_type '", species_raw$dispersal_type[i],
            "' not recognized -- no water rank adjustment applied.")
    water_adjustment <- 0
  }

  landcover_rank_tbl <- landcover_lookup %>%
    filter(LULC_code != 2) %>%
    left_join(sp_rank_by_category %>% dplyr::select(habitat_category, rank), by = "habitat_category") %>%
    mutate(rank = coalesce(fixed_rank, rank)) %>%
    mutate(rank = if_else(LULC_code == 4, pmax(rank + water_adjustment, 0), rank)) %>%
    dplyr::select(cover_class = LULC_code, rank) %>%
    bind_rows(tibble(cover_class = 2, rank = forest_fallback_rank))

  if (any(is.na(landcover_rank_tbl$rank))) {
    warning(sp, ": ", sum(is.na(landcover_rank_tbl$rank)),
            " non-forest land-cover class(es) have no resolvable rank -- check ",
            "landcover_class_lookup.csv's habitat_category/fixed_rank values.")
  }

  landcover_rank_path <- file.path(data_dir, "suitability", paste0(sp, "_landcover_rank.csv"))
  write_csv(landcover_rank_tbl, landcover_rank_path)

  # ---- forest age penalty table (applied only within forest LULC cells) ----
  # species_raw's forest_age column is a comma-separated list of this
  # species' preferred age classes, e.g. "mature, forêt ancienne"
  preferred_ages <- species_raw$forest_age[i] %>%
    str_split(",") %>%
    unlist() %>%
    str_trim()
  
  age_penalty_tbl <- forest_age_lookup %>%
    mutate(
      penalty = if_else(age_category %in% preferred_ages, 0L, age_mismatch_penalty)
    ) %>%
    dplyr::select(forest_age_class, penalty)
  
  age_penalty_path <- file.path(data_dir, "suitability", paste0(sp, "_age_penalty.csv"))
  write_csv(age_penalty_tbl, age_penalty_path)

  # ---- SDM raster path (used only by script 07) ----
  # data/intermediate/SDM_cropped/{Genus species}.tif -- the pre-cropped/
  # reprojected copy 00-PA_preparation.R produces from
  # data/drive_data/SDMs/, named by scientific name (not species code), so
  # derive it from sci_name (genus + species, dropping any subspecies
  # epithet). Species without a matching file (SDMs aren't expected for
  # every species) get NA -- script 07 skips those.
  sdm_genus_species <- stringr::word(species_raw$sci_name[i], 1, 2)
  sdm_path_candidate <- file.path(interm_dir, "SDM_cropped", paste0(sdm_genus_species, ".tif"))
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
    kernel_p_threshold   = 0.05
  )
}

species_metadata <- bind_rows(metadata_rows)
write_csv(species_metadata, paste0(here(data_dir), "/species_metadata.csv"))
