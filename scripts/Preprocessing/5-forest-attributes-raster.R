## Merge forest attributes across datasets
##
## Written by Carina Rauen Firkowski
## 
## Description:
## For the target attributes age and cover type, this script loads the 
## intermediate outputs for each jurisdiction and combines. For age, 
## harmonization across jurisdiction data is done through crosswalks. Lastly, 
## attribute pixels are harmonized with the LULC map so that only forest pixels 
## have an associated age and forest cover type.
##
## Outputs:
## - tif, forest cover type raster
## - tif, forest age raster


# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))



# Cover type -------------------------------------------------------------------

## Load data ----------------------------------------------------

# Quebec
qc_cover_type <- rast(file.path(
  intermediate_dir, 
  "Ecoforestry", 
  paste0("ecoforestry_model_type_couv.tif")
  ))

# NLCD
nlcd_cover_type <- rast(
  file.path(intermediate_dir, "LULC", "nlcd_cover_type.tif")
)

# North American Land Cover
can_land_cover <- rast(
  file.path(
    spatial_dir, 
    "LULC", 
    "can_land_cover_2020v2_30m_tif", 
    "CAN_NALCMS_landcover_2020v2_30m",
    "data",
    "CAN_NALCMS_landcover_2020v2_30m.tif")
)

# Modelling extent
model_area <- st_read(dsn = file.path(output_dir),
                      layer = "modelling_extent")



## Ontario & New Brunswick --------------------------------------

# Mask to modelling extent
can_land_cover_mask <- model_mask(can_land_cover)
rm(can_land_cover)

# Load reclass table
na_reclass_table <- read.csv(
  file.path(tabular_dir, "Crosswalks", "na_lulc_reclass.csv")
  )

# Reclassify raster
can_cover_type <- classify(
  x = can_land_cover_mask, 
  rcl = na_reclass_table[,c(1,4)]
  )

# Project to template raster
can_cover_type_model <- project(
  can_cover_type, template_raster, method = "near"
)

# Save to file
# writeRaster(
#     x = can_cover_type_model,
#     filename = file.path(intermediate_dir, "LULC", "can_cover_type_model.tif"),
#     datatype = "INT4S",
#     overwrite = TRUE)


## Combine data -------------------------------------------------

# Merge data layers
cover_type_merged <- merge(qc_cover_type, nlcd_cover_type, can_cover_type_model) 

# Open LULC data
lulc_model <- rast(file.path(intermediate_dir, "LULC", "lulc_model.tif"))

# Get forest pixels only
lulc_forest <- ifel(lulc_model == 2, 1, NA)

# Mask cover types to forest pixels only data layers
cover_type_model <- mask(cover_type_merged, lulc_forest) 

# Save to file
writeRaster(
  x = cover_type_model,
  filename = file.path(output_dir, "cover_type_model.tif"),
  datatype = "INT4S",
  overwrite = TRUE)



# Age --------------------------------------------------------------------------

## Ontario & New Brunswick --------------------------------------

# Load data
age_can <- rast(file.path(
  spatial_dir, 
  "Ecoforestry", 
  "CA_forest_age_2022",
  "CA_forest_age_2022.tif"))

# Mask to modelling extent
age_can_mask <- model_mask(age_can)

# Write to file
# writeRaster(
#   x = age_can_mask,
#   filename = file.path(intermediate_dir, "Ecoforestry", "age_can_mask.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)

# Load reclass table
age_can_reclass_table <- read.csv(
  file.path(tabular_dir, "Crosswalks", "age_can_reclass.csv")
)

# Reclassify raster
age_can_class <- classify(
  x = age_can_mask, 
  rcl = age_can_reclass_table[,c(1,2,4),],
  right = FALSE,
  include.lowest = TRUE
)
freq(age_can_class)

# Project to template raster
age_can_model <- project(age_can_class, template_raster, method = "near")

# Write to file
# writeRaster(
#   x = age_can_model,
#   filename = file.path(intermediate_dir, "Ecoforestry", "age_can_model.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)


## U.S. ---------------------------------------------------------

# Get list of tiles
us_age_list <- list.files(
    path = file.path(
      spatial_dir, 
      "Ecoforestry", 
      "CMS_Global_Forest_Age_2345_1-20260529_192915"),
    pattern = ".tif$")

# For each tile
for(age_tile in us_age_list){
  
  # Load data
  age_tile_rast <- rast(file.path(
    spatial_dir, 
    "Ecoforestry", 
    "CMS_Global_Forest_Age_2345_1-20260529_192915",
    age_tile))
  
  # Mask to modelling extent
  age_tile_mask <- model_mask(age_tile_rast)
  
  # Combine
  if(age_tile == us_age_list[1]){
    us_age_raster <- age_tile_mask
  } else {
    us_age_raster <- merge(us_age_raster, age_tile_mask)
  }
}

# Project to template raster
us_age_model <- project(us_age_raster, template_raster, method = "near")

# Write to file
# writeRaster(
#   x = us_age_model,
#   filename = file.path(intermediate_dir, "Ecoforestry", "us_age_model.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)

# Load reclass table
us_age_reclass <- read.csv(
  file.path(tabular_dir, "Crosswalks", "age_us_reclass.csv")
  )

# Update age class
us_age_model <- us_age_model*10

# Reclassify raster
us_age_class <- classify(
  x = us_age_model, 
  rcl = us_age_reclass[,c(1,4),],
  right = FALSE,
  include.lowest = TRUE
)

# Write to file
# writeRaster(
#   x = us_age_class,
#   filename = file.path(intermediate_dir, "Ecoforestry", "us_age_class.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)


## Combine data -------------------------------------------------

# Quebec
qc_age_rast <- rast(file.path(
  intermediate_dir, 
  "Ecoforestry", 
  paste0("ecoforestry_model_cl_age.tif")
)) %>%
  mask(model_area)

# Combine
age_raster <- merge(qc_age_rast, age_can_model, us_age_class)

# Open LULC data
lulc_model <- rast(file.path(output_dir, "lulc_road_model.tif"))

# Get forest pixels only
lulc_forest <- ifel(lulc_model == 2, 1, NA)
#rm(lulc_model) ; rm(qc_age_rast) ; rm(age_can_model) ; rm(us_age_class)

# Mask cover types to forest pixels only data layers
age_model <- mask(age_raster, lulc_forest) 

# Write to file
writeRaster(
  x = age_model,
  filename = file.path(output_dir, "age_model.tif"),
  datatype = "INT4S",
  overwrite = TRUE)
