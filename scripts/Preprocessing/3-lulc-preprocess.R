## Pre-process LULC data layers
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This script pre-processes the land use land cover (LULC) data from three 
## data sources: utilisation du territoire for Quebec, AAFC  for Ontario and New
## Brunswick, and NLCD for the United States. Steps involve cropping data to the
## modelling extent, and harmonizing  classes based on crosswalk tables. In the 
## case of the Quebec data, data is also aggregated to the target spatial 
## resolution for the project (30 m). Then, the LULC maps for each jurisdiction
## are merged and NAs in the boundaries are resolved with a focal window 
## approach. 
## Based on the NLCD data, this script also extract the forest cover types for 
## the United States extent.
##
## Intermediate outputs:
## - tif, mask of raw LULC map for Quebec at the modelling extent
## - tif, aggregated and harmonized LULC map for Quebec
## - tif, harmonized LULC map for AAFC
## - tif, mask of raw LULC map for NLCD at the modelling extent
## - tif, harmonized LULC map for NLCD
## - tif, forest cover types for the United States jurisdiction
## - tif, combined LULC map

## Output:
## - tif, water pixels



# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))



# Quebec -----------------------------------------------------------------------

# Quebec LULC
qc_lulc <- rast(file.path(
  spatial_dir, 
  "LULC", 
  "utilisation_territoire_2020", 
  "utilisation_territoire_2020.tif"
))

# Mask to modelling extent
qc_lulc_mask <- model_mask(qc_lulc)
rm(qc_lulc)

# Write to file
writeRaster(
  x = qc_lulc_mask, 
  filename = file.path(intermediate_dir, "LULC", "qc_lulc_mask.tif"),
  datatype = "INT4S",
  overwrite = TRUE)

# Load crosswalk table to simple classes
qc_reclass_table <- read.csv(
  file.path(tabular_dir, "Crosswalks", "qc_lulc_reclass.csv"),
  header = TRUE
)

## Separate forested and open wetlands -----------

# Get pure forest classes
# NOTE: Not those combined with wetland classes (> 1000)
qc_table_forest <- qc_reclass_table %>%
  filter(class_name == "Forestier") %>%
  filter(qc_lulc < 1000)

# Get wetland classes and get source class without wetland type
qc_table_wetland <- qc_reclass_table %>%
  filter(class_name == "Humide") %>%
  mutate(source_lulc = qc_table_wetland$qc_lulc %% 1000)

# For each pure forest class
for(i in qc_table_forest$qc_lulc){
  
  # Reclassify wetland into forested wetland based on forest classes
  qc_table_wetland$class_name[
    qc_table_wetland$source_lulc == i] <- "Humide boisé"
  qc_table_wetland$class_id[
    qc_table_wetland$source_lulc == i] <- 32
}

# Reclassify remaining wetland into open wetland
qc_table_wetland$class_name[
  qc_table_wetland$class_name == "Humide"] <- "Humide ouvert"
qc_table_wetland$class_id[
  qc_table_wetland$class_id == 3] <- 31

# Combine with crosswalk table
qc_reclass_updated <- qc_reclass_table %>%
  filter(class_name != "Humide") %>%
  rbind(qc_table_wetland[,1:3])


## Reclassify ------------------------------------

# Simplify classes
qc_lulc_simple <- classify(x = qc_lulc_mask, rcl = qc_reclass_updated[,1:2])

# Write to file
# writeRaster(
#   x = qc_lulc_simple,
#   filename = file.path(intermediate_dir, "LULC", "qc_lulc_simple.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)

# Aggregate to 30 m resolution
qc_lulc_30m <- aggregate(qc_lulc_simple, fact = 3, fun = "modal")

# Write to file
# writeRaster(
#   x = qc_lulc_30m,
#   filename = file.path(intermediate_dir, "LULC", "qc_lulc_30m.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)


## Complement grassland --------------------------

# Target class to keep (122 = culture pérenne et pâturage)
target_class <- 122

# Reclass matrix
rcl_matrix <- matrix(c(-Inf, target_class - 1, NA,
                       target_class + 1, Inf, NA), 
                     ncol = 3, byrow = TRUE)

# Reclassify
qc_lulc_crop <- classify(qc_lulc_mask, rcl_matrix)

# Load annual crop inventory data
crop_inv <- rast(file.path(
  spatial_dir, 
  "LULC", 
  "aci_2020_qc_v4.tif"
))

# Mask to modelling extent
crop_inv_mask <- model_mask(crop_inv)

# Project to template raster
crop_inv_model <- project(crop_inv_mask, template_raster, method = "near")
qc_crop_model <- project(qc_lulc_crop, template_raster, method = "near")

# Mask to target class pixels
crop_inv_class <- mask(crop_inv_model, qc_crop_model)

# Load crosswalk table
aci_reclass_table <- read.csv(
  file.path(tabular_dir, "Crosswalks", "aci_lulc_reclass.csv"),
  header = TRUE
)

# Get grassland pixels in Quebec
aci_grassland <- classify(x = crop_inv_class, rcl = aci_reclass_table[,c(1,4)])

# Set everything else to NA
qc_grassland <- ifel(aci_grassland == 7, 7, NA)

# Write to file
# writeRaster(
#   x = qc_grassland,
#   filename = file.path(intermediate_dir, "LULC", "qc_grassland.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)

# Merge with LULC
# NOTE: The order of rasters is essential here.
qc_lulc_grass <- merge(qc_grassland, qc_lulc_30m)

# Write to file
# writeRaster(
#   x = qc_lulc_grass,
#   filename = file.path(intermediate_dir, "LULC", "qc_lulc_grass.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)


## Fix NAs ---------------------------------------

# Set NA values (outside modelling extent) to 0 
qc_lulc_grass[is.na(qc_lulc_grass)] <- 0

# Set 9999 values to NA
# NOTE: These represent pixels without classification within the study area
qc_lulc_NA <- classify(
  x = qc_lulc_grass,
  rcl = matrix(c(9999, NA), nrow = 1, byrow = TRUE)
)

# Set condition
condition <- freq(qc_lulc_NA, value = NA)$count[1]

# Repeat until all NA cells are filled
while (condition > 0) {
  
  # Substitute NA values with mode from neighboring cells 
  qc_lulc_NA <- focal(
    qc_lulc_NA, 
    w = 3, 
    fun = "modal", 
    na.policy = "only", 
    na.rm = TRUE
  )
  
  # Update condition value
  condition <- freq(qc_lulc_NA, value = NA)$count[1]
  cat(paste0(condition, " "))
}

# Set 0's to NA
qc_lulc_NA <- classify( 
  x = qc_lulc_NA,
  rcl = matrix(c(0, NA), nrow = 1, byrow = TRUE)
)

# Write to file
writeRaster(
  x = qc_lulc_NA,
  filename = file.path(intermediate_dir, "LULC", "qc_lulc_NA.tif"),
  datatype = "INT4S",
  overwrite = TRUE)

# Project to template raster
qc_lulc_model <- project(qc_lulc_NA, template_raster, method = "near")

# Write to file
writeRaster(
  x = qc_lulc_model, 
  filename = file.path(intermediate_dir, "LULC", "qc_lulc_model.tif"),
  datatype = "INT4S",
  overwrite = TRUE)



# AAFC -------------------------------------------------------------------------

# Load data
aafc_u17 <- rast(file.path(spatial_dir, "LULC", 
                           "LU2020_u17_v5", "LU2020_u17_v5.tif"))
aafc_u18 <- rast(file.path(spatial_dir, "LULC", 
                           "LU2020_u18_v5", "LU2020_u18_v5.tif"))
aafc_u19 <- rast(file.path(spatial_dir, "LULC", 
                           "LU2020_u19_v5", "LU2020_u19_v5.tif"))
aafc_u20 <- rast(file.path(spatial_dir, "LULC", 
                           "LU2020_u20_v5", "LU2020_u20_v5.tif"))

# Mask to modelling extent
aafc_u17_mask <- model_mask(aafc_u17)
aafc_u18_mask <- model_mask(aafc_u18)
aafc_u19_mask <- model_mask(aafc_u19)
aafc_u20_mask <- model_mask(aafc_u20)

# Project to template raster
aafc_u17_model <- project(aafc_u17_mask, template_raster, method = "near")
aafc_u18_model <- project(aafc_u18_mask, template_raster, method = "near")
aafc_u19_model <- project(aafc_u19_mask, template_raster, method = "near")
aafc_u20_model <- project(aafc_u20_mask, template_raster, method = "near")

# Write to file
# writeRaster(
#   x = aafc_u17_model,
#   filename = file.path(intermediate_dir, "LULC", "aafc_u17_model.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)
# writeRaster(
#   x = aafc_u18_model,
#   filename = file.path(intermediate_dir, "LULC", "aafc_u18_model.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)
# writeRaster(
#   x = aafc_u19_model,
#   filename = file.path(intermediate_dir, "LULC", "aafc_u19_model.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)
# writeRaster(
#   x = aafc_u20_model,
#   filename = file.path(intermediate_dir, "LULC", "aafc_u20_model.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)

# Combine AAFC slices
aafc_model <- merge(
  aafc_u17_model, aafc_u18_model, aafc_u19_model, aafc_u20_model
  )

# Write to file
writeRaster(
  x = aafc_model,
  filename = file.path(intermediate_dir, "LULC", "aafc_model.tif"),
  datatype = "INT4S",
  overwrite = TRUE)

# Load crosswalk table to simple classes
aafc_reclass_table <- read.csv(
  file.path(tabular_dir, "Crosswalks", "aafc_lulc_reclass.csv"),
  header = TRUE
)

# Simplify classes
aafc_model_class <- classify(x = aafc_model, rcl = aafc_reclass_table[,c(1,4)])

# Write to file
writeRaster(
  x = aafc_model_class, 
  filename = file.path(intermediate_dir, "LULC", "aafc_model_class.tif"),
  datatype = "INT4S",
  overwrite = TRUE)



# NLCD -------------------------------------------------------------------------

# Load data
nlcd_tile_1 <- rast(file.path(
  spatial_dir, 
  "LULC", 
  "NLCD_tile_1",
  "Annual_NLCD_LndCov_2020_CU_C1V1_3ec067df-979a-4b3c-b17f-3187ae254df0.tiff"
  ))
nlcd_tile_2 <- rast(file.path(
  spatial_dir, 
  "LULC", 
  "NLCD_tile_2",
  "Annual_NLCD_LndCov_2020_CU_C1V1_494f78d7-dff9-43a4-8b72-7865f58a8f56.tiff"
))

# Mask to modelling extent
nlcd_mask_1 <- model_mask(nlcd_tile_1)
nlcd_mask_2 <- model_mask(nlcd_tile_2)

# Merge NLCD tiles
nlcd_mask <- merge(nlcd_mask_1, nlcd_mask_2)

# Project to template raster
nlcd_model <- project(nlcd_mask, template_raster, method = "near")

# Write to file
writeRaster(
  x = nlcd_model,
  filename = file.path(intermediate_dir, "LULC", "nlcd_model.tif"),
  datatype = "INT4S",
  overwrite = TRUE)

# Load crosswalk table to simple classes
nlcd_reclass_table <- read.csv(
  file.path(tabular_dir, "Crosswalks", "nlcd_lulc_reclass.csv"),
  header = TRUE
)

# Simplify classes
nlcd_model_class <- classify(x = nlcd_model, rcl = nlcd_reclass_table[,c(1,4)])

# Write to file
writeRaster(
  x = nlcd_model_class, 
  filename = file.path(intermediate_dir, "LULC", "nlcd_model_class.tif"),
  datatype = "INT4S",
  overwrite = TRUE)


## Forest cover type --------------------------------------------

# Load crosswalk table to get cover types
nlcd_cover_reclass_table <- read.csv(
  file.path(tabular_dir, "Crosswalks", "nlcd_cover_reclass.csv"),
  header = TRUE
)

# Get forest type raster
nlcd_cover_type <- classify(
  x = nlcd_model, 
  rcl = nlcd_cover_reclass_table[,c(1,4)]
  )

# Write to file
writeRaster(
  x = nlcd_cover_type, 
  filename = file.path(intermediate_dir, "LULC", "nlcd_cover_type.tif"),
  datatype = "INT4S",
  overwrite = TRUE)



# Combine ----------------------------------------------------------------------

# Merge data layers
lulc_model_merge <- merge(qc_lulc_model, aafc_model_class, nlcd_model_class)

# Write to file
# writeRaster(
#   x = lulc_model_merge,
#   filename = file.path(intermediate_dir, "LULC", "lulc_model_merge.tif"),
#   datatype = "INT4S",
#   overwrite = TRUE)


## Fix NAs in jurisdiction boundaries ---------------------------

#lulc_model_NA <- rast(file.path(output_dir, "lulc_model_merge.tif"))

# Create copy and set NAs to 0
lulc_model_zero <- lulc_model_merge
lulc_model_zero[is.na(lulc_model_zero)] <- 0 

# Load modelling extent
model_area <- st_read(dsn = file.path(output_dir), layer = "modelling_extent")

# Mask to modelling extent
lulc_inv <- mask(lulc_model_zero, model_area, inverse = TRUE)

# Merge modelling extent set to 0 and LULC map with internal NAs
lulc_fix <- merge(lulc_model_merge, lulc_inv)

# Set condition
condition <- freq(lulc_fix, value = NA)$count[1]

# Repeat until all NA cells are filled
while (condition > 0) {
  
  # Substitute NA values with mode from neighboring cells 
  lulc_fix <- focal(
    lulc_fix, 
    w = 3, 
    fun = "modal", 
    na.policy = "only", 
    na.rm = TRUE
  )
  
  # Update condition value
  condition <- freq(lulc_fix, value = NA)$count[1]
  cat(paste0(condition, " "))
}

# Mask to modelling extent
lulc_model <- mask(lulc_fix, model_area)

# Write to file
writeRaster(
  x = lulc_model, 
  filename = file.path(intermediate_dir, "LULC", "lulc_model.tif"),
  datatype = "INT4S",
  overwrite = TRUE)

# Get water pixels
water_mask <- lulc_model == 4

# Mask to keep only water pixels
water_raster <- mask(lulc_model, water_mask, maskvalues = FALSE)

# Write to file
writeRaster(
  x = water_raster, 
  filename = file.path(output_dir, "water_raster.tif"),
  datatype = "INT4S",
  overwrite = TRUE)


