## Mask selected species SDM to modelling extent
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This script loops through a list of selected species to, one at a time, open
## their Canada-wide SDM and mask it to the modelling extent.
##
## Outputs:
## - tifs, SDM of selected species masked to the modelling extent.



# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))



# Load data --------------------------------------------------------------------

# Modelling extent
model_area <- st_read(
  dsn = file.path(output_dir), 
  layer = "modelling_extent"
  )

# Study area
study_area <- st_read(
  dsn = file.path(output_dir), 
  layer = "study_area"
) %>%
  st_transform(target_crs)



# Mask Eckart SDMs -------------------------------------------------------------

# Selected species
spp_list <- c("Aquila chrysaetos", "Castor canadensis", "Clethrionomys gapperi",
              "Dolichonyx oryzivorus", "Hylatomus pileatus", "Lynx canadensis",
              "Martes americana", "Plethodon cinereus",
              "Asio flammeus")

# Loop through species list
for(i in 1:length(spp_list)){
  
  # Get species name
  species_name <- spp_list[i]
  
  # Load species' SDM for Canada
  spp_sdm_ca <- rast(
    file.path(spatial_dir, "SDMs", paste0(species_name, ".tif"))
    #file.path("/Volumes/KINGSTON/Spatial data/SDMs", paste0(species_name, ".tif"))
    )  
  
  # Transform modeling extent to same CRS as SDM
  model_area_crs <- st_transform(model_area, crs = crs(spp_sdm_ca))
  
  # Clip and mask to modelling extent, and project to target CRS
  spp_sdm <- spp_sdm_ca %>%
    crop(model_area_crs) %>%
    mask(model_area_crs)
  rm(spp_sdm_ca)
  
  # Reproject to template raster
  spp_sdm_model <- project(spp_sdm, target_crs)
  
  # Write to file
  writeRaster(
    x = spp_sdm_model,
    filename = file.path(
      intermediate_dir, 
      "SDMs", 
      paste0(species_name, "_model.tif")),
    overwrite = TRUE,
    datatype = "FLT4S"
    )
  rm(spp_sdm_model)
  
  # Clear temp files
  unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
}



# Mask MELCCFP SDM -------------------------------------------------------------

# Load data
sdm_reptile <- st_read(
  dsn = file.path(spatial_dir, "SDMs", "Aires_repartition_reptiles.gdb"), 
  layer = "Aires_repartition_reptiles"
)

# Filter to "Glyptemys insculpta"
sdm_turtle <- sdm_reptile %>%
  filter(NOM_SCIENT == "Glyptemys insculpta")

# Transform modeling extent to same CRS as SDM
model_area_crs <- st_transform(model_area, crs = crs(sdm_turtle))

# Clip and mask to modelling extent, and project to target CRS
sdm_turtle_model <- sdm_turtle %>%
  st_intersection(model_area_crs)

# Reproject to template raster
sdm_turtle_raster <- rasterize(sdm_turtle_model, template_raster)

# Write to file
writeRaster(
  x = sdm_turtle_raster,
  filename = file.path(
    intermediate_dir, 
    "SDMs", 
    paste0("Glyptemys insculpta", "_model.tif")),
  overwrite = TRUE,
  datatype = "INT4S"
)



