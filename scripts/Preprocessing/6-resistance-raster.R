## Prepare species habitat layers
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This script prepares the resistance layer based on a crosswalk table between
## LULC classes and resistance values. If available, the resistance layer is 
## complemented with the species' habitat model by classify as low resistance 
## pixels of LULC classes representing good habitat that overlap with predicted
## suitable habitat based on the model.
##
## Note: 
## This code currently only works for the short-eared owl, bobolink, and 
## pileated woodpecker. It does not take into consideration forest cover type
## preferences.
##
## Outputs:
## - tif, species-specific resistance raster at modelling extent



# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))

# Load LULC map
lulc_raster <- rast(file.path(output_dir, "lulc_road_model.tif"))

# Load age class map
age_raster <- rast(file.path(output_dir, "age_model.tif"))



# Pre-process habitat data -----------------------------------------------------

# List of species with MQH raster
mqh_list <- c("ASFL")

# Get resistance from habitat models
# NOTE: Code below needs to be generalized before implementing loop 
#for(spp in mqh_list){
  
# Load habitat suitability map for the short-eared owl (Asio flammeus)
spp_mqh <- st_read(
  dsn = file.path(spatial_dir, "MQH", "MQH_HibouDesMarais_SHP"), 
  layer = "MQH_HibouDesMarais_s"
)

# Mask to modelling extent
mqh_model <- model_mask_shp(spp_mqh) %>%
  st_transform(target_crs)

# Rasterize habitat suitability classes
# NOTE: 0 = Hors zone, 1 = Non propice, 2 = Peu propice, 3 = Propice, and
#       4 = Très propice
mqh_raster <- rasterize(
  mqh_model, template_raster, field = "GRIDCODE")

# Keep suitable classes and set rank = 0
habitat_rank <- classify(
  mqh_raster,
  rcl = matrix(c(0, NA,
                 1, NA,
                 2, NA,
                 3, 0,
                 4, 0), ncol = 2, byrow = TRUE))

# Get LULC classes that overlap with suitable habitat
lulc_habitat <- mask(lulc_raster, habitat_rank)

# Reclass overlap of grassland/open wetlands and habitat to top rank
habitat_grassland <- ifel(lulc_habitat == 7, 0, NA)
habitat_wet_open <- ifel(lulc_habitat == 31, 0, NA)

# Combine and convert to resistance (resistance = 2^rank)
habitat_res_model <- 2^merge(habitat_grassland, habitat_wet_open)

#}



# Resistance -------------------------------------------------------------------

# Load LULC to resistance crosswalk
res_reclass_table <- read.csv(
  file.path(tabular_dir, "Crosswalks", "lulc_res_reclass.csv"),
  header = TRUE
)

# Load age penalty crosswalk
age_penalty_table <- read.csv(
  file.path(tabular_dir, "Crosswalks", "age_res_penalty.csv"),
  header = TRUE
)

# Species list
spp_list <- c("ASFL", "DOOR", "DRPI")

# Get resistance layer for each species
for(spp in spp_list){
  
  # Get column for target species
  spp_id <- which(colnames(res_reclass_table) == spp)
  
  # Reclassify LULC map into resistance
  res_model <- classify(
    x = lulc_raster, 
    rcl = res_reclass_table[,c(1, spp_id)]
    )
  
  # If habitat model is available for target species
  if(spp %in% mqh_list){
    # Complement resistance layer with habitat information
    res_model <- merge(habitat_res_model, res_model)
  }
  
  # If age penalty applies to species
  if(spp %in% names(age_penalty_table)){
    
    # Get column for target species in age penalty table
    spp_age_id <- which(colnames(age_penalty_table) == spp)
    
    # Reclassify age class map into penalties
    age_penalty_raster <- classify(
      x = age_raster, 
      rcl = age_penalty_table[,c(1, spp_age_id)]
    )
    
    # Set NAs to zero
    age_penalty_raster[is.na(age_penalty_raster)] <- 0
    
    # Convert resistance back to ranks
    ranks_model <- log2(res_model)
    
    # Add penalties to ranks
    ranks_penalty_model <- ranks_model + age_penalty_raster
    
    # Convert ranks back to resistance
    res_model <- 2^ranks_penalty_model
    
  }
  
  # Write to file
  writeRaster(
    x = res_model, 
    filename = file.path(output_dir, paste0("res_model_", spp, ".tif")),
    datatype = "INT4S",
    overwrite = TRUE)
  
}
