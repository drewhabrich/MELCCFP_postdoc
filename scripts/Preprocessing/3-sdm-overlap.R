## SDM overlap across selected species
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This script loads the masked SDMs prepared in script "3-sdm-mask", creates 
## binary maps by thresholding probability values >= 0.6, and sums all species
## to obtain a map of SDM overlap expressed in number of species. Lastly, it 
## calculates area covered by species' distributions. 
##
## Outputs:
## - tif, binary SDM maps for each selected species
## - tif, SDM overlap across all selected species



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



# SDM overlap ------------------------------------------------------------------

# Selected species
spp_list <- c("Aquila chrysaetos", "Castor canadensis", "Clethrionomys gapperi",
              "Dolichonyx oryzivorus", "Hylatomus pileatus", "Lynx canadensis",
              "Martes americana", "Plethodon cinereus",
              "Asio flammeus")

# Loop through species list
for(i in 1:length(spp_list)){
  
  # Get species name
  species_name <- spp_list[i]
  
  # Load species' masked SDM
  spp_sdm <- rast(
    file.path(intermediate_dir, "SDMs", paste0(species_name, "_model.tif"))
    )  
  
  # Threshold current SDM at 0.6 suitability
  spp_sdm_current <- classify(
    x = spp_sdm$current,
    rcl = matrix(c(0, 0.6, 0,
                   0.6, 1, 1), nrow = 2, byrow = TRUE),
    include.lowest = TRUE
  )
  
  # Write to file
  writeRaster(
    x = spp_sdm_current,
    filename = file.path(
      intermediate_dir, 
      "SDMs", 
      paste0(species_name, "_binary.tif")
    ),
    overwrite = TRUE,
    datatype = "FLT4S"
  )
  
  # Combine results
  if(i == 1){
    sdm_overlap <- spp_sdm_current 
  } else {
    sdm_overlap <- sdm_overlap + spp_sdm_current
  }
  
}
 
# Resample to template raster
sdm_overlap_model <- resample(sdm_overlap, template_raster, method = "near")

# Load wood turtle binary map
sdm_turtle_raster <- rast(
  file.path(intermediate_dir, "SDMs", "Glyptemys insculpta_model.tif"))

# Set NAs to 0
sdm_turtle_raster[is.na(sdm_turtle_raster)] <- 0

# Add wood turtle 
sdm_overlap_spp <- sdm_overlap_model + sdm_turtle_raster

# Mask to modelling extent area
spp_sdm_model <- mask(sdm_overlap_spp, model_area)
# Mask to study area
spp_sdm_study <- mask(sdm_overlap_spp, study_area)

# Write to file
writeRaster(spp_sdm_model, file.path(output_dir, "sdm_overlap_model.tif"),
            overwrite = TRUE)
# Write to file
writeRaster(spp_sdm_study, file.path(output_dir, "sdm_overlap_study.tif"), 
            overwrite = TRUE)



# Coverage ---------------------------------------------------------------------

# Frequency in study area
sdm_study_freq <- freq(spp_sdm_study)

# Coverage min. 1 species
study_coverage_1 <- sum(sdm_study_freq$count[sdm_study_freq$value != 0]) /
  sum(sdm_study_freq$count) * 100
# Coverage min. 2 species
study_coverage_2 <- (sum(sdm_study_freq$count[sdm_study_freq$value != 0 &
                                                sdm_study_freq$value != 1])) /
  sum(sdm_study_freq$count) * 100


