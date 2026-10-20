## Constants for MELCCP 2026-2027 repository
##
## Written by Carina Rauen Firkowski
##
## Description:
## Loads packages, paths and functions needed across scripts.


# Load packages ---------------------------------------------------------------- 

library(dplyr)
library(tidyr)
library(terra)
library(sf)
library(sfheaders)
library(ggplot2)
library(smoothr)
library(spatialEco)
library(rvest)
library(stars)



# Setup paths ------------------------------------------------------------------

# Simple paths
data_dir <- file.path("Data")
intermediate_dir <- file.path("Intermediate")
output_dir <- file.path("Output")

# Composite paths
spatial_dir <- file.path(data_dir, "Spatial")
tabular_dir <- file.path(data_dir, "Tabular")

# Create folder, if needed
if (!dir.exists(intermediate_dir)) { dir.create(intermediate_dir) }
if (!dir.exists(output_dir)) { dir.create(output_dir) }


# Spatial ----------------------------------------------------------------------

# CRS for Quebec
target_crs <- "EPSG:32198"

# Target resolution
target_res <- c(30, 30)

# Reference raster for grid
reference_ext <- ext(c(-830370, 429480, 117940, 744580)) # Quebec LULC
reference_raster <- rast(
  extent = reference_ext, 
  res = target_res, 
  crs = target_crs
) 

# Modelling extent
model_ext <- ext(c(
  -942562.735316926, 
  429467.232992499, 
  5971.82034726957, 
  744578.222380984
  ))

# Target extent
target_ext <- align(model_ext, reference_raster, snap = "near")
  
# Template raster
template_raster <- rast(
  extent = target_ext, 
  res = target_res, 
  crs = target_crs
) 



# Functions --------------------------------------------------------------------

# Mask LULC layer to modelling extent
model_mask <- function(x) {
  
  # Load modelling extent and transform to layer CRS
  model_area <- st_read(
    dsn = file.path(output_dir), 
    layer = "modelling_extent"
  ) %>%
    st_transform(crs(x))
  
  # Crop & Mask
  lulc_masked <- x %>%
    crop(model_area) %>%
    mask(model_area)
  
  return(lulc_masked)
}

# Mask LULC layer to modelling extent
buffer_mask <- function(x) {
  
  # Load modelling extent and transform to layer CRS
  model_area <- st_read(
    dsn = file.path(output_dir), 
    layer = "study_area"
  ) %>%
    st_transform(crs(x))
  
  # Crop & Mask
  lulc_masked <- x %>%
    mask(model_area, inverse = TRUE)
  
  return(lulc_masked)
}

# Mask LULC layer to modelling extent
model_mask_shp <- function(x) {
  
  # Load modelling extent and transform to layer CRS
  model_area <- st_read(
    dsn = file.path(output_dir), 
    layer = "modelling_extent"
  ) %>%
    st_transform(crs(x))
  
  # Crop & Mask
  lulc_masked <- x %>%
    st_intersection(model_area)
  
  return(lulc_masked)
}

# Convert the html table nested within the description attribute for the U.S. PA shapefiles
get_attributes <- function(us_pa, id){

  # Convert html text attribute into attribute datafram
  html_tb <- us_pa$Description[id] %>%
    minimal_html() %>%
    html_table()

  # Keep second table and pivot wide
  attribute_tb <- html_tb[[2]] %>%
    pivot_wider(names_from = X1, values_from = X2)
  
  return(attribute_tb)
}

# Open U.S. state PA FMZ file, mask to modelling extent, and build attribute table off of html table 
get_us_state_pa <- function(state){

  # Unzip KMZ file into temporary directory
  unzip(
  file.path(
    spatial_dir, 
    "PAs", 
    paste0("PADUS4_1_State_", state, "_GDB_KMZ"), 
    paste0("PADUS4_1Combined_State", state, ".kmz")
    ),
  exdir = tempdir()
  )

  # Get path to the extracted KML
  kml_path <- file.path(tempdir(), "doc.kml")

  # Check available layers
  st_layers(kml_path)

  # Load file
  state_pa <- st_read(
    dsn = kml_path,
    layer = "SHP_lyr"
  ) %>%
  st_make_valid() 

  # Mask to modelling extent
  state_pa_mask <- model_mask_shp(state_pa)

  # Loop through each polygon
  for(i in 1:dim(state_pa_mask)[1]){

    # Get attribute table
    attribute_tb <- get_attributes(state_pa_mask, id = i)

    # Merge attributes and geometry
    state_pa_temp <- state_pa_mask[i,] %>%
      select(geometry) %>%
      cbind(attribute_tb)

    # Combine results
    if(i == 1){
      state_pa_attr <- state_pa_temp
    } else {
      state_pa_attr <- rbind(state_pa_attr, state_pa_temp)
    }
  }

  # Write to file
  st_write(state_pa_attr, 
          file.path(intermediate_dir, "PAs", paste0(state, "_pa_attr.shp")),
          append = FALSE)

  # Select target PA types
  state_pa_select <- state_pa_attr %>%
    filter(d_GAP_Sts != "3 - managed for multiple uses - subject to extractive (e.g. mining or logging) or OHV use" &
             d_GAP_Sts != "4 - no known mandate for biodiversity protection")
  
  # Select and rename attributes
  state_pa_model <- state_pa_select %>%
    st_transform(target_crs) %>%
    select(id = FID,
           name = Unit_Nm,
           type = d_Des_Tp,
           iucn = IUCN_Cat)
  
  # Set FID to character
  state_pa_model$id <- as.character(state_pa_model$id)

  # Write to file
  st_write(state_pa_model, 
          file.path(intermediate_dir, "PAs", paste0(state, "_pa_model.shp")),
          append = FALSE)

  return(state_pa_model)
}

# Rasterizes tile based on a defined forest attribute and saves to file
get_tile_raster <- function(shp_file, tile_id, attr) {
  
  # Rasterize attribute
  attr_tile <- rasterize(shp_file, template_raster, field = attr)
  
  # Write to file
  writeRaster(
    x = attr_tile, 
    filename = file.path(
      intermediate_dir, 
      "Ecoforestry", 
      paste0("ecoforestry_tile_", tile_id, "_", attr, ".tif")
    ),
    datatype = "INT4S",
    overwrite = TRUE)
}



