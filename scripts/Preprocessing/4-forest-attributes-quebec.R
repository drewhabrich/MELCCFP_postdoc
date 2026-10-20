## Pre-process forest attributes
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This script loads each ecoforestry tile for Quebec created in script 
## "3-forest-attributes" and creates a raster tile for the target attributes  
## age, cover type, and density. Then for each target attribute, it combines its
## raster tiles into a single tile
##
## Intermediate outputs:
## - tif, age raster tiles for Quebec
## - tif, forest cover type tiles for Quebec
## - tif, forest density tiles for Quebec
## - tif, age raster for Quebec
## - tif, forest cover type for Quebec
## - tif, forest density for Quebec



# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))


# Rasterize tiles --------------------------------------------------------------

# List of forest attributes to rasterize
attr_list <- c("cl_age", "type_couv", "cl_dens")

# Load crosswalk tables
age_reclass <- read.csv(
  file.path(tabular_dir, "Crosswalks", "age_reclass.csv")
  )
dens_reclass <- read.csv(
  file.path(tabular_dir, "Crosswalks", "density_reclass.csv")
  )

# List all gpkg files
gpkg_files <- list.files(
  path = file.path(intermediate_dir, "Ecoforestry"),
  pattern = ".gpkg"
)

# Rasterize tiles
for(gpkg_file in gpkg_files){

  # Load data
  ecoforestry_tile <- st_read(
    dsn = file.path(intermediate_dir, "Ecoforestry", gpkg_file)
    ) %>%
    st_transform(target_crs)
  # Load data
  # ecoforestry_tile <- st_read(
  #   dsn = file.path(
  #     intermediate_dir, 
  #     "Ecoforestry", 
  #     paste0("ecoforestry_tile_", i, ".gpkg"))
  #   ) %>%
  #   st_transform(target_crs)
    
  # Get tile id
  i <- str_sub(gpkg_file, -7, -6)
  
  # If tile is empty, skip
  if(dim(ecoforestry_tile)[1] == 0){
    next
  } else {
    # For each attribute
    for(j in attr_list){
      
      # Select attribute and drop NA geometries
      ecoforestry_tile_attr <- ecoforestry_tile %>%
        select(all_of(j)) %>%
        drop_na(j)
      
      # "Reclassify" categories
      if(j == "cl_age"){
        # Merge with crosswalk
        ecoforestry_tile_int <- merge(
          ecoforestry_tile_attr, 
          age_reclass, 
          by = "cl_age"
          )
        # Assign reclass to attribute column
        ecoforestry_tile_int <- ecoforestry_tile_int %>%
          select(all_of("Code")) %>%
          rename(cl_age = Code)
      }
      if(j == "type_couv"){
        ecoforestry_tile_int <- ecoforestry_tile_attr %>%
          mutate(type_couv = case_when(
            type_couv %in% c("F") ~ 1,
            type_couv %in% c("M") ~ 2,
            type_couv %in% c("R") ~ 3)
          )
      }
      if(j == "cl_dens"){
        # Merge with crosswalk
        ecoforestry_tile_int <- merge(
          ecoforestry_tile_attr, 
          dens_reclass, 
          by = "cl_dens"
        )
        # Assign reclass to attribute column
        ecoforestry_tile_int <- ecoforestry_tile_int %>%
          select(all_of("Code")) %>%
          rename(cl_dens = Code)
      }
      
      # Rasterize and save tile
      get_tile_raster(shp_file = ecoforestry_tile_int, tile_id = i, attr = j)
    }
  }
}



# Combine ----------------------------------------------------------------------

# Load modelling extent
model_area <- st_read(dsn = file.path(output_dir), layer = "modelling_extent")

for(attr in attr_list[3]){
  
  # List all raster files for attribute
  tiles_list <- list.files(
    path = file.path(intermediate_dir, "Ecoforestry"),
    pattern = paste0(attr, ".tif$")
    )
  
  # Set loop counter
  loop_id <- 1
  
  # For each tile
  for(tile_id in tiles_list){
    
    # Load data
    tile_rast <- rast(file.path(intermediate_dir, "Ecoforestry", tile_id))
    
    # Combine
    if(loop_id == 1){
      attr_model <- tile_rast
    } else {
      attr_model <- merge(attr_model, tile_rast)
    }
    
    # Update loop count
    loop_id <- loop_id + 1
  }
  
  # Mask to modelling extent
  attr_model_mask <- mask(attr_model, model_area)
  
  # Write to file
  writeRaster(
    x = attr_model_mask, 
    filename = file.path(
      intermediate_dir, 
      "Ecoforestry", 
      paste0("ecoforestry_model_", attr, ".tif")
    ),
    datatype = "INT4S",
    overwrite = TRUE)
  
}

#plot(attr_model)



