## Add roads to LULC layer
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This script prepares the road layer and merges it with the LULC map. Data for
## each jurisdiction is pre-processed separately and combined at the end. The 
## hierarchy of road overlap is: railway (top), major roads, minor roads, 
## forestry resource roads (bottom).
##
## Outputs:
## - tif, road raster
## - tif, LULC map with roads



# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))

# Road buffer size to align with LULC resolution
road_buffer <- 30



# Quebec -----------------------------------------------------------------------

# Reclass table
qc_roads_reclass <- read.csv(
  file.path(tabular_dir, "Crosswalks", "qc_road_reclass.csv"), header = TRUE
)
qc_roads_forest_reclass <- read.csv(
  file.path(tabular_dir, "Crosswalks", "qc_road_forest_reclass.csv"), header = TRUE
)


## Major and minor roads ----------------------------------------

# List of datasets to combine
data_list <- c("AQreseau_SHP", "ESRI(SHP)")

for(dataset in data_list){
  
  # Load data
  if(dataset == "AQreseau_SHP"){
    qc_roads <- st_read(
      dsn = file.path(spatial_dir, "Roads", "AQreseau_SHP"), 
      layer = "AQ_ROUTES"
    )
  }
  if(dataset == "ESRI(SHP)"){
    qc_roads <- st_read(
      dsn = file.path(spatial_dir, "Roads", "ESRI(SHP)"), 
      layer = "Reseau_routier"
    )
  }
  
  # Mask to modelling extent
  qc_roads_mask <- model_mask_shp(qc_roads) %>%
    st_transform(target_crs)
  rm(qc_roads)
  
  # Write to file
  st_write(qc_roads_mask, 
           file.path(intermediate_dir, 
                     "Roads", 
                     paste0("qc_roads_mask_", dataset, ".shp")),
           append = FALSE)
  
  # Merge reclass table with shapefile
  qc_roads_class <- merge(
    qc_roads_mask, qc_roads_reclass, by.x = "ClsRte", by.y = "qc_class"
  )
  rm(qc_roads_mask)
  
  # Keep target attributes
  qc_roads_target <- qc_roads_class %>%
    select(IdRte, ClsRte, class_name, class_id, geometry)
  
  # Add Cls_CheFor attribute
  if(data_list == "AQreseau_SHP"){
    qc_roads_target$Cls_CheFor <- NA
  }
  if(data_list == "ESRI(SHP)"){
    qc_roads_target <- cbind(qc_roads_target, qc_roads_class$Cls_CheFor)
    names(qc_roads_target)[5] <- "Cls_CheFor"
  }
  rm(qc_roads_class)
  
  # Combine shapefiles
  if (dataset == data_list[1]) {
    qc_all_roads_class <- qc_roads_target
  } else {
    qc_all_roads_class <- rbind(qc_all_roads_class, qc_roads_target)
  }
  rm(qc_roads_target)
  
}

# Write to file
st_write(qc_all_roads_class, 
         file.path(intermediate_dir, "Roads", "qc_all_roads_class.shp"),
         append = FALSE)

# Separate road types
for(road_type in c("Minor", "Major")){
  
  # Filter target road type
  qc_all_roads_target <- qc_all_roads_class %>%
    filter(class_name == road_type)

  # Get number of rows in dataset
  num_rows <- dim(qc_all_roads_target)[1]
  
  # Set number of breaks
  num_breaks <- 3
  
  # Get breaks
  row_breaks <- seq(1, num_rows, by = num_rows/num_breaks)
  
  for (i in 1:num_breaks) {
    
    cat(i)
    
    # Define starting and ending row
    start_row <- row_breaks[i]
    end_row <- row_breaks[i] + num_rows/num_breaks-1
    
    # Add buffer to roads
    qc_all_roads_buffer <- qc_all_roads_target[start_row:end_row,] %>%
      st_buffer(road_buffer)
    
    # Separate and rasterize roads
    qc_roads_sub_raster <- qc_all_roads_buffer %>%
      rasterize(template_raster, field = "class_id")
    
    # Write to file
    writeRaster(
      x = qc_roads_sub_raster, 
      filename = file.path(intermediate_dir, 
                           "Roads", 
                           paste0("qc_road_raster_", i, "_", road_type,".tif")),
      datatype = "INT4S",
      overwrite = TRUE)
  
  }
}


## Forestry roads -----------------------------------------------

# For the forest resource roads, break down classification further

# Filter target road type
qc_all_roads_target <- qc_all_roads_class %>%
  filter(class_name == "Resource")

# Merge reclass table with shapefile
qc_all_roads_forest <- merge(
  qc_all_roads_target, qc_roads_forest_reclass, 
  by.x = "Cls_CheFor", by.y = "qc_class"
  )

# For each forest resource road class
for(res_class in c("Resource A", "Resource B", "Resource C")){
  
  # Filter target road type
  qc_all_roads_class <- qc_all_roads_forest %>%
    filter(res_class_name == res_class)
  
  # Add buffer to roads
  qc_all_roads_class_buffer <- qc_all_roads_class %>%
    st_buffer(road_buffer)
  
  # Separate and rasterize roads
  qc_roads_sub_raster <- qc_all_roads_class_buffer %>%
    rasterize(template_raster, field = "res_class_id")
  
  # Write to file
  writeRaster(
    x = qc_roads_sub_raster, 
    filename = file.path(intermediate_dir, 
                         "Roads", 
                         paste0("qc_road_raster_", res_class,".tif")),
    datatype = "INT4S",
    overwrite = TRUE)
  
}
 

## Railway ------------------------------------------------------

# Load dataset
qc_rail <- st_read(
  dsn = file.path(spatial_dir, "Roads", "ESRI(SHP)"), 
  layer = "Reseau_ferroviaire"
)

# Mask to modelling extent
qc_rail_mask <- model_mask_shp(qc_rail) %>%
  st_transform(target_crs)
rm(qc_rail)

# Remove ferry lines
qc_rail_filter <- qc_rail_mask %>%
  filter(IdSegChfer != 69 &
           IdSegChfer != 12124)

# Add buffer to roads
qc_rail_buffer <- qc_rail_filter %>%
  st_buffer(road_buffer)

# Rasterize roads
qc_rail_raster <- rasterize(
  qc_rail_buffer, template_raster, field = ""
)

# Set class ID
qc_rail_raster <- qc_rail_raster * 6

# Write to file
writeRaster(
  x = qc_rail_raster, 
  filename = file.path(intermediate_dir, 
                       "Roads", 
                       paste0("qc_rail_raster.tif")),
  datatype = "INT4S",
  overwrite = TRUE)



# Ontario ----------------------------------------------------------------------

# Load data
on_roads <- st_read(
  dsn = file.path(spatial_dir, "Roads", "LIO-2026-06-09"), 
  layer = "ORN_ROAD_NET_ELEMENT"
)

# Mask to modelling extent
on_roads_mask <- model_mask_shp(on_roads) %>%
  st_transform(target_crs)
rm(on_roads)

# Load road class information
on_roads_class_df <- read.csv(
  file.path(spatial_dir, "Roads", "LIO-2026-06-09", "ORN_ROAD_CLASS.csv")
  ) %>%
  select(ORN_ROAD_NET_ELEMENT_ID, ROAD_CLASS)

# Merge with shapefile
on_roads_class <- merge(
  on_roads_mask, on_roads_class_df, 
  by.x = "OGF_ID", by.y = "ORN_ROAD_NET_ELEMENT_ID"
  )

# Reclass table
on_roads_reclass <- read.csv(
  file.path(tabular_dir, "Crosswalks", "on_road_reclass.csv")
  )

# Merge with shapefile
on_roads_class <- merge(
  on_roads_class, on_roads_reclass, by.x = "ROAD_CLASS", by.y = "on_class"
)

# Add buffer to roads
on_roads_buffer <- on_roads_class %>%
  st_buffer(road_buffer)

# Separate and rasterize major roads
on_major_raster <- on_roads_buffer %>%
  filter(class_name == "Major") %>%
  rasterize(template_raster, field = "class_id")

# Separate and rasterize minor roads
on_minor_raster <- on_roads_buffer %>%
  filter(class_name == "Minor") %>%
  rasterize(template_raster, field = "class_id")

# Combine roads
on_road_raster <- merge(on_major_raster, on_minor_raster)

# Write to file
writeRaster(
  x = on_road_raster, 
  filename = file.path(intermediate_dir, "Roads", "on_road_raster.tif"),
  datatype = "INT4S",
  overwrite = TRUE)



# New Brunswick ----------------------------------------------------------------

# Load data
nb_roads <- st_read(
  dsn = file.path(spatial_dir, "Roads", "geonb_nbrn-rrnb_shp"), 
  layer = "geonb_nbrn-rrnb_road-route"
)

# Mask to modelling extent
nb_roads_mask <- model_mask_shp(nb_roads) %>%
  st_transform(target_crs)
rm(nb_roads)

# Reclass table
nb_roads_reclass <- read.csv(
  file.path(tabular_dir, "Crosswalks", "nb_road_reclass.csv")
)

# Merge with shapefile
nb_roads_class <- merge(
  nb_roads_mask, nb_roads_reclass, by.x = "ROADCLASS", by.y = "nb_class"
)

# Add buffer to roads
nb_roads_buffer <- nb_roads_class %>%
  st_buffer(road_buffer)

# Separate and rasterize major roads
nb_major_raster <- nb_roads_buffer %>%
  filter(class_name == "Major") %>%
  rasterize(template_raster, field = "class_id")

# Separate and rasterize minor roads
nb_minor_raster <- nb_roads_buffer %>%
  filter(class_name == "Minor") %>%
  rasterize(template_raster, field = "class_id")

# Combine roads
nb_road_raster <- merge(nb_major_raster, nb_minor_raster)

# Write to file
writeRaster(
  x = nb_road_raster, 
  filename = file.path(intermediate_dir, "Roads", "nb_road_raster.tif"),
  datatype = "INT4S",
  overwrite = TRUE)



# U.S. -------------------------------------------------------------------------

# List all folders
folder_list <- basename(list.dirs(
  path = file.path(spatial_dir, "Roads", "US"),
  recursive = FALSE
))

# Loop through all folder to combine shapefiles
for(folder in folder_list){
  
  # Load data
  county_roads <- st_read(
    dsn = file.path(
      spatial_dir, 
      "Roads", "US",
      folder
      ), 
    layer = folder
  )
  
  # Mask to modelling extent
  county_roads_mask <- model_mask_shp(county_roads) %>%
    st_transform(target_crs)
  
  # Combine results
  if(folder == folder_list[1]){
    county_model <- county_roads_mask
  } else {
    county_model <- rbind(county_model, county_roads_mask)
  }
}

# Reclass table
us_roads_reclass <- read.csv(
  file.path(tabular_dir, "Crosswalks", "us_road_reclass.csv")
)

# Merge with shapefile
county_roads_class <- merge(
  county_model, us_roads_reclass, by.x = "MTFCC", by.y = "us_id"
)

# Add buffer to roads
county_roads_buffer <- county_roads_class %>%
  st_buffer(road_buffer)

# Separate and rasterize major roads
us_major_raster <- county_roads_buffer %>%
  filter(class_name == "Major") %>%
  rasterize(template_raster, field = "class_id")

# Separate and rasterize minor roads
us_minor_raster <- county_roads_buffer %>%
  filter(class_name == "Minor") %>%
  rasterize(template_raster, field = "class_id")

# Combine roads
county_road_raster <- merge(us_major_raster, us_minor_raster)

# Write to file
writeRaster(
  x = county_road_raster, 
  filename = file.path(intermediate_dir, "Roads", "county_road_raster.tif"),
  datatype = "INT4S",
  overwrite = TRUE)



# Combine ----------------------------------------------------------------------

## Roads --------------------------------------------------------

# Start loop counter
loop_id <- 1

# Combine resource roads
# NOTE: The res_class loop order matters as it dictates the overlay hierarchy
for(res_class in c("Resource C", "Resource B", "Resource A")){
  
  # Load raster
  qc_road_raster <- rast(
    file.path(intermediate_dir, "Roads", paste0("qc_road_raster_", res_class,".tif"))
  )
  
  # Merge
  if (loop_id == 1) {
    qc_road_all_raster <- qc_road_raster
  } else {
    qc_road_all_raster <- merge(qc_road_raster, qc_road_all_raster)
  }
  
  # Update loop counter
  loop_id <- loop_id + 1 
  
}

# Add Quebec's major and minor roads
# NOTE: The road_type loop order matters as it dictates the overlay hierarchy
for(road_type in c("Minor", "Major")){
  for (i in 1:num_breaks) {
  
    # Load raster
    qc_road_raster <- rast(
      file.path(intermediate_dir, "Roads", 
                paste0("qc_road_raster_", i, "_", road_type,".tif"))
    )
    
    qc_road_all_raster <- merge(qc_road_raster, qc_road_all_raster)
  
  }
}

# Write to file
writeRaster(
  x = qc_road_all_raster, 
  filename = file.path(intermediate_dir, "Roads", "qc_road_all_raster.tif"),
  datatype = "INT4S",
  overwrite = TRUE)

# Add other jurisdictions 
road_raster <- merge(
  qc_rail_raster,
  qc_road_all_raster,
  on_road_raster, 
  nb_road_raster, 
  county_road_raster
  )

# Write to file
writeRaster(
  x = road_raster, 
  filename = file.path(output_dir, "road_raster.tif"),
  datatype = "INT4S",
  overwrite = TRUE)


## LULC ---------------------------------------------------------

# Load data
lulc_model <- rast(file.path(intermediate_dir, "LULC", "lulc_model.tif"))

# Merge
lulc_raster <- merge(road_raster*100, lulc_model)

# Set 0 pixels around modelling extent boundary to NA
lulc_raster[lulc_raster == 0] <- NA

# Write to file
writeRaster(
  x = lulc_raster, 
  filename = file.path(output_dir, "lulc_road_model.tif"),
  datatype = "INT4S",
  overwrite = TRUE)
