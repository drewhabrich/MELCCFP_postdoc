## Prepare study area
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This code loads the Canadian provinces and coastal water shapefiles
## to retain the boundaries for Quebec and the St. Lawrence river. Then, it
## creates the parallel line for latitude 49, which requires working in CRS 4326
## (geodetic coordinate system for world). Next, it defines the buffer size for
## the modelling extent, based on the northern limit of the St. Lawrence river.
## Then, the strict study area is created by defining the spatial extent of the 
## mainland plus St. Lawrence Lowlands south of the 49th parallel. The extent is
## used to clip the Quebec province to define the strict study area. The 
## modelling extent is created by slicing the mainland plus St. Lawrence 
## Lowlands at the 49th parallel plus buffer, adding a buffer around all edges, 
## and slicing the polygon again at the 49th parallel plus buffer. Finally, the 
## study area (which acts as a clipping mask) is created by slicing the mainland
## and St. Lawrence river polygons after adding to the 49th parallel line a 
## small increment northwards to ensure it include the entire Gaspé peninsula.
##
## Outputs:
## - shp, 49th parallel
## - shp, strict study area
## - shp, modelling extent
## - shp, study area (clipping mask)



# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))



# Load data --------------------------------------------------------------------

# Provincial boundaries
cad <- st_read(
  dsn = file.path(spatial_dir, "Canada"),
  layer = "lpr_000b16a_e"
  ) %>%
  st_transform(crs = 4326)

# Coastal waterbodies
water_coast <- st_read(
  dsn = file.path(spatial_dir, "Coastal waterbodies"),
  layer = "ghy_000h11a_e"
) %>%
  st_transform(crs = 4326) %>%
  st_make_valid()



# Study area -------------------------------------------------------------------

## Pre-process spatial layers -----------------------------------

# Get Quebec boundary
qc <- cad[cad$PRENAME == "Quebec",]

# Cast single multipolygon to polygons
qc_poly <- st_cast(qc, "POLYGON")

# Calculate area of polygons & sort
qc_poly$area <- st_area(qc_poly)
area_temp <- sort(qc_poly$area, decreasing = T)

# Calculate extent
qc_ext <- ext(qc_poly)

# Get St. Lawrence river and golf boundary
st_lwrnc <- water_coast[water_coast$HYDROUID == "6858881" | 
                          water_coast$HYDROUID == "7764631",] %>%
  summarize()
rm(water_coast)


## Latitude lines & Buffer size ---------------------------------

# Get northern most extent of St. Lawrence river
ext_n <- round(as.numeric(st_bbox(st_lwrnc)$ymax), 0)

# Calculate buffer size 
# (i.e., difference between 49th and northern extent limit of St. Lawrence)
# In latitude
buffer_lat <- round(abs(49 - ext_n), 1)
# In degrees
buffer_d <- round(terra::distance(
  matrix(c(qc_ext[1], ext_n, qc_ext[1], 49), byrow = TRUE, ncol = 2), 
  lonlat = TRUE), 0)[1]

# Set coordinates for 49th parallel line with buffer
line49b_coords <- st_linestring(
  matrix(c(-180, 49 + buffer_lat, 180, 49 + buffer_lat), ncol = 2, byrow = TRUE)
)
# Set coordinates for 49th parallel line without buffer
line49_coords <- st_linestring(
  matrix(c(-180, 49, 180, 49), ncol = 2, byrow = TRUE)
)

# Convert coordinates to sf object
lat49b_line <- st_sf(geometry = st_sfc(line49b_coords), crs = 4326)
lat49_line <- st_sf(geometry = st_sfc(line49_coords), crs = 4326)

# Densify 49th parallel line
lat49_line_dense <- smoothr::densify(lat49_line, max_distance = 0.1)

# Write to file
st_write(lat49_line_dense, file.path(output_dir, "lat49_line.shp"), append = FALSE)


## Strict study area (south of 49th parallel) -------------------

# Keep mainland and St. Lawrence lowlands
qc_main <- qc_poly[qc_poly$area == area_temp[1] | qc_poly$area == area_temp[2],]

# Slice Quebec at 49th parallel
qc_slice49 <- lwgeom::st_split(qc_main, lat49_line) %>%
  st_collection_extract("POLYGON")

# Calculate area of polygons & sort
qc_slice49$area <- st_area(qc_slice49)
area_temp <- sort(qc_slice49$area, decreasing = T)

# Discard north of 49th parallel
qc_south49 <- qc_slice49[qc_slice49$area != area_temp[1] & 
                         qc_slice49$area != area_temp[4],]

# Define extent of study area
study_ext <- ext(qc_south49)

sf_use_s2(FALSE)
# Crop Quebec by study extent
study_area <- st_crop(qc_poly, study_ext)
sf_use_s2(TRUE)

# Densify 49th parallel line
study_area_dense <- smoothr::densify(study_area, max_distance = 0.1)

# Transform to target CRS
study_area <- study_area_dense %>%
  st_transform(target_crs)

# Write to file
st_write(
  study_area, file.path(output_dir, "strict_study_area.shp"), append = FALSE
  )

# Largest dimension of study area (m)
study_area_width <- terra::distance(
  matrix(c(study_ext[1], study_ext[3], study_ext[2], study_ext[3]), 
         byrow = TRUE, ncol = 2), 
  lonlat = TRUE)[1]

# Calculate buffer proportion
buffer_prop <- (buffer_d * 100) / study_area_width



# Modelling extent -------------------------------------------------------------

# Slice Quebec at 49th parallel + buffer
qc_slice49b <- lwgeom::st_split(qc_main, lat49b_line) %>%
  st_collection_extract("POLYGON")

# Calculate area of polygons & sort
qc_slice49b$area <- st_area(qc_slice49b)
area_temp <- sort(qc_slice49b$area, decreasing = T)

# Discard north of 49th parallel
qc_south49b <- qc_slice49b[qc_slice49b$area != area_temp[1],]

# Dissolve all polygons into one & add smooth buffer to boundary
qc_south49b_buffer <- st_union(qc_south49b) %>%
  st_sf() %>%
  st_buffer(dist = buffer_d, max_cells = 20000) 

# Slice again at northern limit
buffer_slice49b <- lwgeom::st_split(qc_south49b_buffer, lat49b_line) %>%
  st_collection_extract("POLYGON")

# Keep south portion
model_area <- buffer_slice49b[1,]

# Densify 49th parallel line
model_area_dense <- smoothr::densify(model_area, max_distance = 0.1)

# Transform to target CRS
model_area <- model_area_dense %>%
  st_transform(target_crs)

# Load Quebec municipalities
# qc_munic <- st_read(dsn = file.path(spatial_dir, "SDA"), layer = "munic_s") %>%
#   st_transform(target_crs)
# 
# # Keep Anticosti island
# qc_island <- qc_munic[qc_munic$MUS_NM_NMC == "Municipalité de L'Île-d'Anticosti",]
# rm(qc_munic)
# 
# # Remove Anticosti island from modelling extent
# model_area_no_island <- st_difference(model_area, st_geometry(qc_island))

# Write to file
st_write(model_area, 
         file.path(output_dir, "modelling_extent.shp"),
         append = FALSE)



# Study area (clipping mask) ---------------------------------------------------

# Get north most point of peninsula
clip_ymax <- round(as.numeric(st_bbox(qc_main$geometry[2])$ymax), 2) - 0.004

# Set coordinates for 49th parallel line plus Gaspé
line49g_coords <- st_linestring(
  matrix(c(-180, clip_ymax, 180, clip_ymax), ncol = 2, byrow = TRUE)
)

# Convert coordinates to sf object
lat49g_line <- st_sf(geometry = st_sfc(line49g_coords), crs = 4326)

# Slice Quebec at 49th parallel plus Gaspé
qc_slice49g <- lwgeom::st_split(qc_main, lat49g_line) %>%
  st_collection_extract("POLYGON")

# Calculate area of polygons & sort
qc_slice49g$area <- st_area(qc_slice49g)
area_temp <- sort(qc_slice49g$area, decreasing = T)

# Discard north of 49th parallel
qc_south49g <- qc_slice49g[qc_slice49g$area != area_temp[1],]

# Define extent of study area
clip_ext <- ext(qc_south49g)

sf_use_s2(FALSE)
# Crop Quebec by study extent
clip_area <- st_crop(qc_poly, clip_ext)
sf_use_s2(TRUE)

# Densify 49th parallel line
clip_area_dense <- smoothr::densify(clip_area, max_distance = 0.1)

# Slice waterbodies at 49th parallel plus Gaspé
st_lwrnc_slice49g <- lwgeom::st_split(st_lwrnc, lat49g_line) %>%
  st_collection_extract("POLYGON")

# Calculate area of polygons & sort
st_lwrnc_slice49g$area <- st_area(st_lwrnc_slice49g)
area_temp <- sort(st_lwrnc_slice49g$area, decreasing = T)

# Discard north of 49th parallel
st_lwrnc_49g <- st_lwrnc_slice49g[st_lwrnc_slice49g$area == area_temp[2],]

# Combine study area with waterbody, add small buffer (CRS transformation 
# required), and remove holes
clip_river <- clip_area_dense %>%
  st_union(st_lwrnc_49g) %>%
  st_transform(crs = target_crs) %>%
  st_buffer(dist = 1800) %>%
  st_union() %>%
  sf_remove_holes() %>%
  st_transform(crs = 4326)

# Slice again at northern limit
clip_river_slice <- lwgeom::st_split(clip_river, lat49g_line) %>%
  st_collection_extract("POLYGON")

# Keep south portion
clip_mask <- clip_river_slice[1,]

# Densify 49th parallel line
clip_mask_dense <- smoothr::densify(clip_mask, max_distance = 0.1)

# Transform to target CRS
clip_mask <- clip_mask_dense %>%
  st_transform(target_crs)

# Write to file
st_write(clip_mask, 
         file.path(output_dir, "study_area.shp"),
         append = FALSE)


