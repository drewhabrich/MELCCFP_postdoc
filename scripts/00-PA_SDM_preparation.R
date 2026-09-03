## ------------------------------- ##
## Script name: 00-PA_SDM_preparation.R
## Purpose of script: Setup for working with the PAs network form QC and study extent. 
## This is a copy of the file on the github page, modified for the new directory.
#
#  This script combines all the PA file sources into a geospatial dataset and 
#  generates two subsets of PAs, based on area (>5km2, >10km2)
#
#  This script also processes the SDMs from Eckert et al. for each species to fit the study extent

## Author: Andrew Habrich
##
## Notes ------------------------- ##
source(file.path("scripts", "setup_script.R"))

# PAs --------------------------------------------------------------------------
# Mask LULC layer to modelling
model_mask_shp <- function(x) {
  
  # Load modelling extent and transform to layer CRS
  model_area <- st_read(
    dsn = file.path(data_dir), 
    layer = "modelling_extent"
  ) %>%
    st_transform(crs(x))
  
  # Crop & Mask
  lulc_masked <- x %>%
    st_intersection(model_area)
  
  return(lulc_masked)
}
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
get_us_state_pa <- function(state){
  
  # Unzip KMZ file into temporary directory
  unzip(
    file.path(
      data_dir, 
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
      dplyr::select(geometry) %>%
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
           file.path(interm_dir, "PAs", paste0(state, "_pa_attr.shp")),
           append = FALSE)
  
  # Select target PA types
  state_pa_select <- state_pa_attr %>%
    filter(d_GAP_Sts != "3 - managed for multiple uses - subject to extractive (e.g. mining or logging) or OHV use" &
             d_GAP_Sts != "4 - no known mandate for biodiversity protection")
  
  # Select and rename attributes
  state_pa_model <- state_pa_select %>%
    st_transform(target_crs) %>%
    dplyr::select(id = FID,
                  name = Unit_Nm,
                  type = d_Des_Tp,
                  iucn = IUCN_Cat)
  
  # Set FID to character
  state_pa_model$id <- as.character(state_pa_model$id)
  
  # Write to file
  st_write(state_pa_model, 
           file.path(interm_dir, "PAs", paste0(state, "_pa_model.shp")),
           append = FALSE)
  
  return(state_pa_model)
}

## Ontario federal PAs -----------------------------------------

# Load file
on_fed_pa <- st_read(
  dsn = file.path(data_dir, "PAs", "Federal_Protected_Area"), 
  layer = "Federal_Protected_Area"
) %>%
  st_make_valid() 

# Mask to modelling extent
on_fed_pa_mask <- model_mask_shp(on_fed_pa)

# Write to file
st_write(on_fed_pa_mask, 
         file.path(interm_dir, "PAs", "on_fed_pa_mask.shp"),
         append = FALSE)

# Select and rename attributes
on_fed_pa_model <- on_fed_pa_mask %>%
  st_transform(target_crs) %>%
  dplyr::select(id = OGF_ID,
         name = OFFICIAL_N,
         type = PROTECTED_)

# Set ID to character
on_fed_pa_model$id <- as.character(on_fed_pa_model$id)

# Add missing IUCN column
on_fed_pa_model$iucn <- "Unknown"

# Write to file
st_write(on_fed_pa_model, 
         file.path(interm_dir, "PAs", "on_fed_pa_model.shp"),
         append = FALSE)


## Ontario provincial PAs --------------------------------------

# Load file
on_prov_pa <- st_read(
  dsn = file.path(data_dir, "PAs", "Provincial_park_regulated"), 
  layer = "Provincial_park_regulated"
) %>%
  st_make_valid() 

# Mask to modelling extent
on_prov_pa_mask <- model_mask_shp(on_prov_pa)

# Write to file
st_write(on_prov_pa_mask, 
         file.path(interm_dir, "PAs", "on_prov_pa_mask.shp"),
         append = FALSE)

# Exclude select PA types
on_prov_pa_filter <- on_prov_pa_mask %>%
  filter(PROVINCIAL != "Waterway" &
           PROVINCIAL != "Cultural Heritage" 
  )

# Select and rename attributes
on_prov_pa_model <- on_prov_pa_filter %>%
  st_transform(target_crs) %>%
  dplyr::select(id = OGF_ID,
         name = PROTECTE_1,
         type = TYPE_ENG,
         iucn = IUCN_CATEG)

# Set ID to character
on_prov_pa_model$id <- as.character(on_prov_pa_model$id)

# Write to file
st_write(on_prov_pa_model, 
         file.path(interm_dir, "PAs", "on_prov_pa_model.shp"),
         append = FALSE)


## New Brunswick wildlife refuges ------------------------------

# Load file
nb_wild_pa <- st_read(
  dsn = file.path(data_dir, "PAs", "Wildlife_Refuge"), 
  layer = "Wildlife_Refuge"
)

# Mask to modelling extent
nb_wild_pa_mask <- model_mask_shp(nb_wild_pa)

# Select and rename attributes
nb_wild_pa_model <- nb_wild_pa_mask %>%
  st_transform(target_crs) %>%
  dplyr::select(id = OBJECTID,
         name = NAME,
         type = W_CLASS)

# Set ID to character
nb_wild_pa_model$id <- as.character(nb_wild_pa_model$id)

# Add missing IUCN category
nb_wild_pa_model$iucn <- "Unknown"

# Write to file
st_write(nb_wild_pa_model, 
         file.path(interm_dir, "PAs", "nb_wild_pa_model.shp"),
         append = FALSE)


## New Brunswick nature legacy PAs -----------------------------

# Load file
nb_leg_pa <- st_read(
  dsn = file.path(data_dir, "PAs", "Nature Legacy protected areas _ Les Aires protégées de l initiative Patrimoine naturel_20260521"), 
  layer = "geo_export_4695a750-664a-437e-a09f-e9593c4f9b07"
) %>%
  st_make_valid() 

# Mask to modelling extent
nb_leg_pa_mask <- model_mask_shp(nb_leg_pa)

# Select and rename attributes
nb_leg_pa_model <- nb_leg_pa_mask %>%
  st_transform(target_crs) %>%
  dplyr::select(id = objectid,
         name = name) %>%
  mutate(type = "Nature Legacy")

# Set ID to character
nb_leg_pa_model$id <- as.character(nb_leg_pa_model$id)

# Add missing IUCN column
nb_leg_pa_model$iucn <- "Unknown"

# Write to file
st_write(nb_leg_pa_model, 
         file.path(interm_dir, "PAs", "nb_leg_pa_model.shp"),
         append = FALSE)


## New Brunswick protected natural areas -----------------------

# Load file
nb_nat_pa <- st_read(
  dsn = file.path(data_dir, "PAs", "Protected Natural Areas _ Zones naturelles protégées_20260521"), 
  layer = "geo_export_65e5ab37-706b-4eda-9e3d-9cebf617eb51"
) %>%
  st_make_valid() 

# Mask to modelling extent
nb_nat_pa_mask <- model_mask_shp(nb_nat_pa)

# Select and rename attributes
nb_nat_pa_model <- nb_nat_pa_mask %>%
  st_transform(target_crs) %>%
  dplyr::select(id = objectid,
         name = name,
         type = class,
         iucn = iucn_categ)

# Set ID to character
nb_nat_pa_model$id <- as.character(nb_nat_pa_model$id)

# Write to file
st_write(nb_nat_pa_model, 
         file.path(interm_dir, "PAs", "nb_nat_pa_model.shp"),
         append = FALSE)


## Quebec open PAs ---------------------------------------------

# Load file
qc_pa <- st_read(
  dsn = file.path(data_dir, "PAs", "registre_aires_prot"), 
  layer = "AP_REG_S"
)

# Mask to modelling extent
qc_pa_mask <- model_mask_shp(qc_pa)

# Write to file
st_write(qc_pa_mask, 
         file.path(interm_dir, "PAs", "qc_pa_mask.shp"),
         append = FALSE)

# Exclude select PA types
qc_pa_filter <- qc_pa_mask %>%
  filter(DESIGNOM != "Réserve marine" &
           DESIGNOM != "Aire de concentration d'oiseaux aquatiques" &
           DESIGNOM != "Parc marin du Saguenay - Saint-Laurent"
  )

# Select and rename attributes
qc_pa_model <- qc_pa_filter %>%
  st_transform(target_crs) %>%
  dplyr::select(id = MACODE,
         name = TOPONYME,
         type = DESIGNOM,
         iucn = UICN)

# Set ID to character
qc_pa_model$id <- as.character(qc_pa_model$id)

# Write to file
st_write(qc_pa_model, 
         file.path(interm_dir, "PAs", "qc_pa_model.shp"),
         append = FALSE)


## Quebec RMN PAs ----------------------------------------------

# Load file (note: this file has some invalid geometries)
qc_pa_rmn <- st_read(
  dsn = file.path(data_dir, "PAs", "ProtectedConservedArea_2025", "ProtectedConservedArea_2025.gdb"), 
  layer = "ProtectedConservedArea_2025") |> 
  filter(BIOME == "T") |> 
  filter(!TYPE_E %in% c("Water fowl gathering area", "Proposed aquatic reserve", "Marine reserve",
                        "To Be Determined", "National Marine Park"))
  
# 1. Find the problem geometries
bad_idx <- which(sf::st_geometry_type(qc_pa_rmn, by_geometry = TRUE) == "MULTISURFACE")
# 2. Round-trip through Shapefile format to force GDAL to linearize curves into MULTIPOLYGON
tmp <- tempfile(fileext = ".shp")
sf::st_write(qc_pa_rmn[bad_idx, ], tmp, quiet = TRUE)
bad_lin <- sf::st_read(tmp, quiet = TRUE)
# 3. Swap the linearized geometries back into qc_pa_rmn (keep original attributes)
geom <- sf::st_geometry(qc_pa_rmn)
geom[bad_idx] <- sf::st_geometry(bad_lin)
sf::st_geometry(qc_pa_rmn) <- geom
# 4. Now st_make_valid() works on the whole object
qc_pa_rmn <- sf::st_make_valid(qc_pa_rmn)

# Mask to modelling extent
qc_pa_rmn_mask <- model_mask_shp(qc_pa_rmn)

# Write to file
st_write(qc_pa_rmn_mask, 
         file.path(interm_dir, "PAs", "qc_pa_rmn_mask.shp"),
         append = FALSE)

# Select and rename attributes
qc_pa_rmn_model <- qc_pa_rmn_mask %>%
  st_transform(target_crs) %>%
  dplyr::select(id = PARENT_ID,
         name = NAME_E,
         type = TYPE_E)

# Set ID to character
qc_pa_rmn_model$id <- as.character(qc_pa_rmn_model$id)

# Add missing IUCN column
qc_pa_rmn_model$iucn <- "Unknown"

# Write to file
st_write(qc_pa_rmn_model, 
         file.path(interm_dir, "PAs", "qc_pa_rmn_model.shp"),
         append = FALSE)

## U.S. state PAs ----------------------------------------------
# Get and save shapefile at modelling extent with target attributes
get_us_state_pa(state = "NY")
get_us_state_pa(state = "VT")
get_us_state_pa(state = "NH")
get_us_state_pa(state = "ME")

# Merge all the individual PA shapefiles together ###########################
# Load data --------------------------------------------------------------------

# Ontario PAs
on_fed_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "on_fed_pa_model"
)
on_prov_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "on_prov_pa_model"
)

# New Brunswick PAs
nb_wild_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "nb_wild_pa_model"
)
nb_leg_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "nb_leg_pa_model"
)
nb_nat_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "nb_nat_pa_model"
)

# Quebec PAs
qc_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "qc_pa_model"
)
qc_pa_rmn_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "qc_pa_rmn_model"
)

# U.S. PAs
ny_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "NY_pa_model"
)
vt_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "VT_pa_model"
)
nh_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "NH_pa_model"
)
me_pa_model <- st_read(
  dsn = file.path(interm_dir, "PAs"),
  layer = "ME_pa_model"
)

# Combine ----------------------------------------------------------------------

# Combine all PA shapefiles
pa_merged_model <- bind_rows(
  qc_pa_model, 
  qc_pa_rmn_model,
  on_prov_pa_model, 
  on_fed_pa_model,
  nb_wild_pa_model, 
  nb_leg_pa_model, 
  nb_nat_pa_model,
  ny_pa_model,
  nh_pa_model,
  me_pa_model,
  vt_pa_model
)

# Write to file
st_write(pa_merged_model,
         file.path(interm_dir, "PAs", "pa_merged_model.shp"),
         append = FALSE)

# Remove large aquatic PAs
pa_no_water <- pa_merged_model %>%
  filter(id != "42" &
           id != "667" &
           id != "5048" &
           id != "6525" &
           id != "6741" &
           id != "166653" &
           id != "166654" &
           id != "166655" &
           id != "166657" &
           id != "166658" &
           id != "166663" &
           id != "166664" &
           id != "166665" &
           id != "166666" &
           id != "166667" &
           id != "166668" &
           id != "167224" &
           id != "7086" &
           id != "167092")

# Write to file
pa_no_water <- st_write(pa_no_water,
         file.path(interm_dir, "PAs", "pa_no_water.shp"),
         append = FALSE)

# Dissolve boundaries of overlapping polygons
pa_dsslv <- pa_no_water %>%
  st_union() %>%
  st_as_sf() %>%
  st_cast("POLYGON") %>%
  st_set_geometry("geometry")

# Write to file
st_write(pa_dsslv,
         file.path(interm_dir, "PAs", "pa_dsslv.shp"),
         append = FALSE)

# Set buffer for joining neighboring PAs
buffer_dist <- 50

# Join neighboring PAs
pa_model <- pa_dsslv %>%
  st_buffer(dist = buffer_dist) %>%
  st_union() %>%
  st_as_sf() %>%
  st_cast("POLYGON") %>%
  st_set_geometry("geometry") %>%
  st_buffer(dist = -buffer_dist)

# Calculate area of polygons
pa_model$area <- as.numeric(st_area(pa_model))

# Write to file
st_write(pa_model, 
         file.path(interm_dir, "PAs", "pa_model.shp"),
         append = FALSE)

### Filter small polygons #############################
##
# pa_no_water <- st_read(file.path(interm_dir, "PAs", "pa_no_water.shp"))
# pa_model <- st_read(file.path(interm_dir, "PAs", "pa_model.shp"))
##
# NOTE: Cut-off was determined by visually inspecting the histogram
pa_model_10km2 <- pa_model %>%
  filter(area >= (10000*1000)) #10km2 filter
pa_model_5km2 <- pa_model %>%
  filter(area >= (5000*1000)) #5km2 filter

# Write to file
st_write(pa_model_10km2, 
         file.path(interm_dir, "pa_model_10km2.shp"),
         append = FALSE)
st_write(pa_model_5km2, 
         file.path(interm_dir, "pa_model_5km2.shp"),
         append = FALSE)

# Reattach attributes
# Use a spatial join instead of st_intersection so the pa_model_filter
# polygons stay intact (no splitting/duplicating geometry for each
# overlapping input). Where multiple pa_no_water polygons intersect the
# same pa_model_filter polygon, keep only the attributes of the first match.
# pa_model_iucn <- pa_model_filter %>%
#   mutate(feat_id = row_number()) %>%
#   st_join(pa_no_water, join = st_intersects, left = TRUE) %>%
#   group_by(feat_id) %>%
#   slice(1) %>%
#   ungroup() %>%
#   select(-feat_id) %>%
#   st_make_valid()
# 
# # Tidy IUCN categories
# pa_model_iucn$iucn[pa_model_iucn$iucn == "Y"] <- "Unassigned" 
# pa_model_iucn$iucn[pa_model_iucn$iucn == "M"] <- "Multiple" 
# 
# # Write to file
# st_write(pa_model_iucn, 
#          file.path(interm_dir, "pa_model_iucn.shp"),
#          append = FALSE)
# 
# # Generate a histogram of distribution of area in the pa_model object
# ggplot(pa_model, aes(x = log(area))) + 
#   xlab("log(Area)") + ylab("Count") +
#   geom_vline(aes(xintercept = log(5000*1000)), linetype = "dashed") +
#   geom_density(color="black", fill = "steelblue", alpha = 0.2, 
#                  linewidth = 0.4) + 
#   theme_classic() + 
#   theme(axis.title = element_text(size = 14.5),
#         axis.text = element_text(size = 14.5))

# SDMs ---------------------------------------------------------------------------
# Crop + reproject continent-wide SDMs to the modelling extent (one-time)
# SDMs in data/drive_data/SDMs/ are continent-wide, 6-layer in their own CRS 
# select just the "current" layer, crop (in the SDM's own CRS, cheap) to modelling_extent, then reproject 

sdm_dir <- here(data_dir, "SDMs")
sdm_cropped_dir <- here(interm_dir, "SDM_cropped")
if (!dir.exists(sdm_cropped_dir)) dir.create(sdm_cropped_dir)

modelling_extent <- vect(here(data_dir, "modelling_extent.shp"))

sdm_files <- list.files(sdm_dir, pattern = "\\.tif$", full.names = TRUE)
for (sdm_file in sdm_files) {
  out_file <- file.path(sdm_cropped_dir, basename(sdm_file))
  if (file.exists(out_file)) {
    message("Already cropped: ", basename(sdm_file), " -- skipping.")
    next
  }
  message("Cropping/reprojecting SDM: ", basename(sdm_file))

  sdm <- rast(sdm_file)
  if (!"current" %in% names(sdm)) {
    warning(basename(sdm_file), ": no 'current' layer (layers: ",
            paste(names(sdm), collapse = ", "), ") -- skipping.")
    next
  }
  
  #This layer is the T0 SDM from Eckert et al. (not projected forward for climate scenarios)
  sdm <- sdm[["current"]] 

  # fill NA gaps (e.g. reprojection edge cells, or areas the SDM has no
  # prediction for) with 0 -- treated as "no suitability" rather than left
  # undefined, so script 07's corridor x SDM multiplication doesn't punch
  # NA holes into the overlap where corridor density is otherwise valid.
  sdm[is.na(sdm)] <- 0
  
  # crop in the SDM's own CRS first (cheap), using modelling_extent
  # reprojected into that CRS and -- a rectangle in
  # one CRS doesn't map to a rectangle in another, so the buffer has to be
  # generous enough that the true corresponding region stays inside it --
  # then reproject the now-small piece onto the project's CRS.
  extent_in_sdm_crs <- project(modelling_extent, crs(sdm))
  sdm <- crop(sdm, ext(extent_in_sdm_crs) + 10000)
  sdm <- project(sdm, target_crs)

  writeRaster(sdm, out_file, overwrite = TRUE)
  message("  -> saved ", out_file)
}