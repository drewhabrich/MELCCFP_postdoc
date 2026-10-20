## Pre-process other data layers
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This script pre-processes the following data by masking it to the modelling 
## extent: Resistance or movement cost raster from Pither et al. 2023 (PLOS 
## One), digital elevation model for Canada and U.S., and protected area layers 
## across all jurisdictions.
## For the PAs, this script also harmonized the attributes across all layers, 
## and keeps only terrestrial protected areas
##
## Outputs:
## - tif, resistance raster from Pither et al.
## - tif, DEM 
##
## Intermediate outputs:
## - shp, Ontario federal PAs
## - shp, Ontario provincial PAs
## - shp, New Brunswick wildlife refuges
## - shp, New Brunswick nature legacy PAs
## - shp, New Brunswick protected natural areas
## - shp, Quebec open-source PAs
## - shp, Quebec RMN PAs
## - shp, Vermont PAs
## - shp, New York PAs
## - shp, New Hampshire PAs
## - shp, Maine PAs



# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))



# Resistance (Pither et al.) ---------------------------------------------------

# Load data
pither <- rast(file.path(spatial_dir, "Other", "Movement_Cost_Layer.tif"))

# Mask to modelling extent
pither_mask <- model_mask(pither)

# Modify template raster to match resolution of resistance layer
template_raster_res <- template_raster
res(template_raster_res) <- res(pither_mask)

# Project to template raster
pither_model <- project(pither_mask, template_raster_res, method = "near")

# Write to file
writeRaster(
  x = pither_model, 
  filename = file.path(output_dir, "pither_model.tif"),
  datatype = "INT4S",
  overwrite = TRUE)



# DEM --------------------------------------------------------------------------

## Canada ------------------------------------------------------

# Load data
dem_can <- rast(file.path(spatial_dir, "Other", "mrdem-30-dtm.tif"))

# Mask to modelling extent
dem_can_mask <- model_mask(dem_can)

# Write to file
# writeRaster(
#   x = dem_can_mask, 
#   filename = file.path(intermediate_dir, "dem_can_mask.tif"),
#   datatype = "FLT4S",
#   overwrite = TRUE)

# Project to template raster
dem_can_model <- project(dem_can_mask, template_raster)

# Write to file
# writeRaster(
#   x = dem_can_model, 
#   filename = file.path(intermediate_dir, "dem_can_model.tif"),
#   datatype = "FLT4S",
#   overwrite = TRUE)


## U.S. --------------------------------------------------------

# Load data
dem_us <- rast(file.path(spatial_dir, "Other", "output_USGS30m.tif"))

# Mask to modelling extent
dem_us_mask <- model_mask(dem_us)

# Project to template raster
dem_us_model <- project(dem_us_mask, template_raster)

# Write to file
# writeRaster(
#   x = dem_us_model, 
#   filename = file.path(intermediate_dir, "dem_us_model.tif"),
#   datatype = "FLT4S",
#   overwrite = TRUE)


## Combine -----------------------------------------------------

# Merge raster layers
dem_model <- merge(dem_can_model, dem_us_model)

# Write to file
writeRaster(
  x = dem_model, 
  filename = file.path(output_dir, "dem_model.tif"),
  datatype = "FLT4S",
  overwrite = TRUE)



# PAs --------------------------------------------------------------------------

## Ontario federal PAs -----------------------------------------

# Load file
on_fed_pa <- st_read(
  dsn = file.path(spatial_dir, "PAs", "Federal_Protected_Area"), 
  layer = "Federal_Protected_Area"
  ) %>%
  st_make_valid() 

# Mask to modelling extent
on_fed_pa_mask <- model_mask_shp(on_fed_pa)

# Write to file
st_write(on_fed_pa_mask, 
         file.path(intermediate_dir, "PAs", "on_fed_pa_mask.shp"),
         append = FALSE)

# Select and rename attributes
on_fed_pa_model <- on_fed_pa_mask %>%
  st_transform(target_crs) %>%
  select(id = OGF_ID,
         name = OFFICIAL_N,
         type = PROTECTED_)

# Set ID to character
on_fed_pa_model$id <- as.character(on_fed_pa_model$id)

# Add missing IUCN column
on_fed_pa_model$iucn <- "Unknown"

# Write to file
st_write(on_fed_pa_model, 
         file.path(intermediate_dir, "PAs", "on_fed_pa_model.shp"),
         append = FALSE)


## Ontario provincial PAs --------------------------------------

# Load file
on_prov_pa <- st_read(
  dsn = file.path(spatial_dir, "PAs", "Provincial_park_regulated"), 
  layer = "Provincial_park_regulated"
  ) %>%
  st_make_valid() 

# Mask to modelling extent
on_prov_pa_mask <- model_mask_shp(on_prov_pa)

# Write to file
st_write(on_prov_pa_mask, 
         file.path(intermediate_dir, "PAs", "on_prov_pa_mask.shp"),
         append = FALSE)

# Exclude select PA types
on_prov_pa_filter <- on_prov_pa_mask %>%
  filter(PROVINCIAL != "Waterway" &
          PROVINCIAL != "Cultural Heritage" 
  )

# Select and rename attributes
on_prov_pa_model <- on_prov_pa_filter %>%
  st_transform(target_crs) %>%
  select(id = OGF_ID,
         name = PROTECTE_1,
         type = TYPE_ENG,
         iucn = IUCN_CATEG)

# Set ID to character
on_prov_pa_model$id <- as.character(on_prov_pa_model$id)

# Write to file
st_write(on_prov_pa_model, 
         file.path(intermediate_dir, "PAs", "on_prov_pa_model.shp"),
         append = FALSE)


## New Brunswick wildlife refuges ------------------------------

# Load file
nb_wild_pa <- st_read(
  dsn = file.path(spatial_dir, "PAs", "Wildlife_Refuge"), 
  layer = "Wildlife_Refuge"
)

# Mask to modelling extent
nb_wild_pa_mask <- model_mask_shp(nb_wild_pa)

# Select and rename attributes
nb_wild_pa_model <- nb_wild_pa_mask %>%
  st_transform(target_crs) %>%
  select(id = OBJECTID,
         name = NAME,
         type = W_CLASS)

# Set ID to character
nb_wild_pa_model$id <- as.character(nb_wild_pa_model$id)

# Add missing IUCN category
nb_wild_pa_model$iucn <- "Unknown"

# Write to file
st_write(nb_wild_pa_model, 
         file.path(intermediate_dir, "PAs", "nb_wild_pa_model.shp"),
         append = FALSE)


## New Brunswick nature legacy PAs -----------------------------

# Load file
nb_leg_pa <- st_read(
  dsn = file.path(spatial_dir, "PAs", "Nature Legacy protected areas _ Les Aires protégées de l initiative Patrimoine naturel_20260521"), 
  layer = "geo_export_4695a750-664a-437e-a09f-e9593c4f9b07"
) %>%
  st_make_valid() 

# Mask to modelling extent
nb_leg_pa_mask <- model_mask_shp(nb_leg_pa)

# Select and rename attributes
nb_leg_pa_model <- nb_leg_pa_mask %>%
  st_transform(target_crs) %>%
  select(id = objectid,
         name = name) %>%
  mutate(type = "Nature Legacy")

# Set ID to character
nb_leg_pa_model$id <- as.character(nb_leg_pa_model$id)

# Add missing IUCN column
nb_leg_pa_model$iucn <- "Unknown"

# Write to file
st_write(nb_leg_pa_model, 
         file.path(intermediate_dir, "PAs", "nb_leg_pa_model.shp"),
         append = FALSE)


## New Brunswick protected natural areas -----------------------

# Load file
nb_nat_pa <- st_read(
  dsn = file.path(spatial_dir, "PAs", "Protected Natural Areas _ Zones naturelles protégées_20260521"), 
  layer = "geo_export_65e5ab37-706b-4eda-9e3d-9cebf617eb51"
) %>%
  st_make_valid() 

# Mask to modelling extent
nb_nat_pa_mask <- model_mask_shp(nb_nat_pa)

# Select and rename attributes
nb_nat_pa_model <- nb_nat_pa_mask %>%
  st_transform(target_crs) %>%
  select(id = objectid,
         name = name,
         type = class,
         iucn = iucn_categ)

# Set ID to character
nb_nat_pa_model$id <- as.character(nb_nat_pa_model$id)

# Write to file
st_write(nb_nat_pa_model, 
         file.path(intermediate_dir, "PAs", "nb_nat_pa_model.shp"),
         append = FALSE)


## Quebec open-sources PAs -------------------------------------

# Load file
qc_pa <- st_read(
  dsn = file.path(spatial_dir, "PAs", "registre_aires_prot"), 
  layer = "AP_REG_S"
  )

# Mask to modelling extent
qc_pa_mask <- model_mask_shp(qc_pa)

# Write to file
st_write(qc_pa_mask, 
         file.path(intermediate_dir, "PAs", "qc_pa_mask.shp"),
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
  select(id = MACODE,
         name = TOPONYME,
         type = DESIGNOM,
         iucn = UICN)

# Set ID to character
qc_pa_model$id <- as.character(qc_pa_model$id)

# Write to file
st_write(qc_pa_model, 
         file.path(intermediate_dir, "PAs", "qc_pa_model.shp"),
         append = FALSE)


## Quebec RMN PAs ----------------------------------------------

# Load file
qc_pa_rmn <- st_read(
  dsn = file.path(spatial_dir, "PAs", "Rep_Sites_proteges_202511"), 
  layer = "Rep_Sites_proteges_202511"
) %>%
  st_make_valid()

# Mask to modelling extent
qc_pa_rmn_mask <- model_mask_shp(qc_pa_rmn)

# Write to file
st_write(qc_pa_rmn_mask, 
         file.path(intermediate_dir, "PAs", "qc_pa_rmn_mask.shp"),
         append = FALSE)

# Select and rename attributes
qc_pa_rmn_model <- qc_pa_rmn_mask %>%
  st_transform(target_crs) %>%
  select(id = ID_SITE,
         name = Nom_Site,
         type = Ty_Prot)

# Set ID to character
qc_pa_rmn_model$id <- as.character(qc_pa_rmn_model$id)

# Add missing IUCN column
qc_pa_rmn_model$iucn <- "Unknown"

# Write to file
st_write(qc_pa_rmn_model, 
         file.path(intermediate_dir, "PAs", "qc_pa_rmn_model.shp"),
         append = FALSE)



## U.S. state PAs ----------------------------------------------

# Get and save shapefile at modelling extent with target attributes
get_us_state_pa(state = "NY")
get_us_state_pa(state = "VT")
get_us_state_pa(state = "NH")
get_us_state_pa(state = "ME")


