## Combine PA polygons into source raster
##
## Written by Carina Rauen Firkowski
## 
## Description:
## This script loads and combined the PAs layers created as intermediate outputs
## in the script "3-other-preprocess". Remaining aquatic PAs are removed, 
## polygon with overlapping boundaries are merged, and neighboring PAs less than
## 50 m apart are joined. Then, a filtering process is applied to keep only PAs
## larger than a size threshold based on the distribution of PA sizes.
##
## Outputs:
## - tif, filtered PA layer with IUCN categories



# Setup environment ------------------------------------------------------------ 

# Load constants
source(file.path("Scripts", "Preprocessing", "0-constants.R"))



# Load data --------------------------------------------------------------------

# Ontario PAs
on_fed_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "on_fed_pa_model"
  )
on_prov_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "on_prov_pa_model"
  )

# New Brunswick PAs
nb_wild_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "nb_wild_pa_model"
  )
nb_leg_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "nb_leg_pa_model"
  )
nb_nat_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "nb_nat_pa_model"
  )

# Quebec PAs
qc_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "qc_pa_model"
  )
qc_pa_rmn_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "qc_pa_rmn_model"
)

# U.S. PAs
ny_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "NY_pa_model"
)
vt_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "VT_pa_model"
)
nh_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
  layer = "NH_pa_model"
)
me_pa_model <- st_read(
  dsn = file.path(intermediate_dir, "PAs"),
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
         file.path(intermediate_dir, "PAs", "pa_merged_model.shp"),
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
st_write(pa_no_water,
         file.path(intermediate_dir, "PAs", "pa_no_water.shp"),
         append = FALSE)

# Dissolve boundaries of overlapping polygons
pa_dsslv <- pa_no_water %>%
  st_union() %>%
  st_as_sf() %>%
  st_cast("POLYGON") %>%
  st_set_geometry("geometry")

# Write to file
# st_write(pa_dsslv,
#          file.path(intermediate_dir, "PAs", "pa_dsslv.shp"),
#          append = FALSE)

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
# st_write(pa_model,
#          file.path(intermediate_dir, "PAs", "pa_model.shp"),
#          append = FALSE)

# Filter small polygons
# NOTE: Cut-off was determined by visually inspecting the histogram (see below).
pa_model_filter <- pa_model %>%
  filter(area >= (5000*1000))

# Write to file
# st_write(pa_model_filter,
#          file.path(output_dir, "pa_model_filter_500ha.shp"),
#          append = FALSE)

# Reattach attributes
pa_model_iucn <- st_intersection(pa_no_water, pa_model_filter) %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON")

# Tidy IUCN categories
pa_model_iucn$iucn[pa_model_iucn$iucn == "Y"] <- "Unassigned" 
pa_model_iucn$iucn[pa_model_iucn$iucn == "M"] <- "Multiple" 

# Write to file
st_write(pa_model_iucn, 
         file.path(output_dir, "pa_model_500ha_iucn.shp"),
         append = FALSE)



# Plot -------------------------------------------------------------------------

# Load data
#pa_model <- st_read(dsn = file.path(output_dir), layer = "pa_model")
# Get area from modelling extent
pa_area_df <- data.frame(area = pa_model$area, extent = "model")

# Sort area and get cumulative sum
pa_area_df <- pa_area_df %>%
  arrange(area) %>%
  mutate(cumsum = cumsum(area / sum(area)) )

# Load study area
study_area <- st_read(dsn = file.path(output_dir), layer = "study_area") %>%
  st_transform(target_crs)

# Mask to study area
pa_study <- st_intersection(pa_model, study_area)

# Sort area and get cumulative sum
pa_study <- pa_study %>%
  arrange(area) %>%
  mutate(cumsum = cumsum(area / sum(area)) )

# Setup dataframe
pa_study_df <- data.frame(
  area = pa_study$area,
  extent = "study",
  cumsum = pa_study$cumsum
  )

# Filter small polygons
pa_study_filter <- pa_study %>%
  filter(area >= (5000*1000))

# Get represented surface area
surf_area_study <- (sum(pa_study_df$area[pa_study_df$area >= 5000*1000])*100) / 
  sum(pa_study_df$area)
surf_area_model <- (sum(pa_area_df$area[pa_area_df$area >= 5000*1000])*100) /
  sum(pa_area_df$area)

# Get merge with area from study area
pa_areas_df <- rbind(pa_area_df, pa_study_df)

# Sort area
# PA area histogram at the modelling extent only
hist_model <- ggplot(pa_model, aes(x = log(area))) + 
  xlab("log(Area)") + ylab("Count") +
  geom_vline(aes(xintercept = log(5000*1000)), linetype = "dashed") +
  geom_histogram(color="black", fill = "steelblue", alpha = 0.2, 
                 linewidth = 0.4) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 14.5),
        axis.text = element_text(size = 14.5))

# PA area histogram at the study area extent
hist_study <- ggplot(pa_study, aes(x = log(area))) + 
  xlab("log(Area)") + ylab("Count") +
  geom_vline(aes(xintercept = log(5000*1000)), linetype = "dashed") +
  geom_histogram(color="black", fill = "steelblue", alpha = 0.2, 
                 linewidth = 0.4) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 14.5),
        axis.text = element_text(size = 14.5))

# PA area histogram of both extents
hist_both <- ggplot(pa_areas_df, aes(x = log(area), fill = extent)) + 
  labs(x = "log(Superficie)", y = "Fréquence") +
  geom_vline(aes(xintercept = log(5000*1000)), linetype = "dashed") +
  geom_histogram(color="black", alpha = 0.9, linewidth = 0.4) + 
  geom_line(mapping = aes(x = log(area), y = cumsum*3000, colour = extent),
            linewidth = 1.2) +
  scale_fill_manual(values=c("steelblue", "lightgreen"),
                    labels = c("Étendue de modélisation", "Zone d’étude")) +
  scale_colour_manual(values=c("steelblue", "lightgreen"),
                    labels = c("Étendue de modélisation", "Zone d’étude")) +
  scale_y_continuous(name = "Fréquence", 
                     sec.axis = sec_axis(~./3000,
                                         name = "Superficie cumulée")) +
  theme_classic() + 
  theme(axis.title = element_text(size = 14.5),
        axis.text = element_text(size = 14.5),
        legend.text = element_text(size = 14.5),
        legend.position = "bottom",
        legend.title = element_blank())

