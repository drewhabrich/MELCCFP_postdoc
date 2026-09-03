## ------------------------------- ##
## Script name: 01-PA_candidates
##
## Purpose of script: Prefilter node pairs by straight-line (euclidean) distance,
# per species' dispersal capacity. Cheap to compute, and avoids running
# expensive cost-distance on pairs that could never be within range --
# cost-distance is always >= euclidean distance when resistance >= 1.
#
# Saves one .csv per species to
# data/intermediate/candidate_pairs/{species}_PA_candidate_pairs.csv
# for two minimum PA sizes 10km2 and 5km2
#
# Also generates a set of candidate PA pairs using dispersal kernel estimation
# to determine 95% of *possible* movement. 
# Uses a negative exponential dispersal kernel to model probability of movement.
#
## Author: Andrew Habrich
##
## Notes ------------------------- ##
source(file.path("scripts", "setup_script.R"))

## 1. Load the protected area model shapefile ##################################
pa_mod <- vect(here(interm_dir, "pa_model_5km2.shp"),
               crs = target_crs)
# remove duplicated polygons (if any) based on geometry (may occur during merging)
#pa_clean <- tidyterra::distinct(pa_mod, geometry, .keep_all = T)
# save it to file for later use
#writeVector(pa_clean, here(interm_dir, "pa_model_clean.shp"), overwrite = T)

## Extract the centroid of each PA polygon as the nodes for the minimum spanning tree
nodes <- centroids(pa_mod, inside = T) #inside ensures the node is inside the polygon
nodes$site_id <- seq_len(nrow(nodes))

## Quick visual check of the centroids
plot(pa_mod)
plot(nodes, add = T, col = "red")
hist(nodes$area)   # distribution of patch areas (this should be the same as above)

## Save the centroid nodes to a vector for later
writeVector(nodes, here(output_dir, "pa_ctroidnodes.shp"), overwrite = TRUE)

## 2. Eligible PA pairs with distance < than the species-specific median distance #############
#read in the specieslist.xlsx
specieslist <- readxl::read_excel(here(data_dir, "specieslist.xlsx")) 

## Read the node centroids as an sf object
nodes <- read_sf(here(output_dir, "pa_ctroidnodes.shp"))

## Compute the distance matrix (in metres) between all centroids
coords <- st_coordinates(nodes) 
d <- as.matrix(dist(coords))
rownames(d) <- colnames(d) <- nodes$site_id
## Set the upper triangle to 0 to remove duplicate pairs.
d[upper.tri(d, diag = TRUE)] <- NA ## Because the matrix is symmetric, we only need the lower triangle.

## This is the distance between EACH PA centroid in the dataset
d_long <- as_tibble(d, rownames = "from") %>%
  pivot_longer(-from, names_to = "to", values_to = "euclidean_dist") %>%
  filter(!is.na(euclidean_dist))

## Create a loop to do the calculations for all the species with their respective dispersal distance
## NOTE: If the dispersal distance is too short, no eligible pairs will be found. MIN. PA-PAIR DISTANCE = 4.125 km
for (i in seq_len(nrow(specieslist))) {
  sp <- specieslist$species[i]
  disp <- specieslist$dispersal_km[i]*1000 #in meters
  ## filter to PA pairs that are within the dispersal distances
  candidates <- d_long %>% filter(euclidean_dist < disp)
  ## Add a column for the cost‑distance (to be filled later by the cost‑distance script)
  candidates$cost_dist <- NA
  write.csv(candidates, here(interm_dir, paste0(sp, "_PApairs_dist_cutoff.csv")), row.names = F)
  
  message(sp, ": ", nrow(candidates), " candidate pair(s) within ",
          disp, "m")
}

## 3. Eligible PA pairs using dispersal kernel cutoff ##########################
## Alternative approach to finding PA pairs using a dispersal kernel to estimate 95% of *potential* movement
#read in the specieslist.xlsx
specieslist <- readxl::read_excel(here(data_dir, "specieslist.xlsx")) 

## Read the node centroids as an sf object
nodes <- read_sf(here(interm_dir, "pa_ctroidnodes.shp"))

#### MODIFY FROM HERE
kernel_params <- list()

for (i in seq_len(nrow(specieslist))) {
  sp <- specieslist$species[i]
  #d_ref is the dispersal_km of the species, derived from the median literature value
  #p_ref describes what proportion of dispersers the d_ref refers to (median = 50%)
  alpha <- dispersal_kernel_alpha(
    d_ref = specieslist$dispersal_km[i]*1000,
    p_ref = 0.5
  )
  # p_threshold describes the probability cutoff we are looking to describe as potential movement
  p_threshold <- 0.05 
  
  # distance at which the kernel probability equals p_threshold -- the
  # search radius for this species' prefilter
  search_radius <- -log(p_threshold) / alpha
  
  kernel_params[[sp]] <- tibble(
    species             = sp,
    kernel_ref_dist_m = specieslist$dispersal_km[i]*1000,
    kernel_ref_prob     = 0.5,
    alpha               = alpha,
    kernel_p_threshold  = p_threshold,
    search_radius       = search_radius
  )
  
  candidates <- d_long %>%
    filter(euclidean_dist <= search_radius) %>%
    mutate(euclidean_kernel_prob = dispersal_probability(euclidean_dist, alpha))
  
  out_path <- file.path(interm_dir, paste0(sp, "_PApairs_dkernel_cutoff.csv"))
  write_csv(candidates, out_path)
  
  message(sp, ": alpha = ", signif(alpha, 4), ", search radius = ",
          round(search_radius), " map units -> ", nrow(candidates),
          " candidate pair(s) -> ", out_path)
  
  if (nrow(candidates) == 0) {
    warning(sp, ": no candidate pairs at all, even under the kernel's ",
            "long-distance tail (p >= ", p_threshold, ") -- every node is ",
            "isolated for this species at this site spacing. Consider ",
            "lowering kernel_p_threshold if that's ecologically defensible.")
  }
}

kernel_params_tbl <- bind_rows(kernel_params)
write_csv(kernel_params_tbl, file.path(interm_dir, "dispersal_kernel_params.csv"))

