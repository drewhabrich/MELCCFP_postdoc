#TEST SCRIPT#
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))

## CHECKING COST DIST CALCULATION TIMES FROM FIR.ALLIANCE ##
library(ggpubr)
here()
time <- read.csv(here(costdist_dir, "ASFL_costdist_profiling.csv"))
glimpse(time)
## Plot n_cells against total_sec to see if there is a relationship
time |> 
  ggscatter("n_cells", "total_sec") +
  stat_cor()

## compare the mean times for each step
# 1. Reshape the wide data frame into a long format
time_l <- pivot_longer(time, cols = everything(), names_to = "step", values_to = "time")
glimpse(time_l)
# 2. Plot using ggplot
time_l |> filter(!step %in% c("n_cells", "n_dest","origin_id","buffer_mult")) |> 
ggplot(aes(x = step, y = time, fill = step)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Steps in Cost-distance calculation", y = "Time (s)")

#check the pa costs output to see if there are any obvious errors
asfl_cd <- read.csv(here(costdist_dir, "ASFL_PA_costs.csv"))
glimpse(asfl_cd)
ggdensity(asfl_cd$cost_dist)
ggscatter(asfl_cd, "euclidean_dist", "cost_dist")+
  stat_cor()

#how many unique origins are there? There should be 556 origin nodes
asfl_cd %>% 
  summarise(n_nodes = n_distinct(from))

#check raster values
door1 <- rast(here(resistance_dir, "DOOR_res.tif"))
door2 <- rast(here(output_dir, "res_model_DOOR.tif"))

### CHECK SDMS ###
maamsdm <- rast(here(spatial_dir, "SDMs", "Martes americana.tif"))
glimpse(maamsdm)
plot(maamsdm$current)
plot(maamsdm$current_clipped)
plot(maamsdm$current_thresholded)
global(maamsdm[["current_thresholded"]], "min", na.rm=TRUE)
