## ------------------------------- ##
## Script name: setup_script
##
## Purpose of script: Setup for workflow, including loading scripts, file locations etc.
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##

## 1. Load relevant packages--------
# General utility
library(tidyverse)
library(dplyr)
library(rvest)
library(here)
library(patchwork)

# Spatial data manipulation
library(terra)
library(tidyterra)
library(sf)
library(sfheaders)
library(smoothr)
library(spatialEco)
library(gdistance)

## 2. File paths and directories --------
# Same project-root layout as Scripts/Preprocessing/0-constants.R, so the
# analysis reads the preprocessing outputs from where they are written.
data_dir <- here("Data")                          # raw inputs
spatial_dir <- here(data_dir, "Spatial")          # raw spatial data (PAs, SDMs, ...)
tabular_dir <- here(data_dir, "Tabular")          # tabular inputs (species list, lookup tables, species tables)
intermediate_dir <- here("Intermediate")
output_dir <- here("Output")                      # preprocessing outputs + numbered analysis outputs
fig_dir <- here(output_dir, "Figures")
table_dir <- here(output_dir, "Tables")

if (!dir.exists(intermediate_dir)) { dir.create(intermediate_dir) }
if (!dir.exists(output_dir)) { dir.create(output_dir) }

# Per-stage folders, each prefixed with the number of the script that writes
# it. Defined once here because later stages read them too -- the writing
# script creates its own folder(s).
pa_dir          <- here(intermediate_dir, "00_PAs")                # 00: harmonized PA polygons
sdm_cropped_dir <- here(intermediate_dir, "00_SDM_cropped")        # 00: SDMs cropped to the modelling extent
pa_pairs_dir    <- here(intermediate_dir, "01_PA_candidates")      # 01: {sp}_PApairs_dkernel_cutoff.csv / _edge_lines.gpkg
pa_patches_dir  <- here(output_dir, "01_PA_patches")               # 01: pa_patches.gpkg, pa_ctroidnodes.shp
resistance_dir  <- here(output_dir, "02_resistance")               # 02: {sp}_res.tif
costdist_dir    <- here(intermediate_dir, "03_cost_distances")     # 03: {sp}_PA_costs.csv, {sp}_costdist_profiling.csv
mst_dir         <- here(intermediate_dir, "04_mst")                # 04: {sp}_mst_edges.csv, {sp}_mst.rds
mst_lines_dir   <- here(output_dir, "04_mst")                      # 04: {sp}_mst_lines.gpkg
corridor_dir    <- here(intermediate_dir, "05_lcp_corridors")      # 05: {sp}_lcp_corridors.gpkg
density_dir     <- here(intermediate_dir, "06_corridor_density")   # 06: {sp}_corridor_density.tif, all_species_corridor_density.tif
sdm_overlap_dir <- here(output_dir, "07_mst_SDM_overlap")          # 07: {sp}_habitat_overlap.tif
sensitivity_dir <- here(intermediate_dir, "08_sensitivity")        # 08: {sp}/{scenario}_result.csv, _mst_edges.csv, _overlap_1km.tif
sensitivity_out_dir <- here(output_dir, "09_sensitivity")          # 09: {sp}_mst_edge_retention.gpkg, {sp}_overlap_variability.tif

## 3. Spatial ---------------------------
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

#---- negative exponential dispersal kernel helpers ----
# Shared by *PA candidate prefilter* and *MST edge cutoff*,
# so both use exactly the same kernel definition.
#
# p(d) = exp(-alpha * d); calibrated from a known (distance, probability)
# reference point via alpha = -log(p_ref) / d_ref. Same convention as
# graph4lg::kernel_param(mode = "A") if you want to cross-check.
dispersal_kernel_alpha <- function(p_ref, d_ref) {
    -log(p_ref) / d_ref
  }

dispersal_probability <- function(dist, alpha) {
  exp(-alpha * dist)
}

# sentinel value used to mark origin cells for terra::costDist() -- must not
# occur naturally in any resistance raster
sentinel_value <- -1

#---- species-specific corridor-density smoothing bandwidth (script 06) ----
# sigma = the MEAN of the species' negative-exponential dispersal kernel
# (median_dispersal / ln(2), since dispersal_kernel_alpha() is calibrated at
# p_ref = 0.5 i.e. the median) -- deliberately larger than the median
# dispersal distance itself, so the smoothed corridor surface reflects a
# species' broader potential-movement footprint, not just its single
# cheapest least-cost path. Species with dispersal_km > cap_km are treated
# as cap_km (continuous with `sparsify_dispersal_km` in script 01, which
# sparsifies candidate pairs the same way). Floored at 2 cells so
# near-sessile species still get a usable smoothing kernel.
#
# LITERATURE BASIS: this is not an arbitrary inflation of dispersal_km.
# Compton, B.W., McGarigal, K., Cushman, S.A., Gamble, L.R. (2007). "A
# Resistant-Kernel Model of Connectivity for Amphibians that Breed in
# Vernal Pools." Conservation Biology 21(3):788-799 -- the founding paper
# for resistant-kernel connectivity modeling (the same family of technique
# script 06 implements: kernel-smoothed density around dispersal paths,
# not a flat buffer) -- define their kernel's smoothing bandwidth h as the
# STANDARD DEVIATION of a curve fit to the species' actual dispersal-
# distance data (a normal curve there; h = 399.6 m for marbled
# salamanders). A separate, larger parameter -- the kernel's outer
# radius/truncation threshold -- is set to the species' MAXIMUM dispersal
# distance in the related UNICOR/resistant-kernel literature (Landguth et
# al.; Cushman et al.), not the bandwidth itself; that convention is what
# the cap_km/sparsify_dispersal_km treatment above separately mirrors.
# For THIS project's kernel family (negative exponential, not normal),
# mean = standard deviation = 1/alpha = median/ln(2) -- a basic property
# of the exponential distribution. So this formula computes exactly "1 SD
# of the fitted dispersal kernel," i.e. Compton et al.'s own definition of
# kernel bandwidth, just for a different parametric family. sigma
# exceeding the median dispersal distance is a mathematical consequence of
# SD > median for an exponential, not a separate inflation choice.
corridor_sigma <- function(dispersal_km, cap_km = 100,
                            min_sigma = 2 * res(template_raster)[1]) {
  d_m <- min(dispersal_km, cap_km) * 1000
  max(1 / dispersal_kernel_alpha(p_ref = 0.5, d_ref = d_m), min_sigma)
}

## 4. Shared pipeline parameters -----------
# Method settings used by more than one script (the stage scripts AND the
# sensitivity analysis, 08/09) -- defined once here so every script runs
# with the same values. Pure compute knobs (worker counts, batch sizes)
# stay in each script.

# Kernel probability defining the dispersal cutoff (01 prefilter, 04 MST
# edges): pairs are kept while p(d) >= this, i.e. out to a search radius
# of dispersal_km * log2(1 / kernel_p_threshold) (~4.3 x at 0.05).
kernel_p_threshold <- 0.05

# 01: species with dispersal_km above this get their candidate pairs
# sparsified to each PA's `sparsify_knn` nearest partners (edge-to-edge)
# plus the Euclidean-MST edges (see 01's header).
sparsify_dispersal_km <- 100
sparsify_knn <- 20

# 03: how far beyond an origin's candidate-destination bounding box to crop
# the resistance raster, as a fraction of that origin's farthest candidate
# distance (plus a minimum, so origins with few/close destinations still
# get a sane amount of room). Doubled once automatically if any
# destination comes back unreachable within the crop.
cd_buffer_frac <- 0.25
cd_buffer_min_m <- 1000

# 05: how far beyond each edge's straight-line bounding box to crop the
# resistance raster, as a fraction of that edge's straight-line distance
# (plus a minimum, so very short edges still get a sane amount of room).
# Doubled once automatically if the first attempt fails to find a path.
# Shrinking these is a free way to speed up tracing further, at some risk
# of more edges needing the doubled-buffer retry.
lcp_buffer_frac <- 0.25
lcp_buffer_min_m <- 1000
lcp_buffer_max_m <- 5000  # absolute cap regardless of edge length -- keeps very long edges' corridors from also being very wide

# 06: how finely each species' own coarse grid resolves relative to its
# sigma (grid cell ~= sigma / cells_per_sigma): smaller = more detail
# retained, larger = faster/coarser. Cells are capped at common_res_m.
cells_per_sigma <- 10

# 06/08/09: common grid resolution (m) for combining/comparing surfaces
# across species or runs -- 1 km matches the SDMs' native resolution. Also
# the coarsest grid 06 smooths corridor density on (so 07 overlaps are <= 1 km).
common_res_m <- 1000

# 07: fraction of the corridor density raster's lowest values to treat as
# noise (set to NA before rescaling), e.g. 0.1 = drop the bottom 10%.
# Adapted from the worked example 07 was based on, not a validated
# threshold -- check whether it's dropping meaningful low-density
# corridor or genuinely negligible noise for your data before trusting it.
quantile_cutoff <- 0.1

# Max forest-age penalty (in ranks) added to a forest type the species rates
# "Unsuitable". Age can't make an already-unsuitable forest much worse, so
# Unsuitable forest tops out at Unsuitable rank + 1 (= urban's rank, 6).
# Applied in build_resistance(); see forest_unsuitable_rank there.
unsuitable_age_penalty_cap <- 1

## 5. Per-species pipeline functions -----------
# The core of every stage (01-07) as functions, shared with the
# sensitivity analysis -- see pipeline_functions.R.
source(file.path("Scripts", "Analysis", "pipeline_functions.R"))
