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

# Spatial data manipulation
library(terra)
library(tidyterra)
library(sf)
library(sfheaders)
library(smoothr)
library(spatialEco)
library(gdistance)

## 2. File paths and directories --------
data_dir <- here("data/drive_data/")
interm_dir <- here("data/intermediate/")
output_dir <- here("data/output_data/")
resist_dir <- here("data/output_data/resistance/")
fig_dir <- here("output/figures/")
table_dir <- here("output/tables/")

if (!dir.exists(interm_dir)) { dir.create(interm_dir) }
if (!dir.exists(output_dir)) { dir.create(output_dir) }

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
