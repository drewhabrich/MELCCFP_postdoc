## ------------------------------- ##
## Script name: fig_reachable_PAs
##
## Purpose of script: Slide figure showing, for one species, the full set of
## PAs vs. the subset of PAs that are "potentially reachable" -- i.e. have at
## least one candidate partner within that species' dispersal-kernel cutoff
## (Intermediate/01_PA_candidates/{species}_PApairs_dkernel_cutoff.csv, from script
## 01-PA_candidates.R). Two side-by-side panels: all PAs, and reachable PAs
## highlighted against the full set.
##
## Also produces a pairwise-connection-density figure: one straight line per
## candidate pair, edge-to-edge between the two PA polygons' nearest boundary
## points (sf::st_nearest_points(), same convention as scripts 01/03/05 --
## not centroids), rasterized and kernel-smoothed to show where dispersal-
## kernel-eligible connections are spatially dense vs. sparse. NOTE these are
## still straight-line distances, not resistance-weighted least-cost
## corridors -- it is not the corridor density from script 06 (which needs
## cost-distance rasters that don't exist for a species until stage 03/05
## have run).
##
## Also draws every individual candidate-pair line directly (network-style),
## so individual connections stay visible rather than only their smoothed
## density.
##
## Saves Output/Figures/{species}_reachable_PAs.png,
## Output/Figures/{species}_connection_density.png,
## Output/Figures/{species}_connection_network.png, and
## Output/{species}_PApairs_edge_lines.gpkg (the edge-to-edge lines
## themselves, for use in GIS)
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
source(file.path("Scripts", "Analysis", "setup_script.R"))
library(patchwork)

## ---- Species to plot ----
species_code <- "MAAM"  # American marten

## 1. Load PA patches and this species' candidate pairs ##########################
pa_patches <- st_read(here(pa_patches_dir, "pa_patches.gpkg"), quiet = TRUE)

pairs <- read_csv(
  here(pa_pairs_dir, paste0(species_code, "_PApairs_dkernel_cutoff.csv")),
  show_col_types = FALSE
)

reachable_ids <- union(pairs$from, pairs$to)

pa_patches <- pa_patches %>%
  mutate(reachable = site_id %in% reachable_ids)

message(species_code, ": ", sum(pa_patches$reachable), " of ", nrow(pa_patches),
        " PAs are potentially reachable (", nrow(pairs), " candidate pair(s)).")

## 2. Build the two panels ##########################
p_all <- ggplot(pa_patches) +
  geom_sf(fill = "grey40", color = NA) +
  labs(title = "All protected areas") +
  theme_void()

p_reachable <- ggplot(pa_patches) +
  geom_sf(fill = "grey85", color = NA) +
  geom_sf(data = filter(pa_patches, reachable), fill = "forestgreen", color = NA) +
  labs(title = paste0("Potentially reachable PAs (", species_code, ")")) +
  theme_void()

fig <- p_all + p_reachable +
  plot_annotation(
    title = paste0(species_code, ": dispersal-based PA reachability"),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

## 3. Save ##########################
out_path <- here(fig_dir, paste0(species_code, "_reachable_PAs.png"))
ggsave(out_path, fig, width = 10, height = 6, dpi = 300, bg = "white")
message("Saved figure to ", out_path)

## 4. Pairwise connection density (dispersal-kernel-limited) ##########################
## One edge-to-edge line per candidate pair -- the geometrically nearest
## points between the two PA polygons' boundaries (st_nearest_points(),
## pairwise so each pair gets its own line in one vectorized call), not
## centroids. Rasterized as a connection count per cell, then Gaussian-
## smoothed into a continuous density surface. Cropped to a buffered
## bounding box around this species' PAs, not the full modelling extent
## (never operate on the full extent -- see CLAUDE.md), since this is a
## lightweight figure-only raster.
from_polys <- pa_patches[match(pairs$from, pa_patches$site_id), ]
to_polys   <- pa_patches[match(pairs$to,   pa_patches$site_id), ]

pair_lines <- st_nearest_points(from_polys, to_polys, pairwise = TRUE) %>%
  st_as_sf() %>%
  mutate(from = pairs$from, to = pairs$to,
         euclidean_dist = pairs$euclidean_dist,
         kernel_prob = pairs$euclidean_kernel_prob)

## Persist for GIS use (e.g. QGIS/ArcGIS), alongside pa_patches.gpkg.
lines_out_path <- here(output_dir, paste0(species_code, "_PApairs_edge_lines.gpkg"))
st_write(pair_lines, lines_out_path, delete_dsn = TRUE, quiet = TRUE)
message("Saved edge-to-edge candidate-pair lines to ", lines_out_path)

## Coarse (1 km) raster, buffered 20 km around the PAs this species could
## ever reach -- plenty fine for a smoothed density figure, far cheaper than
## the pipeline's 30 m template_raster.
buffer_m <- 20000
pa_bbox <- st_bbox(pa_patches)
analysis_ext <- ext(pa_bbox["xmin"] - buffer_m, pa_bbox["xmax"] + buffer_m,
                     pa_bbox["ymin"] - buffer_m, pa_bbox["ymax"] + buffer_m)
density_raster <- rast(extent = analysis_ext, resolution = c(1000, 1000), crs = target_crs)

conn_count <- rasterizeGeom(vect(pair_lines), density_raster, fun = "count", unit = "m")

## Gaussian smoothing (5 km bandwidth) turns the raw per-cell line count into
## a continuous "connection density" surface.
smooth_w <- focalMat(conn_count, d = 5000, type = "Gauss")
conn_density <- focal(conn_count, w = smooth_w, fun = "sum", na.rm = TRUE)
names(conn_density) <- "density"

p_density <- ggplot() +
  geom_spatraster(data = conn_density) +
  scale_fill_viridis_c(name = "Connection\ndensity", na.value = "transparent") +
  geom_sf(data = pa_patches, fill = NA, color = "white", linewidth = 0.1) +
  labs(title = paste0(species_code, ": pairwise connection density (dispersal-kernel-limited)")) +
  theme_void()

density_out_path <- here(fig_dir, paste0(species_code, "_connection_density.png"))
ggsave(density_out_path, p_density, width = 8, height = 6, dpi = 300, bg = "white")
message("Saved figure to ", density_out_path)

## 5. Individual connection lines (network view) ##########################
## Same candidate-pair lines as above, drawn directly rather than smoothed --
## low alpha so overlapping lines still visually build up into denser-looking
## bundles where many candidate connections run close together.
p_network <- ggplot() +
  geom_sf(data = pa_patches, fill = "grey85", color = NA) +
  geom_sf(data = pair_lines, aes(color = kernel_prob), alpha = 0.35, linewidth = 0.3) +
  scale_color_viridis_c(name = "Kernel\nprobability") +
  labs(title = paste0(species_code, ": candidate PA-pair connections (edge-to-edge, dispersal-kernel-limited)")) +
  theme_void()

network_out_path <- here(fig_dir, paste0(species_code, "_connection_network.png"))
ggsave(network_out_path, p_network, width = 8, height = 6, dpi = 300, bg = "white")
message("Saved figure to ", network_out_path)
