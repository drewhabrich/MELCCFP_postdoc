## ------------------------------- ##
## Script name: fig_deliverable3
##
## Purpose of script: One four-panel figure per example species for
## "Deliverable 3 - MST results summary":
##   A. resistance surface, B. candidate links (01), C. minimum spanning
##   tree by detour factor (04), D. corridor density x SDM overlap (07).
## Labels avoid acronyms and the styling follows fig_results_document.R, to
## match the document text. Reads only finished stage outputs; the 30 m
## resistance raster is only sampled for display (maxcell), never processed
## over the full extent.
## Output: Output/Figures/deliverable3_maps_{species}.png
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))
library(cowplot)

species_all <- c(ASFL = "Short-eared owl", AQCH = "Golden eagle", DOOR = "Bobolink",
                 DRPI = "Pileated woodpecker", MAAM = "American marten")  # by decreasing dispersal
dpi <- 200

pa_patches <- st_read(here(pa_patches_dir, "pa_patches.gpkg"), quiet = TRUE)
study_area <- st_transform(st_read(here(output_dir, "study_area.shp"), quiet = TRUE), st_crs(pa_patches))
map_bbox   <- st_bbox(st_union(st_as_sfc(st_bbox(study_area)), st_as_sfc(st_bbox(pa_patches))))

base_layers <- function() list(
  geom_sf(data = study_area, fill = "grey96", colour = "grey45", linewidth = 0.25),
  geom_sf(data = pa_patches, fill = "grey72", colour = NA)
)
full_extent <- function() coord_sf(xlim = map_bbox[c("xmin", "xmax")], ylim = map_bbox[c("ymin", "ymax")],
                                   expand = FALSE, datum = NA)
map_theme <- theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom",
        legend.title = element_text(size = 9), legend.text = element_text(size = 8),
        legend.key.width = grid::unit(1.4, "cm"), legend.key.height = grid::unit(0.25, "cm"),
        plot.margin = margin(4, 6, 4, 6))
bar_guide <- guide_colourbar(title.position = "top", title.hjust = 0.5)
res_breaks <- 2^(0:6)

for (i in seq_along(species_all)) {
  sp <- names(species_all)[i]
  message("Figure ", i, ": ", sp)

  res_r      <- rast(here(resistance_dir, paste0(sp, "_res.tif")))
  cand_lines <- st_read(here(pa_pairs_dir, paste0(sp, "_PApairs_edge_lines.gpkg")), quiet = TRUE)
  mst_lines  <- st_read(here(mst_lines_dir, paste0(sp, "_mst_lines.gpkg")), quiet = TRUE) %>%
    mutate(detour = if_else(euclidean_dist > 0, weight / euclidean_dist, NA_real_))
  overlap_r  <- rast(here(sdm_overlap_dir, paste0(sp, "_habitat_overlap.tif")))

  ## PAs outside the largest MST component, highlighted in panel C
  g    <- igraph::graph_from_data_frame(st_drop_geometry(mst_lines)[, c("from", "to")], directed = FALSE)
  memb <- igraph::components(g)$membership
  main <- as.integer(names(which.max(table(memb))))
  minor_ids <- names(memb)[memb != main]
  pa_minor  <- filter(pa_patches, as.character(site_id) %in% minor_ids)
  n_comp    <- length(unique(memb))

  pa <- ggplot() +
    geom_spatraster(data = res_r, maxcell = 2e6) +
    geom_sf(data = study_area, fill = NA, colour = "grey30", linewidth = 0.25) +
    scale_fill_viridis_c(option = "rocket", direction = -1, trans = "log2", breaks = res_breaks,
                         labels = res_breaks, na.value = "transparent", name = "Resistance (log scale)",
                         guide = bar_guide) +
    full_extent() + labs(title = "A. Resistance to movement (30 m)") + map_theme

  pb <- ggplot() + base_layers() +
    geom_sf(data = arrange(cand_lines, euclidean_kernel_prob), aes(colour = euclidean_kernel_prob),
            linewidth = 0.15, alpha = 0.6) +
    scale_colour_viridis_c(name = "Dispersal probability", limits = c(kernel_p_threshold, 1), guide = bar_guide) +
    full_extent() +
    labs(title = sprintf("B. Candidate links (%s pairs within range)", format(nrow(cand_lines), big.mark = ","))) +
    map_theme

  pc <- ggplot() + base_layers() +
    ## second component in orange, as in report figure 1D
    geom_sf(data = pa_minor, fill = "#e76f51", colour = "#e76f51", linewidth = 0.6) +
    geom_sf(data = arrange(mst_lines, detour), aes(colour = detour), linewidth = 0.5) +
    scale_colour_viridis_c(option = "plasma", trans = "log10", limits = c(1, 80), oob = scales::squish,
                           breaks = c(1, 2, 5, 10, 20, 50), na.value = "grey40",
                           name = "Detour factor (cost / straight-line length)", guide = bar_guide) +
    full_extent() +
    labs(title = sprintf("C. Minimum spanning tree (%d links, %d component%s)",
                         nrow(mst_lines), n_comp, if (n_comp > 1) "s" else "")) +
    map_theme

  pd <- ggplot() +
    geom_sf(data = study_area, fill = "grey96", colour = "grey45", linewidth = 0.25) +
    geom_spatraster(data = overlap_r) +
    geom_sf(data = pa_patches, fill = NA, colour = "grey35", linewidth = 0.15) +
    scale_fill_viridis_c(option = "inferno", na.value = "transparent", limits = c(0, 1),
                         name = "Corridor density × occurrence (0–1)", guide = bar_guide) +
    full_extent() + labs(title = "D. Overlap of corridors and distribution model") + map_theme

  fig <- plot_grid(pa, pb, pc, pd, ncol = 2, align = "hv", axis = "tblr")
  out <- here(fig_dir, sprintf("deliverable3_maps_%s.png", sp))
  ggsave(out, fig, width = 11, height = 8, dpi = dpi, bg = "white")
  message("Saved ", out)
}
