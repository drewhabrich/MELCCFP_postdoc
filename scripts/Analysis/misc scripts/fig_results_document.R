## ------------------------------- ##
## Script name: fig_results_document
##
## Purpose of script: Figures for "MST results summary.docx". American marten
## (MAAM) is the worked example for every methods step; the results figures
## add MAAM's tree and overlap, a five-species comparison of the trees, and
## MAAM's sensitivity maps. Labels avoid acronyms to match the document text.
##   report_fig1_network_MAAM.png    -- A. dispersal kernel, B. candidate links,
##                                      C. resistance, D. minimum spanning tree
##   report_fig2_corridors_MAAM.png  -- 100 km window: A. straight tree links vs
##                                      least-cost corridors, B. corridor density,
##                                      C. distribution model, D. overlap
##   report_fig3_results_MAAM.png    -- A. tree links by detour factor, B. overlap
##   report_fig4_trees_5species.png  -- each analysed species' tree by detour factor
##   report_fig5_sensitivity_MAAM.png -- A. share of scenarios keeping each link,
##                                      B. coefficient of variation of the overlap
## Reads only finished stage outputs (01, 02, 04-07, 09); nothing is recomputed.
## The 30 m resistance raster is only sampled for display (maxcell), never
## processed over the full extent.
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))
library(cowplot)

## ---- Species ----
species_code <- "MAAM"  # American marten -- the worked example
species_all  <- c(ASFL = "Short-eared owl", AQCH = "Golden eagle", DOOR = "Bobolink",
                  DRPI = "Pileated woodpecker", MAAM = "American marten")  # by decreasing dispersal
window_km <- 100        # side of the zoomed square in figure 2 (densest part of the tree)
dpi <- 300

specieslist  <- readxl::read_excel(here(tabular_dir, "specieslist.xlsx"))
dispersal_km <- specieslist$dispersal_km[specieslist$species == species_code]
species_meta <- read_csv(here(tabular_dir, "species_metadata.csv"), show_col_types = FALSE)
sdm_path     <- species_meta$sdm_raster_path[species_meta$species == species_code]
stopifnot(length(dispersal_km) == 1, file.exists(sdm_path))

## ---- Shared layers and styling ----
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

## detour factor = cost distance / straight-line length (touching patches have
## length 0 and no defined detour)
with_detour <- function(lines) mutate(lines, detour = if_else(euclidean_dist > 0, weight / euclidean_dist, NA_real_))
detour_scale <- function() scale_colour_viridis_c(
  option = "plasma", trans = "log10", limits = c(1, 80), oob = scales::squish,
  breaks = c(1, 2, 5, 10, 20, 50), na.value = "grey40",
  name = "Detour factor (cost / straight-line length)", guide = bar_guide)

save_fig <- function(p, name, width, height) {
  out <- here(fig_dir, name)
  ggsave(out, p, width = width, height = height, dpi = dpi, bg = "white")
  message("Saved ", out)
}

## ---- Stage outputs for the example species ----
cand_lines <- st_read(here(pa_pairs_dir, paste0(species_code, "_PApairs_edge_lines.gpkg")), quiet = TRUE)
mst_lines  <- with_detour(st_read(here(mst_lines_dir, paste0(species_code, "_mst_lines.gpkg")), quiet = TRUE))
lcp_lines  <- st_read(here(corridor_dir, paste0(species_code, "_lcp_corridors.gpkg")), quiet = TRUE)
res_r      <- rast(here(resistance_dir, paste0(species_code, "_res.tif")))
density_r  <- rast(here(density_dir, paste0(species_code, "_corridor_density.tif")))
sdm_r      <- rast(sdm_path)
overlap_r  <- rast(here(sdm_overlap_dir, paste0(species_code, "_habitat_overlap.tif")))

## tree components (the tree is a forest when some patches are out of range)
g <- igraph::graph_from_data_frame(st_drop_geometry(mst_lines)[, c("from", "to")], directed = FALSE)
memb <- igraph::components(g)$membership
mst_lines$component <- factor(memb[as.character(mst_lines$from)])
n_comp <- nlevels(mst_lines$component)

## densest window_km square lying inside the study area, by corridor midpoints
## on a 10 km grid
mid  <- st_coordinates(st_centroid(st_geometry(lcp_lines)))
bb   <- st_bbox(lcp_lines); w <- window_km * 1000
cand <- expand.grid(x0 = seq(bb["xmin"], bb["xmax"] - w, by = 10000),
                    y0 = seq(bb["ymin"], bb["ymax"] - w, by = 10000))
corners <- st_as_sf(data.frame(x = c(cand$x0, cand$x0 + w, cand$x0, cand$x0 + w),
                               y = c(cand$y0, cand$y0, cand$y0 + w, cand$y0 + w)),
                    coords = c("x", "y"), crs = st_crs(pa_patches))
inside <- matrix(lengths(st_intersects(corners, study_area)) > 0, ncol = 4)
cand <- cand[rowSums(inside) == 4, ]
cand$n <- mapply(function(x0, y0) sum(mid[, 1] >= x0 & mid[, 1] < x0 + w & mid[, 2] >= y0 & mid[, 2] < y0 + w),
                 cand$x0, cand$y0)
best <- cand[which.max(cand$n), ]
win  <- ext(best$x0, best$x0 + w, best$y0, best$y0 + w)
win_sf <- st_as_sfc(st_bbox(c(xmin = best$x0, xmax = best$x0 + w, ymin = best$y0, ymax = best$y0 + w),
                            crs = st_crs(pa_patches)))
message("Zoom window: ", best$n, " of ", nrow(lcp_lines), " corridors centred in the densest ", window_km, " km square.")

## 1. Figure 1: building the network ###########################################
alpha_k   <- dispersal_kernel_alpha(p_ref = 0.5, d_ref = dispersal_km * 1000)
cutoff_km <- -log(kernel_p_threshold) / alpha_k / 1000
kern  <- data.frame(d = seq(0, ceiling(cutoff_km * 1.2), length.out = 400))
kern$p <- exp(-alpha_k * kern$d * 1000)
marks <- data.frame(d = c(dispersal_km, cutoff_km), p = c(0.5, kernel_p_threshold),
                    lab = c(sprintf("median dispersal\n%.0f km (p = 0.50)", dispersal_km),
                            sprintf("search radius\n%.0f km (p = %.2f)", cutoff_km, kernel_p_threshold)),
                    hj = c(-0.08, 0.5), ly = c(0.5, 0.13))

p1a <- ggplot(kern, aes(d, p)) +
  geom_area(data = subset(kern, d <= cutoff_km), fill = "#2a9d8f", alpha = 0.15) +
  geom_line(linewidth = 0.9, colour = "grey15") +
  geom_segment(data = marks, aes(x = d, xend = d, y = 0, yend = p), linetype = "dashed", colour = "#2a9d8f") +
  geom_point(data = marks, size = 2.4, colour = "#2a9d8f") +
  geom_text(data = marks, aes(y = ly, label = lab, hjust = hj), vjust = -0.25, size = 3.1, lineheight = 0.9) +
  scale_y_continuous(limits = c(0, 1.02), expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A. Dispersal kernel", x = "Edge-to-edge distance between protected areas (km)",
       y = "Dispersal probability") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11),
        plot.margin = margin(6, 12, 6, 6))

p1b <- ggplot() + base_layers() +
  geom_sf(data = cand_lines, aes(colour = euclidean_kernel_prob), linewidth = 0.2, alpha = 0.7) +
  scale_colour_viridis_c(name = "Dispersal probability", limits = c(kernel_p_threshold, 1), guide = bar_guide) +
  full_extent() + labs(title = sprintf("B. Candidate links (%s pairs within range)", format(nrow(cand_lines), big.mark = ","))) +
  map_theme

res_breaks <- 2^(0:6)
p1c <- ggplot() +
  geom_spatraster(data = res_r, maxcell = 2e6) +
  geom_sf(data = study_area, fill = NA, colour = "grey30", linewidth = 0.25) +
  scale_fill_viridis_c(option = "rocket", direction = -1, trans = "log2", breaks = res_breaks,
                       labels = res_breaks, na.value = "transparent", name = "Resistance (log scale)",
                       guide = bar_guide) +
  full_extent() + labs(title = "C. Resistance to movement (30 m)") + map_theme

## components numbered largest first; the smaller ones are drawn on top, thicker
comp_size <- table(mst_lines$component)
mst_lines$component <- factor(mst_lines$component, levels = names(sort(comp_size, decreasing = TRUE)),
                              labels = seq_len(n_comp))
## protected areas in the smaller component(s) are filled in that component's
## colour too -- their links are short and easy to miss at this scale
comp_of_site <- setNames(as.character(mst_lines$component), mst_lines$from)
comp_of_site <- c(comp_of_site, setNames(as.character(mst_lines$component), mst_lines$to))
pa_minor <- pa_patches %>%
  mutate(component = factor(comp_of_site[as.character(site_id)], levels = levels(mst_lines$component))) %>%
  filter(!is.na(component), component != "1")
p1d <- ggplot() + base_layers() +
  geom_sf(data = pa_minor, aes(fill = component), colour = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("#1d3557", "#e76f51", "#2a9d8f", "#8338ec")[seq_len(n_comp)], drop = FALSE) +
  geom_sf(data = arrange(mst_lines, component), aes(colour = component, linewidth = component != "1")) +
  scale_linewidth_manual(values = c(`FALSE` = 0.45, `TRUE` = 1.3), guide = "none") +
  geom_sf(data = win_sf, fill = NA, colour = "black", linewidth = 0.5) +
  scale_colour_manual(values = c("#1d3557", "#e76f51", "#2a9d8f", "#8338ec")[seq_len(n_comp)],
                      name = "Network component") +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5,
                               override.aes = list(linewidth = 1.5))) +
  full_extent() +
  labs(title = sprintf("D. Minimum spanning tree (%d links, %d components)", nrow(mst_lines), n_comp)) +
  map_theme

fig1 <- plot_grid(p1a, p1b, p1c, p1d, ncol = 2, align = "hv", axis = "tblr")
save_fig(fig1, paste0("report_fig1_network_", species_code, ".png"), width = 11, height = 7.6)

## 2. Figure 2: corridors and distribution (zoom window) #######################
zoom <- function() coord_sf(xlim = c(win[1], win[2]), ylim = c(win[3], win[4]), expand = FALSE, datum = NA)
pa_win   <- st_filter(pa_patches, win_sf)
lcp_win  <- st_filter(lcp_lines, win_sf)
mst_win  <- st_filter(mst_lines, win_sf)
## protected areas: 50% transparent grey fill so they read over every surface
outline  <- function(col = "grey20") geom_sf(data = pa_win, fill = "grey50", alpha = 0.5, colour = col, linewidth = 0.3)
crop_to  <- function(r) crop(r, win)

p2a <- ggplot() +
  geom_spatraster(data = crop_to(res_r), maxcell = 1.2e6) +
  ## greyscale resistance here so the coloured links stand out (figure 1C keeps the colour scale)
  scale_fill_gradient(low = "white", high = "grey10", trans = "log2", breaks = res_breaks,
                      labels = res_breaks, na.value = "transparent", name = "Resistance (log scale)",
                      guide = bar_guide) +
  outline() +
  ## white halo under each line keeps it visible over dark (high-resistance) cells
  geom_sf(data = lcp_win, colour = "white", linewidth = 1.8) +
  geom_sf(data = mst_win, colour = "white", linewidth = 1.4) +
  geom_sf(data = lcp_win, aes(colour = "Least-cost corridor"), linewidth = 0.9) +
  geom_sf(data = mst_win, aes(colour = "Straight tree link"), linetype = "22", linewidth = 0.7) +
  scale_colour_manual(values = c("Least-cost corridor" = "#d7191c", "Straight tree link" = "#0057e7"),
                      name = NULL, guide = guide_legend(override.aes = list(linewidth = 1.2,
                                                                            linetype = c("solid", "22")))) +
  zoom() + labs(title = "A. Tree links traced as least-cost corridors") + map_theme +
  theme(legend.box = "vertical", legend.spacing.y = grid::unit(2, "pt"))

p2b <- ggplot() +
  geom_spatraster(data = crop_to(density_r)) +
  scale_fill_viridis_c(option = "mako", direction = -1, na.value = "transparent",
                       name = expression("Corridor density (km per km"^2*")"), guide = bar_guide) +
  outline() + geom_sf(data = lcp_win, colour = "#e76f51", linewidth = 0.3) +
  zoom() + labs(title = sprintf("B. Smoothed corridor density (standard deviation %.1f km)",
                                corridor_sigma(dispersal_km) / 1000)) + map_theme

p2c <- ggplot() +
  geom_spatraster(data = crop_to(sdm_r)) +
  scale_fill_viridis_c(option = "viridis", na.value = "transparent", limits = c(0, 1),
                       name = "Predicted probability of occurrence", guide = bar_guide) +
  outline() +
  zoom() + labs(title = "C. Species distribution model (Eckert et al. 2023)") + map_theme

p2d <- ggplot() +
  geom_spatraster(data = crop_to(overlap_r)) +
  scale_fill_viridis_c(option = "inferno", na.value = "transparent", limits = c(0, 1),
                       name = "Corridor density × occurrence (0–1)", guide = bar_guide) +
  outline() +
  zoom() + labs(title = "D. Overlap of corridors and distribution") + map_theme

fig2 <- plot_grid(p2a, p2b, p2c, p2d, ncol = 2, align = "hv", axis = "tblr")
save_fig(fig2, paste0("report_fig2_corridors_", species_code, ".png"), width = 10, height = 11)

## 3. Figure 3: example species results ########################################
p3a <- ggplot() + base_layers() +
  geom_sf(data = arrange(mst_lines, detour), aes(colour = detour), linewidth = 0.55) +
  detour_scale() + full_extent() +
  labs(title = "A. Minimum spanning tree links by detour factor") + map_theme

p3b <- ggplot() +
  geom_sf(data = study_area, fill = "grey96", colour = "grey45", linewidth = 0.25) +
  geom_spatraster(data = overlap_r) +
  geom_sf(data = pa_patches, fill = NA, colour = "grey35", linewidth = 0.15) +
  scale_fill_viridis_c(option = "inferno", na.value = "transparent", limits = c(0, 1),
                       name = "Corridor density × occurrence (0–1)", guide = bar_guide) +
  full_extent() + labs(title = "B. Overlap of corridors and distribution") + map_theme

fig3 <- plot_grid(p3a, p3b, ncol = 1, align = "v")
save_fig(fig3, paste0("report_fig3_results_", species_code, ".png"), width = 9, height = 10)

## 4. Figure 4: the five analysed species' trees ###############################
tree_panel <- function(sp) {
  l <- with_detour(st_read(here(mst_lines_dir, paste0(sp, "_mst_lines.gpkg")), quiet = TRUE))
  ggplot() + base_layers() +
    geom_sf(data = arrange(l, detour), aes(colour = detour), linewidth = 0.4) +
    detour_scale() + full_extent() +
    labs(title = sprintf("%s (median dispersal %s km)", species_all[[sp]],
                         format(specieslist$dispersal_km[specieslist$species == sp]))) +
    map_theme + theme(legend.position = "none")
}
panels4 <- lapply(names(species_all), tree_panel)
legend4 <- get_legend(ggplot() + geom_sf(data = mst_lines, aes(colour = detour)) + detour_scale() +
                        map_theme + theme(legend.key.width = grid::unit(2.4, "cm")))
fig4 <- plot_grid(plotlist = c(panels4, list(legend4)), ncol = 2)
save_fig(fig4, "report_fig4_trees_5species.png", width = 11, height = 9)

## 5. Figure 5: example species sensitivity ####################################
retention <- st_read(here(sensitivity_out_dir, paste0(species_code, "_mst_edge_retention.gpkg")), quiet = TRUE)
var_r <- rast(here(sensitivity_out_dir, paste0(species_code, "_overlap_variability.tif")))
## CV is unstable where the mean overlap is near zero -- show it only where the
## mean is at least 1% of its maximum
mean_max <- global(var_r[["mean"]], "max", na.rm = TRUE)[1, 1]
cv_r <- ifel(var_r[["mean"]] >= 0.01 * mean_max, var_r[["cv"]], NA)

p5a <- ggplot() + base_layers() +
  geom_sf(data = arrange(retention, desc(retention_all)), aes(colour = retention_all), linewidth = 0.55) +
  scale_colour_viridis_c(option = "viridis", limits = c(0, 1), labels = scales::percent,
                         name = sprintf("Share of the %d scenarios keeping the link", retention$n_scenarios[1]),
                         guide = bar_guide) +
  full_extent() + labs(title = "A. Robustness of baseline tree links") + map_theme

p5b <- ggplot() +
  geom_sf(data = study_area, fill = "grey96", colour = "grey45", linewidth = 0.25) +
  geom_spatraster(data = cv_r) +
  geom_sf(data = pa_patches, fill = NA, colour = "grey35", linewidth = 0.15) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 2), oob = scales::squish,
                       na.value = "transparent", name = "Coefficient of variation of the overlap (capped at 2)",
                       guide = bar_guide) +
  full_extent() + labs(title = "B. Variability of the overlap across scenarios") + map_theme

fig5 <- plot_grid(p5a, p5b, ncol = 1, align = "v")
save_fig(fig5, paste0("report_fig5_sensitivity_", species_code, ".png"), width = 9, height = 10)
