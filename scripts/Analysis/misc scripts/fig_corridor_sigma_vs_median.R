## ------------------------------- ##
## Script name: fig_corridor_sigma_vs_median
##
## Purpose of script: Explanatory figure for script 06's corridor-smoothing
## bandwidth, for one species: why sigma = corridor_sigma(dispersal_km)
## (= median dispersal / ln(2), the SD = mean of the negative-exponential
## dispersal kernel) rather than the median dispersal distance itself.
## Four panels:
##   A. the dispersal kernel p(d), with the median, sigma, and the p = 0.05
##      search cutoff (scripts 01/04) marked
##   B. the Gaussian cross-section a single traced corridor line is smoothed
##      with, for sigma = median vs sigma = corridor_sigma() (equal total
##      mass, so the wider one is lower)
##   C./D. the species' real traced corridors (script 05) smoothed both ways
##      in the densest `window_km` window, on one shared colour scale
##
## C/D use the production worker (edge_density_worker() in
## pipeline_functions.R) on the same species grid rule as script 06, so the
## sigma = corridor_sigma() panel reproduces the density script 06 wrote
## (checked against {species}_corridor_density.tif below); only sigma differs
## between the two panels.
##
## Saves Output/Figures/{species}_corridor_sigma_vs_median.png
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
source(file.path("Scripts", "Analysis", "setup_script.R"))

## ---- Species to plot ----
species_code <- "MAAM"  # American marten
window_km <- 200        # side of the square map window (densest part of the network)

specieslist <- readxl::read_excel(here(tabular_dir, "specieslist.xlsx"))
dispersal_km <- specieslist$dispersal_km[specieslist$species == species_code]
stopifnot(length(dispersal_km) == 1)

## 1. Sigmas and kernel #########################################################
## (If a species' dispersal_km exceeds corridor_sigma()'s 100 km cap, the
## "median" sigma below is still the uncapped median -- the comparison is
## then median vs the capped production value.)
sigma_med  <- dispersal_km * 1000            # naive choice: sigma = median dispersal distance
sigma_prod <- corridor_sigma(dispersal_km)   # production: SD of the kernel = median / ln(2)
alpha      <- dispersal_kernel_alpha(p_ref = 0.5, d_ref = dispersal_km * 1000)
cutoff_km  <- -log(kernel_p_threshold) / alpha / 1000

sig <- "σ"; alp <- "α"; pm <- "±"; sq <- "²"
col_med  <- "#4C78A8"   # sigma = median dispersal
col_prod <- "#E4572E"   # sigma = corridor_sigma()

## 2. Panel A: the dispersal kernel #############################################
kern <- data.frame(d = seq(0, 62, length.out = 400))
kern$p <- exp(-alpha * kern$d * 1000)

marks <- data.frame(
  d = c(dispersal_km, sigma_prod / 1000, cutoff_km),
  col = c(col_med, col_prod, "grey35"),
  lab = c(sprintf("median dispersal\n%.1f km (p = 0.50)", dispersal_km),
          sprintf("%s = 1/%s = mean = SD\n%.1f km (p = %.2f)", sig, alp, sigma_prod / 1000, exp(-1)),
          sprintf("search cutoff\n%.1f km (p = %.2f)", cutoff_km, kernel_p_threshold))
)
marks$p <- exp(-alpha * marks$d * 1000)

pA <- ggplot(kern, aes(d, p)) +
  geom_line(linewidth = 0.9, colour = "grey20") +
  geom_segment(data = marks, aes(x = d, xend = d, y = 0, yend = p, colour = I(col)),
               linetype = "dashed", linewidth = 0.6) +
  geom_point(data = marks, aes(colour = I(col)), size = 2.6) +
  annotate("text", x = marks$d[1] - 1, y = 0.02, label = marks$lab[1], colour = col_med,
           hjust = 1, vjust = 0, size = 3.3, lineheight = 0.95) +
  annotate("text", x = marks$d[2] + 1, y = 0.02, label = marks$lab[2], colour = col_prod,
           hjust = 0, vjust = 0, size = 3.3, lineheight = 0.95) +
  annotate("text", x = marks$d[3], y = marks$p[3] + 0.08, label = marks$lab[3], colour = "grey35",
           hjust = 0.5, vjust = 0, size = 3.3, lineheight = 0.95) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 1.02), expand = expansion(mult = c(0, 0))) +
  labs(title = "A. Dispersal kernel (negative exponential)",
       subtitle = sprintf("p(d) = exp(-%s d), %s = ln2 / %.0f km", alp, alp, dispersal_km),
       x = "Dispersal distance d (km)", y = "Dispersal probability p(d)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

## 3. Panel B: cross-section of the smoothed corridor ##########################
lat <- seq(-3.2 * sigma_prod / 1000, 3.2 * sigma_prod / 1000, length.out = 600)
prof <- rbind(
  data.frame(x = lat, y = dnorm(lat, 0, sigma_med / 1000),  who = "med"),
  data.frame(x = lat, y = dnorm(lat, 0, sigma_prod / 1000), who = "prod")
)
band <- function(curve, s_km) {
  b <- prof[prof$who == curve & abs(prof$x) <= s_km, ]
  b$ymin <- 0
  b
}

pB <- ggplot() +
  geom_ribbon(data = band("med", sigma_med / 1000),   aes(x, ymin = ymin, ymax = y), fill = col_med,  alpha = 0.18) +
  geom_ribbon(data = band("prod", sigma_prod / 1000), aes(x, ymin = ymin, ymax = y), fill = col_prod, alpha = 0.18) +
  geom_line(data = prof, aes(x, y, colour = who), linewidth = 0.95) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
  annotate("text", x = -2, y = max(prof$y) * 1.03, label = "traced least-cost path (1 pixel wide)",
           size = 3.2, colour = "grey30", hjust = 1, vjust = 0) +
  scale_colour_manual(
    name = NULL,
    values = c(med = col_med, prod = col_prod),
    labels = c(med  = sprintf("%s = median dispersal = %.1f km\n68%% of density within %s%.1f km", sig, sigma_med / 1000, pm, sigma_med / 1000),
               prod = sprintf("%s = kernel SD = %.1f km (production)\n68%% of density within %s%.1f km", sig, sigma_prod / 1000, pm, sigma_prod / 1000))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "B. Smoothing applied across one corridor",
       subtitle = "Same corridor length, spread over a wider swath (shaded = ±1 sigma)",
       x = "Distance from the least-cost path (km)", y = "Share of corridor density (per km)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.spacing.y = grid::unit(8, "pt"),
        legend.text = element_text(size = 9, lineheight = 1.0))

## 4. Panels C/D: real MAAM corridors, smoothed both ways ######################
edges <- st_read(here(corridor_dir, paste0(species_code, "_lcp_corridors.gpkg")), quiet = TRUE)

## Densest window_km x window_km square, by edge-midpoint count on a 25 km grid.
cent <- st_coordinates(st_centroid(st_geometry(edges)))
bb   <- st_bbox(edges)
w    <- window_km * 1000
cand <- expand.grid(x0 = seq(bb["xmin"], bb["xmax"] - w, by = 25000),
                    y0 = seq(bb["ymin"], bb["ymax"] - w, by = 25000))
cand$n <- mapply(function(x0, y0) sum(cent[, 1] >= x0 & cent[, 1] < x0 + w & cent[, 2] >= y0 & cent[, 2] < y0 + w),
                 cand$x0, cand$y0)
best <- cand[which.max(cand$n), ]
win  <- ext(best$x0, best$x0 + w, best$y0, best$y0 + w)
message("Window: ", best$n, " of ", nrow(edges), " edge midpoints in the densest ", window_km, " km square.")

box_sfc <- function(e) {
  e <- unname(e)  # as.vector(<SpatExtent>) is named, which would mangle c(xmin = e[1], ...)
  st_as_sfc(st_bbox(c(xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4]), crs = st_crs(edges)))
}
## Edges within 3 x production sigma of the window all contribute to it.
edges_calc <- st_filter(edges, box_sfc(as.vector(win + 3 * sigma_prod)))
edges_win  <- st_filter(edges, box_sfc(as.vector(win)))

## Smooth `edges_calc` with the production worker at a given sigma (same grid
## rule as corridor_density() in pipeline_functions.R, run serially here).
smooth_edges <- function(edges_sf, sigma) {
  v <- vect(edges_sf)
  edge_coords_list <<- lapply(seq_len(nrow(v)), function(e) terra::geom(v[e, ])[, c("x", "y")])
  template_crs_wkt <<- crs(template_raster)
  res_val          <- min(common_res_m, max(30, round(sigma / cells_per_sigma / 30) * 30))
  sp_grid_res      <<- c(res_val, res_val)
  sp_grid_ext_vec  <<- as.vector(ext(template_raster))
  crop_buffer_m    <<- 3 * sigma
  cell_area_km2    <<- res_val^2 / 1e6
  tmpl <- rast(ext(sp_grid_ext_vec), resolution = sp_grid_res, crs = template_crs_wkt)
  g2 <- focalMat(tmpl, sigma, "Gauss")
  g1 <- g2[(nrow(g2) + 1) / 2, ]
  gauss_1d <<- g1 / sum(g1)

  paths <- vapply(seq_len(nrow(v)), edge_density_worker, character(1))
  r_path <- tempfile(fileext = ".tif")
  writeRaster(mosaic(sprc(lapply(paths, rast)), fun = "sum"), r_path, overwrite = TRUE)
  unlink(paths)
  list(r = rast(r_path), cell_km2 = cell_area_km2, res = res_val)
}

dens_med  <- smooth_edges(edges_calc, sigma_med)
dens_prod <- smooth_edges(edges_calc, sigma_prod)

## Check 1: the sigma = corridor_sigma() panel reproduces script 06's output.
prod_file <- file.path(density_dir, paste0(species_code, "_corridor_density.tif"))
if (file.exists(prod_file)) {
  a <- crop(dens_prod$r, win)
  b <- crop(rast(prod_file), win)
  if (isTRUE(all.equal(res(a), res(b), tolerance = 1e-6))) {
    d <- abs(a - resample(b, a));
    message(sprintf("Check vs %s: max |diff| = %.2e km/km2 (%.3f%% of window peak).",
                    basename(prod_file), global(d, "max", na.rm = TRUE)[1, 1],
                    100 * global(d, "max", na.rm = TRUE)[1, 1] / global(a, "max", na.rm = TRUE)[1, 1]))
  } else message("Check vs production tif skipped: grid differs (", paste(round(res(a)), collapse = "x"),
                 " vs ", paste(round(res(b)), collapse = "x"), " m) -- rerun 06?")
}
## Check 2: smoothing conserves corridor length (same mass either sigma).
len_km <- sum(as.numeric(st_length(edges_calc))) / 1000
message(sprintf("Corridor length in calc set: %.0f km | smoothed mass: sigma=median %.0f km, sigma=prod %.0f km",
                len_km, global(dens_med$r, "sum", na.rm = TRUE)[1, 1] * dens_med$cell_km2,
                global(dens_prod$r, "sum", na.rm = TRUE)[1, 1] * dens_prod$cell_km2))

## Cells no edge's crop reached (beyond 3 sigma of every corridor) are NA in
## the raster; for display they are zero density, so fill them -- otherwise the
## smaller-sigma panel shows white gaps where its crop buffer ends.
fill0  <- function(r) subst(extend(crop(r, win), win), NA, 0)
r_med  <- fill0(dens_med$r)
r_prod <- fill0(dens_prod$r)
dens_max <- max(global(r_med, "max", na.rm = TRUE)[1, 1], global(r_prod, "max", na.rm = TRUE)[1, 1])
peak_med  <- global(r_med,  "max", na.rm = TRUE)[1, 1]
peak_prod <- global(r_prod, "max", na.rm = TRUE)[1, 1]

pa_patches <- st_read(here(pa_patches_dir, "pa_patches.gpkg"), quiet = TRUE)
pa_win <- st_filter(pa_patches, box_sfc(as.vector(win)))

map_panel <- function(r, sigma_m, col, title, peak) {
  ggplot() +
    geom_spatraster(data = r, maxcell = 5e5) +
    geom_sf(data = pa_win, fill = "grey55", colour = NA, alpha = 0.55) +
    geom_sf(data = edges_win, colour = "#E4572E", linewidth = 0.28) +
    ## sigma scale bar on a white card (readable over any fill)
    annotate("rect", xmin = win[1] + 0.02 * w, xmax = win[1] + 0.27 * w,
             ymin = win[3] + 0.02 * w, ymax = win[3] + 0.105 * w, fill = "white", alpha = 0.88) +
    annotate("text", x = win[1] + 0.035 * w, y = win[3] + 0.07 * w,
             label = sprintf("%s = %.1f km", sig, sigma_m / 1000), size = 3.5, fontface = "bold",
             colour = col, hjust = 0, vjust = 0) +
    annotate("segment", x = win[1] + 0.035 * w, xend = win[1] + 0.035 * w + sigma_m, y = win[3] + 0.045 * w,
             yend = win[3] + 0.045 * w, colour = col, linewidth = 1.8) +
    scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, dens_max), na.value = NA,
                         name = expression("Corridor density (km / km"^2*")"),
                         guide = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    coord_sf(xlim = c(win[1], win[2]), ylim = c(win[3], win[4]), expand = FALSE, datum = NA) +
    labs(title = title, subtitle = sprintf("peak density %.3f km/km%s", peak, sq)) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
          legend.key.width = grid::unit(1.6, "cm"), legend.key.height = grid::unit(0.3, "cm"),
          plot.margin = margin(4, 8, 4, 8))
}

pC <- map_panel(r_med,  sigma_med,  col_med,
                sprintf("C. %s = median dispersal distance", sig), peak_med)
pD <- map_panel(r_prod, sigma_prod, col_prod,
                sprintf("D. %s = kernel SD (production)", sig), peak_prod)

## 5. Assemble + save ###########################################################
out_path <- here(fig_dir, paste0(species_code, "_corridor_sigma_vs_median.png"))
png(out_path, width = 14, height = 10.5, units = "in", res = 200, type = "cairo")
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(
  3, 2, heights = grid::unit(c(0.085, 0.38, 0.535), "npc"))))
vp <- function(r, c) grid::viewport(layout.pos.row = r, layout.pos.col = c)
grid::grid.text(
  sprintf("%s: corridor smoothing %s vs. median dispersal distance", species_code, sig),
  x = 0.5, y = 0.62, vp = vp(1, 1:2), gp = grid::gpar(fontsize = 17, fontface = "bold"))
grid::grid.text(
  sprintf("Kernel SD (%s = median / ln 2 = %.1f km) is the resistant-kernel bandwidth of Compton et al. (2007), applied to this project's exponential dispersal kernel (median dispersal = %.0f km)",
          sig, sigma_prod / 1000, dispersal_km),
  x = 0.5, y = 0.2, vp = vp(1, 1:2), gp = grid::gpar(fontsize = 10.5, col = "grey30"))
print(pA, vp = vp(2, 1)); print(pB, vp = vp(2, 2))
print(pC, vp = vp(3, 1)); print(pD, vp = vp(3, 2))
invisible(dev.off())
message("Saved figure to ", out_path)
