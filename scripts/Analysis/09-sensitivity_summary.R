## ------------------------------- ##
## Script name: 09-sensitivity_summary
##
## Purpose of script: Summarise the one-at-a-time sensitivity scenarios run
# by 08-sensitivity_runs.R (one _result.csv row per scenario, each already
# compared with the species' production outputs):
#   Output/Tables/09_sensitivity_results.csv
#       every scenario row, plus changes vs. production's 04_mst_summary.csv
#   Output/Figures/09_{sp}_sensitivity_tornado.png
#       per input and level: share of MST edges changed (1 - edge Jaccard),
#       corridor location changed (1 - share of corridor length within 1 km
#       of production's corridors; Beier et al. 2009), SDM overlap changed
#       (1 - Schoener's D; Warren et al. 2008), change in MST components
#   Output/09_sensitivity/{sp}_mst_edge_retention.gpkg
#       production's MST edges, each with the share of scenarios that keep it
#       (all scenarios / dispersal scenarios / contrast scenarios)
#   Output/09_sensitivity/{sp}_overlap_variability.tif
#       per-cell mean, SD and CV of the scenarios' SDM overlap (1 km)
#       ({sp}_density_variability.tif for species with no SDM)
# The "identity_check" test scenario, if run, stays in the table but is
# left out of the figure and maps.
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))

species_subset <- c("MAAM", "DOOR", "DRPI")  # NULL = every species with scenario results

result_files <- list.files(sensitivity_dir, pattern = "_result[.]csv$", recursive = TRUE, full.names = TRUE)
results <- bind_rows(lapply(result_files, read_csv, show_col_types = FALSE))
if (!is.null(species_subset)) results <- results %>% filter(species %in% species_subset)
if (nrow(results) == 0) stop("No scenario results in ", sensitivity_dir, " -- run 08-sensitivity_runs.R first.")

## ---- Table: scenario rows + change vs. production ----
production <- read_csv(here(table_dir, "04_mst_summary.csv"), show_col_types = FALSE) %>%
  dplyr::select(species, prod_n_components = n_components, prod_euclid_m = total_mst_euclidean_dist_m,
                prod_eca_frac = eca_frac)
results <- results %>%
  left_join(production, by = "species") %>%
  mutate(delta_n_components = n_components - prod_n_components,
         delta_total_mst_euclid_km = (total_mst_euclidean_dist_m - prod_euclid_m) / 1000,
         delta_eca_frac = eca_frac - prod_eca_frac) %>%
  dplyr::select(-starts_with("prod_")) %>%
  mutate(
    input = if_else(parameter == "dispersal_mult", "Dispersal distance", "Resistance contrast"),
    level = if_else(parameter == "dispersal_mult", paste0("x", signif(value, 3)), paste0("m = ", value)),
    ## both inputs are multipliers: < 1 = shorter dispersal / flatter contrast
    direction = if_else(value < 1, "below production", "above production")
  )
write_csv(results, here(table_dir, "09_sensitivity_results.csv"))
results %>%
  dplyr::select(species, scenario, edge_jaccard, corridor_agreement, surface_schoener_d,
                delta_n_components, minutes) %>%
  mutate(across(where(is.double), ~ round(.x, 3))) %>%
  print(n = Inf)

if (!dir.exists(sensitivity_out_dir)) dir.create(sensitivity_out_dir, recursive = TRUE)
scen <- results %>% filter(scenario != "identity_check")

## Species silhouette (Data/Phylopic/{sp}.png) as a grob for the
## figure title; NULL if there is none. Read with terra (no png/magick needed):
## 1-2 bands = grey (+ alpha), 3-4 bands = RGB (+ alpha).
phylopic_grob <- function(sp) {
  f <- file.path(data_dir, "Phylopic", paste0(sp, ".png"))
  if (!file.exists(f)) return(NULL)
  v <- terra::as.array(suppressWarnings(terra::rast(f)))
  nb <- dim(v)[3]
  rgb_v <- if (nb >= 3) v[, , 1:3] else array(v[, , 1], c(dim(v)[1:2], 3))
  alpha <- if (nb %in% c(2, 4)) v[, , nb] else array(255, dim(v)[1:2])
  alpha[is.na(alpha)] <- 0; rgb_v[is.na(rgb_v)] <- 0
  cols <- rgb(rgb_v[, , 1], rgb_v[, , 2], rgb_v[, , 3], alpha, maxColorValue = 255)
  grid::rasterGrob(matrix(cols, nrow = dim(v)[1]), interpolate = TRUE)
}

for (sp in unique(scen$species)) {
  d <- scen %>% filter(species == sp)
  surface <- d$surface[1]
  sp_dir <- file.path(sensitivity_dir, sp)

  ## ---- Tornado-style figure ----
  surface_label <- if (surface == "overlap") "SDM overlap changed (1 - D)" else "Corridor density changed (1 - D)"
  plot_df <- d %>%
    mutate(`MST edges changed (1 - Jaccard)` = 1 - edge_jaccard,
           `Corridor location changed` = 1 - corridor_agreement,
           !!surface_label := 1 - surface_schoener_d,
           `Change in MST components` = delta_n_components) %>%
    pivot_longer(c(`MST edges changed (1 - Jaccard)`, `Corridor location changed`,
                   all_of(surface_label), `Change in MST components`),
                 names_to = "metric", values_to = "change") %>%
    mutate(metric = factor(metric, levels = c("MST edges changed (1 - Jaccard)", "Corridor location changed",
                                              surface_label, "Change in MST components")),
           ## one row per scenario: dispersal levels on top, contrast levels below
           scenario_label = paste(if_else(parameter == "dispersal_mult", "Dispersal", "Contrast"), level))
  label_order <- d %>%
    mutate(scenario_label = paste(if_else(parameter == "dispersal_mult", "Dispersal", "Contrast"), level)) %>%
    arrange(parameter != "dispersal_mult", value) %>% pull(scenario_label) %>% unique()
  plot_df$scenario_label <- factor(plot_df$scenario_label, levels = rev(label_order))
  n_disp <- sum(d$parameter == "dispersal_mult")
  ## 2 x 2 panels (one per measure) on a page-width figure, so the text stays readable
  p <- ggplot(plot_df, aes(change, scenario_label, colour = direction)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_hline(yintercept = length(label_order) - n_disp + 0.5, colour = "grey75", linetype = "dashed") +
    geom_segment(aes(x = 0, xend = change, yend = scenario_label), alpha = 0.6) +
    geom_point(size = 2.2) +
    facet_wrap(vars(metric), ncol = 2, scales = "free_x") +
    scale_colour_manual(values = c("below production" = "#2166ac", "above production" = "#b2182b")) +
    labs(x = "Change vs. production", y = NULL, colour = NULL,
         title = paste(sp, "-- one-at-a-time sensitivity (each scenario changes one input)")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title.position = "plot")
  fig_h <- 2.2 + 0.5 * length(label_order)
  icon <- phylopic_grob(sp)
  if (!is.null(icon)) {
    ## silhouette in the top-left corner, title shifted right to make room
    icon_in <- 0.55  # icon box (in)
    p <- p + theme(plot.title = element_text(size = 11, margin = margin(l = icon_in * 72 + 6, t = 8, b = 10)))
    p <- cowplot::ggdraw(p) +
      cowplot::draw_grob(icon, x = 0.005, y = 1 - 0.01, width = icon_in / 8, height = icon_in / fig_h,
                         hjust = 0, vjust = 1)
  } else {
    p <- p + theme(plot.title = element_text(size = 11))
  }
  ggsave(file.path(fig_dir, paste0("09_", sp, "_sensitivity_tornado.png")), p,
         width = 8, height = fig_h, dpi = 300, bg = "white")

  ## ---- Edge retention on production's MST lines ----
  scenario_keys <- lapply(d$scenario, function(s) {
    unique(mst_edge_keys(read_csv(file.path(sp_dir, paste0(s, "_mst_edges.csv")), show_col_types = FALSE)))
  })
  kept_share <- function(keys, use) {
    if (!any(use)) return(rep(NA_real_, length(keys)))
    rowMeans(do.call(cbind, lapply(scenario_keys[use], function(k) keys %in% k)))
  }
  prod_lines <- st_read(file.path(mst_lines_dir, paste0(sp, "_mst_lines.gpkg")), quiet = TRUE)
  keys <- mst_edge_keys(prod_lines)
  is_disp <- d$parameter == "dispersal_mult"
  prod_lines$retention_all <- kept_share(keys, rep(TRUE, nrow(d)))
  prod_lines$retention_dispersal <- kept_share(keys, is_disp)
  prod_lines$retention_contrast <- kept_share(keys, !is_disp)
  prod_lines$n_scenarios <- nrow(d)
  st_write(prod_lines, file.path(sensitivity_out_dir, paste0(sp, "_mst_edge_retention.gpkg")),
           delete_dsn = TRUE, quiet = TRUE)

  ## ---- Surface variability across scenarios ----
  files <- file.path(sp_dir, paste0(d$scenario, "_", surface, "_1km.tif"))
  files <- files[file.exists(files)]
  if (length(files) >= 2) {
    s <- subst(rast(files), NA, 0)  # no corridor in a scenario = 0 there
    m <- mean(s)
    sdv <- stdev(s, pop = FALSE)
    out <- c(m, sdv, ifel(m > 0, sdv / m, NA))
    names(out) <- c("mean", "sd", "cv")
    writeRaster(out, file.path(sensitivity_out_dir, paste0(sp, "_", surface, "_variability.tif")), overwrite = TRUE)
  }
  message(sp, ": figure, edge-retention map and ", surface, " variability written (", nrow(d), " scenarios)")
}
