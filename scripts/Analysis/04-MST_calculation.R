## ------------------------------- ##
## Script name: 04-MST_calculation
##
## Purpose of script: Re-apply the dispersal cutoff on cost-distance itself
#    (cost-distance can exceed euclidean, so some prefiltered
#    pairs may still fail the real dispersal limit), then build
#    the MST (or minimum spanning FOREST if disconnected) for each species
#    using the cost-distance weights to determine the connected edges of the MST.
#    MST generated using the PA pair candidates and cost distances with Prim's greedy algorithm.
#    (The cutoff + MST + summary metrics are build_mst() in pipeline_functions.R.)
#
#    Produces individual MST gpkg for each species as output
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))

plot_diagnostics <- FALSE   # set FALSE to skip the per-species QA plot() calls (e.g. for an unattended batch run)
# NOTE: unlike script 03, this script doesn't skip already-done species --
# building the MST from an existing {sp}_PA_costs.csv is cheap (igraph, not
# costDist()), and always recomputing keeps every species' row present in
# the aggregated 04_mst_summary.csv written at the end (skipping would silently
# drop that species' row on a partial rerun).

# Load required files
pa_nodes <- st_read(here(pa_patches_dir, "pa_ctroidnodes.shp"))
## pa_patches: the actual PA polygons (site_id-keyed), needed to draw each
## MST edge's line from patch boundary to patch boundary rather than
## centroid to centroid -- see the "Convert to a vector layer" section below.
pa_patches <- read_sf(here(pa_patches_dir, "pa_patches.gpkg"))
pa_area_km2 <- pa_area_lookup(pa_patches)  # for the PA-area columns of 04_mst_summary.csv
dkernel_params <- read_csv(here(table_dir, "01_dispersal_kernel_params.csv"))
summary_rows <- list()
edge_stats_rows <- list()

## Node coordinates, keyed by site_id -- used below for the diagnostic
## plot() layout only (line geometries are now built from pa_patches
## polygon boundaries, not centroids -- see "Convert to a vector layer" below).
node_coords <- st_coordinates(pa_nodes)
rownames(node_coords) <- pa_nodes$site_id

# Load in species list for analysis
specieslist <- readxl::read_excel(here(tabular_dir, "specieslist.xlsx"))

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to quickly test/compare just those species instead (must already
## have {sp}_PA_costs.csv from script 03 for each one).
species_subset <- c("MAAM", "ASFL", "DOOR", "DRPI", "AQCH")  # e.g. NULL for all species
if (!is.null(species_subset)) {
  specieslist <- specieslist %>% filter(species %in% species_subset)
}

## Create the MST output directories (mst_dir/mst_lines_dir are set in setup_script.R)
if (!dir.exists(mst_dir)) dir.create(mst_dir)
if (!dir.exists(mst_lines_dir)) dir.create(mst_lines_dir)

## Calculate MST forest for the desired PAs and species
for (i in seq_len(nrow(specieslist))) {
  sp <- specieslist$species[i]
  lines_out_path <- file.path(mst_lines_dir, paste0(sp, "_mst_lines.gpkg"))

  message("Building MST for: ", sp)

  alpha <- dkernel_params |> filter(species == sp) |> pull(alpha)
  p_threshold <- dkernel_params |> filter(species == sp) |> pull(kernel_p_threshold)

  edges_path <- file.path(costdist_dir, paste0(sp, "_PA_costs.csv"))
  if (!file.exists(edges_path)) {
    warning("  missing ", edges_path, " -- run cost-distance script first. Skipping.")
    next
  }

  pair_costs <- read_csv(edges_path, show_col_types = FALSE)

  ## Cutoff on euclidean-distance kernel probability, MST by cost_dist,
  ## plus this species' mst_summary / mst_edge_stats rows
  mst_result <- build_mst(pair_costs, pa_nodes$site_id, alpha, p_threshold, sp, pa_area_km2 = pa_area_km2)
  mst_g <- mst_result$mst
  mst_edges <- mst_result$edges

  out_edges_path <- file.path(mst_dir, paste0(sp, "_mst_edges.csv"))
  write_csv(mst_edges, out_edges_path)

  out_graph_path <- file.path(mst_dir, paste0(sp, "_mst.rds"))
  saveRDS(mst_g, out_graph_path)

  message("  -> ", nrow(mst_edges), " MST edge(s) saved to ", out_edges_path)

  summary_rows[[sp]] <- mst_result$summary
  edge_stats_rows[[sp]] <- mst_result$edge_stats

  weakest_link_prob <- mst_result$summary$weakest_link_prob
  if (!is.na(weakest_link_prob) && weakest_link_prob < 5 * p_threshold) {
    message("  note: weakest MST edge (p = ", signif(weakest_link_prob, 3),
            ") is close to the inclusion threshold (p = ", p_threshold,
            ") -- this connection is fragile and worth a manual look.")
  }

  ### Diagnostic plots (optional QA) ############################################
  if (plot_diagnostics) {
    ## This is to just check the node attributes
    plot(mst_g, main = sp)
    ## Plot with spatially correct node locations (reorder to match the graph's vertex order)
    layout_matrix <- node_coords[V(mst_g)$name, ]
    plot(mst_g, layout = layout_matrix, vertex.size = 6,
         vertex.label.cex = 0.6, edge.width = 2, main = sp)
  }

  ### Convert to a vector layer for saving #######################################
  ## Boundary-to-boundary (edge-to-edge), not centroid-to-centroid: each
  ## MST edge's line runs between the nearest points on the two patches'
  ## actual boundaries (sf::st_nearest_points(), same endpoint method
  ## script 05 uses). Still a STRAIGHT line, not the routed least-cost
  ## corridor -- that's script 05's job -- this just anchors the line to
  ## the patches' edges instead of their centroids.
  mst_lines <- mst_edge_lines(mst_edges, pa_patches)

  ## save to a gpkg
  st_write(mst_lines, lines_out_path, delete_dsn = TRUE)
  message("  -> saved ", lines_out_path)
}

summary_tbl <- bind_rows(summary_rows)
summary_tbl
write_csv(summary_tbl, file.path(table_dir, "04_mst_summary.csv"))

edge_stats_tbl <- bind_rows(edge_stats_rows)
edge_stats_tbl %>% mutate(across(where(is.double), ~ signif(.x, 3))) %>% print(n = Inf, width = Inf)
write_csv(edge_stats_tbl, file.path(table_dir, "04_mst_edge_stats.csv"))
