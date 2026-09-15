## ------------------------------- ##
## Script name: 04-MST_calculation
##
## Purpose of script: Re-apply the dispersal cutoff on cost-distance itself
#    (cost-distance can exceed euclidean, so some prefiltered
#    pairs may still fail the real dispersal limit), then build
#    the MST (or minimum spanning FOREST if disconnected) for each species
#    using the cost-distance weights to determine the connected edges of the MST.
#    MST generated using the PA pair candidates and cost distances with Prim's greedy algorithm.
#
#    Produces individual MST gpkg for each species as output
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
suppressMessages(suppressWarnings(source(file.path("scripts", "setup_script.R"))))

plot_diagnostics <- TRUE   # set FALSE to skip the per-species QA plot() calls (e.g. for an unattended batch run)
# NOTE: unlike script 03, this script doesn't skip already-done species --
# building the MST from an existing {sp}_PA_costs.csv is cheap (igraph, not
# costDist()), and always recomputing keeps every species' row present in
# the aggregated mst_summary.csv written at the end (skipping would silently
# drop that species' row on a partial rerun).

# Load required files
pa_nodes <- st_read(here(output_dir, "pa_ctroidnodes.shp"))
## pa_patches: the actual PA polygons (site_id-keyed), needed to draw each
## MST edge's line from patch boundary to patch boundary rather than
## centroid to centroid -- see the "Convert to a vector layer" section below.
pa_patches <- read_sf(here(output_dir, "pa_patches.gpkg"))
dkernel_params <- read_csv(here(table_dir, "dispersal_kernel_params.csv"))
summary_rows <- list()

## Node coordinates, keyed by site_id -- used below for the diagnostic
## plot() layout only (line geometries are now built from pa_patches
## polygon boundaries, not centroids -- see "Convert to a vector layer" below).
node_coords <- st_coordinates(pa_nodes)
rownames(node_coords) <- pa_nodes$site_id

# Load in species list for analysis
specieslist <- readxl::read_excel(here(data_dir, "specieslist.xlsx"))

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to quickly test/compare just those species instead (must already
## have {sp}_PA_costs.csv from script 03 for each one).
species_subset <- c("MAAM", "ASFL", "DOOR")  # e.g. NULL for all species
if (!is.null(species_subset)) {
  specieslist <- specieslist %>% filter(species %in% species_subset)
}

## Set the MST output directories
dir_mst <- here(interm_dir, "mst")
if (!dir.exists(dir_mst)) dir.create(dir_mst)
dir_mst_lines <- here(output_dir, "mst")
if (!dir.exists(dir_mst_lines)) dir.create(dir_mst_lines)

## Calculate MST forest for the desired PAs and species
for (i in seq_len(nrow(specieslist))) {
  sp <- specieslist$species[i]
  lines_out_path <- file.path(dir_mst_lines, paste0(sp, "_mst_lines.gpkg"))

  message("Building MST for: ", sp)

  alpha <- dkernel_params |> filter(species == sp) |> pull(alpha)
  p_threshold <- dkernel_params |> filter(species == sp) |> pull(kernel_p_threshold)
  
  edges_path <- file.path(paste0(interm_dir, "/cost_distances/", sp,"_PA_costs.csv"))
  if (!file.exists(edges_path)) {
    warning("  missing ", edges_path, " -- run cost-distance script first. Skipping.")
    next
  }
  
  edges <- read_csv(edges_path, show_col_types = FALSE) %>%
    filter(!is.na(cost_dist)) %>%
    ## SHOULD WE USE THE CUTOFF FROM EUCL. DISTANCE? COST DIST DOESNT REALLY MAKE SENSE...CAN BE EITHER ONE
    mutate(dispersal_probability = dispersal_probability(euclidean_dist, alpha)) %>%
    filter(dispersal_probability >= p_threshold)
  
  g <- graph_from_data_frame(
    edges %>% dplyr::select(from, to, weight = cost_dist, euclidean_dist, dispersal_probability),
    directed = FALSE,
    vertices = pa_nodes$site_id
  )
  
  n_components <- components(g)$no
  if (n_components > 1) {
    message("  ", n_components, " disconnected components -- result is a ",
            "minimum spanning FOREST, not a single tree.")
  }
  
  ## Generate the MST using the graph of PA nodes and the weighted edges (of cost)
  # The minimum spanning tree will have the smallest sum of edge weights (lowest costs)
  mst_g <- mst(g, weights = E(g)$weight, algorithm = "prim")
  mst_edges <- as_data_frame(mst_g, what = "edges") %>% as_tibble()
  
  out_edges_path <- file.path(dir_mst, paste0(sp, "_mst_edges.csv"))
  write_csv(mst_edges, out_edges_path)
  
  out_graph_path <- file.path(dir_mst, paste0(sp, "_mst.rds"))
  saveRDS(mst_g, out_graph_path)
  
  message("  -> ", nrow(mst_edges), " MST edge(s) saved to ", out_edges_path)
  
  weakest_link_prob <- if (nrow(mst_edges) > 0) min(mst_edges$dispersal_probability) else NA_real_
  mean_link_prob <- if (nrow(mst_edges) > 0) mean(mst_edges$dispersal_probability) else NA_real_

  ## detour factor: how much farther the cheapest cost-path is than a
  ## straight line, per MST edge -- flags edges where resistance forces a
  ## long detour relative to geography alone. NA (not Inf) for edges whose
  ## patches already touch (euclidean_dist = 0), since the ratio is
  ## undefined there.
  detour_factor <- ifelse(mst_edges$euclidean_dist > 0,
                           mst_edges$weight / mst_edges$euclidean_dist, NA_real_)
  n_detour_valid <- sum(!is.na(detour_factor))
  mst_degree <- degree(mst_g)
  ## cost-weighted diameter: worst-case cumulative cost between the two
  ## most distant connected PAs (NA, not 0, when there are no MST edges).
  mst_diameter_cost <- if (nrow(mst_edges) > 0) diameter(mst_g, weights = E(mst_g)$weight) else NA_real_

  summary_rows[[sp]] <- tibble(
    species                    = sp,
    kernel_p_threshold         = p_threshold,
    n_nodes                    = vcount(g),
    n_candidate_edges          = ecount(g),
    n_isolated_nodes           = sum(degree(g) == 0),  # nodes with no candidate pair within dispersal range
    n_mst_edges                = ecount(mst_g),
    n_components               = n_components,
    total_mst_cost             = sum(E(mst_g)$weight),
    mean_mst_edge_cost         = if (nrow(mst_edges) > 0) mean(mst_edges$weight) else NA_real_,
    median_mst_edge_cost       = if (nrow(mst_edges) > 0) median(mst_edges$weight) else NA_real_,
    max_mst_edge_cost          = if (nrow(mst_edges) > 0) max(mst_edges$weight) else NA_real_,
    total_mst_euclidean_dist_m = if (nrow(mst_edges) > 0) sum(mst_edges$euclidean_dist) else NA_real_,
    mean_detour_factor         = if (n_detour_valid > 0) mean(detour_factor, na.rm = TRUE) else NA_real_,
    max_detour_factor          = if (n_detour_valid > 0) max(detour_factor, na.rm = TRUE) else NA_real_,
    mst_diameter_cost          = mst_diameter_cost,
    max_node_degree            = if (length(mst_degree) > 0) max(mst_degree) else NA_real_,
    weakest_link_prob          = weakest_link_prob,  # lowest dispersal probability among the MST's own edges
    mean_link_prob             = mean_link_prob
  )
  
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
  edge_lines <- Map(function(from_id, to_id) {
    origin_poly <- pa_patches[as.character(pa_patches$site_id) == as.character(from_id), ]
    dest_poly   <- pa_patches[as.character(pa_patches$site_id) == as.character(to_id), ]
    st_nearest_points(origin_poly, dest_poly)
  }, mst_edges$from, mst_edges$to)

  mst_lines <- mst_edges
  st_geometry(mst_lines) <- do.call(c, edge_lines)
  st_crs(mst_lines) <- st_crs(pa_patches)

  ## save to a gpkg
  st_write(mst_lines, lines_out_path, delete_dsn = TRUE)
  message("  -> saved ", lines_out_path)
}

summary_tbl <- bind_rows(summary_rows)
summary_tbl
write_csv(summary_tbl, file.path(table_dir, "mst_summary.csv"))

