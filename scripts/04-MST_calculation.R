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
source(file.path("scripts", "setup_script.R"))

plot_diagnostics <- TRUE   # set FALSE to skip the per-species QA plot() calls (e.g. for an unattended batch run)
# NOTE: unlike script 03, this script doesn't skip already-done species --
# building the MST from an existing {sp}_PA_costs.csv is cheap (igraph, not
# costDist()), and always recomputing keeps every species' row present in
# the aggregated mst_summary.csv written at the end (skipping would silently
# drop that species' row on a partial rerun).

# Load required files
pa_nodes <- st_read(here(output_dir, "pa_ctroidnodes.shp"))
dkernel_params <- read_csv(here(interm_dir, "dispersal_kernel_params.csv"))
summary_rows <- list()

## Node coordinates, keyed by site_id -- used below to build line geometries
node_coords <- st_coordinates(pa_nodes)
rownames(node_coords) <- pa_nodes$site_id
node_coords_tbl <- node_coords |>
  as.data.frame() |>
  tibble::rownames_to_column("site_id") |>
  as_tibble()

# Load in species list for analysis
specieslist <- readxl::read_excel(here(data_dir, "specieslist.xlsx"))

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to quickly test/compare just those species instead (must already
## have {sp}_PA_costs.csv from script 03 for each one).
species_subset <- c("MAAM", "DOOR")  # e.g. NULL for all species
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
    edges %>% dplyr::select(from, to, weight = cost_dist, dispersal_probability),
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
  
  summary_rows[[sp]] <- tibble(
    species              = sp,
    kernel_p_threshold   = p_threshold,
    n_nodes              = vcount(g),
    n_candidate_edges    = ecount(g),
    n_mst_edges          = ecount(mst_g),
    n_components         = n_components,
    total_mst_cost       = sum(E(mst_g)$weight),
    weakest_link_prob    = weakest_link_prob  # lowest dispersal probability among the MST's own edges
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
  edges_coords <- mst_edges |>
    left_join(node_coords_tbl, by = c("from" = "site_id")) |>
    rename(x_from = X, y_from = Y) |>
    left_join(node_coords_tbl, by = c("to" = "site_id")) |>
    rename(x_to = X, y_to = Y)

  mst_lines <- edges_coords |>
    rowwise() |>
    mutate(geometry = st_sfc(
      st_linestring(matrix(c(x_from, x_to, y_from, y_to), ncol = 2)),
      crs = st_crs(pa_nodes)
    )) |>
    ungroup() |>
    st_as_sf()

  ## save to a gpkg
  st_write(mst_lines, lines_out_path, delete_dsn = TRUE)
  message("  -> saved ", lines_out_path)
}

summary_tbl <- bind_rows(summary_rows)
summary_tbl
write_csv(summary_tbl, file.path(dir_mst, "mst_summary.csv"))

