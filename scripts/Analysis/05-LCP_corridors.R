## ------------------------------- ##
## Script name: 05-LCP_corridors
##
## Purpose of script: Trace the actual least-cost path (LCP) corridor for each MST
# edge, using that species' combined resistance raster -- NOT a straight
# line between PAs, but the routed path a disperser would actually take
# through the landscape. All of one species' valid edge paths are merged
# into a single vector layer, which is what script 06 reads.
#
# Uses gdistance (already a dependency of this project -- see
# setup_script.R -- and, per the comment in 03-costdistance_calculation.R,
# the approach WY's earlier prototype used for cost-distance too). Edges
# are traced in parallel across CPU cores (base R's `parallel` package --
# ships with R, no install needed) since gdistance::transition() is the
# expensive step and every edge's trace is independent of every other.
# The tracing functions (attempt_edge_lcp(), trace_edge_lcp(), lcp_worker(),
# trace_lcp_corridors()) are in pipeline_functions.R; the crop buffer
# settings (lcp_buffer_frac/min/max) are in setup_script.R.
#
# Saves one gpkg per species (one LINESTRING feature per successfully
# traced MST edge, all in one merged layer) to:
#   Intermediate/05_lcp_corridors/{species}_lcp_corridors.gpkg
##
## Author: Andrew Habrich
##
## Notes ------------------------- ##
# WHY THIS BUILDS A TRANSITION SURFACE PER EDGE, NOT ONCE PER SPECIES:
# gdistance::transition() builds an explicit cell-adjacency graph over
# whatever raster you give it. The full modelling-extent resistance raster
# is ~45,700 x 24,600 cells at 30 m -- over a billion cells -- which is not
# feasible to build one graph over (this is why script 03 uses
# terra::costDist() instead, which flood-fills directly on the raster
# without materializing that graph). So for each MST edge, this script
# crops the resistance raster to a small buffered bounding box around just
# that edge's two PA nodes, and builds a fresh (small, fast) transition
# surface on the crop. The buffer (lcp_buffer_frac) gives the path
# room to route around a barrier near the straight line between the two
# nodes -- a pure bounding box with no buffer would force the path
# through it. If tracing fails within the initial buffer (e.g. the crop
# cuts off the only viable route), it's retried once with a doubled
# buffer before being skipped.
#
# RESISTANCE -> CONDUCTANCE: gdistance's transition() takes a
# `transitionFunction` that converts a PAIR of adjacent cells' values into
# the conductance (ease of movement) of the edge between them. This
# uses `function(x) 1 / mean(x)` -- the standard gdistance
# convention for a resistance/cost surface (higher combined resistance ->
# lower conductance -> less preferred route) -- which is the well-
# documented convention gdistance's own vignette recommends, unlike
# leastcostpath's create_cs(), whose cost-vs-conductance direction isn't
# clearly documented. `geoCorrection(type = "c")` then corrects for the
# actual distance between adjacent/diagonal cells, per gdistance's
# standard shortestPath() workflow.
#
# PARALLELIZATION: `terra::SpatRaster` objects wrap an external C++
# pointer and don't survive being sent to a parallel worker process --
# using one on a worker silently produces garbage, not a clean error. So
# `attempt_edge_lcp()` takes `resistance_path` (a plain string) and loads
# the raster itself on whichever process calls it (master or worker);
# `terra::rast()` opens files lazily via GDAL, so this doesn't mean
# re-reading the whole (huge) file each time -- crop() still only reads
# the small window it needs. A single PSOCK cluster (the only option on
# Windows -- no fork()) is built once and reused across every species.
# The very first edge of the very first species is always traced on the
# master process (not the cluster) purely so the sanity-check plot below
# can render to your interactive session -- a worker process has no way
# to pop a plot into RStudio's graphics device.
#
# MEMORY: a plain bounding-box crop grows with the SQUARE of edge length,
# and gdistance::transition() builds an explicit adjacency graph over
# every non-NA cell in it -- so a handful of unusually long MST edges can
# make a worker try to build a graph over millions of cells, which can
# exhaust memory and kill the worker process outright (visible as
# "error reading/writing to connection" at the master -- a dead worker,
# not a normal R error, so the tryCatch() inside attempt_edge_lcp() can't
# catch it). To keep this bounded, the crop is additionally masked down
# to a buffered corridor around the straight line between the two nodes
# (transition() excludes NA cells from its graph), turning the effective
# graph size from ~O(edge_length^2) into ~O(edge_length); lcp_buffer_max_m
# further caps how wide that corridor gets regardless of edge length. If
# a worker still dies on some outlier edge despite this, the per-batch
# dispatch in trace_lcp_corridors() rebuilds the cluster and marks that
# batch failed rather than letting the whole run abort.
# ============================================================
suppressMessages(suppressWarnings(source(file.path("Scripts", "Analysis", "setup_script.R"))))

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to test/compare just those species instead.
species_subset <- c("MAAM", "ASFL", "DOOR", "DRPI", "AQCH")  # e.g. NULL for all species

overwrite <- TRUE  # set TRUE to force recomputation of species already done

## How many MST edges to hand to the cluster per parLapplyLB() call. Batching
## (rather than one call for a whole species) is what lets the "N/total
## edges traced" progress message still print periodically.
lcp_batch_size <- 100

n_cores <- max(1, parallel::detectCores() - 2)  # tune down if you want to keep using your machine for other things

specieslist <- readxl::read_excel(here(tabular_dir, "specieslist.xlsx"))
if (!is.null(species_subset)) {
  specieslist <- specieslist %>% filter(species %in% species_subset)
}

## pa_patches: sf POLYGON object of protected areas -- loaded as sf (not
## terra::vect) specifically because sf objects serialize cleanly to
## parallel workers, unlike a SpatVector (see setup_lcp_cluster()).
pa_patches <- read_sf(here(pa_patches_dir, "pa_patches.gpkg"))

if (!dir.exists(corridor_dir)) dir.create(corridor_dir)

sanity_check_done <- FALSE  # only plot once, for the very first successfully traced path overall

message("Setting up a ", n_cores, "-worker cluster for parallel LCP tracing.")
cl <- setup_lcp_cluster(n_cores)

tryCatch({

  for (i in seq_len(nrow(specieslist))) {
    sp <- specieslist$species[i]
    out_path <- file.path(corridor_dir, paste0(sp, "_lcp_corridors.gpkg"))

    if (file.exists(out_path) && !overwrite) {
      message(sp, ": ", out_path, " already exists -- skipping (set overwrite <- TRUE to redo).")
      next
    }

    edges_path <- here(mst_dir, paste0(sp, "_mst_edges.csv"))
    if (!file.exists(edges_path)) {
      warning(sp, ": missing ", edges_path, " -- run script 04 first. Skipping.")
      next
    }
    mst_edges <- read_csv(edges_path, show_col_types = FALSE)

    resistance_path <- here(resistance_dir, paste0(sp, "_res.tif"))
    if (!file.exists(resistance_path)) {
      warning(sp, ": missing ", resistance_path, " -- run script 02 first. Skipping.")
      next
    }

    n_edges <- nrow(mst_edges)
    message("Tracing LCP corridors for: ", sp, " (", n_edges, " MST edge(s), ", n_cores, " worker(s))")

    edge_paths <- vector("list", n_edges)
    remaining <- seq_len(n_edges)

    ## Trace edge 1 on the master process (not the cluster) so the sanity
    ## check plot, if this is the first successful edge overall, can render.
    if (!sanity_check_done && n_edges > 0) {
      from_id_1 <- as.character(mst_edges$from[1])
      to_id_1   <- as.character(mst_edges$to[1])
      result <- trace_edge_lcp(from_id_1, to_id_1, resistance_path)
      remaining <- remaining[-1]

      if (!is.null(result)) {
        lcp_sf <- result$path
        lcp_sf$from <- from_id_1
        lcp_sf$to <- to_id_1
        lcp_sf$cost_dist <- mst_edges$weight[1]
        lcp_sf$dispersal_probability <- mst_edges$dispersal_probability[1]
        edge_paths[[1]] <- lcp_sf

        plot(result$resistance_crop, main = paste(sp, "edge 1 -- sanity check"))
        plot(st_geometry(lcp_sf), add = TRUE, col = "red", lwd = 2)
        message("SANITY CHECK: confirm the red traced path above hugs low-resistance ",
                "(dark/low-value) terrain rather than cutting through barriers.")
        sanity_check_done <- TRUE
      } else {
        warning("  ", sp, ": edge 1 could not be traced -- no sanity-check plot to show; ",
                "continuing with the rest of the run.")
      }
    }

    ## Remaining edges, dispatched to the cluster in batches (keep the
    ## returned cluster -- it's rebuilt if a worker died).
    lcp <- trace_lcp_corridors(mst_edges, resistance_path, cl, n_cores,
                               edge_ids = remaining, batch_size = lcp_batch_size, sp = sp)
    cl <- lcp$cl
    edge_paths[remaining] <- lcp$paths[remaining]

    n_failed <- sum(vapply(edge_paths, is.null, logical(1)))
    edge_paths <- edge_paths[!vapply(edge_paths, is.null, logical(1))]
    if (length(edge_paths) == 0) {
      warning(sp, ": no MST edges could be traced -- nothing saved.")
      next
    }

    ## Merge every valid edge's line into one vector layer for this species --
    ## this merged file is what script 06 reads.
    corridors <- do.call(rbind, edge_paths)
    st_write(corridors, out_path, delete_dsn = TRUE, delete_layer = TRUE, quiet = TRUE)
    message("  -> ", nrow(corridors), " of ", n_edges, " edge(s) traced and saved to ", out_path,
            if (n_failed > 0) paste0(" (", n_failed, " edge(s) skipped)") else "")
  }

}, finally = {
  parallel::stopCluster(cl)
  message("Cluster stopped.")
})
