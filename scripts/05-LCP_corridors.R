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
#
# Saves one gpkg per species (one LINESTRING feature per successfully
# traced MST edge, all in one merged layer) to:
#   data/intermediate/corridors/{species}_lcp_corridors.gpkg
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
# surface on the crop. The buffer (lcp_buffer_frac below) gives the path
# room to route around a barrier near the straight line between the two
# nodes -- a pure bounding box with no buffer would force the path
# through it. If tracing fails within the initial buffer (e.g. the crop
# cuts off the only viable route), it's retried once with a doubled
# buffer before being skipped.
#
# RESISTANCE -> CONDUCTANCE: gdistance's transition() takes a
# `transitionFunction` that converts a PAIR of adjacent cells' values into
# the conductance (ease of movement) of the edge between them. This
# script uses `function(x) 1 / mean(x)` -- the standard gdistance
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
# dispatch below rebuilds the cluster and marks that batch failed rather
# than letting the whole run abort.
# ============================================================
source(file.path("scripts", "setup_script.R"))

## ---- Sample test: restrict to a couple of example species ----
## Set to NULL to run the full species list; set to a vector of species
## codes to test/compare just those species instead.
species_subset <- c("MAAM", "DOOR")  # e.g. NULL for all species

overwrite <- TRUE  # set TRUE to force recomputation of species already done

## How far beyond each edge's straight-line bounding box to crop the
## resistance raster, as a fraction of that edge's straight-line distance
## (plus a minimum, so very short edges still get a sane amount of room).
## Doubled once automatically if the first attempt fails to find a path.
## Shrinking these is a free way to speed up tracing further, at some risk
## of more edges needing the doubled-buffer retry.
lcp_buffer_frac <- 0.25
lcp_buffer_min_m <- 1000
lcp_buffer_max_m <- 5000  # absolute cap regardless of edge length -- keeps very long edges' corridors from also being very wide

## How many MST edges to hand to the cluster per parLapplyLB() call. Batching
## (rather than one call for a whole species) is what lets the "N/total
## edges traced" progress message below still print periodically.
lcp_batch_size <- 100

n_cores <- max(1, parallel::detectCores() - 10)  # tune down if you want to keep using your machine for other things

specieslist <- readxl::read_excel(here(data_dir, "specieslist.xlsx"))
if (!is.null(species_subset)) {
  specieslist <- specieslist %>% filter(species %in% species_subset)
}

pa_nodes <- st_read(here(output_dir, "pa_ctroidnodes.shp"), quiet = TRUE)
node_coords <- st_coordinates(pa_nodes)
rownames(node_coords) <- pa_nodes$site_id

dir_corridors <- here(interm_dir, "corridors")
if (!dir.exists(dir_corridors)) dir.create(dir_corridors)

sanity_check_done <- FALSE  # only plot once, for the very first successfully traced path overall

## Attempt to trace one edge's LCP on a crop of the raster at `resistance_path`,
## buffered by `buffer_mult` x the edge's normal buffer. Returns NULL (never
## throws) if the crop, transition build, or path tracing fails -- safe to
## call from a parallel worker.
attempt_edge_lcp <- function(from_xy, to_xy, resistance_path, buffer_mult) {
  edge_dist <- sqrt(sum((from_xy - to_xy)^2))
  buffer_dist <- min(max(edge_dist * lcp_buffer_frac, lcp_buffer_min_m), lcp_buffer_max_m) * buffer_mult

  edge_ext <- ext(
    min(from_xy[, 1], to_xy[, 1]) - buffer_dist,
    max(from_xy[, 1], to_xy[, 1]) + buffer_dist,
    min(from_xy[, 2], to_xy[, 2]) - buffer_dist,
    max(from_xy[, 2], to_xy[, 2]) + buffer_dist
  )

  tryCatch({
    resistance <- rast(resistance_path)
    resistance_crop <- crop(resistance, edge_ext)

    ## Mask down to a corridor around the straight line between the two
    ## nodes, not the full rectangle -- transition() excludes NA cells from
    ## the graph it builds, so this keeps graph size proportional to edge
    ## length rather than edge length squared (see MEMORY note above).
    edge_line <- vect(rbind(from_xy, to_xy), type = "lines", crs = crs(resistance))
    edge_corridor <- buffer(edge_line, width = buffer_dist)
    resistance_crop <- mask(resistance_crop, edge_corridor)

    resistance_crop_r <- raster::raster(resistance_crop)

    tr <- transition(resistance_crop_r, transitionFunction = function(x) 1 / mean(x), directions = 8)
    tr <- geoCorrection(tr, type = "c", multpl = FALSE)

    path <- gdistance::shortestPath(tr, from_xy, to_xy, output = "SpatialLines")
    path_sf <- st_as_sf(path)
    st_crs(path_sf) <- crs(resistance)

    list(path = path_sf, resistance_crop = resistance_crop)
  }, error = function(e) NULL)
}

## Trace one MST edge's LCP, retrying once with a doubled buffer on failure.
trace_edge_lcp <- function(from_id, to_id, resistance_path) {
  from_xy <- node_coords[from_id, , drop = FALSE]
  to_xy   <- node_coords[to_id, , drop = FALSE]

  result <- attempt_edge_lcp(from_xy, to_xy, resistance_path, buffer_mult = 1)
  if (is.null(result)) {
    result <- attempt_edge_lcp(from_xy, to_xy, resistance_path, buffer_mult = 2)
  }
  result
}

## One MST edge's full worker task: trace it, then attach its edge metadata.
## References `mst_edges` / `resistance_path` as free variables -- exported
## to the cluster (clusterExport, below) each time they change (i.e. each
## species), and present in the master's own environment for the one edge
## traced serially per species.
lcp_worker <- function(j) {
  from_id <- as.character(mst_edges$from[j])
  to_id   <- as.character(mst_edges$to[j])

  result <- trace_edge_lcp(from_id, to_id, resistance_path)
  if (is.null(result)) return(NULL)

  lcp_sf <- result$path
  lcp_sf$from <- from_id
  lcp_sf$to <- to_id
  lcp_sf$cost_dist <- mst_edges$weight[j]
  lcp_sf$dispersal_probability <- mst_edges$dispersal_probability[j]
  lcp_sf
}

## Builds a fresh cluster with everything that doesn't change across species
## already exported (node_coords, buffer settings, the tracing functions).
## Used for the initial cluster and again if a worker dies mid-run.
setup_cluster <- function() {
  new_cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(new_cl, {
    library(terra); library(sf); library(gdistance); library(raster)
  })
  parallel::clusterExport(new_cl, varlist = c(
    "node_coords", "lcp_buffer_frac", "lcp_buffer_min_m", "lcp_buffer_max_m",
    "attempt_edge_lcp", "trace_edge_lcp", "lcp_worker"
  ))
  new_cl
}

message("Setting up a ", n_cores, "-worker cluster for parallel LCP tracing.")
cl <- setup_cluster()

tryCatch({

  for (i in seq_len(nrow(specieslist))) {
    sp <- specieslist$species[i]
    out_path <- file.path(dir_corridors, paste0(sp, "_lcp_corridors.gpkg"))

    if (file.exists(out_path) && !overwrite) {
      message(sp, ": ", out_path, " already exists -- skipping (set overwrite <- TRUE to redo).")
      next
    }

    edges_path <- here(interm_dir, "mst", paste0(sp, "_mst_edges.csv"))
    if (!file.exists(edges_path)) {
      warning(sp, ": missing ", edges_path, " -- run script 04 first. Skipping.")
      next
    }
    mst_edges <- read_csv(edges_path, show_col_types = FALSE)

    resistance_path <- here(output_dir, "resistance", paste0(sp, "_res.tif"))
    if (!file.exists(resistance_path)) {
      warning(sp, ": missing ", resistance_path, " -- run script 02 first. Skipping.")
      next
    }

    ## `mst_edges` / `resistance_path` change every species -- re-export.
    parallel::clusterExport(cl, varlist = c("mst_edges", "resistance_path"))

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

    ## Remaining edges, dispatched to the cluster in batches (batching keeps
    ## periodic progress messages -- a single call for the whole species
    ## would go silent until every edge in it finished).
    tictoc::tic(paste(sp, "LCP corridor extraction"))
    batches <- split(remaining, ceiling(seq_along(remaining) / lcp_batch_size))
    n_done <- n_edges - length(remaining)

    for (b in batches) {
      ## If a worker dies mid-batch (e.g. an outlier-long edge exhausted its
      ## memory -- see MEMORY note above), parLapplyLB() throws rather than
      ## returning a result. Rebuild the cluster and mark this batch's edges
      ## as failed rather than letting that abort the whole run.
      batch_results <- tryCatch(
        parallel::parLapplyLB(cl, b, lcp_worker),
        error = function(e) {
          warning("  ", sp, ": cluster worker failed on this batch (", conditionMessage(e),
                  ") -- rebuilding the cluster and marking these ", length(b),
                  " edge(s) as failed, then continuing.")
          cl <<- setup_cluster()
          ## setup_cluster() only exports the species-invariant objects --
          ## this species' mst_edges/resistance_path need re-exporting to
          ## the new cluster too, or the next batch fails the same way.
          parallel::clusterExport(cl, varlist = c("mst_edges", "resistance_path"))
          vector("list", length(b))
        }
      )
      edge_paths[b] <- batch_results
      n_done <- n_done + length(b)
      message("  ", sp, ": ", n_done, "/", n_edges, " edge(s) traced")
    }
    tictoc::toc()
    
    n_failed <- sum(vapply(edge_paths, is.null, logical(1)))
    edge_paths <- edge_paths[!vapply(edge_paths, is.null, logical(1))]
    if (length(edge_paths) == 0) {
      warning(sp, ": no MST edges could be traced -- nothing saved.")
      next
    }

    ## Merge every valid edge's line into one vector layer for this species --
    ## this merged file is what script 06 reads.
    corridors <- do.call(rbind, edge_paths)
    st_write(corridors, out_path, delete_dsn = TRUE, quiet = TRUE)
    message("  -> ", nrow(corridors), " of ", n_edges, " edge(s) traced and saved to ", out_path,
            if (n_failed > 0) paste0(" (", n_failed, " edge(s) skipped)") else "")
  }

}, finally = {
  parallel::stopCluster(cl)
  message("Cluster stopped.")
})
