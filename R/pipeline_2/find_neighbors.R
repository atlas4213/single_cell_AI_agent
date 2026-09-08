library(Seurat)

# Find the neighbors for clustering

find_neighbors <- function(
    pbmc,
    dims_high,
    dims_low = 1,
    reduction_type = "pca",
    k_param = 20,
    graph_name = NULL
) {
  
  # Start benchmarking the entire function
  total_start_time <- Sys.time()
  
  # -------------------------------------------------
  # 1. Validate the Seurat object
  # -------------------------------------------------
  
  if (missing(pbmc) || is.null(pbmc)) {
    stop("A Seurat object must be supplied to 'pbmc'.")
  }
  
  if (!inherits(pbmc, "Seurat")) {
    stop(
      "'pbmc' must be a Seurat object. Received: ",
      paste(class(pbmc), collapse = ", ")
    )
  }
  
  if (ncol(pbmc) == 0) {
    stop("The Seurat object does not contain any cells.")
  }
  
  # -------------------------------------------------
  # 2. Validate the reduction
  # -------------------------------------------------
  
  if (
    !is.character(reduction_type) ||
    length(reduction_type) != 1 ||
    is.na(reduction_type) ||
    trimws(reduction_type) == ""
  ) {
    stop(
      "'reduction_type' must be one non-empty character value."
    )
  }
  
  available_reductions <- Reductions(pbmc)
  
  if (!reduction_type %in% available_reductions) {
    stop(
      "The reduction '",
      reduction_type,
      "' was not found. Available reductions: ",
      if (length(available_reductions) == 0) {
        "none"
      } else {
        paste(available_reductions, collapse = ", ")
      },
      "."
    )
  }
  
  # -------------------------------------------------
  # 3. Validate dimensions
  # -------------------------------------------------
  
  validate_dimension <- function(value, argument_name) {
    if (
      !is.numeric(value) ||
      length(value) != 1 ||
      is.na(value) ||
      !is.finite(value) ||
      value <= 0 ||
      value != as.integer(value)
    ) {
      stop(
        "'", argument_name,
        "' must be one positive whole number."
      )
    }
  }
  
  validate_dimension(dims_low, "dims_low")
  validate_dimension(dims_high, "dims_high")
  
  if (dims_low > dims_high) {
    stop(
      "'dims_low' cannot be greater than 'dims_high'."
    )
  }
  
  available_dimension_count <- tryCatch(
    {
      ncol(
        Embeddings(
          object = pbmc,
          reduction = reduction_type
        )
      )
    },
    error = function(e) {
      stop(
        "Unable to retrieve the dimensional-reduction embeddings: ",
        conditionMessage(e)
      )
    }
  )
  
  if (dims_high > available_dimension_count) {
    stop(
      "Dimension ",
      dims_high,
      " was requested, but the '",
      reduction_type,
      "' reduction contains only ",
      available_dimension_count,
      " dimensions."
    )
  }
  
  selected_dims <- seq.int(
    from = as.integer(dims_low),
    to = as.integer(dims_high)
  )
  
  # -------------------------------------------------
  # 4. Validate k_param
  # -------------------------------------------------
  
  if (
    !is.numeric(k_param) ||
    length(k_param) != 1 ||
    is.na(k_param) ||
    !is.finite(k_param) ||
    k_param <= 0 ||
    k_param != as.integer(k_param)
  ) {
    stop("'k_param' must be one positive whole number.")
  }
  
  if (k_param >= ncol(pbmc)) {
    stop(
      "'k_param' must be smaller than the number of cells. ",
      "The object contains ",
      ncol(pbmc),
      " cells."
    )
  }
  
  # -------------------------------------------------
  # 5. Validate optional graph names
  # -------------------------------------------------
  
  if (!is.null(graph_name)) {
    if (
      !is.character(graph_name) ||
      length(graph_name) != 2 ||
      anyNA(graph_name) ||
      any(trimws(graph_name) == "")
    ) {
      stop(
        "'graph_name' must be NULL or a character vector ",
        "containing two names: the NN and SNN graph names."
      )
    }
    
    if (anyDuplicated(graph_name)) {
      stop("'graph_name' must contain two unique graph names.")
    }
  }
  
  graphs_before <- Graphs(pbmc)
  
  # -------------------------------------------------
  # 6. Run and benchmark FindNeighbors()
  # -------------------------------------------------
  
  message(
    "Building neighbor graphs using ",
    reduction_type,
    " dimensions ",
    dims_low,
    " through ",
    dims_high,
    "..."
  )
  
  neighbors_start_time <- Sys.time()
  
  pbmc <- tryCatch(
    {
      if (is.null(graph_name)) {
        FindNeighbors(
          object = pbmc,
          reduction = reduction_type,
          dims = selected_dims,
          k.param = as.integer(k_param),
          verbose = FALSE
        )
      } else {
        FindNeighbors(
          object = pbmc,
          reduction = reduction_type,
          dims = selected_dims,
          k.param = as.integer(k_param),
          graph.name = graph_name,
          verbose = FALSE
        )
      }
    },
    error = function(e) {
      stop(
        "Neighbor graph construction failed: ",
        conditionMessage(e)
      )
    }
  )
  
  neighbors_end_time <- Sys.time()
  
  neighbors_elapsed_seconds <- as.numeric(
    difftime(
      neighbors_end_time,
      neighbors_start_time,
      units = "secs"
    )
  )
  
  # -------------------------------------------------
  # 7. Confirm graph creation
  # -------------------------------------------------
  
  graphs_after <- Graphs(pbmc)
  
  if (length(graphs_after) == 0) {
    stop(
      "FindNeighbors() completed, but no neighbor graphs ",
      "were found in the Seurat object."
    )
  }
  
  if (!is.null(graph_name)) {
    missing_graphs <- setdiff(graph_name, graphs_after)
    
    if (length(missing_graphs) > 0) {
      stop(
        "The following requested graphs were not created: ",
        paste(missing_graphs, collapse = ", "),
        "."
      )
    }
  }
  
  newly_created_graphs <- setdiff(
    graphs_after,
    graphs_before
  )
  
  # -------------------------------------------------
  # 8. Calculate total runtime
  # -------------------------------------------------
  
  total_end_time <- Sys.time()
  
  total_elapsed_seconds <- as.numeric(
    difftime(
      total_end_time,
      total_start_time,
      units = "secs"
    )
  )
  
  benchmark <- data.frame(
    step = c(
      "FindNeighbors",
      "Entire function"
    ),
    elapsed_seconds = round(
      c(
        neighbors_elapsed_seconds,
        total_elapsed_seconds
      ),
      3
    )
  )
  
  # -------------------------------------------------
  # 9. Store parameters and benchmark
  # -------------------------------------------------
  
  pbmc@misc$pipeline_parameters$neighbors <- list(
    reduction = reduction_type,
    dimensions = selected_dims,
    k_param = k_param,
    graphs = graphs_after
  )
  
  pbmc@misc$benchmarks$find_neighbors <- benchmark
  
  # -------------------------------------------------
  # 10. Report results
  # -------------------------------------------------
  
  message("Neighbor graph construction completed successfully.")
  
  if (length(newly_created_graphs) > 0) {
    message(
      "New graphs created: ",
      paste(newly_created_graphs, collapse = ", ")
    )
  } else {
    message(
      "Existing neighbor graphs were updated: ",
      paste(graphs_after, collapse = ", ")
    )
  }
  
  message(
    "FindNeighbors runtime: ",
    round(neighbors_elapsed_seconds, 2),
    " seconds."
  )
  
  message(
    "Total function runtime: ",
    round(total_elapsed_seconds, 2),
    " seconds."
  )
  
  print(benchmark)
  
  # -------------------------------------------------
  # 11. Return the updated Seurat object
  # -------------------------------------------------
  
  return(pbmc)
}