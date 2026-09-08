library(Seurat)

# Cluster cells together

cluster_genes <- function(
    pbmc,
    resolution = 0.5,
    algorithm = 4,
    random_seed = 1234,
    graph_name = NULL
) {
  
  # Start total function benchmark
  total_start_time <- Sys.time()
  
  # -----------------------------------------------
  # 1. Validate the Seurat object
  # -----------------------------------------------
  
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
  
  # -----------------------------------------------
  # 2. Validate resolution
  # -----------------------------------------------
  
  if (
    !is.numeric(resolution) ||
    length(resolution) != 1 ||
    is.na(resolution) ||
    !is.finite(resolution) ||
    resolution <= 0
  ) {
    stop("'resolution' must be one positive numeric value.")
  }
  
  # -----------------------------------------------
  # 3. Validate clustering algorithm
  # -----------------------------------------------
  
  valid_algorithms <- 1:4
  
  if (
    !is.numeric(algorithm) ||
    length(algorithm) != 1 ||
    is.na(algorithm) ||
    !is.finite(algorithm) ||
    algorithm != as.integer(algorithm) ||
    !algorithm %in% valid_algorithms
  ) {
    stop("'algorithm' must be a whole number from 1 through 4.")
  }
  
  # -----------------------------------------------
  # 4. Validate random seed
  # -----------------------------------------------
  
  if (
    !is.numeric(random_seed) ||
    length(random_seed) != 1 ||
    is.na(random_seed) ||
    !is.finite(random_seed) ||
    random_seed != as.integer(random_seed)
  ) {
    stop("'random_seed' must be one whole number.")
  }
  
  # -----------------------------------------------
  # 5. Check for a neighbor graph
  # -----------------------------------------------
  
  available_graphs <- Graphs(pbmc)
  
  if (length(available_graphs) == 0) {
    stop(
      "No neighbor graph was found in the Seurat object. ",
      "Run FindNeighbors() before running cluster_genes()."
    )
  }
  
  # If graph_name is supplied, confirm that it exists
  if (!is.null(graph_name)) {
    
    if (
      !is.character(graph_name) ||
      length(graph_name) != 1 ||
      is.na(graph_name) ||
      trimws(graph_name) == ""
    ) {
      stop(
        "'graph_name' must be NULL or one non-empty character value."
      )
    }
    
    if (!graph_name %in% available_graphs) {
      stop(
        "The graph '",
        graph_name,
        "' was not found. Available graphs: ",
        paste(available_graphs, collapse = ", "),
        "."
      )
    }
  }
  
  # -----------------------------------------------
  # 6. Run and benchmark clustering
  # -----------------------------------------------
  
  message(
    "Clustering cells at resolution ",
    resolution,
    " using algorithm ",
    algorithm,
    "..."
  )
  
  clustering_start_time <- Sys.time()
  
  pbmc <- tryCatch(
    {
      if (is.null(graph_name)) {
        FindClusters(
          object = pbmc,
          resolution = resolution,
          algorithm = as.integer(algorithm),
          random.seed = as.integer(random_seed),
          verbose = FALSE
        )
      } else {
        FindClusters(
          object = pbmc,
          graph.name = graph_name,
          resolution = resolution,
          algorithm = as.integer(algorithm),
          random.seed = as.integer(random_seed),
          verbose = FALSE
        )
      }
    },
    error = function(e) {
      elapsed_before_error <- as.numeric(
        difftime(
          Sys.time(),
          clustering_start_time,
          units = "secs"
        )
      )
      
      stop(
        "Cell clustering failed after ",
        round(elapsed_before_error, 3),
        " seconds: ",
        conditionMessage(e)
      )
    }
  )
  
  clustering_end_time <- Sys.time()
  
  clustering_elapsed_seconds <- as.numeric(
    difftime(
      clustering_end_time,
      clustering_start_time,
      units = "secs"
    )
  )
  
  # -----------------------------------------------
  # 7. Confirm that clusters were created
  # -----------------------------------------------
  
  if (!"seurat_clusters" %in% colnames(pbmc[[]])) {
    stop(
      "FindClusters() finished, but no 'seurat_clusters' ",
      "column was added to the cell metadata."
    )
  }
  
  cluster_counts <- table(pbmc$seurat_clusters)
  
  if (length(cluster_counts) == 0) {
    stop("Clustering completed, but no clusters were identified.")
  }
  
  if (anyNA(pbmc$seurat_clusters)) {
    warning(
      sum(is.na(pbmc$seurat_clusters)),
      " cells were not assigned to a cluster."
    )
  }
  
  # -----------------------------------------------
  # 8. Calculate total runtime
  # -----------------------------------------------
  
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
      "FindClusters",
      "Entire function"
    ),
    elapsed_seconds = round(
      c(
        clustering_elapsed_seconds,
        total_elapsed_seconds
      ),
      3
    )
  )
  
  # -----------------------------------------------
  # 9. Store parameters and benchmark
  # -----------------------------------------------
  
  pbmc@misc$pipeline_parameters$clustering <- list(
    resolution = resolution,
    algorithm = algorithm,
    algorithm_name = switch(
      as.character(algorithm),
      "1" = "Original Louvain",
      "2" = "Louvain with multilevel refinement",
      "3" = "SLM",
      "4" = "Leiden"
    ),
    random_seed = random_seed,
    graph_name = graph_name,
    number_of_clusters = length(cluster_counts)
  )
  
  pbmc@misc$benchmarks$cluster_genes <- benchmark
  
  # -----------------------------------------------
  # 10. Report results
  # -----------------------------------------------
  
  message(
    "Clustering completed successfully. Number of clusters: ",
    length(cluster_counts)
  )
  
  message(
    "FindClusters runtime: ",
    round(clustering_elapsed_seconds, 2),
    " seconds."
  )
  
  message(
    "Total function runtime: ",
    round(total_elapsed_seconds, 2),
    " seconds."
  )
  
  print(cluster_counts)
  print(benchmark)
  
  # -----------------------------------------------
  # 11. Return the updated Seurat object
  # -----------------------------------------------
  
  return(pbmc)
}