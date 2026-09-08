library(Seurat)

# Plot selected gene expression on a UMAP

umap_features <- function(
    pbmc,
    feature_list,
    reduction_type = "umap",
    assay_name = NULL,
    slot_type = "data",
    plot_columns = NULL,
    point_size = NULL,
    plot_high_expression_last = TRUE
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
  
  if (nrow(pbmc) == 0) {
    stop("The Seurat object does not contain any features.")
  }
  
  # -----------------------------------------------
  # 2. Validate the UMAP reduction
  # -----------------------------------------------
  
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
      ". Run RunUMAP() before calling umap_features()."
    )
  }
  
  reduction_embeddings <- tryCatch(
    {
      Embeddings(
        object = pbmc,
        reduction = reduction_type
      )
    },
    error = function(e) {
      stop(
        "Unable to retrieve embeddings from reduction '",
        reduction_type,
        "': ",
        conditionMessage(e)
      )
    }
  )
  
  if (ncol(reduction_embeddings) < 2) {
    stop(
      "The reduction '",
      reduction_type,
      "' must contain at least two dimensions."
    )
  }
  
  if (nrow(reduction_embeddings) != ncol(pbmc)) {
    stop(
      "The number of cells in the '",
      reduction_type,
      "' reduction does not match the Seurat object."
    )
  }
  
  # -----------------------------------------------
  # 3. Validate or select the assay
  # -----------------------------------------------
  
  available_assays <- Assays(pbmc)
  
  if (length(available_assays) == 0) {
    stop("The Seurat object does not contain any assays.")
  }
  
  if (is.null(assay_name)) {
    assay_name <- DefaultAssay(pbmc)
  }
  
  if (
    !is.character(assay_name) ||
    length(assay_name) != 1 ||
    is.na(assay_name) ||
    trimws(assay_name) == ""
  ) {
    stop(
      "'assay_name' must be NULL or one non-empty character value."
    )
  }
  
  if (!assay_name %in% available_assays) {
    stop(
      "The assay '",
      assay_name,
      "' was not found. Available assays: ",
      paste(available_assays, collapse = ", "),
      "."
    )
  }
  
  # -----------------------------------------------
  # 4. Validate the data slot
  # -----------------------------------------------
  
  if (
    !is.character(slot_type) ||
    length(slot_type) != 1 ||
    is.na(slot_type) ||
    trimws(slot_type) == ""
  ) {
    stop("'slot_type' must be one non-empty character value.")
  }
  
  valid_slots <- c("counts", "data", "scale.data")
  
  if (!slot_type %in% valid_slots) {
    stop(
      "'slot_type' must be one of: ",
      paste(valid_slots, collapse = ", "),
      "."
    )
  }
  
  # -----------------------------------------------
  # 5. Validate feature names
  # -----------------------------------------------
  
  if (
    missing(feature_list) ||
    is.null(feature_list) ||
    !is.character(feature_list) ||
    length(feature_list) == 0 ||
    anyNA(feature_list) ||
    any(trimws(feature_list) == "")
  ) {
    stop(
      "'feature_list' must be a character vector containing ",
      "at least one valid feature name."
    )
  }
  
  feature_list <- unique(feature_list)
  
  available_features <- rownames(pbmc[[assay_name]])
  
  missing_features <- setdiff(
    feature_list,
    available_features
  )
  
  if (length(missing_features) > 0) {
    stop(
      "The following features were not found in assay '",
      assay_name,
      "': ",
      paste(missing_features, collapse = ", "),
      "."
    )
  }
  
  # Confirm that the features can be retrieved from the selected slot
  tryCatch(
    {
      FetchData(
        object = pbmc,
        vars = feature_list,
        assay = assay_name,
        slot = slot_type
      )
    },
    error = function(e) {
      stop(
        "The requested features could not be retrieved from the '",
        slot_type,
        "' slot of assay '",
        assay_name,
        "': ",
        conditionMessage(e)
      )
    }
  )
  
  # -----------------------------------------------
  # 6. Validate plotting arguments
  # -----------------------------------------------
  
  if (!is.null(plot_columns)) {
    if (
      !is.numeric(plot_columns) ||
      length(plot_columns) != 1 ||
      is.na(plot_columns) ||
      !is.finite(plot_columns) ||
      plot_columns <= 0 ||
      plot_columns != as.integer(plot_columns)
    ) {
      stop(
        "'plot_columns' must be NULL or one positive whole number."
      )
    }
  }
  
  if (!is.null(point_size)) {
    if (
      !is.numeric(point_size) ||
      length(point_size) != 1 ||
      is.na(point_size) ||
      !is.finite(point_size) ||
      point_size <= 0
    ) {
      stop(
        "'point_size' must be NULL or one positive numeric value."
      )
    }
  }
  
  if (
    !is.logical(plot_high_expression_last) ||
    length(plot_high_expression_last) != 1 ||
    is.na(plot_high_expression_last)
  ) {
    stop(
      "'plot_high_expression_last' must be TRUE or FALSE."
    )
  }
  
  # -----------------------------------------------
  # 7. Create and benchmark the FeaturePlot
  # -----------------------------------------------
  
  message(
    "Plotting ",
    length(feature_list),
    " feature(s) on the '",
    reduction_type,
    "' reduction..."
  )
  
  plot_start_time <- Sys.time()
  
  feature_plot <- tryCatch(
    {
      FeaturePlot(
        object = pbmc,
        features = feature_list,
        reduction = reduction_type,
        assay = assay_name,
        slot = slot_type,
        ncol = plot_columns,
        pt.size = point_size,
        order = plot_high_expression_last,
        combine = TRUE
      )
    },
    error = function(e) {
      elapsed_before_error <- as.numeric(
        difftime(
          Sys.time(),
          plot_start_time,
          units = "secs"
        )
      )
      
      stop(
        "UMAP feature plotting failed after ",
        round(elapsed_before_error, 3),
        " seconds: ",
        conditionMessage(e)
      )
    }
  )
  
  plot_end_time <- Sys.time()
  
  plot_elapsed_seconds <- as.numeric(
    difftime(
      plot_end_time,
      plot_start_time,
      units = "secs"
    )
  )
  
  print(feature_plot)
  
  # -----------------------------------------------
  # 8. Calculate the total runtime
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
      "FeaturePlot creation",
      "Entire function"
    ),
    elapsed_seconds = round(
      c(
        plot_elapsed_seconds,
        total_elapsed_seconds
      ),
      3
    )
  )
  
  # -----------------------------------------------
  # 9. Store parameters and benchmark
  # -----------------------------------------------
  
  pbmc@misc$pipeline_parameters$umap_features <- list(
    features = feature_list,
    reduction = reduction_type,
    assay = assay_name,
    slot = slot_type,
    plot_columns = plot_columns,
    point_size = point_size,
    plot_high_expression_last = plot_high_expression_last
  )
  
  pbmc@misc$benchmarks$umap_features <- benchmark
  
  # -----------------------------------------------
  # 10. Report results
  # -----------------------------------------------
  
  message("UMAP feature plot created successfully.")
  
  message(
    "Plot creation runtime: ",
    round(plot_elapsed_seconds, 2),
    " seconds."
  )
  
  message(
    "Total function runtime: ",
    round(total_elapsed_seconds, 2),
    " seconds."
  )
  
  print(benchmark)
  
  # -----------------------------------------------
  # 11. Return results
  # -----------------------------------------------
  
  return(
    list(
      pbmc = pbmc,
      feature_plot = feature_plot,
      benchmark = benchmark,
      features = feature_list
    )
  )
}