library(Seurat)

# Plot feature expression counts from a selected assay layer

plot_raw_counts <- function(
    pbmc,
    feature_list,
    layer_type = "counts",
    assay_name = NULL,
    log_scale = TRUE,
    pt_size = 0
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
  # 2. Validate or select the assay
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
  # 3. Validate the layer
  # -----------------------------------------------
  
  if (
    !is.character(layer_type) ||
    length(layer_type) != 1 ||
    is.na(layer_type) ||
    trimws(layer_type) == ""
  ) {
    stop("'layer_type' must be one non-empty character value.")
  }
  
  available_layers <- tryCatch(
    {
      Layers(pbmc[[assay_name]])
    },
    error = function(e) {
      stop(
        "Unable to retrieve layers from assay '",
        assay_name,
        "': ",
        conditionMessage(e)
      )
    }
  )
  
  if (length(available_layers) == 0) {
    stop(
      "The assay '",
      assay_name,
      "' does not contain any data layers."
    )
  }
  
  if (!layer_type %in% available_layers) {
    stop(
      "The layer '",
      layer_type,
      "' was not found in assay '",
      assay_name,
      "'. Available layers: ",
      paste(available_layers, collapse = ", "),
      "."
    )
  }
  
  # -----------------------------------------------
  # 4. Validate features
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
  
  assay_features <- rownames(pbmc[[assay_name]])
  
  missing_features <- setdiff(
    feature_list,
    assay_features
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
  
  # Confirm that requested features exist in the selected layer
  layer_features <- tryCatch(
    {
      rownames(
        LayerData(
          object = pbmc,
          assay = assay_name,
          layer = layer_type
        )
      )
    },
    error = function(e) {
      stop(
        "Unable to read layer '",
        layer_type,
        "' from assay '",
        assay_name,
        "': ",
        conditionMessage(e)
      )
    }
  )
  
  missing_from_layer <- setdiff(
    feature_list,
    layer_features
  )
  
  if (length(missing_from_layer) > 0) {
    stop(
      "The following features exist in the assay but not in layer '",
      layer_type,
      "': ",
      paste(missing_from_layer, collapse = ", "),
      "."
    )
  }
  
  # -----------------------------------------------
  # 5. Validate plotting arguments
  # -----------------------------------------------
  
  if (
    !is.logical(log_scale) ||
    length(log_scale) != 1 ||
    is.na(log_scale)
  ) {
    stop("'log_scale' must be either TRUE or FALSE.")
  }
  
  if (
    !is.numeric(pt_size) ||
    length(pt_size) != 1 ||
    is.na(pt_size) ||
    !is.finite(pt_size) ||
    pt_size < 0
  ) {
    stop("'pt_size' must be one non-negative numeric value.")
  }
  
  # -----------------------------------------------
  # 6. Create and benchmark the violin plot
  # -----------------------------------------------
  
  message(
    "Creating violin plot for ",
    length(feature_list),
    " feature(s) using the '",
    layer_type,
    "' layer..."
  )
  
  plot_start_time <- Sys.time()
  
  counts_plot <- tryCatch(
    {
      VlnPlot(
        object = pbmc,
        features = feature_list,
        assay = assay_name,
        layer = layer_type,
        log = log_scale,
        pt.size = pt_size
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
        "Violin-plot creation failed after ",
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
  
  # Display the plot when running interactively
  print(counts_plot)
  
  # -----------------------------------------------
  # 7. Calculate total runtime
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
      "Violin plot creation",
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
  # 8. Store parameters and benchmark
  # -----------------------------------------------
  
  pbmc@misc$pipeline_parameters$raw_count_plot <- list(
    assay = assay_name,
    layer = layer_type,
    features = feature_list,
    log_scale = log_scale,
    pt_size = pt_size
  )
  
  pbmc@misc$benchmarks$plot_raw_counts <- benchmark
  
  # -----------------------------------------------
  # 9. Report results
  # -----------------------------------------------
  
  message("Violin plot created successfully.")
  
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
  # 10. Return results
  # -----------------------------------------------
  
  return(
    list(
      pbmc = pbmc,
      plot = counts_plot,
      benchmark = benchmark,
      assay = assay_name,
      layer = layer_type,
      features = feature_list
    )
  )
}