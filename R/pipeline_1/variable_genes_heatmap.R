# Create Heatmap for Highly variable genes 
#
#
#

variable_genes_heatmap <- function(
    pbmc,
    output_name = "single_cell",
    output_directory = ".",
    reduction_type = "pca",
    dims = 1,
    cell_count = 500,
    width = 12,
    height = 5
) {
  
  # Start total function benchmark
  overall_start_time <- Sys.time()
  
  # -------------------------------------------------
  # 1. Validate the Seurat object
  # -------------------------------------------------
  
  if (missing(pbmc) || is.null(pbmc)) {
    stop("A Seurat object must be provided to 'pbmc'.")
  }
  
  if (!inherits(pbmc, "Seurat")) {
    stop(
      "'pbmc' must be a Seurat object. Received: ",
      paste(class(pbmc), collapse = ", ")
    )
  }
  
  # -------------------------------------------------
  # 2. Validate text arguments
  # -------------------------------------------------
  
  if (
    !is.character(output_name) ||
    length(output_name) != 1 ||
    is.na(output_name) ||
    trimws(output_name) == ""
  ) {
    stop("'output_name' must be one non-empty character value.")
  }
  
  if (
    !is.character(output_directory) ||
    length(output_directory) != 1 ||
    is.na(output_directory) ||
    trimws(output_directory) == ""
  ) {
    stop("'output_directory' must be one valid directory path.")
  }
  
  if (
    !is.character(reduction_type) ||
    length(reduction_type) != 1 ||
    is.na(reduction_type) ||
    trimws(reduction_type) == ""
  ) {
    stop("'reduction_type' must be one non-empty character value.")
  }
  
  # -------------------------------------------------
  # 3. Validate numeric arguments
  # -------------------------------------------------
  
  if (
    !is.numeric(dims) ||
    length(dims) == 0 ||
    anyNA(dims) ||
    any(!is.finite(dims)) ||
    any(dims <= 0) ||
    any(dims != as.integer(dims))
  ) {
    stop("'dims' must contain positive whole numbers.")
  }
  
  if (anyDuplicated(dims)) {
    stop("'dims' cannot contain duplicate dimensions.")
  }
  
  if (
    !is.numeric(cell_count) ||
    length(cell_count) != 1 ||
    is.na(cell_count) ||
    !is.finite(cell_count) ||
    cell_count <= 0 ||
    cell_count != as.integer(cell_count)
  ) {
    stop("'cell_count' must be one positive whole number.")
  }
  
  if (
    !is.numeric(width) ||
    length(width) != 1 ||
    is.na(width) ||
    !is.finite(width) ||
    width <= 0
  ) {
    stop("'width' must be one positive numeric value.")
  }
  
  if (
    !is.numeric(height) ||
    length(height) != 1 ||
    is.na(height) ||
    !is.finite(height) ||
    height <= 0
  ) {
    stop("'height' must be one positive numeric value.")
  }
  
  # -------------------------------------------------
  # 4. Validate the reduction
  # -------------------------------------------------
  
  available_reductions <- Reductions(pbmc)
  
  if (!reduction_type %in% available_reductions) {
    stop(
      "The reduction '",
      reduction_type,
      "' was not found in the Seurat object. Available reductions: ",
      if (length(available_reductions) == 0) {
        "none"
      } else {
        paste(available_reductions, collapse = ", ")
      },
      ". Run RunPCA() before creating a PCA heatmap."
    )
  }
  
  available_dimensions <- tryCatch(
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
        "Unable to inspect the '",
        reduction_type,
        "' reduction: ",
        conditionMessage(e)
      )
    }
  )
  
  if (max(dims) > available_dimensions) {
    stop(
      "Dimension ",
      max(dims),
      " was requested, but the '",
      reduction_type,
      "' reduction contains only ",
      available_dimensions,
      " dimensions."
    )
  }
  
  # -------------------------------------------------
  # 5. Validate the requested cell count
  # -------------------------------------------------
  
  total_cells <- ncol(pbmc)
  
  if (cell_count > total_cells) {
    warning(
      "'cell_count' is greater than the total number of cells. ",
      "Using all ",
      total_cells,
      " cells instead."
    )
    
    cell_count <- total_cells
  }
  
  # -------------------------------------------------
  # 6. Create the output directory
  # -------------------------------------------------
  
  if (!dir.exists(output_directory)) {
    directory_created <- dir.create(
      output_directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    if (!directory_created) {
      stop(
        "The output directory does not exist and could not be created: ",
        output_directory
      )
    }
  }
  
  if (file.access(output_directory, mode = 2) != 0) {
    stop(
      "The output directory is not writable: ",
      output_directory
    )
  }
  
  # -------------------------------------------------
  # 7. Create the output filename
  # -------------------------------------------------
  
  heatmap_filename <- file.path(
    output_directory,
    paste0(
      output_name,
      "_highly_variable_genes_heatmap.pdf"
    )
  )
  
  # -------------------------------------------------
  # 8. Create and benchmark the heatmap
  # -------------------------------------------------
  
  message(
    "Creating heatmap for ",
    length(dims),
    " dimension(s) and ",
    cell_count,
    " cells..."
  )
  
  plot_start_time <- Sys.time()
  
  heatmap_plot <- tryCatch(
    {
      DimHeatmap(
        object = pbmc,
        dims = dims,
        cells = cell_count,
        reduction = reduction_type,
        balanced = TRUE
      )
    },
    error = function(e) {
      stop(
        "Heatmap creation failed: ",
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
  
  # -------------------------------------------------
  # 9. Save and benchmark the heatmap
  # -------------------------------------------------
  
  save_start_time <- Sys.time()
  
  tryCatch(
    {
      ggsave(
        filename = heatmap_filename,
        plot = heatmap_plot,
        width = width,
        height = height,
        units = "in"
      )
    },
    error = function(e) {
      stop(
        "The heatmap was created but could not be saved to '",
        heatmap_filename,
        "': ",
        conditionMessage(e)
      )
    }
  )
  
  save_end_time <- Sys.time()
  
  save_elapsed_seconds <- as.numeric(
    difftime(
      save_end_time,
      save_start_time,
      units = "secs"
    )
  )
  
  # -------------------------------------------------
  # 10. Calculate the total runtime
  # -------------------------------------------------
  
  overall_end_time <- Sys.time()
  
  total_elapsed_seconds <- as.numeric(
    difftime(
      overall_end_time,
      overall_start_time,
      units = "secs"
    )
  )
  
  benchmark <- data.frame(
    step = c(
      "Heatmap creation",
      "Heatmap saving",
      "Entire function"
    ),
    elapsed_seconds = round(
      c(
        plot_elapsed_seconds,
        save_elapsed_seconds,
        total_elapsed_seconds
      ),
      3
    )
  )
  
  message("Heatmap created successfully.")
  message("Heatmap saved to: ", heatmap_filename)
  message(
    "Total runtime: ",
    round(total_elapsed_seconds, 2),
    " seconds."
  )
  
  print(benchmark)
  
  # -------------------------------------------------
  # 11. Return the results
  # -------------------------------------------------
  
  return(
    list(
      pbmc = pbmc,
      heatmap_plot = heatmap_plot,
      heatmap_file = heatmap_filename,
      benchmark = benchmark
    )
  )
}