library(Seurat)

# Run UMAP, create a plot, and save the Seurat object

run_umap <- function(
    pbmc,
    output_file_name = "../output/pbmc_tutorial.rds",
    input_reduction = "pca",
    output_reduction = "umap",
    dims_low = 1,
    dims_high = 20,
    random_seed = 1234
) {
  
  # Start total function benchmark
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
  # 2. Validate reduction names
  # -------------------------------------------------
  
  validate_text_argument <- function(value, argument_name) {
    if (
      !is.character(value) ||
      length(value) != 1 ||
      is.na(value) ||
      trimws(value) == ""
    ) {
      stop(
        "'", argument_name,
        "' must be one non-empty character value."
      )
    }
  }
  
  validate_text_argument(input_reduction, "input_reduction")
  validate_text_argument(output_reduction, "output_reduction")
  validate_text_argument(output_file_name, "output_file_name")
  
  available_reductions <- Reductions(pbmc)
  
  if (!input_reduction %in% available_reductions) {
    stop(
      "The input reduction '",
      input_reduction,
      "' was not found. Available reductions: ",
      if (length(available_reductions) == 0) {
        "none"
      } else {
        paste(available_reductions, collapse = ", ")
      },
      "."
    )
  }
  
  if (input_reduction == output_reduction) {
    stop(
      "'input_reduction' and 'output_reduction' must have ",
      "different names. For example, use 'pca' and 'umap'."
    )
  }
  
  if (output_reduction %in% available_reductions) {
    warning(
      "The reduction '",
      output_reduction,
      "' already exists and will be replaced."
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
      "'dims_low' cannot be greater than 'dims_high'. ",
      "Received dims_low = ",
      dims_low,
      " and dims_high = ",
      dims_high,
      "."
    )
  }
  
  available_dimension_count <- tryCatch(
    {
      ncol(
        Embeddings(
          object = pbmc,
          reduction = input_reduction
        )
      )
    },
    error = function(e) {
      stop(
        "Unable to retrieve embeddings from the '",
        input_reduction,
        "' reduction: ",
        conditionMessage(e)
      )
    }
  )
  
  if (dims_high > available_dimension_count) {
    stop(
      "Dimension ",
      dims_high,
      " was requested, but the '",
      input_reduction,
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
  # 4. Validate random seed
  # -------------------------------------------------
  
  if (
    !is.numeric(random_seed) ||
    length(random_seed) != 1 ||
    is.na(random_seed) ||
    !is.finite(random_seed) ||
    random_seed != as.integer(random_seed)
  ) {
    stop("'random_seed' must be one whole number.")
  }
  
  # -------------------------------------------------
  # 5. Validate the output path
  # -------------------------------------------------
  
  output_extension <- tools::file_ext(output_file_name)
  
  if (tolower(output_extension) != "rds") {
    stop(
      "'output_file_name' must end in '.rds'. Received: ",
      output_file_name
    )
  }
  
  output_directory <- dirname(output_file_name)
  
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
  # 6. Run and benchmark UMAP
  # -------------------------------------------------
  
  message(
    "Running UMAP using ",
    input_reduction,
    " dimensions ",
    dims_low,
    " through ",
    dims_high,
    "..."
  )
  
  umap_start_time <- Sys.time()
  
  pbmc <- tryCatch(
    {
      RunUMAP(
        object = pbmc,
        reduction = input_reduction,
        dims = selected_dims,
        reduction.name = output_reduction,
        reduction.key = "UMAP_",
        seed.use = as.integer(random_seed),
        verbose = FALSE
      )
    },
    error = function(e) {
      elapsed_before_error <- as.numeric(
        difftime(
          Sys.time(),
          umap_start_time,
          units = "secs"
        )
      )
      
      stop(
        "UMAP failed after ",
        round(elapsed_before_error, 3),
        " seconds: ",
        conditionMessage(e)
      )
    }
  )
  
  umap_end_time <- Sys.time()
  
  umap_elapsed_seconds <- as.numeric(
    difftime(
      umap_end_time,
      umap_start_time,
      units = "secs"
    )
  )
  
  # Confirm that UMAP was created
  if (!output_reduction %in% Reductions(pbmc)) {
    stop(
      "RunUMAP() finished, but the '",
      output_reduction,
      "' reduction was not added to the Seurat object."
    )
  }
  
  umap_embeddings <- Embeddings(
    object = pbmc,
    reduction = output_reduction
  )
  
  if (nrow(umap_embeddings) != ncol(pbmc)) {
    stop(
      "The number of UMAP embeddings does not match ",
      "the number of cells in the Seurat object."
    )
  }
  
  # -------------------------------------------------
  # 7. Create and benchmark the UMAP plot
  # -------------------------------------------------
  
  plot_start_time <- Sys.time()
  
  umap_plot <- tryCatch(
    {
      DimPlot(
        object = pbmc,
        reduction = output_reduction
      )
    },
    error = function(e) {
      stop(
        "UMAP was created, but plotting failed: ",
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
  
  print(umap_plot)
  
  # -------------------------------------------------
  # 8. Store parameters before saving
  # -------------------------------------------------
  
  pbmc@misc$pipeline_parameters$umap <- list(
    input_reduction = input_reduction,
    output_reduction = output_reduction,
    dimensions = selected_dims,
    random_seed = random_seed
  )
  
  # -------------------------------------------------
  # 9. Save and benchmark the Seurat object
  # -------------------------------------------------
  
  save_start_time <- Sys.time()
  
  tryCatch(
    {
      saveRDS(
        object = pbmc,
        file = output_file_name
      )
    },
    error = function(e) {
      elapsed_before_error <- as.numeric(
        difftime(
          Sys.time(),
          save_start_time,
          units = "secs"
        )
      )
      
      stop(
        "UMAP completed, but the Seurat object could not be saved after ",
        round(elapsed_before_error, 3),
        " seconds: ",
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
  # 10. Calculate benchmark results
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
      "RunUMAP",
      "UMAP plot creation",
      "Save Seurat object",
      "Entire function"
    ),
    elapsed_seconds = round(
      c(
        umap_elapsed_seconds,
        plot_elapsed_seconds,
        save_elapsed_seconds,
        total_elapsed_seconds
      ),
      3
    )
  )
  
  # Store benchmark in the object
  pbmc@misc$benchmarks$run_umap <- benchmark
  
  # Because the benchmark was added after the first saveRDS(),
  # save the object again so the benchmark is included in the file.
  tryCatch(
    {
      saveRDS(
        object = pbmc,
        file = output_file_name
      )
    },
    error = function(e) {
      stop(
        "The analysis completed, but the final object with benchmark ",
        "metadata could not be saved: ",
        conditionMessage(e)
      )
    }
  )
  
  # -------------------------------------------------
  # 11. Report results
  # -------------------------------------------------
  
  message("UMAP completed successfully.")
  message("Seurat object saved to: ", output_file_name)
  
  message(
    "RunUMAP runtime: ",
    round(umap_elapsed_seconds, 2),
    " seconds."
  )
  
  message(
    "Plot creation runtime: ",
    round(plot_elapsed_seconds, 2),
    " seconds."
  )
  
  message(
    "Object-saving runtime: ",
    round(save_elapsed_seconds, 2),
    " seconds."
  )
  
  message(
    "Total function runtime: ",
    round(total_elapsed_seconds, 2),
    " seconds."
  )
  
  print(benchmark)
  
  # -------------------------------------------------
  # 12. Return the updated Seurat object
  # -------------------------------------------------
  
  return(pbmc)
}