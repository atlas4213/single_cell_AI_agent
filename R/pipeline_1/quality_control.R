# Quality Control Script
# Creator: Dave Hill

library(Seurat)
library(ggplot2)
library(patchwork)


#' Calculate and Visualize Single-Cell Quality-Control Metrics
#'
#' Calculates mitochondrial expression percentages, creates QC plots,
#' saves the plots, and records execution times.
#'
#' @param pbmc A Seurat object.
#' @param file_name Character string used to name the output files.
#' @param output_dir Directory where QC plots will be saved.
#'   Default is `"./test_output"`.
#' @param benchmark Logical. If `TRUE`, print benchmark results.
#'   Default is `TRUE`.
#'
#' @return A Seurat object containing the `percent.mt` metadata column.
#'   Benchmark results are stored in `pbmc@misc$benchmark$quality_control`.
#'
#' @examples
#' \dontrun{
#' pbmc <- qc_data(
#'   pbmc = pbmc,
#'   file_name = "pbmc3k",
#'   output_dir = "./test_output"
#' )
#'
#' pbmc@misc$benchmark$quality_control
#' }
#'
#' @export
qc_data <- function(
    pbmc,
    file_name,
    output_dir = "./test_output",
    benchmark = TRUE
) {
  
  total_start <- proc.time()[["elapsed"]]
  
  # DEBUGGING 
  #pbmc = pbmc_raw
  #file_name = "tester"
  #output_dir = "./test_output"
  #benchmark = TRUE
  
  # 1. Validate inputs
  if (!inherits(pbmc, "Seurat")) {
    stop(
      "`pbmc` must be a Seurat object.",
      call. = FALSE
    )
  }
  
  if (
    missing(file_name) ||
    !is.character(file_name) ||
    length(file_name) != 1L ||
    is.na(file_name) ||
    !nzchar(file_name)
  ) {
    stop(
      "`file_name` must be a single, non-empty character string.",
      call. = FALSE
    )
  }
  
  if (
    !is.character(output_dir) ||
    length(output_dir) != 1L ||
    is.na(output_dir) ||
    !nzchar(output_dir)
  ) {
    stop(
      "`output_dir` must be a single, non-empty character string.",
      call. = FALSE
    )
  }
  
  if (
    !is.logical(benchmark) ||
    length(benchmark) != 1L ||
    is.na(benchmark)
  ) {
    stop(
      "`benchmark` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }
  
  # Create the output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  if (!dir.exists(output_dir)) {
    stop(
      "The output directory could not be created: ",
      output_dir,
      call. = FALSE
    )
  }
  
  # 2. Calculate mitochondrial percentage
  mitochondrial_start <- proc.time()[["elapsed"]]
  
  pbmc[["percent.mt"]] <- tryCatch(
    Seurat::PercentageFeatureSet(
      object = pbmc,
      pattern = "^MT-"
    ),
    error = function(e) {
      stop(
        "Failed to calculate mitochondrial percentages.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  mitochondrial_seconds <-
    proc.time()[["elapsed"]] - mitochondrial_start
  
  # 3. Create the violin plot
  violin_start <- proc.time()[["elapsed"]]
  
  qc_violin_plot <- tryCatch(
    Seurat::VlnPlot(
      object = pbmc,
      features = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt"
      ),
      ncol = 3
    ),
    error = function(e) {
      stop(
        "Failed to create the QC violin plot.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  violin_seconds <- proc.time()[["elapsed"]] - violin_start
  
  # Display the plot during an interactive R session
  if (interactive()) {
    print(qc_violin_plot)
  }
  
  # 4. Save the violin plot
  violin_save_start <- proc.time()[["elapsed"]]
  
  violin_file <- file.path(
    output_dir,
    paste0(file_name, "_single_cell_qc_violin_plot.pdf")
  )
  
  tryCatch(
    ggplot2::ggsave(
      filename = violin_file,
      plot = qc_violin_plot,
      width = 12,
      height = 5,
      units = "in"
    ),
    error = function(e) {
      stop(
        "Failed to save the QC violin plot to: ",
        violin_file,
        "\nOriginal error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  violin_save_seconds <-
    proc.time()[["elapsed"]] - violin_save_start
  
  # 5. Create feature-scatter plots
  scatter_start <- proc.time()[["elapsed"]]
  
  scatter_plots <- tryCatch(
    {
      plot1 <- Seurat::FeatureScatter(
        object = pbmc,
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA"
      )
      
      plot2 <- Seurat::FeatureScatter(
        object = pbmc,
        feature1 = "nCount_RNA",
        feature2 = "percent.mt"
      )
      
      plot1 + plot2
    },
    error = function(e) {
      stop(
        "Failed to create the QC feature-scatter plots.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  scatter_seconds <- proc.time()[["elapsed"]] - scatter_start
  
  if (interactive()) {
    print(scatter_plots)
  }
  
  # 6. Save the feature-scatter plots
  scatter_save_start <- proc.time()[["elapsed"]]
  
  scatter_file <- file.path(
    output_dir,
    paste0(file_name, "_single_cell_qc_feature_plots.pdf")
  )
  
  tryCatch(
    ggplot2::ggsave(
      filename = scatter_file,
      plot = scatter_plots,
      width = 12,
      height = 5,
      units = "in"
    ),
    error = function(e) {
      stop(
        "Failed to save the QC feature plots to: ",
        scatter_file,
        "\nOriginal error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  scatter_save_seconds <-
    proc.time()[["elapsed"]] - scatter_save_start
  
  # 7. Calculate total runtime
  total_seconds <- proc.time()[["elapsed"]] - total_start
  
  benchmark_results <- data.frame(
    step = c(
      "Calculate mitochondrial percentage",
      "Create violin plot",
      "Save violin plot",
      "Create feature-scatter plots",
      "Save feature-scatter plots",
      "Total qc_data function"
    ),
    elapsed_seconds = round(
      c(
        mitochondrial_seconds,
        violin_seconds,
        violin_save_seconds,
        scatter_seconds,
        scatter_save_seconds,
        total_seconds
      ),
      digits = 3
    ),
    row.names = NULL
  )
  
  # Store results and paths inside the Seurat object
  pbmc@misc$benchmark$quality_control <- list(
    benchmark_results
    )
  
  pbmc@misc$plots$quality_control <- list(
    violin_plot = qc_violin_plot,
    scatter_plots = scatter_plots
  )
  
  pbmc@misc$output_files$quality_control <- list(
    violin_plot = normalizePath(
      violin_file,
      mustWork = FALSE
    ),
    scatter_plots = normalizePath(
      scatter_file,
      mustWork = FALSE
    )
  )
  
  message(
    "Quality control completed successfully for ",
    ncol(pbmc),
    " cells."
  )
  
  if (benchmark) {
    message("\nQuality-control benchmark results:")
    print(benchmark_results, row.names = FALSE)
  }
  
  return(pbmc)
}