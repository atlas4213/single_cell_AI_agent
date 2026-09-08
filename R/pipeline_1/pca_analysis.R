#' Run Principal Component Analysis and Create PCA Plots
#'
#' Runs principal component analysis on the highly variable genes stored in a
#' Seurat object, identifies genes that contribute strongly to selected
#' principal components, creates PCA loading and cell-distribution plots,
#' saves the plots, and measures PCA execution time.
#'
#' @description
#' Principal component analysis reduces a high-dimensional gene-expression
#' matrix into a smaller set of principal components. Each principal component
#' captures a portion of the variation observed across cells.
#'
#' The function uses the highly variable genes stored in the Seurat object as
#' the PCA input features. These features must be selected before calling this
#' function.
#'
#' @param pbmc A Seurat object containing normalized and scaled
#'   gene-expression data. Highly variable features must already have been
#'   identified and stored in the object.
#'
#' @param nfeature_count A positive integer specifying the number of
#'   highest-loading genes to print and display for each selected principal
#'   component. This does not control how many genes are used to run PCA.
#'   Default is `10`.
#'
#' @param reduction_type A single non-empty character string specifying the
#'   name used to store the dimensional reduction in the Seurat object.
#'   Default is `"pca"`.
#'
#' @param output_name A single non-empty character string used as the prefix
#'   for saved output files. Default is `"single_cell"`.
#'
#' @param output_directory A single character string specifying the directory
#'   where the PCA plots will be saved. If the directory does not exist, the
#'   function attempts to create it recursively. Default is the current
#'   directory (`"."`).
#'
#' @param width A positive numeric value specifying the width of each saved
#'   plot in inches. Default is `12`.
#'
#' @param height A positive numeric value specifying the height of each saved
#'   plot in inches. Default is `5`.
#'
#' @param dims_start A positive integer specifying the first principal
#'   component whose highest-loading genes should be printed. Default is `1`.
#'
#' @param dims_end A positive integer specifying the final principal
#'   component whose highest-loading genes should be printed. It must be
#'   greater than or equal to `dims_start`. Default is `5`.
#'
#' @param dims_influence_1 A positive integer specifying the first principal
#'   component included in the PCA loading plot. Default is `1`.
#'
#' @param dims_influence_2 A positive integer specifying the final principal
#'   component included in the PCA loading plot. It must be greater than or
#'   equal to `dims_influence_1`. Default is `2`.
#'
#' @details
#' Before calling this function, the Seurat object should normally pass through
#' the following processing steps:
#'
#' 1. Normalization with [Seurat::NormalizeData()].
#' 2. Variable-feature selection with [Seurat::FindVariableFeatures()].
#' 3. Scaling with [Seurat::ScaleData()].
#'
#' The PCA is run using all variable features returned by:
#'
#' `SeuratObject::VariableFeatures(pbmc)`
#'
#' The function saves two PDF files:
#'
#' * `<output_name>_pca_loadings.pdf`
#' * `<output_name>_pca_cells.pdf`
#'
#' `nfeature_count` only controls how many influential genes are printed or
#' displayed. It does not determine the number of variable genes supplied to
#' [Seurat::RunPCA()].
#'
#' @return A named list containing:
#'
#' * `pbmc`: The updated Seurat object containing the PCA reduction.
#' * `loading_plot`: A plot showing genes with the strongest PCA loadings.
#' * `cell_plot`: A PCA plot showing the distribution of cells.
#' * `loading_file`: The path to the saved PCA-loading plot.
#' * `cell_file`: The path to the saved PCA cell plot.
#'
#' @examples
#' \dontrun{
#' # Run PCA with the default settings
#' pca_results <- pca_analysis(
#'   pbmc = pbmc
#' )
#'
#' # Retrieve the updated Seurat object
#' pbmc <- pca_results$pbmc
#'
#' # Access the plots
#' pca_results$loading_plot
#' pca_results$cell_plot
#'
#' # Access the saved file paths
#' pca_results$loading_file
#' pca_results$cell_file
#'
#' # Run PCA with customized dimensions and output settings
#' pca_results <- pca_analysis(
#'   pbmc = pbmc,
#'   nfeature_count = 15,
#'   reduction_type = "pca",
#'   output_name = "pbmc3k",
#'   output_directory = "./results/pca",
#'   width = 14,
#'   height = 6,
#'   dims_start = 1,
#'   dims_end = 10,
#'   dims_influence_1 = 1,
#'   dims_influence_2 = 4
#' )
#'
#' pbmc <- pca_results$pbmc
#' }
#'
#' @seealso
#' [Seurat::RunPCA()]
#' [Seurat::VizDimLoadings()]
#' [Seurat::DimPlot()]
#' [Seurat::Embeddings()]
#' [Seurat::NoLegend()]
#' [SeuratObject::VariableFeatures()]
#' [ggplot2::ggsave()]
#'
#' @export
#'
#' @importFrom Seurat RunPCA VizDimLoadings DimPlot Embeddings NoLegend
#' @importFrom SeuratObject VariableFeatures
#' @importFrom ggplot2 ggsave

pca_analysis <- function(
    pbmc,
    nfeature_count = 10,
    reduction_type = "pca",
    output_name = "single_cell",
    output_directory = ".",
    width = 12,
    height = 5,
    dims_start = 1,
    dims_end = 5,
    dims_influence_1 = 1,
    dims_influence_2 = 2
) {
  # Benchmark timing 
  overall_start_time <- Sys.time()
  
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
  
  # -------------------------------------------------
  # 2. Validate text arguments
  # -------------------------------------------------
  
  if (
    !is.character(reduction_type) ||
    length(reduction_type) != 1 ||
    is.na(reduction_type) ||
    trimws(reduction_type) == ""
  ) {
    stop("'reduction_type' must be one non-empty character value.")
  }
  
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
  
  # -------------------------------------------------
  # 3. Validate numeric arguments
  # -------------------------------------------------
  
  numeric_arguments <- list(
    nfeature_count = nfeature_count,
    width = width,
    height = height,
    dims_start = dims_start,
    dims_end = dims_end,
    dims_influence_1 = dims_influence_1,
    dims_influence_2 = dims_influence_2
  )
  
  for (argument_name in names(numeric_arguments)) {
    value <- numeric_arguments[[argument_name]]
    
    if (
      !is.numeric(value) ||
      length(value) != 1 ||
      is.na(value) ||
      !is.finite(value) ||
      value <= 0
    ) {
      stop(
        "'", argument_name,
        "' must be one positive numeric value."
      )
    }
  }
  
  # Dimensions and feature count must be whole numbers
  integer_arguments <- c(
    "nfeature_count",
    "dims_start",
    "dims_end",
    "dims_influence_1",
    "dims_influence_2"
  )
  
  for (argument_name in integer_arguments) {
    value <- numeric_arguments[[argument_name]]
    
    if (value != as.integer(value)) {
      stop("'", argument_name, "' must be a whole number.")
    }
  }
  
  if (dims_start > dims_end) {
    stop("'dims_start' cannot be greater than 'dims_end'.")
  }
  
  if (dims_influence_1 > dims_influence_2) {
    stop(
      "'dims_influence_1' cannot be greater than ",
      "'dims_influence_2'."
    )
  }
  
  # -------------------------------------------------
  # 4. Validate output directory
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
  
  # -------------------------------------------------
  # 5. Retrieve highly variable genes
  # -------------------------------------------------
  
  variable_features <- tryCatch(
    VariableFeatures(object = pbmc),
    error = function(e) {
      stop(
        "Unable to retrieve highly variable genes: ",
        conditionMessage(e)
      )
    }
  )
  
  if (length(variable_features) == 0) {
    stop(
      "No highly variable genes were found. ",
      "Run FindVariableFeatures() before pca_analysis()."
    )
  }
  
  # -------------------------------------------------
  # 6. Run PCA
  # -------------------------------------------------
  
  # Benchmark PCA
  pca_start_time <- Sys.time()
  
  pbmc <- tryCatch(
    RunPCA(
      object = pbmc,
      features = variable_features,
      reduction.name = reduction_type,
      verbose = FALSE
    ),
    error = function(e) {
      stop("PCA failed: ", conditionMessage(e))
    }
  )
  
  pca_end_time <- Sys.time()
  
  pca_elapsed_seconds <- as.numeric(
    difftime(pca_end_time, pca_start_time, units = "secs")
  )
  
  message(
    "PCA runtime: ",
    round(pca_elapsed_seconds, 2),
    " seconds."
  )
  # -------------------------------------------------
  # 7. Validate requested dimensions
  # -------------------------------------------------
  
  available_dimensions <- ncol(
    Embeddings(
      object = pbmc,
      reduction = reduction_type
    )
  )
  
  requested_dimensions <- c(
    dims_start,
    dims_end,
    dims_influence_1,
    dims_influence_2
  )
  
  if (any(requested_dimensions > available_dimensions)) {
    stop(
      "A requested PCA dimension exceeds the ",
      available_dimensions,
      " dimensions available in the '",
      reduction_type,
      "' reduction."
    )
  }
  
  # Print the genes contributing most strongly to the PCs
  print(
    pbmc[[reduction_type]],
    dims = dims_start:dims_end,
    nfeatures = nfeature_count
  )
  
  # -------------------------------------------------
  # 8. Create plots
  # -------------------------------------------------
  
  loading_plot <- tryCatch(
    VizDimLoadings(
      object = pbmc,
      dims = dims_influence_1:dims_influence_2,
      reduction = reduction_type,
      nfeatures = nfeature_count
    ),
    error = function(e) {
      stop(
        "Unable to create the PCA loading plot: ",
        conditionMessage(e)
      )
    }
  )
  
  cell_plot <- tryCatch(
    DimPlot(
      object = pbmc,
      reduction = reduction_type
    ) + NoLegend(),
    error = function(e) {
      stop(
        "Unable to create the PCA cell plot: ",
        conditionMessage(e)
      )
    }
  )
  
  # -------------------------------------------------
  # 9. Build output filenames
  # -------------------------------------------------
  
  loading_filename <- file.path(
    output_directory,
    paste0(output_name, "_pca_loadings.pdf")
  )
  
  cell_filename <- file.path(
    output_directory,
    paste0(output_name, "_pca_cells.pdf")
  )
  
  # -------------------------------------------------
  # 10. Save plots
  # -------------------------------------------------
  
  tryCatch(
    {
      ggsave(
        filename = loading_filename,
        plot = loading_plot,
        width = width,
        height = height
      )
      
      ggsave(
        filename = cell_filename,
        plot = cell_plot,
        width = width,
        height = height
      )
    },
    error = function(e) {
      stop(
        "The plots were created but could not be saved: ",
        conditionMessage(e)
      )
    }
  )
  
  message("PCA analysis completed successfully.")
  message("Loading plot saved to: ", loading_filename)
  message("Cell plot saved to: ", cell_filename)
  
  # -------------------------------------------------
  # 11. Return results
  # -------------------------------------------------
  
  return(
    list(
      pbmc = pbmc,
      loading_plot = loading_plot,
      cell_plot = cell_plot,
      loading_file = loading_filename,
      cell_file = cell_filename
    )
  )
}