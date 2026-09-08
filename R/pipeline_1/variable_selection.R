# Highly Variable Gene Selection
# Creator: Dave Hill

#' Select Highly Variable Genes in a Seurat Object
#'
#' Identifies highly variable genes using
#' [Seurat::FindVariableFeatures()], retrieves the most variable genes,
#' creates labeled and unlabeled variable-feature plots, saves the combined
#' plot, and records execution benchmarks.
#'
#' @description
#' Highly variable genes exhibit greater cell-to-cell variation than expected
#' based on their average expression. Selecting these genes reduces the number
#' of features used during downstream dimensionality reduction and clustering.
#'
#' The function creates two plots:
#'
#' * A variable-feature plot without gene labels.
#' * A variable-feature plot labeling the requested top variable genes.
#'
#' The two plots are combined and saved as a single file.
#'
#' @param pbmc A Seurat object containing normalized gene-expression data.
#'
#' @param selection_method_type A single character string specifying the
#'   variable-feature selection method. Supported values are:
#'
#'   * `"vst"`: Variance-stabilizing transformation. This is the default and
#'     generally recommended method.
#'   * `"mean.var.plot"`: Selects features based on mean expression and
#'     dispersion.
#'   * `"dispersion"`: Selects features using standardized dispersion.
#'
#' @param n_feature_number A positive integer specifying the number of
#'   variable genes to select. Default is `2000`.
#'
#' @param variable_genes_top A positive integer specifying how many of the
#'   most variable genes should be labeled and returned. Default is `10`.
#'
#' @param output_file A single character string specifying the path and
#'   filename where the combined variable-feature plot will be saved.
#'   The file format is inferred from the extension. Default is
#'   `"./single_cell_variable_genes_plots.pdf"`.
#'
#' @param width A positive numeric value specifying the saved plot width in
#'   inches. Default is `12`.
#'
#' @param height A positive numeric value specifying the saved plot height in
#'   inches. Default is `6`.
#'
#' @param benchmark Logical. If `TRUE`, print execution times for variable
#'   feature selection, gene retrieval, plot creation, plot saving, and the
#'   complete function. Default is `TRUE`.
#'
#' @details
#' The input object should normally be normalized before this function is
#' called. For example, it can be processed using
#' [Seurat::NormalizeData()] or the pipeline's `normalize_matrix()` function.
#'
#' Benchmark results are stored in:
#'
#' `pbmc@misc$benchmark$variable_gene_selection`
#'
#' The variable-gene settings are stored in:
#'
#' `pbmc@misc$variable_gene_settings`
#'
#' If the directory specified in `output_file` does not exist, the function
#' attempts to create it recursively.
#'
#' @return A named list containing:
#'
#' * `pbmc`: The updated Seurat object containing the selected variable
#'   features and stored benchmark results.
#' * `top_variable_genes`: A character vector containing the requested top
#'   variable genes.
#' * `plot`: The combined variable-feature plot.
#' * `plot_file`: The normalized path to the saved plot.
#' * `benchmark`: A data frame containing step names and elapsed execution
#'   times in seconds.
#'
#' @examples
#' \dontrun{
#' # Run variable-feature selection with default settings
#' variable_gene_results <- variable_gene_selection(
#'   pbmc = pbmc
#' )
#'
#' # Retrieve the updated Seurat object
#' pbmc <- variable_gene_results$pbmc
#'
#' # View the most variable genes
#' variable_gene_results$top_variable_genes
#'
#' # View the benchmark results
#' variable_gene_results$benchmark
#'
#' # Run with customized settings
#' variable_gene_results <- variable_gene_selection(
#'   pbmc = pbmc,
#'   selection_method_type = "vst",
#'   n_feature_number = 3000,
#'   variable_genes_top = 20,
#'   output_file = "./results/pbmc_variable_genes.pdf",
#'   width = 14,
#'   height = 7,
#'   benchmark = TRUE
#' )
#'
#' pbmc <- variable_gene_results$pbmc
#' }
#'
#' @seealso
#' [Seurat::FindVariableFeatures()]
#' [Seurat::VariableFeaturePlot()]
#' [Seurat::LabelPoints()]
#' [SeuratObject::VariableFeatures()]
#' [ggplot2::ggsave()]
#'
#' @export
#'
#' @importFrom Seurat FindVariableFeatures VariableFeaturePlot LabelPoints
#' @importFrom SeuratObject VariableFeatures
#' @importFrom ggplot2 ggsave

variable_gene_selection <- function(
    pbmc,
    selection_method_type = "vst",
    n_feature_number = 2000,
    variable_genes_top = 10,
    output_file = "./single_cell_variable_genes_plots.pdf",
    width = 12,
    height = 6,
    benchmark = TRUE
) {
  
  total_start <- proc.time()[["elapsed"]]
  
  # ---------------------------------------------------------------
  # 1. Validate inputs
  # ---------------------------------------------------------------
  
  validation_start <- proc.time()[["elapsed"]]
  
  validate_positive_integer <- function(
    value,
    argument_name
  ) {
    if (
      length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value) ||
      !is.finite(value) ||
      value <= 0 ||
      value %% 1 != 0
    ) {
      stop(
        "`",
        argument_name,
        "` must be a single integer greater than zero.",
        call. = FALSE
      )
    }
  }
  
  if (missing(pbmc) || is.null(pbmc)) {
    stop(
      "A Seurat object must be supplied as `pbmc`.",
      call. = FALSE
    )
  }
  
  if (!inherits(pbmc, "Seurat")) {
    stop(
      "`pbmc` must be a valid Seurat object.",
      call. = FALSE
    )
  }
  
  if (nrow(pbmc) == 0L || ncol(pbmc) == 0L) {
    stop(
      "`pbmc` contains zero genes or zero cells.",
      call. = FALSE
    )
  }
  
  valid_methods <- c(
    "vst",
    "mean.var.plot",
    "dispersion"
  )
  
  if (
    length(selection_method_type) != 1L ||
    !is.character(selection_method_type) ||
    is.na(selection_method_type) ||
    !selection_method_type %in% valid_methods
  ) {
    stop(
      "`selection_method_type` must be one of: ",
      paste(valid_methods, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  
  validate_positive_integer(
    n_feature_number,
    "n_feature_number"
  )
  
  validate_positive_integer(
    variable_genes_top,
    "variable_genes_top"
  )
  
  total_genes <- nrow(pbmc)
  
  if (n_feature_number > total_genes) {
    stop(
      "`n_feature_number` cannot exceed the total number of genes (",
      total_genes,
      ").",
      call. = FALSE
    )
  }
  
  if (variable_genes_top > n_feature_number) {
    stop(
      "`variable_genes_top` cannot be greater than ",
      "`n_feature_number`.",
      call. = FALSE
    )
  }
  
  for (argument in c("width", "height")) {
    value <- get(argument)
    
    if (
      length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value) ||
      !is.finite(value) ||
      value <= 0
    ) {
      stop(
        "`",
        argument,
        "` must be a single number greater than zero.",
        call. = FALSE
      )
    }
  }
  
  if (
    length(output_file) != 1L ||
    !is.character(output_file) ||
    is.na(output_file) ||
    !nzchar(trimws(output_file))
  ) {
    stop(
      "`output_file` must be a valid file path.",
      call. = FALSE
    )
  }
  
  if (
    length(benchmark) != 1L ||
    !is.logical(benchmark) ||
    is.na(benchmark)
  ) {
    stop(
      "`benchmark` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }
  
  # Create the parent output directory if needed
  output_directory <- dirname(output_file)
  
  if (
    !identical(output_directory, ".") &&
    !dir.exists(output_directory)
  ) {
    directory_created <- dir.create(
      output_directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    if (!directory_created && !dir.exists(output_directory)) {
      stop(
        "The output directory could not be created: ",
        output_directory,
        call. = FALSE
      )
    }
  }
  
  validation_seconds <-
    proc.time()[["elapsed"]] - validation_start
  
  # ---------------------------------------------------------------
  # 2. Identify highly variable genes
  # ---------------------------------------------------------------
  
  feature_selection_start <- proc.time()[["elapsed"]]
  
  pbmc <- tryCatch(
    {
      Seurat::FindVariableFeatures(
        object = pbmc,
        selection.method = selection_method_type,
        nfeatures = n_feature_number,
        verbose = FALSE
      )
    },
    error = function(e) {
      stop(
        "Variable-gene selection failed using method `",
        selection_method_type,
        "`.\nOriginal error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  feature_selection_seconds <-
    proc.time()[["elapsed"]] - feature_selection_start
  
  # ---------------------------------------------------------------
  # 3. Retrieve the top variable genes
  # ---------------------------------------------------------------
  
  gene_retrieval_start <- proc.time()[["elapsed"]]
  
  variable_features <- tryCatch(
    {
      SeuratObject::VariableFeatures(pbmc)
    },
    error = function(e) {
      stop(
        "Variable features could not be retrieved.\nOriginal error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  if (length(variable_features) == 0L) {
    stop(
      "No variable genes were identified.",
      call. = FALSE
    )
  }
  
  if (length(variable_features) < variable_genes_top) {
    warning(
      "Only ",
      length(variable_features),
      " variable genes were identified. All identified genes ",
      "will be returned.",
      call. = FALSE
    )
  }
  
  top_genes <- head(
    variable_features,
    n = min(
      variable_genes_top,
      length(variable_features)
    )
  )
  
  gene_retrieval_seconds <-
    proc.time()[["elapsed"]] - gene_retrieval_start
  
  # ---------------------------------------------------------------
  # 4. Create variable-feature plots
  # ---------------------------------------------------------------
  
  plotting_start <- proc.time()[["elapsed"]]
  
  combined_plot <- tryCatch(
    {
      plot1 <- Seurat::VariableFeaturePlot(
        object = pbmc
      )
      
      plot2 <- Seurat::LabelPoints(
        plot = plot1,
        points = top_genes,
        repel = TRUE
      )
      
      plot1 + plot2
    },
    error = function(e) {
      stop(
        "Variable-feature plot creation failed.\nOriginal error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  plotting_seconds <-
    proc.time()[["elapsed"]] - plotting_start
  
  # if (interactive()) {
  #   print(combined_plot)
  # }
  
  # ---------------------------------------------------------------
  # 5. Save the plot
  # ---------------------------------------------------------------
  
  plot_saving_start <- proc.time()[["elapsed"]]
  
  tryCatch(
    {
      ggplot2::ggsave(
        filename = output_file,
        plot = combined_plot,
        width = width,
        height = height,
        units = "in"
      )
    },
    error = function(e) {
      stop(
        "The plot could not be saved to `",
        output_file,
        "`.\nOriginal error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  plot_saving_seconds <-
    proc.time()[["elapsed"]] - plot_saving_start
  
  # ---------------------------------------------------------------
  # 6. Create and store benchmark results
  # ---------------------------------------------------------------
  
  total_seconds <-
    proc.time()[["elapsed"]] - total_start
  
  benchmark_results <- data.frame(
    step = c(
      "Validate inputs",
      "Select variable features",
      "Retrieve top variable genes",
      "Create variable-feature plots",
      "Save variable-feature plots",
      "Total variable_gene_selection function"
    ),
    elapsed_seconds = round(
      c(
        validation_seconds,
        feature_selection_seconds,
        gene_retrieval_seconds,
        plotting_seconds,
        plot_saving_seconds,
        total_seconds
      ),
      digits = 3
    ),
    row.names = NULL
  )
  
  # Convert old benchmark storage into a named list
  if (is.data.frame(pbmc@misc$benchmark)) {
    pbmc@misc$benchmark <- list(
      previous_benchmark = pbmc@misc$benchmark
    )
  }
  
  if (is.null(pbmc@misc$benchmark)) {
    pbmc@misc$benchmark <- list()
  }
  
  if (!is.list(pbmc@misc$benchmark)) {
    stop(
      "`pbmc@misc$benchmark` exists but is not a list.",
      call. = FALSE
    )
  }
  
  pbmc@misc$benchmark[["variable_gene_selection"]] <-
    benchmark_results
  
  # Store parameters used for reproducibility
  pbmc@misc$variable_gene_settings <- list(
    selection_method = selection_method_type,
    requested_features = n_feature_number,
    identified_features = length(variable_features),
    top_gene_count = length(top_genes),
    output_file = normalizePath(
      output_file,
      mustWork = FALSE
    )
  )
  
  message(
    "Variable-gene selection completed successfully: ",
    length(variable_features),
    " variable genes identified using `",
    selection_method_type,
    "`."
  )
  
  if (benchmark) {
    message("\nVariable-gene selection benchmark:")
    print(
      benchmark_results,
      row.names = FALSE
    )
  }
  
  # Store useful results inside the Seurat object
  pbmc@misc$variable_gene_selection <- list(
    top_variable_genes = top_genes,
    selection_method = selection_method_type,
    number_of_features = length(variable_features),
    plot_file = normalizePath(
      output_file,
      mustWork = FALSE
    )
  )
  
  # Benchmark should already be stored here
  pbmc@misc$benchmark[["variable_gene_selection"]] <-
    benchmark_results
  
  return(pbmc)
}