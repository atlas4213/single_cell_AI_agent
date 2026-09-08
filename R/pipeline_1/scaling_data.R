# Scale Gene-Expression Data
# Creator: Dave Hill


suppressPackageStartupMessages({
  library(Seurat)
})


#' Scale Gene-Expression Data in a Seurat Object
#'
#' Centers and scales gene-expression values using [Seurat::ScaleData()]
#' and records execution benchmarks.
#'
#' @description
#' Scaling transforms each selected gene so that its average expression
#' across cells is approximately zero and its variance is approximately one.
#' This prevents highly expressed genes from dominating downstream analyses
#' such as principal component analysis.
#'
#' @param pbmc A Seurat object containing normalized expression data.
#'
#' @param features A character vector containing the genes to scale.
#'   If `NULL`, the function uses the variable features stored in the Seurat
#'   object. If no variable features exist, all genes are used. Default is
#'   `NULL`.
#'
#' @param vars_to_regress Optional character vector containing metadata
#'   variables to regress out during scaling, such as `"percent.mt"` or
#'   `"nCount_RNA"`. Default is `NULL`.
#'
#' @param scale_max Maximum value allowed after scaling. Values greater than
#'   this threshold are clipped. Default is `10`.
#'
#' @param verbose Logical. If `TRUE`, display Seurat progress messages.
#'   Default is `FALSE`.
#'
#' @param benchmark Logical. If `TRUE`, print benchmark results.
#'   Default is `TRUE`.
#'
#' @details
#' The input Seurat object should be normalized and should generally have
#' variable features identified before this function is called.
#'
#' Benchmark results are stored in:
#'
#' `pbmc@misc$benchmark$scale_matrix_data`
#'
#' Scaling settings are stored in:
#'
#' `pbmc@misc$scaling_settings`
#'
#' @return A named list containing:
#'
#' * `pbmc`: The scaled Seurat object.
#' * `features_scaled`: The genes included in scaling.
#' * `benchmark`: A data frame containing execution times.
#'
#' @examples
#' \dontrun{
#' scaling_results <- scale_matrix_data(
#'   pbmc = pbmc
#' )
#'
#' pbmc <- scaling_results$pbmc
#'
#' scaling_results$benchmark
#'
#' # Regress out mitochondrial percentage
#' scaling_results <- scale_matrix_data(
#'   pbmc = pbmc,
#'   vars_to_regress = "percent.mt"
#' )
#'
#' pbmc <- scaling_results$pbmc
#' }
#'
#' @seealso
#' [Seurat::ScaleData()]
#' [Seurat::VariableFeatures()]
#' [Seurat::RunPCA()]
#'
#' @export
#' @importFrom Seurat ScaleData VariableFeatures
scale_matrix_data <- function(
    pbmc,
    features = NULL,
    vars_to_regress = NULL,
    scale_max = 10,
    verbose = FALSE,
    benchmark = TRUE
) {
  
  total_start <- proc.time()[["elapsed"]]
  
  # ---------------------------------------------------------------
  # 1. Validate inputs
  # ---------------------------------------------------------------
  
  validation_start <- proc.time()[["elapsed"]]
  
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
  
  if (!is.null(features)) {
    if (
      !is.character(features) ||
      length(features) == 0L ||
      anyNA(features) ||
      any(!nzchar(features))
    ) {
      stop(
        "`features` must be `NULL` or a character vector of gene names.",
        call. = FALSE
      )
    }
  }
  
  if (!is.null(vars_to_regress)) {
    if (
      !is.character(vars_to_regress) ||
      length(vars_to_regress) == 0L ||
      anyNA(vars_to_regress) ||
      any(!nzchar(vars_to_regress))
    ) {
      stop(
        "`vars_to_regress` must be `NULL` or a character vector ",
        "of metadata column names.",
        call. = FALSE
      )
    }
    
    missing_variables <- setdiff(
      vars_to_regress,
      colnames(pbmc@meta.data)
    )
    
    if (length(missing_variables) > 0L) {
      stop(
        "The following regression variables are missing from the ",
        "Seurat metadata: ",
        paste(missing_variables, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  if (
    !is.numeric(scale_max) ||
    length(scale_max) != 1L ||
    is.na(scale_max) ||
    !is.finite(scale_max) ||
    scale_max <= 0
  ) {
    stop(
      "`scale_max` must be a single number greater than zero.",
      call. = FALSE
    )
  }
  
  if (
    !is.logical(verbose) ||
    length(verbose) != 1L ||
    is.na(verbose)
  ) {
    stop(
      "`verbose` must be either TRUE or FALSE.",
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
  
  validation_seconds <-
    proc.time()[["elapsed"]] - validation_start
  
  # ---------------------------------------------------------------
  # 2. Select features to scale
  # ---------------------------------------------------------------
  
  feature_selection_start <- proc.time()[["elapsed"]]
  
  if (is.null(features)) {
    features <- SeuratObject::VariableFeatures(pbmc)
    
    if (length(features) == 0L) {
      warning(
        "No variable features were found. All genes will be scaled. ",
        "This may require considerably more memory.",
        call. = FALSE
      )
      
      features <- rownames(pbmc)
    }
  }
  
  features <- unique(features)
  
  missing_features <- setdiff(
    features,
    rownames(pbmc)
  )
  
  if (length(missing_features) > 0L) {
    warning(
      length(missing_features),
      " requested features were not found and will be excluded.",
      call. = FALSE
    )
    
    features <- intersect(
      features,
      rownames(pbmc)
    )
  }
  
  if (length(features) == 0L) {
    stop(
      "None of the requested features exist in the Seurat object.",
      call. = FALSE
    )
  }
  
  feature_selection_seconds <-
    proc.time()[["elapsed"]] - feature_selection_start
  
  # ---------------------------------------------------------------
  # 3. Scale the expression data
  # ---------------------------------------------------------------
  
  scaling_start <- proc.time()[["elapsed"]]
  
  scaled_pbmc <- tryCatch(
    {
      Seurat::ScaleData(
        object = pbmc,
        features = features,
        vars.to.regress = vars_to_regress,
        scale.max = scale_max,
        verbose = verbose
      )
    },
    error = function(e) {
      stop(
        "Scaling the expression matrix failed.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  scaling_seconds <-
    proc.time()[["elapsed"]] - scaling_start
  
  # ---------------------------------------------------------------
  # 4. Validate the result
  # ---------------------------------------------------------------
  
  result_validation_start <- proc.time()[["elapsed"]]
  
  if (!inherits(scaled_pbmc, "Seurat")) {
    stop(
      "Scaling did not return a Seurat object.",
      call. = FALSE
    )
  }
  
  if (
    nrow(scaled_pbmc) != nrow(pbmc) ||
    ncol(scaled_pbmc) != ncol(pbmc)
  ) {
    warning(
      "The dimensions of the Seurat object changed during scaling.",
      call. = FALSE
    )
  }
  
  result_validation_seconds <-
    proc.time()[["elapsed"]] - result_validation_start
  
  # ---------------------------------------------------------------
  # 5. Store benchmark results
  # ---------------------------------------------------------------
  
  total_seconds <-
    proc.time()[["elapsed"]] - total_start
  
  benchmark_results <- data.frame(
    step = c(
      "Validate inputs",
      "Select scaling features",
      "Scale expression matrix",
      "Validate scaled object",
      "Total scale_matrix_data function"
    ),
    elapsed_seconds = round(
      c(
        validation_seconds,
        feature_selection_seconds,
        scaling_seconds,
        result_validation_seconds,
        total_seconds
      ),
      digits = 3
    ),
    row.names = NULL
  )
  
  if (is.data.frame(scaled_pbmc@misc$benchmark)) {
    scaled_pbmc@misc$benchmark <- list(
      previous_benchmark = scaled_pbmc@misc$benchmark
    )
  }
  
  if (is.null(scaled_pbmc@misc$benchmark)) {
    scaled_pbmc@misc$benchmark <- list()
  }
  
  if (!is.list(scaled_pbmc@misc$benchmark)) {
    stop(
      "`pbmc@misc$benchmark` exists but is not a list.",
      call. = FALSE
    )
  }
  
  scaled_pbmc@misc$benchmark[["scale_matrix_data"]] <-
    benchmark_results
  
  scaled_pbmc@misc$scaling_settings <- list(
    number_of_features = length(features),
    features = features,
    variables_regressed = vars_to_regress,
    scale_max = scale_max
  )
  
  message(
    "Scaling completed successfully for ",
    length(features),
    " genes in ",
    round(scaling_seconds, 2),
    " seconds."
  )
  
  if (benchmark) {
    message("\nScaling benchmark:")
    print(
      benchmark_results,
      row.names = FALSE
    )
  }
  
  return(
    list(
      pbmc = scaled_pbmc,
      features_scaled = features,
      benchmark = benchmark_results
    )
  )
}