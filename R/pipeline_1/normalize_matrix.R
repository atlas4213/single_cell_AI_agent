# Normalize Dataset
# Creator: Dave Hill


suppressPackageStartupMessages({
  library(Seurat)
})


#' Normalize a Seurat Gene-Expression Matrix
#'
#' Normalizes a Seurat object using one of Seurat's supported normalization
#' methods and records execution times.
#'
#' @param pbmc A Seurat object.
#' @param normalization_method Normalization method. Must be
#'   `"LogNormalize"`, `"CLR"`, or `"RC"`.
#' @param scale_factor_number Numeric scale factor used during normalization.
#'   Default is `10000`.
#' @param benchmark Logical. If `TRUE`, print benchmark results.
#'
#' @return A normalized Seurat object. Benchmark results are stored in
#'   `pbmc@misc$benchmark$normalize_matrix`.
#'
#' @export
normalize_matrix <- function(
    pbmc,
    normalization_method = "LogNormalize",
    scale_factor_number = 10000,
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
  
  valid_methods <- c(
    "LogNormalize",
    "CLR",
    "RC"
  )
  
  if (
    length(normalization_method) != 1L ||
    !is.character(normalization_method) ||
    is.na(normalization_method) ||
    !normalization_method %in% valid_methods
  ) {
    stop(
      "`normalization_method` must be one of: ",
      paste(valid_methods, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  
  if (
    length(scale_factor_number) != 1L ||
    !is.numeric(scale_factor_number) ||
    is.na(scale_factor_number) ||
    !is.finite(scale_factor_number) ||
    scale_factor_number <= 0
  ) {
    stop(
      "`scale_factor_number` must be a single number greater than zero.",
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
  
  validation_seconds <-
    proc.time()[["elapsed"]] - validation_start
  
  # ---------------------------------------------------------------
  # 2. Normalize the count matrix
  # ---------------------------------------------------------------
  
  normalization_start <- proc.time()[["elapsed"]]
  
  normalized_pbmc <- tryCatch(
    {
      Seurat::NormalizeData(
        object = pbmc,
        normalization.method = normalization_method,
        scale.factor = scale_factor_number,
        verbose = FALSE
      )
    },
    error = function(e) {
      stop(
        "Normalization failed using method `",
        normalization_method,
        "`.\nOriginal error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  normalization_seconds <-
    proc.time()[["elapsed"]] - normalization_start
  
  # ---------------------------------------------------------------
  # 3. Validate the normalized result
  # ---------------------------------------------------------------
  
  result_validation_start <- proc.time()[["elapsed"]]
  
  if (!inherits(normalized_pbmc, "Seurat")) {
    stop(
      "Normalization did not return a Seurat object.",
      call. = FALSE
    )
  }
  
  if (
    nrow(normalized_pbmc) != nrow(pbmc) ||
    ncol(normalized_pbmc) != ncol(pbmc)
  ) {
    warning(
      "The dimensions of the Seurat object changed during normalization.",
      call. = FALSE
    )
  }
  
  result_validation_seconds <-
    proc.time()[["elapsed"]] - result_validation_start
  
  total_seconds <-
    proc.time()[["elapsed"]] - total_start
  
  # ---------------------------------------------------------------
  # 4. Create and store benchmark results
  # ---------------------------------------------------------------
  
  benchmark_results <- data.frame(
    step = c(
      "Validate inputs",
      "Normalize matrix",
      "Validate normalized result",
      "Total normalize_matrix function"
    ),
    elapsed_seconds = round(
      c(
        validation_seconds,
        normalization_seconds,
        result_validation_seconds,
        total_seconds
      ),
      digits = 3
    ),
    row.names = NULL
  )
  
  # Convert old benchmark storage into a named list
  if (is.data.frame(normalized_pbmc@misc$benchmark)) {
    normalized_pbmc@misc$benchmark <- list(
      previous_benchmark = normalized_pbmc@misc$benchmark
    )
  }
  
  if (is.null(normalized_pbmc@misc$benchmark)) {
    normalized_pbmc@misc$benchmark <- list()
  }
  
  if (!is.list(normalized_pbmc@misc$benchmark)) {
    stop(
      "`pbmc@misc$benchmark` exists but is not a list.",
      call. = FALSE
    )
  }
  
  normalized_pbmc@misc$benchmark[["normalize_matrix"]] <-
    benchmark_results
  
  # Store the normalization settings
  normalized_pbmc@misc$normalization_settings <- list(
    method = normalization_method,
    scale_factor = scale_factor_number
  )
  
  message(
    "Normalization completed successfully using ",
    normalization_method,
    ": ",
    nrow(normalized_pbmc),
    " genes x ",
    ncol(normalized_pbmc),
    " cells."
  )
  
  if (benchmark) {
    message("\nNormalization benchmark:")
    print(
      benchmark_results,
      row.names = FALSE
    )
  }
  
  return(normalized_pbmc)
}