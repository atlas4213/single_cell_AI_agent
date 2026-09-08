# Filter Cell Matrix Script
# Creator: Dave Hill


# Load required library
library(Seurat)


#' Filter Cells Using Quality-Control Thresholds
#'
#' Filters cells based on detected features and mitochondrial expression.
#' Execution times and filtering statistics are stored in the returned
#' Seurat object.
#'
#' @param pbmc A Seurat object containing `nFeature_RNA` and `percent.mt`
#'   metadata columns.
#' @param nfeature_high Maximum number of detected features allowed per cell.
#' @param nfeature_low Minimum number of detected features allowed per cell.
#' @param percent_mt_number Maximum mitochondrial-expression percentage
#'   allowed per cell.
#' @param benchmark Logical. If `TRUE`, print benchmark and filtering results.
#'   Default is `TRUE`.
#'
#' @return A filtered Seurat object. Benchmark results are stored in
#'   `pbmc@misc$benchmark$filter_cells`.
#'
#' @examples
#' \dontrun{
#' filtered_pbmc <- filter_cells(
#'   pbmc = pbmc,
#'   nfeature_high = 2500,
#'   nfeature_low = 200,
#'   percent_mt_number = 5
#' )
#'
#' filtered_pbmc@misc$benchmark$filter_cells
#' }
#'
#' @export
#' @importFrom Seurat subset
filter_cells <- function(
    pbmc,
    nfeature_high,
    nfeature_low,
    percent_mt_number,
    benchmark = TRUE
) {
  
  # Start total execution timer
  total_start <- proc.time()[["elapsed"]]
  
  # Start validation timer
  validation_start <- proc.time()[["elapsed"]]
  
  # 1. Validate the Seurat object
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
  
  # 2. Validate the benchmark argument
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
  
  # 3. Validate feature thresholds
  validate_positive_integer <- function(value, argument_name) {
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
        "` must be a single integer greater than 0.",
        call. = FALSE
      )
    }
  }
  
  validate_positive_integer(
    nfeature_high,
    "nfeature_high"
  )
  
  validate_positive_integer(
    nfeature_low,
    "nfeature_low"
  )
  
  # Mitochondrial percentage can be a decimal, such as 7.5
  if (
    length(percent_mt_number) != 1L ||
    !is.numeric(percent_mt_number) ||
    is.na(percent_mt_number) ||
    !is.finite(percent_mt_number) ||
    percent_mt_number < 0 ||
    percent_mt_number > 100
  ) {
    stop(
      "`percent_mt_number` must be a single number between 0 and 100.",
      call. = FALSE
    )
  }
  
  if (nfeature_low >= nfeature_high) {
    stop(
      "`nfeature_low` must be less than `nfeature_high`.",
      call. = FALSE
    )
  }
  
  # 4. Confirm that the required metadata columns exist
  required_columns <- c(
    "nFeature_RNA",
    "percent.mt"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(pbmc@meta.data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "The Seurat object's metadata is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  validation_seconds <-
    proc.time()[["elapsed"]] - validation_start
  
  # Record the number of cells before filtering
  cells_before <- ncol(pbmc)
  
  # 5. Filter the cells
  filtering_start <- proc.time()[["elapsed"]]
  
  filtered_pbmc <- tryCatch(
    {
      subset(
        x = pbmc,
        subset =
          nFeature_RNA > nfeature_low &
          nFeature_RNA < nfeature_high &
          percent.mt < percent_mt_number
      )
    },
    error = function(e) {
      stop(
        "Cell filtering failed.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  filtering_seconds <-
    proc.time()[["elapsed"]] - filtering_start
  
  # Record filtering statistics
  cells_after <- ncol(filtered_pbmc)
  cells_removed <- cells_before - cells_after
  
  percent_retained <- if (cells_before > 0) {
    (cells_after / cells_before) * 100
  } else {
    0
  }
  
  percent_removed <- if (cells_before > 0) {
    (cells_removed / cells_before) * 100
  } else {
    0
  }
  
  if (cells_after == 0L) {
    stop(
      "All cells were removed during filtering. ",
      "Review `nfeature_low`, `nfeature_high`, and ",
      "`percent_mt_number`.",
      call. = FALSE
    )
  }
  
  # Calculate total execution time
  total_seconds <-
    proc.time()[["elapsed"]] - total_start
  
  # 6. Create benchmark table
  benchmark_results <- data.frame(
    step = c(
      "Validate inputs",
      "Filter cells",
      "Total filter_cells function"
    ),
    elapsed_seconds = round(
      c(
        validation_seconds,
        filtering_seconds,
        total_seconds
      ),
      digits = 3
    ),
    row.names = NULL
  )
  
  # 7. Create filtering-summary table
  filtering_summary <- data.frame(
    metric = c(
      "Cells before filtering",
      "Cells after filtering",
      "Cells removed",
      "Percent retained",
      "Percent removed"
    ),
    value = c(
      cells_before,
      cells_after,
      cells_removed,
      round(percent_retained, 2),
      round(percent_removed, 2)
    ),
    row.names = NULL
  )
  
  # 8. Prepare the benchmark storage structure
  
  # Convert benchmarks created using the old data-frame format
  if (is.data.frame(filtered_pbmc@misc$benchmark)) {
    filtered_pbmc@misc$benchmark <- list(
      load_data = filtered_pbmc@misc$benchmark
    )
  }
  
  # Initialize the benchmark list if it does not exist
  if (is.null(filtered_pbmc@misc$benchmark)) {
    filtered_pbmc@misc$benchmark <- list()
  }
  
  # Confirm that benchmark storage is a list
  if (!is.list(filtered_pbmc@misc$benchmark)) {
    stop(
      "`pbmc@misc$benchmark` exists but is not a list.",
      call. = FALSE
    )
  }
  
  # 9. Store benchmark results
  filtered_pbmc@misc$benchmark[["filter_cells"]] <-
    benchmark_results
  
  # Store filtering statistics separately
  if (is.null(filtered_pbmc@misc$filtering_summary)) {
    filtered_pbmc@misc$filtering_summary <- list()
  }
  
  filtered_pbmc@misc$filtering_summary[["filter_cells"]] <-
    filtering_summary
  
  message(
    "Cell filtering completed successfully: ",
    cells_before,
    " cells before filtering; ",
    cells_after,
    " cells retained; ",
    cells_removed,
    " cells removed."
  )
  
  if (benchmark) {
    message("\nFiltering summary:")
    print(
      filtering_summary,
      row.names = FALSE
    )
    
    message("\nFiltering benchmark:")
    print(
      benchmark_results,
      row.names = FALSE
    )
  }
  
  return(filtered_pbmc)
}