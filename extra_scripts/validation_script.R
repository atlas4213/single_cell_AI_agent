### Validation Script
### Creator: Dave Hill

library(Seurat)


#' Load 10x Genomics Data into a Seurat Object
#'
#' Reads a filtered 10x Genomics feature-barcode matrix, validates the
#' data, creates a Seurat object, and records execution times.
#'
#' @param file A character string containing the path to a 10x Genomics
#'   data directory.
#' @param benchmark Logical. If `TRUE`, print and store execution times.
#'   Default is `TRUE`.
#'
#' @return A Seurat object. Benchmark results are stored in
#'   `pbmc@misc$benchmark`.
#'
#' @examples
#' \dontrun{
#' pbmc <- load_data("./filtered_gene_bc_matrices/hg19/")
#'
#' # View benchmark results
#' pbmc@misc$benchmark
#' }
#'
#' @seealso
#' [Seurat::Read10X()]
#' [Seurat::CreateSeuratObject()]
#'
#' @export
#' @importFrom Seurat Read10X CreateSeuratObject
load_data <- function(file, benchmark = TRUE) {
  
  # Start total execution timer
  total_start <- proc.time()[["elapsed"]]
  
  # 1. Validate the input path
  if (
    missing(file) ||
    is.null(file) ||
    length(file) != 1L ||
    !is.character(file) ||
    is.na(file) ||
    !nzchar(file)
  ) {
    stop(
      "A single, non-empty 10x data directory must be provided.",
      call. = FALSE
    )
  }
  
  if (!dir.exists(file)) {
    stop(
      "The 10x data directory does not exist: ",
      file,
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
  
  # 2. Read the 10X count matrix
  read_start <- proc.time()[["elapsed"]]
  
  pbmc.data <- tryCatch(
    Seurat::Read10X(data.dir = file),
    error = function(e) {
      stop(
        "Failed to read 10X data from: ",
        file,
        "\nOriginal error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  read_seconds <- proc.time()[["elapsed"]] - read_start
  
  # 3. Check that data was returned
  if (is.null(pbmc.data) || length(pbmc.data) == 0) {
    stop(
      "Read10X returned no data for directory: ",
      file,
      call. = FALSE
    )
  }
  
  # Read10X can return a list when multiple assays are present
  if (is.list(pbmc.data)) {
    if ("Gene Expression" %in% names(pbmc.data)) {
      pbmc.data <- pbmc.data[["Gene Expression"]]
    } else {
      stop(
        "Multiple assays were found, but no 'Gene Expression' assay exists. ",
        "Available assays: ",
        paste(names(pbmc.data), collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  # 4. Validate the count matrix
  if (nrow(pbmc.data) == 0 || ncol(pbmc.data) == 0) {
    stop(
      "The count matrix is empty. Dimensions: ",
      nrow(pbmc.data),
      " genes x ",
      ncol(pbmc.data),
      " cells.",
      call. = FALSE
    )
  }
  
  # 5. Create the Seurat object
  creation_start <- proc.time()[["elapsed"]]
  
  pbmc <- tryCatch(
    Seurat::CreateSeuratObject(
      counts = pbmc.data,
      project = "pbmc3k",
      min.cells = 3,
      min.features = 200
    ),
    error = function(e) {
      stop(
        "The count matrix was read successfully, but the Seurat object ",
        "could not be created.\nOriginal error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  creation_seconds <- proc.time()[["elapsed"]] - creation_start
  
  # 6. Final validation
  if (ncol(pbmc) == 0) {
    stop(
      "The Seurat object contains zero cells after filtering. ",
      "Consider lowering `min.features`.",
      call. = FALSE
    )
  }
  
  total_seconds <- proc.time()[["elapsed"]] - total_start
  
  # Create a benchmark table
  benchmark_results <- data.frame(
    step = c(
      "Read 10X matrix",
      "Create Seurat object",
      "Total load_data function"
    ),
    elapsed_seconds = round(
      c(
        read_seconds,
        creation_seconds,
        total_seconds
      ),
      digits = 3
    ),
    row.names = NULL
  )
  
  # Store the benchmark results in the Seurat object
  pbmc@misc$benchmark <- list( 
    benchmark_results
  )
  
  message(
    "Seurat object created successfully: ",
    nrow(pbmc),
    " genes x ",
    ncol(pbmc),
    " cells."
  )
  
  if (benchmark) {
    message("\nBenchmark results:")
    print(benchmark_results, row.names = FALSE)
  }
  
  return(pbmc)
}




