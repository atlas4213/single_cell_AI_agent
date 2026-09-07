#!/usr/bin/env Rscript


# ---------------------------------------------------------------
# 1. Read command-line arguments
# ---------------------------------------------------------------

arguments <- commandArgs(
  trailingOnly = TRUE
)

if (length(arguments) != 3L) {
  stop(
    paste(
      "Three arguments are required:",
      "\n1. Raw 10X path",
      "\n2. Filtered 10X path",
      "\n3. Output directory",
      "\n\nExample:",
      "\nRscript run_pipeline.R",
      "/data/raw_feature_bc_matrix",
      "/data/filtered_feature_bc_matrix",
      "/results"
    ),
    call. = FALSE
  )
}

raw_path <- arguments[[1L]]
filtered_path <- arguments[[2L]]
output_directory <- arguments[[3L]]


# ---------------------------------------------------------------
# 2. Validate input paths
# ---------------------------------------------------------------

if (!dir.exists(raw_path) && !file.exists(raw_path)) {
  stop(
    "Raw 10X input was not found: ",
    raw_path,
    call. = FALSE
  )
}

if (
  !dir.exists(filtered_path) &&
  !file.exists(filtered_path)
) {
  stop(
    "Filtered 10X input was not found: ",
    filtered_path,
    call. = FALSE
  )
}

if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!dir.exists(output_directory)) {
  stop(
    "The output directory could not be created: ",
    output_directory,
    call. = FALSE
  )
}


# ---------------------------------------------------------------
# 3. Load pipeline functions
# ---------------------------------------------------------------

function_files <- list.files(
  path = "./R/pipeline_1/",
  pattern = "\\.R$",
  full.names = TRUE
)

if (length(function_files) == 0L) {
  stop(
    "No R function files were found in /app/R.",
    call. = FALSE
  )
}

invisible(
  lapply(
    function_files,
    source
  )
)


# ---------------------------------------------------------------
# 4. Run pipeline
# ---------------------------------------------------------------

pipeline_start <- Sys.time()

message("Starting single-cell pipeline.")

pbmc <- run_soupx(
  raw_path = raw_path,
  filtered_path = filtered_path,
  project_name = "pbmc",
  output_dir = file.path(
    output_directory,
    "soupx_corrected_counts"
  )
)

pbmc <- qc_data(
  pbmc = pbmc,
  file_name = "pbmc",
  output_dir = file.path(
    output_directory,
    "quality_control"
  )
)

pbmc <- filter_cells(
  pbmc = pbmc,
  nfeature_low = 200,
  nfeature_high = 2500,
  percent_mt_number = 5
)

pbmc <- run_doublet_finder(
  pbmc = pbmc,
  split_by = NULL,
  doublet_rate = 0.075
)

pbmc <- normalize_matrix(
  pbmc = pbmc,
  normalization_method = "LogNormalize",
  scale_factor_number = 10000
)

# If variable_gene_selection() returns a Seurat object
pbmc <- variable_gene_selection(
  pbmc = pbmc,
  output_file = file.path(
    output_directory,
    "variable_genes.pdf"
  )
)

# If scale_matrix_data() still returns a list
scaling_results <- scale_matrix_data(
  pbmc = pbmc
)

pbmc <- scaling_results$pbmc

# If pca_analysis() still returns a list
pca_results <- pca_analysis(
  pbmc = pbmc,
  output_name = "pbmc",
  output_directory = file.path(
    output_directory,
    "pca"
  )
)

pbmc <- pca_results$pbmc


# ---------------------------------------------------------------
# 5. Save final results
# ---------------------------------------------------------------

final_object_path <- file.path(
  output_directory,
  "pbmc_processed.rds"
)

saveRDS(
  object = pbmc,
  file = final_object_path
)

pipeline_seconds <- as.numeric(
  difftime(
    Sys.time(),
    pipeline_start,
    units = "secs"
  )
)

message(
  "Pipeline completed successfully in ",
  round(pipeline_seconds, 2),
  " seconds."
)

message(
  "Final Seurat object saved to: ",
  final_object_path
)