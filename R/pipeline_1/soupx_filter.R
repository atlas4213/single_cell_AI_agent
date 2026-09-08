# Filter Ambient RNA
# Creator: Dave Hill


# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(DropletUtils)
})


#' Remove Ambient RNA Using SoupX
#'
#' Estimates ambient RNA contamination using raw and filtered 10x count
#' matrices, corrects the counts, and returns a Seurat object containing
#' the corrected count matrix.
#'
#' A preliminary Seurat analysis is performed to provide SoupX with cluster
#' assignments and UMAP coordinates.
#'
#' @param raw_matrix A gene-by-droplet raw count matrix containing all
#'   droplets, including empty droplets.
#' @param filt_matrix A gene-by-cell filtered count matrix containing
#'   cell-associated barcodes.
#' @param project_name Character string used as the Seurat project name.
#' @param dims Integer vector specifying the dimensions used for PCA-based
#'   neighbor detection and UMAP. Default is `1:30`.
#' @param resolution Numeric clustering resolution. Default is `0.5`.
#' @param output_dir Optional directory where corrected counts will be saved
#'   in 10x format. If `NULL`, corrected counts are not written to disk.
#' @param overwrite Logical. If `TRUE`, permit `write10xCounts()` to overwrite
#'   an existing output. Default is `FALSE`.
#' @param benchmark Logical. If `TRUE`, print execution times.
#' @param seed Integer random seed used for reproducible dimensionality
#'   reduction and clustering.
#'
#' @return A Seurat object constructed from SoupX-corrected counts.
#'   Benchmark results are stored in
#'   `pbmc@misc$benchmark$run_soupx`.
#'
#' @examples
#' \dontrun{
#' raw_matrix <- Seurat::Read10X(
#'   "./raw_feature_bc_matrix"
#' )
#'
#' filt_matrix <- Seurat::Read10X(
#'   "./filtered_feature_bc_matrix"
#' )
#'
#' pbmc <- run_soupx(
#'   raw_matrix = raw_matrix,
#'   filt_matrix = filt_matrix,
#'   output_dir = "./soupX_corrected_counts"
#' )
#'
#' pbmc@misc$benchmark$run_soupx
#' }
#'
#' @export
run_soupx <- function(
    raw_path,
    filtered_path,
    project_name = "soupX_corrected",
    dims = 1:30,
    resolution = 0.5,
    output_dir = NULL,
    overwrite = FALSE,
    benchmark = TRUE,
    seed = 1234
) {
  
  total_start <- proc.time()[["elapsed"]]
  
  # ---------------------------------------------------------------
  # Internal function for reading a 10X directory or HDF5 file
  # ---------------------------------------------------------------
  
  read_10x_input <- function(path, argument_name) {
    
    if (
      missing(path) ||
      !is.character(path) ||
      length(path) != 1L ||
      is.na(path) ||
      !nzchar(path)
    ) {
      stop(
        "`",
        argument_name,
        "` must be a single, non-empty path.",
        call. = FALSE
      )
    }
    
    if (!dir.exists(path) && !file.exists(path)) {
      stop(
        "The path supplied for `",
        argument_name,
        "` does not exist: ",
        path,
        call. = FALSE
      )
    }
    
    matrix_object <- tryCatch(
      {
        if (dir.exists(path)) {
          
          # Read a directory containing:
          # matrix.mtx, features.tsv and barcodes.tsv
          Seurat::Read10X(
            data.dir = path,
            unique.features = TRUE,
            strip.suffix = FALSE
          )
          
        } else if (
          grepl(
            pattern = "\\.h5$",
            x = path,
            ignore.case = TRUE
          )
        ) {
          
          # Read a 10X HDF5 file
          Seurat::Read10X_h5(
            filename = path,
            use.names = TRUE,
            unique.features = TRUE
          )
          
        } else {
          stop(
            "The path must point to a 10X matrix directory ",
            "or a 10X `.h5` file."
          )
        }
      },
      error = function(e) {
        stop(
          "Failed to read `",
          argument_name,
          "` from: ",
          path,
          "\nOriginal error: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    
    # Read10X can return a list for multimodal datasets
    if (is.list(matrix_object)) {
      if ("Gene Expression" %in% names(matrix_object)) {
        matrix_object <- matrix_object[["Gene Expression"]]
      } else {
        stop(
          "`",
          argument_name,
          "` contains multiple assays, but no ",
          "'Gene Expression' assay was found. Available assays: ",
          paste(names(matrix_object), collapse = ", "),
          call. = FALSE
        )
      }
    }
    
    if (
      !inherits(matrix_object, "Matrix") &&
      !is.matrix(matrix_object)
    ) {
      stop(
        "`",
        argument_name,
        "` did not produce a valid count matrix.",
        call. = FALSE
      )
    }
    
    if (
      nrow(matrix_object) == 0L ||
      ncol(matrix_object) == 0L
    ) {
      stop(
        "`",
        argument_name,
        "` produced an empty count matrix.",
        call. = FALSE
      )
    }
    
    return(matrix_object)
  }
  
  # ---------------------------------------------------------------
  # 1. Read the raw matrix
  # ---------------------------------------------------------------
  
  raw_read_start <- proc.time()[["elapsed"]]
  
  raw_matrix <- read_10x_input(
    path = raw_path,
    argument_name = "raw_path"
  )
  
  raw_read_seconds <-
    proc.time()[["elapsed"]] - raw_read_start
  
  message(
    "Raw matrix loaded: ",
    nrow(raw_matrix),
    " genes x ",
    ncol(raw_matrix),
    " droplets."
  )
  
  # ---------------------------------------------------------------
  # 2. Read the filtered matrix
  # ---------------------------------------------------------------
  
  filtered_read_start <- proc.time()[["elapsed"]]
  
  filt_matrix <- read_10x_input(
    path = filtered_path,
    argument_name = "filtered_path"
  )
  
  filtered_read_seconds <-
    proc.time()[["elapsed"]] - filtered_read_start
  
  message(
    "Filtered matrix loaded: ",
    nrow(filt_matrix),
    " genes x ",
    ncol(filt_matrix),
    " cells."
  )
  
  # ---------------------------------------------------------------
  # 3. Validate the relationship between the matrices
  # ---------------------------------------------------------------
  
  validation_start <- proc.time()[["elapsed"]]
  
  if (is.null(rownames(raw_matrix))) {
    stop(
      "The raw matrix does not contain gene names.",
      call. = FALSE
    )
  }
  
  if (is.null(rownames(filt_matrix))) {
    stop(
      "The filtered matrix does not contain gene names.",
      call. = FALSE
    )
  }
  
  if (is.null(colnames(raw_matrix))) {
    stop(
      "The raw matrix does not contain droplet barcodes.",
      call. = FALSE
    )
  }
  
  if (is.null(colnames(filt_matrix))) {
    stop(
      "The filtered matrix does not contain cell barcodes.",
      call. = FALSE
    )
  }
  
  # Identify genes shared between the two matrices
  common_genes <- intersect(
    rownames(raw_matrix),
    rownames(filt_matrix)
  )
  
  if (length(common_genes) == 0L) {
    stop(
      "The raw and filtered matrices have no genes in common. ",
      "Confirm that they came from the same dataset.",
      call. = FALSE
    )
  }
  
  # Restrict both matrices to shared genes in identical order
  if (
    length(common_genes) != nrow(raw_matrix) ||
    length(common_genes) != nrow(filt_matrix)
  ) {
    warning(
      "The matrices do not contain identical gene sets. ",
      "Both will be restricted to ",
      length(common_genes),
      " shared genes.",
      call. = FALSE
    )
  }
  
  raw_matrix <- raw_matrix[
    common_genes,
    ,
    drop = FALSE
  ]
  
  filt_matrix <- filt_matrix[
    common_genes,
    ,
    drop = FALSE
  ]
  
  # Filtered barcodes should also exist in the raw matrix
  missing_barcodes <- setdiff(
    colnames(filt_matrix),
    colnames(raw_matrix)
  )
  
  if (length(missing_barcodes) > 0L) {
    stop(
      length(missing_barcodes),
      " filtered-cell barcodes are absent from the raw matrix. ",
      "The matrices may not be from the same 10X capture.",
      call. = FALSE
    )
  }
  
  validation_seconds <-
    proc.time()[["elapsed"]] - validation_start
  
  # ---------------------------------------------------------------
  # 4. Create the SoupX channel
  # ---------------------------------------------------------------
  
  channel_start <- proc.time()[["elapsed"]]
  
  soup_channel <- tryCatch(
    {
      SoupX::SoupChannel(
        tod = raw_matrix,
        toc = filt_matrix
      )
    },
    error = function(e) {
      stop(
        "Failed to create the SoupX channel.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  channel_seconds <-
    proc.time()[["elapsed"]] - channel_start
  
  # Continue with the preliminary Seurat processing from your
  # existing run_soupx() function below this point.
  
  # ------------------------------------------------------------------
  # 3. Create the preliminary Seurat object
  # ------------------------------------------------------------------
  
  seurat_start <- proc.time()[["elapsed"]]
  
  srat <- tryCatch(
    {
      Seurat::CreateSeuratObject(
        counts = filt_matrix,
        project = paste0(project_name, "_preliminary"),
        min.cells = 0,
        min.features = 0
      )
    },
    error = function(e) {
      stop(
        "Failed to create the preliminary Seurat object.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  seurat_seconds <-
    proc.time()[["elapsed"]] - seurat_start
  
  # ------------------------------------------------------------------
  # 4. Run preliminary normalization and clustering
  # ------------------------------------------------------------------
  
  clustering_start <- proc.time()[["elapsed"]]
  
  set.seed(seed)
  
  srat <- tryCatch(
    {
      srat <- Seurat::SCTransform(
        object = srat,
        verbose = FALSE
      )
      
      srat <- Seurat::RunPCA(
        object = srat,
        npcs = max(dims),
        verbose = FALSE,
        seed.use = seed
      )
      
      # Check that the requested dimensions were produced
      available_pcs <- ncol(
        Seurat::Embeddings(srat, reduction = "pca")
      )
      
      valid_dims <- dims[dims <= available_pcs]
      
      if (length(valid_dims) == 0L) {
        stop(
          "No requested dimensions are available. ",
          "PCA produced ",
          available_pcs,
          " components."
        )
      }
      
      if (length(valid_dims) < length(dims)) {
        warning(
          "Some requested dimensions exceed the number of available PCs. ",
          "Using dimensions 1 through ",
          max(valid_dims),
          ".",
          call. = FALSE
        )
      }
      
      srat <- Seurat::RunUMAP(
        object = srat,
        dims = valid_dims,
        verbose = FALSE,
        seed.use = seed
      )
      
      srat <- Seurat::FindNeighbors(
        object = srat,
        dims = valid_dims,
        verbose = FALSE
      )
      
      srat <- Seurat::FindClusters(
        object = srat,
        resolution = resolution,
        random.seed = seed,
        verbose = FALSE
      )
      
      srat
    },
    error = function(e) {
      stop(
        "Preliminary Seurat processing or clustering failed.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  clustering_seconds <-
    proc.time()[["elapsed"]] - clustering_start
  
  # ------------------------------------------------------------------
  # 5. Supply clusters and UMAP coordinates to SoupX
  # ------------------------------------------------------------------
  
  metadata_start <- proc.time()[["elapsed"]]
  
  soup_channel <- tryCatch(
    {
      cluster_assignments <- setNames(
        as.character(Seurat::Idents(srat)),
        names(Seurat::Idents(srat))
      )
      
      umap_coordinates <- Seurat::Embeddings(
        srat,
        reduction = "umap"
      )
      
      soup_channel <- SoupX::setClusters(
        soup_channel,
        cluster_assignments
      )
      
      soup_channel <- SoupX::setDR(
        soup_channel,
        umap_coordinates
      )
      
      soup_channel
    },
    error = function(e) {
      stop(
        "Failed to add cluster or UMAP information to SoupX.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  metadata_seconds <-
    proc.time()[["elapsed"]] - metadata_start
  
  # ------------------------------------------------------------------
  # 6. Estimate ambient RNA contamination
  # ------------------------------------------------------------------
  
  estimation_start <- proc.time()[["elapsed"]]
  
  soup_channel <- tryCatch(
    {
      SoupX::autoEstCont(soup_channel)
    },
    error = function(e) {
      stop(
        "SoupX could not estimate the ambient RNA contamination rate.\n",
        "Check the preliminary clusters and verify that the dataset ",
        "contains informative cell-type markers.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  estimation_seconds <-
    proc.time()[["elapsed"]] - estimation_start
  
  # Capture the most abundant genes in the ambient profile
  soup_profile <- soup_channel$soupProfile
  
  top_ambient_genes <- head(
    soup_profile[
      order(
        soup_profile$est,
        decreasing = TRUE
      ),
      ,
      drop = FALSE
    ],
    n = 20
  )
  
  # ------------------------------------------------------------------
  # 7. Correct the count matrix
  # ------------------------------------------------------------------
  
  correction_start <- proc.time()[["elapsed"]]
  
  adjusted_matrix <- tryCatch(
    {
      SoupX::adjustCounts(
        soup_channel,
        roundToInt = TRUE
      )
    },
    error = function(e) {
      stop(
        "SoupX failed to create the corrected count matrix.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  correction_seconds <-
    proc.time()[["elapsed"]] - correction_start
  
  if (
    is.null(adjusted_matrix) ||
    nrow(adjusted_matrix) == 0L ||
    ncol(adjusted_matrix) == 0L
  ) {
    stop(
      "SoupX returned an empty corrected count matrix.",
      call. = FALSE
    )
  }
  
  # ------------------------------------------------------------------
  # 8. Create the corrected Seurat object
  # ------------------------------------------------------------------
  
  corrected_object_start <- proc.time()[["elapsed"]]
  
  corrected_srat <- tryCatch(
    {
      Seurat::CreateSeuratObject(
        counts = adjusted_matrix,
        project = project_name,
        min.cells = 0,
        min.features = 0
      )
    },
    error = function(e) {
      stop(
        "The corrected matrix was created, but the final Seurat object ",
        "could not be created.\n",
        "Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  corrected_object_seconds <-
    proc.time()[["elapsed"]] - corrected_object_start
  
  # ------------------------------------------------------------------
  # 9. Optionally write corrected counts
  # ------------------------------------------------------------------
  
  write_seconds <- 0
  
  if (!is.null(output_dir)) {
    write_start <- proc.time()[["elapsed"]]
    
    tryCatch(
      {
        DropletUtils::write10xCounts(
          path = output_dir,
          x = adjusted_matrix,
          overwrite = overwrite
        )
      },
      error = function(e) {
        stop(
          "Corrected counts were created, but they could not be written to: ",
          output_dir,
          "\nOriginal error: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    
    write_seconds <-
      proc.time()[["elapsed"]] - write_start
  }
  
  # ------------------------------------------------------------------
  # 10. Create and store benchmark results
  # ------------------------------------------------------------------
  
  total_seconds <-
    proc.time()[["elapsed"]] - total_start
  
  benchmark_results <- data.frame(
    step = c(
      "Validate inputs",
      "Create SoupX channel",
      "Create preliminary Seurat object",
      "Normalize and cluster cells",
      "Add clusters and UMAP to SoupX",
      "Estimate ambient RNA contamination",
      "Adjust count matrix",
      "Create corrected Seurat object",
      "Write corrected counts",
      "Total run_soupx function"
    ),
    elapsed_seconds = round(
      c(
        validation_seconds,
        channel_seconds,
        seurat_seconds,
        clustering_seconds,
        metadata_seconds,
        estimation_seconds,
        correction_seconds,
        corrected_object_seconds,
        write_seconds,
        total_seconds
      ),
      digits = 3
    ),
    row.names = NULL
  )
  
  # Initialize benchmark storage
  corrected_srat@misc$benchmark <- list(
    run_soupx = benchmark_results
  )
  
  # Save additional SoupX information
  corrected_srat@misc$soupx <- list(
    top_ambient_genes = top_ambient_genes,
    raw_droplets = ncol(raw_matrix),
    filtered_cells = ncol(filt_matrix),
    corrected_cells = ncol(adjusted_matrix),
    output_directory = output_dir
  )
  
  message(
    "SoupX correction completed successfully: ",
    ncol(raw_matrix),
    " raw droplets, ",
    ncol(filt_matrix),
    " filtered cells, and ",
    ncol(adjusted_matrix),
    " corrected cells."
  )
  
  if (benchmark) {
    message("\nTop estimated ambient genes:")
    print(top_ambient_genes)
    
    message("\nSoupX benchmark:")
    print(
      benchmark_results,
      row.names = FALSE
    )
  }
  
  return(corrected_srat)
}