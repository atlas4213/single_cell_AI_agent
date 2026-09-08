library(Seurat)

# Plot a heatmap of the top marker genes for each cluster

gene_expression_heatmap <- function(
    pbmc,
    pbmc_markers,
    filter_threshold = 1,
    n = 10,
    assay_name = NULL,
    slot_type = "scale.data",
    group_by = NULL
) {
  
  # Start total function benchmark
  total_start_time <- Sys.time()
  
  # -----------------------------------------------
  # 1. Validate the Seurat object
  # -----------------------------------------------
  
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
  
  if (nrow(pbmc) == 0) {
    stop("The Seurat object does not contain any features.")
  }
  
  # -----------------------------------------------
  # 2. Validate the marker table
  # -----------------------------------------------
  
  if (
    missing(pbmc_markers) ||
    is.null(pbmc_markers)
  ) {
    stop(
      "A marker-gene table must be supplied to 'pbmc_markers'. ",
      "Use the output from FindAllMarkers()."
    )
  }
  
  if (!is.data.frame(pbmc_markers)) {
    stop(
      "'pbmc_markers' must be a data frame. Received: ",
      paste(class(pbmc_markers), collapse = ", ")
    )
  }
  
  if (nrow(pbmc_markers) == 0) {
    stop("The marker table does not contain any rows.")
  }
  
  if (!"cluster" %in% colnames(pbmc_markers)) {
    stop(
      "The marker table must contain a column named 'cluster'."
    )
  }
  
  # Detect the gene-name column
  if ("gene" %in% colnames(pbmc_markers)) {
    gene_column <- "gene"
  } else if (!is.null(rownames(pbmc_markers))) {
    pbmc_markers$gene <- rownames(pbmc_markers)
    gene_column <- "gene"
    
    warning(
      "No 'gene' column was found. Row names will be used ",
      "as the gene names."
    )
  } else {
    stop(
      "The marker table must contain a 'gene' column ",
      "or gene names as row names."
    )
  }
  
  # Detect the fold-change column
  if ("avg_log2FC" %in% colnames(pbmc_markers)) {
    fold_change_column <- "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(pbmc_markers)) {
    fold_change_column <- "avg_logFC"
    
    warning(
      "The marker table contains 'avg_logFC' instead of ",
      "'avg_log2FC'. Check the Seurat version and logarithm base."
    )
  } else {
    stop(
      "The marker table must contain an 'avg_log2FC' ",
      "or 'avg_logFC' column."
    )
  }
  
  if (!is.numeric(pbmc_markers[[fold_change_column]])) {
    stop(
      "The fold-change column '",
      fold_change_column,
      "' must be numeric."
    )
  }
  
  if (anyNA(pbmc_markers$cluster)) {
    stop("The marker table contains missing cluster values.")
  }
  
  if (
    anyNA(pbmc_markers[[gene_column]]) ||
    any(trimws(as.character(
      pbmc_markers[[gene_column]]
    )) == "")
  ) {
    stop("The marker table contains missing or empty gene names.")
  }
  
  # -----------------------------------------------
  # 3. Validate filtering arguments
  # -----------------------------------------------
  
  if (
    !is.numeric(filter_threshold) ||
    length(filter_threshold) != 1 ||
    is.na(filter_threshold) ||
    !is.finite(filter_threshold)
  ) {
    stop(
      "'filter_threshold' must be one finite numeric value."
    )
  }
  
  if (
    !is.numeric(n) ||
    length(n) != 1 ||
    is.na(n) ||
    !is.finite(n) ||
    n <= 0 ||
    n != as.integer(n)
  ) {
    stop("'n' must be one positive whole number.")
  }
  
  # -----------------------------------------------
  # 4. Validate assay and slot
  # -----------------------------------------------
  
  available_assays <- Assays(pbmc)
  
  if (length(available_assays) == 0) {
    stop("The Seurat object does not contain any assays.")
  }
  
  if (is.null(assay_name)) {
    assay_name <- DefaultAssay(pbmc)
  }
  
  if (
    !is.character(assay_name) ||
    length(assay_name) != 1 ||
    is.na(assay_name) ||
    trimws(assay_name) == ""
  ) {
    stop(
      "'assay_name' must be NULL or one non-empty character value."
    )
  }
  
  if (!assay_name %in% available_assays) {
    stop(
      "The assay '",
      assay_name,
      "' was not found. Available assays: ",
      paste(available_assays, collapse = ", "),
      "."
    )
  }
  
  valid_slots <- c("counts", "data", "scale.data")
  
  if (
    !is.character(slot_type) ||
    length(slot_type) != 1 ||
    is.na(slot_type) ||
    !slot_type %in% valid_slots
  ) {
    stop(
      "'slot_type' must be one of: ",
      paste(valid_slots, collapse = ", "),
      "."
    )
  }
  
  # -----------------------------------------------
  # 5. Validate optional grouping column
  # -----------------------------------------------
  
  if (!is.null(group_by)) {
    if (
      !is.character(group_by) ||
      length(group_by) != 1 ||
      is.na(group_by) ||
      trimws(group_by) == ""
    ) {
      stop(
        "'group_by' must be NULL or one non-empty character value."
      )
    }
    
    if (!group_by %in% colnames(pbmc[[]])) {
      stop(
        "The grouping column '",
        group_by,
        "' was not found in the cell metadata."
      )
    }
  } else {
    if (length(unique(Idents(pbmc))) < 2) {
      warning(
        "The Seurat object contains fewer than two active identities. ",
        "The heatmap may not show meaningful cluster groupings."
      )
    }
  }
  
  # -----------------------------------------------
  # 6. Filter and benchmark marker selection
  # -----------------------------------------------
  
  message(
    "Selecting up to ",
    n,
    " marker genes per cluster with ",
    fold_change_column,
    " > ",
    filter_threshold,
    "..."
  )
  
  selection_start_time <- Sys.time()
  
  top_markers <- tryCatch(
    {
      filtered_markers <- pbmc_markers[
        !is.na(pbmc_markers[[fold_change_column]]) &
          pbmc_markers[[fold_change_column]] >
          filter_threshold,
        ,
        drop = FALSE
      ]
      
      if (nrow(filtered_markers) == 0) {
        stop(
          "No markers passed the fold-change threshold of ",
          filter_threshold,
          "."
        )
      }
      
      # Split markers by cluster
      markers_by_cluster <- split(
        filtered_markers,
        filtered_markers$cluster
      )
      
      # Select the top n markers in each cluster
      top_marker_list <- lapply(
        markers_by_cluster,
        function(cluster_markers) {
          
          cluster_markers <- cluster_markers[
            order(
              cluster_markers[[fold_change_column]],
              decreasing = TRUE,
              na.last = NA
            ),
            ,
            drop = FALSE
          ]
          
          head(cluster_markers, n = as.integer(n))
        }
      )
      
      do.call(
        rbind,
        top_marker_list
      )
    },
    error = function(e) {
      elapsed_before_error <- as.numeric(
        difftime(
          Sys.time(),
          selection_start_time,
          units = "secs"
        )
      )
      
      stop(
        "Top-marker selection failed after ",
        round(elapsed_before_error, 3),
        " seconds: ",
        conditionMessage(e)
      )
    }
  )
  
  selection_end_time <- Sys.time()
  
  selection_elapsed_seconds <- as.numeric(
    difftime(
      selection_end_time,
      selection_start_time,
      units = "secs"
    )
  )
  
  # Remove duplicated genes to avoid repeated heatmap rows
  selected_features <- unique(
    as.character(top_markers[[gene_column]])
  )
  
  if (length(selected_features) == 0) {
    stop("No genes remained after top-marker selection.")
  }
  
  # -----------------------------------------------
  # 7. Confirm genes exist in the assay
  # -----------------------------------------------
  
  available_features <- rownames(pbmc[[assay_name]])
  
  missing_features <- setdiff(
    selected_features,
    available_features
  )
  
  if (length(missing_features) > 0) {
    warning(
      length(missing_features),
      " selected marker gene(s) were not found in assay '",
      assay_name,
      "' and will be excluded: ",
      paste(missing_features, collapse = ", ")
    )
    
    selected_features <- setdiff(
      selected_features,
      missing_features
    )
  }
  
  if (length(selected_features) == 0) {
    stop(
      "None of the selected marker genes were found in assay '",
      assay_name,
      "'."
    )
  }
  
  # -----------------------------------------------
  # 8. Confirm genes exist in the selected slot
  # -----------------------------------------------
  
  slot_features <- tryCatch(
    {
      rownames(
        GetAssayData(
          object = pbmc,
          assay = assay_name,
          slot = slot_type
        )
      )
    },
    error = function(e) {
      stop(
        "Unable to retrieve the '",
        slot_type,
        "' slot from assay '",
        assay_name,
        "': ",
        conditionMessage(e)
      )
    }
  )
  
  missing_from_slot <- setdiff(
    selected_features,
    slot_features
  )
  
  if (length(missing_from_slot) > 0) {
    warning(
      length(missing_from_slot),
      " selected marker gene(s) were not found in the '",
      slot_type,
      "' slot and will be excluded: ",
      paste(missing_from_slot, collapse = ", ")
    )
    
    selected_features <- setdiff(
      selected_features,
      missing_from_slot
    )
  }
  
  if (length(selected_features) == 0) {
    stop(
      "None of the selected marker genes are available in the '",
      slot_type,
      "' slot. If using 'scale.data', make sure ScaleData() ",
      "was run for these genes."
    )
  }
  
  # -----------------------------------------------
  # 9. Create and benchmark the heatmap
  # -----------------------------------------------
  
  message(
    "Creating heatmap with ",
    length(selected_features),
    " unique marker genes..."
  )
  
  plot_start_time <- Sys.time()
  
  heatmap_plot <- tryCatch(
    {
      DoHeatmap(
        object = pbmc,
        features = selected_features,
        assay = assay_name,
        slot = slot_type,
        group.by = group_by
      ) + NoLegend()
    },
    error = function(e) {
      elapsed_before_error <- as.numeric(
        difftime(
          Sys.time(),
          plot_start_time,
          units = "secs"
        )
      )
      
      stop(
        "Heatmap creation failed after ",
        round(elapsed_before_error, 3),
        " seconds: ",
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
  
  print(heatmap_plot)
  
  # -----------------------------------------------
  # 10. Calculate total runtime
  # -----------------------------------------------
  
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
      "Top-marker selection",
      "Heatmap creation",
      "Entire function"
    ),
    elapsed_seconds = round(
      c(
        selection_elapsed_seconds,
        plot_elapsed_seconds,
        total_elapsed_seconds
      ),
      3
    )
  )
  
  # -----------------------------------------------
  # 11. Store parameters and benchmark
  # -----------------------------------------------
  
  pbmc@misc$pipeline_parameters$gene_expression_heatmap <- list(
    filter_threshold = filter_threshold,
    top_genes_per_cluster = n,
    assay = assay_name,
    slot = slot_type,
    group_by = group_by,
    fold_change_column = fold_change_column,
    selected_features = selected_features
  )
  
  pbmc@misc$benchmarks$gene_expression_heatmap <- benchmark
  
  # -----------------------------------------------
  # 12. Report results
  # -----------------------------------------------
  
  message("Gene-expression heatmap created successfully.")
  message("Clusters represented: ", length(unique(top_markers$cluster)))
  message("Unique genes plotted: ", length(selected_features))
  
  message(
    "Marker-selection runtime: ",
    round(selection_elapsed_seconds, 2),
    " seconds."
  )
  
  message(
    "Heatmap runtime: ",
    round(plot_elapsed_seconds, 2),
    " seconds."
  )
  
  message(
    "Total function runtime: ",
    round(total_elapsed_seconds, 2),
    " seconds."
  )
  
  print(benchmark)
  
  # -----------------------------------------------
  # 13. Return results
  # -----------------------------------------------
  
  return(
    list(
      pbmc = pbmc,
      top_markers = top_markers,
      selected_features = selected_features,
      heatmap_plot = heatmap_plot,
      benchmark = benchmark
    )
  )
}

# -------------------------------------------------
# Example usage
# -------------------------------------------------

# Use the marker table returned by FindAllMarkers() or by find_markers().
# This performs the original requested workflow:
#   1. filter markers above the fold-change threshold,
#   2. select the top n genes from each cluster, and
#   3. plot their scaled expression with DoHeatmap().

# heatmap_results <- gene_expression_heatmap(
#   pbmc = pbmc,
#   pbmc_markers = pbmc.markers,
#   filter_threshold = 1,
#   n = 10,
#   assay_name = "RNA",
#   slot_type = "scale.data",
#   group_by = "seurat_clusters"
# )

# Preserve the updated Seurat object and access the outputs:
# pbmc <- heatmap_results$pbmc
# heatmap_results$heatmap_plot
# heatmap_results$top_markers
# heatmap_results$benchmark





# library(Seurat)
# 
# # Find marker genes for all clusters and visualize selected genes
# 
# find_markers <- function(
#     pbmc,
#     log2fc_filter_value = 1,
#     feature_names = c("MS4A1", "CD79A"),
#     only_positive = TRUE,
#     min_pct = 0.1,
#     test_method = "wilcox"
# ) {
#   
#   # Start total function benchmark
#   total_start_time <- Sys.time()
#   
#   # -----------------------------------------------
#   # 1. Validate the Seurat object
#   # -----------------------------------------------
#   
#   if (missing(pbmc) || is.null(pbmc)) {
#     stop("A Seurat object must be supplied to 'pbmc'.")
#   }
#   
#   if (!inherits(pbmc, "Seurat")) {
#     stop(
#       "'pbmc' must be a Seurat object. Received: ",
#       paste(class(pbmc), collapse = ", ")
#     )
#   }
#   
#   if (ncol(pbmc) == 0) {
#     stop("The Seurat object does not contain any cells.")
#   }
#   
#   if (nrow(pbmc) == 0) {
#     stop("The Seurat object does not contain any features.")
#   }
#   
#   # -----------------------------------------------
#   # 2. Validate cluster identities
#   # -----------------------------------------------
#   
#   cell_identities <- Idents(pbmc)
#   
#   if (length(cell_identities) == 0) {
#     stop("No active cell identities were found.")
#   }
#   
#   if (anyNA(cell_identities)) {
#     stop(
#       "The active cell identities contain missing values. ",
#       "Assign every cell to a cluster before finding markers."
#     )
#   }
#   
#   number_of_clusters <- length(unique(cell_identities))
#   
#   if (number_of_clusters < 2) {
#     stop(
#       "At least two clusters are required for FindAllMarkers(). ",
#       "The object currently contains ",
#       number_of_clusters,
#       " cluster."
#     )
#   }
#   
#   # -----------------------------------------------
#   # 3. Validate the log2 fold-change threshold
#   # -----------------------------------------------
#   
#   if (
#     !is.numeric(log2fc_filter_value) ||
#     length(log2fc_filter_value) != 1 ||
#     is.na(log2fc_filter_value) ||
#     !is.finite(log2fc_filter_value)
#   ) {
#     stop(
#       "'log2fc_filter_value' must be one finite numeric value."
#     )
#   }
#   
#   if (only_positive && log2fc_filter_value < 0) {
#     warning(
#       "'only_positive' is TRUE, but 'log2fc_filter_value' is negative. ",
#       "This will retain nearly all positive markers."
#     )
#   }
#   
#   # -----------------------------------------------
#   # 4. Validate feature names
#   # -----------------------------------------------
#   
#   if (
#     !is.character(feature_names) ||
#     length(feature_names) == 0 ||
#     anyNA(feature_names) ||
#     any(trimws(feature_names) == "")
#   ) {
#     stop(
#       "'feature_names' must be a character vector containing ",
#       "at least one valid feature name."
#     )
#   }
#   
#   feature_names <- unique(feature_names)
#   
#   available_features <- rownames(pbmc)
#   
#   missing_features <- setdiff(
#     feature_names,
#     available_features
#   )
#   
#   if (length(missing_features) > 0) {
#     stop(
#       "The following requested features were not found in the ",
#       "Seurat object: ",
#       paste(missing_features, collapse = ", "),
#       "."
#     )
#   }
#   
#   # -----------------------------------------------
#   # 5. Validate remaining marker arguments
#   # -----------------------------------------------
#   
#   if (
#     !is.logical(only_positive) ||
#     length(only_positive) != 1 ||
#     is.na(only_positive)
#   ) {
#     stop("'only_positive' must be either TRUE or FALSE.")
#   }
#   
#   if (
#     !is.numeric(min_pct) ||
#     length(min_pct) != 1 ||
#     is.na(min_pct) ||
#     !is.finite(min_pct) ||
#     min_pct < 0 ||
#     min_pct > 1
#   ) {
#     stop("'min_pct' must be one numeric value between 0 and 1.")
#   }
#   
#   if (
#     !is.character(test_method) ||
#     length(test_method) != 1 ||
#     is.na(test_method) ||
#     trimws(test_method) == ""
#   ) {
#     stop("'test_method' must be one non-empty character value.")
#   }
#   
#   # -----------------------------------------------
#   # 6. Find and benchmark markers
#   # -----------------------------------------------
#   
#   message(
#     "Finding marker genes for ",
#     number_of_clusters,
#     " clusters..."
#   )
#   
#   marker_start_time <- Sys.time()
#   
#   all_markers <- tryCatch(
#     {
#       FindAllMarkers(
#         object = pbmc,
#         only.pos = only_positive,
#         min.pct = min_pct,
#         test.use = test_method,
#         verbose = FALSE
#       )
#     },
#     error = function(e) {
#       elapsed_before_error <- as.numeric(
#         difftime(
#           Sys.time(),
#           marker_start_time,
#           units = "secs"
#         )
#       )
#       
#       stop(
#         "Marker detection failed after ",
#         round(elapsed_before_error, 3),
#         " seconds: ",
#         conditionMessage(e)
#       )
#     }
#   )
#   
#   marker_end_time <- Sys.time()
#   
#   marker_elapsed_seconds <- as.numeric(
#     difftime(
#       marker_end_time,
#       marker_start_time,
#       units = "secs"
#     )
#   )
#   
#   # -----------------------------------------------
#   # 7. Validate marker results
#   # -----------------------------------------------
#   
#   if (is.null(all_markers) || nrow(all_markers) == 0) {
#     stop(
#       "FindAllMarkers() completed, but no marker genes were found. ",
#       "Consider lowering min_pct or checking the cluster identities."
#     )
#   }
#   
#   if (!"cluster" %in% colnames(all_markers)) {
#     stop(
#       "The marker results do not contain the expected 'cluster' column."
#     )
#   }
#   
#   # Current Seurat versions generally use avg_log2FC.
#   # Older versions may return avg_logFC.
#   if ("avg_log2FC" %in% colnames(all_markers)) {
#     fold_change_column <- "avg_log2FC"
#   } else if ("avg_logFC" %in% colnames(all_markers)) {
#     fold_change_column <- "avg_logFC"
#     
#     warning(
#       "The marker results contain 'avg_logFC' instead of ",
#       "'avg_log2FC'. Check the Seurat version and logarithm base."
#     )
#   } else {
#     stop(
#       "The marker results do not contain 'avg_log2FC' ",
#       "or 'avg_logFC'."
#     )
#   }
#   
#   # -----------------------------------------------
#   # 8. Filter and benchmark marker results
#   # -----------------------------------------------
#   
#   filter_start_time <- Sys.time()
#   
#   filtered_markers <- all_markers[
#     !is.na(all_markers[[fold_change_column]]) &
#       all_markers[[fold_change_column]] >
#       log2fc_filter_value,
#     ,
#     drop = FALSE
#   ]
#   
#   # Sort results by cluster and decreasing fold change
#   if (nrow(filtered_markers) > 0) {
#     filtered_markers <- filtered_markers[
#       order(
#         filtered_markers$cluster,
#         -filtered_markers[[fold_change_column]]
#       ),
#       ,
#       drop = FALSE
#     ]
#   }
#   
#   filter_end_time <- Sys.time()
#   
#   filter_elapsed_seconds <- as.numeric(
#     difftime(
#       filter_end_time,
#       filter_start_time,
#       units = "secs"
#     )
#   )
#   
#   if (nrow(filtered_markers) == 0) {
#     warning(
#       "Markers were found, but none passed the fold-change cutoff of ",
#       log2fc_filter_value,
#       ". Consider using a lower cutoff."
#     )
#   }
#   
#   # -----------------------------------------------
#   # 9. Create and benchmark violin plots
#   # -----------------------------------------------
#   
#   plot_start_time <- Sys.time()
#   
#   violin_plot <- tryCatch(
#     {
#       VlnPlot(
#         object = pbmc,
#         features = feature_names
#       )
#     },
#     error = function(e) {
#       stop(
#         "Marker detection completed, but violin-plot creation failed: ",
#         conditionMessage(e)
#       )
#     }
#   )
#   
#   plot_end_time <- Sys.time()
#   
#   plot_elapsed_seconds <- as.numeric(
#     difftime(
#       plot_end_time,
#       plot_start_time,
#       units = "secs"
#     )
#   )
#   
#   print(violin_plot)
#   
#   # -----------------------------------------------
#   # 10. Calculate total runtime
#   # -----------------------------------------------
#   
#   total_end_time <- Sys.time()
#   
#   total_elapsed_seconds <- as.numeric(
#     difftime(
#       total_end_time,
#       total_start_time,
#       units = "secs"
#     )
#   )
#   
#   benchmark <- data.frame(
#     step = c(
#       "FindAllMarkers",
#       "Filter markers",
#       "Violin plot",
#       "Entire function"
#     ),
#     elapsed_seconds = round(
#       c(
#         marker_elapsed_seconds,
#         filter_elapsed_seconds,
#         plot_elapsed_seconds,
#         total_elapsed_seconds
#       ),
#       3
#     )
#   )
#   
#   # -----------------------------------------------
#   # 11. Store parameters and benchmark
#   # -----------------------------------------------
#   
#   pbmc@misc$pipeline_parameters$markers <- list(
#     log2fc_filter_value = log2fc_filter_value,
#     fold_change_column = fold_change_column,
#     feature_names = feature_names,
#     only_positive = only_positive,
#     min_pct = min_pct,
#     test_method = test_method,
#     number_of_clusters = number_of_clusters,
#     markers_found = nrow(all_markers),
#     markers_passing_filter = nrow(filtered_markers)
#   )
#   
#   pbmc@misc$benchmarks$find_markers <- benchmark
#   
#   # -----------------------------------------------
#   # 12. Report results
#   # -----------------------------------------------
#   
#   message("Marker detection completed successfully.")
#   message("All markers found: ", nrow(all_markers))
#   message(
#     "Markers passing the fold-change cutoff: ",
#     nrow(filtered_markers)
#   )
#   
#   message(
#     "FindAllMarkers runtime: ",
#     round(marker_elapsed_seconds, 2),
#     " seconds."
#   )
#   
#   message(
#     "Total function runtime: ",
#     round(total_elapsed_seconds, 2),
#     " seconds."
#   )
#   
#   print(benchmark)
#   
#   # -----------------------------------------------
#   # 13. Return results
#   # -----------------------------------------------
#   
#   return(
#     list(
#       pbmc = pbmc,
#       all_markers = all_markers,
#       filtered_markers = filtered_markers,
#       violin_plot = violin_plot,
#       benchmark = benchmark
#     )
#   )
# }