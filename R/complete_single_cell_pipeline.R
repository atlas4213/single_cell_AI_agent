# Single cell pipeine script

source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/soupx_filter.R")
#source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/validation_script.R")
source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/quality_control.R")
source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/cell_filtering.R")
source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/run_doublet.R")
source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/normalize_matrix.R")
source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/variable_selection.R")
source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/scaling_data.R")
source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/pca_analysis.R")
source("/Users/Atlas/Desktop/single_cell_practice_pipeline/R/variable_genes_heatmap.R")

#directory_path <- "/Users/Atlas/Desktop/single_cell_practice_pipeline/filtered_gene_bc_matrices/hg19/"
directory_path_raw <- "/Users/Atlas/Desktop/single_cell_practice_pipeline/data/raw_gene_bc_matrices/hg19/"
directory_path_filtered <- "/Users/Atlas/Desktop/single_cell_practice_pipeline/data/filtered_gene_bc_matrices/hg19/"





pbmc_soup <- run_soupx(directory_path_raw, directory_path_filtered)
#pbmc_raw <- load_data(directory_path)
pbmc_qc <- qc_data(pbmc_soup, "single_cell")
pbmc_filtered <- filter_cells(
  pbmc_qc,
  nfeature_high = 2500,
  nfeature_low = 200,
  percent_mt_number = 5
)
pbmc_doublet <- run_doublet_finder(pbmc_filtered)
pbmc_norm <- normalize_matrix(pbmc_doublet)
pbmc_norm_var_sel <- variable_gene_selection(pbmc_norm)
pbmc_scaled <- scale_matrix_data(pbmc = pbmc_norm_var_sel)
pbmc_pca <- pca_analysis(pbmc_scaled$pbmc)
pbmc_heatmap <-variable_genes_heatmap(pbmc_pca$pbmc)
