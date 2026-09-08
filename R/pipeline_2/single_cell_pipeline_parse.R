setwd("~/Desktop/single_cell_practice_pipeline")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#Load in the pmbc dataset 
pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")

#Create the seurat object with the raw (non-normalized data)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", 
                           min.cells = 3, min.features = 200)
pbmc

#An object of class Seurat 
#13714 features across 2700 samples within 1 assay 
#Active assay: RNA (13714 features, 0 variable features)


#Quality Control 
#Add mitochondiral percentages to the metadata 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#Visualize the QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Create and store the plot
qc_violin_plot <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Display the plot
qc_violin_plot


ggsave(
  filename = "./single_cell_qc_violin_plot.pdf",
  plot = qc_violin_plot,
  width = 12,
  height = 5
)


#Visualize feature to feature relationships 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Save the feature plots
ggsave(
  filename = "./single_cell_qc_feature_plots.pdf",
  plot = plot1 + plot2,
  width = 12,
  height = 5
)

#Subset the data for the proper QC metrics 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# Normalize the data 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify Highly Variable Genes 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the top 10 most highly variable genes 
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable featues with and without labels 
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2 

#Save the Highly Variable Genes plots
ggsave(
  filename = "./single_cell_variable_genes_plots.pdf",
  plot = plot1 + plot2,
  width = 12,
  height = 5
)

# Scale the data (to normalize expression amongst genes)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Perform linear reduction (PCA) and get principal componets from HVG's 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Create a plot showing the genes that most heavily effected the PCA 
VizDimLoadings(pbmc, dims=1:2, reduction="pca")

#Save the PCA plots
ggsave(
  filename = "./single_cell_PC_compare_plot.pdf",
  plot = plot1 + plot2,
  width = 12,
  height = 5
)

#Create graph of the way the cells are effected by the 2 PC's 
DimPlot(pbmc, reduction = "pca") + NoLegend()

#Save the PCA plots
ggsave(
  filename = "./single_cell_pca1_vs_pca2_plots.pdf",
  plot = plot1 + plot2,
  width = 12,
  height = 5
)

# Creates a heatmap to further explore the dataset
DimHeatmap(pbmc, dims =1, cells = 500, balanced =TRUE)

#Save the Highly Variable Genes plots
ggsave(
  filename = "./single_cell_heatmap_plot.pdf",
  plot = plot1 + plot2,
  width = 12,
  height = 5
)

# Rank the principal components  for a cutoff 
ElbowPlot(pbmc)

#Save the Highly Variable Genes plots
ggsave(
  filename = "./single_cell_elbow_plot.pdf",
  plot = plot1 + plot2,
  width = 12,
  height = 5
)

# Cluster the cells 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


# Run UMAP 
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Plot UMAP
DimPlot(pbmc, reduction = "umap")


#Save single cell object
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

# Find differential expressed features 
# Find all markers of cluster 2 
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n =5)


# Find all markers distinguishing cluster 5 from cluster 0 and 3 
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0,3))
head(cluster5.markers)


# Find markers for every cluster compared to all remaining cells 
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


# Visualizing the marker expressions 
VlnPlot(pbmc, features = c("MS4A1","CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), layer = "counts", log = TRUE)

# Shows where in the umap certain genes were expressed 
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# Heatmap to show gene expression across all cells and features 
pbmc.markers %>% 
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


# Cannonical markers for identifying Genes 
# Cluster ID	Markers	Cell Type
# 0	IL7R, CCR7	Naive CD4+ T
# 1	CD14, LYZ	CD14+ Mono
# 2	IL7R, S100A4	Memory CD4+
# 3	MS4A1	B
# 4	CD8A	CD8+ T
# 5	FCGR3A, MS4A7	FCGR3A+ Mono
# 6	GNLY, NKG7	NK
# 7	FCER1A, CST3	DC
# 8	PPBP	Platelet

# Asssinging cell type to clusters 
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Plotting new umap with cell type labeling 
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))

# ggsave(filename = '../output/images/pbmc3k_umap.jpg', height = 7, width = 12, plot = plot,
# quality = 50)

saveRDS(pbmc, file = "../output/pbmc3k_final.rds")

