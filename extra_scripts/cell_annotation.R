library(Seurat)

#Assign cell annotation and save the final RDS object. 


single_cell_annotation <- function(
    pbmc,
    cluster_ids_names,
    reduction_type, 
    pt_size = 0.5,
    title_size = 18,
    legend_size = 18,
    label_size = 4.5,
    guide_legend_size = 10
    ){
  # Assinging cell type to clusters 
  new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                       "NK", "DC", "Platelet")
  names(new.cluster.ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.cluster.ids)
  
  # Plotting new umap with cell type labeling 
  DimPlot(pbmc, reduction = reduction_type, label = TRUE, pt.size = pt_size) + NoLegend()
  
  
  
  plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = label_size) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = title_size), legend.text = element_text(size = legend_size)) + guides(colour = guide_legend(override.aes = list(size = 10)))
  
  # ggsave(filename = '../output/images/pbmc3k_umap.jpg', height = 7, width = 12, plot = plot,
  # quality = 50)
  
  saveRDS(pbmc, file = "../output/pbmc3k_final.rds")
}
