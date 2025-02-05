# Remove clusters with high mt
# cluster each sample, remove clusters with median mt > 10%

library(Seurat)
library(dplyr)
library(ggpubr)

wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/01-Regional_data_processing/VST/"
setwd(wd)

parent_dir <- "/banach2/SCORCH/analysis/240610_pipeline/vst/"
samples <- list.dirs(parent_dir,recursive = F,full.names = F)
for (sample in samples) {
  seurat.obj <- readRDS(paste0(parent_dir,sample,"/label_transfer/labeled_",sample,".rds"))
  metadata <- read.csv(paste0(parent_dir,sample,"/label_transfer/labeled_metadata_",sample,".csv"),row.names = 1)
  seurat.obj@meta.data <- metadata
  
  # note: PCA computed from SCT
  DefaultAssay(seurat.obj) <- "SCT"
  seurat.obj <- FindNeighbors(seurat.obj)
  seurat.obj <- FindClusters(seurat.obj,resolution = 3)
  
  p1 <- FeaturePlot(seurat.obj,features = "percent.mt")
  p2 <- DimPlot(seurat.obj,reduction = "umap",group.by = "seurat_clusters",label = T)+ NoLegend()
  p3 <- VlnPlot(seurat.obj,features = "percent.mt",group.by = "seurat_clusters")+ NoLegend()
  
  median_mt <- seurat.obj@meta.data %>%
    select(percent.mt,seurat_clusters) %>%
    group_by(seurat_clusters) %>%
    summarise(median_mt = median(percent.mt)) %>%
    mutate(cluster2remove = median_mt > 10) %>%
    mutate(cluster_idx = as.character(seurat_clusters))
  
  seurat.obj@meta.data$cell2remove_highMT <- seurat.obj@meta.data$seurat_clusters %in% 
    median_mt[median_mt$cluster2remove,]$cluster_idx
  
  message(paste0("Sample: ",sample,",",sum(seurat.obj@meta.data$cell2remove_highMT)," cells removed in ",sum(median_mt$cluster2remove)," clusters"))
  
  # output
  write.csv(seurat.obj@meta.data,paste0("./data/sample_metadata_withMT/",sample,".csv"))
  
  png(filename = paste0("./figures/sample_MT_check/",sample,".png"),width = 1500,height = 800,res = 100)
  print(ggarrange(p1,p2,p3,nrow = 1,ncol=3))
  dev.off()
}
