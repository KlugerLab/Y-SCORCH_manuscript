# for the glial cells
# remove Endo/Peri
# remove ambiguous clusters

library(Seurat)

wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/01-Regional_data_processing/"
setwd(wd)

seurat.combined <- readRDS("./combined/data/seurat_combined.RDS")

seurat.glial <- subset(seurat.combined,
                       cells = colnames(seurat.combined)[!(seurat.combined@meta.data$consensus_cell_type_coarse_l1 %in% c("Ext","Inh","MSNs"))])

# re-integrate the data, cluster the data, and annotate again
DefaultAssay(seurat.glial) <- "RNA"
seurat.glial <- JoinLayers(seurat.glial)
seurat.glial[["RNA"]] <- split(seurat.glial[["RNA"]], f = seurat.glial$sample_name)
seurat.glial <- NormalizeData(seurat.glial)
seurat.glial <- FindVariableFeatures(seurat.glial)
seurat.glial <- ScaleData(seurat.glial)
seurat.glial <- RunPCA(seurat.glial)
seurat.glial <- IntegrateLayers(
  object = seurat.glial, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
# cluster the data
seurat.glial <- FindNeighbors(seurat.glial, reduction = "harmony", dims = 1:50)
seurat.glial <- FindClusters(seurat.glial, resolution = 2, cluster.name = "seurat_clusters")
# visualization
seurat.glial <-  RunTSNE(
  seurat.glial,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 100,
  fast_tsne_path="/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)

p1 <- DimPlot(seurat.glial,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

p2 <- DimPlot(seurat.glial,reduction = "tsne",group.by = "seurat_clusters",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)+NoLegend()

manual_markers <- c("AQP4","MBP","PDGFRA","P2RY12","SKAP1","CLDN5","GRM8","ACTA2","COL1A2")
d1 <- DotPlot(seurat.glial,features = manual_markers)

f1 <- FeaturePlot(seurat.glial,features = c("CLDN5","GRM8"))

# these are glial clusters that look ambiguous in terms of marker expressions
markers <- readRDS("../reference/fine_grain_markers_INS_PFC.rds")

d2 <- DotPlot(seurat.glial,features = markers[["oligo"]])
f2 <- FeaturePlot(seurat.glial,features = markers[["oligo"]])

f3 <- FeaturePlot(seurat.glial,features = c("COLEC12","F13A1"))


cluster2remove <- c(21,24,26,28,33,37,39,40,42,43,47,48,50,53,54)

seurat.glial.sub <- subset(seurat.glial,cells = colnames(seurat.glial)[!(Idents(seurat.glial) %in% cluster2remove)])

p1_new <- DimPlot(seurat.glial.sub,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
                  raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

seurat.glial.sub <- JoinLayers(seurat.glial.sub)
markers.38 <- FindMarkers(seurat.glial.sub,ident.1 = 38,ident.2 = NULL)
markers.38.top10 <- markers.38 %>%
  slice_max(avg_log2FC,n=10)
# these are all macrophage markers

cells2keep <- colnames(seurat.combined)[!(colnames(seurat.combined) %in%  colnames(seurat.glial)[(Idents(seurat.glial) %in% cluster2remove)])]
macrophage_cells <- colnames(seurat.glial)[Idents(seurat.glial) == 38]

# subset the original data, and rename macrophages in V3
seurat.combined.sub <- subset(seurat.combined,cells = cells2keep)
seurat.combined.sub@meta.data$consensus_annotation_v3 <- as.character(seurat.combined.sub@meta.data$consensus_annotation_v2)
seurat.combined.sub@meta.data$consensus_annotation_v3[colnames(seurat.combined.sub) %in% macrophage_cells] <- "Macrophages"

p4 <- DimPlot(seurat.combined.sub,reduction = "tsne",group.by = "consensus_annotation_v3",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
saveRDS(seurat.combined.sub,"./combined/data/seurat_combined_v2.RDS")


