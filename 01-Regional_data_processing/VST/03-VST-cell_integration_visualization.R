library(Seurat)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(cowplot)
options(future.globals.maxSize = 1e12)



wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/01-Regional_data_processing/VST/"
setwd(wd)

sce.obj <- readRDS("./data/02-cell_annotation/processed_sce_obj.rds")
count_data <- sce.obj@assays@data$counts
colnames(count_data) <- paste0(colnames(sce.obj),"-",sce.obj@colData$sample_ids)
seu.obj <- CreateSeuratObject(counts = count_data,
                              meta.data = as.data.frame(sce.obj@colData))
rm(count_data)
seu.obj <- NormalizeData(seu.obj)


# integration
seu.obj[["RNA"]] <- split(seu.obj[["RNA"]], f = seu.obj$sample_ids)
seu.obj <- NormalizeData(seu.obj)
seu.obj <- FindVariableFeatures(seu.obj)
seu.obj <- ScaleData(seu.obj)
seu.obj <- RunPCA(seu.obj)
seu.obj <- IntegrateLayers(
  object = seu.obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
seu.obj <- FindNeighbors(seu.obj, reduction = "harmony", dims = 1:50)
seu.obj <- FindClusters(seu.obj, resolution = 1, cluster.name = "seurat_clusters")


# visualization
seu.obj <-  RunTSNE(
  seu.obj,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path="/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)


# tocheck:
# 1. are there a lot of predicted MSNs also predicted as ext neurons with high scores
# 2. are these cells enriched in specific samples/conditions
# 3. are these cells more similar to INS than to PFC

# not sure:
# for the neuron population, cluster the data, then check the marker genes to make sure whether 
# each cluster is indeed ext or msn

p0 <- DimPlot(seu.obj,reduction = "tsne",group.by = "condition_new",label=F,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

p1 <- DimPlot(seu.obj,reduction = "tsne",group.by = "sample_ids",label=F,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p2 <- DimPlot(seu.obj,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p3 <- DimPlot(seu.obj,reduction = "tsne",group.by = "seurat_clusters",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)+ NoLegend()
Ext_markers <- c("GRIN1","SLC17A7","SATB2","SV2B")
Inh_markers <- c("GAD1","GAD2")
msn_markers <- c('PPP1R1B', 'PDE1B', 'BCL11B')
msn_sub_markers <- c("DRD2","STXBP6","DRD1","FAM163A","GRIA4","CPNE4","CASZ1","RXFP1","GRIK1")

f1 <- FeaturePlot(seu.obj,reduction = "tsne",features = Ext_markers)
f2 <- FeaturePlot(seu.obj,reduction = "tsne",features = Inh_markers)
f3 <- FeaturePlot(seu.obj,reduction = "tsne",features = msn_markers)

# so there are several clusters that are MSNs instead of Ext
# what are these clusters

# TODO: 1.identify the clusters that are ext based on the marker genes
#         verify by the max scores
#       2. visualize the sample/condition distributions
#       3. remove these clusters
v1 <- DotPlot(seu.obj,features = c(Ext_markers,Inh_markers,msn_markers))
ext_clusters <- c(14,22,26,33,35,38,41,43,45,47,48,50,52,56) # 35 not sure, 54 seem to be a mixture with MSN
seu.obj@meta.data$putative_ext <- seu.obj@meta.data$seurat_clusters %in% ext_clusters
DimPlot(seu.obj,reduction = "tsne",group.by = "putative_ext",label=F,
        raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
FeaturePlot(seu.obj,features = c("ma_prediction_score.max","nhp_prediction_score.max","biccn_prediction.score.max"))
table(seu.obj@meta.data$condition_new[seu.obj@meta.data$putative_ext])
table(seu.obj@meta.data$sample_ids[seu.obj@meta.data$putative_ext])
table(seu.obj@meta.data$sample_name[seu.obj@meta.data$putative_ext])

seu.obj@meta.data$celltypes <- "Others"
seu.obj@meta.data$celltypes[seu.obj@meta.data$putative_ext] <-"Ext"
metadata_all <- seu.obj@meta.data

DotPlot(seu.obj,features = Ext_markers,group.by = "celltypes")

n_counts <- metadata_all %>%
  group_by(condition_new,sample_ids,celltypes) %>%
  summarise(n=n())%>%
  filter(celltypes == "Ext")
cond_col <- c("#0072B2", "darkorange1")
names(cond_col) <- c("CTR","HIV")
ggplot(n_counts, aes(x=sample_ids,y=n,fill=condition_new)) +
  scale_fill_manual(values = cond_col)+
  facet_wrap(~condition_new, scales = "free_x")+
  geom_bar(stat="identity") +
  geom_text(aes(label=n), vjust=-1) +
  ylab("# Cells")+
  theme_minimal()+
  theme(axis.text = element_text(size=10),
        axis.text.x = element_text(angle = 45),
        strip.text = element_text(size=15,face = "bold"),
        legend.title=element_text(size=15,face = "bold"),
        title=element_text(size=15,face = "bold"),
        legend.text=element_text(size=10))+
  theme(axis.ticks.length.x = unit(0, "cm")) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))




DimPlot(seu.obj,reduction = "tsne",group.by = "biccn_predicted.id",label=T,
        raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
DimPlot(seu.obj,reduction = "tsne",group.by = "ma_predicted.id",label=T,
        raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
DimPlot(seu.obj,reduction = "tsne",group.by = "nhp_vst_predicted.id",label=T,
        raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

DotPlot(seu.obj,features = Ext_markers)

saveRDS(seu.obj,"./data/03-cell_integration_vis/integrated_seurat_obj.RDS")
seu.obj <- readRDS("./data/03-cell_integration_vis/integrated_seurat_obj.RDS")
