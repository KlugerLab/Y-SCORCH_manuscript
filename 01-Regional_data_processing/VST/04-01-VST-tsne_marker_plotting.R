library(Seurat)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(ggpubr)
options(future.globals.maxSize = 1e12)

seu.obj <- readRDS("./data/03-cell_integration_vis/integrated_seurat_obj.RDS")

# remove samples with abundant ext 
table(seu.obj@meta.data$sample_ids[seu.obj@meta.data$putative_ext])
sample2remove <- c("1039","1068","1119","7200776574","HCCRS","HCTME","HCTXJ")

seu.obj <- subset(seu.obj,cells = colnames(seu.obj)[!(seu.obj@meta.data$sample_ids %in% sample2remove)])


# remove the remaining exts
seu.obj <- subset(seu.obj,cells = colnames(seu.obj)[!seu.obj@meta.data$putative_ext])

p1 <- DimPlot(seu.obj,reduction = "tsne",group.by = "sample_ids",label=F,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p3 <- DimPlot(seu.obj,reduction = "tsne",group.by = "seurat_clusters",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)+ NoLegend()

Ext_markers <- c("GRIN1","SLC17A7","SATB2","SV2B")
f1 <- DotPlot(seu.obj,features = Ext_markers)

# coarse label the neuron cell types (MSNs, INHs, glials)
seu.obj@meta.data$nhp_vst_predicted.id
Inh_markers <- c("GAD1","GAD2")
msn_markers <- c('PPP1R1B', 'PDE1B', 'BCL11B')
msn_sub_markers <- c("DRD2","STXBP6","DRD1","FAM163A","GRIA4","CPNE4","CASZ1","RXFP1","GRIK1")
f2 <- DotPlot(seu.obj,features = Inh_markers)
f3 <- DotPlot(seu.obj,features = msn_sub_markers)

f4 <- DotPlot(seu.obj,features = c(Ext_markers,Inh_markers,msn_markers))
glial_clusters <- as.character(c(c(0:5),8,10,11,12,13,15,17,18,20,
                    21,24,28,34,36,37,39,40,44,46,49,
                    53,57,58,59,61,62,63))
msn_clusters <- as.character(c(6,7,9,16,25))
inh_clusters <- as.character(c(19,23,27,29,30,31,32,42,51)) #32 seem to be mixed of ext/inh (2 cells only)
names(glial_clusters) <- rep("Glial",length(glial_clusters))
names(msn_clusters) <- rep("MSNs",length(msn_clusters))
names(inh_clusters) <- rep("Inh",length(inh_clusters))
manual_coarse_annotations <- c(inh_clusters,msn_clusters,glial_clusters)
manual_coarse_annotations <- setNames(names(manual_coarse_annotations), manual_coarse_annotations)

seu.obj@meta.data$manual_cell_type_coarse <- manual_coarse_annotations[as.character(seu.obj@meta.data$seurat_clusters)]
p2 <- DimPlot(seu.obj,reduction = "tsne",group.by = "consensus_cell_type_coarse",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

# remove the cells with low prediction scores
seu.obj@meta.data$cell2remove = ifelse((seu.obj@meta.data$nhp_prediction_score.max < 0.5) & 
                                         (seu.obj@meta.data$biccn_prediction.score.max < 0.5) &
                                         (seu.obj@meta.data$ma_prediction.score.max  < 0.5),T,F)

seu.obj <- subset(seu.obj,cells = colnames(seu.obj)[!seu.obj@meta.data$cell2remove])

# check consistency
all_cp_unique_ordered <- c("L2/3 IT","L2-3 IT",
                           "L4 IT","L5 IT","L3-5 IT-1","L3-5 IT-2","L3-5 IT-3",
                           "L5 ET","L5/6 NP","L5-6 NP",
                           "L6 CT","L6 IT","L6 IT Car3","L6 IT-1","L6 IT-2","L6b","L6B",
                           "Lamp5","LAMP5 RELN","Lamp5 Lhx6","LAMP5 LHX6",
                           "Chandelier","PVALB ChC","Pax6","ADARB2 KCNG1","Pvalb","PVALB",
                           "Sncg","Sst","SST","Sst Chodl","SST HGF","SST NPY",
                           "Vip","VIP",
                           "Interneurons","MSNs",
                           "Astro","Astrocytes","Endo","Endothelial","Oligo","Oligos","Oligos_Pre","OPC","Microglia","Micro/PVM","Micro","Immune","VLMC","PC",
                           "SMC","RB","Mural/Fibroblast")
coarse_cp_mapper <- c(rep("Ext",17),rep("Inh",18),"Interneurons","MSNs",rep("Glial",17))
names(coarse_cp_mapper) <- all_cp_unique_ordered


seu.obj@meta.data$consistent_labeled <- FALSE
seu.obj@meta.data$consistent_labeled[(seu.obj@meta.data$manual_cell_type_coarse == "Glial")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$biccn_predicted.id)] == "Glial")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$ma_predicted.id)] == "Glial")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$nhp_vst_predicted.id)] == "Glial")
                                     ] <- TRUE

seu.obj@meta.data$consistent_labeled[(seu.obj@meta.data$manual_cell_type_coarse == "Inh")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$biccn_predicted.id)] == "Inh")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$ma_predicted.id)] == "Inh")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$nhp_vst_predicted.id)] == "Interneurons")
] <- TRUE

seu.obj@meta.data$consistent_labeled[(seu.obj@meta.data$manual_cell_type_coarse == "MSNs")&
                                       (coarse_cp_mapper[as.character(seu.obj@meta.data$nhp_vst_predicted.id)] == "MSNs")
] <- TRUE

seu.obj <- subset(seu.obj,cells = colnames(seu.obj)[seu.obj@meta.data$consistent_labeled])



# Glial cells labeled by Ma et al, inh labeled by biccn, MSNs labeled by non-human primate

# glial
seu.obj@meta.data$consensus_annotation_v2 <-  seu.obj@meta.data$manual_cell_type_coarse
seu.obj@meta.data[seu.obj@meta.data$consensus_annotation_v2 == "Glial","consensus_annotation_v2"] <- 
  as.character(seu.obj@meta.data[seu.obj@meta.data$consensus_annotation_v2 == "Glial","ma_predicted.id"])

  
seu.obj@meta.data$consensus_cell_type_coarse_l1 <- seu.obj@meta.data$consensus_annotation_v2

# inh
seu.obj@meta.data[seu.obj@meta.data$consensus_annotation_v2 == "Inh","consensus_annotation_v2"] <- 
  as.character(seu.obj@meta.data[seu.obj@meta.data$consensus_annotation_v2 == "Inh","biccn_predicted.id"])

# vst 
seu.obj.msn <- subset(seu.obj,cells = colnames(seu.obj)[seu.obj@meta.data$consensus_annotation_v2 == "MSNs"])
saveRDS(seu.obj.msn,"./data/msn_subset.RDS")
# run 04-02-VST-msn_subtyping.R

prediction.df <- readRDS("./data/msn_prediction.RDS")
prediction.df <- bind_rows(prediction.df)
prediction.df$sample_ids <- sapply(rownames(prediction.df), function(sname){strsplit(sname,split = '-')[[1]][3]})

seu.obj@meta.data$msn_subtypes <- NULL
seu.obj@meta.data[rownames(prediction.df),"msn_subtypes"] <- prediction.df[,"predicted.id"]

seu.obj@meta.data[rownames(prediction.df),"consensus_annotation_v2"]<- prediction.df[,"predicted.id"]

# remove the MSNs that are from samples with < 50 MSN cells
seu.obj <- subset(seu.obj,cells = colnames(seu.obj)[seu.obj@meta.data$consensus_annotation_v2 != "MSNs"])


# re-integrate the data
# re-integrate the data
seu.obj <- FindVariableFeatures(seu.obj)
seu.obj <- ScaleData(seu.obj)
seu.obj <- RunPCA(seu.obj)
seu.obj <- IntegrateLayers(
  object = seu.obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
# visualization
seu.obj <-  RunTSNE(
  seu.obj,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path="/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)

p2 <- DimPlot(seu.obj,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
                               

seu.obj <- JoinLayers(seu.obj)

saveRDS(seu.obj,"./data/04-VST-tsne_marker_plotting/integrated_seurat_obj.RDS")



