# REVISION: run on the CellBender-corrected data; paths (including fast_tsne_path) come
# from revision/config/paths.R and outputs go to OUT_ROOT. Changes from the published
# script are marked ## REVISION.

library(Seurat)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(cowplot)
options(future.globals.maxSize = 1e12)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

REGION <- "INS"

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup(file.path("01-Regional_data_processing", REGION),
            subdirs = "data/03-cell_integration_vis")

## REVISION: the published script used sample_meta_filter without defining it (it was
## carried over from the 02-cell_annotation session); it is now loaded explicitly.
sample_meta_filter <- load_sample_meta(REGION)

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
seu.obj <- FindClusters(seu.obj, resolution = 3, cluster.name = "seurat_clusters")

# visualization
seu.obj <-  RunTSNE(
  seu.obj,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path=FAST_TSNE_PATH, ann_not_vptree=FALSE, nthreads=12
)

p1 <- DimPlot(seu.obj,reduction = "tsne",group.by = "sample_ids",label=F,
        raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p2 <- DimPlot(seu.obj,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p3 <- DimPlot(seu.obj,reduction = "tsne",group.by = "seurat_clusters",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

# TODO: remove clusters composed of mostly one sample
sample_purity <- seu.obj@meta.data %>%
  group_by(seurat_clusters,sample_ids) %>%
  summarise(n = n()) %>%
  group_by(seurat_clusters) %>%
  mutate(n_cluster = sum(n)) %>%
  mutate(percent_sample = n/n_cluster) %>%
  group_by(seurat_clusters) %>%
  summarise(percent_sample_max = max(percent_sample))


cluster2remove <- unique(sample_purity[sample_purity$percent_sample_max > 0.5,]$seurat_clusters)


seu.obj_sub <- subset(seu.obj,cells = colnames(seu.obj)[!(Idents(seu.obj) %in% cluster2remove)])

p3_new <- DimPlot(seu.obj,reduction = "tsne",
              cells = colnames(seu.obj)[!(Idents(seu.obj) %in% cluster2remove)],
              group.by = "seurat_clusters",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

p2_new <- DimPlot(seu.obj,reduction = "tsne",
              cells = colnames(seu.obj)[!(Idents(seu.obj) %in% cluster2remove)],
              group.by = "consensus_annotation_v2",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p1_new <- DimPlot(seu.obj,reduction = "tsne",
                  cells = colnames(seu.obj)[!(Idents(seu.obj) %in% cluster2remove)],
                  group.by = "sample_ids",label=F,
                  raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

seu.obj_sub@meta.data$consensus_cell_type_coarse_l1 <- seu.obj_sub@meta.data$consensus_cell_type_coarse
seu.obj_sub@meta.data[seu.obj_sub@meta.data$consensus_cell_type_coarse_l1 == "Glial","consensus_cell_type_coarse_l1"] <- 
  as.character(seu.obj_sub@meta.data[seu.obj_sub@meta.data$consensus_cell_type_coarse_l1 == "Glial","consensus_annotation_v2"])
table(seu.obj_sub@meta.data$consensus_cell_type_coarse_l1)

# coarse markers
markers_coarse <- c("GRIN1","SLC17A7","GAB2","AQP4","MBP",
                    "PDGFRA","P2RY12","CLDN5","ANPEP","CD3D","CD3G","CD3E","COL1A2","ACTA2")
t_markers <- c("CD4","SKAP1","CD8A",
               "CD8B","THEMIS","GNLY","PRF1",
               "CD3D","CD3G","CD3E")
endo_markers <- c("CLDN5","ITM2A","PLTP","LY6A","FLT1","PTPRB")

DotPlot(seu.obj_sub,features = markers_coarse,group.by = "consensus_cell_type_coarse_l1")
DotPlot(seu.obj_sub,features = endo_markers,group.by = "consensus_cell_type_coarse_l1")
DotPlot(seu.obj_sub,features = markers_coarse,group.by = "seurat_clusters")
DotPlot(seu.obj_sub,features = t_markers,group.by = "consensus_cell_type_coarse_l1")

# re-integrate the data
seu.obj_sub <- FindVariableFeatures(seu.obj_sub)
seu.obj_sub <- ScaleData(seu.obj_sub)
seu.obj_sub <- RunPCA(seu.obj_sub)
seu.obj_sub <- IntegrateLayers(
  object = seu.obj_sub, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
# visualization
seu.obj_sub <-  RunTSNE(
  seu.obj_sub,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path=FAST_TSNE_PATH, ann_not_vptree=FALSE, nthreads=12
)

## REVISION (reviewer request): UMAP alongside the tSNE above.
##   Same harmony reduction and dims, so the two embeddings are directly
##   comparable. Stored on the object before it is saved. Not part of the
##   published analysis.
##
##   Threads: Seurat passes n_threads = future::nbrOfWorkers() to uwot, so a
##   direct n_threads argument to RunUMAP is silently dropped. run_umap_threaded()
##   sets the future plan for the call and restores it afterwards.
seu.obj_sub <- run_umap_threaded(seu.obj_sub, reduction = "harmony", dims = 1:50)

## REVISION: write named tSNE and UMAP panels. The published scripts call
##   DimPlot() with no device open, so under Rscript every plot lands unlabelled
##   in Rplots.pdf and the embedding figures cannot be found.
save_embedding_panels(seu.obj_sub, prefix = "03_INS_cell_integration_visualization")
p1 <- DimPlot(seu.obj_sub,reduction = "tsne",group.by = "sample_ids",label=F,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p2 <- DimPlot(seu.obj_sub,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
p3 <- DimPlot(seu.obj_sub,reduction = "tsne",group.by = "seurat_clusters",label=T,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

## REVISION: the published code compared against the strings "True"/"False"; read.csv()
## gives logicals here, so neither branch fired and condition_new stayed TRUE/FALSE.
seu.obj_sub@meta.data[,"condition_new"] <- sapply(seu.obj_sub@meta.data$sample_name,function(sid) sample_meta_filter[sample_meta_filter$sample_name == sid,"hiv_jay_v2"])
seu.obj_sub@meta.data$condition_new[seu.obj_sub@meta.data$condition_new == TRUE] <- "HIV"
seu.obj_sub@meta.data$condition_new[seu.obj_sub@meta.data$condition_new == FALSE] <- "CTR"
stopifnot(all(seu.obj_sub@meta.data$condition_new %in% c("HIV","CTR")))


p4 <- DimPlot(seu.obj_sub,reduction = "tsne",group.by = "condition_new",label=F,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)
DotPlot(seu.obj_sub,features = markers_coarse,group.by = "consensus_cell_type_coarse_l1")

# call the PC population Endo/PC?

# TODO: marker gene plot, TSNE

cp_col2 <- c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1", 
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
cell_types_consensus_levels <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Immune","PC","Endo","SMC","VLMC")
names(cp_col2) <- cell_types_consensus_levels

p2 <- DimPlot(seu.obj_sub,reduction = "tsne",group.by = "consensus_annotation_v2",label=T,
              cols = cp_col2,
              raster = F,shuffle = T,label.size = 5,pt.size = 0.5)

cell_types_coarse_levels <- c("Ext", "Inh" , "Astro" ,"Oligo" ,"OPC","Micro","Immune","Endo","PC","SMC","VLMC" )

seu.obj_sub@meta.data$consensus_cell_type_coarse_l1 <- factor(seu.obj_sub@meta.data$consensus_cell_type_coarse_l1,
                                                              levels = cell_types_coarse_levels)
markers_coarse <- c("GRIN1","SLC17A7","GAB2","AQP4","MBP",
                    "PDGFRA","P2RY12","THEMIS","CLDN5","ANPEP","ACTA2","COL1A2")

DotPlot(seu.obj_sub,features = markers_coarse,group.by = "consensus_cell_type_coarse_l1")

DotPlot(seu.obj_sub,features = t_markers,group.by = "consensus_cell_type_coarse_l1")

seu.obj_sub <- JoinLayers(seu.obj_sub)
saveRDS(seu.obj_sub,"./data/03-cell_integration_vis/integrated_seurat_obj.RDS")
