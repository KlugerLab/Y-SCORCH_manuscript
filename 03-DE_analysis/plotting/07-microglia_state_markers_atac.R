# ==============================================================================
# REVISION: authors' plotting script, added to this repository with the revision; run on
# the CellBender-corrected data; paths come from revision/config/paths.R and
# outputs go to OUT_ROOT. Changes are marked ## REVISION.
#
# Figure S9b (gene activities): microglia state markers across samples, columns grouped by
# condition then region, rows the marker panels (homeostatic, anti-inflammatory,
# interferon-responsive, antigen presentation, proliferative, disease-associated).
# Only samples with more than 30 microglia are shown.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(limma)
library(scater)
library(muscat)
library(purrr)
library(DESeq2)
library(Seurat)
library(ggpubr)
library(clusterProfiler)
library(enrichR)
library(enrichplot)
organism = "org.Hs.eg.db"

library(organism, character.only = TRUE)
library(Homo.sapiens)
library(RColorBrewer)
library(circlize)
library(grid)

#BiocManager::install("GOSim")
## REVISION: GOSim is no longer in Bioconductor and no GOSim function is called
##   in this script, so the import is disabled.
# library(GOSim)
## REVISION: working directory is set by rerun_setup() under OUT_ROOT instead of a hardcoded server path.
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
rerun_setup("03-DE_analysis/cross_region/",
            subdirs = c("data", "figures/microglia_state_markers"))
edgeR_res <- readRDS("./data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")  ## REVISION: only the _withFreq table is written
pb.norm_mean <- readRDS("./data/pseudo_norm_mean_bulk.rds")








pb.norm_mean_act_ins <- readRDS("../../04-ATAC_DE/INS/data_GeneACT/pseudo_norm_mean_bulk.rds")
pb.norm_mean_act_pfc <- readRDS("../../04-ATAC_DE/PFC/data_GeneACT/pseudo_norm_mean_bulk.rds")
pb.norm_mean_act_vst <- readRDS("../../04-ATAC_DE/VST/data_GeneACT/pseudo_norm_mean_bulk.rds")
#only keep microglia
pb.norm_mean_act_ins@assays@data <- SimpleList(list("Micro"=pb.norm_mean_act_ins@assays@data$Micro))
pb.norm_mean_act_pfc@assays@data <- SimpleList(list("Micro"=pb.norm_mean_act_pfc@assays@data$Micro))
pb.norm_mean_act_vst@assays@data <- SimpleList(list("Micro"=pb.norm_mean_act_vst@assays@data$Micro))
pb.norm_mean_act <- cbind(cbind(pb.norm_mean_act_ins,pb.norm_mean_act_pfc),pb.norm_mean_act_vst)

coldf <- data.frame(
  region = pb.norm_mean_act@colData$region,
  condition = pb.norm_mean_act@colData$group_id,
  sample_id = pb.norm_mean_act@colData$sample_name)

rownames(coldf) <- coldf$sample_id
## REVISION: condition is relabelled PWoH / PWH, as in every other revision figure.
##   The factor keeps PWoH first, so the column order matches the published CTR-then-HIV
##   layout (plain alphabetical sorting would put PWH first).
coldf$condition <- factor(as.character(pwh_relabel(coldf$condition)), levels = PWH_LEVELS)
coldf <- coldf %>%
  arrange(condition,region,sample_id)

condition_cols <- c("#0072B2","darkorange1")
names(condition_cols) <- PWH_LEVELS   ## REVISION: PWoH / PWH, same colours

region_cols <- c("#E69F00", "#56B4E9","#CC79A7")
names(region_cols) <-c("INS","PFC","VST")


ann_colors = list(
  condition = condition_cols,
  region = region_cols
)
target_genes <- c("P2RY12","CX3CR1","TMEM119","SALL1",
                  "IL4", "IL13", "IL10", "TGFB1", "TGFB3",
                  "IFITM3","IFIT3","ISG15","CD74","HLA-DRA",
                  "MKI67","PCNA","APOE", "LPL")

pb.norm_mean.data.sub <- pb.norm_mean_act@assays@data[["Micro"]]
colnames(pb.norm_mean.data.sub) <- pb.norm_mean_act@colData$sample_name
pb.norm_mean.data.sub <- pb.norm_mean.data.sub[,rownames(coldf)]

# only keep samples with enough microglia
n_cells <- do.call(rbind,pb.norm_mean_act@int_colData$n_cells)
rownames(n_cells) <- pb.norm_mean_act@colData$sample_name
n_cells <- n_cells[rownames(coldf),]
s2keep <- n_cells[,"Micro"] > 30



mat <- t(scale(t(pb.norm_mean.data.sub[target_genes,s2keep])))
coldf2keep <- coldf[s2keep,]
colors <- c(-4.1,seq(-4,4,by=0.001),4.1)
my_palette <- c("steelblue",colorRampPalette(colors = c("steelblue", "white", "red1"))
                (n = length(colors)-3), "red1")
col_fun = colorRamp2(c(-4, 0, 4), c("steelblue", "white", "red1"))
p1 <- ComplexHeatmap::Heatmap(mat, 
                              col=col_fun,
                              column_names_gp = gpar(fontsize = 25),
                              row_names_gp = gpar(fontsize = 30),
                              cluster_columns = F, cluster_rows = F,
                              show_row_dend = F,
                              name=" ",
                              column_gap = unit(0.1, "mm"),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(
                                "Condition"=coldf2keep$condition,
                                "Region" = coldf2keep$region,
                                col=list(
                                  Region = region_cols,
                                  Condition = condition_cols
                                  
                                ),
                                
                                simple_anno_size = unit(2, "cm"),
                                annotation_name_gp = gpar(fontsize = 30,fontface="bold"),
                                annotation_legend_param = list(
                                  Region = list(
                                    labels_gp = gpar(fontsize = 25),
                                    title_gp = gpar(fontsize = 30,fontface="bold")),
                                  Condition = list(
                                    labels_gp = gpar(fontsize = 25),
                                    title_gp = gpar(fontsize = 30,fontface="bold"))
                                )
                              ),
                              heatmap_legend_param = list(
                                labels_gp = gpar(fontsize = 25),
                                title_gp = gpar(fontsize = 25)
                              ),
                              
                              show_column_names = F)

png(filename = paste0("./figures/microglia_state_markers/microglia_markers_act.png"),width=3500,height = 3500,res = 300)
print(p1)
dev.off()



## REVISION: disabled. Everything below is exploratory and saves nothing; it
##   reloads the 15 GB PFC Signac object and redoes TF-IDF / SVD / Harmony / tSNE /
##   FindAllMarkers (hours) only to draw interactive plots. Not part of S9.
if (FALSE) {
# check the expression of the interferon-response marker expressions on tsne
signac.pfc <- readRDS("../../02-ATAC_processing/sampleQC_redo/PFC/data/signac_integrated_v2.RDS")
signac.pfc.micro <- subset(signac.pfc,cells = colnames(signac.pfc)[signac.pfc@meta.data$consensus_annotation_v3 == "Micro"])


DefaultAssay(signac.pfc.micro) <- "ATAC"
signac.pfc.micro <- RunTFIDF(signac.pfc.micro)
signac.pfc.micro <- FindTopFeatures(signac.pfc.micro)
signac.pfc.micro <- RunSVD(signac.pfc.micro)

# this is not optimal: first component not ruled out
signac.pfc.micro <- RunHarmony(
  object = signac.pfc.micro,
  group.by.vars = 'sample_name',
  reduction.use = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

signac.pfc.micro <-  RunTSNE(
  signac.pfc.micro,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path="/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
DefaultAssay(signac.pfc.micro) <- "GeneACT"
signac.pfc.micro <- FindNeighbors(signac.pfc.micro, reduction = "harmony", dims = 1:50)
signac.pfc.micro <- FindClusters(signac.pfc.micro, graph.name = "ATAC_snn",resolution = 1, 
                                 cluster.name = "seurat_clusters")
cluster.markers <- FindAllMarkers(signac.pfc.micro,assay = "GeneACT",only.pos = T)
cluster.markers.top10 <- cluster.markers %>%
  group_by(cluster) %>%
  dplyr::filter((pct.1 > 0.3) & (pct.2 < 0.1)) %>%
  slice_max(avg_log2FC,n = 10)

FeaturePlot(signac.pfc.micro,features = c("NUPR1","IL15","IFITM3","IFIT3","ISG15"))
FeaturePlot(signac.pfc.micro,features = c("nCount_ATAC"))

DimPlot(signac.pfc.micro,group.by = "sample_name",shuffle = T)+NoLegend()
DimPlot(signac.pfc.micro,group.by = "condition_new",shuffle = T)
DimPlot(signac.pfc.micro,group.by = "seurat_clusters",shuffle = T,label = T)


VlnPlot(signac.pfc.micro,features = c("IL15"),group.by = "sample_name")+
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+NoLegend()

}
