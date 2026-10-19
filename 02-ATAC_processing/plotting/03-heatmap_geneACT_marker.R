# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Seurat)
library(Signac)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggpubr)
library(patchwork)
library(ggrepel)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT by rerun_setup() instead of a
##   hard-coded setwd().
rerun_setup("02-ATAC_processing/plotting/",
            subdirs = c("data", "figures"))
signac.obj.INS <- readRDS("../sampleQC_redo/INS/data/signac_integrated_v2.RDS")
signac.obj.PFC <- readRDS("../sampleQC_redo/PFC/data/signac_integrated_v2.RDS")
## REVISION: VST reads the object carrying the nine MSN subtypes
##   (revision/tools/patch_msn_subtypes.R --target atac). The unpatched object labels every
##   striatal neuron "MSNs", which the palette and consensus_annotation_v3 grouping
##   below do not expect.
signac.obj.VST <- readRDS("../sampleQC_redo/VST/data/signac_integrated_v2_msn.RDS")


cp_col <- c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1",
  "purple1", "lightslateblue","orchid3","pink3","plum1",  "tomato1", "tomato3","sienna2", "burlywood3",
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
cell_types_consensus_levels_old <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                     "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                     "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                     "Astro","Oligo","OPC","Micro","Macrophages","Immune","PC","Endo","SMC","VLMC")
cell_types_consensus_levels <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                 "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Macrophages","Other immune","PC","Endo","SMC","VLMC")
names(cp_col) <- cell_types_consensus_levels


# heatmap visualization of the gene activities
Idents(signac.obj.INS) <- "consensus_annotation_v3"
da_gat.ins <- FindAllMarkers(
  object = signac.obj.INS,
  test.use = 'wilcox',
  only.pos = T
)
saveRDS(da_gat.ins,"./data/differential_geneACT_INS.RDS")
da_gat.ins <- readRDS("./data/differential_geneACT_INS.RDS")
top10.ins <-  da_gat.ins %>%
  group_by(cluster) %>%
  dplyr::filter(pct.1 > 0.5) %>%
  slice_max(avg_log2FC,n = 50) %>%
  ungroup() %>%
  mutate(cluster=factor(cluster,levels = cell_types_consensus_levels_old))%>%
  arrange(cluster)


Idents(signac.obj.PFC) <- "consensus_annotation_v3"
da_gat.pfc <- FindAllMarkers(
  object = signac.obj.PFC,
  test.use = 'wilcox',
  only.pos = T
)
saveRDS(da_gat.pfc,"./data/differential_geneACT_PFC.RDS")
da_gat.pfc <- readRDS("./data/differential_geneACT_PFC.RDS")
top10.pfc <-  da_gat.pfc %>%
  group_by(cluster) %>%
  dplyr::filter(pct.1 > 0.5) %>%
  slice_max(avg_log2FC,n = 50) %>%
  ungroup() %>%
  mutate(cluster=factor(cluster,levels = cell_types_consensus_levels_old))%>%
  arrange(cluster)



Idents(signac.obj.VST) <- "consensus_annotation_v3"
da_gat.vst <- FindAllMarkers(
  object = signac.obj.VST,
  test.use = 'wilcox',
  only.pos = T
)
saveRDS(da_gat.vst,"./data/differential_geneACT_VST.RDS")

da_gat.vst <- readRDS("./data/differential_geneACT_VST.RDS")
top10.vst <-  da_gat.vst %>%
  group_by(cluster) %>%
  dplyr::filter(pct.1 > 0.5) %>%
  slice_max(avg_log2FC,n = 50) %>%
  ungroup() %>%
  mutate(cluster=factor(cluster,levels = cell_types_consensus_levels_old))%>%
  arrange(cluster)



cell2keep <- data.frame("cells"=colnames(signac.obj.INS),
                        "cell_type"=signac.obj.INS@meta.data$consensus_annotation_v3)
cell2keep <- cell2keep %>% group_by(cell_type) %>% slice_sample(n=500)
signac.obj.INS.sub <- subset(signac.obj.INS,cells = cell2keep$cells)
atac_names <- rownames(signac.obj.INS.sub)
heatmap_atac_name_sorted <- atac_names[atac_names %in% unique(top10.ins$gene)]
signac.obj.INS.sub <- ScaleData(signac.obj.INS.sub,
                                features = heatmap_atac_name_sorted)

cell2keep <- data.frame("cells"=colnames(signac.obj.PFC),
                        "cell_type"=signac.obj.PFC@meta.data$consensus_annotation_v3)
cell2keep <- cell2keep %>% group_by(cell_type) %>% slice_sample(n=500)
signac.obj.PFC.sub <- subset(signac.obj.PFC,cells = cell2keep$cells)
atac_names <- rownames(signac.obj.PFC.sub)
heatmap_atac_name_sorted <- atac_names[atac_names %in% unique(top10.pfc$gene)]
signac.obj.PFC.sub <- ScaleData(signac.obj.PFC.sub,
                                features = heatmap_atac_name_sorted)


cell2keep <- data.frame("cells"=colnames(signac.obj.VST),
                        "cell_type"=signac.obj.VST@meta.data$consensus_annotation_v3)
cell2keep <- cell2keep %>% group_by(cell_type) %>% slice_sample(n=500)
signac.obj.VST.sub <- subset(signac.obj.VST,cells = cell2keep$cells)
atac_names <- rownames(signac.obj.VST.sub)
heatmap_atac_name_sorted <- atac_names[atac_names %in% unique(top10.vst$gene)]
signac.obj.VST.sub <- ScaleData(signac.obj.VST.sub,
                                features = heatmap_atac_name_sorted)



cell_types_consensus_levels_1 <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Macrophages","Immune","PC","Endo","SMC","VLMC")
cp_col1 <- c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1",
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
names(cp_col1) <- cell_types_consensus_levels_1

cell_types_consensus_levels_2 <- c("D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Macrophages","Immune","PC","Endo","SMC","VLMC")

cp_col2 <- c(
  "purple1", "lightslateblue","orchid3","pink3","plum1",  "tomato1", "tomato3","sienna2", "burlywood3",
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
names(cp_col2) <- cell_types_consensus_levels_2


signac.obj.INS.sub@meta.data$consensus_annotation_v3 <- factor(signac.obj.INS.sub@meta.data$consensus_annotation_v3,
                                                               levels = cell_types_consensus_levels_1)
signac.obj.PFC.sub@meta.data$consensus_annotation_v3 <- factor(signac.obj.PFC.sub@meta.data$consensus_annotation_v3,
                                                               levels = cell_types_consensus_levels_1)
signac.obj.VST.sub@meta.data$consensus_annotation_v3 <- factor(signac.obj.VST.sub@meta.data$consensus_annotation_v3,
                                                               levels = cell_types_consensus_levels_2)
Idents(signac.obj.INS.sub) <- "consensus_annotation_v3"
d1 <- DoHeatmap(signac.obj.INS.sub, 
                features = unique(top10.ins$gene),
                assay = "GeneACT",
                group.by = "consensus_annotation_v3",
                draw.lines = F,group.colors = cp_col1)+ 
  theme(axis.text.y = element_blank()) +
  scale_fill_gradientn(colors = c("white","black"))

Idents(signac.obj.PFC.sub) <- "consensus_annotation_v3"
d2 <- DoHeatmap(signac.obj.PFC.sub, 
                features = unique(top10.pfc$gene),
                assay = "GeneACT",
                group.by = "consensus_annotation_v3",
                draw.lines = F,group.colors = cp_col1)+ 
  theme(axis.text.y = element_blank()) +
  scale_fill_gradientn(colors = c("white","black"))

Idents(signac.obj.VST.sub) <- "consensus_annotation_v3"
d3 <- DoHeatmap(signac.obj.VST.sub, 
                features = unique(top10.vst$gene),
                assay = "GeneACT",
                group.by = "consensus_annotation_v3",
                draw.lines = F,group.colors = cp_col2)+ 
  theme(axis.text.y = element_blank()) +
  scale_fill_gradientn(colors = c("white","black"))


png("./figures/heatmap_gene_act_ins.png",width = 3000,height = 3000,res = 200)
print(d1)
dev.off()

png("./figures/heatmap_gene_act_pfc.png",width = 3000,height = 3000,res = 200)
print(d2)
dev.off()

png("./figures/heatmap_gene_act_vst.png",width = 3000,height = 3000,res = 200)
print(d3)
dev.off()



