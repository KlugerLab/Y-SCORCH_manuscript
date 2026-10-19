# TSNE visualization of the ATAC data
#
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
## REVISION: ggplot2 is used but was not loaded by the published script;
##   loaded explicitly so the script runs under Rscript.
library(ggplot2)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT by rerun_setup() instead of a
##   hard-coded setwd().
rerun_setup("02-ATAC_processing/plotting",
            subdirs = c("figures/tsne_new"))
signac.obj.INS <- readRDS("../sampleQC_redo/INS/data/signac_integrated_v2.RDS")
signac.obj.PFC <- readRDS("../sampleQC_redo/PFC/data/signac_integrated_v2.RDS")
## REVISION: VST reads the object carrying the nine MSN subtypes
##   (revision/tools/patch_msn_subtypes.R --target atac). The unpatched object labels every
##   striatal neuron "MSNs", which the palette and consensus_annotation_v3 grouping
##   below do not expect.
signac.obj.VST <- readRDS("../sampleQC_redo/VST/data/signac_integrated_v2_msn.RDS")

signac.obj.INS@meta.data$consensus_annotation_v3_relabeled <- signac.obj.INS@meta.data$consensus_annotation_v3
signac.obj.INS@meta.data$consensus_annotation_v3_relabeled[signac.obj.INS@meta.data$consensus_annotation_v3 == "Immune"] <- "Other immune"

signac.obj.PFC@meta.data$consensus_annotation_v3_relabeled <- signac.obj.PFC@meta.data$consensus_annotation_v3
signac.obj.PFC@meta.data$consensus_annotation_v3_relabeled[signac.obj.PFC@meta.data$consensus_annotation_v3 == "Immune"] <- "Other immune"

signac.obj.VST@meta.data$consensus_annotation_v3_relabeled <- signac.obj.VST@meta.data$consensus_annotation_v3
signac.obj.VST@meta.data$consensus_annotation_v3_relabeled[signac.obj.VST@meta.data$consensus_annotation_v3 == "Immune"] <- "Other immune"


# color
cp_col <- c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1",
  "purple1", "lightslateblue","orchid3","pink3","plum1",  "tomato1", "tomato3","sienna2", "burlywood3",
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
cell_types_consensus_levels <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                 "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Macrophages","Other immune","PC","Endo","SMC","VLMC")
names(cp_col) <- cell_types_consensus_levels



region_cols <- c("#E69F00", "#56B4E9","#CC79A7")
names(region_cols) <-c("INS","PFC","VST")

cond_col <- c("#0072B2", "darkorange1")
names(cond_col) <- c("CTR","HIV")


# my legend
legend_df <- as.data.frame(cp_col)
legend_df$cp <- factor(rownames(legend_df),levels=cell_types_consensus_levels)
legend_df$somev <- 1
mylegend <- as_ggplot(ggpubr::get_legend(ggplot(legend_df,aes(x=cp,y=somev,color=cp))+
                                           geom_point()+
                                           guides(color=guide_legend(nrow=3,title=NULL,override.aes = list(size=10)))+
                                           
                                           theme_void()+
                                           theme(legend.title=element_blank(),legend.text = element_text(size=20))+
                                           scale_color_manual(values = cp_col),position="bottom"))


t1 <- DimPlot(signac.obj.INS,reduction = "tsne",group.by = "consensus_annotation_v3_relabeled",
              cols = cp_col,
              label=T,repel = T,
              raster = F,shuffle = T,label.size = 10,pt.size = 0.5) +
  theme_void()+
  NoLegend()+
  theme(title = element_blank())

t1_legend <- as_ggplot(ggpubr::get_legend(DimPlot(signac.obj.INS,reduction = "tsne",group.by = "consensus_annotation_v3_relabeled",
                                                  cols = cp_col,
                                                  label=T,repel = T,
                                                  raster = F,shuffle = T,label.size = 10,pt.size = 0.5) +
                                            theme_void()+
                                            
                                            theme(title = element_blank()),position="bottom"))



t2 <- DimPlot(signac.obj.PFC,reduction = "tsne",group.by = "consensus_annotation_v3_relabeled",
              cols = cp_col,
              label=T,repel = T,
              raster = F,shuffle = T,label.size = 10,pt.size = 0.5) +
  theme_void()+
  NoLegend()+
  theme(title = element_blank())

t2_legend <- as_ggplot(ggpubr::get_legend(DimPlot(signac.obj.PFC,reduction = "tsne",group.by = "consensus_annotation_v3_relabeled",
                                                  cols = cp_col,
                                                  label=T,repel = T,
                                                  raster = F,shuffle = T,label.size = 10,pt.size = 0.5) +
                                            theme_void()+
                                            
                                            theme(title = element_blank()),position="bottom"))


vst_tsne_df <- as.data.frame(signac.obj.VST@reductions$tsne@cell.embeddings)
vst_tsne_df$labels <- signac.obj.VST@meta.data$consensus_annotation_v3_relabeled
t3 <- DimPlot(signac.obj.VST,reduction = "tsne",group.by = "consensus_annotation_v3_relabeled",
              cols = cp_col,
              label=F,repel = T,
              
              raster = F,shuffle = T,label.size = 10,pt.size = 0.5) +
  #geom_text_repel(data=vst_tsne_df,aes(x=tSNE_1,y=tSNE_2,label=labels),max.overlaps = 5)+
  theme_void()+
  NoLegend()+
  theme(title = element_blank())


t3 <- LabelClusters(t3,clusters = unique(signac.obj.VST@meta.data$consensus_annotation_v3_relabeled),
                    id = "consensus_annotation_v3_relabeled",repel = T,max.overlaps=5,
                    size = 10)

t3_legend <- as_ggplot(ggpubr::get_legend(DimPlot(signac.obj.VST,reduction = "tsne",group.by = "consensus_annotation_v3_relabeled",
                                                  cols = cp_col,
                                                  label=T,repel = T,
                                                  raster = F,shuffle = T,label.size = 10,pt.size = 0.5) +
                                            theme_void()+
                                            
                                            theme(title = element_blank()),position="bottom"))

# save the plots

png(file="./figures/tsne_new/tsne_INS_cp.png",width = 3500,height = 3500,res=300)
print(t1)
dev.off()


png(file="./figures/tsne_new/tsne_INS_cp_legend.png",width = 1500,height = 1500,res=300)
print(t1_legend)
dev.off()


png(file="./figures/tsne_new/tsne_PFC_cp.png",width = 3500,height = 3500,res=300)
print(t2)
dev.off()


png(file="./figures/tsne_new/tsne_PFC_cp_legend.png",width = 1500,height = 1500,res=300)
print(t2_legend)
dev.off()



png(file="./figures/tsne_new/tsne_VST_cp.png",width = 3500,height = 3500,res=300)
print(t3)
dev.off()


png(file="./figures/tsne_new/tsne_VST_cp_legend.png",width = 1500,height = 1500,res=300)
print(t3_legend)
dev.off()



p1 <- FeaturePlot(signac.obj.INS,features = c("SATB2", "GAD2","PDGFRA","P2RY12"),
            min.cutoff = "q10",max.cutoff = "q90")&
  theme_void()&theme(title = element_text(size=25), 
                     legend.text=element_text(size=20))

p2 <- FeaturePlot(signac.obj.PFC,features = c("SATB2", "GAD2","PDGFRA","P2RY12"),
                  min.cutoff = "q10",max.cutoff = "q90")&
  theme_void()&theme(title = element_text(size=25), 
                    legend.text=element_text(size=20))

p3 <- FeaturePlot(signac.obj.VST,features = c("BCL11B", "GAD2","PDGFRA","P2RY12"),
                  min.cutoff = "q10",max.cutoff = "q90")&
  theme_void()&theme(title = element_text(size=25), 
                     legend.text=element_text(size=20))

png(file="./figures/tsne_new/tsne_PFC_markers.png",width = 3500,height = 3500,res=300)
print(p2)
dev.off()


png(file="./figures/tsne_new/tsne_INS_markers.png",width = 3500,height = 3500,res=300)
print(p1)
dev.off()

png(file="./figures/tsne_new/tsne_VST_markers.png",width = 3500,height = 3500,res=300)
print(p3)
dev.off()

