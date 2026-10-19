# stacked violin for all markers 
#
# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Seurat)
## REVISION: packages used but never loaded in the published script (they were
##   loaded earlier in the interactive session); under Rscript this is an error.
library(ggplot2)
## REVISION: also used but never loaded.
library(cowplot)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup("01-Regional_data_processing",
            subdirs = c("VST/data/04-VST-tsne_marker_plotting", "VST/figures/marker_and_abundance"))
seurat.vst <- readRDS("./VST/data/04-VST-tsne_marker_plotting/integrated_seurat_obj_v3.RDS")

cp_col1 <- c(
  "purple1", "lightslateblue","orchid3","pink3","plum1",  "tomato1", "tomato3","sienna2", "burlywood3", 
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
cell_types_consensus_levels <- c("D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Macrophages","Other immune","PC","Endo","SMC","VLMC")


names(cp_col1) <- cell_types_consensus_levels


# load markers

marker.list <- readRDS(file.path(REFERENCE_DIR, "fine_grain_markers_INS_PFC.rds"))

# add macrophage markers along with the t markers
marker.list[["t"]] <- c(c("COLEC12","F13A1"),marker.list[["t"]])
marker.list[["micro"]] <- c( "ITGAM","CX3CR1","P2RY12")
msn_markers <- c('PPP1R1B', 'PDE1B', 'BCL11B')
msn_sub_markers <- c("DRD2","STXBP6","DRD1","FAM163A","GRIA4","CPNE4","CASZ1","RXFP1","GRIK1")

markers <- c(msn_markers,msn_sub_markers,marker.list[["lamp5_biccn"]][1:5],
             marker.list[["lamp5_lhx6_biccn"]][1:5],marker.list[["pax6_biccn"]][1:5],
             marker.list[["pvalb_biccn"]][1:5],marker.list[["chandelier_biccn"]][1:5],marker.list[["sncg_biccn"]][1:5],
             marker.list[["sst_biccn"]][1:5],marker.list[["sst_chodl_biccn"]][1:5],
             marker.list[["vip_biccn"]][1:5],marker.list[["astro"]][1:5],marker.list[["oligo"]][1:5],
             marker.list[["OPC"]][1:5],marker.list[["micro"]][1:3],marker.list[["t"]][1:5],
             marker.list[["peri"]][1:5], marker.list[["endo"]][1:5],marker.list[["SMC"]][1:5],marker.list[["VLMC"]][1:5])
markers.unique <- unique(markers)



vst_data <- seurat.vst@assays$RNA@layers$data
rownames(vst_data) <- rownames(seurat.vst)
colnames(vst_data) <- colnames(seurat.vst)
markers.unique <- markers.unique[markers.unique %in% rownames(vst_data)]

stacked_violin_data <- as.data.frame(t(vst_data[markers.unique,]))
stacked_violin_data$Cell <- rownames(stacked_violin_data)
stacked_violin_data$cell_types <- seurat.vst@meta.data$consensus_annotation_v3_relabeled
stacked_violin_data$cell_types[stacked_violin_data$cell_types == "Other immune"] <- "Other\nimmune"
cell_types_consensus_levels <- c("D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Macrophages","Other\nimmune","PC","Endo","SMC","VLMC")

stacked_violin_data$cell_types <- factor(stacked_violin_data$cell_types,levels = cell_types_consensus_levels)
# Use melt to change data.frame format
stacked_violin_data <- reshape2::melt(stacked_violin_data, id.vars = c("Cell","cell_types"), measure.vars = markers.unique,
                                      variable.name = "Feat", value.name = "Expr")
p3 <- ggplot(stacked_violin_data, aes(Expr, factor(Feat), fill = cell_types)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(cols = vars(cell_types), scales = "free") +
  theme_cowplot(font_size = 12) +
  scale_fill_manual(values=cp_col1)+
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size=20,face = "bold")) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=30),
        
  ) + xlab("Expression Level") + ylab("Marker genes")

png(file="./VST/figures/marker_and_abundance/markers_all.png",width = 18000,height = 12000,res=300)
print(p3)
dev.off()


