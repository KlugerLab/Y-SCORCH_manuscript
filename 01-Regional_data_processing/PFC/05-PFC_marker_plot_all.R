# stacked violin for all markers 

library(Seurat)

setwd("/banach2/SCORCH/scripts/combined_analysis_paper/01-Regional_data_processing/")
seurat.pfc <- readRDS("./PFC/data/03-cell_integration_vis/integrated_seurat_obj_v3.RDS")

cp_col1 <- c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1", 
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","gold3",
  "olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)

cell_types_consensus_levels <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Macrophages","Other immune","PC","Endo","SMC","VLMC")


names(cp_col1) <- cell_types_consensus_levels

# load markers

marker.list <- readRDS("../reference/fine_grain_markers_INS_PFC.rds")

# add macrophage markers along with the t markers
marker.list[["t"]] <- c(c("COLEC12","F13A1"),marker.list[["t"]])
marker.list[["micro"]] <- c( "ITGAM","CX3CR1","P2RY12")
markers <- c()
for (i in 1:length(marker.list)) {
  k <- ifelse(length(marker.list[[i]])>5,5,length(marker.list[[i]]))
  markers <- c(markers,marker.list[[i]][1:k])
}
markers.unique <- unique(markers)


pfc_data <- seurat.pfc@assays$RNA@layers$data
rownames(pfc_data) <- rownames(seurat.pfc)
colnames(pfc_data) <- colnames(seurat.pfc)
markers.unique <- markers.unique[markers.unique %in% rownames(pfc_data)]

stacked_violin_data <- as.data.frame(t(pfc_data[markers.unique,]))
stacked_violin_data$Cell <- rownames(stacked_violin_data)
stacked_violin_data$cell_types <- seurat.pfc@meta.data$consensus_annotation_v3_relabeled
stacked_violin_data$cell_types[stacked_violin_data$cell_types == "Other immune"] <- "Other\nimmune"
cell_types_consensus_levels <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
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


png(file="./PFC/figures/marker_and_abundance/markers_all.png",width = 15000,height = 12000,res=300)
print(p3)
dev.off()
