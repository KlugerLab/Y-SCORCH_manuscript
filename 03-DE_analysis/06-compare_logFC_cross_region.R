# compare the DEs between regions

library(dplyr)
library(ggplot2)
library(limma)
library(scater)
library(muscat)
library(purrr)
library(DESeq2)
library(Seurat)
library(ggpubr)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(enrichR)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(Homo.sapiens)
library(RColorBrewer)
library(circlize)
library(grid)
library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(scater)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggrepel)


setwd("/banach2/SCORCH/scripts/combined_analysis_paper/03-DE_analysis/")

INS_edgeR <- readRDS("./INS/results/INS_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
INS_all <- INS_edgeR %>%
  ungroup() %>%
  dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05))%>%
  dplyr::select(gene,cluster_id,logFC,logCPM,p_adj.glb,region) 

PFC_edgeR <- readRDS("./PFC//results/PFC_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
PFC_all <- PFC_edgeR %>%
  ungroup() %>%
  dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05))%>%
  dplyr::select(gene,cluster_id,logFC,logCPM,p_adj.glb,region) 

VST_edgeR <- readRDS("./VST//results/VST_HIV-CTR_res_edgeR_filtered_withFreq.RDS")

VST_all <- VST_edgeR %>%
  ungroup() %>%
  dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05))%>%
  dplyr::select(gene,cluster_id,logFC,logCPM,p_adj.glb,region)



DE_combined <- INS_all %>%
  inner_join(PFC_all,by = c("cluster_id","gene"))

DE_INS_PFC <- INS_all %>%
  inner_join(PFC_all,by = c("cluster_id","gene"))
DE_INS_PFC <- as.data.frame(DE_INS_PFC)

DE_INS_VST <- INS_all %>%
  inner_join(VST_all,by = c("cluster_id","gene"))
DE_INS_VST <- as.data.frame(DE_INS_VST)

DE_PFC_VST <- PFC_all %>%
  inner_join(VST_all,by = c("cluster_id","gene"))
DE_PFC_VST <- as.data.frame(DE_PFC_VST)

# split by cell types
cp_col <- c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1",
  "purple1", "lightslateblue","orchid3","pink3","plum1",  "tomato1", "tomato3","sienna2", "burlywood3",
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
cell_types_consensus_levels <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                 "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Immune","PC","Endo","SMC","VLMC")
neuron_cps <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip")

names(cp_col) <- cell_types_consensus_levels

cp_col2 <-  c("tomato2","royalblue3","dodgerblue","gold2")
names(cp_col2) <- c("Ext","MSNs","Inh","Glial")

cell_types_coarse_levels <- c(rep("Ext",9),rep("MSNs",9),
                              rep("Inh",9),rep("Glial",10))
names(cell_types_coarse_levels) <- cell_types_consensus_levels

# 
# between two regions, compute the pairwise correlation of each cell type
corr_mat_INS_PFC <- sapply(unique(DE_INS_PFC$cluster_id),function(cp){
  cor.test(DE_INS_PFC[DE_INS_PFC$cluster_id == cp,"logFC.x"],DE_INS_PFC[DE_INS_PFC$cluster_id == cp,"logFC.y"])$estimate
}
)
INS_PFC_cps <- unique(DE_INS_PFC$cluster_id)
names(corr_mat_INS_PFC) <- INS_PFC_cps


corr_mat_INS_VST <- sapply(unique(DE_INS_VST$cluster_id),function(cp){
  cor.test(DE_INS_VST[DE_INS_VST$cluster_id == cp,"logFC.x"],DE_INS_VST[DE_INS_VST$cluster_id == cp,"logFC.y"])$estimate
}
)
INS_VST_cps <- unique(DE_INS_VST$cluster_id)
names(corr_mat_INS_VST) <- INS_VST_cps


corr_mat_PFC_VST <- sapply(unique(DE_PFC_VST$cluster_id),function(cp){
  cor.test(DE_PFC_VST[DE_PFC_VST$cluster_id == cp,"logFC.x"],DE_PFC_VST[DE_PFC_VST$cluster_id == cp,"logFC.y"])$estimate
}
)
PFC_VST_cps <- unique(DE_PFC_VST$cluster_id)
names(corr_mat_PFC_VST) <- PFC_VST_cps

combined_cell_types <- unique(c(DE_INS_PFC$cluster_id,DE_INS_VST$cluster_id,DE_PFC_VST$cluster_id))
corr_mat_df <- Matrix::Matrix(data = NA,nrow=3,ncol = length(combined_cell_types))
colnames(corr_mat_df) <- combined_cell_types
for(cp in INS_PFC_cps){
  corr_mat_df[1,cp] <- corr_mat_INS_PFC[cp]
}
for(cp in INS_VST_cps){
  corr_mat_df[2,cp] <- corr_mat_INS_VST[cp]
}
for(cp in PFC_VST_cps){
  corr_mat_df[3,cp] <- corr_mat_PFC_VST[cp]
}
corr_mat_df <- as.data.frame(corr_mat_df)
corr_mat_df$regions <- c("INS-PFC","INS-VST","PFC-VST")
corr_mat_df.long <- corr_mat_df %>%
  pivot_longer(!regions, names_to = "cell_type", values_to = "correlation")
cp_orders <- c("Astro","Micro","Oligo","OPC","Endo","Pvalb","Sst","Vip","Other immune","VLMC",
               "L2/3 IT","L4 IT","L5 IT","L5/6 NP","Lamp5","Lamp5 Lhx6",
               "Chandelier"
)
corr_mat_df.long$regions <- factor(corr_mat_df.long$regions,levels = c("INS-PFC","INS-VST","PFC-VST"))
corr_mat_df.long$cell_type[corr_mat_df.long$cell_type == "Immune"] <- "Other immune"
corr_mat_df.long$cell_type <- factor(corr_mat_df.long$cell_type,levels = cp_orders)


corr_mat_df.long$cell_type_coarse <- cell_types_coarse_levels[as.character(corr_mat_df.long$cell_type)]
corr_mat_df.long$cell_type_coarse <- factor(corr_mat_df.long$cell_type_coarse,levels = c("Ext","Inh","Glial"))
p0 <- ggplot(corr_mat_df.long, aes(x = cell_type, y = regions, fill = correlation)) +
  geom_tile()+
  geom_text(data=corr_mat_df.long,aes(y=regions,label=round(correlation,digits = 2)),size=5)+
  scale_x_discrete(position = "bottom")+
  scale_fill_continuous(low = "white", high = "red",na.value = "white")+
  theme_void()+
  coord_flip()+
  theme(
    axis.text=element_text(size=20),
    #axis.text.x=element_text(hjust = 0),
    axis.text.y=element_text(size=20),
    axis.title=element_text(size=20,face="bold"),
    legend.title=element_text(size=20,face = "bold"), 
    legend.text=element_text(size=20)
  )+
  ylab("")+
  xlab("")+#theme(legend.position='bottom')+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

p0_colorbar <- ggplot(corr_mat_df.long)+
  geom_bar(mapping = aes(x = cell_type, y = 1, fill = cell_type_coarse), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values = cp_col2)+
  theme_void()+
  theme(
    axis.text=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=20),
    axis.title=element_text(size=20,face="bold"),
    legend.title=element_text(size=20,face = "bold"), 
    legend.text=element_text(size=20)
  )+
  ylab("")+
  xlab("")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme(legend.position='bottom')+
  guides(fill = guide_legend(ncol = 3,title = ""))+ coord_flip()


png(file="./figures/corr_logfc_heatmap.png",width = 1800,height = 3000,res=300)
print(p0+NoLegend())
dev.off()

png(file="./figures/corr_logfc_heatmap_colorbar_legend.png",width = 1000,height = 500,res=300)
print(as_ggplot(ggpubr::get_legend(p0_colorbar)))
dev.off()

png(file="./figures/corr_logfc_heatmap_corr_legend.png",width = 1000,height = 2000,res=300)
print(as_ggplot(get_legend(p0)))
dev.off()




# only label the most top ones


DE_INS_PFC <- DE_INS_PFC %>%
  mutate(DE_INS = if_else((abs(logFC.x) > 0.25)  &(p_adj.glb.x < 0.05),TRUE,FALSE)) %>%
  mutate(DE_PFC = if_else((abs(logFC.y) > 0.25)  &(p_adj.glb.y < 0.05),TRUE,FALSE)) %>%
  mutate(nonDE_INS = if_else((abs(logFC.x) < 0.25)  &(p_adj.glb.x > 0.05),TRUE,FALSE)) %>%
  mutate(nonDE_PFC = if_else((abs(logFC.y) < 0.25)  &(p_adj.glb.y > 0.05),TRUE,FALSE)) %>%
  mutate(DE = case_when((DE_INS&DE_PFC)~"Both DE", 
                        (DE_INS&(!DE_PFC))~"INS DE only",
                        ((!DE_INS)& DE_PFC)~"PFC DE only",
                        ((!DE_INS)&(!DE_PFC))~"None DE")) %>%
  mutate(isDE = if_else(DE!="None DE",TRUE,FALSE)) %>%
  mutate(top_DE_INS = ifelse((abs(logFC.x) > 1) &(p_adj.glb.x < 0.05),TRUE,FALSE)) %>%
  mutate(top_DE_PFC = ifelse((abs(logFC.y) > 1) &(p_adj.glb.y < 0.05),TRUE,FALSE)) %>%
  mutate(top_DE = top_DE_INS&top_DE_PFC) %>%
  mutate(top_DE_INS_only = top_DE_INS&nonDE_PFC) %>%
  mutate(top_DE_PFC_only = top_DE_PFC&nonDE_INS) %>%
  mutate(delabel = ifelse(top_DE|top_DE_INS_only|top_DE_PFC_only,gene,NA))

DE_INS_VST <- DE_INS_VST %>%
  mutate(DE_INS = if_else((abs(logFC.x) > 0.25)  &(p_adj.glb.x < 0.05),TRUE,FALSE)) %>%
  mutate(DE_VST = if_else((abs(logFC.y) > 0.25)  &(p_adj.glb.y < 0.05),TRUE,FALSE)) %>%
  mutate(nonDE_INS = if_else((abs(logFC.x) < 0.25)  &(p_adj.glb.x > 0.05),TRUE,FALSE)) %>%
  mutate(nonDE_VST = if_else((abs(logFC.y) < 0.25)  &(p_adj.glb.y > 0.05),TRUE,FALSE)) %>%
  mutate(isDE = if_else((DE_INS)|(DE_VST),TRUE,FALSE)) %>%
  mutate(top_DE_INS = ifelse((abs(logFC.x) > 1)&(p_adj.glb.x < 0.05),TRUE,FALSE)) %>%
  mutate(top_DE_VST = ifelse((abs(logFC.y) > 1)&(p_adj.glb.y < 0.05),TRUE,FALSE)) %>%
  mutate(top_DE = top_DE_INS&top_DE_VST) %>%
  mutate(top_DE_INS_only = top_DE_INS&nonDE_VST) %>%
  mutate(top_DE_VST_only = top_DE_VST&nonDE_INS) %>%
  mutate(delabel = ifelse(top_DE|top_DE_INS_only|top_DE_VST_only,gene,NA))



DE_PFC_VST <- DE_PFC_VST %>%
  mutate(DE_PFC = if_else((abs(logFC.x) > 0.25)  &(p_adj.glb.x < 0.05),TRUE,FALSE)) %>%
  mutate(DE_VST = if_else((abs(logFC.y) > 0.25)  &(p_adj.glb.y < 0.05),TRUE,FALSE)) %>%
  mutate(nonDE_PFC = if_else((abs(logFC.x) < 0.25)  &(p_adj.glb.x > 0.05),TRUE,FALSE)) %>%
  mutate(nonDE_VST = if_else((abs(logFC.y) < 0.25)  &(p_adj.glb.y > 0.05),TRUE,FALSE)) %>%
  mutate(isDE = if_else((DE_PFC)|(DE_VST),TRUE,FALSE)) %>%
  mutate(top_DE_PFC = ifelse((abs(logFC.x) > 1) &(p_adj.glb.x < 0.05),TRUE,FALSE)) %>%
  mutate(top_DE_VST = ifelse((abs(logFC.y) > 1) &(p_adj.glb.y < 0.05),TRUE,FALSE)) %>%
  mutate(top_DE = top_DE_PFC&top_DE_VST) %>%
  mutate(top_DE_PFC_only = top_DE_PFC&nonDE_VST) %>%
  mutate(top_DE_VST_only = top_DE_VST&nonDE_PFC) %>%
  mutate(delabel = ifelse(top_DE|top_DE_PFC_only|top_DE_VST_only,gene,NA)) 




p1 <- ggplot(DE_INS_PFC,aes(x = logFC.x,y=logFC.y,col=cluster_id,
                            alpha=isDE,label=delabel))+
  geom_point(size=2)+
  scale_alpha_discrete(range = c(0.01, 0.9))+
  scale_color_manual(values = cp_col)+
  geom_label_repel( max.overlaps = 40,size=5,show.legend = FALSE,aes(fill=cluster_id),color="black",fontface = "bold") +
  scale_fill_manual(values = cp_col)+
  #scale_x_continuous(alpha=0.1)+
  geom_vline(xintercept=0,linetype="dotted",alpha=0.5) +
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  xlab("logFC - INS")+
  ylab("logFC - PFC")+
  #theme_cowplot(font_size = 5) +
  theme_minimal()+
  theme(
    #axis.line = element_line(linetype = "dashed"),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20,face="bold"),
    legend.title=element_text(size=20,face = "bold"), 
    legend.text=element_text(size=20))+
  
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +guides(color=guide_legend(ncol=2))

png(file="./regional/figures/corr_INS_PFC_logfc_labeled.png",width = 3500,height = 3500,res=300)
print(p1+NoLegend())
dev.off()





p2 <- ggplot(DE_INS_VST,aes(x = logFC.x,y=logFC.y,col=cluster_id,
                            alpha=isDE,label=delabel))+
  geom_point(size=2)+
  scale_alpha_discrete(range = c(0.01, 0.9))+
  scale_color_manual(values = cp_col)+
  geom_label_repel( max.overlaps = 20,size=5,show.legend = FALSE,aes(fill=cluster_id),color="black",fontface = "bold") +
  scale_fill_manual(values = cp_col)+
  #scale_x_continuous(alpha=0.1)+
  geom_vline(xintercept=0,linetype="dotted",alpha=0.5) +
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  xlab("logFC - INS")+
  ylab("logFC - VST")+
  ylim(-7.5,7.5)+
  #theme_cowplot(font_size = 5) +
  theme_minimal()+
  theme(
    #axis.line = element_line(linetype = "dashed"),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20,face="bold"),
    legend.title=element_text(size=20,face = "bold"), 
    legend.text=element_text(size=20))+
  
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +guides(color=guide_legend(ncol=2))

png(file="./regional/figures/corr_INS_VST_logfc_labeled.png",width = 3500,height = 3500,res=300)
print(p2+NoLegend())
dev.off()




p3 <- ggplot(DE_PFC_VST,aes(x = logFC.x,y=logFC.y,col=cluster_id,
                            alpha=isDE,label=delabel))+
  geom_point(size=2)+
  scale_alpha_discrete(range = c(0.01, 0.9))+
  scale_color_manual(values = cp_col)+
  geom_label_repel( max.overlaps = 40,size=5,show.legend = FALSE,aes(fill=cluster_id),color="black",fontface = "bold") +
  scale_fill_manual(values = cp_col)+
  #scale_x_continuous(alpha=0.1)+
  geom_vline(xintercept=0,linetype="dotted",alpha=0.5) +
  geom_hline(yintercept = 0,linetype="dotted",alpha=0.5)+
  xlab("logFC - PFC")+
  ylab("logFC - VST")+
  ylim(-7.5,7.5)+
  #theme_cowplot(font_size = 5) +
  theme_minimal()+
  theme(
    #axis.line = element_line(linetype = "dashed"),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20,face="bold"),
    legend.title=element_text(size=20,face = "bold"), 
    legend.text=element_text(size=20))+
  
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +guides(color=guide_legend(ncol=2))

png(file="./regional/figures/corr_PFC_VST_logfc_labeled.png",width = 3500,height = 3500,res=300)
print(p3+NoLegend())
dev.off()

