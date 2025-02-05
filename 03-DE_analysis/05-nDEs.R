# plot the number of DEs
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/03-DE_analysis/"
setwd(wd)

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


de_edgeR <- rbind(rbind(PFC_all,INS_all),VST_all)




de_edgeR.top <- de_edgeR %>%
  dplyr::filter((abs(logFC) > 0.25) &(p_adj.glb < 0.05)) %>%
  mutate(direction = ifelse(logFC >0,"HIV UP","HIV DOWN")) %>%
  group_by(region,cluster_id,direction) %>%
  summarise(n_DEs = n())

region_cols <- c("#E69F00", "#56B4E9","#CC79A7")
names(region_cols) <-c("INS","PFC","VST")

condition_cols <- c("#0072B2","darkorange1")
names(condition_cols) <- c("HIV DOWN","HIV UP")

cell_types_consensus_levels <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                 "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Immune","Macrophages","PC","Endo","SMC","VLMC")
cp_col <- c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1",
  "purple1", "lightslateblue","orchid3","pink3","plum1",  "tomato1", "tomato3","sienna2", "burlywood3",
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","olivedrab2","olivedrab4","darkorange3","goldenrod","yellow4", "yellow3"
  
)
names(cp_col) <- cell_types_consensus_levels
cp_col2 <-  c("tomato2","royalblue3","dodgerblue","gold2")
names(cp_col2) <- c("Ext","MSNs","Inh","Glial")

cell_types_coarse_levels <- c(rep("Ext",9),rep("MSNs",9),
                                        rep("Inh",9),rep("Glial",10))
names(cell_types_coarse_levels) <- cell_types_consensus_levels

de_edgeR.top$cluster_id_coarse <- factor(cell_types_coarse_levels[de_edgeR.top$cluster_id],
                                         levels = c("Ext","MSNs","Inh","Glial"))

de_edgeR.top$cluster_id <- factor(de_edgeR.top$cluster_id,levels = cell_types_consensus_levels)
de_edgeR.top$direction <- factor(de_edgeR.top$direction,levels = c("HIV UP","HIV DOWN"))

de_edgeR.top.up <- de_edgeR.top %>%
  dplyr::filter(direction == "HIV UP")

de_edgeR.top.down <- de_edgeR.top %>%
  dplyr::filter(direction == "HIV DOWN")

# plot as heatmaps
for (re in c("PFC","INS","VST")) {
  p0 <- ggplot(de_edgeR.top.up[de_edgeR.top.up$region == re,], aes(x = cluster_id, y = region, fill = n_DEs)) +
    #facet_wrap(~region, scales='free',nrow=3)+
    geom_tile()+
    geom_text(data=de_edgeR.top.up[de_edgeR.top.up$region == re,],aes(y=region,label=n_DEs),size=10)+
    scale_fill_continuous(low = "white", high = "red",na.value = "white")+
    theme_void()+
    theme(
      axis.text=element_text(size=25),
      axis.text.x=element_text(angle = 60, vjust = 0.5),
      axis.text.y=element_blank(),
      axis.title=element_text(size=25,face="bold"),
      legend.title=element_text(size=25,face = "bold"), 
      legend.text=element_text(size=25)
    )+
    ylab("")+
    xlab("")+#theme(legend.position='bottom')+
    theme(plot.margin = unit(c(0,0,0,0), "cm"))
  
  p0_colorbar_up <- ggplot(de_edgeR.top.up[de_edgeR.top.up$region == re,])+
    geom_bar(mapping = aes(x = cluster_id, y = 1, fill = cluster_id), 
             stat = "identity", 
             width = 1)+
    scale_fill_manual(values = cp_col)+
    theme_void()+
    theme(
      axis.text=element_text(size=20),
      axis.text.y=element_blank(),
      axis.text.x=element_blank(),
      axis.title=element_text(size=20,face="bold"),
      legend.title=element_text(size=20,face = "bold"), 
      legend.text=element_text(size=20)
    )+
    ylab("")+
    xlab("")+
    theme(plot.margin = unit(c(0,0,0,0), "cm"))+
    theme(legend.position='bottom')+
    guides(fill = guide_legend(ncol = 2))
  
  p0_colorbar_up2 <- ggplot(de_edgeR.top.up[de_edgeR.top.up$region == re,])+
    geom_bar(mapping = aes(x = cluster_id, y = 1, fill = cluster_id_coarse), 
             stat = "identity", 
             width = 1)+
    scale_fill_manual(values = cp_col2)+
    theme_void()+
    theme(
      axis.text=element_text(size=25),
      axis.text.y=element_blank(),
      axis.text.x=element_blank(),
      axis.title=element_text(size=25,face="bold"),
      legend.title=element_text(size=25,face = "bold"), 
      legend.text=element_text(size=25)
    )+
    ylab("")+
    xlab("")+
    theme(plot.margin = unit(c(0,0,0,0), "cm"))+
    theme(legend.position='bottom')+
    guides(fill = guide_legend(ncol = 2))
  
  #p0_all <- plot_grid(p0+NoLegend(), NULL,p0_colorbar_up+NoLegend(),NULL,
  #                    p0_colorbar_up2+NoLegend(), align = "v",
  #                         nrow  = 5, rel_heights = c(2,-0.4,1,-0.4,1))
  
  p0_all <- plot_grid(p0_colorbar_up2+NoLegend(), NULL,p0+NoLegend(), 
                      align = "v",
                      nrow  = 3, rel_heights = c(1,-0.3,3))
  
  png(filename = paste0("./plotting/figures/heatmap_n_DE_",re,"_up.png"),width=6000,height = 1000,res = 300)
  print(p0_all+NoLegend())
  dev.off()
  
  png(filename = paste0("./plotting/figures/heatmap_n_DE_",re,"_up_legend.png"),width=1000,height = 500,res = 300)
  print(as_ggplot(ggpubr::get_legend(p0)))
  dev.off()
  
  p1 <- ggplot(de_edgeR.top.down[de_edgeR.top.down$region == re,], aes(x = cluster_id, y = region, fill = n_DEs)) +
    #facet_wrap(~region, scales='free',nrow=3)+
    geom_tile()+
    geom_text(data=de_edgeR.top.down[de_edgeR.top.down$region == re,],aes(y=region,label=n_DEs),size=10)+
    scale_fill_continuous(low = "white", high = "blue",na.value = "white")+
    theme_void()+
    theme(
      axis.text=element_text(size=25),
      axis.text.x=element_text(angle = 60, vjust = 0.5),
      axis.text.y=element_blank(),
      axis.title=element_text(size=25,face="bold"),
      legend.title=element_text(size=25,face = "bold"), 
      legend.text=element_text(size=25)
    )+
    ylab("")+
    xlab("")+#theme(legend.position='bottom')+
    theme(plot.margin = unit(c(0,0,0,0), "cm"))
  
  p1_colorbar_down <- ggplot(de_edgeR.top.down[de_edgeR.top.down$region == re,])+
    geom_bar(mapping = aes(x = cluster_id, y = 1, fill = cluster_id), 
             stat = "identity", 
             width = 1)+
    scale_fill_manual(values = cp_col)+
    theme_void()+
    theme(
      axis.text=element_text(size=20),
      axis.text.y=element_blank(),
      axis.text.x=element_blank(),
      axis.title=element_text(size=20,face="bold"),
      legend.title=element_text(size=20,face = "bold"), 
      legend.text=element_text(size=20)
    )+
    ylab("")+
    xlab("")+
    theme(plot.margin = unit(c(0,0,0,0), "cm"))+
    theme(legend.position='bottom')+
    guides(fill = guide_legend(ncol = 2))
  
  
  p1_colorbar_down2 <- ggplot(de_edgeR.top.down[de_edgeR.top.down$region == re,])+
    geom_bar(mapping = aes(x = cluster_id, y = 1, fill = cluster_id_coarse), 
             stat = "identity", 
             width = 1)+
    scale_fill_manual(values = cp_col2)+
    theme_void()+
    theme(
      axis.text=element_text(size=20),
      axis.text.y=element_blank(),
      axis.text.x=element_blank(),
      axis.title=element_text(size=20,face="bold"),
      legend.title=element_text(size=20,face = "bold"), 
      legend.text=element_text(size=20)
    )+
    ylab("")+
    xlab("")+
    theme(plot.margin = unit(c(0,0,0,0), "cm"))+
    theme(legend.position='bottom')+
    guides(fill = guide_legend(ncol = 2))
  
  #p1_down_all <- plot_grid(p1+NoLegend(), NULL,p1_colorbar_down+NoLegend(),NULL,
  #                         p1_colorbar_down2+NoLegend(), align = "v",
  #                    nrow  = 5, rel_heights = c(2,-0.4,1,-0.4,1))
  
  p1_down_all <- plot_grid(p1_colorbar_down2+NoLegend(), NULL,p1+NoLegend(), 
                      align = "v",
                      nrow  = 3, rel_heights = c(1,-0.3,3))
  
  png(filename = paste0("./plotting/figures/heatmap_n_DE_",re,"_down.png"),width=6000,height = 1000,res = 300)
  print(p1_down_all)
  dev.off()
  
  png(filename = paste0("./plotting/figures/heatmap_n_DE_",re,"_down_legend.png"),width=1000,height = 500,res = 300)
  print(as_ggplot(ggpubr::get_legend(p1)))
  dev.off()
}


coarse_label_legend <- as_ggplot(ggpubr::get_legend(ggplot(de_edgeR.top.up)+
  geom_bar(mapping = aes(x = cluster_id, y = 1, fill = cluster_id_coarse), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values = cp_col2)+
  theme_void()+
  theme(
    axis.text=element_text(size=25),
    axis.text.y=element_blank(),
    axis.text.x=element_blank(),
    axis.title=element_text(size=25,face="bold"),
    legend.title=element_text(size=25,face = "bold"), 
    legend.text=element_text(size=25)
  )+
  ylab("")+
  xlab("")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  theme(legend.position='bottom')+
  guides(fill = guide_legend(title="",ncol = 1))))

png(filename = paste0("./plotting/figures/heatmap_n_DE_coarse_cp_legend.png"),width=1000,height = 1000,res = 300)
print(coarse_label_legend)
dev.off()
