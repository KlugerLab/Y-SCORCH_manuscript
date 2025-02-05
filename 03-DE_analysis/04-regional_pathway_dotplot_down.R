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
library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(scater)
library(tidyr)
library(dplyr)
library(cowplot)



setwd("/banach2/SCORCH/scripts/combined_analysis_paper/03-DE_analysis/")


INS_edgeR <- readRDS("./INS/results/INS_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
INS_edgeR <- INS_edgeR %>%
  ungroup() %>%
  dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05))%>%
  dplyr::select(gene,cluster_id,logFC,logCPM,p_adj.glb,region) 

PFC_edgeR <- readRDS("./PFC//results/PFC_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
PFC_edgeR <- PFC_edgeR %>%
  ungroup() %>%
  dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05))%>%
  dplyr::select(gene,cluster_id,logFC,logCPM,p_adj.glb,region) 

VST_edgeR <- readRDS("./VST//results/VST_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
VST_edgeR <- VST_edgeR %>%
  ungroup() %>%
  dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05))%>%
  dplyr::select(gene,cluster_id,logFC,logCPM,p_adj.glb,region) 

DE_all <- rbind(rbind(PFC_edgeR,INS_edgeR),VST_edgeR)


top_edgeR_all.down <- DE_all  %>%
  dplyr::filter((logFC < -0.25) &(p_adj.glb < 0.05)) #&(logCPM > 0)

top_edgeR_all.down$cell_type <- top_edgeR_all.down$cluster_id
cell_types <- unique(top_edgeR_all.down$cluster_id)


de_up_list <- list()
for(reg in c("PFC","INS","VST")){
  for(cp in cell_types){
    top_edgeR <- top_edgeR_all.down %>%
      dplyr::filter((region == reg) & (cell_type == cp))
    if(nrow(top_edgeR) > 1){
      # go analysis
      de <- AnnotationDbi::select(Homo.sapiens, keys = unique(top_edgeR$gene), 
                                  keytype = "SYMBOL", columns="ENTREZID")
      degs_converted_top_down <- top_edgeR %>%
        mutate(SYMBOL=gene) %>%
        left_join(de,by="SYMBOL") %>%
        dplyr::filter(!is.na(ENTREZID))
      
      de_up_list[[paste0(reg,'-',cp)]] <- degs_converted_top_down$ENTREZID
      
    }
    
  }
}


ego <- compareCluster(geneClusters = de_up_list,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff  = 1,
                      keyType = "ENTREZID",
                      OrgDb = org.Hs.eg.db,
                      readable = TRUE,
                      ont = "BP",
                      fun = "enrichGO")

ego_df <- as.data.frame(ego@compareClusterResult) 
ego_df <- ego_df %>% extract(Cluster, into = c("region","cell_type"), "^(^[A-Z]*)\\-(.*)")

ego_df.sig <- ego_df %>%
  group_by(cell_type,region) %>%
  dplyr::filter(Count > 2) %>% # have more than 1 gene enriched
  slice_min(pvalue,n=5)

ego_df.list <- lapply(cell_types,function(cp) ego_df[ego_df$cell_type == cp,])
names(ego_df.list) <- cell_types

ego_df.sig.list <- lapply(cell_types,function(cp) ego_df.sig[ego_df.sig$cell_type == cp,])
names(ego_df.sig.list) <- cell_types

ego_df.sig.list.expanded <- list()
for(cell_type in cell_types){
  ego_df.sig.list.expanded[[cell_type]] <- ego_df.list[[cell_type]][ego_df.list[[cell_type]]$ID %in% ego_df.sig.list[[cell_type]]$ID,]
}

ego_df.sig.df.expanded <- do.call(rbind,ego_df.sig.list.expanded)

ego_df.long <- ego_df.sig.df.expanded %>%
  separate_longer_delim(geneID, delim ='/') %>% # split the genes
  mutate(gene = geneID,pathway = Description) %>%
  left_join(top_edgeR_all.down,by = c("region","cell_type","gene"))%>%
  group_by(region,cell_type,pathway) %>%
  mutate(avr_logFC = mean(logFC)) %>%
  ungroup() %>%
  dplyr::select(region,cell_type,pathway,p.adjust,avr_logFC) %>%
  distinct() %>%
  mutate(pathway_short = substr(pathway,1,50)) 


ego_df.long$region <- factor(ego_df.long$region,levels = c("PFC","INS","VST"))
neuron_types <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                  "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                  "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip")

for (cp in unique(ego_df.long$cell_type)) {
  margin_x <- ifelse(cp %in% neuron_types, 2.3,2)
  margin_x <- ifelse(cp == "Astro",2.5,margin_x)
  margin_x <- ifelse(cp == "Micro",1.8,margin_x)
  margin_x <- ifelse(cp == "Endo",2.1,margin_x)
  margin_x <- ifelse(cp == "OPC",2.1,margin_x)
  
  g1 <- ggplot(ego_df.long[ego_df.long$cell_type == cp,], aes(x=region, y=pathway_short,size = -log(p.adjust))) +
    
    geom_point(aes(color = abs(avr_logFC)),alpha=1)+
    scale_size(range=c(1,10))+
    scale_color_continuous(low = "lightgrey", high = "blue",n.breaks=3,limit=c(0,3))+
    theme_void()+
    theme(
      axis.text=element_text(size=20),
      axis.text.y=element_text(hjust = 0),
      axis.title=element_text(size=20,face="bold"),
      legend.title=element_text(size=20,face = "bold"), 
      legend.text=element_text(size=15),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
    scale_y_discrete(position = "right")+
    ylab("")+
    xlab("")+
    theme(legend.position=c(margin_x,1),legend.direction = "horizontal",
          plot.margin = unit(c(0,0,0,0), "cm"))+
    guides(fill = guide_legend(order = 1), 
           size = guide_legend(order = 2))
  
  g2 <- as_ggplot(ggpubr::get_legend(g1))
  png(filename = paste0("./regional/figures_tmp//dotplot_pathway_heatmap_down/",make.names(cp),".png"),width=2500,height = 2000,res = 300)
  print((g1+NoLegend()) / g2 +   plot_layout(heights = c(12, 1)))
  dev.off()
  
  png(filename = paste0("./regional/figures_tmp/dotplot_pathway_heatmap_down/",make.names(cp),"_legend_.png"),width=3000,height = 500,res = 300)
  print(as_ggplot(ggpubr::get_legend(g1)))
  dev.off()
}

