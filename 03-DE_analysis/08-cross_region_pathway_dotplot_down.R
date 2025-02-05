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
library(GOSim)
wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/03-DE_analysis/cross_region//"
setwd(wd)

de_edgeR <- readRDS("./data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")


top_edgeR_all.down <- de_edgeR  %>%
  dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05))%>%
  dplyr::select(gene,cluster_id,logFC,logCPM,p_adj.glb) %>%
  dplyr::filter((logFC < -0.25) &(p_adj.glb < 0.05))
top_edgeR_all.down$method <- "edgeR"

cell_types <- unique(top_edgeR_all.down$cluster_id)
top_edgeR_all.down$cell_type <- top_edgeR_all.down$cluster_id


top_edgeR_list <- lapply(cell_types,function(cp) top_edgeR_all.down[top_edgeR_all.down$cluster_id == cp,])
names(top_edgeR_list) <- cell_types



de_up_list <- list()
for(cp in cell_types){
  top_edgeR <- top_edgeR_all.down %>%
    dplyr::filter((cell_type == cp))
  if(nrow(top_edgeR) > 1){
    # go analysis
    de <- AnnotationDbi::select(Homo.sapiens, keys = unique(top_edgeR$gene), 
                                keytype = "SYMBOL", columns="ENTREZID")
    degs_converted_top_down <- top_edgeR %>%
      mutate(SYMBOL=gene) %>%
      left_join(de,by="SYMBOL") %>%
      dplyr::filter(!is.na(ENTREZID))
    
    de_up_list[[cp]] <- degs_converted_top_down$ENTREZID
    
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

#saveRDS(ego,"./data/all_ego_down.RDS")
ego_df <- as.data.frame(ego@compareClusterResult) %>%
  mutate(cell_type = Cluster)

ego_df.sig <- ego_df %>%
  group_by(cell_type) %>%
  slice_min(pvalue,n=5)

ego_df.list <- lapply(cell_types,function(cp) ego_df[ego_df$cell_type == cp,])
names(ego_df.list) <- cell_types

ego_df.sig.list <- lapply(cell_types,function(cp) ego_df.sig[ego_df.sig$cell_type == cp,])
names(ego_df.sig.list) <- cell_types

ego_df.sig.df.expanded <- ego_df[ego_df$ID %in% ego_df.sig$ID,]

ego_df.long <- ego_df.sig.df.expanded %>%
  separate_longer_delim(geneID, delim ='/') %>% # split the genes
  mutate(gene = geneID,pathway = Description) %>%
  left_join(top_edgeR_all.down,by = c("cell_type","gene"))%>%
  group_by(cell_type,pathway) %>%
  mutate(avr_logFC = mean(logFC)) %>%
  ungroup() %>%
  dplyr::select(cell_type,pathway,p.adjust,avr_logFC) %>%
  distinct() %>%
  mutate(pathway_short = substr(pathway,1,55)) 

ego_df.wide <- ego_df.long %>%
  dplyr::select(cell_type,pathway_short,avr_logFC) %>%
  distinct() %>% # this is to avoid a bug of different pathways shortened to the same names 
  pivot_wider(names_from = pathway_short, values_from = avr_logFC)

#colnames(ego_df.wide) <- sapply(colnames(ego_df.wide),function(cname) substr(cname,1,50))

ego_df.wide.mat <- ego_df.wide[,2:ncol(ego_df.wide)]
ego_df.wide.mat[is.na(ego_df.wide.mat)] <- 0
#rownames(ego_df.wide.mat) <- ego_df.wide[,1]

#biclust_order <- biclust::biclust(ego_df.wide.mat,method=biclust::BCCC())
row.order <- hclust(dist(ego_df.wide.mat))$order
col.order <- hclust(dist(t(ego_df.wide.mat)))$order

#pathway_cluster <- kmeans(t(ego_df.wide.mat[,2:ncol(ego_df.wide.mat)]),10)$cluster

ego_df.long$cell_type <- factor(ego_df.long$cell_type,levels = ego_df.wide$cell_type[row.order])
#ego_df.long$pathway_short <- sapply(ego_df.long$pathway_short,function(cname) substr(cname,1,50))
ego_df.long$pathway_short <- factor(ego_df.long$pathway_short,levels = colnames(ego_df.wide.mat)[col.order])

g1 <- ggplot(ego_df.long, aes(x=cell_type, y=pathway_short,size = -log(p.adjust))) +
  
  geom_point(aes(color = abs(avr_logFC)),alpha=1)+
  scale_size(range=c(0.5,10))+
  scale_color_continuous(low = "lightgrey", high = "blue",n.breaks=5,limit=c(0,3))+
  theme_void()+
  theme(
    axis.text=element_text(size=20),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text.y=element_text(hjust = 0),
    axis.title=element_text(size=20,face="bold"),
    legend.title=element_text(size=20,face = "bold"), 
    legend.text=element_text(size=15),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
  scale_y_discrete(position = "right")+
  ylab("")+
  xlab("")+
  theme(legend.position=c(1.3,1),legend.direction = "horizontal",
        plot.margin = unit(c(0,0,0,0), "cm"))+
  guides(fill = guide_legend(order = 1), 
         size = guide_legend(order = 2))

g2 <- as_ggplot(ggpubr::get_legend(g1))
png(filename = paste0("./figures/cross_region_pathway_dotplot_down.png"),width=4000,height = 4000,res = 300)
print((g1+NoLegend()) / g2 +   plot_layout(heights = c(12, 1)))
dev.off()


