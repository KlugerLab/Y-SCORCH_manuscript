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

top_edgeR_all <- de_edgeR  %>%
  dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05))%>%
  dplyr::select(gene,cluster_id,logFC,logCPM,p_adj.glb) %>%
  dplyr::filter((logFC > 0.25)&(p_adj.glb < 0.05))
top_edgeR_all$method <- "edgeR"

top_edgeR_all.up.microglia <- top_edgeR_all %>%
  dplyr::filter(cluster_id == "Micro")

cell_types <- unique(top_edgeR_all$cluster_id)
top_edgeR_list <- lapply(cell_types,function(cp) top_edgeR_all[top_edgeR_all$cluster_id == cp,])
names(top_edgeR_list) <- cell_types




cp <- "Micro"
top_edgeR <- top_edgeR_list[[cp]]

# go analysis
de <- AnnotationDbi::select(Homo.sapiens, keys = unique(top_edgeR$gene), keytype = "SYMBOL", columns="ENTREZID")
degs_converted_top_up <- top_edgeR %>%
  mutate(SYMBOL=gene) %>%
  left_join(de,by="SYMBOL") %>%
  dplyr::filter(!is.na(ENTREZID))

# 
gene_list <- degs_converted_top_up$logFC
names(gene_list) <- degs_converted_top_up$ENTREZID

ego <- enrichGO(gene          = degs_converted_top_up$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)


edox <- setReadable(ego, 'org.Hs.eg.db', 'ENTREZID')
ego_df <- ego@result
new_gene_list <- as.numeric(degs_converted_top_up$logFC)
names(new_gene_list) <- degs_converted_top_up$gene
new_gene_list = sort(new_gene_list, decreasing = TRUE)

p1 <- cnetplot(edox, color.params = list(foldChange = new_gene_list),
               cex_label_category=1.2,cex_gene=3,cex_category=1.5,node_label="gene",
               showCategory=5)

p1 <- p1+ scale_color_gradientn(colours = c("white", "red"))+NoLegend()

png(filename = paste0("./figures/microglia_pathways/go_netplot_micro_new.png"),width=5000,height = 5000,res = 600)
print(p1)
dev.off()

p1_withlabel <- cnetplot(edox, color.params = list(foldChange = new_gene_list),
                         cex_label_category=1.2,cex_gene=3,cex_category=1.5,
                         showCategory=5)+ scale_color_gradientn(colours = c("white", "red"))
png(filename = paste0("./figures/microglia_pathways/go_netplot_micro_withlabel_new.png"),width=5000,height = 5000,res = 300)
print(p1_withlabel)
dev.off()






# heatmap
top_edgeR_all <- de_edgeR  %>%
  dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05))%>%
  dplyr::select(gene,cluster_id,logFC,logCPM,p_adj.glb) %>%
  dplyr::filter((abs(logFC) > 0.25) &(p_adj.glb < 0.05))
top_edgeR_all$method <- "edgeR"



cell_types <- unique(top_edgeR_all$cluster_id)
top_edgeR_list <- lapply(cell_types,function(cp) top_edgeR_all[top_edgeR_all$cluster_id == cp,])
names(top_edgeR_list) <- cell_types

top_edgeR_target.up <- top_edgeR_list[[cp]] %>%
  dplyr::filter(logFC > 0) %>%
  slice_min(p_adj.glb ,n=25)

top_edgeR_target.down <- top_edgeR_list[[cp]] %>%
  dplyr::filter(logFC < 0) %>%
  slice_min(p_adj.glb ,n=25)

top_edgeR_target <- rbind(top_edgeR_target.down,top_edgeR_target.up)



coldf <- data.frame(
  region = pb@colData$region,
  condition = pb@colData$condition_new,
  sample_id = colnames(pb))
rownames(coldf) <- colnames(pb.norm_mean)
coldf <- coldf %>%
  arrange(condition,region,sample_id)

condition_cols <- c("#0072B2","darkorange1")
names(condition_cols) <- c("CTR","HIV")

region_cols <- c("#E69F00", "#56B4E9","#CC79A7")
names(region_cols) <-c("INS","PFC","VST")


ann_colors = list(
  condition = condition_cols,
  region = region_cols
)

pb.norm_mean.data.sub <- pb.norm_mean@assays@data[[cp]]
pb.norm_mean.data.sub <- pb.norm_mean.data.sub[,rownames(coldf)]



# only keep samples with enough microglia
n_cells <- do.call(rbind,pb@int_colData$n_cells)
n_cells <- n_cells[rownames(coldf),]
s2keep <- n_cells[,"Micro"] > 30

coldf <- coldf[s2keep,]

mat <- t(scale(t(pb.norm_mean.data.sub[top_edgeR_target$gene,s2keep])))
colors <- c(-4.1,seq(-4,4,by=0.001),4.1)
my_palette <- c("steelblue",colorRampPalette(colors = c("steelblue", "white", "red1"))
                (n = length(colors)-3), "red1")
col_fun = colorRamp2(c(-4, 0, 4), c("steelblue", "white", "red1"))
p1 <- ComplexHeatmap::Heatmap(mat[,rownames(coldf)], 
                              col=col_fun,
                              column_names_gp = gpar(fontsize = 25),
                              row_names_gp = gpar(fontsize = 30),
                              cluster_columns = F, cluster_rows = F,
                              show_row_dend = F,
                              name=" ",
                              column_gap = unit(0.1, "mm"),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(
                                "Condition"=coldf$condition,
                                "Region" = coldf$region,
                                
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

png(filename = paste0("./figures/microglia_pathways/heatmap.png"),width=6000,height = 6000,res = 300)
print(p1)
dev.off()

svg(filename = paste0("./figures/heatmaps/",make.names(cp),".svg"),width=20,height = 15)
print(p1)
dev.off()
