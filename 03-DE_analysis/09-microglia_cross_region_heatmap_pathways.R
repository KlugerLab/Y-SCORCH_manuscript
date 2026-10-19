# REVISION: run on the CellBender-corrected data; paths come from
# revision/config/paths.R and outputs go to OUT_ROOT. Changes are marked ## REVISION.

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
## in this script, so the import is disabled.
# library(GOSim)
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION
## REVISION: run against either DE result set.
##   SCORCH_DE_COUNTS = "uncorrected" (default, the primary analysis) | "corrected"
##   Inputs and outputs both carry the suffix, so the two never mix and the
##   corrected figures cannot overwrite the primary ones.
DE_COUNTS <- match.arg(Sys.getenv("SCORCH_DE_COUNTS", unset = "uncorrected"),
                       c("uncorrected", "corrected"))
OUT_SFX <- if (DE_COUNTS == "corrected") "_corrected" else ""
## every png()/pdf()/ggsave() in this script picks this up via revision/config/paths.R
options(scorch.out_suffix = paste0(OUT_SFX, XR_OUT_SFX))   ## REVISION: + "_donor_aware" when SCORCH_XR_DE=donor_aware
message("[fig] DE result set: ", DE_COUNTS)

## REVISION: working directory is set by rerun_setup() under OUT_ROOT instead of a hardcoded server path.
rerun_setup("03-DE_analysis/cross_region/",
            subdirs = c("data", "figures/heatmaps", "figures/microglia_pathways"))
de_edgeR <- readRDS(xr_de_file(OUT_SFX))   ## REVISION: edgeR table by default; the donor-aware (primary) table under SCORCH_XR_DE=donor_aware

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

## REVISION: current enrichplot renamed the cnetplot() arguments (color_category /
##   size_category / color_item / size_item / node_label / foldChange); the
##   published color.params=list(...), cex_gene and cex_category raise "unused
##   arguments". cex_label_category has no equivalent and is dropped; label size
##   is controlled by the theme.
p1 <- cnetplot(edox, foldChange = new_gene_list,
               size_item = 3, size_category = 1.5, node_label = "gene",
               showCategory=5)

p1 <- p1+ scale_color_gradientn(colours = c("white", "red"))+NoLegend()

png(filename = paste0("./figures/microglia_pathways/go_netplot_micro_new.png"),width=5000,height = 5000,res = 600)
print(p1)
dev.off()

p1_withlabel <- cnetplot(edox, foldChange = new_gene_list,
                         size_item = 3, size_category = 1.5,
                         showCategory=5)+ scale_color_gradientn(colours = c("white", "red"))
png(filename = paste0("./figures/microglia_pathways/go_netplot_micro_withlabel_new.png"),width=5000,height = 5000,res = 300)
print(p1_withlabel)
dev.off()






# heatmap
## REVISION: pb / pb.norm_mean were used below but never loaded in the published
##   script. Load the cross-region pseudobulks written by 07-cross_region_edgeR.R.
pb <- readRDS(paste0("./data/pseudo_bulk", OUT_SFX, ".rds"))
pb.norm_mean <- readRDS(paste0("./data/pseudo_norm_mean_bulk", OUT_SFX, ".rds"))
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
