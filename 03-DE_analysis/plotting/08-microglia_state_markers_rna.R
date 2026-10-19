# ==============================================================================
# REVISION: authors' plotting script, added to this repository with the revision; run on
# the CellBender-corrected data; paths come from revision/config/paths.R and
# outputs go to OUT_ROOT. Changes are marked ## REVISION.
#
# Figure S9a (gene expression): microglia state markers across samples, columns grouped by
# condition then region, rows the marker panels (homeostatic, anti-inflammatory,
# interferon-responsive, antigen presentation, proliferative, disease-associated).
# Samples with <= 50 microglia are dropped.
# ==============================================================================

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
##   in this script, so the import is disabled.
# library(GOSim)
## REVISION: working directory is set by rerun_setup() under OUT_ROOT instead of a hardcoded server path.
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
rerun_setup("03-DE_analysis/cross_region/",
            subdirs = c("data", "figures/microglia_state_markers"))
edgeR_res <- readRDS("./data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")  ## REVISION: only the _withFreq table is written
pb.norm_mean <- readRDS("./data/pseudo_norm_mean_bulk.rds")


coldf <- data.frame(
  region = pb.norm_mean@colData$region,
  condition = pb.norm_mean@colData$condition_new,
  sample_id = pb.norm_mean@colData$sample_name)

rownames(coldf) <- coldf$sample_id
## REVISION: condition is relabelled PWoH / PWH, as in every other revision figure.
##   The factor keeps PWoH first, so the column order matches the published CTR-then-HIV
##   layout (plain alphabetical sorting would put PWH first).
coldf$condition <- factor(as.character(pwh_relabel(coldf$condition)), levels = PWH_LEVELS)
coldf <- coldf %>%
  arrange(condition,region,sample_id)

condition_cols <- c("#0072B2","darkorange1")
names(condition_cols) <- PWH_LEVELS   ## REVISION: PWoH / PWH, same colours

region_cols <- c("#E69F00", "#56B4E9","#CC79A7")
names(region_cols) <-c("INS","PFC","VST")


ann_colors = list(
  condition = condition_cols,
  region = region_cols
)
target_genes <- c("P2RY12","CX3CR1","TMEM119","SALL1",
                  "IL4", "IL13", "IL10", "TGFB1", "TGFB3",
                  "IFITM3","IFIT3","ISG15","CD74","HLA-DRA",
                  "MKI67","PCNA","APOE", "LPL")

pb.norm_mean.sub <- pb.norm_mean@assays@data[["Micro"]]
colnames(pb.norm_mean.sub) <- pb.norm_mean@colData$sample_name
pb.norm_mean.sub <- pb.norm_mean.sub[,rownames(coldf)]

# only keep samples with enough microglia
n_cells <- do.call(rbind,pb.norm_mean@int_colData$n_cells)
rownames(n_cells) <- pb.norm_mean@colData$sample_name
n_cells <- n_cells[rownames(coldf),]
s2keep <- n_cells[,"Micro"] > 50



mat <- t(scale(t(pb.norm_mean.sub[target_genes,s2keep])))
coldf2keep <- coldf[s2keep,]
colors <- c(-4.1,seq(-4,4,by=0.001),4.1)
my_palette <- c("steelblue",colorRampPalette(colors = c("steelblue", "white", "red1"))
                (n = length(colors)-3), "red1")
col_fun = colorRamp2(c(-4, 0, 4), c("steelblue", "white", "red1"))
p1 <- ComplexHeatmap::Heatmap(mat, 
                              col=col_fun,
                              column_names_gp = gpar(fontsize = 25),
                              row_names_gp = gpar(fontsize = 30),
                              cluster_columns = F, cluster_rows = F,
                              show_row_dend = F,
                              name=" ",
                              column_gap = unit(0.1, "mm"),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(
                                "Condition"=coldf2keep$condition,
                                "Region" = coldf2keep$region,
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

png(filename = paste0("./figures/microglia_state_markers/microglia_markers_rna.png"),width=3500,height = 3500,res = 300)
print(p1)
dev.off()


# check the significance of these genes
## REVISION: the published script read all_HIV-CTR_res_edgeR_filtered.RDS (a muscat
##   result list, not written here) and indexed $table[[1]]$Micro. The same check
##   is done on the flat _withFreq data frame, and the answer is printed to the log.
de.micro.all <- readRDS("./data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS") %>%
  ungroup() %>% dplyr::filter(cluster_id == "Micro")
de.micro.all.sig <- de.micro.all %>%
  dplyr::filter((abs(logFC) > 0.25) & (logCPM > 0) & (p_adj.glb < 0.05))
message("[S9a] marker genes significant in cross-region microglia DE: ",
        paste(target_genes[target_genes %in% de.micro.all.sig$gene], collapse = ", "))
print(as.data.frame(de.micro.all.sig[de.micro.all.sig$gene %in% target_genes,
                                     c("gene", "logFC", "logCPM", "p_adj.glb")]), row.names = FALSE)
