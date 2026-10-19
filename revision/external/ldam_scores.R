# ==============================================================================
# LDAM signature score in microglia (analysis by Junchen).
# ssGSEA (GSVA) of the four-gene LDAM set (NAMPT, ACSL1, DPYD, CD163; Haney et al.
# 2024) on per-donor microglia pseudobulk, per region; PWH vs PWoH by t-test.
# Output: OUT_ROOT/08-communication/LDAM/LDAM_signature_boxplot.png
# Outlier-donor sensitivity: revision/tools/ldam_ssgsea_sensitivity.R
# Usage: revision/docker/run.sh revision/external/ldam_scores.R
# ==============================================================================
library(Seurat)
#BiocManager::install("GSVA")
library(GSVA)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggpubr)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
LDAM_OUT <- file.path(OUT_ROOT, "08-communication", "LDAM")
dir.create(LDAM_OUT, recursive = TRUE, showWarnings = FALSE)

seurat.obj <- readRDS(file.path(OUT_ROOT, "01-Regional_data_processing/combined/data/seurat_combined_v3_no1040_msn.RDS"))

seurat.obj.micro <- subset(seurat.obj,cells = colnames(seurat.obj)[seurat.obj@meta.data$consensus_annotation_v3 == "Micro"])

# extract the counts, genes, and metadata
seurat.obj.micro <- JoinLayers(seurat.obj.micro)

pseudo_counts <- AggregateExpression(
  object = seurat.obj.micro, 
  assays = "RNA", 
  slot = "counts", 
  group.by = c("region","condition_new", "sample_ids")
)


# ldam_signature: APOE4/4 is linked to damaging lipid droplets in Alzheimer’s disease microglia
ldam_signature <- list("LDAM"=c("NAMPT","ACSL1","DPYD","CD163"))

ssgsea_param <- ssgseaParam(exprData = pseudo_counts$RNA, 
                            geneSets = ldam_signature)
ssgsea_scores <- gsva(ssgsea_param)

# plotting
region_cols <- c("#E69F00", "#56B4E9","#CC79A7")
names(region_cols) <-c("INS","PFC","VST")

cond_col <- c("#0072B2", "darkorange1")
names(cond_col) <-c("PWoH","PWH")

ssgsea_df <- as.data.frame(t(ssgsea_scores)) %>%
  rownames_to_column(var="sample_ids") %>%
  separate("sample_ids",into = c("region","condition","sample_id")) %>%
  mutate(condition = ifelse(condition == "HIV","PWH","PWoH"))

p <- ggplot(ssgsea_df, aes(x = condition, y = LDAM, fill = condition)) +
  facet_grid(~region) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "t.test", label = "p.format") +
  scale_fill_manual(values = cond_col) +
  ylab("LDAM signature score")+
  NoLegend()+
  theme_minimal(base_size = 15)+
  theme(legend.position = "none")

png(filename = file.path(LDAM_OUT, "LDAM_signature_boxplot.png"),width = 2500,height = 1500,res=300)
print(p)
dev.off()
