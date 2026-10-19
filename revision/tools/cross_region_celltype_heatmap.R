# Cross-region top-DEG heatmap for one cell type (astrocytes by default), the same figure as the published microglia heatmap.
#
# Method: ports the heatmap section of the published astrocyte cross-region
#   script, which is identical to 03-DE_analysis/09-microglia_cross_region_heatmap_pathways.R
#   except for cp <- "Astro" (the astrocyte script is not in the repository).
#   The published script uses pb / pb.norm_mean without loading them; this tool
#   loads them explicitly and draws only the heatmap (no cnetplot).
#   Same as published:
#   genes    cross-region DE table, detected in > 5% of nuclei in either group,
#            |logFC| > 0.25, p_adj.glb < 0.05; top 25 down + top 25 up by
#            p_adj.glb (slice_min, ties kept), down block first, rows not clustered
#   values   pseudo_norm_mean_bulk.rds for the cell type, row z-scored
#   samples  > 30 nuclei of the cell type; ordered by condition, region, sample
#   styling  ComplexHeatmap, steelblue-white-red1 on [-4, 4], Condition + Region
#            annotation bars, the published region and condition colours
#   Differences:
#   - condition labels PWoH / PWH (person-first), same colours
#   - the DE table is chosen by SCORCH_XR_DE (revision/config/paths.R): the
#     cross-region edgeR table (with the negative-F fallback; `test` column
#     written to the CSV) or, with SCORCH_XR_DE=donor_aware, the donor-aware
#     limma-voom table; outputs then carry the _donor_aware suffix
#   - --rank fc, --include and --min-frq options (below)
#   - also writes a PDF and a CSV of exactly the genes and samples drawn
#
# Inputs:  OUT_ROOT/03-DE_analysis/cross_region/data/{DE table, pseudo_bulk.rds, pseudo_norm_mean_bulk.rds}
# Outputs: OUT_ROOT/03-DE_analysis/cross_region/figures/celltype_heatmaps/
#   <CT>_top<N>_heatmap.png / .pdf
#   <CT>_top<N>_heatmap_genes.csv     gene, direction, logFC, p_adj.glb, test
#   <CT>_top<N>_heatmap_samples.csv   sample, region, condition, n nuclei, kept
#
# Usage: revision/docker/run.sh revision/tools/cross_region_celltype_heatmap.R [--celltype Astro]
#          [--n 25] [--min-cells 30] [--rank padj|fc] [--include GENE1,GENE2] [--min-frq 0.05]

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({
  library(SingleCellExperiment); library(dplyr)
  library(ComplexHeatmap); library(circlize); library(grid)
})

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
CP        <- getarg("--celltype", "Astro")
N_TOP     <- as.integer(getarg("--n", "25"))
MIN_CELLS <- as.integer(getarg("--min-cells", "30"))

DATA   <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "data")
OUTDIR <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "figures", "celltype_heatmaps")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --rank fc / --include (the Fig 1E ranking rule). The published rule ranks by
##   p_adj.glb; under the donor-aware model that pushes SERPINA3, CHI3L1 (Astro) and CD163
##   (Micro) out of the top 25 although they stay significant. fc = top N significant
##   genes by |logFC| per direction; --include adds named genes to their block when they
##   are significant and not already drawn (never forced in if not significant).
RANK    <- match.arg(getarg("--rank", "padj"), c("padj", "fc"))
INCLUDE <- trimws(strsplit(getarg("--include", ""), ",")[[1]])
## --min-frq: minimum fraction
##   of nuclei expressing the gene in either group for a gene to be ELIGIBLE AS A ROW.
##   Default 0.05 = the DEG rule. 0.10 drops the low-expression genes (often unnamed
##   lncRNAs) that fold-change ranking favours; the DEG calls themselves are unchanged.
MIN_FRQ <- as.numeric(getarg("--min-frq", "0.05"))
RANK_SFX <- paste0(if (RANK == "fc") "_fc" else "", if (length(INCLUDE)) "_plusnamed" else "",
                   if (MIN_FRQ != 0.05) sprintf("_frq%02d", round(100 * MIN_FRQ)) else "")
stem <- file.path(OUTDIR, sprintf("%s_top%d_heatmap%s%s", make.names(CP), N_TOP, RANK_SFX, XR_OUT_SFX))   # + _donor_aware under SCORCH_XR_DE

## ---- inputs (the published script used these without loading them) -----------
de_edgeR     <- readRDS(xr_de_file(dir = DATA))   # edgeR table, or donor-aware under SCORCH_XR_DE=donor_aware (revision/config/paths.R)
pb           <- readRDS(file.path(DATA, "pseudo_bulk.rds"))
pb.norm_mean <- readRDS(file.path(DATA, "pseudo_norm_mean_bulk.rds"))
stopifnot(identical(colnames(pb), colnames(pb.norm_mean)), CP %in% assayNames(pb.norm_mean))

## ---- genes: published selection ------------------------------------------------
top_edgeR_all <- de_edgeR %>% ungroup() %>%
  dplyr::filter((avg_CTR_frq > MIN_FRQ) | (avg_HIV_frq > MIN_FRQ)) %>%
  dplyr::filter((abs(logFC) > 0.25) & (p_adj.glb < 0.05), cluster_id == CP)

pick <- function(x) if (RANK == "fc") slice_max(x, abs(logFC), n = N_TOP) else slice_min(x, p_adj.glb, n = N_TOP)
add_named <- function(sel, pool) {
  extra <- pool %>% dplyr::filter(gene %in% INCLUDE, !gene %in% sel$gene)
  if (nrow(extra)) message("[hm] added named gene(s): ", paste(extra$gene, collapse = ", "))
  bind_rows(sel, extra)
}
up_pool   <- top_edgeR_all %>% dplyr::filter(logFC > 0)
down_pool <- top_edgeR_all %>% dplyr::filter(logFC < 0)
top_edgeR_target.up   <- add_named(pick(up_pool), up_pool)
top_edgeR_target.down <- add_named(pick(down_pool), down_pool)
if (length(INCLUDE)) {
  ns <- setdiff(INCLUDE, top_edgeR_all$gene)
  if (length(ns)) message("[hm] named but NOT significant in ", CP, " (not drawn): ", paste(ns, collapse = ", "))
}
top_edgeR_target <- rbind(top_edgeR_target.down, top_edgeR_target.up)
message(sprintf("[hm] %s: %d up / %d down significant; drawing %d down + %d up",
                CP, sum(top_edgeR_all$logFC > 0), sum(top_edgeR_all$logFC < 0),
                nrow(top_edgeR_target.down), nrow(top_edgeR_target.up)))

## ---- samples: published order and nuclei filter ---------------------------------
n_cells <- do.call(rbind, pb@int_colData$n_cells)
coldf <- data.frame(
  region    = pb@colData$region,
  condition = pwh_relabel(pb@colData$condition_new),
  sample_id = colnames(pb),
  n_nuclei  = n_cells[colnames(pb), CP])
rownames(coldf) <- colnames(pb.norm_mean)
coldf <- coldf %>% arrange(condition, region, sample_id)
coldf$kept <- coldf$n_nuclei > MIN_CELLS
write.csv(coldf, paste0(stem, "_samples.csv"), row.names = FALSE)
message("[hm] samples dropped (<= ", MIN_CELLS, " nuclei): ",
        if (any(!coldf$kept)) paste(coldf$sample_id[!coldf$kept], collapse = ", ") else "none")
coldf <- coldf[coldf$kept, ]

mat <- pb.norm_mean@assays@data[[CP]][top_edgeR_target$gene, rownames(coldf)]
mat <- t(scale(t(as.matrix(mat))))
if (any(!is.finite(mat))) {
  message("[hm] ", sum(apply(!is.finite(mat), 1, any)), " zero-variance row(s) set to 0")
  mat[!is.finite(mat)] <- 0
}

write.csv(top_edgeR_target %>%
            transmute(gene, direction = ifelse(logFC > 0, "up in PWH", "down in PWH"),
                      logFC, logCPM, p_adj.glb,
                      test = if ("test" %in% names(.)) test else NA_character_),
          paste0(stem, "_genes.csv"), row.names = FALSE)

## ---- plot: published styling ----------------------------------------------------
condition_cols <- PWH_COLOURS
region_cols    <- c(INS = "#E69F00", PFC = "#56B4E9", VST = "#CC79A7")
col_fun <- colorRamp2(c(-4, 0, 4), c("steelblue", "white", "red1"))

p1 <- ComplexHeatmap::Heatmap(mat,
  col = col_fun,
  column_names_gp = gpar(fontsize = 25),
  row_names_gp = gpar(fontsize = 30),
  cluster_columns = FALSE, cluster_rows = FALSE,
  show_row_dend = FALSE,
  name = " ",
  column_gap = unit(0.1, "mm"),
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    "Condition" = as.character(coldf$condition),
    "Region"    = coldf$region,
    col = list(Region = region_cols, Condition = condition_cols),
    simple_anno_size = unit(2, "cm"),
    annotation_name_gp = gpar(fontsize = 30, fontface = "bold"),
    annotation_legend_param = list(
      Region    = list(labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 30, fontface = "bold")),
      Condition = list(labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 30, fontface = "bold")))),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 25)),
  show_column_names = FALSE)

png(filename = paste0(stem, ".png"), width = 6000, height = 6000, res = 300)
print(p1)
dev.off()
pdf(paste0(stem, ".pdf"), width = 20, height = 20)
print(p1)
dev.off()

message(sprintf("[hm] %d genes x %d samples (INS %d, PFC %d, VST %d)", nrow(mat), ncol(mat),
                sum(coldf$region == "INS"), sum(coldf$region == "PFC"), sum(coldf$region == "VST")))
message("[hm] wrote ", stem, ".png/.pdf/_genes.csv/_samples.csv")
message("[hm] done")
