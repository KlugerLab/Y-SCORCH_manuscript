# ==============================================================================
# Figure S14b: GO pathways of genes with higher gene activity (snATAC) in the
# Reactive-PWH astrocyte cluster, per region (the cluster was "HIV-Reactive" in
# the published run).
#
# PUBLISHED CODE  06-astrocyte_analysis/astrocyte_analysis.Rmd, section
#   "##S13b" (lines ~691-740), written out for INS only ("it's the same for all the
#   regions"):
#     1 ATAC cell names <sample>_<barcode> -> <barcode>-<sample>, matched to the
#       astrocyte RNA nuclei (multiome: one barcode, both modalities)
#     2 astrocyte subtype labels copied onto those ATAC nuclei
#     3 FindMarkers(ident.1 = reactive cluster, test.use = "wilcox", min.pct = 0.1)
#       against every other astrocyte subtype, on gene activity
#     4 genes with avg_log2FC > 0 and p_val_adj < 0.05 -> enrichGO BP, BH,
#       p < 0.01, q < 0.05 -> dotplot of the top 12 terms
#
# CHANGES FROM THE PUBLISHED CODE
#   * label "HIV-Reactive" -> "Reactive-PWH"
#   * gene-activity assay set explicitly (the Rmd relied on the object's default
#     assay; enrichGO keyType = SYMBOL only works if that was gene activity)
#   * enrichGO universe = genes tested (gene activity detected in >= 10% of nuclei
#     in either group), as for Fig 1E - not the whole genome
#   * all three regions, one job each
#   The per-nucleus Wilcoxon test is kept as published; it is NOT donor-aware.
#
# DONOR-AWARE VERSION. The Fig 1E method, on gene-activity counts:
#   revision/tools/astro_subtype_de_donor.R::astro_subtype_de(assay = <gene activity>) -
#   per donor, pseudobulk the Reactive-PWH nuclei and the other astrocytes
#   separately; a donor enters only if both halves have >= 10 nuclei; >= 5 donors
#   needed; edgeR ~ donor + half (paired, within person). GO BP on genes up at
#   FDR < 0.05, logFC > 0.25, universe = genes kept by filterByExpr.
#   Each region's raw donor x half pseudobulk is also saved, so
#   revision/tools/figS14b_pooled_donor.R can pool the regions without reloading anything.
#
# INPUT   OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled_no1040vst.RDS
#         OUT_ROOT/02-ATAC_processing/sampleQC_redo/<R>/data/signac_integrated_v2[_msn].RDS
# OUTPUT  OUT_ROOT/06-astrocyte_analysis/figures/S14b/
#           S14b_<R>_reactive_pwh_atac_pathways.png/.pdf
#           S14b_<R>_reactive_pwh_geneact_markers.csv   FindMarkers table
#           S14b_<R>_reactive_pwh_GO_BP.csv              enrichGO table
#
#           donor_aware/<R>/  paired edgeR table, donor counts, per-donor heatmap,
#                             S14b_<R>_donor_GO_BP.csv + _donor_atac_pathways.png
#           donor_aware/S14b_<R>_pseudobulk.rds   raw donor x half counts (for pooling)
#
# Usage: revision/docker/run.sh revision/tools/figS14b_reactive_atac_pathways.R --region INS
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(dplyr); library(ggplot2)
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_subtype_de_donor.R"))
suppressPackageStartupMessages(library(edgeR))

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
REGION <- toupper(getarg("--region", NA_character_))
if (is.na(REGION)) stop("--region is required")
TARGET <- "Reactive-PWH"
CLEAN_SUBTYPES <- c("Homeostatic", "Reactive", "Reactive-PWH", "Signaling",
                    "Synaptic", "Micro-Astro", "Mural-Astro")
OUT <- file.path(OUT_ROOT, "06-astrocyte_analysis", "figures", "S14b")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time(); el <- function() format(round(difftime(Sys.time(), t0, units = "mins"), 1))

## ---- astrocyte RNA nuclei of this region, with subtype labels ----------------
astro <- readRDS(file.path(OUT_ROOT, "06-astrocyte_analysis", "data",
                           "astro_all_clusters_labeled_no1040vst.RDS"))
am <- astro@meta.data
am$cell <- rownames(am)
am <- am[as.character(am$subcluster_labels) %in% CLEAN_SUBTYPES &
         grepl(paste0("-", REGION, "$"), am$cell), c("cell", "subcluster_labels")]
am$key <- sub(paste0("-", REGION, "$"), "", am$cell)          # <barcode>-1-<sample>
rm(astro); invisible(gc())
message(sprintf("[S14b] %s: %d astrocyte RNA nuclei in the 7 subtypes (%d %s)",
                REGION, nrow(am), sum(am$subcluster_labels == TARGET), TARGET))

## ---- ATAC object, cells renamed to the RNA key and matched --------------------
fn <- if (REGION == "VST") "signac_integrated_v2_msn.RDS" else "signac_integrated_v2.RDS"
atac <- readRDS(file.path(OUT_ROOT, "02-ATAC_processing", "sampleQC_redo", REGION, "data", fn))
message("[S14b] ATAC loaded (", el(), " min); assays: ", paste(Assays(atac), collapse = ", "))
GA <- intersect(c("GeneACT", "ACTIVITY", "RNA"), Assays(atac))[1]
if (is.na(GA)) stop("no gene-activity assay on the ATAC object")
## <sample>_<barcode>  ->  <barcode>-<sample>   (Rmd line ~696)
akey <- vapply(strsplit(colnames(atac), "_"), function(p) paste0(p[2], "-", p[1]), character(1))
hit  <- akey %in% am$key
message(sprintf("[S14b] ATAC nuclei matched to astrocyte RNA nuclei: %d of %d astrocytes",
                sum(hit), nrow(am)))
if (sum(hit) < 50) stop("too few matched nuclei - check the barcode conversion")
atac <- subset(atac, cells = colnames(atac)[hit])
atac$subtype <- am$subcluster_labels[match(akey[hit], am$key)]
DefaultAssay(atac) <- GA
Idents(atac) <- "subtype"
tab <- table(atac$subtype)
message("[S14b] matched nuclei per subtype: ", paste(names(tab), tab, sep = "=", collapse = ", "))
if (!TARGET %in% names(tab)) stop("no ", TARGET, " nuclei with ATAC in ", REGION)
## a region too small for the per-nucleus test still contributes its donors to the
## pooled donor-aware model, so it must reach the pseudobulk step below
cnt <- GetAssayData(atac, assay = GA, layer = "counts")
WILCOX_OK <- tab[[TARGET]] >= 10
if (!WILCOX_OK) message("[S14b] only ", tab[[TARGET]], " ", TARGET, " nuclei - Wilcoxon version skipped")

## ---- published test: Reactive-PWH vs all other astrocyte subtypes -------------
if (WILCOX_OK) {
mk <- FindMarkers(atac, ident.1 = TARGET, assay = GA, test.use = "wilcox", min.pct = 0.1)
mk$gene <- rownames(mk)
utils::write.csv(mk, file.path(OUT, sprintf("S14b_%s_reactive_pwh_geneact_markers.csv", REGION)),
                 row.names = FALSE)
up <- mk %>% dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05)
message(sprintf("[S14b] %d genes with higher gene activity in %s (of %d tested)", nrow(up), TARGET, nrow(mk)))

## universe: gene activity detected in >= 10% of nuclei in either group
in_t <- atac$subtype == TARGET
pct_t <- Matrix::rowMeans(cnt[, in_t, drop = FALSE] > 0)
pct_o <- Matrix::rowMeans(cnt[, !in_t, drop = FALSE] > 0)
universe <- rownames(cnt)[pct_t >= 0.1 | pct_o >= 0.1]

ego <- enrichGO(gene = up$gene, universe = universe, keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05,
                readable = FALSE)
res <- as.data.frame(ego)
utils::write.csv(res, file.path(OUT, sprintf("S14b_%s_reactive_pwh_GO_BP.csv", REGION)), row.names = FALSE)
message(sprintf("[S14b] %d GO BP terms (universe %d genes)", nrow(res), length(universe)))
if (!nrow(res)) {
  message("[S14b] no significant terms - no Wilcoxon dotplot")
} else {
  message("[S14b] top terms: ", paste(utils::head(res$Description, 5), collapse = " | "))
  p <- dotplot(ego, showCategory = 12) + ggtitle(sprintf("Top %s ATAC Reactive-PWH Pathways", REGION))
  f <- file.path(OUT, sprintf("S14b_%s_reactive_pwh_atac_pathways.png", REGION))
  ggsave(f, p, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(sub("\\.png$", ".pdf", f), p, width = 8, height = 6, device = grDevices::cairo_pdf, bg = "white")
  message("[S14b] Wilcoxon version -> ", f)
}
}   # WILCOX_OK

## ---- donor-aware: paired pseudobulk, the Fig 1E method ---------------------
DOUT <- file.path(OUT, "donor_aware")
dir.create(file.path(DOUT, REGION), recursive = TRUE, showWarnings = FALSE)
atac$condition_new <- as.character(atac$condition_new)

## raw donor x half pseudobulk for EVERY donor (no threshold), saved for pooling
don  <- as.character(atac$sample_ids)
half <- ifelse(atac$subtype == TARGET, "subtype", "other")
key  <- factor(paste(don, half, sep = "__"))
pb   <- as.matrix(cnt %*% Matrix::sparse.model.matrix(~ 0 + key))
colnames(pb) <- sub("^key", "", colnames(pb))
ncell <- as.data.frame(table(donor = don, half = half))
cond  <- tapply(atac$condition_new, don, function(x) x[1])
saveRDS(list(pb = pb, ncell = ncell, condition = cond, region = REGION, assay = GA),
        file.path(DOUT, sprintf("S14b_%s_pseudobulk.rds", REGION)))
message(sprintf("[S14b] saved donor x half pseudobulk: %d genes x %d profiles", nrow(pb), ncol(pb)))

dres <- tryCatch(
  astro_subtype_de(atac, subtype = TARGET, label_col = "subtype", donor_col = "sample_ids",
                   cond_col = "condition_new", assay = GA, outdir = file.path(DOUT, REGION)),
  error = function(e) { message("[S14b] donor-aware ", REGION, " not run: ", conditionMessage(e)); NULL })

if (!is.null(dres)) {
  tt <- dres$de
  up_d <- tt$gene[tt$FDR < 0.05 & tt$logFC > 0.25]
  message(sprintf("[S14b] donor-aware %s: %d donors, %d genes up in %s (FDR < 0.05, logFC > 0.25)",
                  REGION, sum(dres$donors$usable), length(up_d), TARGET))
  if (length(up_d) >= 5) {
    ego_d <- enrichGO(gene = up_d, universe = rownames(dres$dge), keyType = "SYMBOL",
                      OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                      pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    rd <- as.data.frame(ego_d)
    utils::write.csv(rd, file.path(DOUT, REGION, sprintf("S14b_%s_donor_GO_BP.csv", REGION)), row.names = FALSE)
    message(sprintf("[S14b] donor-aware GO: %d terms", nrow(rd)))
    if (nrow(rd)) {
      message("[S14b] donor-aware top terms: ", paste(utils::head(rd$Description, 5), collapse = " | "))
      pd <- dotplot(ego_d, showCategory = 12) +
        ggtitle(sprintf("Top %s ATAC Reactive-PWH Pathways (donor-aware)", REGION))
      fd <- file.path(DOUT, REGION, sprintf("S14b_%s_donor_atac_pathways.png", REGION))
      ggsave(fd, pd, width = 8, height = 6, dpi = 300, bg = "white")
      ggsave(sub("\\.png$", ".pdf", fd), pd, width = 8, height = 6, device = grDevices::cairo_pdf, bg = "white")
    }
  }
}
message("[S14b] ", REGION, " done in ", el(), " min")
