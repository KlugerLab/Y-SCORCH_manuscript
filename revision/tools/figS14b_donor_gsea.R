# ==============================================================================
# Figure S14b, donor-aware PATHWAY test by gene-set enrichment (GSEA) on the
# paired Reactive-PWH gene-activity DE.
#
# WHY  The paired within-donor test (revision/tools/figS14b_reactive_atac_pathways.R, the Fig 1E
#   method on gene activity) leaves 2-27 genes per region and 124 pooled at
#   FDR < 0.05 & logFC > 0.25 - too few for over-representation, 0 GO terms. A
#   threshold throws away the modest, coordinated shifts a pathway test is meant
#   to find. GSEA uses EVERY tested gene, ranked by the paired donor-level
#   statistic, so it stays donor-aware without a gene cut-off.
#
# INPUT  the paired edgeR tables already written (no object is loaded):
#          OUT_ROOT/06-astrocyte_analysis/figures/S14b/donor_aware/<R>/astro_Reactive_PWH_donorDE.csv
#          OUT_ROOT/06-astrocyte_analysis/figures/S14b/donor_aware/pooled/S14b_pooled_donorDE.csv
# RANK   sign(logFC) * sqrt(F) - the signed square root of the QL F statistic,
#        i.e. the paired t-like statistic, Reactive-PWH vs the donor's other
#        astrocytes.
# TEST   clusterProfiler::gseGO, BP, SYMBOL, gene sets of 10-500, BH, seed 1.
# PLOT   top 12 pathways with NES > 0 (higher gene activity in Reactive-PWH - the
#        direction of the published panel); the full table keeps both directions.
#
# OUTPUT  OUT_ROOT/06-astrocyte_analysis/figures/S14b/donor_aware/gsea/
#           S14b_<R>_donor_gsea.csv / .png / .pdf   (R = INS, PFC, VST, pooled)
#
# Usage: revision/docker/run.sh revision/tools/figS14b_donor_gsea.R
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(ggplot2); library(dplyr)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

DIN <- file.path(OUT_ROOT, "06-astrocyte_analysis", "figures", "S14b", "donor_aware")
OUT <- file.path(DIN, "gsea")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SETS <- c(INS = file.path(DIN, "INS", "astro_Reactive_PWH_donorDE.csv"),
          PFC = file.path(DIN, "PFC", "astro_Reactive_PWH_donorDE.csv"),
          VST = file.path(DIN, "VST", "astro_Reactive_PWH_donorDE.csv"),
          pooled = file.path(DIN, "pooled", "S14b_pooled_donorDE.csv"))

for (nm in names(SETS)) {
  f <- SETS[[nm]]
  if (!file.exists(f)) { message("[gsea] ", nm, ": no paired table - skipped"); next }
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d <- d[is.finite(d$F) & is.finite(d$logFC), ]
  stat <- sign(d$logFC) * sqrt(pmax(d$F, 0))
  ## tapply returns a 1-d array; gseaResult's geneList slot must be a plain numeric
  v  <- tapply(stat, d$gene, function(x) x[which.max(abs(x))])
  gl <- sort(setNames(as.numeric(v), names(v)), decreasing = TRUE)
  set.seed(1)
  g <- gseGO(geneList = gl, ont = "BP", keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH",
             eps = 0, seed = TRUE, verbose = FALSE)
  res <- as.data.frame(g)
  utils::write.csv(res, file.path(OUT, sprintf("S14b_%s_donor_gsea.csv", nm)), row.names = FALSE)
  n_up <- sum(res$NES > 0); n_dn <- sum(res$NES < 0)
  message(sprintf("[gsea] %s: %d genes ranked; %d pathways at FDR < 0.05 (%d higher in Reactive-PWH, %d lower)",
                  nm, length(gl), nrow(res), n_up, n_dn))
  if (n_up > 0) {
    message("[gsea]   top up: ", paste(head(res$Description[res$NES > 0], 6), collapse = " | "))
    up_ids <- head(res$ID[res$NES > 0], 12)
    gu <- g; gu@result <- gu@result[gu@result$ID %in% up_ids, ]
    ttl <- if (nm == "pooled") "Reactive-PWH ATAC pathways, donor-aware GSEA (regions pooled)"
           else sprintf("Top %s ATAC Reactive-PWH Pathways (donor-aware GSEA)", nm)
    ## long GO names wrap onto 2-3 lines, so rows get more room and the wrap
    ## width is fixed rather than left to the default
    p <- dotplot(gu, showCategory = 12, x = "NES", label_format = 45) + ggtitle(ttl) +
      theme(axis.text.y = element_text(size = 11, lineheight = 0.9))
    n_terms <- nrow(gu@result)
    h <- max(6, 2.2 + 0.62 * n_terms)
    fp <- file.path(OUT, sprintf("S14b_%s_donor_gsea.png", nm))
    ggsave(fp, p, width = 9, height = h, dpi = 300, bg = "white")
    ggsave(sub("\\.png$", ".pdf", fp), p, width = 9, height = h, device = grDevices::cairo_pdf, bg = "white")
  }
  if (n_dn > 0) message("[gsea]   top down: ", paste(head(res$Description[res$NES < 0], 4), collapse = " | "))
}
message("[gsea] done -> ", OUT)
