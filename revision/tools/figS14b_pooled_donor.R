# ==============================================================================
# Figure S14b, donor-aware and pooled across regions: GO pathways of genes with
# higher gene activity in Reactive-PWH astrocytes than in the donor's other astrocytes.
#
# Reads the raw donor x half gene-activity pseudobulks that
# revision/tools/figS14b_reactive_atac_pathways.R saves for each region
# (figures/S14b/donor_aware/S14b_<R>_pseudobulk.rds), sums each donor over the
# regions it contributes - exactly how Fig 1E pools astrocytes before its paired
# test - and runs the same paired edgeR as revision/tools/astro_subtype_de_donor.R:
#   a donor enters if both halves have >= 10 nuclei (summed over regions),
#   >= 5 donors required, DGEList -> filterByExpr(design) -> calcNormFactors ->
#   estimateDisp(robust) -> glmQLFit(robust) -> glmQLFTest on ~ donor + half.
# GO BP on genes up in Reactive-PWH (FDR < 0.05, logFC > 0.25), universe = genes
# kept by filterByExpr.
#
# This exists because Reactive-PWH is small: a single region may not reach 5
# donors with >= 10 matched nuclei, but the regions together can.
#
# INPUT   OUT_ROOT/06-astrocyte_analysis/figures/S14b/donor_aware/S14b_<R>_pseudobulk.rds
# OUTPUT  OUT_ROOT/06-astrocyte_analysis/figures/S14b/donor_aware/pooled/
#           S14b_pooled_donorDE.csv, S14b_pooled_donor_counts.csv,
#           S14b_pooled_donor_GO_BP.csv, S14b_pooled_donor_atac_pathways.png/.pdf
#
# Usage: revision/docker/run.sh revision/tools/figS14b_pooled_donor.R
# ==============================================================================

suppressPackageStartupMessages({
  library(edgeR); library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

MIN_CELLS_PER_HALF <- 10; MIN_DONORS_PAIRED <- 5     # as astro_subtype_de_donor.R
DIN <- file.path(OUT_ROOT, "06-astrocyte_analysis", "figures", "S14b", "donor_aware")
OUT <- file.path(DIN, "pooled")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

fs <- file.path(DIN, sprintf("S14b_%s_pseudobulk.rds", c("INS", "PFC", "VST")))
fs <- fs[file.exists(fs)]
if (!length(fs)) stop("no regional pseudobulks found in ", DIN)
L <- lapply(fs, readRDS)
message("[S14b-pool] regions: ", paste(vapply(L, `[[`, "", "region"), collapse = ", "))

genes <- Reduce(intersect, lapply(L, function(x) rownames(x$pb)))
cols  <- unique(unlist(lapply(L, function(x) colnames(x$pb))))
pb <- matrix(0, length(genes), length(cols), dimnames = list(genes, cols))
for (x in L) pb[, colnames(x$pb)] <- pb[, colnames(x$pb)] + x$pb[genes, , drop = FALSE]

nc <- do.call(rbind, lapply(L, function(x) x$ncell))
nc$Freq <- as.integer(nc$Freq)
nc <- aggregate(Freq ~ donor + half, data = nc, FUN = sum)
wide <- reshape(nc, idvar = "donor", timevar = "half", direction = "wide")
wide[is.na(wide)] <- 0
names(wide) <- sub("^Freq\\.", "n_", names(wide))
cond <- unlist(lapply(L, function(x) x$condition))
wide$condition <- as.character(pwh_relabel(cond[match(wide$donor, names(cond))]))
wide$usable <- wide$n_subtype >= MIN_CELLS_PER_HALF & wide$n_other >= MIN_CELLS_PER_HALF
utils::write.csv(wide, file.path(OUT, "S14b_pooled_donor_counts.csv"), row.names = FALSE)
use <- wide$donor[wide$usable]
message(sprintf("[S14b-pool] donors with >= %d nuclei in both halves: %d (PWH %d, PWoH %d)",
                MIN_CELLS_PER_HALF, length(use), sum(wide$usable & wide$condition == "PWH"),
                sum(wide$usable & wide$condition == "PWoH")))
if (length(use) < MIN_DONORS_PAIRED)
  stop("only ", length(use), " donors can contribute a pair even after pooling regions")

keep_cols <- as.vector(outer(use, c("subtype", "other"), paste, sep = "__"))
keep_cols <- keep_cols[keep_cols %in% colnames(pb)]
pb <- pb[, keep_cols, drop = FALSE]
meta <- data.frame(sample = colnames(pb), donor = sub("__.*$", "", colnames(pb)),
                   half = factor(sub("^.*__", "", colnames(pb)), levels = c("other", "subtype")))

y <- DGEList(counts = pb, samples = meta)
design <- model.matrix(~ donor + half, data = meta)
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
tt <- topTags(glmQLFTest(fit, coef = "halfsubtype"), n = Inf)$table
tt$gene <- rownames(tt)
utils::write.csv(tt[, c("gene", "logFC", "logCPM", "F", "PValue", "FDR")],
                 file.path(OUT, "S14b_pooled_donorDE.csv"), row.names = FALSE)
up <- tt$gene[tt$FDR < 0.05 & tt$logFC > 0.25]
message(sprintf("[S14b-pool] %d genes tested, %d up in Reactive-PWH (FDR < 0.05, logFC > 0.25); BCV %.3f",
                nrow(tt), length(up), sqrt(y$common.dispersion)))
if (length(up) < 5) { message("[S14b-pool] too few genes for GO"); quit(save = "no") }

ego <- enrichGO(gene = up, universe = rownames(y), keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
res <- as.data.frame(ego)
utils::write.csv(res, file.path(OUT, "S14b_pooled_donor_GO_BP.csv"), row.names = FALSE)
message(sprintf("[S14b-pool] %d GO BP terms", nrow(res)))
if (nrow(res)) {
  message("[S14b-pool] top terms: ", paste(utils::head(res$Description, 5), collapse = " | "))
  p <- dotplot(ego, showCategory = 12) +
    ggtitle("Top ATAC Reactive-PWH Pathways (donor-aware, regions pooled)")
  f <- file.path(OUT, "S14b_pooled_donor_atac_pathways.png")
  ggsave(f, p, width = 8, height = 6, dpi = 300, bg = "white")
  ggsave(sub("\\.png$", ".pdf", f), p, width = 8, height = 6, device = grDevices::cairo_pdf, bg = "white")
}
message("[S14b-pool] done -> ", OUT)
