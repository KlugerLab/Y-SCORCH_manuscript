# Donor-aware cross-region DE (primary cross-region model): the published estimand, with repeated regional measures from the same donor modelled.
#
# Method: the published 03-DE_analysis/07-cross_region_edgeR.R aggregates to
#   region.donor samples and fits
#       ~ 0 + region_HIV + Age + Sex + Race + PMI
#   then contrasts the average within-region HIV effect:
#       (INS.HIV-INS.CTR)/3 + (PFC.HIV-PFC.CTR)/3 + (VST.HIV-VST.CTR)/3
#   A donor contributing three regions appears as three independent samples, so
#   n is inflated up to threefold and the standard errors are too small
#   (reviewer request: whether donor pairing across regions was modelled).
#   A donor fixed effect is not an option: HIV status is constant within a donor,
#   so it is collinear with the contrast and absorbs it. The repetition is
#   handled as a correlation structure instead: limma-voom on the same
#   pseudobulk, design and contrast, with duplicateCorrelation(block = donor).
#   The estimand is unchanged; only the variance is corrected.
#   Cell-type samples need >= 10 nuclei (as pbDS(min_cells = 10)); BH is applied
#   globally across cell types (p_adj.glb), as in the edgeR tables.
#   The consensus correlation answers the reviewer's question directly: near zero
#   means regional measures from one donor are effectively independent; clearly
#   positive means the donor-naive model was anticonservative.
#
# Inputs:  OUT_ROOT/03-DE_analysis/cross_region/data/{pseudo_bulk.rds, mm.rds,
#          all_HIV-CTR_res_edgeR_filtered_withFreq.RDS (donor-naive edgeR, for comparison)}
# Outputs: OUT_ROOT/03-DE_analysis/cross_region/donor_aware/
#          cross_region_donor_aware_DE.csv, intra_donor_correlation.csv,
#          edgeR_vs_donor_aware.csv, edgeR_vs_donor_aware_by_celltype.csv
#
# Usage: revision/docker/run.sh revision/tools/cross_region_donor_aware.R

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(SingleCellExperiment)
  library(limma); library(edgeR); library(dplyr)
})

MIN_CELLS <- 10
FDR       <- 0.05
PB_FILE   <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "data", "pseudo_bulk.rds")
MM_FILE   <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "data", "mm.rds")
## The comparison baseline is the reprocessed cross-region edgeR table
##   (donor-naive), not the published paper result.
EDGER_FILE <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "data",
                        "all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
OUTDIR    <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "donor_aware")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

CONTRAST_STR <- "(INS.HIV-INS.CTR)/3+(PFC.HIV-PFC.CTR)/3+(VST.HIV-VST.CTR)/3"

for (f in c(PB_FILE, MM_FILE)) if (!file.exists(f)) stop("missing input: ", f)

message("[xr] loading pseudobulk and design")
pb <- readRDS(PB_FILE)
mm <- readRDS(MM_FILE)

## sample ids are "<REGION>.<donor>" - the donor is what repeats
stopifnot(!is.null(rownames(mm)))
sid    <- rownames(mm)
donor  <- sub("^[^.]+\\.", "", sid)
region <- sub("\\..*$", "", sid)

message("[xr] ", length(sid), " region-donor samples, ",
        length(unique(donor)), " unique donors, ",
        length(unique(region)), " regions")
rep_tab <- table(table(donor))
message("[xr] donors contributing 1/2/3 regions: ",
        paste(sprintf("%s region(s)=%d", names(rep_tab), rep_tab), collapse = "  "))

contrast <- limma::makeContrasts(contrasts = CONTRAST_STR, levels = mm)

nc <- int_colData(pb)$n_cells          # sample -> named vector of cells per type
cell_types <- assayNames(pb)
message("[xr] ", length(cell_types), " cell types")

# ------------------------------------------------------------------------------
run_one <- function(k) {
  y <- assay(pb, k)
  y <- y[, sid, drop = FALSE]                       # align to the design

  ## same min_cells rule as pbDS(min_cells = 10)
  nck <- vapply(sid, function(s) {
    v <- nc[[s]]; if (is.null(v) || !k %in% names(v)) 0L else as.integer(v[[k]])
  }, integer(1))
  keep_s <- nck >= MIN_CELLS
  if (sum(keep_s) < 8) return(NULL)

  d  <- mm[keep_s, , drop = FALSE]
  d  <- d[, colSums(abs(d)) > 0, drop = FALSE]      # drop collapsed levels
  if (qr(d)$rank < ncol(d)) return(NULL)
  ## the contrast must still be expressible after dropping columns
  ctr <- tryCatch(limma::makeContrasts(contrasts = CONTRAST_STR, levels = d),
                  error = function(e) NULL)
  if (is.null(ctr)) return(NULL)

  dge <- DGEList(counts = y[, keep_s, drop = FALSE])
  keep_g <- filterByExpr(dge, d)
  if (sum(keep_g) < 50) return(NULL)
  dge <- dge[keep_g, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  blk <- donor[keep_s]
  v   <- voom(dge, d)

  ## duplicateCorrelation twice is the documented limma idiom: the first pass
  ## uses voom weights computed without the block, the second re-estimates both
  ## with the correlation in hand.
  dc1 <- tryCatch(duplicateCorrelation(v, d, block = blk), error = function(e) NULL)
  if (is.null(dc1) || !is.finite(dc1$consensus.correlation)) return(NULL)
  v2  <- voom(dge, d, block = blk, correlation = dc1$consensus.correlation)
  dc2 <- duplicateCorrelation(v2, d, block = blk)
  rho <- dc2$consensus.correlation
  if (!is.finite(rho)) return(NULL)

  fit <- lmFit(v2, d, block = blk, correlation = rho)
  fit <- contrasts.fit(fit, ctr)
  fit <- eBayes(fit)
  tt  <- topTable(fit, number = Inf, sort.by = "none")

  data.frame(cluster_id = k, gene = rownames(tt),
             logFC = tt$logFC, p_val = tt$P.Value, p_adj.loc = tt$adj.P.Val,
             rho = rho, n_samples = sum(keep_s),
             n_donors = length(unique(blk)),
             stringsAsFactors = FALSE)
}

message("[xr] fitting donor-aware model per cell type")
res <- list()
for (k in cell_types) {
  out <- tryCatch(run_one(k), error = function(e) {
    message("[xr] ", k, ": FAILED - ", conditionMessage(e)); NULL })
  if (is.null(out)) { message("[xr] ", k, ": not estimable"); next }
  message(sprintf("[xr] %-14s rho = %+.3f   %d samples / %d donors   %d sig",
                  k, out$rho[1], out$n_samples[1], out$n_donors[1],
                  sum(out$p_adj.loc < FDR, na.rm = TRUE)))
  res[[k]] <- out
}
if (!length(res)) stop("no cell type was estimable")
R <- dplyr::bind_rows(res)
## global BH across cell types, the same p_adj.glb the edgeR tables and every
## downstream figure use (topTable's adj.P.Val is per cell type = p_adj.loc)
R$p_adj.glb <- p.adjust(R$p_val, method = "BH")

# ---- the headline: how much does donor repetition matter? --------------------
rho_tab <- R %>% group_by(cluster_id) %>%
  summarise(rho = rho[1], n_donors = n_donors[1], n_samples = n_samples[1],
            n_sig_donor_aware = sum(p_adj.glb < FDR, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(rho))
cat("\n=== INTRA-DONOR CORRELATION ACROSS REGIONS ===\n")
cat("rho > 0 means regional measures from one donor are correlated, i.e. the\n")
cat("published model (which treated them as independent) was anticonservative.\n\n")
print(as.data.frame(rho_tab), row.names = FALSE)
cat(sprintf("\nmedian rho = %+.3f   (n = %d cell types)\n",
            median(rho_tab$rho), nrow(rho_tab)))

# ---- comparison against the donor-naive cross-region edgeR table -------------
## Both sides on p_adj.glb (global BH across cell types).
if (file.exists(EDGER_FILE)) {
  eg <- as.data.frame(readRDS(EDGER_FILE))
  eg <- data.frame(cluster_id = as.character(eg$cluster_id), gene = eg$gene,
                   logFC_edgeR = eg$logFC, p_adj.glb_edgeR = eg$p_adj.glb,
                   test_edgeR = if ("test" %in% colnames(eg)) eg$test else "QL",
                   stringsAsFactors = FALSE)
  cmp <- merge(R, eg, by = c("cluster_id", "gene"), all = TRUE)
  cmp$sig_edgeR <- !is.na(cmp$p_adj.glb_edgeR) & cmp$p_adj.glb_edgeR < FDR
  cmp$sig_donor <- !is.na(cmp$p_adj.glb)       & cmp$p_adj.glb       < FDR
  both <- cmp$sig_edgeR & cmp$sig_donor
  tot <- data.frame(
    n_sig_edgeR       = sum(cmp$sig_edgeR),
    n_sig_donor_aware = sum(cmp$sig_donor),
    n_both            = sum(both),
    pct_edgeR_kept    = round(100 * sum(both) / max(1, sum(cmp$sig_edgeR)), 1),
    only_donor_aware  = sum(cmp$sig_donor & !cmp$sig_edgeR),
    sign_agree_pct    = round(100 * mean(sign(cmp$logFC[both]) == sign(cmp$logFC_edgeR[both]),
                                         na.rm = TRUE), 1))
  per_ct <- cmp %>% group_by(cluster_id) %>%
    summarise(sig_edgeR = sum(sig_edgeR), sig_donor_aware = sum(sig_donor),
              both = sum(sig_edgeR & sig_donor), .groups = "drop") %>%
    mutate(pct_edgeR_kept = round(100 * both / pmax(1, sig_edgeR), 1)) %>%
    arrange(desc(sig_edgeR))
  cat("\n=== cross-region edgeR (donor-naive, rerun) vs DONOR-AWARE, p_adj.glb < ", FDR, " ===\n", sep = "")
  print(as.data.frame(tot), row.names = FALSE)
  print(as.data.frame(per_ct), row.names = FALSE)
  cat("\n=== reactive-astrocyte panel ===\n")
  print(cmp[cmp$cluster_id == "Astro" & cmp$gene %in% c("SERPINA3", "OSMR", "CHI3L1", "CD44", "GFAP", "C3"),
            c("gene", "logFC", "p_adj.glb", "logFC_edgeR", "p_adj.glb_edgeR", "test_edgeR")], row.names = FALSE)
  utils::write.csv(cmp,    file.path(OUTDIR, "edgeR_vs_donor_aware.csv"), row.names = FALSE)
  utils::write.csv(per_ct, file.path(OUTDIR, "edgeR_vs_donor_aware_by_celltype.csv"), row.names = FALSE)
} else message("[xr] no cross-region edgeR table to compare against: ", EDGER_FILE)

utils::write.csv(R,       file.path(OUTDIR, "cross_region_donor_aware_DE.csv"), row.names = FALSE)
utils::write.csv(rho_tab, file.path(OUTDIR, "intra_donor_correlation.csv"), row.names = FALSE)
message("\n[xr] wrote ", OUTDIR)
message("[xr] done")
