# Pooled, donor-aware chromVAR motif activity test across regions, PWH vs PWoH (Figure 2 motif panel).
#
# Why: per region (revision/tools/atac_motif_chromvar.R) no astrocyte motif passed
#   FDR across 746 motifs with ~30 donors, but AP-1, STAT, BACH/NFE2 and CEBP
#   activity was higher in PWH in all three regions (INS strongest). Pooling the
#   regions lets a consistent effect accumulate, as in the donor-aware
#   cross-region RNA DE (revision/tools/cross_region_donor_aware.R).
#
# Model: mirrors revision/tools/cross_region_donor_aware.R, on chromVAR scores
#   data      each donor's mean chromVAR deviation z-score, per region
#             (<R>/<CT>_chromvar_donor_means.csv; donors with >= 20 nuclei)
#   design    the cross-region DE design matrix (03-DE_analysis/cross_region/
#             data/mm.rds): region x condition means + Age + Sex + Race + PMI,
#             rows "<REGION>.<donor>"
#   contrast  the average within-region HIV effect,
#             (INS.HIV-INS.CTR)/3 + (PFC.HIV-PFC.CTR)/3 + (VST.HIV-VST.CTR)/3
#             (over the regions present for that cell type)
#   donors    a donor contributing several regions is one person:
#             duplicateCorrelation(block = donor), then lmFit with that
#             correlation. limma on the scores directly (they are continuous
#             z-scores, so no voom).
#   Per-region effects are read off the same fit, for the heatmap columns.
#   Three versions are fitted - plain, CENTERED (primary: each sample's mean
#   deviation over all motifs removed, so global shifts do not pass as TF
#   activity) and WEIGHTED (limma arrayWeights) - see the note in the loop.
#   BH across the 746 motifs, per cell type.
#
# Inputs:  OUT_ROOT/04-ATAC_DE/motif_chromvar/<REGION>/{<CT>_chromvar_donor_means.csv, motif_names.csv}
#          OUT_ROOT/03-DE_analysis/cross_region/data/mm.rds
# Outputs: OUT_ROOT/04-ATAC_DE/motif_chromvar/crossregion/
#   <CT>_chromvar_crossregion.csv   motif, name, pooled PWH-PWoH, t, P, FDR,
#                                   per-region differences and P, n, rho
#
# Usage: revision/docker/run.sh revision/tools/atac_motif_crossregion.R [--celltypes Astro,Micro,OPC]

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(limma); library(dplyr) })

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
CTS  <- trimws(strsplit(getarg("--celltypes", "Astro,Micro,OPC"), ",")[[1]])
BASE <- file.path(OUT_ROOT, "04-ATAC_DE", "motif_chromvar")
OUT  <- file.path(BASE, "crossregion")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
mm <- readRDS(file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "data", "mm.rds"))

for (ct in CTS) {
  message("\n[xr-motif] ===== ", ct, " =====")
  regs <- c("INS", "PFC", "VST")
  ## file names are sanitised the same way as in atac_motif_chromvar.R, because
  ## names like "L2/3 IT" and "D1/D2-Hybrid" contain a slash
  safe <- gsub("[^A-Za-z0-9_.+-]", "_", ct)
  fs <- setNames(file.path(BASE, regs, paste0(safe, "_chromvar_donor_means.csv")), regs)
  regs <- regs[file.exists(fs)]
  if (length(regs) < 2) { message("[xr-motif] fewer than 2 regions finished - skipped"); next }
  mats <- lapply(regs, function(R) {
    x <- utils::read.csv(fs[[R]], check.names = FALSE, stringsAsFactors = FALSE)
    m <- as.matrix(x[, -(1:2), drop = FALSE]); rownames(m) <- x$motif
    colnames(m) <- paste0(R, ".", colnames(m)); m
  })
  motifs <- Reduce(intersect, lapply(mats, rownames))
  M <- do.call(cbind, lapply(mats, function(m) m[motifs, , drop = FALSE]))
  names_map <- utils::read.csv(file.path(BASE, regs[1], "motif_names.csv"), stringsAsFactors = FALSE)
  motif_name <- setNames(names_map$name, names_map$motif)

  keep <- colnames(M) %in% rownames(mm)
  if (any(!keep)) message("[xr-motif] not in the cross-region design (dropped): ", paste(colnames(M)[!keep], collapse = ", "))
  M <- M[, keep, drop = FALSE]
  design <- mm[colnames(M), , drop = FALSE]
  design <- design[, colSums(design != 0) > 0, drop = FALSE]
  ne <- limma::nonEstimable(design)
  if (!is.null(ne)) { message("[xr-motif] dropping non-estimable: ", paste(ne, collapse = ", "))
                      design <- design[, setdiff(colnames(design), ne), drop = FALSE] }
  donor <- sub("^[^.]+\\.", "", colnames(M))
  message("[xr-motif] regions: ", paste(regs, collapse = ", "), "; ", ncol(M), " region-donor samples, ",
          length(unique(donor)), " donors (",
          paste(sprintf("%s: %d", regs, sapply(regs, function(R) sum(startsWith(colnames(M), paste0(R, "."))))), collapse = ", "), ")")

  ## pooled + per-region contrasts, all from one fit
  pooled <- paste(sprintf("(%s.HIV-%s.CTR)/%d", regs, regs, length(regs)), collapse = "+")
  per    <- setNames(sprintf("%s.HIV-%s.CTR", regs, regs), regs)
  ctr <- limma::makeContrasts(contrasts = c(pooled, per), levels = design)
  colnames(ctr) <- c("pooled", regs)

  ## Three versions of the same model. In the plain fit INS Micro (402 motifs down,
  ##   driven by donor HCTMU) and INS OPC (463 up, 96% of all motifs shifted the
  ##   same way) are global, not motif-specific, shifts:
  ##   plain     as specified above
  ##   centered  each sample's mean deviation over all motifs subtracted first, so
  ##             a donor-wide global shift cannot pass for TF-specific activity -
  ##             tests what is specific to a motif family (primary)
  ##   weighted  limma arrayWeights: data-driven reliability weight per sample,
  ##             down-weighting outlier donors without dropping them by hand
  fit_variant <- function(X, weights = NULL) {
    dc  <- duplicateCorrelation(X, design, block = donor, weights = weights)
    rho <- dc$consensus.correlation
    list(fit = eBayes(contrasts.fit(lmFit(X, design, block = donor, correlation = rho, weights = weights), ctr)),
         rho = rho)
  }
  Mc <- sweep(M, 2, colMeans(M))                       # per-sample global offset removed
  aw <- arrayWeights(M, design)
  message("[xr-motif] lowest sample weights: ",
          paste(sprintf("%s=%.2f", colnames(M)[order(aw)][1:4], sort(aw)[1:4]), collapse = ", "))
  V <- list(plain = fit_variant(M), centered = fit_variant(Mc), weighted = fit_variant(M, weights = aw))
  rho <- V$plain$rho
  message(sprintf("[xr-motif] within-donor correlation across regions: %.3f", rho))

  tp <- topTable(V$centered$fit, coef = "pooled", number = Inf, sort.by = "none")
  res <- data.frame(motif = rownames(tp), name = motif_name[rownames(tp)],
                    diff_pooled = tp$logFC, t = tp$t, P = tp$P.Value, FDR = tp$adj.P.Val,
                    stringsAsFactors = FALSE)                             # primary = centered
  for (R in regs) {
    tr <- topTable(V$centered$fit, coef = R, number = Inf, sort.by = "none")
    res[[paste0("diff_", R)]] <- tr$logFC
    res[[paste0("P_", R)]]    <- tr$P.Value
  }
  for (v in c("plain", "weighted")) {
    tv <- topTable(V[[v]]$fit, coef = "pooled", number = Inf, sort.by = "none")
    res[[paste0("diff_", v)]] <- tv$logFC; res[[paste0("FDR_", v)]] <- tv$adj.P.Val
  }
  res$global_shift_pooled <- mean(topTable(V$plain$fit, coef = "pooled", number = Inf)$logFC)
  res$n_samples <- ncol(M); res$n_donors <- length(unique(donor)); res$rho <- rho; res$cell_type <- ct
  res <- res[order(res$P), ]
  utils::write.csv(res, file.path(OUT, paste0(safe, "_chromvar_crossregion.csv")), row.names = FALSE)

  summ <- function(d, f) sprintf("FDR<0.05 %3d (up %3d, down %3d)", sum(f < 0.05), sum(f < 0.05 & d > 0), sum(f < 0.05 & d < 0))
  message(sprintf("[xr-motif] %s | plain: %s | centered: %s | weighted: %s | mean plain shift %+.2f",
                  ct, summ(res$diff_plain, res$FDR_plain), summ(res$diff_pooled, res$FDR),
                  summ(res$diff_weighted, res$FDR_weighted), res$global_shift_pooled[1]))
  message("[xr-motif] centered, P < 0.05: ", sum(res$P < 0.05), " (chance ~", round(0.05 * nrow(res)), ")")
  print(utils::head(res[, c("name", "diff_pooled", "P", "FDR", paste0("diff_", regs), "FDR_plain", "FDR_weighted")], 20),
        row.names = FALSE, digits = 2)
}
message("\n[xr-motif] done - ", OUT)
