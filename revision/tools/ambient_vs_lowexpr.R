# ==============================================================================
# Overlap between ambient-flagged regional DE calls and low-detection calls.
#
# Method: two separate concerns apply to the regional edgeR DE:
#   (1) low expression - genes detected in <=5% of nuclei in both groups are still
#       tested and carry inflated |logFC| (INS: 30% of genes, 23-26% of significant
#       calls, median |logFC| 2.03 vs 1.06);
#   (2) ambient RNA - genes whose signal CellBender attributed largely to ambient
#       RNA (ambient_frac > 0.10; revision/tools/ambient_fraction.R).
# For all tested, significant (FDR < 0.05) and significant with |logFC| > 5 rows,
# this cross-tabulates the two flags (Fisher's exact test, Jaccard index) and
# summarises |logFC| per flag category. Largely overlapping flags mean one reporting
# filter handles both; disjoint flags are two independent caveats.
#
# Inputs   OUT_ROOT/03-DE_analysis/<REGION>/results/
#            <REGION>_HIV-CTR_res_edgeR_ambientFlagged[_corrected].csv
# Outputs  OUT_ROOT/03-DE_analysis/ambient_vs_lowexpr[_<REGIONS>].csv
#
# Usage: revision/docker/run.sh revision/tools/ambient_vs_lowexpr.R [INS PFC VST]
#        (no arguments = every region present)
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(dplyr) })

AMB_THRESHOLD <- 0.10
LOW_FRQ       <- 0.05     # "low in both" = detected in <=5% of nuclei in each group
FDR           <- 0.05
EXTREME_LFC   <- 5        # inflated |logFC| tail seen in low-expression genes

.args <- commandArgs(trailingOnly = TRUE)
REGIONS <- if (length(.args)) toupper(.args) else c("INS", "PFC", "VST")

flagged_path <- function(region, corrected = FALSE)
  file.path(OUT_ROOT, "03-DE_analysis", region, "results",
            sprintf("%s_HIV-CTR_res_edgeR_ambientFlagged%s.csv",
                    region, if (corrected) "_corrected" else ""))

#' 2x2 agreement between the two flags, plus Fisher's exact test.
.crosstab <- function(d, label) {
  a <- d$ambient_flag
  l <- d$low_in_both
  ok <- !is.na(a) & !is.na(l)
  a <- a[ok]; l <- l[ok]
  tab <- table(ambient = a, low_expr = l)
  ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)

  both <- sum(a & l); either <- sum(a | l)
  data.frame(
    set = label, n = length(a),
    n_ambient = sum(a), n_lowexpr = sum(l), n_both = both,
    pct_ambient_also_low = 100 * both / max(sum(a), 1),
    pct_low_also_ambient = 100 * both / max(sum(l), 1),
    jaccard = both / max(either, 1),
    odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
    p = if (is.null(ft)) NA_real_ else ft$p.value,
    row.names = NULL)
}

out <- list()

for (rg in REGIONS) {
  for (corrected in c(FALSE, TRUE)) {
    p <- flagged_path(rg, corrected)
    if (!file.exists(p)) { message("[x] not found, skipping: ", p); next }
    cnt <- if (corrected) "corrected" else "uncorrected"
    message("\n===== ", rg, " (", cnt, ") =====")

    d <- utils::read.csv(p, stringsAsFactors = FALSE)
    need <- c("ambient_frac", "avg_CTR_frq", "avg_HIV_frq", "p_adj.glb", "logFC")
    miss <- setdiff(need, colnames(d))
    if (length(miss)) { message("[x] lacks: ", paste(miss, collapse = ", ")); next }

    d$ambient_flag <- d$ambient_frac > AMB_THRESHOLD
    d$low_in_both  <- d$avg_CTR_frq <= LOW_FRQ & d$avg_HIV_frq <= LOW_FRQ
    d$sig          <- !is.na(d$p_adj.glb) & d$p_adj.glb < FDR
    d$extreme      <- abs(d$logFC) > EXTREME_LFC

    r <- dplyr::bind_rows(
      .crosstab(d, "all tested"),
      .crosstab(d[d$sig, , drop = FALSE], "significant"),
      .crosstab(d[d$sig & d$extreme, , drop = FALSE], "significant & |logFC|>5"))
    r$region <- rg; r$counts <- cnt
    out[[length(out) + 1]] <- r
    print(r[, c("set", "n", "n_ambient", "n_lowexpr", "n_both",
                "pct_ambient_also_low", "pct_low_also_ambient",
                "jaccard", "odds_ratio", "p")],
          row.names = FALSE, digits = 3)

    ## How much of the significant set does each flag account for, alone and together?
    s <- d[d$sig, , drop = FALSE]
    message(sprintf(
      "[x] of %d significant calls: %d (%.1f%%) ambient-only, %d (%.1f%%) low-only, %d (%.1f%%) both, %d (%.1f%%) neither",
      nrow(s),
      sum(s$ambient_flag & !s$low_in_both, na.rm = TRUE),
      100 * mean(s$ambient_flag & !s$low_in_both, na.rm = TRUE),
      sum(!s$ambient_flag & s$low_in_both, na.rm = TRUE),
      100 * mean(!s$ambient_flag & s$low_in_both, na.rm = TRUE),
      sum(s$ambient_flag & s$low_in_both, na.rm = TRUE),
      100 * mean(s$ambient_flag & s$low_in_both, na.rm = TRUE),
      sum(!s$ambient_flag & !s$low_in_both, na.rm = TRUE),
      100 * mean(!s$ambient_flag & !s$low_in_both, na.rm = TRUE)))

    ## |logFC| by category - the practical question is whether either flag
    ## picks out the inflated effect sizes.
    lf <- s %>%
      dplyr::mutate(cat = dplyr::case_when(
        ambient_flag &  low_in_both ~ "both",
        ambient_flag & !low_in_both ~ "ambient only",
        !ambient_flag &  low_in_both ~ "low-expr only",
        TRUE ~ "neither")) %>%
      dplyr::group_by(cat) %>%
      dplyr::summarise(n = dplyr::n(),
                       median_abs_logFC = stats::median(abs(logFC), na.rm = TRUE),
                       pct_extreme = 100 * mean(abs(logFC) > EXTREME_LFC, na.rm = TRUE),
                       .groups = "drop")
    message("[x] |logFC| by category:")
    print(as.data.frame(lf), row.names = FALSE, digits = 3)
  }
}

if (length(out)) {
  res <- dplyr::bind_rows(out)
  f <- file.path(OUT_ROOT, "03-DE_analysis",
                 if (length(.args)) paste0("ambient_vs_lowexpr_",
                                           paste(REGIONS, collapse = "-"), ".csv")
                 else "ambient_vs_lowexpr.csv")
  utils::write.csv(res, f, row.names = FALSE)
  message("\n[x] wrote ", f)
} else message("\n[x] nothing to summarise")
