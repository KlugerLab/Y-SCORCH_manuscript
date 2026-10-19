# Concordance of regional edgeR DE on uncorrected vs CellBender-corrected counts (ambient RNA correction, reviewer request).
#
# Method: two separate questions.
#   1. Concordance: do the same genes come out, in the same direction, with the
#      same magnitude? Overlap counts, Jaccard, sign agreement, and the logFC
#      correlation among gene x cell type calls significant (p_adj.loc < 0.05)
#      in either run.
#   2. Why uncorrected counts are the primary basis for DE: correction inflates
#      the extreme tail. Significant |logFC| > 5 calls roughly double under
#      correction in every region, a complete-separation artefact (genes driven
#      to zero in one group), not a biological gain.
#   For context, the ambient fraction of the typical gene is ~0.6-0.7%, which is
#   why concordance is high.
#   Significant-call counts move by -17% (PFC) and +17% (VST), so agreement holds
#   at the gene-identity level, not at the count level; the summary reports both.
#   Calls significant in only one run are counted with their p-value in the other
#   run, to show whether they were marginal, and a per-cell-type breakdown
#   (cell types with >= 20 calls) shows any cell type a region average would hide.
#
# Inputs:  OUT_ROOT/03-DE_analysis/<REGION>/results/<REGION>_HIV-CTR_res_edgeR_filtered_withFreq[_corrected].RDS
# Outputs: OUT_ROOT/03-DE_analysis/ambient_concordance/
#          concordance_by_region.csv, discordant_calls_marginality.csv,
#          extreme_tail_inflation.csv, concordance_by_celltype.csv,
#          all_significant_either_run.csv, logFC_concordance.{png,pdf}
#
# Usage: revision/docker/run.sh revision/tools/de_concordance_corrected.R

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

REGIONS <- c("INS", "PFC", "VST")
FDR     <- 0.05
OUTDIR  <- file.path(OUT_ROOT, "03-DE_analysis", "ambient_concordance")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## The paper's condition colours, reused so the figure matches the rest.
COL_BOTH <- "grey30"
COL_UNC  <- unname(PWH_COLOURS["PWH"])    # darkorange1
COL_COR  <- unname(PWH_COLOURS["PWoH"])   # #0072B2

.de_path <- function(region, corrected = FALSE)
  file.path(OUT_ROOT, "03-DE_analysis", region, "results",
            sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq%s.RDS",
                    region, if (corrected) "_corrected" else ""))

#' Read one DE result down to the columns needed, keyed on gene x cluster_id.
.read_de <- function(region, corrected) {
  p <- .de_path(region, corrected)
  if (!file.exists(p)) stop("missing DE result: ", p)
  d <- readRDS(p)
  d <- as.data.frame(d)
  need <- c("gene", "cluster_id", "logFC", "p_adj.loc")
  if (!all(need %in% colnames(d))) stop("unexpected columns in ", p)
  data.frame(key = paste(d$cluster_id, d$gene, sep = "|"),
             gene = d$gene, cluster_id = d$cluster_id,
             logFC = d$logFC, p_adj = d$`p_adj.loc`,
             stringsAsFactors = FALSE)
}

# ------------------------------------------------------------------------------
message("[conc] reading DE results")
all_joined <- list()

for (r in REGIONS) {
  u <- .read_de(r, FALSE); c_ <- .read_de(r, TRUE)
  message(sprintf("[conc] %s: uncorrected %d tests, corrected %d tests",
                  r, nrow(u), nrow(c_)))
  j <- merge(u, c_, by = c("key", "gene", "cluster_id"),
             suffixes = c("_unc", "_cor"), all = TRUE)
  j$region <- r
  all_joined[[r]] <- j
}
J <- dplyr::bind_rows(all_joined)

J$sig_unc <- !is.na(J$p_adj_unc) & J$p_adj_unc < FDR
J$sig_cor <- !is.na(J$p_adj_cor) & J$p_adj_cor < FDR
J$cat <- ifelse(J$sig_unc & J$sig_cor, "both",
         ifelse(J$sig_unc, "uncorrected only",
         ifelse(J$sig_cor, "corrected only", "neither")))

# ---- 1. region-level concordance ---------------------------------------------
region_tab <- J %>% group_by(region) %>% summarise(
  n_tested_both   = sum(!is.na(logFC_unc) & !is.na(logFC_cor)),
  n_sig_unc       = sum(sig_unc),
  n_sig_cor       = sum(sig_cor),
  pct_change_sig  = round(100 * (sum(sig_cor) - sum(sig_unc)) / sum(sig_unc), 1),
  n_both          = sum(cat == "both"),
  n_unc_only      = sum(cat == "uncorrected only"),
  n_cor_only      = sum(cat == "corrected only"),
  jaccard         = round(sum(cat == "both") /
                          sum(sig_unc | sig_cor), 3),
  pct_unc_sig_kept = round(100 * sum(cat == "both") / sum(sig_unc), 1),
  ## agreement among calls significant in BOTH - the claim that matters
  sign_agree_pct  = round(100 * mean(sign(logFC_unc[cat == "both"]) ==
                                     sign(logFC_cor[cat == "both"]), na.rm = TRUE), 2),
  pearson_both    = round(stats::cor(logFC_unc[cat == "both"],
                                     logFC_cor[cat == "both"],
                                     use = "complete.obs"), 3),
  spearman_both   = round(stats::cor(logFC_unc[cat == "both"],
                                     logFC_cor[cat == "both"],
                                     method = "spearman", use = "complete.obs"), 3),
  .groups = "drop")

cat("\n=== 1. CONCORDANCE BY REGION (FDR < ", FDR, ") ===\n", sep = "")
print(as.data.frame(region_tab), row.names = FALSE)

## Were the discordant calls marginal? If a gene significant in only one run sat
## just under the threshold in the other, the disagreement is threshold noise,
## not a changed conclusion.
marg <- J %>% filter(cat %in% c("uncorrected only", "corrected only")) %>%
  mutate(other_p = ifelse(cat == "uncorrected only", p_adj_cor, p_adj_unc)) %>%
  group_by(region, cat) %>%
  summarise(n = dplyr::n(),
            median_p_in_other_run = signif(median(other_p, na.rm = TRUE), 3),
            pct_other_p_under_0.10 = round(100 * mean(other_p < 0.10, na.rm = TRUE), 1),
            pct_absent_from_other  = round(100 * mean(is.na(other_p)), 1),
            .groups = "drop")
cat("\n=== 2. ARE THE DISCORDANT CALLS MARGINAL? ===\n")
print(as.data.frame(marg), row.names = FALSE)

# ---- 3. the justification: extreme-tail inflation -----------------------------
ext <- J %>% group_by(region) %>% summarise(
  extreme_unc = sum(abs(logFC_unc) > 5 & sig_unc, na.rm = TRUE),
  extreme_cor = sum(abs(logFC_cor) > 5 & sig_cor, na.rm = TRUE),
  fold        = round(sum(abs(logFC_cor) > 5 & sig_cor, na.rm = TRUE) /
                      pmax(1, sum(abs(logFC_unc) > 5 & sig_unc, na.rm = TRUE)), 2),
  median_abs_logFC_unc = round(median(abs(logFC_unc[sig_unc]), na.rm = TRUE), 3),
  median_abs_logFC_cor = round(median(abs(logFC_cor[sig_cor]), na.rm = TRUE), 3),
  .groups = "drop")
cat("\n=== 3. EXTREME-TAIL INFLATION UNDER CORRECTION ===\n")
cat("This is WHY uncorrected counts are primary. Correction multiplies the\n")
cat("number of runaway fold-changes without changing the bulk of the data.\n\n")
print(as.data.frame(ext), row.names = FALSE)

# ---- 4. per-cell-type, so a region average cannot hide a moved cell type ------
ct_tab <- J %>% group_by(region, cluster_id) %>%
  filter(sum(sig_unc | sig_cor) >= 20) %>%
  summarise(n_sig_unc = sum(sig_unc), n_sig_cor = sum(sig_cor),
            jaccard = round(sum(cat == "both") / sum(sig_unc | sig_cor), 3),
            pearson = round(stats::cor(logFC_unc[cat == "both"],
                                       logFC_cor[cat == "both"],
                                       use = "complete.obs"), 3),
            .groups = "drop") %>%
  arrange(jaccard)
cat("\n=== 4. PER CELL TYPE (lowest concordance first, >= 20 calls) ===\n")
print(as.data.frame(head(ct_tab, 15)), row.names = FALSE)

# ---- outputs -----------------------------------------------------------------
utils::write.csv(region_tab, file.path(OUTDIR, "concordance_by_region.csv"), row.names = FALSE)
utils::write.csv(marg,       file.path(OUTDIR, "discordant_calls_marginality.csv"), row.names = FALSE)
utils::write.csv(ext,        file.path(OUTDIR, "extreme_tail_inflation.csv"), row.names = FALSE)
utils::write.csv(ct_tab,     file.path(OUTDIR, "concordance_by_celltype.csv"), row.names = FALSE)
utils::write.csv(J[J$sig_unc | J$sig_cor,
                   c("region","cluster_id","gene","logFC_unc","p_adj_unc",
                     "logFC_cor","p_adj_cor","cat")],
                 file.path(OUTDIR, "all_significant_either_run.csv"), row.names = FALSE)

# ---- figure: logFC concordance, one panel per region -------------------------
pd <- J %>% filter(cat != "neither", !is.na(logFC_unc), !is.na(logFC_cor))
pd$cat <- factor(pd$cat, levels = c("both", "uncorrected only", "corrected only"))
lab <- region_tab %>%
  mutate(txt = sprintf("r = %.3f\n%d of %d kept (%.0f%%)",
                       pearson_both, n_both, n_sig_unc, pct_unc_sig_kept))

p <- ggplot(pd, aes(x = logFC_unc, y = logFC_cor, colour = cat)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 0.5, alpha = 0.35) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt),
            inherit.aes = FALSE, hjust = -0.08, vjust = 1.15, size = 3) +
  facet_wrap(~ region, nrow = 1) +
  scale_colour_manual(values = c("both" = COL_BOTH,
                                 "uncorrected only" = COL_UNC,
                                 "corrected only" = COL_COR), name = NULL) +
  labs(title = "Differential expression is concordant with and without ambient RNA correction",
       subtitle = sprintf(paste("pseudobulk edgeR, FDR < %.2f; each point is one gene x cell type",
                                "significant in either run.\nGrey line is y = x."), FDR),
       x = expression(log[2]*"FC, uncorrected counts"),
       y = expression(log[2]*"FC, CellBender-corrected counts")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank())

ggsave(file.path(OUTDIR, "logFC_concordance.png"), p, width = 10, height = 4.2, dpi = 300)
ggsave(file.path(OUTDIR, "logFC_concordance.pdf"), p, width = 10, height = 4.2)

message("\n[conc] wrote ", OUTDIR)
message("[conc] done")
