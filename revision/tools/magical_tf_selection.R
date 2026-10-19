# ==============================================================================
# MAGICAL astrocyte GRN: per-TF eligibility tables for the figure TFs, and the
# TFs that link the reactive markers (SERPINA3, OSMR, CHI3L1).
#
# Method: per region, from the expanded MAGICAL circuits, a TF meets the
# RNA-based rule if
#   1 its own RNA is up in PWH astrocytes: logFC > 0 and p_adj.glb < 0.05 in the
#     cross-region edgeR table (the negative-F fallback result), and
#   2 it has at least --min-up target genes that are up-regulated in that region's
#     own astrocyte DE (logFC > 0, p_adj.glb < 0.05).
# The same rule is also evaluated on the donor-aware cross-region model (limma-voom
# + duplicateCorrelation, block = donor), the primary cross-region model; IRF1 and
# MAFG do not pass there (p_adj.glb 0.33 and 0.20). Column `meets_rule_donor_aware`.
# The TF selection used in the paper is the rule_activity criterion - TF
# expressed in >= 5% of nuclei, >= 5 up-regulated targets, chromVAR motif
# activity up in the donor-aware test - which is not computed by this script.
# Detection frequency (expr_frq) is the PWH frequency (avg_HIV_frq) of the
# cross-region table, used for the reactive-marker rule (>= --min-expr).
# Existing tables in the output folders are compared against the recomputed ones.
#
# Inputs: OUT_ROOT/magical/analyzed_output/<REGION>_new_<celltype>/expanded_circuits.csv
#         OUT_ROOT/03-DE_analysis/cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#         OUT_ROOT/03-DE_analysis/cross_region/donor_aware/cross_region_donor_aware_DE.csv
#         OUT_ROOT/03-DE_analysis/<REGION>/results/<REGION>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
# Output (with --write), per region in OUT_ROOT/magical/analyzed_output/<REGION>_new_<celltype>/:
#           tf_selection_table.csv     TF x (targets, up targets, own DE, rule flags)
#           reactive_marker_TFs.csv    TFs linked to SERPINA3 / OSMR / CHI3L1
#         and a combined OUT_ROOT/magical/analyzed_output/tf_selection_summary.csv
#
# Usage: revision/docker/run.sh revision/tools/magical_tf_selection.R [--celltype Astro]
#          [--min-up 5] [--min-expr 0.05] [--markers SERPINA3,OSMR,CHI3L1] [--write]
#        Without --write it only compares with the stored tables and writes nothing.
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(dplyr) })

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
CT       <- getarg("--celltype", "Astro")
MIN_UP   <- as.integer(getarg("--min-up", "5"))
MIN_EXPR <- as.numeric(getarg("--min-expr", "0.05"))
MARKERS  <- trimws(strsplit(getarg("--markers", "SERPINA3,OSMR,CHI3L1"), ",")[[1]])
WRITE    <- "--write" %in% args
PUBLISHED <- c("STAT3", "STAT1", "FOSL2", "FOS", "JUNB", "IRF1", "KLF4")
DE  <- file.path(OUT_ROOT, "03-DE_analysis")
AN  <- file.path(OUT_ROOT, "magical", "analyzed_output")

## ---- the TF's own differential expression, cross-region -----------------------------
xr <- as.data.frame(readRDS(file.path(DE, "cross_region", "data", "all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")))
xr <- xr[as.character(xr$cluster_id) == CT, ]
own <- data.frame(TF = xr$gene, logFC = xr$logFC, p_adj.glb = xr$p_adj.glb,
                  expr_frq = xr$avg_HIV_frq, stringsAsFactors = FALSE)
da <- utils::read.csv(file.path(DE, "cross_region", "donor_aware", "cross_region_donor_aware_DE.csv"),
                      stringsAsFactors = FALSE)
da <- da[da$cluster_id == CT, c("gene", "logFC", "p_adj.glb")]
names(da) <- c("TF", "logFC_donor_aware", "p_adj.glb_donor_aware")
own <- dplyr::left_join(own, da, by = "TF")

same <- function(a, b, tol = 1e-8) {
  if (is.numeric(a) && is.numeric(b)) all(abs(a - b) <= tol | (is.na(a) & is.na(b))) else all(a == b | (is.na(a) & is.na(b)))
}
summ <- list()
for (R in c("INS", "PFC", "VST")) {
  folder <- file.path(AN, magical_tag(R, CT))
  f_circ <- file.path(folder, "expanded_circuits.csv")
  if (!file.exists(f_circ)) { message("[tfsel] ", R, ": no expanded_circuits.csv - skipped"); next }
  ec <- utils::read.csv(f_circ, stringsAsFactors = FALSE)

  ## up-regulated targets come from the REGION's own DE, not the cross-region table
  reg <- as.data.frame(readRDS(file.path(DE, R, "results", sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq.RDS", R))))
  reg <- reg[as.character(reg$cluster_id) == CT, ]
  up_genes <- unique(reg$gene[reg$logFC > 0 & reg$p_adj.glb < 0.05])

  t <- ec %>% group_by(TF = TFs) %>%
    summarise(up_targets = length(intersect(unique(Gene_symbol), up_genes)),
              targets    = dplyr::n_distinct(Gene_symbol), .groups = "drop") %>%
    arrange(desc(up_targets)) %>%
    mutate(rank_up_targets = seq_len(dplyr::n()))
  t$SERPINA3 <- t$TF %in% ec$TFs[ec$Gene_symbol == "SERPINA3"]
  t <- dplyr::left_join(t, own, by = "TF")
  t$TF_up_in_PWH <- !is.na(t$p_adj.glb) & t$logFC > 0 & t$p_adj.glb < 0.05
  t$meets_rule   <- t$TF_up_in_PWH & t$up_targets >= MIN_UP
  t$meets_rule_donor_aware <- !is.na(t$p_adj.glb_donor_aware) & t$logFC_donor_aware > 0 &
                              t$p_adj.glb_donor_aware < 0.05 & t$up_targets >= MIN_UP
  t$published_set <- t$TF %in% PUBLISHED
  t <- t[, c("TF", "rank_up_targets", "up_targets", "targets", "SERPINA3", "logFC", "p_adj.glb",
             "expr_frq", "TF_up_in_PWH", "meets_rule", "published_set",
             "logFC_donor_aware", "p_adj.glb_donor_aware", "meets_rule_donor_aware")]

  ## TFs that bind a peak linked to one of the reactive markers
  mk <- ec[ec$Gene_symbol %in% MARKERS, ]
  rm_tab <- mk %>% group_by(TF = TFs) %>%
    summarise(n_genes = dplyr::n_distinct(Gene_symbol),
              genes   = paste(sort(unique(Gene_symbol)), collapse = "+"),
              max_prob = max(binding_prob), .groups = "drop") %>%
    dplyr::left_join(own %>% dplyr::select(TF, TF_logFC = logFC, TF_padj = p_adj.glb, TF_expr = expr_frq), by = "TF") %>%
    mutate(expressed = !is.na(TF_expr) & TF_expr >= MIN_EXPR,
           TF_up = !is.na(TF_padj) & TF_logFC > 0 & TF_padj < 0.05) %>%
    arrange(desc(n_genes), desc(expressed), desc(max_prob))

  ## ---- check against the stored tables ------------------------------------------------
  chk <- c()
  f_old <- file.path(folder, "tf_selection_table.csv")
  if (file.exists(f_old)) {
    o <- utils::read.csv(f_old, stringsAsFactors = FALSE)
    j <- dplyr::inner_join(o, t, by = "TF", suffix = c("_old", "_new"))
    ok <- nrow(j) == nrow(o) && same(j$up_targets_old, j$up_targets_new) && same(j$targets_old, j$targets_new) &&
          same(j$meets_rule_old, j$meets_rule_new) && same(round(j$logFC_old, 2), round(j$logFC_new, 2))
    chk <- c(chk, sprintf("tf_selection_table: %s (%d TFs)", if (ok) "reproduced" else "DIFFERS", nrow(o)))
  }
  f_old2 <- file.path(folder, "reactive_marker_TFs.csv")
  if (file.exists(f_old2)) {
    o2 <- utils::read.csv(f_old2, stringsAsFactors = FALSE)
    j2 <- dplyr::inner_join(o2, rm_tab, by = "TF", suffix = c("_old", "_new"))
    ## the stored file rounded max_prob to 3 decimals
    ok2 <- nrow(j2) == nrow(o2) && same(j2$genes_old, j2$genes_new) && same(j2$n_genes_old, j2$n_genes_new) &&
           same(round(j2$max_prob_old, 3), round(j2$max_prob_new, 3))
    chk <- c(chk, sprintf("reactive_marker_TFs: %s (%d TFs)", if (ok2) "reproduced" else "DIFFERS", nrow(o2)))
  }
  message(sprintf("[tfsel] %s: %d TFs in circuits; meet the rule %d (edgeR) / %d (donor-aware). %s",
                  R, nrow(t), sum(t$meets_rule), sum(t$meets_rule_donor_aware), paste(chk, collapse = "; ")))
  drop <- t$TF[t$meets_rule & !t$meets_rule_donor_aware]
  if (length(drop)) message("[tfsel]   eligible under edgeR but NOT donor-aware: ", paste(drop, collapse = ", "))
  if (WRITE) {
    utils::write.csv(t, f_old, row.names = FALSE)
    utils::write.csv(rm_tab, f_old2, row.names = FALSE)
  }
  summ[[R]] <- t %>% mutate(region = R)
}
S <- dplyr::bind_rows(summ)
if (WRITE && nrow(S)) utils::write.csv(S, file.path(AN, "tf_selection_summary.csv"), row.names = FALSE)
if (nrow(S)) {
  cat("\nTFs drawn in the figure (rule-based set), by region:\n")
  print(as.data.frame(S %>% dplyr::filter(TF %in% c("STAT3", "STAT1", "IRF1", "MAFG", "CEBPD")) %>%
    dplyr::select(region, TF, up_targets, targets, logFC, p_adj.glb, meets_rule,
                  logFC_donor_aware, p_adj.glb_donor_aware, meets_rule_donor_aware)), row.names = FALSE, digits = 3)
}
message("[tfsel] done", if (!WRITE) " (no files written; pass --write)" else "")
