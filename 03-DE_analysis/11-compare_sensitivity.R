#!/usr/bin/env Rscript
# ==============================================================================
# Compare two pseudobulk DE results gene by gene, per cell type (revision).
#
# Method: each gene in each cell type is a DEG when p_adj.glb < --padj (0.05)
# and |logFC| > --logfc (0.5). Genes are classed as retained, flipped, lost,
# gained or ns_both between run A and run B. A gene that "loses" significance
# because the QL F-test degenerated (negative F, or p_adj = 1 with |logFC| > 1)
# is classed separately as lost_test_failed: that is not the same as an effect
# that went away, and conflating them would overstate fragility. Also reports
# logFC correlation over genes tested in both runs.
#
# Inputs:  two DE result RDS files, paths relative to OUT_ROOT/03-DE_analysis.
# Outputs: OUT_ROOT/03-DE_analysis/comparisons/<out>/gene_level_fates.csv,
#          celltype_summary.csv, logfc_concordance.png, fate_by_celltype.png.
#
# Usage: revision/docker/run.sh 03-DE_analysis/11-compare_sensitivity.R --a ... --b ...
#
#   VST exclusions (reviewer request: dependence of VST findings on exclusions):
#     --a VST/results/VST_published26_res_edgeR_filtered_withFreq.RDS --a-label published26
#     --b VST/results/VST_criterion29_res_edgeR_filtered_withFreq.RDS --b-label criterion29
#
#   HIVE leave-one-out:
#     --a INS/results/INS_HIV-CTR_res_edgeR_filtered_withFreq.RDS --a-label full
#     --b INS/results/INS_noHIVE_res_edgeR_filtered_withFreq.RDS    --b-label noHIVE
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(ggrepel); library(optparse)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--a",       type = "character"),
  make_option("--b",       type = "character"),
  make_option("--a-label", type = "character", default = "A"),
  make_option("--b-label", type = "character", default = "B"),
  make_option("--out",     type = "character", default = "sensitivity_comparison"),
  make_option("--padj",    type = "double",    default = 0.05),
  make_option("--logfc",   type = "double",    default = 0.5)
)))

A_LAB <- opt$`a-label`; B_LAB <- opt$`b-label`
rerun_setup("03-DE_analysis", subdirs = file.path("comparisons", opt$out))
OUT <- file.path("comparisons", opt$out)

read_res <- function(p) {
  if (!file.exists(p)) stop("missing DE result: ", p)
  r <- readRDS(p)
  # muscat pbDS results: $table[[contrast]] is a list of per-cluster data frames
  if (is.list(r) && "table" %in% names(r)) r <- r$table[[1]]
  if (is.list(r) && !is.data.frame(r)) r <- bind_rows(r)
  as.data.frame(r)
}
a <- read_res(opt$a); b <- read_res(opt$b)
cat(sprintf("%s: %d rows | %s: %d rows\n", A_LAB, nrow(a), B_LAB, nrow(b)))

## A QL F-test that degenerated, not a real loss of signal. p_adj pinned at 1
## with a large fold change, or a negative F, is the failure mode on
## near-separated genes.
flag_test_failure <- function(df) {
  df$test_failed <- !is.na(df$p_adj.glb) & df$p_adj.glb == 1 &
    !is.na(df$logFC) & abs(df$logFC) > 1
  if ("F" %in% colnames(df)) df$test_failed <- df$test_failed | (!is.na(df$F) & df$F < 0)
  df
}
tag_degs <- function(df) {
  df %>% mutate(
    is_deg = !is.na(p_adj.glb) & p_adj.glb < opt$padj & abs(logFC) > opt$logfc,
    direction = case_when(!is_deg ~ "ns", logFC > 0 ~ "up", TRUE ~ "down"))
}
a <- tag_degs(flag_test_failure(a)); b <- tag_degs(flag_test_failure(b))

joined <- full_join(
  a %>% select(cluster_id, gene, logFC_a = logFC, padj_a = p_adj.glb,
               is_deg_a = is_deg, direction_a = direction, tf_a = test_failed),
  b %>% select(cluster_id, gene, logFC_b = logFC, padj_b = p_adj.glb,
               is_deg_b = is_deg, direction_b = direction, tf_b = test_failed),
  by = c("cluster_id", "gene")) %>%
  mutate(across(c(is_deg_a, is_deg_b, tf_a, tf_b), ~ coalesce(.x, FALSE)),
         across(c(direction_a, direction_b), ~ coalesce(.x, "ns")),
         fate = case_when(
           is_deg_a &  is_deg_b & direction_a != direction_b ~ "flipped",
           is_deg_a &  is_deg_b                              ~ "retained",
           is_deg_a & !is_deg_b & tf_b                       ~ "lost_test_failed",
           is_deg_a & !is_deg_b                              ~ "lost",
          !is_deg_a &  is_deg_b                              ~ "gained",
           TRUE                                              ~ "ns_both"))

write.csv(joined, file.path(OUT, "gene_level_fates.csv"), row.names = FALSE)

summ <- joined %>% group_by(cluster_id) %>%
  summarise(n_deg_a = sum(is_deg_a), n_deg_b = sum(is_deg_b),
            retained = sum(fate == "retained"), lost = sum(fate == "lost"),
            lost_test_failed = sum(fate == "lost_test_failed"),
            gained = sum(fate == "gained"), flipped = sum(fate == "flipped"),
            .groups = "drop") %>%
  mutate(pct_retained = ifelse(n_deg_a > 0, round(100 * retained / n_deg_a, 1), NA_real_)) %>%
  arrange(desc(n_deg_a))
names(summ)[2:3] <- paste0("n_deg_", c(A_LAB, B_LAB))
write.csv(summ, file.path(OUT, "celltype_summary.csv"), row.names = FALSE)
print(as.data.frame(summ), row.names = FALSE)

## headline concordance
tot_a <- sum(joined$is_deg_a); tot_b <- sum(joined$is_deg_b)
both  <- sum(joined$fate %in% c("retained", "flipped"))
cat(sprintf("\nDEGs %s: %d | %s: %d | shared: %d (%.1f%% of %s)\n",
            A_LAB, tot_a, B_LAB, tot_b, both,
            ifelse(tot_a > 0, 100 * both / tot_a, NA), A_LAB))
cat(sprintf("flipped direction: %d   lost: %d   lost to test failure: %d   gained: %d\n",
            sum(joined$fate == "flipped"), sum(joined$fate == "lost"),
            sum(joined$fate == "lost_test_failed"), sum(joined$fate == "gained")))

## logFC concordance among genes tested in both — the robustness evidence
shared <- joined %>% filter(!is.na(logFC_a), !is.na(logFC_b))
if (nrow(shared) > 10) {
  r <- cor(shared$logFC_a, shared$logFC_b, use = "complete.obs")
  rd <- cor(shared$logFC_a, shared$logFC_b, method = "spearman", use = "complete.obs")
  cat(sprintf("logFC correlation over %d shared genes: pearson %.3f  spearman %.3f\n",
              nrow(shared), r, rd))
  deg_any <- shared %>% filter(is_deg_a | is_deg_b)
  p <- ggplot(shared, aes(logFC_a, logFC_b)) +
    geom_point(alpha = 0.15, size = 0.5, colour = "grey60") +
    geom_point(data = deg_any, aes(colour = fate), alpha = 0.8, size = 0.9) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30") +
    labs(title = sprintf("logFC concordance: %s vs %s", A_LAB, B_LAB),
         subtitle = sprintf("pearson r = %.3f over %d genes", r, nrow(shared)),
         x = sprintf("logFC (%s)", A_LAB), y = sprintf("logFC (%s)", B_LAB)) +
    theme_classic(base_size = 11)
  ggsave(file.path(OUT, "logfc_concordance.png"), p, width = 7.5, height = 6.5, dpi = 150)

  pb <- summ %>% select(cluster_id, retained, lost, lost_test_failed, gained, flipped) %>%
    pivot_longer(-cluster_id, names_to = "fate", values_to = "n")
  p2 <- ggplot(pb, aes(reorder(cluster_id, -n), n, fill = fate)) +
    geom_col() + coord_flip() +
    labs(title = sprintf("DEG fate by cell type: %s -> %s", A_LAB, B_LAB),
         x = NULL, y = "genes") + theme_classic(base_size = 11)
  ggsave(file.path(OUT, "fate_by_celltype.png"), p2, width = 8, height = 7, dpi = 150)
}
cat("\nwrote", file.path(getwd(), OUT), "\n")
