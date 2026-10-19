# Compares the donor-naive edgeR and donor-aware cross-region DE, both judged by the paper's full DEG rule.
#
# Method: the comparison written by revision/tools/cross_region_donor_aware.R
#   (donor_aware/edgeR_vs_donor_aware.csv) calls a gene significant on adjusted p
#   alone. The paper's DEG definition also requires |logFC| > 0.25 and expression
#   in > 5% of nuclei in either group; this applies the full rule (FDR < 0.05,
#   |logFC| > 0.25, > 5% expressed) to both the donor-aware limma-voom +
#   duplicateCorrelation fit and the donor-naive edgeR fit, so they are compared
#   like with like. The expression filter is a property of the gene in that cell
#   type, not of the model, so the edgeR table's avg_CTR_frq / avg_HIV_frq
#   columns are used for both.
#
# Inputs:  OUT_ROOT/03-DE_analysis/cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#          OUT_ROOT/03-DE_analysis/cross_region/donor_aware/cross_region_donor_aware_DE.csv
# Outputs: OUT_ROOT/03-DE_analysis/cross_region/donor_aware/
#            paper_rule_comparison_by_celltype.csv
#            paper_rule_comparison_overall.csv
#
# Usage: revision/docker/run.sh revision/tools/cross_region_paper_rule_compare.R

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(dplyr) })

D   <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region")
FDR <- 0.05; LFC <- 0.25; FRQ <- 0.05

ed <- readRDS(file.path(D, "data", "all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")) %>%
  ungroup() %>%
  transmute(cluster_id = as.character(cluster_id), gene,
            logFC_edgeR = logFC, p_adj.glb_edgeR = p_adj.glb, avg_CTR_frq, avg_HIV_frq)
da <- utils::read.csv(file.path(D, "donor_aware", "cross_region_donor_aware_DE.csv"),
                      stringsAsFactors = FALSE) %>%
  transmute(cluster_id, gene, logFC_da = logFC, p_adj.glb_da = p_adj.glb)

x <- dplyr::inner_join(ed, da, by = c("cluster_id", "gene"))
message("[cmp] genes tested by both models: ", nrow(x), " across ",
        length(unique(x$cluster_id)), " cell types")

expressed <- (x$avg_CTR_frq > FRQ) | (x$avg_HIV_frq > FRQ)
x$deg_edgeR <- expressed & !is.na(x$p_adj.glb_edgeR) & x$p_adj.glb_edgeR < FDR &
               abs(x$logFC_edgeR) > LFC
x$deg_da    <- expressed & !is.na(x$p_adj.glb_da) & x$p_adj.glb_da < FDR &
               abs(x$logFC_da) > LFC

summ <- function(d) data.frame(
  degs_edgeR          = sum(d$deg_edgeR),
  degs_donor_aware    = sum(d$deg_da),
  both                = sum(d$deg_edgeR & d$deg_da),
  pct_edgeR_retained  = round(100 * sum(d$deg_edgeR & d$deg_da) / max(1, sum(d$deg_edgeR)), 1),
  donor_aware_only    = sum(d$deg_da & !d$deg_edgeR),
  sign_agreement_pct  = round(100 * mean(sign(d$logFC_edgeR[d$deg_edgeR & d$deg_da]) ==
                                         sign(d$logFC_da[d$deg_edgeR & d$deg_da])), 1),
  logFC_r_among_edgeR_degs = round(stats::cor(d$logFC_edgeR[d$deg_edgeR], d$logFC_da[d$deg_edgeR]), 3))

by_ct <- x %>% group_by(cluster_id) %>% group_modify(~ summ(.x)) %>% ungroup() %>%
  arrange(desc(degs_edgeR))
overall <- summ(x)
utils::write.csv(by_ct,   file.path(D, "donor_aware", "paper_rule_comparison_by_celltype.csv"), row.names = FALSE)
utils::write.csv(overall, file.path(D, "donor_aware", "paper_rule_comparison_overall.csv"), row.names = FALSE)
message("\n[cmp] overall (paper DEG rule applied to both models):")
print(overall, row.names = FALSE)
message("\n[cmp] per cell type:")
print(as.data.frame(by_ct), row.names = FALSE)
