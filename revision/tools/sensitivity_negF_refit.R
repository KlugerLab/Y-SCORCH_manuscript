# ==============================================================================
# Negative-F refit of the regional edgeR DE for the outlier-donor sensitivity
# comparison (primary cohort vs. donors 1039 and 1181 excluded).
#
# Method
#   With 1039 and 1181 excluded, six primary DEGs return p = 1 with an unchanged
#   or larger fold change (e.g. PFC astrocyte SERPINA3, log2FC 6.40 -> 6.64).
#   This is the edgeR 4 default quasi-likelihood test failing on near-separated
#   genes (negative F), the same failure handled in the cross-region DE, not a
#   loss of the effect. pbds_negF_fallback() (revision/config/paths.R) keeps the
#   default QL test for every gene and refits only rows with F < 0 with
#   glmQLFit(legacy = TRUE), or glmLRT where the legacy fit is also negative.
#   It is applied to both runs, because the primary tables have failed tests of
#   their own (INS 47, PFC 19, VST 20) and fixing one side only would bias the
#   comparison toward "gained" calls.
#     exclusion  the saved pbDS object (10-sensitivity_DE.R keeps it unrounded,
#                with $data); the fallback is applied directly
#     primary    only resDS() tables were saved, so pbDS is rerun on the saved
#                pseudobulk and design with the published call; it must
#                reproduce the saved primary table (logFC, p_adj.glb) before the
#                fallback is applied, or the script stops
#   The regional primary results used in the paper are not changed; the outputs
#   are separate files used only for this comparison.
#
# Inputs   OUT_ROOT/03-DE_analysis/<R>/data/{pseudo_bulk.rds, mm.rds, mm_noHIVE_noPath.rds}
#          OUT_ROOT/03-DE_analysis/<R>/results/<R>_{HIV-CTR,noHIVE_noPath}_res_edgeR_filtered_withFreq.RDS
# Outputs  OUT_ROOT/03-DE_analysis/<R>/results/<R>_{primary,noHIVE_noPath}_negF_res_edgeR.RDS
#          (then compared with 03-DE_analysis/11-compare_sensitivity.R,
#          --out donor_sensitivity/<R>_noHIVE_noPath_negF)
#
# Usage
#   revision/docker/run.sh revision/tools/sensitivity_negF_refit.R [--regions INS,PFC,VST]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(muscat); library(dplyr); library(limma) })
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
REGIONS <- strsplit(getarg("--regions", "INS,PFC,VST"), ",")[[1]]
DE <- file.path(OUT_ROOT, "03-DE_analysis")

flat <- function(res) bind_rows(lapply(res$table[[1]], as.data.frame))
n_negF <- function(res) sum(flat(res)$F < 0, na.rm = TRUE)

for (R in REGIONS) {
  message("\n===== ", R)
  ## ---- exclusion run: fallback on the saved pbDS object -------------------------------------
  ex  <- readRDS(file.path(DE, R, "results", sprintf("%s_noHIVE_noPath_res_edgeR_filtered_withFreq.RDS", R)))
  if (is.null(ex$data)) stop(R, ": exclusion object has no $data - cannot refit")
  mmx <- readRDS(file.path(DE, R, "data", "mm_noHIVE_noPath.rds"))
  message("[negF] exclusion: ", n_negF(ex), " rows with F < 0")
  ## same contrast 10-sensitivity_DE.R built; its column name matches names(ex$table)
  cx  <- makeContrasts("HIV-CTR" = HIV - CTR, levels = mmx)
  stopifnot(identical(colnames(cx), names(ex$table)))
  ex2 <- pbds_negF_fallback(ex, design = mmx, contrast = cx)
  saveRDS(ex2, file.path(DE, R, "results", sprintf("%s_noHIVE_noPath_negF_res_edgeR.RDS", R)))

  ## ---- primary run: re-derive pbDS from the saved inputs, verify, then fallback -------------
  pb <- readRDS(file.path(DE, R, "data", "pseudo_bulk.rds"))
  mm <- readRDS(file.path(DE, R, "data", "mm.rds"))
  contrast <- makeContrasts("HIV-CTR", levels = mm)
  pr <- pbDS(pb, design = mm, contrast = contrast, filter = "both", min_cells = 10, method = "edgeR")
  saved <- as.data.frame(readRDS(file.path(DE, R, "results",
                                           sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq.RDS", R))))
  chk <- inner_join(flat(pr) %>% mutate(cluster_id = as.character(cluster_id)),
                    saved %>% mutate(cluster_id = as.character(cluster_id)) %>%
                      select(gene, cluster_id, lfc_s = logFC, padj_s = p_adj.glb),
                    by = c("gene", "cluster_id"))
  ## resDS() rounds to 3 significant digits, so the tolerance is relative for logFC:
  ## an absolute 0.01 would fail on rounding alone (a |logFC| ~ 10 rounded to 3
  ## digits differs by up to ~0.04). Rounding error is at most 0.5% of the value;
  ## 1% leaves margin and still catches any real refit change.
  rel <- abs(chk$logFC - chk$lfc_s) / pmax(abs(chk$lfc_s), 0.01)
  d_p <- max(abs(chk$p_adj.glb - chk$padj_s))
  w <- which.max(rel)
  message(sprintf(paste("[negF] primary refit vs saved table: %d of %d rows matched; max relative dlogFC %.2g",
                        "(worst: %s %s, refit %.4f vs saved %.4f); max |dFDR| %.2g"),
                  nrow(chk), nrow(saved), rel[w], chk$cluster_id[w], chk$gene[w],
                  chk$logFC[w], chk$lfc_s[w], d_p))
  if (nrow(chk) != nrow(saved) || rel[w] > 0.01 || d_p > 0.01)
    stop(R, ": re-derived primary fit does not reproduce the saved primary table")
  message("[negF] primary: ", n_negF(pr), " rows with F < 0")
  pr2 <- pbds_negF_fallback(pr, design = mm, contrast = contrast)
  saveRDS(pr2, file.path(DE, R, "results", sprintf("%s_primary_negF_res_edgeR.RDS", R)))

  ## rows that were refit in the exclusion run
  show <- flat(ex2) %>% filter(test != "QL") %>%
    select(cluster_id, gene, logFC, F, p_val, p_adj.glb, test) %>% arrange(p_adj.glb)
  message("[negF] exclusion rows refit (first 15 by FDR):")
  print(head(as.data.frame(show), 15), row.names = FALSE, digits = 3)
}
message("\n[negF] done")
