# ==============================================================================
# Outlier-donor sensitivity of the microglial LDAM signature test (PWH vs PWoH).
#
# Reproduces the LDAM signature analysis of Junchen (revision/external/ldam_scores.R,
# not modified): microglia from the combined
# object, AggregateExpression counts per region x condition x donor, GSVA ssGSEA
# on the 4-gene LDAM set (NAMPT, ACSL1, DPYD, CD163), then per region a
# two-sample Welch t-test PWH vs PWoH (the ggpubr stat_compare_means default).
#
# Adds, per region: the same t-test without donors 1039 and 1181, and with each
# donor left out in turn. ssGSEA scores each sample from its own gene ranks; the
# only cross-sample step is normalize = TRUE, a single linear rescale by the
# score range, which cannot change a t-test. So scores are computed once on all
# donors and the exclusions only change which donors enter the test.
# As in the original analysis, every donor is kept regardless of its microglia
# count (reported per donor as n_micro).
#
# Input:  OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3_no1040_msn.RDS
# Output: OUT_ROOT/03-DE_analysis/donor_sensitivity/ldam_ssgsea/
#           ldam_ssgsea_scores.csv, ldam_ssgsea_sensitivity.csv
# Requires the RLIBS_CCC library (revision/tools/install_ccc_rlibs.R).
#
# Usage: revision/docker/run.sh revision/tools/ldam_ssgsea_sensitivity.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
.libPaths(c(RLIBS_CCC, .libPaths()))
suppressPackageStartupMessages({ library(Seurat); library(GSVA); library(dplyr); library(tidyr) })

OUT <- file.path(OUT_ROOT, "03-DE_analysis", "donor_sensitivity", "ldam_ssgsea")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PATHOLOGY_DONORS <- c("1039", "1181")

obj <- readRDS(file.path(OUT_ROOT, "01-Regional_data_processing", "combined", "data",
                         "seurat_combined_v3_no1040_msn.RDS"))
mic <- subset(obj, cells = colnames(obj)[obj@meta.data$consensus_annotation_v3 == "Micro"])
rm(obj); invisible(gc())
mic <- JoinLayers(mic)
n_micro <- mic@meta.data %>% count(region, sample_ids, name = "n_micro")

pc <- AggregateExpression(mic, assays = "RNA", slot = "counts",
                          group.by = c("region", "condition_new", "sample_ids"))
sets <- list(LDAM = c("NAMPT", "ACSL1", "DPYD", "CD163"))
sc <- gsva(ssgseaParam(exprData = pc$RNA, geneSets = sets))

d <- data.frame(id = colnames(sc), LDAM = as.numeric(sc["LDAM", ]), stringsAsFactors = FALSE) %>%
  separate(id, into = c("region", "condition", "donor"), sep = "_", remove = FALSE, extra = "merge") %>%
  mutate(condition = ifelse(condition == "HIV", "PWH", "PWoH")) %>%
  left_join(n_micro %>% mutate(donor = as.character(sample_ids)) %>% select(region, donor, n_micro),
            by = c("region", "donor"))
if (anyNA(d$n_micro)) message("[ldam] donor ids not matched to cell counts: ",
                              paste(d$id[is.na(d$n_micro)], collapse = ", "))
utils::write.csv(d, file.path(OUT, "ldam_ssgsea_scores.csv"), row.names = FALSE)

tt <- function(x) {
  if (length(unique(x$condition)) < 2 || min(table(x$condition)) < 2) return(c(NA, NA))
  t <- t.test(LDAM ~ condition, data = x)             # Welch, as stat_compare_means
  c(unname(diff(rev(t$estimate))), t$p.value)          # PWH - PWoH
}
res <- bind_rows(lapply(sort(unique(d$region)), function(R) {
  x <- d[d$region == R, ]
  full <- tt(x)
  present <- intersect(PATHOLOGY_DONORS, x$donor)
  ex <- tt(x[!x$donor %in% present, ])
  lo <- t(vapply(x$donor, function(dn) tt(x[x$donor != dn, ]), numeric(2)))
  w <- which.max(lo[, 2])
  data.frame(region = R, n_donors = nrow(x),
             n_PWH = sum(x$condition == "PWH"), n_PWoH = sum(x$condition == "PWoH"),
             diff_full = full[1], p_full = full[2],
             excluded = paste(present, collapse = ","), diff_excl = ex[1], p_excl = ex[2],
             lodo_n_sig = sum(lo[, 2] < 0.05), lodo_n_fits = nrow(lo),
             lodo_p_max = lo[w, 2], most_influential = rownames(lo)[w],
             lodo_diff_min = min(lo[, 1]), lodo_diff_max = max(lo[, 1]),
             stringsAsFactors = FALSE)
}))
utils::write.csv(res, file.path(OUT, "ldam_ssgsea_sensitivity.csv"), row.names = FALSE)
options(width = 200); print(res, row.names = FALSE, digits = 3)
message("[ldam] check p_full against the full-cohort LDAM_signature_boxplot.png (revision/external/ldam_scores.R) before reporting")
message("[ldam] done")
