# ==============================================================================
# Donor-level tests of the microglial state signatures: per region (primary)
# and pooled across regions, with a camera summary table.
#
# Method
#   Driver for microglia_score_all() in revision/tools/microglia_signature_score.R:
#   donor pseudobulk of uncorrected counts within CellBender-defined microglia,
#   per-donor signature scores, and camera (competitive) and fry
#   (self-contained) tests with n = donors. Each set is also tested with
#   ambient-suspect genes removed (fdr_drop). The summary prints sets at
#   FDR < 0.1 in any region.
#
# Inputs   the three regional objects, loaded one at a time (~2.4-3 GB each)
# Outputs  OUT_ROOT/03-DE_analysis/microglia_signatures/
#            microglia_signature_results.rds
#            camera_summary_all_regions.csv
#
# Usage
#   revision/docker/run.sh revision/tools/run_microglia_signatures.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/microglia_signature_score.R"))

res <- microglia_score_all(regions = c("INS", "PFC", "VST"))

saveRDS(res, file.path(MICRO_OUT, "microglia_signature_results.rds"))
message("\n[mgsig] wrote ", file.path(MICRO_OUT, "microglia_signature_results.rds"))

## ---- headline table: camera, full vs drop-ambient, all regions -------------
hd <- do.call(rbind, lapply(names(res$per_region), function(r) {
  x <- res$per_region[[r]]
  cf <- x$tests$camera; cd <- x$tests_drop$camera
  data.frame(region = r, signature = cf$signature,
             n_genes = cf$NGenes, direction = cf$Direction,
             p_full = cf$PValue, fdr_full = cf$FDR,
             fdr_drop = cd$FDR[match(cf$signature, cd$signature)],
             stringsAsFactors = FALSE)
}))
hd <- hd[order(hd$fdr_full), ]
utils::write.csv(hd, file.path(MICRO_OUT, "camera_summary_all_regions.csv"),
                 row.names = FALSE)

message("\n########## CAMERA SUMMARY (FDR < 0.1 in any region) ##########")
sig <- hd[!is.na(hd$fdr_full) & hd$fdr_full < 0.1, ]
if (nrow(sig)) print(sig, row.names = FALSE) else
  message("nothing at FDR < 0.1 - report that plainly rather than lowering the bar")

message("\n[mgsig] done")
