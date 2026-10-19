# ==============================================================================
# Per-region ambient-fraction tables, and the regional edgeR DE results
# annotated with an ambient flag for each gene.
#
# Method
#   ambient_fraction_table() (revision/tools/ambient_fraction.R) computes, per
#   gene and cell type, the fraction of signal CellBender attributed to ambient
#   RNA; annotate_de_with_ambient() joins it onto each DE call with a flag at
#   THRESHOLD (0.10). Both the uncorrected (primary) and the CellBender-corrected
#   DE results are annotated, and the flag rate among tested and significant
#   genes is summarised. The DE results themselves are not modified; annotated
#   copies are written alongside them. A cached ambient table is reused.
#
# Inputs
#   OUT_ROOT/03-DE_analysis/<REGION>/results/
#     <REGION>_HIV-CTR_res_edgeR_filtered_withFreq[_corrected].RDS
# Outputs
#   OUT_ROOT/03-DE_analysis/<REGION>/results/
#     <REGION>_ambient_fraction_by_celltype.csv                gene x cell type
#     <REGION>_HIV-CTR_res_edgeR_ambientFlagged.csv            uncorrected DE + flag
#     <REGION>_HIV-CTR_res_edgeR_ambientFlagged_corrected.csv  corrected DE + flag
#   OUT_ROOT/03-DE_analysis/ambient_flag_summary.csv           cross-region summary
#
# Usage
#   revision/docker/run.sh revision/tools/run_ambient_fraction.R [INS PFC VST]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/ambient_fraction.R"))

## Regions may be restricted from the command line, so INS can be inspected
## before committing to the other two:
##   revision/docker/run.sh revision/tools/run_ambient_fraction.R INS
.args <- commandArgs(trailingOnly = TRUE)
REGIONS <- if (length(.args)) toupper(.args) else c("INS", "PFC", "VST")
THRESHOLD <- 0.10
message("[amb] regions: ", paste(REGIONS, collapse = ", "))

de_path <- function(region, corrected = FALSE)
  file.path(OUT_ROOT, "03-DE_analysis", region, "results",
            sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq%s.RDS",
                    region, if (corrected) "_corrected" else ""))

summ <- list()

for (rg in REGIONS) {
  message("\n===================== ", rg, " =====================")

  ## reuse the cache if it is already there - the object load is the slow part
  cache <- ambient_table_path(rg)
  amb <- if (file.exists(cache)) {
    message("[amb] reusing ", cache)
    utils::read.csv(cache, stringsAsFactors = FALSE)
  } else ambient_fraction_table(rg)

  for (corrected in c(FALSE, TRUE)) {
    p <- de_path(rg, corrected)
    if (!file.exists(p)) { message("[amb] no DE result: ", p); next }
    de <- as.data.frame(readRDS(p))
    de <- annotate_de_with_ambient(de, amb, threshold = THRESHOLD)

    out <- sub("_filtered_withFreq", "_ambientFlagged", sub("\\.RDS$", ".csv", p))
    utils::write.csv(de, out, row.names = FALSE)
    message("[amb] wrote ", out)

    s <- ambient_flag_summary(de)
    s$region <- rg
    s$counts <- if (corrected) "corrected" else "uncorrected"
    summ[[length(summ) + 1]] <- s

    message("--- ", rg, " (", s$counts[1], ") flag rate, tested vs significant ---")
    print(s[, c("threshold", "pct_tested", "pct_significant",
                "n_significant_flagged", "n_significant", "enrichment")],
          row.names = FALSE, digits = 3)

    message("--- by cell type at ", THRESHOLD, " ---")
    print(as.data.frame(ambient_flag_by_celltype(de, THRESHOLD)),
          row.names = FALSE, digits = 3)

    rm(de); invisible(gc(FALSE))
  }
  rm(amb); invisible(gc(FALSE))
}

allsum <- dplyr::bind_rows(summ)
## A partial run must not overwrite the all-region summary with a subset.
f <- file.path(OUT_ROOT, "03-DE_analysis",
               if (length(.args)) paste0("ambient_flag_summary_",
                                         paste(REGIONS, collapse = "-"), ".csv")
               else "ambient_flag_summary.csv")
utils::write.csv(allsum, f, row.names = FALSE)
message("\n[amb] wrote ", f)
print(as.data.frame(allsum), row.names = FALSE, digits = 3)
