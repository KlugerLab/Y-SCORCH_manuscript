# Tabulates donor representation (the effective n) per region x cell type for the regional pseudobulk DE.
#
# Method: reviewer request for the unit of inference, minimum donor/cell
#   thresholds per cell type and region, and donor representation. The regional
#   edgeR DE calls
#     pbDS(pb, design = mm, contrast = contrast, filter = "both",
#          min_cells = 10, method = "edgeR")
#   so a donor contributes to a cell type's test only if it has >= 10 nuclei of
#   that type; the effective n varies by cell type and for rare types falls far
#   below the cohort size. This counts, from the saved pseudobulk objects, how
#   many donors clear the threshold per region x cell type, by condition, with
#   the cell counts behind them.
#   pbDS() also drops a cell type silently when the design loses rank for it
#   (muscat: `qr(d)$rank < ncol(d)` -> return(NULL)), so a cell type with too few
#   donors in one arm reads as "no differential expression" rather than "not
#   estimable". The `estimable` column (>= 2 passing donors per arm) flags these.
#
# Inputs:  OUT_ROOT/03-DE_analysis/<REGION>/data/pseudo_bulk.rds
# Outputs: OUT_ROOT/03-DE_analysis/donor_representation/
#          donor_representation_by_celltype.csv, donor_representation_summary.csv,
#          cells_per_donor_per_celltype.csv
#
# Usage: revision/docker/run.sh revision/tools/donor_representation.R

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(SingleCellExperiment); library(dplyr)
})

MIN_CELLS <- 10        # pbDS(min_cells = 10) in 03-DE_analysis/0[1-3]-*_edgeR.R
REGIONS   <- c("INS", "PFC", "VST")
OUTDIR    <- file.path(OUT_ROOT, "03-DE_analysis", "donor_representation")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

#' Long table of cells per donor x cell type for one region's pseudobulk.
.rep_one <- function(region, pb_file = "pseudo_bulk.rds") {
  p <- file.path(OUT_ROOT, "03-DE_analysis", region, "data", pb_file)
  if (!file.exists(p)) { warning("no pseudobulk for ", region); return(NULL) }
  pb <- readRDS(p)
  nc <- int_colData(pb)$n_cells          # named list: sample -> named int vector
  ei <- S4Vectors::metadata(pb)$experiment_info

  long <- do.call(rbind, lapply(names(nc), function(s) {
    v <- nc[[s]]
    data.frame(region = region, donor = s, cell_type = names(v),
               n_cells = as.integer(v), stringsAsFactors = FALSE)
  }))
  long$condition <- as.character(pwh_relabel(
    ei$group_id[match(long$donor, ei$sample_id)]))
  long
}

message("[rep] reading pseudobulk objects")
LONG <- dplyr::bind_rows(lapply(REGIONS, .rep_one))
if (!nrow(LONG)) stop("no pseudobulk objects found")

## ---- per region x cell type: effective n --------------------------------
tab <- LONG %>%
  group_by(region, cell_type) %>%
  summarise(
    donors_total   = dplyr::n(),
    PWH_pass       = sum(n_cells >= MIN_CELLS & condition == "PWH"),
    PWoH_pass      = sum(n_cells >= MIN_CELLS & condition == "PWoH"),
    donors_pass    = PWH_pass + PWoH_pass,
    pct_donors_pass = round(100 * (PWH_pass + PWoH_pass) / dplyr::n(), 1),
    median_cells   = median(n_cells[n_cells >= MIN_CELLS]),
    total_cells    = sum(n_cells),
    ## A contrast needs both arms populated. Two donors per arm is already
    ## fragile; below that the design will not support the covariates.
    estimable      = PWH_pass >= 2 & PWoH_pass >= 2,
    .groups = "drop") %>%
  arrange(region, desc(donors_pass))

cat("\n=== DONOR REPRESENTATION PER REGION x CELL TYPE (min_cells = ", MIN_CELLS, ") ===\n", sep = "")
cat("donors_pass is the EFFECTIVE n for that cell type's test, not the cohort size.\n\n")
for (r in REGIONS) {
  cat("---", r, "---\n")
  print(as.data.frame(tab[tab$region == r,
        c("cell_type","PWH_pass","PWoH_pass","donors_pass","pct_donors_pass",
          "median_cells","total_cells","estimable")]), row.names = FALSE)
  cat("\n")
}

## ---- summary per region -------------------------------------------------
summ <- tab %>% group_by(region) %>% summarise(
  n_cell_types      = dplyr::n(),
  cohort_donors     = max(donors_total),
  types_all_donors  = sum(donors_pass == donors_total),
  types_ge_20       = sum(donors_pass >= 20),
  types_lt_10       = sum(donors_pass < 10),
  types_not_estimable = sum(!estimable),
  median_effective_n = median(donors_pass),
  .groups = "drop")
cat("=== SUMMARY ===\n")
print(as.data.frame(summ), row.names = FALSE)

cat("\n=== CELL TYPES AT RISK (fewer than 2 donors in an arm) ===\n")
risk <- tab[!tab$estimable, c("region","cell_type","PWH_pass","PWoH_pass","total_cells")]
if (nrow(risk)) print(as.data.frame(risk), row.names = FALSE) else
  cat("none - every cell type has >= 2 donors per arm\n")

utils::write.csv(tab,  file.path(OUTDIR, "donor_representation_by_celltype.csv"), row.names = FALSE)
utils::write.csv(summ, file.path(OUTDIR, "donor_representation_summary.csv"), row.names = FALSE)
utils::write.csv(LONG, file.path(OUTDIR, "cells_per_donor_per_celltype.csv"), row.names = FALSE)
message("\n[rep] wrote ", OUTDIR)
message("[rep] done")
