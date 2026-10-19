#!/usr/bin/env Rscript
# ==============================================================================
# repair_fixed_meta.R: regenerates corrupt percent.mt metadata files.
#
# FIXED_META_DIR holds labeled_metadata_with_pt_mt_<sample>.csv for all 31/31/33
# samples, produced by the pipeline's create_fixed_metadata.R. Two of the 95 files
# contain embedded NUL bytes and are unreadable:
#
#   ins/1119_HIV_INS      44 NUL bytes, first at byte 6291560  (~6 MiB)
#   ins/HCTMQ_CTR_INS   2240 NUL bytes, first at byte 2097547  (~2 MiB)
#
# Both begin at a MiB boundary, which points to an interrupted write rather than
# anything wrong with the underlying data. In R they surface as
#   "EOF within quoted string" / "embedded nul(s) found in input"
# and then as a barcode mismatch, because read.csv silently returns short.
#
# That directory is read-only, so this regenerates the affected samples into
# REPAIRED_META_DIR using the same recipe as the original script: percent.mt from
# the RNA assay of raw_<sample>.rds (CellBender-corrected counts), joined by
# barcode onto labeled_metadata_<sample>.csv. fixed_meta_csv() in
# revision/config/paths.R prefers a repaired file when one exists.
#
# Scans every sample, so re-running it also catches any file corrupted later.
#
# Inputs:  FIXED_META_DIR, and PIPELINE_DIR/<region>/<sample>/ (labeled metadata
#          and raw_<sample>.rds).
# Outputs: REPAIRED_META_DIR (OUT_ROOT/meta/fixed_meta_repaired/<region>/<sample>/...).
#
# Usage: revision/docker/run.sh revision/config/repair_fixed_meta.R [ins pfc vst]
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

REGIONS <- c("ins", "pfc", "vst")
ONLY    <- commandArgs(trailingOnly = TRUE)   # optional: region names to limit to
if (length(ONLY)) REGIONS <- tolower(ONLY)

has_nul <- function(path) {
  n <- file.info(path)$size
  if (is.na(n) || n == 0) return(TRUE)
  any(readBin(path, what = "raw", n = n) == as.raw(0))
}

dir.create(REPAIRED_META_DIR, recursive = TRUE, showWarnings = FALSE)

repaired <- 0L
checked  <- 0L

for (region in REGIONS) {
  samples <- basename(list.dirs(region_dir(region), recursive = FALSE))
  for (sample in samples) {
    checked <- checked + 1L
    src <- file.path(FIXED_META_DIR, region, sample, "label_transfer",
                     paste0("labeled_metadata_with_pt_mt_", sample, ".csv"))

    ok <- file.exists(src) && !has_nul(src)
    if (ok) next

    cat("REPAIRING", region, "/", sample,
        if (!file.exists(src)) "(missing)" else "(NUL bytes)", "\n")

    labeled_csv <- file.path(region_dir(region), sample, "label_transfer",
                             paste0("labeled_metadata_", sample, ".csv"))
    raw_rds     <- file.path(region_dir(region), sample,
                             paste0("raw_", sample, ".rds"))
    if (!file.exists(labeled_csv)) stop("no labeled metadata: ", labeled_csv)
    if (!file.exists(raw_rds))     stop("no raw object: ", raw_rds)

    # percent.mt from the RNA assay of the raw object — the corrected counts,
    # matching what create_fixed_metadata.R did.
    raw_obj <- readRDS(raw_rds)
    if (!"RNA" %in% Assays(raw_obj)) stop(sample, ": no RNA assay in raw object")
    DefaultAssay(raw_obj) <- "RNA"
    raw_obj <- PercentageFeatureSet(raw_obj, pattern = "^MT-", col.name = "percent.mt")
    pct <- data.frame(barcode    = colnames(raw_obj),
                      percent.mt = raw_obj$percent.mt,
                      stringsAsFactors = FALSE)
    rm(raw_obj); gc(verbose = FALSE)

    labeled <- read.csv(labeled_csv, row.names = 1, stringsAsFactors = FALSE)
    labeled$barcode <- rownames(labeled)

    joined <- labeled %>% left_join(pct, by = "barcode")
    rownames(joined) <- joined$barcode

    n_missing <- sum(is.na(joined$percent.mt))
    if (n_missing > 0) {
      stop(sample, ": ", n_missing, " of ", nrow(joined),
           " barcodes absent from raw_<sample>.rds — not a corruption problem")
    }

    # Written barcode-first so read.csv(row.names = 1) recovers the barcodes,
    # matching how the original files are consumed.
    joined <- joined[, c("barcode", setdiff(names(joined), "barcode"))]

    out_dir <- file.path(REPAIRED_META_DIR, region, sample, "label_transfer")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out <- file.path(out_dir, paste0("labeled_metadata_with_pt_mt_", sample, ".csv"))
    write.csv(joined, out, row.names = FALSE)

    if (has_nul(out)) stop("regenerated file is itself corrupt: ", out)
    cat("  wrote", nrow(joined), "cells ->", out, "\n")
    repaired <- repaired + 1L
  }
}

cat("\nchecked ", checked, " samples; repaired ", repaired, "\n", sep = "")
if (repaired == 0L) cat("All percent.mt metadata files are readable.\n")
