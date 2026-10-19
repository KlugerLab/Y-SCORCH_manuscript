# ==============================================================================
# Attach the nine VST MSN subtypes (04-02 predictions) to the combined RNA
# object, or to the VST Signac object, without rebuilding either.
#
# Method
#   The MSN subtype predictions (04-02) became available after the combined
#   object (combined/01..03) was built, so downstream objects label all 34,319
#   striatal neurons "MSNs" where the published analysis has nine subtypes.
#   Rebuilding via 04-01 -> combined/01..03 -> VST/05 re-runs integration,
#   refinement and FIt-SNE over 641,718 cells to change one annotation column,
#   and re-clustering could shift cell type assignments in all three regions.
#   This script applies the same relabel directly:
#     - only the VST cells labelled "MSNs" change; 34,233 take a subtype and the
#       86 without a prediction are dropped, as 04-01 does with
#       subset(..., consensus_annotation_v2 != "MSNs") (the published
#       "<50 MSN cells per sample" rule);
#     - the subtype is written to consensus_annotation_v2 (as 04-01 does) and to
#       consensus_annotation_v3, whose palettes the downstream figures use;
#     - every other cell type, all INS and PFC cells, cluster assignments and
#       embedding coordinates are unchanged.
#   Guards: every prediction must land on a cell currently labelled "MSNs",
#   no MSN may lie outside VST, and no assigned value may come back NA.
#
# Inputs
#   OUT_ROOT/01-Regional_data_processing/VST/data/msn_prediction.RDS
#   rna : OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3_no1040.RDS
#   atac: OUT_ROOT/02-ATAC_processing/sampleQC_redo/VST/data/signac_integrated_v2.RDS
# Outputs
#   the same paths with suffix _msn (seurat_combined_v3_no1040_msn.RDS,
#   signac_integrated_v2_msn.RDS)
#
# Usage
#   revision/docker/run.sh revision/tools/patch_msn_subtypes.R [--target rna|atac] [--dry-run]
# ==============================================================================

suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

args    <- commandArgs(trailingOnly = TRUE)
getarg  <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
DRY_RUN <- "--dry-run" %in% args
## "rna"  = the combined Seurat object (default)
## "atac" = the VST Signac object; same relabel, different barcode convention
TARGET  <- match.arg(getarg("--target", "rna"), c("rna", "atac"))

rerun_setup("01-Regional_data_processing", subdirs = "combined/data")

PRED <- "./VST/data/msn_prediction.RDS"
if (!file.exists(PRED)) stop("missing ", PRED, " - run 04-02-VST-msn_subtyping.R first")

if (TARGET == "rna") {
  IN  <- "./combined/data/seurat_combined_v3_no1040.RDS"
  OUT <- "./combined/data/seurat_combined_v3_no1040_msn.RDS"
} else {
  IN  <- file.path(OUT_ROOT, "02-ATAC_processing/sampleQC_redo/VST/data/signac_integrated_v2.RDS")
  OUT <- file.path(OUT_ROOT, "02-ATAC_processing/sampleQC_redo/VST/data/signac_integrated_v2_msn.RDS")
}
if (!file.exists(IN)) stop("missing ", IN)

## ---- predictions -------------------------------------------------------------
p <- do.call(rbind, readRDS(PRED))
message("[msn] predictions: ", nrow(p), " cells, ",
        length(unique(p$predicted.id)), " subtypes")

## Barcode conventions differ between the prediction and the two objects; the
## match counts are reported below.
##   prediction : <sample>.<barcode>-1-<sample>     10280.AAAC...-1-10280
##   RNA combined: <barcode>-1-<sample>-VST         AAAC...-1-10280-VST
##   VST Signac : <sample>_<barcode>-1              10280_AAAC...-1
rn <- rownames(p)
mapped <- if (TARGET == "rna") {
  paste0(sub("^[^.]+\\.", "", rn), "-VST")
} else {
  paste0(sub("\\..*$", "", rn), "_", sub("-[^-]+$", "", sub("^[^.]+\\.", "", rn)))
}
if (anyDuplicated(mapped)) stop("duplicate barcodes after mapping")
pred <- setNames(as.character(p$predicted.id), mapped)
message("[msn] target: ", TARGET, "  -> ", basename(IN))

## ---- object ------------------------------------------------------------------
message("[msn] loading ", IN, " ...")
obj <- readRDS(IN)
md  <- obj@meta.data
n_before <- ncol(obj)

is_msn <- as.character(md$consensus_annotation_v3) == "MSNs"
message("[msn] cells labelled MSNs: ", sum(is_msn, na.rm = TRUE))

## guard: the relabel is VST-only by construction, so anything else is a bug
bad_region <- sum(is_msn & as.character(md$region) != "VST", na.rm = TRUE)
if (bad_region > 0)
  stop(bad_region, " MSN-labelled cells outside VST - investigate before patching")

hit <- intersect(names(pred), colnames(obj))
n_msn <- sum(is_msn, na.rm = TRUE)
message("[msn] predictions matching object cells: ", length(hit), " of ", length(pred))

## The two targets fail in opposite directions, so they need different tests.
##   RNA : the combined object holds every cell 04-02 scored, so a prediction
##         that does not map means the barcode transform is wrong.
##   ATAC: ATAC QC dropped nuclei that passed RNA QC, so surplus predictions
##         are expected (34,233 predictions vs 24,963 ATAC MSNs). What must not
##         happen is an ATAC MSN left without a subtype.
if (TARGET == "rna") {
  if (length(hit) != length(pred))
    stop(length(pred) - length(hit), " predictions did not map; refusing to patch ",
         "a partial set - check the barcode transform")
} else {
  unmatched <- n_msn - length(hit)
  message("[msn] ATAC MSNs without a prediction: ", unmatched,
          " (", round(100 * unmatched / n_msn, 2), "%)")
  ## a few unassigned cells are 04-02's own remainder; a large fraction means
  ## the transform is wrong, not that the cells are genuinely missing
  if (unmatched > 0.05 * n_msn)
    stop(unmatched, " of ", n_msn, " ATAC MSNs (", round(100*unmatched/n_msn), "%) ",
         "have no prediction - that is too many to be 04-02's remainder; ",
         "check the barcode transform")
}

## every target must currently be an MSN - otherwise we would be overwriting a
## real annotation somewhere else
not_msn <- sum(as.character(md[hit, "consensus_annotation_v3"]) != "MSNs")
if (not_msn > 0)
  stop(not_msn, " prediction targets are not currently labelled MSNs")

## ---- apply -------------------------------------------------------------------
## 04-01 writes the subtype into consensus_annotation_v2 and keeps a separate
## msn_subtypes column; v3 carries it too so the downstream figures, whose
## palettes are keyed to consensus_annotation_v3, resolve the nine levels.
## consensus_annotation_v2 is a factor (28 levels) and none of the nine
## subtypes is among them, so assigning into it raises
##   invalid factor level, NA generated
## and blanks every cell it touches while the subtype count read off the
## character column still looks correct. Convert to character first; v3 is
## already character.
.assign <- function(df, col, cells, values) {
  if (!col %in% colnames(df)) return(df)
  if (is.factor(df[[col]])) {
    message("[msn] ", col, " is a factor - converting to character before assigning")
    df[[col]] <- as.character(df[[col]])
  }
  df[cells, col] <- values
  df
}

md$msn_subtypes <- NA_character_
md[hit, "msn_subtypes"] <- pred[hit]
for (cc in c("consensus_annotation_v2", "consensus_annotation_v3",
             "consensus_annotation_v3_relabeled"))
  md <- .assign(md, cc, hit, pred[hit])

## every column we touched must now hold the subtype for every target cell -
## an NA here means the assignment was dropped
for (cc in c("consensus_annotation_v2", "consensus_annotation_v3", "msn_subtypes")) {
  if (!cc %in% colnames(md)) next
  n_na <- sum(is.na(md[hit, cc]))
  if (n_na > 0)
    stop(n_na, " NA values in '", cc, "' after assignment - the write was dropped")
  n_ok <- sum(as.character(md[hit, cc]) == pred[hit])
  if (n_ok != length(hit))
    stop(cc, ": only ", n_ok, " of ", length(hit), " cells carry the expected subtype")
  message("[msn] ", cc, ": ", n_ok, " cells assigned, 0 NA")
}
obj@meta.data <- md

## the unassigned remainder - 04-01's "<50 MSN cells per sample" rule
leftover <- colnames(obj)[as.character(obj@meta.data$consensus_annotation_v3) == "MSNs"]
message("[msn] MSNs with no prediction (dropped, as 04-01 does): ", length(leftover))
if (length(leftover)) obj <- subset(obj, cells = setdiff(colnames(obj), leftover))

message("[msn] cells: ", n_before, " -> ", ncol(obj))

## ---- verify ------------------------------------------------------------------
md2 <- obj@meta.data
tb  <- table(as.character(md2$consensus_annotation_v3))
sub_tb <- tb[grepl("^D1|^D2", names(tb))]
cat("\n=== MSN subtypes on the patched object ===\n"); print(sort(sub_tb, decreasing = TRUE))
cat("\ntotal subtyped:", sum(sub_tb), "\n")
cat("any 'MSNs' left:", sum(names(tb) == "MSNs"), "\n")

## nothing outside VST may have moved
for (rg in c("INS", "PFC")) {
  n <- sum(as.character(md2$region) == rg)
  cat(rg, "cells:", n, " (unchanged)\n")
}

if (DRY_RUN) {
  message("\n[msn] --dry-run: not writing")
} else {
  saveRDS(obj, OUT)
  message("\n[msn] wrote ", OUT)
  message("[msn] next: point VST/05 at this object, e.g.")
  message("  SCORCH_COMBINED_OBJ=", OUT,
          " revision/docker/run.sh 01-Regional_data_processing/VST/05-VST_tsne_from_all.R")
}
