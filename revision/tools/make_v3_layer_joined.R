# ==============================================================================
# Produces seurat_combined_v3_layer_joined.RDS, the input of 05-communication/NICHES.R.
#
# No published script writes this object; the published run has it at the same
# size as seurat_combined_v3.RDS, i.e. the same object with the RNA layers
# joined. This reconstructs that step (load v3, JoinLayers(), save); it is not a
# new analysis.
#
# Why it is needed: combined/03 leaves the RNA assay split (one counts layer per
# sample) because it re-integrates. NICHES needs a single matrix: RunNICHES()
# pulls assay = "RNA" per sample, and LayerData() on a split assay returns only
# the first layer with a warning rather than an error, which would silently run
# the communication analysis on a fraction of the cells.
#
# Input:  OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3.RDS
# Output: OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3_layer_joined.RDS
#
# Usage: revision/docker/run.sh revision/tools/make_v3_layer_joined.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

rerun_setup("01-Regional_data_processing", subdirs = "combined/data")

IN  <- "./combined/data/seurat_combined_v3.RDS"
OUT <- "./combined/data/seurat_combined_v3_layer_joined.RDS"
stopifnot(file.exists(IN))

message("[joined] loading ", IN, " ...")
obj <- readRDS(IN)
n_before <- ncol(obj)
lyr_before <- length(Layers(obj[["RNA"]]))
message("[joined] ", n_before, " cells, ", lyr_before, " RNA layers")

obj <- JoinLayers(obj)
lyr_after <- length(Layers(obj[["RNA"]]))
message("[joined] after JoinLayers: ", ncol(obj), " cells, ", lyr_after, " layers (",
        paste(Layers(obj[["RNA"]]), collapse = ", "), ")")

## Nothing may be lost or gained by joining.
stopifnot(ncol(obj) == n_before)
if (lyr_after >= lyr_before && lyr_before > 2)
  warning("layer count did not fall (", lyr_before, " -> ", lyr_after,
          ") - check the object was actually split")

## The columns NICHES.R depends on.
need <- c("consensus_annotation_v3", "sample_name", "condition_new")
miss <- setdiff(need, colnames(obj@meta.data))
if (length(miss)) stop("missing metadata for NICHES: ", paste(miss, collapse = ", "))

## NICHES subsets to these six populations and needs >= 10 cells per identity
## per sample, so report the per-sample picture rather than discovering it
## halfway through the run.
CPS <- c("Micro", "Astro", "Endo", "L2/3 IT", "L4 IT", "Pvalb")
have <- intersect(CPS, unique(as.character(obj@meta.data$consensus_annotation_v3)))
message("[joined] NICHES cell types present: ", paste(have, collapse = ", "))
if (length(setdiff(CPS, have)))
  warning("NICHES cell type(s) absent from v3: ", paste(setdiff(CPS, have), collapse = ", "))

tab <- table(as.character(obj@meta.data$sample_name),
             as.character(obj@meta.data$consensus_annotation_v3))
tab <- tab[, intersect(colnames(tab), CPS), drop = FALSE]
ok <- rowSums(tab >= 10)
message("[joined] samples with >=10 cells in all ", ncol(tab), " NICHES types: ",
        sum(ok == ncol(tab)), " of ", nrow(tab))
if (any(ok < ncol(tab)))
  message("[joined] samples short in >=1 type (NICHES skips those identities): ",
          paste(utils::head(rownames(tab)[ok < ncol(tab)], 8), collapse = ", "))

message("[joined] writing ", OUT, " ...")
saveRDS(obj, OUT)
message("[joined] DONE")
