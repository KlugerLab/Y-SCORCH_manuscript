# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Seurat)
library(SingleCellExperiment)
library(tidyr)
library(dplyr)
library(cowplot)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: the published script never set a working directory (it relied on the
##   interactive session already being in the VST stage directory).
rerun_setup(file.path("01-Regional_data_processing", "VST"),
            subdirs = c("data"))

options(future.globals.maxSize = 1e12)

seu.obj.msn <- readRDS("./data/msn_subset.RDS")


## REVISION: reference path comes from MSN_REF_RDS (revision/config/paths.R); fail
## with a clear message if it is missing.
if (!file.exists(MSN_REF_RDS)) {
  stop("Missing primate MSN reference: ", MSN_REF_RDS, "\n",
       "  Needed for the 9 MSN subtypes reported in the paper. ",
       "Set MSN_REF_RDS or place the file at the path above.")
}
msn_ref <- readRDS(MSN_REF_RDS)

# for each sample, run against the reference
Idents(seu.obj.msn) <- "sample_ids"
seurat_query_list <- SplitObject(seu.obj.msn, split.by = "sample_ids")

## REVISION: per-sample tryCatch.
##   TransferData(k.weight = 30) errors when a sample yields fewer than 30
##   anchors, and the ncol <= 50 guard does not prevent that - a 60-cell sample
##   can still fall short. Inside a bare lapply, one such error aborts the whole
##   run and discards every sample already transferred (~15 min each). Failures
##   are recorded and reported instead, so a single thin sample costs that
##   sample rather than the job.
.msn_failures <- character(0)

prediction.df <- lapply(names(seurat_query_list), function(sid){
  seurat_obj <- seurat_query_list[[sid]]
  message("[msn] ", sid, ": ", ncol(seurat_obj), " MSN nuclei")
  if(ncol(seurat_obj) <= 50){
    message("[msn] ", sid, ": skipped (<= 50 nuclei)")
    return(NULL)
  }
  tryCatch({
    seurat_obj <- SCTransform(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)

    transfer_anchors_vst_ref <- FindTransferAnchors(reference = msn_ref, query = seurat_obj,
                                                    normalization.method = "SCT",
                                                    reduction = "pcaproject")
    predictions_msn <- TransferData(anchorset = transfer_anchors_vst_ref, k.weight = 30,
                                    refdata = msn_ref$MSN_type)
    message("[msn] ", sid, ": transferred ", nrow(predictions_msn), " predictions")
    predictions_msn
  }, error = function(e){
    message("[msn] ", sid, ": FAILED - ", conditionMessage(e))
    .msn_failures <<- c(.msn_failures, sid)
    NULL
  })
})
names(prediction.df) <- names(seurat_query_list)

saveRDS(prediction.df,"./data/msn_prediction.RDS")

## Report explicitly which samples produced nothing, and why. A NULL entry is
## otherwise indistinguishable from a sample with no MSNs.
.n_ok <- sum(!vapply(prediction.df, is.null, logical(1)))
message("[msn] samples with predictions: ", .n_ok, " of ", length(prediction.df))
if (length(.msn_failures))
  message("[msn] FAILED samples (kept as NULL): ", paste(.msn_failures, collapse = ", "))
message("[msn] wrote ./data/msn_prediction.RDS")
