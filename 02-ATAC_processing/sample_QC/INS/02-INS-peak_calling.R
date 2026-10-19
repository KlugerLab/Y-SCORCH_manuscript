# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Signac)
library(Seurat)
#install.packages("https://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.21.tar.gz", repos=NULL, type="source")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

#install.packages("harmony")
library(harmony)


source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT by rerun_setup() instead of a
##   hard-coded setwd().
rerun_setup("02-ATAC_processing/sampleQC_redo/INS",
            subdirs = c("data"))
## REVISION: resume from the FeatureMatrix checkpoint if it exists (MACS2 +
##   FeatureMatrix take hours; the steps after it take minutes). Set
##   ATAC_FORCE_RESTART=1 to ignore the checkpoint and recount from scratch.
.ckpt_fm <- "./data/signac_postFeatureMatrix.rds"
.resume  <- file.exists(.ckpt_fm) && !nzchar(Sys.getenv("ATAC_FORCE_RESTART"))

if (.resume) {
  message("[atac] RESUMING from ", .ckpt_fm, " - skipping CallPeaks + FeatureMatrix")
  signac.obj.new <- readRDS(.ckpt_fm)
  message("[atac] resumed object: ", ncol(signac.obj.new), " cells x ",
          nrow(signac.obj.new), " peaks")
} else {

signac.obj.combined <- readRDS("./data/signac_combined.RDS")
#test <- subset(signac.obj.combined,cells=colnames(signac.obj.combined)[signac.obj.combined@meta.data$sample_ids == "HCCCY"])
#test@assays$ATAC@fragments[c(2:31)] <- NULL
peaks <- CallPeaks(
  object = signac.obj.combined,
  ## REVISION: the published hard-coded MACS2 path is not readable when the
  ##   container runs as the invoking user (revision/docker/run.sh passes --user);
  ##   MACS2_PATH is resolved in revision/config/paths.R.
  macs2.path = MACS2_PATH,
  group.by = "consensus_annotation_v3",
  outdir = "./data/"
)

saveRDS(peaks,"./data/peaks.RDS")

new_count_data <- FeatureMatrix(
  fragments = signac.obj.combined@assays$ATAC@fragments,
  features = peaks,
  cells = colnames(signac.obj.combined)
)

new_assay <- CreateChromatinAssay(new_count_data, fragments = signac.obj.combined@assays$ATAC@fragments)
signac.obj.new <- CreateSeuratObject(new_assay, assay = "ATAC", meta.data=signac.obj.combined@meta.data)

# integrate the data
#signac.obj.new[["ATAC"]] <- split(signac.obj.new[["ATAC"]], f = signac.obj.new$sample_ids)

## REVISION: checkpoint after FeatureMatrix, the most expensive step in this
##   script, so a failure in TFIDF, SVD, Harmony or FIt-SNE does not discard it.
saveRDS(signac.obj.new, "./data/signac_postFeatureMatrix.rds")
message("[atac] checkpoint -> ./data/signac_postFeatureMatrix.rds  (",
        ncol(signac.obj.new), " cells x ", nrow(signac.obj.new), " peaks)")

}  ## REVISION: end of the CallPeaks + FeatureMatrix block skipped on resume

signac.obj.new <- RunTFIDF(signac.obj.new)
signac.obj.new <- FindTopFeatures(signac.obj.new)
signac.obj.new <- RunSVD(signac.obj.new)

# this is not optimal: first component not ruled out
signac.obj.new <- RunHarmony(
  object = signac.obj.new,
  group.by.vars = 'sample_name',
  reduction.use = 'lsi',
  ## REVISION: assay.use removed; harmony 2.0.5 rejects it. Harmony runs on the
  ##   'lsi' reduction either way, so the result is unchanged.
  project.dim = FALSE
)


## REVISION: checkpoint after Harmony, before FIt-SNE (an external binary and
##   the most failure-prone step), so a t-SNE failure does not discard the
##   earlier steps. FAST_TSNE_PATH comes from revision/config/paths.R.
saveRDS(signac.obj.new, "./data/signac_postHarmony.rds")
message("[atac] checkpoint -> ./data/signac_postHarmony.rds")

signac.obj.new <-  RunTSNE(
  signac.obj.new,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path=FAST_TSNE_PATH, ann_not_vptree=FALSE, nthreads=12
)


DimPlot(signac.obj.new,reduction = "tsne",group.by = "consensus_annotation_v3")
DimPlot(signac.obj.new,reduction = "tsne",group.by = "sample_ids")

saveRDS(signac.obj.new,"./data/signac_integrated.RDS")