# QC each sample using Signac
#
# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
## REVISION: dplyr is used but was not loaded by the published script.
library(dplyr)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT by rerun_setup() instead of a
##   hard-coded setwd().
rerun_setup("02-ATAC_processing/sampleQC_redo/PFC",
            subdirs = c("data"))
META_PATH = SAMPLE_META   ## REVISION: rebuilt sample sheet


sample_meta <- read.csv(META_PATH)
# only load samples of interest:
# HIV and CTR samples from PFC
sample_meta_filter = as.data.frame(sample_meta) %>%
  dplyr::filter(!is.na(cellranger_path)) %>%
  dplyr::filter((oud_jay_v2 == FALSE) & (region=="PFC") ) 

sample_meta_filter$cellranger_atac_path <- sapply(sample_meta_filter$cellranger_path,function(crpath) gsub("/filtered_feature_bc_matrix/", "/",crpath))

# load the atac data for each sample
## REVISION: read the CellRanger MTX directory instead of the .h5 file. HDF5
##   requests a file lock on every open, which hangs on the NFS data mount; the
##   MTX and .h5 matrices are identical, so `$Peaks` is unchanged. Per-sample
##   progress messages added.
.n_samp <- nrow(sample_meta_filter)
count_list <- lapply(seq_len(.n_samp), function(i) {
  crpath <- sample_meta_filter$cellranger_atac_path[i]
  sid    <- sample_meta_filter$sample_ids[i]
  t0 <- Sys.time()
  x  <- Read10X(paste0(crpath, "filtered_feature_bc_matrix/"))
  message(sprintf("[atac] %2d/%2d  %-16s peaks %d x %d  (%.0fs)",
                  i, .n_samp, sid, nrow(x$Peaks), ncol(x$Peaks),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  x
})
names(count_list) <- sample_meta_filter$sample_ids
message("[atac] loaded ", length(count_list), " samples")

fragpath_list <- lapply(sample_meta_filter$cellranger_atac_path,function(crpath) paste0(crpath,"atac_fragments.tsv.gz"))
names(fragpath_list) <- sample_meta_filter$sample_ids

#counts <- Read10X_h5(paste0(test_sample_path,"outs/filtered_feature_bc_matrix.h5"))
#fragpath <- paste0(test_sample_path,"outs/atac_fragments.tsv.gz")

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create ATAC assay and add it to the object
signac.obj.list <- lapply(sample_meta_filter$sample_ids, function(sid)
  CreateChromatinAssay(
    counts = count_list[[sid]]$Peaks,
    sep = c(":", "-"),
    fragments = fragpath_list[[sid]],
    annotation = annotation
  )
)
names(signac.obj.list) <- sample_meta_filter$sample_ids

signac.obj.list <- lapply(signac.obj.list,function(obj) CreateSeuratObject(
  counts = obj,
  assay = "peaks"
))

signac.obj.list <- lapply(signac.obj.list,function(obj){
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj)
  obj$nucleosome_group <- ifelse(obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  
  return(obj)
})


# TODO: generate the statistics of the nCount_peaks and TSS.enrichment

# TODO: add RNA labels and remove cells by the metrics

TSS_cutoff <- 1
NS_cutoff <- 2
#seurat.pfc <- readRDS("../../../01-Regional_data_processing/PFC/data/03-cell_integration_vis/integrated_seurat_obj.RDS")

seurat.all <- readRDS("../../../01-Regional_data_processing/combined/data/seurat_combined_v2.RDS")
seurat.pfc <- subset(seurat.all,cells = colnames(seurat.all)[seurat.all@meta.data$region == "PFC"])

RNA_metadata <- seurat.pfc@meta.data
# filter the data based on RNA and QC metrics
signac.obj.list_sub <- lapply(sample_meta_filter$sample_ids,function(sid){
  obj <- signac.obj.list[[sid]]
  obj <- subset(obj,cells = colnames(obj)[paste0(colnames(obj),"-",sid,"-","PFC") %in% rownames(RNA_metadata)])
  obj <- subset(
    x = obj,
    subset = nCount_peaks < 100000 &
      nCount_peaks > 1500 &
      nucleosome_signal < 2 &
      TSS.enrichment > 1
  )
  return(obj)
  
  
})
names(signac.obj.list_sub) <- sample_meta_filter$sample_ids

## REVISION: checkpoint the filtered per-sample list (as the INS script already
##   does). Loading the samples and computing NucleosomeSignal/TSSEnrichment take
##   several hours; the steps below are cheap.
message("[atac] samples after RNA intersection + QC: ", length(signac.obj.list_sub))
message("[atac] cells per sample: ",
        paste(sprintf("%s=%d", names(signac.obj.list_sub),
                      vapply(signac.obj.list_sub, ncol, integer(1))),
              collapse = "  "))
saveRDS(signac.obj.list_sub, "./data/signac_list_filtered.RDS")
message("[atac] checkpoint -> ./data/signac_list_filtered.RDS")

# add RNA metadata
for (sid in sample_meta_filter$sample_ids) {
  metadata_sid <- RNA_metadata[paste0(colnames(signac.obj.list_sub[[sid]]),"-",sid,"-","PFC"),]
  rownames(metadata_sid) <- colnames(signac.obj.list_sub[[sid]])
  signac.obj.list_sub[[sid]] <- AddMetaData(signac.obj.list_sub[[sid]],metadata_sid)
}

# peak calling for each cell type
# merge objects first
combined.peaks <- reduce(x = do.call(c,unlist(sapply(sample_meta_filter$sample_ids,function(sid)signac.obj.list_sub[[sid]]@assays$peaks@ranges),use.names=FALSE)))
# todo: get the same rownames of combined.peaks to signac obj
#combined.peaks_df <- as.data.frame(combined.peaks)
#combined.peaks_rownames <- paste0(combined.peaks_df$seqnames,'-',combined.peaks_df$start,'-',combined.peaks_df$end)


#signac.obj.list_sub <- lapply(signac.obj.list_sub,function(obj){
#  obj@assays$peaks <- obj@assays$peaks[combined.peaks,]
#  return(obj)
#})

signac.obj.list.merged <- lapply(sample_meta_filter$sample_ids, function(sid){
  frag_tmp <- CreateFragmentObject(
    path = fragpath_list[[sid]],
    cells = colnames(signac.obj.list_sub[[sid]])
  )
  count_tmp <- FeatureMatrix(
    fragments = frag_tmp,
    features = combined.peaks,
    cells = colnames(signac.obj.list_sub[[sid]])
  )
  tmp_assay <- CreateChromatinAssay(count_tmp, fragments = frag_tmp)
  tmp_obj <- CreateSeuratObject(tmp_assay, assay = "ATAC", meta.data=signac.obj.list_sub[[sid]]@meta.data)
  return(tmp_obj)
})
names(signac.obj.list.merged) <- names(signac.obj.list_sub)

## REVISION: the published script hard-coded the merge range (c(2:31)). If QC or
##   the RNA-barcode intersection drops a sample, merge() fails or silently drops
##   samples; the range is now derived from the list length.
.n_obj <- length(signac.obj.list.merged)
message("[atac] merging ", .n_obj, " sample objects: ",
        paste(names(signac.obj.list.merged), collapse = ", "))
stopifnot(.n_obj >= 2, !any(duplicated(names(signac.obj.list.merged))))

signac.obj.combined <- merge(
  x = signac.obj.list.merged[[1]],
  y = signac.obj.list.merged[2:.n_obj],
  add.cell.ids = names(signac.obj.list.merged)
)
message("[atac] merged: ", ncol(signac.obj.combined), " cells x ",
        nrow(signac.obj.combined), " peaks")


saveRDS(signac.obj.combined,"./data/signac_combined.RDS")


