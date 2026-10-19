# QC each sample using Signac
#
# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Signac)
library(Seurat)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT by rerun_setup() instead of a
##   hard-coded setwd().
rerun_setup("02-ATAC_processing/sampleQC_redo/VST",
            subdirs = c("data", "figures/qc_metrics"))
META_PATH = SAMPLE_META   ## REVISION: rebuilt sample sheet


sample_meta <- read.csv(META_PATH)
# only load samples of interest:
# HIV and CTR samples from PFC
sample_meta_filter = as.data.frame(sample_meta) %>%
  dplyr::filter(!is.na(cellranger_path)) %>%
  dplyr::filter((oud_jay_v2 == FALSE) & (region=="VST") ) 

sample_meta_filter$cellranger_atac_path <- sapply(sample_meta_filter$cellranger_path,function(crpath) gsub("/filtered_feature_bc_matrix/", "/",crpath))

# load the atac data for each sample
## REVISION: read the CellRanger MTX directory instead of the .h5 file.
##   CellRanger writes the same filtered matrix in both formats (checked on one
##   sample: identical Peaks 154,691 x 14,569), so the `$Peaks` indexing below is
##   unchanged. HDF5 >= 1.10 requests a POSIX file lock on every open, including
##   read-only, and on the NFS data mount (no lock service, `hard` mount) that
##   request hangs indefinitely; gzipped MTX takes no lock. Setting
##   HDF5_USE_FILE_LOCKING=FALSE (done in revision/docker/run.sh) only helps jobs
##   launched through that script. MTX parsing costs ~47 s per sample. Per-sample
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

#DensityScatter(test.obj, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)


#FragmentHistogram(object = test.obj, group.by = 'nucleosome_group')

for (sid in sample_meta_filter$sample_ids) {
  v1 <-  VlnPlot(
    object = signac.obj.list[[sid]],
    features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal"),
    ncol = 3,
    pt.size = 0
  )
  png(file=paste0("./figures/qc_metrics/",sid,".png"),width = 2500,height = 1500,res=300)
  print(v1)
  dev.off()
  
}


# TODO: generate the statistics of the nCount_peaks and TSS.enrichment

# TODO: add RNA labels and remove cells by the metrics

TSS_cutoff <- 1
NS_cutoff <- 2
seurat.all <- readRDS("../../../01-Regional_data_processing/combined/data/seurat_combined_v2.RDS")
seurat.vst <- subset(seurat.all,cells = colnames(seurat.all)[seurat.all@meta.data$region == "VST"])


#seurat.pfc <- readRDS("../../../01-Regional_data_processing/VST/data/04-VST-tsne_marker_plotting/integrated_seurat_obj.RDS")
RNA_metadata <- seurat.vst@meta.data
# filter the data based on RNA and QC metrics

RNA_sids <- unique(RNA_metadata$sample_ids)
# only keep the samples passed RNA filter
signac.obj.list_sub <- lapply(RNA_sids,function(sid){
  obj <- signac.obj.list[[sid]]
  #print(sid)
  obj <- subset(obj,cells = colnames(obj)[paste0(colnames(obj),"-",sid,"-","VST") %in% rownames(RNA_metadata)])
  
  cond = (obj@meta.data$nCount_peaks < 100000)&
    (obj@meta.data$nCount_peaks > 1500)&
    (obj@meta.data$nucleosome_signal < 2)&
    (obj@meta.data$TSS.enrichment > 1)
  
  if(sum(cond) == 0){
    return(NULL)
  }
  else{
    obj <- subset(
      x = obj,
      subset = nCount_peaks < 100000 &
        nCount_peaks > 1500 &
        nucleosome_signal < 2 &
        TSS.enrichment > 1
    )
    return(obj)
  }
  
  
  
})
#24 out of 26 samples passed filters
names(signac.obj.list_sub) <- RNA_sids
signac.obj.list_sub <- signac.obj.list_sub[lengths(signac.obj.list_sub)!=0]

## REVISION: report the surviving sample count and checkpoint the filtered
##   per-sample list. Loading the samples and computing NucleosomeSignal/
##   TSSEnrichment take ~5 h; the steps below are comparatively cheap.
.dropped <- setdiff(RNA_sids, names(signac.obj.list_sub))
message("[atac] samples passing RNA intersection + QC: ",
        length(signac.obj.list_sub), " of ", length(RNA_sids))
if (length(.dropped))
  message("[atac] dropped (no cells passed): ", paste(.dropped, collapse = ", "))
message("[atac] cells per surviving sample: ",
        paste(sprintf("%s=%d", names(signac.obj.list_sub),
                      vapply(signac.obj.list_sub, ncol, integer(1))),
              collapse = "  "))

.ckpt <- "./data/signac_obj_list_sub.rds"
saveRDS(signac.obj.list_sub, .ckpt)
message("[atac] checkpoint -> ", .ckpt,
        "  (rerun from here if a later step fails)")

# add RNA metadata
for (sid in names(signac.obj.list_sub)) {
  metadata_sid <- RNA_metadata[paste0(colnames(signac.obj.list_sub[[sid]]),"-",sid,"-","VST"),]
  rownames(metadata_sid) <- colnames(signac.obj.list_sub[[sid]])
  signac.obj.list_sub[[sid]] <- AddMetaData(signac.obj.list_sub[[sid]],metadata_sid)
}


# peak calling for each cell type
# merge objects first
combined.peaks <- GenomicRanges::reduce(x = do.call(c,unlist(sapply(names(signac.obj.list_sub),function(sid)signac.obj.list_sub[[sid]]@assays$peaks@ranges),use.names=FALSE)))
# todo: get the same rownames of combined.peaks to signac obj
#combined.peaks_df <- as.data.frame(combined.peaks)
#combined.peaks_rownames <- paste0(combined.peaks_df$seqnames,'-',combined.peaks_df$start,'-',combined.peaks_df$end)


#signac.obj.list_sub <- lapply(signac.obj.list_sub,function(obj){
#  obj@assays$peaks <- obj@assays$peaks[combined.peaks,]
#  return(obj)
#})

signac.obj.list.merged <- lapply(names(signac.obj.list_sub), function(sid){
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

## REVISION: the published script hard-coded the merge range (c(2:24)). If QC or
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


