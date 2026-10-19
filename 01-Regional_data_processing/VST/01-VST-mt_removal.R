# Remove clusters with high mt
# cluster each sample, remove clusters with median mt > 10%
#
# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION;
# the high-MT cluster-removal logic is unchanged.

library(Seurat)
library(dplyr)
library(ggpubr)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

REGION <- "VST"

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup(file.path("01-Regional_data_processing", REGION),
            subdirs = c("data/sample_metadata_withMT", "figures/sample_MT_check"))

parent_dir <- paste0(region_dir(REGION), "/")           ## REVISION: CellBender-corrected pipeline
samples <- list.dirs(parent_dir, recursive = F, full.names = F)

## REVISION: the corrected pipeline holds only the paper sample set (the published
## pipeline folder also held OUD samples, dropped in 02-cell_annotation). Assert the
## set matches the sample sheet so a missing sample folder fails instead of
## shrinking the cohort.
expected <- load_sample_meta(REGION)$sample_name
stopifnot(setequal(samples, expected))
message("Processing ", length(samples), " ", REGION, " samples")

for (sample in samples) {
  seurat.obj <- readRDS(paste0(parent_dir,sample,"/label_transfer/labeled_",sample,".rds"))

  ## REVISION: the published code read labeled_metadata_<sample>.csv, which in the
  ## published VST pipeline carried percent.mt. The corrected VST metadata has no
  ## percent.mt, so VST uses the regenerated metadata in FIXED_META_DIR (percent.mt
  ## recomputed from raw_<sample>.rds), as INS and PFC do.
  metadata <- read.csv(fixed_meta_csv(REGION, sample), row.names = 1)

  ## REVISION: align metadata to the object by barcode. The published code assigned
  ## it wholesale, relying on row order, where a mismatch would silently corrupt percent.mt.
  missing <- setdiff(colnames(seurat.obj), rownames(metadata))
  if (length(missing) > 0) {
    stop(sample, ": ", length(missing), " cells have no percent.mt row ",
         "(barcode mismatch between labeled object and fixed metadata)")
  }
  metadata <- metadata[colnames(seurat.obj), , drop = FALSE]
  stopifnot(identical(rownames(metadata), colnames(seurat.obj)),
            !any(is.na(metadata$percent.mt)))
  seurat.obj@meta.data <- metadata

  # note: PCA computed from SCT
  DefaultAssay(seurat.obj) <- "SCT"
  seurat.obj <- FindNeighbors(seurat.obj)
  seurat.obj <- FindClusters(seurat.obj,resolution = 3)

  p1 <- FeaturePlot(seurat.obj,features = "percent.mt")
  p2 <- DimPlot(seurat.obj,reduction = "umap",group.by = "seurat_clusters",label = T)+ NoLegend()
  p3 <- VlnPlot(seurat.obj,features = "percent.mt",group.by = "seurat_clusters")+ NoLegend()

  median_mt <- seurat.obj@meta.data %>%
    select(percent.mt,seurat_clusters) %>%
    group_by(seurat_clusters) %>%
    summarise(median_mt = median(percent.mt)) %>%
    mutate(cluster2remove = median_mt > 10) %>%
    mutate(cluster_idx = as.character(seurat_clusters))

  seurat.obj@meta.data$cell2remove_highMT <- seurat.obj@meta.data$seurat_clusters %in%
    median_mt[median_mt$cluster2remove,]$cluster_idx

  message(paste0("Sample: ",sample,",",sum(seurat.obj@meta.data$cell2remove_highMT)," cells removed in ",sum(median_mt$cluster2remove)," clusters"))

  # output
  write.csv(seurat.obj@meta.data,paste0("./data/sample_metadata_withMT/",sample,".csv"))

  png(filename = paste0("./figures/sample_MT_check/",sample,".png"),width = 1500,height = 800,res = 100)
  print(ggarrange(p1,p2,p3,nrow = 1,ncol=3))
  dev.off()
}
