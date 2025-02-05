library(Signac)
library(Seurat)
#install.packages("https://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.21.tar.gz", repos=NULL, type="source")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

#install.packages("harmony")
library(harmony)


setwd("/banach2/SCORCH/scripts/combined_analysis_paper/02-ATAC_processing/sampleQC_redo/INS/")
signac.obj.combined <- readRDS("./data/signac_combined.RDS")
#test <- subset(signac.obj.combined,cells=colnames(signac.obj.combined)[signac.obj.combined@meta.data$sample_ids == "HCCCY"])
#test@assays$ATAC@fragments[c(2:31)] <- NULL
peaks <- CallPeaks(
  object = signac.obj.combined,
  macs2.path = "/home/rstudio/.local/bin/macs2",
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
signac.obj.new <- RunTFIDF(signac.obj.new)
signac.obj.new <- FindTopFeatures(signac.obj.new)
signac.obj.new <- RunSVD(signac.obj.new)

# this is not optimal: first component not ruled out
signac.obj.new <- RunHarmony(
  object = signac.obj.new,
  group.by.vars = 'sample_name',
  reduction.use = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

signac.obj.new <-  RunTSNE(
  signac.obj.new,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path="/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)


DimPlot(signac.obj.new,reduction = "tsne",group.by = "consensus_annotation_v3")
DimPlot(signac.obj.new,reduction = "tsne",group.by = "sample_ids")

saveRDS(signac.obj.new,"./data/signac_integrated.RDS")