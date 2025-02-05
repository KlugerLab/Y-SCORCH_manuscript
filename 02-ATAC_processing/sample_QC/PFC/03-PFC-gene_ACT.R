# add gene act to the object

library(Seurat)
library(Signac)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(EnsDb.Hsapiens.v86)

setwd("/banach2/SCORCH/scripts/combined_analysis_paper/02-ATAC_processing/plotting/")
signac.obj.PFC <- readRDS("../sampleQC_redo/PFC/data/signac_integrated.RDS")

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(signac.obj.PFC) <- annotations
gene.activities <- GeneActivity(signac.obj.PFC)

# add the gene activity matrix to the Seurat object as a new assay
signac.obj.PFC[['GeneACT']] <- CreateAssayObject(counts = gene.activities)
signac.obj.PFC <- NormalizeData(
  object = signac.obj.PFC,
  assay = 'GeneACT',
  normalization.method = 'LogNormalize',
  scale.factor = median(signac.obj.PFC$nCount_GeneACT)
)

DefaultAssay(signac.obj.PFC) <- 'GeneACT'

saveRDS(signac.obj.PFC,"../sampleQC_redo/PFC/data/signac_integrated_v2.RDS")