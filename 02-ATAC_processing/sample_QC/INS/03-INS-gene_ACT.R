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
signac.obj.INS <- readRDS("../sampleQC_redo/INS/data/signac_integrated.RDS")

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(signac.obj.INS) <- annotations
gene.activities <- GeneActivity(signac.obj.INS)

# add the gene activity matrix to the Seurat object as a new assay
signac.obj.INS[['GeneACT']] <- CreateAssayObject(counts = gene.activities)
signac.obj.INS <- NormalizeData(
  object = signac.obj.INS,
  assay = 'GeneACT',
  normalization.method = 'LogNormalize',
  scale.factor = median(signac.obj.INS$nCount_GeneACT)
)

DefaultAssay(signac.obj.INS) <- 'GeneACT'

saveRDS(signac.obj.INS,"../sampleQC_redo/INS/data/signac_integrated_v2.RDS")