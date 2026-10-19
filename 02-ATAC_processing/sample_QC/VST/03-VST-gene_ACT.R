# add gene act to the object
#
# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Seurat)
library(Signac)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(EnsDb.Hsapiens.v86)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT by rerun_setup() instead of a
##   hard-coded setwd().
rerun_setup("02-ATAC_processing/plotting")
signac.obj.VST <- readRDS("../sampleQC_redo/VST/data/signac_integrated.RDS")

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(signac.obj.VST) <- annotations

## REVISION: checkpoint the GeneActivity matrix, which takes hours to compute
##   from the fragment files. Resumes from the file if it exists; set
##   ATAC_FORCE_RESTART=1 to recompute.
.ckpt_ga <- "../sampleQC_redo/VST/data/gene_activities.rds"
if (file.exists(.ckpt_ga) && !nzchar(Sys.getenv("ATAC_FORCE_RESTART"))) {
  message("[atac] RESUMING gene activities from ", .ckpt_ga)
  gene.activities <- readRDS(.ckpt_ga)
} else {
  gene.activities <- GeneActivity(signac.obj.VST)
  saveRDS(gene.activities, .ckpt_ga)
  message("[atac] checkpoint -> gene_activities.rds  (",
          nrow(gene.activities), " genes x ", ncol(gene.activities), " cells)")
}

# add the gene activity matrix to the Seurat object as a new assay
signac.obj.VST[['GeneACT']] <- CreateAssayObject(counts = gene.activities)
signac.obj.VST <- NormalizeData(
  object = signac.obj.VST,
  assay = 'GeneACT',
  normalization.method = 'LogNormalize',
  scale.factor = median(signac.obj.VST$nCount_GeneACT)
)

DefaultAssay(signac.obj.VST) <- 'GeneACT'

saveRDS(signac.obj.VST,"../sampleQC_redo/VST/data/signac_integrated_v2.RDS")