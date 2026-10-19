# compute how often the nearest neighbors are of the same cell type
# in the atac data
#
# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(tidyr)
library(dplyr)
#BiocManager::install("Repitools")
## REVISION: Repitools import disabled. It cannot be installed (its dependency
##   gsmoothr was archived from CRAN) and no Repitools function is called here.
# library(Repitools)
library(ggpubr)
library(cowplot)
## REVISION: Seurat is used but was not loaded by the published script.
library(Seurat)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT by rerun_setup() instead of a
##   hard-coded setwd().
rerun_setup("02-ATAC_processing/plotting/",
            subdirs = c("data"))
signac.obj.INS <- readRDS("../sampleQC_redo/INS/data/signac_integrated_v2.RDS")
signac.obj.PFC <- readRDS("../sampleQC_redo/PFC/data/signac_integrated_v2.RDS")
## REVISION: VST reads the object carrying the nine MSN subtypes
##   (revision/tools/patch_msn_subtypes.R --target atac). The unpatched object labels every
##   striatal neuron "MSNs", which the palette and consensus_annotation_v3 grouping
##   below do not expect.
signac.obj.VST <- readRDS("../sampleQC_redo/VST/data/signac_integrated_v2_msn.RDS")


signac.obj.INS <- FindNeighbors(object = signac.obj.INS,
                                reduction = 'lsi', dims = 2:50,k.param = 100)
signac.obj.PFC <- FindNeighbors(object = signac.obj.PFC,
                                reduction = 'lsi', dims = 2:50,k.param = 100)
signac.obj.VST <- FindNeighbors(object = signac.obj.VST,
                                reduction = 'lsi', dims = 2:50,k.param = 100)


# for each cell, check its label purity in the neighborhood
cell_types.ins <- signac.obj.INS@meta.data$consensus_annotation_v3
cell_types.pfc <- signac.obj.PFC@meta.data$consensus_annotation_v3
cell_types.vst <- signac.obj.VST@meta.data$consensus_annotation_v3



purity.ins <- unlist(lapply(1:ncol(signac.obj.INS),
                            function(i)mean(cell_types.ins[i] == cell_types.ins[as.logical(signac.obj.INS@graphs$ATAC_nn[i,])])))
purity.pfc <- unlist(lapply(1:ncol(signac.obj.PFC),
                            function(i)mean(cell_types.pfc[i] == cell_types.pfc[as.logical(signac.obj.PFC@graphs$ATAC_nn[i,])])))
purity.vst <- unlist(lapply(1:ncol(signac.obj.VST),
                            function(i)mean(cell_types.vst[i] == cell_types.vst[as.logical(signac.obj.VST@graphs$ATAC_nn[i,])])))

# TODO: only include cell types with more than 100 cells
c2remove.ins <- names(table(cell_types.ins))[table(cell_types.ins) <= 100]
purity.ins.sub <- purity.ins[!(cell_types.ins %in% c2remove.ins)]

c2remove.pfc <- names(table(cell_types.pfc))[table(cell_types.pfc) <= 100]
purity.pfc.sub <- purity.pfc[!(cell_types.pfc %in% c2remove.pfc)]

c2remove.vst <- names(table(cell_types.vst))[table(cell_types.vst) <= 100]
purity.vst.sub <- purity.vst[!(cell_types.vst %in% c2remove.vst)]

purity.df <- data.frame("purity"=c(purity.ins.sub,purity.pfc.sub,purity.vst.sub),
                               "cell_type"=c(cell_types.ins[!(cell_types.ins %in% c2remove.ins)],
                                             cell_types.pfc[!(cell_types.pfc %in% c2remove.pfc)],
                                             cell_types.vst[!(cell_types.vst %in% c2remove.vst)]),
                            "region"=c(rep("INS",length(purity.ins.sub)),
                                       rep("PFC",length(purity.pfc.sub)),
                                       rep("VST",length(purity.vst.sub))))
purity.avg.df <- purity.df %>%
  group_by(region) %>%
  summarise(mean_purity = mean(purity))
saveRDS(purity.df,"./data/purity.df.rds")
