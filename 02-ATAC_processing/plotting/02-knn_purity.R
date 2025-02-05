# compute how often the nearest neighbors are of the same cell type
# in the atac data

library(tidyr)
library(dplyr)
#BiocManager::install("Repitools")
library(Repitools)
library(ggpubr)
library(cowplot)

setwd("/banach2/SCORCH/scripts/combined_analysis_paper/02-ATAC_processing/plotting//")

signac.obj.INS <- readRDS("../sampleQC_redo/INS/data/signac_integrated_v2.RDS")
signac.obj.PFC <- readRDS("../sampleQC_redo/PFC/data/signac_integrated_v2.RDS")
signac.obj.VST <- readRDS("../sampleQC_redo/VST/data/signac_integrated_v2.RDS")


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
