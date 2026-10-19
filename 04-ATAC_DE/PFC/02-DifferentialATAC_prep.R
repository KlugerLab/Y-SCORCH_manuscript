# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Seurat)
library(Signac)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(scater)
library(tidyr)
library(dplyr)
library(cowplot)
## REVISION: limma is used but was not loaded by the published script.
library(limma)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT by rerun_setup() instead of a
##   hard-coded setwd().
rerun_setup("04-ATAC_DE",
            subdirs = c("PFC/data_ATAC"))
signac.obj.PFC <- readRDS("../02-ATAC_processing/sampleQC_redo/PFC/data/signac_integrated_v2.RDS")

DefaultAssay(signac.obj.PFC) <- "ATAC"
count_data <- signac.obj.PFC@assays$ATAC$counts
colnames(count_data) <- colnames(signac.obj.PFC)
rownames(count_data) <- rownames(signac.obj.PFC)
sce.obj <- SingleCellExperiment(list(counts = count_data),colData = signac.obj.PFC@meta.data)

sce.obj <- prepSCE(sce.obj, 
                   kid = "consensus_annotation_v3", # subpopulation assignments
                   gid = "condition_new",  # group IDs (ctrl/stim)
                   sid = "sample_ids",   # sample IDs (ctrl/stim.1234)
                   drop = F)  # drop all other colData columns


pb.pfc <- aggregateData(sce.obj,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
sce.obj <- computeLibraryFactors(sce.obj)
sce.obj <- logNormCounts(sce.obj)

pb.pfc.norm_mean <- aggregateData(sce.obj,
                                  assay = "logcounts", fun = "mean",
                                  by = c("cluster_id", "sample_id"))


ei <- metadata(sce.obj)$experiment_info
ei$group_id <- factor(ei$group_id,levels=c("HIV","CTR"))

META_PATH = SAMPLE_META   ## REVISION: rebuilt sample sheet
sample_meta <- read.csv(META_PATH)
sample_meta_filter = as.data.frame(sample_meta) %>%
  dplyr::filter(!is.na(cellranger_path)) %>%
  dplyr::filter((oud_jay_v2 == FALSE) & (region=="PFC") ) 
table(sample_meta_filter$hiv_jay_v2)

sample_meta_filter_selected <- sample_meta_filter %>%
  dplyr::select(sample_ids,age_value_v2,sex_v2,race_v2,postmortem_interval_v2,hiv_jay_v2)
sample_meta_filter_selected$sample_id <- sample_meta_filter_selected$sample_ids
ei <- ei %>%
  left_join(sample_meta_filter_selected,by = "sample_id")

ei$Age <- as.numeric(ei$age_value_v2)
ei$Age <- scale(ei$Age)[,1]
ei$Sex <- factor(ei$sex_v2,levels = c("male","female"))
ei$Race <- make.names(ei$race_v2)
### ~PMI (scaled)
ei$PMI <- as.numeric(ei$postmortem_interval_v2)
ei$PMI <- scale(ei$PMI)[,1]
ei$HIV <- ei$group_id 

mm <- model.matrix(~ 0 +ei$HIV+ei$Age+ei$Sex +ei$Race+ei$PMI )
dimnames(mm) <- list(ei$sample_id, c(levels(ei$HIV),
                                     "Age","Sex","RaceNAAI","RaceNAAIW",
                                     "RaceW","PMI"))
# save the data for other comparisons
saveRDS(pb.pfc,"./PFC/data_ATAC/pseudo_bulk.rds")
saveRDS(pb.pfc.norm_mean,"./PFC/data_ATAC/pseudo_norm_mean_bulk.rds")
saveRDS(mm,"./PFC/data_ATAC/mm.rds")



contrast <- makeContrasts("HIV-CTR", levels = mm)

res.edgeR <- pbDS(pb.pfc, design = mm,contrast=contrast, filter="both",min_cells=10,method = "edgeR")
res.edgeR_result <- resDS(sce.obj, res.edgeR, frq = TRUE)

colnames(res.edgeR_result) <- sapply(colnames(res.edgeR_result),function(cname){
  num <- extract_numeric(cname)
  if(is.na(num)) return(cname)
  else return(gsub(num, paste0("a",num), cname))
})

## REVISION: the published script averaged group frequencies over a positional
##   c_across() column range, which is wrong if sample order changes; groups are
##   now assigned by sample (add_group_frq() in revision/config/paths.R).
res.edgeR_result <- add_group_frq(res.edgeR_result, sce.obj, "PFC")

saveRDS(res.edgeR_result,"./PFC/data_ATAC/PFC_HIV-CTR_res_edgeR_filtered_withFreq.RDS")

