# prepare the data to sce object and run DEG analysis

library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(scater)
library(tidyr)
library(dplyr)
library(cowplot)

wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/03-DE_analysis/"
setwd(wd)

seurat.obj <- readRDS("../01-Regional_data_processing/PFC/data/03-cell_integration_vis/integrated_seurat_obj_v3.RDS")

count_data <- seurat.obj@assays$RNA@layers$counts
colnames(count_data) <- colnames(seurat.obj)
rownames(count_data) <- rownames(seurat.obj)
sce.obj <- SingleCellExperiment(list(counts = count_data),colData = seurat.obj@meta.data)



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

META_PATH = "/banach2/SCORCH/data/metadata/sample_metadata_v3.csv"
sample_meta <- read.csv(META_PATH)
# only load samples of interest:
# HIV and CTR samples from PFC
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

saveRDS(pb.pfc,"./PFC/data/pseudo_bulk.rds")
saveRDS(pb.pfc.norm_mean,"./PFC/data/pseudo_norm_mean_bulk.rds")
saveRDS(mm,"./PFC/data/mm.rds")


contrast <- makeContrasts("HIV-CTR", levels = mm)

res.edgeR <- pbDS(pb.pfc, design = mm,contrast=contrast, filter="both",min_cells=10,method = "edgeR")
res.edgeR_result <- resDS(sce.obj, res.edgeR, frq = TRUE)

saveRDS(res.edgeR_result,"./PFC/results/PFC_HIV-CTR_res_edgeR_filtered_withFreq.RDS")

res.edgeR_result <- readRDS("./PFC/results/PFC_HIV-CTR_res_edgeR_filtered_withFreq.RDS")

colnames(res.edgeR_result) <- sapply(colnames(res.edgeR_result),function(cname){
  num <- extract_numeric(cname)
  if(is.na(num)) return(cname)
  else return(gsub(num, paste0("a",num), cname))
})
res.edgeR_result <- res.edgeR_result %>%
  rowwise() %>%
  mutate(avg_CTR_frq = mean(c_across(HCCPL.frq:HCTXJ.frq), na.rm = TRUE)) %>%
  mutate(avg_HIV_frq = mean(c_across(a6801066772.frq:a1133.frq), na.rm = TRUE)) %>%
  mutate(region = "PFC")

saveRDS(res.edgeR_result,"./PFC/results/PFC_HIV-CTR_res_edgeR_filtered_withFreq.RDS")

