# prepare the data to sce object and run EdgeR using contrast model

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

seurat.obj <- readRDS("../01-Regional_data_processing/combined/data/seurat_combined_v3.RDS")
seurat.obj <- JoinLayers(seurat.obj)

count_data <- seurat.obj@assays$RNA@layers$counts
colnames(count_data) <- colnames(seurat.obj)
rownames(count_data) <- rownames(seurat.obj)
sce.obj <- SingleCellExperiment(list(counts = count_data),colData = seurat.obj@meta.data)

sce.obj@colData$region_samples_ids <- paste0(sce.obj@colData$region,'.',sce.obj@colData$sample_ids)
sce.obj@colData$region_condition <- paste0(sce.obj@colData$region,'.',sce.obj@colData$condition_new)

sce.obj <- prepSCE(sce.obj, 
                   kid = "consensus_annotation_v3", # subpopulation assignments
                   gid = "region_condition",  # group IDs (ctrl/stim)
                   sid = "region_samples_ids",   # sample IDs (ctrl/stim.1234)
                   drop = F)  # drop all other colData columns


pb <- aggregateData(sce.obj,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

sce.obj <- computeLibraryFactors(sce.obj)
sce.obj <- logNormCounts(sce.obj)

pb.norm_mean <- aggregateData(sce.obj,
                              assay = "logcounts", fun = "mean",
                              by = c("cluster_id", "sample_id"))


ei <- metadata(sce.obj)$experiment_info

META_PATH = "/banach2/SCORCH/data/metadata/sample_metadata_v3.csv"
sample_meta <- read.csv(META_PATH)
sample_meta_filter = as.data.frame(sample_meta) %>%
  dplyr::filter(!is.na(cellranger_path)) %>%
  dplyr::filter((oud_jay_v2 == FALSE)) 
table(sample_meta_filter$hiv_jay_v2)

sample_meta_filter_selected <- sample_meta_filter %>%
  dplyr::select(sample_ids,region,age_value_v2,sex_v2,race_v2,postmortem_interval_v2,hiv_jay_v2)
sample_meta_filter_selected$sample_id <- paste0(sample_meta_filter_selected$region,'.',sample_meta_filter_selected$sample_ids)
ei <- ei %>%
  left_join(sample_meta_filter_selected,by = "sample_id")

ei$Age <- as.numeric(ei$age_value_v2)
ei$Age <- scale(ei$Age)[,1]
ei$Sex <- factor(ei$sex_v2,levels = c("male","female"))
ei$Race <- make.names(ei$race_v2)
### ~PMI (scaled)
ei$PMI <- as.numeric(ei$postmortem_interval_v2)
ei$PMI <- scale(ei$PMI)[,1]
ei$region_HIV <- ei$group_id 

mm <- model.matrix(~ 0 +ei$region_HIV+ei$Age+ei$Sex +ei$Race+ei$PMI )
dimnames(mm) <- list(ei$sample_id, c(levels(ei$region_HIV),
                                     "Age","Sex","RaceNAAI","RaceNAAIW",
                                     "RaceW","PMI"))

saveRDS(pb,"./cross_region/data/pseudo_bulk.rds")
saveRDS(pb.norm_mean,"./cross_region/data/pseudo_norm_mean_bulk.rds")
saveRDS(mm,"./cross_region/data/mm.rds")


contrast <- limma::makeContrasts("(INS.HIV-INS.CTR)/3+(PFC.HIV-PFC.CTR)/3+(VST.HIV-VST.CTR)/3", levels = mm)

res.edgeR <- pbDS(pb, design = mm,contrast=contrast, filter="both",min_cells=10,method = "edgeR")
res.edgeR_result <- resDS(sce.obj, res.edgeR, frq = TRUE)

saveRDS(res.edgeR_result,"./cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")

colnames(res.edgeR_result) <- sapply(colnames(res.edgeR_result),function(cname){
  num <- extract_numeric(cname)
  if(is.na(num)) return(cname)
  else return(gsub(num, paste0("a",num), cname))
})

res.edgeR_result <- res.edgeR_result %>%
  rowwise() %>%
  mutate(avg_CTR_INS_frq = mean(c_across(INS.HCCCY.frq:INS.HCTME.frq), na.rm = TRUE)) %>%
  mutate(avg_HIV_INS_frq = mean(c_across(INS.1063.frq:INS.1040.frq), na.rm = TRUE)) %>%
  mutate(avg_CTR_PFC_frq = mean(c_across(PFC.HCCPL.frq:PFC.HCTXJ.frq), na.rm = TRUE)) %>%
  mutate(avg_HIV_PFC_frq = mean(c_across(PFC.6801066772.frq:PFC.1133.frq), na.rm = TRUE)) %>%
  mutate(avg_CTR_VST_frq = mean(c_across(VST.HCTPK.frq:VST.HCTTS.frq), na.rm = TRUE)) %>%
  mutate(avg_HIV_VST_frq = mean(c_across(VST.10280.frq:VST.1059.frq), na.rm = TRUE)) %>%
  mutate(avg_CTR_frq = mean(c(avg_CTR_INS_frq,avg_CTR_PFC_frq,avg_CTR_VST_frq), na.rm = TRUE)) %>%
  
  mutate(avg_HIV_frq = mean(c(avg_HIV_INS_frq,avg_HIV_PFC_frq,avg_HIV_VST_frq), na.rm = TRUE)) 

res.edgeR_result <- res.edgeR_result %>%
  ungroup()
saveRDS(res.edgeR_result,"./cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq.RDS")

