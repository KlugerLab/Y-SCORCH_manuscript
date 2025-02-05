# clean the cells and annotate the data

library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(scater)
library(tidyr)
library(dplyr)
library(cowplot)
wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/01-Regional_data_processing/VST/"
setwd(wd)

PROCESSED_DATA_PATH = "/banach2/SCORCH/analysis/240610_pipeline/vst/"
META_PATH = "/banach2/SCORCH/data/metadata/sample_metadata_v3.csv"
CELL_META_PATH = "/banach2/SCORCH/scripts/combined_analysis_paper/01-Regional_data_processing/VST/data/sample_metadata_withMT/"

sample_folders <- list.dirs(PROCESSED_DATA_PATH,recursive = F)
sample_folders_short <- list.dirs(PROCESSED_DATA_PATH,recursive = F,full.names = F)
# 56 samples in total

sample_meta <- read.csv(META_PATH)
# only load samples of interest:
# HIV and CTR samples from VST
sample_meta_filter = as.data.frame(sample_meta) %>%
  filter(!is.na(cellranger_path)) %>%
  filter((oud_jay_v2 == FALSE) & (region=="VST") ) 
table(sample_meta_filter$hiv_jay_v2)
# 13 CTRs, 20 HIVs


samples2keep <- sample_meta_filter$sample_name
# change the sample names from HIVOUD to HIV
samples2keep <- sapply(samples2keep,function(sname) gsub("OUD","",sname))
samples2keep %in% sample_folders_short
samples2load_paths <- paste0(PROCESSED_DATA_PATH,samples2keep)

# load the data into sce
sce.list <- lapply(samples2keep,function(sid) SingleCellExperiment(list(counts=readRDS(paste0(PROCESSED_DATA_PATH,sid,"/label_transfer/final_counts_",sid,".rds"))),
                                                                   colData = read.csv(paste0(CELL_META_PATH,sid,".csv"),row.names = 1)) )
sce.obj <- do.call(cbind,sce.list)
sce.obj <- computeLibraryFactors(sce.obj)
sce.obj <- logNormCounts(sce.obj)


sce.obj@colData[,"condition_new"] <- sapply(sce.obj@colData$sample_name,function(sid) sample_meta_filter[samples2keep == sid,"hiv_jay_v2"])
sce.obj@colData$condition_new[sce.obj@colData$condition_new == TRUE] <- "HIV"
sce.obj@colData$condition_new[sce.obj@colData$condition_new == FALSE] <- "CTR"



# 1. remove cells labeled with high MT
# double check the median
median_mt <- as.data.frame(sce.obj@colData) %>%
  select(sample_name,percent.mt,seurat_clusters) %>%
  group_by(sample_name,seurat_clusters) %>%
  summarise(median_mt = median(percent.mt)) %>%
  mutate(cluster2remove = median_mt > 10) %>%
  mutate(cluster_idx = as.character(seurat_clusters))
table(sce.obj@colData$sample_name[sce.obj@colData$cell2remove_highMT])
# only HCTKN has a cluster to remove

sce.obj <- sce.obj[,!sce.obj@colData$cell2remove_highMT]


# 2. clean up the cell types
unique(sce.obj@colData$biccn_predicted.id)
unique(sce.obj@colData$ma_predicted.id)
unique(sce.obj@colData$nhp_vst_predicted.id)
sce.obj@colData$biccn_predicted.id <- as.factor(sce.obj@colData$biccn_predicted.id)
sce.obj@colData$ma_predicted.id <- as.factor(sce.obj@colData$ma_predicted.id)
sce.obj@colData$nhp_vst_predicted.id <- as.factor(sce.obj@colData$nhp_vst_predicted.id)


metadata_all <- as.data.frame(sce.obj@colData)
metadata_all <- metadata_all %>%
  dplyr::select(starts_with("biccn_prediction.score")) %>%
  rowwise() %>%
  mutate(biccn_prediction.score_max = max(c_across(starts_with("biccn_prediction.score"))))
sce.obj@colData$biccn_prediction_score.max <-  metadata_all$biccn_prediction.score_max


metadata_all <- as.data.frame(sce.obj@colData)
metadata_all <- metadata_all %>%
  dplyr::select(starts_with("ma_prediction.score")) %>%
  rowwise() %>%
  mutate(ma_prediction.score_max = max(c_across(starts_with("ma_prediction.score"))))
sce.obj@colData$ma_prediction_score.max <-  metadata_all$ma_prediction.score_max

metadata_all <- as.data.frame(sce.obj@colData)
metadata_all <- metadata_all %>%
  dplyr::select(starts_with("nhp_vst_prediction.score")) %>%
  rowwise() %>%
  mutate(nhp_prediction.score_max = max(c_across(starts_with("nhp_vst_prediction.score"))))
sce.obj@colData$nhp_prediction_score.max <-  metadata_all$nhp_prediction.score_max

# remove cells with max scores < 0.5
metadata_all <- as.data.frame(sce.obj@colData) %>%
  mutate(cell2remove = ifelse((biccn_prediction_score.max < 0.5) & (ma_prediction_score.max < 0.5) & (nhp_prediction_score.max < 0.5),T,F))
sce.obj@colData$cell2remove_score <- metadata_all$cell2remove

# next: integrate the data to visualize
saveRDS(sce.obj,"./data/02-cell_annotation/processed_sce_obj.rds")




