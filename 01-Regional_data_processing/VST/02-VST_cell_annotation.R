# clean the cells and annotate the data
#
# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. The sample sheet is rebuilt from the pipeline manifests by
# revision/config/build_sample_metadata.R. Changes are marked ## REVISION; the annotation
# logic is unchanged.

library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(scater)
library(tidyr)
library(dplyr)
library(cowplot)
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

REGION <- "VST"

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup(file.path("01-Regional_data_processing", REGION),
            subdirs = "data/02-cell_annotation")

PROCESSED_DATA_PATH = paste0(region_dir(REGION), "/")   ## REVISION: CellBender-corrected pipeline
CELL_META_PATH = "./data/sample_metadata_withMT/"       ## REVISION: stage-01 output, same wd

sample_folders <- list.dirs(PROCESSED_DATA_PATH,recursive = F)
sample_folders_short <- list.dirs(PROCESSED_DATA_PATH,recursive = F,full.names = F)
## REVISION: the corrected pipeline holds only the paper sample set.

## REVISION: the published code read sample_metadata_v3.csv and filtered
## !is.na(cellranger_path), oud_jay_v2 == FALSE, region == "VST".
## load_sample_meta() applies the same filter to the rebuilt sample sheet.
sample_meta_filter <- load_sample_meta(REGION)
table(sample_meta_filter$hiv_jay_v2)
stopifnot(sum(sample_meta_filter$hiv_jay_v2 == TRUE)  == 20,
          sum(sample_meta_filter$hiv_jay_v2 == FALSE) == 13)   ## REVISION: assert the cohort


samples2keep <- sample_meta_filter$sample_name
# change the sample names from HIVOUD to HIV
samples2keep <- sapply(samples2keep,function(sname) gsub("OUD","",sname))
stopifnot(all(samples2keep %in% sample_folders_short))   ## REVISION: assert all kept samples exist (the gsub above is a no-op on the corrected manifests)
samples2load_paths <- paste0(PROCESSED_DATA_PATH,samples2keep)

# load the data into sce
## REVISION: the corrected pipeline drops empty subject-metadata columns
##   per sample, so the per-sample colData column sets disagree and the
##   cbind below fails. align_coldata() restores a single common set.
cell_meta <- lapply(samples2keep, function(sid)
  read.csv(paste0(CELL_META_PATH, sid, ".csv"), row.names = 1))
cell_meta <- align_coldata(cell_meta)
sce.list <- Map(function(sid, cd) SingleCellExperiment(
  list(counts = readRDS(paste0(PROCESSED_DATA_PATH, sid, "/label_transfer/final_counts_", sid, ".rds"))),
  colData = cd), samples2keep, cell_meta)
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




