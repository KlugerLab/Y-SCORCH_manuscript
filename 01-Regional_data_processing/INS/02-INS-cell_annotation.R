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

REGION <- "INS"

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup(file.path("01-Regional_data_processing", REGION),
            subdirs = "data/02-cell_annotation")

PROCESSED_DATA_PATH = paste0(region_dir(REGION), "/")   ## REVISION: CellBender-corrected pipeline
CELL_META_PATH = "./data/sample_metadata_withMT/"       ## REVISION: stage-01 output, same wd

sample_folders <- list.dirs(PROCESSED_DATA_PATH,recursive = F)
sample_folders_short <- list.dirs(PROCESSED_DATA_PATH,recursive = F,full.names = F)
## REVISION: the corrected pipeline contains only the paper sample set (31 for INS);
## the published pipeline folder held 53 including OUD samples.

## REVISION: the published code read sample_metadata_v3.csv and filtered
## !is.na(cellranger_path), oud_jay_v2 == FALSE, region == "INS".
## load_sample_meta() applies the same filter to the rebuilt sample sheet.
sample_meta_filter <- load_sample_meta(REGION)
table(sample_meta_filter$hiv_jay_v2)
# 12 CTR vs 19 HIV
stopifnot(sum(sample_meta_filter$hiv_jay_v2 == TRUE)  == 19,
          sum(sample_meta_filter$hiv_jay_v2 == FALSE) == 12)   ## REVISION: assert the cohort


samples2keep <- sample_meta_filter$sample_name
samples2keep %in% sample_folders_short
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


sce.obj@colData[,"condition_new"] <- sapply(sce.obj@colData$sample_name,function(sid) sample_meta_filter[sample_meta_filter$sample_name == sid,"hiv_jay_v2"])
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
sce.obj@colData$biccn_predicted.id <- as.factor(sce.obj@colData$biccn_predicted.id)
sce.obj@colData$ma_predicted.id <- as.factor(sce.obj@colData$ma_predicted.id)


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

# remove cells with max scores < 0.5
metadata_all <- as.data.frame(sce.obj@colData) %>%
  mutate(cell2remove = ifelse((biccn_prediction_score.max < 0.5) & (ma_prediction_score.max < 0.5),T,F))
sce.obj@colData$cell2remove_score <- metadata_all$cell2remove


# 3. cell annotation
# add EXT, INH, Glial annotations

cp_colors <-  c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1", "deeppink1", 
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3","yellow1"
  
)
biccn_cp_unique_ordered <- c("L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP",
                             "L6 CT","L6 IT","L6 IT Car3","L6b",
                             "Lamp5","Lamp5 Lhx6","Chandelier","Sncg","Pvalb",
                             "Sst","Sst Chodl","Pax6","Vip",
                             "Astro","Endo","Micro/PVM","Oligo","OPC","VLMC" )

sestan_unique_ordered <- c("L2-3 IT","L3-5 IT-1","L3-5 IT-2","L3-5 IT-3","L5 ET","L5-6 NP",
                           "L6 CT","L6 IT-1","L6 IT-2","L6B",
                           "LAMP5 RELN","LAMP5 LHX6","PVALB ChC","ADARB2 KCNG1","PVALB",
                           "SST","SST HGF","SST NPY","VIP",
                           "Astro","Endo","Micro","PC","Immune","Oligo","OPC",
                           "VLMC","SMC","RB")
all_cp_unique_ordered <- c("L2/3 IT","L2-3 IT",
                           "L4 IT","L5 IT","L3-5 IT-1","L3-5 IT-2","L3-5 IT-3",
                           "L5 ET","L5/6 NP","L5-6 NP",
                           "L6 CT","L6 IT","L6 IT Car3","L6 IT-1","L6 IT-2","L6b","L6B",
                           "Lamp5","LAMP5 RELN","Lamp5 Lhx6","LAMP5 LHX6",
                           "Chandelier","PVALB ChC","Pax6","ADARB2 KCNG1","Pvalb","PVALB",
                           "Sncg","Sst","SST","Sst Chodl","SST HGF","SST NPY",
                           "Vip","VIP",
                           "Astro","Endo","Oligo","OPC","Micro/PVM","Micro","Immune","VLMC","PC",
                           "SMC","RB")
coarse_cp_mapper <- c(rep("Ext",17),rep("Inh",18),rep("Glial",11))
names(coarse_cp_mapper) <- all_cp_unique_ordered

sce.obj@colData$biccn_predicted.id_coarse <- coarse_cp_mapper[as.character(sce.obj@colData$biccn_predicted.id)]
sce.obj@colData$ma_predicted.id_coarse <- coarse_cp_mapper[as.character(sce.obj@colData$ma_predicted.id)]
sce.obj@colData$cell2remove_consistency <- sce.obj@colData$biccn_predicted.id_coarse != sce.obj@colData$ma_predicted.id_coarse

# note: remove cells that are either have low prediction score or
#       are inconsistent between the two set of annotations (coarse)
sce.obj@colData$cell2remove <- sce.obj@colData$cell2remove_score | sce.obj@colData$cell2remove_consistency

sce.obj <- sce.obj[,!sce.obj@colData$cell2remove]
sce.obj@colData$consensus_cell_type_coarse <- sce.obj@colData$biccn_predicted.id_coarse


# combine the annotations
# if neurons, use biccn, else ma et al
metadata_all <- as.data.frame(sce.obj@colData) 
metadata_all <- metadata_all %>%
  mutate(consensus_annotation_v2 = ifelse(consensus_cell_type_coarse %in% c("Ext","Inh"),as.character(biccn_predicted.id), as.character(ma_predicted.id)))
sce.obj@colData$consensus_annotation_v2 <- metadata_all$consensus_annotation_v2


cp_col2 <- c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1", "deeppink1", 
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
cell_types_consensus_levels <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Immune","PC","Endo","SMC","VLMC")


names(cp_col2) <- cell_types_consensus_levels
#remove NA (RB)
sce.obj@colData$consensus_annotation_v2 <- factor(sce.obj@colData$consensus_annotation_v2,levels = cell_types_consensus_levels)
sce.obj <- sce.obj[,!is.na(sce.obj@colData$consensus_annotation_v2)]


saveRDS(sce.obj,"./data/02-cell_annotation/processed_sce_obj.rds")
