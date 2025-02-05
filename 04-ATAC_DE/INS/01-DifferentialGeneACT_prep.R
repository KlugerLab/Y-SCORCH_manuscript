library(Seurat)
library(Signac)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(scater)
library(tidyr)
library(dplyr)
library(cowplot)

wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/04-ATAC_DE/"
setwd(wd)


signac.obj.INS <- readRDS("../02-ATAC_processing/sampleQC_redo/INS/data/signac_integrated_v2.RDS")

DefaultAssay(signac.obj.INS) 
count_data <- signac.obj.INS@assays$GeneACT@counts
colnames(count_data) <- colnames(signac.obj.INS)
rownames(count_data) <- rownames(signac.obj.INS)
sce.obj <- SingleCellExperiment(list(counts = count_data),colData = signac.obj.INS@meta.data)

sce.obj <- prepSCE(sce.obj, 
                   kid = "consensus_annotation_v3", # subpopulation assignments
                   gid = "condition_new",  # group IDs (ctrl/stim)
                   sid = "sample_ids",   # sample IDs (ctrl/stim.1234)
                   drop = F)  # drop all other colData columns


pb.ins <- aggregateData(sce.obj,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
sce.obj <- computeLibraryFactors(sce.obj)
sce.obj <- logNormCounts(sce.obj)

pb.ins.norm_mean <- aggregateData(sce.obj,
                                  assay = "logcounts", fun = "mean",
                                  by = c("cluster_id", "sample_id"))


ei <- metadata(sce.obj)$experiment_info
ei$group_id <- factor(ei$group_id,levels=c("HIV","CTR"))

META_PATH = "/banach2/SCORCH/data/metadata/sample_metadata_v3.csv"
sample_meta <- read.csv(META_PATH)
# only load samples of interest:
# HIV and CTR samples from INS
sample_meta_filter = as.data.frame(sample_meta) %>%
  dplyr::filter(!is.na(cellranger_path)) %>%
  dplyr::filter((oud_jay_v2 == FALSE) & (region=="INS") ) 
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
saveRDS(pb.ins,"./INS/data_GeneACT/pseudo_bulk.rds")
saveRDS(pb.ins.norm_mean,"./INS/data_GeneACT/pseudo_norm_mean_bulk.rds")
saveRDS(mm,"./INS/data_GeneACT/mm.rds")



contrast <- makeContrasts("HIV-CTR", levels = mm)

res.edgeR <- pbDS(pb.ins, design = mm,contrast=contrast, filter="both",min_cells=10,method = "edgeR")
res.edgeR_result <- resDS(sce.obj, res.edgeR, frq = TRUE)

colnames(res.edgeR_result) <- sapply(colnames(res.edgeR_result),function(cname){
  num <- extract_numeric(cname)
  if(is.na(num)) return(cname)
  else return(gsub(num, paste0("a",num), cname))
})

res.edgeR_result <- res.edgeR_result %>%
  rowwise() %>%
  mutate(avg_CTR_frq = mean(c_across(HCCCY.frq:HCTME.frq), na.rm = TRUE)) %>%
  mutate(avg_HIV_frq = mean(c_across(a1063.frq:a1040.frq), na.rm = TRUE)) %>%
  mutate(region = "INS") %>%
  ungroup()

saveRDS(res.edgeR_result,"./INS/data_GeneACT/INS_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
