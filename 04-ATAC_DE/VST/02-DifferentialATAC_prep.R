# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. The published script is restructured with the same
# design, pbDS call and outputs as the INS and PFC scripts. Changes:
#   1. Group frequency means are computed by condition (add_group_frq()), not over a
#      positional c_across() column range, which encodes one sample ordering and
#      averages the wrong samples, without error, if that ordering changes.
#   2. Design column names are derived from the model matrix and their count is
#      checked, instead of a hard-coded list that assumes a fixed set of race levels;
#      a mislabelled design column silently changes what makeContrasts() builds.
#   3. Reads the object carrying the nine MSN subtypes
#      (revision/tools/patch_msn_subtypes.R --target atac). Without them every
#      striatal neuron is one "MSNs" cluster; the script falls back to the unpatched
#      object with a warning.
# Usage: revision/docker/run.sh 04-ATAC_DE/VST/02-DifferentialATAC_prep.R

library(Seurat)
library(Signac)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(scater)
library(tidyr)
library(dplyr)
library(cowplot)
library(limma)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

rerun_setup("04-ATAC_DE", subdirs = c("VST/data_ATAC"))

## ---- input -------------------------------------------------------------------
.cands <- c(
  "../02-ATAC_processing/sampleQC_redo/VST/data/signac_integrated_v2_msn.RDS",
  "../02-ATAC_processing/sampleQC_redo/VST/data/signac_integrated_v2.RDS"
)
IN <- .cands[file.exists(.cands)][1]
if (is.na(IN)) stop("no VST ATAC object found")
if (!grepl("_msn", IN))
  warning("[VST ATAC] ", basename(IN), " has no MSN subtypes - the DE will test ",
          "one 'MSNs' cluster instead of nine. Run ",
          "revision/tools/patch_msn_subtypes.R --target atac first.", immediate. = TRUE)
message("[VST ATAC] input: ", IN)

signac.obj.VST <- readRDS(IN)

.ct <- as.character(signac.obj.VST@meta.data$consensus_annotation_v3)
message("[VST ATAC] ", ncol(signac.obj.VST), " cells, ",
        length(unique(.ct)), " cell types, ",
        sum(grepl("^D1|^D2", unique(.ct))), " MSN subtypes")

DefaultAssay(signac.obj.VST) <- "ATAC"
count_data <- signac.obj.VST@assays$ATAC$counts
colnames(count_data) <- colnames(signac.obj.VST)
rownames(count_data) <- rownames(signac.obj.VST)
sce.obj <- SingleCellExperiment(list(counts = count_data),
                                colData = signac.obj.VST@meta.data)

sce.obj <- prepSCE(sce.obj,
                   kid = "consensus_annotation_v3", # subpopulation assignments
                   gid = "condition_new",  # group IDs (ctrl/stim)
                   sid = "sample_ids",   # sample IDs (ctrl/stim.1234)
                   drop = FALSE)  # drop all other colData columns

pb.vst <- aggregateData(sce.obj, assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
sce.obj <- computeLibraryFactors(sce.obj)
sce.obj <- logNormCounts(sce.obj)
pb.vst.norm_mean <- aggregateData(sce.obj, assay = "logcounts", fun = "mean",
                                  by = c("cluster_id", "sample_id"))

## ---- design ------------------------------------------------------------------
ei <- metadata(sce.obj)$experiment_info
ei$group_id <- factor(ei$group_id, levels = c("HIV", "CTR"))

sample_meta <- read.csv(SAMPLE_META)
# only load samples of interest:
# HIV and CTR samples from INS
sample_meta_filter <- as.data.frame(sample_meta) %>%
  dplyr::filter(!is.na(cellranger_path)) %>%
  dplyr::filter((oud_jay_v2 == FALSE) & (region == "VST"))

sample_meta_filter_selected <- sample_meta_filter %>%
  dplyr::select(sample_ids, age_value_v2, sex_v2, race_v2,
                postmortem_interval_v2, hiv_jay_v2)
sample_meta_filter_selected$sample_id <- sample_meta_filter_selected$sample_ids
ei <- ei %>% left_join(sample_meta_filter_selected, by = "sample_id")

## REVISION: stop if a sample in the object is missing from the sample sheet.
if (anyNA(ei$age_value_v2) || anyNA(ei$postmortem_interval_v2))
  stop("sample(s) in the object missing from the sample sheet: ",
       paste(ei$sample_id[is.na(ei$age_value_v2)], collapse = ", "))

ei$Age  <- scale(as.numeric(ei$age_value_v2))[, 1]
ei$Sex  <- factor(ei$sex_v2, levels = c("male", "female"))
ei$Race <- make.names(ei$race_v2)
### ~PMI (scaled)
ei$PMI <- scale(as.numeric(ei$postmortem_interval_v2))[, 1]
ei$HIV  <- ei$group_id

message("[VST ATAC] design samples: ", nrow(ei),
        "  (", sum(ei$HIV == "HIV"), " PWH / ", sum(ei$HIV == "CTR"), " PWoH)")
message("[VST ATAC] race levels: ", paste(sort(unique(ei$Race)), collapse = ", "))

mm <- model.matrix(~ 0 + ei$HIV + ei$Age + ei$Sex + ei$Race + ei$PMI)

## REVISION: derive the column names rather than hard-coding them (header note 2).
## The order model.matrix() produces is: the HIV levels, then Age, then Sex,
## then the Race dummies, then PMI.
.race_lv <- sort(unique(ei$Race))
.race_dummies <- .race_lv[-1]                      # first level is the reference
.abbrev <- function(x) gsub("[^A-Za-z]", "", gsub("(?<=\\b)([A-Za-z])[a-z]*", "\\1",
                                                  gsub("[._]", " ", x), perl = TRUE))
.nm <- c(levels(ei$HIV), "Age", "Sex", paste0("Race", .abbrev(.race_dummies)), "PMI")
if (length(.nm) != ncol(mm))
  stop("design name count (", length(.nm), ") != model matrix columns (", ncol(mm),
       "): ", paste(colnames(mm), collapse = " | "))
dimnames(mm) <- list(ei$sample_id, .nm)
message("[VST ATAC] design: ", paste(.nm, collapse = " + "))

# save the data for other comparisons
saveRDS(pb.vst,          "./VST/data_ATAC/pseudo_bulk.rds")
saveRDS(pb.vst.norm_mean, "./VST/data_ATAC/pseudo_norm_mean_bulk.rds")
saveRDS(mm,               "./VST/data_ATAC/mm.rds")

## ---- DE ----------------------------------------------------------------------
contrast <- makeContrasts("HIV-CTR", levels = mm)

res.edgeR <- pbDS(pb.vst, design = mm, contrast = contrast,
                  filter = "both", min_cells = 10, method = "edgeR")
res.edgeR_result <- resDS(sce.obj, res.edgeR, frq = TRUE)

colnames(res.edgeR_result) <- sapply(colnames(res.edgeR_result), function(cname) {
  num <- extract_numeric(cname)
  if (is.na(num)) return(cname) else return(gsub(num, paste0("a", num), cname))
})

## ---- group frequency means, by condition not by column position ---------------
## REVISION: add_group_frq() (revision/config/paths.R) maps each per-sample
## <sample>.frq column (with the "a" prefix applied above to ids starting with a
## digit) to its condition, and skips the group-level CTR.frq / HIV.frq columns.
res.edgeR_result <- add_group_frq(res.edgeR_result, sce.obj, "VST")

saveRDS(res.edgeR_result, "./VST/data_ATAC/VST_HIV-CTR_res_edgeR_filtered_withFreq.RDS")

message("\n[VST ATAC] ", nrow(res.edgeR_result), " rows; ",
        length(unique(res.edgeR_result$cluster_id)), " cell types tested")
message("[VST ATAC] significant (p_adj.glb < 0.05): ",
        sum(res.edgeR_result$p_adj.glb < 0.05, na.rm = TRUE))
message("[VST ATAC] wrote ./VST/data_ATAC/VST_HIV-CTR_res_edgeR_filtered_withFreq.RDS")
