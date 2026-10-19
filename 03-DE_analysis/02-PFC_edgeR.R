# prepare the data to sce object and run DEG analysis
#
# REVISION: run on the CellBender-corrected data; paths come from
# revision/config/paths.R and outputs go to OUT_ROOT. Changes are marked ## REVISION.

library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(muscat)
library(scater)
library(tidyr)
library(dplyr)
library(cowplot)
## REVISION: limma (makeContrasts) is used below but was not loaded.
library(limma)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set by rerun_setup() under OUT_ROOT instead of a hardcoded server path.
rerun_setup("03-DE_analysis",
            subdirs = c("PFC/data", "PFC/results"))
seurat.obj <- readRDS("../01-Regional_data_processing/PFC/data/03-cell_integration_vis/integrated_seurat_obj_v3.RDS")

REGION_DE <- "PFC"
## REVISION: choice of count matrix for the test.
##   SCORCH_DE_COUNTS = "uncorrected" (default) | "corrected"
##
##   Primary DE runs on uncorrected counts. CellBender defines the cell
##   populations (QC, clustering, annotation); the counts fed to the test are
##   the CellRanger raw counts for those same cells. CellBender over-corrects
##   the reactive-astrocyte panel (SERPINA3, OSMR, CHI3L1) by zeroing genuine
##   low-level expression in controls, which creates artifactual complete
##   separation: edgeR's QL F-test returns a negative F and the effect size is
##   overstated.
##
##   The corrected run is a concordance check, not the source of reported
##   magnitudes. It writes to *_corrected filenames so it cannot overwrite the
##   primary result, and rounds counts to integers.
DE_COUNTS <- match.arg(Sys.getenv("SCORCH_DE_COUNTS", unset = "uncorrected"),
                       c("uncorrected", "corrected"))
OUT_SFX <- if (DE_COUNTS == "corrected") "_corrected" else ""
OUT_SFX <- paste0(OUT_SFX, if (identical(Sys.getenv("SCORCH_RACE_CODING"), "collapsed")) "_raceCollapsed" else "")
message("[DE] count matrix: ", DE_COUNTS,
        if (DE_COUNTS == "corrected") "  (CONCORDANCE ONLY - not for reported magnitudes)" else "  (PRIMARY)")

if (DE_COUNTS == "uncorrected") {
  count_data <- uncorrected_counts(seurat.obj, REGION_DE)
} else {
  count_data <- seurat.obj@assays$RNA@layers$counts
  colnames(count_data) <- colnames(seurat.obj)
  rownames(count_data) <- rownames(seurat.obj)
  ## CellBender posterior counts can be fractional; the QL F-test requires
  ## integers.
  if (any(count_data@x != floor(count_data@x))) {
    message("[DE] corrected counts are fractional - rounding (", 
            sum(count_data@x != floor(count_data@x)), " non-integer entries)")
    count_data@x <- round(count_data@x)
  }
}
stopifnot(identical(colnames(count_data), colnames(seurat.obj)),
          identical(rownames(count_data), rownames(seurat.obj)))
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

META_PATH = SAMPLE_META   ## REVISION: rebuilt sample sheet
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

## REVISION: race coding.
##   SCORCH_RACE_CODING = "published" (default) | "collapsed"
##
##   The published coding has four categories, and the two rarest are
##   represented by 1 and 2 donors, all PWH. With no controls in either, those
##   coefficients are estimated from donors nested entirely inside the
##   exposure group and are not identified by any within-category contrast.
##
##   This is also what silently removes cell types: pbDS() drops a cluster
##   without warning when the design loses rank
##   (muscat: `qr(d)$rank < ncol(d)` -> return(NULL)), and a one-donor race
##   column collapses to all zeros as soon as that donor falls below
##   min_cells = 10. In PFC this costs 11 cell types including L6 IT,
##   L6 IT Car3 and L6b, in the published analysis as well as here.
##
##   "collapsed" merges the two Native Alaskan/American Indian categories into a
##   single Other level (n = 3), which keeps the adjustment but stops one donor
##   from gating a whole cluster. It deviates from the published model and is
##   reported as a sensitivity analysis, not as the primary result.
RACE_CODING <- match.arg(Sys.getenv("SCORCH_RACE_CODING", unset = "published"),
                         c("published", "collapsed"))
if (RACE_CODING == "collapsed") {
  .rare <- grepl("Native", ei$Race)
  ei$Race[.rare] <- "Other"
  message("[DE] race coding: collapsed — ", sum(.rare),
          " Native Alaskan/American Indian donor(s) merged into 'Other'")
} else {
  message("[DE] race coding: published (4 categories, as in the paper)")
}
message("[DE] race levels: ",
        paste(sprintf("%s=%d", names(table(ei$Race)), table(ei$Race)), collapse = "  "))
### ~PMI (scaled)
ei$PMI <- as.numeric(ei$postmortem_interval_v2)
ei$PMI <- scale(ei$PMI)[,1]
ei$HIV <- ei$group_id 

## REVISION: Age, Sex, Race and PMI are joined from the sample sheet, and a
##   single NA in any of them makes model.matrix() silently drop that donor
##   from the design, so the DE would run on fewer samples than the cohort with
##   no error or warning. Assert completeness here rather than trusting the join.
.cov <- list(Age = ei$Age, Sex = ei$Sex, Race = ei$Race, PMI = ei$PMI, HIV = ei$HIV)
.bad <- vapply(.cov, function(v) sum(is.na(v)), integer(1))
if (any(.bad > 0))
  stop("covariate NA after joining the sample sheet - model.matrix would drop ",
       "these donors from the design: ",
       paste(sprintf("%s=%d", names(.bad)[.bad > 0], .bad[.bad > 0]), collapse = ", "))
if (nlevels(droplevels(factor(ei$Sex))) < 2)
  stop("Sex has <2 levels after the join; the design is rank deficient")
message("[DE] design covariates complete for ", nrow(ei), " donors: ",
        "HIV + Age + Sex + Race + PMI")
mm <- model.matrix(~ 0 +ei$HIV+ei$Age+ei$Sex +ei$Race+ei$PMI )
## REVISION: the published code hardcoded eight column names, which assumes
##   exactly four race categories. Under SCORCH_RACE_CODING=collapsed there are
##   three, so the hardcoded vector would not match. The names are derived from
##   the matrix instead; for the published coding this gives the original names.
colnames(mm) <- sub("^ei\\$HIV", "", colnames(mm))
colnames(mm) <- sub("^ei\\$Race", "Race", colnames(mm))
colnames(mm) <- sub("^ei\\$Sex.*", "Sex", colnames(mm))
colnames(mm) <- sub("^ei\\$Age$", "Age", colnames(mm))
colnames(mm) <- sub("^ei\\$PMI$", "PMI", colnames(mm))
colnames(mm) <- make.names(colnames(mm))
rownames(mm) <- as.character(ei$sample_id)
stopifnot(all(levels(ei$HIV) %in% colnames(mm)),
          !anyNA(rownames(mm)), !any(duplicated(rownames(mm))))
message("[DE] design: ", nrow(mm), " donors x ", ncol(mm), " terms  (",
        paste(colnames(mm), collapse = ", "), ")")

saveRDS(pb.pfc,paste0("./PFC/data/pseudo_bulk", OUT_SFX, ".rds"))
saveRDS(pb.pfc.norm_mean,paste0("./PFC/data/pseudo_norm_mean_bulk", OUT_SFX, ".rds"))
saveRDS(mm,paste0("./PFC/data/mm", OUT_SFX, ".rds"))


contrast <- makeContrasts("HIV-CTR", levels = mm)

res.edgeR <- pbDS(pb.pfc, design = mm,contrast=contrast, filter="both",min_cells=10,method = "edgeR")
## REVISION: pbDS() drops a cell type silently when the design loses rank for
##   it (muscat: `qr(d)$rank < ncol(d)` -> return(NULL), then the NULLs are
##   removed with no message). A one-donor covariate level collapses to all
##   zeros as soon as that donor falls below min_cells, so whole cell types
##   vanish from the results and read as "no DE" rather than "not estimable".
##   In the published race coding this costs PFC 11 cell types. Log them.
.tested  <- names(res.edgeR$table[[1]])
## REVISION: this script names the pseudobulk `pb.pfc`, not `pb`.
.dropped <- setdiff(SummarizedExperiment::assayNames(pb.pfc), .tested)
message("[DE] cell types tested: ", length(.tested), " of ",
        length(SummarizedExperiment::assayNames(pb.pfc)))
if (length(.dropped))
  message("[DE] NOT ESTIMABLE (dropped by pbDS, design rank/min_cells): ",
          paste(.dropped, collapse = ", "))

res.edgeR_result <- resDS(sce.obj, res.edgeR, frq = TRUE)

saveRDS(res.edgeR_result,paste0("./PFC/results/PFC_HIV-CTR_res_edgeR_filtered_withFreq", OUT_SFX, ".RDS"))

res.edgeR_result <- readRDS(paste0("./PFC/results/PFC_HIV-CTR_res_edgeR_filtered_withFreq", OUT_SFX, ".RDS"))

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

saveRDS(res.edgeR_result,paste0("./PFC/results/PFC_HIV-CTR_res_edgeR_filtered_withFreq", OUT_SFX, ".RDS"))

