# prepare the data to sce object and run EdgeR using contrast model
#
# REVISION: run on the CellBender-corrected data; paths come from
# revision/config/paths.R and outputs go to OUT_ROOT. Changes are marked ## REVISION.
# In the revision this published edgeR model is the concordance check; the primary
# cross-region model is donor-aware (revision/tools/cross_region_donor_aware.R).

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

## REVISION: negative-F fallback. edgeR 4.8.2's default glmQLFit returns F < 0
##   (p = 1) for 52 gene x cell type rows in this 90-sample, six-group design,
##   astrocyte SERPINA3 among them (F = -82.9; the published run, older edgeR,
##   F = 108). The counts are not the cause: the Astro pseudobulk is identical to
##   the regional one and matches the published pseudobulk.
##   Every gene keeps the edgeR 4 default QL test; only rows with F < 0 are refit
##   with glmQLFit(legacy = TRUE), or glmLRT where legacy is negative too, and a
##   `test` column records which (pbds_negF_fallback() in revision/config/paths.R).
##   SCORCH_EDGER_QL = fallback (default)
##                   | legacy  (every gene legacy via edger_ql_legacy(); *_legacyQL)
##                   | default (no fallback; *_edgeR4defaultQL)
QL_MODE <- match.arg(tolower(Sys.getenv("SCORCH_EDGER_QL", unset = "fallback")),
                     c("fallback", "legacy", "default"))
edger_ql_legacy(QL_MODE == "legacy")

## REVISION: working directory is set by rerun_setup() under OUT_ROOT instead of a hardcoded server path.
rerun_setup("03-DE_analysis",
            subdirs = c("cross_region/data"))
## REVISION: prefer the combined object with donor 1040's VST library removed
##   and the nine MSN subtypes attached (revision/tools/patch_msn_subtypes.R).
##   The unfiltered object still carries VST.1040 as one of 91 pseudobulk samples.
.combined_candidates <- c(
  Sys.getenv("SCORCH_COMBINED_OBJ", unset = NA_character_),
  "../01-Regional_data_processing/combined/data/seurat_combined_v3_no1040_msn.RDS",
  "../01-Regional_data_processing/combined/data/seurat_combined_v3_no1040.RDS",
  "../01-Regional_data_processing/combined/data/seurat_combined_v3.RDS")
.combined_in <- .combined_candidates[!is.na(.combined_candidates) & file.exists(.combined_candidates)][1]
if (!grepl("no1040", .combined_in))
  warning("[DE] ", basename(.combined_in), " still contains donor 1040's VST cells", immediate. = TRUE)
message("[DE] combined object: ", .combined_in)
seurat.obj <- readRDS(.combined_in)
seurat.obj <- JoinLayers(seurat.obj)

## REVISION: same count-matrix switch as the regional DE (01-03), same
##   default. Primary DE is on uncorrected counts within CellBender-defined
##   populations. Cells here span all three regions, so uncorrected_counts() is
##   called without a region and uses the object's own region column per cell.
DE_COUNTS <- match.arg(Sys.getenv("SCORCH_DE_COUNTS", unset = "uncorrected"),
                       c("uncorrected", "corrected"))
OUT_SFX <- if (DE_COUNTS == "corrected") "_corrected" else ""
## REVISION: the alternative QL modes get their own suffixes, so they never
##   overwrite the default (fallback) results
OUT_SFX <- paste0(OUT_SFX, switch(QL_MODE, fallback = "", legacy = "_legacyQL",
                                  default = "_edgeR4defaultQL"))
## every png()/pdf()/ggsave() in this script picks this up via revision/config/paths.R
options(scorch.out_suffix = OUT_SFX)
message("[DE] count matrix: ", DE_COUNTS,
        if (DE_COUNTS == "uncorrected") "  (PRIMARY)" else "  (CONCORDANCE ONLY)")
if (DE_COUNTS == "uncorrected") {
  count_data <- uncorrected_counts(seurat.obj)
} else {
  count_data <- seurat.obj@assays$RNA@layers$counts
  colnames(count_data) <- colnames(seurat.obj)
  rownames(count_data) <- rownames(seurat.obj)
  if (any(count_data@x != floor(count_data@x))) count_data@x <- round(count_data@x)
}
stopifnot(identical(colnames(count_data), colnames(seurat.obj)),
          identical(rownames(count_data), rownames(seurat.obj)))
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

META_PATH = SAMPLE_META   ## REVISION: rebuilt sample sheet
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
## REVISION: the published names were hardcoded, assuming four race categories.
##   Derived from the matrix so SCORCH_RACE_CODING=collapsed also works; for the
##   published coding the positions are unchanged.
colnames(mm) <- sub("^ei\\$region_HIV", "", colnames(mm))
colnames(mm) <- sub("^ei\\$Race", "Race", colnames(mm))
colnames(mm) <- sub("^ei\\$Sex.*", "Sex", colnames(mm))
colnames(mm) <- sub("^ei\\$Age$", "Age", colnames(mm))
colnames(mm) <- sub("^ei\\$PMI$", "PMI", colnames(mm))
colnames(mm) <- make.names(colnames(mm))
rownames(mm) <- as.character(ei$sample_id)
stopifnot(!anyNA(rownames(mm)), !any(duplicated(rownames(mm))))
message("[DE] design: ", nrow(mm), " samples x ", ncol(mm), " terms")

saveRDS(pb,paste0("./cross_region/data/pseudo_bulk", OUT_SFX, ".rds"))
saveRDS(pb.norm_mean,paste0("./cross_region/data/pseudo_norm_mean_bulk", OUT_SFX, ".rds"))
saveRDS(mm,paste0("./cross_region/data/mm", OUT_SFX, ".rds"))


contrast <- limma::makeContrasts("(INS.HIV-INS.CTR)/3+(PFC.HIV-PFC.CTR)/3+(VST.HIV-VST.CTR)/3", levels = mm)

res.edgeR <- pbDS(pb, design = mm,contrast=contrast, filter="both",min_cells=10,method = "edgeR")
## REVISION: report, then (default mode) apply, the negative-F fallback
.negF <- function(r) sum(vapply(r$table[[1]], function(t) sum(t$F < 0, na.rm = TRUE), numeric(1)))
message("[DE] QL mode: ", QL_MODE, " | rows with F < 0 from pbDS: ", .negF(res.edgeR))
if (QL_MODE == "fallback") {
  res.edgeR <- pbds_negF_fallback(res.edgeR, design = mm, contrast = contrast)
  message("[DE] after fallback: rows with F < 0: ", .negF(res.edgeR))
}
.a <- res.edgeR$table[[1]]$Astro
if (!is.null(.a)) message("[DE] Astro SERPINA3: F = ", signif(.a$F[.a$gene == "SERPINA3"], 3),
                          ", p = ", signif(.a$p_val[.a$gene == "SERPINA3"], 3),
                          if ("test" %in% names(.a)) paste0(", test = ", .a$test[.a$gene == "SERPINA3"]) else "")
res.edgeR_result <- resDS(sce.obj, res.edgeR, frq = TRUE)

saveRDS(res.edgeR_result,paste0("./cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq", OUT_SFX, ".RDS"))

colnames(res.edgeR_result) <- sapply(colnames(res.edgeR_result),function(cname){
  num <- extract_numeric(cname)
  if(is.na(num)) return(cname)
  else return(gsub(num, paste0("a",num), cname))
})

## REVISION: the published six regional frequency averages were positional
##   c_across() ranges (e.g. VST.10280.frq:VST.1059.frq), which silently break when
##   the sample set changes, as VST's did. add_group_frq_by_region() in
##   revision/config/paths.R maps columns by group_id and skips the group-level
##   columns; on the unchanged sample set it reproduces the positional result exactly.
res.edgeR_result <- add_group_frq_by_region(res.edgeR_result, sce.obj)
saveRDS(res.edgeR_result,paste0("./cross_region/data/all_HIV-CTR_res_edgeR_filtered_withFreq", OUT_SFX, ".RDS"))

