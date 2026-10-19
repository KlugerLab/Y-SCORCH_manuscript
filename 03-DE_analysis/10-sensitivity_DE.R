#!/usr/bin/env Rscript
# ==============================================================================
# Sensitivity and regrouping variants of the regional pseudobulk DE.
#
# Added for the revision (not part of the published analysis). Reruns the
# regional edgeR DE under a different sample set or grouping, for two reviewer
# requests:
#
#  VST exclusions: how much the VST findings depend on the excluded samples
#    -> run VST DE with and without the excluded samples and compare.
#
#  Viremia / HIV encephalitis (HIVE): a three-way suppressed vs viremic vs
#  control analysis, excluding the HIVE donor.
#
#    HIVE drop -> --exclude 1039. This is the primary sensitivity analysis.
#    1039 is the only HIVE donor in the cohort (brain_pathology_v2 =
#    "HIV Encephalitis, Progressive multifocal leukoencephalopathy"), so it is a
#    clean leave-one-out defined by neuropathology rather than by a lab value.
#
#    Three-way (--group viremia) is exploratory only. The cohort is clinically
#    suppressed and the recorded plasma titres are not reliable for stratification:
#      - "viremic" tracks ART status, not biology. Of the donors above 10,000
#        copies/mL, 1068, 1039 and 7200776574 are no_art, 1119 is on 3TC
#        monotherapy, and 1133 has no ART recorded.
#      - the measurements are not contemporaneous with death:
#        months_since_last_plasma_viral_measurement ranges 0 to 13.1.
#      - 2 of 20 PWH have no recorded titre at all (1137, 7100518287).
#    Splitting on VL would conflate treatment status, adherence, measurement
#    timing and terminal course, then attribute the result to viremia.
#
# Method: the same recipe as 01/02/03-*_edgeR.R:
# prepSCE -> aggregateData -> model.matrix(~ 0 + <group> + Age + Sex + Race + PMI)
# -> pbDS(method = "edgeR", filter = "both", min_cells = 10). Only the sample set
# or the grouping variable changes.
#
# Inputs:  OUT_ROOT/01-Regional_data_processing/<REGION>/data/.../integrated_seurat_obj_v3.RDS
#          and the rebuilt sample sheet (SAMPLE_META).
# Outputs: OUT_ROOT/03-DE_analysis/<REGION>/results/<REGION>_<label>_res_edgeR_filtered_withFreq.RDS,
#          <REGION>_<label>_runinfo.txt, and pseudobulk / design matrices under
#          OUT_ROOT/03-DE_analysis/<REGION>/data/.
#
# Usage (inside the container: revision/docker/run.sh 03-DE_analysis/10-sensitivity_DE.R [args]):
#   Rscript 10-sensitivity_DE.R --region INS --exclude 1039 --label noHIVE
#   Rscript 10-sensitivity_DE.R --region VST --keep-all-samples --label allVST
#   Rscript 10-sensitivity_DE.R --region INS --group viremia --label threeway
#   Rscript 10-sensitivity_DE.R --region INS --group viremia --exclude 1039 \
#                               --label threeway_noHIVE
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(muscat)
  library(scater);  library(edgeR); library(limma)
  library(dplyr);   library(tidyr); library(optparse)
})

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--region",  type = "character", default = "INS"),
  make_option("--exclude", type = "character", default = "",
              help = "comma-separated sample_ids to drop, e.g. 1039"),
  make_option("--group",   type = "character", default = "HIV",
              help = "HIV (two-group, as published) or viremia (three-group)"),
  make_option("--label",   type = "character", default = "sensitivity",
              help = "suffix for the output files"),
  make_option("--vl-threshold", type = "double", default = 200,
              help = "plasma copies/mL at or above which a donor counts as viremic [default 200]"),
  make_option("--keep-all-samples", action = "store_true", default = FALSE,
              help = "VST only: keep the samples the published run excluded")
)))

REGION  <- toupper(opt$region)
EXCLUDE <- trimws(unlist(strsplit(opt$exclude, ",")))
EXCLUDE <- EXCLUDE[nzchar(EXCLUDE)]
LABEL   <- opt$label

rerun_setup("03-DE_analysis",
            subdirs = c(file.path(REGION, "data"), file.path(REGION, "results")))

# ---- input: the same object the published regional DE reads --------------------
IN <- sprintf("../01-Regional_data_processing/%s/data/%s/integrated_seurat_obj_v3.RDS",
              REGION, if (REGION == "VST") "04-VST-tsne_marker_plotting" else "03-cell_integration_vis")
if (!file.exists(IN)) stop("missing input object: ", IN, "\n  run the regional stage-04 step first")
seurat.obj <- readRDS(IN)
message("[sens] loaded ", IN, ": ", ncol(seurat.obj), " cells")

# ---- grouping ------------------------------------------------------------------
meta <- load_sample_meta(REGION)

parse_vl <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- suppressWarnings(as.numeric(x))
  # "below_det_limit_50_copies" / "below_det_limit_200_copies" -> treat as suppressed
  out[grepl("below_det", x)] <- 0
  # "above_det_limit_750000_copies" -> take the stated bound
  hi <- grepl("above_det", x)
  out[hi] <- suppressWarnings(as.numeric(gsub("\\D+", "", x[hi])))
  out
}

if (opt$group == "viremia") {
  warning("--group viremia is EXPLORATORY. The cohort is clinically suppressed ",
          "and the recorded titres track ART status and measurement timing more ",
          "than biology - see the header. Do not report as a headline analysis.")
  meta$vl_numeric <- parse_vl(meta$plasma_viral_load_measurement_v2)
  meta$group <- ifelse(meta$hiv_jay_v2 %in% FALSE, "CTR",
                ifelse(is.na(meta$vl_numeric), NA_character_,
                ifelse(meta$vl_numeric >= opt$`vl-threshold`, "VIREMIC", "SUPPRESSED")))
  dropped_vl <- meta$sample_ids[is.na(meta$group)]
  if (length(dropped_vl))
    message("[sens] no viral load, dropped: ", paste(dropped_vl, collapse = ", "))
  meta <- meta[!is.na(meta$group), , drop = FALSE]
} else {
  meta$group <- ifelse(meta$hiv_jay_v2 %in% TRUE, "HIV", "CTR")
}

# ---- sample set ----------------------------------------------------------------
if (length(EXCLUDE)) {
  n0 <- nrow(meta)
  meta <- meta[!(meta$sample_ids %in% EXCLUDE), , drop = FALSE]
  message("[sens] excluded ", n0 - nrow(meta), " sample(s): ", paste(EXCLUDE, collapse = ", "))
}
keep_ids <- as.character(meta$sample_ids)
message("[sens] cohort: ", nrow(meta), " samples")
print(table(meta$group))

seurat.obj <- subset(seurat.obj,
                     cells = colnames(seurat.obj)[
                       as.character(seurat.obj@meta.data$sample_ids) %in% keep_ids])
message("[sens] cells after subsetting: ", ncol(seurat.obj))

## carry the grouping onto the cells
grp <- setNames(meta$group, as.character(meta$sample_ids))
seurat.obj@meta.data$condition_new <- unname(grp[as.character(seurat.obj@meta.data$sample_ids)])
stopifnot(!any(is.na(seurat.obj@meta.data$condition_new)))

# ---- identical to the published DE from here ----------------------------------
## Same count-matrix switch as 01-03, and the same default (uncorrected counts).
##   A sensitivity analysis is only interpretable against the primary DE if both
##   use the same counts; otherwise a difference could come from the cohort
##   change or from the count matrix and the two cannot be separated.
DE_COUNTS <- match.arg(Sys.getenv("SCORCH_DE_COUNTS", unset = "uncorrected"),
                       c("uncorrected", "corrected"))
message("[sens] count matrix: ", DE_COUNTS,
        if (DE_COUNTS == "uncorrected") "  (matches the primary DE)" else
          "  (CONCORDANCE ONLY)")
if (DE_COUNTS == "uncorrected") {
  count_data <- uncorrected_counts(seurat.obj, opt$region)
} else {
  count_data <- seurat.obj@assays$RNA@layers$counts
  colnames(count_data) <- colnames(seurat.obj)
  rownames(count_data) <- rownames(seurat.obj)
  if (any(count_data@x != floor(count_data@x))) count_data@x <- round(count_data@x)
}
stopifnot(identical(colnames(count_data), colnames(seurat.obj)),
          identical(rownames(count_data), rownames(seurat.obj)))
sce.obj <- SingleCellExperiment(list(counts = count_data), colData = seurat.obj@meta.data)

sce.obj <- prepSCE(sce.obj, kid = "consensus_annotation_v3",
                   gid = "condition_new", sid = "sample_ids", drop = FALSE)

pb <- aggregateData(sce.obj, assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
sce.obj <- computeLibraryFactors(sce.obj)
sce.obj <- logNormCounts(sce.obj)
pb.norm_mean <- aggregateData(sce.obj, assay = "logcounts", fun = "mean",
                              by = c("cluster_id", "sample_id"))

ei <- metadata(sce.obj)$experiment_info
sel <- meta %>%
  dplyr::select(sample_ids, age_value_v2, sex_v2, race_v2,
                postmortem_interval_v2, group) %>%
  dplyr::mutate(sample_id = as.character(sample_ids))
ei <- ei %>% dplyr::mutate(sample_id = as.character(sample_id)) %>%
  dplyr::left_join(sel, by = "sample_id")

ei$Age   <- as.numeric(ei$age_value_v2)
ei$Sex   <- make.names(ei$sex_v2)
ei$Race  <- make.names(ei$race_v2)
ei$PMI   <- as.numeric(ei$postmortem_interval_v2)
ei$GROUP <- factor(make.names(ei$group))
stopifnot(!any(is.na(ei$Age)), !any(is.na(ei$PMI)))

mm <- model.matrix(~ 0 + ei$GROUP + ei$Age + ei$Sex + ei$Race + ei$PMI)
colnames(mm) <- gsub("^ei\\$GROUP", "", colnames(mm))
colnames(mm) <- make.names(colnames(mm))

## The design needs sample rownames, not just column names. pbDS() subsets it
##   as design[colnames(pb), ] (pseudobulk columns are sample ids), so a design
##   with default integer rownames fails with "subscript out of bounds" inside
##   BiocParallel, where the cause is hard to see.
rownames(mm) <- as.character(ei$sample_id)
stopifnot(!anyNA(rownames(mm)), !any(duplicated(rownames(mm))))
message("[sens] design: ", nrow(mm), " donors x ", ncol(mm), " terms  (",
        paste(colnames(mm), collapse = ", "), ")")

## contrasts: the published HIV-CTR, or all pairwise for the three-way
contrast <- if (opt$group == "viremia") {
  limma::makeContrasts(
    "VIREMIC-CTR"        = VIREMIC - CTR,
    "SUPPRESSED-CTR"     = SUPPRESSED - CTR,
    "VIREMIC-SUPPRESSED" = VIREMIC - SUPPRESSED,
    levels = mm)
} else {
  limma::makeContrasts("HIV-CTR" = HIV - CTR, levels = mm)
}

res.edgeR <- pbDS(pb, design = mm, contrast = contrast,
                  filter = "both", min_cells = 10, method = "edgeR")

out_res <- sprintf("./%s/results/%s_%s_res_edgeR_filtered_withFreq.RDS", REGION, REGION, LABEL)
saveRDS(res.edgeR,      out_res)
saveRDS(pb,             sprintf("./%s/data/pseudo_bulk_%s.rds", REGION, LABEL))
saveRDS(pb.norm_mean,   sprintf("./%s/data/pseudo_norm_mean_bulk_%s.rds", REGION, LABEL))
saveRDS(mm,             sprintf("./%s/data/mm_%s.rds", REGION, LABEL))

## record exactly what this run was, so the comparison script needs no guessing
writeLines(c(
  paste("region       :", REGION),
  paste("label        :", LABEL),
  paste("group        :", opt$group),
  paste("vl_threshold :", opt$`vl-threshold`),
  paste("excluded     :", if (length(EXCLUDE)) paste(EXCLUDE, collapse = ",") else "(none)"),
  paste("n_samples    :", nrow(meta)),
  paste("groups       :", paste(names(table(meta$group)), table(meta$group),
                                sep = "=", collapse = " ")),
  paste("input        :", IN)
), sprintf("./%s/results/%s_%s_runinfo.txt", REGION, REGION, LABEL))

message("[sens] wrote ", out_res)
