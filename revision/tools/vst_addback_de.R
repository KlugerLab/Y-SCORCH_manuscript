# ==============================================================================
# VST exclusion sensitivity: regional VST edgeR DE with the five excluded
# libraries added back (reviewer request: how much the VST findings depend on
# the exclusions).
#
# Method
#   04-01 drops the five excluded libraries (1039, 1040, 1068, HCCRS, HCTXJ)
#   before annotation, so their nuclei carry no cell type label and the regional
#   DE cannot simply be re-run with them kept; re-running 04-01 without the
#   exclusion would re-integrate and re-cluster all of VST (and drop the patched
#   MSN subtypes). Instead:
#   1 take the five libraries' nuclei from the stage-03 object (all 33
#     libraries), minus the putative excitatory nuclei, the same cell-level rule
#     04-01 applies to the libraries it keeps;
#   2 label them by Seurat label transfer from the current VST annotation
#     (consensus_annotation_v3 of integrated_seurat_obj_v3.RDS, MSN subtypes
#     included); reference = up to 5,000 nuclei per label, LogNormalize,
#     2,000 HVGs, 30 PCs. Nuclei with prediction score < 0.5 are left out and
#     counted;
#   3 run the regional VST DE recipe of 03-DE_analysis/03-VST_edgeR.R unchanged
#     (uncorrected counts, prepSCE -> aggregateData(sum) ->
#     ~ 0 + HIV + Age + Sex + Race + PMI, Age and PMI scaled over the donors in
#     the run -> pbDS(edgeR, filter "both", min_cells 10)) twice:
#       rebuild28   the 28 primary libraries only; must reproduce the primary
#                   VST result, and the script stops if logFC or p_adj.glb differ
#       addback33   the same 28 plus the five labelled libraries
#   The difference between the two is therefore the added libraries only.
#   --race collapsed merges the Native Alaskan/American Indian level into
#   "Other" (see the note below).
#
# Inputs   OUT_ROOT/01-Regional_data_processing/VST/data/
#            03-cell_integration_vis/integrated_seurat_obj.RDS
#            04-VST-tsne_marker_plotting/integrated_seurat_obj_v3.RDS
#          OUT_ROOT/03-DE_analysis/VST/results/VST_HIV-CTR_res_edgeR_filtered_withFreq.RDS
# Outputs  OUT_ROOT/03-DE_analysis/VST/results/
#            VST_addback33[_raceCollapsed]_res_edgeR.RDS, VST_rebuild28[_raceCollapsed]_res_edgeR.RDS
#          OUT_ROOT/03-DE_analysis/comparisons/vst_addback/
#            addback_labels.csv            every added nucleus: sample, label, score
#            addback_label_summary.csv     per library x label
#            rebuild28_vs_primary_check.txt
#
# Usage
#   revision/docker/run.sh revision/tools/vst_addback_de.R [--race published|collapsed]
#   then 03-DE_analysis/11-compare_sensitivity.R --a <primary> --b VST/results/VST_addback33_res_edgeR.RDS
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(muscat); library(limma); library(dplyr)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
set.seed(1)
t0 <- Sys.time(); el <- function() format(round(difftime(Sys.time(), t0, units = "mins"), 1))
VD  <- file.path(OUT_ROOT, "01-Regional_data_processing", "VST", "data")
RES <- file.path(OUT_ROOT, "03-DE_analysis", "VST", "results")
CMP <- file.path(OUT_ROOT, "03-DE_analysis", "comparisons", "vst_addback")
dir.create(CMP, recursive = TRUE, showWarnings = FALSE)
## --race collapsed: one added donor carries a race level no primary donor has
##   (Native Alaskan/American Indian), a one-donor design column. Wherever that donor
##   has < 10 nuclei the design loses rank and pbDS silently drops the cell type (with
##   the published coding the 33-library run tests 9 cell types, not 14: D1/D2-Shell/OT,
##   Immune, Macrophages and VLMC are lost). The collapsed coding is run for both 28
##   and 33 libraries so the comparison stays like-for-like.
args <- commandArgs(trailingOnly = TRUE)
RACE <- if ("--race" %in% args) args[match("--race", args) + 1] else "published"
RACE <- match.arg(RACE, c("published", "collapsed"))
RSFX <- if (RACE == "collapsed") "_raceCollapsed" else ""
LAB <- "consensus_annotation_v3"
ADD <- c("1039", "1040", "1068", "HCCRS", "HCTXJ")
MIN_SCORE <- 0.5

## ---- 1. the excluded libraries' non-excitatory nuclei ------------------------------
s3 <- readRDS(file.path(VD, "03-cell_integration_vis", "integrated_seurat_obj.RDS"))
md3 <- s3@meta.data
if (!"putative_ext" %in% colnames(md3)) stop("stage-03 object has no putative_ext column")
tab_all <- table(as.character(md3$sample_ids)[as.character(md3$sample_ids) %in% ADD])
q_cells <- rownames(md3)[as.character(md3$sample_ids) %in% ADD & !md3$putative_ext]
message(sprintf("[addback] stage-03: %d nuclei in the 5 libraries, %d after removing putative excitatory nuclei (%s min)",
                sum(tab_all), length(q_cells), el()))
q <- s3[, q_cells]; rm(s3); invisible(gc())
q[["RNA"]] <- JoinLayers(q[["RNA"]])

## ---- 2. label transfer from the current VST annotation ------------------------------
ref <- readRDS(file.path(VD, "04-VST-tsne_marker_plotting", "integrated_seurat_obj_v3.RDS"))
if (any(ADD %in% as.character(ref$sample_ids))) stop("the reference already holds an excluded library")
message(sprintf("[addback] reference: %d nuclei, %d libraries (%s min)", ncol(ref),
                length(unique(ref$sample_ids)), el()))
rl <- as.character(ref@meta.data[[LAB]])
pick <- unlist(lapply(split(colnames(ref), rl), function(v) if (length(v) > 5000) sample(v, 5000) else v))
rs <- ref[, pick]
rs[["RNA"]] <- JoinLayers(rs[["RNA"]])
rs <- NormalizeData(rs, verbose = FALSE) |> FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE)
q  <- NormalizeData(q, verbose = FALSE)
anc <- FindTransferAnchors(reference = rs, query = q, normalization.method = "LogNormalize",
                           reference.reduction = "pca", dims = 1:30, verbose = FALSE)
pred <- TransferData(anchorset = anc, refdata = as.character(rs@meta.data[[LAB]]), dims = 1:30, verbose = FALSE)
lab_q <- data.frame(cell = colnames(q), sample_ids = as.character(q$sample_ids),
                    condition_new = as.character(q$condition_new),
                    label = pred$predicted.id, score = pred$prediction.score.max, stringsAsFactors = FALSE)
lab_q$kept <- lab_q$score >= MIN_SCORE
utils::write.csv(lab_q, file.path(CMP, "addback_labels.csv"), row.names = FALSE)
lsum <- lab_q %>% group_by(sample_ids, condition_new, label) %>%
  summarise(nuclei = dplyr::n(), kept = sum(kept), median_score = median(score), .groups = "drop")
utils::write.csv(lsum, file.path(CMP, "addback_label_summary.csv"), row.names = FALSE)
message(sprintf("[addback] labelled %d nuclei; %d (%.1f%%) score >= %.1f (%s min)", nrow(lab_q), sum(lab_q$kept),
                100 * mean(lab_q$kept), MIN_SCORE, el()))
print(as.data.frame(lab_q %>% dplyr::filter(kept) %>% dplyr::count(sample_ids, label) %>%
                      tidyr::pivot_wider(names_from = sample_ids, values_from = n, values_fill = 0)), row.names = FALSE)
rm(rs, anc, pred); invisible(gc())

## ---- 3. counts: uncorrected, as the primary DE ---------------------------------------
c28 <- uncorrected_counts(ref, "VST")
qk  <- q[, lab_q$cell[lab_q$kept]]
cadd <- uncorrected_counts(qk, "VST")
if (!identical(rownames(c28), rownames(cadd))) {
  g <- intersect(rownames(c28), rownames(cadd))
  message("[addback] gene sets differ (", nrow(c28), " vs ", nrow(cadd), "); using the ", length(g), " shared")
  c28 <- c28[g, ]; cadd <- cadd[g, ]
}
cd28 <- ref@meta.data[, c(LAB, "condition_new", "sample_ids")]
cdad <- data.frame(lab_q$label[lab_q$kept], qk$condition_new, qk$sample_ids,
                   row.names = colnames(qk), stringsAsFactors = FALSE)
names(cdad) <- names(cd28)
stopifnot(all(c(cd28$condition_new, cdad$condition_new) %in% c("HIV", "CTR")))
rm(ref, q, qk); invisible(gc())
message("[addback] counts assembled (", el(), " min)")

## ---- the regional VST recipe (03-VST_edgeR.R), unchanged ------------------------------
run_de <- function(counts, cd, tag) {
  sce <- SingleCellExperiment(list(counts = counts), colData = cd)
  sce <- prepSCE(sce, kid = LAB, gid = "condition_new", sid = "sample_ids", drop = FALSE)
  pb  <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
  ei  <- metadata(sce)$experiment_info
  ei$group_id <- factor(ei$group_id, levels = c("HIV", "CTR"))
  sm <- utils::read.csv(SAMPLE_META) %>% dplyr::filter(!is.na(cellranger_path)) %>%
    dplyr::filter((oud_jay_v2 == FALSE) & (region == "VST")) %>%
    dplyr::select(sample_ids, age_value_v2, sex_v2, race_v2, postmortem_interval_v2, hiv_jay_v2)
  sm$sample_id <- sm$sample_ids
  ei <- ei %>% left_join(sm, by = "sample_id")
  ei$Age  <- scale(as.numeric(ei$age_value_v2))[, 1]
  ei$Sex  <- factor(ei$sex_v2, levels = c("male", "female"))
  ei$Race <- make.names(ei$race_v2)
  ## the pipeline's SCORCH_RACE_CODING = "collapsed" variant (03-VST_edgeR.R):
  ## Native Alaskan/American Indian levels merged into "Other"
  if (RACE == "collapsed") ei$Race[grepl("Native", ei$Race)] <- "Other"
  ei$PMI  <- scale(as.numeric(ei$postmortem_interval_v2))[, 1]
  ei$HIV  <- ei$group_id
  if (anyNA(ei$Age) || anyNA(ei$PMI) || anyNA(ei$Sex)) stop(tag, ": covariate NA after the sample-sheet join")
  mm <- model.matrix(~ 0 + ei$HIV + ei$Age + ei$Sex + ei$Race + ei$PMI)
  colnames(mm) <- sub("^ei\\$HIV", "", colnames(mm))
  colnames(mm) <- sub("^ei\\$Race", "Race", colnames(mm))
  colnames(mm) <- sub("^ei\\$Sex.*", "Sex", colnames(mm))
  colnames(mm) <- sub("^ei\\$Age$", "Age", colnames(mm))
  colnames(mm) <- sub("^ei\\$PMI$", "PMI", colnames(mm))
  colnames(mm) <- make.names(colnames(mm))
  rownames(mm) <- as.character(ei$sample_id)
  message(sprintf("[addback] %s: %d libraries (HIV %d, CTR %d), %d cell types in the pseudobulk", tag, nrow(ei),
                  sum(ei$HIV == "HIV"), sum(ei$HIV == "CTR"), length(SummarizedExperiment::assayNames(pb))))
  contrast <- makeContrasts("HIV-CTR", levels = mm)
  res <- pbDS(pb, design = mm, contrast = contrast, filter = "both", min_cells = 10, method = "edgeR")
  message(sprintf("[addback] %s: %d cell types tested (%s min)", tag, length(res$table[[1]]), el()))
  resDS(sce, res, frq = FALSE)
}

r28 <- run_de(c28, cd28, "rebuild28")
saveRDS(r28, file.path(RES, sprintf("VST_rebuild28%s_res_edgeR.RDS", RSFX)))

## ---- check: the rebuild must reproduce the primary result -----------------------------
## Under the collapsed coding the primary (published-coding) result is not the right
## reference, and the pipeline's own collapsed-coding result lacks the MSN subtypes,
## so the check is skipped and the 28-library run is compared with that result
## separately.
if (RACE == "collapsed") {
  r33 <- run_de(cbind(c28, cadd), rbind(cd28, cdad), "addback33")
  saveRDS(r33, file.path(RES, sprintf("VST_addback33%s_res_edgeR.RDS", RSFX)))
  message("[addback] done (", el(), " min)")
  quit(save = "no")
}
prim <- readRDS(file.path(RES, "VST_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))
j <- dplyr::inner_join(as.data.frame(prim)[, c("cluster_id", "gene", "logFC", "p_adj.glb")],
                       as.data.frame(r28)[, c("cluster_id", "gene", "logFC", "p_adj.glb")],
                       by = c("cluster_id", "gene"), suffix = c("_prim", "_rebuild"))
chk <- c(sprintf("primary rows %d, rebuild rows %d, joined %d", nrow(prim), nrow(r28), nrow(j)),
         sprintf("max |logFC diff| %.3g, max |p_adj.glb diff| %.3g",
                 max(abs(j$logFC_prim - j$logFC_rebuild)), max(abs(j$p_adj.glb_prim - j$p_adj.glb_rebuild))))
writeLines(chk, file.path(CMP, "rebuild28_vs_primary_check.txt")); message("[addback] ", paste(chk, collapse = " | "))
if (nrow(j) != nrow(prim) || nrow(j) != nrow(r28) || max(abs(j$logFC_prim - j$logFC_rebuild)) > 1e-6 ||
    max(abs(j$p_adj.glb_prim - j$p_adj.glb_rebuild)) > 1e-6)
  stop("the 28-library rebuild does NOT reproduce the primary VST result - the add-back comparison would not be clean")

r33 <- run_de(cbind(c28, cadd), rbind(cd28, cdad), "addback33")
saveRDS(r33, file.path(RES, "VST_addback33_res_edgeR.RDS"))
message("[addback] done (", el(), " min)")
