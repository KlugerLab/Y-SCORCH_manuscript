# Attributes the cells the corrected combined object loses relative to the published one, by sample, region and cell type.
#
# Method: the corrected seurat_combined.RDS holds 661,409 cells, the published
#   one 673,993 (12,584 fewer, -1.9%), about the size of one sample, so this
#   establishes whether the loss is one sample or a thin slice off all of them.
#   Samples are matched on donor|region because donor directory names differ
#   between the pipelines. The loss arises upstream, in the per-cell QC of the
#   regional processing: the >10% mitochondrial cut, the <0.5 label-transfer
#   score cut and the biccn/ma coarse-consistency cut. Ambient correction changes
#   percent.mt and both prediction scores, so cells move across all three
#   thresholds. The script does not attribute the loss to a specific threshold;
#   a loss spread across samples is the expected consequence of rerunning QC on
#   corrected counts, while one dominant sample would need follow-up.
#
# Inputs:  ORIG_ANALYSIS_DIR/01-Regional_data_processing/combined/data/seurat_combined.RDS (published)
#          OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined.RDS
# Outputs: OUT_ROOT/01-Regional_data_processing/combined/data/cmp_published_by_{sample,region,annotation}.csv,
#          cmp_missing_cells.csv (only if barcodes correspond). Read-only otherwise.
#
# Usage: revision/docker/run.sh revision/tools/compare_combined_to_published.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

rerun_setup("01-Regional_data_processing", subdirs = "combined/data")

PUB <- file.path(ORIG_ANALYSIS_DIR,
                 "01-Regional_data_processing/combined/data/seurat_combined.RDS")
OURS <- "./combined/data/seurat_combined.RDS"
stopifnot(file.exists(PUB), file.exists(OURS))

message("[cmp] loading published ", PUB, " ...")
pub <- readRDS(PUB)
message("[cmp] published: ", ncol(pub), " cells")

message("[cmp] loading ours ...")
our <- readRDS(OURS)
message("[cmp] ours     : ", ncol(our), " cells")
message("[cmp] difference: ", ncol(pub) - ncol(our),
        sprintf("  (%.2f%%)", 100 * (ncol(pub) - ncol(our)) / ncol(pub)))

pm <- pub@meta.data
om <- our@meta.data

## ---- do barcodes even correspond? --------------------------------------------
## If they do, an exact set difference is possible. If naming changed between
## pipelines, only aggregate counts are comparable - say so rather than
## silently producing a meaningless "missing" list.
ov <- length(intersect(colnames(pub), colnames(our)))
message("\n[cmp] barcode overlap: ", ov, " of ", ncol(our), " ours (",
        sprintf("%.1f%%", 100 * ov / ncol(our)), ")")
COMPARABLE <- ov > 0.5 * ncol(our)
if (!COMPARABLE) {
  message("[cmp] barcodes DO NOT correspond between the two objects - ",
          "reporting aggregate counts only, not a per-cell missing list.")
  message("[cmp] published example: ", paste(utils::head(colnames(pub), 2), collapse = " | "))
  message("[cmp] ours example     : ", paste(utils::head(colnames(our), 2), collapse = " | "))
}

## ---- pick the columns that exist in both -------------------------------------
pick <- function(cands, m) { h <- intersect(cands, colnames(m)); if (length(h)) h[1] else NA }
SAMP_P <- pick(c("sample_name","sample_ids","sample_id"), pm)
SAMP_O <- pick(c("sample_name","sample_ids","sample_id"), om)
REG_P  <- pick(c("region","brain_region"), pm)
REG_O  <- pick(c("region","brain_region"), om)
ANN_P  <- pick(c("consensus_annotation_v2","consensus_annotation"), pm)
ANN_O  <- pick(c("consensus_annotation_v2","consensus_annotation"), om)
message("[cmp] columns - sample: ", SAMP_P, "/", SAMP_O,
        "  region: ", REG_P, "/", REG_O, "  annotation: ", ANN_P, "/", ANN_O)

cmp_counts <- function(vp, vo, label) {
  tp <- table(as.character(vp)); to <- table(as.character(vo))
  keys <- union(names(tp), names(to))
  d <- data.frame(
    key       = keys,
    published = as.integer(tp[keys]),
    ours      = as.integer(to[keys]),
    stringsAsFactors = FALSE
  )
  d$published[is.na(d$published)] <- 0L
  d$ours[is.na(d$ours)] <- 0L
  d$diff <- d$ours - d$published
  d$pct  <- round(100 * d$diff / pmax(d$published, 1), 2)
  d <- d[order(d$diff), ]
  names(d)[1] <- label
  d
}

## ---- 1. by sample - the question that matters --------------------------------
## Donor directories were renamed between pipelines: 1039_HIVOUD_INS became
## 1039_HIV_INS. Comparing sample_name literally reports seven renamed INS
## donors as "entirely absent" and their new names as "new", which is a pure
## naming artifact and reads as catastrophic. Key on donor|region instead.
norm_sample <- function(x) {
  p <- strsplit(as.character(x), "_", fixed = TRUE)
  vapply(p, function(v) paste0(v[1], "|", v[length(v)]), character(1))
}
if (!is.na(SAMP_P) && !is.na(SAMP_O)) {
  by_s <- cmp_counts(norm_sample(pm[[SAMP_P]]), norm_sample(om[[SAMP_O]]), "sample")
  write.csv(by_s, "./combined/data/cmp_published_by_sample.csv", row.names = FALSE)
  message("\n[cmp] === by sample: 15 largest losses ===")
  print(utils::head(by_s, 15), row.names = FALSE)

  gone <- by_s$sample[by_s$ours == 0 & by_s$published > 0]
  if (length(gone)) {
    message("\n[cmp] *** SAMPLES ENTIRELY ABSENT FROM OURS: ",
            paste(gone, collapse = ", "), " ***")
  } else {
    message("\n[cmp] no sample is entirely absent - every published sample is present")
  }
  new_s <- by_s$sample[by_s$published == 0 & by_s$ours > 0]
  if (length(new_s))
    message("[cmp] samples present only in ours: ", paste(new_s, collapse = ", "))

  ## is the loss concentrated or spread?
  tot_loss <- sum(pmin(by_s$diff, 0))
  worst <- utils::head(by_s[by_s$diff < 0, ], 1)
  if (nrow(worst))
    message(sprintf("[cmp] largest single-sample loss: %s  %d cells = %.1f%% of the total loss",
                    worst$sample[1], worst$diff[1], 100 * worst$diff[1] / tot_loss))
  n_losing <- sum(by_s$diff < 0)
  message("[cmp] samples losing cells: ", n_losing, " of ", nrow(by_s))
}

## ---- 2. by region -------------------------------------------------------------
if (!is.na(REG_P) && !is.na(REG_O)) {
  by_r <- cmp_counts(pm[[REG_P]], om[[REG_O]], "region")
  message("\n[cmp] === by region ===")
  print(by_r, row.names = FALSE)
  write.csv(by_r, "./combined/data/cmp_published_by_region.csv", row.names = FALSE)
}

## ---- 3. by cell type ----------------------------------------------------------
if (!is.na(ANN_P) && !is.na(ANN_O)) {
  by_a <- cmp_counts(pm[[ANN_P]], om[[ANN_O]], "annotation")
  message("\n[cmp] === by annotation ===")
  print(by_a, row.names = FALSE)
  write.csv(by_a, "./combined/data/cmp_published_by_annotation.csv", row.names = FALSE)
}

## ---- 4. exact missing cells, if barcodes correspond ---------------------------
if (COMPARABLE) {
  missing <- setdiff(colnames(pub), colnames(our))
  gained  <- setdiff(colnames(our), colnames(pub))
  message("\n[cmp] cells in published but not ours: ", length(missing))
  message("[cmp] cells in ours but not published: ", length(gained))

  if (length(missing)) {
    mm <- pm[missing, , drop = FALSE]
    message("[cmp] missing cells by annotation:")
    print(sort(table(as.character(mm[[ANN_P]])), decreasing = TRUE))
    if (!is.na(REG_P)) {
      message("[cmp] missing cells by region:")
      print(table(as.character(mm[[REG_P]])))
    }
    ## percent.mt of the lost cells, if recorded - speaks to the MT threshold
    if ("percent.mt" %in% colnames(mm)) {
      message("[cmp] percent.mt of missing cells: ",
              sprintf("median %.2f, mean %.2f, >10%%: %d (%.1f%%)",
                      median(mm$percent.mt, na.rm = TRUE),
                      mean(mm$percent.mt, na.rm = TRUE),
                      sum(mm$percent.mt > 10, na.rm = TRUE),
                      100 * mean(mm$percent.mt > 10, na.rm = TRUE)))
      message("[cmp]   for comparison, all published cells: ",
              sprintf("median %.2f, >10%%: %.1f%%",
                      median(pm$percent.mt, na.rm = TRUE),
                      100 * mean(pm$percent.mt > 10, na.rm = TRUE)))
    }
    out <- data.frame(cell = missing,
                      sample = if (!is.na(SAMP_P)) as.character(mm[[SAMP_P]]) else NA,
                      region = if (!is.na(REG_P)) as.character(mm[[REG_P]]) else NA,
                      annotation = if (!is.na(ANN_P)) as.character(mm[[ANN_P]]) else NA,
                      stringsAsFactors = FALSE)
    write.csv(out, "./combined/data/cmp_missing_cells.csv", row.names = FALSE)
    message("[cmp] full missing-cell list -> combined/data/cmp_missing_cells.csv")
  }
}

message("\n[cmp] DONE")
