# ==============================================================================
# Supplementary table of nuclei with HIV-aligned reads at >= 1, >= 2 (the Wei et
# al. criterion), >= 3 (the original criterion) and >= 5 reads per nucleus.
#
# Method
#   One row per nucleus, not per cell type: donor, condition, region, assay, the
#   HIV count, and both annotations side by side:
#     annotation_published          the published combined object (654,918 nuclei)
#     annotation_after_cellbender   the ambient-corrected reprocessing (641,923 nuclei)
#   The HIV counts are identical under both: viral mapping ran on the raw FASTQs
#   and was not repeated. Only cell calling and labelling differ, so a barcode can
#   be a labelled nucleus in one column and absent from the other.
#   Each row carries its exact count, so any threshold can be re-derived; the
#   highest_threshold_met column is a convenience for the caption.
#   "Nuclei" and "barcodes" are reported separately: the supplementary tables
#   hold barcodes that are a nucleus in at least one annotation, while
#   viral_all_barcodes.csv keeps every barcode with HIV reads, including those
#   that are not nuclei in either run. At the original threshold 63% of the HIV
#   RNA signal sits in these non-nuclei, which bears on the ambient-background
#   interpretation and is therefore kept.
#
# Inputs   OUT_ROOT/07-Viral_Detection/viral_barcode_counts.csv.gz
#          (per-barcode table cached by revision/tools/viral_threshold_table.R)
# Outputs  OUT_ROOT/07-Viral_Detection/
#            viral_nuclei_supplementary_{RNA,ATAC}.csv   nuclei in either annotation
#            viral_nuclei_supplementary.md
#            viral_nuclei_cellbender_{RNA,ATAC}.csv      CellBender nuclei only
#            viral_nuclei_cellbender_by_threshold.csv, viral_nuclei_cellbender_by_celltype.csv
#            viral_all_barcodes.csv                      every barcode, with provenance
#            viral_nuclei_counts_by_threshold.csv
#
# Usage
#   revision/docker/run.sh revision/tools/viral_nuclei_table.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(data.table) })

OUT <- file.path(OUT_ROOT, "07-Viral_Detection")
THRESH <- c(1, 2, 3, 5)
THRESH_NOTE <- c("1" = "", "2" = "Wei et al. criterion", "3" = "original criterion", "5" = "")

## na.strings matters: fwrite wrote unmatched annotations as empty fields, and
## without this fread returns "" rather than NA, so !is.na() is TRUE for every
## barcode and every barcode counts as a nucleus.
tx <- data.table::fread(file.path(OUT, "viral_barcode_counts.csv.gz"),
                        na.strings = c("NA", ""))
message("[nuc] ", nrow(tx), " barcodes with HIV-aligned reads")

## ---- condition, from the sample sheet -----------------------------------------------------
meta <- load_sample_meta()
cond <- unique(meta[, c("sample_ids", "hiv_jay_v2")])
## Donor IDs are inconsistently capitalised in the raw filenames (HCtME, HCcPL,
## HCtNZ against HCTME, HCCPL, HCTNZ in the sheet), so the join is case-insensitive
## and the donor is rewritten to the sheet's spelling for the published table.
idx <- match(toupper(tx$donor), toupper(cond$sample_ids))
tx[, condition := ifelse(toupper(as.character(cond$hiv_jay_v2[idx])) == "TRUE", "PWH", "PWoH")]
tx[, donor := ifelse(is.na(idx), donor, cond$sample_ids[idx])]
if (anyNA(tx$condition)) stop("no condition for donors: ",
                              paste(unique(tx$donor[is.na(tx$condition)]), collapse = ", "))

setnames(tx, c("ct_published", "ct_corrected"),
         c("annotation_published", "annotation_after_cellbender"))
tx[, `:=`(hiv_counts = n,
          nucleus_published = !is.na(annotation_published),
          nucleus_after_cellbender = !is.na(annotation_after_cellbender),
          highest_threshold_met = fifelse(n >= 5, ">=5",
                                  fifelse(n >= 3, ">=3",
                                  fifelse(n >= 2, ">=2", ">=1"))))]

## The supplementary tables are split into one file per assay (RNA, ATAC), so
## `assay` is dropped from them. `variant` names the pipeline output each row
## came from; it is kept in viral_all_barcodes.csv, the provenance file, and
## dropped from the supplementary tables.
COLS <- c("barcode", "donor", "condition", "region", "assay", "hiv_counts",
          "annotation_published", "annotation_after_cellbender",
          "highest_threshold_met", "variant")
SPLIT_COLS <- setdiff(COLS, c("assay", "variant"))
data.table::setorder(tx, -hiv_counts, donor, region)
data.table::fwrite(tx[, ..COLS], file.path(OUT, "viral_all_barcodes.csv"))

nuc <- tx[nucleus_published | nucleus_after_cellbender]
## "not retained" (rather than "not called as a nucleus"): the barcode is absent
## from that object's cell set. Every cell in both objects carries a label
## (0 NAs in 654,918 and 641,923), so a blank can only mean the barcode was not
## retained as a nucleus there. The cause is not established: CellBender may not
## have called it a cell, or a downstream QC filter (mitochondrial fraction,
## doublets, cluster removal) may have dropped it. The wording claims no
## mechanism and reads correctly in both directions; some barcodes are retained
## only after reprocessing (1119 VST).
NOT_RETAINED <- "not retained"
nuc[, annotation_published := fifelse(is.na(annotation_published),
                                      NOT_RETAINED, annotation_published)]
nuc[, annotation_after_cellbender := fifelse(is.na(annotation_after_cellbender),
                                             NOT_RETAINED, annotation_after_cellbender)]
## one file per assay; "DNA" is the provirus/ATAC pipeline, named ATAC in the outputs
for (a in c("RNA", "DNA")) {
  lab <- ifelse(a == "DNA", "ATAC", a)
  d <- nuc[assay == a]
  data.table::fwrite(d[, ..SPLIT_COLS],
                     file.path(OUT, sprintf("viral_nuclei_supplementary_%s.csv", lab)))
  message("[nuc] ", lab, ": ", nrow(d), " nuclei -> viral_nuclei_supplementary_", lab, ".csv")
}
## ---- CellBender-only version ----------------------------------------------------------------
## The same rows, but judged solely by the ambient-corrected reprocessing: a
## barcode appears only if it is a nucleus there, and there is one annotation
## column rather than two. No "not retained" rows, because under this view a
## barcode that the corrected run did not keep is simply not a nucleus and has no
## place in a table of nuclei. The published annotation is not consulted at all.
cb <- tx[nucleus_after_cellbender == TRUE]
cb[, cell_type := annotation_after_cellbender]
CB_COLS <- c("barcode", "donor", "condition", "region", "hiv_counts",
             "cell_type", "highest_threshold_met")
for (a in c("RNA", "DNA")) {
  lab <- ifelse(a == "DNA", "ATAC", a)
  d <- cb[assay == a]
  data.table::fwrite(d[, ..CB_COLS],
                     file.path(OUT, sprintf("viral_nuclei_cellbender_%s.csv", lab)))
  message("[nuc] CellBender-only ", lab, ": ", nrow(d), " nuclei")
}

## counts per threshold, CellBender nuclei only
cb_summ <- data.table::rbindlist(lapply(THRESH, function(t)
  data.table::rbindlist(lapply(c("RNA", "DNA"), function(a) {
    d <- cb[assay == a & n >= t]
    data.table(threshold = paste0(">=", t), note = unname(THRESH_NOTE[as.character(t)]),
               assay = ifelse(a == "DNA", "ATAC", a),
               nuclei = nrow(d), hiv_counts = sum(d$n),
               donors_PWH = uniqueN(d$donor[d$condition == "PWH"]),
               donors_PWoH = uniqueN(d$donor[d$condition == "PWoH"]),
               nuclei_in_PWoH = sum(d$condition == "PWoH"))
  }))))
data.table::fwrite(cb_summ, file.path(OUT, "viral_nuclei_cellbender_by_threshold.csv"))
message("\n[nuc] ===== CellBender nuclei per threshold =====")
print(as.data.frame(cb_summ), row.names = FALSE)

## per cell type, CellBender nuclei only
cb_ct <- data.table::rbindlist(lapply(THRESH, function(t)
  data.table::rbindlist(lapply(c("RNA", "DNA"), function(a) {
    d <- cb[assay == a & n >= t]
    if (!nrow(d)) return(NULL)
    d[, .(threshold = paste0(">=", t), assay = ifelse(a == "DNA", "ATAC", a),
          nuclei = .N, hiv_counts = sum(n), n_donors = uniqueN(donor)), by = cell_type]
  }))))
data.table::setorder(cb_ct, assay, threshold, -nuclei)
data.table::fwrite(cb_ct, file.path(OUT, "viral_nuclei_cellbender_by_celltype.csv"))
message("\n[nuc] CellBender nuclei by cell type at the original threshold (>=3):")
print(as.data.frame(cb_ct[threshold == ">=3"]), row.names = FALSE)

message("[nuc] nuclei with HIV-aligned reads: ", nrow(nuc),
        " (of ", nrow(tx), " barcodes)")
if (nrow(nuc) == nrow(tx))
  stop("every barcode counted as a nucleus - the annotation join is not working ",
       "(empty fields read as \"\" instead of NA?)")

## ---- counts per threshold, for the sentence ------------------------------------------------
summ <- data.table::rbindlist(lapply(THRESH, function(t)
  data.table::rbindlist(lapply(c("RNA", "DNA"), function(a) {
    d <- tx[assay == a & n >= t]
    data.table(threshold = paste0(">=", t), note = unname(THRESH_NOTE[as.character(t)]),
               assay = a,
               barcodes = nrow(d), hiv_counts = sum(d$n),
               nuclei_published = sum(d$nucleus_published),
               nuclei_after_cellbender = sum(d$nucleus_after_cellbender),
               nuclei_either = sum(d$nucleus_published | d$nucleus_after_cellbender),
               donors_PWH = uniqueN(d$donor[d$condition == "PWH"]),
               donors_PWoH = uniqueN(d$donor[d$condition == "PWoH"]),
               ## the specificity check: viral barcodes in control donors
               barcodes_in_PWoH = sum(d$condition == "PWoH"))
  }))))
data.table::fwrite(summ, file.path(OUT, "viral_nuclei_counts_by_threshold.csv"))
message("\n[nuc] ===== nuclei per threshold =====")
print(as.data.frame(summ), row.names = FALSE)

## ---- the supplementary table as markdown ----------------------------------------------------
md <- c("# Supplementary Table X. Nuclei with HIV-aligned reads",
        "",
        paste0("One row per nucleus, split by assay. HIV counts are unique reads (RNA) or unique ",
               "fragments (ATAC) aligned to the HIV genome; they are identical under both ",
               "annotations because viral mapping was performed on the raw sequencing reads and ",
               "not repeated. `annotation_published` is the cell type in the originally published ",
               "object; `annotation_after_cellbender` is the cell type after ambient RNA correction ",
               "and reprocessing. \"not retained\" means the barcode is not part of that run's ",
               "nucleus set: it was either not called as a nucleus or removed by quality ",
               "control during that processing. Counts are unchanged in either case."))
for (a in c("RNA", "DNA")) {
  lab <- ifelse(a == "DNA", "ATAC", a)
  d <- nuc[assay == a]
  md <- c(md, "", paste0("## ", lab), "",
          "| barcode | donor | condition | region | HIV counts | annotation (published) | annotation (after CellBender) |",
          "|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(d))) with(d[i], md <<- c(md, sprintf(
    "| %s | %s | %s | %s | %d | %s | %s |",
    barcode, donor, condition, region, hiv_counts,
    annotation_published, annotation_after_cellbender)))
}
md <- c(md, "",
        "## Counts per threshold", "",
        "| threshold | assay | barcodes | nuclei (published) | nuclei (after CellBender) | PWH donors | PWoH donors |",
        "|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(summ))) with(summ[i], md <<- c(md, sprintf(
  "| %s%s | %s | %d | %d | %d | %d | %d |",
  threshold, ifelse(nzchar(note), paste0(" (", note, ")"), ""), assay,
  barcodes, nuclei_published, nuclei_after_cellbender, donors_PWH, donors_PWoH)))
writeLines(md, file.path(OUT, "viral_nuclei_supplementary.md"))
message("[nuc] wrote ", OUT)
