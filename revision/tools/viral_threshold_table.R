# ==============================================================================
# HIV-aligned reads per cell type, and how the counts change with the detection
# threshold (1, 2, 3, 5, 10 unique reads per barcode), under both the published
# and the CellBender-corrected annotation.
#
# Method
#   Published threshold: the published viral-mapping analysis keeps barcodes with
#   `V2 > 2`, i.e. at least three unique viral counts, for both RNA and DNA (the
#   Venus config sets VIRUS_THRESHOLD=2 independently). ">= 3" is therefore the
#   published rule; the other rows show what a different choice would give.
#   Cell names follow the combined object: <barcode>-1-<donor>-<REGION>.
#   Both annotations are reported. The corrected run keeps ~6% fewer nuclei, so
#   a barcode can be a labelled cell in one and absent from the other.
#   Most viral-read barcodes are not called cells. They are counted as "not a
#   called cell" rather than dropped, because the ratio of in-cell to
#   out-of-cell viral reads measures how much of the signal is ambient.
#   Empty count files are data: 127 of 152 RNA files are 0 bytes because the
#   pipeline ran and found no HIV-mapped barcode. They count as zeros, and
#   coverage is reported, so "13 libraries" means libraries with any HIV RNA
#   read, not libraries processed.
#   Protocol variants: a sample can have several count files. The PFC 09232023
#   run holds ho_protocol2 and ho_protocol3 for each RNA sample (same 93
#   barcodes for 1039, different counts) and hiv_atac and hiv_atac_mixed for
#   each DNA sample; 32 of 301 sample x region x assay combinations are
#   affected. Summing them would double-count, so one file per sample x region x
#   assay is chosen by a stated precedence (RNA: ho_protocol > ho_protocol2 >
#   ho_protocol3; DNA: hiv_atac > hiv_atac_mixed, the variants the published
#   script read) and every variant is reported side by side in
#   viral_file_variants.csv. The published script read an ho_protocol file for
#   the PFC 09232023 run that no longer exists (only ho_protocol2 and
#   ho_protocol3 do); which of the two reproduces the published numbers cannot
#   be determined from the files.
#   DNA barcode space: 10x multiome ATAC and GEX barcodes differ, and whether
#   the hiv_atac files were already translated to GEX space is not documented.
#   The match rate against the object is reported; a near-zero rate means they
#   are still in ATAC space and need the 737K-arc-v1 translation.
#
# Inputs   per sample, headerless CSV of <barcode>,<unique count>:
#            RNA  RAW_DIR/<run>/ho_protocol*/<sample>_uniqbarcode_count.txt
#            DNA  RAW_DIR/<run>/hiv_atac*/<sample>_ATAC_MAH/filtered_barcode_viral_unique.txt
#          ORIG_ANALYSIS_DIR/01-Regional_data_processing/combined/data/seurat_combined_v3.RDS
#          OUT_ROOT/01-Regional_data_processing/combined/data/seurat_combined_v3.RDS
# Outputs  OUT_ROOT/07-Viral_Detection/
#            viral_barcode_counts.csv.gz        every barcode, its counts, both labels
#            viral_threshold_by_celltype.csv    cell type x threshold x assay x annotation
#            viral_threshold_summary.csv        totals per threshold, both annotations
#            viral_threshold_table.md           summary table for the manuscript
#            viral_file_variants.csv, viral_file_coverage.csv
#            annotation_{published,corrected}.csv.gz, viral_files.txt  (caches)
#
# Usage
#   revision/docker/run.sh revision/tools/viral_threshold_table.R [--refresh]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(data.table) })
args    <- commandArgs(trailingOnly = TRUE)
REFRESH <- "--refresh" %in% args

RAW <- RAW_DIR
OUT <- file.path(OUT_ROOT, "07-Viral_Detection")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
THRESH <- c(1, 2, 3, 5, 10)          # 3 = the published rule (> 2)
PUBLISHED_T <- 3

OBJ <- list(
  published = file.path(ORIG_ANALYSIS_DIR, "01-Regional_data_processing/combined/data/seurat_combined_v3.RDS"),
  corrected = file.path(OUT_ROOT, "01-Regional_data_processing", "combined", "data",
                        "seurat_combined_v3.RDS"))

## ---- annotation lookup, cached -----------------------------------------------------------
ann_file <- function(k) file.path(OUT, paste0("annotation_", k, ".csv.gz"))
for (k in names(OBJ)) {
  if (!REFRESH && file.exists(ann_file(k))) next
  if (!file.exists(OBJ[[k]])) { message("[viral] MISSING object: ", OBJ[[k]]); next }
  suppressPackageStartupMessages(library(Seurat))
  message("[viral] loading ", k, " object (this is the slow step)")
  o  <- readRDS(OBJ[[k]])
  md <- o@meta.data
  col <- intersect(c("consensus_annotation_v3", "consensus_annotation", "cell_types_level0"),
                   colnames(md))[1]
  if (is.na(col)) stop("no annotation column in ", OBJ[[k]])
  message("[viral] ", k, ": ", ncol(o), " cells; annotation column ", col)
  message("[viral] ", k, ": example cell name ", colnames(o)[1])
  data.table::fwrite(data.table(cell = colnames(o), cell_type = as.character(md[[col]])),
                     ann_file(k))
  rm(o, md); invisible(gc())
}
ANN <- lapply(names(OBJ), function(k)
  if (file.exists(ann_file(k))) data.table::fread(ann_file(k)) else NULL)
names(ANN) <- names(OBJ)

## ---- the paper cohort ----------------------------------------------------------------------
meta <- load_sample_meta()
cohort <- unique(toupper(paste(meta$sample_ids, meta$region, sep = "_")))
message("[viral] paper cohort: ", length(cohort), " sample x region combinations")

## ---- find and parse the viral files ---------------------------------------------------------
## The recursive listing over RAW_DIR is slow (network storage), so it is cached.
FILE_LIST <- file.path(OUT, "viral_files.txt")
if (REFRESH || !file.exists(FILE_LIST)) {
  files <- c(list.files(RAW, pattern = "uniqbarcode_count\\.txt$", recursive = TRUE, full.names = TRUE),
             list.files(RAW, pattern = "^filtered_barcode_viral_unique\\.txt$", recursive = TRUE,
                        full.names = TRUE))
  writeLines(files, FILE_LIST)
} else files <- readLines(FILE_LIST)
message("[viral] ", length(files), " count files under ", RAW)

## RNA prefers the plain protocol folder; DNA prefers the unmixed one. Lower
## number = earlier variant. Precedence is stated, not inferred from the data.
RNA_ORDER <- c("ho_protocol", "ho_protocol2", "ho_protocol3")
DNA_ORDER <- c("hiv_atac", "hiv_atac_mixed")

parse_one <- function(f) {
  assay <- if (grepl("hiv_atac", f)) "DNA" else "RNA"
  ## RNA: <donor>_<cond>_<REGION>[12]_MAH_uniqbarcode_count.txt
  ## DNA: .../hiv_atac/<donor>_<cond>_<REGION>_ATAC_MAH/filtered_barcode_viral_unique.txt
  stem  <- if (assay == "DNA") basename(dirname(f)) else basename(f)
  donor <- sub("_.*$", "", stem)
  reg   <- regmatches(toupper(stem), regexpr("INS|PFC|VST", toupper(stem)))
  if (!length(reg)) return(NULL)
  d <- try(data.table::fread(f, header = FALSE, col.names = c("barcode", "n"),
                             showProgress = FALSE), silent = TRUE)
  if (inherits(d, "try-error") || !nrow(d)) return(NULL)
  rel <- sub(paste0("^", RAW, "/"), "", f)
  d[, `:=`(donor = donor, region = reg, assay = assay, file = f,
           run = strsplit(rel, "/")[[1]][1], variant = strsplit(rel, "/")[[1]][2])]
  d
}
tx <- data.table::rbindlist(lapply(files, parse_one), fill = TRUE)
message("[viral] ", nrow(tx), " barcode rows before cohort filtering")

## ---- coverage: which libraries were PROCESSED, empty files included ----------------------
cov <- data.table::rbindlist(lapply(files, function(f) {
  assay <- if (grepl("hiv_atac", f)) "DNA" else "RNA"
  stem  <- if (assay == "DNA") basename(dirname(f)) else basename(f)
  reg   <- regmatches(toupper(stem), regexpr("INS|PFC|VST", toupper(stem)))
  if (!length(reg)) return(NULL)
  rel <- sub(paste0("^", RAW, "/"), "", f)
  data.table(donor = sub("_.*$", "", stem), region = reg, assay = assay,
             variant = strsplit(rel, "/")[[1]][2], bytes = file.size(f))
}))
cov[, key := toupper(paste(donor, region, sep = "_"))]
cov <- cov[key %in% cohort]

tx[, key := toupper(paste(donor, region, sep = "_"))]
in_cohort <- tx$key %in% cohort
message("[viral] rows in the paper cohort: ", sum(in_cohort), " of ", nrow(tx),
        " (samples outside the cohort: ",
        paste(utils::head(sort(unique(tx$key[!in_cohort])), 8), collapse = ", "), " ...)")
tx <- tx[in_cohort]
## ---- one file per sample x region x assay -------------------------------------------------
variants <- tx[, .(n_barcodes = .N, total_counts = sum(n),
                   barcodes_ge3 = sum(n >= PUBLISHED_T), counts_ge3 = sum(n[n >= PUBLISHED_T])),
               by = .(donor, region, assay, run, variant)]
data.table::setorder(variants, donor, region, assay, variant)
data.table::fwrite(variants, file.path(OUT, "viral_file_variants.csv"))
dup <- variants[, .N, by = .(donor, region, assay)][N > 1]
if (nrow(dup)) {
  message("[viral] ", nrow(dup), " sample x region x assay combinations have MORE THAN ONE file; ",
          "they are NOT summed. Side-by-side values: viral_file_variants.csv")
  print(utils::head(as.data.frame(
    variants[dup, on = .(donor, region, assay)][, .(donor, region, assay, variant,
                                                    n_barcodes, total_counts, counts_ge3)]), 12),
        row.names = FALSE)
}
rank_of <- function(assay, variant) {
  o <- ifelse(assay == "RNA", match(variant, RNA_ORDER), match(variant, DNA_ORDER))
  ifelse(is.na(o), 99L, o)
}
tx[, vrank := rank_of(assay, variant)]
tx[, keep := vrank == min(vrank), by = .(donor, region, assay)]
drop_n <- sum(!tx$keep)
if (drop_n) message("[viral] dropped ", drop_n, " rows from non-preferred variants")
tx <- tx[keep == TRUE]
## a sample sequenced twice under the SAME variant would still be duplicated
resq <- tx[, uniqueN(run), by = .(donor, region, assay)][V1 > 1]
if (nrow(resq)) message("[viral] NOTE ", nrow(resq), " combinations appear in more than one ",
                        "sequencing run under the same variant: ",
                        paste(resq$donor, resq$region, resq$assay, collapse = "; "))
message("[viral] variants used: ",
        paste(sprintf("%s=%s", tx[, .N, by = variant]$variant, tx[, .N, by = variant]$N),
              collapse = ", "))

tx[, cell := paste0(barcode, "-1-", donor, "-", region)]
message("[viral] samples represented: ", length(unique(tx$key)),
        "  (RNA ", length(unique(tx$key[tx$assay == "RNA"])),
        ", DNA ", length(unique(tx$key[tx$assay == "DNA"])), ")")

## ---- join both annotations --------------------------------------------------------------
for (k in names(ANN)) {
  if (is.null(ANN[[k]])) { tx[[paste0("ct_", k)]] <- NA_character_; next }
  tx[[paste0("ct_", k)]] <- ANN[[k]]$cell_type[match(tx$cell, ANN[[k]]$cell)]
}
for (a in c("RNA", "DNA")) for (k in names(ANN)) {
  col <- paste0("ct_", k)
  m <- mean(!is.na(tx[[col]][tx$assay == a]))
  message(sprintf("[viral] %s barcodes matching a %s cell: %.1f%%", a, k, 100 * m))
  if (a == "DNA" && m < 0.01)
    message("[viral] WARNING DNA barcodes essentially never match - they are probably still ",
            "in ATAC barcode space and need the 737K-arc-v1 translation before use")
}
data.table::fwrite(tx, file.path(OUT, "viral_barcode_counts.csv.gz"))

## ---- the threshold tables ------------------------------------------------------------------
lab <- function(x) ifelse(is.na(x), "not a called cell", x)
per_ct <- data.table::rbindlist(lapply(THRESH, function(t)
  data.table::rbindlist(lapply(c("RNA", "DNA"), function(a)
    data.table::rbindlist(lapply(names(ANN), function(k) {
      d <- tx[assay == a & n >= t]
      if (!nrow(d)) return(NULL)
      d[, .(threshold = t, assay = a, annotation = k,
            n_barcodes = .N, total_viral_counts = sum(n),
            n_donors = uniqueN(donor), n_samples = uniqueN(key)),
        by = .(cell_type = lab(get(paste0("ct_", k))))]
    }))))))
data.table::setorder(per_ct, annotation, assay, threshold, -n_barcodes)
data.table::fwrite(per_ct, file.path(OUT, "viral_threshold_by_celltype.csv"))

summ <- data.table::rbindlist(lapply(THRESH, function(t) {
  d <- tx[n >= t]
  rna <- d[assay == "RNA"]; dna <- d[assay == "DNA"]
  both <- length(intersect(rna$cell, dna$cell))
  data.table(
    threshold = t,
    published_rule = ifelse(t == PUBLISHED_T, "yes", ""),
    rna_barcodes = nrow(rna), rna_counts = sum(rna$n),
    rna_in_published_cells = sum(!is.na(rna$ct_published)),
    rna_in_corrected_cells = sum(!is.na(rna$ct_corrected)),
    dna_barcodes = nrow(dna), dna_counts = sum(dna$n),
    dna_in_published_cells = sum(!is.na(dna$ct_published)),
    dna_in_corrected_cells = sum(!is.na(dna$ct_corrected)),
    rna_and_dna_same_cell = both,
    donors_with_any = uniqueN(d$donor),
    rna_libraries_with_any = uniqueN(rna$key), dna_libraries_with_any = uniqueN(dna$key))
}))
data.table::fwrite(summ, file.path(OUT, "viral_threshold_summary.csv"))
## ---- coverage report -----------------------------------------------------------------------
cov[, vrank := rank_of(assay, variant)]
cov <- cov[cov[, .I[vrank == min(vrank)], by = .(donor, region, assay)]$V1]
cov_summ <- cov[, .(libraries_processed = .N,
                    empty_no_viral_read = sum(bytes == 0),
                    with_any_viral_barcode = sum(bytes > 0)), by = assay]
utils::write.csv(cov, file.path(OUT, "viral_file_coverage.csv"), row.names = FALSE)
message("\n[viral] ===== coverage of the ", length(cohort), "-library cohort =====")
print(as.data.frame(cov_summ), row.names = FALSE)
miss <- setdiff(cohort, unique(cov$key))
if (length(miss)) message("[viral] libraries with NO count file at all: ", length(miss), " (",
                          paste(utils::head(miss, 6), collapse = ", "), ")")

message("\n[viral] ===== summary by threshold =====")
print(as.data.frame(summ), row.names = FALSE)

## ---- a small markdown table for the manuscript ----------------------------------------------
md <- c("# HIV read detection by threshold",
        "",
        paste0("Unique HIV-mapped counts per nucleus barcode, from the RNA (`ho_protocol`) and ",
               "DNA (`hiv_atac`) pipelines, restricted to the ", length(cohort),
               " sample-by-region libraries of the published cohort. The published analysis used ",
               "a threshold of >2, i.e. at least 3 unique counts."),
        "",
        "| threshold (unique counts) | RNA barcodes | RNA in a called cell (published / corrected) | DNA barcodes | DNA in a called cell (published / corrected) | same cell in both |",
        "|---|---|---|---|---|---|")
for (i in seq_len(nrow(summ))) with(summ[i], md <<- c(md, sprintf(
  "| >= %d%s | %d | %d / %d | %d | %d / %d | %d |",
  threshold, ifelse(published_rule == "yes", " **(published)**", ""),
  rna_barcodes, rna_in_published_cells, rna_in_corrected_cells,
  dna_barcodes, dna_in_published_cells, dna_in_corrected_cells, rna_and_dna_same_cell)))
md <- c(md, "", "Per cell type: `viral_threshold_by_celltype.csv`.")
writeLines(md, file.path(OUT, "viral_threshold_table.md"))
message("[viral] wrote ", OUT)
