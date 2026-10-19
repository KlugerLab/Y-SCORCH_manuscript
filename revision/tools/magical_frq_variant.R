# ==============================================================================
# Builds a MAGICAL input variant with the paper's expression filter applied to
# the candidate genes and peaks.
#
# Method: MAGICAL's candidates are |logFC| >= 0.5 & p_adj.glb <= 0.05 with no
#   expression filter (prepare_MAGICAL_input.R, as published, and as the GRN
#   methods state). The paper's DEG definition also requires expression in > 5%
#   of nuclei in the PWH or the PWoH group, so e.g. CP, SOCS3 and CIITA sit in the
#   PFC astrocyte network (FDR 0.0008 / 0.044 / 0.003) but are in no DEG table.
#   This variant applies the same > 5% rule to the candidate genes and, by
#   default, the candidate peaks (the paper filters DA peaks "with the same
#   cutoffs"), keeping MAGICAL's own |logFC| >= 0.5.
#
#   Everything else MAGICAL loads - read counts, cell metadata, peak list, motif
#   prior - does not depend on which features are candidates, so it is symlinked
#   from the original input folder and only the six candidate files are
#   rewritten. Any difference in the circuits therefore comes from the filter
#   (plus MAGICAL's ~10% sampling noise between runs).
#
#   Self-check: the unfiltered lists are regenerated from the same tables first
#   and must equal the original candidate files exactly; otherwise the script
#   stops, because it would be filtering a different table from the one the
#   original run used.
#
# Inputs: OUT_ROOT/03-DE_analysis/<REGION>/results/<REGION>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#         OUT_ROOT/04-ATAC_DE/<REGION>/data_ATAC/<REGION>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#         MAGICAL_INPUT_DIR/<REGION>_new_<CT>/
# Output: MAGICAL_INPUT_DIR/<REGION>_new_<CT>_frq05/  (then: run_magical.R --variant frq05)
#         .../filter_report.csv   what the filter removed
#
# Usage: revision/docker/run.sh revision/tools/magical_frq_variant.R --region PFC --celltype Astro
#        add --peaks no to filter genes only (variant frq05genes)
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(dplyr); library(data.table) })
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
REGION  <- getarg("--region", "PFC"); CT <- getarg("--celltype", "Astro")
FRQ     <- as.numeric(getarg("--min-frq", "0.05"))
## The 5% rule removes few candidate genes but most candidate peaks (snATAC is
## sparse - INS up-peaks 754 -> 32), a far larger change to MAGICAL than the
## DEG-table mismatch that motivates it. --peaks no filters genes only (variant
## "frq05genes"); the default filters both.
APPLY_PEAKS <- getarg("--peaks", "yes") != "no"
VARIANT <- getarg("--variant", if (APPLY_PEAKS) "frq05" else "frq05genes")

src <- file.path(MAGICAL_INPUT_DIR, magical_tag(REGION, CT))
dst <- paste0(src, "_", VARIANT)
if (!file.exists(file.path(src, "_DONE"))) stop("original MAGICAL input incomplete: ", src)
dir.create(dst, showWarnings = FALSE)

rna  <- as.data.frame(readRDS(file.path(OUT_ROOT, "03-DE_analysis", REGION, "results",
                      paste0(REGION, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))))
atac <- as.data.frame(readRDS(file.path(OUT_ROOT, "04-ATAC_DE", REGION, "data_ATAC",
                      paste0(REGION, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))))
for (nm in c("rna", "atac")) {
  x <- get(nm)
  if (!all(c("avg_CTR_frq", "avg_HIV_frq") %in% names(x)))
    stop(nm, " table has no avg_CTR_frq / avg_HIV_frq - cannot apply the expression filter")
}
rna  <- rna[as.character(rna$cluster_id) == CT, ]
atac <- atac[as.character(atac$cluster_id) == CT, ]

## MAGICAL's own rule, exactly as prepare_MAGICAL_input.R writes it
sel <- function(d, dir) switch(dir,
  all  = which(abs(d$logFC) >= 0.5 & d$p_adj.glb <= 0.05),
  up   = which(d$logFC >= 0.5 & d$p_adj.glb <= 0.05),
  down = which(d$logFC <= -0.5 & d$p_adj.glb <= 0.05))
expressed <- function(d) (d$avg_CTR_frq > FRQ) | (d$avg_HIV_frq > FRQ)
split_peaks <- function(ids) {
  s <- strsplit(as.character(ids), "-")
  rbindlist(lapply(s, function(x) as.data.table(list(Chromosome = x[1], Start = x[2], End = x[3]))))
}
wr <- function(x, f) write.table(x, file = file.path(dst, f), quote = FALSE,
                                 row.names = FALSE, col.names = FALSE, sep = "\t")
rd <- function(f) { p <- file.path(src, f); if (file.size(p) == 0) character(0) else readLines(p) }

GENE_F <- c(all = "Cell type candidate genes_DEGs.txt", up = "Cell type candidate genes_DEGs_up.txt",
            down = "Cell type candidate genes_DEGs_down.txt")
PEAK_F <- c(all = "all_de_peak_Cell type candidate peaks.txt", up = "up_de_peak_Cell type candidate peaks.txt",
            down = "down_de_peak_Cell type candidate peaks.txt")

## ---- self-check: the unfiltered selection must reproduce the files on disk ---------------
for (k in names(GENE_F)) {
  mine <- rna$gene[sel(rna, k)]
  if (!identical(mine, rd(GENE_F[[k]])))
    stop("unfiltered ", k, " genes do not match ", GENE_F[[k]], " (", length(mine), " vs ",
         length(rd(GENE_F[[k]])), ") - not the table the original run used")
}
for (k in names(PEAK_F)) {
  mine <- if (length(sel(atac, k))) apply(as.matrix(split_peaks(atac$gene[sel(atac, k)])), 1, paste, collapse = "\t") else character(0)
  if (!identical(unname(mine), rd(PEAK_F[[k]])))
    stop("unfiltered ", k, " peaks do not match ", PEAK_F[[k]])
}
message("[frq] self-check passed: unfiltered selection reproduces the original candidate files")

## ---- the filtered candidates ---------------------------------------------------------------
rep <- list()
for (k in names(GENE_F)) {
  i <- sel(rna, k); keep <- i[expressed(rna[i, ])]
  wr(rna$gene[keep], GENE_F[[k]])
  rep[[paste0("genes_", k)]] <- data.frame(set = paste0("genes_", k), before = length(i),
    after = length(keep), removed = paste(setdiff(rna$gene[i], rna$gene[keep]), collapse = ";"))
}
for (k in names(PEAK_F)) {
  i <- sel(atac, k); keep <- if (APPLY_PEAKS) i[expressed(atac[i, ])] else i
  wr(if (length(keep)) split_peaks(atac$gene[keep]) else data.table(), PEAK_F[[k]])
  rep[[paste0("peaks_", k)]] <- data.frame(set = paste0("peaks_", k), before = length(i),
    after = length(keep), removed = length(i) - length(keep))
}
report <- bind_rows(lapply(rep, function(x) { x$removed <- as.character(x$removed); x }))
utils::write.csv(report, file.path(dst, "filter_report.csv"), row.names = FALSE)
print(report[, c("set", "before", "after")], row.names = FALSE)

## ---- everything else: link to the original inputs ----------------------------------------
for (f in setdiff(list.files(src), c(GENE_F, PEAK_F, "_DONE"))) {
  to <- file.path(dst, f)
  if (!file.exists(to)) file.symlink(file.path(src, f), to)
}
file.create(file.path(dst, "_DONE"))
up_removed <- strsplit(report$removed[report$set == "genes_up"], ";")[[1]]
message(sprintf("[frq] %s %s: %d up-regulated candidate genes removed by the > %.0f%% filter",
                REGION, CT, length(up_removed), 100 * FRQ))
message("[frq] of note: ", paste(intersect(up_removed, c("CP", "SOCS3", "CIITA", "SERPINA3", "OSMR", "CHI3L1")),
                                 collapse = ", "))
message("[frq] wrote ", dst)
