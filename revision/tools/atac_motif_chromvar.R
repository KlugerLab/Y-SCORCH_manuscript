# Donor-aware chromVAR motif activity test, PWH vs PWoH, per region (Figure 2 motif panel).
#
# Method: the published panel ran chromVAR on each regional ATAC object, selected
#   motifs by over-representation in PWH-up DA peaks (FindMotifs, p.adj < 1e-10,
#   top 30 by fold enrichment) and plotted cell-level mean scores. That treats
#   nuclei as replicates (the reviewer objection to Fig 1E), and the selection is
#   an over-representation result, not a test of motif activity; with the
#   reprocessed DA peaks it is also weak (INS Astro has 32 up peaks, INS OPC 1).
#   Here, per region and cell type:
#   - peaks on standard chromosomes; JASPAR2020 CORE vertebrates, latest versions
#     (as published; the same set MAGICAL uses)
#   - AddMotifs on the region's peaks, then chromVAR (addGCBias ->
#     getBackgroundPeaks -> computeDeviations) on each cell type's nuclei, so the
#     GC/accessibility background is cell-type specific
#   - each donor's mean deviation z-score (donors with >= --min-cells nuclei),
#     then limma on the donor means with the ATAC DA design matrix
#     (HIV, CTR, Age, Sex, Race, PMI), contrast HIV - CTR, BH within cell type x
#     region. All motifs are tested; none is pre-selected by the same contrast.
#   - published-style over-representation (FindMotifs on up- and down-DA peaks)
#     is written alongside for comparison.
#   With --celltypes all, a cell type runs only if >= 3 PWH and >= 3 PWoH donors
#   have >= MIN_CELLS nuclei.
#
# Inputs:  OUT_ROOT/02-ATAC_processing/sampleQC_redo/<REGION>/data/signac_integrated_v2[_msn].RDS
#          OUT_ROOT/04-ATAC_DE/<REGION>/data_ATAC/{<REGION>_HIV-CTR_res_edgeR_filtered_withFreq.RDS, mm.rds}
# Outputs: OUT_ROOT/04-ATAC_DE/motif_chromvar/<REGION>/
#   <CT>_chromvar_donor_test.csv   motif, name, PWH-PWoH difference, t, P, FDR,
#                                  mean PWoH / PWH, n donors
#   <CT>_chromvar_donor_means.csv  donor x motif means (the data behind the test)
#   <CT>_chromvar_cells.rds        motif x cell deviation z-scores + cell meta
#   <CT>_findmotifs_{up,down}.csv  over-representation in DA peaks
#   motif_names.csv
# The ATAC objects are not modified or re-saved.
#
# Usage: revision/docker/run.sh revision/tools/atac_motif_chromvar.R --region INS
#          [--celltypes Astro,Micro,OPC | all] [--min-cells 20] [--threads 8] [--refresh]

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(JASPAR2020); library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38); library(chromVAR); library(motifmatchr)
  library(limma); library(dplyr); library(GenomeInfoDb)
})

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
REGION    <- toupper(getarg("--region", NA_character_))
CTS       <- trimws(strsplit(getarg("--celltypes", "Astro,Micro,OPC"), ",")[[1]])
MIN_CELLS <- as.integer(getarg("--min-cells", "20"))
REFRESH   <- "--refresh" %in% args
THREADS   <- as.integer(getarg("--threads", "8"))
if (is.na(REGION)) stop("--region is required")
BiocParallel::register(BiocParallel::MulticoreParam(THREADS))
set.seed(1)

OUT <- file.path(OUT_ROOT, "04-ATAC_DE", "motif_chromvar", REGION)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()

## ---- object: the one the ATAC DA used ----------------------------------------
fn <- if (REGION == "VST") "signac_integrated_v2_msn.RDS" else "signac_integrated_v2.RDS"
f  <- file.path(OUT_ROOT, "02-ATAC_processing/sampleQC_redo", REGION, "data", fn)
message("[motif] ", REGION, ": loading ", f)
obj <- readRDS(f)
DefaultAssay(obj) <- "ATAC"
md <- obj@meta.data
for (col in c("consensus_annotation_v3", "sample_ids"))
  if (!col %in% colnames(md)) stop("no '", col, "' column on the ATAC object")
## --celltypes all: every annotated type that could actually be tested. A type is kept only
## if >= 3 PWH and >= 3 PWoH donors have >= MIN_CELLS nuclei - below that the limma
## step is skipped anyway, and chromVAR on it is ~10 min of wasted compute.
## File names are sanitised because names like "L2/3 IT" and "D1/D2-Hybrid" contain
## a slash; names without special characters are unchanged.
safe_name <<- function(x) gsub("[^A-Za-z0-9_.+-]", "_", x)
all_ct <- unique(as.character(md$consensus_annotation_v3))
if (identical(tolower(CTS[1]), "all")) CTS <- sort(all_ct) else CTS <- intersect(CTS, all_ct)
cond_of <- setNames(as.character(pwh_relabel(md$condition_new)), rownames(md))
ct_ok <- vapply(CTS, function(ct) {
  i <- as.character(md$consensus_annotation_v3) == ct
  tb <- table(as.character(md$sample_ids)[i])
  don <- names(tb)[tb >= MIN_CELLS]
  cnd <- vapply(don, function(d) cond_of[rownames(md)[as.character(md$sample_ids) == d][1]], character(1))
  sum(cnd == "PWH") >= 3 && sum(cnd == "PWoH") >= 3
}, logical(1))
if (any(!ct_ok)) message("[motif] too few donors with >= ", MIN_CELLS, " nuclei, skipped: ",
                         paste(CTS[!ct_ok], collapse = ", "))
CTS <- CTS[ct_ok]
if (!REFRESH) {
  done <- vapply(CTS, function(ct) file.exists(file.path(OUT, paste0(safe_name(ct), "_chromvar_donor_means.csv"))),
                 logical(1))
  if (any(done)) message("[motif] already computed, kept as is (--refresh to redo): ",
                         paste(CTS[done], collapse = ", "))
  CTS <- CTS[!done]
}
if (!length(CTS)) { message("[motif] nothing left to compute"); quit(save = "no") }
message("[motif] cell types: ", paste(CTS, collapse = ", "))

## keep only the nuclei needed, then standard chromosomes (as published)
obj <- subset(obj, cells = colnames(obj)[as.character(md$consensus_annotation_v3) %in% CTS])
main_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep <- as.character(seqnames(granges(obj))) %in% main_chroms
obj[["ATAC"]] <- subset(obj[["ATAC"]], features = rownames(obj[["ATAC"]])[keep])
message("[motif] ", ncol(obj), " nuclei, ", sum(keep), " peaks on standard chromosomes")

## ---- motifs: JASPAR2020 CORE vertebrates, latest versions (as published) ------
pfm <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates",
                                            all_versions = FALSE))
obj <- AddMotifs(obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
if (!"GC.percent" %in% colnames(obj[["ATAC"]][[]]))
  obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)
motif_names <- unlist(obj[["ATAC"]]@motifs@motif.names)
utils::write.csv(data.frame(motif = names(motif_names), name = unname(motif_names)),
                 file.path(OUT, "motif_names.csv"), row.names = FALSE)
message("[motif] ", length(motif_names), " motifs matched to peaks (",
        format(round(difftime(Sys.time(), t0, units = "mins"), 1)), ")")

## ---- ATAC DA results and design ------------------------------------------------
da <- as.data.frame(readRDS(file.path(OUT_ROOT, "04-ATAC_DE", REGION, "data_ATAC",
                                      paste0(REGION, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))))
da <- da[(da$avg_CTR_frq > 0.05) | (da$avg_HIV_frq > 0.05), ]
mm <- readRDS(file.path(OUT_ROOT, "04-ATAC_DE", REGION, "data_ATAC", "mm.rds"))
stopifnot(all(c("HIV", "CTR") %in% colnames(mm)))

for (ct in CTS) {
  message("\n[motif] ===== ", REGION, " ", ct, " =====")
  sub <- subset(obj, cells = colnames(obj)[as.character(obj$consensus_annotation_v3) == ct])

  ## published-style over-representation in DA peaks (for comparison only)
  for (dir in c("up", "down")) {
    sel <- da$gene[da$cluster_id == ct & da$p_adj.glb < 0.05 &
                   (if (dir == "up") da$logFC > 0.25 else da$logFC < -0.25)]
    sel <- intersect(sel, rownames(sub))
    fm <- if (length(sel) >= 10) tryCatch(FindMotifs(sub, features = sel), error = function(e) NULL) else NULL
    if (is.null(fm)) fm <- data.frame()
    utils::write.csv(fm, file.path(OUT, paste0(safe_name(ct), "_findmotifs_", dir, ".csv")), row.names = FALSE)
    message("[motif] ", ct, " ", dir, "-DA peaks: ", length(sel), "; motifs over-represented at p.adj < 0.05: ",
            if (nrow(fm)) sum(fm$p.adjust < 0.05) else 0)
  }

  ## chromVAR on this cell type's nuclei.
  ## Signac 1.17 no longer exports RunChromVAR(), so this is the same pipeline written out with chromVAR
  ## itself, as Signac's wrapper did it: counts + peak ranges -> GC bias ->
  ## GC/accessibility-matched background peaks -> expectations -> deviations;
  ## keep the deviation z-scores. Peaks with no reads in this cell type are
  ## dropped first (chromVAR's background sampling needs non-zero peaks).
  ## seeded per cell type, not once per run: getBackgroundPeaks draws from the RNG,
  ## so with one global seed a cell type's z-scores would depend on the loop order
  set.seed(1)
  cnt <- GetAssayData(sub, assay = "ATAC", layer = "counts")
  nz  <- Matrix::rowSums(cnt) > 0
  se  <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = cnt[nz, ]),
                                                    rowRanges = granges(sub[["ATAC"]])[nz])
  se  <- chromVAR::addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg38)
  rd_ <- SummarizedExperiment::rowData(se); rd_$bias[is.na(rd_$bias)] <- 0
  SummarizedExperiment::rowData(se) <- rd_
  bg  <- chromVAR::getBackgroundPeaks(se, bias = rd_$bias)
  mot <- Motifs(sub[["ATAC"]])@data
  stopifnot(all(rownames(cnt)[nz] %in% rownames(mot)))
  mot <- mot[rownames(cnt)[nz], , drop = FALSE]                  # aligned by peak name
  dev <- chromVAR::computeDeviations(se, annotations = mot, background_peaks = bg,
                                     expectation = chromVAR::computeExpectations(se))
  z <- as.matrix(SummarizedExperiment::assays(dev)[["z"]])
  z[!is.finite(z)] <- 0
  message("[motif] ", ct, ": chromVAR on ", ncol(z), " nuclei, ", sum(nz), " peaks with reads, ",
          nrow(z), " motifs")
  cm <- sub@meta.data[, intersect(c("sample_ids", "sample_name", "condition_new", "consensus_annotation_v3"),
                                  colnames(sub@meta.data))]
  saveRDS(list(z = z, meta = cm), file.path(OUT, paste0(safe_name(ct), "_chromvar_cells.rds")))

  ## donor means, donors with >= MIN_CELLS nuclei that are in the DA design
  don <- as.character(sub$sample_ids)
  n_d <- table(don)
  use <- intersect(names(n_d)[n_d >= MIN_CELLS], rownames(mm))
  M <- sapply(use, function(d) rowMeans(z[, don == d, drop = FALSE]))
  utils::write.csv(data.frame(motif = rownames(M), name = motif_names[rownames(M)], M, check.names = FALSE),
                   file.path(OUT, paste0(safe_name(ct), "_chromvar_donor_means.csv")), row.names = FALSE)

  design <- mm[use, , drop = FALSE]
  design <- design[, colSums(design != 0) > 0, drop = FALSE]      # covariate levels absent here
  ne <- limma::nonEstimable(design)
  if (!is.null(ne)) { message("[motif] dropping non-estimable: ", paste(ne, collapse = ", "))
                      design <- design[, setdiff(colnames(design), ne), drop = FALSE] }
  n_hiv <- sum(design[, "HIV"] == 1); n_ctr <- sum(design[, "CTR"] == 1)
  message("[motif] donors with >= ", MIN_CELLS, " nuclei: ", length(use), " (PWH ", n_hiv, ", PWoH ", n_ctr, ")")
  ## A cell type with too few residual df (e.g. VST Vip, 6 donors, 3/3) makes
  ## eBayes fail with "No finite residual standard deviations"; it is skipped
  ## here, and its donor means (written above) can still go into the pooled model.
  rdf <- nrow(design) - qr(design)$rank
  if (n_hiv < 3 || n_ctr < 3 || rdf < 2) {
    message("[motif] too few donors per group (or ", rdf, " residual df) - per-region test skipped"); next }

  fit <- lmFit(M, design)
  fit <- eBayes(contrasts.fit(fit, makeContrasts(HIV - CTR, levels = design)))
  tt  <- topTable(fit, number = Inf, sort.by = "P")
  res <- data.frame(motif = rownames(tt), name = motif_names[rownames(tt)],
                    diff_PWH_minus_PWoH = tt$logFC, t = tt$t, P = tt$P.Value, FDR = tt$adj.P.Val,
                    mean_PWoH = rowMeans(M[rownames(tt), design[, "CTR"] == 1, drop = FALSE]),
                    mean_PWH  = rowMeans(M[rownames(tt), design[, "HIV"] == 1, drop = FALSE]),
                    n_PWH = n_hiv, n_PWoH = n_ctr, region = REGION, cell_type = ct)
  utils::write.csv(res, file.path(OUT, paste0(safe_name(ct), "_chromvar_donor_test.csv")), row.names = FALSE)
  message("[motif] ", ct, ": motifs at FDR < 0.05: ", sum(res$FDR < 0.05),
          " (higher in PWH ", sum(res$FDR < 0.05 & res$diff_PWH_minus_PWoH > 0), ", lower ",
          sum(res$FDR < 0.05 & res$diff_PWH_minus_PWoH < 0), ")")
  print(utils::head(res[, c("name", "diff_PWH_minus_PWoH", "t", "FDR")], 12), row.names = FALSE)
  rm(sub, z); invisible(gc())
}
message("\n[motif] ", REGION, " done in ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)),
        " -> ", OUT)
