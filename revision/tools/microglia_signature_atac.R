# ==============================================================================
# Microglial signatures in chromatin accessibility: an orthogonal, donor-level,
# cell-type-specific check of the microglial activation/metabolic program
# (reviewer request for donor-level validation).
#
# Method: same nuclei (multiome), independent molecular layer. Gene activity =
#   ATAC fragments over each gene body + 2 kb upstream (Signac GeneActivity),
#   aggregated per donor within microglia by 04-ATAC_DE/<R>/01-DifferentialGeneACT.
#   No RNA counts enter the score.
#   Statistics are identical to the RNA signatures (revision/tools/microglia_signature_score.R,
#   unchanged, on this pseudobulk): per-donor score (mean z of logCPM over the
#   set), camera (competitive), fry (self-contained), design ~ condition + Age +
#   Sex + Race + PMI per region, and the cross-region mixed model with
#   (1 | donor). Gene sets: revision/config/microglia_signatures.R. Only the
#   modality differs. The *_dropAmbient outputs are written by the shared
#   function but are not meaningful for ATAC (ambient RNA is an RNA artefact).
#   Donors: >= 10 microglia nuclei in ATAC (same min_cells as the RNA
#   pseudobulk); fewer than RNA: INS 21 (RNA 28), PFC 29 (29), VST 20 (26).
#   Compared against the saved RNA result:
#     per region   camera direction/FDR in each modality, PWH-PWoH score
#                  difference, and Spearman of per-donor ATAC vs RNA score
#     cross-region the two mixed-model estimates side by side
#     per gene     microglia DE logFC in RNA vs gene-activity DE logFC
#
# Inputs: OUT_ROOT/04-ATAC_DE/<R>/data_GeneACT/pseudo_bulk.rds and
#         <R>_HIV-CTR_res_edgeR_filtered_withFreq.RDS;
#         OUT_ROOT/03-DE_analysis/<R>/results/<R>_HIV-CTR_res_edgeR_filtered_withFreq.RDS;
#         OUT_ROOT/03-DE_analysis/microglia_signatures/microglia_signature_results.rds
# Output: OUT_ROOT/03-DE_analysis/microglia_signatures/atac_gene_activity/
#
# Usage: revision/docker/run.sh revision/tools/microglia_signature_atac.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(SingleCellExperiment); library(muscat); library(dplyr)
})
source(file.path(SCORCH_REPO, "revision/tools/microglia_signature_score.R"))

REGIONS <- c("INS", "PFC", "VST")
OUTDIR  <- file.path(MICRO_OUT, "atac_gene_activity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
RNA     <- readRDS(file.path(MICRO_OUT, "microglia_signature_results.rds"))
SIGS    <- names(MICROGLIA_SIGNATURES)

## ---- ATAC gene-activity signatures, same functions as RNA ------------------------
res <- list()
for (R in REGIONS) {
  message("\n########## ", R, " - microglia, ATAC gene activity ##########")
  pbs  <- readRDS(file.path(OUT_ROOT, "04-ATAC_DE", R, "data_GeneACT", "pseudo_bulk.rds"))
  m    <- assay(pbs, MICRO_LABEL)
  nc   <- muscat:::.n_cells(pbs)[MICRO_LABEL, colnames(m)]
  keep <- colnames(m)[nc >= 10]
  message("[atac] ", R, ": ", length(keep), " of ", ncol(m),
          " donors with >= 10 microglia nuclei in ATAC")
  meta <- .microglia_attach_meta(data.frame(donor = keep, n_cells = as.integer(nc[keep]),
                                            region = R, stringsAsFactors = FALSE))
  pb <- list(counts = as.matrix(m[, keep, drop = FALSE]), meta = meta, region = R)
  res[[R]] <- microglia_score_region(pb, outdir = OUTDIR)
}
message("\n########## CROSS-REGION (ATAC) ##########")
cr_atac <- microglia_cross_region(res, outdir = OUTDIR)

## ---- concordance with RNA, per region ----------------------------------------------
cam_get <- function(tab, s, col) if (!is.null(tab) && s %in% rownames(tab)) tab[s, col] else NA
score_diff <- function(sc, s) {
  if (!s %in% colnames(sc)) return(NA_real_)
  mean(sc[[s]][sc$condition == "PWH"], na.rm = TRUE) - mean(sc[[s]][sc$condition == "PWoH"], na.rm = TRUE)
}
conc <- bind_rows(lapply(REGIONS, function(R) {
  a <- res[[R]]; r <- RNA$per_region[[R]]
  bind_rows(lapply(SIGS, function(s) {
    d <- merge(a$scores[, c("donor", s)], r$scores[, c("donor", s)], by = "donor",
               suffixes = c("_atac", "_rna"))
    d <- d[stats::complete.cases(d), , drop = FALSE]
    ct <- if (nrow(d) >= 5)
      suppressWarnings(stats::cor.test(d[[paste0(s, "_atac")]], d[[paste0(s, "_rna")]],
                                       method = "spearman")) else NULL
    data.frame(region = R, signature = s,
               rna_camera_dir  = cam_get(r$tests$camera, s, "Direction"),
               rna_camera_FDR  = cam_get(r$tests$camera, s, "FDR"),
               rna_score_diff  = round(score_diff(r$scores, s), 3),
               atac_camera_dir = cam_get(a$tests$camera, s, "Direction"),
               atac_camera_FDR = cam_get(a$tests$camera, s, "FDR"),
               atac_fry_FDR    = cam_get(a$tests$fry, s, "FDR"),
               atac_score_diff = round(score_diff(a$scores, s), 3),
               atac_genes_tested = if (s %in% names(a$tests$n_genes_tested)) a$tests$n_genes_tested[[s]] else NA,
               n_shared_donors = nrow(d),
               donor_spearman_rho = if (!is.null(ct)) round(unname(ct$estimate), 3) else NA,
               donor_spearman_p   = if (!is.null(ct)) signif(ct$p.value, 3) else NA,
               stringsAsFactors = FALSE)
  }))
}))
utils::write.csv(conc, file.path(OUTDIR, "rna_vs_atac_per_region.csv"), row.names = FALSE)

## ---- cross-region, side by side ----------------------------------------------------
cr_both <- merge(RNA$cross_region[, c("signature", "estimate", "p", "FDR", "n_donors")],
                 if (is.null(cr_atac)) data.frame(signature = character(0)) else
                   cr_atac[, c("signature", "estimate", "p", "FDR", "n_donors")],
                 by = "signature", all = TRUE, suffixes = c("_rna", "_atac"))
cr_both <- cr_both[order(cr_both$p_rna), ]
utils::write.csv(cr_both, file.path(OUTDIR, "rna_vs_atac_cross_region.csv"), row.names = FALSE)

## ---- per gene: RNA DE vs gene-activity DE in microglia --------------------------
sig_genes <- unique(unlist(lapply(MICROGLIA_SIGNATURES, function(x) x$genes)))
set_of <- vapply(sig_genes, function(g)
  paste(SIGS[vapply(MICROGLIA_SIGNATURES, function(x) g %in% x$genes, logical(1))], collapse = "; "),
  character(1))
gl <- bind_rows(lapply(REGIONS, function(R) {
  pick <- function(f, pre) {
    if (!file.exists(f)) return(data.frame(gene = character(0)))
    x <- as.data.frame(readRDS(f))
    x <- x[as.character(x$cluster_id) == MICRO_LABEL & x$gene %in% sig_genes, c("gene", "logFC", "p_adj.glb")]
    names(x)[2:3] <- paste0(pre, c("logFC", "p_adj.glb")); x
  }
  rn <- pick(file.path(OUT_ROOT, "03-DE_analysis", R, "results", paste0(R, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS")), "rna_")
  ga <- pick(file.path(OUT_ROOT, "04-ATAC_DE", R, "data_GeneACT", paste0(R, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS")), "atac_")
  out <- merge(rn, ga, by = "gene", all = TRUE)
  if (nrow(out)) out$region <- R
  out
}))
gl$sets <- set_of[gl$gene]
gl <- gl[order(gl$region, gl$rna_p_adj.glb), c("region", "gene", "sets", "rna_logFC", "rna_p_adj.glb", "atac_logFC", "atac_p_adj.glb")]
utils::write.csv(gl, file.path(OUTDIR, "rna_vs_atac_gene_level.csv"), row.names = FALSE)

## ---- report ------------------------------------------------------------------------
cat("\n=== PER REGION: RNA vs ATAC gene activity (camera = competitive; score diff = PWH - PWoH mean z) ===\n")
print(conc, row.names = FALSE)
cat("\n=== CROSS-REGION mixed model, condition PWH (donor random effect) ===\n")
print(cr_both, row.names = FALSE)
both <- gl[!is.na(gl$rna_logFC) & !is.na(gl$atac_logFC), ]
sig  <- both[both$rna_p_adj.glb < 0.05, ]
cat(sprintf("\n=== GENE LEVEL: %d signature gene x region pairs in both modalities; among %d RNA-significant, same sign in ATAC: %d (%.0f%%); ATAC also significant: %d ===\n",
            nrow(both), nrow(sig), sum(sign(sig$rna_logFC) == sign(sig$atac_logFC)),
            100 * mean(sign(sig$rna_logFC) == sign(sig$atac_logFC)), sum(sig$atac_p_adj.glb < 0.05, na.rm = TRUE)))
print(sig, row.names = FALSE)
message("\n[atac] wrote ", OUTDIR)
message("[atac] done")
