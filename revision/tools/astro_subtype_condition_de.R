# ==============================================================================
# Donor-aware DE of PWH vs PWoH within each astrocyte subtype (input to Fig S15
# with astro_manuscript_figures.R --s15-de subtype).
#
# Method
#   The published S15 read `deg_results`, cell-level Wilcoxon (FindMarkers) tables
#   built in an interactive session and not defined in the Rmd; cell-level tests
#   are not donor-aware. This rebuilds the DEGs like the other revision DE:
#   donor-level pseudobulk on uncorrected counts (DE magnitudes are never taken
#   from CellBender-corrected counts).
#   The model is the primary cross-region model of
#   revision/tools/cross_region_donor_aware.R (limma-voom + duplicateCorrelation,
#   block = donor), run on one astrocyte subtype at a time instead of one cell type:
#     counts    uncorrected_counts() for the subtype's nuclei (revision/config/paths.R)
#     samples   pseudobulk per region x donor, rows "<REGION>.<donor>"; a sample
#               needs >= 10 nuclei (pbDS min_cells = 10)
#     design    the cross-region DE design (03-DE_analysis/cross_region/data/mm.rds):
#               region x condition means + Age + Sex + Race + PMI
#     contrast  average within-region HIV effect,
#               (INS.HIV-INS.CTR)/3 + (PFC.HIV-PFC.CTR)/3 + (VST.HIV-VST.CTR)/3
#     fit       filterByExpr -> TMM -> voom -> duplicateCorrelation(block = donor),
#               twice (limma idiom) -> lmFit(block, correlation) -> eBayes
#     FDR       p_adj.loc within subtype; p_adj.glb = BH across all subtypes jointly
#               (the p_adj.glb convention of every other table)
#   Detection (pct_PWH / pct_PWoH) is the fraction of the subtype's nuclei with a
#   non-zero uncorrected count, so the S15 code's pct.1 > 0.05 filter still applies.
#   Subtypes with fewer than 8 qualifying samples, or a missing region/condition,
#   are skipped.
#
# Inputs   OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled_no1040vst.RDS
#          OUT_ROOT/03-DE_analysis/cross_region/data/mm.rds
# Outputs  OUT_ROOT/06-astrocyte_analysis/subtype_condition_de/
#            astro_subtype_PWH_vs_PWoH.csv          every gene x subtype
#            astro_subtype_PWH_vs_PWoH_summary.csv
#
# Usage: revision/docker/run.sh revision/tools/astro_subtype_condition_de.R
#          [--subtypes Homeostatic,Reactive,Signaling,Synaptic]
# ==============================================================================

suppressPackageStartupMessages({ library(Seurat); library(limma); library(edgeR); library(dplyr) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
SUBTYPES <- trimws(strsplit(getarg("--subtypes", "Homeostatic,Reactive,Signaling,Synaptic"), ",")[[1]])
MIN_CELLS <- 10
CONTRAST_STR <- "(INS.HIV-INS.CTR)/3+(PFC.HIV-PFC.CTR)/3+(VST.HIV-VST.CTR)/3"
OUT <- file.path(OUT_ROOT, "06-astrocyte_analysis", "subtype_condition_de")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(OUT_ROOT, "06-astrocyte_analysis", "data", "astro_all_clusters_labeled_no1040vst.RDS"))
md <- obj@meta.data
if (!"region" %in% colnames(md))
  obj$region <- sub("^.*-([A-Z]+)$", "\\1", colnames(obj))    # cell names end in -<REGION>
obj <- obj[, as.character(obj$subcluster_labels) %in% SUBTYPES]
message("[sub-de] nuclei in ", paste(SUBTYPES, collapse = "/"), ": ", ncol(obj))

cts <- uncorrected_counts(obj)                                  # genes x cells, obj's cell order
stopifnot(identical(colnames(cts), colnames(obj)))
md  <- obj@meta.data
mm  <- readRDS(file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "data", "mm.rds"))
cond <- as.character(pwh_relabel(md$condition_new))

run_one <- function(st) {
  i   <- as.character(md$subcluster_labels) == st
  key <- paste0(toupper(as.character(md$region[i])), ".", as.character(md$sample_ids[i]))
  f   <- factor(key)
  pb  <- as.matrix(cts[, i, drop = FALSE] %*% Matrix::sparse.model.matrix(~ 0 + f))
  colnames(pb) <- levels(f)
  nck <- table(key)[colnames(pb)]
  keep_s <- as.integer(nck) >= MIN_CELLS & colnames(pb) %in% rownames(mm)
  sid <- colnames(pb)[keep_s]
  if (length(sid) < 8) { message("[sub-de] ", st, ": only ", length(sid), " samples - skipped"); return(NULL) }

  d <- mm[sid, , drop = FALSE]
  d <- d[, colSums(abs(d)) > 0, drop = FALSE]
  if (qr(d)$rank < ncol(d)) {
    ne <- limma::nonEstimable(d); d <- d[, setdiff(colnames(d), ne), drop = FALSE]
  }
  ctr <- tryCatch(limma::makeContrasts(contrasts = CONTRAST_STR, levels = d), error = function(e) NULL)
  if (is.null(ctr)) { message("[sub-de] ", st, ": contrast not estimable (a region/condition missing)"); return(NULL) }

  dge <- DGEList(counts = pb[, sid, drop = FALSE])
  keep_g <- filterByExpr(dge, d)
  dge <- calcNormFactors(dge[keep_g, , keep.lib.sizes = FALSE])
  blk <- sub("^[^.]+\\.", "", sid)
  v   <- voom(dge, d)
  dc1 <- duplicateCorrelation(v, d, block = blk)
  v2  <- voom(dge, d, block = blk, correlation = dc1$consensus.correlation)
  rho <- duplicateCorrelation(v2, d, block = blk)$consensus.correlation
  fit <- eBayes(contrasts.fit(lmFit(v2, d, block = blk, correlation = rho), ctr))
  tt  <- topTable(fit, number = Inf, sort.by = "none")

  g   <- rownames(tt)
  pct <- function(sel) Matrix::rowMeans(cts[g, i & sel, drop = FALSE] > 0)
  n_cond <- tapply(cond[i], sub("^[^.]+\\.", "", key), function(x) x[1])
  message(sprintf("[sub-de] %s: %d samples, %d donors (PWH %d, PWoH %d), %d genes tested, rho %.3f",
                  st, length(sid), length(unique(blk)),
                  sum(n_cond[unique(blk)] == "PWH", na.rm = TRUE),
                  sum(n_cond[unique(blk)] == "PWoH", na.rm = TRUE), length(g), rho))
  data.frame(cluster_id = st, gene = g, logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
             p_val = tt$P.Value, p_adj.loc = tt$adj.P.Val,
             pct_PWH = pct(cond == "PWH"), pct_PWoH = pct(cond == "PWoH"),
             rho = rho, n_samples = length(sid), n_donors = length(unique(blk)),
             stringsAsFactors = FALSE)
}

R <- bind_rows(lapply(SUBTYPES, run_one))
R$p_adj.glb <- p.adjust(R$p_val, method = "BH")                  # across subtypes jointly
utils::write.csv(R, file.path(OUT, "astro_subtype_PWH_vs_PWoH.csv"), row.names = FALSE)
S <- R %>% group_by(cluster_id) %>%
  summarise(genes_tested = dplyr::n(),
            up_FDR05_lfc0.5   = sum(p_adj.glb < 0.05 & logFC >  0.5),
            down_FDR05_lfc0.5 = sum(p_adj.glb < 0.05 & logFC < -0.5),
            n_donors = n_donors[1], rho = rho[1], .groups = "drop")
utils::write.csv(S, file.path(OUT, "astro_subtype_PWH_vs_PWoH_summary.csv"), row.names = FALSE)
print(as.data.frame(S), row.names = FALSE)
message("[sub-de] done -> ", OUT)
