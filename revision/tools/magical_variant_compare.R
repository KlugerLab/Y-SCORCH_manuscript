# ==============================================================================
# How much the MAGICAL astrocyte network changes with the expression-filtered
# candidate variants (revision/tools/magical_frq_variant.R).
#
# Compares each variant's astrocyte circuits with the primary run, per region:
#   size       genes, peaks, TFs
#   overlap    gene and TF-gene pair Jaccard vs primary. Benchmark: re-running
#              the primary on identical inputs with another seed gives Jaccard
#              0.91 (PFC), so ~0.9 is "unchanged within noise".
#   TF rule    the five figure TFs (STAT3, STAT1, IRF1, MAFG, CEBPD): up-regulated
#              targets, using the same definition as revision/tools/magical_tf_selection.R
#              (region's astrocyte DE, logFC > 0 & p_adj.glb < 0.05), and whether
#              each still has >= 5. (The TF's own cross-region DE does not change.)
#   network    the pills figure's content: edges from those five TFs to
#              up-regulated genes - retained, lost, gained
#   markers    which TFs bind SERPINA3 and OSMR
#
# Inputs: MAGICAL_OUTPUT_DIR/[variants/<variant>/]<R>_AstroMAGICAL_selected_regulatory_circuits.txt
#         OUT_ROOT/03-DE_analysis/<R>/results/<R>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
# Output: MAGICAL_OUTPUT_DIR/variants/compare_variants.csv, compare_variants_network_genes.csv
#
# Usage: revision/docker/run.sh revision/tools/magical_variant_compare.R [--variants frq05genes,frq05]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(dplyr) })
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
VARS <- strsplit(getarg("--variants", "frq05genes,frq05"), ",")[[1]]
CT   <- "Astro"
FIG_TFS <- c("STAT3", "STAT1", "IRF1", "MAFG", "CEBPD")
MIN_UP  <- 5

read_circ <- function(f) {
  if (!file.exists(f)) return(NULL)
  L <- readLines(f)[-1]; L <- L[nzchar(trimws(L))]
  bind_rows(lapply(strsplit(L, "\t"), function(x) {
    tf <- strsplit(x[8], ",")[[1]]; tf <- tf[!grepl("^\\(", tf)]
    data.frame(gene = x[1], peak = paste(x[4], x[5], x[6], sep = ":"), tf = tf,
               stringsAsFactors = FALSE)
  }))
}
jac <- function(a, b) round(length(intersect(a, b)) / max(1, length(union(a, b))), 2)

MO  <- MAGICAL_OUTPUT_DIR
out <- list(); net_rows <- list()
for (R in c("PFC", "INS", "VST")) {
  base <- sprintf("%s_%sMAGICAL_selected_regulatory_circuits.txt", R, CT)
  prim <- read_circ(file.path(MO, base))
  reg  <- as.data.frame(readRDS(file.path(OUT_ROOT, "03-DE_analysis", R, "results",
                                sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq.RDS", R))))
  reg  <- reg[as.character(reg$cluster_id) == CT, ]
  up_genes <- unique(reg$gene[reg$logFC > 0 & reg$p_adj.glb < 0.05])

  sets <- c(list(primary = prim), setNames(lapply(VARS, function(v)
    read_circ(file.path(MO, "variants", v, base))), VARS))
  pairs_of <- function(d) unique(paste(d$tf, d$gene))
  net_of   <- function(d) unique(d[d$tf %in% FIG_TFS & d$gene %in% up_genes, c("tf", "gene")])

  for (v in names(sets)) {
    d <- sets[[v]]
    if (is.null(d)) { message("[cmp] ", R, " ", v, ": not finished"); next }
    upt <- vapply(FIG_TFS, function(t) length(intersect(unique(d$gene[d$tf == t]), up_genes)), integer(1))
    net <- net_of(d); pnet <- net_of(prim)
    out[[paste(R, v)]] <- data.frame(
      region = R, version = v,
      genes = n_distinct(d$gene), peaks = n_distinct(d$peak), tfs = n_distinct(d$tf),
      gene_jaccard_vs_primary = jac(unique(d$gene), unique(prim$gene)),
      pair_jaccard_vs_primary = jac(pairs_of(d), pairs_of(prim)),
      figure_edges = nrow(net), figure_genes = n_distinct(net$gene),
      figure_genes_kept = length(intersect(unique(net$gene), unique(pnet$gene))),
      t(setNames(upt, paste0(FIG_TFS, "_up_targets"))),
      tfs_passing_min5 = paste(FIG_TFS[upt >= MIN_UP], collapse = ";"),
      SERPINA3_tfs = paste(sort(unique(d$tf[d$gene == "SERPINA3"])), collapse = ";"),
      SERPINA3_fig_tfs = paste(sort(intersect(unique(d$tf[d$gene == "SERPINA3"]), FIG_TFS)), collapse = ";"),
      OSMR_fig_tfs = paste(sort(intersect(unique(d$tf[d$gene == "OSMR"]), FIG_TFS)), collapse = ";"),
      check.names = FALSE)
    if (v != "primary") {
      g_now <- unique(net$gene); g_was <- unique(pnet$gene)
      net_rows[[paste(R, v)]] <- data.frame(region = R, version = v,
        lost = paste(sort(setdiff(g_was, g_now)), collapse = ";"),
        gained = paste(sort(setdiff(g_now, g_was)), collapse = ";"))
    }
  }
}
res <- bind_rows(out); nets <- bind_rows(net_rows)
OUTD <- file.path(MO, "variants"); dir.create(OUTD, showWarnings = FALSE)
utils::write.csv(res,  file.path(OUTD, "compare_variants.csv"), row.names = FALSE)
utils::write.csv(nets, file.path(OUTD, "compare_variants_network_genes.csv"), row.names = FALSE)
options(width = 220)
message("\n[cmp] size and overlap (seed-replicate benchmark: Jaccard ~0.91)")
print(res[, c("region", "version", "genes", "peaks", "tfs", "gene_jaccard_vs_primary",
              "pair_jaccard_vs_primary", "figure_edges", "figure_genes", "figure_genes_kept")], row.names = FALSE)
message("\n[cmp] figure TFs: up-regulated targets (rule needs >= ", MIN_UP, ")")
print(res[, c("region", "version", paste0(FIG_TFS, "_up_targets"), "tfs_passing_min5")], row.names = FALSE)
message("\n[cmp] reactive markers")
print(res[, c("region", "version", "SERPINA3_fig_tfs", "OSMR_fig_tfs", "SERPINA3_tfs")], row.names = FALSE)
message("\n[cmp] figure genes lost / gained vs primary")
print(nets, row.names = FALSE)
