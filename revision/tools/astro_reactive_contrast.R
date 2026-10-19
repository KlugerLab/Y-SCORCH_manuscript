# ==============================================================================
# Donor-paired DE and GO enrichment of Reactive-PWH against Reactive astrocytes,
# compared with the Fig1E contrast (Reactive-PWH vs all other astrocytes).
#
# Method
#   Fig1E compares Reactive-PWH against every other astrocyte, so its gene list
#   mixes what makes a cell reactive at all with what is specific to the PWH state.
#   Reactive astrocytes exist in both groups, so using them as the comparator
#   subtracts the generic reactivity programme and leaves the part particular to
#   Reactive-PWH.
#   DE: astro_subtype_de() (revision/tools/astro_subtype_de_donor.R) pseudobulks
#   each donor's Reactive-PWH and Reactive nuclei and tests them paired within
#   donor (edgeR QL). Donors need >= 10 cells in both halves: all 20 PWH donors and
#   4 of 13 PWoH donors (HCTMU 141, HCTMQ 40, HCTKN 19, HCTNZ 11 Reactive-PWH
#   cells); the other PWoH donors have 0-8, the population's defining feature.
#   Two donor sets:
#   all   every donor with both halves (24). The PWoH donors act as a check: if the
#         signature holds where the state occurs without HIV, it is a cell state
#         rather than a disease effect.
#   pwh   the 20 PWH donors only; cleaner interpretation, unaffected by the few
#         PWoH cells whose label is least certain.
#   Each is the sensitivity analysis for the other.
#   GO is run as for Fig1E - up-regulated genes, FDR < 0.05, logFC > 0.5,
#   background = the genes tested in this contrast, BP, BH, p < 0.01 - so the only
#   difference from Fig1E is the comparator.
#
# Inputs   OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled_no1040vst.RDS
#          OUT_ROOT/06-astrocyte_analysis/astro_Reactive_PWH_donorDE_no1040vst.csv (Fig1E DE)
# Outputs  OUT_ROOT/06-astrocyte_analysis/
#            astro_Reactive_PWH_vs_Reactive_donorDE_no1040vst.csv  (+ _counts, _heatmap;
#            PWH-only run under pwh_only/)
#          OUT_ROOT/06-astrocyte_analysis/figures/Fig1E_vs_reactive/
#            go_<set>.csv, up_genes_<set>.csv  enrichment and up genes per donor set
#            gene_overlap_vs_published.csv     gene list vs Fig1E's
#            term_comparison.csv               terms gained, kept, lost vs Fig1E
#            fig1e_vs_reactive.png/.pdf
#
# Usage: revision/docker/run.sh revision/tools/astro_reactive_contrast.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
## The input object already excludes donor 1040's VST cells, but .sfxname() takes
## the suffix from this global - unset, the tables would be named as if they were
## the 1040-included run. Set it before sourcing the DE helper.
ASTRO_OUT_SUFFIX <- "_no1040vst"
TOOLS <- file.path(SCORCH_REPO, "revision/tools")
for (f in c("astro_pub_theme", "astro_subtype_de_donor")) source(file.path(TOOLS, paste0(f, ".R")))
suppressPackageStartupMessages({
  library(Seurat); library(clusterProfiler); library(org.Hs.eg.db)
  library(ggplot2); library(dplyr)
})

D   <- file.path(OUT_ROOT, "06-astrocyte_analysis")
OUT <- file.path(D, "figures", "Fig1E_vs_reactive")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CLEAN <- c("Homeostatic", "Reactive", "Reactive-PWH", "Signaling",
           "Synaptic", "Micro-Astro", "Mural-Astro")
PCUT <- 0.01; QCUT <- 0.05; FDR <- 0.05; LFC <- 0.5

IN <- file.path(D, "data", "astro_all_clusters_labeled_no1040vst.RDS")
if (!file.exists(IN)) stop("no astrocyte object at ", IN)
message("[vsR] loading ", IN)
obj <- readRDS(IN)
Idents(obj) <- "subcluster_labels"
obj <- subset(obj, idents = intersect(CLEAN, unique(as.character(obj$subcluster_labels))))
message("[vsR] clean subset: ", ncol(obj), " cells")

## ---- the two paired runs ----------------------------------------------------------------
runs <- list()
runs$all <- astro_subtype_de(obj, subtype = "Reactive-PWH", vs_subtype = "Reactive")
runs$pwh <- astro_subtype_de(obj, subtype = "Reactive-PWH", vs_subtype = "Reactive",
                             restrict_to_condition = "HIV",
                             outdir = file.path(D, "pwh_only"))

## ---- enrichment, identical to Fig1E apart from the comparator ---------------------------
go_of <- function(de) {
  g <- de$gene[de$FDR < FDR & de$logFC > LFC]
  message("[vsR]   up genes: ", length(g), " of ", nrow(de), " tested")
  if (length(g) < 10) return(list(genes = g, res = NULL))
  e <- try(enrichGO(gene = g, universe = unique(de$gene), keyType = "SYMBOL",
                    OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = PCUT, qvalueCutoff = QCUT, readable = TRUE), silent = TRUE)
  list(genes = g, res = if (inherits(e, "try-error")) NULL else as.data.frame(e))
}

pub <- utils::read.csv(file.path(D, "astro_Reactive_PWH_donorDE_no1040vst.csv"),
                       stringsAsFactors = FALSE)
pub_go <- go_of(pub)                      # the current Fig1E, recomputed identically
message("[vsR] current Fig1E (vs all other astrocytes): ",
        if (is.null(pub_go$res)) 0 else nrow(pub_go$res), " terms")

out <- list()
for (nm in names(runs)) {
  de <- runs[[nm]]$de
  message("\n[vsR] === ", nm, " donors: Reactive-PWH vs Reactive ===")
  g <- go_of(de)
  if (!is.null(g$res)) {
    utils::write.csv(g$res, file.path(OUT, paste0("go_", nm, ".csv")), row.names = FALSE)
    message("[vsR]   top terms: ", paste(head(g$res$Description, 5), collapse = " | "))
  }
  out[[nm]] <- g
  utils::write.csv(data.frame(gene = g$genes),
                   file.path(OUT, paste0("up_genes_", nm, ".csv")), row.names = FALSE)
}

## ---- how the gene lists relate to the published one --------------------------------------
ov <- dplyr::bind_rows(lapply(names(out), function(nm) {
  a <- out[[nm]]$genes; b <- pub_go$genes
  data.frame(donor_set = nm, n_genes = length(a), n_fig1e = length(b),
             shared = length(intersect(a, b)),
             pct_of_new_in_fig1e = round(100 * length(intersect(a, b)) / max(length(a), 1), 1),
             jaccard = round(length(intersect(a, b)) / length(union(a, b)), 3))
}))
print(ov, row.names = FALSE)
utils::write.csv(ov, file.path(OUT, "gene_overlap_vs_published.csv"), row.names = FALSE)

## ---- terms gained / kept / lost -----------------------------------------------------------
pub_terms <- if (is.null(pub_go$res)) character(0) else pub_go$res$Description
cmp <- dplyr::bind_rows(lapply(names(out), function(nm) {
  t2 <- if (is.null(out[[nm]]$res)) character(0) else out[[nm]]$res$Description
  data.frame(donor_set = nm,
             term = union(pub_terms, t2),
             in_fig1e = union(pub_terms, t2) %in% pub_terms,
             in_vs_reactive = union(pub_terms, t2) %in% t2)
}))
cmp$status <- ifelse(cmp$in_fig1e & cmp$in_vs_reactive, "kept",
                     ifelse(cmp$in_vs_reactive, "gained", "lost"))
utils::write.csv(cmp, file.path(OUT, "term_comparison.csv"), row.names = FALSE)
message("\n[vsR] terms by status:"); print(table(cmp$donor_set, cmp$status))

## ---- figure: top terms of each contrast, side by side --------------------------------------
top_of <- function(r, lab, n = 12) {
  if (is.null(r) || !nrow(r)) return(NULL)
  r <- r[order(r$p.adjust), ][seq_len(min(n, nrow(r))), ]
  data.frame(contrast = lab, Description = r$Description, p.adjust = r$p.adjust, Count = r$Count)
}
pf <- dplyr::bind_rows(
  top_of(pub_go$res, "vs all other astrocytes\n(current Fig1E)"),
  top_of(out$all$res, "vs Reactive\n(24 donors)"),
  top_of(out$pwh$res, "vs Reactive\n(20 PWH donors)"))
if (nrow(pf)) {
  lv <- pf %>% group_by(Description) %>% summarise(p = min(p.adjust)) %>% arrange(p) %>% pull(Description)
  pf$Description <- factor(pf$Description, levels = rev(lv))
  p <- ggplot(pf, aes(contrast, Description, size = -log10(p.adjust), colour = Count)) +
    geom_point() +
    scale_colour_viridis_c(name = "genes", option = "D") +
    scale_size_continuous(name = expression(-log[10]~adj.~italic(p)), range = c(1, 4.5)) +
    labs(x = NULL, y = NULL) + theme_astro_pub(7) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")
  ggsave(file.path(OUT, "fig1e_vs_reactive.png"), p, width = 8.4, height = 6.6, dpi = 600, bg = "white")
  ggsave(file.path(OUT, "fig1e_vs_reactive.pdf"), p, width = 8.4, height = 6.6,
         device = grDevices::cairo_pdf, bg = "white")
}
message("[vsR] wrote ", OUT)
