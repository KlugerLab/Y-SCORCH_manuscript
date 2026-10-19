# ==============================================================================
# Figure 2F: TF expression vs TF motif activity, PWH vs PWoH.
#
# QUESTION  Is a TF whose motifs open up in PWH also induced at the RNA level in
#   the same cell type? Both axes are pooled across regions and donor-aware:
#     x  RNA log2FC of the TF gene - cross-region donor-aware DE
#        (03-DE_analysis/cross_region/donor_aware/cross_region_donor_aware_DE.csv,
#        limma-voom + duplicateCorrelation, average within-region HIV effect)
#     y  PWH - PWoH chromVAR deviation z-score of the TF's motif - pooled model,
#        CENTERED version (motif_chromvar/crossregion/<CT>_chromvar_crossregion.csv),
#        the same values as the Pooled column of fig2D
#   one point per TF tested in both
#
# MOTIF -> GENE (JASPAR2020 CORE vertebrates, 746 motifs)
#   heterodimer motifs ("FOS::JUN") are dropped - they cannot be assigned to one gene
#   "(var.N)" suffixes stripped and names upper-cased (mouse "Nfe2l2" -> NFE2L2)
#   a gene with several motifs keeps its canonical one (no var suffix, human name)
#
# SIGNIFICANCE  RNA: p_adj.glb < 0.05 and |log2FC| > 0.25 (the DEG rule used
#   throughout); motif: centered-model FDR < 0.05. Labels: TFs significant on BOTH
#   axes - a fixed rule, no forced labels.
#
# INPUT   OUT_ROOT/03-DE_analysis/cross_region/donor_aware/cross_region_donor_aware_DE.csv
#         OUT_ROOT/04-ATAC_DE/motif_chromvar/crossregion/<CT>_chromvar_crossregion.csv
# OUTPUT  OUT_ROOT/04-ATAC_DE/figures/fig2/fig2F_tf_rna_vs_motif<suffix>.png/.pdf/_caption.txt
#         + fig2F_tf_rna_vs_motif<suffix>.csv (every point)
#
# Usage: revision/docker/run.sh revision/tools/fig2_tf_rna_vs_motif.R
#          [--celltypes Astro,Micro] [--suffix <tag>]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }

MOT <- file.path(OUT_ROOT, "04-ATAC_DE", "motif_chromvar", "crossregion")
RNA <- file.path(OUT_ROOT, "03-DE_analysis", "cross_region", "donor_aware", "cross_region_donor_aware_DE.csv")
OUT <- file.path(OUT_ROOT, "04-ATAC_DE", "figures", "fig2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
## OPCs are not drawn: their hits were mostly FOX motifs (the AT/GC composition
## axis, see revision/tools/atac_motif_crossregion.R). Excitatory neurons are
## cortex-only and have no cross-region RNA fit; the four main interneuron types
## can be drawn as a supplement panel (--celltypes). Main panel = Astro + Micro.
## Pvalb/Sst/Lamp5 have no TF significant on both axes and no motif at FDR < 0.05
## in the pooled model, and Vip has one (TEF).
## A cell type is drawn only if its pooled chromVAR table exists; a few get a
## longer display name. File names are sanitised as in atac_motif_chromvar.R.
safe_name <- function(x) gsub("[^A-Za-z0-9_.+-]", "_", x)
NICE <- c(Astro = "Astrocytes", Micro = "Microglia", Oligo = "Oligodendrocytes",
          OPC = "OPCs", Endo = "Endothelial", PC = "Pericytes")
args_cts <- getarg("--celltypes", "Astro,Micro")
SUFFIX   <- getarg("--suffix", "")
cts <- trimws(strsplit(args_cts, ",")[[1]])
CTS <- setNames(ifelse(cts %in% names(NICE), NICE[cts], cts), cts)
CTS <- CTS[file.exists(file.path(MOT, paste0(safe_name(names(CTS)), "_chromvar_crossregion.csv")))]
message("[fig2F] cell types: ", paste(CTS, collapse = ", "))

rna <- utils::read.csv(RNA, stringsAsFactors = FALSE)
d <- bind_rows(lapply(names(CTS), function(ct) {
  m <- utils::read.csv(file.path(MOT, paste0(safe_name(ct), "_chromvar_crossregion.csv")), stringsAsFactors = FALSE)
  m <- m[!grepl("::", m$name), ]
  m$is_var   <- grepl("\\(var\\.[0-9]+\\)", m$name)
  m$is_human <- m$name == toupper(m$name)
  m$gene     <- toupper(sub("\\(var\\.[0-9]+\\)", "", m$name))
  m <- m %>% arrange(gene, is_var, !is_human, motif) %>% distinct(gene, .keep_all = TRUE)
  r <- rna[rna$cluster_id == ct, c("gene", "logFC", "p_adj.glb")]
  inner_join(m[, c("gene", "name", "motif", "diff_pooled", "FDR", "FDR_plain", "FDR_weighted")], r, by = "gene") %>%
    mutate(cell_type = ct, n_motif_genes = nrow(m))
}))
d <- d %>% mutate(
  rna_sig   = p_adj.glb < 0.05 & abs(logFC) > 0.25,
  motif_sig = FDR < 0.05,
  class = factor(case_when(rna_sig & motif_sig ~ "RNA and motif",
                           rna_sig ~ "RNA only", motif_sig ~ "Motif only", TRUE ~ "Neither"),
                 levels = c("RNA and motif", "RNA only", "Motif only", "Neither")),
  same_dir = sign(logFC) == sign(diff_pooled),
  lab = ifelse(class == "RNA and motif", gene, NA_character_))
utils::write.csv(d %>% dplyr::select(cell_type, gene, motif, motif_name = name, rna_log2FC = logFC,
                                     rna_FDR = p_adj.glb, motif_diff = diff_pooled, motif_FDR = FDR,
                                     motif_FDR_plain = FDR_plain, motif_FDR_weighted = FDR_weighted, class, same_dir),
                 file.path(OUT, paste0("fig2F_tf_rna_vs_motif", SUFFIX, ".csv")), row.names = FALSE)

st <- d %>% group_by(cell_type) %>%
  summarise(n = dplyr::n(), both = sum(class == "RNA and motif"),
            both_same = sum(class == "RNA and motif" & same_dir),
            rho = cor(logFC, diff_pooled, method = "spearman"),
            rho_rnasig = cor(logFC[rna_sig], diff_pooled[rna_sig], method = "spearman"),
            n_rnasig = sum(rna_sig), .groups = "drop")
print(as.data.frame(st), row.names = FALSE, digits = 2)
for (ct in names(CTS)) {
  x <- d[d$cell_type == ct & d$class == "RNA and motif", ]
  message("[fig2F] ", ct, " RNA and motif: ",
          paste(sprintf("%s (RNA %+.2f, motif %+.2f)", x$gene, x$logFC, x$diff_pooled)[order(-x$diff_pooled)], collapse = ", "))
}

strip_lab <- setNames(unname(CTS[st$cell_type]), st$cell_type)   # no counts in titles
d$cell_type <- factor(d$cell_type, levels = names(CTS))
d <- d %>% arrange(desc(class))                                   # significant points drawn on top
CLASS_COLS <- c(`RNA and motif` = "#B2182B", `RNA only` = "#4393C3", `Motif only` = "#F4A582", Neither = "grey82")
p <- ggplot(d, aes(x = logFC, y = diff_pooled, colour = class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey40") +
  geom_point(size = 1.1, alpha = 0.9, stroke = 0) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.1, fontface = "italic", colour = "black", na.rm = TRUE,
                           max.overlaps = Inf, box.padding = 0.3, min.segment.length = 0.1, force = 2,
                           segment.size = 0.2, segment.colour = "grey50", seed = 1, show.legend = FALSE) +
  scale_colour_manual(values = CLASS_COLS, name = "FDR < 0.05", drop = FALSE) +
  facet_wrap(~ cell_type, ncol = min(3, length(CTS)), scales = "free", labeller = as_labeller(strip_lab)) +
  labs(x = "log2 FC, TF gene expression (RNA)", y = "PWH − PWoH, TF motif activity (chromVAR z)") +
  theme_astro_pub(base_size = 7) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.3, "cm"), strip.text = element_text(size = 7.5, lineheight = 0.9)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2)))
f <- file.path(OUT, paste0("fig2F_tf_rna_vs_motif", SUFFIX, ".png"))
if (length(CTS) <= 2) save_pub(p, f, width = 6.0, height = 3.2) else
  save_pub(p, f, width = 7.2, height = 0.6 + 2.6 * ceiling(length(CTS) / 3))
astro_caption(f, paste(
  "Transcription-factor expression and motif activity in PWH vs PWoH. Each point is a TF with a JASPAR2020 motif",
  "that was tested for differential expression in that cell type: x, log2 fold change of the TF gene (donor-aware",
  "pseudobulk limma-voom, the three regions modelled together); y, PWH minus PWoH chromVAR deviation z-score of its",
  "motif (donor-aware limma, regions together, each sample's mean over all motifs removed). Heterodimer motifs are",
  "omitted. Colour: significance (RNA, FDR < 0.05 and |log2FC| > 0.25; motif, FDR < 0.05). TFs significant on",
  "both axes are labelled."))
message("[fig2F] done")
