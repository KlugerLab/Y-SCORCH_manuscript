# ==============================================================================
# Figure 2D: TF motif activity (chromVAR) in PWH vs PWoH, per cell type and region.
#
# Reports motif activity as effect sizes; significance is marked but not required
# for a row to be shown. Reads the pooled, donor-aware chromVAR results
#   OUT_ROOT/04-ATAC_DE/motif_chromvar/crossregion/<CT>_chromvar_crossregion.csv
# (revision/tools/atac_motif_crossregion.R). No object is loaded.
#
#   rows     PRE-SPECIFIED TF families, chosen from the RNA / GRN results before
#            looking at motif significance: AP-1, STAT, CNC/small-MAF, CEBP,
#            NF-kB (Hallmark TNFa/NF-kB and IL6/JAK/STAT3; rule-based GRN TFs
#            STAT1, STAT3, IRF1, MAFG, CEBPD), plus SOX9 (astrocyte identity).
#            One or a few representative JASPAR2020 motifs per family.
#   columns  Astro / Micro / OPC / Endo x INS, PFC, VST (per-region effect from the same
#            fit) and Pooled (average within-region PWH effect)
#   colour   donor-level PWH - PWoH difference in chromVAR deviation z-score,
#            CENTERED model (each sample's mean over all motifs removed)
#   marks    Pooled column only: filled dot = FDR < 0.05 in all three model
#            versions (plain, centered, weighted); open dot = centered only
#
# OUTPUT  OUT_ROOT/04-ATAC_DE/figures/fig2/fig2D_motif_activity.png/.pdf/_caption.txt
#         + fig2D_motif_activity.csv (every value drawn)
#
# Usage: revision/docker/run.sh revision/tools/fig2_motif_heatmap.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

IN  <- file.path(OUT_ROOT, "04-ATAC_DE", "motif_chromvar", "crossregion")
OUT <- file.path(OUT_ROOT, "04-ATAC_DE", "figures", "fig2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
## Endothelium is included because it carries the strongest RNA-matched TF signal
## after Astro/Micro (fig2F candidates: 8 TFs significant on both axes vs 2 for Oligo).
CTS <- c("Astro", "Micro", "OPC", "Endo")
REGIONS <- c("INS", "PFC", "VST")

## family -> JASPAR2020 motif names (as in motif_names.csv) -> display label
ROWS <- tibble::tribble(
  ~family,        ~motif_name,  ~label,
  "AP-1",         "FOSL2",      "FOSL2",
  "AP-1",         "FOS::JUN",   "FOS::JUN",
  "AP-1",         "JUNB",       "JUNB",
  "STAT",         "STAT1",      "STAT1",
  "STAT",         "STAT3",      "STAT3",
  "CNC / MAF",    "BACH1",      "BACH1",
  "CNC / MAF",    "NFE2",       "NFE2",
  "CNC / MAF",    "Nfe2l2",     "NFE2L2",
  "CNC / MAF",    "MAFG",       "MAFG",
  "CEBP",         "CEBPB",      "CEBPB",
  "CEBP",         "CEBPD",      "CEBPD",
  "NF-κB",   "RELA",       "RELA",
  "NF-κB",   "REL",        "REL",
  "NF-κB",   "NFKB1",      "NFKB1",
  "IRF",          "IRF1",       "IRF1",
  "Identity",     "SOX9",       "SOX9")

d <- bind_rows(lapply(CTS, function(ct) {
  x <- utils::read.csv(file.path(IN, paste0(ct, "_chromvar_crossregion.csv")), stringsAsFactors = FALSE)
  miss <- setdiff(ROWS$motif_name, x$name)
  if (length(miss)) stop(ct, ": motifs not found: ", paste(miss, collapse = ", "))
  x <- x[match(ROWS$motif_name, x$name), ]
  robust <- x$FDR < 0.05 & x$FDR_plain < 0.05 & x$FDR_weighted < 0.05
  mark <- ifelse(robust, "all three", ifelse(x$FDR < 0.05, "centered only", NA))
  bind_rows(
    lapply(REGIONS, function(R) data.frame(ROWS, cell_type = ct, column = R,
                                           diff = x[[paste0("diff_", R)]], mark = NA_character_)),
    list(data.frame(ROWS, cell_type = ct, column = "Pooled", diff = x$diff_pooled, mark = mark,
                    FDR_centered = x$FDR, FDR_plain = x$FDR_plain, FDR_weighted = x$FDR_weighted)))
}))
utils::write.csv(d, file.path(OUT, "fig2D_motif_activity.csv"), row.names = FALSE)

d$label  <- factor(d$label, levels = rev(unique(ROWS$label)))
d$family <- factor(d$family, levels = unique(ROWS$family))
d$column <- factor(d$column, levels = c(REGIONS, "Pooled"))
d$cell_type <- factor(d$cell_type, levels = CTS,
                      labels = c("Astrocytes", "Microglia", "OPCs", "Endothelial")[seq_along(CTS)])
lim <- max(abs(d$diff), na.rm = TRUE)

p <- ggplot(d, aes(x = column, y = label, fill = diff)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_point(data = d[!is.na(d$mark), ], aes(shape = mark), size = 1.6, stroke = 0.5,
             colour = "black", fill = "black") +
  scale_shape_manual(values = c(`all three` = 16, `centered only` = 1),
                     name = "Pooled FDR < 0.05") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), name = "PWH − PWoH\nchromVAR z") +
  facet_grid(family ~ cell_type, scales = "free_y", space = "free_y", switch = "y") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_astro_pub(base_size = 7) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 6.5, face = "bold"),
        axis.text.y = element_text(size = 7, face = "italic"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
        strip.text.x = element_text(size = 8),
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.25, "cm"),
        panel.spacing.x = unit(0.5, "lines"), panel.spacing.y = unit(0.2, "lines"))
f <- file.path(OUT, "fig2D_motif_activity.png")
save_pub(p, f, width = 7.2 + 1.1 * (length(CTS) - 3), height = 3.6)
astro_caption(f, paste(
  "Transcription-factor motif activity in PWH vs PWoH. chromVAR deviation z-scores (JASPAR2020 CORE vertebrate",
  "motifs) were averaged per donor within each cell type and region and compared between PWH and PWoH across",
  "donors with limma, adjusting for age, sex, race and PMI; the three regions were modelled together, with the",
  "within-donor correlation across regions estimated by duplicateCorrelation. Each sample's mean deviation over",
  "all motifs was removed first, so global shifts in accessibility are not read as TF activity. Colour: PWH minus",
  "PWoH, per region and pooled (average within-region effect). Rows: TF families pre-specified from the RNA and",
  "gene-regulatory-network results. Pooled column: filled dot, FDR < 0.05 in the centered, uncentered and",
  "sample-weighted (arrayWeights) models; open dot, centered model only."))
message("[fig2D] robust (all three): ", paste(unique(paste(d$cell_type, d$label)[d$mark %in% "all three"]), collapse = ", "))
message("[fig2D] centered only: ", paste(unique(paste(d$cell_type, d$label)[d$mark %in% "centered only"]), collapse = ", "))
message("[fig2D] done")
