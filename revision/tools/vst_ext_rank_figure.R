# ==============================================================================
# Figure of cortical excitatory contamination per VST library, with the sample
# exclusion criterion and the published removal list marked.
#
# Method
#   All 33 VST libraries ordered by the percentage of nuclei in cortical
#   excitatory clusters; fill = condition (PWoH / PWH); dashed line = the
#   exclusion criterion (VST_SAMPLE_EXT_THRESHOLD, 20%). Libraries removed in
#   the published analysis are marked with a black outline (bar style) or an
#   open circle (rank style, where the excluded libraries are also labelled).
#   --style bar (default) draws one bar per library; --style rank draws the
#   ranked-points version.
#
# Inputs   OUT_ROOT/01-Regional_data_processing/VST/data/04-VST-tsne_marker_plotting/
#            sample_exclusion_decision.csv   ext_fraction, excluded, published_removal
#          condition from the sample sheet (SAMPLE_META, revision/config/paths.R)
# Outputs  OUT_ROOT/01-Regional_data_processing/VST/figures/vst_ext/
#            vst_excitatory_by_library[_rank].png/.pdf + _caption.txt
#
# Usage
#   revision/docker/run.sh revision/tools/vst_ext_rank_figure.R [--style bar|rank]
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))

VD  <- file.path(OUT_ROOT, "01-Regional_data_processing", "VST")
OUT <- file.path(VD, "figures", "vst_ext")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- utils::read.csv(file.path(VD, "data", "04-VST-tsne_marker_plotting", "sample_exclusion_decision.csv"),
                     stringsAsFactors = FALSE)
m <- utils::read.csv(SAMPLE_META, stringsAsFactors = FALSE, colClasses = "character")
m <- m[!duplicated(m$sample_ids), c("sample_ids", "hiv_jay_v2")]
d$condition <- factor(ifelse(toupper(m$hiv_jay_v2[match(d$sample_ids, m$sample_ids)]) == "TRUE", "PWH", "PWoH"),
                      levels = PWH_LEVELS)
if (anyNA(d$condition)) stop("no condition for: ", paste(d$sample_ids[is.na(d$condition)], collapse = ", "))
d <- d[order(-d$ext_fraction), ]
d$rank <- seq_len(nrow(d))
d$pct  <- 100 * d$ext_fraction
d$excluded <- toupper(d$excluded) == "TRUE"
d$published_removal <- toupper(d$published_removal) == "TRUE"
message(sprintf("[rank] %d libraries; excluded now %d; published list %d; overlap %d",
                nrow(d), sum(d$excluded), sum(d$published_removal), sum(d$excluded & d$published_removal)))

## Default: one bar per library, ordered; the sample ids go on the axis so no labels
##   are needed, and the published removals are marked by a black outline.
##   --style rank draws the ranked-points version.
STYLE <- if ("--style" %in% commandArgs(TRUE))
  commandArgs(TRUE)[match("--style", commandArgs(TRUE)) + 1] else "bar"
d$sample_ids <- factor(d$sample_ids, levels = d$sample_ids[order(-d$pct)])
thr_lab <- sprintf("exclusion criterion, %.0f%%", 100 * VST_SAMPLE_EXT_THRESHOLD)

if (STYLE == "bar") {
  p <- ggplot(d, aes(sample_ids, pct, fill = condition)) +
    geom_col(aes(colour = published_removal), width = 0.75, linewidth = 0.35) +
    geom_hline(yintercept = 100 * VST_SAMPLE_EXT_THRESHOLD, linetype = "dashed",
               colour = "grey35", linewidth = 0.3) +
    annotate("text", x = nrow(d) - 0.2, y = 100 * VST_SAMPLE_EXT_THRESHOLD + 1.8, hjust = 1,
             size = 2.5, label = thr_lab) +
    scale_fill_manual(values = PWH_COLOURS, name = NULL) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = NA), na.value = NA,
                        breaks = "TRUE", labels = "removed in the published analysis", name = NULL) +
    scale_y_continuous(limits = c(0, max(d$pct) * 1.08), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "nuclei in cortical excitatory clusters (%)") +
    theme_astro_pub(8) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
          legend.position = "top", legend.margin = margin(b = -4),
          legend.spacing.x = unit(0.2, "cm")) +
    guides(colour = guide_legend(override.aes = list(fill = NA, linewidth = 0.4)))
} else {
  p <- ggplot(d, aes(rank, pct)) +
    geom_line(colour = "grey75", linewidth = 0.3) +
    geom_hline(yintercept = 100 * VST_SAMPLE_EXT_THRESHOLD, linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    geom_point(data = d[d$published_removal, ], shape = 21, size = 3.4, stroke = 0.5,
               colour = "grey25", fill = NA) +
    geom_point(aes(fill = condition), shape = 21, size = 2.1, stroke = 0.2, colour = "white") +
    ggrepel::geom_text_repel(data = d[d$excluded, ], aes(label = sample_ids), size = 2.4,
                             nudge_x = 1.6, direction = "y", hjust = 0, segment.size = 0.2,
                             segment.colour = "grey60", seed = 1) +
    annotate("text", x = nrow(d), y = 100 * VST_SAMPLE_EXT_THRESHOLD + 1.6, hjust = 1, size = 2.5,
             label = thr_lab) +
    scale_fill_manual(values = PWH_COLOURS, name = NULL) +
    scale_x_continuous(breaks = c(1, seq(5, nrow(d), 5))) +
    scale_y_continuous(limits = c(0, max(d$pct) * 1.08), expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = "VST library, ranked", y = "nuclei in cortical excitatory clusters (%)") +
    theme_astro_pub(8) + theme(legend.position = c(0.88, 0.78), legend.background = element_blank())
}

f <- file.path(OUT, sprintf("vst_excitatory_by_library%s.png", if (STYLE == "bar") "" else "_rank"))
W <- if (STYLE == "bar") 5.6 else 4.6
ggsave(f, p, width = W, height = 3.4, dpi = 600, bg = "white")
ggsave(sub("\\.png$", ".pdf", f), p, width = W, height = 3.4, device = grDevices::cairo_pdf, bg = "white")
writeLines(strwrap(sprintf(paste(
  "Cortical excitatory contamination of the VST libraries. Each bar is one library (n = %d), ordered by the",
  "percentage of its nuclei assigned to cortical excitatory clusters and coloured by condition. The dashed line is",
  "the exclusion criterion used here (%.0f%%); the %d libraries above it were excluded. A black outline marks",
  "the %d libraries removed in the original analysis: %d of them fall below the criterion and are retained here,",
  "while donor 1040 (%.0f%%) was not on that list and is now excluded. The distribution is a step rather than a",
  "gradient - excluded libraries contain %.0f-%.0f%% excitatory nuclei and every retained library at most %.1f%% -",
  "which indicates variable inclusion of adjacent cortex during dissection rather than sampling of bordering tissue."),
  nrow(d), 100 * VST_SAMPLE_EXT_THRESHOLD, sum(d$excluded), sum(d$published_removal),
  sum(d$published_removal & !d$excluded), d$pct[d$sample_ids == "1040"],
  min(d$pct[d$excluded]), max(d$pct[d$excluded]), max(d$pct[!d$excluded])), 110),
  sub("\\.png$", "_caption.txt", f))
message("[rank] wrote ", f)
