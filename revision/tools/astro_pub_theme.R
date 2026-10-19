# ==============================================================================
# Shared publication styling for the revision figures (function library, sourced
# by the astrocyte and other revision figure tools).
#
#   theme_astro_pub()      one theme for every figure: classic axes, black text,
#                          no grid, bold facet labels without grey boxes, small
#                          legend keys, no subtitle
#   save_pub()             PNG at 600 dpi plus a vector PDF, white background,
#                          sizes in inches at journal widths (7.2 in = 2 columns)
#   astro_caption()        writes the reading instructions to
#                          <figure>_caption.txt beside the PNG, for the legend,
#                          instead of an in-figure subtitle
#   astro_subtype_palette() the Fig1D tSNE colours (ASTRO_SUBTYPE_COLOURS in
#                          revision/tools/astro_reviewer_figures.R) in tSNE order,
#                          grey for anything unexpected
#
# Usage: source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

theme_astro_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text             = element_text(colour = "black"),
      axis.text        = element_text(colour = "black", size = base_size),
      axis.title       = element_text(colour = "black", size = base_size + 1),
      axis.line        = element_line(colour = "black", linewidth = 0.3),
      axis.ticks       = element_line(colour = "black", linewidth = 0.3),
      strip.background = element_blank(),
      strip.text       = element_text(colour = "black", face = "bold", size = base_size + 1),
      legend.title     = element_text(size = base_size + 1),
      legend.text      = element_text(size = base_size),
      legend.key.size  = unit(0.35, "cm"),
      plot.title       = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle    = element_blank(),
      panel.spacing    = unit(0.8, "lines"),
      plot.margin      = margin(4, 6, 4, 4))
}

save_pub <- function(p, png_file, width, height, dpi = 600) {
  dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(png_file, p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  ## cairo_pdf: the base pdf device cannot encode Greek letters (e.g. "TNFα
  ##   signaling via NF-κB" fails with a mbcsToSbcs conversion error)
  ggsave(sub("\\.png$", ".pdf", png_file), p, width = width, height = height,
         device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf",
         bg = "white", limitsize = FALSE)
  message("[pub] ", basename(png_file), " (+ .pdf), ", width, " x ", height, " in")
  invisible(png_file)
}

astro_caption <- function(png_file, text) {
  f <- sub("\\.png$", "_caption.txt", png_file)
  writeLines(strwrap(text, width = 90), f)
  invisible(f)
}

astro_subtype_palette <- function(labs) {
  ref <- if (exists("ASTRO_SUBTYPE_COLOURS")) ASTRO_SUBTYPE_COLOURS else character(0)
  present <- unique(as.character(labs))
  ord <- c(intersect(names(ref), present), sort(setdiff(present, names(ref))))
  setNames(ifelse(ord %in% names(ref), ref[ord], "grey60"), ord)
}
