# ==============================================================================
# Figure 3 per-region pathway panels with up- and down-regulated terms together.
#
# It does NOT redo the enrichment or the term selection. It reads the rows that
# revision/tools/fig3_pathway_panels.R already chose and drew
#   <run>/fig3_<class>_<up|down>_top<N>_plotted.csv
# and redraws them as ONE horizontal panel per class:
#   x  pathways, in two blocks - up-regulated | gap | down-regulated
#   y  cell types, the three regions stacked as coloured row bands
#   colour  mean logFC on a diverging scale, so the blocks read blue vs red
#   size    -log10 adjusted p
# A vertical orientation (pathways on y, up block above down) is also written.
#
# INPUT   OUT_ROOT/03-DE_analysis/figure_prototypes/fig3_runs/<run>/
# OUTPUT  <run>/updown/fig3_<class>_updown_top<N>_<horizontal|vertical>[_readable].png/.pdf
#         + _gg.rds (the ggplot object, for re-theming at print size)
#
# Usage: revision/docker/run.sh revision/tools/fig3_updown_panels.R [--run <folder under fig3_runs>]
#          [--top 2] [--classes glial,excitatory,inhibitory] [--per-col 0.26] [--wrap 44]
#          [--orient horizontal,vertical] [--readable] [--p-cap 10] [--lfc-cap 1.5]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(cowplot) })

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
RUNS   <- file.path(OUT_ROOT, "03-DE_analysis", "figure_prototypes", "fig3_runs")
## Default: the run using the published DEG rule (|logFC| > 0.25); a run with
## logFC > 0 is selectable with --run.
RUN    <- getarg("--run", "glia-core_simplify-0.6_top2_compact_draw1gene_empty-tested_lfc0.25")
TOP_N  <- as.integer(getarg("--top", "2"))
CLASSES <- trimws(strsplit(getarg("--classes", "glial,excitatory,inhibitory"), ",")[[1]])
PER_COL <- as.numeric(getarg("--per-col", "0.26"))
## Both orientations are written each run. horizontal: pathways on x, up | gap | down blocks side by side,
## regions as row bands. vertical: pathways on y, up block above the down block,
## regions across the top - narrow and tall instead of wide and short.
ORIENTS <- trimws(strsplit(getarg("--orient", "horizontal,vertical"), ",")[[1]])
WRAP_V  <- as.integer(getarg("--wrap", "44"))
## --readable: one outlier (-log10 adj. p ~ 20, cytoplasmic translation in L6 IT Car3)
##   sets the size scale, so every other dot sits at the bottom of it; and the colour
##   scale runs to +/-3 while most terms average +/-0.5-1, so dots are pastel on white.
##   Caps significance at --p-cap (10) and colour at --lfc-cap (1.5; beyond is drawn
##   at the cap, legend says so), enlarges the dot range, and outlines every dot.
##   Files get a _readable suffix; the default panels are unchanged.
READABLE <- "--readable" %in% args
P_CAP    <- as.numeric(getarg("--p-cap", "10"))
LFC_CAP  <- as.numeric(getarg("--lfc-cap", "1.5"))
SFX_R    <- if (READABLE) "_readable" else ""
## Angled labels must stay on ONE line: wrapped two- and three-line labels at 45
## degrees are wide blocks and run into each other. Instead the
## standard GO shorthand of fig3_pathway_panels.R --abbrev keeps them short.
abbreviate_go <- function(x, cap = TRUE) {
  rules <- c(
    "immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains" =
      "immune response (Ig-superfamily receptors)",
    "positive regulation of " = "pos. reg. of ", "negative regulation of " = "neg. reg. of ",
    "regulation of adaptive immune response based on somatic recombination of\u2026" = "adaptive immune response (recombination)",
    "adaptive immune response based on somatic recombination of\u2026" = "adaptive immune response (recombination)",
    "adaptive immune response based on somatic recombination of..." = "adaptive immune response (recombination)",
    "transposable element silencing by heterochromatin formation" = "transposable element silencing (heterochr.)",
    "homophilic cell adhesion via plasma membrane adhesion molecules" = "homophilic cell adhesion",
    "regulation of " = "reg. of ", "biosynthetic process" = "biosynthesis",
    "catabolic process" = "catabolism", "metabolic process" = "metabolism",
    "signaling pathway" = "signaling", "signal transduction" = "signaling")
  for (k in names(rules)) x <- gsub(k, rules[[k]], x, fixed = TRUE)
  ## a hard cap keeps one runaway term from setting the height of the whole figure -
  ## only needed for the angled labels of the horizontal layout; vertical labels are
  ## horizontal text and simply wrap, so they keep the full name
  if (cap) {
    long <- nchar(x) > 46
    x[long] <- paste0(substr(x[long], 1, 45), "\u2026")
  }
  x
}
IN  <- file.path(RUNS, RUN)
OUT <- file.path(IN, "updown")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

REGIONS     <- c("PFC", "INS", "VST")
REGION_FILL <- c(PFC = "#E69F00", INS = "#56B4E9", VST = "#CC79A7")
CT_ORDER <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
              "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT",
              "D2-Striosome","D1/D2-Hybrid","Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg",
              "Sst","Sst Chodl","Vip","Astro","Oligo","OPC","Micro","Macrophages","Other immune",
              "Immune","PC","Endo","SMC","VLMC")
DIR_LAB <- c(up = "Up in PWH", down = "Down in PWH")

read_dir <- function(cls, dir) {
  f <- file.path(IN, sprintf("fig3_%s_%s_top%d_plotted.csv", tolower(cls), dir, TOP_N))
  if (!file.exists(f)) { message("[updown] missing ", basename(f)); return(NULL) }
  utils::read.csv(f, stringsAsFactors = FALSE) %>% mutate(direction = dir)
}

for (cls in CLASSES) {
  d <- bind_rows(lapply(c("up", "down"), read_dir, cls = cls))
  if (is.null(d) || !nrow(d)) next
  d$direction <- factor(DIR_LAB[d$direction], levels = unname(DIR_LAB))
  d$region    <- factor(d$region, levels = REGIONS)
  ct_lv <- intersect(CT_ORDER, unique(d$cell_type))
  ## x order: within each block, the pathway with the strongest dot first, so the
  ## two blocks each read from most to least significant
  ord <- d %>% group_by(direction, label) %>%
    summarise(best = min(p.adjust), .groups = "drop") %>% arrange(direction, best)
  d$label_wrapped <- factor(abbreviate_go(d$label), levels = unique(abbreviate_go(ord$label)))
  lab_chars <- max(nchar(levels(d$label_wrapped)))

  for (ORIENT in ORIENTS) {
  d$label_v <- factor(stringr::str_wrap(abbreviate_go(d$label, cap = FALSE), WRAP_V),
                      levels = unique(stringr::str_wrap(abbreviate_go(ord$label, cap = FALSE), WRAP_V)))
  ## y reads bottom-up, x reads left-right, so the level order flips with orientation
  d$cell_type <- factor(as.character(d$cell_type),
                        levels = if (ORIENT == "horizontal") rev(ct_lv) else ct_lv)
  p <- ggplot(d, aes(x = label_wrapped, y = cell_type,
                     size = -log10(p.adjust), colour = avr_logFC)) +
    geom_point(na.rm = TRUE) +
    scale_size(range = c(0.9, 4.2), name = expression(-log[10]~adj.~p)) +
    ## diverging here, unlike the single-direction panels: with both blocks in one
    ## figure the sign of mean logFC is exactly what the colour should carry
    scale_colour_gradient2(low = "#2166AC", mid = "grey93", high = "#B2182B", midpoint = 0,
                           name = "mean logFC") +
    facet_grid(region ~ direction, scales = "free", space = "free", switch = "y") +
    labs(title = paste0(toupper(substring(cls, 1, 1)), substring(cls, 2)), x = NULL, y = NULL) +
    theme_cowplot(font_size = 9) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank(),
          panel.border = element_rect(colour = "grey78", fill = NA, linewidth = 0.5),
          panel.spacing.x = unit(1.4, "lines"),      # the gap between up and down
          panel.spacing.y = unit(0.35, "lines"),
          strip.text.y.left = element_text(face = "bold", size = 9, colour = "white", angle = 90),
          strip.text.x = element_text(face = "bold", size = 9.5),
          strip.placement = "outside",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
          axis.text.y = element_text(size = 8.5),
          legend.position = "right", legend.title = element_text(size = 8.5, face = "bold"),
          legend.text = element_text(size = 7.5), legend.key.height = unit(0.75, "lines"),
          plot.title = element_text(size = 11, face = "bold"), plot.title.position = "plot",
          plot.margin = margin(6, 10, 4, 26))   # room for the left-most angled label

  if (READABLE) {
    d$size_v <- pmin(-log10(d$p.adjust), P_CAP)
    d$fill_v <- pmax(pmin(d$avr_logFC, LFC_CAP), -LFC_CAP)
    br <- pretty(c(-log10(0.05), P_CAP), n = 3); br <- br[br > 0 & br <= P_CAP]
    p$layers <- list()
    p <- p + geom_point(data = d, aes(size = size_v, fill = fill_v), shape = 21, colour = "grey30",
                        stroke = 0.25, na.rm = TRUE, inherit.aes = TRUE) +
      scale_size(range = c(1.6, 5), limits = c(0, P_CAP), breaks = br,
                 labels = ifelse(br == P_CAP, paste0("\u2265", br), br),
                 name = expression(-log[10]~adj.~p)) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                           limits = c(-LFC_CAP, LFC_CAP),
                           breaks = c(-LFC_CAP, 0, LFC_CAP),
                           labels = c(paste0("\u2264", -LFC_CAP), "0", paste0("\u2265", LFC_CAP)),
                           name = "mean logFC") +
      guides(colour = "none", size = guide_legend(override.aes = list(fill = "grey60")))
  }
  if (ORIENT == "vertical") {
    ## pathways on y (read horizontally), cell types on x, up block above down
    p <- p + aes(x = cell_type, y = label_v) +
      facet_grid(direction ~ region, scales = "free", space = "free", switch = "y") +
      theme(strip.text.y.left = element_text(face = "bold", size = 9, colour = "black", angle = 90),
            strip.text.x = element_text(face = "bold", size = 9.5, colour = "white"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
            axis.text.y = element_text(size = 9, lineheight = 0.85),
            panel.spacing.x = unit(0.35, "lines"),
            panel.spacing.y = unit(1.0, "lines"),      # the gap between up and down
            plot.margin = margin(6, 10, 4, 6))
  }

  ## colour the region strips - left band in horizontal, top band in vertical
  g <- ggplot2::ggplotGrob(p)
  sl <- which(grepl(if (ORIENT == "horizontal") "^strip-l" else "^strip-t", g$layout$name))
  regs <- levels(droplevels(d$region))
  for (i in seq_along(sl)) {
    gr <- g$grobs[[sl[i]]]
    j <- which(vapply(gr$grobs[[1]]$children, function(x) inherits(x, "rect"), logical(1)))
    fill <- if (i <= length(regs)) REGION_FILL[[regs[i]]] else "grey40"
    if (length(j)) gr$grobs[[1]]$children[[j[1]]]$gp$fill <- fill
    g$grobs[[sl[i]]] <- gr
  }

  n_path <- dplyr::n_distinct(d$label_wrapped)
  n_cells <- d %>% distinct(region, cell_type) %>% nrow()
  if (ORIENT == "horizontal") {
    w <- 1.1 + 1.5 + PER_COL * n_path                   # y labels + legend + columns
    ## a 7.5 pt character is ~0.062 in wide; at 45 degrees a label eats 0.7 of its length
    h <- 1.0 + 0.26 * n_cells + 0.7 * 0.062 * lab_chars
  } else {
    ## A discrete axis spaces rows evenly, so one three-line label would inflate
    ## EVERY row. The 90th percentile sets the row height instead.
    v_lines <- lengths(regmatches(levels(d$label_v), gregexpr("\n", levels(d$label_v)))) + 1
    q_lines <- max(1, as.numeric(stats::quantile(v_lines, 0.9)))
    ## The name column (min 2.5 in) follows the longest wrapped line, so a wider
    ## --wrap does not squeeze the panel (9 pt is ~0.064 in per character).
    max_line <- max(nchar(unlist(strsplit(levels(d$label_v), "\n"))))
    w <- max(2.5, 0.064 * max_line + 0.2) + 1.4 + 0.30 * n_cells   # pathway names + legend + cell types
    h <- 1.05 + 0.185 * q_lines * n_path + 0.55
  }
  f <- file.path(OUT, sprintf("fig3_%s_updown_top%d_%s%s.png", tolower(cls), TOP_N, ORIENT, SFX_R))
  ggsave(f, g, width = w, height = h, dpi = 300, limitsize = FALSE, bg = "white")
  ggsave(sub("\\.png$", ".pdf", f), g, width = w, height = h, device = grDevices::cairo_pdf, bg = "white")
  ## the ggplot itself, so a combined figure can re-theme this panel at print
  ## size instead of shrinking a finished image
  saveRDS(list(plot = p, regions = regs, region_fill = REGION_FILL, orient = ORIENT,
               class = cls, top = TOP_N, run = RUN),
          sub("\\.png$", "_gg.rds", f))
  message(sprintf("[updown] %s %s: %d pathways (%d up, %d down), %d region x cell-type slots -> %.1f x %.1f in",
                  cls, ORIENT, n_path, dplyr::n_distinct(d$label_wrapped[d$direction == DIR_LAB[["up"]]]),
                  dplyr::n_distinct(d$label_wrapped[d$direction == DIR_LAB[["down"]]]),
                  n_cells, w, h))
  }
}
message("[updown] done")
