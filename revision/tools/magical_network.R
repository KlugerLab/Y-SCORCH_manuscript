# ==============================================================================
# MAGICAL regulatory network tables and figures for one region x cell type.
#
# Method: reproduces the tables of
# 06-astrocyte_analysis/gene-regulatory-networks/analyze_magical_output.Rmd
# (the "filter the TF first" path used for the published figure): expand the
# MAGICAL circuits, keep the chosen TFs, then (--genes up) keep circuits whose
# target is an up-regulated DE gene (|logFC| >= 0.5, p_adj.glb <= 0.05, the
# MAGICAL candidate cutoffs). The igraph figure reproduces the published
# plotting script (layout_with_gem, seed 2; gold circles = up-regulated genes,
# beige squares = TFs, red / blue rhombi = DA peaks up / down). ggraph versions
# add a TF -> gene view with peaks folded into the edges, a compact page-size
# version and a "pill" version with gene names as the nodes.
#
# Differences from the Rmd / plotting script, needed for the current data:
#   - DE inputs are the flat *_withFreq tables, not pbDS objects with $table;
#     filtered on cluster_id with the same cutoffs.
#   - current MAGICAL writes the TF column as "TFs(binding prob)" =
#     "JUNB,(0.953),NRF1,(0.941),...". A plain comma split would turn every
#     probability into a "TF", so it is parsed as TF/probability pairs;
#     binding_prob and Looping_prob are kept as extra columns.
#   - a TF that is also a target gene (e.g. PFC: STAT1) is one igraph vertex;
#     node_types keeps it as a TF, avoiding a duplicate row name.
#   - the output folder is created.
#
# Inputs: OUT_ROOT/03-DE_analysis/<R>/results/<R>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#         OUT_ROOT/04-ATAC_DE/<R>/data_ATAC/<R>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#         MAGICAL_OUTPUT_DIR/[variants/<variant>/]<R>_<CT>MAGICAL_selected_regulatory_circuits.txt
# Output: MAGICAL_ANALYZED/<R>_new_<CT>[_<variant>][/<tag>]/ - expanded circuits,
#         Cytoscape edge/node files, tf_summary.csv, tf_gene_edges.csv, figures
#
# Usage: revision/docker/run.sh revision/tools/magical_network.R --region INS [--celltype Astro]
#          [--tfs STAT3,STAT1,FOSL2,FOS,JUNB,IRF1,KLF4] [--genes up|all]
#          [--style both|igraph|ggraph] [--variant <v>] [--tag <t>] [--highlight G1,G2]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(readr); library(igraph)
})

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
region   <- getarg("--region", NA_character_)
celltype <- getarg("--celltype", "Astro")
## the published PFC astrocyte TF set (analyze_magical_output.Rmd)
chosen_TF <- trimws(strsplit(getarg("--tfs", "STAT3,STAT1,FOSL2,FOS,JUNB,IRF1,KLF4"), ",")[[1]])
GENES    <- match.arg(getarg("--genes", "up"), c("up", "all"))
SEED     <- as.integer(getarg("--seed", "2"))
## igraph = the published plotting script's look; ggraph = the ggplot2 versions
STYLE    <- match.arg(getarg("--style", "both"), c("both", "igraph", "ggraph"))
if (is.na(region)) stop("--region is required")

folder <- file.path(MAGICAL_ANALYZED, magical_tag(region, celltype))
## --tag: an alternative TF set writes to its own subfolder, so it never
##   overwrites the published-TF outputs or their tf_gene_edges.csv.
## --variant (e.g. the expression-filtered MAGICAL run) reads
##   output/variants/<variant>/ and writes analyzed_output/<folder>_<variant>/,
##   so the primary figures and tables are never overwritten.
VARIANT <- getarg("--variant", "")
if (nzchar(VARIANT)) folder <- paste0(folder, "_", VARIANT)
TAG <- getarg("--tag", "")
if (nzchar(TAG)) folder <- file.path(folder, TAG)
dir.create(folder, recursive = TRUE, showWarnings = FALSE)
message("[net] ", region, " ", celltype, " -> ", folder)

## ---- DE (flat withFreq tables) -------------------------------------------------
read_de <- function(f) {
  x <- as.data.frame(readRDS(f))
  x[as.character(x$cluster_id) == celltype, , drop = FALSE]
}
atac_de_micro <- read_de(file.path(OUT_ROOT, "04-ATAC_DE", region, "data_ATAC",
                                   paste0(region, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS")))
atac_de_micro_up   <- atac_de_micro[which(atac_de_micro$logFC >= 0.5 & atac_de_micro$p_adj.glb <= 0.05), ]
atac_de_micro_down <- atac_de_micro[which(atac_de_micro$logFC <= -0.5 & atac_de_micro$p_adj.glb <= 0.05), ]

rna_de_micro <- read_de(file.path(OUT_ROOT, "03-DE_analysis", region, "results",
                                  paste0(region, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS")))
rna_de_micro_up   <- rna_de_micro[which(rna_de_micro$logFC >= 0.5 & rna_de_micro$p_adj.glb <= 0.05), ]
rna_de_micro_down <- rna_de_micro[which(rna_de_micro$logFC <= -0.5 & rna_de_micro$p_adj.glb <= 0.05), ]

## ---- circuits (TF / binding-probability pairs) --------------------------------
circuits_file <- file.path(MAGICAL_OUTPUT_DIR, if (nzchar(VARIANT)) file.path("variants", VARIANT) else "",
                           paste0(region, "_", gsub("[/ ]", "-", celltype),
                                  "MAGICAL_selected_regulatory_circuits.txt"))
message("[net] circuits: ", circuits_file)
circuits <- read_tsv(circuits_file, col_types = cols(.default = "c"), show_col_types = FALSE)
tf_col <- grep("^TFs", names(circuits), value = TRUE)[1]
expanded_circuits <- circuits %>%
  mutate(.id = row_number(), tok = str_split(.data[[tf_col]], ",")) %>%
  unnest(tok) %>% mutate(tok = str_trim(tok)) %>% filter(tok != "") %>%
  group_by(.id) %>%
  mutate(is_prob = str_detect(tok, "^\\("), pair = cumsum(!is_prob)) %>%
  group_by(.id, pair) %>%
  summarise(TFs = tok[!is_prob][1],
            binding_prob = suppressWarnings(as.numeric(str_remove_all(tok[is_prob][1], "[()]"))),
            Gene_symbol = first(Gene_symbol), Peak_chr = first(Peak_chr),
            Peak_start = first(Peak_start), Peak_end = first(Peak_end),
            Looping_prob = as.numeric(first(Looping_prob)), .groups = "drop") %>%
  filter(!is.na(TFs), TFs != "") %>%
  mutate(peakname = paste0(Peak_chr, "-", Peak_start, "-", Peak_end)) %>%
  select(TFs, peakname, Gene_symbol, binding_prob, Looping_prob) %>%
  arrange(TFs)
stopifnot(!any(str_detect(expanded_circuits$TFs, "^\\(")))
write.csv(expanded_circuits, file.path(folder, "expanded_circuits.csv"), row.names = FALSE)

cnt <- function(d, lab) message(sprintf("[net] %-28s TFs %4d | peaks %4d | genes %4d", lab,
                                        n_distinct(d$TFs), n_distinct(d$peakname), n_distinct(d$Gene_symbol)))
cnt(expanded_circuits, "all circuits")

## ---- filter the TF first (the published path) ---------------------------------
expanded_circuits_filtered <- expanded_circuits[which(expanded_circuits$TFs %in% chosen_TF), ]
cnt(expanded_circuits_filtered, "after TF filter")
missing_tf <- setdiff(chosen_TF, expanded_circuits_filtered$TFs)
if (length(missing_tf)) message("[net] chosen TFs absent from this network: ", paste(missing_tf, collapse = ", "))
if (GENES == "up") {
  expanded_circuits_filtered <- expanded_circuits_filtered[which(expanded_circuits_filtered$Gene_symbol %in% rna_de_micro_up$gene), ]
  cnt(expanded_circuits_filtered, "after up-regulated genes")
}
if (!nrow(expanded_circuits_filtered)) stop("nothing left to draw for ", region, " ", celltype)
write.csv(expanded_circuits_filtered, file.path(folder, "expanded_circuits_filtered.csv"), row.names = FALSE)

tf_summary <- expanded_circuits_filtered %>% group_by(TFs) %>%
  summarise(peaks = n_distinct(peakname), genes = n_distinct(Gene_symbol), .groups = "drop") %>%
  arrange(desc(genes))
write.csv(tf_summary, file.path(folder, "tf_summary.csv"), row.names = FALSE)
print(as.data.frame(tf_summary), row.names = FALSE)

## ---- network files (as the Rmd) ------------------------------------------------
edges <- function(data) rbind(
  data %>% dplyr::select(Source = TFs, Target = peakname) %>% mutate(Interaction = "regulates"),
  data %>% dplyr::select(Source = peakname, Target = Gene_symbol) %>% mutate(Interaction = "associated_with"))
network <- edges(expanded_circuits_filtered)
write.table(network, file.path(folder, "cyto_df.csv"), sep = "\t", row.names = FALSE, quote = FALSE)
network <- unique(network)
write.table(network, file.path(folder, "cyto_df_no_overlap.csv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(network[, c("Source", "Interaction", "Target")], file.path(folder, "cyto_df_no_overlap.csv.sif"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## ---- node types (as the Rmd; a TF wins over its own gene) ---------------------
f <- expanded_circuits_filtered
node_types <- bind_rows(
  tibble(node = unique(f$TFs), type = "TFs"),
  tibble(node = unique(f$peakname[f$peakname %in% atac_de_micro_up$gene]), type = "peak_up"),
  tibble(node = unique(f$peakname[f$peakname %in% atac_de_micro_down$gene]), type = "peak_down"),
  tibble(node = unique(f$Gene_symbol[f$Gene_symbol %in% rna_de_micro_up$gene]), type = "de_up"),
  tibble(node = unique(f$Gene_symbol[f$Gene_symbol %in% rna_de_micro_down$gene]), type = "de_down"))
dup <- unique(node_types$node[duplicated(node_types$node)])
if (length(dup)) message("[net] TF and target gene at once (one vertex, kept as TF): ", paste(dup, collapse = ", "))
node_types <- node_types[!duplicated(node_types$node), ]
untyped <- setdiff(unique(c(network$Source, network$Target)), node_types$node)
if (length(untyped)) {
  message("[net] ", length(untyped), " node(s) with no DE/DA type (drawn grey): ", paste(head(untyped, 8), collapse = ", "))
  node_types <- bind_rows(node_types, tibble(node = untyped,
                          type = ifelse(grepl("^chr", untyped), "peak_other", "gene_other")))
}
write.table(node_types, file.path(folder, "node_types.csv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(table(node_types$type))

## ---- figure: the published plotting script -------------------------------------
myrhombus <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) vertex.color <- vertex.color[v]
  ## igraph >= 2 hands custom shapes the size already scaled by 1/200
  ## (vertex.size = 10 arrives as 0.05, and igraph's own circle uses it
  ## directly). The published plotting script was written for igraph 1.x and
  ## divided again, which draws every peak rhombus ~200x too small.
  vertex.size <- params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) vertex.size <- vertex.size[v]
  symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color,
          stars = cbind(1.2 * vertex.size, vertex.size, 1.2 * vertex.size, vertex.size),
          add = TRUE, inches = FALSE)
}
add_shape("rhombus", clip = shapes("circle")$clip, plot = myrhombus)

shape_mapper <- c(de_up = "circle", TFs = "square", peak_up = "rhombus", peak_down = "rhombus",
                  de_down = "circle", peak_other = "rhombus", gene_other = "circle")
col_mapper   <- c(peak_down = "#78B7C5", peak_up = "#F2300F", de_up = "gold1", TFs = "beige",
                  de_down = "#78B7C5", peak_other = "grey80", gene_other = "grey80")

ig <- graph_from_data_frame(network[, c("Source", "Target")], directed = FALSE)
nt <- setNames(node_types$type, node_types$node)
V(ig)$node_type <- unname(nt[V(ig)$name])
V(ig)$labels <- ifelse(V(ig)$node_type %in% c("de_up"), V(ig)$name, "")
V(ig)$shape  <- unname(shape_mapper[V(ig)$node_type])
V(ig)$color  <- unname(col_mapper[V(ig)$node_type])
node.size_mapper <- c(de_up = 17, TFs = 10, peak_up = 10, peak_down = 10, de_down = 17,
                      peak_other = 10, gene_other = 17)
node.sizes <- unname(node.size_mapper[V(ig)$node_type])

set.seed(SEED)
layout <- layout_with_gem(ig)

draw <- function(file, label, dev = c("png", "svg")) {
  dev <- match.arg(dev)
  if (dev == "png") png(file, height = 3500, width = 3500, res = 300) else svg(file, height = 10, width = 10)
  plot(ig, layout = layout, vertex.label.dist = 0,
       vertex.size = if (dev == "png" && label) node.sizes else 10,
       vertex.label = if (label) V(ig)$labels else rep("", vcount(ig)),
       vertex.shape = V(ig)$shape,
       vertex.label.cex = if (dev == "png" && label) 1.2 else 1,
       vertex.label.color = "black",
       edge.width = if (dev == "png" && label) 2 else 1)
  invisible(dev.off())
}
if (STYLE %in% c("igraph", "both")) {
  draw(file.path(folder, "igraph_layout_label.png"), TRUE, "png")
  draw(file.path(folder, "igraph_layout_nolabel.png"), FALSE, "png")
  if (isTRUE(capabilities("cairo"))) {
    draw(file.path(folder, "igraph_layout_label.svg"), TRUE, "svg")
    draw(file.path(folder, "igraph_layout_nolabel.svg"), FALSE, "svg")
  } else message("[net] no cairo - SVGs skipped")
}

## ---- figure 2: ggraph -----------------------------------------------------------
## Same encoding as the published figure (colours, shapes: circle = gene,
## square = TF, diamond = peak), but ggplot2: TFs labelled, repelled labels, a
## legend, and a PDF for editing. Two layouts to choose between:
##   stress   - force-directed, closest to the published look
##   layered  - TF -> peak -> gene, left to right (sugiyama)
if (STYLE %in% c("ggraph", "both")) {
  suppressPackageStartupMessages({ library(ggraph); library(tidygraph); library(ggplot2) })
  type_lab   <- c(TFs = "Transcription factor", peak_up = "Peak, more accessible in PWH",
                  peak_down = "Peak, less accessible in PWH", de_up = "Gene, up in PWH",
                  de_down = "Gene, down in PWH", peak_other = "Peak (not DA)", gene_other = "Gene (not DE)")
  type_shape <- c(TFs = 22, peak_up = 23, peak_down = 23, de_up = 21, de_down = 21,
                  peak_other = 23, gene_other = 21)
  type_size  <- c(TFs = 7, peak_up = 3, peak_down = 3, de_up = 4.5, de_down = 4.5,
                  peak_other = 3, gene_other = 4.5)
  present <- intersect(names(type_lab), unique(V(ig)$node_type))
  is_gene <- V(ig)$node_type %in% c("de_up", "de_down", "gene_other")
  layer   <- ifelse(V(ig)$node_type == "TFs", 1, ifelse(grepl("^peak", V(ig)$node_type), 2, 3))

  tg <- as_tbl_graph(ig) %>% activate(nodes) %>%
    mutate(type     = factor(node_type, levels = present),
           gene_lab = ifelse(node_type %in% c("de_up", "de_down", "gene_other"), name, NA_character_),
           tf_lab   = ifelse(node_type == "TFs", name, NA_character_))

  scales_common <- list(
    scale_shape_manual(values = type_shape[present], labels = type_lab[present], name = NULL, drop = FALSE),
    scale_fill_manual(values = col_mapper[present], labels = type_lab[present], name = NULL, drop = FALSE),
    scale_size_manual(values = type_size[present], labels = type_lab[present], name = NULL, drop = FALSE),
    labs(title = paste(region, celltype)),
    theme_void(base_size = 12),
    theme(legend.position = "bottom", legend.text = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 14),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(12, 12, 12, 12)))
  save_gg <- function(p, stub, w, h) {
    ggsave(file.path(folder, paste0(stub, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
    ggsave(file.path(folder, paste0(stub, ".pdf")), p, width = w, height = h,
           device = if (isTRUE(capabilities("cairo"))) cairo_pdf else "pdf")
  }

  ## stress (force-directed). TFs are drawn as their label box, with no square
  ## point under the name; the legend keeps the TF key through drop = FALSE.
  ## Legend on two rows - on one row it runs off both sides of the INS figure.
  set.seed(SEED)
  side <- max(6, min(14, 4 + 0.08 * vcount(ig)))
  p1 <- ggraph(tg, layout = "stress") +
    geom_edge_link(colour = "grey72", edge_width = 0.35, edge_alpha = 0.9) +
    geom_node_point(aes(filter = node_type != "TFs", shape = type, fill = type, size = type),
                    colour = "grey15", stroke = 0.4) +
    geom_node_text(aes(label = gene_lab), size = 3, repel = TRUE, max.overlaps = Inf,
                   segment.colour = "grey60", segment.size = 0.2, na.rm = TRUE) +
    geom_node_label(aes(filter = node_type == "TFs", label = name), size = 3.8,
                    fontface = "bold", fill = col_mapper[["TFs"]], colour = "black",
                    label.r = unit(0.08, "lines"), label.padding = unit(0.35, "lines")) +
    scales_common +
    guides(shape = guide_legend(nrow = 2, byrow = TRUE),
           fill  = guide_legend(nrow = 2, byrow = TRUE),
           size  = guide_legend(nrow = 2, byrow = TRUE))
  save_gg(p1, "network_ggraph_stress", side, side)

  ## layered TF -> peak -> gene, left to right
  n_rows <- max(table(layer))
  p2 <- ggraph(tg, layout = "sugiyama", layers = layer) +
    geom_edge_link(colour = "grey72", edge_width = 0.35, edge_alpha = 0.9) +
    geom_node_point(aes(shape = type, fill = type, size = type), colour = "grey15", stroke = 0.4) +
    geom_node_text(aes(label = gene_lab), size = 3, hjust = 0, nudge_y = 0.12, na.rm = TRUE) +
    geom_node_text(aes(label = tf_lab), size = 3.8, fontface = "bold", hjust = 1,
                   nudge_y = -0.12, na.rm = TRUE) +
    coord_flip(clip = "off") +
    scale_y_reverse(expand = expansion(mult = 0.18)) +
    scales_common
  save_gg(p2, "network_ggraph_layered", 8.5, max(5, min(20, 1.5 + 0.2 * n_rows)))
  ## TF -> gene with the peaks folded into the edges.
  ## One edge per TF-peak-gene circuit, so a TF reaching a gene through two
  ## peaks draws two fanned lines. Colour = accessibility of the linked peak in
  ## PWH; width = MAGICAL's TF-peak binding probability.
  e2 <- expanded_circuits_filtered %>%
    group_by(TFs, peakname, Gene_symbol) %>%
    summarise(binding_prob = max(binding_prob), Looping_prob = max(Looping_prob), .groups = "drop") %>%
    mutate(peak_dir = case_when(peakname %in% atac_de_micro_up$gene   ~ "More accessible in PWH",
                                peakname %in% atac_de_micro_down$gene ~ "Less accessible in PWH",
                                TRUE ~ "Not differentially accessible")) %>%
    transmute(from = TFs, to = Gene_symbol, peak_dir, binding_prob, Looping_prob, peakname)
  n2 <- tibble(name = unique(c(e2$from, e2$to))) %>%
    mutate(node_type = ifelse(name %in% e2$from, "TFs", unname(nt[name])),
           node_type = ifelse(is.na(node_type), "gene_other", node_type))
  write.csv(e2, file.path(folder, "tf_gene_edges.csv"), row.names = FALSE)
  tg2 <- tbl_graph(nodes = n2, edges = e2, directed = TRUE)
  dir_cols <- c("More accessible in PWH" = col_mapper[["peak_up"]],
                "Less accessible in PWH" = col_mapper[["peak_down"]],
                "Not differentially accessible" = "grey70")
  gene_types <- intersect(c("de_up", "de_down", "gene_other"), unique(n2$node_type))

  ## Binding probabilities only span ~0.8-1.0, so the width range is thin
  ## (0.15-1.5) and edges semi-transparent; gene names get a white halo
  ## (ggrepel bg.colour) and extra repulsion for legibility.
  set.seed(SEED)
  side2 <- max(7, min(16, 5 + 0.1 * nrow(n2)))
  p3 <- ggraph(tg2, layout = "stress") +
    geom_edge_fan(aes(colour = peak_dir, width = binding_prob), alpha = 0.6, strength = 1.3) +
    geom_edge_loop(aes(colour = peak_dir, width = binding_prob), alpha = 0.6) +
    scale_edge_colour_manual(values = dir_cols, name = "Linked peak") +
    scale_edge_width(range = c(0.15, 1.5), name = "TF binding probability") +
    geom_node_point(aes(filter = node_type != "TFs", fill = node_type), shape = 21, size = 3.2,
                    colour = "grey15", stroke = 0.35) +
    scale_fill_manual(values = col_mapper[gene_types], labels = type_lab[gene_types], name = NULL) +
    geom_node_text(aes(filter = node_type != "TFs", label = name), size = 3.4, repel = TRUE,
                   max.overlaps = Inf, bg.colour = "white", bg.r = 0.15,
                   box.padding = 0.45, point.padding = 0.25, min.segment.length = 0.3,
                   segment.colour = "grey55", segment.size = 0.25, seed = SEED) +
    geom_node_label(aes(filter = node_type == "TFs", label = name), size = 4, fontface = "bold",
                    fill = col_mapper[["TFs"]], colour = "black",
                    label.r = unit(0.08, "lines"), label.padding = unit(0.35, "lines")) +
    labs(title = paste(region, celltype)) +
    theme_void(base_size = 12) +
    theme(legend.position = "right", legend.text = element_text(size = 10),
          legend.title = element_text(size = 10, face = "bold"),
          plot.title = element_text(face = "bold", size = 14),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(12, 12, 12, 12))
  save_gg(p3, "network_ggraph_tf_gene", side2 + 2.5, side2)

  ## Compact publication version. The tf_gene figure above is drawn 12-16 in
  ##   wide with ~10 pt labels, so at page width the names shrink to ~4 pt. This
  ##   one is drawn at page size (--width, default 7.2 in = two columns), so the
  ##   sizes are the printed sizes: gene names ~7.5 pt, TF names 9 pt bold,
  ##   edges 0.12-0.6, small gene dots, legend in one block underneath.
  W <- as.numeric(getarg("--width", "7.2")); H <- as.numeric(getarg("--height", "6.8"))
  set.seed(SEED)
  p4 <- ggraph(tg2, layout = "stress") +
    geom_edge_fan(aes(colour = peak_dir, width = binding_prob), alpha = 0.5, strength = 1.2) +
    geom_edge_loop(aes(colour = peak_dir, width = binding_prob), alpha = 0.5) +
    scale_edge_colour_manual(values = dir_cols, name = "Linked peak") +
    scale_edge_width(range = c(0.12, 0.6), name = "TF binding prob.") +
    geom_node_point(aes(filter = node_type != "TFs", fill = node_type), shape = 21, size = 1.3,
                    colour = "grey20", stroke = 0.2) +
    scale_fill_manual(values = col_mapper[gene_types], labels = type_lab[gene_types], name = NULL) +
    geom_node_text(aes(filter = node_type != "TFs", label = name), size = 2.65, repel = TRUE,
                   max.overlaps = Inf, bg.colour = "white", bg.r = 0.12,
                   box.padding = 0.18, point.padding = 0.1, min.segment.length = 0.25, force = 1.5,
                   segment.colour = "grey60", segment.size = 0.15, seed = SEED) +
    geom_node_label(aes(filter = node_type == "TFs", label = name), size = 3.2, fontface = "bold",
                    fill = col_mapper[["TFs"]], colour = "black",
                    label.r = unit(0.06, "lines"), label.padding = unit(0.22, "lines")) +
    theme_void(base_size = 8) +
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.text = element_text(size = 7), legend.title = element_text(size = 7, face = "bold"),
          legend.key.size = unit(0.3, "cm"), legend.spacing.x = unit(0.2, "cm"),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(4, 4, 4, 4)) +
    guides(edge_colour = guide_legend(order = 1, nrow = 2, override.aes = list(edge_width = 0.8)),
           edge_width  = guide_legend(order = 2, nrow = 2),
           fill        = guide_legend(order = 3))
  save_gg(p4, "network_tf_gene_compact", W, H)
  message("[net] compact figure: ", W, " x ", H, " in, TFs: ", paste(chosen_TF, collapse = ", "))

  ## Pill version: each gene is its label - the name sits inside a rounded gold
  ##   pill (8 pt), each TF inside its box (9.5 pt bold), and the edges run
  ##   underneath. Because the labels are the nodes they must not overlap, so the
  ##   stress layout is scaled to inches and nodes are pushed apart until no two
  ##   label boxes intersect. The canvas
  ##   starts at --pill-width x --pill-height and grows only as far as the
  ##   de-overlapped layout needs, so the text stays at the printed size.
  PILL_PT <- as.numeric(getarg("--pill-pt", "8")); TF_PT <- PILL_PT + 1.5
  PW <- as.numeric(getarg("--pill-width", "5.6")); PH <- as.numeric(getarg("--pill-height", "5.2"))
  ## --highlight SERPINA3,OSMR,CHI3L1: these gene names in bold, for the
  ##   reactive-marker version. Empty by default = no highlighting.
  HL <- trimws(strsplit(getarg("--highlight", ""), ",")[[1]])
  LEG_IN <- 0.75                                     # room for the legend underneath
  lay <- { set.seed(SEED); create_layout(tg2, layout = "stress") }
  is_tf <- lay$node_type == "TFs"
  lay$face <- ifelse(is_tf | lay$name %in% HL, "bold", "plain")
  grDevices::pdf(NULL)                               # a device, so text can be measured
  txt_w <- vapply(seq_len(nrow(lay)), function(i)
    grid::convertWidth(grid::grobWidth(grid::textGrob(lay$name[i], gp = grid::gpar(
      fontsize = if (is_tf[i]) TF_PT else PILL_PT, fontface = lay$face[i]))),
      "inches", valueOnly = TRUE), numeric(1))
  grDevices::dev.off()
  ## label padding + a small gap; TF boxes are padded more so adjacent TF
  ## borders do not touch
  box_w <- txt_w + ifelse(is_tf, 0.2, 0.12)
  box_h <- ifelse(is_tf, TF_PT, PILL_PT) * 1.25 / 72 + ifelse(is_tf, 0.13, 0.09)
  ## scale the layout into the starting panel, in inches
  sx <- diff(range(lay$x)); sy <- diff(range(lay$y))
  sc <- min((PW - max(box_w) - 0.3) / sx, (PH - LEG_IN - max(box_h) - 0.3) / sy)
  x <- (lay$x - min(lay$x)) * sc; y <- (lay$y - min(lay$y)) * sc
  ## Disconnected pieces (INS: MAFG -> HERC6 on its own) are placed by stress far
  ## from the main network, leaving empty space. Keep the largest piece where it
  ## is and pack the others in rows directly underneath it, largest first.
  memb <- igraph::components(igraph::as.undirected(igraph::as.igraph(tg2)))$membership
  if (length(unique(memb)) > 1) {
    sizes <- sort(table(memb), decreasing = TRUE)
    main  <- memb == as.integer(names(sizes)[1])
    left  <- min(x[main] - box_w[main] / 2); right <- max(x[main] + box_w[main] / 2)
    cursor_x <- left; row_top <- min(y[main] - box_h[main] / 2) - 0.25; row_h <- 0
    for (k in as.integer(names(sizes)[-1])) {
      m <- memb == k
      cw <- diff(range(x[m])) + max(box_w[m]); ch <- diff(range(y[m])) + max(box_h[m])
      if (cursor_x > left && cursor_x + cw > right) {           # start a new row
        cursor_x <- left; row_top <- row_top - row_h - 0.2; row_h <- 0 }
      x[m] <- x[m] - min(x[m] - box_w[m] / 2) + cursor_x
      y[m] <- y[m] - max(y[m] + box_h[m] / 2) + row_top
      cursor_x <- cursor_x + cw + 0.25; row_h <- max(row_h, ch)
    }
    message("[net] pill layout: ", length(sizes) - 1, " disconnected piece(s) packed under the main network")
  }
  ## push overlapping label boxes apart along the axis of least overlap
  for (it in seq_len(400)) {
    moved <- FALSE
    for (i in seq_len(nrow(lay) - 1)) for (j in (i + 1):nrow(lay)) {
      dx <- x[j] - x[i]; dy <- y[j] - y[i]
      ox <- (box_w[i] + box_w[j]) / 2 - abs(dx); oy <- (box_h[i] + box_h[j]) / 2 - abs(dy)
      if (ox > 0 && oy > 0) {
        moved <- TRUE
        if (ox / (box_w[i] + box_w[j]) < oy / (box_h[i] + box_h[j])) {
          s <- if (dx == 0) 1 else sign(dx); x[i] <- x[i] - s * ox / 2; x[j] <- x[j] + s * ox / 2
        } else {
          s <- if (dy == 0) 1 else sign(dy); y[i] <- y[i] - s * oy / 2; y[j] <- y[j] + s * oy / 2
        }
      }
    }
    if (!moved) break
  }
  message("[net] pill layout: overlaps resolved after ", it, " passes", if (moved) " (NOT fully resolved)" else "")
  lay$x <- x; lay$y <- y
  pad_x <- max(box_w) / 2 + 0.05; pad_y <- max(box_h) / 2 + 0.05
  ## never narrower than the legend (at 3.3-3.8 in it is clipped on both
  ## sides); 5.2 in = the PFC figure, so the three regions share one width
  W5 <- max(diff(range(x)) + 2 * pad_x + 0.2, as.numeric(getarg("--pill-min-width", "5.2")))
  H5 <- diff(range(y)) + 2 * pad_y + 0.2 + LEG_IN
  p5 <- ggraph(lay) +
    geom_edge_fan(aes(colour = peak_dir, width = binding_prob), alpha = 0.45, strength = 1.1) +
    geom_edge_loop(aes(colour = peak_dir, width = binding_prob), alpha = 0.45) +
    scale_edge_colour_manual(values = dir_cols, name = "Linked peak") +
    scale_edge_width(range = c(0.1, 0.5), name = "TF binding prob.") +
    geom_node_label(aes(filter = node_type != "TFs", label = name, fill = node_type, fontface = face),
                    size = PILL_PT / .pt, colour = "black", label.r = unit(0.45, "lines"),
                    label.padding = unit(0.14, "lines"), key_glyph = "point") +
    scale_fill_manual(values = col_mapper[gene_types], labels = type_lab[gene_types], name = NULL) +
    geom_node_label(aes(filter = node_type == "TFs", label = name), size = TF_PT / .pt,
                    fontface = "bold", fill = col_mapper[["TFs"]], colour = "black",
                    label.r = unit(0.05, "lines"), label.padding = unit(0.2, "lines")) +
    scale_x_continuous(expand = expansion(add = pad_x)) +
    scale_y_continuous(expand = expansion(add = pad_y)) +
    coord_fixed(clip = "off") +
    theme_void(base_size = 8) +
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.text = element_text(size = 7), legend.title = element_text(size = 7, face = "bold"),
          legend.key.size = unit(0.3, "cm"), legend.spacing.x = unit(0.2, "cm"),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(3, 3, 3, 3)) +
    guides(edge_colour = guide_legend(order = 1, nrow = 2, override.aes = list(edge_width = 0.8)),
           edge_width  = guide_legend(order = 2, nrow = 2),
           fill        = guide_legend(order = 3, override.aes = list(size = 3, shape = 21)))
  save_gg(p5, "network_tf_gene_pills", W5, H5)
  message(sprintf("[net] pill figure: %.1f x %.1f in (gene names %g pt, TF %g pt)", W5, H5, PILL_PT, TF_PT))
  message("[net] ggraph figures written (stress, layered, tf_gene)")
}
message(sprintf("[net] drew %d nodes, %d edges", vcount(ig), ecount(ig)))
message("[net] done")
