# ==============================================================================
# Figure 3 regional GO pathway dot plots, as THREE class panels (excitatory,
# inhibitory, glial) per direction instead of one file per cell type.
#
# WHY
#   The published 04-regional_pathway_dotplot_up.R computes one enrichment table,
#   then loops over cell types writing a separate PNG each (43 files). Assembled
#   into a figure that becomes ~16 mini heatmaps, none large enough to read.
#
# THE LAYOUT HERE
#   Three panels - Excitatory, Inhibitory, Glial/non-neuronal - each a single
#   dot plot:
#       x      = cell type (within that class)
#       y      = pathway (top N per cell type, unioned within the class)
#       facet  = region (PFC | INS | VST)
#       size   = -log10(p.adjust)
#       colour = mean logFC of the pathway's DE genes
#   Same enrichGO call as the published script, with these differences:
#     - DE genes: p_adj.glb < 0.05, |logFC| > --lfc (default 0.25, the published
#       rule), expressed in > 5% of nuclei in either group.
#     - rows are chosen by raw p from terms with >= 3 DE genes (the published
#       Count > 2), but a chosen term is drawn for every cell type whose DE genes
#       appear in it, significant or not (--min-select-genes / --min-draw-genes).
#       --padj < 1 adds a significance cutoff.
#     - top 2 per cell type per region, where the published script took 5.
#     - near-duplicate GO terms are merged (GOSemSim Wang similarity > 0.6).
#
# ENRICHMENT IS CACHED
#   compareCluster(fun = "enrichGO") over ~22 cell types x 3 regions is the slow
#   part and does not change when the layout is tweaked. It is cached to
#   fig3_pathway_enrichment_<direction>_allcounts*.rds; pass --refresh to
#   recompute. Always --refresh after the DE results change.
#
# INPUT   OUT_ROOT/03-DE_analysis/<region>/results/<region>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
# OUTPUT  OUT_ROOT/03-DE_analysis/figure_prototypes/ (enrichment cache + full CSV)
#         OUT_ROOT/03-DE_analysis/figure_prototypes/fig3_runs/<run tag>/
#           fig3_<class>_<dir>_top<N>_<horizontal|vertical>.png, _plotted.csv,
#           _omitted.csv, _simplify_merged.csv, fig3_panel_contents_<dir>.csv,
#           settings_<dir>.txt
#
# Usage: revision/docker/run.sh revision/tools/fig3_pathway_panels.R
#          [--direction up|down] [--top 2] [--lfc 0.25] [--universe genome|tested] [--refresh]
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
  library(clusterProfiler); library(org.Hs.eg.db); library(Homo.sapiens)
})

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }

DIRECTION <- match.arg(getarg("--direction", "up"), c("up", "down"))
TOP_N     <- as.integer(getarg("--top", "2"))
## DEG logFC threshold. The published DEG rule is |logFC| > 0.25
## (04-regional_pathway_dotplot_up.R:79), as in every other figure; --lfc 0
## takes any logFC > 0 / < 0.
LFC      <- as.numeric(getarg("--lfc", "0.25"))
## Enrichment background:
##   genome  compareCluster/enrichGO with no universe = every GO-annotated gene
##           (what the published script did)
##   tested  each region x cell type against ITS OWN tested genes - those that
##           passed the > 5% expression filter and so could have been called DE
UNIVERSE <- match.arg(getarg("--universe", "genome"), c("genome", "tested"))
## cache/run suffix encodes the logFC threshold and background (empty for
## --lfc 0 with the genome background)
ENRICH_SFX <- paste0(if (LFC > 0) paste0("_lfc", LFC) else "",
                     if (UNIVERSE == "tested") "_universe-tested" else "")
REFRESH   <- "--refresh" %in% args

## Keep only pathways picked up by at least this many cell types WITHIN the
## panel's class. Default 1 (off). At --top 3 with no filter, ~20 of 35 rows are
## a single dot - a pathway that is top-3 for one cell type in one region and
## absent everywhere else; --min-celltypes 2 removes them.
MIN_CT    <- as.integer(getarg("--min-celltypes", "1"))

## SIGNIFICANCE: off by default. The plot shows which pathways the DEGs belong
##   to, so, as published: top N terms by raw p per cell type per region, no
##   cutoff, every measurement of those terms drawn. --padj 0.05 removes every
##   cell type whose DE genes give no term at adj. p < 0.05 (e.g. INS Astro and
##   Micro down, PFC Vip/Endo/OPC up).
PADJ      <- as.numeric(getarg("--padj", "1"))

## GENES PER TERM. The published script kept only GO terms with > 2 of the cell
##   type's DE genes. Applied to every dot, that empties whole columns (INS Oligo;
##   VST Astro up, VST OPC, VST Oligo down: 3-13 DE genes, never 3 in one term)
##   and drops small dots elsewhere. The cache keeps every term and the rule is
##   split in two:
##     --min-select-genes 3  a term must have >= 3 DE genes in SOME cell type to
##                           become a row (the published Count > 2)
##     --min-draw-genes 1    a row is then drawn for every cell type whose DE
##                           genes appear in that term at all
##   --min-draw-genes 3 applies the published filter to every dot.
MIN_SEL_GENES  <- as.integer(getarg("--min-select-genes", "3"))
MIN_DRAW_GENES <- as.integer(getarg("--min-draw-genes", "1"))

## EMPTY COLUMNS. Default tested: every cell type pbDS tested gets a column,
##   labelled when it has no dots, so every major cell type is shown or its
##   absence explained. Cell types pbDS never tested (race rank deficiency)
##   cannot have dots; they are listed with the reason in
##   fig3_<class>_<dir>_top<N>_omitted.csv. none = no empty columns;
##   all = untested cell types as empty columns too.
SHOW_EMPTY <- match.arg(getarg("--show-empty", "tested"), c("all", "tested", "none"))

## A. GLIAL PANEL. "core" = Astro, Oligo, OPC, Micro.
##   Endo, VLMC, PC, SMC (vascular) and Macrophages, Immune are not glia, and
##   they contribute rows nothing else shares (VST Endo IL-5 terms, INS
##   Macrophages pH terms). "all" = every glial and non-neuronal cell type.
GLIA      <- match.arg(getarg("--glia", "core"), c("core", "all"))
GLIA_CORE <- c("Astro", "Oligo", "OPC", "Micro")

## B. MERGE NEAR-DUPLICATE GO TERMS. clusterProfiler's
##   simplify() rule (GOSemSim Wang similarity; among similar terms keep the
##   best) applied to each cell type's --simplify-pool best terms, then again
##   across the panel's selected terms. 0 = off.
##   Cutoff 0.6, not simplify()'s default 0.7: the duplicates that prompted this
##   sit just under 0.7 - response to virus / defense response to virus 0.63,
##   cilium assembly / cilium organization 0.65.
SIMPLIFY      <- as.numeric(getarg("--simplify", "0.6"))
SIMPLIFY_POOL <- as.integer(getarg("--simplify-pool", "20"))

## Layout. "horizontal" puts pathways on x and cell types on y, with the three
## regions stacked - a wide, short figure. "vertical" is the original
## orientation (pathways on y, regions side by side).
ORIENT    <- match.arg(getarg("--orient", "horizontal"), c("horizontal", "vertical"))

## SIZE. compact (default): pathway names wrap at ~26-28 characters so most sit
##   on two short lines (up to three), smaller fonts and dots, tighter spacing,
##   figure sized from the rows/columns actually drawn, 300 dpi, written to a
##   *_compact folder. large = bigger fonts and dots, 200 dpi.
SIZE <- match.arg(getarg("--size", "compact"), c("compact", "large"))
H_COL <- suppressWarnings(as.numeric(getarg("--h-col", NA)))   # in per pathway, horizontal
SZ <- if (SIZE == "compact") list(
  ## 28, not 24: at 24 the longest GO names need four lines and are cut
  ## with an ellipsis, and the extra third lines make every row taller
  lab_v = 28L, lab_h = 26L, max_lines = 3L, lineheight = 0.85,
  font = 11, axis_y = 9.5, axis_x = 9.5, strip = 10.5, leg_title = 9.5, leg_text = 8.5,
  title = 13, dot = c(1.2, 6), spacing = 0.8, key_h = 0.8, leg_space = 0.5,
  legend_in = 1.3, ylab_v = 1.9, ylab_h = 1.3, col_in = 0.36, dpi = 300
) else list(
  lab_v = 46L, lab_h = 34L, max_lines = 2L, lineheight = 0.95,
  font = 15, axis_y = 13, axis_x = 13, strip = 15, leg_title = 13, leg_text = 12,
  title = 19, dot = c(2, 11), spacing = 1.6, key_h = 1.0, leg_space = 0.8,
  legend_in = 2.4, ylab_v = 3.6, ylab_h = 2.6, col_in = 0.70, dpi = 200
)

## --abbrev: standard GO shorthand so
##   nearly every label fits on two lines. Display only - the plotted CSVs keep
##   the full term next to its label, and a label that would merge two
##   different terms onto one row falls back to the full names.
ABBREV <- "--abbrev" %in% args
abbreviate_go <- function(x) {
  rules <- c(
    "immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains" =
      "immune response (Ig-superfamily receptors)",
    "positive regulation of " = "pos. reg. of ",
    "negative regulation of " = "neg. reg. of ",
    "regulation of "          = "reg. of ",
    "biosynthetic process"    = "biosynthesis",
    "catabolic process"       = "catabolism",
    "metabolic process"       = "metabolism",
    "signaling pathway"       = "signaling",
    "signal transduction"     = "signaling")
  for (k in names(rules)) x <- gsub(k, rules[[k]], x, fixed = TRUE)
  x
}

## Pathway labels are long; horizontal needs them shorter to keep the angled
## axis from eating the panel.
LAB_CHARS <- if (ORIENT == "horizontal") SZ$lab_h else SZ$lab_v

OUT <- file.path(OUT_ROOT, "03-DE_analysis", "figure_prototypes")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
## _allcounts: the cache keeps terms with any number of DE genes.
CACHE <- file.path(OUT, sprintf("fig3_pathway_enrichment_%s_allcounts%s.rds", DIRECTION, ENRICH_SFX))

## Every panel variant gets its OWN folder, named from its settings.
## The enrichment cache and the full enrichment CSV are shared inputs and stay
## in OUT. Up and down write the same folder, so settings are per direction.
## NB eight arguments need eight conversions (count them when adding one): with
## a "%s" missing, sprintf silently drops an argument and runs overwrite each other.
RUN_TAG <- getarg("--tag", sprintf("glia-%s_simplify-%s_top%d%s%s%s%s%s", GLIA,
                  if (SIMPLIFY > 0) SIMPLIFY else "off", TOP_N,
                  if (PADJ < 1) paste0("_padj", PADJ) else "",
                  if (SIZE == "compact") "_compact" else "",
                  if (ABBREV) "_abbrev" else "",
                  if (MIN_DRAW_GENES < 3) paste0("_draw", MIN_DRAW_GENES, "gene") else "",
                  if (SHOW_EMPTY != "none") paste0("_empty-", SHOW_EMPTY) else ""))
if (!"--tag" %in% args) RUN_TAG <- paste0(RUN_TAG, ENRICH_SFX)
RUN_OUT <- file.path(OUT, "fig3_runs", RUN_TAG)
dir.create(RUN_OUT, recursive = TRUE, showWarnings = FALSE)
writeLines(c(paste("run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
             paste("args:", paste(args, collapse = " ")),
             paste("DE genes: p_adj.glb < 0.05, |logFC| >", LFC, "; expressed in > 5% of nuclei"),
             paste("enrichment background:", if (UNIVERSE == "tested")
               "genes tested for that region x cell type" else "all GO-annotated genes (no universe)"),
             paste("glia:", GLIA, if (GLIA == "core") paste0("(", paste(GLIA_CORE, collapse = " "), ")") else ""),
             paste("simplify:", if (SIMPLIFY > 0) paste0("Wang similarity > ", SIMPLIFY, ", pool ", SIMPLIFY_POOL) else "off"),
             paste("top_n:", TOP_N, "per cell type per region by raw p;",
                   if (PADJ < 1) paste("BH <", PADJ) else "no significance cutoff"),
             paste("show_empty:", SHOW_EMPTY),
             paste("size:", SIZE),
             paste("abbreviated GO names:", ABBREV),
             paste("genes per term: >=", MIN_SEL_GENES, "to select a row; >=", MIN_DRAW_GENES,
                   "to draw a dot"),
             paste0("DE genes: p_adj.glb < 0.05, |logFC| > ", LFC, ", expressed in > 5% of cells")),
           file.path(RUN_OUT, sprintf("settings_%s.txt", DIRECTION)))

## ---- cell-type classes -------------------------------------------------------
## MSNs are GABAergic projection neurons, so they sit with the inhibitory panel.
## They are also the defining striatal population and VST-only - if that reads
## oddly next to cortical interneurons, split them into a fourth panel.
##
## The nine MSN SUBTYPES are listed as well as "MSNs": with the subtypes
##   attached (revision/tools/patch_msn_subtypes.R) the VST DE reports
##   D1-NUDAP, D1-Shell/OT, D2-Shell/OT etc. instead of one "MSNs" row.
##   Without these names class_of() returns NA and those rows are silently
##   dropped from every panel.
MSN_SUBTYPES <- c("D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP",
                  "D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid")
CLASS <- list(
  Excitatory = c("L2/3 IT","L4 IT","L5 IT","L5 ET","L5/6 NP","L6 IT","L6 IT Car3",
                 "L6b","L6 CT"),
  Inhibitory = c("Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst",
                 "Sst Chodl","Vip","MSNs", MSN_SUBTYPES),
  Glial      = c("Astro","Oligo","OPC","Micro","Macrophages","Immune","Endo",
                 "PC","SMC","VLMC")
)
## CLASS_ALL keeps every name so types dropped by --glia core are still
## reported in omitted.csv and are not mistaken for unmapped labels.
CLASS_ALL <- CLASS
if (GLIA == "core") CLASS$Glial <- GLIA_CORE
class_of <- function(ct) {
  for (k in names(CLASS)) if (ct %in% CLASS[[k]]) return(k)
  NA_character_
}

## ---- DE input ----------------------------------------------------------------
load_de <- function() {
  regs <- c(INS = "INS", PFC = "PFC", VST = "VST")
  out <- lapply(names(regs), function(r) {
    f <- file.path(OUT_ROOT, "03-DE_analysis", r, "results",
                   sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq.RDS", r))
    if (!file.exists(f)) { message("[fig3] missing ", f); return(NULL) }
    x <- readRDS(f) %>% ungroup() %>%
      ## same expression filter as the published script
      dplyr::filter((avg_CTR_frq > 0.05) | (avg_HIV_frq > 0.05)) %>%
      dplyr::select(gene, cluster_id, logFC, logCPM, p_adj.glb, region)
    x$region <- r
    x
  })
  bind_rows(out)
}

## ---- enrichment (cached) -----------------------------------------------------
compute_enrichment <- function() {
  de_all <- load_de()
  message("[fig3] DE rows: ", nrow(de_all))

  sig <- de_all %>%
    dplyr::filter(p_adj.glb < 0.05,
                  if (DIRECTION == "up") logFC > LFC else logFC < -LFC) %>%
    mutate(cell_type = as.character(cluster_id))
  message("[fig3] significant ", DIRECTION, "-regulated (|logFC| > ", LFC, "): ", nrow(sig))

  gene_lists <- list()
  for (reg in unique(sig$region)) {
    for (cp in unique(sig$cell_type[sig$region == reg])) {
      sub <- sig %>% dplyr::filter(region == reg, cell_type == cp)
      if (nrow(sub) <= 1) next
      conv <- suppressMessages(
        AnnotationDbi::select(Homo.sapiens, keys = unique(sub$gene),
                              keytype = "SYMBOL", columns = "ENTREZID"))
      ent <- conv$ENTREZID[!is.na(conv$ENTREZID)]
      if (length(ent) > 1) gene_lists[[paste0(reg, "-", cp)]] <- ent
    }
  }
  message("[fig3] gene lists: ", length(gene_lists))

  if (UNIVERSE == "genome") {
    ego <- compareCluster(geneClusters = gene_lists, pAdjustMethod = "BH",
                          pvalueCutoff = 1, qvalueCutoff = 1,
                          keyType = "ENTREZID", OrgDb = org.Hs.eg.db,
                          readable = TRUE, ont = "BP", fun = "enrichGO")
    ego_df <- as.data.frame(ego@compareClusterResult)
  } else {
    ## compareCluster takes ONE universe for every cluster, but each region x
    ## cell type was tested on its own genes. So enrichGO runs per cluster with
    ## that cluster's universe - the same test compareCluster performs
    ## internally, which also adjusts p within each cluster.
    ego_df <- dplyr::bind_rows(lapply(names(gene_lists), function(k) {
      reg <- sub("-.*$", "", k); cp <- sub("^[A-Z]*-", "", k)
      uni_sym <- unique(de_all$gene[de_all$region == reg & as.character(de_all$cluster_id) == cp])
      uni <- suppressMessages(AnnotationDbi::select(Homo.sapiens, keys = uni_sym,
                                                    keytype = "SYMBOL", columns = "ENTREZID"))
      uni <- unique(uni$ENTREZID[!is.na(uni$ENTREZID)])
      e <- tryCatch(enrichGO(gene = gene_lists[[k]], universe = uni, keyType = "ENTREZID",
                             OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                             pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE),
                    error = function(err) NULL)
      if (is.null(e) || !nrow(as.data.frame(e))) return(NULL)
      d <- as.data.frame(e); d$Cluster <- k
      d
    }))
    message("[fig3] per-cluster tested-gene universes: ", length(gene_lists), " clusters")
  }
  ego_df <- ego_df %>%
    tidyr::extract(Cluster, into = c("region", "cell_type"), "^(^[A-Z]*)\\-(.*)")

  ## mean logFC of the DE genes behind each pathway - the colour channel
  ## No Count > 2 filter here: filtering before the cache would make 1-2 gene
  ##   terms unrecoverable. The published threshold is applied when rows are
  ##   SELECTED (MIN_SEL_GENES), not when dots are drawn.
  long <- ego_df %>%
    separate_longer_delim(geneID, delim = "/") %>%
    mutate(gene = geneID, pathway = Description) %>%
    left_join(sig %>% dplyr::select(region, cell_type, gene, logFC),
              by = c("region", "cell_type", "gene")) %>%
    group_by(region, cell_type, pathway, ID, p.adjust, pvalue) %>%
    summarise(avr_logFC = mean(logFC, na.rm = TRUE), n_genes = dplyr::n(),
              .groups = "drop")

  saveRDS(long, CACHE)
  message("[fig3] cached enrichment -> ", CACHE)
  long
}

long <- if (!REFRESH && file.exists(CACHE)) {
  message("[fig3] using cached enrichment (pass --refresh to recompute): ", CACHE)
  readRDS(CACHE)
} else compute_enrichment()

## Full enrichment table as CSV - every GO term for every region x cell type.
## The RDS cache is not something to hand to a co-author or a supplement.
utils::write.csv(long, file.path(OUT, sprintf("fig3_pathway_enrichment_%s_allcounts%s.csv", DIRECTION, ENRICH_SFX)),
                 row.names = FALSE)
message("[fig3] wrote fig3_pathway_enrichment_", DIRECTION, "_allcounts", ENRICH_SFX, ".csv (", nrow(long), " rows)")

## ---- which cell types exist / were tested, per region ------------------------
## The pseudobulk's assays are every cell type in the region; the DE result's
## cluster_ids are the ones pbDS actually tested (it silently drops the rest,
## e.g. cell types whose design is rank-deficient for race).
ct_status <- bind_rows(lapply(c("INS", "PFC", "VST"), function(r) {
  def <- file.path(OUT_ROOT, "03-DE_analysis", r, "results",
                   sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq.RDS", r))
  pbf <- file.path(OUT_ROOT, "03-DE_analysis", r, "data", "pseudo_bulk.rds")
  if (!file.exists(def)) return(NULL)
  tested <- unique(as.character(readRDS(def)$cluster_id))
  all_ct <- if (file.exists(pbf)) SummarizedExperiment::assayNames(readRDS(pbf)) else tested
  data.frame(region = r, cell_type = union(all_ct, tested), stringsAsFactors = FALSE) %>%
    mutate(tested = cell_type %in% tested)
}))
message("[fig3] cell types per region (tested/all): ",
        paste(tapply(ct_status$tested, ct_status$region,
                     function(x) paste0(sum(x), "/", length(x))), collapse = "  "))

## ---- assemble panels ---------------------------------------------------------
long$class <- vapply(long$cell_type, class_of, character(1))
unmapped <- sort(setdiff(unique(long$cell_type[is.na(long$class)]), unlist(CLASS_ALL)))
if (length(unmapped))
  message("[fig3] NOT in any class (add them to CLASS): ", paste(unmapped, collapse = ", "))
long <- long[!is.na(long$class), ]

long$region <- factor(long$region, levels = c("PFC", "INS", "VST"))

GO_SEM <- if (SIMPLIFY > 0)
  suppressMessages(GOSemSim::godata(annoDb = "org.Hs.eg.db", ont = "BP", computeIC = FALSE)) else NULL

## clusterProfiler:::simplify_internal's rule on a short list of GO IDs, so it
## runs in seconds on the cached table: for every term, among the terms whose
## Wang similarity to it exceeds the cutoff (itself included), only the best
## (lowest) score survives. The overall best term is never removed. Returns the
## kept IDs, with attr "merged" = removed term, the kept term it is most
## similar to, and that similarity.
simplify_ids <- function(ids, score, cutoff) {
  if (length(ids) < 2) return(structure(ids, merged = NULL))
  sim <- GOSemSim::termSim(ids, ids, semData = GO_SEM, method = "Wang")
  sim[is.na(sim)] <- 0                      # obsolete / unknown IDs never merge
  remove <- character(0)
  for (i in seq_along(ids)) {
    ii <- which(sim[, i] > cutoff)
    if (length(ii) < 2) next
    best <- ii[score[ii] == min(score[ii])]
    remove <- union(remove, ids[setdiff(ii, best)])
  }
  kept <- setdiff(ids, remove)
  merged <- if (length(remove)) data.frame(
    removed    = remove,
    into       = vapply(remove, function(r) kept[which.max(sim[r, kept])], character(1)),
    similarity = vapply(remove, function(r) max(sim[r, kept]), numeric(1)),
    stringsAsFactors = FALSE) else NULL
  structure(kept, merged = merged)
}

make_panel <- function(cls) {
  d <- long[long$class == cls, ]
  if (!nrow(d)) { message("[fig3] no data for ", cls); return(NULL) }

  ## top N per cell type per region by raw p, then keep every (cell type,
  ## region) measurement of those pathways. With --padj < 1, selection and
  ## drawing are restricted to terms under that BH cutoff.
  if (PADJ < 1) d <- d[d$p.adjust < PADJ, ]
  ## d  = every measurement that may be DRAWN;  ds = those that may be SELECTED
  d  <- d[d$n_genes >= MIN_DRAW_GENES, ]
  ds <- d[d$n_genes >= MIN_SEL_GENES, ]
  if (!nrow(ds)) { message("[fig3] no terms left for ", cls); return(NULL) }
  had_terms <- unique(paste(ds$region, ds$cell_type))
  if (SIMPLIFY > 0) {
    ## B, step 1 - per cell type: merge near-duplicates among its best
    ## SIMPLIFY_POOL terms (by raw p), then take the top N survivors
    merged_log <- list()
    sel <- ds %>% group_by(cell_type, region) %>% slice_min(pvalue, n = SIMPLIFY_POOL) %>%
      group_modify(function(g, k) {
        kept <- simplify_ids(g$ID, g$pvalue, SIMPLIFY)
        m <- attr(kept, "merged")
        if (!is.null(m))
          merged_log[[length(merged_log) + 1]] <<- cbind(level = paste(k$region, k$cell_type), m)
        dplyr::slice_min(g[g$ID %in% kept, ], pvalue, n = TOP_N)
      }) %>% ungroup()
    ## step 2 - across the panel, so two cell types' near-identical picks become
    ## one row: the term picked by more cell types wins, then the smaller p
    rank_tab <- sel %>% group_by(ID) %>%
      summarise(n_sel = dplyr::n(), best = min(pvalue), .groups = "drop") %>%
      arrange(desc(n_sel), best)
    kept_ids <- simplify_ids(rank_tab$ID, seq_len(nrow(rank_tab)), SIMPLIFY)
    if (!is.null(attr(kept_ids, "merged")))
      merged_log[[length(merged_log) + 1]] <- cbind(level = "panel", attr(kept_ids, "merged"))
    message("[fig3] ", cls, ": simplify ", SIMPLIFY, " - ", nrow(rank_tab),
            " selected terms -> ", length(kept_ids), " after merging across cell types")
    if (length(merged_log)) {
      ml <- do.call(rbind, merged_log)
      nm <- setNames(d$pathway, d$ID)
      ml$removed_name <- unname(nm[ml$removed]); ml$into_name <- unname(nm[ml$into])
      utils::write.csv(ml[, c("level", "removed", "removed_name", "into", "into_name", "similarity")],
                       file.path(RUN_OUT, sprintf("fig3_%s_%s_top%d_simplify_merged.csv",
                                                  tolower(cls), DIRECTION, TOP_N)),
                       row.names = FALSE)
    }
    keep <- data.frame(ID = as.character(kept_ids))
  } else {
    keep <- ds %>% group_by(cell_type, region) %>% slice_min(pvalue, n = TOP_N) %>%
      ungroup() %>% dplyr::select(ID) %>% distinct()
  }
  d <- d[d$ID %in% keep$ID, ]

  ## shared-pathway filter: how many distinct cell types select this pathway?
  breadth <- d %>% group_by(ID) %>%
    summarise(n_ct = dplyr::n_distinct(cell_type), .groups = "drop")
  before <- dplyr::n_distinct(d$ID)
  d <- d[d$ID %in% breadth$ID[breadth$n_ct >= MIN_CT], ]
  message("[fig3] ", cls, ": ", before, " pathways -> ",
          dplyr::n_distinct(d$ID), " after >= ", MIN_CT, " cell types")
  if (!nrow(d)) { message("[fig3] nothing left for ", cls); return(NULL) }

  ## Wrap at word boundaries rather than substr() - a hard character cut leaves
  ## labels like "protein localization to microtubule cytoskelet" severed mid-word.
  ## Long names wrap onto a second line; only genuinely enormous GO terms are
  ## then elided, and with an ellipsis so it is visibly a truncation.
  .wrap <- function(x, width) {
    out <- vapply(x, function(s) {
      w <- strwrap(s, width = width, simplify = TRUE)
      kept <- if (length(w) > SZ$max_lines) w[seq_len(SZ$max_lines)] else w   # line cap from --size
      ## Truncation test compares WORD CONTENT, not character count. strwrap
      ## drops the space at each break, so a character comparison makes every
      ## wrapped label look shortened and appends a false ellipsis.
      truncated <- length(w) > SZ$max_lines
      paste0(paste(kept, collapse = "\n"), if (truncated) "…" else "")
    }, character(1))
    unname(out)
  }
  lab_src <- if (ABBREV) abbreviate_go(d$pathway) else d$pathway
  ## an abbreviation must never merge two different terms onto one row
  clash <- tapply(d$ID, lab_src, function(z) length(unique(z)) > 1)
  if (any(clash)) lab_src <- ifelse(lab_src %in% names(clash)[clash], d$pathway, lab_src)
  d$pathway_short <- .wrap(lab_src, LAB_CHARS)
  ## order pathways by how broadly they come up, so the busiest sit together
  ord <- d %>% group_by(pathway_short) %>%
    summarise(n = dplyr::n(), best = min(p.adjust), .groups = "drop") %>%
    arrange(n, -best)
  d$pathway_short <- factor(d$pathway_short, levels = ord$pathway_short)
  ## One column per (region, cell type), INCLUDING cell types with no dots.
  ## The plotted variable is "region||cell type" so the same name can carry a
  ## different label per region (INS Astro has dots, another region's may not).
  d$cell_type <- as.character(d$cell_type)
  d$region    <- as.character(d$region)
  have <- dplyr::distinct(d[, c("region", "cell_type")])
  ## what this class has no dots for, and why - written for the caption
  ## whatever --show-empty is
  omit <- ct_status[ct_status$cell_type %in% CLASS_ALL[[cls]] &
                    !paste(ct_status$region, ct_status$cell_type) %in%
                      paste(have$region, have$cell_type), ]
  omit$reason <- ifelse(!omit$cell_type %in% CLASS[[cls]], "not core glia (--glia core)",
                 ifelse(!omit$tested, "not tested by DE (dropped by pbDS)",
                 ifelse(paste(omit$region, omit$cell_type) %in% had_terms,
                        "its top terms merged into near-duplicates it lacks (--simplify)",
                        if (PADJ < 1) "no GO term below the BH cutoff"
                        else "too few DE genes: none of them fall in any plotted pathway")))
  omit <- omit[order(match(omit$region, c("PFC", "INS", "VST")), match(omit$cell_type, CLASS_ALL[[cls]])), ]
  utils::write.csv(omit[, c("region", "cell_type", "reason")],
                   file.path(RUN_OUT, sprintf("fig3_%s_%s_top%d_omitted.csv", tolower(cls), DIRECTION, TOP_N)),
                   row.names = FALSE)
  st <- if (SHOW_EMPTY == "none") have else
    ct_status[ct_status$cell_type %in% CLASS[[cls]] &
              (SHOW_EMPTY == "all" | ct_status$tested), c("region", "cell_type")]
  st <- dplyr::distinct(bind_rows(st, have))
  st <- left_join(st, ct_status, by = c("region", "cell_type"))
  st$tested[is.na(st$tested)] <- TRUE
  st$has_dots <- paste(st$region, st$cell_type) %in% paste(have$region, have$cell_type)
  ## NB the n.s./too-few choice is a scalar if(), not ifelse(): ifelse() with a
  ## length-1 test returns ONE value, which gets recycled - every tested-but-empty
  ## column would be labelled with the class's first cell type.
  empty_tag <- if (PADJ < 1) " (n.s.)" else " (too few DEGs)"
  st$label <- ifelse(st$has_dots, st$cell_type,
              ifelse(!st$tested, paste0(st$cell_type, " (not tested)"),
                                 paste0(st$cell_type, empty_tag)))
  st <- st[order(match(st$region, c("PFC", "INS", "VST")), match(st$cell_type, CLASS[[cls]])), ]
  st$key <- paste(st$region, st$cell_type, sep = "||")
  d$key  <- paste(d$region, d$cell_type, sep = "||")
  miss <- st[!st$has_dots, ]
  if (nrow(miss)) {
    ## invisible placeholder (NA size and colour) parked on an existing pathway
    ## row, only so the facet keeps the column
    blank <- data.frame(class = cls, region = miss$region, cell_type = miss$cell_type,
                        key = miss$key, stringsAsFactors = FALSE)
    blank$pathway_short <- factor(levels(d$pathway_short)[1], levels = levels(d$pathway_short))
    d <- bind_rows(d, blank)
  }
  message("[fig3] ", cls, ": ", sum(st$has_dots), " columns with dots, ",
          sum(!st$has_dots), " empty columns drawn | omitted: ",
          paste(names(table(omit$reason)), table(omit$reason), sep = " = ", collapse = "; "))
  d$region <- factor(d$region, levels = c("PFC", "INS", "VST"))
  d$key    <- factor(d$key, levels = st$key)
  key_lab  <- setNames(st$label, st$key)
  lab_fun  <- function(k) unname(key_lab[k])

  ## SEQUENTIAL colour, not diverging. With --direction up every avr_logFC is
  ## positive, so a blue-to-red diverging ramp reads as "blue = downregulated"
  ## when it only means "smaller positive fold change". Pale-to-red for
  ## up, pale-to-blue for down, no midpoint, no sign implied.
  ramp_hi <- if (DIRECTION == "up") "#B2182B" else "#2166AC"

  base <- ggplot(d, aes(size = -log10(p.adjust), colour = avr_logFC)) +
    geom_point(na.rm = TRUE) +   # placeholder rows for empty columns are NA
    scale_size(range = SZ$dot, name = expression(-log[10]~adj.~p)) +
    ## For DOWN every value is negative, so the most strongly down-regulated
    ## pathways sit at the LOW end. The ramp must run dark-at-low for down, or the
    ## strongest effects render palest.
    scale_colour_gradient(low  = if (DIRECTION == "up") "grey88" else ramp_hi,
                          high = if (DIRECTION == "up") ramp_hi else "grey88",
                          name = "mean logFC") +
    ## Just the class name. Selection rules (top N per cell type, shared by
    ## >= MIN_CT, direction) belong in the figure caption, not burned into the
    ## image - they are also recorded in the filename and in fig3_panel_contents_<dir>.csv.
    labs(title = cls, x = NULL, y = NULL) +
    theme_cowplot(font_size = SZ$font) +
    theme(## nothing but dots inside the panel - no grid, no ticks, no axis line
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          ## a light box per region: enough to separate, too faint to compete
          ## with the dots for attention
          panel.border = element_rect(colour = "grey78", fill = NA, linewidth = 0.6),
          panel.spacing = unit(SZ$spacing, "lines"),
          ## pathway names are what the reader navigates by
          strip.text = element_text(face = "bold", size = SZ$strip, colour = "white"),
          legend.position = "right",
          legend.title = element_text(size = SZ$leg_title, face = "bold"),
          legend.text = element_text(size = SZ$leg_text),
          legend.key.height = unit(SZ$key_h, "lines"),
          legend.spacing.y = unit(SZ$leg_space, "lines"),
          plot.title = element_text(size = SZ$title, face = "bold"),
          plot.title.position = "plot",   # title aligns to the canvas, not the panel
          plot.margin = margin(10, 14, 6, 10))

  if (ORIENT == "horizontal") {
    ## facet_grid + space="free_y" drops the cell types a region has no dots for
    ## AND sizes each region band to the rows it actually uses - that is what
    ## removes the blank Immune/Macrophages strips.
    ## horizontal: pathways on x, so THEY get the bold treatment
    base + aes(x = pathway_short, y = key) + scale_y_discrete(labels = lab_fun) +
      facet_grid(region ~ ., scales = "free_y", space = "free_y", switch = "y") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = SZ$axis_x,
                                       face = "bold", lineheight = SZ$lineheight),
            axis.text.y = element_text(size = SZ$axis_y + 1),
            strip.placement = "outside")
  } else {
    ## vertical: pathways on y
    base + aes(x = key, y = pathway_short) + scale_x_discrete(labels = lab_fun) +
      facet_grid(. ~ region, scales = "free_x", space = "free_x") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = SZ$axis_x),
            axis.text.y = element_text(size = SZ$axis_y, face = "bold", lineheight = SZ$lineheight))
  }
}

## Region strip colours, matching the palette the combined figures already use
## (03-integration_tsne_refined.R:77). Applied after the plot is built because
## ggplot has no per-strip fill aesthetic.
REGION_FILL <- c(PFC = "#E69F00", INS = "#56B4E9", VST = "#CC79A7")

colour_strips <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  strips <- which(grepl("^strip-", g$layout$name))
  ## strip order follows the facet levels (PFC, INS, VST)
  labs <- levels(droplevels(factor(p$data$region,
                                   levels = c("PFC", "INS", "VST"))))
  for (i in seq_along(strips)) {
    k <- strips[i]
    lab <- if (i <= length(labs)) labs[i] else NA
    fill <- if (!is.na(lab) && lab %in% names(REGION_FILL)) REGION_FILL[[lab]] else "grey40"
    gr <- g$grobs[[k]]
    j <- which(vapply(gr$grobs[[1]]$children, function(x) inherits(x, "rect"), logical(1)))
    if (length(j)) gr$grobs[[1]]$children[[j[1]]]$gp$fill <- fill
    g$grobs[[k]] <- gr
  }
  g
}

## Both orientations are written each run, so they can be compared directly
## rather than by rerunning with a different flag.
ORIENTS <- if (identical(getarg("--orient", ""), "")) c("horizontal", "vertical") else ORIENT

for (cls in names(CLASS)) for (o in ORIENTS) {
  ORIENT    <<- o
  LAB_CHARS <<- if (o == "horizontal") SZ$lab_h else SZ$lab_v
  p <- make_panel(cls)
  if (is.null(p)) next

  ## size from what actually survived the filters, not the unfiltered table
  pd     <- p$data
  n_path <- dplyr::n_distinct(pd$pathway_short)
  ## one column (vertical) / row (horizontal) per region x cell type, empty ones included
  n_ct   <- dplyr::n_distinct(pd$key)
  n_rows <- n_ct

  ## Width must cover the y-axis labels and the legend, not just the panel -
  ## sizing to the data alone truncates the title and squeezes the key.
  LEGEND_IN <- SZ$legend_in
  YLAB_IN   <- if (o == "horizontal") SZ$ylab_h else SZ$ylab_v   # wrapped pathway names / "(not tested)" labels are wide
  ## wrapped labels take two lines, so rows need more height
  lab_lines <- lengths(regmatches(as.character(unique(pd$pathway_short)),
                                  gregexpr("\n", as.character(unique(pd$pathway_short))))) + 1
  n_lines <- max(1, mean(lab_lines))
  if (SIZE == "large") {
    if (o == "horizontal") {
      w <- max(11, YLAB_IN + LEGEND_IN + 0.46 * n_path)
      h <- max(6,  3.0 + 0.46 * n_rows + 0.5 * n_lines)
    } else {
      w <- max(11, YLAB_IN + LEGEND_IN + 0.70 * n_ct)
      h <- max(7,  4.0 + 0.40 * n_path * n_lines)   # +1 in for angled "(not tested)" labels
    }
  } else {
    ## compact: a discrete axis spaces rows EVENLY, so each row must fit the
    ## TALLEST label. At 0.115 in per line (90th percentile) two adjacent
    ## three-line labels run into each other; a rendered line of 9.5 pt is ~0.14 in.
    q_lines <- as.numeric(stats::quantile(lab_lines, 0.9))
    if (o == "horizontal") {
      ## A label wrapped onto k lines and turned 45 degrees spans about 0.22 in
      ## of x per line (at 0.19 adjacent 3-line labels touch), so a fixed
      ## 0.34 in per pathway lets 2- and 3-line labels run into their neighbours.
      ## The pitch follows the longest label actually drawn; --h-col sets it by hand.
      pitch <- if (!is.na(H_COL)) H_COL else max(0.34, 0.10 + 0.22 * max(lab_lines))
      w <- max(6, YLAB_IN + LEGEND_IN + pitch * n_path)
      h <- max(4, 1.6 + 0.24 * n_rows + 0.55 * q_lines)
      message(sprintf("[fig3] horizontal: %d pathways x %.2f in = %.1f in wide (labels up to %d lines)",
                      n_path, pitch, w, max(lab_lines)))
    } else {
      w <- max(5, YLAB_IN + LEGEND_IN + SZ$col_in * n_ct + 0.3)
      h <- max(4, 1.6 + (0.05 + 0.14 * max(lab_lines)) * n_path)
    }
  }

  g <- colour_strips(p)
  f <- file.path(RUN_OUT, sprintf("fig3_%s_%s_top%d_%s.png",
                              tolower(cls), DIRECTION, TOP_N, o))
  ## bg = "white": theme_cowplot leaves the background transparent, and viewers
  ## that paint transparency black show black text on black (panels look blank).
  ggsave(f, g, width = w, height = h, dpi = SZ$dpi, limitsize = FALSE, bg = "white")
  ## the exact rows drawn in this panel, once per class (same for both orientations)
  if (o == ORIENTS[1])
    utils::write.csv(pd[!is.na(pd$ID), ] %>%
                       dplyr::mutate(label = gsub("\n", " ", as.character(pathway_short))) %>%
                       dplyr::select(dplyr::any_of(c("class", "region", "cell_type",
                       "pathway", "label", "ID", "p.adjust", "pvalue", "avr_logFC", "n_genes"))),
                     file.path(RUN_OUT, sprintf("fig3_%s_%s_top%d_plotted.csv",
                                            tolower(cls), DIRECTION, TOP_N)),
                     row.names = FALSE)
  message("[fig3] wrote ", basename(f), "  (", n_path, " pathways x ",
          n_ct, " cell types, ", n_rows, " filled rows)")
}

## ---- what the reader is actually looking at ----------------------------------
summ <- long %>% group_by(class, region) %>%
  summarise(cell_types = dplyr::n_distinct(cell_type),
            pathways = dplyr::n_distinct(pathway), .groups = "drop")
message("\n########## panel contents ##########")
print(as.data.frame(summ), row.names = FALSE)

## Everything the caption needs, written out rather than drawn on the figure.
summ$direction      <- DIRECTION
summ$top_n_per_ct   <- TOP_N
summ$min_cell_types <- MIN_CT
summ$padj_threshold <- PADJ
summ$gene_universe  <- "GO BP, enrichGO, BH-adjusted"
summ$de_filter      <- "p_adj.glb < 0.05; expressed in > 5% of cells in either group"
summ$glia           <- GLIA
summ$simplify       <- SIMPLIFY
utils::write.csv(summ, file.path(RUN_OUT, sprintf("fig3_panel_contents_%s.csv", DIRECTION)),
                 row.names = FALSE)

message("\nCaption facts: top ", TOP_N, " ", DIRECTION, "-regulated GO BP terms per ",
        "cell type per region",
        if (PADJ < 1) paste0(" among terms with BH-adjusted p < ", PADJ,
                             "; only significant enrichments are shown")
        else " by raw p with no significance cutoff",
        paste0(", chosen from terms containing >= ", MIN_SEL_GENES, " of the cell type's DE genes;",
               " every measurement of a chosen term is drawn, significant or not",
               if (MIN_DRAW_GENES < MIN_SEL_GENES) paste0(", including cell types with only ",
                 paste(seq(MIN_DRAW_GENES, MIN_SEL_GENES - 1), collapse = "-"), " genes in the term") else ""),
        if (SHOW_EMPTY == "none")
          paste0("; cell types not tested by DE or with too few DE genes for any GO term are omitted",
                 " (listed in fig3_<class>_", DIRECTION, "_top", TOP_N, "_omitted.csv)")
        else if (SHOW_EMPTY == "tested")
          paste0("; every cell type tested by DE has a column, labelled (too few DEGs) when none of its",
                 " DE genes fall in a plotted pathway; untested cell types are listed in",
                 " fig3_<class>_", DIRECTION, "_top", TOP_N, "_omitted.csv")
        else "; every cell type is shown, empty columns labelled (not tested) or (too few DEGs)",
        if (SIMPLIFY > 0) paste0("; redundant GO terms merged (Wang semantic similarity > ", SIMPLIFY,
                                 ", clusterProfiler simplify rule)") else "",
        if (GLIA == "core") "; glial panel limited to astrocytes, oligodendrocytes, OPCs and microglia" else "",
        if (MIN_CT > 1) paste0(", retained where shared by >= ", MIN_CT, " cell types") else "",
        ". Dot size -log10 BH-adjusted p; colour mean logFC of the term's DE genes. ",
        "Plotted rows in fig3_<class>_", DIRECTION, "_top", TOP_N, "_plotted.csv")
message("\n[fig3] done - ", RUN_OUT)
