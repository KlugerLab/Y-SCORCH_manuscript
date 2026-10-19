# ==============================================================================
# Donor counts and per-donor nuclei for every figure panel that makes a
# statistical claim (reviewer request for per-donor sample sizes and donor counts
# in figure legends).
#
# Every number is read from the object or table the analysis itself used, with
# that analysis's own inclusion rule - nothing is re-derived from the cohort size.
# Inputs (under OUT_ROOT):
#   regional RNA DE      03-DE_analysis/<R>/data/pseudo_bulk.rds, pbDS min_cells 10;
#                        "tested" = cell type present in <R>_HIV-CTR_res_edgeR_filtered_withFreq.RDS
#   cross-region RNA DE  03-DE_analysis/cross_region/data/pseudo_bulk.rds, >= 10 nuclei per
#                        region x donor, rows of mm.rds; checked against n_samples / n_donors
#                        in cross_region_donor_aware_DE.csv
#   ATAC DA / gene act.  04-ATAC_DE/<R>/data_{ATAC,GeneACT}/pseudo_bulk.rds, pbDS min_cells 10
#   chromVAR             donors in <CT>_chromvar_donor_means.csv (>= 20 nuclei), per region
#                        and pooled
#   microglia            PCA > 20 nuclei (09-microglia_pca.R); heatmap > 30
#                        (cross_region_celltype_heatmap.R); state markers RNA > 50 / ATAC > 30
#                        (plotting/07, 08); signature scores >= 10 (microglia_signatures/*_scores.csv)
#   astrocytes           composition (astro_per_donor_counts), paired Reactive-PWH DE and
#                        S14b (usable donors), S15 subtype DE (n_donors + log), MAGICAL samples
#
# RNAscope and NICHES panels are not covered.
#
# Output: OUT_ROOT/legend_donor_counts/
#           panel_counts_long.csv   one row per figure panel x region x cell type
#           panel_summary.csv       one row per figure panel x region (the legend numbers)
#           legend_text.md          legend-ready sentences, numbers filled in from the tables
#
# Usage: revision/docker/run.sh revision/tools/legend_donor_counts.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(SingleCellExperiment); library(dplyr) })

REGIONS <- c("INS", "PFC", "VST")
OUT <- file.path(OUT_ROOT, "legend_donor_counts")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DE  <- file.path(OUT_ROOT, "03-DE_analysis")
AT  <- file.path(OUT_ROOT, "04-ATAC_DE")
AS  <- file.path(OUT_ROOT, "06-astrocyte_analysis")
cond_chr <- function(x) as.character(pwh_relabel(x))

## ---- helpers -----------------------------------------------------------------
## muscat keeps nuclei per sample x cluster in int_colData(pb)$n_cells
nc_long <- function(pb) {
  nc <- int_colData(pb)$n_cells
  if (is.null(names(nc))) names(nc) <- colnames(pb)
  do.call(rbind, lapply(names(nc), function(s)
    data.frame(sample = s, cell_type = names(nc[[s]]), n = as.integer(nc[[s]]),
               stringsAsFactors = FALSE)))
}
ei_cond <- function(pb) {
  ei <- S4Vectors::metadata(pb)$experiment_info
  setNames(cond_chr(ei$group_id), as.character(ei$sample_id))
}
tested_types <- function(f) {
  x <- readRDS(f)
  if (!is.data.frame(x) && !is.null(x$table)) x <- dplyr::bind_rows(lapply(x$table[[1]], identity))
  unique(as.character(x$cluster_id))
}
## one row per (region, cell type): donors passing a nuclei rule, split by condition
summarise_pass <- function(long, rule, min_n) {
  long %>% dplyr::filter(if (rule == ">=") n >= min_n else n > min_n) %>%
    group_by(region, cell_type) %>%
    summarise(donors_PWH = sum(condition == "PWH"), donors_PWoH = sum(condition == "PWoH"),
              donors_total = dplyr::n(), nuclei_total = sum(n),
              nuclei_per_donor_median = median(n), nuclei_per_donor_min = min(n),
              nuclei_per_donor_max = max(n), .groups = "drop")
}
ROWS <- list()
add <- function(df, figure, panel, analysis, unit, rule) {
  df$figure <- figure; df$panel <- panel; df$analysis <- analysis
  df$unit <- unit; df$inclusion_rule <- rule
  ROWS[[length(ROWS) + 1]] <<- df
  invisible(df)
}

## ---- 1. regional RNA DE (Fig 3, Fig 4, S6, S10 claims) --------------------------
rna_long <- do.call(rbind, lapply(REGIONS, function(R) {
  pb <- readRDS(file.path(DE, R, "data", "pseudo_bulk.rds"))
  L <- nc_long(pb); L$region <- R; L$condition <- ei_cond(pb)[L$sample]; L
}))
rna_tested <- do.call(rbind, lapply(REGIONS, function(R)
  data.frame(region = R, cell_type = tested_types(
    file.path(DE, R, "results", sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq.RDS", R))))))
reg_rna <- summarise_pass(rna_long, ">=", 10) %>%
  dplyr::inner_join(rna_tested, by = c("region", "cell_type"))
message("[leg] regional RNA DE: ", paste(sprintf("%s %d types", REGIONS,
        sapply(REGIONS, function(R) sum(reg_rna$region == R))), collapse = ", "))

## cohort, all nuclei (Fig 1)
cohort_rna <- rna_long %>% group_by(region, sample, condition) %>% summarise(n = sum(n), .groups = "drop") %>%
  group_by(region) %>% summarise(cell_type = "all", donors_PWH = sum(condition == "PWH"),
    donors_PWoH = sum(condition == "PWoH"), donors_total = dplyr::n(), nuclei_total = sum(n),
    nuclei_per_donor_median = median(n), nuclei_per_donor_min = min(n),
    nuclei_per_donor_max = max(n), .groups = "drop")
add(cohort_rna, "Fig 1", "a-f", "snRNA cohort after QC", "donor (one library per donor per region)", "all nuclei")
add(reg_rna, "Fig 3", "a-c", "regional pseudobulk edgeR DE", "donor", ">= 10 nuclei of the cell type (pbDS min_cells)")
add(reg_rna, "Fig 4", "a-c", "regional pseudobulk edgeR DE (per region)", "donor", ">= 10 nuclei of the cell type")
add(reg_rna[reg_rna$cell_type == "Endo", ], "Fig S6", "", "Endo regional DE (violins show nuclei)", "donor", ">= 10 nuclei")
add(reg_rna[reg_rna$cell_type == "Micro", ], "Fig S10", "", "Micro regional DE (violins show nuclei)", "donor", ">= 10 nuclei")

## ---- 2. cross-region RNA DE (Fig 5b/c, Fig 6a/b, Fig 2F x-axis, S7) -------------
xpb <- readRDS(file.path(DE, "cross_region", "data", "pseudo_bulk.rds"))
mm  <- readRDS(file.path(DE, "cross_region", "data", "mm.rds"))
xcd <- as.data.frame(colData(xpb))
x_long <- nc_long(xpb)
x_long$region    <- as.character(xcd[x_long$sample, "region"])
x_long$donor     <- as.character(xcd[x_long$sample, "sample_ids"])
x_long$condition <- cond_chr(xcd[x_long$sample, "condition_new"])
donor_cond <- unique(x_long[, c("donor", "condition")]); donor_cond <- setNames(donor_cond$condition, donor_cond$donor)

da <- utils::read.csv(file.path(DE, "cross_region", "donor_aware", "cross_region_donor_aware_DE.csv"),
                      stringsAsFactors = FALSE)
da_n <- unique(da[, c("cluster_id", "n_samples", "n_donors", "rho")])
xr <- x_long %>% dplyr::filter(n >= 10, sample %in% rownames(mm)) %>%
  group_by(cell_type) %>%
  summarise(region = "cross-region",
            samples_INS = sum(region == "INS"), samples_PFC = sum(region == "PFC"),
            samples_VST = sum(region == "VST"), samples_total = dplyr::n(),
            donors_PWH = dplyr::n_distinct(donor[condition == "PWH"]),
            donors_PWoH = dplyr::n_distinct(donor[condition == "PWoH"]),
            donors_total = dplyr::n_distinct(donor), nuclei_total = sum(n),
            nuclei_per_donor_median = median(n), nuclei_per_donor_min = min(n),
            nuclei_per_donor_max = max(n), .groups = "drop") %>%
  dplyr::inner_join(da_n, by = c("cell_type" = "cluster_id"))
bad <- xr[xr$samples_total != xr$n_samples | xr$donors_total != xr$n_donors, ]
if (nrow(bad)) {
  message("[leg] WARNING cross-region counts differ from the DE table for: ",
          paste(bad$cell_type, collapse = ", "))
  print(as.data.frame(bad[, c("cell_type", "samples_total", "n_samples", "donors_total", "n_donors")]))
} else message("[leg] cross-region counts match n_samples / n_donors of the DE table for all ",
               nrow(xr), " cell types")
xr$check_vs_DE_table <- ifelse(xr$samples_total == xr$n_samples & xr$donors_total == xr$n_donors, "match", "DIFFERS")
add(xr, "Fig S7", "a-c", "cross-region donor-aware limma-voom + duplicateCorrelation(donor)",
    "region x donor sample; donor as correlated block", ">= 10 nuclei per region x donor")
add(xr[xr$cell_type == "Micro", ], "Fig 5", "c", "cross-region donor-aware DE, microglia GO", "region x donor; donor block", ">= 10 nuclei")
add(xr[xr$cell_type == "Astro", ], "Fig 6", "b", "cross-region donor-aware DE, astrocyte GO", "region x donor; donor block", ">= 10 nuclei")
add(xr[xr$cell_type %in% c("Astro", "Micro", "Endo"), ], "Fig 2", "F (x-axis)", "TF RNA logFC, cross-region donor-aware DE",
    "region x donor; donor block", ">= 10 nuclei")

## sample-column panels (heatmaps / PCA): per region, nuclei rule of that script
col_panel <- function(ct, min_n, rule = ">") {
  x_long %>% dplyr::filter(cell_type == ct, if (rule == ">") n > min_n else n >= min_n) %>%
    group_by(region, cell_type) %>%
    summarise(donors_PWH = sum(condition == "PWH"), donors_PWoH = sum(condition == "PWoH"),
              donors_total = dplyr::n(), nuclei_total = sum(n), nuclei_per_donor_median = median(n),
              nuclei_per_donor_min = min(n), nuclei_per_donor_max = max(n), .groups = "drop")
}
add(col_panel("Micro", 20), "Fig 5", "a", "microglia PCA (DESeq2 vst, top 300 genes)", "donor (one point per sample)", "> 20 microglia")
add(col_panel("Micro", 30), "Fig 5", "b", "microglia top-DEG heatmap (columns = samples)", "region x donor", "> 30 microglia")
add(col_panel("Astro", 30), "Fig 6", "a", "astrocyte top-DEG heatmap (columns = samples)", "region x donor", "> 30 astrocytes")
add(col_panel("Micro", 50), "Fig S9", "a", "microglia state markers, RNA (columns = samples)", "region x donor", "> 50 microglia")

## ---- 3. ATAC: DA peaks and gene activity (Fig 2B/C/E, S8) ------------------------
atac_long <- list(); atac_tested <- list()
for (assay in c("ATAC", "GeneACT")) {
  atac_long[[assay]] <- do.call(rbind, lapply(REGIONS, function(R) {
    pb <- readRDS(file.path(AT, R, paste0("data_", assay), "pseudo_bulk.rds"))
    L <- nc_long(pb); L$region <- R; L$condition <- ei_cond(pb)[L$sample]; L
  }))
  atac_tested[[assay]] <- do.call(rbind, lapply(REGIONS, function(R)
    data.frame(region = R, cell_type = tested_types(
      file.path(AT, R, paste0("data_", assay), sprintf("%s_HIV-CTR_res_edgeR_filtered_withFreq.RDS", R))))))
}
same_nuclei <- isTRUE(all.equal(atac_long$ATAC[order(atac_long$ATAC$region, atac_long$ATAC$sample, atac_long$ATAC$cell_type), c("region", "sample", "cell_type", "n")],
                                atac_long$GeneACT[order(atac_long$GeneACT$region, atac_long$GeneACT$sample, atac_long$GeneACT$cell_type), c("region", "sample", "cell_type", "n")],
                                check.attributes = FALSE))
message("[leg] ATAC and gene-activity pseudobulks hold the same nuclei: ", same_nuclei)
cohort_atac <- atac_long$ATAC %>% group_by(region, sample, condition) %>% summarise(n = sum(n), .groups = "drop") %>%
  group_by(region) %>% summarise(cell_type = "all", donors_PWH = sum(condition == "PWH"),
    donors_PWoH = sum(condition == "PWoH"), donors_total = dplyr::n(), nuclei_total = sum(n),
    nuclei_per_donor_median = median(n), nuclei_per_donor_min = min(n),
    nuclei_per_donor_max = max(n), .groups = "drop")
add(cohort_atac, "Fig 2", "A", "snATAC cohort after QC", "donor", "all nuclei")
reg_da <- summarise_pass(atac_long$ATAC, ">=", 10) %>% dplyr::inner_join(atac_tested$ATAC, by = c("region", "cell_type"))
reg_ga <- summarise_pass(atac_long$GeneACT, ">=", 10) %>% dplyr::inner_join(atac_tested$GeneACT, by = c("region", "cell_type"))
add(reg_da, "Fig 2", "B", "DA peaks, regional pseudobulk edgeR", "donor", ">= 10 nuclei (pbDS min_cells)")
add(reg_da, "Fig S8", "a", "DA peaks, regional pseudobulk edgeR", "donor", ">= 10 nuclei")
add(reg_ga, "Fig S8", "b", "differential gene activity, regional pseudobulk edgeR", "donor", ">= 10 nuclei")
both <- dplyr::inner_join(reg_rna[, c("region", "cell_type", "donors_total")],
                          reg_ga[, c("region", "cell_type", "donors_total")],
                          by = c("region", "cell_type"), suffix = c("_RNA", "_GA"))
add(dplyr::inner_join(reg_ga, both[, c("region", "cell_type", "donors_total_RNA")], by = c("region", "cell_type")),
    "Fig 2", "C, E", "RNA DE vs gene-activity DE (both tested; GA donors shown, RNA donors in donors_total_RNA)",
    "donor", ">= 10 nuclei in each modality")
mg_atac <- atac_long$GeneACT %>% dplyr::filter(cell_type == "Micro", n > 30) %>% group_by(region, cell_type) %>%
  summarise(donors_PWH = sum(condition == "PWH"), donors_PWoH = sum(condition == "PWoH"), donors_total = dplyr::n(),
            nuclei_total = sum(n), nuclei_per_donor_median = median(n), nuclei_per_donor_min = min(n),
            nuclei_per_donor_max = max(n), .groups = "drop")
add(mg_atac, "Fig S9", "b", "microglia state markers, gene activity (columns = samples)", "region x donor", "> 30 microglia (ATAC)")

## ---- 4. chromVAR (Fig 2D/F, S11) ----------------------------------------------------
CV <- file.path(AT, "motif_chromvar")
cv_one <- function(R, ct) {
  f <- file.path(CV, R, paste0(gsub("[^A-Za-z0-9_.+-]", "_", ct), "_chromvar_donor_means.csv"))
  if (!file.exists(f)) return(NULL)
  don <- setdiff(names(utils::read.csv(f, check.names = FALSE, nrows = 1)), c("motif", "name"))
  a <- atac_long$ATAC[atac_long$ATAC$region == R & atac_long$ATAC$cell_type == ct & atac_long$ATAC$sample %in% don, ]
  data.frame(region = R, cell_type = ct, donor = don, condition = donor_cond[don],
             n = a$n[match(don, a$sample)], stringsAsFactors = FALSE)
}
cv_types <- sub("_chromvar_crossregion\\.csv$", "", list.files(file.path(CV, "crossregion"), "_chromvar_crossregion\\.csv$"))
cv_names <- vapply(cv_types, function(s) utils::read.csv(file.path(CV, "crossregion", paste0(s, "_chromvar_crossregion.csv")),
                                                          nrows = 1, stringsAsFactors = FALSE)$cell_type, "")
cv_long <- do.call(rbind, unlist(lapply(REGIONS, function(R) lapply(cv_names, function(ct) cv_one(R, ct))), recursive = FALSE))
cv_reg <- cv_long %>% group_by(region, cell_type) %>%
  summarise(donors_PWH = sum(condition == "PWH"), donors_PWoH = sum(condition == "PWoH"), donors_total = dplyr::n(),
            nuclei_total = sum(n, na.rm = TRUE), nuclei_per_donor_median = median(n, na.rm = TRUE),
            nuclei_per_donor_min = min(n, na.rm = TRUE), nuclei_per_donor_max = max(n, na.rm = TRUE), .groups = "drop")
cv_pool <- cv_long %>% dplyr::filter(paste0(region, ".", donor) %in% rownames(mm)) %>% group_by(cell_type) %>%
  summarise(region = "pooled", samples_total = dplyr::n(),
            donors_PWH = dplyr::n_distinct(donor[condition == "PWH"]), donors_PWoH = dplyr::n_distinct(donor[condition == "PWoH"]),
            donors_total = dplyr::n_distinct(donor), nuclei_total = sum(n, na.rm = TRUE),
            nuclei_per_donor_median = median(n, na.rm = TRUE), nuclei_per_donor_min = min(n, na.rm = TRUE),
            nuclei_per_donor_max = max(n, na.rm = TRUE), .groups = "drop")
xrn <- do.call(rbind, lapply(cv_types, function(s) {
  d <- utils::read.csv(file.path(CV, "crossregion", paste0(s, "_chromvar_crossregion.csv")), nrows = 1, stringsAsFactors = FALSE)
  data.frame(cell_type = d$cell_type, n_samples = d$n_samples, n_donors = d$n_donors, rho = d$rho)
}))
cv_pool <- dplyr::left_join(cv_pool, xrn, by = "cell_type")
cv_pool$check_vs_DE_table <- ifelse(cv_pool$samples_total == cv_pool$n_samples & cv_pool$donors_total == cv_pool$n_donors, "match", "DIFFERS")
if (any(cv_pool$check_vs_DE_table != "match"))
  message("[leg] WARNING pooled chromVAR counts differ for: ", paste(cv_pool$cell_type[cv_pool$check_vs_DE_table != "match"], collapse = ", "))
cv_all <- dplyr::bind_rows(cv_reg, cv_pool)
FIG2D <- c("Astro", "Micro", "OPC", "Endo"); FIG2F <- c("Astro", "Micro", "Endo")
add(cv_all[cv_all$cell_type %in% FIG2D, ], "Fig 2", "D", "chromVAR motif activity, donor-mean limma (per region) and pooled cross-region model",
    "donor (per region); region x donor with donor block (pooled)", ">= 20 nuclei per donor")
add(cv_all[cv_all$cell_type %in% FIG2F, ], "Fig 2", "F (y-axis)", "chromVAR pooled centered model", "region x donor; donor block", ">= 20 nuclei")
add(cv_all, "Fig 2", "F supplement", "chromVAR, all cell types", "region x donor; donor block", ">= 20 nuclei")
add(cv_all[cv_all$cell_type == "Astro", ], "Fig S11", "a-b", "astrocyte chromVAR (heatmap = region x condition means; marker = pooled donor-aware FDR)",
    "donor", ">= 20 nuclei")

## ---- 5. microglia signatures and modules ------------------------------------------
sig <- function(dir, lab) do.call(rbind, lapply(REGIONS, function(R) {
  s <- utils::read.csv(file.path(DE, "microglia_signatures", dir, paste0(R, "_scores.csv")), stringsAsFactors = FALSE)
  data.frame(region = R, cell_type = "Micro", donors_PWH = sum(s$condition == "PWH"), donors_PWoH = sum(s$condition == "PWoH"),
             donors_total = nrow(s), nuclei_total = sum(s$n_cells), nuclei_per_donor_median = median(s$n_cells),
             nuclei_per_donor_min = min(s$n_cells), nuclei_per_donor_max = max(s$n_cells))
}))
add(sig(".", "RNA"), "Fig S10", "LDAM signature claim", "microglia signature score, donor pseudobulk (camera / cross-region mixed model)",
    "donor", ">= 10 microglia")
add(sig("atac_gene_activity", "ATAC"), "Response letter", "LDAM in chromatin", "microglia signature score, gene activity",
    "donor", ">= 10 microglia (ATAC)")

mod <- do.call(rbind, lapply(c("CTR", "HIV"), function(g) {
  o <- readRDS(file.path(OUT_ROOT, "microglia_module_analysis", "input", sprintf("scorch_microglia_%s.rds", g)))
  md <- o@meta.data
  dcol <- intersect(c("sample_ids", "orig.ident", "sample"), colnames(md))[1]
  if (is.na(dcol)) stop("no donor column in the module object; columns: ", paste(colnames(md), collapse = ", "))
  rcol <- intersect(c("region", "Region"), colnames(md))[1]
  reg <- if (is.na(rcol)) sub("^.*-([A-Z]+)$", "\\1", rownames(md)) else as.character(md[[rcol]])
  data.frame(region = reg, donor = as.character(md[[dcol]]), condition = cond_chr(g), stringsAsFactors = FALSE)
}))
mod_s <- mod %>% group_by(region, donor, condition) %>% summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(region) %>% summarise(cell_type = "Micro", donors_PWH = sum(condition == "PWH"), donors_PWoH = sum(condition == "PWoH"),
    donors_total = dplyr::n(), nuclei_total = sum(n), nuclei_per_donor_median = median(n),
    nuclei_per_donor_min = min(n), nuclei_per_donor_max = max(n), .groups = "drop")
mod_all <- mod %>% group_by(donor, condition) %>% summarise(n = dplyr::n(), .groups = "drop") %>%
  summarise(region = "all", cell_type = "Micro", donors_PWH = sum(condition == "PWH"), donors_PWoH = sum(condition == "PWoH"),
    donors_total = dplyr::n(), nuclei_total = sum(n), nuclei_per_donor_median = median(n),
    nuclei_per_donor_min = min(n), nuclei_per_donor_max = max(n))
add(dplyr::bind_rows(mod_s, mod_all), "Fig 5", "d-g", "gene modules (gene-gene correlation over nuclei; descriptive, no PWH-vs-PWoH test)",
    "nucleus", "all microglia in the module-analysis input objects")

## ---- 6. astrocytes ---------------------------------------------------------------
comp <- utils::read.csv(file.path(AS, "astro_per_donor_counts_no1040vst.csv"), check.names = FALSE, stringsAsFactors = FALSE)
comp$condition <- cond_chr(comp$HIV)
add(data.frame(region = "all (regions pooled per donor)", cell_type = "Astro", donors_PWH = sum(comp$condition == "PWH"),
               donors_PWoH = sum(comp$condition == "PWoH"), donors_total = nrow(comp), nuclei_total = sum(comp$total),
               nuclei_per_donor_median = median(comp$total), nuclei_per_donor_min = min(comp$total),
               nuclei_per_donor_max = max(comp$total)),
    "Fig S12", "a-b + composition test", "subtype fraction per donor, quasibinomial GLM ~ HIV + Age + Sex + Race + PMI",
    "donor", "every donor with astrocytes (zero counts kept)")

paired <- function(f, lab_col) {
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  cnd <- if (!is.null(d[[lab_col]])) cond_chr(d[[lab_col]]) else cond_chr(d$condition)
  u <- d$usable %in% TRUE
  data.frame(region = "all (regions pooled per donor)", cell_type = "Reactive-PWH vs other astrocytes",
             donors_PWH = sum(u & cnd == "PWH"), donors_PWoH = sum(u & cnd == "PWoH"), donors_total = sum(u),
             nuclei_total = sum(d$n_subtype[u]), nuclei_per_donor_median = median(d$n_subtype[u]),
             nuclei_per_donor_min = min(d$n_subtype[u]), nuclei_per_donor_max = max(d$n_subtype[u]),
             donors_excluded = sum(!u))
}
add(paired(file.path(AS, "astro_Reactive_PWH_donor_counts_no1040vst.csv"), "condition_label"),
    "Fig 6", "e", "paired within-donor pseudobulk edgeR ~ donor + half (nuclei columns = Reactive-PWH nuclei)",
    "donor (paired)", ">= 10 nuclei in both halves")
add(paired(file.path(AS, "figures", "S14b", "donor_aware", "pooled", "S14b_pooled_donor_counts.csv"), "condition"),
    "Fig S14", "b", "paired within-donor gene-activity pseudobulk, GSEA (nuclei columns = Reactive-PWH ATAC nuclei)",
    "donor (paired)", ">= 10 matched ATAC nuclei in both halves")

sub_de <- utils::read.csv(file.path(AS, "subtype_condition_de", "astro_subtype_PWH_vs_PWoH.csv"), stringsAsFactors = FALSE)
sub_de <- unique(sub_de[, c("cluster_id", "n_samples", "n_donors", "rho")])
logs <- sort(Sys.glob(file.path(SCORCH_REPO, "logs", "astro_subtype_condition_de.*.log")), decreasing = TRUE)
ll <- if (length(logs)) grep("\\[sub-de\\] .*samples, .*donors \\(PWH", readLines(logs[1]), value = TRUE) else character(0)
split <- do.call(rbind, lapply(ll, function(l) {
  m <- regmatches(l, regexec("\\[sub-de\\] ([^:]+): (\\d+) samples, (\\d+) donors \\(PWH (\\d+), PWoH (\\d+)\\)", l))[[1]]
  data.frame(cluster_id = m[2], donors_PWH = as.integer(m[5]), donors_PWoH = as.integer(m[6]))
}))
sub_de <- dplyr::left_join(sub_de, split, by = "cluster_id")
add(data.frame(region = "cross-region", cell_type = sub_de$cluster_id, donors_PWH = sub_de$donors_PWH,
               donors_PWoH = sub_de$donors_PWoH, donors_total = sub_de$n_donors, samples_total = sub_de$n_samples,
               rho = sub_de$rho),
    "Fig S15", "a-b", "PWH vs PWoH within astrocyte subtype, cross-region limma-voom + duplicateCorrelation(donor)",
    "region x donor; donor block", ">= 10 nuclei of the subtype per region x donor")

mst <- utils::read.csv(file.path(OUT_ROOT, "magical", "output", "MAGICAL_run_status.csv"), stringsAsFactors = FALSE)
mst <- mst[mst$cell_type == "Astro", ]
add(data.frame(region = mst$region, cell_type = "Astro", samples_total = mst$samples, nuclei_total = mst$cells),
    "Fig 6", "c", "MAGICAL (candidate genes/peaks from the donor-level DE / DA)", "donor-level pseudobulk inputs",
    "all astrocyte samples of the region")

## ---- write -------------------------------------------------------------------------
LONG <- dplyr::bind_rows(ROWS)
lead <- c("figure", "panel", "analysis", "unit", "inclusion_rule", "region", "cell_type",
          "donors_PWH", "donors_PWoH", "donors_total")
LONG <- LONG[, c(lead, setdiff(names(LONG), lead))]
utils::write.csv(LONG, file.path(OUT, "panel_counts_long.csv"), row.names = FALSE)

rng <- function(x) { x <- x[is.finite(x)]; if (!length(x)) return(NA_character_)
  if (min(x) == max(x)) as.character(min(x)) else sprintf("%s-%s", min(x), max(x)) }
SUMM <- LONG %>% group_by(figure, panel, analysis, unit, inclusion_rule, region) %>%
  ## new names for the range columns: summarise() sees a column it has just
  ## replaced, so reusing "donors_total" would break the median below it
  summarise(cell_types = dplyr::n(),
            median_donors = if (all(is.na(donors_total))) NA_real_ else as.numeric(median(donors_total, na.rm = TRUE)),
            donors_range = rng(donors_total), donors_PWH_range = rng(donors_PWH), donors_PWoH_range = rng(donors_PWoH),
            nuclei_per_donor_median = if (all(is.na(nuclei_per_donor_median))) NA_real_ else median(nuclei_per_donor_median, na.rm = TRUE),
            nuclei_per_donor_range = rng(c(nuclei_per_donor_min, nuclei_per_donor_max)),
            .groups = "drop")
utils::write.csv(SUMM, file.path(OUT, "panel_summary.csv"), row.names = FALSE)

## ---- legend text ------------------------------------------------------------------
g <- function(fig, pan) LONG[LONG$figure == fig & LONG$panel == pan, ]
per_region <- function(d) paste(vapply(unique(c(intersect(c(REGIONS, "cross-region", "pooled"), d$region), d$region)), function(R) {
  x <- d[d$region == R, ]
  if (nrow(x) == 1)
    sprintf("%s %d donors (%d PWH, %d PWoH)%s", R, x$donors_total, x$donors_PWH, x$donors_PWoH,
            if (!is.null(x$nuclei_per_donor_median) && is.finite(x$nuclei_per_donor_median))
              sprintf(", median %s nuclei per donor (range %s)", format(round(x$nuclei_per_donor_median), big.mark = ","),
                      rng(c(x$nuclei_per_donor_min, x$nuclei_per_donor_max))) else "")
  else sprintf("%s: %d cell types, %s donors per cell type (median %g; PWH %s, PWoH %s), median %s nuclei per donor",
               R, nrow(x), rng(x$donors_total), median(x$donors_total), rng(x$donors_PWH), rng(x$donors_PWoH),
               format(round(median(x$nuclei_per_donor_median, na.rm = TRUE)), big.mark = ","))
}, ""), collapse = "; ")
ct_line <- function(d) paste(sprintf("%s %d (%d PWH, %d PWoH)", d$cell_type, d$donors_total, d$donors_PWH, d$donors_PWoH), collapse = ", ")
rho_txt <- function(d) paste(sprintf("%s %.2f", d$cell_type, d$rho), collapse = ", ")

L <- c("# Figure legend sample sizes - drafted from the saved analysis objects",
       sprintf("Generated %s by revision/tools/legend_donor_counts.R. Every number comes from panel_counts_long.csv;", format(Sys.Date())),
       "per-cell-type counts too long for a legend belong in a supplementary table built from that file.", "",
       "Wording convention: *donor* = one person; a donor contributes one library per region.", "",
       "## Methods sentence (unit of inference)",
       paste0("For every differential test the unit of inference is the donor. Nuclei were summed per donor (pseudobulk) ",
              "within each cell type; a donor contributed to a cell type's test only if it had at least 10 nuclei of that ",
              "type (20 for chromVAR). Regional models have one sample per donor. In the cross-region models each donor ",
              "contributes up to three region samples; these were modelled as correlated observations of the same person with ",
              "limma duplicateCorrelation (block = donor), so the donor, not the region sample, carries the replication."), "",
       "## Fig 1", paste0("Cohort after QC (snRNA): ", per_region(g("Fig 1", "a-f")), "."), "",
       "## Fig 2",
       paste0("A (snATAC cohort): ", per_region(g("Fig 2", "A")), "."),
       paste0("B, differential accessibility (pseudobulk edgeR, >= 10 nuclei per donor): ", per_region(g("Fig 2", "B")), "."),
       paste0("C, E (genes tested in both RNA and gene activity; gene-activity donors): ", per_region(g("Fig 2", "C, E")), "."),
       paste0("D, chromVAR (>= 20 nuclei per donor), pooled model: ", ct_line(g("Fig 2", "D")[g("Fig 2", "D")$region == "pooled", ]),
              ". Per region: ", per_region(g("Fig 2", "D")[g("Fig 2", "D")$region %in% REGIONS, ]), "."),
       paste0("F, x-axis (cross-region RNA DE): ", ct_line(g("Fig 2", "F (x-axis)")), "; y-axis (pooled chromVAR): ",
              ct_line(g("Fig 2", "F (y-axis)")[g("Fig 2", "F (y-axis)")$region == "pooled", ]), "."), "",
       "## Fig 3 and Fig 4 (regional DE; donors per cell type)", paste0(per_region(g("Fig 3", "a-c")), "."), "",
       "## Fig 5 (microglia)",
       paste0("a, PCA (samples with > 20 microglia): ", per_region(g("Fig 5", "a")), "."),
       paste0("b, heatmap columns (> 30 microglia): ", per_region(g("Fig 5", "b")), "."),
       paste0("b-c, cross-region donor-aware DE: ", with(g("Fig 5", "c"), sprintf("%d region samples from %d donors (%d PWH, %d PWoH); intra-donor correlation %.2f",
              samples_total, donors_total, donors_PWH, donors_PWoH, rho)), "."),
       paste0("d-g, modules (nuclei; descriptive): ", per_region(g("Fig 5", "d-g")), "."),
       "h, feature selection: NOT COVERED - the training data for this panel is not among the rerun outputs; state its unit from the original code.", "",
       "## Fig 6 (astrocytes)",
       paste0("a, heatmap columns (> 30 astrocytes): ", per_region(g("Fig 6", "a")), "."),
       paste0("a-b, cross-region donor-aware DE: ", with(g("Fig 6", "b"), sprintf("%d region samples from %d donors (%d PWH, %d PWoH); intra-donor correlation %.2f",
              samples_total, donors_total, donors_PWH, donors_PWoH, rho)), "."),
       paste0("c, MAGICAL inputs: ", paste(sprintf("%s %d samples", g("Fig 6", "c")$region, g("Fig 6", "c")$samples_total), collapse = ", "), "."),
       paste0("e, paired within-donor test (>= 10 nuclei in both Reactive-PWH and the other astrocytes): ",
              with(g("Fig 6", "e"), sprintf("%d donors (%d PWH, %d PWoH); median %g Reactive-PWH nuclei per donor (range %s); %d donors did not qualify",
                   donors_total, donors_PWH, donors_PWoH, nuclei_per_donor_median, rng(c(nuclei_per_donor_min, nuclei_per_donor_max)), donors_excluded)), "."), "",
       "## Supplementary",
       paste0("S6 (Endo; violins show nuclei, the test is the pseudobulk DE): ", per_region(g("Fig S6", "")), "."),
       paste0("S7 (cross-region, per cell type): ", paste(sprintf("%s %d samples/%d donors (%d PWH, %d PWoH)", g("Fig S7", "a-c")$cell_type,
              g("Fig S7", "a-c")$samples_total, g("Fig S7", "a-c")$donors_total, g("Fig S7", "a-c")$donors_PWH, g("Fig S7", "a-c")$donors_PWoH), collapse = "; "),
              ". Intra-donor correlation: ", rho_txt(g("Fig S7", "a-c")), "."),
       paste0("S8a (DA peaks): ", per_region(g("Fig S8", "a")), "."),
       paste0("S8b (gene activity): ", per_region(g("Fig S8", "b")), "."),
       paste0("S9a (RNA, > 50 microglia): ", per_region(g("Fig S9", "a")), "."),
       paste0("S9b (gene activity, > 30 microglia): ", per_region(g("Fig S9", "b")), "."),
       paste0("S10 (Micro regional DE, >= 10 nuclei): ", per_region(g("Fig S10", "")), ". Signature scores (>= 10 microglia): ",
              per_region(g("Fig S10", "LDAM signature claim")), "."),
       paste0("S11 (astrocyte chromVAR, >= 20 nuclei): ", per_region(g("Fig S11", "a-b")), "."),
       paste0("S12 (composition; one observation per donor, regions pooled, zero counts kept): ",
              per_region(g("Fig S12", "a-b + composition test")), "."),
       paste0("S14b (paired gene-activity test, regions pooled): ",
              with(g("Fig S14", "b"), sprintf("%d donors (%d PWH, %d PWoH); median %g Reactive-PWH nuclei per donor (range %s); %d donors did not qualify",
                   donors_total, donors_PWH, donors_PWoH, nuclei_per_donor_median, rng(c(nuclei_per_donor_min, nuclei_per_donor_max)), donors_excluded)), "."),
       paste0("S15 (within-subtype PWH vs PWoH, cross-region): ", paste(with(g("Fig S15", "a-b"),
              sprintf("%s %d region samples from %d donors (%d PWH, %d PWoH)", cell_type, samples_total, donors_total, donors_PWH, donors_PWoH)), collapse = "; "), "."),
       "", "Not covered here: RNAscope (S13, S17) and NICHES (S16).")
writeLines(L, file.path(OUT, "legend_text.md"))
message("[leg] wrote ", OUT, " (", nrow(LONG), " rows)")
message("[leg] done")
