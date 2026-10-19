# ==============================================================================
# Supplementary HIV-alignment figure (HIV-aligned reads per nucleus), remade with
# the CellBender-corrected nuclei and cell type annotation.
#
# Method
#   Same layout as the published figure:
#     rows      one nucleus each, labelled "<cell type>\n<barcode>", grouped by donor
#     panels    HIV DNA (Wei et al.; snATAC) | HIV RNA (Wei et al.) | HIV RNA (Venus),
#               each split by region
#     dots      size = HIV-aligned reads (PCR-collapsed unique counts)
#     inclusion nuclei with >= MIN_READS reads (default 3, i.e. the published > 2)
#               in a Wei et al. method (DNA or RNA); Venus counts are drawn for
#               those nuclei only and cannot admit a nucleus alone
#   Nuclei and labels come from the CellBender-corrected object. A barcode that
#   CellBender did not retain as a nucleus has no cell type and is not drawn; such
#   nuclei from the published figure are listed in viral_figure_not_retained.csv.
#
# Inputs
#   Wei DNA / RNA  OUT_ROOT/07-Viral_Detection/viral_barcode_counts.csv.gz
#                  (from revision/tools/viral_threshold_table.R; hiv_atac /
#                  ho_protocol, the variants the published script read)
#   Venus RNA      RAW_DIR/<run>/venus/detection/<sample>/output.txt
#                  ("<barcode>, <count>"; reproduces the published Venus column
#                  for 1039 PFC: 24 / 20 / 9 / 3)
#   annotations    OUT_ROOT/07-Viral_Detection/annotation_{corrected,published}.csv.gz
# Outputs  OUT_ROOT/07-Viral_Detection/figure/
#           hiv_alignment_cellbender.png/.pdf   the figure
#           viral_figure_rows.csv               every dot drawn
#           viral_figure_not_retained.csv       published-figure nuclei CellBender dropped
#
# Usage
#   revision/docker/run.sh revision/tools/viral_figure.R [--min-reads 3] [--immune-label Immune]
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2) })
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
MIN_READS <- as.integer(getarg("--min-reads", "3"))            # published: > 2
## The published figure called this group "Myeloid"; it is labelled "Immune", as
## in the CellBender annotation. --immune-label overrides.
IMMUNE_LAB <- getarg("--immune-label", "Immune")
CT_LAB <- c(Astro = "Astrocyte", Micro = "Microglia", Immune = IMMUNE_LAB,
            `Other immune` = IMMUNE_LAB, Macrophages = "Macrophage", Oligo = "Oligodendrocyte",
            OPC = "OPC", Endo = "Endothelial")

VD  <- file.path(OUT_ROOT, "07-Viral_Detection")
OUT <- file.path(VD, "figure"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RAW <- RAW_DIR

## ---- cohort and annotations ------------------------------------------------------------------
meta   <- load_sample_meta()
cohort <- unique(toupper(paste(meta$sample_ids, meta$region, sep = "_")))
cond   <- unique(meta[, c("sample_ids", "hiv_jay_v2")])
ann <- lapply(c(corrected = "corrected", published = "published"), function(k)
  data.table::fread(file.path(VD, paste0("annotation_", k, ".csv.gz"))))

## ---- Wei DNA and Wei RNA: the cached per-barcode table --------------------------------------
wei <- data.table::fread(file.path(VD, "viral_barcode_counts.csv.gz"), na.strings = c("NA", ""))
## the table holds every protocol variant; `keep` marks the one the precedence rule chose
wei <- wei[keep == TRUE, .(barcode, donor, region, cell, n,
               method = ifelse(assay == "DNA", "HIV DNA\nWei et al", "HIV RNA\nWei et al"))]

## ---- Venus RNA ---------------------------------------------------------------------------------
vf <- Sys.glob(file.path(RAW, "*", "venus", "detection", "*", "output.txt"))
message("[viral-fig] Venus output files: ", length(vf))
venus <- rbindlist(lapply(vf, function(f) {
  stem  <- basename(dirname(f))
  reg   <- regmatches(toupper(stem), regexpr("INS|PFC|VST", toupper(stem)))
  if (!length(reg) || file.size(f) == 0) return(NULL)
  d <- tryCatch(fread(f, header = FALSE, col.names = c("barcode", "n"), sep = ","),
                error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d[, `:=`(barcode = trimws(barcode), donor = sub("_.*$", "", stem), region = reg, file = f)]
}), fill = TRUE)
venus <- venus[toupper(paste(donor, region, sep = "_")) %in% cohort]
## one output per sample x region; a sample sequenced in two runs would double here
dup <- venus[, uniqueN(file), by = .(donor, region)][V1 > 1]
if (nrow(dup)) stop("Venus: more than one output for ", paste(dup$donor, dup$region, collapse = ", "),
                    " - decide which run before plotting")
venus[, cell := paste0(barcode, "-1-", donor, "-", region)]
venus <- venus[, .(barcode, donor, region, cell, n, method = "HIV RNA\nVenus")]
message("[viral-fig] Venus barcodes in the cohort: ", nrow(venus))

## ---- combine, label, select ----------------------------------------------------------------------
x <- rbind(wei, venus)
## donor IDs are inconsistently capitalised in raw file names (HCtME vs HCTME)
idx <- match(toupper(x$donor), toupper(cond$sample_ids))
x[, donor := ifelse(is.na(idx), donor, cond$sample_ids[idx])]
x[, condition := ifelse(toupper(as.character(cond$hiv_jay_v2[idx])) == "TRUE", "PWH", "PWoH")]
x[, cell := paste0(barcode, "-1-", donor, "-", region)]
x[, ct_cb  := ann$corrected$cell_type[match(cell, ann$corrected$cell)]]
x[, ct_pub := ann$published$cell_type[match(cell, ann$published$cell)]]

## a nucleus enters on a Wei et al. hit (DNA or RNA); Venus is shown alongside but
## cannot admit a nucleus alone, as in the published figure (so 6801066772 INS
## AAATGCCTCGGTAACT, Venus 16 / Wei 0, is not drawn)
hit_any <- unique(x[n >= MIN_READS & grepl("Wei", method)]$cell)
venus_only <- setdiff(unique(x[n >= MIN_READS]$cell), hit_any)
message("[viral-fig] Venus-only nuclei left out: ", paste(venus_only, collapse = ", "))
## nuclei the PUBLISHED figure would have drawn but CellBender did not retain
lost <- unique(x[cell %in% hit_any & !is.na(ct_pub) & is.na(ct_cb),
                 .(cell, donor, region, published_cell_type = ct_pub)])
fwrite(lost, file.path(OUT, "viral_figure_not_retained.csv"))

d <- x[cell %in% hit_any & !is.na(ct_cb) & n >= MIN_READS]
if (!nrow(d)) stop("no CellBender-retained nucleus has >= ", MIN_READS, " reads in any method")
d[, cell_type := ifelse(ct_cb %in% names(CT_LAB), CT_LAB[ct_cb], ct_cb)]
d[, row_lab := paste0(cell_type, "\n", barcode)]
d[, donor_lab := paste0(condition, "\n", donor)]
d[, method := factor(method, levels = c("HIV DNA\nWei et al", "HIV RNA\nWei et al", "HIV RNA\nVenus"))]
REG <- intersect(c("INS", "PFC", "VST"), unique(d$region))
d[, region := factor(region, levels = REG)]
## rows: by donor, then by cell type, then by strongest signal
ord <- d[, .(mx = max(n)), by = .(donor_lab, row_lab, cell_type)][order(donor_lab, cell_type, -mx)]
d[, row_lab := factor(row_lab, levels = rev(unique(ord$row_lab)))]
fwrite(d[order(donor_lab, row_lab)], file.path(OUT, "viral_figure_rows.csv"))
message("[viral-fig] drawn: ", uniqueN(d$cell), " nuclei, ", nrow(d), " dots; ",
        "published-figure nuclei not retained by CellBender: ", nrow(lost))
print(as.data.frame(d[order(donor_lab, row_lab), .(donor, region, barcode, cell_type, method = gsub("\n", " ", method), n)]),
      row.names = FALSE)

## ---- the figure, in the published style ----------------------------------------------------
## every row x region x method gets a faint guide line, as in the original
grid_rows <- CJ(row_lab = levels(d$row_lab), region = REG, method = levels(d$method))
grid_rows <- merge(grid_rows, unique(d[, .(row_lab = as.character(row_lab), donor_lab)]), by = "row_lab")
grid_rows[, `:=`(row_lab = factor(row_lab, levels = levels(d$row_lab)),
                 region = factor(region, levels = REG),
                 method = factor(method, levels = levels(d$method)))]
FILL <- c("HIV DNA\nWei et al" = "#7E4FA0", "HIV RNA\nWei et al" = "#F39120",
          "HIV RNA\nVenus" = "#E8702A")
p <- ggplot(d, aes(region, row_lab)) +
  ## thin tiles on the discrete region axis (a numeric-x segment clashes with it)
  geom_tile(data = grid_rows, aes(region, row_lab), width = 0.9, height = 0.06,
            fill = "grey88", inherit.aes = FALSE) +
  geom_point(aes(size = n, fill = method), shape = 21, colour = "black", stroke = 0.5) +
  scale_fill_manual(values = FILL, guide = "none") +
  scale_size_area(max_size = 9, breaks = c(3, 10, 20), name = "# of reads") +
  facet_grid(donor_lab ~ method, scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        strip.text.x = element_text(size = 15, lineheight = 0.95),
        strip.text.y.left = element_text(size = 13, angle = 0, hjust = 1),
        strip.placement = "outside",
        axis.text.y = element_text(size = 8, lineheight = 0.9, hjust = 1),
        axis.text.x = element_text(size = 13),
        panel.spacing.x = unit(1.4, "lines"), panel.spacing.y = unit(0.9, "lines"),
        legend.position = "right", plot.background = element_rect(fill = "white", colour = NA))
H <- 1.4 + 0.52 * uniqueN(d$row_lab) + 0.35 * uniqueN(d$donor_lab)
W <- 3.2 + 1.35 * length(REG) * 3
f <- file.path(OUT, "hiv_alignment_cellbender.png")
ggsave(f, p, width = W, height = H, dpi = 300, bg = "white")
ggsave(sub("\\.png$", ".pdf", f), p, width = W, height = H, device = grDevices::cairo_pdf, bg = "white")
message(sprintf("[viral-fig] wrote %s (%.1f x %.1f in)", f, W, H))
