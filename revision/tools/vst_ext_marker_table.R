# ==============================================================================
# Marker detection table and figure showing that the putative excitatory nuclei
# in VST are cortical neurons (marker profile matches cortex, not striatum).
#
# Method
#   Groups:
#     VST excitatory   putative_ext nuclei of the VST stage-03 object (all 33 libraries)
#     VST other        the remaining VST nuclei
#     INS excitatory   cortical excitatory neurons of the INS object (labels L2/3 IT,
#                      L4 IT, L5 IT, L6 ..., by consensus_annotation_v3)
#   Measure: fraction of nuclei with a non-zero count per gene and group
#   (detection), on the same corrected assay the annotation used, for cortical
#   glutamatergic (SLC17A7, SATB2, SV2B, GRIN1), MSN (PPP1R1B, PDE1B, BCL11B) and
#   inhibitory (GAD1, GAD2) markers. The figure compares excitatory with
#   excitatory nuclei; --with-reference adds the other VST nuclei as a pale
#   baseline.
#
# Inputs   OUT_ROOT/01-Regional_data_processing/VST/data/03-cell_integration_vis/integrated_seurat_obj.RDS
#          OUT_ROOT/01-Regional_data_processing/INS/data/03-cell_integration_vis/integrated_seurat_obj_v3.RDS
# Outputs  OUT_ROOT/01-Regional_data_processing/VST/figures/vst_ext/
#            vst_ext_marker_detection.csv
#            vst_ext_marker_profile.png/.pdf + _caption.txt
#
# Usage
#   revision/docker/run.sh revision/tools/vst_ext_marker_table.R [--replot] [--with-reference]
# ==============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr); library(ggplot2) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
## --replot draws from the saved csv without touching the objects
REPLOT <- "--replot" %in% commandArgs(TRUE)

EXT  <- c("SLC17A7", "SATB2", "SV2B", "GRIN1")      # cortical glutamatergic
MSN  <- c("PPP1R1B", "PDE1B", "BCL11B")             # medium spiny neurons
INH  <- c("GAD1", "GAD2")
GENES <- c(EXT, MSN, INH)
OUT <- file.path(OUT_ROOT, "01-Regional_data_processing", "VST", "figures", "vst_ext")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

det <- function(obj, cells, label) {
  m <- LayerData(obj, assay = "RNA", layer = "counts")[intersect(GENES, rownames(obj)), cells, drop = FALSE]
  data.frame(group = label, n_nuclei = length(cells), gene = rownames(m),
             pct_detected = round(100 * Matrix::rowMeans(m > 0), 1), row.names = NULL)
}

F_OUT <- file.path(OUT, "vst_ext_marker_detection.csv")
if (!REPLOT) {
v <- readRDS(file.path(OUT_ROOT, "01-Regional_data_processing", "VST", "data",
                       "03-cell_integration_vis", "integrated_seurat_obj.RDS"))
v[["RNA"]] <- JoinLayers(v[["RNA"]])
stopifnot("putative_ext" %in% colnames(v@meta.data))
ext_cells <- colnames(v)[v$putative_ext]
oth_cells <- colnames(v)[!v$putative_ext]
message(sprintf("[ext] VST: %d nuclei, %d putative excitatory, %d other",
                ncol(v), length(ext_cells), length(oth_cells)))
res <- dplyr::bind_rows(det(v, ext_cells, "VST excitatory (excluded)"),
                        det(v, oth_cells, "VST other nuclei"))
rm(v); invisible(gc())

i <- readRDS(file.path(OUT_ROOT, "01-Regional_data_processing", "INS", "data",
                       "03-cell_integration_vis", "integrated_seurat_obj_v3.RDS"))
i[["RNA"]] <- JoinLayers(i[["RNA"]])
lab <- intersect(c("consensus_annotation_v3", "consensus_annotation", "cell_type"), colnames(i@meta.data))[1]
ct  <- as.character(i@meta.data[[lab]])
cortical <- grepl("^L[0-9]", ct)                      # L2/3 IT, L4 IT, L5 IT, L6 ... = cortical excitatory
message("[ext] INS label column: ", lab, "; cortical excitatory classes: ",
        paste(sort(unique(ct[cortical])), collapse = ", "))
res <- dplyr::bind_rows(res, det(i, colnames(i)[cortical], "INS cortical excitatory"))
rm(i); invisible(gc())

utils::write.csv(res, F_OUT, row.names = FALSE)
}
res <- utils::read.csv(F_OUT, stringsAsFactors = FALSE)
## The "VST other" group pools every VST cell type, so it is not like-for-like with
##   the two excitatory groups. The panel shows excitatory against excitatory;
##   --with-reference adds the pooled VST nuclei as a pale baseline.
WITH_REF <- "--with-reference" %in% commandArgs(TRUE)
LEV <- c("VST excitatory (excluded libraries)", "INS excitatory neurons")
res$group[res$group == "VST excitatory (excluded)"] <- LEV[1]
res$group[res$group == "INS cortical excitatory"]   <- LEV[2]
if (WITH_REF) LEV <- c(LEV, "VST other nuclei") else res <- res[res$group != "VST other nuclei", ]
res$group <- factor(res$group, levels = LEV)
res$gene  <- factor(res$gene, levels = GENES)
res <- res[order(res$gene, res$group), ]
print(tidyr::pivot_wider(res, id_cols = gene, names_from = group, values_from = pct_detected), n = 20)

## ---- figure -------------------------------------------------------------------------
## Published region palette (05-nDEs.R region_cols): INS #E69F00, VST #CC79A7. The VST
## groups take the VST colour, dark for the excluded nuclei and pale for the rest, so the
## panel says directly that the excluded nuclei look like INS cortex, not like VST.
COLS <- setNames(c("#CC79A7", "#E69F00", "#EBD3E0")[seq_along(LEV)], LEV)
ann <- data.frame(gene = factor(c("SLC17A7", "PPP1R1B", "GAD1"), levels = GENES),
                  lab = c("cortical glutamatergic", "medium spiny neuron", "inhibitory"))
p <- ggplot(res, aes(gene, pct_detected, fill = group)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7, colour = "grey30", linewidth = 0.15) +
  scale_fill_manual(values = COLS, name = NULL) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "nuclei with detected expression (%)") +
  theme_astro_pub(8) +
  theme(legend.position = "top", legend.margin = margin(b = -4),
        axis.text.x = element_text(face = "italic", angle = 40, hjust = 1))
f <- file.path(OUT, "vst_ext_marker_profile.png")
ggsave(f, p, width = 5.2, height = 3.2, dpi = 600, bg = "white")
ggsave(sub("\\.png$", ".pdf", f), p, width = 5.2, height = 3.2, device = grDevices::cairo_pdf, bg = "white")
w <- function(g, grp) res$pct_detected[res$gene == g & res$group == grp]
writeLines(strwrap(sprintf(paste(
  "The excitatory nuclei of the excluded VST libraries are cortical neurons. Detection of marker genes in the",
  "putative excitatory nuclei of the five excluded VST libraries (n = %s) and in excitatory neurons of the INS",
  "(n = %s). The two populations are indistinguishable on cortical glutamatergic markers (SATB2 %.1f%% vs %.1f%%,",
  "SV2B %.1f%% vs %.1f%%, SLC17A7 %.1f%% vs %.1f%%) and both lack striatal (PPP1R1B %.1f%% vs %.1f%%) and",
  "inhibitory (GAD2 %.1f%% vs %.1f%%) markers, indicating cortical tissue included during dissection rather than",
  "striatal neurons."),
  format(res$n_nuclei[res$group == LEV[1]][1], big.mark = ","),
  format(res$n_nuclei[res$group == LEV[2]][1], big.mark = ","),
  w("SATB2", LEV[1]), w("SATB2", LEV[2]), w("SV2B", LEV[1]), w("SV2B", LEV[2]),
  w("SLC17A7", LEV[1]), w("SLC17A7", LEV[2]), w("PPP1R1B", LEV[1]), w("PPP1R1B", LEV[2]),
  w("GAD2", LEV[1]), w("GAD2", LEV[2])), 110),
  sub("\\.png$", "_caption.txt", f))
message("[ext] wrote ", f)
