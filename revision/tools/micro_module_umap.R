# ==============================================================================
# Fig 5e microglia module scores drawn on a UMAP instead of the tSNE.
#
# Method: same scoring and styling as Ruiqi's module-score figure, only the
#   embedding differs: AddModuleScore(features = modules, ctrl = 5,
#   name = "Module"), modules in split() order (hiv_down_module, hiv_up_module1,
#   hiv_up_module2 = Module1-3), FeaturePlot order = TRUE, RdBu ramp, limits
#   -0.5..0.75 squished, NoAxes, PWH row over PWoH row, 13 x 8 in.
#   UMAP: RunUMAP on the objects' own "harmony" reduction, dims 1:30 (the dims
#   their tSNE used), seed 3. A tSNE panel set is drawn from the same scores as
#   a check against the published tSNE figure.
#
# Inputs: OUT_ROOT/microglia_module_analysis/input/scorch_microglia_{HIV,CTR}.rds
#         (verbatim, md5-checked copies of the module analysis's per-condition
#         microglia objects) and selected_gene_modules.csv
# Output: OUT_ROOT/microglia_module_analysis/umap/
#           fig5e_modulescore_umap.png/.pdf
#           fig5e_modulescore_tsne_check.png
#           umap_coords_<PWH|PWoH>.csv            cached embedding
#
# Usage: revision/docker/run.sh revision/tools/micro_module_umap.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(RColorBrewer)
})
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

BASE <- file.path(OUT_ROOT, "microglia_module_analysis")
OUT  <- file.path(BASE, "umap")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

df <- utils::read.csv(file.path(BASE, "selected_gene_modules.csv"), stringsAsFactors = FALSE)
selected_modules <- split(df$gene, df$module)
message("[modumap] modules: ", paste(sprintf("%s (%d)", names(selected_modules), lengths(selected_modules)),
                                     collapse = ", "))

objs <- list(PWH  = readRDS(file.path(BASE, "input", "scorch_microglia_HIV.rds")),
             PWoH = readRDS(file.path(BASE, "input", "scorch_microglia_CTR.rds")))

for (cond in names(objs)) {
  o <- objs[[cond]]
  stopifnot("harmony" %in% Reductions(o))
  o <- AddModuleScore(o, features = selected_modules, ctrl = 5, name = "Module")
  o <- RunUMAP(o, reduction = "harmony", dims = 1:30, seed.use = 3, verbose = FALSE)
  emb <- Embeddings(o, "umap")
  utils::write.csv(data.frame(cell = rownames(emb), emb, o@meta.data[rownames(emb), paste0("Module", 1:3)]),
                   file.path(OUT, sprintf("umap_coords_%s.csv", cond)), row.names = FALSE)
  message(sprintf("[modumap] %s: %d nuclei, %d donors", cond, ncol(o), dplyr::n_distinct(o$sample_ids)))
  objs[[cond]] <- o
}

plot_module_scores <- function(srat, reduction) {
  pl <- lapply(paste0("Module", 1:3), function(module) {
    FeaturePlot(srat, reduction = reduction, features = module, order = TRUE) +
      scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),
                            limits = c(-0.5, 0.75), oob = scales::squish) +
      NoAxes() + ggtitle(module)
  })
  wrap_plots(pl, nrow = 1) + plot_layout(guides = "collect")
}

for (red in c("umap", "tsne")) {
  p <- plot_module_scores(objs$PWH, red) / plot_module_scores(objs$PWoH, red)
  f <- if (red == "umap") file.path(OUT, "fig5e_modulescore_umap.png")
       else file.path(OUT, "fig5e_modulescore_tsne_check.png")
  ggsave(f, p, width = 13, height = 8, dpi = 300, bg = "white")
  if (red == "umap") ggsave(sub("\\.png$", ".pdf", f), p, width = 13, height = 8,
                            device = grDevices::cairo_pdf, bg = "white")
  message("[modumap] wrote ", basename(f))
}
message("[modumap] done -> ", OUT)
