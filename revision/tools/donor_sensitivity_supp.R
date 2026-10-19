# Supplementary figure and table for the outlier-donor sensitivity analyses (exclusion of donors 1039 and 1181, and leave-one-donor-out).
#
# Method: answers the reviewer requests to exclude the HIVE donor (1039) and the
#   donor with multiple reported pathologies (1181), and for leave-one-donor-out
#   analyses. Nothing is refitted; every number is read from existing outputs:
#   a. DE, all donors vs excluding 1039 + 1181: gene_level_fates.csv from
#      03-DE_analysis/11-compare_sensitivity.R (DEG = p_adj.glb < 0.05 and
#      |logFC| > 0.5, that script's default, not the paper's 0.25 + >5% rule).
#      --cmp-suffix _negF reads the pair with the negative-F fallback applied to
#      both runs (revision/tools/sensitivity_negF_refit.R). r is Pearson over
#      every gene tested in both runs (r = 0.96 / 0.96 / 0.90 for INS / PFC / VST).
#   b. microglial signatures, LODO + exclusion (revision/tools/donor_lodo_microglia.R)
#   c. Reactive-PWH astrocyte composition, LODO + exclusion
#      (revision/tools/donor_lodo_astro.R, _no1040vst; drawn as in
#      revision/tools/astro_lodo_forest.R)
#   Donors actually removed differ by analysis, and the caption says so:
#   DE           INS/PFC 31 -> 29 (1039, 1181); VST 28 -> 27 (1039's VST library is
#                already out under the contamination filter, so only 1181 drops)
#   microglia    INS drops 1039 only (1181 has < 10 INS microglia); VST 1181 only
#
# Inputs:  OUT_ROOT/03-DE_analysis/comparisons/donor_sensitivity/<R>_noHIVE_noPath[<cmp-suffix>]/gene_level_fates.csv
#          OUT_ROOT/03-DE_analysis/donor_sensitivity/{microglia_lodo_summary,<R>_microglia_lodo_allfits}[_modules].csv
#          OUT_ROOT/03-DE_analysis/donor_sensitivity/astro_composition_lodo_{summary,allfits}_no1040vst.csv
#          OUT_ROOT/06-astrocyte_analysis/astro_per_donor_counts_no1040vst.csv
# Outputs: OUT_ROOT/03-DE_analysis/donor_sensitivity/supp/
#            supp_donor_sensitivity.png/.pdf + _caption.txt
#            supp_donor_sensitivity_table.csv (+ .md)
#
# Usage: revision/docker/run.sh revision/tools/donor_sensitivity_supp.R [--cmp-suffix _negF]

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

DE  <- file.path(OUT_ROOT, "03-DE_analysis")
DS  <- file.path(DE, "donor_sensitivity")
CMP <- file.path(DE, "comparisons", "donor_sensitivity")
AS  <- file.path(OUT_ROOT, "06-astrocyte_analysis")
OUT <- file.path(DS, "supp"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
REGIONS <- c("INS", "PFC", "VST")
args   <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
CMP_SFX <- getarg("--cmp-suffix", "")
N_DONORS <- list(INS = c(31, 29), PFC = c(31, 29), VST = c(28, 27))   # donors per DE run (all / excluding), from the runs' run-info output
EXCL_OUT <- "#D55E00"   # "without 1039 + 1181" marker, distinct from PWH orange

## ---- a: gene-level DE, primary vs excluding both ---------------------------------------------
de_rows <- list(); sc <- list()
for (R in REGIONS) {
  g <- data.table::fread(file.path(CMP, paste0(R, "_noHIVE_noPath", CMP_SFX), "gene_level_fates.csv"),
                         na.strings = c("NA", ""))
  shared <- g[!is.na(logFC_a) & !is.na(logFC_b)]
  r <- cor(shared$logFC_a, shared$logFC_b)
  prim <- g[is_deg_a == TRUE]
  rt <- prim[!is.na(logFC_b)]
  same <- sum(sign(rt$logFC_a) == sign(rt$logFC_b))
  de_rows[[R]] <- data.frame(
    block = "Differential expression", region = R, finding = "All DEGs (per cell type)",
    donors = sprintf("%d -> %d", N_DONORS[[R]][1], N_DONORS[[R]][2]),
    primary = sprintf("%s DEGs", format(nrow(prim), big.mark = ",")),
    excluding_1039_1181 = sprintf("%s DEGs; %.1f%% of primary DEGs still significant",
                                  format(sum(g$is_deg_b), big.mark = ","),
                                  100 * sum(g$fate == "retained") / nrow(prim)),
    lodo = "not run (per-gene LODO not performed)",
    detail = sprintf("logFC r = %.3f over %s genes; same sign %d/%d re-tested primary DEGs (%.2f%%); %d significant in both with opposite sign",
                     r, format(nrow(shared), big.mark = ","), same, nrow(rt), 100 * same / nrow(rt),
                     sum(g$fate == "flipped")),
    stringsAsFactors = FALSE)
  sc[[R]] <- list(d = shared, r = r, deg = shared[is_deg_a == TRUE | is_deg_b == TRUE],
                  pct = 100 * same / nrow(rt))
  message(sprintf("[sens] %s: r = %.3f, sign %d/%d", R, r, same, nrow(rt)))
}
lim <- max(abs(unlist(lapply(sc, function(x) quantile(abs(c(x$d$logFC_a, x$d$logFC_b)), 0.999)))))
lim <- ceiling(lim)
pa <- lapply(REGIONS, function(R) {
  x <- sc[[R]]
  ## bins over the DRAWN window only: a few |logFC| > 15 genes in VST otherwise
  ## stretch the bin grid and every bin becomes a large block
  inwin <- x$d[abs(logFC_a) <= lim & abs(logFC_b) <= lim]
  ggplot(x$d, aes(logFC_a, logFC_b)) +
    ## density of every tested gene (~150-300k per region) as bins, so the PDF stays light
    geom_bin2d(data = inwin, binwidth = c(2 * lim / 90, 2 * lim / 90),
               aes(fill = after_stat(log10(count)))) +
    scale_fill_gradient(low = "grey88", high = "grey35", guide = "none") +
    geom_point(data = x$deg, colour = PWH_COLOURS[["PWH"]], size = 0.15, alpha = 0.35) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed") +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey50") +
    annotate("text", x = -lim * 0.95, y = lim * 0.95, hjust = 0, vjust = 1, size = 2.3,
             label = sprintf("r = %.2f\nsame sign %.1f%%", x$r, floor(x$pct * 10) / 10)) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(title = R, x = "log2FC, all donors", y = "log2FC, without 1039 + 1181") +
    theme_astro_pub(7) + theme(plot.title = element_text(hjust = 0.5, size = 8))
})

## ---- b: microglial signatures, LODO ---------------------------------------------------------------
SIG_LAB <- c(LDAM_Haney2024 = "LDAM (Haney 2024)",
             DAM_stage1_KerenShaul2017 = "DAM stage 1 (Keren-Shaul 2017)",
             IRM_core = "Interferon-response (IRM)",
             Homeostatic = "Homeostatic (down)",
             MG4_lipid_processing_Sun2023 = "Lipid processing, MG4 (Sun 2023)",
             ## the paper's own CS-CORE modules (Fig 5), revision/tools/donor_lodo_microglia.R --sets modules
             Module1_homeostatic_down = "Module 1, homeostatic (down)",
             Module2_inflammatory_up  = "Module 2, inflammatory",
             Module3_hypoxia_lipid_up = "Module 3, hypoxia/lipid")
ms <- bind_rows(utils::read.csv(file.path(DS, "microglia_lodo_summary.csv"), stringsAsFactors = FALSE),
                utils::read.csv(file.path(DS, "microglia_lodo_summary_modules.csv"), stringsAsFactors = FALSE))
## show every signature significant in the primary analysis of at least one region
keep_sig <- intersect(names(SIG_LAB), unique(ms$signature[ms$fdr_full < 0.05]))
mf <- bind_rows(lapply(REGIONS, function(R)
  bind_rows(utils::read.csv(file.path(DS, paste0(R, "_microglia_lodo_allfits.csv")), stringsAsFactors = FALSE),
            utils::read.csv(file.path(DS, paste0(R, "_microglia_lodo_allfits_modules.csv")), stringsAsFactors = FALSE)) %>%
    mutate(region = R))) %>% filter(signature %in% keep_sig)
ms_k <- ms %>% filter(signature %in% keep_sig)
lev <- rev(SIG_LAB[keep_sig])
mf$sig <- factor(SIG_LAB[mf$signature], levels = lev)
refs <- bind_rows(
  ms_k %>% transmute(region, sig = SIG_LAB[signature], fdr = fdr_full, what = "All donors"),
  ms_k %>% transmute(region, sig = SIG_LAB[signature], fdr = fdr_excl_pathology,
                     what = "Without 1039 + 1181")) %>%
  mutate(sig = factor(sig, levels = lev))
xl <- -log10(0.05)
XCAP <- 8   # module FDRs reach 1e-41; beyond the cap points are drawn AT the cap
mf$x <- pmin(-log10(mf$fdr), XCAP); refs$x <- pmin(-log10(refs$fdr), XCAP)
pb <- ggplot(mf, aes(x, sig)) +
  geom_vline(xintercept = xl, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_point(position = position_jitter(height = 0.12, width = 0, seed = 1),
             colour = "grey60", size = 0.6, alpha = 0.8) +
  geom_point(data = refs, aes(shape = what, colour = what, fill = what), size = 1.9, stroke = 0.5) +
  scale_shape_manual(values = c(`All donors` = 23, `Without 1039 + 1181` = 24), name = NULL) +
  scale_colour_manual(values = c(`All donors` = "black", `Without 1039 + 1181` = "black"), name = NULL) +
  scale_fill_manual(values = c(`All donors` = "white", `Without 1039 + 1181` = EXCL_OUT), name = NULL) +
  facet_wrap(~ region, nrow = 1) +
  scale_x_continuous(limits = c(0, XCAP), breaks = c(0, 2, 4, 6, 8), labels = c("0", "2", "4", "6", "\u22658")) +
  labs(x = expression(-log[10]~FDR), y = NULL) +
  theme_astro_pub(7) + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))

## ---- c: Reactive-PWH astrocytes, LODO forest (as revision/tools/astro_lodo_forest.R) -------------
ST <- "Reactive-PWH"
af <- utils::read.csv(file.path(DS, "astro_composition_lodo_allfits_no1040vst.csv"), stringsAsFactors = FALSE)
as_ <- utils::read.csv(file.path(DS, "astro_composition_lodo_summary_no1040vst.csv"), stringsAsFactors = FALSE)
cnt <- utils::read.csv(file.path(AS, "astro_per_donor_counts_no1040vst.csv"), check.names = FALSE,
                       stringsAsFactors = FALSE)
cond <- setNames(as.character(pwh_relabel(cnt$HIV)), as.character(cnt$donor))
f <- af[af$subtype == ST, ]; s <- as_[as_$subtype == ST, ]
f$group <- cond[as.character(f$dropped)]
f <- f[order(f$odds_ratio), ]
Dc <- rbind(data.frame(label = c("All donors", "Without 1039 + 1181"),
                       odds_ratio = c(s$odds_ratio_full, s$or_excl_pathology),
                       FDR = c(s$fdr_full, s$fdr_excl_pathology), group = c("All donors", "Without 1039 + 1181")),
            data.frame(label = paste0("−", f$dropped), odds_ratio = f$odds_ratio, FDR = f$FDR, group = f$group))
Dc$label <- factor(Dc$label, levels = rev(Dc$label))
Dc$sig <- ifelse(Dc$FDR < 0.05, "FDR < 0.05", "FDR ≥ 0.05")
pc <- ggplot(Dc, aes(odds_ratio, label)) +
  geom_vline(xintercept = 1, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = s$odds_ratio_full, colour = "grey40", linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(colour = group, shape = sig), size = 1.3, stroke = 0.5) +
  scale_colour_manual(values = c(`All donors` = "black", `Without 1039 + 1181` = "black", PWH_COLOURS),
                      breaks = c("PWoH", "PWH"), name = "Dropped donor") +
  scale_shape_manual(values = c(`FDR < 0.05` = 16, `FDR ≥ 0.05` = 1), name = NULL) +
  scale_x_log10(breaks = c(1, 3, 10, 30)) +
  labs(x = "Reactive-PWH odds ratio,\nPWH vs PWoH (log scale)", y = NULL) +
  theme_astro_pub(7) +
  guides(colour = guide_legend(nrow = 1, order = 1, title.position = "top"),
         shape = guide_legend(nrow = 1, order = 2)) +
  theme(axis.text.y = element_text(size = 4.8), legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 6), legend.title = element_text(size = 6.5),
        legend.spacing.x = unit(2, "pt"), legend.margin = margin(0, 0, 0, 0))

## ---- assemble --------------------------------------------------------------------------------------
## the top row is not axis-aligned with b, whose long signature labels would push it right
top <- wrap_plots(pa, nrow = 1)
bottom <- wrap_plots(pb, pc, widths = c(2.7, 1))
## wrap_elements makes the top row one tagged element ("a") and frees its alignment
fig <- wrap_plots(wrap_elements(full = top + plot_annotation(theme = theme())), bottom,
                  ncol = 1, heights = c(1, 1.45)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))
fp <- file.path(OUT, "supp_donor_sensitivity.png")
save_pub(fig, fp, width = 7.2, height = 6.6)

## ---- table ------------------------------------------------------------------------------------------
fmt_p <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x))
mic_rows <- ms_k %>% mutate(
  block = "Microglial signatures (camera, donor-level)",
  finding = paste0(SIG_LAB[signature], ", ", tolower(direction)),
  donors = sprintf("%d; excluded: %s", n_donors_full, gsub("\\s+", ", ", trimws(excluded_donors))),
  primary = paste0("FDR ", fmt_p(fdr_full)),
  excluding_1039_1181 = paste0("FDR ", fmt_p(fdr_excl_pathology)),
  lodo = sprintf("%d/%d refits FDR < 0.05 (FDR %s-%s)", lodo_n_sig, lodo_n_fits, fmt_p(lodo_fdr_min), fmt_p(lodo_fdr_max)),
  detail = paste0("most influential donor: ", most_influential,
                  ifelse(direction_stable, "", "; direction NOT stable across refits"))) %>%
  arrange(match(signature, keep_sig), match(region, REGIONS)) %>%
  select(block, region, finding, donors, primary, excluding_1039_1181, lodo, detail)
ast_rows <- as_ %>% filter(fdr_full < 0.05) %>% mutate(
  block = "Astrocyte subtype abundance (quasibinomial, donor-level)", region = "INS+PFC+VST",
  finding = paste0(subtype, " (", mean_pct_CTR, "% PWoH vs ", mean_pct_HIV, "% PWH)"),
  donors = sprintf("%d; excluded: %s", n_donors_full, excluded_donors),
  primary = sprintf("OR %.1f, FDR %s", odds_ratio_full, fmt_p(fdr_full)),
  excluding_1039_1181 = sprintf("OR %.1f, FDR %s", or_excl_pathology, fmt_p(fdr_excl_pathology)),
  lodo = sprintf("%d/%d refits FDR < 0.05 (OR %.1f-%.1f)", lodo_n_sig, lodo_n_fits, lodo_or_min, lodo_or_max),
  detail = paste0("most influential donor: ", most_influential)) %>%
  select(block, region, finding, donors, primary, excluding_1039_1181, lodo, detail)
tab <- bind_rows(bind_rows(de_rows), mic_rows, ast_rows)
utils::write.csv(tab, file.path(OUT, "supp_donor_sensitivity_table.csv"), row.names = FALSE)
md <- c("| Analysis | Region | Finding | Donors | All donors | Without 1039 + 1181 | Leave-one-donor-out | Notes |",
        "|---|---|---|---|---|---|---|---|",
        apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
writeLines(md, file.path(OUT, "supp_donor_sensitivity_table.md"))
options(width = 250); print(tab[, -1], row.names = FALSE)

## ---- caption -------------------------------------------------------------------------------------
st <- setNames(lapply(REGIONS, function(R) sc[[R]]), REGIONS)
cap <- paste0(
  "Sensitivity of the principal findings to donors 1039 (HIV encephalitis) and 1181 (additional pathologies). ",
  "(a) Pseudobulk log2 fold changes (PWH vs PWoH) for every gene and cell type tested in both the primary analysis ",
  "and the analysis excluding both donors (grey, density; orange, significant in either). Pearson r: ",
  paste(sprintf("%s %.2f", REGIONS, sapply(st, `[[`, "r")), collapse = ", "), ". ",
  "INS and PFC, 31 to 29 donors; VST, 28 to 27 donors (the VST library of 1039 is already excluded in the primary analysis). ",
  "(b) Microglial gene-signature scores tested at the donor level (camera; published signatures and the three ",
  "co-expression modules of Fig. 5, the latter tested separately); grey points, leave-one-donor-out refits; ",
  "open diamond, all donors; filled triangle, excluding both donors (INS excludes 1039 only and VST 1181 only, ",
  "because the other donor contributes fewer than 10 microglia). Dashed line, FDR = 0.05. ",
  "The modules were defined from PWH-vs-PWoH genes in these data, so their full-cohort FDRs are not independent ",
  "evidence; values beyond 8 on the axis are drawn at 8 (exact values in Supplementary Table X). The refits show whether the separation depends on any one donor. ",
  "(c) Reactive-PWH astrocyte abundance, donor-level composition model: odds ratio (PWH vs PWoH) with all donors, ",
  "without both donors, and with each donor left out in turn (point colour, condition of the dropped donor; ",
  "open points, FDR >= 0.05). Dashed line, all-donor odds ratio.")
astro_caption(fp, cap)
message("[sens] done: ", OUT)
