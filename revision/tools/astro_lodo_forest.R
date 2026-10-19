# ==============================================================================
# Leave-one-donor-out forest plot for the donor-level abundance test of an
# astrocyte subtype (default Reactive-PWH; reviewer request).
#
# Method
#   Plots refits already made by revision/tools/donor_lodo_astro.R; nothing is
#   refitted. One row per dropped donor: odds ratio (PWH vs PWoH) of the donor-level
#   quasibinomial composition model (~ HIV + Age + Sex + Race + PMI) refitted without
#   that donor, filled when FDR < 0.05 (BH across the seven subtypes). Top row = all
#   donors; second = pathology donors 1039 + 1181 removed. Point colour = the dropped
#   donor's condition. Any refit that loses FDR < 0.05 is explained in the caption,
#   written from the tables.
#
# Inputs   OUT_ROOT/03-DE_analysis/donor_sensitivity/astro_composition_lodo_allfits_no1040vst.csv
#          OUT_ROOT/03-DE_analysis/donor_sensitivity/astro_composition_lodo_summary_no1040vst.csv
#          OUT_ROOT/06-astrocyte_analysis/astro_per_donor_counts_no1040vst.csv
# Outputs  OUT_ROOT/06-astrocyte_analysis/figures/ask2_lodo_forest_<subtype>.png/.pdf/_caption.txt
#
# Usage: revision/docker/run.sh revision/tools/astro_lodo_forest.R [--subtype Reactive-PWH]
# ==============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
source(file.path(SCORCH_REPO, "revision/tools/astro_pub_theme.R"))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
ST <- getarg("--subtype", "Reactive-PWH")
DS <- file.path(OUT_ROOT, "03-DE_analysis", "donor_sensitivity")
AS <- file.path(OUT_ROOT, "06-astrocyte_analysis")
OUT <- file.path(AS, "figures")

fits <- utils::read.csv(file.path(DS, "astro_composition_lodo_allfits_no1040vst.csv"), stringsAsFactors = FALSE)
summ <- utils::read.csv(file.path(DS, "astro_composition_lodo_summary_no1040vst.csv"), stringsAsFactors = FALSE)
cnt  <- utils::read.csv(file.path(AS, "astro_per_donor_counts_no1040vst.csv"), check.names = FALSE, stringsAsFactors = FALSE)
cond <- setNames(as.character(pwh_relabel(cnt$HIV)), as.character(cnt$donor))

f <- fits[fits$subtype == ST, ]
s <- summ[summ$subtype == ST, ]
if (!nrow(f) || !nrow(s)) stop("no leave-one-donor-out rows for ", ST)
f$label <- paste0("without ", f$dropped)
f$group <- cond[as.character(f$dropped)]
top <- data.frame(label = c("all donors", "without 1039 + 1181"),
                  odds_ratio = c(s$odds_ratio_full, s$or_excl_pathology),
                  FDR = c(s$fdr_full, s$fdr_excl_pathology), group = "reference")
D <- rbind(top, f[order(f$odds_ratio), c("label", "odds_ratio", "FDR", "group")])
D$label <- factor(D$label, levels = rev(D$label))
D$sig <- ifelse(D$FDR < 0.05, "FDR < 0.05", "FDR ≥ 0.05")

p <- ggplot(D, aes(odds_ratio, label)) +
  geom_vline(xintercept = 1, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = s$odds_ratio_full, colour = "grey40", linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(colour = group, shape = sig), size = 1.8, stroke = 0.6) +
  scale_colour_manual(values = c(reference = "black", PWH_COLOURS),
                      breaks = c("PWoH", "PWH"), name = "Dropped donor") +
  scale_shape_manual(values = c(`FDR < 0.05` = 16, `FDR ≥ 0.05` = 1), name = NULL) +
  scale_x_log10(name = sprintf("%s odds ratio, PWH vs PWoH (log scale)", ST)) + ylab(NULL) +
  theme_astro_pub() + theme(axis.text.y = element_text(size = 6))
fp <- file.path(OUT, sprintf("ask2_lodo_forest_%s.png", gsub("[^A-Za-z0-9]+", "_", ST)))
ggsave(fp, p, width = 3.6, height = 5.2, dpi = 600, bg = "white")
ggsave(sub("\\.png$", ".pdf", fp), p, width = 3.6, height = 5.2, device = grDevices::cairo_pdf, bg = "white")
cap <- sprintf(paste(
  "Leave-one-donor-out refits of the donor-level composition model (quasibinomial, one observation per donor,",
  "adjusted for age, sex, race and PMI). Each row drops one donor; point colour = that donor's condition; filled",
  "= FDR < 0.05 across the seven subtypes. Dashed line: all-donor odds ratio (%.1f, FDR %.3f). Odds ratios range",
  "%.1f-%.1f; %d of %d refits have FDR < 0.05."), s$odds_ratio_full, s$fdr_full, s$lodo_or_min, s$lodo_or_max,
  s$lodo_n_sig, s$lodo_n_fits)
## Explain any refit that loses FDR < 0.05 in the caption itself. Written from the
##   tables, so a changed result changes the sentence.
pct <- setNames(100 * cnt[[ST]] / cnt$total, as.character(cnt$donor))
for (i in which(f$FDR >= 0.05)) {
  d <- as.character(f$dropped[i]); cd <- cond[d]
  same <- names(cond)[cond == cd & names(cond) != d]
  bigger <- f$odds_ratio[i] > s$odds_ratio_full
  cap <- paste(cap, sprintf(paste(
    "The only refit above FDR 0.05 drops %s, a %s donor with %.1f%% %s nuclei (%d nuclei; the next-highest %s donor",
    "has %.1f%%). Without it the odds ratio %s to %.1f (p = %.3f; FDR %.3f after adjustment across the seven",
    "subtypes)%s"),
    d, cd, pct[d], ST, cnt[[ST]][match(d, cnt$donor)], cd, max(pct[same]),
    if (bigger) "rises" else "falls", f$odds_ratio[i], f$p[i], f$FDR[i],
    if (bigger && cd == "PWoH") paste0(": the effect grows, and the loss of significance reflects a control rate ",
                                       "estimated from fewer nuclei, not a weaker difference.") else "."))
}
writeLines(strwrap(cap, 110), sub("\\.png$", "_caption.txt", fp))
message("[lodo] ", fp)
message("[lodo] done")
