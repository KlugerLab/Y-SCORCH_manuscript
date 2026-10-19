# ==============================================================================
# Astrocyte nuclei per region and condition, under each population definition in use.
#
# Method: two definitions count different populations:
#   all       every nucleus in the astrocyte object, including the Oligo-like and
#             Endo-Astro contamination clusters (the 74,890 total)
#   subtypes  the seven astrocyte subtypes only - the "clean" subset that Fig1D,
#             the composition tests and the reviewer figures use
# Both are reported for the current object (donor 1040's VST library removed) and,
# for reference, the object before that removal. A manuscript-ready sentence with
# the per-region counts is printed for each definition.
#
# Inputs   OUT_ROOT/06-astrocyte_analysis/data/astro_all_clusters_labeled[_no1040vst].RDS
# Outputs  OUT_ROOT/06-astrocyte_analysis/astro_nuclei_counts.csv
#
# Usage: revision/docker/run.sh revision/tools/astro_nuclei_counts.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
suppressPackageStartupMessages({ library(dplyr) })
D <- file.path(OUT_ROOT, "06-astrocyte_analysis")
CLEAN <- c("Homeostatic", "Reactive", "Reactive-PWH", "Signaling",
           "Synaptic", "Micro-Astro", "Mural-Astro")
OBJS <- c(current = "astro_all_clusters_labeled_no1040vst.RDS",
          before_1040_removal = "astro_all_clusters_labeled.RDS")

res <- list()
for (k in names(OBJS)) {
  f <- file.path(D, "data", OBJS[[k]])
  if (!file.exists(f)) { message("[count] missing ", f); next }
  md <- readRDS(f)@meta.data
  reg_col <- intersect(c("region", "Region", "brain_region"), colnames(md))[1]
  region <- if (!is.na(reg_col)) as.character(md[[reg_col]]) else
    regmatches(md$orig.ident, regexpr("INS|PFC|VST", md$orig.ident))
  cond <- as.character(md$condition_new)
  cond <- ifelse(toupper(cond) %in% c("HIV", "TRUE", "PWH"), "PWH", "PWoH")
  sub  <- as.character(md$subcluster_labels)
  message(sprintf("[count] %s: %d nuclei; subtypes present: %s", k, nrow(md),
                  paste(sort(unique(sub)), collapse = ", ")))
  for (def in c("all", "subtypes")) {
    keep <- if (def == "all") rep(TRUE, nrow(md)) else sub %in% CLEAN
    t <- as.data.frame(table(region = region[keep], condition = cond[keep]))
    t$object <- k; t$definition <- def
    res[[paste(k, def)]] <- t
  }
}
tab <- bind_rows(res) %>% tidyr::pivot_wider(names_from = condition, values_from = Freq) %>%
  mutate(total = PWH + PWoH) %>% select(object, definition, region, total, PWoH, PWH)
utils::write.csv(tab, file.path(D, "astro_nuclei_counts.csv"), row.names = FALSE)
options(width = 150); print(as.data.frame(tab), row.names = FALSE)

fmt <- function(x) format(x, big.mark = ",")
for (def in c("subtypes", "all")) {
  t <- tab %>% filter(object == "current", definition == def)
  s <- paste(vapply(c("INS", "PFC", "VST"), function(r) {
    x <- t[t$region == r, ]
    sprintf("%s nuclei (%s PWoH; %s PWH) from the %s", fmt(x$total), fmt(x$PWoH), fmt(x$PWH), r)
  }, character(1)), collapse = ", ")
  s <- sub(", ([^,]*)$", ", and \\1", s)
  message("\n[count] current object, ", def, ":\n  we analyzed ", s, ".")
}
