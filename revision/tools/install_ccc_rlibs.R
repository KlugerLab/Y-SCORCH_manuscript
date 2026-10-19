# ==============================================================================
# Installs GSVA, nichenetr and multinichenetr for the outlier-donor sensitivity
# reruns of the cell-cell communication (MultiNicheNet) and LDAM (ssGSEA)
# analyses, which follow Junchen's analysis.
#
# Packages go into a separate library, RLIBS_CCC, not the shared library that
# revision/config/paths.R prepends for every script: multinichenetr pulls in many
# dependencies that would otherwise shadow the image's versions for every
# analysis. Only tools that opt in (.libPaths(c(RLIBS_CCC, .libPaths()))) see
# these packages.
#
# upgrade = "never": anything the image already has is used as is; only missing
# packages are installed. Versions are written to RLIBS_CCC/INSTALLED_VERSIONS.txt
# because they may differ from the original analysis environment, which matters
# for MultiNicheNet's prioritisation scores.
#
# Usage: run in the image with RLIBS_CCC mounted read-write and network access:
#   docker run --rm --user $(id -u):$(id -g) -v <RLIBS_CCC>:<RLIBS_CCC> \
#     -v <repo>:<repo> scorch_seurat_v5:latest Rscript <repo>/revision/tools/install_ccc_rlibs.R
# ==============================================================================

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))
LIB <- RLIBS_CCC
stopifnot(dir.exists(LIB))
.libPaths(c(LIB, .libPaths()))
options(repos = BiocManager::repositories(), timeout = 600, Ncpus = 4)

BiocManager::install("GSVA", lib = LIB, update = FALSE, ask = FALSE)
remotes::install_github("saeyslab/nichenetr", lib = LIB, upgrade = "never",
                        dependencies = TRUE, build_vignettes = FALSE)
remotes::install_github("saeyslab/multinichenetr", lib = LIB, upgrade = "never",
                        dependencies = TRUE, build_vignettes = FALSE)

pk <- c("GSVA", "nichenetr", "multinichenetr")
ok <- vapply(pk, requireNamespace, logical(1), quietly = TRUE)
v <- vapply(pk, function(p) if (requireNamespace(p, quietly = TRUE))
  as.character(packageVersion(p)) else "NOT INSTALLED", character(1))
added <- setdiff(rownames(installed.packages(lib.loc = LIB)), character(0))
writeLines(c(sprintf("installed %s", format(Sys.time())), R.version.string,
             sprintf("Bioconductor %s", BiocManager::version()),
             sprintf("%s %s", pk, v), "",
             sprintf("all packages in %s (%d):", LIB, length(added)),
             sprintf("  %s %s", added, installed.packages(lib.loc = LIB)[added, "Version"])),
           file.path(LIB, "INSTALLED_VERSIONS.txt"))
print(v)
if (!all(ok)) stop("not installed: ", paste(pk[!ok], collapse = ", "))
message("[install] done")
