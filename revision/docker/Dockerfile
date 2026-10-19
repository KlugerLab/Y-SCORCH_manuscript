# ==============================================================================
# Y-SCORCH revision — single image covering stages 01-07.
#
# Base: the official Bioconductor image, which ships a coherent R + Bioconductor
# pair plus the system libraries Bioconductor packages need, and installs from
# binaries rather than source.
#
# Why not the original 2024 stack (rocker/rstudio:4.4.1 + Bioconductor 3.19): Bioc
# 3.19's ggtree 3.12.0 calls getFromNamespace("check_linewidth", "ggplot2"), an
# internal removed from current ggplot2, so ggtree -> enrichplot -> clusterProfiler
# fail to install; pinning ggplot2 3.5.1 fixes that and breaks fgsea instead. One
# coherent current release avoids both. It also brings edgeR 4, whose quasi-likelihood
# test is the one the revision uses (with a fallback for negative F statistics, see
# revision/config/paths.R).
#
# Paths below are load-bearing — the scripts hardcode them:
#   /bin/fast_tsne                 15 RunTSNE(tsne.method = "FIt-SNE") calls
#   /home/rstudio/.local/bin/macs2 3 Signac CallPeaks() calls
#
# Each install block verifies itself, so a failure names the block that broke
# rather than surfacing 40 packages later. install.packages() only warns on
# failure — an unverified image is not a working image.
# ==============================================================================
FROM bioconductor/bioconductor_docker:RELEASE_3_21

ARG DEBIAN_FRONTEND=noninteractive

# ---- System libraries --------------------------------------------------------
# cmake is needed by nloptr (-> car -> rstatix -> ggpubr); libglpk-dev supplies
# the headers igraph builds against (the runtime libglpk40 alone is not enough).
RUN apt-get update && apt-get install -y \
      git wget cmake \
      libhdf5-dev libglpk-dev libglpk40 libgsl-dev zlib1g-dev libcairo2-dev \
      libx11-dev libxt-dev liblzma-dev libbz2-dev libcurl4-openssl-dev \
      libfftw3-dev libxml2-dev libpng-dev libjpeg-dev \
      libharfbuzz-dev libfribidi-dev libfreetype6-dev libtiff5-dev \
      libmagick++-dev libgeos-dev libudunits2-dev libgdal-dev \
      python3 python3-pip python3-dev build-essential \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Some packages still declare CXX_STD = CXX11, which current Boost headers reject
# ('is_final' has not been declared in 'std'). fgsea is the one that bites here.
RUN printf 'CXX11STD = -std=gnu++17\nCXX14STD = -std=gnu++17\nCXX17STD = -std=gnu++17\n' \
      >> /usr/local/lib/R/etc/Makevars.site

# ---- Align BiocManager with the packages the base image actually ships -------
# The RELEASE_3_21 image is mislabelled: its BiocManager reports 3.21, but the
# installed packages are the 3.22 set (GenomeInfoDb 1.46.2, where Bioc 3.21 ships
# 1.44.3). BiocManager then refuses to install the older 3.21 versions over the
# newer ones ("version(s) same as or greater than current"), so anything it does
# install is built against a GenomeInfoDb from a different release. The result is
# S4 dispatch failure at load time, not install time:
#
#   BSgenome.Hsapiens.UCSC.hg38 installs, prints DONE, then fails to load with
#   "unable to find an inherited method for 'seqinfo' for 'TwobitNamedSequences'"
#
# Declaring 3.22 makes the repo match what is on disk. Verified: BSgenome,
# TFBSTools and EnsDb.Hsapiens.v86 all load cleanly afterwards. This runs before
# the install blocks so the entire stack is resolved against one release.
RUN R -e "BiocManager::install(version = '3.22', ask = FALSE); \
    stopifnot(as.character(BiocManager::version()) == '3.22'); \
    cat('BiocManager aligned to', as.character(BiocManager::version()), '\n')"

# ---- FIt-SNE -> /bin/fast_tsne ----------------------------------------------
# Every tSNE in the paper goes through this binary.
WORKDIR /opt
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git \
    && g++ -std=c++11 -O3 \
         FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp \
         -o /bin/fast_tsne -pthread -lfftw3 -lm -Wno-address-of-packed-member \
    && test -x /bin/fast_tsne

# ---- CRAN --------------------------------------------------------------------
RUN R -e "install.packages(c( \
      'remotes','devtools', \
      'ggplot2','cowplot','patchwork','ggpubr','ggrepel','RColorBrewer','viridis','scales','circlize', \
      'dplyr','tidyr','purrr','tibble','forcats','stringr','reshape2','magrittr','data.table','rlang', \
      'Matrix','optparse','enrichR','R.utils','hdf5r','harmony' \
    )); \
    p <- c('ggplot2','cowplot','patchwork','ggpubr','ggrepel','RColorBrewer','circlize', \
           'dplyr','tidyr','purrr','stringr','reshape2','data.table','Matrix','optparse', \
           'enrichR','harmony','remotes'); \
    m <- p[!sapply(p, requireNamespace, quietly=TRUE)]; \
    if (length(m)) stop('CRAN BLOCK FAILED: ', paste(m, collapse=', ')); \
    cat('CRAN block OK\n')"

# ---- Seurat v5 stack ---------------------------------------------------------
RUN R -e "install.packages(c('Seurat','SeuratObject')); \
    setRepositories(ind=1:3, addURLs=c('https://satijalab.r-universe.dev','https://bnprks.r-universe.dev/')); \
    install.packages(c('BPCells','presto','glmGamPoi','Signac')); \
    p <- c('Seurat','SeuratObject','Signac','presto','glmGamPoi'); \
    m <- p[!sapply(p, requireNamespace, quietly=TRUE)]; \
    if (length(m)) stop('SEURAT BLOCK FAILED: ', paste(m, collapse=', ')); \
    cat('Seurat block OK\n')"

# ---- Bioconductor ------------------------------------------------------------
# Version is set above, once, for the whole image.
# GOSim and Repitools are deliberately absent. GOSim was removed from
# Bioconductor; Repitools depends on gsmoothr, archived from CRAN. Every script
# that library()s either one calls no function from it — both are dead imports.
RUN R -e "BiocManager::install(c( \
      'edgeR','limma','DESeq2','SingleCellExperiment','scater','muscat','scDblFinder', \
      'biomaRt','org.Hs.eg.db','fgsea','ggtree','clusterProfiler','enrichplot', \
      'EnsDb.Hsapiens.v86','BSgenome.Hsapiens.UCSC.hg38','Homo.sapiens', \
      'JASPAR2020','TFBSTools','motifmatchr','chromVAR','ComplexHeatmap' \
    ), ask=FALSE, update=FALSE); \
    p <- c('edgeR','limma','DESeq2','SingleCellExperiment','scater','muscat', \
           'org.Hs.eg.db','fgsea','ggtree','clusterProfiler','enrichplot', \
           'EnsDb.Hsapiens.v86','BSgenome.Hsapiens.UCSC.hg38','Homo.sapiens', \
           'JASPAR2020','TFBSTools'); \
    m <- p[!sapply(p, requireNamespace, quietly=TRUE)]; \
    if (length(m)) stop('BIOC BLOCK FAILED: ', paste(m, collapse=', ')); \
    cat('Bioconductor block OK\n')"

# ---- GitHub-only packages ----------------------------------------------------
# None of these are on CRAN; install.packages() only warns, so they were silently
# absent from the earlier image.
#
# NICHES needs BisRNA, which has been archived from CRAN — so it must come from
# the archive first, or NICHES fails with "dependency 'BisRNA' is not available".
# Verified: all four install and load with this ordering.
RUN R -e "install.packages( \
      'https://cran.r-project.org/src/contrib/Archive/BisRNA/BisRNA_0.2.2.tar.gz', \
      repos = NULL, type = 'source'); \
    stopifnot(requireNamespace('BisRNA', quietly=TRUE)); \
    remotes::install_github('neurorestore/Libra', upgrade='never'); \
    remotes::install_github('GreenleafLab/chromVARmotifs', upgrade='never'); \
    remotes::install_github('msraredon/NICHES', ref='master', upgrade='never'); \
    p <- c('BisRNA','Libra','chromVARmotifs','NICHES'); \
    m <- p[!sapply(p, requireNamespace, quietly=TRUE)]; \
    if (length(m)) stop('GITHUB BLOCK FAILED: ', paste(m, collapse=', ')); \
    cat('GitHub block OK\n')"

# ---- MAGICAL (stage 06 gene-regulatory networks) -----------------------------
# run_magical.R sources MAGICAL_functions.R; config/paths.R reads MAGICAL_DIR.
# The repo is xichensf/magical, and MAGICAL_functions.R sits in its R/ subdir.
# Symlinked to the repo root so config/paths.R's MAGICAL_DIR resolves unchanged.
# Not tolerant of failure: a silent `|| echo WARNING` here previously produced an
# image with no MAGICAL at all, which only surfaced on inspection.
RUN git clone --depth 1 https://github.com/xichensf/magical.git /opt/MAGICAL \
    && ln -sf /opt/MAGICAL/R/MAGICAL_functions.R /opt/MAGICAL/MAGICAL_functions.R \
    && test -f /opt/MAGICAL/MAGICAL_functions.R \
    && echo 'MAGICAL OK'
ENV MAGICAL_DIR=/opt/MAGICAL

# ---- MACS2 for Signac CallPeaks ----------------------------------------------
# MACS2 2.2.9.1 (the last MACS2 release, and the version the published analysis
# used — see software/MACS2-2.2.9.1.dist-info in the analysis backup) cannot
# build on this image's Python 3.12: its Cython sources fail with
# "CompileError: MACS2/IO/CallPeakUnit.pyx". Pinning Cython<3 does not help.
#
# It gets its own isolated Python 3.10 environment rather than being replaced by
# MACS3. MACS3 would pip-install cleanly, but it calls peaks differently, and the
# premise of the revision is that ambient RNA removal is the only intended change
# from the published analysis — a different peak caller would add a second
# variable to every ATAC result.
#
# Signac shells out to this binary from inside R, so it has to live in this
# image; the peak-calling scripts hardcode /home/rstudio/.local/bin/macs2.
# curl comes from the base image; wget is the fallback (installed above) so this
# does not require touching the apt layer.
#
# Two things this block has to get right:
#
# 1. The bioconda BINARY is unusable on this image's glibc. Its Cython
#    extensions were linked against a glibc that still exported the finite-math
#    aliases, so loading MACS2/IO/ScoreTrack.*.so dies with
#        undefined symbol: __log_finite
#    glibc dropped those in 2.28+. `macs2 --version` does NOT import ScoreTrack,
#    so it happily reports 2.2.9.1 while peak calling is broken - which is how
#    this reached a real run. Rebuilding from source inside the same Python 3.10
#    env links against the current glibc and fixes it. (The Cython build failure
#    noted above is specific to Python 3.12; it compiles fine on 3.10.)
#
# 2. The /home/rstudio symlink is unreachable for non-root users. /home/rstudio
#    is 0750 rstudio:rstudio, and docker/run.sh passes --user <uid>, so
#    CallPeaks reports "MACS2 not found". config/paths.R::MACS2_PATH resolves to
#    the first EXECUTABLE candidate and prefers /opt/macs2env/bin/macs2; the
#    symlink is kept only for anything still hardcoding the old path, and the
#    directory is made traversable so it actually works.
#
# The verification below imports ScoreTrack and runs a real callpeak on a
# planted peak. Checking --version is not sufficient and never was.
RUN (curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
     || wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest) \
      | tar -xj -C /usr/local bin/micromamba \
    && MAMBA_ROOT_PREFIX=/opt/mamba /usr/local/bin/micromamba create -y \
         -p /opt/macs2env -c conda-forge -c bioconda python=3.10 macs2 \
    && /opt/macs2env/bin/pip install --no-cache-dir --no-binary :all: \
         --force-reinstall "MACS2==2.2.9.1" \
    && mkdir -p /home/rstudio/.local/bin \
    && ln -sf /opt/macs2env/bin/macs2 /home/rstudio/.local/bin/macs2 \
    && chown -R rstudio:rstudio /home/rstudio/.local \
    && chmod 0755 /home/rstudio /home/rstudio/.local /home/rstudio/.local/bin \
    && /opt/macs2env/bin/python -c "from MACS2.IO import ScoreTrack" \
    && printf 'chr1\t%s\t%s\tr\t0\t+\n' $(seq 10000 10400) > /tmp/_m.bed \
    && /opt/macs2env/bin/macs2 callpeak -t /tmp/_m.bed -f BED -n _m \
         --outdir /tmp/_mout -g 2.5e8 --nomodel >/dev/null 2>&1 \
    && test -s /tmp/_mout/_m_peaks.narrowPeak \
    && rm -rf /tmp/_m.bed /tmp/_mout \
    && echo 'MACS2 OK (imports + calls peaks)'

# ---- Python: the microglia feature-selection notebook ------------------------
# In a venv, not --user. Installing into the user site leaves the distro's own
# pandas at /usr/lib/python3/dist-packages ahead on sys.path; that copy is built
# against numpy 1.x, so importing it alongside a user-site numpy 2.x fails with
#   "ValueError: numpy.dtype size changed ... Expected 96 from C header, got 88"
# A venv without system site-packages keeps the distro copies out entirely.
RUN python3 -m venv /opt/pyenv \
    && /opt/pyenv/bin/pip install --no-cache-dir \
         numpy scipy pandas matplotlib scikit-learn scanpy anndata jupyter \
    && /opt/pyenv/bin/python -c "import scanpy, anndata, sklearn, pandas, matplotlib" \
    && ln -sf /opt/pyenv/bin/jupyter /usr/local/bin/jupyter \
    && ln -sf /opt/pyenv/bin/python /usr/local/bin/pyenv-python
ENV PYENV_PYTHON=/opt/pyenv/bin/python
#
# NOTE: 03-DE_analysis/09-micro_feature_selection-allRegions.ipynb is NOT covered
# by this venv. Its saved output shows it ran under a conda env named
# bipca-experiment on Python 3.10, and it imports `bipca`, a Kluger lab package
# that is not on PyPI. Run that notebook in the existing bipca-experiment image
# instead. It produces the 642-gene classifier and the 0.990/0.935/0.873
# balanced accuracies, so it does need to run — just not here.

# ---- Final verification ------------------------------------------------------
RUN R -e "pkgs <- c( \
      'Seurat','SeuratObject','Signac','harmony','Matrix', \
      'SingleCellExperiment','scater','muscat','DESeq2','edgeR','limma', \
      'clusterProfiler','enrichplot','enrichR','org.Hs.eg.db', \
      'EnsDb.Hsapiens.v86','BSgenome.Hsapiens.UCSC.hg38','Homo.sapiens', \
      'JASPAR2020','TFBSTools','chromVARmotifs','NICHES','Libra', \
      'ggplot2','cowplot','ggpubr','ggrepel','patchwork','circlize','RColorBrewer', \
      'dplyr','tidyr','purrr','stringr','data.table','optparse','reshape2'); \
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
    if (length(missing)) stop('MISSING PACKAGES: ', paste(missing, collapse=', ')); \
    cat('All', length(pkgs), 'packages present\n'); \
    for (p in c('edgeR','Seurat','ggplot2','clusterProfiler')) \
      cat(sprintf('  %-16s %s\n', p, as.character(packageVersion(p))))"

WORKDIR /work
