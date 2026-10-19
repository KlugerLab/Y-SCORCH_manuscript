# REVISION: run on the CellBender-corrected data; paths come from
# revision/config/paths.R and outputs go to OUT_ROOT. Changes are marked ## REVISION.


## REVISION: provides MAGICAL_DIR, the other MAGICAL_* paths and OUT_ROOT; it
##   must be sourced before MAGICAL_functions.R, whose location depends on MAGICAL_DIR.
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

library(Seurat)
library(Matrix)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(Signac)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARmotifs)
library(optparse)

## REVISION: MAGICAL itself lives in the image at /opt/MAGICAL (MAGICAL_DIR).
## Fail with a useful message rather than a bare "cannot open connection".
.magical_fn <- file.path(MAGICAL_DIR, "MAGICAL_functions.R")
if (!file.exists(.magical_fn))
  stop("MAGICAL_functions.R not found at ", .magical_fn,
       "\n  Set MAGICAL_DIR, or check the image.")
source(.magical_fn)

## REVISION: MAGICAL_functions.R needs hg38_Refseq.txt (a tab table of chr,
##   strand, start, end, gene symbol) to link peaks to genes, and it is not at
##   the expected path in a plain MAGICAL clone (see MAGICAL_REFSEQ in
##   revision/config/paths.R). Check up front rather than failing deep in the run.
if (!file.exists(MAGICAL_REFSEQ))
  stop("missing MAGICAL reference: ", MAGICAL_REFSEQ,
       "\n  This ships in the MAGICAL repo under input_data/. Copy it there,",
       "\n  or point MAGICAL_REFSEQ at a copy.")

option_list <- list(
  make_option(c("-r", "--region"), type = "character", default = "PFC",
              help = "Input region", metavar = "character"),
  make_option(c("-c", "--celltype"), type = "character", default = "Micro",
              help = "Input cell type", metavar = "character"),
  ## REVISION: a seed replicate on identical inputs, to measure how much of the
  ##   difference from the published circuits is MAGICAL's own sampling
  ##   variability. 123 = the primary run.
  make_option("--seed", type = "integer", default = 123L, dest = "seed",
              help = "RNG seed; anything but 123 writes to output/seed_replicates/"),
  ## REVISION: an alternative candidate set, prepared by
  ##   revision/tools/magical_frq_variant.R into input_data/<tag>_<variant>/. Circuits go to
  ##   output/variants/<variant>/ so the primary circuits are never overwritten.
  make_option("--variant", type = "character", default = "", dest = "variant",
              help = "read input_data/<tag>_<variant>/, write output/variants/<variant>/")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt$region)
print(opt$celltype)

## REVISION: timing, and a fixed seed so a rerun reproduces the circuits
##   (MAGICAL_estimation samples; the published script set no seed).
t0 <- Sys.time()
set.seed(opt$seed)
message("[magical] seed ", opt$seed)

## REVISION: "L2/3 IT" and "D1-Shell/OT" turn paste0(region, "_new_",
##   celltype) into nested folders. magical_tag() (revision/config/paths.R) is shared
##   with prepare_MAGICAL_input.R, which writes _DONE when a folder is complete.
folder_name <- magical_tag(opt$region, opt$celltype)
input_dir <- paste0(paste0(MAGICAL_INPUT_DIR, "/"),folder_name,
                    if (nzchar(opt$variant)) paste0("_", opt$variant) else "", "/")
if (!file.exists(file.path(input_dir, "_DONE")))
  stop("MAGICAL input for ", opt$region, " ", opt$celltype,
       " is missing or incomplete (no _DONE in ", input_dir, ")")


Candidate_gene_file_path = paste0(input_dir,'Cell type candidate genes_DEGs.txt')
Candidate_peak_file_path = paste0(input_dir,"all_de_peak_Cell type candidate peaks.txt")

# filtered scRNA data of the cell type
scRNA_readcount_file_path = paste0(input_dir,'Cell type scRNA read count.txt')
scRNA_gene_file_path = paste0(input_dir,'scRNA genes.txt')
scRNA_cellmeta_file_path = paste0(input_dir,'Cell type scRNA cell meta.txt')


# filtered scATAC data of the cell type
scATAC_readcount_file_path = paste0(input_dir,'Cell type scATAC read count.txt')
scATAC_peak_file_path = paste0(input_dir,'scATAC peaks.txt')
scATAC_cellmeta_file_path = paste0(input_dir,'Cell type scATAC cell meta.txt')

# TF motif prior on all ATAC peaks
Motif_mapping_file_path = paste0(input_dir,'Motif mapping prior.txt')
Motif_name_file_path = paste0(input_dir,'Motifs.txt')

Ref_seq_file_path = MAGICAL_REFSEQ

loaded_data <- Data_loading(Candidate_gene_file_path, Candidate_peak_file_path,
                            scRNA_readcount_file_path, scRNA_gene_file_path, scRNA_cellmeta_file_path,
                            scATAC_readcount_file_path, scATAC_peak_file_path, scATAC_cellmeta_file_path,
                            Motif_mapping_file_path, Motif_name_file_path, Ref_seq_file_path)


distance_control=5e5
Candidate_circuits <- Candidate_circuits_construction_without_TAD(loaded_data, distance_control)


#Model initialization
Initial_model<-MAGICAL_initialization(loaded_data, Candidate_circuits)

Circuits_linkage_posterior<-MAGICAL_estimation(loaded_data, Candidate_circuits, Initial_model, iteration_num = 1000)

## REVISION: create the output folder, and replace "/" and spaces in the cell type name
dir.create(MAGICAL_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
Output_path <- paste0(paste0(MAGICAL_OUTPUT_DIR, "/"),opt$region,"_",gsub("[/ ]", "-", opt$celltype),"MAGICAL_selected_regulatory_circuits.txt")
## REVISION: a seed replicate never overwrites the primary circuits. It goes to
##   output/seed_replicates/ with the seed in the name, and keeps the posterior
##   probabilities, so circuits near MAGICAL's 0.8 / 0.95 cutoffs can be examined.
if (opt$seed != 123L) {
  rep_dir <- file.path(MAGICAL_OUTPUT_DIR, "seed_replicates")
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  Output_path <- file.path(rep_dir, sub("\\.txt$", paste0("_seed", opt$seed, ".txt"), basename(Output_path)))
  saveRDS(list(posterior = Circuits_linkage_posterior, circuits = Candidate_circuits),
          sub("\\.txt$", "_posterior.rds", Output_path))
}

if (nzchar(opt$variant)) {
  var_dir <- file.path(MAGICAL_OUTPUT_DIR, "variants", opt$variant)
  dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
  Output_path <- file.path(var_dir, basename(Output_path))
}
MAGICAL_circuits_output(Output_file_path = Output_path,
                        Candidate_circuits, Circuits_linkage_posterior)
message("[magical] ", opt$region, " ", opt$celltype, " done in ",
        format(round(difftime(Sys.time(), t0, units = "mins"), 1)), " -> ", Output_path)