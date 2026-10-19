# REVISION: run on the CellBender-corrected data; paths come from
# revision/config/paths.R and outputs go to OUT_ROOT. Changes are marked ## REVISION.
#
## REVISION: one job per region prepares every requested cell type. The
##   published script took a single --celltype and reloaded a 9 GB RNA object
##   and a 9-17 GB ATAC object for each one. The per-cell-type steps are the
##   published ones.
##
##   --celltype  one name, a comma list, or "auto" = every cell type with
##               >= --min-degs DE genes and >= --min-peaks DA peaks at this
##               script's own cutoffs (|logFC| >= 0.5, p_adj.glb <= 0.05).
##               The 25 / 25 cutoff gives INS 15, PFC 9, VST 7 cell types.
##               Written to <MAGICAL_INPUT_DIR>/<REGION>_celltype_eligibility.csv.
##   Each finished cell-type folder gets a _DONE marker (run_magical.R checks it);
##   <REGION>_prepared.csv records what was written.
##
##   revision/docker/run.sh 06-astrocyte_analysis/gene-regulatory-networks/prepare_MAGICAL_input.R \
##       --region PFC --celltype auto

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
library(JASPAR2020)
## REVISION: provides the MAGICAL_* paths, OUT_ROOT and uncorrected_counts().
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))

option_list <- list(
  make_option(c("-r", "--region"), type = "character", default = "PFC",
              help = "Input region", metavar = "character"),
  make_option(c("-c", "--celltype"), type = "character", default = "Astro",
              help = "Input cell type, a comma-separated list, or 'auto'", metavar = "character"),
  ## REVISION: thresholds for --celltype auto. dest is explicit because optparse
  ##   keeps the hyphen ("min-degs"), which would leave opt$min_degs NULL.
  make_option("--min-degs", type = "integer", default = 25L, dest = "min_degs",
              help = "auto: minimum DE genes (|logFC| >= 0.5, p_adj.glb <= 0.05)"),
  make_option("--min-peaks", type = "integer", default = 25L, dest = "min_peaks",
              help = "auto: minimum DA peaks at the same cutoffs")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt$region)
print(opt$celltype)
dir.create(MAGICAL_INPUT_DIR, recursive = TRUE, showWarnings = FALSE)


#*****************DE results (read first: cheap, and they decide the cell types)*******************

## REVISION: reads the flat resDS tables (*_withFreq.RDS: one row per feature
##   per cell type, with the same logFC and p_adj.glb as the pbDS object's
##   $table plus frequency columns) and filters on cluster_id, instead of the
##   pbDS object the published script indexed with $table$`HIV-CTR`[[celltype]].
rna_de_path <- file.path(OUT_ROOT, "03-DE_analysis", opt$region, "results",
                         paste0(opt$region, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))
atac_de_path <- ## REVISION: the published code read results_ATAC/INS_*.RDS for every
##   region (the published files were all named INS_*). 04-ATAC_DE writes
##   <REGION>_HIV-CTR_res_edgeR_filtered_withFreq.RDS in data_ATAC/.
file.path(OUT_ROOT, "04-ATAC_DE", opt$region, "data_ATAC",
          paste0(opt$region, "_HIV-CTR_res_edgeR_filtered_withFreq.RDS"))
rna_de_all  <- as.data.frame(readRDS(rna_de_path))
atac_de_all <- as.data.frame(readRDS(atac_de_path))

sig_n <- function(de) tapply(abs(de$logFC) >= 0.5 & de$p_adj.glb <= 0.05,
                             as.character(de$cluster_id), sum, na.rm = TRUE)
n_deg <- sig_n(rna_de_all); n_peak <- sig_n(atac_de_all)
all_ct <- union(names(n_deg), names(n_peak))
elig <- data.frame(region = opt$region, cell_type = all_ct,
                   DEGs = as.integer(n_deg[all_ct]), DA_peaks = as.integer(n_peak[all_ct]),
                   stringsAsFactors = FALSE)
elig$DEGs[is.na(elig$DEGs)] <- 0L; elig$DA_peaks[is.na(elig$DA_peaks)] <- 0L
elig$eligible <- elig$DEGs >= opt$min_degs & elig$DA_peaks >= opt$min_peaks
elig <- elig[order(!elig$eligible, -elig$DEGs), ]
write.csv(elig, file.path(MAGICAL_INPUT_DIR, paste0(opt$region, "_celltype_eligibility.csv")), row.names = FALSE)
print(elig, row.names = FALSE)

celltypes <- if (identical(opt$celltype, "auto")) elig$cell_type[elig$eligible] else
  trimws(strsplit(opt$celltype, ",")[[1]])
not_elig <- setdiff(celltypes, elig$cell_type[elig$eligible])
if (length(not_elig))
  warning("below the ", opt$min_degs, "/", opt$min_peaks, " cutoff but requested: ",
          paste(not_elig, collapse = ", "), immediate. = TRUE)
## astrocytes first, so the figure cell type's inputs are ready soonest
celltypes <- c(intersect("Astro", celltypes), setdiff(celltypes, "Astro"))
message("[magical-prep] ", opt$region, ": preparing ", length(celltypes), " cell types: ",
        paste(celltypes, collapse = ", "))


#*****************load data, subset to region and celltype, and match cells*******************

## REVISION: seurat_combined_v2 holds donor 1040's VST cells and a single
##   "MSNs" class. seurat_combined_v3_no1040_msn has 1040 removed and the nine
##   MSN subtypes in consensus_annotation_v3.
scRNA_all <- readRDS(file.path(OUT_ROOT, "01-Regional_data_processing/combined/data/seurat_combined_v3_no1040_msn.RDS"))
# subset to region + the cell types being prepared
scRNA_all <- scRNA_all[,which(scRNA_all$region == opt$region)]
scRNA_all <- scRNA_all[,which(as.character(scRNA_all$consensus_annotation_v3) %in% celltypes)]

## REVISION: the RNA side uses uncorrected counts, like the DE that supplies
##   MAGICAL's candidate genes. CellBender zeroes genuine low PWoH expression of
##   reactive-astrocyte genes (SERPINA3, OSMR), which would fake condition
##   differences. uncorrected_counts() (revision/config/paths.R) pulls the raw
##   CellRanger assay for exactly these cells, in the object's cell order.
raw_counts <- uncorrected_counts(scRNA_all, region = opt$region)
rna_meta <- scRNA_all@meta.data
rm(scRNA_all); invisible(gc())
scRNA_all <- CreateSeuratObject(counts = raw_counts, meta.data = rna_meta)
rm(raw_counts); invisible(gc())

# rename RNA cell names to make it the same format as atac
rna_cell_names <- colnames(scRNA_all)
new_rna_cell_names <- sapply(rna_cell_names, function(name) {
  parts <- unlist(strsplit(name, split = "-"))
  new_name <- paste(parts[3], paste(parts[1], parts[2], sep = "-"), sep = "_")
  return(new_name)
})
stopifnot(!anyDuplicated(new_rna_cell_names))
scRNA_all <- RenameCells(scRNA_all, new.names = unname(new_rna_cell_names))

## REVISION: MAGICAL must read the object and assay the ATAC DE was run on,
##   because it inner-joins the candidate DA peaks to this peak list on exact
##   coordinates. That is signac_integrated_v2 (ATAC assay) under sampleQC_redo/
##   for INS and PFC (the published INS read the v1 object) and
##   signac_integrated_v2_msn for VST.
atac_file <- if (opt$region == "VST") "signac_integrated_v2_msn.RDS" else "signac_integrated_v2.RDS"
atac_seurat_path <- file.path(OUT_ROOT, "02-ATAC_processing/sampleQC_redo", opt$region, "data", atac_file)
scATAC_all <- readRDS(atac_seurat_path)
DefaultAssay(scATAC_all) <- "ATAC"

matched <- intersect(colnames(scRNA_all), Cells(scATAC_all))
message("[magical-prep] ", length(matched), " of ", ncol(scRNA_all), " RNA cells found in ", atac_file)
if (length(matched) < 0.1 * ncol(scRNA_all))
  stop("barcodes do not line up - RNA example ", colnames(scRNA_all)[1],
       ", ATAC example ", Cells(scATAC_all)[1])
scATAC_all <- scATAC_all[, matched]
scRNA_all  <- scRNA_all[, matched]

## REVISION: MAGICAL pairs the two modalities only through subject_ID
##   (intersect of the RNA and ATAC meta columns). The published script wrote
##   sample_name from each object; take it from the shared barcode prefix
##   (<sample_ids>_<barcode>) so the two files cannot disagree.
scRNA_all$subject  <- sub("_[^_]*$", "", colnames(scRNA_all))
scATAC_all$subject <- sub("_[^_]*$", "", Cells(scATAC_all))

main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)

#TF motif set (identical for every cell type, so built once)
opts <- list(
  species = 9606,         # NCBI Taxonomy ID for Homo sapiens
  collection = "CORE",     # CORE collection of JASPAR
  all_versions = FALSE
)

# Retrieve the matrix set for Homo sapiens
pfm <- getMatrixSet(JASPAR2020, opts)
motif_names <- sapply(pfm@listData, function(x) x@name)
names(pfm) <- motif_names

extract_tf_name <- function(name) {
  # Case 1: Remove '(var.x)'
  name <- sub("\\(var\\.\\d+\\)", "", name)

  # Case 2: Keep only the first part before '::'
  name <- sub("::.*$", "", name)

  # Case 3: Return the processed name (unchanged for other cases)
  return(name)
}

# Apply the function to all motif names
tf_names <- sapply(motif_names, extract_tf_name)
names(tf_names) <- NULL

print("finish load and matched assay")

prepared <- list()
for (Cell_type_selected in celltypes) {
  t0 <- Sys.time()
  message("\n########## ", opt$region, " ", Cell_type_selected, " ##########")
  ## REVISION: magical_tag() replaces paste0(region, "_new_", celltype), which
  ##   turns "L2/3 IT" and "D1-Shell/OT" into nested folders; the folder is
  ##   also created here.
  folder_name <- magical_tag(opt$region, Cell_type_selected)
  save_folder <- paste0(paste0(MAGICAL_INPUT_DIR, "/"),folder_name,"/")
  dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(save_folder, "_DONE"))

  res <- tryCatch({
    scRNA <- scRNA_all[,which(scRNA_all$consensus_annotation_v3 == Cell_type_selected)]
    ## REVISION: same cells in the same order in both modalities
    scATAC <- scATAC_all[, colnames(scRNA)]
    stopifnot(identical(colnames(scRNA), Cells(scATAC)))

    #*****************Deal with RNA assay*******************
    sample_exclude <- names(which(table(scRNA$subject) <= 1))
    scRNA <- scRNA[,which(!scRNA$subject %in% sample_exclude)]
    scATAC <- scATAC[,which(!scATAC$subject %in% sample_exclude)]

    scRNA <- FindTopFeatures(scRNA, min.cutoff = 5)
    scRNA <- scRNA[VariableFeatures(scRNA),]

    write.table(rownames(scRNA), file = paste0(save_folder,'scRNA genes.txt'), quote = FALSE, row.names = TRUE, col.names = FALSE, sep = '\t')

    #RNA assay cell read counts
    cell_index=which(scRNA$consensus_annotation_v3==Cell_type_selected)
    Cell_type_scRNA_counts = as(LayerData(scRNA, assay = "RNA", layer = "counts")[, cell_index], "dgTMatrix")
    write.table(summary(Cell_type_scRNA_counts), file=paste0(save_folder,'Cell type scRNA read count.txt'),quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')

    #RNA assay cell meta including barcode, cell type annotation, sample names/IDs and group/condition
    Cell_type_scRNA_meta = scRNA@meta.data[cell_index,]
    write.table(data.frame(rownames(Cell_type_scRNA_meta), Cell_type_scRNA_meta$consensus_annotation_v3, Cell_type_scRNA_meta$subject, Cell_type_scRNA_meta$condition_new),
                file = paste0(save_folder,'Cell type scRNA cell meta.txt'),
                quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = '\t')

    print("finish RNA assay")

    #*****************Deal with ATAC assay*******************
    keep.peaks <- as.logical(seqnames(granges(scATAC)) %in% main.chroms)
    scATAC <- scATAC[keep.peaks, ]

    scATAC <- FindTopFeatures(scATAC, min.cutoff = 5)
    scATAC <- scATAC[VariableFeatures(scATAC),]

    #ATAC assay peaks
    write.table(scATAC@assays$ATAC@ranges, file = paste0(save_folder,'scATAC peaks.txt'), quote = FALSE, col.names = FALSE,  sep = '\t')

    #ATAC assay cell read counts
    Cell_type_scATAC_counts = as(scATAC@assays$ATAC@counts, "dgTMatrix")
    write.table(summary(Cell_type_scATAC_counts), file=paste0(save_folder,'Cell type scATAC read count.txt'),
                quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')

    #ATAC assay cell meta including barcode, cell type annotation, and sample names/IDs
    Cell_type_scATAC_meta = scATAC@meta.data
    write.table(data.frame(rownames(Cell_type_scATAC_meta), Cell_type_scRNA_meta$consensus_annotation_v3, Cell_type_scATAC_meta$subject),
                file = paste0(save_folder,'Cell type scATAC cell meta.txt'),
                quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = '\t')

    #TF motif mapping
    scATAC <- AddMotifs(object = scATAC, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
    dense_matrix <- as.matrix(scATAC@assays$ATAC@motifs@data)
    Peak_motif_mapping <- as(dense_matrix, "dgTMatrix")
    write.table(summary(Peak_motif_mapping), file=paste0(save_folder,'Motif mapping prior.txt'),
                quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')
    rm(dense_matrix, Peak_motif_mapping)

    write.table(tf_names, file = paste0(save_folder,'Motifs.txt'), quote = FALSE, col.names = FALSE,  sep = '\t')

    print("finish ATAC assay")

    #*****************prepare gene list*******************
    all_gene <- rownames(scRNA)

    write.table(all_gene,
                file = paste0(save_folder,'Cell type candidate genes_allgene.txt'),
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = '\t')
    rna_de_micro <- rna_de_all[which(as.character(rna_de_all$cluster_id) == Cell_type_selected),]
    rna_de_micro_up <- rna_de_micro[which(rna_de_micro$logFC >= 0.5 & rna_de_micro$p_adj.glb <= 0.05),]
    rna_de_micro_down <- rna_de_micro[which(rna_de_micro$logFC <= -0.5 & rna_de_micro$p_adj.glb <= 0.05),]
    rna_de_micro_all <- rna_de_micro[which(abs(rna_de_micro$logFC) >= 0.5 & rna_de_micro$p_adj.glb <= 0.05),]

    write.table(rna_de_micro_all$gene,
                file = paste0(save_folder,'Cell type candidate genes_DEGs.txt'),
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = '\t')

    write.table(rna_de_micro_up$gene,
                file = paste0(save_folder,'Cell type candidate genes_DEGs_up.txt'),
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = '\t')

    write.table(rna_de_micro_down$gene,
                file = paste0(save_folder,'Cell type candidate genes_DEGs_down.txt'),
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = '\t')
    print("finish RNA list")

    #**************** Prepare atac list*****************
    atac_de_micro <- atac_de_all[which(as.character(atac_de_all$cluster_id) == Cell_type_selected),]
    atac_de_micro_all <- atac_de_micro[which(abs(atac_de_micro$logFC) >= 0.5 & atac_de_micro$p_adj.glb <= 0.05),]
    atac_de_micro_up <- atac_de_micro[which(atac_de_micro$logFC >= 0.5 & atac_de_micro$p_adj.glb <= 0.05),]
    atac_de_micro_down <- atac_de_micro[which(atac_de_micro$logFC <= -0.5 & atac_de_micro$p_adj.glb <= 0.05),]

    split_peaks <- function(ids) {
      split_data <- strsplit(as.character(ids), "-")
      rbindlist(lapply(split_data, function(x) {
        as.data.table(list(Chromosome = x[1], Start = x[2], End = x[3]))
      }))
    }

    write.table(split_peaks(atac_de_micro_all$gene),
                file = paste0(save_folder,'all_de_peak_Cell type candidate peaks.txt'),
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = '\t')

    write.table(split_peaks(atac_de_micro_up$gene),
                file = paste0(save_folder,'up_de_peak_Cell type candidate peaks.txt'),
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = '\t')

    write.table(split_peaks(atac_de_micro_down$gene),
                file = paste0(save_folder,'down_de_peak_Cell type candidate peaks.txt'),
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE,
                sep = '\t')
    print("Finish atac list")

    ## how many candidate DA peaks survive into the peak list MAGICAL loads
    peak_ids <- paste(seqnames(granges(scATAC)), start(granges(scATAC)), end(granges(scATAC)), sep = "-")
    writeLines(format(Sys.time()), file.path(save_folder, "_DONE"))
    data.frame(region = opt$region, cell_type = Cell_type_selected, folder = folder_name,
               cells = ncol(scRNA), samples = length(unique(scRNA$subject)),
               rna_genes = nrow(scRNA), atac_peaks = nrow(scATAC),
               candidate_genes = nrow(rna_de_micro_all),
               candidate_genes_in_rna = sum(rna_de_micro_all$gene %in% rownames(scRNA)),
               candidate_peaks = nrow(atac_de_micro_all),
               candidate_peaks_in_atac = sum(as.character(atac_de_micro_all$gene) %in% peak_ids),
               minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
               error = NA_character_)
  }, error = function(e) {
    message("[magical-prep] FAILED ", opt$region, " ", Cell_type_selected, ": ", conditionMessage(e))
    data.frame(region = opt$region, cell_type = Cell_type_selected, folder = folder_name,
               error = conditionMessage(e))
  })
  prepared[[Cell_type_selected]] <- res
  print(res)
  write.csv(bind_rows(prepared), file.path(MAGICAL_INPUT_DIR, paste0(opt$region, "_prepared.csv")), row.names = FALSE)
  invisible(gc())
}
message("[magical-prep] done ", opt$region)
