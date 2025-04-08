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

option_list <- list(
  make_option(c("-r", "--region"), type = "character", default = "PFC",
              help = "Input region", metavar = "character"),
  make_option(c("-c", "--celltype"), type = "character", default = "Astro",
              help = "Input cell type", metavar = "character")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt$region)
print(opt$celltype)


#*****************load data, subset to region and celltype, and match cells*******************

scRNA <- readRDS("/banach2/SCORCH/scripts/combined_analysis_paper/01-Regional_data_processing/combined/data/seurat_combined_v2.RDS")
# subset to region + celltype
scRNA <- scRNA[,which(scRNA$region == opt$region)]
scRNA <- scRNA[,which(scRNA$consensus_annotation_v3 == opt$celltype)]
if (opt$region == "INS"){
  atac_seurat_path <- paste0("/banach2/SCORCH/scripts/combined_analysis_paper/02-ATAC_processing/sampleQC/",opt$region,"/data/signac_integrated.RDS")
} else{
  atac_seurat_path <- paste0("/banach2/SCORCH/scripts/combined_analysis_paper/02-ATAC_processing/sampleQC/",opt$region,"/data/signac_integrated_v2.RDS")
}
scATAC <- readRDS(atac_seurat_path)

# rename RNA cell names to make it the same format as atac
rna_cell_names <- colnames(scRNA)
new_rna_cell_names <- sapply(rna_cell_names, function(name) {
  parts <- unlist(strsplit(name, split = "-"))
  new_name <- paste(parts[3], paste(parts[1], parts[2], sep = "-"), sep = "_")
  return(new_name)
})
colnames(scRNA) <- new_rna_cell_names

scATAC <- scATAC[,which(Cells(scATAC) %in% new_rna_cell_names)]
scRNA <- scRNA[,which(Cells(scRNA) %in% Cells(scATAC))]
table(Cells(scRNA) == Cells(scATAC))

folder_name <- paste0(opt$region,"_new_",opt$celltype)
save_folder <- paste0("/banach2/chang/new_SCORCH_analysis/MAGICAL/input_data/",folder_name,"/")

print("finish load and matched assay")

#*****************Deal with RNA assay*******************
Cell_type_selected = opt$celltype
sample_exclude <- names(which(table(scRNA$sample_name) <= 1))
scRNA <- scRNA[,which(!scRNA$sample_name %in% sample_exclude)]
scATAC <- scATAC[,which(!scATAC$sample_name %in% sample_exclude)]

scRNA <- FindTopFeatures(scRNA, min.cutoff = 5)
scRNA <- scRNA[VariableFeatures(scRNA),]

write.table(scRNA@assays$RNA$counts@Dimnames[[1]], file = paste0(save_folder,'scRNA genes.txt'), quote = FALSE, row.names = TRUE, col.names = FALSE, sep = '\t')


#RNA assay cell read counts
cell_index=which(scRNA$consensus_annotation_v3==Cell_type_selected)
Cell_type_scRNA_counts = as(scRNA@assays$RNA$counts[, cell_index], "dgTMatrix")
write.table(summary(Cell_type_scRNA_counts), file=paste0(save_folder,'Cell type scRNA read count.txt'),quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')


#RNA assay cell meta including barcode, cell type annotation, sample names/IDs and group/condition
Cell_type_scRNA_meta = scRNA@meta.data[cell_index,]
write.table(data.frame(rownames(Cell_type_scRNA_meta), Cell_type_scRNA_meta$consensus_annotation_v3, Cell_type_scRNA_meta$sample_name, Cell_type_scRNA_meta$condition_new),
            file = paste0(save_folder,'Cell type scRNA cell meta.txt'),
            quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = '\t')

print("finish RNA assay")

#*****************Deal with ATAC assay*******************

Cell_type_selected <- opt$celltype

DefaultAssay(scATAC) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
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
write.table(data.frame(rownames(Cell_type_scATAC_meta), Cell_type_scRNA_meta$consensus_annotation_v3, Cell_type_scATAC_meta$sample_name),
            file = paste0(save_folder,'Cell type scATAC cell meta.txt'),
            quote = FALSE, row.names = TRUE, col.names = FALSE,  sep = '\t')


#TF motif mapping
opts <- list(
  species = 9606,         # NCBI Taxonomy ID for Homo sapiens
  collection = "CORE",     # CORE collection of JASPAR
  all_versions = FALSE
)

# Retrieve the matrix set for Homo sapiens
pfm <- getMatrixSet(JASPAR2020, opts)
motif_names <- sapply(pfm@listData, function(x) x@name)
names(pfm) <- motif_names

scATAC <- AddMotifs(object = scATAC, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
dense_matrix <- as.matrix(scATAC@assays$ATAC@motifs@data)
Peak_motif_mapping <- as(dense_matrix, "dgTMatrix")
write.table(summary(Peak_motif_mapping), file=paste0(save_folder,'Motif mapping prior.txt'),
            quote = FALSE, row.names = FALSE, col.names = FALSE,  sep = '\t')

#motif_name <-Peak_motif_mapping@Dimnames[[2]]
#extract_tf_names <- function(motif_names) {
#  sapply(motif_names, function(x) {
#    parts <- strsplit(x, "_")[[1]]
#    tf_name <- parts[3]  # Assuming TF name is always in the second position
#    return(tf_name)
#  })
#}

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
#TF_name <- extract_tf_names(motif_name)

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
rna_de_path <- paste0("/banach2/SCORCH/scripts/combined_analysis_paper/03-DE_analysis/",opt$region,"/results/",opt$region,"_HIV-CTR_res_edgeR_filtered.RDS")
rna_de <-  readRDS(rna_de_path)
celltype <- opt$celltype
rna_de_micro <- rna_de$table$`HIV-CTR`[[celltype]]
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
atac_de_path <- paste0('/banach2/SCORCH/scripts/combined_analysis_paper/04-ATAC_DE/',opt$region,"/results_ATAC/INS_HIV-CTR_res_edgeR_filtered.RDS")
atac_de <- readRDS(atac_de_path)
atac_de_micro <- atac_de$table$`HIV-CTR`[[celltype]]
atac_de_micro_all <- atac_de_micro[which(abs(atac_de_micro$logFC) >= 0.5 & atac_de_micro$p_adj.glb <= 0.05),]
atac_de_micro_up <- atac_de_micro[which(atac_de_micro$logFC >= 0.5 & atac_de_micro$p_adj.glb <= 0.05),]
atac_de_micro_down <- atac_de_micro[which(atac_de_micro$logFC <= -0.5 & atac_de_micro$p_adj.glb <= 0.05),]



split_data <- strsplit(as.character(atac_de_micro_all$gene), "-")
split_data <- rbindlist(lapply(split_data, function(x) {
  as.data.table(list(Chromosome = x[1], Start = x[2], End = x[3]))
}))

write.table(split_data, 
            file = paste0(save_folder,'all_de_peak_Cell type candidate peaks.txt'),
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = '\t')

split_data <- strsplit(as.character(atac_de_micro_up$gene), "-")
split_data <- rbindlist(lapply(split_data, function(x) {
  as.data.table(list(Chromosome = x[1], Start = x[2], End = x[3]))
}))

write.table(split_data, 
            file = paste0(save_folder,'up_de_peak_Cell type candidate peaks.txt'),
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = '\t')

split_data <- strsplit(as.character(atac_de_micro_down$gene), "-")
split_data <- rbindlist(lapply(split_data, function(x) {
  as.data.table(list(Chromosome = x[1], Start = x[2], End = x[3]))
}))

write.table(split_data, 
            file = paste0(save_folder,'down_de_peak_Cell type candidate peaks.txt'),
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = '\t')
print("Finish atac list")
