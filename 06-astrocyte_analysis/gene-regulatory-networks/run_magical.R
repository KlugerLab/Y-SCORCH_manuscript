
source('/banach2/chang/new_SCORCH_analysis/MAGICAL/MAGICAL_functions.R')

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

option_list <- list(
  make_option(c("-r", "--region"), type = "character", default = "PFC",
              help = "Input region", metavar = "character"),
  make_option(c("-c", "--celltype"), type = "character", default = "Micro",
              help = "Input cell type", metavar = "character")
)

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
print(opt$region)
print(opt$celltype)

folder_name <- paste0(opt$region,"_new_",opt$celltype)
input_dir <- paste0("/banach2/chang/new_SCORCH_analysis/MAGICAL/input_data/",folder_name,"/")


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

Ref_seq_file_path = '/banach2/chang/new_SCORCH_analysis/MAGICAL/input_data/hg38_Refseq.txt'

loaded_data <- Data_loading(Candidate_gene_file_path, Candidate_peak_file_path,
                            scRNA_readcount_file_path, scRNA_gene_file_path, scRNA_cellmeta_file_path,
                            scATAC_readcount_file_path, scATAC_peak_file_path, scATAC_cellmeta_file_path,
                            Motif_mapping_file_path, Motif_name_file_path, Ref_seq_file_path)


distance_control=5e5
Candidate_circuits <- Candidate_circuits_construction_without_TAD(loaded_data, distance_control)


#Model initialization
Initial_model<-MAGICAL_initialization(loaded_data, Candidate_circuits)

Circuits_linkage_posterior<-MAGICAL_estimation(loaded_data, Candidate_circuits, Initial_model, iteration_num = 1000)

Output_path <- paste0("/banach2/chang/new_SCORCH_analysis/MAGICAL/output/",opt$region,"_",opt$celltype,"MAGICAL_selected_regulatory_circuits.txt")

MAGICAL_circuits_output(Output_file_path = Output_path, 
                        Candidate_circuits, Circuits_linkage_posterior)