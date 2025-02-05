library(Seurat)
library(dplyr)
library(ggplot2)
library(NICHES)
setwd("/banach2/SCORCH/scripts/combined_analysis_paper/")


seurat.combined <- readRDS("./01-Regional_data_processing/combined/data/seurat_combined_v3_layer_joined.RDS")



seurat.combined.sub <- subset(seurat.combined,
                              cells = colnames(seurat.combined)[(seurat.combined@meta.data$consensus_annotation_v3 %in% c("L2/3 IT","L4 IT","Pvalb", "Micro","Astro","Endo"))])

cps <- c("Micro","Astro","Endo","L2/3 IT","L4 IT","Pvalb")

test_sample_obj_list <- list()
for (sample_id in unique(seurat.combined.sub@meta.data$sample_name)) {
  test_sample_obj_list[[sample_id]] <- subset(seurat.combined.sub,cells = colnames(seurat.combined.sub)[seurat.combined.sub@meta.data$sample_name == sample_id])
}

# run niches
niches_obj_list <- list()
# only run on the ones with more than 500 cells
sample2keep <- names(test_sample_obj_list)
for(sample_id in names(test_sample_obj_list)){
  if(ncol(test_sample_obj_list[[sample_id]]) > 500){
    Idents(test_sample_obj_list[[sample_id]]) <- "consensus_annotation_v3"
    niches_obj_list[[sample_id]] <- RunNICHES(test_sample_obj_list[[sample_id]],
                                              assay = "RNA",LR.database = "fantom5",min.cells.per.ident = 10,
                                              species = "human",meta.data.to.map = c("astro_subtype","condition_new","region"),
                                              cell_types ="consensus_annotation_v3",CellToCell = T )
    niches_obj_list[[sample_id]] <- niches_obj_list[[sample_id]][["CellToCell"]]
    if(length(unique(niches_obj_list[[sample_id]]@meta.data$VectorType)) > 1){
      
      if(ncol(niches_obj_list[[sample_id]]) > 200){
        sample2keep <- c(sample2keep,sample_id)
        niches_obj_list[[sample_id]] <- FindVariableFeatures(niches_obj_list[[sample_id]],selection.method = "disp")
        niches_obj_list[[sample_id]] <- ScaleData(niches_obj_list[[sample_id]])
        niches_obj_list[[sample_id]] <- RunPCA(niches_obj_list[[sample_id]])
      }
      
      
    }}
  
}

# merge the two and do DE
sample2keep <- names(niches_obj_list)
count1 <- niches_obj_list[[sample2keep[1]]]@assays$CellToCell@layers$counts
colnames(count1) <- colnames(niches_obj_list[[sample2keep[1]]])
rownames(count1) <- rownames(niches_obj_list[[sample2keep[1]]])

for(sample_id in sample2keep[2:length(sample2keep)]){
  count2 <- niches_obj_list[[sample_id]]@assays$CellToCell@layers$counts
  colnames(count2) <- colnames(niches_obj_list[[sample_id]])
  rownames(count2) <- rownames(niches_obj_list[[sample_id]])
  count1 <- cbind(count1,count2)
  
}


merged_meta1 <- do.call(rbind,lapply(sample2keep,function(sid)niches_obj_list[[sid]]@meta.data))


merge_niches_obj <- CreateSeuratObject(counts = count1,
                                       meta.data = merged_meta1)
merge_niches_obj <- NormalizeData(merge_niches_obj)
merge_niches_obj@assays$RNA@layers$data <- merge_niches_obj@assays$RNA@layers$counts

merge_niches_obj@meta.data$condition <- as.character(merge_niches_obj@meta.data$condition_new.Receiving)
merge_niches_obj@meta.data$vec_condition <- paste0(merge_niches_obj@meta.data$VectorType,'|',merge_niches_obj@meta.data$condition)

region_cols <- c("#E69F00", "#56B4E9","#CC79A7")
names(region_cols) <-c("INS","PFC","VST")

cond_col <- c("#0072B2", "darkorange1")
names(cond_col) <- c("CTR","HIV")


Idents(merge_niches_obj) <- "vec_condition"
p_list_up <- list()
p_list_down <- list()

merge_niches_obj_sub_list <- list()

for(vec in unique(merge_niches_obj@meta.data$VectorType)){
  
  markers <- FindMarkers(merge_niches_obj,
                         ident.1 = paste0(vec,"|HIV"),ident.2 = paste0(vec,"|CTR")
  )
  markers.up.sig <- markers %>%
    filter(p_val_adj < 0.01)%>%
    slice_max(avg_log2FC,n=20) 
  markers.down.sig <- markers %>%
    filter(p_val_adj < 0.01)%>%
    slice_max(-avg_log2FC,n=20) 
  merge_niches_obj_sub <- subset(merge_niches_obj,
                                 cells = colnames(merge_niches_obj)[Idents(merge_niches_obj) %in% c(paste0(vec,"|HIV"),paste0(vec,"|CTR"))])
  merge_niches_obj_sub <- ScaleData(merge_niches_obj_sub,features = rownames(merge_niches_obj_sub))
  
  
  p_list_up[[vec]] <- DotPlot(merge_niches_obj_sub,
                              features = rownames(markers.up.sig),
                              group.by = "region.Receiving",cols = "RdBu",
                              split.by = "condition")+ coord_flip()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p_list_down[[vec]] <- DotPlot(merge_niches_obj_sub,
                                features = rownames(markers.down.sig),
                                group.by = "region.Receiving",cols = "RdBu",
                                split.by = "condition")+ coord_flip()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
}

for(vec in unique(merge_niches_obj@meta.data$VectorType)){
  
  png(file=paste0("./08-communication/NICHES/figures_new/",make.names(vec),".png"),width = 2000,height = 2000,res=300)
  print(p_list_up[[vec]])
  dev.off()
}

