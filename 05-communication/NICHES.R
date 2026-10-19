# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Seurat)
library(dplyr)
library(ggplot2)
library(NICHES)
source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT by rerun_setup() instead of a
##   hard-coded setwd().
rerun_setup("",
            subdirs = c("01-Regional_data_processing/combined/data", "08-communication/NICHES/figures_new"))
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
                                              ## REVISION: "astro_subtype" removed from meta.data.to.map.
                                              ##   The column does not exist on the object, so RunNICHES()
                                              ##   errors, and nothing downstream reads it.
                                              species = "human",meta.data.to.map = c("condition_new","region"),
                                              cell_types ="consensus_annotation_v3",CellToCell = T )
    niches_obj_list[[sample_id]] <- niches_obj_list[[sample_id]][["CellToCell"]]

    ## REVISION: carry donor identity onto the interaction vectors before the
    ##   merge. RunNICHES() maps only the meta.data.to.map columns, none of which
    ##   identify the donor, so donor-level statistics would otherwise be impossible.
    .md <- test_sample_obj_list[[sample_id]]@meta.data
    niches_obj_list[[sample_id]]@meta.data$sample_name <- sample_id
    if ("sample_ids" %in% colnames(.md))
      niches_obj_list[[sample_id]]@meta.data$sample_ids <- as.character(.md$sample_ids)[1]
    else
      warning("no sample_ids on ", sample_id, " - donor-level NICHES DE will ",
              "have to fall back to sample_name, which splits a donor by region")
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

## REVISION: person-first labels (PWH / PWoH instead of HIV / CTR). Everything
##   downstream that keys on them (vec_condition, the FindMarkers idents, the
##   DotPlot split) is updated together.
merge_niches_obj@meta.data$condition <- as.character(pwh_relabel(merge_niches_obj@meta.data$condition_new.Receiving))
merge_niches_obj@meta.data$vec_condition <- paste0(merge_niches_obj@meta.data$VectorType,'|',merge_niches_obj@meta.data$condition)

## REVISION: save the merged object (the published script kept only the figures),
##   so donor-level reanalysis is possible.
dir.create("./08-communication/NICHES", recursive = TRUE, showWarnings = FALSE)
saveRDS(merge_niches_obj, "./08-communication/NICHES/merge_niches_obj.RDS")
message("[niches] wrote merge_niches_obj.RDS: ", ncol(merge_niches_obj),
        " interaction vectors, ",
        length(unique(merge_niches_obj@meta.data$VectorType)), " VectorTypes, ",
        if ("sample_ids" %in% colnames(merge_niches_obj@meta.data))
          paste0(length(unique(merge_niches_obj@meta.data$sample_ids)), " donors")
        else "NO donor column")

region_cols <- c("#E69F00", "#56B4E9","#CC79A7")
names(region_cols) <-c("INS","PFC","VST")

cond_col <- PWH_COLOURS   ## REVISION: person-first labels; PWoH = #0072B2, PWH = darkorange1


Idents(merge_niches_obj) <- "vec_condition"
p_list_up <- list()
p_list_down <- list()

merge_niches_obj_sub_list <- list()
## REVISION: keep the FindMarkers output for every VectorType. The published
##   script used it only to pick 20 interactions per DotPlot and did not save the
##   ranked ligand-receptor pairs.
marker_list <- list()

for(vec in unique(merge_niches_obj@meta.data$VectorType)){

  markers <- FindMarkers(merge_niches_obj,
                         ident.1 = paste0(vec,"|PWH"),ident.2 = paste0(vec,"|PWoH")
  )
  marker_list[[vec]] <- tibble::rownames_to_column(markers, "interaction")
  marker_list[[vec]]$VectorType <- vec
  markers.up.sig <- markers %>%
    filter(p_val_adj < 0.01)%>%
    slice_max(avg_log2FC,n=20) 
  markers.down.sig <- markers %>%
    filter(p_val_adj < 0.01)%>%
    slice_max(-avg_log2FC,n=20) 
  merge_niches_obj_sub <- subset(merge_niches_obj,
                                 cells = colnames(merge_niches_obj)[Idents(merge_niches_obj) %in% c(paste0(vec,"|PWH"),paste0(vec,"|PWoH"))])
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

## REVISION: write the ranked interactions (per VectorType RDS and one combined
##   CSV) so the full result is available, not only the dot plots.
saveRDS(marker_list, "./08-communication/NICHES/niches_markers_by_vectortype.RDS")
all_markers <- do.call(rbind, marker_list)
all_markers <- all_markers[order(all_markers$p_val_adj, -abs(all_markers$avg_log2FC)), ]
write.csv(all_markers, "./08-communication/NICHES/niches_markers_all.csv", row.names = FALSE)
message("[niches] wrote ", nrow(all_markers), " interaction tests across ",
        length(marker_list), " VectorTypes; ",
        sum(all_markers$p_val_adj < 0.01, na.rm = TRUE), " at p_adj < 0.01")

for(vec in unique(merge_niches_obj@meta.data$VectorType)){

  png(file=paste0("./08-communication/NICHES/figures_new/",make.names(vec),".png"),width = 2000,height = 2000,res=300)
  print(p_list_up[[vec]])
  dev.off()

  ## REVISION: the published script built p_list_down but never wrote it; it is
  ##   now saved alongside the up panel.
  png(file=paste0("./08-communication/NICHES/figures_new/",make.names(vec),"_down.png"),width = 2000,height = 2000,res=300)
  print(p_list_down[[vec]])
  dev.off()
}

