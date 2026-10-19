# REVISION: run on the CellBender-corrected data; paths (including fast_tsne_path) come
# from revision/config/paths.R and outputs go to OUT_ROOT. Changes from the published
# script are marked ## REVISION.

library(Seurat)
## REVISION: packages used but never loaded in the published script (they were
##   loaded earlier in the interactive session); under Rscript this is an error.
library(ggplot2)
library(dplyr)
## REVISION: also used but never loaded.
library(cowplot)
library(ggpubr)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup("01-Regional_data_processing",
            subdirs = c("combined/data", "combined/figures", "combined/figures/tsne"))
seurat.combined <- readRDS("./combined/data/seurat_combined_v2.RDS")

# re-integrate the data
## REVISION: seurat_combined_v2.RDS arrives with the RNA assay still split
##   (183 layers - combined/01 saved it split and combined/02 only subsets it),
##   and Seurat v5 refuses to split an already-split assay:
##     "Please join before splitting"
##   The published script ran this line against a joined object. JoinLayers()
##   is a no-op on a joined assay, so this is correct either way.
seurat.combined <- JoinLayers(seurat.combined)
seurat.combined[["RNA"]] <- split(seurat.combined[["RNA"]], f = seurat.combined$sample_name)
seurat.combined <- NormalizeData(seurat.combined)
seurat.combined <- FindVariableFeatures(seurat.combined)
seurat.combined <- ScaleData(seurat.combined)
seurat.combined <- RunPCA(seurat.combined)
seurat.combined <- IntegrateLayers(
  object = seurat.combined, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
# visualization
seurat.combined <-  RunTSNE(
  seurat.combined,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path=FAST_TSNE_PATH, ann_not_vptree=FALSE, nthreads=36
)

## REVISION (reviewer request): UMAP alongside the tSNE above.
##   Same harmony reduction and dims, so the two embeddings are directly
##   comparable. Stored on the object before it is saved. Not part of the
##   published analysis.
##
##   Threads: Seurat passes n_threads = future::nbrOfWorkers() to uwot, so a
##   direct n_threads argument to RunUMAP is silently dropped. run_umap_threaded()
##   sets the future plan for the call and restores it afterwards.
seurat.combined <- run_umap_threaded(seurat.combined, reduction = "harmony", dims = 1:50)

## REVISION: checkpoint the object as soon as the expensive work is done.
##   Everything above this line is ~2 hours (JoinLayers, Harmony, FIt-SNE,
##   UMAP); everything below is plotting. Writing the object before any plot is
##   drawn means a plotting failure costs only the plots. The published saveRDS
##   further down is left in place; it rewrites the same path with the
##   annotation added later in the script.
saveRDS(seurat.combined, "./combined/data/seurat_combined_v3.RDS")
message("[combined/03] checkpoint written: seurat_combined_v3.RDS (",
        ncol(seurat.combined), " cells, reductions: ",
        paste(names(seurat.combined@reductions), collapse = ", "), ")")

## REVISION: write named tSNE and UMAP panels. The published scripts call
##   DimPlot() with no device open, so under Rscript every plot lands unlabelled
##   in Rplots.pdf and the embedding figures cannot be found.
save_embedding_panels(seurat.combined, prefix = "03_integration_tsne_refined")


region_cols <- c("#E69F00", "#56B4E9","#CC79A7")
names(region_cols) <-c("INS","PFC","VST")
#region_cols <- factor(region_cols,levels = c("#56B4E9","#E69F00","#CC79A7"))
t1 <- DimPlot(seurat.combined,reduction = "tsne",group.by = "region",
              cols = region_cols,
              label=F,repel = T,
              raster = F,shuffle = T,label.size = 10,pt.size = 0.5) +
  theme_void()+
  NoLegend()+
  theme(title = element_blank())



t1_legend <- as_ggplot(get_legend(DimPlot(seurat.combined,reduction = "tsne",group.by = "region",
                                          cols = region_cols,
                                          order = c("VST","INS","PFC"),
                                          label=F,
                                          raster = F,shuffle = T,label.size = 5,pt.size = 0.5) +
                                    theme_void()+
                                    guides(color = guide_legend(ncol = 1,override.aes = list(size=5)))+
                                    
                                    theme(title = element_blank(),legend.text = element_text(size=20))))


cond_col <- c("#0072B2", "darkorange1")
cond_plot <- DimPlot(seurat.combined,reduction = "tsne",group.by = "condition_new",
                     label=F,
                     cols = cond_col,
                     raster = F,shuffle = T,label.size = 10,pt.size = 0.5)+
  theme_void()+
  NoLegend()+
  theme(title = element_blank())

cond_plot_legend <- as_ggplot(get_legend(DimPlot(seurat.combined,reduction = "tsne",group.by = "condition_new",
                                                 label=F,
                                                 cols = cond_col,
                                                 raster = F,shuffle = T,label.size = 10,pt.size = 0.5)+
                                           theme_void()+
                                           guides(color = guide_legend(ncol = 1,override.aes = list(size=5)))+
                                           theme(title = element_blank(),legend.text = element_text(size=20))))


t1_legend <- as_ggplot(get_legend(DimPlot(seurat.combined,reduction = "tsne",group.by = "region",
                                          cols = region_cols,
                                          order = c("VST","INS","PFC"),
                                          label=F,
                                          raster = F,shuffle = T,label.size = 5,pt.size = 0.5) +
                                    theme_void()+
                                    guides(color = guide_legend(ncol = 1,override.aes = list(size=5)))+
                                    
                                    theme(title = element_blank(),legend.text = element_text(size=20))))


cond_col <- c("#0072B2", "darkorange1")
cond_plot <- DimPlot(seurat.combined,reduction = "tsne",group.by = "condition_new",
                     label=F,
                     cols = cond_col,
                     raster = F,shuffle = T,label.size = 10,pt.size = 0.5)+
  theme_void()+
  NoLegend()+
  theme(title = element_blank())

cond_plot_legend <- as_ggplot(get_legend(DimPlot(seurat.combined,reduction = "tsne",group.by = "condition_new",
                                                 label=F,
                                                 cols = cond_col,
                                                 raster = F,shuffle = T,label.size = 10,pt.size = 0.5)+
                                           theme_void()+
                                           guides(color = guide_legend(ncol = 1,override.aes = list(size=5)))+
                                           theme(title = element_blank(),legend.text = element_text(size=20))))

png(file="./combined/figures/tsne/tsne_condition_new.png",width = 3500,height = 3500,res=300)
print(cond_plot)
dev.off()

png(file="./combined/figures/tsne/tsne_condition_legend_new.png",width = 300,height = 300,res=300)
print(cond_plot_legend)
dev.off()

png(file="./combined/figures/tsne/tsne_region_new.png",width = 3500,height = 3500,res=300)
print(t1)
dev.off()

png(file="./combined/figures/tsne/tsne_region_legend_new.png",width = 300,height = 300,res=300)
print(t1_legend)
dev.off()




metadata_all <- seurat.combined@meta.data
# number of cells in each condition in each region
n_counts <- metadata_all %>%
  group_by(region,condition_new) %>%
  summarise(n=n()) 

n_count_sample <- metadata_all %>%
  group_by(sample_ids,region) %>%
  summarise(n=n()) 

n_counts$condition_new <- factor(n_counts$condition_new,levels = c("HIV","CTR"))
n_counts$region <- factor(n_counts$region,levels = c("VST","INS","PFC"))
names(cond_col) <- c("CTR","HIV")
p1 <- ggplot(n_counts, aes(x=region,y=n,fill=condition_new)) +
  scale_fill_manual(values = cond_col)+
  #facet_wrap(~group_id, scales = "free_x")+
  geom_bar(stat="identity",position=position_dodge()) +
  geom_text(aes(label=n), hjust=-0.05,position = position_dodge(width=1),size = 20) +
  ylab("# Cells")+
  xlab("")+
  theme_cowplot(font_size = 30) +
  #scale_x_discrete(position = "top") +
  #scale_y_continuous(trans='sqrt',limits = c(0,150000))+
  ylim(c(0,180000))+
  theme(axis.text = element_text(size=60),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        
        strip.text = element_text(size=15),
        legend.title=element_text(size=20),
        title=element_text(size=60),
        legend.text=element_text(size=10)) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+NoLegend() + coord_flip()

png(file="./combined/figures/n_counts_new.png",width = 4500,height = 4500,res=300)
print(p1)
dev.off()

saveRDS(seurat.combined,"./combined/data/seurat_combined_v3.RDS")


seurat.combined <- readRDS("./combined/data/seurat_combined_v3.RDS")

# get legend


cp_col1 <- c(
  "indianred1", "tomato1","sienna1","sienna3","maroon",  "#6A3D9A", "#FB9A99", "#CAB2D6", "orchid1", 
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "purple1", "lightslateblue","orchid3","pink3","plum1",  "tomato1", "tomato3","sienna2", "burlywood3", 
  
  "#FF7F00", "#FDBF6F","khaki2", "gold1","gold3",
  "olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
cell_types_consensus_levels <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                 "Astro","Oligo","OPC","Micro","Macrophages","Other immune","PC","Endo","SMC","VLMC")


names(cp_col1) <- cell_types_consensus_levels
seurat.combined@meta.data$consensus_annotation_v3_relabeled <- seurat.combined@meta.data$consensus_annotation_v3
seurat.combined@meta.data$consensus_annotation_v3_relabeled[seurat.combined@meta.data$consensus_annotation_v3 == "Immune"] = "Other immune"
seurat.combined@meta.data$consensus_annotation_v3_relabeled <- factor(seurat.combined@meta.data$consensus_annotation_v3_relabeled,
                                                                      levels = cell_types_consensus_levels)
t1_legend <- as_ggplot(get_legend(DimPlot(seurat.combined,reduction = "tsne",group.by = "consensus_annotation_v3_relabeled",
                                          cols = cp_col1,
                                          label=F,
                                          raster = F,shuffle = T,label.size = 5,pt.size = 0.5) +
                                    theme_void()+
                                    guides(color = guide_legend(nrow = 3,override.aes = list(size=10)))+
                                    theme(title = element_blank(),legend.text = element_text(size=20))))
png(file="./combined/figures/legend_all_cp.png",width = 8000,height = 2000,res=300)
print(t1_legend)
dev.off()
