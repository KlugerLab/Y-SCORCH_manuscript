# get the tsne from the new combined object
#
# REVISION: run on the CellBender-corrected data; paths come from revision/config/paths.R
# and outputs go to OUT_ROOT. Changes from the published script are marked ## REVISION.

library(Seurat)
## REVISION: packages used but never loaded in the published script (they were
##   loaded earlier in the interactive session); under Rscript this is an error.
library(ggplot2)
library(patchwork)
## REVISION: also used but never loaded.
library(cowplot)
library(dplyr)
library(ggpubr)

source(file.path(Sys.getenv("SCORCH_REPO", "."), "revision/config/paths.R"))   ## REVISION

## REVISION: working directory is set under OUT_ROOT instead of a hard-coded server path.
rerun_setup("01-Regional_data_processing",
            subdirs = c("VST/data/04-VST-tsne_marker_plotting", "VST/figures/marker_and_abundance", "VST/figures/tsne_new", "combined/data"))
## REVISION: prefer the combined object with donor 1040's VST cells removed.
##   1040's VST library is excluded under the post-ambient-correction rule
##   (38.1% excitatory, the smallest VST library); see VST_SAMPLE_KEEP_OVERRIDE
##   in revision/config/paths.R.
##
##   Resolution order: SCORCH_COMBINED_OBJ env override, then the filtered
##   object, then the unfiltered one. The unfiltered fallback lets the script
##   run on a fresh checkout, but it raises a warning, since a silent fallback
##   to the object containing 1040 is the failure this guards against.
## Order matters. The _msn object is the only one carrying the nine MSN
## subtypes (revision/tools/patch_msn_subtypes.R); without it every striatal
## neuron is a single "MSNs" class, cp_col1 below has no entry for it, and the
## fine-grained tSNE renders the largest VST population unmapped.
.combined_candidates <- c(
  Sys.getenv("SCORCH_COMBINED_OBJ", unset = NA_character_),
  "./combined/data/seurat_combined_v3_no1040_msn.RDS",
  "./combined/data/seurat_combined_v3_no1040.RDS",
  "./combined/data/seurat_combined_v3.RDS"
)
.combined_in <- .combined_candidates[!is.na(.combined_candidates) &
                                     file.exists(.combined_candidates)][1]
if (is.na(.combined_in)) stop("no combined object found; tried: ",
                              paste(na.omit(.combined_candidates), collapse = ", "))
if (!grepl("no1040", .combined_in))
  warning("[VST] using ", .combined_in, " - this object still contains donor ",
          "1040's VST cells. Run revision/tools/combined_drop_1040.R first unless that ",
          "is deliberate.", immediate. = TRUE)
message("[VST] combined object: ", .combined_in)
seurat.combined <- readRDS(.combined_in)

seurat.vst <- subset(seurat.combined,cells = colnames(seurat.combined)[seurat.combined@meta.data$region == "VST"])
message("[VST] VST cells: ", ncol(seurat.vst), " from ",
        length(unique(as.character(seurat.vst@meta.data$sample_ids))), " samples")

## Report whether the subtypes are present. cp_col1 is keyed to the nine of
## them; a collapsed "MSNs" class silently drops out of the palette instead of
## erroring, so the figure looks finished and is not.
.ann <- as.character(seurat.vst@meta.data$consensus_annotation_v3)
.n_msn_sub <- length(unique(.ann[grepl("^D1|^D2", .ann)]))
message("[VST] MSN subtypes present: ", .n_msn_sub, " of 9",
        if (any(.ann == "MSNs")) paste0("  (", sum(.ann == "MSNs"),
                                        " cells still labelled 'MSNs')") else "")
if (.n_msn_sub == 0)
  warning("[VST] no MSN subtypes on this object - the fine-grained tSNE will ",
          "render every striatal neuron unmapped. Expected ",
          "seurat_combined_v3_no1040_msn.RDS.", immediate. = TRUE)

seurat.vst <- NormalizeData(seurat.vst)
seurat.vst <- FindVariableFeatures(seurat.vst)
seurat.vst <- ScaleData(seurat.vst)
seurat.vst <- RunPCA(seurat.vst)
seurat.vst <- IntegrateLayers(
  object = seurat.vst, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat.vst <-  RunTSNE(
  seurat.vst,reduction = "harmony", tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:50, perplexity = 200,
  fast_tsne_path=FAST_TSNE_PATH, ann_not_vptree=FALSE, nthreads=36#,initialization=seurat.combined@reductions$pca@cell.embeddings[colnames(seurat.vst),]
)

## REVISION (reviewer request): UMAP alongside the tSNE above.
##   Same harmony reduction and dims, so the two embeddings are directly
##   comparable. Stored on the object before it is saved. Not part of the
##   published analysis.
##
##   Threads: Seurat passes n_threads = future::nbrOfWorkers() to uwot, so a
##   direct n_threads argument to RunUMAP is silently dropped. run_umap_threaded()
##   sets the future plan for the call and restores it afterwards.
seurat.vst <- run_umap_threaded(seurat.vst, reduction = "harmony", dims = 1:50)

## REVISION: write named tSNE and UMAP panels. The published scripts call
##   DimPlot() with no device open, so under Rscript every plot lands unlabelled
##   in Rplots.pdf and the embedding figures cannot be found.
save_embedding_panels(seurat.vst, prefix = "05_VST_tsne_from_all")

# cell types (fine grained)
cp_col1 <- c(
  "purple1", "lightslateblue","orchid3","pink3","plum1",  "tomato1", "tomato3","sienna2", "burlywood3", 
  "green4","palegreen2", "darkseagreen","dodgerblue2","skyblue3", "darkturquoise", "skyblue1", "steelblue4","deepskyblue",
  "#FF7F00", "#FDBF6F","khaki2", "gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3"
  
)
cell_types_consensus_levels <- c("D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid",
                                 "Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip",
                                 "Astro","Oligo","OPC","Micro","Macrophages","Other immune","PC","Endo","SMC","VLMC")


names(cp_col1) <- cell_types_consensus_levels
seurat.vst@meta.data$consensus_annotation_v3_relabeled <- seurat.vst@meta.data$consensus_annotation_v3
seurat.vst@meta.data$consensus_annotation_v3_relabeled[seurat.vst@meta.data$consensus_annotation_v3 == "Immune"] = "Other immune"
seurat.vst <- JoinLayers(seurat.vst)
saveRDS(seurat.vst,"./VST/data/04-VST-tsne_marker_plotting/integrated_seurat_obj_v3.RDS")

seurat.vst <- readRDS("./VST/data/04-VST-tsne_marker_plotting/integrated_seurat_obj_v3.RDS")


t1 <- DimPlot(seurat.vst,reduction = "tsne",group.by = "consensus_annotation_v3_relabeled",
              cols = cp_col1,
              label=T,repel = T,
              raster = F,shuffle = T,label.size = 10,pt.size = 0.5) +
  theme_void()+
  NoLegend()+
  theme(title = element_blank())

t1_legend <- as_ggplot(get_legend(DimPlot(seurat.vst,reduction = "tsne",group.by = "consensus_annotation_v3_relabeled",
                                          cols = cp_col1,
                                          label=F,
                                          raster = F,shuffle = T,label.size = 5,pt.size = 0.5) +
                                    theme_void()+
                                    guides(color = guide_legend(ncol = 3))+
                                    
                                    theme(title = element_blank(),legend.text = element_text(size=20))))

# cell types (coarse grained)
cp_col2 <-  c("tomato2","royalblue3","dodgerblue","gold2")
names(cp_col2) <- c("Ext","MSNs","Inh","Glial")

t1_coarse <- DimPlot(seurat.vst,reduction = "tsne",group.by = "consensus_cell_type_coarse",
                     cols = cp_col2,
                     label=F,
                     raster = F,shuffle = T,label.size = 5,pt.size = 0.5) +
  theme_void()+
  NoLegend()+
  theme(title = element_blank())
t1_coarse_legend <- as_ggplot(get_legend(DimPlot(seurat.vst,reduction = "tsne",group.by = "consensus_cell_type_coarse",
                                                 cols = cp_col2,
                                                 label=F,
                                                 raster = F,shuffle = T,label.size = 5,pt.size = 0.5) +
                                           theme_void()+
                                           
                                           theme(title = element_blank(),legend.text = element_text(size=20))))

manual_markers <- c("GAD2","BCL11B","AQP4","MBP","PDGFRA","P2RY12","F13A1",
                    "SKAP1","CLDN5","GRM8","ACTA2","COL1A2")

# stacked violin
vst_data <- seurat.vst@assays$RNA@layers$data
rownames(vst_data) <- rownames(seurat.vst)
colnames(vst_data) <- colnames(seurat.vst)
stacked_violin_data <- as.data.frame(t(vst_data[manual_markers[manual_markers %in% rownames(vst_data)],]))
stacked_violin_data$Cell <- rownames(stacked_violin_data)

stacked_violin_data$cell_types <- seurat.vst@meta.data$consensus_annotation_v3_relabeled
ext_cps <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L6b","L5 ET","L5/6 NP","L6 CT")
inh_cps <- c("Lamp5","Lamp5 Lhx6","Pax6","Pvalb","Chandelier","Sncg","Sst","Sst Chodl","Vip")
msns_cps <- c("D1-Matrix","D1-Shell/OT","D1-Striosome","D1-ICj","D1-NUDAP","D2-Matrix","D2-Shell/OT","D2-Striosome","D1/D2-Hybrid")
stacked_violin_data$cell_types[stacked_violin_data$cell_types %in% ext_cps] <- "Ext"
stacked_violin_data$cell_types[stacked_violin_data$cell_types %in% inh_cps] <- "Inh"
stacked_violin_data$cell_types[stacked_violin_data$cell_types %in% msns_cps] <- "MSNs"

stacked_violin_data$cell_types <- factor(stacked_violin_data$cell_types,
                                         levels = c("Ext","Inh","MSNs",
                                                    "Astro","Oligo","OPC",
                                                    "Micro","Macrophages","Other immune",
                                                    "Endo","PC","SMC","VLMC"))
seurat.vst@meta.data$consensus_cell_type_coarse_relabeled <- stacked_violin_data$cell_types
#stacked_violin_data$cell_types <- seurat.ins@meta.data$consensus_cell_type_coarse_l1
cp_col3 <- c("tomato2","dodgerblue","royalblue3","#FF7F00", "#FDBF6F","khaki2", 
             "gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3")

names(cp_col3) <- levels(stacked_violin_data$cell_types)

# Use melt to change data.frame format
stacked_violin_data <- reshape2::melt(stacked_violin_data, id.vars = c("Cell","cell_types"), 
                                      measure.vars = manual_markers[manual_markers %in% rownames(vst_data)],
                                      variable.name = "Feat", value.name = "Expr")

# add a fake cell from Ext
#stacked_violin_data[nrow(stacked_violin_data)+1,] <- list("tmp","Ext","SATB2",0)
v1 <- ggplot(stacked_violin_data, aes(Expr, factor(Feat), fill = cell_types)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(cols = vars(cell_types), scales = "free") +
  theme_cowplot(font_size = 12) +
  scale_fill_manual(values=cp_col3)+
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        axis.text.y = element_text(size=25), 
        axis.title = element_text(size=25), 
        
        legend.title=element_text(size=20)
  )  + xlab("Expression Level") + ylab("Marker genes")


# add on top the proportions of samples in each cell type

metadata_all <- seurat.vst@meta.data
n_counts <- metadata_all %>%
  group_by(consensus_cell_type_coarse_relabeled,sample_ids) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  group_by(consensus_cell_type_coarse_relabeled) %>%
  mutate(n_cp = sum(n)) %>%
  mutate(percent = n/n_cp) %>%
  mutate(cell_types = consensus_cell_type_coarse_relabeled)

n_counts_total <- n_counts %>%
  dplyr::select(cell_types,n_cp) %>%
  distinct() 
new_levels <- levels(stacked_violin_data$cell_types)
new_levels[new_levels == "Macrophages"] <- "Macro"
new_levels[new_levels == "Other immune"] <- "Other\nimmune"

n_counts_total$cell_types <- as.character(n_counts_total$cell_types)
n_counts_total$cell_types[n_counts_total$cell_types == "Macrophages"] <- "Macro"
n_counts_total$cell_types[n_counts_total$cell_types == "Other immune"] <- "Other\nimmune"

n_counts_total$cell_types <- factor(n_counts_total$cell_types,levels = new_levels)
#n_counts_total[13,] <- list("Ext","Ext",0)

cp_col4 <- c("tomato2","dodgerblue","royalblue3","#FF7F00", "#FDBF6F","khaki2", 
             "gold1","gold3","olivedrab2","darkorange3","goldenrod","yellow4", "yellow3")

names(cp_col4) <- new_levels


p2 <- ggplot(n_counts_total, aes(x=cell_types,y=n_cp,fill=cell_types)) +
  scale_fill_manual(values = cp_col4)+
  #facet_wrap(~group_id, scales = "free_x")+
  geom_bar(stat="identity") +
  ylab("# Cells")+
  theme_cowplot(font_size = 15) +
  scale_x_discrete(position = "top") +
  scale_y_continuous(trans='sqrt')+
  theme(axis.text.y = element_text(size=25),
        axis.text.x = element_text(angle = 0,size=20),
        axis.title.x = element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size=20),
        legend.title=element_text(size=25),
        title=element_text(size=25),
        legend.text=element_text(size=15)) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+NoLegend()


#ggpubr::ggarrange(p2,v1,ncol=1)
marker_and_ncells <- p2 / v1 +   plot_layout(heights = c(1, 3))


png(file="./VST/figures/tsne_new/tsne_cp_coarse.png",width = 3500,height = 3500,res=300)
print(t1_coarse)
dev.off()

png(file="./VST/figures/tsne_new/tsne_cp_coarse_legend.png",width = 500,height = 500,res=300)
print(t1_coarse_legend)
dev.off()


png(file="./VST/figures/tsne_new/tsne_cp.png",width = 3500,height = 3500,res=300)
print(t1)
dev.off()

png(file="./VST/figures/tsne_new/tsne_cp_legend.png",width = 1500,height = 1500,res=300)
print(t1_legend)
dev.off()

png(file="./VST/figures/marker_and_abundance/marker_and_ncells_new.png",width = 4000,height = 4000,res=300)
print(marker_and_ncells)
dev.off()

