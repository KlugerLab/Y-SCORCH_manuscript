library(DESeq2)
library(muscat)
wd <- "/banach2/SCORCH/scripts/combined_analysis_paper/03-DE_analysis/"
setwd(wd)


META_PATH = "/banach2/SCORCH/data/metadata/sample_metadata_v3.csv"
sample_meta <- read.csv(META_PATH)

sample_source <- sample_meta %>%
  dplyr::select(sample_ids,subject_source_v2) %>%
  distinct()
rownames(sample_source) <- sample_source$sample_ids



pb <- readRDS("./cross_region/data/pseudo_bulk.rds")
mm <- readRDS("./cross_region/data/mm.rds")
cell_types <- names(pb@assays@data)

pb@colData$sample_source <- sample_source[pb@colData$sample_ids,]$subject_source_v2

condition_cols=c("#0072B2","darkorange1")

p_list_all <- list()
for(cp in  c("Micro")){
  p_list_all[[cp]] <- list()
  
  for(region in c("INS","PFC","VST")){
    
    n_cells_df <- as.data.frame(do.call(rbind,pb@int_colData$n_cells))
    s2keep <- rownames(n_cells_df)[(n_cells_df[,cp]>20)]
    pb_sub <- pb[,s2keep]
    
    s2keep <- colSums(pb_sub@assays@data[[cp]]) > 20
    pb_sub <- pb_sub[,s2keep]
    
    pb_sub <- pb_sub[,pb_sub@colData$region == region]
    
    dds <- DESeqDataSetFromMatrix(pb_sub@assays@data[[cp]], 
                                  colData = pb_sub@colData, 
                                  design = ~condition_new)
    vsd <- vst(dds, blind=F)
    
    pcaData <- plotPCA(vsd, intgroup=c("condition_new","sample_source"), returnData=TRUE,ntop=300)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    p1 <- ggplot(pcaData, aes(PC1, PC2, color=condition_new,label=name)) +
      geom_point(size=5) +
      scale_color_manual(values = condition_cols)+
      #geom_text(size=3)+
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      guides(color=guide_legend(title="Condition"))
    p1 <- p1+
      theme_minimal()+
      theme(legend.title=element_text(size=20),
            axis.text = element_text(size=20),
            title=element_text(size=20),
            legend.text=element_text(size=20))
    
    
    p_list_all[[cp]][[region]] <- p1
  }
  
}

library(patchwork)

png(filename = paste0("./plotting/figures/micro_stacked_condition.png"),width=1500,height = 3000,res = 300)
print(p_list_all[[1]][["PFC"]]+NoLegend()+
        p_list_all[[1]][["INS"]]+NoLegend()+
        p_list_all[[1]][["VST"]]+NoLegend()+ plot_layout(ncol = 1)
)
dev.off()





p_list_all <- list()
for(cp in  c("Micro")){
  p_list_all[[cp]] <- list()
  
  for(region in c("INS","PFC","VST")){
    
    n_cells_df <- as.data.frame(do.call(rbind,pb@int_colData$n_cells))
    s2keep <- rownames(n_cells_df)[(n_cells_df[,cp]>20)]
    
    pb_sub <- pb[,s2keep]
    
    s2keep <- colSums(pb_sub@assays@data[[cp]]) > 20
    pb_sub <- pb_sub[,s2keep]
    
    pb_sub <- pb_sub[,pb_sub@colData$region == region]
    
    dds <- DESeqDataSetFromMatrix(pb_sub@assays@data[[cp]], 
                                  colData = pb_sub@colData, 
                                  design = ~condition_new)
    vsd <- vst(dds, blind=F)
    
    pcaData <- plotPCA(vsd, intgroup=c("condition_new","sample_source"), returnData=TRUE,ntop=300)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    p1 <- ggplot(pcaData, aes(PC1, PC2, color=sample_source,label=name)) +
      geom_point(size=5) +
      #scale_color_manual(values = condition_cols)+
      #geom_text(size=3)+
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      guides(color=guide_legend(title="Sources"))
    p1 <- p1+
      theme_minimal()+
      theme(legend.title=element_text(size=20),
            axis.text = element_text(size=20),
            title=element_text(size=20),
            legend.text=element_text(size=20))
    
    
    p_list_all[[cp]][[region]] <- p1
  }
  
}

library(patchwork)

png(filename = paste0("./plotting/figures/micro_stacked_sources.png"),width=1500,height = 3000,res = 300)
print(p_list_all[[1]][["PFC"]]+NoLegend()+
        p_list_all[[1]][["INS"]]+NoLegend()+
        p_list_all[[1]][["VST"]]+NoLegend()+ plot_layout(ncol = 1)
)
dev.off()

