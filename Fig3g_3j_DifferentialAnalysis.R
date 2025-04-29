rm(list = ls())
library(plotly)
library(tidyr)
library(ggrepel)
library(gtools)
library(data.table)
library(gplots)
library(useful)
library(locfit)
library(Matrix)
library(dplyr)
library(ggplot2)
library(Seurat)
library(Rmagic)
library(readr)
library(phateR)
library(viridis)
library(magrittr)
library(harmony)
library(RColorBrewer)
library(tidyverse)
library(future)
library(GSVA)
library(GSEABase)
library(parallel)
library(sigPathway)
library(limma)
library(gskb)
library(Biobase)
library(genefilter)
library(GSVAdata)
data(c2BroadSets)
library(heatmap3)
plan("multicore", workers = (availableCores()-1))
options(future.globals.maxSize = 3000 * 1024^2)


#Set directory of dataset to analyse
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/Mapping_WTonly/Non_integrated_2")

project_name <- "RETINA_CD31_ECs_Integrated_WTonly"

Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))


#### For DGE between conditions globally

##Analysis global
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/DifferentialExpression/Non_integrated/Genes/OIRvsNORM") 

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Dataset %in% c("NORM_WT","OIR_WT")))

Seurat_object_Subset <- subset(Seurat_object, cells = Subseted_cells)


###Subset low cell count clusters

Seurat_object_Subset[["Cell_SubType_Condition"]] <- paste(Seurat_object_Subset@meta.data$Cell_SubType,Seurat_object_Subset@meta.data$Condition,sep="_")

Seurat_object_Subset <- SetIdent(Seurat_object_Subset, value = "Cell_SubType_Condition")

Seurat_object_Subset <- subset(Seurat_object_Subset, ident = c("Tuft_NORM_P14_WT","Proliferative ECs_NORM_P17_WT"), invert=T)

#Seurat_object_Subset <- SetIdent(Seurat_object_Subset, value = "Cell_SubType")


##Subset based on Age

Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Dataset %in% c("OIR_WT")))

Seurat_object_Subset_specific <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)



##Find DEGs in whole ECs

cells_in_this_cluster <- SetIdent(Seurat_object_Subset, value = "Condition")

  
#Perform DGE analysis using one of the model above

DGE_test_P14 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P14_WT", ident.2 = "NORM_P14_WT",logfc.threshold = 0.5, min.pct=0.2)
DGE_test_P14 <- subset(DGE_test_P14, p_val < 0.05)
write.table(DGE_test_P14, paste("wholeEC_OIRvsNORM_P14_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


top.1 <- rownames(head(DGE_test_P14[order(DGE_test_P14$avg_logFC),], 5))
top.2 <- rownames(head(DGE_test_P14[order(-DGE_test_P14$avg_logFC),], 5))

top10 <- c(top.1,top.2)  

dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("wholeEC_top10_DE_genes.Wilcox.OIRvsNORM_P14.png", sep=""), dpi=300, width = 10, height=5)


DGE_test_P17 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P17_WT", ident.2 = "NORM_P17_WT",logfc.threshold = 0.5, min.pct=0.2)
DGE_test_P17 <- subset(DGE_test_P17, p_val < 0.05)
write.table(DGE_test_P17, paste("wholeEC_OIRvsNORM_P17_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


top.1 <- rownames(head(DGE_test_P17[order(DGE_test_P17$avg_logFC),], 5))
top.2 <- rownames(head(DGE_test_P17[order(-DGE_test_P17$avg_logFC),], 5))
top10 <- c(top.1,top.2)  
dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))


ggsave(dotplot, filename=paste("wholeEC_top10_DE_genes.Wilcox.OIRvsNORM_P17.png", sep=""), dpi=300, width = 10, height=5)


#Time point confounded

cells_in_this_cluster <- SetIdent(Seurat_object_Subset, value = "Treatment")

DGE_test <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR", ident.2 = "NORM",logfc.threshold = -Inf, min.pct=-Inf)
DGE_test <- subset(DGE_test, p_val < 0.05)
write.table(DGE_test, paste("wholeEC_OIRvsNORM_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


top.1 <- rownames(head(DGE_test[order(DGE_test$avg_logFC),], 5))
top.2 <- rownames(head(DGE_test[order(-DGE_test$avg_logFC),], 5))

top10 <- c(top.1,top.2)  

dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("wholeEC_top10_DE_genes.Wilcox.OIRvsNORM.png", sep=""), dpi=300, width = 10, height=5)

#Make volcano plot


theme_Publication <- function(base_size=14, base_family="helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

volcanoplot <- ggplot(DGE_test,aes(avg_logFC,-log(p_val))) +
      geom_point(size=2) +
      geom_point(data = subset(DGE_test,(avg_logFC > 0.5) & (p_val < 0.05) & (p_val_adj < 0.05)),col="red") +
      geom_point(data = subset(DGE_test,(avg_logFC < -0.5 ) & (p_val < 0.05) & (p_val_adj < 0.05)),col="green") +
      geom_text_repel(
        data = subset(DGE_test, (avg_logFC > 0.5 | avg_logFC < -0.5) & (p_val < 0.05) & (p_val_adj < 0.05)),
        aes(label = rownames(subset(DGE_test, (avg_logFC > 0.5 | avg_logFC < -0.5) & (p_val < 0.05) & (p_val_adj < 0.05)))),
        size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) +
      theme_Publication() 
ggsave(volcanoplot, filename=paste("volcanoplot_wholeEC.Wilcox.OIRvsNORM.png", sep=""), dpi=300, width = 7, height=7)


##Find DEGs per subtype
#cell_type_list <- sort(unique(Seurat_object@ident))

#cell_type_list <- c("Tip ECs")                 

cell_type_list <- unique(Seurat_object@active.ident)

for(this_cluster in cell_type_list){

  #####this_cluster <- "Tip ECs"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetIdent(cells_in_this_cluster, value = "Condition")


        if(length(grep("OIR_P14_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("NORM_P14_WT", cells_in_this_cluster@meta.data$Condition))>2) {
  
        #Perform DGE analysis using one of the model above

        DGE_test_P14 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P14_WT", ident.2 = "NORM_P14_WT",logfc.threshold = 0.5, min.pct=0.2)
        DGE_test_P14 <- subset(DGE_test_P14, p_val < 0.05)
        write.table(DGE_test_P14, paste(this_cluster,"OIRvsNORM_P14_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


        top.1 <- rownames(head(DGE_test_P14[order(DGE_test_P14$avg_logFC),], 5))
        top.2 <- rownames(head(DGE_test_P14[order(-DGE_test_P14$avg_logFC),], 5))

        top10 <- c(top.1,top.2)  

        dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))

        ggsave(dotplot, filename=paste("cluster_",this_cluster,"_top10_DE_genes.Wilcox.OIRvsNORM_P14.png", sep=""), dpi=300, width = 10, height=5)

        }

        if(length(grep("OIR_P17_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("NORM_P17_WT", cells_in_this_cluster@meta.data$Condition))>2) {

        DGE_test_P17 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P17_WT", ident.2 = "NORM_P17_WT",logfc.threshold = 0.5, min.pct=0.2)
        DGE_test_P17 <- subset(DGE_test_P17, p_val < 0.05)
        write.table(DGE_test_P17, paste(this_cluster,"OIRvsNORM_P17_WT_DGE_test_p0.05.txt", sep="_"), sep="\t", row.names=TRUE, col.names=TRUE)


        top.1 <- rownames(head(DGE_test_P17[order(DGE_test_P17$avg_logFC),], 5))
        top.2 <- rownames(head(DGE_test_P17[order(-DGE_test_P17$avg_logFC),], 5))

        top10 <- c(top.1,top.2)  

        dotplot <- DotPlot(
        cells_in_this_cluster,
        assay = NULL,
        features=unique(top10),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA) + theme(axis.text.x = element_text(angle = 90))


        ggsave(dotplot, filename=paste("cluster_",this_cluster,"_top10_DE_genes.Wilcox.OIRvsNORM_P17.png", sep=""), dpi=300, width = 10, height=5)

        }




    ## Make venn diag for all differentially expressed genes containing testing results

        library(VennDiagram)
        flog.threshold(ERROR)

        if(length(grep("OIR_P14_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("NORM_P14_WT", cells_in_this_cluster@meta.data$Condition))>2
          & length(grep("OIR_P17_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("NORM_P17_WT", cells_in_this_cluster@meta.data$Condition))>2) {

        venn.diagram(x = list(P14=rownames(DGE_test_P14), P17=rownames(DGE_test_P17)), filename= paste("cluster_",this_cluster,"_VennDiagram_OIRvsNORM.png"), 
        height = 2000, width = 2000, resolution =500, 
        imagetype = "png", units = "px", compression ="lzw", 
        na = "stop", main = NULL, sub = NULL, main.pos= c(0.5, 1.05), main.fontface = "plain",
        main.fontfamily = "serif", main.col = "black",
        main.cex = 2, main.just = c(0.5, 1), sub.pos = c(0.5,
        1.05), sub.fontface = "plain", sub.fontfamily ="serif", sub.col = "blue", sub.cex = 1)
        }
}


###Check expression of specific genes

png(paste("FeaturePlot_ANGPT2_ESM1.png", sep=""), res=300, width = 4000, height=1000)
FeaturePlot(
  Seurat_object_Subset,
  c("ESM1", "ANGPT2"),
  dims = c(1, 2),
  cells = NULL,
  cols = c("red", "green"),
  pt.size = NULL,
  order = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  shape.by = NULL,
  slot = "data",
  blend = TRUE,
  blend.threshold = 0.5,
  label = FALSE,
  label.size = 4,
  repel = FALSE,
  ncol = NULL,
  coord.fixed = FALSE,
  by.col = TRUE,
  sort.cell = FALSE,
  combine = TRUE
)
dev.off()



png(paste("DotPlot_REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION_OF_SATURATED_FATTY_ACIDS.bycondition.png", sep=""), res=300, width = 4500, height=2000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("ACADL",
"ACADM",
"ACADS",
"ACADVL",
"ACSM3",
"ACSM6",
"ECHS1",
"HADH",
"HADHA",
"HADHB",
"MECR", "DECR1", "ECI1", "CPT1A", "SIRT3"),
  cols = c("darkblue", "orange", "purple", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = "Condition",
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)
dev.off()


png(paste("DotPlot_PARD.bycondition.png", sep=""), res=300, width = 2500, height=2000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = rownames(Seurat_object_Subset@assays$RNA[grep("CYP", rownames(Seurat_object_Subset@assays$RNA)),])
,
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_SubType_Condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)
dev.off()

dotplot <- DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c(rownames(Seurat_object_Subset@assays$RNA[grep("GPR", rownames(Seurat_object_Subset@assays$RNA)),]), rownames(Seurat_object_Subset@assays$RNA[grep("CYP", rownames(Seurat_object_Subset@assays$RNA)),])),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_SubType_Condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = 0,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
ggsave(dotplot, filename=paste("PaulPark/DotPlot_GPR_CYP.bycondition.pdf", sep=""),width = 37, height=6)

dotplot <- DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = "SOX9" ,
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_SubType_Condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = 0,
  scale.max = NA
)
ggsave(dotplot, filename=paste("Dubrac/DotPlot_SOX9.bycondition.pdf", sep=""),width = 7, height=6)


dotplot <- DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = rownames(Seurat_object)[grep("SLC25A20", rownames(Seurat_object))],
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_SubType_Condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = 0,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
ggsave(dotplot, filename=paste("Ketone/DotPlot_SLC25A20.bycondition.pdf", sep=""),width = 8, height=6)


dotplot <- DotPlot(
  Seurat_object_Subset_specific,
  assay = NULL,
  features = c("HMGCS2","HMGCL", "BDH1", "SLC16A6", "SLC16A1", "OXCT1", "ACAT1"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_SubType_Condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = 0,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
ggsave(dotplot, filename=paste("DotPlot_Ketogenesis_Ketolysis.Cell_SubType_Condition_P14.eps", sep=""),width = 8, height=8)

dotplot <- DotPlot(
  Seurat_object_Subset_specific,
  assay = NULL,
  features = c("ACAT1", "HMGCS2","HMGCL", "BDH1"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_SubType_Condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = 0,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
ggsave(dotplot, filename=paste("Ketone/DotPlot_Ketogenesis_Ketogenesis.Cell_SubType_Condition_OIR.P14_P17.eps", sep=""),width = 8, height=6)


dotplot <- DotPlot(
  Seurat_object_Subset_specific,
  assay = NULL,
  features = c("HMGCS2","HMGCL", "BDH1", "SLC16A6"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_SubType",
  split.by = NULL,
  scale.by = "radius",
  scale.min = 0,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
ggsave(dotplot, filename=paste("Ketone/DotPlot_Ketogenesis_only.Cell_SubType_OIR.P14_P17.eps", sep=""),width =6, height=6)


rownames(Seurat_object_Subset@assays$RNA[grep("PARD", rownames(Seurat_object_Subset@assays$RNA)),])

png(paste("DotPlot_tuft_markers.bycondition.png", sep=""), res=300, width = 3000, height=2000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("AQP1",
"PLVAP",
"HADHA", "HBEGF","TFPI"),
  cols = c("darkblue", "orange", "purple", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = "Condition",
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)
dev.off()


png(paste("DotPlot_Tipcells_glycolysis_markers_OIRvsNORM_P14_P17.png", sep=""), res=300, width = 3000, height=2000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("PKM",
"ALDOA",
"ENO1", "TPI1","PFKP", "LDHA", "CHST1", "CXCR4", "MIF", "PLOD1"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()

png(paste("DotPlot_Tipcells_glyco_ppp_hexo_markers_OIRvsNORM_P14_P17.png", sep=""), res=300, width = 3000, height=2000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = c("FBP1", "EXT1", "PFKP", "IDH1", "EGLN3", "SLC16A3", "PGK1", "HK1"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()

###Get genes from gene list

geneset1 <- importGeneSets("~/projects/def-jsjoyal/jhowa105/projects/Mike/MALLETTE_SENESCENCE_UP.gmx", verbose = TRUE)

geneset1 <- GeneSet(geneset1[[1]]$probes, geneIdType=SymbolIdentifier(), setName="MALLETTE_SENESCENCE_UP")

geneset2 <- importGeneSets("~/projects/def-jsjoyal/jhowa105/projects/Mike/GLOBAL_SENESCENCE_LITERATURE_CURATED_2020.gmx", verbose = TRUE)

geneset2 <- GeneSet(geneset2[[1]]$probes, geneIdType=SymbolIdentifier(), setName="GLOBAL_SENESCENCE_LITERATURE_CURATED")

geneset3 <- importGeneSets("~/projects/def-jsjoyal/jhowa105/projects/Mike/SAWCHYN_UNBIASED_UP_v4.gmx", verbose = TRUE)

geneset3 <- GeneSet(geneset3[[1]]$probes, geneIdType=SymbolIdentifier(), setName="SAWCHYN_UNBIASED_UP")

geneset4 <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/Pathways/GO_HEXOSE_CATABOLIC_PROCESS.gmt", geneIdType=SymbolIdentifier())

gmt.MsigDB.h <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MsigDB/V7.1/h.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier())

gmt.MsigDB.c2 <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MsigDB/V7.1/c2.cp.v7.1.symbols.gmt", geneIdType=SymbolIdentifier())

gmt.MsigDB.c5 <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MsigDB/V7.1/c5.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier())

gsc_final <- GeneSetCollection(c(gmt.MsigDB.h, gmt.MsigDB.c2,gmt.MsigDB.c5,geneset1,geneset2,geneset3,geneset4))

###

 geneset
 [1] "REACTOME_KETONE_BODY_METABOLISM"                   
 [2] "REACTOME_SYNTHESIS_OF_KETONE_BODIES"               
 [3] "GO_CELLULAR_KETONE_METABOLIC_PROCESS"              
 [4] "GO_RESPONSE_TO_KETONE"                             
 [5] "GO_REGULATION_OF_CELLULAR_KETONE_METABOLIC_PROCESS"
 [6] "GO_KETONE_BIOSYNTHETIC_PROCESS"                    
 [7] "GO_CELLULAR_RESPONSE_TO_KETONE"                    
 [8] "GO_REGULATION_OF_KETONE_BIOSYNTHETIC_PROCESS"      
 [9] "GO_KETONE_CATABOLIC_PROCESS"                       
[10] "GO_KETONE_BODY_BIOSYNTHETIC_PROCESS"               
[11] "GO_KETONE_BODY_METABOLIC_PROCESS"       



geneset <- names(gsc_final)[unique(c(
	grep("GO_CELLULAR_KETONE_METABOLIC_PROCESS", names(gsc_final))
	))]

unique(geneIds(gsc_final[[geneset]]))

dotplot <- DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = unique(geneIds(gsc_final[[geneset]])),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_GO_CELLULAR_KETONE_METABOLIC_PROCESS.byTreatment.pdf", sep=""),width = 44, height=8)



Pathway <- c("OXCT1", "SLC16A6", "BDH1", "BDH2", "HMGCL", "HMGCS2", "HMGCLL1")
Pathway <- c("OXCT1", "SLC16A1", "SLC16A6", "BDH1", "BDH2", "HMGCL", "HMGCS2", "HMGCLL1")


dotplot <- DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = Pathway,
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_Ketones_Genes_2.byTreatment.pdf", sep=""), width = 10, height=8)


###FInd markers using gene list

geneset <- names(gsc_final)[unique(c(
	grep("HALLMARK_GLYCOLYSIS", names(gsc_final)),
	grep("KEGG_PENTOSE_PHOSPHATE_PATHWAY", names(gsc_final)),
	grep("REACTOME_PENTOSE_PHOSPHATE_PATHWAY", names(gsc_final)),
	grep("GO_PENTOSE_PHOSPHATE_SHUNT_NON_OXIDATIVE_BRANCH", names(gsc_final)),
	grep("GO_HEXOSE_CATABOLIC_PROCESS", names(gsc_final))
	))]


genelist <- unique(c(geneIds(gsc_final[["HALLMARK_GLYCOLYSIS"]]),
		 geneIds(gsc_final[["KEGG_PENTOSE_PHOSPHATE_PATHWAY"]]),
		 geneIds(gsc_final[["REACTOME_PENTOSE_PHOSPHATE_PATHWAY"]]),
		 geneIds(gsc_final[["GO_PENTOSE_PHOSPHATE_SHUNT_NON_OXIDATIVE_BRANCH"]]),
		 geneIds(gsc_final[["GO_HEXOSE_CATABOLIC_PROCESS"]])
	   ))


genelist <- intersect(genelist,rownames(Seurat_object_Subset))


###Subset condition of interest

cond4subset <- c("OIR_WT")

Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Dataset %in% c(cond4subset)))

Seurat_object_Subset_Condition <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)

DGE_test <- FindMarkers(Seurat_object_Subset_Condition, ident.1 = "Tuft", ident.2 = "Tip ECs", logfc.threshold = -Inf, min.pct=0.01, features = NULL, only.pos=FALSE)

DGE_test <- subset(DGE_test, p_val < 0.05)
write.table(DGE_test, "TuftVSTip_OIR_P17_DEGs_p05.txt", sep="\t")

dotplot <- DotPlot(
  Seurat_object_Subset_Condition,
  assay = NULL,
  features = rownames(DGE_test),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_Tip_P14_Markers_GLYCO_PPP_HEXOSE.",cond4subset,".png", sep=""), dpi=150, width = 20, height=6)


dotplot <- DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  features = rownames(DGE_test),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_Tip_P14_Markers_GLYCO_PPP_HEXOSE.byTreatment.png", sep=""), dpi=150, width = 20, height=8)

dotplot <- DotPlot(
  Seurat_object_Subset_Condition,
  assay = NULL,
  features = c("SIRT3", "HADHA", "HADH", "HDLBP", "ACADVL", "SLC25A20", "CPT2", "FAR2", "NAMPT", "HMOX1", "NT5C3", "TEAD2", "NKD1", "INSR"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_Sirt3_FAO.",cond4subset,".png", sep=""), dpi=150, width = 8, height=4)

Seurat_object_Subset_Condition <- ReorderIdent(Seurat_object_Subset_Condition, 
var = c("HADHA", "ACADVL", "NAMPT", "PKM", "ALDOA", "TPI1", "GPI1"), reverse = FALSE)

dotplot <- DotPlot(
  Seurat_object_Subset_Condition,
  assay = NULL,
  features = c("AQP1", "HADHA", "ACADVL", "SLC25A20", "ESM1", "PKM", "ALDOA", "TPI1", "GPI1", "PFKFB3", "EGLN3"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_HADHA_AQP1_ESM1_EGLN3.",cond4subset,".eps", sep=""), dpi=300, width = 8, height=4)

dotplot <- DotPlot(
  Seurat_object_Subset_Condition,
  assay = NULL,
  features = c("PKM", "ALDOA", "TPI1", "GPI1", "HADHA", "ACADVL", "NAMPT", "NT5C3"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_TipvsTuft_metabolism.",cond4subset,".eps", sep=""), dpi=150, width = 8, height=4)


dotplot <- DotPlot(
  Seurat_object_Subset_Condition,
  assay = NULL,
  features = c("NAMPT"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_NAMPT.",cond4subset,".eps", sep=""), dpi=150, width = 8, height=4)


dotplot <- DotPlot(
  Seurat_object_Subset_Condition,
  assay = NULL,
  features = c("AQP1","ESM1"),
  cols = c("darkblue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "Cell_SubType_Condition",
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))

ggsave(dotplot, filename=paste("DotPlot_AQP1_ESM1.",cond4subset,".eps", sep=""), dpi=150, width = 8, height=4)


###Make volcano plot with enhanceVolcano

library(EnhancedVolcano)


keyvals.colour <- ifelse(
    DGE_test$avg_logFC < -0.5 & DGE_test$p_val < 0.05, 'royalblue',
      ifelse(DGE_test$avg_logFC > 0.5 & DGE_test$p_val < 0.05, 'red',
        'grey'))
  keyvals.colour[is.na(keyvals.colour)] <- 'grey'
  names(keyvals.colour)[keyvals.colour == 'red'] <- 'high'
  names(keyvals.colour)[keyvals.colour == 'grey'] <- 'mid'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'low'


genetohi <- unique(c('AQP1', "TRP53I11", "APLN", "APLNR", "ANGPT2", "TPI1", "PDGFB", "GPI1",
      "COX6B1", "VWF", "ALDOA", "ATPIF1" ,"PKM", "COL18A1", "SLC2A1", "INSR", "ELOVL7", "ACADVL", "NT5C3", "NAMPT", "CACNA1A",
      "PIK3R3", "SGMS1", "TEAD2", "SREBF2", "BPGM", "CFH","HTRA1","CLU", "HMOX1", "MYLK", "ID1", "NLK", "CASP2", "IFITM3",
      "CTHRC1", "SERPINE2", "SLC6A6", "HEY1", "HEY2", "PAPSS1", "SMAD6", "SIX3", "SIX3OS1"))

genetohi <- unique(c('AQP1', "TRP53I11", "APLN", "APLNR", "ANGPT2", "TPI1", "GPI1",
      "VWF", "ALDOA","PKM", "COL18A1", "SLC2A1", "INSR", "ELOVL7", "ACADVL", "NT5C3", "NAMPT",
      "SGMS1", "TEAD2", "SREBF2", "CFH","HTRA1","HMOX1", "CASP2",
      "CTHRC1", "SERPINE2", "SIX3", "SIX3OS1", "HADHA"))


volcanoplot <- EnhancedVolcano(DGE_test,
    lab = rownames(DGE_test),
    selectLab = genetohi,
    x = 'avg_logFC',
    y = 'p_val',
    xlim = c(-5, 5),
    title = 'Tip versus Tuft',
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 1.0,
    labSize = 4.0,
    col=c('grey', 'pink', 'pink', 'red3'),
    colAlpha = 1,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    colCustom = keyvals.colour,
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey30'
    )

ggsave(volcanoplot, filename=paste("volcanoplot_TipVSTuft.",cond4subset,".eps", sep=""), width = 8, height=8)




q("no")
