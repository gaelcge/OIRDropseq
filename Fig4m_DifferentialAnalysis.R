
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
library(limma)
library(future)
library(dendextend)
library(VennDiagram)
flog.threshold(ERROR)
plan("multiprocess", workers = (availableCores()-1))
options(future.globals.maxSize = 3000 * 1024^3)


#Set directory of dataset to analyse
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/Mapping")

project_name <- "RETINA_CD31_ECs_Integrated"

Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

Seurat_object[["Cell_Type_Condition"]] <- paste(Seurat_object@meta.data$Cell_Type,Seurat_object@meta.data$Condition, sep="_")
Seurat_object[["Cell_Type_Dataset"]] <- paste(Seurat_object@meta.data$Cell_Type,Seurat_object@meta.data$Dataset, sep="_")


###Add GO function to Seurat object


GO_score=fread("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/GSVA/May2020/gsva.exprs_data.msigdb_May2020.csv", sep=",",data.table=FALSE)
rownames(GO_score) <- GO_score$V1
GO_score$V1 <- NULL

Seurat_object[["GO"]] <- CreateAssayObject(counts = GO_score)


Seurat_object <- NormalizeData(Seurat_object, assay = "GO", normalization.method = "CLR", margin = 2)

Seurat_object <- ScaleData(Seurat_object, assay = "GO", verbose = FALSE, model.use = "negbinom")


#Seurat_object <- SetAllIdent(Seurat_object, id="cell_subtype")


###Modify ident and cell type to use for analysis

Seurat_object <- RenameIdents(object = Seurat_object, 'Rods' = 'Photoreceptors', 'Cones' = 'Photoreceptors')
Seurat_object <- RenameIdents(object = Seurat_object, 'Amacrine cells' = 'Neuronal cells', 'Bipolar cells' = 'Neuronal cells')
Seurat_object <- RenameIdents(object = Seurat_object, 'Muller Glial cells' = 'Glial cells', 'Astrocytes' = 'Glial cells')


Seurat_object <-SubsetData(
  Seurat_object,
  assay = NULL,
  cells = NULL,
  subset.name = NULL,
  ident.use = NULL,
  ident.remove = c("RBCs", "Lens cells"))

Seurat_object[["Cell_Type"]] <- Idents(Seurat_object)

Seurat_object[["Cell_Type_Dataset"]] <- paste(Seurat_object@meta.data$Cell_Type,Seurat_object@meta.data$Dataset,sep="_")

##Filter on crystal genes


CRYA.genes <- grep(pattern = "^CRYA", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
CRYB.genes <- grep(pattern = "^CRYB", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
CRYG.genes <- grep(pattern = "^CRYG", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
crystal.genes <- c(CRYA.genes, CRYB.genes)
percent.crystal <- Matrix::colSums(Seurat_object@assays$RNA[crystal.genes, ])/Matrix::colSums(Seurat_object@assays$RNA)*100
Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.crystal, col.name = "percent.crystal")

Seurat_object <- subset(Seurat_object, subset = percent.crystal < 0.25)




#### For DGE between conditions globally

##Analysis global
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/DifferentialExpression/Functions/OIRKOvsOIRWT") 

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Dataset %in% c("OIR_WT", "OIR_Sirt3KO")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells = Subseted_cells)

Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Age %in% c("P17")))

Seurat_object_Subset <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)


png(filename="umap.sct_Condition_Dataset.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object_Subset, reduction = "umap", label=FALSE, group.by = "Dataset", cols =c("green", "red"))
dev.off()

setEPS(bg = "white", family = "Times", width=14, height=7)
postscript("umap.sct_Condition_split_Dataset.eps")
DimPlot(Seurat_object_Subset, reduction = "umap", label=FALSE, split.by = "Dataset", group.by = "Cell_Type")
dev.off()


###Make distribution plots


#cell_type_list <- sort(unique(Seurat_object@ident))

#cell_type_list <- c("Muller Glial cells")                 

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


        if(length(grep("OIR_P14_Sirt3KO", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("OIR_P14_WT", cells_in_this_cluster@meta.data$Condition))>2) {
  
        #Perform DGE analysis using one of the model above

        DGE_test_P14 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P14_Sirt3KO", ident.2 = "OIR_P14_WT",logfc.threshold = 0, min.pct=0.2, assay="GO")
        DGE_test_P14 <- subset(DGE_test_P14, p_val < 0.05)
        write.table(DGE_test_P14,paste("cluster_",this_cluster,"_DE_function.Wilcox.OIRKOvsWT_P14.txt", sep=""),sep="\t")

        top.1 <- rownames(head(DGE_test_P14[order(DGE_test_P14$avg_logFC),], 5))
        top.2 <- rownames(head(DGE_test_P14[order(-DGE_test_P14$avg_logFC),], 5))

        top10 <- c(top.1,top.2)  



        ridgeplot <- RidgePlot(
        cells_in_this_cluster,
        features=unique(top10),
        cols = viridis(4, dir=-1),
        idents = NULL,
        sort = TRUE,
        assay = "GO",
        group.by = NULL,
        y.max = NULL,
        same.y.lims = FALSE,
        log = FALSE,
        ncol = NULL,
        slot = "data",
        combine = TRUE
        )

        ggsave(ridgeplot, filename=paste("cluster_",this_cluster,"_top10_DE_functions.Wilcox.OIR_Sirt3KOvsWT_P14.png", sep=""), dpi=150, width = 15, height=7)

        
        Vlnplot <- VlnPlot(
        cells_in_this_cluster,
        features=unique(top10),
        cols = viridis(4, dir=-1),
        pt.size = 1,
        idents = NULL,
        sort = TRUE,
        assay = "GO",
        group.by = NULL,
        split.by = NULL,
        adjust = 1,
        y.max = NULL,
        same.y.lims = FALSE,
        log = FALSE,
        ncol = NULL,
        slot = "data",
        combine = TRUE
        )

        ggsave(Vlnplot, filename=paste("cluster_",this_cluster,"_top10_DE_functions.Wilcox.OIR_Sirt3KOvsWT_P14.Vlnplot.png", sep=""), dpi=150, width = 25, height=15)



        }

        if(length(grep("OIR_P17_Sirt3KO", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("OIR_P17_WT", cells_in_this_cluster@meta.data$Condition))>2) {

        DGE_test_P17 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_P17_Sirt3KO", ident.2 = "OIR_P17_WT",logfc.threshold = 0, min.pct=0.5, assay="GO")
        DGE_test_P17 <- subset(DGE_test_P17, p_val < 0.05)
        #DGE_test_P17_test <- subset(DGE_test_P17, myAUC < 0.4 | myAUC > 0.6)
        write.table(DGE_test_P17,paste("cluster_",this_cluster,"_DE_function.Wilcox.OIRKOvsWT_P17.txt", sep=""),sep="\t")

        top.1 <- rownames(head(DGE_test_P17[order(DGE_test_P17$avg_logFC),], 5))
        top.2 <- rownames(head(DGE_test_P17[order(-DGE_test_P17$avg_logFC),], 5))

        top10 <- c(top.1,top.2)  

        ridgeplot <- RidgePlot(
        cells_in_this_cluster,
        features=unique(top10),
        cols = viridis(4, dir=-1),
        idents = NULL,
        sort = TRUE,
        assay = "GO",
        group.by = NULL,
        y.max = NULL,
        same.y.lims = FALSE,
        log = FALSE,
        ncol = NULL,
        slot = "data",
        combine = TRUE
        )

        ggsave(ridgeplot, filename=paste("cluster_",this_cluster,"_top10_DE_functions.Wilcox.OIR_Sirt3KOvsWT_P17.png", sep=""), dpi=150, width = 15, height=7)


        Vlnplot <- VlnPlot(
        cells_in_this_cluster,
        features=unique(top10),
        cols = viridis(4, dir=-1),
        pt.size = 1,
        idents = NULL,
        sort = TRUE,
        assay = "GO",
        group.by = NULL,
        split.by = NULL,
        adjust = 1,
        y.max = NULL,
        same.y.lims = FALSE,
        log = FALSE,
        ncol = NULL,
        slot = "data",
        combine = TRUE
        )

        ggsave(Vlnplot, filename=paste("cluster_",this_cluster,"_top10_DE_functions.Wilcox.OIR_Sirt3KOvsWT_P17.Vlnplot.png", sep=""), dpi=150, width = 25, height=15)

        }




    ## Write table for all differentially expressed genes containing testing results



        if(length(grep("OIR_P14_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("OIR_P14_Sirt3KO", cells_in_this_cluster@meta.data$Condition))>2
          & length(grep("OIR_P17_WT", cells_in_this_cluster@meta.data$Condition))>2 & length(grep("OIR_P17_Sirt3KO", cells_in_this_cluster@meta.data$Condition))>2) {

        venn.diagram(x = list(P14=rownames(DGE_test_P14), P17=rownames(DGE_test_P17)), filename= paste("cluster_",this_cluster,"_VennDiagram_OIR_Sirt3KOvsWT.png"), 
        height = 2000, width = 2000, resolution =500, 
        imagetype = "png", units = "px", compression ="lzw", 
        na = "stop", main = NULL, sub = NULL, main.pos= c(0.5, 1.05), main.fontface = "plain",
        main.fontfamily = "serif", main.col = "black",
        main.cex = 2, main.just = c(0.5, 1), sub.pos = c(0.5,
        1.05), sub.fontface = "plain", sub.fontfamily ="serif", sub.col = "blue", sub.cex = 1)
        }
}


#########
########



###Make heatmap of KO vs WT for specific functions acrross cell types

rowf <- unique(c(grep("HALLMARK-MYC-", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-FOXO", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-REGULATION-OF-LOCALIZATION-OF-FOXO", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-REGULATION-OF-FOXO", rownames(Seurat_object_Subset[["GO"]])),
          grep("GO-HYPOXIA-INDUCIBLE-FACTOR-1ALPHA-SIGNALING-PATHWAY", rownames(Seurat_object_Subset[["GO"]]))
          ))

rowf <- unique(c(grep("HALLMARK-MYC-TARGETS-V1", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-MYC-TARGETS-V2", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-FOXO-MEDIATED-TRANSCRIPTION$", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-REGULATION-OF-LOCALIZATION-OF-FOXO-TRANSCRIPTION-FACTORS", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-REGULATION-OF-FOXO-TRANSCRIPTIONAL-ACTIVITY-BY-ACETYLATION", rownames(Seurat_object_Subset[["GO"]])),
          grep("GO-HYPOXIA-INDUCIBLE-FACTOR-1ALPHA-SIGNALING-PATHWAY", rownames(Seurat_object_Subset[["GO"]]))
          ))

rowf <- unique(c(grep("REACTOME-MITOCHONDRIAL-FATTY-ACID-BETA-OXIDATION-OF-SATURATED-FATTY-ACIDS", rownames(Seurat_object_Subset[["GO"]])),
      grep("REACTOME-REGULATION-OF-FOXO-TRANSCRIPTIONAL-ACTIVITY-BY-ACETYLATION", rownames(Seurat_object_Subset[["GO"]])),
          grep("GO-HYPOXIA-INDUCIBLE-FACTOR-1ALPHA-SIGNALING-PATHWAY", rownames(Seurat_object_Subset[["GO"]])),
          grep("HALLMARK-GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]]))
          ))

rowf <- unique(c(grep("REACTOME-MITOCHONDRIAL-FATTY-ACID-BETA-OXIDATION-OF-SATURATED-FATTY-ACIDS", rownames(Seurat_object_Subset[["GO"]])),
       grep("HALLMARK-GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]]))
          ))

rowf <- unique(c(grep("REACTOME-REGULATION-OF-FOXO-TRANSCRIPTIONAL-ACTIVITY-BY-ACETYLATION", rownames(Seurat_object_Subset[["GO"]])),
       grep("MYC-TARGETS-V1", rownames(Seurat_object_Subset[["GO"]])),
       grep("GO-REGULATION-OF-FATTY-ACID-BETA-OXIDATION", rownames(Seurat_object_Subset[["GO"]])),
       grep("KEGG-GLYCOLYSIS-GLUCONEOGENESIS", rownames(Seurat_object_Subset[["GO"]]))
          ))


rowf <- rev(unique(c(
      grep("GO-AEROBIC-RESPIRATION", rownames(Seurat_object_Subset[["GO"]])),
      grep("REACTOME-MITOCHONDRIAL-FATTY-ACID-BETA-OXIDATION-OF-SATURATED-FATTY-ACIDS", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-HYPOXIA", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-SPROUTING-ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]]))
          )))

rowf <- rev(unique(c(grep("-TCF", rownames(Seurat_object_Subset[["GO"]])),
      grep("FOXO", rownames(Seurat_object_Subset[["GO"]])),
      grep("PPARG", rownames(Seurat_object_Subset[["GO"]])),
      grep("HYPOXIA", rownames(Seurat_object_Subset[["GO"]])),
      grep("ESTROGEN", rownames(Seurat_object_Subset[["GO"]])),
      grep("AUTOPHAGY", rownames(Seurat_object_Subset[["GO"]]))
          )))

rowf <- rev(unique(c(
      grep("REACTOME-OPSINS", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-REGULATION-OF-LONG-TERM-SYNAPTIC-DEPRESSION", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-PROTEIN-CHROMOPHORE-LINKAGE", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-AEROBIC-ELECTRON-TRANSPORT-CHAIN", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-PHOTOTRANSDUCTION", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-MAIN-AXON", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-AXON-INITIAL-SEGMENT", rownames(Seurat_object_Subset[["GO"]])),
      grep("REACTOME-MITOCHONDRIAL-FATTY-ACID-BETA-OXIDATION-OF-SATURATED-FATTY-ACIDS", rownames(Seurat_object_Subset[["GO"]]))
          )))

rowf <- rev(unique(c(
      grep("REACTOME-MITOCHONDRIAL-FATTY-ACID-BETA-OXIDATION", rownames(Seurat_object_Subset[["GO"]]))
          )))

rowf <- rev(unique(c(
      grep("HALLMARK-GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("REACTOME-GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("GLUCOSE", rownames(Seurat_object_Subset[["GO"]]))
          )))


rowf <- rev(unique(c(
      grep("HALLMARK-ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-REGULATION-OF-SPROUTING-ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-BLOOD-VESSEL-ENDOTHELIAL-CELL-PROLIFERATION-INVOLVED-IN-SPROUTING-ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-SPROUTING-ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-CELL-MIGRATION-INVOLVED-IN-SPROUTING-ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("REACTOME-MITOCHONDRIAL-FATTY-ACID-BETA-OXIDATION", rownames(Seurat_object_Subset[["GO"]]))
          )))


rowf <- rev(unique(c(
      grep("HALLMARK-ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-HYPOXIA", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-OXIDATIVE-PHOSPHORYLATION", rownames(Seurat_object_Subset[["GO"]]))
          )))

rowf <- rev(unique(c(
      grep("HALLMARK-ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-HYPOXIA", rownames(Seurat_object_Subset[["GO"]])),
      grep("REACTOME-MITOCHONDRIAL-FATTY-ACID-BETA-OXIDATION$", rownames(Seurat_object_Subset[["GO"]]))
          )))

selected_functions <- rownames(Seurat_object_Subset[["GO"]])[rowf]

dataframe.all <- data.frame()

dataframe.all <- dataframe.all[1:length(selected_functions),]

rownames(dataframe.all) <- selected_functions

Seurat_object_loop <- SetIdent(Seurat_object_Subset, value = "Age")

timepoint <- c("P14", "P17")

Seurat_object_loop <- SubsetData(Seurat_object_loop,
                                      ident.use=timepoint)

Seurat_object_loop <- SetIdent(Seurat_object_loop, value = "Cell_Type")  

###
loop_list <- unique(Seurat_object_loop@active.ident)

#loop_list <- loop_list[loop_list != "Tuft ECs"]

for(this_cluster in loop_list){

  #####this_cluster <- "Tip_ECs"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_loop,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetIdent(cells_in_this_cluster, value = "Dataset")

  DGE_test_P14 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR_Sirt3KO", ident.2 = "OIR_WT", logfc.threshold = -Inf, min.pct=-Inf, 
    features=selected_functions, assay = "GO")

  #DGE_test_P14 <- subset(DGE_test_P14, p_val < 0.05)

  DGE_test_P14_avg_logFC <- as.data.frame(DGE_test_P14$avg_logFC)

  rownames(DGE_test_P14_avg_logFC) <- rownames(DGE_test_P14)

  colnames(DGE_test_P14_avg_logFC) <- paste(this_cluster, "avg_logFC_OIR_KO_vs_OIR_WT", "P14_P17_Comb", sep="_")

  dataframe.all <- merge(dataframe.all,DGE_test_P14_avg_logFC, by=0, all=TRUE)

  rownames(dataframe.all) <- dataframe.all$Row.names
 
  dataframe.all$Row.names <- NULL
}

dataframe <- dataframe.all

##Select for specific subtypes
#dataframe <- dataframe[complete.cases(dataframe[ ,1]),]

#Remove row containing only NA
ind <- apply(dataframe, 1, function(x) all(is.na(x)))

dataframe <- dataframe[ !ind, ]

#### make matrix for heatmap

my_palette <- colorRampPalette(c("royalblue4", "royalblue", "white", "orangered", "red4"))

dataframe <- as.matrix(dataframe)


Colv  <- dataframe %>% t %>% dist(method = "euclidean") %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", k=2) %>% set("branches_lwd", 2) %>%
           ladderize #%>% rotate_DendSer

        #png("Dend.png")
        #plot(Colv)
        #dev.off()


Rowv  <- dataframe %>% dist %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", 2) %>% set("branches_lwd", 2) %>%
           ladderize


dataframe[is.na(dataframe)] <- 0
dataframe[is.nan(dataframe)] <- 0

dim(dataframe)

dataframe <- dataframe[selected_functions,]

setEPS(bg = "white", family = "Times", width=35, height=15)
postscript(paste("OIRSirt3KO_vs_OIRWT_","_P14_P17_Comb_","_FAbO_GLYCOLYSIS_Hypoxia_Angiogenesis_AllCells_heatmap.eps",sep=""))
heatmap.2(dataframe, col=my_palette(100), symbreak=TRUE, trace='none', cexRow=2, cexCol= 2,
          Rowv=Rowv,
          #Rowv=TRUE,
          Colv=Colv,
          #Colv=TRUE,
          srtCol=45,
          revC=FALSE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(2,5), 
          key.par=list(mar=c(5, 8, 8, 8)), 
          margins=c(45,100), key.ylab=NA, main="OIR_KO vs OIR_WT")
dev.off()





##Make heatmap plot using limma

cell_type_list <- unique(Seurat_object@active.ident)

for(this_cluster in cell_type_list){

  #####this_cluster <- "Tip ECs"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetIdent(cells_in_this_cluster, value = "Treatment")

  
        #Perform DGE analysis using one of the model above

        DGE_test <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR", ident.2 = "NORM",logfc.threshold = 0, min.pct=0.2, assay="GO")
        DGE_test <- subset(DGE_test, p_val < 0.05)

        top.1 <- rownames(head(DGE_test[order(DGE_test$avg_logFC),], 10))
        top.2 <- rownames(head(DGE_test[order(-DGE_test$avg_logFC),], 10))

        top10 <- c(top.1,top.2)  

        Selected.genes <- unique(top10)

        ##Make dataframe

        dataframe <- data.frame()

        dataframe <- dataframe[1:length(Selected.genes),]

        rownames(dataframe) <- Selected.genes

        for(this_cluster_2 in cell_type_list){

          #####this_cluster <- "Pericytes"

          ## Print status for which identity is being processed
          print(paste("Working on cluster #",this_cluster_2,sep=""))

          ## Subset Seurat object to only contain cells from this cluster
          cells_in_this_cluster <- SubsetData(Seurat_object_Subset,
                                              ident.use=this_cluster_2)

          ## Check whether there are cells in both groups, otherwise skip this cluster
                
          cells_in_this_cluster <- SetIdent(cells_in_this_cluster, value = "Treatment")
          
          #Perform DGE analysis using one of the model above
          
          average_cells_in_this_cluster <- AverageExpression(
          cells_in_this_cluster,
          assays = "GO",
          features = Selected.genes,
          return.seurat = FALSE,
          add.ident = NULL,
          slot = "scale.data",
          use.scale = FALSE,
          use.counts = FALSE,
          verbose = TRUE)

          average_cells_in_this_cluster <- as.data.frame(average_cells_in_this_cluster)
        
          design <- model.matrix(~0+colnames(average_cells_in_this_cluster))

          colnames(design) <- c("NORM","OIR")

          contrast_dir <- paste(colnames(design)[2], "vs", colnames(design)[1], sep="")

          contrast <- makeContrasts(OIR - NORM, levels = design)

          fit <- lmFit(average_cells_in_this_cluster, design)

          fit <- contrasts.fit(fit, contrast)

          sel_FC <- as.data.frame(fit$coefficients[Selected.genes,])

          colnames(sel_FC) <- paste(this_cluster_2)

          dataframe <- cbind(dataframe,sel_FC)
        }

        my_palette <- colorRampPalette(c("blue", "white", "red"))

        dataframe <- as.matrix(dataframe)


        Colv  <- dataframe %>% t %>% dist(method = "euclidean") %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", k=2) %>% set("branches_lwd", 2) %>%
           ladderize %>% rotate_DendSer

        #png("Dend.png")
        #plot(Colv)
        #dev.off()


        Rowv  <- dataframe %>% dist %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", k = 2) %>% set("branches_lwd", 2) %>%
           ladderize


        png(filename=paste("cluster_",this_cluster,"_top10_DE_functions.Wilcox.OIRvsNORM_heatmap.png", sep=""), width=2000, height=1500, res = 150)  
        heatmap.2(dataframe, col=my_palette, symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, 
          Rowv=Rowv,
          #Rowv=TRUE,
          Colv=Colv,
          #Colv=TRUE,
          srtCol=45,
          revC=TRUE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(1,6), 
          key.par=list(mar=c(3, 1, 1, 0)), 
          margins=c(10,35), key.ylab=NA, main="top10_DE_functions in OIR vs NORM retina")
        dev.off()

}



###Check expression of specific functions

##On p17 time point

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Age %in% c("P17")))

Seurat_object_P17 <- SubsetData(Seurat_object, cells = Subseted_cells)




functionOfInterest <- "GO-RECEPTOR-TRANSACTIVATION"

functions <- rownames(Seurat_object_Subset$GO)[grep(functionOfInterest, rownames(Seurat_object_Subset$GO))]

png(paste(functionOfInterest, "FeaturePlot.png", sep=""), res=100, width = 2000, height=2000)
FeaturePlot(
  Seurat_object_Subset,
  functions,
  dims = c(1, 2),
  cells = NULL,
  cols = c("lawngreen","yellow","red"),
  pt.size = 0.3,
  order = TRUE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  shape.by = NULL,
  slot = "data",
  blend = FALSE,
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

functionOfInterest <- "Sirt3"

functions <- rownames(Seurat_object_Subset$GO)[grep(functionOfInterest, rownames(Seurat_object_Subset$GO))]

pdf(paste(functionOfInterest,"Cell_Type_Dataset_RidgePlot.pdf", sep=""), width = 10, height=8)
RidgePlot(
        Seurat_object_Subset,
        features=functions,
        cols = viridis(18, dir=-1),
        idents = NULL,
        sort = TRUE,
        assay = "GO",
        group.by = "Cell_Type_Dataset",
        y.max = NULL,
        same.y.lims = FALSE,
        log = FALSE,
        ncol = NULL,
        slot = "scale.data",
        combine = TRUE
        )
dev.off()


png(paste(functionOfInterest,"RidgePlot.ECs.byconditions.png", sep=""), res=100, width = 3000, height=3000)
RidgePlot(
        Seurat_object_Subset,
        features=functions,
        cols = viridis(4, dir=-1),
        idents = "ECs",
        sort = TRUE,
        assay = "GO",
        group.by = "Condition",
        y.max = NULL,
        same.y.lims = FALSE,
        log = FALSE,
        ncol = NULL,
        slot = "scale.data",
        combine = TRUE
        )
dev.off()


functionOfInterest <- "BIOCARTA-ERBB4-PATHWAY"

functions <- rownames(Seurat_object_Subset$GO)[grep(functionOfInterest, rownames(Seurat_object_Subset$GO))]

setEPS()
postscript(paste("MullerGlia/",functionOfInterest,"VlnPlot.MullerGlia.byconditions.eps", sep=""), width = 10, height=10)
VlnPlot(
        Seurat_object_P17,
        features=functions,
        cols = viridis(3, dir=-1),
        idents = "Muller Glial cells",
        sort = TRUE,
        assay = "GO",
        group.by = "Condition",
        y.max = NULL,
        same.y.lims = FALSE,
        log = FALSE,
        ncol = NULL,
        slot = "data",
        combine = TRUE
        )
dev.off()

###Make heatmap for specific functions

functionOfInterest <- "SYNAPTIC"

Selected.genes <- rownames(Seurat_object_Subset$GO)[grep(functionOfInterest, rownames(Seurat_object_Subset$GO))]



dataframe.all <- data.frame()

dataframe.all <- dataframe.all[1:length(Selected.genes),]

rownames(dataframe.all) <- Selected.genes

##Make dataframe

timepoint_list <- rev(unique(Seurat_object_Subset@meta.data$Age))

for(this_timepoint in timepoint_list){

    print(paste("Working on data #",this_timepoint,sep=""))

    Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Age %in% this_timepoint))

    cells_in_this_timepoint <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)


    dataframe <- data.frame()

    dataframe <- dataframe[1:length(Selected.genes),]

    rownames(dataframe) <- Selected.genes


    cell_type_list <- unique(Seurat_object@active.ident)

for(this_cluster_2 in cell_type_list){

          #####this_cluster <- "Pericytes"

          ## Print status for which identity is being processed
          print(paste("Working on cluster #",this_cluster_2,sep=""))

          ## Subset Seurat object to only contain cells from this cluster
          cells_in_this_cluster <- SubsetData(cells_in_this_timepoint,
                                              ident.use=this_cluster_2)

          ## Check whether there are cells in both groups, otherwise skip this cluster
                
          cells_in_this_cluster <- SetIdent(cells_in_this_cluster, value = "Genotype")
          
          #Perform DGE analysis using one of the model above
          
          average_cells_in_this_cluster <- AverageExpression(
          cells_in_this_cluster,
          assays = "GO",
          features = Selected.genes,
          return.seurat = FALSE,
          add.ident = NULL,
          slot = "scale.data",
          use.scale = FALSE,
          use.counts = FALSE,
          verbose = TRUE)

          average_cells_in_this_cluster <- as.data.frame(average_cells_in_this_cluster)
        
          design <- model.matrix(~0+colnames(average_cells_in_this_cluster))

          colnames(design) <- c("Sirt3KO","WT")

          contrast_dir <- paste(colnames(design)[2], "vs", colnames(design)[1], sep="")

          contrast <- makeContrasts(Sirt3KO - WT, levels = design)

          fit <- lmFit(average_cells_in_this_cluster, design)

          fit <- contrasts.fit(fit, contrast)

          sel_FC <- as.data.frame(fit$coefficients[Selected.genes,])

          colnames(sel_FC) <- paste(this_cluster_2, this_timepoint)

          dataframe <- cbind(dataframe,sel_FC)
        }

    dataframe.all <- cbind(dataframe.all,dataframe)
    dataframe.all[,this_timepoint] <- NA

    my_palette <- colorRampPalette(c("blue", "white", "red"))

    dataframe <- as.matrix(dataframe)


    Colv  <- dataframe %>% t %>% dist(method = "euclidean") %>% hclust(method = "average") %>% as.dendrogram %>%
               set("branches_k_color", k=2) %>% set("branches_lwd", 2) %>%
               ladderize %>% rotate_DendSer

            #png("Dend.png")
            #plot(Colv)
            #dev.off()


    Rowv  <- dataframe %>% dist %>% hclust(method = "average") %>% as.dendrogram %>%
               set("branches_k_color", k = 2) %>% set("branches_lwd", 2) %>%
               ladderize


    png(filename=paste(this_timepoint,functionOfInterest,"functions.KOvsWT_OIR_heatmap.png", sep="_"), width=2000, height=1500, res = 150)  
    heatmap.2(dataframe, col=my_palette(100), symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, 
              Rowv=Rowv,
              #Rowv=TRUE,
              Colv=Colv,
              #Colv=TRUE,
              srtCol=45,
              revC=TRUE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(2,6), 
              key.par=list(mar=c(4, 4, 4, 3)), 
              margins=c(10,45), key.ylab=NA, main="Specific functions in 
              OIR vs NORM retina")
    dev.off()

}

my_palette <- colorRampPalette(c("blue", "white", "red"))

dataframe.all <- as.matrix(dataframe.all)


Colv  <- dataframe.all %>% t %>% dist(method = "euclidean") %>% hclust(method = "average") %>% as.dendrogram %>%
               set("branches_k_color", k=2) %>% set("branches_lwd", 2) %>%
               ladderize %>% rotate_DendSer

            #png("Dend.png")
            #plot(Colv)
            #dev.off()


Rowv  <- dataframe.all %>% dist %>% hclust(method = "average") %>% as.dendrogram %>%
               set("branches_k_color", k = 2) %>% set("branches_lwd", 2) %>%
               ladderize


png(filename=paste(functionOfInterest,"functions.OIRvsNORM_p14_P17_heatmap.png", sep="_"), width=2400, height=1500, res = 150)  
heatmap.2(dataframe.all, col=my_palette(100), symbreak=TRUE, trace='none', cexRow=1, cexCol= 1, 
              Rowv=Rowv,
              #Rowv=TRUE,
              Colv=FALSE,
              #Colv=TRUE,
              srtCol=45,
              revC=TRUE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(2,6), 
              key.par=list(mar=c(4, 4, 4, 4)), 
              margins=c(30,45), key.ylab=NA, main="Specific functions in 
              OIR vs NORM retina")
dev.off()



q("no")
