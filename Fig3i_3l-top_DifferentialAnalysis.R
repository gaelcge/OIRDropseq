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
plan("multisession", workers = (availableCores()-1))
options(future.globals.maxSize = 3000 * 1024^2)


#Set directory of dataset to analyse
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/Mapping_WTonly/Non_integrated_2")

project_name <- "RETINA_CD31_ECs_Integrated_WTonly"

Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

###R 4.0
library("Nebulosa")
png("plot_density_Tuft.png", width = 4000, height =2000, res = 150)
plot_density(Seurat_object, 
  c("MAL", "NKD1", "PLVAP", "SAT1", "COL18A1", "EML1", "COL15A1", "PRSS23", "SIX3", "SIX3OS1",  "AQP1"), 
  joint = TRUE, combine = TRUE, reduction = "umap", method = c("wkde"), adjust = 2) + 
plot_layout(ncol = 4)
dev.off()

#### For Diff functions between conditions globally

##Analysis global
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/DifferentialExpression/Non_integrated/Functions/OIRvsNORM") 

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Dataset %in% c("NORM_WT","OIR_WT")))

Seurat_object_Subset <- subset(Seurat_object, cells = Subseted_cells)

png(filename="umap.sct_filtered.annotated.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object_Subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()



###Add GO function to Seurat object

GO_score=fread("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/GSVA/WTonly/gsva.exprs.MsigDB_senescence.csv", sep=",",data.table=FALSE)
rownames(GO_score) <- GO_score$V1
GO_score$V1 <- NULL

GO_score_subset = GO_score[,rownames(Seurat_object_Subset@meta.data)]

Seurat_object_Subset[["GO"]] <- CreateAssayObject(counts = GO_score_subset)

Seurat_object_Subset <- NormalizeData(Seurat_object_Subset, assay = "GO", normalization.method = "CLR")

Seurat_object_Subset <- ScaleData(Seurat_object_Subset, assay = "GO", verbose = FALSE, model.use = "negbinom")

Seurat_object_Subset@meta.data$Cell_SubType <- factor(Seurat_object_Subset@meta.data$Cell_SubType, 
                                                  levels = c("Arterial ECs", "Capillary ECs", "Vein ECs", "Proliferative ECs", "Tip ECs", "Tuft")
                                                  )

### Make ridge plot for specific functions acrross cell types
pathways <- "KETO"
pathways <- gsub("_", "-", pathways)

rowf <- grep(pathways, rownames(Seurat_object_Subset[["GO"]]))

selected_functions <- rownames(Seurat_object_Subset[["GO"]])[rowf]

Seurat_object_loop <- SetIdent(Seurat_object_Subset, value = "Condition")

loop_list <- unique(Seurat_object_loop@active.ident)

for(this_cluster in loop_list){

    cells_in_this_cluster <- SubsetData(Seurat_object_loop,
                                      ident.use=this_cluster)

    nbcelltype <- length(unique(cells_in_this_cluster@meta.data$Cell_SubType))

    ridgeplot <- RidgePlot(
    cells_in_this_cluster,
    assay = "GO",
    features=selected_functions,
    cols = viridis(nbcelltype),
    idents = NULL,
    sort = "decreasing",
    group.by = "Cell_SubType",
    y.max = NULL,
    same.y.lims = FALSE,
    log = FALSE,
    ncol = NULL,
    slot = "scale.data",
    combine = TRUE)
    
    ggsave(ridgeplot, filename=paste("Ridgeplot/Ketone/",this_cluster,"_",pathways,"_ridgeplot.new.eps",sep=""), device = "eps", width = 16, height=16)

}

###Make heatmap of celltype specific functions


Subseted_cells <- rownames(subset(Seurat_object_Subset@meta.data, Age %in% c("P17")))

Seurat_object_Subset_Age <- SubsetData(Seurat_object_Subset, cells = Subseted_cells)

DGE_test_P14 <- FindMarkers(Seurat_object_Subset_Age, ident.1="Tuft", logfc.threshold = 0.25, min.pct=0.5, 
    features= NULL, assay = "GO")

DGE_test_P14 <- FindAllMarkers(Seurat_object_Subset_Age, logfc.threshold = 0.20, min.pct=0.1, 
    features= NULL, assay = "GO")

top10 <- DGE_test_P14 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)



png(filename="DoHeatmap_top10Markers.ECs_subtypes_OIRandNORM_P17.png", width=2000, height=2000, bg = "white", res = 150)
DoHeatmap(Seurat_object_Subset, features = top10$gene,assay = "GO") #+ NoLegend()
dev.off()

selected_functions <- unique(top10$gene)

###Compare using doheatmap all cell type in separate condtions



rowf <- unique(c(grep("HALLMARK-MYC-", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-FOXO", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-REGULATION-OF-LOCALIZATION-OF-FOXO", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-REGULATION-OF-FOXO", rownames(Seurat_object_Subset[["GO"]]))
          ))

rowf <- unique(c(grep("HALLMARK-MYC-TARGETS-V1", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-MYC-TARGETS-V2", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-FOXO-MEDIATED-TRANSCRIPTION$", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-REGULATION-OF-LOCALIZATION-OF-FOXO-TRANSCRIPTION-FACTORS", rownames(Seurat_object_Subset[["GO"]])),
          grep("REACTOME-REGULATION-OF-FOXO-TRANSCRIPTIONAL-ACTIVITY-BY-ACETYLATION", rownames(Seurat_object_Subset[["GO"]])),
          grep("GO-HYPOXIA-INDUCIBLE-FACTOR-1ALPHA-SIGNALING-PATHWAY", rownames(Seurat_object_Subset[["GO"]]))
          ))

rowf <- rev(unique(c(grep("GO-AEROBIC-RESPIRATION", rownames(Seurat_object_Subset[["GO"]])),
      grep("REACTOME-MITOCHONDRIAL-FATTY-ACID-BETA-OXIDATION-OF-SATURATED-FATTY-ACIDS", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("HALLMARK-HYPOXIA", rownames(Seurat_object_Subset[["GO"]])),
      grep("GO-SPROUTING-ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]]))
          )))

rowf <- rev(unique(c(grep("ANGIOGENESIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("VASCULAR", rownames(Seurat_object_Subset[["GO"]])),
      grep("SENESCENCE", rownames(Seurat_object_Subset[["GO"]]))
           )))

rowf <- rev(unique(c(grep("BARRIER", rownames(Seurat_object_Subset[["GO"]])),
      grep("PERMEABILITY", rownames(Seurat_object_Subset[["GO"]])),
      grep("JUNCTION", rownames(Seurat_object_Subset[["GO"]]))
          )))

rowf <- rev(unique(c(grep("PROLIFERATION", rownames(Seurat_object_Subset[["GO"]])),
      grep("APOPTOSIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("HYPOXIA", rownames(Seurat_object_Subset[["GO"]]))
          )))

rowf <- rev(unique(c(grep("RESPIRATION", rownames(Seurat_object_Subset[["GO"]])),
      grep("OXIDATIVE-PHOSPHO", rownames(Seurat_object_Subset[["GO"]])),
      grep("TCA", rownames(Seurat_object_Subset[["GO"]])),
      grep("GLYCOLYSIS", rownames(Seurat_object_Subset[["GO"]])),
      grep("FATTY-ACID", rownames(Seurat_object_Subset[["GO"]])),
      grep("GLUCOSE", rownames(Seurat_object_Subset[["GO"]]))
          )))


rowf <- rev(unique(c(grep(gsub("−", "-","REACTOME−TRANSPORT−OF−FATTY−ACIDS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−LONG−CHAIN−FATTY−ACID−TRANSPORT"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−GLUCOSE−CATABOLIC−PROCESS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","KEGG−GLYCOLYSIS−GLUCONEOGENESIS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","HALLMARK−GLYCOLYSIS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−FATTY−ACID−ELONGATION"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","KEGG−BIOSYNTHESIS−OF−UNSATURATED−FATTY−ACIDS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−RESPONSE−TO−FATTY−ACID"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−UNSATURATED−FATTY−ACID−METABOLIC−PROCESS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","KEGG−BIOSYNTHESIS−OF−UNSATURATED−FATTY−ACIDS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−FATTY−ACID−LIGASE−ACTIVITY"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","REACTOME−CITRIC−ACID−CYCLE−TCA−CYCLE"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−AEROBIC−RESPIRATION"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","HALLMARK−OXIDATIVE−PHOSPHORYLATION"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","REACTOME−THE−CITRIC−ACID−TCA−CYCLE−AND−RESPIRATORY−ELECTRON−TRANSPORT"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","REACTOME−PYRUVATE−METABOLISM−AND−CITRIC−ACID−TCA−CYCLE"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","KEGG−CITRATE−CYCLE−TCA−CYCLE"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−ESTABLISHMENT−OF−ENDOTHELIAL−BLOOD−BRAIN−BARRIER"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−ESTABLISHMENT−OF−ENDOTHELIAL−BARRIER"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−REGULATION−OF−VASCULAR−PERMEABILITY"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−REGULATION−OF−MEMBRANE−PERMEABILITY"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","REACTOME−VEGFR2−MEDIATED−VASCULAR−PERMEABILITY"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−ESTABLISHMENT−OF−BLOOD−BRAIN−BARRIER"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−TRANSPORT−ACROSS−BLOOD−BRAIN−BARRIER"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−SPROUTING−ANGIOGENESIS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","REACTOME−CELLULAR−SENESCENCE"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−VASCULAR−ENDOTHELIAL−GROWTH−FACTOR−RECEPTOR−BINDING"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−POSITIVE−REGULATION−OF−CELL−MIGRATION−INVOLVED−IN−SPROUTING−ANGIOGENESIS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","REACTOME−OXIDATIVE−STRESS−INDUCED−SENESCENCE"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−NEGATIVE−REGULATION−OF−CELL−MIGRATION−INVOLVED−IN−SPROUTING−ANGIOGENESIS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−ANGIOGENESIS−INVOLVED−IN−WOUND−HEALING"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−GLIAL−CELL−PROLIFERATION"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","HALLMARK−APOPTOSIS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","KEGG−APOPTOSIS"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","HALLMARK−HYPOXIA"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","REACTOME−VEGFR2−MEDIATED−CELL−PROLIFERATION"), rownames(Seurat_object_Subset[["GO"]])),
      grep(gsub("−", "-","GO−HYPOXIA−INDUCIBLE−FACTOR−1ALPHA−SIGNALING−PATHWAY"), rownames(Seurat_object_Subset[["GO"]]))
          )))

selected_functions <- rownames(Seurat_object_Subset[["GO"]])[rowf]

Seurat_object_loop <- SetIdent(Seurat_object_Subset, value = "Age")

timepoint <- "P17"

Seurat_object_loop <- SubsetData(Seurat_object_loop,
                                      ident.use=timepoint)

#Subseted_cells <- rownames(subset(Seurat_object_loop@meta.data, Cell_SubType %in% c("Tip ECs", "Tuft ECs")))

#Seurat_object_loop <- SubsetData(Seurat_object_loop, cells = Subseted_cells)

Seurat_object_loop <- SetIdent(Seurat_object_loop, value = "Cell_SubType")                                      


AverageExpressionseurat <- AverageExpression(
  Seurat_object_loop, , return.seurat = FALSE,
  #add.ident = "Condition", 
  features = selected_functions,
  assay = "GO",
  slot = "scale.data")

my_palette <- colorRampPalette(c("royalblue4", "royalblue", "white", "orangered", "red4"))

#my_palette <- colorRampPalette(c("blue", "yellow", "orangered", "red"))

dataframe <- as.matrix(AverageExpressionseurat$GO)

#dataframe <- dataframe[,grep("OIR",colnames(dataframe))]


pdf(paste("OIRNORMWT_",timepoint,"Selection_Large_GO_ECs.scale.data_heatmap.pdf"), width=45, heigh=40)  
heatmap.2(dataframe, col=my_palette(100), symbreak=TRUE, trace='none', cexRow=3, cexCol= 4,
          #Rowv=Rowv,
          Rowv=TRUE,
          #Colv=Colv,
          Colv=TRUE,
          srtCol=45,
          revC=TRUE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(2,5), 
          key.par=list(mar=c(5, 8, 8, 8)), 
          margins=c(25,155), key.ylab=NA, main="Markers GO")
dev.off()


###Make heatmap for specific functions acrross cell types

rowf <- unique(c(grep("SPROUTING", rownames(Seurat_object_Subset[["GO"]])),
          grep("RESPIRATION", rownames(Seurat_object_Subset[["GO"]]))
          ))

selected_functions <- rownames(Seurat_object_Subset[["GO"]])[rowf]

selected_functions <- rownames(Seurat_object_Subset[["GO"]])

dataframe.all <- data.frame()

dataframe.all <- dataframe.all[1:length(selected_functions),]

rownames(dataframe.all) <- selected_functions

Seurat_object_loop <- SetIdent(Seurat_object_Subset, value = "Age")

timepoint <- "P17"

Seurat_object_loop <- SubsetData(Seurat_object_loop,
                                      ident.use=timepoint)

Seurat_object_loop <- SetIdent(Seurat_object_loop, value = "Dataset")

Dataset <- "OIR_WT"

Seurat_object_loop <- SubsetData(Seurat_object_loop,
                                      ident.use=Dataset)

Seurat_object_loop <- SetIdent(Seurat_object_loop, value = "Cell_SubType")

for(this_cluster in unique(Seurat_object_loop@active.ident)) {

  if(length(grep(this_cluster, Seurat_object_loop@active.ident,fixed = TRUE)) > 3) {

    print(paste("Working on cluster #",this_cluster,sep=""))

    DGE_test_P14 <- FindMarkers(Seurat_object_loop, ident.1=this_cluster, logfc.threshold = -Inf, min.pct=0.1, 
          features=selected_functions, assay = "GO", only.pos=FALSE)

    DGE_test_P14 <- subset(DGE_test_P14, p_val_adj < 0.05)

    #DGE_test_P14 <- subset(DGE_test_P14, avg_logFC > 0.4 | avg_logFC < -0.4 )

    DGE_test_P14_avg_logFC <- as.data.frame(DGE_test_P14$avg_logFC)

    rownames(DGE_test_P14_avg_logFC) <- rownames(DGE_test_P14)

    colnames(DGE_test_P14_avg_logFC) <- paste(this_cluster, "avg_logFC", timepoint, sep="_")

    dataframe.all <- merge(dataframe.all,DGE_test_P14_avg_logFC, by=0, all=TRUE)

    rownames(dataframe.all) <- dataframe.all$Row.names
   
    dataframe.all$Row.names <- NULL
 
  }
   
}

dataframe <- dataframe.all

##Select for specific subtypes
dataframe <- dataframe[complete.cases(dataframe[ ,6]),]
dataframe <- dataframe[complete.cases(dataframe[ ,3]),]


#Remove row containing only NA
ind <- apply(dataframe, 1, function(x) all(is.na(x)))

dataframe <- dataframe[ !ind, ]

#### make matrix for heatmap

my_palette <- colorRampPalette(c("royalblue4", "royalblue", "white", "orangered", "red4"))

dataframe <- as.matrix(dataframe)


Colv  <- dataframe %>% t %>% dist(method = "euclidean") %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", k=2) %>% set("branches_lwd", 2) %>%
           ladderize %>% rotate_DendSer

        #png("Dend.png")
        #plot(Colv)
        #dev.off()


Rowv  <- dataframe %>% dist %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", 2) %>% set("branches_lwd", 2) %>%
           ladderize


dataframe[is.na(dataframe)] <- 0
dataframe[is.nan(dataframe)] <- 0

dim(dataframe)

png(filename=paste("EC_markers_",timepoint,Dataset,"_GLYCOLYSIS_ECs_heatmap.png", sep=""), width=4500, height=4500, res = 300)  
heatmap.2(dataframe, col=my_palette(100), symbreak=TRUE, trace='none', cexRow=1, cexCol= 1,
          #Rowv=Rowv,
          Rowv=TRUE,
          #Colv=Colv,
          Colv=TRUE,
          srtCol=45,
          revC=TRUE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(2,5), 
          key.par=list(mar=c(5, 8, 8, 8)), 
          margins=c(65,45), key.ylab=NA, main="OIR_WT P17")
dev.off()


###Make heatmap of OIR vs NORM for specific functions acrross cell types

rowf <- unique(c(grep("BETA-OX", rownames(Seurat_object_Subset[["GO"]])),
          grep("BETA-OX", rownames(Seurat_object_Subset[["GO"]]))
          ))

selected_functions <- rownames(Seurat_object_Subset[["GO"]])[rowf]

dataframe.all <- data.frame()

dataframe.all <- dataframe.all[1:length(selected_functions),]

rownames(dataframe.all) <- selected_functions

Seurat_object_loop <- SetIdent(Seurat_object_Subset, value = "Age")

timepoint <- "P17"

Seurat_object_loop <- SubsetData(Seurat_object_loop,
                                      ident.use=timepoint)

Seurat_object_loop <- SetIdent(Seurat_object_loop, value = "Cell_SubType")                                      

loop_list <- unique(Seurat_object_loop@active.ident)

loop_list <- c("Arterial ECs", "Capillary II", "Capillary I", "Vein ECs", "Tip ECs")


#loop_list <- loop_list[-5]

for(this_cluster in loop_list){

  #####this_cluster <- "Tip_ECs"

  ## Print status for which identity is being processed
  print(paste("Working on cluster #",this_cluster,sep=""))

  ## Subset Seurat object to only contain cells from this cluster
  cells_in_this_cluster <- SubsetData(Seurat_object_loop,
                                      ident.use=this_cluster)

  ## Check whether there are cells in both groups, otherwise skip this cluster
        
  cells_in_this_cluster <- SetIdent(cells_in_this_cluster, value = "Treatment")

  DGE_test_P14 <- FindMarkers(cells_in_this_cluster, ident.1 = "OIR", ident.2 = "NORM", logfc.threshold = -Inf, min.pct=-Inf, 
    features=selected_functions, assay = "GO")

  #DGE_test_P14 <- subset(DGE_test_P14, p_val < 0.05)

  DGE_test_P14_avg_logFC <- as.data.frame(DGE_test_P14$avg_logFC)

  rownames(DGE_test_P14_avg_logFC) <- rownames(DGE_test_P14)

  colnames(DGE_test_P14_avg_logFC) <- paste(this_cluster, "avg_logFC_OIR_vs_NORM", timepoint, sep="_")

  dataframe.all <- merge(dataframe.all,DGE_test_P14_avg_logFC, by=0, all=TRUE)

  rownames(dataframe.all) <- dataframe.all$Row.names
 
  dataframe.all$Row.names <- NULL
}

dataframe <- dataframe.all

##Select for specific subtypes
#dataframe <- dataframe[complete.cases(dataframe[ ,5]),]

#Remove row containing only NA
#ind <- apply(dataframe, 1, function(x) all(is.na(x)))

#dataframe <- dataframe[ !ind, ]

#### make matrix for heatmap

my_palette <- colorRampPalette(c("royalblue4", "royalblue", "white", "orangered", "red4"))

dataframe <- as.matrix(dataframe)


Colv  <- dataframe %>% t %>% dist(method = "euclidean") %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", k=2) %>% set("branches_lwd", 2) %>%
           ladderize %>% rotate_DendSer

        #png("Dend.png")
        #plot(Colv)
        #dev.off()


Rowv  <- dataframe %>% dist %>% hclust(method = "average") %>% as.dendrogram %>%
           set("branches_k_color", 2) %>% set("branches_lwd", 2) %>%
           ladderize


dataframe[is.na(dataframe)] <- 0
dataframe[is.nan(dataframe)] <- 0

png(filename=paste("OIR_vs_NORM_",timepoint,"_RESPI_ECs_heatmap.png", sep=""), width=3500, height=3000, res = 200)  
heatmap.2(dataframe, col=my_palette(100), symbreak=TRUE, trace='none', cexRow=1, cexCol= 1,
          #Rowv=Rowv,
          Rowv=TRUE,
          #Colv=Colv,
          Colv=TRUE,
          srtCol=45,
          revC=TRUE, key.title = NA, density.info="none", lhei=c(1,6), lwid=c(2,5), 
          key.par=list(mar=c(5, 8, 8, 8)), 
          margins=c(25,55), key.ylab=NA, main="OIR vs NORM P17")
dev.off()







q("no")
