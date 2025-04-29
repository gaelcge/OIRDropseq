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
library(DoubletFinder)
library(harmony)
library(sctransform)
library(future)
library(Hmisc)
plan("multiprocess", workers = (availableCores()-1))
options(future.globals.maxSize = 3000 * 1024^2)



# Load the dataset

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/NORM_WT/Subclustering/ECs/No_CC_Correction")

project_name <- "RETINA_CD31_NORM_WT_ECs"

Seurat_object_NORM_WT <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

Seurat_object_NORM_WT[["Dataset"]] <- "NORM_WT"

Seurat_object_NORM_WT <- SubsetData(
  Seurat_object_NORM_WT,
  assay = NULL,
  cells = NULL,
  subset.name = NULL,
  ident.use = NULL,
  ident.remove = c("Rod contam"),
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1)


setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/OIR_WT/Subclustering/ECs/No_CC_Correction")

project_name <- "RETINA_CD31_OIR_WT_ECs"

Seurat_object_OIR_WT <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

Seurat_object_OIR_WT[["Dataset"]] <- "OIR_WT"

Seurat_object_OIR_WT <- SubsetData(
  Seurat_object_OIR_WT,
  assay = NULL,
  cells = NULL,
  subset.name = NULL,
  ident.use = NULL,
  ident.remove = c("Mural contam"),
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1)


setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/Mapping_WTonly/Non_integrated_2")

project_name <- "RETINA_CD31_ECs_Integrated_WTonly"

Seurat_object <- merge(Seurat_object_NORM_WT, y = Seurat_object_OIR_WT, add.cell.ids = NULL, project = project_name)

Seurat_object@meta.data$Condition <- factor(Seurat_object@meta.data$Condition, levels = c("NORM_P14_WT", "NORM_P17_WT",
                                              "OIR_P14_WT", "OIR_P17_WT"))

#Normalizing the data
Seurat_object <- NormalizeData(Seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

#Scaling the data
Seurat_object$CC.Difference <- Seurat_object$S.Score - Seurat_object$G2M.Score

Seurat_object <- ScaleData(Seurat_object, vars.to.regress = c("nFeature_RNA", "percent.mt", "CC.Difference"))

#Feature Selection
Seurat_object<- FindVariableFeatures(Seurat_object, selection.method = "vst", nfeatures = 2000)



###

Seurat_object <- RunPCA(Seurat_object, features = VariableFeatures(object = Seurat_object))

# Examine and visualize PCA results a few different ways
print(Seurat_object[["pca"]], dims = 1:5, nfeatures = 5)


png(filename="VizDimLoadings.sct.png", width=1500, height=900, bg = "white", res = 150)
VizDimLoadings(Seurat_object, dims = 1:2, reduction = "pca")
dev.off()

png(filename="DimPlot.Condition.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "pca", group.by ="Condition")
dev.off()

png(filename="DimHeatmap.sct.png", width=1500, height=3000, bg = "white", res = 150)
DimHeatmap(Seurat_object, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

png(filename="ElbowPlot.sct.png", width=1500, height=1000, bg = "white", res = 150)
ElbowPlot(Seurat_object)
dev.off()


png(filename="DimPlot.cellcycle.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "pca", group.by ="Phase")
dev.off()

##Run non-linear dimensional reduction (UMAP/tSNE)


Seurat_object <- RunUMAP(Seurat_object, dims = 1:20,
                            n.components = 2L
                            ,min.dist = 0.3
                            ,n.neighbors = 30L
                            )


##Cluster the cells

Seurat_object <- FindNeighbors(Seurat_object, 
								dims = 1:2, 
                reduction = "umap"
								)

Seurat_object <- FindClusters(Seurat_object, 
								resolution = 1.8, 
                reduction = "umap"
								)

png(filename="umap.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, group.by = "seurat_clusters")
dev.off()

png(filename="umap.sct_Cell_Subtype.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Cell_SubType")
dev.off()


png(filename="umap.sct.splitedbyCondition.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, group.by = "seurat_clusters", split.by = "Condition")
dev.off()

png(filename="umap.sct_Condition.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Condition")
dev.off()

png(filename="umap.sct_Dataset.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Dataset")
dev.off()


png(filename="umap.sct_Batch.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Batch")
dev.off()

png(filename="umap.sct_Age.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Age")
dev.off()

png(filename="umap.sct_Phase.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Phase")
dev.off()


##Run Tsne

Seurat_object <- RunTSNE(Seurat_object, dims = 1:17)



png(filename="tsne.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "tsne", label=TRUE)
dev.off()

png(filename="tsne.sct_phase.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "tsne", label=FALSE, group.by = "Phase")
dev.off()

png(filename="tsne.sct_Condition.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "tsne", label=FALSE, group.by = "Condition")
dev.off()

png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


##Remove unwanted cluster


Seurat_object <- SubsetData(
  Seurat_object,
  assay = NULL,
  cells = NULL,
  subset.name = NULL,
  ident.use = NULL,
  ident.remove = c(35),
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1)


#doublet detection

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
#sweep.res.list_Seurat_object <- paramSweep_v3(Seurat_object, PCs = 1:10, sct = FALSE)
#sweep.stats_Seurat_object <- summarizeSweep(sweep.res.list_Seurat_object, GT = FALSE)
#bcmvn_Seurat_object <- find.pK(sweep.stats_Seurat_object)

#

# find markers for every cluster compared to all remaining cells, report only the positive ones
###Detect know marker

markers.retina.dotplot <- rev(c("RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "PAX6", "KCNJ8", "FBN1", "CLDN5", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A"))


png(filename="DotPlot_markers.retina.sct.png", res = 150, width=1000, height=1500)
DotPlot(
  Seurat_object,
  assay = NULL,
  markers.retina.dotplot,
  cols = c("blue", "red"),
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

markers.ECs_Contam.dotplot <- rev(c("RHO","RPE65", "OPN1SW", "TRPM1", "SNHG11", "PAX6", "KCNJ8", "FBN1", "LYZ2", "RLBP1", "GFAP", "OPTC"))


png(filename="DotPlot_markers.ECs_Contam.sct.png", res = 150, width=1000, height=1500)
DotPlot(
  Seurat_object,
  assay = NULL,
  markers.ECs_Contam.dotplot,
  cols = c("blue", "red"),
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


png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1"))
dev.off()


png("LargeArterialFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("BMX", "GKN3", "EFNB2", "MGP", "CYTL1", "FBLN5"), cols=c("blue", "red"))+RotatedAxis()
dev.off()


png("ArterialFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("GKN3", "HEY1", "EDN3"), cols=c("blue", "red"))+RotatedAxis()
dev.off()

png("ArterialCapillaryFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("GLUL", "SLC26A10"), cols=c("blue", "red"))+RotatedAxis()
dev.off()

png("CapillaryFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("MFSD2A", "TFRC", "RGCC", "KDR", "CXCL12", "SPOCK2"), cols=c("blue", "red"))+RotatedAxis()
dev.off()


png("CapillaryVeinFeaturePlot.1.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("MFSD2A", "TFRC", "RGCC", "KDR", "CAR4", "ITM2A", "GATM"), cols=c("blue", "red"))+RotatedAxis()
dev.off()

png("CapillaryVeinFeaturePlot.2.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("HMCN1", "SLC7A1", "NKD1", "VCAM1", "IER3", "PLTP"), cols=c("blue", "red"))+RotatedAxis()
dev.off()


png("VeinsFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("NR2F2", "VWF", "EPHB4", "TMSB10", "ICAM1", "DARC"), cols=c("blue", "red"))+RotatedAxis()
dev.off()


png("StalkFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("APLNR", "JAG1", "FLT1"), cols=c("blue", "red"))+RotatedAxis()
dev.off()

png("TipFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("APLN", "ESM1", "ANGPT2", "TRP53I11"), cols=c("blue", "red"))+RotatedAxis()
dev.off()

png("ECSubtypeFeaturePlot.png", width=1500, height =1500)
FeaturePlot(Seurat_object, features=c("AQP1", "BMX", "GJA4", "FBLN2", "PTN", "MGP", "CXCL12", "SLC6A6", "GLUL", "APOD", "FOS", "FOSB", "ATF3", "BTG2", "APOLD1", "EGR1", "APLN", "SERPINE1", "ESM1", "ANGPT2", "PRSS23", "TRP53I11", "TOP2A", "MKI67",
"BIRC5"), cols=c("blue", "red")) + RotatedAxis()
dev.off()


##Find markers

Seurat_object.markers <- FindAllMarkers(Seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Seurat_object.markers, "Seurat_object.cluster.markers.txt")

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


png(filename="DotPlot_markers.top10.sct.png", res = 150, width=4500, height=1500)
DotPlot(
  Seurat_object,
  assay = NULL,
  unique(top10$gene),
  cols = c("blue", "red"),
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


png(filename="VlnPlot_top2Markers.sct.png", width=3500, height=4500, bg = "white", res = 150)
VlnPlot(Seurat_object, features = unique(top2$gene))
dev.off()

png(filename="FeaturePlot_top2Markers.sct.png", width=3000, height=4000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c(unique(top2$gene), "RPL13"))
dev.off()

png(filename="DoHeatmap_top10Markers.sct.png", width=2000, height=3000, bg = "white", res = 150)
DoHeatmap(Seurat_object, features = top10$gene) + NoLegend()
dev.off()

png(filename="DoHeatmap_top2Markers.sct.png", width=900, height=700, bg = "white", res = 150)
DoHeatmap(Seurat_object, features = top2$gene) + NoLegend()
dev.off()



### Assigning annotation
Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")

new.cluster.ids <- c("Vein ECs",      # cluster 0
                      "Vein ECs",     # cluster 1
                      "Vein ECs",     # cluster 2
                      "Capillary ECs",     # cluster 3
                      "Capillary ECs",     # cluster 4
                      "Arterial ECs",     # cluster 5
                      "Capillary ECs",     # cluster 6
                      "Vein ECs",     # cluster 7
                      "Arterial ECs",     # cluster 8
                      "Vein ECs",     # cluster 9
                      "Tip ECs",     # cluster 10
                      "Capillary ECs",     # cluster 11
                      "Vein ECs",     # cluster 12
                      "Vein ECs",     # cluster 13
                      "Tip ECs",     # cluster 14
                      "Vein ECs",     # cluster 15
                      "Capillary ECs",     # cluster 16
                      "Capillary ECs",     # cluster 17
                      "Arterial ECs",     # cluster 18
                      "Vein ECs",     # cluster 19
                      "Proliferative ECs",     # cluster 20
                      "Proliferative ECs",     # cluster 21
                      "Vein ECs" ,     # cluster 22
                      "Vein ECs",     # cluster 23
                      "Vein ECs",     # cluster 24
                      "Capillary ECs",     # cluster 25
                      "Capillary ECs",     # cluster 26
                      "Vein ECs",     # cluster 27
                      "Vein ECs",     # cluster 28
                      "Vein ECs",     # cluster 29
                      "Vein ECs" ,     # cluster 30
                      "Capillary ECs",     # cluster 31
                      "Tip ECs",     # cluster 32
                      "Proliferative ECs",     # cluster 33
                      "Tuft"     # cluster 34
                        )    



names(new.cluster.ids) <- levels(Seurat_object)
Seurat_object <- RenameIdents(Seurat_object, new.cluster.ids)

Seurat_object[["Cell_SubType"]] <- Idents(object = Seurat_object)

png(filename="umap.sct_filtered.annotated.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

setEPS(bg = "white", family = "Times", width=7, height=4)
postscript("umap.sct_filtered.annotated.eps")
DimPlot(Seurat_object, reduction = "umap", label=FALSE)
dev.off()

setEPS(bg = "white", family = "Times", width=7, height=5)
postscript("umap.sct_filtered.condition.annotated.eps")
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Condition")
dev.off()



png(filename="umap.sct.splitedbyCondition.annotated.png", width=3000, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Condition")
dev.off()

setEPS(bg = "white", family = "Times", width=15, height=4)
postscript("umap.sct.splitedbyCondition.annotated.nolabel.eps")
DimPlot(Seurat_object, reduction = "umap", label=FALSE, split.by = "Condition")
dev.off()


png(filename="umap.sct.splitedbyDataset.annotated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Dataset")+ NoLegend()
dev.off()

png(filename="umap.sct.splitedbyDataset.annotated.nolabel.png", width=2500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, split.by = "Dataset")
dev.off()

levels <- c("Arterial ECs", "Capillary ECs", "Vein ECs", "Proliferative ECs", "Tip ECs", "Tuft")

Seurat_object@active.ident <- factor(Seurat_object@active.ident, level =levels)

Seurat_object.markers <- FindAllMarkers(Seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

setEPS(bg = "white", family = "Times", width=20, height=5)
postscript("DotPlot_markers.top10.sct.renamed.eps")
DotPlot(
  Seurat_object,
  assay = NULL,
  c("CLDN5", "CDH5", "PECAM1",unique(top10$gene)),
  cols = c("blue", "red"),
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


phase_counts <- Seurat_object@meta.data %>%
  group_by(Condition) %>%
  count(Cell_SubType) %>%
  mutate(freq = n / sum(n)*100)  


png(filename=paste("Cell_Type_proportions_Condition.png", sep="."), width=3000, height=2000, res = 300) 
ggplot(phase_counts,aes(x=reorder(Cell_SubType,-freq),y=freq,fill=Cell_SubType)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Condition, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


phase_counts <- phase_counts %>% group_by(Condition) %>% mutate(freq_sum = sum(freq))

phase_counts$Cell_SubType  = factor(phase_counts$Cell_SubType, levels=rev(c("Arterial ECs", "Capillary ECs", "Vein ECs", "Proliferative ECs", "Tip ECs", "Tuft")))

ggplot(phase_counts,aes(x=Condition,y=freq,fill=Cell_SubType)) +
  geom_bar(stat="identity",col="black") +
  geom_text(aes(label = paste0(round(freq), "%")), 
              position = position_stack(vjust = 0.5)) +
  theme_light() +
  labs(y='Relative Proportion',x='Condition') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename=paste("Cell_Type_proportions_Condition.stacked.eps", sep="."), width=12, height=12) 



WT_prop <- filter(phase_counts, Condition%in% c("NORM_P14_WT", "NORM_P17_WT"))

WT_prop$Timepoint <- WT_prop$Condition

WT_prop$Timepoint <- gsub("NORM_P14_WT", "14", WT_prop$Timepoint)
WT_prop$Timepoint <- gsub("NORM_P17_WT", "17", WT_prop$Timepoint)

WT_prop$Timepoint <- as.integer(WT_prop$Timepoint)


ggplot(WT_prop,aes(x=Timepoint,y=freq,fill=Cell_SubType)) + 
  geom_area(size=1, colour="black") +
ggsave(filename=paste("Cell_Type_proportions_Condition.stacked_continous_NORM.eps", sep="."), width=12, height=12) 

OIR_prop <- filter(phase_counts, Condition%in% c("OIR_P14_WT", "OIR_P17_WT"))

OIR_prop$Timepoint <- OIR_prop$Condition

OIR_prop$Timepoint <- gsub("OIR_P14_WT", "14", OIR_prop$Timepoint)
OIR_prop$Timepoint <- gsub("OIR_P17_WT", "17", OIR_prop$Timepoint)

OIR_prop$Timepoint <- as.integer(OIR_prop$Timepoint)


ggplot(OIR_prop,aes(x=Timepoint,y=freq,fill=Cell_SubType)) + 
  geom_area(size=1, colour="black") +
ggsave(filename=paste("Cell_Type_proportions_Condition.stacked_continous_OIR.eps", sep="."), width=12, height=12) 




phase_counts <- Seurat_object@meta.data %>%
  group_by(Condition) %>%
  count(Cell_SubType) %>%
  mutate(freq = n / sum(n)*100)  

phase_counts$Condition  = factor(phase_counts$Condition, levels=c("NORM_P14_WT", "OIR_P14_WT", "NORM_P17_WT", "OIR_P17_WT"))


setEPS(bg = "white", family = "Times", width=8, height=8)
postscript(paste("Cell_Type_proportions_Condition_perCell_SubType.YlGnBu.eps", sep=".")) 
ggplot(phase_counts,aes(x=Condition,y=freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_SubType, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette="YlGnBu", direction=1)
dev.off()

ggplot(phase_counts,aes(x="", y=freq, fill=reorder(Cell_SubType,freq))) +
geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  facet_wrap(~ Condition) +
    geom_text(aes(label = paste0(round(freq), "%")), 
              position = position_stack(vjust = 0.5)) +
  coord_polar("y")+
  theme_light()
ggsave(filename=paste("Cell_Type_proportions_Condition_perNew_Cell_SubType.YlGnBu.pie.eps", sep="."), width=12, height=12) 



 
exprs <- GetAssayData(object = Seurat_object, assay = "RNA", slot = "data")


write.csv(exprs, "./GSVA/NotImputed/exprs.data.NonImputed.csv")


saveRDS(Seurat_object, paste(project_name, "Seurat_object.rds", sep="."))



###Reload object

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/Mapping_WTonly/Non_integrated_2")

project_name <- "RETINA_CD31_ECs_Integrated_WTonly"

Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

Seurat_object <- BuildClusterTree(
  Seurat_object,
  assay = NULL,
  features = VariableFeatures(Seurat_object),
  dims = NULL,
  graph = NULL,
  slot = "data",
  reorder = FALSE,
  reorder.numeric = FALSE,
  verbose = TRUE
)

setEPS(bg = "white", family = "Times", width=8, height=8)
postscript(paste("BuildClusterTree.eps", sep=".")) 
PlotClusterTree(Seurat_object)
dev.off()


##Find narker for tuft in OIR P17

Seurat_object_P17_OIR <- SetIdent(Seurat_object, value="Condition")

Seurat_object_P17_OIR <- SubsetData(
  Seurat_object_P17_OIR,
  assay = NULL,
  cells = NULL,
  subset.name = NULL,
  ident.use = "OIR_P17_WT",
  ident.remove = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1)

Seurat_object_P17_OIR <- SetIdent(Seurat_object_P17_OIR, value="Cell_SubType")
 
Seurat_object.markers <- FindMarkers(Seurat_object_P17_OIR, ident.1="Tuft", only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)

library(tibble)
top <- Seurat_object.markers %>% rownames_to_column('gene') %>% top_n(n = 50, wt = avg_logFC)

png(filename="DotPlot_tuft_markers.top50.sct.renamed.P17_OIR.png", res = 150, width=2500, height=700)
DotPlot(
  Seurat_object_P17_OIR,
  assay = NULL,
  unique(top$gene),
  cols = c("blue", "red"),
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



png(filename="DotPlot_tuft_markers.top50.sct.renamed.splitbycondition.png", res = 150, width=2500, height=1700)
DotPlot(
  Seurat_object,
  assay = NULL,
  unique(top$gene),
  cols = c("limegreen", "limegreen","magenta4", "magenta4"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = NULL,
  split.by = "Condition",
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()


selectedmarker <- c("HEY1", "SIX3", "NSD1",   # TransFact MSX1-NRF2-SMAD1-MTA1
	 "HADHA" , "SAT1", "SLC2A1", "NSD1", "KMT2A", "MRPL20", "UQCRFS1", #Metabolis
	 "INSR", "JKAMP", "CPE", #Insulin signalling
	 "HBEGF", #Hypertrophy, clathrin
	 "VASP", "PHLDB2",  ## Cell adhesion
	 "COL18A1", # ECM
	 "CTSA",  #Angiotensin, peptidase activity
	 "NOSTRIN", #eNOS
	 "LITAF", #IntraC Sign
	 "AQP1", "PLXND1", #Vasc Dev
	 "EML1", #Cytosk
	 "PRSS23", "MAL", # Interleukine
	 "GIMAP4", "GIMAP5", #GTPase
	 "GM23935" #miRNA
	)


png(filename="DotPlot_tuft_selected_markers.sct.renamed.P17_OIR.png", res = 150, width=2500, height=700)
DotPlot(
  Seurat_object_P17_OIR,
  assay = NULL,
  unique(selectedmarker),
  cols = c("blue", "red"),
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



library('org.Mm.eg.db')
library('org.Hs.eg.db')

GeneID_entrez <- as.data.frame(mapIds(org.Mm.eg.db, capitalize(tolower(rownames(Seurat_object.markers))), 'ENTREZID', 'SYMBOL'))
colnames(GeneID_entrez) <- "GeneID"
rownames(GeneID_entrez) <- toupper(rownames(GeneID_entrez))

Seurat_object.markers_GeneID <- merge(Seurat_object.markers, GeneID_entrez, by="row.names",all.x=TRUE)

write.table(Seurat_object.markers_GeneID, "Seurat_object.markers_GeneID.txt", row.names=FALSE, sep="\t")

q("no")
