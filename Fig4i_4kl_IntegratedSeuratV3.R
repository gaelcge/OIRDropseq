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
plan("multiprocess", workers = (availableCores()-1))
options(future.globals.maxSize = 3000 * 1024^2)

# Load the dataset

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/NORM_WT/Clustering")

project_name <- "RETINA_CD31_NORM_WT"

Seurat_object_NORM_WT <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

Seurat_object_NORM_WT[["Dataset"]] <- "NORM_WT"

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/OIR_WT/Clustering")

project_name <- "RETINA_CD31_OIR_WT"

Seurat_object_OIR_WT <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

Seurat_object_OIR_WT[["Dataset"]] <- "OIR_WT"

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/OIR_KO/Clustering")

project_name <- "RETINA_CD31_OIR_Sirt3KO"

Seurat_object_OIR_Sirt3KO <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

Seurat_object_OIR_Sirt3KO[["Dataset"]] <- "OIR_Sirt3KO"


setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/Mapping")

project_name <- "RETINA_CD31_ECs_Integrated"

Seurat_object <- merge(Seurat_object_NORM_WT, y = Seurat_object_OIR_WT, add.cell.ids = NULL, project = project_name)

Seurat_object <- merge(Seurat_object, y = Seurat_object_OIR_Sirt3KO, add.cell.ids = NULL, project = project_name)

Seurat_object@meta.data$Condition <- factor(Seurat_object@meta.data$Condition, levels = c("NORM_P14_WT", "NORM_P17_WT",
																							"OIR_P14_WT", "OIR_P17_WT",
																							"OIR_P14_Sirt3KO", "OIR_P17_Sirt3KO"))

#DefaultAssay(Seurat_object) <- "SCT"

Seurat_object.list <- SplitObject(Seurat_object, split.by = "Dataset")
for (i in names(Seurat_object.list)) {
    Seurat_object.list[[i]] <- SCTransform(Seurat_object.list[[i]], verbose = FALSE)
}

Seurat_object.features <- SelectIntegrationFeatures(object.list = Seurat_object.list, nfeatures = 3000)
Seurat_object.list <- PrepSCTIntegration(object.list = Seurat_object.list, anchor.features = Seurat_object.features)


Seurat_object.anchors <- FindIntegrationAnchors(object.list = Seurat_object.list, normalization.method = "SCT", 
    anchor.features = Seurat_object.features)
Seurat_object.integrated <- IntegrateData(anchorset = Seurat_object.anchors, normalization.method = "SCT")


###

Seurat_object <- Seurat_object.integrated

DefaultAssay(Seurat_object) <- "integrated"

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
                            n.components = 2L)

Seurat_object <- RunTSNE(Seurat_object, dims = 1:20)


##Cluster the cells

Seurat_object <- FindNeighbors(Seurat_object, 
								dims = 1:2, 
								reduction = "umap"
								)

Seurat_object <- FindClusters(Seurat_object, 
								resolution = 0.3, 
                reduction = "umap"
								)

png(filename="umap.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE)
dev.off()

png(filename="umap.sct.splitedbyCondition.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Condition")
dev.off()

png(filename="umap.sct_Condition.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Condition")
dev.off()

png(filename="umap.sct_Dataset.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Dataset")
dev.off()

png(filename="umap.sct_Cell_Type.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Cell_Type")
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


### look at gene markers



png(filename="umap.sct.2.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE)
dev.off()


markers.retina.dotplot <- rev(c("LHX1", "PAX6", "SLC16A6", "VSX2", "SLC5A7", "RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "FBN1", "CLDN5", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A"))


png(filename="DotPlot_markers.retina.sct.initial.png", res = 150, width=1000, height=1500)
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

png(filename="FeaturePlot_retina_markers_sct.png", width=2500, height=2500, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = markers.retina.dotplot)
dev.off()


##Remove unwanted cluster


Seurat_object <- SubsetData(
  Seurat_object,
  assay = NULL,
  cells = NULL,
  subset.name = NULL,
  ident.use = NULL,
  ident.remove = c(13, 16,33, 34),
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = NULL,
  max.cells.per.ident = Inf,
  random.seed = 1)


### Assigning annotation
Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")

new.cluster.ids <- c("ECs",      # cluster 0
                      "Immune cells",     # cluster 1
                      "ECs",     # cluster 2
                      "Rods",     # cluster 3
                      "ECs",     # cluster 4
                      "ECs",     # cluster 5
                      "Mural cells",     # cluster 6
                      "ECs",     # cluster 7
                      "ECs",     # cluster 8
                      "ECs",     # cluster 9
                      "Rods",     # cluster 10
                      "ECs",     # cluster 11
                      "Bipolar cells",     # cluster 12
                      "ECs",     # cluster 13
                      "Amacrine cells",     # cluster 14
                      "Muller Glial cells",     # cluster 15
                      "Rods",     # cluster 16
                      "ECs",     # cluster 17
                      "Bipolar cells",     # cluster 18
                      "Cones",     # cluster 19
                      "Bipolar cells",     # cluster 20
                      "RBCs",     # cluster 21
                      "Mural cells",     # cluster 22
                      "Lens cells",     # cluster 23
                      "Astrocytes"     # cluster 24 
                        )    



names(new.cluster.ids) <- levels(Seurat_object)
Seurat_object <- RenameIdents(Seurat_object, new.cluster.ids)

Seurat_object[["Cell_Type"]] <- Idents(object = Seurat_object)

png(filename="umap.sct_filtered.annotated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

png(filename="umap.sct.splitedbyCondition.annotated.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Condition")
dev.off()



#doublet detection

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
#sweep.res.list_Seurat_object <- paramSweep_v3(Seurat_object, PCs = 1:10, sct = FALSE)
#sweep.stats_Seurat_object <- summarizeSweep(sweep.res.list_Seurat_object, GT = FALSE)
#bcmvn_Seurat_object <- find.pK(sweep.stats_Seurat_object)

#

# find markers for every cluster compared to all remaining cells, report only the positive ones

DefaultAssay(Seurat_object) <- "RNA"

Seurat_object <- NormalizeData(Seurat_object, verbose = FALSE)

Seurat_object <- ScaleData(Seurat_object, verbose = FALSE, vars.to.regress = c("nFeature_RNA", "percent.mt", "Batch", "S.Score", "G2M.Score"))


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

markers.ECs_Contam.dotplot <- rev(c("RPE65", "OPN1SW", "TRPM1", "SNHG11", "PAX6", "KCNJ8", "FBN1", "LYZ2", "RLBP1", "GFAP", "OPTC"))


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


png("NORM.integrated.LargeArterialFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("BMX", "GKN3", "EFNB2", "MGP", "CYTL1", "FBLN5"), cols=c("blue", "red"))+RotatedAxis()
dev.off()


png("NORM.integrated.ArterialFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("GKN3", "HEY1", "EDN3"), cols=c("blue", "red"))+RotatedAxis()
dev.off()

png("NORM.integrated.ArterialCapillaryFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("GLUL", "SLC26A10"), cols=c("blue", "red"))+RotatedAxis()
dev.off()

png("NORM.integrated.CapillaryFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("MFSD2A", "TFRC", "RGCC", "KDR", "CXCL12", "SPOCK2"), cols=c("blue", "red"))+RotatedAxis()
dev.off()


png("NORM.integrated.CapillaryVeinFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("MFSD2A", "TFRC", "RGCC", "KDR", "CAR4", "ITM2A"), cols=c("blue", "red"))+RotatedAxis()
dev.off()

png("NORM.integrated.VeinsFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("NR2F2", "VWF", "EPHB4", "TMSB10", "ICAM1", "DARC"), cols=c("blue", "red"))+RotatedAxis()
dev.off()


png("NORM.integrated.StalkFeaturePlot.png", width=1500, height =1000)
FeaturePlot(Seurat_object, features=c("APLNR", "JAG1", "FLT1"), cols=c("blue", "red"))+RotatedAxis()
dev.off()

png("NORM.integrated.TipFeaturePlot.png", width=1500, height =1000)
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

png("NORM.integrated.Tree.png", width=500, height =500)
Seurat_object <- BuildClusterTree(
  Seurat_object,
  assay = NULL,
  features = VariableFeatures.Seurat.object,
  dims = NULL,
  graph = NULL,
  slot = "data",
  reorder = FALSE,
  reorder.numeric = FALSE,
  verbose = TRUE
)
PlotClusterTree(object = Seurat_object)
dev.off()





phase_counts <- Seurat_object@meta.data %>%
  group_by(Condition) %>%
  count(Cell_Type) %>%
  mutate(freq = n / sum(n)*100)  

#phase_counts$Condition  = factor(phase_counts$Condition, levels=c("day0", "day3", "day10", "day24", "day31"))

png(filename=paste("Cell_Type_proportions_Condition.png", sep="."), width=3000, height=2000, res = 300) 
ggplot(phase_counts,aes(x=reorder(Cell_Type,-freq),y=freq,fill=Cell_Type)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Condition, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


phase_counts <- phase_counts %>% group_by(Condition) %>% mutate(freq_sum = sum(freq))

png(filename=paste("Cell_Type_proportions_Condition.stacked.png", sep="."), width=2000, height=1500, res=300)
ggplot(phase_counts,aes(x=Condition,y=freq,fill=Cell_Type)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  labs(y='Relative Proportion',x='Condition') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



phase_counts <- Seurat_object@meta.data %>%
  group_by(Condition) %>%
  count(Cell_Type) %>%
  mutate(freq = n / sum(n)*100)  


png(filename=paste("Cell_Type_proportions_Condition_perCell_SubType.png", sep="."), width=3500, height=3000, res = 300) 
ggplot(phase_counts,aes(x=Condition,y=freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_Type, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

 
exprs <- GetAssayData(object = Seurat_object, assay = "RNA", slot = "data")


write.csv(exprs, "./GSVA/NotImputed/exprs.data.NonImputed.csv")


saveRDS(Seurat_object, paste(project_name, "Seurat_object.rds", sep="."))

###

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/Mapping")

project_name <- "RETINA_CD31_ECs_Integrated"

Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))


table(Seurat_object[["Cell_Type"]])


Seurat_object[["Cell_Type_Condition"]] <- paste(Seurat_object@meta.data$Cell_Type,Seurat_object@meta.data$Condition, sep="_")
Seurat_object[["Cell_Type_Dataset"]] <- paste(Seurat_object@meta.data$Cell_Type,Seurat_object@meta.data$Dataset, sep="_")


table(Seurat_object[["Cell_Type_Condition"]])

q("no")
