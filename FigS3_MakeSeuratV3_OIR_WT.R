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

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/OIR_WT/Clustering")

project_name <- "RETINA_CD31_OIR_WT"


DGE.data=data.frame(fread("/home/gaelcge/projects/def-jsjoyal/gaelcge/Sequencing/CD31_retina/Retina_mix_P14_P17_NORMOIR_WTSirt3KO_CD31only.txt", sep="\t", header=TRUE), row.names=1)

# Subset dataset and matrix

DGE.data <- DGE.data %>% dplyr:: select(grep("OIR.P14.CD31.WT", names(DGE.data)), grep("OIR.P17.CD31.WT", names(DGE.data)))


colnames(DGE.data) = gsub(".", "_", colnames(DGE.data), fixed = TRUE)

# Look at the data matrix
corner(DGE.data)
dim(DGE.data)


###Optional = downsize to 2000

#DGE.data <- DGE.data[ ,sample(names(DGE.data), 2000)]

#Initialize the Seurat object with the raw (non-normalized data)
# Keep all genes expressed in >= 3 cells, keep all cells with >= 100 genes
Seurat_object <- CreateSeuratObject(DGE.data, project = paste0(project_name), assay = "RNA",
  min.cells = 3, min.features = 100)



#nGene and nUMI are automatically calculated for every object by Seurat. For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData. The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.


Seurat_object[["percent.mt"]] <- PercentageFeatureSet(Seurat_object, pattern = "^MT-")


CRYA.genes <- grep(pattern = "^CRYA", x = rownames(x = Seurat_object[['RNA']]), value = TRUE)
CRYB.genes <- grep(pattern = "^CRYB", x = rownames(x = Seurat_object[['RNA']]), value = TRUE)
CRYG.genes <- grep(pattern = "^CRYG", x = rownames(x = Seurat_object[['RNA']]), value = TRUE)
crystal.genes <- c(CRYA.genes, CRYB.genes)
percent.crystal <- Matrix::colSums(Seurat_object[['RNA']][crystal.genes, ])/Matrix::colSums(Seurat_object[['RNA']])

Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.crystal, col.name = "percent.crystal")

#Seurat_object[["percent.crystal"]] <- PercentageFeatureSet(Seurat_object, pattern = crystal.genes)

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats


#Assigning batch to data
Colname_long <- data.frame("merged_name"=colnames(DGE.data))

batch_assigned <- Colname_long %>%
  tidyr::separate(merged_name,c("Treatment", "Age", "Sorting", "Genotype", "Strain", "Lab", "Replicate", "Cell"),sep="_")

batch_assigned <- batch_assigned %>%
  mutate("Batch"=paste(Treatment,Age,Genotype,Replicate,sep="_"))

batch_assigned <- batch_assigned %>%
  mutate("Condition"=paste(Treatment,Age,Genotype,sep="_"))

rownames(batch_assigned) <- Colname_long$merged_name

Seurat_object <- AddMetaData(object = Seurat_object, metadata = batch_assigned)

#Seurat_object$Condition <- factor(Seurat_object$Condition, levels = c('NORM_P14_WT', 'OIR_P14_WT', "OIR_P14_Sirt3KO", 'NORM_P17_WT', 'OIR_P17_WT', "OIR_P17_Sirt3KO"))



png(filename="VlnPlotQC.Batch.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.crystal"), group.by = "Batch")
dev.off()

png(filename="VlnPlotQC.Conditions.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.crystal"), group.by = "Condition")
dev.off()

#GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage, and also low UMI content, we filter these as well

plot1 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="percentMitoPlot.nGenePlot.png", width=1300, height=600, bg = "white", res = 150)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#We filter out cells that have unique gene counts over 6000
#Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nGene but on anything in the object (as in GenePlot above)


Seurat_object <- subset(Seurat_object, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA < 10000 & percent.mt < 10 & percent.crystal < 0.01)



png(filename="VlnPlotQC.Batch.subset.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.crystal"), group.by = "Batch")
dev.off()

png(filename="VlnPlotQC.Conditions.subset.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.crystal"), group.by = "Condition")
dev.off()



#Cell Cycle scoring

cc.genes <- readLines(con="/home/gaelcge/projects/def-jsjoyal/gaelcge/Seurat_ressource/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]


Seurat_object <- CellCycleScoring(Seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

phase_counts <- Seurat_object@meta.data %>%
  group_by(Condition) %>%
  count(Phase) %>%
  mutate(freq = n / sum(n))  

phase_counts$Phase  = factor(phase_counts$Phase, levels=rev(c("G1", "S", "G2M")))

png(filename=paste("CellCycle_proportions_Condition.png", sep="."), width=1200, height=1000, res = 150) 
ggplot(phase_counts,aes(x=Phase,freq,fill=Phase)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Condition, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


phase_counts <- phase_counts %>% group_by(Condition) %>% mutate(freq_sum = sum(freq))

png(filename=paste("CellCycle_proportions_Condition.stacked.png", sep="."), width=1000, height=800, res=300)
ggplot(phase_counts,aes(x=Condition,y=freq,fill=Phase)) +
  geom_bar(stat="identity",col="black") +
  theme_light() +
  labs(y='Relative Proportion',x='Condition') + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_fill_brewer(palette="PuBuGn")
dev.off()



### Subset for testing

#Seurat_object <- SubsetData(Seurat_object, assay = NULL, cells = NULL,
#  subset.name = NULL, ident.use = NULL, ident.remove = NULL,
#  low.threshold = -Inf, high.threshold = Inf, accept.value = NULL,
#  max.cells.per.ident = 50, random.seed = 1)

###

Seurat_object <- SCTransform(Seurat_object, verbose = FALSE, vars.to.regress = c("nFeature_RNA", "percent.mt", "Batch", "S.Score", "G2M.Score"))


plot1 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="percentMitoPlot.nGenePlot.SCT.png", width=900, height=600, bg = "white", res = 50)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#saveRDS(Seurat_object, paste(project_name, "Seurat_object.rds", sep="."))

#Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(Seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(filename="VariableFeaturePlot.png", width=2500, height=600, bg = "white", res = 150)
CombinePlots(plots = list(plot1, plot2))
dev.off()

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
								resolution = 0.5, 
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
DimPlot(Seurat_object, reduction = "tsne", label=FALSE, group.by = "umap")
dev.off()

png(filename="tsne.sct_Condition.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "tsne", label=FALSE, group.by = "Condition")
dev.off()


png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1"))
dev.off()

#doublet detection

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
#sweep.res.list_Seurat_object <- paramSweep_v3(Seurat_object, PCs = 1:10, sct = FALSE)
#sweep.stats_Seurat_object <- summarizeSweep(sweep.res.list_Seurat_object, GT = FALSE)
#bcmvn_Seurat_object <- find.pK(sweep.stats_Seurat_object)

#
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(Seurat_object@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*length(rownames(Seurat_object@meta.data)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## Run DoubletFinder with varying classification stringencies ---------------------

Seurat_object <- doubletFinder_v3(Seurat_object, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

summary(Seurat_object@meta.data)

png(filename="FeaturePlot_QC_doublet.sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "pANN_0.25_0.09_328"))
dev.off()

png(filename="FeaturePlot_VascularCells.sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1", "KCNJ8", "RGS5"))
dev.off()


png(filename="umap.corrected.doublet.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "DF.classifications_0.25_0.09_328")
dev.off()

Seurat_object <- SubsetData(
  Seurat_object,
  assay = NULL,
  cells = NULL,
  subset.name = "DF.classifications_0.25_0.09_328",
  ident.use = NULL,
  ident.remove = NULL,
  low.threshold = -Inf,
  high.threshold = Inf,
  accept.value = "Singlet",
  max.cells.per.ident = Inf,
  random.seed = 1)

png(filename="umap.corrected.doublet_filtered.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "DF.classifications_0.25_0.09_328")
dev.off()


# find markers for every cluster compared to all remaining cells, report only the positive ones

DefaultAssay(Seurat_object) <- "RNA"

Seurat_object <- NormalizeData(Seurat_object, verbose = FALSE)

Seurat_object <- ScaleData(Seurat_object, verbose = FALSE, vars.to.regress = c("nFeature_RNA", "percent.mt", "Batch", "S.Score", "G2M.Score"))


markers.retina.dotplot <- rev(c("RHO", "OPN1SW", "TRPM1", "SNHG11", "PAX6", "KCNJ8", "FBN1", "CLDN5", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A"))


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


Seurat_object.markers <- FindAllMarkers(Seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Seurat_object.markers, "Seurat_object.cluster.markers.txt")

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

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

png(filename="DotPlot_markers.top10.sct.png", res = 100, width=4000, height=1500)
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


### Assigning annotation
Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="SCT_snn_res.0.5")

new.cluster.ids <- c("ECs",      # cluster 0
                      "Rods",     # cluster 1
                      "ECs",     # cluster 2
                      "Mural cells",     # cluster 3
                      "ECs",     # cluster 4
                      "Immune cells",     # cluster 5
                      "ECs",     # cluster 6
                      "ECs",     # cluster 7
                      "Rods",     # cluster 8
                      "Neuronal cells",     # cluster 9
                      "Immune cells",     # cluster 10
                      "Neuronal cells",     # cluster 11
                      "Neuronal cells",     # cluster 12
                      "Neuronal cells",     # cluster 13
                      "ECs",     # cluster 14
                      "ECs",     # cluster 15
                      "Glial cells",     # cluster 16
                      "Cones",     # cluster 17
                      "ECs",     # cluster 18
                      "RBCs",     # cluster 19
                      "Neuronal cells",     # cluster 20
                      "Lens cells",     # cluster 21
                      "Cones",     # cluster 22
                      "Glial cells"     # cluster 23
                      
                      )    



names(new.cluster.ids) <- levels(Seurat_object)
Seurat_object <- RenameIdents(Seurat_object, new.cluster.ids)

Seurat_object[["Cell_Type"]] <- Idents(object = Seurat_object)

png(filename="umap.sct_filtered.annotated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

png(filename="umap.sct.splitedbyCondition.annotated.png", width=2500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Condition")
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


png(filename=paste("Cell_Type_proportions_Condition_percelltype.png", sep="."), width=3500, height=3000, res = 300) 
ggplot(phase_counts,aes(x=Condition,y=freq,fill=Condition)) +
  geom_bar(stat="identity",col="black") +
  facet_wrap(~ Cell_Type, scale="free") +
  theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

 
exprs <- GetAssayData(object = Seurat_object, assay = "RNA", slot = "data")


write.csv(exprs, "./GSVA/NotImputed/exprs.data.NonImputed.csv")


saveRDS(Seurat_object, paste(project_name, "Seurat_object.rds", sep="."))

Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))


q("no")
