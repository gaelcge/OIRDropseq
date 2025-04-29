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
#library(DoubletFinder)
library(harmony)
library(sctransform)
library(future)
plan("multiprocess", workers = (availableCores()-1))
options(future.globals.maxSize = 50000 * 1024^2)

# Load the dataset

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/OIR_NORM_TimeCourse/Aligned/SeuratV3/Aligned_Sorting/Clustering")

project_name <- "RETINA_RYTVELA_OIRTimeCourse_AlignedBySorting"

res = 1

DIM_nb <- c(1:20)

Dim <- 20

perp = 30


DGE.data=data.frame(fread("/home/gaelcge/projects/def-jsjoyal/gaelcge/Sequencing/Merging/TimeCourseOIR/Retina_NORM-OIRTimeCourse_WR_Cd73ft.txt", sep="\t", header=TRUE), row.names=1)

# Subset dataset and matrix

#retina.data <- retina.data[ ,sample(names(retina.data), 500)]

DGE.data.Mccarrol <- DGE.data %>% dplyr:: select(grep("Mccarrol_r3", names(DGE.data)),grep("Mccarrol_r5", names(DGE.data)))

#write.table(DGE.data.Mccarrol, "Dropseq_p14_retina_Mccarroll.txt", quote=FALSE, row.names=TRUE, sep="\t", col.names=TRUE)


#DGE.data.Mccarrol <- DGE.data.Mccarrol[ ,sample(names(DGE.data.Mccarrol), 6000)]

corner(DGE.data.Mccarrol)
dim(DGE.data.Mccarrol)

DGE.data.CHUSJ <- DGE.data %>% dplyr:: select(grep("Joyal", names(DGE.data)))

#write.table(DGE.data.CHUSJ, "Dropseq_p14_retina_Joyal.txt", quote=FALSE, row.names=TRUE, sep="\t", col.names=TRUE)

#DGE.data.CHUSJ_OIR <- DGE.data %>% dplyr:: select(grep("OIR.P14.WR.Joyal.r1", names(DGE.data)))

#write.table(DGE.data.CHUSJ_OIR, "Dropseq_p14_retina_retinopathy_Joyal.txt", quote=FALSE, row.names=TRUE, sep="\t", col.names=TRUE)

corner(DGE.data.CHUSJ)
dim(DGE.data.CHUSJ)

DGE.data <- merge(DGE.data.Mccarrol, DGE.data.CHUSJ, by=0, all=TRUE)

rownames(DGE.data) <- DGE.data$Row.names


DGE.data <- DGE.data %>% dplyr:: select(grep("NORM", names(DGE.data)), grep("OIR", names(DGE.data)))

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

mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
percent.mito <- Matrix::colSums(Seurat_object@assays$RNA[mito.genes, ])/Matrix::colSums(Seurat_object@assays$RNA)

CRYA.genes <- grep(pattern = "^CRYA", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
CRYB.genes <- grep(pattern = "^CRYB", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
CRYG.genes <- grep(pattern = "^CRYG", x = rownames(x = Seurat_object@assays$RNA), value = TRUE)
crystal.genes <- c(CRYA.genes, CRYB.genes)
percent.crystal <- Matrix::colSums(Seurat_object@assays$RNA[crystal.genes, ])/Matrix::colSums(Seurat_object@assays$RNA)


Seurat_object[["percent.crystal"]] <- PercentageFeatureSet(Seurat_object, pattern = crystal.genes)

#AddMetaData adds columns to object@data.info, and is a great place to stash QC stats

Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.mito, col.name = "percent.mito")

Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.crystal, col.name = "percent.crystal")

#Assigning batch to data
batch_long <- data.frame("merged_name"=colnames(DGE.data))

batch_assigned <- batch_long %>% separate(merged_name,c("Condition", "TimePoint", "Sorting", "Labo", "Replicate", "CellBarCode"),sep="_")

batch_assigned <- batch_assigned %>%
  mutate("Batch"=paste(Condition,TimePoint,Sorting,Replicate,sep="_"))

batch_assigned <- batch_assigned %>%
  mutate("Cond_Sorting"=paste(Condition,Sorting,sep="_"))

batch_assigned <- batch_assigned %>%
  mutate("Cond_TimePoint"=paste(Condition,TimePoint,sep="_"))
 
rownames(batch_assigned) <- batch_long$merged_name

Seurat_object <- AddMetaData(object = Seurat_object, metadata = batch_assigned)

png(filename="VlnPlotQC.Cond_Sorting.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal"), group.by = "Cond_Sorting")
dev.off()

png(filename="VlnPlotQC.Conditions.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal"), group.by = "Condition")
dev.off()

png(filename="VlnPlotQC.Sorting.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal"), group.by = "Sorting")
dev.off()

#GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage, and also low UMI content, we filter these as well

plot1 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="percentMitoPlot.nGenePlot.png", width=900, height=600, bg = "white", res = 50)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#We filter out cells that have unique gene counts over 6000
#Note that accept.high and accept.low can be used to define a 'gate', and can filter cells not only based on nGene but on anything in the object (as in GenePlot above)


Seurat_object <- subset(Seurat_object, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA < 10000 & percent.mito < 0.10 & percent.crystal < 0.025)


#Cell Cycle scoring

cc.genes <- readLines(con="/home/gaelcge/projects/def-jsjoyal/gaelcge/Seurat_ressource/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]


Seurat_object <- CellCycleScoring(Seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

### Subset for testing

#Seurat_object <- SubsetData(Seurat_object, assay = NULL, cells = NULL,
#  subset.name = NULL, ident.use = NULL, ident.remove = NULL,
#  low.threshold = -Inf, high.threshold = Inf, accept.value = NULL,
#  max.cells.per.ident = 50, random.seed = 1)

###

Seurat_object.list <- SplitObject(Seurat_object, split.by = "Sorting")

for (i in 1:length(Seurat_object.list)) {
    Seurat_object.list[[i]] <- SCTransform(Seurat_object.list[[i]], verbose = FALSE, vars.to.regress = c("nFeature_RNA", "percent.mito", "Batch", "S.Score", "G2M.Score"))
}

#Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.

Seurat_object.features <- SelectIntegrationFeatures(object.list = Seurat_object.list, nfeatures = 3000)

Seurat_object.list <- PrepSCTIntegration(object.list = Seurat_object.list, anchor.features = Seurat_object.features, 
    verbose = FALSE)

#Next, identify anchors and integrate the datasets. Commands are identical to the standard workflow, but make sure to set  normalization.method = 'SCT':

Seurat_object.anchors <- FindIntegrationAnchors(object.list = Seurat_object.list, normalization.method = "SCT", 
    anchor.features = Seurat_object.features, verbose = FALSE, dims = DIM_nb)

Seurat_object.integrated <- IntegrateData(anchorset = Seurat_object.anchors, normalization.method = "SCT", 
    verbose = FALSE)

plot1 <- FeatureScatter(Seurat_object.integrated, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(Seurat_object.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="percentMitoPlot.nGenePlot.integrated.png", width=900, height=600, bg = "white", res = 50)
CombinePlots(plots = list(plot1, plot2))
dev.off()


saveRDS(Seurat_object.integrated, paste(project_name, "Seurat_object.integrated.rds", sep="."))



q("no")
