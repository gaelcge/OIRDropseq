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
library(doubletFinder)
library(harmony)
library(sctransform)
library(future)


Sys.unsetenv("GITHUB_PAT")

x <- c("Augur", "Seurat", "dplyr", "tibble", "magrittr", "purrr", "tester", "Matrix", 
  "sparseMatrixStats", "parsnip", "recipes", "rsample", "yardstick",
   "pbmcapply", "lmtest", "rlang", "glmnet", "randomForest", "future", "parallel")
lapply(x, library, character.only = TRUE)


plan("multiprocess", workers = (availableCores()-1))
options(future.globals.maxSize = 3000 * 1024^2)
set.seed(42)


#Set directory of dataset to analyse
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/Mapping")

project_name <- "RETINA_CD31_ECs_Integrated"

Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))

Seurat_object <- RenameIdents(object = Seurat_object, "Astrocytes" = "Macroglia")
Seurat_object <- RenameIdents(object = Seurat_object, "Muller Glial cells" = "Macroglia")

Seurat_object[["cell_type"]] <- Seurat_object@active.ident 
Seurat_object[["label"]] <- Seurat_object@meta.data$Dataset




##SUbset condition and cells of interest and set default assay

Seurat_object_Subset <- subset(Seurat_object, 
  idents = c("RBCs", "Lens cells"), invert = T)

Seurat_object_Subset <- subset(Seurat_object_Subset, subset = Dataset == "NORM_WT", invert = T)

Seurat_object_Subset <- subset(Seurat_object_Subset, subset = Age == "P17")


table(Seurat_object_Subset@active.ident)

Seurat_object_Subset[["cell_type_label"]] <- paste(Seurat_object_Subset@meta.data$cell_type,Seurat_object_Subset@meta.data$label,sep="_") 


DefaultAssay(Seurat_object_Subset) <- "RNA"


#Seurat_object_Subset <- subset(Seurat_object_Subset, 
#  idents = c("Pericytes", "Endothelial_cells", "Astrocytes", "Immune_cells"))


setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/Augur")

colorsplot <- DiscretePalette(15, palette = "polychrome")


setEPS(bg = "white", family = "Times", width=12, height=6)
postscript("umap.annotated.Cell_Type_Condition_P17_filtered.eps")
DimPlot(Seurat_object_Subset, reduction = "umap", label = FALSE, pt.size = 0.1, 
  group.by = "cell_type", split.by="Condition", cols = colorsplot)
dev.off()




###RUn augur

augur = calculate_auc(Seurat_object_Subset, n_threads=5, classifier = "rf")

pdf("AUGUR_P17_RNA_Sirt3KOvsWT_OIR_lr_lollipop.pdf")
plot_lollipop(augur)
dev.off()


pdf("augur_P17_RNA_Sirt3KOvsWT_OIR_umap.YlGnBu.pdf", width=3, height=3)
plot_umap(augur, Seurat_object_Subset, palette = "YlGnBu")
dev.off()




top2 <- augur$feature_importance %>% group_by(cell_type) %>% top_n(n = 5, wt = importance)


png(filename="DotPlot_markers.top2.rna.P17_RNA_Sirt3KOvsWT_OIR.png", res = 150, width=2500, height=1000)
DotPlot(
  Seurat_object_Subset,
  assay = NULL,
  unique(top2$gene),
  cols = c("blue", "red"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  group.by = "cell_type_label",
  split.by = NULL,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) + theme(axis.text.x = element_text(angle = 90))
dev.off()




q("no")
