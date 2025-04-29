rm(list = ls())
library(plotly)
library(ggrepel)
library(gtools)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(useful)
library(Matrix)
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(scico)
library(DropletUtils)
library(devtools)
library(velocyto.R)
library(SeuratWrappers)
theme_set(theme_bw())
library(future)
plan("multicore", workers = (availableCores()-1))
options(future.globals.maxSize = 3000 * 1024^2)



setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/References/kb_python_velo/mouse/")
tr2g <- read_tsv("t2g.txt", col_names = c("transcript", "gene", "gene_name"))
tr2g <- distinct(tr2g[, c("gene", "gene_name")])


read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}



#Create seurat object with spliced and unspliced matrices

####P14


# Seurat_object_P14_CD31_WT_Rep1

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P14_S129_WT_CD31/Sample_Gael_190812_1200st_P14_CD31_SIRT3KO_Retina_1200st_N706")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))

dim(res_mat_unspliced)

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)

Seurat_object_P14_CD31_WT_Rep1 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P14_CD31_WT_Rep1)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

Seurat_object_P14_CD31_WT_Rep1[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P14_CD31_WT_Rep1[["Rep"]] <- "rep1"
Seurat_object_P14_CD31_WT_Rep1[["TimePoint"]] <- "P14"
Seurat_object_P14_CD31_WT_Rep1[["Genotype"]] <- "WT"
Seurat_object_P14_CD31_WT_Rep1[["Condition"]] <- "OIR"
Seurat_object_P14_CD31_WT_Rep1[["Sample"]] <- "P14_OIR_WT"

summary(Seurat_object_P14_CD31_WT_Rep1@meta.data)


# Seurat_object_P14_CD31_WT_Rep2

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P14_S129_WT_CD31/Sample_Gael_190822_1400st_P14_S129_OIR_CD31_N701")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))
dim(res_mat_unspliced)

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)

Seurat_object_P14_CD31_WT_Rep2 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P14_CD31_WT_Rep2)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

Seurat_object_P14_CD31_WT_Rep2[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P14_CD31_WT_Rep2[["Rep"]] <- "rep2"
Seurat_object_P14_CD31_WT_Rep2[["TimePoint"]] <- "P14"
Seurat_object_P14_CD31_WT_Rep2[["Genotype"]] <- "WT"
Seurat_object_P14_CD31_WT_Rep2[["Condition"]] <- "OIR"
Seurat_object_P14_CD31_WT_Rep2[["Sample"]] <- "P14_OIR_WT"

summary(Seurat_object_P14_CD31_WT_Rep2@meta.data)

# Seurat_object_P14_CD31_SIRT3KO_Rep1

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P14_S129_Sirt3KO_CD31/Sample_Gael_20190830_2000st_SIRT3_CD31_P14_n703")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))
dim(res_mat_unspliced)

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)



res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)

res_mat_spliced <- res_mat_spliced[,colnames(res_mat_unspliced)]

Seurat_object_P14_CD31_SIRT3KO_Rep1 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P14_CD31_SIRT3KO_Rep1)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

Seurat_object_P14_CD31_SIRT3KO_Rep1[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P14_CD31_SIRT3KO_Rep1[["Rep"]] <- "rep1"
Seurat_object_P14_CD31_SIRT3KO_Rep1[["TimePoint"]] <- "P14"
Seurat_object_P14_CD31_SIRT3KO_Rep1[["Genotype"]] <- "SIRT3KO"
Seurat_object_P14_CD31_SIRT3KO_Rep1[["Condition"]] <- "OIR"
Seurat_object_P14_CD31_SIRT3KO_Rep1[["Sample"]] <- "P14_OIR_SIRT3KO"

summary(Seurat_object_P14_CD31_SIRT3KO_Rep1@meta.data)



# Seurat_object_P14_CD31_SIRT3KO_Rep2

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P14_S129_Sirt3KO_CD31/Sample_Gael_20190911_1000st_P14_OIR_SIRT3_N704")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

#bc_unspl <- readLines("./counts_unfiltered/spliced.barcodes.txt")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))
dim(res_mat_unspliced)

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)

res_mat_spliced <- res_mat_spliced[,colnames(res_mat_unspliced)]

Seurat_object_P14_CD31_SIRT3KO_Rep2 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P14_CD31_SIRT3KO_Rep2)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

Seurat_object_P14_CD31_SIRT3KO_Rep2[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P14_CD31_SIRT3KO_Rep2[["Rep"]] <- "rep2"
Seurat_object_P14_CD31_SIRT3KO_Rep2[["TimePoint"]] <- "P14"
Seurat_object_P14_CD31_SIRT3KO_Rep2[["Genotype"]] <- "SIRT3KO"
Seurat_object_P14_CD31_SIRT3KO_Rep2[["Condition"]] <- "OIR"
Seurat_object_P14_CD31_SIRT3KO_Rep2[["Sample"]] <- "P14_OIR_SIRT3KO"

summary(Seurat_object_P14_CD31_SIRT3KO_Rep2@meta.data)



###P17


# Seurat_object_P17_CD31_WT_Rep1

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P17_S129_WT_CD31/Sample_Gael_181005_R4_wt_1000st_N705")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))
dim(res_mat_unspliced)

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)

Seurat_object_P17_CD31_WT_Rep1 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P17_CD31_WT_Rep1)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

Seurat_object_P17_CD31_WT_Rep1[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P17_CD31_WT_Rep1[["Rep"]] <- "rep1"
Seurat_object_P17_CD31_WT_Rep1[["TimePoint"]] <- "P17"
Seurat_object_P17_CD31_WT_Rep1[["Genotype"]] <- "WT"
Seurat_object_P17_CD31_WT_Rep1[["Condition"]] <- "OIR"
Seurat_object_P17_CD31_WT_Rep1[["Sample"]] <- "P17_OIR_WT"

summary(Seurat_object_P17_CD31_WT_Rep1@meta.data)


# Seurat_object_P17_CD31_WT_Rep2

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P17_S129_WT_CD31/Sample_Gael_181115_R1_S129_CD31pos_P17_OIR_2000ST_N701")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))
dim(res_mat_unspliced)

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)

Seurat_object_P17_CD31_WT_Rep2 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P17_CD31_WT_Rep2)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

Seurat_object_P17_CD31_WT_Rep2[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P17_CD31_WT_Rep2[["Rep"]] <- "rep2"
Seurat_object_P17_CD31_WT_Rep2[["TimePoint"]] <- "P17"
Seurat_object_P17_CD31_WT_Rep2[["Genotype"]] <- "WT"
Seurat_object_P17_CD31_WT_Rep2[["Condition"]] <- "OIR"
Seurat_object_P17_CD31_WT_Rep2[["Sample"]] <- "P17_OIR_WT"

summary(Seurat_object_P17_CD31_WT_Rep2@meta.data)




# Seurat_object_P17_CD31_WT_Rep3

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P17_S129_WT_CD31/Sample_Gael_R1_cd31pos_3000STAMPS_12032018_N701")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))
dim(res_mat_unspliced)

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)

Seurat_object_P17_CD31_WT_Rep3 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P17_CD31_WT_Rep3)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)


Seurat_object_P17_CD31_WT_Rep3[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P17_CD31_WT_Rep3[["Rep"]] <- "rep3"
Seurat_object_P17_CD31_WT_Rep3[["TimePoint"]] <- "P17"
Seurat_object_P17_CD31_WT_Rep3[["Genotype"]] <- "WT"
Seurat_object_P17_CD31_WT_Rep3[["Condition"]] <- "OIR"
Seurat_object_P17_CD31_WT_Rep3[["Sample"]] <- "P17_OIR_WT"

summary(Seurat_object_P17_CD31_WT_Rep3@meta.data)




# Seurat_object_P17_CD31_SIRT3KO_Rep1

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P17_S129_Sirt3KO_CD31/Sample_Gael_181004_R3_ko_1000st_N704")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))
dim(res_mat_unspliced)

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)

Seurat_object_P17_CD31_SIRT3KO_Rep1 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P17_CD31_SIRT3KO_Rep1)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)
Seurat_object_P17_CD31_SIRT3KO_Rep1[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P17_CD31_SIRT3KO_Rep1[["Rep"]] <- "rep1"
Seurat_object_P17_CD31_SIRT3KO_Rep1[["TimePoint"]] <- "P17"
Seurat_object_P17_CD31_SIRT3KO_Rep1[["Genotype"]] <- "SIRT3KO"
Seurat_object_P17_CD31_SIRT3KO_Rep1[["Condition"]] <- "OIR"
Seurat_object_P17_CD31_SIRT3KO_Rep1[["Sample"]] <- "P17_OIR_SIRT3KO"

summary(Seurat_object_P17_CD31_SIRT3KO_Rep1@meta.data)


# Seurat_object_P17_CD31_SIRT3KO_Rep2

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P17_S129_Sirt3KO_CD31/Sample_Gael_181114_R1_CD31pos_P17_OIR_2000ST_N702")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))
dim(res_mat_unspliced)

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)


Seurat_object_P17_CD31_SIRT3KO_Rep2 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P17_CD31_SIRT3KO_Rep2)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)

Seurat_object_P17_CD31_SIRT3KO_Rep2[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P17_CD31_SIRT3KO_Rep2[["Rep"]] <- "rep2"
Seurat_object_P17_CD31_SIRT3KO_Rep2[["TimePoint"]] <- "P17"
Seurat_object_P17_CD31_SIRT3KO_Rep2[["Genotype"]] <- "SIRT3KO"
Seurat_object_P17_CD31_SIRT3KO_Rep2[["Condition"]] <- "OIR"
Seurat_object_P17_CD31_SIRT3KO_Rep2[["Sample"]] <- "P17_OIR_SIRT3KO"

summary(Seurat_object_P17_CD31_SIRT3KO_Rep2@meta.data)



# Seurat_object_P17_CD31_SIRT3KO_Rep3

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Kb_python_Velocity/OIR_P17_S129_Sirt3KO_CD31/Sample_Gael_R2_cd31pos_3000STAMPS_14032018_N702")
list.files(".", recursive = TRUE)
list.files(".", recursive = FALSE)

res_mat_unspliced <- read_count_output("./counts_unfiltered", name = "unspliced")

res_mat_unspliced <- res_mat_unspliced[Matrix::rowSums(res_mat_unspliced) > 0,] ##genes with at least one count
dim(res_mat_unspliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_unspliced) <- tr2g$gene_name[match(rownames(res_mat_unspliced), tr2g$gene)]
rownames(res_mat_unspliced) <- toupper(rownames(res_mat_unspliced))
dim(res_mat_unspliced)


res_mat_spliced <- read_count_output("./counts_unfiltered", name = "spliced")

res_mat_spliced <- res_mat_spliced[Matrix::rowSums(res_mat_spliced) > 0,] ##genes with at least one count
dim(res_mat_spliced)

# Convert from Ensembl gene ID to gene symbol

rownames(res_mat_spliced) <- tr2g$gene_name[match(rownames(res_mat_spliced), tr2g$gene)]
rownames(res_mat_spliced) <- toupper(rownames(res_mat_spliced))
dim(res_mat_spliced)
dim(res_mat_unspliced)

res_mat_spliced <- res_mat_spliced[,colnames(res_mat_unspliced)]

Seurat_object_P17_CD31_SIRT3KO_Rep3 <- CreateSeuratObject(counts = res_mat_spliced,
  min.cells = 3, min.features = 100)

res_mat_unspliced <- res_mat_unspliced[,colnames(Seurat_object_P17_CD31_SIRT3KO_Rep3)]

res_mat_unspliced_assay <- CreateAssayObject(counts = res_mat_unspliced)


Seurat_object_P17_CD31_SIRT3KO_Rep3[["unspliced"]] <- res_mat_unspliced_assay

Seurat_object_P17_CD31_SIRT3KO_Rep3[["Rep"]] <- "rep3"
Seurat_object_P17_CD31_SIRT3KO_Rep3[["TimePoint"]] <- "P17"
Seurat_object_P17_CD31_SIRT3KO_Rep3[["Genotype"]] <- "SIRT3KO"
Seurat_object_P17_CD31_SIRT3KO_Rep3[["Condition"]] <- "OIR"
Seurat_object_P17_CD31_SIRT3KO_Rep3[["Sample"]] <- "P17_OIR_SIRT3KO"

summary(Seurat_object_P17_CD31_SIRT3KO_Rep3@meta.data)



Seurat_object_P14_CD31 <- merge(x=Seurat_object_P14_CD31_WT_Rep1 , y=c(Seurat_object_P14_CD31_WT_Rep2, Seurat_object_P14_CD31_SIRT3KO_Rep1, Seurat_object_P14_CD31_SIRT3KO_Rep2))

Seurat_object_P17_CD31 <- merge(x=Seurat_object_P17_CD31_WT_Rep1 , y=c(Seurat_object_P17_CD31_WT_Rep2, Seurat_object_P17_CD31_WT_Rep3, Seurat_object_P17_CD31_SIRT3KO_Rep1, Seurat_object_P17_CD31_SIRT3KO_Rep2, Seurat_object_P17_CD31_SIRT3KO_Rep3))



Seurat_object <- merge(x=Seurat_object_P14_CD31 , y=c(Seurat_object_P17_CD31))



Seurat_object[["percent.mt"]] <- PercentageFeatureSet(Seurat_object, pattern = "^MT-")


CRYA.genes <- grep(pattern = "^CRYA", x = rownames(x = Seurat_object[['RNA']]), value = TRUE)
CRYB.genes <- grep(pattern = "^CRYB", x = rownames(x = Seurat_object[['RNA']]), value = TRUE)
CRYG.genes <- grep(pattern = "^CRYG", x = rownames(x = Seurat_object[['RNA']]), value = TRUE)
crystal.genes <- c(CRYA.genes, CRYB.genes)
percent.crystal <- Matrix::colSums(Seurat_object[['RNA']][crystal.genes, ])/Matrix::colSums(Seurat_object[['RNA']])

Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.crystal, col.name = "percent.crystal")

summary(Seurat_object@meta.data)



setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/Integration")



png(filename="VlnPlotQC.Sample.Init.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.crystal"), group.by = "Sample", pt.size=-1)
dev.off()



Seurat_object <- subset(Seurat_object, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & nCount_RNA < 10000 & percent.mt < 10 & percent.crystal < 0.01)



png(filename="VlnPlotQC.Sample.filter.png", width=1500, height=1000, bg = "white", res = 50)
VlnPlot(object = Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.crystal"), group.by = "Sample", pt.size=-1)
dev.off()




#Cell Cycle scoring

cc.genes <- readLines(con="/home/gaelcge/projects/def-jsjoyal/gaelcge/Seurat_ressource/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]


Seurat_object <- CellCycleScoring(Seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)


Seurat_object <- NormalizeData(Seurat_object, verbose = FALSE)

Seurat_object <- ScaleData(Seurat_object, verbose = FALSE, vars.to.regress = c("nFeature_RNA", "percent.mt", "Sample", "S.Score", "G2M.Score"))



#GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object, i.e. columns in object@data.info, PC scores etc.
#Since there is a rare subset of cells with an outlier level of high mitochondrial percentage, and also low UMI content, we filter these as well

plot1 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="percentMitoPlot.nGenePlot.png", width=1300, height=600, bg = "white", res = 150)
CombinePlots(plots = list(plot1, plot2))
dev.off()


Seurat_object <- SCTransform(Seurat_object, verbose = FALSE, vars.to.regress = c("nFeature_RNA", "percent.mt", "Sample", "S.Score", "G2M.Score"))



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

png(filename="umap.sct.splitedbySample.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Sample")
dev.off()

png(filename="umap.sct_Sample.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Sample")
dev.off()

png(filename="umap.sct_TimePoint.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "TimePoint")
dev.off()

png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


png(filename="umap.sct_Phase.integrated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Phase")
dev.off()

png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1", "ESM1", "AQP1"))
dev.off()




markers.retina.dotplot <- rev(c("LHX1", "PAX6", "SLC16A6", "VSX2", "SLC5A7", "RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "FBN1", "CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A"))


png(filename="DotPlot_markers.retina.sct.initial.png", res = 150, width=1000, height=1500)
DotPlot(
  Seurat_object,
  assay = "SCT",
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



png(filename="DotPlot_markers.retina.rna.initial.png", res = 150, width=1000, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
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




##Subset ECs

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/ECs")


Seurat_object <- subset(
  Seurat_object, ident= c(0,2,4,6,7,8,9,10,11,12,13,14,17,21,25,27,28)
  )

#Filter out contam twice

Seurat_object <- subset(
  Seurat_object, ident= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
  )

Seurat_object <- subset(
  Seurat_object, ident= c(0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17)
  )

###


DefaultAssay(Seurat_object) <- "RNA"


Seurat_object <- NormalizeData(Seurat_object, verbose = FALSE)


#Scaling the data
all.genes <- rownames(Seurat_object)

Seurat_object$CC.Difference <- Seurat_object$S.Score - Seurat_object$G2M.Score

Seurat_object <- ScaleData(Seurat_object, features = all.genes, vars.to.regress = c("nFeature_RNA", "percent.mt", "CC.Difference", "Sample"))


Seurat_object<- FindVariableFeatures(Seurat_object, selection.method = "vst", nfeatures = 2000)

top10_VG <- head(VariableFeatures(Seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10_VG, repel = TRUE)
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

png(filename="umap.sct.splitedbySample.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Sample")
dev.off()

png(filename="umap.sct_Sample.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Sample")
dev.off()

png(filename="umap.sct_TimePoint.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "TimePoint")
dev.off()

png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


png(filename="umap.sct_Phase.integrated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Phase")
dev.off()

png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "APLN", "BMX"))
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



markers.retina.dotplot <- unique(c("LHX1", "PAX6", "SLC16A6", "VSX2", "SLC5A7", "RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "FBN1", "CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A",top10_VG))


png(filename="DotPlot_markers.retina.sct.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "SCT",
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



png(filename="DotPlot_markers.retina.rna.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
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


png(filename="DotPlot_markers.Tuft.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("AQP1", "SIX3", "MAL", "COL18A1"),
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



png(filename="DotPlot_markers.Tip.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("ESM1", "ANGPT2", "APLN", "SERPINE1"),
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

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


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


saveRDS(Seurat_object, "ECs.rds")


##Subset tip and tuft

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/ECs")

Seurat_object <- readRDS("ECs.rds")

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/Tiptufts")

Seurat_object <- subset(
  Seurat_object, ident= c(1,4,5,8,9,11,15)
  )


Seurat_object <- subset(
  Seurat_object, ident= c(0,1,2,3,5,6,7,8,9,10,11)
  )

DefaultAssay(Seurat_object) <- "RNA"


Seurat_object<- FindVariableFeatures(Seurat_object, selection.method = "vst", nfeatures = 2000)

top10_VG <- head(VariableFeatures(Seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10_VG, repel = TRUE)
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

png(filename="umap.sct.splitedbySample.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Sample")
dev.off()

png(filename="umap.sct_Sample.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Sample")
dev.off()

png(filename="umap.sct_TimePoint.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "TimePoint")
dev.off()

png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


png(filename="umap.sct_Phase.integrated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Phase")
dev.off()

png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "APLN", "BMX"))
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



markers.retina.dotplot <- unique(c("LHX1", "PAX6", "SLC16A6", "VSX2", "SLC5A7", "RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "FBN1", "CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A",top10_VG))


png(filename="DotPlot_markers.retina.sct.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "SCT",
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



png(filename="DotPlot_markers.retina.rna.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
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


png(filename="DotPlot_markers.Tuft.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("AQP1", "SIX3", "MAL", "COL18A1"),
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



png(filename="DotPlot_markers.Tip.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("ESM1", "ANGPT2", "APLN", "SERPINE1"),
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

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


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


##Annotate

Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")

cell_type_assigned <- c("Tip",      # cluster 0
                      "Tip",     # cluster 1
                      "Tip",     # cluster 2
                      "Tip",     # cluster 3
                      "Tip",     # cluster 4
                      "Tip",     # cluster 5
                      "Tip",     # cluster 6
                      "Tip",     # cluster 7
                      "Tip",     # cluster 8
                      "Tuft",     # cluster 9
                      "Tip",     # cluster 10
                      "Tip"     # cluster 11
                       )    



names(cell_type_assigned) <- levels(Seurat_object)

Seurat_object <- RenameIdents(Seurat_object, cell_type_assigned)

Seurat_object[["Cell_Type_TipTuft"]] <- Idents(object = Seurat_object)

png(filename="umap.sct_filtered.annotated.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()



library(pals)
tab = prop.table(table(Seurat_object$Cell_Type_TipTuft, Seurat_object$Sample), margin = 2)*100
pdf(paste("prop.table","Cond_TimePoint","plot.pdf", sep="."), width=5, height=6)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
dev.off()



#Backup for reusing the object

Seurat_object_bup <- Seurat_object

##Subset WT tip and tuft
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/TiptuftsWTSampled")

Seurat_object <- Seurat_object_bup

Idents(Seurat_object) <- "Genotype"

Seurat_object <- subset(
  Seurat_object, ident="WT")


DefaultAssay(Seurat_object) <- "RNA"


Seurat_object<- FindVariableFeatures(Seurat_object, selection.method = "vst", nfeatures = 2000)

top10_VG <- head(VariableFeatures(Seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10_VG, repel = TRUE)
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


Seurat_object <- RunUMAP(Seurat_object, dims = 1:10,
                            n.components = 2L)


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

png(filename="umap.sct.splitedbySample.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Sample")
dev.off()

png(filename="umap.sct_Sample.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Sample")
dev.off()

png(filename="umap.sct_TimePoint.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "TimePoint")
dev.off()

png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


png(filename="umap.sct_Phase.integrated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Phase")
dev.off()

png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "APLN", "BMX"))
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



markers.retina.dotplot <- unique(c("LHX1", "PAX6", "SLC16A6", "VSX2", "SLC5A7", "RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "FBN1", "CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A",top10_VG))


png(filename="DotPlot_markers.retina.sct.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "SCT",
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



png(filename="DotPlot_markers.retina.rna.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
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


png(filename="DotPlot_markers.Tuft.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("AQP1", "SIX3", "MAL", "COL18A1"),
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



png(filename="DotPlot_markers.Tip.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("ESM1", "ANGPT2", "APLN", "SERPINE1"),
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

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


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


##Run Velocity
library(velocyto.R)
library(SeuratWrappers)


Seurat_object <- RunVelocity(object = Seurat_object, spliced ="RNA", unspliced = "unspliced", 
  deltaT = 1, kCells = 25, fit.quantile = 0.02, ncores = availableCores()-1)

# RunVelocity <- function(
#   object,
#   spliced = 'spliced',
#   unspliced = 'unspliced',
#   ambiguous = NULL,
#   spliced.average = 0.2,
#   unspliced.average = 0.05,
#   reduction = 'pca',
#   group.by = 'ident',
#   cells = NULL,
#   graph = NULL,
#   ncores = 1,
#   verbose = TRUE




ident.colors <- (scales::hue_pal())(n = length(x = levels(x = Seurat_object)))
names(x = ident.colors) <- levels(x = Seurat_object)
cell.colors <- ident.colors[Idents(object = Seurat_object)]
names(x = cell.colors) <- colnames(x = Seurat_object)


png("show.velocity.on.embedding.cor_Seurat_object.png", height=1000, width=1000)
show.velocity.on.embedding.cor(
  emb = Embeddings(object = Seurat_object, reduction = "umap"), 
  vel = Tool(object = Seurat_object, 
  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
  cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, 
  arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1, n.cores = availableCores()-1)
dev.off()



##Subset SIRT3KO tip and tuft
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/TiptuftsKO")

Seurat_object <- Seurat_object_bup

Idents(Seurat_object) <- "Genotype"

Seurat_object <- subset(
  Seurat_object, ident="SIRT3KO")


DefaultAssay(Seurat_object) <- "RNA"


Seurat_object<- FindVariableFeatures(Seurat_object, selection.method = "vst", nfeatures = 2000)

top10_VG <- head(VariableFeatures(Seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10_VG, repel = TRUE)
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


Seurat_object <- RunUMAP(Seurat_object, dims = 1:10,
                            n.components = 2L)


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

png(filename="umap.sct.splitedbySample.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Sample")
dev.off()

png(filename="umap.sct_Sample.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Sample")
dev.off()

png(filename="umap.sct_TimePoint.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "TimePoint")
dev.off()

png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


png(filename="umap.sct_Phase.integrated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Phase")
dev.off()

png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "APLN", "BMX"))
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



markers.retina.dotplot <- unique(c("LHX1", "PAX6", "SLC16A6", "VSX2", "SLC5A7", "RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "FBN1", "CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A",top10_VG))


png(filename="DotPlot_markers.retina.sct.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "SCT",
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



png(filename="DotPlot_markers.retina.rna.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
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


png(filename="DotPlot_markers.Tuft.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("AQP1", "SIX3", "MAL", "COL18A1"),
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



png(filename="DotPlot_markers.Tip.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("ESM1", "ANGPT2", "APLN", "SERPINE1"),
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

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


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


##Run Velocity
library(velocyto.R)
library(SeuratWrappers)


Seurat_object <- RunVelocity(object = Seurat_object, spliced ="RNA", unspliced = "unspliced", 
  deltaT = 1, kCells = 25, fit.quantile = 0.02, ncores = availableCores()-1)

# RunVelocity <- function(
#   object,
#   spliced = 'spliced',
#   unspliced = 'unspliced',
#   ambiguous = NULL,
#   spliced.average = 0.2,
#   unspliced.average = 0.05,
#   reduction = 'pca',
#   group.by = 'ident',
#   cells = NULL,
#   graph = NULL,
#   ncores = 1,
#   verbose = TRUE




ident.colors <- (scales::hue_pal())(n = length(x = levels(x = Seurat_object)))
names(x = ident.colors) <- levels(x = Seurat_object)
cell.colors <- ident.colors[Idents(object = Seurat_object)]
names(x = cell.colors) <- colnames(x = Seurat_object)


png("show.velocity.on.embedding.cor_Seurat_object.png", height=1000, width=1000)
show.velocity.on.embedding.cor(
  emb = Embeddings(object = Seurat_object, reduction = "umap"), 
  vel = Tool(object = Seurat_object, 
  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
  cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, 
  arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1, n.cores = availableCores()-1)
dev.off()





##SUbset tip and tuft but at same cell counts between sample

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/ECs")

Seurat_object <- readRDS("ECs.rds")

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/TiptuftsSampled")

Seurat_object <- subset(
  Seurat_object, ident= c(1,4,5,8,9,11)
  )

#Down sample per sample

Idents(Seurat_object) <- "Sample"

Seurat_object <- subset(
  Seurat_object, downsample = 500)

##REmove unwanted clusters

Seurat_object <- subset(
  Seurat_object, ident= c(0,1,2,3,5,6,7,8,9,10,11)
  )

DefaultAssay(Seurat_object) <- "RNA"


Seurat_object<- FindVariableFeatures(Seurat_object, selection.method = "vst", nfeatures = 2000)

top10_VG <- head(VariableFeatures(Seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10_VG, repel = TRUE)
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


Seurat_object <- RunUMAP(Seurat_object, dims = 1:15,
                            n.components = 2L,
                            n.neighbors = 30L, # 30L /n.neighbors: This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
                            min.dist = 0.2, #  0.3 / This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are moreevenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
                            spread = 5 # 1 / The effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
          )


##Cluster the cells

Seurat_object <- FindNeighbors(Seurat_object, 
                dims = 1:2, 
                reduction = "umap"
                )

Seurat_object <- FindClusters(Seurat_object, 
                resolution = 0.8, 
                reduction = "umap"
                )

png(filename="umap.sct.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE)
dev.off()

png("ECSubtypeFeaturePlot.png", width=1500, height =1500)
FeaturePlot(Seurat_object, features=c("AQP1", "BMX", "GJA4", "FBLN2", "PTN", "MGP", "CXCL12", "SLC6A6", "GLUL", "APOD", "FOS", "FOSB", "ATF3", "BTG2", "APOLD1", "EGR1", "APLN", "SERPINE1", "ESM1", "ANGPT2", "PRSS23", "TRP53I11", "TOP2A", "MKI67",
"BIRC5"), cols=c("blue", "red")) + RotatedAxis()
dev.off()


png(filename="umap.sct.splitedbySample.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Sample")
dev.off()

png(filename="umap.sct_Sample.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Sample")
dev.off()

png(filename="umap.sct_TimePoint.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "TimePoint")
dev.off()

png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


png(filename="umap.sct_Phase.integrated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Phase")
dev.off()

png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "APLN", "BMX"))
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


markers.retina.dotplot <- unique(c("LHX1", "PAX6", "SLC16A6", "VSX2", "SLC5A7", "RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "FBN1", "CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A",top10_VG))


png(filename="DotPlot_markers.retina.sct.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "SCT",
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



png(filename="DotPlot_markers.retina.rna.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
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


png(filename="DotPlot_markers.Tuft.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("AQP1", "SIX3", "MAL", "COL18A1"),
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



png(filename="DotPlot_markers.Tip.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("ESM1", "ANGPT2", "APLN", "SERPINE1"),
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

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


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


##Annotate

Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")

cell_type_assigned <- c("Tip",      # cluster 0
                      "Tip",     # cluster 1
                      "Tip",     # cluster 2
                      "Tip",     # cluster 3
                      "Tip",     # cluster 4
                      "Tip",     # cluster 5
                      "Tuft",     # cluster 6
                      "Tip",     # cluster 7
                      "Tip",     # cluster 8
                      "Tip",     # cluster 9
                      "Tip",     # cluster 10
                      "Tip",     # cluster 11
                      "Tuft",     # cluster 12
                      "Tip"     # cluster 13
                       )    



names(cell_type_assigned) <- levels(Seurat_object)

Seurat_object <- RenameIdents(Seurat_object, cell_type_assigned)

Seurat_object[["Cell_Type_TipTuft"]] <- Idents(object = Seurat_object)

png(filename="umap.sct_filtered.annotated.png", width=1500, height=1500, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()



library(pals)
tab = prop.table(table(Seurat_object$Cell_Type_TipTuft, Seurat_object$Sample), margin = 2)*100
pdf(paste("prop.table","Cond_TimePoint","plot.pdf", sep="."), width=5, height=6)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
dev.off()



library(pals)
tab = prop.table(table(Seurat_object$Phase, Seurat_object$Cell_Type_TipTuft), margin = 2)*100
pdf(paste("prop.table","Cell_Type_TipTuft_Phase","plot.pdf", sep="."), width=5, height=6)
ggplot(as.data.frame(tab),aes(x=Var2,y=Freq,fill=Var1)) + 
geom_col() +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
dev.off()


saveRDS(Seurat_object, "TiptuftsSampled.rds")


setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/TiptuftsSampled")

Seurat_object <- readRDS("TiptuftsSampled.rds")



##Run Velocity
library(velocyto.R)
library(SeuratWrappers)


Seurat_object <- RunVelocity(object = Seurat_object, spliced ="RNA", unspliced = "unspliced", 
  deltaT = 1, kCells = 25, fit.quantile = 0.02, ncores = availableCores()-1)

# RunVelocity <- function(
#   object,
#   spliced = 'spliced',
#   unspliced = 'unspliced',
#   ambiguous = NULL,
#   spliced.average = 0.2,
#   unspliced.average = 0.05,
#   reduction = 'pca',
#   group.by = 'ident',
#   cells = NULL,
#   graph = NULL,
#   ncores = 1,
#   verbose = TRUE




ident.colors <- (scales::hue_pal())(n = length(x = levels(x = Seurat_object)))
names(x = ident.colors) <- levels(x = Seurat_object)
cell.colors <- ident.colors[Idents(object = Seurat_object)]
names(x = cell.colors) <- colnames(x = Seurat_object)


png("show.velocity.on.embedding.cor_Seurat_object.png", height=1000, width=1000)
show.velocity.on.embedding.cor(
  emb = Embeddings(object = Seurat_object, reduction = "umap"), 
  vel = Tool(object = Seurat_object, 
  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
  cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, 
  arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1, n.cores = availableCores()-1)
dev.off()



#Backup for reusing the object

Seurat_object_bup <- Seurat_object

##Subset WT tip and tuft
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/TiptuftsWTSampled")

Seurat_object <- Seurat_object_bup

Idents(Seurat_object) <- "Genotype"

Seurat_object <- subset(
  Seurat_object, ident="WT")


DefaultAssay(Seurat_object) <- "RNA"

Idents(Seurat_object) <- "Cell_Type_TipTuft"

##Run Velocity on non recluster cells

library(velocyto.R)
library(SeuratWrappers)



Seurat_object <- RunVelocity(object = Seurat_object, spliced ="RNA", unspliced = "unspliced", 
  reduction = 'umap', spliced.average = 0.2 , unspliced.average = 0.05, 
  ncores = availableCores()-1)

# RunVelocity <- function(
#   object,
#   spliced = 'spliced',
#   unspliced = 'unspliced',
#   ambiguous = NULL,
#   spliced.average = 0.2,
#   unspliced.average = 0.05,
#   reduction = 'pca',
#   group.by = 'ident',
#   cells = NULL,
#   graph = NULL,
#   ncores = 1,
#   verbose = TRUE




ident.colors <- (scales::hue_pal())(n = length(x = levels(x = Seurat_object)))
names(x = ident.colors) <- levels(x = Seurat_object)
cell.colors <- ident.colors[Idents(object = Seurat_object)]
names(x = cell.colors) <- colnames(x = Seurat_object)


pdf("show.velocity.on.embedding.cor_Seurat_object_WT_Native.pdf", height=5, width=5)
show.velocity.on.embedding.cor(
  emb = Embeddings(object = Seurat_object, reduction = "umap"), 
  vel = Tool(object = Seurat_object, 
  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 1), 
  cex = 1.5, arrow.scale = 8, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, 
  arrow.lwd = 0.8, do.par = FALSE, cell.border.alpha = 0.01, n.cores = availableCores()-1)
dev.off()


## Re cluster
Seurat_object<- FindVariableFeatures(Seurat_object, selection.method = "vst", nfeatures = 2000)

top10_VG <- head(VariableFeatures(Seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10_VG, repel = TRUE)
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


Seurat_object <- RunUMAP(Seurat_object, dims = 1:10,
                            n.components = 2L)


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

png(filename="umap.sct.splitedbySample.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Sample")
dev.off()

png(filename="umap.sct_Sample.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Sample")
dev.off()

png(filename="umap.sct_TimePoint.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "TimePoint")
dev.off()

png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


png(filename="umap.sct_Phase.integrated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Phase")
dev.off()

png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "APLN", "BMX"))
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



markers.retina.dotplot <- unique(c("LHX1", "PAX6", "SLC16A6", "VSX2", "SLC5A7", "RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "FBN1", "CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A",top10_VG))


png(filename="DotPlot_markers.retina.sct.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "SCT",
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



png(filename="DotPlot_markers.retina.rna.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
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


png(filename="DotPlot_markers.Tuft.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("AQP1", "SIX3", "MAL", "COL18A1"),
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



png(filename="DotPlot_markers.Tip.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("ESM1", "ANGPT2", "APLN", "SERPINE1"),
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

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


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



Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")

cell_type_assigned <- c("Tip",      # cluster 0
                      "Tuft",     # cluster 1
                      "Tip",     # cluster 2
                      "Tip",     # cluster 3
                      "Tip",     # cluster 4
                      "Tip",     # cluster 5
                      "Tip",     # cluster 6
                      "Tip",     # cluster 7
                      "Tip"     # cluster 8
                       )    



names(cell_type_assigned) <- levels(Seurat_object)

Seurat_object <- RenameIdents(Seurat_object, cell_type_assigned)

Seurat_object[["Cell_Type_TipTuft"]] <- Idents(object = Seurat_object)

Idents(Seurat_object) <- "Cell_Type_TipTuft"



##Run Velocity
library(velocyto.R)
library(SeuratWrappers)



Seurat_object <- RunVelocity(object = Seurat_object, spliced ="RNA", unspliced = "unspliced", 
  reduction = 'umap', spliced.average = 0.2 , unspliced.average = 0.05, 
  ncores = availableCores()-1)

# RunVelocity <- function(
#   object,
#   spliced = 'spliced',
#   unspliced = 'unspliced',
#   ambiguous = NULL,
#   spliced.average = 0.2,
#   unspliced.average = 0.05,
#   reduction = 'pca',
#   group.by = 'ident',
#   cells = NULL,
#   graph = NULL,
#   ncores = 1,
#   verbose = TRUE




ident.colors <- (scales::hue_pal())(n = length(x = levels(x = Seurat_object)))
names(x = ident.colors) <- levels(x = Seurat_object)
cell.colors <- ident.colors[Idents(object = Seurat_object)]
names(x = cell.colors) <- colnames(x = Seurat_object)


pdf("show.velocity.on.embedding.cor_Seurat_object_WT_reclustered.pdf", height=5, width=5)
show.velocity.on.embedding.cor(
  emb = Embeddings(object = Seurat_object, reduction = "umap"), 
  vel = Tool(object = Seurat_object, 
  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 1), 
  cex = 1.5, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, 
  arrow.lwd = 0.8, do.par = FALSE, cell.border.alpha = 0.01, n.cores = availableCores()-1)
dev.off()




##Subset SIRT3KO tip and tuft
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Sirt3_Velocity/TiptuftsKOSampled")

Seurat_object <- Seurat_object_bup

Idents(Seurat_object) <- "Genotype"

Seurat_object <- subset(
  Seurat_object, ident="SIRT3KO")


DefaultAssay(Seurat_object) <- "RNA"


Idents(Seurat_object) <- "Cell_Type_TipTuft"

##Run Velocity on non recluster cells

library(velocyto.R)
library(SeuratWrappers)


Seurat_object <- RunVelocity(object = Seurat_object, spliced ="RNA", unspliced = "unspliced", 
  reduction = 'umap', spliced.average = 0.2 , unspliced.average = 0.05, 
  ncores = availableCores()-1)

# RunVelocity <- function(
#   object,
#   spliced = 'spliced',
#   unspliced = 'unspliced',
#   ambiguous = NULL,
#   spliced.average = 0.2,
#   unspliced.average = 0.05,
#   reduction = 'pca',
#   group.by = 'ident',
#   cells = NULL,
#   graph = NULL,
#   ncores = 1,
#   verbose = TRUE




ident.colors <- (scales::hue_pal())(n = length(x = levels(x = Seurat_object)))
names(x = ident.colors) <- levels(x = Seurat_object)
cell.colors <- ident.colors[Idents(object = Seurat_object)]
names(x = cell.colors) <- colnames(x = Seurat_object)


pdf("show.velocity.on.embedding.cor_Seurat_object_KO_Native.pdf", height=5, width=5)
show.velocity.on.embedding.cor(
  emb = Embeddings(object = Seurat_object, reduction = "umap"), 
  vel = Tool(object = Seurat_object, 
  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 1), 
  cex = 1.5, arrow.scale = 8, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, 
  arrow.lwd = 0.8, do.par = FALSE, cell.border.alpha = 0.01, n.cores = availableCores()-1)
dev.off()


##Recluster

Seurat_object<- FindVariableFeatures(Seurat_object, selection.method = "vst", nfeatures = 2000)

top10_VG <- head(VariableFeatures(Seurat_object), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10_VG, repel = TRUE)
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


Seurat_object <- RunUMAP(Seurat_object, dims = 1:10,
                            n.components = 2L)


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

png(filename="umap.sct.splitedbySample.png", width=3500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by = "Sample")
dev.off()

png(filename="umap.sct_Sample.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Sample")
dev.off()

png(filename="umap.sct_TimePoint.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "TimePoint")
dev.off()

png(filename="FeaturePlot_sct.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()


png(filename="umap.sct_Phase.integrated.png", width=1500, height=900, bg = "white", res = 150)
DimPlot(Seurat_object, reduction = "umap", label=FALSE, group.by = "Phase")
dev.off()

png(filename="FeaturePlot_sct_ECsmarkers.png", width=1500, height=1000, bg = "white", res = 150)
FeaturePlot(Seurat_object, features = c("CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "APLN", "BMX"))
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



markers.retina.dotplot <- unique(c("LHX1", "PAX6", "SLC16A6", "VSX2", "SLC5A7", "RPE65", "RHO", "OPN1SW", "TRPM1", "SNHG11", "KCNJ8", "FBN1", "CLDN5", "PECAM1", "ESM1", "AQP1", "CXCL12", "LYZ2", "RLBP1", "GFAP", "OPTC", "TOP2A",top10_VG))


png(filename="DotPlot_markers.retina.sct.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "SCT",
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



png(filename="DotPlot_markers.retina.rna.initial.png", res = 150, width=1500, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
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


png(filename="DotPlot_markers.Tuft.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("AQP1", "SIX3", "MAL", "COL18A1"),
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



png(filename="DotPlot_markers.Tip.rna.initial.png", res = 150, width=600, height=1500)
DotPlot(
  Seurat_object,
  assay = "RNA",
  c("ESM1", "ANGPT2", "APLN", "SERPINE1"),
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

top10 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top2 <- Seurat_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


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


Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")

cell_type_assigned <- c("Tip",      # cluster 0
                      "Tip",     # cluster 1
                      "Tip",     # cluster 2
                      "Tuft",     # cluster 3
                      "Tip",     # cluster 4
                      "Tuft",     # cluster 5
                      "Tip",     # cluster 6
                      "Tip"     # cluster 7
                       )    



names(cell_type_assigned) <- levels(Seurat_object)

Seurat_object <- RenameIdents(Seurat_object, cell_type_assigned)

Seurat_object[["Cell_Type_TipTuft"]] <- Idents(object = Seurat_object)

Idents(Seurat_object) <- "Cell_Type_TipTuft"


##Run Velocity
library(velocyto.R)
library(SeuratWrappers)


Seurat_object <- RunVelocity(object = Seurat_object, spliced ="RNA", unspliced = "unspliced", 
  reduction = 'umap', spliced.average = 0.2 , unspliced.average = 0.05, 
  ncores = availableCores()-1)

# RunVelocity <- function(
#   object,
#   spliced = 'spliced',
#   unspliced = 'unspliced',
#   ambiguous = NULL,
#   spliced.average = 0.2,
#   unspliced.average = 0.05,
#   reduction = 'pca',
#   group.by = 'ident',
#   cells = NULL,
#   graph = NULL,
#   ncores = 1,
#   verbose = TRUE




ident.colors <- (scales::hue_pal())(n = length(x = levels(x = Seurat_object)))
names(x = ident.colors) <- levels(x = Seurat_object)
cell.colors <- ident.colors[Idents(object = Seurat_object)]
names(x = cell.colors) <- colnames(x = Seurat_object)


pdf("show.velocity.on.embedding.cor_Seurat_object_KO_reclustered.pdf", height=5, width=5)
show.velocity.on.embedding.cor(
  emb = Embeddings(object = Seurat_object, reduction = "umap"), 
  vel = Tool(object = Seurat_object, 
  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 1), 
  cex = 1.5, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, 
  arrow.lwd = 0.8, do.par = FALSE, cell.border.alpha = 0.01, n.cores = availableCores()-1)
dev.off()

















q("no")
