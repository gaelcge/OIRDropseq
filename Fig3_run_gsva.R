library(GSVA)
library(GSEABase)
library(parallel)
library(sigPathway)
library(Seurat)
library(limma)
library(gskb)
library(future)
library(Biobase)
library(genefilter)
library(RColorBrewer)
library(GSVAdata)
data(c2BroadSets)
library(gplots)
library(heatmap3)
library(viridis)


##Import matrix
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/Mapping")

project_name <- "RETINA_CD31_ECs_Integrated"

Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))


Seurat_object[["CelltypeDataset"]] <- paste(Seurat_object@meta.data$Cell_Type,Seurat_object@meta.data$Dataset,sep="_")
table(Seurat_object[["CelltypeDataset"]])

exp_mat <- GetAssayData(object = Seurat_object, assay = "RNA", slot = "data")

#exp_mat <- exp_mat[,1:500]

###Import database from gskb

data(mm_GO)
data(mm_metabolic)
data(mm_pathway)


# #Check pathway of interest in GSKB
# grep_pathway <- grep("GLYCOLYSIS", names(mm_GO)) 
# mm_GO[grep_pathway]


# grep_pathway <- grep("GLYCOLYSIS", names(mm_metabolic))
# mm_metabolic[grep_pathway]


# grep_pathway <- grep("GLYCOLYSIS", names(mm_pathway)) 
# mm_pathway[grep_pathway]

# ###Import database from gskb gmt files

# gmt.Metabolic <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MousePath_Metabolic_gmt.gmt", geneIdType=SymbolIdentifier())

# gmt.Co_expression <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MousePath_Co-expression_gmt.gmt", geneIdType=SymbolIdentifier())

# gmt.GO <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MousePath_GO_gmt.gmt", geneIdType=SymbolIdentifier())

# gmt.curated <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MousePath_Pathway_gmt.gmt", geneIdType=SymbolIdentifier())

# gmt.TF <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MousePath_TF_gmt.gmt", geneIdType=SymbolIdentifier())

# gmt.miRNA <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MousePath_miRNA_gmt.gmt", geneIdType=SymbolIdentifier())

# gmt.Location <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MousePath_Location_gmt.gmt", geneIdType=SymbolIdentifier())

# gmt.Other <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/MousePath_Other_gmt.gmt", geneIdType=SymbolIdentifier())

# gmt.TFBS_ChipSeq <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/TFBS_ChipSeq_multiple_sites/", geneIdType=SymbolIdentifier())

# gmt.gskb_all <- append(gmt.Metabolic, gmt.Co_expression, gmt.GO, gmt.curated, gmt.TF, gmt.miRNA, gmt.Location, gmt.Other, gmt.TFBS_ChipSeq)

# ###Import database from MsigDB

gmt.MsigDB_c2 <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/c2.cp.v7.1.symbols.gmt", geneIdType=SymbolIdentifier())

gmt.MsigDB_c5 <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/c5.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier())

gmt.MsigDB_h <- getGmt("~/projects/def-jsjoyal/gaelcge/Gene_list/Gmt.file/h.all.v7.1.symbols.gmt", geneIdType=SymbolIdentifier())

# ##Import home made dataset

# geneset <- importGeneSets("~/projects/def-jsjoyal/jhowa105/projects/Mike/MALLETTE_SENESCENCE_UP.gmx", verbose = TRUE)

# geneset <- GeneSet(geneset[[1]]$probes, geneIdType=SymbolIdentifier(), setName="MALLETTE_SENESCENCE_UP")

gsc <- GeneSetCollection(c(gmt.MsigDB_c2,gmt.MsigDB_c5,gmt.MsigDB_h))


# ###Merge GetSet Collection

# gsc_final <- GeneSetCollection(c(gsc, gmt.MsigDB))


###Selected gmt to run in GSVA

gmt <- gsc

###Run GSVA
setwd("/home/gaelcge/scratch")

parallel.sz <- availableCores()-1

parallel.sz <- 14

out.file <- "gsva.exprs_data.msigdb_May2020.csv"

print(paste0('parallel.sz = ', parallel.sz))

gsva.mat <- gsva(as.matrix(exp_mat), gmt,
				method="gsva", #"ssgsea", "zscore", "plage"),
				kcdf="Gaussian", # "Poisson", "none"),
				min.sz=1,
				max.sz=Inf,
				verbose=TRUE)

#gsva.mat <- gsva(as.matrix(exp_mat), gmt, method='gsva',kcdf='Gaussian', parallel.sz = parallel.sz, verbose = TRUE, parallel.type="SOCK")

message(gsva.mat[1:5,1:5])

setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/GSVA")

write.table(gsva.mat, file = out.file, sep = ",", row.names = TRUE, col.names = TRUE)

q("no")


