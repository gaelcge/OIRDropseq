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
library(monocle3)
library(DDRTree)
library(pheatmap)
library(reshape)
#plan("multiprocess", workers = (availableCores()-1))
#options(future.globals.maxSize = 3000 * 1024^2)

# install.packages(c("plotly",
# 					"tidyr",
# 					"ggrepel",
# 					"gtools",
# 					"data.table",
# 					"gplots",
# 					"useful",
# 					"locfit",
# 					"Matrix",
# 					"dplyr",
# 					"ggplot2",
# 					"Seurat",
# 					"Rmagic",
# 					"readr",
# 					"phateR",
# 					"viridis",
# 					"magrittr",
# 					"harmony",
# 					"RColorBrewer",
# 					"tidyverse",
# 					"future"))

#Set directory of dataset to analyse
setwd("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/Mapping_OIRonly/Subclustering/Tip_tuft_ECs/Non_integrated")

project_name <- "RETINA_CD31_ECs_Integrated_OIRonly"

Seurat_object <- readRDS(paste(project_name, "Seurat_object.rds", sep="."))



###Subset dataset

Subseted_cells <- rownames(subset(Seurat_object@meta.data, Dataset %in% c("OIR_KO", "OIR_WT")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells = Subseted_cells)
#
Subseted_cells <- rownames(subset(Seurat_object@meta.data, Dataset %in% c("OIR_WT")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells = Subseted_cells)
#
Subseted_cells <- rownames(subset(Seurat_object@meta.data, Dataset %in% c("OIR_KO")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells = Subseted_cells)
#
Subseted_cells <- rownames(subset(Seurat_object@meta.data, Condition %in% c("OIR_P17_WT", "OIR_P17_Sirt3KO")))

Seurat_object_Subset <- SubsetData(Seurat_object, cells = Subseted_cells)


###Add GO function to Seurat object

GO_score = read.table("/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/GSVA/OIRonly/gsva.exprs.MsigDB.csv", sep = ",", header = T, row.names=1, stringsAsFactors=F)

GO_score_subset = GO_score[,rownames(Seurat_object_Subset@meta.data)]

Seurat_object_Subset[["GO"]] <- CreateAssayObject(counts = GO_score_subset)

Seurat_object_Subset <- NormalizeData(Seurat_object_Subset, assay = "GO", normalization.method = "CLR")

Seurat_object_Subset <- ScaleData(Seurat_object_Subset, assay = "GO", verbose = FALSE, model.use = "negbinom")

Seurat_object_Subset@meta.data$New_Cell_SubType <- factor(Seurat_object_Subset@meta.data$New_Cell_SubType, 
                                                  levels = c("Tip U ECs", "Tip S ECs", "Tip D ECs")
                                                  )



###Re run umap

Seurat_object_Subset<- FindVariableFeatures(Seurat_object_Subset, selection.method = "vst", nfeatures = 2000)

Seurat_object_Subset <- RunPCA(Seurat_object_Subset, features = VariableFeatures(object = Seurat_object_Subset))

Seurat_object_Subset <- RunUMAP(Seurat_object_Subset, dims = 1:10,
                            n.components = 2L)




####Create moncocle3 object

dir = "/home/gaelcge/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/Mapping_OIRonly/Subclustering/Tip_tuft_ECs/Monocle3/WTonly"
setwd(dir)

### Create inputs for new_cell_data_set object
### expression matrix should be non normalized UMI counts
exp_mat <- as.matrix(GetAssayData(Seurat_object_Subset, slot = "data", assay = "RNA"))
#exp_mat <- exp_mat[,colnames(seurat_obj@data)]
cell_metadata <- Seurat_object_Subset@meta.data
gene_annotation <- as.data.frame(rownames(exp_mat), stringsAsFactors=F)
rownames(gene_annotation) <- gene_annotation[,1]
colnames(gene_annotation) <- "gene_short_name"

### Create monocle object
cds <- new_cell_data_set(as(exp_mat,"sparseMatrix"),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


cds = preprocess_cds(cds, num_dim = 2, method = "PCA",
  norm_method = c("none"))

#cds <- align_cds(cds, alignment_group = "Batch")

### Elbow plot to vizualize variance described by each PC
#png(filename="elbowplot.cca.png", units="in",height = 10, width = 10, res = 480)
#plot_pc_variance_explained(cds)
#dev.off()


### Run UMAP 
#cds <- reduce_dimension(cds, reduction_method="UMAP", cores = 14, umap.metric = "cosine", umap.min_dist = 0.3,
#  umap.n_neighbors = 30L, umap.fast_sgd = FALSE,
#  umap.nn_method = "annoy")
#cds <- reduce_dimension(cds)

cds@int_colData@listData$reducedDims@listData$UMAP <- Embeddings(object = Seurat_object_Subset, reduction= "umap")


png(filename="umap.Subtype.png", units="in",height = 3, width = 3, res = 300)
plot_cells(cds, color_cells_by="New_Cell_SubType")
dev.off()

png(filename="umap.Condition.png", units="in",height = 3, width = 3, res = 300)
plot_cells(cds, color_cells_by="Condition")
dev.off()


cds <- cluster_cells(cds, reduction_method = c("UMAP"), k = 50,
  louvain_iter = 1, partition_qval = 1, weight = FALSE,
  resolution = 0.1, random_seed = 0L, verbose = TRUE)

png(filename="umap.partition.png", units="in",height = 5, width = 5, res = 150)
plot_cells(cds, color_cells_by = "partition", group_cells_by="partition")
dev.off()


png(filename="umap.cluster.png", units="in",height = 5, width = 5, res = 150)
plot_cells(cds)
dev.off()

marker_test_res <- top_markers(cds, group_cells_by="New_Cell_SubType", reference_cells=1000, cores=12)

#marker_test_res <- top_markers(cds, reference_cells=1000, cores=12)

top_specific_markers <- marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

png(filename="dotplot_top_specific_markers.Subtype.png", units="in",height = 5, width = 5, res = 150)
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="New_Cell_SubType",
                    ordering_type="maximal_on_diag",
                    max.size=3)
dev.off()

top_specific_markers <- marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

png(filename="dotplot_top_specific_markers.Subtype_maximal_on_diag.png", units="in",height = 5, width = 8, res = 150)
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="New_Cell_SubType",
                    ordering_type="maximal_on_diag",
                    max.size=5)
dev.off()


png(filename="umap.top_specific_marker_ids.png", units="in", height = 7, width = 10, res = 150)
plot_cells(cds,
           genes=top_specific_marker_ids,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
dev.off()

### Automated annotation with Garnett

# Require that markers have at least JS specificty score > 0.5 and
# be significant in the logistic test for identifying their cell type:
garnett_markers <- marker_test_res %>%
                        filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
                        group_by(cell_group) %>%
                        top_n(5, marker_score)
# Exclude genes that are good markers for more than one cell type:
garnett_markers <- garnett_markers %>% 
                        group_by(gene_short_name) %>%
                        filter(n() == 1)


generate_garnett_marker_file(garnett_markers, file="./marker_file.txt")


library(garnett)

colData(cds)$garnett_cluster <- clusters(cds)
garnett_classifier <- train_cell_classifier(cds = cds,
                                         marker_file = "./marker_file.txt", 
                                         db=org.Hs.eg.db::org.Hs.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL",
                                         cores=8)


cds <- classify_cells(cds, garnett_classifier,
                      db = org.Hs.eg.db::org.Hs.eg.db,
                      cluster_extend = TRUE,
                      cds_gene_id_type = "SYMBOL")

png(filename="dotplot_top_specific_markers.cluster_ext_type.png", units="in",height = 5, width = 8, res = 150)
plot_cells(cds,
           group_cells_by="cluster",
           color_cells_by="cluster_ext_type")
dev.off()



#######

#########

#Learn the trajectory graph

cds <- learn_graph(cds, use_partition = TRUE, close_loop = TRUE,
 verbose = FALSE, learn_graph_control = list(minimal_branch_len=10, euclidean_distance_ratio=1, geodesic_distance_ratio=1/3))

png(filename="umap.learn_graph.Subtype.png", units="in",height = 2, width = 4, res = 300)
plot_cells(cds, 
           color_cells_by = "New_Cell_SubType", 
           label_cell_groups=FALSE, 
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
dev.off()


png(filename="umap.learn_graph_Age.png", units="in",height = 2, width = 3, res = 300)
plot_cells(cds, 
           color_cells_by = "Age", 
           label_cell_groups=FALSE, 
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
dev.off()


png(filename="umap.learn_graph.Cell_type.png", units="in",height = 2, width = 3, res = 300)
plot_cells(cds, 
           color_cells_by = "Cell_SubType", 
           label_cell_groups=FALSE, 
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
dev.off()



##If UMAP reduc done in 3D
#plot_cells_3d(cds,color_cells_by="New_Cell_SubType")



# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="P14"){
  cell_ids <- which(colData(cds)[, "Age"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

setEPS(bg = "white", family = "Times", width=3, height=2)
postscript("umap.learn_graph_ordered_pseudotime_New_Cell_SubTyperoot.eps")
plot_cells(cds, 
           color_cells_by = "pseudotime", 
           label_cell_groups=FALSE, 
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()


setEPS(bg = "white", family = "Times", width=3, height=2)
postscript("umap.learn_graph_ordered_Age_New_Cell_SubTyperoot.eps")
plot_cells(cds, 
           color_cells_by = "Age", 
           label_cell_groups=FALSE, 
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()


setEPS(bg = "white", family = "Times", width=4, height=2)
postscript("umap.learn_graph_ordered_Cell_SubType_New_Cell_SubTyperoot.eps")
plot_cells(cds, 
           color_cells_by = "New_Cell_SubType", 
           label_cell_groups=FALSE, 
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()


setEPS(bg = "white", family = "Times", width=3, height=2)
postscript("umap.learn_graph_ordered_Cell_Type_SubTyperoot.eps")
plot_cells(cds, 
           color_cells_by = "Cell_SubType", 
           label_cell_groups=FALSE, 
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()


saveRDS(cds, "cds_OIR_WT_Tip_and_tuft.rds")

# cds_3d <- reduce_dimension(cds, max_components = 3)
# cds_3d <- cluster_cells(cds_3d)
# cds_3d <- learn_graph(cds_3d)
# cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

# png(filename="umap.learn_graph_Subtype_3D.png", units="in",height = 5, width = 8, res = 150)
# plot_cells_3d(cds_3d, color_cells_by="New_Cell_SubType")
# dev.off()

###########
#######

######Regression analysis


gene_fits <- fit_models(cds, model_formula_str = "~Age")

fit_coefs <- coefficient_table(gene_fits)

emb_time_terms <- fit_coefs %>% filter(term == "AgeP17")

SignTerms <- emb_time_terms %>% filter (q_value < 0.00001) %>%
         select(gene_short_name, term, q_value, estimate) 

cds_subset <- cds[rowData(cds)$gene_short_name %in% SignTerms$gene_short_name,]

png(filename="plot_genes_violin_fit_model_Cell_SubType.png", units="in",height = 10, width = 10, res = 100)
plot_genes_violin(cds_subset, group_cells_by="Cell_SubType", ncol=5) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

evaluate_fits_cds <- evaluate_fits(gene_fits)



#######
######   GGraph-autocorrelation analysis for comparing clusters

##Using Knn

pr_graph_test_res = graph_test(cds, neighbor_graph="knn", cores=12)

pr_deg_ids = row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

##Finding modules of co-regulated genes
gene_module_df = find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)

png(filename="umap_gene_modules.knn.png", units="in", height = 6, width = 10, res = 150)
plot_cells(cds, genes=gene_module_df, show_trajectory_graph=FALSE, label_cell_groups=FALSE)
dev.off()


###

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$Cell_SubType)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))


png(filename="module_heatmap.Subtype.knn.png", units="in",height = 8, width = 5, res = 300)
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
dev.off()


png(filename="umap_modules_specific_cluster.knn.png", units="in",height = 10, width = 15, res = 150)
plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(18,13,7,21,24,1,33,9,35)), 
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           cell_size = 1)
dev.off()


## Finding genes that change as a function of pseudotime
##Using principal_graph

cds_pr_test_res = graph_test(cds, neighbor_graph="principal_graph", cores=12)

pr_deg_ids = row.names(subset(cds_pr_test_res, q_value < 0.05))

moran_genes <- rownames(head(cds_pr_test_res[order(cds_pr_test_res$morans_test_statistic,decreasing=TRUE),], 50))

png(filename="umap_Markers_moran_genes.png", units="in", height = 10, width = 10, res = 200)
plot_cells(cds, genes=moran_genes, 
           show_trajectory_graph=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 1)
dev.off()


gene_module_df <- monocle3::find_gene_modules(cds[pr_deg_ids,], resolution=1e-2, random_seed = 40L)

png(filename="umap_gene_modules.principal_graph.png", units="in", height = 6, width = 10, res = 150)
plot_cells(cds, genes=gene_module_df, show_trajectory_graph=TRUE, label_cell_groups=FALSE)
dev.off()


cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$New_Cell_SubType)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

png(filename="module_heatmap.Subtype.principal_graph.png", units="in",height = 8, width = 5, res = 300)
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
dev.off()

pal <- colorRampPalette(c("red", "blue"))


png(filename="umap_modules_specific_cluster.principal_graph.png", units="in",height = 8, width = 13, res = 300)
plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(1,2,3,4,5)), 
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           cell_size = 1) +
           scale_colour_gradient2(low="gray100", mid = "lightcyan2",high="red2")
        
dev.off()


png(filename="umap_modules_specific_nb2_cluster.principal_graph.png", units="in",height = 3, width = 6, res = 300)
plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(2,7)), 
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           cell_size = 1) +
           scale_colour_gradient2(low="gray100", mid = "lightcyan2",high="red2")
        
dev.off()


png(filename="umap_modules_specific_nb3_cluster.principal_graph.png", units="in",height = 3, width = 6, res = 300)
plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(3,7)), 
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           cell_size = 1) +
           scale_colour_gradient2(low="gray100", mid = "lightcyan2",high="red2")
        
dev.off()


png(filename="umap_modules_specific_nb4_cluster.principal_graph.png", units="in",height = 3, width = 6, res = 300)
plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(4,7)), 
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           cell_size = 1) +
           scale_colour_gradient2(low="gray100", mid = "lightcyan2",high="red2")
        
dev.off()


png(filename="umap_modules_specific_nb5_cluster.principal_graph.png", units="in",height = 3, width = 6, res = 300)
plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(5,7)), 
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           cell_size = 1) +
           scale_colour_gradient2(low="gray100", mid = "lightcyan2",high="red2")
        
dev.off()



png(filename="umap_modules_specific_nb1_cluster.principal_graph.png", units="in",height = 3, width = 6, res = 300)
plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(1,7)), 
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           cell_size = 1) +
           scale_colour_gradient2(low="gray100", mid = "lightcyan2",high="red2")
        
dev.off()





##Plot genes in module

genes_in_module_1 <- cds_pr_test_res[which(cds_pr_test_res$gene_short_name%in%as.data.frame(gene_module_df)[which(gene_module_df$module==1),1]),2:4]

mod_genes_module_1 <- rownames(genes_in_module_1[order(genes_in_module_1$morans_test_statistic,decreasing=T),])[1:20]


mod_genes_module_1_lineage_cds = cds[rowData(cds)$gene_short_name %in% mod_genes_module_1,
                      colData(cds)$Cell_SubType %in% c("Tip ECs")]

png(filename="module_1.Pseudotime.moran_genes.pseudotime.png", units="in", height = 20, width = 5, res = 300)
plot_genes_in_pseudotime(mod_genes_module_1_lineage_cds,
                         min_expr=0.5)
dev.off()

png(filename="module_1.Pseudotime.moran_genes.New_Cell_SubType.png", units="in", height = 20, width = 5, res = 300)
plot_genes_in_pseudotime(mod_genes_module_1_lineage_cds, 
                         color_cells_by="New_Cell_SubType",
                         min_expr=0.5)
dev.off()

png(filename="module_1.Pseudotime.moran_genes.Age.png", units="in", height = 20, width = 5, res = 300)
plot_genes_in_pseudotime(mod_genes_module_1_lineage_cds, 
                         color_cells_by="Age",
                         min_expr=0.5)
dev.off()


##Plot moran genes


moran_genes_lineage_cds = cds[rowData(cds)$gene_short_name %in% moran_genes,
                      colData(cds)$Cell_SubType %in% c("Tip ECs")]



png(filename="moran_genes.Pseudotime.moran_genes.pseudotime.png", units="in", height = 40, width = 3, res = 150)
plot_genes_in_pseudotime(moran_genes_lineage_cds,
                         min_expr=0.5)
dev.off()

png(filename="moran_genes.Pseudotime.moran_genes.New_Cell_SubType.png", units="in", height = 40, width = 3, res = 150)
plot_genes_in_pseudotime(moran_genes_lineage_cds, 
                         color_cells_by="New_Cell_SubType",
                         min_expr=0.5)
dev.off()

png(filename="moran_genes.Pseudotime.moran_genes.Age.png", units="in", height = 40, width = 3, res = 150)
plot_genes_in_pseudotime(moran_genes_lineage_cds, 
                         color_cells_by="Age",
                         min_expr=0.5)
dev.off()


#########




cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$Subtype)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

png(filename="module_heatmap.Subtype.png", units="in",height = 20, width = 10, res = 480)
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
dev.off()

png(filename="Umap_modules_specific_cluster.png", units="in",height = 10, width = 15, res = 150)
plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(3,9,10,23,8,14,20,6,5,16)), 
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           cell_size = 1)
dev.off()


cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=cds@clusters$UMAP$clusters)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))



png(filename="module_heatmap.cluster.png", units="in",height = 20, width = 10, res = 480)
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
dev.off()

png(filename="Umap_modules_specific_cluster_2.png", units="in",height = 5, width = 7, res = 150)
plot_cells(cds, genes=gene_module_df %>% filter(module %in% c(23)), 
           show_trajectory_graph=TRUE,
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           cell_size = 1)
dev.off()


png(filename="Umap_cluster.png", units="in",height = 5, width = 5, res = 150)
plot_cells(cds, show_trajectory_graph=TRUE)
dev.off()


ECProlif_genes = moran_genes

ECProlif_lineage_cds = cds[rowData(cds)$gene_short_name %in% ECProlif_genes,
                      clusters(cds) %in% c(12, 4, 3, 7, 11)]

png(filename="ECProlif.Pseudotime.moran_genes.pseudotime.png", units="in", height = 5, width = 5, res = 150)
plot_genes_in_pseudotime(ECProlif_lineage_cds, 
                         color_cells_by="pseudotime",
                         min_expr=0.5)
dev.off()

png(filename="ECProlif.Pseudotime.moran_genes.Subtype.png", units="in", height = 5, width = 5, res = 150)
plot_genes_in_pseudotime(ECProlif_lineage_cds, 
                         color_cells_by="Subtype",
                         min_expr=0.5)
dev.off()

whole_lineage_cds = cds[rowData(cds)$gene_short_name %in% ECProlif_genes]

png(filename="cds.Pseudotime.moran_genes.pseudotime.png", units="in", height = 5, width = 5, res = 150)
plot_genes_in_pseudotime(whole_lineage_cds, 
                         color_cells_by="pseudotime",
                         min_expr=0.5)
dev.off()

png(filename="cds.Pseudotime.moran_genes.Subtype.png", units="in", height = 5, width = 5, res = 150)
plot_genes_in_pseudotime(whole_lineage_cds, 
                         color_cells_by="Subtype",
                         min_expr=0.5)
dev.off()


Cluster2_genes = mod_genes_module_23

Cluster2_genes_lineage_cds = cds[rowData(cds)$gene_short_name %in% Cluster2_genes,
                      clusters(cds) %in% c(12, 4, 3, 7, 8, 2)]

png(filename="Cluster2.Pseudotime.moran_genes.pseudotime.png", units="in", height = 5, width = 5, res = 150)
plot_genes_in_pseudotime(Cluster2_genes_lineage_cds,
                         min_expr=0.5)
dev.off()

png(filename="Cluster2.Pseudotime.moran_genes.Cell_type.png", units="in", height = 5, width = 5, res = 150)
plot_genes_in_pseudotime(Cluster2_genes_lineage_cds, 
                         color_cells_by="Cell_type",
                         min_expr=0.5)
dev.off()

png(filename="Cluster2.Pseudotime.moran_genes.Clusters.png", units="in", height = 5, width = 5, res = 150)
plot_genes_in_pseudotime(Cluster2_genes_lineage_cds, 
                         color_cells_by="cluster",
                         min_expr=0.5)
dev.off()


###Plot gene of interest


Gene_of_inter_cds = cds[rowData(cds)$gene_short_name %in% c("AQP1", "ESM1")]



pdf("Aqp1_Esm1.Pseudotime_genes.Clusters.pdf", height = 5, width = 5)
plot_genes_in_pseudotime(Gene_of_inter_cds#, 
                         #color_cells_by="Size_Factor",
                         #min_expr=0.01
                         )
dev.off()



#### use library(ReductionWrappers)
library(ReductionWrappers)

Seurat_paga <- PAGA(Seurat_object_Subset)


q("no")
