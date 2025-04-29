library(Seurat)
library(tidyverse)
library(nichenetr)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)


#packageurl <- "https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.5.1.tar.gz" 
#install.packages(packageurl, repos=NULL, type="source")

##Set some parameters

Log2FC_var <- 0.25
pct_var <-0.4

##load nichenet data bases
setwd ("~/projects/def-jsjoyal/mchen251/NicheNet_Tutorial/Data")
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))


-----
#convert to mouse
lr_network = lr_network %>% mutate(from = toupper(convert_human_to_mouse_symbols(from)), to = toupper(convert_human_to_mouse_symbols(to))) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() %>% toupper()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() %>% toupper()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
weighted_networks_lr = weighted_networks_lr %>% mutate(from = toupper(convert_human_to_mouse_symbols(from)), to = toupper(convert_human_to_mouse_symbols(to))) %>% drop_na()

#all ligand and all receptor from network
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

-----
##load all the resources needed and baseline
setwd("~/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Subclustering/ECs/Mapping_OIRonly/Non_integrated_WT_KO_CCregression_2")

##this will be genesetoi and the receiver
##this object has all subtypes of ecs and all ecs
seuratObj_receiver <- readRDS("RETINA_CD31_ECs_Integrated_OIRonly.Seurat_object.rds")
unique(seuratObj_receiver$Cell_Type)
# [1] "ECs"
unique(seuratObj_receiver$Cell_SubType)
# [1] "Vein ECs"          "Capillary ECs"     "Tip ECs"          
# [4] "Proliferative ECs" "Arterial ECs"      "Tuft ECs" 
table(seuratObj_receiver$Condition)
# [1] OIR_P14_Sirt3KO OIR_P17_Sirt3KO OIR_P14_WT      OIR_P17_WT 

seuratObj_receiver = SetIdent(seuratObj_receiver, value = seuratObj_receiver[["Condition"]])
seuratObj_receiver= subset(seuratObj_receiver, idents = c("OIR_P14_Sirt3KO","OIR_P14_WT"))

### Define receiver cell type = "Tip ECs"

seuratObj_receiver= subset(seuratObj_receiver, subset = Cell_SubType == "Tuft ECs")

------
##this is the larger object, will be senders/micro enviro
setwd("~/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/Clustering/Mapping")
seuratObj_sender <- readRDS("RETINA_CD31_ECs_Integrated.Seurat_object.rds")
table(seuratObj_sender$Cell_Type)
 # [1] ECs                Mural cells        Immune cells       Astrocytes        
 # [5] Lens cells         Muller Glial cells RBCs               Cones             
 # [9] Bipolar cells      Amacrine cells     Rods   

table(seuratObj_sender$Dataset)
# [1] NORM_P17_WT     NORM_P14_WT     OIR_P14_WT      OIR_P17_WT     
# [5] OIR_P14_Sirt3KO OIR_P17_Sirt3KO

seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[["Dataset"]])
seuratObj_sender= subset(seuratObj_sender, idents = c("OIR_Sirt3KO","OIR_WT"))
seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[["Cell_Type"]])
seuratObj_sender= subset(seuratObj_sender, idents = c("ECs","Lens cells","RBCs","Amacrine cells","Cones","Rods","Immune cells","Bipolar cells"), invert=T)
table(seuratObj_sender$Cell_Type)
#[1] Mural cells        Muller Glial cells           


##Set wd

setwd("~/projects/def-jsjoyal/gaelcge/Retina/Retina_Mix/Test14_CD31_NORM_OIR_P14_P17_Sirt3/Aligned/SeuratV3/Aligned_Condition/Integrated/NicheNet/TuftCells/WT_SenderP14_P17")


png(filename="umap.P14_Senders.png", width=2500, height=900, bg = "white", res = 150)
DimPlot(seuratObj_sender, group.by="Cell_Type", split.by="Dataset", , reduction = "umap", label=TRUE)
dev.off()

png(filename="umap.P14_Receivers.png", width=2500, height=900, bg = "white", res = 150)
DimPlot(seuratObj_receiver, group.by="New_Cell_SubType", split.by="Condition", reduction = "umap", label=TRUE)
dev.off()
-----
###reciver
###DE between tuft and tip is the genesetoi
###find genesetoi in the larger ec group 



###sirt3ko effect on larger EC group

seuratObj_receiver = SetIdent(seuratObj_receiver, value = seuratObj_receiver[["Cell_Type"]])
receiver = "ECs"
##will be used to define the background, need background genes for function predict_ligand_activities
expressed_genes_receiver = get_expressed_genes(receiver, seuratObj_receiver, pct = 0.10,assay_oi="RNA")

####DE between condition of interest

seuratObj_receiver = SetIdent(seuratObj_receiver, value = seuratObj_receiver[["Cell_SubType"]])
#DE_table_tip_tuft = FindMarkers(object = seuratObj_receiver, ident.1 = "Tuft", ident.2 = "Tip ECs",min.pct = -Inf, logfc.threshold = -Inf,only.pos = TRUE) %>% rownames_to_column("gene")

#seuratObj_receiver= subset(seuratObj_receiver, idents = c("Tip ECs"))

condition_oi_receive = "OIR_P14_WT"
condition_reference_receive = "OIR_P14_Sirt3KO" 

seuratObj_receiver = SetIdent(seuratObj_receiver, value = seuratObj_receiver[["Condition"]])

DE_table_receiver = FindMarkers(object = seuratObj_receiver, ident.1 = condition_oi_receive, ident.2 = condition_reference_receive,min.pct = 0.1, logfc.threshold = Log2FC_var,only.pos = TRUE, assay_oi="RNA") %>% rownames_to_column("gene")
#3000
geneset_oi=DE_table_receiver$gene


###now setting very specific receiver/receptors and background
background_expressed_genes = setdiff(expressed_genes_receiver,geneset_oi)
#1109
expressed_receptors=intersect(receptors,geneset_oi)
list(expressed_receptors)

# [1] "CXCR4"    "LIFR"     "FLT1"     "KDR"      "IFNGR2"   "IL10RB"  
#  [7] "CD40"     "ICAM1"    "VCAM1"    "CDH2"     "OCLN"     "F11R"    
# [13] "GLG1"     "PTPRM"    "ITGA1"    "EDNRB"    "SCARF1"   "EPHA4"   
# [19] "AMFR"     "BAMBI"    "NRP1"     "NCSTN"    "NRP2"     "TNFRSF21"
# [25] "FZD8"     "RYK"      "ICAM2"    "PTK7"     "ITPR3"    "ITPR1"   
# [31] "ENG"      "GPR4"     "RARG"     "PTPRK"    "RRBP1"    "PRNP"    
# [37] "TRIM27"   "PCDH7"    "PKD1"     "XPR1"     "CNGB1"   





------

#sender
#find de genes of senders now

seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[["Cell_Type"]])

#setting parameters
condition_oi = "OIR_WT"
condition_reference = "OIR_Sirt3KO" 

# seurat_obj_ec= subset(seuratObj, idents = "ECs")
# seurat_obj_ec = SetIdent(seurat_obj_ec, value = seurat_obj_ec[["Condition"]])
# DE_table_ec = FindMarkers(object = seurat_obj_ec, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = -Inf, logfc.threshold = Log2FC_var,only.pos = TRUE) %>% rownames_to_column("gene")

seurat_obj_mural= subset(seuratObj_sender, idents = "Mural cells")
seurat_obj_mural = SetIdent(seurat_obj_mural, value = seurat_obj_mural[["Dataset"]])
DE_table_mural = FindMarkers(object = seurat_obj_mural, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = pct_var, logfc.threshold = Log2FC_var,only.pos = TRUE, assay_oi="RNA") %>% rownames_to_column("gene")

# seurat_obj_immune= subset(seuratObj, idents = "Immune cells")
# seurat_obj_immune = SetIdent(seurat_obj_immune, value = seurat_obj_immune[["Condition"]])
# DE_table_immune = FindMarkers(object = seurat_obj_immune, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = pct_var, logfc.threshold = Log2FC_var,only.pos = TRUE, assay_oi="RNA") %>% rownames_to_column("gene")

seurat_obj_astro= subset(seuratObj_sender, idents = "Astrocytes")
seurat_obj_astro = SetIdent(seurat_obj_astro, value = seurat_obj_astro[["Dataset"]])
DE_table_astro = FindMarkers(object = seurat_obj_astro, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = pct_var, logfc.threshold = Log2FC_var,only.pos = TRUE, assay_oi="RNA") %>% rownames_to_column("gene")

# seurat_obj_lens= subset(seuratObj, idents = "Lens cells")
# seurat_obj_lens = SetIdent(seurat_obj_lens, value = seurat_obj_lens[["Condition"]])
# DE_table_lense = FindMarkers(object = seurat_obj_lens, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = -Inf, logfc.threshold = Log2FC_var,only.pos = TRUE) %>% rownames_to_column("gene")

seurat_obj_muller= subset(seuratObj_sender, idents = "Muller Glial cells")
seurat_obj_muller = SetIdent(seurat_obj_muller, value = seurat_obj_muller[["Dataset"]])
DE_table_muller = FindMarkers(object = seurat_obj_muller, ident.1 = condition_oi, ident.2 = condition_reference,
  min.pct = pct_var, logfc.threshold = Log2FC_var ,only.pos = TRUE, assay_oi="RNA") %>% rownames_to_column("gene")



# seurat_obj_rbc= subset(seuratObj, idents = "RBCs")
# seurat_obj_rbc = SetIdent(seurat_obj_rbc, value = seurat_obj_rbc[["Condition"]])
# DE_table_rbc = FindMarkers(object = seurat_obj_rbc, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = pct_var, logfc.threshold = Log2FC_var,only.pos = TRUE, assay_oi="RNA") %>% rownames_to_column("gene")

# seurat_obj_cones= subset(seuratObj, idents = "Cones")
# seurat_obj_cones = SetIdent(seurat_obj_cones, value = seurat_obj_cones[["Condition"]])
# DE_table_cones = FindMarkers(object = seurat_obj_cones, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = -Inf, logfc.threshold = Log2FC_var,only.pos = TRUE) %>% rownames_to_column("gene")
##cell group 1 has fewer than 3 cells

# seurat_obj_bi= subset(seuratObj, idents = "Bipolar cells")
# seurat_obj_bi = SetIdent(seurat_obj_bi, value = seurat_obj_bi[["Condition"]])
# DE_table_bi = FindMarkers(object = seurat_obj_bi, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = pct_var, logfc.threshold = Log2FC_var,only.pos = TRUE, assay_oi="RNA") %>% rownames_to_column("gene")

# seurat_obj_am= subset(seuratObj, idents = "Amacrine cells")
# seurat_obj_am = SetIdent(seurat_obj_am, value = seurat_obj_am[["Condition"]])
# DE_table_am = FindMarkers(object = seurat_obj_am, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = -Inf, logfc.threshold = Log2FC_var,only.pos = TRUE) %>% rownames_to_column("gene")
##cell group 1 has fewer than 3 cells

# seurat_obj_rods= subset(seuratObj, idents = "Rods")
# seurat_obj_rods = SetIdent(seurat_obj_rods, value = seurat_obj_rods[["Condition"]])
# DE_table_rods = FindMarkers(object = seurat_obj_rods, ident.1 = condition_oi, ident.2 = condition_reference,min.pct = pct_var, logfc.threshold = Log2FC_var,only.pos = TRUE, assay_oi="RNA") %>% rownames_to_column("gene")

DE_table_astro[grep("VEGFA", DE_table_astro$gene),]
#       gene     p_val avg_log2FC pct.1 pct.2 p_val_adj
#3 VEGFA 0.003608842   2.044399   0.6   0.1         1

DE_table_muller[grep("VEGFA", DE_table_muller$gene),]
#      gene p_val avg_log2FC pct.1 pct.2 p_val_adj
# 1742 CALR     1  0.9934302 0.429 0.444         1

DE_table_mural[grep("VEGFA", DE_table_mural$gene),]
#     gene      p_val avg_log2FC pct.1 pct.2 p_val_adj
# 101 APOE 0.03541451  0.9670743 0.762 0.806         1

#expressed_ligands_ec= intersect(ligands,DE_table_ec$gene)
expressed_ligands_mural= intersect(ligands,DE_table_mural$gene)
#expressed_ligands_immune= intersect(ligands,DE_table_immune$gene)
#expressed_ligands_lens= intersect(ligands,DE_table_lense$gene)
expressed_ligands_muller= intersect(ligands,DE_table_muller$gene)
#expressed_ligands_rbc= intersect(ligands,DE_table_rbc$gene)
#expressed_ligands_bi= intersect(ligands,DE_table_bi$gene)
#expressed_ligands_rods= intersect(ligands,DE_table_rods$gene)
expressed_ligands_astro= intersect(ligands,DE_table_astro$gene)


expressed_ligands = unique(c(#expressed_ligands_ec,
expressed_ligands_mural,
#expressed_ligands_immune,
#expressed_ligands_lens,
expressed_ligands_muller,
#expressed_ligands_rbc,
#expressed_ligands_bi,
expressed_ligands_astro
#expressed_ligands_rods
))
expressed_ligands[grep("VEGFA", expressed_ligands)]
expressed_ligands[grep("IGF", expressed_ligands)]


-------


###ligand focus

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
# # A tibble: 6 × 4
#   from  to     source           database      
#   <chr> <chr>  <chr>            <chr>         
# 1 CDH2  CDH2   kegg_cams        kegg          
# 2 CALR  SCARF1 ramilowski_known ramilowski    
# 3 APP   IGF1R  ppi_lr           ppi_prediction
# 4 APP   LPAR4  ppi_lr           ppi_prediction
# 5 VEGFA FLT4   ppi_lr           ppi_prediction
# 6 CALR  CD40   ppi_lr           ppi_prediction

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

ligand_activities_oir = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities_oir = ligand_activities_oir %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
best_upstream_ligands_oir = ligand_activities_oir%>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

best_upstream_ligands_oir
#[1] "APP"   "CALR"  "VEGFA" "CDH2"  "ITGB1"


pdf ("DotPlot_best_upstream_ligands_RdYlGn_2.pdf",width=10, height=5)
DotPlot(seuratObj_sender, features = best_upstream_ligands_oir %>% rev(), cols = "RdYlGn", split.by="Dataset") + RotatedAxis()
dev.off()



----------
###ligand gene target

active_ligand_target_links_df_oir = best_upstream_ligands_oir %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links_oir = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df_oir, ligand_target_matrix = ligand_target_matrix, cutoff = 0.9)

order_ligands_oir = intersect(best_upstream_ligands_oir, colnames(active_ligand_target_links_oir)) %>% rev() %>% make.names()
order_targets_oir = active_ligand_target_links_df_oir$target %>% unique() %>% intersect(rownames(active_ligand_target_links_oir)) %>% make.names()
rownames(active_ligand_target_links_oir) = rownames(active_ligand_target_links_oir) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links_oir) = colnames(active_ligand_target_links_oir) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target_oir = active_ligand_target_links_oir[order_targets_oir,order_ligands_oir] %>% t()
p_ligand_target_network_oir = vis_ligand_target_oir 

dim(p_ligand_target_network_oir)
#[1]  3 16

#top 20 ranked ligands and their predicted target genes downstream

pdf ("heatmap_ligands_predictedGeneTarget_downstream_dif_0.9.pdf",width=7, height=3)
make_heatmap_ggplot(p_ligand_target_network_oir,"Prioritized ligands","Predicted target genes", color = "tomato",legend_position = "top", 
  x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
dev.off()


---------
###receptor focus


lr_network_top_oir = lr_network %>% filter(from %in% best_upstream_ligands_oir & to %in% expressed_receptors) %>% distinct(from,to)

best_upstream_receptors_oir = lr_network_top_oir %>% pull(to) %>% unique()

#weighted consideration based on lit
lr_network_top_df_large_oir = weighted_networks_lr %>% filter(from %in% best_upstream_ligands_oir & to %in% best_upstream_receptors_oir)

lr_network_top_df_oir = lr_network_top_df_large_oir %>% spread("from","weight",fill = 0)
lr_network_top_matrix_oir = lr_network_top_df_oir %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_oir$to)

##receptors
dist_receptors_oir = dist(lr_network_top_matrix_oir, method = "binary")
hclust_receptors_oir = hclust(dist_receptors_oir, method = "ward.D2")
order_receptors_oir= hclust_receptors_oir$labels[hclust_receptors_oir$order]


##ligands and ligand receptor
dist_ligands_oir = dist(lr_network_top_matrix_oir %>% t(), method = "binary")
hclust_ligands_oir = hclust(dist_ligands_oir, method = "ward.D2")
order_ligands_receptor_oir = hclust_ligands_oir$labels[hclust_ligands_oir$order]


##find genes of ligands and receptor
order_receptors_oir = order_receptors_oir %>% intersect(rownames(lr_network_top_matrix_oir))
order_ligands_receptor_oir = order_ligands_receptor_oir %>% intersect(colnames(lr_network_top_matrix_oir))


vis_ligand_receptor_network_oir = lr_network_top_matrix_oir[order_receptors_oir, order_ligands_receptor_oir]
rownames(vis_ligand_receptor_network_oir) = order_receptors_oir %>% make.names()
colnames(vis_ligand_receptor_network_oir) = order_ligands_receptor_oir %>% make.names()


p_ligand_receptor_network_oir = vis_ligand_receptor_network_oir %>% t()
dim(p_ligand_receptor_network_oir)
#[1]  6 11

pdf("heatmap_ligand_receptors.pdf", height=4, width=5)
make_heatmap_ggplot(p_ligand_receptor_network_oir,"Ligands","Receptors", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Prior interaction potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple")
dev.off()


---------
#overall ligand activity heatmap

ligand_pearson_matrix_oir = ligand_activities_oir %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_oir$test_ligand)
rownames(ligand_pearson_matrix_oir) = rownames(ligand_pearson_matrix_oir) %>% make.names()
colnames(ligand_pearson_matrix_oir) = colnames(ligand_pearson_matrix_oir) %>% make.names()


vis_ligand_pearson_oir = ligand_pearson_matrix_oir[order_ligands_oir, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson_oir = vis_ligand_pearson_oir 


pdf ("ligand_activity_heatmap.pdf")
make_heatmap_ggplot(p_ligand_pearson_oir,"Prioritized ligands","Ligand activity", color = "purple",legend_position = "bottom", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
dev.off()


###dot plot of ligands from cell type
pdf ("ligand_cell_type.pdf",width=15,height=5)
DotPlot(seuratObj_sender, features = best_upstream_ligands_oir %>% rev(), cols = "RdYlBu") + RotatedAxis()
dev.off()

-----


########################
###circos
########################

---------

#best_upstream_ligands_ec = best_upstream_ligands_oir %>% intersect(expressed_ligands_ec) 
#best_upstream_ligands_mural = best_upstream_ligands_oir %>% intersect(expressed_ligands_mural) 
#best_upstream_ligands_immune = best_upstream_ligands_oir %>% intersect(expressed_ligands_immune) 
#best_upstream_ligands_lens = best_upstream_ligands_oir %>% intersect(expressed_ligands_lens) 
#best_upstream_ligands_muller = best_upstream_ligands_oir %>% intersect(expressed_ligands_muller) 
#best_upstream_ligands_rbc = best_upstream_ligands_oir %>% intersect(expressed_ligands_rbc) 
#best_upstream_ligands_bi = best_upstream_ligands_oir %>% intersect(expressed_ligands_bi) 
#best_upstream_ligands_rods = best_upstream_ligands_oir %>% intersect(expressed_ligands_rods) 


expression = t(as.matrix(GetAssayData(seuratObj_sender, assay="RNA", slot="data")))

seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[["Cell_Type"]])

#id_ecs= subset(seuratObj, idents = "ECs")%>% colnames()
id_mural= subset(seuratObj_sender, idents = "Mural cells")%>% colnames()
#id_immune= subset(seuratObj, idents = "Immune cells")%>% colnames()
id_astro= subset(seuratObj_sender, idents = "Astrocytes")%>% colnames()
#id_lens= subset(seuratObj, idents = "Lens cells")%>% colnames()
id_muller= subset(seuratObj_sender, idents = "Muller Glial cells")%>% colnames()
#id_rbc= subset(seuratObj, idents = "RBCs")%>% colnames()
#id_bi= subset(seuratObj, idents = "Bipolar cells")%>% colnames()
#id_rods= subset(seuratObj, idents = "Rods")%>% colnames()

ligand_expression_tbl = tibble(
  ligand = best_upstream_ligands_oir, 
#  ECs = expression[id_ecs,best_upstream_ligands_oir] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
  Mural_cells = expression[id_mural,best_upstream_ligands_oir] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
#  Immune_cells = expression[id_immune,best_upstream_ligands_oir] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
  Astrocytes = expression[id_astro,best_upstream_ligands_oir] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
#  Lens_cells = expression[id_lens,best_upstream_ligands_oir] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
  Muller_Glial_cells = expression[id_muller,best_upstream_ligands_oir] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)})
#  RBCs = expression[id_rbc,best_upstream_ligands_oir] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
#  Bipolar_cells = expression[id_bi,best_upstream_ligands_oir] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)})
#  Rods = expression[id_rods,best_upstream_ligands_oir] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)})
 )

###Ligand more expressed compared to other cell types

#Seurat_Markers_Mural <- FindAllMarkers(seuratObj_sender, only.pos=TRUE,logfc.threshold = Log2FC_var, min.pct = pct_var)
##only take wt

#EC_specific_ligands <- intersect(ligand_expression_tbl$ligand,subset(Seurat_Markers, subset = cluster =="ECs")$gene)



mural_specific_ligands <- intersect(ligand_expression_tbl$ligand,DE_table_mural$gene)
#immune_specific_ligands <- intersect(ligand_expression_tbl$ligand,subset(Seurat_Markers, subset = cluster =="Immune cells")$gene)
astro_specific_ligands <- intersect(ligand_expression_tbl$ligand,DE_table_astro$gene)
#lens_specific_ligands <- intersect(ligand_expression_tbl$ligand,subset(Seurat_Markers, subset = cluster =="Lens cells")$gene)
muller_specific_ligands <- intersect(ligand_expression_tbl$ligand,DE_table_muller$gene)
#rbc_specific_ligands <- intersect(ligand_expression_tbl$ligand,subset(Seurat_Markers, subset = cluster =="RBCs")$gene)
#bi_specific_ligands <- intersect(ligand_expression_tbl$ligand,subset(Seurat_Markers, subset = cluster =="Bipolar cells")$gene)
#rod_specific_ligands <- intersect(ligand_expression_tbl$ligand,subset(Seurat_Markers, subset = cluster =="Rods")$gene)


general_ligands = setdiff(best_upstream_ligands_oir,
    c(#EC_specific_ligands,
    mural_specific_ligands, 
   # immune_specific_ligands,
    astro_specific_ligands,
    #lens_specific_ligands,
    muller_specific_ligands
  #rbc_specific_ligands,
#  bi_specific_ligands
  #rod_specific_ligands
  ))



ligand_type_indication_df = tibble(
  ligand_type = c(#rep("ECs", times = EC_specific_ligands %>% length()),
                  rep("Mural cells", times = mural_specific_ligands %>% length()),
#                  rep("Immune cells", times = immune_specific_ligands %>% length()),
                rep("Astrocytes", times = astro_specific_ligands %>% length()),
#                rep("Lens cells", times = lens_specific_ligands %>% length()),
                rep("Muller Glial cells", times = muller_specific_ligands %>% length()),
#                rep("RBCs", times = rbc_specific_ligands %>% length()),
#                rep("Bipolar cells", times = bi_specific_ligands %>% length()),
#                rep("Rods", times = rod_specific_ligands %>% length()),
                rep("General", times = general_ligands %>% length())

                ),
  ligand = c(#EC_specific_ligands,
    mural_specific_ligands, 
#    immune_specific_ligands,
    astro_specific_ligands,
#    lens_specific_ligands,
    muller_specific_ligands,
#    rbc_specific_ligands,
#    bi_specific_ligands,
#    rod_specific_ligands,
    general_ligands))

#37 rows

#targets genes (downstream effects) of ligands
active_ligand_target_links_df = best_upstream_ligands_oir %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "Tip ECs") %>% inner_join(ligand_type_indication_df) 
active_ligand_target_links_df$ligand <- paste(active_ligand_target_links_df$ligand_type,active_ligand_target_links_df$ligand,sep="_")


----
library(RColorBrewer)
grid_col_ligand =c("ECs" = "indianred", 
                    "Mural cells"="blue4", 
                    "Immune cells"="blue",
                    "Astrocytes"="mediumseagreen",
                    "Lens cells"="darkseagreen3",
                    "Muller Glial cells"="lightsalmon",
                    "RBCs"="darkgreen",
                    "Bipolar cells"="darkorchid1",
                    "Rods"="brown4",
                    "General"="red")
#set your targets = p-emt gene pathway in malignant receiver cells
grid_col_target =c(
            "Tip ECs" = "indianred")
grid_col_receptor =c(
            "Tip ECs" = "indianred")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)


---------
#########################

##ligand target

#########################

--------

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0,na.rm = TRUE)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)
#78
write.csv(circos_links,file="p14_tip_LT_nothres.csv")



------
#join according to colour
circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)

#makes the links/lines between the ligand and target 
links_circle = circos_links %>% select(ligand,target, weight)

##aesthetic stuff of the link line
ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)


# give the option that links in the circos plot will be transparant ~ ligand-target potential score
                    ####check with gael

transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

target_order = unique(circos_links$target)
ligand_order <- sort(unique(circos_links$ligand))
order = c(ligand_order,target_order)



pdf("circos_p14_tip_Lt_nothres_grad.pdf",width =50,height=55)
circos.par(track.height = 0.1)
chordDiagram(links_circle, directional = 1,order=order,
    link.sort = TRUE, 
    link.decreasing = FALSE, 
    grid.col = grid_col,
    diffHeight = 0.005, 
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow", 
    link.visible = TRUE,
    annotationTrack = "grid", 
    transparency = transparency,
    scale=FALSE,
    preAllocateTracks = list(track.height = 0.075))

# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 3)
}, bg.border = NA) 


#draw(lgd_cells, x = unit(20, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))
dev.off()

pdf("circos_p14_tip_Lt_nothres.pdf",width =50,height=55)
circos.par(track.height = 0.1)
chordDiagram(links_circle, directional = 1,order=order,
    link.sort = TRUE, 
    link.decreasing = FALSE, 
    grid.col = grid_col,
    diffHeight = 0.005, 
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow", 
    link.visible = TRUE,
    annotationTrack = "grid", 
    transparency = 0,
    scale=FALSE,
    preAllocateTracks = list(track.height = 0.075))

# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 3)
}, bg.border = NA) 


#draw(lgd_cells, x = unit(20, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))
dev.off()




---------
#########################

##ligand receptor

#########################

--------
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands_oir & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df = weighted_networks_lr%>% filter(from %in% best_upstream_ligands_oir & to %in% best_upstream_receptors) %>% rename(ligand = from, receptor = to)

lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "Tip ECs") %>% inner_join(ligand_type_indication_df)

circos_links = lr_network_top_df
circos_links$ligand <- paste(circos_links$ligand_type,circos_links$ligand,sep="_")


cutoff_include_all_ligands = circos_links$weight %>% quantile(0,na.rm = TRUE)

circos_links = circos_links %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(circos_links$ligand %>% unique(), circos_links$ligand %>% unique())
receptors_to_remove = setdiff(circos_links$receptor %>% unique(), circos_links$receptor %>% unique())

circos_links = circos_links %>% filter(!receptor %in% receptors_to_remove &!ligand %in% ligands_to_remove)
#11

write.csv(circos_links,file="p14_tip_LR_nothres.csv")


circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)

#makes the links/lines between the ligand and target 
links_circle = circos_links %>% select(ligand,receptor, weight)


ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

receptor_order = circos_links$receptor %>% unique()
ligand_order <- sort(unique(circos_links$ligand))
order = c(ligand_order,receptor_order)


pdf("circos_p14_tip_LR_nothres_grad.pdf",width =50,height=60)
circos.par(track.height = 0.1)
chordDiagram(links_circle, directional = 1,order=order,
    link.sort = TRUE, 
    link.decreasing = FALSE, 
    grid.col = grid_col,
    diffHeight = 0.005, 
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow", 
    link.visible = TRUE,
    annotationTrack = "grid", 
    transparency = transparency,
    scale=FALSE,
    preAllocateTracks = list(track.height = 0.075))

# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 3)
}, bg.border = NA) 


#draw(lgd_cells, x = unit(20, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))
dev.off()



pdf("circos_p14_tip_LR_nothres.pdf",width =50,height=60)
circos.par(track.height = 0.1)
chordDiagram(links_circle, directional = 1,order=order,
    link.sort = TRUE, 
    link.decreasing = FALSE, 
    grid.col = grid_col,
    diffHeight = 0.005, 
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow", 
    link.visible = TRUE,
    annotationTrack = "grid", 
    transparency = 0,
    scale=FALSE,
    preAllocateTracks = list(track.height = 0.075))

# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 3)
}, bg.border = NA) 


#draw(lgd_cells, x = unit(20, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))
dev.off()


