#!/usr/bin/env Rscript

#################################################################
# Description: Clustering and evaluation of subclusters using CDI
#################################################################

# Following script contains the full clusteringpipeline for the Ubels26_HairCycle
# publication. Renv/Conda environments are dynamically set. Pytorch compatability
# has to be set by user and will not be supported.
#################################################################
# LIBRARY LOADING
#################################################################

library(Seurat)
library(SeuratDisk)
library(Nebulosa)
library(reticulate)
library(XML)
library(GenomicFeatures)
library(ggplot2)
library(scCustomize)

options(future.globals.maxSize = 5 * 1024^3)
set.seed(43)

# There may be an incompatability caused by XML due to libxml2
# In this case simply force install of libxml2 to the latest version
# and run renv::install("XML", type = "source", rebuild = TRUE) to
# force reinstall into the renv.

#################################################################
# SETUP PROJECT PARAMETERS
#################################################################
project <- "scCD83_Mallia"
main_folder <- "./"

source("./scripts/helper_functions.R")
source("./scripts/gene_lists.R")
source("./scripts/setup_py_environment.R")

#################################################################
# SETTING UP DATA
#################################################################

object_list <- readRDS(paste0(main_folder, "obj_doublet_filtered.rds"))
object_list <- filter_doublets(object_list)
unique(object_list[[1]]$scDblFinder.class)

#################################################################
# HARMONY INTEGRATION
#################################################################

obj <- merge_and_preprocess(object_list)
obj <- run_harmony_integration(obj)
obj$treatment <- ifelse(grepl("^PBS", obj$orig.ident), "PBS", "sCD83")

plot_integration_comparison(obj, output_path = paste0(main_folder, "/integration_plots/harmony_comparison.png"))
saveRDS(obj, paste0(main_folder, "obj_harmony_integrated.rds"))

rm(object_list)
gc()

#obj <- readRDS("obj_harmony_integrated.rds")

vis_obj <- readRDS(paste0("../ubels26_hair_cycle/annotated_vis_obj.rds"))
ImageFeaturePlot(vis_obj, features = c("NR4A2"))

#################################################################
# RUNNING BROAD MARKER GENES FOR INITIAL CLUSTERIZATION
#################################################################

# Easy visualization through Nebulosa to get a better overview of gene expression particularly
# for when low cell count has high gene expression in a particular cluster
obj <- cluster_subcluster(obj, output_dir = "./marker_genes/")

Idents(obj) <- "SCT_snn_res.0.8"
broad_markers <- FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.5)

plot_marker_genes(obj = obj, 
                  genes = broad_gene_list, 
                  cluster_col = "SCT_snn_res.0.8",
                  reduction = "umap", 
                  output_dir = "./marker_genes/broad_markers", 
                  pt_size = 1,
                  show_labels = TRUE
                  )

################################################################
# SHOW AND REMOVE CLUSTER SPECIFIC TO PBS_72H 
################################################################

sample_colors <- c(
  "PBS_48h"   = "#4E79A7",
  "PBS_72h"   = "#A0CBE8",
  "sCD83_48h" = "#E15759",
  "sCD83_72h" = "#F1A2A0"
)

df <- data.frame(
  cluster = obj$SCT_snn_res.0.8,
  sample = obj$orig.ident
)

ggplot(df, aes(x = cluster, fill = sample)) +
  geom_bar(position = "stack") +
  theme_minimal() +
  scale_fill_manual(values = sample_colors) +
  labs(x = "Cluster", y = "Cell count", fill = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(df, aes(x = cluster, fill = sample)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  scale_fill_manual(values = sample_colors) +
  labs(x = "Cluster", y = "Proportion", fill = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

obj <- subset(obj, idents = "0", invert = TRUE)

#################################################################
# RECLUSTER WITHOUT ABBARENT CELLTYPE
#################################################################

# Easy visualization through Nebulosa to get a better overview of gene expression particularly
# for when low cell count has high gene expression in a particular cluster
obj <- cluster_subcluster(obj, output_dir = "./marker_genes/")

Idents(obj) <- "SCT_snn_res.0.8"
broad_markers <- FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.5)

marker_genes <- broad_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 200) -> broad_markers

write.csv(marker_genes, file = "marker_genes.csv")

plot_marker_genes(obj = obj, 
                  genes = broad_gene_list, 
                  cluster_col = "SCT_snn_res.0.8",
                  reduction = "umap", 
                  output_dir = "./marker_genes/broad_markers", 
                  pt_size = 1,
                  show_labels = TRUE
)



#################################################################
# ASSIGN BROAD MARKER IDENTIFICATION TO CLUSTERS
#################################################################

broad_cluster_identification <- list(
  `0` = "Fibroblasts.1",          #PDGFRA/FAP/BICC1             
  `1` = "Keratinocytes.1",        #KRT5/KRT14/KRT17             
  `2` = "Endothelial.1",          #PECAM1/VWF/KDR/PLVAP/CLDN5
  `3` = "Fibroblasts.2",          #COL1A1/COL3A1/DCN/LUM/POSTN/TNC 
  `4` = "Immune.2",               #PTPRC/CD63/CD96
  `5` = "Fibroblasts.3",          #RGS5/KCNJ8/ABCC9/MYOCD/ANGPT1
  `6` = "Fibroblasts.4",          #MMP1/MMP3/MMP10/CXCL5/CXCL6/PDGFRA
  `7` = "Fibroblasts.5",          #VCAM1/ACHE/NPR3/HES4/LRRC32/PER1   
  `8` = "Immune.3",               #PTPRC/CD63/CD96
  `9` = "Remaining.1",            #MLANA/SOX5/FMN1        
  `10` = "Keratinocytes.2",       #KRT6A/KRT16/KRT17        
  `11` = "Fibroblasts.6",         #DDX21/WDR43/NOP16/SNHG15/POSTN  
  `12` = "Endothelial.2",         #PECAM1/ERG/NOTCH4/CDH5
  `13` = "Keratinocytes.3",       #KRT85/KRT35/KRT32
  `14` = "Fibroblasts.7",         #CXCL1/2/3/8/IL6/FGF7/WNT2   
  `15` = "Immune.4",              #MRC1/HLA-DRA/CD14
  `16` = "Keratinocytes.4",       #EGFR/COL17A1/CDH13                  
  `17` = "Fibroblasts.8",         #CEMIP/TWIST2/RSPO3/FGF2/SPON1
  `18` = "Immune.5",              #PTPRC/CD96/TOX
  `19` = "Immune.6",              #TPSAB1/TPSB2/CPA3/CMA1/HDC
  `20` = "Keratinocytes.5",       #KRT5/KRT14/S100A9        
  `21` = "Keratinocytes.6",       #KRT19/SAA1/MUCL1           
  `22` = "Fibroblasts.9",         #PDGFRA/FAP/BICC1
  `23` = "Fibroblasts.10",        #PDGFRA/FAP/BICC1
  `24` = "Immune.7",              #HLA-DRA/CD74/CD86/LAMP3
  `25` = "Remaining.2",           #CDH19/SOX10/
  `26` = "Immune.8"               #IGKC/IGHG1/MZB1
)

obj$broad_cluster <- unname(unlist(broad_cluster_identification[as.character(obj$SCT_snn_res.0.8)]))
Idents(obj) <- "broad_cluster"

p <- DimPlot(obj, label = TRUE, repel = TRUE, label.size = 2)
ggsave(filename = "./marker_genes/broad_cluster_annotation.png", p, width = 15, height = 10)

####################################################################################
# PLOT BROAD MARKER GENES DOTPLOT
####################################################################################

group_def <- split(
  unique(obj$broad_cluster),
  sub("\\.\\d+$", "", unique(obj$broad_cluster))
)

# Natural sort within each group
group_def <- lapply(group_def, function(x) {
  nums <- as.numeric(sub(".*\\.", "", x))
  x[order(nums)]
})

group_def <- group_def[c("Keratinocytes","Fibroblasts","Endothelial","Immune","Remaining")]
names(group_def) <- c("Keratinocytes","Fibroblasts","Endo","Immune","Other")

p <- grouped_dotplot(
  obj,
  group_def = group_def,
  gene_groups = list(
    Keratinocytes = c("KRT5","KRT14","KRT17","KRT6A","S100A2","S100A9",
                      "COL17A1","TP63","KRT85","KRT19"),
    Fibroblasts   = c("PDGFRA","COL1A1","COL3A1","DCN","LUM","FAP","MME",
                      "RGS5","PRRX1","MMP1","MMP3","VCAM1"),
    Endo          = c("PECAM1","VWF","CDH5","CLDN5","ERG","DLL4"),
    Immune        = c("CD3D","CD3E","PTPRC","CD14","C1QA",
                      "TPSAB1","CPA3","LAMP3","CCR7","IDO1","IGKC","MZB1"),
    Other         = c("MLANA","DCT","TYRP1","SOX10","CDH19","PMP2")
  ),
  ident_col = "broad_cluster",
  bar_height = 1.2
)

p

ggsave("./marker_genes/grouped_dotplot.pdf", p, width = 14, height = 10)

saveRDS(obj, file = "./broad_annotated_obj.rds")

####################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO KERATINOCYTE CLUSTERS
####################################################################

obj$fine_clust <- as.character(obj$broad_cluster)

krt_obj <- subset(obj, idents = grep("^Keratinocytes", levels(Idents(obj)), value = TRUE))
krt_obj <- cluster_subcluster(krt_obj, output_dir = "./marker_genes/keratinocytes/", dims = 1:30)
Idents(krt_obj) <- "SCT_snn_res.0.4"

krt_markers <- FindAllMarkers(krt_obj, min.pct = 0.1, logfc.threshold = 0.5)

krt_cluster_identification <- list(
  `0` = "ORS.1",              #CXCL8/MMP3/CXCL1                   
  `1` = "ORS.2",              #CXCL14/S100A2   
  `2` = "ORS.Suprabasal",     #KRT6A/KRT16/KRT17       
  `3` = "ORS.3",              #KRT5/KRT14/GAS5
  `4` = "Bulge",              #EGFR/CDH13/RUNX1
  `5` = "FBs.Reticular",      #COL1A2/IGFBP5/VIM
  `6` = "Eccrine.cells",      #KRT19/AQP5/MUCL1  
  `7` = "Ishtmus",            #KLK10/MXD1/FOXO3/PRDM1
  `8` = "Matrix"              #KRT35/KRT85
)

new_ids <- unlist(krt_cluster_identification)

krt_obj$fine_clust <- plyr::mapvalues(
  x = as.character(Idents(krt_obj)),  
  from = names(new_ids),
  to = new_ids
)

obj$fine_clust[colnames(krt_obj)] <- krt_obj$fine_clust

####################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO FIBROBLAST CLUSTERS
####################################################################

Idents(obj) <- "broad_cluster"
fb_obj <- subset(obj, idents = grep("^Fibroblasts", levels(Idents(obj)), value = TRUE))
fb_all_barcodes <- colnames(fb_obj)

fb_obj <- cluster_subcluster(fb_obj, output_dir = "./marker_genes/fibroblasts/", n_genes = 2000, features = 2000, dims = 1:20, conserve.memory = FALSE)
Idents(fb_obj) <- "SCT_snn_res.0.5"
fb_markers <- FindAllMarkers(fb_obj, min.pct = 0.2, logfc.threshold = 1, only.pos = TRUE)

g1 <- DimPlot(fb_obj, group.by = "SCT_snn_res.0.5", label = TRUE)
g2 <- DimPlot(fb_obj, group.by = "Phase")
(g1 + g2)

fb_genes <- fb_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 300) -> fb_markers

write.csv(fb_markers, "fb_marker_genes.csv")

fb_obj <- subset(fb_obj, idents = "8", invert = TRUE)

fb_cluster_identification <- list(
  `0` = "FBs.Activated",           #CEMIP/MT1X/IGFBP2
  `1` = "FBs.Interstitial",        #COL1A1/COL6A3/GAS1
  `2` = "FBs.Infl.Myofibroblast",  #MMP3/MMP1/CXCL5 — maps to Steele F6
  `3` = "FBs.Pericytes",           #PI15/ESAM/ABCC9/RGS5
  `4` = "FBs.DP-like",             #VCAM1/NGFR/OXTR
  `5` = "FBs.Perifolicular",       #NTN1/PLXDC2/LAMA2
  `6` = "FBs.Cycling",             #DDX21/WDR43/ODC1
  `7` = "FBs.NF-kB",               #CXCL2/CXCL3/RCSD1/NFKB
  `8` = "FBs.HF.Progenitor",       #HSPA6/HES5/CRYAB
  `9` = "FBs.Dermal.Sheath"        #COL11A1/TAGLN/PPP1R14A
)

new_ids <- unlist(fb_cluster_identification)
fb_obj$fine_clust <- plyr::mapvalues(
  x    = as.character(Idents(fb_obj)),
  from = names(new_ids),
  to   = new_ids
)

obj$fine_clust[colnames(fb_obj)] <- fb_obj$fine_clust

###############################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TOIMMUNE CLUSTERS
###############################################################################

Idents(obj) <- "broad_cluster"
im_obj <- subset(obj, idents = grep("^Immune", levels(Idents(obj)), value = TRUE))
im_all_barcodes <- colnames(im_obj)

im_obj <- cluster_subcluster(im_obj, output_dir = "./marker_genes/immune/", n_genes = 2000, features = 2000, dims = 1:20, conserve.memory = FALSE)
Idents(im_obj) <- "SCT_snn_res.0.7"
im_markers <- FindAllMarkers(im_obj, min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)

g1 <- DimPlot(im_obj, group.by = "SCT_snn_res.0.7", label = TRUE)
g2 <- FeaturePlot(im_obj, features = "CCR7")
(g1 + g2)

sub_markers <- FindMarkers(im_obj, ident.1 = c("0", "1", "2", "4", "10"))

im_genes <- im_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 300) -> im_markers

write.csv(im_markers, "im_marker_genes.csv")

im_cluster_identification <- list(
  `0` = "T.Naive-ANK3",            #PTPRC/CD96/IL7R/ANK3
  `1` = "T.Naive",                 #PTPRC/CD96/IL7R
  `2` = "T.Naive",                 #PTPRC/CD96/IL7R
  `3` = "M2.Macrophages",          #MRC1/CCL3/CD83
  `4` = "T.Naive-ANK3",            #PTPRC/CD96/ANK3
  `5` = "MAIT.cells",              #KLRB1/BLK
  `6` = "Mast.cells",              #CPE3/KIT 
  `7` = "T.reg",                   #FOXP3/IKZF2
  `8` = "T.cyto",                  #CD8A/NKG7
  `9` = "Dendritic.cells",         #HLA-DRA/LAMP3
  `10` = "T.Naive",                #PTPRC/CD96/IL7R
  `11` = "FBs.Dermal.Sheath",      #TAGLN/COL1A1/CALD
  `12` = "Plasma.cells",           #IGKC/JCHAIN/MZB1
  `13` = "B.cells"                 #CD19/MS4A1/PAX5
)

new_ids <- unlist(im_cluster_identification)
im_obj$fine_clust <- plyr::mapvalues(
  x    = as.character(Idents(im_obj)),
  from = names(new_ids),
  to   = new_ids
)

obj$fine_clust[colnames(im_obj)] <- im_obj$fine_clust

DimPlot(obj, group.by = "fine_clust", label = TRUE)

###############################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION TO ENDOTHELIAL CLUSTERS
###############################################################################

Idents(obj) <- "broad_cluster"
en_obj <- subset(obj, idents = grep("^Endothelial", levels(Idents(obj)), value = TRUE))
en_all_barcodes <- colnames(en_obj)

en_obj <- cluster_subcluster(en_obj, output_dir = "./marker_genes/endothelial/", n_genes = 2000, features = 2000, dims = 1:20, conserve.memory = FALSE)
Idents(en_obj) <- "SCT_snn_res.0.1"
en_markers <- FindAllMarkers(en_obj, min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)

en_genes <- en_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 300) -> en_markers

write.csv(en_markers, "endo_marker_genes.csv")

endo_cluster_identification <- list(
  `0` = "Endo.General",            #PLVAP/PECAM1/VWF    
  `1` = "Endo.Vascular",           #VEGFA/PLVAP/PECAM1
  `2` = "Endo.Angiogenic",         #BICC1/PLVAP/PECAM1
  `3` = "Endo.Lymphatic"           #CCL21
)

new_ids <- unlist(endo_cluster_identification)
en_obj$fine_clust <- plyr::mapvalues(
  x    = as.character(Idents(en_obj)),
  from = names(new_ids),
  to   = new_ids
)

obj$fine_clust[colnames(en_obj)] <- en_obj$fine_clust

###############################################################################
# ASSIGN FINE GRAINED MARKER IDENTIFICATION FOR REMAINING CLUSTERS
###############################################################################

Idents(obj) <- "broad_cluster"
re_obj <- subset(obj, idents = grep("^Remaining", levels(Idents(obj)), value = TRUE))
re_all_barcodes <- colnames(re_obj)

re_obj <- cluster_subcluster(re_obj, output_dir = "./marker_genes/remaining/", n_genes = 2000, features = 2000, dims = 1:20, conserve.memory = FALSE)
Idents(re_obj) <- "SCT_snn_res.0.1"
re_markers <- FindAllMarkers(re_obj, min.pct = 0.1, logfc.threshold = 0.5, only.pos = TRUE)

re_genes <- re_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 300) -> re_markers

write.csv(re_markers, "remaining_marker_genes.csv")

re_cluster_identification <- list(
  `0` = "Melanocytes",            #MLANA/SOX5/DCT
  `1` = "Schwann.cells"           #NGFR/PMP2/SOX10
)

new_ids <- unlist(re_cluster_identification)
re_obj$fine_clust <- plyr::mapvalues(
  x    = as.character(Idents(re_obj)),
  from = names(new_ids),
  to   = new_ids
)

obj$fine_clust[colnames(re_obj)] <- re_obj$fine_clust
obj$mapping_cell_type <- sub("\\.\\d+$", "", obj$broad_cluster)

Idents(obj) <- "fine_clust"
saveRDS(obj, file = "./fine_annotated_obj.rds")

#################################################################
# SETUP PY ENVIRONMENT
#################################################################

# Please note that for GPU support you need to manually change
# parameters in setup_py_env.R Due to this being highly user 
# dependent, questions regarding setting up appropriate pytorch
# compatibility will not be supported. CellBender can run 
# without GPU support but this will take a very long time.

py_location <- "/home/uvictor/miniconda3/bin/conda"
conda_info_env <- setup_py_env(project, py_location)

obj <- readRDS(file = "./fine_annotated_obj.rds")

reticulate::py_install(
  packages = c("anndata", "scipy", "h5py"),
  pip = TRUE
)

reticulate::py_module_available("anndata")

as.anndata(x = obj, file_path = "./", file_name = "obj_anndata.h5ad")
