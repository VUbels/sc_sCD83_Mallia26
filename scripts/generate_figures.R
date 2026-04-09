library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(DESeq2)
library(tidyr)
library(speckle)

source("./scripts/helper_functions.R")


obj <- readRDS("./fine_annotated_obj.rds")

###########################################################################
# PLOT SAMPLE BALANCE
###########################################################################

print(colnames(obj@meta.data))
table(obj$orig.ident, obj$treatment, obj$mapping_cell_type)

# Check balance BEFORE running
show_sample_balance(
  obj,
  condition_col = "treatment",
  celltype_col = "mapping_cell_type",
  sample_col = "orig.ident",        
  subcluster_col = "fine_clust"
)

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
    Endothelial          = c("PECAM1","VWF","CDH5","CLDN5","ERG","DLL4"),
    Immune        = c("CD3D","CD3E","PTPRC","CD14","C1QA",
                      "TPSAB1","CPA3","LAMP3","CCR7","IDO1","IGKC","MZB1"),
    Other         = c("MLANA","DCT","TYRP1","SOX10","CDH19","PMP2")
  ),
  ident_col = "broad_cluster",
  bar_height = 1.2
)

p

ggsave("./marker_genes/grouped_dotplot.pdf", p, width = 10, height = 10)

####################################################################################
# PLOT BROAD MARKER GENES DOTPLOT
####################################################################################

# Run the improved pipeline
results <- run_full_de_pipeline(
  seurat_obj = obj,
  condition_col = "treatment",
  celltype_col = "mapping_cell_type",
  sample_col = "orig.ident",            
  ident.1 = "sCD83",
  ident.2 = "PBS",
  blacklist_genes = NULL,               
  clusters_to_remove = NULL,
  cluster_col = "fine_clust",
  cells_to_remove = NULL,
  top_n_genes = 2000,
  logfc_thresh = 0.5,
  pval_thresh = 0.05,
  pval_col = "p_val_BH",               
  run_pseudobulk = FALSE,                
  run_mast = TRUE,                      
  min_cells_per_sample = 10,
  min_samples_per_condition = 2,
  output_dir = "de_results"
)

# View results
print(results$plots$Global$volcano)
draw(results$plots$Global$heatmap)

# View concordance between methods
print(results$concordance_plots$Global)

# Access primary DE results
head(results$de_primary$global)

# View summary across methods
print(results$summary)

# gene_groups = list(Keratinocytes = c("KRT14","KRT17","KRT6A","S100A2","COL17A1","TP63","KRT85","KRT19"),
#                    Fibroblasts = c("COL1A1", "RGS5","MMP1","MMP3","VCAM1", "CRYAB", "TAGLN", "NFKB", "GAS1"),
#                    Endothelial = c("PECAM1","VWF","VEGFA","BICC1","CCL21","DLL4"),
#                    Immune = c("CD3D", "CD96", "KLRB1", "MRC1", "CD83", "KIT", "FOXP3", "CD8A", "LAMP3", "IGKC", "PAX5"),
#                    Remaining = c("MLANA","SOX5","DCT","NGFR","PMP2","SOX10"))


# plot_gene_bars_seurat(
#   seurat_obj = obj,
#   gene_groups = gene_groups,
#   large_cluster_col = "mapping_cell_type",
#   fine_cluster_col = "fine_clust",
#   outdir = "./marker_genes"
# )

# plot_umap_highlight_clusters(
#     seurat_obj = obj,
#     large_cluster_col = "mapping_cell_type",
#     pt_size = 2,
#     alpha_bg = 0.4,
#     outdir = "./marker_genes"
# )

missing_cells <- setdiff(Cells(full_obj), Cells(obj))

# Subset full_obj to just those cells
ors_r_obj <- subset(full_obj, cells = missing_cells)

# Set the metadata on the subset
ors_r_obj$fine_clust <- "ORS.R"
ors_r_obj$treatment  <- "PBS"
ors_r_obj$orig.ident <- "PBS_72h"

# Merge back into obj
obj <- merge(obj, ors_r_obj)

# Run enrichment test on full object
# ORS.2 is excluded from baseline computation but still tested and reported
# res <- test_cluster_enrichment(
#   obj              = obj,
#   treatment_col    = "treatment",
#   cluster_col      = "fine_clust",
#   sample_col       = "orig.ident",
#   condition_a      = "PBS",
#   condition_b      = "sCD83",
#   exclude_clusters = NULL,
#   alpha            = 0.001
# )
# 
# # View results table
# print(res, row.names = FALSE)
# 
# # Plot
# p <- plot_enrichment(
#   res,
#   title      = "PBS vs sCD83 cluster enrichment",
#   condA_color = "green4",
#   condB_color = "red3"
# )
# print(p)
# 
# # Save
# ggsave("cluster_enrichment.pdf", p, width = 14, height = 6)
