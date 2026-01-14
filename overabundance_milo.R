#!/usr/bin/env Rscript


####################
# 1. LIBRARY LOADING
####################

library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(scuttle)
library(irlba)
library(BiocParallel)
library(ggplot2)
library(sparseMatrixStats)
library(igraph) 

source("./helper_functions.R")

####################
# 2. PARAMETERS
####################

main_folder <- "/mnt/d/scRNA_output/Mallia_25/"
input_file <- paste0(main_folder, "Mallia_25_LogN_Annotated.rds")

# Create output directory
milo_dir <- file.path(main_folder, "Milo_results/")
dir.create(milo_dir, showWarnings = FALSE, recursive = TRUE)

####################
# 3. DATA LOADING
####################

# Load full Seurat object
obj_full <- readRDS(input_file)
DefaultAssay(obj_full) <- "RNA"

# Normalize if needed
obj_full <- NormalizeData(obj_full, normalization.method = "LogNormalize", scale.factor = 10000)

# Get unique cell types from mapping_cell_type
cell_types <- unique(obj_full@meta.data$mapping_cell_type)
cell_types <- cell_types[!is.na(cell_types)]  # Remove NAs if any

cat("Cell types to analyze:\n")
print(cell_types)

####################
# 4. LOOP THROUGH CELL TYPES
####################

for (cell_type in cell_types) {

  cat("Processing:", cell_type, "\n")

  cell_type_dir <- file.path(milo_dir, cell_type)
  dir.create(cell_type_dir, showWarnings = FALSE, recursive = TRUE)
  obj <- subset(obj_full, subset = mapping_cell_type == cell_type)
  
  cat("Number of cells:", ncol(obj), "\n")

  ####################
  # 5. MILOR SETUP
  ####################
  
  # Convert to SingleCellExperiment
  obj.sce <- as.SingleCellExperiment(obj)
  rm(obj)
  gc()
  
  # Set up parallel processing
  bpparam <- SerialParam()
  register(bpparam)
  set.seed(42)
  
  # Perform normalization if not already done
  if(!"logcounts" %in% assayNames(obj.sce)) {
    obj.sce <- logNormCounts(obj.sce)
  }
  
  # Feature selection
  dec <- modelGeneVar(obj.sce)
  hvgs <- getTopHVGs(dec, n = 2000)
  
  # PCA
  obj.sce <- runPCA(obj.sce, subset_row = hvgs, ncomponents = 50)
  
  # UMAP
  obj.sce <- runUMAP(obj.sce, dimred = "PCA", n_neighbors = 30, min_dist = 0.3)
  
  # Save UMAP plot
  png(file.path(cell_type_dir, "UMAP_treatment.png"), width = 8, height = 6, units = "in", res = 300)
  print(plotUMAP(obj.sce, colour_by = "treatment") + ggtitle(cell_type))
  dev.off()
  
  ####################
  # 6. CREATE MILO OBJECT
  ####################
  
  cat("Creating Milo object...\n")
  obj.milo <- Milo(obj.sce)
  
  # Build KNN graph
  obj.milo <- buildGraph(obj.milo, k = 50, d = 50, reduced.dim = "PCA")
  
  # Make neighborhoods
  obj.milo <- makeNhoods(obj.milo, prop = 0.5, k = 50, d = 50, 
                         refined = TRUE, reduced_dims = "PCA")
  
  # Save neighborhood size histogram
  png(file.path(cell_type_dir, "neighborhood_size_hist.png"), width = 6, height = 5, units = "in", res = 300)
  print(plotNhoodSizeHist(obj.milo) +
          labs(title = paste(cell_type, "- Neighborhood Sizes"),
               subtitle = paste("Total neighborhoods:", ncol(nhoods(obj.milo)))))
  dev.off()
  
  cat("Number of neighborhoods:", ncol(nhoods(obj.milo)), "\n")
  
  ####################
  # 7. COUNT CELLS AND CALC DISTANCES
  ####################
  
  cat("Counting cells in neighborhoods...\n")
  obj.milo <- countCells(obj.milo, 
                         meta.data = data.frame(colData(obj.milo)), 
                         samples = "orig.ident")
  
  cat("Calculating neighborhood distances...\n")
  obj.milo <- calcNhoodDistance(obj.milo, d = 30, reduced.dim = "PCA")
  
  ####################
  # 8. CREATE DESIGN MATRIX
  ####################
  
  sample_meta <- colData(obj.milo) %>%
    as.data.frame() %>%
    dplyr::select(orig.ident, treatment) %>%
    distinct()
  
  design_matrix <- data.frame(
    Sample = sample_meta$orig.ident,
    Treatment = sample_meta$treatment
  )
  
  rownames(design_matrix) <- design_matrix$Sample
  design_matrix <- design_matrix[colnames(nhoodCounts(obj.milo)), ]
  
  cat("Design matrix:\n")
  print(design_matrix)
  
  ####################
  # 9. TEST FOR DA
  ####################
  
  cat("Testing for differential abundance...\n")
  da_results <- testNhoods(obj.milo, 
                           design = ~ Treatment,
                           design.df = design_matrix)
  
  cat("Significant neighborhoods (P < 0.05):", sum(da_results$PValue < 0.05, na.rm = TRUE), "\n")
  cat("  - Upregulated:", sum(da_results$PValue < 0.05 & da_results$logFC > 0, na.rm = TRUE), "\n")
  cat("  - Downregulated:", sum(da_results$PValue < 0.05 & da_results$logFC < 0, na.rm = TRUE), "\n")
  
  ####################
  # 10. BUILD NHOOD GRAPH
  ####################
  
  cat("Building neighborhood graph...\n")
  obj.milo <- buildNhoodGraph(obj.milo)
  
  ####################
  # 11. VISUALIZATIONS
  ####################
  
  cat("Generating visualizations...\n")
  
  # P-value histogram
  png(file.path(cell_type_dir, "pvalue_histogram.png"), width = 6, height = 5, units = "in", res = 300)
  print(ggplot(da_results, aes(PValue)) + 
          geom_histogram(bins = 50, fill = "#272E6A", color = "white") +
          theme_classic() +
          labs(title = paste(cell_type, "- P-value Distribution"), 
               x = "P-value", y = "Count"))
  dev.off()
  
  # Volcano plot
  png(file.path(cell_type_dir, "volcano_plot.png"), width = 7, height = 6, units = "in", res = 300)
  print(ggplot(da_results, aes(logFC, -log10(PValue))) + 
          geom_point(aes(color = PValue < 0.05), size = 2) +
          scale_color_manual(values = c("grey60", "#D51F26"), 
                             labels = c("NS", "P < 0.05"),
                             name = "Significance") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
          theme_classic() +
          labs(title = paste(cell_type, "- Volcano Plot"), 
               x = "Log Fold Change", y = "-log10(P-value)"))
  dev.off()
  
  # Neighborhood graph with P-value filtering
  png(file.path(cell_type_dir, "nhood_graph_pval.png"), width = 8, height = 8, units = "in", res = 300)
  print(plotNhoodGraphDA_pval(obj.milo, da_results, alpha = 0.05, 
                              use_pvalue = TRUE, layout = "UMAP"))
  dev.off()
  
  ####################
  # 12. ANNOTATE BY FINE_CLUST
  ####################
  
  # Check if fine_clust exists
  if("fine_clust" %in% colnames(colData(obj.milo))) {
    cat("Annotating neighborhoods by fine_clust...\n")
    da_results <- annotateNhoods(obj.milo, da_results, coldata_col = "fine_clust")
    
    # Beeswarm plot
    png(file.path(cell_type_dir, "DA_beeswarm.png"), width = 10, height = 6, units = "in", res = 300)
    print(plotDAbeeswarm(da_results, group.by = "fine_clust") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                  axis.text.y = element_text(size = 10),
                  axis.title = element_text(size = 12)) +
            ggtitle(cell_type))
    dev.off()
  }
  
  ####################
  # 13. EXTRACT DA CELLS
  ####################
  
  cat("Extracting cells from significant neighborhoods...\n")
  
  # Upregulated neighborhoods
  cells_up <- extract_DA_cells(obj.milo, da_results, alpha = 0.05, 
                               direction = "up", use_pvalue = TRUE)
  cat("Cells in upregulated neighborhoods:", length(cells_up), "\n")
  
  # Downregulated neighborhoods
  cells_down <- extract_DA_cells(obj.milo, da_results, alpha = 0.05, 
                                 direction = "down", use_pvalue = TRUE)
  cat("Cells in downregulated neighborhoods:", length(cells_down), "\n")
  
  # Save cell lists
  write.csv(data.frame(barcode = cells_up), 
            file.path(cell_type_dir, "DA_cells_upregulated.csv"), row.names = FALSE)
  write.csv(data.frame(barcode = cells_down), 
            file.path(cell_type_dir, "DA_cells_downregulated.csv"), row.names = FALSE)
  
  ####################
  # 14. SUMMARY STATS
  ####################
  
  # Summary by fine_clust if available
  if("fine_clust" %in% colnames(da_results)) {
    summary_stats <- da_results %>%
      group_by(fine_clust) %>%
      summarise(
        n_nhoods = n(),
        n_sig_p0.05 = sum(PValue < 0.05, na.rm = TRUE),
        n_up = sum(PValue < 0.05 & logFC > 0, na.rm = TRUE),
        n_down = sum(PValue < 0.05 & logFC < 0, na.rm = TRUE),
        mean_logFC = mean(logFC, na.rm = TRUE),
        median_logFC = median(logFC, na.rm = TRUE)
      ) %>%
      arrange(desc(n_sig_p0.05))
    
    write.csv(summary_stats, file.path(cell_type_dir, "DA_summary_by_fineclust.csv"), row.names = FALSE)
  }
  
  ####################
  # 15. SAVE RESULTS
  ####################
  
  write.csv(da_results, file.path(cell_type_dir, "milo_da_results_full.csv"), row.names = FALSE)
  saveRDS(obj.milo, file.path(cell_type_dir, "milo_object.rds"))
  
  cat("Completed:", cell_type, "\n")
  
  # Clean up
  rm(obj.sce, obj.milo, da_results)
  gc()
}
