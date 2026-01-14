#!/usr/bin/env Rscript

####################
# 1. LIBRARY LOADING
####################

library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(patchwork)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("glmGamPoi")


options(future.globals.maxSize = 2 * 1024^3)

####################
# 2. PARAMETERS
####################

input_folder <- "/mnt/d/scRNA_output/Mallia_25/Mallia_aRNA_corrected/"
output_folder <- "/mnt/d/scRNA_output/Mallia_25/"
objects <- list.files(input_folder, recursive = FALSE, include.dirs = FALSE, pattern = ".h5")
variables <- c("PBS_48h", "PBS_72h", "sCD83_48h", "sCD83_72h")

# Create output directory for DoubletFinder results
dir.create(paste0(output_folder, "DoubletFinder/"), showWarnings = FALSE)

####################
# 3. DATA LOADING AND INITIAL QC
####################

object.list <- list()

for (i in seq_along(objects)) {
  object <- objects[[i]]  
  stage <- variables[[i]]
  
  # Load CellBender-corrected data
  data.arna_corrected <- Read10X_h5(filename = paste0(input_folder, object), use.names = TRUE)  
  obj <- CreateSeuratObject(counts = data.arna_corrected, project = stage)
  obj$orig.ident <- stage
  
  # Calculate QC metrics
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # QC visualization
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  
  object.list[[i]] <- obj
  
  rm(data.arna_corrected)
  rm(obj)
}

####################
# 4. INITIAL QC FILTERING
####################

for (i in seq_along(object.list)) {
  obj <- object.list[[i]]
  
  # Standard QC thresholds - adjust based on your data
  obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 30)
  
  object.list[[i]] <- obj
  cat("Remaining cells after initial QC for", unique(obj$orig.ident), "is", ncol(obj), "cells\n")
}

####################
# 5. DOUBLETFINDER PIPELINE
####################

for (i in seq_along(object.list)) {
  obj <- object.list[[i]]
  
  # REMOVE any existing DoubletFinder columns (in case of re-running)
  obj@meta.data <- obj@meta.data[, !grepl("^pANN|^DF.classifications", colnames(obj@meta.data))]
  
  # DOUBLET RATE CALCULATION
  # Standard 10X Genomics doublet rate: ~0.8% per 1000 cells loaded
  doublet_assumption <- (ncol(obj)/1000) * 0.008
  
  cat("Processing:", unique(obj$orig.ident), "\n")
  cat("Cells:", ncol(obj), "\n")
  cat("Expected doublet rate:", round(doublet_assumption * 100, 2), "%\n")
  
  ####################
  # PREPROCESSING
  ####################
  
  obj <- SCTransform(obj, seed.use = 1)
  obj <- RunPCA(obj)
  
  # Use consistent PC dimensions throughout
  n_pcs <- 30
  
  obj <- RunUMAP(obj, dims = 1:n_pcs)
  obj <- FindNeighbors(obj, dims = 1:n_pcs, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 1.0)  # Higher resolution for better homotypic estimation
  
  # pK PARAMETER OPTIMIZATION
  
  cat("Running parameter sweep for pK optimization...\n")
  sweep.res.list <- paramSweep(obj, PCs = 1:n_pcs, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # Get the optimal pK value (highest BCmetric)
  pK_optimal <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  cat("Optimal pK:", pK_optimal, "\n")
  
  ####################
  # HOMOTYPIC DOUBLET PROPORTION ESTIMATE
  ####################
  
  annotations <- obj$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  # Calculate expected doublets
  nExp_poi <- round(doublet_assumption * nrow(obj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  cat("Expected doublets (raw):", nExp_poi, "\n")
  cat("Expected doublets (homotypic-adjusted):", nExp_poi.adj, "\n")
  
  pN <- 0.25
  
  ####################
  # RUN DOUBLETFINDER - FIRST PASS
  ####################
  
  cat("Running DoubletFinder (first pass)...\n")
  obj <- doubletFinder(obj, PCs = 1:n_pcs, pN = pN, pK = pK_optimal, 
                       nExp = nExp_poi, sct = TRUE)
  
  # CRITICAL: Get pANN column correctly (starts with "pANN", not "SCT_snn")
  pANN_col <- grep("^pANN", colnames(obj@meta.data), value = TRUE)
  pANN_col <- pANN_col[length(pANN_col)]
  
  cat("Using pANN column:", pANN_col, "\n")
  
  ####################
  # RUN DOUBLETFINDER - SECOND PASS (WITH HOMOTYPIC ADJUSTMENT)
  ####################
  
  cat("Running DoubletFinder (second pass with homotypic adjustment)...\n")
  obj <- doubletFinder(obj, PCs = 1:n_pcs, pN = pN, pK = pK_optimal, 
                       nExp = nExp_poi.adj, reuse.pANN = pANN_col, sct = TRUE)
  
  ####################
  # VISUALIZATION - pK OPTIMIZATION PLOT
  ####################
  
  pk_visualized <- ggplot(bcmvn, aes(x = as.numeric(as.character(pK)), y = BCmetric)) +
    geom_point(color = "#272E6A", size = 2) +
    geom_line(color = "#272E6A") +
    geom_vline(xintercept = pK_optimal, 
               color = "#D51F26", linetype = "dashed", size = 1) +
    geom_text(aes(x = pK_optimal, 
                  y = max(BCmetric), 
                  label = pK_optimal),
              vjust = -0.5, hjust = -0.2, color = "#D51F26", size = 3) +
    labs(x = "pK", y = "BCmetric", 
         title = paste0(unique(obj$orig.ident), " - pK Optimization")) +
    theme_classic() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12))
  
  figure_output <- paste0(output_folder, "DoubletFinder/", unique(obj$orig.ident), "_pKvalue.png")
  ggsave(filename = figure_output, plot = pk_visualized, width = 5, height = 4)
  
  ####################
  # VISUALIZATION - DOUBLET CLASSIFICATION
  ####################
  
  # Get the final classification column (most recent one)
  doublet_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)
  doublet_col <- doublet_col[length(doublet_col)]
  
  doublet_plot <- DimPlot(obj, group.by = doublet_col, pt.size = 0.5) +
    ggtitle(paste0(unique(obj$orig.ident), " - Doublet Classification"))
  
  figure_output <- paste0(output_folder, "DoubletFinder/", unique(obj$orig.ident), "_df_classification.png")
  ggsave(filename = figure_output, plot = doublet_plot, width = 6, height = 5)
  
  # Create standardized column name for downstream use
  obj@meta.data$doublet <- obj@meta.data[[doublet_col]]
  
  # Print summary
  doublet_table <- table(obj$doublet)
  cat("  Singlets:", doublet_table["Singlet"], "\n")
  cat("  Doublets:", doublet_table["Doublet"], "\n")
  cat("  Doublet %:", round(doublet_table["Doublet"] / ncol(obj) * 100, 2), "%\n")
  
  # UPDATE OBJECT LIST
  object.list[[i]] <- obj
  rm(obj)
  gc()
}

####################
# 6. FILTER DOUBLETS
####################

# Filter singlets only
object.list_filtered <- lapply(object.list, function(obj) {
  stage <- unique(obj$orig.ident)
  obj_filtered <- subset(obj, subset = doublet == "Singlet")
  cat("Remaining cells after doublet removal for", stage, ":", ncol(obj_filtered), "cells\n")
  return(obj_filtered)
})

####################
# 7. SAVE OUTPUT
####################

saveRDS(object.list, paste0(output_folder, "seurat_objects_with_doublets.rds"))
saveRDS(object.list_filtered, paste0(output_folder, "seurat_objects_doublet_filtered.rds"))

cat("Results saved to:", output_folder, "\n")

rm(object.list)
gc()