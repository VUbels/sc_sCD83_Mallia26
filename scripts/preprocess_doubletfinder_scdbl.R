#!/usr/bin/env Rscript

####################
# 1. LIBRARY LOADING
####################

library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = 2 * 1024^3)

####################
# 2. PARAMETERS
####################

input_folder <- "./arna_corrected/"
output_folder <- "./"


objects <- list.files(input_folder, recursive = FALSE, include.dirs = FALSE, pattern = ".h5")
variables <- c("PBS_48h", "PBS_72h", "sCD83_48h", "sCD83_72h")

# Create plot output directory
plot_folder <- "./arna_corrected/plots/"
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

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
  
  # QC violin plots — save and print
  vln <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(filename = paste0(plot_folder, stage, "_QC_violin.png"),
         plot = vln, width = 12, height = 5, dpi = 300)
  print(vln)
  
  # QC scatter plots — save and print
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  scatter <- plot1 + plot2
  ggsave(filename = paste0(plot_folder, stage, "_QC_scatter.png"),
         plot = scatter, width = 12, height = 5, dpi = 300)
  print(scatter)
  
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
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 30, nCount_RNA > 1000, nCount_RNA < 50000)
  
  object.list[[i]] <- obj
  cat("Remaining cells after initial QC for", unique(obj$orig.ident), "is", ncol(obj), "cells\n")
}

###############################################################
# 5. DOUBLET REMOVAL USING SCDBLFINDER
###############################################################

doublet_identification <- function(object.list, assumed_doublet_rate = 0.004) {
  for (i in seq_along(object.list)) {
    obj <- object.list[[i]]
    
    cat("\n###################################################\n")
    cat("Running scDblFinder on", unique(obj$orig.ident), "dataset ...\n")
    cat("###################################################\n")
    
    sce <- as.SingleCellExperiment(obj)
    
    set.seed(123)
    sce <- scDblFinder(sce, clusters = FALSE, dbr.per1k = assumed_doublet_rate)
    
    obj$scDblFinder.class <- sce$scDblFinder.class
    obj$scDblFinder.score <- sce$scDblFinder.score
    
    doublet_table <- table(sce$scDblFinder.class)
    cat("Doublets detected:", doublet_table["doublet"],
        "(", round(doublet_table["doublet"]/ncol(obj)*100, 2), "%)\n")
    cat("Singlets:", doublet_table["singlet"],
        "(", round(doublet_table["singlet"]/ncol(obj)*100, 2), "%)\n")
    
    object.list[[i]] <- obj
    rm(sce)
  }
  
  return(object.list)
}

###############################################################
# 6. VISUALIZATION OF DOUBLET DETECTION
###############################################################

doublet_visualization <- function(object.list, plot_folder) {
  
  dir.create(paste0(plot_folder, "doublet_detection/"), showWarnings = FALSE)
  
  for (i in seq_along(object.list)) {
    obj <- object.list[[i]]
    stage <- unique(obj$orig.ident)
    
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
    
    p1 <- DimPlot(obj, group.by = "scDblFinder.class", pt.size = 0.1) +
      ggtitle(paste(stage, "- Doublet Classification"))
    
    p2 <- FeaturePlot(obj, features = "scDblFinder.score", pt.size = 0.1) +
      ggtitle(paste(stage, "- Doublet Score"))
    
    combined_plot <- p1 + p2
    ggsave(filename = paste0(plot_folder, "doublet_detection/", stage, "_doublet_detection.png"),
           plot = combined_plot, width = 12, height = 5, dpi = 300)
    
    print(combined_plot)
    
    object.list[[i]] <- obj
  }
}

###############################################################
# 7. FILTER DOUBLETS
###############################################################

filter_doublets <- function(object.list) {
  object.list_filtered <- lapply(object.list, function(obj) {
    stage <- unique(obj$orig.ident)
    obj_filtered <- subset(obj, subset = scDblFinder.class == "singlet")
    cat("Remaining cells after doublet removal for", stage, ":", ncol(obj_filtered), "cells\n")
    return(obj_filtered)
  })
}

####################
# 8. EXECUTE PIPELINE
####################

object.list <- doublet_identification(object.list, assumed_doublet_rate = 0.004)
doublet_visualization(object.list, plot_folder)
object.list_filtered <- filter_doublets(object.list)

####################
# 9. SAVE OUTPUT
####################

saveRDS(object.list, paste0(output_folder, "obj_with_doublets.rds"))
saveRDS(object.list_filtered, paste0(output_folder, "obj_doublet_filtered.rds"))

cat("Results saved to:", output_folder, "\n")

rm(object.list)
gc()