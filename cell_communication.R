#!/usr/bin/env Rscript

####################
# 1. LIBRARY LOADING
####################

library(CellChat)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)
library(future)

####################
# 2. PARAMETERS
####################

future::plan("sequential")

main_folder <- "/mnt/d/scRNA_output/Mallia_25/"
output_folder <- paste0("./cell_communication_results/")
dir.create(output_folder, showWarnings = FALSE)

# Cell types to exclude from analysis (Population descrepancy to high between treatments)
exclude_cell_types <- c("ORS.1", "ORS.4")

####################
# 3. DATA PREPARATION
####################


obj <- readRDS(paste0(main_folder, "Mallia_25_LogN_Annotated.rds"))

if (length(exclude_cell_types) > 0) {
  cells_to_keep <- !obj$fine_clust %in% exclude_cell_types
  obj <- subset(obj, cells = which(cells_to_keep))
  cat("Removed", sum(!cells_to_keep), "cells from excluded cell types\n")
}


DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)

####################
# 4. CREATE CELLCHAT OBJECTS PER CONDITION
####################

# Split by treatment condition
conditions <- unique(obj$treatment)
cat("Treatments found:", paste(conditions, collapse = ", "), "\n")

cellchat.list <- list()

for (cond in conditions) {
  cat("Processing treatment:", cond, "\n")
  
  # Subset to condition
  obj_subset <- subset(obj, subset = treatment == cond)
  
  # Create CellChat object
  cellchat <- createCellChat(object = obj_subset, group.by = "fine_clust")
  
  # Set the ligand-receptor database
  CellChatDB <- CellChatDB.human
  
  # Use only Secreted Signaling for cleaner results (optional)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  
  # Preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probabilities
  cellchat <- computeCommunProb(cellchat, 
                                type = "truncatedMean", 
                                trim = 0.1,
                                nboot = 100)
  
  # Filter communications
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  # Compute pathway-level communication
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Aggregate network
  cellchat <- aggregateNet(cellchat)
  
  # Compute centrality
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  cellchat.list[[cond]] <- cellchat
  
  cat("Completed processing for", cond, "\n")
}

####################
# 5. MERGE CELLCHAT OBJECTS FOR COMPARISON
####################

cellchat <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))

####################
# 6. COMPARE NUMBER AND STRENGTH OF INTERACTIONS
####################

# Compare total number of interactions
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")

pdf(paste0(output_folder, "Interaction_comparison.pdf"), width = 10, height = 5)
print(gg1 + gg2)
dev.off()

####################
# 7. COMPARE INFORMATION FLOW BY PATHWAY
####################

# Pathway ranking comparison
gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)

pdf(paste0(output_folder, "pathway_ranking.pdf"), width = 12, height = 8)
print(gg1 + gg2)
dev.off()

####################
# 8. IDENTIFY SPECIFIC L-R PAIR CHANGES IN HFSCs/HFDSCs/
####################

cell_types <- levels(cellchat@idents$joint)
cell_types

# Increased signaling
result_increased <- netVisual_bubble(cellchat,
                                     sources.use = c(24),
                                     targets.use = c(9,10,14,15,16),
                                     comparison = c(1, 2), 
                                     max.dataset = 2, 
                                     title.name = "Increased signaling after sCD83", 
                                     angle.x = 45, 
                                     remove.isolate = TRUE,
                                     show.legend = FALSE,
                                     dot.size.min = 3.5,
                                     dot.size.max = 3.5,
                                     font.size.title = 15,
                                     return.data = TRUE
)

# Decreased signaling
result_decreased <- netVisual_bubble(cellchat,
                                     sources.use = c(24),
                                     targets.use = c(9,10,14,15,16),
                                     comparison = c(1, 2), 
                                     max.dataset = 1, 
                                     title.name = "Decreased signaling after sCD83", 
                                     angle.x = 45, 
                                     remove.isolate = TRUE,
                                     dot.size.min = 3.5,
                                     dot.size.max = 3.5,
                                     font.size.title = 15,
                                     return.data = TRUE
)

# Modify font sizes (adjust values as needed)
gg_increased <- result_increased$gg.obj + 
  theme(
    axis.text.x = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

gg_decreased <- result_decreased$gg.obj + 
  theme(
    axis.text.x = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.05, 'npc')
  )

pdf(paste0(output_folder, "LR_pair_changes_sCD83.pdf"), width = 15, height = 10)
print(gg_increased + gg_decreased)
dev.off()

# Extract data for downstream analysis
lr_increased_data <- result_increased$data
lr_decreased_data <- result_decreased$data

write.csv(lr_increased_data, paste0(output_folder, "LR_increased_sCD83_treatment.csv"), row.names = FALSE)
write.csv(lr_decreased_data, paste0(output_folder, "LR_decreased_sCD83_treatment.csv"), row.names = FALSE)

####################
# 9. SIGNALING PLOTS
####################

# Identify cell types with altered signaling patterns
# Exclude common large pathways that dominate the plot
exclude_pathways <- c("MIF", "LAMININ", "COLLAGEN")

pdf(paste0(output_folder, "Signaling_role_scatter.pdf"), width = 6, height = 6)
for (ct in cell_types) {
  tryCatch({
    gg <- netAnalysis_signalingChanges_scatter(cellchat, 
                                               idents.use = ct, 
                                               signaling.exclude = exclude_pathways,
                                               font.size.title = 20,
                                               font.size = 10
    ) +
      ggtitle(ct) +
      theme(aspect.ratio = 1) 
    print(gg)
  }, error = function(e) {
    cat("Could not plot", ct, ":", e$message, "\n")
  })
}
dev.off()

####################
# 10. PATHWAY-SPECIFIC ANALYSIS
####################

# Analyze specific pathways of interest
pathways_of_interest <- c("SPP1", "WNT", "PDGF", "BMP", "FGF", "EDN",
                          "EGF", "HGF", "IGF", "IL1", "IL2", "IL6", "IL10",
                          "IL16", "LIFR", "PTN")

dir.create(paste0(output_folder, "/circle_plots"), showWarnings = TRUE, recursive = FALSE)

for (i in seq_along(pathways_of_interest)) {
  
  pathway = pathways_of_interest[i]
  pdf(paste0(output_folder, "/circle_plots/", pathway, "_circleplot.pdf"), width = 8, height = 8)
  
    tryCatch({
      # Check if pathway exists
      pathway.union <- union(cellchat.list[[1]]@netP$pathways, 
                             cellchat.list[[2]]@netP$pathways)
      
      if (pathway %in% pathway.union) {
        # Circle plot comparison
        par(mfrow = c(1, 2))
        netVisual_aggregate(cellchat.list[[1]], signaling = pathway, 
                            layout = "circle", title = paste0(pathway, " - ", names(cellchat.list)[1]))
        netVisual_aggregate(cellchat.list[[2]], signaling = pathway, 
                            layout = "circle", title = paste0(pathway, " - ", names(cellchat.list)[2]))
      }
    }, error = function(e) {
      cat("Could not visualize", pathway, ":", e$message, "\n")
    })
  
  dev.off()

}

dir.create(paste0(output_folder, "/vln_plots"), showWarnings = TRUE, recursive = FALSE)

for (i in seq_along(pathways_of_interest)) {
  
  pathway = pathways_of_interest[i]
  pdf(paste0(output_folder, "/vln_plots/", pathway, "_vlnplot.pdf"), width = 8, height = 8)
  
  tryCatch({
    
    cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("PBS", "sCD83"))
    p <- plotGeneExpression(cellchat, signaling = pathway, split.by = "datasets", colors.ggplot = T, type = "violin")
    print(p)
    
    }
   , error = function(e) {
    cat("Could not visualize", pathway, ":", e$message, "\n")
  })
  
  dev.off()
  
}

dir.create(paste0(output_folder, "/heatmap_plots"), showWarnings = TRUE, recursive = FALSE)

for (i in seq_along(pathways_of_interest)) {
  
  pathway = pathways_of_interest[i]
  pdf(paste0(output_folder, "/heatmap_plots/", pathway, "_heatmap.pdf"), width = 12, height = 8)
  
  tryCatch({
    
    par(mfrow = c(1,2), xpd=TRUE)
    ht <- list()
    for (i in 1:length(cellchat.list)) {
      ht[[i]] <- netVisual_heatmap(cellchat.list[[i]], signaling = pathway, color.heatmap = "Reds",title.name = paste(pathway, "signaling ",names(cellchat.list)[i]))
    }
    ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
    
  }
  , error = function(e) {
    cat("Could not visualize", pathway, ":", e$message, "\n")
  })
  
  dev.off()
  
}

####################
# 11. OUTGOING/INCOMING SIGNALING HEATMAPS
####################

# Get union of pathways across conditions
pathway.union <- union(cellchat.list[[1]]@netP$pathways, cellchat.list[[2]]@netP$pathways)

# Outgoing signaling
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.list[[1]], 
                                         pattern = "outgoing", 
                                         signaling = pathway.union, 
                                         title = paste0(names(cellchat.list)[1]), 
                                         width = 11, height = 15,
                                         color.heatmap = "Blues")

ht2 <- netAnalysis_signalingRole_heatmap(cellchat.list[[2]], 
                                         pattern = "outgoing", 
                                         signaling = pathway.union, 
                                         title = paste0(names(cellchat.list)[2]), 
                                         width = 11, height = 15,
                                         color.heatmap = "Blues")

pdf(paste0(output_folder, "Outgoing_heatmaps.pdf"), width = 11, height = 9)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

# Incoming signaling
ht3 <- netAnalysis_signalingRole_heatmap(cellchat.list[[1]], 
                                         pattern = "incoming", 
                                         signaling = pathway.union, 
                                         title = paste0(names(cellchat.list)[1]), 
                                         width = 11, height = 15,
                                         color.heatmap = "Blues")

ht4 <- netAnalysis_signalingRole_heatmap(cellchat.list[[2]], 
                                         pattern = "incoming", 
                                         signaling = pathway.union, 
                                         title = paste0(names(cellchat.list)[2]), 
                                         width = 11, height = 15,
                                         color.heatmap = "Blues")

pdf(paste0(output_folder, "Incoming_heatmaps.pdf"), width = 11, height = 9)
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
dev.off()

####################
# 12. EXTRACT SPECIFIC COMMUNICATIONS
####################

# Functional analysis 

####################
# 13. SAVE OBJECTS
####################

saveRDS(cellchat, paste0(output_folder, "cellchat_merged.rds"))
saveRDS(cellchat.list, paste0(output_folder, "cellchat_list.rds"))
