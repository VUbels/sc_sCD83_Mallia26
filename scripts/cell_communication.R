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
library(reticulate)

source("./scripts/helper_functions.R")

####################
# 2. PARAMETERS
####################

future::plan("sequential")

main_folder <- "./"
output_folder <- paste0("./cell_communication_results/")
dir.create(output_folder, showWarnings = FALSE)

#cellchat <- readRDS(paste0(output_folder, "cellchat_merged.rds"))
#cellchat.list <- readRDS(paste0(output_folder, "cellchat_list.rds"))

# Cell types to exclude from analysis (Population descrepancy to high between treatments)
exclude_cell_types <- c("ORS.1")

#######################
# CONDA SETUP
#######################

env_name <- "cell_communication"
options(reticulate.conda_binary = "/home/uvictor/miniconda3/bin/conda")

# Check existing environments
reticulate::conda_list()

# Create if doesn't exist, otherwise just use it
if (!(env_name %in% reticulate::conda_list()$name)) {
  # Environment doesn't exist - create it
  conda_create(env_name,
               python_version = "3.9",
               packages = c("pip", "umap-learn"))
  use_condaenv(paste0("/home/uvictor/miniconda3/envs/", env_name), required = TRUE)
} else {
  # Environment exists
  use_condaenv(paste0("/home/uvictor/miniconda3/envs/", env_name), required = TRUE)
}

# Verify configuration
py_config()


####################
# 3. DATA PREPARATION
####################


obj <- readRDS(paste0(main_folder, "fine_annotated_obj.rds"))

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
  CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor"))
  cellchat@DB <- CellChatDB.use
  
  # Preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute communication probabilities
  cellchat <- computeCommunProb(cellchat, 
                                type = "triMean", 
                                trim = 0.1,
                                nboot = 100)
  
  # Filter communications
  cellchat <- filterCommunication(cellchat, min.cells = 25)
  
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

pdf(paste0(output_folder, "pathway_ranking_stacked.pdf"), width = 4, height = 8)
print(gg1)
dev.off()

pdf(paste0(output_folder, "pathway_ranking.pdf"), width = 4, height = 8)
print(gg2)
dev.off()

####################
# 8. IDENTIFY SPECIFIC L-R PAIR CHANGES IN HFSCs/HFDSCs/
####################

cell_types <- levels(cellchat@idents$joint)
cell_types

# Increased signaling
result_increased <- netVisual_bubble(cellchat,
                                     sources.use = c(21),
                                     targets.use = c(2,9,10,11,12,13,17,20,26,27,28),
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
                                     sources.use = c(22,23),
                                     targets.use = c(2,9,10,11,12,13,17,20,26,27,28),
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

ggsave(paste0(output_folder, "LR_pair_changes_up_sCD83.png"), plot = gg_increased, width = 6, height = 10)
ggsave(paste0(output_folder, "LR_pair_changes_down_sCD83.png"), plot = gg_decreased, width = 7, height = 10)


pdf(paste0(output_folder, "LR_pair_changes_sCD83.pdf"), width = 10.5, height = 10)
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
exclude_pathways <- c("MIF", "LAMININ", "COLLAGEN", "VISFATIN")

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
pathways_of_interest <- unique(c(cellchat@netP$PBS$pathways, cellchat@netP$sCD83$pathways))

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

#############################################################################
# 10B. VIOLINPLOT SPECIFIC GENES
#############################################################################

pdf(paste0(output_folder, "/vln_plots/", "specific_genes_vlnplot.pdf"), width = 6, height = 6)

p <- plot_custom_expression(
  obj,
  features = c("SPP1", "CD44", "ITGB1", "ITGB3", "TGFB1", "TGFB3", "ACVR1B", "CXCL12", "AREG", "FGF2", "FGF7", "FGFR1"),
  idents = c("M2.Macrophages", "Bulge", "ORS.2", "ORS.3", "ORS.Suprabasal", "FBs.DP-like","FBs.Dermal.Sheath", "FBs.Cycling", "FBs.Activated", "Mast.cells"),
  split.by = "treatment"
)

print(p)

dev.off()

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
# 12. FUNCTIONAL ANALYSIS
####################

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

pdf(paste0(output_folder, "functional_pathway_analysis.pdf"), width = 11, height = 9)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
dev.off()

pdf(paste0(output_folder, "pathway_distance_analysis.pdf"), width = 11, height = 9)
rankSimilarity(cellchat, type = "functional")
dev.off()


########################################################
# FURTHER DIFFERENTIAL PATHWAYS
########################################################
pos.dataset = "sCD83"
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.fc = 0.05, thresh.pc = 0.1, thresh.p = 0.05, group.DE.combined = FALSE) 

net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat, net = net, datasets = "sCD83", ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "PBS", ligand.logFC = -0.01, receptor.logFC = NULL)

shared <- intersect(net.up$interaction_name, net.down$interaction_name)
net.up <- net.up[!net.up$interaction_name %in% shared, ]
net.down <- net.down[!net.down$interaction_name %in% shared, ]

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 21, targets.use = c(2,9,10,11,12,13,17,20,26,27,28), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(cellchat.list)[2]), show.legend = FALSE, font.size = 12)
#> Comparing communications on a merged object

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(22,23), targets.use = c(1:35), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(cellchat.list)[2]), font.size = 12)

gg1 <- clean_bubble_xaxis(gg1)
gg2 <- clean_bubble_xaxis(gg2)

# Comparing communications on a merged object
pdf(paste0(output_folder, "differential_pairs_analysis.pdf"), width = 12, height = 7)
wrap_plots(gg1, gg2, nrow = 1)
dev.off()

####################
# 13. SAVE OBJECTS
####################

saveRDS(cellchat, paste0(output_folder, "cellchat_merged.rds"))
saveRDS(cellchat.list, paste0(output_folder, "cellchat_list.rds"))
