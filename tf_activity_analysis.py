#!/usr/bin/env python3

import os
import re
import pickle
import torch

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad

from collections import Counter
from sklearn.preprocessing import minmax_scale
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from scipy.cluster.hierarchy import linkage, dendrogram

import warnings
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)
warnings.filterwarnings('ignore')

################
# CONFIGURATION
################

# Paths
main_folder = "/mnt/d/scrna_output/mallia_25/"
save_dir = os.path.join("./mallia_tfactivity")
os.makedirs(save_dir, exist_ok=True)

# Input files
adata_path = os.path.join(main_folder, "Mallia_25_LogN.h5ad")
grn_path = os.path.join(main_folder, "Greenleaf23_Skin_GRN_dataframe.parquet")

# Cell types to exclude (optional)
cells_to_remove = ['ORS.1', 'ORS.4', 'ORS.5', 'ORS.6'] #This is to prevent anomalies introduced due to severe population difference size between conditions.

# Weighting parameters
MIN_WEIGHT = 0.1        # Minimum weight for unexpressed TFs (protects niche populations)
MAX_WEIGHT = 1.0        # Maximum weight cap
SCALE_TO_MAX = 10.0     # Scale final weights to [0, 10] range (99th percentile normalization)

################
# COMMAND LINE ARGUMENTS
################

import argparse

def parse_arguments():
    """
    Parse command line arguments for model loading vs training.
    """
    parser = argparse.ArgumentParser(
        description='scRegulate: Weighted GRN + Differential TF + Ontology Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--load-model', '-l',
        type=str,
        default=None,
        metavar='PATH',
        help='Path to directory containing pre-trained model outputs.'
    )
    
    parser.add_argument(
        '--train', '-t',
        action='store_true',
        help='Train a new model (default behavior if --load-model not specified)'
    )
    
    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        default=None,
        metavar='PATH',
        help='Output directory for results (overrides config)'
    )
    
    # Individual file specification arguments
    parser.add_argument(
        '--adata-file',
        type=str,
        default="processed_adata.h5ad",
        metavar='FILENAME',
        help='Filename for processed AnnData (default: processed_adata.h5ad)'
    )
    
    parser.add_argument(
        '--tf-file',
        type=str,
        default="tf_activities.h5ad",
        metavar='FILENAME',
        help='Filename for TF activities (default: tf_activities.h5ad)'
    )
    
    parser.add_argument(
        '--grn-file',
        type=str,
        default="GRN.pkl",
        metavar='FILENAME',
        help='Filename for GRN pickle (default: GRN.pkl)'
    )
    
    parser.add_argument(
        '--model-file',
        type=str,
        default="fine_model.pt",
        metavar='FILENAME',
        help='Filename for PyTorch model (default: fine_model.pt)'
    )
    
    parser.add_argument(
        '--skip-enrichment',
        action='store_true',
        help='Skip ontology enrichment analysis'
    )
    
    parser.add_argument(
        '--skip-plots',
        action='store_true',
        help='Skip visualization generation'
    )
    
    return parser.parse_args()

# scRegulate training parameters
TRAINING_PARAMS = {
    'encoder_dims': [2048, 256, 64],  
    'z_dim': 40,                       
    'decoder_dims': [256],            
    'epochs': 10000,                   
    'batch_size': 1024,                
    'lr': 1e-4,                        
    'patience': 2500,                  
    'scheduler_patience': 200,         
    'min_lr': 1e-6,                    
    'min_targets': 10,                 
    'min_TFs': 5                       
}

# Fine-tuning parameters (two-stage)
FINETUNE_PARAMS = {
    'initial_finetune_epochs': 2000,
    'cluster_finetune_epochs_max': 10000,  
    'cluster_finetune_epochs_min': 5000,   
    'tf_mapping_lr': 4e-4,      
    'fc_output_lr': 2e-7,       
    'other_layers_lr': 3.5e-7   
}

# Differential analysis parameters
DIFF_PARAMS = {
    'condition_col': 'treatment',
    'condition_A': 'PBS',        
    'condition_B': 'sCD83',      
    'cluster_col': 'FineClust',
    'min_cells_per_condition': 3,  # ≥3 cells per condition required
    'p_threshold': 0.05,           # FDR α=0.05
    'min_abs_difference': 0.5      # |difference| ≥ 0.5 for significance
}

# Enrichment parameters
ENRICHMENT_PARAMS = {
    'databases': [
        'GO_Biological_Process_2023'  
    ],
    'qval_cutoff': 0.15,              
    'min_overlap': 2,                 
    'min_tfs_for_enrichment': 3,       
    'min_abs_difference_enrichment': 1.0  # |difference| ≥ 1.0 for enrichment selection
}

# Visualization parameters for differential TF heatmap
HEATMAP_PARAMS = {
    'n_top_per_celltype': 20,          # Top TFs per cell type
    'min_score_threshold': 30,         # Minimum combined score
    'selection_criterion': 'combined', # 'combined', 'fold_change', 'significance'
    'min_cell_types_appearance': 3,    # TF must appear in N cell types
    'high_score_threshold': 80,        # OR have score above this
    'final_top_n_tfs': 50,             # Final TFs for heatmap
    'cluster_rows': True,
    'cluster_columns': False,          # Use custom order instead
    'clustering_method': 'complete',
    'clustering_metric': 'correlation',
    'colormap': 'bwr',
    'show_significance_stars': True,
    'dpi_resolution': 330
}

# Custom cell type order for visualizations
CUSTOM_CELLTYPE_ORDER = [
    "Activated.EDs", "Arterial.EDs", "Capillary.EDs",
    "Lymphatic.EDs", "PID1.EDs", "Venous.EDs", "APM", "Pericytes",
    "DP-like.cells.1", "DP-like.cells.2", "Dermal.Sheath",
    "Fibroblasts.1", "Fibroblasts.2", "Fibroblasts.3", "HFDSCs",
    "IGFBP2.FBs", "Immune-associated.FBs",
    "Mast.cells", "Dendritic.cells", "M2.Macrophages", "NK.cells",
    "T.cells.1", "T.cells.2", "T.cyto", "T.reg", "T.residual", "B.cells",
    "HFSC.1", "HFSC.2", "IRS.1", "IRS.2",
    "Immune-associated.KCs", "Infindibulum", "Matrix", "ORS.2", "ORS.3",
    "Proliferating.KCs", "Sebaceous.Gland", "Eccrine.Gland", "Melanocytes",
    "Schwann.cells"
]

# Background colors for enrichment plots (range-based)
ENRICHMENT_RANGE_COLORS = {
    (0, 6): 'lightblue',      # Endothelial
    (6, 17): 'lightgreen',    # Fibroblasts 
    (17, 27): 'lightcoral',   # Immune
    (27, 37): 'lightyellow',  # Keratinocytes
    (37, 41): 'lightgray',    # Glandular/Neural-Crest 
}

print("scRegulate: WEIGHTED GRN + DIFFERENTIAL TF + ONTOLOGY PIPELINE")


################
# DATA PREPARATION
################

def prepare_data(adata, cells_to_remove=None):
    """
    Prepare data for scRegulate analysis.
    
    - Creates treatment column from cell index
    - Filters unwanted cell types
    - Extracts raw counts and normalizes
    """
    print("DATA PREPARATION")
    # Create treatment column from index if not present
    if 'treatment' not in adata.obs.columns:
        print("\nCreating treatment column from cell index...")
        adata.obs['treatment'] = adata.obs.index.str.split('_').str[0]
        adata.obs['treatment'] = adata.obs['treatment'].map({
            'pbs': 'PBS',
            'sCD83': 'sCD83'
        })
        print(f"  Treatment distribution: {adata.obs['treatment'].value_counts().to_dict()}")
    
    # Filter cell types if specified
    if cells_to_remove and len(cells_to_remove) > 0:
        print(f"\nRemoving cell types: {cells_to_remove}")
        print(f"  Before: {adata.shape[0]} cells")
        adata = adata[~adata.obs['FineClust'].isin(cells_to_remove)].copy()
        print(f"  After: {adata.shape[0]} cells")
    
    # Extract raw data and normalize
    print("\nExtracting and normalizing data...")
    if adata.raw is not None:
        raw_X = adata.raw.X.copy()
        obs = adata.obs.copy()
        var = adata.raw.var.copy()
        rna_data = sc.AnnData(X=raw_X, obs=obs, var=var)
        print("  Using raw counts from adata.raw")
    else:
        rna_data = adata.copy()
        print("  Using adata.X (no raw layer found)")
    
    # Normalize
    sc.pp.normalize_total(rna_data)
    sc.pp.log1p(rna_data)
    print("  Applied normalize_total and log1p")
    
    # Copy obsm if present
    if 'X_umap' in adata.obsm:
        rna_data.obsm['X_umap'] = adata.obsm['X_umap']
    
    print(f"\nPrepared data: {rna_data.shape}")
    
    return rna_data, adata


################
# BUILD WEIGHTED GRN
################

def build_weighted_grn(adata, skin_grn_raw, save_dir, 
                       min_weight=0.1, max_weight=1.0, scale_to_max=10.0):
    """
    Build expression-weighted GRN from Greenleaf ATAC data.
    
    Weighting logic:
    - TFs that are highly expressed get higher weights
    - TFs that are not expressed get minimum weight (not zero, to allow learning)
    - Log-scale transformation handles log-normal expression distribution
    """
    print("BUILDING WEIGHTED GRN")
    # Get TF columns
    tf_columns = [col for col in skin_grn_raw.columns 
                  if col not in ['peak_id', 'gene_short_name']]
    print(f"\nNumber of TFs in GRN: {len(tf_columns)}")
    
    # Extract expression matrix
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X
    
    gene_names = adata.var_names.tolist()
    
    # Calculate TF expression weights
    print("\nCalculating TF expression weights...")
    
    # Get mean expression per TF
    tf_expression = {}
    for tf in tf_columns:
        if tf in gene_names:
            tf_idx = gene_names.index(tf)
            tf_expression[tf] = X[:, tf_idx].mean()
        else:
            tf_expression[tf] = 0.0
    
    # Find max expression for normalization
    max_expr = max(tf_expression.values()) if tf_expression.values() else 1.0
    print(f"Max TF expression: {max_expr:.4f}")
    
    # Calculate weights using log-scale transformation
    tf_weights = {}
    for tf, expr in tf_expression.items():
        if expr <= 0:
            weight = min_weight
        else:
            # Log-scale weight
            weight = np.log1p(expr) / np.log1p(max_expr)
            weight = max(min_weight, min(max_weight, weight))
        tf_weights[tf] = weight
    
    # Show weight distribution
    weight_values = list(tf_weights.values())
    print(f"\nWeight distribution:")
    print(f"  Min:    {min(weight_values):.4f}")
    print(f"  25%:    {np.percentile(weight_values, 25):.4f}")
    print(f"  Median: {np.median(weight_values):.4f}")
    print(f"  75%:    {np.percentile(weight_values, 75):.4f}")
    print(f"  Max:    {max(weight_values):.4f}")
    
    # Show example TF weights
    print("\nExample TF weights:")
    example_tfs = ['MITF', 'SOX18', 'PAX5', 'TCF7', 'SOX9', 'MYF5', 'SOX10', 'TP63', 'KLF4']
    for tf in example_tfs:
        if tf in tf_weights:
            print(f"  {tf:10s}: expr={tf_expression[tf]:7.4f}, weight={tf_weights[tf]:.4f}")
    
    # Convert to long format with weights
    print("\nConverting to long format with expression weights...")
    tf_target_pairs = []
    
    for idx, row in skin_grn_raw.iterrows():
        target_gene = row['gene_short_name']
        peak_id = row['peak_id']
        
        for tf in tf_columns:
            if row[tf] == 1:
                tf_target_pairs.append({
                    'source': tf,
                    'target': target_gene,
                    'peak_id': peak_id,
                    'weight': tf_weights[tf]
                })
        
        if (idx + 1) % 10000 == 0:
            print(f"  Processed {idx + 1}/{len(skin_grn_raw)} peaks...")
    
    grn_long = pd.DataFrame(tf_target_pairs)
    print(f"\nTotal TF-peak-gene edges: {len(grn_long)}")
    
    # Aggregate by TF-target pairs
    print("\nAggregating by TF-target pairs...")
    grn_aggregated = grn_long.groupby(['source', 'target'], as_index=False).agg({
        'weight': 'sum',
        'peak_id': lambda x: ';'.join(x)
    })
    
    print(f"Unique TF-target pairs: {len(grn_aggregated)}")
    
    # Scale weights to [0, scale_to_max] range
    max_weight_99 = grn_aggregated['weight'].quantile(0.99)
    grn_aggregated['weight'] = (grn_aggregated['weight'] / max_weight_99) * scale_to_max
    grn_aggregated['weight'] = grn_aggregated['weight'].clip(upper=scale_to_max)
    
    print(f"\nFinal weight distribution (scaled to [0, {scale_to_max}]):")
    print(f"  Mean:   {grn_aggregated['weight'].mean():.4f}")
    print(f"  Median: {grn_aggregated['weight'].median():.4f}")
    print(f"  Max:    {grn_aggregated['weight'].max():.4f}")
    
    # Create CollecTRI-compatible format with PMID column
    net_weighted = grn_aggregated[['source', 'target', 'weight']].copy()
    net_weighted['PMID'] = 'CellOracle_' + grn_aggregated['peak_id'].astype(str)
    
    # Reorder columns to match CollecTRI format
    net_weighted = net_weighted[['source', 'target', 'weight', 'PMID']]
    
    # Filter TFs with ≥10 target genes as per methods
    min_targets = TRAINING_PARAMS.get('min_targets', 10)
    tf_target_counts = net_weighted.groupby('source').size()
    valid_tfs = tf_target_counts[tf_target_counts >= min_targets].index
    net_weighted = net_weighted[net_weighted['source'].isin(valid_tfs)]
    
    # Final statistics
    print(f"\nGRN Statistics (after filtering TFs with ≥{min_targets} targets):")
    print(f"  Weighted edges: {len(net_weighted):,}")
    print(f"  TFs: {net_weighted['source'].nunique():,}")
    print(f"  Target genes: {net_weighted['target'].nunique():,}")
    
    # Save both formats
    output_path = os.path.join(save_dir, "net_weighted_grn.parquet")
    net_weighted.to_parquet(output_path, index=False)
    print(f"\nSaved weighted GRN (CollecTRI format) to: {output_path}")
    
    # Also save without PMID for scRegulate (it only needs source, target, weight)
    net_for_scregulate = net_weighted[['source', 'target', 'weight']].copy()
    
    return net_for_scregulate, tf_weights


################
# TRAIN scRegulate MODEL
################

def train_scregulate_model(adata, net_weighted, save_dir, training_params, finetune_params):
    """
    Train scRegulate model with weighted GRN.
    
    Two-stage fine-tuning approach:
    1. Initial fine-tuning on all cells (2,000 epochs)
    2. Cluster-specific fine-tuning (5,000-10,000 epochs per cluster)
    
    Differential learning rates:
    - TF mapping layer: 4×10⁻⁴
    - Final output layer: 2×10⁻⁷  
    - All other layers: 3.5×10⁻⁷
    """
    print("TRAINING scRegulate MODEL")
    try:
        import scregulate as reg
    except ImportError:
        print("ERROR: scregulate not installed. Install with: pip install scregulate")
        return None, None
    
    # Prepare RNA data
    rna_data = adata.to_df()
    print(f"\nRNA data shape: {rna_data.shape}")
    print(f"GRN edges: {len(net_weighted)}")
    
    # Train model
    print("\nTraining scRegulate model...")
    print(f"Encoder: {training_params['encoder_dims']}, z_dim={training_params['z_dim']}")
    print(f"Decoder: {training_params['decoder_dims']}")
    print(f"Epochs: {training_params['epochs']}, Batch size: {training_params['batch_size']}")
    print(f"Learning rate: {training_params['lr']}, Patience: {training_params['patience']}")
    
    model, processed_adata, GRN = reg.train_model(
        rna_data=rna_data,
        net=net_weighted,
        encoder_dims=training_params['encoder_dims'],
        z_dim=training_params['z_dim'],
        decoder_dims=training_params['decoder_dims'],
        epochs=training_params['epochs'],
        batch_size=training_params['batch_size'],
        lr=training_params['lr'],
        patience=training_params['patience'],
        min_targets=training_params['min_targets'],
        min_TFs=training_params['min_TFs'],
        verbose=True
    )
    
    print("\nInitial training complete!")
    
    # Transfer metadata
    for col in adata.obs.columns:
        if col not in processed_adata.obs.columns:
            processed_adata.obs[col] = adata.obs[col].values
    
    # Stage 1: Initial fine-tuning on all cells
    print(f"\nStage 1: Initial fine-tuning ({finetune_params['initial_finetune_epochs']} epochs)...")
    
    processed_adata = reg.fine_tune_global(
        processed_adata=processed_adata,
        model=model,
        finetune_epochs=finetune_params['initial_finetune_epochs'],
        tf_mapping_lr=finetune_params['tf_mapping_lr'],
        fc_output_lr=finetune_params['fc_output_lr'],
        other_layers_lr=finetune_params['other_layers_lr'],
        verbose=True
    )
    
    # Stage 2: Cluster-specific fine-tuning
    print(f"\nStage 2: Cluster-specific fine-tuning...")
    print(f"  Epochs per cluster: {finetune_params['cluster_finetune_epochs_min']}-{finetune_params['cluster_finetune_epochs_max']}")
    print(f"  TF mapping LR: {finetune_params['tf_mapping_lr']}")
    print(f"  Output layer LR: {finetune_params['fc_output_lr']}")
    print(f"  Other layers LR: {finetune_params['other_layers_lr']}")
    
    cluster_col = DIFF_PARAMS['cluster_col']
    
    processed_adata = reg.fine_tune_clusters(
        processed_adata=processed_adata,
        model=model,
        cluster_column=cluster_col,
        finetune_epochs=finetune_params['cluster_finetune_epochs_max'],
        min_epochs=finetune_params['cluster_finetune_epochs_min'],
        tf_mapping_lr=finetune_params['tf_mapping_lr'],
        fc_output_lr=finetune_params['fc_output_lr'],
        other_layers_lr=finetune_params['other_layers_lr'],
        verbose=True
    )
    
    print("\nFine-tuning complete!")
    
    # Save model and data
    reg.save_model(model, os.path.join(save_dir, "scregulate_model.pt"))
    processed_adata.write_h5ad(os.path.join(save_dir, "processed_adata.h5ad"))
    
    return model, processed_adata


################
# DIFFERENTIAL TF ACTIVITY ANALYSIS
################

def differential_tf_activity(processed_adata, save_dir, diff_params):
    """
    Compute differential TF activity between conditions per cell type.
    
    Statistics computed per TF:
    - Mean activity in PBS and sCD83 cells
    - Difference in mean activity (sCD83_mean - PBS_mean)
    - Log2 fold change with pseudocount (0.0001)
    - Mann-Whitney U test p-values (two-sided)
    - Cohen's d effect size
    """
    print("DIFFERENTIAL TF ACTIVITY ANALYSIS")
    from scipy import stats
    
    condition_col = diff_params['condition_col']
    condition_A = diff_params['condition_A']
    condition_B = diff_params['condition_B']
    cluster_col = diff_params['cluster_col']
    p_threshold = diff_params['p_threshold']
    min_cells = diff_params['min_cells_per_condition']
    min_abs_diff = diff_params['min_abs_difference']
    
    # Pseudocount for log2 fold change
    PSEUDOCOUNT = 0.0001
    
    # Get TF activity matrix
    tf_activity = processed_adata.obsm['TF_activity']
    tf_names = processed_adata.uns['TF_names']
    
    print(f"\nTF activity matrix: {tf_activity.shape}")
    print(f"Conditions: {condition_A} vs {condition_B}")
    print(f"Cluster column: {cluster_col}")
    print(f"Min cells per condition: {min_cells}")
    
    # Get unique clusters
    clusters = processed_adata.obs[cluster_col].unique()
    print(f"Number of clusters: {len(clusters)}")
    
    results = []
    
    for cluster in clusters:
        cluster_mask = processed_adata.obs[cluster_col] == cluster
        cluster_data = tf_activity[cluster_mask]
        cluster_conditions = processed_adata.obs.loc[cluster_mask, condition_col]
        
        mask_A = cluster_conditions == condition_A
        mask_B = cluster_conditions == condition_B
        
        n_A = mask_A.sum()
        n_B = mask_B.sum()
        
        # Require ≥3 cells per condition as per methods
        if n_A < min_cells or n_B < min_cells:
            print(f"  Skipping {cluster}: insufficient cells (A={n_A}, B={n_B})")
            continue
        
        print(f"  Processing {cluster}: {n_A} {condition_A}, {n_B} {condition_B}")
        
        for i, tf in enumerate(tf_names):
            tf_A = cluster_data[mask_A, i]
            tf_B = cluster_data[mask_B, i]
            
            # Calculate means
            mean_A = tf_A.mean()
            mean_B = tf_B.mean()
            
            # Difference (B - A)
            difference = mean_B - mean_A
            
            # Log2 fold change with pseudocount
            log2fc = np.log2((mean_B + PSEUDOCOUNT) / (mean_A + PSEUDOCOUNT))
            
            # Wilcoxon rank-sum test (Mann-Whitney U, two-sided)
            try:
                stat, pval = stats.ranksums(tf_B, tf_A)
            except:
                pval = 1.0
            
            # Cohen's d effect size
            pooled_std = np.sqrt(((n_A - 1) * tf_A.std()**2 + (n_B - 1) * tf_B.std()**2) / (n_A + n_B - 2))
            if pooled_std > 0:
                cohens_d = difference / pooled_std
            else:
                cohens_d = 0.0
            
            results.append({
                'cell_type': cluster,
                'TF': tf,
                'mean_A': mean_A,
                'mean_B': mean_B,
                'difference': difference,
                'log2fc': log2fc,
                'pvalue': pval,
                'cohens_d': cohens_d,
                'n_A': n_A,
                'n_B': n_B
            })
    
    # Create DataFrame
    results_df = pd.DataFrame(results)
    
    # Multiple testing correction (Benjamini-Hochberg per cell type)
    from statsmodels.stats.multitest import multipletests
    
    results_df['p_adjusted'] = 1.0
    for cluster in results_df['cell_type'].unique():
        mask = results_df['cell_type'] == cluster
        pvals = results_df.loc[mask, 'pvalue'].values
        _, padj, _, _ = multipletests(pvals, method='fdr_bh')
        results_df.loc[mask, 'p_adjusted'] = padj
    
    # Mark significant (FDR p < 0.05 AND |difference| >= 0.5)
    results_df['significant'] = (
        (results_df['p_adjusted'] < p_threshold) & 
        (np.abs(results_df['difference']) >= min_abs_diff)
    )
    
    # Summary
    print(f"\nDifferential analysis complete!")
    print(f"  Total tests: {len(results_df)}")
    print(f"  Significant (FDR<{p_threshold}, |diff|>={min_abs_diff}): {results_df['significant'].sum()}")
    
    # Save
    output_path = os.path.join(save_dir, "differential_tf_activity_full.csv")
    results_df.to_csv(output_path, index=False)
    print(f"  Saved to: {output_path}")
    
    return results_df


################
# HELPER FUNCTIONS
################

def is_valid_gene_name(gene_name):
    """
    Check if gene name is a valid/known transcription factor symbol.
    Returns False for ENSG IDs, AC/AL contigs, LINC, etc.
    """
    gene_name = str(gene_name)
    
    invalid_patterns = [
        r'^ENSG\d+',           # Ensembl human genes
        r'^ENSMUSG\d+',        # Ensembl mouse genes
        r'^AC\d+\.\d+',        # Chromosome contigs
        r'^AL\d+\.\d+',        # Chromosome contigs
        r'^LINC\d+',           # Long intergenic non-coding RNAs
        r'^IRC\d+',            # IRC codes
        r'^LOC\d+',            # Uncharacterized loci
        r'^FP\d+',             # Predicted genes
        r'^RP\d+-\d+',         # Ribosomal pseudogenes
        r'^CTD-\d+',           # Clone-based gene names
        r'^KB-\d+',            # Clone-based gene names
        r'^\d+',               # Starts with number
    ]
    
    for pattern in invalid_patterns:
        if re.match(pattern, gene_name):
            return False
    
    if '.' in gene_name:
        return False
    
    if len(gene_name) < 2 or not any(c.isalpha() for c in gene_name):
        return False
    
    return True


def get_top_tfs_per_celltype_filtered(differential_results, tf_names, n_top=10, 
                                       criterion='combined', min_score_threshold=None):
    """
    Get top differentially active TFs per cell type with gene name filtering.
    """
    top_tfs_dict = {}
    
    for cell_type, results in differential_results.items():
        if criterion == 'combined':
            score = np.abs(results['fold_change']) * -np.log10(results['p_adjusted'] + 1e-300)
        elif criterion == 'fold_change':
            score = np.abs(results['fold_change'])
        elif criterion == 'significance':
            score = -np.log10(results['p_adjusted'] + 1e-300)
        elif criterion == 'effect_size':
            score = np.abs(results['difference'])
        
        top_indices = np.argsort(score)[-n_top:][::-1]
        
        if min_score_threshold is not None:
            top_indices = [i for i in top_indices if score[i] >= min_score_threshold]
        
        # Filter to valid gene names
        valid_tfs = [tf_names[i] for i in top_indices if is_valid_gene_name(tf_names[i])]
        
        # If filtered too many, get more candidates
        if len(valid_tfs) < n_top and len(top_indices) < len(tf_names):
            extended_indices = np.argsort(score)[-(n_top*3):][::-1]
            if min_score_threshold is not None:
                extended_indices = [i for i in extended_indices if score[i] >= min_score_threshold]
            valid_tfs = [tf_names[i] for i in extended_indices if is_valid_gene_name(tf_names[i])]
            valid_tfs = valid_tfs[:n_top]
        
        top_tfs_dict[cell_type] = valid_tfs
    
    return top_tfs_dict


################
# CLUSTER SIMILARITY ANALYSIS
################

def analyze_cluster_similarity(processed_adata, save_dir, subset_clusters=None):
    """
    Analyze similarity of regulatory networks between cell types using 
    cosine similarity of GRN posterior matrices.
    """
    print("CLUSTER SIMILARITY ANALYSIS")
    W_posteriors_per_cluster = processed_adata.uns.get("W_posteriors_per_cluster", None)
    
    if W_posteriors_per_cluster is None:
        print("No W_posteriors_per_cluster found - run fine-tuning first")
        return None
    
    cell_type_columns = list(W_posteriors_per_cluster.keys())
    
    # Build average W matrices
    print("\nBuilding average W matrices...")
    average_W_matrices = {
        cell_type: minmax_scale(np.abs(W_posteriors_per_cluster[cell_type]).ravel()).reshape(
            W_posteriors_per_cluster[cell_type].shape).mean(axis=0)
        for cell_type in cell_type_columns
    }
    
    # Combine into DataFrame
    combined_average_W = pd.DataFrame(average_W_matrices).T
    
    # Compute cosine similarity
    cosine_sim = np.clip(cosine_similarity(combined_average_W), 0, 1)
    similarity_matrix = cosine_sim**8  # Power transform for contrast
    
    similarity_df = pd.DataFrame(
        similarity_matrix, 
        index=cell_type_columns, 
        columns=cell_type_columns
    )
    
    # Plot full heatmap
    print("\nCreating similarity heatmap...")
    fig, ax = plt.subplots(figsize=(12, 12))
    sns_plot = sns.clustermap(
        similarity_df, 
        figsize=(12, 12),
        annot=False, 
        fmt=".1f", 
        annot_kws={"size": 8},
        cmap="RdBu_r",
        center=0,
        cbar=False
    )
    sns_plot.cax.set_visible(False)
    sns_plot.ax_heatmap.set_xlabel("", fontsize=16)
    sns_plot.ax_heatmap.set_ylabel("", fontsize=16)
    sns_plot.ax_heatmap.set_xticklabels(
        sns_plot.ax_heatmap.get_xticklabels(), fontsize=9, rotation=45, ha='right')
    sns_plot.ax_heatmap.set_yticklabels(
        sns_plot.ax_heatmap.get_yticklabels(), fontsize=9)
    
    output_path = os.path.join(save_dir, 'tf_grn_similarity_heatmap.png')
    sns_plot.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved to: {output_path}")
    
    # Subset analysis if requested
    if subset_clusters is not None:
        print(f"\nCreating subset heatmap for {len(subset_clusters)} clusters...")
        subset_W = {ct: average_W_matrices[ct] for ct in subset_clusters if ct in average_W_matrices}
        subset_combined_W = pd.DataFrame(subset_W).T
        subset_combined_W_rescaled = subset_combined_W.apply(lambda x: minmax_scale(x), axis=0)
        
        subset_cosine_sim = np.clip(cosine_similarity(subset_combined_W_rescaled), 0, 1)
        subset_similarity = subset_cosine_sim**8
        
        subset_similarity_df = pd.DataFrame(
            subset_similarity, 
            index=list(subset_W.keys()), 
            columns=list(subset_W.keys())
        )
        
        sns_plot_subset = sns.clustermap(
            subset_similarity_df, 
            figsize=(12, 10),
            annot=True, 
            fmt=".1f", 
            annot_kws={"size": 8},
            cmap="RdBu_r",
            center=0,
            cbar=False
        )
        sns_plot_subset.cax.set_visible(False)
        
        output_path_zoom = os.path.join(save_dir, 'tf_grn_similarity_heatmap_zoom.png')
        sns_plot_subset.savefig(output_path_zoom, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved subset to: {output_path_zoom}")
    
    return similarity_df


################
# DIFFERENTIAL TF HEATMAP
################

def create_differential_heatmap(processed_adata, save_dir, heatmap_params, 
                                 custom_order=None):
    """
    Create differential TF activity heatmap with significance stars.
    """
    print("CREATING DIFFERENTIAL TF HEATMAP")
    # Get TF activity matrix
    if 'TF_finetuned' in processed_adata.obsm:
        W_matrix = processed_adata.obsm['TF_finetuned']
    elif 'TF_activity' in processed_adata.obsm:
        W_matrix = processed_adata.obsm['TF_activity']
    else:
        print("No TF activity matrix found")
        return None
    
    if hasattr(W_matrix, 'toarray'):
        W_matrix = W_matrix.toarray()
    
    tf_names = processed_adata.uns['TF_names']
    conditions = processed_adata.obs['treatment'].values
    cell_types = processed_adata.obs['FineClust'].values
    cell_type_columns = np.unique(cell_types)
    
    # Perform differential analysis
    print("\nPerforming differential analysis...")
    differential_results = {}
    
    for cell_type in cell_type_columns:
        cell_type_mask = (cell_types == cell_type)
        
        pbs_mask = cell_type_mask & (conditions == 'PBS')
        scd83_mask = cell_type_mask & (conditions == 'sCD83')
        
        n_pbs = pbs_mask.sum()
        n_scd83 = scd83_mask.sum()
        
        if n_pbs >= 3 and n_scd83 >= 3:
            W_pbs_cells = W_matrix[pbs_mask, :]
            W_scd83_cells = W_matrix[scd83_mask, :]
            
            W_pbs_mean = W_pbs_cells.mean(axis=0)
            W_scd83_mean = W_scd83_cells.mean(axis=0)
            
            difference = W_scd83_mean - W_pbs_mean
            fold_change = np.log2((W_scd83_mean + 0.0001) / (W_pbs_mean + 0.0001))
            
            p_values = []
            for tf_idx in range(W_matrix.shape[1]):
                try:
                    stat, p_val = mannwhitneyu(
                        W_pbs_cells[:, tf_idx], 
                        W_scd83_cells[:, tf_idx], 
                        alternative='two-sided'
                    )
                    p_values.append(p_val)
                except:
                    p_values.append(1.0)
            
            p_values = np.array(p_values)
            reject, p_adjusted, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
            
            differential_results[cell_type] = {
                'PBS_mean': W_pbs_mean,
                'sCD83_mean': W_scd83_mean,
                'difference': difference,
                'fold_change': fold_change,
                'p_values': p_values,
                'p_adjusted': p_adjusted,
                'significant': reject
            }
    
    # Get top TFs
    print("\nSelecting top TFs...")
    top_tfs_per_celltype = get_top_tfs_per_celltype_filtered(
        differential_results, 
        tf_names, 
        n_top=heatmap_params['n_top_per_celltype'],
        criterion=heatmap_params['selection_criterion'],
        min_score_threshold=heatmap_params['min_score_threshold']
    )
    
    # Count frequency
    tf_frequency = Counter()
    tf_max_scores = {}
    
    for cell_type, tfs in top_tfs_per_celltype.items():
        results = differential_results[cell_type]
        for tf in tfs:
            tf_idx = tf_names.index(tf)
            tf_frequency[tf] += 1
            score = np.abs(results['fold_change'][tf_idx]) * -np.log10(results['p_adjusted'][tf_idx] + 1e-300)
            if tf not in tf_max_scores or score > tf_max_scores[tf]:
                tf_max_scores[tf] = score
    
    # Filter TFs
    filtered_tfs = [
        tf for tf in tf_frequency.keys() 
        if tf_frequency[tf] >= heatmap_params['min_cell_types_appearance'] or 
           tf_max_scores[tf] >= heatmap_params['high_score_threshold']
    ]
    
    # Select most variable
    analyzed_cell_types = [ct for ct in cell_type_columns if ct in differential_results]
    
    variance_scores = []
    for tf in filtered_tfs:
        tf_idx = tf_names.index(tf)
        tf_values = [differential_results[ct]['difference'][tf_idx] for ct in analyzed_cell_types]
        variance_scores.append(np.var(tf_values))
    
    if len(filtered_tfs) > heatmap_params['final_top_n_tfs']:
        most_variable_idx = np.argsort(variance_scores)[-heatmap_params['final_top_n_tfs']:]
        all_top_tfs = [filtered_tfs[i] for i in most_variable_idx]
    else:
        all_top_tfs = filtered_tfs
    
    print(f"Selected {len(all_top_tfs)} TFs for heatmap")
    
    # Build heatmap data
    heatmap_data_list = []
    for cell_type in analyzed_cell_types:
        results = differential_results[cell_type]
        row_data = []
        for tf in all_top_tfs:
            tf_idx = tf_names.index(tf)
            row_data.append(results['difference'][tf_idx])
        heatmap_data_list.append(row_data)
    
    heatmap_df = pd.DataFrame(heatmap_data_list, index=analyzed_cell_types, columns=all_top_tfs)
    heatmap_df_T = heatmap_df.T
    
    # Apply custom order if provided
    if custom_order is not None:
        ordered_cell_types = [ct for ct in custom_order if ct in heatmap_df_T.columns]
        heatmap_df_T = heatmap_df_T[ordered_cell_types]
    else:
        ordered_cell_types = list(heatmap_df_T.columns)
    
    # Cluster rows
    if heatmap_params['cluster_rows']:
        row_linkage = linkage(heatmap_df_T, method='complete', metric='correlation')
        row_dendrogram = dendrogram(row_linkage, no_plot=True)
        row_order = row_dendrogram['leaves']
        ordered_tfs = [all_top_tfs[i] for i in row_order]
        heatmap_df_T = heatmap_df_T.loc[ordered_tfs]
    else:
        ordered_tfs = all_top_tfs
    
    # Build significance matrix
    sig_matrix = np.zeros_like(heatmap_df_T.values, dtype=str)
    if heatmap_params['show_significance_stars']:
        for i, tf in enumerate(ordered_tfs):
            for j, cell_type in enumerate(ordered_cell_types):
                tf_idx = tf_names.index(tf)
                p_adj = differential_results[cell_type]['p_adjusted'][tf_idx]
                if p_adj < 0.001:
                    sig_matrix[i, j] = '***'
                elif p_adj < 0.01:
                    sig_matrix[i, j] = '**'
                elif p_adj < 0.05:
                    sig_matrix[i, j] = '*'
    
    # Plot
    fig_height = max(8, len(ordered_tfs) * 0.3)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    
    sns.heatmap(
        heatmap_df_T,
        cmap=heatmap_params['colormap'],
        center=0,
        robust=True,
        cbar_kws={'label': 'TF Activity Change\n(sCD83 - PBS)'},
        yticklabels=True,
        xticklabels=True,
        linewidths=0.5,
        linecolor='lightgray',
        ax=ax
    )
    
    # Add significance stars
    if heatmap_params['show_significance_stars']:
        for i in range(len(ordered_tfs)):
            for j in range(len(ordered_cell_types)):
                text = sig_matrix[i, j]
                if text:
                    ax.text(j + 0.5, i + 0.5, text,
                           ha='center', va='center',
                           color='black', fontsize=8, fontweight='bold')
    
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=10, rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=9)
    ax.set_title('Differential TF Activity (sCD83 vs PBS)', fontsize=14, pad=20)
    
    plt.tight_layout()
    
    output_path = os.path.join(save_dir, 'differential_tf_activity_heatmap.png')
    plt.savefig(output_path, dpi=heatmap_params['dpi_resolution'], bbox_inches='tight')
    plt.close()
    
    print(f"\nSaved heatmap to: {output_path}")
    
    return heatmap_df_T, differential_results


################
# MODEL SAVE/LOAD FUNCTIONS
################

def save_model_outputs(model, processed_adata, fine_tuned_tf_activities, GRN, save_dir, prefix=""):
    """
    Save all model outputs for later loading.
    """
    print("SAVING MODEL OUTPUTS")
    os.makedirs(save_dir, exist_ok=True)
    
    # Save PyTorch model
    model_path = os.path.join(save_dir, f"{prefix}model.pt")
    torch.save(model.state_dict(), model_path)
    print(f"Saved model to: {model_path}")
    
    # Convert W_posteriors_per_cluster keys to strings for h5ad compatibility
    if 'W_posteriors_per_cluster' in processed_adata.uns:
        processed_adata.uns['W_posteriors_per_cluster'] = {
            str(k): v for k, v in processed_adata.uns['W_posteriors_per_cluster'].items()
        }
    
    # Save processed AnnData
    adata_path = os.path.join(save_dir, f"{prefix}processed_adata.h5ad")
    processed_adata.write(adata_path)
    print(f"Saved processed_adata to: {adata_path}")
    
    # Save TF activities
    if fine_tuned_tf_activities is not None:
        tf_path = os.path.join(save_dir, f"{prefix}tf_activities.h5ad")
        fine_tuned_tf_activities.write(tf_path)
        print(f"Saved tf_activities to: {tf_path}")
    
    # Save GRN
    grn_path = os.path.join(save_dir, f"{prefix}GRN.pkl")
    with open(grn_path, "wb") as f:
        pickle.dump(GRN, f)
    print(f"Saved GRN to: {grn_path}")
    
    print("\nAll outputs saved successfully")


def load_model_outputs(load_dir, adata_file="processed_adata.h5ad", 
                       tf_file="tf_activities.h5ad", grn_file="GRN.pkl",
                       model_file="fine_model.pt"):
    """
    Load previously saved model outputs.
    
    Parameters:
    -----------
    load_dir : str
        Directory containing saved outputs
    adata_file : str
        Filename for processed AnnData
    tf_file : str
        Filename for TF activities
    grn_file : str
        Filename for GRN pickle
    model_file : str
        Filename for PyTorch model
    
    Returns:
    --------
    model, processed_adata, tf_activities, GRN
    """
    print("LOADING MODEL OUTPUTS")
    print(f"  Directory: {load_dir}")
    print(f"  AnnData:   {adata_file}")
    print(f"  TF file:   {tf_file}")
    print(f"  GRN file:  {grn_file}")
    print(f"  Model:     {model_file}")
    
    import scregulate as reg
    
    # Load processed AnnData
    adata_path = os.path.join(load_dir, adata_file)
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"AnnData file not found: {adata_path}")
    
    processed_adata = ad.read_h5ad(adata_path)
    print(f"\nLoaded processed_adata: {processed_adata.shape}")
    
    # Load TF activities (optional)
    tf_path = os.path.join(load_dir, tf_file)
    if os.path.exists(tf_path):
        tf_activities = ad.read_h5ad(tf_path)
        print(f"Loaded tf_activities: {tf_activities.shape}")
    else:
        tf_activities = None
        print(f"TF activities file not found (optional): {tf_file}")
    
    # Load GRN
    grn_path = os.path.join(load_dir, grn_file)
    if not os.path.exists(grn_path):
        raise FileNotFoundError(f"GRN file not found: {grn_path}")
    
    with open(grn_path, "rb") as f:
        GRN = pickle.load(f)
    
    # Handle different GRN formats
    if isinstance(GRN, dict):
        print(f"Loaded GRN: dict with {len(GRN)} keys")
        print(f"  Keys: {list(GRN.keys())[:10]}{'...' if len(GRN) > 10 else ''}")
        
        # Extract dimensions from dict or from processed_adata
        if 'matrix' in GRN:
            grn_matrix = GRN['matrix']
            input_dim = grn_matrix.shape[0]
            tf_dim = grn_matrix.shape[1]
        elif 'shape' in GRN:
            input_dim, tf_dim = GRN['shape']
        else:
            # Infer from processed_adata
            input_dim = processed_adata.n_vars
            if tf_activities is not None:
                tf_dim = tf_activities.n_vars
            elif 'TF_names' in processed_adata.uns:
                tf_dim = len(processed_adata.uns['TF_names'])
            else:
                raise ValueError("Cannot determine TF dimension from GRN dict. "
                                "Keys available: " + str(list(GRN.keys())))
    elif isinstance(GRN, (np.ndarray, pd.DataFrame)):
        if isinstance(GRN, pd.DataFrame):
            print(f"Loaded GRN: DataFrame {GRN.shape}")
        else:
            print(f"Loaded GRN: array {GRN.shape}")
        input_dim = GRN.shape[0]
        tf_dim = GRN.shape[1]
    else:
        print(f"Loaded GRN: {type(GRN)}")
        # Fallback to processed_adata dimensions
        input_dim = processed_adata.n_vars
        if tf_activities is not None:
            tf_dim = tf_activities.n_vars
        elif 'TF_names' in processed_adata.uns:
            tf_dim = len(processed_adata.uns['TF_names'])
        else:
            raise ValueError(f"Unknown GRN type: {type(GRN)} and cannot infer dimensions")
    
    print(f"  Model dimensions: input_dim={input_dim}, tf_dim={tf_dim}")
    
    # Set up modality structure
    processed_adata.modality = {}
    processed_adata.modality['RNA'] = processed_adata.copy()
    if tf_activities is not None:
        processed_adata.modality['TF'] = tf_activities.copy()
    
    # Copy TF activity to obsm if not already present
    if 'TF_activity' not in processed_adata.obsm:
        if tf_activities is not None:
            # TF activities stored in separate AnnData - copy X matrix to obsm
            if hasattr(tf_activities.X, 'toarray'):
                processed_adata.obsm['TF_activity'] = tf_activities.X.toarray()
            else:
                processed_adata.obsm['TF_activity'] = tf_activities.X.copy()
            print(f"Copied TF activity to obsm: {processed_adata.obsm['TF_activity'].shape}")
            
            # Also copy TF names if available
            if 'TF_names' not in processed_adata.uns:
                processed_adata.uns['TF_names'] = tf_activities.var_names.tolist()
                print(f"Copied TF names: {len(processed_adata.uns['TF_names'])} TFs")
        else:
            print("WARNING: No TF activity data found in obsm or separate file")
    else:
        print(f"TF_activity already in obsm: {processed_adata.obsm['TF_activity'].shape}")
    
    # Initialize and load model
    model = reg.scRNA_VAE(
        input_dim=input_dim,
        encode_dims=TRAINING_PARAMS['encoder_dims'],
        decode_dims=TRAINING_PARAMS['decoder_dims'],
        z_dim=TRAINING_PARAMS['z_dim'],
        tf_dim=tf_dim
    )
    
    model_path = os.path.join(load_dir, model_file)
    if os.path.exists(model_path):
        model.load_state_dict(torch.load(model_path, map_location='cpu'))
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = model.to(device)
        print(f"Loaded model to device: {device}")
    else:
        model = None
        print(f"Model file not found (optional): {model_file}")
        print("  Analysis can still proceed with processed_adata")
    
    return model, processed_adata, tf_activities, GRN


################
# ONTOLOGY ENRICHMENT ANALYSIS
################

def run_ontology_enrichment(results_df, save_dir, enrichment_params):
    """
    Run ontology enrichment on up/down regulated TFs per cell type.
    """
    print("ONTOLOGY ENRICHMENT ANALYSIS")
    try:
        import gseapy as gp
    except ImportError:
        print("ERROR: gseapy not installed. Install with: pip install gseapy")
        return None, None
    
    p_threshold = DIFF_PARAMS['p_threshold']
    min_abs_diff_enrichment = enrichment_params['min_abs_difference_enrichment']  # Stricter threshold for enrichment
    min_tfs = enrichment_params['min_tfs_for_enrichment']
    databases = enrichment_params['databases']
    qval_cutoff = enrichment_params['qval_cutoff']
    min_overlap = enrichment_params['min_overlap']
    
    # Separate UP and DOWN regulated TFs per cell type
    # Using stricter |difference| ≥ 1.0 threshold for enrichment as per methods
    upregulated_tfs = {}
    downregulated_tfs = {}
    
    print("\nSeparating UP vs DOWN regulated TFs per cell type...")
    print(f"Criteria: p_adjusted < {p_threshold}, |difference| >= {min_abs_diff_enrichment}")
    
    for cell_type in results_df['cell_type'].unique():
        ct_data = results_df[results_df['cell_type'] == cell_type]
        
        # Filter significant with STRICTER threshold for enrichment
        sig_up = ct_data[
            (ct_data['p_adjusted'] < p_threshold) &
            (ct_data['difference'] > 0) &
            (np.abs(ct_data['difference']) >= min_abs_diff_enrichment)
        ]
        
        sig_down = ct_data[
            (ct_data['p_adjusted'] < p_threshold) &
            (ct_data['difference'] < 0) &
            (np.abs(ct_data['difference']) >= min_abs_diff_enrichment)
        ]
        
        up_tfs = sig_up['TF'].tolist()
        down_tfs = sig_down['TF'].tolist()
        
        if len(up_tfs) >= min_tfs:
            upregulated_tfs[cell_type] = up_tfs
        if len(down_tfs) >= min_tfs:
            downregulated_tfs[cell_type] = down_tfs
    
    print(f"\nCell types with ≥{min_tfs} UP TFs: {len(upregulated_tfs)}")
    print(f"Cell types with ≥{min_tfs} DOWN TFs: {len(downregulated_tfs)}")
    
    # Run enrichment for UP-regulated TFs
    print("RUNNING ENRICHMENT FOR UP-REGULATED TFs")
    up_enrichment_results = {}
    
    for i, (cell_type, tf_list) in enumerate(upregulated_tfs.items(), 1):
        print(f"\n[{i}/{len(upregulated_tfs)}] {cell_type}: {len(tf_list)} UP TFs", end=' ')
        
        try:
            enr = gp.enrichr(
                gene_list=tf_list,
                gene_sets=databases,
                organism='human',
                outdir=None,
                no_plot=True,
                cutoff=qval_cutoff
            )
            
            if enr.results is not None and len(enr.results) > 0:
                results = enr.results.copy()
                results['cell_type'] = cell_type
                results['direction'] = 'UP'
                
                # Filter by significance
                results = results[
                    (results['Adjusted P-value'] < qval_cutoff) &
                    (results['Overlap'].str.split('/').str[0].astype(int) >= min_overlap)
                ]
                
                if len(results) > 0:
                    up_enrichment_results[cell_type] = results
                    print(f"{len(results)} terms")
                else:
                    print("No significant enrichment")
            else:
                print("No results")
                
        except Exception as e:
            print(f"Error: {str(e)[:50]}")
    
    # Run enrichment for DOWN-regulated TFs
    print("RUNNING ENRICHMENT FOR DOWN-REGULATED TFs")
    down_enrichment_results = {}
    
    for i, (cell_type, tf_list) in enumerate(downregulated_tfs.items(), 1):
        print(f"\n[{i}/{len(downregulated_tfs)}] {cell_type}: {len(tf_list)} DOWN TFs", end=' ')
        
        try:
            enr = gp.enrichr(
                gene_list=tf_list,
                gene_sets=databases,
                organism='human',
                outdir=None,
                no_plot=True,
                cutoff=qval_cutoff
            )
            
            if enr.results is not None and len(enr.results) > 0:
                results = enr.results.copy()
                results['cell_type'] = cell_type
                results['direction'] = 'DOWN'
                
                # Filter by significance
                results = results[
                    (results['Adjusted P-value'] < qval_cutoff) &
                    (results['Overlap'].str.split('/').str[0].astype(int) >= min_overlap)
                ]
                
                if len(results) > 0:
                    down_enrichment_results[cell_type] = results
                    print(f"{len(results)} terms")
                else:
                    print("No significant enrichment")
            else:
                print("No results")
                
        except Exception as e:
            print(f"Error: {str(e)[:50]}")
    
    # Combine results
    print("COMBINING RESULTS")
    all_results = []
    
    for cell_type, df in up_enrichment_results.items():
        all_results.append(df)
    for cell_type, df in down_enrichment_results.items():
        all_results.append(df)
    
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        
        # Save
        output_path = os.path.join(save_dir, "enrichment_up_down_all.csv")
        combined_df.to_csv(output_path, index=False)
        print(f"\nSaved {len(combined_df)} enriched terms to: {output_path}")
        
        return up_enrichment_results, down_enrichment_results
    else:
        print("\nNo enrichment results to save")
        return {}, {}


################
# VISUALIZATION
################

def create_enrichment_dotplot(up_results, down_results, save_dir,
                              max_terms=50, max_celltypes=35, 
                              custom_celltype_order=None,
                              range_colors=None,
                              show_all_celltypes=True):
    """
    Create directional enrichment dot plot showing UP (red) vs DOWN (blue) terms.
    
    Parameters:
    -----------
    up_results : dict
        UP-regulated enrichment results
    down_results : dict
        DOWN-regulated enrichment results
    save_dir : str
        Output directory
    max_terms : int
        Maximum number of terms to display
    max_celltypes : int
        Maximum number of cell types to display
    custom_celltype_order : list or None
        Custom order for cell types
    range_colors : dict or None
        Background colors by range, e.g., {(0, 6): 'lightblue', (6, 12): 'lightgreen'}
    show_all_celltypes : bool
        If True, show all cell types even if they have no enrichment data
    """
    print("CREATING ENRICHMENT VISUALIZATIONS")
    if not up_results and not down_results:
        print("No enrichment results to plot")
        return
    
    # Combine all results for plotting
    plot_data = []
    
    for cell_type, df in up_results.items():
        for _, row in df.iterrows():
            plot_data.append({
                'cell_type': cell_type,
                'term': row['Term'][:60],  # Truncate long terms
                'pval': row['Adjusted P-value'],
                'overlap': int(row['Overlap'].split('/')[0]),
                'direction': 'UP'
            })
    
    for cell_type, df in down_results.items():
        for _, row in df.iterrows():
            plot_data.append({
                'cell_type': cell_type,
                'term': row['Term'][:60],
                'pval': row['Adjusted P-value'],
                'overlap': int(row['Overlap'].split('/')[0]),
                'direction': 'DOWN'
            })
    
    if not plot_data:
        print("No data to plot")
        return
    
    plot_df = pd.DataFrame(plot_data)
    
    # Select top terms by frequency and significance
    term_scores = plot_df.groupby('term').agg({
        'pval': lambda x: -np.log10(x).mean(),
        'cell_type': 'count'
    }).reset_index()
    term_scores.columns = ['term', 'avg_significance', 'frequency']
    term_scores['combined_score'] = term_scores['avg_significance'] * term_scores['frequency']
    term_scores = term_scores.nlargest(max_terms, 'combined_score')
    selected_terms = term_scores['term'].tolist()
    
    # Filter plot data
    plot_df = plot_df[plot_df['term'].isin(selected_terms)]
    
    # Get cell types with custom order support
    if custom_celltype_order is not None:
        if show_all_celltypes:
            # Show ALL cell types from custom order
            celltypes = custom_celltype_order[:max_celltypes]
            celltypes_with_data = plot_df['cell_type'].unique()
            missing = [ct for ct in celltypes if ct not in celltypes_with_data]
            if missing:
                print(f"\nNote: {len(missing)} cell types have no enrichment data")
        else:
            celltypes = [ct for ct in custom_celltype_order if ct in plot_df['cell_type'].unique()]
            celltypes = celltypes[:max_celltypes]
    else:
        celltypes = sorted(plot_df['cell_type'].unique())[:max_celltypes]
    
    # Create figure
    fig_height = max(8, len(selected_terms) * 0.3)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    
    # Create position mappings
    celltype_pos = {ct: i for i, ct in enumerate(celltypes)}
    term_pos = {term: i for i, term in enumerate(selected_terms)}
    
    # Add background colors if specified
    if range_colors is not None:
        for (start, end), color in range_colors.items():
            ax.axvspan(start - 0.5, end - 0.5, alpha=0.15, color=color)
    
    # Plot dots
    for _, row in plot_df.iterrows():
        if row['cell_type'] not in celltype_pos:
            continue
        
        x = celltype_pos[row['cell_type']]
        y = term_pos[row['term']]
        
        # Color by direction
        color = 'darkred' if row['direction'] == 'UP' else 'darkblue'
        
        ax.scatter(
            x, y,
            s=row['overlap'] * 25,
            c=color,
            alpha=0.7,
            edgecolors='black',
            linewidth=0.5
        )
    
    # Customize axes
    ax.set_xticks(range(len(celltypes)))
    ax.set_xticklabels(celltypes, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(selected_terms)))
    ax.set_yticklabels(selected_terms, fontsize=8)
    
    ax.set_xlabel('Cell Type', fontsize=12, fontweight='bold')
    ax.set_ylabel('Enriched Term', fontsize=12, fontweight='bold')
    ax.set_title('Directional Pathway Enrichment\n(UP-regulated vs DOWN-regulated TFs)',
                fontsize=14, fontweight='bold', pad=20)
    
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Set axis limits
    ax.set_xlim(-0.5, len(celltypes) - 0.5)
    ax.set_ylim(-0.5, len(selected_terms) - 0.5)
    
    # Add legends
    from matplotlib.lines import Line2D
    
    # Size legend
    sizes = [3, 6, 9]
    size_legend = [plt.scatter([], [], s=s*25, c='gray', alpha=0.7, 
                               edgecolors='black', linewidth=0.5) for s in sizes]
    legend1 = ax.legend(size_legend, [f'{s} TFs' for s in sizes], 
                       title='TF Count', loc='upper left', 
                       bbox_to_anchor=(1.02, 1), frameon=True)
    
    # Direction legend
    up_marker = plt.scatter([], [], s=100, c='darkred', alpha=0.7,
                           edgecolors='black', linewidth=0.5, label='UP-regulated')
    down_marker = plt.scatter([], [], s=100, c='darkblue', alpha=0.7,
                             edgecolors='black', linewidth=0.5, label='DOWN-regulated')
    ax.legend(handles=[up_marker, down_marker], title='Direction',
             loc='upper left', bbox_to_anchor=(1.02, 0.7), frameon=True)
    ax.add_artist(legend1)
    
    plt.tight_layout()
    
    # Save
    output_path = os.path.join(save_dir, "enrichment_dotplot_directional.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"\nSaved dot plot to: {output_path}")


################
# MAIN EXECUTION
################

if __name__ == "__main__":
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Override output directory if specified
    if args.output_dir:
        save_dir = args.output_dir
        os.makedirs(save_dir, exist_ok=True)
        print(f"Output directory set to: {save_dir}")
    
    print("scRegulate: WEIGHTED GRN + DIFFERENTIAL TF + ONTOLOGY PIPELINE")
    
    ################
    # LOAD EXISTING MODEL (if --load-model specified)
    ################
    if args.load_model:
        print(f"\n{'='*60}")
        print(f"LOADING PRE-TRAINED MODEL FROM: {args.load_model}")
        print(f"{'='*60}")
        
        model, processed_adata, tf_activities, GRN = load_model_outputs(
            load_dir=args.load_model,
            adata_file=args.adata_file,
            tf_file=args.tf_file,
            grn_file=args.grn_file,
            model_file=args.model_file
        )
        
    ################
    # TRAIN NEW MODEL
    ################
    else:
        print(f"\n{'='*60}")
        print("TRAINING NEW MODEL")
        print(f"{'='*60}")
        
        # Load data
        print("\nLoading data...")
        adata = sc.read_h5ad(adata_path)
        print(f"Loaded adata: {adata.shape}")
        
        # Prepare data
        rna_data, adata_filtered = prepare_data(adata, cells_to_remove=cells_to_remove)
        
        # Load GRN
        skin_grn_raw = pd.read_parquet(grn_path)
        print(f"Loaded GRN: {skin_grn_raw.shape}")
        
        # Build weighted GRN
        net_weighted, tf_weights = build_weighted_grn(
            adata=adata_filtered,
            skin_grn_raw=skin_grn_raw,
            save_dir=save_dir,
            min_weight=MIN_WEIGHT,
            max_weight=MAX_WEIGHT,
            scale_to_max=SCALE_TO_MAX
        )
        
        if args.train:
            # Train scRegulate model
            model, processed_adata = train_scregulate_model(
                adata=rna_data,
                net_weighted=net_weighted,
                save_dir=save_dir,
                training_params=TRAINING_PARAMS,
                finetune_params=FINETUNE_PARAMS
            )
        else:
            # Load existing processed data (original behavior)
            processed_adata_path = os.path.join(save_dir, "processed_adata_Mallia.h5ad")
            if not os.path.exists(processed_adata_path):
                processed_adata_path = os.path.join(save_dir, "processed_adata.h5ad")
            
            if os.path.exists(processed_adata_path):
                print(f"\nLoading existing processed data from: {processed_adata_path}")
                processed_adata = sc.read_h5ad(processed_adata_path)
            else:
                print(f"\nNo processed data found. Use --train to train a new model.")
                exit(1)
    
    ################
    # DOWNSTREAM ANALYSIS
    ################
    print(f"\n{'='*60}")
    print("RUNNING DOWNSTREAM ANALYSIS")
    print(f"{'='*60}")
    
    # Cluster similarity analysis
    if not args.skip_plots:
        similarity_df = analyze_cluster_similarity(
            processed_adata=processed_adata,
            save_dir=save_dir,
            subset_clusters=CUSTOM_CELLTYPE_ORDER[:25]
        )
    
    # Differential TF activity analysis
    results_df = differential_tf_activity(
        processed_adata=processed_adata,
        save_dir=save_dir,
        diff_params=DIFF_PARAMS
    )
    
    # Create differential heatmap
    if not args.skip_plots:
        heatmap_df, diff_results = create_differential_heatmap(
            processed_adata=processed_adata,
            save_dir=save_dir,
            heatmap_params=HEATMAP_PARAMS,
            custom_order=CUSTOM_CELLTYPE_ORDER
        )
    
    # Ontology enrichment
    if not args.skip_enrichment:
        up_results, down_results = run_ontology_enrichment(
            results_df=results_df,
            save_dir=save_dir,
            enrichment_params=ENRICHMENT_PARAMS
        )
        
        # Visualization
        if not args.skip_plots:
            create_enrichment_dotplot(
                up_results=up_results,
                down_results=down_results,
                save_dir=save_dir,
                custom_celltype_order=CUSTOM_CELLTYPE_ORDER
            )
    
    print(f"\n{'='*60}")
    print("PIPELINE COMPLETE")
    print(f"{'='*60}")
