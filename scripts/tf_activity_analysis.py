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
main_folder = "./"
save_dir = os.path.join("./tfactivity")
os.makedirs(save_dir, exist_ok=True)

# Input files
adata_path = os.path.join(main_folder, "obj_anndata.h5ad")
grn_path = os.path.join(main_folder, "Skin_GRN_dataframe.parquet")

# Cell types to exclude (optional)
cells_to_remove = []

# GRN weight scaling (applied to peak multiplicity counts)
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
    
    parser.add_argument(
        '--adata-file',
        type=str,
        default="obj_annadata.h5ad",
        metavar='FILENAME',
        help='Filename for processed AnnData (default: obj_annadata.h5ad)'
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

    parser.add_argument(
        '--minimal', '-m',
        action='store_true',
        help=(
            'Minimal test-run mode: subsets to a small number of cells, genes, '
            'and GRN edges, and drastically reduces training epochs so the full '
            'pipeline can be validated quickly without GPU resources. '
            'Not intended for scientific results.'
        )
    )
    
    return parser.parse_args()

# Minimal mode: limits used when --minimal is passed.
#
# Only CELLS and EPOCHS are reduced. Gene set and GRN are intentionally
# kept full because scRegulate's ULM initialisation requires a sufficient
# gene-to-network overlap; subsampling either dimension collapses that
# overlap and produces NaN ULM estimates. The real training bottleneck
# is cells x epochs -- the GRN is a sparse prior and costs almost nothing.

MINIMAL_PARAMS = {
    'n_cells':          500,    # Cells to subsample (stratified by cell type)
    'epochs':           50,
    'freeze_epochs':    20,
    'batch_size':       128,
    'initial_finetune_epochs':          20,
    'cluster_finetune_epochs_max':      50,
    'cluster_finetune_epochs_min':      20,
    'early_stopping_patience':          50,
    'log_interval':                     10,
}


# scRegulate training parameters
TRAINING_PARAMS = {
    'encode_dims': [2048, 512, 128],
    'z_dim': 64,
    'decode_dims': [512, 2048],
    'epochs': 20000,
    'freeze_epochs': 12000,
    'batch_size': 2048,
    'learning_rate': 1.5e-4,
    'alpha_start': 0,
    'alpha_max': 0.8,
    'alpha_scale': 0.04,
    'beta_start': 0,
    'beta_max': 0.3,
    'gamma_start': 0,
    'gamma_max': 3.0,
    'log_interval': 200,
    'early_stopping_patience': 3000,
    'train_val_split_ratio': 0.9,
    'min_targets': 10,
    'min_TFs': 5
}

# Fine-tuning parameters (two-stage)
#   Stage 1: Global fine-tuning across all cells (no cluster_key)
#   Stage 2: Cluster-specific fine-tuning (with cluster_key)
FINETUNE_PARAMS = {
    'initial_finetune_epochs': 2000,        # Stage 1: global
    'cluster_finetune_epochs_max': 10000,   # Stage 2: per-cluster max
    'cluster_finetune_epochs_min': 5000,    # Stage 2: per-cluster min
    'tf_mapping_lr': 4e-4,      
    'fc_output_lr': 2e-7,       
    'other_layers_lr': 3.5e-7   
}

# Differential analysis parameters
DIFF_PARAMS = {
    'condition_col': 'treatment',
    'condition_A': 'PBS',        
    'condition_B': 'sCD83',      
    'cluster_col': 'fine_clust',
    'min_cells_per_condition': 3,
    'p_threshold': 0.05,
    'min_abs_difference': 0.5
}

# Enrichment parameters
ENRICHMENT_PARAMS = {
    'databases': [
        'GO_Biological_Process_2023'  
    ],
    'qval_cutoff': 0.15,              
    'min_overlap': 2,                 
    'min_tfs_for_enrichment': 3,       
    'min_abs_difference_enrichment': 1.0
}

# Visualization parameters for differential TF heatmap
HEATMAP_PARAMS = {
    'n_top_per_celltype': 20,
    'min_score_threshold': 30,
    'selection_criterion': 'combined',
    'min_cell_types_appearance': 3,
    'high_score_threshold': 80,
    'final_top_n_tfs': 50,
    'cluster_rows': True,
    'cluster_columns': False,
    'clustering_method': 'complete',
    'clustering_metric': 'correlation',
    'colormap': 'bwr',
    'show_significance_stars': True,
    'dpi_resolution': 330
}

# Custom cell type order for visualizations
CUSTOM_CELLTYPE_ORDER = [
    # Endothelial
    "Endo.General", "Endo.Vascular", "Endo.Angiogenic", "Endo.Lymphatic",
    # Fibroblasts / Stromal
    "FBs.Activated", "FBs.Interstitial", "FBs.Infl.Myofibroblast",
    "FBs.Pericytes", "FBs.DP-like", "FBs.Perifolicular",
    "FBs.Cycling", "FBs.NF-kB", "FBs.HF.Progenitor",
    "FBs.Dermal.Sheath", "FBs.Reticular",
    # Immune
    "Mast.cells", "Dendritic.cells", "M2.Macrophages",
    "MAIT.cells", "T.Naive-ANK3", "T.Naive",
    "T.cyto", "T.reg", "Plasma.cells", "B.cells",
    # Keratinocytes / HF
    "Bulge", "Ishtmus", "ORS.1", "ORS.2",
    "ORS.Suprabasal", "ORS.3", "Matrix", "Eccrine.cells",
    # Neural-crest
    "Melanocytes", "Schwann.cells"
]

# Background colors for enrichment plots (range-based)
ENRICHMENT_RANGE_COLORS = {
    (0, 4): 'lightblue',       # Endothelial
    (4, 15): 'lightgreen',     # Fibroblasts/Stromal
    (15, 26): 'lightcoral',    # Immune
    (26, 34): 'lightyellow',   # Keratinocytes/HF
    (34, 36): 'lightgray',     # Neural-crest
}

print("scRegulate: GRN + DIFFERENTIAL TF + ONTOLOGY PIPELINE")


################
# MINIMAL MODE HELPERS
################

def apply_minimal_overrides(training_params, finetune_params, minimal_params):
    """
    Return copies of training / fine-tune param dicts with epoch counts and
    batch size replaced by the fast smoke-test values from minimal_params.
    Original dicts are not mutated.
    """
    tp = training_params.copy()
    fp = finetune_params.copy()

    tp['epochs']                    = minimal_params['epochs']
    tp['freeze_epochs']             = minimal_params['freeze_epochs']
    tp['batch_size']                = minimal_params['batch_size']
    tp['early_stopping_patience']   = minimal_params['early_stopping_patience']
    tp['log_interval']              = minimal_params['log_interval']

    fp['initial_finetune_epochs']       = minimal_params['initial_finetune_epochs']
    fp['cluster_finetune_epochs_max']   = minimal_params['cluster_finetune_epochs_max']
    fp['cluster_finetune_epochs_min']   = minimal_params['cluster_finetune_epochs_min']

    return tp, fp


def subsample_cells(adata, n_cells, cluster_col='fine_clust', seed=0):
    """
    Stratified cell subsample for minimal smoke-test runs.

    Genes and GRN are intentionally left untouched: scRegulate's ULM
    initialisation needs a healthy gene-to-network overlap, which
    collapses if either the gene set or the GRN edges are thinned.
    The training bottleneck is cells x epochs, not GRN size.
    """
    if adata.n_obs <= n_cells:
        print(f"[MINIMAL]   Cell count ({adata.n_obs}) already <= {n_cells}; no cell filter applied")
        return adata.copy()

    print(f"\n[MINIMAL] Subsampling cells: {adata.n_obs} -> target {n_cells}")
    rng = np.random.default_rng(seed)

    if cluster_col in adata.obs.columns:
        clusters      = adata.obs[cluster_col].values
        unique_cl     = np.unique(clusters)
        n_per_cluster = max(3, n_cells // len(unique_cl))   # keep >=3 per cluster

        keep_idx = []
        for cl in unique_cl:
            idx = np.where(clusters == cl)[0]
            chosen = rng.choice(idx, size=min(n_per_cluster, len(idx)), replace=False)
            keep_idx.extend(chosen.tolist())

        if len(keep_idx) > n_cells:
            keep_idx = rng.choice(keep_idx, size=n_cells, replace=False).tolist()
    else:
        keep_idx = rng.choice(adata.n_obs, size=n_cells, replace=False).tolist()

    adata_sub = adata[keep_idx, :].copy()
    print(f"[MINIMAL]   Cells after subsample: {adata_sub.n_obs}")
    return adata_sub


################
# NORMALIZATION DETECTION
################

def detect_normalization_status(X, n_sample_cells=500, n_sample_genes=200):
    """
    Infer whether a count matrix has already been normalized / transformed,
    and recommend what (if anything) still needs to be applied.
    """
    import scipy.sparse as sp

    n_cells = X.shape[0]
    n_genes = X.shape[1]
    cell_idx = np.random.default_rng(42).choice(n_cells, size=min(n_sample_cells, n_cells), replace=False)
    gene_idx = np.random.default_rng(42).choice(n_genes, size=min(n_sample_genes, n_genes), replace=False)

    X_sample = X[np.ix_(cell_idx, gene_idx)]
    if sp.issparse(X_sample):
        X_sample = X_sample.toarray()
    X_sample = X_sample.astype(float)

    has_negative   = bool((X_sample < 0).any())
    global_max     = float(X_sample.max())
    global_min     = float(X_sample.min())
    global_mean    = float(X_sample.mean())
    global_std     = float(X_sample.std())

    nonzero        = X_sample[X_sample != 0]
    frac_integer   = float((np.abs(nonzero - np.round(nonzero)) < 1e-4).mean()) if len(nonzero) else 1.0
    looks_integer  = frac_integer > 0.99

    cell_sums      = X_sample.sum(axis=1)
    cell_sum_cv    = float(cell_sums.std() / (cell_sums.mean() + 1e-9))

    stats = {
        'has_negative':   has_negative,
        'global_min':     round(global_min, 4),
        'global_max':     round(global_max, 4),
        'global_mean':    round(global_mean, 4),
        'global_std':     round(global_std, 4),
        'frac_integer':   round(frac_integer, 4),
        'cell_sum_cv':    round(cell_sum_cv, 4),
    }

    if has_negative:
        return dict(
            status='sct_residuals', needs_norm=False, needs_log=False,
            reason="Negative values detected — data appears to be SCTransform Pearson residuals. Renormalization would corrupt the signal.",
            stats=stats
        )

    if not looks_integer and global_max < 30:
        return dict(
            status='log_normalized', needs_norm=False, needs_log=False,
            reason=f"Non-integer values with max={global_max:.2f} — data appears to be already log-normalized (SCTransform corrected counts or normalize_total + log1p). Skipping renormalization.",
            stats=stats
        )

    if looks_integer and global_max > 30:
        return dict(
            status='raw_counts', needs_norm=True, needs_log=True,
            reason=f"Integer-valued matrix with max={global_max:.0f} — data appears to be raw counts. Applying normalize_total + log1p.",
            stats=stats
        )

    return dict(
        status='ambiguous', needs_norm=True, needs_log=True,
        reason=f"Ambiguous matrix (non-integer={not looks_integer}, max={global_max:.2f}, cv={cell_sum_cv:.2f}). Applying normalize_total + log1p as a precaution.",
        stats=stats
    )


def print_normalization_report(norm_info, source_label):
    """Pretty-print the normalization detection report."""
    s = norm_info['stats']
    status_labels = {
        'raw_counts':    'RAW COUNTS',
        'log_normalized': 'ALREADY LOG-NORMALIZED',
        'sct_residuals': 'SCTransform RESIDUALS',
        'ambiguous':     'AMBIGUOUS — normalizing as precaution',
    }
    label = status_labels.get(norm_info['status'], norm_info['status'].upper())

    print(f"\n  Normalization check ({source_label}):")
    print(f"    Status       : {label}")
    print(f"    Has negatives: {s['has_negative']}")
    print(f"    Value range  : [{s['global_min']}, {s['global_max']}]")
    print(f"    Mean / Std   : {s['global_mean']} / {s['global_std']}")
    print(f"    Frac integer : {s['frac_integer']:.1%}")
    print(f"    Cell-sum CV  : {s['cell_sum_cv']:.3f}")
    print(f"    Decision     : {norm_info['reason']}")


################
# DATA PREPARATION
################

def prepare_data(adata, cells_to_remove=None):
    """
    Prepare data for scRegulate analysis.
    """
    print("DATA PREPARATION")
    if 'treatment' not in adata.obs.columns:
        print("\nCreating treatment column from cell index...")
        adata.obs['treatment'] = adata.obs.index.str.split('_').str[0]
        adata.obs['treatment'] = adata.obs['treatment'].map({
            'pbs': 'PBS',
            'sCD83': 'sCD83'
        })
        print(f"  Treatment distribution: {adata.obs['treatment'].value_counts().to_dict()}")

    if cells_to_remove and len(cells_to_remove) > 0:
        print(f"\nRemoving cell types: {cells_to_remove}")
        print(f"  Before: {adata.shape[0]} cells")
        adata = adata[~adata.obs['FineClust'].isin(cells_to_remove)].copy()
        print(f"  After: {adata.shape[0]} cells")

    print("\nExtracting and normalizing data...")

    if adata.raw is not None:
        raw_X = adata.raw.X.copy()
        obs   = adata.obs.copy()
        var   = adata.raw.var.copy()
        rna_data = sc.AnnData(X=raw_X, obs=obs, var=var)
        source_label = "adata.raw"
        print("  Source: adata.raw")
    else:
        rna_data = adata.copy()
        source_label = "adata.X"
        print("  Source: adata.X (no raw layer found)")

    norm_info = detect_normalization_status(rna_data.X)
    print_normalization_report(norm_info, source_label)

    if norm_info['needs_norm'] and norm_info['needs_log']:
        sc.pp.normalize_total(rna_data)
        sc.pp.log1p(rna_data)
        print("\n  Applied: normalize_total + log1p")
    elif norm_info['needs_log']:
        sc.pp.log1p(rna_data)
        print("\n  Applied: log1p only (already library-size normalized)")
    else:
        print("\n  Skipped normalization — data is already appropriately transformed")

    if 'X_umap' in adata.obsm:
        rna_data.obsm['X_umap'] = adata.obsm['X_umap']

    print(f"\nPrepared data: {rna_data.shape}")
    return rna_data, adata


################
# BUILD GRN (UNIFORM WEIGHTS — PEAK MULTIPLICITY ONLY)
################

def build_grn(skin_grn_raw, save_dir, scale_to_max=10.0):
    """
    Build GRN from Greenleaf ATAC data with uniform per-edge weights.
    """
    print("BUILDING GRN (UNIFORM WEIGHTS — PEAK MULTIPLICITY)")
    
    tf_columns = [col for col in skin_grn_raw.columns 
                  if col not in ['peak_id', 'gene_short_name']]
    print(f"\nNumber of TFs in GRN: {len(tf_columns)}")
    
    print("\nConverting to long format...")
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
                    'weight': 1.0
                })
        
        if (idx + 1) % 10000 == 0:
            print(f"  Processed {idx + 1}/{len(skin_grn_raw)} peaks...")
    
    grn_long = pd.DataFrame(tf_target_pairs)
    print(f"\nTotal TF-peak-gene edges: {len(grn_long)}")
    
    print("\nAggregating by TF-target pairs...")
    grn_aggregated = grn_long.groupby(['source', 'target'], as_index=False).agg({
        'weight': 'sum',
        'peak_id': lambda x: ';'.join(x)
    })
    
    print(f"Unique TF-target pairs: {len(grn_aggregated)}")
    
    max_weight_99 = grn_aggregated['weight'].quantile(0.99)
    grn_aggregated['weight'] = (grn_aggregated['weight'] / max_weight_99) * scale_to_max
    grn_aggregated['weight'] = grn_aggregated['weight'].clip(upper=scale_to_max)
    
    print(f"\nWeight distribution (peak multiplicity, scaled to [0, {scale_to_max}]):")
    print(f"  Mean:   {grn_aggregated['weight'].mean():.4f}")
    print(f"  Median: {grn_aggregated['weight'].median():.4f}")
    print(f"  75%:    {grn_aggregated['weight'].quantile(0.75):.4f}")
    print(f"  Max:    {grn_aggregated['weight'].max():.4f}")
    
    net = grn_aggregated[['source', 'target', 'weight']].copy()
    net['PMID'] = 'CellOracle_' + grn_aggregated['peak_id'].astype(str)
    net = net[['source', 'target', 'weight', 'PMID']]
    
    min_targets = TRAINING_PARAMS.get('min_targets', 10)
    tf_target_counts = net.groupby('source').size()
    valid_tfs = tf_target_counts[tf_target_counts >= min_targets].index
    net = net[net['source'].isin(valid_tfs)]
    
    print(f"\nGRN Statistics (after filtering TFs with >={min_targets} targets):")
    print(f"  Edges: {len(net):,}")
    print(f"  TFs: {net['source'].nunique():,}")
    print(f"  Target genes: {net['target'].nunique():,}")
    
    output_path = os.path.join(save_dir, "net_grn.parquet")
    net.to_parquet(output_path, index=False)
    print(f"\nSaved GRN (CollecTRI format) to: {output_path}")
    
    net_for_scregulate = net[['source', 'target', 'weight']].copy()
    
    return net_for_scregulate


################
# TRAIN scRegulate MODEL
################

def train_scregulate_model(adata, net, save_dir, training_params, finetune_params):
    """
    Train scRegulate model with GRN, then perform two-stage fine-tuning.

    Stage 1: Global fine-tuning across all cells (no cluster_key).
             This refines TF activities jointly across the full dataset.
    Stage 2: Cluster-specific fine-tuning (with cluster_key).
             This learns cell-type-resolved regulatory programmes.

    Both stages use reg.fine_tuning.fine_tune_clusters(), which returns
    four values: (processed_adata, fine_tuned_tf_activities, model, GRN).
    Omitting cluster_key triggers global mode (per official tutorial).
    """
    print("TRAINING scRegulate MODEL")
    try:
        import scregulate as reg
    except ImportError:
        print("ERROR: scregulate not installed. Install with: pip install scregulate")
        return None, None

    # ── Device selection ─────────────────────────────────────────────────
    if torch.cuda.is_available():
        device = torch.device("cuda")
        gpu_name = torch.cuda.get_device_name(0)
        print(f"\nGPU detected: {gpu_name}")
        print(f"  VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    else:
        device = torch.device("cpu")
        print("\nNo GPU detected — training on CPU (this will be slow)")
    print(f"Device: {device}")

    rna_data = adata
    print(f"\nRNA data shape: {rna_data.shape}")
    print(f"GRN edges: {len(net)}")
    
    print("\nTraining scRegulate model...")
    print(f"Encoder: {training_params['encode_dims']}, z_dim={training_params['z_dim']}")
    print(f"Decoder: {training_params['decode_dims']}")
    print(f"Epochs: {training_params['epochs']}, Batch size: {training_params['batch_size']}")
    print(f"Learning rate: {training_params['learning_rate']}, Patience: {training_params['early_stopping_patience']}")
    
    # ── Initial VAE training ─────────────────────────────────────────────
    model, processed_adata, GRN = reg.train_model(
        rna_data=rna_data,
        net=net,
        encode_dims=training_params['encode_dims'],
        z_dim=training_params['z_dim'],
        decode_dims=training_params['decode_dims'],
        epochs=training_params['epochs'],
        freeze_epochs=training_params['freeze_epochs'],
        batch_size=training_params['batch_size'],
        learning_rate=training_params['learning_rate'],
        alpha_start=training_params['alpha_start'],
        alpha_max=training_params['alpha_max'],
        alpha_scale=training_params['alpha_scale'],
        beta_start=training_params['beta_start'],
        beta_max=training_params['beta_max'],
        gamma_start=training_params['gamma_start'],
        gamma_max=training_params['gamma_max'],
        log_interval=training_params['log_interval'],
        early_stopping_patience=training_params['early_stopping_patience'],
        train_val_split_ratio=training_params['train_val_split_ratio'],
        min_targets=training_params['min_targets'],
        min_TFs=training_params['min_TFs'],
        device=device,
        verbose=True
    )
    
    print("\nInitial training complete!")
    
    # Transfer metadata from original adata to processed_adata
    for col in adata.obs.columns:
        if col not in processed_adata.obs.columns:
            processed_adata.obs[col] = adata.obs[col].values
    
    # ── Stage 1: Global fine-tuning (no cluster_key) ─────────────────────
    # Per the official scRegulate tutorial, calling fine_tune_clusters
    # without cluster_key performs global fine-tuning on all cells.
    print(f"\nStage 1: Global fine-tuning ({finetune_params['initial_finetune_epochs']} epochs)...")

    processed_adata, fine_tuned_tf_activities, fine_model, GRN_global = reg.fine_tuning.fine_tune_clusters(
        model=model,
        processed_adata=processed_adata,
        epochs=finetune_params['initial_finetune_epochs'],
        tf_mapping_lr=finetune_params['tf_mapping_lr'],
        fc_output_lr=finetune_params['fc_output_lr'],
        default_lr=finetune_params['other_layers_lr'],
        device=device,
        verbose=True
    )

    # ── Stage 2: Cluster-specific fine-tuning ────────────────────────────
    cluster_col = DIFF_PARAMS['cluster_col']
    print(f"\nStage 2: Cluster-specific fine-tuning (cluster_key='{cluster_col}')...")
    print(f"  Epochs per cluster: {finetune_params['cluster_finetune_epochs_min']}-{finetune_params['cluster_finetune_epochs_max']}")
    print(f"  TF mapping LR: {finetune_params['tf_mapping_lr']}")
    print(f"  Output layer LR: {finetune_params['fc_output_lr']}")
    print(f"  Other layers LR: {finetune_params['other_layers_lr']}")

    processed_adata, fine_tuned_tf_activities, final_model, GRN_final = reg.fine_tuning.fine_tune_clusters(
        processed_adata=processed_adata,
        model=fine_model,
        cluster_key=cluster_col,
        epochs=finetune_params['cluster_finetune_epochs_max'],
        min_epochs=finetune_params['cluster_finetune_epochs_min'],
        tf_mapping_lr=finetune_params['tf_mapping_lr'],
        fc_output_lr=finetune_params['fc_output_lr'],
        default_lr=finetune_params['other_layers_lr'],
        device=device,
        verbose=True
    )

    print("\nFine-tuning complete!")

    # ── Save all outputs ─────────────────────────────────────────────────
    # Model weights
    torch.save(final_model.state_dict(), os.path.join(save_dir, "scregulate_model.pt"))
    print(f"Saved model to: {os.path.join(save_dir, 'scregulate_model.pt')}")

    # Processed AnnData (convert W_posteriors keys to strings for h5ad)
    if 'W_posteriors_per_cluster' in processed_adata.uns:
        processed_adata.uns['W_posteriors_per_cluster'] = {
            str(k): v for k, v in processed_adata.uns['W_posteriors_per_cluster'].items()
        }
    processed_adata.write_h5ad(os.path.join(save_dir, "processed_adata.h5ad"))
    print(f"Saved processed_adata to: {os.path.join(save_dir, 'processed_adata.h5ad')}")

    # TF activities as separate AnnData
    fine_tuned_tf_activities.write_h5ad(os.path.join(save_dir, "tf_activities.h5ad"))
    print(f"Saved tf_activities to: {os.path.join(save_dir, 'tf_activities.h5ad')}")

    # GRN as pickle
    with open(os.path.join(save_dir, "GRN.pkl"), "wb") as f:
        pickle.dump(GRN_final, f)
    print(f"Saved GRN to: {os.path.join(save_dir, 'GRN.pkl')}")
    
    return final_model, processed_adata


################
# DIFFERENTIAL TF ACTIVITY ANALYSIS
################

def differential_tf_activity(processed_adata, save_dir, diff_params):
    """
    Compute differential TF activity between conditions per cell type.
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
    
    PSEUDOCOUNT = 0.0001
    
    tf_activity = processed_adata.obsm['TF_activity']
    tf_names = processed_adata.uns['TF_names']
    
    print(f"\nTF activity matrix: {tf_activity.shape}")
    print(f"Conditions: {condition_A} vs {condition_B}")
    print(f"Cluster column: {cluster_col}")
    print(f"Min cells per condition: {min_cells}")
    
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
        
        if n_A < min_cells or n_B < min_cells:
            print(f"  Skipping {cluster}: insufficient cells (A={n_A}, B={n_B})")
            continue
        
        print(f"  Processing {cluster}: {n_A} {condition_A}, {n_B} {condition_B}")
        
        for i, tf in enumerate(tf_names):
            tf_A = cluster_data[mask_A, i]
            tf_B = cluster_data[mask_B, i]
            
            mean_A = tf_A.mean()
            mean_B = tf_B.mean()
            difference = mean_B - mean_A
            log2fc = np.log2((mean_B + PSEUDOCOUNT) / (mean_A + PSEUDOCOUNT))
            
            try:
                stat, pval = stats.ranksums(tf_B, tf_A)
            except:
                pval = 1.0
            
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
    
    results_df = pd.DataFrame(results)
    
    results_df['p_adjusted'] = 1.0
    for cluster in results_df['cell_type'].unique():
        mask = results_df['cell_type'] == cluster
        pvals = results_df.loc[mask, 'pvalue'].values
        _, padj, _, _ = multipletests(pvals, method='fdr_bh')
        results_df.loc[mask, 'p_adjusted'] = padj
    
    results_df['significant'] = (
        (results_df['p_adjusted'] < p_threshold) & 
        (np.abs(results_df['difference']) >= min_abs_diff)
    )
    
    print(f"\nDifferential analysis complete!")
    print(f"  Total tests: {len(results_df)}")
    print(f"  Significant (FDR<{p_threshold}, |diff|>={min_abs_diff}): {results_df['significant'].sum()}")
    
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
        r'^ENSG\d+',
        r'^ENSMUSG\d+',
        r'^AC\d+\.\d+',
        r'^AL\d+\.\d+',
        r'^LINC\d+',
        r'^IRC\d+',
        r'^LOC\d+',
        r'^FP\d+',
        r'^RP\d+-\d+',
        r'^CTD-\d+',
        r'^KB-\d+',
        r'^\d+',
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
        
        valid_tfs = [tf_names[i] for i in top_indices if is_valid_gene_name(tf_names[i])]
        
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
    
    print("\nBuilding average W matrices...")
    average_W_matrices = {
        cell_type: minmax_scale(np.abs(W_posteriors_per_cluster[cell_type]).ravel()).reshape(
            W_posteriors_per_cluster[cell_type].shape).mean(axis=0)
        for cell_type in cell_type_columns
    }
    
    combined_average_W = pd.DataFrame(average_W_matrices).T
    
    cosine_sim = np.clip(cosine_similarity(combined_average_W), 0, 1)
    similarity_matrix = cosine_sim**8
    
    similarity_df = pd.DataFrame(
        similarity_matrix, 
        index=cell_type_columns, 
        columns=cell_type_columns
    )
    
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
    
    if subset_clusters is not None:
        print(f"\nCreating subset heatmap for {len(subset_clusters)} clusters...")
        subset_W = {ct: average_W_matrices[ct] for ct in subset_clusters if ct in average_W_matrices}
        if len(subset_W) >= 2:
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
        else:
            print(f"  Not enough matching clusters for subset heatmap ({len(subset_W)} found)")
    
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
    if 'TF_finetuned' in processed_adata.obsm:
        W_matrix = processed_adata.obsm['TF_finetuned']
    elif 'TF_activity' in processed_adata.obsm:
        W_matrix = processed_adata.obsm['TF_activity']
    else:
        print("No TF activity matrix found")
        return None, None
    
    if hasattr(W_matrix, 'toarray'):
        W_matrix = W_matrix.toarray()
    
    tf_names = processed_adata.uns['TF_names']
    conditions = processed_adata.obs['treatment'].values
    cell_types = processed_adata.obs['fine_clust'].values
    cell_type_columns = np.unique(cell_types)
    
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
    
    if not differential_results:
        print("No cell types with sufficient cells in both conditions")
        return None, None
    
    print("\nSelecting top TFs...")
    top_tfs_per_celltype = get_top_tfs_per_celltype_filtered(
        differential_results, 
        tf_names, 
        n_top=heatmap_params['n_top_per_celltype'],
        criterion=heatmap_params['selection_criterion'],
        min_score_threshold=heatmap_params['min_score_threshold']
    )
    
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
    
    filtered_tfs = [
        tf for tf in tf_frequency.keys() 
        if tf_frequency[tf] >= heatmap_params['min_cell_types_appearance'] or 
           tf_max_scores[tf] >= heatmap_params['high_score_threshold']
    ]
    
    if not filtered_tfs:
        print("No TFs passed filtering criteria")
        return None, None
    
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
    
    if custom_order is not None:
        ordered_cell_types = [ct for ct in custom_order if ct in heatmap_df_T.columns]
        heatmap_df_T = heatmap_df_T[ordered_cell_types]
    else:
        ordered_cell_types = list(heatmap_df_T.columns)
    
    if heatmap_params['cluster_rows'] and len(all_top_tfs) >= 2:
        row_linkage = linkage(heatmap_df_T, method='complete', metric='correlation')
        row_dendrogram = dendrogram(row_linkage, no_plot=True)
        row_order = row_dendrogram['leaves']
        ordered_tfs = [all_top_tfs[i] for i in row_order]
        heatmap_df_T = heatmap_df_T.loc[ordered_tfs]
    else:
        ordered_tfs = all_top_tfs
    
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
    
    model_path = os.path.join(save_dir, f"{prefix}model.pt")
    torch.save(model.state_dict(), model_path)
    print(f"Saved model to: {model_path}")
    
    if 'W_posteriors_per_cluster' in processed_adata.uns:
        processed_adata.uns['W_posteriors_per_cluster'] = {
            str(k): v for k, v in processed_adata.uns['W_posteriors_per_cluster'].items()
        }
    
    adata_path = os.path.join(save_dir, f"{prefix}processed_adata.h5ad")
    processed_adata.write(adata_path)
    print(f"Saved processed_adata to: {adata_path}")
    
    if fine_tuned_tf_activities is not None:
        tf_path = os.path.join(save_dir, f"{prefix}tf_activities.h5ad")
        fine_tuned_tf_activities.write(tf_path)
        print(f"Saved tf_activities to: {tf_path}")
    
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
    """
    print("LOADING MODEL OUTPUTS")
    print(f"  Directory: {load_dir}")
    print(f"  AnnData:   {adata_file}")
    print(f"  TF file:   {tf_file}")
    print(f"  GRN file:  {grn_file}")
    print(f"  Model:     {model_file}")
    
    import scregulate as reg
    
    adata_path = os.path.join(load_dir, adata_file)
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"AnnData file not found: {adata_path}")
    
    processed_adata = ad.read_h5ad(adata_path)
    print(f"\nLoaded processed_adata: {processed_adata.shape}")
    
    tf_path = os.path.join(load_dir, tf_file)
    if os.path.exists(tf_path):
        tf_activities = ad.read_h5ad(tf_path)
        print(f"Loaded tf_activities: {tf_activities.shape}")
    else:
        tf_activities = None
        print(f"TF activities file not found (optional): {tf_file}")
    
    grn_path = os.path.join(load_dir, grn_file)
    if not os.path.exists(grn_path):
        raise FileNotFoundError(f"GRN file not found: {grn_path}")
    
    with open(grn_path, "rb") as f:
        GRN = pickle.load(f)
    
    if isinstance(GRN, dict):
        print(f"Loaded GRN: dict with {len(GRN)} keys")
        if 'matrix' in GRN:
            input_dim = GRN['matrix'].shape[0]
            tf_dim = GRN['matrix'].shape[1]
        elif 'shape' in GRN:
            input_dim, tf_dim = GRN['shape']
        else:
            input_dim = processed_adata.n_vars
            tf_dim = len(processed_adata.uns.get('TF_names', []))
    elif isinstance(GRN, (np.ndarray, pd.DataFrame)):
        if isinstance(GRN, pd.DataFrame):
            print(f"Loaded GRN: DataFrame {GRN.shape}")
        else:
            print(f"Loaded GRN: array {GRN.shape}")
        input_dim = GRN.shape[0]
        tf_dim = GRN.shape[1]
    else:
        print(f"Loaded GRN: {type(GRN)}")
        input_dim = processed_adata.n_vars
        if tf_activities is not None:
            tf_dim = tf_activities.n_vars
        elif 'TF_names' in processed_adata.uns:
            tf_dim = len(processed_adata.uns['TF_names'])
        else:
            raise ValueError(f"Unknown GRN type: {type(GRN)} and cannot infer dimensions")
    
    print(f"  Model dimensions: input_dim={input_dim}, tf_dim={tf_dim}")
    
    # Set up modality structure (as per tutorial)
    processed_adata.modality = {}
    processed_adata.modality['RNA'] = processed_adata.copy()
    if tf_activities is not None:
        processed_adata.modality['TF'] = tf_activities.copy()
    
    if 'TF_activity' not in processed_adata.obsm:
        if tf_activities is not None:
            if hasattr(tf_activities.X, 'toarray'):
                processed_adata.obsm['TF_activity'] = tf_activities.X.toarray()
            else:
                processed_adata.obsm['TF_activity'] = tf_activities.X.copy()
            print(f"Copied TF activity to obsm: {processed_adata.obsm['TF_activity'].shape}")
            
            if 'TF_names' not in processed_adata.uns:
                processed_adata.uns['TF_names'] = tf_activities.var_names.tolist()
                print(f"Copied TF names: {len(processed_adata.uns['TF_names'])} TFs")
        else:
            print("WARNING: No TF activity data found in obsm or separate file")
    else:
        print(f"TF_activity already in obsm: {processed_adata.obsm['TF_activity'].shape}")
    
    # Reinitialize model with same architecture
    model = reg.scRNA_VAE(
        input_dim=input_dim,
        encode_dims=TRAINING_PARAMS['encode_dims'],
        decode_dims=TRAINING_PARAMS['decode_dims'],
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
    min_abs_diff_enrichment = enrichment_params['min_abs_difference_enrichment']
    min_tfs = enrichment_params['min_tfs_for_enrichment']
    databases = enrichment_params['databases']
    qval_cutoff = enrichment_params['qval_cutoff']
    min_overlap = enrichment_params['min_overlap']
    
    upregulated_tfs = {}
    downregulated_tfs = {}
    
    print("\nSeparating UP vs DOWN regulated TFs per cell type...")
    print(f"Criteria: p_adjusted < {p_threshold}, |difference| >= {min_abs_diff_enrichment}")
    
    for cell_type in results_df['cell_type'].unique():
        ct_data = results_df[results_df['cell_type'] == cell_type]
        
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
    
    print(f"\nCell types with >={min_tfs} UP TFs: {len(upregulated_tfs)}")
    print(f"Cell types with >={min_tfs} DOWN TFs: {len(downregulated_tfs)}")
    
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
    
    print("COMBINING RESULTS")
    all_results = []
    
    for cell_type, df in up_enrichment_results.items():
        all_results.append(df)
    for cell_type, df in down_enrichment_results.items():
        all_results.append(df)
    
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        
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
    """
    print("CREATING ENRICHMENT VISUALIZATIONS")
    if not up_results and not down_results:
        print("No enrichment results to plot")
        return
    
    plot_data = []
    
    for cell_type, df in up_results.items():
        for _, row in df.iterrows():
            plot_data.append({
                'cell_type': cell_type,
                'term': row['Term'][:60],
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
    
    term_scores = plot_df.groupby('term').agg({
        'pval': lambda x: -np.log10(x).mean(),
        'cell_type': 'count'
    }).reset_index()
    term_scores.columns = ['term', 'avg_significance', 'frequency']
    term_scores['combined_score'] = term_scores['avg_significance'] * term_scores['frequency']
    term_scores = term_scores.nlargest(max_terms, 'combined_score')
    selected_terms = term_scores['term'].tolist()
    
    plot_df = plot_df[plot_df['term'].isin(selected_terms)]
    
    if custom_celltype_order is not None:
        if show_all_celltypes:
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
    
    fig_height = max(8, len(selected_terms) * 0.3)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    
    celltype_pos = {ct: i for i, ct in enumerate(celltypes)}
    term_pos = {term: i for i, term in enumerate(selected_terms)}
    
    if range_colors is not None:
        for (start, end), color in range_colors.items():
            ax.axvspan(start - 0.5, end - 0.5, alpha=0.15, color=color)
    
    for _, row in plot_df.iterrows():
        if row['cell_type'] not in celltype_pos:
            continue
        
        x = celltype_pos[row['cell_type']]
        y = term_pos[row['term']]
        
        color = 'darkred' if row['direction'] == 'UP' else 'darkblue'
        
        ax.scatter(
            x, y,
            s=row['overlap'] * 25,
            c=color,
            alpha=0.7,
            edgecolors='black',
            linewidth=0.5
        )
    
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
    
    ax.set_xlim(-0.5, len(celltypes) - 0.5)
    ax.set_ylim(-0.5, len(selected_terms) - 0.5)
    
    from matplotlib.lines import Line2D
    
    sizes = [3, 6, 9]
    size_legend = [plt.scatter([], [], s=s*25, c='gray', alpha=0.7, 
                               edgecolors='black', linewidth=0.5) for s in sizes]
    legend1 = ax.legend(size_legend, [f'{s} TFs' for s in sizes], 
                       title='TF Count', loc='upper left', 
                       bbox_to_anchor=(1.02, 1), frameon=True)
    
    up_marker = plt.scatter([], [], s=100, c='darkred', alpha=0.7,
                           edgecolors='black', linewidth=0.5, label='UP-regulated')
    down_marker = plt.scatter([], [], s=100, c='darkblue', alpha=0.7,
                             edgecolors='black', linewidth=0.5, label='DOWN-regulated')
    ax.legend(handles=[up_marker, down_marker], title='Direction',
             loc='upper left', bbox_to_anchor=(1.02, 0.7), frameon=True)
    ax.add_artist(legend1)
    
    plt.tight_layout()
    
    output_path = os.path.join(save_dir, "enrichment_dotplot_directional.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(output_path.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"\nSaved dot plot to: {output_path}")


################
# MAIN EXECUTION
################

if __name__ == "__main__":
    
    args = parse_arguments()

    # ── Minimal-mode banner ──────────────────────────────────────────────
    if args.minimal:
        print("\n" + "="*60)
        print("  MINIMAL TEST-RUN MODE  (-m / --minimal)")
        print("  Dataset and epochs severely reduced for smoke-testing.")
        print("  Results are NOT scientifically meaningful.")
        print("="*60 + "\n")

    if args.output_dir:
        save_dir = args.output_dir
        os.makedirs(save_dir, exist_ok=True)
        print(f"Output directory set to: {save_dir}")
    
    print("scRegulate: GRN + DIFFERENTIAL TF + ONTOLOGY PIPELINE")

    # ── Resolve effective training parameters ────────────────────────────
    if args.minimal:
        effective_training_params, effective_finetune_params = apply_minimal_overrides(
            TRAINING_PARAMS, FINETUNE_PARAMS, MINIMAL_PARAMS
        )
        print("[MINIMAL] Effective training epochs  :", effective_training_params['epochs'])
        print("[MINIMAL] Effective finetune epochs  :", effective_finetune_params['initial_finetune_epochs'])
        print("[MINIMAL] Effective cluster ft epochs:", effective_finetune_params['cluster_finetune_epochs_max'])
    else:
        effective_training_params = TRAINING_PARAMS
        effective_finetune_params = FINETUNE_PARAMS
    
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
        
        print("\nLoading data...")
        adata = sc.read_h5ad(adata_path)
        print(f"Loaded adata: {adata.shape}")
        
        rna_data, adata_filtered = prepare_data(adata, cells_to_remove=cells_to_remove)

        # ── Minimal: subsample cells only (genes and GRN kept full) ─────
        if args.minimal:
            rna_data = subsample_cells(
                rna_data,
                n_cells=MINIMAL_PARAMS['n_cells'],
                cluster_col=DIFF_PARAMS['cluster_col']
            )
        
        skin_grn_raw = pd.read_parquet(grn_path)
        print(f"Loaded GRN: {skin_grn_raw.shape}")
        
        # Build GRN with uniform weights (peak multiplicity only)
        net = build_grn(
            skin_grn_raw=skin_grn_raw,
            save_dir=save_dir,
            scale_to_max=SCALE_TO_MAX
        )
        
        if args.train:
            model, processed_adata = train_scregulate_model(
                adata=rna_data,
                net=net,
                save_dir=save_dir,
                training_params=effective_training_params,
                finetune_params=effective_finetune_params
            )
        else:
            processed_adata_path = os.path.join(save_dir, "obj_annadata.h5ad")
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
    
    if not args.skip_plots:
        similarity_df = analyze_cluster_similarity(
            processed_adata=processed_adata,
            save_dir=save_dir,
            subset_clusters=CUSTOM_CELLTYPE_ORDER[:25]
        )
    
    results_df = differential_tf_activity(
        processed_adata=processed_adata,
        save_dir=save_dir,
        diff_params=DIFF_PARAMS
    )
    
    if not args.skip_plots:
        heatmap_result = create_differential_heatmap(
            processed_adata=processed_adata,
            save_dir=save_dir,
            heatmap_params=HEATMAP_PARAMS,
            custom_order=CUSTOM_CELLTYPE_ORDER
        )
        if heatmap_result is not None:
            heatmap_df, diff_results = heatmap_result
    
    if not args.skip_enrichment:
        up_results, down_results = run_ontology_enrichment(
            results_df=results_df,
            save_dir=save_dir,
            enrichment_params=ENRICHMENT_PARAMS
        )
        
        if not args.skip_plots and up_results is not None:
            create_enrichment_dotplot(
                up_results=up_results,
                down_results=down_results,
                save_dir=save_dir,
                custom_celltype_order=CUSTOM_CELLTYPE_ORDER
            )
    
    print(f"\n{'='*60}")
    if args.minimal:
        print("MINIMAL TEST-RUN COMPLETE  ✓")
        print("Re-run without -m / --minimal for a full scientific analysis.")
    else:
        print("PIPELINE COMPLETE")
    print(f"{'='*60}")