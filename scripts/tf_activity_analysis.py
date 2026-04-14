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

main_folder = "./"
save_dir = os.path.join("./tfactivity")
os.makedirs(save_dir, exist_ok=True)

adata_path = os.path.join(main_folder, "obj_anndata.h5ad")
grn_path = os.path.join(main_folder, "Skin_GRN_dataframe.parquet")

cells_to_remove = []

# GRN weighting parameters
MIN_WEIGHT = 0.1
MAX_WEIGHT = 1.0
SCALE_TO_MAX = 10.0

################
# COMMAND LINE ARGUMENTS
################

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='scRegulate: Weighted GRN + Differential TF + Ontology Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--load-model', '-l', type=str, default=None, metavar='PATH')
    parser.add_argument('--train', '-t', action='store_true')
    parser.add_argument('--output-dir', '-o', type=str, default=None, metavar='PATH')
    parser.add_argument('--adata-file', type=str, default="processed_adata.h5ad", metavar='FILENAME')
    parser.add_argument('--tf-file', type=str, default="tf_activities.h5ad", metavar='FILENAME')
    parser.add_argument('--grn-file', type=str, default="GRN.pkl", metavar='FILENAME')
    parser.add_argument('--model-file', type=str, default="scregulate_model.pt", metavar='FILENAME')
    parser.add_argument('--skip-enrichment', action='store_true')
    parser.add_argument('--skip-plots', action='store_true')
    parser.add_argument('--minimal', '-m', action='store_true',
        help='Minimal test-run mode: reduced cells/epochs for smoke-testing.')
    return parser.parse_args()

MINIMAL_PARAMS = {
    'n_cells': 500, 'epochs': 50, 'freeze_epochs': 20, 'batch_size': 128,
    'initial_finetune_epochs': 20, 'cluster_finetune_epochs_max': 50,
    'cluster_finetune_epochs_min': 20, 'early_stopping_patience': 50, 'log_interval': 10,
}

TRAINING_PARAMS = {
    'encode_dims': [2048, 512, 64], 'decode_dims': [512],
    'epochs': 8000, 'freeze_epochs': 4500, 'batch_size': 4096,
    'learning_rate': 2.5e-4,
    'alpha_start': 0, 'alpha_max': 1, 'alpha_scale': 0.003,
    'beta_start': 0, 'beta_max': 0.4,
    'gamma_start': 0, 'gamma_max': 3.0, 'z_dim': 40,
    'log_interval': 250, 'early_stopping_patience': 4500,
    'train_val_split_ratio': 0.85, 'min_targets': 10, 'min_TFs': 1
}

FINETUNE_PARAMS = {
    'initial_finetune_epochs': 2000,
    'cluster_finetune_epochs_max': 5000, 'cluster_finetune_epochs_min': 2000,
    'tf_mapping_lr': 2.5e-4, 'fc_output_lr': 2e-7, 'other_layers_lr': 3.5e-7
}

DIFF_PARAMS = {
    'condition_col': 'treatment', 'condition_A': 'PBS', 'condition_B': 'sCD83',
    'cluster_col': 'fine_clust', 'min_cells_per_condition': 3,
    'p_threshold': 0.05, 'min_abs_difference': 0.5
}

ENRICHMENT_PARAMS = {
    'databases': ['GO_Biological_Process_2023'],
    'qval_cutoff': 0.15, 'min_overlap': 2,
    'min_tfs_for_enrichment': 2, 'min_abs_difference_enrichment': 0.5,
}
HEATMAP_PARAMS = {
    'n_top_per_celltype': 20, 'min_score_threshold': 20,
    'selection_criterion': 'combined', 'min_cell_types_appearance': 1,
    'high_score_threshold': 50, 'final_top_n_tfs': 30,
    'cluster_rows': True, 'cluster_columns': False,
    'clustering_method': 'complete', 'clustering_metric': 'correlation',
    'colormap': 'bwr', 'show_significance_stars': True, 'dpi_resolution': 330
}

CUSTOM_CELLTYPE_ORDER = [
    "Endo.General", "Endo.Vascular", "Endo.Angiogenic", "Endo.Lymphatic",
    "FBs.Activated", "FBs.Interstitial", "FBs.Infl.Myofibroblast",
    "FBs.Pericytes", "FBs.DP-like", "FBs.Perifolicular",
    "FBs.Cycling", "FBs.NF-kB", "FBs.HF.Progenitor",
    "FBs.Dermal.Sheath", "FBs.Reticular",
    "Mast.cells", "Dendritic.cells", "M2.Macrophages",
    "MAIT.cells", "T.Naive-ANK3", "T.Naive",
    "T.cyto", "T.reg", "Plasma.cells", "B.cells",
    "Bulge", "Ishtmus", "ORS.1", "ORS.2",
    "ORS.Suprabasal", "ORS.3", "Matrix", "Eccrine.cells",
    "Melanocytes", "Schwann.cells"
]

ENRICHMENT_RANGE_COLORS = {
    (0, 4): 'lightblue', (4, 15): 'lightgreen',
    (15, 25): 'lightcoral', (25, 33): 'lightyellow', (33, 35): 'lightgray',
}

print("scRegulate: WEIGHTED GRN + DIFFERENTIAL TF + ONTOLOGY PIPELINE")


################
# MINIMAL MODE HELPERS
################

def apply_minimal_overrides(training_params, finetune_params, minimal_params):
    tp = training_params.copy()
    fp = finetune_params.copy()
    tp['epochs'] = minimal_params['epochs']
    tp['freeze_epochs'] = minimal_params['freeze_epochs']
    tp['batch_size'] = minimal_params['batch_size']
    tp['early_stopping_patience'] = minimal_params['early_stopping_patience']
    tp['log_interval'] = minimal_params['log_interval']
    fp['initial_finetune_epochs'] = minimal_params['initial_finetune_epochs']
    fp['cluster_finetune_epochs_max'] = minimal_params['cluster_finetune_epochs_max']
    fp['cluster_finetune_epochs_min'] = minimal_params['cluster_finetune_epochs_min']
    return tp, fp


def subsample_cells(adata, n_cells, cluster_col='fine_clust', seed=0):
    if adata.n_obs <= n_cells:
        print(f"[MINIMAL]   Cell count ({adata.n_obs}) already <= {n_cells}; no filter")
        return adata.copy()
    print(f"\n[MINIMAL] Subsampling cells: {adata.n_obs} -> target {n_cells}")
    rng = np.random.default_rng(seed)
    if cluster_col in adata.obs.columns:
        clusters = adata.obs[cluster_col].values
        unique_cl = np.unique(clusters)
        n_per = max(3, n_cells // len(unique_cl))
        keep_idx = []
        for cl in unique_cl:
            idx = np.where(clusters == cl)[0]
            chosen = rng.choice(idx, size=min(n_per, len(idx)), replace=False)
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
    import scipy.sparse as sp
    n_cells, n_genes = X.shape
    cell_idx = np.random.default_rng(42).choice(n_cells, size=min(n_sample_cells, n_cells), replace=False)
    gene_idx = np.random.default_rng(42).choice(n_genes, size=min(n_sample_genes, n_genes), replace=False)
    X_sample = X[np.ix_(cell_idx, gene_idx)]
    if sp.issparse(X_sample): X_sample = X_sample.toarray()
    X_sample = X_sample.astype(float)
    has_negative = bool((X_sample < 0).any())
    global_max = float(X_sample.max())
    global_min = float(X_sample.min())
    global_mean = float(X_sample.mean())
    global_std = float(X_sample.std())
    nonzero = X_sample[X_sample != 0]
    frac_integer = float((np.abs(nonzero - np.round(nonzero)) < 1e-4).mean()) if len(nonzero) else 1.0
    looks_integer = frac_integer > 0.99
    cell_sums = X_sample.sum(axis=1)
    cell_sum_cv = float(cell_sums.std() / (cell_sums.mean() + 1e-9))
    stats = {'has_negative': has_negative, 'global_min': round(global_min, 4),
        'global_max': round(global_max, 4), 'global_mean': round(global_mean, 4),
        'global_std': round(global_std, 4), 'frac_integer': round(frac_integer, 4),
        'cell_sum_cv': round(cell_sum_cv, 4)}
    if has_negative:
        return dict(status='sct_residuals', needs_norm=False, needs_log=False,
            reason="Negative values — SCTransform Pearson residuals.", stats=stats)
    if not looks_integer and global_max < 30:
        return dict(status='log_normalized', needs_norm=False, needs_log=False,
            reason=f"Non-integer, max={global_max:.2f} — already log-normalized.", stats=stats)
    if looks_integer and global_max > 30:
        return dict(status='raw_counts', needs_norm=True, needs_log=True,
            reason=f"Integer-valued, max={global_max:.0f} — raw counts.", stats=stats)
    return dict(status='ambiguous', needs_norm=True, needs_log=True,
        reason="Ambiguous matrix. Applying normalize_total + log1p.", stats=stats)


def print_normalization_report(norm_info, source_label):
    s = norm_info['stats']
    label = {'raw_counts': 'RAW COUNTS', 'log_normalized': 'ALREADY LOG-NORMALIZED',
        'sct_residuals': 'SCTransform RESIDUALS', 'ambiguous': 'AMBIGUOUS'
    }.get(norm_info['status'], norm_info['status'].upper())
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
    print("DATA PREPARATION")
    if 'treatment' not in adata.obs.columns:
        print("\nCreating treatment column from cell index...")
        adata.obs['treatment'] = adata.obs.index.str.split('_').str[0]
        adata.obs['treatment'] = adata.obs['treatment'].map({'pbs': 'PBS', 'sCD83': 'sCD83'})
        print(f"  Treatment distribution: {adata.obs['treatment'].value_counts().to_dict()}")
    if cells_to_remove and len(cells_to_remove) > 0:
        print(f"\nRemoving cell types: {cells_to_remove}")
        adata = adata[~adata.obs['FineClust'].isin(cells_to_remove)].copy()
    print("\nExtracting and normalizing data...")
    if adata.raw is not None:
        rna_data = sc.AnnData(X=adata.raw.X.copy(), obs=adata.obs.copy(), var=adata.raw.var.copy())
        source_label = "adata.raw"
        print("  Source: adata.raw")
    else:
        rna_data = adata.copy()
        source_label = "adata.X"
        print("  Source: adata.X (no raw layer found)")
    norm_info = detect_normalization_status(rna_data.X)
    print_normalization_report(norm_info, source_label)
    if norm_info['needs_norm'] and norm_info['needs_log']:
        sc.pp.normalize_total(rna_data); sc.pp.log1p(rna_data)
        print("\n  Applied: normalize_total + log1p")
    elif norm_info['needs_log']:
        sc.pp.log1p(rna_data); print("\n  Applied: log1p only")
    else:
        print("\n  Skipped normalization — already transformed")
    if 'X_umap' in adata.obsm:
        rna_data.obsm['X_umap'] = adata.obsm['X_umap']
    print(f"\nPrepared data: {rna_data.shape}")
    return rna_data, adata


################
# BUILD EXPRESSION-WEIGHTED GRN
################

def build_weighted_grn(adata, skin_grn_raw, save_dir,
                       min_weight=0.1, max_weight=1.0, scale_to_max=10.0,
                       cluster_col='fine_clust'):
    """
    Build expression-weighted GRN from Greenleaf ATAC data.

    Weights are based on max-cluster-mean expression: for each TF, we take
    the maximum of per-cluster mean expression values. This captures TFs
    that are highly active in at least one cell type (e.g. SOX10 in
    melanocytes) rather than penalizing cell-type-specific TFs with low
    global mean expression. scRegulate's cluster-specific fine-tuning
    then learns which TFs matter in which cell types.
    """
    print("BUILDING EXPRESSION-WEIGHTED GRN (max-cluster-mean)")
    tf_columns = [col for col in skin_grn_raw.columns if col not in ['peak_id', 'gene_short_name']]
    print(f"\nNumber of TFs in GRN: {len(tf_columns)}")

    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = adata.X
    gene_names = adata.var_names.tolist()

    # Get cluster labels for per-cluster mean calculation
    if cluster_col in adata.obs.columns:
        clusters = adata.obs[cluster_col].values
        unique_clusters = np.unique(clusters)
        print(f"  Using {len(unique_clusters)} clusters from '{cluster_col}' for per-cluster expression")
    else:
        print(f"  WARNING: '{cluster_col}' not found in obs — falling back to global mean")
        clusters = None

    print("\nCalculating TF expression weights (max-cluster-mean)...")
    tf_expression = {}
    for tf in tf_columns:
        if tf in gene_names:
            tf_idx = gene_names.index(tf)
            if clusters is not None:
                # Max of per-cluster means: captures TFs important in ANY cell type
                cluster_means = [X[clusters == c, tf_idx].mean() for c in unique_clusters]
                tf_expression[tf] = max(cluster_means)
            else:
                tf_expression[tf] = X[:, tf_idx].mean()
        else:
            tf_expression[tf] = 0.0

    max_expr = max(tf_expression.values()) if tf_expression.values() else 1.0
    print(f"Max TF expression (max-cluster-mean): {max_expr:.4f}")

    tf_weights = {}
    for tf, expr in tf_expression.items():
        if expr <= 0:
            weight = min_weight
        else:
            weight = (expr / max_expr) * (max_weight - min_weight) + min_weight
        tf_weights[tf] = weight

    weight_values = list(tf_weights.values())
    print(f"\nWeight distribution:")
    print(f"  Min:    {min(weight_values):.4f}")
    print(f"  25%:    {np.percentile(weight_values, 25):.4f}")
    print(f"  Median: {np.median(weight_values):.4f}")
    print(f"  75%:    {np.percentile(weight_values, 75):.4f}")
    print(f"  Max:    {max(weight_values):.4f}")

    print("\nExample TF weights:")
    for tf in ['MITF', 'SOX18', 'PAX5', 'TCF7', 'SOX9', 'MYF5', 'SOX10', 'TP63', 'KLF4', 'FOXP3', 'STAT3']:
        if tf in tf_weights:
            print(f"  {tf:10s}: expr={tf_expression[tf]:7.4f}, weight={tf_weights[tf]:.4f}")

    print("\nConverting to long format with expression weights...")
    tf_target_pairs = []
    for idx, row in skin_grn_raw.iterrows():
        target_gene = row['gene_short_name']
        peak_id = row['peak_id']
        for tf in tf_columns:
            if row[tf] == 1:
                tf_target_pairs.append({
                    'source': tf, 'target': target_gene,
                    'peak_id': peak_id, 'weight': tf_weights[tf]
                })
        if (idx + 1) % 10000 == 0:
            print(f"  Processed {idx + 1}/{len(skin_grn_raw)} peaks...")

    grn_long = pd.DataFrame(tf_target_pairs)
    print(f"\nTotal TF-peak-gene edges: {len(grn_long)}")

    print("\nAggregating by TF-target pairs...")
    grn_aggregated = grn_long.groupby(['source', 'target'], as_index=False).agg({
        'weight': 'sum', 'peak_id': lambda x: ';'.join(x)
    })
    print(f"Unique TF-target pairs: {len(grn_aggregated)}")

    max_weight_99 = grn_aggregated['weight'].quantile(0.99)
    grn_aggregated['weight'] = (grn_aggregated['weight'] / max_weight_99) * scale_to_max
    grn_aggregated['weight'] = grn_aggregated['weight'].clip(upper=scale_to_max)

    print(f"\nFinal weight distribution (scaled to [0, {scale_to_max}]):")
    print(f"  Mean:   {grn_aggregated['weight'].mean():.4f}")
    print(f"  Median: {grn_aggregated['weight'].median():.4f}")
    print(f"  Max:    {grn_aggregated['weight'].max():.4f}")

    net_weighted = grn_aggregated[['source', 'target', 'weight']].copy()
    net_weighted['PMID'] = 'CellOracle_' + grn_aggregated['peak_id'].astype(str)
    net_weighted = net_weighted[['source', 'target', 'weight', 'PMID']]

    min_targets = TRAINING_PARAMS.get('min_targets', 10)
    tf_target_counts = net_weighted.groupby('source').size()
    valid_tfs = tf_target_counts[tf_target_counts >= min_targets].index
    net_weighted = net_weighted[net_weighted['source'].isin(valid_tfs)]

    print(f"\nGRN Statistics (after filtering TFs with >={min_targets} targets):")
    print(f"  Weighted edges: {len(net_weighted):,}")
    print(f"  TFs: {net_weighted['source'].nunique():,}")
    print(f"  Target genes: {net_weighted['target'].nunique():,}")

    output_path = os.path.join(save_dir, "net_weighted_grn.parquet")
    net_weighted.to_parquet(output_path, index=False)
    print(f"\nSaved weighted GRN to: {output_path}")

    return net_weighted[['source', 'target', 'weight']].copy(), tf_weights


################
# TRAIN scRegulate MODEL
################

def train_scregulate_model(adata, net, save_dir, training_params, finetune_params):
    """
    Train scRegulate model, then two-stage fine-tuning per official tutorial:
    Stage 1: Global (no cluster_key)  Stage 2: Per-cluster (with cluster_key)
    """
    print("TRAINING scRegulate MODEL")
    try:
        import scregulate as reg
    except ImportError:
        print("ERROR: scregulate not installed."); return None, None

    if torch.cuda.is_available():
        device = torch.device("cuda")
        gpu_name = torch.cuda.get_device_name(0)
        print(f"\nGPU detected: {gpu_name}")
        print(f"  VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    else:
        device = torch.device("cpu")
        print("\nNo GPU detected — training on CPU")
    print(f"Device: {device}")

    rna_data = adata
    print(f"\nRNA data shape: {rna_data.shape}")
    print(f"GRN edges: {len(net)}")
    print(f"\nTraining scRegulate model...")
    print(f"Encoder: {training_params['encode_dims']}, z_dim={training_params['z_dim']}")
    print(f"Decoder: {training_params['decode_dims']}")
    print(f"Epochs: {training_params['epochs']}, Batch size: {training_params['batch_size']}")
    print(f"Learning rate: {training_params['learning_rate']}, Patience: {training_params['early_stopping_patience']}")

    model, processed_adata, GRN = reg.train_model(
        rna_data=rna_data, net=net,
        encode_dims=training_params['encode_dims'], z_dim=training_params['z_dim'],
        decode_dims=training_params['decode_dims'], epochs=training_params['epochs'],
        freeze_epochs=training_params['freeze_epochs'], batch_size=training_params['batch_size'],
        learning_rate=training_params['learning_rate'],
        alpha_start=training_params['alpha_start'], alpha_max=training_params['alpha_max'],
        alpha_scale=training_params['alpha_scale'],
        beta_start=training_params['beta_start'], beta_max=training_params['beta_max'],
        gamma_start=training_params['gamma_start'], gamma_max=training_params['gamma_max'],
        log_interval=training_params['log_interval'],
        early_stopping_patience=training_params['early_stopping_patience'],
        train_val_split_ratio=training_params['train_val_split_ratio'],
        min_targets=training_params['min_targets'], min_TFs=training_params['min_TFs'],
        device=device, verbose=True)

    print("\nInitial training complete!")
    for col in adata.obs.columns:
        if col not in processed_adata.obs.columns:
            processed_adata.obs[col] = adata.obs[col].values

    # Stage 1: Global fine-tuning (no cluster_key)
    print(f"\nStage 1: Global fine-tuning ({finetune_params['initial_finetune_epochs']} epochs)...")
    processed_adata, fine_tuned_tf_activities, fine_model, GRN_global = reg.fine_tuning.fine_tune_clusters(
        model=model, processed_adata=processed_adata,
        epochs=finetune_params['initial_finetune_epochs'],
        tf_mapping_lr=finetune_params['tf_mapping_lr'],
        fc_output_lr=finetune_params['fc_output_lr'],
        default_lr=finetune_params['other_layers_lr'],
        device=device, verbose=True)

    # Stage 2: Cluster-specific fine-tuning
    cluster_col = DIFF_PARAMS['cluster_col']
    print(f"\nStage 2: Cluster-specific fine-tuning (cluster_key='{cluster_col}')...")
    processed_adata, fine_tuned_tf_activities, final_model, GRN_final = reg.fine_tuning.fine_tune_clusters(
        processed_adata=processed_adata, model=fine_model,
        cluster_key=cluster_col,
        epochs=finetune_params['cluster_finetune_epochs_max'],
        min_epochs=finetune_params['cluster_finetune_epochs_min'],
        tf_mapping_lr=finetune_params['tf_mapping_lr'],
        fc_output_lr=finetune_params['fc_output_lr'],
        default_lr=finetune_params['other_layers_lr'],
        device=device, verbose=True)

    print("\nFine-tuning complete!")

    # Save outputs
    torch.save(final_model.state_dict(), os.path.join(save_dir, "scregulate_model.pt"))
    if 'W_posteriors_per_cluster' in processed_adata.uns:
        processed_adata.uns['W_posteriors_per_cluster'] = {
            str(k): v for k, v in processed_adata.uns['W_posteriors_per_cluster'].items()}
    if 'TF_finetuned' in processed_adata.obsm and 'TF_activity' not in processed_adata.obsm:
        processed_adata.obsm['TF_activity'] = processed_adata.obsm['TF_finetuned']
    if 'TF_names' not in processed_adata.uns and fine_tuned_tf_activities is not None:
        processed_adata.uns['TF_names'] = fine_tuned_tf_activities.var_names.tolist()

    processed_adata.write_h5ad(os.path.join(save_dir, "processed_adata.h5ad"))
    fine_tuned_tf_activities.write_h5ad(os.path.join(save_dir, "tf_activities.h5ad"))
    with open(os.path.join(save_dir, "GRN.pkl"), "wb") as f:
        pickle.dump(GRN_final, f)
    print("All outputs saved.")

    return final_model, processed_adata


################
# DIFFERENTIAL TF ACTIVITY ANALYSIS
################

def differential_tf_activity(processed_adata, save_dir, diff_params):
    print("DIFFERENTIAL TF ACTIVITY ANALYSIS")
    from scipy import stats

    condition_col, condition_A, condition_B = diff_params['condition_col'], diff_params['condition_A'], diff_params['condition_B']
    cluster_col = diff_params['cluster_col']
    p_threshold, min_cells = diff_params['p_threshold'], diff_params['min_cells_per_condition']
    min_abs_diff = diff_params['min_abs_difference']
    PSEUDOCOUNT = 0.0001

    if 'TF_activity' in processed_adata.obsm:
        tf_activity = processed_adata.obsm['TF_activity']
    elif 'TF_finetuned' in processed_adata.obsm:
        tf_activity = processed_adata.obsm['TF_finetuned']
    else:
        raise KeyError("No TF activity found in obsm")

    tf_names = list(processed_adata.uns['TF_names'])
    print(f"\nTF activity matrix: {tf_activity.shape}")
    print(f"Conditions: {condition_A} vs {condition_B}")
    clusters = processed_adata.obs[cluster_col].unique()
    print(f"Number of clusters: {len(clusters)}")

    results = []
    for cluster in clusters:
        mask = processed_adata.obs[cluster_col] == cluster
        conds = processed_adata.obs.loc[mask, condition_col]
        mA, mB = conds == condition_A, conds == condition_B
        nA, nB = mA.sum(), mB.sum()
        if nA < min_cells or nB < min_cells:
            print(f"  Skipping {cluster}: (A={nA}, B={nB})"); continue
        print(f"  Processing {cluster}: {nA} {condition_A}, {nB} {condition_B}")
        data = tf_activity[mask]
        for i, tf in enumerate(tf_names):
            tA, tB = data[mA, i], data[mB, i]
            mA_val, mB_val = tA.mean(), tB.mean()
            diff = mB_val - mA_val
            l2fc = np.log2((mB_val + PSEUDOCOUNT) / (mA_val + PSEUDOCOUNT))
            try: _, pval = stats.ranksums(tB, tA)
            except: pval = 1.0
            ps = np.sqrt(((nA-1)*tA.std()**2 + (nB-1)*tB.std()**2) / (nA+nB-2))
            cd = diff / ps if ps > 0 else 0.0
            results.append({'cell_type': cluster, 'TF': tf, 'mean_A': mA_val, 'mean_B': mB_val,
                'difference': diff, 'log2fc': l2fc, 'pvalue': pval, 'cohens_d': cd, 'n_A': nA, 'n_B': nB})

    results_df = pd.DataFrame(results)
    results_df['p_adjusted'] = 1.0
    for cluster in results_df['cell_type'].unique():
        m = results_df['cell_type'] == cluster
        _, padj, _, _ = multipletests(results_df.loc[m, 'pvalue'].values, method='fdr_bh')
        results_df.loc[m, 'p_adjusted'] = padj
    results_df['significant'] = (results_df['p_adjusted'] < p_threshold) & (np.abs(results_df['difference']) >= min_abs_diff)

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
    gene_name = str(gene_name)
    for p in [r'^ENSG\d+', r'^ENSMUSG\d+', r'^AC\d+\.\d+', r'^AL\d+\.\d+',
              r'^LINC\d+', r'^IRC\d+', r'^LOC\d+', r'^FP\d+',
              r'^RP\d+-\d+', r'^CTD-\d+', r'^KB-\d+', r'^\d+']:
        if re.match(p, gene_name): return False
    if '.' in gene_name: return False
    if len(gene_name) < 2 or not any(c.isalpha() for c in gene_name): return False
    return True


def get_top_tfs_per_celltype_filtered(differential_results, tf_names, n_top=10,
                                       criterion='combined', min_score_threshold=None):
    top_tfs_dict = {}
    for cell_type, results in differential_results.items():
        if criterion == 'combined':
            score = np.abs(results['fold_change']) * -np.log10(results['p_adjusted'] + 1e-300)
        elif criterion == 'fold_change':
            score = np.abs(results['fold_change'])
        elif criterion == 'significance':
            score = -np.log10(results['p_adjusted'] + 1e-300)
        else:
            score = np.abs(results['difference'])
        top_indices = np.argsort(score)[-n_top:][::-1]
        if min_score_threshold is not None:
            top_indices = [i for i in top_indices if score[i] >= min_score_threshold]
        valid_tfs = [tf_names[i] for i in top_indices if is_valid_gene_name(tf_names[i])]
        if len(valid_tfs) < n_top:
            ext = np.argsort(score)[-(n_top*3):][::-1]
            if min_score_threshold is not None:
                ext = [i for i in ext if score[i] >= min_score_threshold]
            valid_tfs = [tf_names[i] for i in ext if is_valid_gene_name(tf_names[i])][:n_top]
        top_tfs_dict[cell_type] = valid_tfs
    return top_tfs_dict


################
# CLUSTER SIMILARITY ANALYSIS
################

def analyze_cluster_similarity(processed_adata, save_dir, subset_clusters=None):
    print("CLUSTER SIMILARITY ANALYSIS")
    W_post = processed_adata.uns.get("W_posteriors_per_cluster", None)
    if W_post is None:
        print("No W_posteriors_per_cluster found"); return None
    ct_cols = list(W_post.keys())
    print("\nBuilding average W matrices...")
    avg_W = {ct: minmax_scale(np.abs(W_post[ct]).ravel()).reshape(W_post[ct].shape).mean(axis=0) for ct in ct_cols}
    combined = pd.DataFrame(avg_W).T
    sim = np.clip(cosine_similarity(combined), 0, 1)**8
    sim_df = pd.DataFrame(sim, index=ct_cols, columns=ct_cols)

    print("\nCreating similarity heatmap...")
    sp = sns.clustermap(sim_df, figsize=(12, 12), annot=False, cmap="RdBu_r", center=0, cbar=False,
                        dendrogram_ratio=(0.05, 0.05))
    sp.cax.set_visible(False)
    sp.ax_heatmap.set_xticklabels(sp.ax_heatmap.get_xticklabels(), fontsize=11, fontweight='bold', rotation=90, ha='right')
    sp.ax_heatmap.set_yticklabels(sp.ax_heatmap.get_yticklabels(), fontsize=11, fontweight='bold')
    sp.savefig(os.path.join(save_dir, 'tf_grn_similarity_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()

    if subset_clusters is not None:
        sub_W = {ct: avg_W[ct] for ct in subset_clusters if ct in avg_W}
        if len(sub_W) >= 2:
            sub_df = pd.DataFrame(sub_W).T
            sub_df = sub_df.apply(lambda x: minmax_scale(x), axis=0)
            sub_sim = np.clip(cosine_similarity(sub_df), 0, 1)**8
            sub_sim_df = pd.DataFrame(sub_sim, index=list(sub_W.keys()), columns=list(sub_W.keys()))
            sp2 = sns.clustermap(sub_sim_df, figsize=(12, 10), annot=True, fmt=".1f",
                annot_kws={"size": 8}, cmap="RdBu_r", center=0, cbar=False,
                dendrogram_ratio=(0.08, 0.08))
            sp2.cax.set_visible(False)
            sp2.ax_heatmap.set_xticklabels(sp2.ax_heatmap.get_xticklabels(), fontsize=11, fontweight='bold', rotation=90, ha='right')
            sp2.ax_heatmap.set_yticklabels(sp2.ax_heatmap.get_yticklabels(), fontsize=11, fontweight='bold')
            sp2.savefig(os.path.join(save_dir, 'tf_grn_similarity_heatmap_zoom.png'), dpi=300, bbox_inches='tight')
            plt.close()
    return sim_df


################
# DIFFERENTIAL TF HEATMAP
################

def create_differential_heatmap(processed_adata, save_dir, heatmap_params, custom_order=None):
    print("CREATING DIFFERENTIAL TF HEATMAP")
    if 'TF_finetuned' in processed_adata.obsm: W_matrix = processed_adata.obsm['TF_finetuned']
    elif 'TF_activity' in processed_adata.obsm: W_matrix = processed_adata.obsm['TF_activity']
    else: print("No TF activity found"); return None, None
    if hasattr(W_matrix, 'toarray'): W_matrix = W_matrix.toarray()

    tf_names = list(processed_adata.uns['TF_names'])
    conditions = processed_adata.obs['treatment'].values
    cell_types = processed_adata.obs['fine_clust'].values
    ct_cols = np.unique(cell_types)

    print("\nPerforming differential analysis...")
    diff_res = {}
    for ct in ct_cols:
        cm = cell_types == ct
        pm, sm = cm & (conditions == 'PBS'), cm & (conditions == 'sCD83')
        if pm.sum() >= 3 and sm.sum() >= 3:
            Wp, Ws = W_matrix[pm, :], W_matrix[sm, :]
            diff = Ws.mean(0) - Wp.mean(0)
            fc = np.log2((Ws.mean(0) + 0.0001) / (Wp.mean(0) + 0.0001))
            pv = np.array([mannwhitneyu(Wp[:, j], Ws[:, j], alternative='two-sided')[1]
                if Wp[:, j].std() + Ws[:, j].std() > 0 else 1.0 for j in range(W_matrix.shape[1])])
            _, padj, _, _ = multipletests(pv, alpha=0.05, method='fdr_bh')
            diff_res[ct] = {'difference': diff, 'fold_change': fc, 'p_adjusted': padj, 'significant': padj < 0.05}

    if not diff_res: print("No cell types with data"); return None, None

    print("\nSelecting top TFs...")
    top_tfs = get_top_tfs_per_celltype_filtered(diff_res, tf_names,
        n_top=heatmap_params['n_top_per_celltype'], criterion=heatmap_params['selection_criterion'],
        min_score_threshold=heatmap_params['min_score_threshold'])

    tf_freq, tf_max = Counter(), {}
    for ct, tfs in top_tfs.items():
        for tf in tfs:
            ti = tf_names.index(tf)
            tf_freq[tf] += 1
            s = np.abs(diff_res[ct]['fold_change'][ti]) * -np.log10(diff_res[ct]['p_adjusted'][ti] + 1e-300)
            if tf not in tf_max or s > tf_max[tf]: tf_max[tf] = s

    filtered = [tf for tf in tf_freq if tf_freq[tf] >= heatmap_params['min_cell_types_appearance']
                or tf_max[tf] >= heatmap_params['high_score_threshold']]
    if not filtered: print("No TFs passed filtering"); return None, None

    analyzed = [ct for ct in ct_cols if ct in diff_res]
    var_scores = [np.var([diff_res[ct]['difference'][tf_names.index(tf)] for ct in analyzed]) for tf in filtered]
    if len(filtered) > heatmap_params['final_top_n_tfs']:
        idx = np.argsort(var_scores)[-heatmap_params['final_top_n_tfs']:]
        all_top = [filtered[i] for i in idx]
    else:
        all_top = filtered
    print(f"Selected {len(all_top)} TFs for heatmap")

    data = [[diff_res[ct]['difference'][tf_names.index(tf)] for tf in all_top] for ct in analyzed]
    hm = pd.DataFrame(data, index=analyzed, columns=all_top).T
    if custom_order:
        oc = [ct for ct in custom_order if ct in hm.columns]; hm = hm[oc]
    else: oc = list(hm.columns)

    if heatmap_params['cluster_rows'] and len(all_top) >= 2:
        rl = linkage(hm, method='complete', metric='correlation')
        ro = dendrogram(rl, no_plot=True)['leaves']
        ot = [all_top[i] for i in ro]; hm = hm.loc[ot]
    else: ot = all_top

    sig = np.zeros_like(hm.values, dtype=str)
    if heatmap_params['show_significance_stars']:
        for i, tf in enumerate(ot):
            for j, ct in enumerate(oc):
                pa = diff_res[ct]['p_adjusted'][tf_names.index(tf)]
                if pa < 0.001: sig[i, j] = '***'
                elif pa < 0.01: sig[i, j] = '**'
                elif pa < 0.05: sig[i, j] = '*'

    fh = max(8, len(ot) * 0.3)
    fig, ax = plt.subplots(figsize=(11, fh))
    
    from matplotlib.colors import LinearSegmentedColormap
    
    custom_cmap = LinearSegmentedColormap.from_list('pbs_scd83', ['#F8766C', '#FFFFFF', '#00BEF9'])
    sns.heatmap(hm, cmap=custom_cmap, center=0, robust=True,
        cbar_kws={'label': 'TF Activity Change\n(sCD83 - PBS)'},
        yticklabels=True, xticklabels=True, linewidths=0.5, linecolor='lightgray', ax=ax)
    if heatmap_params['show_significance_stars']:
        for i in range(len(ot)):
            for j in range(len(oc)):
                if sig[i, j]:
                    ax.text(j+0.5, i+0.5, sig[i, j], ha='center', va='center', color='black', fontsize=8, fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=12, fontweight='bold', rotation=90, ha='center', va='top')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=11, fontweight='bold')
    ax.set_title('Differential TF Activity (sCD83 vs PBS)', fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'differential_tf_activity_heatmap.png'), dpi=heatmap_params['dpi_resolution'], bbox_inches='tight')
    plt.close()
    print(f"\nSaved heatmap")
    return hm, diff_res


################
# MODEL LOAD
################

def load_model_outputs(load_dir, adata_file="processed_adata.h5ad",
                       tf_file="tf_activities.h5ad", grn_file="GRN.pkl",
                       model_file="scregulate_model.pt"):
    print("LOADING MODEL OUTPUTS")
    import scregulate as reg

    processed_adata = ad.read_h5ad(os.path.join(load_dir, adata_file))
    print(f"Loaded processed_adata: {processed_adata.shape}")

    tf_path = os.path.join(load_dir, tf_file)
    tf_activities = ad.read_h5ad(tf_path) if os.path.exists(tf_path) else None

    with open(os.path.join(load_dir, grn_file), "rb") as f:
        GRN = pickle.load(f)

    # Ensure keys exist
    if 'TF_activity' not in processed_adata.obsm:
        if 'TF_finetuned' in processed_adata.obsm:
            processed_adata.obsm['TF_activity'] = processed_adata.obsm['TF_finetuned']
        elif tf_activities is not None:
            processed_adata.obsm['TF_activity'] = tf_activities.X.toarray() if hasattr(tf_activities.X, 'toarray') else tf_activities.X.copy()
    if 'TF_names' not in processed_adata.uns and tf_activities is not None:
        processed_adata.uns['TF_names'] = tf_activities.var_names.tolist()
    if 'TF_names' in processed_adata.uns:
        processed_adata.uns['TF_names'] = list(processed_adata.uns['TF_names'])

    # Determine dimensions
    if isinstance(GRN, (np.ndarray, pd.DataFrame)):
        input_dim, tf_dim = GRN.shape
    else:
        input_dim = processed_adata.n_vars
        tf_dim = len(processed_adata.uns.get('TF_names', []))

    model = reg.scRNA_VAE(input_dim=input_dim, encode_dims=TRAINING_PARAMS['encode_dims'],
        decode_dims=TRAINING_PARAMS['decode_dims'], z_dim=TRAINING_PARAMS['z_dim'], tf_dim=tf_dim)
    mp = os.path.join(load_dir, model_file)
    if os.path.exists(mp):
        model.load_state_dict(torch.load(mp, map_location='cpu'))
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = model.to(device)
        print(f"Loaded model to device: {device}")
    else:
        model = None

    processed_adata.modality = {}
    processed_adata.modality['RNA'] = processed_adata.copy()
    if tf_activities is not None:
        processed_adata.modality['TF'] = tf_activities.copy()

    return model, processed_adata, tf_activities, GRN


################
# ONTOLOGY ENRICHMENT
################

def run_ontology_enrichment(results_df, save_dir, enrichment_params):
    """
    Enrichment on significant TFs directly, split by cell type and direction.
    
    Uses the TFs that passed differential activity thresholds
    (FDR < 0.05, |diff| >= 0.5) as gene lists for Enrichr.
    """
    print("ONTOLOGY ENRICHMENT ANALYSIS (significant TFs)")
    try:
        import gseapy as gp
    except ImportError:
        print("gseapy not installed"); return None, None

    p_thr = DIFF_PARAMS['p_threshold']
    min_diff = enrichment_params['min_abs_difference_enrichment']
    min_tfs = enrichment_params['min_tfs_for_enrichment']
    dbs = enrichment_params['databases']
    qcut = enrichment_params['qval_cutoff']
    min_ov = enrichment_params['min_overlap']

    # ---- collect significant TFs per cell type per direction ----
    up_tfs, dn_tfs = {}, {}
    for ct in results_df['cell_type'].unique():
        d = results_df[results_df['cell_type'] == ct]
        sig = d[(d['p_adjusted'] < p_thr) & (d['difference'].abs() >= min_diff)]

        su = sig[sig['difference'] > 0].sort_values('difference', ascending=False)
        sd = sig[sig['difference'] < 0].sort_values('difference', ascending=True)

        if len(su) >= min_tfs:
            up_tfs[ct] = su['TF'].tolist()
        if len(sd) >= min_tfs:
            dn_tfs[ct] = sd['TF'].tolist()

    print(f"\nCell types with >={min_tfs} UP TFs: {len(up_tfs)}")
    for ct, tfs in up_tfs.items():
        print(f"  {ct}: {len(tfs)} TFs")
    print(f"\nCell types with >={min_tfs} DOWN TFs: {len(dn_tfs)}")
    for ct, tfs in dn_tfs.items():
        print(f"  {ct}: {len(tfs)} TFs")

    # ---- run enrichment ----
    def _run(tf_dict, direction):
        res = {}
        for i, (ct, tl) in enumerate(tf_dict.items(), 1):
            print(f"\n[{i}/{len(tf_dict)}] {ct}: {len(tl)} {direction} TFs", end=' ')
            try:
                enr = gp.enrichr(
                    gene_list=tl,
                    gene_sets=dbs,
                    organism='human',
                    outdir=None,
                    no_plot=True,
                    cutoff=1.0)  # don't pre-filter, we handle it

                if enr.results is not None and len(enr.results) > 0:
                    r = enr.results.copy()
                    r['cell_type'] = ct
                    r['direction'] = direction
                    r['n_input_tfs'] = len(tl)
                    r['input_TFs'] = ';'.join(tl)

                    # Handle column names across gseapy versions
                    pval_col = next((c for c in r.columns
                                    if c in ['Adjusted P-value', 'FDR q-val',
                                             'Adjusted_P-value', 'adjusted_p_value']), None)
                    overlap_col = next((c for c in r.columns
                                       if c in ['Overlap', 'overlap']), None)

                    if pval_col is None:
                        pval_candidates = [c for c in r.columns
                                          if 'p' in c.lower() or 'fdr' in c.lower()]
                        pval_col = pval_candidates[0] if pval_candidates else None

                    if pval_col is None:
                        print(f"-> Cannot find p-value column: {list(r.columns)}")
                        continue

                    # Standardize
                    r['Adjusted P-value'] = r[pval_col].astype(float)

                    # Filter by significance
                    r = r[r['Adjusted P-value'] < qcut]

                    # Filter by overlap if available
                    if overlap_col and len(r) > 0:
                        overlap_str = r[overlap_col].astype(str)
                        if overlap_str.str.contains('/').any():
                            overlap_count = overlap_str.str.split('/').str[0].astype(int)
                        else:
                            overlap_count = overlap_str.astype(int)
                        r = r[overlap_count >= min_ov]

                    if len(r) > 0:
                        res[ct] = r
                        print(f"-> {len(r)} terms")
                    else:
                        print("-> No enrichment")
                else:
                    print("-> No results")
            except Exception as e:
                print(f"-> Error: {str(e)[:120]}")
        return res

    print("\n" + "="*50)
    print("RUNNING ENRICHMENT FOR UP-REGULATED TFs")
    print("="*50)
    up_res = _run(up_tfs, 'UP')

    print("\n" + "="*50)
    print("RUNNING ENRICHMENT FOR DOWN-REGULATED TFs")
    print("="*50)
    dn_res = _run(dn_tfs, 'DOWN')

    # ---- save ----
    all_r = list(up_res.values()) + list(dn_res.values())
    if all_r:
        combined = pd.concat(all_r, ignore_index=True)
        combined.to_csv(os.path.join(save_dir, "enrichment_up_down_all.csv"), index=False)
        print(f"\nSaved {len(combined)} enriched terms to enrichment_up_down_all.csv")
        return up_res, dn_res

    print("\nNo enrichment results")
    return {}, {}

######################################
# MODEL FITTING ANALYSIS
######################################


def post_hoc_model_quality(model, processed_adata, tf_activities, GRN, save_dir,
                           cluster_col='fine_clust', condition_col='treatment'):
    """
    Post-hoc quality assessment of a trained scRegulate model.

    Assesses quality from saved outputs: TF activities, W_posteriors.
    Does NOT rely on model forward pass or marker-based validation
    (TF activity ≠ TF expression).

    Generates:
      1. TF activity vs expression correlation (sanity check)
      2. TF activity space structure (silhouette, UMAP)
      3. W_posteriors analysis (learned TF importance, cleaned)
      4. Effect size analysis (distribution of differential TF effects)
      5. Per-cluster differential TF summary
    """
    from scipy.stats import spearmanr
    from sklearn.metrics import silhouette_score, silhouette_samples

    qc_dir = os.path.join(save_dir, "model_quality")
    os.makedirs(qc_dir, exist_ok=True)

    clusters = processed_adata.obs[cluster_col].values
    unique_clusters = np.unique(clusters)
    conditions = processed_adata.obs[condition_col].values if condition_col in processed_adata.obs.columns else None

    tf_names = list(processed_adata.uns.get('TF_names', []))

    if 'TF_activity' in processed_adata.obsm:
        tf_matrix = processed_adata.obsm['TF_activity']
    elif 'TF_finetuned' in processed_adata.obsm:
        tf_matrix = processed_adata.obsm['TF_finetuned']
    else:
        print("ERROR: No TF activity matrix found"); return None
    if hasattr(tf_matrix, 'toarray'):
        tf_matrix = tf_matrix.toarray()

    print("MODEL QUALITY ASSESSMENT")
    print(f"  Cells: {processed_adata.n_obs}")
    print(f"  Genes: {processed_adata.n_vars}")
    print(f"  Clusters: {len(unique_clusters)}")
    print(f"  TFs: {len(tf_names)}")
    print(f"  TF activity shape: {tf_matrix.shape}")

    summary = {
        'n_cells': processed_adata.n_obs,
        'n_genes': processed_adata.n_vars,
        'n_clusters': len(unique_clusters),
        'n_tfs': len(tf_names),
    }

    # ================================================================
    # 1. TF ACTIVITY vs EXPRESSION CORRELATION
    # ================================================================
    print("\n1. TF ACTIVITY vs EXPRESSION CORRELATION")
    print("   (Sanity check: do learned activities relate to expression at all?)")

    if hasattr(processed_adata.X, 'toarray'):
        X_expr = processed_adata.X.toarray()
    else:
        X_expr = np.array(processed_adata.X)

    gene_names = processed_adata.var_names.tolist()

    tf_expr_corrs = []
    for i, tf in enumerate(tf_names):
        if tf in gene_names:
            gene_idx = gene_names.index(tf)
            tf_act = tf_matrix[:, i]
            tf_expr = X_expr[:, gene_idx]

            cluster_act_means = [tf_act[clusters == ct].mean() for ct in unique_clusters]
            cluster_expr_means = [tf_expr[clusters == ct].mean() for ct in unique_clusters]

            if np.std(cluster_act_means) > 0 and np.std(cluster_expr_means) > 0:
                r, p = spearmanr(cluster_act_means, cluster_expr_means)
                tf_expr_corrs.append({'TF': tf, 'spearman_r': r, 'pvalue': p,
                                     'mean_activity': np.mean(tf_act),
                                     'mean_expression': np.mean(tf_expr)})

    if tf_expr_corrs:
        corr_df = pd.DataFrame(tf_expr_corrs)
        corr_df.to_csv(os.path.join(qc_dir, "tf_activity_vs_expression.csv"), index=False)

        median_r = corr_df['spearman_r'].median()
        pos_frac = (corr_df['spearman_r'] > 0).mean()
        sig_pos = ((corr_df['spearman_r'] > 0.3) & (corr_df['pvalue'] < 0.05)).sum()
        summary['tf_expr_corr_median'] = float(median_r)
        summary['tf_expr_corr_positive_frac'] = float(pos_frac)

        print(f"  TFs with expression data: {len(corr_df)}")
        print(f"  Median Spearman r (cluster means): {median_r:.4f}")
        print(f"  Fraction positive: {pos_frac:.1%}")
        print(f"  Significantly positive (r>0.3, p<0.05): {sig_pos}")
        print(f"  Note: Low correlation is expected — TF activity reflects")
        print(f"  regulatory impact on targets, not expression of the TF itself.")

        fig, axes = plt.subplots(1, 2, figsize=(10, 6))

        ax = axes[0]
        ax.hist(corr_df['spearman_r'], bins=50, color='#3498db', edgecolor='black',
               linewidth=0.5, alpha=0.8, density=True)
        from scipy.stats import gaussian_kde
        kde_x = np.linspace(corr_df['spearman_r'].min() - 0.1, corr_df['spearman_r'].max() + 0.1, 200)
        kde_y = gaussian_kde(corr_df['spearman_r'])(kde_x)
        ax.plot(kde_x, kde_y, color='black', linewidth=2, label='KDE')
        ax.axvline(x=0, color='red', linestyle='--', alpha=0.7)
        ax.axvline(x=median_r, color='green', linestyle='-', linewidth=2,
                  label=f'Median = {median_r:.3f}')
        ax.set_xlabel('Spearman r (cluster-level)')
        ax.set_ylabel('Density')
        ax.set_title('TF Activity vs Expression Correlation')
        ax.legend()

        ax = axes[1]
        top_pos = corr_df.nlargest(15, 'spearman_r')
        top_neg = corr_df.nsmallest(15, 'spearman_r')
        show_tfs = pd.concat([top_pos, top_neg]).sort_values('spearman_r')
        colors = ['#F8766C' if v < 0 else '#00BEC4' for v in show_tfs['spearman_r']]
        ax.barh(range(len(show_tfs)), show_tfs['spearman_r'], color=colors,
               edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(show_tfs)))
        ax.set_yticklabels(show_tfs['TF'], fontsize=8)
        ax.set_xlabel('Spearman r')
        ax.set_title('Top/Bottom TF Activity-Expression Correlations')
        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

        plt.tight_layout()
        plt.savefig(os.path.join(qc_dir, "tf_activity_vs_expression.png"), dpi=200, bbox_inches='tight')
        plt.close()

    # ================================================================
    # 2. TF ACTIVITY SPACE STRUCTURE
    # ================================================================
    print("\n2. TF ACTIVITY SPACE STRUCTURE")

    n_sil = min(5000, tf_matrix.shape[0])
    sil_idx = np.random.choice(tf_matrix.shape[0], n_sil, replace=False)

    sil_score = silhouette_score(tf_matrix[sil_idx], clusters[sil_idx])
    summary['silhouette_clusters_tf_space'] = float(sil_score)
    print(f"  Silhouette (clusters in TF space, n={n_sil}): {sil_score:.4f}")

    if conditions is not None:
        sil_cond = silhouette_score(tf_matrix[sil_idx], conditions[sil_idx])
        summary['silhouette_treatment_tf_space'] = float(sil_cond)
        print(f"  Silhouette (treatment in TF space, n={n_sil}): {sil_cond:.4f}")

    sil_samples_arr = silhouette_samples(tf_matrix[sil_idx], clusters[sil_idx])
    cluster_sil = {}
    for ct in unique_clusters:
        ct_mask = clusters[sil_idx] == ct
        if ct_mask.sum() > 0:
            cluster_sil[ct] = sil_samples_arr[ct_mask].mean()

    sil_series = pd.Series(cluster_sil).sort_values(ascending=True)
    sil_series.to_csv(os.path.join(qc_dir, "silhouette_per_cluster.csv"))

    fig, ax = plt.subplots(figsize=(6, max(8, len(sil_series) * 0.3)))
    colors_sil = ['#e74c3c' if v < 0 else '#f39c12' if v < 0.1 else '#2ecc71'
                  for v in sil_series.values]
    ax.barh(range(len(sil_series)), sil_series.values, color=colors_sil,
           edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(sil_series)))
    ax.set_yticklabels(sil_series.index, fontsize=8)
    ax.set_xlabel('Silhouette Score')
    ax.set_title('Per-Cluster Silhouette in TF Activity Space')
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.axvline(x=sil_score, color='blue', linestyle='--', alpha=0.5,
              label=f'Global mean = {sil_score:.3f}')
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(qc_dir, "silhouette_per_cluster.png"), dpi=330, bbox_inches='tight')
    plt.close()

    # UMAP of TF activity space
    print("  Computing UMAP of TF activity space...")
    tf_adata = ad.AnnData(X=tf_matrix, obs=processed_adata.obs.copy())
    sc.pp.neighbors(tf_adata, use_rep='X', n_neighbors=30)
    sc.tl.umap(tf_adata)

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    sc.pl.umap(tf_adata, color=cluster_col, ax=axes[0], show=False,
               title='TF Activity Space — Clusters', legend_loc='on data',
               legend_fontsize=6, legend_fontoutline=1)
    if conditions is not None:
        sc.pl.umap(tf_adata, color=condition_col, ax=axes[1], show=False,
                   title='TF Activity Space — Treatment',
                   palette={'PBS': '#3498db', 'sCD83': '#e74c3c'})
    plt.tight_layout()
    plt.savefig(os.path.join(qc_dir, "tf_activity_umap.png"), dpi=200, bbox_inches='tight')
    plt.close()

    if conditions is not None:
        print("\n  Per-cluster treatment composition:")
        mixing = {}
        for ct in unique_clusters:
            ct_mask = clusters == ct
            n_a = (conditions[ct_mask] == 'PBS').sum()
            n_b = (conditions[ct_mask] == 'sCD83').sum()
            total = n_a + n_b
            p_a = n_a / total if total > 0 else 0
            p_b = n_b / total if total > 0 else 0
            entropy = 0
            if p_a > 0: entropy -= p_a * np.log2(p_a)
            if p_b > 0: entropy -= p_b * np.log2(p_b)
            mixing[ct] = {'PBS': n_a, 'sCD83': n_b, 'frac_PBS': round(p_a, 3),
                          'mixing_entropy': round(entropy, 3)}
        mix_df = pd.DataFrame(mixing).T
        mix_df.to_csv(os.path.join(qc_dir, "treatment_mixing_per_cluster.csv"))
        for ct, row in mix_df.iterrows():
            flag = " << imbalanced" if row['mixing_entropy'] < 0.5 else ""
            print(f"    {ct:30s}: PBS={int(row['PBS']):5d}, sCD83={int(row['sCD83']):5d}, "
                  f"entropy={row['mixing_entropy']:.3f}{flag}")

    # ================================================================
    # 3. W_POSTERIORS ANALYSIS (cleaned)
    # ================================================================
    print("\n3. W_POSTERIORS ANALYSIS")

    W_post = processed_adata.uns.get('W_posteriors_per_cluster', None)
    if W_post is not None:
        print(f"  Found W_posteriors for {len(W_post)} clusters")

        cluster_tf_importance = {}
        for ct, W in W_post.items():
            W_abs = np.abs(W)
            mean_importance = W_abs.mean(axis=0) if W.ndim > 1 else W_abs
            cluster_tf_importance[ct] = mean_importance

        print("\n  Top 5 TFs by learned importance per cluster:")
        for ct, importance in cluster_tf_importance.items():
            if len(importance) != len(tf_names):
                continue
            valid_mask = [is_valid_gene_name(tf_names[i]) for i in range(len(tf_names))]
            valid_idx = [i for i in range(len(tf_names)) if valid_mask[i]]
            valid_importance = importance[valid_idx]
            top_local_idx = np.argsort(valid_importance)[-5:][::-1]
            top_tfs = [(tf_names[valid_idx[j]], valid_importance[j]) for j in top_local_idx]
            top_str = ', '.join([f"{t[0]}({t[1]:.3f})" for t in top_tfs])
            print(f"    {ct:30s}: {top_str}")

        if len(cluster_tf_importance) >= 2:
            valid_tf_mask = [is_valid_gene_name(tf) for tf in tf_names]
            valid_tf_indices = [i for i, v in enumerate(valid_tf_mask) if v]
            valid_tf_labels = [tf_names[i] for i in valid_tf_indices]

            sorted_cts = sorted([ct for ct in cluster_tf_importance.keys()
                                if len(cluster_tf_importance[ct]) == len(tf_names)])
            imp_matrix = np.array([cluster_tf_importance[ct][valid_tf_indices] for ct in sorted_cts])

            if imp_matrix.shape[0] >= 2 and imp_matrix.shape[1] >= 2:
                tf_variance = imp_matrix.var(axis=0)
                n_show = min(40, len(valid_tf_labels))
                top_var_local = np.argsort(tf_variance)[-n_show:]
                top_var_labels = [valid_tf_labels[i] for i in top_var_local]

                sub_matrix = imp_matrix[:, top_var_local]
                sub_df = pd.DataFrame(sub_matrix, index=sorted_cts, columns=top_var_labels)

                fig_h = max(8, sub_df.shape[0] * 0.3)
                sp = sns.clustermap(sub_df, figsize=(14, fig_h), cmap='viridis',
                                   standard_scale=1, method='complete',
                                   cbar_kws={'label': 'Learned TF Importance\n(normalized)'})
                sp.ax_heatmap.set_xticklabels(sp.ax_heatmap.get_xticklabels(),
                                              fontsize=8, rotation=45, ha='right')
                sp.ax_heatmap.set_yticklabels(sp.ax_heatmap.get_yticklabels(), fontsize=9)
                sp.fig.suptitle('Top Variable TFs by Learned Importance (W_posteriors)',
                               fontsize=13, y=1.02)
                plt.savefig(os.path.join(qc_dir, "W_posteriors_top_variable_tfs.png"),
                           dpi=200, bbox_inches='tight')
                plt.close()
                print(f"  Saved W_posteriors heatmap ({n_show} TFs x {len(sorted_cts)} clusters)")

    # ================================================================
    # 4. EFFECT SIZE ANALYSIS
    # ================================================================
    print("\n4. EFFECT SIZE ANALYSIS")

    diff_path = os.path.join(save_dir, "differential_tf_activity_full.csv")
    if os.path.exists(diff_path):
        diff_df = pd.read_csv(diff_path)
        sig_df = diff_df[diff_df['significant'] == True]

        print(f"  Total significant TF-cell type pairs: {len(sig_df)}")
        print(f"  Upregulated: {(sig_df['difference'] > 0).sum()}")
        print(f"  Downregulated: {(sig_df['difference'] < 0).sum()}")
        print(f"  Mean |difference|: {sig_df['difference'].abs().mean():.4f}")
        print(f"  Mean |Cohen's d|: {sig_df['cohens_d'].abs().mean():.4f}")

        summary['n_significant_pairs'] = len(sig_df)
        summary['n_upregulated'] = int((sig_df['difference'] > 0).sum())
        summary['n_downregulated'] = int((sig_df['difference'] < 0).sum())

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # Volcano
        ax = axes[0]
        nonsig = diff_df[~diff_df['significant']]
        sig_up = sig_df[sig_df['difference'] > 0]
        sig_dn = sig_df[sig_df['difference'] < 0]
        ax.scatter(nonsig['difference'], -np.log10(nonsig['p_adjusted'] + 1e-300),
                  s=2, alpha=0.1, color='gray', rasterized=True)
        ax.scatter(sig_up['difference'], -np.log10(sig_up['p_adjusted'] + 1e-300),
                  s=8, alpha=0.5, color='#e74c3c', label=f'UP ({len(sig_up)})')
        ax.scatter(sig_dn['difference'], -np.log10(sig_dn['p_adjusted'] + 1e-300),
                  s=8, alpha=0.5, color='#3498db', label=f'DOWN ({len(sig_dn)})')
        ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
        ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(x=-0.5, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('TF Activity Difference (sCD83 - PBS)')
        ax.set_ylabel('-log10(FDR)')
        ax.set_title('TF Activity Volcano (all cell types)')
        ax.legend(fontsize=9)

        # Per cell type counts
        ax = axes[1]
        ct_counts = sig_df.groupby('cell_type').agg(
            up=('difference', lambda x: (x > 0).sum()),
            down=('difference', lambda x: (x < 0).sum())
        ).sort_values('up', ascending=True)
        y_pos = range(len(ct_counts))
        ax.barh(y_pos, ct_counts['up'], color='#e74c3c', alpha=0.8, label='UP')
        ax.barh(y_pos, -ct_counts['down'], color='#3498db', alpha=0.8, label='DOWN')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(ct_counts.index, fontsize=7)
        ax.set_xlabel('Number of Significant TFs')
        ax.set_title('Significant TFs per Cell Type')
        ax.axvline(x=0, color='black', linewidth=0.5)
        ax.legend(fontsize=8)

        # Cohen's d distribution
        ax = axes[2]
        ax.hist(sig_df['cohens_d'], bins=40, color='#9b59b6', edgecolor='black',
               linewidth=0.5, alpha=0.8)
        ax.axvline(x=0, color='gray', linestyle='--')
        ax.axvline(x=sig_df['cohens_d'].median(), color='red', linewidth=2,
                  label=f"Median d = {sig_df['cohens_d'].median():.3f}")
        ax.set_xlabel("Cohen's d")
        ax.set_ylabel('Count')
        ax.set_title('Effect Size Distribution (significant TFs)')
        ax.legend(fontsize=9)

        plt.tight_layout()
        plt.savefig(os.path.join(qc_dir, "effect_size_analysis.png"), dpi=200, bbox_inches='tight')
        plt.close()

    # ================================================================
    # 5. PER-CLUSTER DIFFERENTIAL SUMMARY
    # ================================================================
    print("\n5. PER-CLUSTER DIFFERENTIAL SUMMARY")

    if os.path.exists(diff_path):
        diff_df = pd.read_csv(diff_path)
        sig_df = diff_df[diff_df['significant'] == True]

        cluster_summary = []
        for ct in unique_clusters:
            ct_sig = sig_df[sig_df['cell_type'] == ct]
            ct_all = diff_df[diff_df['cell_type'] == ct]
            n_up = (ct_sig['difference'] > 0).sum()
            n_down = (ct_sig['difference'] < 0).sum()

            # Top upregulated TF
            top_up = ct_sig[ct_sig['difference'] > 0].nlargest(1, 'difference')
            top_up_name = top_up['TF'].values[0] if len(top_up) > 0 else '-'
            top_up_diff = top_up['difference'].values[0] if len(top_up) > 0 else 0

            # Top downregulated TF
            top_dn = ct_sig[ct_sig['difference'] < 0].nsmallest(1, 'difference')
            top_dn_name = top_dn['TF'].values[0] if len(top_dn) > 0 else '-'
            top_dn_diff = top_dn['difference'].values[0] if len(top_dn) > 0 else 0

            # Mean absolute effect among significant
            mean_effect = ct_sig['difference'].abs().mean() if len(ct_sig) > 0 else 0

            # Get cell counts
            ct_mask = clusters == ct
            if conditions is not None:
                n_pbs = (conditions[ct_mask] == 'PBS').sum()
                n_scd83 = (conditions[ct_mask] == 'sCD83').sum()
            else:
                n_pbs = n_scd83 = 0

            cluster_summary.append({
                'cell_type': ct, 'n_PBS': n_pbs, 'n_sCD83': n_scd83,
                'n_sig_up': n_up, 'n_sig_down': n_down,
                'top_up_TF': top_up_name, 'top_up_diff': round(top_up_diff, 3),
                'top_down_TF': top_dn_name, 'top_down_diff': round(top_dn_diff, 3),
                'mean_abs_effect': round(mean_effect, 3)
            })

        cs_df = pd.DataFrame(cluster_summary)
        cs_df = cs_df.sort_values('mean_abs_effect', ascending=False)
        cs_df.to_csv(os.path.join(qc_dir, "per_cluster_differential_summary.csv"), index=False)

        print(f"  {'Cell Type':30s} {'UP':>4s} {'DN':>4s} {'|Eff|':>6s} {'Top UP':>12s} {'Top DN':>12s}")
        print(f"  {'-'*80}")
        for _, row in cs_df.iterrows():
            print(f"  {row['cell_type']:30s} {row['n_sig_up']:4d} {row['n_sig_down']:4d} "
                  f"{row['mean_abs_effect']:6.3f} {row['top_up_TF']:>12s} {row['top_down_TF']:>12s}")

    # ================================================================
    # SUMMARY
    # ================================================================
    print(f"\n{'='*50}")
    print("QUALITY SUMMARY")
    print(f"{'='*50}")

    summary_df = pd.Series(summary)
    summary_df.to_csv(os.path.join(qc_dir, "quality_summary.csv"))

    for k, v in summary.items():
        print(f"  {k}: {v}")

    print(f"\n  All outputs saved to: {qc_dir}/")
    return summary



################
# VISUALIZATION
################

def create_enrichment_dotplot(up_results, down_results, save_dir,
                              max_terms=30, max_celltypes=35,
                              custom_celltype_order=None, range_colors=None,
                              show_all_celltypes=True):
    """
    Dotplot for enrichment results. Filters generic GO terms,
    caps at max_terms, handles both Enrichr and GSEA output.
    """
    print("CREATING ENRICHMENT VISUALIZATIONS")
    if not up_results and not down_results:
        print("No results to plot"); return

    # Generic GO terms to exclude — these appear for any TF list and add no insight
    EXCLUDE_PATTERNS = [
        'regulation of transcription',
        'regulation of dna-templated transcription',
        'regulation of nucleic acid-templated transcription',
        'regulation of rna metabolic process',
        'regulation of gene expression',
        'positive regulation of transcription',
        'negative regulation of transcription',
        'positive regulation of dna-templated',
        'negative regulation of dna-templated',
        'positive regulation of nucleic acid',
        'negative regulation of nucleic acid',
        'regulation of transcription by rna polymerase',
        'positive regulation of transcription by rna polymerase',
        'negative regulation of transcription by rna polymerase',
        'positive regulation of macromolecule metabolic',
        'positive regulation of nitrogen compound',
        'regulation of macromolecule metabolic',
        'positive regulation of biosynthetic process',
        'negative regulation of biosynthetic process',
        'positive regulation of cellular biosynthetic',
        'negative regulation of cellular biosynthetic',
        'regulation of biosynthetic process',
        'positive regulation of macromolecule biosynthetic',
        'negative regulation of macromolecule biosynthetic',
    ]

    def _is_generic(term_str):
        t = term_str.lower()
        return any(pat in t for pat in EXCLUDE_PATTERNS)

    plot_data = []

    def _extract_rows(results_dict, default_direction):
        rows = []
        for ct, df in results_dict.items():
            for _, r in df.iterrows():
                # Resolve term
                term = str(r.get('Term', r.get('Name', r.get('term', 'unknown'))))
                if '__' in term:
                    term = term.split('__', 1)[-1]
                term = term[:80]

                if _is_generic(term):
                    continue

                # Direction
                direction = default_direction
                if 'direction' in r.index and pd.notna(r['direction']):
                    direction = r['direction']

                # P-value
                pval = 0.05
                for col in ['Adjusted P-value', 'FDR', 'Adjusted_P-value']:
                    if col in r.index and pd.notna(r[col]):
                        pval = float(r[col])
                        break

                # Overlap count for dot size
                overlap = 1
                if 'n_lead_genes' in r.index and pd.notna(r['n_lead_genes']) and int(r['n_lead_genes']) > 0:
                    overlap = int(r['n_lead_genes'])
                elif 'Overlap' in r.index and pd.notna(r['Overlap']):
                    ov_str = str(r['Overlap'])
                    if '/' in ov_str:
                        try: overlap = int(ov_str.split('/')[0])
                        except: pass
                if overlap <= 1 and 'n_input_tfs' in r.index and pd.notna(r['n_input_tfs']):
                    overlap = int(r['n_input_tfs'])
                overlap = max(overlap, 1)

                rows.append({
                    'cell_type': ct, 'term': term,
                    'pval': pval, 'overlap': overlap,
                    'direction': direction
                })
        return rows

    plot_data.extend(_extract_rows(up_results, 'UP'))
    plot_data.extend(_extract_rows(down_results, 'DOWN'))

    if not plot_data:
        print("No data to plot after filtering generic terms"); return

    pdf = pd.DataFrame(plot_data)
    print(f"  {len(pdf)} data points after removing generic GO terms")

    # Rank terms by significance × frequency across cell types
    ts = pdf.groupby('term').agg(
        avg_sig=('pval', lambda x: -np.log10(x + 1e-300).mean()),
        freq=('cell_type', 'nunique')
    ).reset_index()
    ts['score'] = ts['avg_sig'] * ts['freq']
    sel = ts.nlargest(max_terms, 'score')['term'].tolist()
    pdf = pdf[pdf['term'].isin(sel)]

    if custom_celltype_order:
        cts = (custom_celltype_order if show_all_celltypes else
               [c for c in custom_celltype_order if c in pdf['cell_type'].unique()])[:max_celltypes]
    else:
        cts = sorted(pdf['cell_type'].unique())[:max_celltypes]

    # Figure sizing
    n_terms = len(sel)
    fig_w = max(6, len(cts) * 0.4 + 4)
    fig_h = max(3, n_terms * 0.3 + 2)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cp = {c: i for i, c in enumerate(cts)}
    tp = {t: i for i, t in enumerate(sel)}

    if range_colors:
        for (s, e), col in range_colors.items():
            if s < len(cts):
                ax.axvspan(s - 0.5, min(e, len(cts)) - 0.5, alpha=0.12, color=col)

    # Dot sizing: scale so largest dot is clearly visible
    max_ov = pdf['overlap'].max() if len(pdf) > 0 else 1
    size_scale = 300 / max(max_ov, 1)
    min_dot = 30

    for _, r in pdf.iterrows():
        if r['cell_type'] not in cp:
            continue
        dot_size = max(r['overlap'] * size_scale, min_dot)
        ax.scatter(cp[r['cell_type']], tp[r['term']],
                   s=dot_size,
                   c='#00BEC4' if r['direction'] == 'UP' else '#F8766C',
                   alpha=0.75, edgecolors='black', linewidth=0.5)

    ax.set_xticks(range(len(cts)))
    ax.set_xticklabels(cts, rotation=90, ha='center', fontsize=14, va='top', fontweight="bold")
    ax.set_yticks(range(len(sel)))
    ax.set_yticklabels(sel, fontsize=12, fontweight="bold")
    ax.set_xlabel('Cell Type', fontsize=16, fontweight='bold')
    ax.set_ylabel('Enriched GO Term', fontsize=12, fontweight='bold')
    ax.set_title('TF Enrichment (sCD83 vs PBS)\nFiltered: top 30 specific terms',
                 fontsize=16, fontweight='bold', pad=20)
    ax.grid(True, alpha=0.2, linestyle='--')
    ax.set_axisbelow(True)
    ax.set_xlim(-0.5, len(cts) - 0.5)
    ax.set_ylim(-0.5, len(sel) - 0.5)

    # Unified legend combining direction and overlap size
    size_vals = [2, 4, 6]

    legend_handles = []
    legend_labels = []

    # Direction entries
    um = plt.scatter([], [], s=120, c='#00BEC4', alpha=0.75,
                     edgecolors='black', linewidth=0.5)
    dm = plt.scatter([], [], s=120, c='#F8766C', alpha=0.75,
                     edgecolors='black', linewidth=0.5)
    legend_handles += [um, dm]
    legend_labels += ['Upregulated', 'Downregulated']

    # Blank spacer between the two groups
    spacer = plt.scatter([], [], s=0, c='none', edgecolors='none')
    legend_handles.append(spacer)
    legend_labels.append('')

    # Overlap / size entries
    for s in size_vals:
        h = plt.scatter([], [], s=max(s * size_scale, min_dot) * 1.4,
                        c='gray', alpha=0.7, edgecolors='black', linewidth=0.5)
        legend_handles.append(h)
        legend_labels.append(f'{s}+ TFs' if s == size_vals[-1] else f'{s} TFs')

    ax.legend(legend_handles, legend_labels,
              title='Direction / Overlap',
              loc='upper left', bbox_to_anchor=(1.02, 1),
              fontsize=12, title_fontsize=14,
              frameon=True, framealpha=0.9,
              borderpad=1, labelspacing=1.4,
              handletextpad=1.0)

    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "enrichment_dotplot.png"),
                dpi=330, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(save_dir, "enrichment_dotplot.pdf"),
                bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved enrichment dotplot ({n_terms} terms x {len(cts)} cell types)")


################
# MAIN
################

if __name__ == "__main__":
    args = parse_arguments()

    if args.minimal:
        print("\n" + "="*60)
        print("  MINIMAL TEST-RUN MODE"); print("="*60 + "\n")
    if args.output_dir:
        save_dir = args.output_dir; os.makedirs(save_dir, exist_ok=True)

    print("scRegulate: WEIGHTED GRN + DIFFERENTIAL TF + ONTOLOGY PIPELINE")

    if args.minimal:
        effective_tp, effective_fp = apply_minimal_overrides(TRAINING_PARAMS, FINETUNE_PARAMS, MINIMAL_PARAMS)
    else:
        effective_tp, effective_fp = TRAINING_PARAMS, FINETUNE_PARAMS

    # Track the GRN dataframe for enrichment use downstream
    grn_for_enrichment = None

    if args.load_model:
        print(f"\n{'='*60}\nLOADING PRE-TRAINED MODEL FROM: {args.load_model}\n{'='*60}")
        model, processed_adata, tf_activities, GRN = load_model_outputs(
            load_dir=args.load_model, adata_file=args.adata_file,
            tf_file=args.tf_file, grn_file=args.grn_file, model_file=args.model_file)

        # Resolve GRN for enrichment from loaded outputs
        if isinstance(GRN, pd.DataFrame) and {'source', 'target'}.issubset(GRN.columns):
            grn_for_enrichment = GRN
        else:
            grn_path_parquet = os.path.join(
                args.load_model, "net_weighted_grn.parquet")
            if os.path.exists(grn_path_parquet):
                grn_for_enrichment = pd.read_parquet(grn_path_parquet)
                print(f"Loaded GRN for enrichment from: {grn_path_parquet}")
            else:
                # Fallback: check save_dir
                grn_path_parquet_alt = os.path.join(save_dir, "net_weighted_grn.parquet")
                if os.path.exists(grn_path_parquet_alt):
                    grn_for_enrichment = pd.read_parquet(grn_path_parquet_alt)
                    print(f"Loaded GRN for enrichment from: {grn_path_parquet_alt}")
                else:
                    print("WARNING: No usable GRN found for target-gene enrichment. "
                          "Enrichment will be skipped.")
    else:
        print(f"\n{'='*60}\nTRAINING NEW MODEL\n{'='*60}")
        adata = sc.read_h5ad(adata_path); print(f"Loaded adata: {adata.shape}")
        rna_data, adata_filtered = prepare_data(adata, cells_to_remove=cells_to_remove)
        if args.minimal:
            rna_data = subsample_cells(rna_data, n_cells=MINIMAL_PARAMS['n_cells'],
                                       cluster_col=DIFF_PARAMS['cluster_col'])
        skin_grn_raw = pd.read_parquet(grn_path); print(f"Loaded GRN: {skin_grn_raw.shape}")

        # BUILD EXPRESSION-WEIGHTED GRN
        net, tf_weights = build_weighted_grn(
            adata=adata_filtered, skin_grn_raw=skin_grn_raw, save_dir=save_dir,
            min_weight=MIN_WEIGHT, max_weight=MAX_WEIGHT, scale_to_max=SCALE_TO_MAX)

        grn_for_enrichment = net  # source/target/weight DataFrame

        if args.train:
            model, processed_adata = train_scregulate_model(
                adata=rna_data, net=net, save_dir=save_dir,
                training_params=effective_tp, finetune_params=effective_fp)
        else:
            pp = os.path.join(save_dir, "processed_adata.h5ad")
            if os.path.exists(pp):
                processed_adata = sc.read_h5ad(pp)
            else:
                print("No processed data found. Use --train."); exit(1)

    print(f"\n{'='*60}\nRUNNING DOWNSTREAM ANALYSIS\n{'='*60}")

    if not args.skip_plots:
        analyze_cluster_similarity(processed_adata=processed_adata, save_dir=save_dir,
            subset_clusters=CUSTOM_CELLTYPE_ORDER[:25])

    results_df = differential_tf_activity(processed_adata=processed_adata,
                                          save_dir=save_dir, diff_params=DIFF_PARAMS)

    if not args.skip_plots:
        create_differential_heatmap(processed_adata=processed_adata,
            save_dir=save_dir, heatmap_params=HEATMAP_PARAMS,
            custom_order=CUSTOM_CELLTYPE_ORDER)

    if not args.skip_plots:
        quality = post_hoc_model_quality(
        model=model if model is not None else None,
        processed_adata=processed_adata,
        tf_activities=tf_activities if 'tf_activities' in dir() else None,
        GRN=GRN if 'GRN' in dir() else None,
        save_dir=save_dir,
        cluster_col=DIFF_PARAMS['cluster_col'],
        condition_col=DIFF_PARAMS['condition_col'])

    if not args.skip_enrichment:
        up_r, dn_r = run_ontology_enrichment(
            results_df=results_df, save_dir=save_dir,
            enrichment_params=ENRICHMENT_PARAMS)
    if not args.skip_plots and up_r is not None:
        create_enrichment_dotplot(up_results=up_r, down_results=dn_r, save_dir=save_dir,
            custom_celltype_order=CUSTOM_CELLTYPE_ORDER, range_colors=ENRICHMENT_RANGE_COLORS)

    print(f"\n{'='*60}")
    print("MINIMAL TEST-RUN COMPLETE" if args.minimal else "PIPELINE COMPLETE")
    print(f"{'='*60}")
