#!/usr/bin/env python3
"""
plot_marker_genes.py

Recreates per-cell expression bar plots (one row per gene) with cluster-mean
step-line overlays and a top colour bar, following the layout logic of the
original notebook code.

Inputs
------
--h5ad    : path to an AnnData .h5ad file (must contain raw counts or
            normalised expression in adata.X or a named layer).
--layer   : layer to pull expression from (default: uses adata.X).
--outdir  : root output directory.

The script:
1.  Derives ``mapping_cell_type`` from ``adata.obs["broad_cluster"]``.
2.  For each mapping_cell_type that has ≥1 cell, produces a PDF with one
    row per marker gene defined in ``gene_groups``.
"""

import argparse
import os
import sys
from collections import Counter, OrderedDict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scanpy as sc

# ── font setup: prefer Arial, fall back to DejaVu Sans ──────────────────────
from matplotlib.font_manager import fontManager
_available_families = {f.name for f in fontManager.ttflist}
FONT_FAMILY = "Arial" if "Arial" in _available_families else "DejaVu Sans"


# ── gene groups ──────────────────────────────────────────────────────────────

GENE_GROUPS: dict[str, list[str]] = {
    "Keratinocytes": [
        "KRT14", "CXCL8", "KRT6A", "S100A2", "COL17A1", "TP63", "KRT85", "KRT19", "KLK10",
    ],
    "Fibroblasts": [
        "PLXDC2", "RGS5", "MMP1", "VCAM1", "CRYAB", "TAGLN", "NFKBIA", "TXN", "CXCL2",
    ],
    "Endothelial": [
        "PECAM1", "VWF", "VEGFA", "BICC1", "CCL21", "DLL4",
    ],
    "Immune": [
        "CD3D", "KLRB1", "MRC1", "CD83", "KIT", "FOXP3", "CD8A", "LAMP3",
        "IGKC", "PAX5", "ANK3",
    ],
    "Remaining": [
        "MLANA", "SOX5", "DCT", "NGFR", "PMP2", "SOX10",
    ],
}


# ── helper functions ─────────────────────────────────────────────────────────

def return_unique(seq):
    """Return unique values of *seq* in first-occurrence order."""
    seen = set()
    out = []
    for v in seq:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def get_cluster_ticks(cl):
    """Return x-positions at the midpoint of each contiguous cluster block."""
    ticks = []
    cnt = 0
    for i in return_unique(cl):
        l = Counter(cl)[i]
        ticks.append(cnt + l / 2)
        cnt += l
    return ticks


def clean_axis(ax):
    """Remove all spines and ticks from an axis."""
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)


def remove_ticks(ax, linewidth=1):
    """Remove tick marks but keep spines."""
    ax.tick_params(
        left=False, bottom=False, top=False, right=False,
        labelleft=False, labelbottom=False, labeltop=False, labelright=False,
    )


def build_cluster_cmap(unique_clusters):
    """
    Assign a distinct colour to each cluster label using a fixed 20-colour
    palette (matching the standard ArchR / custom R colour map).

    Clusters are sorted alphabetically and assigned colours positionally:
    the first cluster alphabetically receives colour 1 (#D51F26), the
    second receives colour 2 (#272E6A), and so on.  If there are more
    than 20 clusters the palette cycles.
    """
    _PALETTE_COLOURS = [
        "#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B",
        "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4",
        "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8",
        "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", "#3D3D3D",
    ]

    sorted_clusters = sorted(unique_clusters, key=lambda x: str(x).casefold())
    cmap = {
        cl: _PALETTE_COLOURS[i % len(_PALETTE_COLOURS)]
        for i, cl in enumerate(sorted_clusters)
    }
    return cmap


def order_genes_and_clusters(dataset, cell_groups, genes):
    """
    Jointly reorder genes (rows) and fine clusters (x-axis blocks) to produce
    a diagonal-like pattern: genes whose expression peaks in early clusters
    appear at the top, cluster-specific genes sink to the bottom-right.

    Algorithm
    ---------
    1.  Build a (genes × clusters) mean-expression matrix.
    2.  For each gene, find the cluster with peak mean expression.
    3.  Determine cluster order: walk through genes in the input list; the
        first time a cluster appears as a peak, it is appended to the
        cluster order.  Clusters that never peak for any gene are appended
        at the end alphabetically.
    4.  Re-sort genes by (peak-cluster position, -peak mean) so that within
        the same peak cluster, the gene with the highest mean comes first.
    5.  Re-sort cells by the new cluster order.

    Returns
    -------
    genes_ordered   : list[str]
    cell_groups_ordered : pd.Series  (sorted so clusters are contiguous in
                         the new order)
    """
    unique_clusters = return_unique(cell_groups)

    # mean expression matrix: genes × clusters
    mean_mat = pd.DataFrame(index=genes, columns=unique_clusters, dtype=float)
    for cl in unique_clusters:
        ix = cell_groups[cell_groups == cl].index
        mean_mat[cl] = dataset.loc[genes, ix].mean(axis=1)

    # peak cluster per gene
    peak_cluster = mean_mat.idxmax(axis=1)
    peak_value = mean_mat.max(axis=1)

    # derive cluster order from gene order: first-seen peak cluster
    cluster_order = []
    seen = set()
    for g in genes:
        pc = peak_cluster[g]
        if pc not in seen:
            cluster_order.append(pc)
            seen.add(pc)
    # append any clusters that were never a peak
    for cl in sorted(unique_clusters):
        if cl not in seen:
            cluster_order.append(cl)
            seen.add(cl)

    # map cluster → position in the new order
    cl_pos = {cl: i for i, cl in enumerate(cluster_order)}

    # sort genes by (cluster position of peak, -peak mean value)
    gene_sort_key = [(cl_pos[peak_cluster[g]], -peak_value[g], g) for g in genes]
    gene_sort_key.sort()
    genes_ordered = [t[2] for t in gene_sort_key]

    # re-sort cells so clusters appear in cluster_order
    cl_rank = cell_groups.map(cl_pos)
    cell_groups_ordered = cell_groups.iloc[cl_rank.argsort(kind="mergesort").values]

    return genes_ordered, cell_groups_ordered


# ── core plotting function ───────────────────────────────────────────────────

def plot_genes(
    dataset: pd.DataFrame,
    cell_groups: pd.Series,
    genes: list[str],
    cmap: dict,
    outpath: str,
):
    """
    Produce the multi-gene bar plot PDF.

    Parameters
    ----------
    dataset     : DataFrame  (genes × cells), expression values.
    cell_groups : Series     (index = cell barcodes, values = cluster labels),
                  already sorted by cluster.
    genes       : list of gene names to plot (rows).
    cmap        : dict mapping cluster label → colour.
    outpath     : full path for the output PDF.
    """
    if not genes:
        return

    # ── colour list per cell ──
    clist_bar = [cmap[val] for val in cell_groups]

    # ── figure dimensions ──
    bar_size = 0.5
    bar_pad = 0.25
    width = 15
    height = max(len(genes) * 1.1, 1.5)

    # layout margins (in figure-fraction)
    label_left = 3.5 / width      # space for gene names on the left
    legend_right = 4.0 / width    # space for the legend on the right
    plot_left = label_left
    plot_right = 1.0 - legend_right

    fig = plt.figure(facecolor="w", figsize=(width, height))

    # GridSpec for gene rows — single column for the bar plots
    gs1 = plt.GridSpec(
        len(genes), 1,
        hspace=0.00, wspace=0.0,
        top=1 - (bar_pad + bar_size) / height,
        bottom=0,
        left=plot_left,
        right=plot_right,
    )

    # GridSpec for the cluster colour bar at the top
    gs0 = plt.GridSpec(
        1, 1,
        left=plot_left,
        right=plot_right,
        top=1,
        bottom=1 - bar_size / height,
        hspace=0.0, wspace=0.0,
    )

    n_cells = len(dataset.columns)

    for row_idx, g in enumerate(genes):
        # ── cluster means ──
        unique_groups = return_unique(cell_groups)
        mean_tmp = pd.Series(index=unique_groups, dtype=float)
        for gr in unique_groups:
            ix_tmp = cell_groups[cell_groups == gr].index
            mean_tmp[gr] = dataset.loc[g, ix_tmp].mean()

        # main expression axis
        ax = plt.subplot(gs1[row_idx, 0])
        ax.axvspan(0, n_cells, color="#FFFFFF", zorder=0)

        # x-axis
        ax.set_xlim(0, n_cells)
        ax.xaxis.set_ticks(get_cluster_ticks(cell_groups))
        ax.xaxis.set_ticklabels([])

        # y-axis
        y_max = np.ceil(np.nanmax(mean_tmp) * 2)
        if y_max == 0 or np.isnan(y_max):
            y_max = 1
        ax.set_ylim(0, y_max)

        # gene label on the LEFT
        ax.set_ylabel(
            g,
            family=FONT_FAMILY,
            fontsize=36,
            fontweight="bold",
            rotation="horizontal",
            va="center",
            ha="right",
        )
        ax.yaxis.set_label_coords(-0.12, 0.5)
        ax.locator_params(axis="y", nbins=2)

        for tick_pos, tick in enumerate(ax.yaxis.get_major_ticks()):
            if tick_pos % 2 == 1:
                tick.label1.set_family(FONT_FAMILY)
                tick.label1.set_fontsize(25)
            else:
                tick.label1.set_visible(False)

        # per-cell bars
        expr_vals = dataset.loc[g].values
        ax.bar(
            np.arange(n_cells),
            expr_vals,
            color=clist_bar,
            linewidth=0,
            width=1.0,
            rasterized=True,
        )

        # cluster-mean step line
        mean_per_cell = [mean_tmp[val] for val in cell_groups]
        ax.step(
            range(len(cell_groups)),
            mean_per_cell,
            where="mid",
            color="black",
            linewidth=2,
        )

    # ── top colour bar (cluster identity) ──
    ax_top = plt.subplot(gs0[0, 0])
    ax_top.set_xlim(0, n_cells)
    for pos, barcode in enumerate(cell_groups.index):
        ax_top.axvspan(
            xmin=pos, xmax=pos + 1, color=cmap[cell_groups[barcode]]
        )
    remove_ticks(ax_top, linewidth=1)
    clean_axis(ax_top)

    # ── legend axis on the right ──
    unique_groups = return_unique(cell_groups)
    n_clusters = len(unique_groups)

    ax_leg = fig.add_axes([
        plot_right + 0.02,           # x origin: just right of plot area
        max(0.0, 1 - (bar_pad + bar_size) / height - n_clusters * 0.06),  # y origin
        legend_right - 0.04,         # width
        min(1.0, n_clusters * 0.06), # height: scale with number of clusters
    ])
    clean_axis(ax_leg)

    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=cmap[cl], edgecolor="none")
        for cl in unique_groups
    ]
    ax_leg.legend(
        legend_handles,
        unique_groups,
        loc="upper left",
        frameon=False,
        prop={"family": FONT_FAMILY, "size": 22},
        handlelength=1.4,
        handleheight=1.4,
        labelspacing=0.6,
    )

    plt.savefig(
        outpath,
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=150,
    )
    plt.close(fig)
    print(f"  saved → {outpath}")


# ── UMAP highlight plot ──────────────────────────────────────────────────────

def plot_umap(
    adata,
    cell_type: str,
    outpath: str,
    highlight_colour: str = "#E64B35",
    obsm_key: str = "X_umap",
    point_size: float = 5,
    figsize: tuple = (10, 10),
):
    """
    UMAP scatter: cells belonging to *cell_type* are highlighted in a single
    colour; all other cells are drawn in silver behind them.

    Parameters
    ----------
    adata            : full AnnData object (all cells).
    cell_type        : the mapping_cell_type to highlight.
    outpath          : output PDF path.
    highlight_colour : colour for the highlighted cell type.
    obsm_key         : key in adata.obsm holding the 2-D embedding.
    point_size       : marker size for scatter.
    figsize          : (width, height) in inches.
    """
    if obsm_key not in adata.obsm:
        print(f"  WARNING: '{obsm_key}' not found in adata.obsm — skipping UMAP.")
        return

    coords = pd.DataFrame(
        adata.obsm[obsm_key][:, :2],
        index=adata.obs_names,
        columns=["x", "y"],
    )

    fg_mask = adata.obs["mapping_cell_type"] == cell_type
    fg_idx = adata.obs_names[fg_mask]
    bg_idx = adata.obs_names[~fg_mask]

    # ── figure ──
    fig, ax = plt.subplots(figsize=figsize, facecolor="w")

    # axis limits with equal aspect & padding
    x_min, x_max = coords["x"].min(), coords["x"].max()
    y_min, y_max = coords["y"].min(), coords["y"].max()
    x_diff = x_max - x_min
    y_diff = y_max - y_min
    pad = 2.0

    if x_diff >= y_diff:
        mid_y = (y_min + y_max) / 2
        half = x_diff / 2
        xlim = (x_min - pad, x_max + pad)
        ylim = (mid_y - half - pad, mid_y + half + pad)
    else:
        mid_x = (x_min + x_max) / 2
        half = y_diff / 2
        xlim = (mid_x - half - pad, mid_x + half + pad)
        ylim = (y_min - pad, y_max + pad)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")

    # ── background cells (silver, shuffled) ──
    bg_order = bg_idx.to_numpy().copy()
    np.random.seed(42)
    np.random.shuffle(bg_order)
    ax.scatter(
        coords.loc[bg_order, "x"],
        coords.loc[bg_order, "y"],
        s=point_size,
        c="silver",
        edgecolors="none",
        linewidth=0,
        rasterized=True,
    )

    # ── foreground cells (single highlight colour, shuffled) ──
    fg_order = fg_idx.to_numpy().copy()
    np.random.shuffle(fg_order)

    ax.scatter(
        coords.loc[fg_order, "x"],
        coords.loc[fg_order, "y"],
        s=point_size,
        c=highlight_colour,
        edgecolors="none",
        linewidth=0,
        rasterized=True,
    )

    remove_ticks(ax)
    clean_axis(ax)

    plt.savefig(
        outpath,
        format="pdf",
        transparent=True,
        bbox_inches="tight",
        pad_inches=0.1,
        dpi=150,
    )
    plt.close(fig)
    print(f"  saved → {outpath}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Per-cell marker-gene expression bar plots."
    )
    parser.add_argument("--h5ad", required=True, help="Path to .h5ad file.")
    parser.add_argument(
        "--layer", default=None,
        help="Layer to use for expression values (default: adata.X).",
    )
    parser.add_argument(
        "--outdir", default="marker_plots",
        help="Root output directory (default: marker_plots).",
    )
    parser.add_argument(
        "--fine-clust-col", default="fine_clust",
        help="obs column with fine sub-cluster labels (default: fine_clust).",
    )
    parser.add_argument(
        "--no-reorder", action="store_true",
        help="Disable diagonal reordering of genes and clusters. "
             "Keeps gene order as defined in GENE_GROUPS and clusters sorted "
             "alphabetically.",
    )
    args = parser.parse_args()

    # ── load data ──
    print(f"Reading {args.h5ad} ...")
    adata = sc.read_h5ad(args.h5ad)

    # ── derive mapping_cell_type ──
    if "broad_cluster" not in adata.obs.columns:
        sys.exit("ERROR: adata.obs must contain a 'broad_cluster' column.")

    fine_col = args.fine_clust_col
    if fine_col not in adata.obs.columns:
        sys.exit(f"ERROR: adata.obs must contain a '{fine_col}' column.")

    adata.obs["mapping_cell_type"] = (
        adata.obs["broad_cluster"]
        .astype(str)
        .str.replace(r"\.\d+$", "", regex=True)
    )

    # ── check gene availability per group ──
    for grp, glist in GENE_GROUPS.items():
        missing = [g for g in glist if g not in adata.var_names]
        if missing:
            print(f"WARNING [{grp}]: genes not found and will be skipped: {missing}")

    # ── build expression DataFrame (genes × cells) for all available genes ──
    all_genes = list({g for gl in GENE_GROUPS.values() for g in gl if g in adata.var_names})
    if args.layer is not None:
        expr_matrix = adata[:, all_genes].layers[args.layer]
    else:
        expr_matrix = adata[:, all_genes].X

    if hasattr(expr_matrix, "toarray"):
        expr_matrix = expr_matrix.toarray()

    dataset_full = pd.DataFrame(
        expr_matrix.T,
        index=all_genes,
        columns=adata.obs_names,
    )

    # ── iterate over mapping_cell_types ──
    cell_types = sorted(adata.obs["mapping_cell_type"].unique())
    print(f"Found {len(cell_types)} mapping cell types: {cell_types}")

    for cell_type in cell_types:

        # ── pick the matching gene group ──
        if cell_type not in GENE_GROUPS:
            print(f"\nSkipping '{cell_type}' — no gene group defined.")
            continue

        genes = [g for g in GENE_GROUPS[cell_type] if g in dataset_full.index]
        if not genes:
            print(f"\nSkipping '{cell_type}' — none of its marker genes are in the data.")
            continue

        mask = adata.obs["mapping_cell_type"] == cell_type
        barcodes = adata.obs_names[mask]
        if len(barcodes) == 0:
            continue

        print(f"\nProcessing '{cell_type}' ({len(barcodes)} cells) ...")

        # ── use fine_clust for sub-cluster grouping & colouring ──
        cell_groups = adata.obs.loc[barcodes, fine_col].astype(str)

        # sort cells by fine cluster so blocks are contiguous
        cell_groups = cell_groups.sort_values()
        barcodes_sorted = cell_groups.index

        dataset_sub = dataset_full.loc[genes, barcodes_sorted]

        # ── diagonal reordering: genes & clusters ──
        if not args.no_reorder:
            genes, cell_groups = order_genes_and_clusters(
                dataset_sub, cell_groups, genes
            )
            barcodes_sorted = cell_groups.index
            dataset_sub = dataset_full.loc[genes, barcodes_sorted]

        print(f"  gene order : {genes}")
        print(f"  cluster order: {return_unique(cell_groups)}")

        # cluster colour map based on fine sub-clusters
        unique_clusters = return_unique(cell_groups)
        cmap = build_cluster_cmap(unique_clusters)

        # output path
        outdir = os.path.join(args.outdir, cell_type.lower())
        os.makedirs(outdir, exist_ok=True)
        figname = f"{cell_type.lower()}_marker_genes.pdf"
        outpath = os.path.join(outdir, figname)

        plot_genes(dataset_sub, cell_groups, genes, cmap, outpath)

        # ── UMAP highlight plot ──
        umap_figname = f"{cell_type.lower()}_umap.pdf"
        umap_outpath = os.path.join(outdir, umap_figname)
        plot_umap(adata, cell_type, umap_outpath)

    print("\nDone.")


if __name__ == "__main__":
    main()
