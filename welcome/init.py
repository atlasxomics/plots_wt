import anndata
import math
import numpy as np
import pandas as pd
import plotly.express as px
import scanpy as sc
import scipy.sparse as sp


from anndata import AnnData
from pathlib import Path
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from typing import List

from lplots import palettes, submit_widget_state
from lplots.reactive import Signal
from lplots.widgets.button import w_button
from lplots.widgets.checkbox import w_checkbox
from lplots.widgets.column import w_column
from lplots.widgets.grid import w_grid
from lplots.widgets.h5 import w_h5
from lplots.widgets.ldata import w_ldata_picker
from lplots.widgets.multiselect import w_multi_select
from lplots.widgets.plot import w_plot
from lplots.widgets.row import w_row
from lplots.widgets.select import w_select
from lplots.widgets.table import w_table
from lplots.widgets.text import w_text_input, w_text_output

w_text_output(content="# **AnnData H5AD Viewer**")
w_text_output(content="""

This notebook provides interactive tools for **exploratory data analysis** and **figure generation** from AnnData `.h5ad` files. Plotting modules are organized into tabs at the **top of this window**--move between tabs to explore results.

""")

DEFAULT_H5_CATEGORICAL_PALETTE = [
    "#C33530", "#282E66", "#43884A", "#7E2F8A", "#E48341",
    "#FAE64D", "#8E9ECD", "#B570A8", "#E0C3DA", "#9FD3E2",
    "#96C56C", "#E38180", "#9584B9", "#C25434", "#63B9A8",
    "#694D99", "#33707A", "#731F1C", "#D0A970", "#3D3D3D",
]
DEFAULT_CATEGORICAL_PALETTE_NAME = "Default H5 Viewer Palette"
DEFAULT_CONTINUOUS_PALETTE = "PuBu_r"
DEFAULT_CONTINUOUS_PALETTE_NAME = "Default Continuous Palette"
MAX_DEFAULT_CATEGORIES = 30
MAX_PLOT_CATEGORIES = 100
OBS_ID_KEYWORDS = (
    "barcode",
    "barcodes",
    "cell_id",
    "cellid",
    "cell_name",
    "cellname",
    "spot_id",
    "spotid",
    "uuid",
)
CELL_ID_INDEX_KEYWORDS = ("barcode", "cell", "spot", "aaac", "tttg")
OBS_METADATA_HINTS = ("sample", "condition", "cluster")
VAR_FEATURE_HINTS = (
    "gene",
    "gene_id",
    "gene_ids",
    "gene_name",
    "gene_names",
    "gene_symbol",
    "gene_symbols",
    "feature",
    "feature_name",
    "feature_names",
    "symbol",
)
GENE_NAME_PREFIXES = ("ENSG", "ENSMUSG", "ENSGALG", "LOC", "MIR", "RPL", "RPS")
KNOWN_GENE_NAMES = {
    "ACTB",
    "GAPDH",
    "MALAT1",
    "B2M",
    "TPT1",
    "ACTG1",
    "EEF1A1",
    "FTL",
    "FTH1",
    "RPLP0",
    "RPS18",
    "MS4A1",
    "CD3D",
    "CD3E",
    "LYZ",
    "PPBP",
    "NKG7",
}


# Globals ------------------------------------------------------------------

if "new_data_signal" not in globals():
    new_data_signal = Signal(False)
if "choose_subset_signal" not in globals():
    choose_subset_signal = Signal(False)
if "gene_score_done_signal" not in globals():
    gene_score_done_signal = Signal(False)
if "refresh_h5_signal" not in globals():
    refresh_h5_signal = Signal(False)

na_keys = ['barcode', 'on_off', 'row', 'col', 'xcor', 'ycor', 'score']

adata = None
adata_g = None


# Functions ----------------------------------------------------------------

def create_proportion_dataframe(
    adata, group_by, stack_by, return_type="proportion"
):
    """Create a DataFrame for stacked count/proportion plots."""
    count_df = pd.crosstab(adata.obs[group_by], adata.obs[stack_by])

    if return_type == "proportion":
        result_df = count_df.div(count_df.sum(axis=1), axis=0)
    elif return_type == "counts":
        result_df = count_df
    else:
        raise ValueError("Invalid return_type. Use 'proportion' or 'counts'.")

    long_df = result_df.reset_index().melt(
      id_vars=group_by,
      value_name="value",
      var_name=stack_by
    )
    long_df.columns = ["group_by", "stack_by", "value"]

    return long_df


def create_violin_data(adata, group_by, plot_data, data_type="obs"):
    """Create a DataFrame for violin or box plots."""
    if data_type == "obs":
        values = adata.obs[plot_data]
    elif data_type == "feature":
        feature_matrix = adata[:, plot_data].X
        if hasattr(feature_matrix, "toarray"):
            feature_matrix = feature_matrix.toarray()
        values = np.asarray(feature_matrix).ravel()
    else:
        raise ValueError("data_type must be either 'obs' or 'feature'.")

    return pd.DataFrame({
        "group": adata.obs[group_by],
        "value": values,
    })


def drop_obs_column(adata_objects, col_to_drop="orig.ident"):
    """Remove a specified obs column from AnnData objects when present."""
    for adata_obj in adata_objects:
        if col_to_drop in adata_obj.obs.columns:
            adata_obj.obs = adata_obj.obs.drop(columns=[col_to_drop])


def empty_notebook_palettes():
    return {"categorical": [], "continuous": []}


async def get_notebook_palettes():
    try:
        palette_data = await palettes.get()
    except RuntimeError as e:
        if "Event loop is closed" in str(e):
            return empty_notebook_palettes()
        print(f"Unable to load notebook palettes: {e}")
        return empty_notebook_palettes()
    except Exception as e:
        print(f"Unable to load notebook palettes: {e}")
        return empty_notebook_palettes()

    if not isinstance(palette_data, dict):
        return empty_notebook_palettes()

    return {
        kind: [
            palette for palette in palette_data.get(kind, [])
            if palette.get("display_name") and palette.get("colors")
        ]
        for kind in ("categorical", "continuous")
    }


def get_palette_selector_options(
    palette_data,
    kind="categorical",
    fallback_name=DEFAULT_CATEGORICAL_PALETTE_NAME,
):
    palette_names = []
    seen = {fallback_name}

    for palette in palette_data.get(kind, []):
        display_name = palette["display_name"]
        if display_name not in seen:
            palette_names.append(display_name)
            seen.add(display_name)

    return tuple([fallback_name] + palette_names)


def get_selected_palette_colors(
    palette_data,
    selected_name,
    kind="categorical",
    fallback_colors=None,
    fallback_name=DEFAULT_CATEGORICAL_PALETTE_NAME,
):
    if fallback_colors is None:
        fallback_colors = DEFAULT_H5_CATEGORICAL_PALETTE

    if not selected_name or selected_name == fallback_name:
        return fallback_colors

    for palette in palette_data.get(kind, []):
        if palette["display_name"] == selected_name:
            return palette["colors"]

    return fallback_colors


def get_selected_continuous_palette(
    palette_data,
    selected_name,
    fallback_colors=DEFAULT_CONTINUOUS_PALETTE,
    fallback_name=DEFAULT_CONTINUOUS_PALETTE_NAME,
):
    return get_selected_palette_colors(
        palette_data,
        selected_name,
        kind="continuous",
        fallback_colors=fallback_colors,
        fallback_name=fallback_name,
    )


def build_discrete_color_map(categories, colors):
    if not colors:
        colors = DEFAULT_H5_CATEGORICAL_PALETTE

    return {category: colors[i % len(colors)] for i, category in enumerate(categories)}


def get_groups(adata: anndata.AnnData) -> List[str]:
    """Return preferred grouping columns with more than one observed value."""
    groups = []
    for group in ["cluster", "sample", "condition"]:
        if group in adata.obs and len(adata.obs[group].dropna().unique()) > 1:
            groups.append(group)

    return groups


def get_categorical_obs_keys(adata: anndata.AnnData) -> List[str]:
    """Return all categorical-like obs columns."""
    return [
        key for key in adata.obs_keys()
        if (
            pd.api.types.is_object_dtype(adata.obs[key])
            or pd.api.types.is_categorical_dtype(adata.obs[key])
            or pd.api.types.is_bool_dtype(adata.obs[key])
        )
    ]


def get_obs_category_summary(adata: anndata.AnnData) -> List[dict]:
    """Summarize categorical obs columns for plot usability."""
    summary = []
    for key in get_categorical_obs_keys(adata):
        values = adata.obs[key].dropna()
        n_unique = int(values.nunique())
        key_lower = key.lower()
        is_id_like = any(token in key_lower for token in OBS_ID_KEYWORDS)
        default_ok = (
            not is_id_like
            and 2 <= n_unique <= MAX_DEFAULT_CATEGORIES
        )
        plot_ok = (
            not is_id_like
            and 2 <= n_unique <= MAX_PLOT_CATEGORIES
        )

        summary.append({
            "key": key,
            "n_unique": n_unique,
            "is_id_like": is_id_like,
            "default_ok": default_ok,
            "plot_ok": plot_ok,
        })

    return summary


def get_groupable_obs_keys(
    adata: anndata.AnnData,
    *,
    max_categories: int = MAX_PLOT_CATEGORIES,
) -> List[str]:
    """Return categorical obs columns that should be safe to plot/group by."""
    return [
        row["key"]
        for row in get_obs_category_summary(adata)
        if row["plot_ok"] and row["n_unique"] <= max_categories
    ]


def get_excluded_groupable_obs_keys(adata: anndata.AnnData) -> List[str]:
    """Return categorical obs columns hidden from grouping controls."""
    return [
        f"{row['key']} ({row['n_unique']} categories)"
        for row in get_obs_category_summary(adata)
        if not row["plot_ok"]
    ]


def get_numeric_obs_keys(adata: anndata.AnnData) -> List[str]:
    """Return numeric obs columns suitable for continuous plots."""
    return [
        key for key in adata.obs_keys()
        if pd.api.types.is_numeric_dtype(adata.obs[key])
    ]


def _sample_index_values(index, n=200) -> List[str]:
    if len(index) == 0:
        return []
    sample_size = min(n, len(index))
    positions = np.linspace(0, len(index) - 1, sample_size, dtype=int)
    return [str(index[i]) for i in positions]


def _fraction_matching(values, predicate) -> float:
    if not values:
        return 0.0
    return sum(1 for value in values if predicate(value)) / len(values)


def _looks_like_cell_id(value: str) -> bool:
    value_upper = value.upper()
    value_lower = value.lower()
    if any(token in value_lower for token in CELL_ID_INDEX_KEYWORDS):
        return True
    if "-" in value and len(value) >= 10:
        return True
    if len(value) >= 12 and value_upper == value and any(c.isdigit() for c in value):
        return True
    return False


def _looks_like_gene_name(value: str) -> bool:
    value_upper = value.upper()
    if value_upper in KNOWN_GENE_NAMES:
        return True
    return any(value_upper.startswith(prefix) for prefix in GENE_NAME_PREFIXES)


def check_anndata_orientation(adata: anndata.AnnData) -> List[str]:
    """Return non-fatal warnings when an AnnData object may be transposed."""
    obs_cols = set(adata.obs_keys())
    var_cols = set(adata.var_keys())
    obs_hints = [key for key in OBS_METADATA_HINTS if key in obs_cols]
    var_obs_hints = [key for key in OBS_METADATA_HINTS if key in var_cols]
    var_feature_hints = [key for key in VAR_FEATURE_HINTS if key in var_cols]
    obs_feature_hints = [key for key in VAR_FEATURE_HINTS if key in obs_cols]

    obs_name_sample = _sample_index_values(adata.obs_names)
    var_name_sample = _sample_index_values(adata.var_names)
    obs_gene_fraction = _fraction_matching(obs_name_sample, _looks_like_gene_name)
    var_cell_fraction = _fraction_matching(var_name_sample, _looks_like_cell_id)
    obs_cell_fraction = _fraction_matching(obs_name_sample, _looks_like_cell_id)
    var_gene_fraction = _fraction_matching(var_name_sample, _looks_like_gene_name)

    evidence = []
    if var_obs_hints and not obs_hints:
        evidence.append(
            "expected cell metadata columns found in `.var` instead of `.obs`: "
            + ", ".join(var_obs_hints)
        )
    if obs_feature_hints and not var_feature_hints:
        evidence.append(
            "feature annotation columns found in `.obs` instead of `.var`: "
            + ", ".join(obs_feature_hints)
        )
    if obs_gene_fraction >= 0.25 and var_cell_fraction >= 0.25:
        evidence.append(
            "observation names look gene-like while variable names look cell/barcode-like"
        )
    elif var_cell_fraction >= 0.5 and obs_cell_fraction < 0.1:
        evidence.append("variable names look cell/barcode-like")
    elif obs_gene_fraction >= 0.5 and var_gene_fraction < 0.1:
        evidence.append("observation names look gene-like")

    if not evidence:
        return []

    return [
        "This AnnData object may be transposed (genes in `.obs` and cells in `.var`). "
        + "Evidence: "
        + "; ".join(evidence)
        + ". Downstream plots will still run on the loaded orientation, but "
        + "feature plots and gene scoring may not behave as expected."
    ]


def choose_default_option(options, preferred=None, fallback=None):
    """Choose a stable widget default from available options."""
    options = list(options or [])
    if preferred in options:
        return preferred
    if fallback in options:
        return fallback
    return options[0] if options else None


def choose_group_default(
    options,
    preferred=("cluster", "condition", "sample"),
    fallback=None,
):
    """Choose a grouping default from cardinality-filtered obs options."""
    options = list(options or [])
    for key in preferred:
        if key in options:
            return key
    if fallback in options:
        return fallback
    return options[0] if options else None


def process_matrix_layout(
    adata_all,
    n_rows: int,
    n_cols: int,
    sample_key: str = "sample",
    spatial_key: str = "spatial",
    new_obsm_key: str = "X_dataset",
    tile_spacing: float = 100.0,
    flipy: bool = False,
    sample_order_mode: str = "original",
    condition_key: str = "condition",
):
    """Add an obsm with spatial coordinates offset into a sample grid."""
    if sample_order_mode == "original":
        samples = list(pd.unique(adata_all.obs[sample_key]))
    elif sample_order_mode == "sample":
        samples = sorted(adata_all.obs[sample_key].astype(str).unique().tolist())
    elif sample_order_mode == "condition":
        obs_tmp = adata_all.obs[[sample_key, condition_key]].copy()
        cond_per_sample = (
            obs_tmp
            .assign(_i=np.arange(len(obs_tmp)))
            .sort_values("_i")
            .groupby(sample_key, sort=False)[condition_key]
            .first()
        )
        samples = (
            cond_per_sample.reset_index()
            .sort_values([condition_key, sample_key], kind="stable")
            [sample_key]
            .astype(str)
            .tolist()
        )
    else:
        raise ValueError("sample_order_mode must be one of {'original','sample','condition'}")

    n_samples = len(samples)
    total_positions = n_rows * n_cols
    if n_samples > total_positions:
        raise ValueError(
            f"Not enough grid positions ({n_rows}x{n_cols}={total_positions}) for {n_samples} samples"
        )

    X_new = np.empty_like(adata_all.obsm[spatial_key], dtype=float)
    grid_bounds = {}
    sample_positions = {}

    for idx, sample_name in enumerate(samples):
        row = idx // n_cols
        col = idx % n_cols
        sample_positions[sample_name] = (row, col)

        mask = (adata_all.obs[sample_key].astype(str) == str(sample_name))
        xspa = adata_all.obsm[spatial_key][mask]

        l_max = xspa.max(axis=0)
        l_min = xspa.min(axis=0)
        grid_bounds[(row, col)] = {
            "width": float(l_max[0] - l_min[0]),
            "height": float(l_max[1] - l_min[1]),
            "min_x": float(l_min[0]),
            "min_y": float(l_min[1]),
            "max_x": float(l_max[0]),
            "max_y": float(l_max[1]),
        }

    row_heights = [
        max(
            (
                grid_bounds[(r, c)]["height"]
                for c in range(n_cols)
                if (r, c) in grid_bounds
            ),
            default=0.0,
        )
        for r in range(n_rows)
    ]
    col_widths = [
        max(
            (
                grid_bounds[(r, c)]["width"]
                for r in range(n_rows)
                if (r, c) in grid_bounds
            ),
            default=0.0,
        )
        for c in range(n_cols)
    ]

    row_y_offsets = [0.0]
    for i in range(n_rows - 1):
        row_y_offsets.append(row_y_offsets[-1] - row_heights[i] - tile_spacing)

    col_x_offsets = [0.0]
    for i in range(n_cols - 1):
        col_x_offsets.append(col_x_offsets[-1] + col_widths[i] + tile_spacing)

    for sample_name in samples:
        row, col = sample_positions[sample_name]
        mask = (adata_all.obs[sample_key].astype(str) == str(sample_name))
        xspa = adata_all.obsm[spatial_key][mask].copy().astype(float)

        bounds = grid_bounds[(row, col)]
        target_x = col_x_offsets[col]
        target_y = row_y_offsets[row]

        if flipy:
            center_y = (bounds["min_y"] + bounds["max_y"]) / 2.0
            xspa[:, 1] = 2.0 * center_y - xspa[:, 1]

        xspa[:, 0] += target_x - bounds["min_x"]
        xspa[:, 1] += target_y - (bounds["max_y"] if not flipy else bounds["min_y"])
        X_new[mask] = xspa

    adata_all.obsm[new_obsm_key] = X_new


def rename_obs_keys(adata: anndata.AnnData) -> anndata.AnnData:
    """Add preferred obs columns by copying common alternate names."""
    key_map = {
        "Sample": "sample",
        "nFrags": "n_fragment",
        "Condition": "condition",
        "Clusters": "cluster",
    }

    keys = adata.obs_keys()
    for src, dest in key_map.items():
        if src in keys and dest not in keys:
            adata.obs[dest] = adata.obs[src]

    return adata


def reorder_obs_columns(adata, first_col="cluster"):
    """Move a column to the first obs position when present."""
    if first_col not in adata.obs.columns:
        return
    new_order = [first_col] + [c for c in adata.obs.columns if c != first_col]
    adata.obs = adata.obs[new_order]


def sort_group_categories(values):
    """Sort group labels numerically when possible, otherwise alphabetically."""
    num_vals = []
    all_numeric = True
    for value in values:
        try:
            num_vals.append(float(value))
        except (ValueError, TypeError):
            all_numeric = False
            break

    if all_numeric:
        sorted_pairs = sorted(zip(values, num_vals), key=lambda x: x[1])
        return [value for value, _ in sorted_pairs]

    return sorted(map(str, values))


def cluster_marker_to_dataframe(raw_value, key):
    """Convert an AnnData `.uns` value into a DataFrame."""
    if isinstance(raw_value, pd.DataFrame):
        return raw_value.copy()
    try:
        return pd.DataFrame(raw_value)
    except Exception as err:
        raise ValueError(f"Could not convert `adata.uns['{key}']` to a table: {err}")


def parse_cluster_marker_int(value, default, label, minimum=1):
    """Parse an integer plotting control value."""
    try:
        parsed = int(str(value).strip())
    except Exception:
        raise ValueError(f"{label} must be an integer.")
    if parsed < minimum:
        raise ValueError(f"{label} must be at least {minimum}.")
    return parsed


def parse_cluster_marker_float(value, default, label, minimum=None):
    """Parse a floating point plotting control value."""
    raw = str(value).strip()
    if raw == "":
        return default
    try:
        parsed = float(raw)
    except Exception:
        raise ValueError(f"{label} must be a number.")
    if minimum is not None and parsed < minimum:
        raise ValueError(f"{label} must be at least {minimum}.")
    return parsed


def cluster_marker_sort_key(cluster):
    """Sort cluster labels numerically when possible, otherwise alphabetically."""
    cluster_str = str(cluster)
    return (0, int(cluster_str)) if cluster_str.isdigit() else (1, cluster_str)


def cluster_marker_zscore_heatmap(values_df):
    """Z-score heatmap columns and clip to the marker heatmap display range."""
    return values_df.apply(
        lambda col: (col - col.mean()) / (col.std() + 1e-9),
        axis=0,
    ).clip(-3, 3)


def get_cached_cluster_marker_heatmap(
    adata,
    heatmap_key="cluster_marker_heatmap",
):
    """Read a precomputed cluster marker heatmap from AnnData if present."""
    if heatmap_key not in adata.uns:
        return None
    heatmap_df = cluster_marker_to_dataframe(adata.uns[heatmap_key], heatmap_key)
    heatmap_df = heatmap_df.apply(pd.to_numeric, errors="coerce")
    heatmap_df = heatmap_df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    heatmap_df.index = heatmap_df.index.map(str)
    heatmap_df.columns = heatmap_df.columns.map(str)
    if heatmap_df.empty:
        return None
    return heatmap_df


def compute_cluster_marker_heatmap_from_degs(
    adata,
    deg_df,
    top_n,
    pval_col,
    pval_cutoff,
    log2fc_cutoff,
    order_mode,
    user_order,
    deg_key="cluster_marker_degs",
    groupby="cluster",
    layer="log1p",
):
    """Build the dynamic marker heatmap from stored cluster DEG results."""
    required_cols = {groupby, "names", "logfoldchanges", pval_col}
    missing = sorted(required_cols - set(deg_df.columns))
    if missing:
        raise ValueError(f"`adata.uns['{deg_key}']` is missing columns: {missing}.")
    if layer not in adata.layers:
        raise ValueError(f"This plot requires `adata.layers['{layer}']`.")
    if groupby not in adata.obs:
        raise ValueError(f"This plot requires `adata.obs['{groupby}']`.")

    deg_df = deg_df.copy()
    deg_df[groupby] = deg_df[groupby].astype(str)
    deg_df["names"] = deg_df["names"].astype(str)
    deg_df["logfoldchanges"] = pd.to_numeric(
        deg_df["logfoldchanges"],
        errors="coerce",
    )
    deg_df[pval_col] = pd.to_numeric(deg_df[pval_col], errors="coerce")
    deg_df = deg_df.dropna(subset=[groupby, "names", "logfoldchanges", pval_col])
    deg_df = deg_df[
        (deg_df[pval_col] <= pval_cutoff)
        & (deg_df["logfoldchanges"] >= log2fc_cutoff)
    ]

    if deg_df.empty:
        raise ValueError("No DEGs remain after the selected filters.")

    clusters = sorted(deg_df[groupby].unique().tolist(), key=cluster_marker_sort_key)
    top_genes_per_cluster = {}
    for cluster in clusters:
        cluster_df = deg_df[deg_df[groupby] == cluster].head(top_n)
        top_genes_per_cluster[cluster] = cluster_df["names"].tolist()

    seen = set()
    all_top_genes = []
    for cluster in clusters:
        for gene in top_genes_per_cluster[cluster]:
            if gene not in seen:
                all_top_genes.append(gene)
                seen.add(gene)

    if len(all_top_genes) == 0:
        raise ValueError("No marker genes are available for the selected filters.")

    gene_idx = adata.var_names.get_indexer(all_top_genes)
    valid = gene_idx >= 0
    genes = [gene for gene, keep in zip(all_top_genes, valid) if keep]
    gene_idx = gene_idx[valid]
    if len(genes) == 0:
        raise ValueError("None of the selected marker genes are present in `adata.var_names`.")

    obs_clusters = adata.obs[groupby].astype(str)
    X = adata.layers[layer]
    mean_expr = pd.DataFrame(index=clusters, columns=genes, dtype=float)
    for cluster in clusters:
        mask = obs_clusters == cluster
        if int(mask.sum()) == 0:
            mean_expr.loc[cluster] = np.nan
            continue
        sub = X[mask.to_numpy(), :][:, gene_idx]
        if sp.issparse(sub):
            sub = sub.toarray()
        mean_expr.loc[cluster] = np.asarray(sub).mean(axis=0)

    mean_expr = mean_expr.dropna(axis=0, how="all").dropna(axis=1, how="all")
    if mean_expr.empty:
        raise ValueError("Unable to compute mean expression for the selected genes.")

    clusters = mean_expr.index.tolist()
    if order_mode == "DEG similarity" and len(clusters) > 1:
        values = mean_expr.to_numpy(dtype=float)
        scaled = (values - values.mean(axis=0)) / (values.std(axis=0) + 1e-9)
        cluster_order = [
            clusters[i]
            for i in leaves_list(linkage(pdist(scaled), method="ward"))
        ]
    elif order_mode == "user-selected":
        requested = [
            item.strip()
            for item in str(user_order).replace(";", ",").split(",")
            if item.strip()
        ]
        requested = [cluster for cluster in requested if cluster in clusters]
        cluster_order = requested + [
            cluster for cluster in clusters if cluster not in requested
        ]
    else:
        cluster_order = sorted(clusters, key=cluster_marker_sort_key)

    valid_gene_set = set(mean_expr.columns)
    seen_ordered = set()
    ordered_genes = []
    for cluster in cluster_order:
        for gene in top_genes_per_cluster.get(cluster, []):
            if gene in valid_gene_set and gene not in seen_ordered:
                ordered_genes.append(gene)
                seen_ordered.add(gene)

    if len(ordered_genes) == 0:
        raise ValueError("No ordered marker genes are available for the heatmap.")

    heatmap_values = mean_expr.loc[cluster_order, ordered_genes].astype(float)
    return cluster_marker_zscore_heatmap(heatmap_values)
