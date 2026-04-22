import anndata
import math
import numpy as np
import pandas as pd
import plotly.express as px
import scanpy as sc

from anndata import AnnData
from pathlib import Path
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


async def get_notebook_palettes():
    try:
        palette_data = await palettes.get()
    except Exception as e:
        print(f"Unable to load notebook palettes: {e}")
        return {"categorical": [], "continuous": []}

    if not isinstance(palette_data, dict):
        return {"categorical": [], "continuous": []}

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
    """Return obs columns suitable for grouping/filtering widgets."""
    return [
        key for key in adata.obs_keys()
        if (
            pd.api.types.is_object_dtype(adata.obs[key])
            or pd.api.types.is_categorical_dtype(adata.obs[key])
            or pd.api.types.is_bool_dtype(adata.obs[key])
        )
    ]


def get_numeric_obs_keys(adata: anndata.AnnData) -> List[str]:
    """Return numeric obs columns suitable for continuous plots."""
    return [
        key for key in adata.obs_keys()
        if pd.api.types.is_numeric_dtype(adata.obs[key])
    ]


def choose_default_option(options, preferred=None, fallback=None):
    """Choose a stable widget default from available options."""
    options = list(options or [])
    if preferred in options:
        return preferred
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
