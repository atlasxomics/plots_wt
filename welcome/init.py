import anndata
import json
import math
import matplotlib.pyplot as plt
import matplotlib.path as mplpath
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import scanpy as sc

from anndata import AnnData
from functools import lru_cache
from pathlib import Path
from plotly.subplots import make_subplots
from typing import Any, Dict, List, Optional

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
SEQLOGO_JSON_FILENAMES = {
    "hg38": "seqlogo_hg38.json",
    "mm10": "seqlogo_mm10.json",
    "rn6": "seqlogo_rn6.json",
}
SEQLOGO_BASES = ["A", "C", "G", "T"]
SEQLOGO_BASE_COLORS = {
    "A": "#2E8B57",
    "C": "#3B6FB6",
    "G": "#E39C23",
    "T": "#C53A2C",
}


# Globals ------------------------------------------------------------------

if "new_data_signal" not in globals():
    new_data_signal = Signal(False)
if "choose_group_signal" not in globals():
    choose_group_signal = Signal(False)
if "groupselect_signal" not in globals():
    groupselect_signal = Signal(False)
if "barcodes_signal" not in globals():
    barcodes_signal = Signal(False)
if "wf_ready_signal" not in globals():
    wf_ready_signal = Signal(False)
if "wf_exe_signal" not in globals():
    wf_exe_signal = Signal(False)
if "wf_results_signal" not in globals():
    wf_results_signal = Signal(False)
if "wf_bigwigs_signal" not in globals():
    wf_bigwigs_signal = Signal(False)
if "h5_viewer_signal" not in globals():
    h5_viewer_signal = Signal(False)
if "compare_signal" not in globals():
    compare_signal = Signal(False)
if "heatmap_signal" not in globals():
    heatmap_signal = Signal(False)
if "tracks_signal" not in globals():
    tracks_signal = Signal(False)
if "choose_subset_signal" not in globals():
    choose_subset_signal = Signal(False)
if "gene_score_done_signal" not in globals():
    gene_score_done_signal = Signal(False)
if "refresh_h5_signal" not in globals():
    refresh_h5_signal = Signal(False)

obsm_keys = ("X_umap", "spatial")
na_keys = ['barcode', 'on_off', 'row', 'col', 'xcor', 'ycor', 'score']
all_colors = (
    'Paired', 'Paired_r',
    'Set1', 'Set1_r',
    'Set2', 'Set2_r',
    'tab10', 'tab10_r',
    'tab20', 'tab20_r',
    'tab20b', 'tab20b_r',
    'tab20c', 'tab20c_r',
    'deep', 'muted', 'pastel', 'bright', 'dark', 'colorblind',
    "Alphabet", "Dark24", "Light24"
)

adata = None
adata_g = None
adata_m = None

# Functions ----------------------------------------------------------------
def adjust_pvals(
    df: pd.DataFrame,
    pval_col: str = "pvals_adj",
    threshold: float = 0.05,
    display_pval: bool = True
) -> pd.DataFrame:
    """
    Prepare adjusted p‐values for plotting or filtering.

    If `display_pval=True`, zeros and NaNs in `pval_col` are replaced by
    a small positive number so that –log10(p) or color scales don’t blow up.
    If `display_pval=False`, rows where `pval_col` is zero or NaN are dropped.

    Parameters
    ----------
    df
        Input DataFrame (will be copied).
    pval_col
        Name of the adjusted‐p‐value column.
    threshold
        Any nonzero pval above this is treated as too big; we’ll use eps instead.
    display_pval
        If True, replace 0/NaN with eps or fraction of the minimum nonzero pval;
        if False, drop 0/NaN rows entirely.

    Returns
    -------
    A new DataFrame with `pval_col` massaged as above.
    """
    df = df.copy()
    if pval_col not in df.columns:
        raise ValueError(f"Column '{pval_col}' not found in DataFrame.")

    if display_pval:
        # fill NaNs with 0
        df[pval_col] = df[pval_col].fillna(0)

        # find smallest nonzero
        nonzero = df.loc[df[pval_col] > 0, pval_col]
        min_nonzero = nonzero.min() if not nonzero.empty else None

        # choose replacement
        if min_nonzero is None or np.isnan(min_nonzero) or min_nonzero > threshold:
            replacement = np.finfo(float).eps
        else:
            replacement = min_nonzero / 10

        # replace zeros (and any negatives, just in case)
        df.loc[df[pval_col] <= 0, pval_col] = replacement

    else:
        # drop any rows where pval is zero or NaN
        df = df[df[pval_col].notna() & (df[pval_col] > 0)].copy()

    return df


def convert_feature_expression(anndata, feature_name):
    """Create a new column in .obs with feature value from .X.
    """
    try:
        feature_matrix = anndata[:, feature_name].X
        if hasattr(feature_matrix, "toarray"):
            feature_matrix = feature_matrix.toarray()
        feature_value = np.asarray(feature_matrix).ravel()
        anndata.obs[feature_name] = list(feature_value)
    except KeyError as e:
        print(f"Error {e}: No feature {feature_name} found in .var_names")


def create_proportion_dataframe(
    adata, group_by, stack_by, return_type="proportion"
):
    """Create a DataFrame for a proportion or raw count plot (stacked bar
    graph) from an AnnData object.
    """

    # Create a cross-tabulation (contingency table) to get the counts
    count_df = pd.crosstab(adata.obs[group_by], adata.obs[stack_by])

    if return_type == 'proportion':
        # Calculate proportions for each group
        result_df = count_df.div(count_df.sum(axis=1), axis=0)
    elif return_type == 'counts':
        # Return raw counts
        result_df = count_df
    else:
        raise ValueError("Invalid return_type. Use 'proportion' or 'counts'.")

    # Reshape the DataFrame to long format for easier plotting (stacked format)
    long_df = result_df.reset_index().melt(
      id_vars=group_by,
      value_name="value",
      var_name=stack_by
    )
    long_df.columns = ["group_by", "stack_by", "value"]

    return long_df


def create_violin_data(adata, group_by, plot_data, data_type="obs"):
    """Create a DataFrame from an AnnData object to be used for violin plots;
    returns pandas DataFrame with columns: 'group', 'value', and 'type'
    (either 'obs' or feature values from '.X').
    """

    # Check if the data to plot is from .obs
    if data_type == "obs":

        # Extract the data from .obs
        values = adata.obs[plot_data]
        df = pd.DataFrame({
            'group': adata.obs[group_by],
            'value': values
        })

    # Check if the data to plot is a feature from .X
    elif data_type in ["gene", "motif", "feature"]:

        # Extract feature values.
        feature_matrix = adata[:, plot_data].X
        if hasattr(feature_matrix, "toarray"):
            feature_matrix = feature_matrix.toarray()
        values = np.asarray(feature_matrix).ravel()

        df = pd.DataFrame({
            "group": adata.obs[group_by],
            "value": values
        })

    else:
        raise ValueError("data_type must be either 'obs', 'gene', 'motif', or 'feature'.")

    return df


def custom_plotly(
    snap_fig, color_scheme='bright', width=400, height=500, hide_axes=True
):

    orig_data = snap_fig.data
    unique_groups = sorted([d.name for d in orig_data if d.name is not None])
    n_groups = len(unique_groups)

    # Generate new colors
    colors = generate_color_palette(n_groups, color_scheme)

    group_color_map = {unique_groups[i]: colors[i] for i in range(n_groups)}

    # Create new figure
    new_fig = go.Figure()

    new_fig.update_layout(
        snap_fig.layout
    )
    # Add traces with new colors
    for i, trace in enumerate(orig_data):
        if trace.name is not None:  # Skip traces without names
            new_trace = go.Scatter(
                x=trace.x,
                y=trace.y,
                mode=trace.mode,
                name=trace.name,
                marker=dict(
                    size=trace.marker.size,
                    color=group_color_map[trace.name],
                    opacity=trace.marker.opacity if trace.marker.opacity else 0.7
                ),
                showlegend=trace.showlegend
            )
            new_fig.add_trace(new_trace)

    # Update only width and height
    new_fig.update_layout(
        width=width,
        height=height
    )

    if hide_axes:
        new_fig.update_layout(
            xaxis=dict(showticklabels=False, showline=False, zeroline=False, showgrid=False),
            yaxis=dict(showticklabels=False, showline=False, zeroline=False, showgrid=False)
        )

    return new_fig


def drop_obs_column(adata_objects, col_to_drop="orig.ident"):
    """Remove specified column from obs DataFrame of multiple AnnData objects."""
    for adata_obj in adata_objects:
        if col_to_drop in adata_obj.obs.columns:
            adata_obj.obs = adata_obj.obs.drop(columns=[col_to_drop])


def filter_adata_by_groups(adata, group, group_a, group_b="All"):
    """Filter adata to two values in obs."""
    assert group_a != group_b, "Groups must be different."
    return adata[adata.obs[group].isin([group_a, group_b])]


def filter_anndata(
    adata: AnnData, group: str, subgroup: List[str], mem=False
) -> AnnData:
    if mem:
        return adata[adata.obs[group] == subgroup].to_memory()
    return adata[adata.obs[group] == subgroup]


def generate_color_palette(length, scheme="bright"):
    """
    Generate a list of hex color codes  with the specified number of colors.
    
    length : int
        Number of colors in the palette
    scheme : str, default="bright"
        Name of the color scheme to use. Can be one of:
        - Matplotlib palettes: 'Paired', 'Set1', 'Set2', 'tab10', 'tab20', etc.
        - Seaborn palettes: 'deep', 'muted', 'pastel', 'bright', 'dark', 'colorblind'
        - Plotly palettes: 'Alphabet', 'Dark24', 'Light24'
        
    """

    
    # Helper function to convert RGB to hex
    def rgb_to_hex(rgb):
        return "#{:02x}{:02x}{:02x}".format(
            int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255)
        )
    
    matplotlib_palettes = (
        "Paired", "Paired_r",
        "Set1", "Set1_r",
        "Set2", "Set2_r",
        "tab10", "tab10_r",
        "tab20", "tab20_r",
        "tab20b", "tab20b_r",
        "tab20c", "tab20c_r",
    )
    
    seaborn_palettes = (
        "deep", "muted", "pastel", "bright", "dark", "colorblind"
    )
    
    plotly_palettes = {
        "Alphabet": [
            "#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", 
            "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F", "#C4451C", "#DEA0FD", 
            "#FE00FA", "#325A9B", "#FEAF16", "#F8A19F", "#90AD1C", "#F6222E", 
            "#1CFFCE", "#2ED9FF", "#B10DA1", "#C075A6", "#FC1CBF", "#B00068", 
            "#FBE426", "#FA0087"
        ],
        "Dark24": [
            "#2E91E5", "#E15F99", "#1CA71C", "#FB0D0D", "#DA16FF", "#222A2A", 
            "#B68100", "#750D86", "#EB663B", "#511CFB", "#00A08B", "#FB00D1", 
            "#FC0080", "#B2828D", "#6C7C32", "#778AAE", "#862A16", "#A777F1", 
            "#620042", "#1616A7", "#DA60CA", "#6C4516", "#0D2A63", "#AF0038"
        ],
        "Light24": [
            "#FD3216", "#00FE35", "#6A76FC", "#FED4C4", "#FE00CE", "#0DF9FF", 
            "#F6F926", "#FF9616", "#479B55", "#EEA6FB", "#DC587D", "#D626FF", 
            "#6E899C", "#00B5F7", "#B68E00", "#C9FBE5", "#FF0092", "#22FFA7", 
            "#E3EE9E", "#86CE00", "#BC7196", "#7E7DCD", "#FC6955", "#E48F72"
        ]
    }

    if scheme in matplotlib_palettes:
        if scheme == 'rainbow':
            cm = plt.cm.rainbow
        else:
            cm = plt.get_cmap(scheme)
        colors = [rgb_to_hex(cm(i)[:3]) for i in np.linspace(0, 1, length)]
    elif scheme in seaborn_palettes:
        colors = sns.color_palette(scheme, length)
        colors = [rgb_to_hex(c) for c in colors]
    elif scheme in plotly_palettes:
        # Get Plotly palette and cycle if needed
        plotly_colors = plotly_palettes[scheme]
        colors = [plotly_colors[i % len(plotly_colors)] for i in range(length)]
    else:
        cm = plt.get_cmap('viridis')
        colors = [rgb_to_hex(cm(i)[:3]) for i in np.linspace(0, 1, length)]

    return colors


async def get_notebook_palettes():
    try:
        palette_data = await palettes.get()
    except Exception as e:
        print(f"Unable to load notebook palettes: {e}")
        return {"categorical": [], "continuous": []}

    if not isinstance(palette_data, dict):
        return {"categorical": [], "continuous": []}

    cleaned_palettes = {}
    for kind in ("categorical", "continuous"):
        cleaned_palettes[kind] = [
            palette for palette in palette_data.get(kind, [])
            if palette.get("display_name") and palette.get("colors")
        ]

    return cleaned_palettes


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


def normalize_plotly_colorscale(colorscale):
    if isinstance(colorscale, list) and len(colorscale) > 0:
        if isinstance(colorscale[0], str):
            if len(colorscale) == 1:
                return [[0, colorscale[0]], [1, colorscale[0]]]
            return [
                [i / (len(colorscale) - 1), color]
                for i, color in enumerate(colorscale)
            ]

    return colorscale


def build_discrete_color_map(categories, colors):
    if not colors:
        colors = DEFAULT_H5_CATEGORICAL_PALETTE

    return {category: colors[i % len(colors)] for i, category in enumerate(categories)}


def get_current_igv_genome(default: str = "hg38") -> str:
    if "coverages_genome" in globals() and coverages_genome.value is not None:
        return coverages_genome.value
    return default


def resolve_seqlogo_json_path(genome: str) -> Optional[Path]:
    filename = SEQLOGO_JSON_FILENAMES.get(genome)
    if filename is None:
        return None

    path = Path.cwd() / filename
    if path.exists():
        return path

    return None


@lru_cache(maxsize=8)
def load_seqlogos(path: str) -> Dict[str, pd.DataFrame]:
    with open(path) as f:
        raw = json.load(f)

    seqlogos: Dict[str, pd.DataFrame] = {}
    for motif, rows in raw.items():
        normalized_rows = [
            row[0] if isinstance(row, list) else row
            for row in rows
        ]
        df = pd.DataFrame(normalized_rows)
        missing_bases = [base for base in SEQLOGO_BASES if base not in df.columns]
        if missing_bases:
            raise ValueError(
                f"Motif '{motif}' is missing base columns: {', '.join(missing_bases)}"
            )

        df = df[SEQLOGO_BASES].astype(float)
        seqlogos[motif] = df

    return seqlogos


def get_seqlogos_for_genome(genome: str) -> Dict[str, pd.DataFrame]:
    seqlogo_path = resolve_seqlogo_json_path(genome)
    if seqlogo_path is None:
        raise FileNotFoundError(
            f"No seqlogo JSON found for genome '{genome}'. "
            f"Expected `{SEQLOGO_JSON_FILENAMES.get(genome, genome)}` "
            "in the notebook working directory."
        )

    return load_seqlogos(str(seqlogo_path))


def normalize_motif_logo_name(name: str) -> str:
    return name.replace("-", "_")


def resolve_motif_logo_name(
    motif_name: str,
    seqlogos: Dict[str, pd.DataFrame],
) -> Optional[str]:
    if motif_name in seqlogos:
        return motif_name

    normalized_name = normalize_motif_logo_name(motif_name)
    if normalized_name in seqlogos:
        return normalized_name

    return None


def get_available_motif_logos(
    motif_names: List[str],
    genome: str,
) -> Dict[str, pd.DataFrame]:
    seqlogos = get_seqlogos_for_genome(genome)
    matched_logos: Dict[str, pd.DataFrame] = {}
    for motif in motif_names:
        resolved_name = resolve_motif_logo_name(motif, seqlogos)
        if resolved_name is not None:
            matched_logos[motif] = seqlogos[resolved_name]

    return matched_logos


def _draw_logo_letter(ax, letter: str, x: float, y: float, height: float, width: float) -> None:
    if height <= 0:
        return

    text_path = TextPath((0, 0), letter, size=1)
    bbox = text_path.get_extents()
    if bbox.width == 0 or bbox.height == 0:
        return

    transform = (
        Affine2D()
        .translate(-bbox.x0, -bbox.y0)
        .scale(width / bbox.width, height / bbox.height)
        .translate(x, y)
    )
    patch = PathPatch(
        text_path,
        transform=transform + ax.transData,
        facecolor=SEQLOGO_BASE_COLORS.get(letter, "black"),
        edgecolor="none",
    )
    ax.add_patch(patch)


def plot_motif_logo(
    prob_df: pd.DataFrame,
    title: str = "",
    information_content: bool = True,
) -> plt.Figure:
    df = prob_df.copy()[SEQLOGO_BASES].astype(float)
    df = df.div(df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    if information_content:
        entropy = -(df * np.log2(df.replace(0, np.nan))).sum(axis=1).fillna(0.0)
        total_height = 2.0 - entropy
    else:
        total_height = pd.Series(1.0, index=df.index)

    letter_heights = df.mul(total_height, axis=0)

    fig_width = max(6, min(16, 0.55 * len(df) + 1.5))
    fig, ax = plt.subplots(figsize=(fig_width, 3.0))

    for i, (_, row) in enumerate(letter_heights.iterrows()):
        y_offset = 0.0
        for base, height in sorted(row.items(), key=lambda item: item[1]):
            _draw_logo_letter(ax, base, i + 0.1, y_offset, float(height), 0.8)
            y_offset += float(height)

    if information_content:
        ymax = float(letter_heights.sum(axis=1).max()) + 0.1
    else:
        ymax = max(1.05, float(letter_heights.sum(axis=1).max()) + 0.02)
    ax.set_xlim(0, len(df))
    ax.set_ylim(0, ymax)
    ax.set_xticks(np.arange(len(df)) + 0.5)
    ax.set_xticklabels(np.arange(1, len(df) + 1))
    ax.set_xlabel("Position")
    ax.set_ylabel("Bits" if information_content else "Probability")
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    fig.tight_layout()

    return fig


def get_barcodes_by_lasso(
    adata: anndata.AnnData,
    obsm_key: str,
    lasso_points: list[tuple[int, int]],
) -> List:
    """Function to retrieve cell barcodes based on lasso-select"""

    embedding = adata.obsm[obsm_key][:, :2]  # type: ignore

    if len(lasso_points) < 3:
        return []

    polygon = mplpath.Path(lasso_points)
    mask = polygon.contains_points(embedding)
    return adata.obs_names[mask]


def get_feature_heatmap(df, features, rank_by="scores"):
    """Returns a heatmap DataFrame by selecting the selected features per group
    based on a ranking column.
    """
    groups = df["group"].unique()

    # Extract gene scores from each group
    feats_dfs = []
    for group in groups:

        # Filter the main DataFrame to only include the selected top genes
        feats_df = df[df["names"].isin(features)]

        feats_dfs.append(feats_df)

    # Combine all top n genes into a single DataFrame
    combined_df = pd.concat(feats_dfs, ignore_index=True)

    # Prepare the DataFrame for heatmap plotting
    pivot_df = combined_df[["group", "names", rank_by]]
    pivot_df["names"] = pd.Categorical(
      pivot_df["names"], categories=pivot_df["names"].unique(), ordered=True
    )

    heatmap_df = pivot_df.pivot_table(
      index='group', columns='names', values=rank_by, aggfunc='mean'
    )

    return heatmap_df


def get_groups(adata: anndata.AnnData) -> List[str]:
    """Set 'groups' list for differential analysis."""

    groups = []
    for group in ["cluster", "sample", "condition"]:
        length = 1
        try:
            length = len(adata.obs[group].unique())
        except KeyError as e:
            print(f"{e}")

        if length > 1:
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


def get_top_n_heatmap(df, rank_by="scores", n_top=5):
    """Returns a heatmap DataFrame by selecting the top n features per group
    based on a ranking column.
    """

    groups = df["group"].unique()

    # Loop through each group and select the top n genes
    topn_gene_dfs = []
    for group in groups:

        group_df = df[df["group"] == group]

        if rank_by in ["logfoldchanges", "scores"]:
            topn = group_df.nlargest(n_top, columns=rank_by)["names"]
        elif rank_by in ["pvals", "pvals_adj"]:
            topn = group_df.nsmallest(n_top, columns=rank_by)["names"]

        # Filter the main DataFrame to only include the selected top genes
        gene_df = df[df["names"].isin(topn)]

        topn_gene_dfs.append(gene_df)

    # Combine all top n genes into a single DataFrame
    combined_df = pd.concat(topn_gene_dfs, ignore_index=True)

    # Prepare the DataFrame for heatmap plotting
    pivot_df = combined_df[["group", "names", rank_by]]
    pivot_df["names"] = pd.Categorical(
      pivot_df["names"], categories=pivot_df["names"].unique(), ordered=True
    )

    heatmap_df = pivot_df.pivot_table(
      index='group', columns='names', values=rank_by, aggfunc='mean'
    )

    return heatmap_df


def coerce_uns_to_df(obj):
    """Best-effort conversion of .uns payloads into DataFrame."""
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    if isinstance(obj, dict):
        try:
            return pd.DataFrame(obj)
        except Exception:
            return None
    if isinstance(obj, np.ndarray) and obj.dtype.names is not None:
        try:
            return pd.DataFrame.from_records(obj)
        except Exception:
            return None
    try:
        return pd.DataFrame(obj)
    except Exception:
        return None


def resolve_heatmap_stats_table(adata_hm, hm_feats, hm_group):
    """Resolve the stats table used by heatmap plotting."""
    if hm_feats == "gene":
        feature_label = "gene"
        key_map = {
            "cluster": ["ranked_genes_per_cluster", "marker_genes_per_cluster"],
            "sample": ["ranked_genes_per_sample", "marker_genes_per_sample"],
            "condition": [
                "ranked_genes_per_condition",
                "marker_genes_per_condition",
                "ranked_genes_per_conditions1",
                "marker_genes_per_conditions1",
                "ranked_genes_per_condition_1",
                "marker_genes_per_condition_1",
            ],
        }
    elif hm_feats == "motif":
        feature_label = "motif"
        key_map = {
            "cluster": ["enrichedMotifs_cluster"],
            "sample": ["enrichedMotifs_sample"],
            "condition": [
                "enrichedMotifs_condition",
                "enrichedMotifs_conditions1",
                "enrichedMotifs_condition_1",
            ],
        }
    else:
        raise ValueError(f"Unsupported heatmap feature type: {hm_feats}")

    key_candidates = key_map.get(hm_group, [])
    for key in key_candidates:
        if key in adata_hm.uns:
            stats_df = coerce_uns_to_df(adata_hm.uns[key])
            if stats_df is not None and not stats_df.empty:
                return feature_label, key, stats_df, key_candidates

    raise ValueError(
        f"No stats table found in `.uns` for {feature_label} and group "
        f"'{hm_group}'. Tried keys: {', '.join(key_candidates)}"
    )


def get_heatmap_stats_columns(stats_df, stats_key, preferred_sig_metric="FDR"):
    """Detect required columns in the stats table."""
    if "group_name" in stats_df.columns:
        group_col = "group_name"
    elif "group" in stats_df.columns:
        group_col = "group"
    else:
        raise ValueError(
            f"Could not identify group column in stats table '{stats_key}'."
        )

    feature_col = None
    for c in ["name", "names", "feature", "features", "motif", "gene"]:
        if c in stats_df.columns:
            feature_col = c
            break
    if feature_col is None:
        raise ValueError(
            f"Could not identify feature column in stats table '{stats_key}'."
        )

    if "Log2FC" not in stats_df.columns:
        raise ValueError(
            f"Log2FC is required for heatmap value/ranking but is missing in "
            f"stats table '{stats_key}'."
        )

    sig_aliases = {
        "FDR": ["FDR", "p_val_adj", "pvals_adj", "pval_adj", "padj"],
        "Pval": ["Pval", "p_val", "pvals", "pval", "PValue", "p_value"],
    }
    sig_cols = {}
    for metric, aliases in sig_aliases.items():
        for c in aliases:
            if c in stats_df.columns:
                sig_cols[metric] = c
                break

    available_sig_metrics = tuple(
        metric for metric in ("FDR", "Pval") if metric in sig_cols
    )
    fallback_order = [preferred_sig_metric, "FDR", "Pval"]
    selected_sig_metric = next(
        (metric for metric in fallback_order if metric in sig_cols),
        None,
    )
    sig_col = sig_cols.get(selected_sig_metric)

    return group_col, feature_col, sig_col, selected_sig_metric, available_sig_metrics


def prepare_heatmap_work_df(
    stats_df,
    group_col,
    feature_col,
    value_metric="Log2FC",
    rank_metric="Log2FC",
    sig_col=None,
    sig_threshold=0.01,
):
    """Prepare and filter the long-form stats table for heatmap generation."""
    keep_cols = [group_col, feature_col, value_metric, rank_metric]
    if sig_col is not None:
        keep_cols.append(sig_col)
    keep_cols = list(dict.fromkeys(keep_cols))

    work_df = stats_df[keep_cols].copy()
    work_df[group_col] = work_df[group_col].astype(str)
    work_df[feature_col] = work_df[feature_col].astype(str)
    work_df[value_metric] = pd.to_numeric(work_df[value_metric], errors="coerce")
    work_df[rank_metric] = pd.to_numeric(work_df[rank_metric], errors="coerce")
    if sig_col is not None:
        work_df[sig_col] = pd.to_numeric(work_df[sig_col], errors="coerce")

    work_df = work_df.dropna(subset=[group_col, feature_col, value_metric])
    if sig_col is not None:
        work_df = work_df[work_df[sig_col] <= sig_threshold]

    if work_df.empty:
        raise ValueError("No rows passed the selected significance filter.")

    return work_df


def select_archr_like_heatmap_features(
    work_df,
    group_col,
    feature_col,
    rank_metric,
    top_n,
    effect_threshold,
    effect_direction,
    feature_input="",
):
    """Select heatmap features from user input or top-N directional ranking."""
    selected_features = [x.strip() for x in feature_input.split(",") if x.strip()]
    selected_features = list(dict.fromkeys(selected_features))

    if len(selected_features) > 0:
        filt_df = work_df[work_df[feature_col].isin(selected_features)]
        if filt_df.empty:
            raise ValueError(
                "None of the requested features were found after filtering."
            )
        return filt_df, selected_features

    rank_df = work_df.dropna(subset=[rank_metric]).copy()
    if rank_df.empty:
        raise ValueError(
            f"No numeric values found for ranking metric '{rank_metric}'."
        )

    if effect_direction == "positive":
        rank_df = rank_df[rank_df[rank_metric] >= effect_threshold]
    elif effect_direction == "negative":
        rank_df = rank_df[rank_df[rank_metric] <= -effect_threshold]
    else:
        rank_df = rank_df[rank_df[rank_metric].abs() >= effect_threshold]

    if rank_df.empty:
        raise ValueError(
            "No features passed current effect threshold/direction settings."
        )

    selected = []
    group_order = sort_group_categories(rank_df[group_col].unique().tolist())
    for grp in group_order:
        group_df = rank_df[rank_df[group_col] == grp]
        if group_df.empty:
            continue

        if effect_direction == "positive":
            group_feats = group_df.nlargest(top_n, rank_metric)[feature_col].tolist()
        elif effect_direction == "negative":
            group_feats = group_df.nsmallest(top_n, rank_metric)[feature_col].tolist()
        else:
            group_feats = (
                group_df.assign(_abs_rank=group_df[rank_metric].abs())
                .nlargest(top_n, "_abs_rank")[feature_col]
                .tolist()
            )

        for feat in group_feats:
            if feat not in selected:
                selected.append(feat)

    filt_df = work_df[work_df[feature_col].isin(selected)]
    return filt_df, selected


def build_archr_like_heatmap_df(
    work_df,
    hm_feats,
    group_col,
    feature_col,
    value_metric="Log2FC",
    sig_col=None,
    z_clip=2.0,
):
    """Build ArchR-like heatmap matrix and legend text from filtered stats."""
    value_df = work_df.pivot_table(
        index=feature_col,
        columns=group_col,
        values=value_metric,
        aggfunc="max",
    )
    if value_df.empty:
        raise ValueError("No values available to plot after matrix construction.")

    ordered_cols = sort_group_categories(value_df.columns.tolist())
    value_df = value_df[[c for c in ordered_cols if c in value_df.columns]]

    if hm_feats == "motif" and sig_col is not None:
        p_df = work_df.pivot_table(
            index=feature_col,
            columns=group_col,
            values=sig_col,
            aggfunc="min",
        )
        p_df = p_df[[c for c in ordered_cols if c in p_df.columns]]
        p_df = p_df.clip(lower=np.finfo(float).tiny)
        score_df = -np.log10(p_df)
        row_min = score_df.min(axis=1)
        row_max = score_df.max(axis=1)
        denom = (row_max - row_min).replace(0, np.nan)
        heatmap_df = score_df.sub(row_min, axis=0).div(denom, axis=0).fillna(0) * 100.0
        legend_title = "row-normalized -log10(adj p-value) [0-100]"
    else:
        row_mean = value_df.mean(axis=1)
        row_std = value_df.std(axis=1, ddof=0).replace(0, np.nan)
        heatmap_df = value_df.sub(row_mean, axis=0).div(row_std, axis=0).fillna(0)
        heatmap_df = heatmap_df.clip(-z_clip, z_clip)
        legend_title = f"row z-score ({value_metric}, clip ±{z_clip:g})"

    return heatmap_df, legend_title


def is_motif_wide_stats_table(stats_df):
    """Detect ArchR motif enrichment table stored in wide format."""
    if not isinstance(stats_df, pd.DataFrame):
        return False
    if "group_name" not in stats_df.columns:
        return False

    group_name_vals = set(stats_df["group_name"].astype(str).unique())
    has_mlog = ("mlog10Padj" in group_name_vals) or ("mlog10p" in group_name_vals)
    has_feature = "feature" in group_name_vals
    if not (has_mlog and has_feature):
        return False

    id_col = stats_df.columns[0]
    meta_cols = {id_col, "group", "group_name"}
    value_cols = [c for c in stats_df.columns if c not in meta_cols]
    return len(value_cols) > 0


def _infer_motif_wide_value_cols(stats_df):
    """Infer value columns for wide motif tables."""
    id_col = stats_df.columns[0]
    meta_cols = {id_col, "group", "group_name"}
    value_cols = [c for c in stats_df.columns if c not in meta_cols]

    if len(value_cols) == 0:
        raise ValueError("Could not infer value columns for motif stats table.")

    def _x_key(c):
        c_str = str(c)
        if c_str.startswith("X") and c_str[1:].isdigit():
            return (0, int(c_str[1:]))
        return (1, c_str)

    value_cols = sorted(value_cols, key=_x_key)
    return id_col, value_cols


def build_motif_metric_matrix_from_wide(
    stats_df,
    metric_name,
    group_labels=None,
):
    """Build a feature x group matrix for one motif metric from wide-format stats."""
    id_col, value_cols = _infer_motif_wide_value_cols(stats_df)

    metric_rows = stats_df[stats_df["group_name"].astype(str) == metric_name].copy()
    if metric_rows.empty:
        raise ValueError(f"Metric '{metric_name}' not found in motif stats table.")

    metric_df = metric_rows[[id_col] + value_cols].copy()
    metric_df[value_cols] = metric_df[value_cols].apply(pd.to_numeric, errors="coerce")
    metric_df[id_col] = metric_df[id_col].astype(str)
    metric_df = metric_df.set_index(id_col)

    # Prefer canonical motif names from the `feature` block when aligned.
    feature_rows = stats_df[stats_df["group_name"].astype(str) == "feature"].copy()
    if len(feature_rows) == len(metric_df):
        feature_names = []
        for _, row in feature_rows[value_cols].iterrows():
            vals = [
                str(v).strip()
                for v in row.tolist()
                if str(v).strip() not in {"", "nan", "None"}
            ]
            feature_names.append(vals[0] if len(vals) > 0 else None)
        if any(name is not None for name in feature_names):
            fallback_idx = metric_df.index.tolist()
            resolved = [
                feature_names[i] if feature_names[i] is not None else fallback_idx[i]
                for i in range(len(fallback_idx))
            ]
            metric_df.index = resolved

    metric_df = metric_df.groupby(metric_df.index).max()

    if group_labels is not None and len(group_labels) == len(metric_df.columns):
        metric_df.columns = list(group_labels)

    return metric_df


def build_motif_archr_like_heatmap_from_wide(
    stats_df,
    *,
    sig_threshold=0.01,
    top_n=25,
    feature_input="",
    group_labels=None,
):
    """Build ArchR-like motif heatmap from wide motif enrichment stats."""
    metric_for_sig = "mlog10Padj"
    if metric_for_sig not in stats_df["group_name"].astype(str).unique():
        metric_for_sig = "mlog10p"

    mlog10_df = build_motif_metric_matrix_from_wide(
        stats_df,
        metric_name=metric_for_sig,
        group_labels=group_labels,
    )

    cutoff = -np.log10(max(sig_threshold, np.finfo(float).tiny))

    selected_features = [x.strip() for x in feature_input.split(",") if x.strip()]
    selected_features = list(dict.fromkeys(selected_features))

    if len(selected_features) > 0:
        sel = [f for f in selected_features if f in mlog10_df.index]
        if len(sel) == 0:
            raise ValueError("None of the requested motif features were found.")
    else:
        sel = []
        for col in mlog10_df.columns:
            scores = mlog10_df[col].dropna()
            scores = scores[scores >= cutoff]
            if scores.empty:
                continue
            top_feats = scores.sort_values(ascending=False).head(top_n).index.tolist()
            for feat in top_feats:
                if feat not in sel:
                    sel.append(feat)

        if len(sel) == 0:
            raise ValueError(
                "No motifs passed the selected significance/cutoff settings."
            )

    selected_mlog10 = mlog10_df.loc[sel]
    row_min = selected_mlog10.min(axis=1)
    row_max = selected_mlog10.max(axis=1)
    denom = (row_max - row_min).replace(0, np.nan)
    heatmap_df = (
        selected_mlog10.sub(row_min, axis=0).div(denom, axis=0).fillna(0) * 100.0
    )
    legend_title = "row-normalized -log10(adj p-value) [0-100]"

    meta = {
        "metric": metric_for_sig,
        "mlog10_cutoff": cutoff,
        "selected_features": sel,
    }
    return heatmap_df, legend_title, meta


def make_volcano_df(
    adata,
    group,
    group_a,
    group_b,
    feature,
    threshold=0.01,
    display_gm=True
):
    """Using sc.get.rank_genes_groups_df, make dataframe for volcano plot.
    Optionally filter out genes with the "Gm" prefix.
    """

    assert group_a != group_b, "Groups must be different."
    assert group in adata.obs.columns, f"No group {group} for in AnnData."

    subgroups = adata.obs[group].unique()
    assert group_a in subgroups, f"Group A {group_a} not found in subgroups."

    key = f"{group_a}_{group_b}_genes"

    if group_b == "All":
        df = sc.get.rank_genes_groups_df(
            adata, group=group_a, key=f"{group}_{feature}"
        )
    else:
        assert group_b in subgroups, f"Group B {group_b} not found in subgroups."
        adata = filter_adata_by_groups(adata, group, group_a, group_b)
        sc.tl.rank_genes_groups(
            adata,
            groupby=group,
            method="t-test",
            key_added=key,
            use_raw=False
        )
        df = sc.get.rank_genes_groups_df(adata, group=group_a, key=key)

    if not display_gm:
        df = df[~df["names"].str.startswith("Gm")]

    return df


def plot_neighborhood_groups(
  group_adatas: Dict["str", anndata.AnnData],
  title: str,
  key: str = "cluster",
  uns_key: Optional[str] = None,
  mode: str = 'zscore',
  vmin: Optional[float] = None,
  vmax: Optional[float] = None,
  colorscale: Any = None,
):
  groups = list(group_adatas.keys())
  num_groups = len(groups)

  # Dynamically determine optimal grid layout based on number of groups
  if num_groups <= 2:
    num_cols = num_groups
    num_rows = 1
    base_width = 550 * num_cols
    base_height = 506 * num_rows
  elif num_groups <= 4:
    num_cols = 2
    num_rows = (num_groups + 1) // 2
    base_width = 550 * num_cols
    base_height = 506 * num_rows
  elif num_groups <= 9:
    num_cols = 3
    num_rows = (num_groups + 2) // 3
    base_width = 400 * num_cols
    base_height = 360 * num_rows
  else:
    num_cols = 4
    num_rows = (num_groups + 3) // 4
    base_width = 336 * num_cols
    base_height = 302 * num_rows

  # Create subplots with calculated rows and columns
  combined_fig = make_subplots(
    rows=num_rows,
    cols=num_cols,
    subplot_titles=groups,
    horizontal_spacing=0.05,  # Reduced horizontal spacing
    vertical_spacing=0.05     # Reduced vertical spacing
  )

  # Create a shared colorscale range for all subplots if not provided
  if vmin is None or vmax is None:
    all_mins = []
    all_maxs = []
    for group in groups:
      if uns_key in group_adatas[group].uns:
        data = group_adatas[group].uns[uns_key][mode]
        all_mins.append(np.nanmin(data))
        all_maxs.append(np.nanmax(data))

    if vmin is None and all_mins:
      vmin = min(all_mins)
    if vmax is None and all_maxs:
      vmax = max(all_maxs)

  if colorscale is not None:
    custom_colorscale = normalize_plotly_colorscale(colorscale)
  else:
    # Create dynamic colorscale with white at zero
    abs_max = max(abs(vmin) if vmin is not None else 0, abs(vmax) if vmax is not None else 0)

    # If all values are positive or all negative, create appropriate one-sided colorscale
    if vmin is not None and vmax is not None:
      if vmin >= 0:  # All positive values
        custom_colorscale = [
          [0, 'white'],
          [1, 'red']
        ]
      elif vmax <= 0:  # All negative values
        custom_colorscale = [
          [0, 'blue'],
          [1, 'white']
        ]
      else:  # Mixed positive and negative values
        # Calculate the midpoint (0) in the normalized scale
        midpoint = abs(vmin) / (abs(vmin) + abs(vmax))

        # Create a colorscale with white at the midpoint
        custom_colorscale = [
          [0, 'blue'],
          [midpoint, 'white'],
          [1, 'red']
        ]
    else:
      # Default to RdBu_r if bounds not determined
      custom_colorscale = "RdBu_r"

  # Loop through each sample
  data_tables = {}
  for i, group in enumerate(groups):
    row = (i // num_cols) + 1
    col = (i % num_cols) + 1

    # Get data directly rather than creating a full figure
    adata = group_adatas[group]

    # Validate input
    if uns_key not in adata.uns:
      continue

    # Get the data
    data = adata.uns[uns_key][mode]
    categories = adata.obs[key].cat.categories
    row_labels = categories.copy()
    col_labels = categories.copy()

    # Sort categories and data numerically
    numeric_like = np.array([
        float(c) if isinstance(c, str) and c.replace('.', '', 1).isdigit() else np.inf
        for c in categories
    ])
    sorted_idx = np.argsort(numeric_like)
    data = data[sorted_idx][:, sorted_idx]
    row_labels = row_labels[sorted_idx]
    col_labels = col_labels[sorted_idx]
    data_tables[group] = data

    # Create heatmap - only the first subplot will show a colorbar
    heatmap = go.Heatmap(
      z=data,
      x=row_labels,
      y=col_labels,
      colorscale=custom_colorscale,
      showscale=(i == 0),  # Only show colorbar for the first subplot
      textfont={"size": 12},  # Smaller font for dense plots
      hoverongaps=False,
      zmin=vmin,
      zmax=vmax,
      colorbar=dict(
        title=mode,
        tickformat='.2f',
        len=0.9,
        x=1.02,  # Position colorbar slightly to the right
        yanchor="middle"
      ) if i == 0 else None
    )

    # Add trace to combined figure
    combined_fig.add_trace(heatmap, row=row, col=col)

    # Update axes for each subplot - make more compact
    combined_fig.update_xaxes(
      showgrid=False,
      title=None,
      side='bottom',
      tickfont=dict(size=10),  # Smaller font size for tick labels
      row=row,
      col=col
    )
    combined_fig.update_yaxes(
      showgrid=False,
      title=None,
      autorange='reversed',
      tickfont=dict(size=10),  # Smaller font size for tick labels
      row=row,
      col=col
    )

  # For very large numbers of samples, increase the base size
  if num_groups > 16:
    base_width = max(base_width, 1000)
    base_height = max(base_height, 900)

  # Update layout once for all subplots
  combined_fig.update_layout(
    title={
      'text': title,
      'x': 0.5,  # Center the title
      'xanchor': 'center',
      'yanchor': 'top',
      'font': {'size': 18}  # Slightly larger font size for title
    },
    plot_bgcolor='rgba(0,0,0,0)',
    autosize=False,  # Explicitly set size
    width=base_width,
    height=base_height,
    margin=dict(l=80, r=80, t=100, b=40),  # Tighter margins
    showlegend=False
  )

  # Update subplot titles with smaller font
  for i in range(len(combined_fig.layout.annotations)):
    combined_fig.layout.annotations[i].font.size = 14

  combined_data = pd.concat(
      {k: pd.DataFrame(v) for k, v in data_tables.items()}, axis=1
    )

  return combined_fig, combined_data


def plot_umap_for_samples(
  adata,
  samples,
  color_by='cluster',
  pt_size=3,
  coords="spatial",
  flipY=True,
  color_scheme='bright',
  show_cluster=None,
  vmin=None,
  vmax=None
):
    # Determine if color_by is discrete or continuous
    obs_values = adata.obs[color_by]
    is_discrete = pd.api.types.is_categorical_dtype(obs_values) or obs_values.dtype.name == 'category' or obs_values.nunique() < 30

    print(f"Coloring by: {color_by} (Discrete: {is_discrete})")

    # Automatically set vmin and vmax if not provided and data is continuous
    if not is_discrete:
        if vmin is None:
            vmin = obs_values.min()
        if vmax is None:
            vmax = obs_values.max()

    if is_discrete:
        obs_groups = sorted(obs_values.unique())
        colors = generate_color_palette(len(obs_groups), color_scheme)
        group_color_map = {obs_groups[i]: colors[i] for i in range(len(obs_groups))}
    else:
        group_color_map = None

    num_rows = (len(samples) - 1) // 3 + 1
    num_cols = min(len(samples), 3)

    combined_fig = make_subplots(
        rows=num_rows,
        cols=num_cols,
        subplot_titles=samples,
        horizontal_spacing=0,
        vertical_spacing=0.05
    )

    shown_clusters = set() if is_discrete else None

    if flipY:
        flipped_y = adata.obsm['spatial'].copy()
        flipped_y[:, 1] = -flipped_y[:, 1]
        adata.obsm['spatial_flippedY'] = flipped_y
        coords = 'spatial_flippedY'

    for i, sample in enumerate(samples):
        row = i // 3 + 1
        col = i % 3 + 1

        sample_data = adata[adata.obs['sample'] == sample]

        def apply_trace_settings(trace):
            if is_discrete:
                trace['marker']['color'] = group_color_map.get(trace.name, 'black')
                if trace.name not in shown_clusters:
                    shown_clusters.add(trace.name)
                else:
                    trace.showlegend = False
            else:
                if vmin is not None:
                    trace['marker']['cmin'] = vmin
                if vmax is not None:
                    trace['marker']['cmax'] = vmax

        if show_cluster is not None:
            highlight_mask = sample_data.obs['cluster'] == show_cluster
            highlight_data = sample_data[highlight_mask]
            background_data = sample_data[~highlight_mask]

            bg_fig = snap.pl.umap(
                background_data,
                color=color_by,
                use_rep=coords,
                marker_size=pt_size,
                show=False
            )
            for trace in bg_fig.data:
                trace['marker']['color'] = 'lightgrey'
                trace['marker']['opacity'] = 0.3
                trace.showlegend = False
                combined_fig.add_trace(trace, row=row, col=col)

            fg_fig = snap.pl.umap(
                highlight_data,
                color=color_by,
                use_rep=coords,
                marker_size=pt_size,
                show=False
            )
            for trace in fg_fig.data:
                apply_trace_settings(trace)
                combined_fig.add_trace(trace, row=row, col=col)

        else:
            fig = snap.pl.umap(
                sample_data,
                color=color_by,
                use_rep=coords,
                marker_size=pt_size,
                show=False
            )
            for trace in fig.data:
                apply_trace_settings(trace)
                combined_fig.add_trace(trace, row=row, col=col)

    subplot_width = 500
    subplot_height = 500
    combined_fig.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        autosize=False,
        width=subplot_width * num_cols,
        height=subplot_height * num_rows,
        margin=dict(l=0, r=0, t=25, b=0, pad=0),
        legend=dict(
            font=dict(size=18),
            itemsizing='constant',
            itemwidth=30,
            xanchor='right',
            yanchor='top',
            x=1.1
        )
    )

    if not is_discrete:
        combined_fig.update_coloraxes(
            colorscale='Spectral_r',
            colorbar_title=color_by,
            cmin=vmin,
            cmax=vmax
        )

    for i in range(1, num_rows + 1):
        for j in range(1, num_cols + 1):
            combined_fig.update_xaxes(
                scaleratio=1,
                showticklabels=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                row=i,
                col=j
            )
            combined_fig.update_yaxes(
                showticklabels=False,
                showline=False,
                zeroline=False,
                showgrid=False,
                row=i,
                col=j
            )

    return combined_fig


def plot_volcano(
  vol_df,
  pvals_adj_threshold,
  log2fc_threshold,
  group_a,
  group_b,
  pval_key="pvals_adj",
  l2fc_key="logfoldchanges",
  names_key="name",
  title="",
  plot_width=1000,
  plot_height=425,
  top_n=2
):
    """Creates a volcano plot using Plotly with labels for the top n points by
    p-value, highest log2 fold change, and lowest log2 fold change, minimizing
    label overlap and alternating label positions.
    """
    fig_volcano_plot = go.Figure()

    # Add scatter plot
    fig_volcano_plot.add_trace(go.Scattergl(
        x=vol_df[l2fc_key],
        y=-np.log10(vol_df[pval_key]),
        mode='markers',
        marker=dict(
            size=5,
            color=np.where(
              (vol_df[pval_key] < pvals_adj_threshold)
              & (abs(vol_df[l2fc_key]) > log2fc_threshold),
              np.where(vol_df[l2fc_key] > 0, 'red', 'blue'),
              'grey'
            )
        ),
          text=vol_df[names_key],
        hoverinfo='text',
        hovertemplate='<b>%{text}</b><br>Log2FC: %{x:.2f}<br>-Log10(Adj P-value): %{y:.2f}<br>Adj P-value: %{customdata:.2e}<extra></extra>',
        customdata=vol_df[pval_key]
    ))

    # Add labels for top n points by p-value, highest and lowest log2 fold change
    top_pvals = vol_df.nsmallest(top_n, pval_key)
    top_logfc = vol_df.nlargest(top_n, l2fc_key)
    top_logfc_neg = vol_df.nsmallest(top_n, l2fc_key)
    lowest_logfc = vol_df.nlargest(top_n, l2fc_key)
    top_points = pd.concat([top_pvals, top_logfc, top_logfc_neg, lowest_logfc]).drop_duplicates()

    annotations = []
    for i, row in top_points.iterrows():
        annotations.append(
          dict(
              x=row[l2fc_key],
              y=-np.log10(row[pval_key]),
              text=row[names_key],
              showarrow=True,
              arrowhead=2,
              ax=0,
              ay=(-40 if i % 2 == 0 else 40),  # Alternate label placement
              xanchor='auto',
              yanchor='auto',
              textangle=0,
              align='center'
          )
      )

    # Adjust annotations to minimize overlap
    fig_volcano_plot.update_layout(annotations=annotations)

    # Add horizontal line for p-value threshold
    fig_volcano_plot.add_hline(
        y=-np.log10(pvals_adj_threshold), line_dash="dash", line_color="grey"
    )

    # Add vertical lines for log2 fold change thresholds
    fig_volcano_plot.add_vline(x=-log2fc_threshold, line_dash="dash", line_color="grey")
    fig_volcano_plot.add_vline(x=log2fc_threshold, line_dash="dash", line_color="grey")

    # Update layout
    fig_volcano_plot.update_layout(
        title=title,
        xaxis_title=l2fc_key,
        yaxis_title=f"-Log10({pval_key})",
        showlegend=False,
        width=plot_width,
        height=plot_height
    )

    return fig_volcano_plot


def plot_ranked_feature_plotly(
    df: pd.DataFrame,
    y_col: str,
    x_col: Optional[str] = None,
    label_col: Optional[str] = None,
    color_col: Optional[str] = None,
    n_labels: int = 30,
    colorscale: str = "Viridis",
    marker_size: int = 10,
    text_font_size: int = 10,
    ascending: bool = False,
    title: Optional[str] = None,
    y_label: Optional[str] = None,
) -> go.Figure:
    """
    Interactive Plotly scatter of y_col vs. x_col (or auto-ranked), colored by color_col or y_col,
    with top and bottom n_labels points labeled using annotations, no legend, optional title, 
    and custom axis labels.
    """
    if label_col is None:
        raise ValueError("`label_col` must be provided to draw text labels.")
    if color_col is not None and color_col not in df.columns:
        raise ValueError(f"`color_col` {color_col!r} not in DataFrame columns.")

    df_plot = df.copy()
    # Determine x-axis or compute rank
    if x_col is None:
        df_plot["rank"] = df_plot[y_col].rank(method="first", ascending=ascending)
        x = "rank"
    else:
        if x_col not in df_plot.columns:
            raise ValueError(f"`x_col` {x_col!r} not in DataFrame columns.")
        x = x_col

    # Decide which column to color by
    ccol = color_col if color_col is not None else y_col
    cbar_title = ccol.replace("_", " ").title()
    
    df_plot[x] = df_plot[x].astype(float)
    df_plot[y_col] = df_plot[y_col].astype(float)    
    df_plot[ccol] = df_plot[ccol].astype(float)

    # Create the single scatter trace for all points
    fig = go.Figure(
        go.Scatter(
            x=df_plot[x],
            y=df_plot[y_col],
            mode='markers',
            showlegend=False,
            marker=dict(
                size=marker_size,
                color=df_plot[ccol],
                colorscale=colorscale,
                showscale=True,
                opacity=0.8,
                colorbar=dict(title=cbar_title),
            ),
            hovertemplate=(
                f"{x}: %{{x}}<br>"
                f"{y_col}: %{{y}}<br>"
                f"{label_col}: %{{customdata}}<br>"
                f"{cbar_title}: %{{marker.color}}"
                "<extra></extra>"
            ),
            customdata=df_plot[label_col].astype(str).values,
        )
    )
    
    # Sort by y value to get top/bottom points
    df_sorted = df_plot.sort_values(y_col, ascending=not ascending)
    
    # Select top and bottom points for labeling
    top_points = df_sorted.head(n_labels)
    bottom_points = df_sorted.tail(n_labels) if n_labels > 0 else pd.DataFrame()

    # Function to add annotations with alternating left-right positions
    def add_point_annotations(points_df, is_top=True):
        if points_df.empty:
            return

        # Define annotation positions that alternate left and right
        # Using different offsets for top vs bottom points for better spacing
        positions = [
            {'ax': 40, 'ay': 10},    # right
            {'ax': -40, 'ay': 10},   # left
            {'ax': 40, 'ay': -10},   # left
            {'ax': -40, 'ay': -10},    # right
        ]

        for i, (_, row) in enumerate(points_df.iterrows()):
            position = positions[i % len(positions)]

            fig.add_annotation(
                x=row[x],
                y=row[y_col],
                text=str(row[label_col]),
                showarrow=True,
                arrowhead=2,
                arrowsize=1,
                arrowwidth=1,
                arrowcolor="#636363",
                ax=position['ax'],
                ay=position['ay'],
                font=dict(size=text_font_size),
                bgcolor="rgba(255, 255, 255, 0.7)",
                bordercolor="#c7c7c7",
                borderwidth=1,
                borderpad=2,
                standoff=2,
            )

    # Add annotations for top and bottom points
    add_point_annotations(top_points, is_top=True)
    add_point_annotations(bottom_points, is_top=False)

    # Calculate padding for x-axis
    x_min, x_max = df_plot[x].min(), df_plot[x].max()
    pad = (x_max - x_min) * 0.05

    # Determine axis titles
    xaxis_title = "Rank" if x_col is None else x.replace("_", " ").title()
    yaxis_title = y_label if y_label is not None else y_col.replace("_", " ").title()

    # Final layout
    layout_kwargs = dict(
        showlegend=False,
        xaxis=dict(
            title=xaxis_title,
            range=[x_min - pad, x_max + pad],
            automargin=True,
        ),
        yaxis=dict(
            title=yaxis_title,
            automargin=True,
        ),
        template="plotly_white",
        margin=dict(l=80, r=80, t=60, b=80),  # Increased margins for labels
    )
    if title:
        layout_kwargs["title"] = dict(text=title, x=0.5)

    fig.update_layout(**layout_kwargs)
    return fig


def plotly_heatmap(
  adata: AnnData,
  uns_key: str,
  key: str = "cluster",
  title: str = "",
  colorscale: str = "RdBu_r",
  width: Optional[int] = 700,
  height: Optional[int] = 700,
  mode: str = 'zscore',
  vmin: Optional[float] = None,
  vmax: Optional[float] = None,
  **kwargs: Any,
) -> go.Figure:

    # Validate input
    if key not in adata.obs:
        raise ValueError(f"Key '{key}' not found in adata.obs")
    if uns_key not in adata.uns:
        raise ValueError(f"Key '{uns_key}' not found in adata.uns")
    # Ensure the key column is categorical
    if not pd.api.types.is_categorical_dtype(adata.obs[key]):
        # Convert to categorical if it's not already
        adata.obs[key] = adata.obs[key].astype('category')

    data = adata.uns[uns_key][mode]
    categories = adata.obs[key].cat.categories

    row_labels = categories.copy()
    col_labels = categories.copy()

    # Set color scale range
    if vmin is None:
        vmin = np.nanmin(data)
    if vmax is None:
        vmax = np.nanmax(data)

    numeric_like = np.array([
        float(c) if isinstance(c, str) and c.replace('.', '', 1).isdigit() else np.inf
        for c in categories
    ])
    sorted_idx = np.argsort(numeric_like)
    data = data[sorted_idx][:, sorted_idx]
    row_labels = row_labels[sorted_idx]
    col_labels = col_labels[sorted_idx]

    fig = go.Figure()

    # Create heatmap
    heatmap = go.Heatmap(
        z=data,
        x=col_labels,
        y=row_labels,
        colorscale=colorscale,
        showscale=True,
        textfont={"size": 12},
        hoverongaps=False,
        zmin=vmin,  # Set minimum color scale value
        zmax=vmax,  # Set maximum color scale value
        colorbar=dict(
            title=key,
            tickformat='.2f',
            len=0.9
        ),
        **kwargs
    )
    fig.add_trace(heatmap)
    # Update layout
    fig.update_layout(
        title={
            'text': title,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        width=width,
        height=height,
        showlegend=False,
        xaxis=dict(
            showgrid=False,
            side='bottom'
        ),
        yaxis=dict(
            showgrid=False,
            autorange='reversed'
        ),
    )

    return fig, data


def process_matrix_layout(
    adata_all,
    n_rows: int,
    n_cols: int,
    sample_key: str = "sample",
    spatial_key: str = "spatial",
    new_obsm_key: str = "X_dataset",
    tile_spacing: float = 100.0,
    flipy: bool = False,
    sample_order_mode: str = "original",  # "original", "sample", or "condition"
    condition_key: str = "condition",
):
    """
    Add new obsm with offset spatial coords to display all samples at once.

    Parameters
    ----------
    adata_all : AnnData
    n_rows, n_cols : int
        Grid size for placing samples.
    sample_key : str
        Column in .obs with sample names.
    spatial_key : str
        Key in .obsm for input spatial coordinates (N x 2).
    new_obsm_key : str
        Key in .obsm to store transformed coordinates.
    tile_spacing : float
        Spacing between tiles.
    flipy : bool
        If True, vertically flip each sample across the y axis.
    sample_order_mode : str
        "original" (pd.unique order), "sample" (alphabetical by sample),
        or "condition" (condition A..Z, then sample A..Z).
    condition_key : str
        Column in .obs with condition labels (only used when sample_order_mode="condition").
    """

    # --- Decide sample placement order ---
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

    # --- Validate grid dims vs number of samples ---
    n_samples = len(samples)
    if n_cols is not None and n_rows is not None:
        total_positions = n_rows * n_cols
        if n_samples > total_positions:
            raise ValueError(
                f"Not enough grid positions ({n_rows}x{n_cols}={total_positions}) for {n_samples} samples"
            )
    elif n_cols is not None and n_rows is None:
        n_rows = (n_samples + n_cols - 1) // n_cols
    elif n_rows is not None and n_cols is None:
        n_cols = (n_samples + n_rows - 1) // n_rows

    total_positions = n_rows * n_cols
    if n_samples > total_positions:
        raise ValueError(f"Not enough grid positions ({total_positions}) for {n_samples} samples")

    # --- Prepare containers ---
    X_new = np.empty_like(adata_all.obsm[spatial_key], dtype=float)

    grid_bounds = {}
    sample_positions = {}

    # First pass: bounds
    for idx, sample_name in enumerate(samples):
        row = idx // n_cols
        col = idx % n_cols
        sample_positions[sample_name] = (row, col)

        mask = (adata_all.obs[sample_key].astype(str) == str(sample_name))
        xspa = adata_all.obsm[spatial_key][mask]

        l_max = xspa.max(axis=0)
        l_min = xspa.min(axis=0)
        width = float(l_max[0] - l_min[0])
        height = float(l_max[1] - l_min[1])

        grid_bounds[(row, col)] = {
            "width": width,
            "height": height,
            "min_x": float(l_min[0]),
            "min_y": float(l_min[1]),
            "max_x": float(l_max[0]),
            "max_y": float(l_max[1]),
        }

    # Row/col offsets
    row_heights = [max((grid_bounds[(r, c)]["height"] for c in range(n_cols) if (r, c) in grid_bounds), default=0.0) for r in range(n_rows)]
    col_widths = [max((grid_bounds[(r, c)]["width"] for r in range(n_rows) if (r, c) in grid_bounds), default=0.0) for c in range(n_cols)]
    row_y_offsets = [0.0]
    for i in range(n_rows - 1):
        row_y_offsets.append(row_y_offsets[-1] - row_heights[i] - tile_spacing)
    col_x_offsets = [0.0]
    for i in range(n_cols - 1):
        col_x_offsets.append(col_x_offsets[-1] + col_widths[i] + tile_spacing)

    # Second pass: transform coords
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

        dif_x = target_x - bounds["min_x"]
        dif_y = target_y - (bounds["max_y"] if not flipy else bounds["min_y"])

        xspa[:, 0] += dif_x
        xspa[:, 1] += dif_y

        X_new[mask] = xspa

    adata_all.obsm[new_obsm_key] = X_new


def rename_obs_keys(adata: anndata.AnnData) -> anndata.AnnData:
    """Add obs columns for Plots by renaming keys, if needed."""
    key_map = {
        "Sample": "sample",
        "nFrags": "n_fragment",
        "Condition": "condition",
        "Clusters": "cluster"
    }

    keys = adata.obs_keys()
    for src, dest in key_map.items():
        if src in keys:
            if dest not in keys:
                adata.obs[dest] = adata.obs[src]
            else:
                print(
                    f"Target key '{dest}' already exists in obs; skipping \
                    copy from '{src}'."
                )

    return adata


def rgb_to_hex(rgb):
    """Convert RGB tuple to hex color code"""
    return '#{:02x}{:02x}{:02x}'.format(
        int(rgb[0] * 255),
        int(rgb[1] * 255),
        int(rgb[2] * 255)
    )


def reorder_obs_columns(adata, first_col="cluster"):
    """Move specified column to first position in obs DataFrame."""
    if first_col not in adata.obs.columns:
        return
    new_order = [first_col] + [c for c in adata.obs.columns if c != first_col]
    adata.obs = adata.obs[new_order]


def safe_float(val, warn_msg=None):
  """Validate for umap plots custom max/mins."""
  if val == "":
      return None
  try:
      return float(val)
  except (TypeError, ValueError):
      if warn_msg:
        w_text_output(
            content=warn_msg,
            appearance={"message_box": "warning"}
        )
      else:
        pass
      return None


def sort_group_categories(values):
    """Sort group labels numerically if possible, else alphabetically."""
    # Try to convert to numbers
    num_vals = []
    all_numeric = True
    for v in values:
        try:
            num_vals.append(float(v))
        except (ValueError, TypeError):
            all_numeric = False
            break

    if all_numeric:
        # Sort by numeric value
        sorted_pairs = sorted(zip(values, num_vals), key=lambda x: x[1])
        return [v for v, _ in sorted_pairs]
    else:
        # Fall back to string sort
        return sorted(map(str, values))


def squidpy_analysis(
    adata: anndata.AnnData,
    cluster_key: str = "cluster",
    sample_key: Optional[str] = None
) -> anndata.AnnData:
    """Perform squidpy Neighbors enrichment analysis.
    """
    from squidpy.gr import nhood_enrichment, spatial_neighbors

    if not adata.obs[cluster_key].dtype.name == "category":
        adata.obs[cluster_key] = adata.obs["cluster"].astype("category")

    if sample_key:
        if not adata.obs[sample_key].dtype.name == "category":
            adata.obs[sample_key] = adata.obs[sample_key].astype("category")

    spatial_neighbors(
        adata, coord_type="grid", n_neighs=4, n_rings=1, library_key=sample_key
    )
    nhood_enrichment(
        adata, cluster_key=cluster_key, library_key=sample_key, seed=42
    )

    return adata


def sync_obs_metadata(
    adata1: AnnData,
    adata2: AnnData,
    *,
    reconcile_shared: bool = True,
    reconcile_direction: str = "adata1_to_adata2",
) -> None:
    """
    Reciprocally copy any *non-numeric* obs columns so both AnnData objects end up
    with the same set of obs columns. Safe to differing cell order: aligns by obs_names.
    Operates in-place and does NOT modify columns that already exist in a target.

    Optionally detects differences within shared columns and reconciles them.

    Args:
        adata1, adata2: AnnData objects to synchronize.
        reconcile_shared: If True, compare shared non-numeric columns and overwrite
            one side when values differ. If False, only copies missing columns.
        reconcile_direction:
            - "adata1_to_adata2": overwrite adata2 from adata1 when different
            - "adata2_to_adata1": overwrite adata1 from adata2 when different

    Raises:
        ValueError: if the two objects don't contain exactly the same cells (by name).
        ValueError: if reconcile_direction is invalid.
    """
    idx1 = pd.Index(adata1.obs_names)
    idx2 = pd.Index(adata2.obs_names)

    if not idx1.equals(idx2):
        if set(idx1) != set(idx2):
            raise ValueError("AnnData objects must contain exactly the same cells (obs_names).")

    if reconcile_direction not in {"adata1_to_adata2", "adata2_to_adata1"}:
        raise ValueError("reconcile_direction must be 'adata1_to_adata2' or 'adata2_to_adata1'.")

    obs1 = adata1.obs
    obs2 = adata2.obs

    # Identify non-numeric columns in each .obs
    nonnum1 = [c for c in obs1.columns if not pd.api.types.is_numeric_dtype(obs1[c])]
    nonnum2 = [c for c in obs2.columns if not pd.api.types.is_numeric_dtype(obs2[c])]

    # Determine which columns are missing in the other object
    missing_in_2 = [c for c in nonnum1 if c not in obs2.columns]
    missing_in_1 = [c for c in nonnum2 if c not in obs1.columns]

    # Copy from adata1 -> adata2 (aligned by cell names)
    for col in missing_in_2:
        adata2.obs[col] = obs1[col].reindex(adata2.obs_names)

    # Copy from adata2 -> adata1 (aligned by cell names)
    for col in missing_in_1:
        adata1.obs[col] = obs2[col].reindex(adata1.obs_names)

    # Optionally reconcile differences in shared non-numeric columns
    if not reconcile_shared:
        return

    shared_nonnum = [c for c in nonnum1 if c in nonnum2]

    # Compare on a common order to avoid order-related mismatches
    common_idx = pd.Index(adata1.obs_names)  # order doesn't matter as long as consistent

    for col in shared_nonnum:
        s1 = obs1[col].reindex(common_idx)
        s2 = obs2[col].reindex(common_idx)

        differs = False
        try:
            differs = not s1.equals(s2)  # equals handles NaNs
        except Exception:
            differs = True

        if differs:
            if reconcile_direction == "adata1_to_adata2":
                adata2.obs[col] = obs1[col].reindex(adata2.obs_names)
            else:  # "adata2_to_adata1"
                adata1.obs[col] = obs2[col].reindex(adata1.obs_names)
