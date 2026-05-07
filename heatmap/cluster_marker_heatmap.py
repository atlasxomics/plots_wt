new_data_signal()

w_text_output(content="""

# Cluster Marker Heatmap

Create a marker-gene heatmap from cluster DEG results stored in optimized
whole-transcriptome `combined.h5ad` files.

""")

if adata_g is None or not isinstance(adata_g, AnnData):
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

deg_key = "cluster_marker_degs"
deg_params_key = "cluster_marker_degs_params"
heatmap_key = "cluster_marker_heatmap"
heatmap_params_key = "cluster_marker_heatmap_params"

notebook_palettes = await get_notebook_palettes()

hm_palette = w_select(
    key="cluster_marker_heatmap_palette",
    label="Colorscale",
    default="Default Cluster Marker Heatmap Colorscale",
    options=get_palette_selector_options(
        notebook_palettes,
        kind="continuous",
        fallback_name="Default Cluster Marker Heatmap Colorscale",
    ),
    appearance={
        "help_text": "Use a continuous palette saved from the H5 Viewer or fall back to the default marker heatmap colors."
    }
)

heatmap_params = adata_g.uns.get(heatmap_params_key, {})
deg_params = adata_g.uns.get(deg_params_key, {})
default_top_n = 50
if isinstance(heatmap_params, dict) and heatmap_params.get("marker_top_n") is not None:
    default_top_n = int(heatmap_params.get("marker_top_n"))
elif isinstance(deg_params, dict) and deg_params.get("marker_top_n") is not None:
    default_top_n = int(deg_params.get("marker_top_n"))

top_n_input = w_text_input(
    key="cluster_marker_heatmap_top_n",
    label="Top genes per cluster",
    default=str(default_top_n),
    appearance={"help_text": "Number of top DEG rows per cluster to include."},
)

pval_options = ["pvals"]
if deg_key in adata_g.uns:
    candidate_degs = cluster_marker_to_dataframe(adata_g.uns[deg_key], deg_key)
    if "pvals_adj" in candidate_degs.columns:
        pval_options.append("pvals_adj")
else:
    candidate_degs = None

pval_col_select = w_select(
    key="cluster_marker_heatmap_pval_col",
    label="P-value column",
    default="pvals",
    options=tuple(pval_options),
)

pval_cutoff_input = w_text_input(
    key="cluster_marker_heatmap_pval_cutoff",
    label="P-value cutoff",
    default="0.05",
)

log2fc_cutoff_input = w_text_input(
    key="cluster_marker_heatmap_log2fc_cutoff",
    label="Log2FC cutoff",
    default="0.25",
)

order_select = w_select(
    key="cluster_marker_heatmap_order",
    label="Cluster order",
    default="DEG similarity",
    options=("DEG similarity", "numeric", "user-selected"),
)

default_order = ""
if "cluster" in adata_g.obs:
    default_order = ", ".join(
        sorted(
            adata_g.obs["cluster"].astype(str).unique().tolist(),
            key=cluster_marker_sort_key,
        )
    )

user_order_input = w_text_input(
    key="cluster_marker_heatmap_user_order",
    label="User cluster order",
    default=default_order,
    appearance={
        "help_text": "Comma-separated cluster order. Used only when Cluster order is user-selected."
    },
)

w_row(items=[hm_palette, top_n_input, pval_col_select])
w_row(items=[pval_cutoff_input, log2fc_cutoff_input])
w_row(items=[user_order_input, order_select])

try:
    if candidate_degs is not None:
        top_n = parse_cluster_marker_int(
            top_n_input.value,
            default_top_n,
            "Top genes per cluster",
        )
        pval_cutoff = parse_cluster_marker_float(
            pval_cutoff_input.value,
            0.05,
            "P-value cutoff",
            minimum=0.0,
        )
        log2fc_cutoff = parse_cluster_marker_float(
            log2fc_cutoff_input.value,
            0.25,
            "Log2FC cutoff",
        )
        marker_heatmap_df = compute_cluster_marker_heatmap_from_degs(
            adata_g,
            candidate_degs,
            top_n=top_n,
            pval_col=pval_col_select.value,
            pval_cutoff=pval_cutoff,
            log2fc_cutoff=log2fc_cutoff,
            order_mode=order_select.value,
            user_order=user_order_input.value,
            deg_key=deg_key,
        )
        title = f"Top {top_n} Marker Genes per Cluster"
    else:
        marker_heatmap_df = get_cached_cluster_marker_heatmap(
            adata_g,
            heatmap_key=heatmap_key,
        )
        if marker_heatmap_df is None:
            raise ValueError(
                f"No dynamic DEG table (`adata.uns['{deg_key}']`) or cached heatmap "
                f"(`adata.uns['{heatmap_key}']`) was found."
            )
        w_text_output(
            content=(
                f"Differential stats not found, defauling to cached heatmap. "
                f"Some options will not work. "
            ),
            appearance={"message_box": "warning"}
        )
        title = "Cached Cluster Marker Heatmap"
except Exception as e:
    w_text_output(
        content=str(e),
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

if isinstance(deg_params, dict) and deg_params:
    fixed_filters = []
    if deg_params.get("excluded_prefixes") is not None:
        fixed_filters.append(f"excluded prefixes: {deg_params.get('excluded_prefixes')}")
    if deg_params.get("expression_layer") is not None:
        fixed_filters.append(f"expression layer: {deg_params.get('expression_layer')}")
    if fixed_filters:
        w_text_output(
            content="Stored DEG table filters: " + "; ".join(fixed_filters),
            appearance={"message_box": "info"}
        )

gene_labels = marker_heatmap_df.columns.tolist()
label_every_n = max(1, math.ceil(len(gene_labels) / 80))
visible_gene_labels = [
    gene
    for i, gene in enumerate(gene_labels)
    if i % label_every_n == 0
]

cluster_marker_heatmap = px.imshow(
    marker_heatmap_df,
    color_continuous_scale=get_selected_continuous_palette(
        notebook_palettes,
        hm_palette.value,
        fallback_colors="RdYlBu_r",
        fallback_name="Default Cluster Marker Heatmap Colorscale",
    ),
    aspect="auto",
    origin="lower",
    zmin=-3,
    zmax=3,
)

cluster_marker_heatmap.update_layout(
    title=title,
    xaxis_title="Marker gene",
    yaxis_title="Cluster",
    coloraxis_colorbar=dict(
        title="Z-score",
        title_side="right"
    )
)

cluster_marker_heatmap.update_xaxes(
    side="bottom",
    tickmode="array",
    tickvals=visible_gene_labels,
    ticktext=visible_gene_labels,
    tickangle=90
)
cluster_marker_heatmap.update_yaxes(
    autorange="reversed",
    tickmode="array",
    tickvals=marker_heatmap_df.index.tolist(),
    ticktext=marker_heatmap_df.index.tolist()
)

cluster_hm_plot = w_plot(source=cluster_marker_heatmap)

cluster_hm_table = w_table(
    label="Cluster marker heatmap data",
    source=marker_heatmap_df,
)
