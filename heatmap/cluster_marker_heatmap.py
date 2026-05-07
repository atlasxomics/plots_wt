new_data_signal()

w_text_output(content="""

# Cluster Marker Heatmap

Display the marker-gene heatmap stored in optimized whole-transcriptome
`combined.h5ad` files.

""")

if adata_g is None or not isinstance(adata_g, AnnData):
    w_text_output(
        content="No data loaded...",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

heatmap_key = "cluster_marker_heatmap"
heatmap_params_key = "cluster_marker_heatmap_params"

if heatmap_key not in adata_g.uns:
    w_text_output(
        content=(
            "No cluster marker heatmap was found in this AnnData object. "
            f"Expected `adata.uns['{heatmap_key}']`, which is added by "
            "`optimize_wt` when cluster marker computation succeeds."
        ),
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

raw_heatmap = adata_g.uns[heatmap_key]

try:
    if isinstance(raw_heatmap, pd.DataFrame):
        marker_heatmap_df = raw_heatmap.copy()
    else:
        marker_heatmap_df = pd.DataFrame(raw_heatmap)
except Exception as e:
    w_text_output(
        content=f"Could not convert `adata.uns['{heatmap_key}']` to a table: {e}",
        appearance={"message_box": "danger"}
    )
    submit_widget_state()
    exit()

marker_heatmap_df = marker_heatmap_df.apply(pd.to_numeric, errors="coerce")
marker_heatmap_df = marker_heatmap_df.dropna(axis=0, how="all").dropna(axis=1, how="all")
marker_heatmap_df.index = marker_heatmap_df.index.map(str)
marker_heatmap_df.columns = marker_heatmap_df.columns.map(str)

if marker_heatmap_df.empty:
    w_text_output(
        content=f"`adata.uns['{heatmap_key}']` is empty after removing missing values.",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

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

show_table = w_checkbox(
    key="cluster_marker_heatmap_table",
    label="Display heatmap data table",
    default=False,
)

w_row(items=[hm_palette, show_table])

heatmap_params = adata_g.uns.get(heatmap_params_key, {})
legend_title = "Z-score"
title = "Cluster Marker Heatmap"

if isinstance(heatmap_params, dict):
    marker_top_n = heatmap_params.get("marker_top_n")
    values_label = heatmap_params.get("values")
    if marker_top_n is not None:
        title = f"Top {marker_top_n} Marker Genes per Cluster"
    if isinstance(values_label, str) and values_label:
        w_text_output(
            content=f"Values: {values_label}",
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
        title=legend_title,
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

w_plot(source=cluster_marker_heatmap)

if show_table.value:
    w_table(
        label="Cluster marker heatmap data",
        source=marker_heatmap_df,
    )
