new_data_signal()

w_text_output(content="""

# Heatmaps

Display heatmaps of differential statistics across project groupings (e.g., **Clusters**, **Samples**, or **Conditions**) for either genes or motifs.

Heatmaps are generated dynamically from statistics tables stored in `.uns` (for example, `ranked_genes_per_*`, `marker_genes_per_*`, and `enrichedMotifs_*`), so you can adjust filters on demand.

You can:
- set a significance threshold (FDR/adjusted p-value when available),
- use ArchR-like marker filtering:
  Log2FC is used for both feature ranking and heatmap values, with directional effect-size filtering,
- use ArchR-like display scaling:
  genes are row z-scored and clipped (default ±2), motifs use row-normalized `-log10(adj p-value)`,
- select top features per group,
- or provide an explicit comma-separated feature list.

""")

# Abort if no data loaded
if not adata_g:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

# Choose whether to display gene or motif data
choose_heatmap_data = w_select(
    label="Select Data for Heatmap Plots",
    default="gene",
    options=["gene", "motif"],
    appearance={
        "help_text": "Select which features to display in the heatmap."
    }
)

if choose_heatmap_data.value is not None:
  heatmap_signal(True)
