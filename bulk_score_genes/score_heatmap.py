new_data_signal()
if "gene_score_done_signal" not in globals():
  gene_score_done_signal = Signal(False)
gene_score_done_signal()

if adata_g is None or not isinstance(adata_g, AnnData):
  w_text_output(content=" ")
  exit()

if gene_score_done_signal.sample() is not True:
  w_text_output(content=" ")
  exit()

if "score_cols" not in globals() or not score_cols:
  w_text_output(
      content=" ",
  )
  exit()

available_score_cols = [
    col for col in score_cols
    if col in adata_g.obs.columns
  ]

categorical_obs = get_categorical_obs_keys(adata_g)

if not categorical_obs:
  w_text_output(
      content="No categorical metadata columns available to group by.",
      appearance={"message_box": "warning"}
  )
  exit()

notebook_palettes = await get_notebook_palettes()
score_palette = w_select(
    label="Colorscale",
    default="Default Gene Score Heatmap Colorscale",
    options=get_palette_selector_options(
        notebook_palettes,
        kind="continuous",
        fallback_name="Default Gene Score Heatmap Colorscale",
    ),
    appearance={
        "help_text": "Use a continuous palette saved from the H5 Viewer or fall back to the default gene score heatmap colors."
    },
    key="score_heatmap_palette"
)

score_select = w_multi_select(
    label="Score columns",
    options=tuple(available_score_cols),
    appearance={"help_text": "Select one or more gene score columns for the heatmap."},
    key="score_heatmap_cols"
)

group_select_gs = w_select(
    label="Group by metadata column",
    default=choose_default_option(categorical_obs, preferred="cluster"),
    options=tuple(categorical_obs),
    appearance={"help_text": "Select a categorical obs column to aggregate scores by."},
    key="score_heatmap_group"
)

gs_row = w_row(items=[score_select, group_select_gs, score_palette])

selected_scores = list(score_select.value or [])
selected_group_gs = group_select_gs.value

if not selected_scores or selected_group_gs is None:
  w_text_output(
      content="Select at least one score column and a grouping column to view the heatmap.",
      appearance={"message_box": "info"}
  )
  exit()

obs_df = adata_g.obs[selected_scores + [selected_group_gs]].copy()
obs_df = obs_df.dropna(subset=[selected_group_gs])

if obs_df.empty:
  w_text_output(
      content="No observations available after filtering missing group values.",
      appearance={"message_box": "warning"}
  )
  exit()

group_series = adata_g.obs[selected_group_gs]
if pd.api.types.is_categorical_dtype(group_series):
  category_order = [str(cat) for cat in group_series.cat.categories]
  obs_df[selected_group_gs] = group_series.astype(str)
else:
  category_order = None
  obs_df[selected_group_gs] = obs_df[selected_group_gs].astype(str)

heatmap_df_gs = (
    obs_df.groupby(selected_group_gs)[selected_scores]
    .mean()
    .T
)

if heatmap_df_gs.empty:
  w_text_output(
      content="Unable to compute aggregated scores for the selected inputs.",
      appearance={"message_box": "warning"}
  )
  exit()

if category_order:
  ordered_cols = [cat for cat in category_order if cat in heatmap_df_gs.columns]
  remaining_cols = [col for col in heatmap_df_gs.columns if col not in ordered_cols]
  heatmap_df_gs = heatmap_df_gs[ordered_cols + remaining_cols]

heatmap_title_gs = f"Gene Score Heatmap by {selected_group_gs}"

heatmap_gs = px.imshow(
    heatmap_df_gs,
    color_continuous_scale=get_selected_continuous_palette(
        notebook_palettes,
        score_palette.value,
        fallback_colors="Spectral_r",
        fallback_name="Default Gene Score Heatmap Colorscale",
    ),
    aspect="auto",
    origin="lower"
)

heatmap_gs.update_layout(
    title=heatmap_title_gs,
    xaxis_title=selected_group_gs,
    yaxis_title="Gene Score",
    coloraxis_colorbar=dict(
        title="Score",
        title_side="right"
    )
)

heatmap_gs.update_xaxes(
    side="bottom",
    tickmode="array",
    tickvals=list(range(len(heatmap_df_gs.columns))),
    ticktext=heatmap_df_gs.columns,
    tickangle=45
)
heatmap_gs.update_yaxes(
    autorange="reversed",
    tickmode="array",
    tickvals=list(range(len(heatmap_df_gs.index))),
    ticktext=heatmap_df_gs.index
)

w_plot(source=heatmap_gs)
