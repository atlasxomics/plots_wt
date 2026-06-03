w_text_output(content="""

# Violin Plot

Generate violin plots to visualize the distribution of values across a selected grouping (e.g., **Clusters**, **Samples**, or **Conditions**).  
You can display numeric cell-level metrics or feature values from `.X`.

> Feature names come from `.var_names`.

""")

new_data_signal()

if adata_g is None:
    w_text_output(
        content="No data selected...",
        appearance={"message_box": "warning"}
    )
    exit()

notebook_palettes = await get_notebook_palettes()
categorical_palette = w_select(
  label="palette",
  default=DEFAULT_CATEGORICAL_PALETTE_NAME,
  options=get_palette_selector_options(notebook_palettes),
  appearance={
    "help_text": "Use a palette saved from the H5 Viewer or fall back to the default palette."
  }
)

violin_groups = get_groupable_obs_keys(adata_g)

available_metadata = tuple(
  key for key in adata_g.obs_keys()
  if key not in na_keys
)

numeric_metadata = [
  data for data in get_numeric_obs_keys(adata_g)
  if data in available_metadata
]
has_features = adata_g.n_vars > 0

if not violin_groups:
  w_text_output(
    content="No low-cardinality categorical metadata columns are available for grouping.",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()

if not numeric_metadata and not has_features:
  w_text_output(
    content="No numeric metadata columns or features are available to plot.",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()

data_source_options = []
if numeric_metadata:
  data_source_options.append("metadata")
if has_features:
  data_source_options.append("feature")

violin_data_source = w_select(
  label="data source",
  default=choose_default_option(data_source_options, preferred="metadata"),
  options=tuple(data_source_options),
  appearance={
    "help_text": "Choose whether to plot obs metadata or a feature from `.X`."
  }
)

violin_metadata = None
violin_feature = None
if violin_data_source.value == "metadata":
  if not numeric_metadata:
    w_text_output(
      content="No numeric metadata columns are available to plot.",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()
  violin_metadata = w_select(
    label="metadata",
    default=choose_default_option(numeric_metadata, preferred="tsse"),
    options=tuple(numeric_metadata),
    appearance={
      "detail": "numeric obs",
      "help_text": "Select numeric observation metadata."
    }
  )
else:
  violin_feature = w_text_input(
    label="feature",
    key="violin_feature",
    default="",
    appearance={
      "placeholder": "Enter feature name",
      "help_text": "Type an exact `.var_names` feature name."
    }
  )

violin_group_by = w_select(
  label="group",
  default=choose_group_default(
    violin_groups,
    preferred=("cluster", "condition", "sample"),
  ),
  options=tuple(violin_groups),
  appearance={
    "detail": "categorical obs",
    "help_text": "Select group to display on x-axis."
  }
)

violin_type = w_select(
  label="plot type",
  default="box",
  options=tuple(["box", "violin"]),
  appearance={
    "detail": "(box, violin)",
    "help_text": "Use box to speed up large datasets."
  }
)
violin_max_points = w_text_input(
  label="max plotted cells",
  key="violin_max_points",
  default=str(DEFAULT_VIOLIN_MAX_POINTS),
  appearance={
    "help_text": "Caps plotting data with group-stratified sampling."
  }
)

data_widget = violin_metadata if violin_metadata is not None else violin_feature
violin_row = w_row(items=[
  violin_data_source,
  data_widget,
  violin_group_by,
  violin_type,
  violin_max_points,
  categorical_palette
])

data_type = None
plot_data = None
if violin_data_source.value == "metadata" and violin_metadata is not None:
  data_type = "obs"
  plot_data = violin_metadata.value
elif violin_data_source.value == "feature" and violin_feature is not None:
  data_type = "feature"
  plot_data = (violin_feature.value or "").strip()

if not data_type or not plot_data:
  w_text_output(
    content="Select metadata or enter a feature name to plot.",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()

if data_type == "feature" and plot_data not in adata_g.var_names:
  w_text_output(
    content=f"Feature `{plot_data}` was not found in `.var_names`.",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()

if data_type is not None:

  v_adata = adata_g
  group_col = violin_group_by.value
  
  try:
      max_points = parse_positive_int(
          violin_max_points.value,
          DEFAULT_VIOLIN_MAX_POINTS,
          "Max plotted cells",
      )
      violin_obs_positions = stratified_obs_positions(
          v_adata.obs[group_col].astype(str),
          max_points=max_points,
      )
      violin_df = create_violin_data(
          v_adata,
          group_col,
          plot_data,
          data_type=data_type,
          obs_positions=violin_obs_positions,
      )
  except KeyError as e:
      w_text_output(
          content=(
              f"Could not find {e} in the AnnData object."
          ),
          appearance={"message_box": "warning"},
      )
      submit_widget_state()
      exit()
  except ValueError as e:
      w_text_output(
          content=str(e),
          appearance={"message_box": "warning"},
      )
      submit_widget_state()
      exit()

  if len(violin_obs_positions) < v_adata.n_obs:
    w_text_output(
      content=(
        f"Plotting {len(violin_obs_positions):,} group-stratified cells "
        f"from {v_adata.n_obs:,} total cells."
      ),
      appearance={"message_box": "info"}
    )

  if violin_type.value == "box":
    violin_categories = sort_group_categories(violin_df['group'].unique().tolist())
    selected_colors = get_selected_palette_colors(
        notebook_palettes,
        categorical_palette.value,
        fallback_colors=DEFAULT_H5_CATEGORICAL_PALETTE,
    )
    violin_color_map = build_discrete_color_map(violin_categories, selected_colors)

    violin_fig = px.box(
        violin_df,
        x='group',
        y='value',
        points=False,
        color='group',
        color_discrete_map=violin_color_map,
        category_orders={"group": violin_categories},
    )

  elif violin_type.value == "violin":
    violin_categories = sort_group_categories(violin_df['group'].unique().tolist())
    selected_colors = get_selected_palette_colors(
        notebook_palettes,
        categorical_palette.value,
        fallback_colors=DEFAULT_H5_CATEGORICAL_PALETTE,
    )
    violin_color_map = build_discrete_color_map(violin_categories, selected_colors)

    violin_fig = px.violin(
        violin_df,
        x='group',
        y='value',
        box=True,
        points=False,
        color='group',
        color_discrete_map=violin_color_map,
        category_orders={"group": violin_categories},
      )

  else:
    w_text_output(
      content="Plot type not recognized",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

  violin_fig.update_layout(
      title=f"Distribution of {plot_data} by {violin_group_by.value}",
      xaxis_title=violin_group_by.value,
      yaxis_title=plot_data,
      plot_bgcolor='rgba(0,0,0,0)',
      showlegend=False
  )
  violin_fig.update_xaxes(showgrid=False)
  violin_fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')

  # Determine ordering: numeric sort if possible, else alphabetical
  v_category_order = violin_categories
  
  violin_fig.update_xaxes(categoryorder='array', categoryarray=v_category_order)

  w_plot(source=violin_fig)

  violin_data_button = w_checkbox(
    label="Display Violin Data",
    key="violin_data_button",
    default=False,
  )

  if violin_data_button.value:
    violin_table = w_table(
      label=f"Distribution of {plot_data} by {violin_group_by.value}",
      source=violin_df
    )

else:
  w_output_text(content="  ")
  sumbit_widget_state()
