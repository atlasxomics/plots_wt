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

violin_groups = get_categorical_obs_keys(adata_g)

available_metadata = tuple(
  key for key in adata_g.obs_keys()
  if key not in na_keys
)

numeric_metadata = [
  data for data in get_numeric_obs_keys(adata_g)
  if data in available_metadata
]
feature_options = available_features if "available_features" in globals() else available_genes
violin_data_options = list(dict.fromkeys(numeric_metadata + feature_options))

if not violin_groups:
  w_text_output(
    content="No categorical metadata columns available for grouping.",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()

if not violin_data_options:
  w_text_output(
    content="No numeric metadata columns or features are available to plot.",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()

violin_data = w_select(
  label="data",
  default=choose_default_option(violin_data_options, preferred="tsse"),
  options=tuple(violin_data_options),
  appearance={
    "help_text": "Select values to plot.",
    "description": "Metadata columns come from `.obs`; features come from `.var_names`."
  }
)

violin_group_by = w_select(
  label="group",
  default=choose_default_option(violin_groups, preferred="cluster"),
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
violin_row = w_row(items=[violin_data, violin_group_by, violin_type, categorical_palette])

data_type = None
if violin_data.value in numeric_metadata:
  data_type = "obs"
elif violin_data.value in feature_options:
  data_type = "feature"

if not data_type:
  w_text_output(
    content="Selected data not found in AnnData object...",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()

if data_type is not None:

  v_adata = adata_g
  group_col = violin_group_by.value
  
  try:
      violin_df = create_violin_data(
          v_adata, group_col, violin_data.value, data_type=data_type
      )
  except KeyError:
      w_text_output(
          content=(
              f"Could not find {group_col} in the AnnData object."
          ),
          appearance={"message_box": "warning"},
      )
      submit_widget_state()
      exit()

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
      title=f"Distribution of {violin_data.value} by {violin_group_by.value}",
      xaxis_title=violin_group_by.value,
      yaxis_title=violin_data.value,
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
      label=f"Distribution of {violin_data.value} by {violin_group_by.value}",
      source=violin_df
    )

else:
  w_output_text(content="  ")
  sumbit_widget_state()
