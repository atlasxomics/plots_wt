new_data_signal()

w_text_output(content="""

# Proportion Plot

Generate stacked bar plots where the **x-axis** represents a primary grouping (e.g., conditions) and the **stacked segments** represent sub-groups (e.g., clusters).  
The **y-axis** can be displayed as either raw counts or proportions of cells.

""")

# Abort if no data loaded
if adata_g is None:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
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

prop_groups = get_categorical_obs_keys(adata_g)

if not prop_groups:
    w_text_output(
        content="No categorical metadata columns available for proportion plotting.",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

default_group_by = choose_default_option(prop_groups, preferred="sample")
default_stack_by = choose_default_option(
    prop_groups,
    preferred="cluster",
    fallback=next((group for group in prop_groups if group != default_group_by), None),
)

group_by = w_select(
    label="group by",
    default=default_group_by,
    options=tuple(prop_groups),
    appearance={
        "detail": "categorical obs",
        "help_text": "Select group to display on x-axis."
    }
)

stack_by = w_select(
    label="stack by",
    default=default_stack_by,
    options=tuple(prop_groups),
    appearance={
        "detail": "categorical obs",
        "help_text": "Select group to stack the bars by."
    }
)
return_type = w_select(
    label="return type",
    default="proportion",
    options=("proportion", "counts"),
    appearance={
        "help_text": "Display cell counts or proportions per grouping."
    }
)

prop_row = w_row(items=[group_by, stack_by, return_type, categorical_palette])

stacked_df = create_proportion_dataframe(
    adata_g, group_by.value, stack_by.value, return_type=return_type.value
)

stack_categories = sort_group_categories(stacked_df["stack_by"].unique().tolist())
selected_colors = get_selected_palette_colors(
    notebook_palettes,
    categorical_palette.value,
    fallback_colors=DEFAULT_H5_CATEGORICAL_PALETTE,
)
stack_color_map = build_discrete_color_map(stack_categories, selected_colors)

# Create proportion plot from the stacked dataframe created above
proportion_plot = px.bar(
    stacked_df,
    x="group_by",
    y="value",
    color="stack_by",
    barmode="stack",
    color_discrete_map=stack_color_map,
    category_orders={"stack_by": stack_categories},
    title=f"Distribution of {stack_by.value} by {group_by.value}"
)

# Update layout
proportion_plot.update_layout(
    xaxis_title=group_by.value,
    yaxis_title="Proportion" if return_type.value == "proportion" else "Count",
    plot_bgcolor='rgba(0,0,0,0)',
    showlegend=True,
    legend_title=stack_by.value
)

# Update axes
proportion_plot.update_xaxes(showgrid=False)
proportion_plot.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')

# Determine ordering: numeric sort if possible, else alphabetical
p_cats = stacked_df['group_by'].unique().tolist()
p_category_order = sort_group_categories(p_cats)

proportion_plot.update_xaxes(categoryorder='array', categoryarray=p_category_order)

w_plot(source=proportion_plot)
