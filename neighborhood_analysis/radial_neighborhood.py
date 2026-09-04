new_data_signal()

# The selector is rendered by the heatmap cell above. Reading its value makes
# this cell rerun when the user switches display type.
_neigh_display_value = (
    neigh_display.value if "neigh_display" in globals() else "heatmap"
)
if _neigh_display_value != "radial":
    exit()

w_text_output(content="""
## Radial Neighborhood Enrichment

This is a spoke view of the same enrichment shown in the heatmap. Each subplot is one **focal cluster**, with spokes pointing to neighboring clusters around the circle. **Red** indicates positive enrichment (co-localization), **blue** indicates negative enrichment (avoidance), and spoke thickness scales with absolute z-score. Hover over a spoke to see its value.
""")

if adata_g is None:
    w_text_output(
        content="No AnnData file selected...",
        appearance={"message_box": "warning"},
    )
    exit()

if "cluster" not in adata_g.obs:
    w_text_output(
        content="Neighborhood analysis requires `cluster` in `adata.obs`.",
        appearance={"message_box": "warning"},
    )
    exit()

if "filtered_groups" not in globals():
    filtered_groups = {}


def ordered_enrichment_matrix(adata_source, uns_key, mode):
    """Return z-score and selected-mode matrices with the same label order."""
    zscore, labels = _ordered_neighborhood_data(
        adata_source, uns_key, "zscore", "cluster"
    )
    selected, selected_labels = _ordered_neighborhood_data(
        adata_source, uns_key, mode, "cluster"
    )
    if labels != selected_labels:
        raise ValueError("Neighborhood matrix cluster labels do not align.")
    return zscore, selected, labels


def build_radial_neighborhood_fig(
    zscore_matrix,
    mode_matrix,
    labels,
    title,
    threshold=0.0,
    show_self=False,
    radius_max=None,
    zscore_max=None,
):
    """Build a grid of fixed-radius polar spoke plots, one per focal cluster."""
    del mode_matrix  # The selected metric is reported in the table below.
    n_clusters = len(labels)
    if n_clusters == 0:
        raise ValueError("No cluster categories are available to plot.")

    if radius_max is None or not np.isfinite(radius_max) or radius_max <= 0:
        radius_max = 1.0

    if zscore_max is None:
        off_diagonal = np.abs(zscore_matrix).copy()
        np.fill_diagonal(off_diagonal, 0.0)
        finite = off_diagonal[np.isfinite(off_diagonal)]
        zscore_max = float(np.nanmax(finite)) if finite.size else 1.0
    if not np.isfinite(zscore_max) or zscore_max <= 0:
        zscore_max = 1.0

    min_width, max_width = 0.04, 0.55
    num_cols = min(4, n_clusters)
    num_rows = math.ceil(n_clusters / num_cols)
    specs = [
        [{"type": "polar"} for _ in range(num_cols)]
        for _ in range(num_rows)
    ]
    fig = make_subplots(
        rows=num_rows,
        cols=num_cols,
        specs=specs,
        subplot_titles=[f"C{label}" for label in labels],
        horizontal_spacing=0.04,
        vertical_spacing=0.11,
    )

    for i, focal_cluster in enumerate(labels):
        row, col = i // num_cols + 1, i % num_cols + 1
        radii, angles, colors, widths, zscores = [], [], [], [], []
        for j, neighbor_cluster in enumerate(labels):
            if not show_self and neighbor_cluster == focal_cluster:
                continue
            zscore = float(zscore_matrix[i, j])
            if not np.isfinite(zscore) or abs(zscore) < threshold:
                continue

            radii.append(radius_max)
            angles.append(str(neighbor_cluster))
            colors.append("#C33530" if zscore >= 0 else "#282E66")
            denominator = zscore_max - threshold if zscore_max > threshold else 1.0
            fraction = max(0.0, min(1.0, (abs(zscore) - threshold) / denominator))
            widths.append(min_width + fraction * (max_width - min_width))
            zscores.append(zscore)

        fig.add_trace(
            go.Barpolar(
                theta=angles,
                r=radii,
                marker=dict(color=colors, line=dict(width=0)),
                width=widths,
                customdata=zscores,
                showlegend=False,
                hovertemplate=(
                    f"focal C{focal_cluster}"
                    "<br>neighbor C%{theta}"
                    "<br>z-score %{customdata:.2f}<extra></extra>"
                ),
            ),
            row=row,
            col=col,
        )

    fig.update_polars(
        radialaxis=dict(
            range=[0, radius_max],
            showticklabels=False,
            ticks="",
            showline=False,
        ),
        angularaxis=dict(
            direction="clockwise",
            rotation=90,
            categoryorder="array",
            categoryarray=[str(label) for label in labels],
            tickfont=dict(size=7),
        ),
        bgcolor="white",
    )
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center", font=dict(size=16)),
        width=250 * num_cols,
        height=270 * num_rows,
        margin=dict(l=20, r=20, t=80, b=20),
    )
    for annotation in fig.layout.annotations:
        annotation.font.size = 12
        annotation.yshift = 10
    return fig


radial_precomputed = (
    precomputed_neighborhoods
    if "precomputed_neighborhoods" in globals()
    else load_precomputed_neighborhood_groups(adata_g)
)

radial_group_dict = (
    group_dict
    if "group_dict" in globals()
    else {
        group: list(pd.unique(adata_g.obs[group].dropna()))
        for group in get_groupable_obs_keys(adata_g)
        if group != "cluster" and group not in na_keys
    }
)

radial_group_by = w_select(
    label="subplot groups",
    key="radial_group_by",
    default="all",
    options=tuple(["all"] + list(radial_group_dict)),
    appearance={
        "detail": "(all, categorical observation)",
        "help_text": "Facet radial plots by a categorical observation.",
    },
)

radial_mode = w_select(
    label="value metric",
    key="radial_mode",
    default="zscore",
    options=("zscore", "count"),
    appearance={
        "help_text": "Metric reported in the table; spoke color and thickness always use z-score."
    },
)

radial_threshold = w_text_input(
    label="z-score threshold",
    key="radial_threshold",
    default="0",
    appearance={
        "help_text": "Only connections at or above this absolute z-score are drawn."
    },
)

radial_show_self = w_checkbox(
    label="Show self-enrichment spoke",
    key="radial_show_self",
    default=False,
    appearance={"description": "Include the focal cluster's diagonal value."},
)

w_row(items=[radial_group_by, radial_mode, radial_threshold, radial_show_self])

if (
    radial_group_by.value not in (None, "all")
    and radial_group_by.value not in radial_precomputed
):
    w_text_output(
        content=(
            "This annotation was not precomputed by the workflow. Its spatial "
            "neighborhoods will be computed when first displayed and may take longer."
        ),
        appearance={"message_box": "warning"},
    )

radial_button = w_button(
    label="Update Radial Plots",
    key="radial_button",
)

if radial_group_by.value is not None and radial_button.value:
    try:
        threshold_value = (
            float(radial_threshold.value)
            if radial_threshold.value not in (None, "") else 0.0
        )
        if not np.isfinite(threshold_value) or threshold_value < 0:
            raise ValueError("z-score threshold must be a non-negative number.")

        mode_value = radial_mode.value
        show_self_value = bool(radial_show_self.value)
        radial_figures = []
        radial_table_rows = []

        if radial_group_by.value == "all":
            sample_key = "sample" if "sample" in adata_g.obs else None
            if "cluster_nhood_enrichment" not in adata_g.uns:
                w_text_output(
                    content="Computing neighborhoods for all cells...",
                    appearance={"message_box": "info"},
                )
                submit_widget_state()
                squidpy_analysis(adata_g, sample_key=sample_key)
            else:
                w_text_output(
                    content="Using workflow-precomputed neighborhood enrichment for all cells...",
                    appearance={"message_box": "info"},
                )
                submit_widget_state()

            zscore_matrix, mode_matrix, labels = ordered_enrichment_matrix(
                adata_g, "cluster_nhood_enrichment", mode_value
            )
            radial_figures.append(
                build_radial_neighborhood_fig(
                    zscore_matrix,
                    mode_matrix,
                    labels,
                    "All cells: Radial Neighborhood Enrichment",
                    threshold=threshold_value,
                    show_self=show_self_value,
                )
            )
            for i, focal_cluster in enumerate(labels):
                for j, neighbor_cluster in enumerate(labels):
                    radial_table_rows.append({
                        "subgroup": "all",
                        "focal_cluster": focal_cluster,
                        "neighbor_cluster": neighbor_cluster,
                        "zscore": zscore_matrix[i, j],
                        mode_value: mode_matrix[i, j],
                    })

        else:
            group = radial_group_by.value
            subgroup_values = radial_group_dict[group]
            group_adatas = filtered_groups.setdefault(group, {})
            for subgroup_key, subgroup_adata in radial_precomputed.get(group, {}).items():
                group_adatas.setdefault(subgroup_key, subgroup_adata)

            displayed_adatas = {}
            for subgroup in subgroup_values:
                subgroup_key = str(subgroup)
                if subgroup_key not in group_adatas:
                    group_adatas[subgroup_key] = make_lightweight_neighborhood_adata(
                        adata_g, group, subgroup
                    )
                subgroup_adata = group_adatas[subgroup_key]
                if "cluster_nhood_enrichment" not in subgroup_adata.uns:
                    w_text_output(
                        content=f"Computing spatial neighborhoods for {subgroup_key}...",
                        appearance={"message_box": "info"},
                    )
                    submit_widget_state()
                    sample_key = "sample" if "sample" in subgroup_adata.obs else None
                    squidpy_analysis(subgroup_adata, sample_key=sample_key)
                displayed_adatas[subgroup_key] = subgroup_adata

            matrices = {}
            shared_zscore_max = 0.0
            for subgroup_key, subgroup_adata in displayed_adatas.items():
                zscore_matrix, mode_matrix, labels = ordered_enrichment_matrix(
                    subgroup_adata, "cluster_nhood_enrichment", mode_value
                )
                matrices[subgroup_key] = (zscore_matrix, mode_matrix, labels)
                off_diagonal = np.abs(zscore_matrix).copy()
                np.fill_diagonal(off_diagonal, 0.0)
                finite = off_diagonal[np.isfinite(off_diagonal)]
                if finite.size:
                    shared_zscore_max = max(
                        shared_zscore_max, float(np.nanmax(finite))
                    )
            if shared_zscore_max <= 0:
                shared_zscore_max = 1.0

            for subgroup_key in sort_group_categories(list(matrices)):
                zscore_matrix, mode_matrix, labels = matrices[subgroup_key]
                radial_figures.append(
                    build_radial_neighborhood_fig(
                        zscore_matrix,
                        mode_matrix,
                        labels,
                        f"{group} = {subgroup_key}: Radial Neighborhood Enrichment",
                        threshold=threshold_value,
                        show_self=show_self_value,
                        zscore_max=shared_zscore_max,
                    )
                )
                for i, focal_cluster in enumerate(labels):
                    for j, neighbor_cluster in enumerate(labels):
                        radial_table_rows.append({
                            "subgroup": subgroup_key,
                            "focal_cluster": focal_cluster,
                            "neighbor_cluster": neighbor_cluster,
                            "zscore": zscore_matrix[i, j],
                            mode_value: mode_matrix[i, j],
                        })

    except Exception as error:
        w_text_output(
            content=f"Unable to create radial neighborhood plots: {error}",
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()

    # Index the list directly so every widget has a distinct source expression.
    # Passing each entry through one loop variable makes Latch treat every plot
    # as the same live global and causes later subgroup plots to overwrite it.
    for figure_index in range(len(radial_figures)):
        w_plot(
            source=radial_figures[figure_index],
            key=f"radial_plot_{figure_index}",
        )

    radial_enrichment_df = pd.DataFrame(radial_table_rows)
    w_table(
        label="Radial neighborhood enrichment values",
        source=radial_enrichment_df,
        key="radial_enrichment_table",
    )
