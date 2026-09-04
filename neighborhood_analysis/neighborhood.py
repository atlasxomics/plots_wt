new_data_signal()

w_text_output(content="""

# Neighborhood Analysis

Explore spatial neighborhood enrichment among clusters, either for **all cells** or split by a subgroup such as **sample** or **condition**.

<details>
<summary><i>details</i></summary>

Each heatmap cell shows how often cells from **cluster A** neighbor cells from **cluster B** compared with chance. Values can be shown as standardized **z-scores** (recommended) or raw neighborhood **counts**.

Workflow-generated sample and condition results are loaded from the precomputed neighborhood matrices in the selected AnnData file. A categorical annotation added later can also be selected; its results are computed when first displayed and cached for this notebook session.

### Controls

1. **display type** switches between this heatmap and the radial spoke view below.
2. **subplot groups** displays all cells together or facets results by a categorical observation.
3. **displayed data** selects z-scores or counts.
4. **colorscale maximum / minimum** optionally fixes the shared heatmap range.

</details>

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


def load_precomputed_neighborhood_groups(adata):
    """Rebuild lightweight AnnData objects from workflow-precomputed results."""
    root = adata.uns.get("cluster_nhood_enrichment_by_group")
    if not isinstance(root, dict) or root.get("schema_version") != 1:
        return {}
    if str(root.get("cluster_key", "cluster")) != "cluster":
        return {}

    result = {}
    for group_entry in root.get("groups", {}).values():
        group_key = str(group_entry["group_key"])
        subgroups = {}
        for subgroup_entry in group_entry.get("subgroups", {}).values():
            group_value = str(subgroup_entry["group_value"])
            categories = pd.Index(
                np.asarray(subgroup_entry["cluster_categories"]).astype(str)
            )
            obs = pd.DataFrame(
                {
                    "cluster": pd.Categorical(
                        categories,
                        categories=categories,
                        ordered=True,
                    )
                },
                index=[f"cluster-{i}" for i in range(len(categories))],
            )
            neighborhood_adata = anndata.AnnData(obs=obs)
            neighborhood_adata.uns["cluster_nhood_enrichment"] = {
                "zscore": np.asarray(subgroup_entry["zscore"], dtype=float),
                "count": np.asarray(subgroup_entry["count"], dtype=float),
            }
            subgroups[group_value] = neighborhood_adata
        result[group_key] = subgroups
    return result


def make_lightweight_neighborhood_adata(adata, group, subgroup):
    """Subset metadata and coordinates without copying the feature matrix."""
    spatial_key = get_spatial_layout_key(adata)
    if spatial_key is None:
        raise KeyError(
            "Spatial coordinates were not found in `.obsm`; expected "
            "`spatial` or `spatial_offset`."
        )

    mask = (adata.obs[group] == subgroup).to_numpy()
    obs_keys = ["cluster"]
    if "sample" in adata.obs and group != "sample":
        obs_keys.append("sample")
    subset_obs = adata.obs.loc[mask, obs_keys].copy()
    for obs_key in obs_keys:
        if isinstance(subset_obs[obs_key].dtype, pd.CategoricalDtype):
            subset_obs[obs_key] = subset_obs[obs_key].cat.remove_unused_categories()

    lightweight_adata = anndata.AnnData(obs=subset_obs)
    lightweight_adata.obsm[spatial_key] = np.asarray(adata.obsm[spatial_key])[mask].copy()
    return lightweight_adata


def parse_optional_neighborhood_float(value, label):
    if value is None or str(value).strip() == "":
        return None
    try:
        return float(str(value).strip())
    except ValueError as error:
        raise ValueError(f"{label} must be numeric.") from error


precomputed_neighborhoods = load_precomputed_neighborhood_groups(adata_g)

neighbor_groups = [
    key for key in get_groupable_obs_keys(adata_g)
    if key != "cluster" and key not in na_keys
]
for key in precomputed_neighborhoods:
    if key in adata_g.obs and key not in neighbor_groups:
        neighbor_groups.append(key)

group_dict = {
    group: list(pd.unique(adata_g.obs[group].dropna()))
    for group in neighbor_groups
}

neigh_display = w_select(
    label="display type",
    key="neigh_display",
    default="heatmap",
    options=("heatmap", "radial"),
    appearance={
        "help_text": "Switch between heatmap and radial spoke views."
    },
)

if neigh_display.value != "heatmap":
    exit()

neigh_group_by = w_select(
    label="subplot groups",
    key="neigh_group_by",
    default="all",
    options=tuple(["all"] + neighbor_groups),
    appearance={
        "detail": "(all, categorical observation)",
        "help_text": "Facet neighborhood plots by a categorical observation.",
    },
)

mode = w_select(
    label="displayed data",
    key="neigh_mode",
    default="zscore",
    options=("zscore", "count"),
    appearance={"help_text": "Neighborhood values to plot."},
)

scale_max = w_text_input(
    label="colorscale maximum",
    key="neigh_scale_max",
    default=None,
    appearance={"help_text": "Optional maximum value of the shared colorscale."},
)

scale_min = w_text_input(
    label="colorscale minimum",
    key="neigh_scale_min",
    default=None,
    appearance={"help_text": "Optional minimum value of the shared colorscale."},
)

w_row(items=[neigh_group_by, mode, scale_max, scale_min])

if (
    neigh_group_by.value not in (None, "all")
    and neigh_group_by.value not in precomputed_neighborhoods
):
    w_text_output(
        content=(
            "This annotation was not precomputed by the workflow. Its spatial "
            "neighborhoods will be computed when first displayed and may take longer."
        ),
        appearance={"message_box": "warning"},
    )

neigh_button = w_button(
    label="Update Neighborhood Plots",
    key="neigh_button",
)

if neigh_group_by.value is not None and neigh_button.value:
    try:
        vmax = parse_optional_neighborhood_float(
            scale_max.value, "colorscale maximum"
        )
        vmin = parse_optional_neighborhood_float(
            scale_min.value, "colorscale minimum"
        )
        if vmin is not None and vmax is not None and vmin >= vmax:
            raise ValueError("colorscale minimum must be less than its maximum.")

        if neigh_group_by.value == "all":
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

            neigh_heatmap, neigh_data = plotly_heatmap(
                adata_g,
                uns_key="cluster_nhood_enrichment",
                title="All cells: Neighborhood Enrichment",
                mode=mode.value,
                vmax=vmax,
                vmin=vmin,
            )

        elif neigh_group_by.value in group_dict:
            group = neigh_group_by.value
            subgroup_values = group_dict[group]
            filtered_adatas = filtered_groups.setdefault(group, {})
            for subgroup_key, subgroup_adata in precomputed_neighborhoods.get(group, {}).items():
                filtered_adatas.setdefault(subgroup_key, subgroup_adata)

            displayed_adatas = {}
            for subgroup in subgroup_values:
                subgroup_key = str(subgroup)
                if subgroup_key not in filtered_adatas:
                    filtered_adatas[subgroup_key] = make_lightweight_neighborhood_adata(
                        adata_g, group, subgroup
                    )

                filtered_adata = filtered_adatas[subgroup_key]
                if "cluster_nhood_enrichment" not in filtered_adata.uns:
                    w_text_output(
                        content=f"Computing spatial neighborhoods for {subgroup_key}...",
                        appearance={"message_box": "info"},
                    )
                    submit_widget_state()
                    sample_key = "sample" if "sample" in filtered_adata.obs else None
                    squidpy_analysis(filtered_adata, sample_key=sample_key)
                else:
                    w_text_output(
                        content=f"Using existing neighborhood enrichment for {subgroup_key}...",
                        appearance={"message_box": "info"},
                    )
                    submit_widget_state()
                displayed_adatas[subgroup_key] = filtered_adata

            neigh_heatmap, neigh_data = plot_neighborhood_groups(
                displayed_adatas,
                f"Neighborhoods by {group}",
                uns_key="cluster_nhood_enrichment",
                mode=mode.value,
                vmax=vmax,
                vmin=vmin,
            )
        else:
            raise KeyError("Unexpected neighborhood grouping selection.")

    except Exception as error:
        w_text_output(
            content=f"Unable to create neighborhood plots: {error}",
            appearance={"message_box": "danger"},
        )
        submit_widget_state()
        exit()

    w_plot(source=neigh_heatmap, key="neigh_heatmap")
    w_table(source=neigh_data, key="neigh_table")
