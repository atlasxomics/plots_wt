new_data_signal()

w_text_output(content="""

# Neighborhood Analysis

Explore spatial neighborhood enrichment among clusters, either for **all cells** or split by a subgroup (e.g., **sample** or **condition**).  See squipy [neighbors enrichment analysis](https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_nhood_enrichment.html) for more information.


<details>
<summary><i>details</i></summary>

Each heatmap cell reflects how often cells from **cluster A** neighbor cells from **cluster B** compared with chance.  You can view values as **z-scores** (standardized enrichment; recommended) or **counts** (raw neighborhood counts).  Optionally, you can facet plots by any categorical observation to compare neighborhood structure across custom annotations.

### Controls

1. **subplot groups** 
   - Options: **all** plus any categorical observation (for example **sample**, **condition**, or custom H5 Viewer annotations)  
   - **all**: one heatmap using all cells.
   - **categorical observation**: one heatmap per subgroup (faceted).

2. **displayed data**
   - Options: **zscore**, **count**  
   - **zscore**: standardized neighborhood enrichment (best for comparisons). 
   - **count**: raw neighbor counts (scale depends on dataset size).

3. **hierarchical clustering method** 
   - Options: **None**, **single**, **complete**, **average**, **weighted**, **centroid**, **median**, **ward**  
   - Choose **None** to keep the original cluster order, or a method to cluster rows/columns and group similar patterns.

4. **colorscale maximum / minimum** 
   - Optional numeric limits for the heatmap color range (e.g., max = `5`, min = `-2`).  
   - Leave blank to auto-scale.
</details>

""")

if not adata_g:
    w_text_output(
        content="No data gene activity data selected...",
        appearance={"message_box": "warning"}
    )
    exit()

notebook_palettes = await get_notebook_palettes()

neighbor_groups = [
  key for key in adata_g.obs_keys()
  if key != "cluster" and (
    pd.api.types.is_object_dtype(adata_g.obs[key]) or
    pd.api.types.is_categorical_dtype(adata_g.obs[key])
  )
]
group_dict = {g: adata_g.obs[g].dropna().unique() for g in neighbor_groups}

neigh_group_by = w_select(
  label="subplot groups",
  default="all",
  options=tuple(["all"] + list(group_dict.keys())),
  appearance={
    "detail": "(all, categorical observation)",
    "help_text": "Facet neighborhood plots by any categorical observation."
  }
)

mode = w_select(
  label="displayed data",
  default="zscore",
  options=("zscore", "count"),
  appearance={
    "help_text": "Data to be plotted"
  }
)

scale_max = w_text_input(
  label="colorscale maximum",
  default=None,
  appearance={
    "help_text": "Maximum value of colorscale"
  }
)

scale_min = w_text_input(
  label="colorscale minimum",
  default=None,
  appearance={
  "help_text": "Minimum value of colorscale"
  }
)

neigh_palette = w_select(
  label="colorscale palette",
  default="Default Neighborhood Colorscale",
  options=get_palette_selector_options(
    notebook_palettes,
    kind="continuous",
    fallback_name="Default Neighborhood Colorscale",
  ),
  appearance={
    "help_text": "Use a continuous palette saved from the H5 Viewer or fall back to the default neighborhood colors."
  }
)

w_row(items=[neigh_group_by, mode, scale_max, scale_min, neigh_palette])

neigh_button = w_button(label="Update Neighborhood Plots")

if neigh_group_by.value is not None and neigh_button.value:
  neigh_colorscale = get_selected_continuous_palette(
    notebook_palettes,
    neigh_palette.value,
    fallback_colors="RdBu_r",
    fallback_name="Default Neighborhood Colorscale",
  )
  
  vmax = int(scale_max.value) if scale_max.value and scale_max.value.strip().isdigit() else None
  try:  # Handle negative values
    vmin = int(scale_min.value) if scale_min.value.strip() != '' else None
  except ValueError:
    vmin = None  # Fallback if the value can't be converted to an integer
  
  # --------------------------------------------------------------------------------
  
  if neigh_group_by.value == "all":
    sample_key = "sample" if "sample" in groups else None
    if "cluster_nhood_enrichment" not in adata_g.uns:
      w_text_output(
        content="Computing neighborhoods for all cells...",
        appearance={"message_box": "info"}
      )
      submit_widget_state()
      squidpy_analysis(adata_g, sample_key=sample_key)
    else:
      w_text_output(
        content="Using existing neighborhood enrichment for all cells...",
        appearance={"message_box": "info"}
      )
      submit_widget_state()

    neigh_heatmap, neigh_data = plotly_heatmap(
      adata_g,
      uns_key="cluster_nhood_enrichment",
      title=f"{neigh_group_by.value} cells: Neighborhood Enrichment",
      mode=mode.value,
      colorscale=normalize_plotly_colorscale(neigh_colorscale),
      vmax=vmax,
      vmin=vmin
    )

    neigh_data = pd.DataFrame(neigh_data)
  
  
  elif neigh_group_by.value in group_dict:
  
    group = neigh_group_by.value
    sub_groups = group_dict[group]
    if group not in filtered_groups:
      filtered_adatas: dict[str, anndata.AnnData] = {}

      filtered_groups[group] = filtered_adatas

    filtered_adatas = filtered_groups[group]

    for sg in sub_groups:
      if sg not in filtered_adatas:
        filtered_adatas[sg] = filter_anndata(adata_g, group, sg)

      filtered_adata = filtered_adatas[sg]
      if "cluster_nhood_enrichment" not in filtered_adata.uns:
        w_text_output(
          content=f"Computing spatial neighborhoods for {sg}...",
          appearance={"message_box": "info"}
        )
        submit_widget_state()
        sample_key = "sample" if "sample" in groups else None
        squidpy_analysis(filtered_adata, sample_key=sample_key)
      else:
        w_text_output(
          content=f"Using existing neighborhood enrichment for {sg}...",
          appearance={"message_box": "info"}
        )
        submit_widget_state()

    neigh_heatmap, neigh_data = plot_neighborhood_groups(
      filtered_groups[group],
      f"Neighborhoods by {group}",
      uns_key="cluster_nhood_enrichment",
      mode=mode.value,
      colorscale=neigh_colorscale,
      vmax=vmax,
      vmin=vmin
    )

  else:
    raise KeyError("Group by not expected value")
  
  w_plot(source=neigh_heatmap)
  w_table(source=neigh_data)
