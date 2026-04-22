w_text_output(content="## H5 Viewer")

w_text_output(content="""
<details>
<summary><i>H5 Viewer Instructions</i></summary>

Once initialized, the viewer will display your project cells plotted in UMAP space.

## Navigating the Viewer
- **Change coordinates:**
  - `X_umap`: UMAP coordinates
  - `X_dataset`: samples arranged in spatial layout
  - `spatial`: *not recommended* (overlays all samples at once)
- **Coloring cells:**
  - In the left panel, click the **paint bucket** icon next to any annotation to color by that variable.
  - In the *Genes of Interest* section (bottom of left panel), click **+**, type a gene name, select it, and click **Add**. Then click the paint bucket icon next to the feature to color cells by that gene or motif.

---

## Selecting Cells for the Compare Clusters Workflow

The H5 Viewer can be used to create custom annotations in the AnnData object. These annotations can then be used as input to the **Compare Clusters Workflow**.

- **Add a new annotation:**  
  In the left panel, click the **+** next to either *Continuous* or *Categorical* observations.
- **Select cells with the lasso tool:**
  - Hover over the AnnData figure, click the **lasso icon** in the Plotly toolbar (top-right).
  - Use your cursor to manually outline cells of interest.
  - To select multiple regions, **hold Shift**.
- Once cells are lasso-selected, a new **+** will appear next to the observation fields. Click it to add your selected cells to the annotation.

---

## Using Filters to Define New Observations

- **Filter cells:**  
  Click the **ellipsis (…)** next to an observation → select **Filter**, then enter your criteria.  
- **Capture filtered cells:** 
  Use the lasso tool to select the filtered cells and add them to a custom annotation as described above.  
- **Clear filters:**  
  Click the **filter icon** in the bottom-right of the display panel to remove active filters.  
- Repeat this process to create additional filter-defined groups in your AnnData object.

</details>
""")

new_data_signal()
if not adata_g:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    exit()

if "adata_h5" not in globals() or adata_h5 is None:
    adata_h5 = adata_g
    loaded_h5_data_key = "gene"
elif adata_h5 is not adata_g and adata_h5 is not adata_m:
    adata_h5 = adata_g
    loaded_h5_data_key = "gene"
elif "loaded_h5_data_key" not in globals():
    if adata_h5 is adata_m:
        loaded_h5_data_key = "motif"
    else:
        loaded_h5_data_key = "gene"

refresh_h5_signal()
viewer = w_h5(ann_data=adata_h5)

h5_obs_button = w_checkbox(
    label="Display Cell Metadata Table",
    key="h5_obs_button",
    default=False,
)

h5_obs = adata_h5.obs
if h5_obs_button.value:
  h5_table = w_table(label="Metadata (adata.obs)", source=h5_obs)

h5_viewer_signal(True)
