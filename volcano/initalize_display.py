new_data_signal()

w_text_output(content="""
# Compare Conditions

Visualize pre-computed differential statistics for project **conditions** in either genes or motifs with a **volcano plot** and a **ranked feature plot**.

To initialize the figures, select the features (genes or motifs) you’d like to display from the drop-down menu below.

<details>
<summary><i>Instructions</i></summary>

### Steps

1. **Select Condition**
   - Use **Condition** to pick the condition you want to test (e.g., *treated* vs *untreated*).

2. **Optional: Filter by Cluster**
   - Use **Cluster** to limit the analysis to a single cluster, or choose **All** to use all cells.

3. **Set Significance Thresholds**
   - **pval adjust Threshold**: maximum adjusted p-value (e.g., `0.05`).
   - **Difference Metric Threshold**: minimum effect size. 
     - For **genes**: Log2FC 
     - For **motifs**: MeanDiff
   - Points below these thresholds are de-emphasized in the plots.

4. **Rank & Color Controls**
   - **Rank By**: metric used to order features in the ranked plot (e.g., `-log10(p_val_adj)`, `Log2FC`, `MeanDiff`).  
     - If you choose a p-value metric, the app automatically converts it to **−log10(p)** for readability.
   - **Color By**: which metric determines point color in the ranked plot.

5. **Review the Outputs**
   - **Volcano Plot**: interactive plot highlighting top features (labels) that pass your thresholds.
   - **Ranked Feature Plot**: interactive dot plot ordered by *Rank By*, colored by *Color By*.
   - **Results Table**: sortable table with all features and statistics used by the plots.

---

### Notes & Warnings
- **Imbalanced clusters**: If a cluster contains **>90%** of one condition, a volcano plot is not generated for that cluster.  
  - You’ll see: *“There is no volcano plot for cluster X because it contains more than 90% of one of the conditions. Please check Proportion plot.”*
- **All vs. single cluster**: Start with **All** to get a global view; drill down by selecting a specific cluster.
- **Performance**: Large datasets may take time to render, especially when re-ranking or re-coloring features.

""")

# Abort if no data loaded
if not adata_g:
    w_text_output(
        content="No data loaded…",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

try:
  conditions = group_options["condition"]
except KeyError:
  w_text_output(
    content="No conditions found in experiment...",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  compare_signal(False)
  exit()    

# Choose whether to display gene or motif data
choose_compare_data = w_select(
    label="Select Data for Comparison Plots",
    default="gene",
    options=["gene", "motif"],
    appearance={
        "help_text": "Select which features to display in the comparison plots."
    }
)

if choose_compare_data.value == "gene":
  rankby_opts = ['Log2FC', 'p_val_adj']
  rankby_default = 'Log2FC'
  diff_metric = 'Log2FC'
elif choose_compare_data.value == "motif":
  rankby_opts = ['p_val_adj', 'MeanDiff']
  rankby_default = 'MeanDiff'
  diff_metric = 'MeanDiff'

compare_signal(True)
