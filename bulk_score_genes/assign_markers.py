w_text_output(content="# Gene Scoring and Cell Type Assignment")

if "gene_score_done_signal" not in globals():
    gene_score_done_signal = Signal(False)
if "choose_subset_signal" not in globals():
    choose_subset_signal = Signal(False)
if "bulk_score_source_key" not in globals():
    bulk_score_source_key = None

new_data_signal()
current_bulk_score_source_key = None
if "adata_path" in globals() and adata_path is not None:
    current_bulk_score_source_key = adata_path.path

if current_bulk_score_source_key != bulk_score_source_key:
    bulk_score_source_key = current_bulk_score_source_key
    gene_score_done_signal(False)
    choose_subset_signal(False)
# Ensure AnnData is loaded
if adata_g is None or not isinstance(adata_g, AnnData):
    w_text_output(
        content="No AnnData object loaded...",
        appearance={"message_box": "warning"}
    )
    exit()

w_text_output(content="""
Use this tab to score cells with the average expression of a set of genes (see [scanpy.tl.score_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html)).
These scores can be used to guide the assignment of putative cell types with the H5 Viewer at the bottom of the tab. 
"""
)

# --- Step 1: Text input for number of cell types ---
w_text_output(content="""
## Select Number of Cell Types
Specify how many cell types you want to predict, along with the marker gene sets for each. For each cell, the tool calculates a gene set score for every cell type and assigns the cell
type with the highest score as the _predicted cell type_, stored in the Categorical Observation column **pred_cell_type**.
""")

w_text_output(content="""
Gene scores can vary across conditions. Subseting your dataset to a single condition before computing scores could improve the data clarity.
"""
)
# Checkbox to show/hide layout controls
subset_button = w_checkbox(
    label="Filter dataset",
    default=False,
    appearance={"description": "Toggle to turn on/off subsetting."}
)

if subset_button.value:
  w_text_output(content="""
  **Optionally**, use this section to subset your data.  Start by selecting the metadata column that
  defines your groups, then choose the specific group you'd like to focus on.
  """)

  filter_groups = [
    key for key in get_groupable_obs_keys(adata_g)
    if key != "cluster"
  ]
  filter_groups = filter_groups + ["None"]
  
  filter_col = w_select(
      label="Filter by metadata column",
      default=None,
      options=tuple(filter_groups),
      key="filter_col",
      appearance={"help_text": "Select an categorical observation to filter cells."}
  )
  
  use_filter = False if filter_col.value == None or filter_col.value == "None" else True

  # Widget to choose the value within that column
  if use_filter:
      vals = adata_g.obs[filter_col.value].dropna().unique().tolist()
      filter_val = w_select(
          label=f"Filter value for {filter_col.value}",
          default=None,
          options=tuple(vals),
          key="filter_val",
          appearance={"help_text": f"Select a value in {filter_col.value} to subset cells."}
      )
  else:
      filter_val = None

  # Require a filter value if a column is chosen
  if use_filter and filter_val.value is None:
      w_text_output(
          content="Please select a filter value to subset cells.",
          appearance={"message_box": "warning"}
      )
      submit_widget_state()
      exit()
  
# Create the subset AnnData
if subset_button.value and use_filter and filter_val.value is not None:
    adata_subset = adata_g[adata_g.obs[filter_col.value] == filter_val.value].copy()
else:
    adata_subset = adata_g.copy()

n_types_input = w_text_input(
    label="Number of cell types",
    default="",
    key="no_cell_types_key",
    appearance={"help_text": "Enter an integer number of cell types"}
)

# --- Step 2: Parse the input into an integer ---
try:
    n_types = int(n_types_input.value)
except (ValueError, TypeError):
    n_types = 0

if n_types == 1:
  w_text_output(
    content="Only one cell types specified; all cells will be assigned to this cell type...",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()

# Prepare lists to store widgets
label_inputs = []
feature_selects = []

# --- Step 3: Dynamically create widgets for each cell type ---
for idx in range(n_types):
    label_input = w_text_input(
        label=f"Cell type {idx+1} label",
        default="",
        appearance={"help_text": "Enter a name for this cell type"},
        key=f"Cell type {idx+1} label"
    )
    feature_select = w_multi_select(
        label=f"Select features for cell type {idx+1}",
        options=tuple(adata_g.var_names),
        appearance={"help_text": "Choose one or more features for this cell type"},
        key=f"Select features for cell type {idx+1}"
    )
    label_inputs.append(label_input)
    feature_selects.append(feature_select)
    w_row(items=[label_input, feature_select])

confirm_celltypes = w_button(label="Confirm Inputs")

if confirm_celltypes.value:
  
  # --- Validation and gating for downstream steps ---
  validation_errors = []
  
  if n_types <= 0:
      validation_errors.append("Number of cell types must be greater than 0.")
  
  for idx in range(n_types):
      label_value = (label_inputs[idx].value or "").strip()
      if not label_value:
          validation_errors.append(f"Cell type {idx + 1} label cannot be empty.")

      features = feature_selects[idx].value
      if not features:
          validation_errors.append(f"Select at least one feature for cell type {idx + 1}.")

  if validation_errors:
        w_text_output(
          content="Please resolve the following before continuing:\n" + "\n".join(
              f"- {msg}" for msg in validation_errors
          ),
          appearance={"message_box": "warning"}
        )
        choose_subset_signal(False)
  else:
      choose_subset_signal(True)
      submit_widget_state()
