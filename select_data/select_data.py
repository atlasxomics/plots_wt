w_text_output(content="""
## Select Data
<details>
<summary><i>Instructions</i></summary>

**Loading data**
- Click the **Select File** icon and choose a `.h5ad` file from Latch Data.
- The file should contain a valid AnnData object.
- The notebook will use `sample`, `condition`, and `cluster` from `.obs` when present, but these columns are optional.
- Loading large datasets into memory may take several minutes.
- If the notebook becomes frozen, try refreshing the browser tab or clicking the Run All In Tab button in the Run All dropdown menu.

</details>
""")

# Select input data -----------------------------------------------------------

data_path = w_ldata_picker(
  label="AnnData H5AD file",
  key="data_path",
  appearance={
    "placeholder": "Select a .h5ad file"
  }
)

# Load .h5ad file -------------------------------------------------------------

if data_path.value is not None:

  if data_path.value.is_dir():
      w_text_output(
          content="Selected resource must be a `.h5ad` file, not a directory.",
          appearance={"message_box": "danger"}
      )
      submit_widget_state()
      exit()

  adata_path = data_path.value

  adata_filename = adata_path.name()
  if not adata_filename.endswith(".h5ad"):
      w_text_output(
          content="Selected file must end with `.h5ad`.",
          appearance={"message_box": "danger"}
      )
      submit_widget_state()
      exit()

  local_adata_path = Path(adata_filename)

  # Download file -------------------------------------------------------------

  w_text_output(
    content="Downloading and reading the AnnData file; this may take a few minutes...",
    appearance={"message_box": "info"}
  )
  submit_widget_state()

  try:
    adata_path.download(local_adata_path, cache=True)
  except Exception as e:
    w_text_output(
      content=f"Error downloading input file: {e}",
      appearance={"message_box": "danger"}
    )
    submit_widget_state()
    exit()

  # Load file -----------------------------------------------------------------

  try:
    adata = sc.read_h5ad(local_adata_path)
  except Exception as e:
    w_text_output(
      content=f"Error loading AnnData object: {e}\nPlease check input file.",
      appearance={"message_box": "danger"}
    )
    submit_widget_state()
    exit()

  if not isinstance(adata, AnnData):
    w_text_output(
      content="Selected file did not load as an AnnData object.",
      appearance={"message_box": "danger"}
    )
    submit_widget_state()
    exit()

  # Compatibility alias for tabs that still refer to the primary object as adata_g.
  adata_g = adata

  # Normalize common alternate metadata names without requiring them.
  adata = rename_obs_keys(adata)
  adata_g = adata

  for col in ["n_fragment", "n_counts", "total_counts"]:
    if col in adata.obs_keys():
      adata.obs[col] = pd.to_numeric(adata.obs[col], errors="ignore")

  for group in get_groups(adata):
      if adata.obs[group].dtype != object:
          adata.obs[group] = adata.obs[group].astype(str)

  available_features = list(adata.var_names)

  missing_recommended = [
    col for col in ["sample", "condition", "cluster"]
    if col not in adata.obs
  ]
  if missing_recommended:
    w_text_output(
      content=(
        "Loaded without optional obs column(s): "
        + ", ".join(missing_recommended)
        + ". Some default groupings will use other available metadata columns."
      ),
      appearance={"message_box": "warning"}
    )
    submit_widget_state()

  if "sample" in adata.obs and "spatial" in adata.obsm_keys() and "X_dataset" not in adata.obsm_keys():
      n_samples = adata.obs["sample"].nunique()
      if n_samples > 0:
        n_cols = min(2, max(1, n_samples))
        n_rows = math.ceil(n_samples / n_cols)
        try:
          process_matrix_layout(
            adata,
            n_rows=n_rows,
            n_cols=n_cols,
            tile_spacing=300,
            new_obsm_key="X_dataset",
          )
        except Exception as e:
          w_text_output(
            content=f"Could not create combined sample spatial layout: {e}",
            appearance={"message_box": "warning"}
          )
          submit_widget_state()

  reorder_obs_columns(adata)
  drop_obs_column([adata], col_to_drop="orig.ident")

  w_text_output(
    content=f"Data successfully loaded: {adata.n_obs} observations x {adata.n_vars} features.",
    appearance={"message_box": "success"}
  )
  submit_widget_state()

  choose_subset_signal(False)
  gene_score_done_signal(False)
  refresh_h5_signal(False)

  new_data_signal(True)
else:
  # Reset dynamic globals when no data path is selected.
  adata = None
  adata_g = None
  adata_path = None
  available_features = []
  choose_subset_signal(False)
  gene_score_done_signal(False)
  refresh_h5_signal(False)

  new_data_signal(True)
