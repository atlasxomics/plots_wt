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
  adata_g_path = adata_path
  adata_m_path = None

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
  adata_m = None

  # Normalize common alternate metadata names without requiring them.
  adata = rename_obs_keys(adata)
  adata_g = adata

  for col in ["n_fragment", "n_counts", "total_counts"]:
    if col in adata.obs_keys():
      adata.obs[col] = pd.to_numeric(adata.obs[col], errors="ignore")

  groups = get_groups(adata)
  for group in groups:
      if adata.obs[group].dtype != object:
          adata.obs[group] = adata.obs[group].astype(str)

  available_features = list(adata.var_names)
  available_genes = available_features
  available_motifs = []
  available_metadata = tuple(key for key in adata.obs_keys()
                             if key not in na_keys)

  filtered_groups: dict[str, dict[str, anndata.AnnData]] = {}

  group_options = {
    group: list(adata.obs[group].dropna().unique())
    for group in groups
  }
  samples = adata.obs["sample"].dropna().unique() if "sample" in adata.obs else []
  clusters = group_options.get("cluster", [])
  conditions = group_options.get("condition", [])

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

  coverages_dict = {}
  archrproj_dir = None
  workspace_account_id = None
  genome_dict = {"hg38": Genome.hg38, "mm10": Genome.mm10}

  groupA_cells = []
  groupB_cells = []

  h5data_dict = {
    "adata": adata
  }

  adata_h5 = adata
  loaded_h5_data_key = "adata"

  results_dict = {}
  feats = ["feature"]

  w_text_output(
    content=f"Data successfully loaded: {adata.n_obs} observations x {adata.n_vars} features.",
    appearance={"message_box": "success"}
  )
  submit_widget_state()

  choose_group_signal(False)
  groupselect_signal(False)
  barcodes_signal(False)
  wf_ready_signal(False)
  wf_exe_signal(False)
  wf_results_signal(False)
  wf_bigwigs_signal(False)

  # Other signals ------------------------------------------------------
  h5_viewer_signal(False)
  compare_signal(False)
  heatmap_signal(False)
  tracks_signal(False)
  choose_subset_signal(False)
  gene_score_done_signal(False)
  refresh_h5_signal(False)

  new_data_signal(True)
else:
  # Reset dynamic globals when no data path is selected.
  adata = None
  adata_g = None
  adata_m = None
  adata_path = None
  adata_g_path = None
  adata_m_path = None
  local_adata_path = None
  available_features = []
  available_genes = []
  available_motifs = []
  samples = []
  groups = []
  group_options = {}
  clusters = []
  conditions = []
  available_metadata = ()
  coverages_dict = {}
  archrproj_dir = None
  workspace_account_id = None
  h5data_dict = {}
  adata_h5 = None
  loaded_h5_data_key = None
  results_dict = {}
  feats = []

  choose_group_signal(False)
  groupselect_signal(False)
  barcodes_signal(False)
  wf_ready_signal(False)
  wf_exe_signal(False)
  wf_results_signal(False)
  wf_bigwigs_signal(False)

  h5_viewer_signal(False)
  compare_signal(False)
  heatmap_signal(False)
  tracks_signal(False)
  choose_subset_signal(False)
  gene_score_done_signal(False)
  refresh_h5_signal(False)

  new_data_signal(True)
