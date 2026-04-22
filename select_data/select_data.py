w_text_output(content="""
## Select Data
<details>
<summary><i>Instructions</i></summary>

**Loading data**  
- Click the **Select File** icon and choose a directory containing AnnData objects from the Latch Data module.  
- The directory should contain at least one of the following files:  
  - `adata_ge.h5ad`: a SnapATAC2 `AnnData` object with `.X` as a gene accessibility matrix.  
  - `adata_motifs.h5ad`: a SnapATAC2 `AnnData` object with `.X` as a motif deviation matrix.  
- BigWig files for cluster-, sample-, and condition-level groups should be saved in the output directory under subfolders named `[group]_coverages`.
- Loading large datasets into memory may take several minutes.  
- By default, compatible files are located in `latch:///snap_outs/[project_name]/`.
- If the notebook becomes frozen, try refreshing the browser tab or clicking the Run All In Tab b button in the Run All dropdown menu.

</details>
""")

# Select input data -----------------------------------------------------------

data_path = w_ldata_picker(
  label="atx_snap output folder",
  key="data_path",
  appearance={
    "placeholder": "Placeholder…"
  }
)

# Get .h5ad files -------------------------------------------------------------

if data_path.value is not None:

  if not data_path.value.is_dir():
      w_text_output(
          content="Selected resource must be a directory...",
          appearance={"message_box": "danger"}
      )
      submit_widget_state()
      exit()

  adata_g_paths = [f for f in data_path.value.iterdir() if "sm_ge.h5ad" in f.name()]
  adata_m_paths = [f for f in data_path.value.iterdir() if "sm_motifs.h5ad" in f.name()]

  if len(adata_g_paths) == 1:
      adata_g_path = adata_g_paths[0]
  elif len(adata_g_paths) == 0:
      adata_g_path = None
      adata_g = None
      w_text_output(
          content="No file with suffix 'sm_ge.h5ad' (gene data) found in \
            selected folder; selected folder MUST contain a \
            file ending in '_ge.h5ad'",
          appearance={"message_box": "danger"}
      )
      submit_widget_state()
      exit()
  elif len(adata_g_paths) > 1:
      adata_g_path = None
      adata_g = None
      w_text_output(
          content="Multiple files with suffix 'sm_ge.h5ad' (gene data) found \
            in selected folder; please ensure the output folder contains only \
            one file ending in '_ge.h5ad'",
          appearance={"message_box": "danger"}
      )
      submit_widget_state()
      exit()

  if len(adata_m_paths) == 1:
      adata_m_path = adata_m_paths[0]
  elif len(adata_m_paths) == 0:
      adata_m_path = None
      adata_m = None
      w_text_output(
          content="No file with suffix 'sm_motifs.h5ad' (motif data) found in \
            selected folder; selcted folder MUST contain a file \
            ending in '_motifs.h5ad'",
          appearance={"message_box": "danger"}
      )
      submit_widget_state()
      exit()
  elif len(adata_m_paths) > 1:
      adata_m_path = None  
      adata_m = None
      w_text_output(
          content="Multiple files with suffix 'sm_motifs.h5ad' (motif data) \
            found in selected folder; please ensure the output folder \
            contains only one file ending in '_motifs.h5ad'",
          appearance={"message_box": "danger"}
      )
      submit_widget_state()
      exit()
  
  if adata_g_path is None or adata_m_path is None:
      exit()
  
  # Download files ------------------------------------------------------------
  
  w_text_output(
    content="Downloading files and reading files; this may take a few minutes...",
    appearance={"message_box": "info"}
  )
  submit_widget_state()

  for path in [adata_g_path, adata_m_path]:
    path.download(Path(path.name()), cache=True)

  # Load files ----------------------------------------------------------------

  try:
    adata_g = sc.read(Path(adata_g_path.name()))
  except Exception as e:
    w_text_output(
      content=f"Error loading gene data: {e}\nPlease check input files.",
      appearance={"message_box": "danger"}
    )
    submit_widget_state()
    exit()
  
  available_genes = list(adata_g.var_names)

  # Ensure essential obs keys from ArchR
  adata_g = rename_obs_keys(adata_g)

  # Make obsm with all spatials offset
  if "spatial_offset" not in adata_g.obsm_keys():
      n_samples = adata_g.obs["sample"].nunique()
      n_cols = min(2, max(1, n_samples))
      n_rows = math.ceil(n_samples / n_cols)
      process_matrix_layout(adata_g, n_rows=n_rows, n_cols=n_cols, tile_spacing=300, new_obsm_key="spatial_offset")

  # Convert n_fragment to float for plotting
  if "n_fragment" in adata_g.obs_keys():
    adata_g.obs["n_fragment"] = adata_g.obs["n_fragment"].astype(float)

  try:
    adata_m = sc.read(Path(adata_m_path.name()))
  except Exception as e:
    w_text_output(
      content=f"Error loading motif data: {e}\nPlease check input files.",
      appearance={"message_box": "danger"}
    )
    submit_widget_state()
    exit()
    
  available_motifs = list(adata_m.var_names)

  # Ensure essential obs keys from ArchR
  adata_m = rename_obs_keys(adata_m)

  # Make obsm with all spatials offset
  if "spatial_offset" not in adata_m.obsm_keys():
    n_samples = adata_m.obs["sample"].nunique()
    n_cols = min(2, max(1, n_samples))
    n_rows = math.ceil(n_samples / n_cols)
    process_matrix_layout(adata_m, n_rows=n_rows, n_cols=n_cols, tile_spacing=300, new_obsm_key="spatial_offset")

  w_text_output(
    content=f"Data successfully loaded!",
    appearance={"message_box": "success"}
  )
  submit_widget_state()
  
  # Set default values --------------------------------------------------------

  samples = adata_g.obs["sample"].unique()
  groups = get_groups(adata_g)

  for data in [adata_g, adata_m]:
      for group in groups:
          if data.obs[group].dtype != object:  # Ensure groups are str
              data.obs[group] = data.obs[group].astype(str)

  available_metadata = tuple(key for key in adata_g.obs_keys()
                             if key not in na_keys)

  filtered_groups: dict[str, dict[str, anndata.AnnData]] = {}

  group_options = dict()
  for group in groups:
      group_options[group] = list(adata_g.obs[group].unique())

  clusters = group_options["cluster"]
  if "condition" in groups:
    conditions = group_options["condition"]

  # Reorder columns for H5 Viewer  ---------------------------------------------
  reorder_obs_columns(adata_g)
  reorder_obs_columns(adata_m) 

  drop_obs_column([adata_g, adata_m], col_to_drop="orig.ident") # remove orig.idents
  # Stuff for IGV  ------------------------------------------------------------

  coverages_dict = {}
  coverage_groups = groups if "sample" in groups else groups + ["sample"]
  for group in coverage_groups:
      for file in data_path.value.iterdir():
          if file.path.endswith(f"{group}_coverages"):
              coverages_dict[group] = file

  if len(coverages_dict) == 0:
      w_text_output(
          content="No coverage folders were found for project...",
          appearance={"message_box": "warning"}
      )
      submit_widget_state()

  # Stuff for Compare wf  ------------------------------------------------------

  # Get the current workspace account
  account = Account.current()
  account.load()
  workspace_account_id = account.id

  archrproj_dirs = [
    f for f in data_path.value.iterdir() if f.name().endswith("_ArchRProject")
  ]
  if len(archrproj_dirs) > 0:
    archrproj_dir = archrproj_dirs[0]
  else:
    w_text_output(
      content="No ArchRProject found for project...",
      appearance={"message_box": "warning"}
    )
    archrproj_dir = None

  genome_dict = {"hg38": Genome.hg38, "mm10": Genome.mm10}

  groupA_cells = []
  groupB_cells = []

  h5data_dict = {
    "gene": adata_g,
    "motif": adata_m
  }

  adata_h5 = None

  results_dict = {}
  feats = ["gene", "motif"]

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

  new_data_signal(True)
else:
  # Reset dynamic globals when no data path is selected.
  adata_g = None
  adata_m = None
  adata_g_path = None
  adata_m_path = None
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
  h5data_dict = {}
  adata_h5 = None
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
  refresh_h5_signal(False)

  new_data_signal(True)
