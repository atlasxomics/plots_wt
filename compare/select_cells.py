new_data_signal()

w_text_output(content="""
  # Compare Workflow
  
  Use the annotations from the AnnData object above to define cell groups for pairwise comparison with [ArchR::getMarkerFeatures](https://www.archrproject.com/reference/getMarkerFeatures.html).  
  Custom annotations can be added with the lasso-select tool or filters; see the instructions in the **H5 Viewer** section.
  
  <details>
  <summary><i>Compare Clusters Instructions</i></summary>
  
  1. In the **Select Annotation** input below, choose a **categorical** annotation from the AnnData object.  
     > If you do not see a recently added custom annotation, refresh the data by selecting *any* annotation value.  
     > Clicking the drop-down menu again will display your custom annotation.
  
  2. Select annotation values for **Group A** and **Group B**.
  
  3. Click **Confirm Cells** to ensure the selections are sufficient for the workflow.
  
  4. In the **Set Inputs for Compare Workflow** section:  
     - Enter a project name in the **Output Directory Name** field.  
     - Select the reference genome.
  
  5. Once the inputs have been correctly entered, click **Launch Workflow**.  
     This will trigger a workflow execution in your Latch workspace.  
     - You can monitor the execution status directly from this notebook.
  
  6. When the execution status is **Success**, click **Fetch Workflow Results** to load the outputs into interactive plots for data exploration.  
     - The output files are saved to **Latch Data** in the folder `compare_outs`.  
     - If the execution returns a status **Failed** or **Aborted**, please contact an ATX Support Scientist.
  
  </details>
""")


w_text_output(content="""## Select Annotation for Compare Workflow""")

if not adata_g:
  w_text_output(
    content="No data loaded...",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()

if not archrproj_dir:
  w_text_output(
    content="No ArchRProject found in input directory; directory MUST contain an ArchRProject for this module.",
    appearance={"message_box": "warning"}
  )
  submit_widget_state()
  exit()

choose_obs = w_select(
      label="Select Annotation",
      default=None,
      options=adata_g.obs_keys(),
      appearance={
        "help_text": "Selection with Annotation (.obs column) to define groups by."
      }
    )

if choose_obs.value is not None:

  obs_values = adata_g.obs[choose_obs.value]
  is_continuous = not (
    pd.api.types.is_categorical_dtype(obs_values)
    or obs_values.dtype.name == "category"
    or obs_values.nunique() < 30
  )
  if is_continuous:
    w_text_output(content="Continous data detected in Annotation; please select a Categorical Annotation.", appearance={"message_box": "warning"})
    submit_widget_state()
    groupselect_signal(False)
    barcodes_signal(False)
    exit()

  ann_keys = [ob for ob in obs_values.unique() if not pd.isna(ob)]

  if len(ann_keys) < 2:
    w_text_output(content="Fewer than two unique values detected in Annotation; please select an Annotation with two or greater values.", appearance={"message_box": "warning"})
    submit_widget_state()
    groupselect_signal(False)
    barcodes_signal(False)
    exit()

  groupA_ann = w_select(
    label="Group A Value",
    default=None,
    options=ann_keys,
    appearance={
      "help_text": "Select value for Group A."
    }
  )

  groupB_ann = w_select(
    label="Group B Value",
    default=None,
    options=ann_keys,
    appearance={
      "help_text": "Select value for Group B."
    }
  )

  groups_row = w_row(items=[groupA_ann, groupB_ann])

  if groupA_ann.value is not None and groupB_ann.value is not None:
    groupA_val = groupA_ann._signal.sample()
    groupB_val = groupB_ann._signal.sample()
    choose_group_signal(True)
