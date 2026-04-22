new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

groupselect_signal()
if groupselect_signal.sample() == True:
  confirm_cells = w_button(label="Confirm Cell Selections")
  
  if confirm_cells.value:
  
    try:
     wf_name._signal(None)
    except NameError:
      pass
    
    if len(groupA_cells) > 0 and len(groupB_cells) > 0:
      shared = set(groupA_cells) & set(groupB_cells)
      if shared:
        w_text_output(
          content=f"Selections share {len(shared)} cells; please confirm selections do not overlap.",
          appearance={"message_box": "danger"}
        )
        submit_widget_state()
        barcodes_signal(False)
      else:
        remote_bcs = None
        try:
          cells_payload = {"groupA": groupA_cells, "groupB": groupB_cells}
          with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
              json.dump(cells_payload, f)
              local_cfg = f.name
          dataset_name = Path(data_path.value.path.rstrip("/")).name or "dataset"
          dataset_slug = re.sub(r"[^A-Za-z0-9._-]+", "_", dataset_name).strip("._-")
          if len(dataset_slug) == 0:
            dataset_slug = "dataset"
          timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
          remote_bcs = LPath(
            f"{data_path.value.path}/compare_config_{dataset_slug}_A{len(groupA_cells)}_B{len(groupB_cells)}_{timestamp}.json"
          )
          remote_bcs.upload_from(Path(local_cfg))
        except:
          bcs_fail = w_text_output(
            content="Failed to upload barcodes to remote.",
            appearance={"message_box": "danger"}
          )
        
        barcodes_signal(True)
        bcs_success = w_text_output(
          content=(
            "Selections are ready for the Workflow! "
            f"Config saved to `{remote_bcs.path}`."
          ),
          appearance={"message_box": "success"}
        )
        submit_widget_state()
  
    else:
        w_text_output(
          content=f"Group A contains {len(groupA_cells)} cells, Group B contains {len(groupB_cells)}; please ensure selections for both groups are greater than 0.",
          appearance={"message_box": "danger"}
        )
        submit_widget_state()
        barcodes_signal(False)
