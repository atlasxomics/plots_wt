new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

choose_group_signal()
# Prevent user selecting the same group twice

# if "groupA_val" in globals() and "groupB_val" in globals():
if choose_group_signal.sample() == True:
  if groupA_val == groupB_val:
    w_text_output(
        content="Please ensure different values are selected for Group A and Group B.",
        appearance={"message_box": "warning"}
    )
    submit_widget_state()
    groupselect_signal(False)
    barcodes_signal(False)
    exit()

  # Valid selections: update signals and compute barcodes lists
  groupselect_signal(True)
  groupA_cells = list(adata_h5[adata_h5.obs[choose_obs.value] == groupA_val].obs_names)
  groupB_cells = list(adata_h5[adata_h5.obs[choose_obs.value] == groupB_val].obs_names)
  
  # Display collapsible lists of barcodes for each group
  groupA_text = w_text_output(
      content=f"<details><summary><i>Group A barcodes</i></summary>{','.join(groupA_cells)}</details>"
  )
  groupB_text = w_text_output(
      content=f"<details><summary><i>Group B barcodes</i></summary>{','.join(groupB_cells)}</details>"
  )
  
  with w_grid(columns=2) as g2:
      g2.add(item=groupA_text, col_span=1)
      g2.add(item=groupB_text, col_span=1)
  
  # If both groups have at least one cell, show success message
  if len(groupA_cells) > 0 and len(groupB_cells) > 0:
      w_text_output(
          content=f"Group A ({groupA_val}): {len(groupA_cells)} cells; "
                  f"Group B ({groupB_val}): {len(groupB_cells)} cells",
          appearance={"message_box": "success"}
      )
      
else:
  w_text_output(content="   ")
  submit_widget_state()
  exit()
