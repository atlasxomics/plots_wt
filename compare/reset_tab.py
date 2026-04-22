new_data_signal()

# Abort if no data loaded
if not adata_g:
  w_text_output(content="   ")
  exit()

reset_tab2 = w_button(label="Reset Tab")

if reset_tab2.value:
    # Reset core signals
    choose_group_signal(False)
    groupselect_signal(False)
    barcodes_signal(False)
    wf_exe_signal(False)
    wf_results_signal(False)
    wf_bigwigs_signal(False)

    # Reset metadata widgets to their defaults
    if "choose_h5_data" in globals():
      choose_h5_data._signal(None)
    if "sample_layout_button" in globals():
      sample_layout_button._signal(False)

    if "adata_h5" in globals():
      del adata_h5

    if "h5_cols" in globals():
      h5_cols._signal(None)
    if "h5_rows" in globals():
      h5_rows._signal(None)
    if "h5_spacing" in globals():
      h5_spacing._signal("100.0")    

    if "choose_obs" in globals():
      choose_obs._signal(None)

    if "groupA_ann" in globals():
      groupA_ann._signal(None)
    if "groupB_ann" in globals():
      groupB_ann._signal(None)

    if "wf_name" in globals():
      wf_name._signal("")

    if "groupA_cells" in globals():
      groupA_cells = []
    if "groupB_cells" in globals():
      groupB_cells = []

    # Reset the update plot button
    if "h5_button" in globals():
      h5_button._signal(False)

    if "load_compare_box" in globals():
      load_compare_box._signal(False)
    
    if "load_compare_button" in globals():
      load_compare_button._signal(False)

    if "compare_genome" in globals():
      compare_genome._signal(None)
      
    if "compare_path" in globals():
      compare_path._signal(None)
    
    reset_tab2._signal(False)

    # Ensure all cells initalize
    new_data_signal(True)