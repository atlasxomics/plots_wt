# w_text_output(content="## H5 Viewer Settings")

new_data_signal()

# Abort if no data loaded
if not adata_g:
    w_text_output(
        content=" ",
    )
    exit()

w_text_output(content="""
H5 Viewer Advanced Options
<ul>
<li><b>Refresh H5 Viewer</b>: reload the selected gene or motif AnnData object and apply any spatial layout changes.</li>
<li><b>Synch Metadata</b>: copy missing non-numeric <code>obs</code> columns between the gene and motif AnnData objects.</li>
<li><b>Save H5 Data</b>: write the currently displayed H5 Viewer AnnData object to disk and upload it back to Latch Data.</li>
</ul>
""")

# Choose whether to display gene or motif data
choose_h5_data = w_select(
    label="Select Data for H5 Viewer",
    default="gene",
    options=["gene", "motif"],
    appearance={
        "help_text": "Select which features to display in the H5 Viewer."
    }
)

# Checkbox to show/hide layout controls
sample_layout_button = w_checkbox(
    label="Change H5 Viewer spatial arrangement",
    default=False,
    appearance={"description": "Toggle to specify rows, columns, and spacing."}
)

if sample_layout_button.value:

    h5_sortby_opts = ["original", "sample"]
    if "condition" in groups:
      h5_sortby_opts.append("condition")

    h5_cols = w_text_input(
        label="Number of Columns",
        key="h5_cols",
        default=None,
        appearance={"help_text": "Specify the number of columns in the layout."}
    )
    h5_rows = w_text_input(
        label="Number of Rows",
        key="h5_rows",
        default=None,
        appearance={"help_text": "Specify the number of rows in the layout."}
    )
    h5_spacing = w_text_input(
        label="Spacing Between Samples",
        key="h5_spacing",
        default="100.0",
        appearance={"help_text": "Specify the spacing between samples."}
    )

    h5_flipy = w_checkbox(
        label="Flip Y Axis",
        key="h5_flipy",
        default=False,
        appearance={"description": "Rotate samples around the y axis."}
    )

    h5_sortby = w_select(
        label="Sort Samples By",
        key="h5_sortby",
        default="original",
        options=tuple(h5_sortby_opts),
        appearance={
          "help_text": "Sort samples alphabetically or by condition.",
          "description": "'original' maintains original order."
        }
    )

    layout_col1 = w_column(items=[h5_cols, h5_rows])
    layout_col2 = w_column(items=[h5_spacing, h5_sortby, h5_flipy])

    with w_grid(columns=5) as layout_grid:
        layout_grid.add(item=layout_col1, col_span=1)
        layout_grid.add(item=layout_col2, col_span=4)

# Button to start the H5 viewer
h5_button = w_button(label="Refresh H5 Viewer")
h5_data_col = w_column(items=[h5_button])

sync_button = w_button(label="Synch Metadata")

save_button = w_button(label="Save H5 Data")

right_col = w_column(items=[sync_button, save_button])

with w_grid(columns=2) as h5_grid:
    h5_grid.add(item=h5_data_col, col_span=1)
    h5_grid.add(item=right_col, col_span=1)

if sync_button.value:
    if adata_m is None:
        w_text_output(
            content="Motif AnnData is not loaded, so metadata sync could not run.",
            appearance={"message_box": "warning"}
        )
        submit_widget_state()
        exit()

    try:
        sync_obs_metadata(
            adata_g,
            adata_m,
            reconcile_shared=False,
        )
    except ValueError as e:
        w_text_output(
            content=f"Metadata sync failed: {e}",
            appearance={"message_box": "warning"}
        )
        submit_widget_state()
        exit()

    w_text_output(
        content="Gene and motif metadata synchronized for non-numeric obs columns.",
        appearance={"message_box": "success"}
    )
    refresh_h5_signal(True)
    submit_widget_state()

if save_button.value:
    loaded_key = loaded_h5_data_key if "loaded_h5_data_key" in globals() else None
    if loaded_key is None:
      if adata_h5 is adata_m:
        loaded_key = "motif"
      else:
        loaded_key = "gene"

    if choose_h5_data.value != loaded_key:
      w_text_output(
        content=(
          f"The H5 Viewer currently has {loaded_key} data loaded. "
          f"Saving the loaded object instead of the selected {choose_h5_data.value} object."
        ),
        appearance={"message_box": "warning"}
      )
      submit_widget_state()

    if loaded_key == "gene":
      save_path = adata_g_path
    elif loaded_key == "motif":
      save_path = adata_m_path
    else:
      w_text_output(
        content="Could not determine which H5 object is currently loaded.",
        appearance={"message_box": "warning"}
      )
      submit_widget_state()
      exit()

    w_text_output(
      content="Writing data to disk...",
      key="writing_message",
      appearance={"message_box": "info"}
    )
    submit_widget_state()
    try:
      adata_h5.write(save_path.name())
    except:
      w_text_output(
        content="Write to disk failed...",
        key="writing_failed",
        appearance={"message_box": "warning"}
      )
      submit_widget_state()
      exit()

    upload_message = w_text_output(
      content="Uploading data to Latch Data...",
      key="upload_message",
      appearance={"message_box": "info"}
    )
    try:
      save_path.upload_from(Path(save_path.name()))
    except:
      w_text_output(
        content="Upload failed...",
        key="upload_failed",
        appearance={"message_box": "warning"}
      )
      submit_widget_state()
      exit()

    w_text_output(
      content="Upload success!",
      key="upload_success",
      appearance={"message_box": "success"}
    )
    submit_widget_state()

if h5_button.value:

    adata_h5 = h5data_dict[choose_h5_data.value]
    loaded_h5_data_key = choose_h5_data.value
    
    if sample_layout_button.value:
    
        proceed = True
    
        # Ensure fields are present and non-empty after stripping
        cols_val = (h5_cols.value or "").strip() if h5_cols is not None else ""
        rows_val = (h5_rows.value or "").strip() if h5_rows is not None else ""
        spacing_val = (h5_spacing.value or "").strip() if h5_spacing is not None else ""

        if h5_flipy is not None and isinstance(h5_flipy.value, bool):
            flipy_val = h5_flipy.value
        else:
            flipy_val = False
        
        valid_sort_modes = {"original", "sample", "condition"}
        if h5_sortby is not None and (h5_sortby.value in valid_sort_modes):
            sort_val = h5_sortby.value
        else:
            sort_val = "original"
    
        if cols_val and rows_val and spacing_val:
    
            try:
                n_cols = int(cols_val)
            except (TypeError, ValueError):
                proceed = False
                w_text_output(
                    content="Cannot convert 'Number of Columns' input into an integer; ignoring...",
                    appearance={"message_box": "warning"}
                )
    
            try:
                n_rows = int(rows_val)
            except (TypeError, ValueError):
                proceed = False
                w_text_output(
                    content="Cannot convert 'Number of Rows' input into an integer; ignoring...",
                    appearance={"message_box": "warning"}
                )
    
            try:
                spacing = float(spacing_val)
            except (TypeError, ValueError):
                proceed = False
                w_text_output(
                    content="Cannot convert 'Spacing Between Samples' input into a float; ignoring...",
                    appearance={"message_box": "warning"}
                )
    
            # Optional: basic range checks
            if proceed:
                if n_cols < 1 or n_rows < 1 or spacing <= 0:
                    proceed = False
                    w_text_output(
                        content="Rows/Columns must be ≥ 1 and Spacing must be > 0; ignoring...",
                        appearance={"message_box": "warning"}
                    )

            if proceed:
              total_positions = n_rows * n_cols
              if len(samples) > total_positions:
                proceed = False
                w_text_output(
                  content=f"Not enough grid positions ({n_rows}x{n_cols}={total_positions}) for {len(samples)} samples",
                  appearance={"message_box": "warning"}
                )
    
            if proceed:
                new_obsm = f"spatial_offset_{n_rows}x{n_cols}-{spacing}-{('FlipY' if flipy_val else 'noFlipY')}-{sort_val}"
                process_matrix_layout(
                    adata_h5,
                    n_rows=n_rows,
                    n_cols=n_cols,
                    tile_spacing=spacing,
                    flipy=flipy_val,
                    sample_order_mode=sort_val,
                    new_obsm_key=new_obsm
                )
                refresh_h5_signal(True)
        else:
            w_text_output(
                content="Please complete all fields in 'Change H5 Viewer spatial arrangement' to specify layout.",
                appearance={"message_box": "warning"}
            )
    refresh_h5_signal(True)
