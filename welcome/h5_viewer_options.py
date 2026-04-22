# w_text_output(content="## H5 Viewer Settings")

new_data_signal()

# Abort if no data loaded
if adata_g is None:
    w_text_output(
        content=" ",
    )
    exit()

if "adata_h5" not in globals() or adata_h5 is None:
    adata_h5 = adata_g
    loaded_h5_data_key = "adata"

w_text_output(content="""
H5 Viewer Advanced Options
<ul>
<li><b>Refresh H5 Viewer</b>: reload the selected AnnData object and apply any spatial layout changes.</li>
<li><b>Save H5 Data</b>: write the currently displayed AnnData object to disk and upload it back to Latch Data.</li>
</ul>
""")

sample_layout_available = (
    "sample" in adata_g.obs
    and "spatial" in adata_g.obsm_keys()
)

sample_layout_button = None
if sample_layout_available:
    sample_layout_button = w_checkbox(
        label="Change H5 Viewer spatial arrangement",
        default=False,
        appearance={"description": "Toggle to specify rows, columns, and spacing."}
    )
else:
    w_text_output(
        content="Spatial arrangement controls are unavailable because this AnnData object does not contain both `obs['sample']` and `obsm['spatial']`.",
        appearance={"message_box": "info"}
    )

if sample_layout_button is not None and sample_layout_button.value:

    h5_sortby_opts = ["original", "sample"]
    if "condition" in adata_g.obs:
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
save_button = w_button(label="Save H5 Data")

with w_grid(columns=2) as h5_grid:
    h5_grid.add(item=h5_button, col_span=1)
    h5_grid.add(item=save_button, col_span=1)

if save_button.value:
    save_path = adata_path if "adata_path" in globals() and adata_path is not None else adata_g_path
    if save_path is None:
      w_text_output(
        content="Could not determine the source H5AD path for saving.",
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
    except Exception as e:
      w_text_output(
        content=f"Write to disk failed: {e}",
        key="writing_failed",
        appearance={"message_box": "warning"}
      )
      submit_widget_state()
      exit()

    w_text_output(
      content="Uploading data to Latch Data...",
      key="upload_message",
      appearance={"message_box": "info"}
    )
    try:
      save_path.upload_from(Path(save_path.name()))
    except Exception as e:
      w_text_output(
        content=f"Upload failed: {e}",
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

    adata_h5 = adata_g
    loaded_h5_data_key = "adata"

    if sample_layout_button is not None and sample_layout_button.value:

        proceed = True

        cols_val = (h5_cols.value or "").strip() if h5_cols is not None else ""
        rows_val = (h5_rows.value or "").strip() if h5_rows is not None else ""
        spacing_val = (h5_spacing.value or "").strip() if h5_spacing is not None else ""

        if h5_flipy is not None and isinstance(h5_flipy.value, bool):
            flipy_val = h5_flipy.value
        else:
            flipy_val = False

        valid_sort_modes = set(h5_sortby_opts)
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

            if proceed:
                if n_cols < 1 or n_rows < 1 or spacing <= 0:
                    proceed = False
                    w_text_output(
                        content="Rows/Columns must be >= 1 and Spacing must be > 0; ignoring...",
                        appearance={"message_box": "warning"}
                    )

            if proceed:
              n_samples = adata_h5.obs["sample"].nunique()
              total_positions = n_rows * n_cols
              if n_samples > total_positions:
                proceed = False
                w_text_output(
                  content=f"Not enough grid positions ({n_rows}x{n_cols}={total_positions}) for {n_samples} samples",
                  appearance={"message_box": "warning"}
                )

            if proceed:
                new_obsm = f"X_dataset_{n_rows}x{n_cols}-{spacing}-{('FlipY' if flipy_val else 'noFlipY')}-{sort_val}"
                try:
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
                except Exception as e:
                  w_text_output(
                    content=f"Could not create spatial layout: {e}",
                    appearance={"message_box": "warning"}
                  )
        else:
            w_text_output(
                content="Please complete all fields in 'Change H5 Viewer spatial arrangement' to specify layout.",
                appearance={"message_box": "warning"}
            )
    refresh_h5_signal(True)
