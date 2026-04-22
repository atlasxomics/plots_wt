new_data_signal()

if adata_g is None:
  w_text_output(content="  ")
  submit_widget_state()
  exit()

if "gene_score_done_signal" not in globals():
  gene_score_done_signal = Signal(False)

gene_score_done_signal()
if gene_score_done_signal.sample() is True:

  w_text_output(content="""
  Click **Save H5AD Data** to save your custom H5 Viewer annotations; new annotations will be available in your next session.
  """)

  save_button_gs = w_button(label="Save H5AD Data")
  save_warning_gs = w_text_output(content="""_This operation may take a couple minutes._""")

  save_col_gs = w_column(items=[save_button_gs, save_warning_gs])

  if save_button_gs.value:
      save_path_gs = adata_path if "adata_path" in globals() and adata_path is not None else None
      if save_path_gs is None:
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
        adata_g.write(save_path_gs.name())
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
        save_path_gs.upload_from(Path(save_path_gs.name()))
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
