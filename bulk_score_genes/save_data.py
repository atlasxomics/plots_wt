new_data_signal()

if not adata_g:
  w_text_output(content="  ")
  submit_widget_state()
  exit()

if "gene_score_done_signal" not in globals():
  gene_score_done_signal = Signal(False)

gene_score_done_signal()
if gene_score_done_signal.sample() is True:

  w_text_output(content="""
  Click **Save H5AD Data** to save your custom H5 Viewer annotations; new annotations with be available in your next session. <br>  Click **Copy Annotations to motif data** to add new annotations to motif AnnData object.
  """)

  save_button_gs = w_button(label="Save H5AD Data")
  copy_to_motif_button = w_button(label="Copy Annotations to motif data")
  save_warning_gs = w_text_output(content="""_This operation may take a couple minutes._""")

  save_col_gs = w_column(items=[save_button_gs, copy_to_motif_button, save_warning_gs])

  if copy_to_motif_button.value:
      if adata_m is None:
        w_text_output(
          content="Motif AnnData is not loaded, so annotations could not be copied.",
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
          content=f"Annotation copy failed: {e}",
          appearance={"message_box": "warning"}
        )
        submit_widget_state()
        exit()

      w_text_output(
        content="Copied non-numeric annotations from gene to motif AnnData.",
        appearance={"message_box": "success"}
      )
      submit_widget_state()

  if save_button_gs.value:
      save_path_gs = adata_g_path

      w_text_output(
        content="Writing data to disk...",
        key="writing_message",
        appearance={"message_box": "info"}
      )
      submit_widget_state()
      try:
        adata_g.write(save_path_gs.name())
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
        save_path_gs.upload_from(Path(save_path_gs.name()))
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
