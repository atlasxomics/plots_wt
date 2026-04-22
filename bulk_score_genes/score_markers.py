new_data_signal()

if adata_g is None:
  w_text_output(content="   ")
  exit()

if "choose_subset_signal" not in globals():
  choose_subset_signal = Signal(False)

choose_subset_signal()
if choose_subset_signal.sample():
  
  # --- Step 4: Button to trigger scoring ---
  construct_button = w_button(label="Compute Gene Set Score", key="gene_set_button")
  skip_unassigned_checkbox = w_checkbox(
      label="Do not label low-confidence cells as 'Unassigned'",
      default=False,
      appearance={"help_text": "If checked, low-confidence cells keep their predicted type instead of 'Unassigned'"}
  )
  
  if construct_button.value:
      # Build marker dictionary
      marker_dict = {}
      for lbl_w, feat_w in zip(label_inputs, feature_selects):
          lbl = lbl_w.value
          feats = feat_w.value or []
          if lbl and isinstance(feats, list) and feats:
              marker_dict[lbl] = feats
  
      w_text_output(
          content=f"Constructed cell type-feature set dictionary: {marker_dict}",
          appearance={"message_box": "info"}
      )
      submit_widget_state()
  
      # --- Compute gene set scores on the subset ---
      for cell_type, genes in marker_dict.items():
          valid_genes = [g for g in genes if g in adata_subset.var_names]
          w_text_output(
              content=f"Starting feature scoring for {cell_type}...",
              appearance={"message_box": "info"}
          )
          if len(valid_genes) > 0:
            submit_widget_state()
            if valid_genes:
                sc.tl.score_genes(
                    adata_subset,
                    gene_list=valid_genes,
                    score_name=f"{cell_type}_score"
                )
                w_text_output(
                    content=f"Feature scoring finished for {cell_type}.",
                    appearance={"message_box": "success"}
                )
                submit_widget_state()
            else:
                w_text_output(
                    content=f"Cell type {cell_type} contains no valid features.",
                    appearance={"message_box": "warning"}
                )

      # --- Assign best-matching cell type per cell in the subset ---
      score_cols = [
          f"{ct}_score"
          for ct in marker_dict.keys()
          if f"{ct}_score" in adata_subset.obs
      ]
      if score_cols:
          w_text_output(
              content="Assigning predicted cell types...",
              appearance={"message_box": "info"}
          )
          submit_widget_state()
          scores = adata_subset.obs[score_cols].to_numpy()
          best_idx = np.argmax(scores, axis=1)
          best_types = np.array(score_cols)[best_idx]
          adata_subset.obs["pred_cell_type"] = [
              bt.replace("_score", "") for bt in best_types
          ]
          # Compute confidence margin
          if len(marker_dict) > 1:
            second_best = np.partition(scores, -2, axis=1)[:, -2]
            margin = scores[np.arange(scores.shape[0]), best_idx] - second_best
            adata_subset.obs["pred_cell_type_conf"] = margin

            # Optionally label low-confidence cells as Unassigned
            if not skip_unassigned_checkbox.value:
                low_conf_mask = adata_subset.obs["pred_cell_type_conf"] < 0.01
                adata_subset.obs.loc[low_conf_mask, "pred_cell_type"] = "Unassigned"
          else:
              adata_subset.obs["pred_cell_type_conf"] = 1.0
              w_text_output(
                  content="Only one cell type used, all cells are assigned to that cell type.",
                  appearance={"message_box": "warning"}
              )
              submit_widget_state()
  
      # --- Merge scores and predictions back into adata_g.obs, only for subset cells ---
      # Ensure columns exist
      for col in score_cols:
          if col not in adata_g.obs.columns:
              adata_g.obs[col] = np.nan
      if "pred_cell_type" not in adata_g.obs.columns:
          adata_g.obs["pred_cell_type"] = None
      if "pred_cell_type_conf" not in adata_g.obs.columns:
          adata_g.obs["pred_cell_type_conf"] = np.nan
  
      # Populate only subset cells
      for col in score_cols:
          adata_g.obs.loc[adata_subset.obs_names, col] = adata_subset.obs[col]
      adata_g.obs["pred_cell_type"] = adata_g.obs["pred_cell_type"].astype(str)
      adata_g.obs.loc[adata_subset.obs_names, "pred_cell_type"] = adata_subset.obs["pred_cell_type"]
      adata_g.obs.loc[adata_subset.obs_names, "pred_cell_type_conf"] = adata_subset.obs["pred_cell_type_conf"]
  
      w_text_output(
          content="Finished assigning cell types based on best score",
          appearance={"message_box": "success"}
      )
      submit_widget_state()
      gene_score_done_signal(True)

else:
  w_text_output(content="   ")
  submit_widget_state()
  exit()
  
