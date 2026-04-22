new_data_signal()
if "gene_score_done_signal" not in globals():
    gene_score_done_signal = Signal(False)
gene_score_done_signal()


# Ensure gene activity AnnData is loaded
if not adata_g or not isinstance(adata_g, AnnData):
    w_text_output(content=" ")
    exit()

if gene_score_done_signal.sample() == True:
  w_text_output(content="""### Manually add cell tpe annotations <br>
  Use the H5 Viewer to guide cell type assignment. All cell type scores are stored as Continuous Observations with the suffix '_score.'"""
  )

  w_h5(ann_data=adata_g)
  
  adata_g_obs_df = adata_g.obs
  w_table(label="Cell Metadata", source=adata_g_obs_df)
