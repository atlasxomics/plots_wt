new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

wf_results_signal()
if wf_results_signal.sample() == True:

  if "gene" in results_dict.keys():
    
    w_text_output(content="""## Differential Gene Accessibility""")
    
    g_group = w_select(
        label="Group",
        default=results_dict["gene"]["group_name"].unique()[0],
        options=tuple(results_dict["gene"]["group_name"].unique()),
    )
    
    g_pvals_adj_threshold = w_text_input(
      label="pval adjust threshold",
      default="0.05",
    )
    
    g_log2fc_threshold = w_text_input(
      label="log2fc threshold",
      default="0.01",
    )
    
    gcompare_rankby = w_select(
        label="Rank By",
        default="Log2FC",
        options=tuple(['Log2FC', 'FDR', 'MeanDiff']),
    )
    
    gcompare_colorby = w_select(
        label="Color By",
        default="FDR",
        options=tuple(['Log2FC', 'FDR', 'MeanDiff']),
    )
    
    
    gvol_opts = w_row(items=([g_pvals_adj_threshold, g_log2fc_threshold]))
    grank_opts = w_row(items=[gcompare_rankby, gcompare_colorby])
    
    # ----------------------------------------------------------------------
    
    if g_group.value is not None:
      
      
      group_g = g_group.value
      df_g = results_dict["gene"]
      df_g = df_g[df_g["group_name"] == group_g]
      
      gvol = plot_volcano(
        df_g,
        float(g_pvals_adj_threshold.value),
        float(g_log2fc_threshold.value),
        "GroupA",
        "GroupB",
        pval_key="FDR",
        l2fc_key="Log2FC",
        names_key="name",
        plot_width=750,
        plot_height=640,
        top_n=2
      )
      gvol_plot = w_plot(source=gvol)
      gvol_col = w_column(items=[gvol_plot, gvol_opts])
      
      grank = plot_ranked_feature_plotly(
          df_g,
          y_col=gcompare_rankby.value,
          x_col=None,
          n_labels=4,
          label_col="name",
          color_col=gcompare_colorby.value,
          colorscale="PuBu_r",
          marker_size=6,
          title="",
          y_label=gcompare_rankby.value
      )
      grank_plot = w_plot(source=grank)
      grank_col = w_column(items=[grank_plot, grank_opts])
      
      with w_grid(columns=2) as grid:
          grid.add(item=gvol_col, col_span=1)
          grid.add(item=grank_col, col_span=1)

      g_table = w_table(source=df_g)

  else:
    w_text_output(
      content="No differential gene analysis found; please check Execution logs.",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()
else:
  w_text_output(
    content="   ",
  )
  submit_widget_state()
