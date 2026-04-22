new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

wf_results_signal()
if wf_results_signal.sample() == True:
  
  if "motif" in results_dict.keys():
  
    w_text_output(content="""## Differential Motif Enrichment""")
    
    m_group = w_select(
        label="Group",
        default=results_dict["motif"]["group_name"].unique()[0],
        options=tuple(results_dict["motif"]["group_name"].unique()),
    )
    
    m_pvals_adj_threshold = w_text_input(
      label="pval adjust threshold",
      default="0.05",
    )
    
    m_meandiff_threshold = w_text_input(
      label="MeanDiff threshold",
      default="0.01",
    )
    
    mcompare_rankby = w_select(
        label="Rank By",
        default="MeanDiff",
        options=tuple(['FDR', 'MeanDiff']),
    )
    
    mcompare_colorby = w_select(
        label="Color By",
        default="FDR",
        options=tuple(['FDR', 'MeanDiff']),
    )
    
    
    mvol_opts = w_row(items=([m_pvals_adj_threshold, m_meandiff_threshold]))
    mrank_opts = w_row(items=[mcompare_rankby, mcompare_colorby])
    
    # ----------------------------------------------------------------------
    
    if m_group.value is not None:
      
      
      group_m = m_group.value
      df_m = results_dict["motif"]
      df_m = df_m[df_m["group_name"] == group_m]
      
      mvol = plot_volcano(
        df_m,
        float(m_pvals_adj_threshold.value),
        float(m_meandiff_threshold.value),
        "GroupA",
        "GroupB",
        pval_key="FDR",
        l2fc_key="MeanDiff",
        names_key="name",
        plot_width=750,
        plot_height=640,
        top_n=2
      )
      mvol_plot = w_plot(source=mvol)
      mvol_col = w_column(items=[mvol_plot, mvol_opts])
      
      mrank = plot_ranked_feature_plotly(
          df_m,
          y_col=mcompare_rankby.value,
          x_col=None,
          n_labels=4,
          label_col="name",
          color_col=mcompare_colorby.value,
          colorscale="PuBu_r",
          marker_size=6,
          title="",
          y_label=mcompare_rankby.value
      )
      mrank_plot = w_plot(source=mrank)
      mrank_col = w_column(items=[mrank_plot, mrank_opts])
      
      with w_grid(columns=2) as m_grid:
          m_grid.add(item=mvol_col, col_span=1)
          m_grid.add(item=mrank_col, col_span=1)

      m_table = w_table(source=df_m)
  
  else:
    w_text_output(
      content="No differential motif analysis found; please check Execution logs.",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()

else:
  w_text_output(
    content="   ",
  )
  submit_widget_state()
