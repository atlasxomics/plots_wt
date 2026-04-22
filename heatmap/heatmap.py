new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

heatmap_signal()
if heatmap_signal.sample() is True and choose_heatmap_data.value is not None:
  notebook_palettes = await get_notebook_palettes()
  hm_palette = w_select(
    key="hm_palette",
    label="colorscale",
    default="Default Heatmap Colorscale",
    options=get_palette_selector_options(
      notebook_palettes,
      kind="continuous",
      fallback_name="Default Heatmap Colorscale",
    ),
    appearance={
      "help_text": "Use a continuous palette saved from the H5 Viewer or fall back to the default heatmap colors."
    }
  )

  adata_hm = h5data_dict[choose_heatmap_data.value]
  hm_feats = choose_heatmap_data.value

  hm_group_widget = w_select(
    key="hm_group_widget",
    label="group",
    default="cluster",
    options=tuple(groups),
    appearance={
      "detail": "(cluster, sample, condition)",
      "help_text": "Select grouping for heatmap."
    }
  )
  hm_group = hm_group_widget.value
  hm_top_row = w_row(items=[hm_group_widget, hm_palette])

  try:
    feature_label, stats_key, stats_df, _ = resolve_heatmap_stats_table(
      adata_hm, hm_feats, hm_group
    )
  except ValueError as e:
    w_text_output(
      content=str(e),
      appearance={"message_box": "warning"}
    )
    exit()

  if stats_key.startswith("marker_genes_per_"):
    w_text_output(
      content=(
        "Using a prefiltered marker table for this heatmap. "
        "Some genes may be omitted from the available statistics."
      ),
      appearance={"message_box": "warning"}
    )

  # Motif stats can be emitted in ArchR's wide format:
  # each row is a metric (mlog10Padj, Enrichment, etc.), columns are groups.
  if hm_feats == "motif" and is_motif_wide_stats_table(stats_df):
    hm_sig_threshold = w_text_input(
      key="hm_sig_threshold",
      label="Adjusted p-value threshold",
      default="0.01",
      appearance={"help_text": "Converted to -log10(padj) internally for motif filtering."}
    )
    hm_top_n = w_text_input(
      key="hm_top_n",
      label="Top motifs per group",
      default="25",
      appearance={"help_text": "Used only when feature list is empty."}
    )
    hm_feature_list = w_text_input(
      key="hm_feature_list",
      label="Motif list (optional, comma-separated)",
      default="",
      appearance={"help_text": "If provided, this list overrides Top motifs per group."}
    )
    controls_row1 = w_row(items=[hm_sig_threshold, hm_top_n])
    controls_row2 = w_row(items=[hm_feature_list])

    sig_threshold = safe_float(
      hm_sig_threshold.value,
      warn_msg="Adjusted p-value threshold must be numeric; defaulting to 0.01."
    )
    if sig_threshold is None or sig_threshold <= 0:
      sig_threshold = 0.01

    try:
      top_n = int(hm_top_n.value)
      if top_n < 1:
        raise ValueError
    except (TypeError, ValueError):
      w_text_output(
        content="Top motifs per group must be a positive integer; defaulting to 25.",
        appearance={"message_box": "warning"}
      )
      top_n = 25

    feature_input = hm_feature_list.value if hm_feature_list.value is not None else ""
    group_labels = None
    if hm_group in adata_hm.obs.columns:
      group_labels = sort_group_categories(
        adata_hm.obs[hm_group].astype(str).unique().tolist()
      )

    try:
      heatmap_df, legend_title, motif_meta = build_motif_archr_like_heatmap_from_wide(
        stats_df,
        sig_threshold=sig_threshold,
        top_n=top_n,
        feature_input=feature_input,
        group_labels=group_labels,
      )
    except ValueError as e:
      w_text_output(
        content=str(e),
        appearance={"message_box": "warning"}
      )
      exit()

    selected = motif_meta["selected_features"]
    if len(selected) > 0:
      ordered_rows = [f for f in selected if f in heatmap_df.index]
      ordered_rows += [f for f in heatmap_df.index if f not in ordered_rows]
      heatmap_df = heatmap_df.reindex(ordered_rows)

    title = f"Motif heatmap by {hm_group}"
    title += f" | metric={motif_meta['metric']}"
    title += f" | cutoff>={motif_meta['mlog10_cutoff']:.2f}"
    title += f" | top_n={top_n}"
  else:
    try:
      (
        group_col,
        feature_col,
        sig_col,
        selected_sig_metric,
        available_sig_metrics,
      ) = get_heatmap_stats_columns(
        stats_df, stats_key
      )
    except ValueError as e:
      w_text_output(
        content=str(e),
        appearance={"message_box": "warning"}
      )
      exit()

    if len(available_sig_metrics) > 0:
      hm_sig_metric = w_select(
        key="hm_sig_metric",
        label="Significance statistic",
        default=selected_sig_metric,
        options=available_sig_metrics,
        appearance={
          "help_text": "Choose whether the significance cutoff uses FDR or Pval."
        }
      )
      selected_sig_metric = hm_sig_metric.value
      _, _, sig_col, selected_sig_metric, _ = get_heatmap_stats_columns(
        stats_df,
        stats_key,
        preferred_sig_metric=selected_sig_metric,
      )
    else:
      selected_sig_metric = None
      sig_col = None

    sig_threshold_label = "Significance threshold"
    sig_threshold_help = "Marker significance cutoff using FDR."
    if selected_sig_metric == "Pval":
      sig_threshold_label = "Significance threshold)"
      sig_threshold_help = "Marker significance cutoff using unadjusted p-values."
    elif selected_sig_metric is None:
      sig_threshold_label = "Significance threshold"
      sig_threshold_help = "No significance column found; this cutoff will not be applied."

    hm_sig_threshold = w_text_input(
      key="hm_sig_threshold",
      label=sig_threshold_label,
      default="0.01",
      appearance={"help_text": sig_threshold_help}
    )

    hm_effect_threshold = w_text_input(
      key="hm_effect_threshold",
      label="Log2FoldChange threshold",
      default="0.5" if hm_feats == "gene" else "0.1",
      appearance={
        "help_text": "Directional Log2FoldChange cutoff for selecting markers."
      }
    )

    hm_effect_direction = w_select(
      key="hm_effect_direction",
      label="Log2FoldChange direction",
      default="positive",
      options=("positive", "negative", "absolute"),
      appearance={
        "help_text": "Direction used when selecting top markers per group; use 'positive' for upregulated, 'negative' for downregulated."
      }
    )

    hm_z_clip = w_text_input(
      key="hm_z_clip",
      label="Gene z-score clip (max/min)",
      default="2.0",
      appearance={
        "help_text": "Max/min of color scale, only applies to gene heatmaps."
      }
    )

    hm_top_n = w_text_input(
      key="hm_top_n",
      label="Top features per group",
      default="25",
      appearance={"help_text": "Used only when feature list is empty."}
    )

    hm_feature_list = w_text_input(
      key="hm_feature_list",
      label="Feature list (optional, comma-separated)",
      default="",
      appearance={"help_text": "If provided, this list overrides Top features per group."}
    )

    row2_items = [hm_sig_threshold]
    if len(available_sig_metrics) > 0:
      row2_items = [hm_sig_metric] + row2_items
    controls_row2 = w_row(items=row2_items)
    controls_row3 = w_row(items=[hm_effect_threshold, hm_effect_direction])
    controls_row4 = w_row(items=[hm_z_clip, hm_top_n, hm_feature_list])

    value_metric = "Log2FC"
    rank_metric = "Log2FC"
    effect_direction = hm_effect_direction.value

    if sig_col is not None:
      sig_threshold = safe_float(
        hm_sig_threshold.value,
        warn_msg="Significance threshold must be numeric; defaulting to 0.01."
      )
      if sig_threshold is None:
        sig_threshold = 0.01
    else:
      sig_threshold = None

    effect_threshold = safe_float(
      hm_effect_threshold.value,
      warn_msg="Effect-size threshold must be numeric; using default."
    )
    if effect_threshold is None:
      effect_threshold = 0.5 if hm_feats == "gene" else 0.1

    z_clip = safe_float(
      hm_z_clip.value,
      warn_msg="Row z-score clip must be numeric; defaulting to 2.0."
    )
    if z_clip is None or z_clip <= 0:
      z_clip = 2.0

    try:
      top_n = int(hm_top_n.value)
      if top_n < 1:
        raise ValueError
    except (TypeError, ValueError):
      w_text_output(
        content="Top features per group must be a positive integer; defaulting to 25.",
        appearance={"message_box": "warning"}
      )
      top_n = 25

    feature_input = hm_feature_list.value if hm_feature_list.value is not None else ""
    try:
      work_df = prepare_heatmap_work_df(
        stats_df=stats_df,
        group_col=group_col,
        feature_col=feature_col,
        value_metric=value_metric,
        rank_metric=rank_metric,
        sig_col=sig_col,
        sig_threshold=sig_threshold,
      )
      work_df, selected = select_archr_like_heatmap_features(
        work_df=work_df,
        group_col=group_col,
        feature_col=feature_col,
        rank_metric=rank_metric,
        top_n=top_n,
        effect_threshold=effect_threshold,
        effect_direction=effect_direction,
        feature_input=feature_input,
      )
      heatmap_df, legend_title = build_archr_like_heatmap_df(
        work_df=work_df,
        hm_feats=hm_feats,
        group_col=group_col,
        feature_col=feature_col,
        value_metric=value_metric,
        sig_col=sig_col,
        z_clip=z_clip,
      )
    except ValueError as e:
      w_text_output(
        content=str(e),
        appearance={"message_box": "warning"}
      )
      exit()

    if len(selected) > 0:
      ordered_rows = [f for f in selected if f in heatmap_df.index]
      ordered_rows += [f for f in heatmap_df.index if f not in ordered_rows]
      heatmap_df = heatmap_df.reindex(ordered_rows)

    title = f"{feature_label.capitalize()} heatmap by {hm_group}"
    if sig_col is not None and sig_threshold is not None:
      title += f" | {selected_sig_metric}<={sig_threshold}"
    if effect_direction == "positive":
      title += f" | {rank_metric}>={effect_threshold}"
    elif effect_direction == "negative":
      title += f" | {rank_metric}<=-{effect_threshold}"
    else:
      title += f" | |{rank_metric}|>={effect_threshold}"

  heatmap = px.imshow(
    heatmap_df,
    color_continuous_scale=get_selected_continuous_palette(
      notebook_palettes,
      hm_palette.value,
      fallback_colors="Spectral_r",
      fallback_name="Default Heatmap Colorscale",
    ),
    aspect="auto",
    origin="lower"
  )

  heatmap.update_layout(
    title=title,
    xaxis_title=hm_group,
    yaxis_title=f"{feature_label}s",
    coloraxis_colorbar=dict(
      title=legend_title,
      title_side="right"
    )
  )

  heatmap.update_xaxes(
    side="bottom",
    tickmode="array",
    tickvals=list(range(len(heatmap_df.columns))),
    ticktext=heatmap_df.columns,
    tickangle=45
  )
  heatmap.update_yaxes(autorange="reversed")

  w_plot(source=heatmap)
