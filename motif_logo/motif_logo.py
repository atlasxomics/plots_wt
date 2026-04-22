new_data_signal()

w_text_output(content="""
# Motif Logo

Display CisBP motif logos for motifs present in the motif AnnData object.
""")

if not adata_g or adata_m is None:
    w_text_output(
        content="Motif data is not loaded...",
        appearance={"message_box": "warning"}
    )
    exit()

default_logo_genome = get_current_igv_genome()
motif_logo_genome = w_select(
    label="genome",
    options=("hg38", "mm10", "rn6"),
    default=default_logo_genome,
    appearance={
        "help_text": "Defaults to the Track Browser genome when available."
    }
)

logo_genome = motif_logo_genome.value
seqlogo_path = resolve_seqlogo_json_path(logo_genome)

if seqlogo_path is None:
    expected_filename = SEQLOGO_JSON_FILENAMES.get(logo_genome, f"seqlogo_{logo_genome}.json")
    w_text_output(
        content=(
            f"No motif logo JSON found for genome `{logo_genome}`. "
            f"Expected `{expected_filename}` in the notebook working directory."
        ),
        appearance={"message_box": "warning"}
    )
    exit()

try:
    motif_logos = get_available_motif_logos(available_motifs, logo_genome)
except Exception as e:
    w_text_output(
        content=f"Failed to load motif logos: {e}",
        appearance={"message_box": "warning"}
    )
    exit()

if len(motif_logos) == 0:
    w_text_output(
        content=(
            f"No overlapping motif logos were found for genome `{logo_genome}` "
            f"in `{seqlogo_path.name}`."
        ),
        appearance={"message_box": "warning"}
    )
    exit()

motif_logo_names = sorted(motif_logos.keys())

motif_logo_select = w_select(
    label="motif",
    options=tuple(motif_logo_names),
    default=motif_logo_names[0],
    appearance={
        "help_text": "Select a motif to display its sequence logo."
    }
)

motif_logo_mode = w_select(
    label="height mode",
    options=("information", "probability"),
    default="information",
    appearance={
        "help_text": "Use information content (bits) or raw probabilities for logo height."
    }
)

motif_logo_show_matrix = w_checkbox(
    label="Display motif matrix",
    default=False,
    appearance={
        "description": "Show the probability matrix used to render the logo."
    }
)

w_row(items=[motif_logo_select, motif_logo_mode, motif_logo_show_matrix])

selected_motif = motif_logo_select.value
selected_logo_df = motif_logos[selected_motif].copy()
logo_title = f"{selected_motif} ({logo_genome})"

motif_logo_fig = plot_motif_logo(
    selected_logo_df,
    title=logo_title,
    information_content=(motif_logo_mode.value == "information"),
)

w_plot(source=motif_logo_fig)

if motif_logo_show_matrix.value:
    matrix_df = selected_logo_df.copy()
    matrix_df.index = np.arange(1, len(matrix_df) + 1)
    matrix_df.index.name = "position"
    w_table(
        label=f"Probability matrix for {selected_motif}",
        source=matrix_df.reset_index()
    )
