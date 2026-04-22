new_data_signal()

if not adata_g:
    w_text_output(
        content=" ",
    )
    exit()

if coverages_group.value is not None and coverages_genome.value is not None:

    coverages_dir = coverages_dict[coverages_group.value]

    bedgraphs = []
    bigwigs = []
    for file in coverages_dir.iterdir():
        suffix = file.path.split(".")[-1]
        if suffix == "gz":
            bedgraphs.append(file)
        elif suffix == "bw":
            bigwigs.append(file)

    bedgraphs = sorted(bedgraphs, key=lambda x: x.path.split('/')[-1].split('_')[0])
    bigwigs = sorted(bigwigs, key=lambda x: x.path.split('/')[-1].split('_')[0])

    files = []
    if len(bedgraphs) + len(bigwigs) == 0:
        w_text_output(
            content="No bedgraph (.gz) or bigwig (.bw) files found in coverages directory...",
            appearance={"message_box": "warning"}
        )
        submit_widget_state()
        exit()
    elif len(bedgraphs) > 0 and len(bigwigs) == 0:
        w_text_output(
            content="Using bigwigs instead of bedgraph could make IGV much faster...",
            appearance={"message_box": "warning"}
        )
        submit_widget_state()
        format = "bedgraph"
        files = bedgraphs
    elif len(bedgraphs) == 0 and len(bigwigs) > 0:
        format = "bigwig"
        files = bigwigs
    elif len(bedgraphs) > 0 and len(bigwigs) > 0:
        submit_widget_state()
        format = "bigwig"
        files = bigwigs

    tracks = []
    for f in files:
        tracks.append({
            "name": f.path.split("/")[-1],
            "type": "wig",
            "format": format,
            "url": f.path,
            "autoscale": False,
            "visibilityWindow": 1000000000,
        })

    opts: IGVOptions = {
        "genome": coverages_genome.value,
        "tracks": tracks,
    }

    igviewer = w_igv(options=opts)
