new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

wf_bigwigs_signal()
  
if wf_bigwigs_signal.sample() == True:

  igv_genome = None

  if "load_compare_box" in globals() and load_compare_box.value:
      igv_genome = compare_genome.value if "compare_genome" in globals() else None
  elif "coverages_genome" in globals():
      igv_genome = coverages_genome.value

  if igv_genome is None:
      w_text_output(
          content="Please ensure genome is selected; set genome on Track Browser tab.",
          appearance={"message_box": "warning"},
      )
      submit_widget_state()
      exit()
  
  if igv_genome is None:
    w_text_output(
      content="Please ensure genome is selected.",
      appearance={"message_box": "warning"}
    )
    submit_widget_state()
    exit()

  tracks = []
  for f in files:
      tracks.append({
          "name": f.path.split("/")[-1],
          "type": "wig",
          "format": "bigwig",
          "url": f.path,
          "autoscale": False,
          "visibilityWindow": 1000000000,
      })
  
  opts: IGVOptions = {
      "genome": igv_genome,
      "tracks": tracks,
  }
  
  igviewer = w_igv(options=opts)

else:
  w_text_output(
    content="   ",
  )
  submit_widget_state()