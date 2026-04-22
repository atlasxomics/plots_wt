new_data_signal()

if not adata_g:
  w_text_output(content="   ")
  exit()

w_text_output(content="""## Set Inputs for Compare Workflow""")

barcodes_signal()
if barcodes_signal.sample() == True:

  wf_name = w_text_input(
    key="wf_name",
    label="Output Directory Name",
    default="",
    appearance={
      "help_text": "Name of output dirctory results save in Latch Data in the directory `compare_outs`."
    }
  )

  show_workflow_params = w_checkbox(
    key="show_workflow_params",
    label="Show optional workflow parameters",
    default=False,
    appearance={
      "description": "Expand to configure max_cells and use_max_possible_cells."
    }
  )
  max_cells_value = 500
  use_max_possible_cells_value = False

  if show_workflow_params.value:
    max_cells = w_text_input(
      key="max_cells",
      label="Maximum cells per group",
      default="500",
      appearance={
        "help_text": "Upper bound on the number of cells used per group in the compare workflow."
      }
    )

    use_max_possible_cells = w_checkbox(
      key="use_max_possible_cells",
      label="Use maximum possible cells",
      default=False,
      appearance={
        "description": "Ignore the explicit max_cells cap and use the largest possible matched group size."
      }
    )
    workflow_controls = w_row(items=[max_cells, use_max_possible_cells])

    if use_max_possible_cells.value:
      w_text_output(
        content="`Use maximum possible cells` is enabled. The workflow will ignore `max_cells` and use the largest possible matched group size.",
        appearance={"message_box": "info"}
      )

    max_cells_input = (
      max_cells.value.strip()
      if max_cells.value is not None else ""
    )
    if len(max_cells_input) > 0:
      try:
        max_cells_value = int(max_cells_input)
        if max_cells_value < 1:
          raise ValueError
      except (TypeError, ValueError):
        w_text_output(
          content="Maximum cells per group must be a positive integer; defaulting to 500.",
          appearance={"message_box": "warning"}
        )
        max_cells_value = 500

    use_max_possible_cells_value = use_max_possible_cells.value

  if (wf_name.value is not None and
      len(wf_name.value) > 0 and
      archrproj_dir is not None and
      remote_bcs is not None
    ):

    params = {
        "project_name": wf_name.value,
        "groupings": LatchFile(remote_bcs.path),
        "archrproject": LatchDir(archrproj_dir.path),
        "max_cells": max_cells_value,
        "use_max_possible_cells": use_max_possible_cells_value,
    }

    w = w_workflow(
      wf_name="wf.__init__.compare_workflow",
      version="0.10.4-81336e-wip-d3e423",
      params=params,
      label="Launch Workflow"
    )

    wf_exe_signal(True)
    execution = w.value
    
  else:
    w_text_output(
      content="Please set Workflow inputs.",
      appearance={"message_box": "info"}
    )
    submit_widget_state()
    wf_exe_signal(False)    
  
else:
  w_text_output(
    content="Please ensure cells are selected for Group A and Group B.",
    appearance={"message_box": "neutral"}
  )
  submit_widget_state()
  wf_exe_signal(False)    
  exit()
