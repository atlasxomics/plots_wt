new_data_signal()

w_text_output(content="""## Compare Workflow Results""")

if not adata_g:
  w_text_output(content="   ")
  exit()

# Checkbox to load in previous compare clusters data
load_compare_box = w_checkbox(
    label="Load in previous compare Workflow data",
    default=False,
    appearance={"description": "Select a compare_outs folder from Latch Data to plot."}
)

if not load_compare_box.value:
  
  wf_exe_signal()
  if wf_exe_signal.sample() == True:

    get_results = w_button(label="Fetch Workflow Results")

    w_text_output(content="Click to load in Execution results.")
    
    
    if get_results.value:
      if execution is not None:
  
        list_resp = post(
            url=config.api.execution.list,
            headers = {"Authorization": get_auth_header()},
            json={"ws_account_id": f"{workspace_account_id}"},
        ).json()
        target_execution = list_resp[str(execution.id)]
        status = target_execution.get('status')
      
        if status == 'SUCCEEDED':
  
            w_text_output(
                content=f"Compare Clusters successfully finished running; loading results...",
                appearance={"message_box": "success"}
            )
            submit_widget_state()
            
            res = LPath(f"latch:///compare_outs/{wf_name.value}")

            # Check if dir exists
            try:
              res.is_dir()
            except LatchPathError:
              w_text_output(
                content=f"Cannot find compare Workflow outputs, please check Execution logs...",
                appearance={"message_box": "success"}
              )
              submit_widget_state()
              exit()

            res_contents = {p.name(): p for p in res.iterdir()}
            
    
            for feat in feats:
              feat_dir = f"{feat}_results"
              
              if feat_dir in res_contents.keys():
                feat_csv = LPath(f"{res_contents[feat_dir].path}/all_{feat}s.csv")
                
                try:
                    feat_csv.download(Path(f"all_{feat}s.csv"), cache=True)
                    results_dict[feat] = pd.read_csv(f"all_{feat}s.csv")
                except ValueError:
                  w_text_output(
                    content=f"Could not find results csv for {feat}s in {feat_dir}.",
                    appearance={"message_box": "warning"}
                  )
                  submit_widget_state()
        
              else:
                w_text_output(
                  content=f"Could not find results for {feat}s in {res.name}.",
                  appearance={"message_box": "warning"}
                )
                submit_widget_state()
    
            if len(results_dict.keys()) > 0:
              w_text_output(
                content=f"Loaded {', '.join(list(results_dict.keys()))} results",
                appearance={"message_box": "success"}
              )
              submit_widget_state()
              wf_results_signal(True)
            else:
              w_text_output(
                content="Could not find output tables for Execution; please check Execution logs.",
                appearance={"message_box": "warning"}
              )
              submit_widget_state()
              wf_results_signal(False)
  
            # Coverages
            coverages_dir = LatchDir(f"{res.path}/coverages")
  
            files = []
            try:
              for file in coverages_dir.iterdir():
                  suffix = file.path.split(".")[-1]
                  if suffix == "bw":
                      files.append(file)
              if len(files) > 0:
                wf_bigwigs_signal(True)
                w_text_output(
                  content="Loaded coverage bigwigs for Group A and Group B...",
                  appearance={"message_box": "success"}
                )
                submit_widget_state()
              else:
                w_text_output(
                  content="No coverage files found...",
                  appearance={"message_box": "warning"}
                )
                submit_widget_state()
            except ValueError:
                w_text_output(
                  content="Could not find remote coverages folder; please check Execution logs.",
                  appearance={"message_box": "warning"}
                )
                submit_widget_state()
  
        elif status in ["ABORTED", "FAILED"]:
          w_text_output(
            content="Workflow Execution failed or aborted, please check Execution logs.",
            appearance={"message_box": "warning"}
          )
          submit_widget_state()
          wf_results_signal(False)
          exit()
    
        elif execution.status in ["UNDEFINED", "RUNNING"]:
          w_text_output(
            content="Workflow still running, click button again once it has completed.",
            appearance={"message_box": "warning"}
          )
          submit_widget_state()
          wf_results_signal(False)
          exit()
        else:
          w_text_output(
            content="Unknown Execution status...",
            appearance={"message_box": "neutral"}
          )
          submit_widget_state()
          
      else:
        w_text_output(
          content="""Awaiting Workflow launch...""",
          appearance={"message_box": "neutral"}
        )
        submit_widget_state()
        wf_results_signal(False)
  
  else:
      w_text_output(
        content="""Awaiting Workflow launch...""",
        appearance={"message_box": "neutral"}
      )
      submit_widget_state()
      wf_results_signal(False)
else: 
  compare_path = w_ldata_picker(
    label="compare_outs folder",
    key="compare_path",
    appearance={
      "placeholder": "Select output directory from Latch Data `/compare_outs/`."
    }
  )
  compare_genome = w_select(
    label="Reference Geneome",
    default=None,
    key="compare_genome",
    options=tuple(["hg38", "mm10", "rn6"]),
    appearance={
      "help_text": "Reference genome for experiment."
    }
  )

  input_compare_row = w_row(items=[compare_path, compare_genome])

  load_compare_button = w_button(label="Load Data")
  
  if load_compare_button.value:
    if compare_genome.value is not None:
      if compare_path.value is not None:
        res = compare_path.value
    
        if not res.is_dir():
          w_text_output(
            content="Selected resource is not a directory; please select a directory.",
            appearance={"message_box": "warning"}
          )
          submit_widget_state()
          exit()
    
        res_contents = {p.name(): p for p in res.iterdir()}
        
        for feat in feats:
          feat_dir = f"{feat}_results"
          
          if feat_dir in res_contents.keys():
            feat_csv = LPath(f"{res_contents[feat_dir].path}/all_{feat}s.csv")
            
            try:
                feat_csv.download(Path(f"all_{feat}s.csv"), cache=True)
                results_dict[feat] = pd.read_csv(f"all_{feat}s.csv")
            except ValueError:
              w_text_output(
                content=f"Could not find results csv for {feat}s in {feat_dir}.",
                appearance={"message_box": "warning"}
              )
              submit_widget_state()
    
          else:
            w_text_output(
              content=f"Could not find results for {feat}s in {res.path}.",
              appearance={"message_box": "warning"}
            )
            submit_widget_state()
      
        if len(results_dict.keys()) > 0:
            w_text_output(
              content=f"Loaded {', '.join(list(results_dict.keys()))} results",
              appearance={"message_box": "success"}
            )
            submit_widget_state()
            wf_results_signal(True)
        else:
          w_text_output(
            content="Could not find expected outputs in selected folder.",
            appearance={"message_box": "warning"}
          )
          submit_widget_state()
          wf_results_signal(False)
      
        # Coverages
        coverages_dir = LatchDir(f"{res.path}/coverages")
      
        files = []
        try:
          for file in coverages_dir.iterdir():
              suffix = file.path.split(".")[-1]
              if suffix == "bw":
                  files.append(file)
          if len(files) > 0:
            wf_bigwigs_signal(True)
            w_text_output(
              content="Loaded coverage bigwigs for Group A and Group B...",
              appearance={"message_box": "success"}
            )
            submit_widget_state()
          else:
            w_text_output(
              content="No coverage files found in selected folder.",
              appearance={"message_box": "warning"}
            )
            submit_widget_state()
        except ValueError:
            w_text_output(
              content="Could not find remote coverages folder; please check Execution logs.",
              appearance={"message_box": "warning"}
            )
            submit_widget_state()
    else:
      w_text_output(
        content="Please select reference genome.",
        appearance={"message_box": "warning"}
      )
      submit_widget_state()
      exit()
