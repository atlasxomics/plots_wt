This code is for a notebook platform developed by LatchBio called Plots.  See this link for information on Plots: https://wiki.latch.bio/plots/overview.  The GitHub for Plots can be found here: https://wiki.latch.bio/plots/overview.

Each Plot notebook runs in a virtual machine (AWS EC2).  Users interact with Plots from console.latch.bio via a GUI.  The Plots GUI is organized into tabs.  Each of the top directories in the repository corresponds to a tab in the Plots UI.  Within a tab, the page is subdivided into 'Cells'; in this repository, each file within a directory corresponds to a Cell.  

We organize this notebook with the tab 'welcome' first and 'select_data' last.  welcome/init.py initializes the notebooks, importing packages, setting important global variables, and defining functions.  The remaining Cells in welcome are then run.  The next tabs are then initialized by Plots in the order that they occur.  select_data is initalized last.  Without a data path specified, select_data/select_data.py does not change the state of the other Cells.  Once a data path is specified, select_data/select_data.py reinitalizes all other cells (except welcome/init.py) with the `new_data_signal()` call.

Users can then interact with an navigate the Plots GUI. 

- A notebook is shutdown after use. When it is turned back on, no global variables persist in memory.  However, previously download data exsists on disk.
- Data consists of two h5ad files containing AnnData objects from a DBiT-seq experiemnt, and folders containing bigwigs. The files contain gene accessiblity (*_sm_ge.h5ad) and motif deviation data (*_sm_motifs.h5ad).  The data was generated on the Latch Workflow platform via either of these two repositories: https://github.com/atlasxomics/ATX_snap, https://github.com/atlasxomics/archrproject_latch.
- Plots is run in a mamba env.  Python package versions can be found in the file ~/mamba.txt.

Tabs and cells are initalized in the following order:
- welcome
    - init.py
- compare
    - select_cells.py
    - check_selection.py
    - confirm_selection.py
    - launch_workflow.py
    - fetch_results.py
    - gene_results.py
    - motif_results.py
    - track_browser.py
    - reset_tab.py
- track_browser
    - select_coverage.py
    - track_browser.py
- proportion_plot
    - proportion_plot.py
- violin_plots
    - violin_plot.py
- volcano
    - initalize_display.py
    - compare_plots.py
- heatmap
    - heatmap.py
- bulk_score_genes
    - assign_markers.py
    - score_markers.py
    - h5_viewer.py
    - score_heatmap.py
    - save_data.py
- neighborhood_analysis
    - neighborhood.py
- select_data
    - select_data.py
