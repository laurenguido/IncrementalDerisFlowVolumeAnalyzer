# IncrementalDerisFlowVolumeAnalyzer

This repo is for code associated with the IncrementalDebrisFlowVolumeAnalyzer. This code is fully functional. Additional features may be added, and the repo updated, as dissertation work using this code is ongoing. All code in the repo is written in Python, and developed in Spyder. Please email laurenmiller@mines.edu with questions.

For analysis the user will need:

- .LAS or .TIF change detection
- A DEM of the AOI
- A flow path of interest

The processes for the IncrementalDebrisFlowVolumeAnalyzer have been subdivided into 6 scripts to help with run times, trouble shooting, and batch processing. These should be run in the following order:

1) 01_Las2Ras.py
2) 02_Process_Watersheds.py
3) 03_Generate_Query_Points.py
4) 04_Prepare_Path.py
5) 05_Apply_Modified_Voronoi.py
6) 06_Estimate_Volume_and_Plot.py

The files in this repo have been developed to handle simple (non-branching) as well as complex, multi-strahler-order branching debris flow paths.

Current objectives and ongoing work include standardizing function calls and improving efficiency to reduce the number of scripts. 
