# IncrementalDerisFlowVolumeAnalyzer

This repo is for code associated with the IncrementalDebrisFlowVolumeAnalyzer (idfva). This code is fully functional; additional features may be added, and the repo updated, as research using this code is ongoing. All code in the repo is written in Python. Please email laurenmiller@mines.edu with questions.

For analysis the user will need:

- .LAS or .TIF change detection
- A DEM of the AOI
- A flow path of interest

A packaged version of the code is available at https://pypi.org/project/idfva/, with the associated source code on this repo: https://github.com/laurenguido/IncrementalDerisFlowVolumeAnalyzer/tree/main/idfva.

## Installation

You can install **idfva** directly from [PyPI](https://pypi.org/project/idfva/):

```bash
pip install idfva
```
## Additional Documentation

All six idfva modules are provided in the scripts folder, which allows direct access to individual scripts and serves as a built-in test suite. Run the scripts in the following order:

1) 01_Las2Ras.py
2) 02_Process_Watersheds.py
3) 03_Generate_Query_Points.py
4) 04_Prepare_Path.py
5) 05_Apply_Modified_Voronoi.py
6) 06_Estimate_Volume_and_Plot.py

It is recomended that for the most up-to-date version of the package, users download via pip/PyPI. 

The files in this repo have been developed to handle simple (non-branching) as well as complex, multi-strahler-order branching debris flow paths.

Current objectives and ongoing work include standardizing function calls and improving efficiency to reduce the number of scripts. 
