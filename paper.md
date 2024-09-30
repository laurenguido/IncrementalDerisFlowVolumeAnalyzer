title: 'IncrementalDebrisFlowVolumeAnalyzer: A Python tool for estimating volumes along a debris flow path'
tags:
  - Python
  - geohazards
  - debris flow
  - erosion and deposition
  - geomorphology
authors:
  - name: Lauren E. Guido
    orcid: 0000-0003-4449-560X
    equal-contrib: true
    affiliation: 1
  - name: Paul M. Santi
  - equal-contrib: true
  - affiliation: 1

affiliations:
 - name: Colorado School of Mines, Department of Geology and Geological Engineering, Golden, CO, USA
   index: 1

date: 29 September 2024
bibliography: paper.bib

# Summary

Debris flows are common, costly, and deadly post wildfire hazards. As wildfire activity 
is expected to increase in the coming decades, a precipitous increase post-fire debris 
flow hazards is expected. The tool presented here allows for the semi-automated 
estimation of erosional and depositional volumes along the path of a debris flow for use
in ongoing research which aims to better understand the intra-channel mechanisms of 
bulking and deposition. Understanding these field-scale mechanics will play an important
role in better focusing mitigation efforts, improving modeling for runout and hazard 
prediction, and predicting locations of volume growth or decline. 

# Statement of need

The `IncrementalDebrisFlowVolumeAnalyzer` is a Python tool for estimating intra-channel 
volume of erosion and deposition along a flow path. The series of semi-automated 
scripts which make up `IncrementalDebrisFlowVolumeAnalyzer` can be looped using a master 
text file and/or the `glob` package depending on the user file structure, to efficiently 
run volume estimations across an area of interest. The `IncrementalDebrisFlowVolumeAnalyzer`
relies heavily on the `Fiona`, `Rasterio`, and `Shapely` packages for reading and writing 
geospatial data. 

The `IncrementalDebrisFlowVolumeAnalyzer` was designed to be used by geohazard and 
geomorphology researchers in tandem with external data, field work, and 
geomorphometric analyses. Students may also use this tool to gain familiarity with 
debris flow hazards, change detection data types, and geospatial data manipulation
in Python with hands-on application to real world problems. The 
`IncrementalDebrisFlowVolumeAnalyzer` has been used in multiple forthcoming 
scientific publications, and has been previously presented on [@Guido:2024]. This work 
builds off the concepts of volume estimation demonstrated by [@Scheip & Wegmann:2023], 
with improvemnts in standardized sample areas, ability to handle complex tributary 
systems, referencing to hydrological locations in the watershed, and increasing 
the capacity for automation. These improvements, combined with the growing need
to investigate intra-channel debris flow dynamics and increasing availability of
high-resolution topographic data, will enable exciting scientific exploration and 
developemnts in mitigation, modeling, and hazard prediction in debris flow and 
geohazard research. 


# Graphical Summary 



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
