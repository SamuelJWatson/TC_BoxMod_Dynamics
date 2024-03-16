## Low-Order Tropical Cyclone Model Dynamics
This repository contains MATLAB code and data structures for analysis the dynamics of a low-order tropical cyclone model developed by Schönemann and Frisius [[1]](#1). Each folder is self contained. 

Bifurcation analysis using the Continuation Core and Toolboxes (COCO) is provided in `Curve_Tracking`. 
COCO is avaliable at https://sourceforge.net/projects/cocotools/. 
`TCCurveTracking_beta.m` runs a bifurcation analysis in $\beta$ and `TCCurveTracking_SST.m` in SST. 

Basin Instability analysis is provided in `Basin_Instability_Analysis` and run with `BasinInstability.m`.

Time dependent sech forcing of $\beta$ and SST is provided in `Parameter_Forcing` and run with `S_TCModRun.m`.


# References
<a id="1">[1]</a> D. Schönemann and T. Frisius, “Dynamical system analysis of a low-order tropical cyclone model,” Tellus A:305 Dynamic Meteorology Oceanography 64, 15817 (2012).
