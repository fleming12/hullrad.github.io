# hullrad.github.io
HullRad is an algorithm for calculating hydrodynamic properties of a macromolecule from a structure file.
A webserver to run HullRad on your macromolecule can be found at (http://52.14.70.9/).

To run HullRad on your local machine as a Python script, download the file HullRadV9.py and read the header for instructions on installing the necessary libraries. 

HullRadSAS is similar to the original HullRad but it also provides the amount of two types of macromolecular hydration, shell and entrained. Version 1.1 fixes a bug that overestimated Dmax and affected rotational diffusion calculations.

If you publish work that uses HullRad please cite one of the following references: 

Fleming, P.J. and Fleming, K. G., HullRad: Fast Calculations of Folded and Disordered Protein and Nucleic Acid Hydrodynamic Properties. Biophysical Journal, Volume 114, Issue 4, p856–869, 27 February 2018. 

Fleming, P.J., Correia, J.J, and Fleming, K.G., Revisiting Macromolecular Hydration with HullRadSAS.
(2023) Eur. Biophys. J. Online Jan 5. https://doi.org/10.1007/s00249-022-01627-8

Note: The EBJ paper also describes updates to the original HullRad incorporated in HullRadV9.

Fleming, P.J., Correia, J.J, and Fleming, K.G., The Molecular Basis for Hydrodynamic Properties of PEGylated Human Serum Albumin.
Biophysical Journal, Online June 21, 2024  DOI: 10.1016/j.bpj.2024.05.019

Download fleming_bj_2018.pdf for a PDF copy of the paper describing the original HullRad algorithm.

Download the Python script, Display_hull.py (Python2) or Display_Hull3.py (Python3), to draw convex hull around a displayed structure in PyMOL.

Other scripts to run HullRad in batch mode are available at the HullRad website.
