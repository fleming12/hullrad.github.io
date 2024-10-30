# hullrad.github.io
HullRad is an algorithm for calculating hydrodynamic properties of a macromolecule from a structure file.
A webserver to run HullRad on your macromolecule can be found at (http://52.14.70.9/).

To run HullRad on your local machine as a Python script, download the file HullRadV10.1.py and read the header for instructions on installing the necessary libraries. 

HullRadSAS is similar to the original HullRad but it also provides the amount of two types of macromolecular hydration, shell and entrained. The current version is HullRadSASV3.1.py.

If you publish work that uses HullRad please cite one of the following references: 

Fleming, P.J. and Fleming, K. G., HullRad: Fast Calculations of Folded and Disordered Protein and Nucleic Acid Hydrodynamic Properties. Biophysical Journal, Volume 114, Issue 4, p856–869, 27 February 2018. 

Fleming, P.J., Correia, J.J, and Fleming, K.G., Revisiting Macromolecular Hydration with HullRadSAS.
(2023) Eur. Biophys. J. Online Jan 5. https://doi.org/10.1007/s00249-022-01627-8

Note: The EBJ paper also describes updates to the original HullRad incorporated in HullRadV9.

Fleming, P.J., Correia, J.J, and Fleming, K.G., The Molecular Basis for Hydrodynamic Properties of PEGylated Human Serum Albumin.
Biophysical Journal, Online June 21, 2024  DOI: 10.1016/j.bpj.2024.05.019

Note: The above BJ paper describes the addition of non-ideality constants to the output of HullRad.

Download fleming_bj_2018.pdf for a PDF copy of the paper describing the original HullRad algorithm.

Download the Python script, Display_hull.py (Python2) or Display_Hull3.py (Python3), to draw convex hull around a displayed structure in PyMOL.

Other scripts to run HullRad in batch mode are available at the HullRad website.

Version 5 (The program as described in Fleming and Fleming, BJ, 2018)
Proteins
Nucleic Acids

Version 6 (November, 2019)
Includes many saccharides for glycosylation
Rg is calculated from all atom model, not reduced model as in version 5
Includes option to use numpy and scipy
Asphericity is calculated from gyration tensor (requires numpy)
Option to use either qconvex from qhull or ConvexHull from scipy
       (Many thanks to Chad Brautigam for implementation)

Version 7 (July, 2020)
Ported to Python 3 (Python 2 edition of version 7 also available)
Addition of chain ID
Includes several detergents for analysis of protein/detergent micelles.
        (Note that hydrodynamic calculations of detergents
       is sensitive to partial specific volume of the detergent micelle.
       Literature values of detergent partial specific volumes are not consistent.
       Solution conditions may also affect partial specific volumes - read the notes in the
       code and use with caution.)
Calculation of intrinsic viscosity was removed.
       The original implementation was appropriate only for very sperical particles.
       This parameter is sensitive to the model used for axial ratio determination and needs further development.

Version 8 (November, 2021)
Accepts mmCIF file format as well as PDB format
Additional error checking for missing backbone atoms

Version 8.1 (March, 2022)
Fixed bugs introduced with mmCIF reader

Version 9 (October, 2022)
Addition of Intrinsic Viscosity
Addition of Non-idealtiy Constant

Version 10 (July, 2024)
Update Non-idealtiy Constant equation for ks (See DOI:https://doi.org/10.1016/j.bpj.2024.05.019)
Note that the sedimentation non-ideality constant (ks) is now calculated using a modified version of the Rowe equation and gives a different (corrected) value compared to HullRad Version 9
Add kd and second virial coeficient

Version 10.1 (Oct, 2024)
Update nucleic acid partial specific volumes.
Many thanks to Brad Chaires and Rob Monsen for pointing out that experimentally the partial specific volumes of all nucleic acids is approximately 0.55.
This is different from that calculated from structure using Voronoi volumes, which HullRad previously used.
Curiously, for proteins the experimentally determined partial specific volumes agree pretty well with those calculated from strucucture using Voronoi volumes. Why this doesn't work for nucleic acids is a mystery.
