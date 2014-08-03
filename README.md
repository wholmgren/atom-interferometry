Introduction
============
The code contained in this repository was used during the course of my Ph.D. research in Prof. Alex Cronin's atom interferometry lab at the University of Arizona. The research culiminated in the following publications:
* Holmgren et. al., *Measurement of a Wavelength of Light for Which the Energy Shift for an Atom Vanishes*, Physical Review Letters **109**, 243004 (2012). [link](http://dx.doi.org/10.1103/PhysRevLett.109.243004)
* Holmgren et. al., *Atom beam velocity measurements using phase choppers*, New Journal of Physics, New Journal of Physics **13**, 115007 (2011). [link](http://dx.doi.org/10.1088/1367-2630/13/11/115007)
* Holmgren et. al., *Absolute and ratio measurements of the polarizability of Na, K, and Rb with an atom interferometer*, Physical Review A **81** 053607 (2010). [link](http://dx.doi.org/10.1103/PhysRevA.81.053607)
* Holmgren, *Polarizability and magic-zero wavelength measurements of alkali atoms*, Ph.D. thesis, University of Arizona (2013). [link](http://www.atomwave.org/otherarticles/Holmgren%202013%20thesis.pdf)

I did not use a good version control system during my Ph.D. research, so I cannot promise that the code in this repository is the exact code that led to the publications above. I hope that future students and researchers do not repeat my mistake.


The code
========
99% of this code is for Igor Pro.

The most developed and important procedures are in the 'User Procedures' folder.
The most important procedures are listed here:

``Fringe Fit.ipf``: fits the interferometer data.

``Fit Diffraction with Convolution.ipf``: fit diffraction data to find atom beam velocity and velocity distribution.

``TOW procedures.ipf``: used for analyzing the TOW/MZW/lambda_zero experiments.

``Pol2Electrodes.ipf``: used for analyzing static polarizability experiments. All other 'Pol' files are older, but may contain useful stuff still.
    
``Phase Choppers.ipf``: my version of the phase choppers code. Superseded by Ivan and Maxwell's work, but there are still useful helper functions here.

``PhysicalConstants.ipf``: eps0, masses, conversions, polarizabilities, etc

``LabConstants.ipf``: Not much at this point, but I recommend a common file for this type of data.
