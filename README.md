# PyEmittance

PyEmittance is a tool for the adaptive measurement of beam emittance (*e.g.* of an electron beam at OTR/YAG beam profile monitors) using a single quadrupole scan developed at SLAC National Accelerator Laboratory. The Twiss parameters and the 'Bmag' match parameter can also be obtained. It can be used in different beamlines/machines by defining machine-specific configuration files. 

This tool was designed for robustness during online machine learning optimizations where each measurement needs to be reliable. It is still a work in progress, with version 1.0 to be released soon.


Clone the repository and run the command below to install: 

    pip install -e .

Or install by running:
    
    pip install pyemittance
    
