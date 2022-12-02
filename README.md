# PyEmittance

PyEmittance is a tool for the adaptive measurement of beam emittance (*e.g.* of an electron beam at OTR/YAG beam profile monitors) using a single quadrupole scan developed at SLAC National Accelerator Laboratory. The Twiss parameters and the 'Bmag' match parameter can also be obtained. It can be used in different beamlines/machines by defining machine-specific configuration files. 

This tool was designed for robustness during online machine learning optimizations where each measurement needs to be reliable. It is still a work in progress, with version 1.0 to be released soon.

For more information, please see the IPAC coneference proceedings [DOI: 10.18429/JACoW-IPAC2022-TUPOST059](https://accelconf.web.cern.ch/ipac2022/doi/JACoW-IPAC2022-TUPOST059.html)

## Installation

Clone the repository and run the command below to install: 

    pip install -e .

Or install by running:
    
    pip install pyemittance
    
## Setting the correct configuration parameters for the machine/simulation

For examples on how to set the correct configuration and how to use the new API, please see [`examples/tutorial.ipynb`](https://github.com/pluflou/PyEmittance/blob/main/examples/tutorial.ipynb).

## Citations
For citing this work, please use:
```
{@article{Miskovich:2022kqg,
    author = "Miskovich, Sara and Edelen, Auralee and Mayes, Christopher",
    title = "{PyEmittance: A General Python Package for Particle Beam Emittance Measurements with Adaptive Quadrupole Scans}",
    doi = "10.18429/JACoW-IPAC2022-TUPOST059",
    journal = "JACoW",
    volume = "IPAC2022",
    pages = "TUPOST059",
    year = "2022"
}
```
