# pyshdom

This is a python wrapper for SHDOM created by Aviad Levis and Amit Aides, Technion Inst. of Technology.
The purpose of this wrapper is to develop 3D remote sensing metoglogies. 

SHDOM: Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer

The documentation of the source Fortran code by Frank Evans can be found at
http://nit.colorado.edu/shdom/shdomdoc

Information about the source code taken from the documentation page:
This program computes unpolarized monochromatic or spectral band radiative transfer in a one, two,
or three-dimensional medium for either collimated solar and/or thermal emission sources of radiation.
The properties of the medium can be specified completely generally, i.e. the extinction, single 
scattering albedo, Legendre coefficients of the scattering phase function, and temperature for
the particular wavelength or spectral band may be specified at each input grid point. SHDOM is
superior to Monte Carlo radiative transfer methods when many radiative quantities are desired,
e.g. the radiance field across the domain top or the 3D distribution of heating. Radiances at
any angle, hemispheric fluxes, net fluxes, mean radiances, and net flux convergence (related
to heating rates) may be output anywhere in the domain. For highly peaked phase functions the 
delta-M method may be chosen, in which case the radiance is computed with an untruncated phase
function single scattering correction. A correlated k-distribution approach is used for the
integration over a spectral band. There may be uniform or spatially variable Lambertian
reflection and emission from the ground surface. Several types of bidirectional reflection
distribution functions (BRDF) for the surface are implemented, and more may be added easily.
SHDOM may be run on a single processor or on multiple processors (e.g. an SMP machine or a
cluster) using the Message Passing Interface (MPI).


## Requirements (python 2.7):
(required) numpy, scipy, matplotlib, dill, joblib, opencv

(optional - summary and display) pillow, tensorflow, tensorboard, tensorboardX

## Installation (using anaconda package management)

Start a clean virtual enviroment
```
conda create -n pyshdom
source activate pyshdom
```

Install required packages
```
conda install numpy scipy matplotlib dill joblib opencv
```

Install jupyter notebook and extentions
```
conda install jupyter
```

Install optional packages
```
conda install tensorflow tensorboard pillow
pip install tensorboardX
```

Install pyshdom distribution with (either install or develop flag)
```
python setup.py install
```

## Basic usage
For basic usage follow the jupyter notebook tutorials
 - notebooks/Make Mie Table.ipynb
 - notebooks/Forward Rendering.ipynb 

## Main scripts
For rendering and optimization scripts see
  - scripts/render_radiance_toa.py
  - scripts/optimize_extinction.py
