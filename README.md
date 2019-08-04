# pyshdom
---
This is a python wrapper for SHDOM (Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer).
The purpose of this wrapper is to develop 3D remote sensing methodologies.
It was created by Aviad Levis, Amit Aides (Technion - Israel Institute of Technology) and Jesse Loveridge (University of Illinois).



The documentation of the source Fortran SHDOM code by Frank K. Evans can be found at http://nit.colorado.edu/shdom/shdomdoc

Information about the source code taken from the documentation page:
This program computes unpolarized monochromatic or spectral band radiative transfer in a one, two, or three-dimensional medium for either collimated solar and/or thermal emission sources of radiation. The properties of the medium can be specified completely generally, i.e. the extinction, single scattering albedo, Legendre coefficients of the scattering phase function, and temperature for the particular wavelength or spectral band may be specified at each input grid point. SHDOM is superior to Monte Carlo radiative transfer methods when many radiative quantities are desired, e.g. the radiance field across the domain top or the 3D distribution of heating. Radiances at any angle, hemispheric fluxes, net fluxes, mean radiances, and net flux convergence (related to heating rates) may be output anywhere in the domain. For highly peaked phase functions the delta-M method may be chosen, in which case the radiance is computed with an untruncated phase function single scattering correction. A correlated k-distribution approach is used for the integration over a spectral band. There may be uniform or spatially variable Lambertian reflection and emission from the ground surface. Several types of bidirectional reflection distribution functions (BRDF) for the surface are implemented, and more may be added easily. SHDOM may be run on a single processor or on multiple processors (e.g. an SMP machine or a cluster) using the Message Passing Interface (MPI).

---
### Updates in pyshdom3.0
 - Code migration to python3
 - Multispectral rendering and optimization
 - Microphysical optimization
 - Main changes to SHDOM core code: 
     1. Removal of global varibles (property array)
     2. Particle mixing is done at runtime and not as an a-proiri computation 
     3. SOLVE_RTE is broken down into initialization and solution iterations
     4. Mie computations are broken down to monodisperse and polydisperse
---
### Installation 
Installation using using anaconda package management

Start a clean virtual enviroment
```
conda create -n pyshdom python=3
source activate pyshdom
```

Install required packages
```
conda install anaconda dill tensorflow tensorboard pillow
```

Install pyshdom distribution with (either install or develop flag)
```
python setup.py install
```
---
### Basic usage
For basic usage follow the following jupyter notebook tutorials
- notebooks/Single Image Rendering.ipynb
- notebooks/Multiview Rendering.ipynb
- notebooks/Multispectral Rendering.ipynb
 
---
### Main scripts
For generating rendering and optimization scripts see the list below. 
The scripts folder contains another readme file with examples of how to run each script.
  - scripts/generate_mie_tables.py
  - scripts/render_monochromatic_radiance_toa.py
  - scripts/render_polychromatic_radiance_toa.py
  - scripts/optimize_extinction.py
  - scripts/optimize_microphysics.py
