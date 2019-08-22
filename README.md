# pyshdom
---
Pyshdom performs 3D reconstruction of cloud microphysical properties from multi-angule, multi-spectral solar reflected radiation using a non-linear optimization procedure [[1],[2]]. The core radiative transfer routines are sourced from the Fortran SHDOM (Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer) code by Frank K. Evans [[3]]. The python package was created by Aviad Levis, Amit Aides (Technion - Israel Institute of Technology) and Jesse Loveridge (University of Illinois).

[1]: http://openaccess.thecvf.com/content_iccv_2015/html/Levis_Airborne_Three-Dimensional_Cloud_ICCV_2015_paper.html
[2]: http://openaccess.thecvf.com/content_cvpr_2017/html/Levis_Multiple-Scattering_Microphysics_Tomography_CVPR_2017_paper.html
[3]: http://coloradolinux.com/~evans/shdom.html

---

At present pyshdom has the following features:

* Mie & Rayleigh scattering optical property calculations. Optical properties of other species (e.g. non-spherical ice or aerosol) can be included but must be calculated externally.
* Cloud data of varying complexity can be generated or read from LES output.
* Scalar radiative transfer from SHDOM with Perspective or Orthographic sensor geometries and Lambertian Surface.
* Each SHDOM solution is serial but independent wavelengths and pixel radiance calculations are
parallelised in a shared memory framework.
* Local & Global optimization procedures for recovery of cloud microphysical properties (liquid water content,
droplet effective radius, droplet effective variance) on 3D (or reduced order) grids from simulated (or observed) radiances.
* The calculation of optical properties and each SHDOM solution have been streamlined to minimize computational resources
necessary at each iteration of the optimization routines.

Future Improvements:

* Implement vector SHDOM to include polarisation information in the optimization.
* Add additional sensor geometries (cross-track scan, push-broom) & expand to other surface types in SHDOM.
* Include a more flexible parallelisation scheme.
* Add useful regularisation options in the optimization procedure.
* Add further accelerations for computational efficiency.
* Include gaseous absorption for greater realism.
---

The documentation of the source Fortran SHDOM code by Frank K. Evans can be found at http://nit.colorado.edu/shdom/shdomdoc

Information about the source code taken from the documentation page:
This program computes unpolarized monochromatic or spectral band radiative transfer in a one, two, or three-dimensional medium for either collimated solar and/or thermal emission sources of radiation. The properties of the medium can be specified completely generally, i.e. the extinction, single scattering albedo, Legendre coefficients of the scattering phase function, and temperature for the particular wavelength or spectral band may be specified at each input grid point. SHDOM is superior to Monte Carlo radiative transfer methods when many radiative quantities are desired, e.g. the radiance field across the domain top or the 3D distribution of heating. Radiances at any angle, hemispheric fluxes, net fluxes, mean radiances, and net flux convergence (related to heating rates) may be output anywhere in the domain. For highly peaked phase functions the delta-M method may be chosen, in which case the radiance is computed with an untruncated phase function single scattering correction. A correlated k-distribution approach is used for the integration over a spectral band. There may be uniform or spatially variable Lambertian reflection and emission from the ground surface. Several types of bidirectional reflection distribution functions (BRDF) for the surface are implemented, and more may be added easily. SHDOM may be run on a single processor or on multiple processors (e.g. an SMP machine or a cluster) using the Message Passing Interface (MPI).

---
### Updates in pyshdom3.0
 - Code migration to python3
 - Multispectral rendering and optimization
 - Microphysical optimization
 - Main changes to SHDOM core code: 
     1. Removal of global variables (property array)
     2. Particle mixing is done at runtime and not as an a-priori computation 
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
