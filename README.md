# pyshdom

Pyshdom performs 3D reconstruction of cloud/aerosol microphysical properties from multi-angle, multi-spectral solar reflected radiation using a non-linear optimization procedure [[1],[2],[3]]. The core radiative transfer routines are sourced from the Fortran SHDOM (Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer) code by Frank K. Evans [[4]]. The python package was created by Aviad Levis, Amit Aides (Technion - Israel Institute of Technology) and Jesse Loveridge (University of Illinois).

[1]: http://openaccess.thecvf.com/content_iccv_2015/html/Levis_Airborne_Three-Dimensional_Cloud_ICCV_2015_paper.html
[2]: http://openaccess.thecvf.com/content_cvpr_2017/html/Levis_Multiple-Scattering_Microphysics_Tomography_CVPR_2017_paper.html
[3]: https://www.mdpi.com/2072-4292/12/17/2831
[4]: http://coloradolinux.com/~evans/shdom.html

&nbsp;

## Features

* The key features of polarized SHDOM are included
  1. Solar/Thermal/Combined sources
  2. A variety of (spatially variable) surface BRDFs
  3. Vector or scalar radiative transfer.
  4. Open or periodic boundary conditions.

* Local optimization procedures are included for recovery of the microphysical or optical properties of atmospheric constituents.
  1. The linearization used in the optimization employs an approximation to the Frechet derivatives of the RTE developed by Levis et al. [[3]].
* Each RTE solution is serial (**unlike SHDOM**) but independent wavelengths and pixel radiance calculations are parallelized using either MPI or a multi-threading shared memory framework.
* Wide field-of-view radiances can be calculated and arbitrary observation geometries are supported. Defaults for both Perspective and Orthographic sensor geometries are included.
* Mie & Rayleigh scattering optical property calculations. Optical properties of other species (e.g. non-spherical ice or aerosol or absorbing gases) can be included but must be calculated externally.
* Microphysical/optical properties can be generated or be read from netCDF or the SHDOM/i3rc file format.

### Future Improvements:

* Add additional sensor geometries (cross-track scan, push-broom).
* Parallelize RTE solution with MPI.
* Include retrieval of surface BRDF.


&nbsp;

## Updates in pyshdom 4.0
 - Data is represented using xarray objects.
 - Wide field of fiew radiances are now modeled.
 - Multiple 

&nbsp;

## Installation
Compilation of this package requires Fortran & C compilers (e.g. GCC 9.3.0_1) to be available and correctly linked. Installation has been tested on Mac and Linux.

Installation using using anaconda package management

Start a clean virtual environment
```
conda create -n pyshdom python=3
source activate pyshdom
```

Install [xarray](http://xarray.pydata.org/) its dependencies 
```
conda install -c conda-forge xarray dask netCDF4 bottleneck
```

Install other required packages
```
pip install -r requirements.txt
```

Install pyshdom distribution. This should be run from within the folder containing setup.py.
```
pip install .
```

&nbsp;

## Running Tests
After successful installation, run the tests using Python's [nosetests](https://nose.readthedocs.io/en/latest/index.html) package
and make sure they all succeed:
```
cd tests
nosetests
```
This command will execute all files starting with *test_\*.py*.
To execute only one specific test file, use
```
nosetests test_XYZ.py
```

&nbsp;

## Basic usage
For basic usage follow the following jupyter notebook tutorials
TODO

&nbsp;

## Main scripts
For generating rendering and optimization scripts see the list below.
The scripts folder contains another readme file with examples of how to run each script.
TODO


&nbsp;

## Usage and Contact
If you find this package useful please let me know at aviad.levis@gmail.com, I am interested.
If you use this package in an academic publication please acknowledge the appropriate publications (see LICENSE file).
