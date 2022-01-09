# pyshdom

Pyshdom performs 3D reconstruction of cloud/aerosol microphysical properties from multi-angle, multi-spectral solar reflected radiation using a non-linear optimization procedure [[1],[2],[3]].
The core radiative transfer routines are sourced from the Fortran SHDOM (Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer) code by Frank K. Evans [[4]].
The python package was created by [Aviad Levis](https://www.aviadlevis.com/3d-remote-sensing), Amit Aides (Technion - Israel Institute of Technology) and Jesse Loveridge (University of Illinois). Code contribution were made by Linda Forster and Vadim Holodovsky.

[1]: http://openaccess.thecvf.com/content_iccv_2015/html/Levis_Airborne_Three-Dimensional_Cloud_ICCV_2015_paper.html
[2]: http://openaccess.thecvf.com/content_cvpr_2017/html/Levis_Multiple-Scattering_Microphysics_Tomography_CVPR_2017_paper.html
[3]: https://www.mdpi.com/2072-4292/12/17/2831
[4]: http://coloradolinux.com/~evans/shdom.html

&nbsp;

## Features

#### Forward (RTE solver):
Pyshdom is a python wrapper for [polarized SHDOM](https://coloradolinux.com/~evans/shdom.html) and can be used to compute radiative quantities for a variety of atmospheric configurations.
The key features of polarized SHDOM are included
  1. Solar/Thermal/Combined sources
  2. A variety of (spatially variable) surface BRDFs
  3. Vector or scalar radiative transfer.
  4. Open or periodic boundary conditions.

Note that each RTE solution is serial (**unlike SHDOM**) but independent wavelengths and pixel radiance calculations are parallelized using either MPI or a multi-threading shared memory framework.
Other key features that are implemented are:
  * Several sensor configurations (e.g. Perspective, Orthographic) and arbitrary observation geometries.
  * Mie & Rayleigh scattering optical property calculations. Optical properties of other species (e.g. non-spherical ice or aerosol or absorbing gases) can be included but must be calculated externally.
  * Microphysical/optical properties can be generated or be read from netCDF or the SHDOM/[I3RC](https://i3rc.gsfc.nasa.gov/) file format.

#### Inverse (remote-sensing):
Pyshdom recovers microphysical/optical properties of atmospheric constituents that fit measured radiances.
In contrast to most availble codes, pyshdom can recover **3D variable** atmospheric properties. This is achieved by local optimization procedures which employ an approximation to the [Frechet derivatives](https://en.wikipedia.org/wiki/Fr%C3%A9chet_derivative) of the RTE developed by Levis et al. [[3]].

#### Future Improvements
Future improvement include:
* Add additional sensor geometries (cross-track scan, push-broom).
* Parallelize RTE solution with MPI.
* Include retrieval of surface BRDF.

To contribute to the development effort, contact us! see `Usage and Contact` section below.

&nbsp;

## Updates in pyshdom 4.0
 - Data is represented using xarray objects.
 - Wide field of fiew radiances are now modeled.

&nbsp;

## Installation
Compilation of this package requires Fortran & C compilers (e.g. GCC 9.3.0_1) to be available and correctly linked. Installation has been tested on Mac and Linux using using [anaconda](https://www.anaconda.com/) package management.

Clone the repository into your local machine
```
git clone https://github.com/CloudTomography/pyshdom.git
cd pyshdom
```

Start a clean virtual environment and setup enviroment variables
```
conda create -n pyshdom python=3
conda activate pyshdom
conda env config vars set PYSHDOM_DIR=$(pwd)
conda activate pyshdom
```

Install [xarray](http://xarray.pydata.org/) and its dependencies
```
conda install -c conda-forge xarray dask netCDF4 bottleneck
```

Install other required packages
```
pip install -r requirements.txt
```

Install pyshdom distribution. This should be run from within the folder containing setup.py. For development mode add the flag `-e`.
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
For basic usage follow the following jupyter notebook tutorials under the notebooks directory:

* MakeMieTablesExample
* MakePolydisperseMie
* MakeOpticalProperties
* MakeSensors
* SolveRTE
* SimulatingRadiances
* SimpleInverseProblem

&nbsp;

## Main scripts
For generating rendering and optimization scripts see the list below.
The scripts folder contains another readme file with examples of how to run each script.
TODO


&nbsp;

## Usage and Contact
If you find this package useful and/or would like to contribute code please let us know: aviad.levis@gmail.com; jesserl2@illinois.edu.  
If you use this package in an academic publication please acknowledge the appropriate publications (see LICENSE file).
