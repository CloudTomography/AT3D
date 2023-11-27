# Atmospheric Tomography with 3D Radiative Transfer (AT3D)

AT3D performs 3D reconstruction of cloud/aerosol microphysical properties from multi-angle, multi-spectral solar reflected radiation using a non-linear optimization procedure [[1],[2],[3]].
The core radiative transfer routines are sourced from the Fortran SHDOM (Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer) code by Frank K. Evans [[4],[5]].

The python package was created by [Aviad Levis](https://www.aviadlevis.info), Amit Aides (Technion - Israel Institute of Technology) and [Jesse Loveridge](https://cloud-radiation.atmos.colostate.edu/jesse-loveridge/) (University of Illinois). Code contributions have been made so far by Linda Forster and Vadim Holodovsky.

## Usage

The AT3D software is built around SHDOM, which is freely distributed online [[4]]. Please contact Frank Evans if you have concerns about the licensing of his code as it appears in this package. This package (AT3D) is distributed under the GNU General Public License (see the `LICENSE` file).
If you want to acknowledge the use of this repository in a publication (e.g. scientific journal article), then please cite the appropriate release, or the most recent release, which is available at the following DOI. See the `CITATION.cff` file for how to reference this repository.

[![DOI](https://zenodo.org/badge/342386439.svg)](https://zenodo.org/badge/latestdoi/342386439)

If you want to acknowledge the scientific origin of a particular feature of this software in a publication, then please cite the appropriate journal or conference articles in which the feature originates [[1],[2],[3]]. This work relies on the generosity of Frank Evans in making his code publicly available, for which we are very grateful. Please acknowledge his work appropriately. In particular, use of the SHDOM solver as a part of the AT3D software in a scientific work should cite [[5]].

Any publications using the synthetic les clouds in the ./data/synthetic_cloud_fields/jpl_les directory which is included in the distribution must cite the following work [[7]].

## Contact

If you find this package useful and/or would like to contribute code please let us know: aviad.levis@gmail.com; Jesse.Loveridge@colostate.edu. 

[1]: http://openaccess.thecvf.com/content_iccv_2015/html/Levis_Airborne_Three-Dimensional_Cloud_ICCV_2015_paper.html
[2]: http://openaccess.thecvf.com/content_cvpr_2017/html/Levis_Multiple-Scattering_Microphysics_Tomography_CVPR_2017_paper.html
[3]: https://doi.org/10.3390/rs12172831
[4]: http://coloradolinux.com/~evans/shdom.html
[5]: https://doi.org/10.1175/1520-0469(1998)055<0429:TSHDOM>2.0.CO;2
[6]: https://doi.org/10.1175/1520-0477(1998)079<0831:OPOAAC>2.0.CO;2
[7]: https://doi.org/10.1175/JAS-D-13-0306.1

&nbsp;

## Features

#### Forward (RTE solver):
AT3D is a python wrapper for [polarized SHDOM](https://coloradolinux.com/~evans/shdom.html) and can be used to compute radiative quantities for a variety of atmospheric configurations.
The key features of polarized SHDOM are included
  1. Solar/Thermal/Combined sources
  2. A variety of (spatially variable) surface BRDFs
  3. Vector or scalar radiative transfer.
  4. Open or periodic boundary conditions.

Note that each RTE solution is serial (**unlike SHDOM**) but independent wavelengths and pixel radiance calculations are parallelized using either MPI or a multi-threading shared memory framework.
Other key features that are implemented are:
  * Several sensor configurations (e.g. Perspective, Orthographic) and arbitrary observation geometries.
  * Mie & Rayleigh scattering optical property calculations including [OPAC aerosols](6). Optical properties of other species (e.g. non-spherical ice or aerosol) can be included but must be calculated externally.
  * Microphysical/optical properties can be generated or be read from netCDF or the SHDOM/[I3RC](https://i3rc.gsfc.nasa.gov/) file format.

#### Inverse (remote-sensing):
AT3D recovers microphysical/optical properties of atmospheric constituents that fit measured radiances.
In contrast to most availble codes, AT3D can recover **3D variable** atmospheric properties. This is achieved by local optimization procedures which employ an approximation to the [Frechet derivatives](https://en.wikipedia.org/wiki/Fr%C3%A9chet_derivative) of the RTE developed by Levis et al. [[3]].

#### Future Improvements
Future improvement include:
* Add additional sensor geometries (cross-track scan, push-broom).
* Parallelize RTE solution with MPI.
* Include retrieval of surface BRDF.

To contribute to the development effort, contact us! see `Contact` section above.

&nbsp;

## Updates in AT3D 4.0
 - Data is represented using xarray objects.
 - Wide field of fiew radiances are now modeled.

&nbsp;

## Installation
Compilation of this package requires Fortran & C compilers (e.g. GCC 9.3.0_1) to be available and correctly linked. Installation has been tested on Mac and Linux using using [anaconda](https://www.anaconda.com/) package management. 

The treatment of legacy Fortran code has changed from GCC 9.X to 10.X+ so currently there is a flag in the setup.py script which needs to be commented if trying to install using GCC 9.X or earlier versions. The flag is `extra_f77_compile_args=["-fallow-argument-mismatch"]`. There is additional discussion of this point written as comments in setup.py but please feel raise an issue or discussion if you run into issues with the compiler version. The default version of the package should compile with GCC v11.3.

Clone the repository into your local machine
```
git clone https://github.com/CloudTomography/AT3D.git
cd AT3D
```

Start a clean virtual environment and setup environment variables
```
conda create -n at3d python=3.10.4
conda activate at3d
```

Install required packages
```
pip install -r requirements.txt
```

Install AT3D distribution. This should be run from within the folder containing setup.py. For development mode add the flag `-e`.
```
pip install .
```

&nbsp;

## Running Tests
After successful installation, run the tests using Python's [nosetests](https://nose.readthedocs.io/en/latest/index.html) package
and make sure they all succeed:
```
cd tests
nose2 -v
```
This command will execute all files starting with *test_\*.py*.
To execute only one specific test file, `test.py` use
```
nose2 -v test
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
