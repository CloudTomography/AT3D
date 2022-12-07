"""
Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer

This is a python wrapper for SHDOM created by Aviad Levis and Amit Aides, Technion Inst. of Technology
and Jesse Loveridge, University of Illinois at Urbana-Champaign.
The purpose of this wrapper is to develop 3D remote sensing metoglogies.

The documentation of the source Fortran code by Frank Evans can be found at
https://nit.coloradolinux.com/~evans/shdom.html

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
"""
import at3d.aerosol
import at3d.rayleigh
import at3d.medium
import at3d.mie
import at3d.size_distribution
import at3d.surface
import at3d.source
import at3d.configuration
import at3d.solver
import at3d.sensor
import at3d.grid
import at3d.gas_absorption
import at3d.checks
import at3d.gradient
import at3d.optimize
import at3d.space_carve
import at3d.stereo
import at3d.util
import at3d.uncertainties
import at3d.callback
import at3d.containers
import at3d.parallel
import at3d.exceptions
import at3d.regularization
import at3d.generate
import at3d.transforms
import at3d.initialization
import at3d.preprocessing
import at3d.visualization
import at3d.initialization_workflow
