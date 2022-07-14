"""This module contains functions to generate Rayleigh scattering properties.
The Wigner d-function coefficients for either the Rayleigh scalar phase function
or polarized phase matrix for molecular scattering by air. A molecular Rayleigh
extinction profile is be computed from an input temperature profile (in Kelvin).

Computations are done using SHDOM fortran routines. From SHDOM docstring:
Rayleigh scattering including the wavelength depolarization factor.
From Mishchenko's book "Multiple scattering of light by particles:
Radiative Transfer and Coherent Backscattering", Cambridge, 2006. Thanks to Adrian Doicu.
"""
from collections import OrderedDict
import numpy as np
import xarray as xr

import at3d.core
import at3d.checks

def to_grid(wavelengths, atmosphere, rte_grid):
    """
    Interpolate atmosphere to rte_grid and compute Rayleigh optical properties profile
    for each input wavelength.

    Rayleigh optical properties include the scattering phase matrix (for polarized light),
    the single scattering albedo, and the extinction profile. The assumption is of a
    linear lapse rate between levels to compute the pressure at each level.
    The Rayleigh extinction is proportional to air density.

    Parameters
    ----------
    wavelengths: list of floats or numpy array or scalar
        A list/numpy array or scalar with wavelengths in [micron].
    atmosphere: xr.Dataset
        Dataset containing temperature [Kelvin] and pressure [mbar] variables and z coordinate [km].
    rte_grid: xr.Dataset
        Dataset containing at least z coordinate [km] and data_vars 'delx' [km] (resolution in x)
        and 'dely' [km] (resolution in y direction). Rayleigh properties will be interpolated
        onto this grid.

    Returns
    -------
    output: OrderedDict
        Dictionary with wavelengths as keys and corresponding Rayleigh xr.Datasets as values.

    Notes
    -----
    single scattering albedo is 1.0.
    """
    at3d.checks.check_grid(rte_grid)
    wavelengths = np.atleast_1d(wavelengths)
    atmosphere_on_rte_grid = atmosphere.interp({'z': rte_grid.z})
    rayleigh_poly_tables = compute_table(wavelengths).rename('legcoef')
    rayleigh_extinction = compute_extinction(
        wavelengths,
        atmosphere_on_rte_grid.temperature,
        surface_pressure=atmosphere_on_rte_grid.pressure.data[0]
    )

    rayleigh_ssalb = xr.DataArray(
        name='ssalb',
        data=np.ones(rayleigh_extinction.extinction.shape, dtype=np.float32),
        dims=['wavelength', 'z']
    )

    rayleigh = xr.merge([rayleigh_extinction, rayleigh_ssalb]).broadcast_like(rte_grid)

    rayleigh['phase_weights'] = (['num_micro', 'x', 'y', 'z', 'wavelength'],
                                 np.ones((1,) +rayleigh.extinction.shape, dtype=np.float32)
                                )
    rayleigh['table_index'] = (['num_micro', 'x', 'y', 'z', 'wavelength'],
                               np.ones((1,) +rayleigh.extinction.shape, dtype=np.int32)
                               )

    rayleigh = rayleigh.set_coords('table_index')
    rayleigh_final = xr.merge(
        [rayleigh, rayleigh_poly_tables.expand_dims(dim='table_index', axis=-2)]
        )


    output = OrderedDict()
    for wavelength in wavelengths:
        temp = rayleigh_final.sel({'wavelength': wavelength})
        #add a 'wavelength_center' attribute as this is what solver.RTE needs.
        temp.attrs['wavelength_center'] = wavelength
        temp = temp.assign_attrs(atmosphere.attrs)
        temp['delx'] = rte_grid.delx
        temp['dely'] = rte_grid.dely
        output[wavelength] = temp

    return output

def compute_table(wavelengths):
    """
    Retrieve the Rayleigh Wigner d-function coefficients of the phase matrix at
    each specified wavelength.

    Parameters
    ----------
    wavelengths: list of floats or numpy array or scalar
        A list/numpy array or scalar with wavelengths in [micron].

    Returns
    -------
    table: xr.Dataset
        A dataset of Legendre coefficients for each wavelength specified.

    Notes
    -----
    Includes the wavelength depolarization factor in accordance to [1].
    This was implemented in SHDOM tanks to Adrian Doicu.

    References
    ----------
    .. [1] Mishchenko, Michael I., Larry D. Travis, and Andrew A. Lacis.
    Multiple scattering of light by particles: radiative transfer and coherent
    backscattering. Cambridge University Press, 2006.
    """
    wavelengths = np.atleast_1d(wavelengths)
    legcoefs, table_types = zip(*[at3d.core.rayleigh_phase_function(wvl) for wvl in wavelengths])
    data_arrays = [
        xr.DataArray(name='rayleigh_table',
                     data=legcoef[..., None],
                     dims=['stokes_index', 'legendre_index', 'wavelength'],
                     coords={'stokes_index': (['stokes_index'], ['P11','P22','P33','P44','P12','P34']),
                             'wavelength': np.array([wavelength])})
        for legcoef, wavelength in zip(legcoefs, wavelengths) 
    ]
    table = xr.concat(data_arrays, dim='wavelength')
    table.attrs = {'table_type': 'vector',
                   'units': 'wavelength [micron]'}

    return table

def compute_extinction(wavelengths, temperature_profile, surface_pressure=1013.0):
    """
    Compute the Rayleigh extinction profile (as a function of altitude).

    Parameters
    ----------
    wavelengths: list / numpy array or scalar
        A list/numpy array or scalar with wavelengths in [micron].
    temperature_profile: xarray.DataArray
        A DataArray containing the altitude grid points in [km] and temperatures in [Kelvin].
    surface_pressure: float
        Surface pressure in units [mbar] (default value is 1013.0 mbar).

    Returns
    -------
    rayleigh_profile: xarray.Dataset
        A DataArray containing the altitude grid points in [km], temperatures in [Kelvin]
        and molecular extinction [km^-1].

    Notes
    -----
    Use the parameterization of [1]_ Eq. 30 for tau_R at
    sea level and convert to Rayleigh density coefficient:
    k = 0.03370*(p_sfc/1013.25)*tau_R  for k in K/(mbar km).

    References
    ----------
    .. [1] Bodhaine, B. A., Wood, N. B., Dutton, E. G., & Slusser, J. R. (1999).
        On Rayleigh Optical Depth Calculations, Journal of Atmospheric and Oceanic Technology,
        16(11), 1854-1861.
        URL: https://journals.ametsoc.org/view/journals/atot/16/11/1520-0426_1999_016_1854_orodc_2_0_co_2.xml.
    """

    wavelengths = np.atleast_1d(wavelengths)
    raylcoefs = 0.03370 * (surface_pressure / 1013.25) * 0.0021520 * (
            1.0455996 - 341.29061 / wavelengths ** 2 - 0.90230850 * wavelengths ** 2) / \
                (1 + 0.0027059889 / wavelengths ** 2 - 85.968563 * wavelengths ** 2)
    ext_profile = [xr.DataArray(
        dims='z',
        name='rayleigh_extinction_profile',
        coords={'z': temperature_profile.z.values, 'wavelength': wavelength},
        data=at3d.core.rayleigh_extinct(nzt=temperature_profile.size,
                                   zlevels=temperature_profile.z,
                                   temp=temperature_profile.data,
                                   raysfcpres=surface_pressure,
                                   raylcoef=raylcoef)
    ) for raylcoef, wavelength in zip(raylcoefs, wavelengths)]

    # Create a Rayleigh profile xarray.Dataset
    rayleigh_profile = temperature_profile.copy().to_dataset(name='temperature')
    rayleigh_profile['extinction'] = xr.concat(ext_profile, dim='wavelength')
    rayleigh_profile.attrs['units'] = ['wavelength [micron]', 'temperature [K]', 'extinction [km^-1]']
    return rayleigh_profile
