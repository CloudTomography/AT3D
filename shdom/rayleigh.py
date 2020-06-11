"""
Rayleigh scattering for temperature profile.

Description taken from cloudprp.f:
     Computes the molecular Rayleigh extinction profile EXTRAYL [/km]
     from the temperature profile TEMP [K] at ZLEVELS [km].  Assumes
     a linear lapse rate between levels to compute the pressure at
     each level.  The Rayleigh extinction is proportional to air
     density, with the coefficient RAYLCOEF in [K/(mb km)].

Parameters
----------
wavelengths: float or list / numpy array.
    The wavelengths in [microns].
"""
from shdom import core
import numpy as np
import xarray as xr
from collections import OrderedDict

def to_grid(wavelengths, atmosphere, rte_grid):
    """
    TODO
    """
    atmosphere_on_rte_grid = atmosphere.interp({'z':rte_grid.z})
    rayleigh_poly_tables = compute_table(wavelengths).rename('legcoef')
    rayleigh_extinction = compute_extinction(wavelengths,atmosphere_on_rte_grid.Temperature,
                                                           surface_pressure=atmosphere_on_rte_grid.Pressure.data[0])

    rayleigh_ssalb = xr.DataArray(
                        name='ssalb',
                        data = np.ones(rayleigh_extinction.extinction.shape,dtype=np.float32),
                        dims=['wavelength','z']
    )

    rayleigh_table_index = xr.DataArray(
                        name='table_index',
                        data = np.ones(rayleigh_extinction.extinction.shape, dtype=np.int),
                        dims=['wavelength','z']
    )

    rayleigh = xr.merge([rayleigh_extinction,rayleigh_ssalb, rayleigh_table_index]).broadcast_like(rte_grid)

    rayleigh = rayleigh.set_coords('table_index')
    rayleigh_final = xr.merge([rayleigh, rayleigh_poly_tables.expand_dims(dim='table_index', axis=-2)])

    #TODO need to decide whether this should return a large Dataset with a wavelength dimension
    #or a list/dict. For consistency with the way that microphysical scatterers are determined at the script level,
    #the output is currently an OrderedDict. compute_table/rayleigh.compute_extinction COULD be reformaulted to be
    #mononcrhomatic for consistency with microphysics but that is not really necessary.

    output = OrderedDict()
    for wavelength in wavelengths:
        temp = rayleigh_final.sel({'wavelength': wavelength})
        #add a 'wavelength_center' attribute as this is what solver.RTE needs.
        temp.attrs['wavelength_center'] = wavelength
        output[wavelength] = temp

    return output

def compute_table(wavelengths):
    """
    Retrieve the Rayleigh phase function (Legendre table).

    Parameters
    ----------
    wavelengths: list or scalar
        A list/numpy array or scalar with wavelengths in [micron]

    Returns
    -------
    table: xarray.Dataset
        A dataset of Legendre coefficients for each wavelength specified.
    """
    wavelengths = np.atleast_1d(wavelengths)
    legcoefs, table_types = zip(*[core.rayleigh_phase_function(wavelen=wavelen) for wavelen in wavelengths])
    data_arrays = [
        xr.DataArray(name='rayleigh_table',
                     data=legcoef[..., None],
                     dims=['stokes_index', 'legendre_index', 'wavelength'],
                     coords={'stokes_index': (['stokes_index'], ['P11','P22','P33','P44','P12','P34']),
                             'wavelength': np.array([wavelength])})
        for legcoef, table_type, wavelength in zip(legcoefs, table_types, wavelengths)
    ]
    table = xr.concat(data_arrays, dim='wavelength')
    table.attrs = {'table_type': 'vector',
                   'units': 'wavelength [micron]'}

    return table

def compute_extinction(wavelengths, temperature_profile, surface_pressure=1013.0):
    """
    Retrieve the Rayleigh extinction profile (as a function of altitude).

    Parameters
    ----------
    wavelengths: list or scalar
        A list/numpy array or scalar with wavelengths in [micron]
    temperature_profile: xarray.DataArray
        a DataArray containing the altitude grid points in [km] and temperatures in [kelvin].
    surface_pressure: float
        Surface pressure in units [mb] (default value is 1013.0)

    Returns
    -------
    rayleigh_profile: xarray.Dataset
        a DataArray containing the altitude grid points in [km], temperatures in [kelvin] and molecular extinction [km^-1].
    """

    # Use the parameterization of Bodhaine et al. (1999) eq 30 for tau_R at
    # sea level and convert to Rayleigh density coefficient:
    #   k = 0.03370*(p_sfc/1013.25)*tau_R  for k in K/(mb km)
    wavelengths = np.atleast_1d(wavelengths)
    raylcoefs = 0.03370 * (surface_pressure / 1013.25) * 0.0021520 * (
            1.0455996 - 341.29061 / wavelengths ** 2 - 0.90230850 * wavelengths ** 2) / \
                (1 + 0.0027059889 / wavelengths ** 2 - 85.968563 * wavelengths ** 2)
    ext_profile = [xr.DataArray(
        dims='z',
        name='rayleigh_extinction_profile',
        coords={'z': temperature_profile.z.values, 'wavelength': wavelength},
        data=core.rayleigh_extinct(nzt=temperature_profile.size,
                                   zlevels=temperature_profile.z,
                                   temp=temperature_profile.data,
                                   raysfcpres=surface_pressure,
                                   raylcoef=raylcoef)
    ) for raylcoef, wavelength in zip(raylcoefs, wavelengths)]

    # Create a Rayleigh profile xarray.Dataset
    rayleigh_profile = temperature_profile.copy().to_dataset(name='temperature')#.swap_dims({'Altitude': 'z'}).drop('Altitude')
    rayleigh_profile['extinction'] = xr.concat(ext_profile, dim='wavelength')
    rayleigh_profile.attrs['units'] = ['wavelength [micron]', 'temperature [k]', 'extinction [km^-1]']
    return rayleigh_profile
