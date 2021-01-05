"""
TODO module docstring
"""
import xarray as xr
import numpy as np

import pyshdom.core

def lambertian(albedo, ground_temperature=298.15, delx=None, dely=None):
    """
    TODO
    """
    ground_temperature = np.atleast_2d(ground_temperature)
    albedo = np.atleast_2d(albedo)

    if (albedo.size == 1) & (ground_temperature.size == 1):
        dataset = xr.Dataset(
            data_vars={
                'name': 'fixed_lambertian_surface',
                'sfctype':'FL',
                'gndalbedo':albedo[0, 0],
                'gndtemp':ground_temperature[0, 0],
                'maxsfcpars':4,
                'nxsfc': 0,
                'nysfc': 0,
                'delxsfc': 0,
                'delysfc': 0,
                'nsfcpar': 1,
                'sfcparms': [],
            }
            )
    elif ground_temperature.size == 1:
        ground_temperature = np.full_like(albedo, fill_value=ground_temperature[0, 0])
        if dely is None or delx is None:
            raise ValueError('dely and delx must be defined for variable surface parameters.')
        dataset = _make_surface_dataset('variable_lambertian', ground_temperature,
                                        delx, dely, albedo=albedo)
    else:
        raise ValueError('ground temperature must have a compatible shape.')
    return dataset

def wave_fresnel(real_refractive_index, imaginary_refractive_index, surface_wind_speed,
                 ground_temperature=298.15, delx=None, dely=None):
    """
    TODO
    The mean square surface slope (s^2) or SIGMA2 in Mishchenko and Travis
    is obtained from the Cox and Munk formula in the paper.
       For all formulas and definitions, see the paper:
    M. I. Mishchenko and L. D. Travis, 1997: Satellite retrieval
    of aerosol properties over the ocean using polarization as well as
    intensity of reflected sunlight.  J. Geophys. Res. 102, 16989-17013.
    """
    ground_temperature = np.atleast_2d(ground_temperature)
    real_refractive_index = np.atleast_2d(real_refractive_index)
    imaginary_refractive_index = np.atleast_2d(imaginary_refractive_index)
    surface_wind_speed = np.atleast_2d(surface_wind_speed)

    if not ((real_refractive_index.shape == imaginary_refractive_index.shape) &
            (real_refractive_index.shape == surface_wind_speed.shape)):
        raise ValueError('All surface brdf parameters must have the same shape.')

    if ground_temperature.size == 1:
        ground_temperature = np.full_like(real_refractive_index,
                                          fill_value=ground_temperature[0, 0])
    else:
        raise ValueError('ground temperature must have a compatible shape')

    if real_refractive_index.size == 1:
        if delx is None:
            delx = 0.02
        if dely is None:
            dely = 0.02
    else:
        if dely is None or delx is None:
            raise ValueError('dely and delx must be defined for variable surface parameters.')


    return _make_surface_dataset('wave_fresnel', ground_temperature, delx, dely,
                                 real_refractive_index=real_refractive_index,
                                 imaginary_refractive_index=imaginary_refractive_index,
                                 surface_wind_speed=surface_wind_speed)

def diner(A, K, B, ZETA, SIGMA, ground_temperature=298.15, delx=None, dely=None):
    """
    TODO
    For all formulas and definitions, see the paper:
        Diner, D. J., F. Xu, J. V. Martonchik, B. E. Rheingans, S. Geier,
        V. M. Jovanovic, A. Davis, R. A. Chipman, S. C. McClain, 2012:
        Exploration of a polarized surface bidirectional reflectance model
        using the ground-based multiangle spectropolarimetric imager.
        Atmosphere 2012, 3, 591-619; doi:10.3390/atmos3040591.
    """
    ground_temperature = np.atleast_2d(ground_temperature)
    A = np.atleast_2d(A)
    B = np.atleast_2d(B)
    K = np.atleast_2d(K)
    ZETA = np.atleast_2d(ZETA)
    SIGMA = np.atleast_2d(SIGMA)

    if not ((A.shape == K.shape) &
            (A.shape == B.shape) & (A.shape == ZETA.shape) & (A.shape == SIGMA.shape)):
        raise ValueError('All surface brdf parameters must have the same shape.')

    if ground_temperature.size == 1:
        ground_temperature = np.full_like(A, fill_value=ground_temperature[0, 0])
    else:
        raise ValueError('ground temperature must have a compatible shape')

    if A.size == 1:
        if delx is None:
            delx = 0.02
        if dely is None:
            dely = 0.02
    else:
        if dely is None or delx is None:
            raise ValueError('dely and delx must be defined for variable surface parameters.')

    return _make_surface_dataset('diner', ground_temperature, delx, dely,
                                 A=A, K=K, B=B, ZETA=ZETA, SIGMA=SIGMA)


def ocean_unpolarized(surface_wind_speed, pigmentation, ground_temperature=298.15,
                      delx=None, dely=None):
    """
    TODO
    """
    ground_temperature = np.atleast_2d(ground_temperature)
    pigmentation = np.atleast_2d(pigmentation)
    surface_wind_speed = np.atleast_2d(surface_wind_speed)

    if not surface_wind_speed .shape == pigmentation.shape:
        raise ValueError('All surface brdf parameters must have the same shape.')

    if ground_temperature.size == 1:
        ground_temperature = np.full_like(pigmentation, fill_value=ground_temperature[0, 0])
    else:
        raise ValueError('ground temperature must have a compatible shape')

    if pigmentation.size == 1:
        if delx is None:
            delx = 0.02
        if dely is None:
            dely = 0.02
    else:
        if dely is None or delx is None:
            raise ValueError('dely and delx must be defined for variable surface parameters.')

    return _make_surface_dataset('ocean_unpolarized', ground_temperature, delx, dely,
                                 surface_wind_speed=surface_wind_speed,
                                 pigmentation=pigmentation)


def RPV_unpolarized(RHO0, K, THETA, ground_temperature=298.15, delx=None, dely=None):
    """
    TODO
    The reference is:
       Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere
       Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable
       With NOAA Advanced Very High Resolution Radiometer Data,
       J. Geophys. Res., 98, 20791-20801.
    """

    ground_temperature = np.atleast_2d(ground_temperature)
    RHO0 = np.atleast_2d(RHO0)
    THETA = np.atleast_2d(THETA)
    K = np.atleast_2d(K)


    if not ((RHO0.shape == THETA.shape) &
            (RHO0.shape == K.shape)):
        raise ValueError('All surface brdf parameters must have the same shape.')

    if ground_temperature.size == 1:
        ground_temperature = np.full_like(K, fill_value=ground_temperature[0, 0])
    else:
        raise ValueError('ground temperature must have a compatible shape')

    if K.size == 1:
        if delx is None:
            delx = 0.02
        if dely is None:
            dely = 0.02
    else:
        if dely is None or delx is None:
            raise ValueError('dely and delx must be defined for variable surface parameters.')

    return _make_surface_dataset('rpv_unpolarized', ground_temperature, delx, dely,
                                 RHO0=RHO0, K=K, THETA=THETA)

def _make_surface_dataset(surface_type, ground_temperature, delx, dely, **kwargs):
    """
    TODO
    A hidden method which does the hard work for preparing surfaces.
    """
    nxsfc, nysfc = ground_temperature.shape
    maxsfcpts = (nxsfc+1)*(nysfc+1)

    if surface_type == 'variable_lambertian':
        maxsfcpars = 2
        sfctype = 'VL'
    elif surface_type == 'wave_fresnel':
        maxsfcpars = 4
        sfctype = 'VW'
    elif surface_type == 'diner':
        maxsfcpars = 6
        sfctype = 'VD'
    elif surface_type == 'ocean_unpolarized':
        maxsfcpars = 3
        sfctype = 'VO'
    elif surface_type == 'rpv_unpolarized':
        maxsfcpars = 4
        sfctype = 'VR'
    else:
        raise NotImplementedError

    delxsfc, delysfc = delx, dely
    grid_coords_temp = np.meshgrid(np.arange(1, nxsfc+1), np.arange(1, nysfc+1), indexing='ij')
    grid_coords = np.stack([grid_coords_temp[0].ravel(order='F'),
                            grid_coords_temp[1].ravel(order='F')], axis=0)
    list_of_params = [ground_temperature.ravel(order='F')] + [val.ravel(order='F')
                                                              for val in kwargs.values()]

    parms_in = np.stack(list_of_params, axis=0)
    nsfcpar, sfcparms, gndtemp, gndalbedo = pyshdom.core.prep_surface(maxsfcpts=maxsfcpts,
                                                                      maxsfcpars=maxsfcpars,
                                                                      sfctype=sfctype,
                                                                      nxsfc=nxsfc,
                                                                      nysfc=nysfc,
                                                                      delxsfc=delxsfc,
                                                                      delysfc=delysfc,
                                                                      parms_in=parms_in,
                                                                      grid_coords=grid_coords)
    return xr.Dataset(
        data_vars={
            'name': surface_type,
            'sfctype': sfctype,
            'gndalbedo': gndalbedo,
            'gndtemp': gndtemp,
            'maxsfcpars': maxsfcpars,
            'nxsfc': nxsfc,
            'nysfc': nysfc,
            'delxsfc': delxsfc,
            'delysfc': delysfc,
            'nsfcpar': nsfcpar,
            'sfcparms': ('nparameters', sfcparms),
        }
    )
