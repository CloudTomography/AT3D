"""
This module contains functions for generating the parameters that describe
the surface BRDF, the bottom boundary condition used in solver.RTE (SHDOM).
Rather than reading the surface parameters from a file as in SHDOM, they are
directly generated and are formatted by a fortran subroutine at3d.core.prep_surface
which is a minor modification of shdomsub3.f READ_SURFACE.

Each of these surfaces have corresponding fortran subroutines which evaluate
the BRDF based on the prescribed parameters.
Implementing a new surface would require writing a new subroutine that evaluates
the BRDF to shdomsub2.f (model after RPV_REFLECTION). Then add a block to
SURFACE_BRDF in shdomsub2.f that calls this model, and at3d.core.prep_surface
(PREP_SURFACE in surface.f) that prepares the input parameters. Additionally,
a python wrapper should be added here, and the chosen `sfctype` should be added
to the list of valid surface types in solver.RTE._setup_surface.

The spatially varying surface parameters are defined using a linear interpolation
kernel on a uniformly spaced grid in each of the 'x' and 'y' directions. This surface
grid does not have to correspond with the grid used to define the volumetric
optical properties in solver.RTE. Linear interpolation will be used to map from
the surface grid to solver.RTE grid within SHDOM.
"""
import xarray as xr
import numpy as np

import at3d.core

def lambertian(albedo, ground_temperature=298.15, delx=None, dely=None):
    """
    Defines either a fixed or spatially variable Lambertian surface for use
    in solver.RTE (SHDOM).

    A Lambertian surface has an isotropic reflected radiance. so that the
    upwelling (reflected) radiance is given by the ratio downwelling_flux / pi.
    Naturally, the emissivity is also isotropic and is given by 1 - `albedo`.

    Parameters
    ----------
    albedo : float or np.ndarray
        values should be in range [0, 1]. Can be of shape=(nxsfc, nysfc)
        where nxsfc and nysfc are the number of points in the x and y directions.
    ground_temperature : float or np.ndarray
        The blackbody temperature of the surface emission with units of [K].
        Can be of shape=(nxsfc, nysfc)  where nxsfc and nysfc are the number of
        points in the x and y directions.
    delx : float or None, optional
        The spacing of the surface parameters in their 1st dimension, in units of [km].
        If size of `albedo`/`ground_temperature` is not equal to 1, then this must be
        specified.
    dely : float or None, optional.
        The spacing of the surface parameters in their 2nd dimension, in units of [km].
        If size of `albedo`/`ground_temperature` is not equal to 1, then this must be
        specified.

    Returns
    -------
    dataset : xr.Dataset
        The dataset contains the parameters for the surface along with the correct
        flags for use by solver.RTE (SHDOM).

    Raises
    ------
    ValueError
        If the `albedo` is not in [0, 1].
        If `albedo` and `ground_temperature` don't have matching shapes.
    """
    if np.any(albedo > 1.0) or np.any(albedo < 0.0):
        raise ValueError("surface albedo should be in [0, 1] not '{}'".format(albedo))

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
    Defines either a fixed or spatially varying polarized surface BRDF for
    rough (oceanic) surfaces.

    Models the surface as a polarized Fresnel reflection from a dielectric
    interface with a Gaussian distribution of slopes, including shadowing.
    The properties of the interface can be fixed or spatially varying. The
    surface roughness (mean square slope) is calculated from the parameterization
    of Cox and Munk.

    Parameters
    ----------
    real_refractive_index : float or np.ndarray
        The real refractive index of the surface medium. Can be of shape=(nxsfc, nysfc)
        where nxsfc and nysfc are the number of points in the x and y directions.
    imaginary_refractive_index : float or np.ndarray
        The imaginary refractive index of the surface medium. This should be
        negative for absorbing surfaces. Can be of shape=(nxsfc, nysfc)
        where nxsfc and nysfc are the number of points in the x and y directions.
    surface_wind_speed : float or np.ndarray
        The surface wind speed in m/s. Can be of shape=(nxsfc, nysfc)
        where nxsfc and nysfc are the number of points in the x and y directions.
    ground_temperature : float or np.ndarray
        The blackbody temperature of the surface emission units of [K].
        Can be of shape=(nxsfc, nysfc) where nxsfc and nysfc are the number
        of points in the x and y directions.
    delx : float or None, optional
        The spacing of the surface parameters in their 1st dimension, in units of [km].
        If size of `ground_temperature` or other parameters is not equal to 1,
        then this must be specified.
    dely : float or None, optional.
        The spacing of the surface parameters in their 2nd dimension, in units of [km].
        If size of `ground_temperature` or other parameters is not equal to 1,
        then this must be specified.

    Returns
    -------
    dataset : xr.Dataset
        The dataset contains the parameters for the surface along with the correct
        flags for use by solver.RTE (SHDOM).

    Raises
    ------
    ValueError
        If the `albedo` is not in [0, 1].
        If all surface parameters (arguments) don't have matching shapes.
        If size of surface parameters is larger than 1 but delx, dely are
        not specified.

    Notes
    -----
    See WAVE_FRESNEL_REFLECTION in shdomsub2.f for the implementation of the
    BRDF calculation. When the windspeed becomes small high angular resolution (NMU, NPHI)
    is needed to resolve the specular reflection peak. For this reason,
    tests should be always be performed to ensure convergence.
    No checks are performed on the input values, users should read and understand the
    reference.

    Reference
    ---------
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


    dataset = _make_surface_dataset('wave_fresnel', ground_temperature, delx, dely,
                                 real_refractive_index=real_refractive_index,
                                 imaginary_refractive_index=imaginary_refractive_index,
                                 surface_wind_speed=surface_wind_speed)
    return dataset

def diner(A, K, B, ZETA, SIGMA, ground_temperature=298.15, delx=None, dely=None):
    """
    Defines either a fixed or spatially varying polarized surface BRDF for
    generic rough surfaces including a volumetric scattering term.

    Computes polarized reflection from the Diner et al. (2012) model
    which two terms: 1) completely depolarizing volumetric scattering with
    the modified RPV model, and 2) Fresnel reflection from randomly
    oriented microfacets with either a uniform or Gaussian distribution
    of slopes. This BRDF model also includes the hotspot factor which
    was not included in Diner et al. (2012).

    Parameters
    ----------
    ground_temperature : float or np.ndarray
        The blackbody temperature of the surface emission with units of [K].
        Can be of shape=(nxsfc, nysfc) where nxsfc and nysfc are the number of
        points in the x and y directions.
    SIGMA : float or np.ndarray
        The standard deviation of the gaussian distribution of microfacets of
        the surface (roughness). If SIGMA <= 0.0 the uniform orientation distribution
        is used. Can be of shape=(nxsfc, nysfc) where nxsfc and nysfc are the number of
        points in the x and y directions.
    ZETA : float or np.ndarray
        A scaling parameter which sets the magnitude of the microfacet term of the
        BRDF.
    A : float or np.ndarray
        A scaling parameter which sets the magnitude of the depolarizing volumetric
        scattering term of the BRDF.
    K : float or np.ndarray
        Describes the strength of the azimuthally independent component of the
        depolarizing volumetric scattering (modified RPV model) term in the BRDF.
        Isotropic volume component has K=1.0.
    B : float or np.ndarray
        controls the scattering angle dependence, (forward/backward asymmetry) and
        hence all azimuthal dependence) of the depolarizing volumetric scattering
        (modified RPV model) term in the BRDF. Isotropic volume component has B=0.0.
    delx : float or None, optional
        The spacing of the surface parameters in their 1st dimension, in units of [km].
        If the size of `ground_temperature` or other parameters is not equal to 1, then
        this must be specified.
    dely : float or None, optional.
        The spacing of the surface parameters in their 2nd dimension, in units of [km].
        If the size of `ground_temperature` or other parameters is not equal to 1, then
        this must be specified.

    Returns
    -------
    dataset : xr.Dataset
        The dataset contains the parameters for the surface along with the correct
        flags for use by solver.RTE (SHDOM).

    Raises
    ------
    ValueError
        If all surface parameters (arguments) don't have matching shapes.
        If size of surface parameters is larger than 1 but delx, dely are
        not specified.

    Notes
    -----
    See the fortran subroutine DINER_REFLECTION (shdomsub2.f) for the implementation.
    The complex refractive index is hardcoded to  1.5 - 0.0i.
    No checks are performed on the input values, users should read and understand the
    reference.

    Reference
    ---------
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

    dataset = _make_surface_dataset('diner', ground_temperature, delx, dely,
                                 A=A, K=K, B=B, ZETA=ZETA, SIGMA=SIGMA)
    return dataset


def ocean_unpolarized(surface_wind_speed, pigmentation, ground_temperature=298.15,
                      delx=None, dely=None):
    """
    Defines either a fixed or spatially varying scalar surface BRDF for
    rough (oceanic) surfaces.

    This ocean model is extracted and modified from the 6S model (reference below).
    The backscattered reflectance from the seawater is modeled as as
    lambertian, parameterized by the `pigmentation` and salinity.
    Sunglint is modeled taking into account surface roughness based on
    a Gaussian distribution of slopes parameterized according to Cox and Munk.
    The effect of whitecaps is also modeled as Lambertian.

    Parameters
    ----------
    surface_wind_speed : float or np.ndarray
        Minimum value of 0.25 m/s
        The surface wind speed in m/s. Can be of shape=(nxsfc, nysfc)
        where nxsfc and nysfc are the number of points in the x and y directions.
    pigmentation : float or np.ndarray
        The pigment concentration [mg/m^3] used to set the index of refraction.
        Can be of shape=(nxsfc, nysfc) where nxsfc and nysfc are the number of
        points in the x and y directions.
    ground_temperature : float or np.ndarray
        The blackbody temperature of the surface emission units of [K]. Can be
        of shape=(nxsfc, nysfc) where nxsfc and nysfc are the number of points
        in the x and y directions.
    delx : float or None, optional
        The spacing of the surface parameters in their 1st dimension, in units of [km].
        If size of `ground_temperature` or other parameters is not equal to 1,
        then this must be specified.
    dely : float or None, optional.
        The spacing of the surface parameters in their 2nd dimension, in units of [km].
        If size of `ground_temperature` or other parameters is not equal to 1,
        then this must be specified.

    Returns
    -------
    dataset : xr.Dataset
        The dataset contains the parameters for the surface along with the correct
        flags for use by solver.RTE (SHDOM).

    Raises
    ------
    ValueError
        If all surface parameters (arguments) don't have matching shapes.
        If size of surface parameters is larger than 1 but delx, dely are
        not specified.

    Notes
    -----
    See ocean_brdf.f for all subroutines.
    Salinity is hardcoded to a default value of 34.3 ppt in
    SURFACE_BRDF (shdomsub2.f). The azimuthal orientation of the wind
    is always along 0.0.
    When the windspeed becomes small high angular resolution (NMU, NPHI)
    is needed to resolve the specular reflection peak. For this reason,
    tests should be always be performed to ensure convergence.
    No checks are performed on the input values, users should read and understand the
    reference.

    Reference
    ---------
    E. F. Vermote, D. Tanre, J. L. Deuze, M. Herman and J. -. Morcette,
    "Second Simulation of the Satellite Signal in the Solar Spectrum, 6S:
    an overview," in IEEE Transactions on Geoscience and Remote Sensing,
    vol. 35, no. 3, pp. 675-686, May 1997, doi: 10.1109/36.581987.
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

    dataset = _make_surface_dataset('ocean_unpolarized', ground_temperature, delx, dely,
                                 surface_wind_speed=surface_wind_speed,
                                 pigmentation=pigmentation)
    return dataset


def RPV_unpolarized(RHO0, K, THETA, ground_temperature=298.15, delx=None, dely=None):
    """
    Defines either a fixed or spatially varying scalar surface BRDF for
    generic natural surfaces.

    The empirical RPV model is used appropriate for scalar radiative transfer.
    This BRDF model also includes a hotspot factor.

    Parameters
    ----------
    ground_temperature : float or np.ndarray
        The blackbody temperature of the surface emission with units of [K].
        Can be of shape=(nxsfc, nysfc) where nxsfc and nysfc are the number of
        points in the x and y directions.
    RHO0 : float or np.ndarray
        A scaling parameter which sets the magnitude of the BRDF function.
    K : float or np.ndarray
        Describes the strength of the azimuthally independent component of the
        term in the BRDF. Larger K increases anisotropy, K=1.0 is required for
        isotropy.
    THETA : float or np.ndarray
        The scattering angle dependence is modeled as a henyey greenstein phase
        function where THETA is the asymmetry parameter. THETA should be in
        range [-1, 1], where more positive indicates more forward scattering.
        THETA=0.0 is required for isotropy.
    delx : float or None, optional
        The spacing of the surface parameters in their 1st dimension, in units of [km].
        If the size of `ground_temperature` or other parameters is not equal to 1, then
        this must be specified.
    dely : float or None, optional.
        The spacing of the surface parameters in their 2nd dimension, in units of [km].
        If the size of `ground_temperature` or other parameters is not equal to 1, then
        this must be specified.

    Returns
    -------
    dataset : xr.Dataset
        The dataset contains the parameters for the surface along with the correct
        flags for use by solver.RTE (SHDOM).

    Raises
    ------
    ValueError
        If all surface parameters (arguments) don't have matching shapes.
        If size of surface parameters is larger than 1 but delx, dely are
        not specified.

    Notes
    -----
    See the fortran subroutine RPV_REFLECTION (shdomsub2.f) for the implementation.
    No checks are performed on the input values, users should read and understand the
    reference.

    Reference
    ---------
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

    dataset = _make_surface_dataset('rpv_unpolarized', ground_temperature, delx, dely,
                                 RHO0=RHO0, K=K, THETA=THETA)
    return dataset

def _make_surface_dataset(surface_type, ground_temperature, delx, dely, **kwargs):
    """
    This function prepares surface inputs generated by the python functions in
    this module for use in solver.RTE (SHDOM).

    This method reorganizes the inputs into the correct format
    and sets the array size parameters. The filling of the horizontal boundary
    conditions and other preparations is done by at3d.core.prep_surface
    (PREP_SURFACE in surface.f). Only `surface_type` implemented in SHDOM
    are allowed. This function should not be used directly, instead, use the
    individual functions which are defined in this module.

    Parameters
    ----------
    surface_type : string
        The type of the surface to be generated.
    ground_temperature : float or np.ndarray
        The array of surface blackbody emission temperatures.
    delx : float
        The 'x' spacing of gridpoints in variable surfaces in units of [km].
    dely : float
        The 'y' spacing of gridpoints in variable surfaces in units of [km].
    kwargs : np.ndarray
        Each entry in kwargs should be a parameter describing the surface
        (other than `ground_temperature`). Note that the order of these is
        significant. See the use of this function within surface.py for
        examples.

    Returns
    -------
    surface_dataset : xr.Dataset
        The dataset contains the parameters for the surface along with the correct
        flags for use by solver.RTE (SHDOM).

    Notes
    -----
    There are no checks performed to make sure that the correct kwargs for each
    surface type are passed, or that they are in the correct order. For this reason,
    do not use this function in isolation.
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
    nsfcpar, sfcparms, gndtemp, gndalbedo, ierr, errmsg = at3d.core.prep_surface(maxsfcpts=maxsfcpts,
                                                                      maxsfcpars=maxsfcpars,
                                                                      sfctype=sfctype,
                                                                      nxsfc=nxsfc,
                                                                      nysfc=nysfc,
                                                                      delxsfc=delxsfc,
                                                                      delysfc=delysfc,
                                                                      parms_in=parms_in,
                                                                      grid_coords=grid_coords)
    at3d.checks.check_errcode(ierr, errmsg)
    surface_dataset = xr.Dataset(
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
    return surface_dataset
