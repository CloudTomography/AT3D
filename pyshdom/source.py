"""This module contains methods for defining the three different types
of sources used in SHDOM; `solar`, `thermal`, and `combined`, which
produce collimated flux from above, thermal emission or a combination of the two.
The source also defines the top boundary condition, which may be a
downwelling isotropic radiation field.
Note that only monochromatic sources are currently supported,
so spectral integration schemes must be implemented 'manually' at the script level.
Gradient calculations are only available for `solar` sources, not
`thermal` or `combined`.
"""
import xarray as xr
import numpy as np

def solar(wavelength, solarmu, solar_azimuth, solarflux=1.0, skyrad=0.0):
    """Defines a solar source for use in solver.RTE (SHDOM).

    Parameters
    ----------
    wavelength : float
        The monochromatic (center) wavelength for use in the solar source in [microns]
    solarmu : float
        cosine of solar zenith angle. Points in the direction of propagation of
        the solar beam. Is required to be negative, if positive values are supplied
        they are set to negative.
    solar_azimuth : float
        The azimuthal angle of the solar beam. Points in the direction of
        propagation of the solar beam.
    solarflux : float
        The flux of the solar beam on a horizontal surface.
    skyrad : float
        The RADIANCE of the downwelling isotropic radiance field at domain top.

    Returns
    -------
    source_dataset : xr.Dataset
        The dataset contains the flags, units and above parameters that define
        the source for solver.RTE (SHDOM).

    Raises
    ------
    ValueError
        If `solarmu` is outside of the range [-1.0, 0.0). Note that `solarmu`=0.0
        (purely horizontal solar beam) is not permitted.

    Notes
    -----
    Several of the numerical parameters used in SHDOM such as the splitting accuracy
    and spherical harmonics accuracy are absolute measures and so they should scale
    with the magnitude of the radiation field and therefore the prescribed source.
    To avoid having to retune the numerical parameters, it is advisible to set
    `solarflux`=1.0 for all purely solar problems.
    """
    solarmu = -1*np.abs(solarmu)
    if not (-1.0 <= solarmu) & (solarmu < 0.0):
        raise ValueError("solarmu must be in the range -1.0 <= solarmu < 0.0 not '{}'. "
                         "The SHDOM convention for solar direction is that it points"
                         "in the direction of the propagation of radiance.".format(solarmu))

    source_dataset = xr.Dataset(
        data_vars={
            'name': 'solar_source',
            'wavelength': wavelength,
            'solarflux': solarflux,
            'solarmu': solarmu,
            'solaraz': np.deg2rad(solar_azimuth),
            'srctype': 'S',
            'units': 'R',
            'wavenumber': [10000, 10001], #only used for CKD
            'skyrad': skyrad #isotropic diffuse radiance from above
        }
    )
    return source_dataset

def thermal(wavelength, skyrad=0.0, units='radiance'):
    """Defines a thermal emission source for solver.RTE (SHDOM).

    Parameters
    ----------
    wavelength : float
        The monochromatic (center) wavelength for use in the solar source in [microns]
    skyrad : float
        The BRIGHTNESS TEMPERATURE of the downwelling isotropic radiance field at domain top.
    units : string
        Determines whether to set units to spectral radiance [W/m^2/steradian/micron] or
        to units of brightness temperature [K].

    Returns
    -------
    source_dataset : xr.Dataset
        The dataset contains the flags, units and above parameters that define
        the source for solver.RTE (SHDOM).
    Raises
    ------
    ValueError
        If the prescribed units are not either 'radiance' or 'brightness_temperature'

    Notes
    -----
    skyrad is now in units of BRIGHTNESS TEMPERATURE.
    When using a thermal source, it is not normalizable like the solar source.
    This means that care must be taken with the prescription of the accompanying
    numerical parameters such as the splitting accuracy.
    A good rule of thumb is to proportionally scale splitting accuracy and
    spherical harmonics accuracy based on the largest flux in the domain,
    which can be approximated as the flux emitted by the surface.
    """
    if units == 'radiance':
        units_flag = 'R'
    elif units == 'brightness_temperature':
        units_flag = 'T'
    else:
        raise ValueError("`units` should be either 'radiance' or 'brightness_temperature'.")

    source_dataset = xr.Dataset(
        data_vars={
            'name': 'solar_source',
            'wavelength': wavelength,
            'solarflux': 0.0,
            'solarmu': -0.5,
            'solaraz': 0.0,
            'srctype': 'T',
            'units': units_flag,
            'wavenumber': [10000, 10001], #only used for CKD
            'skyrad': skyrad #in thermal only this is brightness temperature of the isotropic diffuse radiance.
        }
    )
    return source_dataset

def combined(wavelength, solarmu, solar_azimuth, solarflux=1.0, skyrad=0.0):
    """Defines a combined (both solar and thermal) source for use in solver.RTE (SHDOM).

    Parameters
    ----------
    wavelength : float
        The monochromatic (center) wavelength for use in the solar source in [microns]
    solarmu : float
        cosine of solar zenith angle. Points in the direction of propagation of
        the solar beam. Is required to be negative, if positive values are supplied
        they are set to negative.
    solar_azimuth : float
        The azimuthal angle of the solar beam. Points in the direction of
        propagation of the solar beam.
    solarflux : float
        The flux of the solar beam on a horizontal surface.
    skyrad : float
        The RADIANCE of the downwelling isotropic radiance field at domain top.

    Returns
    -------
    source_dataset : xr.Dataset
        The dataset contains the flags, units and above parameters that define
        the source for solver.RTE (SHDOM).

    Raises
    ------
    ValueError
        If `solarmu` is outside of the range [-1.0, 0.0). Note that `solarmu`=0.0
        (purely horizontal solar beam) is not permitted.

    Notes
    -----
    Units of output (including skyrad) are in radiance [W/m^2/steradian/micron].
    When using a combined source, the `solarflux` should be set to its non-normalized
    value so that its relative magnitude with the thermal emission is correct.
    In this case, the split accuracy and spherical harmonics accuracy should again be
    correspondingly scaled, these could be chosen as a proportionally scaling by
    the sum of `solarflux` and the flux emitted from the surface, for example.
    """
    solarmu = -1*np.abs(solarmu)
    if not (-1.0 <= solarmu) & (solarmu < 0.0):
        raise ValueError("solarmu must be in the range -1.0 <= solarmu < 0.0 not '{}'. "
                         "The SHDOM convention for solar direction is that it points"
                         "in the direction of the propagation of radiance.".format(solarmu))

    source_dataset = xr.Dataset(
        data_vars={
            'name': 'solar_source',
            'wavelength': wavelength,
            'solarflux': solarflux,
            'solarmu': solarmu,
            'solaraz': np.deg2rad(solar_azimuth),
            'srctype': 'B',
            'units': 'R',
            'wavenumber': [10000, 10001], #only used for CKD
            'skyrad': skyrad #In combined this is the same as solar (radiance)
        }
    )
    return source_dataset
