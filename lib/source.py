"""
TODO module docstring
"""

import xarray as xr
import numpy as np

def solar(wavelength, solarmu, solar_azimuth, solarflux=1.0, skyrad=0.0):
    """
    TODO
    """
    if not (-1.0 <= solarmu) & (solarmu < 0.0):
        raise ValueError("solarmu must be in the range -1.0 <= solarmu < 0.0 not '{}'. "
                         "The SHDOM convention for solar direction is that it points"
                         "in the direction of the propagation of radiance.".format(solarmu))

    return xr.Dataset(
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

#useful comment: proportionally scale the split_accuracy by the
#planck function at a typical temperature (as Frank Evans did)
#set spherical_harmonics_accuracy = 0.0, unless you have tested.
#when using both combined source: use the solar spectrum to calculate
#solar source at wavelength. scale split_accuracy by the sum.
#output units of radiance.
#thermal or combined source is not yet supported for gradient calculations.

def thermal(wavelength, skyrad=0.0, units='radiance'):
    """
    TODO
    """
    if units == 'radiance':
        units_flag = 'R'
    elif units == 'brightness_temperature':
        units_flag = 'T'

    return xr.Dataset(
        data_vars={
            'name': 'solar_source',
            'wavelength': wavelength,
            'solarflux': 0.0,
            'solarmu': np.nan,
            'solaraz': np.nan,
            'srctype': 'T',
            'units': units_flag,
            'wavenumber': [10000, 10001], #only used for CKD
            'skyrad': skyrad #in thermal only this is brightness temperature of the isotropic diffuse radiance.
        }
    )

def combined(wavelength, solarmu, solar_azimuth, solarflux=1.0, skyrad=0.0):
    """
    TODO
    """
    if not (-1.0 <= solarmu) & (solarmu < 0.0):
        raise ValueError("solarmu must be in the range -1.0 <= solarmu < 0.0 not '{}'. "
                         "The SHDOM convention for solar direction is that it points"
                         "in the direction of the propagation of radiance.".format(solarmu))

    return xr.Dataset(
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
