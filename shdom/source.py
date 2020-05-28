import xarray as xr
import numpy as np

def solar_source(solar_zenith, solar_azimuth,solarflux=1.0,skyrad=0.0):
    """
    TODO
    """
    return xr.Dataset(
        data_vars={
            'name': 'solar_source'
            'solarflux': solarflux,
            'solar_zenith': solar_zenith,
            'solar_azimuth': solar_azimuth,
            'solarmu': np.cos(np.deg2rad(solar_zenith)),
            'solaraz': np.deg2rad(solar_azimuth),
            'srctype': 'S'
            'units': 'R', #only used for thermal
            'wavenumber': [10000,10001], #only used for thermal.
            'skyrad': skyrad #TODO figure out if skyrad is thermal only.
        }
    )
