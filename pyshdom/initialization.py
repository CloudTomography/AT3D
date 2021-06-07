"""
This module holds routines for initialization and preprocessing
based on measurements.
"""
import xarray as xr
import numpy as np
import pyshdom.space_carve

def mean_ext_estimate(rte_grid, sensors, solar_mu, solar_azimuth,
                     chi=2/3, g=0.86, sun_distance_reflect=0.1,
                     sun_distance_transmit=0.1):
    """
    Estimate the extinction of a cloud using diffusion theory.

    Given a masked volume `space_carved_volume`, the geometric distance
    between each point and the sun is calculated. The value of the geometric distance
    from the sun through the cloud at the first intersection of a sensor ray
    with the cloud volume is used to classify whether sensor pixels are
    observing shadowed or directly illuminated portions of the cloud.

    The mean of all 'shadowed' and 'illuminated' pixels is used to derive an
    optical diameter using diffusion theory and the extrapolation length `chi`
    and an asymmetry factor. This optical diameter is converted to an extinction
    using the length scale of the maximum chord length through the cloud in the solar
    direction. This length scale is chosen because it collapses to the relevant case
    for several geometries.
    """
    space_carver = pyshdom.space_carve.SpaceCarver(rte_grid)
    if isinstance(sensors, xr.Dataset):
        sensor_list = [sensors]
    elif isinstance(sensors, type([])):
        sensor_list = sensors
    elif isinstance(sensors, pyshdom.containers.SensorsDict):
        sensor_list = []
        for instrument in sensors:
            sensor_list.extend(sensors[instrument]['sensor_list'])

    volume = space_carver.carve(sensor_list, agreement=(0.0, 1.0), linear_mode=False)
    sundistance = space_carver.shadow_mask(volume.mask, sensor_list, solar_mu, solar_azimuth)

    reflected = []
    transmitted = []
    for sensor in sensor_list:
        reflected.extend(sensor.I.data[np.where((sensor.sun_distance.data < sun_distance_reflect) &
                                                (sensor.cloud_mask.data == 1))])
        transmitted.extend(sensor.I.data[np.where((sensor.sun_distance.data >= sun_distance_transmit) &
                                                  (sensor.cloud_mask.data == 1))])

    sundistance_radius = sundistance.sun_distance.data[np.where(sundistance.sun_distance > 0.0)].max()

    tau_estimate = 2*chi*np.mean(reflected)/np.mean(transmitted)/(1.0-g)
    ext_estimate = tau_estimate/sundistance_radius

    extinction = np.zeros(volume.mask.shape)
    extinction[np.where(volume.mask == 1.0)] = ext_estimate
    extinction = xr.Dataset(
        data_vars={
            'extinction': (['x', 'y', 'z'], extinction)
        },
        coords={
            'x': rte_grid.x,
            'y': rte_grid.y,
            'z': rte_grid.z,
        },
        attrs={
            'tau_estimate': tau_estimate,
            'chi': chi,
            'g': g,
            'radius': sundistance_radius
        }
    )
    return extinction
