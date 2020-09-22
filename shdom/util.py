"""
Utility functions
"""
import numpy as np
from collections import OrderedDict
import shdom
import netCDF4 as nc
import xarray as xr

def float_round(x):
    """Round a float or np.float32 to a 3 digits float"""
    if type(x) == np.float32:
        x = x.item()
    return round(x,3)

def int_round(x):
    """Round a float or np.float32 to a 3 digits integer by 1000x scaling"""
    return int(np.round(x*1000))

def find_nearest(array, value):
    """Find the nearest element index in an array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def set_pyshdom_path():
    """set path to pyshdom parent directory"""
    import os, shdom
    from pathlib import Path
    os.chdir(str(Path(shdom.__path__[0]).parent))

def planck_function(temperature, wavelength, c=2.99792458e8,h=6.62606876e-34,k=1.3806503e-23):
    """
    temperature
        units, Kelvin
    wavelength
        units, micrometers
    radiance
        units, Watts/m^2/micrometer/steradian (SHDOM units)
    """
    wavelength = wavelength*1e-6
    radiance = 2*h*c**2/ wavelength**5 * 1.0 / (np.exp((h*c) / (wavelength*k*temperature)) - 1.0) * 1e-6
    return radiance


def load_forward_model(file_name):
    """
    TODO
    """
    dataset = nc.Dataset(file_name)

    groups = dataset.groups
    sensors = groups['sensors'].groups
    solvers = groups['solvers'].groups
    sensor_dict = OrderedDict()
    solver_dict = OrderedDict()

    for key,sensor in sensors.items():
        sensor_list = []
        for i, image in sensor.groups.items():
            sensor_list.append(xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                'sensors/'+str(key)+'/'+str(i)])))
        sensor_dict[key] = {'sensor_list':sensor_list}

    for key, solver in solvers.items():
        numerical_params = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                        'solvers/'+str(key)+'/numerical_parameters']))
        num_stokes = numerical_params.num_stokes.data
        surface = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/surface']))
        source = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/source']))

        medium_list = []
        for i, med in solver['medium'].groups.items():
            medium_list.append(xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                        'solvers/'+str(key)+'/medium/'+str(i)])))

        if 'atmosphere' in solver.groups.keys():
            atmosphere = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                'solvers/'+str(key)+'/atmosphere']))
        else:
            atmosphere=None

        solver_dict[float(key)] = shdom.solver.RTE(numerical_params=numerical_params,
                                            medium=medium_list,
                                           source=source,
                                           surface=surface,
                                            num_stokes=num_stokes,
                                            name=None
                                           )
        rte_grid = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/grid']))

    return sensor_dict, solver_dict, rte_grid

def save_forward_model(file_name,sensors, solvers):
    """
    TODO
    """
    for i, (key,sensor) in enumerate(sensors.items()):
        for j, image in enumerate(sensor['sensor_list']):
            if (i==0) & (j==0):
                image.to_netcdf(file_name, 'w', group = 'sensors/'+str(key)+'/'+str(j), format='NETCDF4',engine='netcdf4')
            else:
                image.to_netcdf(file_name, 'a', group = 'sensors/'+str(key)+'/'+str(j), format='NETCDF4',engine='netcdf4')

    for i, (key, solver) in enumerate(solvers.items()):

        numerical_params = solver.numerical_params
        numerical_params['num_stokes'] = solver._nstokes

        solver.numerical_params.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'numerical_parameters', format='NETCDF4',engine='netcdf4')
        solver.surface.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'surface', format='NETCDF4',engine='netcdf4')
        solver.source.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'source', format='NETCDF4',engine='netcdf4')
        solver._grid.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'grid', format='NETCDF4',engine='netcdf4')
        if solver.atmosphere is not None:
            solver.atmosphere.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'atmosphere', format='NETCDF4',engine='netcdf4')
        for j, med in enumerate(solver.medium):
            med.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'medium/'+str(j), format='NETCDF4',engine='netcdf4')
