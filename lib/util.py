"""
Utility functions
"""
import numpy as np
from collections import OrderedDict
import netCDF4 as nc
import xarray as xr
import pandas as pd

import pyshdom.core
import pyshdom.solver
import pyshdom.grid

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
    import os
    from pathlib import Path
    os.chdir(str(Path(pyshdom.__path__[0]).parent))

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

def cell_average_comparison(reference, other, variable_name):
    """
    calculates average values of 'variable name' in the cells
    of reference's grid for both reference and other (other is on a different grid.)
    """
    xgrid1 = reference.x.data
    ygrid1 = reference.y.data
    zgrid1 = reference.z.data
    values1 = reference[variable_name].data

    xgrid2 = other.x.data
    ygrid2 = other.y.data
    zgrid2 = other.z.data
    values2 = other[variable_name].data

    ref_vol,ref_val, other_vol, other_val = pyshdom.core.cell_average(xgrid1=xgrid1,ygrid1=ygrid1,zgrid1=zgrid1,
                        xgrid2=xgrid2,ygrid2=ygrid2,zgrid2=zgrid2,
                        values1=values1,
                        values2=values2)
    cell_average_ref = np.zeros(ref_vol.shape)
    cell_average_ref[np.where(ref_vol > 0.0)] = ref_val[np.where(ref_vol > 0.0)] /ref_vol[np.where(ref_vol > 0.0)]
    cell_average_other= np.zeros(other_vol.shape)
    cell_average_other[np.where(other_vol > 0.0)] = other_val[np.where(other_vol > 0.0)] /other_vol[np.where(other_vol > 0.0)]
    return cell_average_ref, cell_average_other

def load_2parameter_lwc_file(file_name, density='lwc'):
    """
    TODO
    Function that loads a scatterer from the '2 parameter lwc file' format used by
    SHDOM and i3rc monte carlo model.
    """
    header = pd.read_csv(file_name, nrows=4)
    nx, ny, nz = np.fromstring(header['2 parameter LWC file'][0], sep=' ').astype(np.int)
    dx, dy = np.fromstring(header['2 parameter LWC file'][1], sep=' ').astype(np.float)
    z = np.fromstring(header['2 parameter LWC file'][2], sep=' ').astype(np.float)
    temperature = np.fromstring(header['2 parameter LWC file'][3], sep=' ').astype(np.float)
    dset = pyshdom.grid.make_grid(dx, nx, dy, ny, z)

    data = np.genfromtxt(file_name, skip_header=5)

    lwc = np.zeros((nx, ny, nz))*np.nan
    reff = np.zeros((nx, ny, nz))*np.nan

    i, j, k = data[:, 0].astype(np.int)-1, data[:, 1].astype(np.int)-1, data[:, 2].astype(np.int)-1
    lwc[i, j, k] = data[:, 3]
    reff[i, j, k] = data[:, 4]

    dset['density'] = xr.DataArray(
        data=lwc,
        dims=['x', 'y', 'z']
    )

    dset['reff'] = xr.DataArray(
        data=reff,
        dims=['x', 'y', 'z']
    )

    dset['temperature'] = xr.DataArray(
        data=temperature,
        dims=['z']
    )

    dset.attrs['density_name'] = density
    dset.attrs['file_name'] = file_name

    return dset

def to_2parameter_lwc_file(file_name, cloud_scatterer, atmosphere=None, fill_temperature=280.0):
    """
    TODO
    Write lwc & reff to the '2 parameter lwc' file format used by i3rc MonteCarlo model and SHDOM.
    atmosphere should contain the temperature. It is interpolated to the specified z grid.
    If no atmosphere is included then a fill_temperature is used (Temperature is required
    in the file).
    """

    nx, ny, nz = cloud_scatterer.density.shape
    dx, dy = (cloud_scatterer.x[1]-cloud_scatterer.x[0]).data, (cloud_scatterer.y[1] - cloud_scatterer.y[0]).data
    z = cloud_scatterer.z.data

    if atmosphere is not None:
        temperature = atmosphere.interp({'z': cloud_scatterer.z}).temperature.data
    else:
        temperature = np.ones(z.shape)*fill_temperature

    i, j, k = np.meshgrid(np.arange(1, nx+1), np.arange(1, ny+1), np.arange(1, nz+1), indexing='ij')

    lwc = cloud_scatterer.density.data.ravel()
    reff = cloud_scatterer.reff.data.ravel()

    z_string = ''
    for z_value in z:
        if z_value == z[-1]:
            z_string += '{}'.format(z_value)
        else:
            z_string += '{} '.format(z_value)

    t_string = ''
    for index, temp_value in enumerate(temperature):
        if index == len(temperature) - 1:
            t_string += '{:5.2f}'.format(temp_value)
        else:
            t_string += '{:5.2f} '.format(temp_value)

    with open(file_name, "w") as f:
        f.write('2 parameter LWC file\n')
        f.write(' {} {} {}\n'.format(nx, ny, nz))
        f.write('{} {}\n'.format(dx, dy))
        f.write('{}\n'.format(z_string))
        f.write('{}\n'.format(t_string))
        for x, y, z, l, r in zip(i.ravel(), j.ravel(), k.ravel(), lwc.ravel(), reff.ravel()):
            f.write('{} {} {} {:5.4f} {:3.2f}\n'.format(x, y, z, l, r))

def load_from_csv(path, density=None, origin=(0.0,0.0)):
    """
    TODO
    """
    df = pd.read_csv(path, comment='#', skiprows=4, index_col=['x', 'y', 'z'])
    nx, ny, nz = np.genfromtxt(path, max_rows=1, dtype=int, delimiter=',')
    dx, dy = np.genfromtxt(path, max_rows=1, dtype=float, skip_header=2, delimiter=',')
    z = xr.DataArray(np.genfromtxt(path, max_rows=1, dtype=float, skip_header=3, delimiter=','), coords=[range(nz)], dims=['z'])

    dset = pyshdom.grid.make_grid(dx, nx, dy, ny, z)
    i, j, k = zip(*df.index)

    for name in df.columns:
        #initialize with np.nans so that empty data is np.nan
        variable_data = np.zeros((dset.sizes['x'], dset.sizes['y'], dset.sizes['z']))*np.nan
        variable_data[i, j, k] = df[name]
        dset[name] = (['x', 'y', 'z'], variable_data)

    if density is not None:
        assert density in dset.data_vars, \
        "density variable: '{}' must be in the file".format(density)

        dset = dset.rename_vars({density: 'density'})
        dset.attrs['density_name'] = density

    dset.attrs['file_name'] = path

    return dset

def load_from_netcdf(path, density=None):
    """
        TODO
    A shallow wrapper around open_dataset that sets the density_name.
    """
    dset = xr.open_dataset(path)

    if density is not None:
        if density not in dset.data_vars:
            raise ValueError("density variable: '{}' must be in the file".format(density))
        dset = dset.rename_vars({density: 'density'})
        dset.attrs['density_name'] = density

    dset.attrs['file_name'] = path

    return dset


def load_forward_model(file_name):
    """
    TODO
    """
    dataset = nc.Dataset(file_name)

    groups = dataset.groups
    sensors = groups['sensors'].groups
    solvers = groups['solvers'].groups
    sensor_dict = pyshdom.containers.SensorsDict()
    solver_dict = pyshdom.containers.SolversDict()

    for key,sensor in sensors.items():
        sensor_list = []
        for i, image in sensor.groups.items():

            sensor_dict.add_sensor(key, xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                'sensors/'+str(key)+'/'+str(i)])))
        #     sensor_list.append(xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
        #         'sensors/'+str(key)+'/'+str(i)])))
        # sensor_dict[key] = {'sensor_list':sensor_list}

    for key, solver in solvers.items():
        numerical_params = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                        'solvers/'+str(key)+'/numerical_parameters']))
        num_stokes = numerical_params.num_stokes.data
        surface = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/surface']))
        source = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/source']))

        mediums = OrderedDict()
        for name, med in solver['medium'].groups.items():
            mediums[name] = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                        'solvers/'+str(key)+'/medium/'+str(name)]))

        if 'atmosphere' in solver.groups.keys():
            atmosphere = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                'solvers/'+str(key)+'/atmosphere']))
        else:
            atmosphere=None

        solver_dict.add_solver(float(key), pyshdom.solver.RTE(numerical_params=numerical_params,
                                            medium=mediums,
                                           source=source,
                                           surface=surface,
                                            num_stokes=num_stokes,
                                            name=None
                                           )
                                           )
        rte_grid = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/grid']))

    return sensor_dict, solver_dict, rte_grid
#TODO add checks here for if file exists etc.
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
        for name, med in solver.medium.items():
            med.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'medium/'+str(name), format='NETCDF4',engine='netcdf4')
