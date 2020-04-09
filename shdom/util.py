"""
Utility functions
"""
import numpy as np

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

def combine_1d_grids(z1, z2):
    """
    Find a common grid (union of the grids) which maintains the high resolution grid.

    Parameters
    ----------
    z1: np.array(dtype=np.float32)
        Input 1D grid
    z2: np.array(dtype=np.float32)
        Input 1D grid

    Returns
    -------
    output_grid: np.array(dtype=np.float32)
        The common grid.
    """
    print('Combining multi-resolution grids...')
    # Bottom part of the atmosphere (no grid intersection)
    z_bottom = z1[z1 < z2[0]] if z1[0] < z2[0] else z2[z2 < z1[0]]

    # Top part of the atmosphere (no grid intersection)
    z_top = z2[z2 > z1[-1]] if z1[-1] < z2[-1] else z1[z1 > z2[-1]]

    # Middle part of the atmosphere (grids intersect)
    z1_middle = z1
    z2_middle = z2
    if z_bottom.any():
        z1_middle = z1[z1 > z_bottom[-1]]
        z2_middle = z2[z2 > z_bottom[-1]]
    if z_top.any():
        z1_middle = z1[z1 < z_top[0]]
        z2_middle = z2[z2 < z_top[0]]

    z_middle = z1_middle if len(z1_middle) > len(z2_middle) else z2_middle

    # Check if an extra point is necessary at the bottom
    if z_bottom.any() & len(z_middle)>2:
        extra_zlevel = 2*z_middle[0] - z_middle[1]
        if extra_zlevel > z_bottom[-1]:
            z_middle = np.append(extra_zlevel, z_middle)

    # Check if an extra point is necessary at the top
    if z_top.any() & len(z_middle)>2:
        extra_zlevel = 2*z_middle[-1] - z_middle[-2]
        if extra_zlevel < z_top[0]:
            z_middle = np.append(z_middle, extra_zlevel)

    return np.concatenate((z_bottom, z_middle, z_top))


#WORKING BELOW HERE
import xarray as xr
import numpy as np
import warnings
import shdom
import pandas as pd

@xr.register_dataarray_accessor("GridData")
class GridData(object):
    def __init__(self, data_array, grid, method='linear'):

        temp_data = np.zeros((len(grid.x),len(grid.y),len(grid.z)))
        temp_grid_data = xr.DataArray(data=temp_data,coords=[grid.x,grid.y,grid.z],
                                     dims=['x','y','z'])

        for key,coord in grid.items():
            if key in data_array.coords:
                if coord.max() > data_array.coords[key].max():
                    warnings.warn('{} coordinate exceeded the grid maximum. Data has been truncated.'.format(key))
                if coord.min() < data_array.coords[key].max():
                    warnings.warn('{} coordinate exceeded the grid minimum. Data has been truncated.'.format(key))

        interped = data_array.interp_like(temp_grid_data,method=method)
        resampled, temp = xr.broadcast(interped, temp_grid_data)

        self._obj = resampled

    def nx(self):
        return len(self._obj.coords['x'])

    def ny(self):
        return len(self._obj.coords['y'])

    def nz(self):
        return len(self._obj.coords['z'])

#make some data
xs = np.linspace(0,1.0,15)
ys = np.linspace(0,1.0,16)
zs = np.linspace(0,1.0,17)

zs2 = np.linspace(0,1.2,32)

data = np.random.normal(size=(len(xs),len(ys),len(zs)))
data_1d = np.arange(len(zs2))

lwc = xr.DataArray(name='lwc',data=data,coords=[xs,ys,zs], dims=['x','y','z'])
reff = xr.DataArray(name='reff',data=data_1d,coords=[zs2],dims=['z'])

# Air grid
z_a = np.arange(0, 20, 1.0)
z_c = np.arange(0.5, 1.5, 0.04)
# Cloud grid
x = np.linspace(0.0,1.0,15)
y = np.linspace(0.0,1.0,16)
z = shdom.combine_1d_grids(z_c, z_a)
# Atmosphere grid
grid = pd.Series(index=['x','y','z'], data=[x, y, z])

grid_data = xr.DataArray.GridData(lwc,grid)
