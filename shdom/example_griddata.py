import xarray as xr
import numpy as np
import warnings
import shdom
import pandas as pd

@xr.register_dataarray_accessor("GridData")
class GridData(object):
    def __init__(self, data, coords, dims, grid, method='linear',name=None):

        data_array = xr.DataArray(data=data,coords=coords,dims=dims,name=name)

        temp_data = np.zeros((len(grid['x']),len(grid['y']),len(grid['z'])))
        temp_grid_data = xr.DataArray(data=temp_data,coords=[grid['x'],grid['y'],grid['z']],
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

#lwc = xr.DataArray(name='lwc',data=data,coords=[xs,ys,zs], dims=['x','y','z'])
#reff = xr.DataArray(name='reff',data=data_1d,coords=[zs2],dims=['z'])

# Air grid
z_a = np.arange(0, 20, 1.0)
z_c = np.arange(0.5, 1.5, 0.04)
# Cloud grid
x = np.linspace(0.0,1.0,15)
y = np.linspace(0.0,1.0,16)
z = shdom.combine_1d_grids(z_c, z_a)
# Atmosphere grid
grid = pd.Series(index=['x','y','z'], data=[x, y, z])

lwc = xr.DataArray.GridData(data=lwc,grid=grid,coords=[xs,ys,zs],dims=['x','y','z'], name='lwc')
reff = xr.DataArray.GridData(data=reff,grid=grid,coords=[zs2],dims=['z'],name='reff')
