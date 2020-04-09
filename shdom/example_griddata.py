import xarray as xr
import numpy as np
import warnings
import shdom
import pandas as pd

@xr.register_dataarray_accessor("shdomSampling")
class shdomSampling(object):
    def __init__(self, data_array):
        self._obj = data_array

    def resample(self, grid, inplace=True, method='linear'):
        temp_data = np.zeros((len(grid['x']), len(grid['y']), len(grid['z'])))
        temp_grid_data = xr.DataArray(
            data=temp_data, coords=[grid['x'], grid['y'], grid['z']], dims=['x', 'y', 'z']
        )

        for key, coord in grid.items():
            if key in self._obj.coords:
                if coord.max() > self._obj.coords[key].max():
                    warnings.warn('{} coordinate exceeded the grid maximum. Data has been truncated.'.format(key))
                if coord.min() < self._obj.coords[key].max():
                    warnings.warn('{} coordinate exceeded the grid minimum. Data has been truncated.'.format(key))

        interped = self._obj.interp_like(temp_grid_data, method=method)
        resampled, temp = xr.broadcast(interped, temp_grid_data)

        if inplace:
            self._obj = resampled
        else:
            return resampled


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

lwc = xr.DataArray(data=data,coords=[xs,ys,zs], dims=['x','y','z'], name='lwc')
reff = xr.DataArray(data=data_1d, coords=[zs2], dims=['z'], name='reff')

# Resampling inplace
lwc.shdomSampling.resample(grid)

# Return new object
reff = reff.shdomSampling.resample(grid, inplace=False)