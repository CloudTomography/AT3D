import xarray as xr
import numpy as np
import warnings
import shdom
import pandas as pd
import core

@xr.register_dataarray_accessor("shdomSampling")
class shdomSampling(object):
    def __init__(self, data_array):
        self._obj = data_array

    def resample(self, grid, method='linear'):
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


### Load microphysical scatterer
path = '../synthetic_cloud_fields/jpl_les/rico32x37x26.txt'
df = pd.read_csv(path, comment='#', skiprows=4, index_col=['i', 'j', 'k'])
nx, ny, nz = np.genfromtxt(path, max_rows=1, dtype=int, delimiter=',')
dx, dy = np.genfromtxt(path, max_rows=1, dtype=float, skip_header=2, delimiter=',')
z = xr.DataArray(np.genfromtxt(path, max_rows=1, dtype=float, skip_header=3, delimiter=','), coords=[range(nz)], dims=['z'])
scatterer = xr.Dataset.from_dataframe(df)
scatterer = scatterer.rename({'i': 'x', 'j':'y', 'k':'z'})
scatterer = scatterer.assign_coords({'x': dx * scatterer.x, 'y': dy * scatterer.y, 'z': z[scatterer.z]})

scatterer.lwc.isel(x=10).plot()
scatterer.reff.isel(x=10).plot()

###New microphysics checks and _to_size_distribution here.

import itertools

@xr.register_dataset_accessor("shdomMicrophysics")
class shdomMicrophysics(object):
    def __init__(self, data_array):
        assert len(data_array.items()) > 0, 'Dataset must be nonempty.'
        self._obj = data_array

    def _grid_checks(self):
        """
        TODO
        """
        #checks of DataArrays being 'GridData':
        #'x','y','z' coords must be common to
        #all variables and also consistent between variables.

        grid_test = True
        for data_var in self._obj.items():
            for coord in ('x', 'y','z'):
                if coord not in data_var[1].coords:
                    grid_test = False
                    warnings.warn('\'{}\' coord is not in {}. (x, y, z) coords must be present.'.format(data_var[0]))
            for data_var2 in self._obj.items():
                if not all(itertools.product(data_var[1].coords,data_var2[1].coords)):
                    grid_test = False
                    warnings.warn('Variables \'{}\' and \'{}\' do not have consistent coords. All coords must be consistent.'.format(data_var[0],data_var2[0]))
        assert grid_test,'Grid Check Failed.'

    #THIS IS REALLY A GENERATOR FOR DATA SHOULD COME FROM SOMEWHERE ELSE I GUESS.
    def _add_homogeneous_like(self, reference_variable_name, variable_name, value):
        """
        TODO
        """
        self._grid_checks()

        #copy first variable's coordinates etc but not any attributes.
        to_copy = self._obj.data_vars[reference_variable_name]
        new_data = np.full_like(to_copy, fill_value=value)
        new_data_array = xr.DataArray(coords=to_copy.coords,data=new_data)
        new_data_array = new_data_array.where(np.bitwise_not(np.isnan(to_copy)),np.nan)
        self._obj[variable_name] = new_data_array


    def microphysics_checks(self, replace_defaults=False, default_distribution_type=None, default_veff=None):
        """
        TODO
        """
        try:
            distribution_type = self._obj.attrs['distribution_type']
        except KeyError:
            if replace_defaults:
                if default_distribution_type in ('gamma','lognormal'):
                    self._obj.attrs['distribution_type'] = default_distribution_type
                    distribution_type = self._obj.attrs['distribution_type']
                    warnings.warn('No \'distribution_type\' attribute specified. Replacing with default \'{}\''.format(default_distribution_type))
                else:
                    raise ValueError("A valid \'default_distribution_type\' must be specified. Supported values are \'gamma\' and \'lognormal\'")
            else:
                raise KeyError('A \'distribution_type\' attribute must be specified. Supported values are \'gamma\' and \'lognormal\'.')

        if distribution_type in ('gamma','lognormal'):
            try:
                lwc = self._obj['lwc']
                assert np.all(lwc.where(np.bitwise_not(np.isnan(lwc)),0.0) >=0.0), 'lwc data must be positive semidefinite.'
            except:
                raise

            try:
                reff = self._obj['reff']
                assert np.all(reff.where(np.bitwise_not(np.isnan(reff)),0.0) >=0.0), 'reff data must be positive definite.'
            except:
                raise

            try:
                veff = self._obj['veff']
                assert np.all(veff.where(np.bitwise_not(np.isnan(veff)),0.0) >=0.0), 'veff data must be positive definite.'
            except KeyError:
                try:
                    alpha = self._obj['alpha']
                    assert np.all(alpha.where(np.bitwise_not(np.isnan(alpha)),0.0) >=0.0), 'alpha data must be positive definite.'

                except KeyError:
                    if replace_defaults:
                        assert default_veff > 0.0, 'default_veff must be positive definite.'
                        self._add_homogeneous_like('lwc','veff', default_veff)
                        warnings.warn('No \'alpha\' or \'veff\' specified for distribution_type \'{}\'. Replacing with homogeneous default value \'veff\'={}'.format(distribution_type,default_veff))

                    else:
                        raise
        else:
            raise NotImplementedError("Supported '\distribution_type'\ are \'gamma\' and \'lognormal\'")



    def _prepare_size_distribution(self, particle_density=1.0, size_distribution_mode='GRID',max_grid_phase=1000,Ntable=50,
                              default_distribution_type='gamma', default_veff=0.1, replace_defaults=False):

        self._grid_checks()
        self.microphysics_checks(replace_defaults=replace_defaults, default_distribution_type=default_distribution_type, default_veff=default_veff)

        #in case more types are implemented I explicitly separate.
        if self._obj.attrs['distribution_type'] in ('gamma','lognormal'):
            reff = self._obj['reff'].values.ravel()
            gamma=np.zeros(reff.shape)

            if self._obj.attrs['distribution_type'] == 'gamma':
                if 'veff' in self._obj:
                    alpha = (1.0/self._obj['veff'] - 3.0).values.ravel()
                else:
                    alpha = self._obj['alpha'].values.ravel()
            elif self._obj.attrs['distribution_type'] == 'lognormal':
                if 'veff' in self._obj:
                    alpha = np.sqrt(np.log(self._obj['veff'] + 1.0)).values.ravel()
                else:
                    alpha = self._obj['alpha'].values.ravel()

            """
            TODO BELOW
            """
            #Get Reff & Alpha coordinates for table depending whether its 'GRID' or 'TABULATED'.
            if size_distribution_mode == 'GRID':

                unique_pairs = list(set([(i,j) for i,j in zip(reff,alpha)])) #smarter way needed.

                if len(unique_pairs) > max_grid_phase:
                    warnings.warn('Number of phase functions \'{}\' is larger than suggested number {}'.format(len(unique_pairs),max_grid_phase))

                #only take the unique elements.
                reff_table, alpha_table = zip(*unique_pairs)


            elif size_distribution_mode == 'TABULATED':
                assert Ntable > 0
                reff_table, alpha_table = np.meshgrid(np.arange(reff.min(),reff.max(),Ntable),
                                                   np.arange(alpha.min(),alpha.max(),Ntable))

        else:
            raise NotImplementedError("Supported '\distribution_type'\ are \'gamma\' and \'lognormal\'")

        self._obj.assign_coords({'table_reff':reff_table,
                                 'table_alpha': alpha_table})

#checks can be run independently.
scatterer.shdomMicrophysics.microphysics_checks(replace_defaults=True,default_distribution_type='gamma',default_veff=0.1)

#this definitely doesn't work. One issue is deciding on masking convention - while masked values always be np.nan?
#what about missing, or invalid data?
scatterer.shdomMicrophysics._prepare_size_distribution()
