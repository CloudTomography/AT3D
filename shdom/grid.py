"""
Utility functions to create grids and resample scatterers onto them.
TODO the behaviour of ALL functions still needs to be verified.
All functions other than 'resample_onto_grid' and 'make_grid' are
designed purely to reproduce 'pre-refactoring' grid merging.

In general, RTE Grids can be arbitrarily defined using the 'make_grid' function
and scatterers forced to resample onto the specified grid.

"""
import pandas as pd
import numpy as np
import xarray as xr
import typing

def load_2parameter_lwc_file(file_name, density='lwc', origin=(0.0,0.0)):
    """
    Function that loads a scatterer from the '2 parameter lwc file' format used by
    SHDOM and i3rc monte carlo model.
    """
    header = pd.read_csv(file_name, nrows=4)
    nx,ny,nz = np.fromstring(header['2 parameter LWC file'][0],sep=' ').astype(np.int)
    dx,dy = np.fromstring(header['2 parameter LWC file'][1],sep=' ').astype(np.float)
    z = np.fromstring(header['2 parameter LWC file'][2],sep=' ').astype(np.float)
    temperature = np.fromstring(header['2 parameter LWC file'][3],sep=' ').astype(np.float)

    dset = make_grid(origin[0],(nx-1)*dx,nx,origin[1],(ny-1)*dy,ny,z)

    data = np.genfromtxt(file_name,skip_header=5)

    lwc = np.zeros((nx,ny,nz))*np.nan
    reff = np.zeros((nx,ny,nz))*np.nan

    lwc[data[:,0].astype(np.int)-1,data[:,1].astype(np.int)-1,data[:,2].astype(np.int)-1] = data[:,3]
    reff[data[:,0].astype(np.int)-1,data[:,1].astype(np.int)-1,data[:,2].astype(np.int)-1] = data[:,4]

    dset['density'] = xr.DataArray(
                        data=lwc,
                        dims=['x','y','z']
    )

    dset['reff'] = xr.DataArray(
                        data=reff,
                        dims=['x','y','z']
    )

    dset['temperature'] = xr.DataArray(
                        data=temperature,
                        dims=['z']
    )

    dset.attrs['density_name'] = 'lwc'
    dset.attrs['file_name'] = file_name

    return dset

def to_2parameter_lwc_file(file_name, cloud_scatterer,atmosphere=None,fill_temperature=280.0):
    """
    Write lwc & reff to the '2 parameter lwc' file format used by i3rc MonteCarlo model and SHDOM.
    atmosphere should contain the temperature. It is interpolated to the specified z grid.
    If no atmosphere is included then a fill_temperature is used (Temperature is required
    in the file).
    """

    nx,ny,nz = cloud_scatterer.density.shape
    dx,dy = (cloud_scatterer.x[1]-cloud_scatterer.x[0]).data, (cloud_scatterer.y[1] - cloud_scatterer.y[0]).data
    z = cloud_scatterer.z.data

    if atmosphere is not None:
        temperature = atmosphere.interp({'z': cloud_scatterer.z}).temperature.data
    else:
        temperature = np.ones(z.shape)*fill_temperature

    i,j,k = np.meshgrid(np.arange(1,nx+1),np.arange(1,ny+1),np.arange(1,nz+1), indexing='ij')

    lwc = cloud_scatterer.density.data.ravel()
    reff = cloud_scatterer.reff.data.ravel()

    z_string = ''
    for a in z:
        if a == z[-1]:
            z_string += '{}'.format(a)
        else:
            z_string += '{} '.format(a)

    t_string = ''
    for iter,a in enumerate(temperature):
        if iter == len(temperature) - 1:
            t_string += '{:5.2f}'.format(a)
        else:
            t_string +='{:5.2f} '.format(a)

    with open(file_name, "w") as f:
        f.write('2 parameter LWC file\n')
        f.write(' {} {} {}\n'.format(nx,ny,nz))
        f.write('{} {}\n'.format(dx,dy))
        f.write('{}\n'.format(z_string))
        f.write('{}\n'.format(t_string))
        for x,y,z,l,r in zip(i.ravel(),j.ravel(),k.ravel(),lwc.ravel(),reff.ravel()):
            f.write('{} {} {} {:5.4f} {:3.2f}\n'.format(x, y, z,l,r))

def load_from_csv(path, density=None,origin=(0.0,0.0)):
    """
    TODO
    """
    df = pd.read_csv(path, comment='#', skiprows=4, index_col=['x', 'y', 'z'])
    nx, ny, nz = np.genfromtxt(path, max_rows=1, dtype=int, delimiter=',')
    dx, dy = np.genfromtxt(path, max_rows=1, dtype=float, skip_header=2, delimiter=',')
    z = xr.DataArray(np.genfromtxt(path, max_rows=1, dtype=float, skip_header=3, delimiter=','), coords=[range(nz)], dims=['z'])

    dset = make_grid(origin[0],(nx-1)*dx,nx,origin[1],(ny-1)*dy,ny,z)
    i,j,k = zip(*df.index)

    for name in df.columns:
        #initialize with np.nans so that empty data is np.nan
        variable_data = np.zeros((dset.sizes['x'],dset.sizes['y'],dset.sizes['z']))*np.nan
        variable_data[i,j,k] = df[name]
        dset[name] = (['x', 'y', 'z'], variable_data)

    if density is not None:
        assert density in dset.data_vars, "density variable: '{}' must be in the file".format(density)
        dset = dset.rename_vars({density: 'density'})
        dset.attrs['density_name'] = density

    dset.attrs['file_name'] = path

    return dset

def load_from_netcdf(path, density=None):
    """
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

def resample_onto_grid(grid, data):
    """
    TODO
    """
    if not isinstance(grid, (xr.core.coordinates.DatasetCoordinates,xr.core.coordinates.DataArrayCoordinates, xr.Dataset, xr.DataArray)):
        raise ValueError("'grid' should be an xr.Dataset, xr.DataArray or xarray coordinates object, not '{}'".format(type(grid)))
    if not isinstance(data, (xr.Dataset, xr.DataArray)):
        raise ValueError("'grid' should be an xr.Dataset, xr.DataArray, not '{}'".format(type(grid)))
    #if coordinates are passed, then make a dataset.
    if isinstance(grid, (xr.core.coordinates.DatasetCoordinates,
                                xr.core.coordinates.DataArrayCoordinates)):
        grid = xr.Dataset(
                    coords=grid
        )
    #linearly interpolate onto the grid.
    #data is broadcasted to 3D.
    data_copy = data.copy(deep=True)
    if 'density' in data:
        data_copy['density'] = data.density.fillna(0.0)
    for name,var in data.data_vars.items():
        if name != 'density':
            if ('x' in var.coords) & ('y' in var.coords) & ('z' in var.coords):
                data_copy[name] = var.bfill(dim='x').ffill(dim='x').bfill(dim='y').ffill(dim='y').bfill(dim='z').ffill(dim='z')
            elif ('x' in var.coords) & ('y' in var.coords):
                data_copy[name] = var.bfill(dim='x').ffill(dim='x').bfill(dim='y').ffill(dim='y')
            elif ('x' in var.coords) & ('z' in var.coords):
                data_copy[name] = var.bfill(dim='x').ffill(dim='x').bfill(dim='z').ffill(dim='z')
            elif ('y' in var.coords) & ('z' in var.coords):
                data_copy[name] = var.bfill(dim='y').ffill(dim='y').bfill(dim='z').ffill(dim='z')
            elif ('x' in var.coords):
                data_copy[name] = var.bfill(dim='x').ffill(dim='x')
            elif ('y' in var.coords):
                data_copy[name] = var.bfill(dim='y').ffill(dim='y')
            elif ('z' in var.coords):
                data_copy[name] = var.bfill(dim='z').ffill(dim='z')

    resampled_data = data_copy.interp_like(grid, method='linear').broadcast_like(grid)
    filled = resampled_data.bfill(dim='x').ffill(dim='x').bfill(dim='y').ffill(dim='y').bfill(dim='z').ffill(dim='z')

    #overwrite density values so missing data is filled with 0.0
    if 'density' in filled:
        filled['density'] = resampled_data.density.fillna(0.0)
    for name, datavar in filled.data_vars: #consistency check.
        assert np.bitwise_not(np.all(np.isnan(datavar.data))), "Unexpected NaN in '{}'".format(name)
    return filled

def make_grid(xmax,nx,ymax,ny,z):
    """
    TODO
    """
    #TODO checks on z to make sure its monotonic.
    z = np.asarray(z)
    if (not np.all(np.sort(z) == z)) or (not np.all(z >=0.0)) or (np.unique(z1).size != z.size) or (z.ndim != 1):
        raise ValueError('z must be >= 0 and strictly increasing.')
    return xr.Dataset(
        coords = {
            'x': np.linspace(xmin,xmax,nx),
            'y': np.linspace(ymin,ymax,ny),
            'z': z,
        }
    )
def combine_z_coordinates(scatterer_list):
    """
    A wrapper around merge_two_z_coordinates.
    """
    if not isinstance(scatterer_list, (typing.List, typing.Tuple)):
        raise TypeError("scatterer_list should be a Tuple or List not '{}''".format(type(scatterer_list)))
    for item in scatterer_list:
        if not isinstance(item, (xr.Dataset, xr.DataArray)):
            raise TypeError("Elements of 'scatterer_list' should be xr.Dataset or xr.DataArray not '{}'".format(type(item)))

    z_coordinate_list = [scatterer.coords['z'].data for scatterer in scatterer_list]
    if len(z_coordinate_list) == 1:
        combined = z_coordinate_list[0]
    else:
        combined = merge_two_z_coordinates(z_coordinate_list[0],z_coordinate_list[1])
        for z_coord in z_coordinate_list[2:]:
            combined = merge_two_z_coordinates(z_coord,combined)
    assert np.unique(combined).size == combined.size, 'unexpected repeated elements.'
    assert np.all(np.sort(combined) == combined), 'unexpectedly not strictly increasing.'
    return combined

def merge_two_z_coordinates(z1,z2):
    """
    TODO
    """
    z1 = np.asarray(z1)
    z2 = np.asarray(z2)

    if (not np.all(np.sort(z1) == z1)) or (not np.all(z1 >=0.0)) or (np.unique(z1).size != z1.size) or (z1.ndim != 1):
        raise ValueError('z1 must be >= 0 and strictly increasing.')
    if (not np.all(np.sort(z2) == z2)) or (not np.all(z2 >=0.0)) or (np.unique(z2).size != z2.size) or (z2.ndim != 1):
        raise ValueError('z1 must be >= 0 and strictly increasing.')

    # Bottom part of the atmosphere (no grid intersection)
    z_bottom = z1[z1 < z2[0]] if z1[0] < z2[0] else z2[z2 < z1[0]]

    # Top part of the atmosphere (no grid intersection)
    z_top = z2[z2 > z1[-1]] if z1[-1] < z2[-1] else z1[z1 > z2[-1]]

    # Middle part of the atmosphere (grids intersect)
    z1_middle = z1
    z2_middle = z2
    if (z_bottom.size > 0) & (z_top.size > 0):
        z1_middle = z1[(z1 > z_bottom[-1]) & (z1 < z_top[0])]
        z2_middle = z2[(z2 > z_bottom[-1]) & (z2 < z_top[0])]
    elif z_top.size > 0:
        z1_middle = z1[z1 < z_top[0]]
        z2_middle = z2[z2 < z_top[0]]
    elif z_bottom.size > 0:
        z1_middle = z1[z1 > z_bottom[-1]]
        z2_middle = z2[z2 > z_bottom[-1]]
    #pick the higher resolution middle based no number of points in
    #region of grid intersection.
    #
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

    combined = np.concatenate((z_bottom, z_middle, z_top))
    assert np.unique(combined).size == combined.size, 'unexpected repeated elements.'
    assert np.all(np.sort(combined) == combined), 'unexpectedly not strictly increasing.'

    return combined


#
# def find_horizontal_union(scatterer_list):
#     """
#     TODO
#     """
#     x_mins = []
#     x_maxs = []
#     y_mins = []
#     y_maxs = []
#     for scatterer in scatterer_list:
#         try: #TODO INCORRECT need to test for existence of different coords.
#             x_mins.append(scatterer.coords['x'].data.min())
#             y_mins.append(scatterer.coords['y'].data.min())
#             x_maxs.append(scatterer.coords['x'].data.max())
#             y_maxs.append(scatterer.coords['y'].data.max())
#         except:
#             continue
#
#     assert len(x_mins) > 0, 'At least one scatterer must have an x dimension'
#     assert len(y_mins) > 0, 'At least one scatterer must have a y dimension'
#
#     return min(x_mins), max(x_maxs), min(y_mins),max(y_maxs)

# def find_max_horizontal_resolution(scatterer_list):
#     """
#     TODO
#     """
#     dx_mins = []
#     dy_mins = []
#     for scatterer in scatterer_list:
#         try: #TODO INCORRECT need to test for existence of different coords.
#             dx_mins.append(scatterer.coords['x'].diff(dim='x').data.min())
#             dy_mins.append(scatterer.coords['y'].diff(dim='y').data.min())
#         except:
#             continue
#
#     assert len(dx_mins) > 0, 'At least one scatterer must have an x dimension'
#     assert len(dy_mins) > 0, 'At least one scatterer must have a y dimension'
#
#     return min(dx_mins), min(dy_mins)

# def merge_scatterer_grids(scatterer_list):
#     """
#     A function that produces the grid addition behaviour
#     of scatterer merging from the 'pre-refactoring' shdom.
#     TODO
#     """
#     assert len(scatterer_list) > 1,'need more than 1 scatterer to combine.'
#
#     xmin,xmax,ymin,ymax = find_horizontal_union(scatterer_list)
#     dx,dy = find_max_horizontal_resolution(scatterer_list)
#
#     #defined so that linspace includes xmax as the last point.
#     #supports unevenly spaced x,y in input data.
#     nx = int(1+(xmax-xmin)/dx)
#     ny = int(1+(ymax-ymin)/dy)
#
#     z = combine_z_coordinates(scatterer_list)
#     grid = make_grid(xmin,xmax,nx,ymin,ymax,ny,z)
#     merged_scatterers = [resample_onto_grid(grid,scatterer) for scatterer in scatterer_list]
#
#     return merged_scatterers
