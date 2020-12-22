"""
Utility functions to create grids and resample scatterers onto them.
TODO the behaviour of ALL functions still needs to be verified.
All functions other than 'resample_onto_grid' and 'make_grid' are
designed purely to reproduce 'pre-refactoring' grid merging.

In general, RTE Grids can be arbitrarily defined using the 'make_grid' function
and scatterers forced to resample onto the specified grid.

"""
import typing
import numpy as np
import xarray as xr

def make_grid(delx: float, nx: int, dely: float, ny: int, z: np.ndarray) -> xr.Dataset:
    """
    TODO
    """
    #checks on z to make sure its monotonic.
    z = np.asarray(z)
    if (not np.all(np.sort(z) == z)) or (not np.all(z >= 0.0)) or \
       (np.unique(z).size != z.size) or (z.ndim != 1) or (z.size < 2):
        raise ValueError('z must be >= 0, strictly increasing, 1-D and contain at least 2 points.')

    grid = xr.Dataset(
        coords={
            'x': np.linspace(0.0, delx*(nx-1), nx),
            'y': np.linspace(0.0, dely*(ny-1), ny),
            'z': z,
        }
    )
    grid['delx'] = delx
    grid['dely'] = dely
    return grid

def resample_onto_grid(grid, data):
    """
    TODO
    """
    if not isinstance(grid, (xr.core.coordinates.DatasetCoordinates,
                             xr.core.coordinates.DataArrayCoordinates, xr.Dataset, xr.DataArray)):
        raise ValueError("'grid' should be an xr.Dataset, "
                         "xr.DataArray or xarray coordinates object, not '{}'".format(type(grid)))
    if not isinstance(data, (xr.Dataset, xr.DataArray)):
        raise ValueError("'grid' should be an xr.Dataset, "
                         "xr.DataArray, not '{}'".format(type(grid)))
    #if coordinates are passed, then make a dataset.
    if isinstance(grid, (xr.core.coordinates.DatasetCoordinates,
                         xr.core.coordinates.DataArrayCoordinates)):
        grid = xr.Dataset(
            coords=grid
        )
    data_copy = data.copy(deep=True)
    if 'density' in data:
        data_copy['density'] = data.density.fillna(0.0)

    # all variables except 'density'shouldn't be filled with zero at undefined points
    # as this can cause very different values from expected in the following linear interpolation.
    # Note that the following choice of backward and forward filling in (x,y,z) order is subjective.
    # The validity of this method relies on the assumption that microphysics don't decay towards cloud edge and instead maintain
    # a typical value. 'z' is always filled last as microphysics are expected to vary strongly in the vertical.
    for name, var in data.data_vars.items():
        if name != 'density':
            if ('x' in var.coords) & ('y' in var.coords) & ('z' in var.coords):
                data_copy[name] = var.bfill(dim='x').ffill(dim='x').bfill(
                    dim='y').ffill(dim='y').bfill(dim='z').ffill(dim='z')
            elif ('x' in var.coords) & ('y' in var.coords):
                data_copy[name] = var.bfill(dim='x').ffill(dim='x').bfill(dim='y').ffill(dim='y')
            elif ('x' in var.coords) & ('z' in var.coords):
                data_copy[name] = var.bfill(dim='x').ffill(dim='x').bfill(dim='z').ffill(dim='z')
            elif ('y' in var.coords) & ('z' in var.coords):
                data_copy[name] = var.bfill(dim='y').ffill(dim='y').bfill(dim='z').ffill(dim='z')
            elif 'x' in var.coords:
                data_copy[name] = var.bfill(dim='x').ffill(dim='x')
            elif 'y' in var.coords:
                data_copy[name] = var.bfill(dim='y').ffill(dim='y')
            elif 'z' in var.coords:
                data_copy[name] = var.bfill(dim='z').ffill(dim='z')

    resampled_data = data_copy.interp_like(grid, method='linear').broadcast_like(grid)
    #for variables which weren't defined at every point in the rte_grid, perform filling.
    filled = resampled_data.bfill(dim='x').ffill(dim='x').bfill(
                dim='y').ffill(dim='y').bfill(dim='z').ffill(dim='z')

    #overwrite density values so missing data is filled with 0.0
    if 'density' in filled:
        filled['density'] = resampled_data.density.fillna(0.0)
    for name, datavar in filled.data_vars.items(): #consistency check.
        assert np.bitwise_not(np.all(np.isnan(datavar.data))), "Unexpected NaN in '{}'".format(name)

    filled['delx'] = grid.delx
    filled['dely'] = grid.dely
    return filled

def combine_z_coordinates(scatterer_list):
    """
    A wrapper around merge_two_z_coordinates.
    """
    if not isinstance(scatterer_list, (typing.List, typing.Tuple)):
        raise TypeError("scatterer_list should be a Tuple "
                        "or List not '{}''".format(type(scatterer_list)))
    for item in scatterer_list:
        if not isinstance(item, (xr.Dataset, xr.DataArray)):
            raise TypeError("Elements of 'scatterer_list' should be "
                            "xr.Dataset or xr.DataArray not '{}'".format(type(item)))

    z_coordinate_list = [scatterer.coords['z'].data for scatterer in scatterer_list]
    if len(z_coordinate_list) == 1:
        combined = z_coordinate_list[0]
    else:
        combined = merge_two_z_coordinates(z_coordinate_list[0], z_coordinate_list[1])
        for z_coord in z_coordinate_list[2:]:
            combined = merge_two_z_coordinates(z_coord, combined)
    assert np.unique(combined).size == combined.size, 'unexpected repeated elements.'
    assert np.all(np.sort(combined) == combined), 'unexpectedly not strictly increasing.'
    return combined

def merge_two_z_coordinates(z1, z2):
    """
    TODO change variable names?
    """
    z1 = np.asarray(z1)
    z2 = np.asarray(z2)

    if (not np.all(np.sort(z1) == z1)) or (not np.all(z1 >= 0.0)) or \
        (np.unique(z1).size != z1.size) or (z1.ndim != 1):
        raise ValueError('z1 must be >= 0,strictly increasing and 1-D')
    if (not np.all(np.sort(z2) == z2)) or (not np.all(z2 >= 0.0)) or \
        (np.unique(z2).size != z2.size) or (z2.ndim != 1):
        raise ValueError('z2 must be >= 0,strictly increasing and 1-D')

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
    if z_bottom.any() & len(z_middle) > 2:
        extra_zlevel = 2*z_middle[0] - z_middle[1]
        if extra_zlevel > z_bottom[-1]:
            z_middle = np.append(extra_zlevel, z_middle)

    # Check if an extra point is necessary at the top
    if z_top.any() & len(z_middle) > 2:
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
