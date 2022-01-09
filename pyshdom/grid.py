"""This module contains functions to generate grids that follow the
spatial discretization scheme of SHDOM and to resample data on a different set
of coordinates onto an SHDOM style grid.

SHDOM grids must be equispaced in the 'x' and 'y' coordinates and monotonic and
increasing in the 'z' coordinate. To ensure that these conditions are met, the
'make_grid' function should be used to generate grids, which takes as input a
spacing in the 'x' and 'y' coordinate and number of grid points just as in SHDOM.

When resampling data between two different grids there are many different ways to do this,
the provided functions: 'resample_onto_grid' & 'combine_z_coordinates'/'merge_two_z_coordinates'
assume certain behaviour which is described in their docstrings. If users desire
different behaviour they should write their own functions and add them to this module.
"""
import typing
import numpy as np
import xarray as xr
import pyshdom.checks

import scipy.spatial as ss


def make_grid(delx: float, npx: int, dely: float, npy: int, z: np.ndarray,
              nx=None, ny=None, nz=None) -> xr.Dataset:
    """
    Defines a 3D grid of 'x', 'y', 'z' coordinates according to SHDOM's
    conventions.

    Defines an equispaced grid in the 'x' and 'y' dimensions with spacing `delx`/'dely'
    and number of points `npx`/`npy`, respectively. The starting value of 'x' and 'y'
    is always 0.0. The vertical coordinate, `z` is directly input.

    Parameters
    ----------
    delx, dely : float
        grid spacing in the 'x' and 'y' dimensions, respectively.
    npx, npy : int
        number of grid points in the 'x' and 'y' dimensions, respectively.
    z : array_like of floats
        1D strictly increasing and positive values of the vertical grid points.
    nx, ny, nz : int
        number of grid points for the SHDOM grid. If nz is specified it overrides
        the property grid with an equispaced vertical grid. These variables
        control the numerical resolution used in the SHDOM solution. By default
        these are set to equal npx, npy and the possibly uneven `z` are used in
        the numerical solution.

    Returns
    -------
    grid : xr.Dataset
        Contains coordinates of 'x', 'y' and 'z'. It also contains the 'delx'/'dely'
        variables that are the spacing of the 'x' and 'y' coordinates. Optionally
        contains 'nx', 'ny', 'nz'.

    Raises
    ------
    ValueError
        If the z coordinate is not 1D, positive, strictly increasing and contains
        at least two points.

    Notes
    -----
    More than a single point should generally be defined in the horizontal dimensions
    unless operating in 2D or 1D radiative transfer mode. When operating in those modes,
    the radiative transfer is set to periodic and therefore a single point can be supplied.
    """
    #checks on z to make sure its monotonic.
    z = np.asarray(z)
    if (not np.all(np.sort(z) == z)) or (not np.all(z >= 0.0)) or \
       (np.unique(z).size != z.size) or (z.ndim != 1) or (z.size < 2):
        raise ValueError('z must be >= 0, strictly increasing, 1-D and contain at least 2 points.')

    grid = xr.Dataset(
        coords={
            'x': np.linspace(0.0, delx*(npx-1), npx),
            'y': np.linspace(0.0, dely*(npy-1), npy),
            'z': z,
        }
    )
    grid['delx'] = delx
    grid['dely'] = dely
    if nx is not None:
        grid['nx'] = nx
    if ny is not None:
        grid['ny'] = ny
    if nz is not None:
        grid['nz'] = nz

    return grid

def resample_onto_grid(grid, data):
    """
    Resamples multi-variate data in an xarray Dataset onto the grid
    provided by another xarray Dataset using linear N-dimensional interpolation.

    Missing data (np.nan) are first filled in all variables. If the variable
    is named 'density' then missing data are filled with 0.0. Other variables
    are assumed to not necessarily decay to 0.0 at cloud edge but are instead
    better described by their nearest value with non-zero 'density' so np.nans
    are filled by backwards and forward filled in the 'x', 'y', and then 'z' directions.
    Then all variables are linearly interpolated onto the coordinates in `grid`
    and broadcast to 3D.

    Parameters
    ----------
    grid : xr.Dataset / xr.DataArray
        This contains the coordinates named 'x', 'y', 'z' that data will be
        interpolated onto. It also needs the 'delx' and 'dely' variables that
        correspond to the equispacing of 'x' and 'y'. This input should be
        generated by the 'grid.make_grid' function to avoid errors.
    data : xr.Dataset
        This dataset should contain variables that are to be interpolated onto
        the grid. They should be defined on coordinates named 'x', 'y', 'z'.

    Returns
    -------
    filled : xr.Dataset
        Contains variables from `data` interpolated onto the coordinates of `grid`,
        with prescribed filling of NaN values. Also propagates the 'delx' and 'dely'
        variables from `grid` and attributes from `data`.

    Raises
    ------
    TypeError
        Will be raised if inputs are of the wrong type.
    AssertionError
        If there are (unexpected) np.nan in `filled` before returning.
    exceptions.GridError
        If `grid` is not a valid SHDOM grid. See checks.check_grid.
    """
    if not isinstance(grid, (xr.Dataset, xr.DataArray)):
        raise TypeError("'grid' should be an xr.Dataset, "
                         "xr.DataArray  object, not '{}'".format(type(grid)))
    pyshdom.checks.check_grid(grid)

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

    filled = add_grid_variables(grid, filled)
    return filled

def add_grid_variables(grid, dataset):
    dataset['delx'] = grid.delx
    dataset['dely'] = grid.dely
    if 'nx' in grid.data_vars:
        dataset['nx'] = grid.nx
    if 'ny' in grid.data_vars:
        dataset['ny'] = grid.ny
    if 'nz' in grid.data_vars:
        dataset['nz'] = grid.nz

    return dataset

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
    Merges two arrays of vertical coordinates to provide a combined coordinate.

    Takes two monotonically increasing and positive input z coordinates and combines
    them according to the following rules. In the regions where the coordinates do
    not overlap, the respective coordinate that is defined is used. In the overlapping
    region, the coordinate that has the higher resolution (greater number of points)
    is used.

    Parameters
    ----------
    z1, z2 : array_like of floats.
        1D, monotonically increasing, positive, vertical coordinate.

    Returns
    -------
    combined : np.array of floats.
        The combined coordinate.

    Raises
    ------
    ValueError
        If `z1` or `z2` are not monotonically increasing, positive and 1D.
    AssertionError
        If `combined` has repeated elements or is increasing.
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

def from_scatterer(scatterer):
    """
    Defines an RTE grid from a scatterer object.
    Only covers the simplest case.
    """
    grid = make_grid(scatterer.x.diff('x')[0].data,
                          scatterer.x.size,
                          scatterer.y.diff('y')[0].data,
                          scatterer.y.size,
                          scatterer.z.data)
    resampled_scatterer = resample_onto_grid(grid, scatterer)
    return resampled_scatterer

@xr.register_dataset_accessor("grid")
@xr.register_dataarray_accessor("grid")
class _GridAccessor(object):
    """
    Register a custom accessor for Grid properties particular to xarray Datasets
    and DataArrays used in pyshdom.
    """
    def __init__(self, xarray_obj):
        pyshdom.checks.check_grid(xarray_obj)
        self._obj = xarray_obj

    @property
    def shape(self):
        return (self._obj.x.size, self._obj.y.size, self._obj.z.size)

    def distance_to_clear(self, mask):

        grid = self._obj

        if grid.grid.shape != mask.shape:
            raise ValueError(
                "`grid` and `mask` must be of consistent shape."
            )

        from scipy.spatial import cKDTree
        xs = grid.x.data
        ys = grid.y.data
        zs = grid.z.data
        xs, ys, zs = np.meshgrid(xs, ys, zs, indexing='ij')
        positions = np.stack((xs.ravel(), ys.ravel(), zs.ravel()), axis=1)

        clear_positions = positions[np.where(~mask.ravel()), :][0]
        cloud_positions = positions[np.where(mask.ravel()), :][0]
        tree = ss.cKDTree(clear_positions)
        minimum_distance_to_clear, indices = tree.query(cloud_positions)

        out = np.zeros(mask.shape)
        out[np.where(mask)] = minimum_distance_to_clear

        return out
