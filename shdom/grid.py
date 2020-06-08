"""
Utility functions to create grids and resample scatterers onto them.
TODO the behaviour of ALL functions still needs to be verified.
All functions other than 'resample_onto_grid' and 'make_grid' are
designed purely to reproduce 'pre-refactoring' grid merging.

In general, RTE Grids can be arbitrarily defined using the 'make_grid' function
and scatterers forced to resample onto the specified grid.

"""

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
        assert density in dset.data_vars, "density variable: '{}' must be in the file".format(density)
        dset = dset.rename_vars({density: 'density'})
        dset.attrs['density_name'] = density

    dset.attrs['file_name'] = path

    return dset


def resample_onto_grid(grid, data):
    """
    TODO
    """
    #if coordinates are passed, then make a dataset.
    if isinstance(grid, (xr.core.coordinates.DatasetCoordinates,
                                xr.core.coordinates.DataArrayCoordinates)):
        grid = xr.Dataset(
                    coords=grid
        )
    #linearly interpolate onto the grid.
    #data is broadcasted to 3D.
    resampled_data = data.interp_like(grid).broadcast_like(grid)

    #fill all missing values with unlimited forward and backward filling.
    filled = resampled_data.bfill(dim='x').ffill(dim='x').bfill(dim='y').ffill(dim='y').bfill(dim='z').ffill(dim='z')

    #overwrite density values so missing data is filled with 0.0
    filled['density'] = resampled_data.density.fillna(0.0)

    return filled

def make_grid(xmin,xmax,nx,ymin,ymax,ny,z):
    """
    TODO
    """
    #TODO checks on z to make sure its monotonic.

    return xr.Dataset(
        coords = {
            'x': np.linspace(xmin,xmax,nx),
            'y': np.linspace(ymin,ymax,ny),
            'z': z,
        }
    )

def find_horizontal_union(scatterer_list):
    """
    TODO
    """
    x_mins = []
    x_maxs = []
    y_mins = []
    y_maxs = []
    for scatterer in scatterer_list:
        try: #TODO INCORRECT need to test for existence of different coords.
            x_mins.append(scatterer.coords['x'].data.min())
            y_mins.append(scatterer.coords['y'].data.min())
            x_maxs.append(scatterer.coords['x'].data.max())
            y_maxs.append(scatterer.coords['y'].data.max())
        except:
            continue

    assert len(x_mins) > 0, 'At least one scatterer must have an x dimension'
    assert len(y_mins) > 0, 'At least one scatterer must have a y dimension'

    return min(x_mins), max(x_maxs), min(y_mins),max(y_maxs)

def find_max_horizontal_resolution(scatterer_list):
    """
    TODO
    """
    dx_mins = []
    dy_mins = []
    for scatterer in scatterer_list:
        try: #TODO INCORRECT need to test for existence of different coords.
            dx_mins.append(scatterer.coords['x'].diff(dim='x').data.min())
            dy_mins.append(scatterer.coords['y'].diff(dim='y').data.min())
        except:
            continue

    assert len(dx_mins) > 0, 'At least one scatterer must have an x dimension'
    assert len(dy_mins) > 0, 'At least one scatterer must have a y dimension'

    return min(dx_mins), min(dy_mins)

def merge_scatterer_grids(scatterer_list):
    """
    A function that produces the grid addition behaviour
    of scatterer merging from the 'pre-refactoring' shdom.
    TODO
    """
    assert len(scatterer_list) > 1,'need more than 1 scatterer to combine.'

    xmin,xmax,ymin,ymax = find_horizontal_union(scatterer_list)
    dx,dy = find_max_horizontal_resolution(scatterer_list)

    #defined so that linspace includes xmax as the last point.
    #supports unevenly spaced x,y in input data.
    nx = int(1+(xmax-xmin)/dx)
    ny = int(1+(ymax-ymin)/dy)

    z = combine_z_coordinates(scatterer_list)
    grid = make_grid(xmin,xmax,nx,ymin,ymax,ny,z)
    merged_scatterers = [resample_onto_grid(grid,scatterer) for scatterer in scatterer_list]

    return merged_scatterers


def combine_z_coordinates(scatterer_list):
    """
    A wrapper around merge_two_z_coordinates.
    """
    assert len(scatterer_list) > 1,'need more than 1 z coordinate to combine.'

    z_coordinate_list = [scatterer.coords['z'].data for scatterer in scatterer_list]

    combined = merge_two_z_coordinate(z_coordinate_list[0],z_coordinate_list[1])
    for scatterer in scatterer_list[2:]:
        combined = merge_two_z_coordinates(scatterer,combined)

    return combined

def merge_two_z_coordinates(z1,z2):
    """
    TODO
    """
    # Bottom part of the atmosphere (no grid intersection)
    z_bottom = z1[z1 < z2[0]] if z1[0] < z2[0] else z2[z2 < z1[0]]

    # Top part of the atmosphere (no grid intersection)
    z_top = z2[z2 > z1[-1]] if z1[-1] < z2[-1] else z1[z1 > z2[-1]]

    # Middle part of the atmosphere (grids intersect)
    z1_middle = z1
    z2_middle = z2
    if z_bottom.any() & z_top.any():
        z1_middle = z1[(z1 > z_bottom[-1]) & (z1 < z_top[0])]
        z2_middle = z2[(z2 > z_bottom[-1]) & (z2 < z_top[0])]
    elif z_top.any():
        z1_middle = z1[z1 < z_top[0]]
        z2_middle = z2[z2 < z_top[0]]
    elif z_bottom.any():
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

    return np.concatenate((z_bottom, z_middle, z_top))
