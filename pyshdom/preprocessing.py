"""
Data pre-processing tools
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import netCDF4 as nc
from scipy import optimize

def project_to_altitude(theta, phi, altitude):
    """
    Get horizontal grid offsets for projection of a ray to a specified altitude.

    Parameters
    ----------
    theta : float or np.array,
        Spherical coordinates zenith angle in [deg]
    phi : float or np.array,
        Spherical coordinates azimuth angle in [deg]
    altitude : float or np.array,
        Altitude (e.g. top of the atmosphere)

    Returns
    -------
    delta_x: float or np.array,
        offsets in x-direction
    delta_y: float or np.array,
        offsets in y-direction
    """
    tan_theta = np.tan(np.deg2rad(theta)) * altitude
    delta_x = tan_theta * np.cos(np.deg2rad(phi))
    delta_y = tan_theta * np.sin(np.deg2rad(phi))
    return delta_x, delta_y

def spherical_coords_to_vector(theta, phi):
    """
    Transform spherical coordinate angles into a direction vector

    Parameters
    ----------
    theta : float,
        Spherical coordinates zenith angle in [deg]
    phi : float,
        Spherical coordinates azimuth angle in [deg]

    Returns
    -------
    vector: np.array,
        3D direction vector
    """
    theta = np.deg2rad(np.mod(theta, 180.0))
    phi = np.deg2rad(phi)
    vector = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)]).T
    return vector.squeeze()

def find_intersection_line(rays, reference_view='000N', x_bounds=(-np.inf, np.inf),
                           y_bounds=(-np.inf, np.inf), z_bounds=(-np.inf, np.inf),
                           vx_bounds=(-np.inf, np.inf), vy_bounds=(-np.inf, np.inf),
                           vz_bounds=(-np.inf, np.inf)):
    """Find the intersection of rays with multiple views with a line of constant velocity.


    Parameters
    ----------
    rays: dict,
        Dictionary with multiple viewing angles each containing following fields:
            'x': float,
            'y': float,
            'z': float,
            'direction': np.array(shape=(3), dtype=float)
            'time': numpy.datetime64
    reference_view: string, default='000N'
        Reference view should be  key of the rays dictionary.
    x_bounds: (float, float), default=(-np.inf, np.inf),
        Limit the optimization search space for 3D COM output coordinate 'x'
    y_bounds: (float, float), default=(-np.inf, np.inf),
        Limit the optimization search space for 3D COM output coordinate 'y'
    z_bounds: (float, float), default=(-np.inf, np.inf),
        Limit the optimization search space for 3D COM output coordinate 'z'
    vx_bounds: (float, float), default=(-np.inf, np.inf),
        Limit the optimization search space for 3D COM output coordinate 'vx'
    vy_bounds: (float, float), default=(-np.inf, np.inf),
        Limit the optimization search space for 3D COM output coordinate 'vy'
    vz_bounds: (float, float), default=(-np.inf, np.inf),
        Limit the optimization search space for 3D COM output coordinate 'vz'

    Returns
    -------
    com: np.array(shape=(3), dtype=float)
        Center of mass in 3D space.
    velocity: np.array(shape=(3), dtype=float)
        Velocity vectors (in grid_units/second) in 3D space.

    Notes
    -----
    Bounds should be in units of the grid.
    """
    ndim = rays[reference_view]['direction'].size
    R = np.zeros((ndim * 2, ndim * 2))
    q = np.zeros(ndim * 2)
    for ray in rays.values():
        direction = ray['direction']
        com_coords = np.array([ray['x'], ray['y'], ray['z']])
        dt = np.timedelta64(ray['time'] - rays[reference_view]['time'], 's').astype('d')

        A = np.concatenate((np.eye(ndim), dt * np.eye(ndim)), axis=1)
        InnT = np.eye(ndim) - np.outer(direction, direction)
        R += np.matmul(A.T, np.matmul(InnT, A))
        q += np.matmul(A.T, np.matmul(InnT, com_coords))

    bounds = list(zip(x_bounds, y_bounds, z_bounds, vx_bounds, vy_bounds, vz_bounds))
    x = optimize.lsq_linear(R, q, bounds=bounds).x
    com, velocity = x[:ndim], x[ndim:]
    return com, velocity

def get_com_indices(image, thresh):
    """
    Compute the center of mass indices for the largest connected component object.

    Parameters
    ----------
    image : xr.DataArray or np.array,
        A 2D image
    thresh : float,
        A threshold used to make a binary `cloud mask`

    Returns
    -------
    comy, comx: int, int
        COM image coordinates in x,y axes.
    """
    mask_ = image > thresh
    mask = ndimage.binary_opening(ndimage.binary_closing(mask_))

    # label objects
    lbl, nlbl = ndimage.label(mask.astype('int'))

    # get largest object
    object_sizes = ndimage.sum(mask, lbl, range(nlbl + 1))
    max_lbl = np.argmax(object_sizes)

    # compute COMs in pixel coordinates
    coms = ndimage.measurements.center_of_mass(mask.astype('int'), lbl, np.unique(lbl))

    # COM for largest object
    comy, comx = coms[max_lbl]
    return int(comy), int(comx)

def get_roi(image, figsize=(10, 10)):
    """
    Interactively mark a Region of Interest (RIO) on a 2D image.
    Polygon points are selected using left click. Right click completes the selection and closes the image.
    See [1] for more detail.

    Parameters
    ----------
    image : xr.DataArray or np.array
        A 2D image
    figsize : tuple,
        Figure size.

    Returns
    -------
    roi: roipoly.RoiPoly
        Region of interest polygon.

    Notes
    -----
    Requires installing the package roipoly [1]. Use `pip install roipoly` if not installed.

    References
    ----------
    [1] https://github.com/jdoepfert/roipoly.py
    """
    try:
        from roipoly import RoiPoly
    except ModuleNotFoundError as e:
        print(e)

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.imshow(image)
    roi = RoiPoly(fig=fig, ax=ax, color='r', close_fig=True)
    return roi

def load_airmspi_data(path, bands=['355nm', '380nm', '445nm', '470nm', '555nm', '660nm', '865nm', '935nm']):
    """
    Utility function to load AirMSPI hdf files into xr.Dataset

    Parameters
    ----------
    path : str,
        path to file
    bands : list or string,
        AirMSPI bands to load. Default is all bands.

    Returns
    -------
    output : xr.Dataset or tuple of datasets
        A dataset with radiance or polarized fields or a tuple of datasets with both radiance and polarization fields.

    Notes
    -----
    AirMSPI is an 8-band (355, 380, 445, 470, 555, 660, 865, 935 nm) pushbroom camera, measuring polarization in the
    470, 660, and 865 nm bands, mounted on a gimbal to acquire multiangular observations over a ±67° along-track range.
    Two principal observing modes are employed: step-and-stare, in which 11 km x 11 km targets are observed at a
    discrete set of view angles with a spatial resolution of ~10 m; and continuous sweep, in which the camera
    slews back and forth along the flight track between ±67° to acquire wide area coverage
    (11 km swath at nadir, target length 108 km) with ~25 m spatial resolution. See references for more information.

    References
    ----------
    https://asdc.larc.nasa.gov/project/AIRMSPI

    Raises
    ------
    AttributeError if band not supported.
    Supported types are: '355nm', '380nm', '445nm', '470nm', '555nm', '660nm', '865nm', '935nm'
    """
    ncf = nc.Dataset(path, diskless=True, persist=False)
    pol_bands = ['470nm', '660nm', '865nm']
    rad_bands = ['355nm', '380nm', '445nm', '555nm', '935nm']
    bands = np.atleast_1d(bands)

    for band in bands:
        if (band not in pol_bands) and (band not in rad_bands):
            raise AttributeError('Band: {} not recognized'.format(band))

    dims, rad_dataset, pol_dataset = {}, [], []
    for band in bands:
        data = []
        data_fields = ncf['HDFEOS/GRIDS/{}_band/Data Fields/'.format(band)]
        for name, var in data_fields.variables.items():
            if name in dims.keys():
                continue
            attrs_keys = [attr for attr in var.ncattrs() if not attr.startswith('_')]
            attrs_vals = [var.getncattr(attr) for attr in attrs_keys]
            attrs = dict(zip(attrs_keys, attrs_vals))

            use_dims = True
            for dim in var.dimensions:
                if dim in data_fields.variables.keys():
                    if dim not in dims:
                        dim_attrs_keys = [attr for attr in data_fields[dim].ncattrs() if not attr.startswith('_')]
                        dim_attrs_vals = [data_fields[dim].getncattr(attr) for attr in dim_attrs_keys]
                        dim_attrs = dict(zip(dim_attrs_keys, dim_attrs_vals))
                        dims[dim] = xr.DataArray(data_fields[dim], dims=dim, attrs=dim_attrs)
                else:
                    use_dims = False

            if use_dims:
                current_dims = {k: dims[k] for k in var.dimensions}
                data_field = xr.DataArray(var[:], attrs=attrs, name=name, dims=list(current_dims.keys()),
                                          coords=current_dims)
                # Longitude and Latitude
                lon = xr.DataArray(ncf['HDFEOS/GRIDS/Ancillary/Data Fields/Longitude'], dims=['YDim', 'XDim'],
                                   coords=current_dims)
                lat = xr.DataArray(ncf['HDFEOS/GRIDS/Ancillary/Data Fields/Latitude'], dims=['YDim', 'XDim'],
                                   coords=current_dims)
                data_field = data_field.assign_coords({'Latitude': lat, 'Longitude': lon})
            else:
                data_field = xr.DataArray(var[:], attrs=attrs, name=name)

            data.append(data_field.expand_dims(band=[band]))

        data = xr.merge(data).squeeze()
        if band in pol_bands:
            pol_dataset.append(data)
        else:
            rad_dataset.append(data)

    output = []
    if len(rad_dataset) > 0:
        rad_dataset = xr.concat(rad_dataset, dim='band')
        output.append(rad_dataset)
    if len(pol_dataset) > 0:
        pol_dataset = xr.concat(pol_dataset, dim='band')
        output.append(pol_dataset)

    output = output[0] if len(output) == 1 else output
    return output