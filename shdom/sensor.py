import xarray as xr
import numpy as np
import itertools

def make_sensor_dataset(x, y, z, mu, phi, stokes, wavelength):
    """
    TODO
    """
    stokes_bool = [component in stokes for component in ['I', 'Q', 'U', 'V']]
    dataset = xr.Dataset(
            data_vars={
                'wavelength': wavelength,
                'stokes': ('stokes_index', np.array(stokes_bool)),
                'cam_x': (['total_pixels'], x.ravel().astype(np.float64)),
                'cam_y': (['total_pixels'], y.ravel().astype(np.float64)),
                'cam_z': (['total_pixels'], z.ravel().astype(np.float64)),
                'cam_mu'  : (['total_pixels'], mu.ravel()),
                'cam_phi' : (['total_pixels'], phi.ravel()),
                },
        coords = {'stokes_index': np.array(['I', 'Q', 'U', 'V'])}
    )

    #Various other information could be added
    #about the grouping and relation of pixels to one another,
    #but this is the bare minimum for rendering.
    return dataset

def merge_sensor_list(list_of_sensor_datasets):
    """
#merging sensors
#concatenates, camx/y/z/mu/phi/superpix_i/superpix_wt along total_pixels dimension.
#concatenates npixels/observable_list along nimage dimension.
    """
    var_list = ['cam_x','cam_y','cam_z','cam_mu','cam_phi','super_pixel_index','super_pixel_weight']

    #make sure super_pixel_indices are unique.
    for i,sensor in enumerate(list_of_sensor_datasets):
        if i>=1:
            sensor['super_pixel_index'] += list_of_sensor_datasets[i-1]['super_pixel_index'].max() + 1

    concat_observable = xr.concat([data.stokes for data in list_of_sensor_datasets],dim='nimage')
    concat_image_shape = xr.concat([data.image_shape for data in list_of_sensor_datasets],dim='nimage')

    concatenated_pixels = xr.concat(list_of_sensor_datasets,data_vars=var_list,dim='total_pixels')
    merge_list = [concatenated_pixels.data_vars[name] for name in var_list]
    merge_list += [concat_observable,concat_image_shape]
    merged = xr.merge(merge_list)
    return merged

def split_sensors(combined_render_output):
    """
    TODO
    This function both splits the sensors, does the averaging of rays over pixels
    and subsets which observables are specified by the original sensor.
    These latter functions could be done at the sensor level in a 'make measurements' step.
    """

    #averaging over rays in each super_pixel
    averaged_over_super_pixels = (combined_render_output['super_pixel_weight']*combined_render_output).groupby('super_pixel_index').mean()
    #drop the super_pixel_dimension when irrelevant broadcasting occurs.
    #TODO tidy this up. only apply to necessary variables.
    averaged_over_super_pixels['image_shape'] = averaged_over_super_pixels.image_shape[0].astype(np.int64)
    averaged_over_super_pixels['observable_list'] = averaged_over_super_pixels.observable_list[0].astype(np.int64)

    list_of_unmerged = []
    count = 0
    for i in range(averaged_over_super_pixels.sizes['nimage']):
        split_index = averaged_over_super_pixels.image_shape[i].prod()
        split = averaged_over_super_pixels.sel({'super_pixel_index': slice(count,count+split_index),
                                               'nimage':i})

        #only_take the I,Q,U,V designated as observables
        for i,observable in enumerate(split.coords['observables'].data):
            if (split.observable_list[i]==0) and (observable in split):
                split = split.drop_vars(observable)

        list_of_unmerged.append(split)
        count+=split_index
    return list_of_unmerged

def add_subpixel_rays(sensor,IFOV, n_rays=1, weights='gaussian', sampling='deterministic'):
    """
    TODO
    """
    if IFOV == 0.0:
        sensor['super_pixel_index'] = ('total_pixels', range(sensor.sizes['total_pixels']))
        sensor['super_pixel_weight'] = ('total_pixels', np.ones(sensor.sizes['total_pixels']))
        return sensor
    else:
        raise NotImplementedError


def _homography_projection(projection_matrix, point_array):
    """
    Project 3D coordinates according to a projection matrix

    Parameters
    ----------
    projection matrix: np.array(shape=(3,4), dtype=float)
        The projection matrix.
    point_array: np.array(shape=(3, num_points), dtype=float)
        An array of num_points 3D points (x,y,z) [km]
    """
    homogenic_point_array = np.pad(point_array, ((0, 1), (0, 0)), 'constant', constant_values=1)
    return np.dot(projection_matrix, homogenic_point_array)


def orthographic_projection(wavelength, bounding_box, x_resolution, y_resolution, azimuth, zenith,
                            altitude='TOA', stokes='I'):
    """
    A parallel ray projection.

    Parameters
    ----------
    wavelength: float,
        Wavelength in [micron]
    bounding_box: xarray.Dataset/DataArray
        An xarray which contains the x,y,z coordinates which will be used
        to compute the bounding_box for the domain.
    x_resolution: float
        Pixel resolution [km] in x axis (North)
    y_resolution: float
        Pixel resolution [km] in y axis (East)
    azimuth: float
        Azimuth angle [deg] of the measurements (direction of photons)
    zenith: float
        Zenith angle [deg] of the measurements (direction of photons)
    altitude: float or 'TOA' (default)
       1. 'TOA': Top of the atmosphere.
       2. float: Altitude of the  measurements.
    stokes: list or string
       list or string of stokes components to observe ['I', 'Q', 'U', 'V']
    """
    mu = np.cos(np.deg2rad(zenith))
    phi = np.deg2rad(azimuth)
    xmin, ymin, zmin = bounding_box.x.data.min(),bounding_box.y.data.min(),bounding_box.z.data.min()
    xmax, ymax, zmax = bounding_box.x.data.max(),bounding_box.y.data.max(),bounding_box.z.data.max()

    altitude = zmax if altitude == 'TOA' else altitude

    # Project the bounding box onto the image plane
    alpha = np.sqrt(1 - mu ** 2) * np.cos(phi) / mu
    beta = np.sqrt(1 - mu ** 2) * np.sin(phi) / mu
    projection_matrix = np.array([
        [1, 0, -alpha, alpha * altitude],
        [0, 1, -beta, beta * altitude],
        [0, 0, 0, altitude]
    ])
    bounding_box_8point_array = np.array(list(itertools.product([xmin, xmax], [ymin, ymax], [zmin, zmax]))).T
    projected_bounding_box = _homography_projection(projection_matrix, bounding_box_8point_array)

    # Use projected bounding box to define image sampling
    x_s, y_s = projected_bounding_box[:2, :].min(axis=1)
    x_e, y_e = projected_bounding_box[:2, :].max(axis=1)
    x = np.arange(x_s, x_e + 1e-6, x_resolution)
    y = np.arange(y_s, y_e + 1e-6, y_resolution)
    z = altitude
    image_shape = [x.size, y.size]

    x, y, z, mu, phi = np.meshgrid(x, y, z, mu, phi)
    sensor = make_sensor_dataset(x.ravel(), y.ravel(), z.ravel(), mu.ravel(), phi.ravel(), stokes, wavelength)
    sensor['bounding_box'] = xr.DataArray(np.array([xmin,ymin,zmin,xmax,ymax,zmax]),
    coords={'bbox': ['xmin','ymin','zmin','xmax','ymax','zmax']},dims='bbox')
    sensor['image_shape'] = xr.DataArray(image_shape, coords={'image_dims': ['nx', 'ny']}, dims='image_dims')
    sensor.attrs = {
        'projection': 'Orthographic',
        'altitude': altitude,
        'x_resolution': x_resolution,
        'y_resolution': y_resolution,
        'projection_azimuth': azimuth,
        'projection_zenith': zenith
    }
    return sensor
