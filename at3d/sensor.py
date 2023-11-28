"""This module contains functions for creating 'sensor' datasets that contain
the geometric information necessary to render/calculate observables
(e.g. Stokes Vector) given an RTE solution from a solver.RTE object.

A generic method, `make_sensor_dataset` is first defined as well as two specific
geometric projections: `orthographic_projection` and `perspective_projection`.

A sensor dataset contains both 'pixel' and 'ray' variables.
Each pixel has a single pointing direction and position which define the line of
sight. For each pixel, there is a single 'observable' quantity; the Stokes Vector
at a specified wavelength. A pixel may also have a non-zero field of view.
As such, the field of view of each pixel is described by a set of ray variables.
Each pixel may have several rays, each of which has its own pointing direction,
position and weight. The observables are calculated at the specified ray geometries
and then a weighted sum is performed over the rays to calculate the pixel level
observable quantity.
As the rays are the variables that are actually passed to the solver.RTE object for
rendering, a sensor MUST have ray variables. However, there is no generic method
for generating sub-pixel ray geometry, as it depends on the assumed sensor geometry.
Users should add their own generating functions for specialized sensors.
"""
import itertools
import inspect
import sys
import xarray as xr
import numpy as np
from collections import OrderedDict

import at3d.checks

def make_sensor_dataset(x, y, z, mu, phi, stokes, wavelength, fill_ray_variables=False):
    """
    A generic method to generate an xr.Dataset which specifies the measurement
    geometry of a sensor.

    The position (`x`, `y`, `z`) and angle (`mu`, `phi`) of each pixel in the sensor is
    specified along with the components of the Stokes Vector to be simulated at each pixel
    and the monochromatic wavelength of the measurements. `x`, `y`, `z`, `mu`, `phi`
    should all be 1D and have the same length.

    Parameters
    ----------
    x, y, z : array_like of floats
        The spatial positions of the sensor's pixels.
    mu : array_like of floats
        The cosine of zenith angle pointing TOWARDS the sensor's pixels.
    phi : array_like of floats
        The azimuthal angle angle (in radians) pointing TOWARDS the sensor's pixels.
    stokes : string or list of strings
        The Stokes components that will be calculated for the sensor at each pixel.
        Valid values ('I', 'Q', 'U', 'V')
    wavelength : float
        monochromatic wavelength [microns] that the `stokes` observables will be calculated at.
    fill_ray_variables : {'False', 'True'}, optional
        If True, then 'ray' variables are created corresponding to the specified
        pixel variables. If this function is being used in isolation without a geometric projection
        (e.g. `orthographic_projection`, `perspective_projection`) to generate subpixel
        rays then `fill_ray_variables` should be set to 'True'.

    Returns
    -------
    dataset : xr.Dataset
        A dataset containing the sensor pixel geometries, wavelength and
        Stokes observables. May also contain ray variables.

    Raises
    ------
    ValueError
        If `x`, `y`, `z`,`mu`, `phi` are not the same shape (and 1D) or if
        any `mu` = 0.0 (directly horizontal radiances are not allowed by SHDOM) or
        if the specified pixel altitudes `z` are less than or equal to zero.
        Or if the stokes components are not valid.
    """
    x, y, z, mu, phi, wavelength, stokes = np.asarray(x), np.asarray(y), np.asarray(z), \
                                           np.asarray(mu), np.asarray(phi), np.asarray(wavelength),\
                                           np.atleast_1d(stokes)
    for totest, name in zip((x, y, z, mu, phi), ('x', 'y', 'z', 'mu', 'phi')):
        if not totest.ndim == 1:
            raise ValueError("'{}' should be a 1-D np.ndarray".format(name))
    if not all([x.size == i.size for i in (y, z, mu, phi)]):
        raise ValueError("All of x, y, z, mu, phi should have the same size.")

    if not np.all(z >= 0.0):
        raise ValueError("All altitudes (z) must be >= 0.0")
    if np.any(mu == 0.0):
        raise ValueError("'mu' values of 0.0 are not allowed.")

    for i in stokes:
        if i not in ('I', 'Q', 'U', 'V'):
            raise ValueError("Valid Stokes components are 'I', 'Q', 'U', 'V' not '{}'".format(i))

    stokes_bool = [component in stokes for component in ('I', 'Q', 'U', 'V')]
    dataset = xr.Dataset(
        data_vars={
            'wavelength': wavelength,
            'stokes': ('stokes_index', np.array(stokes_bool)),
            'cam_x': (['npixels'], x.ravel().astype(np.float64)),
            'cam_y': (['npixels'], y.ravel().astype(np.float64)),
            'cam_z': (['npixels'], z.ravel().astype(np.float64)),
            'cam_mu'  : (['npixels'], mu.ravel()),
            'cam_phi' : (['npixels'], phi.ravel()),
            },
        coords={'stokes_index': np.array(['I', 'Q', 'U', 'V'])}
    )
    if fill_ray_variables:
        dataset = _add_null_subpixel_rays(dataset)
    return dataset

def stochastic(npixels, nrays, seed=None):
    """
    Generates random position perturbations and weights for sub-pixel rays.

    Random numbers are drawn from a uniform distribution in the range (-1.0,1.0]
    with shape=(`npixels`, `nrays`).

    Parameters
    ----------
    npixels : int
        The number of pixels for which sub-pixel rays should be generated.
    nrays : int
        The number of sub-pixel rays to generate per pixel.
    seed : {None, int}
        The random seed to use to generate

    Returns
    -------
    position_perturbations : np.ndarray
        Has shape=(`npixels`, `nrays`). perturbations to the position of a subpixel-ray
        in the IMAGE PLANE.
    weights : np.ndarray
        Has shape=(`npixels`, `nrays`) with a uniform value of 1.0/`nrays` so that they
        sum to unity.

    Notes
    -----
    These random numbers are designed to be position perturbations in pixel
    coordinates in the IMAGE PLANE of a sensor though this function could be
    used for other purposes.
    The weights are normalized to sum to unity.
    """
    if seed is not None:
        np.random.seed(seed)
    position_perturbations = np.random.uniform(low=-1.0, high=1.0, size=(npixels, nrays))
    weights = np.ones((npixels, nrays))/nrays
    return position_perturbations, weights

def uniform(npixels, nrays):
    position_perturbations = np.linspace(-1.0+1/nrays,1.0-1/nrays,nrays)
    weights = np.ones((npixels, nrays))/nrays
    return position_perturbations, weights

def gaussian(npixels, degree):
    """
    Generates gauss-legendre weights and positions for integration of sub-pixel
    rays.

    Parameters
    ----------
    npixels : int
        The number of pixels for which sub-pixel rays should be generated.
    degree : int
        The degree of the gauss-legendre quadrature.

    Returns
    -------
    position_perturbations : np.ndarray
        shape=(`npixels`, `degree`). perturbations to the position of a subpixel-ray
        in the IMAGE PLANE.
    weights : np.ndarray
        shape=(`npixels`, `degree`). The corresponding weights of each subpixel-ray,
        normalized to sum to unity.

    See Also
    --------
    np.polynomial.legendre.leggauss

    Notes
    -----
    This integration scheme treats every pixel independently so the integration
    of radiance in a pixel in the IMAGE PLANE does not utilize information from the
    surrounding pixels to interpolate the radiance at the edge of the pixel.
    """
    out_mu, out_mu_weights = np.polynomial.legendre.leggauss(degree)
    position_perturbations = np.repeat(out_mu[np.newaxis, :], npixels, axis=0)
    weights = np.repeat(out_mu_weights[np.newaxis, :], npixels, axis=0)/np.sum(out_mu_weights)
    return position_perturbations, weights

def orthographic_projection(wavelength, bounding_box, x_resolution, y_resolution, azimuth, zenith,
                            altitude='TOA', stokes='I', sub_pixel_ray_args={'method': None}):
    """
    Generates a sensor dataset which views a bounding box with an orthographic projection at
    a given orientation and resolution.

    This function will use the geometric information from the `bounding_box` and viewing orientation
    (`azimuth`, `zenith`) and pixel spacing (`x_resolution`, `y_resolution`) to generate the positions
    and angles of individual pixels in the sensor. In addition, as the geometric projection
    is well defined and simple, sub-pixel rays can be generated by passing a method to the
    `sub_pixel_ray_args`. In addition to viewing geometry, a sensor dataset has some
    observables. These are defined by a monochromatic `wavelength` and a list of
    desired stokes components `stokes`.

    Parameters
    ----------
    wavelength : float
        The monochromatic wavelength that the sensor will observe at [micron].
    bounding_box : xr.Dataset/xr.DataArray
        A grid object which contains 'x', 'y', 'z' coordinates whose maximum and minimum
        will set the span of the domain that is observed by the sensor.
    x_resolution : float
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
    sub_pixel_ray_args : dict
        dictionary defining the method for generating sub-pixel rays. The callable
        which generates the position_perturbations and weights (e.g. at3d.sensor.gaussian)
        should be set as the 'method', while arguments to that callable, should be set as
        other entries in the dict. Each argument have two values, one for each of the
        x and y axes of the image plane, respectively.
        E.g. sub_pixel_ray_args={'method':at3d.sensor.gaussian, 'degree': (2, 3)}

    Returns
    -------
    sensor : xr.Dataset
        A dataset containing all of the information required to define a sensor
        for which synthetic measurements can be simulated;
        positions and angles of all pixels, sub-pixel rays and their associated weights,
        and the sensor's observables.

    See Also
    --------
    at3d.sensor.make_sensor_dataset

    Notes
    -----
    An orthographic projection has parallel rays, so each pixel will view the bounding box
    at the same angle. The `x_resolution` and `y_resolution` is the spacing of the
    observed rays in the horizontal coordinates of the `bounding box`. This means that
    as the viewing angle of the bounding box becomes more oblique, the orthogonal distance
    between rays decreases.
    """
    mu = np.cos(np.deg2rad(zenith))
    phi = np.deg2rad(azimuth)
    at3d.checks.check_grid(bounding_box)
    xmin, ymin, zmin = bounding_box.x.data.min(), bounding_box.y.data.min(), bounding_box.z.data.min()
    xmax, ymax, zmax = bounding_box.x.data.max(), bounding_box.y.data.max(), bounding_box.z.data.max()

    altitude = zmax if altitude == 'TOA' else altitude

    # Project the bounding box onto the image plane
    alpha = np.sqrt(1 - mu ** 2) * np.cos(phi) / mu
    beta = np.sqrt(1 - mu ** 2) * np.sin(phi) / mu
    projection_matrix = np.array([
        [1, 0, -alpha, alpha * altitude],
        [0, 1, -beta, beta * altitude],
        [0, 0, 0, altitude]
    ])
    bounding_box_8point_array = np.array(list(itertools.product([xmin, xmax],
                                                                [ymin, ymax],
                                                                [zmin, zmax]))).T
    projected_bounding_box = _homography_projection(projection_matrix, bounding_box_8point_array)

    # Use projected bounding box to define image sampling
    x_s, y_s = projected_bounding_box[:2, :].min(axis=1)
    x_e, y_e = projected_bounding_box[:2, :].max(axis=1)

    # Always pad this so that there is a point at the domain edge
    # even if it is past it.
    # Note that this introduces an asymmetry between a forward and backward
    # sampling sensor i.e. between a nadir sensor with azimuths of 0 and 180.0
    # even though they should be identical.
    x = np.arange(x_s, x_e + x_resolution, x_resolution)
    y = np.arange(y_s, y_e + y_resolution, y_resolution)

    z = altitude
    image_shape = [x.size, y.size]

    x, y, z, mu, phi = np.meshgrid(x, y, z, mu, phi)
    sensor = make_sensor_dataset(x.ravel(), y.ravel(), z.ravel(), mu.ravel(), phi.ravel(), stokes, wavelength)
    sensor['bounding_box'] = xr.DataArray(np.array([xmin, ymin, zmin, xmax, ymax, zmax]),
    coords={'bbox': ['xmin', 'ymin', 'zmin', 'xmax', 'ymax', 'zmax']}, dims='bbox')
    sensor['image_shape'] = xr.DataArray(image_shape, coords={'image_dims': ['nx', 'ny']}, dims='image_dims')
    sensor.attrs = {
        'projection': 'Orthographic',
        'altitude': altitude,
        'x_resolution': x_resolution,
        'y_resolution': y_resolution,
        'projection_azimuth': azimuth,
        'projection_zenith': zenith
    }
    if sub_pixel_ray_args['method'] is not None:

        #generate the weights and perturbations to the pixel positions in the image plane.
        sub_pixel_ray_method, subpixel_ray_kwargs_x, subpixel_ray_kwargs_y = \
            _parse_sub_pixel_ray_args(sub_pixel_ray_args)
        position_perturbations_x, weights_x = sub_pixel_ray_method(x.size, **subpixel_ray_kwargs_x)
        position_perturbations_y, weights_y = sub_pixel_ray_method(y.size, **subpixel_ray_kwargs_y)

        #merge the two dimensions
        perturbations_x = np.repeat(position_perturbations_x[..., np.newaxis]*x_resolution/2.0,
                                    position_perturbations_y.shape[-1], axis=-1)
        perturbations_y = np.repeat(position_perturbations_y[..., np.newaxis, :]*y_resolution/2.0,
                                    position_perturbations_x.shape[-1], axis=-2)
        big_weightx = np.repeat(weights_x[..., np.newaxis], weights_y.shape[-1], axis=-1)
        big_weighty = np.repeat(weights_y[..., np.newaxis, :], weights_x.shape[-1], axis=-2)

        #apply perturbations to original image plane coordinates.
        x_ray = (x.ravel()[:, np.newaxis, np.newaxis] + perturbations_x).ravel()
        y_ray = (y.ravel()[:, np.newaxis, np.newaxis] + perturbations_y).ravel()
        z_ray = np.repeat(np.repeat(z.ravel()[:, np.newaxis, np.newaxis],
                                    perturbations_x.shape[-2], axis=-2),
                          perturbations_y.shape[-1], axis=-1).ravel()
        mu_ray = np.repeat(np.repeat(mu.ravel()[:, np.newaxis, np.newaxis],
                                     perturbations_x.shape[-2], axis=-2),
                           perturbations_y.shape[-1], axis=-1).ravel()
        phi_ray = np.repeat(np.repeat(phi.ravel()[:, np.newaxis, np.newaxis],
                                      perturbations_x.shape[-2], axis=-2),
                            perturbations_y.shape[-1], axis=-1).ravel()
        #make the pixel indices and ray weights.
        pixel_index = np.repeat(np.repeat(range(len(sensor.cam_mu.data)),
                                          weights_x.shape[-1]), weights_y.shape[-1])
        ray_weight = (big_weightx*big_weighty).ravel()
        #update ray variables to sensor dataset.
        sensor['ray_mu'] = ('nrays', mu_ray)
        sensor['ray_phi'] = ('nrays', phi_ray)
        sensor['ray_x'] = ('nrays', x_ray)
        sensor['ray_y'] = ('nrays', y_ray)
        sensor['ray_z'] = ('nrays', z_ray)
        sensor['pixel_index'] = ('nrays', pixel_index)
        sensor['ray_weight'] = ('nrays', ray_weight)
        sensor['use_subpixel_rays'] = True
        sub_pixel_ray_args['method'] = sub_pixel_ray_args['method'].__name__
        for attribute in sub_pixel_ray_args:
            sensor.attrs['sub_pixel_ray_args_{}'.format(attribute)] = sub_pixel_ray_args[attribute]

    else:
        #duplicate ray variables to sensor dataset.
        sensor = _add_null_subpixel_rays(sensor)
    return sensor

def _homography_projection(projection_matrix, point_array):
    """
    Project 3D coordinates according to a projection matrix

    Parameters
    ----------
    projection matrix: np.array(shape=(3,4), dtype=float)
        The projection matrix.
    point_array: np.array(shape=(3, num_points), dtype=float)
        An array of num_points 3D points (x,y,z) [km]

    Returns
    -------
    The projection of the 3D points.
    """
    homogenic_point_array = np.pad(point_array, ((0, 1), (0, 0)), 'constant', constant_values=1)
    return np.dot(projection_matrix, homogenic_point_array)

def perspective_projection(wavelength, fov, x_resolution, y_resolution,
                           position_vector, lookat_vector, up_vector,
                           stokes='I', sub_pixel_ray_args={'method':None}):
    """
    Generates a sensor dataset that observes a target location with
    a perspective (pinhole camera) projection.

    Parameters
    ----------
    wavelength: float,
        Wavelength in [micron]
    fov: float
        Field of view [deg]
    x_resolution: int
        Number of pixels in camera x axis
    y_resolution: int
        Number of pixels in camera y axis
    position_vector: list of 3 float elements
        [x , y , z] which are:
        Location in global x coordinates [km] (North)
        Location in global y coordinates [km] (East)
        Location in global z coordinates [km] (Up)
    lookat_vector: list of 3 float elements
        [x , y , z] which are:
        Point in global x coordinates [km] (North) where the camera is pointing at
        Point in global y coordinates [km] (East) where the camera is pointing at
        Point in global z coordinates [km] (Up) where the camera is pointing at
    up_vector: list of 3 float elements
        The up vector determines the roll of the camera.
    stokes: list or string
       list or string of stokes components to observe ['I', 'Q', 'U', 'V'].
    sub_pixel_ray_args : dict
        dictionary defining the method for generating sub-pixel rays. The callable
        which generates the position_perturbations and weights (e.g. at3d.sensor.gaussian)
        should be set as the 'method', while arguments to that callable, should be set as
        other entries in the dict. Each argument have two values, one for each of the
        x and y axes of the image plane, respectively.
        E.g. sub_pixel_ray_args={'method':at3d.sensor.gaussian, 'degree': (2, 3)}

    Returns
    -------
    sensor : xr.Dataset
        A dataset containing all of the information required to define a sensor
        for which synthetic measurements can be simulated;
        positions and angles of all pixels, sub-pixel rays and their associated weights,
        and the sensor's observables.

    """
    norm = lambda x: x / np.linalg.norm(x, axis=0)

    #assert samples>=1, "Sample per pixel is an integer >= 1"
    #assert int(samples) == samples, "Sample per pixel is an integer >= 1"

    assert int(x_resolution) == x_resolution, "x_resolution is an integer >= 1"
    assert int(y_resolution) == y_resolution, "y_resolution is an integer >= 1"

    # The bounding_box is not nessesary in the prespactive projection, but we still may consider
    # to use if we project the rays on the bounding_box when the differences in mu , phi angles are below certaine precision.
    #     if(bounding_box is not None):

    #         xmin, ymin, zmin = bounding_box.x.data.min(),bounding_box.y.data.min(),bounding_box.z.data.min()
    #         xmax, ymax, zmax = bounding_box.x.data.max(),bounding_box.y.data.max(),bounding_box.z.data.max()

    nx = x_resolution
    ny = y_resolution
    position = np.array(position_vector, dtype=np.float32)
    lookat = np.array(lookat_vector, dtype=np.float32)
    up = np.array(up_vector)
    direction = lookat - position

    zaxis = norm(direction)
    xaxis = norm(np.cross(up, zaxis))
    yaxis = np.cross(zaxis, xaxis)
    rotation_matrix = np.stack((xaxis, yaxis, zaxis), axis=1)

    M = max(nx, ny)
    npix = nx*ny
    R = np.array([nx, ny])/M # R will be used to scale the sensor meshgrid.
    dy = 2*R[1]/ny # pixel length in y direction in the normalized image plane.
    dx = 2*R[0]/nx # pixel length in x direction in the normalized image plane.
    x_s, y_s, z_s = np.meshgrid(np.linspace(-R[0]+dx/2, R[0]-dx/2, nx),
                                np.linspace(-R[1]+dy/2, R[1]-dy/2, ny), 1.0)

    # Here x_c, y_c, z_c coordinates on the image plane before transformation to the requaired observation angle
    focal = 1.0 / np.tan(np.deg2rad(fov) / 2.0) # focal (normalized) length when the sensor size is 2 e.g. r in [-1,1).
    fov_x = np.rad2deg(2*np.arctan(R[0]/focal))
    fov_y = np.rad2deg(2*np.arctan(R[1]/focal))

    k = np.array([[focal, 0, 0],
                  [0, focal, 0],
                  [0, 0, 1]], dtype=np.float32)
    inv_k = np.linalg.inv(k)

    homogeneous_coordinates = np.stack([x_s.ravel(), y_s.ravel(), z_s.ravel()])

    x_c, y_c, z_c = norm(np.matmul(
        rotation_matrix, np.matmul(inv_k, homogeneous_coordinates)))
    # Here x_c, y_c, z_c coordinates on the image plane after transformation to the requaired observation

    # x,y,z mu, phi in the global coordinates:
    mu = -z_c.astype(np.float64)
    phi = (np.arctan2(y_c, x_c) + np.pi).astype(np.float64)
    x = np.full(npix, position[0], dtype=np.float32)
    y = np.full(npix, position[1], dtype=np.float32)
    z = np.full(npix, position[2], dtype=np.float32)

    image_shape = [nx,ny]
    sensor = make_sensor_dataset(x.ravel(), y.ravel(), z.ravel(),
                                 mu.ravel(), phi.ravel(), stokes, wavelength)
    # compare to orthographic projection, prespective projection may not have bounding box.
    #     if(bounding_box is not None):
    #         sensor['bounding_box'] = xr.DataArray(np.array([xmin,ymin,zmin,xmax,ymax,zmax]),
    #                                               coords={'bbox': ['xmin','ymin','zmin','xmax','ymax','zmax']},dims='bbox')

    sensor['image_shape'] = xr.DataArray(image_shape,
                                         coords={'image_dims': ['nx', 'ny']},
                                         dims='image_dims')
    sensor.attrs = {
        'projection': 'Perspective',
        'fov_deg': fov,
        'fov_x_deg': fov_x,
        'fov_y_deg': fov_y,
        'x_resolution': x_resolution,
        'y_resolution': y_resolution,
        'position': position,
        'lookat': lookat,
        'rotation_matrix': rotation_matrix.ravel(),
        'sensor_to_camera_transform_matrix':k.ravel()

    }

    if sub_pixel_ray_args['method'] is not None:

        #generate the weights and perturbations to the pixel positions in the image plane.
        sub_pixel_ray_method, subpixel_ray_kwargs_x, subpixel_ray_kwargs_y =  \
                            _parse_sub_pixel_ray_args(sub_pixel_ray_args)
        position_perturbations_x, weights_x = sub_pixel_ray_method(x_s.size,
                                                                   **subpixel_ray_kwargs_x)
        position_perturbations_y, weights_y = sub_pixel_ray_method(y_s.size,
                                                                   **subpixel_ray_kwargs_y)

        #merge the two dimensions
        perturbations_x = np.repeat(position_perturbations_x[..., np.newaxis]*dx/2.0,
                                    position_perturbations_y.shape[-1], axis=-1)
        perturbations_y = np.repeat(position_perturbations_y[..., np.newaxis, :]*dy/2.0,
                                    position_perturbations_x.shape[-1], axis=-2)
        big_weightx = np.repeat(weights_x[..., np.newaxis], weights_y.shape[-1], axis=-1)
        big_weighty = np.repeat(weights_y[..., np.newaxis, :], weights_x.shape[-1], axis=-2)

        #apply perturbations to original image plane coordinates.
        x_ray = (x_s.ravel()[:, np.newaxis, np.newaxis] + perturbations_x).ravel()
        y_ray = (y_s.ravel()[:, np.newaxis, np.newaxis] + perturbations_y).ravel()
        z_ray = np.repeat(np.repeat(z_s.ravel()[:, np.newaxis, np.newaxis],
                                    perturbations_x.shape[-2], axis=-2),
                          perturbations_y.shape[-1], axis=-1).ravel()
        ray_homogeneous = np.stack([x_ray, y_ray, z_ray])

        x_c, y_c, z_c = norm(np.matmul(
            rotation_matrix, np.matmul(inv_k, ray_homogeneous)))
        # Here x_c, y_c, z_c coordinates on the image plane after transformation to the requaired observation

        # x,y,z mu, phi in the global coordinates:
        mu = -z_c.astype(np.float64)
        phi = (np.arctan2(y_c, x_c) + np.pi).astype(np.float64)
        x = np.full(x_c.size, position[0], dtype=np.float32)
        y = np.full(x_c.size, position[1], dtype=np.float32)
        z = np.full(x_c.size, position[2], dtype=np.float32)

        #make the pixel indices and ray weights.
        pixel_index = np.repeat(np.repeat(range(len(sensor.cam_mu.data)),
                                          weights_x.shape[-1]), weights_y.shape[-1])
        ray_weight = (big_weightx*big_weighty).ravel()
        #update ray variables to sensor dataset.
        sensor['ray_mu'] = ('nrays', mu)
        sensor['ray_phi'] = ('nrays', phi)
        sensor['ray_x'] = ('nrays', x)
        sensor['ray_y'] = ('nrays', y)
        sensor['ray_z'] = ('nrays', z)
        sensor['pixel_index'] = ('nrays', pixel_index)
        sensor['ray_weight'] = ('nrays', ray_weight)
        sensor['use_subpixel_rays'] = True

        sub_pixel_ray_args['method'] = sub_pixel_ray_args['method'].__name__
        for attribute in sub_pixel_ray_args:
            sensor.attrs['sub_pixel_ray_args_{}'.format(attribute)] = sub_pixel_ray_args[attribute]
    else:
            #duplicate ray variables to sensor dataset.
        sensor = _add_null_subpixel_rays(sensor)
    return sensor

def domaintop_projection(wavelength, bounding_box, x_resolution, y_resolution, azimuth, zenith,
               x_offset=0.0, y_offset=0.0, stokes='I',sub_pixel_ray_args={'method': None}):
    """
    This projection generator is like the default SHDOM radiance calculation method. A single uniform spacing
    of radiances is specified across the domain top as well as a global offset of the grid.
    This is useful for 1D simulations, where the registration of oblique angles to the domain top is important.
    In 3D the orthographic_projection calculates the correct offset for you so it is preferred.

    Parameters
    ----------
    wavelength : float
        The monochromatic wavelength that the sensor will observe at [micron].
    bounding_box : xr.Dataset/xr.DataArray
        A grid object which contains 'x', 'y', 'z' coordinates whose maximum and minimum
        will set the span of the domain that is observed by the sensor.
    x_resolution : float
        Pixel resolution [km] in x axis (North)
    y_resolution: float
        Pixel resolution [km] in y axis (East)
    azimuth: float
        Azimuth angle [deg] of the measurements (direction of photons)
    zenith: float
        Zenith angle [deg] of the measurements (direction of photons)
    stokes: list or string
       list or string of stokes components to observe ['I', 'Q', 'U', 'V']
    sub_pixel_ray_args : dict
        dictionary defining the method for generating sub-pixel rays. The callable
        which generates the position_perturbations and weights (e.g. at3d.sensor.gaussian)
        should be set as the 'method', while arguments to that callable, should be set as
        other entries in the dict. Each argument have two values, one for each of the
        x and y axes of the image plane, respectively.
        E.g. sub_pixel_ray_args={'method':at3d.sensor.gaussian, 'degree': (2, 3)}

    Returns
    -------
    sensor : xr.Dataset
        A dataset containing all of the information required to define a sensor
        for which synthetic measurements can be simulated;
        positions and angles of all pixels, sub-pixel rays and their associated weights,
        and the sensor's observables.

    See Also
    --------
    at3d.sensor.make_sensor_dataset
    at3d.sensor.orthographic_projection
    """
    sensor = at3d.sensor.orthographic_projection(wavelength, bounding_box, x_resolution, y_resolution, 0.0,
                                           0.0, altitude='TOA', stokes=stokes,
                                                   sub_pixel_ray_args=sub_pixel_ray_args)
    sensor.cam_x[:] += x_offset
    sensor.ray_x[:] += x_offset
    sensor.cam_y[:] += y_offset
    sensor.ray_y[:] += y_offset

    sensor.cam_mu[:] = np.cos(np.deg2rad(zenith))
    sensor.cam_phi[:] = np.deg2rad(azimuth)
    sensor.ray_mu[:] = np.cos(np.deg2rad(zenith))
    sensor.ray_phi[:] = np.deg2rad(azimuth)

    sensor.attrs['projection'] = 'DomainTop'
    sensor.attrs['projection_azimuth'] = azimuth
    sensor.attrs['projection_zenith'] = zenith

    return sensor


def _parse_sub_pixel_ray_args(sub_pixel_ray_args):
    """

    """
    try:
        subpixel_ray_params = inspect.signature(sub_pixel_ray_args['method'])
    except TypeError as err:
        raise TypeError("sub_pixel_ray_args 'method' must be a"
                        "callable object not '{}' of type '{}'".format(
                            sub_pixel_ray_args['method'],
                            type(sub_pixel_ray_args['method']))) from err
    subpixel_ray_kwargs_x = {}
    subpixel_ray_kwargs_y = {}
    subpixel_ray_params = inspect.signature(sub_pixel_ray_args['method']).parameters
    # Find only the kwargs that are in the signature of the subpixel ray method.
    # Bad kwargs raise a KeyError.
    for name, value in sub_pixel_ray_args.items():
        if name != 'method':
            try:
                if subpixel_ray_params[name].kind in (subpixel_ray_params[name].POSITIONAL_OR_KEYWORD,
                                                      subpixel_ray_params[name].KEYWORD_ONLY,
                                                      subpixel_ray_params[name].VAR_KEYWORD):
                    if isinstance(value, tuple):
                        subpixel_ray_kwargs_x[name] = value[0]
                        subpixel_ray_kwargs_y[name] = value[1]
                    else:
                        subpixel_ray_kwargs_x[name] = value
                        subpixel_ray_kwargs_y[name] = value
            except KeyError as err:
                raise type(err)(str(err).replace('"', "") + \
                    "Invalid kwarg '{}' passed to the "
                    "sub_pixel_ray_args['method'] callable. '{}'".format(
                        name, sub_pixel_ray_args['method'].__name__)
                    ).with_traceback(sys.exc_info()[2])

    return sub_pixel_ray_args['method'], subpixel_ray_kwargs_x, subpixel_ray_kwargs_y

def _add_null_subpixel_rays(sensor):
    """
    Fills in missing ray variables with copies of pixel (cam) variables in a sensor.

    Notes
    -----
    Ray variables are those that are used directly to calculate radiances
    so they must be present in a sensor for it to be valid. In the simplest case,
    there are no sub-pixel rays and the ray variables are just copies of the pixel
    variables.
    """
    for name in ('cam_mu', 'cam_phi', 'cam_x', 'cam_y', 'cam_z', 'cam_mu'):
        if name not in sensor.data_vars:
            raise ValueError("'{}' is missing from sensor. This is not a valid sensor.".format(name))
    sensor['ray_mu'] = ('nrays', sensor.cam_mu.data)
    sensor['ray_phi'] = ('nrays', sensor.cam_phi.data)
    sensor['ray_x'] = ('nrays', sensor.cam_x.data)
    sensor['ray_y'] = ('nrays', sensor.cam_y.data)
    sensor['ray_z'] = ('nrays', sensor.cam_z.data)
    sensor['pixel_index'] = ('nrays', range(len(sensor.cam_mu.data)))
    sensor['ray_weight'] = ('nrays', np.ones(len(sensor.cam_mu.data)))
    sensor['use_subpixel_rays'] = False
    return sensor

class BandModel:

    def __init__(self, identifier, wavelengths, weights):
        self._identifier = identifier
        self.weights = OrderedDict()
        for wavelength, weight in zip(np.atleast_1d(wavelengths), np.atleast_1d(weights)):
            self.weights[wavelength] = weight

    def weight(self, wavelength):
        return self.weights[wavelength]
    
    @property
    def wavelengths(self):
        return np.array(list(self.weights.keys()))

    @property
    def id(self):
        return self._identifier
        
class Monochromatic(BandModel):

    def __init__(self, wavelength):
        BandModel.__init__(self, wavelength, wavelength, 1.0)
        
class Satellite(BandModel):

    def __init__(self, instrument, band, parameterization='reptran'):

        if parameterization != 'reptran':
            raise ValueError(
                "Only the `reptran` satellite band model is currently supported"
            )
            
        if not instrument in ('terra-modis', 'aqua-modis'):
            raise ValueError(
                "Currently only supports Terra and Aqua MODIS not other instruments."
            )

        test_atmosphere = at3d.gas_absorption.load_standard_atmosphere()
        try:
            reptran = at3d.gas_absorption.Reptran()
            gas_abs = reptran.get_absorption_data(test_atmosphere,instrument=instrument.split('-')[0],band=band)
            weights = gas_abs.weights*gas_abs.irradiance
            weights /= weights.sum('wavelength')
            wavelengths = gas_abs.wavelength
        except ValueError:
            reptran = at3d.gas_absorption.Reptran(parameterization_name='thermal_modis')
            gas_abs = reptran.get_absorption_data(test_atmosphere,instrument=instrument.split('-')[0],band=band)
            weights = gas_abs.weights
            weights /= weights.sum('wavelength')
            wavelengths = gas_abs.wavelength

        BandModel.__init__(self, instrument+'-'+str(band), wavelengths.values, weights.values)

# def gaussian_cone(npixels,degree, FOV):
#     """
#     TODO
#     Number of sampling points scales as some quadratic function
#     of degree so caution should be taken when using this.
#     degree & FOV are scalars.
#     THIS ASSUMES FOV OF A PIXEL IS A CONE
#     """
#     out_mu,out_mu_weights = np.polynomial.legendre.leggauss(2*degree)
#     mus = out_mu[len(out_mu)//2:]
#     mu_weights = out_mu_weights[len(out_mu)//2:]
#
#     nphis = [int(0.9+(2*degree)*np.sqrt(1.0-mu**2)) for mu in mus]
#     phis = np.concatenate([np.cumsum([2*np.pi/nphi]*nphi)  for nphi in nphis],axis=-1)
#     big_mus = np.concatenate([[mu]*nphi for mu,nphi in zip(mus,nphis)],axis=-1)
#     weights = np.concatenate([[mu_weight/(nphi)]*nphi for mu_weight,nphi in zip(mu_weights,nphis)],axis=-1)
#
#     #scale to the correct angle interval
#     thetas = np.arccos(1.0- (1.0 - big_mus)*(1.0 - np.cos(np.deg2rad(FOV))))
#
#     #cosine weighting.
#     weights *= np.cos(thetas)
#     weights /= np.sum(weights)
#
#     thetas = np.repeat(np.expand_dims(thetas,-1),npixels,axis=-1)
#     phis = np.repeat(np.expand_dims(phis,-1),npixels,axis=-1)
#     weights = np.repeat(np.expand_dims(weights,-1),npixels,axis=-1)
#
#     return thetas,phis,weights
#
# def stochastic_cone(npixels,nrays, FOV, seed):
#     """
#     TODO
#     Generate rays that have 'equal energy' contribution (and therefore weight)
#     so that each ray is as equally worthwhile to simulate.
#     THIS ASSUMES FOV OF A PIXEL IS A CONE
#     """
#     if seed is not None:
#         np.random.seed(seed)
#
#     mu_rand = np.random.uniform(low=0.0,high=1.0,size=(nrays,npixels))
#     rand_phi = np.random.uniform(low=-np.pi,high=np.pi,size=(nrays,npixels))
#     rand_theta = np.arccos(1.0- (1.0 - mu_rand)*(1.0 - np.cos(np.deg2rad(FOV))))
#     weights = np.cos(rand_theta)/np.sum(np.cos(rand_theta),axis=0)
#     return rand_theta, rand_phi, weights
#
# def make_new_rays(theta_ref,phi_ref,theta_prime,phi_prime):
#     """
#     TODO
#     """
#     T_inv = np.array([[np.cos(theta_ref)*np.cos(phi_ref), -np.sin(phi_ref), np.sin(theta_ref)*np.cos(phi_ref)],
#                             [np.cos(theta_ref)*np.sin(phi_ref), np.cos(phi_ref), np.sin(theta_ref)*np.sin(phi_ref)],
#                              [-np.sin(theta_ref), np.zeros(theta_ref.shape), np.cos(theta_ref)]])
#
#     prime = np.array([np.sin(theta_prime)*np.cos(phi_prime),np.sin(theta_prime)*np.sin(phi_prime), np.cos(theta_prime)])
#
#     #hacks to ensure that broadcasting works properly.
#     if prime.ndim == 2:
#         prime = np.expand_dims(prime,-1)
#     prime = prime.transpose(-1,0,1)
#     T_inv = T_inv.transpose(-1,0,1)
#
#     new_ray = np.matmul(T_inv,prime)
#     new_mu = new_ray[:,-1,:]
#     new_phi = np.arctan2(new_ray[:,1,:],new_ray[:,0,:])
#     return new_mu,new_phi
#
# def add_sub_pixel_rays(sensor,FOV,inplace=True,sampling='gaussian_cone',seed=None,degree=None,
#                       nrays=None):
#     """
#     TODO
#     This is the main interface that makes the sub-pixel 'ray' variables that are actually integrated
#     by RTE solver. Therefore it needs to ALWAYS be called.
#     If FOV=0 then this just copies the 'cam' variables regardless of sampling method.
#     The new 'ray' variables can be added to the sensor object
#     which modifies it inplace or returned as a separate object.
#     TODO
#     New methods for subsampling may have a tuple of FOV to allow rectangular FOVs.
#     These new methods may be applied to arbitrary sensors.
#     Or the sub_pixel_ray generation methods may be integrated into the sensor definition
#     if they are unique to that sensor.
#     """
#
#     cam_mu = sensor.cam_mu.data
#     cam_phi = sensor.cam_phi.data
#     cam_x = sensor.cam_x.data
#     cam_y = sensor.cam_y.data
#     cam_z = sensor.cam_z.data
#     output = {}
#
#     FOV = np.atleast_1d(FOV)
#     if all(FOV==0.0):
#
#         output['ray_mu'] = ('nrays',cam_mu)
#         output['ray_phi'] = ('nrays',cam_phi)
#         output['ray_x'] = ('nrays',cam_x)
#         output['ray_y'] = ('nrays',cam_y)
#         output['ray_z'] = ('nrays',cam_z)
#         output['pixel_index'] = ('nrays', range(len(cam_mu)))
#         output['ray_weight'] = ('nrays', np.ones(len(cam_mu)))
#         output['pixel_fov'] = FOV
#
#     else:
#         if sampling == 'gaussian_cone':
#             assert degree is not None, 'degree must be defined for gaussian_cone sampling.'
#             assert len(FOV) == 1, 'FOV is a single scalar for conical sub pixel.'
#             theta_prime,phi_prime,weights = gaussian_cone(npixels=len(cam_mu),FOV=FOV,degree=degree)
#
#         elif sampling =='stochastic_cone':
#             assert nrays is not None,'seed and nrays must be defined for stochastic_cone sampling.'
#             assert len(FOV) == 1, 'FOV is a single scalar for conical sub pixel.'
#             theta_prime,phi_prime,weights = stochastic_cone(npixels=len(cam_mu),FOV=FOV,nrays=nrays,seed=seed)
#         new_mu,new_phi = make_new_rays(np.arccos(cam_mu),cam_phi,theta_prime,phi_prime)
#
#         weights_all = weights.ravel(order='F')
#         mu_all = new_mu.ravel()
#         phi_all = new_phi.ravel()
#         xs_all = np.repeat(np.expand_dims(cam_x,0),new_mu.shape[-1]).ravel()
#         ys_all = np.repeat(np.expand_dims(cam_y,0),new_mu.shape[-1]).ravel()
#         zs_all = np.repeat(np.expand_dims(cam_z,0),new_mu.shape[-1]).ravel()
#         indices = np.repeat(np.expand_dims(np.arange(len(cam_x)),1),new_mu.shape[-1]).ravel()
#
#         output['pixel_fov'] = FOV
#         output['ray_mu'] = ('nrays',mu_all)
#         output['ray_phi'] = ('nrays',phi_all)
#         output['ray_x'] = ('nrays',xs_all)
#         output['ray_y'] = ('nrays',ys_all)
#         output['ray_z'] = ('nrays',zs_all)
#         output['pixel_index'] = ('nrays', indices)
#         output['ray_weight'] = ('nrays', weights_all)
#
#     #make outputs
#     dset = xr.Dataset(data_vars=output)
#     dset.attrs['units'] = ['pixel_fov [degrees]']
#     dset.attrs['sub_pixel_sampling'] = sampling
#
#     if sampling == 'gaussian_cone':
#         dset.attrs['gaussian_sub_pixel_degree'] = degree
#     elif sampling == 'stochastic_cone':
#         dset.attrs['seed'] = seed
#         dset.attrs['n_rays'] = nrays
#
#     if inplace:
#         final = xr.merge([sensor,dset])
#         final=final.assign_attrs(dset.attrs)
#         final=final.assign_attrs(sensor.attrs)
#     else:
#         final = dset
#
#     return final
