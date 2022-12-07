import xarray as xr
import numpy as np


def test_stereo_applicability(sensors, instrument):
    """
    Test if the sensors contained in the given instrument's sensor list are
    the right format for use in the stereo pipeline (epipolar rectified).

    Parameters
    ----------
    sensors : at3d.containers.SensorsDict
        Contains the rendered imagery used as input to calculate stereo correspondences between images.
    instrument : str
        The key in `sensors` that contains the images to use.

    Returns
    -------
    x_resolution : float
        The common resolution of all of the sensors in the x direction.
    y_resolution : float
        The common resolution of all of the sensors in the y direction.
    orientation : bool
        Whether the sensor data is oriented in the x or y azimuthal directions.
        True for x.

    Raises
    ------
    ValueError
        If the group of sensors are not appropriate.
    """

    # TODO Need to add a sorting test to make sure views are sorted by zenith and azimuth before stereo.
    # As well as a sorting function to make the views ordered.

    orientations = []
    x_resolutions = []
    y_resolutions = []
    altitudes = []

    for sensor in sensors[instrument]['sensor_list']:

        if not sensor.attrs['projection'] in ('Orthographic', 'DomainTop'):
            raise ValueError(
            "All sensor projections should be Orthographic or DomainTop"
            )

        x_resolutions.append(sensor.attrs['x_resolution'])
        y_resolutions.append(sensor.attrs['y_resolution'])
        altitudes.append(sensor.cam_z.data.mean())

        y_oriented = (np.allclose(sensor.cam_phi, np.pi/2) or np.allclose(sensor.cam_phi, -1*np.pi/2))
        x_oriented = (np.allclose(sensor.cam_phi, 0) or np.allclose(sensor.cam_phi, np.pi))
        nadir = np.allclose(sensor.cam_mu,1.0)

        if nadir:
            orientations.append(2)
        elif y_oriented:
            orientations.append(1)
        elif x_oriented:
            orientations.append(0)

        if not (y_oriented or x_oriented):
            raise ValueError(
            "Sensors should be oriented in either the exact 'x' or 'y' direction."
            )
    orientation = set(orientations) - set([2])
    if len(orientation) > 1:
        raise ValueError(
        "All sensors should have the same orientation."
        )
    if len(set(x_resolutions)) > 1 or len(set(y_resolutions)) > 1:
        raise ValueError(
        "All matched images should have the same resolution"
        )

    if len(set(altitudes)) > 1:
        raise ValueError(
        "All sensor altitudes should be the same."
        )

    x_resolution = list(set(x_resolutions))[0]
    y_resolution = list(set(y_resolutions))[0]
    orientation = bool(list(orientation)[0])
    altitude = list(set(altitudes))[0]

    return x_resolution, y_resolution, orientation, altitude


def register_to_ground(sensors, instrument, x_resolution, y_resolution, altitude, fill_value=0.0):
    """
    Register all views to the ground because they will have a
    small overlapping region that can be used for finding correspondences.

    Parameters
    ----------
    sensors : at3d.containers.SensorsDict
        Contains the rendered imagery used as input to calculate stereo correspondences between images.
    instrument : str
        The key in `sensors` that contains the images to use.
    x_resolution : float
        The common x resolution of the sensors
    y_resolution : float
        The common y resolution of the sensors
    altitude : float
        The common altitude of all of the sensors.

    Returns
    -------
    registered_views : list
        The list of all of the views registered to a common ground coordinate system.
    """

    images = sensors.get_images(instrument)
    projected_views = []

    min_x = []
    max_x = []
    min_y = []
    max_y = []

    for view in images:

        projected_sensor = view.copy(deep=True)

        new_y = (view.y[0] - altitude*np.tan(np.arccos(view.mu[0,0]))*np.sin(view.phi[0,0])).data
        new_x = (view.x[:,0] - altitude*np.tan(np.arccos(view.mu[0,0]))*np.cos(view.phi[0,0])).data

        padded_I = view.I.pad({'imgdim0':(1,1), 'imgdim1': (1,1)}, mode='edge')

        projected_sensor = xr.Dataset(
            data_vars={
                'I': (['x','y'], padded_I.data)
            },
            coords={
                'x': np.append(np.append(new_x[0]-x_resolution, new_x),new_x[-1]+x_resolution),
                'y': np.append(np.append(new_y[0]-y_resolution, new_y),new_y[-1]+y_resolution),
                'mu': view.mu[0,0].data,
                'phi': view.phi[0,0].data
            }
        )
        projected_views.append(projected_sensor)
        min_x.append(projected_sensor.x.data.min())
        max_x.append(projected_sensor.x.data.max())
        min_y.append(projected_sensor.y.data.min())
        max_y.append(projected_sensor.y.data.max())

    min_x = min(min_x)
    max_x = max(max_x)

    min_y = min(min_y)
    max_y = max(max_y)

    common_x_grid = np.arange(min_x, max_x+x_resolution, x_resolution)
    common_y_grid = np.arange(min_y, max_y+y_resolution, y_resolution)

    registered_views = [
        view.interp(y=common_y_grid, x=common_x_grid, method='linear'
                   ).fillna(fill_value) for view in projected_views
    ]

    return registered_views

def compute_point_cloud(disparity, view1, view2, x_resolution, y_resolution, orientation):
    """
    Compute the positions of the point cloud from the disparities
    and sensor characteristics.

    Parameters
    ----------
    disparity : np.ndarray
        The map of disparities (view1 perspective).
    view1 : xr.Dataset
        Image data with common x, y coordinates registered to ground.
    view2 : xr.Dataset
        Image data with common x, y coordinates registered to ground.
    x_resolution : float
        The spacing of the common x coordinate
    y_resolution : float
        The spacing of the common y coordinate
    orientation : bool
        True if sensors are oriented along the y coordinate.
    """
    big_x, big_y = np.meshgrid(view1.x.data, view1.y.data, indexing='ij')

    if orientation:
        factor1 = np.sin(view1.phi.data.mean())
        factor2 = np.sin(view2.phi.data.mean())
        resolution = y_resolution
    else:
        factor1 = np.cos(view1.phi.data.mean())
        factor2 = np.cos(view2.phi.data.mean())
        resolution = x_resolution

        #big_x = big_x.T
        #big_y = big_y.T

    intersect_z = disparity*resolution/(factor1*np.tan(np.arccos(view1.mu.data.mean())) - factor2*np.tan(np.arccos(view2.mu.data.mean())))
    intersect_y = big_y + intersect_z*np.tan(np.arccos(view1.mu.data.mean()))*np.sin(view1.phi.data.mean())
    intersect_x = big_x + intersect_z*np.tan(np.arccos(view1.mu.data.mean()))*np.cos(view1.phi.data.mean())

    return intersect_x, intersect_y, intersect_z

def compute_disparities(registered_views, orientation, matcher, spacing=1, min_disparity=-500, max_disparity=500):

    disparities = []
    computed_view_pairs = []
    backflows = []
    costs = []

    for j in range(len(registered_views)):
        i = j + spacing
        if i >= len(registered_views)-spacing:
            i = j - spacing

        view1 = registered_views[j]
        view2 = registered_views[i]

        computed_view_pairs.append((view1, view2))

        image1 = view1.I.data[:,:]
        image2 = view2.I.data[:,:]

        if not orientation:
            image1 = image1.T
            image2 = image2.T

        disparity, cost, backflow = matcher.match_stereorectified_images(image1, image2,
                                                                        min_disparity=min_disparity, max_disparity=max_disparity)

        if not orientation:
            disparity = disparity.T

        disparities.append(disparity)
        backflows.append(backflow)
        costs.append(cost)

    return disparities, costs, backflows, computed_view_pairs


def compute_filtered_point_cloud(disparities, computed_view_pairs, masking_function, x_resolution, y_resolution, orientation):

    x = []
    y = []
    z = []

    for disparity, (view1, view2) in zip(disparities, computed_view_pairs):

        # Compute threshold here.
        mask = masking_function(view1.I.data, view1.mu.data, view1.phi.data)
        disparity[np.where(mask==0)] = 0.0

        intersect_x, intersect_y, intersect_z = compute_point_cloud(disparity, view1, view2, x_resolution,
                                                                   y_resolution, orientation)

        x.append(intersect_x)
        y.append(intersect_y)
        z.append(intersect_z)

    x = np.stack(x, axis=0)
    y = np.stack(y, axis=0)
    z = np.stack(z, axis=0)

    return x, y, z


def make_sensor_masks(sensors, computed_view_pairs, begin, end, subset, nray_per_pixel,
                        altitude, z_points, cloud_masker_clear,cloud_masker_conservative):

    sensor_masks = []
    no_distance_masks = []

    for i, original_sensor in enumerate(sensors[begin:end:subset]):

        mask_sensor = original_sensor.copy(deep=True)

        mask = cloud_masker_clear(original_sensor.I.data.reshape(original_sensor.image_shape.data,order='F'),
                           original_sensor.cam_mu.data[0], original_sensor.cam_phi.data[0])
        mask_sensor['cloud_mask'] = ('nrays', np.repeat(mask.ravel(order='F'),nray_per_pixel,axis=0))
        no_distance_masks.append(mask_sensor)
        mask_sensor2 = mask_sensor.copy(deep=True)

        mask = cloud_masker_conservative(original_sensor.I.data.reshape(original_sensor.image_shape.data,order='F'),
                           original_sensor.cam_mu.data[0], original_sensor.cam_phi.data[0])
        mask_sensor2['cloud_mask'] = ('nrays', np.repeat(mask.ravel(order='F'),nray_per_pixel,axis=0))

        # Find the sensors in computed_view_pairs that match.
        first_set, second_set = zip(*computed_view_pairs)

        #views = []
        distance_pair = []
        for view_set in (list(first_set), list(second_set)):

            mus = np.array([view.mu.data for view in view_set])
            phis = np.array([view.phi.data for view in view_set])
            index = np.where((mus==original_sensor.cam_mu.data[0])&(phis==original_sensor.cam_phi.data[0]))#[0][0] # all sensor.cam_mu should be the same.
            index = index[0]

            if len(index > 0):

                view = view_set[index[0]]
                distances = (altitude - z_points[i])/view.mu.data

                new_y = view.y.data + (altitude*np.tan(np.arccos(view.mu))*np.sin(view.phi)).data
                new_x = view.x.data + (altitude*np.tan(np.arccos(view.mu))*np.cos(view.phi)).data

                distance_limit = xr.DataArray(
                    data=distances,
                    dims=['x','y'],
                    coords={
                        'x': new_x,
                        'y': new_y
                    },
                )

                distance_limit = distance_limit.interp(x=original_sensor.cam_x, y=original_sensor.cam_y, method='nearest')
                distance_limit = distance_limit.fillna(0.0)
                distance_pair.append(distance_limit)

        if len(distance_pair) > 1:
            max_distance = np.maximum(distance_pair[0], distance_pair[1])
        else:
            max_distance = distance_pair[0]
        mask_sensor2['distance_limits'] = ('nrays', np.repeat(max_distance.data,nray_per_pixel,axis=0))
        sensor_masks.append(mask_sensor2)

    return sensor_masks, no_distance_masks

def calculate_reflected_transmitted(sensor_masks, nray_per_pixel, sun_distance_reflect=0.1,
                                   sun_distance_transmit=0.2):


    reflecteds = []
    transmitteds = []
    for sensor in sensor_masks:
        clear = sensor.I.data[np.where(sensor.cloud_mask.data[::nray_per_pixel] == 0)].mean()

        reflecteds.append((sensor.I.data)[np.where((sensor.sun_distance.data < sun_distance_reflect) &
                                                (sensor.cloud_mask.data[::nray_per_pixel] == 1))])

        transmitteds.append((sensor.I.data)[np.where((sensor.sun_distance.data >= sun_distance_transmit) &
                                                  (sensor.cloud_mask.data[::nray_per_pixel] == 1))])

    reflect_means = np.zeros(len(sensor_masks))
    transmit_means = np.zeros(len(sensor_masks))

    for i,(reflected, transmitted) in enumerate(zip(reflecteds, transmitteds)):

        if len(transmitted) > 1:
            transmit_mean = np.percentile(transmitted, 50)
            transmit_means[i] = transmit_mean
        else:
            transmit_means[i] = np.nan

        if (len(reflected) > 1) & (len(transmitted) > 1):
            reflect_means[i] = np.percentile((np.array(reflected)[np.where(reflected>transmit_mean)]),50)

        else:
            reflect_means[i] = np.nan

    return reflect_means, transmit_means

def initialize_extinction(carved, carved_no_distance, reflect_means, transmit_means,
                          sundistance, chi=2.0/3.0, g=0.85, fill_ext=0.01):

    condition = np.where(sundistance.sun_distance > 0.0)
    if len(condition[0]) > 0:

        sundistance_radius = np.percentile(sundistance.sun_distance.data[np.where(sundistance.sun_distance > 0.0)],90)
        condition = np.where(~np.isnan(transmit_means))
        tau_estimate = 2*chi*(np.nanmean(reflect_means[condition])/np.nanmean(transmit_means[condition]))/(1.0-g)
        tau2 = (2*chi/(1-g))*np.pi*np.nanmean(reflect_means)/(1.0-np.pi*np.nanmean(reflect_means))

        central_ext_estimate = tau_estimate/sundistance_radius
    else:
        central_ext_estimate = fill_ext

    zs = carved.z.where(carved.mask==1)
    zmin = zs.min()
    zmax = zs.max()

    adiabatic_profile = ((zs-zmin)**(2.0/3.0))
    scale_factor = central_ext_estimate/adiabatic_profile.mean()
    initial_profile = adiabatic_profile*scale_factor

    initial_extinction = initial_profile.where(~((carved_no_distance.mask==1) & (np.isnan(initial_profile))), 0.01).transpose('x','y','z')

    return initial_extinction
