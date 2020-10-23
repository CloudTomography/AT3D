import xarray as xr
import numpy as np
import itertools
import inspect

def make_sensor_dataset(x, y, z, mu, phi, stokes, wavelength):
    """
    TODO
    """
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
        coords = {'stokes_index': np.array(['I', 'Q', 'U', 'V'])}
    )

    return dataset

def merge_sensors_forward(sensors):
    """
    TODO
    """

    var_list = ['ray_x','ray_y','ray_z','ray_mu','ray_phi','ray_weight','pixel_index']

    output = {}
    for var in var_list:
        concatenated = xr.concat([sensor[var] for sensor in sensors], dim='nrays')
        output[var] = ('nrays', concatenated)


    output['stokes'] = xr.concat([sensor.stokes for sensor in sensors],dim='nimage')
    output['rays_per_image'] = ('nimage', np.array([sensor.sizes['nrays'] for sensor in sensors]))
    output['rays_per_pixel'] = ('npixels', np.concatenate([np.unique(sensor.pixel_index,return_counts=True)[1]\
                                                          for sensor in sensors]))

    merged_dataset = xr.Dataset(data_vars=output)

    return merged_dataset

def split_sensor_rays(merged):
    """
    TODO
    """
    list_of_unmerged = []
    count = 0

    for i in range(merged.sizes['nimage']):
        split_index = merged.rays_per_image[i].data
        split = merged.sel({'nrays': slice(count, count+split_index), 'nimage':i})

        list_of_unmerged.append(split)
        count += split_index
    return list_of_unmerged

def split_sensor_pixels(merged):
    """
    TODO
    """
    pixel_inds = np.cumsum(np.concatenate([np.array([0]), merged.rays_per_pixel.data])).astype(np.int)
    pixels_per_image = []
    for ray_image in merged.rays_per_image.data:
        pixels_per_image.append(np.where(pixel_inds == ray_image)[0][0])

    list_of_unmerged = []
    count = 0
    for i in range(merged.sizes['nimage']):
        split_index = pixels_per_image[i]#merged.rays_per_image[i].data
        split = merged.sel({'npixels': slice(count, count+split_index), 'nimage':i})
        list_of_unmerged.append(split)
        count += split_index
    return list_of_unmerged

def get_image(sensor):
    """
    TODO
    Make a 2D sensor for visualization purposes.
    """
    if not hasattr(sensor, 'image_shape'):
        raise ValueError("Sensor dataset does not have an 'image_shape' variable. A 2D image cannot be formed.")
    img_data = xr.Dataset(
        data_vars={
            'x': (['imgdim0', 'imgdim1'],sensor.cam_x.data.reshape(sensor.image_shape.data,order='F')),
            'y': (['imgdim0', 'imgdim1'],sensor.cam_y.data.reshape(sensor.image_shape.data,order='F')),
            'mu': (['imgdim0', 'imgdim1'],sensor.cam_mu.data.reshape(sensor.image_shape.data,order='F')),
            'phi': (['imgdim0', 'imgdim1'],sensor.cam_phi.data.reshape(sensor.image_shape.data,order='F'))
        }
    )
    for stokes in sensor.stokes_index:
        if sensor.stokes.sel({'stokes_index': stokes}):
            img_data[str(stokes.data)] = (['imgdim0','imgdim1'],sensor[str(stokes.data)].data.reshape(
                                        sensor.image_shape.data,order='F'))
    return img_data

def get_observables(sensor, rendered_rays):
    """
    TODO
    Currently averages over the sub_pixel rays and only takes
    the observables that are needed.
    """
    if not sensor.use_subpixel_rays:
        merge_list = []
    else:
        merge_list = [sensor.pixel_index]
    for stokes in sensor.stokes_index:
        if rendered_rays['stokes'].sel({'stokes_index':stokes}):
            weighted_stokes = (sensor.ray_weight*rendered_rays[str(stokes.data)])
            weighted_stokes.name = '{}'.format(str(stokes.data))
            merge_list.append(weighted_stokes)

    temp = xr.merge(merge_list)

    if not sensor.use_subpixel_rays:
        #avoid expensive groupby if no subpixel rays are used.
        for stokes in temp.data_vars:
            sensor[stokes] = ('npixels', temp[stokes].data)
    else:
        #only do the group_by once, not for every stokes.
        observed = temp.groupby('pixel_index').sum()

        for stokes in observed.data_vars:
            sensor[stokes] = ('npixels', observed[stokes].data)

    return sensor

def merge_sensors_inverse(sensors, solver):
    """
    TODO
    """
    var_list = ['ray_x','ray_y','ray_z','ray_mu','ray_phi','ray_weight','pixel_index']

    output = {}
    for var in var_list:
        concatenated = xr.concat([sensor[var] for sensor in sensors], dim='nrays')
        output[var] = ('nrays', concatenated)

    output['stokes'] = xr.concat([sensor.stokes for sensor in sensors],dim='nimage')
    output['rays_per_image'] = ('nimage', np.array([sensor.sizes['nrays'] for sensor in sensors]))
    output['rays_per_pixel'] = ('npixels', np.concatenate([np.unique(sensor.pixel_index,return_counts=True)[1]\
                                                          for sensor in sensors]))
    stokes_weights = []
    stokes_datas = []
    for sensor in sensors:
        stokes_weight = np.zeros((solver._nstokes, sensor.sizes['npixels']))
        stokes_data = np.zeros((solver._nstokes, sensor.sizes['npixels']))
        for i,stokes in enumerate(sensor.stokes_index):
            if stokes in list(sensor.data_vars):
                stokes_weight[i,:] = 1.0
                stokes_data[i,:] = sensor[str(stokes.data)]

        stokes_weights.append(stokes_weight)
        stokes_datas.append(stokes_data)

    output['stokes_weights'] = (['nstokes','npixels'], np.concatenate(stokes_weights,axis=-1))
    output['measurement_data'] = (['nstokes','npixels'], np.concatenate(stokes_datas,axis=-1))

    merged_sensors = xr.Dataset(data_vars=output)

    #TODO HARDCODED UNCERTAINTIES
    merged_sensors['uncertainties'] = xr.DataArray(data= np.ones((merged_sensors.sizes['nstokes'],
                                merged_sensors.sizes['nstokes'], merged_sensors.sizes['npixels'])),
                                dims=['nstokes', 'nstokes2', 'npixels']) #NB Be aware that repeated dims cause errors. so the second dim is set to 'nstokes2' Issue 1378 on xarray.
    return merged_sensors

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


def parse_sub_pixel_ray_args(sub_pixel_ray_args):
    """
    TODO
    """
    try:
        subpixel_ray_params = inspect.signature(sub_pixel_ray_args['method'])
    except TypeError as err:
        raise TypeError("sub_pixel_ray_args 'method' must be a"
                        "callable object not '{}' of type '{}'".format(
                            sub_pixel_ray_args['method'],type(sub_pixel_ray_args['method']))) from err
    subpixel_ray_kwargs_x = {}
    subpixel_ray_kwargs_y = {}
    subpixel_ray_params = inspect.signature(sub_pixel_ray_args['method']).parameters
    for name,value in sub_pixel_ray_args.items():
        if name is not 'method':
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
            except:
                raise ValueError("Invalid kwarg '{}' passed to the "
                                 "sub_pixel_ray_args['method'] callable. '{}'".format(
                                     name, sub_pixel_ray_args['method'].__name__))
    return sub_pixel_ray_args['method'], subpixel_ray_kwargs_x, subpixel_ray_kwargs_y

def add_null_subpixel_rays(sensor):
    """
    TODO
    """
    sensor['ray_mu'] = ('nrays',sensor.cam_mu.data)
    sensor['ray_phi'] = ('nrays',sensor.cam_phi.data)
    sensor['ray_x'] = ('nrays',sensor.cam_x.data)
    sensor['ray_y'] = ('nrays',sensor.cam_y.data)
    sensor['ray_z'] = ('nrays',sensor.cam_z.data)
    sensor['pixel_index'] = ('nrays', range(len(sensor.cam_mu.data)))
    sensor['ray_weight'] = ('nrays', np.ones(len(sensor.cam_mu.data)))
    sensor['use_subpixel_rays'] = False
    return sensor

def stochastic(npixels,nrays, seed=None):
    """
    TODO
    """
    if seed is not None:
        np.random.seed(seed)
    position_perturbations = np.random.uniform(low=-1.0,high=1.0, size=(npixels, nrays))
    weights = np.ones((npixels, nrays))/nrays
    return position_perturbations, weights

def gaussian(npixels, degree):
    """
    TODO
    """
    out_mu,out_mu_weights = np.polynomial.legendre.leggauss(degree)
    position_perturbations = np.repeat(out_mu[np.newaxis,:], npixels, axis=0)
    weights = np.repeat(out_mu_weights[np.newaxis,:], npixels, axis=0)/np.sum(out_mu_weights)
    return position_perturbations, weights

def orthographic_projection(wavelength, bounding_box, x_resolution, y_resolution, azimuth, zenith,
                            altitude='TOA', stokes='I', sub_pixel_ray_args={'method': None}):
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

    if np.isclose(x_s,xmin) or np.isclose(x_s,xmax):
        x = np.arange(x_s, x_e + 1e-6, x_resolution)
    elif np.isclose(x_e,xmin) or np.isclose(x_e,xmax):
        x = -1*np.arange(-1*x_e,-1*x_s+1e-6, x_resolution)[::-1]
    if np.isclose(y_s,ymin) or np.isclose(y_s,ymax):
        y = np.arange(y_s, y_e + 1e-6, y_resolution)
    elif np.isclose(y_e,ymin) or np.isclose(y_e,ymax):
        y = -1*np.arange(-1*y_e,-1*y_s+1e-6, y_resolution)[::-1]
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
    if sub_pixel_ray_args['method'] is not None:

        #generate the weights and perturbations to the pixel positions in the image plane.
        sub_pixel_ray_method, subpixel_ray_kwargs_x,subpixel_ray_kwargs_y = parse_sub_pixel_ray_args(sub_pixel_ray_args)
        position_perturbations_x, weights_x = sub_pixel_ray_method(x.size, **subpixel_ray_kwargs_x)
        position_perturbations_y, weights_y = sub_pixel_ray_method(y.size, **subpixel_ray_kwargs_y)

        #merge the two dimensions
        perturbations_x = np.repeat(position_perturbations_x[...,np.newaxis]*x_resolution/2.0,
                                    position_perturbations_y.shape[-1], axis=-1)
        perturbations_y = np.repeat(position_perturbations_y[...,np.newaxis,:]*y_resolution/2.0,
                                    position_perturbations_x.shape[-1], axis=-2)
        big_weightx = np.repeat(weights_x[...,np.newaxis], weights_y.shape[-1],axis=-1)
        big_weighty = np.repeat(weights_y[...,np.newaxis,:], weights_x.shape[-1],axis=-2)

        #apply perturbations to original image plane coordinates.
        x_ray = (x.ravel()[:,np.newaxis,np.newaxis] + perturbations_x).ravel()
        y_ray = (y.ravel()[:,np.newaxis,np.newaxis] + perturbations_y).ravel()
        z_ray = np.repeat(np.repeat(z.ravel()[:,np.newaxis,np.newaxis], perturbations_x.shape[-2], axis=-2),
                          perturbations_y.shape[-1], axis=-1).ravel()
        mu_ray = np.repeat(np.repeat(mu.ravel()[:,np.newaxis,np.newaxis], perturbations_x.shape[-2], axis=-2),
                          perturbations_y.shape[-1], axis=-1).ravel()
        phi_ray = np.repeat(np.repeat(phi.ravel()[:,np.newaxis,np.newaxis], perturbations_x.shape[-2], axis=-2),
                          perturbations_y.shape[-1], axis=-1).ravel()
        #make the pixel indices and ray weights.
        pixel_index = np.repeat(np.repeat(range(len(sensor.cam_mu.data)),
                                          weights_x.shape[-1]), weights_y.shape[-1])
        ray_weight = (big_weightx*big_weighty).ravel()
        #update ray variables to sensor dataset.
        sensor['ray_mu'] = ('nrays',mu_ray)
        sensor['ray_phi'] = ('nrays',phi_ray)
        sensor['ray_x'] = ('nrays',x_ray)
        sensor['ray_y'] = ('nrays',y_ray)
        sensor['ray_z'] = ('nrays',z_ray)
        sensor['pixel_index'] = ('nrays', pixel_index)
        sensor['ray_weight'] = ('nrays', ray_weight)
        sensor['use_subpixel_rays'] = True
        sensor['subpixel_ray_args'] = sub_pixel_ray_args

    else:
        #duplicate ray variables to sensor dataset.
        sensor = add_null_subpixel_rays(sensor)
    return sensor

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
