import xarray as xr
import numpy as np

def make_sensor_dataset(x,y,z,mu,phi):
    """
    TODO
    """
    dataset = xr.Dataset(
            data_vars={
                'observable_list': (['nimage','observables'],np.array([[True,True,False,False]])),
                'cam_x': (['total_pixels'], x.ravel().astype(np.float64)),
                'cam_y': (['total_pixels'], y.ravel().astype(np.float64)),
                'cam_z': (['total_pixels'], z.ravel().astype(np.float64)),
                'cam_mu'  : (['total_pixels'], mu.ravel()),
                'cam_phi' : (['total_pixels'], phi.ravel()),
                'super_pixel_index': (['total_pixels'], np.append(np.arange(mu.ravel().shape[0]//2),
                                                                  np.arange(mu.ravel().shape[0]//2)).astype(np.int64)) ,
                'super_pixel_weight': (['total_pixels'], np.ones(mu.ravel().shape)),
                'image_shape': (['nimage','image_dims'], [[scatterer.x.shape[0], scatterer.y.shape[0]]])
            },
        coords = {'observables': np.array(['I', 'Q', 'U', 'V']),
                 }
    )
    #TODO remove the 'nimage' dimension.

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

    concat_observable = xr.concat([data.observable_list for data in list_of_sensor_datasets],dim='nimage')
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
