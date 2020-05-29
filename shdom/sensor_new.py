import xarray as xr
import numpy as np
import shdom

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

    #Various other information could be added
    #about the grouping and relation of pixels to one another,
    #but this is the bare minimum for rendering.
    return dataset


def merge_sensor_list(list_of_sensor_datasets):
    """
    TODO
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
