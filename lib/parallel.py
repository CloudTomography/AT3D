"""
TODO
"""
import inspect
import warnings
import numpy as np
from collections import OrderedDict
from joblib import Parallel, delayed

import pyshdom.gradient

def parallel_gradient(solvers, rte_sensors, sensor_mappings, forward_sensors, gradient_fun,
                     mpi_comm=None, n_jobs=1, **kwargs):
    """
    TODO
    """
    #organize **kwargs safely.
    grad_kwargs = {}
    grad_params = inspect.signature(gradient_fun).parameters
    for name,value in kwargs.items():
        if name in grad_params.keys():
            if grad_params[name].kind in (grad_params[name].POSITIONAL_OR_KEYWORD, grad_params[name].KEYWORD_ONLY,
                    grad_params[name].VAR_KEYWORD):
                grad_kwargs[name] = value
            else:
                warnings.warn("kwarg '{}' passed to pyshdom.gradient.parallel_gradient is unused by \
                                gradient_fun '{}''".format(name, gradient_fun.__name__))
        else:
            warnings.warn("kwarg '{}' passed to pyshdom.gradient.parallel_gradient is unused by \
                            gradient_fun '{}''".format(name, gradient_fun.__name__))

    if mpi_comm is not None:
        raise NotImplementedError

    else:
        if n_jobs == 1 or n_jobs >= forward_sensors.npixels:
            out = [gradient_fun(solvers[key],rte_sensors[key], **grad_kwargs) for key in solvers]
            keys = list(solvers.keys())
        else:
            #decide on the division of n_jobs among solvers based on total number of rays.
            keys, ray_start_end, pixel_start_end = subdivide_raytrace_jobs(rte_sensors, n_jobs)

            out = Parallel(n_jobs=n_jobs, backend='threading')(
                delayed(gradient_fun)(solvers[key],rte_sensors[key].sel(nrays=slice(ray_start,ray_end),
                                        npixels=slice(pix_start, pix_end)),**grad_kwargs)
                                for key, (ray_start,ray_end),(pix_start,pix_end) in zip(keys, ray_start_end, pixel_start_end))

        gradient = np.stack([i[0] for i in out],axis=-1)
        loss = np.array([i[1] for i in out])
        forward_model_output = [i[2] for i in out]
        #modify forward sensors in place to contain updated forward model estimates.
        forward_sensors.add_measurements_inverse(sensor_mappings, forward_model_output, keys)

    #special treatment of jacobian out.
    if gradient_fun == pyshdom.gradient.jacobian:
        jacobian_list = [i[-1] for i in out]
        return loss, gradient, jacobian_list
    else:
        return loss, gradient

def subdivide_raytrace_jobs(rte_sensors, n_jobs):
    """
    TODO
    """
    #loose distribution of workers by wavelength based on number
    #of rays at each wavelength.
    render_jobs = OrderedDict()
    ray_count = 0
    for key, merged_sensor in rte_sensors.items():
        ray_count+= merged_sensor.sizes['nrays']
        render_jobs[key] = ray_count
    for key,render_job in render_jobs.items():
        render_jobs[key] = max(np.round(render_job/ray_count * n_jobs).astype(np.int), 1)

    #find the ray indices to split each sensor at.
    keys = []
    pixel_start_end = []
    ray_start_end = []
    for key, merged_sensor in rte_sensors.items():
        split = np.array_split(np.arange(merged_sensor.sizes['nrays']),render_jobs[key])
        start_end = [(i.min(),i.max()) for i in split]
        #print(start_end)
        #adjust start and end indices so that rays are grouped by their parent pixel.
        index_diffs = np.append(merged_sensor.pixel_index.diff(dim='nrays').data,1) #add last end pixel.
        ends = np.where(index_diffs==1)[0] + 1
        updated_start_end = []
        updated_start_end.append((0, ends[np.abs(ends - start_end[0][1]).argmin()]))

        for i in range(1,len(start_end)):
            new_start = updated_start_end[i-1][1]
            if i < len(start_end) - 1:
                new_end = ends[np.abs(ends - start_end[i][1]).argmin()]
            else:
                new_end = start_end[i][1] + 1

            updated_start_end.append((new_start, new_end))

        ray_start_end.extend(updated_start_end)
        pixel_inds = np.cumsum(np.concatenate([np.array([0]), merged_sensor.rays_per_pixel.data])).astype(np.int)
        #return updated_start_end, pixel_inds
        for start,end in updated_start_end:
            a = np.where(pixel_inds == start)[0][0]
            b=np.where(pixel_inds == end)[0][0]

            pixel_start_end.append((a,b))
            keys.append(key)

    return keys, ray_start_end, pixel_start_end
