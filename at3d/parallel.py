"""
This  module contains functions that are used to parallelize
internal components of the processing. This includes parallelizing
the gradient calculation and also the subdivision of a sensor
into `n_jobs` for parallelization.
"""
import inspect
import warnings
from collections import OrderedDict
from joblib import Parallel, delayed
import numpy as np

def parallel_gradient(solvers, rte_sensors, sensor_mappings, forward_sensors, gradient_fun,
                      mpi_comm=None, n_jobs=1, **kwargs):
    """
    Parallelizes the evaluation of a gradient across several solvers and sensors.

    Evaluates `gradient_fun` for solvers by subdividing `rte_sensors` into `n_jobs`
    distributed to `n_jobs` parallel workers.

    Parameters
    ----------
    solvers : at3d.containers.SolversDict
        Contains the solvers for which `gradient_fun` should be evaluated.
    rte_sensors : OrderedDict
        Contains the sensors, grouped by solver key for which `gradient_fun` should
        be evavluated.
    sensor_mappings : OrderedDict
        Contains the mapping between of pixels between the `rte_sensors` which are
        grouped by solver, and the `forward_sensors`, which are grouped by instrument
        name and sensor index.
    forward_sensors : at3d.containers.SensorsDict
        Container for the synthetic measurements that will be stored by instrument
        key and sensor index.
    gradient_fun : callable
        When evaluated, this function will return the loss, gradient and
        synthetic measurements.
    mpi_comm : None
        TODO
        MPI-based parallelization for the gradient is not yet implemented.
    n_jobs : int
        The number of parallel workers for which multi-threading will be used
        to parallelize the evaluation of `gradient_fun`.
    kwargs : dict
        Arguments to `gradient_fun`.

    Returns
    -------
    loss : float
        The loss evaluated by `gradient_fun`
    gradient : np.ndarray, float
        The gradient of a specified cost function (determined by `gradient_fun`).
    jacobian : np.ndarray, optional
        Only returned if `gradient_fun`=`at3d.gradient.jacobian`
    """
    #organize **kwargs safely.
    grad_kwargs = {}
    grad_params = inspect.signature(gradient_fun).parameters
    for name, value in kwargs.items():
        if name in grad_params.keys():
            if grad_params[name].kind in (
                    grad_params[name].POSITIONAL_OR_KEYWORD,
                    grad_params[name].KEYWORD_ONLY,
                    grad_params[name].VAR_KEYWORD):
                grad_kwargs[name] = value
            else:
                warnings.warn(
                    "kwarg '{}' passed to at3d.gradient.parallel_gradient is unused by"
                    "'{}'".format(name, gradient_fun.__name__))
        else:
            warnings.warn("kwarg '{}' passed to at3d.gradient.parallel_gradient is unused by"
                            "'{}'".format(name, gradient_fun.__name__))

    if mpi_comm is not None:
        raise NotImplementedError

    else:
        if n_jobs == 1 or n_jobs >= forward_sensors.npixels:
            out = [gradient_fun(solvers[key], rte_sensors[key], **grad_kwargs) for key in solvers]
            keys = list(solvers.keys())
        else:
            #decide on the division of n_jobs among solvers based on total number of rays.
            keys, ray_start_end, pixel_start_end = subdivide_raytrace_jobs(rte_sensors, n_jobs)
            #for (key, sensor) in rte_sensors.items():
                #print(key, sensor.npixels.size, sensor.nrays.size, sensor.pixel_index.data.min(), sensor.pixel_index.data.max())
            #for key, (ray_start, ray_end), (pix_start, pix_end) in zip(keys, ray_start_end, pixel_start_end):
                #print('ray', key, ray_start, ray_end)
                #print('pixel', key, pix_start, pix_end)
            out = Parallel(n_jobs=n_jobs, backend='threading')(
                delayed(gradient_fun)(
                    solvers[key],
                    rte_sensors[key].sel(
                        nrays=slice(ray_start, ray_end),
                        npixels=slice(pix_start, pix_end)),
                    **grad_kwargs)
                for key, (ray_start, ray_end), (pix_start, pix_end) in
                zip(keys, ray_start_end, pixel_start_end)
                )

        gradient = np.stack([i[0] for i in out], axis=-1)
        loss = np.array([i[1] for i in out])
        forward_model_output = [i[2] for i in out]

        other_output = []
        for i in range(3, len(out[0])):
            if out[0][i] is not None:
                other_output.append([entry[i] for entry in out])

        #modify forward sensors in place to contain updated forward model estimates.
        forward_sensors.add_measurements_inverse(sensor_mappings, forward_model_output, keys)

    return loss, gradient, other_output

def subdivide_raytrace_jobs(rte_sensors, n_jobs, job_factor=1):
    """
    Subdivides a sensor for parallelized processing.

    The pixels and rays in `rte_sensors` are subdivided into roughly
    `n_jobs` pieces, taking care to make sure that all rays are grouped
    with their corresponding pixels.

    Parameters
    ----------
    rte_sensors : OrderedDict
        A group of merged sensors corresponding to a particular solver.RTE's key.
        See at3d.containers.sort_sensors()
    n_jobs : int
        The number of groups to subdivide each sensor into.
    job_factor : int
        A factor which divides up the jobs into more pieces which would improve
        load balancing by making tasks similar sizes with a trade off against
        overhead of workers calling the function multiple times.
    """
    #loose distribution of workers by sensor key based on number
    #of rays at each sensor key.
    render_jobs = OrderedDict()
    ray_count = 0
    for key, merged_sensor in rte_sensors.items():
        ray_count += merged_sensor.sizes['nrays']
        render_jobs[key] = merged_sensor.sizes['nrays']
    for key, render_job in render_jobs.items():
        render_jobs[key] = max(np.ceil(render_job/ray_count * n_jobs * job_factor).astype(int), 1)
    #find the ray indices to split each sensor at.
    keys = []
    pixel_start_end = []
    ray_start_end = []
    for key, merged_sensor in rte_sensors.items():
        split = np.array_split(np.arange(merged_sensor.sizes['nrays'] + 1), render_jobs[key])
        start_end = [(i.min(), i.max()) for i in split]
        #adjust start and end indices so that rays are grouped by their parent pixel.
        index_diffs = np.append(
            merged_sensor.pixel_index.diff(dim='nrays').data,
            1
            ) #add last end pixel.
        ends = np.where(index_diffs == 1)[0] + 1
        updated_start_end = []
        new_start = 0
        for i, (start, end) in enumerate(start_end):
            updated_start_end.append((new_start, ends[np.abs(ends-end).argmin()]))
            new_start = updated_start_end[i][1]
        ray_start_end.extend(updated_start_end)

        pixel_inds = np.cumsum(np.concatenate(
            [np.array([0]), merged_sensor.rays_per_pixel.data])).astype(int)
        for i, (start, end) in enumerate(updated_start_end):
            final_start = np.where(pixel_inds == start)[0][0]
            final_end = np.where(pixel_inds == end)[0][0]
            pixel_start_end.append((final_start, final_end))
            keys.append(key)

        assert ray_start_end[-1][1] == merged_sensor.nrays.size, "Unexpectedly, number of rays after parallelization is not equivalent to original"
        assert pixel_start_end[-1][1] == merged_sensor.npixels.size, "Unexpectedly, number of pixels after parallelization is not equivalent to original"

    return keys, ray_start_end, pixel_start_end
