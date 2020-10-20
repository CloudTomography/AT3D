
from joblib import Parallel, delayed
import shdom
import numpy as np
from collections import OrderedDict
import xarray as xr

def make_forward_sensors(sensors):

    forward_sensors = OrderedDict()
    for key,sensor in sensors.items():
        forward_sensor= OrderedDict()
        forward_sensor['sensor_list'] = [single_sensor.copy(deep=True) for single_sensor in sensor['sensor_list']]
        forward_sensors[key] = forward_sensor

    return forward_sensors

def get_unique_wavelengths(sensor_dict):
    """
    TODO
    """
    wavelength_list = []
    for key,instrument in sensor_dict.items():
        for sensor in instrument['sensor_list']:
            wavelength_list.append(sensor.wavelength)

    return np.unique(wavelength_list)

def parallel_solve(solvers, n_jobs=1, mpi_comm=None, maxiter=100, verbose=True, init_solution=True,
                    setup_grid=True):
    """
    TODO
    """
    if mpi_comm is not None:
        raise NotImplementedError

    else:
        if n_jobs==1:
            for solver in solvers.values():
                solver.solve(maxiter=maxiter, init_solution=init_solution, verbose=verbose,setup_grid=setup_grid)
        else:
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(solver.solve, check_pickle=False)(maxiter, init_solution, verbose,setup_grid=setup_grid) for solver in solvers.values())

def render_one_solver(solver,merged_sensor, sensor_mapping, sensors, maxiter=100, verbose=True, n_jobs=1, init_solution=True,
                        setup_grid=True):
    """
    TODO
    """
    solver.solve(maxiter=maxiter,init_solution=True,verbose=verbose,setup_grid=setup_grid)
    split = np.array_split(np.arange(merged_sensor.sizes['nrays']),n_jobs)
    start_end = [(a.min(),a.max()) for a in split]
    out = Parallel(n_jobs=n_jobs)(delayed(solver.integrate_to_sensor, check_pickle=False)(merged_sensor.sel(nrays=slice(start,end+1))) for start,end in start_end)

    var_list_nray = [str(name) for name in out[0].data_vars if str(name) not in ('rays_per_image', 'stokes',
                                                                                'rays_per_pixel')]
    merged = {}
    for var in var_list_nray:
        concatenated= xr.concat([data[var] for data in out], dim='nrays')
        merged[var] = concatenated
    merged['stokes'] = out[0].stokes
    merged['rays_per_image'] = out[0].rays_per_image
    integrated_rays = xr.Dataset(merged)

    rendered_rays = shdom.sensor.split_sensor_rays(integrated_rays)
    for i,rendered_ray in enumerate(rendered_rays):
        mapping = sensor_mapping[i]
        sensor = sensors[mapping[0]]['sensor_list'][mapping[1]]
        unused = shdom.sensor.get_observables(sensor, rendered_ray)

def sort_sensors(sensors, solvers, merge):

    rte_sensors = OrderedDict()
    sensor_mapping = OrderedDict()

    for key, solver in solvers.items():
        sensor_list = []
        mapping_list = []
        for instrument,instrument_data in sensors.items():
            for i,sensor in enumerate(instrument_data['sensor_list']):
                if key == sensor.wavelength:
                    sensor_list.append(sensor)
                    mapping_list.append((instrument,i))
        if merge == 'forward':
            merged = shdom.sensor.merge_sensors_forward(sensor_list)
        elif merge == 'inverse':
            merged = shdom.sensor.merge_sensors_inverse(sensor_list, solver)
        rte_sensors[key] = merged
        sensor_mapping[key] = mapping_list
    return rte_sensors, sensor_mapping

def combine_to_medium(scatterers):
    """
    TODO
    """
    mediums = OrderedDict()
    for key in np.atleast_1d(scatterers)[0].keys():
        scatterer_list = []
        for scatterer in scatterers:
            scatterer_list.append(scatterer[key])
        mediums[key] = scatterer_list
    return mediums

def subdivide_raytrace_jobs(rte_sensors, solvers, n_jobs):
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


def get_measurements(solvers,sensors, n_jobs=1, mpi_comm=None, destructive=False, maxiter=100,verbose=True, init_solution=True,
                    setup_grid=True):
    """
    In-place modification of sensors and solvers.
    """
    rte_sensors, sensor_mapping = sort_sensors(sensors, solvers, 'forward')
    parallel_solve(solvers, n_jobs=n_jobs, mpi_comm=mpi_comm, maxiter=maxiter, verbose=verbose, init_solution=init_solution,
                        setup_grid=setup_grid)
    if mpi_comm is not None:
        raise NotImplementedError

    else:

        if destructive:
            for key in solvers.keys():

                render_one_solver(solvers.pop(key), rte_sensors[key], sensor_mapping[key],
                                sensors, maxiter, verbose, n_jobs)

        else:
            if n_jobs==1:
                out = [solvers[key].integrate_to_sensor(rte_sensors[key] for key in solvers.keys()]

            else:
                #decide on the division of n_jobs among solvers based on total number of pixels.
                #Note that the number of n_jobs here doesn't have to be the number of workers but can instead
                #be the number of subdivided job. This could be modified to ensure that all jobs are the correct size.
                #as is, this will make slightly more tasks than there are workers (n_jobs).
                keys, ray_start_end, pixel_start_end = subdivide_raytrace_jobs(rte_sensors, solvers, n_jobs)

                out= Parallel(n_jobs=n_jobs)(
                                delayed(solvers[key].integrate_to_sensor)(rte_sensors[key].sel(
                                nrays=slice(ray_start,ray_end),npixels=slice(pix_start, pix_end)))
                            for key, (ray_start,ray_end),(pix_start,pix_end) in zip(keys, ray_start_end, pixel_start_end))

            #re arrange output so it goes in the original sensor objects.
            for key in np.unique(keys):

                #group output by solver.
                indices = np.where(key == np.array(keys))[0]
                out_key = []
                for index in indices:
                    out_key.append(out[index])

                var_list_nray = [str(name) for name in out_key[0].data_vars if str(name) not in ('rays_per_image', 'stokes',
                                                                                            'rays_per_pixel')]
                merged = {}
                for var in var_list_nray:
                    concatenated= xr.concat([data[var] for data in out_key], dim='nrays')
                    merged[var] = concatenated
                merged['stokes'] = out_key[0].stokes
                merged['rays_per_image'] = out_key[0].rays_per_image
                integrated_rays = xr.Dataset(merged)
                sensor_mapping = sensor_mappings[key]
                rendered_rays = shdom.sensor.split_sensor_rays(integrated_rays) #splits into individual images.
                for i,rendered_ray in enumerate(rendered_rays):
                    mapping = sensor_mapping[i]
                    sensor = sensors[mapping[0]]['sensor_list'][mapping[1]]
                    shdom.sensor.get_observables(sensor, rendered_ray)     #modifies sensor in-place. Also averages over sub-pixel rays.
