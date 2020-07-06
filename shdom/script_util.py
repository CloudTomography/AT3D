
from joblib import Parallel, delayed
import shdom
import numpy as np
from collections import OrderedDict
import xarray as xr

def get_unique_wavelengths(sensor_dict):
    """
    TODO
    """
    wavelength_list = []
    for key,instrument in sensor_dict.items():
        for sensor in instrument['sensor_list']:
            wavelength_list.append(sensor.wavelength)

    return np.unique(wavelength_list)

def parallel_solve(solvers, n_jobs=1, mpi_comm=None, maxiter=100, verbose=True):
    """
    TODO
    """
    if mpi_comm is not None:
        raise NotImplementedError

    else:
        if n_jobs==1:
            for solver in solvers.values():
                solver.solve(maxiter=maxiter, verbose=verbose)
        else:
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(solver.solve, check_pickle=False)(maxiter, verbose) for solver in solvers.values())

def render_one_solver(solver,merged_sensor, sensor_mapping, sensors, maxiter=100, verbose=True, n_jobs=1):
    """
    TODO
    """
    solver.solve(maxiter=maxiter,verbose=verbose)

    split = np.array_split(np.arange(merged_sensor.sizes['nrays']),n_jobs)
    start_end = [(a.min(),a.max()) for a in split]
    out = Parallel(n_jobs=n_jobs, backend="threading")(delayed(solver.integrate_to_sensor)(merged_sensor.sel(nrays=slice(start,end+1))) for start,end in start_end)

    var_list_nray = [str(name) for name in out[0].data_vars if str(name) not in ('rays_per_image', 'stokes')]
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


def get_measurements(solvers,sensors, n_jobs=1, mpi_comm=None, destructive=False, maxiter=100,verbose=True):
    """
    In-place modification of sensors and solvers.
    """
    rte_sensors, sensor_mapping = sort_sensors(sensors, solvers, 'forward')

    if mpi_comm is not None:
        raise NotImplementedError

    else:

        if destructive:
            for key in solvers.keys():

                render_one_solver(solvers.pop(key), rte_sensors[key], sensor_mapping[key],
                                sensors, maxiter, verbose, n_jobs)

        else:
            #decide on the division of n_jobs among solvers based on total number of pixels.
            render_jobs = OrderedDict()
            pixel_count = 0
            for key, merged_sensor in rte_sensors.items():
                pixel_count += merged_sensor.sizes['nrays']
                render_jobs[key] = pixel_count

            for key,render_job in render_jobs.items():
                render_jobs[key] = np.round(render_job/pixel_count * n_jobs).astype(np.int)

            Parallel(n_jobs=len(list(solvers.keys())), backend="threading")(
                delayed(render_one_solver, check_pickle=False)(solvers[key],rte_sensors[key], sensor_mapping[key],
                sensors,maxiter, verbose, render_jobs[key]) for key in solvers.keys())
