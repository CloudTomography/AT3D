
from joblib import Parallel, delayed
import shdom


def get_unique_wavelengths(sensordict):
    """
    TODO
    """
    wavelength_list = []
    for key,instrument in sensordict.items():
        for sensor in instrument['sensor_list']:
            wavelength_list.append(sensor.wavelength)

    return np.unique(wavelength_list)

def parallel_solve(solvers, n_jobs=1, mpi_comm=None, maxiter=100, verbose=True):
    """
    TODO
    See 'get_measurements' for more TODOs
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
    integrated_rays = solver.integrate_to_sensor(merged_sensor, n_jobs=n_jobs)
    rendered_rays = shdom.sensor.split_sensor_rays(integrated_rays)
    for i,rendered_ray in enumerate(rendered_rays):
        mapping = sensor_mapping[i]
        sensor = sensors[mapping[0]]['sensor_list'][mapping[1]]
        unused = shdom.sensor.get_observables(sensor, rendered_ray)

def sort_sensors(solvers, sensors):
    """
    TODO
    """
    #sort sensors
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
        merged_sensor = shdom.sensor.merge_sensor_rays(sensor_list)
        rte_sensors[key] = merged_sensor
        sensor_mapping[key] = mapping_list
    return rte_sensors, sensor_mapping

def get_measurements(solvers,sensors, n_jobs=1, mpi_comm=None, destructive=False, maxiter=100,verbose=True):
    """
    TODO
    """
    rte_sensors, sensor_mapping = sort_sensors(solvers, sensors)

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
