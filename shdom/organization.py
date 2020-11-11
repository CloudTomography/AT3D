
from collections import OrderedDict
import numpy as np
import xarray as xr
from joblib import Parallel, delayed
from shdom import core, gradient, medium

class SensorsDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        super(OrderedDict, self).__init__(*args, **kwargs)

    def add_sensor(self, key, value):
        if not key in self:
            self._add_instrument(key)
        self[key]['sensor_list'].append(value)

    def _add_instrument(self, key):
        self[key] = {'sensor_list': []}

    def make_forward_sensors(self):
        """
        TODO
        """
        forward_sensors = SensorsDict()
        for key,instrument in self.items():
            forward_instrument= OrderedDict()
            forward_instrument['sensor_list'] = [single_sensor.copy(deep=True) for
                                             single_sensor in instrument['sensor_list']]
            forward_sensors[key] = forward_instrument

        return forward_sensors

    @property
    def nmeasurements(self):
        nmeasurements = 0
        for instrument in self.values():
            for sensor in instrument['sensor_list']:
                nmeasurements += sensor.npixels.size*np.sum(sensor.stokes.data)
        return nmeasurements

    def sort_sensors(self, solvers, inverse=False):
        """
        TODO
        """
        rte_sensors = OrderedDict()
        sensor_mappings = OrderedDict()

        for key, solver in solvers.items():
            sensor_list = []
            mapping_list = []
            for instrument,instrument_data in self.items():
                for i,sensor in enumerate(instrument_data['sensor_list']):
                    if key == sensor.wavelength:
                        sensor_list.append(sensor)
                        mapping_list.append((instrument,i))

            var_list = ['ray_x','ray_y','ray_z','ray_mu','ray_phi','ray_weight','pixel_index']
            output = {}
            for var in var_list:
                concatenated = xr.concat([sensor[var] for sensor in sensor_list], dim='nrays')
                output[var] = ('nrays', concatenated)

            output['stokes'] = xr.concat([sensor.stokes for sensor in sensor_list],dim='nimage')
            output['rays_per_image'] = ('nimage', np.array([sensor.sizes['nrays'] for sensor in sensor_list]))
            output['rays_per_pixel'] = ('npixels', np.concatenate([np.unique(sensor.pixel_index,return_counts=True)[1]\
                                                                  for sensor in sensor_list]))
            if inverse:
                stokes_weights = []
                stokes_datas = []
                for sensor in sensor_list:
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

            if inverse:
                #TODO HARDCODED UNCERTAINTIES.
                #If we add an uncertainty model to sensors then this should be processed here.
                merged_sensors['uncertainties'] = xr.DataArray(data= np.ones((merged_sensors.sizes['nstokes'],
                                            merged_sensors.sizes['nstokes'], merged_sensors.sizes['npixels'])),
                                            dims=['nstokes', 'nstokes2', 'npixels'])
                    #NB Be aware that repeated dims cause errors. so the second dim is set to 'nstokes2' even though they are identical.
                    #Issue 1378 on xarray.
            rte_sensors[key] = merged_sensors
            sensor_mappings[key] = mapping_list

        return rte_sensors, sensor_mappings

    def get_unique_solvers(self):
        """
        scans the sensors and finds the unique solvers required
        and assigns each one a key. At the moment, each
        TODO
        """
        wavelength_list = []
        for key,instrument in self.items():
            for sensor in instrument['sensor_list']:
                wavelength_list.append(sensor.wavelength)
        return np.unique(wavelength_list)

    def get_minimum_stokes(self):
        """
        TODO
        """
        min_stokes = OrderedDict({key:0 for key in self.get_unique_solvers()})
        for key in min_stokes:
            for instrument in self.values():
                for sensor in instrument['sensor_list']:
                    if key == float(sensor.wavelength.data):
                        if np.all(sensor.stokes.data):
                            nstoke = 4
                        else:
                            nstoke = np.where(sensor.stokes.data == False)[0][0]
                            nstoke = nstoke if nstoke != 2 else 3
                        if min_stokes[float(sensor.wavelength.data)] < nstoke:
                            min_stokes[float(sensor.wavelength.data)] = nstoke
        return min_stokes

    def get_images(self, key):
        """
        TODO
        """
        return [self.get_image(key, index) for index in range(len(self[key]['sensor_list']))]


    def get_image(self, key, index):
        """
        TODO
        Make a 2D sensor for visualization purposes.
        """
        sensor = self[key]['sensor_list'][index]
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

    def add_measurements_inverse(self, sensor_mappings, measurements, measurement_keys):
        """
        TODO
        """

        for key in sensor_mappings:
            #group output from different (possibly parallel) workers by unique solver key (sensor_mappings.keys() == solvers.keys()).
            #if not parallel then measurement_keys == solvers.keys().
            indices = np.where(key == np.array(measurement_keys))[0]
            measurements_by_key = [measurements[i] for i in indices]

            variable_list_nray = [str(name) for name in measurements_by_key[0].data_vars if str(name) not in ('rays_per_image', 'stokes',
                                                                                           'rays_per_pixel', 'uncertainties',
                                                                                           'measurement_data','stokes_weights',
                                                                                           'I','Q','U','V')]
            variable_list_npixel = [str(name) for name in measurements_by_key[0].data_vars if str(name) in ('rays_per_pixel', 'uncertainties',
                                                                                           'measurement_data','stokes_weights',
                                                                                           'I','Q','U','V')]
            merged = {}
            for name in variable_list_nray:
                merged[name] = xr.concat([data[name] for data in measurements_by_key], dim='nrays')
            #process pixel variables for inverse problem.
            merged_pixel_variables = {}
            for name in variable_list_npixel:
                merged_pixel_variables[name] = xr.concat([data[name] for data in measurements_by_key], dim='npixels')
            merged.update(merged_pixel_variables)
            merged['stokes'] = measurements_by_key[0].stokes
            merged['rays_per_image'] = measurements_by_key[0].rays_per_image
            merged_measurements = xr.Dataset(merged)

            pixel_inds = np.cumsum(np.concatenate([np.array([0]), merged_measurements.rays_per_pixel.data])).astype(np.int)
            pixels_per_image = [np.where(pixel_inds == ray_image)[0][0] for ray_image in merged_measurements.rays_per_image.data]

            count = 0
            for i in range(merged_measurements.sizes['nimage']):
                split_index = pixels_per_image[i]
                sensor_measurements = merged_measurements.sel({'npixels': slice(count, count+split_index), 'nimage':i})
                forward_sensor = self[sensor_mappings[key][i][0]]['sensor_list'][sensor_mappings[key][i][1]]
                for stokes in forward_sensor.stokes_index:
                    if forward_sensor['stokes'].sel({'stokes_index':stokes}):
                        forward_sensor[str(stokes.data)] = sensor_measurements[str(stokes.data)]
                count += split_index


    def add_measurements_forward(self, sensor_mappings, measurements, measurement_keys):
        """
        TODO
        output_keys
        """
        for key in sensor_mappings:
            #group output from different (possibly parallel) workers by unique solver key (sensor_mappings.keys() == solvers.keys()).
            #if not parallel then measurement_keys == solvers.keys().
            indices = np.where(key == np.array(measurement_keys))[0]
            measurements_by_key = [measurements[i] for i in indices]
            # for i in indices:
            #     measurements_by_key.append(measurements[i])

            #Separately treat inverse and forward problem as inverse has different variables
            #and pixel variables have already been created.
            variable_list_nray = [str(name) for name in measurements_by_key[0].data_vars if str(name) not in ('rays_per_image', 'stokes',
                                                                                        'rays_per_pixel')]
            merged = {}
            for name in variable_list_nray:
                merged[name] = xr.concat([data[name] for data in measurements_by_key], dim='nrays')
            merged['stokes'] = measurements_by_key[0].stokes
            merged['rays_per_image'] = measurements_by_key[0].rays_per_image
            merged_measurements = xr.Dataset(merged)

            #split merged_measurements back into a list of individual images.
            #update each stored sensor with observables from each image. Use sensor mappings for this.
            count = 0
            for i in range(merged_measurements.sizes['nimage']):
                split_index = merged_measurements.rays_per_image[i].data
                sensor_measurements = merged_measurements.sel({'nrays': slice(count, count+split_index), 'nimage':i})
                #calculate observables.
                self._calculate_observables(sensor_mappings[key][i], sensor_measurements)
                count += split_index

    def _calculate_observables(self, mapping, rendered_rays):
        """
        TODO TESTS NEEDED.
        Currently averages over the sub_pixel rays and only takes
        the observables that are needed.
        """
        sensor = self[mapping[0]]['sensor_list'][mapping[1]]
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
            for stokes in temp.data_vars:
                sensor[stokes] = ('npixels', temp[stokes].data)
        else:
            observed_stokes = np.stack([var.data for name,var in temp.data_vars.items() if name != 'pixel_index'],axis=0)
            stokes_names = [name for name,var in temp.data_vars.items() if name != 'pixel_index']
            observables = core.average_subpixel_rays(pixel_index = temp.pixel_index.data, nstokes=observed_stokes.shape[0],
                                                weighted_stokes=observed_stokes,nrays=temp.nrays.size,
                                                npixels=sensor.npixels.size)
            for i, name in enumerate(stokes_names):
                sensor[name] = ('npixels', observables[i])

class SolversDict(OrderedDict):

    def __init__(self, *args, **kwargs):
        super(OrderedDict, self).__init__(*args, **kwargs)

    def add_solver(self, key, value):
        #TODO check that value is a valid solver and whether or not
        #the key already exists.
        self[key] = value

    def parallel_solve(self, n_jobs=1,mpi_comm=None,maxiter=100, verbose=True, init_solution=True,
                    setup_grid=True):
        """
        TODO
        """
        if mpi_comm is not None:
            for i in range(0, len(solvers), mpi_comm.Get_size()):
                index = i + mpi_comm.Get_rank()
                if index < len(solvers):
                    key = list(solvers)[index]
                    solvers[key].solve(maxiter=maxiter, verbose=verbose, init_solution=init_solution,setup_grid=setup_grid)

        if n_jobs==1:
            for solver in self.values():
                solver.solve(maxiter=maxiter, init_solution=init_solution, verbose=verbose,setup_grid=setup_grid)
        else:
            Parallel(n_jobs=n_jobs, backend="threading")(
                delayed(solver.solve, check_pickle=False)(maxiter=maxiter, init_solution=init_solution, verbose=verbose,setup_grid=setup_grid)
                for solver in self.values())

    def add_direct_beam_derivatives(self):
        """
        Calculate the geometry of the direct beam at each point and solver.
        Solver is modified in-place.
        TODO
        """
        #Only use the first solver under the assumption that all of them have the same
        #base grid.
        rte_solver = list(self.values())[0]
        rte_solver._direct_beam_derivative()

        for solver in self.values():
            solver._direct_derivative_path = rte_solver._direct_derivative_path
            solver._direct_derivative_ptr = rte_solver._direct_derivative_ptr

    def add_microphysical_partial_derivatives(self, table_to_grid_method, table_data):
        """
        TODO checks on inputs.
        """
        for key in self:
            self[key].calculate_microphysical_partial_derivatives(table_to_grid_method, table_data[key])

class UnknownScatterers(OrderedDict):
    def __init__(self, *args, **kwargs):
        super(OrderedDict, self).__init__(*args, **kwargs)

        def regular_grid(scatterer, table_data):
            return medium.table_to_grid(scatterer, table_data, inverse_mode=True)

        self._table_to_grid_method = regular_grid

    def add_unknown(self, scatterer_name, variable_name_list, table_data):
        if self._table_to_grid_method.__name__ == 'regular_grid':
            self[scatterer_name] = {'variable_name_list': np.atleast_1d(variable_name_list), 'table_data_input': table_data}
        else:
            raise NotImplementedError

    def create_derivative_tables(self):

        inputs = [(scatterer_name, self[scatterer_name]['table_data_input'],
                   self[scatterer_name]['variable_name_list'])
                 for scatterer_name in self]
        self.table_data = gradient.create_derivative_tables(inputs)

    @property
    def table_to_grid_method(self):
        return self._table_to_grid_method


def get_measurements(solvers, sensors, n_jobs=1, mpi_comm=None, maxiter=100,verbose=True, init_solution=True,
                    setup_grid=True, destructive=False):

    rte_sensors, sensor_mappings = sensors.sort_sensors(solvers)

    if mpi_comm is not None:
        out = []
        keys = []
        for i in range(0, len(solvers), mpi_comm.Get_size()):
            index = i + mpi_comm.Get_rank()
            if index < len(solvers):
                key = list(solvers)[index]
                solvers[key].solve(maxiter=maxiter, verbose=verbose, init_solution=init_solution,setup_grid=setup_grid)
                out.append(solvers[key].integrate_to_sensor(rte_sensors[key]))
                keys.append(key)
                if destructive: #memory management. After rendering the large arrays are released to ensure that the largest
                                #max_total_mb is not exceeded.
                    solvers[key]._release_big_arrays()

        out = mpi_comm.gather(out, root=0)
        keys = mpi_comm.gather(keys, root=0)
        if mpi_comm.Get_rank() == 0:
            out_flat = [item for sublist in out for item in sublist]
            keys_flat = [item for sublist in keys for item in sublist]
            organized_out = [out_flat[keys_flat.index(key)] for key in solvers]
            sensors.add_measurements_forward(sensor_mappings, organized_out, list(solvers))

    else:
        solvers.parallel_solve(n_jobs=n_jobs, mpi_comm=mpi_comm,maxiter=maxiter, verbose=verbose, init_solution=init_solution,
                            setup_grid=setup_grid)
        if n_jobs==1:
            out = [solvers[key].integrate_to_sensor(rte_sensors[key]) for key in solvers]
            keys = list(solvers)
        else:
            #decide on the division of n_jobs among solvers based on total number of pixels.
            #Note that the number of n_jobs here doesn't have to be the number of workers but can instead
            #be the number of subdivided job. This could be modified to ensure that all jobs are the correct size.
            #as is, this will make slightly more tasks than there are workers (n_jobs).
            keys, ray_start_end, pixel_start_end = subdivide_raytrace_jobs(rte_sensors, n_jobs)

            out= Parallel(n_jobs=n_jobs, backend='threading')(
                            delayed(solvers[key].integrate_to_sensor, check_pickle=False)(rte_sensors[key].sel(
                            nrays=slice(ray_start,ray_end),npixels=slice(pix_start, pix_end)))
                        for key, (ray_start,ray_end),(pix_start,pix_end) in zip(keys, ray_start_end, pixel_start_end))

        sensors.add_measurements_forward(sensor_mappings, out, keys)

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
