"""
Utility functions
"""
import numpy as np
from collections import OrderedDict
import netCDF4 as nc
import xarray as xr
from joblib import Parallel, delayed
import pandas as pd

import pyshdom.core
import pyshdom.gradient
import pyshdom.medium
import pyshdom.solver
import pyshdom.grid

class SensorsDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        super(OrderedDict, self).__init__(*args, **kwargs)

    def add_sensor(self, key, value):
        if not key in self:
            self._add_instrument(key)
        self[key]['sensor_list'].append(value)

    def _add_instrument(self, key):
        self[key] = {'sensor_list': [], 'uncertainty_model': None}

    def add_uncertainty_model(self, key, uncertainty_model):
        self[key]['uncertainty_model'] = uncertainty_model

    def add_noise(self, key, nstokes):
        for sensor in self[key]['sensor_list']:
            self[key]['uncertainty_model'].add_noise(sensor, nstokes)

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

    @property
    def npixels(self):
        npixels = 0
        for instrument in self.values():
            for sensor in instrument['sensor_list']:
                npixels += sensor.npixels.size
        return npixels

    def sort_sensors(self, solvers, measurements=None):
        """
        TODO check that measurements matches self.
        """
        rte_sensors = OrderedDict()
        sensor_mappings = OrderedDict()

        if measurements is not None:
            for key, solver in solvers.items():
                for instrument, instrument_data in measurements.items():
                    for i,sensor in enumerate(instrument_data['sensor_list']):
                        if key == sensor.wavelength:
                            if instrument_data['uncertainty_model'] is None:
                                uncertainty_data = np.ones((solver._nstokes, solver._nstokes, sensor.sizes['npixels']))
                                sensor['uncertainties'] = (['nstokes', 'nstokes2', 'npixels'], uncertainty_data)
                            else:
                                if 'uncertainties' not in sensor:
                                    instrument_data['uncertainty_model'].calculate_uncertainties(sensor, solver._nstokes)

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

            if measurements is not None:

                sensor_list_measure = []
                uncertainty_list_measure = []
                for instrument, instrument_data in measurements.items():
                    for sensor in instrument_data['sensor_list']:
                        if key == sensor.wavelength:
                            sensor_list_measure.append(sensor)
                            #no need for mapping as the sensors should match.

                stokes_weights = []
                stokes_datas = []
                for sensor in sensor_list_measure:
                    stokes_weight = np.zeros((solver._nstokes, sensor.sizes['npixels']))
                    stokes_data = np.zeros((solver._nstokes, sensor.sizes['npixels']))
                    for i,stokes in enumerate(sensor.stokes_index):
                        if stokes in list(sensor.data_vars):
                            stokes_weight[i,:] = 1.0
                            stokes_data[i,:] = sensor[str(stokes.data)]

                    stokes_weights.append(stokes_weight)
                    stokes_datas.append(stokes_data)
                output['uncertainties'] = xr.concat([sensor.uncertainties for sensor in sensor_list_measure], dim='npixels')
                output['stokes_weights'] = (['nstokes','npixels'], np.concatenate(stokes_weights,axis=-1))
                output['measurement_data'] = (['nstokes','npixels'], np.concatenate(stokes_datas,axis=-1))


            merged_sensors = xr.Dataset(data_vars=output)

            # if measurements is not None:
            #     #TODO HARDCODED UNCERTAINTIES.
            #     #If we add an uncertainty model to sensors then this should be processed here.
            #
            #
            #     error_variance = (np.repeat(np.concatenate(stokes_datas,axis=-1)[np.newaxis,...],merged_sensors.sizes['nstokes'],axis=0)*0.02)**2
            #     uncertainties = np.ones(error_variance.shape)/((0.1*0.001)**2) #error covariance for nothing just in case.
            #     #uncertainties[np.where(error_variance > 0.0)] = 1.0/error_variance[np.where(error_variance > 0.0)]
            #     #1.0/(np.repeat(np.concatenate(stokes_datas,axis=-1)[np.newaxis,...],merged_sensors.sizes['nstokes'],axis=0)*0.02)**2
            #     #uncertainties[np.where(~np.isfinite(uncertainties))] = 0.0
            #     merged_sensors['uncertainties'] = xr.DataArray(data= uncertainties,
            #                                                     dims=['nstokes', 'nstokes2', 'npixels'])
            #
            #     # np.ones((merged_sensors.sizes['nstokes'],
            #     #                             merged_sensors.sizes['nstokes'], merged_sensors.sizes['npixels'])),
            #     #                             dims=['nstokes', 'nstokes2', 'npixels'])
            #         #NB Be aware that repeated dims cause errors. so the second dim is set to 'nstokes2' even though they are identical.
            #         #Issue 1378 on xarray.
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
            observables = pyshdom.core.average_subpixel_rays(pixel_index = temp.pixel_index.data, nstokes=observed_stokes.shape[0],
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
                delayed(solver.solve)(maxiter=maxiter, init_solution=init_solution, verbose=verbose,setup_grid=setup_grid)
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
            return pyshdom.medium.table_to_grid(scatterer, table_data, inverse_mode=True)

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
        self.table_data = pyshdom.gradient.create_derivative_tables(inputs)

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
        if n_jobs==1 or n_jobs >= sensors.npixels:
            out = [solvers[key].integrate_to_sensor(rte_sensors[key]) for key in solvers]
            keys = list(solvers)
        else:
            #decide on the division of n_jobs among solvers based on total number of pixels.
            #Note that the number of n_jobs here doesn't have to be the number of workers but can instead
            #be the number of subdivided job. This could be modified to ensure that all jobs are the correct size.
            #as is, this will make slightly more tasks than there are workers (n_jobs).
            keys, ray_start_end, pixel_start_end = subdivide_raytrace_jobs(rte_sensors, n_jobs)

            out= Parallel(n_jobs=n_jobs, backend='threading')(
                            delayed(solvers[key].integrate_to_sensor)(rte_sensors[key].sel(
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

def float_round(x):
    """Round a float or np.float32 to a 3 digits float"""
    if type(x) == np.float32:
        x = x.item()
    return round(x,3)

def int_round(x):
    """Round a float or np.float32 to a 3 digits integer by 1000x scaling"""
    return int(np.round(x*1000))

def find_nearest(array, value):
    """Find the nearest element index in an array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def set_pyshdom_path():
    """set path to pyshdom parent directory"""
    import os
    from pathlib import Path
    os.chdir(str(Path(pyshdom.__path__[0]).parent))

def planck_function(temperature, wavelength, c=2.99792458e8,h=6.62606876e-34,k=1.3806503e-23):
    """
    temperature
        units, Kelvin
    wavelength
        units, micrometers
    radiance
        units, Watts/m^2/micrometer/steradian (SHDOM units)
    """
    wavelength = wavelength*1e-6
    radiance = 2*h*c**2/ wavelength**5 * 1.0 / (np.exp((h*c) / (wavelength*k*temperature)) - 1.0) * 1e-6
    return radiance

def cell_average_comparison(reference, other, variable_name):
    """
    calculates average values of 'variable name' in the cells
    of reference's grid for both reference and other (other is on a different grid.)
    """
    xgrid1 = reference.x.data
    ygrid1 = reference.y.data
    zgrid1 = reference.z.data
    values1 = reference[variable_name].data

    xgrid2 = other.x.data
    ygrid2 = other.y.data
    zgrid2 = other.z.data
    values2 = other[variable_name].data

    ref_vol,ref_val, other_vol, other_val = pyshdom.core.cell_average(xgrid1=xgrid1,ygrid1=ygrid1,zgrid1=zgrid1,
                        xgrid2=xgrid2,ygrid2=ygrid2,zgrid2=zgrid2,
                        values1=values1,
                        values2=values2)
    cell_average_ref = np.zeros(ref_vol.shape)
    cell_average_ref[np.where(ref_vol > 0.0)] = ref_val[np.where(ref_vol > 0.0)] /ref_vol[np.where(ref_vol > 0.0)]
    cell_average_other= np.zeros(other_vol.shape)
    cell_average_other[np.where(other_vol > 0.0)] = other_val[np.where(other_vol > 0.0)] /other_vol[np.where(other_vol > 0.0)]
    return cell_average_ref, cell_average_other

def load_2parameter_lwc_file(file_name, density='lwc'):
    """
    TODO
    Function that loads a scatterer from the '2 parameter lwc file' format used by
    SHDOM and i3rc monte carlo model.
    """
    header = pd.read_csv(file_name, nrows=4)
    nx, ny, nz = np.fromstring(header['2 parameter LWC file'][0], sep=' ').astype(np.int)
    dx, dy = np.fromstring(header['2 parameter LWC file'][1], sep=' ').astype(np.float)
    z = np.fromstring(header['2 parameter LWC file'][2], sep=' ').astype(np.float)
    temperature = np.fromstring(header['2 parameter LWC file'][3], sep=' ').astype(np.float)
    dset = pyshdom.grid.make_grid(dx, nx, dy, ny, z)

    data = np.genfromtxt(file_name, skip_header=5)

    lwc = np.zeros((nx, ny, nz))*np.nan
    reff = np.zeros((nx, ny, nz))*np.nan

    i, j, k = data[:, 0].astype(np.int)-1, data[:, 1].astype(np.int)-1, data[:, 2].astype(np.int)-1
    lwc[i, j, k] = data[:, 3]
    reff[i, j, k] = data[:, 4]

    dset['density'] = xr.DataArray(
        data=lwc,
        dims=['x', 'y', 'z']
    )

    dset['reff'] = xr.DataArray(
        data=reff,
        dims=['x', 'y', 'z']
    )

    dset['temperature'] = xr.DataArray(
        data=temperature,
        dims=['z']
    )

    dset.attrs['density_name'] = density
    dset.attrs['file_name'] = file_name

    return dset

def to_2parameter_lwc_file(file_name, cloud_scatterer, atmosphere=None, fill_temperature=280.0):
    """
    TODO
    Write lwc & reff to the '2 parameter lwc' file format used by i3rc MonteCarlo model and SHDOM.
    atmosphere should contain the temperature. It is interpolated to the specified z grid.
    If no atmosphere is included then a fill_temperature is used (Temperature is required
    in the file).
    """

    nx, ny, nz = cloud_scatterer.density.shape
    dx, dy = (cloud_scatterer.x[1]-cloud_scatterer.x[0]).data, (cloud_scatterer.y[1] - cloud_scatterer.y[0]).data
    z = cloud_scatterer.z.data

    if atmosphere is not None:
        temperature = atmosphere.interp({'z': cloud_scatterer.z}).temperature.data
    else:
        temperature = np.ones(z.shape)*fill_temperature

    i, j, k = np.meshgrid(np.arange(1, nx+1), np.arange(1, ny+1), np.arange(1, nz+1), indexing='ij')

    lwc = cloud_scatterer.density.data.ravel()
    reff = cloud_scatterer.reff.data.ravel()

    z_string = ''
    for z_value in z:
        if z_value == z[-1]:
            z_string += '{}'.format(z_value)
        else:
            z_string += '{} '.format(z_value)

    t_string = ''
    for i, temp_value in enumerate(temperature):
        if i == len(temperature) - 1:
            t_string += '{:5.2f}'.format(temp_value)
        else:
            t_string += '{:5.2f} '.format(temp_value)

    with open(file_name, "w") as f:
        f.write('2 parameter LWC file\n')
        f.write(' {} {} {}\n'.format(nx, ny, nz))
        f.write('{} {}\n'.format(dx, dy))
        f.write('{}\n'.format(z_string))
        f.write('{}\n'.format(t_string))
        for x, y, z, l, r in zip(i.ravel(), j.ravel(), k.ravel(), lwc.ravel(), reff.ravel()):
            f.write('{} {} {} {:5.4f} {:3.2f}\n'.format(x, y, z, l, r))

def load_from_csv(path, density=None, origin=(0.0,0.0)):
    """
    TODO
    """
    df = pd.read_csv(path, comment='#', skiprows=4, index_col=['x', 'y', 'z'])
    nx, ny, nz = np.genfromtxt(path, max_rows=1, dtype=int, delimiter=',')
    dx, dy = np.genfromtxt(path, max_rows=1, dtype=float, skip_header=2, delimiter=',')
    z = xr.DataArray(np.genfromtxt(path, max_rows=1, dtype=float, skip_header=3, delimiter=','), coords=[range(nz)], dims=['z'])

    dset = pyshdom.grid.make_grid(dx, nx, dy, ny, z)
    i, j, k = zip(*df.index)

    for name in df.columns:
        #initialize with np.nans so that empty data is np.nan
        variable_data = np.zeros((dset.sizes['x'], dset.sizes['y'], dset.sizes['z']))*np.nan
        variable_data[i, j, k] = df[name]
        dset[name] = (['x', 'y', 'z'], variable_data)

    if density is not None:
        assert density in dset.data_vars, \
        "density variable: '{}' must be in the file".format(density)

        dset = dset.rename_vars({density: 'density'})
        dset.attrs['density_name'] = density

    dset.attrs['file_name'] = path

    return dset

def load_from_netcdf(path, density=None):
    """
        TODO
    A shallow wrapper around open_dataset that sets the density_name.
    """
    dset = xr.open_dataset(path)

    if density is not None:
        if density not in dset.data_vars:
            raise ValueError("density variable: '{}' must be in the file".format(density))
        dset = dset.rename_vars({density: 'density'})
        dset.attrs['density_name'] = density

    dset.attrs['file_name'] = path

    return dset


def load_forward_model(file_name):
    """
    TODO
    """
    dataset = nc.Dataset(file_name)

    groups = dataset.groups
    sensors = groups['sensors'].groups
    solvers = groups['solvers'].groups
    sensor_dict = pyshdom.util.SensorsDict()
    solver_dict = pyshdom.util.SolversDict()

    for key,sensor in sensors.items():
        sensor_list = []
        for i, image in sensor.groups.items():

            sensor_dict.add_sensor(key, xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                'sensors/'+str(key)+'/'+str(i)])))
        #     sensor_list.append(xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
        #         'sensors/'+str(key)+'/'+str(i)])))
        # sensor_dict[key] = {'sensor_list':sensor_list}

    for key, solver in solvers.items():
        numerical_params = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                        'solvers/'+str(key)+'/numerical_parameters']))
        num_stokes = numerical_params.num_stokes.data
        surface = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/surface']))
        source = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/source']))

        mediums = OrderedDict()
        for name, med in solver['medium'].groups.items():
            mediums[name] = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                        'solvers/'+str(key)+'/medium/'+str(name)]))

        if 'atmosphere' in solver.groups.keys():
            atmosphere = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset[
                'solvers/'+str(key)+'/atmosphere']))
        else:
            atmosphere=None

        solver_dict.add_solver(float(key), pyshdom.solver.RTE(numerical_params=numerical_params,
                                            medium=mediums,
                                           source=source,
                                           surface=surface,
                                            num_stokes=num_stokes,
                                            name=None
                                           )
                                           )
        rte_grid = xr.open_dataset(xr.backends.NetCDF4DataStore(dataset['solvers/'+str(key)+'/grid']))

    return sensor_dict, solver_dict, rte_grid
#TODO add checks here for if file exists etc.
def save_forward_model(file_name,sensors, solvers):
    """
    TODO
    """
    for i, (key,sensor) in enumerate(sensors.items()):
        for j, image in enumerate(sensor['sensor_list']):
            if (i==0) & (j==0):
                image.to_netcdf(file_name, 'w', group = 'sensors/'+str(key)+'/'+str(j), format='NETCDF4',engine='netcdf4')
            else:
                image.to_netcdf(file_name, 'a', group = 'sensors/'+str(key)+'/'+str(j), format='NETCDF4',engine='netcdf4')

    for i, (key, solver) in enumerate(solvers.items()):

        numerical_params = solver.numerical_params
        numerical_params['num_stokes'] = solver._nstokes

        solver.numerical_params.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'numerical_parameters', format='NETCDF4',engine='netcdf4')
        solver.surface.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'surface', format='NETCDF4',engine='netcdf4')
        solver.source.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'source', format='NETCDF4',engine='netcdf4')
        solver._grid.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'grid', format='NETCDF4',engine='netcdf4')
        if solver.atmosphere is not None:
            solver.atmosphere.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'atmosphere', format='NETCDF4',engine='netcdf4')
        for name, med in solver.medium.items():
            med.to_netcdf(file_name,'a', group='solvers/'+str(key)+'/'+'medium/'+str(name), format='NETCDF4',engine='netcdf4')
