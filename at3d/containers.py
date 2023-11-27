"""
This module defines containers for groups of RTE solvers and sensors
and scattering tables to streamline retrievals or the generation of
synthetic measurements in scripts.

The three objects defined here are `SensorsDict`, `SolversDict` and
`UnknownScatterers` and they are all wrappers around OrderedDict. Data
are not stored as attributes meaning they can be modified without
utilizing the additional 'add_X' methods.

`SensorsDict` holds measurements, either real or synthetic. The `SensorsDict`
groups individual sensor xr.Datasets (see sensor.py) by their instrument
which has a common uncertainty model. The `SensorsDict` object contains the methods
for reorganizing measurements from this sensor-centric organization to match
them with the appropriate solver.RTE objects.

The solver.RTE objects are stored in `SolversDict`. This object contains method
for solving a set of solver.RTE objects in parallel or calling a couple of methods
on all solvers in the `SolversDict.`

`UnknownScatterers` is only utilized during the optimization and it stores the
data and methods for calculating the partial derivatives of optical properties
with respect to the unknowns.
"""

from collections import OrderedDict
from joblib import Parallel, delayed

import numpy as np
import xarray as xr
import warnings
import typing

import at3d.core
import at3d.gradient
import at3d.medium
import at3d.parallel
import at3d.checks

class SensorsDict(OrderedDict):
    """
    Stores measurement geometry and has methods for rearranging the measurement
    geometry to send to solver.RTE for calculation and then postprocessing the
    measurements and stores them.

    Sensor data (see sensor.py) is stored by instrument, each of which has a
    list of sensors and an accompanying uncertainty model.

    """
    def add_sensor(self, instrument, sensor):
        """
        Adds a sensor Dataset to a given instrument's sensor list.

        Parameters
        ----------
        instrument : Any Type.
            The name of the instrument for this sensor.
        sensor : xr.Dataset
            A valid sensor.
        """
        if not instrument in self:
            self._add_instrument(instrument)
        at3d.checks.check_sensor(sensor)
        self[instrument]['sensor_list'].append(sensor)

    def _add_instrument(self, key):
        """Initializes the key for an instrument"""
        self[key] = {'sensor_list': [], 'uncertainty_model': None}

    def add_uncertainty_model(self, instrument, uncertainty_model):
        """Updates the uncertainty model for a given instrument.

        Paramaters
        ----------
        instrument : Any Type
            The key for the instrument.
        uncertainty_model : callable
            The model that will be called on an xr.Dataset to calculate
            the error-covariance matrix from the Stokes Vector.
        """
        if not instrument in self:
            self._add_instrument(instrument)
        self[instrument]['uncertainty_model'] = uncertainty_model

    def calculate_uncertainties(self, instrument):
        """Evaluate the uncertainty matrix for each sensor using
        the instrument's uncertainty model.

        Parameters
        ----------
        instrument : Any Type
            The key for the instrument.
        """
        if 'uncertainty_model' not in self[instrument]:
            raise KeyError("There is no uncertainty model for instrument '{}'."
                " please add one.")
        for sensor in self[instrument]['sensor_list']:
            self[instrument]['uncertainty_model'].calculate_uncertainties(sensor)

    def add_noise(self, instrument):
        """Add noise to the observables in the sensors for the
        specified instrument according to the prescribed noise model.

        Parameters
        ----------
        instrument : Any Type
            The key for the instrument.
        """
        for sensor in self[instrument]['sensor_list']:
            self[instrument]['uncertainty_model'].add_noise(sensor)

    def make_forward_sensors(self, instrument_list=None):
        """Make a deep copy of self.

        Used when wanting to exactly replicate the sensor geometry to store the
        output of the forward model during optimization.
        """
        forward_sensors = SensorsDict()
        if instrument_list is None:
            instrument_list = self
        for key in instrument_list:
            if not key in self:
                raise KeyError("Instrument '{}' is not in SensorsDict")
            instrument = self[key]
        # for key, instrument in self.items():
            forward_instrument = OrderedDict()
            forward_instrument['sensor_list'] = [single_sensor.copy(deep=True) for
                                                 single_sensor in instrument['sensor_list']]
            forward_instrument['uncertainty_model'] = instrument['uncertainty_model']
            forward_sensors[key] = forward_instrument

        return forward_sensors

    def get_measurements(self, solvers, n_jobs=1, mpi_comm=None, maxiter=100, verbose=True,
                         init_solution=True, setup_grid=True, destructive=False,
                         overwrite_solver=False):
        """
        Calculates the observables specified in the sensors in self using
        RTE solvers.

        If necessary, the SHDOM solutions are performed and the StokesVector at
        each specified ray is calculated. Rays are then averaged to 'pixel'
        observables that are then stored matching the sensor they came from in
        self.

        Parameters
        ----------
        solvers : at3d.containers.SolversDict
            The solvers that will be used to calculate the measurements.
        n_jobs : 1
            The number of workers for shared memory parallelization of the RTE
            solutions and the
        mpi_comm : mpi4py.MPI.Intracomm
            An MPI communicator for the MPI parallelization of different RTE
            solutions. MPI is prioritized over other parallelization methods.
        destructive : bool
            If True when MPI is used then the large memory is released after measurements
            are simulated. This means that solvers will have to be resolved to
            recalculate measurements but can prevent running out of memory if a large
            number of sequential RTE solutions are calculated by each worker.
        overwrite_solver : bool
            If True then forces all solvers to resolve the RTE even if they have
            already converged.
        maxiter: integer
            Maximum number of iterations for the iterative solution to SHDOM.
        setup_grid: boolean
            If True then a new grid is initialized. If False then the grid
            (including adaptive grid points)
        init_solution: boolean, default=True
            If False then a solution is initialized. This is overwritten
            to True if no existing Radiance/Source function fields exist.
            The solution initialization will use a Radiance/Source field provided
            by RTE.load_solution or use a 1D two-stream model.
        verbose: boolean
            True will output solution iteration information into stdout.
        """
        if not isinstance(solvers, SolversDict):
            raise TypeError(
                "`solvers` should be of type '{}' not '{}'".format(
                    SolversDict, type(solvers))
                )
        if not isinstance(destructive, bool):
            raise TypeError("`destructive` should be a boolean.")

        rte_sensors, sensor_mappings = self.sort_sensors(solvers)
        keys, to_solve = solvers.to_solve(overwrite_solver)

        if mpi_comm is not None:
            out = []
            keys = []
            for i in range(0, len(to_solve), mpi_comm.Get_size()):
                index = i + mpi_comm.Get_rank()
                if index < len(to_solve):
                    key = keys[index]
                    to_solve[index].solve(
                        maxiter=maxiter, verbose=verbose, init_solution=init_solution,
                        setup_grid=setup_grid)
                    out.append(solvers[key].integrate_to_sensor(rte_sensors[key]))
                    keys.append(key)
                    if destructive:
                        # memory management. After rendering, the large arrays are
                        # released to ensure that the largest
                        # max_total_mb is not exceeded.
                        solvers[key]._release_big_arrays()

            out = mpi_comm.gather(out, root=0)
            keys = mpi_comm.gather(keys, root=0)
            if mpi_comm.Get_rank() == 0:
                out_flat = [item for sublist in out for item in sublist]
                keys_flat = [item for sublist in keys for item in sublist]
                organized_out = [out_flat[keys_flat.index(key)] for key in solvers]
                self.add_measurements_forward(sensor_mappings, organized_out, list(solvers))

        else:
            solvers.solve(n_jobs=n_jobs, mpi_comm=mpi_comm, maxiter=maxiter,
                                   verbose=verbose, init_solution=init_solution,
                                   setup_grid=setup_grid, overwrite_solver=overwrite_solver)
            if n_jobs == 1 or n_jobs >= self.npixels:
                out = [solvers[key].integrate_to_sensor(rte_sensors[key]) for key in solvers]
                keys = list(solvers)
            else:
                #decide on the division of n_jobs among solvers based on total number of pixels.
                #Note that the number of n_jobs here doesn't have to be the number of workers but can instead
                #be the number of subdivided job. This could be modified to ensure that all jobs are the correct size.
                #as is, this will make slightly more tasks than there are workers (n_jobs).
                keys, ray_start_end, pixel_start_end = at3d.parallel.subdivide_raytrace_jobs(rte_sensors, n_jobs)

                out = Parallel(n_jobs=n_jobs, backend='threading')(
                    delayed(solvers[key].integrate_to_sensor)(rte_sensors[key].sel(
                        nrays=slice(ray_start, ray_end), npixels=slice(pix_start, pix_end)))
                    for key, (ray_start, ray_end), (pix_start, pix_end) in
                    zip(keys, ray_start_end, pixel_start_end))

            self.add_measurements_forward(sensor_mappings, out, keys)

    def sort_sensors(self, solvers, measurements=None):
        """Groups sensors by RTE solver for evaluation of observables.
        Also prepares measurements for cost/gradient calculation in inverse
        problem.

        `solvers` contains the solver.RTE objects that will be used to
        calculate the observables. This function groups all sensors in
        `self` by the keys of `solvers` and the inverse mapping.
        If `measurements` is not None then the inverse error-covariance
        is evaluated in preparation for evaluation of the cost function.

        Parameters
        ----------
        solvers : at3d.containers.SolversDict
            The SolversDict containing the solver.RTE objects that will be
            used to evaluate the observables.
        measurements : at3d.containers.SensorsDict
            This contains the actual measurement data that is used to constrain
            the retrieval.

        Returns
        -------
        rte_sensors : OrderedDict
            Keys match `solvers` and each entry is an xr.Dataset containing the
            concatenated sensors from `self`. Also contains the start/end
            pixels/rays for mapping back to a list of xr.Datasets.
        sensor_mappings : OrderedDict
            Keys are from `solvers` and the entries are the instrument name and
            index for mapping each entry in the sensor list back into `self`.

        Raises
        ------
        TypeError
            If `solvers` is not at3d.containers.SolversDict or `measurements` is
            not None or at3d.containers.SensorsDict.
        """
        if not isinstance(solvers, SolversDict):
            raise TypeError("`solvers` should be of type '{}' not '{}'".format(SolversDict, type(solvers)))
        if measurements is not None:
            if not isinstance(measurements, SensorsDict):
                raise TypeError("`measurements` should be of type '{}' not '{}'".format(SensorsDict, type(measurements)))

        rte_sensors = OrderedDict()
        sensor_mappings = OrderedDict()

#         if measurements is not None:
#             for key, solver in solvers.items():
#                 for instrument, instrument_data in measurements.items():
#                     for i, sensor in enumerate(instrument_data['sensor_list']):
#                         if key == sensor.wavelength:
#                             if instrument_data['uncertainty_model'] is None:
#                                 uncertainty_data = np.ones((solver._nstokes,
#                                                             solver._nstokes,
#                                                             sensor.sizes['npixels']))
#
#                                 sensor['uncertainties'] = (['nstokes', 'nstokes2', 'npixels'],
#                                                            uncertainty_data)
#                             else:
#                                 if 'uncertainties' not in sensor:
#                                     instrument_data['uncertainty_model'].calculate_uncertainties(
#                                         sensor, solver._nstokes)
        var_list = ['ray_x', 'ray_y', 'ray_z', 'ray_mu', 'ray_phi', 'ray_weight', 'pixel_index']
        for key, solver in solvers.items():
            sensor_list = []
            mapping_list = []
            for instrument, instrument_data in self.items():
                for i, sensor in enumerate(instrument_data['sensor_list']):
                    if key == sensor.wavelength:
                        sensor_list.append(sensor)
                        mapping_list.append((instrument, i))
            output = {}
            if not sensor_list:
                warnings.warn("No sensors found matching solver with key '{}'".format(key))
            else:
                for var in var_list:
                    concatenated = xr.concat([sensor[var] for sensor in sensor_list], dim='nrays')
                    output[var] = concatenated

                output['stokes'] = xr.concat([sensor.stokes for sensor in sensor_list], dim='nimage')
                output['rays_per_image'] = ('nimage', np.array([sensor.sizes['nrays']
                                                                for sensor in sensor_list]))
                output['rays_per_pixel'] = ('npixels', np.concatenate([np.unique(sensor.pixel_index,
                                                                                 return_counts=True)[1]
                                                                       for sensor in sensor_list]))

                if measurements is not None:

                    sensor_list_measure = []
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
                        for i, stokes in enumerate(sensor.stokes_index):
                            if stokes in list(sensor.data_vars):
                                stokes_weight[i, :] = 1.0
                                stokes_data[i, :] = sensor[str(stokes.data)]

                        stokes_weights.append(stokes_weight)
                        stokes_datas.append(stokes_data)
                    output['uncertainties'] = xr.concat([sensor.uncertainties
                                                         for sensor in sensor_list_measure],
                                                        dim='npixels')
                    output['stokes_weights'] = (['nstokes', 'npixels'],
                                                np.concatenate(stokes_weights, axis=-1))
                    output['measurement_data'] = (['nstokes', 'npixels'],
                                                  np.concatenate(stokes_datas, axis=-1))

            merged_sensors = xr.Dataset(data_vars=output)
            rte_sensors[key] = merged_sensors
            sensor_mappings[key] = mapping_list

        return rte_sensors, sensor_mappings

    def get_unique_solvers(self):
        """
        Finds the set of unique wavelengths among all sensors.
        """
        wavelength_list = []
        for instrument in self.values():
            for sensor in instrument['sensor_list']:
                wavelength_list.append(sensor.wavelength)
        return np.unique(wavelength_list)

    def get_minimum_stokes(self):
        """
        Find the minimum required number of stokes parameters
        for each solver to simulate the required observables.

        Notes
        -----
        Users may want to not use this if they want more Stokes
        components for accuracy reasons.
        """
        min_stokes = OrderedDict({key:0 for key in self.get_unique_solvers()})
        for key in min_stokes:
            for instrument in self.values():
                for sensor in instrument['sensor_list']:
                    if key == float(sensor.wavelength.data):
                        if np.all(sensor.stokes.data):
                            nstoke = 4
                        else:
                            nstoke = np.where(~sensor.stokes.data)[0][0]
                            nstoke = nstoke if nstoke != 2 else 3
                        if min_stokes[float(sensor.wavelength.data)] < nstoke:
                            min_stokes[float(sensor.wavelength.data)] = nstoke
        return min_stokes

    def get_images(self, instrument):
        """Returns a list of images for a specified instrument.

        Parameters
        ----------
        instrument : Any
            The key in `self` corresponding to the given instrument. Typically
            a string.

        Returns
        -------
        image_list : List
            A list of xr.Datasets which are the xr.Datasets in `instrument`'s
            sensor_list rearranged to form 2D images if applicable.

        Raises
        ------
        ValueError
            If any of the xr.Datasets in `self`[`instrument`] does not have a
            valid image shape.
        """
        image_list = [self.get_image(instrument, index) for index in range(len(self[instrument]['sensor_list']))]
        return image_list


    def get_image(self, instrument, sensor_index):
        """Make a 2D sensor for visualization purposes.

        Parameters
        ----------
        instrument : Any
            The key corresponding to the desired instrument.
        sensor_index : int
            The index of the sensor in `instrument`'s sensor_list.

        Returns
        -------
        img_data : xr.Dataset
            Contains the given sensor dataset reshaped into a 2D image.

        Raises
        ------
        ValueError
            If the specified sensor doesn't have a 2D image.
        KeyError
            If the observables (e.g. I, Q, U, V) have not been calculated
            for the sensor yet.
        """
        sensor = self[instrument]['sensor_list'][sensor_index]
        if not hasattr(sensor, 'image_shape'):
            raise ValueError("Sensor dataset does not have an "
                             "'image_shape' variable."
                             " A 2D image cannot be formed.")
        img_data = xr.Dataset(
            data_vars={
                'x': (['imgdim0', 'imgdim1'],
                      sensor.cam_x.data.reshape(sensor.image_shape.data, order='F')),
                'y': (['imgdim0', 'imgdim1'],
                      sensor.cam_y.data.reshape(sensor.image_shape.data, order='F')),
                'mu': (['imgdim0', 'imgdim1'],
                       sensor.cam_mu.data.reshape(sensor.image_shape.data, order='F')),
                'phi': (['imgdim0', 'imgdim1'],
                        sensor.cam_phi.data.reshape(sensor.image_shape.data, order='F'))
            }
        )
        for stokes in sensor.stokes_index:
            if sensor.stokes.sel({'stokes_index': stokes}):
                img_data[str(stokes.data)] = (['imgdim0', 'imgdim1'],
                                              sensor[str(stokes.data)].data.reshape(
                                                  sensor.image_shape.data, order='F'))
        return img_data

    def add_measurements_inverse(self, sensor_mappings, measurements, measurement_keys):
        """Takes the output of a (possibly parallel) rendering of synthetic measurements
        by at3d.gradient.gradient_l2 during the optimization and updates
        the sensors in `self` with the observables (e.g. I, Q, etc).

        Parameters
        ----------
        sensor_mappings : OrderedDict
            Keys are from the SolversDict used to calculate the
            synthetic measurements. Values are are list of tuples of
            the instrument names and sensor indices.
        measurements : List
            A list of xr.Datasets containing the Stokes observables calculated
            by at3d.gradient.gradient_l2.
        measurement_keys : List
            The solver key for each element in `measurements`. These keys
            may be duplicated as `measurements` may have been calculated
            in parallel.

        Notes
        -----
        This method is similar to SensorsDict.add_measurements_forward but as the output
        of evaluating at3d.gradient.gradient_l2 is different to
        at3d.solver.RTE.integrate_to_sensor there are differences. Particularly,
        pixel averaged quantities are output by at3d.gradient.gradient_l2
        while ray quantities are output by input to self.add_measurements_forward.
        See sensor.py for the difference between ray and pixel variables.
        """

        for key in sensor_mappings:
            #group output from different (possibly parallel) workers by
            #unique solver key (sensor_mappings.keys() == solvers.keys()).
            #if not parallel then measurement_keys == solvers.keys().
            indices = np.where(key == np.array(measurement_keys))[0]
            measurements_by_key = [measurements[i] for i in indices]

            variable_list_nray = [str(name) for name in measurements_by_key[0].data_vars
                                  if str(name) not in ('rays_per_image', 'stokes',
                                                       'rays_per_pixel', 'uncertainties',
                                                       'measurement_data', 'stokes_weights',
                                                       'I', 'Q', 'U', 'V')
                                  ]
            variable_list_npixel = [str(name) for name in measurements_by_key[0].data_vars
                                    if str(name) in ('rays_per_pixel', 'uncertainties',
                                                     'measurement_data', 'stokes_weights',
                                                     'I', 'Q', 'U', 'V')
                                    ]
            merged = {}
            for name in variable_list_nray:
                merged[name] = xr.concat([data[name] for data in measurements_by_key], dim='nrays')
            #process pixel variables for inverse problem.
            merged_pixel_variables = {}
            for name in variable_list_npixel:
                merged_pixel_variables[name] = xr.concat([data[name]
                                                          for data in measurements_by_key],
                                                         dim='npixels')
            merged.update(merged_pixel_variables)
            merged['stokes'] = measurements_by_key[0].stokes
            merged['rays_per_image'] = measurements_by_key[0].rays_per_image
            merged_measurements = xr.Dataset(merged)

            pixel_inds = np.cumsum(np.concatenate([np.array([0]),
                                                   merged_measurements.rays_per_pixel.data
                                                   ])).astype(int)
            pixels_per_image = [np.where(pixel_inds == ray_image)[0][0]
                                for ray_image in merged_measurements.rays_per_image.data]

            count = 0
            for i in range(merged_measurements.sizes['nimage']):
                split_index = pixels_per_image[i]
                sensor_measurements = merged_measurements.sel({
                    'npixels': slice(count, count+split_index),
                    'nimage':i})
                forward_sensor = self[sensor_mappings[key][i][0]]['sensor_list'][sensor_mappings[key][i][1]]
                for stokes in forward_sensor.stokes_index:
                    if forward_sensor['stokes'].sel({'stokes_index':stokes}):
                        forward_sensor[str(stokes.data)] = sensor_measurements[str(stokes.data)]
                count += split_index


    def add_measurements_forward(self, sensor_mappings, measurements, measurement_keys):
        """Takes the output of a (possibly parallel) rendering of synthetic measurements
        by at3d.solver.integrate_to_sensor and updates the sensors in `self`
        with the observables (e.g. pixel-averaged I, Q, etc).

        Parameters
        ----------
        sensor_mappings : OrderedDict
            Keys are from the SolversDict used to calculate the
            synthetic measurements. Values are are list of tuples of
            the instrument names and sensor indices.
        measurements : List
            A list of xr.Datasets containing the Stokes observables calculated
            by at3d.solver.integrate_to_sensor
        measurement_keys : List
            The solver key for each element in `measurements`. These keys
            may be duplicated as `measurements` may have been calculated
            in parallel.

        Notes
        -----
        Calls self._calculate_observables to evaluate the pixel-averaged observables.
        See also SensorsDict.add_measurements_inverse and sensor.py for definition
        of ray and pixel quantities.
        """
        for key in sensor_mappings:
            #group output from different (possibly parallel) workers by unique solver key (sensor_mappings.keys() == solvers.keys()).
            #if not parallel then measurement_keys == solvers.keys().
            indices = np.where(key == np.array(measurement_keys))[0]
            measurements_by_key = [measurements[i] for i in indices]

            #Separately treat inverse and forward problem as inverse has different variables
            #and pixel variables have already been created.
            variable_list_nray = [str(name) for name in measurements_by_key[0].data_vars
                                  if str(name) not in ('rays_per_image',
                                                       'stokes',
                                                       'rays_per_pixel')
                                  ]
            merged = {}
            for name in variable_list_nray:
                merged[name] = xr.concat([data[name]
                                          for data in measurements_by_key],
                                         dim='nrays')
            merged['stokes'] = measurements_by_key[0].stokes
            merged['rays_per_image'] = measurements_by_key[0].rays_per_image
            merged_measurements = xr.Dataset(merged)

            #split merged_measurements back into a list of individual images.
            #update each stored sensor with observables from each image. Use sensor mappings for this.
            count = 0
            for i in range(merged_measurements.sizes['nimage']):
                split_index = merged_measurements.rays_per_image[i].data
                sensor_measurements = merged_measurements.sel({
                    'nrays': slice(count, count+split_index),
                    'nimage':i})
                #calculate observables.
                self._calculate_observables(sensor_mappings[key][i], sensor_measurements)
                count += split_index

    def _calculate_observables(self, mapping, rendered_rays):
        """
        Calculates the observables (pixel averaged Stokes components)
        required by the sensor using the ray Stokes variables output
        from at3d.solver.RTE.integrate_to_sensor.

        Parameters
        ----------
        mapping : Tuple
            First component is the instrument key in `self` and
            second component is the index of the sensor in that
            instrument's sensor_list.
        rendered_rays : xr.Dataset
            Contains the ray variables including Stokes components
            calculated by at3d.solver.RTE.integrate_to_sensor.

        Notes
        -----
        Calls at3d.core.average_subpixel_rays to do the averaging.
        See src/util.f90.
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

        merged = xr.merge(merge_list)

        if not sensor.use_subpixel_rays:
            for stokes in merged.data_vars:
                sensor[stokes] = ('npixels', merged[stokes].data)
        else:
            observed_stokes = np.stack([var.data for name, var in
                                        merged.data_vars.items()
                                        if name != 'pixel_index'], axis=0)
            stokes_names = [name for name, var in merged.data_vars.items() if name != 'pixel_index']
            observables = at3d.core.average_subpixel_rays(
                pixel_index=merged.pixel_index.data,
                nstokes=observed_stokes.shape[0],
                weighted_stokes=observed_stokes,
                nrays=merged.nrays.size,
                npixels=sensor.npixels.size)
            for i, name in enumerate(stokes_names):
                sensor[name] = ('npixels', observables[i])

    @property
    def nmeasurements(self):
        """Calculates the total number of observed quantities
        (pixel-averaged Stokes components) in the sensors in `self`.
        """
        nmeasurements = 0
        for instrument in self.values():
            for sensor in instrument['sensor_list']:
                nmeasurements += sensor.npixels.size*np.sum(sensor.stokes.data)
        return nmeasurements

    @property
    def npixels(self):
        """
        Calculates the total number of pixels in the sensors in `self`.
        """
        npixels = 0
        for instrument in self.values():
            for sensor in instrument['sensor_list']:
                npixels += sensor.npixels.size
        return npixels


class SolversDict(OrderedDict):
    """
    Stores multiple solver.RTE objects and has methods for solving in parallel
    as well as pre-processing for the evaluation of the cost/gradient.
    """
    def add_solver(self, key, solver):
        """Adds a at3d.solver.RTE object to self.

        Parameters
        ----------
        key : Any
            The key used to uniquely identify the solver which is typically
            the monochromatic wavelength as a float.
        solver : at3d.solver.RTE
            The solver ojbect to add to `self`.

        Raises
        ------
        TypeError
            If `solver` is not of the type at3d.solver.RTE
        """
        if not isinstance(solver, at3d.solver.RTE):
            raise TypeError("solver should be of type '{}'".format(at3d.solver.RTE))
        self[key] = solver

    def solve(self, n_jobs=1, mpi_comm=None, overwrite_solver=False, maxiter=100,
                       verbose=True, init_solution=True, setup_grid=True):
        """
        Solves in parallel all solver.RTE objects using MPI or multi-threading.

        Parameters
        ----------
        n_jobs : int
            The number of workers if using multi-threaded parallelization.
        mpi_comm : mpi4py.MPI.Intracomm
            The MPI communicator forces the MPI parallelization to run.
            Note that there is no type check on mpi_comm.
        overwrite_solver : bool
            If True then forces all solvers to resolve the RTE even if they have
            already converged.
        maxiter: integer
            Maximum number of iterations for the iterative solution to SHDOM.
        setup_grid: boolean
            If True then a new grid is initialized. If False then the grid
            (including adaptive grid points)
        init_solution: boolean, default=True
            If False then a solution is initialized. This is overwritten
            to True if no existing Radiance/Source function fields exist.
            The solution initialization will use a Radiance/Source field provided
            by RTE.load_solution or use a 1D two-stream model.
        verbose: boolean
            True will output solution iteration information into stdout.
        """
        key_list, to_solve = self.to_solve(overwrite_solver)
        if mpi_comm is not None:
            for i in range(0, len(to_solve), mpi_comm.Get_size()):
                index = i + mpi_comm.Get_rank()
                if index < len(to_solve):
                    to_solve[index].solve(
                        maxiter=maxiter, verbose=verbose, init_solution=init_solution,
                        setup_grid=setup_grid)
        else:
            if n_jobs == 1:
                for solver in to_solve:
                    solver.solve(maxiter=maxiter, init_solution=init_solution,
                                 verbose=verbose, setup_grid=setup_grid)
            else:
                Parallel(n_jobs=n_jobs, backend="threading")(
                    delayed(solver.solve)(maxiter=maxiter, init_solution=init_solution,
                                          verbose=verbose, setup_grid=setup_grid)
                    for solver in to_solve)

    def to_solve(self, overwrite_solver):
        """
        Returns the keys and solvers of solvers that have not already
        converged to the prescribed solution accuracy.

        This is used to avoid repeating possibly expensive RTE solutions
        if not necessary.

        Parameters
        ----------
        overwrite_solver : bool
            If True then returns all keys and solvers in self. Else, only those
            which have not converged.
        """
        if overwrite_solver:
            solver_list = self.values()
            key_list = self.keys()
        else:
            solver_list = [solver for solver in self.values() if not solver.check_solved(verbose=False)]
            key_list = [key for key, solver in self.items() if not solver.check_solved(verbose=False)]
        return key_list, solver_list

    #def calculate_beam_derivatives(self):
    def calculate_direct_beam_derivative(self):
        """
        Calculate the contributions of each voxel to the direct beam
        derivative for each solver. Solvers are modified in-place.
        Used in the inverse problem.
        """
        for solver in self.values():
            solver.calculate_direct_beam_derivative()

    def calculate_microphysical_partial_derivatives(self, unknown_scatterers):
        """
        Calculates the partial derivatives of optical properties with respect
        to microphysical variables.

        solver.RTE objects in `self` are modified in place.

        Parameters
        ----------
        unknown_scatterers : at3d.containers.UnknownScatterers
            The unknown_scatterers which supply the look up tables of optical properties
            as a function of microphysical variables and the methods that should be used
            to calculate the derivatives.

        Raises
        ------
        TypeError
            If `unknown_scatterers` is not of the correct type.
        """
        if not isinstance(unknown_scatterers, UnknownScatterers):
            raise TypeError(
                "`unknown_scatterers` should be of type 'at3d.containers.UnknownScatterers'"
                " not '{}'".format(type(unknown_scatterers))
                )
        wavelength_ordered_derivatives = OrderedDict()
        for key in self:
            wavelength_ordered_derivatives[key] = OrderedDict()

        for scatterer_name, unknown_scatterer in unknown_scatterers.items():
            medium_data = list(self.values())[0].medium[scatterer_name]
            derivatives_by_wavelength = unknown_scatterer.grid_to_optical_properties.calculate_derivatives(
                list(unknown_scatterer.variables.keys()), medium_data
            )
            for key in self:
                wavelength_ordered_derivatives[key][scatterer_name] = derivatives_by_wavelength[key]

        for key, solver in self.items():
            solver.calculate_microphysical_partial_derivatives(
                wavelength_ordered_derivatives[key]
            )

class UnknownScatterers(OrderedDict):
    """
    Holds the information about which microphysical or optical
    variables to calculate derivatives for and the data and methods
    to use.

    See Also
    --------
    at3d.solver.RTE.calculate_microphysical_partial_derivatives
    at3d.medium.table_to_grid
    """
    def __init__(self, *args, global_transform=None):
        for unknown_scatterer in args:
            self.add_unknowns(unknown_scatterer)


        self.add_global_transform(global_transform)

    def add_unknowns(self, unknown_scatterer):

        if not isinstance(unknown_scatterer, at3d.medium.UnknownScatterer):
            raise TypeError(
                "`unknown_scatterer` argument should be of type '{}'".format(
                    at3d.medium.UnknownScatterer
                )
            )
        # check for consistency of unknown_scatterer with other existing unknown_scatterers
        # ie is the grid consistent.
        if self: # Test for at least one member
            reference_grid = list(self.values())[0].grid_to_optical_properties._rte_grid
            at3d.checks.check_grid_consistency(
            reference_grid,
            unknown_scatterer.grid_to_optical_properties._rte_grid
            )

        self[unknown_scatterer.scatterer_name] = unknown_scatterer

    def add_global_transform(self, transform):

        if transform is None:
            transform = at3d.transforms.CoordinateTransform()
        if not isinstance(transform, at3d.transforms.CoordinateTransform):
            raise TypeError(
                "`transform` is of an invalid type."
            )
        self.global_transform = transform
