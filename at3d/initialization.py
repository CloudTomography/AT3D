"""
This module holds routines for initialization and preprocessing
based on measurements.
"""
import xarray as xr
import numpy as np

import scipy.interpolate as si
import scipy.optimize as so
import at3d

class LUTRetrieval:
    """
    Perform LUT based retrievals of homogeneous clouds
    using 1D radiative transfer. Mono-angle but arbitrary number of spectral channels.

    Currently only supports retrieval of up to two unknowns. Though Look Up Tables
    can be formed based on arbitrary number of unknowns.
    """
    def __init__(self, view_zenith, view_azimuth, sources, optical_property_generator,
                 solver_config=None,
                 surfaces=at3d.surface.lambertian(0.0),
                 background_optical_scatterers={},
                 num_stokes=1,
                 **coordinate_data):

        self.view_zenith = view_zenith
        self.view_azimuth = view_azimuth

        self.coordinate_data = coordinate_data

        if solver_config is None:
            solver_config = at3d.configuration.get_config()
            solver_config['adapt_grid_factor'] = 800
            #This should be large due to lack of base vertical resolution.
            solver_config['split_accuracy'] = 1e-3

        solver_config['ip_flag'] = 3
        self.solver_config = solver_config
        self.num_stokes = num_stokes

        self.optical_property_generator = optical_property_generator

        self.background_optical_scatterers = background_optical_scatterers

        self.optical_properties, self.rte_grid = self._prepare_optical_properties()
        if isinstance(surfaces, xr.Dataset):
            surfaces = {wavelength: surfaces for wavelength in self.optical_properties}
        if isinstance(sources, xr.Dataset):
            sources = {wavelength: sources for wavelength in self.optical_properties}

        self.surfaces = surfaces
        self.sources = sources

        self.saved_solutions = {wavelength: None for wavelength in self.optical_properties}

        self.make_lut()


    def set_surface(self, surfaces):

        if isinstance(surfaces, xr.Dataset):
            surfaces = {wavelength: surfaces for wavelength in self.optical_properties}

        self.surfaces = surfaces
        self.make_lut()

    def _prepare_optical_properties(self):

        big_meshgrid = np.meshgrid(*self.coordinate_data.values(), indexing='ij')

        # reshape into 1D grid for the RTE solution
        reshaped_data = [meshgrid.ravel() for meshgrid in big_meshgrid]

        rte_grid = at3d.grid.make_grid(0.05, reshaped_data[0].size,
                                          0.05, 1, np.linspace(0.0, 1.0, 21))

        for name, data in zip(self.coordinate_data.keys(), big_meshgrid):
            rte_grid[name] = (
                ['x', 'y', 'z'],
                np.repeat(data.ravel()[:, None, None], 21, axis=-1)
                )

        optical_properties = self.optical_property_generator(rte_grid)

        return optical_properties, rte_grid


    def rte_model(self, maxiter=800, verbose=False, n_jobs=4):
        """
        Run the SHDOM solution to form the Look Up Table for the retrieval.
        """
        sensors = at3d.containers.SensorsDict()
        for wavelength in self.optical_properties:
            sensors.add_sensor(
                wavelength,
                at3d.sensor.domaintop_projection(
                    wavelength,
                    self.rte_grid,
                    float(self.rte_grid.delx),
                    float(self.rte_grid.dely),
                    self.view_azimuth,
                    self.view_zenith,
                )
            )

        solvers_dict = at3d.containers.SolversDict()

        for wavelength, optical_property_data in self.optical_properties.items():

            medium = {'retrieved': optical_property_data}
            if wavelength in self.background_optical_scatterers:
                medium[wavelength] = self.background_optical_scatterers[wavelength]

            solver = at3d.solver.RTE(
                numerical_params=self.solver_config,
                medium=medium,
                source=self.sources[wavelength],
                surface=self.surfaces[wavelength],
                num_stokes=self.num_stokes,
                name=None
            )
            if self.saved_solutions[wavelength] is not None:
                solver.load_solution(self.saved_solutions[wavelength])

            solvers_dict.add_solver(wavelength, solver)


        sensors.get_measurements(solvers_dict, maxiter=maxiter, verbose=verbose, n_jobs=n_jobs)

        for key, solver in solvers_dict.items():
            self.saved_solutions[key] = solver.save_solution()

        return sensors, solvers_dict


    def make_lut(self):

        # we discard the direct solver output and
        sensors, solvers = self.rte_model()

        shape = [coordinate_data.size for coordinate_data in self.coordinate_data.values()]
        data_vars = {}
        for name, instrument in sensors.items():
            sensor = instrument['sensor_list'][0]
            data_vars[name] = (list(self.coordinate_data.keys()), sensor.I.data.reshape(shape))

        lut_dataset = xr.Dataset(data_vars=data_vars,
                                 coords=self.coordinate_data)

        self.lut = lut_dataset

    def retrieve(self, observed_data, initial_vector,
                minimize_options={'maxiter': 100, 'maxls': 40, 'disp': False,
                'gtol': 1e-8, 'ftol': 1e-8},
                **fixed_variables):
        """
        This only supports up to two variables.
        """
        subset_lut = self.lut
        for name, fixed_variable in fixed_variables.items():
            if self.lut[name].size > 1:
                subset_lut = subset_lut.interp({name: fixed_variable})
            else:
                subset_lut = subset_lut.sel({name: fixed_variable})

        valid_coordinates = {}
        for key in self.coordinate_data:
            if subset_lut[key].size > 1:
                valid_coordinates[key] = subset_lut[key].data

        if len(valid_coordinates) > 2:
            raise NotImplementedError(
                "This method does not currently support more than two unknowns."
            )

        retrieved_quantities_size = len(valid_coordinates)
        retrievals = np.zeros((retrieved_quantities_size, len(list(observed_data.values())[0])))
        costs = np.zeros(len(list(observed_data.values())[0]))

        for j, wavelength_data in enumerate(zip(*observed_data.values())):

            cost = 0.0
            for wavelength, data in zip(self.optical_properties, wavelength_data):
                cost += (np.log(data) - np.log(subset_lut[wavelength]))**2

            if retrieved_quantities_size == 1:
                coordinate1 = list(valid_coordinates.values())[0]
                cost_spline = si.CubicSpline(
                    coordinate1,
                    cost,
                    extrapolate=None
                    )

                def cost_function(state):

                    cost = cost_spline(state)
                    grad = cost_spline(state, nu=1)
                    return cost, np.array([grad])
                result = so.minimize(
                    cost_function,
                    x0=initial_vector,
                    jac=True,
                    method='L-BFGS-B',
                    bounds=[(coordinate1.min(), coordinate1.max())],
                    options=minimize_options
                )


            elif retrieved_quantities_size == 2:
                coordinate1 = list(valid_coordinates.values())[0]
                coordinate2 = list(valid_coordinates.values())[1]

                cost_spline = si.RectBivariateSpline(coordinate1, coordinate2,
                                                     cost,
                                                     bbox=[coordinate1.min(),
                                                           coordinate1.max(),
                                                           coordinate2.min(),
                                                           coordinate2.max(),],
                                                     kx=3, ky=3, s=0)
                def cost_function(state):

                    cost = cost_spline(state[0], state[1])
                    gradx = cost_spline(state[0], state[1], dx=1)
                    grady = cost_spline(state[0], state[1], dy=1)
                    return cost, np.array([gradx[0, 0], grady[0, 0]])

                result = so.minimize(
                    cost_function,
                    x0=initial_vector,
                    jac=True,
                    method='L-BFGS-B',
                    bounds=[(coordinate1.min(), coordinate1.max()),
                            (coordinate2.min(), coordinate2.max())],
                    options=minimize_options
                )

            for i, retrieved in enumerate(result.x):
                retrievals[i, j] = retrieved
            costs[j] = result.fun

        return retrievals, costs


def mean_ext_estimate(rte_grid, sensors, solar_mu, solar_azimuth,
                     chi=2/3, g=0.86, sun_distance_reflect=0.1,
                     sun_distance_transmit=0.1):
    """
    Estimate the extinction of a cloud using diffusion theory.

    Given a masked volume `space_carved_volume`, the geometric distance
    between each point and the sun is calculated. The value of the geometric distance
    from the sun through the cloud at the first intersection of a sensor ray
    with the cloud volume is used to classify whether sensor pixels are
    observing shadowed or directly illuminated portions of the cloud.

    The mean of all 'shadowed' and 'illuminated' pixels is used to derive an
    optical diameter using diffusion theory and the extrapolation length `chi`
    and an asymmetry factor. This optical diameter is converted to an extinction
    using the length scale of the maximum chord length through the cloud in the solar
    direction. This length scale is chosen because it collapses to the relevant case
    for several geometries.
    """
    space_carver = at3d.space_carve.SpaceCarver(rte_grid)
    if isinstance(sensors, xr.Dataset):
        sensor_list = [sensors]
    elif isinstance(sensors, type([])):
        sensor_list = sensors
    elif isinstance(sensors, at3d.containers.SensorsDict):
        sensor_list = []
        for instrument in sensors:
            sensor_list.extend(sensors[instrument]['sensor_list'])

    volume = space_carver.carve(sensor_list, agreement=(0.0, 1.0), linear_mode=False)
    sundistance = space_carver.shadow_mask(volume.mask, sensor_list, solar_mu, solar_azimuth)

    reflected = []
    transmitted = []
    for sensor in sensor_list:
        reflected.extend(sensor.I.data[np.where((sensor.sun_distance.data < sun_distance_reflect) &
                                                (sensor.cloud_mask.data == 1))])
        transmitted.extend(sensor.I.data[np.where((sensor.sun_distance.data >= sun_distance_transmit) &
                                                  (sensor.cloud_mask.data == 1))])

    sundistance_radius = sundistance.sun_distance.data[np.where(sundistance.sun_distance > 0.0)].max()

    tau_estimate = 2*chi*np.mean(reflected)/np.mean(transmitted)/(1.0-g)
    ext_estimate = tau_estimate/sundistance_radius

    extinction = np.zeros(volume.mask.shape)
    extinction[np.where(volume.mask == 1.0)] = ext_estimate
    extinction = xr.Dataset(
        data_vars={
            'extinction': (['x', 'y', 'z'], extinction)
        },
        coords={
            'x': rte_grid.x,
            'y': rte_grid.y,
            'z': rte_grid.z,
        },
        attrs={
            'tau_estimate': tau_estimate,
            'chi': chi,
            'g': g,
            'radius': sundistance_radius
        }
    )
    return extinction
