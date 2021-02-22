from unittest import TestCase
from collections import OrderedDict
import numpy as np
import xarray as xr
import pyshdom

class CostFunctionL2(TestCase):
    @classmethod
    def setUpClass(cls):
        cost = 0.0
        gradout = np.zeros((10, 1, 1))
        raygrad_pixel = np.ones((4, 10, 1))
        uncertainties = np.ones((4, 4))*5
        costfunc = 'L2'
        stokesout = np.ones(4)*10.0
        stokesout[3] = 0.0
        measurement = np.ones(4)*13.0
        gradout, cost = pyshdom.core.update_costfunction(
            cost=cost,
            gradout=gradout,
            stokesout=stokesout,
            measurement=measurement,
            raygrad_pixel=raygrad_pixel,
            uncertainties=uncertainties,
            costfunc=costfunc,
        )
        cls.gradout = gradout
        cls.cost = cost
    def test_cost(self):
        self.assertAlmostEqual(self.cost, 1960.0, places=5)
    def test_gradient(self):
        self.assertAlmostEqual(self.gradout[0, 0, 0], -440.0, places=5)

class CostFunctionLL(TestCase):
    @classmethod
    def setUpClass(cls):
        cost = 0.0
        gradout = np.zeros((10, 1, 1))
        raygrad_pixel = np.ones((3, 10, 1))
        uncertainties = np.zeros((2, 2))
        uncertainties[0, 0] = (1.0/0.03)**2
        uncertainties[1, 1] = (1.0/0.005)**2
        costfunc = 'LL'
        stokesout = np.ones(3)
        stokesout[1] = 0.5
        stokesout[2] = 0.0
        measurement = np.ones(3)*1.25
        measurement[1] = 0.25
        measurement[2] = 0.25

        gradout, cost = pyshdom.core.update_costfunction(
            cost=cost,
            gradout=gradout,
            stokesout=stokesout,
            measurement=measurement,
            raygrad_pixel=raygrad_pixel,
            uncertainties=uncertainties,
            costfunc=costfunc,
        )
        cls.gradout = gradout
        cls.cost = cost
    def test_cost(self):
        self.assertAlmostEqual(self.cost, 6519.21, places=2)
    def test_gradient(self):
        self.assertAlmostEqual(self.gradout[0, 0, 0], 45329.43, places=2)


class Microphysical_Derivatives(TestCase):
    @classmethod
    def setUpClass(cls):

        reff = 0.5
        ext = 20.0
        veff = 0.1
        rte_grid = pyshdom.grid.make_grid(0.05, 13, 0.05, 13, np.linspace(0.1,0.7,13))
        grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
        rte_grid['density'] = (['x','y','z'], np.ones(grid_shape))
        rte_grid['reff'] = (['x','y','z'], np.zeros(grid_shape) + reff)
        rte_grid['veff'] = (['x','y','z'], np.zeros(grid_shape) + veff)

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, rte_grid)

        #define any necessary variables for microphysics here.
        size_distribution_function = pyshdom.size_distribution.gamma

        #define sensors.
        Sensordict = pyshdom.containers.SensorsDict()
        sensor_zenith_list = [75.0,60.0,45.6,26.1]*2 + [0.0]
        sensor_azimuth_list = [90]*4 + [-90]*4 +[0.0]
        wavelengths = [0.86,0.86,0.86,1.38,1.38,2.2,2.2,3.4,3.4]
        for zenith,azimuth,wavelength in zip(sensor_zenith_list,sensor_azimuth_list,wavelengths):
            Sensordict.add_sensor('MISR', pyshdom.sensor.orthographic_projection(
                                                                wavelength,
                                                                cloud_scatterer_on_rte_grid,
                                                                0.04, 0.04,
                                                                azimuth, zenith,
                                                                altitude='TOA', stokes=['I']
                                                                )
                                                                )

        #Define the RTE solvers needed to model the measurements and
        #calculate optical properties.
        wavelengths = Sensordict.get_unique_solvers()

        cloud_poly_tables = OrderedDict()
        solvers = pyshdom.containers.SolversDict()


        for wavelength in wavelengths:

            #optical properties from mie calculations.
            mie_mono_table = pyshdom.mie.get_mono_table('Water',(wavelength,wavelength),
                                                      max_integration_radius=10.0,
                                                      minimum_effective_radius=0.1,
                                                      relative_dir='../mie_tables',
                                                      verbose=False)
            cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
                                                                    mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function,particle_density=1.0,
                                reff={'coord_min':0.2, 'coord_max': 1.0, 'npoints': 10,
                                'spacing': 'logarithmic', 'units': 'micron'},
                                veff={'coord_min':0.09, 'coord_max': 0.11, 'npoints': 12,
                                'spacing': 'linear', 'units': 'unitless'}
                                )
            poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            optical_properties = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                             exact_table=False)
            optical_properties['ssalb'][:,:,:] = 1.0
            extinction = np.zeros(optical_properties.extinction.shape)
            np.random.seed(1)
            extinction[1:-1,1:-1,1:-1] = ext + np.random.uniform(low=0.0,high=10.0,size=(11,11,11))
            #extinction[a,b,c] += step
            optical_properties['legcoef'][:,1:,:] = 0.0
            optical_properties['extinction'][:,:,:] = extinction
            cloud_poly_tables[wavelength] = poly_table
            config = pyshdom.configuration.get_config('../default_config.json')
            config['num_mu_bins'] = 4
            config['num_phi_bins'] = 8
            config['split_accuracy'] = 0.1
            config['spherical_harmonics_accuracy'] = 0.0
            config['solution_accuracy'] = 1e-4
            config['deltam'] = True
            solver = pyshdom.solver.RTE(
                            numerical_params=config,
                            medium={'cloud': optical_properties},
                            source=pyshdom.source.solar(wavelength,-1*np.cos(np.deg2rad(60.0)), 0.0, solarflux=1.0),
                            surface=pyshdom.surface.ocean_unpolarized(surface_wind_speed=10.0, pigmentation=0.0),
                            num_stokes=1,
                            name=None
                            )

            solvers.add_solver(wavelength, solver)

        solvers.parallel_solve(maxiter=100, verbose=False)
        cls.solvers=solvers
        cls.cloud_poly_tables = cloud_poly_tables

    def test_ssalb(self):
        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['ssalb'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        self.solvers.add_microphysical_partial_derivatives(unknown_scatterers)
        solvers = self.solvers
        self.assertTrue(all([all([np.all(solvers[key]._dalb==1.0) for key in solvers]),
            all([np.all(solvers[key]._dext==0.0) for key in solvers]),
            all([np.all(solvers[key]._dleg==0.0)for key in solvers]),
            all([np.all(solvers[key]._dphasetab==0.0)for key in solvers])]))

    def test_extinction(self):
        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['extinction'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers = self.solvers
        solvers.add_microphysical_partial_derivatives(unknown_scatterers)

        self.assertTrue(all([all([np.all(solvers[key]._dalb==0.0) for key in solvers]),
            all([np.all(solvers[key]._dext==1.0) for key in solvers]),
            all([np.all(solvers[key]._dleg==0.0)for key in solvers]),
            all([np.all(solvers[key]._dphasetab==0.0)for key in solvers])]))

    def test_density(self):
        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['density'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers = self.solvers
        solvers.add_microphysical_partial_derivatives(unknown_scatterers)

        self.assertTrue(all([all([np.all(solvers[key]._dalb==0.0) for key in solvers]),
            all([np.allclose(solvers[key]._dext, self.cloud_poly_tables[key].extinction.interp(
                    {'reff':0.5,'veff':0.1}, method='linear').data) for key in solvers]),
                all([np.all(solvers[key]._dleg==0.0)for key in solvers]),
                all([np.all(solvers[key]._dphasetab==0.0)for key in solvers])]))

    def test_legendre(self):
        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['legendre_0_10'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers = self.solvers
        solvers.add_microphysical_partial_derivatives(unknown_scatterers)

        self.assertTrue(all([all([np.all(solvers[key]._dalb==0.0) for key in solvers]),
            all([np.all(solvers[key]._dext==0.0) for key in solvers]),
                all([np.all(solvers[key]._dleg[0,10]==1.0/(2*10.0 + 1.0))for key in solvers]),
                ]))


class Verify_Jacobian(TestCase):
    @classmethod
    def setUpClass(cls):

        ext = 0.1
        veff = 0.1
        reff=10.0
        rte_grid = pyshdom.grid.make_grid(0.05, 3, 0.05, 3, np.arange(0.1, 0.25, 0.05))
        grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
        rte_grid['density'] = (['x','y','z'], np.ones(grid_shape))
        rte_grid['reff'] = (['x','y','z'], np.zeros(grid_shape) + reff)
        rte_grid['veff'] = (['x','y','z'] ,np.zeros(grid_shape) + veff)

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, rte_grid)

        #define any necessary variables for microphysics here.
        size_distribution_function = pyshdom.size_distribution.gamma

        #define sensors.
        Sensordict = pyshdom.containers.SensorsDict()

        sensor = pyshdom.sensor.make_sensor_dataset(np.array([0.05]),np.array([0.05]),
                                                  np.array([0.7]),np.array([1.0]),
                                                  np.array([0.0]),
                                                 wavelength=0.86,stokes=['I'],fill_ray_variables=True)
        Sensordict.add_sensor('MISR', sensor)

        #Define the RTE solvers needed to model the measurements and
        #calculate optical properties.
        wavelengths = Sensordict.get_unique_solvers()

        cloud_poly_tables = OrderedDict()
        solvers = pyshdom.containers.SolversDict()

        for wavelength in wavelengths:

            #optical properties from mie calculations.
            mie_mono_table = pyshdom.mie.get_mono_table('Water',(wavelength,wavelength),
                                                      max_integration_radius=65.0,
                                                      minimum_effective_radius=0.1,
                                                      relative_dir='../mie_tables',
                                                      verbose=False)
            cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
                                                                    mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function,particle_density=1.0,
                                reff={'coord_min':9.0, 'coord_max': 11.0, 'npoints': 100,
                                'spacing': 'logarithmic', 'units': 'micron'},
                                veff={'coord_min':0.09, 'coord_max': 0.11, 'npoints': 12,
                                'spacing': 'linear', 'units': 'unitless'}
                                )
            poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            optical_properties = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                             exact_table=False)
            optical_properties['ssalb'][:,:,:] = 1.0
            extinction = np.zeros(optical_properties.extinction.shape)
            extinction[1:-1,1:-1,1:-1] = ext
            #extinction[a,b,c] += step
            optical_properties['legcoef'][:,1:,:] = 0.0
            optical_properties['extinction'][:,:,:] = extinction
            cloud_poly_tables[wavelength] = poly_table
            config = pyshdom.configuration.get_config('../default_config.json')
            config['num_mu_bins'] = 16
            config['num_phi_bins'] = 32
            config['split_accuracy'] = 0.0003
            config['spherical_harmonics_accuracy'] = 0.0
            config['solution_accuracy'] = 1e-5
            config['deltam'] = True
            solver = pyshdom.solver.RTE(
                                numerical_params=config,
                                medium={'cloud': optical_properties},
                                source=pyshdom.source.solar(wavelength,-1,0.0,solarflux=1.0),
                                surface=pyshdom.surface.lambertian(albedo=0.0),
                                num_stokes=1,
                                name=None
                                )

            solvers.add_solver(wavelength, solver)
        Sensordict.get_measurements(solvers, maxiter=100, n_jobs=8, verbose=False)

        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['extinction'], cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers.add_microphysical_partial_derivatives(unknown_scatterers)

        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = pyshdom.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'n_jobs':4, 'maxiter': 100, 'setup_grid':True, 'verbose':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': ([1],[1],[1])}, uncertainty_kwargs={'add_noise': False})
        out, gradient, jacobian_exact = gradient_call()

        cls.jacobian_exact = jacobian_exact

    def test_jacobian(self):
        self.assertAlmostEqual(self.jacobian_exact['jacobian_0.860'][0,0,0,0].data, 0.00396336, places=5)

class ThermalJacobian(TestCase):
    @classmethod
    def setUpClass(cls):

        ext = 30
        veff = 0.1
        reff = 10.0
        rte_grid = pyshdom.grid.make_grid(0.05, 15, 0.05, 15, np.arange(0.1, 0.75, 0.05))
        grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
        rte_grid['density'] = (['x','y','z'], np.ones(grid_shape))
        rte_grid['reff'] = (['x','y','z'], np.zeros(grid_shape) + reff)
        rte_grid['veff'] = (['x','y','z'] ,np.zeros(grid_shape) + veff)

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, rte_grid)

        #define any necessary variables for microphysics here.
        size_distribution_function = pyshdom.size_distribution.gamma

        #define sensors.
        Sensordict = pyshdom.containers.SensorsDict()

        sensor = pyshdom.sensor.make_sensor_dataset(np.array([0.25]),np.array([0.25]),
                                                  np.array([0.7]),np.array([1.0]),
                                                  np.array([0.0]),
                                                 wavelength=11.0,stokes=['I'],fill_ray_variables=True)
        Sensordict.add_sensor('MISR', sensor)

        #Define the RTE solvers needed to model the measurements and
        #calculate optical properties.
        wavelengths = Sensordict.get_unique_solvers()

        cloud_poly_tables = OrderedDict()
        solvers = pyshdom.containers.SolversDict()

        atmosphere = xr.Dataset(data_vars={
            'temperature': (['x','y','z'], np.ones(grid_shape)*280.0)
        }, coords={'x':rte_grid.x, 'y': rte_grid.y, 'z': rte_grid.z})
        atmosphere2 = pyshdom.grid.resample_onto_grid(rte_grid, atmosphere)

        mie_mono_table = pyshdom.mie.get_mono_table('Water',(11.0,11.0),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='../mie_tables',
                                                  verbose=False)
        for wavelength in wavelengths:

            cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
                                                                    mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function,particle_density=1.0,
                                reff={'coord_min':9.0, 'coord_max': 11.0, 'npoints': 100,
                                'spacing': 'logarithmic', 'units': 'micron'},
                                veff={'coord_min':0.09, 'coord_max': 0.11, 'npoints': 12,
                                'spacing': 'linear', 'units': 'unitless'}
                                )
            poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            optical_properties = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                             exact_table=False)
            optical_properties['ssalb'][:,:,:] = 0.0#9999
            extinction = np.zeros(optical_properties.extinction.shape)
            extinction[3:-3,3:-3,3:-3] = ext
            extinction[3:-3,3:-3,-3:] = 0.0
            #extinction[a,b,c] += step
            #optical_properties['legcoef'][:,1:,:] = 0.0
            optical_properties['extinction'][:,:,:] = extinction
            cloud_poly_tables[wavelength] = poly_table
            config = pyshdom.configuration.get_config('../default_config.json')
            config['num_mu_bins'] = 16
            config['num_phi_bins'] = 32
            config['split_accuracy'] = 0.1
            config['spherical_harmonics_accuracy'] = 0.0
            config['adapt_grid_factor'] = 100.0
            config['solution_accuracy'] = 1e-5
            config['deltam'] = True
            solver = pyshdom.solver.RTE(
                                numerical_params=config,
                                medium={'cloud': optical_properties},
                                source=pyshdom.source.thermal(wavelength),
                                surface=pyshdom.surface.lambertian(albedo=0.0),
                                num_stokes=1,
                                atmosphere=atmosphere2,
                                name=None
                                )

            solvers.add_solver(wavelength,solver)
        Sensordict.get_measurements(solvers, maxiter=100, n_jobs=8, verbose=False)

        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['extinction'], cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers.add_microphysical_partial_derivatives(unknown_scatterers)

        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = pyshdom.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':100,'n_jobs':1, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': np.where(extinction > 0.0)}, uncertainty_kwargs={'add_noise': False})

        cost, gradient, jacobian = gradient_call()
        derivs = np.zeros(extinction.shape)
        derivs[np.where(extinction > 0.0)] = jacobian['jacobian_11.000'][0,0,:,0].data
        cls.jacobian = derivs[5, 5]

    def test_jacobian(self):
        self.assertAlmostEqual(self.jacobian.sum(), 0.8255725, places=7)
