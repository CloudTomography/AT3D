from unittest import TestCase
from collections import OrderedDict
import numpy as np
import xarray as xr
import pyshdom


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
                            source=pyshdom.source.solar(-1*np.cos(np.deg2rad(60.0)), 0.0, solarflux=1.0),
                            surface=pyshdom.surface.ocean_unpolarized(surface_wind_speed=10.0, pigmentation=0.0),
                            num_stokes=1,
                            name=None
                            )

            solvers.add_solver(wavelength,solver)

        solvers.parallel_solve(maxiter=100,verbose=False)
        cls.solvers=solvers
        cls.cloud_poly_tables = cloud_poly_tables

    def test_ssalb(self):
        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['ssalb'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        self.solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method, unknown_scatterers.table_data)
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
        solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method, unknown_scatterers.table_data)

        self.assertTrue(all([all([np.all(solvers[key]._dalb==0.0) for key in solvers]),
            all([np.all(solvers[key]._dext==1.0) for key in solvers]),
            all([np.all(solvers[key]._dleg==0.0)for key in solvers]),
            all([np.all(solvers[key]._dphasetab==0.0)for key in solvers])]))

    def test_density(self):
        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['density'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers = self.solvers
        solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method, unknown_scatterers.table_data)

        self.assertTrue(all([all([np.all(solvers[key]._dalb==0.0) for key in solvers]),
            all([np.all(solvers[key]._dext==self.cloud_poly_tables[key].extinction.interp(
                    {'reff':0.5,'veff':0.1}, method='linear').data) for key in solvers]),
                all([np.all(solvers[key]._dleg==0.0)for key in solvers]),
                all([np.all(solvers[key]._dphasetab==0.0)for key in solvers])]))

    def test_legendre(self):
        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['legendre_0_10'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers = self.solvers
        solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method, unknown_scatterers.table_data)

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
                                source=pyshdom.source.solar(-1,0.0,solarflux=1.0),
                                surface=pyshdom.surface.lambertian(albedo=0.0),
                                num_stokes=1,
                                name=None
                                )

            solvers.add_solver(wavelength,solver)
        Sensordict.get_measurements(solvers, maxiter=100, n_jobs=8, verbose=False)

        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['extinction'],cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method, unknown_scatterers.table_data)

        forward_sensors = Sensordict.make_forward_sensors()

        out,gradient,jacobian_exact = pyshdom.gradient.levis_approx_jacobian(Sensordict, solvers, forward_sensors, unknown_scatterers,
                    ([1],[1],[1]), n_jobs=8,mpi_comm=None,verbose=False,
                    maxiter=100, init_solution=False, setup_grid=False,
                    exact_single_scatter=True)

        cls.jacobian_exact = jacobian_exact

    def test_jacobian(self):
        self.assertAlmostEqual(self.jacobian_exact['jacobian_0.860'][0,0,0,0].data, 0.00396336, places=5)
