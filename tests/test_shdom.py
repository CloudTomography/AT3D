from unittest import TestCase
from collections import OrderedDict
import numpy as np
import xarray as xr
import shdom

def parse_shdom_output(filename, comment='*'):
    """
    Parse output from original SHDOM files
    """
    output = []
    with open(filename) as file:
        data = file.readlines()
        for line in data:
            if not comment in line:
                output.append(np.fromstring(line, sep=' '))
    return np.array(output)

def get_basic_state_for_surface():
    config = shdom.configuration.get_config('../default_config.json')
    config['split_accuracy'] = 0.001
    config['spherical_harmonics_accuracy'] = 0.0
    config['num_mu_bins'] = 16
    config['num_phi_bins'] = 32
    config['solution_accuracy'] = 1e-5
    config['x_boundary_condition'] = 'periodic'
    config['y_boundary_condition'] = 'periodic'
    config['ip_flag'] = 3

    rte_grid = shdom.grid.make_grid(0.02, 50, 0.02, 1,
                               np.array([0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0]))

    atmosphere = xr.Dataset(
        data_vars = {
            'temperature': ('z', np.array([288.0,269.0,249.0,230.0,217.0,217.0,217.0,218.0,221.0,224.0,227.0])),
            'pressure': ('z', np.ones(rte_grid.z.size)*1013.25)
                        },
        coords = {'z': rte_grid.z.data}
    )
    wavelengths = np.atleast_1d(0.85)
    rayleigh= shdom.rayleigh.to_grid(wavelengths,atmosphere,rte_grid)
    rayleigh[0.85]['extinction'] = (['x','y','z'],np.round(rayleigh[0.85].extinction.data,4))

    x = np.linspace(0,1.0-1.0/50,50)
    y = np.zeros(50)
    z = np.ones(50)*30.0
    mu = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]+[1.0]+[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9][::-1])
    phi = np.array([180.0]*9+[0.0]*10)

    x2 = np.tile(x,19)
    y2 = np.tile(y,19)
    z2 = np.tile(z,19)
    mu2 = np.repeat(mu,50)
    phi2 = np.repeat(phi,50)

    sensor = shdom.sensor.make_sensor_dataset(x2.ravel(),y2.ravel(),z2.ravel(),mu2.ravel(),np.deg2rad(phi2.ravel()),['I'],
                                             0.85, fill_ray_variables=True)
    return sensor, rayleigh, config

class Verify_Solver(TestCase):
    @classmethod
    def setUpClass(cls):
        #shdom_polarized = xr.open_dataset('../shdom_verification/shdomout_rico32x36x26w672_polarized.nc')

        cloud_scatterer = shdom.grid.load_2parameter_lwc_file('../shdom_verification/rico32x36x26.txt',
                                                   density='lwc',origin=(0.0,0.0))

        rte_grid = shdom.grid.make_grid(cloud_scatterer.x.data[1] -cloud_scatterer.x.data[0],cloud_scatterer.x.data.size,
                                   cloud_scatterer.y.data[1] -cloud_scatterer.y.data[0],cloud_scatterer.y.data.size,
                                   np.append(np.array([0.0]),cloud_scatterer.z.data))

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = shdom.grid.resample_onto_grid(rte_grid, cloud_scatterer)

        #define any necessary variables for microphysics here.
        size_distribution_function = shdom.size_distribution.gamma
        #We choose a gamma size distribution and therefore need to define a 'veff' variable.
        cloud_scatterer_on_rte_grid['alpha'] = (cloud_scatterer_on_rte_grid.reff.dims,
                                               np.full_like(cloud_scatterer_on_rte_grid.reff.data, fill_value=7))

        wavelengths = np.atleast_1d(0.672)
        x,y = np.meshgrid(np.arange(0.0,0.62,0.0133), np.arange(0.0,0.70,0.0133))

        x=x.ravel()
        y=y.ravel()

        mu = np.array([1.0,0.5,0.2,0.5,0.2])
        phi = np.array([0.0,0.0,0.0,90.0,90.0])
        z = np.ones(x.size)*1.0
        x2 = np.tile(x,5)
        y2 = np.tile(y,5)
        z2 = np.tile(z,5)
        mu2 = np.repeat(mu,x.size)
        phi2 = np.repeat(phi,x.size)

        sensor = shdom.sensor.make_sensor_dataset(x2.ravel(),y2.ravel(),z2.ravel(),mu2.ravel(),np.deg2rad(phi2.ravel()),['I'],
                                                 0.672, fill_ray_variables=True)
        Sensordict = shdom.util.SensorsDict()
        Sensordict.add_sensor('Sensor0', sensor)

        config = shdom.configuration.get_config('../default_config.json')
        config['split_accuracy'] = 0.1
        config['spherical_harmonics_accuracy'] = 0.01
        config['num_mu_bins'] = 8
        config['num_phi_bins'] = 16
        config['solution_accuracy'] = 1e-4

        solvers = shdom.util.SolversDict()
        for wavelength in wavelengths:

            #mie table and medium doesn't matter here as it is overwritten by property file.
            #just need it to initialize a solver.
            mie_mono_table = shdom.mie.get_mono_table('Water',(wavelength,wavelength), relative_path='../mie_tables')
            cloud_size_distribution = shdom.size_distribution.get_size_distribution_grid(
                                mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function,
                                particle_density=1.0,
                                reff=[5.0, 20.0, 50, 'linear', 'micron'],
                                alpha=[6.9, 7.1, 2, 'linear', 'unitless']
                                )
            poly_table = shdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            cloud_optical_scatterer = shdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)

            solvers.add_solver(wavelength, shdom.solver.RTE(
                                            numerical_params=config,
                                            medium={'cloud': cloud_optical_scatterer},
                                            source=shdom.source.solar(-0.5, 0.0, solarflux=1.0),
                                            surface=shdom.surface.lambertian(albedo=0.05),
                                            num_stokes=3,
                                            name=None,
                                            atmosphere=None
                                            )
                                        )

        solver = solvers[0.672]
        solver._solve_prop(filename = '../shdom_verification/rico32x36x26w672.prp')
        shdom_source = parse_shdom_output('../shdom_verification/shdom_verification_source_out.out', comment='*')
        cls.truth = shdom_source.T
        cls.testing = solver._source[:, :solver._npts]
        integrated_rays = solver.integrate_to_sensor(sensor)
        cls.integrated_rays = integrated_rays
        cls.radiances = parse_shdom_output('../shdom_verification/rico32x36x26w672ar.out', comment='!')

    def test_solver(self):
        self.assertTrue(np.allclose(self.testing, self.truth))

    def test_radiance(self):
        self.assertTrue(np.allclose(self.integrated_rays.I.data, self.radiances[:,2].data, atol=3e-3))

    def test_Q(self):
        self.assertTrue(np.allclose(self.integrated_rays.Q.data, self.radiances[:,3].data, atol=2e-4))

    def test_U(self):
        self.assertTrue(np.allclose(self.integrated_rays.U.data, self.radiances[:,4].data, atol=7e-5))


class Parallelization_No_SubpixelRays(TestCase):
    @classmethod
    def setUpClass(cls):

        reff = 0.5
        ext = 20.0
        rte_grid = shdom.grid.make_grid(0.05,13,
                                   0.05,13,
                                   np.linspace(0.1,0.7,13))

        rte_grid['density'] = (['x','y','z'], np.ones((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)))
        rte_grid['reff'] = (['x','y','z'], np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+reff)
        rte_grid['veff'] = (['x','y','z'] ,np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+0.1)

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = shdom.grid.resample_onto_grid(rte_grid, rte_grid)

        #define any necessary variables for microphysics here.
        size_distribution_function = shdom.size_distribution.gamma

        #define sensors.
        Sensordict = shdom.util.SensorsDict()
        Sensordict2 = shdom.util.SensorsDict()
        misr_list = []
        sensor_zenith_list = [75.0,60.0,45.6,26.1]*2 + [0.0]
        sensor_azimuth_list = [90]*4 + [-90]*4 +[0.0]
        wavelengths = [0.86,0.86,0.86,1.38,1.38,2.2,2.2,3.4,3.4]
        for zenith,azimuth,wavelength in zip(sensor_zenith_list,sensor_azimuth_list,wavelengths):
            Sensordict.add_sensor('MISR',
                            shdom.sensor.orthographic_projection(wavelength, cloud_scatterer_on_rte_grid,0.04,0.04, azimuth, zenith,
                                                     altitude='TOA', stokes=['I'])
                                 )
            Sensordict2.add_sensor('MISR',
                            shdom.sensor.orthographic_projection(wavelength, cloud_scatterer_on_rte_grid,0.04,0.04, azimuth, zenith,
                                                     altitude='TOA', stokes=['I'])
            )

        #Define the RTE solvers needed to model the measurements and
        #calculate optical properties.
        wavelengths = Sensordict.get_unique_solvers()

        cloud_poly_tables = OrderedDict()
        solvers = shdom.util.SolversDict()

        for wavelength in wavelengths:

            #optical properties from mie calculations.
            mie_mono_table = shdom.mie.get_mono_table('Water',(wavelength,wavelength),
                                                      max_integration_radius=10.0,
                                                      minimum_effective_radius=0.1,
                                                      relative_path='../mie_tables',
                                                      verbose=False)
            #mie_mono_table.to_netcdf('./mie_tables/MieTable_{}.nc'.format(reff))
            cloud_size_distribution = shdom.size_distribution.get_size_distribution_grid(
                                                                    mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function,particle_density=1.0,
                                reff=[0.2,1.0,10,'logarithmic','micron'],
                                veff=[0.09,0.11,12,'linear','unitless'],
                                )
            poly_table = shdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            optical_properties = shdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                             exact_table=False)
            optical_properties['ssalb'][:,:,:] = 1.0
            extinction = np.zeros(optical_properties.extinction.shape)
            extinction[1:-1,1:-1,1:-1] = ext + np.random.uniform(low=0.0,high=10.0,size=(11,11,11))
            #extinction[a,b,c] += step
            optical_properties['legcoef'][:,1:,:] = 0.0
            optical_properties['extinction'][:,:,:] = extinction
            cloud_poly_tables[wavelength] = poly_table
            config = shdom.configuration.get_config('../default_config.json')
            config['num_mu_bins'] = 4
            config['num_phi_bins'] = 8
            config['split_accuracy'] = 0.1
            config['spherical_harmonics_accuracy'] = 0.0
            config['solution_accuracy'] = 1e-4
            solver = shdom.solver.RTE(numerical_params=config,
                                            medium={'cloud': optical_properties},
                                           source=shdom.source.solar(-1*np.cos(np.deg2rad(60.0)),0.0,solarflux=1.0),
                                           surface=shdom.surface.ocean_unpolarized(surface_wind_speed=10.0,
                                                                                  pigmentation=0.0),
                                            num_stokes=1,
                                            name=None
                                           )
            solvers.add_solver(wavelength,solver)
            shdom.util.get_measurements(solvers, Sensordict, maxiter=100, n_jobs=8, verbose=False)
            shdom.util.get_measurements(solvers, Sensordict2, maxiter=100, n_jobs=1, verbose=False)
            #Sensordict['MISR']['sensor_list'][0].to_netcdf('../shdom_verification/RenderedSensorReference_nosubpixel.nc')

            cls.solvers = solvers
            cls.Sensordict = Sensordict
            cls.Sensordict2 = Sensordict2

    def test_radiance(self):
        self.assertTrue(all([np.allclose(self.Sensordict2['MISR']['sensor_list'][i].I, self.Sensordict['MISR']['sensor_list'][i].I) for i in range(9)]))

    def test_no_subpixel_args_reference(self):
        test = xr.open_dataset('../shdom_verification/RenderedSensorReference_nosubpixel.nc')
        self.assertTrue(test.equals(self.Sensordict['MISR']['sensor_list'][0]))


class Parallelization_SubpixelRays(TestCase):
    @classmethod
    def setUpClass(cls):

        reff = 0.5
        ext = 20.0
        rte_grid = shdom.grid.make_grid(0.05,13,
                                   0.05,13,
                                   np.linspace(0.1,0.7,13))

        rte_grid['density'] = (['x','y','z'], np.ones((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)))
        rte_grid['reff'] = (['x','y','z'], np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+reff)
        rte_grid['veff'] = (['x','y','z'] ,np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+0.1)

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = shdom.grid.resample_onto_grid(rte_grid, rte_grid)

        #define any necessary variables for microphysics here.
        size_distribution_function = shdom.size_distribution.gamma

        #define sensors.
        Sensordict = shdom.util.SensorsDict()
        Sensordict2 = shdom.util.SensorsDict()
        misr_list = []
        sensor_zenith_list = [75.0,60.0,45.6,26.1]*2 + [0.0]
        sensor_azimuth_list = [90]*4 + [-90]*4 +[0.0]
        wavelengths = [0.86,0.86,0.86,1.38,1.38,2.2,2.2,3.4,3.4]
        for zenith,azimuth,wavelength in zip(sensor_zenith_list,sensor_azimuth_list,wavelengths):
            Sensordict.add_sensor('MISR',
                            shdom.sensor.orthographic_projection(wavelength, cloud_scatterer_on_rte_grid,0.04,0.04, azimuth, zenith,
                                                     altitude='TOA', stokes=['I'],
                                                     sub_pixel_ray_args={'method': shdom.sensor.gaussian,
                                                     'degree':4})
                                 )
            Sensordict2.add_sensor('MISR',
                            shdom.sensor.orthographic_projection(wavelength, cloud_scatterer_on_rte_grid,0.04,0.04, azimuth, zenith,
                                                     altitude='TOA', stokes=['I'],
                                                  sub_pixel_ray_args={'method': shdom.sensor.gaussian,
                                                  'degree':4})
            )

        #Define the RTE solvers needed to model the measurements and
        #calculate optical properties.
        wavelengths = Sensordict.get_unique_solvers()

        cloud_poly_tables = OrderedDict()
        solvers = shdom.util.SolversDict()

        for wavelength in wavelengths:

            #optical properties from mie calculations.
            mie_mono_table = shdom.mie.get_mono_table('Water',(wavelength,wavelength),
                                                      max_integration_radius=10.0,
                                                      minimum_effective_radius=0.1,
                                                      relative_path='../mie_tables',
                                                      verbose=False)

            cloud_size_distribution = shdom.size_distribution.get_size_distribution_grid(
                                                                    mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function,particle_density=1.0,
                                reff=[0.2,1.0,10,'logarithmic','micron'],
                                veff=[0.09,0.11,12,'linear','unitless'],
                                )
            poly_table = shdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            optical_properties = shdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                             exact_table=False)
            optical_properties['ssalb'][:,:,:] = 1.0
            extinction = np.zeros(optical_properties.extinction.shape)
            extinction[1:-1,1:-1,1:-1] = ext + np.random.uniform(low=0.0,high=10.0,size=(11,11,11))
            #extinction[a,b,c] += step
            optical_properties['legcoef'][:,1:,:] = 0.0
            optical_properties['extinction'][:,:,:] = extinction
            cloud_poly_tables[wavelength] = poly_table
            config = shdom.configuration.get_config('../default_config.json')
            config['num_mu_bins'] = 4
            config['num_phi_bins'] = 8
            config['split_accuracy'] = 0.1
            config['spherical_harmonics_accuracy'] = 0.0
            config['solution_accuracy'] = 1e-4
            solver = shdom.solver.RTE(numerical_params=config,
                                            medium={'cloud': optical_properties},
                                           source=shdom.source.solar(-1*np.cos(np.deg2rad(60.0)),0.0,solarflux=1.0),
                                           surface=shdom.surface.ocean_unpolarized(surface_wind_speed=10.0,
                                                                                  pigmentation=0.0),
                                            num_stokes=1,
                                            name=None
                                           )
            solvers.add_solver(wavelength,solver)
            shdom.util.get_measurements(solvers, Sensordict, maxiter=100, n_jobs=8, verbose=False)
            shdom.util.get_measurements(solvers, Sensordict2, maxiter=100, n_jobs=1, verbose=False)

            #Sensordict['MISR']['sensor_list'][0].to_netcdf('../shdom_verification/RenderedSensorReference_subpixelargs.nc')
            cls.solvers = solvers
            cls.Sensordict = Sensordict
            cls.Sensordict2 = Sensordict2

    def test_subpixel_args(self):
        self.assertTrue(all([np.allclose(self.Sensordict2['MISR']['sensor_list'][i].I, self.Sensordict['MISR']['sensor_list'][i].I) for i in range(9)]))

    def test_subpixel_args_reference(self):
        test = xr.open_dataset('../shdom_verification/RenderedSensorReference_subpixelargs.nc')
        self.assertTrue(test.equals(self.Sensordict['MISR']['sensor_list'][0]))


class Verify_Lambertian_Surfaces(TestCase):
    @classmethod
    def setUpClass(cls):

        sensor, rayleigh, config = get_basic_state_for_surface()
        variable_lambertian = shdom.surface.lambertian(albedo=np.linspace(0.0,0.3- (0.3-0.0)/50.0,50)[:,np.newaxis],
                                                       ground_temperature=288.0,delx=0.02,dely=0.04)
        # variable_ocean = shdom.surface.ocean_unpolarized(surface_wind_speed = np.linspace(4.0, 12.0-(12.0-4.0)/50.0,50)[:,np.newaxis],
        #                                                 pigmentation = np.zeros((50,1)), ground_temperature=288.0,delx=0.02,dely=0.04)
        solver = shdom.solver.RTE(numerical_params=config,
                                        medium={'rayleigh': rayleigh[0.85]},
                                       source=shdom.source.solar(-0.707, 0.0, solarflux=1.0),
                                       surface=variable_lambertian,
                                        num_stokes=1,
                                        name=None,
                                        atmosphere=None)

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver
        cls.fluxes = parse_shdom_output('../shdom_verification/brdf_L1f.out', comment='!')
        cls.radiances = parse_shdom_output('../shdom_verification/brdf_L1r.out', comment='!')

    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_up(self):
        self.assertTrue(np.allclose(self.fluxes[:,2],self.solver.fluxes.flux_up[:,0,0].data, atol=4e-6))

    def test_radiance(self):
        self.assertTrue(np.allclose(self.radiances[:,-1], self.integrated_rays.I.data, atol=9e-6))


class Verify_Ocean_Unpolarized_Surfaces(TestCase):
    @classmethod
    def setUpClass(cls):

        sensor, rayleigh, config = get_basic_state_for_surface()
        variable_ocean = shdom.surface.ocean_unpolarized(
                                surface_wind_speed = np.linspace(4.0, 12.0-(12.0-4.0)/50.0,50)[:,np.newaxis],
                                pigmentation = np.zeros((50,1)),
                                ground_temperature=288.0,
                                delx=0.02, dely=0.04
                                )
        solver = shdom.solver.RTE(
                                numerical_params=config,
                                medium={'rayleigh': rayleigh[0.85]},
                                source=shdom.source.solar(-0.707, 0.0, solarflux=1.0),
                                surface=variable_ocean,
                                num_stokes=1,
                                name=None,
                                atmosphere=None
                                )

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver
        cls.fluxes = parse_shdom_output('../shdom_verification/brdf_O1f.out', comment='!')
        cls.radiances = parse_shdom_output('../shdom_verification/brdf_O1r.out', comment='!')

    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_up(self):
        self.assertTrue(np.allclose(self.fluxes[:,2],self.solver.fluxes.flux_up[:,0,0].data, atol=4e-6))

    def test_radiance(self):
        self.assertTrue(np.allclose(self.radiances[:,-1], self.integrated_rays.I.data, atol=9e-6))


class Verify_RPV_Surfaces(TestCase):
    @classmethod
    def setUpClass(cls):

        sensor, rayleigh, config = get_basic_state_for_surface()
        k = np.linspace(0.5, 1.0 - 0.5/50, 50)[:,np.newaxis]
        theta = np.ones(50)[:,np.newaxis]*-0.24
        rho0 = np.ones(50)[:,np.newaxis]*0.1
        surface = shdom.surface.RPV_unpolarized(rho0,k,theta, ground_temperature=288.0, delx=0.02,dely=0.02)
        solver = shdom.solver.RTE(numerical_params=config,
                                                medium={'rayleigh': rayleigh[0.85]},
                                               source=shdom.source.solar(-0.707, 0.0, solarflux=1.0),
                                               surface=surface,
                                                num_stokes=1,
                                                name=None,
                                                atmosphere=None)

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver
        cls.fluxes =  parse_shdom_output('../shdom_verification/brdf_R1f.out', comment='!')
        cls.radiances = parse_shdom_output('../shdom_verification/brdf_R1r.out', comment='!')

    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_up(self):
        self.assertTrue(np.allclose(self.fluxes[:,2],self.solver.fluxes.flux_up[:,0,0].data, atol=4e-6))

    def test_radiance(self):
        self.assertTrue(np.allclose(self.radiances[:,-1], self.integrated_rays.I.data, atol=9e-6))


class Verify_WaveFresnel_Surfaces(TestCase):
    @classmethod
    def setUpClass(cls):

        sensor, rayleigh, config = get_basic_state_for_surface()
        real_refractive_index =np.ones(50)[:,np.newaxis]*1.33
        imaginary_refractive_index = np.zeros(50)[:,np.newaxis]
        surface_wind_speed = np.linspace(4.0, 12.0-8.0/50, 50)[:,np.newaxis]
        surface = shdom.surface.wave_fresnel(real_refractive_index, imaginary_refractive_index, surface_wind_speed,
                                            ground_temperature=288.0, delx=0.02,dely=0.02)
        solver = shdom.solver.RTE(numerical_params=config,
                                                medium={'rayleigh': rayleigh[0.85]},
                                               source=shdom.source.solar(-0.707, 0.0, solarflux=1.0),
                                               surface=surface,
                                                num_stokes=3,
                                                name=None,
                                                atmosphere=None)

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver
        cls.fluxes = parse_shdom_output('../shdom_verification/brdf_W1f.out', comment='!')
        cls.radiances = parse_shdom_output('../shdom_verification/brdf_W1r.out', comment='!')

    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_up(self):
        self.assertTrue(np.allclose(self.fluxes[:,2],self.solver.fluxes.flux_up[:,0,0].data, atol=4e-6))

    def test_radiance(self):
        self.assertTrue(np.allclose(self.radiances[:,2], self.integrated_rays.I.data, atol=9e-6))

    def test_Q(self):
        self.assertTrue(np.allclose(self.radiances[:,3], self.integrated_rays.Q.data, atol=9e-6))

    def test_U(self):
        self.assertTrue(np.allclose(self.radiances[:,4], self.integrated_rays.U.data, atol=9e-6))


class Verify_Diner_Surfaces(TestCase):
    @classmethod
    def setUpClass(cls):

        sensor, rayleigh, config = get_basic_state_for_surface()
        a = np.ones(50)[:,np.newaxis]*0.2
        k = np.ones(50)[:,np.newaxis]*0.8
        b = np.ones(50)[:,np.newaxis]*0.3
        zeta = np.linspace(0., 1.0 - 1.0/50, 50)[:, np.newaxis]
        sigma = np.ones(50)[:,np.newaxis] * -1.0

        surface = shdom.surface.diner(a, k, b, zeta, sigma, ground_temperature=288.0, delx=0.02,dely=0.02)
        solver = shdom.solver.RTE(
                    numerical_params=config,
                    medium={'rayleigh': rayleigh[0.85]},
                    source=shdom.source.solar(-0.707, 0.0, solarflux=1.0),
                    surface=surface,
                    num_stokes=3,
                    name=None,
                    atmosphere=None
                    )

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver
        cls.fluxes = parse_shdom_output('../shdom_verification/brdf_D1f.out', comment='!')
        cls.radiances = parse_shdom_output('../shdom_verification/brdf_D1r.out', comment='!')

    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_up(self):
        self.assertTrue(np.allclose(self.fluxes[:,2],self.solver.fluxes.flux_up[:,0,0].data, atol=4e-6))

    def test_radiance(self):
        self.assertTrue(np.allclose(self.radiances[:,2], self.integrated_rays.I.data, atol=9e-6))

    def test_Q(self):
        self.assertTrue(np.allclose(self.radiances[:,3], self.integrated_rays.Q.data, atol=9e-6))

    def test_U(self):
        self.assertTrue(np.allclose(self.radiances[:,4], self.integrated_rays.U.data, atol=9e-6))
