import os
from unittest import TestCase
from collections import OrderedDict
import numpy as np
import xarray as xr
import pyshdom

def make_test_cloud():
    """
    Create test cloud and atmosphere Datasets using random values

    Returns
    -------
    scatterer : xarray.Dataset
        Dataset containing density [] and reff [micron] with values between [0, 1)
        along x, y, z grid [km]
    atmosphere : xarray.Dataset
        Dataset containing temperature profile [K] with values between [0, 200) along z [km]
    """
    np.random.seed(1)
    scatterer = xr.Dataset(
                    data_vars={
                     'density':(['x','y','z'],np.random.random(size = (20,26,32))),
                     'reff':(['x','y','z'],np.random.random(size = (20,26,32))),
                    },
        coords = {
            'x':np.linspace(0.0,5.5,20),
            'y': np.linspace(0.0,2.6,26),
            'z': np.linspace(0.3,2.6,32)
        }
    )
    atmosphere = xr.Dataset(
                data_vars = {
                    'temperature': ('z', (np.random.random(size=11)*100)+200)
                },
        coords = {
            'z': np.linspace(0.3,4.0,11)
        }
    )
    return scatterer, atmosphere

class Load_from_csv_test(TestCase):
    def test_load_from_csv(self, path= '../data/synthetic_cloud_fields/jpl_les/rico32x37x26.txt'):
        scatterer = pyshdom.util.load_from_csv(path, density='lwc', origin=(0.0,0.0))

class Microphysics_load_from_csv_test(TestCase):
    def setUp(self):
        scatterer = pyshdom.util.load_from_csv('../data/synthetic_cloud_fields/jpl_les/rico32x37x26.txt', density='lwc', origin=(0.0,0.0))
        x = scatterer.x.data
        y = scatterer.y.data
        z = scatterer.z.data
        rte_grid = pyshdom.grid.make_grid(x[1]-x[0], x.size,y[1]-y[0], y.size, z)
        self.droplets = pyshdom.grid.resample_onto_grid(rte_grid, scatterer)
        self.droplets['veff'] = (self.droplets.reff.dims, np.full_like(self.droplets.reff.data, fill_value=0.1))
    def test_load_from_csv_lwc(self):
            self.assertAlmostEqual(self.droplets.density.data[self.droplets.density.data>0].mean(), 0.265432, places=6)
    def test_load_from_csv_reff(self):
        self.assertAlmostEqual(self.droplets.reff.data[self.droplets.reff.data>0].mean(), 15.972615, places=6)
    def test_load_from_csv_veff(self):
        self.assertAlmostEqual(self.droplets.veff.data[self.droplets.veff.data>0].mean(), 0.100000, places=6)

class Save_2parameter_lwc(TestCase):
    def test_save_2parameter_lwc_atmos(self, path='test_2paramlwc_atmos.lwc'):
        scatterer, atmosphere = make_test_cloud()
        pyshdom.util.to_2parameter_lwc_file(path, scatterer, atmosphere=atmosphere)
        os.remove(path)
    def test_save_2parameter_lwc_filltemp(self, path='test_2paramlwc_filltemp.lwc'):
        scatterer, atmosphere = make_test_cloud()
        pyshdom.util.to_2parameter_lwc_file(path, scatterer)
        os.remove(path)

class Load_from_parameter_lwc_file(TestCase):
    def setUp(self):
        np.random.seed(1)
        scatterer, atmosphere = make_test_cloud()
        pyshdom.util.to_2parameter_lwc_file('test_2paramlwc_atmos.lwc', scatterer, atmosphere=atmosphere)
        pyshdom.util.to_2parameter_lwc_file('test_2paramlwc_filltemp.lwc', scatterer)
    def test_load_2parameter_atmos(self):
        self.scatterer_atmos = pyshdom.util.load_2parameter_lwc_file('test_2paramlwc_atmos.lwc')
        os.remove('test_2paramlwc_atmos.lwc')
    def test_load_2parameter_filltemp(self):
        self.scatterer = pyshdom.util.load_2parameter_lwc_file('test_2paramlwc_filltemp.lwc')
        os.remove('test_2paramlwc_filltemp.lwc')

class Microphysics_2parameter_lwc(TestCase):

    def setUp(self):
        scatterer, atmosphere = make_test_cloud()
        self.scatterer_original = scatterer
        self.atmosphere = atmosphere

        self.scatterer_original['temperature'] = self.atmosphere.interp({'z':self.scatterer_original.z}).temperature

        pyshdom.util.to_2parameter_lwc_file('test_2paramlwc_atmos.lwc', scatterer, atmosphere=atmosphere)
        pyshdom.util.to_2parameter_lwc_file('test_2paramlwc_filltemp.lwc', scatterer)

        self.scatterer_atmos = pyshdom.util.load_2parameter_lwc_file('test_2paramlwc_atmos.lwc')
        self.scatterer = pyshdom.util.load_2parameter_lwc_file('test_2paramlwc_filltemp.lwc')

    def test_load_2parameter_atmos_lwc(self):
        self.assertAlmostEqual(self.scatterer_atmos.density.data.mean(), self.scatterer_original.density.data.mean(),places=6)
    def test_load_2parameter_atmos_reff(self):
        self.assertAlmostEqual(self.scatterer_atmos.reff.data.mean(), self.scatterer_original.reff.data.mean(),places=4)
    def test_load_2parameter_atmos_temp(self):
        self.assertAlmostEqual(self.scatterer_atmos.temperature.data.mean(), self.scatterer_original.temperature.data.mean(),places=2)

    def test_load_2parameter_lwc(self):
        self.assertAlmostEqual(self.scatterer.density.data.mean(), self.scatterer_original.density.data.mean(),places=6)
    def test_load_2parameter_reff(self):
        self.assertAlmostEqual(self.scatterer.reff.data.mean(), self.scatterer_original.reff.data.mean(),places=4)
    def test_load_2parameter_temp(self):
        self.assertAlmostEqual(self.scatterer.temperature.data.mean(), 280.0,places=6)

    def tearDown(self):
        os.remove('test_2paramlwc_atmos.lwc')
        os.remove('test_2paramlwc_filltemp.lwc')

class Merge_two_z_coordinates(TestCase):
    def setUp(self):
        np.random.seed(1)
        #inputs must be sorted numpy arrays and not have repeated values (strictly increasing).
        z1 = np.sort(np.random.random(size=np.random.randint(1,high=1000))*np.random.randint(1,high=10))
        z2 = np.sort(np.random.random(size=np.random.randint(1,high=1000))*np.random.randint(1,high=10))
        assert np.unique(z1).size == z1.size
        assert np.unique(z2).size == z2.size
        combine = pyshdom.grid.merge_two_z_coordinates(z1,z2)
        self.z1 = z1
        self.z2 = z2
        self.combine = combine
    def test_equivalent(self):
        self.assertTrue(np.all(pyshdom.grid.merge_two_z_coordinates(self.z1,self.z2)==pyshdom.grid.merge_two_z_coordinates(self.z2,self.z1)))
    def test_combine_grid_min(self):
        self.assertEqual(min(self.z1.min(), self.z2.min()), self.combine.min())
    def test_combine_grid_max(self):
        self.assertEqual(max(self.z1.max(), self.z2.max()), self.combine.max())
    def test_combine_grid_overlap(self):
        #tests that the number of points in the overlap region is equal
        #to the maximum of the number of points in either of the two grids in the overlap region.
        combine_overlap = len(self.combine[np.where((self.combine <= min(self.z1.max(),self.z2.max())) & (self.combine >= max(self.z1.min(), self.z2.min())))])
        z1_overlap = len(self.z1[np.where((self.z1 <= min(self.z1.max(),self.z2.max())) & (self.z1 >= max(self.z1.min(), self.z2.min())))])
        z2_overlap = len(self.z2[np.where((self.z2 <= min(self.z1.max(),self.z2.max())) & (self.z2 >= max(self.z1.min(), self.z2.min())))])
        self.assertEqual(combine_overlap, max(z1_overlap, z2_overlap))
    def test_combine_grid_unique(self):
        self.assertEqual(np.unique(self.combine).size,self.combine.size)
    def test_combine_grid_sorted(self):
        self.assertTrue(np.allclose(np.sort(self.combine), self.combine))

class Combine_z_coordinates(TestCase):
    def setUp(self):
        np.random.seed(1)
        self.zlist = []
        for i in range(10):
            data = np.sort(np.random.random(size=np.random.randint(1,high=1000))*np.random.randint(1,high=10))
            if np.unique(data).size == data.size:
                dset = xr.Dataset(
                            coords={
                            'z':data
                            }
                )
            self.zlist.append(dset)
        self.combine = pyshdom.grid.combine_z_coordinates(self.zlist)
    def test_combine(self):
         combined = pyshdom.grid.combine_z_coordinates(self.zlist)
    #same tests should apply for Combine_z_coordinates as Merge_two_z_coordinates.
    def test_combine_grid_unique(self):
        self.assertEqual(np.unique(self.combine).size,self.combine.size)
    def test_combine_grid_sorted(self):
        self.assertTrue(np.allclose(np.sort(self.combine), self.combine))
    def test_combine_grid_min(self):
        minimum = min([dset.z.min() for dset in self.zlist])
        self.assertEqual(minimum, self.combine.min())
    def test_combine_grid_max(self):
        maximum = max([dset.z.max() for dset in self.zlist])
        self.assertEqual(maximum, self.combine.max())

class Resample_onto_grid(TestCase):

    def setUp(self):
        cloud, atmosphere = make_test_cloud()
        self.cloud = cloud
        self.atmosphere = atmosphere
        self.rte_grid = pyshdom.grid.make_grid(2.5/29,30, 1.6/27 ,28,np.linspace(0.0,4.6,19))
    def test_resample_3d_onto_grid(self):
        pyshdom.grid.resample_onto_grid(self.rte_grid, self.cloud)
    def test_resample_1d_onto_grid(self):
        pyshdom.grid.resample_onto_grid(self.rte_grid, self.atmosphere)

class Resample_onto_grid_values(TestCase):
    def setUp(self):
        cloud, atmosphere = make_test_cloud()
        self.cloud = cloud
        self.atmosphere = atmosphere
        self.rte_grid = pyshdom.grid.make_grid(2.5/29,30, 1.6/27 ,28,np.linspace(0.0,4.6,19))
        self.grid_1d = pyshdom.grid.resample_onto_grid(self.rte_grid, atmosphere)
        self.grid_3d = pyshdom.grid.resample_onto_grid(self.rte_grid, cloud)

    def test_resample_3d_onto_grid_valid_density(self):
        self.assertTrue(np.bitwise_not(np.all(np.isnan(self.grid_3d.density.data))))
    def test_resample_3d_onto_grid_valid_reff(self):
        self.assertTrue(np.bitwise_not(np.all(np.isnan(self.grid_3d.reff.data))))
    def test_resample_1d_onto_grid_valid_temp(self):
        self.assertTrue(np.bitwise_not(np.all(np.isnan(self.grid_1d.temperature.data))))


class VerifySolver(TestCase):

    @classmethod
    def setUpClass(cls):
        #shdom_polarized = xr.open_dataset('data/shdomout_rico32x36x26w672_polarized.nc')

        cloud_scatterer = pyshdom.util.load_2parameter_lwc_file('data/rico32x36x26.txt', density='lwc')

        rte_grid = pyshdom.grid.make_grid(cloud_scatterer.x.data[1] -cloud_scatterer.x.data[0],cloud_scatterer.x.data.size,
                                   cloud_scatterer.y.data[1] -cloud_scatterer.y.data[0],cloud_scatterer.y.data.size,
                                   np.append(np.array([0.0]),cloud_scatterer.z.data))

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, cloud_scatterer)

        #define any necessary variables for microphysics here.
        size_distribution_function = pyshdom.size_distribution.gamma
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

        sensor = pyshdom.sensor.make_sensor_dataset(x2.ravel(),y2.ravel(),z2.ravel(),mu2.ravel(),np.deg2rad(phi2.ravel()),['I'],
                                                 0.672, fill_ray_variables=True)
        Sensordict = pyshdom.util.SensorsDict()
        Sensordict.add_sensor('Sensor0', sensor)

        config = pyshdom.configuration.get_config('../default_config.json')
        config['split_accuracy'] = 0.1
        config['spherical_harmonics_accuracy'] = 0.01
        config['num_mu_bins'] = 8
        config['num_phi_bins'] = 16
        config['solution_accuracy'] = 1e-4


        solvers = pyshdom.util.SolversDict()

        for wavelength in wavelengths:

            #mie table and medium doesn't matter here as it is overwritten by property file.
            #just need it to initialize a solver.
            mie_mono_table = pyshdom.mie.get_mono_table('Water',(wavelength,wavelength), relative_dir='../mie_tables')
            cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
                                                                    mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function,particle_density=1.0,
                                reff=[5.0,20.0,50,'linear','micron'],
                                alpha=[6.9,7.1,2,'linear','unitless'],
                                )
            poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            cloud_optical_scatterer = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)

            solvers.add_solver(wavelength, pyshdom.solver.RTE(numerical_params=config,
                                            medium={'cloud':cloud_optical_scatterer},
                                           source=pyshdom.source.solar(-0.5, 0.0, solarflux=1.0),
                                           surface=pyshdom.surface.lambertian(albedo=0.05),
                                            num_stokes=3,
                                            name=None,
                                            atmosphere=None
                                           )
                                           )

        solver = solvers[0.672]
        solver._solve_prop(filename = 'data/rico32x36x26w672.prp')

        shdom_source = []
        with open('data/shdom_verification_source_out.out') as f:
            data = f.readlines()
            for line in data:
                if not '*' in line:
                    shdom_source.append(np.fromstring(line, sep=' '))
        shdom_source = np.array(shdom_source)

        cls.testing = solver._source[:,:solver._npts]
        cls.truth = shdom_source.T
        integrated_rays = solver.integrate_to_sensor(sensor)
        cls.integrated_rays = integrated_rays

        radiances = []
        with open('data/rico32x36x26w672ar.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    radiances.append(np.fromstring(line, sep=' '))
        radiances = np.array(radiances)
        cls.radiances = radiances

    def test_solver(self):
        self.assertTrue(np.allclose(self.testing, self.truth))

    def test_radiance(self):
        self.assertTrue(np.allclose(self.integrated_rays.I.data, self.radiances[:,2].data, atol=3e-3))

    def test_Q(self):
        self.assertTrue(np.allclose(self.integrated_rays.Q.data, self.radiances[:,3].data, atol=2e-4))

    def test_U(self):
        self.assertTrue(np.allclose(self.integrated_rays.U.data, self.radiances[:,4].data, atol=7e-5))

def get_basic_state_for_surface():
    config = pyshdom.configuration.get_config('../default_config.json')
    config['split_accuracy'] = 0.001
    config['spherical_harmonics_accuracy'] = 0.0
    config['num_mu_bins'] = 16
    config['num_phi_bins'] = 32
    config['solution_accuracy'] = 1e-5
    config['x_boundary_condition'] = 'periodic'
    config['y_boundary_condition'] = 'periodic'
    config['ip_flag'] = 3

    rte_grid = pyshdom.grid.make_grid(0.02, 50, 0.02, 1,
                               np.array([0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0]))

    atmosphere = xr.Dataset(
                        data_vars = {
                            'temperature': ('z', np.array([288.0,269.0,249.0,230.0,217.0,217.0,217.0,
                                                           218.0,221.0,224.0,227.0])),
                            'pressure': ('z', np.ones(rte_grid.z.size)*1013.25)
                        },
        coords = {
            'z': rte_grid.z.data
        }
    )
    wavelengths = np.atleast_1d(0.85)
    rayleigh= pyshdom.rayleigh.to_grid(wavelengths,atmosphere,rte_grid)
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

    sensor = pyshdom.sensor.make_sensor_dataset(x2.ravel(),y2.ravel(),z2.ravel(),mu2.ravel(),np.deg2rad(phi2.ravel()),['I'],
                                             0.85, fill_ray_variables=True)
    return sensor, rayleigh, config



class Verify_Lambertian_Surfaces(TestCase):
    @classmethod
    def setUpClass(cls):

        sensor, rayleigh, config = get_basic_state_for_surface()
        variable_lambertian = pyshdom.surface.lambertian(albedo=np.linspace(0.0,0.3- (0.3-0.0)/50.0,50)[:,np.newaxis],
                                                       ground_temperature=288.0,delx=0.02,dely=0.04)
        # variable_ocean = pyshdom.surface.ocean_unpolarized(surface_wind_speed = np.linspace(4.0, 12.0-(12.0-4.0)/50.0,50)[:,np.newaxis],
        #                                                 pigmentation = np.zeros((50,1)), ground_temperature=288.0,delx=0.02,dely=0.04)
        solver = pyshdom.solver.RTE(numerical_params=config,
                                        medium={'rayleigh': rayleigh[0.85]},
                                       source=pyshdom.source.solar(-0.707, 0.0, solarflux=1.0),
                                       surface=variable_lambertian,
                                        num_stokes=1,
                                        name=None,
                                        atmosphere=None)

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver

        fluxes = []
        with open('data/brdf_L1f.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    fluxes.append(np.fromstring(line, sep=' '))
        fluxes = np.array(fluxes)
        cls.fluxes = fluxes

        radiances = []
        with open('data/brdf_L1r.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    radiances.append(np.fromstring(line, sep=' '))
        radiances = np.array(radiances)
        cls.radiances = radiances

    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,2],self.solver.fluxes.flux_up[:,0,0].data, atol=4e-6))

    def test_radiance(self):
        self.assertTrue(np.allclose(self.radiances[:,-1], self.integrated_rays.I.data, atol=9e-6))


class Verify_Ocean_Unpolarized_Surfaces(TestCase):

    @classmethod
    def setUpClass(cls):

        sensor, rayleigh, config = get_basic_state_for_surface()
        variable_ocean = pyshdom.surface.ocean_unpolarized(surface_wind_speed = np.linspace(4.0, 12.0-(12.0-4.0)/50.0,50)[:,np.newaxis],
                                                        pigmentation = np.zeros((50,1)), ground_temperature=288.0,delx=0.02,dely=0.04)
        solver = pyshdom.solver.RTE(numerical_params=config,
                                        medium={'rayleigh': rayleigh[0.85]},
                                       source=pyshdom.source.solar(-0.707, 0.0, solarflux=1.0),
                                       surface=variable_ocean,
                                        num_stokes=1,
                                        name=None,
                                        atmosphere=None)

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver

        fluxes = []
        with open('data/brdf_O1f.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    fluxes.append(np.fromstring(line, sep=' '))
        fluxes = np.array(fluxes)
        cls.fluxes = fluxes

        radiances = []
        with open('data/brdf_O1r.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    radiances.append(np.fromstring(line, sep=' '))
        radiances = np.array(radiances)
        cls.radiances = radiances


    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
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
        surface = pyshdom.surface.RPV_unpolarized(rho0,k,theta, ground_temperature=288.0, delx=0.02,dely=0.02)
        solver = pyshdom.solver.RTE(numerical_params=config,
                                                medium={'rayleigh': rayleigh[0.85]},
                                               source=pyshdom.source.solar(-0.707, 0.0, solarflux=1.0),
                                               surface=surface,
                                                num_stokes=1,
                                                name=None,
                                                atmosphere=None)

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver

        fluxes = []
        with open('data/brdf_R1f.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    fluxes.append(np.fromstring(line, sep=' '))
        fluxes = np.array(fluxes)
        cls.fluxes = fluxes

        radiances = []
        with open('data/brdf_R1r.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    radiances.append(np.fromstring(line, sep=' '))
        radiances = np.array(radiances)
        cls.radiances = radiances

    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
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
        surface = pyshdom.surface.wave_fresnel(real_refractive_index, imaginary_refractive_index, surface_wind_speed,
                                            ground_temperature=288.0, delx=0.02,dely=0.02)
        solver = pyshdom.solver.RTE(numerical_params=config,
                                                medium={'rayleigh': rayleigh[0.85]},
                                               source=pyshdom.source.solar(-0.707, 0.0, solarflux=1.0),
                                               surface=surface,
                                                num_stokes=3,
                                                name=None,
                                                atmosphere=None)

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver

        fluxes = []
        with open('data/brdf_W1f.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    fluxes.append(np.fromstring(line, sep=' '))
        fluxes = np.array(fluxes)
        cls.fluxes = fluxes

        radiances = []
        with open('data/brdf_W1r.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    radiances.append(np.fromstring(line, sep=' '))
        radiances = np.array(radiances)
        cls.radiances = radiances

    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
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
        zeta = np.linspace(0., 1.0 - 1.0/50, 50)[:,np.newaxis]
        sigma = np.ones(50)[:,np.newaxis]*-1.0

        surface = pyshdom.surface.diner(a,k,b,zeta,sigma, ground_temperature=288.0, delx=0.02,dely=0.02)
        solver = pyshdom.solver.RTE(numerical_params=config,
                                                        medium={'rayleigh': rayleigh[0.85]},
                                                       source=pyshdom.source.solar(-0.707, 0.0, solarflux=1.0),
                                                       surface=surface,
                                                        num_stokes=3,
                                                        name=None,
                                                        atmosphere=None)

        solver.solve(maxiter=100, verbose=False)
        integrated_rays = solver.integrate_to_sensor(sensor)

        cls.integrated_rays = integrated_rays
        cls.solver = solver

        fluxes = []
        with open('data/brdf_D1f.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    fluxes.append(np.fromstring(line, sep=' '))
        fluxes = np.array(fluxes)
        cls.fluxes = fluxes

        radiances = []
        with open('data/brdf_D1r.out') as f:
            data = f.readlines()
            for line in data:
                if not '!' in line:
                    radiances.append(np.fromstring(line, sep=' '))
        radiances = np.array(radiances)
        cls.radiances = radiances

    def test_flux_direct(self):
        self.assertTrue(np.allclose(self.fluxes[:,4],self.solver.fluxes.flux_direct[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,3],self.solver.fluxes.flux_down[:,0,0].data, atol=4e-6))

    def test_flux_down(self):
        self.assertTrue(np.allclose(self.fluxes[:,2],self.solver.fluxes.flux_up[:,0,0].data, atol=4e-6))

    def test_radiance(self):
        self.assertTrue(np.allclose(self.radiances[:,2], self.integrated_rays.I.data, atol=9e-6))

    def test_Q(self):
        self.assertTrue(np.allclose(self.radiances[:,3], self.integrated_rays.Q.data, atol=9e-6))

    def test_U(self):
        self.assertTrue(np.allclose(self.radiances[:,4], self.integrated_rays.U.data, atol=9e-6))

class ParallelizationNoSubpixelRays(TestCase):
    @classmethod
    def setUpClass(cls):

        reff = 0.5
        ext = 20.0
        rte_grid = pyshdom.grid.make_grid(0.05,13,
                                   0.05,13,
                                   np.linspace(0.1,0.7,13))

        rte_grid['density'] = (['x','y','z'], np.ones((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)))
        rte_grid['reff'] = (['x','y','z'], np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+reff)
        rte_grid['veff'] = (['x','y','z'] ,np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+0.1)

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, rte_grid)

        #define any necessary variables for microphysics here.
        size_distribution_function = pyshdom.size_distribution.gamma

        #define sensors.
        Sensordict = pyshdom.util.SensorsDict()
        Sensordict2 = pyshdom.util.SensorsDict()
        misr_list = []
        sensor_zenith_list = [75.0,60.0,45.6,26.1]*2 + [0.0]
        sensor_azimuth_list = [90]*4 + [-90]*4 +[0.0]
        wavelengths = [0.86,0.86,0.86,1.38,1.38,2.2,2.2,3.4,3.4]
        for zenith,azimuth,wavelength in zip(sensor_zenith_list,sensor_azimuth_list,wavelengths):
            Sensordict.add_sensor('MISR',
                            pyshdom.sensor.orthographic_projection(wavelength, cloud_scatterer_on_rte_grid,0.04,0.04, azimuth, zenith,
                                                     altitude='TOA', stokes=['I'])
                                 )
            Sensordict2.add_sensor('MISR',
                            pyshdom.sensor.orthographic_projection(wavelength, cloud_scatterer_on_rte_grid,0.04,0.04, azimuth, zenith,
                                                     altitude='TOA', stokes=['I'])
            )

        #Define the RTE solvers needed to model the measurements and
        #calculate optical properties.
        wavelengths = Sensordict.get_unique_solvers()

        cloud_poly_tables = OrderedDict()
        solvers = pyshdom.util.SolversDict()

        for wavelength in wavelengths:

            #optical properties from mie calculations.
            mie_mono_table = pyshdom.mie.get_mono_table('Water',(wavelength,wavelength),
                                                      max_integration_radius=10.0,
                                                      minimum_effective_radius=0.1,
                                                      relative_dir='../mie_tables',
                                                      verbose=False)
            #mie_mono_table.to_netcdf('./mie_tables/MieTable_{}.nc'.format(reff))
            cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
                                                                    mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function,particle_density=1.0,
                                reff=[0.2,1.0,10,'logarithmic','micron'],
                                veff=[0.09,0.11,12,'linear','unitless'],
                                )
            poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            optical_properties = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                             exact_table=False)
            optical_properties['ssalb'][:,:,:] = 1.0
            extinction = np.zeros(optical_properties.extinction.shape)
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
            solver = pyshdom.solver.RTE(numerical_params=config,
                                            medium={'cloud': optical_properties},
                                           source=pyshdom.source.solar(-1*np.cos(np.deg2rad(60.0)),0.0,solarflux=1.0),
                                           surface=pyshdom.surface.ocean_unpolarized(surface_wind_speed=10.0,
                                                                                  pigmentation=0.0),
                                            num_stokes=1,
                                            name=None
                                           )
            solvers.add_solver(wavelength,solver)
            pyshdom.util.get_measurements(solvers, Sensordict, maxiter=100, n_jobs=8, verbose=False)
            pyshdom.util.get_measurements(solvers, Sensordict2, maxiter=100, n_jobs=1, verbose=False)
            #Sensordict['MISR']['sensor_list'][0].to_netcdf('data/RenderedSensorReference_nosubpixel.nc')

            cls.solvers = solvers
            cls.Sensordict = Sensordict
            cls.Sensordict2 = Sensordict2

    def test_radiance(self):
        self.assertTrue(all([np.allclose(self.Sensordict2['MISR']['sensor_list'][i].I, self.Sensordict['MISR']['sensor_list'][i].I) for i in range(9)]))

    def test_no_subpixel_args_reference(self):
        test = xr.open_dataset('data/RenderedSensorReference_nosubpixel.nc')
        self.assertTrue(test.equals(self.Sensordict['MISR']['sensor_list'][0]))

class ParallelizationSubpixelRays(TestCase):
    @classmethod
    def setUpClass(cls):

        reff = 0.5
        ext = 20.0
        rte_grid = pyshdom.grid.make_grid(0.05,13,
                                   0.05,13,
                                   np.linspace(0.1,0.7,13))

        rte_grid['density'] = (['x','y','z'], np.ones((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)))
        rte_grid['reff'] = (['x','y','z'], np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+reff)
        rte_grid['veff'] = (['x','y','z'] ,np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+0.1)

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, rte_grid)

        #define any necessary variables for microphysics here.
        size_distribution_function = pyshdom.size_distribution.gamma

        #define sensors.
        Sensordict = pyshdom.util.SensorsDict()
        Sensordict2 = pyshdom.util.SensorsDict()
        misr_list = []
        sensor_zenith_list = [75.0,60.0,45.6,26.1]*2 + [0.0]
        sensor_azimuth_list = [90]*4 + [-90]*4 +[0.0]
        wavelengths = [0.86,0.86,0.86,1.38,1.38,2.2,2.2,3.4,3.4]
        for zenith,azimuth,wavelength in zip(sensor_zenith_list,sensor_azimuth_list,wavelengths):
            Sensordict.add_sensor('MISR',
                            pyshdom.sensor.orthographic_projection(wavelength, cloud_scatterer_on_rte_grid,0.04,0.04, azimuth, zenith,
                                                     altitude='TOA', stokes=['I'],
                                                     sub_pixel_ray_args={'method': pyshdom.sensor.gaussian,
                                                     'degree':4})
                                 )
            Sensordict2.add_sensor('MISR',
                            pyshdom.sensor.orthographic_projection(wavelength, cloud_scatterer_on_rte_grid,0.04,0.04, azimuth, zenith,
                                                     altitude='TOA', stokes=['I'],
                                                  sub_pixel_ray_args={'method': pyshdom.sensor.gaussian,
                                                  'degree':4})
            )

        #Define the RTE solvers needed to model the measurements and
        #calculate optical properties.
        wavelengths = Sensordict.get_unique_solvers()

        cloud_poly_tables = OrderedDict()
        solvers = pyshdom.util.SolversDict()

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
                                reff=[0.2,1.0,10,'logarithmic','micron'],
                                veff=[0.09,0.11,12,'linear','unitless'],
                                )
            poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            optical_properties = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                             exact_table=False)
            optical_properties['ssalb'][:,:,:] = 1.0
            extinction = np.zeros(optical_properties.extinction.shape)
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
            solver = pyshdom.solver.RTE(numerical_params=config,
                                            medium={'cloud': optical_properties},
                                           source=pyshdom.source.solar(-1*np.cos(np.deg2rad(60.0)),0.0,solarflux=1.0),
                                           surface=pyshdom.surface.ocean_unpolarized(surface_wind_speed=10.0,
                                                                                  pigmentation=0.0),
                                            num_stokes=1,
                                            name=None
                                           )
            solvers.add_solver(wavelength,solver)
            pyshdom.util.get_measurements(solvers, Sensordict, maxiter=100, n_jobs=8, verbose=False)
            pyshdom.util.get_measurements(solvers, Sensordict2, maxiter=100, n_jobs=1, verbose=False)

            #Sensordict['MISR']['sensor_list'][0].to_netcdf('data/RenderedSensorReference_subpixelargs.nc')
            cls.solvers = solvers
            cls.Sensordict = Sensordict
            cls.Sensordict2 = Sensordict2

    def test_subpixel_args(self):
        self.assertTrue(all([np.allclose(self.Sensordict2['MISR']['sensor_list'][i].I, self.Sensordict['MISR']['sensor_list'][i].I) for i in range(9)]))

    def test_subpixel_args_reference(self):
        test = xr.open_dataset('data/RenderedSensorReference_subpixelargs.nc')
        self.assertTrue(test.equals(self.Sensordict['MISR']['sensor_list'][0]))

class MicrophysicalDerivatives(TestCase):

    @classmethod
    def setUpClass(cls):

        reff = 0.5
        ext = 20.0
        rte_grid = pyshdom.grid.make_grid(0.05,13,
                                   0.05,13,
                                   np.linspace(0.1,0.7,13))

        rte_grid['density'] = (['x','y','z'], np.ones((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)))
        rte_grid['reff'] = (['x','y','z'], np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+reff)
        rte_grid['veff'] = (['x','y','z'] ,np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+0.1)

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, rte_grid)

        #define any necessary variables for microphysics here.
        size_distribution_function = pyshdom.size_distribution.gamma

        #define sensors.
        Sensordict = pyshdom.util.SensorsDict()
        misr_list = []
        sensor_zenith_list = [75.0,60.0,45.6,26.1]*2 + [0.0]
        sensor_azimuth_list = [90]*4 + [-90]*4 +[0.0]
        wavelengths = [0.86,0.86,0.86,1.38,1.38,2.2,2.2,3.4,3.4]
        for zenith,azimuth,wavelength in zip(sensor_zenith_list,sensor_azimuth_list,wavelengths):
            Sensordict.add_sensor('MISR',
                            pyshdom.sensor.orthographic_projection(wavelength, cloud_scatterer_on_rte_grid,0.04,0.04, azimuth, zenith,
                                                     altitude='TOA', stokes=['I'])
                                 )

        #Define the RTE solvers needed to model the measurements and
        #calculate optical properties.
        wavelengths = Sensordict.get_unique_solvers()

        cloud_poly_tables = OrderedDict()
        solvers = pyshdom.util.SolversDict()


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
                                reff=[0.2,1.0,10,'logarithmic','micron'],
                                veff=[0.09,0.11,12,'linear','unitless'],
                                )
            poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            optical_properties = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                             exact_table=False)
            optical_properties['ssalb'][:,:,:] = 1.0
            extinction = np.zeros(optical_properties.extinction.shape)
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
            solver = pyshdom.solver.RTE(numerical_params=config,
                                            medium={'cloud': optical_properties},
                                           source=pyshdom.source.solar(-1*np.cos(np.deg2rad(60.0)),0.0,solarflux=1.0),
                                           surface=pyshdom.surface.ocean_unpolarized(surface_wind_speed=10.0,
                                                                                  pigmentation=0.0),
                                            num_stokes=1,
                                            name=None
                                           )

            solvers.add_solver(wavelength,solver)

        solvers.parallel_solve(maxiter=100,verbose=False)
        cls.solvers=solvers
        cls.cloud_poly_tables = cloud_poly_tables

    def test_ssalb(self):

        unknown_scatterers = pyshdom.util.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['ssalb'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        self.solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method, unknown_scatterers.table_data)
        solvers = self.solvers
        self.assertTrue(all([all([np.all(solvers[key]._dalb==1.0) for key in solvers]),
            all([np.all(solvers[key]._dext==0.0) for key in solvers]),
            all([np.all(solvers[key]._dleg==0.0)for key in solvers]),
            all([np.all(solvers[key]._dphasetab==0.0)for key in solvers])]))

    def test_extinction(self):

        unknown_scatterers = pyshdom.util.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['extinction'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers = self.solvers
        solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method, unknown_scatterers.table_data)

        self.assertTrue(all([all([np.all(solvers[key]._dalb==0.0) for key in solvers]),
            all([np.all(solvers[key]._dext==1.0) for key in solvers]),
            all([np.all(solvers[key]._dleg==0.0)for key in solvers]),
            all([np.all(solvers[key]._dphasetab==0.0)for key in solvers])]))

    def test_density(self):

        unknown_scatterers = pyshdom.util.UnknownScatterers()
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

        unknown_scatterers = pyshdom.util.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['legendre_0_10'],self.cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers = self.solvers
        solvers.add_microphysical_partial_derivatives(unknown_scatterers.table_to_grid_method, unknown_scatterers.table_data)

        self.assertTrue(all([all([np.all(solvers[key]._dalb==0.0) for key in solvers]),
            all([np.all(solvers[key]._dext==0.0) for key in solvers]),
                all([np.all(solvers[key]._dleg[0,10]==1.0/(2*10.0 + 1.0))for key in solvers]),
                ]))

class VerifyJacobian(TestCase):

    @classmethod
    def setUpClass(cls):

        ext = 0.1
        reff=10.0

        rte_grid = pyshdom.grid.make_grid(0.05,3,
                           0.05,3,
                           np.arange(0.1,0.25,0.05))

        rte_grid['density'] = (['x','y','z'], np.ones((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)))
        rte_grid['reff'] = (['x','y','z'], np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+reff)
        rte_grid['veff'] = (['x','y','z'] ,np.zeros((rte_grid.x.size, rte_grid.y.size, rte_grid.z.size))+0.1)

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, rte_grid)

        #define any necessary variables for microphysics here.
        size_distribution_function = pyshdom.size_distribution.gamma

        #define sensors.
        Sensordict = pyshdom.util.SensorsDict()

        sensor = pyshdom.sensor.make_sensor_dataset(np.array([0.05]),np.array([0.05]),
                                                  np.array([0.7]),np.array([1.0]),
                                                  np.array([0.0]),
                                                 wavelength=0.86,stokes=['I'],fill_ray_variables=True)
        Sensordict.add_sensor('MISR', sensor)

        #Define the RTE solvers needed to model the measurements and
        #calculate optical properties.
        wavelengths = Sensordict.get_unique_solvers()

        cloud_poly_tables = OrderedDict()
        solvers = pyshdom.util.SolversDict()

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
                                reff=[9.0,11.0,100,'logarithmic','micron'],
                                veff=[0.09,0.11,12,'linear','unitless'],
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
            solver = pyshdom.solver.RTE(numerical_params=config,
                                            medium={'cloud': optical_properties},
                                           source=pyshdom.source.solar(-1,0.0,solarflux=1.0),
                                           surface=pyshdom.surface.lambertian(albedo=0.0),
                                            num_stokes=1,
                                            name=None
                                           )

            solvers.add_solver(wavelength,solver)
        pyshdom.util.get_measurements(solvers, Sensordict, maxiter=100, n_jobs=8, verbose=False)

        unknown_scatterers = pyshdom.util.UnknownScatterers()
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
        self.assertAlmostEqual(self.jacobian_exact['jacobian_0.860'][0,0,0,0].data,0.00396336,places=5)
