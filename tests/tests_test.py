
import shdom
import os
import numpy as np
import xarray as xr
import pylab as py
from unittest import TestCase
import typing

def make_test_cloud():

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
    def test_load_from_csv(self, path= '../synthetic_cloud_fields/jpl_les/rico32x37x26.txt'):
        scatterer = shdom.grid.load_from_csv(path, density='lwc', origin=(0.0,0.0))

class Microphysics_load_from_csv_test(TestCase):
    def setUp(self):
        scatterer = shdom.grid.load_from_csv('../synthetic_cloud_fields/jpl_les/rico32x37x26.txt', density='lwc', origin=(0.0,0.0))
        x = scatterer.x.data
        y = scatterer.y.data
        z = scatterer.z.data
        rte_grid = shdom.grid.make_grid(x.min(), x.max(), x.size, y.min(), y.max(), y.size, z)
        self.droplets = shdom.grid.resample_onto_grid(rte_grid, scatterer)
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
        shdom.grid.to_2parameter_lwc_file(path, scatterer, atmosphere=atmosphere)
        os.remove(path)
    def test_save_2parameter_lwc_filltemp(self, path='test_2paramlwc_filltemp.lwc'):
        scatterer, atmosphere = make_test_cloud()
        shdom.grid.to_2parameter_lwc_file(path, scatterer)
        os.remove(path)

class Load_from_parameter_lwc_file(TestCase):
    def setUp(self):
        np.random.seed(1)
        scatterer, atmosphere = make_test_cloud()
        shdom.grid.to_2parameter_lwc_file('test_2paramlwc_atmos.lwc', scatterer, atmosphere=atmosphere)
        shdom.grid.to_2parameter_lwc_file('test_2paramlwc_filltemp.lwc', scatterer)
    def test_load_2parameter_atmos(self):
        self.scatterer_atmos = shdom.grid.load_2parameter_lwc_file('test_2paramlwc_atmos.lwc',origin=(0.0,0.0))
        os.remove('test_2paramlwc_atmos.lwc')
    def test_load_2parameter_filltemp(self):
        self.scatterer = shdom.grid.load_2parameter_lwc_file('test_2paramlwc_filltemp.lwc',origin=(0.0,0.0))
        os.remove('test_2paramlwc_filltemp.lwc')

class Microphysics_2parameter_lwc(TestCase):

    def setUp(self):
        scatterer, atmosphere = make_test_cloud()
        self.scatterer_original = scatterer
        self.atmosphere = atmosphere

        self.scatterer_original['temperature'] = self.atmosphere.interp({'z':self.scatterer_original.z}).temperature

        shdom.grid.to_2parameter_lwc_file('test_2paramlwc_atmos.lwc', scatterer, atmosphere=atmosphere)
        shdom.grid.to_2parameter_lwc_file('test_2paramlwc_filltemp.lwc', scatterer)

        self.scatterer_atmos = shdom.grid.load_2parameter_lwc_file('test_2paramlwc_atmos.lwc',origin=(0.1,0.1))
        self.scatterer = shdom.grid.load_2parameter_lwc_file('test_2paramlwc_filltemp.lwc',origin=(0.1,0.1))

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
        combine = shdom.grid.merge_two_z_coordinates(z1,z2)
        self.z1 = z1
        self.z2 = z2
        self.combine = combine
    def test_equivalent(self):
        self.assertTrue(np.all(shdom.grid.merge_two_z_coordinates(self.z1,self.z2)==shdom.grid.merge_two_z_coordinates(self.z2,self.z1)))
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
        self.combine = shdom.grid.combine_z_coordinates(self.zlist)
    def test_combine(self):
         combined = shdom.grid.combine_z_coordinates(self.zlist)
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
        self.rte_grid = shdom.grid.make_grid(0.0, 2.5,30, 0.0,1.6,28,np.linspace(0.0,4.6,19))
    def test_resample_3d_onto_grid(self):
        shdom.grid.resample_onto_grid(self.rte_grid, self.cloud)
    def test_resample_1d_onto_grid(self):
        shdom.grid.resample_onto_grid(self.rte_grid, self.atmosphere)

class Resample_onto_grid_values(TestCase):
    def setUp(self):
        cloud, atmosphere = make_test_cloud()
        self.cloud = cloud
        self.atmosphere = atmosphere
        self.rte_grid = shdom.grid.make_grid(0.0, 2.5,30, 0.0,1.6,28,np.linspace(0.0,4.6,19))
        self.grid_1d = shdom.grid.resample_onto_grid(self.rte_grid, atmosphere)
        self.grid_3d = shdom.grid.resample_onto_grid(self.rte_grid, cloud)

    def test_resample_3d_onto_grid_valid_density(self):
        self.assertTrue(np.bitwise_not(np.all(np.isnan(self.grid_3d.density.data))))
    def test_resample_3d_onto_grid_valid_reff(self):
        self.assertTrue(np.bitwise_not(np.all(np.isnan(self.grid_3d.reff.data))))
    def test_resample_1d_onto_grid_valid_temp(self):
        self.assertTrue(np.bitwise_not(np.all(np.isnan(self.grid_1d.temperature.data))))


class VerifySolver(TestCase):
    def setUp(self):
        shdom_polarized = xr.open_dataset('../shdom_verification/shdomout_rico32x36x26w672_polarized.nc')

        cloud_scatterer = shdom.grid.load_2parameter_lwc_file('../shdom_verification/rico32x36x26.txt',
                                                   density='lwc',origin=(0.0,0.0))

        #simple 'union' horizontal grid merging for 3D and 1D needs to be fixed.
        rte_grid = shdom.grid.make_grid(cloud_scatterer.x.data.min(),cloud_scatterer.x.data.max(),cloud_scatterer.x.data.size,
                                   cloud_scatterer.y.data.min(),cloud_scatterer.y.data.max(),cloud_scatterer.y.data.size,
                                   np.append(np.array([0.0]),cloud_scatterer.z.data))

        #resample the cloud onto the rte_grid
        cloud_scatterer_on_rte_grid = shdom.grid.resample_onto_grid(rte_grid, cloud_scatterer)

        #define any necessary variables for microphysics here.
        size_distribution_function = shdom.size_distribution.gamma
        #We choose a gamma size distribution and therefore need to define a 'veff' variable.
        cloud_scatterer_on_rte_grid['alpha'] = (cloud_scatterer_on_rte_grid.reff.dims,
                                               np.full_like(cloud_scatterer_on_rte_grid.reff.data, fill_value=7))

        wavelengths = np.atleast_1d(0.672)

        x = np.repeat(np.repeat(shdom_polarized.radiance_x.data[np.newaxis,:,np.newaxis],shdom_polarized.radiance_y.size, axis=2),
                     shdom_polarized.radiance_mu.data.size,axis=0).ravel()
        y = np.repeat(np.repeat(shdom_polarized.radiance_y.data[np.newaxis,np.newaxis,:],shdom_polarized.radiance_x.size, axis=1),
                     shdom_polarized.radiance_mu.data.size,axis=0).ravel()
        mu = np.repeat(np.repeat(shdom_polarized.radiance_mu.data[:,np.newaxis,np.newaxis],shdom_polarized.radiance_x.size, axis=0),
                       shdom_polarized.radiance_y.size, axis=-1).ravel()
        phi = np.repeat(np.repeat(shdom_polarized.radiance_phi.data[:,np.newaxis,np.newaxis],shdom_polarized.radiance_x.size, axis=-1),
                       shdom_polarized.radiance_y.size, axis=-1).ravel()
        sensor = shdom.sensor.make_sensor_dataset(
            x.astype(np.float32),
            y.astype(np.float32),
            np.full_like(x,fill_value=1.4).astype(np.float32),
            mu.astype(np.float32),
            np.deg2rad(phi).astype(np.float32),
            stokes=['I','Q','U'],wavelength = wavelengths[0], fill_ray_variables=True)
        Sensordict = shdom.organization.SensorsDict()
        Sensordict.add_sensor('Sensor0', sensor)

        config = shdom.configuration.get_config('../default_config.json')
        config['split_accuracy'] = 0.1
        config['spherical_harmonics_accuracy'] = 0.01
        config['num_mu_bins'] = 8
        config['num_phi_bins'] = 16
        config['solution_accuracy'] = 1e-4
        config['x_boundary_condition'] = 'periodic'
        config['y_boundary_condition'] = 'periodic'

        solvers = shdom.organization.SolversDict()

        for wavelength in wavelengths:
            #print('making mie_table. . . may take a while.')

            #mie table and medium doesn't matter here as it is overwritten by property file.
            #just need it to initialize a solver.
            mie_mono_table = shdom.mie.get_mono_table('Water',(wavelength,wavelength), relative_path='../mie_tables')
            cloud_size_distribution = shdom.size_distribution.get_size_distribution_grid(
                                                                    mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function,particle_density=1.0,
                                reff=[5.0,20.0,50,'linear','micron'],
                                alpha=[7,7.1,2,'linear','unitless'],
                                )
            poly_table = shdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
            cloud_optical_scatterer = shdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)

            solvers.add_solver(wavelength, shdom.solver.RTE(numerical_params=config,
                                            medium={'cloud':cloud_optical_scatterer},
                                           source=shdom.source.solar(-0.5, 0.0, solarflux=1.0),
                                           surface=shdom.surface.lambertian(albedo=0.05),
                                            num_stokes=3,
                                            name=None,
                                            atmosphere=None
                                           )
                                           )

        solver = solvers[0.672]
        solver._solve_prop(filename = '../shdom_verification/rico32x36x26w672.prp')

        shdom_source = []
        with open('../shdom_verification/shdom_verification_source_out.out') as f:
            data = f.readlines()
            for line in data:
                if not '*' in line:
                    shdom_source.append(np.fromstring(line, sep=' '))
        shdom_source = np.array(shdom_source)

        self.testing = solver._source[:,:solver._npts]
        self.truth = shdom_source.T

    def test_solver(self):
        self.assertTrue(np.allclose(self.testing, self.truth))



if __name__ == '__main__':
    import unittest
    unittest.main()

# cloud_scatterer = shdom.grid.load_from_csv('./synthetic_cloud_fields/jpl_les/rico32x37x26.txt',
#                                            density='lwc',origin=(0.5,0.1))

# shdom.grid.to_2parameter_lwc_file('test_2paramlwc.lwc', cloud_scatterer)
# cloud_scatterer2 = shdom.grid.load_2parameter_lwc_file('test_2paramlwc.lwc',origin=(0.1,0.1))

# print(np.allclose(cloud_scatterer, cloud_scatterer2, atol=5e-5,equal_nan=True))
# os.remove('test_2paramlwc.lwc')
