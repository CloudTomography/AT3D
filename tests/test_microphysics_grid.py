from unittest import TestCase
import os
import numpy as np
import xarray as xr
import at3d
import importlib

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
    def test_load_from_csv(self, path=os.path.join(importlib.resources.files('at3d'),'data/synthetic_cloud_fields/jpl_les/rico32x37x26.txt')):
        scatterer = at3d.util.load_from_csv(path, density='lwc', origin=(0.0,0.0))

class Microphysics_load_from_csv_test(TestCase):
    def setUp(self):
        self.droplets = at3d.util.load_from_csv(os.path.join(importlib.resources.files('at3d'),'data/synthetic_cloud_fields/jpl_les/rico32x37x26.txt'), density='lwc', origin=(0.0,0.0))
        self.droplets['veff'] = (self.droplets.reff.dims, np.full_like(self.droplets.reff.data, fill_value=0.1))
    def test_load_from_csv_lwc(self):
        self.assertAlmostEqual(self.droplets.density.data[self.droplets.density.data>0].mean(), 0.265432, places=6)
    def test_load_from_csv_reff(self):
        self.assertAlmostEqual(self.droplets.reff.data[self.droplets.reff.data>0].mean(), 15.713033, places=6)
    def test_load_from_csv_veff(self):
        self.assertAlmostEqual(self.droplets.veff.data[self.droplets.veff.data>0].mean(), 0.100000, places=6)


class Test_resample_onto_grid(TestCase):
    def setUp(self):
        scatterer = at3d.util.load_from_csv(os.path.join(importlib.resources.files('at3d'),'data/synthetic_cloud_fields/jpl_les/rico32x37x26.txt'), density='lwc', origin=(0.0,0.0))
        scatterer['veff'] = (scatterer.reff.dims, np.full_like(scatterer.reff.data, fill_value=0.1))
        x = scatterer.x.data
        y = scatterer.y.data
        z = scatterer.z.data
        rte_grid = at3d.grid.make_grid(x[1]-x[0], x.size,y[1]-y[0], y.size, z)
        self.droplets = at3d.grid.resample_onto_grid(scatterer, scatterer)
    def test_resampled_lwc(self):
        # Do not use 0 as the threshold for comparison as the interpolation routine can create some
        # very small non-negative values that would otherwise bias the comparison. Those changes
        # are unstable, and change depending on the version of xarray (which controls interpolation) that
        # is used.
        self.assertAlmostEqual(self.droplets.density.data[self.droplets.density.data > 1e-12].mean(), 0.265432, places=6)
    def test_resampled_reff(self):
        self.assertAlmostEqual(self.droplets.reff.data[self.droplets.reff.data>0].mean(), 15.972615, places=6)
    def test_resampled_veff(self):
        self.assertAlmostEqual(self.droplets.veff.data[self.droplets.veff.data>0].mean(), 0.100000, places=6)


class Save_2parameter_lwc(TestCase):
    def test_save_2parameter_lwc_atmos(self, path='test_2paramlwc_atmos.lwc'):
        scatterer, atmosphere = make_test_cloud()
        at3d.util.to_2parameter_lwc_file(path, scatterer, atmosphere=atmosphere)
        os.remove(path)
    def test_save_2parameter_lwc_filltemp(self, path='test_2paramlwc_filltemp.lwc'):
        scatterer, atmosphere = make_test_cloud()
        at3d.util.to_2parameter_lwc_file(path, scatterer)
        os.remove(path)

class Load_from_parameter_lwc_file(TestCase):
    def setUp(self):
        np.random.seed(1)
        scatterer, atmosphere = make_test_cloud()
        at3d.util.to_2parameter_lwc_file('test_2paramlwc_atmos.lwc', scatterer, atmosphere=atmosphere)
        at3d.util.to_2parameter_lwc_file('test_2paramlwc_filltemp.lwc', scatterer)
    def test_load_2parameter_atmos(self):
        self.scatterer_atmos = at3d.util.load_2parameter_lwc_file('test_2paramlwc_atmos.lwc')
        os.remove('test_2paramlwc_atmos.lwc')
    def test_load_2parameter_filltemp(self):
        self.scatterer = at3d.util.load_2parameter_lwc_file('test_2paramlwc_filltemp.lwc')
        os.remove('test_2paramlwc_filltemp.lwc')

class Microphysics_2parameter_lwc(TestCase):

    def setUp(self):
        scatterer, atmosphere = make_test_cloud()
        self.scatterer_original = scatterer
        self.atmosphere = atmosphere

        self.scatterer_original['temperature'] = self.atmosphere.interp({'z':self.scatterer_original.z}).temperature

        at3d.util.to_2parameter_lwc_file('test_2paramlwc_atmos.lwc', scatterer, atmosphere=atmosphere)
        at3d.util.to_2parameter_lwc_file('test_2paramlwc_filltemp.lwc', scatterer)

        self.scatterer_atmos = at3d.util.load_2parameter_lwc_file('test_2paramlwc_atmos.lwc')
        self.scatterer = at3d.util.load_2parameter_lwc_file('test_2paramlwc_filltemp.lwc')

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
        combine = at3d.grid.merge_two_z_coordinates(z1,z2)
        self.z1 = z1
        self.z2 = z2
        self.combine = combine
    def test_equivalent(self):
        self.assertTrue(np.all(at3d.grid.merge_two_z_coordinates(self.z1,self.z2)==at3d.grid.merge_two_z_coordinates(self.z2,self.z1)))
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
        self.combine = at3d.grid.combine_z_coordinates(self.zlist)
    def test_combine(self):
         combined = at3d.grid.combine_z_coordinates(self.zlist)
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
        self.rte_grid = at3d.grid.make_grid(2.5/29,30, 1.6/27 ,28,np.linspace(0.0,4.6,19))
    def test_resample_3d_onto_grid(self):
        at3d.grid.resample_onto_grid(self.rte_grid, self.cloud)
    def test_resample_1d_onto_grid(self):
        at3d.grid.resample_onto_grid(self.rte_grid, self.atmosphere)

class Resample_onto_grid_values(TestCase):
    def setUp(self):
        cloud, atmosphere = make_test_cloud()
        self.cloud = cloud
        self.atmosphere = atmosphere
        self.rte_grid = at3d.grid.make_grid(2.5/29,30, 1.6/27 ,28,np.linspace(0.0,4.6,19))
        self.grid_1d = at3d.grid.resample_onto_grid(self.rte_grid, atmosphere)
        self.grid_3d = at3d.grid.resample_onto_grid(self.rte_grid, cloud)

    def test_resample_3d_onto_grid_valid_density(self):
        self.assertTrue(np.bitwise_not(np.all(np.isnan(self.grid_3d.density.data))))
    def test_resample_3d_onto_grid_valid_reff(self):
        self.assertTrue(np.bitwise_not(np.all(np.isnan(self.grid_3d.reff.data))))
    def test_resample_1d_onto_grid_valid_temp(self):
        self.assertTrue(np.bitwise_not(np.all(np.isnan(self.grid_1d.temperature.data))))

class Adjoint_Linear_Interpolation(TestCase):
    def dotproduct_test(self):

        np.random.seed(1)
        xinput = np.random.uniform(size=(2,1,1))
        yinput = np.random.uniform(size=(30,1,1))

        test = xr.Dataset(data_vars={
            'h': (['x','y','z'],xinput)
        },
                         coords={'x': np.linspace(0.0,1.0,2),
                                'y': np.array([0.0]),
                                'z': np.array([0.0])})
        test2 = test.interp(x = np.linspace(0.0,1.0,30))
        ytilde = test2.h.data

        xtilde = at3d.core.adjoint_linear_interpolation(
            xgrid=test.x.data,
            ygrid=test.y.data,
            zgrid=test.z.data,
            inxgrid=test2.x.data,
            inygrid=test2.y.data,
            inzgrid=test2.z.data,
            field=yinput
                                                 )

        dot1 = np.dot(xtilde.ravel(),xinput.ravel())
        dot2 = np.dot(ytilde.ravel(),yinput.ravel())
        self.assertTrue(np.allclose(dot1,dot2))
