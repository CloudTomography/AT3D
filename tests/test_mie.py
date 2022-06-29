from unittest import TestCase
import os
import numpy as np
import xarray as xr
import pathlib

import at3d.mie

class Mie_tables(TestCase):
    @classmethod
    def setUpClass(cls, path='data/test_load_mie_table.nc'):
        cls.mie_water_table = at3d.mie.get_mono_table('Water', (0.6, 0.6),
                                                     minimum_effective_radius=4.0,
                                                     max_integration_radius=45.0,
                                                     wavelength_averaging=False,
                                                     verbose=False)

        cls.mie_aerosol_table = at3d.mie.get_mono_table('Aerosol', (0.75, 0.8),
                                                           minimum_effective_radius=5.0,
                                                           max_integration_radius=20.0,
                                                           wavelength_averaging=True,
                                                           wavelength_resolution=0.01,
                                                           refractive_index=1.53-0.08j,
                                                           verbose=False)
        cls.mie_water_table.to_netcdf(path)
        cls.path = path
        cls.relative_dir = str(pathlib.Path(path).parent)
        size_distribution = at3d.size_distribution.get_size_distribution_grid(
            cls.mie_water_table.radius,
            size_distribution_function=at3d.size_distribution.gamma,
            particle_density=1.0,
            reff={'coord_min':4.0, 'coord_max': 25.0, 'npoints': 25,
            'spacing': 'logarithmic', 'units': 'micron'},
            veff={'coord_min':0.09, 'coord_max': 0.11, 'npoints': 2,
            'spacing': 'linear', 'units': 'unitless'}
        )
        cls.poly_table = at3d.mie.get_poly_table(size_distribution, cls.mie_water_table)

    def test_load_table(self):
        loaded_mie_table = at3d.mie._load_table(relative_dir=self.relative_dir,
                                                   particle_type='Water',
                                                   wavelength_band=(0.6, 0.6),
                                                   minimum_effective_radius=4.0,
                                                   max_integration_radius=45.0,
                                                   wavelength_averaging=False,
                                                   wavelength_resolution=0.001,
                                                   refractive_index=None)
        self.assertIsNotNone(loaded_mie_table)
        self.assertTrue(np.allclose(self.mie_water_table.extinction.data, loaded_mie_table.extinction.data))
        self.assertTrue(np.allclose(self.mie_water_table.scatter.data, loaded_mie_table.scatter.data))
        self.assertTrue(np.allclose(self.mie_water_table.nleg.data, loaded_mie_table.nleg.data))
        self.assertTrue(np.allclose(self.mie_water_table.legendre.data, loaded_mie_table.legendre.data))

    def test_water_table_extinction(self, path='data/water600nm_extinction.nc'):
        self.assertTrue(np.allclose(self.mie_water_table.extinction.data, xr.load_dataarray(path).data))
    def test_water_table_scatter(self, path='data/water600nm_scatter.nc'):
        self.assertTrue(np.allclose(self.mie_water_table.scatter.data, xr.load_dataarray(path).data))
    def test_water_table_nleg(self, path='data/water600nm_nleg.nc'):
        self.assertTrue(np.allclose(self.mie_water_table.nleg.data, xr.load_dataarray(path).data))
    def test_water_table_legnedre_mean(self, path='data/water600nm_legendre_mean.nc'):
        self.assertTrue(np.allclose(self.mie_water_table.legendre.mean('legendre_index').data,
                                    xr.load_dataarray(path).data))

    def test_aerosol_table_extinction(self, path='data/aerosol750-800nm_extinction.nc'):
        self.assertTrue(np.allclose(self.mie_aerosol_table.extinction.data, xr.load_dataarray(path).data))
    def test_aerosol_table_scatter(self, path='data/aerosol750-800nm_scatter.nc'):
        self.assertTrue(np.allclose(self.mie_aerosol_table.scatter.data, xr.load_dataarray(path).data))
    def test_aerosol_table_nleg(self, path='data/aerosol750-800nm_nleg.nc'):
        self.assertTrue(np.allclose(self.mie_aerosol_table.nleg.data, xr.load_dataarray(path).data))
    def test_aerosol_table_legnedre_mean(self, path='data/aerosol750-800nm_legendre_mean.nc'):
        self.assertTrue(np.allclose(self.mie_aerosol_table.legendre.mean('legendre_index').data,
                                    xr.load_dataarray(path).data))

    def test_poly_table_extinction(self, path='data/poly_table_test.nc'):
        with xr.open_dataset(path) as dataset:
            ref_extinction = dataset['extinction']
        self.assertTrue(np.allclose(self.poly_table.extinction.data, ref_extinction.data))
    def test_poly_table_ssalb(self, path='data/poly_table_test.nc'):
        with xr.open_dataset(path) as dataset:
            ref_ssalb = dataset['ssalb']
        self.assertTrue(np.allclose(self.poly_table.ssalb.data, ref_ssalb.data))
    def test_poly_table_legendre(self, path='data/poly_table_test.nc'):
        with xr.open_dataset(path) as dataset:
            ref_legcoef = dataset['legcoef_mean']
        self.assertTrue(np.allclose(self.poly_table.legcoef.mean('legendre_index').data, ref_legcoef.data))

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.path)
