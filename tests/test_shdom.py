from unittest import TestCase
import numpy as np
import os
import shdom
import pandas as pd
import xarray as xr
from collections import OrderedDict

def calculate_mie(mie_dir, particle_type, wavelength):
    mie_fname = '{}_{:.3f}.nc'.format(particle_type.lower(), wavelength)
    mie_path = os.path.join(mie_dir, mie_fname)
    if not os.path.exists(mie_path):
        mie_table = shdom.mie.get_mono_table(particle_type, (wavelength,wavelength), relative_path=mie_path)
    else:
        mie_table = xr.load_dataset(mie_path)
    number_density_grid = shdom.size_distribution.get_size_distribution_grid(
                                                        mie_table['radius'],
                                                        size_distribution_function=shdom.size_distribution.gamma,
                                                        reff=[4.5, 20.0, 10, 'linear', 'micron'],
                                                        veff=[0.01, 0.1, 5, 'linear', 'unitless'],
                                                        radius_units=mie_table.attrs['units'][0]
                                                        )
    poly_table = shdom.mie.get_poly_table(number_density_grid, mie_table)
    return poly_table

class Load_from_csv_test(TestCase):
    def test_load_from_csv(self, path='../synthetic_cloud_fields/jpl_les/rico32x37x26.txt'):
        scatterer = shdom.grid.load_from_csv(path, density='lwc', origin=(0.0,0.0))

class Microphysics_load_from_csv_test(TestCase):
    def setUp(self):
        scatterer = shdom.grid.load_from_csv('../synthetic_cloud_fields/jpl_les/rico32x37x26.txt', density='lwc', origin=(0.0,0.0))
        self.droplets = shdom.grid.resample_onto_grid(scatterer, scatterer)
        self.droplets['veff'] = (self.droplets.reff.dims, np.full_like(self.droplets.reff.data, fill_value=0.1))

    def test_load_from_csv_lwc(self):
            self.assertAlmostEqual(self.droplets.density.data[self.droplets.density.data>0].mean(), 0.265432, places=6)

    def test_load_from_csv_reff(self):
        self.assertAlmostEqual(self.droplets.reff.data[self.droplets.reff.data>0].mean(), 15.713033, places=6)

    def test_load_from_csv_veff(self):
        self.assertAlmostEqual(self.droplets.veff.data[self.droplets.veff.data>0].mean(), 0.100000, places=6)


class Microphysics_save_to_csv_test(TestCase):
    def setUp(self):
        scatterer = shdom.grid.load_from_csv('../synthetic_cloud_fields/jpl_les/rico32x37x26.txt', density='lwc', origin=(0.0,0.0))
        self.test_fname = 'test_rico32x37x26.nc'
        self.droplets = shdom.grid.resample_onto_grid(scatterer, scatterer)
        self.droplets['veff'] = (self.droplets.reff.dims, np.full_like(self.droplets.reff.data, fill_value=0.1))
        self.droplets.to_netcdf(test_fname)

    def test_saved_microphysics(self):
        new_droplets = xr.load_dataset(self.test_fname)
        self.assertEqual(np.sum(self.droplets.density.data - new_droplets.density.data), 0)
        self.assertEqual(np.sum(self.droplets.reff.data - new_droplets.reff.data), 0)
        self.assertEqual(np.sum(self.droplets.veff.data - new_droplets.veff.data), 0)

    def tearDown(self):
        os.remove(self.test_fname)

class Shdom_rt_simulation_test(TestCase):
    def setUp(self):
        self.wavelength = 0.672
        self.droplets = shdom.grid.load_from_csv('../synthetic_cloud_fields/jpl_les/rico32x37x26.txt', density='lwc', origin=(0.0,0.0))
        scatterer_grid = shdom.grid.resample_onto_grid(self.droplets, self.droplets)
        scatterer_grid['veff'] = (scatterer_grid.reff.dims, np.full_like(scatterer_grid.reff.data, fill_value=0.1))
        poly_table = calculate_mie('../mie_tables', 'Water', self.wavelength)
        optical_properties = shdom.medium.table_to_grid(scatterer_grid, poly_table)

        config_dataset = shdom.configuration.get_config('../orig_config.json')
        surface = shdom.surface.lambertian(albedo=0.0)
        source = shdom.source.solar(solar_zenith=0, solar_azimuth=0.0)
        self.solvers = OrderedDict()
        self.solvers[self.wavelength] = shdom.solver.RTE(
                                        config_dataset,
                                        optical_properties,
                                        source,
                                        surface,
                                        num_stokes=3,
                                        name=None
                                        )

    def test_forward_solver(self):
        sensors = OrderedDict()
        camera = shdom.sensor.orthographic_projection(
                        wavelength=self.wavelength,
                        bounding_box=self.droplets,
                        x_resolution=0.02,
                        y_resolution=0.02,
                        azimuth=0,
                        zenith=0,
                        altitude="TOA",
                        stokes=["I", "Q", "U"]
                        )
        camera_rays = shdom.sensor.add_sub_pixel_rays(camera, 0, degree=2)
        sensors['MISR'] = {'sensor_list': [camera_rays]}
        shdom.script_util.get_measurements(self.solvers, sensors, maxiter=4, n_jobs=40, verbose=True)

        sensor = sensors['MISR']['sensor_list'][0]
        I = sensor.I.data.reshape(sensor.image_shape.data, order='F')
        Q = sensor.Q.data.reshape(sensor.image_shape.data, order='F')
        U = sensor.U.data.reshape(sensor.image_shape.data, order='F')
        self.assertAlmostEqual(I[I>0].mean(), 0.08, places=2)
        self.assertAlmostEqual(I[I>0].std(), 0.05, places=2)
        self.assertAlmostEqual(Q[Q>0].mean(), 0.000131, places=6)
        self.assertAlmostEqual(Q[Q>0].std(), 0.000111, places=6)
        self.assertAlmostEqual(U[U>0].mean(), 0.000136, places=6)
        self.assertAlmostEqual(U[U>0].std(), 0.000135, places=6)

    def test_forward_solver_oblique(self):
        sensors = OrderedDict()
        camera = shdom.sensor.orthographic_projection(
                        wavelength=self.wavelength,
                        bounding_box=self.droplets,
                        x_resolution=0.02,
                        y_resolution=0.02,
                        azimuth=0,
                        zenith=26.1,
                        altitude="TOA",
                        stokes=["I", "Q", "U"]
                        )
        camera_rays = shdom.sensor.add_sub_pixel_rays(camera, 0, degree=2)
        sensors['MISR'] = {'sensor_list': [camera_rays]}
        shdom.script_util.get_measurements(self.solvers, sensors, maxiter=4, n_jobs=40, verbose=True)

        sensor = sensors['MISR']['sensor_list'][0]
        I = sensor.I.data.reshape(sensor.image_shape.data, order='F')
        Q = sensor.Q.data.reshape(sensor.image_shape.data, order='F')
        U = sensor.U.data.reshape(sensor.image_shape.data, order='F')
        self.assertAlmostEqual(I[I>0].mean(), 0.08, places=2)
        self.assertAlmostEqual(I[I>0].std(), 0.05, places=2)
        self.assertAlmostEqual(Q[Q>0].mean(), 0.000131, places=6)
        self.assertAlmostEqual(Q[Q>0].std(), 0.000111, places=6)
        self.assertAlmostEqual(U[U>0].mean(), 0.000136, places=6)
        self.assertAlmostEqual(U[U>0].std(), 0.000135, places=6)
