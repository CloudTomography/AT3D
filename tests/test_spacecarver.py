import xarray as xr
import numpy as np
import at3d

from unittest import TestCase

class VerifySpaceCarver(TestCase):
    @classmethod
    def setUpClass(cls):
        np.random.seed(1)

        grid = at3d.grid.make_grid(1.0, 10, 1.2, 12, z=np.linspace(0,1.2,13))

        space_carver = at3d.space_carve.SpaceCarver(grid, 3)

        density = np.zeros((grid.x.size, grid.y.size, grid.z.size))
        xinput = np.random.uniform(0.1, 100.0, size=density.shape)
        density[:,:,:] = xinput
        weights = xr.Dataset(
            data_vars={
                'density': (['x', 'y', 'z'], density)
            },
            coords=grid.coords
        )

        sensors = at3d.containers.SensorsDict()

        zeniths = [88.0, 0.0]
        azimuths = [45.0,45.0]
        for zenith, azimuth in zip(zeniths, azimuths):
            sensors.add_sensor('TEST',
                              at3d.sensor.orthographic_projection(0.672,
                                                                    grid,
                                                                    0.5,0.5,azimuth,zenith,
                                                                    stokes=['I'])
                              )
        space_carver.project(weights, sensors)

        sensors2 = sensors.make_forward_sensors()
        yinput1 = np.random.uniform(0.1, 20.0,size=sensors2['TEST']['sensor_list'][0].nrays.size)
        yinput2 = np.random.uniform(0.1, 20.0,size=sensors2['TEST']['sensor_list'][1].nrays.size)


        sensors2['TEST']['sensor_list'][0]['weights'] = (['nrays'], yinput1)
        sensors2['TEST']['sensor_list'][0]['cloud_mask'] = (['nrays'], np.ones(yinput1.size))
        sensors2['TEST']['sensor_list'][1]['weights'] = (['nrays'], yinput2)
        sensors2['TEST']['sensor_list'][1]['cloud_mask'] = (['nrays'], np.ones(yinput2.size))

        space_carver2 = at3d.space_carve.SpaceCarver(grid, 3)

        output = space_carver2.carve(sensors2, linear_mode=True)
        xtilde = output.weights.data

        ytilde1 = sensors['TEST']['sensor_list'][0].weights.data
        ytilde2 = sensors['TEST']['sensor_list'][1].weights.data
        cls.oneside = np.dot(xtilde.sum(axis=0).ravel(),xinput.ravel())
        cls.otherside = np.dot(ytilde1,yinput1) + np.dot(ytilde2,yinput2)
    def dotproduct_test1(self):
        self.assertEqual(83123068.05097, np.round(self.oneside,5))
    def dotproduct_test2(self):
        self.assertEqual(83123068.05097, np.round(self.oneside,5))
