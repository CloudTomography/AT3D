"""
A quick code for 1D bispectral retrievals.
This is not validated, use at your own risk.
In particular, it uses a black box local minimization method
which does not have the same retrieval logic as e.g. MODIS operational
retrieval in the presence of multiple solutions at small
optical depths.

TODO - add more angles, varying surface,
     - add documentation.
"""

import scipy.interpolate as si
import scipy.optimize as so
import numpy as np
from collections import OrderedDict

import pyshdom

class BispectralRetrieval:
    """
    Single view/solar/angle bispectral retrieval using a look up table.
    """
    def __init__(self, nir_wavelength=0.865, swir_wavelength=2.265, view_zenith=0.0,
                 view_azimuth=0.0,
                 solar_mu=0.5, solar_azimuth=0.0, tau_min=0.1, tau_max=100.0, ntau=42,
                 reff_min=5.0, reff_max=35.0, nreff=32, veff=0.1,
                 mie_save_directory='../mie_tables', ext0=4.0, reff0=10.0,
                 solver_config=None):

        self.ext = np.logspace(np.log10(tau_min), np.log10(tau_max), ntau)
        self.reff = np.linspace(reff_min, reff_max, nreff)
        self.veff = veff

        self.nir_wavelength = nir_wavelength
        self.swir_wavelength = swir_wavelength
        self.mie_save_directory = mie_save_directory

        self.view_zenith = view_zenith
        self.view_azimuth = view_azimuth
        self.solar_mu = solar_mu
        self.solar_azimuth = solar_azimuth

        self.ext0 = ext0
        self.reff0 = reff0

        if solver_config is None:
            config = pyshdom.configuration.get_config('../default_config.json')
            config['num_mu_bins'] = 8
            config['num_phi_bins'] = 16
            config['split_accuracy'] = 0.03
            config['spherical_harmonics_accuracy'] = 0.0
            config['adapt_grid_factor'] = 100.0
            config['solution_accuracy'] = 1e-6
            self.solver_config = config
        else:
            self.solver_config = solver_config

        self._prepare_mie_tables()
        self._make_lut()


    def _prepare_mie_tables(self):

        self.mie_tables = OrderedDict()
        self.mie_tables[0.865] = pyshdom.mie.get_mono_table(
            'Water',
            (self.nir_wavelength, self.nir_wavelength),
            max_integration_radius=65.0,
            minimum_effective_radius=0.1,
            relative_dir=self.mie_save_directory,
            verbose=True
            )
        self.mie_tables[2.265] = pyshdom.mie.get_mono_table(
            'Water',
            (self.swir_wavelength, self.swir_wavelength),
            max_integration_radius=65.0,
            minimum_effective_radius=0.1,
            relative_dir=self.mie_save_directory,
            verbose=True
            )

        self.cloud_poly_tables = OrderedDict()
        size_distribution_function = pyshdom.size_distribution.gamma
        for wavelength in self.mie_tables.keys():
            mie_mono_table = self.mie_tables[wavelength]
            cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
                                mie_mono_table.radius.data,
                                size_distribution_function=size_distribution_function, particle_density=1.0,
                                reff={'coord_min':self.reff.min()-0.1, 'coord_max': self.reff.max()+0.1, 'npoints': 100,
                                'spacing': 'logarithmic', 'units': 'micron'},
                                veff={'coord_min':0.09, 'coord_max': 0.11, 'npoints': 2,
                                'spacing': 'linear', 'units': 'unitless'}
                                )
            poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution, mie_mono_table)
            self.cloud_poly_tables[wavelength] = poly_table

        scale_table = self.cloud_poly_tables[self.nir_wavelength].copy(deep=True)
        for wavelength, table in self.cloud_poly_tables.items():
            table['extinction'][:] /= scale_table.extinction.data

    def forward_model(self, taus, reffs):
        # put taus in increasing order for maximum efficiency
        taus = np.atleast_1d(taus)
        reffs = np.atleast_1d(reffs)

        output_nir = np.zeros(taus.size)-1
        output_swir = np.zeros(taus.size)-1
        nir_solution = None
        swir_solution = None
        for i, (tau, reff) in enumerate(zip(taus.ravel(), reffs.ravel())):
            ext = np.array([tau, tau])
            reffgrid = np.array([reff, reff])
            sensor_dict, solver, nir_solution, swir_solution = self.rte_model(
                ext, reffgrid, nir_solution=nir_solution, swir_solution=swir_solution
                )
            nir = sensor_dict['NIR']['sensor_list'][0].I.data.max()
            swir = sensor_dict['SWIR']['sensor_list'][0].I.data.max()
            output_nir[i] = nir
            output_swir[i] = swir
        output_nirs = output_nir.reshape(taus.shape)
        output_swirs = output_swir.reshape(reffs.shape)
        return output_nirs, output_swirs

    def rte_model(self, ext, reff, nir_solution=None, swir_solution=None):

        rte_grid = pyshdom.grid.make_grid(0.05, len(ext), 0.05, len(reff), np.linspace(0.1, 1.1,4))
        grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
        extinction, reff_grid = np.meshgrid(ext, reff, indexing='ij')
        rte_grid['density'] = (['x', 'y', 'z'],
            np.repeat(extinction[..., np.newaxis], grid_shape[-1], axis=-1))
        rte_grid['reff'] = (['x', 'y', 'z'], np.repeat(reff_grid[...,
            np.newaxis], grid_shape[-1], axis=-1))
        rte_grid['veff'] = (['x', 'y', 'z'] , np.zeros(grid_shape) + self.veff)
        cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, rte_grid)

        x = rte_grid.x.data
        y = rte_grid.y.data
        x2, y2 = np.meshgrid(x,y, indexing='ij')
        xs = x2.ravel()
        ys = y2.ravel()
        zs = np.array([rte_grid.z.data[-1]]*len(xs))
        mus = np.array([np.cos(np.deg2rad(self.view_zenith))]*len(xs))
        phis = np.array([self.view_azimuth]*len(xs))
        Sensordict = pyshdom.containers.SensorsDict()
        Sensordict.add_sensor('NIR',
                pyshdom.sensor.make_sensor_dataset(x=xs, y=ys, z=zs, mu=mus, phi=phis,
                                                  wavelength=self.nir_wavelength, stokes=['I'],
                                                  fill_ray_variables=True)
                             )
        Sensordict.add_sensor('SWIR',
                pyshdom.sensor.make_sensor_dataset(x=xs, y=ys, z=zs, mu=mus, phi=phis,
                                                  wavelength=self.swir_wavelength, stokes=['I'],
                                                  fill_ray_variables=True)
                             )

        wavelengths = Sensordict.get_unique_solvers()

        solvers_dict = pyshdom.containers.SolversDict()
        for wavelength in wavelengths:
            poly_table = self.cloud_poly_tables[wavelength]
            optical_properties = pyshdom.medium.table_to_grid(
                cloud_scatterer_on_rte_grid, poly_table, exact_table=False
                )
            config = self.solver_config
            config['ip_flag'] = 3 # make sure 1D mode.

            solver = pyshdom.solver.RTE(
                numerical_params=config,
                medium={'cloud': optical_properties},
                source=pyshdom.source.solar(wavelength, self.solar_mu, self.solar_azimuth),
                surface=pyshdom.surface.lambertian(albedo=0.03),
                num_stokes=1,
                name=None
            )
            solvers_dict.add_solver(wavelength, solver)
        if nir_solution is not None:
            solvers_dict[wavelengths[0]].load_solution(nir_solution)
        if swir_solution is not None:
            solvers_dict[wavelengths[1]].load_solution(swir_solution)

        Sensordict.get_measurements(solvers_dict, maxiter=800)

        nir_saved = solvers_dict[wavelengths[0]].save_solution()
        swir_saved = solvers_dict[wavelengths[1]].save_solution()

        return Sensordict, solvers_dict, nir_saved, swir_saved

    def _make_lut(self):
        print('Making Look-Up-Table, this will take a little while. . .')
        sensor_dict, solvers, nir_solution, swir_solution = self.rte_model(self.ext, self.reff)
        self.nirs_lut = sensor_dict['NIR']['sensor_list'][0].I.data.reshape(
            (self.ext.size, self.reff.size)
        )
        self.swirs_lut = sensor_dict['SWIR']['sensor_list'][0].I.data.reshape(
            (self.ext.size, self.reff.size)
        )
        self.nir_solution = nir_solution
        self.swir_solution = swir_solution

    def retrieve(self, nirs_to_test, swirs_to_test):
        nirs_to_test = np.atleast_1d(nirs_to_test)
        swirs_to_test = np.atleast_1d(swirs_to_test)
        retrieved_reff = np.zeros(nirs_to_test.ravel().shape)
        retrieved_tau = np.zeros(retrieved_reff.shape)
        costs = np.zeros(retrieved_reff.shape)
        for i, (nir_to_test, swir_to_test) in enumerate(
                zip(nirs_to_test.ravel(), swirs_to_test.ravel())
            ):
            cost = 0.5*((np.log(nir_to_test) - np.log(self.nirs_lut))**2 + \
                        (np.log(swir_to_test) - np.log(self.swirs_lut))**2)

            cost_spline = si.RectBivariateSpline(self.ext, self.reff, cost,
                                                 bbox=[self.ext.min(), self.ext.max(),
                                                       self.reff.min(), self.reff.max()],
                                                 kx=3, ky=3, s=0)
            def cost_function(state):

                cost = cost_spline(state[0], state[1])
                gradx = cost_spline(state[0], state[1], dx=1)
                grady = cost_spline(state[0], state[1], dy=1)
                return cost, np.array([gradx[0, 0], grady[0, 0]])

            result = so.minimize(
                cost_function,
                x0=np.array([self.ext0, self.reff0]),
                jac=True,
                method='L-BFGS-B',
                bounds=[(self.ext.min(), self.ext.max()),
                        (self.reff.min(), self.reff.max())]
            )
            retrieved_reff[i] = result.x[1]
            retrieved_tau[i] = result.x[0]
            costs[i] = result.fun

        return retrieved_tau, retrieved_reff, costs
