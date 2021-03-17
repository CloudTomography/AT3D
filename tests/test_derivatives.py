from unittest import TestCase
from collections import OrderedDict
import numpy as np
import xarray as xr
import pyshdom

import warnings
warnings.filterwarnings('ignore')

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


# class Verify_Jacobian(TestCase):
#     @classmethod
#     def setUpClass(cls):
#
#         ext = 0.1
#         veff = 0.1
#         reff=10.0
#         rte_grid = pyshdom.grid.make_grid(0.05, 3, 0.05, 3, np.arange(0.1, 0.25, 0.05))
#         grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
#         rte_grid['density'] = (['x','y','z'], np.ones(grid_shape))
#         rte_grid['reff'] = (['x','y','z'], np.zeros(grid_shape) + reff)
#         rte_grid['veff'] = (['x','y','z'] ,np.zeros(grid_shape) + veff)
#
#         #resample the cloud onto the rte_grid
#         cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, rte_grid)
#
#         #define any necessary variables for microphysics here.
#         size_distribution_function = pyshdom.size_distribution.gamma
#
#         #define sensors.
#         Sensordict = pyshdom.containers.SensorsDict()
#
#         sensor = pyshdom.sensor.make_sensor_dataset(np.array([0.05]),np.array([0.05]),
#                                                   np.array([0.7]),np.array([1.0]),
#                                                   np.array([0.0]),
#                                                  wavelength=0.86,stokes=['I'],fill_ray_variables=True)
#         Sensordict.add_sensor('MISR', sensor)
#
#         #Define the RTE solvers needed to model the measurements and
#         #calculate optical properties.
#         wavelengths = Sensordict.get_unique_solvers()
#
#         cloud_poly_tables = OrderedDict()
#         solvers = pyshdom.containers.SolversDict()
#
#         for wavelength in wavelengths:
#
#             #optical properties from mie calculations.
#             mie_mono_table = pyshdom.mie.get_mono_table('Water',(wavelength,wavelength),
#                                                       max_integration_radius=65.0,
#                                                       minimum_effective_radius=0.1,
#                                                       relative_dir='../mie_tables',
#                                                       verbose=False)
#             cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
#                                                                     mie_mono_table.radius.data,
#                                 size_distribution_function=size_distribution_function,particle_density=1.0,
#                                 reff={'coord_min':9.0, 'coord_max': 11.0, 'npoints': 100,
#                                 'spacing': 'logarithmic', 'units': 'micron'},
#                                 veff={'coord_min':0.09, 'coord_max': 0.11, 'npoints': 12,
#                                 'spacing': 'linear', 'units': 'unitless'}
#                                 )
#             poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
#             optical_properties = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
#                                                                              exact_table=False)
#             optical_properties['ssalb'][:,:,:] = 1.0
#             extinction = np.zeros(optical_properties.extinction.shape)
#             extinction[1:-1,1:-1,1:-1] = ext
#             #extinction[a,b,c] += step
#             optical_properties['legcoef'][:,1:,:] = 0.0
#             optical_properties['extinction'][:,:,:] = extinction
#             cloud_poly_tables[wavelength] = poly_table
#             config = pyshdom.configuration.get_config('../default_config.json')
#             config['num_mu_bins'] = 16
#             config['num_phi_bins'] = 32
#             config['split_accuracy'] = 0.0003
#             config['spherical_harmonics_accuracy'] = 0.0
#             config['solution_accuracy'] = 1e-5
#             config['deltam'] = True
#             solver = pyshdom.solver.RTE(
#                                 numerical_params=config,
#                                 medium={'cloud': optical_properties},
#                                 source=pyshdom.source.solar(wavelength,-1,0.0,solarflux=1.0),
#                                 surface=pyshdom.surface.lambertian(albedo=0.0),
#                                 num_stokes=1,
#                                 name=None
#                                 )
#
#             solvers.add_solver(wavelength, solver)
#         Sensordict.get_measurements(solvers, maxiter=100, n_jobs=8, verbose=False)
#
#         unknown_scatterers = pyshdom.containers.UnknownScatterers()
#         unknown_scatterers.add_unknown('cloud', ['extinction'], cloud_poly_tables)
#         unknown_scatterers.create_derivative_tables()
#         solvers.add_microphysical_partial_derivatives(unknown_scatterers)
#
#         forward_sensors = Sensordict.make_forward_sensors()
#
#         gradient_call = pyshdom.gradient.LevisApproxGradientUncorrelated(Sensordict,
#         solvers, forward_sensors, unknown_scatterers,
#         parallel_solve_kwargs={'n_jobs':4, 'maxiter': 100, 'setup_grid':True, 'verbose':False},
#         gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
#         'indices_for_jacobian': ([1],[1],[1])}, uncertainty_kwargs={'add_noise': False})
#         out, gradient, jacobian_exact = gradient_call()
#
#         cls.jacobian_exact = jacobian_exact
#
#     def test_jacobian(self):
#         self.assertAlmostEqual(self.jacobian_exact['jacobian_0.860'][0,0,0,0].data, 0.00396336, places=5)

#A function for computing clouds for the thermal jacobian reference case.
def cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=0.0, index=(1,1,1),
          nmu=16, split=0.03, load_solution=None, resolution=1):

    rte_grid = pyshdom.grid.make_grid(0.05/resolution, 7*resolution, 0.05/resolution, 7*resolution,
                                      np.arange(0.1, 0.5, 0.05/resolution))

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
    zeniths =  [75.0,60.0,45.6,26.1]*2 + [0.0] +  [75.0,60.0,45.6,26.1]*6
    azimuths =  [90.0]*4 + [-90]*4 + [0.0] + [0]*4 + [180]*4 + [-45]*4 + [135]*4+ [45]*4 + [-135]*4

    target = [rte_grid.x.data[-1]/2, rte_grid.y.data[-1]/2, rte_grid.z.data[rte_grid.z.size//2]]
    R = 10.0
    xs = target[0] + R*np.sin(np.deg2rad(zeniths))*np.cos(np.deg2rad(azimuths))
    ys = target[1] + R*np.sin(np.deg2rad(zeniths))*np.sin(np.deg2rad(azimuths))
    zs = target[2] + R*np.cos(np.deg2rad(zeniths))
    for x,y,z in zip(xs,ys,zs):
        Sensordict.add_sensor('MISR',
        pyshdom.sensor.perspective_projection(11.0, 5.0, 26, 26,[x,y,z], target,
                                              [0,1,0],stokes=['I'])
                             )

    wavelengths = Sensordict.get_unique_solvers()

    cloud_poly_tables = OrderedDict()
    solvers = pyshdom.containers.SolversDict()
    temp = np.linspace(280,300, grid_shape[-1])[np.newaxis,np.newaxis,:]
    atmosphere = xr.Dataset(data_vars={
        'temperature': (['x','y','z'], np.repeat(np.repeat(temp,grid_shape[0], axis=0),grid_shape[1],axis=1))
    }, coords={'x':rte_grid.x, 'y': rte_grid.y, 'z': rte_grid.z})
    atmosphere2 = pyshdom.grid.resample_onto_grid(rte_grid, atmosphere)

    for wavelength in wavelengths:

        cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
                                                                mie_mono_table.radius.data,
                            size_distribution_function=size_distribution_function,particle_density=1.0,
                            reff={'coord_min':0.1, 'coord_max': 11.0, 'npoints': 100,
                            'spacing': 'logarithmic', 'units': 'micron'},
                            veff={'coord_min':0.09, 'coord_max': 0.11, 'npoints': 12,
                            'spacing': 'linear', 'units': 'unitless'}
                            )
        poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
        optical_properties = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                         exact_table=False)
        optical_properties['ssalb'][:,:,:] = ssalb
        extinction = np.zeros(optical_properties.extinction.shape)
        x,y,z = np.meshgrid(np.arange(extinction.shape[0]),
                       np.arange(extinction.shape[1]),
                       np.arange(extinction.shape[2]),indexing='ij')
        R = np.sqrt((x-x.mean())**2+ (y-y.mean())**2 + (z-z.mean())**2)
        extinction[1:-1,1:-1,1:-1] = ext*np.exp(-1*R/(0.2*R.mean()))[1:-1,1:-1,1:-1]
        extinction[index[0],index[1],index[2]] += step
        optical_properties['extinction'][:,:,:] = extinction

        cloud_poly_tables[wavelength] = poly_table
        config = pyshdom.configuration.get_config('../default_config.json')
        config['num_mu_bins'] = nmu
        config['num_phi_bins'] = nmu*2
        config['split_accuracy'] = split
        config['spherical_harmonics_accuracy'] = 0.0
        config['adapt_grid_factor'] = 1000.0
        config['solution_accuracy'] = 1e-7
        config['high_order_radiance'] = False
        config['acceleration_flag'] = True
        config['deltam'] = False
        solver = pyshdom.solver.RTE(
                            numerical_params=config,
                            medium={'cloud': optical_properties},
                            source=pyshdom.source.thermal(wavelength),
                            surface=pyshdom.surface.lambertian(albedo=surfacealb,
                            ground_temperature=ground_temperature),
                            num_stokes=1,
                            atmosphere=atmosphere2,
                            name=None
                            )
        if load_solution is not None:
            solver.load_solution(load_solution)
        solvers.add_solver(wavelength,solver)

    Sensordict.get_measurements(solvers, maxiter=200, n_jobs=1, verbose=False)
    return solvers, Sensordict, cloud_poly_tables, step

class ThermalJacobianNoSurface(TestCase):
    @classmethod
    def setUpClass(cls):
        tautol=0.2
        surfacealb=0.0
        ssalb = 0.0
        solarmu = 1.0
        reff=10.0
        resolutionfactor = 1
        nmu=2
        split=0.5
        ext=10.0
        veff=0.1
        step=1e-4
        ground_temperature = 0.0
        mie_mono_table = pyshdom.mie.get_mono_table('Water',(11.0,11.0),
                                              max_integration_radius=65.0,
                                              minimum_effective_radius=0.1,
                                              relative_dir='../mie_tables',
                                              verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor)
        Sensordict.add_uncertainty_model('MISR', pyshdom.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[11.0]
        solver = solvers[11.0]

        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['extinction'], cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers.add_microphysical_partial_derivatives(unknown_scatterers)
        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = pyshdom.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':200,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': np.where(solver.medium['cloud'].extinction.data > 0.0)}, uncertainty_kwargs={'add_noise': False})
        cost, gradient, jacobian = gradient_call()
        cls.jacobian = jacobian['jacobian_11.000'][0,0].data

        # CODE THAT IS USED TO GENERATE THE FINITE DIFFERENCE REFERENCE
        # indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        # saved = solver.save_solution()
        # print("calculating finite differencing lineariztion. This will take a while. . .")
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #     data = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb, ground_temperature,step=step,index=(a,b,c),nmu=nmu,load_solution=saved, split=0.0,
        #                  resolution=resolutionfactor)
        #     data[1].add_uncertainty_model('MISR', pyshdom.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #
        #     data = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=-1*step,index=(a,b,c),nmu=nmu,load_solution=saved, split=0.0,
        #                  resolution=resolutionfactor)
        #     data[1].add_uncertainty_model('MISR', pyshdom.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_low, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #     out.append((rte_sensor_high[11.0].measurement_data[0].data - rte_sensor_low[11.0].measurement_data[0].data)/(2*step))#[a,b,c])
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/reference_noscat_nosurface_jacobian.npy', finite_jacobian)
        cls.jacobian_reference = np.load('./data/reference_noscat_nosurface_jacobian.npy')

    def test_jacobian(self):
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=8e-3))


class ThermalJacobianWithSurface(TestCase):
    @classmethod
    def setUpClass(cls):
        tautol=0.2
        surfacealb=0.0
        ssalb = 0.0
        solarmu = 1.0
        reff=10.0
        resolutionfactor = 1
        nmu=2
        split=0.5
        ext=10.0
        veff=0.1
        step=1e-4
        ground_temperature = 200.0
        mie_mono_table = pyshdom.mie.get_mono_table('Water',(11.0,11.0),
                                              max_integration_radius=65.0,
                                              minimum_effective_radius=0.1,
                                              relative_dir='../mie_tables',
                                              verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor)
        Sensordict.add_uncertainty_model('MISR', pyshdom.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[11.0]
        solver = solvers[11.0]

        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['extinction'], cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers.add_microphysical_partial_derivatives(unknown_scatterers)
        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = pyshdom.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':200,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': np.where(solver.medium['cloud'].extinction.data > 0.0)}, uncertainty_kwargs={'add_noise': False})
        cost, gradient, jacobian = gradient_call()
        cls.jacobian = jacobian['jacobian_11.000'][0,0].data

        # CODE THAT IS USED TO GENERATE THE FINITE DIFFERENCE REFERENCE
        # indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        # saved = solver.save_solution()
        # print("calculating finite differencing lineariztion. This will take a while. . .")
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #     data = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb, ground_temperature,step=step,index=(a,b,c),nmu=nmu,load_solution=saved, split=0.0,
        #                  resolution=resolutionfactor)
        #     data[1].add_uncertainty_model('MISR', pyshdom.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #
        #     data = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=-1*step,index=(a,b,c),nmu=nmu,load_solution=saved, split=0.0,
        #                  resolution=resolutionfactor)
        #     data[1].add_uncertainty_model('MISR', pyshdom.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_low, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #     out.append((rte_sensor_high[11.0].measurement_data[0].data - rte_sensor_low[11.0].measurement_data[0].data)/(2*step))#[a,b,c])
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/reference_noscat_surface_jacobian.npy', finite_jacobian)
        cls.jacobian_reference = np.load('./data/reference_noscat_surface_jacobian.npy')

    def test_jacobian(self):
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=8e-3))


def cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb, ground_temperature, step=0.0, index=(1,1,1),
          nmu=16, split=0.03, load_solution=None, resolution=1):
    rte_grid = pyshdom.grid.make_grid(0.05/resolution, 9*resolution, 0.05/resolution, 9*resolution,
                                      np.arange(0.1, 0.65, 0.05/resolution))

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
    zeniths =  [75.0,60.0,45.6,26.1]*2 + [0.0] +  [75.0,60.0,45.6,26.1]*6
    azimuths =  [90.0]*4 + [-90]*4 + [0.0] + [0]*4 + [180]*4 + [-45]*4 + [135]*4+ [45]*4 + [-135]*4

    target = [rte_grid.x.data[-1]/2, rte_grid.y.data[-1]/2, rte_grid.z.data[rte_grid.z.size//2]]
    R = 10.0
    xs = target[0] + R*np.sin(np.deg2rad(zeniths))*np.cos(np.deg2rad(azimuths))
    ys = target[1] + R*np.sin(np.deg2rad(zeniths))*np.sin(np.deg2rad(azimuths))
    zs = target[2] + R*np.cos(np.deg2rad(zeniths))
    for x,y,z in zip(xs,ys,zs):
        Sensordict.add_sensor('MISR',
        pyshdom.sensor.perspective_projection(0.86, 4.0, 13, 13,[x,y,z], target,
                                              [0,1,0],stokes=['I'])
                             )

    sensor = pyshdom.sensor.make_sensor_dataset(np.array([0.15]),np.array([0.055]),
                                              np.array([0.44999998807907104]),np.array([1.0]),
                                             np.array([0.0]),
                                             wavelength=0.86,stokes=['I'],fill_ray_variables=True)

    Sensordict.add_sensor('MISR', sensor)
    wavelengths = Sensordict.get_unique_solvers()

    cloud_poly_tables = OrderedDict()
    solvers = pyshdom.containers.SolversDict()
    temp = np.linspace(280,300, grid_shape[-1])[np.newaxis,np.newaxis,:]
    atmosphere = xr.Dataset(data_vars={
        'temperature': (['x','y','z'], np.repeat(np.repeat(temp,grid_shape[0], axis=0),grid_shape[1],axis=1))
    }, coords={'x':rte_grid.x, 'y': rte_grid.y, 'z': rte_grid.z})
    atmosphere2 = pyshdom.grid.resample_onto_grid(rte_grid, atmosphere)

    for wavelength in wavelengths:

        cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
                                                                mie_mono_table.radius.data,
                            size_distribution_function=size_distribution_function,particle_density=1.0,
                            reff={'coord_min':0.1, 'coord_max': 11.0, 'npoints': 100,
                            'spacing': 'logarithmic', 'units': 'micron'},
                            veff={'coord_min':0.09, 'coord_max': 0.11, 'npoints': 12,
                            'spacing': 'linear', 'units': 'unitless'}
                            )
        poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
        optical_properties = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table,
                                                                         exact_table=False)
        optical_properties['ssalb'][:,:,:] = ssalb
        extinction = np.zeros(optical_properties.extinction.shape)
        x,y,z = np.meshgrid(np.arange(extinction.shape[0]),
                       np.arange(extinction.shape[1]),
                       np.arange(extinction.shape[2]),indexing='ij')
        R = np.sqrt((x-x.mean())**2+ (y-y.mean())**2 + (z-z.mean())**2)
        #extinction[1:-1,1:-1,1:-1] = 0.0#*np.exp(-1*R/(0.2*R.mean()))[1:-1,1:-1,1:-1]
        extinction[2:-2,2:-2,2:-2] = ext
        extinction[index[0],index[1],index[2]] += step
        optical_properties['extinction'][:,:,:] = extinction
        cloud_poly_tables[wavelength] = poly_table
        config = pyshdom.configuration.get_config('../default_config.json')
        config['num_mu_bins'] = nmu
        config['num_phi_bins'] = nmu*2
        config['split_accuracy'] = split
        config['spherical_harmonics_accuracy'] = 0.0
        config['adapt_grid_factor'] = 10.0
        config['solution_accuracy'] = 1e-8
        config['high_order_radiance'] = True
        config['acceleration_flag'] = False
        config['deltam'] = True
        solver = pyshdom.solver.RTE(
                            numerical_params=config,
                            medium={'cloud': optical_properties},
                            source=pyshdom.source.solar(wavelength, solarmu, 0.0),
                            surface=pyshdom.surface.lambertian(albedo=surfacealb, ground_temperature=ground_temperature),
                            num_stokes=1,
                            atmosphere=atmosphere2,
                            name=None
                            )

        if load_solution is not None:
            solver.load_solution(load_solution)
        solvers.add_solver(wavelength,solver)

        Sensordict.get_measurements(solvers, maxiter=200, n_jobs=4, verbose=False)
        return solvers, Sensordict, cloud_poly_tables, step

class SolarJacobianThinNoSurface(TestCase):
    @classmethod
    def setUpClass(cls):
        tautol=0.2
        surfacealb=0.0
        ssalb = 1.0
        solarmu = 1.0
        reff=10.0
        resolutionfactor = 1
        nmu=16
        split=0.0
        ext=0.1
        veff=0.1
        step=1e-3
        ground_temperature=200.0
        mie_mono_table = pyshdom.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='../mie_tables',
                                                  verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor)

        Sensordict.add_uncertainty_model('MISR', pyshdom.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[0.86]
        solver = solvers[0.86]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        # CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        out = []
        for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
            print(i)
            data = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
                         resolution=resolutionfactor)
            data[1].add_uncertainty_model('MISR', pyshdom.uncertainties.NullUncertainty('L2'))
            for sensor in data[1]['MISR']['sensor_list']:
                data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
            rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
            # note only forward difference not central difference here, because I was testing
            # derivatives at ext=0.0 as well.
            out.append((rte_sensor_high[0.86].measurement_data[0].data - rte_sensor_ref.measurement_data[0].data)/step)
        finite_jacobian = np.stack(out, axis=0)
        np.save('./data/reference_{}_{}_{}_sfcalbedo_jacobian.npy'.format(surfacealb, step,ext), finite_jacobian)
        cls.jacobian_reference = np.load('./data/reference_{}_{}_{}_sfcalbedo_jacobian.npy'.format(surfacealb, step,ext))

        unknown_scatterers = pyshdom.containers.UnknownScatterers()
        unknown_scatterers.add_unknown('cloud', ['extinction'], cloud_poly_tables)
        unknown_scatterers.create_derivative_tables()
        solvers.add_microphysical_partial_derivatives(unknown_scatterers)
        solvers[0.86]._dext[np.where(solvers[0.86].medium['cloud'].extinction==0.0),0] = 0.0
        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = pyshdom.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':200,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        cost, gradient, jacobian = gradient_call()
        cls.jacobian = jacobian['jacobian_0.860'][0,0].data

    def test_jacobian(self):
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=5.52e-6))
