# This is an independent test to ensure the direct beam derivatives are correct.
# not yet complete.
from unittest import TestCase
from collections import OrderedDict
import numpy as np
import xarray as xr
import at3d
import sys

import warnings
warnings.filterwarnings('ignore')

import builtins as __builtin__
def print(*args, **kwargs):
    if '-vv' in sys.argv:
        return __builtin__.print(*args, **kwargs)

def cloud_direct_beam(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb, ground_temperature, step=0.0, index=(1,1,1),
          nmu=16, split=0.03, load_solution=None, resolution=1, boundary='open',deltam=True):
    rte_grid = at3d.grid.make_grid(0.05/resolution, 8*resolution, 0.05/resolution, 3*resolution,
                                      np.arange(0.1, 0.4, 0.05/resolution))

    grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
    rte_grid['density'] = (['x','y','z'], np.ones(grid_shape))
    rte_grid['reff'] = (['x','y','z'], np.zeros(grid_shape) + reff)
    rte_grid['veff'] = (['x','y','z'] ,np.zeros(grid_shape) + veff)

    #resample the cloud onto the rte_grid
    cloud_scatterer_on_rte_grid = at3d.grid.resample_onto_grid(rte_grid, rte_grid)

    #define any necessary variables for microphysics here.
    size_distribution_function = at3d.size_distribution.gamma


    wavelengths = np.array([0.86])

    cloud_poly_tables = OrderedDict()
    solvers = at3d.containers.SolversDict()
    temp = np.linspace(280,300, grid_shape[-1])[np.newaxis,np.newaxis,:]
    atmosphere = xr.Dataset(data_vars={
        'temperature': (['x','y','z'], np.repeat(np.repeat(temp,grid_shape[0], axis=0),grid_shape[1],axis=1))
    }, coords={'x':rte_grid.x, 'y': rte_grid.y, 'z': rte_grid.z})
    atmosphere2 = at3d.grid.resample_onto_grid(rte_grid, atmosphere)

    for wavelength in wavelengths:

        cloud_size_distribution = at3d.size_distribution.get_size_distribution_grid(
                                                                mie_mono_table.radius.data,
                            size_distribution_function=size_distribution_function,particle_density=1.0,
                            reff={'coord_min':0.1, 'coord_max': 11.0, 'npoints': 100,
                            'spacing': 'logarithmic', 'units': 'micron'},
                            veff={'coord_min':0.09, 'coord_max': 0.11, 'npoints': 12,
                            'spacing': 'linear', 'units': 'unitless'}
                            )
        poly_table = at3d.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
        optical_properties = at3d.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)
        optical_properties['ssalb'][:,:,:] = ssalb
        extinction = np.zeros(optical_properties.extinction.shape)
        x,y,z = np.meshgrid(np.arange(extinction.shape[0]),
                       np.arange(extinction.shape[1]),
                       np.arange(extinction.shape[2]),indexing='ij')
        R = np.sqrt((x-x.mean())**2+ (y-y.mean())**2 + (z-z.mean())**2)

        np.random.seed(1)
        ext = np.random.uniform(1e-6, 100, size=extinction.shape)
        extinction[:,:,:] = ext
        extinction[index[0],index[1],index[2]] += step
        optical_properties['extinction'][:,:,:] = extinction
        cloud_poly_tables[wavelength] = poly_table
        config = at3d.configuration.get_config('../default_config.json')
        config['num_mu_bins'] = nmu
        config['num_phi_bins'] = nmu*2
        config['split_accuracy'] = split
        config['spherical_harmonics_accuracy'] = 0.0
        config['adapt_grid_factor'] = 1.01
        config['solution_accuracy'] = 1e-8
        config['high_order_radiance'] = False
        config['acceleration_flag'] = False
        config['deltam'] = deltam
        config['x_boundary_condition'] = boundary
        config['y_boundary_condition'] = boundary
        solver = at3d.solver.RTE(
                            numerical_params=config,
                            medium={'cloud': optical_properties},
                            source=at3d.source.solar(wavelength, solarmu, 10.0),
                            surface=at3d.surface.lambertian(albedo=surfacealb, ground_temperature=ground_temperature),
                            num_stokes=1,
                            atmosphere=atmosphere2,
                            name=None
                            )

        if load_solution is not None:
            solver.load_solution(load_solution)
        solvers.add_solver(wavelength, solver)
        Sensordict = at3d.containers.SensorsDict()
        return solvers, Sensordict, cloud_poly_tables, step, rte_grid



class DirectBeamDerivativeDeltaMOpen(TestCase):
    @classmethod
    def setUpClass(cls):
        tautol=0.2
        surfacealb=1.0
        ssalb = 0.0
        solarmu = 0.3
        reff=10.0
        resolutionfactor = 1
        nmu=2
        split=0.0
        ext=10.0
        veff=0.1
        step=1e-3
        ground_temperature=200.0
        deltam = True
        boundary='open'
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data/',
                                                  verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_direct_beam(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor,deltam=True,
                                                                boundary=boundary)

        solver = solvers[0.86]
        solver._init_solution()
        solver._make_direct()
        reference = solver._dirflux[:solver._npts]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        #CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #
        #     data = cloud_direct_beam(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #                  resolution=resolutionfactor)
        #     solver2 = data[0][0.86]
        #     solver2._init_solution()
        #     solver2._make_direct()
        #     upper = solver2._dirflux[:solver._npts]
        #
        #     out.append((upper - reference)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/dirflux_gradient_{}_{}.npy'.format(deltam, boundary), finite_jacobian)
        cls.finite_jacobian = np.load('./data/dirflux_gradient_{}_{}.npy'.format(deltam, boundary))


        solver.calculate_direct_beam_derivative()
        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'extinction'
            )
        )
        forward_sensors = Sensordict.make_forward_sensors()


        # gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        # solvers, forward_sensors, unknown_scatterers,
        # parallel_solve_kwargs={'maxiter':800,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        # gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        # 'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        # cost, gradient, jacobian = gradient_call()

        #jacobian = jacobian['jacobian_0.860'][0,0].data
        solvers.calculate_microphysical_partial_derivatives(unknown_scatterers)
        #solver.calculate_direct_beam_derivative()

        solver.calculate_direct_beam_derivative()
        raygrad = np.zeros((solver._nstokes, solver._maxpg, 1))
        jacobian = []
        for i in range(solver._npts):
            raygrad = np.zeros((solver._nstokes, solver._maxpg, 1))
            raygrad = at3d.core.compute_direct_beam_deriv(
                dpath=solver._direct_derivative_path[:,i],
                dptr=solver._direct_derivative_ptr[:,i],
                numder=1,
                dextm=solver._dextm,
                transmit=1.0,
                abscell=1.0,
                inputweight=solver._dirflux[i],
                nstokes=1,
                longest_path_pts=solver._longest_path_pts,
                maxpg=solver._maxpg,
                raygrad=raygrad,
            )
            jacobian.append(raygrad[0,:,0])
        cls.jacobian = np.stack(jacobian, axis=1)

    def test_direct_beam_derivative(self):
        print('DirectBeamDerivativeDeltaMOpen', np.all(np.isfinite(self.jacobian)), np.all(np.isfinite(self.finite_jacobian)),
        np.max(np.abs(self.jacobian.ravel()-self.finite_jacobian.ravel())))
        print(self.jacobian.ravel())
        print('reference below:')
        print(self.finite_jacobian.ravel())
        # import pylab as py
        # py.figure()
        # py.plot(self.jacobian.ravel(), self.finite_jacobian.ravel(), 'x')
        # py.show()
        self.assertTrue(np.allclose(self.jacobian, self.finite_jacobian, atol=1e-5))

class DirectBeamDerivativeDeltaMPeriodic(TestCase):
    @classmethod
    def setUpClass(cls):
        tautol=0.2
        surfacealb=1.0
        ssalb = 0.0
        solarmu = 0.3
        reff=10.0
        resolutionfactor = 1
        nmu=2
        split=0.0
        ext=10.0
        veff=0.1
        step=1e-3
        ground_temperature=200.0
        deltam = True
        boundary='periodic'
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data/',
                                                  verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_direct_beam(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor,deltam=True,
                                                                boundary=boundary)

        solver = solvers[0.86]
        solver._init_solution()
        solver._make_direct()
        reference = solver._dirflux[:solver._npts]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        #CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #
        #     data = cloud_direct_beam(
        #         mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
        #         step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #         resolution=resolutionfactor, deltam=deltam, boundary=boundary)
        #     solver2 = data[0][0.86]
        #     solver2._init_solution()
        #     solver2._make_direct()
        #     upper = solver2._dirflux[:solver._npts]
        #
        #     out.append((upper - reference)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/dirflux_gradient_{}_{}.npy'.format(deltam, boundary), finite_jacobian)
        cls.finite_jacobian = np.load('./data/dirflux_gradient_{}_{}.npy'.format(deltam, boundary))


        solver.calculate_direct_beam_derivative()

        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'extinction'
            )
        )
        forward_sensors = Sensordict.make_forward_sensors()


        # gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        # solvers, forward_sensors, unknown_scatterers,
        # parallel_solve_kwargs={'maxiter':800,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        # gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        # 'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        #cost, gradient, jacobian = gradient_call()
        #self.solvers.calculate_microphysical_partial_derivatives(self.unknown_scatterers)
        #jacobian = jacobian['jacobian_0.860'][0,0].data

        solvers.calculate_microphysical_partial_derivatives(unknown_scatterers)
        solver.calculate_direct_beam_derivative()
        raygrad = np.zeros((solver._nstokes, solver._maxpg, 1))
        jacobian = []
        for i in range(solver._npts):
            raygrad = np.zeros((solver._nstokes, solver._maxpg, 1))
            raygrad = at3d.core.compute_direct_beam_deriv(
                dpath=solver._direct_derivative_path[:,i],
                dptr=solver._direct_derivative_ptr[:,i],
                numder=1,
                dextm=solver._dextm,
                transmit=1.0,
                abscell=1.0,
                inputweight=solver._dirflux[i],
                nstokes=1,
                longest_path_pts=solver._longest_path_pts,
                maxpg=solver._maxpg,
                raygrad=raygrad,
            )
            jacobian.append(raygrad[0,:,0])
        cls.jacobian = np.stack(jacobian, axis=1)

    def test_direct_beam_derivative(self):
        self.assertTrue(np.allclose(self.jacobian, self.finite_jacobian, atol=1e-5))

class DirectBeamDerivativePeriodic(TestCase):
    @classmethod
    def setUpClass(cls):
        tautol=0.2
        surfacealb=1.0
        ssalb = 0.0
        solarmu = 0.3
        reff=10.0
        resolutionfactor = 1
        nmu=2
        split=0.0
        ext=10.0
        veff=0.1
        step=1e-3
        ground_temperature=200.0
        deltam = False
        boundary='periodic'
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data/',
                                                  verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_direct_beam(
            mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
            step=0.0,nmu=nmu,split=split,
            resolution=resolutionfactor,deltam=deltam,
            boundary=boundary)

        solver = solvers[0.86]
        solver._init_solution()
        solver._make_direct()
        reference = solver._dirflux[:solver._npts]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        #CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #
        #     data = cloud_direct_beam(
        #         mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
        #         step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #         resolution=resolutionfactor, deltam=deltam, boundary=boundary)
        #     solver2 = data[0][0.86]
        #     solver2._init_solution()
        #     solver2._make_direct()
        #     upper = solver2._dirflux[:solver._npts]
        #
        #     out.append((upper - reference)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/dirflux_gradient_{}_{}.npy'.format(deltam, boundary), finite_jacobian)
        cls.finite_jacobian = np.load('./data/dirflux_gradient_{}_{}.npy'.format(deltam, boundary))


        solver.calculate_direct_beam_derivative()
        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'extinction'
            )
        )
        forward_sensors = Sensordict.make_forward_sensors()


        # gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        # solvers, forward_sensors, unknown_scatterers,
        # parallel_solve_kwargs={'maxiter':800,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        # gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        # 'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        # cost, gradient, jacobian = gradient_call()
        # jacobian = jacobian['jacobian_0.860'][0,0].data
        solvers.calculate_microphysical_partial_derivatives(unknown_scatterers)
        solver.calculate_direct_beam_derivative()
        raygrad = np.zeros((solver._nstokes, solver._maxpg, 1))
        jacobian = []
        for i in range(solver._npts):
            raygrad = np.zeros((solver._nstokes, solver._maxpg, 1))
            raygrad = at3d.core.compute_direct_beam_deriv(
                dpath=solver._direct_derivative_path[:,i],
                dptr=solver._direct_derivative_ptr[:,i],
                numder=1,
                dextm=solver._dextm,
                transmit=1.0,
                abscell=1.0,
                inputweight=solver._dirflux[i],
                nstokes=1,
                longest_path_pts=solver._longest_path_pts,
                maxpg=solver._maxpg,
                raygrad=raygrad,
            )
            jacobian.append(raygrad[0,:,0])
        cls.jacobian = np.stack(jacobian, axis=1)

    def test_direct_beam_derivative(self):
        self.assertTrue(np.allclose(self.jacobian, self.finite_jacobian, atol=1e-5))

class DirectBeamDerivativeOpen(TestCase):
    @classmethod
    def setUpClass(cls):
        tautol=0.2
        surfacealb=1.0
        ssalb = 0.0
        solarmu = 0.3
        reff=10.0
        resolutionfactor = 1
        nmu=2
        split=0.0
        ext=10.0
        veff=0.1
        step=1e-3
        ground_temperature=200.0
        deltam = False
        boundary='open'
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data/',
                                                  verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_direct_beam(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor,deltam=True,
                                                                boundary=boundary)

        solver = solvers[0.86]
        solver._init_solution()
        solver._make_direct()
        reference = solver._dirflux[:solver._npts]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        #CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #
        #     data = cloud_direct_beam(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #                  resolution=resolutionfactor)
        #     solver2 = data[0][0.86]
        #     solver2._init_solution()
        #     solver2._make_direct()
        #     upper = solver2._dirflux[:solver._npts]
        #
        #     out.append((upper - reference)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/dirflux_gradient_{}_{}.npy'.format(deltam, boundary), finite_jacobian)
        cls.finite_jacobian = np.load('./data/dirflux_gradient_{}_{}.npy'.format(deltam, boundary))


        solver.calculate_direct_beam_derivative()
        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'extinction'
            )
        )
        forward_sensors = Sensordict.make_forward_sensors()


        # gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        # solvers, forward_sensors, unknown_scatterers,
        # parallel_solve_kwargs={'maxiter':800,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        # gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        # 'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        # cost, gradient, jacobian = gradient_call()
        # jacobian = jacobian['jacobian_0.860'][0,0].data

        solvers.calculate_microphysical_partial_derivatives(unknown_scatterers)
        solver.calculate_direct_beam_derivative()
        raygrad = np.zeros((solver._nstokes, solver._maxpg, 1))
        jacobian = []
        for i in range(solver._npts):
            raygrad = np.zeros((solver._nstokes, solver._maxpg, 1))
            raygrad = at3d.core.compute_direct_beam_deriv(
                dpath=solver._direct_derivative_path[:,i],
                dptr=solver._direct_derivative_ptr[:,i],
                numder=1,
                dextm=solver._dextm,
                transmit=1.0,
                abscell=1.0,
                inputweight=solver._dirflux[i],
                nstokes=1,
                longest_path_pts=solver._longest_path_pts,
                maxpg=solver._maxpg,
                raygrad=raygrad,
            )
            jacobian.append(raygrad[0,:,0])
        cls.jacobian = np.stack(jacobian, axis=1)

    def test_direct_beam_derivative(self):
        print('DirectBeamDerivativeOpen', np.all(np.isfinite(self.jacobian)), np.all(np.isfinite(self.finite_jacobian)),
        np.max(np.abs(self.jacobian.ravel()-self.finite_jacobian.ravel())))
        self.assertTrue(np.allclose(self.jacobian, self.finite_jacobian, atol=1e-5))
