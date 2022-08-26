from unittest import TestCase
from collections import OrderedDict
import numpy as np
import xarray as xr
import at3d
from scipy import stats
import sys
import warnings
warnings.filterwarnings('ignore')

import builtins as __builtin__
def print(*args, **kwargs):
    if '-vv' in sys.argv:
        return __builtin__.print(*args, **kwargs)

class PlanckDerivative(TestCase):

    @classmethod
    def setUpClass(cls):
        np.random.seed(1)
        units_all = ('B', 'T', 'R')
        size = 1000
        temp_test = np.random.uniform(0.0, 900.0, size=size)
        wavelens = np.random.uniform(0.0, 20.0, size=size)
        wavenos = np.zeros((size, 2))
        wavenos[:,0] = np.random.uniform(0.0, 10000.0, size=size)
        wavenos[:,1] = wavenos[:,0]+ np.random.uniform(0.0, 300.0 ,size)

        finite_diff_step = 0.15
        cls.all_ref = []
        cls.all_finite_diff = []
        cls.all_ref_grad = []

        for units in units_all:
            finite_diff_grad = []
            planck_ref = []
            planck_upper = []
            ref_grad = []
            for temp, wavelen, waveno in zip(temp_test,wavelens, wavenos):
                planck_ref.append(at3d.core.planck_function(temp=temp, units=units, waveno=waveno, wavelen=wavelen))
                planck_upper.append(at3d.core.planck_function(temp=temp+finite_diff_step, units=units, waveno=waveno, wavelen=wavelen))

                ref_grad.append(at3d.core.planck_derivative(temp=temp, units=units, waveno=waveno, wavelen=wavelen))
            ref_grad = np.array(ref_grad)
            planck_ref = np.array(planck_ref)
            planck_upper = np.array(planck_upper)
            finite_diff_grad = ((planck_upper - planck_ref)/finite_diff_step)
            cls.all_ref.append(planck_ref)
            cls.all_finite_diff.append(finite_diff_grad)
            cls.all_ref_grad.append(ref_grad)

    def testTemperatureUnits(self):
        self.assertTrue(np.allclose(self.all_ref_grad[1], 1.0))

    def testBandUnits(self):
        good_data = self.all_ref_grad[0][np.where(~np.isnan(self.all_ref_grad[0]))]
        good_finite = self.all_finite_diff[0][np.where(~np.isnan(self.all_ref_grad[0]))]
        maxerror = np.max(np.abs(good_data - good_finite))
        self.assertTrue(maxerror < 1.7e-3)

    def testBadBandUnits(self):
        self.assertTrue(np.all(self.all_ref[0][np.where(np.isnan(self.all_ref_grad[0]))] < 1e-9))

    def testBadRadianceUnits(self):
        self.assertTrue(np.all(self.all_ref[-1][np.where(np.isnan(self.all_ref_grad[-1]))] < 1e-9))

    def testRadianceUnits(self):
        good_data = self.all_ref_grad[-1][np.where(~np.isnan(self.all_ref_grad[-1]))]
        good_finite = self.all_finite_diff[-1][np.where(~np.isnan(self.all_ref_grad[-1]))]
        maxerror = np.max(np.abs(good_data - good_finite))
        self.assertTrue(maxerror < 1.23e-2)



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
        gradout, cost, ierr, errmsg = at3d.core.update_costfunction(
            cost=cost,
            gradout=gradout,
            stokesout=stokesout,
            measurement=measurement,
            raygrad_pixel=raygrad_pixel,
            uncertainties=uncertainties,
            costfunc=costfunc,
        )
        at3d.checks.check_errcode(ierr, errmsg)
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

        gradout, cost, ierr, errmsg = at3d.core.update_costfunction(
            cost=cost,
            gradout=gradout,
            stokesout=stokesout,
            measurement=measurement,
            raygrad_pixel=raygrad_pixel,
            uncertainties=uncertainties,
            costfunc=costfunc,
        )
        at3d.checks.check_errcode(ierr, errmsg)
        cls.gradout = gradout
        cls.cost = cost
    def test_cost(self):
        self.assertAlmostEqual(self.cost, 6519.21, places=2)
    def test_gradient(self):
        self.assertAlmostEqual(self.gradout[0, 0, 0], 45329.43, places=2)

def legacy_perspective_projection(wavelength, fov, x_resolution, y_resolution,
                           position_vector, lookat_vector, up_vector,
                           stokes='I', sub_pixel_ray_args={'method':None}):
    """
    A SEPARATE IMPLEMENTATION AFTER A BUGFIX WHICH PRODUCES THE EXACT
    GEOMETRY USED FOR THE FINITE DIFFERENCE LINEARIZATIONS HERE.

    Generates a sensor dataset that observes a target location with
    a perspective (pinhole camera) projection.

    Parameters
    ----------
    wavelength: float,
        Wavelength in [micron]
    fov: float
        Field of view [deg]
    x_resolution: int
        Number of pixels in camera x axis
    y_resolution: int
        Number of pixels in camera y axis
    position_vector: list of 3 float elements
        [x , y , z] which are:
        Location in global x coordinates [km] (North)
        Location in global y coordinates [km] (East)
        Location in global z coordinates [km] (Up)
    lookat_vector: list of 3 float elements
        [x , y , z] which are:
        Point in global x coordinates [km] (North) where the camera is pointing at
        Point in global y coordinates [km] (East) where the camera is pointing at
        Point in global z coordinates [km] (Up) where the camera is pointing at
    up_vector: list of 3 float elements
        The up vector determines the roll of the camera.
    stokes: list or string
       list or string of stokes components to observe ['I', 'Q', 'U', 'V'].
    sub_pixel_ray_args : dict
        dictionary defining the method for generating sub-pixel rays. The callable
        which generates the position_perturbations and weights (e.g. at3d.sensor.gaussian)
        should be set as the 'method', while arguments to that callable, should be set as
        other entries in the dict. Each argument have two values, one for each of the
        x and y axes of the image plane, respectively.
        E.g. sub_pixel_ray_args={'method':at3d.sensor.gaussian, 'degree': (2, 3)}

    Returns
    -------
    sensor : xr.Dataset
        A dataset containing all of the information required to define a sensor
        for which synthetic measurements can be simulated;
        positions and angles of all pixels, sub-pixel rays and their associated weights,
        and the sensor's observables.

    """
    norm = lambda x: x / np.linalg.norm(x, axis=0)

    #assert samples>=1, "Sample per pixel is an integer >= 1"
    #assert int(samples) == samples, "Sample per pixel is an integer >= 1"

    assert int(x_resolution) == x_resolution, "x_resolution is an integer >= 1"
    assert int(y_resolution) == y_resolution, "y_resolution is an integer >= 1"

    # The bounding_box is not nessesary in the prespactive projection, but we still may consider
    # to use if we project the rays on the bounding_box when the differences in mu , phi angles are below certaine precision.
    #     if(bounding_box is not None):

    #         xmin, ymin, zmin = bounding_box.x.data.min(),bounding_box.y.data.min(),bounding_box.z.data.min()
    #         xmax, ymax, zmax = bounding_box.x.data.max(),bounding_box.y.data.max(),bounding_box.z.data.max()

    nx = x_resolution
    ny = y_resolution
    position = np.array(position_vector, dtype=np.float32)
    lookat = np.array(lookat_vector, dtype=np.float32)
    up = np.array(up_vector)
    direction = lookat - position

    zaxis = norm(direction)
    xaxis = norm(np.cross(up, zaxis))
    yaxis = np.cross(zaxis, xaxis)
    rotation_matrix = np.stack((xaxis, yaxis, zaxis), axis=1)

    M = max(nx, ny)
    npix = nx*ny
    R = np.array([nx, ny])/M # R will be used to scale the sensor meshgrid.
    dy = 2*R[1]/ny # pixel length in y direction in the normalized image plane.
    dx = 2*R[0]/nx # pixel length in x direction in the normalized image plane.

    # THIS DIFFERS FROM sensor.perspective_projection HERE.
    x_s, y_s, z_s = np.meshgrid(np.linspace(-R[0], R[0]-dx, nx),
                                np.linspace(-R[1], R[1]-dy, ny), 1.0)

    # Here x_c, y_c, z_c coordinates on the image plane before transformation to the requaired observation angle
    focal = 1.0 / np.tan(np.deg2rad(fov) / 2.0) # focal (normalized) length when the sensor size is 2 e.g. r in [-1,1).
    fov_x = np.rad2deg(2*np.arctan(R[0]/focal))
    fov_y = np.rad2deg(2*np.arctan(R[1]/focal))

    k = np.array([[focal, 0, 0],
                  [0, focal, 0],
                  [0, 0, 1]], dtype=np.float32)
    inv_k = np.linalg.inv(k)

    homogeneous_coordinates = np.stack([x_s.ravel(), y_s.ravel(), z_s.ravel()])

    x_c, y_c, z_c = norm(np.matmul(
        rotation_matrix, np.matmul(inv_k, homogeneous_coordinates)))
    # Here x_c, y_c, z_c coordinates on the image plane after transformation to the requaired observation

    # x,y,z mu, phi in the global coordinates:
    mu = -z_c.astype(np.float64)
    phi = (np.arctan2(y_c, x_c) + np.pi).astype(np.float64)
    x = np.full(npix, position[0], dtype=np.float32)
    y = np.full(npix, position[1], dtype=np.float32)
    z = np.full(npix, position[2], dtype=np.float32)

    image_shape = [nx,ny]
    sensor = at3d.sensor.make_sensor_dataset(x.ravel(), y.ravel(), z.ravel(),
                                 mu.ravel(), phi.ravel(), stokes, wavelength)
    # compare to orthographic projection, prespective projection may not have bounding box.
    #     if(bounding_box is not None):
    #         sensor['bounding_box'] = xr.DataArray(np.array([xmin,ymin,zmin,xmax,ymax,zmax]),
    #                                               coords={'bbox': ['xmin','ymin','zmin','xmax','ymax','zmax']},dims='bbox')

    sensor['image_shape'] = xr.DataArray(image_shape,
                                         coords={'image_dims': ['nx', 'ny']},
                                         dims='image_dims')
    sensor.attrs = {
        'projection': 'Perspective',
        'fov_deg': fov,
        'fov_x_deg': fov_x,
        'fov_y_deg': fov_y,
        'x_resolution': x_resolution,
        'y_resolution': y_resolution,
        'position': position,
        'lookat': lookat,
        'rotation_matrix': rotation_matrix,
        'sensor_to_camera_transform_matrix':k

    }

    if sub_pixel_ray_args['method'] is not None:

        #generate the weights and perturbations to the pixel positions in the image plane.
        sub_pixel_ray_method, subpixel_ray_kwargs_x, subpixel_ray_kwargs_y =  \
                            at3d.sensor._parse_sub_pixel_ray_args(sub_pixel_ray_args)
        position_perturbations_x, weights_x = sub_pixel_ray_method(x_s.size,
                                                                   **subpixel_ray_kwargs_x)
        position_perturbations_y, weights_y = sub_pixel_ray_method(y_s.size,
                                                                   **subpixel_ray_kwargs_y)

        #merge the two dimensions
        perturbations_x = np.repeat(position_perturbations_x[..., np.newaxis]*dx/2.0,
                                    position_perturbations_y.shape[-1], axis=-1)
        perturbations_y = np.repeat(position_perturbations_y[..., np.newaxis, :]*dy/2.0,
                                    position_perturbations_x.shape[-1], axis=-2)
        big_weightx = np.repeat(weights_x[..., np.newaxis], weights_y.shape[-1], axis=-1)
        big_weighty = np.repeat(weights_y[..., np.newaxis, :], weights_x.shape[-1], axis=-2)

        #apply perturbations to original image plane coordinates.
        x_ray = (x_s.ravel()[:, np.newaxis, np.newaxis] + perturbations_x).ravel()
        y_ray = (y_s.ravel()[:, np.newaxis, np.newaxis] + perturbations_y).ravel()
        z_ray = np.repeat(np.repeat(z_s.ravel()[:, np.newaxis, np.newaxis],
                                    perturbations_x.shape[-2], axis=-2),
                          perturbations_y.shape[-1], axis=-1).ravel()
        ray_homogeneous = np.stack([x_ray, y_ray, z_ray])

        x_c, y_c, z_c = norm(np.matmul(
            rotation_matrix, np.matmul(inv_k, ray_homogeneous)))
        # Here x_c, y_c, z_c coordinates on the image plane after transformation to the requaired observation

        # x,y,z mu, phi in the global coordinates:
        mu = -z_c.astype(np.float64)
        phi = (np.arctan2(y_c, x_c) + np.pi).astype(np.float64)
        x = np.full(x_c.size, position[0], dtype=np.float32)
        y = np.full(x_c.size, position[1], dtype=np.float32)
        z = np.full(x_c.size, position[2], dtype=np.float32)

        #make the pixel indices and ray weights.
        pixel_index = np.repeat(np.repeat(range(len(sensor.cam_mu.data)),
                                          weights_x.shape[-1]), weights_y.shape[-1])
        ray_weight = (big_weightx*big_weighty).ravel()
        #update ray variables to sensor dataset.
        sensor['ray_mu'] = ('nrays', mu)
        sensor['ray_phi'] = ('nrays', phi)
        sensor['ray_x'] = ('nrays', x)
        sensor['ray_y'] = ('nrays', y)
        sensor['ray_z'] = ('nrays', z)
        sensor['pixel_index'] = ('nrays', pixel_index)
        sensor['ray_weight'] = ('nrays', ray_weight)
        sensor['use_subpixel_rays'] = True

        sub_pixel_ray_args['method'] = sub_pixel_ray_args['method'].__name__
        for attribute in sub_pixel_ray_args:
            sensor.attrs['sub_pixel_ray_args_{}'.format(attribute)] = sub_pixel_ray_args[attribute]
    else:
            #duplicate ray variables to sensor dataset.
        sensor = at3d.sensor._add_null_subpixel_rays(sensor)
    return sensor


#A function for computing clouds for the thermal jacobian reference case.
def cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=0.0, index=(1,1,1),
          nmu=16, split=0.03, load_solution=None, resolution=1):

    rte_grid = at3d.grid.make_grid(0.05/resolution, 7*resolution, 0.05/resolution, 7*resolution,
                                      np.arange(0.1, 0.5, 0.05/resolution))

    grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
    rte_grid['density'] = (['x','y','z'], np.ones(grid_shape))
    rte_grid['reff'] = (['x','y','z'], np.zeros(grid_shape) + reff)
    rte_grid['veff'] = (['x','y','z'] ,np.zeros(grid_shape) + veff)

    #resample the cloud onto the rte_grid
    cloud_scatterer_on_rte_grid = at3d.grid.resample_onto_grid(rte_grid, rte_grid)

    #define any necessary variables for microphysics here.
    size_distribution_function = at3d.size_distribution.gamma

    #define sensors.
    Sensordict = at3d.containers.SensorsDict()
    zeniths =  [75.0,60.0,45.6,26.1]*2 + [0.0] +  [75.0,60.0,45.6,26.1]*6
    azimuths =  [90.0]*4 + [-90]*4 + [0.0] + [0]*4 + [180]*4 + [-45]*4 + [135]*4+ [45]*4 + [-135]*4

    target = [rte_grid.x.data[-1]/2, rte_grid.y.data[-1]/2, rte_grid.z.data[rte_grid.z.size//2]]
    R = 10.0
    xs = target[0] + R*np.sin(np.deg2rad(zeniths))*np.cos(np.deg2rad(azimuths))
    ys = target[1] + R*np.sin(np.deg2rad(zeniths))*np.sin(np.deg2rad(azimuths))
    zs = target[2] + R*np.cos(np.deg2rad(zeniths))
    for x,y,z in zip(xs,ys,zs):
        Sensordict.add_sensor('MISR',
        legacy_perspective_projection(11.0, 5.0, 26, 26,[x,y,z], target,
                                              [0,1,0],stokes=['I'])
                             )

    wavelengths = Sensordict.get_unique_solvers()

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
        extinction[1:-1,1:-1,1:-1] = ext*np.exp(-1*R/(0.2*R.mean()))[1:-1,1:-1,1:-1]
        extinction[index[0],index[1],index[2]] += step
        optical_properties['extinction'][:,:,:] = extinction

        cloud_poly_tables[wavelength] = poly_table
        config = at3d.configuration.get_config('../default_config.json')
        config['num_mu_bins'] = nmu
        config['num_phi_bins'] = nmu*2
        config['split_accuracy'] = split
        config['spherical_harmonics_accuracy'] = 0.0
        config['adapt_grid_factor'] = 1000.0
        config['solution_accuracy'] = 1e-7
        config['high_order_radiance'] = False
        config['acceleration_flag'] = True
        config['deltam'] = False
        config['tautol'] = 0.2
        config['transcut'] = 5e-5
        solver = at3d.solver.RTE(
                            numerical_params=config,
                            medium={'cloud': optical_properties},
                            source=at3d.source.thermal(wavelength),
                            surface=at3d.surface.lambertian(albedo=surfacealb,
                            ground_temperature=ground_temperature),
                            num_stokes=1,
                            atmosphere=atmosphere2,
                            name=None
                            )
        if load_solution is not None:
            solver.load_solution(load_solution)
        solvers.add_solver(wavelength,solver)

    Sensordict.get_measurements(solvers, maxiter=200, n_jobs=1, verbose=False)
    return solvers, Sensordict, cloud_poly_tables, step, rte_grid

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
        mie_mono_table = at3d.mie.get_mono_table('Water',(11.0,11.0),
                                              max_integration_radius=65.0,
                                              minimum_effective_radius=0.1,
                                              relative_dir='../mie_tables',
                                              verbose=False)
        mie_mono_table.to_netcdf('./data/mie_table_11micron.nc')
        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor)
        Sensordict.add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[11.0]
        solver = solvers[11.0]

        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 11.0)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'extinction'
            )
        )

        # deriv_gen = at3d.medium.OpticalGenerator(rte_grid,'cloud', 11.0)
        # deriv_info = OrderedDict()
        #
        # unknown_scatterers = at3d.containers.UnknownScatterers()
        # unknown_scatterers.add_unknowns(['extinction'], deriv_gen)
        # unknown_scatterers = at3d.containers.UnknownScatterers()
        # unknown_scatterers.add_unknown('cloud', ['extinction'], cloud_poly_tables)
        # unknown_scatterers.create_derivative_tables()
        # solvers.add_microphysical_partial_derivatives(unknown_scatterers)
        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
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
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #
        #     data = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=-1*step,index=(a,b,c),nmu=nmu,load_solution=saved, split=0.0,
        #                  resolution=resolutionfactor)
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_low, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #     out.append((rte_sensor_high[11.0].measurement_data[0].data - rte_sensor_low[11.0].measurement_data[0].data)/(2*step))#[a,b,c])
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/reference_noscat_nosurface_jacobian.npy', finite_jacobian)
        cls.jacobian_reference = np.load('./data/reference_noscat_nosurface_jacobian.npy')

    def test_jacobian(self):
        print('thermal no surface', np.all(np.isfinite(self.jacobian)), np.all(np.isfinite(self.jacobian_reference)),
        np.max(np.abs(self.jacobian.ravel()-self.jacobian_reference.ravel())))
        print(np.max(np.abs(self.jacobian.ravel()-self.jacobian_reference.ravel())/np.abs(self.jacobian_reference.ravel())))
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
        mie_mono_table = at3d.mie.get_mono_table('Water',(11.0,11.0),
                                              max_integration_radius=65.0,
                                              minimum_effective_radius=0.1,
                                              relative_dir='./data',
                                              verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor)
        Sensordict.add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[11.0]
        solver = solvers[11.0]

        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 11.0)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'extinction'
            )
        )

        # deriv_gen = at3d.medium.OpticalGenerator(rte_grid,'cloud', 11.0)
        # deriv_info = OrderedDict()
        #
        # unknown_scatterers = at3d.containers.UnknownScatterers()
        # unknown_scatterers.add_unknowns(['extinction'], deriv_gen)

        # unknown_scatterers = at3d.containers.UnknownScatterers()
        # unknown_scatterers.add_unknown('cloud', ['extinction'], cloud_poly_tables)
        # unknown_scatterers.create_derivative_tables()
        # solvers.add_microphysical_partial_derivatives(unknown_scatterers)
        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
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
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #
        #     data = cloud(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=-1*step,index=(a,b,c),nmu=nmu,load_solution=saved, split=0.0,
        #                  resolution=resolutionfactor)
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_low, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #     out.append((rte_sensor_high[11.0].measurement_data[0].data - rte_sensor_low[11.0].measurement_data[0].data)/(2*step))#[a,b,c])
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/reference_noscat_surface_jacobian.npy', finite_jacobian)
        cls.jacobian_reference = np.load('./data/reference_noscat_surface_jacobian.npy')

    def test_jacobian(self):
        print('thermal with surface', np.all(np.isfinite(self.jacobian)), np.all(np.isfinite(self.jacobian_reference)),
        np.max(np.abs(self.jacobian.ravel()-self.jacobian_reference.ravel())))
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=8e-3))


def cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb, ground_temperature, step=0.0, index=(1,1,1),
          nmu=16, split=0.03, load_solution=None, resolution=1, random=False, random_ssalb=False,
          random_reff=False,
          perturb='extinct', solve=True, deltam=True):
    np.random.seed(1)
    rte_grid = at3d.grid.make_grid(0.05/resolution, 9*resolution, 0.05/resolution, 9*resolution,
                                      np.arange(0.1, 0.65, 0.05/resolution))
    if random_reff:
        reff = np.random.uniform(0.1, 11.0, size=grid_shape)
    grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
    rte_grid['density'] = (['x','y','z'], np.ones(grid_shape))
    rte_grid['reff'] = (['x','y','z'], np.zeros(grid_shape) + reff)
    rte_grid['veff'] = (['x','y','z'] ,np.zeros(grid_shape) + veff)

    #resample the cloud onto the rte_grid
    cloud_scatterer_on_rte_grid = at3d.grid.resample_onto_grid(rte_grid, rte_grid)

    #define any necessary variables for microphysics here.
    size_distribution_function = at3d.size_distribution.gamma

    #define sensors.
    Sensordict = at3d.containers.SensorsDict()
    zeniths =  [75.0,60.0,45.6,26.1]*2 + [0.0] +  [75.0,60.0,45.6,26.1]*6
    azimuths =  [90.0]*4 + [-90]*4 + [0.0] + [0]*4 + [180]*4 + [-45]*4 + [135]*4+ [45]*4 + [-135]*4

    target = [rte_grid.x.data[-1]/2, rte_grid.y.data[-1]/2, rte_grid.z.data[rte_grid.z.size//2]]
    R = 10.0
    xs = target[0] + R*np.sin(np.deg2rad(zeniths))*np.cos(np.deg2rad(azimuths))
    ys = target[1] + R*np.sin(np.deg2rad(zeniths))*np.sin(np.deg2rad(azimuths))
    zs = target[2] + R*np.cos(np.deg2rad(zeniths))
    for x,y,z in zip(xs,ys,zs):
        Sensordict.add_sensor('MISR',
        legacy_perspective_projection(0.86, 4.0, 13, 13,[x,y,z], target,
                                              [0,1,0],stokes=['I'])
                             )

    sensor = at3d.sensor.make_sensor_dataset(np.array([0.15]),np.array([0.055]),
                                              np.array([0.44999998807907104]),np.array([1.0]),
                                             np.array([0.0]),
                                             wavelength=0.86,stokes=['I'],fill_ray_variables=True)

    Sensordict.add_sensor('MISR', sensor)
    wavelengths = Sensordict.get_unique_solvers()

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
        if not deltam:
            optical_properties.legcoef[:,1:] = 0.0
        ssalb = np.zeros(optical_properties.extinction.shape) + ssalb
        if random_ssalb:
            ssalb = np.random.uniform(2e-3,1.0, size=optical_properties.extinction.shape)
        if perturb == 'ssalb':
            ssalb[index[0],index[1],index[2]] += step
        ssalb[np.where(ssalb < 0.0)] = 0.0
        ssalb[np.where(ssalb > 1.0)] = 1.0
        optical_properties['ssalb'][:] = ssalb
        extinction = np.zeros(optical_properties.extinction.shape)
        x,y,z = np.meshgrid(np.arange(extinction.shape[0]),
                       np.arange(extinction.shape[1]),
                       np.arange(extinction.shape[2]),indexing='ij')
        R = np.sqrt((x-x.mean())**2+ (y-y.mean())**2 + (z-z.mean())**2)
        #extinction[1:-1,1:-1,1:-1] = 0.0#*np.exp(-1*R/(0.2*R.mean()))[1:-1,1:-1,1:-1]
        if random:

            ext = np.random.uniform(1e-6, ext, size=extinction[2:-2,2:-2,2:-2].shape)
        extinction[2:-2,2:-2,2:-2] = ext
        if perturb == 'extinct':
            extinction[index[0],index[1],index[2]] += step
        optical_properties['extinction'][:,:,:] = extinction

        if perturb == 'g':
            new_legcoef = xr.concat([optical_properties.legcoef.copy(deep=True), optical_properties.legcoef.copy(deep=True)], dim='table_index')
            new_legcoef[0,1,1] += step
            optical_properties['legcoef'] = (['stokes_index', 'legendre_index', 'table_index'],new_legcoef.data)
            optical_properties['table_index'][:,index[0],index[1],index[2]] = 2
        #print('after',optical_properties.legcoef[0,1,:])
        cloud_poly_tables[wavelength] = poly_table
        config = at3d.configuration.get_config('../default_config.json')
        config['num_mu_bins'] = nmu
        config['num_phi_bins'] = nmu*2
        config['split_accuracy'] = split
        config['spherical_harmonics_accuracy'] = 0.0
        config['adapt_grid_factor'] = 10.0
        config['solution_accuracy'] = 1e-8
        config['high_order_radiance'] = True
        config['acceleration_flag'] = False
        config['deltam'] = deltam
        config['tautol'] = 0.2
        config['transcut'] = 5e-5
        solver = at3d.solver.RTE(
                            numerical_params=config,
                            medium={'cloud': optical_properties},
                            source=at3d.source.solar(wavelength, solarmu, 0.0),
                            surface=at3d.surface.lambertian(albedo=surfacealb, ground_temperature=ground_temperature),
                            num_stokes=1,
                            atmosphere=atmosphere2,
                            name=None
                            )

        if load_solution is not None:
            solver.load_solution(load_solution)
        solvers.add_solver(wavelength,solver)
        if solve:
            Sensordict.get_measurements(solvers, maxiter=200, n_jobs=4, verbose=False)
        return solvers, Sensordict, cloud_poly_tables, step, rte_grid


class SolarJacobianThinNoSurfaceExtinction(TestCase):
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
        step=-1e-3
        ground_temperature=200.0
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data',
                                                  verbose=False)
        mie_mono_table.to_netcdf('./data/mie_table_860nm.nc')
        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor, random=True,
                                                                               random_ssalb=True, perturb='extinct')

        Sensordict.add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[0.86]
        solver = solvers[0.86]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        # CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #     data = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #                  resolution=resolutionfactor, random=True, random_ssalb=True, perturb='extinct')
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #     # note only forward difference not central difference here, because I was testing
        #     # derivatives at ext=0.0 as well.
        #     out.append((rte_sensor_high[0.86].measurement_data[0].data - rte_sensor_ref.measurement_data[0].data)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/thin_reference_ext_{}_{}_sfcalbedo_jacobian.npy'.format(surfacealb, step), finite_jacobian)
        cls.jacobian_reference = np.load('./data/thin_reference_ext_{}_{}_sfcalbedo_jacobian.npy'.format(surfacealb, step))

        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'extinction'
            )
        )


        # deriv_gen = at3d.medium.OpticalGenerator(rte_grid,'cloud', 0.86)
        #
        # unknown_scatterers = at3d.containers.UnknownScatterers()
        # unknown_scatterers.add_unknowns(['extinction'], deriv_gen)

        # unknown_scatterers = at3d.containers.UnknownScatterers()
        # unknown_scatterers.add_unknown('cloud', ['extinction'], cloud_poly_tables)
        # unknown_scatterers.create_derivative_tables()
        # solvers.add_microphysical_partial_derivatives(unknown_scatterers)
        # solvers[0.86]._dext[np.where(solvers[0.86].medium['cloud'].extinction==0.0),0] = 0.0
        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':200,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        cost, gradient, jacobian = gradient_call()
        cls.jacobian = jacobian['jacobian_0.860'][0,0].data

    def test_jacobian(self):
        cond = np.where(self.jacobian != 0.0)
        res = stats.linregress(self.jacobian_reference[cond],self.jacobian[cond])
        print(np.max(np.abs(self.jacobian_reference.ravel()-self.jacobian.ravel())),
        np.sqrt(np.sum((self.jacobian_reference.ravel()-self.jacobian.ravel())**2))/np.sqrt(np.sum(self.jacobian_reference**2)),
        res.intercept, res.slope, res.rvalue**2)
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=9.2e-6))#8.32e-6))

class SolarJacobianThinNoSurfaceAlbedo(TestCase):
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
        step=-1e-3
        ground_temperature=200.0
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data',
                                                  verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor, random=True,
                                                                random_ssalb=True, perturb='ssalb')

        Sensordict.add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[0.86]
        solver = solvers[0.86]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        # CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #     data = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #                  resolution=resolutionfactor, random=True, random_ssalb=True, perturb='ssalb')
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #     # note only forward difference not central difference here, because I was testing
        #     # derivatives at ext=0.0 as well.
        #     out.append((rte_sensor_high[0.86].measurement_data[0].data - rte_sensor_ref.measurement_data[0].data)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/thin_reference_ssalb_0.0_-0.001_sfcalbedo_jacobian.npy'.format(surfacealb, step), finite_jacobian)
        cls.jacobian_reference = np.load('./data/thin_reference_ssalb_0.0_-0.001_sfcalbedo_jacobian.npy'.format(surfacealb, step))

        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'ssalb'
            )
        )

        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':200,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        cost, gradient, jacobian = gradient_call()
        cls.jacobian = jacobian['jacobian_0.860'][0,0].data

    def test_jacobian(self):
        cond = np.where(self.jacobian != 0.0)
        res = stats.linregress(self.jacobian_reference[cond],self.jacobian[cond])
        print(np.max(np.abs(self.jacobian_reference.ravel()-self.jacobian.ravel())),
        np.sqrt(np.sum((self.jacobian_reference.ravel()-self.jacobian.ravel())**2))/np.sqrt(np.sum(self.jacobian_reference**2)),
        res.intercept, res.slope, res.rvalue**2)
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=5.3e-7))#5e-7))

class SolarJacobianThinNoSurfaceAsymmetry(TestCase):
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
        step=-1e-1
        ground_temperature=200.0
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data',
                                                  verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor, random=True,
                                                                random_ssalb=True, perturb='g')

        Sensordict.add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[0.86]
        solver = solvers[0.86]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        # CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #     data = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #                  resolution=resolutionfactor, random=True, random_ssalb=True, perturb='g')
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #     # note only forward difference not central difference here, because I was testing
        #     # derivatives at ext=0.0 as well.
        #     out.append((rte_sensor_high[0.86].measurement_data[0].data - rte_sensor_ref.measurement_data[0].data)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/thin_reference_g_0.0_-0.1_sfcalbedo_jacobian.npy'.format(surfacealb, step), finite_jacobian)
        cls.jacobian_reference = np.load('./data/thin_reference_g_{}_{}_sfcalbedo_jacobian.npy'.format(surfacealb, step))

        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'legendre_0_1'
            )
        )

        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':200,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        cost, gradient, jacobian = gradient_call()
        cls.jacobian = jacobian['jacobian_0.860'][0,0].data

    def test_jacobian(self):
        cond = np.where(self.jacobian != 0.0)
        res = stats.linregress(self.jacobian_reference[cond],self.jacobian[cond])
        print(np.max(np.abs(self.jacobian_reference.ravel()-self.jacobian.ravel())),
        np.sqrt(np.sum((self.jacobian_reference.ravel()-self.jacobian.ravel())**2))/np.sqrt(np.sum(self.jacobian_reference**2)),
        res.intercept, res.slope, res.rvalue**2)
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=1.3e-6))#3.8e-7))



class SolarJacobianThinNoSurfaceExtinctionNoDeltaM(TestCase):
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
        step=-1e-3
        ground_temperature=200.0
        deltam = False
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data',
                                                  verbose=False)
        mie_mono_table.to_netcdf('./data/mie_table_860nm.nc')
        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor, random=True,
                                                                               random_ssalb=True, perturb='extinct',
                                                                               deltam=deltam)

        Sensordict.add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[0.86]
        solver = solvers[0.86]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        # CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #     data = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #                  resolution=resolutionfactor, random=True, random_ssalb=True, perturb='extinct',deltam=deltam)
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #
        #     out.append((rte_sensor_high[0.86].measurement_data[0].data - rte_sensor_ref.measurement_data[0].data)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/thin_reference_ext_{}_{}_sfcalbedo_jacobian_nodeltam.npy'.format(surfacealb, step), finite_jacobian)
        cls.jacobian_reference = np.load('./data/thin_reference_ext_{}_{}_sfcalbedo_jacobian_nodeltam.npy'.format(surfacealb, step))


        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'extinction'
            )
        )

        # deriv_gen = at3d.medium.OpticalGenerator(rte_grid,'cloud', 0.86)
        #
        # unknown_scatterers = at3d.containers.UnknownScatterers()
        # unknown_scatterers.add_unknowns(['extinction'], deriv_gen)

        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':200,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        cost, gradient, jacobian = gradient_call()
        cls.jacobian = jacobian['jacobian_0.860'][0,0].data

    def test_jacobian(self):
        cond = np.where(self.jacobian != 0.0)
        res = stats.linregress(self.jacobian_reference[cond],self.jacobian[cond])
        print(np.max(np.abs(self.jacobian_reference.ravel()-self.jacobian.ravel())),
        np.sqrt(np.sum((self.jacobian_reference.ravel()-self.jacobian.ravel())**2))/np.sqrt(np.sum(self.jacobian_reference**2)),
        res.intercept, res.slope, res.rvalue**2)
        # import pylab as py
        # py.figure()
        # py.plot(self.jacobian.ravel(), self.jacobian_reference.ravel(), 'x')
        # py.show()
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=2.3e-5))#9.2e-6))

class SolarJacobianThinNoSurfaceAlbedoNoDeltaM(TestCase):
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
        step=-1e-3
        ground_temperature=200.0
        deltam = False
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data',
                                                  verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor, random=True,
                                                                random_ssalb=True, perturb='ssalb',deltam=deltam)

        Sensordict.add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[0.86]
        solver = solvers[0.86]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        # CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #     data = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #                  resolution=resolutionfactor, random=True, random_ssalb=True, perturb='ssalb',deltam=deltam)
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #     # note only forward difference not central difference here, because I was testing
        #     # derivatives at ext=0.0 as well.
        #     out.append((rte_sensor_high[0.86].measurement_data[0].data - rte_sensor_ref.measurement_data[0].data)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/thin_reference_ssalb_{}_{}_sfcalbedo_jacobian_nodeltam.npy'.format(surfacealb, step), finite_jacobian)
        cls.jacobian_reference = np.load('./data/thin_reference_ssalb_{}_{}_sfcalbedo_jacobian_nodeltam.npy'.format(surfacealb, step))

        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'ssalb'
            )
        )

        # deriv_gen = at3d.medium.OpticalGenerator(rte_grid,'cloud', 0.86)
        #
        # unknown_scatterers = at3d.containers.UnknownScatterers()
        # unknown_scatterers.add_unknowns(['ssalb'], deriv_gen)
        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':200,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        cost, gradient, jacobian = gradient_call()
        cls.jacobian = jacobian['jacobian_0.860'][0,0].data

    def test_jacobian(self):
        cond = np.where(self.jacobian != 0.0)
        res = stats.linregress(self.jacobian_reference[cond],self.jacobian[cond])
        print(np.max(np.abs(self.jacobian_reference.ravel()-self.jacobian.ravel())),
        np.sqrt(np.sum((self.jacobian_reference.ravel()-self.jacobian.ravel())**2))/np.sqrt(np.sum(self.jacobian_reference**2)),
        res.intercept, res.slope, res.rvalue**2)
        # import pylab as py
        # py.figure()
        # py.plot(self.jacobian.ravel(), self.jacobian_reference.ravel(), 'x')
        # py.show()
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=1.25e-6))#5e-7))

class SolarJacobianThinNoSurfaceAsymmetryNoDeltaM(TestCase):
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
        step=-1e-1
        ground_temperature=200.0
        deltam=False
        mie_mono_table = at3d.mie.get_mono_table('Water',(0.86,0.86),
                                                  max_integration_radius=65.0,
                                                  minimum_effective_radius=0.1,
                                                  relative_dir='./data',
                                                  verbose=False)

        solvers, Sensordict,cloud_poly_tables,final_step,rte_grid = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature,
                                                                 step=0.0,nmu=nmu,split=split,
                                                                resolution=resolutionfactor, random=True,
                                                                random_ssalb=True, perturb='g',deltam=deltam)

        Sensordict.add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        for sensor in Sensordict['MISR']['sensor_list']:
            Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        rte_sensor_ref, mapping_ref = Sensordict.sort_sensors(solvers, measurements=Sensordict)
        rte_sensor_ref = rte_sensor_ref[0.86]
        solver = solvers[0.86]

        indices_for_jacobian = np.where(solver.medium['cloud'].extinction.data > 0.0)
        # CODE FOR GENERATING THE FINITE DIFFERENCE REFERENCE.
        # out = []
        # for i,(a,b,c) in enumerate(zip(*indices_for_jacobian)):
        #     print(i)
        #     data = cloud_solar(mie_mono_table,ext,veff,reff,ssalb,solarmu,surfacealb,ground_temperature, step=step,index=(a,b,c),nmu=nmu,load_solution=None, split=split,
        #                  resolution=resolutionfactor, random=True, random_ssalb=True, perturb='g',deltam=deltam)
        #     data[1].add_uncertainty_model('MISR', at3d.uncertainties.NullUncertainty('L2'))
        #     for sensor in data[1]['MISR']['sensor_list']:
        #         data[1]['MISR']['uncertainty_model'].calculate_uncertainties(sensor)
        #     rte_sensor_high, mapping = data[1].sort_sensors(solvers, measurements=data[1])
        #     # note only forward difference not central difference here, because I was testing
        #     # derivatives at ext=0.0 as well.
        #     out.append((rte_sensor_high[0.86].measurement_data[0].data - rte_sensor_ref.measurement_data[0].data)/step)
        # finite_jacobian = np.stack(out, axis=0)
        # np.save('./data/thin_reference_g_{}_{}_sfcalbedo_jacobian_nodeltam.npy'.format(surfacealb, step), finite_jacobian)
        cls.jacobian_reference = np.load('./data/thin_reference_g_{}_{}_sfcalbedo_jacobian_nodeltam.npy'.format(surfacealb, step))

        deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
        unknown_scatterers = at3d.containers.UnknownScatterers(
            at3d.medium.UnknownScatterer(
                deriv_gen, 'legendre_0_1'
            )
        )

        # deriv_gen = at3d.medium.OpticalGenerator(rte_grid,'cloud', 0.86)
        #
        # unknown_scatterers = at3d.containers.UnknownScatterers()
        # unknown_scatterers.add_unknowns(['legendre_0_1'], deriv_gen)
        forward_sensors = Sensordict.make_forward_sensors()

        gradient_call = at3d.gradient.LevisApproxGradientUncorrelated(Sensordict,
        solvers, forward_sensors, unknown_scatterers,
        parallel_solve_kwargs={'maxiter':200,'n_jobs':4, 'setup_grid':False, 'verbose': False, 'init_solution':False},
        gradient_kwargs={'exact_single_scatter': True, 'cost_function': 'L2',
        'indices_for_jacobian': indices_for_jacobian}, uncertainty_kwargs={'add_noise': False})
        cost, gradient, jacobian = gradient_call()
        cls.jacobian = jacobian['jacobian_0.860'][0,0].data

    def test_jacobian(self):
        cond = np.where(self.jacobian != 0.0)
        res = stats.linregress(self.jacobian_reference[cond],self.jacobian[cond])
        print(np.max(np.abs(self.jacobian_reference.ravel()-self.jacobian.ravel())),
        np.sqrt(np.sum((self.jacobian_reference.ravel()-self.jacobian.ravel())**2))/np.sqrt(np.sum(self.jacobian_reference**2)),
        res.intercept, res.slope, res.rvalue**2)
        self.assertTrue(np.allclose(self.jacobian.ravel(), self.jacobian_reference.ravel(), atol=1.3e-6))#3.8e-7))





# class AdjointSource(TestCase):
#
#     @classmethod
#     def setUpClass(cls):
#
#         np.random.seed(1)
#
#         nx=npx=10
#         ny=npy=24
#         npz=nz=30
#         xstart=0.0
#         ystart=0.0
#         bcflag=0
#         gridtype='P'
#         ipflag=0
#         delx=0.02
#         dely=0.02
#         zlevels=np.linspace(0.0,1.0,nz)
#
#         def ibits(val, bit, ret_val):
#             if val & 2 ** bit:
#                 return ret_val
#             return 0
#
#         nx1, ny1 = nx + 1, ny + 1
#         if bcflag & 5 or ibits(ipflag, 0, 1):
#             nx1 -= 1
#         if bcflag & 7 or ibits(ipflag, 1, 1):
#             ny1 -= 1
#         nbpts = nx1 * ny1 * nz
#
#         xgrid, ygrid, zgrid = at3d.core.new_grids(
#             bcflag=bcflag,
#             gridtype=gridtype,
#             npx=npx,
#             npy=npy,
#             nx=nx,
#             ny=ny,
#             nz=nz,
#             xstart=xstart,
#             ystart=ystart,
#             delxp=delx,
#             delyp=dely,
#             zlevels=zlevels
#         )
#
#         npts, ncells, gridpos, gridptr, neighptr, \
#         treeptr, cellflags = at3d.core.init_cell_structure(
#             maxig=2.0*nbpts,
#             maxic=2.0*2*nbpts,
#             bcflag=bcflag,
#             ipflag=ipflag,
#             nx=nx,
#             ny=ny,
#             nz=nz,
#             nx1=nx1,
#             ny1=ny1,
#             xgrid=xgrid,
#             ygrid=ygrid,
#             zgrid=zgrid
#         )
#         number_tests = 500
#         lefts = np.zeros(number_tests)
#         rights = np.zeros(number_tests)
#         for i in range(number_tests):
#             x0=np.random.uniform(low=0.0, high=xgrid.max())
#             y0=np.random.uniform(low=0.0, high=ygrid.max())
#             z0=np.random.uniform(low=0.2, high=zgrid.max())
#             mu=np.random.uniform(low=-1.0, high=1.0)
#             phi=np.random.uniform(low=0.0, high=np.pi)
#             adjoint_magnitude=np.random.uniform(low=0.0,high=1.0)
#             total_ext = np.random.uniform(low=0.0,high=40.0,size=npts)
#             field = np.random.uniform(low=0.0,high=40.0, size=npts)
#             adjoint_source=np.zeros(total_ext.shape)
#             transmit=1.0
#             xe=ye=ze=0.0
#
#             adjoint_source, side, xe,ye,ze, transmit, ierr, errmsg = at3d.core.pencil_beam_prop(x0=x0,y0=y0,z0=z0,
#              transmit=transmit,
#              xe=xe,
#              ye=ye,
#              ze=ze,
#              dirflux=adjoint_source,
#              tautol=0.2,
#              bcflag=bcflag,
#              ipflag=ipflag,
#              magnitude=adjoint_magnitude,
#              mu2=mu,
#              phi2=phi,
#              total_ext=total_ext,
#              nx=nx,
#              ny=ny,
#              nz=nz,
#              ncells=ncells,
#              npts=npts,
#              cellflags=cellflags[:ncells],
#              xgrid=xgrid,
#              ygrid=ygrid,
#              zgrid=zgrid,
#              gridpos=gridpos[:,:npts],
#              gridptr=gridptr[:,:ncells],
#              neighptr=neighptr[:,:ncells],
#              treeptr=treeptr[:,:ncells]
#              )
#
#             magnitude = at3d.core.transmission_integral(x0=x0,y0=y0,z0=z0,
#                                          bcflag=bcflag,
#                                          ipflag=ipflag,
#                                          mu2=mu,
#                                          phi2=phi,
#                                          field=field,
#                                          nx=nx,
#                                          ny=ny,
#                                          nz=nz,
#                                          ncells=ncells,
#                                          total_ext=total_ext,
#                                          npts=npts,
#                                          cellflags=cellflags[:ncells],
#                                          xgrid=xgrid,
#                                          ygrid=ygrid,
#                                          zgrid=zgrid,
#                                          gridpos=gridpos[:,:npts],
#                                          gridptr=gridptr[:,:ncells],
#                                          neighptr=neighptr[:,:ncells],
#                                          treeptr=treeptr[:,:ncells]
#                                               )
#             lefts[i] = magnitude*adjoint_magnitude
#             rights[i] = np.dot(adjoint_source, field*total_ext)
#         cls.lefts = lefts
#         cls.rights = rights
#     def dot_product_test(self):
#         self.assertTrue(np.allclose(self.lefts, self.rights))
