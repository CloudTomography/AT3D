#imports
import shdom
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import OrderedDict
import shutil
import datetime as dt

os.chdir('/Users/jesserl2/Documents/Code/aviad_pyshdom_dev/pyshdom_dev/')
#shdom.util.set_pyshdom_path()
import sys

file_name = sys.argv[1]
save_name = sys.argv[2]

#Load forward model. Some inputs from the forward problem are used.
#but when using real measurements nothing from 'solvers' will be available.
#'sensor_dict' would also only supply 'cam' variables not 'ray' variables
#in that case.
sensor_dict, solvers, rte_grid =shdom.util.load_forward_model(file_name)

#define an approximate forward_sensor that doesn't have any subpixel rays.
# forward_sensors = shdom.organization.SensorsDict()
# for instrument in sensor_dict:
#     sensor_list = sensor_dict[instrument]['sensor_list']
#     for sensor in sensor_list:
#         approximated_sensor = shdom.sensor.make_sensor_dataset(sensor.cam_x,
#                                                                sensor.cam_y,
#                                                                sensor.cam_z,
#                                                                sensor.cam_mu,
#                                                                sensor.cam_phi,
#                                                                stokes=['I'],
#                                                                wavelength=sensor.wavelength,
#                                                                fill_ray_variables=True)
#         forward_sensors.add_sensor(instrument, sensor)

forward_sensors = sensor_dict.make_forward_sensors() # make a perfect copy of the sensor with the same subpixel rays.
wavelengths = forward_sensors.get_unique_solvers()

#define all scattering properties. In general, these will be defined here rather than read from forward problem.
size_distribution_function = shdom.size_distribution.gamma

cloud_poly_tables = OrderedDict() # keeping poly_tables for later use.
for wavelength in wavelengths:
    #mie_table path is hardcoded for now. . .
    mie_mono_table = shdom.mie.get_mono_table('Water',(wavelength,wavelength), relative_path='./mie_tables')
    cloud_size_distribution = shdom.size_distribution.get_size_distribution_grid(
                                                            mie_mono_table.radius.data,
                        size_distribution_function=size_distribution_function,particle_density=1.0,
                        reff=[9.0,11.0,13,'linear','micron'], #reff/veff are not going to be optimized for
                        veff=[0.04,0.6,13,'linear','unitless'],#so only small number of entries in the table.
                        )
    poly_table = shdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
    cloud_poly_tables[wavelength] = poly_table

#make an rte_grid
# rte_grid = shdom.grid.make_grid(0.0, _.x.data.max(), 25, 0.0, _.y.data.max(), 25,
#                     np.linspace(0.4, 1.45, 25)) #simple z coordinate as no atmosphere

density = np.zeros(26)
density[1:-1] = 0.001

initial = xr.Dataset(
        data_vars={
        'density': (['z'], density),
        'reff': (['z'], np.ones(26)*10.0),
        'veff': ('z', np.ones(26)*0.05)
        },
        coords={
        'z': rte_grid.z.data
        }
)
initial_on_grid = shdom.grid.resample_onto_grid(rte_grid, initial)

#apply a mask here if desired.


#initialize containers for solvers
#This 'globals' that will be modified by set_state_fn during the iterations
#of the reconstruction.
solvers_reconstruct = shdom.organization.SolversDict()

def set_state_fn(state):
    print(state.min(), state.max())
    print(state)
    #update microphysics
    cloud_scatterer_on_rte_grid = initial_on_grid.copy(deep=True)
    #cloud_scatterer_on_rte_grid.density[:,:,:] = state.reshape(cloud_scatterer_on_rte_grid.density.shape)

    #update optical properties
    for wavelength in wavelengths:
        poly_table = cloud_poly_tables[wavelength] #tables are used here.
        cloud_optical_scatterer = shdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)
        extinction = np.zeros((cloud_optical_scatterer.extinction.shape))

        reshaped_state = state.reshape((extinction.shape[0]-2,
                                                    extinction.shape[1]-2,
                                                    extinction.shape[2]-2))
        #reshaped_state[np.where(reshaped_state <=0.0)] = 0.001
        cloud_optical_scatterer.extinction[1:-1,1:-1,1:-1] = reshaped_state
        #decide on numerical accuracy. Here a reduced accuracy is used compared
        #to the forward problem (as is typical when nature is supplying the input measurements.)
        config = shdom.configuration.get_config('./default_config.json')
        config['num_mu_bins'] = 8
        config['num_phi_bins'] = 16
        config['split_accuracy'] = 0.1
        config['spherical_harmonics_accuracy'] = 0.0
        config['solution_accuracy'] = 1e-4
        solvers_reconstruct.add_solver( wavelength,
                                shdom.solver.RTE(numerical_params=config,
                                                medium={'cloud': cloud_optical_scatterer},
                                                source=solvers[wavelength].source,
                                                surface=solvers[wavelength].surface,
                                                num_stokes=solvers[wavelength]._nstokes,
                                                name=None
                                                )
                                )

unknown_scatterers = shdom.organization.UnknownScatterers()
unknown_scatterers.add_unknown('cloud', 'extinction', cloud_poly_tables)
unknown_scatterers.create_derivative_tables()

def project_gradient_to_state(gradient):
    #for help writing this run the gradient e.g. shdom.gradient.levis_approx_uncorrelated_l2
    #on the initial state and examine the output and code of that function.
    #This is based on unknown scatterers and set_state_fn.

    #find the component of the gradient that is relevant.
    #apply any change of spatial basis (in this case there is none.)
    out_gradient = gradient.gradient.sel(variable_name='extinction', scatterer_name='cloud').data[1:-1,1:-1,1:-1]#[mask[0],mask[1],mask[2]]
    print(out_gradient.ravel())
    return out_gradient.ravel()

#at the moment the bounds are specified unknown-by-unknown.
#this allows for the fact we might want coordinate dependent bounds on e.g. LWC. (adiabatic maximum)
#to stop some of the run away values that can happen.
min_bound = np.zeros(initial_on_grid.density[1:-1,1:-1,1:-1].size) + 0.001
max_bound = np.zeros(initial_on_grid.density[1:-1,1:-1,1:-1].size) + 300

obj_fun = shdom.optimize.ObjectiveFunction.LevisApproxUncorrelatedL2(sensor_dict, solvers_reconstruct,
                                                                     forward_sensors, unknown_scatterers,
                                                                    set_state_fn, project_gradient_to_state,n_jobs=8,mpi_comm=None,
                                                                     verbose=True, maxiter=100, min_bounds=min_bound,
                                                                    max_bounds=max_bound)

optimizer = shdom.optimize.Optimizer(obj_fun)
optimizer._options['maxls'] = 30
optimizer._options['maxiter'] = 10 #maxiter to 1 to debug the saving of result.
#optimize for the specified initial condition.
result = optimizer.minimize(0.1*np.ones(initial_on_grid.density.data[1:-1,1:-1,1:-1].size))

#save the final synthetic forward model and measurements.
shdom.util.save_forward_model(save_name, forward_sensors, solvers_reconstruct)

to_save = xr.Dataset(
                data_vars={
                'final_state': ('state_dim',result.x),
                'intial_state': ('state_dim',0.1*np.ones(initial_on_grid.density.data[1:-1,1:-1,1:-1].size)),
                'input_name': file_name,
                'file_name': sys.argv[0], #optimization result is not reproducible without the script so save the script.
                'time': dt.datetime.now() #and the date.
                }
)
#append some extra stuff.
#This file coulbe be opened earlier and have callback output appended to it as well.
to_save.to_netcdf(save_name,'a',group='optimization_result')
shutil.copy(__file__, 'copied_script.py')
