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

#use ground truth solvers and rte_grid because lazy.
#these are loaded but could be independently defined if desired.
sensor_dict, solvers, rte_grid = shdom.util.load_forward_model(file_name)

#define a copy of sensor_dict for storing synthetic measurements.
forward_sensors = shdom.script_util.make_forward_sensors(sensor_dict)
wavelengths = np.array(list(solvers.keys())) #strictly speaking these are keys may not always be wavelength.
                                             #if not then loop through solvers like a normal person.

#define all scattering properties. These could also be parsed directly from the medium objects in solvers.
#but lazy.
size_distribution_function = shdom.size_distribution.gamma

cloud_poly_tables = OrderedDict()
for wavelength in wavelengths:
    #mie_table path is hardcoded for now. . .
    mie_mono_table = shdom.mie.get_mono_table('Water',(wavelength,wavelength), relative_path='./mie_tables')
    cloud_size_distribution = shdom.size_distribution.get_size_distribution_grid(
                                                            mie_mono_table.radius.data,
                        size_distribution_function=size_distribution_function,particle_density=1.0,
                        reff=[4.0,25.0,25,'logarithmic','micron'],
                        veff=[0.09,0.11,2,'linear','unitless'],
                        )
    poly_table = shdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
    cloud_poly_tables[wavelength] = poly_table

#use a ground truth mask because lazy.
ground_truth = list(solvers.values())[0].medium[0]
mask = np.where(ground_truth.density.data > 0.0)

#make initial condition by copying some things from ground truth.
initial = ground_truth.copy(deep=True)
initial = initial.drop_vars(['table_index', 'extinction', 'ssalb','legcoef','stokes_index'])
initial = initial.reset_coords(['reff', 'veff'])

#ground truth reff/veff but a linear profile of density to initialize.
initial['density'][:,:,:] = np.repeat(np.repeat(np.append(np.append(np.zeros(4),np.linspace(0.0,1.0,21)),
                                                          np.zeros(2))[np.newaxis,np.newaxis,:],
                                                initial.density.shape[0],axis=0),initial.density.shape[1],axis=1)

#initialize containers for solvers and cloud optical scatterers.
#These are 'globals' that will be modified by set_state_fn
solvers_reconstruct = OrderedDict()
cloud_optical_scatterers = OrderedDict()

def set_state_fn(state):

    #update microphysics
    cloud_scatterer_on_rte_grid = shdom.grid.resample_onto_grid(rte_grid, initial)
    state_on_grid = np.zeros(cloud_scatterer_on_rte_grid.density.shape)
    state_on_grid[mask] = state
    cloud_scatterer_on_rte_grid.density[:,:,:] = state_on_grid

    #update optical properties
    for wavelength in wavelengths:
        poly_table = cloud_poly_tables[wavelength]
        cloud_optical_scatterers[wavelength] = shdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)

    #make new solver
    mediums = shdom.script_util.combine_to_medium([cloud_optical_scatterers])

    #lazily add ground truth rayleigh scatter.
    for key, medium in mediums.items():
        medium.append(solvers[key].medium[1])

    for key in solvers.keys():
        solvers_reconstruct[key] = shdom.solver.RTE(numerical_params=solvers[key].numerical_params,
                                        medium=mediums[key],
                                       source=solvers[key].source,
                                       surface=solvers[key].surface,
                                        num_stokes=solvers[key]._nstokes,
                                        name=None
                                       )

#SIMPLIFIED definition of unknown scatterers. These must match names in the scatterers.
unknown_scatterers = [('cloud', cloud_poly_tables,['density', 'reff'])]

def project_gradient_to_state(gradient):
    #for help writing this run the gradient e.g. shdom.gradient.levis_approx_uncorrelated_l2
    #on the initial state and examine the output and code of that function.
    #This is based on unknown scatterers and set_state_fn.

    #find the component of the gradient that is relevant.
    #NEW just select the names here based on what was chosen in unknown scatterers.
    #apply any change of spatial basis (in this case just a mask.)
    out_gradient = gradient.gradient.sel(variable_name='density', scatterer_name='cloud').data[mask[0],mask[1],mask[2]]
    return out_gradient.ravel(order='F')

#at the moment the bounds are specified unknown-by-unknown.
#this allows for the fact we might want coordinate dependent bounds on e.g. LWC. (adiabatic maximum)
#to stop some of the run away values that can happen.
min_bound = np.zeros(initial.density.data[mask].shape)
max_bound = np.zeros(initial.density.data[mask].shape) + 2

obj_fun = shdom.optimize.ObjectiveFunction.LevisApproxUncorrelatedL2(sensor_dict, solvers_reconstruct,
                                                                     forward_sensors, unknown_scatterers,
                                                                    set_state_fn, project_gradient_to_state,n_jobs=8,mpi_comm=None,
                                                                     verbose=True, maxiter=100, min_bounds=min_bound,
                                                                    max_bounds=max_bound)

optimizer = shdom.optimize.Optimizer(obj_fun)
optimizer._options['maxls'] = 30
optimizer._options['maxiter'] = 1 #maxiter to 1 to debug the saving of result.

#optimize for the specified initial condition.
result = optimizer.minimize(initial.density.data[mask])

#save the final synthetic forward model and measurements.
shdom.util.save_forward_model(save_name, forward_sensors, solvers_reconstruct)

to_save = xr.Dataset(
                data_vars={
                'final_state': ('state_dim',result.x),
                'intial_state': ('state_dim',initial.density.data[mask]),
                'input_name': file_name,
                'file_name': sys.argv[0], #optimization result is not reproducible without the script so save the script.
                'time': dt.datetime.now() #and the date.
                }
)
#append some extra stuff.
#This file coulbe be opened earlier and have callback output appended to it as well.
to_save.to_netcdf(save_name,'a',group='optimization_result')
shutil.copy(__file__, 'copied_script.py')
