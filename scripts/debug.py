#imports
import shdom
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import OrderedDict

os.chdir('/Users/jesserl2/Documents/Code/aviad_pyshdom_dev/pyshdom_dev/')
#shdom.util.set_pyshdom_path()
import sys

#load a cloud.
#locate the 'origin' of the cloud at 0.0,0.0 for simplicity.
#this option allows us to easily move individual clouds with respect to each other
#and even if overlapping they will be merged onto the RTE grid.
cloud_scatterer = shdom.grid.load_from_csv('./synthetic_cloud_fields/jpl_les/rico32x37x26.txt',
                                           density='lwc',origin=(0.0,0.0))

#load atmosphere file for rayleigh. (and eventually gases)
#'Altitude' coordinate is renamed to 'z'.
atmosphere = xr.open_dataset('./ancillary_data/AFGL_summer_mid_lat.nc')#.rename(Altitude='z')

#extract a chosen temperature_profile and the surface_pressure.
#only model atmosphere below 20 km. This choice needs to be made here so
#that an RTE grid can be defined.

#try where with drop.na
reduced_atmosphere = atmosphere.sel({'z': atmosphere.coords['z'].data[atmosphere.coords['z'].data <= 10.0]})

# -----  make the RTE grid ---------------------------

#make RTE grid just using cloud_scatterer for horizontal grid and 'merged' z coordinates.
merged_z_coordinate = shdom.grid.combine_z_coordinates([reduced_atmosphere,cloud_scatterer])

#simple 'union' horizontal grid merging for 3D and 1D needs to be fixed.
rte_grid = shdom.grid.make_grid(cloud_scatterer.x.data.min(),cloud_scatterer.x.data.max(),cloud_scatterer.x.data.size,
                           cloud_scatterer.y.data.min(),cloud_scatterer.y.data.max(),cloud_scatterer.y.data.size,
                           merged_z_coordinate)
#resample the cloud onto the rte_grid
cloud_scatterer_on_rte_grid = shdom.grid.resample_onto_grid(rte_grid, cloud_scatterer)

#define any necessary variables for microphysics here.
size_distribution_function = shdom.size_distribution.gamma

#We choose a gamma size distribution and therefore need to define a 'veff' variable.
cloud_scatterer_on_rte_grid['veff'] = xr.full_like(cloud_scatterer_on_rte_grid.reff, 0.1)

#shdom.grid.to_2parameter_lwc_file('./synthetic_cloud_fields/jpl_les/fixed_rico32x37x26.txt', cloud_scatterer_on_rte_grid, reduced_atmosphere)

#new = shdom.grid.load_2parameter_lwc_file('/Users/jesserl2/fixed_rico32x37x26.txt')

#This is modified by the user as needed.

#idealized monochromatic orthographic sensors at different wavelengths.
#9 'MISR-like' VIS cameras
#1 'MODIS-like' nadir multi-spectral sensor.

#TODO make this defined as a dict with entries for each 'instrument' to group sensors
Sensordict = OrderedDict()

misr_list = []

#add MISR-like sensors
sensor_zenith_list = [75.0,70.6]#,65.0,60.0,55.0,50.0,45.6,40.0,35.0,30.0,26.1,20.0,15.0,10.0,5.0]*2 + [0.0]
sensor_azimuth_list = [90]*2 #+ [-90]*15 +[0.0]

for zenith,azimuth in zip(sensor_zenith_list,sensor_azimuth_list):
    misr_list.append(
        shdom.sensor.add_sub_pixel_rays(shdom.sensor.orthographic_projection(1.65, cloud_scatterer,0.021,0.021, azimuth, zenith,
                                             altitude='TOA', stokes='I'
                                            ),FOV=0.0,degree=2)

    )

Sensordict['MISR'] = {'sensor_list': misr_list}

modis_list = []
#add MODIS-like sensors
wavelength_list = [1.65]
for wavelength in wavelength_list:
    modis_list.append(
       shdom.sensor.add_sub_pixel_rays(shdom.sensor.orthographic_projection(wavelength,cloud_scatterer,0.049,0.049,0.0,0.0,
                                            altitude='TOA',
                                            stokes='I'
                                            ),0.0,degree=2)#.sel(nrays=slice(6894,6901),npixels=1149)
    )

Sensordict['MODIS'] = {'sensor_list': modis_list}

def make_forward_sensors(sensors):

    forward_sensors = OrderedDict()
    for key,sensor in sensors.items():
        forward_sensor= OrderedDict()
        forward_sensor['sensor_list'] = [single_sensor.copy(deep=True) for single_sensor in sensor['sensor_list']]
        forward_sensors[key] = forward_sensor

    return forward_sensors


forward_sensors = make_forward_sensors(Sensordict)

#num_stokes should be set to choose whether to use num_stokes=USER_SPECIFIED
#even if only radiance needs to be simulated for accuracy reasons.
num_stokes_override_flag = False
num_stokes_override=3

#TODO hand specify num_stokes as in the configuration file, for example (or manually)
#Define the unique wavelengths from all sensors.

#extract all unique_wavelengths
#this treats even very slightly different wavelengths as unique.
wavelengths = shdom.script_util.get_unique_wavelengths(Sensordict)

names = OrderedDict()
surfaces = OrderedDict()
sources = OrderedDict()
numerical_parameters = OrderedDict()
num_stokes = OrderedDict()

config = shdom.configuration.get_config('./default_config.json')
config['spherical_harmonics_accuracy'] = 0.01

for wavelength in wavelengths:
    num_stokes[wavelength] = 3
    names[wavelength] = None
    surfaces[wavelength] = shdom.surface.lambertian(albedo=0.01) #surface is wavelength independent.
    sources[wavelength] = shdom.source.solar(145.0,0.0,solarflux=1.0)

    numerical_parameters[wavelength] = config#shdom.configuration.get_config('./default_config.json') #all use defaults.

#resample the cloud onto the rte_grid
cloud_scatterer_on_rte_grid = shdom.grid.resample_onto_grid(rte_grid, cloud_scatterer)

#define any necessary variables for microphysics here.
size_distribution_function = shdom.size_distribution.gamma
#We choose a gamma size distribution and therefore need to define a 'veff' variable.
cloud_scatterer_on_rte_grid['veff'] = (cloud_scatterer_on_rte_grid.reff.dims,
                                       np.full_like(cloud_scatterer_on_rte_grid.reff.data, fill_value=0.1))

#calculate the optical properties for this scatterer.
#All wavelengths use consistent settings.
cloud_poly_tables = OrderedDict()
cloud_optical_scatterers = OrderedDict()
for wavelength in wavelengths:
    print('making mie_table. . . may take a while.')
    mie_mono_table = shdom.mie.get_mono_table('Water',(wavelength,wavelength))
    cloud_size_distribution = shdom.size_distribution.get_size_distribution_grid(
                                                            mie_mono_table.radius.data,
                        size_distribution_function=size_distribution_function,particle_density=1.0,
                        reff=[4.0,25.0,25,'logarithmic','micron'],
                        veff=[0.09,0.11,2,'linear','unitless'],
                        )
    poly_table = shdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
    cloud_optical_scatterers[wavelength] = shdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)
    cloud_poly_tables[wavelength] = poly_table

#get rayleigh.
#This is self contained due to its simplicity.
rayleigh_scatterer_list = shdom.rayleigh.to_grid(wavelengths,atmosphere,rte_grid)


#TODO make this a function (combine scatterers to medium)
mediums = shdom.script_util.combine_to_medium([cloud_optical_scatterers, rayleigh_scatterer_list])
#group properties
# mediums = OrderedDict()
# for key,optical in cloud_optical_scatterers.items():

#     rayleigh = rayleigh_scatterer_list[key]
#     mediums[key] = [optical]

#make solver list
solvers = OrderedDict()

for key,name in names.items():
    solvers[key] = shdom.solver.RTE(numerical_params=numerical_parameters[key],
                                    medium=mediums[key],
                                   source=sources[key],
                                   surface=surfaces[key],
                                    num_stokes=num_stokes[key],
                                    name=name
                                   )
shdom.script_util.get_measurements(solvers, Sensordict, maxiter=2, n_jobs=10)

#resample the cloud onto the rte_grid
cloud_scatterer_on_rte_grid = shdom.grid.resample_onto_grid(rte_grid, cloud_scatterer)
cloud_scatterer_on_rte_grid['density'] = cloud_scatterer_on_rte_grid['density'] * 0.1

#define any necessary variables for microphysics here.
size_distribution_function = shdom.size_distribution.gamma
#We choose a gamma size distribution and therefore need to define a 'veff' variable.
cloud_scatterer_on_rte_grid['veff'] = (cloud_scatterer_on_rte_grid.reff.dims,
                                       np.full_like(cloud_scatterer_on_rte_grid.reff.data, fill_value=0.1))

#calculate the optical properties for this scatterer.
#All wavelengths use consistent settings.
cloud_poly_tables = OrderedDict()
cloud_optical_scatterers = OrderedDict()
for wavelength in wavelengths:
    print('making mie_table. . . may take a while.')
    mie_mono_table = shdom.mie.get_mono_table('Water',(wavelength,wavelength))
    cloud_size_distribution = shdom.size_distribution.get_size_distribution_grid(
                                                            mie_mono_table.radius.data,
                        size_distribution_function=size_distribution_function,particle_density=1.0,
                        reff=[4.0,25.0,25,'logarithmic','micron'],
                        veff=[0.09,0.11,2,'linear','unitless'],
                        )
    poly_table = shdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
    cloud_optical_scatterers[wavelength] = shdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)
    cloud_poly_tables[wavelength] = poly_table


#TODO make this a function (combine scatterers to medium)
mediums = shdom.script_util.combine_to_medium([cloud_optical_scatterers, rayleigh_scatterer_list])
#group properties
# mediums = OrderedDict()
# for key,optical in cloud_optical_scatterers.items():

#     rayleigh = rayleigh_scatterer_list[key]
#     mediums[key] = [optical]


#make solver list
solvers_reconstruct = OrderedDict()

for key,name in names.items():
    solvers_reconstruct[key] = shdom.solver.RTE(numerical_params=numerical_parameters[key],
                                    medium=mediums[key],
                                   source=sources[key],
                                   surface=surfaces[key],
                                    num_stokes=num_stokes[key],
                                    name=name
                                   )
cloud_unknowns = OrderedDict()
for key, scatterer in cloud_optical_scatterers.items():

    cloud_unknowns[key] = (scatterer, cloud_poly_tables[key], 'reff')

unordered_unknown_scatterers = [cloud_unknowns]
unknown_scatterers = shdom.script_util.combine_to_medium(unordered_unknown_scatterers)

table_derivatives = shdom.gradient.create_derivative_tables(solvers_reconstruct, unknown_scatterers)

n_jobs=10
verbose=True
init_solution=True
maxiter=2
exact_single_scatter=True
mpi_comm=None

for solver in solvers_reconstruct.values():
    if solver._srctype != 'S':
        raise NotImplementedError('Only Solar Source is supported for gradient calculations.')

shdom.script_util.parallel_solve(solvers_reconstruct, n_jobs=n_jobs, mpi_comm=mpi_comm,verbose=verbose,maxiter=maxiter,
    init_solution=init_solution)

    #These are called after the solution because they require at least _init_solution to be run.
if not hasattr(list(solvers_reconstruct.values())[0], '_direct_derivative_path'):
        #only needs to be called once.
        #If this doesn't work after 'solve' and requires only '_init_solution' then it will need to be
        #called inside parallel_solve. . .
    shdom.gradient.calculate_direct_beam_derivative(solvers_reconstruct)

solvers_reconstruct  = shdom.gradient.get_derivatives(solvers_reconstruct, table_derivatives)  #adds the _dext/_dleg/_dalb etc to the solvers.

    #prepare the sensors for the fortran subroutine for calculating gradient.
rte_sensors, sensor_mapping = shdom.script_util.sort_sensors(Sensordict, solvers_reconstruct, 'inverse')

gradient_function = shdom.gradient.grad_l2_old
out = shdom.gradient.parallel_gradient(solvers_reconstruct, rte_sensors, sensor_mapping, forward_sensors,
                                   gradient_fun=gradient_function,
                    mpi_comm=mpi_comm, n_jobs=n_jobs, exact_single_scatter=False)

print(out[-1].min(), out[-1].max())
