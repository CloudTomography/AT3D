import shdom
import numpy as np
import os
import sys
os.chdir('/Users/jesserl2/Documents/Code/aviad_pyshdom_dev/pyshdom_dev/')

file_name = sys.argv[-1]
#load the cloud.
cloud_scatterer = pyshdom.grid.load_from_csv('./synthetic_cloud_fields/jpl_les/rico32x37x26.txt',
                                           density='lwc',origin=(0.0,0.0))
#load atmosphere
#atmosphere = xr.open_dataset('./ancillary_data/AFGL_summer_mid_lat.nc')
#subset the atmosphere, choose only the bottom kilometre.
#reduced_atmosphere = atmosphere.sel({'z': atmosphere.coords['z'].data[atmosphere.coords['z'].data <= 1.0]})
#merge the atmosphere and cloud z coordinates
merged_z_coordinate = pyshdom.grid.combine_z_coordinates([cloud_scatterer])

#make a merged grid for the rte.
rte_grid = pyshdom.grid.make_grid(cloud_scatterer.x[1]-cloud_scatterer.x[0],cloud_scatterer.x.size,
                           cloud_scatterer.y[1]-cloud_scatterer.y[0],cloud_scatterer.y.size,
                           merged_z_coordinate)

#resample the cloud onto the rte_grid
cloud_scatterer_on_rte_grid = pyshdom.grid.resample_onto_grid(rte_grid, cloud_scatterer)

#define any necessary variables for microphysics here.
size_distribution_function = pyshdom.size_distribution.gamma
#We choose a gamma size distribution and therefore need to define a 'veff' variable.
cloud_scatterer_on_rte_grid['veff'] = (cloud_scatterer_on_rte_grid.reff.dims,
                                       np.full_like(cloud_scatterer_on_rte_grid.reff.data, fill_value=0.1))

#define sensors.
sensors_dict = pyshdom.util.SensorsDict()

sensor_zenith_list =  [75.0,60.0,45.6,26.1]*2 + [0.0]
sensor_azimuth_list = [90]*4 + [-90]*4 +[0.0]
for zenith,azimuth in zip(sensor_zenith_list,sensor_azimuth_list):
    sensors_dict.add_sensor('MISR',
                    pyshdom.sensor.orthographic_projection(0.86, cloud_scatterer,0.02,0.02,
                                                         azimuth, zenith,
                                                         altitude='TOA', stokes=['I'],
                                                         sub_pixel_ray_args={'method':pyshdom.sensor.gaussian,
                                                         'degree':2})
    )
for wavelength in (1.65, 2.2, 3.7):
    sensors_dict.add_sensor('MODIS',
                    pyshdom.sensor.orthographic_projection(wavelength,
                                                         cloud_scatterer,
                                                         0.04,0.04, 0.0,0.0,altitude='TOA',
                                                         stokes=['I'],
                                                         sub_pixel_ray_args={'method':pyshdom.sensor.gaussian,
                                                         'degree':4})
    )

wavelengths = sensors_dict.get_unique_solvers()
min_stokes = sensors_dict.get_minimum_stokes()
solvers = pyshdom.util.SolversDict()

#rayleigh optical properties if desired.
rayleigh_scatterer_list = pyshdom.rayleigh.to_grid(wavelengths,atmosphere,rte_grid)

for wavelength in wavelengths:
    #optical properties from mie calculations.
    mie_mono_table = pyshdom.mie.get_mono_table('Water', (wavelength, wavelength), relative_dir='./mie_tables')
    cloud_size_distribution = pyshdom.size_distribution.get_size_distribution_grid(
                                                            mie_mono_table.radius.data,
                        size_distribution_function=size_distribution_function,particle_density=1.0,
                        reff=[4.0,25.0,25,'logarithmic','micron'],
                        veff=[0.09,0.11,2,'linear','unitless'],
                        )
    poly_table = pyshdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
    cloud_optical_scatterer = pyshdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)
    config = pyshdom.configuration.get_config('./default_config.json')
    # config['num_mu_bins'] = 16
    # config['num_phi_bins'] = 32
    # config['split_accuracy'] = 0.03
    # config['spherical_harmonics_accuracy'] = 0.0
    # config['solution_accuracy'] = 1e-5
    solvers.add_solver(wavelength, pyshdom.solver.RTE(numerical_params=config,
                                    medium={'cloud': cloud_optical_scatterer },#, 'rayleigh': rayleigh_scatterer_list[wavelength]},
                                   source=pyshdom.source.solar(-1*np.cos(np.deg2rad(40.0)),0.0,solarflux=1.0),
                                   surface=pyshdom.surface.lambertian(albedo=0.04),
                                    num_stokes=min_stokes[wavelength],
                                    name=None
                                   ))

#import mpi4py.MPI as MPI
#mpi_comm = MPI.COMM_WORLD
#print('I am {} of {}'.format(mpi_comm.Get_rank(), mpi_comm.Get_size()))

pyshdom.util.get_measurements(solvers, sensors_dict, n_jobs=4)

#if mpi_comm.Get_rank() == 0:
pyshdom.util.save_forward_model(file_name, sensors_dict, solvers)
