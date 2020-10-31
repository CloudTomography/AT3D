#imports
import shdom
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import OrderedDict
import sys
import shutil

file_name = sys.argv[-1]

os.chdir('/Users/jesserl2/Documents/Code/aviad_pyshdom_dev/pyshdom_dev/')

#load the cloud.
cloud_scatterer = shdom.grid.load_from_csv('./synthetic_cloud_fields/jpl_les/rico32x37x26.txt',
                                           density='lwc',origin=(0.0,0.0))
#load atmosphere
atmosphere = xr.open_dataset('./ancillary_data/AFGL_summer_mid_lat.nc')
#subset the atmosphere, choose only the bottom kilometre.
reduced_atmosphere = atmosphere.sel({'z': atmosphere.coords['z'].data[atmosphere.coords['z'].data <= 1.0]})
#merge the atmosphere and cloud z coordinates
merged_z_coordinate = shdom.grid.combine_z_coordinates([reduced_atmosphere,cloud_scatterer])

#make a merged grid for the rte.
rte_grid = shdom.grid.make_grid(cloud_scatterer.x.data.min(),cloud_scatterer.x.data.max(),cloud_scatterer.x.data.size,
                           cloud_scatterer.y.data.min(),cloud_scatterer.y.data.max(),cloud_scatterer.y.data.size,
                           merged_z_coordinate)

#resample the cloud onto the rte_grid
cloud_scatterer_on_rte_grid = shdom.grid.resample_onto_grid(rte_grid, cloud_scatterer)

#define any necessary variables for microphysics here.
size_distribution_function = shdom.size_distribution.gamma
#We choose a gamma size distribution and therefore need to define a 'veff' variable.
cloud_scatterer_on_rte_grid['veff'] = (cloud_scatterer_on_rte_grid.reff.dims,
                                       np.full_like(cloud_scatterer_on_rte_grid.reff.data, fill_value=0.1))

#define sensors.
Sensordict = OrderedDict()
misr_list = []
sensor_zenith_list = [75.0,60.0,45.6,26.1]*2 + [0.0]
sensor_azimuth_list = [90]*4 + [-90]*4 +[0.0]
for zenith,azimuth in zip(sensor_zenith_list,sensor_azimuth_list):
    misr_list.append(
                    shdom.sensor.orthographic_projection(0.86, cloud_scatterer,0.02,0.02, azimuth, zenith,
                                             altitude='TOA', stokes='I')#, sub_pixel_ray_args={'method': shdom.sensor.gaussian,
                                                                    #                        'degree': (3,2)} #non-square pixel so add more
                                                                                                            #subpixel rays in larger direction.
                                            #)
    )
Sensordict['MISR'] = {'sensor_list': misr_list}

#Define the RTE solvers needed to model the measurements and
#calculate optical properties.
wavelengths = shdom.script_util.get_unique_wavelengths(Sensordict)

names = OrderedDict()
surfaces = OrderedDict()
sources = OrderedDict()
numerical_parameters = OrderedDict()
num_stokes = OrderedDict()
cloud_poly_tables = OrderedDict()
cloud_optical_scatterers = OrderedDict()

config = shdom.configuration.get_config('./default_config.json')
config['spherical_harmonics_accuracy'] = 0.01

for wavelength in wavelengths:
    num_stokes[wavelength] = 1
    names[wavelength] = None
    surfaces[wavelength] = shdom.surface.lambertian(albedo=0.0) #surface is wavelength independent.
    sources[wavelength] = shdom.source.solar(-1*np.cos(np.deg2rad(10.0)),0.0,solarflux=1.0)
    numerical_parameters[wavelength] = config

    #optical properties from mie calculations.
    mie_mono_table = shdom.mie.get_mono_table('Water',(wavelength,wavelength), relative_path='./mie_tables')
    cloud_size_distribution = shdom.size_distribution.get_size_distribution_grid(
                                                            mie_mono_table.radius.data,
                        size_distribution_function=size_distribution_function,particle_density=1.0,
                        reff=[4.0,25.0,25,'logarithmic','micron'],
                        veff=[0.09,0.11,2,'linear','unitless'],
                        )
    poly_table = shdom.mie.get_poly_table(cloud_size_distribution,mie_mono_table)
    cloud_optical_scatterers[wavelength] = shdom.medium.table_to_grid(cloud_scatterer_on_rte_grid, poly_table)
    cloud_poly_tables[wavelength] = poly_table

#calculate rayleigh properties.
rayleigh_scatterer_list = shdom.rayleigh.to_grid(wavelengths,atmosphere,rte_grid)

#just regroup optical properties to be grouped by wavelength for solver.

#SCATTERERS CAN (and should) BE NAMED.
mediums = shdom.script_util.combine_to_medium({'cloud':cloud_optical_scatterers, 'rayleigh':rayleigh_scatterer_list})

#define all solvers.
solvers = OrderedDict()
for key,name in names.items():
    solvers[key] = shdom.solver.RTE(numerical_params=numerical_parameters[key],
                                    medium=mediums[key],
                                   source=sources[key],
                                   surface=surfaces[key],
                                    num_stokes=num_stokes[key],
                                    name=name
                                   )

#solve RTE and calculate measurements.
shdom.script_util.get_measurements(solvers, Sensordict, maxiter=100, n_jobs=8, verbose=True)
shutil.copy(__file__, 'copied_script.py')
#save measurements, and also document all inputs.
shdom.util.save_forward_model(file_name, Sensordict, solvers)
