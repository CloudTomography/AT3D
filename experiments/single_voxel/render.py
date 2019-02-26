""" 
Forward Rendering: Single Voxel
-------------------------------

Define an atmospheric volume of 3x3x3 voxels with a single cloudy voxel (Mie scattering) in the center. 
The cloudy voxel is surrounded by air (Rayleigh scattering).

For a tutorial overview of how to operate the forward rendering model see the following notebooks:
 - notebooks/Make Mie Table.ipynb
 - notebooks/Forward Rendering.ipynb
"""

import os 
import matplotlib.pyplot as plt
import numpy as np
import shdom

# Mie scattering for water droplets
mie = shdom.Mie()
mie.read_table(file_path='../../mie_tables/Water_672nm.scat')

# Rayleigh scattering for air molecules
temperatures = np.array([292.220, 292.040, 291.860, 291.680, 291.500, 291.320, 291.140, 290.960, 290.780, 
                         290.600, 290.420, 290.240, 290.060, 289.880, 289.700, 289.920, 290.140, 290.360, 
                         290.580, 290.800, 291.020, 291.240, 291.460, 291.680, 291.900])
temp_grid = shdom.Grid(z=np.linspace(0.0, 20.0, len(temperatures)))
temperature_profile = shdom.GridData(temp_grid, temperatures)
rayleigh = shdom.Rayleigh(wavelength=0.672, temperature_profile=temperature_profile)


# Extinction of the center cloudy voxel elemet 
cloud_extinction = 1.34

# Define the atmospheric grid
bounding_box = shdom.BoundingBox(-0.5, -0.5, 0.0, 0.5, 0.5, 1.0)
grid = shdom.Grid(bounding_box=bounding_box, nx=3, ny=3, nz=3)

# Define the optical fields: extinction, albedo, phase
ext_data = np.zeros(shape=(3, 3, 3), dtype=np.float32)
alb_data = np.zeros(shape=(3, 3, 3), dtype=np.float32)
reff_data = np.zeros(shape=(3, 3, 3), dtype=np.float32)

# Define center voxel parameters
ext_data[1, 1, 1] = cloud_extinction
alb_data[1, 1, 1] = 1.0
reff_data[1, 1, 1] = 10.0

# Define a Medium object and inilize it's optical properties
extinction_c = shdom.GridData(grid, ext_data)
albedo_c = shdom.GridData(grid, alb_data)
phase_c = mie.get_grid_phase(shdom.GridData(grid, reff_data))

atmosphere = shdom.Medium()
atmosphere.set_optical_properties(extinction_c, albedo_c, phase_c)

# Solve the Radiative Transfer for the domain
rte_solver = shdom.RteSolver()
scene_params = shdom.SceneParameters()
numerical_params = shdom.NumericalParameters()
scene_params.source.azimuth = 35
scene_params.source.zenith = 165
rte_solver.set_params(scene_params, numerical_params)
rte_solver.init_medium(atmosphere)
rte_solver.solve(maxiter=100)

# Define a sensor and render an image of the domain. 
sensor = shdom.OrthographicSensor(bounding_box=atmosphere.bounding_box, 
                                  x_resolution=0.1, 
                                  y_resolution=0.1, 
                                  azimuth=0.0, 
                                  zenith=0.0,
                                  altitude='TOA')
image = sensor.render(rte_solver)

# Save measurements, groud-truth medium and rte solver parameters
np.save('measurements.npy', image)
atmosphere.save('ground_truth_medium')
sensor.save('sensor')
rte_solver.save_params('solver_params')


# extinction_a, albedo_a, phase_a = rayleigh.get_scattering_field(grid.z_grid, phase_type='Grid')
# air.set_optical_properties(extinction_a, albedo_a, phase_a)
#scene_params.source.zenith = 160
#scene_params.source.azimuth = 30
#scene_params.source.flux = 1.0