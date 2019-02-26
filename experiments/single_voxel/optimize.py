import numpy as np
import shdom
import os


# Load the ground truth atmosphere for error analysis
ground_truth_atmosphere = shdom.Medium()
ground_truth_atmosphere.load('ground_truth_medium')

# Load the ground truth images for reconstruction
measurements = np.load('measurements.npy')

# Inverse crime: use the same grid for reconstruction
grid = ground_truth_atmosphere.grid

# Inverse crime: use the ground-truth phase/albedo as a prescribed phase function
phase_gt = ground_truth_atmosphere.phase
albedo_gt = ground_truth_atmosphere.albedo

# Define the optical fields: extinction, albedo, phase
ext_init = shdom.GridData(grid, np.full(shape=(3, 3, 3), fill_value=1e-5, dtype=np.float32))
atmosphere = shdom.Medium()
atmosphere.set_optical_properties(ext_init, albedo_gt, phase_gt)

# Load RteSolver according to parameters and scene parameters
rte_solver = shdom.RteSolver()
rte_solver.load_params('solver_params')
rte_solver.init_medium(atmosphere)

# Load Sensor according to it's parameters
sensor = shdom.Sensor()
sensor.load('sensor')

# Setup an optimizer object
optimizer = shdom.Optimizer()
optimizer.set_measurements(measurements)
optimizer.set_sensor(sensor)
optimizer.set_rte_solver(rte_solver)

result = optimizer.minimize()
print(result)


pass


