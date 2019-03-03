""" 
Inverse: optimize extinction (simulation)
-----------------------------------------
Optimize for the extinction coefficient with a homogeneous field initialization. 
Measurements should be generated using a forward script (e.g. scripts/forward_single_voxel).
The phase function, albedo and rayleigh scattering are assumed known.

Example usage:
  TODO
  
For information about the command line flags see:
  python render.py --help

For a tutorial overview of how to operate the forward rendering model see the following notebooks:
 - notebooks/Make Mie Table.ipynb
 - notebooks/Forward Rendering.ipynb
"""

import os 
import numpy as np
import argparse
import shdom

parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', 
                    help='Path to an input directory where the forward modeling parameters are be saved. \
                          This directory will be used to save the optimization results and progress.')

parser.add_argument('--initial_value', 
                    default=1e-6,
                    type=np.float32, 
                    help='(default value: %(default)s) Initial value for a homogeneous extinction field [km^-1]. \
                          Note: initial_extinction=0.0 is an unstable value.')

parser.add_argument('--use_forward_grid',
                    action='store_true',
                    help='Use the same grid for the reconstruction. This is a sort of inverse crime which is \
                          usefull for debugging/development.')

args = parser.parse_args()


def load_forward_model(directory):
    """
    Save the forward model parameters for reconstruction.
    
    Parameters
    ----------
    directory: str
        Directory path where the forward modeling parameters are saved. 
        If the folder doesnt exist it will be created.
    
    Notes
    -----
    medium: shdom.Medium object
        The atmospheric medium. This ground-truth medium will be used to 
    sensor: shdom.Sensor object
        The sensor used to image the medium.
    solver: shdom.RteSolver object
        The solver and the parameters used. This includes the scene parameters (such as solar and surface parameters)
        and the numerical parameters.
    
    Notes
    -----
    The ground-truth atmosphere is used for evaulation of the recovery.
    """  
    
    # Load the ground truth atmosphere for error analysis and ground-truth known phase and albedo
    medium = shdom.Medium()
    medium.load(path=os.path.join(directory, 'ground_truth_medium'))   
    
    # Load Sensor according to it's parameters
    sensor = shdom.Sensor()
    sensor.load(path=os.path.join(directory, 'sensor'))
    
    # Load RteSolver according to numerical and scene parameters
    solver = shdom.RteSolver()
    solver.load_params(path=os.path.join(directory, 'solver_params'))   

    # Load measurements for reconstruction
    measurements = np.load(os.path.join(directory, 'measurements.npy')) 
    
    return medium, sensor, solver, measurements
    
if __name__ == "__main__":
    
    medium_gt, sensor, rte_solver, measurements = load_forward_model(directory)
    
    # Define the grid for recovery
    if args.use_forward_grid:
        grid = ground_truth_atmosphere.grid
    else:
        raise NotImplementedError #TODO
        
    # Define the optical fields: extinction, albedo, phase
    extinction_init = shdom.GridData(grid, np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=args.initial_value))
    medium_init = shdom.Medium()
    medium_init.set_optical_properties(extinction_init, medium_gt.albedo, medium_gt.phase)    

    # Setup an optimizer object
    optimizer = shdom.Optimizer()
    optimizer.set_measurements(measurements)
    optimizer.set_sensor(sensor)
    optimizer.set_rte_solver(rte_solver)
    
    # Define L-BFGS-B options
    options = {
        'maxiter': 100,
        'maxls': 100,
        'disp': True,
        'gtol': 1e-12,
        'ftol': 1e-12 
    }
    
    # Start optimization 
    medium_rec, optimizer_result = optimizer.minimize(init=medium_init, 
                                                      options=options, 
                                                      method='L-BFGS-B')
    
    
    


