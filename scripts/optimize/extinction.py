""" 
Optimize: Extinction
--------------------

Optimize for the extinction coefficient based on radiance measurements.
Measurements are either:
  1. Simulated measurements using a forward rendering script (e.g. in scripts/render).
  2. Real radiance measurements (e.g. AirMSPI multi-view)

The phase function, albedo and rayleigh scattering are assumed known.

Example usage:
  TODO
  
For information about the command line flags see:
  python scripts/optimize/extinction.py --help
"""

import os 
import numpy as np
import argparse
import shdom
from shdom import parameters
import scripts.generate.air as generate_air

parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', 
                    help='Path to an input directory where the forward modeling parameters are be saved. \
                          This directory will be used to save the optimization results and progress.')

parser.add_argument('--initial_value', 
                    default=0.1,
                    type=np.float32, 
                    help='(default value: %(default)s) Initial value for a homogeneous extinction field [km^-1]. \
                          Note: initial_extinction=0.0 is an unstable value.')

parser.add_argument('--use_forward_grid',
                    action='store_true',
                    help='Use the same grid for the reconstruction. This is a sort of inverse crime which is \
                          usefull for debugging/development.')

parser.add_argument('--gt_cloud_mask',
                    action='store_true',
                    help='Use the ground-truth cloud masks for reconstruction. This is an inverse crime which is \
                          usefull for debugging/development.')

# Additional arguments to the parser
subparser = argparse.ArgumentParser(add_help=False)
subparser.add_argument('--add_rayleigh', action='store_true')
parser.add_argument('--add_rayleigh',
                    action='store_true',
                    help='Overlay the atmosphere with (known) Rayleigh scattering due to air molecules. \
                          Temperature profile is taken from AFGL measurements of summer mid-lat.')

add_rayleigh = subparser.parse_known_args()[0].add_rayleigh 

if add_rayleigh:
    parser = generate_air.update_parser(parser)

args = parser.parse_args()


if __name__ == "__main__":
    
    medium_gt, rte_solver, measurements = shdom.load_forward_model(args.input_dir)
    
    # Multi-band optimization is not supported
    if measurements.sensors.type == 'SensorArray':
        wavelengths = map(lambda x: x.wavelength, measurements.sensors.sensor_list)
        assert np.all(np.isclose(np.diff(wavelengths), 0.0)), 'Multi-band optimization is not supported. S'
        args.wavelength = wavelengths[0]
    else:
        args.wavelength = measurements.sensors.wavelength

    # Define medium grid for recovery
    if args.use_forward_grid:
        grid = medium_gt.grid
    else:
        raise NotImplementedError #TODO
    
    
    # Define cloud mask for recovery
    if args.gt_cloud_mask:
        cloud_mask = medium_gt.get_cloud_mask()
    else:
        raise NotImplementedError #TODO
    
    
    # Define the unknown extinction field 
    extinction = parameters.Extinction()
    ext_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz))
    ext_data[cloud_mask.data] = args.initial_value
    extinction_init = shdom.GridData(grid, ext_data)
    extinction.initialize(extinction_init, min_bound=0.0)
    
    
    # Initilize the RTE solver to the medium to ground-truth albedo and phase
    medium_init = shdom.Medium()
    medium_init.set_optical_properties(extinction, medium_gt.albedo, medium_gt.phase) 

    if args.add_rayleigh:
        air = generate_air.summer_mid_lat(args)
        medium_init = medium_init + air
    rte_solver.init_medium(medium_init)
    
    
    # Setup an optimizer object
    optimizer = shdom.Optimizer()
    optimizer.set_measurements(measurements)
    optimizer.set_rte_solver(rte_solver)
    optimizer.set_cloud_mask(cloud_mask)
    optimizer.add_parameter(extinction)
    
    if add_rayleigh:
        optimizer.set_known_medium(air)
    
    
    # Define L-BFGS-B options
    options = {
        'maxiter': 100,
        'maxls': 100,
        'disp': True,
        'gtol': 1e-12,
        'ftol': 1e-12 
    }
    
    
    # Start optimization 
    optimizer_result = optimizer.minimize(options=options)
    print(optimizer_result)
    
    


