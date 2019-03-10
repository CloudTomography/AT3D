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
import importlib
import scripts.generate.air as generate_air

parser = argparse.ArgumentParser()

parser.add_argument('--input_dir', 
                    help='Path to an input directory where the forward modeling parameters are be saved. \
                          This directory will be used to save the optimization results and progress.')

parser.add_argument('--summary',
                    help='Write intermediate TensorBoardX results. \
                          The provided string is added as a comment to the specific run.')

parser.add_argument('--use_forward_grid',
                    action='store_true',
                    help='Use the same grid for the reconstruction. This is a sort of inverse crime which is \
                          usefull for debugging/development.')

parser.add_argument('--use_forward_albedo',
                    action='store_true',
                    help='Use the same grid for the reconstruction. This is a sort of inverse crime which is \
                          usefull for debugging/development.')

parser.add_argument('--use_forward_phase',
                    action='store_true',
                    help='Use the same grid for the reconstruction. This is a sort of inverse crime which is \
                          usefull for debugging/development.')

parser.add_argument('--cloud_mask',
                    choices=['spacecarve', 'groundtruth'],
                    help='(default value: %(default)s) spacecarve - find the mask using space carving (geometric considirations). \
                                                       groundtruth - use an extinction based mask cfrom the forward model. \
                                                       This is an inverse crime which is usefull for debugging/development.')

# Additional arguments to the parser
subparser = argparse.ArgumentParser(add_help=False)
subparser.add_argument('--init')
subparser.add_argument('--add_rayleigh', action='store_true')
parser.add_argument('--init', 
                    help='Name of the generator used to initialize the atmosphere. \
                          or additional generator arguments: python scripts/optimize/extinction.py --generator=GENERATOR -h. \
                          A generator should have the following methods to generate the medium: \
                          update_parser, grid, extinction, albedo, phase. \
                          See scripts/generate/__init__.py for more documentation.')

parser.add_argument('--add_rayleigh',
                    action='store_true',
                    help='Overlay the atmosphere with (known) Rayleigh scattering due to air molecules. \
                          Temperature profile is taken from AFGL measurements of summer mid-lat.')

init = subparser.parse_known_args()[0].init 
add_rayleigh = subparser.parse_known_args()[0].add_rayleigh 
if init:
    init = importlib.import_module('scripts.generate.' + init)
    parser = init.update_parser(parser)
    
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

    # Define a Grid
    if args.use_forward_grid:
        grid = medium_gt.grid
    else:
        grid = init.grid(args)
    
    
    # Define a cloud mask 
    if args.cloud_mask is None:
        cloud_mask = None
    elif args.cloud_mask == 'groundtruth':
        cloud_mask = medium_gt.get_cloud_mask()
    elif args.cloud_mask == 'spacecarve':
        raise NotImplementedError
    
    
    # Define the unknown extinction field 
    extinction = parameters.Extinction()
    extinction_init = init.extinction(grid, args, cloud_mask) 
    extinction.initialize(extinction_init, min_bound=0.0)
    
    
    # Define the known albedo field (cloud be ground-truth or specified, but it is not optimized)
    if args.use_forward_albedo:
        albedo = medium_gt.albedo
    else:
        albedo = init.albedo(grid, args, cloud_mask)
    
    
    # Define the known phase field (cloud be ground-truth or specified, but it is not optimized)
    if args.use_forward_phase:
        phase = medium_gt.phase
    else:
        phase = init.phase(grid, args, cloud_mask)   
        
    
    # Initilize the Medium and RTE solver    
    medium = shdom.Medium()
    medium.set_optical_properties(extinction, albedo, phase) 
    if args.add_rayleigh:
        air = generate_air.summer_mid_lat(args)
        medium += air
    rte_solver.init_medium(medium)
    
    
    # Setup an optimizer object
    optimizer = shdom.Optimizer()
    optimizer.set_measurements(measurements)
    optimizer.set_rte_solver(rte_solver)
    optimizer.set_cloud_mask(cloud_mask)
    optimizer.add_parameter(extinction)
    
    if add_rayleigh:
        optimizer.set_known_medium(air)
    
    if args.summary is not None:
        from tensorboardX import SummaryWriter
        log_dir = os.path.join(args.input_dir, 'log')
        writer = SummaryWriter(log_dir, comment=args.summary)      
        optimizer.set_writer(writer)
        
        
    # Define L-BFGS-B options
    options = {
        'maxiter': 100,
        'maxls': 100,
        'disp': False,
        'gtol': 1e-10,
        'ftol': 1e-10 
    }
    
    
    # Start optimization
    optimizer_result = optimizer.minimize(options=options)
    print(optimizer_result)
    
    


