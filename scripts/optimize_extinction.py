""" 
Optimize: Extinction
--------------------

Optimize for the extinction coefficient based on radiance measurements.
Measurements are either:
  1. Simulated measurements using a forward rendering script (e.g. in scripts/render/).
  2. Real radiance measurements

The phase function, albedo and rayleigh scattering are assumed known.

Example usage:
  python scripts/optimize_extinction.py --input_dir DIRECTORY --mie_table_path TABLE_PATH \
         --use_forward_grid --use_forward_albedo --use_forward_phase \
         --init Homogeneous --extinction 0.0 --add_rayleigh --wavelength 0.672 --log test
  
For information about the command line flags see:
  python scripts/optimize/extinction.py --help
"""

import os, time
import numpy as np
import argparse
import shdom


def argument_parsing():
    """
    Handle all the argument parsing needed for this script.
    
    Returns
    -------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this script.
    CloudGenerator: a shdom.Generator class object.
        Creates the cloudy medium.
    AirGenerator: a shdom.Air class object
        Creates the scattering due to air molecules
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_dir', 
                        help='Path to an input directory where the forward modeling parameters are be saved. \
                              This directory will be used to save the optimization results and progress.')
    parser.add_argument('--log',
                        help='Write intermediate TensorBoardX results. \
                              The provided string is added as a comment to the specific run.')
    parser.add_argument('--use_forward_grid',
                        action='store_true',
                        help='Use the same grid for the reconstruction. This is a sort of inverse crime which is \
                              usefull for debugging/development.')
    parser.add_argument('--use_forward_albedo',
                        action='store_true',
                        help='Use the ground truth albedo.')
    parser.add_argument('--use_forward_phase',
                        action='store_true',
                        help='Use the ground-truth phase reconstruction.')
    parser.add_argument('--use_forward_mask',
                        action='store_true',
                        help='Use the ground-truth cloud mask. This is an inverse crime which is \
                              usefull for debugging/development.')
    
    # Additional arguments to the parser
    subparser = argparse.ArgumentParser(add_help=False)
    subparser.add_argument('--init')
    subparser.add_argument('--add_rayleigh', action='store_true')
    parser.add_argument('--init',
                        default='Homogeneous',
                        help='(default value: %(default)s) Name of the generator used to initialize the atmosphere. \
                              for additional generator arguments: python scripts/optimize/extinction.py --generator=GENERATOR -h. \
                              See scripts/generate.py for more documentation.')
    parser.add_argument('--add_rayleigh',
                        action='store_true',
                        help='Overlay the atmosphere with (known) Rayleigh scattering due to air molecules. \
                              Temperature profile is taken from AFGL measurements of summer mid-lat.')
    
    init = subparser.parse_known_args()[0].init 
    add_rayleigh = subparser.parse_known_args()[0].add_rayleigh 
    
    if init:
        CloudGenerator = getattr(shdom.generate, init)
        parser = CloudGenerator.update_parser(parser)
     
    AirGenerator = None  
    if add_rayleigh:
        AirGenerator = shdom.generate.AFGLSummerMidLatAir
        parser = AirGenerator.update_parser(parser)
    
    args = parser.parse_args()
    
    return args, CloudGenerator, AirGenerator


def init_atmosphere(args, CloudGenerator, AirGenerator):
    """
    Load forward modeling setup parameters.
    
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this script.
    CloudGenerator: a shdom.Generator class object.
        Creates the cloudy medium.
    AirGenerator: a shdom.Air class object
        Creates the scattering due to air molecules
    measurements: shdom.Measurements object
        Used for space-carving a mask if args.use_forward_mask is False
        
    Returns
    -------
    medium: shdom.OpticalMediumEstimator
        An OpticalMediumEstimator object.
    """
    
    cloud_generator = CloudGenerator(args)

    # Define the grid for reconstruction
    if args.use_forward_grid:
        cloud_generator.set_grid(medium_gt.grid)
    else:
        cloud_generator.init_grid()

    cloud_generator.init_extinction()
    
    # Define the known albedo and phase (either ground-truth or specified, but it is not optimized)
    if args.use_forward_albedo:
        cloud_generator.set_albedo(medium_gt.get_scatterer('cloud').albedo)
    else:
        cloud_generator.init_albedo()
        
    if args.use_forward_phase:
        cloud_generator.set_phase(medium_gt.get_scatterer('cloud').phase)
    else:
        cloud_generator.init_phase()
           
    cloud = shdom.ScattererEstimator(
        extinction=shdom.GridDataEstimator(cloud_generator.extinction, min_bound=0.0),
        albedo=cloud_generator.albedo,
        phase=cloud_generator.phase
    )
    
    # Set a cloud mask for non-cloudy voxels
    if args.use_forward_mask:
        mask = medium_gt.get_scatterer('cloud').get_mask(threshold=1.0)
    else:
        carver = shdom.SpaceCarver(measurements)
        mask = carver.carve(cloud.grid, agreement=0.7, thresholds=0.05)
    cloud.set_mask(mask)
    
    medium = shdom.OpticalMediumEstimator(cloud.grid)
    medium.add_scatterer(cloud, name='cloud estimator')
    
    if args.add_rayleigh:
        air_generator = AirGenerator(args)
        medium.add_scatterer(air_generator.get_scatterer(), 'air')        

    return medium


if __name__ == "__main__":
    
    args, CloudGenerator, AirGenerator = argument_parsing()
    
    # Load forward model
    medium_gt, rte_solver, measurements = shdom.load_forward_model(args.input_dir)
    
    # Init medium and solver
    medium = init_atmosphere(args, CloudGenerator, AirGenerator)
    
    # Define a summary writer
    writer = None
    if args.log is not None:
        log_dir = os.path.join(args.input_dir, 'logs', args.log + time.strftime("%Y%m%d-%H%M%S"))
        writer = shdom.SummaryWriter(log_dir)
        writer.monitor_loss()
        writer.monitor_images(acquired_images=measurements.images)
        writer.monitor_parameter_error(ground_truth=medium_gt.get_scatterer('cloud').extinction)
        
    optimizer = shdom.Optimizer()
        
    # Define L-BFGS-B options
    options = {
        'maxiter': 1000,
        'maxls': 100,
        'disp': True,
        'gtol': 1e-16,
        'ftol': 1e-16 
    }
    
    optimizer.set_measurements(measurements)
    optimizer.set_rte_solver(rte_solver)
    optimizer.set_medium_estimator(medium)
    optimizer.set_writer(writer)

    # Optimization process
    result = optimizer.minimize(ckpt_period=30*60, options=options)
    print('\n------------------ Optimization Finished ------------------\n')
    print(result)
    print('Success: {}'.format(result.success))
    print('Message: {}'.format(result.message))
    print('Final loss: {}'.format(result.fun))
    print('Number iterations: {}'.format(result.nit))    

    optimizer.save(os.path.join(args.input_dir, 'optimizer'))


