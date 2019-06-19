""" 
Optimize: Microphysics
----------------------
 
Optimize for the microphysics based on radiance measurements.
Measurements are either:
  1. Simulated measurements using a forward rendering script (e.g. in scripts/render_radiance_toa.py).
  2. Real radiance measurements

For example usage see the README.md

For information about the command line flags see:
  python scripts/optimize_microphysics.py --help
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
    parser.add_argument('--use_forward_veff',
                        action='store_true',
                        help='Use the ground-truth effective variance.')
    parser.add_argument('--use_forward_mask',
                        action='store_true',
                        help='Use the ground-truth cloud mask. This is an inverse crime which is \
                              usefull for debugging/development.')
    parser.add_argument('--radiance_threshold',
                        default=[0.05],
                        nargs='+',
                        type=np.float32,
                        help='(default value: %(default)s) Threshold for the radiance to create a cloud mask.' \
                        'Threshold is either a scalar or a list of length of measurements.')    
    parser.add_argument('--n_jobs',
                        default=1,
                        type=int,
                        help='(default value: %(default)s) Number of jobs for parallel rendering. n_jobs=1 uses no parallelization')
    
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


def init_atmosphere_estimation():
    """
    Initilize the atmosphere for optimization.
    """
    cloud_generator = CloudGenerator(args)
    
    # Define the grid for reconstruction
    if args.use_forward_grid:
        lwc_grid = cloud_gt.lwc.grid
        reff_grid = cloud_gt.reff.grid
        veff_grid = cloud_gt.reff.grid
    else:
        lwc_grid = reff_grid = veff_grid = cloud_generator.grid
   
    # Define the Microphysical parameters
    # Effective variance is either initialized (optimized) or ground-truth (not optimized)    
    lwc = shdom.GridDataEstimator(cloud_generator.get_lwc(lwc_grid), 
                                  min_bound=0.0)
    reff = shdom.GridDataEstimator(cloud_generator.get_reff(reff_grid),
                                   max_bound=cloud_gt.max_reff, 
                                   min_bound=cloud_gt.min_reff)
    if args.use_forward_veff:
        veff = cloud_gt.veff
    else:
        veff = shdom.GridDataEstimator(cloud_generator.get_veff(veff_grid), 
                                       max_bound=cloud_gt.max_veff,
                                       min_bound=cloud_gt.min_veff)
        
    cloud_estimator = shdom.MicrophysicalScattererEstimator(lwc, reff, veff)
    cloud_estimator.add_mie(cloud_gt.mie)
    
    # Set a cloud mask for non-cloudy voxels
    if args.use_forward_mask:
        mask = cloud_gt.get_mask(threshold=0.001)
    else:
        carver = shdom.SpaceCarver(measurements)
        mask = carver.carve(grid, agreement=0.95, thresholds=args.radiance_threshold)
    cloud_estimator.set_mask(mask)

    # Create a medium estimator object (optional Rayleigh scattering)
    medium_estimator = shdom.MediumEstimator()
    if args.add_rayleigh:
        air_generator = AirGenerator(args)
        air = air_generator.get_scatterer(cloud_estimator.wavelength)   
        medium_estimator.set_grid(cloud_estimator.grid + air.grid)
        medium_estimator.add_scatterer(air, 'air') 
    else:
        medium_estimator.set_grid(cloud_estimator.grid)
    medium_estimator.add_scatterer(cloud_estimator, name='cloud estimator')
    
    return medium_estimator


if __name__ == "__main__":
    
    args, CloudGenerator, AirGenerator = argument_parsing()
    
    # Load forward model
    medium_gt, rte_solver, measurements = shdom.load_forward_model(args.input_dir)
    
    # Get micro-physical medium ground-truth and mie tables
    cloud_gt = medium_gt.get_scatterer('cloud')
    
    # Init medium estimator
    medium_estimator = init_atmosphere_estimation()
    
    # Define a summary writer
    writer = None
    if args.log is not None:
        log_dir = os.path.join(args.input_dir, 'logs', args.log + '-' + time.strftime("%Y%m%d-%H%M%S"))
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
        'gtol': 1e-18,
        'ftol': 1e-18 
    }
    
    optimizer.set_measurements(measurements)
    optimizer.set_rte_solver(rte_solver)
    optimizer.set_medium_estimator(medium_estimator)
    optimizer.set_writer(writer)

    # Optimization process
    result = optimizer.minimize(ckpt_period=30*60, options=options, n_jobs=args.n_jobs)
    print('\n------------------ Optimization Finished ------------------\n')
    print('Success: {}'.format(result.success))
    print('Message: {}'.format(result.message))
    print('Final loss: {}'.format(result.fun))
    print('Number iterations: {}'.format(result.nit))    

    optimizer.save(os.path.join(args.input_dir, 'optimizer'))


