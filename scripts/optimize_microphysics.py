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
        Arguments required for this script.
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
    parser.add_argument('--use_forward_lwc',
                        action='store_true',
                        help='Use the ground-truth LWC.')
    parser.add_argument('--use_forward_reff',
                            action='store_true',
                            help='Use the ground-truth effective radius.')    
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
                        help='(default value: %(default)s) Threshold for the radiance to create a cloud mask.'
                        'Threshold is either a scalar or a list of length of measurements.')    
    parser.add_argument('--n_jobs',
                        default=1,
                        type=int,
                        help='(default value: %(default)s) Number of jobs for parallel rendering. n_jobs=1 uses no parallelization')
    parser.add_argument('--globalopt',
                        action='store_true',
                        help='Global optimization with basin-hopping.'
                             'For more info see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html')

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
        lwc_grid = reff_grid = veff_grid = cloud_generator.get_grid()
   
    # Define Microphysical parameters: either optimize or use ground-truth
    if args.use_forward_lwc:
        lwc = cloud_gt.lwc
    else:
        lwc = shdom.GridDataEstimator(cloud_generator.get_lwc(lwc_grid),
                                      min_bound=0.0,
                                      random_step=0.15)
    
    if args.use_forward_reff:
        reff = cloud_gt.reff
    else:
        reff = shdom.GridDataEstimator(cloud_generator.get_reff(reff_grid),
                                       min_bound=cloud_gt.min_reff,
                                       max_bound=cloud_gt.max_reff,
                                       random_step=5.0)
    if args.use_forward_veff:
        veff = cloud_gt.veff
    else:
        veff = shdom.GridDataEstimator(cloud_generator.get_veff(veff_grid), 
                                       max_bound=cloud_gt.max_veff,
                                       min_bound=cloud_gt.min_veff,
                                       random_step=0.1)
        
    cloud_estimator = shdom.MicrophysicalScattererEstimator(cloud_gt.mie, lwc, reff, veff)
    
    # Set a cloud mask for non-cloudy voxels
    if args.use_forward_mask:
        mask = cloud_gt.get_mask(threshold=0.01)
    else:
        carver = shdom.SpaceCarver(measurements)
        mask = carver.carve(cloud_estimator.grid, agreement=0.95, thresholds=args.radiance_threshold)
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
    medium_estimator.add_scatterer(cloud_estimator, name='cloud')
    
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
    log_dir = args.input_dir
    if args.log is not None:
        log_dir = os.path.join(args.input_dir, 'logs', args.log + '-' + time.strftime("%d-%b-%Y-%H:%M:%S"))
        writer = shdom.SummaryWriter(log_dir)
        writer.monitor_loss()
        writer.save_checkpoints(ckpt_period=30*60)
        writer.monitor_images(acquired_images=measurements.images, ckpt_period=3*60)
        writer.monitor_scatterer_error(estimator_name='cloud', ground_truth=cloud_gt)
        writer.monitor_domain_mean(estimator_name='cloud', ground_truth=cloud_gt)

    # Define an L-BFGS-B local optimizer
    options = {
        'maxiter': 1000,
        'maxls': 100,
        'disp': True,
        'gtol': 1e-18,
        'ftol': 1e-18
    }
    optimizer = shdom.LocalOptimizer(options, n_jobs=args.n_jobs)
    optimizer.set_measurements(measurements)
    optimizer.set_rte_solver(rte_solver)
    optimizer.set_medium_estimator(medium_estimator)
    optimizer.set_writer(writer)

    n_global_iter = 1
    if args.globalopt:
        global_optimizer = shdom.GlobalOptimizer(local_optimizer=optimizer)
        result = global_optimizer.minimize(niter_success=20, T=1e-3)
        n_global_iter = result.nit
        result = result.lowest_optimization_result
        optimizer.set_state(result.x)
    else:
        result = optimizer.minimize()
    optimizer.save(os.path.join(log_dir, 'optimizer'))

    print('\n------------------ Optimization Finished ------------------\n')
    print('Number global iterations: {}'.format(n_global_iter))
    print('Success: {}'.format(result.success))
    print('Message: {}'.format(result.message))
    print('Final loss: {}'.format(result.fun))
    print('Number local iterations: {}'.format(result.nit))
    print(result)



