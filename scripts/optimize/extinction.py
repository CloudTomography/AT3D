""" 
Optimize: Extinction
--------------------

Optimize for the extinction coefficient based on radiance measurements.
Measurements are either:
  1. Simulated measurements using a forward rendering script (e.g. in scripts/render/).
  2. Real radiance measurements

The phase function, albedo and rayleigh scattering are assumed known.

Example usage:
  python scripts/optimize/extinction.py --input_dir DIRECTORY --log LOGDIR --use_forward_grid --use_forward_albedo --use_forward_phase --init Homogeneous --extinction 0.0 --add_rayleigh --mie_table_path TABLE_PATH
  
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
        
    Returns
    -------
    measurments: shdom.Measurements object
        Encapsulates the radiances and the sensor geometry for later optimization
    rte_solver: shdom.RteSolver
        A solver initialized to the initial medium
    mask: shdom.GridData
        A mask of cloudy voxels, either by space-carving or from ground-truth
    extinction: shdom.Parameters.Extinction object
        The unknown parameters to recover
    air: shdom.Medium
        A known medium of air which overlays the atmosphere.
    """
    medium_gt, rte_solver, measurements = shdom.load_forward_model(args.input_dir)
    
    # Multi-band optimization currently not supported
    if measurements.sensors.type == 'SensorArray':
        wavelengths = map(lambda x: x.wavelength, measurements.sensors.sensor_list)
        assert np.all(np.isclose(np.diff(wavelengths), 0.0)), 'Multi-band optimization currently not supported.'
        args.wavelength = wavelengths[0]
    else:
        args.wavelength = measurements.sensors.wavelength
    
        
    cloud_generator = CloudGenerator(args)

    # Define the grid for reconstruction
    if args.use_forward_grid:
        cloud_generator.set_grid(medium_gt.grid)
    else:
        cloud_generator.init_grid()

    # Define the known albedo and phase (either ground-truth or specified, but it is not optimized)
    if args.use_forward_albedo:
        cloud_generator.set_albedo(medium_gt.albedo)
    else:
        cloud_generator.init_albedo()
        
    if args.use_forward_phase:
        cloud_generator.set_phase(medium_gt.phase)
    else:
        cloud_generator.init_phase()
        
    # Define an unknown extinction field to recover
    cloud_generator.init_extinction()
    extinction = shdom.Parameters.Extinction()
    extinction.initialize(cloud_generator.extinction, 
                          min_bound=0.0,
                          max_bound=None)
    
    # Initilize the Medium and RTE solver    
    medium = shdom.Medium()
    medium.set_optical_properties(extinction, 
                                  cloud_generator.albedo, 
                                  cloud_generator.phase)
    
    # Apply a cloud mask for non-cloudy voxels 
    if args.use_forward_mask:
        mask = medium_gt.get_mask(threshold=1.0)
    else:
        raise NotImplementedError 
    medium.apply_mask(mask)
    
    air = None
    if args.add_rayleigh:
        air_generator = AirGenerator(args)
        air = air_generator.get_medium(phase_type=medium.phase.type)   
        medium += air
    rte_solver.init_medium(medium)    
    
    return measurements, rte_solver, mask, extinction, air


def optimize(measurements, rte_solver, mask, extinction, air):
    """
    Optimize over the unknown extinction.
    
    Parameters
    ----------
    measurments: shdom.Measurements object
        Encapsulates the radiances and the sensor geometry for later optimization
    rte_solver: shdom.RteSolver
        A solver initialized to the initial medium
    mask: shdom.GridData
        A mask of cloudy voxels, either by space-carving or from ground-truth
    extinction: shdom.Parameters.Extinction object
        The unknown parameters to recover
    air: shdom.Medium
        A known medium of air which overlays the atmosphere
        
    Returns
    -------
    optimizer: shdom.Optimizer object
        The optimizer after the termination of the optimization process
    """
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
    optimizer.set_cloud_mask(mask)
    optimizer.add_parameter(extinction)
    optimizer.set_known_medium(air) 
    
    if args.log is not None:
        from tensorboardX import SummaryWriter
        timestr = time.strftime("%Y%m%d-%H%M%S")
        log_dir = os.path.join(args.input_dir, 'logs', args.log + timestr)
        writer = SummaryWriter(log_dir)      
        optimizer.set_writer(writer)    

    # Start optimization process
    result = optimizer.minimize(ckpt_period=1*60*60, options=options)
    print('\n------------------ Optimization Finished ------------------\n')
    print('Status: {}'.format(result.status))
    print('Message: {}'.format(result.message))
    print('Final loss: {}'.format(result.fun))
    print('Number iterations: {}'.format(result.nit))
    return optimizer

if __name__ == "__main__":
    
    args, CloudGenerator, AirGenerator = argument_parsing()
    measurements, rte_solver, mask, extinction, air = init_atmosphere(args, CloudGenerator, AirGenerator)
    optimizer = optimize(measurements, rte_solver, mask, extinction, air)
    optimizer.save(os.path.join(writer.log_dir, 'final.ckpt'))


