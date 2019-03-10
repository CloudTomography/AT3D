""" 
Render: Single wavelength radiance at Top of the Atmosphere (TOA)
-----------------------------------------------------------------

Forward rendering of an atmospheric medium with an orthographic sensor measuring exitting radiance at the top of the the domain.
This sensor is an (somewhat crude) approximation for far observing satellites where the rays are parallel. 

As with all `render` scripts a `generator` needs to be specified using the generator flag (--generator). 
The Generator defines the medium parameters: Grid, Extinction, Single Scattering Albedo and Phase function with it's own set of command-line flags.

Example usage:
  python scripts/render/radiance_toa.py --output_dir experiments/single_voxel/ --wavelength 0.672 \
          --generator single_voxel --extinction 10.0 --reff 10.0 --domain_size 1.0 --x_res 0.1 --y_res=0.1 \
          --mie_table_path mie_tables/Water_672nm.scat 

For information about the command line flags see:
  python scripts/render/render_radiance_toa.py --help
  
For a tutorial overview of how to operate the forward rendering model see the following notebooks:
 - notebooks/Make Mie Table.ipynb
 - notebooks/Forward Rendering.ipynb
"""

import os 
import numpy as np
import argparse
import shdom
import importlib
import scripts.generate.air as generate_air


parser = argparse.ArgumentParser()

parser.add_argument('output_dir', 
                    help='Path to an output directory where the measureents and model parameters will be saved. \
                          If the folder doesnt exist, it will be created.')

parser.add_argument('--solar_zenith', 
                    default=180.0,
                    type=np.float32, 
                    help='(default value: %(default)s) Solar zenith [deg]. This is the direction of the photons in range (90, 180]')

parser.add_argument('--solar_azimuth', 
                    default=0.0,
                    type=np.float32, 
                    help='(default value: %(default)s) Solar azimuth [deg]. This is the direction of the photons')

parser.add_argument('--wavelength', 
                    default=0.672,
                    type=np.float32, 
                    help='(default value: %(default)s) Wavelength [micron]. It is used to compute a Mie table if one \
                          is not specified and to compute rayleigh scattering if --add_rayleigh flag is specified')

parser.add_argument('--x_res',
                    default=0.01,
                    type=np.float32,
                    help='Radiance sampling resolution in x axis (North)')

parser.add_argument('--y_res',
                    default=0.01,
                    type=np.float32,
                    help='Radiance sampling resolution in y axis (East)')

parser.add_argument('--azimuth',
                    default=[0.0],
                    nargs='+',
                    type=np.float32,
                    help='Azimuth angles for the radiance measurements [deg]'\
                         '90 is for measuring radiance exiting along Y axis (East)')

parser.add_argument('--zenith',
                    default=[0.0],
                    nargs='+',
                    type=np.float32,
                    help='Zenith angles for the radiance measurements [deg].' \
                         '0 is for measuring radiance exiting directly up.')

parser.add_argument('--n_jobs',
                    default=1,
                    type=int,
                    help='Number of jobs for parallel rendering. n_jobs=1 uses no parallelization')

# Additional arguments to the parser
subparser = argparse.ArgumentParser(add_help=False)
subparser.add_argument('--generator')
subparser.add_argument('--add_rayleigh', action='store_true')
parser.add_argument('--generator', 
                    help='Name of the generator used to generate the atmosphere. \
                          or additional generator arguments: python forward.py --generator=GENERATOR -h. \
                          A generator should have the following methods to generate the medium: \
                          update_parser, grid, extinction, albedo, phase. \
                          See scripts/generate/__init__.py for more documentation.')
parser.add_argument('--add_rayleigh',
                    action='store_true',
                    help='Overlay the atmosphere with (known) Rayleigh scattering due to air molecules. \
                          Temperature profile is taken from AFGL measurements of summer mid-lat.')

generator = subparser.parse_known_args()[0].generator 
add_rayleigh = subparser.parse_known_args()[0].add_rayleigh 
if generator:
    generator = importlib.import_module('scripts.generate.' + generator)
    parser = generator.update_parser(parser)

if add_rayleigh:
    parser = generate_air.update_parser(parser)
    
args = parser.parse_args()

assert len(args.azimuth) == len(args.zenith), 'Length of azimuth and zenith should be equal'
    

if __name__ == "__main__":

    # Define a Medium and initialize optical properties according to the Generator
    grid = generator.grid(args)
    extinction_c = generator.extinction(grid, args)
    albedo_c = generator.albedo(grid, args)
    phase_c = generator.phase(grid, args)
    
    cloud = shdom.Medium()
    cloud.set_optical_properties(extinction_c, albedo_c, phase_c)
    
    
    # Rayleigh scattering for air molecules
    atmosphere = cloud
    if args.add_rayleigh:
        atmosphere += generate_air.summer_mid_lat(args)
        
    # Define an RteSolver object and solve the Radiative Transfer for the domain
    scene_params = shdom.SceneParameters(source=shdom.SolarSource(args.solar_azimuth, args.solar_zenith))
    numerical_params = shdom.NumericalParameters()
    rte_solver = shdom.RteSolver(scene_params, numerical_params)
    rte_solver.init_medium(atmosphere)
    rte_solver.solve(maxiter=100)
    
    
    # Define a sensor and render an image of the domain.
    sensors = [
        shdom.OrthographicMonochromeSensor(
            wavelength=args.wavelength,
            bounding_box=cloud.bounding_box, 
            x_resolution=args.x_res, 
            y_resolution=args.y_res, 
            azimuth=azimuth, 
            zenith=zenith,
            altitude='TOA'
        ) for azimuth, zenith in zip(args.azimuth, args.zenith) 
    ]
    
    if len(sensors) > 1:
        sensors = shdom.SensorArray(sensors)
    else:
        sensors = sensors[0]
    
    if args.n_jobs == 1:
        images = sensors.render(rte_solver)
    elif args.n_jobs > 1:
        images = sensors.par_render(rte_solver, n_jobs=args.n_jobs)
    else:
        raise AssertionError('[assert] Number of jobs is less than 1.')
    
    measurements = shdom.Measurements(sensors, images=images)
    
    # Save measurements, groud-truth (cloud) medium and solver parameters
    shdom.save_forward_model(args.output_dir, cloud, rte_solver, measurements)
    
    
    
