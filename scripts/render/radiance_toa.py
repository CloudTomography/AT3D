""" 
Render: Radiance at top-of-atmosphere
-------------------------------------

python scripts/render/render_radiance_toa.py --generator generate_single_voxel --mie_table_path mie_tables/Water_672nm.scat experiments/single_voxel/ --add_rayleigh --air_max_altitude 1.0 --air_num_grid_points 5

TODO

Define and render an atmospheric volume of 1x1x1 km^3 with a single cloudy voxel (Mie scattering) in the center. 
All parameters are saved for later recovery and analysis. 

Example usage:
  python scripts/render/.py --extinction=15.0 --albedo=0.8 --reff=14 --solar_zenith=170 \
                           --mie_table_path='mie_tables/Water_672nm.scat' --output_dir=experiments/single_voxel
  
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
import importlib
from generate import generate_air


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

# Additional arguments of the generator
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
    generator = importlib.import_module('generate.' + generator)
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
    if args.add_rayleigh:
        air = shdom.Medium()
        extinction_a, albedo_a, phase_a = generate_air.summer_mid_lat(args)
        air.set_optical_properties(extinction_a, albedo_a, phase_a) 
        atmosphere = cloud + air
    else:
        atmosphere = cloud
        
    # Define an RteSolver object and solve the Radiative Transfer for the domain
    scene_params = shdom.SceneParameters(source=shdom.SolarSource(args.solar_azimuth, args.solar_zenith))
    numerical_params = shdom.NumericalParameters()
    rte_solver = shdom.RteSolver(scene_params, numerical_params)
    rte_solver.init_medium(atmosphere)
    rte_solver.solve(maxiter=100)
    
    
    # Define a sensor and render an image of the domain.
    sensors = [
        shdom.OrthographicSensor(bounding_box=cloud.bounding_box, 
                                 x_resolution=args.x_res, 
                                 y_resolution=args.y_res, 
                                 azimuth=azimuth, 
                                 zenith=zenith,
                                 altitude='TOA') 
        for azimuth, zenith in zip(args.azimuth, args.zenith)
    ]
    
    if len(sensors) > 1:
        sensors = shdom.SensorArray(sensors)
    else:
        sensors = sensors[0]
    
    if args.n_jobs == 1:
        radiances = sensors.render(rte_solver)
    elif args.n_jobs > 1:
        radiances = sensors.par_render(rte_solver, n_jobs=args.n_jobs)
    else:
        raise AssertionError('[assert] Number of jobs is less than 1.')
    
    measurements = shdom.Measurements(sensors, radiances)
    
    # Save measurements, groud-truth medium and rte solver parameters
    shdom.save_forward_model(args.output_dir, atmosphere, rte_solver, measurements)
    
    
    
