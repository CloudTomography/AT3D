""" 
Render: [Monochromatic] Radiance at Top of the Atmosphere (TOA)
---------------------------------------------------------------
Forward rendering of an atmospheric medium with an monochromatic orthographic sensor measuring exitting radiance at the top of the the domain.
This sensor is an (somewhat crude) approximation for far observing satellites where the rays are parallel. 

As with all `render` scripts a `generator` needs to be specified using the generator flag (--generator). 
The Generator defines the medium parameters: Grid, Extinction, Single Scattering Albedo and Phase function with it's own set of command-line flags.

For example usage see the README.md
  
For information about the command line flags see:
  python scripts/render/render_monochromatic_radiance_toa.py --help
  
For a tutorial overview of how to operate the forward rendering model see the following notebooks:
 - notebooks/Make Mie Table.ipynb
 - notebooks/Single Image Rendering.ipynb
 - notebooks/Multiview Rendering.ipynb
"""

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
    
    parser.add_argument('output_dir', 
                        help='Path to an output directory where the measureents and model parameters will be saved. \
                              If the folder doesnt exist, it will be created.')
    parser.add_argument('--solar_zenith', 
                        default=165.0,
                        type=np.float32, 
                        help='(default value: %(default)s) Solar zenith [deg]. This is the direction of the photons in range (90, 180]')
    parser.add_argument('--solar_azimuth', 
                        default=0.0,
                        type=np.float32,
                        help='(default value: %(default)s) Solar azimuth [deg]. This is the direction of the photons')
    parser.add_argument('--x_res',
                        default=0.01,
                        type=np.float32,
                        help='(default value: %(default)s) Radiance sampling resolution in x axis (North)')
    parser.add_argument('--y_res',
                        default=0.01,
                        type=np.float32,
                        help='(default value: %(default)s) Radiance sampling resolution in y axis (East)')
    parser.add_argument('--azimuth',
                        default=[0.0],
                        nargs='+',
                        type=np.float32,
                        help='(default value: %(default)s) Azimuth angles for the radiance measurements [deg]'\
                             '90 is for measuring radiance exiting along Y axis (East)')
    parser.add_argument('--zenith',
                        default=[0.0],
                        nargs='+',
                        type=np.float32,
                        help='(default value: %(default)s) Zenith angles for the radiance measurements [deg].' \
                             '0 is for measuring radiance exiting directly up.')   
    parser.add_argument('--n_jobs',
                        default=1,
                        type=int,
                        help='(default value: %(default)s) Number of jobs for parallel rendering. n_jobs=1 uses no parallelization')
    
    # Additional arguments to the parser
    subparser = argparse.ArgumentParser(add_help=False)
    subparser.add_argument('--generator')
    subparser.add_argument('--add_rayleigh', action='store_true')
    parser.add_argument('--generator', 
                        help='Name of the generator used to generate the atmosphere. \
                              or additional generator arguments: python scripts/render/radiace_toa.py --generator GENERATOR -h. \
                              See shdom/generate.py for more documentation.')
    parser.add_argument('--add_rayleigh',
                        action='store_true',
                        help='Overlay the atmosphere with (known) Rayleigh scattering due to air molecules. \
                              Temperature profile is taken from AFGL measurements of summer mid-lat.')
    
    add_rayleigh = subparser.parse_known_args()[0].add_rayleigh 
    generator = subparser.parse_known_args()[0].generator
    
    if generator:
        CloudGenerator = getattr(shdom.Generate, generator)
        parser = CloudGenerator.update_parser(parser)
    
    AirGenerator = None
    if add_rayleigh:
        AirGenerator = shdom.Generate.AFGLSummerMidLatAir
        parser = AirGenerator.update_parser(parser)
        
    args = parser.parse_args()
    assert len(args.azimuth) == len(args.zenith), 'Length of azimuth and zenith should be equal'
    
    return args, CloudGenerator, AirGenerator

    
def generate_atmosphere():
    """
    Generate an atmospheric domain for rendering.
    """
    cloud_generator = CloudGenerator(args)
    cloud_generator.add_mie(args.mie_table_path)
    cloud = cloud_generator.get_scatterer()

    atmosphere = shdom.Medium()
    if args.add_rayleigh:
        air_generator = AirGenerator(args)
        air = air_generator.get_scatterer(cloud.wavelength)    
        atmosphere.set_grid(cloud.grid + air.grid)
        atmosphere.add_scatterer(air, 'air')
    else:
        atmosphere.set_grid(cloud.grid)

    atmosphere.add_scatterer(cloud, 'cloud')

    return atmosphere


def solve_rte(atmosphere):
    """
    Define an RteSolver object and solve the Radiative Transfer for the domain.
    
    Parameters
    ----------
    atmosphere: shdom.Medium
        The atmospheric Medium (used for rendering)
        
    Returns
    -------
    rte_solver: shdom.RteSolver object
        A solver with the solution for the RTE
    """
    scene_params = shdom.SceneParameters(
        source=shdom.SolarSource(args.solar_azimuth, args.solar_zenith)
    )
    numerical_params = shdom.NumericalParameters()
    rte_solver = shdom.RteSolver(scene_params, numerical_params)
    rte_solver.set_medium(atmosphere)
    rte_solver.solve(maxiter=100) 
    return rte_solver
    
    
def render(bounding_box, rte_solver):
    """
    Define a sensor and render an orthographic image at the top domain.
    
    Parameters
    ----------
    bounding_box: shdom.BoundingBox object
        Used to compute the projection that will see the entire bounding box.
    rte_solver: shdom.RteSolver
        A solver that contains the solution to the RTE in the medium
        
    Returns
    -------
    measurments: shdom.Measurements object
        Encapsulates the radiances and the sensor geometry for later optimization
    """  
    projection = shdom.MultiViewProjection()
    for azimuth, zenith in zip(args.azimuth, args.zenith):
        projection.add_projection(
            shdom.OrthographicProjection(
                bounding_box=bounding_box, 
                x_resolution=args.x_res, 
                y_resolution=args.y_res, 
                azimuth=azimuth, 
                zenith=zenith,
                altitude='TOA'))

    camera = shdom.Camera(shdom.RadianceSensor(), projection)
    images = camera.render(rte_solver, args.n_jobs)
    measurements = shdom.Measurements(camera, images=images)
    return measurements


if __name__ == "__main__":
    args, CloudGenerator, AirGenerator = argument_parsing()
    atmosphere = generate_atmosphere()
    rte_solver = solve_rte(atmosphere)
    measurements = render(atmosphere.get_scatterer('cloud').bounding_box, rte_solver)

    # Save measurements, medium and solver parameters
    shdom.save_forward_model(args.output_dir, atmosphere, rte_solver, measurements)
    
    
    
