""" 
Forward Rendering: single voxel
-------------------------------
Define and render an atmospheric volume of 1x1x1 km^3 with a single cloudy voxel (Mie scattering) in the center. 
All parameters are saved for later recovery and analysis. 

Example usage:
  python scripts/render.py --extinction=15.0 --albedo=0.8 --reff=14 --solar_zenith=170 \
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

parser = argparse.ArgumentParser()

parser.add_argument('--output_dir', 
                    help='Path to an output directory where the forward modeling parameters will be saved. \
                          If the folder doesnt exist, it will be created.')

parser.add_argument('--nx', 
                    default=3,
                    type=int, 
                    help='(default value: %(default)s) Number of grid cell in x (North) direction')

parser.add_argument('--ny', 
                    default=3,
                    type=int, 
                    help='(default value: %(default)s) Number of grid cell in y (East) direction')

parser.add_argument('--nz', 
                    default=3,
                    type=int, 
                    help='(default value: %(default)s) Number of grid cell in z (Up) direction')

parser.add_argument('--extinction', 
                    default=10.0,
                    type=np.float32, 
                    help='(default value: %(default)s) Extinction of the unknown voxel [km^-1]')

parser.add_argument('--albedo', 
                    default=1.0,
                    type=np.float32, 
                    help='(default value: %(default)s) Albedo of the unknown voxel in range [0, 1]')

parser.add_argument('--reff', 
                    default=10.0,
                    type=np.float32, 
                    help='(default value: %(default)s) Effective radius of the unknown voxel [micron]')

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

parser.add_argument('--mie_table_path', 
                    help='Path to a precomputed Mie scattering table. \
                          See notebooks/Make Mie Table.ipynb for more details')

parser.add_argument('--add_rayleigh',
                    action='store_true',
                    help='Overlay the atmosphere with (known) Rayleigh scattering due to air molecules.')

args = parser.parse_args()

def get_grid(nx, ny, nz):
    """
    Define an atmospheric grid of 1x1x1 km^3. 
    
    Parameters
    ----------
    nx, int
        The number of grid points in x (North)
    ny, int
        The number of grid points in y (East)
    nz, int
        The number of grid points in z (Up)
    """
    bb = shdom.BoundingBox(-0.5, -0.5, 0.0, 0.5, 0.5, 1.0)
    grid = shdom.Grid(bounding_box=bb, nx=nx, ny=ny, nz=nz)
    return grid


def get_phase(grid, reff, table_path, wavelength):
    """
    Get the Mie scattering phase function for the unknown voxel.
    If the table_path doesnt exist or is not specified, a table is computed for the wavelength requiered. 
    
    Parameters
    ----------
    grid: shdom.Grid object
        A shdom.Grid object specifing the atmospheric grid.
    table_path: str
        Path to a precomputed Mie scattering table. See notebooks/Make Mie Table.ipynb for more details.
    reff: float 
        The effective radius for which to retrieve the phase function [microns].
    wavelength: float
        The wavelength [microns]. It is used to create Mie scattering talbe if a path is unspecified.
    
    Returns
    -------
    phase: GridPhase object
         The phase function specified for the unknown voxel on a 3D grid.
    """
    
    mie = shdom.Mie()
    if os.path.exists(table_path):
        mie.read_table(table_path)
    else:
        mie.set_parameters((wavelength, wavelength), 'Water', 'gamma', 7.0)
        mie.compute_table(1, reff, reff, 65.0)   
    
    reff_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32) 
    reff_data[int((grid.nx+1)/2), int((grid.ny+1)/2), int((grid.nz+1)/2)] = reff
    phase = mie.get_grid_phase(shdom.GridData(grid, reff_data))
    return phase
 
 
def get_extinction(grid, extinction):
    """
    Parameters
    ----------
    grid: shdom.Grid object
        The atmospheric grid.
    
    Returns
    -------
    extinction: shdom.GridData object
        A GridData object with the optical extinction on a grid.
    """    
    ext_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
    ext_data[int((grid.nx+1)/2), int((grid.ny+1)/2), int((grid.nz+1)/2)] = extinction
    extinction = shdom.GridData(grid, ext_data)
    return extinction


def get_albedo(grid, albedo):
    """
    Parameters
    ----------
    grid: shdom.Grid object
        The atmospheric grid.
    
    Returns
    -------
    albedo: shdom.GridData object
        A GridData object with the single scattering albedo on a grid.
    """    
    alb_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
    alb_data[int((grid.nx+1)/2), int((grid.ny+1)/2), int((grid.nz+1)/2)] = albedo
    albedo = shdom.GridData(grid, alb_data)
    return albedo

def save_model(directory, medium, sensor, solver, measurements):
    """
    Save the forward model parameters for reconstruction.
    
    Parameters
    ----------
    directory: str
        Directory path where the forward modeling parameters are saved. 
        If the folder doesnt exist it will be created.
    medium: shdom.Medium object
        The atmospheric medium. This ground-truth medium will be used to 
    sensor: shdom.Sensor object
        The sensor used to image the medium.
    solver: shdom.RteSolver object
        The solver and the parameters used. This includes the scene parameters (such as solar and surface parameters)
        and the numerical parameters.
    
    Notes
    -----
    The ground-truth atmosphere is later used for evaulation of the recovery.
    """  
    if not os.path.isdir(directory):
        os.makedirs(directory)    
    np.save(file=os.path.join(directory, 'measurements.npy'), arr=measurements)
    atmosphere.save(path=os.path.join(directory, 'ground_truth_medium'))
    sensor.save(path=os.path.join(directory, 'sensor'))
    rte_solver.save_params(path=os.path.join(directory, 'solver_parameters'))   
    
    
if __name__ == "__main__":

    # Define a Medium and initialize it's optical properties
    grid = get_grid(args.nx, args.ny, args.nz)
    extinction_c = get_extinction(grid, args.extinction)
    albedo_c = get_albedo(grid, args.albedo)
    phase_c = get_phase(grid, args.reff, args.mie_table_path, args.wavelength)
    atmosphere = shdom.Medium()
    atmosphere.set_optical_properties(extinction_c, albedo_c, phase_c)
    
    
    # Define an RteSolver object and solve the Radiative Transfer for the domain
    rte_solver = shdom.RteSolver()
    scene_params = shdom.SceneParameters(source=shdom.SolarSource(args.solar_azimuth, args.solar_zenith))
    numerical_params = shdom.NumericalParameters()
    rte_solver.set_params(scene_params, numerical_params)
    rte_solver.init_medium(atmosphere)
    rte_solver.solve(maxiter=100)
    
    
    # Define a sensor and render an image of the domain. 
    sensor = shdom.OrthographicSensor(bounding_box=atmosphere.bounding_box, 
                                      x_resolution=0.1, 
                                      y_resolution=0.1, 
                                      azimuth=0.0, 
                                      zenith=0.0,
                                      altitude='TOA')
    measurements = sensor.render(rte_solver)
    
    
    # Save measurements, groud-truth medium and rte solver parameters
    save_model(args.output_dir, atmosphere, sensor, rte_solver, measurements)
    
    
    # Rayleigh scattering for air molecules
    if args.add_rayleigh:
        temperatures = np.array([292.220, 292.040, 291.860, 291.680, 291.500, 291.320, 291.140, 290.960, 290.780, 
                                 290.600, 290.420, 290.240, 290.060, 289.880, 289.700, 289.920, 290.140, 290.360, 
                                 290.580, 290.800, 291.020, 291.240, 291.460, 291.680, 291.900])
        temp_grid = shdom.Grid(z=np.linspace(0.0, 20.0, len(temperatures)))
        temperature_profile = shdom.GridData(temp_grid, temperatures)
        rayleigh = shdom.Rayleigh(wavelength=args.wavelength, temperature_profile=temperature_profile)
        extinction_a, albedo_a, phase_a = rayleigh.get_scattering_field(grid.z_grid, phase_type='Grid')
        air.set_optical_properties(extinction_a, albedo_a, phase_a)
