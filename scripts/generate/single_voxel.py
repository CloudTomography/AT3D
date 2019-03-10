"""
Generate: Single Voxel
----------------------

This script defines a Medium with a single voxel in center of the grid. 
It is useful for developing and debugging of derivatives and gradients and sensitivity analysis.
"""

import os 
import numpy as np
import argparse
import shdom


def update_parser(parser): 
    """
    Update the argument parser with parameters relavent to this generation script. 
    
    Parameters
    ----------
    parser: argparse.ArgumentParser()
        The main parser to update.

    Returns
    -------
    parser: argparse.ArgumentParser()
        The updated parser.

    """    
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
    
    parser.add_argument('--domain_size', 
                        default=1.0,
                        type=float, 
                        help='(default value: %(default)s) Cubic domain size [km]')
    
    parser.add_argument('--extinction', 
                        default=10.0,
                        type=np.float32, 
                        help='(default value: %(default)s) Extinction of the center voxel [km^-1]')
    
    parser.add_argument('--albedo', 
                        default=1.0,
                        type=np.float32, 
                        help='(default value: %(default)s) Albedo of the center voxel in range [0, 1]')
    
    parser.add_argument('--reff', 
                        default=10.0,
                        type=np.float32, 
                        help='(default value: %(default)s) Effective radius of the center voxel [micron]')
    
    parser.add_argument('--mie_table_path', 
                        help='Path to a precomputed Mie scattering table. \
                              See notebooks/Make Mie Table.ipynb for more details')
    
    return parser


def grid(args):
    """
    Define an atmospheric grid of 1x1x1 km^3. 
    
    Parameters
    ----------
    args.nx, int
        The number of grid points in x (North)
    args.ny, int
        The number of grid points in y (East)
    args.nz, int
        The number of grid points in z (Up)
    args.domain_size, float
       The length of the cubic domain along one axis [km].
    """
    bb = shdom.BoundingBox(0.0, 0.0, 0.0, args.domain_size, args.domain_size, args.domain_size)
    grid = shdom.Grid(bounding_box=bb, nx=args.nx, ny=args.ny, nz=args.nz)
    return grid


def phase(grid, args):
    """
    Get the Mie scattering phase function for the unknown voxel.
    If the table_path doesnt exist or is not specified, a table is computed for the wavelength requiered. 
    
    Parameters
    ----------
    grid: shdom.Grid object
        A shdom.Grid object specifing the atmospheric grid.
    args.mie_table_path: str
        Path to a precomputed Mie scattering table. See notebooks/Make Mie Table.ipynb for more details.
    args.reff: float 
        The effective radius for which to retrieve the phase function [microns].
    args.wavelength: float
        The wavelength [microns]. It is used to create Mie scattering talbe if a path is unspecified.
    
    Returns
    -------
    phase: GridPhase object
         The phase function specified for the unknown voxel on a 3D grid.
    """
    mie = shdom.Mie()
    if args.mie_table_path:
        mie.read_table(args.mie_table_path)
    else:
        mie.set_parameters((args.wavelength, args.wavelength), 'Water', 'gamma', 7.0)
        mie.compute_table(1, args.reff, args.reff, 65.0)   
    
    reff_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32) 
    reff_data[grid.nx/2, grid.ny/2, grid.nz/2] = args.reff
    phase = mie.get_grid_phase(shdom.GridData(grid, reff_data))
    return phase
 
 
def extinction(grid, args):
    """
    Parameters
    ----------
    grid: shdom.Grid object
        The atmospheric grid.
    args.extinction: float
        Extinction of the center voxel [km^-1]
        
    Returns
    -------
    extinction: shdom.GridData object
        A GridData object with the optical extinction on a grid.
    """    
    ext_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
    ext_data[grid.nx/2, grid.ny/2, grid.nz/2] = args.extinction
    extinction = shdom.GridData(grid, ext_data)
    return extinction


def albedo(grid, args):
    """
    Parameters
    ----------
    grid: shdom.Grid object
        The atmospheric grid.
    args.albedo: float
        Albedo of the center voxel in range [0, 1].
    
    Returns
    -------
    albedo: shdom.GridData object
        A GridData object with the single scattering albedo on a grid.
    """    
    alb_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
    alb_data[grid.nx/2, grid.ny/2, grid.nz/2] = args.albedo
    albedo = shdom.GridData(grid, alb_data)
    return albedo

