"""
Generate: Homogeneous Cloud
---------------------------

This script defines a homogeneous Medium. A cloud mask can be provided so that the medium is only homogeneous inside the cloud.
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
    
    parser.add_argument('--extinction', 
                        default=1.0,
                        type=np.float32, 
                        help='(default value: %(default)s) Extinction of a homogeneous cloud [km^-1]')
    
    parser.add_argument('--albedo', 
                        default=1.0,
                        type=np.float32, 
                        help='(default value: %(default)s) Albedo of a homogeneous cloud in range [0, 1]')
    
    parser.add_argument('--reff', 
                        default=10.0,
                        type=np.float32, 
                        help='(default value: %(default)s) Effective radius of a homogeneous cloud [micron]')
    
    parser.add_argument('--mie_table_path', 
                        help='Path to a precomputed Mie scattering table. \
                              See notebooks/Make Mie Table.ipynb for more details')    
    return parser


def grid(args):
    """
    Define an atmospheric grid. 
    
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


def phase(grid, args, cloud_mask=None):
    """
    Get the Mie scattering phase function for the homogeneous cloud.
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
    cloud_mask: shdom.GridData object, optional
        A boolean mask on a grid specifying the cloudy voxels with True and non-cloudy with False.
        
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
    
    reff_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=args.reff, dtype=np.float32)
    
    if cloud_mask is not None:
        mask_grid = grid + cloud_mask.grid
        mask = np.array(cloud_mask.resample(mask_grid, method='nearest'), dtype=np.bool)
        reff_data[np.bitwise_not(mask)] = 0.0
        
    phase = mie.get_grid_phase(shdom.GridData(grid, reff_data))
    return phase
 
 
def extinction(grid, args, cloud_mask=None):
    """
    Parameters
    ----------
    grid: shdom.Grid object
        The atmospheric grid.
    args.extinction: float
        Extinction of the center voxel [km^-1]
    cloud_mask: shdom.GridData object, optional
        A boolean mask on a grid specifying the cloudy voxels with True and non-cloudy with False.
        
    Returns
    -------
    extinction: shdom.GridData object
        A GridData object with the optical extinction on a grid.
    """    
    
    ext_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=args.extinction, dtype=np.float32)
    
    if cloud_mask is not None:
        mask_grid = grid + cloud_mask.grid
        mask = np.array(cloud_mask.resample(mask_grid, method='nearest'), dtype=np.bool)
        ext_data[np.bitwise_not(mask)] = 0.0

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
    cloud_mask: shdom.GridData object, optional
        A boolean mask on a grid specifying the cloudy voxels with True and non-cloudy with False.
    
    Returns
    -------
    albedo: shdom.GridData object
        A GridData object with the single scattering albedo on a grid.
    """    
    alb_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=args.albedo, dtype=np.float32)
    
    if cloud_mask is not None:
        mask_grid = grid + cloud_mask.grid
        mask = np.array(cloud_mask.resample(mask_grid, method='nearest'), dtype=np.bool)
        alb_data[np.bitwise_not(mask)] = 0.0

    albedo = shdom.GridData(grid, alb_data)
    return albedo

