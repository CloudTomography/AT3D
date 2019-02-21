"""
A medium object used for atmospheric rendering and inversion.
"""

import core
import numpy as np
from enum import Enum
import warnings
from shdom import Grid, GridData, BoundingBox, Rayleigh


class Medium(object):
    """
    The Medium object encapsulates an atmospheric optical medium. 
    This means specifying extinction, single scattering albedo and phase function at every point in the domain.
    TODO: phase type and sum of phases

    Parameters
    ----------
    extinction: GridData
    albedo: GridData
    phase: Phase
    
    Notes
    -----
    Different grids for extinction and albedo is not supported.
    """
    
    def __init__(self, extinction, albedo, phase):
        assert extinction.grid == albedo.grid == phase.grid, \
               'Different grids for phase, albedo and extinction is not supported.'        
        self.extinction = extinction
        self._albedo = albedo
        self._phase = phase
        self._type = 'Medium'
        
 
    def __add__(self, other):
        extinction = self.extinction + other.extinction
        scat = self.albedo*self.extinction + other.albedo*other.extinction
        albedo = scat / extinction        
        if other.type == 'Medium':
            phase = (self.albedo*self.extinction*self.phase + other.albedo*other.extinction*other.phase) / scat
        if other.type == 'AmbientMedium':
            phase = self.phase.add_ambient(other.phase)  
        return Medium(extinction, albedo, phase)
        
    @property
    def extinction(self):
        return self._ext
    
    @extinction.setter
    def extinction(self, val):
        self._ext = val
        self.grid = self.extinction.grid
    
    @property
    def albedo(self):
        return self._albedo    
    
    @property
    def phase(self):
        return self._phase    
        
    @property
    def grid(self):
        return self._grid
    
    @grid.setter
    def grid(self, val):
        self._grid = val
        self._bounding_box = self.grid.bounding_box
    
    @property
    def bounding_box(self):
        return self._bounding_box
    
    @property
    def type(self):
        return self._type
    
    
    
    
class AmbientMedium(Medium):
    """TODO"""
    def __init__(self, extinction, albedo, phase):

        super(AmbientMedium, self).__init__(extinction, albedo, phase)
        self._type = 'AmbientMedium'
        
        
    def __add__(self, other):
        if other.type == 'AmbientMedium':
            extinction = self.extinction + other.extinction
            albedo = (self.albedo*self.extinction + other.albedo*other.extinction) / extinction
            medium = AmbientMedium(extinction, albedo, self.phase)
        elif other.type == medium:
            medium = other + self
        else:
            raise NotImplementedError
        return medium
    
    
def load_les_from_csv(path_to_csv):
    """ 
    A utility function to load Large Eddy Simulated clouds.
    
    Parameters
    ----------
    path_to_csv: str
         Path to file. 

    Returns
    -------
    lwc: GridData
         a GridData object contatining the liquid water content of the LES cloud.
    reff: GridData
          a GridData object contatining the effective radius of the LES cloud.
    
    Notes
    -----
    CSV format should be as follows:
    
    #name=name of les file
    #original_cloud_data=path to original 
    #resampled_cloud_data_grid_size=grid resolution in meters
    nx ny nz
    dz dy dz     z_levels[0]     z_levels[1] ...  z_levels[nz-1]
    ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
    .
    .
    .
    ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
    """ 
    
    nx, ny, nz = np.genfromtxt(path_to_csv, max_rows=1, dtype=int) 
    dx, dy = np.genfromtxt(path_to_csv, max_rows=1, usecols=(0, 1), dtype=float, skip_header=4)
    z_grid = np.genfromtxt(path_to_csv, max_rows=1, usecols=range(2, 2 + nz), dtype=float, skip_header=4)
    grid_index = np.genfromtxt(path_to_csv, usecols=(0, 1, 2), dtype=int, skip_header=5)
    lwc = np.genfromtxt(path_to_csv, usecols=3, dtype=float, skip_header=5)
    reff = np.genfromtxt(path_to_csv, usecols=4, dtype=float, skip_header=5)
    
    particle_levels = np.array([z in grid_index[:, 2] for z in range(nz)], dtype=int)
    lwc_data  = np.zeros(shape=(nx, ny, nz), dtype=np.float32)
    reff_data = np.zeros(shape=(nx, ny, nz), dtype=np.float32)
    lwc_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]]  = lwc
    reff_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]] = reff
    
    x_grid = np.linspace(0.0, (nx - 1)*dx, nx, dtype=np.float32)
    y_grid = np.linspace(0.0, (ny - 1)*dy, ny, dtype=np.float32)
    grid = Grid(x=x_grid, y=y_grid, z=z_grid)
    return GridData(grid, lwc_data), GridData(grid, reff_data)