"""
Medium and Medium related objects used for atmospheric rendering.
A Medium objects is defined by its optical properties: extinction, albedo and phase function.
Two Medium objects can be added to create a combined Medium. 
This will result in a unified grid with a summation of the optical extinction and a weighted average of the albedo and phase function.
"""

import core
import numpy as np
from enum import Enum
import warnings
from shdom import Grid, GridData, BoundingBox, Rayleigh
import dill as pickle


class Medium(object):
    """
    The Medium object encapsulates an atmospheric optical medium. 
    This means specifying extinction, single scattering albedo and phase function at every point in the domain.
    """
    def __init__(self):
        self._type = 'Medium'
        self._extinction = None
        self._albedo = None
        self._phase = None        
        
        
    def set_optical_properties(self, extinction, albedo, phase):
        """
        Set Medium optical properties: extinction, single scattering albedo and phase function at every point in the domain.
    
        Parameters
        ----------
        extinction: GridData
        albedo: GridData
        phase: Phase (TabulatedPhase or GridPhase)
        
        Notes
        -----
        Different grids for extinction and albedo and phase is not supported.
        """        
        assert extinction.grid == albedo.grid == phase.grid, \
               'Different grids for phase, albedo and extinction is not supported.'        
        self.extinction = extinction
        self.albedo = albedo
        self.phase = phase
        
 
    def __add__(self, other):
        """Adding two Medium objects or a Medium and AmbientMedium."""
        extinction = self.extinction + other.extinction
        scat = self.albedo*self.extinction + other.albedo*other.extinction
        albedo = scat / extinction        
        if other.type == 'Medium':
            if self.phase.type == 'Tabulated' or other.phase.type == 'Tabulated':
                raise NotImplementedError('Medium adding of tabulated phases not implemented.')
            phase = (self.albedo*self.extinction*self.phase + other.albedo*other.extinction*other.phase) / scat
        if other.type == 'AmbientMedium':
            if self.phase.type == 'Grid' or other.phase.type == 'Grid':
                raise NotImplementedError('AmbientMedium adding of grid phases not implemented.')
            phase = self.phase.add_ambient(other.phase)
        medium = Medium()
        medium.set_optical_properties(extinction, albedo, phase)
        return medium
    
      
    def save(self, path):
        """
        Save Medium to file.
        
        Parameters
        ----------
        path: str,
            Full path to file. 
        """
        file = open(path,'w')
        file.write(pickle.dumps(self.__dict__, -1))
        file.close()
        
    
    def load(self, path):
        """
        Load Medium from file.
        
        Parameters
        ----------
        path: str,
            Full path to file. 
        """        
        file = open(path, 'r')
        data = file.read()
        file.close()
        self.__dict__ = pickle.loads(data)    


    def get_mask(self, threshold):
        """
        Get a cloud mask based on the optical extinction.
        
        Parameters
        ----------
        threshold: float
            A threshold which above this value it is considered a cloudy voxel.
        
        Returns
        -------
        mask: shdom.GridData object
            A boolean mask with True making cloudy voxels and False marking non-cloud region.
        """
        data = self.extinction.data > threshold
        return GridData(self.grid, data)
    
    
    def apply_mask(self, mask):
        """
        Zero down the medium properties where the cloud mask is False.
        
        Parameters
        ----------
        mask: shdom.GridData object
            A boolean mask with True making cloudy voxels and False marking non-cloud region.    
        """
        mask_data = np.array(mask.resample(self.grid, method='nearest'), dtype=np.float)
        mask = GridData(self.grid, mask_data)
        self.extinction *= mask
        self.albedo *= mask
        if self.phase.type == 'Tabulated':
            self.phase._index *= mask
        elif self.phase.type == 'Grid':
            self.phase * mask
    
    @property
    def extinction(self):
        return self._extinction
    
    @extinction.setter
    def extinction(self, val):
        self._extinction = val
        self.grid = self.extinction.grid
    
    @property
    def albedo(self):
        return self._albedo    
    
    @albedo.setter
    def albedo(self, val):
        assert (val.max_value <= 1.0 and  val.min_value >= 0.0), 'Single scattering albedo should be in the range [0, 1]'
        self._albedo = val
        
    @property
    def phase(self):
        return self._phase    
      
    @phase.setter
    def phase(self, val):
        self._phase = val
        
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
    """
    An AmbientMedium is defined in the same way a Medium is defined
    by its optical properties: extinction, albedo and phase function. 
    When it is added to a Medium it only `fills in the holes` where there is no medium density. 
    This approximation gives efficiency in memory and computations.
    """
    def __init__(self):
        super(AmbientMedium, self).__init__()
        self._type = 'AmbientMedium'
        
        
    def __add__(self, other):
        if other.type == 'AmbientMedium':
            extinction = self.extinction + other.extinction
            albedo = (self.albedo*self.extinction + other.albedo*other.extinction) / extinction
            medium = AmbientMedium()
            medium.set_optical_properties(extinction, albedo, self.phase)
        elif other.type == medium:
            medium = other + self
        else:
            raise NotImplementedError
        return medium
    
 
class MicrophysicalMedium(object):
    """
    A MicrophysicalMedium encapsulates microphysical properties on a grid
    
    Notes
    -----
    Currently effective variance is a scalar (homogeneous veff).
    """
    def __init__(self, veff=0.1):
        self._lwc = None
        self._reff = None
        self._veff = None
        self._grid = None
    
    def get_grid(self, path):
        """
        A utility function to load Large Eddy Simulated clouds.
        
        Parameters
        ----------
        path: str
             Path to file.
             
        Returns
        -------
        grid: Grid object
            The 3D grid of the LES generated data.
            
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
        nx, ny, nz = np.genfromtxt(path, max_rows=1, dtype=int) 
        dx, dy = np.genfromtxt(path, max_rows=1, usecols=(0, 1), dtype=float, skip_header=4)        
        z_grid = np.genfromtxt(path, max_rows=1, usecols=range(2, 2 + nz), dtype=float, skip_header=4)
        x_grid = np.linspace(0.0, (nx - 1)*dx, nx, dtype=np.float32)
        y_grid = np.linspace(0.0, (ny - 1)*dy, ny, dtype=np.float32)    
        grid = Grid(x=x_grid, y=y_grid, z=z_grid)
        return grid
        
        
    def load_from_csv(self, path):
        """ 
        A utility function to load Large Eddy Simulated clouds.
        
        Parameters
        ----------
        path: str
             Path to file. 
    
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
        
        self._grid = self.get_grid(path)
        grid_index = np.genfromtxt(path, usecols=(0, 1, 2), dtype=int, skip_header=5)
        lwc = np.genfromtxt(path, usecols=3, dtype=float, skip_header=5)
        reff = np.genfromtxt(path, usecols=4, dtype=float, skip_header=5)
        
        particle_levels = np.array([z in grid_index[:, 2] for z in range(self.grid.nz)], dtype=int)
        lwc_data  = np.zeros(shape=(self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float32)
        reff_data = np.zeros(shape=(self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float32)
        lwc_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]]  = lwc
        reff_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]] = reff
        self._lwc = GridData(self.grid, lwc_data)
        self._reff = GridData(self.grid, reff_data)

    @property
    def grid(self):
        return self._grid
    
    @property
    def reff(self):
        return self._reff
    
    @property
    def veff(self):
        return self._veff
    
    @property
    def lwc(self):
        return self._lwc