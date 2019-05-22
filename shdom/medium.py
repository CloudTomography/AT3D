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
import shdom
import dill as pickle
from collections import OrderedDict
 

class Scatterer(object):
    """TODO"""
    def __init__(self, extinction=None, albedo=None, phase=None):
        
        assert extinction.grid == albedo.grid == phase.grid, \
               'Different grids for phase, albedo and extinction is not supported.'        
        
        self.grid = extinction.grid
        self.extinction = extinction
        self.albedo = albedo
        self.phase = phase
        
        
    def resample(self, grid):
        """TODO"""
        if self.grid == grid:
            return self
        else:
            extinction = self.extinction.resample(grid)
            albedo = self.albedo.resample(grid)
            phase = self.phase.resample(grid)            
            return shdom.Scatterer(extinction, albedo, phase)
        
        
    @property
    def extinction(self):
        return self._extinction
    
    @extinction.setter
    def extinction(self, val):
        if val is not None:
            assert val.min_value.astype(np.float32) >= 0.0, 'Extinction should be larger than 0.0'
        self._extinction = val
    
    @property
    def albedo(self):
        return self._albedo    
    
    @albedo.setter
    def albedo(self, val):
        if val is not None: 
            assert (val.max_value.astype(np.float32) <= 1.0 and  val.min_value.astype(np.float32) >= 0.0), 'Single scattering albedo should be in the range [0, 1]'
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
    
    @property
    def bounding_box(self):
        return self.grid.bounding_box
    
    
class OpticalMedium(object):
    """
    The OpticalMedium object encapsulates an atmospheric optical medium with multiple Scatterers.
    This means specifying extinction, single scattering albedo and phase function at every point in the domain for every scatterer.
    """
    def __init__(self, grid):
        self._scatterers = OrderedDict()
        self._num_scatterers = 0
        self.set_grid(grid)
        self._extinctp = None
        self._albedop = None
        self._iphasep = None
        self._legendre_table = None

        
    def set_grid(self, grid):
        """
        TODO
        """
        self._grid = grid
        
    
    def get_scatterer(self, name):
        """
        TODO
        """
        return self.scatterers[name]
    
    
    def add_scatterer(self, scatterer, name=None):
        """
        Add an optical scatterer to the medium
        The extinction, single scattering albedo and phase function at every point in the domain 
        are provided by the scatterer model.
    
        Parameters
        ----------
        scatterer: shdom.Scatterer,
            A scattering particle distribution model
        """
        
        first_particle = True if self.num_scatterers==0 else False
        
        self._num_scatterers += 1
        if name is None:
            name = 'scatterer{:d}'.format(self._num_scatterers)
        self.scatterers[name] = scatterer.resample(self.grid)
        
        extinctp = self.scatterers[name].extinction.data[...,np.newaxis]
        albedop = self.scatterers[name].albedo.data[...,np.newaxis]
        legendre_table = self.scatterers[name].phase.legendre_table  
        iphasep = self.scatterers[name].phase.iphasep[...,np.newaxis]            

        if first_particle: 
            self._extinctp = extinctp
            self._albedop = albedop
            self._legendre_table = legendre_table
            self._iphasep = iphasep

        else:
            self._extinctp = np.append(self.extinctp, extinctp, axis=-1)
            self._albedop = np.append(self.albedop, albedop, axis=-1)
            self._iphasep = np.append(self.iphasep, iphasep + self.legendre_table.numphase, axis=-1)
            self._legendre_table.append(legendre_table)
               
        
    
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
        mask = shdom.GridData(self.grid, mask_data)
        self.extinction *= mask
        self.albedo *= mask
        if self.phase.type == 'Tabulated':
            self.phase._index *= mask
        elif self.phase.type == 'Grid':
            self.phase * mask
    

    @property
    def grid(self):
        return self._grid
    
    @property
    def extinctp(self):
        return self._extinctp
    
    @property
    def albedop(self):
        return self._albedop
    
    @property
    def iphasep(self):
        return self._iphasep    

    @property
    def legendre_table(self):
        return self._legendre_table
    
    @property
    def scatterers(self):
        return self._scatterers
    
    @property
    def num_scatterers(self):
        return self._num_scatterers
    
    
class MicrophysicalMedium(object):
    """
    A MicrophysicalMedium encapsulates microphysical properties on a grid
    
    Parameters
    ----------
    grid: shdom.Grid
        A grid for the microphysical medium.
    """
    def __init__(self, grid=None):
        self.set_grid(grid)
        self._lwc = None
        self._reff = None
        self._veff = None
        
        
    def set_lwc(self, lwc):
        """
        TODO
        """
        self._lwc = lwc.resample(self.grid)


    def set_reff(self, reff):
        """
        TODO
        """
        self._reff = reff.resample(self.grid)
        
        
    def set_veff(self, veff):
        """
        TODO
        """
        self._veff = veff.resample(self.grid)


    def set_grid(self, grid):
        """
        TODO
        """
        self._grid = grid  


    def load_grid(self, path):
        """
        A utility function to load a microphysical medium from file
        
        Parameters
        ----------
        path: str
             Path to file.
             
        Returns
        -------
        grid: shdom.Grid object
            The 3D grid of the medium.
            
        Notes
        -----
        CSV format should be as follows:
        
        # comment line (description)
        nx ny nz
        dz dy dz     z_levels[0]     z_levels[1] ...  z_levels[nz-1]
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
        .
        .
        .
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
        """
        nx, ny, nz = np.genfromtxt(path, max_rows=1, dtype=int) 
        dx, dy = np.genfromtxt(path, max_rows=1, usecols=(0, 1), dtype=float, skip_header=2)        
        z_grid = np.genfromtxt(path, max_rows=1, usecols=range(2, 2 + nz), dtype=float, skip_header=2)
        x_grid = np.linspace(0.0, (nx - 1)*dx, nx, dtype=np.float32)
        y_grid = np.linspace(0.0, (ny - 1)*dy, ny, dtype=np.float32)    
        grid = shdom.Grid(x=x_grid, y=y_grid, z=z_grid)
        return grid
        
    
    def save_to_csv(self, path, comment_line=''):
        """
        A utility function to save a microphysical medium.
        
        Parameters
        ----------
        path: str
            Path to file.
        comment_line: str, optional
            A comment line describing the file.
            
            
        Notes
        -----
        CSV format is as follows:
        
        # comment line (description)
        nx ny nz
        dz dy dz     z_levels[0]     z_levels[1] ...  z_levels[nz-1]
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
        .
        .
        .
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
        """
        np.savetxt(path, X=np.array([self.grid.shape]), fmt='%d', header=comment_line)
        f = open(path, 'ab')
        np.savetxt(f, X=np.concatenate((np.array([self.grid.dx, self.grid.dy]), self.grid.z)).reshape(1,-1), fmt='%2.3f')
        y, x, z = np.meshgrid(range(self.grid.nx), range(self.grid.ny), range(self.grid.nz))
        data = np.vstack((x.ravel(), y.ravel(), z.ravel(), self.lwc.data.ravel(), self.reff.data.ravel())).T
        np.savetxt(f, X=data, fmt='%d %d %d %.5f %.3f')        
        f.close()
        
        
    def load_from_csv(self, path, veff=0.1):
        """ 
        A utility function to load a microphysical medium.
        
        Parameters
        ----------
        path: str
            Path to file. 
        veff: float
            If effective variance is not specified in the csv file as a 6th column,
            this value is used as a homogeneous value. 
            Default value is veff=0.1
    
        Notes
        -----
        CSV format should be as follows:
        
        # comment line (description)
        nx ny nz
        dz dy dz     z_levels[0]     z_levels[1] ...  z_levels[nz-1]
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]  veff[ix, iy, iz](optional)
        .
        .
        .
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]  veff[ix, iy, iz](optional)
        """ 
        self.set_grid(self.load_grid(path))
        data = np.genfromtxt(path, skip_header=3)
        
        grid_index = data[:, :3].astype(int)
        lwc = data[:, 3]
        reff = data[:, 4]
        if data.shape[1] == 6:
            veff = data[:, 5]
        else:
            veff = veff * np.ones_like(reff)
        
        particle_levels = np.array([z in grid_index[:, 2] for z in range(self.grid.nz)], dtype=int)
        lwc_data  = np.zeros(shape=(self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float32)
        reff_data = np.zeros(shape=(self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float32)
        veff_data = np.zeros(shape=(self.grid.nx, self.grid.ny, self.grid.nz), dtype=np.float32)
        lwc_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]]  = lwc
        reff_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]] = reff
        veff_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]] = veff
        
        self.set_lwc(shdom.GridData(self.grid, lwc_data))
        self.set_reff(shdom.GridData(self.grid, reff_data))
        self.set_veff(shdom.GridData(self.grid, veff_data))

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