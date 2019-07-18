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
import copy
 

class Scatterer(object):
    """TODO"""
    def __init__(self):
        self.grid = None
        
    @property
    def grid(self):
        return self._grid
    
    @grid.setter
    def grid(self, val):
        self._grid = val    
    
    @property
    def bounding_box(self):
        return self.grid.bounding_box
    
    
class OpticalScatterer(Scatterer):
    """TODO"""
    def __init__(self, wavelength, extinction=None, albedo=None, phase=None):
        super(OpticalScatterer, self).__init__()
        self._wavelength = round(wavelength, 3)
        if (extinction is not None) and (albedo is not None) and (phase is not None):
            self.grid = extinction.grid + albedo.grid + phase.grid
        self.extinction = extinction
        self.albedo = albedo
        self.phase = phase
        
    def resample(self, grid):
        """TODO"""
        extinction = self.extinction.resample(grid)
        albedo = self.albedo.resample(grid)
        phase = self.phase.resample(grid)            
        return shdom.OpticalScatterer(self.wavelength, extinction, albedo, phase)
        

    def get_mask(self, threshold):
        """
        Get a mask based on the optical extinction.
        
        Parameters
        ----------
        threshold: float
            A threshold which above this value it is considered a populated voxel.
        
        Returns
        -------
        mask: shdom.GridData object
            A boolean mask with True for dense voxels and False for optically thin regions.
        """
        data = self.extinction.data > threshold
        return shdom.GridData(self.grid, data)
    
    @property
    def wavelength(self):
        return self._wavelength    

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
           

class MicrophysicalScatterer(Scatterer):
    """TODO"""
    def __init__(self, lwc=None, reff=None, veff=None):
        super(MicrophysicalScatterer, self).__init__()
        self._mie = OrderedDict()
        self._wavelength = []
        self._min_reff = np.Inf
        self._max_reff = -np.Inf
        self._min_veff = np.Inf
        self._max_veff = -np.Inf           
        self.set_microphysics(lwc, reff, veff)
            
    def get_optical_scatterer(self, wavelength):
        """TODO"""
        wavelength = round(wavelength, 3)
        scatterer = shdom.OpticalScatterer(
            wavelength, 
            extinction=self.mie[wavelength].get_extinction(self.lwc, self.reff, self.veff), 
            albedo=self.mie[wavelength].get_albedo(self.reff, self.veff),
            phase=self.mie[wavelength].get_phase(self.reff, self.veff))
        return scatterer
    
    def resample(self, grid):
        """TODO"""
        lwc = self.lwc.resample(grid)
        reff = self.reff.resample(grid)
        veff = self.veff.resample(grid)            
        self.set_microphysics(lwc, reff, veff)
    
    def get_mask(self, threshold):
        """
        Get a mask based on the liquid water content.
    
        Parameters
        ----------
        threshold: float
            A threshold which above this value it is considered a populated voxel.
    
        Returns
        -------
        mask: shdom.GridData object
            A boolean mask with True for dense voxels and False for thin voxels.
        """
        data = self.lwc.data > threshold
        return shdom.GridData(self.grid, data)    

    def set_microphysics(self, lwc, reff, veff):
        """TODO"""
        self.lwc = lwc
        self.reff = reff
        self.veff = veff
        if (lwc is not None) and (reff is not None) and (veff is not None):
            self._grid = lwc.grid + reff.grid + veff.grid
 
    def add_mie(self, mie):
        """TODO"""
        if isinstance(mie, shdom.MiePolydisperse):
            mie_list = [mie]
        elif isinstance(mie, dict):
            mie_list = mie.values()
            
        for mie in mie_list:
            self._mie[mie.wavelength] = mie
            self._wavelength.append(mie.wavelength)
            self._min_reff = min(self.min_reff, mie.size_distribution.reff.min())
            self._max_reff = max(self.max_reff, mie.size_distribution.reff.max())
            self._min_veff = min(self.min_veff, mie.size_distribution.veff.min())
            self._max_veff = max(self.max_veff, mie.size_distribution.veff.max())        
       
        if self.reff is not None:            
            min_val = self.reff.data[self.reff.data>0.0].min() 
            assert  self.reff.max_value < self.max_reff+1e-3, \
                               'Maximum medium effective radius [{:2.2f}] is larger than the pre-computed table maximum radius [{:2.2f}]. ' \
                               'Recompute Mie table with larger maximum radius.'.format(self.reff.max_va, self.max_reff)
            assert  min_val > self.min_reff-1e-3, \
                    'Minimum medium effective radius [{:2.2f}] is smaller than the pre-computed table minimum radius [{:2.2f}]. ' \
                    'Recompute Mie table with smaller minimum radius.'.format(min_val, self.min_reff)   
        
        if self.veff is not None:
            min_val = self.veff.data[self.veff.data>0.0].min()            
            assert  self.veff.max_value < self.max_veff+1e-3, \
                    'Maximum medium effective radius [{:2.2f}] is larger than the pre-computed table maximum variance [{:2.2f}]. ' \
                    'Recompute Mie table with larger maximum radius.'.format(max_val, self.max_veff)
            assert  min_val > self.min_veff-1e-3, \
                    'Minimum medium effective radius [{:2.2f}] is smaller than the pre-computed table minimum variance [{:2.2f}]. ' \
                    'Recompute Mie table with smaller minimum radius.'.format(min_val, self.min_veff)              

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
        grid = self.load_grid(path)
        data = np.genfromtxt(path, skip_header=3)

        grid_index = data[:, :3].astype(int)
        lwc = data[:, 3]
        reff = data[:, 4]
        if data.shape[1] == 6:
            veff = data[:, 5]
        else:
            veff = veff * np.ones_like(reff)

        particle_levels = np.array([z in grid_index[:, 2] for z in range(grid.nz)], dtype=int)
        lwc_data  = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=np.nan)
        reff_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=np.nan)
        veff_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=np.nan)
        lwc_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]]  = lwc
        reff_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]] = reff
        veff_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]] = veff

        self.set_microphysics(
            lwc=shdom.GridData(grid, lwc_data).squeeze_dims(),
            reff=shdom.GridData(grid, reff_data).squeeze_dims(), 
            veff=shdom.GridData(grid, veff_data).squeeze_dims()
        )

    
    @property
    def lwc(self):
        return self._lwc
    
    @lwc.setter
    def lwc(self, val):
        if val is not None:
            assert val.min_value.astype(np.float32) >= 0.0, 'LWC should be larger than 0.0'
        self._lwc = val
    
    @property
    def reff(self):
        return self._reff    
    
    @reff.setter
    def reff(self, val):
        if val is not None and self.mie: 
            max_val = val.max_value
            min_val = val.data[val.data>0.0].min()
            assert  max_val < self.max_reff+1e-3, \
                    'Maximum medium effective radius [{:2.2f}] is larger than the pre-computed table maximum radius [{:2.2f}]. ' \
                    'Recompute Mie table with larger maximum radius.'.format(max_val, self.max_reff)
            assert  min_val > self.min_reff-1e-3, \
                    'Minimum medium effective radius [{:2.2f}] is smaller than the pre-computed table minimum radius [{:2.2f}]. ' \
                    'Recompute Mie table with smaller minimum radius.'.format(min_val, self.min_reff)            
        self._reff = val
     
    @property
    def veff(self):
        return self._veff    
        
    @veff.setter
    def veff(self, val):
        if val is not None and self.mie: 
            max_val = val.max_value
            min_val = val.data[val.data>0.0].min()            
            assert  max_val < self.max_veff+1e-3, \
                    'Maximum medium effective radius [{:2.2f}] is larger than the pre-computed table maximum variance [{:2.2f}]. ' \
                    'Recompute Mie table with larger maximum radius.'.format(max_val, self.max_veff)
            assert  min_val > self.min_veff-1e-3, \
                    'Minimum medium effective radius [{:2.2f}] is smaller than the pre-computed table minimum variance [{:2.2f}]. ' \
                    'Recompute Mie table with smaller minimum radius.'.format(min_val, self.min_veff)            
        self._veff = val      

    @property
    def mie(self):
        return self._mie
    
    @property
    def min_reff(self):
        return self._min_reff
    
    @property
    def max_reff(self):
        return self._max_reff 
    
    @property
    def min_veff(self):
        return self._min_veff 
    
    @property
    def max_veff(self):
        return self._max_veff

    @property
    def wavelength(self):
        if len(self._wavelength) == 0:
            return None
        if len(self._wavelength) == 1:
            return self._wavelength[0]
        else:
            return self._wavelength

class MultispectralScatterer(object):
    """TODO"""
    def __init__(self, scatterer_list=None):
        self._scatterer = OrderedDict()
        self._num_bands = 0
        self._grid = None
        self._wavelength = []
        if scatterer_list is not None:
            for scatterer in scatterer_list:
                self.add_scatterer(scatterer)
                
    def get_optical_scatterer(self, wavelength):
        """TODO"""
        wavelength = round(wavelength, 3)
        return self.scatterer[wavelength]    

    def resample(self, grid):
        """TODO"""
        for wavelength in self.scatterer.iterkeys():
            self.scatterer[wavelength] = self.scatterer[wavelength].resample(grid)    
        return self
    
    def add_scatterer(self, scatterer):
        if self.num_bands == 0:
            self._grid = scatterer.grid
        else:
            self._grid += scatterer.grid
        self._num_bands += 1
        self.scatterer[scatterer.wavelength] = scatterer.resample(self.grid)
        self._wavelength.append(scatterer.wavelength)
        
    @property
    def scatterer(self):
        return self._scatterer
    
    @property
    def num_bands(self):
        return self._num_bands    

    @property
    def grid(self):
        return self._grid    
    
    @property
    def wavelength(self):
        if len(self._wavelength) == 0:
            return None
        if len(self._wavelength) == 1:
            return self._wavelength[0]
        else:
            return self._wavelength    
    
    
class Medium(object):
    """
    The Medium object encapsulates an atmospheric optical medium with multiple Scatterers.
    This means specifying extinction, single scattering albedo and phase function at every point in the domain for every scatterer.
    """
    def __init__(self, grid=None):
        self._scatterers = OrderedDict()
        self._num_scatterers = 0
        self._wavelength = None
        self.set_grid(grid)
        
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
        Add an optical scatterer at some wavelength to the medium
        The extinction, single scattering albedo and phase function at every point in the domain 
        are provided by the scatterer model.
    
        Parameters
        ----------
        scatterer: shdom.Scatterer,
            A scattering particle distribution model
        """
        first_scatterer = True if self.num_scatterers==0 else False
        
        if first_scatterer:
            self._wavelength = scatterer.wavelength
        else:
            assert self.wavelength == scatterer.wavelength, ' medium wavelength {} differs from scatterer wavelength {}'.format(medium, scatterer)
        self._num_scatterers += 1
        name = 'scatterer{:d}'.format(self._num_scatterers) if name is None else name
        self.scatterers[name] = scatterer
   
   
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

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def grid(self):
        return self._grid       
    
    @property
    def scatterers(self):
        return self._scatterers
    
    @property
    def num_scatterers(self):
        return self._num_scatterers
    