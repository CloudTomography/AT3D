"""
Parameters are used by the Optimizer to minimize the loss function.
These are the parameters we seek to recover. To recover a parameter is has to depends on (at least one):
Extinction, Single Scattering Albedo or Phase Function.
"""

from shdom import GridData
import numpy as np


class Parameter(GridData):
    def __init__(self):
        self._type = 'AbstractParameter'
        
        # Dependency on medium optical parameters
        self._extinction_dependency = None
        self._albedo_dependency = None
        self._phase_dependency = None
        self._min_bound = None
        self._max_bound = None
        self._mask = None
        self._num_parameters = None
        
    def initialize(self, init, min_bound=None, max_bound=None):
        super(Parameter, self).__init__(init.grid, init.data)
        self._min_bound = min_bound
        self._max_bound = max_bound
        self._num_parameters = self.data.size
        
    def get_extinction(self):
        return None
    
    def get_albedo(self):
        return None    

    def get_phase(self):
        return None 
    
    def set_mask(self, mask):
        self._mask = np.array(mask.resample(self.grid, method='nearest'), dtype=np.bool)
        self._num_parameters = self.mask.sum()
    
    def set_data(self, data):
        if self.mask is None:
            self._data = data
        else:
            self._data[self.mask] = data
                    
    @property
    def type(self):
        return self._type
    
    @property
    def extinction_dependency(self):
        return self._extinction_dependency    
    
    @property
    def albedo_dependency(self):
        return self._albedo_dependency

    @property
    def phase_dependency(self):
        return self._phase_dependency
    
    
    @property
    def min_bound(self):
        return self._min_bound    

    @property
    def max_bound(self):
        return self._max_bound
    
    @property 
    def bounds(self):
        return [(self.min_bound, self.max_bound)] * self.num_parameters
    
    @property
    def mask(self):
        return self._mask
    
    @property
    def num_parameters(self):
        return self._num_parameters
    
class Extinction(Parameter):
    def __init__(self):
        super(Extinction, self).__init__()
        self._type = 'Extinction'
        self._extinction_dependency = True
        self._albedo_dependancy = False
        self._phase_dependency = False
        
    def get_extinction(self):
        return self
    