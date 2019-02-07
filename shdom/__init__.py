"""
Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer

This is a python wrapper for SHDOM created by Aviad Levis and Amit Aides, Technion Inst. of Technology.
The purpose of this wrapper is to develop 3D remote sensing metoglogies. 

The documentation of the source Fortran code by Frank Evans can be found at
http://nit.colorado.edu/shdom/shdomdoc

Information about the source code taken from the documentation page:
This program computes unpolarized monochromatic or spectral band radiative transfer in a one, two,
or three-dimensional medium for either collimated solar and/or thermal emission sources of radiation.
The properties of the medium can be specified completely generally, i.e. the extinction, single 
scattering albedo, Legendre coefficients of the scattering phase function, and temperature for
the particular wavelength or spectral band may be specified at each input grid point. SHDOM is
superior to Monte Carlo radiative transfer methods when many radiative quantities are desired,
e.g. the radiance field across the domain top or the 3D distribution of heating. Radiances at
any angle, hemispheric fluxes, net fluxes, mean radiances, and net flux convergence (related
to heating rates) may be output anywhere in the domain. For highly peaked phase functions the 
delta-M method may be chosen, in which case the radiance is computed with an untruncated phase
function single scattering correction. A correlated k-distribution approach is used for the
integration over a spectral band. There may be uniform or spatially variable Lambertian
reflection and emission from the ground surface. Several types of bidirectional reflection
distribution functions (BRDF) for the surface are implemented, and more may be added easily.
SHDOM may be run on a single processor or on multiple processors (e.g. an SMP machine or a
cluster) using the Message Passing Interface (MPI).
"""


class ScalarField(object):
    """
    TODO: add documentation

    Parameters
    ----------
    grid: Grid 
    
    Notes
    -----
    """
    def __init__(self, grid, field):
        self._grid = grid
        self._field = field
    
    @property
    def grid(self):
        return self._grid
    
    @property
    def field(self):
        return self._field
    
    
class VectorField(ScalarField):
    """
    TODO: add documentation

    Parameters
    ----------
    grid: Grid 
    
    Notes
    -----
    """
    def __init__(self, grid, field):
        super(VectorField, self).__init__(grid, field)
        self._depth = field.shape[3] 
    
    @property
    def depth(self):
        return self._depth    
    
class Grid(object):
    """ 
    An grid objects. x and y must have even spacings.
    z can have uneven spacing by specifying the z_levels parameter.
    
    Parameters
    ----------
    nx: integer
        Number of grid point in x axis.
    ny: integer
        Number of grid point in z axis.
    nz: integer
        Number of grid point in z axis.
    z_levels: np.array(dtype=float, shape=(nz,)), optional
        grid along z axis
    Notes
    -----
    """
    def __init__(self, bounding_box, nx, ny, nz, z_levels=None):
        self._bb = bounding_box
        self._nx = nx
        self._ny = ny
        self._nz = nz
        if z_levels is None:
            self._z_levels = np.linspace(bounding_box.zmin, bounding_box.zmax, nz) 
        else:
            self._z_levels = z_levels
        
    @property 
    def nx(self):
        return self._nx
    
    @property 
    def ny(self):
        return self._ny
    
    @property 
    def nz(self):
        return self._nz
    
    @property 
    def z_levels(self):
        return self._z_levels
    
    @property 
    def bounding_box(self):
        return self._bb
    
    @property 
    def dx(self):
        if hasattr(self, '_dx') is False:
            self._dx = (self._bb.xmax - self._bb.xmin) / self._nx
        return self._dx  
    
    @property 
    def dy(self):
        if hasattr(self, '_dy') is False:
            self._dy = (self._bb.ymax - self._bb.ymin) / self._ny
        return self._dy    

    @property 
    def dz(self):
        if hasattr(self, '_dz') is False:
            self._dz = (self._bb.zmax - self._bb.zmin) / self._nz  
        return self._dz 
    
            
class BoundingBox(object):
    """ 
    A bounding box object.
    
    Parameters
    ----------
    xmin: float
         Minimum x (North).
    ymin: float
         Minimum y (East)
    zmin: float
         Minimum z (Up).
    xmax: float
         Maximum x (North).
    ymax: float
         Maximum y (East).
    zmax: float
         Maximum z (Up).

    
    Notes
    -----
    All values are in [km] units
    """ 
    def __init__(self, xmin, ymin, zmin, xmax, ymax, zmax):
        self.xmin = xmin
        self.ymin = ymin
        self.zmin = zmin
        self.xmax = xmax
        self.ymax = ymax
        self.ymin = zmax
        
        
        

from mie import *
from medium import *
from rte_solver import *    