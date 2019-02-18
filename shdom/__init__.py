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
    def __init__(self, grid, data):
        self._grid = grid
        self._data = data
        
        points = []
        if grid.nx > 1:
            points.append(grid.x_grid)          
        if grid.ny > 1:
            points.append(grid.y_grid)        
        points.append(grid.z_levels)
            
        self._data_interpolator = RegularGridInterpolator(points,
                                                          values=np.squeeze(data), 
                                                          bounds_error=False, 
                                                          fill_value=0.0)
  
    def resample(self, new_grid):
        """TODO"""
        
        points = []
        if self.grid.nx > 1:
            points.append(new_grid.x_grid)          
        if self.grid.ny > 1:
            points.append(new_grid.y_grid)
        points.append(new_grid.z_levels)
        data = self._data_interpolator(np.stack(np.meshgrid(*points, indexing='ij'), axis=-1))
        
        if (self.grid.ny == 1) or (new_grid.ny == 1):
            data = data[np.newaxis,...] 
        if (self.grid.nx == 1) or (new_grid.nx == 1):
            data = data[np.newaxis,...] 
        
        return data

        
    def __add__(self, other):
        grid = self.grid + other.grid
        data = self.resample(grid) + other.resample(grid)
        return ScalarField(grid, data)
    
    
    def __mul__(self, other):  
        grid = self.grid + other.grid

        if other.__class__ is VectorField:
            data = self.resample(grid)[..., np.newaxis] * other.resample(grid)
            return VectorField(grid, data)
        
        elif other.__class__ is ScalarField:
            data = self.resample(grid) * other.resample(grid)
            return ScalarField(grid, data)
        else:
            print('Error multiplying {} * {}'.format(self.__class__, other.__class__))        


    def __div__(self, other):
        grid = self.grid + other.grid
        
        if other.__class__ is VectorField:
            data = self.resample(grid)[..., np.newaxis] / other.resample(grid)
            return VectorField(grid, data)
        elif other.__class__ is ScalarField:
            data = self.resample(grid) / other.resample(grid)
            return ScalarField(grid, data)        
        else:
            print('Error dividing {} / {}'.format(self.__class__, other.__class__))
    
    
    @property
    def grid(self):
        return self._grid
    
    @property
    def data(self):
        return self._data
    
    @property
    def shape(self):
        return self._data.shape
        
        
class VectorField(ScalarField):
    """
    TODO: add documentation

    Parameters
    ----------
    grid: Grid 
    
    Notes
    -----
    """
    def __init__(self, grid, data):
        super(VectorField, self).__init__(grid, data)
        self._depth = data.shape[3] - 1

    
    def __add__(self, other):
        grid = self.grid + other.grid
        self_data = self.resample(grid)
        other_data = other.resample(grid)
        depth_diff = self.depth - other.depth
        if depth_diff > 0:
            self_data = self_data
            other_data = np.pad(other_data, ((0, 0), (0, 0), (0, 0), (0, depth_diff)), 'constant')
        elif depth_diff < 0:
            self_data = np.pad(self_data, ((0, 0), (0, 0), (0, 0), (0, -depth_diff)), 'constant')
            other_data = other_data

        return VectorField(grid, self_data +  other_data)

    
    def __mul__(self, other): 
        grid = self.grid + other.grid
        other_data = other.resample(grid)
        if other.__class__ is ScalarField:
            other_data = other_data[..., np.newaxis]
        elif other.__class__ is not VectorField:
            print('Error multiplying {} * {}'.format(self.__class__, other.__class__))        
        data = self.resample(grid) * other_data
        return VectorField(grid, data)
    
    def __div__(self, other):  
        grid = self.grid + other.grid
        other_data = other.resample(grid)        
        if other.__class__ is ScalarField:
            other_data = other_data[..., np.newaxis]
        elif other.__class__ is not VectorField:
            print('Error dividing {} / {}'.format(self.__class__, other.__class__))        
        data = self.resample(grid) / other_data
        return VectorField(grid, data)
    
    @property
    def depth(self):
        return self._depth 

class Grid(object):
    """ 
    An grid objects. x and y must have even spacings.
    z can have uneven spacing by specifying the z_levels parameter.
    
    Parameters
    ----------
    bounding_box: BoundingBox.
        A BoundingBox object. Cloud be unbound in x and y.
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
        self._num_grid_points = self.nx * self.ny *self.nz
        if z_levels is None:
            self._z_levels = np.linspace(bounding_box.zmin, bounding_box.zmax, nz) 
        else:
            self._z_levels = z_levels
        self._dx = (self.bounding_box.xmax - self.bounding_box.xmin) / self.nx
        self._dy = (self.bounding_box.ymax - self.bounding_box.ymin) / self.ny
        self._dz = (self.bounding_box.zmax - self.bounding_box.zmin) / self.nz
        self._x_grid = np.linspace(self.bounding_box.xmin, self.bounding_box.xmax, self.nx) 
        self._y_grid = np.linspace(self.bounding_box.ymin, self.bounding_box.ymax, self.ny) 
            

    def __add__(self, other):
        """ TODO: keeps dx, dy constant"""
    
        bounding_box = self.bounding_box + other.bounding_box
    
        x_size = bounding_box.xmax - bounding_box.xmin
        y_size = bounding_box.ymax - bounding_box.ymin
        dx = min(self.dx, other.dx)
        dy = min(self.dy, other.dy)
    
        if np.isfinite(x_size):
            nx = int(x_size / dx)
        else:
            nx = 1
    
        if np.isfinite(y_size):
            ny = int(y_size / dy)
        else:
            ny = 1 
    
        # Bottom part of the atmosphere (no grid intersection)
        if self.z_levels[0] < other.z_levels[0]:
            z_bottom = self.z_levels[self.z_levels < other.z_levels[0]]
        else:
            z_bottom = other.z_levels[other.z_levels < self.z_levels[0]]
    
        # Top part of the atmosphere (no grid intersection)
        if self.z_levels[-1] < other.z_levels[-1]:
            z_top = other.z_levels[other.z_levels > self.z_levels[-1]]
        else:
            z_top = self.z_levels[self.z_levels > other.z_levels[-1]]
    
        # Middle part of the atmosphere (grids intersect)
        z_middle_self = self.z_levels
        z_middle_other = other.z_levels
        if z_bottom.any():
            z_middle_self = self.z_levels[self.z_levels > z_bottom[-1]]
            z_middle_other = other.z_levels[other.z_levels > z_bottom[-1]]
        if z_top.any():
            z_middle_self = self.z_levels[self.z_levels < z_top[0]]
            z_middle_other = other.z_levels[other.z_levels < z_top[0]]
    
        z_middle = z_middle_self if len(z_middle_self) > len(z_middle_other) else z_middle_other
    
        # Check if an extra point is necessary at the bottom 
        if z_bottom.any() & len(z_middle)>2:
            extra_zlevel = 2*z_middle[0] - z_middle[1]
            if extra_zlevel > z_bottom[-1]:
                z_middle = np.append(extra_zlevel, z_middle)
    
        # Check if an extra point is necessary at the top 
        if z_top.any() & len(z_middle)>2:
            extra_zlevel = 2*z_middle[-1] - z_middle[-2]
            if extra_zlevel < z_top[0]:
                z_middle = np.append(z_middle, extra_zlevel)
    
        z_levels = np.concatenate((z_bottom, z_middle, z_top))
        nz = len(z_levels)
        
        return Grid(bounding_box, nx, ny, nz, z_levels)        

    def __eq__(self, other) : 
        for item1, item2 in zip(self.__dict__.itervalues(), other.__dict__.itervalues()):
            if not np.array_equal(np.nan_to_num(item1), np.nan_to_num(item2)):
                return False
        return True
    
    
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
    def num_grid_points(self):
        return self._num_grid_points
    
    @property 
    def z_levels(self):
        return self._z_levels
    
    @property 
    def bounding_box(self):
        return self._bb
    
    @property 
    def dx(self):
        return self._dx  
    
    @property 
    def dy(self):
        return self._dy    

    @property 
    def dz(self):  
        return self._dz 
    
    @property 
    def x_grid(self):
        return self._x_grid   
    
    @property 
    def y_grid(self):
        return self._y_grid      

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
        assert xmin < xmax, 'Zero area bounding_box along x axis.'  
        assert ymin < ymax, 'Zero area bounding_box along y axis.' 
        assert zmin < zmax, 'Zero area bounding_box along z axis.'       
        self.xmin = xmin
        self.ymin = ymin
        self.zmin = zmin
        self.xmax = xmax
        self.ymax = ymax
        self.zmax = zmax
        
    def __eq__(self, other) : 
        return self.__dict__ == other.__dict__    
    
    
    def __add__(self, other):
        
        xmin = self.xmin
        if np.isfinite(other.xmin):
            if np.isfinite(xmin):
                xmin = min(xmin, other.xmin)
            else:
                xmin = other.xmin
            
        ymin = self.ymin
        if np.isfinite(other.ymin):
            if np.isfinite(ymin):
                ymin = min(ymin, other.ymin)
            else:
                ymin = other.ymin
                
        zmin = self.zmin
        if np.isfinite(other.zmin):
            if np.isfinite(zmin):
                zmin = min(zmin, other.zmin)
            else:
                zmin = other.zmin
        
        xmax = self.xmax
        if np.isfinite(other.xmax):
            if np.isfinite(xmax):
                xmax = max(xmax, other.xmax)
            else:
                xmax = other.xmax
            
        ymax = self.ymax
        if np.isfinite(other.ymax):
            if np.isfinite(ymax):
                ymax = max(ymax, other.ymax)
            else:
                ymax = other.ymax
                
        zmax = self.zmax
        if np.isfinite(other.zmax):
            if np.isfinite(zmax):
                zmax = min(zmax, other.zmax)
            else:
                zmax = other.zmax
                
        return BoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)

 
from phase import *
from medium import *
from sensor import *
from rte_solver import *