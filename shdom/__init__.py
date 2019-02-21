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
from scipy.interpolate import interp1d, RegularGridInterpolator


class Grid(object):
    """ 
    TODO 
    
    Parameters
    ----------
    data: TODO
    bounding_box: BoundingBox.
        A BoundingBox object. Cloud be unbound in x and y.
    nx: integer
        Number of grid point in x axis.
    ny: integer
        Number of grid point in z axis.
    z: np.array(dtype=float, shape=(nz,)), optional
        grid along z axis
    """    
    def __init__(self, **kwargs):
        
        # 3D grid with grids
        if kwargs.has_key('x') and kwargs.has_key('y') and kwargs.has_key('z'):
            self._type = '3D'
            self.x = kwargs['x']
            self.y = kwargs['y']
            self.z = kwargs['z']
            self._bounding_box = BoundingBox(self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax)
            
        # 3D grid with bounding box
        elif kwargs.has_key('bounding_box') and kwargs.has_key('nx') and kwargs.has_key('ny'):
            self._type = '3D'
            bb = kwargs['bounding_box']
            self._bounding_box = bb
            nx, ny = kwargs['nx'], kwargs['ny']
            self.x = np.linspace(bb.xmin, bb.xmax, nx)
            self.y = np.linspace(bb.ymin, bb.ymax, ny)
            if kwargs.has_key('z'):
                self.z = kwargs['z']
            elif kwargs.has_key('nz'):
                nz = kwargs['nx']
                self.z = np.linspace(bb.zmin, bb.zmax, nz)
            else:
                raise AttributeError('z or nz are missing')
            
        # 1D grid 
        elif kwargs.has_key('z'):
            self._type = '1D'
            self.z = kwargs['z']
            self._nx = self._ny = 1
            self._x = self._y = self._bounding_box = None
         
        else:
            raise AttributeError('kwargs in Grid initialization are not defined')
    
    def get_common_x_grid(self, other):
        """
        TODO
        """
        if self.type == '1D' and other.type == '1D':
            return None       
        if self.type == '3D' and other.type == '1D':
            return self.x
        if self.type == '1D' and other.type == '3D':
            return other.x
        
        if np.allclose(self.x, other.x):
            return self.x

        xmax = max(self.xmax, other.xmax)
        xmin = min(self.xmin, other.xmin) 
        x_size = xmax - xmin      
        dx = min(self.dx, other.dx)
        nx = int(y_size / dx)        
        return np.linspace(xmin, xmax, nx, dtype=np.float32)  
    
    
    def get_common_y_grid(self, other):
        """
        TODO
        """
        if self.type == '1D' and other.type == '1D':
            return None       
        if self.type == '3D' and other.type == '1D':
            return self.y
        if self.type == '1D' and other.type == '3D':
            return other.y
        
        if np.allclose(self.y, other.y):
            return self.y
        
        ymax = max(self.ymax, other.ymax)
        ymin = min(self.ymin, other.ymin) 
        y_size = ymax - ymin      
        dy = min(self.dy, other.dy)
        ny = int(y_size / dy)        
        return np.linspace(ymin, ymax, ny, dtype=np.float32)              
       
        
    def get_common_z_grid(self, other):

        """TODO"""
        
        if np.allclose(self.z, other.z):
            return self.z
        
        # Bottom part of the atmosphere (no grid intersection)
        if self.zmin < other.zmin:
            z_bottom = self.z[self.z < other.zmin]
        else:
            z_bottom = other.z[other.z < self.zmin]
    
        # Top part of the atmosphere (no grid intersection)
        if self.zmax < other.zmax:
            z_top = other.z[other.z > self.zmax]
        else:
            z_top = self.z[self.z > other.zmax]
    
        # Middle part of the atmosphere (grids intersect)
        z_middle_self = self.z
        z_middle_other = other.z
        if z_bottom.any():
            z_middle_self = self.z[self.z > z_bottom[-1]]
            z_middle_other = other.z[other.z > z_bottom[-1]]
        if z_top.any():
            z_middle_self = self.z[self.z < z_top[0]]
            z_middle_other = other.z[other.z < z_top[0]]
    
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
    
        return np.concatenate((z_bottom, z_middle, z_top))      

    def __add__(self, other):
        """ TODO: keeps dx, dy constant"""
        x_grid = self.get_common_x_grid(other)
        y_grid = self.get_common_y_grid(other)
        z_grid = self.get_common_z_grid(other)
        if x_grid is not None and y_grid is not None:
            grid = Grid(x=x_grid, y=y_grid, z=z_grid)
        else:
            grid = Grid(z=z_grid)
        return grid
    
    
    def __eq__(self, other) : 
        for item1, item2 in zip(self.__dict__.itervalues(), other.__dict__.itervalues()):
            if not np.array_equal(np.nan_to_num(item1), np.nan_to_num(item2)):
                return False
        return True      

    @property
    def type(self):
        return self._type

    @property
    def x(self):
        return self._x
    
    @x.setter
    def x(self, val):
        val = np.array(val, dtype=np.float32)
        spacing = np.diff(val)
        assert np.all(np.isclose(spacing, spacing[0])), 'x grid supoprt equally spacing only'
        self._x = val
        self._dx = spacing[0]  
        self._nx = len(val)
        self._xmin, self._xmax = val[0], val[-1]
    
    @property
    def y(self):
        return self._y
    
    @y.setter
    def y(self, val):
        val = np.array(val, dtype=np.float32)
        spacing = np.diff(val)
        assert np.all(np.isclose(spacing, spacing[0])), 'y grid supoprt equally spacing only'
        self._y = val
        self._dy = spacing[0] 
        self._ny = len(val)
        self._ymin, self._ymax = val[0], val[-1]

    @property
    def z(self):
        return self._z                 
    
    @z.setter
    def z(self, val):
        val = np.array(val, dtype=np.float32)
        self._z = val
        self._nz = len(val)
        self._zmin, self._zmax = val[0], val[-1]

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
    def shape(self):
        if self.type == '1D':
            return (self.nz,)
        else:
            return (self.nx, self.ny, self.nz)
    
    @property
    def num_points(self):
        return self.nx * self.ny * self.nz
    
    @property 
    def dx(self):
        return self._dx
    
    @property 
    def dy(self):
        return self._dy    
    
    @property 
    def xmin(self):
        return self._xmin
    
    @property 
    def ymin(self):
        return self._ymin

    @property 
    def zmin(self):
        return self._zmin
    
    @property 
    def xmax(self):
        return self._xmax
    
    @property 
    def ymax(self):
        return self._ymax    

    @property 
    def zmax(self):
        return self._zmax    

    @property
    def bounding_box(self):
        return self._bounding_box
    
    
class GridData(object):
    """ 
    TODO 
    
    Parameters
    ----------
    data: TODO
    bounding_box: BoundingBox.
        A BoundingBox object. Cloud be unbound in x and y.
    nx: integer
        Number of grid point in x axis.
    ny: integer
        Number of grid point in z axis.
    z_levels: np.array(dtype=float, shape=(nz,)), optional
        grid along z axis
    """    
    def __init__(self, grid, data):
        self._type = grid.type
        self._grid = grid
        self._data = data
        self._shape = self._data.shape[:3]
        self._ndim = self._data.ndim        
        if self.type == '1D' and self.ndim is not 1:
            raise AttributeError('Grid is 1D but data dimension is:{}'.format(self.ndim))
        if self.type == '3D' and self.ndim < 3:
            raise AttributeError('Grid is 3D but data dimension is:{}'.format(self.ndim))
        
        assert self.shape == grid.shape, 'Data shape is {}, grid shape is {}'.format(self.shape, grid.shape)

        self._linear_interpolator1d = interp1d(grid.z, self.data, assume_sorted=True, copy=False, bounds_error=False, fill_value=0.0) 
        self._nearest_interpolator1d = interp1d(grid.z, self.data, assume_sorted=True, kind='nearest', copy=False, bounds_error=False, fill_value=0)
        if self.type == '3D':
            self._linear_interpolator3d = RegularGridInterpolator((grid.x, grid.y, grid.z), self.data, bounds_error=False, fill_value=0.0)
            self._nearest_interpolator3d = RegularGridInterpolator((grid.x, grid.y, grid.z), self.data, method='nearest', bounds_error=False, fill_value=0)
    
    
    def __add__(self, other):
        """ TODO: keeps dx, dy constant"""
        grid = self.grid + other.grid
        data = self.resample(grid) + other.resample(grid)
        return GridData(grid, data)
    
    
    def __sub__(self, other):
        """ TODO: keeps dx, dy constant"""
        grid = self.grid + other.grid
        data = self.resample(grid) - other.resample(grid)
        return GridData(grid, data)       
    
    
    def __mul__(self, other):
        """ TODO: keeps dx, dy constant"""
        if other.__class__ is GridPhase:
            result = other * self
        else:
            grid = self.grid + other.grid
            data = self.resample(grid) * other.resample(grid)
            result = GridData(grid, data)  
        return result
    
    
    def __div__(self, other):
        """ TODO: keeps dx, dy constant"""
        if other.__class__ is GridPhase:
            result = other * self
        else:
            grid = self.grid + other.grid
            data = self.resample(grid) / other.resample(grid)
            result = GridData(grid, data) 
        return result
    
    
    def resample(self, grid, method='linear'):
        """TODO"""
        if self.type == '1D' or (np.allclose(self.grid.x, grid.x) and np.allclose(self.grid.y, grid.y)):
            if method == 'linear':
                data = self._linear_interpolator1d(grid.z)
            elif method == 'nearest':
                data = self._nearest_interpolator1d(grid.z)
        else:
            if method == 'linear':
                data = self._linear_interpolator3d(np.stack(np.meshgrid(grid.x, grid.y, grid.z, indexing='ij'), axis=-1))
            elif method == 'nearest':
                data = self._nearest_interpolator3d(np.stack(np.meshgrid(grid.x, grid.y, grid.z, indexing='ij'), axis=-1)) 
        return data
    
    @property
    def grid(self):
        return self._grid
    
    @property
    def data(self):
        return self._data
    
    @property
    def shape(self):
        return self._shape
    
    @property
    def ndim(self):
        return self._ndim      
    
    @property
    def type(self):
        return self._type    


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