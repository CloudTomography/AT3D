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
import warnings
import numpy as np

def float_round(x):
    """Round a float or np.float32 to a 3 digits float"""
    if type(x) == np.float32:
        x = x.item()
    return round(x,3) 

def int_round(x):
    """Round a float or np.float32 to a 3 digits integer by 1000x scaling"""
    return int(np.round(x*1000))

def find_nearest(array, value):
    """Find the nearest element index in an array"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


class Grid(object):
    """ 
    A Grid object defining the 3D or 1D grid of the atmopshere. 
    
    A 3D Grid can be defined with:
      1. x, y, z grids.
      2. A BoundingBox and grid resolution (nx, ny, nz).
    
    A 1D grid is defined with a z grid.
    
    Parameters
    ----------
    bounding_box: BoundingBox, optional
        A BoundingBox object for 3D grid. If specified nx, ny and nz must be specified as well.
    nx: integer, optional
        Number of grid point in x axis. Must be specified with a BoundingBox object.
    ny: integer, optional
        Number of grid point in y axis. Must be specified with a BoundingBox object.
    nz: integer, optional
        Number of grid point in z axis. Must be specified with a BoundingBox object.
    x: np.array(dtype=float, shape=(nx,)), optional
        Grid along x axis. Must be specified with y,z grids.
    y: np.array(dtype=float, shape=(ny,)), optional
        Grid along y axis. Must be specified with x,z grids. 
    z: np.array(dtype=float, shape=(nz,)), optional
        Grid along z axis. Either specified with x,y grids (3D grid) or by itself (1D grid).
    """    
    
    def __init__(self, **kwargs):
        self._type = self.get_grid_type(kwargs)
        
        if 'bounding_box' in kwargs:
            self._bounding_box = kwargs['bounding_box']
        else: 
            self._bounding_box = None
            
        if self.type == '3D':
            if 'z' in kwargs:
                self.z = kwargs['z']
            else:
                self.z = np.linspace(kwargs['bounding_box'].zmin, kwargs['bounding_box'].zmax, kwargs['nz'])
                
            if 'x' in kwargs and 'y' in kwargs:
                self.x = kwargs['x']
                self.y = kwargs['y']
                self._bounding_box = BoundingBox(self.xmin, self.ymin, self.zmin, self.xmax, self.ymax, self.zmax)
            else:
                self.x = np.linspace(kwargs['bounding_box'].xmin, kwargs['bounding_box'].xmax, kwargs['nx'])
                self.y = np.linspace(kwargs['bounding_box'].ymin, kwargs['bounding_box'].ymax, kwargs['ny'])

        elif self.type == '1D':
            self.z = kwargs['z']
            self._nx = self._ny = 1
            self._x = self._y = None
        
        elif self.type == 'Homogeneous':
            self._nx = self._ny = self._nz = 1
            self._x = self._y = self._z = None         
            

    
    def get_grid_type(self, kwargs):
        """
        Retrieve the grid type.
        
        Parameters
        ----------
        kwargs: dict
           1. kwargs = {'x', 'y', 'z'} --> grid_type = '3D'
           2. kwargs = {'nx', 'ny', 'z'} --> grid_type = '3D'
           3. kwargs = {'nx', 'ny', 'nz', 'bounding_box'} --> grid_type = '3D'
           4. kwargs = {'z'} --> grid_type = '1D'
           5. else --> grid_type = 'Homogeneous'
        
        Returns
        -------
        grid_type: str
            type could be one of the following: '3D', '1D', 'Homogeneous'
        """
        if 'x' in kwargs and 'y' in kwargs and 'z' in kwargs or \
           'nx' in kwargs and 'ny' in kwargs and 'nz' in kwargs and 'bounding_box' in kwargs or \
           'nx' in kwargs and 'ny' in kwargs and 'z' in kwargs:
            grid_type = '3D'
        
        elif 'z' in kwargs: 
            grid_type = '1D'
            
        else: 
            grid_type = 'Homogeneous'
        return grid_type
            

    def get_common_x(self, other):
        """
        Find the common x which maintains a minimum dx (distance between two grid points).
        
        Parameters
        ----------
        other: shdom.Grid object
           The other for which to find a common x.
           
        Returns
        -------
        x: np.array(dtype=np.float32)
            The common x Grid.
        """
        if other.type == '1D' or other.type == 'Homogeneous' or np.array_equiv(self.x, other.x):
            return self.x       
        elif self.type == '1D' or self.type == 'Homogeneous':
            return other.x

        xmax = max(self.xmax, other.xmax)
        xmin = min(self.xmin, other.xmin) 
        x_size = xmax - xmin      
        dx = min(self.dx, other.dx)
        nx = int(x_size / dx)        
        return np.linspace(xmin, xmax, nx, dtype=np.float32)  
    
    
    def get_common_y(self, other):
        """
        Find the common y which maintains a minimum dy (distance between two grid points).
        
        Parameters
        ----------
        other: shdom.Grid object
           The other grid for which to find a common y grid.
           
        Returns
        -------
        y: np.array(dtype=np.float32)
            The common y.
        """        
        if other.type == '1D' or other.type == 'Homogeneous' or np.array_equiv(self.y, other.y):
            return self.y       
        elif self.type == '1D' or self.type == 'Homogeneous':
            return other.y
        
        if np.array_equiv(self.y, other.y):
            return self.y
        
        ymax = max(self.ymax, other.ymax)
        ymin = min(self.ymin, other.ymin) 
        y_size = ymax - ymin      
        dy = min(self.dy, other.dy)
        ny = int(y_size / dy)        
        return np.linspace(ymin, ymax, ny, dtype=np.float32)              
       
        
    def get_common_z(self, other):
        """
        Find the common z which maintains the high resolution z grid.
        
        Parameters
        ----------
        other: shdom.Grid object
           The other grid for which to find a common y grid.
           
        Returns
        -------
        z: np.array(dtype=np.float32)
            The common z.
        """        
        if other.type == 'Homogeneous' or np.array_equiv(self.z, other.z):
            return self.z  
        elif self.type == 'Homogeneous':
            return other.z
        
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
        """Add two grids by finding the common grid which maintains the higher resolution grid."""
        x = self.get_common_x(other)
        y = self.get_common_y(other)
        z = self.get_common_z(other)
        
        if x is not None and y is not None and z is not None:
            grid = Grid(x=x, y=y, z=z)
        else:
            if self.bounding_box is None:
                bounding_box = other.bounding_box
            elif other.bounding_box is None:
                bounding_box = self.bounding_box
            else:
                bounding_box = self.bounding_box + other.bounding_box
            if x is None and y is None and z is not None:
                grid = Grid(bounding_box=bounding_box, z=z)
            elif x is None and y is None and z is None:
                grid = Grid(bounding_box=bounding_box)
        return grid
    
    
    def __eq__(self, other):
        """Compare two grid objects."""
        for key, item in self.__dict__.items():
            other_item = None
            if key in other.__dict__:
                other_item = other.__dict__[key]
            if not np.array_equiv(np.nan_to_num(item), np.nan_to_num(other_item)):
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
        if self.type == 'Homogeneous':
            return ()
        elif self.type == '1D':
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
    A container for scalar fields which are defined on a Grid. 
    
    Parameters
    ----------
    grid: shdom.Grid object
        A Grid object of type '1D' or '3D'.
    data: np.array
        data contains a scalar field.
    """    
    def __init__(self, grid, data):
        self._type = grid.type
        self._grid = grid
        self._data = np.array(data).squeeze()
        self._shape = self._data.shape[:3]
        self._ndim = self._data.ndim    
        if self.type == 'Homogeneous' and self.ndim is not 0:
            raise AttributeError('Grid is Homogeneous but data dimension is:{}'.format(self.ndim))
        if self.type == '1D' and self.ndim is not 1:
            raise AttributeError('Grid is 1D but data dimension is:{}'.format(self.ndim))
        if self.type == '3D' and self.ndim < 3:
            raise AttributeError('Grid is 3D but data dimension is:{}'.format(self.ndim))
        
        assert self.shape == grid.shape, 'Data shape is {}, grid shape is {}'.format(self.shape, grid.shape)

        if self.type != 'Homogeneous':
            self._linear_interpolator1d = interp1d(grid.z, self.data, assume_sorted=True, copy=False, bounds_error=False, fill_value=0.0) 
            self._nearest_interpolator1d = interp1d(grid.z, self.data, assume_sorted=True, kind='nearest', copy=False, bounds_error=False, fill_value=0)
        if self.type == '3D':
            self._linear_interpolator3d = RegularGridInterpolator((grid.x, grid.y, grid.z), self.data, bounds_error=False, fill_value=0.0)
            self._nearest_interpolator3d = RegularGridInterpolator((grid.x, grid.y, grid.z), self.data, method='nearest', bounds_error=False, fill_value=0)
       
    def __eq__(self, other):
        """check if two GridData objects are equal"""
        if np.allclose(self.data, other.data) and self.grid == other.grid:
            return True
        else:
            return False
        
    def __add__(self, other):
        """Add two GridData objects by resampling to a common grid."""
        if self.grid == other.grid:
            grid = self.grid
            data = self.data * other.data
        else:
            grid = self.grid + other.grid
            data = self.resample(grid).data + other.resample(grid).data
        return GridData(grid, data)
    
    
    def __sub__(self, other):
        """Subtract two GridData objects by resampling to a common grid."""
        if self.grid == other.grid:
            grid = self.grid
            data = self.data - other.data
        else:
            grid = self.grid + other.grid
            data = self.resample(grid).data - other.resample(grid).data
        return GridData(grid, data)       
    
    
    def __mul__(self, other):
        """Multiply two GridData objects by resampling to a common grid."""
        if self.grid == other.grid:
            grid = self.grid
            data = self.data * other.data
        else:
            grid = self.grid + other.grid
            data = self.resample(grid).data * other.resample(grid).data
         
        return GridData(grid, data) 
    
    
    def squeeze_dims(self):
        """Squeezes grid dimensions for which the data is constant"""
        grid = self.grid
        data = self.data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)        
            if self.type == '3D':
                std_x = np.nanstd(data, axis=0, ddof=-1) 
                std_y = np.nanstd(data, axis=1, ddof=-1)
                if (np.nanmax(std_x) < 1e-5 and np.nanmax(std_y) < 1e-5):
                    std_z = np.nanstd(data, axis=2, ddof=-1)
                    if np.nanmax(std_z) < 1e-5:
                        data = np.nanmean(data)
                        grid = shdom.Grid(bounding_box=grid.bounding_box)
                    else:
                        data = np.nanmean(np.nanmean(data, axis=0), axis=0)
                        grid = shdom.Grid(bounding_box=grid.bounding_box, z=grid.z)
            elif self.type == '1D':
                std = np.nanstd(data, ddof=-1)
                if np.nanmax(std) < 1e-5: 
                    data = np.nanmean(data)
                    grid = shdom.Grid(bounding_box=grid.bounding_box)   
        return GridData(grid, np.nan_to_num(data))

        
    def resample(self, grid, method='linear'):
        """
        Resample data to a new Grid
        
        Parameters
        ----------
        grid: shdom.Grid
            The grid to which the data will be resampled
        method: str, default='linear'
            Options are: 'linear', 'nearest'
        
        Returns
        -------
        grid_data: shdom.GridData
            GridData sampled onto the new grid.
            
        Notes
        -----
        1. A 1D/3D grid sampled onto a Homogeneous grid yields the mean of the data.
        2. A 1D grid sampled onto a 3D grid tiles the data along the horizonal axes and reamples along the vertical axis.
        3. A 3D grid sampled onto 1D grid reamples along the vertical axis and averages along the horizontal axes
        """
        if self.grid == grid:
            return self   
        
        else:
            if self.type =='Homogeneous':
                data = np.full(shape=grid.shape, fill_value=self.data)
            
            elif self.type == '1D':
                if grid.type == 'Homogeneous':
                    data = np.mean(self.data)
                else:
                    if method == 'linear':
                        data = self._linear_interpolator1d(grid.z)
                    elif method == 'nearest':
                        data = self._nearest_interpolator1d(grid.z)
                    if grid.type == '3D':
                        data = np.tile(data[np.newaxis, np.newaxis, :], (grid.nx, grid.ny, 1))
                    
            elif self.type == '3D':
                if grid.type == 'Homogeneous':
                    data = np.mean(self.data)
                elif grid.type == '1D':
                    data = self._linear_interpolator1d(grid.z)
                    data = np.mean(np.mean(data, axis=0), axis=0)
                else:
                    if method == 'linear':
                        data = self._linear_interpolator3d(np.stack(np.meshgrid(grid.x, grid.y, grid.z, indexing='ij'), axis=-1))
                    elif method == 'nearest':
                        data = self._nearest_interpolator3d(np.stack(np.meshgrid(grid.x, grid.y, grid.z, indexing='ij'), axis=-1)) 
        return GridData(grid, data.astype(self.data.dtype))
    
    
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
    def max_value(self):
        return self.data.max()
    
    @property
    def min_value(self):
        return self.data.min()    
    
    @property
    def type(self):
        return self._type    


class BoundingBox(object):
    """ 
    A bounding box object.
    
    Parameters
    ----------
    xmin: float
         Minimum x (North) [km].
    ymin: float
         Minimum y (East) [km].
    zmin: float
         Minimum z (Up) [km].
    xmax: float
         Maximum x (North) [km].
    ymax: float
         Maximum y (East) [km].
    zmax: float
         Maximum z (Up) [km].
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
        xmin = min(self.xmin, other.xmin)
        ymin = min(self.ymin, other.ymin)
        zmin = min(self.zmin, other.zmin)
        xmax = max(self.xmax, other.xmax)
        ymax = max(self.ymax, other.ymax)
        zmax = max(self.zmax, other.zmax)
        return BoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)


from shdom.phase import *
from shdom.medium import *
from shdom.sensor import *
from shdom.rte_solver import *
from shdom.optimize import *
import shdom.generate as Generate


def save_forward_model(directory, medium, solver, measurements):
    """
    Save the forward model parameters for reconstruction.
    
    Parameters
    ----------
    directory: str
        Directory path where the forward modeling parameters are saved. 
        If the folder doesnt exist it will be created.
    medium: shdom.Medium object
        The atmospheric medium. This ground-truth medium will be used for comparisons.
    solver: shdom.RteSolver object
        The solver and the parameters used. This includes the scene parameters (such as solar and surface parameters)
        and the numerical parameters.
    measurements: shdom.Measurements
        Contains the camera used and the measurements acquired. 
        
    Notes
    -----
    The ground-truth medium is later used for evaulation of the recovery.
    """  
    if not os.path.isdir(directory):
        os.makedirs(directory)  
    measurements.save(os.path.join(directory, 'measurements'))
    medium.save(os.path.join(directory, 'ground_truth_medium'))
    solver.save_params(os.path.join(directory, 'solver_parameters'))   


def load_forward_model(directory):
    """
    Save the forward model parameters for reconstruction.
    
    Parameters
    ----------
    directory: str
        Directory path where the forward modeling parameters are saved. 
    
    Returns
    -------
    medium: shdom.Medium object
        The ground-truth atmospheric medium. 
    solver: shdom.RteSolver object
        The solver and the parameters used. This includes the scene parameters (such as solar and surface parameters)
        and the numerical parameters.
    measurements: shdom.Measurements
        Contains the sensor used to image the mediu and the radiance measurements. 
        
    Notes
    -----
    The ground-truth medium is used for evaulation of the recovery.
    """  
    # Load the ground truth medium for error analysis and ground-truth known phase and albedo
    medium_path = os.path.join(directory, 'ground_truth_medium')
    if os.path.exists(medium_path):
        medium = shdom.Medium()
        medium.load(path=medium_path)   
    else: 
        medium = None
        
    # Load shdom.Measurements object (sensor geometry and radiances)
    measurements = shdom.Measurements()
    measurements_path = os.path.join(directory, 'measurements')
    assert os.path.exists(measurements_path), 'No measurements file in directory: {}'.format(directory)
    measurements.load(path=measurements_path)
    
    # Load RteSolver according to numerical and scene parameters
    solver_path = os.path.join(directory, 'solver_parameters')
    if np.array(medium.wavelength).size == 1:
        solver = shdom.RteSolver()
    else:
        solver = shdom.RteSolverArray()
    if os.path.exists(solver_path):
        solver.load_params(path=os.path.join(directory, 'solver_parameters'))   
        
    return medium, solver, measurements

class SolarSpectrum(object):
    """
    Loads and interpolates the solar spectrum from ancillary_data/SpectralSolar_MODWehrli_1985_WMO.npz.
    Returns the solar spectral irradiance at specified monochromatic wavelengths. Database is valid for  
    wavelengths in the range 0.2 - 200.0 micrometers. Note that thermal emission becomes substantial beyond 3.0 micrometers
    and dominates at longer wavelengths.
    
    The database has high <.04 nm resolution in the visible range so monochromatic solar irradiances
    may differ substantially from band averaged quantities where the solarspectrum is more variable (at shorter wavelengths).
    
    Parameters
    ----------
    filename: str
        Directory path to the Solar Spectral Irradiance database.
        Default: 'ancillary_data/SpectralSolar_MODWehrli_1985_WMO.npz'
    """
    def __init__(self, filename='ancillary_data/SpectralSolar_MODWehrli_1985_WMO.npz'):
        
        self.filename = filename
        self._load()
        self._interpolate()
    
    def _load(self):
        """
        Read Spectral Irradiance Data from .npz file
        """
        dataset = np.load(self.filename)
        self.wavelengths = dataset['Wavelengths']
        self.Spectral_Irradiance = dataset['SolarSpectralIrradiance']
        self.description = dataset['description']
        self.source = dataset['source']
        self.wavelength_units = dataset['wavelength_units']
        self.Spectral_Irradiance_units = dataset['solarspectralirradiance_units']
    
    
    def _interpolate(self):
        """
        Interpolates the Spectrum so it can be retrieved at any wavelength
        """
        self.interpolater = interp1d(self.wavelengths, self.Spectral_Irradiance)
    
    def get_monochrome_solar_flux(self, wavelengths, Units = 'Micrometers'):
        """
        Returns the Solar Spectral Irradiance interpolated to the provided wavelengths.
        
        Parameters
        ----------
        wavelengths: array_like
            Wavelengths in units of micrometers. Must be in the range 0.2 - 200.0 micrometers.
            Note that thermal emission is significant beyond 3.0 micrometers and dominates at
            longer wavelengths.
            
        Returns
        -------
        Interpolated_Spectral_Irradiance: np.array 
            The Solar Flux interpolated to each of the input wavelengths.
        """
        
        if type(wavelengths) == list:
            wavelengths = np.array(wavelengths)
        
        if Units == 'Micrometers':
            wavelengths = wavelengths * 1000.0
        else:
            print('Wavelength units not supported')
        
        if np.any(wavelengths <= np.min(self.wavelengths)) or np.any(wavelengths >=np.max(self.wavelengths)):
            print('Wavelengths must be in Range {} - {} Micrometer'.format(np.round(np.min(self.wavelengths)/1000.0,6),np.max(self.wavelengths)/1000.0))
            
            return
        
        elif np.any((wavelengths > 3000.0) & (wavelengths <np.max(self.wavelengths))):
            print('Warning: Accurate Modelling of TOA radiances (> 3.0 micrometer) requires consideration of terrestrial thermal emission')
        
        Interpolated_Spectral_Irradiance = self.interpolater(wavelengths)
        
        return Interpolated_Spectral_Irradiance
