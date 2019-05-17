"""
Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer.
Here all the necessary functions used to solve the RT using SHDOM are wrapped.
"""

import core
import numpy as np
from enum import Enum
import warnings
import sys, os
import dill as pickle
from joblib import Parallel, delayed

class BoundaryCondition(Enum):
    open = 1       # open boundary conditions mean that exiting radiance is lost.
    periodic = 2   # periodic boundary conditions mean that exiting radiance returns from the opposite side.
    

class SolarSource(object):
    """ 
    Solar source object.
    
    Parameters
    ----------
    flux: float
        Top of medium solar flux on a horizontal surface (any units).
        For k-distribution this is a multiplier for the solar flux in the CKD file (i.e. normally should be 1).
    azimuth: float,
        Solar beam azimuthal angle (photon direction); specified in degrees but immediately converted to radians for use in code. 
        0 is beam going in positive X direction (North), 90 is positive Y (East).
    zenith: 
        Solar beam zenith angle in range (90,180]; Specified in degrees but immediately converted to the cosine of the angle (mu).
        This angle represents the direction of travel of the solar beam, so is forced to be negative although it can be specified positive.
    
    Notes
    -----
    Paramters: skyrad, units, wavenumber, are used for Thermal source and are not supported. 
    """
    def __init__(self, azimuth, zenith, flux=1.0):
        self.type = 'Solar'
        self.flux = flux
        self.azimuth = azimuth
        self.zenith = zenith
        self.units = 'R'
        self.skyrad = 0.0
        self.wavenumber = [10000, 10001]
        
    @property
    def zenith(self):
        return self._zenith
    
    @zenith.setter
    def zenith(self, val):
        assert 90.0 < val <= 180.0, 'Solar zenith:{} is not in range (90, 180] (photon direction in degrees)'.format(val)
        self._zenith = val

    @property 
    def info(self):
        """
        Print out all the source parameters.        
        """        
        return '{}, flux: {}, azimuth: {}deg, zenith: {}deg'.format(self.type, 
                                                                    self.flux, 
                                                                    self.azimuth, 
                                                                    self.zenith)  
        
class Surface(object):
    """ 
    An abstract sufrace class to be inhirted by different surface types.
    """
    def __init__(self):
        self.type = 'AbstractSurfaceClass'
    
    @property 
    def info(self):
        return '{}'.format(self.type)


class LambertianSurface(Surface):
    """
    A Lambertian surface, defined by a single parameter.
    
    Parameters
    ----------
    albedo: float, optional
        Bottom surface Lambertian albedo the range [0, 1] where 0 means completly absorbing.
        
    Notes
    -----
    Parameter: ground_temperature, is used for Thermal source and is not supported. 
    """
    def __init__(self, albedo):
        super(LambertianSurface, self).__init__()
        self.type = 'Lambertian'
        self.albedo = albedo
        self.ground_temperature = 298.15
    
    @property 
    def info(self):
        """
        Print out all the surface parameters.        
        """        
        return super(LambertianSurface, self).info + ', albedo: {}'.format(self.albedo)
    
        
class NumericalParameters(object):
    """
    This object bundles up together all the numerical parameters requiered by the RteSolver object.
    An instantitation will have all parameters intialized to their default values defined below.
    Below is the description for each parameter as taken from the SHDOM documentation (http://nit.colorado.edu/shdom/shdomdoc/)
    In capitals is the parameter name in the source fortran code.
    
    Properties
    ----------
    num_mu_bins(NMU): number of discrete ordinates covering -1 < mu < 1.
    num_phi_bins(NPHI): number of discrete ordinates covering 0 < phi < 2pi
    split_accuracy(SPLITACC): cell splitting accuracy; grid cells that have the adaptive splitting criterion above this value are split. 
                              This is an absolute measure, but cannot be easily associated with the resulting radiometric accuracy. 
                              Set to zero or negative for no adaptive cell splitting.
    deltam(DELTAM): True for delta-M scaling of medium and Nakajima and Tanaka method of computing radiances.
    spherical_harmonics_accuracy(SHACC): adaptive spherical harmonic truncation accuracy; the spherical harmonic source function series is truncated after 
                                         the terms are below this level. Truncation can still happens if SHACC=0 (for 0 source terms). 
                                         This is also an absolute measure, and is approximately the level of accuracy.
    acceleration_flag(ACCELFLAG): True to do the sequence acceleration. 
                                  An acceleration extrapolation of the source function is done every other iteration.
    solution_accuracy(SOLACC): solution accuracy - tolerance for solution criterion.
    max_total_mb(MAX_TOTAL_MB): approximate maximum memory to use (MB for 4 byte reals)
    adapt_grid_factor(ADAPT_GRID_FACTOR): ratio of total grid points to base grid points
    num_sh_term_factor(NUM_SH_TERM_FACTOR): ratio of average number of spherical harmonic terms to total possible (NLM)
    cell_to_point_ratio(CELL_TO_POINT_RATIO): ratio of number of grid cells to grid points
    high_order_radiance: True to keep the high order radiance field in memory.
    
    Notes
    -----
    deltam is a crucial parameter that should be set to True used with highly peaked (mie) phase function.
    """    
    def __init__(self,
                 num_mu_bins=8,
                 num_phi_bins=16,
                 split_accuracy=0.1,
                 deltam=True, 
                 spherical_harmonics_accuracy=0.01,
                 solution_accuracy=0.0001,
                 acceleration_flag=True,
                 max_total_mb=10000.0,
                 adapt_grid_factor=5,
                 num_sh_term_factor=5,
                 cell_to_point_ratio=1.5,
                 high_order_radiance=True):
        
        self.num_mu_bins = num_mu_bins
        self.num_phi_bins = num_phi_bins
        self.split_accuracy = split_accuracy
        self.deltam = deltam
        self.spherical_harmonics_accuracy = spherical_harmonics_accuracy
        self.solution_accuracy = solution_accuracy
        self.acceleration_flag = acceleration_flag
        self.max_total_mb = max_total_mb
        self.adapt_grid_factor = adapt_grid_factor
        self.num_sh_term_factor = num_sh_term_factor
        self.cell_to_point_ratio = cell_to_point_ratio
        self.high_order_radiance = high_order_radiance
        
    @property 
    def info(self):
        """
        Print out all the numerical parameters.        
        """
        info = 'Numerical Parameters: {}'.format(os.linesep)        
        for item in self.__dict__.iteritems():
            info += '   {}: {}{}'.format(item[0], item[1], os.linesep)
        return info     

class SceneParameters(object):
    """
    This object bundles up together all the scene related parameters requiered by the RteSolver object.
    An instantitation will have all parameters intialized to their default values defined below.
    
    Properties
    ----------
    wavelength: wavelength in microns for 'R' units; WAVELEN not needed for solar sources.
    surface: A Surface object.
    source: A Source object.
    boundary_conditions: a dictionary with BoundaryCondition.open/periodic for 'x' and 'y'.
    info: prints out all the properties as a string.
    
    Notes
    -----
    Currently supports LambertianSurface, SolarSource.
    """  
    def __init__(self, 
                 wavelength=0.672, 
                 surface=LambertianSurface(albedo=0.05), 
                 source=SolarSource(azimuth=0.0, zenith=180.0),
                 boundary_conditions={'x': BoundaryCondition.open,
                                      'y': BoundaryCondition.open}
                 ): 
        self.wavelength = wavelength
        self.surface = surface
        self.source = source
        self.boundary_conditions = boundary_conditions
    
    @property
    def info(self):
        return 'Scene Parameters: {}'.format(os.linesep) + \
               '   Wavelength: [{} micron]{}'.format(self.wavelength, os.linesep) + \
               '   Surface: [{}]{}'.format(self.surface.info, os.linesep) + \
               '   Source: [{}]{}'.format(self.source.info, os.linesep) + \
               '   Boundary Conditions: [x:{}, y:{}]{}'.format(self.boundary_conditions['x'].name, 
                                                               self.boundary_conditions['y'].name, 
                                                               os.linesep)
        


class ShdomPropertyArrays(object):
    """
    Shdom property array module.
    Contain the parameters that are in the original SHDOM_PROPERTY_ARRAYS module 
    in shdom90.f90 .
    
    Parameters
    ----------
    npx: Number of x grid points
    npy: Number of y grid points  
    npz: Number of z grid points 
    numphase: Number of phase function enteries in the table
    delx: Delta-x spacing
    dely: Delta-y spacing 
    xstart: Starting position of x coordinates
    ystart: Starting position of y coordinates 
    zlevels: Altitude grid points
    tempp: Temperatures at altitude grid points
    extinctp: Extinction grid
    albedop: Single scattering albedo grid
    legenp: legendre phase function table
    extdirp: optical depth integrated along sun rays
    iphasep: pointer to phase function table
    nzckd: Number of correlated-k distribution points
    zckd: correlated-k distribution grid
    gasabs: Gas absorption
    """
    def __init__(self):
        self.npx = None
        self.npy = None
        self.npz = None
        self.numphase = None
        self.delx = None
        self.dely = None
        self.xstart = None
        self.ystart = None
        self.zlevels = None
        self.tempp = None
        self.extinctp = None
        self.albedop = None
        self.legenp = None
        self.extdirp = None
        self.iphasep = None
        self.nzckd = None
        self.zckd = None
        self.gasabs = None
        
        
class RteSolver(object):
    """
    Radiative Trasnfer solver object. 
    This object contains the interface to SHDOM internal structures and methods.
    To solve the RTE:
      1. Set the solver parameters with RteSolver.set_parameters(scene_parameters, numerical_parameters)
      2. Attach the solver to a medium object RteSolver.init_medium(Medium)
      3. Run solution iterations with RteSolver.solve(maxiter)
    
    Notes
    -----
    k-distribution not supported.
    """
    
    def __init__(self, scene_params=None, numerical_params=None):
        
        self._type = 'Radiance'
        
        # Start mpi (if available).
        self._masterproc = core.start_mpi()

        # Link to the properties array module.
        self._pa = ShdomPropertyArrays()
        
        if scene_params:
            self.set_scene(scene_params)
        if numerical_params:
            self.set_numerics(numerical_params)
        
    
    def save_params(self, path):
        """
        Save RteSolver parameters from file.
    
        Parameters
        ----------
        path: str,
            Full path to file. 
        """
        param_dict = {'scene_params': self._scene_parameters,
                      'numerical_params': self._numerical_parameters}
        file = open(path,'w')
        file.write(pickle.dumps(param_dict, -1))
        file.close()    
        

    def load_params(self, path):
        """
        Load RteSolver parameters from file.

        Parameters
        ----------
        path: str,
            Full path to file. 
        """        
        file = open(path, 'r')
        data = file.read()
        file.close()
        params = pickle.loads(data)
        self.set_numerics(params['numerical_params'])
        self.set_scene(params['scene_params'])
        

    def set_scene(self, scene_params):
        """
        Set the scene related parameters: 
          wavelength, source, surface, boundary conditions, k-distribution 
    
        Parameters
        ----------
        scene_params : SceneParameters
             An object which encapsulate scene parameters such as solar and surface parameters.
    
        Notes
        -----
        k-distribution not supported. Only solar source supported.
        """  
        
        self._scene_parameters = scene_params
        
        # Wavelength
        self._wavelen = scene_params.wavelength
        
        # Source parameters
        if scene_params.source.type == 'Solar':
            self._srctype = 'S'
        else:
            raise NotImplementedError('Not implemented source type {}'.format(scene_params.source.type))
        self._solarflux = scene_params.source.flux
        self._solarmu = np.cos(np.deg2rad(scene_params.source.zenith))
        self._solaraz = np.deg2rad(scene_params.source.azimuth)
        self._skyrad = scene_params.source.skyrad
        self._units = scene_params.source.units
        self._waveno = scene_params.source.wavenumber
        
        # Surface parameters
        self._maxsfcpars = 4
        if scene_params.surface.type == 'Lambertian':
            self._sfctype = 'FL'  # Fixed Lambertian
            self._gndalbedo = scene_params.surface.albedo
            self._gndtemp = scene_params.surface.ground_temperature
            self._nxsfc = 0
            self._nysfc = 0
            self._delxsfc = 0
            self._delysfc = 0
            self._nsfcpar = 1
            self._sfcparms = []
            self._sfcgridparms = []  
        else:
            raise NotImplementedError('Only Lambertian surface is currently supported')
    
        # Boundary conditions
        self._bcflag = 0
        if scene_params.boundary_conditions['x'] == BoundaryCondition.open:
            self._bcflag += 1
        if scene_params.boundary_conditions['y'] == BoundaryCondition.open:
            self._bcflag += 2
        
        # Independent pixel not supported yet
        self._ipflag = 0
        
        # k-distribution not supported yet
        self._kdist = False
        self._ng = 1
        self._delg = np.ones(1, order='F')
        self._pa.nzckd = 0
        self._baseout = False


    def set_numerics(self, numerical_params):
        """
        Set the numerical parameters of the SHDOM forward solver.
    
        Parameters
        ----------
        numerical_params : NumericalParameters
            An object which encapsulate numerical parameters such as number of azimuthal and zenith bins. 
    
        Notes
        -----
        Curently not delta-m is not suppoted.
        """ 
        self._numerical_parameters = numerical_params
        self._nmu                 = max(2, 2 * int((numerical_params.num_mu_bins + 1) / 2))
        self._nphi                = max(1, numerical_params.num_phi_bins)
        self._deltam              = numerical_params.deltam
        self._max_total_mb        = numerical_params.max_total_mb
        self._adapt_grid_factor   = numerical_params.adapt_grid_factor
        self._accelflag           = numerical_params.acceleration_flag
        self._num_sh_term_factor  = numerical_params.num_sh_term_factor
        self._cell_to_point_ratio = numerical_params.cell_to_point_ratio
        self._splitacc            = numerical_params.split_accuracy
        self._shacc               = numerical_params.spherical_harmonics_accuracy
        self._solacc              = numerical_params.solution_accuracy
        self._highorderrad        = numerical_params.high_order_radiance
    
    
    def init_medium(self, medium):
        """
        Initilize SHDOM internal grid structures according to the Medium.
    
        Parameters
        ----------
        medium: shdom.Medium
            Initilize the RTE solver to a Medium object.

        """
        self._npart = medium.num_scatterers
        self.set_grid(medium.grid)
        self.set_extinction(medium)
        self.set_albedo(medium)
        self.set_phase(medium)
        
        # Temperature is used for thermal radiation. Not supported yet.
        self._pa.tempp = np.zeros(shape=(self._maxpg,), dtype=np.float32)
        
        # Zero itertions so far
        self._iters = 0
        self._solcrit = 1.0
        self._inradflag = False 
        self._oldnpts = 0


    def set_phase(self, medium):
        """
        set the phase function internal SHDOM parameters
        
        Parameters
        ----------
        medium: shdom.OpticalMedium
            an OpticalMedium object contains legenp,iphasep properties
            
        """
        self._pa.iphasep = medium.iphasep.astype(np.int32)
        self._pa.numphase = medium.numphase
        
        # Determine the number of legendre coefficient for a given angular resolution
        if self._deltam:
            self._nleg = self._mm+1
        else:
            self._nleg = self._mm
        self._nleg = self._maxleg = max(medium.maxleg, self._nleg)
        
        # Legenp is without the zero order term which is 1.0 for normalized phase function
        self._pa.legenp = medium.get_legenp(self._nleg).astype(np.float32)
        self._maxasym = medium.maxasym
        self._maxpgl = medium.grid.num_points * medium.maxleg         

        if medium.numphase > 0:
            self._maxigl = medium.numphase*(medium.maxleg + 1)
        else:
            self._maxigl = self._maxig*(medium.maxleg + 1)
    
        self._iphase = np.empty(shape=(self._maxig, self._npart), dtype=np.int32, order='F')
        self._legen = np.empty(shape=(self._maxigl,), dtype=np.float32, order='F')        
        self._ylmsun = np.empty(shape=(self._nlm, ), dtype=np.float32, order='F') 
                   

    def set_albedo(self, medium):
        """
        set the single scattering albedo
        
        Parameters
        ----------
        medium: shdom.OpticalMedium
            an OpticalMedium object contains albedop property
        """
        self._pa.albedop = medium.albedop.astype(np.float32)
        
        
    def set_extinction(self, medium):
        """
        set the optical extinction.
        
        Parameters
        ----------
        medium: shdom.OpticalMedium
            an OpticalMedium object contains extinctp property
        """
        self._pa.extinctp = medium.extinctp.astype(np.float32)
    
    
    def set_grid(self, grid):
        """
        set the base grid for SHDOM.
        
        Parameters
        ----------
        grid: shdom.Grid
            The grid.
        """
        
        def ibits(val, bit, ret_val):
            if val & 2**bit:
                return ret_val
            return 0         

        # Set shdom property array
        self._pa.npx = grid.nx
        self._pa.npy = grid.ny    
        self._pa.npz = grid.nz
        self._pa.xstart = grid.xmin
        self._pa.ystart = grid.ymin
        self._pa.delx = grid.dx
        self._pa.dely = grid.dy
        self._pa.zlevels = grid.z        
        
        # Initialize shdom internal grid sizes to property array grid
        self._nx = grid.nx
        self._ny = grid.ny
        self._nz = grid.nz
        self._maxpg  = grid.num_points
        self._maxnz = grid.nz
        
        # Set the full domain grid sizes (variables end in t)
        self._nxt, self._nyt, self._npxt, self._npyt = self._nx, self._ny, self._pa.npx, self._pa.npy
        self._nx, self._ny, self._nz = max(1, self._nx), max(1, self._ny), max(2, self._nz)
        
        # gridtype = 'P': Z levels taken from property file
        self._gridtype = 'P'        

        # Set up base grid point actual size (NX1xNY1xNZ)
        self._nx1, self._ny1 = self._nx+1, self._ny+1
        if self._bcflag & 5:
            self._nx1 -= 1
        if self._bcflag & 7:
            self._ny1 -= 1
        self._nbpts = self._nx1 * self._ny1 * self._nz 
        
        # Calculate the number of base grid cells depending on the BCFLAG
        self._nbcells = (self._nz - 1) * (self._nx + ibits(self._bcflag, 0, 1) - \
                        ibits(self._bcflag, 2, 1)) * (self._ny + ibits(self._bcflag, 1, 1) - ibits(self._bcflag, 3, 1))
        
        self._xgrid, self._ygrid, self._zgrid = core.new_grids(
            bcflag=self._bcflag,
            gridtype=self._gridtype,
            npx=self._pa.npx,
            npy=self._pa.npy,
            nx=self._nx,
            ny=self._ny,
            nz=self._nz,
            xstart=self._pa.xstart,
            ystart=self._pa.ystart,
            delxp=self._pa.delx,
            delyp=self._pa.dely,
            zlevels=self._pa.zlevels
        )        
        
        self.init_memory()
        
        self._npts, self._ncells, self._gridpos, self._gridptr, self._neighptr, \
            self._treeptr, self._cellflags = core.init_cell_structure(
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                nx=self._nx,
                ny=self._ny,
                nz=self._nz,
                nx1=self._nx1,
                ny1=self._ny1,
                xgrid=self._xgrid,
                ygrid=self._ygrid,
                zgrid=self._zgrid,
                gridpos=self._gridpos,
                gridptr=self._gridptr,
                neighptr=self._neighptr,
                treeptr=self._treeptr,
                cellflags=self._cellflags
            )
        self._nbcells = self._ncells
        

    def init_memory(self):
        """A utility function to initialize internal memory structures and parameters."""
        
        # Make ml and mm from nmu and nphi
        # ML is the maximum meridional mode, MM is the maximum azimuthal mode,
        # and NCS is the azimuthal mode flag (|NCS|=1 for cosine only, |NCS|=2 for 
        # sines and cosines).
        # nphi0max: The maximum number of azimuth angles actually used;
        # for NCS=1 (cosine modes only) NPHI0=INT((NPHI+2)/2),
        # otherwise NPHI0=NPHI.
        self._ncs = 2
        self._nstokes = 1
        self._nstleg = 1
        self._ml = self._nmu - 1
        self._mm = max(0, int(self._nphi / 2) - 1)
        self._nlm = (2 * self._mm + 1) * (self._ml + 1) - self._mm * (self._mm + 1)
    
        self._nphi0max = self._nphi
        self._memword = self._nmu*(2 + 2 * self._nphi + 2 * self._nlm + 2 * 33 *32)
    
        # Guess maximum number of grid points, cells, SH vector size needed
        # but don't let MAX_TOTAL_MB be exceeded
        if self._max_total_mb*1024**2 > 1.75*sys.maxsize:
            self._max_total_mb > 1.75*sys.maxsize/1024**2 
            print 'MAX_TOTAL_MB reduced to fit memory model: ', self._max_total_mb
    
        if self._splitacc < 0.0:
            self._adapt_grid_factor = 1.0
    
        self._big_arrays = 3
        if not self._accelflag:
            self._big_arrays = 2        
    
        wantmem = self._adapt_grid_factor * self._nbpts * (
            28 + 16.5 * self._cell_to_point_ratio + self._nphi0max + self._num_sh_term_factor * self._nlm * self._big_arrays)                
    
        REDUCE = min(1.0, ((self._max_total_mb * 1024**2) / 4 - self._memword) / wantmem)
        self._adapt_grid_factor *= REDUCE
        assert self._adapt_grid_factor > 1.0, 'max_total_mb memory limit exceeded with just base grid.'
    
        if REDUCE < 1.0:
            print 'adapt_grid_factor reduced to ', self._adapt_grid_factor
    
        wantmem *= REDUCE
        assert wantmem < sys.maxsize, 'number of words of memory exceeds max integer size: %d'% wantmem
        self._maxig = int(self._adapt_grid_factor*self._nbpts)
        self._maxic = int(self._cell_to_point_ratio*self._maxig)
        self._maxiv = int(self._num_sh_term_factor*self._nlm*self._maxig)
        self._maxido = self._maxig*self._nphi0max
    
        assert 4.0*(self._maxiv+self._maxig) < sys.maxsize, 'size of big sh arrays (maxiv) probably exceeds max integer number of bytes: %d' % self._maxiv
        assert 4.0*8.0*self._maxic <= sys.maxsize, 'size of gridptr array (8*maxic) probably exceeds max integer number of bytes: %d' % 8*self._maxic
    
        self._maxnbc = int(self._maxig*3/self._nz)
        self._maxsfcpars = 4
        if self._sfctype.endswith('L'):
            self._maxbcrad = 2*self._maxnbc
        else:
            self._maxbcrad = int((2+self._nmu*self._nphi0max/2)*self._maxnbc)
    
        self._sfcgridparms = np.empty(self._maxsfcpars*self._maxnbc,dtype=np.float32)               
        self._bcptr = np.empty(shape=(self._maxnbc, 2), dtype=np.int32, order='F')
        self._bcrad = np.empty(shape=(self._maxbcrad,), dtype=np.float32, order='F')
    
        # Array allocation
        self._total_ext = np.empty(shape=(self._maxig,), dtype=np.float32, order='F')
        self._extinct = np.empty(shape=(self._maxig, self._npart), dtype=np.float32, order='F')
        self._albedo = np.empty(shape=(self._maxig, self._npart), dtype=np.float32, order='F')
        self._mu = np.empty(shape=(self._nmu,), dtype=np.float32, order='F')
        self._wtdo = np.empty(shape=(self._nmu*self._nphi,), dtype=np.float32, order='F')
        self._phi = np.empty(shape=(self._nmu*self._nphi,), dtype=np.float32, order='F')
        self._phi0 = np.empty(shape=(self._nmu,), dtype=np.int32, order='F')
        self._temp = np.empty(shape=(self._maxig,), dtype=np.float32, order='F')
        self._planck = np.empty(shape=(self._maxig, self._npart), dtype=np.float32, order='F')
        self._gridptr = np.empty(shape=(8, self._maxic), dtype=np.int32, order='F')
        self._neighptr = np.empty(shape=(6, self._maxic), dtype=np.int32, order='F')
        self._treeptr = np.empty(shape=(2, self._maxic), dtype=np.int32, order='F')
        self._cellflags = np.empty(shape=(self._maxic,), dtype=np.int16, order='F')
        self._gridpos = np.empty(shape=(3, self._maxig), dtype=np.float32, order='F')
        self._rshptr = np.empty(shape=(self._maxig+2,), dtype=np.int32, order='F')
        self._shptr = np.empty(shape=(self._maxig+1,), dtype=np.int32, order='F')
        self._oshptr = np.empty(shape=(self._maxig+1,), dtype=np.int32, order='F')
        self._radiance = np.empty(shape=(self._maxiv+self._maxig,), dtype=np.float32, order='F')
        self._source = np.empty(shape=(self._maxiv,), dtype=np.float32, order='F')
        self._delsource = np.empty(shape=(self._maxiv,), dtype=np.float32, order='F')
        self._fluxes = np.empty(shape=(2, self._maxig,), dtype=np.float32, order='F')
        self._fluxes1 = np.empty(shape=(4, self._maxig,), dtype=np.float32, order='F')
        self._dirflux = np.empty(shape=(self._maxig,), dtype=np.float32, order='F')
        self._pa.extdirp = np.empty(shape=(self._maxpg,), dtype=np.float32, order='F')
        self._work = np.empty(shape=(self._maxido,), dtype=np.float32, order='F')
        self._work1 = np.empty(shape=(8*self._maxig,), dtype=np.int32, order='F')
        self._work2 = np.empty(shape=(self._maxig), dtype=np.float32, order='F')        
           
    def solve(self, maxiter, verbose=True):
        """
        Main solver routine.

        Performs the SHDOM solution procedure.

        Parameters
        ----------
        maxiter: integer
            Maximum number of iterations for the iterative solution.
        verbose: verbosity
            True will output solution iteration information into stdout.
            
        Returns
        -------
        None
        """
        self._nang, self._nphi0, self._mu, self._phi, self._wtdo, self._sfcgridparms, self._solcrit, \
            iters, self._temp, self._planck, self._extinct, self._albedo, self._legen, self._iphase, \
            self._ntoppts, self._nbotpts, self._bcptr, self._bcrad, self._npts, self._gridpos, self._ncells, self._gridptr, \
            self._neighptr, self._treeptr, self._cellflags, self._rshptr, self._shptr, self._oshptr, self._source, \
            self._delsource, self._radiance, self._fluxes, self._dirflux, self._ylmsun, self._oldnpts, self._total_ext = \
            core.solve_rte(
                npart=self._npart,
                extinctp=self._pa.extinctp,
                albedop=self._pa.albedop,
                planck=self._planck,
                iphase=self._iphase,
                iphasep=self._pa.iphasep,                
                nstokes=self._nstokes,
                nstleg=self._nstleg,
                npx=self._pa.npx,
                npy=self._pa.npy,
                npz=self._pa.npz,
                delx=self._pa.delx,
                dely=self._pa.dely,                
                xstart=self._pa.xstart,
                ystart=self._pa.ystart,
                zlevels=self._pa.zlevels,             
                tempp=self._pa.tempp,          
                legenp=self._pa.legenp,
                extdirp=self._pa.extdirp,                            
                nzckd=self._pa.nzckd,
                zckd=self._pa.zckd,
                gasabs=self._pa.gasabs,                
                solcrit=self._solcrit,
                nx=self._nx,
                ny=self._ny,
                nx1=self._nx1,
                ny1=self._ny1,
                nz=self._nz,
                ml=self._ml,
                mm=self._mm,
                ncs=self._ncs,
                nlm=self._nlm,
                nmu=self._nmu,
                nphi=self._nphi,
                numphase=self._pa.numphase,
                mu=self._mu,
                phi=self._phi,
                wtdo=self._wtdo,
                inradflag=self._inradflag,
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                deltam=self._deltam,
                srctype=self._srctype,
                highorderrad=self._highorderrad,
                solarflux=self._solarflux,
                solarmu=self._solarmu,
                solaraz=self._solaraz,
                skyrad=self._skyrad,
                sfctype=self._sfctype,
                gndtemp=self._gndtemp,
                gndalbedo=self._gndalbedo,
                nxsfc=self._nxsfc,
                nysfc=self._nysfc,
                delxsfc=self._delxsfc,
                delysfc=self._delysfc,
                nsfcpar=self._nsfcpar,
                sfcparms=self._sfcparms,
                sfcgridparms=self._sfcgridparms,
                units=self._units,
                waveno=self._waveno,
                wavelen=self._wavelen,
                accelflag=self._accelflag,
                solacc=self._solacc,
                maxiter=maxiter,
                splitacc=self._splitacc,
                shacc=self._shacc,
                xgrid=self._xgrid,
                ygrid=self._ygrid,
                zgrid=self._zgrid,
                temp=self._temp,
                maxbcrad=self._maxbcrad,
                bcptr=self._bcptr,
                bcrad=self._bcrad,
                npts=self._npts,
                gridpos=self._gridpos,
                ncells=self._ncells,
                gridptr=self._gridptr,
                neighptr=self._neighptr,
                treeptr=self._treeptr,
                cellflags=self._cellflags,
                rshptr=self._rshptr,
                shptr=self._shptr,
                oshptr=self._oshptr,
                work=self._work,
                work1=self._work1,
                work2=self._work2,
                source=self._source,
                delsource=self._delsource,
                radiance=self._radiance,
                fluxes=self._fluxes,
                dirflux=self._dirflux,
                nleg=self._nleg,
                maxiv=self._maxiv,
                maxic=self._maxic,
                maxig=self._maxig,
                maxido=self._maxido,
                verbose=verbose,
                oldnpts=self._oldnpts,
                ylmsun=self._ylmsun
            )
        self._iters += iters
        self._inradflag = True
    
    @property
    def type(self):
        return self._type
    
    @property
    def num_iterations(self):
        return self._iters
    
    @property
    def info(self):
        return self._scene_parameters.info + os.linesep + self._numerical_parameters.info
    
    
    
class RteSolverPolarized(RteSolver):
    """
    Polarized Radiative Trasnfer solver object. 
    This object contains the interface to SHDOM internal structures and methods.
    To solve the RTE:
      1. Set the solver parameters with RteSolverPolarized.set_parameters(scene_parameters, numerical_parameters)
      2. Attach the solver to a medium object RteSolverPolarized.init_medium(Medium)
      3. Run solution iterations with RteSolverPolarized.solve(maxiter)
    
    Parameters
    ----------
    num_stokes: int
        The number of stokes for which to solve the RTE can be 1, 3, or 4. 
        num_stokes=1 means unpolarized.
        num_stokes=3 means linear polarization.
        num_stokes=4 means full polarization.
        
    Notes
    -----
    k-distribution not supported.
    """
    
    def __init__(self, num_stokes, scene_params=None, numerical_params=None):
        super(RteSolverPolarized, self).__init__(scene_params, numerical_params)
        self._type = 'Polarization'
        assert num_stokes in [1, 3, 4], 'num_stokes should be {1, 3, 4}'
        self._nstokes = num_stokes
        

    def set_phase(self, phase):
        """
        set the phase function internal SHDOM parameters
        
        Parameters
        ----------
        phase: shdom.Phase
          TabulatedPhase or GridPhase object.
            
        """
        self._pa.iphasep = phase.iphasep.ravel()
        self._pa.numphase = phase.numphase   
        
        # Determine the number of legendre coefficient for a given angular resolution
        if self._deltam:
            self._nleg = self._mm+1
        else:
            self._nleg = self._mm
        self._nleg = self._maxleg = max(phase.maxleg, self._nleg)

        # Legenp is without the zero order term which is 1.0 for normalized phase function
        self._pa.legenp = phase.get_legenp(self._nleg)

        self._maxasym = phase.maxasym
        self._maxpgl = phase.grid.num_points * phase.maxleg         

        if phase.numphase > 0:
            self._maxigl = phase.numphase*(phase.maxleg + 1)
        else:
            self._maxigl = self._maxig*(phase.maxleg + 1)
    
        self._iphase = np.empty(shape=(self._maxig, self._npart), dtype=np.int32, order='F')
        self._legen = np.empty(shape=(phase.nstleg, self._maxigl,), dtype=np.float32, order='F')        

        if self._nstokes == 1:
            self._nstleg = 1
        else:
            self._nstleg = phase.nstleg
            
        self._ylmsun = np.empty(shape=(self._nstleg, self._nlm), dtype=np.float32, order='F') 
           


    def init_memory(self):
        """A utility function to initialize internal memory structures and parameters."""
        
        # Make ml and mm from nmu and nphi
        # ML is the maximum meridional mode, MM is the maximum azimuthal mode,
        # nphi0max: The maximum number of azimuth angles actually used
        self._ncs = 2
        self._ml = self._nmu - 1
        self._mm = max(0, int(self._nphi / 2) - 1)
        self._nlm = (2 * self._mm + 1) * (self._ml + 1) - self._mm * (self._mm + 1)
    
        self._nphi0max = self._nphi
        self._memword = self._nmu*(2 + 2 * self._nphi + 2 * self._nlm + 2 * 33 *32)
    
        # Guess maximum number of grid points, cells, SH vector size needed
        # but don't let MAX_TOTAL_MB be exceeded
        if self._max_total_mb*1024**2 > 1.75*sys.maxsize:
            self._max_total_mb > 1.75*sys.maxsize/1024**2 
            print 'MAX_TOTAL_MB reduced to fit memory model: ', self._max_total_mb
    
        if self._splitacc < 0.0:
            self._adapt_grid_factor = 1.0
    
        self._big_arrays = 3
        if not self._accelflag:
            self._big_arrays = 2        
    
        wantmem = self._adapt_grid_factor * self._nbpts * (
            28 + 16.5 * self._cell_to_point_ratio + \
            self._nphi0max*self._nstokes + self._num_sh_term_factor*self._nstokes*self._nlm*self._big_arrays)                
    
        REDUCE = min(1.0, ((self._max_total_mb * 1024**2) / 4 - self._memword) / wantmem)
        self._adapt_grid_factor *= REDUCE
        assert self._adapt_grid_factor > 1.0, 'max_total_mb memory limit exceeded with just base grid.'
    
        if REDUCE < 1.0:
            print 'adapt_grid_factor reduced to ', self._adapt_grid_factor
    
        wantmem *= REDUCE
        assert wantmem < sys.maxsize, 'number of words of memory exceeds max integer size: %d'% wantmem
        self._maxig = int(self._adapt_grid_factor*self._nbpts)
        self._maxic = int(self._cell_to_point_ratio*self._maxig)
        self._maxiv = int(self._num_sh_term_factor*self._nlm*self._maxig)
        self._maxido = self._maxig*self._nphi0max
    
        assert 4.0*(self._maxiv+self._maxig)*self._nstokes < sys.maxsize, 'size of big sh arrays (maxiv) probably exceeds max integer number of bytes: %d' % self._maxiv
        assert 4.0*8.0*self._maxic <= sys.maxsize, 'size of gridptr array (8*maxic) probably exceeds max integer number of bytes: %d' % 8*self._maxic
    
        self._maxnbc = int(self._maxig*3/self._nz)
        self._maxsfcpars = 4
        if self._sfctype.endswith('L'):
            self._maxbcrad = 2*self._maxnbc
        else:
            self._maxbcrad = int((2+self._nmu*self._nphi0max/2)*self._maxnbc)
    
        self._sfcgridparms = np.empty(self._maxsfcpars*self._maxnbc,dtype=np.float32)               
        self._bcptr = np.empty(shape=(self._maxnbc, 2), dtype=np.int32, order='F')
        self._bcrad = np.empty(shape=(self._nstokes,self._maxbcrad), dtype=np.float32, order='F')

        # Array allocation
        self._extinct = np.empty(shape=(self._maxig, self._npart), dtype=np.float32, order='F')
        self._albedo = np.empty(shape=(self._maxig, self._npart), dtype=np.float32, order='F')
        self._mu = np.empty(shape=(self._nmu,), dtype=np.float32, order='F')
        self._wtdo = np.empty(shape=(self._nmu*self._nphi,), dtype=np.float32, order='F')
        self._phi = np.empty(shape=(self._nmu*self._nphi,), dtype=np.float32, order='F')
        self._phi0 = np.empty(shape=(self._nmu,), dtype=np.int32, order='F')
        self._temp = np.empty(shape=(self._maxig,), dtype=np.float32, order='F')
        self._planck = np.empty(shape=(self._maxig, self._npart), dtype=np.float32, order='F')
        self._gridptr = np.empty(shape=(8, self._maxic), dtype=np.int32, order='F')
        self._neighptr = np.empty(shape=(6, self._maxic), dtype=np.int32, order='F')
        self._treeptr = np.empty(shape=(2, self._maxic), dtype=np.int32, order='F')
        self._cellflags = np.empty(shape=(self._maxic,), dtype=np.int16, order='F')
        self._gridpos = np.empty(shape=(3, self._maxig), dtype=np.float32, order='F')
        self._rshptr = np.empty(shape=(self._maxig+2,), dtype=np.int32, order='F')
        self._shptr = np.empty(shape=(self._maxig+1,), dtype=np.int32, order='F')
        self._oshptr = np.empty(shape=(self._maxig+1,), dtype=np.int32, order='F')
        self._radiance = np.empty(shape=(self._nstokes, self._maxiv+self._maxig), dtype=np.float32, order='F')
        self._source = np.empty(shape=(self._nstokes, self._maxiv), dtype=np.float32, order='F')
        self._delsource = np.empty(shape=(self._nstokes, self._maxiv), dtype=np.float32, order='F')
        self._fluxes = np.empty(shape=(2, self._maxig,), dtype=np.float32, order='F')
        self._fluxes1 = np.empty(shape=(4, self._maxig,), dtype=np.float32, order='F')
        self._dirflux = np.empty(shape=(self._maxig,), dtype=np.float32, order='F')
        self._pa.extdirp = np.empty(shape=(self._maxpg,), dtype=np.float32, order='F')
        self._work = np.empty(shape=(self._nstokes*self._maxido,), dtype=np.float32, order='F')
        self._work1 = np.empty(shape=(8*self._maxig,), dtype=np.int32, order='F')
        self._work2 = np.empty(shape=(self._nstokes*self._maxig), dtype=np.float32, order='F')        
        

    
class RteSolverArray(object):
    """
    An RteSolverArray object encapsulate several solvers e.g. for multiple spectral imaging
    
    Parameters
    ----------
    solver_list: list, optional
        A list of Sensor objects
    """
    def __init__(self, solver_list=None):
        self._num_solvers = 0
        self._solver_list = []
        self._names = []
        self._solver_type = None
        if solver_list:
            for solver in solver_list:
                self.add_solver(solver)
        
    def __getitem__(self, val):
        return self.solver_list[val]
    
    def add_solver(self, rte_solver, name=None):
        """
        Add a rte_solver to the RteSolverArray
        
        Parameters
        ----------
        rte_solver: RteSolver object
            A RteSolver object to add to the RteSolverArray
        name: str, optional
            An ID for the solver. 
        """
        
        if self.solver_type is None:
            self._solver_type = rte_solver.type
        else:
            assert self.solver_type == rte_solver.type, \
                   '[add_solver] Assert: RteSolverArray is of type {} and new solver is of type {}'.format(self.solver_type, rte_solver.type)
            
        self._solver_list.append(rte_solver)
        self._num_solvers += 1
        if name is None:
            self._names.append('Solver{:d}'.format(self.num_solvers))
        else:
            self._names.append(name)


    def solve(self, maxiter, verbose=True):
        """
        Parallel solving of all solvers.

        Parameters
        ----------
        maxiter: integer
            Maximum number of iterations for the iterative solution.
        verbose: verbosity
            True will output solution iteration information into stdout.
            
        Returns
        -------
        None
        """
        output_arguments = \
            Parallel(n_jobs=self.num_solvers, backend="threading")(
                delayed(core.solve_rte, check_pickle=False)(
                    nstleg=rte_solver._nstleg,
                    nstokes=rte_solver._nstokes,
                    npx=rte_solver._pa.npx,
                    npy=rte_solver._pa.npy,
                    npz=rte_solver._pa.npz,
                    delx=rte_solver._pa.delx,
                    dely=rte_solver._pa.dely,                
                    xstart=rte_solver._pa.xstart,
                    ystart=rte_solver._pa.ystart,
                    zlevels=rte_solver._pa.zlevels,             
                    tempp=rte_solver._pa.tempp,
                    extinctp=rte_solver._pa.extinctp,
                    albedop=rte_solver._pa.albedop,
                    legenp=rte_solver._pa.legenp,
                    extdirp=rte_solver._pa.extdirp,
                    iphasep=rte_solver._pa.iphasep,
                    nzckd=rte_solver._pa.nzckd,
                    zckd=rte_solver._pa.zckd,
                    gasabs=rte_solver._pa.gasabs,
                    solcrit=rte_solver._solcrit,
                    nx=rte_solver._nx,
                    ny=rte_solver._ny,
                    nx1=rte_solver._nx1,
                    ny1=rte_solver._ny1,
                    nz=rte_solver._nz,
                    ml=rte_solver._ml,
                    mm=rte_solver._mm,
                    ncs=rte_solver._ncs,
                    nlm=rte_solver._nlm,
                    nmu=rte_solver._nmu,
                    nphi=rte_solver._nphi,
                    numphase=rte_solver._pa.numphase,
                    mu=rte_solver._mu,
                    phi=rte_solver._phi,
                    wtdo=rte_solver._wtdo,
                    inradflag=rte_solver._inradflag,
                    bcflag=rte_solver._bcflag,
                    ipflag=rte_solver._ipflag,
                    deltam=rte_solver._deltam,
                    srctype=rte_solver._srctype,
                    highorderrad=rte_solver._highorderrad,
                    solarflux=rte_solver._solarflux,
                    solarmu=rte_solver._solarmu,
                    solaraz=rte_solver._solaraz,
                    skyrad=rte_solver._skyrad,
                    sfctype=rte_solver._sfctype,
                    gndtemp=rte_solver._gndtemp,
                    gndalbedo=rte_solver._gndalbedo,
                    nxsfc=rte_solver._nxsfc,
                    nysfc=rte_solver._nysfc,
                    delxsfc=rte_solver._delxsfc,
                    delysfc=rte_solver._delysfc,
                    nsfcpar=rte_solver._nsfcpar,
                    sfcparms=rte_solver._sfcparms,
                    sfcgridparms=rte_solver._sfcgridparms,
                    units=rte_solver._units,
                    waveno=rte_solver._waveno,
                    wavelen=rte_solver._wavelen,
                    accelflag=rte_solver._accelflag,
                    solacc=rte_solver._solacc,
                    maxiter=maxiter,
                    splitacc=rte_solver._splitacc,
                    shacc=rte_solver._shacc,
                    xgrid=rte_solver._xgrid,
                    ygrid=rte_solver._ygrid,
                    zgrid=rte_solver._zgrid,
                    temp=rte_solver._temp,
                    planck=rte_solver._planck,
                    iphase=rte_solver._iphase,
                    maxbcrad=rte_solver._maxbcrad,
                    bcptr=rte_solver._bcptr,
                    bcrad=rte_solver._bcrad,            
                    npts=rte_solver._npts,
                    gridpos=rte_solver._gridpos,
                    ncells=rte_solver._ncells,
                    gridptr=rte_solver._gridptr,
                    neighptr=rte_solver._neighptr,
                    treeptr=rte_solver._treeptr,
                    cellflags=rte_solver._cellflags,
                    rshptr=rte_solver._rshptr,
                    shptr=rte_solver._shptr,
                    oshptr=rte_solver._oshptr,
                    work=rte_solver._work,
                    work1=rte_solver._work1,
                    work2=rte_solver._work2,
                    source=rte_solver._source,
                    delsource=rte_solver._delsource,
                    radiance=rte_solver._radiance,
                    fluxes=rte_solver._fluxes,
                    dirflux=rte_solver._dirflux,
                    nleg=rte_solver._nleg,
                    maxiv=rte_solver._maxiv,
                    maxic=rte_solver._maxic,
                    maxig=rte_solver._maxig,
                    maxido=rte_solver._maxido,
                    verbose=verbose,
                    oldnpts=rte_solver._oldnpts,
                    ylmsun=rte_solver._ylmsun
                    ) for rte_solver in self.solver_list)
     
        # Update solvers internal structures
        for i, solver in enumerate(self.solver_list):
            solver._nang, solver._nphi0, solver._mu, solver._phi, solver._wtdo, solver._sfcgridparms, solver._solcrit, \
                iters, solver._temp, solver._planck, solver._extinct, solver._albedo, solver._legen, solver._iphase, \
                solver._ntoppts, solver._nbotpts, solver._bcptr, solver._bcrad, solver._npts, solver._gridpos, solver._ncells, solver._gridptr, \
                solver._neighptr, solver._treeptr, solver._cellflags, solver._rshptr, solver._shptr, solver._oshptr,\
                solver._source, solver._delsource, solver._radiance, solver._fluxes, solver._dirflux, solver._ylmsun, solver._oldnpts = \
                output_arguments[i]
            solver._iters += iters
            solver._inradflag = True 
            
            
    @property
    def solver_list(self):
        return self._solver_list
    
    @property
    def num_solvers(self):
        return self._num_solvers
    
    @property
    def names(self):
        return self._names  
    
    @property
    def solver_type(self):
        return self._solver_type