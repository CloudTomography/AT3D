"""
Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer.
Here all the necessary functions used to solve the RT using SHDOM are wrapped.
"""
import numpy as np
from enum import Enum
import sys, os, copy
import dill as pickle
from joblib import Parallel, delayed
import shdom 
from shdom import core, float_round


class BC(Enum):
    """
    Two types of boundary conditions:
      1. open - exiting radiance is lost.
      2. periodic - exiting radiance returns from the opposite side.
    """
    open = 1 
    periodic = 2 
    

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
        return '{}, flux: {}, azimuth: {}deg, zenith: {}deg'.format(self.type, self.flux, self.azimuth, self.zenith)


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
    This object bundles up together all the numerical parameters required by the RteSolver object.
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
        for item in self.__dict__.items():
            info += '   {}: {}{}'.format(item[0], item[1], os.linesep)
        return info     


class SceneParameters(object):
    """
    This object bundles up together all the scene related parameters required by the RteSolver object.
    An instantitation will have all parameters intialized to their default values defined below.
    
    Properties
    ----------
    wavelength: wavelength in microns for 'R' units; WAVELEN not needed for solar sources.
    surface: A Surface object.
    source: A Source object.
    boundary_conditions: a dictionary with BC.open/periodic for 'x' and 'y'.
    info: prints out all the properties as a string.
    
    Notes
    -----
    Currently supports LambertianSurface, SolarSource.
    """  
    def __init__(self, 
                 wavelength=0.672, 
                 surface=LambertianSurface(albedo=0.05), 
                 source=SolarSource(azimuth=0.0, zenith=180.0),
                 boundary_conditions={'x': BC.open, 'y': BC.open}): 
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
      2. Attach the solver to a medium object RteSolver.set_medium(Medium)
      3. Run solution iterations with RteSolver.solve(maxiter)
    
    Parameters
    ----------
    scene_params: shdom.SceneParameters
        An object specifying scene parameters (solar direction, ground albedo...)
    numerical_params: shdom.NumericalParameters
        An object specifying numerical parameters (number of discrete ordinates, solution accuracy...)
    num_stokes: int
        The number of stokes for which to solve the RTE can be 1, 3, or 4. 
        num_stokes=1 means unpolarized.
        num_stokes=3 means linear polarization.
        num_stokes=4 means full polarization.
    name: str, optional 
       The name for the solver. Will be used when printing solution iteration messages. 
       If non specified a default <type> <wavelength> is given.
       
    Notes
    -----
    k-distribution not supported.
    """
    
    def __init__(self, scene_params=None, numerical_params=None, num_stokes=1, name=None):
        
        self._name = name 
        
        assert num_stokes in [1, 3, 4], 'num_stokes should be {1, 3, 4}'
        if num_stokes==1:
            self._type = 'Radiance'
        else:
            self._type = 'Polarization'
        self._nstokes = num_stokes
        
        
        # Start mpi (if available).
        self._masterproc = core.start_mpi()

        # Link to the properties array module.
        self._pa = ShdomPropertyArrays()
        
        if scene_params:
            self.set_scene(scene_params)
        if numerical_params:
            self.set_numerics(numerical_params)
            
    def get_param_dict(self):
        """
        Retrieve a dictionary with the solver parameters
        
        Returns
        -------
        param_dict: dict
            The parameters dictionary contains: type, scene_params, numerical_params
        """
        param_dict = {
            '_type': self.type,
            'solver_parameters': [
                {'scene_params': self._scene_parameters,
                 'numerical_params': self._numerical_parameters}
            ]
        }
        return param_dict

    def set_param_dict(self, param_dict):
        """
        Set the solver parameters from a parameter dictionary
        
        Parameters
        ----------
        param_dict: dict
            The parameters dictionary contains: type, scene_params, numerical_params
        """ 
        params = param_dict['solver_parameters'][0]
        self._type = param_dict['_type']
        self._nstokes = 1 if self.type == 'Radiance' else 3
        self.set_numerics(params['numerical_params'])
        self.set_scene(params['scene_params'])

    def save_params(self, path):
        """
        Save RteSolver parameters to file.
    
        Parameters
        ----------
        path: str,
            Full path to file. 
        """
        file = open(path,'wb')
        file.write(pickle.dumps(self.get_param_dict(), -1))
        file.close()

    def load_params(self, path):
        """
        Load RteSolver parameters from file.

        Parameters
        ----------
        path: str,
            Full path to file. 
        """        
        file = open(path,'rb')
        data = file.read()
        file.close()
        params = pickle.loads(data)
        self.set_param_dict(params)

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
        
        # Set a default name for the solver
        if self._name is None:
            self._name = '{} {:1.3f} micron'.format(self._type, self.wavelength)        
            
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
        if scene_params.boundary_conditions['x'] == BC.open:
            self._bcflag += 1
        if scene_params.boundary_conditions['y'] == BC.open:
            self._bcflag += 2
        
        # Independent pixel not supported yet
        self._ipflag = 0
        
        # k-distribution not supported yet
        self._kdist = False
        self._ng = 1
        self._delg = np.ones(1, order='F')
        self._pa.nzckd = 0
        self._baseout = False
        self._npart = 1

        # No iterations have taken place
        self._iters = 0
        
    def set_numerics(self, numerical_params):
        """
        Set the numerical parameters of the SHDOM forward solver.
    
        Parameters
        ----------
        numerical_params : NumericalParameters
            An object which encapsulate numerical parameters such as number of azimuthal and zenith bins. 
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

    def set_medium(self, medium):
        """
        Set the optical medium properties. 
        If a scatterer within the medium is a MultispectralScatterer or MicrophysicalScatterer object then the OpticalScatterer 
        at the solver's wavelength is extracted out of the medium (see get_optical_scatterer method within these classes).
        In addition this method initializes an internal grid and transfers the property array.
        
        Parameters
        ----------
        medium: shdom.Medium
            a Medium object conatining the optical properties.
            
        Notes
        -----
        The solution iteration criteria is set to 1.0 (restarted)
        Temperature (used for thermal radiation) is not supported (set to zero).
        """
        self.set_grid(medium.grid)
        
        # Temperature is used for thermal radiation. Not supported yet.
        self._pa.tempp = np.zeros(shape=(self._nbpts,), dtype=np.float32)     
        
        self._pa.extinctp = np.zeros(shape=[self._nbpts, medium.num_scatterers], dtype=np.float32)
        self._pa.albedop = np.zeros(shape=[self._nbpts, medium.num_scatterers], dtype=np.float32)
        self._pa.iphasep = np.zeros(shape=[self._nbpts, medium.num_scatterers], dtype=np.float32)
        
        for i, scatterer in enumerate(medium.scatterers.values()):
            
            if isinstance(scatterer, shdom.MicrophysicalScatterer) or isinstance(scatterer, shdom.MultispectralScatterer):
                scatterer = scatterer.get_optical_scatterer(self.wavelength)
            resampled_scatterer = scatterer.resample(medium.grid)
            
            self._pa.extinctp[:, i] = resampled_scatterer.extinction.data.ravel()
            self._pa.albedop[:, i] = resampled_scatterer.albedo.data.ravel()
            self._pa.iphasep[:, i] = resampled_scatterer.phase.iphasep.ravel() + self._pa.iphasep.max()

            scat_table = copy.deepcopy(resampled_scatterer.phase.legendre_table)
            if i == 0:
                legendre_table = scat_table
            else:
                legendre_table.append(scat_table)

        self._pa.numphase = legendre_table.numphase
                 
        # Determine the number of legendre coefficient for a given angular resolution
        if self._deltam:
            self._nleg = self._mm+1
        else:
            self._nleg = self._mm
        self._nleg = self._maxleg = max(legendre_table.maxleg, self._nleg)
        self._nscatangle = max(36, min(721, 2*self._nleg))

        # Legenp is without the zero order term which is 1.0 for normalized phase function
        self._pa.legenp = legendre_table.get_legenp(self._nleg).astype(np.float32)
        self._maxasym = legendre_table.maxasym
        self._maxpgl = medium.grid.num_points * legendre_table.maxleg         

        if legendre_table.numphase > 0:
            self._maxigl = legendre_table.numphase*(legendre_table.maxleg + 1)
        else:
            self._maxigl = self._maxig*(legendre_table.maxleg + 1)
    
        self._npart = medium.num_scatterers
        self._nstleg = legendre_table.nstleg
        self._nstphase = min(self._nstleg, 2)

        self._temp, self._planck, self._extinct, self._albedo, self._legen, self._iphase, \
            self._total_ext, self._extmin, self._scatmin, self._albmax = core.transfer_pa_to_grid(
                nstleg=self._nstleg,
                npart=self._npart,
                extinctp=self._pa.extinctp,
                albedop=self._pa.albedop,
                iphasep=self._pa.iphasep,                
                delx=self._pa.delx,
                dely=self._pa.dely,                
                xstart=self._pa.xstart,
                ystart=self._pa.ystart,
                zlevels=self._pa.zlevels,             
                tempp=self._pa.tempp,          
                legenp=self._pa.legenp,                           
                nzckd=self._pa.nzckd,
                zckd=self._pa.zckd,
                gasabs=self._pa.gasabs, 
                ml=self._ml,
                mm=self._mm,
                numphase=self._pa.numphase,
                deltam=self._deltam,
                units=self._units,
                waveno=self._waveno,
                wavelen=self._wavelen,
                gridpos=self._gridpos,
                nleg=self._nleg,
                maxig=self._maxig,
                npx=self._pa.npx,
                npy=self._pa.npy,
                npz=self._pa.npz,
                srctype=self._srctype,
                npts=self._npts)
        
        # Restart solution criteria
        self._solcrit = 1.0

    def set_grid(self, grid):
        """
        Set the base grid and related grid structures for SHDOM.
        
        Parameters
        ----------
        grid: shdom.Grid
            The grid
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
        self._nbpts = self._nx * self._ny * self._nz
        
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
                maxig=self._maxig, 
                maxic=self._maxic,
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                nx=self._nx,
                ny=self._ny,
                nz=self._nz,
                nx1=self._nx1,
                ny1=self._ny1,
                xgrid=self._xgrid,
                ygrid=self._ygrid,
                zgrid=self._zgrid
            )
        self._nbcells = self._ncells

    def init_memory(self):
        """A utility function to initialize internal memory parameters."""
        
        # Make ml and mm from nmu and nphi
        # ML is the maximum meridional mode, MM is the maximum azimuthal mode,
        # and NCS is the azimuthal mode flag (|NCS|=1 for cosine only, |NCS|=2 for 
        # sines and cosines).
        # nphi0max: The maximum number of azimuth angles actually used;
        # for NCS=1 (cosine modes only) NPHI0=INT((NPHI+2)/2),
        # otherwise NPHI0=NPHI.
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
            print('MAX_TOTAL_MB reduced to fit memory model: ', self._max_total_mb)
    
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
            print('adapt_grid_factor reduced to ', self._adapt_grid_factor)
    
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

    def init_solution(self):
        """
        Initilize the solution (I, J fields) from the direct transmition and a simple layered model. 
        Many arrays are initialized including the dirflux (the direct solar transmition).
        """
        self._oldnpts = 0
        self._solcrit = 1.0
        self._nang, self._nphi0, self._mu, self._phi, self._wtdo,  \
            self._sfcgridparms, self._ntoppts, self._nbotpts, self._bcptr, self._rshptr, self._shptr, self._oshptr, \
            self._source, self._delsource, self._radiance, self._fluxes, self._dirflux, self._ylmsun, self._uniformzlev, \
            self._pa.extdirp, self._bcrad, self._cx, self._cy, self._cz, self._cxinv, self._cyinv, self._czinv, \
            self._ipdirect, self._di, self._dj, self._dk, self._epss, self._epsz, self._xdomain, self._ydomain,  \
            self._delxd, self._delyd, self._deljdot, self._deljold, self._deljnew, self._jnorm, self._fftflag, \
            self._cmu1, self._cmu2, self._wtmu, self._cphi1, self._cphi2, self._wphisave, self._work, self._work1, \
            self._work2, self._uniform_sfc_brdf, self._sfc_brdf_do = core.init_solution(
                nstleg=self._nstleg,
                nstokes=self._nstokes,
                nx=self._nx,
                ny=self._ny,
                nz=self._nz,
                nx1=self._nx1,
                ny1=self._ny1,
                ml=self._ml,
                mm=self._mm, 
                maxiv=self._maxiv,
                maxic=self._maxic,
                maxig=self._maxig,
                maxido=self._maxido, 
                shacc=self._shacc,
                ncs=self._ncs,
                nlm=self._nlm,
                nmu=self._nmu,
                nphi=self._nphi,
                nleg=self._nleg,
                numphase=self._pa.numphase,  
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                deltam=self._deltam,
                srctype=self._srctype,
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
                units=self._units,
                waveno=self._waveno,
                wavelen=self._wavelen,
                accelflag=self._accelflag,
                xgrid=self._xgrid,
                ygrid=self._ygrid,
                zgrid=self._zgrid,
                temp=self._temp,
                planck=self._planck,
                extinct=self._extinct,
                albedo=self._albedo,
                legen=self._legen,
                iphase=self._iphase,
                maxnbc=self._maxnbc,
                maxbcrad=self._maxbcrad,
                npts=self._npts,
                gridpos=self._gridpos,
                ncells=self._ncells,
                gridptr=self._gridptr,
                treeptr=self._treeptr,                
                npx=self._pa.npx,
                npy=self._pa.npy,
                npz=self._pa.npz,   
                delx=self._pa.delx,
                dely=self._pa.dely,                
                xstart=self._pa.xstart,
                ystart=self._pa.ystart,
                zlevels=self._pa.zlevels,             
                tempp=self._pa.tempp,
                extinctp=self._pa.extinctp,
                nbpts=self._nbpts,
                albedop=self._pa.albedop,
                legenp=self._pa.legenp,
                iphasep=self._pa.iphasep,
                nzckd=self._pa.nzckd,
                zckd=self._pa.zckd,
                gasabs=self._pa.gasabs,                 
                npart=self._npart,
                total_ext=self._total_ext,
                nbcells=self._nbcells,
                maxsfcpars=self._maxsfcpars,
                nphi0max=self._nphi0max
            )
                
    def solution_iterations(self, maxiter, verbose=True):
        """
        Run SHDOM solution iterations process.

        Parameters
        ----------
        maxiter: integer
            Maximum number of iterations for the iterative solution.
        verbose: verbosity
            True will output solution iteration information into stdout.
            
        Returns
        -------
        output_arguments: tuple,
            A tuple containing all the output arguments, see: update_solution_arguments().
        """  
        output_arguments = core.solution_iterations(
            uniform_sfc_brdf=self._uniform_sfc_brdf,
            sfc_brdf_do=self._sfc_brdf_do,
            work=self._work,
            work1=self._work1,
            work2=self._work2,
            bcrad=self._bcrad,
            fluxes=self._fluxes,
            nang=self._nang,
            nphi0=self._nphi0,
            maxnbc=self._maxnbc,
            ntoppts=self._ntoppts,
            nbotpts=self._nbotpts,
            uniformzlev=self._uniformzlev,
            extmin=self._extmin, 
            scatmin=self._scatmin, 
            cx=self._cx, 
            cy=self._cy, 
            cz=self._cz, 
            cxinv=self._cxinv,
            cyinv=self._cyinv, 
            czinv=self._czinv, 
            ipdirect=self._ipdirect, 
            di=self._di, 
            dj=self._dj, 
            dk=self._dk, 
            nphi0max=self._nphi0max, 
            epss=self._epss, 
            epsz=self._epsz,
            xdomain=self._xdomain, 
            ydomain=self._ydomain, 
            delxd=self._delxd, 
            delyd=self._delyd, 
            albmax=self._albmax, 
            deljdot=self._deljdot, 
            deljold=self._deljold, 
            deljnew=self._deljnew, 
            jnorm=self._jnorm,                
            fftflag=self._fftflag,
            cmu1=self._cmu1,
            cmu2=self._cmu2,
            wtmu=self._wtmu,
            cphi1=self._cphi1,
            cphi2=self._cphi2,
            wphisave=self._wphisave,
            nbpts=self._nbpts,
            npart=self._npart,
            extinct=self._extinct,
            albedo=self._albedo, 
            legen=self._legen, 
            total_ext=self._total_ext,
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
            source=self._source,
            delsource=self._delsource,
            radiance=self._radiance,
            dirflux=self._dirflux,
            nleg=self._nleg,
            maxiv=self._maxiv,
            maxic=self._maxic,
            maxig=self._maxig,
            maxido=self._maxido,
            verbose=verbose,
            oldnpts=self._oldnpts,
            ylmsun=self._ylmsun,
            runname=self._name
        )
        return output_arguments
    
    def make_direct(self):
        """
        Compute the direct transmition from the solar beam.
        """
        self._dirflux, self._pa.extdirp, self._cx, self._cy, self._cz, self._cxinv, self._cyinv, self._czinv, \
            self._di, self._dj, self._dk, self._ipdirect, self._delxd, self._delyd, self._xdomain, self._ydomain, \
            self._epss, self._epsz, self._uniformzlev = \
            core.make_direct(
                nstleg=self._nstleg,
                dirflux=self._dirflux, 
                extdirp=self._pa.extdirp,
                npts=self._npts,
                bcflag=self._bcflag,
                ipflag=self._ipflag,
                deltam=self._deltam,
                ml=self._ml,
                nleg=self._nleg,
                solarflux=self._solarflux,
                solarmu=self._solarmu,
                solaraz=self._solaraz,
                gridpos=self._gridpos,
                npx=self._pa.npx,
                npy=self._pa.npy,
                npz=self._pa.npz,
                numphase=self._pa.numphase,
                delx=self._pa.delx,
                dely=self._pa.dely,
                xstart=self._pa.xstart,
                ystart=self._pa.ystart,
                zlevels=self._pa.zlevels,
                tempp=self._pa.tempp,
                extinctp=self._pa.extinctp,
                albedop=self._pa.albedop,
                legenp=self._pa.legenp,
                iphasep=self._pa.iphasep,
                nzckd=self._pa.nzckd,
                zckd=self._pa.zckd,
                gasabs=self._pa.gasabs
            )
        
    def update_solution_arguments(self, solution_arguments):
        """
        Update the solution arguments from the output of the solution_iteration method.
        """
        self._sfcgridparms, self._solcrit, \
            self._iters, self._temp, self._planck, self._extinct, self._albedo, self._legen, self._iphase, \
            self._ntoppts, self._nbotpts, self._bcptr, self._bcrad, self._npts, self._gridpos, self._ncells, self._gridptr, \
            self._neighptr, self._treeptr, self._cellflags, self._rshptr, self._shptr, self._oshptr, self._source, \
            self._delsource, self._radiance, self._fluxes, self._dirflux, self._uniformzlev, self._pa.extdirp, \
            self._oldnpts, self._total_ext, self._deljdot, self._deljold, self._deljnew, self._jnorm, \
            self._work, self._work1, self._work2 = solution_arguments 
          
    def solve(self, maxiter, init_solution=False, verbose=True):
        """
        Main solver routine. This routine is comprised of two parts:
          1. Initialization (init_solution method), optional
          2. Solution iterations (solution_iterations method followed by update_solution_arguments method)

        Parameters
        ----------
        maxiter: integer
            Maximum number of iterations for the iterative solution.
        init_solution: boolean, default=False
            If False the solution is initialized according to the existing radiance and source function saved within the RteSolver object (previously computed)
            If True or no prior solution (I,J fields) exists then an initialization is preformed (part 1.).
        verbose: boolean
            True will output solution iteration information into stdout.
        """
        if self.num_iterations == 0 or init_solution:
            self.init_solution()
        solution_arguments = self.solution_iterations(maxiter, verbose)
        self.update_solution_arguments(solution_arguments)
    
    @property
    def name(self):
        return self._name    

    @property
    def wavelength(self):
        return float_round(self._wavelen)

    @property
    def type(self):
        return self._type
    
    @property
    def num_iterations(self):
        return self._iters
    
    @property
    def info(self):
        return self._scene_parameters.info + os.linesep + self._numerical_parameters.info
    
    
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
        self._name = []
        self._type = None
        if solver_list is not None:
            for solver in solver_list:
                self.add_solver(solver)
                
    def set_medium(self, medium):
        """
        Set the optical medium properties for all rte_solvers within the list.
        medium.wavelength must match solver.wavelengths
        
        Parameters
        ----------
        medium: shdom.Medium
            a Medium object conatining the optical properties.
            
        Notes
        -----
        medium wavelengths must match solver wavelengths
        """
        solver_wavelengths = [solver.wavelength for solver in self.solver_list]
        assert np.allclose(medium.wavelength, solver_wavelengths), 'medium wavelength {} differs from solver wavelengh {}'.format(medium.wavelength, solver_wavelengths)
        for solver in self.solver_list:
            solver.set_medium(medium)

    def get_param_dict(self):
        """
        Retrieve a dictionary with the solver array parameters
        
        Returns
        -------
        param_dict: dict
            The parameters dictionary contains: type and a list of scene_params, numerical_params matching all the solvers
        """
        param_dict = {}
        for key, param in self.__dict__.items():
            if key != '_solver_list':
                param_dict[key] = param
            else:
                param_dict['solver_parameters'] = []
                for solver in param:
                    solver_param_dict = {'scene_params': solver._scene_parameters,
                                         'numerical_params': solver._numerical_parameters}                    
                    param_dict['solver_parameters'].append(solver_param_dict)        
        return param_dict
    
    def set_param_dict(self, param_dict):
        """
        Set the solver array parameters from a parameter dictionary
        
        Parameters
        ----------
        param_dict: dict
            The parameters dictionary contains: type and a list of scene_params, numerical_params matching all the solvers
        """ 
        num_stokes = 1 if param_dict['_type'] == 'Radiance' else 3
        for solver_params in param_dict['solver_parameters']:
            rte_solver = shdom.RteSolver(scene_params=solver_params['scene_params'], 
                                         numerical_params=solver_params['numerical_params'],
                                         num_stokes=num_stokes)
            self.add_solver(rte_solver)
            
    def save_params(self, path):
        """
        Save RteSolverArray parameters to file.
    
        Parameters
        ----------
        path: str,
            Full path to file. 
        """
        file = open(path,'wb')
        file.write(pickle.dumps(self.get_param_dict(), -1))
        file.close()
        
    def load_params(self, path):
        """
        Load RteSolverArray parameters from file.

        Parameters
        ----------
        path: str,
            Full path to file. 
        """        
        file = open(path,'rb')
        data = file.read()
        file.close()
        params = pickle.loads(data)
        self.set_param_dict(params)

    def __getitem__(self, val):
        return self.solver_list[val]
    
    def add_solver(self, rte_solver):
        """
        Add a rte_solver to the RteSolverArray
        
        Parameters
        ----------
        rte_solver: RteSolver object
            A RteSolver object to add to the RteSolverArray
        """
        
        if self.type is None:
            self._type = rte_solver.type
        else:
            assert self.type == rte_solver.type, \
                   '[add_solver] Assert: RteSolverArray is of type {} and new solver is of type {}'.format(self.type, rte_solver.type)
            
        self._solver_list.append(rte_solver)
        self._name.append(rte_solver.name)
        self._num_solvers += 1
      
    def init_solution(self):
        """
        Initilize the solution (I, J fields) from the direct transmition and a simple layered model. 
        Many arrays are initialized including the dirflux (the direct solar transmition).
        """
        for solver in self.solver_list:
            solver.init_solution() 
            
    def make_direct(self):
        """
        Compute the direct transmition from the solar beam.
        """
        for solver in self.solver_list:
            solver.make_direct()    
            
    def solve(self, maxiter, init_solution=False, verbose=True):
        """
        Parallel solving of all solvers.
        
        Main solver routine. This routine is comprised of two parts:
          1. Initialization (init_solution method), optional
          2. Parallel solution iterations (solution_iterations method followed by update_solution_arguments method)

        Parameters
        ----------
        maxiter: integer
            Maximum number of iterations for the iterative solution.
        init_solution: boolean, default=False
            If False the solution is initialized according to the existing radiance and source function saved within the RteSolver object (previously computed)
            If True or no prior solution (I,J fields) exists then an initialization is preformed (part 1.).
        verbose: boolean
            True will output solution iteration information into stdout.
        """
        for solver in self.solver_list:
            if init_solution or solver.num_iterations == 0:
                solver.init_solution()        

        output_arguments = \
            Parallel(n_jobs=self.num_solvers, backend="threading")(
                delayed(rte_solver.solution_iterations, check_pickle=False)(
                    maxiter, verbose) for rte_solver in self.solver_list)
     
        for solver, arguments in zip(self.solver_list, output_arguments):
            solver.update_solution_arguments(arguments)

    @property
    def name(self):
        return self._name    

    @property
    def solver_list(self):
        return self._solver_list
    
    @property
    def num_solvers(self):
        return self._num_solvers

    @property
    def type(self):
        return self._type