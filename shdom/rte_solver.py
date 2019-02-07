"""
Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer.
Here all the necessary functions used to solve the RT using SHDOM are wrapped.
"""

import core
import numpy as np
from enum import Enum
import warnings


class BoundaryCondition(Enum):
    open = 1       # open boundary conditions mean that exiting radiance is lost.
    periodic = 2   # periodic boundary conditions mean that exiting radiance returns from the opposite side.
    

class PhaseFunctionType(Enum):
    # TODO: Add documentation
    tabulated = 1 
    
    
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
        Solar beam zenith angle; Specified in degrees but immediately converted to the cosine of the angle (mu).
        This angle represents the direction of travel of the solar beam, so is forced to be negative although it can be specified positive.
    
    Notes
    -----
    Paramters: skyrad, units, are used for Thermal source and are not supported. 
    """
    def __init__(self, flux, azimuth, zenith):
        self.flux = flux
        self.azimuth = azimuth
        self.zenith = zenith
        self.units = 'R'
        self.skyrad = 0.0
        
    
class Surface(object):
    """ 
    An abstract sufrace class to be inhirted by different surface types.
    """
    def __init__(self):
        self.type = 'AbstractSurfaceClass'
       
        
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
        self.type = 'Lambertian'
        self.albedo = albedo
        self.ground_temperature = 298.15
        
class NumericalParameters(object):
    """
    TODO: add description

    Parameters
    ----------
    """    
    def __init__(self):
        self.maxnewphase = 200
        self.phase_function_type = PhaseFunctionType.tabulated
        
        # TODO: check these values
        self.asymetric_tolerance = 0.1
        self.frational_phase_tolerance = 0.1
        self.num_mu_bins = 8
        self.num_phi_bins = 16
        self.split_accuracy = 0.1
        self.delta_m = True
        self.spherical_harmonics_accuracy = 0.003
        self.solution_accuracy = 0.0001
        self.acceleration_flag = True
        self.max_total_mb = 100000.0
        self.adapt_grid_factor = 5
        self.num_sh_term_factor = 5
        self.cell_to_point_ratio = 1.5


class SceneParameters(object):
    """
    TODO: add description

    Parameters
    ----------
    """    
    def __init__(self): 
        self.wavelength = 0.670
        self.surface = LambertianSurface(albedo=0.05)
        self.source = SolarSource(flux=np.pi, azimuth=0.0, zenith=180.0)
        self.boundary_conditions = {'x': BoundaryCondition.open, 
                                    'y': BoundaryCondition.open}

        
class RteSolver(object):
    """
    TODO: add documentation

    Parameters
    ----------
    scene_params: SceneParameters
    numerical_params: NumericalParameters
    
    Notes
    -----
    k-distribution not supported.
    """
    
    def __init__(self, scene_params, numerical_params):

        # Start mpi (if available).
        self._masterproc = core.start_mpi()

        # Link to the properties array module.
        self._pa = core.shdom_property_arrays

        # Assign scene parameters to shdom internal structure
        self.set_scene_parameters(scene_params)
        
        # Assign numerical parameters to shdom internal structure
        self.set_numerical_parameters(numerical_params)

    
    
    def __del__(self):
        
        try:
            core.end_shdom_mpi(
                gridpos=self._gridpos,
                npx=self._pa.npx,
                npy=self._pa.npy,
                xstart=self._pa.xstart,
                ystart=self._pa.ystart,
                delx=self._pa.delx,
                dely=self._pa.dely,
                npxt=self._npxt,
                npyt=self._npyt,
                propfile=self._propfile
            )
        except Exception as e:
            warnings.warn(repr(e))

    def set_scene_parameters(self, scene_params):
        """
        Set the scene related parameters: 
          wavelength, source, surface, boundary conditions, k-distribution 
    
        Parameters
        ----------
        scene_params : SceneParameters
        
        Returns
        -------
        None
    
        Notes
        -----
        k-distribution not supported.
    
        """  
        
        # Wavelength
        self._wavelen = scene_params.wavelength
        
        # Source parameters
        self._solarflux = scene_params.source.flux
        self._solarmu = np.cos(np.deg2rad(scene_params.source.zenith))
        self._solaraz = np.deg2rad(scene_params.source.azimuth)
        self._skyrad = scene_params.source.skyrad
        self._units = scene_params.source.units
        
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
        
        # TODO: k-distribution not supported yet
        self._kdist = False
        self._ng = 1
        self._delg = np.ones(1, order='F')
        self._pa.nzckd = 0
        self._baseout = False


    def set_numerical_parameters(self, numerical_params):
        """
        Set the numerical parameters of the SHDOM forward solver.
    
        Parameters
        ----------
        numerical_params : NumericalParameters

        Returns
        -------
        None
    
        Notes
        -----
        """     
        self._maxnewphase         = numerical_params.maxnewphase
        self._proptype            = numerical_params.phase_function_type
        self._asymtol             = numerical_params.asymetric_tolerance   
        self._fracphasetol        = numerical_params.frational_phase_tolerance
        self._nmu                 = max(2, 2 * int((numerical_params.num_mu_bins + 1) / 2))
        self._nphi                = max(1, numerical_params.num_phi_bins)
        self._deltam              = numerical_params.delta_m
        self._max_total_mb        = numerical_params.max_total_mb
        self._adapt_grid_factor   = numerical_params.adapt_grid_factor
        self._accelflag           = numerical_params.acceleration_flag
        self._num_sh_term_factor  = numerical_params.num_sh_term_factor
        self._cell_to_point_ratio = numerical_params.cell_to_point_ratio
        self._splitacc            = numerical_params.split_accuracy
        self._shacc               = numerical_params.spherical_harmonics_accuracy
        self._solacc              = numerical_params.solution_accuracy
           
        
        # Make ml and mm from nmu and nphi
        # ML is the maximum meridional mode, MM is the maximum azimuthal mode,
        # and NCS is the azimuthal mode flag (|NCS|=1 for cosine only, |NCS|=2 for 
        # sines and cosines).
        # nphi0max: The maximum number of azimuth angles actually used;
        # for NCS=1 (cosine modes only) NPHI0=INT((NPHI+2)/2),
        # otherwise NPHI0=NPHI.
        self._ml = self._nmu-1
        self._mm = max(0, int(self._nphi/2)-1)
        self._nlm = (2*self._mm+1)*(self._ml+1) - self._mm*(self._mm+1)
        self._ncs = 2
        self._nphi0max = self._nphi
        self._nleg = self._ml
    
        if self._deltam:
            self._nleg += 1
            
        self._memword = self._nmu*(2+2*self._nphi+2*self._nlm+2*33*32)