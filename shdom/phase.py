"""
This is a python wrapper to handle phase function related computations. 

It includes a wrapper for mie and rayleigh computations. 
The python code was created by Aviad Levis Technion Inst. of Technology February 2019.

Source Fortran files were created by Frank Evans University of Colorado May 2003.
"""

import core
import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator
from shdom import Grid, GridData, BoundingBox


class Phase(object):
    """
    An abstract Phase object that is inherited by two types of Phase functions:
      1. GridPhase: The phase function is specified at every point on the grid.
      2. TabulatedPhase: A table of phase functions is saved with a pointer at every point on the grid.
    
    The TabulatedPhase is more efficient in memory and runtime.
    """
    def __init__(self):
        self._type = 'AbstractPhaseObject'
        self._maxleg = None
        self._iphasep = None
        self._numphase = None
        self._legenp = None
        self._maxasym = None
        self._grid = None
        
    @property
    def type(self):
        return self._type
    
    @property 
    def maxleg(self):
        return self._maxleg
    
    @property 
    def numphase(self):
        return self._numphase  
    
    @property
    def iphasep(self):
        return self._iphasep
    
    @property
    def legenp(self):
        return self._legenp
    
    @property
    def maxasym(self):
        return self._maxasym    
    
    @property
    def grid(self):
        return self._grid
    
    
class GridPhase(Phase):
    """
    A GridPhase object spefies the phase function at every point in the grid. 
    This is the equivalent of the GridData object for the phase function.
    
    Parameters
    ----------
    grid: Grid object
        A Grid object of type '1D' or '3D'.
    data: np.array
        data contains the vector field (legendre coefficients).
    """
    def __init__(self, grid, data):
        self._type = 'Grid'
        self._data = data
        self._grid = grid 
        self._linear_interpolator1d = interp1d(grid.z, self.data, assume_sorted=True, copy=False, bounds_error=False, fill_value=0.0) 
        if grid.type == '3D':
            self._linear_interpolator3d = RegularGridInterpolator((grid.x, grid.y, grid.z), self.data.transpose([1, 2, 3, 0]), bounds_error=False, fill_value=0.0)
    
        self._numphase = grid.num_points
        self._iphasep = np.arange(1, grid.num_points+1, dtype=np.int32)
        
        # Legenp is without the zero order term which is 1.0 for normalized phase function
        self._legenp = data[1:].ravel(order='F')       
        self._maxleg = data.shape[0] - 1
        
        # Asymetry parameter is proportional to the legendre series first coefficient 
        self._maxasym = data[1,...].max() / 3.0


    def __add__(self, other):
        """Adding two GridPhase objects by resampling to a common grid and padding with zeros to the larger legendre series."""
        assert other.__class__ is GridPhase, 'Only GridPhase can be added to GridPhase object'
        grid = self.grid + other.grid
        self_data = self.resample(grid)
        other_data = other.resample(grid)
        depth_diff = self.maxleg - other.maxleg
        if depth_diff > 0:
            data = self_data + np.pad(other_data, ((0, depth_diff), (0, 0), (0, 0), (0, 0)), 'constant')
        elif depth_diff < 0:
            data = other_data + np.pad(self_data, ((0, -depth_diff), (0, 0), (0, 0), (0, 0)), 'constant')
        return GridPhase(grid, data)

    
    def __mul__(self, other):
        """Multiplying a GridPhase objects by a GridData object."""
        assert other.__class__ is GridData, 'Only scalar field (GridData) can multiply a GridPhase object'
        grid = self.grid + other.grid
        other_data = other.resample(grid)
        data = self.resample(grid) * other_data[np.newaxis,...]  
        return GridPhase(grid, data)
    
    
    def __div__(self, other):  
        """Dividing a GridPhase objects by a GridData object."""
        assert other.__class__ is GridData, 'Only scalar field (GridData) can divide a GridPhase object'
        grid = self.grid + other.grid
        other_data = other.resample(grid)        
        data = self.resample(grid) / other_data[np.newaxis,...] 
        return GridPhase(grid, data)    
    
    
    def resample(self, grid):
        """Resample data to a new Grid."""
        if self.grid.type == '1D':
            if np.array_equiv(self.grid.z, grid.z):
                return self.data
            data = self._linear_interpolator1d(grid.z)
        else:
            if self.grid == grid:
                return self.data
            data = self._linear_interpolator3d(np.stack(np.meshgrid(grid.x, grid.y, grid.z, indexing='ij'), axis=-1)).transpose([3, 0, 1, 2])
        return data
    
    
    @property
    def data(self):
        return self._data

    
class TabulatedPhase(Phase):
    """
    The TabulatedPhase internally keeps a phase function table and a pointer array for each grid point.
    This could potentially be more efficient than GridPhase in memory and computation if the number of table entery is much smaller than number of grid points.
    
    Parameters
    ----------
    legendre_table: np.array(shape=(nleg, numphase), dtype=np.float32)
       The Legendre table.
    index: GridData object
       A GridData object with dtype=int. This is a pointer to the enteries in the legendre_table.
    """
    def __init__(self, legendre_table, index):
        self._type = 'Tabulated'
        self._legendre_table = np.array(legendre_table, dtype=np.float32)
        self._index = index
        self._grid = index.grid

        self._numphase = legendre_table.shape[1]
        self._iphasep = index.data
        
        # Legenp is without the zero order term which is 1.0 for normalized phase function
        self._legenp = legendre_table[1:].ravel(order='F')        
        self._maxleg = legendre_table.shape[0] - 1
        
        # Asymetry parameter is proportional to the legendre series first coefficient 
        self._maxasym = legendre_table[1].max() / 3.0

    
    def add_ambient(self, other):
        """
        Adding an ambient medium means that the AmientMedium will only `fill the holes`. 
        This is an approximation which speeds up computations.
        """
        # Join the two tables
        self_legendre_table = self.legendre_table
        other_legendre_table = other.legendre_table
        if self.maxleg > other.maxleg:
            other_legendre_table = np.pad(other_legendre_table, ((0, self.maxleg - other.maxleg), (0, 0)), 'constant')
        elif other.maxleg > self.maxleg:
            self_legendre_table = np.pad(self_legendre_table, ((0, other.maxleg - self.maxleg), (0, 0)), 'constant')
        legendre_table = np.concatenate((other_legendre_table, self_legendre_table), axis=1)        
        
        # Join the indices
        grid = self.grid + other.grid 
        self_index = self.index.resample(grid, method='nearest')
        other_index = other.index.resample(grid, method='nearest')
        self_index[self_index>0] += other_index.max()
        self_index[self_index==0] = (self_index + other_index)[self_index==0]
        index = GridData(grid, self_index.astype(np.int32))
        return TabulatedPhase(legendre_table, index)  

    @property
    def legendre_table(self):
        return self._legendre_table
    
    @property
    def index(self):
        return self._index
    
    
class Mie(object):
    """
    Mie scattering for a particle size distribution. 
    Scattering coefficients are averaged over a range of particle radii and wavelengths. 
    A Mie scattering table can be computed by using:
      1. Mie.set_parameters(...)
      2. Mie.compute_talbe(...)
    Alternatively a pre-computed scattering table can be loaded using:
    Mie.read_table(table_path)
    See: /notebooks/Make Mie Table.ipynb for usage examples. 
    Once a scattering table is computed/loaded, Mie.interpolate_scattering_field(lwc, reff) can be used
    to compute the extinction, singe scattering albedo and phase function on a grid. 
    
    Description taken from make_mie_table.f90:
      Does Mie computations to create a scattering table as a function of
      effective radius for gamma or lognormal size distributions of spherical
      particles.  The particles may be water or ice (in which case the 
      program provides the index of refraction depending on wavelength) or
      "aerosols" (in which case the index of refraction is user specified).
      For water or ice particles the scattering properties may be averaged 
      over the desired spectral range with Planck function weighting.  
      The phase functions in the output scattering table are represented 
      with Legendre series.   The effective radii in the table may be evenly
      or logarithmically spaced.
   """
    def __init__(self):
        self._partype = None
        self._rindex = None
        self._pardens = None
        self._distflag = None
        self._wavelen1 = None
        self._wavelen2 = None
        self._avgflag = None  
        self._deltawave = None
        self._alpha = None
        self._wavelencen = None
        self._reff = None
        self._veff = None
        self._extinct = None
        self._ssalb = None
        self._nleg = None
        self._legcoef = None
    
    def set_parameters(self,
                       wavelength_band,
                       particle_type, 
                       distribution,
                       alpha, 
                       wavelength_averaging=False,
                       wavelength_resolution=0.001):
        """
        Set the Mie parameters to compute a new scattering table.
        
        Parameters
        ----------
        wavelength_band: (float, float)
            (minimum, maximum) wavelength in microns. 
            This defines the spectral band over which to integrate, if both are equal monochrome quantities are computed. 
        particle_type: string
            Options are 'Water' or 'Aerosol'.
        distribution: string
            Particle size-distribution. Options are 'Gamma' or 'Log-normal'. 
            Gamma:
              n(r) = a r^alpha exp(-b*r).
              r - droplet radius.
              a, b, alpha - gamma distribution parameters. 
            Log-normal:
              n(r) = a/r exp( -[ln(r/r0)]^2 / (2*alpha^2) ).
              r0 - logarithmic mode.
              alpha - standard deviation of the log. 
        alpha: float
            Shape parameter for the size distribution. 
            Gamma:
              N = a Gamma(alpha+1)/ b^(alpha+1) - number concentration.
              r_eff = (alpha+3)/b - effective radius.
              v_eff = 1/(alpha+3) - effective variance.
            Log-normal:
              N = sqrt(2*pi)*alpha*a - number concentration. 
              r_eff = r0*exp(2.5*alpha^2) - effective radius.
              v_eff = exp(alpha^2)-1  - effective variance.
        wavelength_averaging: bool
            True - average scattering properties over the wavelength_band.
            False - scattering properties of the central wavelength. 
        wavelength_resolution: float
            The distance between two wavelength samples in the band. Used only if wavelength_averaging is True.
        
        Returns
        -------
        None
            
        Notes
        -----
        Aerosol particle type not supported yet.          
        """
        # Particle type 
        if particle_type == 'Water':
            self._partype = 'W'
            self._rindex = 1.0
            self._pardens = 1.0
        else:
            raise NotImplementedError('Particle type {} not supported'.format(particle_type))
        
        # Size distribution pdf 
        if distribution == 'gamma':
            self._distflag = 'G'
        elif distribution == 'lognormal':
            self._distflag = 'L'
        else:
            raise NotImplementedError('Distribution type {} not supported'.format(distibution))
            
        # Averaging over spectral band or monochrome
        self._wavelen1, self._wavelen2 = wavelength_band
        assert self._wavelen1 <= self._wavelen2, 'Minimum wavelength is smaller than maximum'
        avgflag = 'C'
        if self._wavelen1 == self._wavelen2:
            deltawave = -1
        elif wavelength_averaging:
            avgflag = 'A'
            deltawave = wavelength_resolution
            
        self._avgflag = avgflag  
        self._deltawave = deltawave
        self._alpha = alpha
        
        self._wavelencen = core.get_center_wavelen(
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2)
        
    
    def compute_table(self,
                      num_effective_radii,
                      start_effective_radius,
                      end_effective_radius,
                      max_integration_radius): 
        """
        Compute a scattering table where for each effective radius:
          1. Extinction-cross section per 1 unit of mass content [g/m^3](liquid water content for water clouds)  
          2. Single scattering albedo, unitless in the range [0, 1].
          3. Legendre expansion coefficients of the normalized scattering phase function (first coefficient is always 1.0)
          4. Number of Legendre coefficients for each scattering phase function. 
    
        Parameters
        ----------
        num_effective_radii: int
            Number of effective radii for which to compute the table.
        start_effective_radius: int
            The starting (lowest) effective radius in the table.
        end_effective_radius: int
            The ending (highest) effective radius in the table.
        max_integration_radius: int
            The maximum radius for which to integrate over the size-distribution.
            max_integration_radius > end_effective_radius.
        
    
        Notes
        -----
        Running this function may take some time.
        """
        
        assert None not in (self._partype, self._rindex, self._pardens,
                            self._distflag, self._wavelen1, self._wavelen2,
                            self._avgflag, self._deltawave, self._alpha, self._wavelencen), \
               'Mie parameters were not set. Set them using set_parameters() setter.'
        
        self._nretab = num_effective_radii
        self._sretab = start_effective_radius
        self._eretab = end_effective_radius
        self._maxradius = max_integration_radius
        
        
        # Calculate the maximum size parameter and the max number of Legendre terms
        if self._avgflag == 'A':
            xmax = 2 * np.pi * max_integration_radius / self._wavelen1
        else:
            xmax = 2 * np.pi * max_integration_radius / self._wavelencen
        self._maxleg = int(np.round(2.0 * (xmax + 4.0*xmax**0.3334 + 2.0)))
        
        print('Computing mie table...')
        self._reff, self._extinct, self._ssalb, self._nleg, self._legcoef = core.get_mie_table(
            nretab=self._nretab, 
            maxleg=self._maxleg,
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2, 
            wavelencen=self._wavelencen,
            deltawave=self._deltawave, 
            pardens=self._pardens, 
            sretab=self._sretab, 
            eretab=self._eretab, 
            alpha=self._alpha, 
            maxradius=self._maxradius, 
            rindex=self._rindex, 
            partype=self._partype, 
            avgflag=self._avgflag, 
            distflag=self._distflag)
        
        # A simple hack: duplicate the table for two grid points (effectinve radii). 
        # This is because two or more grid points are requiered for interpolation. 
        if self._nretab == 1:
            self._reff = np.tile(self._reff, 2)
            self._extinct = np.tile(self._extinct, 2)
            self._ssalb = np.tile(self._ssalb, 2)
            self._nleg = np.tile(self._nleg, 2)
            self._legcoef = np.tile(self._legcoef, 2)
            
        self.init_intepolators()
        
        print('Done.')
        
    def write_table(self, file_path): 
        """
        Write a pre-computed table to <file_path>. 
    
        Parameters
        ----------
        file_path: str
            Path to file.
      
        Returns
        -------
        None

        Notes
        -----
        This function must be ran after pre-computing a scattering table with compute_table().
        """        
        print('Writing mie table to file: {}'.format(file_path))
        core.write_mie_table(
            mietabfile=file_path,
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2,
            deltawave=self._deltawave,
            partype=self._partype,
            pardens=self._pardens, 
            rindex=self._rindex,
            distflag=self._distflag,
            alpha=self._alpha, 
            nretab=self._nretab,
            sretab=self._sretab, 
            eretab=self._eretab,             
            reff=self._reff,
            extinct=self._extinct,
            ssalb=self._ssalb,
            nleg=self._nleg,
            legcoef=self._legcoef,
            maxleg=self._maxleg)
        print('Done.')
     
    def read_table(self, file_path): 
        """
        Read a pre-computed table from <file_path>. 
    
        Parameters
        ----------
        file_path: str
            Path to file.
      
        Returns
        -------
        None

        """   
        
        def read_table_header(file_path):
            wavelen1, wavelen2, deltawave = np.genfromtxt(file_path, max_rows=1, skip_header=1, usecols=(0, 1, 2), dtype=float)
            pardens = np.genfromtxt(file_path, max_rows=1, skip_header=2, usecols=(0), dtype=float)
            partype = np.asscalar(np.genfromtxt(file_path, max_rows=1, skip_header=2, usecols=(1), dtype=str))
            rindex  = np.complex(np.genfromtxt(file_path, max_rows=1, skip_header=3, usecols=(0), dtype=float), 
                                 np.genfromtxt(file_path, max_rows=1, skip_header=3, usecols=(1), dtype=float))
            alpha = np.genfromtxt(file_path, max_rows=1, skip_header=4, usecols=(0), dtype=float)
            
            distribution = np.asscalar(np.genfromtxt(file_path, max_rows=1, skip_header=4, usecols=(1), dtype=str))
            if distribution == 'gamma':
                distflag = 'G'
            elif distribution == 'lognormal':
                distflag = 'L'
            else:
                raise NotImplementedError('Distribution type {} not supported'.format(distibution))
            
            nretab = np.genfromtxt(file_path, max_rows=1, skip_header=5, usecols=(0), dtype=int)
            sretab, eretab = np.genfromtxt(file_path, max_rows=1, skip_header=5, usecols=(1, 2), dtype=float)
            maxleg = np.genfromtxt(file_path, max_rows=1, skip_header=6, usecols=(0), dtype=int)
            
            return wavelen1, wavelen2, deltawave, pardens, partype, rindex, alpha, distflag, nretab, sretab, eretab, maxleg        


        
        print('Reading mie table from file: {}'.format(file_path))
        self._wavelen1, self._wavelen2, self._deltawave, self._pardens, \
            self._partype, self._rindex, self._alpha, self._distflag, \
            self._nretab, self._sretab, self._eretab, self._maxleg = read_table_header(file_path)
        
        self._reff, self._extinct, self._ssalb, self._nleg, self._legcoef = \
            core.read_mie_table(mietabfile=file_path, 
                                nretab=self._nretab, 
                                maxleg=self._maxleg)
        
        self.init_intepolators()
        print('Done.')

    
    def init_intepolators(self):
        """
        Initialize interpolators for the extinction, single scattering albedo and phase function.
    
        Notes
        -----
        This function is ran after pre-computing/loading a scattering tables.
        """        
        assert True not in (self._reff is None, self._extinct is None, self._ssalb is None, self._nleg is None, self._legcoef is None), \
                       'Mie scattering table was not computed or read from file. Using compute_table() or read_table().'   
        
        if self._nretab == 1:
            kind = 'nearest'
        elif self._nretab > 1:
            kind = 'linear'
        else:
            raise AttributeError
        
        self._ext_interpolator = interp1d(self.reff, self.extinct, kind=kind, assume_sorted=True, bounds_error=False, fill_value=0.0)
        self._ssalb_interpolator = interp1d(self.reff, self.ssalb, kind=kind, assume_sorted=True, bounds_error=False, fill_value=1.0) 
        self._legcoef_interpolator = interp1d(self.reff, self.legcoeff, kind=kind, assume_sorted=True, bounds_error=False, fill_value=0.0)          
        self._legen_index_interpolator = interp1d(self.reff, range(1, len(self.reff)+1), kind='nearest', assume_sorted=True, copy=False, bounds_error=False, fill_value=0)
        self._nleg_interpolator = interp1d(self.reff, self.nleg, kind='nearest', assume_sorted=True, copy=False, bounds_error=False, fill_value=0)
     
     
    def interpolate_extinction(self, lwc, reff):
        """
        Interpolate the extinciton coefficient function over a grid.
        
        Parameters
        ----------
        extinction: GridData object
            A GridData3D object containting the extinction (1/km) on a 3D grid
        
        Returns
        -------
        extinction: GridData object
            A GridData3D object containting the extinction (1/km) on a 3D grid
        """
        data = self._ext_interpolator(reff.data)
        extinction = lwc * GridData(reff.grid, data)     
        return extinction     

        
    def interpolate_albedo(self, reff):
        """
        Interpolate the single scattering albedo over a grid.
        
        Parameters
        ----------
        reff: GridData 
            A GridData object containting the effective radii (micron) on a 3D grid.
        
        Returns
        -------
        albedo: GridData object
            A GridData object containting the single scattering albedo unitless in range [0, 1] on a 3D grid
        """
        data = self._ssalb_interpolator(reff.data)
        albedo = GridData(reff.grid, data)
        return albedo


    def interpolate_phase(self, reff, phase_type):
        """
        Interpolate the phase function over a grid.
        
        Parameters
        ----------
        reff: GridData 
            A GridData object containting the effective radii (micron) on a 3D grid.
        phase_type: 'Tabulated' or 'Grid'
            Return either a TabulatedPhase or GridPhase object.
        
        Returns
        -------
        phase: Phase object
            A Phase object containting the phase function legendre coeffiecients on a grid.
        """
        if phase_type == 'Tabulated':
            phase = self.get_tabulated_phase(reff)
        elif phase_type == 'Grid':
            phase = self.get_grid_phase(reff)
        else:
            raise AttributeError('Phase type not recognized')        
        return phase

    def interpolate_scattering_field(self, lwc, reff, phase_type='Tabulated'):
        """
        Interpolate the optical quantities (extinction, ssa, phase function) over a grid.
        All optical quantities are a function of the effective radius. The optical extinction also scales
        with the liquid water content.
    
        Parameters
        ----------
        lwc: GridData
            A GridData object containting the liquid water content (g/m^3) on a 3D grid
        reff: GridData 
            A GridData object containting the effective radii (micron) on a 3D grid.
        phase_type: 'Tabulated' or 'Grid'
            Return either a TabulatedPhase or GridPhase object.
            
        Returns
        -------
        extinction: GridData object
            A GridData3D object containting the extinction (1/km) on a 3D grid
        albedo: GridData object
            A GridData object containting the single scattering albedo unitless in range [0, 1] on a 3D grid
        phase: Phase object
            A Phase object containting the phase function legendre coeffiecients on a grid.
        """   
        extinction = self.interpolate_extinction(lwc, reff)
        albedo = self.interpolate_albedo(reff)
        phase = self.interpolate_phase(reff, phase_type)
        return extinction, albedo, phase
    
    
    def get_tabulated_phase(self, reff):
        """
        Interpolate the phase function for an effective radius grid.
    
        Parameters
        ----------
        reff: float
            The effective radius for which to retrieve the phase function [microns].
            
        Returns
        -------
        phase: TabulatedPhase
            A Phase object containting the phase function legendre coeffiecients as a table.
        """   
        maxleg = self.nleg.max()
        phase_table = self.legcoeff[:maxleg, :]
        index_data = self._legen_index_interpolator(reff.data).astype(np.int32)
        index = GridData(reff.grid, index_data)
        phase = TabulatedPhase(phase_table, index)        
        return phase    

    def get_grid_phase(self, reff):
        """
        Interpolate the phase function for an effective radius grid.
    
        Parameters
        ----------
        reff: float
            The effective radius for which to retrieve the phase function [microns].
            
        Returns
        -------
        phase: GridPhase
            A Phase object containting the phase function legendre coeffiecients on a 3D grid
        """   
        maxleg = self._nleg_interpolator(reff.data).max().astype(np.int)
        phase_data = self._legcoef_interpolator(reff.data)[:maxleg + 1,...]
        phase = GridPhase(reff.grid, phase_data) 
        return phase


    @property
    def reff(self):
        if hasattr(self, '_reff'):
            return self._reff
        else:
            print('Mie table was not computed or loaded')    

    @property
    def veff(self):
        if hasattr(self, '_alpha'):
            if self._distflag == 'G':
                return  1.0/(self._alpha+3.0) 
            elif self._distflag == 'L':
                return np.exp(self._alpha**2) - 1.0             
        else:
            print('Mie table was not computed or loaded')  

    @property
    def extinct(self):
        if hasattr(self, '_extinct'):
            return self._extinct
        else:
            print('Mie table was not computed or loaded') 
            
    @property
    def ssalb(self):
        if hasattr(self, '_ssalb'):
            return self._ssalb
        else:
            print('Mie table was not computed or loaded') 
            
    @property
    def nleg(self):
        if hasattr(self, '_nleg'):
            return self._nleg
        else:
            print('Mie table was not computed or loaded') 
            
    @property
    def legcoeff(self):
        if hasattr(self, '_legcoef'):
            return self._legcoef
        else:
            print('Mie table was not computed or loaded')
            
    @property
    def maxleg(self):
        if hasattr(self, '_maxleg'):
            return self._maxleg
        else:
            print('Mie table was not computed or loaded')            
            
    @property
    def distribution(self):
        if hasattr(self, '_distflag'):
            if self._distflag == 'G':
                return 'Gamma'
            elif self._distflag == 'L':
                return 'Lognormal'
        else:
            print('Mie table was not computed or loaded')            

         
class Rayleigh(object):
    """
    Rayleigh scattering for temperature profile.
    
    Description taken from cloudprp.f:
     Computes the molecular Rayleigh extinction profile EXTRAYL [/km]
     from the temperature profile TEMP [K] at ZLEVELS [km].  Assumes
     a linear lapse rate between levels to compute the pressure at
     each level.  The Rayleigh extinction is proportional to air
     density, with the coefficient RAYLCOEF in [K/(mb km)].
     
    Parameters
    ----------
    wavelength: float 
        The wavelength in [microns].
    temperature_profile: TemperatureProfile 
        A TemperatureProfile object containing temperatures and altitudes on a 1D grid.
        
    Notes
    -----
    The ssa of air is 1.0 and the phase function legendre coefficients are [1.0, 0.0, 0.5].
    """
    def __init__(self, wavelength, temperature_profile):
        self._wavelength = wavelength
        self._temperature_profile = temperature_profile
        self._raylcoeff = (2.97e-4) * wavelength**(-4.15 + 0.2 * wavelength)
        self._ssalb = np.array([1.0], dtype=np.float32)
        self._phase = np.array([1.0, 0.0, 0.5], dtype=np.float32)
       
        
    def get_scattering_field(self, grid, phase_type='Tabulated'):
        """
        Interpolate the rayleigh extinction over a given 1D grid (altitude).
         
        Parameters
        ----------
        grid: Grid  
            a Grid object containing the altitude grid points in [km].
        phase_type: 'Tabulated' or 'Grid'
            Return either a TabulatedPhase or GridPhase object.
            
        Returns
        -------
        extinction: GridData
            A GridData object containting the extinction (1/km) on a 1D grid
        albedo: GridData
            A GridData object containting the single scattering albedo unitless in range [0, 1] on a 3D grid
        phase: Phase 
            A Phase object containting the phase function legendre table and index grid.
        """
        temperature_profile = self.temperature_profile.resample(grid)

        extinction_profile = core.rayleigh_extinct(
            nzt=grid.nz,
            zlevels=grid.z,
            temp=temperature_profile,
            raylcoef=self.rayleigh_coefficient
        )
         
        extinction = GridData(grid, extinction_profile)
        albedo = GridData(grid, np.full(shape=(grid.nz,), fill_value=self.ssalb, dtype=np.float32))
        if phase_type == 'Tabulated':
            phase_indices = GridData(grid, np.ones(shape=(grid.nz,), dtype=np.int32))
            phase = TabulatedPhase(self.phase[:, np.newaxis], phase_indices)
        elif phase_type == 'Grid':
            phase = GridPhase(grid, np.tile(self.phase[:, np.newaxis, np.newaxis, np.newaxis], (1, 1, 1, grid.nz))) 
        else:
            raise AttributeError('Phase type not recognized')
      
        return extinction, albedo, phase
        
        
    @property
    def extinction_profile(self):
        if hasattr(self, '_extinction_profile'):
            return self._extinction_profile
        else:
            print('scattering profile was not computed') 
    
    @property
    def temperature_profile(self):
        return self._temperature_profile 
    
    @property
    def wavelength(self):
        return self._wavelength    
    
    @property
    def rayleigh_coefficient(self):
        return self._raylcoeff 
    
    @property
    def ssalb(self):
        return self._ssalb      
    
    @property
    def phase(self):
        return self._phase      