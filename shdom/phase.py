"""
This is a python wrapper to handle phase function related computations. 

It includes a wrapper for mie and rayleigh computations. 
The python code was created by Aviad Levis Technion Inst. of Technology February 2019.

Source Fortran files were created by Frank Evans University of Colorado May 2003.
"""
import core
import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator
import shdom


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


class LegendreTable(object):
    """TODO"""
    def __init__(self, table,  table_type='SCALAR'):
        self._table_type = table_type
        
        # If table is comprised of a single phase entery
        if ((table.ndim==1) and (table_type=='SCALAR')) or \
           ((table.ndim==2) and (table_type=='VECTOR')):
            table = table[...,np.newaxis]        

        self._data = np.array(table, dtype=np.float32)
        self._numphase = self.data.shape[-1]
        self._maxleg = self.data.shape[-2] - 1
        
        # Asymetry parameter is proportional to the legendre series first coefficient 
        if table_type=='SCALAR':
            self._maxasym = self.data[1].max() / 3.0
            self._nstleg = 1
            
        elif table_type=='VECTOR':
            self._maxasym = self.data[0, 1].max() / 3.0        
            self._nstleg = self.data.shape[0]


    def append(self, table):
        """
        TODO
        """
        assert self.table_type==table.table_type, 'Cannot append new table of type {} '\
               'to current table of type {}'.format(table.table_type, self.table_type)
        
        nleg = table.maxleg
        table_data = table.data
        curr_table_data = self.data
        if self.table_type == 'SCALAR':
            if self.maxleg > nleg :
                table_data = np.pad(table_data, ((0, self.maxleg - nleg), (0,0)), 'constant')            
            elif nleg > self.maxleg:
                curr_table_data = np.pad(curr_table_data, ((0, nleg - self.maxleg), (0,0)), 'constant')            
        
        elif self.table_type == 'VECTOR':
            if self.maxleg > nleg :
                table_data = np.pad(table_data, ((0,0), (0, self.maxleg - nleg), (0,0)), 'constant')            
            elif nleg > self.maxleg:
                curr_table_data = np.pad(curr_table_data, ((0,0), (0, nleg - self.maxleg), (0,0)), 'constant')  
        
        self._data = np.append(curr_table_data, table_data, axis=-1)
        self._maxasym = max(self.maxasym, table.maxasym)
        self._maxleg = max(self.maxleg, table.maxleg)
        self._numphase += table.numphase 
        

    def get_legenp(self, nleg):
        """ 
        TODO
        legenp is without the zero order term which is 1.0 for normalized phase function
        """
        legenp = self.data
        
        # Scalar (unpolarized) table
        if self.table_type=='SCALAR':
            legenp = legenp[1:]
            if nleg > self.maxleg:
                legenp = np.pad(legenp, ((0, nleg - self.maxleg), (0,0)), 'constant')
        
        if (self.table_type=='VECTOR') and (nleg>self.maxleg):
            legenp = np.pad(legenp, ((0,0), (0, nleg - self.maxleg), (0,0)) , 'constant')
        
        return legenp.ravel(order='F') 
    
    
    @property
    def data(self):
        return self._data
    
    @property
    def table_type(self):
        return self._table_type
    
    @property
    def numphase(self):
        return self._numphase
    
    @property
    def maxleg(self):
        return self._maxleg
    
    @property
    def maxasym(self):
        return self._maxasym
    
    @property
    def nstleg(self):
        return self._nstleg      
    

class GridPhase(object):
    """
    The GridPhase internally keeps a phase function table and a pointer array for each grid point.

    Parameters
    ----------
    legendre_table: shdom.LegendreTable
       An object containing the Legendre table.
    index: shdom.GridData object
       A shdom.GridData object with dtype=int. This is a pointer to the enteries in the legendre_table.
    """
    def __init__(self, legendre_table, index):
        self._legendre_table = legendre_table
        self._index = index
        self._grid = index.grid
        self._iphasep = index.data

    
    def resample(self, grid):
        """Resample data to a new shdom.Grid."""
        if self.grid == grid:
            grid_phase =  self
        else:
            index = self.index.resample(grid, method='nearest')
            index._data = index.data.clip(1)
            grid_phase = shdom.GridPhase(self.legendre_table, index)
        return grid_phase     
    
    
    @property
    def iphasep(self):
        return self._iphasep
    
    @property
    def grid(self):
        return self._grid
    
    @property
    def legendre_table(self):
        return self._legendre_table
    
    @property
    def index(self):
        return self._index
    


class MieMonodisperse(object):
    """
    Mie scattering for spherical particles.
    
    Parameters
    ----------
    particle_type: string
        Options are 'Water' or 'Aerosol'.
        
        
    Notes
    -----
    Aerosol particle type not supported.  
    """    
    def __init__(self, particle_type='Water'):
        
        if particle_type == 'Water':
            self._partype = 'W'
            self._rindex = 1.0
            self._pardens = 1.0
        else:
            raise NotImplementedError('Particle type {} not supported'.format(particle_type))
        
        self._table_type = None        
        self._partype = particle_type
        self._wavelen1 = None
        self._wavelen2 = None
        self._avgflag = None  
        self._deltawave = None
        self._wavelencen = None
        self._nsize = None
        self._radii = None
        self._maxleg = None
        self._nleg = None
        self._legcoef = None
        self._extinct = None
        self._scatter = None


    def set_wavelength_integration(self,
                                   wavelength_band,
                                   wavelength_averaging=False,
                                   wavelength_resolution=0.001):
        """
        Set the wavelength integration parameters to compute a scattering table.
        
        Parameters
        ----------
        wavelength_band: (float, float)
            (minimum, maximum) wavelength in microns. 
            This defines the spectral band over which to integrate, if both are equal monochrome quantities are computed. 
        wavelength_averaging: bool
            True - average scattering properties over the wavelength_band.
            False - scattering properties of the central wavelength. 
        wavelength_resolution: float
            The distance between two wavelength samples in the band. Used only if wavelength_averaging is True.
        """
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

        self._wavelencen = core.get_center_wavelen(
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2
        )
        
        self._rindex = core.get_refract_index(
            partype=self._partype, 
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2
        )
        
        
    def set_radius_integration(self,
                               minimum_effective_radius,
                               max_integration_radius):
        """
        Set the radius integration parameters to compute a scattering table.
        
        Parameters
        ----------
        minimum_effective_radius: float
            Minimum effective radius in microns. Used to compute minimum radius for integration.
        max_integration_radius: float
            Maximum radius in microns - cutoff for the size distribution integral
        """
        
        self._nsize = core.get_nsize(
            sretab=minimum_effective_radius,
            maxradius=max_integration_radius,
            wavelen=self._wavelencen
        )
        
        self._radii = core.get_sizes(
            sretab=minimum_effective_radius, 
            maxradius=max_integration_radius, 
            wavelen=self._wavelencen,
            nsize=self._nsize
        )
        
        # Calculate the maximum size parameter and the max number of Legendre terms
        if self._avgflag == 'A':
            xmax = 2 * np.pi * max_integration_radius / self._wavelen1
        else:
            xmax = 2 * np.pi * max_integration_radius / self._wavelencen
        self._maxleg = int(np.round(2.0 * (xmax + 4.0*xmax**0.3334 + 2.0)))


    def compute_table(self):
        """
        Compute monodisperse Mie scattering per radius.
        
        Notes
        -----
        This is a time consuming method.
        """    
        self._extinct, self._scatter, self._nleg, self._legcoef, self._table_type = \
            core.compute_mie_all_sizes(
                nsize=self._nsize,
                maxleg=self._maxleg,
                wavelen1=self._wavelen1,
                wavelen2=self._wavelen2,
                deltawave=self._deltawave,
                wavelencen=self._wavelencen,
                radii=self._radii,
                rindex=self._rindex,
                avgflag=self._avgflag,
                partype=self._partype
            )
        

    def write_table(self, file_path): 
        """
        Write a pre-computed table to <file_path>. 
    
        Parameters
        ----------
        file_path: str
            Path to file.
      

        Notes
        -----
        This function must be ran after pre-computing a scattering table with compute_table().
        """        
        print('Writing Mie monodisperse table to file: {}'.format(file_path))
        core.write_mono_table(
            mietabfile=file_path,
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2,
            deltawave=self._deltawave,
            partype=self._partype,
            pardens=self._pardens, 
            rindex=self._rindex,  
            radii=self._radii,
            extinct=self._extinct,
            scatter=self._scatter,
            nleg=self._nleg,
            maxleg=self._maxleg,
            legcoef=self._legcoef
        )      

    
    def read_table_header(self, file_path):
        """
        TODO
        """        
        wavelen1, wavelen2, deltawave = np.genfromtxt(file_path, max_rows=1, skip_header=1, usecols=(0, 1, 2), dtype=float)
        pardens = np.genfromtxt(file_path, max_rows=1, skip_header=2, usecols=(0), dtype=float)
        partype = np.asscalar(np.genfromtxt(file_path, max_rows=1, skip_header=2, usecols=(1), dtype=str))
        rindex  = np.complex(np.genfromtxt(file_path, max_rows=1, skip_header=3, usecols=(0), dtype=float), 
                             np.genfromtxt(file_path, max_rows=1, skip_header=3, usecols=(1), dtype=float))
        nsize = np.genfromtxt(file_path, max_rows=1, skip_header=4, usecols=(0), dtype=int)
        maxleg = np.genfromtxt(file_path, max_rows=1, skip_header=5, usecols=(0), dtype=int)

        return wavelen1, wavelen2, deltawave, pardens, partype, rindex, nsize, maxleg
        
        
    def read_table(self, file_path): 
        """
        Read a pre-computed table from <file_path>. 
    
        Parameters
        ----------
        file_path: str
            Path to file.
        """   
        
        print('Reading mie table from file: {}'.format(file_path))
        self._wavelen1, self._wavelen2, self._deltawave, self._pardens, \
            self._partype, self._rindex, self._nsize, self._maxleg = self.read_table_header(file_path)
    
        self._radii, self._extinct, self._scatter, self._nleg, self._legcoef, self._table_type = \
            core.read_mono_table(
                mietabfile=file_path,
                nrtab=self._nsize,
                maxleg=self._maxleg
            )
    

    @property
    def maxleg(self):
        return self._maxleg
    
    @property 
    def pardens(self):
        return self._pardens
    
    @property 
    def radii(self):
        return self._radii
    
    @property 
    def extinct(self):
        return self._extinct
    
    @property 
    def scatter(self):
        return self._scatter
    
    @property 
    def nleg(self):
        return self._nleg
    
    @property 
    def legcoef(self):
        return self._legcoef    
    
    @property 
    def table_type(self):
        return self._table_type     


class SizeDistribution(object):
    """
    Size distribution object to compute Polydisprese Mie scattering.
    
    Parameters
    ----------
    type: string
        Particle size-distribution type options are 'gamma' or 'log-normal'.

    Notes
    -----
    gamma:
      n(r) = a r^alpha exp(-b*r).
      r - droplet radius.
      a, b, alpha - gamma distribution parameters. 
    
    log-normal:
      n(r) = a/r exp( -[ln(r/r0)]^2 / (2*alpha^2) ).
      r0 - logarithmic mode.
      alpha - standard deviation of the log. 
      
    Modified gamma not supported yet.
    """
    def __init__(self, type='gamma'):
        
        if type == 'gamma':
            self._distflag = 'G'
        elif type == 'lognormal':
            self._distflag = 'L'
        else:
            raise NotImplementedError('Distribution type {} not supported'.format(distibution))
        
        self._alpha = None
        self._reff = None
        self._veff = None
        self._nd = None
        self._nd_interpolator = None
        
        # For modified gamma, currently unused
        self._gamma = 0.0 
    
    
    def set_parameters(self, reff, **kwargs):
        """
        Set size-distribution parameters. 
        
        Parameters
        ----------
        reff: np.array(shape=(nretab,), dtype=float)
            Effective radius array in microns
        veff: np.array(shape=(nvetab,), dtype=float), optional
            Effective variance array for the size distribution
        alpha: np.array(shape=(nvetab,), dtype=float), optional
            Shape parameter array for the size distribution
            
        Notes
        -----
        Either veff or alpha must be specified (not both).
        The reff, veff, N, alpha relationships are as follows:
          Gamma:
            N = a Gamma(alpha+1)/ b^(alpha+1) - number concentration.
            r_eff = (alpha+3)/b - effective radius.
            v_eff = 1/(alpha+3) - effective variance.
          
          Log-normal:
            N = sqrt(2*pi)*alpha*a - number concentration. 
            r_eff = r0*exp(2.5*alpha^2) - effective radius.
            v_eff = exp(alpha^2)-1  - effective variance.
        """
        self.reff = reff
        if kwargs.has_key('alpha'):
            self.alpha = kwargs['alpha']
        elif kwargs.has_key('veff'):
            self.veff = kwargs['veff']
        else:
            raise AttributeError('Neither veff or alpha were specified.')
        
    
    def compute_nd(self, radii, particle_density=1.0):
        """
        Compute the size-distribution and initialize interpolator.
        
        Parameters
        ----------
        radii: np.array(shape=(nsize,), dtype=float)
            The radii for which to compute Nd in microns
        particle_density: np.float,
            Particle density [g/cm^3] is used to normalize the size-distribution for 1.0 [g/cm^3]
            liquid water content or mass content. particle_density=1.0 is used for water particles.
        """
        self._radii = radii
        self._nsize = radii.shape[0]           
        self._pardens = particle_density        

        reff, alpha = np.meshgrid(self.reff, self.alpha)
        self._nd = core.make_multi_size_dist(
            distflag=self.distflag,
            pardens=self.pardens,
            nsize=self.nsize,
            radii=self.radii,
            reff=reff.ravel(),
            alpha=alpha.ravel(),
            gamma=self.gamma,
            ndist=reff.size)
        
        nd = self.nd.T.reshape((self.nretab, self.nvetab, self.nsize), order='F')
        self._nd_interpolator = RegularGridInterpolator(
                (self.reff, self.veff), nd, bounds_error=False, fill_value=0.0)
    
    
    def get_nd(self, reff, veff):
        return self._nd_interpolator((reff, veff)).T
    
    
    @property
    def radii(self):
        return self._radii
    
    @property
    def pardens(self):
        return self._pardens
    
    @property 
    def distflag(self):
        return self._distflag

    @property
    def nsize(self):
        return self._nsize

    @property
    def ndist(self):
        return self.nretab * self.nvetab
    
    @property
    def nretab(self):
        return self._nretab  
    
    @property
    def nvetab(self):
        return self._nvetab    
    
    @property 
    def alpha(self):
        return self._alpha
    
    @alpha.setter
    def alpha(self, val):
        self._alpha = val
        self._nvetab = 1 if np.isscalar(val) else val.shape[0]
        if self.distflag == 'G':
            self._veff = 1.0 / (val + 3.0)
        if self.distflag == 'L':
            self._veff = np.exp(val**2) - 1.0
        
    @property 
    def gamma(self):
        return self._gamma    
    
    @property 
    def reff(self):
        return self._reff   
    
    @reff.setter
    def reff(self, val):
        self._reff = val
        self._nretab = 1 if np.isscalar(val) else val.shape[0] 
    
    @property 
    def veff(self):
        return self._veff 
    
    @veff.setter
    def veff(self, val):
        self._veff = val
        self._nvetab = 1 if np.isscalar(val) else val.shape[0]
        if self.distflag == 'G':
            self._alpha = 1.0 / val - 3.0
        if self.distflag == 'L':
            self._alpha = np.sqrt(np.log(val + 1.0))

    @property
    def nd(self):
        return self._nd
    
    
    
class MiePolydisperse(object):
    """
    Polydisperse Mie scattering for spherical particles with size distribution.
    Scattering coefficients are averaged over a range of particle radii and wavelengths.
    
    Parameters
    ----------
    mono_disperse: shdom.MieMonodisperse
        A monodisperse object containing Mie scattering as a function of radii.
    size_distribution: shdom.SizeDistribution
        A size distribution object
        
    Notes
    -----
    See: 
       - notebooks/Make Mie Table.ipynb for usage examples. 
       - notebooks/Make Mie Table Polarized.ipynb for usage examples. 
    """    
    def __init__(self, mono_disperse=MieMonodisperse(), size_distribution=SizeDistribution()):
        self._table_type = None
        self.set_mono_disperse(mono_disperse)
        self.set_size_distribution(size_distribution)
        self._microphysical_medium = None
        self._extinct = None
        self._ssalb = None
        self._nleg = None
        self._maxleg = None
        self._legendre_table = None    


    def set_mono_disperse(self, mono_disperse):
        """TODO"""
        self._mono_disperse = mono_disperse
        self._table_type = mono_disperse.table_type


    def set_size_distribution(self, size_distribution):
        """TODO"""
        self._size_distribution = size_distribution
        
        
    def compute_table(self):
        """
        Compute a scattering table where for each effective radius:
          1. Extinction-cross section per 1 unit of mass content [g/m^3](liquid water content for water clouds)  
          2. Single scattering albedo, unitless in the range [0, 1].
          3. Legendre expansion coefficients of the normalized scattering phase function (first coefficient is always 1.0)
          4. Number of Legendre coefficients for each scattering phase function. 
        """
        if self.size_distribution.nd is None:
            self.size_distribution.compute_nd(self.mono_disperse.radii, 
                                              self.size_distribution.pardens)

        self._extinct, self._ssalb, self._nleg, self.legcoef = \
            core.get_poly_table(
                nd=self.size_distribution.nd,
                ndist=self.size_distribution.ndist,
                nsize=self.size_distribution.nsize,
                maxleg=self.mono_disperse.maxleg,
                nleg1=self.mono_disperse.nleg,
                extinct1=self.mono_disperse.extinct,
                scatter1=self.mono_disperse.scatter,
                legcoef1=self.mono_disperse.legcoef)
        self.init_intepolators()
        

    def get_legendre(self, reff, veff):
        """
        TODO
        """
        reff_index = find_nearest(self.size_distribution.reff, reff)
        veff_index = find_nearest(self.size_distribution.veff, veff)

        # Fortran style indexing
        index = veff_index*self.size_distribution.nvetab + reff_index
        return self.legcoef[..., index]


    def get_angular_scattering(self, reff, veff, angles, phase_element=1):
        """
        Transfrom the spectral representation into angular representation.

        Parameters
        ----------
        reff: float
            The radius for which to compute the angular scattering
        angles: np.array(dtype=float)
            An array of angles for which to compute the phase function
        phase_element: int (only used if compiled with polarization)
            An integer in the range [1,6] where:
            phase_element=1: P11
            phase_element=2: P22
            phase_element=3: P33
            phase_element=4: P44
            phase_element=5: P12
            phase_element=6: P34

        Returns
        -------
        phase: np.array(dtype=float, shape=(len(angles),))
            The phase element as a function of angles
        """
        reff_index = find_nearest(self.size_distribution.reff, reff)
        veff_index = find_nearest(self.size_distribution.veff, veff)
        # Fortran style indexing
        index = veff_index*self.size_distribution.nvetab + reff_index        

        phase = core.transform_leg_to_phase(
            maxleg=self.maxleg,
            nphasepol=6,
            pelem=phase_element,
            nleg=self.nleg[index],
            legcoef=self.legcoef[...,index],
            nangle=len(angles),
            angle=angles
        )
        return phase   


    def write_table(self, file_path): 
        """
        Write a pre-computed table to <file_path>. 
    
        Parameters
        ----------
        file_path: str
            Path to file.

        Notes
        -----
        This function must be ran after pre-computing a scattering table with compute_table().
        """        
        print('Writing mie table to file: {}'.format(file_path))
        core.write_poly_table(
            mietabfile=file_path,
            wavelen1=self.mono_disperse._wavelen1, 
            wavelen2=self.mono_disperse._wavelen2,
            deltawave=self.mono_disperse._deltawave,
            partype=self.mono_disperse._partype,
            pardens=self.mono_disperse._pardens, 
            rindex=self.mono_disperse._rindex,
            distflag=self.size_distribution.distflag,
            alpha=self.size_distribution.alpha, 
            nretab=self.size_distribution.nretab,
            nvetab=self.size_distribution.nvetab,           
            reff=self.size_distribution.reff,
            veff=self.size_distribution.veff,
            extinct=self.extinct,
            ssalb=self.ssalb,
            nleg=self.nleg,
            legcoef=self.legendre_table.data,
            ndist=self.size_distribution.ndist,
            maxleg=self.legendre_table.maxleg,
            gamma=self.size_distribution.gamma)
    
    
    def read_table_header(self, file_path):
        """
        TODO
        """
        wavelen1, wavelen2, deltawave = np.genfromtxt(file_path, max_rows=1, skip_header=1, usecols=(0, 1, 2), dtype=float)
        pardens = np.genfromtxt(file_path, max_rows=1, skip_header=2, usecols=(0), dtype=float)
        partype = np.asscalar(np.genfromtxt(file_path, max_rows=1, skip_header=2, usecols=(1), dtype=str))
        rindex  = np.complex(np.genfromtxt(file_path, max_rows=1, skip_header=3, usecols=(0), dtype=float), 
                             np.genfromtxt(file_path, max_rows=1, skip_header=3, usecols=(1), dtype=float))
        distribution = np.asscalar(np.genfromtxt(file_path, max_rows=1, skip_header=4, usecols=(0), dtype=str))
        if distribution == 'gamma':
            distflag = 'G'
        elif distribution == 'lognormal':
            distflag = 'L'
        else:
            raise NotImplementedError('Distribution type {} not supported'.format(distribution))
    
        nretab = np.genfromtxt(file_path, max_rows=1, skip_header=5, usecols=(0), dtype=int)
        nvetab = np.genfromtxt(file_path, max_rows=1, skip_header=6, usecols=(0), dtype=int)
        maxleg = np.genfromtxt(file_path, max_rows=1, skip_header=8, usecols=(0), dtype=int)
    
        return wavelen1, wavelen2, deltawave, pardens, partype, rindex, distflag, nretab, nvetab, maxleg        

    
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
        print('Reading mie table from file: {}'.format(file_path))
        self.mono_disperse._wavelen1, self.mono_disperse._wavelen2, self.mono_disperse._deltawave, \
            self.mono_disperse._pardens, self.mono_disperse._partype, self.mono_disperse._rindex, \
            self.mono_disperse._distflag, self.size_distribution._nretab, self.size_distribution._nvetab, \
            self._maxleg = self.read_table_header(file_path)

        self.size_distribution.reff, self.size_distribution.veff, self._extinct, self._ssalb, \
            self._nleg, self.legcoef, self._table_type = core.read_poly_table(
                mietabfile=file_path, 
                nretab=self.size_distribution.nretab, 
                nvetab=self.size_distribution.nvetab,
                ndist=self.size_distribution.ndist,
                maxleg=self.maxleg)           
        self.init_intepolators()


    def init_intepolators(self):
        """
        Initialize interpolators for the extinction, single scattering albedo and phase function.
    
        Notes
        -----
        This function is ran after pre-computing/loading a scattering tables.
        """        
        
        # Truncate the legcoef (or wigcoef if polarized)
        self._maxleg = int(self.nleg.max())
        if self.table_type == 'SCALAR':
            self._legcoef = self.legcoef[:self._maxleg+1, :]
        elif self.table_type == 'VECTOR':
            self._legcoef = self.legcoef[:, :self._maxleg+1, :]

        self._legendre_table = shdom.LegendreTable(self._legcoef, self.table_type)
        
        method = 'nearest' if self.size_distribution.ndist==1 else 'linear'
        
        reff, veff = self.size_distribution.reff, self.size_distribution.veff
        extinct = self.extinct.reshape((self.size_distribution.nretab, self.size_distribution.nvetab), order='F')
        ssalb = self.ssalb.reshape((self.size_distribution.nretab, self.size_distribution.nvetab), order='F')
        nleg = self.nleg.reshape((self.size_distribution.nretab, self.size_distribution.nvetab), order='F')
        legen_index = np.transpose(np.array(np.meshgrid(range(self.size_distribution.nretab),
                                                        range(self.size_distribution.nvetab), 
                                                        indexing='ij')), [1, 2, 0])
        
        self._ext_interpolator = RegularGridInterpolator((reff, veff), extinct, method=method, bounds_error=False, fill_value=0.0)
        self._ssalb_interpolator = RegularGridInterpolator((reff, veff), ssalb, method=method, bounds_error=False, fill_value=1.0)
        self._nleg_interpolator = RegularGridInterpolator((reff, veff), nleg, method='nearest', bounds_error=False, fill_value=0)
        self._legen_index_interpolator = RegularGridInterpolator((reff, veff), legen_index, method='nearest', bounds_error=False, fill_value=0)


    def get_extinction(self, lwc, reff, veff):
        """
        Interpolate the extinciton coefficient function over a grid.
    
        Parameters
        ----------
        lwc: shdom.GridData 
            A GridData object containting liquid water content (g/m^3) on a 3D grid.
        reff: shdom.GridData 
            A GridData object containting effective radii (micron) on a 3D grid.
        veff: shdom.GridData 
            A GridData object containting effective variances on a 3D grid.
        
        Returns
        -------
        extinction: shdom.GridData object
            A shdom.GridData3D object containting the extinction (1/km) on a 3D grid
        """
        grid = veff.grid + reff.grid
        data = self._ext_interpolator((reff.resample(grid).data, veff.resample(grid).data))
        extinction = lwc.resample(grid) * shdom.GridData(grid, data)     
        return extinction  
        
   
    def get_albedo(self, reff, veff):
        """
        Interpolate the single scattering albedo over a grid.
        
        Parameters
        ----------
        reff: shdom.GridData 
            A shdom.GridData object containting the effective radii (micron) on a 3D grid.
        veff: shdom.GridData 
            A GridData object containting effective variances on a 3D grid.
            
        Returns
        -------
        albedo: shdom.GridData object
            A shdom.GridData object containting the single scattering albedo unitless in range [0, 1] on a 3D grid
        """
        grid = veff.grid + reff.grid
        data = self._ssalb_interpolator((reff.resample(grid).data, veff.resample(grid).data))
        albedo = shdom.GridData(grid, data)
        return albedo


    def get_phase(self, reff, veff, squeeze_table=True):
        """
        Interpolate the phase function over a grid.
    
        Parameters
        ----------
        reff: shdom.GridData 
            A shdom.GridData object containting the effective radii (micron) on a 3D grid.
        veff: shdom.GridData 
            A GridData object containting effective variances on a 3D grid.
        squeeze_table: boolean
            True: return compact table containing the range of provided reff, veff. 
            False will return the current table.
            
        Returns
        -------
        phase: GridPhase
            A GridPhase object containting the phase function legendre coeffiecients as a table.
        """
        grid = veff.grid + reff.grid
        
        index = self._legen_index_interpolator((reff.resample(grid).data, veff.resample(grid).data)).astype(np.int32)
        nre, nve = self.size_distribution.nretab, self.size_distribution.nvetab 
        
        # Clip table to make it compact
        if squeeze_table:
            max_re_idx = min(nre-1, find_nearest(self.size_distribution.reff, reff.data.max()))
            max_ve_idx = min(nve-1, find_nearest(self.size_distribution.veff, veff.data.max()))
            min_re_idx = max(0, find_nearest(self.size_distribution.reff, reff.data[reff.data>0.0].min())-1)
            min_ve_idx = max(0, find_nearest(self.size_distribution.veff, veff.data[veff.data>0.0].min())-1)
            legcoef = self.legcoef_2d[..., min_re_idx:max_re_idx+1, min_ve_idx:max_ve_idx+1]
            nre, nve = legcoef.shape[-2:]
            legcoef = legcoef.reshape(self.legcoef.shape[:-1] + (-1,), order='F')
            index[...,0] -= min_re_idx
            index[...,1] -= min_ve_idx
            
        else:
            legcoef = self.legcoef

        legen_table = shdom.LegendreTable(legcoef, self.table_type)
        index = np.ravel_multi_index(np.rollaxis(index, axis=-1), dims=(nre, nve), order='F',  mode='clip') + 1
        legen_index = shdom.GridData(grid, index.astype(np.int32))
        phase = GridPhase(legen_table, legen_index)        
        return phase

        
    def set_microphysical_medium(self, microphysical_medium):
        max_re = microphysical_medium.reff.max_value
        max_ve = microphysical_medium.veff.max_value
        min_re = microphysical_medium.reff.data[microphysical_medium.reff.data>0.0].min()
        min_ve = microphysical_medium.veff.data[microphysical_medium.veff.data>0.0].min()
        
        assert  max_re < self.size_distribution.reff.max()+1e-3, \
               'Maximum medium effective radius [{:2.2f}] is larger than the pre-computed table maximum radius [{:2.2f}]. ' \
               'Recompute Mie table with larger maximum radius.'.format(max_re, self.size_distribution.reff.max())
        assert  min_re > self.size_distribution.reff.min()-1e-3, \
               'Minimum medium effective radius [{:2.2f}] is smaller than the pre-computed table minimum radius [{:2.2f}]. ' \
               'Recompute Mie table with smaller minimum radius.'.format(min_re, self.size_distribution.reff.min())
        assert  max_ve < self.size_distribution.veff.max()+1e-3, \
               'Maximum medium effective variance [{:2.2f}] is larger than the pre-computed table maximum variance [{:2.2f}]. ' \
               'Recompute Mie table with larger maximum variance.'.format(max_ve, self.size_distribution.veff.max())
        assert  min_ve > self.size_distribution.veff.min()-1e-3, \
               'Minimum medium effective variance [{:2.2f}] is smaller than the pre-computed table minimum variance [{:2.2f}]. ' \
               'Recompute Mie table with smaller minimum variance.'.format(min_ve, self.size_distribution.veff.min())
        
        self._microphysical_medium = microphysical_medium


    def get_scatterer(self, squeeze_table=True):
        """
        Interpolate optical quantities (extinction, ssa, phase function) over a grid.
        All optical quantities are a function of the effective radius.
        The optical extinction also scales with liquid water content.
    
        Returns
        -------
        medium: shdom.OpticalMedium object
            An OpticalMedium object encapsulating the the extinction, single scattering albedo and phase function on a 3D grid
        """   
        lwc = self.microphysical_medium.lwc
        reff = self.microphysical_medium.reff
        veff = self.microphysical_medium.veff
        
        extinction = self.get_extinction(lwc, reff, veff)
        albedo = self.get_albedo(reff, veff)
        phase = self.get_phase(reff, veff, squeeze_table)
        return shdom.Scatterer(extinction, albedo, phase)
    
    
    @property
    def nleg(self):
        return self._nleg
    
    @property
    def maxleg(self):
        return self._maxleg      
    
    @property
    def legcoef(self):
        return self._legcoef 
    
    @legcoef.setter
    def legcoef(self, val):
        self._legcoef = val
        self._legcoef_2d = val.reshape(val.shape[:-1] + (self.size_distribution.nretab, self.size_distribution.nvetab), order='F')
    
    @property
    def legcoef_2d(self):
        return self._legcoef_2d    
    
    @property
    def legendre_table(self):
        return self._legendre_table 
            
    @property
    def extinct(self):
        return self._extinct 
    
    @property
    def ssalb(self):
        return self._ssalb 
    
    @property
    def table_type(self):
        return self._table_type
    
    @property
    def mono_disperse(self):
        return self._mono_disperse
    
    @property
    def size_distribution(self):
        return self._size_distribution
    
    @property
    def microphysical_medium(self):
        return self._microphysical_medium    
    
    
            
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
    surface_pressure: float
        Surface pressure in units [mb]
        
        
    Notes
    -----
    The ssa of air is 1.0 and the phase function legendre coefficients are [1.0, 0.0, 0.5].
    """
    def __init__(self, wavelength):
        self._wavelength = wavelength
        self._temperature_profile = None
        self._surface_pressure = None
        self._raylcoeff = None
      
    def get_extinction(self, grid):
        ext_profile = core.rayleigh_extinct(
            nzt=grid.nz,
            zlevels=grid.z,
            temp=self.temperature_profile.data,
            raysfcpres=self.surface_pressure,
            raylcoef=self.raylcoef
        )        
        return shdom.GridData(grid, ext_profile)
    
        
    def get_albedo(self, grid):
        """TODO"""
        return shdom.GridData(grid, data=np.full(shape=grid.shape, fill_value=1.0, dtype=np.float32))
        

    def get_phase(self, grid):
        """TODO"""
        index = shdom.GridData(grid, data=np.ones(shape=grid.shape, dtype=np.int32))
        table, table_type = core.rayleigh_phase_function(wavelen=self._wavelength)
        table = LegendreTable(table.astype(np.float32), table_type) 
        self._phase = GridPhase(table, index)        
        return GridPhase(table, index)
    
    
    def set_profile(self, temperature_profile, surface_pressure=1013.0):
        """        
        Set the tempertature profile and surface pressure and initialize optical fields.
        
        Parameters
        ----------
        temperature_profile: shdom.GridData 
            a shdom.GridData object containing the altitude grid points in [km] and temperatures in [Kelvin].
        surface_pressure: float
            Surface pressure in units [mb] (default value is 1013.0)
        """ 
        self._temperature_profile = temperature_profile
        self._surface_pressure = surface_pressure  
        self._grid = temperature_profile.grid
        
        # Use the parameterization of Bodhaine et al. (1999) eq 30 for tau_R at 
        # sea level and convert to Rayleigh density coefficient:
        #   k = 0.03370*(p_sfc/1013.25)*tau_R  for k in K/(mb km)
        self._raylcoef = 0.03370 * (surface_pressure/1013.25) * 0.0021520 * \
            (1.0455996 - 341.29061/self.wavelength**2 - 0.90230850*self.wavelength**2) / (1 + 0.0027059889/self.wavelength**2 - 85.968563*self.wavelength**2)
        
        self._extinction = self.get_extinction(self.grid)
        self._albedo = self.get_albedo(self.grid)
        self._phase = self.get_phase(self.grid)
        
        
    def get_scatterer(self):
        """
        TODO
    
        Returns
        -------
        rayleigh: shdom.Scatterer
            A shdom.Scatterer object containting the Rayleigh optical properties on a 1D grid.
        """
        rayleigh = shdom.Scatterer(self.extinction, self.albedo, self.phase)
        return rayleigh


    @property
    def wavelength(self):
        return self._wavelength    
    
    @property
    def raylcoef(self):
        return self._raylcoef
    
    @property
    def extinction(self):
        return self._extinction   
    
    @property
    def albedo(self):
        return self._albedo   
    
    @property
    def phase(self):
        return self._phase 
    
    @property
    def temperature_profile(self):
        return self._temperature_profile 
    
    @property
    def grid(self):
        return self._grid  
    
    @property
    def surface_pressure(self):
        return self._surface_pressure  
    
    
    
class RayleighPolarized(Rayleigh):
    """TODO"""
    def __init__(self, wavelength):
        super(RayleighPolarized, self).__init__(wavelength)
                
    def init_phase(self):
        """TODO"""
        self._phase = core.rayleigh_phase_function(wavelen=self._wavelength).astype(np.float32)   

    def get_scattering_field(self, grid, phase_type='Tabulated'):
        """
        TODO
         
        Parameters
        ----------
        grid: Grid  
            a Grid object containing the altitude grid points in [km].
        phase_type: 'Tabulated' or 'Grid'
            Return either a TabulatedPhase or GridPhase object.
            
        Returns
        -------
        extinction: shdom.GridData
            A shdom.GridData object containting the extinction (1/km) on a 1D grid
        albedo: shdom.GridData
            A shdom.GridData object containting the single scattering albedo unitless in range [0, 1] on a 3D grid
        phase: Phase 
            A Phase object containting the phase function legendre table and index grid.
        """
        
        temperature_profile = self.temperature_profile.resample(grid)
        
        extinction_profile = core.rayleigh_extinct(
            nzt=grid.nz,
            zlevels=grid.z,
            temp=temperature_profile,
            raysfcpres=self._surface_pressure,
            raylcoef=self.rayleigh_coefficient
        )
        
        
        if grid.type == '1D':
            extinction = shdom.GridData(grid, extinction_profile)

        elif grid.type == '3D':
            extinction = shdom.GridData(grid, np.tile(extinction_profile, (grid.nx, grid.ny,1)))
    
        albedo = shdom.GridData(grid, np.full(shape=grid.shape, fill_value=self.ssalb, dtype=np.float32))
        
        if phase_type == 'Tabulated':
            phase_indices = shdom.GridData(grid, np.ones(shape=grid.shape, dtype=np.int32))
            phase = TabulatedPhaseMatrix(self.phase[...,np.newaxis], phase_indices)
        elif phase_type == 'Grid':
            phase = GridPhaseMatrix(grid, np.tile(self.phase[..., np.newaxis, np.newaxis, np.newaxis], (1, 1, 1, grid.nz))) 
        else:
            raise AttributeError('Phase type not recognized')
      
        return extinction, albedo, phase    