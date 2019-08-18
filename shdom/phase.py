"""
This is a python wrapper to handle phase function related computations. 

It includes a wrapper for mie and rayleigh computations. 
The python code was created by Aviad Levis Technion Inst. of Technology February 2019.

Source Fortran files were created by Frank Evans University of Colorado May 2003.
"""
import numpy as np
import shdom
from shdom import core, find_nearest
from scipy.interpolate import RegularGridInterpolator


class LegendreTable(object):
    """
    A LegendreTable object to encapsulate phase function tables.

    Parameters
    ----------
    table: np.array(dtype=no.float32)
        A table containing the phase function enteries.
        For a SCALAR table the first dimension is the legendre coeffiecients and the second dimension is the talbe entries.
        For a VECTOR table the first dimension is the stokes components (6), the second dimension is the legendre coeffiecients and the third dimension is the talbe entries.
    """
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
        Appending another LegenreTable to the current.

        Parameters
        ----------
        table: shdom.LegendreTable
            A LegendreTable to append to the current table.
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
        Retrieves a flattened table to be used by the shdom.RteSolver.

        Parameters
        ----------
        nleg: int
            The number of legendre coefficients. If nleg > self.maxleg then the table is padded with zeros to match nleg.

        Returns
        -------
        legenp: np.array(dtype=np.float32)
            A flattened legendre table

        Notes
        -----
        For a SCALAR table the zero order term, which is 1.0 for normalized phase function, is removed.
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
    The GridPhase internally keeps a legendre table and a pointer array for each grid point.

    Parameters
    ----------
    legendre_table: shdom.LegendreTable
       An object encapsulating the Legendre table.
    index: shdom.GridData object
       A shdom.GridData object with dtype=int. This is a pointer to the enteries in the legendre_table.
    """
    def __init__(self, legendre_table, index):
        self._legendre_table = legendre_table
        self._index = index
        self._grid = index.grid
        self._iphasep = index.data
    
    def resample(self, grid):
        """
        The resample method resamples the GridPhase (i.e. the index GridData).

        Parameters
        ----------
        grid: shdom.Grid
            The new grid to which the data will be resampled

        Returns
        -------
        GridPhase: shdom.GridPhase
            A GridPhase resampled onto the input grid
        """
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
    Mie monodisperse scattering for spherical particles.
    
    Parameters
    ----------
    particle_type: string
        Options are 'Water' or 'Aerosol'. Default is 'Water'.
        
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
        self._extinct, self._scatter, self._nleg, self._legcoef, table_type = \
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
        self._table_type = table_type.decode()

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
        Read the table header from <file_path>.

        Parameters
        ----------
        file_path: str
            Path to file.

        Notes
        -----
        The header contains the following information:
            wavelen1, wavelen2, deltawave, pardens, partype, rindex, nsize, maxleg
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
        
        self._wavelencen = core.get_center_wavelen(
            wavelen1=self._wavelen1, 
            wavelen2=self._wavelen2
        )        

        self._radii, self._extinct, self._scatter, self._nleg, self._legcoef, table_type = \
            core.read_mono_table(
                mietabfile=file_path,
                nrtab=self._nsize,
                maxleg=self._maxleg
            )
        self._table_type = table_type.decode()

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
    Size distribution object to compute Polydisprese scattering.
    
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
    
    lognormal:
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
        if 'alpha' in kwargs:
            self.alpha = kwargs['alpha']
        elif 'veff' in kwargs:
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
    Scattering coefficients are averaged over a range of particle radii and wavelengths (optionaly) .
    
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
        self._extinct = None
        self._ssalb = None
        self._nleg = None
        self._maxleg = None
        self._legendre_table = None 

    def set_mono_disperse(self, mono_disperse):
        """
        Set the monodisperse Mie scattering.

        Parameters
        ----------
        mono_disperse: shdom.MieMonodisperse
            The monodisperse Mie scattering.
        """
        self._mono_disperse = mono_disperse
        self._table_type = mono_disperse.table_type

    def set_size_distribution(self, size_distribution):
        """
        Set the size-distribution.

        Parameters
        ----------
        size_distribution: shdom.SizeDistribution
            The size-distribution.
        """
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
            self.size_distribution.compute_nd(self.mono_disperse.radii, self.mono_disperse.pardens)

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
        Retrieve the Legendre coeffiecients at a specific effective radius and variance.

        Parameters
        ----------
        reff: float
            The effective radius [microns]
        veff: float
            The effective variance

        Returns
        -------
        legcoef: np.array(dtype=np.float32)
            The Legendre coefficients

        Notes
        -----
        The legcoef of the nearest reff, veff that is found in the table is retrieved.
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
        Read the table header from <file_path>.

        Parameters
        ----------
        file_path: str
            Path to file.

        Notes
        -----
        The header contains the following information:
            wavelen1, wavelen2, deltawave, pardens, partype, rindex, distflag, nretab, nvetab, maxleg
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
        """   
        print('Reading mie table from file: {}'.format(file_path))
        self._mono_disperse._wavelen1, self._mono_disperse._wavelen2, self._mono_disperse._deltawave, \
            self._mono_disperse._pardens, self._mono_disperse._partype, self._mono_disperse._rindex, \
            self._mono_disperse._distflag, self.size_distribution._nretab, self.size_distribution._nvetab, \
            self._maxleg = self.read_table_header(file_path)
        
        self._mono_disperse._wavelencen = core.get_center_wavelen(
            wavelen1=self.mono_disperse._wavelen1, 
            wavelen2=self.mono_disperse._wavelen2
        )
        
        self.size_distribution.reff, self.size_distribution.veff, self._extinct, self._ssalb, \
            self._nleg, self.legcoef, table_type = core.read_poly_table(
                mietabfile=file_path, 
                nretab=self.size_distribution.nretab, 
                nvetab=self.size_distribution.nvetab,
                ndist=self.size_distribution.ndist,
                maxleg=self.maxleg) 
        self._table_type = table_type.decode()
        self.init_intepolators()

    def init_intepolators(self):
        """
        Initialize interpolators for the extinction, single scattering albedo and phase function.
    
        Notes
        -----
        This function is ran after pre-computing/loading a table.
        """        
        
        # Truncate the legcoef (or wigcoef if polarized)
        self._maxleg = int(self.nleg.max())
        if self.table_type == 'SCALAR':
            self.legcoef = self.legcoef[:self.maxleg+1, :]
        elif self.table_type == 'VECTOR':
            self.legcoef = self.legcoef[:, :self.maxleg+1, :]

        self._legendre_table = shdom.LegendreTable(self.legcoef, self.table_type)
        
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
        Retrieve the extinction coefficient over a grid.
    
        Parameters
        ----------
        lwc: shdom.GridData 
            A GridData object containing liquid water content (g/m^3) on a grid.
        reff: shdom.GridData 
            A GridData object containing effective radii (micron) on a grid.
        veff: shdom.GridData 
            A GridData object containing effective variances on a grid.
        
        Returns
        -------
        extinction: shdom.GridData object
            A GridData object containing the extinction (1/km) on a grid
        """
        grid = lwc.grid + veff.grid + reff.grid
        data = self._ext_interpolator((reff.resample(grid).data, veff.resample(grid).data))
        extinction = lwc.resample(grid) * shdom.GridData(grid, data)     
        return extinction  

    def get_albedo(self, reff, veff):
        """
        Interpolate the single scattering albedo over a grid.
        
        Parameters
        ----------
        reff: shdom.GridData 
            A shdom.GridData object containing the effective radii (micron) on a 3D grid.
        veff: shdom.GridData 
            A GridData object containing effective variances on a 3D grid.
            
        Returns
        -------
        albedo: shdom.GridData object
            A shdom.GridData object containing the single scattering albedo unitless in range [0, 1] on a 3D grid
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
            A shdom.GridData object containing the effective radii (micron) on a 3D grid.
        veff: shdom.GridData 
            A GridData object containing effective variances on a 3D grid.
        squeeze_table: boolean
            True: return compact table containing the range of provided reff, veff. 
            False will return the current table.
            
        Returns
        -------
        phase: GridPhase
            A GridPhase object containing the phase function legendre coeffiecients as a table.
        """
        grid = veff.grid + reff.grid
        
        index = self._legen_index_interpolator((reff.resample(grid).data, veff.resample(grid).data)).astype(np.int32)
        nre, nve = self.size_distribution.nretab, self.size_distribution.nvetab 
        
        # Clip table to make it compact
        if squeeze_table:
            max_re_idx = min(nre-1, find_nearest(self.size_distribution.reff, reff.data.max()))
            max_ve_idx = min(nve-1, find_nearest(self.size_distribution.veff, veff.data.max()))
            min_re_idx = max(0, find_nearest(self.size_distribution.reff, reff.data[reff.data>0.0].min()))
            min_ve_idx = max(0, find_nearest(self.size_distribution.veff, veff.data[veff.data>0.0].min()))
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
    def wavelength(self):
        if self.mono_disperse._wavelencen is None:
            return None
        else:
            return round(self.mono_disperse._wavelencen, 3)


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
    wavelength: float or list of floats
        The wavelengths in [microns].
    temperature_profile: TemperatureProfile 
        A TemperatureProfile object containing temperatures and altitudes on a 1D grid.
    surface_pressure: float
        Surface pressure in units [mb]

    Notes
    -----
    The ssa of air is 1.0 and the phase function legendre coefficients are [1.0, 0.0, 0.5].
    """
    def __init__(self, wavelength):
        if isinstance(wavelength, list):
            self._wavelength = wavelength
        else:
            self._wavelength = [wavelength]
        self._num_wavelengths = len(self._wavelength)
        self._temperature_profile = None
        self._surface_pressure = None
        self._raylcoeff = None

    def get_extinction(self, grid):
        """
        Retrieve the Rayleigh extinction profile (as a function of altitude).

        Parameters
        ----------
        grid: shdom.Grid
            The new grid to which the data will be resampled

        Returns
        -------
        extinction_profile: list of shdom.GridData
            A list of GridData objects containing the extinction on a 1D grid.
            The length of the list is the number of wavelengths.
        """
        ext_profile = [core.rayleigh_extinct(
            nzt=grid.nz,
            zlevels=grid.z,
            temp=self.temperature_profile.data,
            raysfcpres=self.surface_pressure,
            raylcoef=raylcoef
        )  for raylcoef in self._raylcoef]

        return [shdom.GridData(grid, ext) for ext in ext_profile]
        
    def get_albedo(self, grid):
        """
        Retrieve the Rayleigh single scattering albedo.

        Parameters
        ----------
        grid: shdom.Grid
            The new grid to which the data will be resampled

        Returns
        -------
        albedo: list of shdom.GridData
            A list of GridData objects containing the single scattering albedo [0,1] on a grid.
            The length of the list is the number of wavelengths.
        """
        return [shdom.GridData(grid, data=np.full(shape=grid.shape, fill_value=1.0, dtype=np.float32)) for wavelength in self._wavelength]

    def get_phase(self, grid):
        """
        Retrieve the Rayleigh phase function.

        Parameters
        ----------
        grid: shdom.Grid
            The new grid to which the data will be resampled

        Returns
        -------
        phase: list of shdom.GridPhase
            A list of GridPhase objects containing the phase function on a grid.
            The length of the list is the number of wavelengths.
        """
        phase = []
        for wavelength in self._wavelength:
            index = shdom.GridData(grid, data=np.ones(shape=grid.shape, dtype=np.int32))
            table, table_type = core.rayleigh_phase_function(wavelen=wavelength)
            table = LegendreTable(table.astype(np.float32), table_type.decode())
            phase.append(GridPhase(table, index))
        return phase
    
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
        self._raylcoef = [
            0.03370 * (surface_pressure/1013.25) * 0.0021520 * (1.0455996 - 341.29061/wl**2 - 0.90230850*wl**2) / (1 + 0.0027059889/wl**2 - 85.968563*wl**2) for wl in self._wavelength
        ]
        self._extinction = self.get_extinction(self.grid)
        self._albedo = self.get_albedo(self.grid)
        self._phase = self.get_phase(self.grid)

    def get_scatterer(self):
        """
        Retrieve a Scatterer from the medium.
    
        Returns
        -------
        scatterer: shdom.Scatterer
            A Scatterer object containing the Rayleigh optical properties on a 1D grid.

        Notes
        -----
        For a single band a shdom.OpticalScatterer is returned and for multiple wavelengths a shdom.MultispectralScatterer object.
        """
        scatterer_list = [
            shdom.OpticalScatterer(wavelength, extinction, albedo, phase) for \
            wavelength, extinction, albedo, phase in zip(self._wavelength, self._extinction, self._albedo, self._phase)
        ]
        if self.num_wavelengths == 1:
            scatterer = scatterer_list[0]
        else:
            scatterer = shdom.MultispectralScatterer(scatterer_list)
        return scatterer

    @property
    def num_wavelengths(self):
        return self._num_wavelengths

    @property
    def wavelength(self):
        if self.num_wavelengths == 1:
            return self._wavelength[0]
        else:
            return self._wavelength

    @property
    def raylcoef(self):
        if self.num_wavelengths == 1:
            return self._raylcoef[0]
        else:
            return self._raylcoef
    
    @property
    def extinction(self):
        if self.num_wavelengths == 1:
            return self._extinction[0]
        else:
            return self._extinction
    
    @property
    def albedo(self):
        if self.num_wavelengths == 1:
            return self._albedo[0]
        else:
            return self._albedo
    
    @property
    def phase(self):
        if self.num_wavelengths == 1:
            return self._phase[0]
        else:
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