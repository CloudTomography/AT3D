"""
This is a python wrapper to handle phase function related computations.

It includes a wrapper for mie and rayleigh computations.
The python code was created by Aviad Levis Technion Inst. of Technology February 2019.

Source Fortran files were created by Frank Evans University of Colorado May 2003.
"""
import numpy as np
import shdom
import xarray as xr
import pandas as pd
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
        table.pad(self.maxleg)
        self.pad(table.maxleg)
        self._data = np.append(self.data, table.data, axis=-1)
        self._maxasym = max(self.maxasym, table.maxasym)
        self._maxleg = max(self.maxleg, table.maxleg)
        self._numphase += table.numphase

    def pad(self, nleg):
        """
        Pads the table to a nleg number of coeffcients.

        Parameters
        ----------
        nleg: int
            The number of legendre coefficients. If nleg > self.maxleg then the table is padded with zeros to match nleg.
        """
        if self.table_type == 'SCALAR':
            if nleg > self.maxleg:
                self._data = np.pad(self.data, ((0, nleg - self.maxleg), (0, 0)), 'constant')

        if (self.table_type == 'VECTOR') and (nleg > self.maxleg):
            self._data = np.pad(self.data, ((0, 0), (0, nleg - self.maxleg), (0, 0)), 'constant')

    def get_legenp(self, num_stokes):
        """
        Retrieves a flattened table to be used by the shdom.RteSolver.

        Returns
        -------
        legenp: np.array(dtype=np.float32)
            A flattened legendre table

        Notes
        -----
        For a SCALAR table the zero order term, which is 1.0 for normalized phase function, is removed.
        """
        legenp = self.data
        if self.table_type == 'SCALAR':
            assert num_stokes ==1, 'A SCALAR LegendreTable cannot be used for polarized RTE'
            legenp = legenp[1:]
        elif self.table_type == 'VECTOR':
            if num_stokes == 1:
                legenp = legenp[0]
        return legenp.ravel(order='F').astype(np.float32)

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

class Mie(object):
    """
    Mie monodisperse scattering for spherical particles.
    This object computes Legendre coefficients and outputs Mie tables as xarrays.

    Parameters
    ----------
    particle_type: string
        Options are 'Water' or 'Aerosol'.
    wavelength_band: (float, float)
        (minimum, maximum) wavelength in microns.
        This defines the spectral band over which to integrate, if both are equal monochrome quantities are computed.
    minimum_effective_radius: float
        Minimum effective radius in microns. Used to compute minimum radius for integration.
    max_integration_radius: float
        Maximum radius in microns - cutoff for the size distribution integral
    wavelength_averaging: bool
        True - average scattering properties over the wavelength_band.
        False - scattering properties of the central wavelength.
    wavelength_resolution: float
        The distance between two wavelength samples in the band. Used only if wavelength_averaging is True.
    refractive_index: complex number or path to refractive index table (CSV)
        The refractive index should have a negative imaginary part. ri = n - ik
    """

    def __init__(self, particle_type, wavelength_band, minimum_effective_radius=4.0, max_integration_radius=65.0,
                 wavelength_averaging=False, wavelength_resolution=0.001, refractive_index=None):

        self._partype = particle_type

        # Set wavelength integration parameters
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

        # Set particle type properties
        if (particle_type == 'Water'):
            self._rindex = core.get_refract_index(
                partype=self._partype,
                wavelen1=self._wavelen1,
                wavelen2=self._wavelen2
            )

        elif (particle_type == 'Aerosol'):
            self._partype = 'A'
            assert refractive_index is not None, "Refractive index is not specified. \
            This could be a complex number or string path to csv (a function of wavelength)"
            if isinstance(refractive_index, str):
                rindex_df = pd.read_csv('../ancillary_data/dust_volz_1972.ri', comment='#', sep=' ',
                                        names=['wavelength', 'n', 'k'], index_col='wavelength')
                refractive_index = xr.Dataset.from_dataframe(rindex_df).interp({'wavelength': self._wavelencen})
                self._rindex = np.complex(refractive_index['n'], - refractive_index['k'])
            else:
                self._rindex = refractive_index

        else:
            raise AttributeError('Particle type note implemented')

        # Set radius integration parameters
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
        self._maxleg = int(np.round(2.0 * (xmax + 4.0 * xmax ** 0.3334 + 2.0)))

    def compute_table(self):
        """
        Compute monodisperse Mie scattering per radius.

        Notes
        -----
        This is a time consuming method.
        """
        extinct, scatter, nleg, legcoef, table_type = \
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
                partype=self._partype[0]
            )

        table = xr.Dataset(
            data_vars={
                'extinction': (['radius'], extinct),
                'scatter': (['radius'], scatter),
                'nleg': (['radius'], nleg),
                'legendre': (['stokes_index', 'legendre_index', 'radius'], legcoef)
            },
            coords={'radius': self._radii},
            attrs={
                'Particle type': self._partype,
                'Refractive index': self._rindex,
                'Table type': table_type.decode(),
                'units': ['Radius [micron]'],
                'Wavelength band': '[{}, {}] micron'.format(self._wavelen1, self._wavelen2),
                'Wavelength center': '{} micron'.format(self._wavelencen),
                'Wavlegnth averging': self._avgflag,
                'Wavelegnth resolution': self._deltawave,
                'Maximum legendre': self._maxleg
            },
        )
        return table


class SizeDistribution(object):

    def __init__(self, reff, distribution_type='gamma',particle_density=1.0,
                veff=None,alpha=None):

        self.distribution_type = distribution_type

        if distribution_type == 'gamma':
            self._distflag = 'G'
        elif distribution_type == 'lognormal':
            self._distflag = 'L'
        else:
            raise NotImplementedError('Distribution type {} not supported'.format(type))

        self.gamma=0.0

        self.reff = np.atleast_1d(reff)
        self._nretab = self.reff.shape[0]

        assert (veff is not None) or (alpha is not None), 'One of veff and alpha must be specified.'

        if alpha is not None:
            self.alpha = np.atleast_1d(alpha)
            self._nvetab = self.alpha.shape[0]
            if self._distflag == 'G':
                self.veff = 1.0 / (alpha + 3.0)
            if self._distflag == 'L':
                self.veff = np.exp(alpha**2) - 1.0
        elif veff is not None:

            self.veff = np.atleast_1d(veff)
            self._nvetab = self.veff.shape[0]
            if self._distflag == 'G':
                self.alpha = 1.0 / veff - 3.0
            if self._distflag == 'L':
                self.alpha = np.sqrt(np.log(veff + 1.0))

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
        self.radii = radii
        self._nsize = radii.shape[0]
        self._pardens = particle_density

        reff, alpha, gamma = np.meshgrid(self.reff,self.alpha,self.gamma)

        nd = core.make_multi_size_dist(
                    distflag=self._distflag,
                    pardens=self._pardens,
                    nsize=self._nsize,
                    radii=self.radii,
                    reff=reff.ravel(),
                    alpha=alpha.ravel(),
                    gamma=gamma.ravel(),
                    ndist=reff.size)
        nd = nd.T.reshape((self._nretab, self._nvetab, self._nsize), order='F')


        dataset = xr.Dataset(
                data_vars = {
                    'numberdensity': (['reff','veff','radii'], nd)

                },
                coords={'reff': self.reff,
                       'veff': self.veff,
                       'radii': self.radii},
                attrs={
                    'distribution_type': self.distribution_type,
                    'particle_density': particle_density
                })
        return dataset

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
    def __init__(self, mono_disperse=None, size_distribution=None):
        self._table_type = None
        mono_disperse = shdom.MieMonodisperse() if mono_disperse is None else mono_disperse
        size_distribution = shdom.SizeDistribution() if size_distribution is None else size_distribution
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
    wavelengths: float or list / numpy array.
        The wavelengths in [microns].
    """

    def __init__(self, wavelengths):
        self._wavelengths = np.atleast_1d(wavelengths)

    def compute_table(self):
        """
        Retrieve the Rayleigh phase function (Legendre table).

        Parameters
        ----------
        grid: shdom.Grid
            The new grid to which the data will be resampled

        Returns
        -------
        table: xarray.Dataset
            A dataset of Legendre coeffiecients for each wavelength specified.
        """
        legcoefs, table_types = zip(*[core.rayleigh_phase_function(wavelen=wavelen) for wavelen in self.wavelengths])
        arrs = [xr.DataArray(
            name='{:1.3} micron'.format(wavelength),
            data=legcoef, dims=['stokes_index', 'legendre_index'],
            attrs={'Table type': table_type.decode(), 'Wavelength': '{} [micron]'.format(wavelength)}
        ) for legcoef, table_type, wavelength in zip(legcoefs, table_types, self.wavelengths)]
        table = xr.merge(arrs)
        return table

    def compute_extinction(self, temperature_profile, surface_pressure=1013.0):
        """
        Retrieve the Rayleigh extinction profile (as a function of altitude).

        Parameters
        ----------
        temperature_profile: xarray.DataArray
            a DataArray containing the altitude grid points in [km] and temperatures in [Kelvin].
        surface_pressure: float
            Surface pressure in units [mb] (default value is 1013.0)

        Returns
        -------
        rayleigh_profile: xarray.Dataset
            a DataArray containing the altitude grid points in [km], temperatures in [kelvin] and molecular extinction [km^-1].
        """

        # Use the parameterization of Bodhaine et al. (1999) eq 30 for tau_R at
        # sea level and convert to Rayleigh density coefficient:
        #   k = 0.03370*(p_sfc/1013.25)*tau_R  for k in K/(mb km)
        raylcoefs = 0.03370 * (surface_pressure / 1013.25) * 0.0021520 * (
                    1.0455996 - 341.29061 / self.wavelengths ** 2 - 0.90230850 * self.wavelengths ** 2) / \
                    (1 + 0.0027059889 / self.wavelengths ** 2 - 85.968563 * self.wavelengths ** 2)
        ext_profile = [xr.DataArray(
            dims='Altitude',
            name='{:1.3} micron'.format(wavelength),
            coords={'wavelength': wavelength, 'Altitude': temperature_profile.Altitude},
            data=core.rayleigh_extinct(nzt=temperature_profile.size,
                                       zlevels=temperature_profile.Altitude,
                                       temp=temperature_profile.data,
                                       raysfcpres=surface_pressure,
                                       raylcoef=raylcoef)
        ) for raylcoef, wavelength in zip(raylcoefs, self.wavelengths)]

        # Create a Rayleigh profile xarray.Dataset
        rayleigh_profile = temperature_profile.copy().to_dataset(name='Temperature')
        rayleigh_profile['Extinction'] = xr.concat(ext_profile, dim='wavelength')
        rayleigh_profile.attrs['units'] = temperature_profile.attrs['units']
        rayleigh_profile.attrs['units'].append('Extinction [Km^-1]')
        return rayleigh_profile

    @property
    def wavelengths(self):
        return self._wavelengths

    @property
    def surface_pressure(self):
        return self._surface_pressure
