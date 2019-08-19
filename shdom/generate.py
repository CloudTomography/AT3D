"""
Generation objects generate the atmospheric properties and grid.
"""
from shdom import float_round
import numpy as np
import argparse
import shdom
import copy
from collections import OrderedDict

import scipy.stats as st
import numpy.fft as fft
import random as rand
from scipy.spatial import cKDTree
import numpy.ma as ma

class Generator(object):
    """
    The abstract generator class object which is inherited by specific atmosphere generators. 
       
    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        Arguments required for this generator.
    """
    def __init__(self, args):
        self._args = args
        
    @classmethod
    def update_parser(self, parser): 
        """
        Update the argument parser with parameters relevant to this generator.

        Parameters
        ----------
        parser: argparse.ArgumentParser()
            The main parser to update.

        Returns
        -------
        parser: argparse.ArgumentParser()
            The updated parser.

        Notes
        -----
        Dummy method which will be overwritten.
        """
        return parser

    def get_grid(self):
        """
        Retrieve the atmospheric grid.

        Notes
        -----
        Dummy method which will be overwritten.
        """
        return None

    def get_extinction(self, grid=None):
        """
        Retrieve the atmospheric extinction.

        Notes
        -----
        Dummy method which will be overwritten.
        """
        return None

    def get_albedo(self, grid=None):
        """
        Retrieve the atmospheric albedo.

        Notes
        -----
        Dummy method which will be overwritten.
        """
        self._albedo = albedo

    def get_phase(self, grid=None):
        """
        Retrieve the atmospheric phase.

        Notes
        -----
        Dummy method which will be overwritten.
        """
        return None


    def get_scatterer(self):
        """
        Retrieve a shdom.Scatterer object.

        Notes
        -----
        Dummy method which will be overwritten.
        """
        return None

    @property
    def args(self):
        return self._args


class CloudGenerator(Generator):
    """
    An abstract CloudGenerator class which is inherited by specific generators to generate cloudy voxels.

    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        Arguments required for this generator.
    """
    def __init__(self, args):
        super(CloudGenerator, self).__init__(args)
        self._mie = OrderedDict()

    def add_mie(self, table_path):
        """
        Add a Mie table to the generator.

        Parameters
        ----------
        table_path: string
            A path to read the Mie table from file.

        Notes
        -----
        Multiple Mie tables are added for generation of MultispectralScatterer or MicrophysicalScatterer objects.
        """
        mie = shdom.MiePolydisperse()
        mie.read_table(table_path)
        self._mie[mie.wavelength] = copy.deepcopy(mie)

    def get_lwc(self, grid=None):
        """
        Retrieve the liquid water content.

        Notes
        -----
        Dummy method which will be overwritten.
        """
        return None

    def get_reff(self, grid=None):
        """
        Retrieve the effective radius.

        Notes
        -----
        Dummy method which will be overwritten.
        """
        return None

    def get_veff(self, grid=None):
        """
        Retrieve the effective variance.

        Notes
        -----
        Dummy method which will be overwritten.
        """
        return None

    def get_albedo(self, wavelength, grid=None):
        """
        Retrieve the single scattering albedo at a single wavelength on a grid.

        Parameters
        ----------
        wavelength: float
            Wavelength in microns. A Mie table at this wavelength should be added prior (see add_mie method).
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        albedo: shdom.GridData
            A GridData object containing the single scattering albedo [0,1] on a grid

        Notes
        -----
        The input wavelength is rounded to three decimals.
        """
        if grid is None:
            grid = self.get_grid()
        reff = self.get_reff(grid)
        veff = self.get_veff(grid)
        return self.mie[float_round(wavelength)].get_albedo(reff, veff)

    def get_phase(self, wavelength, grid=None):
        """
        Retrieve the phase function at a single wavelength on a grid.

        Parameters
        ----------
        wavelength: float
            Wavelength in microns. A Mie table at this wavelength should be added prior (see add_mie method).
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        phase: shdom.GridPhase
            A GridPhase object containing the phase function on a grid

        Notes
        -----
        The input wavelength is rounded to three decimals.
        """
        if grid is None:
            grid = self.get_grid()
        reff = self.get_reff(grid)
        veff = self.get_veff(grid)
        return self.mie[float_round(wavelength)].get_phase(reff, veff)

    def get_scatterer(self):
        """
        Return a shdom.Scatterer object
        """
        lwc = self.get_lwc()
        reff = self.get_reff()
        veff = self.get_veff()
        if (lwc is not None) and (reff is not None) and (veff is not None):
            scatterer = shdom.MicrophysicalScatterer(lwc, reff, veff)
            for mie in self.mie.values():
                scatterer.add_mie(mie)

        else:
            if len(self.mie) == 1:
                wavelength = next(iter(self.mie.keys()))
                scatterer = shdom.OpticalScatterer(
                    wavelength,
                    self.get_extinction(wavelength),
                    self.get_albedo(wavelength),
                    self.get_phase(wavelength)
                )

            else:
                scatterer = shdom.MultispectralScatterer()
                for wavelength in self.mie.keys():
                    scatterer.add_scatterer(shdom.OpticalScatterer(
                        wavelength,
                        self.get_extinction(wavelength),
                        self.get_albedo(wavelength),
                        self.get_phase(wavelength)
                    ))

        return scatterer

    @property
    def mie(self):
        return self._mie


class AirGenerator(Generator):
    """
    Generate air optical properties according to rayleigh molecular scattering.

    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        Arguments required for this generator.
    """
    def __init__(self, args):
        super(AirGenerator, self).__init__(args)

    @classmethod
    def update_parser(self, parser):
        """
        Update the argument parser with parameters relevant to this generator.

        Parameters
        ----------
        parser: argparse.ArgumentParser()
            The main parser to update.

        Returns
        -------
        parser: argparse.ArgumentParser()
            The updated parser.
        """
        parser.add_argument('--air_max_alt',
                            default=20.0,
                            type=np.float32,
                            help='(default value: %(default)s) Maximum altitude for the air volume.')

        parser.add_argument('--air_num_points',
                            default=20,
                            type=int,
                            help='(default value: %(default)s) Number of altitude grid points for the air volume.')
        return parser

    def set_temperature_profile(self, temperature_profile):
        """
        Set the temperature profile

        Parameters
        ----------
        temperature_profile: shdom.GridData
            A GridData object specifying the temperature (K) on a 1D grid (altitude)
        """
        self._temperature_profile = temperature_profile

    def get_scatterer(self, wavelength):
        """
        Parameters
        -----------
        wavelength: float or list,
            The wavelength or wavelength list in [microns]

        Returns
        -------
        scatterer: shdom.Scatterer or shdom.MultispectralScatterer
        """
        altitudes = shdom.Grid(z=np.linspace(0.0, self.args.air_max_alt, self.args.air_num_points))
        rayleigh = shdom.Rayleigh(wavelength)
        rayleigh.set_profile(self.temperature_profile.resample(altitudes))
        return rayleigh.get_scatterer()

    @property
    def temperature_profile(self):
        return self._temperature_profile


class AFGLSummerMidLatAir(AirGenerator):
    """
    Temperature profile is taken from AFGL measurements of summer mid-lat.

    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        Arguments required for this generator.
    """
    def __init__(self, args):
        super(AFGLSummerMidLatAir, self).__init__(args)
        temperatures = np.array([292.220, 292.040, 291.860, 291.680, 291.500, 291.320, 291.140, 290.960, 290.780,
                                 290.600, 290.420, 290.240, 290.060, 289.880, 289.700, 289.920, 290.140, 290.360,
                                 290.580, 290.800, 291.020, 291.240, 291.460, 291.680, 291.900])
        temp_grid = shdom.Grid(z=np.linspace(0.0, 20.0, len(temperatures)))
        temperature_profile = shdom.GridData(temp_grid, temperatures)
        self.set_temperature_profile(temperature_profile)


class SingleVoxel(CloudGenerator):
    """
    Define a Medium with a single voxel in center of the grid.
    It is useful for developing and debugging of derivatives and sensitivity analysis.

    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        Arguments required for this generator.
    """
    def __init__(self, args):
        super(SingleVoxel, self).__init__(args)

    @classmethod
    def update_parser(self, parser):
        """
        Update the argument parser with parameters relevant to this generator.

        Parameters
        ----------
        parser: argparse.ArgumentParser()
            The main parser to update.

        Returns
        -------
        parser: argparse.ArgumentParser()
            The updated parser.
        """
        parser.add_argument('--nx',
                            default=10,
                            type=int,
                            help='(default value: %(default)s) Number of grid cell in x (North) direction')
        parser.add_argument('--ny',
                            default=10,
                            type=int,
                            help='(default value: %(default)s) Number of grid cell in y (East) direction')
        parser.add_argument('--nz',
                            default=10,
                            type=int,
                            help='(default value: %(default)s) Number of grid cell in z (Up) direction')
        parser.add_argument('--domain_size',
                            default=1.0,
                            type=float,
                            help='(default value: %(default)s) Cubic domain size [km]')
        parser.add_argument('--extinction',
                            default=10.0,
                            type=np.float32,
                            help='(default value: %(default)s) Extinction of the center voxel [km^-1]')
        parser.add_argument('--lwc',
                            default=None,
                            type=np.float32,
                            help='(default value: %(default)s) Liquid water content of the center voxel [g/m^3]. If specified, extinction argument is ignored.')
        parser.add_argument('--reff',
                            default=10.0,
                            type=np.float32,
                            help='(default value: %(default)s) Effective radius of the center voxel [micron]')
        parser.add_argument('--veff',
                            default=0.1,
                            type=np.float32,
                            help='(default value: %(default)s) Effective variance of the center voxel')
        parser.add_argument('--mie_table_path',
                            help='Path to a precomputed Mie scattering table. \
                                  See notebooks/Make Mie Table.ipynb for more details')
        return parser

    def get_grid(self):
        """
        Retrieve the scatterer grid.

        Returns
        -------
        grid: shdom.Grid
            A Grid object for this scatterer
        """
        bb = shdom.BoundingBox(0.0, 0.0, 0.0, self.args.domain_size, self.args.domain_size, self.args.domain_size)
        return shdom.Grid(bounding_box=bb, nx=self.args.nx, ny=self.args.ny, nz=self.args.nz)

    def get_extinction(self, wavelength=None, grid=None):
        """
        Retrieve the optical extinction at a single wavelength on a grid.

        Parameters
        ----------
        wavelength: float
            Wavelength in microns. A Mie table at this wavelength should be added prior (see add_mie method).
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        extinction: shdom.GridData
            A GridData object containing the optical extinction on a grid

        Notes
        -----
        If the LWC is specified then the extinction is derived using (lwc,reff,veff). If not the extinction needs to be directly specified.
        The input wavelength is rounded to three decimals.
        """
        if grid is None:
            grid = self.get_grid()

        if self.args.lwc is None:
            ext_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
            ext_data[int(grid.nx/2), int(grid.ny/2), int(grid.nz/2)] = self.args.extinction
            extinction = shdom.GridData(grid, ext_data)
        else:
            assert wavelength is not None, 'No wavelength provided'
            lwc = self.get_lwc(grid)
            reff = self.get_reff(grid)
            veff = self.get_veff(grid)
            extinction = self.mie[float_round(wavelength)].get_extinction(lwc, reff, veff)
        return extinction

    def get_lwc(self, grid=None):
        """
        Retrieve the liquid water content on a grid.

        Parameters
        ----------
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        lwc: shdom.GridData
            A GridData object containing liquid water content (g/m^3) on a 3D grid.
        """
        if grid is None:
            grid = self.get_grid()

        lwc = self.args.lwc
        if lwc is not None:
            lwc_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
            lwc_data[int(grid.nx/2), int(grid.ny/2), int(grid.nz/2)] = self.args.lwc
            lwc = shdom.GridData(grid, lwc_data)
        return lwc

    def get_reff(self, grid=None):
        """
        Retrieve the effective radius on a grid.

        Parameters
        ----------
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        reff: shdom.GridData
            A GridData object containing the effective radius [microns] on a grid
        """
        if grid is None:
            grid = self.get_grid()
        reff_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
        reff_data[int(grid.nx/2), int(grid.ny/2), int(grid.nz/2)] = self.args.reff
        return shdom.GridData(grid, reff_data)

    def get_veff(self, grid=None):
        """
        Retrieve the effective variance on a grid.

        Parameters
        ----------
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        veff: shdom.GridData
            A GridData object containing the effective variance on a grid
        """
        if grid is None:
            grid = self.get_grid()
        veff_data = np.zeros(shape=(grid.nx, grid.ny, grid.nz), dtype=np.float32)
        veff_data[int(grid.nx/2), int(grid.ny/2), int(grid.nz/2)] = self.args.veff
        return shdom.GridData(grid, veff_data)


class Homogeneous(CloudGenerator):
    """
    Define a homogeneous Medium.

    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        Arguments required for this generator.
    """
    def __init__(self, args):
        super(Homogeneous, self).__init__(args)

    @classmethod
    def update_parser(self, parser):
        """
        Update the argument parser with parameters relevant to this generator.

        Parameters
        ----------
        parser: argparse.ArgumentParser()
            The main parser to update.

        Returns
        -------
        parser: argparse.ArgumentParser()
            The updated parser.
        """
        parser.add_argument('--nx',
                            default=10,
                            type=int,
                            help='(default value: %(default)s) Number of grid cell in x (North) direction')
        parser.add_argument('--ny',
                            default=10,
                            type=int,
                            help='(default value: %(default)s) Number of grid cell in y (East) direction')
        parser.add_argument('--nz',
                            default=10,
                            type=int,
                            help='(default value: %(default)s) Number of grid cell in z (Up) direction')
        parser.add_argument('--domain_size',
                            default=1.0,
                            type=float,
                            help='(default value: %(default)s) Cubic domain size [km]')
        parser.add_argument('--extinction',
                            default=1.0,
                            type=np.float32,
                            help='(default value: %(default)s) Extinction [km^-1]')
        parser.add_argument('--lwc',
                            default=None,
                            type=np.float32,
                            help='(default value: %(default)s) Liquid water content of the center voxel [g/m^3]. If specified, extinction argument is ignored.')
        parser.add_argument('--reff',
                            default=10.0,
                            type=np.float32,
                            help='(default value: %(default)s) Effective radius [micron]')
        parser.add_argument('--veff',
                            default=0.1,
                            type=np.float32,
                            help='(default value: %(default)s) Effective variance')
        parser.add_argument('--mie_table_path',
                            help='Path to a precomputed polydisperse Mie scattering table. \
                                  See notebooks/Make Mie Table.ipynb for more details')
        return parser

    def get_grid(self):
        """
        Retrieve the scatterer grid.

        Returns
        -------
        grid: shdom.Grid
            A Grid object for this scatterer
        """
        bb = shdom.BoundingBox(0.0, 0.0, 0.0, self.args.domain_size, self.args.domain_size, self.args.domain_size)
        return shdom.Grid(bounding_box=bb, nx=self.args.nx, ny=self.args.ny, nz=self.args.nz)

    def get_extinction(self, wavelength=None, grid=None):
        """
        Retrieve the optical extinction at a single wavelength on a grid.

        Parameters
        ----------
        wavelength: float
            Wavelength in microns. A Mie table at this wavelength should be added prior (see add_mie method).
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        extinction: shdom.GridData
            A GridData object containing the optical extinction on a grid

        Notes
        -----
        If the LWC is specified then the extinction is derived using (lwc,reff,veff). If not the extinction needs to be directly specified.
        The input wavelength is rounded to three decimals.
        """
        if grid is None:
            grid = self.get_grid()

        if self.args.lwc is None:
            ext_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=self.args.extinction, dtype=np.float32)
            extinction = shdom.GridData(grid, ext_data)
        else:
            assert wavelength is not None, 'No wavelength provided'
            lwc = self.get_lwc(grid)
            reff = self.get_reff(grid)
            veff = self.get_veff(grid)
            extinction = self.mie[float_round(wavelength)].get_extinction(lwc, reff, veff)
        return extinction

    def get_lwc(self, grid=None):
        """
        Retrieve the liquid water content.

        Parameters
        ----------
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        lwc: shdom.GridData
            A GridData object containing liquid water content (g/m^3) on a 3D grid.
        """
        if grid is None:
            grid = self.get_grid()

        lwc = self.args.lwc
        if lwc is not None:
            lwc_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=lwc, dtype=np.float32)
            lwc = shdom.GridData(grid, lwc_data)
        return lwc

    def get_reff(self, grid=None):
        """
        Retrieve the effective radius on a grid.

        Parameters
        ----------
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        reff: shdom.GridData
            A GridData object containing the effective radius [microns] on a grid
        """
        if grid is None:
            grid = self.get_grid()
        reff_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=self.args.reff, dtype=np.float32)
        return shdom.GridData(grid, reff_data)

    def get_veff(self, grid=None):
        """
        Retrieve the effective variance on a grid.

        Parameters
        ----------
        grid: shdom.Grid, optional
            A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)

        Returns
        -------
        veff: shdom.GridData
            A GridData object containing the effective variance on a grid
        """
        if grid is None:
            grid = self.get_grid()
        veff_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=self.args.veff, dtype=np.float32)
        return shdom.GridData(grid, veff_data)


class LesFile(CloudGenerator):
    """
    Generate an optical medium from an LES file which contains a table with: lwc, reff, veff (optional).

    Parameters
    ----------
    args: arguments from argparse.ArgumentParser()
        Arguments required for this generator.

    Notes
    -----
    If veff is not specified in the file a single value will be used for the entire domain.
    """
    def __init__(self, args):
        super(LesFile, self).__init__(args)
        self._droplets = shdom.MicrophysicalScatterer()
        self._droplets.load_from_csv(args.path, args.veff)

    @classmethod
    def update_parser(self, parser):
        """
        Update the argument parser with parameters relevant to this generator.

        Parameters
        ----------
        parser: argparse.ArgumentParser()
            The main parser to update.

        Returns
        -------
        parser: argparse.ArgumentParser()
            The updated parser.
        """
        parser.add_argument('--path',
                            help='Path to the LES generated file')
        parser.add_argument('--veff', 
                            default=0.1,
                            type=np.float32, 
                            help='(default value: %(default)s) Effective variance (if not provided as a 3D field)')         
        parser.add_argument('--mie_table_path', 
                            help='Path to a precomputed Mie scattering table. \
                                  See notebooks/Make Mie Table.ipynb for more details')    
        return parser
    
    def get_grid(self):
        """
        Retrieve the scatterer grid.

        Returns
        -------
        grid: shdom.Grid
            A Grid object for this scatterer
        """        
        return self._droplets.grid
    
    
    def get_scatterer(self):
        """
        Retrieve a shdom.MicrophysicalScatterer object.

        Returns
        -------
        droplets: shdom.MicrophysicalScatterer
            A microphysical scatterer object.
        """        
        for mie in self.mie.values():
            self._droplets.add_mie(mie)
        return self._droplets
    
class GaussianFieldGenerator(object):
    """
        An object that stochastically generates 3D fields in a cubic domain drawn from
        truncated normal fields with an isotropic power law spectrum.
        
        Parameters
        ----------
        
        nx, ny, nz: int,
        The number of points in the first, second, and thid dimensions of the field.
        
        beta: float,
        The exponent of the power-law controlling the isotropic spectral slope.
        
        domain_size: list of floats,
        The size of each of dimension of the domain in physical units.
        
        field_min, field_max: float,
        The minimum and maximum values of the 'standard' normal distribution at
        which to truncate.
        
        inner_scale: float,
        The length scale at which fourier components at isotropic frequencies
        above 1.0/inner_scale will be set to zero. Units are the same as domain_size.
        
        outer_scale: float,
        The length scale at which fourier components at isotropic frequencies
        below 1.0/outer_scale are set to a constant value. Units are the same as
        domain_size.
        
        horizontal_direction: array_like of floats with len(3)
        A vector in 3D space specifying the preferred anisotropic direction.
        
        direction_strength: float,
        A parameter describing the strength of the anisotropy.
        
        seed, int (< 2**32)
        A seed for the random number generator to ensure reproducibility of the
        results.
        """
    def __init__(self, nx, ny, nz, beta, domain_size, field_min = None,field_max = None,
                 inner_scale = None, outer_scale = None, seed = None):
        
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.beta = beta
        self.domain_size = domain_size
        
        
        if field_min is not None:
            self.field_min = field_min
        elif field_max is not None:
            self.field_min = -10**-9
        else:
            self.field_min = None

        if field_max is not None:
            self.field_max = field_max
        elif field_min is not None:
            self.field_max = 10**9
        else:
            self.field_max = None
        
        
        if (outer_scale is not None) and (inner_scale is not None):
            assert outer_scale > inner_scale,'power_spectrum outer_scale {} must be \
            larger (greater) than power_spectrum inner_scale {}'.format(outer_scale,inner_scale)
        
        self.outer_scale = outer_scale
        self.inner_scale = inner_scale
        
        if seed is not None:
            np.random.seed(seed)
            self.seed = seed
        else:
            self.seed = None

    def _update_args(self,kwargs):
        """
        Updates the attributes of the generator.
        
        Parameters
        ---------
        kwargs: dictionary
        """
        for name in kwargs:
            setattr(self,name,kwargs[name])
        
    def generate_field(self,return_spectra = False, **kwargs):
        """
        Generates the field with the specified parameters.
        Any parameters of the generator can be updated by supplying them as kwargs.
            
        Parameters
        ---------
        return_spectra: Boolean
            Decides whether the spectra and frequencies of the field are
            also returned by the function.
            
        kwargs: dictionary
            A dictionary of updated parameters for the generator.
            See __init__ for a list of possiblities.
            
        Returns
        -------
        output_field: float, shape=(nx, ny, nz)
            The generated field.
            
        if return_spectra is True this function also returns:
            
        power_spectra: float, shape=(nx, ny, nz)
            The fourier transform of output_field
            
        X,Y,Z: float, shape=(nx, ny, nz)
            The frequencies in the first, second, and third dimensions
            respectively.
            """
        
        self._update_args(kwargs)
        
        # Generate fields, ensuring that they are sufficiently gaussian.
        skew = 3
        kurtosis = 5
        optimal_field = None
        min_skew = 10**9
        min_kurtosis = 10**9
        maxiter = 500
        counter = 0
        
        while ((np.abs(skew)>0.1) or (np.abs(kurtosis-0.0)>0.5)) and counter<maxiter:
            data_inv = self._initialize_field()
            skew = st.skew(data_inv.ravel())
            kurtosis = st.kurtosis(data_inv.ravel())
            if np.abs(skew) < np.abs(min_skew) and np.abs(kurtosis)<np.abs(min_kurtosis):
                min_skew=skew
                min_kurtosis=kurtosis
                optimal_field = data_inv
            counter += 1
        
        # Remap normal distribution to truncated dist by rank ordering.
        if self.field_min is not None and self.field_max is not None:
            assert self.field_min < self.field_max, 'field_min {} must be less than field_max {}'.format(
                    self.field_min,self.field_max)
            trunc_field = st.truncnorm.rvs(self.field_min,self.field_max,loc=0.0,scale=1.0,size=(self.nx,self.ny,self.nz))
            
            sorted_trunc_values = np.sort(trunc_field.ravel())
            sorted_white_args = np.argsort(optimal_field.ravel())
            
            output_field = np.zeros(sorted_trunc_values.shape)
            output_field[sorted_white_args] = sorted_trunc_values
            output_field = np.reshape(output_field,(self.nx,self.ny,self.nz))
        
        else:
            output_field = optimal_field

        if return_spectra:
            fourier = fft.fftshift(fft.fftn(output_field))
    
            x,y,z = fft.fftfreq(self.nx,d=self.domain_size[0]/self.nx),fft.fftfreq(self.ny,d=self.domain_size[1]/self.ny), \
            fft.fftfreq(self.nz,d=self.domain_size[2]/self.nz)
            X,Y,Z = np.meshgrid(fft.fftshift(y),fft.fftshift(x),fft.fftshift(z))
            return output_field,fourier,X,Y,Z
        else:
            return output_field

    def _initialize_field(self):
        """
        Generates a white gaussian field and applies the spectral scaling specified by the generator.
        
        Returns
        -------
        scaled_field: np.array shape=(nx, ny, nz)
        """
        
        if self.nx % 2 == 1:
            nx = self.nx+1
        else:
            nx = self.nx
        if self.ny %2 == 1:
            ny = self.ny+1
        else:
            ny = self.ny
        if self.nz % 2 ==1:
            nz = self.nz+1
        else:
            nz=self.nz
        
        white_field = np.random.normal(loc=0.0,scale=1.0,size=(nx,ny,nz))
        white_spectra = fft.fftshift(fft.fftn(white_field))
        x,y,z = fft.fftfreq(nx,d=self.domain_size[0]/nx),fft.fftfreq(ny,d=self.domain_size[1]/ny), \
        fft.fftfreq(nz,d=self.domain_size[2]/nz)
            
            
        X,Y,Z = np.meshgrid(fft.fftshift(x),fft.fftshift(y),fft.fftshift(z),indexing='ij')
        isotropic_wavenumber = np.sqrt(X**2+Y**2+Z**2)
        #np.where(isotropic_wavenumber==0.0)
        
        #if self.ny % 2 == 1 and self.nx % 2 == 0:
        #    isotropic_wavenumber[self.nx//2,]
            
        isotropic_wavenumber[nx//2,ny//2,nz//2] = 10**-6
        spectral_scaling = np.abs(isotropic_wavenumber)**(self.beta/2.0)
        spectral_scaling[nx//2,ny//2,nz//2] = 0.0         
        scaled_spectra = white_spectra*spectral_scaling

        if self.inner_scale is not None:
            scaled_spectra[np.where(isotropic_wavenumber>1.0/self.inner_scale)] = 0.0
                                
        if self.outer_scale is not None:
        #Power above scale-break is an azimuthally averaged power around the scale break.
            value = np.mean(np.abs(scaled_spectra[np.where(np.abs(isotropic_wavenumber-1.0/self.outer_scale)<0.3)]))
            scaled_spectra[np.where(isotropic_wavenumber<1.0/self.outer_scale)] = value
                               
        scaled_field = np.real(fft.ifftn(fft.fftshift(scaled_spectra)))[:self.nx,:self.ny,:self.nz]
        scaled_field = (scaled_field - scaled_field.mean())/scaled_field.std()
        
        return scaled_field
    
class StochasticCloudGenerator(object):
    """
    Generates clouds with either idealised geometry or a cloud mask from thresholding a stochastic field.
    Cloud microphysical variables are generated based on some heuristics of the statistics 
    of clouds. Deterministic scalings can also be applied which allow pseudo realistic
    as well as anti-realistic cloud structures to be generated.
    
    LWC is drawn from a log-normal distribution while Reff and Veff are drawn from truncated 
    gaussian distributions. Each of these fields have an isotropic power spectrum described by a power law
    S(k) = a k^beta.
    This spectrum can be set to zero at frequencies above a specified scale or set to a constant value
    at frequencies below another scale.
    
    A power-law correlation between LWC and Reff of the form Reff = a*LWC^b may be specified.
    This value is set by default to 1/3 following a constant droplet number concentration and Veff.
    
    Polynomial dependences of LWC, Reff and Veff may be specified on both the height above
    cloud base and the minimum distance to a cloud free voxel.
        
    Parameters
    ----------
    
    nx, ny, nz: int
        The number of voxels in first, second and third (vertical) directions.
        Note the resolution in the vertical is set by aspect_ratio*domain_size/nz.
        
    domain_size: float,
        The horizontal length of the domain in physical units.
        
    aspect_ratio: float,
        The ratio of the vertical dimension of the domain to the horizontal scale.
        This parameter also sets the aspect ratio of the idealised geometries.
        
    lwc_mean, reff_mean, veff_mean: float,
        The mean liquid water content [g/m^3], mean droplet effective radius [micrometer]
        and mean droplet effective variance [micrometer^2] in the cloud.
        
    lwc_snr, reff_snr, veff_snr: float,
        The strength of the random fluctuations of liquid water content,
        droplet effective radius, and droplet effective variance relative to the mean.
        lwc_mean.
        
    beta: float,
        The slope parameter of the power spectrum. More negative values means less
        variance at small scales.
        
    cloud_geometry: string,
        Cloud geometry options include ['Sphere','Cuboid','Cylinder', 'Blob','Stochastic']
        where Blob is a filtered Stochastic cloud to make it appear more cohesive and
        focused in the centre of the domain. Idealised (Sphere, Cuboid & Cylinder) geometries
        are centered in the domain centre.
        
    cloud_fraction: float,
        The volumetric fraction of the domain occupied by cloudy voxels. [0,1]
        
    lwc_vertical_dependence, reff_vertical_dependence, veff_vertical_dependence: list of floats,
        Each of these lists contains the coefficients (starting from the linear coefficient)
        of the polynomial that describes the dependence of the microphysical variable on
        altitude above the cloud base height. Units are in lwc_mean/domain_size, lwc_mean/domain_size**2 etc.
        
    lwc_interior_dependence, reff_interior_dependence, veff_interior_dependence: list of floats,
        Each of these lists contains the coefficients (starting from the linear coefficient)
        of the polynomial that describes the dependence of the microphysical variable on the
        minimum distance to a cloud free voxel. Units are in lwc_mean/domain_size, lwc_mean/domain_size**2 etc.
        
    preset: String,
        A string which chooses from among some useful preset combinations of the vertical and 
        interior dependences. This overwrites SOME of the specified dependences. See _process_presets(self) 
        for more details and for inspiration in manually setting the dependences.
        
    spectrum_outer_scale: float,
        Sets the scale (wavelength) in physical units (same as domain_size) above which the power spectrum
        is set to a constant value (the value obtained at the specified scale).
    
    spectrum_inner_scale: float,
        Sets the scale (wavelength) in physical units (same as domain_size) below which the power spectrum
        is set to zero.
        
    cloud_base_height: float,
        A parameter of the Stochastic and Blob cloud geometry generation which specifies the minimum altitude
        of the cloud above the ground in the same units as domain_size.
        
    cloud_mask_beta: float,
        A parameter of the Stochastic and Blob cloud geometry generation which specifies
        how the power-law slope of the field used to generate the cloud mask.
        More negative values make the cloud mask smoother.
        
    cloud_mask_vertical_correlation: float,
        A parameter of the Stochastic and Blob cloud geometry generation which controls 
        the vertical overlap of the cloud mask.
        
    blob_parameter: float,
        A parameter of the Blob cloud geometry generation which controls how focused the blob 
        is in the centre of the domain. Values above 1.0 are the same as the Stochastic 
        geometry generator.
        
    reff_min, reff_max, veff_min, veff_max: float,
        These are the min and max or reff and veff, respectively based on the default
        specifications of a mie table. These variables are truncated to these limits
        either by hard thresholding (reff) or rank mapping to a truncated gaussian (veff).
        This prevents errors when calculating the optical properties later but can be 
        relaxed if larger tables are generated.
     
    Notes
    -----
    For simple use cases simply vary the mean and snr of each variable.
 
    """
    def __init__(self, **kwargs):
        
        self.default_args = {
                  'nx': 50,
                  'ny': 50,
                  'nz': 50,
                  'domain_size': 1.0,
                  'aspect_ratio': 1.0,
                  'lwc_mean': 0.1,
                  'lwc_snr': 0.1,
                  'reff_mean': 10.0,
                  'reff_snr': 5.0,
                  'veff_mean': 0.1,
                  'veff_snr': 0.01,
                  'beta': -5.0/3.0,
                  'cloud_geometry': 'Blob',
                  'preset': 'None',
                  'spectrum_outer_scale': None,
                  'spectrum_inner_scale': None,
                  'lwc_vertical_dependence': [0.0],
                  'lwc_interior_dependence': [0.0],
                  'reff_vertical_dependence': [0.0],
                  'reff_interior_dependence': [0.0],
                  'veff_vertical_dependence': [0.0],
                  'veff_interior_dependence': [0.0],
                  'lwc_reff_correlation': 0.0,
                  'lwc_reff_correlation_slope': 1.0/3.0,
                  'cloud_base_height': 0.1,
                  'reff_min': 4.0,
                  'reff_max': 25.0,
                  'veff_min': 0.01,
                  'veff_max': 0.2,
                  'blob_parameter': 0.2,
                  'cloud_mask_beta': -3.0,
                  'cloud_fraction': 0.4,
                  'cloud_mask_vertical_correlation': 0.8
                  }
        
        for name in self.default_args:
            setattr(self,name,self.default_args[name])
        
        #Hack so unused keyword arguments are ignored
        #instead of filtering all args provided to the generator (such as output_dir)
        unused_kwargs = {}
        for name in kwargs:
            if name not in self.default_args:
                unused_kwargs[name] = kwargs[name]
            else:
                setattr(self,name,kwargs[name])
        
        
        self.list_of_presets = ['Hollow','Inverted','Pseudo_real','None']
        assert self.preset in self.list_of_presets, 'Preset: {}  is not supported. Please select one of {}'.format(
                self.preset,self.list_of_presets)
        self._process_presets()
            
        self.list_of_cloud_geometries = ['Sphere','Cuboid','Cylinder', 'Blob','Stochastic']        
        assert self.cloud_geometry in self.list_of_cloud_geometries, 'Geometry: {}  is not supported. Please select one of {}'.format(
                self.cloud_geometry,self.list_of_cloud_geometries)
        
     
        #calculate lower and upper bounds of truncated normal distributions to avoid 
        #exceeding mie_table bounds.
        self.reff_lower = (self.reff_min - self.reff_mean) / (self.reff_snr*self.reff_mean)
        self.reff_upper = (self.reff_max - self.reff_mean) / (self.reff_snr*self.reff_mean)
        self.veff_lower = (self.veff_min - self.veff_mean) / (self.veff_snr*self.veff_mean)
        self.veff_upper = (self.veff_max - self.veff_mean) / (self.veff_snr*self.veff_mean)

        #Initialize stochastic generator. Specific arguments are updated later.
        self.gaussian_field_generator = GaussianFieldGenerator(self.nx,self.ny,self.nz,self.beta,
                                                               [self.domain_size,self.domain_size,
                                                                self.aspect_ratio*self.domain_size])
          
        self.vertical_resolution = self.aspect_ratio*self.domain_size/self.nz         
    def _update_args(self,**kwargs):
        """Updates parameters of the cloud generator."""
        for name in kwargs:
            assert name in self.default_args, 'Invalid Keyword Argument'
            setattr(self, name, kwargs[name])

    def _process_presets(self):
        """Specifies the vertical and interior dependencies.
            This overwrites the user-specified dependencies."""
        if self.preset == 'Hollow':
            self.lwc_vertical_dependence = [0.0]
            self.lwc_interior_dependence = [100.0,-400.0]
        elif self.preset == 'Inverted':
            self.lwc_vertical_dependence = [-2.0]
            self.lwc_interior_dependence = [-5.0]      
        elif self.preset == 'Pseudo_real':
            self.lwc_vertical_dependence = [2.0]
            self.lwc_interior_dependence = [100.0,-200.0]
            self.reff_vertical_dependence = [2.0]
            self.reff_interior_dependence = [100.0,-1200.0]
        else:
            pass

    def generate_cloud_field(self,**kwargs):
        """
        TODO
        """        
        seed = rand.randrange(2**32)
        np.random.seed(seed)
        self._update_args(**kwargs)
        
        self.cloud_mask = self._generate_cloud_mask()     
        # Generate random fields for microphysics.
        #Note that the ordering of the calls of .generate_field() is important.
        LWC_field = ma.masked_array(data = self.gaussian_field_generator.generate_field(beta=self.beta,
                                                          inner_scale=self.spectrum_inner_scale,
                                                          outer_scale = self.spectrum_outer_scale),
                                    mask=self.cloud_mask)

        Reff_field = ma.masked_array(data = self.gaussian_field_generator.generate_field(
                                             field_min=self.reff_lower,
                                             field_max=self.reff_upper),
                                    mask = self.cloud_mask)

        Veff_field = ma.masked_array(data = self.gaussian_field_generator.generate_field(
                                             field_min=self.veff_lower,
                                             field_max=self.veff_upper),
                                    mask = self.cloud_mask)
        
        # Induce correlation between Reff and LWC
        LWC_Reff_mix = LWC_field*self.lwc_reff_correlation + \
                            np.sqrt(1.0-self.lwc_reff_correlation**2)*Reff_field
        LWC_Reff_mix = LWC_Reff_mix - LWC_Reff_mix.min()
        Reff_correlated = LWC_Reff_mix**self.lwc_reff_correlation_slope/ \
                        np.mean(LWC_Reff_mix**self.lwc_reff_correlation_slope)
        Reff_correlated = (Reff_correlated - 1.0)*self.reff_mean

        LWC_deterministic =  self._get_deterministic_scaling(self.lwc_vertical_dependence,
                                                    self.lwc_interior_dependence,
                                                    self.lwc_mean)
        Reff_deterministic =  self._get_deterministic_scaling(self.reff_vertical_dependence,
                                                    self.reff_interior_dependence,
                                                    self.reff_mean)        
        Veff_deterministic =  self._get_deterministic_scaling(self.veff_vertical_dependence,
                                                    self.veff_interior_dependence,
                                                    self.veff_mean)  
        #Gaussian statistics
        Veff_total = Veff_field*self.veff_snr*self.veff_mean + Veff_deterministic
        Reff_total = ma.masked_array(data=Reff_correlated*self.reff_snr*self.reff_mean + Reff_deterministic,
                                     mask=self.cloud_mask)
        
        #log-normal statistics
        LWC_field = (LWC_field - LWC_field.mean())/LWC_field.std()
        #hack to avoid a divide by zero warning in log calculation for masked values.
        calculation_deterministic = ma.filled(LWC_deterministic,1.0)
        mu = np.log(calculation_deterministic**2/np.sqrt((self.lwc_snr*self.lwc_mean)**2 + calculation_deterministic**2))
        sigma = np.sqrt(np.log(1 + (self.lwc_snr*self.lwc_mean)**2/calculation_deterministic**2))

        LWC_total = np.exp(LWC_field*sigma + mu)
        
        #Exclude values where deterministic is negative.
        #Deterministic can cause large negative values. These values end up as extremely high positive values after
        #log-normal scaling.
        LWC_total[ma.where(LWC_deterministic < 0.0 )] = ma.masked
      
        # Set any negative values caused by deterministic scaling
        # masked.
        LWC_total[ma.where((LWC_total < 0.0))] = ma.masked
        Reff_total[ma.where((LWC_total < 0.0))] = ma.masked
        Veff_total[ma.where((LWC_total < 0.0))] = ma.masked    

        Reff_total[ma.where(Reff_total<=self.reff_min)] = self.reff_min+10**-6
        Reff_total[ma.where(Reff_total>=self.reff_max)] = self.reff_max-10**-6
        Reff_total.mask = LWC_total.mask
        
        #Remap Veff onto truncated gaussians by rank ordering
        sorted_values_veff = ma.sort(Veff_field.ravel())*self.veff_snr*self.veff_mean+self.veff_mean
        sorted_args_veff = ma.argsort(Veff_total.ravel())
        
        Veff = ma.masked_array(data=np.zeros(sorted_values_veff.shape))
        Veff[sorted_args_veff] = sorted_values_veff
        Veff = ma.reshape(Veff,Veff_field.shape)   
        Veff.mask = LWC_total.mask
        
        
        return LWC_total,Reff_total,Veff,seed
    
    def _generate_cloud_mask(self):
        """Chooses the appropriate cloud mask generator."""
        if self.cloud_geometry == 'Sphere':
            cloud_mask = self._generate_spherical_cloud_mask()
        elif self.cloud_geometry == 'Cuboid':
            cloud_mask = self._generate_cuboidal_cloud_mask()
        elif self.cloud_geometry == 'Cylinder':
            cloud_mask = self._generate_cylindrical_cloud_mask()
        elif self.cloud_geometry == 'Stochastic':
            cloud_mask = self._generate_stochastic_cloud_mask(blob_cloud=False)
        elif self.cloud_geometry == 'Blob':
            cloud_mask = self._generate_stochastic_cloud_mask(blob_cloud=True)
            
        return cloud_mask
   
    def _generate_spherical_cloud_mask(self):
        """Generates a spherical cloud mask"""
        radius = (3.0*self.cloud_fraction/(4.0*np.pi))**(1.0/3.0)
        cloud_mask = np.ones((self.nx,self.ny,self.nz))
        
        X,Y,Z = np.meshgrid(np.linspace(-1*self.domain_size/2,self.domain_size/2,self.nx),
                            np.linspace(-1*self.domain_size/2,self.domain_size/2,self.ny),
                            np.linspace(-1*self.aspect_ratio*self.domain_size/2,self.aspect_ratio*self.domain_size/2,self.nz),
                            indexing='ij')
        R = np.sqrt(X**2+Y**2+Z**2)
        cloud_mask[np.where(R<radius)] = 0
        
        # Ensure there is no cloud at domain edges
        # Especially horizontal edges as these are used in plane-parallel calculations and
        # strongly affect images.
        cloud_mask[0,:,:] = cloud_mask[-1,:,:] = cloud_mask[:,0,:] = \
        cloud_mask[:,-1,:] = cloud_mask[:,:,-1] = cloud_mask[:,:,0] = 1
        
        self.surface_gap = int((self.aspect_ratio**2*radius/2)/self.vertical_resolution)        
        
        return cloud_mask
    
    def _generate_cuboidal_cloud_mask(self):
        """Generates a cuboidal cloud mask."""
        dimension = (self.cloud_fraction/self.aspect_ratio)**(1.0/3.0)
        cloud_mask = np.ones((self.nx,self.ny,self.nz))
        
        X,Y,Z = np.meshgrid(np.linspace(-1*self.domain_size/2,self.domain_size/2,self.nx),
                            np.linspace(-1*self.domain_size/2,self.domain_size/2,self.ny),
                            np.linspace(-1*self.aspect_ratio*self.domain_size/2,self.aspect_ratio*self.domain_size/2,self.nz),
                            indexing='ij')
        
        cloud_mask[np.where((np.abs(X)<dimension/2)&(np.abs(Y)<dimension/2)& \
                (np.abs(Z)<self.aspect_ratio**2*dimension/2))] = 0.0

        # Ensure there is no cloud at domain edges
        # Especially horizontal edges as these are used in plane-parallel calculations and
        # strongly affect images.
        cloud_mask[0,:,:] = cloud_mask[-1,:,:] = cloud_mask[:,0,:] = \
        cloud_mask[:,-1,:] = cloud_mask[:,:,-1] = cloud_mask[:,:,0] = 1
        
        self.surface_gap = int((self.aspect_ratio**2*dimension/2)/self.vertical_resolution)

        return cloud_mask
    
    def _generate_cylindrical_cloud_mask(self):
        """Generates a cylindrical cloud mask."""
        dimension = (self.cloud_fraction/(self.aspect_ratio*np.pi))**(1.0/3.0)
        cloud_mask = np.ones((self.nx,self.ny,self.nz))
        
        X,Y,Z = np.meshgrid(np.linspace(-1*self.domain_size/2,self.domain_size/2,self.nx),
                            np.linspace(-1*self.domain_size/2,self.domain_size/2,self.ny),
                            np.linspace(-1*self.aspect_ratio*self.domain_size/2,self.aspect_ratio*self.domain_size/2,self.nz),
                            indexing='ij')
        R = np.sqrt(X**2+Y**2)
        cloud_mask[np.where((R<dimension)& \
                (np.abs(Z)<self.aspect_ratio**2*dimension/2))] = 0.0
        
        # Ensure there is no cloud at domain edges
        # Especially horizontal edges as these are used in plane-parallel calculations and
        # strongly affect images.
        cloud_mask[0,:,:] = cloud_mask[-1,:,:] = cloud_mask[:,0,:] = \
        cloud_mask[:,-1,:] = cloud_mask[:,:,-1] = cloud_mask[:,:,0] = 1
        
        self.surface_gap = int((self.aspect_ratio**2*dimension/2)/self.vertical_resolution)        
    
        return cloud_mask        
    
    def _generate_stochastic_cloud_mask(self,blob_cloud=True):
        """
        Generates a stochastic cloud mask by thresholding a vertically correlated 
        3D field with lognormal statistics and a power-law scaling of power spectrum 
        to a specified cloud fraction.
        If blob_cloud is True then this field is filtered to preferentially
        create a blob shape in the centre of the domain.
        
        Parameters
        ----------
        blob_cloud: Boolean,
            Decides whether filtering is applied.
            
        Notes
        -----
        The character of the stochastic cloud mask is affected by many parameters
        passed to the Stochastic Blob Generator.
        """     
        
        
        cloud_mask_field_temp = self.gaussian_field_generator.generate_field(beta=self.cloud_mask_beta)
        
        self.surface_gap = int(self.cloud_base_height/self.vertical_resolution)
        
        #Add vertical correlation
        cloud_mask_field = cloud_mask_field_temp.copy()
        for i in range(cloud_mask_field.shape[-1]-1):
            cloud_mask_field[:,:,i+1] = self.cloud_mask_vertical_correlation*cloud_mask_field[:,:,i] \
            + np.sqrt(1.0 - self.cloud_mask_vertical_correlation)*cloud_mask_field_temp[:,:,i+1]
        
        cloud_mask_field = np.exp(cloud_mask_field)
        
        if blob_cloud:
            #Apply blob filter (fourth order butterworth)
            X,Y,Z = np.meshgrid(np.linspace(-1*self.domain_size/2,self.domain_size/2,self.nx),
                                np.linspace(-1*self.domain_size/2,self.domain_size/2,self.ny),
                                np.linspace(-1*self.aspect_ratio*self.domain_size/2,self.aspect_ratio*self.domain_size/2,self.nz),
                                indexing='ij')
            R = np.sqrt(X**2+Y**2)/(self.blob_parameter*self.domain_size)
            butter_worth_filter = 1.0/np.sqrt(1 + R**8)
            cloud_mask_field = cloud_mask_field*butter_worth_filter      
        
        threshold = np.percentile(cloud_mask_field,(1.0-self.cloud_fraction)*100)
        cloud_mask_final = np.ones(cloud_mask_field.shape)
        cloud_mask_final[np.where(cloud_mask_field > threshold)] = 0
        cloud_mask_final[:,:,:self.surface_gap] = 1
        
        # Ensure there is no cloud at domain edges
        # Especially horizontal edges as these are used in plane-parallel calculations and
        # strongly affect images.
        cloud_mask_final[0,:,:] = cloud_mask_final[-1,:,:] = cloud_mask_final[:,0,:] = \
        cloud_mask_final[:,-1,:] = cloud_mask_final[:,:,-1] = cloud_mask_final[:,:,0] = 1
        
        return cloud_mask_final

    def _get_deterministic_scaling(self,vertical_dependence,interior_dependence, field_mean):
        """
        Calculates the polynomial dependence of a microphysical variable on 
        altitude above cloud_base_height and minimum distance from a cloud free
        voxel.
        
        Parameters
        ----------
        vertical_dependence: list of floats,
            Coefficients of the polynomial controlling the dependence on height
            above cloud_base_height
        interior_dependence: list of floats,
            Coefficients of the polynomial controlling the dependence on minimum
            distance to cloud free voxels.
        field_mean: float,
            The mean of the microphysical variable over the entire domain.
        """
        #=============== Vertical Scaling =========================
        z = np.linspace(0,1.0,self.nz)
        X,Y,Z = np.meshgrid(np.ones(self.nx),np.ones(self.ny),\
                            np.linspace(0,1.0*self.aspect_ratio,self.nz))
                              
        mid_point = (z[self.surface_gap] + z[-1])/2.0
        vertical_polynomial = []
        for i,coefficient in enumerate(vertical_dependence):
            vertical_polynomial.append(coefficient*(Z-mid_point)**(i+1))
        vertical_polynomial = np.array(vertical_polynomial).sum(axis=0)
        
        data = (vertical_polynomial+1.0)*field_mean
        vertical_field = ma.masked_array(data=data,mask=self.cloud_mask)
        
        #============= Interior Scaling =========================
        
        # Calculate minimum distance from cloud to cloud_free using KDTree
        self.distance_matrix = np.zeros(self.cloud_mask.shape)

        a,b,c = np.where(self.cloud_mask==1)
        a,b,c = (a/self.nx)*self.domain_size, b/self.ny*self.domain_size, c/self.nz*self.domain_size*self.aspect_ratio
        no_cloud_positions = np.stack((a,b,c),axis=1)
        tree = cKDTree(no_cloud_positions)
        
        q,w,e = np.where(self.cloud_mask==0)
        q2,w2,e2 = (q/self.nx)*self.domain_size, w/self.ny*self.domain_size, e/self.nz*self.domain_size*self.aspect_ratio
        cloud_positions = np.stack((q2,w2,e2),axis=1)
        d,i = tree.query(cloud_positions)
        
        self.distance_matrix[q,w,e] = d

        interior_polynomial = []
        for i,coefficient in enumerate(interior_dependence):
            interior_polynomial.append(coefficient*(self.distance_matrix)**(i+1))
        interior_polynomial = np.array(interior_polynomial).sum(axis=0)
    
        deterministic_field = ma.masked_array(data=(interior_polynomial*vertical_field.mean(axis=(0,1))),
                                              mask=self.cloud_mask)
        deterministic_field = deterministic_field - (deterministic_field.mean(axis=(0,1))-vertical_field.mean(axis=(0,1)))
        
        return deterministic_field

class StochasticCloud(CloudGenerator):
    """
    Generate an optical medium using stochastically and deterministically generated cloud microphysical properties.
    This generator is designed to be flexible enough to generate clouds with pseudo-real or anti-real structure.
    While focused on stochastic 'Blob' clouds, this generator can also produce clouds with idealised geometry
    such as Spheres, Cuboids & Cylinders with generated microphysics.
    
    ----------
    Parameters
    
    args: arguments from argparse.ArgumentParser()
    Arguments required for this generator.
    
    Notes
    -----
    The stochastic algorithm for generating the cloud is based on the following heuristics:
    
    Following from the scaling of variance in a turbulent field, cloud properties (lwc, reff, veff)
    are scale invariant with a power spectrum that follows a S(k) = a k^(b) where b is a constant and k
    is the isotropic wavenumber. The default is b = -5/3 from observations of clouds.
    
    lwc has log-normal statistics while reff and veff have normal statistics.
    The truncation of reff and veff is to ensure that values outside the bounds of the default mie-table generation
    are not generated. These bounds are reff in [4,25] and veff in [0.01,0.2].
    
    Power law correlations can be induced between lwc and reff so that reff = a lwc^(b).
    A default of b=1/3 is chosen to mimic the physical case of a constant droplet number concentration and veff assumption.
    
    Polynomial dependencies of lwc, reff or veff on height above cloud base or minimum distance from cloud edge can also be applied allowing
    the creation of more realistic or deliberately unrealistic cloud structures (internal holes).
    """
    
    def __init__(self, args):
        super(StochasticCloud, self).__init__(args)
        
        self.stochastic_cloud_generator = StochasticCloudGenerator()
        filtered_args = self._filter_args()        
        lwc_data,reff_data,veff_data,seed = self.stochastic_cloud_generator.generate_cloud_field(**filtered_args)
        
        self.lwc_data = lwc_data.astype(np.float32)
        self.reff_data = reff_data.astype(np.float32)
        self.veff_data = veff_data.astype(np.float32)
        self.seed = seed
        
    def _filter_args(self):
        filtered_args = {}
        default_args = self.stochastic_cloud_generator.default_args
        
        for name in vars(self.args):
            if name in default_args:
                filtered_args[name] = vars(self.args)[name]
        
        if filtered_args['spectrum_inner_scale'] == 0.0:
            filtered_args['spectrum_inner_scale'] = None
        if filtered_args['spectrum_outer_scale'] == 0.0:
            filtered_args['spectrum_outer_scale'] = None
        return filtered_args   
        

    @classmethod
    def update_parser(self,parser):
        parser.add_argument('--nx',
                            default=50,
                            type=int,
                            help='(default value: %(default)s) Number of grid cell in x (North) direction')
        parser.add_argument('--ny',
                            default=50,
                            type=int,
                            help='(default value: %(default)s) Number of grid cell in y (East) direction')
        parser.add_argument('--nz',
                            default=50,
                            type=int,
                            help='(default value: %(default)s) Number of grid cell in z (Up) direction')
        parser.add_argument('--domain_size',
                            default=1.0,
                            type=float,
                            help='(default value: %(default)s) Horizontal dimension of domain [km]')
        parser.add_argument('--aspect_ratio',
                            default=1.0,
                            type=float,
                            help='(default value: %(default)s) The vertical:horizontal aspect ratio of the domain \
                            (and idealised cloud geometries).')
        parser.add_argument('--cloud_geometry',
                            default='Blob',
                            type=str,
                            help='(default value: %(default)s) The choice of cloud geometries. Choose from \
                            Blob, Stochastic, Sphere, Cylinder & Cuboid.')
        parser.add_argument('--mie_table_path',
                            help='Path to a precomputed polydisperse Mie scattering table. \
                            See notebooks/Make Mie Table.ipynb for more details')
        parser.add_argument('--lwc_mean',
                            default=0.1,
                            type=np.float32,
                            help='(default value: %(default)s) Mean Liquid water content [g/m^3]')
        parser.add_argument('--lwc_snr',
                            default=0.1,
                            type=np.float32,
                            help='(default value: %(default)s) Magnitude of random variations of lwc relative to lwc_mean.')
        parser.add_argument('--reff_mean',
                            default=10.0,
                            type=np.float32,
                            help='(default value: %(default)s) Mean Effective radius [micron]')
        parser.add_argument('--reff_snr',
                            default=1.0,
                            type=np.float32,
                            help='(default value: %(default)s) Magnitude of random variations of reff relative to reff_mean.')
        parser.add_argument('--veff_mean',
                            default=0.1,
                            type=np.float32,
                            help='(default value: %(default)s) Mean Effective variance')
        parser.add_argument('--veff_snr',
                            default=0.1,
                            type=np.float32,
                            help='(default value: %(default)s) Magnitude of random variations of veff relative to veff_mean.')
        parser.add_argument('--beta',
                            default=-5.0/3,
                            type=np.float32,
                            help='(default value: %(default)s) Isotropic Scaling Exponent of Power Spectrum of \
                            microphysical properties (lwc, reff, veff)')
        parser.add_argument('--lwc_reff_correlation',
                            default=0.0,
                            type=np.float32,
                            help='(default value: %(default)s) A coefficient controlling the (power law) correlation \
                            between lwc and reff. Allowable Range: [0,1]')
        parser.add_argument('--lwc_reff_correlation_slope',
                            default=1.0/3.0,
                            type=np.float32,
                            help='(default value: %(default)s) A coefficient controlling the slope of the power law \
                            correlation between lwc and reff.')
        parser.add_argument('--cloud_mask_beta',
                            default=-3,
                            type=np.float32,
                            help='(default value: %(default)s) Isotropic Scaling Exponent of Power Spectrum used to \
                            generate the volumetric cloud mask.')
        parser.add_argument('--cloud_mask_vertical_correlation',
                            default=0.8,
                            type=np.float32,
                            help='(default value: %(default)s) Correlation coefficient used to \
                            control vertical overlap of the cloud mask. Allowable Range [0,1].')
        parser.add_argument('--cloud_fraction',
                            default=0.3,
                            type=np.float32,
                            help='(default value: %(default)s) Volumetric cloud fraction over domain in the cloud mask generation. \
                            Allowable Range [0,1].')
        parser.add_argument('--cloud_base_height',
                            default=0.1,
                            type=np.float32,
                            help='(default value: %(default)s) Cloud Base Height [km]')
        parser.add_argument('--lwc_vertical_dependence',
                            default=[0.0],
                            nargs='+',
                            type=np.float32,
                            help='(default value: %(default)s) Coefficients of the polynomial describing the functional dependence of layer-mean lwc  \
                            on altitude above cloud base height (units of lwc_mean/domain_size, lwc_mean/domain_size^2 etc). The first entry \
                            is the linear coefficient.')
        parser.add_argument('--reff_vertical_dependence',
                            default=[0.0],
                            nargs='+',
                            type=np.float32,
                            help='(default value: %(default)s) Coefficients of the polynomial describing the functional dependence of layer-mean reff  \
                            on altitude above cloud base height (units of reff_mean/domain_size, reff_mean/domain_size^2 etc). The first entry \
                            is the linear coefficient.')
        parser.add_argument('--lwc_interior_dependence',
                            default=[0.0],
                            nargs='+',
                            type=np.float32,
                            help='(default value: %(default)s) Coefficients of the polynomial describing the functional dependence of lwc \
                            on distance from nearest clear sky voxel (units of lwc_mean/domain_size, lwc_mean/domain_size^2 etc). The first entry \
                            is the linear coefficient.')
        parser.add_argument('--reff_interior_dependence',
                            default=[0.0],
                            nargs='+',
                            type=np.float32,
                            help='(default value: %(default)s) Coefficients of the polynomial describing the functional dependence of reff \
                            on distance from nearest clear sky voxel (units of reff_mean/domain_size, reff_mean/domain_size^2 etc). The first entry \
                            is the linear coefficient.')
        parser.add_argument('--spectrum_inner_scale',
                            default=0.0,
                            type=np.float32,
                            help='(default value% %(default)s) The minimum (isotropic) scale over which variability is created [km] \
                                 when cloud microphysics are generated [km]. This scale break is not applied if the scale is set zero. \
                                 This scale must be smaller than the power_spectrum_outer_scale')
        parser.add_argument('--spectrum_outer_scale',
                            default=0.0,
                            type=np.float32,
                            help='(default value% %(default)s) The (isotropic) scale beyond which larger scales have white spectra [km]. \
                                 This scale break is not applied if the scale is set zero. This scale must be \
                                 larger than the power_spectrum_inner_scale')
        parser.add_argument('--reff_min',
                            default=4.0,
                            type=np.float32,
                            help='(default value: %default)s) The minimum of droplet effective radius. Values below \
                            this are set to this value.')
        parser.add_argument('--reff_max',
                            default=25.0,
                            type=np.float32,
                            help='(default value: %default)s) The maximum of droplet effective radius. Values above this are \
                            set to this value.')
        parser.add_argument('--veff_min',
                            default=0.01,
                            type=np.float32,
                            help='(default value: %default)s) The minimum of droplet effective variance. Values below \
                            this are set to this value.')
        parser.add_argument('--veff_max',
                            default=0.2,
                            type=np.float32,
                            help='(default value: %default)s) The maximum of droplet effective variance. Values above this are \
                            set to this value.')
        parser.add_argument('--blob_parameter',
                            default=0.1,
                            help='(default value: %default)s) A parameter of the Blob cloud geometry choice. \
                            Larger values broaden the size of the blob in the centre while smaller values will \
                            localize the cloud more strongly.')
        return parser
    
    def get_grid(self):
        bb = shdom.BoundingBox(0.0, 0.0, 0.0, self.args.domain_size, self.args.domain_size, self.args.domain_size*self.args.aspect_ratio)
        return shdom.Grid(bounding_box=bb, nx=self.args.nx, ny=self.args.ny, nz=self.args.nz)
    
    
    def get_lwc(self, grid=None):
        """
        Retrieve the liquid water content.
            
        Parameters
        ----------
        grid: shdom.Grid, optional
        A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)
            
        Returns
        -------
        lwc: shdom.GridData
        A GridData object containing liquid water content (g/m^3) on a 3D grid.
        """
        if grid is None:
            grid = self.get_grid()
        
        
        
        lwc = shdom.GridData(grid,ma.filled(self.lwc_data,0.0))
        return lwc
    
    def get_reff(self, grid=None):
        """
        Retrieve the effective radius on a grid.
            
        Parameters
        ----------
        grid: shdom.Grid, optional
        A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)
            
        Returns
        -------
        reff: shdom.GridData
        A GridData object containing the effective radius [microns] on a grid
        """
        if grid is None:
            grid = self.get_grid()
        
        reff = shdom.GridData(grid,ma.filled(self.reff_data,0.0))
        return reff
    
    
    def get_veff(self, grid=None):
        """
        Retrieve the effective variance on a grid.
            
        Parameters
        ----------
        grid: shdom.Grid, optional
        A shdom.Grid object. If None is specified than a grid is created from Arguments given to the generator (get_grid method)
            
        Returns
        -------
        veff: shdom.GridData
        A GridData object containing the effective variance on a grid
        """
        if grid is None:
            grid = self.get_grid()
        
        veff = shdom.GridData(grid,ma.filled(self.veff_data,0.0))
        return veff
