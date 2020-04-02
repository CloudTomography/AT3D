"""
Medium and Medium related objects used for atmospheric rendering.
A Medium objects is defined by its optical properties: extinction, albedo and phase function.
Two Medium objects can be added to create a combined Medium.
This will result in a unified grid with a summation of the optical extinction and a weighted average of the albedo and phase function.
"""

import numpy as np
import dill as pickle
from collections import OrderedDict
import shdom
from shdom import float_round


class Scatterer(object):
    """
    A Scatterer class to be inherited by specific scatterer types (e.g. OpticalScatter, MicrophysicalScatterer etc...)
    """
    def __init__(self):
        self.grid = None

    @property
    def grid(self):
        return self._grid

    @grid.setter
    def grid(self, val):
        self._grid = val

    @property
    def bounding_box(self):
        return self.grid.bounding_box


class OpticalScatterer(Scatterer):
    """
    An OpticalScatterer class defined in terms of the optical quantities for a given wavelength.

    Parameters
    ----------
    wavelength: float
        A wavelength in microns
    extinction: shdom.GridData object
        A GridData object containing the extinction (1/km) on a grid
    albedo: shdom.GridData
        A GridData object containing the single scattering albedo [0,1] on a grid
    phase: shdom.GridPhase
        A GridPhase object containing the phase function on a grid
    """
    def __init__(self, wavelength, extinction=None, albedo=None, phase=None):
        self._wavelength = wavelength
        super().__init__()
        if (extinction is not None) and (albedo is not None) and (phase is not None):
            self.grid = extinction.grid + albedo.grid + phase.grid
        self.extinction = extinction
        self.albedo = albedo
        self.phase = phase

    def resample(self, grid):
        """
        The resample method resamples the OpticalScatterer (extinction, albedo, phase).

        Parameters
        ----------
        grid: shdom.Grid
            The new grid to which the data will be resampled

        Returns
        -------
        scatterer: shdom.OpticalScatterer
            An optical scatterer resampled onto the input grid
        """
        extinction = self.extinction.resample(grid)
        albedo = self.albedo.resample(grid)
        phase = self.phase.resample(grid)
        return shdom.OpticalScatterer(self.wavelength, extinction, albedo, phase)

    def get_mask(self, threshold):
        """
        Get a mask based on the optical extinction.

        Parameters
        ----------
        threshold: float
            A threshold which above this value it is considered a populated voxel.

        Returns
        -------
        mask: shdom.GridData object
            A boolean mask with True for dense voxels and False for optically thin regions.
        """
        data = self.extinction.data > threshold
        return shdom.GridData(self.grid, data)

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def extinction(self):
        return self._extinction

    @extinction.setter
    def extinction(self, val):
        if val is not None:
            assert val.min_value.astype(np.float32) >= 0.0, 'Extinction min value {} and should be larger than 0.0'.format(val.min_value)
        self._extinction = val

    @property
    def albedo(self):
        return self._albedo

    @albedo.setter
    def albedo(self, val):
        if val is not None:
            assert (val.max_value.astype(np.float32) <= 1.0 and val.min_value.astype(np.float32) >= 0.0), 'Single scattering albedo should be in the range [0, 1]'
        self._albedo = val

    @property
    def phase(self):
        return self._phase

    @phase.setter
    def phase(self, val):
        self._phase = val


class MicrophysicalScatterer(Scatterer):
    """
    A MicrophysicalScatterer class defined in terms of the microphysical quantities.

    Parameters
    ----------
    lwc: shdom.GridData
        A GridData object containing liquid water content (g/m^3) on a 3D grid.
    reff: shdom.GridData
        A GridData object containing effective radii (micron) on a 3D grid.
    veff: shdom.GridData
        A GridData object containing effective variances on a 3D grid.
    """
    def __init__(self, lwc=None, reff=None, veff=None):
        super().__init__()
        self._mie = OrderedDict()
        self._wavelength = []
        self._min_reff = np.Inf
        self._max_reff = -np.Inf
        self._min_veff = np.Inf
        self._max_veff = -np.Inf
        self.set_microphysics(lwc, reff, veff)

    def set_microphysics(self, lwc, reff, veff):
        """
        Set the microphyscal properties and define the scatterer grid.

        Parameters
        ----------
        lwc: shdom.GridData
            A GridData object containing liquid water content (g/m^3) on a 3D grid.
        reff: shdom.GridData
            A GridData object containing effective radii (micron) on a 3D grid.
        veff: shdom.GridPhase
            A GridData object containing effective variances on a 3D grid.

        Notes
        -----
        The scatterer grid contains all three grids of lwc, reff and veff
        (see __add__ method of shdom.Grid for more info)
        """
        self.lwc = lwc
        self.reff = reff
        self.veff = veff
        if (lwc is not None) and (reff is not None) and (veff is not None):
            self._grid = lwc.grid + reff.grid + veff.grid

    def get_extinction(self, wavelength):
        """
        Get the optical extinction from microphysics and Mie scattering model for a given wavelength.

        Parameters
        ----------
        wavelength: float
            Wavelength in microns.

        Returns
        -------
        extinction: shdom.GridData
            A GridData with the optical extinction

        Notes
        -----
        A Mie scattering model must be defined at the input wavelength by use of the add_mie method.
        The input wavelength is rounded to three decimals.
        """
        extinction = self.mie[float_round(wavelength)].get_extinction(self.lwc, self.reff, self.veff)
        return extinction

    def get_albedo(self, wavelength):
        """
        Get the optical single scattering albedo from microphysics and Mie scattering model for a given wavelength.

        Parameters
        ----------
        wavelength: float
            Wavelength in microns.

        Returns
        -------
        albedo: shdom.GridData
            A GridData with the optical single scattering albedo

        Notes
        -----
        A Mie scattering model must be defined at the input wavelength by use of the add_mie method.
        The input wavelength is rounded to three decimals.
        """
        albedo = self.mie[float_round(wavelength)].get_albedo(self.reff, self.veff)
        return albedo

    def get_phase(self, wavelength):
        """
        Get the optical phase function from microphysics and Mie scattering model for a given wavelength.

        Parameters
        ----------
        wavelength: float
            Wavelength in microns.

        Returns
        -------
        phase: shdom.GridPhase
            A GridPhase with the optical phase function

        Notes
        -----
        A Mie scattering model must be defined at the input wavelength by use of the add_mie method.
        The input wavelength is rounded to three decimals.
        """
        phase = self.mie[float_round(wavelength)].get_phase(self.reff, self.veff)
        return phase

    def get_optical_scatterer(self, wavelength):
        """
        Get the optical scatterer out of the microphysics and Mie scattering model for a given wavelength.

        Parameters
        ----------
        wavelength: float of list of floats
            Wavelength in microns.

        Returns
        -------
        scatterer: shdom.OpticalScatterer or shdom.MultispectralScatterer
            An OpticalScatterer object is returned for an input float and a MultispectralScatterer for a list.

        Notes
        -----
        A Mie scattering model must be defined at the input wavelength by use of the add_mie method.
        The input wavelength is rounded to three decimals.
        """
        if isinstance(wavelength, list):
            scatterer_list = [
                shdom.OpticalScatterer(wl, self.get_extinction(wl), self.get_albedo(wl), self.get_phase(wl)) for wl in wavelength
            ]
            scatterer = shdom.MultispectralScatterer(scatterer_list)
        else:
            scatterer = shdom.OpticalScatterer(
                wavelength, self.get_extinction(wavelength), self.get_albedo(wavelength),  self.get_phase(wavelength)
            )
        return scatterer

    def resample(self, grid):
        """
        The resample method resamples the MicrophysicalScatterer in place(lwc, reff, veff).

        Parameters
        ----------
        grid: shdom.Grid
            The new grid to which the data will be resampled

        Notes
        -----
        Unlike the OpticalScatterer resample method, resampling is done in-place since the Mie structures can be large and are left unchanged.
        """
        lwc = self.lwc.resample(grid)
        reff = self.reff.resample(grid)
        veff = self.veff.resample(grid)
        self.set_microphysics(lwc, reff, veff)
        return self

    def get_mask(self, threshold):
        """
        Get a mask based on the liquid water content.

        Parameters
        ----------
        threshold: float
            A threshold which above this value it is considered a populated voxel.

        Returns
        -------
        mask: shdom.GridData object
            A boolean mask with True for dense voxels and False for thin voxels.
        """
        data = self.lwc.data > threshold
        return shdom.GridData(self.grid, data)

    def add_mie(self, mie):
        """
        Add a Mie scattering model at a given wavelength.

        mie: shdom.MiePolydisperse or list of shdom.MiePolydisperse
            Using Mie model microphyical properties are transformed into optical properties (see get_optical_scatterer method)
        """
        mie_list = mie
        if isinstance(mie_list, shdom.MiePolydisperse):
            mie_list = [mie_list]
        elif isinstance(mie, dict):
            mie_list = mie_list.values()

        for mie in mie_list:
            self._mie[mie.wavelength] = mie
            self._wavelength.append(mie.wavelength)
            self._min_reff = min(self.min_reff, mie.size_distribution.reff.min())
            self._max_reff = max(self.max_reff, mie.size_distribution.reff.max())
            self._min_veff = min(self.min_veff, mie.size_distribution.veff.min())
            self._max_veff = max(self.max_veff, mie.size_distribution.veff.max())

        if self.reff is not None:
            min_val = self.reff.data[self.reff.data > 0.0].min()
            assert self.reff.max_value < self.max_reff+1e-3, \
                               'Maximum medium effective radius [{:2.2f}] is larger than the pre-computed table maximum radius [{:2.2f}]. ' \
                               'Recompute Mie table with larger maximum radius.'.format(self.reff.max_value, self.max_reff)
            assert min_val > self.min_reff-1e-3, \
                    'Minimum medium effective radius [{:2.2f}] is smaller than the pre-computed table minimum radius [{:2.2f}]. ' \
                    'Recompute Mie table with smaller minimum radius.'.format(min_val, self.min_reff)

        if self.veff is not None:
            min_val = self.veff.data[self.veff.data>0.0].min()
            assert self.veff.max_value < self.max_veff+1e-3, \
                    'Maximum medium effective variance [{:2.2f}] is larger than the pre-computed table maximum variance [{:2.2f}]. ' \
                    'Recompute Mie table with larger maximum variance.'.format(self.veff.max_value, self.max_veff)
            assert min_val > self.min_veff-1e-3, \
                    'Minimum medium effective variance [{:2.2f}] is smaller than the pre-computed table minimum variance [{:2.2f}]. ' \
                    'Recompute Mie table with smaller minimum variance.'.format(min_val, self.min_veff)

    def load_grid(self, path):
        """
        A utility function to load a microphysical medium from file

        Parameters
        ----------
        path: str
             Path to file.

        Returns
        -------
        grid: shdom.Grid object
            The 3D grid of the medium.

        Notes
        -----
        CSV format should be as follows:

        # comment line (description)
        nx ny nz
        dz dy dz     z_levels[0]     z_levels[1] ...  z_levels[nz-1]
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
        .
        .
        .
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
        """
        nx, ny, nz = np.genfromtxt(path, max_rows=1, dtype=int)
        dx, dy = np.genfromtxt(path, max_rows=1, usecols=(0, 1), dtype=float, skip_header=2)
        z_grid = np.genfromtxt(path, max_rows=1, usecols=range(2, 2 + nz), dtype=float, skip_header=2)
        x_grid = np.linspace(0.0, (nx - 1)*dx, nx, dtype=np.float32)
        y_grid = np.linspace(0.0, (ny - 1)*dy, ny, dtype=np.float32)
        grid = shdom.Grid(x=x_grid, y=y_grid, z=z_grid)
        return grid

    def save_to_csv(self, path, comment_line=''):
        """
        A utility function to save a microphysical medium.

        Parameters
        ----------
        path: str
            Path to file.
        comment_line: str, optional
            A comment line describing the file.


        Notes
        -----
        CSV format is as follows:

        # comment line (description)
        nx ny nz
        dz dy dz     z_levels[0]     z_levels[1] ...  z_levels[nz-1]
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
        .
        .
        .
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]
        """
        np.savetxt(path, X=np.array([self.grid.shape]), fmt='%d', header=comment_line)
        f = open(path, 'ab')
        np.savetxt(f, X=np.concatenate((np.array([self.grid.dx, self.grid.dy]), self.grid.z)).reshape(1,-1), fmt='%2.3f')
        x, y, z = np.meshgrid(range(self.grid.nx), range(self.grid.ny), range(self.grid.nz), indexing='ij')
        idx = self.lwc.data.ravel() > 1e-8
        data = np.vstack((x.ravel()[idx], y.ravel()[idx], z.ravel()[idx], self.lwc.data.ravel()[idx], self.reff.data.ravel()[idx])).T
        np.savetxt(f, X=data, fmt='%d %d %d %.5f %.3f')
        f.close()

    def load_from_csv(self, path, veff=0.1):
        """
        A utility function to load a microphysical medium.

        Parameters
        ----------
        path: str
            Path to file.
        veff: float
            If effective variance is not specified in the csv file as a 6th column,
            this value is used as a homogeneous value.
            Default value is veff=0.1

        Notes
        -----
        CSV format should be as follows:

        # comment line (description)
        nx ny nz
        dz dy dz     z_levels[0]     z_levels[1] ...  z_levels[nz-1]
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]  veff[ix, iy, iz](optional)
        .
        .
        .
        ix iy iz     lwc[ix, iy, iz]    reff[ix, iy, iz]  veff[ix, iy, iz](optional)
        """
        grid = self.load_grid(path)
        data = np.genfromtxt(path, skip_header=3)

        grid_index = data[:, :3].astype(int)
        lwc = data[:, 3]
        reff = data[:, 4]
        if data.shape[1] == 6:
            veff = data[:, 5]
        else:
            veff = veff * np.ones_like(reff)

        lwc_data  = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=np.nan)
        reff_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=np.nan)
        veff_data = np.full(shape=(grid.nx, grid.ny, grid.nz), fill_value=np.nan)
        lwc_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]]  = lwc
        reff_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]] = reff
        veff_data[grid_index[:, 0], grid_index[:, 1], grid_index[:, 2]] = veff

        self.set_microphysics(
            lwc=shdom.GridData(grid, lwc_data).squeeze_dims(),
            reff=shdom.GridData(grid, reff_data).squeeze_dims(),
            veff=shdom.GridData(grid, veff_data).squeeze_dims()
        )

    @property
    def lwc(self):
        return self._lwc

    @lwc.setter
    def lwc(self, val):
        if val is not None:
            assert val.min_value.astype(np.float32) >= 0.0, 'LWC min value is {} and should be larger than 0.0'.format(val.min_value)
        self._lwc = val

    @property
    def reff(self):
        return self._reff

    @reff.setter
    def reff(self, val):
        if val is not None and self.mie:
            max_val = val.max_value
            min_val = val.data[val.data>0.0].min()
            assert  max_val < self.max_reff+1e-3, \
                    'Maximum medium effective radius [{:2.2f}] is larger than the pre-computed table maximum radius [{:2.2f}]. ' \
                    'Recompute Mie table with larger maximum radius.'.format(max_val, self.max_reff)
            assert  min_val > self.min_reff-1e-3, \
                    'Minimum medium effective radius [{:2.2f}] is smaller than the pre-computed table minimum radius [{:2.2f}]. ' \
                    'Recompute Mie table with smaller minimum radius.'.format(min_val, self.min_reff)
        self._reff = val

    @property
    def veff(self):
        return self._veff

    @veff.setter
    def veff(self, val):
        if val is not None and self.mie:
            max_val = val.max_value
            min_val = val.data[val.data>0.0].min()
            assert  max_val < self.max_veff+1e-3, \
                    'Maximum medium effective radius [{:2.2f}] is larger than the pre-computed table maximum variance [{:2.2f}]. ' \
                    'Recompute Mie table with larger maximum radius.'.format(max_val, self.max_veff)
            assert  min_val > self.min_veff-1e-3, \
                    'Minimum medium effective radius [{:2.2f}] is smaller than the pre-computed table minimum variance [{:2.2f}]. ' \
                    'Recompute Mie table with smaller minimum radius.'.format(min_val, self.min_veff)
        self._veff = val

    @property
    def mie(self):
        return self._mie

    @property
    def min_reff(self):
        return self._min_reff

    @property
    def max_reff(self):
        return self._max_reff

    @property
    def min_veff(self):
        return self._min_veff

    @property
    def max_veff(self):
        return self._max_veff

    @property
    def num_wavelengths(self):
        return len(self._wavelength)

    @property
    def wavelength(self):
        if self.num_wavelengths == 0:
            return None
        if self.num_wavelengths == 1:
            return self._wavelength[0]
        else:
            return self._wavelength


class MultispectralScatterer(object):
    """
    A MultispectralScatterer class defined in terms of the optical quantities (extinction, albedo, phase) at several spectral bands.

    Parameters
    ----------
    scatterer_list: list of shdom.OpticalScatterer
        A list of optical scatterers at different wavelengths.
    """
    def __init__(self, scatterer_list=None):
        self._scatterer = OrderedDict()
        self._num_wavelengths = 0
        self._grid = None
        self._wavelength = []
        if scatterer_list is not None:
            for scatterer in scatterer_list:
                self.add_scatterer(scatterer)

    def get_optical_scatterer(self, wavelength):
        """
        Get the optical scatterer at a given wavelength.

        Parameters
        ----------
        wavelength: float
            Wavelength in microns.

        Notes
        -----
        The input wavelength is rounded to three decimals.
        """
        return self.scatterer[float_round(wavelength)]

    def resample(self, grid):
        """
        The resample method resamples the MutlispectralScatterer (extinction, albedo, phase) at all wavelengths.

        Parameters
        ----------
        grid: shdom.Grid
            The new grid to which the data will be resampled

        Returns
        -------
        scatterer: shdom.MultispectralScatterer
            A multi-spectral scatterer resampled onto the input grid
        """
        for wavelength in self.scatterer.keys():
            self.scatterer[float_round(wavelength)] = self.scatterer[float_round(wavelength)].resample(grid)
        return self

    def add_scatterer(self, scatterer):
        """
        Add an optical scatterer to the list.

        Parameters
        ----------
        scatterer: shdom.OpticalScatterer
            An optical scatterer at a given wavelength.
        """
        if self.num_wavelengths == 0:
            self._grid = scatterer.grid
        else:
            self._grid += scatterer.grid
        self._num_wavelengths += 1
        self.scatterer[float_round(scatterer.wavelength)] = scatterer.resample(self.grid)
        self._wavelength.append(scatterer.wavelength)

    @property
    def scatterer(self):
        return self._scatterer

    @property
    def num_wavelengths(self):
        return self._num_wavelengths

    @property
    def grid(self):
        return self._grid

    @property
    def bounding_box(self):
        return self.grid.bounding_box

    @property
    def wavelength(self):
        if len(self._wavelength) == 0:
            return None
        if len(self._wavelength) == 1:
            return self._wavelength[0]
        else:
            return self._wavelength


class Medium(object):
    """
    The Medium object encapsulates an atmospheric optical medium with multiple Scatterers.

    Parameters
    ----------
    grid: shdom.Grid, optional
        A grid for the Medium object. All scatterers will be resampled to this grid.

    Notes
    -----
    The original SHDOM code combines particles at every grid point apriori, to define effective optical properties.
    In contrast, pyshdom uses all scatterers in real-time, to compute the source function and radiance fields during the solution process.
    """
    def __init__(self, grid=None):
        self._scatterers = OrderedDict()
        self._num_scatterers = 0
        self._wavelength = None
        self.set_grid(grid)

    def set_grid(self, grid):
        """
        Set the medium grid.

        Parameters
        ----------
        grid: shdom.Grid
            A grid for the Medium object. All scatterers will be resampled to this grid.
        """
        self._grid = grid

    def add_scatterer(self, scatterer, name=None):
        """
        Add a Scatterer to the medium.

        Parameters
        ----------
        scatterer: shdom.Scatterer
            A scattering particle distribution (e.g. MicrophysicalScatterer, OpticalScatterer)
        name: string, optional
            A name for the scatterer that will be used to retrieve it (see get_scatterer method).
            If no name is specified the default name is scatterer# where # is the number in which it was input (i.e. scatterer1 for the first scatterer).
        """
        first_scatterer = True if self.num_scatterers==0 else False

        if first_scatterer:
            self._wavelength = scatterer.wavelength
        else:
            assert np.allclose(self.wavelength, scatterer.wavelength), ' medium wavelength {} differs from scatterer wavelength {}'.format(self.wavelength, scatterer.wavelength)
        self._num_scatterers += 1
        name = 'scatterer{:d}'.format(self._num_scatterers) if name is None else name
        self.scatterers[name] = scatterer

    def get_scatterer(self, name):
        """
        Retrieve a Scatterer from the medium.

        Parameters
        ----------
        name: string
            The name used for the scatterer on input (see add_scatterer method)

        Returns
        -------
        scatterer: shdom.Scatterer
            The scatterer with the name matching the query
        """
        if isinstance(name, str):
            scatterer = self.scatterers[name]
        elif isinstance(name, list):
            scatterer = {n: self.scatterers[n] for n in name}
        return scatterer

    def save(self, path):
        """
        Save Medium to file.

        Parameters
        ----------
        path: str,
            Full path to file.
        """
        file = open(path,'wb')
        file.write(pickle.dumps(self.__dict__, -1))
        file.close()

    def load(self, path):
        """
        Load Medium from file.

        Parameters
        ----------
        path: str,
            Full path to file.
        """
        file = open(path, 'rb')
        data = file.read()
        file.close()
        self.__dict__ = pickle.loads(data)

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def grid(self):
        return self._grid

    @property
    def scatterers(self):
        return self._scatterers

    @property
    def num_scatterers(self):
        return self._num_scatterers

    @property
    def num_wavelengths(self):
        if isinstance(self._wavelength, list):
            return len(self._wavelength)
        else:
            return 1
