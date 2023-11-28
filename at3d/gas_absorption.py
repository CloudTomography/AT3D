"""
An interface to the REPTRAN absorption parameterization distributed with Libradtran (v2.0.3)
from http://www.libradtran.org/doku.php?id=download.

J. Gasteiger, C. Emde, B. Mayer, R. Buras, S.A. Buehler, O. Lemke,
Representative wavelengths absorption parameterization applied to satellite channels and spectral bands,
Journal of Quantitative Spectroscopy and Radiative Transfer, Volume 148, 2014, Pages 99-115, ISSN 0022-4073,
https://doi.org/10.1016/j.jqsrt.2014.06.024.

This module contains `load_standard_atmosphere` as a starting point to provide atmospheric data.
Other input to the `Reptran` object which is the primary interface should follow the same
data conventions. Currently, the `Reptran` interface only works with MODIS band models.
Gas absorption data are found in ./data/reptran/. The coarse band data are also distributed
but an interface is not yet implemented.
"""

import numpy as np
import xarray as xr
import os
import importlib
from datetime import datetime

def load_standard_atmosphere(search_directory=None, atmosphere_name='tropical'):
    """
    Reads gas concentrations from AFGL standard atmospheres.

    Standard trace gases are supplemented by extra gas concentrations from the US
    standard atmosphere (N20, CH4, CO, N2).

    Parameters
    ----------
    search_directory : str
        The directory where the atmosphere files are.
    atmosphere_name : str
        The name of the standard atmosphere to use.
        Options include: 'midlatitude_summer', 'midlatitude_winter', 'tropical',
        'subarctic_summer', 'subarctic_winter'

    Returns
    -------
    atmosphere : xr.Dataset
        Contains the gas concentrations in units of cm^-3 as well as
        temperature in Kelvin, pressure in millibars and altitude ('z')
        in kilometers.
    """

    file_name_dict = {
        'midlatitude_summer': 'afglms.txt',
        'midlatitude_winter': 'afglmw.txt',
        'tropical': 'afglt.txt',
        'subarctic_summer': 'afglss.txt',
        'subarctic_winter': 'afglsw.txt'
    }
    try:
        file_name = file_name_dict[atmosphere_name]
    except KeyError:
        raise KeyError("Invalid name for standard atmosphere. `atmosphere_name` should be one of {}".format(tuple(file_name_dict.keys())))

    if search_directory is None:
        search_directory = os.path.join(importlib.resources.files('at3d'),'data/ancillary/')

    main_path = os.path.join(search_directory, file_name)

    atmosphere_data = np.loadtxt(main_path)
    atmosphere = xr.Dataset(
        data_vars={
            'temperature': ('pressure', atmosphere_data[:, 2]),
            'air': ('pressure', atmosphere_data[:, 3]),
            'O3': ('pressure', atmosphere_data[:, 4]),
            'O2': ('pressure', atmosphere_data[:, 5]),
            'H2O': ('pressure', atmosphere_data[:, 6]),
            'CO2': ('pressure', atmosphere_data[:, 7]),
            'NO2': ('pressure', atmosphere_data[:, 8]),
        },
        coords={
            'pressure': ('pressure', atmosphere_data[:, 1]*100),
            'z': ('z', atmosphere_data[:, 0]),
        }
    )
    extra_names = ('N2O', 'CH4', 'CO', 'N2')
    for extra_name in extra_names:
        name = os.path.join(search_directory, 'afglus_{}_vmr.dat'.format(extra_name.lower()))
        atmosphere[extra_name] = (
            'pressure',
            np.loadtxt(name)[:, -1]*atmosphere.air.data
        )
    return atmosphere

class Reptran:
    """
    Loads gas properties from the REPTRAN gas absorption parameterization.

    Parameters
    ----------
    search_directory : str
        The directory where the reptran files are.
    parameterization_name : str
        The name of the parameterization. May be a band model in the thermal or solar
        (e.g. thermal_fine) or may be an instrument model (e.g. solar_modis).

    Reference
    ---------
    https://doi.org/10.1016/j.jqsrt.2014.06.024
    """
    def __init__(self, parameterization_name='solar_modis', search_directory=None):

        if search_directory is None:
            search_directory = os.path.join(importlib.resources.files('at3d'),'data/reptran/')

        self.search_directory = search_directory

        try:
            self.reference = xr.open_dataset(os.path.join(self.search_directory, 'reptran_{}.cdf'.format(parameterization_name)))
        except FileNotFoundError:
            raise FileNotFoundError("No REPTRAN data found for `reference_name` = {}".format(parameterization_name))

        self.parameterization_name = parameterization_name

        self.gas_property_files = {}
        for name in self.reference.species_name:
            decoded = str(name.data.astype(str)).strip()
            gas_file = xr.load_dataset(
                os.path.join(self.search_directory,
                             'reptran_{}.lookup.{}.cdf'.format(parameterization_name, decoded))
            )

            # we change dimension names and assign coordinates so we can use xarray interpolation routines.
            dimension_dict = {'n_pressure': 'pressure', 'n_t_pert': 't_pert', 'n_vmrs': 'wvmr'}
            gas_file = gas_file.swap_dims(dimension_dict)
            coord_dims = {'pressure': gas_file.pressure, 't_pert': gas_file.t_pert, 'wvmr': gas_file.vmrs,
                         }
            gas_file = gas_file.assign_coords(coord_dims)
            self.gas_property_files[decoded] = gas_file

    def get_absorption_data(self, atmosphere, wavelength=None, instrument=None, band=None):
        """
        Calculates the volume extinction coefficient at all representative wavelengths
        for the specified band.

        If a band model is used rather than a particular instrument response function
        then specify the wavelength and you will be given the gas properties of the
        associated band. The resulting calculation then becomes pseudo-monochromatic.

        Parameters
        ---------
        atmosphere : xr.Dataset
            Contains the required gas concentrations in units of cm^-3.
            See `load_standard_atmosphere` for details.
        wavelength : float
            The wavelength in micrometers for which to calculate gas properties.
        instrument : str
            The name of the instrument (e.g. 'terra' for MODIS).
        band : int
            The band index (e.g. 26 for band 26 from MODIS).

        Returns
        -------
        atmosphere : xr.Dataset
            A deep copy of `atmosphere` with the associated volume extinction
            coefficients added for each representative wavelength along with
            associated weights.
        """

        atmosphere = atmosphere.copy(deep=True)

        instrument_model = True
        for name in ('coarse', 'medium', 'fine'):
            if name in self.parameterization_name:
                instrument_model = False

        if instrument_model:
            if band is None:
                raise ValueError(
                    "The band index must be specified for this gas absorption parameterization."
                )

            if ('modis' in self.parameterization_name) & (instrument is None):
                raise ValueError(
                    "The MODIS instrument name (`terra` or `aqua`) must be specified "
                    "for the MODIS band parameterization."
                )
            
            if not isinstance(band, int):
                raise TypeError(
                    "`band` should be an integer."
                )
            absorption_data = self._get_band_absorption_data(atmosphere, instrument, band)

        else:
            if wavelength is None:
                raise ValueError(
                    "The desired monochromatic wavelength around which to get absorption data "
                    "must be specified for '{}'.".format(self.parameterization_name)
                )
            absorption_data = self._get_wavelength_absorption_data(atmosphere, wavelength)

        return absorption_data


    def _get_band_absorption_data(self, atmosphere, instrument, band):
        """
        Returns the gas absorption data as volume coefficient for the specified
        atmosphere.

        Parameters
        ----------
        band : int
            The index of MODIS's bands. Valid indices are (1-20) and 26.
        atmosphere : str
            The filename for an atmosphere (following the format of afglt.dat)

        Returns : xr.Dataset
            Contains the `atmophere` data as well as the gas absorption properties
            in units of 1/km at each of the wavelength quadrature points as well
            as the corresponding weights.
        """
        short_name = 'modis_'+instrument +'_b{}'.format(str(band).zfill(2))
        band_name = (short_name).encode('utf-8').ljust(500)
        try:
            band_index = np.where(self.reference.band_name == band_name)[0][0]
        except IndexError:
            raise ValueError(
                "`band` {} is not a valid name".format(band)
            )

        return self._process_data(atmosphere, short_name, band_index)

    def _get_wavelength_absorption_data(self, atmosphere, wavelength):

        # convert to nanometers from micrometers.
        wavelength *= 1e3

        if (wavelength < self.reference.wvlmin.data.min()) or (wavelength > self.reference.wvlmax.data.max()):
            raise ValueError("`wavelength` '{}' is out of bounds of the reptran parameterization.".format(wavelength))

        band_index = np.digitize(wavelength,
                                 bins=np.append(self.reference.wvlmin.data, self.reference.wvlmax[-1])
                                ) - 1
        short_name = str(self.reference.band_name[band_index].data.astype(str)).strip()

        return self._process_data(atmosphere, short_name, band_index)

    def _process_data(self, atmosphere, short_name, band_index):
        """
        Queries the Look Up Tables (LUT) of gas absorption cross sections
        for each species and adds their contributions to the extinction
        at each representative wavelength.
        """
        iwvl_band = self.reference.iwvl[:, band_index]
        iwvl_band = iwvl_band[np.where(iwvl_band > 0)]

        cross_section_sources = self.reference.cross_section_source[iwvl_band-1, :]
        band_wvls = self.reference.wvl[iwvl_band-1]

        volume_extinctions = []
        # Loop through the different wavelength quadrature points.

        for wvl_index, wavelength_source in zip(iwvl_band, cross_section_sources):

            volume_extinction = np.zeros(atmosphere.z.size)

            for species_index, species_flag in enumerate(wavelength_source):
                if species_flag: # If this species is needed at this wavelength quadrature point.
                    species_name = str(self.reference.species_name[species_index].data.astype(str)).strip()
                    species_file = self.gas_property_files[species_name]

                    wavelength_index = np.where(species_file.wvl_index == wvl_index)[0]

                    # coordinates for the gas absorption cross sections is temperature perturbation (tpert)
                    # from a reference profile as well as pressure (and mixing ratio for H2O).
                    tpert = atmosphere.temperature - species_file.t_ref.interp(pressure=atmosphere.pressure)

                    interp_coords = {'t_pert':tpert, 'pressure':atmosphere.pressure}
                    to_interpolate = species_file.xsec[:, :, wavelength_index[0], :]

                    # special case for H2O which also is sensitive to mixing ratio
                    if species_name == 'H2O':
                        interp_coords['wvmr'] = atmosphere[species_name]/atmosphere.air
                    else:
                        to_interpolate = to_interpolate[:, 0] # remove the useless 'wvmr' dimension

                    cross_section = to_interpolate.interp(interp_coords)

                    #1e5 converts to km # everything here in /cm^3 and cm^2
                    volume_extinction += 1e5* (atmosphere[species_name]*cross_section*1e-16)

            volume_extinctions.append(volume_extinction)

        dimensions = ['wavelength', 'z']
        other_dims = list(atmosphere.air.dims).remove('pressure')
        if other_dims is not None:
            dimensions += other_dims

        atmosphere['gas_absorption'] = (dimensions, np.stack(volume_extinctions, axis=0))
        atmosphere['weights'] = ('wavelength', self.reference.iwvl_weight[:len(iwvl_band), band_index].data)

        if 'solar' in self.parameterization_name:
            atmosphere['irradiance'] = ('wavelength', self.reference.extra[iwvl_band-1].data)
            atmosphere.irradiance.attrs['description'] = 'Extraterrestrial solar flux at representative wavelengths'
            atmosphere.irradiance.attrs['units'] = 'mW / (m2 nm)'
        atmosphere['spectral_integral'] = self.reference.wvl_integral[band_index].data
        atmosphere.spectral_integral.attrs['description'] = 'Integral of band response function over wavelength'
        atmosphere.spectral_integral.attrs['units'] = 'nm'
        atmosphere = atmosphere.assign_coords({'wavelength': band_wvls.data[:len(iwvl_band)]*1e-3})

        atmosphere.gas_absorption.attrs['description'] = 'Volume Extinction Coefficient in units of 1/km derived from REPTRAN'
        atmosphere['band_name'] = short_name
        atmosphere.attrs['parameterization_name'] = self.parameterization_name
        atmosphere.attrs['description'] = 'Gas absorption derived from REPTRAN (https://doi.org/10.1016/j.jqsrt.2014.06.024) as distributed from http://www.meteo.physik.uni-muenchen.de/~libradtran/lib/exe/fetch.php?media=download:reptran_2017_all.tar.gz.'
        atmosphere['created'] = datetime.now()

        return atmosphere
