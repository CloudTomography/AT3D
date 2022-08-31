"""This module loads aerosol properties.

The aerosol properties available are from OPAC. They are generated from
translated scattering tables from libRadtran.
See the ../data/OPAC/aerosol/convert_libradtran_phase_to_at3d.py script for
details on the conversion.
Several different default mixtures of aerosol components are available with
default vertical profiles. Some of the species respond to humidity.

The reference for OPAC is:

Hess, M., Koepke, P., & Schult, I. (1998). Optical Properties of Aerosols
and Clouds: The Software Package OPAC, Bulletin of the American Meteorological
Society, 79(5), 831-844. Retrieved Aug 21, 2022, from
https://journals.ametsoc.org/view/journals/bams/79/5/1520-0477_1998_079_0831_opoaac_2_0_co_2.xml

Currently these aerosol species are not easily optimizeable. In principle,
the `density` field could be retrieved after assuming a certain aerosol component
but that is not a flexible way to handle mixtures of aerosol during retrieval.

See data/OPAC/aerosol/README.md for more info.
"""
import glob
import os

import numpy as np
import xarray as xr

import at3d.medium
import at3d.checks

class MassToNumberConverter:
    """
    Translates mass concentration profiles for the aerosol
    components in each mixture used in libRadtran to number concentration
    profiles so that the relative abundance of each component is correctly
    maintained when humidity is varied from the reference amount.

    Parameters
    ----------
    directory : str
        The directory to the file size_distr.cfg which states the
        size distribution parameters used to generate each aerosol component.
        The default assumes the code is called from the same directory
        as this file.
    """
    def __init__(self, directory='../data/OPAC/aerosol/size_distr.cfg'):

        self._directory = directory
        self._size_config = np.loadtxt(self._directory, skiprows=2)

        self._config_indices = {
            'inso': 1,
            'waso': 2,
            'soot': 3,
            'ssam': 4,
            'sscm': 5,
            'minm': 6,
            'miam': 7,
            'micm': 8,
            'mitr': 9,
            'suso': 10
        }

        self._humidities = np.array([0.0, 50.0, 70.0, 80.0, 90.0, 95.0, 98.0, 99.0])/100.0

        self._look_up_tables = {key: None for key in self._config_indices}

    def __call__(self, species_name, humidity):
        """Calculates a conversion factor for the
        """

        if species_name not in self._look_up_tables:
            raise ValueError(
                "Unknown species name `{}`".format(species_name)
            )

        if self._look_up_tables[species_name] is None:
            self.make_look_up_table(species_name)

        conversion_factors = self._look_up_tables[species_name].interp(
            humidity=humidity,
            method='linear',
            kwargs={'fill_value':'extrapolate'}
        )
        return conversion_factors


    def make_look_up_table(self, species_name):
        """Prepares the conversion factors for all humidities for
        `species_name`.

        Parameters
        ----------
        species_name : str
            The shortened version of the species name that matches the file names
            in ../data/OPAC/aerosol/optical_properties/
            E.g. inso, waso. See also self._config_indices
        """

        conversion_factors = np.array(
            [self._get_conversion_factor(species_name, humidity)
             for humidity in self._humidities]
        )

        # divide by the reference at 50% humidity
        conversion_factors /= self._get_conversion_factor(species_name, 0.5)

        dset = xr.DataArray(
            data=conversion_factors,
            dims=['humidity'],
            coords={
                'humidity': self._humidities
            }
        )
        self._look_up_tables[species_name] = dset

    def _get_conversion_factor(self, species_name, humidity):
        """Calculates the mass to number concentration ratio
        for the size distribution parameters specified in size_distr.cfg.

        Parameters
        ----------
        species_name : str
            The shortened version of the species name that matches the file names
            in ../data/OPAC/aerosol/optical_properties/
            E.g. inso, waso. See also self._config_indices
        humidity : float
            The relative humidity as a fraction (e.g. 0.5).

        Returns
        -------
        conversion_factor : float
            The number concentration in particles per cubic centimeter for
            the particular size distribution normalized to unit mass concentration
            (1 g/m^3).
        """

        config = self._size_config[np.where(self._size_config[:, 0] == self._config_indices[species_name])[0], 1:]

        if config.shape[0] != 1:
            humidity_index = np.where(config[:, 0] == humidity*100)
            config = config[humidity_index]

        rmin, rmax, rmod, rho, sigma = config[0, 1:]

        r = np.linspace(rmin, rmax, 1000)
        delr = np.sqrt(r*np.append(r[1:], r[-1])) - np.sqrt(r*np.append(r[0], r[:-1]))

        prefactor = 1.0/(np.sqrt(2*np.pi)*r*np.log(sigma))
        lognormal_density = prefactor*np.exp(-0.5*((np.log(r) - np.log(rmod))/(np.log(sigma)))**2)

        lognormal_density /= np.sum(delr*lognormal_density)
        vol = (4/3.0)*np.pi*np.sum(delr*lognormal_density*(r)**3)

        conversion_factor = (rho*vol*(1e-4**3)*1e6)

        return conversion_factor

def interp_optical_property_table(table, method='cubic', **coordinates):
    """Interpolates scattering properties in a table.

    Interpolates the optical property table by interpolating extinction efficiency,
    scattering and also the scattering x phase function phase.

    Parameters
    ----------
    table : xr.Dataset
        A valid scattering table. See at3d.checks.check_legendre.
    coordinates : np.array
        Coordinates in the `table` xr.Dataset that are used for interpolation.
        Should be 1D arrays. For example wavelength or humidity.

    Returns
    -------
    new_table : xr.Dataset
        The interpolated table.

    Notes
    -----
    This is only designed to interpolate the OPAC tables in wavelength
    but might be useful elsewhere but hasn't been tested.
    """
    for name in coordinates:
        if name not in table.coords:
            raise KeyError(
                "Interpolation coordinate, '{}' is not a coordinate in `table`.".format(
                    name
                )
            )

    extinct = table.extinction.interp(**coordinates, method=method)
    scat = (table.extinction*table.ssalb).interp(**coordinates, method=method)
    scat_phase = (table.extinction*table.ssalb*table.legcoef).interp(**coordinates, method=method)

    ssalb = scat/extinct
    ssalb.name = 'ssalb'

    legcoef = scat_phase/scat
    legcoef = legcoef.transpose('stokes_index', 'legendre_index', ...)
    legcoef.name = 'legcoef'

    new_table = xr.merge([extinct, ssalb, legcoef])
    new_table = new_table.assign_attrs(table.attrs)

    at3d.checks.check_legendre(new_table)

    return new_table

class OPACMixture:
    """
    Calculates the optical properties for different mixtures of OPAC
    aerosol types.

    Parameters
    ----------
    directory : str
        The directory to the aerosol data which includes both the
        standard profiles and also the scattering tables.
        The default assumes the code is called from the same directory
        as this file.
    """
    def __init__(self, directory='../data/OPAC/aerosol/'):

        self._directory = directory

        self._expected_types = (
            'antarctic',
            'continental_average',
            'continental_clean',
            'continental_polluted',
            #'desert_spheroids',
            'desert',
            'maritime_clean',
            'maritime_polluted',
            'maritime_tropical',
            'urban'
        )

        self._types = glob.glob(os.path.join(self._directory, 'standard_aerosol_files/*.dat'))

        dat_types = [t+'.dat' for t in self._expected_types]

        for name in dat_types:
            not_found = True
            for file_name in self._types:
                name_long = os.path.join(self._directory, 'standard_aerosol_files', name)
                if name_long == file_name:
                    not_found = False
            if not_found:
                raise IOError(
                    "Could not find OPAC aerosol file with name '{}' in '{}'. Please respecify "
                    "the directory.".format(
                        name,
                        os.path.join(self._directory, 'standard_aerosol_files')
                    )
                )

        self._component_files = {
            'water_soluble': 'waso.mie.cdf',
            'sea_salt_accumulation': 'ssam.mie.cdf',
            'sea_salt_coarse': 'sscm.mie.cdf',
            'soot': 'soot.mie.cdf',
            'sulfate': 'suso.mie.cdf',
            'insoluble': 'inso.mie.cdf',
            'mineral_nucleation': 'minm.mie.cdf',
            'mineral_accumulation': 'miam.mie.cdf',
            'mineral_coarse': 'micm.mie.cdf',
            'mineral_transport': 'mitr.mie.cdf'
        }

        self._aerosol_mass_converter = MassToNumberConverter(directory=os.path.join(self._directory,'size_distr.cfg'))

        self._optprop_tables = {}

    def __call__(self, atmosphere, wavelengths, mixture_type='maritime_clean', reference_aod=None, reference_wavelength=0.55):
        """Calculates the optical properties for the desired default aerosol mixture type.

        Parameters
        ----------
        atmosphere : xr.Dataset
            A valid at3d.grid dataset with a 'humidity' variable.
            See `at3d.grid.make_rte_grid`.
        wavelengths : np.array
            A 1D array of wavelengths in micrometers to calculate the optical properties at.
        mixture_type : str
            The type of OPAC aerosol mixture.
        reference_aod : float
            The aerosol optical depth to scale the optical properties to have.
        reference_wavelength : float
            The wavelength in micrometers at which the `reference_aod` will be prescribed.

        Returns
        -------
        aerosol_optical_properties : dict
            The aerosol optical properties on the grid of `atmosphere` with wavelengths
            as keys.
        """
        if any(wavelengths < 0.25) or any(wavelengths > 40.0):
            raise ValueError(
                "`wavelengths` are outside the expected range of 0.25 to 40 micrometers. Check units of input."
            )

        at3d.checks.check_grid(atmosphere)

        # `atmosphere` should have a humidity variable on the rte_grid.
        at3d.checks.check_hasdim(atmosphere, humidity=('x', 'y', 'z'))

        # Load mixture profile and load any required optical property tables for aerosol types.
        standard_profile = self.get_profile(mixture_type)

        humidified_profile = standard_profile.interp(z=atmosphere.z)

        for aerosol_type in standard_profile:
            conversion_factors = self._aerosol_mass_converter(aerosol_type, atmosphere.humidity)
            humidified_profile[aerosol_type] = conversion_factors*humidified_profile[aerosol_type]

        humidified_profile = humidified_profile.drop(['humidity'])
        aerosol = atmosphere.copy(deep=True)
        for aerosol_type in humidified_profile:
            aerosol[aerosol_type] = humidified_profile[aerosol_type]

        aerosol_optical_properties = {}

        tables = {
            aerosol_type: interp_optical_property_table(
                self._optprop_tables[aerosol_type], wavelength=wavelengths, method='linear'
                )
            for aerosol_type in humidified_profile}

        for wavelength in wavelengths:

            all_species_optical_props = []

            for aerosol_type in humidified_profile:

                aerosol['density'] = aerosol[aerosol_type]
                aerosol['humidity'] = atmosphere.humidity*100.0

                table = tables[aerosol_type].sel(
                    wavelength=wavelength,
                    method='nearest'
                    ).drop(['wavelength'])

                table.coords['table_index'] = ('humidity', np.arange(table.humidity.size))

                if table.humidity.size == 1:
                    table = table.sel({'humidity': 0.0}).drop(['humidity'])

                species_opt_prop = at3d.medium.table_to_grid(aerosol, table)

                all_species_optical_props.append(species_opt_prop)

            merged_optical_properties = at3d.medium.merge_optical_properties(
                *all_species_optical_props
                )
            merged_optical_properties.attrs['wavelength_center'] = wavelength
            merged_optical_properties.attrs['description'] = "OPAC Aerosol with mixture type {}. See at3d.aerosol.AerosolMixture for details.".format(mixture_type)

            aerosol_optical_properties[wavelength] = merged_optical_properties

        if reference_aod is not None:
            self.scale_aod(
                aerosol_optical_properties,
                atmosphere,
                mixture_type,
                reference_aod=reference_aod,
                reference_wavelength=reference_wavelength
            )

        return aerosol_optical_properties

    def scale_aod(self, aerosol_optical_properties, atmosphere, mixture_type, reference_aod, reference_wavelength=0.55):
        """
        Scale the aerosol optical properties so that they have a specified aod at a particular wavelength.

        Parameters
        ----------
        aerosol_optical_properties : dict
            Created by OPACMixture.__call__
        atmosphere : xr.Dataset
            A valid at3d.grid dataset with a 'humidity' variable.
            See `at3d.grid.make_rte_grid`.
        mixture_type : str
            The type of OPAC aerosol mixture.
        reference_aod : float
            The aerosol optical depth to scale the optical properties to have.
        reference_wavelength : float
            The wavelength in micrometers at which the `reference_aod` will be prescribed.
        """
        if reference_wavelength in aerosol_optical_properties:
            reference = aerosol_optical_properties[reference_wavelength]
        else:
            reference = self(atmosphere, np.atleast_1d(reference_wavelength), mixture_type)[reference_wavelength]

        aod = (0.5*(reference.extinction[..., 1:].data + reference.extinction[..., :-1].data)*(reference.z[..., 1:].data - reference.z[..., :-1].data)).mean(axis=(0, 1)).sum(axis=-1)
        ratio = reference_aod/aod

        for aerosol_opt in aerosol_optical_properties.values():
            aerosol_opt.extinction.values *= ratio

        return aerosol_optical_properties

    def get_profile(self, mixture_type):
        """
        Load the mass concentrations (g/m^3) of a given standard aerosol profile.

        Profiles are at the fixed humidity of 50%.

        Parameters
        ----------
        mixture_type : str
            The name of the standard aerosol profile mixture.

        Returns
        -------
        profile : xr.Dataset
            The aerosol profile data loaded from the file with abbreviated names
            for each species. See self._component_files for names.

        Notes
        -----
        Long names of aerosol components are TODO as these profiles are
        expected to mostly used internally by OPACMixture.__call__
        """

        if mixture_type not in self._expected_types:
            raise  ValueError(
                "'{}' is not a valid `mixture_type`. Supported values are"
                " {}".format(mixture_type, self._expected_types)
            )

        if not mixture_type.endswith('.dat'):
            mixture_type += '.dat'

        mixture_long = os.path.join(self._directory, 'standard_aerosol_files', mixture_type)

        profile_file = np.loadtxt(mixture_long)

        with open(mixture_long, "r") as f:
            lines = f.readlines(500)

        types = [string.strip() for string in lines[8][:].split(' ')
                 if string.strip() not in ('', '#')
                 ]

        profile = xr.Dataset(
            data_vars={
                atype: (
                ['z'], mass_conc) for atype, mass_conc
                in zip(types[1:], profile_file[:, 1:].T
                )
            },
            coords={
                'z': profile_file[:, 0]
            }
        )

        for aerosol_type in types[1:]:
            if aerosol_type not in self._optprop_tables:
                self._optprop_tables[aerosol_type] = xr.load_dataset(
                    os.path.join(self._directory, 'optical_properties', aerosol_type+'.mie.nc')
                )

        return profile
