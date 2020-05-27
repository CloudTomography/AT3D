
from shdom import core
import os
import xarray as xr
import numpy as np

def get_mono_table(particle_type, wavelength_band, minimum_effective_radius=4.0, max_integration_radius=65.0,
             wavelength_averaging=False, wavelength_resolution=0.001, refractive_index=None,
             relative_path=None):
    """
    Mie monodisperse scattering for spherical particles.
    This function will search a given directory to load the requested mie table or will compute it.
    The table is returned as an xarray Dataset.

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
    relative_path: string
        The path to a directory which contains saved mie_table netcdf files. If there is a file with
        'mie_table' in the name that matches the input parameters this file is loaded.
    """

    table_attempt=None
    if relative_path is not None:
        table_attempt = load_table(relative_path,particle_type, wavelength_band,
                    minimum_effective_radius, max_integration_radius,
                     wavelength_averaging, wavelength_resolution,
                     refractive_index)

    if table_attempt is not None:
        table=table_attempt
    else:
        table = compute_table(particle_type, wavelength_band,
                    minimum_effective_radius, max_integration_radius,
                     wavelength_averaging, wavelength_resolution,
                     refractive_index)

    return table

def compute_table(particle_type, wavelength_band,
            minimum_effective_radius, max_integration_radius,
             wavelength_averaging, wavelength_resolution,
             refractive_index):
    """
    This function does the hard work of computing a monomodal mie table.
    See 'get_mono_table' for more details.
    """

    #wavelength band
    wavelen1, wavelen2 = wavelength_band
    assert wavelen1 <= wavelen2 , 'Minimum wavelength is smaller than maximum'

    avgflag = 'C'
    if wavelen1 == wavelen2:
        deltawave = -1
    elif wavelength_averaging:
        avgflag = 'A'
        deltawave = wavelength_resolution

    wavelencen = core.get_center_wavelen(
            wavelen1=wavelen1,
            wavelen2=wavelen2
        )

    #set particle type properties
    if (particle_type == 'Water'):
        rindex = core.get_refract_index(
                partype=particle_type,
                wavelen1=wavelen1,
                wavelen2=wavelen2
            )
        partype = 'W'

    elif (particle_type == 'Aerosol'):
        partype = 'A'
        assert refractive_index is not None, "Refractive index is not specified. \
        This could be a complex number or string path to csv (a function of wavelength)"
        if isinstance(refractive_index, str):
            rindex_df = pd.read_csv('../ancillary_data/dust_volz_1972.ri', comment='#', sep=' ',
                                    names=['wavelength', 'n', 'k'], index_col='wavelength')
            refractive_index = xr.Dataset.from_dataframe(rindex_df).interp({'wavelength': self._wavelencen})
            rindex = np.complex(refractive_index['n'], - refractive_index['k'])
        else:
            rindex = refractive_index
    else:
        raise AttributeError('Particle type note implemented')

    #set mie integration parameters
    if avgflag =='A':
        xmax = 2 * np.pi * max_integration_radius / wavelen1
    else:
        xmax = 2 * np.pi * max_integration_radius / wavelencen
    maxleg = int(np.round(2.0 * (xmax + 4.0 * xmax ** 0.3334 + 2.0)))

    # Set radius integration parameters
    nsize = core.get_nsize(
            sretab=minimum_effective_radius,
            maxradius=max_integration_radius,
            wavelen=wavelencen
    )

    radii = core.get_sizes(
        sretab=minimum_effective_radius,
        maxradius=max_integration_radius,
        wavelen=wavelencen,
        nsize=nsize
    )
    #compute mie properties
    extinct, scatter, nleg, legcoef, table_type = \
        core.compute_mie_all_sizes(
            nsize=nsize,
            maxleg=maxleg,
            wavelen1=wavelen1,
            wavelen2=wavelen2,
            deltawave=deltawave,
            wavelencen=wavelencen,
            radii=radii,
            rindex=rindex,
            avgflag=avgflag,
            partype=partype
        )
    #return data as an xarray
    table = xr.Dataset(
        data_vars={
            'extinction': (['radius'], extinct),
            'scatter': (['radius'], scatter),
            'nleg': (['radius'], nleg),
            'legendre': (['stokes_index', 'legendre_index', 'radius'], legcoef)
        },
        coords={'radius': radii},
        attrs={
            'particle type': particle_type,
            'refractive index': (rindex.real,rindex.imag),
            'refractive index source': str(refractive_index),
            'table type': table_type.decode(),
            'units': ['Radius [micron]','Wavelength [micron]'],
            'wavelength band': wavelength_band,
            'wavelength center': wavelencen,
            'wavelength averaging': str(wavelength_averaging),
            'wavelength resolution': wavelength_resolution,
            'maximum legendre': maxleg,
            'minimum effective radius':minimum_effective_radius,
            'maximum integration radius':max_integration_radius
        },
    )
    return table


def load_table(relative_path,particle_type, wavelength_band,
            minimum_effective_radius=4.0, max_integration_radius=65.0,
             wavelength_averaging=False, wavelength_resolution=0.001,
             refractive_index=None):
    """
    Function that tests whether there is an existing
    mie table that has the specified properties.
    """
    file_list = os.listdir(relative_path)
    table = None

    list_of_input = [particle_type, str(refractive_index), wavelength_band,
                    minimum_effective_radius, max_integration_radius,
                    str(wavelength_averaging), wavelength_resolution,
                    ]
    names = ['particle type', 'refractive index source', 'wavelength band',
            'minimum effective radius', 'maximum integration radius',
            'wavelength averaging', 'wavelength resolution',
            ]

    for file in file_list:

        if file.endswith('.nc') and 'mie_table' in file:
            with xr.open_dataset(os.path.join(relative_path,file)) as dataset:
                #open_dataset reads data lazily.
                test=True

                for name, attr in zip(names, list_of_input):
                    attr_file = dataset.attrs[name]

                    if isinstance(attr,tuple):
                        temp = np.all(attr==attr_file)
                    else:
                        temp = attr==attr_file

                    if not temp:
                        return name, attr, attr_file
                        test=False
                if test:
                    table = file
    if table is not None:
        table = xr.load_dataset(os.path.join(relative_path,table))
    return table
