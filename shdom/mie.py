from shdom import core
import os
import xarray as xr
import numpy as np
import pandas as pd

def get_mono_table(particle_type, wavelength_band, minimum_effective_radius=4.0, max_integration_radius=65.0,
                   wavelength_averaging=False, wavelength_resolution=0.001, refractive_index=None,
                   relative_path=None, verbose=True):
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
        table_attempt = _load_table(relative_path,particle_type, wavelength_band,
                    minimum_effective_radius, max_integration_radius,
                     wavelength_averaging, wavelength_resolution,
                     refractive_index)

    if table_attempt is not None:
        table=table_attempt
    else:
        if verbose:
            print('making mie_table. . . may take a while.')
        table = _compute_table(particle_type, wavelength_band,
                    minimum_effective_radius, max_integration_radius,
                     wavelength_averaging, wavelength_resolution,
                     refractive_index, verbose=verbose)

    return table

def _compute_table(particle_type, wavelength_band,
                  minimum_effective_radius, max_integration_radius,
                  wavelength_averaging, wavelength_resolution,
                  refractive_index, verbose=True):
    """
    This function does the hard work of computing a monomodal mie table.
    See 'get_mono_table' for more details.
    """
    #wavelength band
    wavelen1, wavelen2 = wavelength_band
    if wavelen1 > wavelen2:
        raise ValueError('wavelen1 must be <= wavelen2')

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
        # TODO : debug aerosol particle tpye
        partype = 'A'
        assert refractive_index is not None, "Refractive index is not specified. \
        This could be a complex number or string path to csv (a function of wavelength)"
        if isinstance(refractive_index, str):
            rindex_df = pd.read_csv('../ancillary_data/dust_volz_1972.ri', comment='#', sep=' ',
                                    names=['wavelength', 'n', 'k'], index_col='wavelength')
            refractive_index = xr.Dataset.from_dataframe(rindex_df).interp({'wavelength': wavelencen})
            rindex = np.complex(refractive_index['n'], - refractive_index['k'])
        else:
            rindex = refractive_index
    else:
        raise AttributeError('Particle type not implemented')

    #set integration parameters
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
            partype=partype,
            verbose=verbose
        )
    #return data as an xarray
    table = xr.Dataset(
        data_vars={
            'extinction': (['radius'], extinct),
            'scatter': (['radius'], scatter),
            'nleg': (['radius'], nleg),
            'legendre': (['stokes_index', 'legendre_index', 'radius'], legcoef)
            },
        coords={
            'radius': radii,
            'stokes_index': (['stokes_index'], ['P11','P22','P33','P44','P12','P34'])
            },
        attrs={
            'particle_type': particle_type,
            'refractive_index': (rindex.real, rindex.imag),
            'refractive_index_source': str(refractive_index),
            'table_type': table_type.decode(),
            'units': ['Radius [micron]', 'Wavelength [micron]'],
            'wavelength_band': wavelength_band,
            'wavelength_center': wavelencen,
            'wavelength_averaging': str(wavelength_averaging),
            'wavelength_resolution': wavelength_resolution,
            'maximum_legendre': maxleg,
            'minimum_effective_radius':minimum_effective_radius,
            'maximum_integration_radius':max_integration_radius
            },
        )
    return table


def _load_table(relative_path,particle_type, wavelength_band,
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
    names = ['particle_type', 'refractive_index_source', 'wavelength_band',
            'minimum_effective_radius', 'maximum_integration_radius',
            'wavelength_averaging', 'wavelength_resolution',
            ]

    for file in file_list:
        try:
            if file.endswith('.nc'):
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
                            test=False
                    if test:
                        table = file
        except:
            continue
    if table is not None:
        table = xr.load_dataset(os.path.join(relative_path,table))
    return table

def get_poly_table(size_distribution, mie_mono_table):
    """
    TODO
    Calculates mie scattering table for a polydisperse size distribution.
    """
    nd = size_distribution['number_density'].values.reshape((len(size_distribution['radius'])),-1)
    #TODO
    #add interpolation onto radius
    assert np.all(size_distribution.coords['radius'] == mie_mono_table.coords['radius']), 'radii should be consistent between size distribution and mie_mono_table'
    extinct, ssalb, nleg, legcoef = \
        core.get_poly_table(
            nd=nd,
            ndist=nd.shape[-1],
            nsize=mie_mono_table.coords['radius'].size,
            maxleg=mie_mono_table.attrs['maximum_legendre'],
            nleg1=mie_mono_table['nleg'],
            extinct1=mie_mono_table['extinction'],
            scatter1=mie_mono_table['scatter'],
            legcoef1=mie_mono_table['legendre'])

    grid_shape = size_distribution['number_density'].shape[1:]

    #all coords except radius
    coords = {name:coord for name, coord in size_distribution.coords.items() if name not in ('radius', 'stokes_index') }
    coord_lengths = [np.arange(coord.size) for name, coord in coords.items()]
    legen_index = np.meshgrid(*coord_lengths, indexing='ij')

    #TODO: Does this need + 1 here?
    table_index = np.ravel_multi_index(legen_index, dims=[coord.size for coord in coord_lengths])
    coords['table_index'] = (list(size_distribution.coords.keys())[1:], table_index)
    coords['stokes_index'] = mie_mono_table.coords['stokes_index']

    poly_table = xr.Dataset(
            data_vars = {
                'extinction': (list(size_distribution.coords.keys())[1:], extinct.reshape(grid_shape)),
                'ssalb': (list(size_distribution.coords.keys())[1:], ssalb.reshape(grid_shape)),
                'legcoef': (['stokes_index','legendre_index'] + list(size_distribution.coords.keys())[1:],
                            legcoef.reshape(legcoef.shape[:2] + grid_shape)),},
        coords=coords
    )
    poly_table = poly_table.assign_attrs(size_distribution.attrs)
    poly_table = poly_table.assign_attrs(mie_mono_table.attrs)
    return poly_table
