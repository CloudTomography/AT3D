"""This module contains functions to generate Mie scattering properties.

Python binding are implemented for the native SHDOM Mie module for computations of
monodisperse, i.e., single particle scattering properties. Supported particle types are
'water', 'aerosol' (in which case the index of refraction is user specified) or
'ice' (in which case index of refraction depends on wavelength).
The module also has additional support for integration over arbitrary particle distributions.

For water or ice particles the scattering properties may be averaged over the
desired spectral range with Planck function weighting. The phase functions in the
output scattering table are represented with Legendre series. For polarized output,
the six unique elements of the phase matrix are represented with Wigner d-function
expansions (Doicu et al., 2013, JQSRT, http://dx.doi.org/10.1016/j.jqsrt.2012.12.009).
"""
import os
import xarray as xr
import numpy as np

import at3d.core
import at3d.checks

def get_mono_table(particle_type, wavelength_band, minimum_effective_radius=4.0,
                   max_integration_radius=65.0, wavelength_averaging=False,
                   wavelength_resolution=0.001, refractive_index=None,
                   relative_dir=None, verbose=True):
    """
    Mie monodisperse scattering for spherical particles.
    This function will search a given directory to load the requested mie table or will compute it.
    The table is returned as an xr Dataset.

    Parameters
    ----------
    particle_type: string
        Options are 'Water' or 'Aerosol'.
    wavelength_band: (float, float)
        (minimum, maximum) wavelength in microns.
        This defines the spectral band over which to integrate, if both are equal
        monochrome quantities are computed.
    minimum_effective_radius: float
        Minimum effective radius in microns. Used to compute minimum radius for integration.
    max_integration_radius: float
        Maximum radius in microns - cutoff for the size distribution integral
    wavelength_averaging: bool
        True - average scattering properties over the wavelength_band.
        False - scattering properties of the central wavelength.
    wavelength_resolution: float
        The distance between two wavelength samples in the band. Used only if
        wavelength_averaging is True.
    refractive_index: complex number or None
        For 'Water' the refractive index is ignored and loaded from a table
        in src/polarized/indexwatice. For 'Aerosol' the refractive index should
        have a negative imaginary part. ri = n - ik
        For water with custom refractive index use 'Aerosol' particle type.
    relative_dir: string
        The path to a directory which contains saved mie_table netcdf files. If there is a file with
        'mie_table' in the name that matches the input parameters this file is loaded.
    verbose: bool
        True for progress prints from the fortran computations.

    Returns
    -------
    table: xr.Dataset
        A Dataset containing tabulated Legendre coefficients per radius per wavelength.
        Each entry has 6 six Wigner d-function elements of the Mie phase matrix
        indexed by: P11, P22, P33, P44, P12, P34.

    Notes
    -----
    Currently supports only 'Water' and 'Aerosol' particles ('Ice' is not yet implemented).

    References
    ----------
    .. [1] Doicu, Adrian, Dmitry Efremenko, and Thomas Trautmann. "A multi-dimensional
        vector spherical harmonics discrete ordinate method for atmospheric radiative transfer."
        Journal of Quantitative Spectroscopy and Radiative Transfer 118 (2013): 121-131.
        URL:http://dx.doi.org/10.1016/j.jqsrt.2012.12.009.
    """
    table_attempt = None
    refractive_index = None if particle_type == 'Water' else refractive_index
    if relative_dir is not None:
        table_attempt = _load_table(relative_dir, particle_type, wavelength_band,
                                    minimum_effective_radius, max_integration_radius,
                                    wavelength_averaging, wavelength_resolution,
                                    refractive_index)

    if table_attempt is not None:
        table = table_attempt
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
    It has python binding for the SHDOM fortran code contained within
    src/polarized/make_mie_table.f90. See 'get_mono_table' for more details.

    Parameters
    ----------
    particle_type: string
        Options are 'Water' or 'Aerosol'.
    wavelength_band: (float, float)
        (minimum, maximum) wavelength in microns.
        This defines the spectral band over which to integrate, if both are
        equal monochrome quantities are computed.
    minimum_effective_radius: float
        Minimum effective radius in microns. Used to compute minimum radius for integration.
    max_integration_radius: float
        Maximum radius in microns - cutoff for the size distribution integral
    wavelength_averaging: bool
        True - average scattering properties over the wavelength_band.
        False - scattering properties of the central wavelength.
    wavelength_resolution: float
        The distance between two wavelength samples in the band.
        Used only if wavelength_averaging is True.
    refractive_index: complex number or None
        For 'Water' the refractive index is ignored and loaded from a table
        in src/polarized/indexwatice. For 'Aerosol' the refractive index should
        have a negative imaginary part. ri = n - ik
        For water with costume refractive index use 'Aerosol' particle type.
    verbose: bool
        True for progress prints from the fortran computations.

    Returns
    -------
    table: xr.Dataset
        A Dataset containing tabulated Legendre coefficients per radius per wavelength.
        Each entry has 6 six Wigner d-function elements of the Mie phase matrix indexed by:
        P11, P22, P33, P44, P12, P34.

    Notes
    -----
    Currently supports only 'Water' and 'Aerosol' particles ('Ice' is not yet implemented).
    """
    # wavelength band
    wavelen1, wavelen2 = wavelength_band
    if wavelen1 > wavelen2:
        raise ValueError('wavelen1 must be <= wavelen2')

    avgflag = 'C'
    deltawave = -1
    if wavelength_averaging:
        avgflag = 'A'
        deltawave = wavelength_resolution

    wavelencen = at3d.core.get_center_wavelen(
        wavelen1=wavelen1,
        wavelen2=wavelen2
    )

    # set particle type properties
    if particle_type == 'Water':
        partype = 'W'
        refractive_index = at3d.core.get_refract_index(
            partype=particle_type,
            wavelen1=wavelen1,
            wavelen2=wavelen2
        )
        refractive_index_source = 'src/polarized/indexwatice.f'

    elif particle_type == 'Aerosol':
        partype = 'A'
        if refractive_index is None:
            raise AttributeError('A refractive index should be supplied for Aerosol particles')
        if refractive_index.imag > 0.0:
            raise AttributeError(
                'The refractive index should have a negative imaginary part. ri = n - ik')
        refractive_index_source = 'user_input'
    else:
        raise AttributeError('Particle type not implemented')

    # set integration parameters
    if avgflag == 'A':
        xmax = 2 * np.pi * max_integration_radius / wavelen1
    else:
        xmax = 2 * np.pi * max_integration_radius / wavelencen
    maxleg = int(np.round(2.0 * (xmax + 4.0 * xmax ** 0.3334 + 2.0)))

    # set radius integration parameters
    nsize = at3d.core.get_nsize(
        sretab=minimum_effective_radius,
        maxradius=max_integration_radius,
        wavelen=wavelencen
    )

    radii = at3d.core.get_sizes(
        sretab=minimum_effective_radius,
        maxradius=max_integration_radius,
        wavelen=wavelencen,
        nsize=nsize
    )
    #compute mie properties
    extinct, scatter, nleg, legcoef, ierr, errmsg = \
        at3d.core.compute_mie_all_sizes(
            nsize=nsize,
            maxleg=maxleg,
            wavelen1=wavelen1,
            wavelen2=wavelen2,
            deltawave=deltawave,
            wavelencen=wavelencen,
            radii=radii,
            rindex=refractive_index,
            avgflag=avgflag,
            partype=partype,
            verbose=verbose
        )
    at3d.checks.check_errcode(ierr, errmsg)
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
            'stokes_index': (['stokes_index'], ['P11', 'P22', 'P33', 'P44', 'P12', 'P34'])
            },
        attrs={
            'particle_type': particle_type,
            'refractive_index': (refractive_index.real, refractive_index.imag),
            'refractive_index_source': refractive_index_source,
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

def _load_table(relative_dir, particle_type, wavelength_band,
                minimum_effective_radius=4.0, max_integration_radius=65.0,
                wavelength_averaging=False, wavelength_resolution=0.001,
                refractive_index=None):
    """
    This methods tests whether there is an existing mie table within the given directory
    with the specified properties.

    Parameters
    ----------
    relative_dir: string
        The path to a directory which contains saved mie_table netcdf files. If there is a file with
        'mie_table' in the name that matches the input parameters this file is loaded.
    particle_type: string
        Options are 'Water' or 'Aerosol'.
    wavelength_band: (float, float)
        (minimum, maximum) wavelength in microns.
        This defines the spectral band over which to integrate, if both are equal
        monochrome quantities are computed.
    minimum_effective_radius: float
        Minimum effective radius in microns. Used to compute minimum radius for integration.
    max_integration_radius: float
        Maximum radius in microns - cutoff for the size distribution integral
    wavelength_averaging: bool
        True - average scattering properties over the wavelength_band.
        False - scattering properties of the central wavelength.
    wavelength_resolution: float
        The distance between two wavelength samples in the band.
        Used only if wavelength_averaging is True.
    refractive_index: complex number or None
        For 'Water' the refractive index is ignored and loaded from a table
        in src/polarized/indexwatice. For 'Aerosol' the refractive index should
        have a negative imaginary part. ri = n - ik
        For water with costume refractive index use 'Aerosol' particle type.

    Returns
    -------
    table: xr.Dataset or None
        If the table exists within the specified path it is loaded and returned.
        Alternatively, returns None.
    """
    table = None
    if os.path.exists(relative_dir):
        file_list = os.listdir(relative_dir)

        list_of_input = [particle_type, wavelength_band,
                         minimum_effective_radius, max_integration_radius,
                         str(wavelength_averaging), wavelength_resolution]
        names = ['particle_type', 'wavelength_band',
                 'minimum_effective_radius', 'maximum_integration_radius',
                 'wavelength_averaging', 'wavelength_resolution']

        if refractive_index is not None:
            list_of_input.append((refractive_index.real, refractive_index.imag))
            names.append('refractive_index')

        for file in file_list:
            try:
                if file.endswith('.nc'):
                    with xr.open_dataset(os.path.join(relative_dir, file)) as dataset:
                        # open_dataset reads data lazily.
                        if ('extinction' not in dataset.data_vars) or \
                            ('scatter' not in dataset.data_vars) or \
                            ('nleg' not in dataset.data_vars) or \
                            ('legendre' not in dataset.data_vars) or \
                            ('radius' not in dataset.coords) or \
                            ('stokes_index' not in dataset.coords):
                            continue
                        test = True
                        for name, attr in zip(names, list_of_input):
                            try:
                                attr_file = dataset.attrs[name]
                            except KeyError:
                                test = False
                                break

                            if isinstance(attr, tuple):
                                temp = np.allclose(attr, attr_file)
                            else:
                                temp = attr == attr_file

                            if not temp:
                                test = False
                                break
                        if test:
                            table = file
            except IOError:
                continue
        if table is not None:
            table = xr.load_dataset(os.path.join(relative_dir, table))
    return table

def get_poly_table(size_distribution, mie_mono_table):
    """
    This methods calculates Mie scattering table for a polydisperse size distribution.
    For more information about the size_distribution see: lib/size_distribution.py.

    This function integrates the Mie table over radii for each entry in the
    size distribution Dataset. This dataset could be parameterized arbitrarily,
    however, a common use-case is according to effective radius and variace.

    Parameters
    ----------
    size_distribution: xr.Dataset
        A Dataset of number_density variable as a function of radius and the table parameterization.
        The shape of the dataset is ('radius', 'param1', 'param2',...'paramN').
        A common case is ('radius', 'reff', 'veff').
    mie_mono_table: xr.Dataset
        A Dataset of Mie legendre coefficients as a function of radius.
        See mie.get_mono_table function for more details.

    Returns
    -------
    poly_table: xr.Dataset
        A Dataset with the polydisperse Mie effective scattering properties:
        extinction, ssalb, legcoef. Each of these are a function of the parameterization
        defined by the size_distribution.

    Raises
    ------
    AssertionError
        If the mie_mono_table are not within the range of the size_distribution radii.

    Notes
    -----
    The radius in size_distribution is interpolated onto the mie_mono_table radii grid.
    This is to avoid interpolation of the Mie table coefficients.
    """
    at3d.checks.check_range(
        mie_mono_table,
        radius=(size_distribution.radius.min(), size_distribution.radius.max())
        )

    if (size_distribution.radius.size != mie_mono_table.radius.size) or \
            np.any(size_distribution.radius.data != mie_mono_table.radius.data):
        print('Warning: size_distribution radii differ to mie_mono_table radii. '
              'Interpolating the size distribution onto the Mie table grid.')
        size_distribution = size_distribution.interp(radius=mie_mono_table.radius)

    number_density = size_distribution['number_density'].values.reshape(
        (len(size_distribution['radius'])), -1
        )

    extinct, ssalb, nleg, legcoef = \
        at3d.core.get_poly_table(
            nd=number_density,
            ndist=number_density.shape[-1],
            nsize=mie_mono_table.coords['radius'].size,
            maxleg=mie_mono_table.attrs['maximum_legendre'],
            nleg1=mie_mono_table['nleg'],
            extinct1=mie_mono_table['extinction'],
            scatter1=mie_mono_table['scatter'],
            legcoef1=mie_mono_table['legendre'])

    grid_shape = size_distribution['number_density'].shape[1:]

    # all coords except radius
    coords = {name:coord for name, coord in size_distribution.coords.items()
              if name not in ('radius', 'stokes_index')}
    microphysics_names = list(coords.keys())
    coord_lengths = [np.arange(coord.size) for name, coord in coords.items()]
    legen_index = np.meshgrid(*coord_lengths, indexing='ij')

    table_index = np.ravel_multi_index(legen_index, dims=[coord.size for coord in coord_lengths])
    coords['table_index'] = (microphysics_names, table_index)
    coords['stokes_index'] = mie_mono_table.coords['stokes_index']

    poly_table = xr.Dataset(
        data_vars={
            'extinction': (microphysics_names, extinct.reshape(grid_shape)),
            'ssalb': (microphysics_names, ssalb.reshape(grid_shape)),
            'legcoef': (['stokes_index', 'legendre_index'] + microphysics_names,
                        legcoef.reshape(legcoef.shape[:2] + grid_shape)),},
        coords=coords
    )
    poly_table = poly_table.assign_attrs(size_distribution.attrs)
    poly_table = poly_table.assign_attrs(mie_mono_table.attrs)
    return poly_table
