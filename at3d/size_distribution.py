"""
This module contains functions to define particle size distributions
and to sample these distributions at radii for a grid/Look-up-Table of
size distribution parameters.

In at3d, the mapping of microphysical properties to optical properties is
done using a look-up-table approach. Bulk optical properties used in the
radiative transfer solver are determined by integrating over single scattering
properties (calculated from mie.py or other methods) using assumed distributions
of particle properties.
The parameterization of the size distributions used, and the spacing and range of the
distribution (microphysical) parameters are specified here.
Two size distributions are currently supported (`gamma` and `lognormal`), but
`get_size_distribution_grid` can produce look-up-tables of size distributions
based on arbitrary distribution parameterizations (e.g. bins, or multi-modes),
all that is required is to add a callable similar to `gamma`/`lognormal` that
generates the weight of the distribution at a set range of radii given the
input parameters.
"""
import typing
from collections import OrderedDict
import xarray as xr
import numpy as np
import at3d.core

def gamma(radii, reff, veff=None, alpha=None, particle_density=1.0,
          normalization='density'):
    """
    Generate a Gamma size distribution.
    Provide either effective variance `veff` or shape parameter `alpha`.

    Parameters
    ----------
    radii: scalar, list/numpy array
        Radii of precomputed Mie tables in `radius units`.
    reff: scalar, list/numpy array
        Effective radius for which to compute the size distribution.
    veff: scalar, list/numpy array
        Effective variance of the size distribution in `radius units`.
        The larger this parameter, the broader the size distribution.
    alpha: scalar, list/numpy array
        Shape parameter.
    particle_density: float
        Particle density in [g/m^3]. Default 1 g/m^3 for Water.
    normalization: str
        Choice of size distribution normalization from
        `density` [g/m^3], `geometric_extinction` [1/km]
        (extinction calculated using a fixed extinction efficiency of 2),
        `number_concentration` [#/cm^3]

    Returns
    -------
    number_density: ndarray
        Number density of the shape (len(radii), len(reff)).

    Notes
    -----
    Given effective variance `veff`, the shape parameter `alpha` is computed according to:
    alpha = 1.0 / veff - 3.0

    """
    reff = np.atleast_1d(reff)
    if veff is not None:
        veff = np.atleast_1d(veff)
        if len(reff) != len(veff):
            raise ValueError("'reff' and 'veff' must have the same length")
        alpha = 1.0/ veff - 3.0
    elif alpha is not None:
        alpha = np.atleast_1d(alpha)
        if len(reff) != len(alpha):
            raise ValueError("'reff' and 'alpha' must have the same length")
    else:
        raise ValueError("One of 'veff' of 'alpha' must be specified.")

    #fortran subroutine requires a 'gamma' variable to be passed.
    gamma_values = np.zeros(alpha.shape)

    number_density, ierr, errmsg = at3d.core.make_multi_size_dist(
                distflag='G',
                pardens=particle_density,
                nsize=len(radii),
                radii=radii,
                reff=reff,
                alpha=alpha,
                gamma=gamma_values,
                ndist=reff.size)
    at3d.checks.check_errcode(ierr, errmsg)

    if normalization == 'geometric_extinction':
        number_density /= 1e-3*(2*number_density*np.pi*np.atleast_1d(radii)[:, None]**2).sum(axis=0)
    elif normalization == 'number_concentration':
        number_density /= (number_density).sum(axis=0)
    elif normalization == 'density':
        # explicit default case.
        pass
    else:
        raise ValueError(
            "Invalid `normalization` argument '{}'".format(normalization)
        )

    return number_density


def lognormal(radii, reff, veff=None, alpha=None, particle_density=1.0,
              normalization='density'):
    """
    Generate a Log-normal size distribution.
    Provide either effective variance `veff` or shape parameter `alpha`.

    Parameters
    ----------
    radii: scalar, list/numpy array
        Radii of precomputed Mie tables.
    reff: scalar, list/numpy array
        Effective radius for which to compute the size distribution.
    veff: scalar, list/numpy array
        Effective variance of the size distribution.
        The larger this parameter, the broader the size distribution.
    alpha: scalar, list/numpy array
        Log-normal standard deviation.
    particle_density: float
        Particle density in [g/m^3]. Default 1 g/m^3 for Water.
    normalization: str
        Choice of size distribution normalization from
        `density` [g/m^3], `geometric_extinction` [1/km]
        (extinction calculated using a fixed extinction efficiency of 2),
        `number_concentration` [#/cm^3]

    Returns
    -------
    number_density: ndarray
        Number density of the shape (len(radii), len(reff)).

    Notes
    -----
    Given effective variance `veff`, the shape parameter `alpha` is computed according to:
    alpha = sqrt(log(veff + 1)

    """
    reff = np.atleast_1d(reff)
    if veff is not None:
        veff = np.atleast_1d(veff)
        if len(reff) != len(veff):
            raise ValueError("'reff' and 'veff' must have the same length")
        alpha = np.sqrt(np.log(np.atleast_1d(veff) + 1.0))
    elif alpha is not None:
        alpha = np.atleast_1d(alpha)
        if len(reff) != len(alpha):
            raise ValueError("'reff' and 'alpha' must have the same length")
    else:
        raise ValueError("One of 'veff' of 'alpha' must be specified.")

    #fortran subroutine requires a 'gamma' attribute to be passed.
    gamma_values = np.zeros(alpha.shape)

    number_density, ierr, errmsg = at3d.core.make_multi_size_dist(
                distflag='L',
                pardens=particle_density,
                nsize=len(radii),
                radii=radii,
                reff=reff,
                alpha=alpha,
                gamma=gamma_values,
                ndist=reff.size)
    at3d.checks.check_errcode(ierr, errmsg)

    if normalization == 'geometric_extinction':
        number_density /= 1e-3*(2*number_density*np.pi*np.atleast_1d(radii)[:, None]**2).sum(axis=0)
    elif normalization == 'number_concentration':
        number_density /= (number_density).sum(axis=0)
    elif normalization == 'density':
        # explicit default case.
        pass
    else:
        raise ValueError(
            "Invalid `normalization` argument '{}'".format(normalization)
        )

    return number_density


def get_size_distribution_grid(radii, size_distribution_function=gamma,
                               particle_density=1.0, radius_units='micron',
                               **size_distribution_parameters):
    """
    A generic interface to get number density (in cm^-3) of a size distribution
    on a grid of size distribution parameters.

    Parameters
    ----------
    radii: list/numpy array
        Radii of precomputed Mie tables in `radius units`
    size_distribution_function: callable
        Predefined function that computes size distribution.
        Implemented options here are `gamma` (default) and `lognormal`.
    particle_density: float
        Particle density in [g/m^3]. Default 1 g/m^3 for Water.
    radius_units: string
        Unit of radii, default [microns].
    size_distribution_parameters: dict
        Each size_distribution_parameter dictionary should contain the following keys:
           'coord_min': float
               The minimum value for that coordinate.
           'coord_max': float
               The maximum value for that coordinate.
           'npoints': integer
               The number of points sampling this dimension.
           'spacing': string
               The type of spacing. Either 'logarithmic' or 'linear'.
           'coord': 1D array
                This overrides the above arguments and specifies the exact
                points to sample along this dimension.
           'units': string
               The units of the microphysical dimension.
        Alternatively, if a 1D numpy array is specified it will be interpreted as
        the 'coord' argument.
    Returns
    -------
    size_dist_grid: xarray.Dataset
        Dataset containing `number_density` in [cm^-3].

    Examples
    --------
    >>> import numpy as np
    >>> radii = np.arange(100)
    >>> size_dist_grid = at3d.size_distribution.get_size_distribution_grid(
                                    radii,
                                    size_distribution_function=at3d.size_distribution.gamma,
                                    particle_density=1.0,
                                    radius_units='micron',
                                    reff={'coord_min':4.0,'coord_max':25.0,
                                          'npoints':25,'spacing':'logarithmic',
                                          'units':'micron'},
                                    veff={'coord_min':0.09,'coord_max':0.11,
                                          'npoints':2,'spacing':'linear',
                                          'units':'unitless'},
                                    )
    >>> size_dist_grid
        Dimensions:         (radius: 100, reff: 25, veff: 2)
        Coordinates:
          * radius          (radius) int64 0 1 2 3 4 5 6 7 8 ... 92 93 94 95 96 97 98 99
          * reff            (reff) float64 4.0 4.317 4.66 5.03 ... 21.46 23.16 25.0
          * veff            (veff) float64 0.09 0.11
        Data variables:
            number_density  (radius, reff, veff) float32 0.0 0.0 ... 1.483e-10 3.284e-09
        Attributes:
            reff_coord_min:     4.0
            reff_coord_max:     25.0
            reff_npoints:       25
            reff_spacing:       logarithmic
            reff_units:         micron
            veff_coord_min:     0.09
            veff_coord_max:     0.11
            veff_npoints:       2
            veff_spacing:       linear
            veff_units:         unitless
            distribution_type:  gamma
            radius_units:       radius units [micron]

    """

    coord_list = []
    name_list = []
    for name, parameter in size_distribution_parameters.items():
        if isinstance(parameter, np.ndarray):
            if parameter.ndim == 1:
                coord_list.append(parameter)
        elif 'coords' in parameter:
            coord_list.append(parameter['coords'])
        else:
            if parameter['spacing'] == 'logarithmic':
                coord = np.logspace(np.log10(parameter['coord_min']),
                                    np.log10(parameter['coord_max']),
                                    parameter['npoints'])
            elif parameter['spacing'] == 'linear':
                coord = np.linspace(parameter['coord_min'],
                                    parameter['coord_max'],
                                    parameter['npoints'])
            else:
                raise NotImplementedError
            coord_list.append(coord)
        name_list.append(name)

    meshgrid = np.meshgrid(*coord_list, indexing='ij')
    parameter_dict = OrderedDict()
    coord_dict = OrderedDict()
    coord_dict['radius'] = radii
    for name, grid, coord in zip(name_list, meshgrid, coord_list):
        parameter_dict[name] = grid.ravel()
        coord_dict[name] = coord

    grid_shape = [len(coord) for name, coord in coord_dict.items()]

    number_density_raveled = size_distribution_function(
        radii,
        **parameter_dict,
        particle_density=particle_density
    )
    if np.any(np.isnan(number_density_raveled)) or np.any(number_density_raveled < 0.0):
        raise ValueError(
            "size_distribution_function: {} produced NaN or negative values "
            "when being evaluated on grid of parameters."
            )

    number_density = number_density_raveled.reshape(grid_shape)

    # create "flat" attrs dictionary to enable saving to netCDF
    size_dist_attrs = OrderedDict()
    for name, parameter in size_distribution_parameters.items():
        if isinstance(parameter, typing.Dict):
            for _name, _param in parameter.items():
                size_dist_attrs[f"{name}_{_name}"] = _param
    size_dist_attrs['distribution_type'] = size_distribution_function.__name__
    size_dist_attrs['radius_units'] = 'radius units [{}]'.format(radius_units)

    size_dist_grid = xr.Dataset(
        data_vars={'number_density': (list(coord_dict.keys()), number_density)},
        coords=coord_dict,
        attrs=size_dist_attrs
    )
    return size_dist_grid

# def get_size_distribution_exact(radii, microphysics, size_distribution_function=gamma,
#                                 particle_density=1.0, radius_units='micron'):
#     """TODO"""
#     names = {name for name in microphysics.data_vars if name != 'density'}
#
#     unique_values = [np.unique(microphysics[name]) for name in names]
#     all_combinations = np.meshgrid(*unique_values, indexing='ij')
#     parameter_dict = {name:variable.ravel() for name, variable in zip(names, all_combinations)}
#
#     number_density_raveled = size_distribution_function(radii, **parameter_dict,
#                                                         particle_density=particle_density)
#     #TODO this can fail silently in some cases (produce nans), add checks.
#     number_density = number_density_raveled.reshape(radii.shape+all_combinations[0].shape)
#
#     coord_dict = OrderedDict()
#     coord_dict['radius'] = radii
#     for name, unique in zip(names, unique_values):
#         coord_dict[name] = unique
#
#     size_dist_grid = xr.Dataset(
#             data_vars = {
#                 'number_density': (list(coord_dict.keys()), number_density),
#             },
#             coords=dict(coord_dict),
#             attrs={
#                   'distribution_type': size_distribution_function.__name__,
#                   },
#                   #TODO add coordmin/max/n and spacing for full traceability.
#     )
#     return size_dist_grid
