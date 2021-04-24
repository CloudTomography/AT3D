"""
This module contains functions to map microphysical variables to
optical properties. In the original SHDOM, these procedures would be performed
by the `propgen` program.
Currently only a Look-Up-Table approach is supported defined by `table_to_grid`
which uses linear interpolation and nearest interpolation for phase functions.
This method is also used to map partial derivatives of optical properties
from the Look-Up-Table onto the SHDOM grid.
"""
import itertools
from collections import OrderedDict
import pandas as pd
import numpy as np
import xarray as xr
import pyshdom.checks

def table_to_grid(microphysics, poly_table, inverse_mode=False):
    """
    Calculates optical properties from microphysical properties using a Look-Up-Table.

    Optical properties are calculated on the 3D grid defined in `microphysics`
    using microphysical parameters defined in `microphysics` and a Look-Up-Table to
    map between microphysical variables and optical variables in `poly_table`
    using linear interpolation.

    Parameters
    ----------
    microphysics : xr.Dataset
        Should contain the microphysical parameters on a valid grid.
        See grid.py for details. The microphysical parameters should match
        those in `poly_table`.
    poly_table : xr.Dataset
        Should contain extinction efficiency, single scatter albedo and phase
        functions as a function microphysical parameters (e.g. effective radius).
    exact_table : bool
        Sets the interpolation method for the calculation of extinction and ssalb.
        linear if False, and nearest if True.
    inverse_mode : bool
        A flag for whether optical properties or their derivatives with respect
        to microphysical properties are being interpolated.
        The only difference is that extinction_efficiency passed instead of
        extinction.

    Returns
    -------
    optical_properties : xr.Dataset
        A dataset containing optical properties defined on an SHDOM grid
        ready to be used as input to solver.RTE.

    Notes
    -----
    The linear interpolation of extinction and single scatter albedo and
    nearest neighbor interpolation of phase functions used here differs
    from SHDOM's propgen program, which uses a different interpolation scheme
    and creates new phase functions if a specified accuracy is not met.
    This is not implemented as the Look-Up-Table approach implemented here
    is simple to apply to the inverse problem without the need to recompute
    any mie properties. To ensure high accuracy, the spacing of the table
    should be fine, see size_distribution.py.
    """
    pyshdom.checks.check_positivity(microphysics, 'density')
    pyshdom.checks.check_grid(microphysics)
    if not inverse_mode:
        pyshdom.checks.check_legendre(poly_table)

    interp_names = set([name for name in poly_table.coords
                        if name not in ('table_index', 'stokes_index')])
    microphysics_names = set([name for name in microphysics.variables.keys()
                              if name not in 'density'])
    missing = interp_names - microphysics_names
    if missing:
        raise KeyError(
            "microphysics dataset is missing variables "
            "for interpolation of table onto grid.", *list(missing)
            )
    interp_coords = {name:microphysics[name] for name in poly_table.coords
                     if name not in ('table_index', 'stokes_index')}
    for interp_coord in interp_coords:
        if np.any(microphysics[interp_coord] <= poly_table[interp_coord].min()) or \
            np.any(microphysics[interp_coord] >= poly_table[interp_coord].max()):
            raise ValueError(
                "Microphysical coordinate '{}' is not"
                " within the range of the mie table.".format(interp_coord)
                )
    interp_method = 'linear' #this is hardcoded as it is suitable for N-dimensions.
    ssalb = poly_table.ssalb.interp(interp_coords, method=interp_method)
    extinction_efficiency = poly_table.extinction.interp(interp_coords, method=interp_method)

    if not inverse_mode:
        extinction = extinction_efficiency * microphysics.density
    else:
        extinction = extinction_efficiency

    assert not np.any(np.isnan(ssalb.data)), 'Unexpected NaN in ssalb'
    assert not np.any(np.isnan(extinction.data)), 'Unexpected NaN in extinction'

    # Different method for phase functions / indices.
    extinction.name = 'extinction'
    table_index = poly_table.coords['table_index'].interp(
        coords=interp_coords, method='nearest'
        ).round().astype(int)
    unique_table_indices, inverse = np.unique(table_index.data, return_inverse=True)
    subset_table_index = xr.DataArray(name=table_index.name,
                                      data=inverse.reshape(table_index.shape) + 1,
                                      dims=table_index.dims,
                                      coords=table_index.coords,
                                      )

    legendre_table_stack = poly_table['legcoef'].stack(table_index=interp_coords)
    subset_legcoef = legendre_table_stack.isel({'table_index':unique_table_indices})
    subset_legcoef = xr.DataArray(name='legcoef',
                                  data=subset_legcoef.data,
                                  dims=['stokes_index', 'legendre_index', 'table_index'],
                                  coords={'stokes_index':subset_legcoef.coords['stokes_index']}
                                 )

    optical_properties = xr.merge([extinction, ssalb, subset_legcoef])
    optical_properties['density'] = microphysics.density
    table_coords = {'table_index': (['num_micro', 'x', 'y', 'z'], subset_table_index.data[np.newaxis])}

    optical_properties = optical_properties.assign_coords(table_coords)
    optical_properties['phase_weights'] = (['num_micro', 'x', 'y', 'z'], np.ones(optical_properties.table_index.shape))
    assert not np.any(np.isnan(optical_properties.table_index.data)), 'Unexpected NaN in table_index'
    assert not np.any(np.isnan(optical_properties.legcoef.data)), 'Unexpected NaN in legcoef'

    #inherit attributes.
    optical_properties = optical_properties.assign_attrs(poly_table.attrs)
    optical_properties = optical_properties.assign_attrs(microphysics.attrs)

    #transfer the grid variables. NOTE that delx, dely exist and be passed.
    # while nx/ny/nz are optional. delx/dely are checked for by check_grid.
    grid_variables = ('delx', 'dely', 'nx', 'ny', 'nz')
    for grid_variable in grid_variables:
        if grid_variable in microphysics.data_vars:
            optical_properties[grid_variable] = microphysics[grid_variable]

    return optical_properties

def get_optical_properties(microphysics, mie_mono_tables, size_distribution_function,
                           size_distribution_grid_parameters, particle_density=1.0,
                           maxnphase=None):
    """
    Calculates optical properties from microphysical properties either exactly
    or by specifying linear mixtures that are combined 'just-in-time' during the
    RTE solution.

    Optical properties are calculated on the 3D grid defined in `microphysics`
    using microphysical parameters defined in `microphysics`. Each grid point
    has a pointer to a single or multiple phase functions and their corresponding
    weights with a preference for a single phase function unless maxnphase is
    exceeded.

    Parameters
    ----------
    microphysics : xr.Dataset
        Should contain the microphysical parameters on a valid grid.
        See grid.py for details. The microphysical parameters should match
        those in specified in size_distribution_grid_parameters.
    mie_mono_tables: Dict
        A dictionary of Datasets of Mie legendre coefficients as a function of radius.
        See mie.get_mono_table function for more details.
    size_distribution_function: callable
        Predefined function that computes size distribution.
        Implemented options here are `gamma` (default) and `lognormal`.
    size_distribution_grid_parameters: dict
        Sets the spacing (and hence accuracy) of the microphysical parameters
        for the computation of optical properties.
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
    maxnphase : int
        Sets the maximum number of unique phase functions that can be used.
        Default is None, in which case the maximum size is the total number
        of points specified in the size_distribution_grid.

    Returns
    -------
    optical_properties : OrderedDict
        A dictionary of datasets containing optical properties defined on an SHDOM grid
        ready to be used as input to solver.RTE.

    Notes
    -----
    This scheme was implemented to ensure continuous variation of
    phase functions with changing microphysics which `table_to_grid` does not
    support. It is more similar to SHDOM's Propgen program but does not limit
    phase functions based on sufficient accuracy. If too many phase functions
    are required it will switch to just-in-time mixing of phase functions
    which will increase computation time but will always be at least as
    accurate as Propgen for a given spacing of microphysical parameters
    (unlike `table_to_grid`) and possibly more so at the expense of
    computation time.
    """
    pyshdom.checks.check_positivity(microphysics, 'density')
    pyshdom.checks.check_grid(microphysics)

    interp_names = set([name for name in size_distribution_grid_parameters])
    microphysics_names = set([name for name in microphysics.variables.keys()
                              if name not in 'density'])
    missing = interp_names - microphysics_names
    if missing:
        raise KeyError(
            "microphysics dataset is missing variables "
            "for interpolation of table onto grid.", *list(missing)
            )

    # make the grid of size distributions. This is only utilized for the computation of the
    # coordinates so this could be sped up by isolating that code instead.
    number_density_grid = pyshdom.size_distribution.get_size_distribution_grid(
        list(mie_mono_tables.values())[0].radius,
        size_distribution_function=size_distribution_function,
        particle_density=particle_density,
        radius_units='micron',
        **size_distribution_grid_parameters
    )
    # if no max number of phase functions is specified then we set it to the total number of table
    # entries specified.
    if maxnphase is None:
        maxnphase = number_density_grid.number_density.size / number_density_grid.radius.size

    interp_coords = {name:microphysics[name] for name in number_density_grid.coords
                     if name not in ('radius',)}

    for interp_coord in interp_coords:
        if np.any(microphysics[interp_coord] <= number_density_grid[interp_coord].min()) or \
            np.any(microphysics[interp_coord] >= number_density_grid[interp_coord].max()):
            raise ValueError(
                "Microphysical coordinate '{}' is not"
                " within the range of the size distribution parameters.".format(interp_coord)
                )

    # use arrays rather than dictionaries to make use of np.unique
    microphysics_data = np.stack([data.data.ravel() for data in interp_coords.values()], axis=0)

    # digitize to find size distribution coordinate bin indices of each variable.
    digitized_microphysics = np.stack(
        [np.digitize(micro.data.ravel(), bins=number_density_grid[name].data)
         for name, micro in interp_coords.items()], axis=0
        )

    # rule for forming phase functions:
    #   1. maximize exact phase functions up to MAXNPHASE
    #   2. if maxnphase exceeded then force to table mode.
    #       treat table vs exact on a variable by variable basis
    #       rather than globally or case by case.

    # First check if any table is needed.
    unique_values = np.unique(microphysics_data, axis=1)
    num_unique = unique_values.shape[1]
    parameter_dict = {}
    if num_unique < maxnphase:
        # make the set of phase indices and interpolation weights.
        # when only exact phase functinos are used.
        # and specify the exact microphysical variables to calculate optical properties for.
        for i, (name, variable_data) in enumerate(zip(interp_coords, microphysics_data)):
            parameter_dict[name] = variable_data
        phase_indices = np.arange(microphysics_data.shape[1])[np.newaxis, :]
        interpolation_weights = np.ones(phase_indices.shape)
    else:
        # lower_upper_flags contains -1,0 for finding the indices of lower and upper bounds
        # based on the bin indices (given the default behaviour of np.digitize).
        lower_upper_flags = np.array(
            [i for i in itertools.product(*[np.arange(-1, 1)]*(microphysics_data.shape[0]))]
            )
        lower_upper_combinations = digitized_microphysics + lower_upper_flags[..., np.newaxis]

        # find the exact microphysical values that we would use if we don't convert to table.
        combinations = []
        exact_points = []
        for variable_data in microphysics_data:
            unique_points, inverses = np.unique(variable_data, return_inverse=True)
            combinations.append(inverses)
            exact_points.append(unique_points)
        combinations = np.repeat(np.stack(combinations, axis=0)[np.newaxis, ...],
                                 lower_upper_flags.shape[0], axis=0)

        # select which variables to first turn into table to most quickly reduce
        # the number of phase functions.
        unique_combinations = np.unique(lower_upper_combinations, axis=2) #number of unique bins.

        # `test` is the metric for choosing which variable to turn into table.
        #  This is not globally optimal.
        test = np.array(
            [np.unique(unique_combinations[:, i, :]).size/np.unique(microphysics_data[i]).size
             for i in range(unique_combinations.shape[1])]
            )

        # convert variables to table until maxnphase is no longer exceeded.
        variables_to_turn_to_table = np.argsort(test)
        table_counter = 0
        while num_unique > maxnphase:
            index = variables_to_turn_to_table[table_counter]
            combinations[:, index, :] = lower_upper_combinations[:, index, :]
            unique = np.unique(
                combinations.transpose([0, -1, 1]).reshape(
                    combinations.shape[0]*combinations.shape[-1], -1),
                axis=0
            )
            num_unique = unique.shape[0]
            table_counter += 1
            if (table_counter > microphysics_data.shape[0]) & (num_unique > maxnphase):
                raise ValueError(
                    "All variables are represented using tables but `maxnphase`='{}'' is "
                    "still exceeded num_unique='{}'. Please increase `maxnphase` or set to "
                    "`None` so that all of the table may be used.".format(maxnphase, num_unique)
                )

        # redo the table with the smaller number of required bounding variables.
        # This is could be more efficiently.
        # taking a subset of `combinations` but would be a little complicated when there
        # are more than two table variables so I just redo.
        # Also this re-computes the table when all variables are table even though that is now
        # in `combinations`. There is a shape difference.
        lower_upper_flags = np.array(
            [i for i in itertools.product(*[np.arange(-1, 1)]*table_counter)]
            )
        table_combinations = digitized_microphysics[variables_to_turn_to_table[:table_counter]] + \
            lower_upper_flags[..., np.newaxis]
        table_unique, table_inverse = np.unique(
            table_combinations.transpose([0, -1, 1]).reshape(
                table_combinations.shape[0]*table_combinations.shape[-1], -1),
            axis=0, return_inverse=True
            )

        # make the final set of phase indices
        phase_indices = table_inverse.reshape(
            table_combinations.shape[0], table_combinations.shape[-1], -1
            ).transpose([0, -1, 1])
        phase_indices = phase_indices[:, 0, :]

        # make the interpolation weights.
        interpolation_weights = np.ones(table_combinations.shape)
        for i in range(table_counter):
            index = variables_to_turn_to_table[i]
            name = list(interp_coords)[index]
            # make the values to generate while we are at it.
            parameter_dict[name] = number_density_grid[name].data[table_unique[:, i]]

            lower = number_density_grid[name].data[table_unique[phase_indices][0, :, i]]
            upper = number_density_grid[name].data[table_unique[phase_indices][-1, :, i]]
            step = upper - lower
            interp_weight = (microphysics_data[index] - lower)/step

            condition = np.where(lower_upper_flags[:, i] == -1)
            condition2 = np.where(lower_upper_flags[:, i] == 0)
            interpolation_weights[condition, i] = (1.0 - interp_weight)
            interpolation_weights[condition2, i] = interp_weight

        interpolation_weights = np.prod(interpolation_weights, axis=1)

        # add all the exact optical properties to generate.
        for i in range(table_counter, microphysics_data.shape[0]):
            index = variables_to_turn_to_table[i]
            name = list(interp_coords)[index]
            parameter_dict[name] = exact_points[index][unique[:, index]]

    assert np.allclose(interpolation_weights.sum(axis=0), 1.0)
    optical_properties = OrderedDict()

    for key, mie_mono_table in mie_mono_tables.items():

        number_density_raveled = size_distribution_function(
            mie_mono_table.radius,
            **parameter_dict,
            particle_density=particle_density
        )
        size_dist_attrs = OrderedDict()
        size_dist_attrs['distribution_type'] = size_distribution_function.__name__
        size_dist_attrs['radius_units'] = number_density_grid.radius_units
        coords = {'radius': mie_mono_table.radius.data}
        coords['microphysics_index'] = pd.MultiIndex.from_arrays(
            [parameter_dict[name] for name in parameter_dict],
            names=list(parameter_dict)
            )

        size_dist_grid = xr.Dataset(
            data_vars={'number_density': (['radius', 'microphysics_index'], number_density_raveled)},
            coords=coords,
            attrs=size_dist_attrs
        )
        poly_table = pyshdom.mie.get_poly_table(size_dist_grid, mie_mono_table)
        # make sure this worked.
        pyshdom.checks.check_legendre(poly_table)

        extinct_efficiency = np.sum(
            interpolation_weights*poly_table.extinction.data[phase_indices],
            axis=0
            )
        ssalb = np.sum(interpolation_weights*poly_table.ssalb.data[phase_indices], axis=0)

        extinction = microphysics.density*extinct_efficiency.reshape(microphysics.density.shape)

        assert not np.any(np.isnan(ssalb)), 'Unexpected NaN in ssalb'
        assert not np.any(np.isnan(extinction)), 'Unexpected NaN in extinction'

        poly_table['extinction'] = extinction
        poly_table['ssalb'] = (['x', 'y', 'z'], ssalb.reshape(microphysics.density.shape))

        poly_table['table_index'] = (
            ['num_micro', 'x', 'y', 'z'],
            1+phase_indices.reshape([-1] + list(microphysics.density.shape))
        )
        poly_table['phase_weights'] = (
            ['num_micro', 'x', 'y', 'z'],
            interpolation_weights.reshape([-1] + list(microphysics.density.shape))
        )

        #transfer the grid variables. NOTE that delx, dely exist and be passed.
        # while nx/ny/nz are optional. delx/dely are checked for by check_grid.
        grid_variables = ('delx', 'dely', 'nx', 'ny', 'nz')
        for grid_variable in grid_variables:
            if grid_variable in microphysics.data_vars:
                poly_table[grid_variable] = microphysics[grid_variable]

        poly_table = poly_table.assign_attrs(microphysics.attrs)
        optical_properties[key] = poly_table

    return optical_properties
