"""
This module contains functions to map microphysical variables to
optical properties and to update solvers with new optical properties
in the inverse problem.

In the original SHDOM, these procedures would be performed
by the `PROPGEN` program.
`mix_optical_properties` performs a similar function to `PROPGEN` but it is
not appropriate for the inverse problem as optical properties do not vary
smoothly with changes in microphysics.

OpticalPropertyGenerator is the improved object which generates exact
optical properties or applies linear interpolation from a table. It also
evaluates the partial derivatives of optical properties with respect to
microphysical unknowns.

DataGenerator objects are used in the inverse problem to hold fixed variables
for each scatterer to produce updated optical properties as the state changes.

StateGenerator is used in the inverse problem to update the solvers as the
state changes (see at3d.containers.SolversDict).
"""
import itertools
import typing
import warnings
import time

from collections import OrderedDict
import pandas as pd
import numpy as np
import xarray as xr
import at3d.checks

def make_optical_properties(grid, extinction, ssalb, phase_function_indices,
                            phase_function_table, phase_function_weights=None):
    """
    Construct an xr.Dataset optical properties by specifying extinction, single scatter albedo
    and gridpoint phase functions.

    This method is provided for the direct specification of optical properties
    rather than using a Mie or other microphysical model which is useful when constructing
    idealized media. This leaves it to the useful to specify the `phase_function_weights`
    and `phase_function_indices` and therefore is best utilized in only simple scenarios.

    Parameters
    ----------
    grid : xr.Dataset
        A valid grid object with x, y, z coordinates. See at3d.grid for details.
    extinction : np.ndarray, float, ndim=3
        An array of volume extinction coefficient values in [1/km].
        This should have the same shape as the `grid` object.
    ssalb : np.ndarray, float ndim=3
        An array of single scatter albedos in [0, 1]. This should have the same shape as the `grid` object.
    phase_function_indices : np.ndarray, int, ndim=4
        An array of indices which point to entries in `phase_function_table` to indicate which entries are used
        to form each grid point function. The first dimension denotes the set of phase functions that are
        combined according to `phase_function_weights` to form the phase function at each grid point. The
        next 3 dimensions should match the shape of `grid`.
        Indices are assumed to start from 0 though in the output they start from 1 according to
        Fortran conventions.
    phase_function_table : np.ndarray, ndim=3
        The legendre / Wigner-d function expansion of the phase functions of
        shape=(StokesIndex, legendreIndex, PhaseFunctionIndex)
    phase_function_weights : np.ndarray, float, ndim=4
        Weights that are used to represent the grid point phase functions as a sum over the table entries.
        If `phase_function_indices` is only specified with one pointer per grid point and this is None
        then uniform weights are specified.

    Return : optical_properties, xr.Dataset
        Valid optical properties
    """
    at3d.checks.check_grid(grid)
    if extinction.shape != grid.grid.shape:
        raise ValueError(
            "`extinction` shape is not consistent with `grid`."
        )
    if ssalb.shape != grid.grid.shape:
        raise ValueError(
            "`extinction` shape is not consistent with `grid`."
        )
    if phase_function_indices.ndim == 3:
        phase_function_indices = phase_function_indices[None, ...]

    if phase_function_indices.shape[0] == 1:
        phase_function_weights = np.ones(phase_function_indices.shape)
    elif phase_function_weights is None:
        raise ValueError(
            "phase_function_weights must be specified."
        )

    if phase_function_indices.shape[1:] != grid.grid.shape:
        raise ValueError(
            "`phase_function_indices` shape is not consistent with `grid`."
        )
    if phase_function_weights.shape != phase_function_indices.shape:
        raise ValueError(
            "`phase_function_indices` and `phase_function_weights` do not have consistent shape."
        )

    optical_properties = grid.copy(deep=True)
    optical_properties['ssalb'] = (['x', 'y', 'z'], ssalb)
    optical_properties['extinction'] = (['x', 'y', 'z'], extinction)
    optical_properties['phase_weights'] = (
        ['num_micro', 'x', 'y', 'z'],
        phase_function_weights
    )
    optical_properties['legcoef'] = (['stokes_index', 'legendre_index', 'table_index'],
                                    phase_function_table)
    optical_properties['table_index'] = (
        ['num_micro', 'x', 'y', 'z'],
        1+phase_function_indices.astype(np.int32)
    )
    at3d.checks.check_optical_properties(optical_properties)
    return optical_properties

def gas_to_scatterer(gas_absorption):

    phase_function_table = np.zeros((6,2,1))
    phase_function_table[0,0,0] = 1.0
    phase_function_table[0,1,0] = 0.0
    gas_scatterer = at3d.medium.make_optical_properties(gas_absorption, gas_absorption.gas_absorption.data,
                                       ssalb=np.zeros(gas_absorption.gas_absorption.shape),
                                       phase_function_indices=np.zeros(gas_absorption.gas_absorption.shape,
                                                                     dtype=int),
                                       phase_function_table=phase_function_table)
    return gas_scatterer

class OpticalPropertyGenerator:
    """
    Transforms microphysical properties into optical properties and calculates
    derivatives of optical properties with respect to microphysical properties
    according to Mie theory.

    Parameters
    ----------
    scatterer_name : str
        The name of the scatterer for reference. e.g. cloud, aerosol.
    monodisperse_tables : Dict
        A dictionary of xr.Dataset objects which contain mie scattering properties as
        generated by at3d.mie.get_mono_table.
    size_distribution_function : callable
        The function which maps microphysical properties to a number density distribution.
        See at3d.size_distribution.gamma as an example.
    particle_density : float
        The density used in the calculation of number density distribution for normalizing
        size_distribution_functions. This is a kwarg to `size_distribution_function`.
    maxnphase : int
        The maximum number of unique phase functions to generate in the optical properties.
        If exceeded, then `interpolation_mode` will be forced to 'table' mode.
    interpolation_mode : str
        The interpolation mode to use when calculating phase functions from microphysical properties.
        If 'exact' then all phase functions will be generated up to `maxnphase`.
    density_normalization : str
        A kwarg for `size_distribution_function` that sets the normalization of the number density
        distribution.
    size_distribution_parameters : np.ndarray
        The
    """

    # TODO.
    # turn this into two different objects.
    # One does does 'linear/table'
    # The other does exact and inherits from the other one.
    # and will use some of the methods defined for 'linear/table'
    # The 'linear/table' takes poly_tables like 'nearest' which
    # opens up for generalization to optical properties for which
    # we don't have single scattering properties.

    def __init__(self, scatterer_name, monodisperse_tables, size_distribution_function,
                 particle_density=1.0, maxnphase=None,
                 interpolation_mode='exact', density_normalization='density', **size_distribution_parameters):

        if not isinstance(monodisperse_tables, typing.Dict):
            raise TypeError(
                "`monodisperse_tables` should be a dictionary of mie monodisperse "
                "tables."
            )

        self._valid_interpolation_modes = ('exact', 'nearest')
        if interpolation_mode not in self._valid_interpolation_modes:
            raise ValueError(
                "`interpolation_mode` must be in '{}' not '{}'".format(
                    self._valid_interpolation_modes, interpolation_mode)
                )
        self._interpolation_mode = interpolation_mode
        self._particle_density = particle_density
        self._maxnphase = maxnphase
        self._scatterer_name = scatterer_name

        self._density_normalization = density_normalization
        import functools
        self._size_distribution_function = functools.update_wrapper(
            functools.partial(size_distribution_function, normalization=self._density_normalization),
            size_distribution_function
        )
        if density_normalization == 'density':
            self._density_bounds = (0.0, 1e2)
        elif density_normalization == 'geometric_extinction':
            self._density_bounds = (0.0, 1e3)
        elif density_normalization == 'number_concentration':
            self._density_bounds = (0.0, 1e4)
        else:
            warnings.warn(
                "No support for default bounds for the specified `density_normalization`."
                )
            self._density_bounds = (0.0, None)


        for variable_name, parameters in size_distribution_parameters.items():
            if isinstance(parameters, np.ndarray):
                if parameters.ndim == 1:
                    size_distribution_parameters[variable_name] = {'coords': parameters}

        self._monodisperse_tables = monodisperse_tables
        self._poly_tables = OrderedDict()
        self._diff_poly_tables = OrderedDict()

        self._cached_microphysics = None
        self._cached_optical_properties = OrderedDict()

        self._size_distribution_parameters = size_distribution_parameters

        self._size_distribution_grids = OrderedDict()
        for key, monodisperse_table in self._monodisperse_tables.items():
            self._size_distribution_grids[key] = at3d.size_distribution.get_size_distribution_grid(
                monodisperse_table.radius,
                size_distribution_function=self._size_distribution_function,
                particle_density=self._particle_density,
                radius_units='micron',
                **self._size_distribution_parameters
            )

        # cached values for exact size distributions to facilitate efficient derivative calculations.
        self._table_variable_names = self._exact_size_distribution_parameters = \
        self._interpolation_weights = self._phase_indices = None

        # if no max number of phase functions is specified then we set it to the total number of table
        # entries specified.
        if self._maxnphase is None:
            self._maxnphase = list(self._size_distribution_grids.values())[0].number_density.size / \
                list(self._size_distribution_grids.values())[0].radius.size

    @property
    def scatterer_name(self):
        return self._scatterer_name

    @property
    def size_distribution_parameters(self):
        return self._size_distribution_parameters

    @property
    def size_distribution_grids(self):
        return self._size_distribution_grids

    @property
    def interpolation_mode(self):
        return self._interpolation_mode

    @property
    def size_distribution_function(self):
        return self._size_distribution_function

    @property
    def cached_optical_properties(self):
        copied = OrderedDict()
        for key,value in self._cached_optical_properties.items():
            copied[key] = value.copy(deep=True)
        return copied

    def __call__(self, microphysics):

        interp_coords = self._check_input(microphysics)
        optical_properties = OrderedDict()

        # if we are called again on the same microphysics we can just return the same
        # cached optical properties. (Maybe need to perform deep copying? Need to test.)
        # really we cache to save time when calculating derivatives afterwards.
        if microphysics.equals(self._cached_microphysics):
            optical_properties = self.cached_optical_properties#self._cached_optical_properties.copy(deep=True)

        else:
            # cache poly tables for the nearest interpolation mode to save time
            # for subsequent calculations of forward / inverse derivatives.
            if self._interpolation_mode == 'nearest' and not self._poly_tables:
                for key, table in self._monodisperse_tables.items():
                    # make sure this worked.
                    poly_table = at3d.mie.get_poly_table(self._size_distribution_grids[key], table)
                    at3d.checks.check_legendre(poly_table)
                    self._poly_tables[key] = poly_table

            if self._interpolation_mode == 'nearest':
                for key, poly_table in self._poly_tables.items():
                    single_optical_properties = self._nearest_optical_properties(
                        microphysics, poly_table, interp_coords, inverse_mode=False, interp_method='linear'
                        )
                    optical_properties[key] = self._postprocess_optical_properties(microphysics, single_optical_properties)

            elif self._interpolation_mode == 'exact':
                self._exact_size_distribution_parameters, self._interpolation_weights, self._phase_indices, \
                self._table_variable_names = self._exact_forward_size_distributions(interp_coords)
                for key, table in self._monodisperse_tables.items():
                    single_optical_properties, self._poly_tables[key] = self._exact_forward_optical_properties(
                        microphysics,
                        table,
                        self._exact_size_distribution_parameters,
                        self._interpolation_weights,
                        self._phase_indices,
                        interp_coords
                    )
                    optical_properties[key] = self._postprocess_optical_properties(microphysics, single_optical_properties)

        self._cached_optical_properties = optical_properties
        self._cached_microphysics = microphysics
        return optical_properties

    def test_valid_names(self, name):
        valid = name in self._size_distribution_parameters
        if name == 'density':
            valid = True
        return valid

    def calculate_derivatives(self, variable_names, microphysics, rel_step=0.01):

        interp_coords = self._check_input(microphysics)
        optical_derivatives = OrderedDict()

        if self._interpolation_mode == 'nearest' and not self._poly_tables:
            for key, table in self._monodisperse_tables.items():
                poly_table = at3d.mie.get_poly_table(self._size_distribution_grids[key], table)
                # make sure this worked.
                at3d.checks.check_legendre(poly_table)
                self._poly_tables[key] = poly_table

        # cache differentiated poly tables for the nearest interpolation mode to save time
        # for subsequent calculations of forward / inverse derivatives.
        if self._interpolation_mode == 'nearest' and not self._diff_poly_tables:
            for key in self._poly_tables:
                self._diff_poly_tables[key] = self._differentiate_nearest(
                    self._poly_tables[key], interp_coords, variable_names
                    )

        if self._interpolation_mode == 'nearest':
            for key, poly_table in self._diff_poly_tables.items():
                variable_optical_properties = OrderedDict()
                for name in variable_names:
                    single_optical_properties = self._nearest_optical_properties(
                        microphysics, poly_table[name],
                        interp_coords,
                        inverse_mode=name == 'density',
                        interp_method='linear'
                    )
                    variable_optical_properties[name] = self._postprocess_optical_properties(
                        microphysics, single_optical_properties
                        )
                optical_derivatives[key] = variable_optical_properties

        elif self._interpolation_mode == 'exact':
            # if we have previously processed the microphysics (e.g. for the forward problem)
            # we can make use of the exactly processed optical properties to
            # calculate the derivatives without too much extra fuss.
            # Otherwise, we have to redo the possibly expensive calculations.
            optical_properties = self._cached_optical_properties
            if not microphysics.equals(self._cached_microphysics):
                optical_properties = OrderedDict()
                self._exact_size_distribution_parameters, self._interpolation_weights,\
                self._phase_indices, self._table_variable_names = \
                    self._exact_forward_size_distributions(interp_coords)

                for key, table in self._monodisperse_tables.items():
                    single_optical_properties, self._poly_tables[key] = \
                    self._exact_forward_optical_properties(
                        microphysics,
                        table,
                        self._exact_size_distribution_parameters,
                        self._interpolation_weights,
                        self._phase_indices,
                        interp_coords
                    )
                    optical_properties[key] = self._postprocess_optical_properties(
                        microphysics, single_optical_properties
                        )
                self._cached_optical_properties = optical_properties

            for key in self._monodisperse_tables:
                optical_derivatives[key] = self._differentiate_exact(
                    microphysics, key, variable_names, interp_coords, rel_step
                    )

        return optical_derivatives

    def _differentiate_exact(self, microphysics, key, variable_names, interp_coords, rel_step=0.01):

        optical_derivatives = OrderedDict()
        for var_name in variable_names:
            # for an exact optical to microphysical relation we perform finite differencing
            # along the variable.
            if var_name == 'density':
                extinct_efficiency = np.sum(
                    self._interpolation_weights*self._poly_tables[key].extinction.data[self._phase_indices],
                    axis=0
                    ).reshape(microphysics.density.shape)
                single_optical_derivative = self._cached_optical_properties[key].copy(deep=True)
                single_optical_derivative['extinction'] = (['x','y','z'], extinct_efficiency)
                single_optical_derivative['ssalb'][:] *= 0.0
                single_optical_derivative['legcoef'][:] *= 0.0
                single_optical_derivative['derivative_method'] = 'exact'
            else:
                if var_name not in self._table_variable_names:
                    exact_size_distribution_parameters = {}
                    # step size for forward difference.
                    step = rel_step * np.sign(self._exact_size_distribution_parameters[var_name]) * \
                        np.maximum(0.05, np.abs(self._exact_size_distribution_parameters[var_name]))
                    exact_size_distribution_parameters[var_name] = step + \
                        self._exact_size_distribution_parameters[var_name]
                    for name in self._exact_size_distribution_parameters:
                        if name != var_name:
                            exact_size_distribution_parameters[name] = \
                                self._exact_size_distribution_parameters[name]

                    optical_properties_above, poly_table_above = self._exact_forward_optical_properties(
                        microphysics,
                        self._monodisperse_tables[key],
                        exact_size_distribution_parameters,
                        self._interpolation_weights,
                        self._phase_indices,
                        interp_coords
                    )

                    step_big = np.sum(
                        self._interpolation_weights*step[self._phase_indices], axis=0
                        ).reshape(microphysics.density.shape)

                    dext = (optical_properties_above.extinction -
                            self._cached_optical_properties[key].extinction)/step_big
                    dalb = (optical_properties_above.ssalb -
                            self._cached_optical_properties[key].ssalb)/step_big
                    dleg = (optical_properties_above.legcoef -
                            self._cached_optical_properties[key].legcoef)/step[None, None, :]
                    dleg[0, 0] = 0.0

                    single_optical_derivative = xr.merge(
                        [dext, dalb, dleg, self._cached_optical_properties[key].phase_weights]
                        )
                    single_optical_derivative['derivative_method'] = 'exact'

                else:
                    # this is a table variable so we have to find the derivatives
                    # based on the table nodes.
                    micro_data = interp_coords[var_name].data.ravel()
                    micro_data_bin = self._size_distribution_grids[key][var_name].data
                    digitized_micro_data = np.digitize(micro_data, bins=micro_data_bin)

                    step = micro_data_bin[digitized_micro_data] - micro_data_bin[digitized_micro_data-1]
                    diff = micro_data - micro_data_bin[digitized_micro_data-1]
                    interp1 = diff/step


                    lower_upper_flags = np.array(
                        [i for i in itertools.product(*[np.arange(-1, 1)]* \
                         len(self._table_variable_names))]
                        )
                    flag_index = np.where(np.array(self._table_variable_names) == var_name)[0][0]

                    derivative_interpolation_indices = np.zeros(
                        shape=self._interpolation_weights.shape,
                        dtype=np.int32
                    )
                    derivative_interpolation_weights = np.zeros(
                        self._interpolation_weights.shape,
                        dtype=np.float32
                    )

                    derivative_interpolation_weights[np.where(lower_upper_flags[:, flag_index] == 0), :] = \
                    self._interpolation_weights[np.where(lower_upper_flags[:, flag_index] == 0), :] / \
                    (step*interp1)
                    derivative_interpolation_weights[np.where(lower_upper_flags[:, flag_index] == -1), :] = \
                    -1*self._interpolation_weights[np.where(lower_upper_flags[:, flag_index] == 0), :] / \
                    (step*interp1)

                    derivative_interpolation_indices[:] = self._phase_indices[:]

                    # dleg has to be mixed at each point using derivative_interpolation_weights and indices.
                    # rather than calculated here.
                    dext = np.sum(self._poly_tables[key].extinction.data[derivative_interpolation_indices]*
                                  derivative_interpolation_weights, axis=0)
                    dalb = np.sum(self._poly_tables[key].ssalb.data[derivative_interpolation_indices]*
                                  derivative_interpolation_weights, axis=0)

                    single_optical_derivative = xr.Dataset(
                        data_vars={
                            'extinction': microphysics.density * dext.reshape(microphysics.density.shape),
                            'ssalb': (['x', 'y', 'z'], dalb.reshape(microphysics.density.shape)),
                            'legcoef': (['stokes_index', 'legendre_index', 'table_index'], # no differentiated legcoef required.
                                        np.zeros((6, 0, 0))),
                            'phase_weights': (['num_micro', 'x', 'y', 'z'],
                                derivative_interpolation_weights.reshape(
                                    (-1,)+microphysics.density.shape)),
                            'table_index': (['num_micro', 'x', 'y', 'z'],
                                1+derivative_interpolation_indices.reshape(
                                    (-1,)+microphysics.density.shape))
                        }
                    )
                    single_optical_derivative['derivative_method'] = 'table'

            single_optical_derivative['derivative_variable'] = var_name
            optical_derivatives[var_name] = single_optical_derivative

        return optical_derivatives

    def _check_input(self, microphysics):

        at3d.checks.check_positivity(microphysics, 'density')
        at3d.checks.check_grid(microphysics)

        interp_names = set([name for name in self._size_distribution_parameters])
        microphysics_names = set([name for name in microphysics.variables.keys()
                                  if name not in 'density'])
        missing = interp_names - microphysics_names
        if missing:
            raise KeyError(
                "microphysics dataset is missing variables "
                "for interpolation of table onto grid.", *list(missing)
                )

        number_density_grid = list(self._size_distribution_grids.values())[0]
        interp_coords = {name:microphysics[name] for name in number_density_grid.coords
                         if name not in ('radius',)}
        for interp_coord in interp_coords:
            if np.any(microphysics[interp_coord] < number_density_grid[interp_coord].min()) or \
                np.any(microphysics[interp_coord] > number_density_grid[interp_coord].max()):
                raise ValueError(
                    "Microphysical coordinate '{}' is not"
                    " within the range of the size distribution parameters.".format(interp_coord)
                    )
        return interp_coords

    def _postprocess_optical_properties(self, microphysics, optical_properties):

        # these are internal consistency checks distinct from the validity checks used for the RTE.
        assert not np.any(np.isnan(optical_properties.ssalb.data)), 'Unexpected NaN in ssalb'
        assert not np.any(np.isnan(optical_properties.extinction.data)), 'Unexpected NaN in extinction'
        assert not np.any(np.isnan(optical_properties.table_index.data)), 'Unexpected NaN in table_index'
        assert not np.any(np.isnan(optical_properties.legcoef.data)), 'Unexpected NaN in legcoef'
        assert np.allclose(optical_properties.phase_weights.sum(axis=0), 1.0), 'Interpolation weights dont add to 1.0'

        #transfer the grid variables. NOTE that delx, dely exist and be passed.
        # while nx/ny/nz are optional. delx/dely are checked for by check_grid.

        optical_properties = at3d.grid.add_grid_variables(microphysics, optical_properties)

        optical_properties = optical_properties.assign_attrs(microphysics.attrs)
        optical_properties['interp_method'] = self._interpolation_mode
        return optical_properties

    def _differentiate_nearest(self, table, interp_coords, variable_names):
        derivatives_tables = OrderedDict()
        for variable_name in variable_names:

            if variable_name in interp_coords.keys():
                differentiated = table.differentiate(coord=variable_name)

            elif variable_name == 'reff_phase_only':
                differentiated = table.differentiate(coord='reff')
                differentiated['extinction'][:] = np.zeros(table.extinction.shape)
                differentiated['ssalb'][:] = np.zeros(table.ssalb.shape)
            elif variable_name == 'veff_phase_only':
                differentiated = table.differentiate(coord='veff')
                differentiated['extinction'][:] = np.zeros(table.extinction.shape)
                differentiated['ssalb'][:] = np.zeros(table.ssalb.shape)
            elif variable_name == 'reff_no_ext':
                differentiated = table.differentiate(coord='reff')
                differentiated['extinction'][:] = np.zeros(table.extinction.shape)
            elif variable_name == 'veff_no_ext':
                differentiated = table.differentiate(coord='veff')
                differentiated['extinction'][:] = np.zeros(table.extinction.shape)
            elif variable_name == 'extinction':
                differentiated = table.copy(data={
                    'extinction': np.ones(table.extinction.shape),
                    'legcoef': np.zeros(table.legcoef.shape),
                    'ssalb': np.zeros(table.ssalb.shape)
                })
            elif variable_name == 'ssalb':
                differentiated = table.copy(data={
                    'extinction': np.zeros(table.extinction.shape),
                    'legcoef': np.zeros(table.legcoef.shape),
                    'ssalb': np.ones(table.ssalb.shape)
                })
            elif 'legendre_' in variable_name:
                leg_index = int(variable_name[len('legendre_X_'):])
                stokes_index = int(variable_name[len('legendre_')])
                legcoef = np.zeros(table.legcoef.shape)
                legcoef[stokes_index, leg_index, ...] = 1.0
                differentiated = table.copy(data={
                    'extinction': np.zeros(table.extinction.shape),
                    'legcoef': legcoef,
                    'ssalb': np.zeros(table.ssalb.shape)
                })
            elif variable_name == 'density':
                differentiated = table.copy(data={
                    'extinction': table.extinction,
                    'legcoef': np.zeros(table.legcoef.shape),
                    'ssalb': np.zeros(table.ssalb.shape)
                })
            else:
                raise ValueError("Invalid variable name to differentiate '{}''".format(variable_name))

            derivatives_tables[variable_name] = differentiated
        return derivatives_tables

    def _nearest_optical_properties(self, microphysics, poly_table, interp_coords, inverse_mode,
                              interp_method):

        # interp_method controls interpolation of optical efficiencies just for ssalb and extinction.
        ssalb = poly_table.ssalb.interp(interp_coords, method=interp_method)
        extinction_efficiency = poly_table.extinction.interp(interp_coords, method=interp_method)

        if not inverse_mode:
            extinction = extinction_efficiency * microphysics.density
        else:
            extinction = extinction_efficiency

        # Nearest neighbor for phase functions / indices.
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

        return optical_properties

    def _exact_forward_optical_properties(self, microphysics, monodisperse_table,
                                          exact_size_distribution_parameters,
                                          interpolation_weights, phase_indices, interp_coords):

        # Make the size distribution for each table's radii.
        number_density_raveled = self._size_distribution_function(
            monodisperse_table.radius,
            **exact_size_distribution_parameters,
            particle_density=self._particle_density
        )
        size_dist_attrs = OrderedDict()
        size_dist_attrs['distribution_type'] = self._size_distribution_function.__name__
        size_dist_attrs['radius_units'] = list(self._size_distribution_grids.values())[0].radius_units
        coords = {'radius': monodisperse_table.radius.data}
        coords['table_index'] = pd.MultiIndex.from_arrays(
            [exact_size_distribution_parameters[name] for name in exact_size_distribution_parameters],
            names=list(exact_size_distribution_parameters)
            )

        size_dist_grid = xr.Dataset(
            data_vars={'number_density': (['radius', 'table_index'], number_density_raveled)},
            coords=coords,
            attrs=size_dist_attrs
        )

        # Use the poly_table as the Dataset to add the main optical properties to.
        poly_table = at3d.mie.get_poly_table(size_dist_grid, monodisperse_table)


        # make sure this worked.
        at3d.checks.check_legendre(poly_table)
        optical_properties = poly_table.copy(deep=True)

        extinct_efficiency = np.sum(
            interpolation_weights*optical_properties.extinction.data[phase_indices],
            axis=0
            )
        ssalb = np.sum(interpolation_weights*optical_properties.ssalb.data[phase_indices], axis=0)
        extinction = microphysics.density*extinct_efficiency.reshape(microphysics.density.shape)
        optical_properties['extinction'] = extinction
        optical_properties['ssalb'] = (['x', 'y', 'z'], ssalb.reshape(microphysics.density.shape))

        optical_properties['table_index'] = (
            ['num_micro', 'x', 'y', 'z'],
            1+phase_indices.reshape([-1] + list(microphysics.density.shape))
        )
        optical_properties['phase_weights'] = (
            ['num_micro', 'x', 'y', 'z'],
            interpolation_weights.reshape([-1] + list(microphysics.density.shape))
        )
        optical_properties['density'] = microphysics.density
#         optical_properties.assign_coords(interp_coords)
        for name, micro_var in interp_coords.items():
            optical_properties.coords[name] = micro_var

        return optical_properties, poly_table

    def _exact_forward_size_distributions(self, interp_coords):

        # use arrays rather than dictionaries to make use of np.unique
        microphysics_data = np.stack([data.data.ravel() for data in interp_coords.values()], axis=0)
        number_density_grid = list(self._size_distribution_grids.values())[0]
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

        variable_names_corresponding_to_table_dims = []

        if num_unique <= self._maxnphase:
            # make the set of phase indices and interpolation weights.
            # when only exact phase functinos are used.
            # and specify the exact microphysical variables to calculate optical properties for.
            unique_exact, phase_indices = np.unique(microphysics_data, axis=1, return_inverse=True)

            for name, variable_data in zip(interp_coords, unique_exact):
                parameter_dict[name] = variable_data
            phase_indices = phase_indices[np.newaxis, ...]
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
            if self._maxnphase <= 0:
                table_counter = microphysics_data.shape[0]
            else:
                while num_unique > self._maxnphase:
                    index = variables_to_turn_to_table[table_counter]
                    combinations[:, index, :] = lower_upper_combinations[:, index, :]
                    unique = np.unique(
                        combinations.transpose([0, -1, 1]).reshape(
                            combinations.shape[0]*combinations.shape[-1], -1),
                        axis=0
                    )
                    num_unique = unique.shape[0]
                    table_counter += 1
                    if (table_counter > microphysics_data.shape[0]) & (num_unique > self._maxnphase):
                        raise ValueError(
                            "All variables are represented using tables but `maxnphase`='{}'' is "
                            "still exceeded num_unique='{}'. Please increase `maxnphase` or set to "
                            "`None` so that all of the table may be used.".format(self._maxnphase, num_unique)
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
                variable_names_corresponding_to_table_dims.append(name)
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

        return parameter_dict, interpolation_weights, phase_indices, variable_names_corresponding_to_table_dims

    def estimate_interpolation_error(self, nsamples=1000):

        if not self._poly_tables:
            for key, table in self._monodisperse_tables.items():
                # make sure this worked.
                poly_table = at3d.mie.get_poly_table(self._size_distribution_grids[key], table)
                at3d.checks.check_legendre(poly_table)
                self._poly_tables[key] = poly_table

        reference = OrderedDict()
        interpolated = OrderedDict()
        for key, table in self._poly_tables.items():

            monodisperse_table = self._monodisperse_tables[key]
            interp_coords = {name: table[name] for name in table.coords if name not in ('table_index', 'stokes_index')}

            random_perturbed = {name: np.random.uniform(low=interp_coords[name].min(),
                                               high=interp_coords[name].max(), size=nsamples) for name in interp_coords}

            number_density_raveled = self._size_distribution_function(
                monodisperse_table.radius,
                **random_perturbed,
                particle_density=self._particle_density
            )
            size_dist_attrs = OrderedDict()
            size_dist_attrs['distribution_type'] = self._size_distribution_function.__name__
            size_dist_attrs['radius_units'] = list(self._size_distribution_grids.values())[0].radius_units
            coords = {'radius': monodisperse_table.radius.data}
            coords['table_index'] = pd.MultiIndex.from_arrays(
                [random_perturbed[name] for name in random_perturbed],
                names=list(random_perturbed)
                )

            size_dist_grid = xr.Dataset(
                data_vars={'number_density': (['radius', 'table_index'], number_density_raveled)},
                coords=coords,
                attrs=size_dist_attrs
            )

            # Use the poly_table as the Dataset to add the main optical properties to.
            optical_properties = at3d.mie.get_poly_table(size_dist_grid, monodisperse_table)

            to_interp = {name: xr.DataArray(random_perturbed[name], dims='table_index') for name in random_perturbed}
            ext_interped = table.extinction.interp(to_interp, method='linear')
            ssalb_interped = table.ssalb.interp(to_interp, method='linear')
            legcoef_interped = table.legcoef.interp(to_interp, method='linear')

            interped = xr.merge([ext_interped, ssalb_interped, legcoef_interped])
            interpolated[key] = interped
            reference[key] = optical_properties

        return reference, interpolated


def mix_optical_properties(*scatterers, asymmetry_tol=0.002, phase_tol=0.01, linear_polarization_tol=0.01,
                           nangles=180, maxnphase=None, maxiter=1):
    """
    Mixes optical properties from different scattering species.

    Extinction and single scatter albedo are mixed exactly between species. A greedy algorithm
    is used to choose a small set of phase functions such that the relative error in the asymmetry
    parameter is less than `asym_tol` and that the relative error in the P11 component of the
    phase function doesn't exceed `phase_tol` at any of `nangles` equispaced scattering angles.
    Note that at present, all traceability information / metadata about the scatterers is lost
    upon mixing.
    Note that error tolerances are of the mixture, not absolute accuracy if the individual scatterers
    are represented approximately (e.g. by linear interpolation.)

    Parameters
    ----------
    scatterers : xr.Dataset
        Each dataset should contain the optical properties (e.g. extinction, ssalb, phase function table)
        of a particle species on a spatial grid.
        All scatterers should be on the same spatial grid and these grids should conform to
        the requirements (e.g. at3d.checks.check_grid).
    asymmetry_tol : float
        The maximum absolute error in the asymmetry parameter allowable for the subset of phase functions
        used to represent the mixture. Smaller values increase accuracy and increase number of phase
        functions required. This parameter affects flux accuracy in RTE solutions.
    phase_tol : float
        The maximum relative error in the P11 component of the phase function allowed for the subset of
        phase functions used to represent the mixture. Smaller values increase accuracy and increase number of phase
        functions required. This parameter affects radiance accuracy in RTE solutions.
    linear_polarization_tol : float
        The maximum relative error in the ratio of P12/P11 phase function componenets allowed for the subset of
        phase functions used to represent the mixture. Smaller values increase accuracy and increase number of phase
        functions required. Under the approximation that the polarization signal is dominated by single scatter,
        this
    nangles : int
        The number of scattering angles to sample the phase function at when testing `phase_tol`.
    maxnphase : int
        The maximum number of phase functions to use in the mixture.
        If None then this is set to the total number of phase functions across all scatterers.
    maxiter : int
        Number of attempts to find the subset of phase functions. If greater than 1 and maxnphase is
        exceeded then maxnphase will be extrapolated based on the number of phase functions per gridpoint
        and a new mixture will be attempted.

    Returns
    -------
    mixed_optical_properties

    Notes
    -----
    This is very similar to SHDOM's PROPGEN program though it additionally constrains the accuracy of the
    single scattering polarization signal.


    """
    #at3d.checks.check_grid_consistency(scatterers)
    for i, scatterer in enumerate(scatterers):
        at3d.checks.check_optical_properties(scatterer, i)

    max_num_micro = max(
        [scatterer.num_micro.size for scatterer in scatterers]
    )

    max_legendre = max([scatterer.sizes['legendre_index'] for
                        scatterer in scatterers])
    padded_legcoefs = [scatterer.legcoef.pad(
        {'legendre_index':
         (0, max_legendre - scatterer.legcoef.sizes['legendre_index'])
         }, constant_values=0.0)
                       for scatterer in scatterers]

    legendre_table = xr.concat(padded_legcoefs, dim='table_index')
    maxpg = scatterers[0].sizes['x']*scatterers[0].sizes['y']*scatterers[0].sizes['z']
    phase_indices = np.zeros(
        shape=[max_num_micro, maxpg, len(scatterers)],
        dtype=np.int32
    )
    phase_weights = np.zeros(
        shape=[max_num_micro, maxpg, len(scatterers)],
        dtype=np.float32
    )

    scatter_coefficients = np.zeros(
        shape=(maxpg, len(scatterers)),
        dtype=np.float32
    )
    for i, scatterer in enumerate(scatterers):
        phase_indices[..., i] = scatterer.table_index.pad(
            {'num_micro': (0, max_num_micro-scatterer.num_micro.size)},
            constant_values=1
            ).data.reshape((max_num_micro, -1), order='F') + phase_indices.max()
        phase_weights[..., i] = scatterer.phase_weights.pad(
            {'num_micro': (0, max_num_micro-scatterer.num_micro.size)},
            constant_values=0.0
            ).data.reshape((max_num_micro, -1), order='F')
        scatter_coefficients[..., i] = scatterer.extinction.data.ravel(order='F')*scatterer.ssalb.data.ravel(order='F')



    # num_micro = opt_prop.sizes['num_micro']
    nparticles = len(scatterers)
    ierr = 1
    if maxnphase is None:
        maxnphase = legendre_table.sizes['table_index']

    itr = 0
    while ierr == 1 and itr < maxiter:
        mixed_phase_table, mixed_phase_indices, nphase, ierr, errmsg = at3d.core.phase_function_mixing(
            scatter_coefficients=scatter_coefficients,
            phase_tables=legendre_table,
            phase_indices=phase_indices,
            interpolation_weights=phase_weights,
            asym_tol=asymmetry_tol,
            phase_tol=phase_tol,
            dolp_tol=linear_polarization_tol,
            nangles=nangles,
            maxnphase=maxnphase
        )

        itr += 1
        try:
            at3d.checks.check_errcode(ierr, errmsg)
        except at3d.exceptions.SHDOMError as err:
            if itr == maxiter:
                raise err

            # if we failed then try again with a larger maxnphase.
            gridpoint = int(errmsg.decode('utf-8')[100:112])
            # extrapolate based on current (nphase/ngridpoints) plus a little bit extra.
            maxnphase = int(nphase/gridpoint * scatterers[0].sizes['x']*scatterers[0].sizes['y']*scatterers[0].sizes['z']) + 10
    extinction = np.sum(np.stack([scatterer.extinction.data for scatterer in scatterers], axis=0), axis=0)
    albedo = np.sum(scatter_coefficients, axis=-1).reshape(extinction.shape,order='F')/extinction
    albedo[np.where(extinction == 0.0)] = 0.0
    mixed_phase_indices = mixed_phase_indices.reshape(extinction.shape,order='F')[np.newaxis, :]
    mixed_phase_table[0, 0] = 1.0
    mixed_optical_properties = xr.Dataset(
       data_vars={
           'legcoef': (['stokes_index', 'legendre_index', 'table_index'], mixed_phase_table[..., :nphase]),
           'phase_weights': (['num_micro', 'x', 'y', 'z'], np.ones(mixed_phase_indices.shape, dtype=np.float32)),
           'extinction': (['x', 'y', 'z'], extinction),
           'ssalb': (['x', 'y', 'z'], albedo),
           'mixture': True
       },
        coords={
            'x': scatterers[0].x,
            'y': scatterers[0].y,
            'z': scatterers[0].z,

        }
    )
    #transfer the grid variables. NOTE that delx, dely exist and be passed.
    # while nx/ny/nz are optional. delx/dely are checked for by check_grid.
    mixed_optical_properties = mixed_optical_properties.assign_coords(
        {'table_index': (['num_micro', 'x', 'y', 'z'], mixed_phase_indices)}
    )

    grid_variables = ('delx', 'dely', 'nx', 'ny', 'nz')
    for grid_variable in grid_variables:
        if grid_variable in scatterers[0].data_vars:
            mixed_optical_properties[grid_variable] = scatterers[0][grid_variable]

    mixed_optical_properties.attrs['asymmetry_tolerance'] = asymmetry_tol
    mixed_optical_properties.attrs['phase_tolerance'] = phase_tol
    mixed_optical_properties.attrs['linear_polarization_tolerance'] = linear_polarization_tol
    mixed_optical_properties.attrs['description'] = "This set of optical properties is a"\
                                                " mixture of different scattering species."\
        " All information about the species that were combined has been discarded."

    return mixed_optical_properties

def merge_optical_properties(*scatterers):
    """
    Exactly mixes the optical properties from several different scatterers.

    This method uses interpolation weights to mix the species so that the total
    number of phase functions is simply the sum of all of the species.

    This is much faster to run than the method in `medium.mix_optical_properties`,
    doesn't make any approximations and can often produce much fewer phase functions.
    The tradeoff is increased expense in the SHDOM solver. See Notes.

    Parameters
    ----------
    scatterers : xr.Dataset
        xarray datasets with valid optical properties. See `medium.make_optical_properties`.

    Returns
    -------
    merged_optical_properties : xr.Dataset
        The merged optical properties.

    Notes
    -----
    This method is useful for forward computations but has not yet been integrated into the
    inverse pipeline. This is TODO.
    The increase in computation time in the SHDOM solve is smaller than passing
    each `scatterer` to a solver.RTE object separately for legacy reasons in
    the inverse pipeline (gradient calculation) that have yet to be resolved.

    See Also
    --------
    `medium.make_optical_properties`
    `medium.mix_optical_properties`
    """
    for i, scatterer in enumerate(scatterers):
        at3d.checks.check_optical_properties(scatterer, i)

    max_num_micro = max(
        [scatterer.num_micro.size for scatterer in scatterers]
    )

    max_legendre = max([scatterer.sizes['legendre_index'] for
                        scatterer in scatterers])
    padded_legcoefs = [scatterer.legcoef.pad(
        {'legendre_index':
         (0, max_legendre - scatterer.legcoef.sizes['legendre_index'])
         }, constant_values=0.0)
                       for scatterer in scatterers]

    legendre_table = xr.concat(padded_legcoefs, dim='table_index')
    maxpg = scatterers[0].sizes['x']*scatterers[0].sizes['y']*scatterers[0].sizes['z']

    table_sizes = [padded_legcoef.table_index.size for padded_legcoef in padded_legcoefs]

    cumsum_table = 0
    padded_table_index = []
    total_scatter_coef = []
    phase_weights = []

    for table_size, scatterer in zip(table_sizes, scatterers):
        padded_table_index.append(scatterer.table_index[:].drop('table_index') + cumsum_table)
        cumsum_table += table_size
        phase_weights.append(scatterer.phase_weights.drop('table_index')*scatterer.extinction*scatterer.ssalb)
        total_scatter_coef.append(scatterer.extinction*scatterer.ssalb)

    phase_weights = xr.concat(phase_weights, dim='num_micro')
    scatter_coefs = xr.concat(total_scatter_coef,dim='num_micro').sum('num_micro')
    new_weights = phase_weights/scatter_coefs
    new_weights = new_weights.fillna(0.0)
    no_phase = np.where(new_weights.sum('num_micro').data==0.0)
    new_weights.values[0,no_phase[0],no_phase[1],no_phase[2]] = 1.0

    new_table_index = xr.concat(padded_table_index, dim='num_micro')

    extinction = sum([scatterer.extinction for scatterer in scatterers])
    albedo = sum(total_scatter_coef)/extinction
    albedo.values[np.where(extinction.values == 0.0)] = 0.0

    merged_optical_properties = xr.Dataset(
        data_vars={
            'extinction': extinction,
            'ssalb': albedo,
            'phase_weights': new_weights,
            'legcoef': legendre_table,
            'table_index': new_table_index,
        },
        coords={
            'x': extinction.x,
            'y': extinction.y,
            'z': extinction.z,
        }
    )

    merged_optical_properties = at3d.grid.add_grid_variables(scatterer, merged_optical_properties)

    at3d.checks.check_optical_properties(merged_optical_properties)

    merged_optical_properties.attrs['description'] = "This set of optical properties is a"\
                                            " merger of different scattering species."\
    " All detailed information about the species that were combined has been discarded."

    return merged_optical_properties

class GridToOpticalProperties:

    def __init__(self, rte_grid, scatterer_name, wavelength, fixed_dataset=None, *fixed_data_arrays,
                 **variable_data_bounds):

        self.scatterer_name = scatterer_name
        self.wavelength = wavelength
        at3d.checks.check_grid(rte_grid)
        self._rte_grid = rte_grid
        self.grid_shape = rte_grid.grid.shape
        self._variable_data_bounds = self._process_bounds(**variable_data_bounds)
        if fixed_dataset is None:
            for data_array in fixed_data_arrays:
                if not isinstance(data_array, (xr.Dataset, xr.DataArray)):
                    raise TypeError(
                        "`fixed_data_arrays` should be of type 'xr.Dataset' "
                        "or xr.DataArray"
                    )
            dataset = xr.merge(fixed_data_arrays)
            fixed_dataset = at3d.grid.resample_onto_grid(self._rte_grid, dataset)

        self._fixed_dataset = fixed_dataset

    def _process_bounds(self, **variable_data_bounds):

        bounds = OrderedDict()
        names = list(variable_data_bounds)
        names.extend(['extinction', 'ssalb'])
        for name in names:
            if name in variable_data_bounds:
                bound = variable_data_bounds[name]
            elif name == 'extinction':
                bound = (np.zeros(self.grid_shape)+ 0.0,
                         np.zeros(self.grid_shape) + 1e3)
            elif name == 'ssalb':
                bound = (np.zeros(self.grid_shape)+ 1e-9,
                         np.ones(self.grid_shape))
            # no supported bounds on legendre as we don't know
            # how big the table is.
            self._check_bound(bound)
            bounds[name] = bound

        return bounds

    def _check_bound(self, bounds):
        if not isinstance(bounds, typing.Tuple):
            raise TypeError(
            "Each `bound` argument should be of type '{}'"
            "".format(typing.Tuple)
            )
        if not ((len(bounds) == 2) & isinstance(bounds[0], np.ndarray) & isinstance(bounds[1], np.ndarray)):
            raise TypeError(
            "Each `bounds` should be a Tuple of two np.ndarrays."
            )
        if (bounds[0].shape != self.grid_shape) | (bounds[1].shape != self.grid_shape):
            raise ValueError(
            "Each `bound` should be of the "
            "same shape as the `rte_grid`'s spatial coordinates."
            )

    def get_bounds(self, variable_name):
        if variable_name not in self._variable_data_bounds:
            raise ValueError(
                "No bounds for '{}'".format(variable_name)
            )
        return self._variable_data_bounds[variable_name]

    def make_full_dataset(self, **variable_data):

        dataset = self._fixed_dataset.copy(deep=True)
        for name, data in variable_data.items():
            if not isinstance(data, np.ndarray):
                raise TypeError(
                    "values of argument `variable_data` should be of type '{}'"
                    "".format(np.ndarray)
                )
            if data.shape != self.grid_shape:
                raise ValueError(
                    "values of argument `variable_data` should have shape matching "
                    "the specified rte_grid '{}' not '{}'".format(
                        self.grid_shape, data.shape
                    )
                )
            # Fill any invalid data at bounds.
            # This is relevant if a mask is used.
            data[np.where(np.isnan(data))] = self._variable_data_bounds[name][0][np.where(np.isnan(data))]
            data = np.maximum(data, self._variable_data_bounds[name][0])
            data = np.minimum(data, self._variable_data_bounds[name][1])
            dataset[name] = (['x', 'y', 'z'], data)
        self._checks(dataset)
        return dataset

    def _checks(self, dataset):
        # checks for optical properties are easy.
        at3d.checks.check_optical_properties(dataset)

    def calculate_optical_properties(self, **variable_data):
        optical_properties = self.make_full_dataset(**variable_data)
        out = OrderedDict()
        out[self.wavelength] = optical_properties
        return out

    def calculate_derivatives(self, variable_names, optical_properties):

        derivatives = OrderedDict()

        for variable_name in variable_names:
            if variable_name == 'extinction':
                differentiated = xr.Dataset(
                    data_vars={
                        'extinction': (['x', 'y', 'z'],
                            np.ones(optical_properties.extinction.shape)),
                        'ssalb': (['x', 'y', 'z'],
                            np.zeros(optical_properties.ssalb.shape)),
                        'legcoef': (['stokes_index', 'legendre_index', 'table_index'],
                            np.zeros(optical_properties.legcoef.shape)),
                        'table_index': optical_properties.table_index,
                        'phase_weights': optical_properties.phase_weights,
                        'delx': optical_properties.delx,
                        'dely': optical_properties.dely
                    },
                    coords={
                        'x': optical_properties.x,
                        'y': optical_properties.y,
                        'z': optical_properties.z,
                        'stokes_index': optical_properties.stokes_index
                    }
                )
            elif variable_name == 'ssalb':
                differentiated = xr.Dataset(
                    data_vars={
                        'extinction': (['x', 'y', 'z'],
                            np.zeros(optical_properties.extinction.shape)),
                        'ssalb': (['x', 'y', 'z'],
                            np.ones(optical_properties.ssalb.shape)),
                        'legcoef': (['stokes_index', 'legendre_index', 'table_index'],
                            np.zeros(optical_properties.legcoef.shape)),
                        'table_index': optical_properties.table_index,
                        'phase_weights': optical_properties.phase_weights,
                        'delx': optical_properties.delx,
                        'dely': optical_properties.dely
                    },
                    coords={
                        'x': optical_properties.x,
                        'y': optical_properties.y,
                        'z': optical_properties.z,
                        'stokes_index': optical_properties.stokes_index
                    }
                )
            elif 'legendre_' in variable_name:
                leg_index = int(variable_name[len('legendre_X_'):])
                stokes_index = int(variable_name[len('legendre_')])
                legcoef = np.zeros(optical_properties.legcoef.shape)
                if not ((stokes_index >= 0) & (stokes_index <= 5) &
                    (leg_index <= legcoef.shape[1]) & (leg_index >= 0)):
                    raise ValueError(
                        "Invalid phase component and legendre index for variable name "
                        "'{}'".format(variable_name))
                legcoef[stokes_index, leg_index, ...] = 1.0
                differentiated = xr.Dataset(
                    data_vars={
                        'extinction': (['x', 'y', 'z'],
                            np.zeros(optical_properties.extinction.shape)),
                        'ssalb': (['x', 'y', 'z'],
                            np.zeros(optical_properties.ssalb.shape)),
                        'legcoef': (['stokes_index', 'legendre_index', 'table_index'],
                            legcoef),
                        'table_index': optical_properties.table_index,
                        'phase_weights': optical_properties.phase_weights,
                        'delx': optical_properties.delx,
                        'dely': optical_properties.dely
                    },
                    coords={
                        'x': optical_properties.x,
                        'y': optical_properties.y,
                        'z': optical_properties.z,
                        'stokes_index': optical_properties.stokes_index
                    }
                )
            else:
                raise ValueError(
                    "variable name is not supported for derivative calculation "
                    "by this generator.'{}'".format(variable_name)
                    )
            differentiated['derivative_method'] = 'exact'
            derivatives[variable_name] = differentiated

        wavelength_organized_derivatives = OrderedDict()
        wavelength_organized_derivatives[self.wavelength] = derivatives
        return wavelength_organized_derivatives


class MicrophysicsGridToOpticalProperties(GridToOpticalProperties):

    def __init__(self, rte_grid, optical_property_generator, fixed_dataset=None,
                 *fixed_data_arrays, **variable_data_bounds):
        
        self._optical_property_generator = optical_property_generator
        GridToOpticalProperties.__init__(
            self,
            rte_grid,
            optical_property_generator.scatterer_name,
            None, # a filler wavelength.
            fixed_dataset=fixed_dataset,
            *fixed_data_arrays,
            **variable_data_bounds
            )

    def calculate_optical_properties(self, **variable_data):

        microphysics_dataset = self.make_full_dataset(**variable_data)
        optical_properties = self.optical_property_generator(microphysics_dataset)
        return optical_properties

    def _process_bounds(self, **variable_data_bounds):

        # use the optical_property_generator to decide the bounds on variables.
        bounds = OrderedDict()
        opt_gen = self.optical_property_generator
        coords = list(opt_gen._size_distribution_grids.values())[0]
        variable_names = list(opt_gen._size_distribution_parameters)
        variable_names.append('density')
        for name in variable_names:
            if name in variable_data_bounds:
                bound = variable_data_bounds[name]
            elif name == 'density':
                default_bounds = opt_gen._density_bounds
                bound = (np.zeros(self.grid_shape) + default_bounds[0],
                         np.zeros(self.grid_shape) + default_bounds[1])
            else:
                bound = (np.zeros(self.grid_shape) + coords[name].min().data,
                         np.zeros(self.grid_shape) + coords[name].max().data)
            self._check_bound(bound)
            bounds[name] = bound
        return bounds

    def _checks(self, dataset):

        at3d.checks.check_grid(dataset)
        # use the optical_property_generator to verify that all the
        # necssary variables are present for the generation of optical properties.
        # this is done anyway when the same optical property generator is used.
        if 'density' not in dataset.data_vars:
            raise ValueError(
                "variable name 'density' is required to use this to generate "
                "optical properties. Please specify it at initialization or "
                "when calling this MicrophysicsGenerator."
            )
        for variable in self._optical_property_generator.size_distribution_parameters:
            if variable not in dataset:
                raise ValueError(
                    "variable '{}' is needed in order to use this to generate "
                    "optical properties. Please specify it at initialization or "
                    "when calling this MicrophysicsGenerator.".format(variable)
                )

    def calculate_derivatives(self, variable_names, microphysics):
        return self.optical_property_generator.calculate_derivatives(variable_names, microphysics)

    @property
    def optical_property_generator(self):
        return self._optical_property_generator

class TransformSet(tuple):
    """
    A simple object whose only purpose is to remove the need to remember the ordering
    of the two transforms in a Tuple.
    """
    @property
    def coordinate_transform(self):
        return self[0]

    @property
    def state_to_grid(self):
        return self[1]

class UnknownScatterer:
    """
    Specifies volumetric unknowns associated with a particular particle species
    (e.g. cloud water / aerosol / cloud ice) for retrieval.

    The unknown variables and their associated transforms can be supplied
    through the `variable_names` args and `variable_transforms` kwargs they can be
    left unspecified and added later using the UnknownScatterer.add_variable method
    which has an easier to understand signature.

    Parameters
    ----------
    grid_to_optical_properties : at3d.medium.GridToOpticalProperties
        Specifies how to map from a particular set of gridded unknowns to a complete
        set of optical properties for the particle species (ie for input as one entry
        in the `medium` argument to solver.RTE).
    variable_names : str
        The names of variables to set as unknowns to retrieve for this scatterer. This
        argument is used when no transforms are associated with the unknown variable.
    variable_transforms : Tuple
        Each kwarg key should be the variable name and the value should be a Tuple of
        a Cooordinate Transform and a StateToGrid transform. See at3d.transformsfor
        details.

    Notes
    -----
    There is some redundancy between `grid_to_optical_properties` and `variable_names`
    as `grid_to_optical_properties` should contain sufficient information to specify
    every missing variable that NEEDS to be retrieved in order to provide valid optical
    properties. The unique aspect of this object is that it is the point where transforms
    are added. Names need to be respecified to link transforms to variables. This
    explicitness also prevents mistakes happening due to misspecification of the
    `fixed_data_arrays` arguments to the `grid_to_optical_properties` object.
    """
    def __init__(self, grid_to_optical_properties, *variable_names, **variable_transforms):

        self.scatterer_name = grid_to_optical_properties.scatterer_name
        self.grid_to_optical_properties = grid_to_optical_properties
        self.variables = OrderedDict()

        for variable_name, (coordinate_transform, state_to_grid) in variable_transforms.items():
            self.add_variable(
                variable_name,
                coordinate_transform=coordinate_transform,
                state_to_grid_transform=state_to_grid
                )
        for variable_name in variable_names:
            self.add_variable(
                variable_name,
                coordinate_transform=None,
                state_to_grid_transform=None
                )

    def add_variable(self, variable_name, coordinate_transform=None,
                     state_to_grid_transform=None):
        """
        Add an unknown variable to retrieve associated with this scatterer.

        The `state_to_grid_transform` argument sets the spatial representation
        of the unknowns to reconstruct by enabling a space of possibly reduced dimension
        for the reconstruction such as a single vertical profile of unknowns or a
        single unknown per column or only a subset of grid points making up the unknowns.
        This simplifies the optimization problem by reducing its dimensionality.
        The `coordinate_transform` enables the use of an abstract set of coordinates for
        the optimization which may provide better conditioned cost function curvature
        and therefore better convergence.

        Parameters
        ----------
        variable_name : str
            The name of the variable. Should be a valid optical property name
            or the name of a microphysical variable associated with the
            `grid_to_optical_properties` object.
        coordinate_transform : at3d.transforms.CoordinateTransform
            Performs coordinate transforms to and from abstract and physical coordinates
            for this variable.
        state_to_grid_transform : at3d.transforms.StateToGridMask
            Performs a transform from gridded data to the 1D state vector of unknowns

        """
        if coordinate_transform is None:
            coordinate_transform = at3d.transforms.CoordinateTransform()
        if state_to_grid_transform is None:
            state_to_grid_transform = at3d.transforms.StateToGridMask(
                self.grid_to_optical_properties.grid_shape
                )
        self.variables[variable_name] = TransformSet(
            [coordinate_transform, state_to_grid_transform]
            )

    def get_grid_data(self, variable_name, state_subset):
        """
        Maps the subset of the abstract state vector corresponding to this unknown variable
        to gridded data in physical coordinates.

        Parameters
        ----------
        variable_name : str
            The name of the variable to form the gridded data for
        state_subset : np.ndarray, ndim=1
            The subset of the 1D state vector corresponding to this variable.

        Returns
        -------
        state_subset_on_grid : np.ndarray, ndim=3
            The gridded data in physical coordinates for this retrieved variable.
        """
        coordinate_transform, state_to_grid_transform = self.variables[variable_name]
        state_in_physical_coordinates = coordinate_transform(state_subset)
        state_subset_on_grid = state_to_grid_transform(state_in_physical_coordinates)
        return state_subset_on_grid

    def calculate_optical_properties(self, state, state_representation):
        """
        Calculate optical properties for this scatterer from the abstract state.

        This is done by identifying subsets of the state vector corresponding to
        each variable and then transforming them to gridded data which is then either
        directly used in the optical properties or used to produce optical properties.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The full abstract state vector
        state_representation : at3d.medium.StateRepresentation
            Knows the organization of different variables into the 1D state vector.

        Returns
        -------
        optical_properties : collections.OrderedDict
            Optical properties for this scatterer organized by solver key
            which is typically wavelength.
        """
        gridded_data = {}
        for variable_name in self.variables.keys():
            state_subset = state_representation.select_variable(
                state,
                self.scatterer_name,
                variable_name
                )
            gridded_data[variable_name] = self.get_grid_data(variable_name, state_subset)

        optical_properties = self.grid_to_optical_properties.calculate_optical_properties(
            **gridded_data
            )
        return optical_properties

class StateRepresentation:
    """
    This object stores information about the structure of the 1D state vector
    that is built from possibly several variables of unequal size.
    This structure is built based on the state_to_grid transforms contained in
    each `UnknownScatterer` object.

    Parameters
    ----------
    unknown_scatterers : at3d.containers.UnknownScatterers

    Notes
    -----
    This object will need to be updated to include surface unknowns.
    """
    def __init__(self, unknown_scatterers):

        grid_shape = list(unknown_scatterers.values())[0].grid_to_optical_properties.grid_shape
        # based on rte_grid we can produce faux input of the correct shape
        # and from state_to_grid.inverse we can define how large each contribution
        # to the state vector is.
        self._total_length = 0

        self._start_end_points = OrderedDict()
        for scatterer_name, unknown_scatterer in unknown_scatterers.items():
            start_end_scatterer = OrderedDict()
            for variable_name, (coordinate_transform, state_to_grid) in unknown_scatterer.variables.items():
                test_gridded_data = np.zeros(grid_shape)
                abstract_state = state_to_grid.inverse_transform(test_gridded_data)
                start_end_scatterer[variable_name] = (self._total_length, self._total_length+len(abstract_state))
                self._total_length += len(abstract_state)
            self._start_end_points[scatterer_name] = start_end_scatterer

    def update_state_vector(self, state, scatterer_name, variable_name, data_to_update):
        start, end = self._start_end_points[scatterer_name][variable_name]
        state[start:end] = data_to_update
        return state

    def select_variable(self, state, scatterer_name, variable_name):
        start, end = self._start_end_points[scatterer_name][variable_name]
        state_out = np.zeros(state[start:end].shape)
        state_out[:] = state[start:end]
        return state_out

    @property
    def start_end_points(self):
        return self._start_end_points

    @property
    def number_of_unknowns(self):
        return self._total_length


class StateGenerator:
    """
    Updates the solver.RTE objects to reflect changes in the retrieval state vector
    and projects gradients calculated by at3d.gradient.LevisApproxGradient onto
    the retrieval state vector.

    Parameters
    ----------
    `solvers_dict` : at3d.containers.SolversDict
        Contains the solver.RTE objects required to simulate the measurements.
    `unknown_scatterers` : at3d.containers.UnknownScatterers
        Contains the at3d.Medium.UnkonwnScatterer objects for all unknowns.
    `surfaces` : Dict
        The fixed surface properties used to define the solver.RTE objects during the retrieval
        organized as a function of solver key.
    `numerical_params` : Dict
        The fixed numerical properties used to define sthe solver.RTE objects  during the retrieval
        organized as a function of solver key.
    `sources` : Dict
        The fixed sources used to define the solver.RTE objects  during the retrieval
        organized as a function of solver key.
    `background_optical_scatterers` : Dict
        A nested dictionary of optical properties organized first by solver key and then
        by scatterer name which forms part of the `medium` argument to the solver.RTE
        objects.
    `num_stokes` : Dict
        The number of Stokes parameters to use in the solver.RTE objects organized
        by solver key.

    Notes
    -----
    Most of the arguments correspond to those in solver.RTE objects.
    """
    def __init__(self, solvers_dict, unknown_scatterers, surfaces,
                 numerical_parameters, sources, background_optical_scatterers,
                 num_stokes, names=None):

        # check compatibility between UnknownScatterers and all other containers.
        # if other containers don't exist then generate them from UnknownScatterers.

        if not isinstance(unknown_scatterers, at3d.containers.UnknownScatterers):
            raise TypeError(
                "`unknown_scatterers` arguments should be of type '{}'"
                " ".format(type(at3d.containers.UnknownScatterers)))
        self._unknown_scatterers = unknown_scatterers

        if not isinstance(solvers_dict, at3d.containers.SolversDict):
            raise TypeError(
                "`solvers_dict` arguments should be of type '{}'"
                " ".format(type(at3d.containers.SolversDict)))
        self._solvers_dict = solvers_dict

        self._rte_grid = list(unknown_scatterers.values())[0].grid_to_optical_properties._rte_grid
        self._grid_shape = self._rte_grid.grid.shape

        # This object defines how the state is concenated as a 1D vector.
        self._state_representation = StateRepresentation(unknown_scatterers)

        # check for compatibility in terms of keys for all solver.RTE inputs.
        # don't check the typing of the inputs, that is done when the solver is formed
        # ie when self.__call__() is executed.
        variables_to_test = [
            sources, num_stokes, surfaces, background_optical_scatterers,
            numerical_parameters
        ]
        names_to_test = [
            'sources', 'num_stokes', 'surfaces', 'background_optical_scatterers',
            'numerical_parameters'
        ]
        if names is None:
            names = OrderedDict()
            for key in sources:
                names[key] = None
        else:
            variables_to_test.append(names)
            names_to_test.append('names')

        ref_keys = None
        for variable, name in zip(variables_to_test, names_to_test):
            if not isinstance(variable, typing.Dict):
                raise TypeError(
                    "`{}` argument should be of type Dict".format(name)
                )
            if ref_keys is None:
                ref_keys = variable.keys()

            if ref_keys != variable.keys():
                raise ValueError(
                    "All of the solver.RTE arguments should have a consistent "
                    "set of keys. '{}'".format(names_to_test)
                )
            if name == 'background_optical_scatterers':
                for opt_scat_dict in background_optical_scatterers.values():
                    if not isinstance(opt_scat_dict, typing.Dict):
                        raise TypeError(
                            "Each entry in `background_optical_scatterers` argument "
                            "should be a dictionary."
                        )
                    for opt_scat in opt_scat_dict.values():
                        at3d.checks.check_optical_properties(opt_scat)
            elif name == 'num_stokes':
                for value in variable.values():
                    if not value in (1, 3, 4):
                        raise ValueError(
                            "`num_stokes` should be an integer from (1, 3, 4) "
                            "not {}".format(variable)
                        )
            else:
                for value in variable.values():
                    if value is not None:
                        if not isinstance(value, xr.Dataset):
                            raise TypeError(
                                "Each entry in `{}` should be of type '{}'".format(
                                    name, xr.Dataset
                                )
                            )

        self._sources = sources
        self._num_stokes = num_stokes
        self._surfaces = surfaces
        self._numerical_parameters = numerical_parameters
        self._background_optical_scatterers = background_optical_scatterers
        self._names = names

        for scatterer_name, unknown_scatterer in self._unknown_scatterers.items():
            if isinstance(unknown_scatterer.grid_to_optical_properties, at3d.medium.GridToOpticalProperties):
                if len(self._sources) != 1:
                    raise ValueError(
                    "Optical property unknowns are not supported for multi-spectral data. "
                    "for similar effect - try forming a 'microphysics' representation. "
                    )
        self._old_solutions = OrderedDict()
        for key in self._sources:
            self._old_solutions[key] = None

        self._total_solver_time = 0.0
        self._total_overhead_time = 0.0

    def __call__(self, state):
        """
        Update the solvers to reflect the new state vector. This does not include
        SOLVING the solver objects. That is performed elsewhere.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The 1D abstract state vector.
        """
        # first extract the previous cpu_time from the solver.
        # To get the total, we have to remember to also do this at the "end" of a retrieval.
        # To get the time for the final call. otherwise it will be biased low.
        for solver in self._solvers_dict.values():
            if hasattr(solver, '_cpu_time'):
                self._total_solver_time += solver._cpu_time

        time1 = time.process_time()
        state = self._unknown_scatterers.global_transform(state)
        new_optical_properties = OrderedDict()
        for scatterer_name, unknown_scatterer in self._unknown_scatterers.items():
            new_optical_properties[scatterer_name] = unknown_scatterer.calculate_optical_properties(
                state,
                self._state_representation
                )

        for wavelength in self._num_stokes:

            # group the optical properties for this wavelength.
            # and add the background_optical_scatterer
            background_optical_scatterer = self._background_optical_scatterers[wavelength]
            for scatterer_name, optical_properties in new_optical_properties.items():
                background_optical_scatterer[scatterer_name] = optical_properties[wavelength]

            solver = at3d.solver.RTE(
                numerical_params=self._numerical_parameters[wavelength],
                medium=background_optical_scatterer,
                source=self._sources[wavelength],
                surface=self._surfaces[wavelength],
                num_stokes=self._num_stokes[wavelength],
                name=self._names[wavelength]
            )
            if self._old_solutions[wavelength] is not None:
                solver.load_solution(self._old_solutions[wavelength])
            self._solvers_dict.add_solver(wavelength, solver)
        self._total_overhead_time += time.process_time() - time1

    def get_state(self):
        """
        Extract the current state vector from the solvers_dict.

        Returns
        -------
        state : np.ndarray, ndim=1
            The abstract state vector.
        """
        if not self._solvers_dict:
            raise ValueError(
                "State must first be set using the call method."
            )
        state = np.zeros(self._state_representation.number_of_unknowns)
        solver = list(self._solvers_dict.values())[0]
        for scatterer_name, unknown_scatterer in self._unknown_scatterers.items():
            for variable_name, (coordinate_transform, state_to_grid) in unknown_scatterer.variables.items():
                unknown_data = solver.medium[scatterer_name][variable_name].data
                state_vector = state_to_grid.inverse_transform(unknown_data)
                state_vector = coordinate_transform.inverse_transform(state_vector)
                state = self._state_representation.update_state_vector(
                    state, scatterer_name, variable_name, state_vector
                )
        state = self._unknown_scatterers.global_transform.inverse_transform(state)
        return state

    def project_gradient_to_state(self, state, gradient_dset):
        """
        Project the gradient to the abstract state.

        Parameters
        ----------
        state : np.array, ndim=1
            The abstract state.
        gradient_dset : xr.Dataset
            gridded derivatives with respect to physical unknowns.
            See at3d.gradient.make_gradient_dataset.

        Returns
        -------
        gradient : np.array
            Same shape as `state`. This is the gradient for the
            abstract state instead of the gridded physical ones.
        """
        gradient = np.zeros(self._state_representation.number_of_unknowns)
        for scatterer_name, unknown_scatterer in self._unknown_scatterers.items():
            for variable_name, (coordinate_transform, state_to_grid) in unknown_scatterer.variables.items():
                gradient_variable = gradient_dset.gradient.sel(
                    scatterer_name=scatterer_name, variable_name=variable_name
                    ).data
                gradient_state = state_to_grid.gradient_transform(gradient_variable)
                transformed_gradient = coordinate_transform.gradient_transform(state, gradient_state)
                gradient = self._state_representation.update_state_vector(
                    gradient, scatterer_name, variable_name, transformed_gradient
                )
        state = self._unknown_scatterers.global_transform.gradient_transform(state, gradient)
        return gradient

    def transform_bounds(self):
        """
        Transforms physical bounds to bounds for the 1D state vector.

        Returns
        -------
        lower_bounds_out : np.ndarray, ndim=1
            An array of lower bounds of the same shape as the state_vector
        upper_bounds_out : np.ndarray, ndim=1
            An array of upper bounds of the same shape as the state_vector
        """
        total_number_unknowns = self._state_representation.number_of_unknowns
        lower_bounds = np.zeros(total_number_unknowns)
        upper_bounds = np.zeros(total_number_unknowns)
        for scatterer_name, unknown_scatterer in self._unknown_scatterers.items():
            data_generator = unknown_scatterer.grid_to_optical_properties
            for variable_name, (coordinate_transform, state_to_grid) in unknown_scatterer.variables.items():
                lower, upper = data_generator.get_bounds(variable_name)
                for data, big_data in zip((lower, upper), (lower_bounds, upper_bounds)):
                    vector = state_to_grid.inverse_bounds_transform(data)
                    coordinate_vector = coordinate_transform.inverse_transform(vector)

                    big_data = self._state_representation.update_state_vector(
                        big_data, scatterer_name, variable_name, coordinate_vector
                    )
        lower_bounds = self._unknown_scatterers.global_transform.inverse_transform(lower_bounds)
        upper_bounds = self._unknown_scatterers.global_transform.inverse_transform(upper_bounds)

        # hack in case the transforms reverse the signs of the bounds.
        lower_bounds_out = np.minimum(lower_bounds, upper_bounds)
        upper_bounds_out = np.maximum(lower_bounds, upper_bounds)

        return lower_bounds_out, upper_bounds_out


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
    at3d.checks.check_positivity(microphysics, 'density')
    at3d.checks.check_grid(microphysics)
    if not inverse_mode:
        at3d.checks.check_legendre(poly_table)

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
    extinction_efficiency = poly_table.extinction.drop(['table_index']).interp(interp_coords, method=interp_method)

    if not inverse_mode:
        extinction = extinction_efficiency * microphysics.density
    else:
        extinction = extinction_efficiency

    extinction.name = 'extinction'

    assert not np.any(np.isnan(extinction.data)), 'Unexpected NaN in extinction'


    if not interp_coords:
        # This means that the table has no microphysical dimensions, so just
        # includes a single phase function etc so we don't do the interpolation (which will break)
        # Under these assumptions, the assertion should always pass.
        assert set(poly_table.legcoef.dims) == set(['stokes_index', 'legendre_index'])

        ssalb = xr.DataArray(
            name='ssalb',
            data=poly_table.ssalb.data*np.ones(microphysics.density.shape),
            dims=['x','y','z'],
            coords={
                'x': microphysics.x,
                'y': microphysics.y,
                'z': microphysics.z,
            }
        )
        subset_legcoef= xr.DataArray(
            name='legcoef',
            data=poly_table.legcoef.expand_dims('table_index', -1),
            dims=['stokes_index', 'legendre_index', 'table_index'],
            coords={
                'stokes_index': poly_table.stokes_index,
            }
        )

        subset_table_index = xr.DataArray(
            data=np.ones(microphysics.density.shape, dtype=int),
            dims=ssalb.dims,
            coords=ssalb.coords,
        )

    else:


        ssalb = poly_table.ssalb.drop(['table_index']).interp(interp_coords, method=interp_method)

        # Different method for phase functions / indices.
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

    assert not np.any(np.isnan(ssalb.data)), 'Unexpected NaN in ssalb'
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



# class OpticalDerivativeGenerator:
#     """
#
#     """
#     def __init__(self, scatterer_name, wavelength=None):
#         self.scatterer_name = scatterer_name
#         self.wavelength = wavelength
#
#     def test_valid_names(self, name):
#
#         valid = False
#         if name in ('extinction', 'ssalb'):
#             valid = True
#         elif 'legendre_' in name:
#             leg_index = int(name[len('legendre_X_'):])
#             stokes_index = int(name[len('legendre_')])
#             if ((stokes_index >= 0) & (stokes_index <= 5) &
#                 (leg_index >= 0)):
#                 valid = True
#         return valid
#
#     def __call__(self, data):
#         at3d.checks.check_optical_properties(data)
#         output = OrderedDict()
#         output[self.wavelength] = data
#         return output
#
#     def calculate_derivatives(self, variable_names, optical_properties):
#
#         derivatives = OrderedDict()
#
#         for variable_name in variable_names:
#             if variable_name == 'extinction':
#                 differentiated = xr.Dataset(
#                     data_vars={
#                         'extinction': (['x', 'y', 'z'],
#                             np.ones(optical_properties.extinction.shape)),
#                         'ssalb': (['x', 'y', 'z'],
#                             np.zeros(optical_properties.ssalb.shape)),
#                         'legcoef': (['stokes_index', 'legendre_index', 'table_index'],
#                             np.zeros(optical_properties.legcoef.shape)),
#                         'table_index': optical_properties.table_index,
#                         'phase_weights': optical_properties.phase_weights,
#                         'delx': optical_properties.delx,
#                         'dely': optical_properties.dely
#                     },
#                     coords={
#                         'x': optical_properties.x,
#                         'y': optical_properties.y,
#                         'z': optical_properties.z,
#                         'stokes_index': optical_properties.stokes_index
#                     }
#                 )
#             elif variable_name == 'ssalb':
#                 differentiated = xr.Dataset(
#                     data_vars={
#                         'extinction': (['x', 'y', 'z'],
#                             np.zeros(optical_properties.extinction.shape)),
#                         'ssalb': (['x', 'y', 'z'],
#                             np.ones(optical_properties.ssalb.shape)),
#                         'legcoef': (['stokes_index', 'legendre_index', 'table_index'],
#                             np.zeros(optical_properties.legcoef.shape)),
#                         'table_index': optical_properties.table_index,
#                         'phase_weights': optical_properties.phase_weights,
#                         'delx': optical_properties.delx,
#                         'dely': optical_properties.dely
#                     },
#                     coords={
#                         'x': optical_properties.x,
#                         'y': optical_properties.y,
#                         'z': optical_properties.z,
#                         'stokes_index': optical_properties.stokes_index
#                     }
#                 )
#             elif 'legendre_' in variable_name:
#                 leg_index = int(variable_name[len('legendre_X_'):])
#                 stokes_index = int(variable_name[len('legendre_')])
#                 legcoef = np.zeros(optical_properties.legcoef.shape)
#                 if not ((stokes_index >= 0) & (stokes_index <= 5) &
#                     (leg_index <= legcoef.shape[1]) & (leg_index >= 0)):
#                     raise ValueError(
#                         "Invalid phase component and legendre index for variable name "
#                         "'{}'".format(variable_name))
#                 legcoef[stokes_index, leg_index, ...] = 1.0
#                 differentiated = xr.Dataset(
#                     data_vars={
#                         'extinction': (['x', 'y', 'z'],
#                             np.zeros(optical_properties.extinction.shape)),
#                         'ssalb': (['x', 'y', 'z'],
#                             np.zeros(optical_properties.ssalb.shape)),
#                         'legcoef': (['stokes_index', 'legendre_index', 'table_index'],
#                             legcoef),
#                         'table_index': optical_properties.table_index,
#                         'phase_weights': optical_properties.phase_weights,
#                         'delx': optical_properties.delx,
#                         'dely': optical_properties.dely
#                     },
#                     coords={
#                         'x': optical_properties.x,
#                         'y': optical_properties.y,
#                         'z': optical_properties.z,
#                         'stokes_index': optical_properties.stokes_index
#                     }
#                 )
#             else:
#                 raise ValueError(
#                     "variable name is not supported for derivative calculation "
#                     "by this generator.'{}'".format(variable_name)
#                     )
#             differentiated['derivative_method'] = 'exact'
#             derivatives[variable_name] = differentiated
#
#         wavelength_organized_derivatives = OrderedDict()
#         wavelength_organized_derivatives[self.wavelength] = derivatives
#         return wavelength_organized_derivatives

# NB This is redundant with optical property generator.
# def get_optical_properties(microphysics, mie_mono_tables, size_distribution_function,
#                            particle_density=1.0, maxnphase=None,
#                            **size_distribution_grid_parameters):
#     """
#     Calculates optical properties from microphysical properties either exactly
#     or by specifying linear mixtures that are combined 'just-in-time' during the
#     RTE solution.
#
#     Optical properties are calculated on the 3D grid defined in `microphysics`
#     using microphysical parameters defined in `microphysics`. Each grid point
#     has a pointer to a single or multiple phase functions and their corresponding
#     weights with a preference for a single phase function unless maxnphase is
#     exceeded.
#
#     Parameters
#     ----------
#     microphysics : xr.Dataset
#         Should contain the microphysical parameters on a valid grid.
#         See grid.py for details. The microphysical parameters should match
#         those in specified in size_distribution_grid_parameters.
#     mie_mono_tables: Dict
#         A dictionary of Datasets of Mie legendre coefficients as a function of radius.
#         See mie.get_mono_table function for more details.
#     size_distribution_function: callable
#         Predefined function that computes size distribution.
#         Implemented options here are `gamma` (default) and `lognormal`.
#     size_distribution_grid_parameters: dict
#         Sets the spacing (and hence accuracy) of the microphysical parameters
#         for the computation of optical properties.
#         Each size_distribution_parameter dictionary should contain the following keys:
#            'coord_min': float
#                The minimum value for that coordinate.
#            'coord_max': float
#                The maximum value for that coordinate.
#            'npoints': integer
#                The number of points sampling this dimension.
#            'spacing': string
#                The type of spacing. Either 'logarithmic' or 'linear'.
#            'coord': 1D array
#                 This overrides the above arguments and specifies the exact
#                 points to sample along this dimension.
#            'units': string
#                The units of the microphysical dimension.
#         Alternatively, if a 1D numpy array is specified it will be interpreted as
#         the 'coord' argument.
#     maxnphase : int
#         Sets the maximum number of unique phase functions that can be used.
#         Default is None, in which case the maximum size is the total number
#         of points specified in the size_distribution_grid.
#
#     Returns
#     -------
#     optical_properties : OrderedDict
#         A dictionary of datasets containing optical properties defined on an SHDOM grid
#         ready to be used as input to solver.RTE.
#
#     Notes
#     -----
#     This scheme was implemented to ensure continuous variation of
#     phase functions with changing microphysics which `table_to_grid` does not
#     support. It is more similar to SHDOM's Propgen program but does not limit
#     phase functions based on sufficient accuracy. If too many phase functions
#     are required it will switch to just-in-time mixing of phase functions
#     which will increase computation time but will always be at least as
#     accurate as Propgen for a given spacing of microphysical parameters
#     (unlike `table_to_grid`) and possibly more so at the expense of
#     computation time.
#     """
#     at3d.checks.check_positivity(microphysics, 'density')
#     at3d.checks.check_grid(microphysics)
#
#     interp_names = set([name for name in size_distribution_grid_parameters])
#     microphysics_names = set([name for name in microphysics.variables.keys()
#                               if name not in 'density'])
#     missing = interp_names - microphysics_names
#     if missing:
#         raise KeyError(
#             "microphysics dataset is missing variables "
#             "for interpolation of table onto grid.", *list(missing)
#             )
#     for variable_name, parameters in size_distribution_grid_parameters.items():
#         if isinstance(parameters, np.ndarray):
#             if parameters.ndim == 1:
#                 size_distribution_grid_parameters[variable_name] = {'coords': parameters}
#     # make the grid of size distributions. This is only utilized for the computation of the
#     # coordinates so this could be sped up by isolating that code instead.
#     number_density_grid = at3d.size_distribution.get_size_distribution_grid(
#         list(mie_mono_tables.values())[0].radius,
#         size_distribution_function=size_distribution_function,
#         particle_density=particle_density,
#         radius_units='micron',
#         **size_distribution_grid_parameters
#     )
#     # if no max number of phase functions is specified then we set it to the total number of table
#     # entries specified.
#     if maxnphase is None:
#         maxnphase = number_density_grid.number_density.size / number_density_grid.radius.size
#
#     interp_coords = {name:microphysics[name] for name in number_density_grid.coords
#                      if name not in ('radius',)}
#
#     for interp_coord in interp_coords:
#         if np.any(microphysics[interp_coord] <= number_density_grid[interp_coord].min()) or \
#             np.any(microphysics[interp_coord] >= number_density_grid[interp_coord].max()):
#             raise ValueError(
#                 "Microphysical coordinate '{}' is not"
#                 " within the range of the size distribution parameters.".format(interp_coord)
#                 )
#
#     # use arrays rather than dictionaries to make use of np.unique
#     microphysics_data = np.stack([data.data.ravel() for data in interp_coords.values()], axis=0)
#
#     # digitize to find size distribution coordinate bin indices of each variable.
#     digitized_microphysics = np.stack(
#         [np.digitize(micro.data.ravel(), bins=number_density_grid[name].data)
#          for name, micro in interp_coords.items()], axis=0
#         )
#
#     # rule for forming phase functions:
#     #   1. maximize exact phase functions up to MAXNPHASE
#     #   2. if maxnphase exceeded then force to table mode.
#     #       treat table vs exact on a variable by variable basis
#     #       rather than globally or case by case.
#
#     # First check if any table is needed.
#     unique_values = np.unique(microphysics_data, axis=1)
#     num_unique = unique_values.shape[1]
#     parameter_dict = {}
#     if num_unique < maxnphase:
#         # make the set of phase indices and interpolation weights.
#         # when only exact phase functinos are used.
#         # and specify the exact microphysical variables to calculate optical properties for.
#         unique_exact, phase_indices = np.unique(microphysics_data, axis=1, return_inverse=True)
#
#         for name, variable_data in zip(interp_coords, unique_exact):
#             parameter_dict[name] = variable_data
#         phase_indices = phase_indices[np.newaxis, ...]
#         interpolation_weights = np.ones(phase_indices.shape)
#     else:
#         # lower_upper_flags contains -1,0 for finding the indices of lower and upper bounds
#         # based on the bin indices (given the default behaviour of np.digitize).
#         lower_upper_flags = np.array(
#             [i for i in itertools.product(*[np.arange(-1, 1)]*(microphysics_data.shape[0]))]
#             )
#         lower_upper_combinations = digitized_microphysics + lower_upper_flags[..., np.newaxis]
#
#         # find the exact microphysical values that we would use if we don't convert to table.
#         combinations = []
#         exact_points = []
#         for variable_data in microphysics_data:
#             unique_points, inverses = np.unique(variable_data, return_inverse=True)
#             combinations.append(inverses)
#             exact_points.append(unique_points)
#         combinations = np.repeat(np.stack(combinations, axis=0)[np.newaxis, ...],
#                                  lower_upper_flags.shape[0], axis=0)
#
#         # select which variables to first turn into table to most quickly reduce
#         # the number of phase functions.
#         unique_combinations = np.unique(lower_upper_combinations, axis=2) #number of unique bins.
#
#         # `test` is the metric for choosing which variable to turn into table.
#         #  This is not globally optimal.
#         test = np.array(
#             [np.unique(unique_combinations[:, i, :]).size/np.unique(microphysics_data[i]).size
#              for i in range(unique_combinations.shape[1])]
#             )
#
#         # convert variables to table until maxnphase is no longer exceeded.
#         variables_to_turn_to_table = np.argsort(test)
#         table_counter = 0
#         while num_unique > maxnphase:
#             index = variables_to_turn_to_table[table_counter]
#             combinations[:, index, :] = lower_upper_combinations[:, index, :]
#             unique = np.unique(
#                 combinations.transpose([0, -1, 1]).reshape(
#                     combinations.shape[0]*combinations.shape[-1], -1),
#                 axis=0
#             )
#             num_unique = unique.shape[0]
#             table_counter += 1
#             if (table_counter > microphysics_data.shape[0]) & (num_unique > maxnphase):
#                 raise ValueError(
#                     "All variables are represented using tables but `maxnphase`='{}'' is "
#                     "still exceeded num_unique='{}'. Please increase `maxnphase` or set to "
#                     "`None` so that all of the table may be used.".format(maxnphase, num_unique)
#                 )
#
#         # redo the table with the smaller number of required bounding variables.
#         # This is could be more efficiently.
#         # taking a subset of `combinations` but would be a little complicated when there
#         # are more than two table variables so I just redo.
#         # Also this re-computes the table when all variables are table even though that is now
#         # in `combinations`. There is a shape difference.
#         lower_upper_flags = np.array(
#             [i for i in itertools.product(*[np.arange(-1, 1)]*table_counter)]
#             )
#         table_combinations = digitized_microphysics[variables_to_turn_to_table[:table_counter]] + \
#             lower_upper_flags[..., np.newaxis]
#         table_unique, table_inverse = np.unique(
#             table_combinations.transpose([0, -1, 1]).reshape(
#                 table_combinations.shape[0]*table_combinations.shape[-1], -1),
#             axis=0, return_inverse=True
#             )
#
#         # make the final set of phase indices
#         phase_indices = table_inverse.reshape(
#             table_combinations.shape[0], table_combinations.shape[-1], -1
#             ).transpose([0, -1, 1])
#         phase_indices = phase_indices[:, 0, :]
#
#         # make the interpolation weights.
#         interpolation_weights = np.ones(table_combinations.shape)
#         for i in range(table_counter):
#             index = variables_to_turn_to_table[i]
#             name = list(interp_coords)[index]
#             # make the values to generate while we are at it.
#             parameter_dict[name] = number_density_grid[name].data[table_unique[:, i]]
#
#             lower = number_density_grid[name].data[table_unique[phase_indices][0, :, i]]
#             upper = number_density_grid[name].data[table_unique[phase_indices][-1, :, i]]
#             step = upper - lower
#             interp_weight = (microphysics_data[index] - lower)/step
#
#             condition = np.where(lower_upper_flags[:, i] == -1)
#             condition2 = np.where(lower_upper_flags[:, i] == 0)
#             interpolation_weights[condition, i] = (1.0 - interp_weight)
#             interpolation_weights[condition2, i] = interp_weight
#
#         interpolation_weights = np.prod(interpolation_weights, axis=1)
#
#         # add all the exact optical properties to generate.
#         for i in range(table_counter, microphysics_data.shape[0]):
#             index = variables_to_turn_to_table[i]
#             name = list(interp_coords)[index]
#             parameter_dict[name] = exact_points[index][unique[:, index]]
#
#     assert np.allclose(interpolation_weights.sum(axis=0), 1.0)
#     optical_properties = OrderedDict()
#
#     for key, mie_mono_table in mie_mono_tables.items():
#
#         number_density_raveled = size_distribution_function(
#             mie_mono_table.radius,
#             **parameter_dict,
#             particle_density=particle_density
#         )
#         size_dist_attrs = OrderedDict()
#         size_dist_attrs['distribution_type'] = size_distribution_function.__name__
#         size_dist_attrs['radius_units'] = number_density_grid.radius_units
#         coords = {'radius': mie_mono_table.radius.data}
#         coords['table_index'] = pd.MultiIndex.from_arrays(
#             [parameter_dict[name] for name in parameter_dict],
#             names=list(parameter_dict)
#             )
#
#         size_dist_grid = xr.Dataset(
#             data_vars={'number_density': (['radius', 'table_index'], number_density_raveled)},
#             coords=coords,
#             attrs=size_dist_attrs
#         )
#         poly_table = at3d.mie.get_poly_table(size_dist_grid, mie_mono_table)
#         # make sure this worked.
#         at3d.checks.check_legendre(poly_table)
#
#         extinct_efficiency = np.sum(
#             interpolation_weights*poly_table.extinction.data[phase_indices],
#             axis=0
#             )
#         ssalb = np.sum(interpolation_weights*poly_table.ssalb.data[phase_indices], axis=0)
#
#         extinction = microphysics.density*extinct_efficiency.reshape(microphysics.density.shape)
#
#         assert not np.any(np.isnan(ssalb)), 'Unexpected NaN in ssalb'
#         assert not np.any(np.isnan(extinction)), 'Unexpected NaN in extinction'
#
#         poly_table['extinction'] = extinction
#         poly_table['ssalb'] = (['x', 'y', 'z'], ssalb.reshape(microphysics.density.shape))
#
#         poly_table['table_index'] = (
#             ['num_micro', 'x', 'y', 'z'],
#             1+phase_indices.reshape([-1] + list(microphysics.density.shape))
#         )
#         poly_table['phase_weights'] = (
#             ['num_micro', 'x', 'y', 'z'],
#             interpolation_weights.reshape([-1] + list(microphysics.density.shape))
#         )
#
#         #transfer the grid variables. NOTE that delx, dely exist and be passed.
#         # while nx/ny/nz are optional. delx/dely are checked for by check_grid.
#         grid_variables = ('delx', 'dely', 'nx', 'ny', 'nz')
#         for grid_variable in grid_variables:
#             if grid_variable in microphysics.data_vars:
#                 poly_table[grid_variable] = microphysics[grid_variable]
#
#         poly_table = poly_table.assign_attrs(microphysics.attrs)
#         optical_properties[key] = poly_table
#
#     return optical_properties
