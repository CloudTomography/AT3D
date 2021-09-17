"""
A module for storing regularization routines.
"""
import numpy as np
import pandas as pd
import xarray as xr
import pyshdom

class Regularization:

    def __init__(self, state_generator, scatterer_name, variable_name,
                 regularization_strength):

        if (not isinstance(regularization_strength, np.float)) or (regularization_strength < 0.0):
            raise ValueError("`regularization_strength` should be a positive float.")
        self._regularization_strength = regularization_strength

        if not isinstance(state_generator, pyshdom.medium.StateGenerator):
            raise TypeError(
                "`state_generator` should be of type {} not {}".format(
                    pyshdom.medium.StateGenerator,
                    type(state_generator)
                    )
                )
        self._state_generator = state_generator
        if scatterer_name not in self.state_generator._unknown_scatterers:
            raise KeyError(
                "Invalid unknown scatterer name `{}`. Valid names are: {}".format(
                    scatterer_name, self.state_generator._unknown_scatterers.keys()
                )
            )
        if variable_name not in \
        self.state_generator._unknown_scatterers[scatterer_name]['variable_name_list']:
            raise KeyError(
                "Invalid unknown variable name `{}`. Valid names are: {}".format(
                    scatterer_name,
                    self.state_generator._unknown_scatterers[scatterer_name]['variable_name_list']
                )
            )

        self._scatterer_name = scatterer_name
        self._variable_name = variable_name

    def _get_gridded_variable(self, state):
        # This is predicated on the idea that the measurement misfit term
        # has always been called first so that the up to date gridded state is
        # available from the state_generator, which holds the solvers.
        # TO BE SAFE we could call the state_generator here again.
        # e.g. self.state_generator(state)
        solver = list(self.state_generator._solvers_dict.values())[0]
        gridded_data = solver.medium[self._scatterer_name][self._variable_name].data
        rte_grid = self.state_generator._rte_grid
        return gridded_data, rte_grid

    def _make_gradient_dset(self, state, gradient):
        # Form a gradient dataset analagous to the observation misfit so the
        # same project_gradient_to_state can be used.
        rte_grid = self.state_generator._rte_grid
        derivative_index = pd.MultiIndex.from_arrays(
            [np.atleast_1d(self._scatterer_name), np.atleast_1d(self._variable_name)],
            names=("scatterer_name", "variable_name"))

        gradient_dataset = xr.Dataset(
            data_vars={
                'gradient': (
                    ['x', 'y', 'z', 'derivative_index'],
                    gradient.reshape(
                        (rte_grid.sizes['x'], rte_grid.sizes['y'], rte_grid.sizes['z'], -1))
                )
            },
            coords={
                'x': rte_grid.x,
                'y': rte_grid.y,
                'z': rte_grid.z,
                'derivative_index': derivative_index
            }
        )
        out_gradient = self.state_generator.project_gradient_to_state(
            state, gradient_dataset
        )
        return out_gradient

    @property
    def state_generator(self):
        return self._state_generator

class WeightedRegularization(Regularization):

    def __init__(
        self, state_generator, scatterer_name, variable_name,
        regularization_strength, spatial_weights='uniform'):

        Regularization.__init__(self, state_generator, scatterer_name, variable_name,
        regularization_strength)

        if spatial_weights == 'uniform':
            spatial_weights = np.ones(self.state_generator._grid_shape)

        if not isinstance(spatial_weights, np.ndarray):
            raise TypeError("`spatial_weights` should be a numpy array.")
        if spatial_weights.shape != self.state_generator._grid_shape:
            raise ValueError(
                "`spatial_weights` should be of the same shape as the property grid."
            )
        if np.any(spatial_weights < 0.0):
            raise ValueError(
                "`spatial_weights` should be positive floats."
            )
        self._spatial_weights = spatial_weights

class Sparsity(WeightedRegularization):

    def __init__(
        self, state_generator, scatterer_name, variable_name,
        regularization_strength, spatial_weights='uniform'):
        WeightedRegularization.__init__(self, state_generator, scatterer_name, variable_name,
        regularization_strength, spatial_weights)

    def __call__(self, state):

        gridded_data, rte_grid = self._get_gridded_variable(state)

        # should this just be a sum or should it be weighted by cell_volume
        # to account for non-uniform grid sizes?
        # If the latter then it will be easiest to add this as a case to the
        # fortran grid_smoothing code in util.f90.
        gradient = self._regularization_strength*np.sign(gridded_data)
        cost = self._regularization_strength*np.sum(np.abs(gridded_data))

        out_gradient = self._make_gradient_dset(state, gradient)

        return cost, out_gradient

class Tikhonov(WeightedRegularization):

    def __init__(
        self, state_generator, scatterer_name, variable_name,
        regularization_strength, spatial_weights='uniform'):
        WeightedRegularization.__init__(self, state_generator, scatterer_name, variable_name,
        regularization_strength, spatial_weights)

    def __call__(self, state):

        gridded_data, rte_grid = self._get_gridded_variable(state)

        # should this just be a sum or should it be weighted by cell_volume
        # to account for non-uniform grid sizes?
        # If the latter then it will be easiest to add this as a case to the
        # fortran grid_smoothing code in util.f90.
        gradient = self._regularization_strength*2*gridded_data
        cost = self._regularization_strength*np.sum(gridded_data**2)

        out_gradient = self._make_gradient_dset(state, gradient)

        return cost, out_gradient

# class CorrelationSmoothing(WeightedRegularization):
#
#     def __init__(
#         self, state_generator, scatterer_name, variable_name,
#         regularization_strength, spatial_weights='uniform'):
#         WeightedRegularization.__init__(self, state_generator, scatterer_name, variable_name,
#         regularization_strength, spatial_weights)
#


class SpatialSmoothing(WeightedRegularization):
    """
    Calculates the L1 or L2 norm of the spatial gradient of a selected variable.
    """
    def __init__(
            self, state_generator, scatterer_name, variable_name,
            regularization_strength, mode='l2', direction_weights=[1.0, 1.0, 1.0],
            spatial_weights='uniform'):
        WeightedRegularization.__init__(self, state_generator, scatterer_name, variable_name,
        regularization_strength, spatial_weights)

        valid_modes = ('l1', 'l2')
        if mode not in valid_modes:
            raise ValueError(
                "Spatial Smoothing `mode` should be in {}".format(valid_modes)
            )
        self._mode = mode

        direction_weights = np.atleast_1d(direction_weights)
        self._direction_weights = direction_weights

    def __call__(self, state):

        gridded_data, rte_grid = self._get_gridded_variable(state)

#       The x,y,z conventions are transposed in this function compared
#       to what is used by pyshdom, that is why the names don't match up.
        cost, gradient, ierr, errmsg = pyshdom.core.grid_smoothing(
            field=gridded_data,
            weights=self._spatial_weights,
            direction_weights=self._direction_weights,
            nx=gridded_data.shape[2],
            ny=gridded_data.shape[1],
            nz=gridded_data.shape[0],
            zgrid=rte_grid.x.data,
            ygrid=rte_grid.y.data,
            xgrid=rte_grid.z.data,
            mode=self._mode
        )
        pyshdom.checks.check_errcode(ierr, errmsg)
        cost *= self._regularization_strength
        gradient *= self._regularization_strength

        # transform the gradient numpy array to an xarray dset that is then
        # projected to state.
        out_gradient = self._make_gradient_dset(state, gradient)
        return cost, out_gradient


# The function below is outdated.
# class LocalCorrelation:
#
#     def __init__(self, grid, state_on_grid_indicators, correlation_length, state_start=None, state_end=None):
#
#         from scipy.spatial import distance_matrix
#         #TODO properly check inputs.
#         self._state_start, self._state_end = state_start, state_end
#
#         a,b,c = state_on_grid_indicators
#         xs = grid.x.data[a]
#         ys = grid.y.data[b]
#         zs = grid.z.data[c]
#         positions = np.stack((xs,ys,zs),axis=1)
#         dist_mat = distance_matrix(positions, positions)
#
#         self._grid = grid
#         self._dist_mat = dist_mat
#         self._correl_matrix = np.exp(-dist_mat/correlation_length)
#
#     def __call__(self, state):
#
#         grad_vect = np.zeros(state.shape)
#         if self._state_start is None and self._state_end is None:
#             state_start = 0
#             state_end = len(state)
#         else:
#             state_start, state_end = self._state_start, self._state_end
#         temp = state[state_start:state_end]
#         residual = temp - np.dot(self._correl_matrix, temp)
#         cost = np.dot(residual, np.dot(np.diag(np.ones(temp.shape)), residual)) / (temp.size**2)
#         grad = 2*np.dot(residual, np.diag(np.ones(temp.shape))) / temp.size
#         grad_vect[state_start:state_end] = grad
#         return cost, grad_vect
