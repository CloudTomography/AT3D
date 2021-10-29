"""
This module contains objects for handling coordinate transforms for
the inverse problem.
While this can be handled at the script level, these objects cover
some useful cases and are integrated into the streamlined
inverse framework using pyshdom.medium.StateGenerator rather
than user specified 'set state functions'.

These objects apply the transformations from the abstract state level
to the spatial grid, any transformations of the state vector and also
the concatentation of the different variables in the state vector.
"""
import pyshdom
import typing
import xarray as xr
import numpy as np
from collections import OrderedDict

class CoordinateTransformNull:
    """
    The base for all coordinate transforms for the state vector.
    This performs no transformation.

    Coordinate Transforms should have the same methods with the same
    signatures as described here.
    `gradient_transform` has both `state` and `gradient` arguments to allow for
    nonlinearity in the transforms.

    Coordinate Transforms should be square, ie the length in and out
    are the same.
    """
    def __call__(self, state):
        return state

    def inverse_transform(self, state):
        return state

    def gradient_transform(self, state, gradient):
        return gradient

class CoordinateTransformLog(CoordinateTransformNull):

    def __call__(self, state):
        return np.exp(state)

    def inverse_transform(self, state):
        return np.log(state)

    def gradient_transform(self, state, gradient):
        return gradient*np.exp(state)

class CoordinateTransformScaling(CoordinateTransformNull):

    def __init__(self, offset, scaling_factor):

        self._offset = offset
        self._scaling_factor = scaling_factor

    def __call__(self, state):
        return state/self._scaling_factor + self._offset

    def inverse_transform(self, state):
        return (state - self._offset)*self._scaling_factor

    def gradient_transform(self, state, gradient):
        return self.inverse_transform(gradient)

class CoordinateTransformExp(CoordinateTransformNull):

    def __init__(self, scaling):

        self._scaling = scaling

    def __call__(self, state):
        return -np.log(1.0-state)*self._scaling

    def inverse_transform(self, state):
        return 1.0 - np.exp(-state/self._scaling)

    def gradient_transform(self, state, gradient):
        return -self._scaling*gradient/(state-1.0)

class CoordinateTransformHyperBol(CoordinateTransformNull):

    def __init__(self, scaling):

        self._scaling = scaling

    def __call__(self, state):
        return (state/(1.0-state))/self._scaling

    def inverse_transform(self, state):
        return self._scaling*state/(1.0 + self._scaling*state)

    def gradient_transform(self, state, gradient):
        return gradient/(self._scaling* ((1.0 - state)**2))

class StateToGridMask:
    """
    Transforms from gridded unknowns to a (possibly reduced) set of
    unknowns in a 1D state vector based on a grid-point mask.

    Parameters
    ----------
    grid_shape : tuple, optional
        The shape of the gridded data.

    """
    def __init__(self, grid_shape=None, mask=None):

        # add checks
        if ((mask is None) & (grid_shape is None)):
            raise ValueError(
                "At least one of `grid_shape` or `mask` arguments must be provided."
            )
        elif mask is None:
            mask = np.ones(grid_shape)
        elif grid_shape is None:
            grid_shape = mask.shape
        else:
            if grid_shape != mask.shape:
                raise ValueError(
                    "Both `grid_shape` and `mask` arguments were provided."
                    " The shape of `mask` is not consistent with `grid_shape`."
                )
        self._mask = mask
        self._grid_shape = grid_shape

    def __call__(self, state):
        gridded_state = np.zeros(self._grid_shape)
        gridded_state[np.where(self._mask)] = state
        return gridded_state

    def inverse_transform(self, gridded_data):
        return gridded_data[np.where(self._mask)]

    def gradient_transform(self, gridded_gradient):
        return gridded_gradient[np.where(self._mask)]

    def inverse_bounds_transform(self, gridded_bounds):
        return gridded_bounds[np.where(self._mask)]

class StateToGridProfile(StateToGridMask):
    """
    A single unknown per vertical level.
    Allows for a fully 3D mask.
    """
    def __call__(self, state):
        gridded_state = np.zeros(self._grid_shape)*np.nan
        for i in range(self._mask.shape[-1]):
            gridded_state[np.where(self._mask[..., i]), i] = state[i]
        return gridded_state

    def inverse_transform(self, gridded_data):
        gridded_state = np.zeros(self._grid_shape)*np.nan
        gridded_state[np.where(self._mask)] = gridded_data[self._mask]
        return np.nanmean(gridded_state, axis=(0, 1))

    def gradient_transform(self, gradient):
        return self.inverse_transform(gradient)

    def inverse_bounds_transform(self, bounds):
        # strictly this isn't true if the bounds are non-uniform
        # in space. But in the simplest case that that is true
        # then this will work.
        if np.size(np.unique(bounds)) == 1:
            out = self.inverse_transform(bounds)
        else:
            raise NotImplementedError(
                "Inverse Transform for non-uniform bounds for single variable"
                " have not yet been implemented."
            )
        return out


class StateToGrid2D(StateToGridMask):
    """
    A single unknown per column.
    Allows for a fully 3D mask.
    """
    def __call__(self, state):
        gridded_state = np.zeros(self._grid_shape)*np.nan
        for i in range(self._mask.shape[0]):
            for j in range(self._mask.shape[1]):
                gridded_state[i, j, np.where(self._mask[i, j, :])] = state.reshape(self._mask.shape[:2])[i, j]
        return gridded_state

    def inverse_transform(self, gridded_data):
        gridded_state = np.zeros(self._grid_shape)*np.nan
        gridded_state[np.where(self._mask)] = gridded_data[np.where(self._mask)]
        return np.nanmean(gridded_state, axis=-1).ravel()

    def gradient_transform(self, gradient):
        return self.inverse_transform(gradient)

    def inverse_bounds_transform(self, bounds):
        # strictly this isn't true if the bounds are non-uniform
        # in space. But in the simplest case that that is true
        # then this will work.
        if np.size(np.unique(bounds)) == 1:
            out = self.inverse_transform(bounds)
        else:
            raise NotImplementedError(
                "Inverse Transform for non-uniform bounds for single variable"
                " have not yet been implemented."
            )
        return out


# from collections import OrderedDict
# import numpy as np
#
# class Transform:
#
#     def __init__(self):
#
#         self._transforms = OrderedDict()
#         self._inverse_transforms = OrderedDict()
#         self._gradient_transforms = OrderedDict()
#         self._bounds_inverse_transforms = OrderedDict()
#
#     def _add_scatterer(self, scatterer_name):
#
#         self._transforms[scatterer_name] = OrderedDict()
#         self._inverse_transforms[scatterer_name] = OrderedDict()
#         self._gradient_transforms[scatterer_name] = OrderedDict()
#         self._bounds_inverse_transforms[scatterer_name] = OrderedDict()
#
#     @property
#     def transforms(self):
#         return self._transforms
#
#     @property
#     def inverse_transforms(self):
#         return self._transforms
#
#     @property
#     def gradient_transforms(self):
#         return self._gradient_transforms
#
#     @property
#     def bounds_inverse_transforms(self):
#         return self._bounds_inverse_transforms
#
#
# class IndependentTransform(Transform):
#     """
#     Handles transform of each variable from the grid to the 1D vector.
#     All in physical coordinates for the unknown data
#     """
#
#     def add_transform(self, scatterer_name, variable_name, transform, inverse_transform,
#                       gradient_transform, bounds_inverse_transform):
#
#         if scatterer_name not in self._transforms:
#             self._add_scatterer(scatterer_name)
#
#         self._transforms[scatterer_name][variable_name] = transform
#         self._inverse_transforms[scatterer_name][variable_name] = inverse_transform
#         self._gradient_transforms[scatterer_name][variable_name] = gradient_transform
#         self._bounds_inverse_transforms[scatterer_name][variable_name] = bounds_inverse_transform
#
#     def __call__(self, state, scatterer_name, variable_name):
#         return self._transforms[scatterer_name][variable_name](state)
#
#     def inverse(self, gridded_data, scatterer_name, variable_name):
#         return self._inverse_transforms[scatterer_name][variable_name](gridded_data)
#
#     def calc_derivative(self, gridded_gradient, scatterer_name, variable_name):
#         return self._gradient_transforms[scatterer_name][variable_name](gridded_gradient)
#
#     def inverse_bounds(self, gridded_bounds, scatterer_name, variable_name):
#         return self._bounds_inverse_transforms[scatterer_name][variable_name](gridded_bounds)
#
# class StateToGridMask(IndependentTransform):
#
#     def add_transform(self, scatterer_name, variable_name, mask, fill_value):
#         mask = mask.astype(np.bool)
#         if scatterer_name not in self._transforms:
#             self._add_scatterer(scatterer_name)
#
#         def transform(state):
#             state_gridded = np.zeros(mask.shape) + fill_value
#             state_gridded[np.where(mask)] = state
#             return state_gridded
#         self._transforms[scatterer_name][variable_name] = transform
#
#         def inverse_transform(gridded_data):
#             return gridded_data[mask]
#
#         self._inverse_transforms[scatterer_name][variable_name] = inverse_transform
#         self._gradient_transforms[scatterer_name][variable_name] = inverse_transform
#         self._bounds_inverse_transforms[scatterer_name][variable_name] = inverse_transform
#
# class StateToGridNull(IndependentTransform):
#
#     def add_transform(self, scatterer_name, variable_name, rte_grid):
#
#         grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
#         mask = np.ones(grid_shape, dtype=np.bool)
#
#         if scatterer_name not in self._transforms:
#             self._add_scatterer(scatterer_name)
#
#         def transform(state):
#             state_gridded = np.zeros(mask.shape)
#             state_gridded[mask] = state
#             return state_gridded
#         self._transforms[scatterer_name][variable_name] = transform
#
#         def inverse_transform(gridded_data):
#             return gridded_data[mask]
#
#         self._inverse_transforms[scatterer_name][variable_name] = inverse_transform
#         self._gradient_transforms[scatterer_name][variable_name] = inverse_transform
#         self._bounds_inverse_transforms[scatterer_name][variable_name] = inverse_transform
#
# class StateRepresentation:
#     """
#     This only represents the physical coordinates.
#     In abstract/general coordinates, there is not necessarily any
#     separation by physical variable name.
#     """
#     def __init__(self, state_to_grid, rte_grid):
#
#         # based on rte_grid we can produce faux input of the correct shape
#         # and from state_to_grid.inverse we can define how large each contribution
#         # to the state vector is.
#         grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
#         self._total_length = 0
#
#         self._start_end_points = OrderedDict()
#         for scatterer_name, variable_name_data in state_to_grid.transforms.items():
#             start_end_scatterer = OrderedDict()
#             for variable_name in variable_name_data:
#                 test_gridded_data = np.zeros(grid_shape)
#                 abstract_state = state_to_grid.inverse(test_gridded_data, scatterer_name, variable_name)
#                 start_end_scatterer[variable_name] = (self._total_length, self._total_length+len(abstract_state))
#                 self._total_length += len(abstract_state)
#             self._start_end_points[scatterer_name] = start_end_scatterer
#
#     def update_state_vector(self, state, scatterer_name, variable_name, data_to_update):
#         start, end = self._start_end_points[scatterer_name][variable_name]
#         state[start:end] = data_to_update
#         return state
#
#     def select_variable(self, state, scatterer_name, variable_name):
#         start, end = self._start_end_points[scatterer_name][variable_name]
#         state_out = np.zeros(state[start:end].shape)
#         state_out[:] = state[start:end]
#         return state_out
#
#     @property
#     def start_end_points(self):
#         return self._start_end_points
#
#     @property
#     def number_of_unknowns(self):
#         return self._total_length
#
# class StateTransform:
#     """
#     Generalized coordinate transforms between 1D state and gradient vectors
#     and the abstract state and its gradient. This is assumed to be non-destructive.
#     ie length of the state before and after is the same.
#     """
#     def __init__(self, state_representation, transform, inverse_transform, gradient_transform):
#         self._state_representation = state_representation
#         self._transform = transform
#         self._inverse_transform = inverse_transform
#         self._gradient_transform = gradient_transform
#
#     def __call__(self, state):
#         return self._transform(self._state_representation, state)
#
#     def inverse(self, state):
#         return self._inverse_transform(self._state_representation, state)
#
#     def calc_derivative(self, state, gradient):
#         return self._gradient_transform(self._state_representation, state, gradient)
#
#     @property
#     def transforms(self):
#         return self._transform
#
#     @property
#     def inverse_transforms(self):
#         return self._inverse_transform
#
#     @property
#     def gradient_transforms(self):
#         return self._gradient_transform
#
# class IndependentStateTransform(IndependentTransform):
#     """
#     We assume that the different variable names don't interact in the transform.
#     Only selected variables are transformed. Each variable may have a select
#     transform.
#     """
#     def __init__(self, state_representation):
#         super().__init__()
#
#         self._state_representation = state_representation
#
#     def __call__(self, state):
#
#         state_copy = np.zeros(state.shape)
#         state_copy[:] = state[:]
#         for scatterer_name, variable_name_data in self._transforms.items():
#             for variable_name in variable_name_data:
#                 state_vector_portion = self._state_representation.select_variable(state, scatterer_name, variable_name)
#                 state_copy = self._state_representation.update_state_vector(
#                     state_copy,
#                     scatterer_name,
#                     variable_name,
#                     self._transforms[scatterer_name][variable_name](state_vector_portion)
#                 )
#         return state_copy
#
#     def inverse(self, state):
#         state_copy = np.zeros(state.shape)
#         state_copy[:] = state[:]
#         for scatterer_name, variable_name_data in self._transforms.items():
#             for variable_name in variable_name_data:
#                 state_vector_portion = self._state_representation.select_variable(state, scatterer_name, variable_name)
#                 state_copy = self._state_representation.update_state_vector(
#                     state_copy,
#                     scatterer_name,
#                     variable_name,
#                     self._inverse_transforms[scatterer_name][variable_name](state_vector_portion)
#                 )
#         return state_copy
#
#     def calc_derivative(self, state, gradient):
#         gradient_copy = np.zeros(gradient.shape)
#         gradient_copy[:] = gradient[:]
#         for scatterer_name, variable_name_data in self._transforms.items():
#             for variable_name in variable_name_data:
#                 state_vector_portion = self._state_representation.select_variable(state, scatterer_name, variable_name)
#                 gradient_vector_portion = self._state_representation.select_variable(gradient, scatterer_name, variable_name)
#                 gradient_copy = self._state_representation.update_state_vector(
#                     gradient_copy,
#                     scatterer_name,
#                     variable_name,
#                     self._gradient_transforms[scatterer_name][variable_name](state_vector_portion, gradient_vector_portion)
#                 )
#         return gradient_copy
#
# class NullStateTransform(IndependentStateTransform):
#
#     def add_transform(self, scatterer_name, variable_name):
#
#         if scatterer_name not in self._transforms:
#             self._add_scatterer(scatterer_name)
#         def null(state):
#             return state
#         def gradient_null(state, gradient):
#             return gradient
#
#         self._transforms[scatterer_name][variable_name] = null
#         self._inverse_transforms[scatterer_name][variable_name] = null
#         self._gradient_transforms[scatterer_name][variable_name] = gradient_null
#
# class LogStateTransform(IndependentStateTransform):
#     """
#     A special case where all variables are log scaled.
#     """
#     def add_transform(self, scatterer_name, variable_name):
#
#         if scatterer_name not in self._transforms:
#             self._add_scatterer(scatterer_name)
#
#         self._transforms[scatterer_name][variable_name] = np.exp
#         self._inverse_transforms[scatterer_name][variable_name] = np.log
#         def gradient_transform(state, gradient):
#             return gradient*np.exp(state)
#         self._gradient_transforms[scatterer_name][variable_name] = gradient_transform
