"""

"""
from collections import OrderedDict
import numpy as np

class Transform:

    def __init__(self):

        self._transforms = OrderedDict()
        self._inverse_transforms = OrderedDict()
        self._gradient_transforms = OrderedDict()
        self._bounds_inverse_transforms = OrderedDict()

    def _add_scatterer(self, scatterer_name):

        self._transforms[scatterer_name] = OrderedDict()
        self._inverse_transforms[scatterer_name] = OrderedDict()
        self._gradient_transforms[scatterer_name] = OrderedDict()
        self._bounds_inverse_transforms[scatterer_name] = OrderedDict()

    @property
    def transforms(self):
        return self._transforms

    @property
    def inverse_transforms(self):
        return self._transforms

    @property
    def gradient_transforms(self):
        return self._gradient_transforms

    @property
    def bounds_inverse_transforms(self):
        return self._bounds_inverse_transforms


class IndependentTransform(Transform):
    """
    Handles transform of each variable from the grid to the 1D vector.
    All in physical coordinates for the unknown data
    """

    def add_transform(self, scatterer_name, variable_name, transform, inverse_transform,
                      gradient_transform, bounds_inverse_transform):

        if scatterer_name not in self._transforms:
            self._add_scatterer(scatterer_name)

        self._transforms[scatterer_name][variable_name] = transform
        self._inverse_transforms[scatterer_name][variable_name] = inverse_transform
        self._gradient_transforms[scatterer_name][variable_name] = gradient_transform
        self._bounds_inverse_transforms[scatterer_name][variable_name] = bounds_inverse_transform

    def __call__(self, state, scatterer_name, variable_name):
        return self._transforms[scatterer_name][variable_name](state)

    def inverse(self, gridded_data, scatterer_name, variable_name):
        return self._inverse_transforms[scatterer_name][variable_name](gridded_data)

    def calc_derivative(self, gridded_gradient, scatterer_name, variable_name):
        return self._gradient_transforms[scatterer_name][variable_name](gridded_gradient)

    def inverse_bounds(self, gridded_bounds, scatterer_name, variable_name):
        return self._bounds_inverse_transforms[scatterer_name][variable_name](gridded_bounds)

class StateToGridMask(IndependentTransform):

    def add_transform(self, scatterer_name, variable_name, mask, fill_value):

        if scatterer_name not in self._transforms:
            self._add_scatterer(scatterer_name)

        def transform(state):
            state_gridded = np.zeros(mask.shape) + fill_value
            state_gridded[mask] = state
            return state_gridded
        self._transforms[scatterer_name][variable_name] = transform

        def inverse_transform(gridded_data):
            return gridded_data[mask]

        self._inverse_transforms[scatterer_name][variable_name] = inverse_transform
        self._gradient_transforms[scatterer_name][variable_name] = inverse_transform
        self._bounds_inverse_transforms[scatterer_name][variable_name] = inverse_transform

class StateToGridNull(IndependentTransform):

    def add_transform(self, scatterer_name, variable_name, rte_grid):

        grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
        mask = np.ones(grid_shape, dtype=np.bool)

        if scatterer_name not in self._transforms:
            self._add_scatterer(scatterer_name)

        def transform(state):
            state_gridded = np.zeros(mask.shape)
            state_gridded[mask] = state
            return state_gridded
        self._transforms[scatterer_name][variable_name] = transform

        def inverse_transform(gridded_data):
            return gridded_data[mask]

        self._inverse_transforms[scatterer_name][variable_name] = inverse_transform
        self._gradient_transforms[scatterer_name][variable_name] = inverse_transform
        self._bounds_inverse_transforms[scatterer_name][variable_name] = inverse_transform

class StateRepresentation:
    """
    This only represents the physical coordinates.
    In abstract/general coordinates, there is not necessarily any
    separation by physical variable name.
    """
    def __init__(self, state_to_grid, rte_grid):

        # based on rte_grid we can produce faux input of the correct shape
        # and from state_to_grid.inverse we can define how large each contribution
        # to the state vector is.
        grid_shape = (rte_grid.x.size, rte_grid.y.size, rte_grid.z.size)
        self._total_length = 0

        self._start_end_points = OrderedDict()
        for scatterer_name, variable_name_data in state_to_grid.transforms.items():
            start_end_scatterer = OrderedDict()
            for variable_name in variable_name_data:
                test_gridded_data = np.zeros(grid_shape)
                abstract_state = state_to_grid.inverse(test_gridded_data, scatterer_name, variable_name)
                start_end_scatterer[variable_name] = (self._total_length, self._total_length+len(abstract_state)+1)
                self._total_length += len(abstract_state)
            self._start_end_points[scatterer_name] = start_end_scatterer

    def update_state_vector(self, state, scatterer_name, variable_name, data_to_update):
        start, end = self._start_end_points[scatterer_name][variable_name]
        state[start:end] = data_to_update
        return state

    def select_variable(self, state, scatterer_name, variable_name):
        start, end = self._start_end_points[scatterer_name][variable_name]
        return state[start:end]

    @property
    def start_end_points(self):
        return self._start_end_points

    @property
    def number_of_unknowns(self):
        return self._total_length

class StateTransform:
    """
    Generalized coordinate transforms between 1D state and gradient vectors
    and the abstract state and its gradient. This is assumed to be non-destructive.
    ie length of the state before and after is the same.
    """
    def __init__(self, state_representation, transform, inverse_transform, gradient_transform):
        self._state_representation = state_representation
        self._transform = transform
        self._inverse_transform = inverse_transform
        self._gradient_transform = gradient_transform

    def __call__(self, state):
        return self._transform(self._state_representation, state)

    def inverse(self, state):
        return self._inverse_transform(self._state_representation, state)

    def calc_derivative(self, state, gradient):
        return self._gradient_transform(self._state_representation, state, gradient)

    @property
    def transforms(self):
        return self._transform

    @property
    def inverse_transforms(self):
        return self._inverse_transform

    @property
    def gradient_transforms(self):
        return self._gradient_transform

class IndependentStateTransform(IndependentTransform):
    """
    We assume that the different variable names don't interact in the transform.
    Only selected variables are transformed. Each variable may have a select
    transform.
    """
    def __init__(self, state_representation):
        super().__init__()

        self._state_representation = state_representation

    def __call__(self, state):

        for scatterer_name, variable_name_data in self._transforms.items():
            for variable_name in variable_name_data:
                state_vector_portion = self._state_representation.select_variable(state, scatterer_name, variable_name)
                state = self._state_representation.update_state_vector(
                    state,
                    scatterer_name,
                    variable_name,
                    self._transforms[scatterer_name][variable_name](state_vector_portion)
                )
        return state

    def inverse(self, state):

        for scatterer_name, variable_name_data in self._transforms.items():
            for variable_name in variable_name_data:
                state_vector_portion = self._state_representation.select_variable(state, scatterer_name, variable_name)
                state = self._state_representation.update_state_vector(
                    state,
                    scatterer_name,
                    variable_name,
                    self._inverse_transforms[scatterer_name][variable_name](state_vector_portion)
                )
        return state

    def calc_derivative(self, state, gradient):

        for scatterer_name, variable_name_data in self._transforms.items():
            for variable_name in variable_name_data:
                state_vector_portion = self._state_representation.select_variable(state, scatterer_name, variable_name)
                gradient_vector_portion = self._state_representation.select_variable(gradient, scatterer_name, variable_name)
                gradient = self._state_representation.update_state_vector(
                    gradient,
                    scatterer_name,
                    variable_name,
                    self._gradient_transforms[scatterer_name][variable_name](state_vector_portion, gradient_vector_portion)
                )
        return gradient

class NullStateTransform(IndependentStateTransform):

    def add_transform(self, scatterer_name, variable_name):

        if scatterer_name not in self._transforms:
            self._add_scatterer(scatterer_name)
        def null(state):
            return state
        def gradient_null(state, gradient):
            return gradient

        self._transforms[scatterer_name][variable_name] = null
        self._inverse_transforms[scatterer_name][variable_name] = null
        self._gradient_transforms[scatterer_name][variable_name] = gradient_null

class LogStateTransform(IndependentStateTransform):
    """
    A special case where all variables are log scaled.
    """
    def add_transform(self, scatterer_name, variable_name):

        if scatterer_name not in self._transforms:
            self._add_scatterer(scatterer_name)

        self._transforms[scatterer_name][variable_name] = np.exp
        self._inverse_transforms[scatterer_name][variable_name] = np.log
        def gradient_transform(state, gradient):
            return gradient*np.exp(state)
        self._gradient_transforms[scatterer_name][variable_name] = gradient_transform
