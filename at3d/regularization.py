"""
A module for regularization routines.
"""
import numpy as np
import pandas as pd
import xarray as xr
import scipy.special as ss
import at3d

class Regularization:
    """
    Base class for regularization routines.

    Parameters
    ----------
    state_generator : at3d.medium.StateGenerator
        This object contains the methods to map from the abstract state to the
        physical coordinates.
    scatterer_name : str
        The name of the scatterer that the regularization term is being applied to. This
        should match one of the keys in the solver.RTE.medium that is within the
        `state_generator`.
    variable_name : str
        The physical variable name that the regularization term is being applied to.
    regularization_strength : float
        A hyperparameter that scales the cost function of the regularization term.
    relaxation_parameter : float
        A hyperparameter controls the geometric sequence of `regularization_strength`
        as a function of optimization iteration number. If this is unity the
        effective regularization strength remains constant throughout the optimization.
        The closer this parameter is to zero the faster the regularization decreases.
    """
    def __init__(self, state_generator, scatterer_name, variable_name,
                 regularization_strength, relaxation_parameter=1.0):

        if (not isinstance(regularization_strength, float)) or (regularization_strength < 0.0):
            raise ValueError("`regularization_strength` should be a positive float.")
        self._regularization_strength = regularization_strength

        if (not isinstance(relaxation_parameter, float)) or (relaxation_parameter < 0.0):
            raise ValueError("`relaxation_parameter` should be a positive float.")
        self._relaxation_parameter = relaxation_parameter

        if not isinstance(state_generator, at3d.medium.StateGenerator):
            raise TypeError(
                "`state_generator` should be of type {} not {}".format(
                    at3d.medium.StateGenerator,
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
        self.state_generator._unknown_scatterers[scatterer_name].variables.keys():
            raise KeyError(
                "Invalid unknown variable name `{}`. Valid names are: {}".format(
                    scatterer_name,
                    self.state_generator._unknown_scatterers[scatterer_name].variables.keys()
                )
            )

        self._scatterer_name = scatterer_name
        self._variable_name = variable_name

    def _get_gridded_variable(self, state):
        """
        Return gridded_data from abstract_state for this variable.

        Parameters
        ----------
        state : np.array, ndim=1
            The abstract 1D state vector.

        Returns
        -------
        gridded_data : np.ndarray, ndim=3
            The gridded data correpsonding for the scatterer and variable
            of this regularization term.
        rte_grid : xr.Dataset
            The dataset containing the coordinate information for the gridded data.
        """
        # This is predicated on the idea that the measurement misfit term
        # has always been called first so that the up to date gridded state is
        # available from the state_generator, which holds the solvers.
        #
        # TO BE SAFE we could instead calculate directly from `state`
        # using transform information contained in self.state_generator._unknown_scatterers.
        solver = list(self.state_generator._solvers_dict.values())[0]
        gridded_data = solver.medium[self._scatterer_name][self._variable_name].data
        rte_grid = self.state_generator._rte_grid
        return gridded_data, rte_grid

    def _make_gradient_dset(self, state, gradient):
        """
        Projects gridded gradient to 1D state vector gradient.

        Form an xr.Dataset for the gridded_gradient so the state_generator.project_gradient_to_state
        method can be used same as for the observational misfit.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The abstract 1D state vector
        gradient : np.ndarray, ndim=3
            The gradient in physical coordinates

        Returns
        -------
        out_gradient : np.ndarray, ndim=1
            The gradient for the abstract 1D state vector

        """
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

    def regularization_strength(self, iteration_number):
        return self._regularization_strength*self._relaxation_parameter**iteration_number

class WeightedRegularization(Regularization):
    """
    A generalization of `Regularization` methods that have spatially varying
    `regularization_strength` hyperparameters

    Parameters
    ----------
    state_generator : at3d.medium.StateGenerator
        This object contains the methods to map from the abstract state to the
        physical coordinates.
    scatterer_name : str
        The name of the scatterer that the regularization term is being applied to. This
        should match one of the keys in the solver.RTE.medium that is within the
        `state_generator`.
    variable_name : str
        The physical variable name that the regularization term is being applied to.
    regularization_strength : float
        A hyperparameter that scales the cost function of the regularization term.
    relaxation_parameter : float
        A hyperparameter controls the geometric sequence of `regularization_strength`
        as a function of optimization iteration number. If this is unity the
        effective regularization strength remains constant throughout the optimization.
        The closer this parameter is to zero the faster the regularization decreases.
    spatial_weights : float, np.ndarray, ndim=3
        These weights scale the `regularization_strength` at each grid point. No
        normalization is imposed so there is redundancy with the `regularization_strength`
        parameter. Both are used. If None then uniform weights are used.
    """
    def __init__(
        self, state_generator, scatterer_name, variable_name,
        regularization_strength, relaxation_parameter=1.0, spatial_weights=None):

        Regularization.__init__(self, state_generator, scatterer_name, variable_name,
        regularization_strength, relaxation_parameter=relaxation_parameter)

        if spatial_weights is None:
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
    """
    Adds a cost term based on the L1 norm of a physical variable thereby
    encouraging sparsity in that variable.

    Parameters
    ----------
    state_generator : at3d.medium.StateGenerator
        This object contains the methods to map from the abstract state to the
        physical coordinates.
    scatterer_name : str
        The name of the scatterer that the regularization term is being applied to. This
        should match one of the keys in the solver.RTE.medium that is within the
        `state_generator`.
    variable_name : str
        The physical variable name that the regularization term is being applied to.
    regularization_strength : float
        A hyperparameter that scales the cost function of the regularization term.
    relaxation_parameter : float
        A hyperparameter controls the geometric sequence of `regularization_strength`
        as a function of optimization iteration number. If this is unity the
        effective regularization strength remains constant throughout the optimization.
        The closer this parameter is to zero the faster the regularization decreases.
    spatial_weights : float, np.ndarray, ndim=3
        These weights scale the `regularization_strength` at each grid point. No
        normalization is imposed so there is redundancy with the `regularization_strength`
        parameter. Both are used. If None then uniform weights are used.
    """
    def __init__(
        self, state_generator, scatterer_name, variable_name,
        regularization_strength, relaxation_parameter=1.0, spatial_weights=None):
        WeightedRegularization.__init__(self, state_generator, scatterer_name, variable_name,
        regularization_strength, spatial_weights=spatial_weights, relaxation_parameter=relaxation_parameter)

    def __call__(self, state, iteration_number):

        gridded_data, rte_grid = self._get_gridded_variable(state)

        # should this just be a sum or should it be weighted by cell_volume
        # to account for non-uniform grid sizes?
        # If the latter then it will be easiest to add this as a case to the
        # fortran grid_smoothing code in util.f90.
        gradient = self.regularization_strength(iteration_number)*np.sign(gridded_data)
        cost = self.regularization_strength(iteration_number)*np.sum(np.abs(gridded_data))

        out_gradient = self._make_gradient_dset(state, gradient)

        return cost, out_gradient

class Tikhonov(WeightedRegularization):
    """
    Adds a cost term based on the L2 norm of a physical variable thereby encouraging
    small vectors, and therefore regularity in the presence of linear observational
    misfits.

    Parameters
    ----------
    state_generator : at3d.medium.StateGenerator
        This object contains the methods to map from the abstract state to the
        physical coordinates.
    scatterer_name : str
        The name of the scatterer that the regularization term is being applied to. This
        should match one of the keys in the solver.RTE.medium that is within the
        `state_generator`.
    variable_name : str
        The physical variable name that the regularization term is being applied to.
    regularization_strength : float
        A hyperparameter that scales the cost function of the regularization term.
    relaxation_parameter : float
        A hyperparameter controls the geometric sequence of `regularization_strength`
        as a function of optimization iteration number. If this is unity the
        effective regularization strength remains constant throughout the optimization.
        The closer this parameter is to zero the faster the regularization decreases.
    spatial_weights : float, np.ndarray, ndim=3
        These weights scale the `regularization_strength` at each grid point. No
        normalization is imposed so there is redundancy with the `regularization_strength`
        parameter. Both are used. If None then uniform weights are used.
    """
    def __init__(
        self, state_generator, scatterer_name, variable_name,
        regularization_strength, relaxation_parameter=1.0, spatial_weights=None):
        WeightedRegularization.__init__(self, state_generator, scatterer_name, variable_name,
        regularization_strength, spatial_weights=spatial_weights, relaxation_parameter=relaxation_parameter)

    def __call__(self, state, iteration_number):

        gridded_data, rte_grid = self._get_gridded_variable(state)

        # should this just be a sum or should it be weighted by cell_volume
        # to account for non-uniform grid sizes?
        # If the latter then it will be easiest to add this as a case to the
        # fortran grid_smoothing code in util.f90.
        gradient = self.regularization_strength(iteration_number)*2*gridded_data
        cost = self.regularization_strength(iteration_number)*np.sum(gridded_data**2)

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
    Adds a cost term based on the L1 or L2 norm of the spatial gradient of a
    physical variable, thereby encouraging smooth state vectors.

    Parameters
    ----------
    state_generator : at3d.medium.StateGenerator
        This object contains the methods to map from the abstract state to the
        physical coordinates.
    scatterer_name : str
        The name of the scatterer that the regularization term is being applied to. This
        should match one of the keys in the solver.RTE.medium that is within the
        `state_generator`.
    variable_name : str
        The physical variable name that the regularization term is being applied to.
    regularization_strength : float
        A hyperparameter that scales the cost function of the regularization term.
    relaxation_parameter : float
        A hyperparameter controls the geometric sequence of `regularization_strength`
        as a function of optimization iteration number. If this is unity the
        effective regularization strength remains constant throughout the optimization.
        The closer this parameter is to zero the faster the regularization decreases.
    spatial_weights : float, np.ndarray, ndim=3
        These weights scale the `regularization_strength` at each grid point. No
        normalization is imposed so there is redundancy with the `regularization_strength`
        parameter. Both are used. If None then uniform weights are used.
    mode : str
        Either 'l2' or 'l1', which chooses which type of norm to apply to the gradient vectors
    direction_weights : list
        A list of three floats corresponding to weights on the gradients in the x,y,z
        directions, respectively.
    """
    def __init__(
            self, state_generator, scatterer_name, variable_name,
            regularization_strength, relaxation_parameter=1.0, mode='l2', direction_weights=[1.0, 1.0, 1.0],
            spatial_weights=None, huber_parameter=1.0):
        WeightedRegularization.__init__(self, state_generator, scatterer_name, variable_name,
        regularization_strength, spatial_weights=spatial_weights, relaxation_parameter=relaxation_parameter)

        valid_modes = ('l1', 'l2', 'ph')
        if mode not in valid_modes:
            raise ValueError(
                "Spatial Smoothing `mode` should be in {}".format(valid_modes)
            )
        self._mode = mode
        self._huber_parameter = huber_parameter

        direction_weights = np.atleast_1d(direction_weights)
        self._direction_weights = direction_weights

    def __call__(self, state, iteration_number):

        gridded_data, rte_grid = self._get_gridded_variable(state)

#       The x,y,z conventions are transposed in this function compared
#       to what is used by at3d, that is why the names don't match up.
        cost, gradient, ierr, errmsg = at3d.core.grid_smoothing(
            field=gridded_data,
            weights=self._spatial_weights,
            direction_weights=self._direction_weights,
            nx=gridded_data.shape[2],
            ny=gridded_data.shape[1],
            nz=gridded_data.shape[0],
            zgrid=rte_grid.x.data,
            ygrid=rte_grid.y.data,
            xgrid=rte_grid.z.data,
            mode=self._mode,
            huber_parameter=self._huber_parameter
        )
        at3d.checks.check_errcode(ierr, errmsg)
        cost *= self.regularization_strength(iteration_number)
        gradient *= self.regularization_strength(iteration_number)

        # transform the gradient numpy array to an xarray dset that is then
        # projected to state.
        out_gradient = self._make_gradient_dset(state, gradient)
        return cost, out_gradient

class PriorCovariance:
    """
    the matrix and vector should be supplied in the coordinates
    of the abstract state.
    """
    def __init__(self, inverse_covariance, central_vector):

        self._inverse_covariance = inverse_covariance
        self._central_vector = central_vector

    def __call__(self, state):
        gradient = np.dot(state - self._central_vector, self._inverse_covariance)
        return gradient


def sigmoidal_spatial_weights(distance_to_clear, scaling_distance=0.1):
    """
    A useful function for computing spatial weights for regularization based
    on distance to the edge of the cloud at each grid point.
    """
    return 2*ss.expit(-1*distance_to_clear/scaling_distance)

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
