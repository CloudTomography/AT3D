"""
This module contains objects for handling coordinate transforms for
the inverse problem.
While this can be handled at the script level, these objects cover
some useful cases and are integrated into the streamlined
inverse framework using at3d.medium.StateGenerator rather
than user specified 'set state functions'.

These objects apply the transformations from the abstract state level
to the spatial grid, any transformations of the state vector.
"""
import numpy as np

class CoordinateTransform:
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
    def __call__(self, abstract_state):
        """
        Transform from abstract state to physical coordinates.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.

        Returns
        -------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.
        """
        physical_state = abstract_state
        return physical_state

    def inverse_transform(self, physical_state):
        """
        Transform from physical coordinates to abstract coordinates.

        Parameters
        ----------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.

        Returns
        -------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        """
        abstract_state = physical_state
        return abstract_state

    def gradient_transform(self, abstract_state, physical_gradient):
        """
        Transform gradient from physical coordinates to abstract coordinates
        applying the chain rule.

        Parameters
        ----------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        physical_gradient : np.ndarray, ndim=1
            The 1D gradient vector in physical coordinates

        Returns
        -------
        abstract_gradient : np.ndarray, ndim=1
            The 1D gradient vector in abstract coordinates.
        """
        abstract_gradient = physical_gradient
        return abstract_gradient

class CoordinateTransformLog(CoordinateTransform):
    """
    State vector is log(physical_coordinates)
    """
    def __call__(self, abstract_state):
        """
        Transform from abstract state to physical coordinates.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.

        Returns
        -------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.
        """
        physical_state = np.exp(abstract_state)
        return physical_state

    def inverse_transform(self, physical_state):
        """
        Transform from physical coordinates to abstract coordinates.

        Parameters
        ----------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.

        Returns
        -------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        """
        abstract_state = np.log(physical_state)
        return abstract_state

    def gradient_transform(self, abstract_state, physical_gradient):
        """
        Transform gradient from physical coordinates to abstract coordinates
        applying the chain rule.

        Parameters
        ----------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        physical_gradient : np.ndarray, ndim=1
            The 1D gradient vector in physical coordinates

        Returns
        -------
        abstract_gradient : np.ndarray, ndim=1
            The 1D gradient vector in abstract coordinates.
        """
        abstract_gradient = physical_gradient*np.exp(abstract_state)
        return abstract_gradient

class CoordinateTransformScaling(CoordinateTransform):
    """
    A linear scaling for all unknowns of a single variable.

    This has little benefit if applied globally, but can be useful
    in multi-variable optimization if applied separately to each variable.

    State coordinates are given by (physical_coordinates - offset)*scaling_factor.

    Parameters
    ----------
    offset : float
        Constant offset parameter in linear transform.
    scaling : float
        Scaling factor in linear transform.
    """
    def __init__(self, offset, scaling_factor):

        self._offset = offset
        self._scaling_factor = scaling_factor

    def __call__(self, abstract_state):
        """
        Transform from abstract state to physical coordinates.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.

        Returns
        -------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.
        """
        physical_state = abstract_state/self._scaling_factor + self._offset
        return physical_state

    def inverse_transform(self, physical_state):
        """
        Transform from physical coordinates to abstract coordinates.

        Parameters
        ----------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.

        Returns
        -------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        """
        abstract_state = (physical_state - self._offset)*self._scaling_factor
        return abstract_state

    def gradient_transform(self, abstract_state, physical_gradient):
        """
        Transform gradient from physical coordinates to abstract coordinates
        applying the chain rule.

        Parameters
        ----------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        physical_gradient : np.ndarray, ndim=1
            The 1D gradient vector in physical coordinates

        Returns
        -------
        abstract_gradient : np.ndarray, ndim=1
            The 1D gradient vector in abstract coordinates.
        """
        abstract_gradient = self.inverse_transform(physical_gradient)
        return abstract_gradient

class CoordinateTransformExp(CoordinateTransform):
    """
    An exponential scaling for all variables.

    This maps unknowns into the interval [0, 1] and compresses the large
    values into a compressed interval.

    Parameters
    ----------
    scaling : float
        The e-folding length used in the exponential transform.
    """
    def __init__(self, scaling):

        self._scaling = scaling

    def __call__(self, abstract_state):
        """
        Transform from abstract state to physical coordinates.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.

        Returns
        -------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.
        """
        physical_state = -np.log(1.0-abstract_state)*self._scaling
        return -physical_state

    def inverse_transform(self, physical_state):
        """
        Transform from physical coordinates to abstract coordinates.

        Parameters
        ----------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.

        Returns
        -------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        """
        abstract_state = 1.0 - np.exp(-physical_state/self._scaling)
        return abstract_state

    def gradient_transform(self, abstract_state, physical_gradient):
        """
        Transform gradient from physical coordinates to abstract coordinates
        applying the chain rule.

        Parameters
        ----------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        physical_gradient : np.ndarray, ndim=1
            The 1D gradient vector in physical coordinates

        Returns
        -------
        abstract_gradient : np.ndarray, ndim=1
            The 1D gradient vector in abstract coordinates.
        """
        abstract_gradient = -self._scaling*physical_gradient/(abstract_state-1.0)
        return abstract_gradient

class CoordinateTransformHyperBol(CoordinateTransform):
    """
    A hyperbolic transform.

    This maps unknowns into the interval [0, 1] and compresses the large
    values into a smaller interval. This transform somewhat linearizes
    the dependence of radiance on extinction by accounting for the
    reduced sensitivity to extinction as extinction increases.

    Parameters
    ----------
    scaling : float
        The e-folding length used in the exponential transform.
    """
    def __init__(self, scaling):

        self._scaling = scaling

    def __call__(self, abstract_state):
        """
        Transform from abstract state to physical coordinates.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.

        Returns
        -------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.
        """
        physical_state = (abstract_state/(1.0-abstract_state))/self._scaling
        return physical_state

    def inverse_transform(self, physical_state):
        """
        Transform from physical coordinates to abstract coordinates.

        Parameters
        ----------
        physical_state : np.ndarray, ndim=1
            The 1D state vector in physical coordinates.

        Returns
        -------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        """
        abstract_state = self._scaling*physical_state/(1.0 + self._scaling*physical_state)
        return abstract_state

    def gradient_transform(self, abstract_state, physical_gradient):
        """
        Transform gradient from physical coordinates to abstract coordinates
        applying the chain rule.

        Parameters
        ----------
        abstract_state : np.ndarray, ndim=1
            The 1D state vector in abstract coordinates.
        physical_gradient : np.ndarray, ndim=1
            The 1D gradient vector in physical coordinates

        Returns
        -------
        abstract_gradient : np.ndarray, ndim=1
            The 1D gradient vector in abstract coordinates.
        """
        abstract_gradient = physical_gradient/(self._scaling* ((1.0 - abstract_state)**2))
        return abstract_gradient

class StateToGridMask:
    """
    Transforms gridded unknowns to and from a (possibly reduced) set of
    unknowns in a 1D state vector based on a grid-point mask.

    This covers the default Null case of simply reshaping into a 1D vector.

    Parameters
    ----------
    grid_shape : tuple, optional
        The shape of the gridded data.
    mask : np.ndarray, optional
        An array that indicates cloud at each element when tested to be True.
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
            if any(grid_shape != mask.shape):
                raise ValueError(
                    "Both `grid_shape` and `mask` arguments were provided."
                    " The shape of `mask` is not consistent with `grid_shape`."
                )
        self._mask = mask
        self._grid_shape = grid_shape

    def __call__(self, state):
        """
        Map from 1D state vector to gridded data using a volumetric mask.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The 1D state vector.

        Returns
        -------
        gridded_state : np.ndarray, ndim=3
            The gridded data of which a subset are unknowns that form the state vector.
        """
        gridded_state = np.zeros(self._grid_shape)
        gridded_state[np.where(self._mask)] = state
        return gridded_state

    def inverse_transform(self, gridded_data):
        """
        Map from gridded data to 1D state vector using a volumetric mask.

        Parameters
        ----------
        gridded_state : np.ndarray, ndim=3
            The gridded data of which a subset are unknowns that form the state vector.

        Returns
        -------
        state : np.ndarray, ndim=1
            The 1D state vector.
        """
        state = gridded_data[np.where(self._mask)]
        return state

    def gradient_transform(self, gridded_gradient):
        """
        Projects gradients of gridded quantities onto the 1D state vector.

        Parameters
        ----------
        gridded_gradient : np.ndarray, ndim=3
            The gridded gradient values of which a subset correspond to entries
            in the state vector.

        Returns
        -------
        gradient_for_state : np.ndarray, ndim=1
            The gradient vector corresponding to the 1D state vector.
        """
        gradient_for_state = gridded_gradient[np.where(self._mask)]
        return gradient_for_state

    def inverse_bounds_transform(self, gridded_bounds):
        """
        Projects bounds of gridded quantities onto the 1D state vector.

        Parameters
        ----------
        gridded_bounds : np.ndarray, ndim=3
            The gridded gradient values of which a subset correspond to entries
            in the state vector. May be either upper or lower bounds.

        Returns
        -------
        state_bounds : np.ndarray, ndim=1
            The 1D vector of bounds corresponding to the 1D state vector.
        """
        state_bounds = gridded_bounds[np.where(self._mask)]
        return state_bounds

class StateToGridProfile(StateToGridMask):
    """
    Maps gridded data to and from a 1D vector based on horizontal averaging
    conditioned on a 3D mask.

    This transform gives one unknown per vertical level.
    """
    def __call__(self, state):
        """
        Map from 1D state vector to gridded data using a volumetric mask.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The 1D state vector.

        Returns
        -------
        gridded_state : np.ndarray, ndim=3
            The gridded data of which a subset are unknowns that form the state vector.
        """
        gridded_state = np.zeros(self._grid_shape)*np.nan
        for i in range(self._mask.shape[-1]):
            gridded_state[np.where(self._mask[..., i]), i] = state[i]
        return gridded_state

    def inverse_transform(self, gridded_data):
        """
        Map from gridded data to 1D state vector using a volumetric mask.

        Parameters
        ----------
        gridded_state : np.ndarray, ndim=3
            The gridded data of which a subset are unknowns that form the state vector.

        Returns
        -------
        state : np.ndarray, ndim=1
            The 1D state vector.
        """
        gridded_state = np.zeros(self._grid_shape)*np.nan
        gridded_state[np.where(self._mask)] = gridded_data[self._mask]
        return np.nanmean(gridded_state, axis=(0, 1))

    def gradient_transform(self, gridded_gradient):
        """
        Projects gradients of gridded quantities onto the 1D state vector.

        Parameters
        ----------
        gridded_gradient : np.ndarray, ndim=3
            The gridded gradient values of which a subset correspond to entries
            in the state vector.

        Returns
        -------
        gradient_for_state : np.ndarray, ndim=1
            The gradient vector corresponding to the 1D state vector.
        """
        gradient_for_state = self.inverse_transform(gridded_gradient)
        return gradient_for_state

    def inverse_bounds_transform(self, gridded_bounds):
        """
        Projects bounds of gridded quantities onto the 1D state vector.

        Parameters
        ----------
        gridded_bounds : np.ndarray, ndim=3
            The gridded gradient values of which a subset correspond to entries
            in the state vector. May be either upper or lower bounds.

        Returns
        -------
        state_bounds : np.ndarray, ndim=1
            The 1D vector of bounds corresponding to the 1D state vector.
        """
        # strictly this isn't true if the bounds are non-uniform
        # in space. But in the simplest case that that is true
        # then this will work.
        if np.size(np.unique(gridded_bounds)) == 1:
            state_bounds = self.inverse_transform(gridded_bounds)
        else:
            raise NotImplementedError(
                "Inverse Transform for non-uniform bounds for single variable"
                " have not yet been implemented."
            )
        return state_bounds


class StateToGrid2D(StateToGridMask):
    """
    Maps gridded data to and from a 1D vector based on horizontal averaging
    conditioned on a 3D mask.

    This transform gives one unknown per vertical column.
    """
    def __call__(self, state):
        """
        Map from 1D state vector to gridded data using a volumetric mask.

        Parameters
        ----------
        state : np.ndarray, ndim=1
            The 1D state vector.

        Returns
        -------
        gridded_state : np.ndarray, ndim=3
            The gridded data of which a subset are unknowns that form the state vector.
        """
        gridded_state = np.zeros(self._grid_shape)*np.nan
        for i in range(self._mask.shape[0]):
            for j in range(self._mask.shape[1]):
                gridded_state[i, j, np.where(self._mask[i, j, :])] = state.reshape(self._mask.shape[:2])[i, j]
        return gridded_state

    def inverse_transform(self, gridded_data):
        """
        Map from gridded data to 1D state vector using a volumetric mask.

        Parameters
        ----------
        gridded_state : np.ndarray, ndim=3
            The gridded data of which a subset are unknowns that form the state vector.

        Returns
        -------
        state : np.ndarray, ndim=1
            The 1D state vector.
        """
        gridded_state = np.zeros(self._grid_shape)*np.nan
        gridded_state[np.where(self._mask)] = gridded_data[np.where(self._mask)]
        state = np.nanmean(gridded_state, axis=-1).ravel()
        return state

    def gradient_transform(self, gridded_gradient):
        """
        Projects gradients of gridded quantities onto the 1D state vector.

        Parameters
        ----------
        gridded_gradient : np.ndarray, ndim=3
            The gridded gradient values of which a subset correspond to entries
            in the state vector.

        Returns
        -------
        gradient_for_state : np.ndarray, ndim=1
            The gradient vector corresponding to the 1D state vector.
        """
        gradient_for_state = self.inverse_transform(gridded_gradient)
        return gradient_for_state

    def inverse_bounds_transform(self, gridded_bounds):
        """
        Projects bounds of gridded quantities onto the 1D state vector.

        Parameters
        ----------
        gridded_bounds : np.ndarray, ndim=3
            The gridded gradient values of which a subset correspond to entries
            in the state vector. May be either upper or lower bounds.

        Returns
        -------
        state_bounds : np.ndarray, ndim=1
            The 1D vector of bounds corresponding to the 1D state vector.
        """
        # strictly this isn't true if the bounds are non-uniform
        # in space. But in the simplest case that that is true
        # then this will work.
        if np.size(np.unique(gridded_bounds)) == 1:
            state_bounds = self.inverse_transform(gridded_bounds)
        else:
            raise NotImplementedError(
                "Inverse Transform for non-uniform bounds for single variable"
                " have not yet been implemented."
            )
        return state_bounds
