"""
A module for storing regularization routines.
"""
import numpy as np

class SpatialGradientL2:
    """
    Calculates the L2 norm of the gradient vector of a given
    gridded unknown.

    Minimizing this norm enforces spatial smoothness of the unknown.
    This class is not fully generalized.
    """
    def __init__(self, grid, state_to_grid_projection,
                 grid_to_state_projection):
        """
        Parameters
        ----------
        grid : xr.Dataset
            A valid pyshdom grid object e.g. the SHDOM grid (see grid.py)
        state_to_grid_projection : callable
            Maps the state vector onto the grid. May only map a specific part
            corresponding to a single variable e.g. extinction.
        grid_to_state_projection : callable
            Maps the spatial gradients back to gradients in the
            state vector.
        """
        self._grid = grid
        self._state_to_grid_projection = state_to_grid_projection
        self._grid_to_state_projection = grid_to_state_projection

        self._x_matrix, self._x_integral_weights = self.get_derivative_matrix('x')
        self._y_matrix, self._y_integral_weights = self.get_derivative_matrix('y')
        self._z_matrix, self._z_integral_weights = self.get_derivative_matrix('z')

    def __call__(self, state):
        """
        Calculates the cost function and gradient vector for a given state vector.
        """

        data = self._state_to_grid_projection(state)
        gradient = np.zeros((3,)+ data.shape)
        cost = 0.0

        #z gradients
        for i in range(self._grid.x.size):
            for j in range(self._grid.y.size):
                gradient[2, i, j] = self._z_integral_weights*np.dot(
                    data[i, j, :], np.dot(self._z_matrix.T, self._z_matrix)
                    )
                cost += np.dot(gradient[2, i, j, :], data[i, j, :])
        #y gradients
        for i in range(self._grid.x.size):
            for j in range(self._grid.z.size):
                gradient[1, i, :, j] = self._y_integral_weights*np.dot(
                    data[i, :, j], np.dot(self._y_matrix.T, self._y_matrix)
                    )
                cost += np.dot(gradient[1, i, :, j], data[i, :, j])
        #x gradients
        for i in range(self._grid.y.size):
            for j in range(self._grid.z.size):
                gradient[0, :, i, j] = self._x_integral_weights*np.dot(
                    data[:, i, j], np.dot(self._x_matrix.T, self._x_matrix)
                    )
                cost += np.dot(gradient[0, :, i, j], data[:, i, j])

        gradient_total = np.sum(gradient, axis=0)
        gradient_for_state = self._grid_to_state_projection(state, gradient_total)
        return cost, gradient_for_state

    def get_derivative_matrix(self, coord):
        """
        Calculates the sensitivity of the gradient norm to each point
        for a given grid accounting for non-uniform grid spacing
        according to the formula used by np.gradient.

        Parameters
        ----------
        coord : string
            One of 'x', 'y', 'z'.
        Returns
        -------
        matrix : np.ndarray
            The matrix for calculating the gradient derivative
            for a column along the coordinate direction.
        integral_weights : The grid spacing that are used to weight the
            cost function.
        """
        diff = self._grid[coord].diff(coord).data
        hd = diff[1:]
        hs = diff[:-1]
        matrix = np.diag(np.append(0, hs**2/(hs*hd*(hs+hd))), k=1) +  \
                -1*np.diag(np.append(hd**2/(hs*hd*(hs+hd)), 0), k=-1) +\
                np.diag(np.append(np.append(0, (hd**2 - hs**2)/(hs*hd*(hs+hd))), 0))
        matrix[-1, -2:] = np.array([-2, 2]) / (2*diff[-1])
        matrix[0, :2] = np.array([2, -2]) / (2*diff[0])

        integral_weights = np.append(np.append(0.5*diff[0], 0.5*(diff[1:] + diff[:-1])),
                                     0.5*diff[-1])
        return matrix, integral_weights


class LocalCorrelation:

    def __init__(self, grid, state_on_grid_indicators, correlation_length, state_start=None, state_end=None):

        from scipy.spatial import distance_matrix
        #TODO properly check inputs.
        self._state_start, self._state_end = state_start, state_end

        a,b,c = state_on_grid_indicators
        xs = grid.x.data[a]
        ys = grid.y.data[b]
        zs = grid.z.data[c]
        positions = np.stack((xs,ys,zs),axis=1)
        dist_mat = distance_matrix(positions, positions)

        self._grid = grid
        self._dist_mat = dist_mat
        self._correl_matrix = np.exp(-dist_mat/correlation_length)

    def __call__(self, state):

        grad_vect = np.zeros(state.shape)
        if self._state_start is None and self._state_end is None:
            state_start = 0
            state_end = len(state)
        else:
            state_start, state_end = self._state_start, self._state_end
        temp = state[state_start:state_end]
        residual = temp - np.dot(self._correl_matrix, temp)
        cost = np.dot(residual, np.dot(np.diag(np.ones(temp.shape)), residual)) / (temp.size**2)
        grad = 2*np.dot(residual, np.diag(np.ones(temp.shape))) / temp.size
        grad_vect[state_start:state_end] = grad
        return cost, grad_vect
