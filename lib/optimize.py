import scipy.optimize
import numpy as np
import time

import pyshdom.gradient

class ObjectiveFunction:
    """
    TODO
    """

    def __init__(self, measurements, loss_fn, min_bounds=None, max_bounds=None):
        self.measurements = measurements
        self.loss_fn = loss_fn
        self._bounds = list(zip(np.atleast_1d(min_bounds), np.atleast_1d(max_bounds)))
        self._loss = None

    def __call__(self, state):
        loss, gradient = self.loss_fn(state, self.measurements)
        self.loss = loss
        return loss, gradient

    @classmethod
    def LevisApproxUncorrelatedL2(cls, measurements, solvers, forward_sensors, unknown_scatterers, set_state_fn,
                                  project_gradient_to_state, n_jobs=1, mpi_comm=None, verbose=True,
                                  maxiter=100, init_solution=True, exact_single_scatter=True,
                                  min_bounds=None, max_bounds=None):
        """
        NB The passed set_state_fn must be defined using the
        solvers/unknown_scatterers defined at the script level.
        """
        def loss_function(state, measurements):

            set_state_fn(state)
            loss, gradient = pyshdom.gradient.levis_approx_uncorrelated_l2(measurements,
                                                                         solvers,
                                                                         forward_sensors,
                                                                         unknown_scatterers,
                                                                         n_jobs=n_jobs,
                                                                         mpi_comm=mpi_comm,
                                                                         verbose=verbose,
                                                                         maxiter=maxiter,
                                                                         init_solution=init_solution,
                                                                         exact_single_scatter=exact_single_scatter)
            state_gradient = project_gradient_to_state(state, gradient)
            return loss, state_gradient

        return cls(measurements, loss_function, min_bounds=min_bounds, max_bounds=max_bounds)

    @property
    def loss(self):
        return self._loss

    @property
    def bounds(self):
        return self._bounds

class PriorFunction:
    """
    TODO
    """
    def __init__(self, prior_fn, scale=1.0):
        self._prior_fn = prior_fn
        self._loss = None
        self._scale = scale

    def __call__(self, state):
        loss, gradient = self._prior_fn(state)
        self._loss = self._scale * loss
        return np.array(self._loss), self._scale * np.array(gradient)

    @classmethod
    def mahalanobis(cls, mean=0, cov=None, cov_inverse=None, scale=1.0e-4):
        """
        TODO
        """
        cov_inverse = np.linalg.inv(cov) if cov is not None else cov_inverse
        def mahalanobis_loss_fn(state):
            gradient = np.matmul(cov_inverse, state-mean)
            loss = np.matmul(state-mean, gradient)
            return np.array(loss), 2*np.array(gradient).ravel()
        return cls(prior_fn=mahalanobis_loss_fn, scale=scale)

    @property
    def loss(self):
        return self._loss

    @property
    def scale(self):
        return self._scale

class Optimizer:
    """
    Optmizer wrapps the scipy optimization methods.
    Notes
    -----
    For documentation:
        https://docs.scipy.org/doc/scipy/reference/optimize.html
    """
    def __init__(self,
                 objective_fn,
                 prior_fn=None,
                 callback_fn=None,
                 method='L-BFGS-B',
                 options={'maxiter': 100, 'maxls': 10, 'disp': True, 'gtol': 1e-16, 'ftol': 1e-8}
                 ):

        self._method = method
        self._options = options
        self._objective_fn = objective_fn
        self._prior_fn = np.atleast_1d(prior_fn) if prior_fn is not None else None
        self._callback_fn = np.atleast_1d(callback_fn)
        self._callback = None if callback_fn is None else self.callback
        self._iteration = None

    def callback(self, state): #TODO check whether the callback function below should call state.
        """
        The callback function invokes the callbacks defined by the writer (if any).
        Additionally it keeps track of the iteration number.
        """
        self._iteration += 1
        [function() for function in self._callback_fn]

    def objective(self, state):
        """
        The callback function invokes the callbacks defined by the writer (if any).
        Additionally it keeps track of the iteration number.
        """
        loss, gradient = self._objective_fn(state)
        if self._prior_fn is not None:
            p_loss, p_gradient = np.array([prior(state) for prior in self._prior_fn]).sum(axis=0)
            loss, gradient = loss + p_loss, gradient + p_gradient
        return loss, gradient

    def minimize(self, initial_state, iteration_step=0):
        """
        Local minimization with respect to the parameters defined.
        """
        self._iteration = iteration_step
        args = {
            'fun': self.objective,
            'x0': initial_state,
            'method': self._method,
            'jac': True,
            'options': self._options,
            'callback': self._callback
        }
        if self.method not in ['CG', 'Newton-CG']:
            if (len(self._objective_fn.bounds) == 1 ) & (self._objective_fn.bounds[0] == (None, None)):
                args['bounds'] = list(zip(np.repeat(np.atleast_1d(None), len(initial_state)),
                                          np.repeat(np.atleast_1d(None), len(initial_state))))
            elif len(self._objective_fn.bounds) == len(initial_state):
                args['bounds'] = self._objective_fn.bounds
            else:
                print('No bounds used.')
        result = scipy.optimize.minimize(**args)
        return result

    #TODO add save method

    @property
    def objective_fn(self):
        return self._objective_fn

    @property
    def iteration(self):
        return self._iteration

    @property
    def method(self):
        return self._method

    @property
    def options(self):
        return self._options


class CallbackFn:
    def __init__(self, callback_fn, ckpt_period=-1):
        self._ckpt_period = ckpt_period
        self._ckpt_time = time.time()
        self._callback_fn = callback_fn

    def __call__(self):
        time_passed = time.time() - self._ckpt_time
        if time_passed > self._ckpt_period:
            self._ckpt_time = time.time()
            self._callback_fn()
