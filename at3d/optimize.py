"""
Defines the high level abstract objects for handling solving
the inverse problem. 

This includes a generic `ObjectiveFunction` object
which specifies the data misfit with a specific case that uses the 
Levis Approximation to compute the gradient of the data misfit.

The `Optimizer` is a wrapper around scipy.optimize.minimize that
includes data misfit terms as well as priors from regularization.py,
for example.
"""

import time
import scipy.optimize
import numpy as np

import at3d.gradient

class ObjectiveFunction:
    """
    Specifies the type of data misfit term to use in the cost function
    in the optimization and concurrently the linearization of that term.

    Parameters
    ----------
    measurements : at3d.containers.SensorsDict
        The observations that constrain the optimization problem.
    loss_fn : callable
        The function which evaluates the data misfit term when supplied with
        a state and `measurements`
    min_bounds :

    max_bounds :
    """
    def __init__(self, measurements, loss_fn, min_bounds=None, max_bounds=None):
        self.measurements = measurements
        self.loss_fn = loss_fn
        self._bounds = list(zip(np.atleast_1d(min_bounds), np.atleast_1d(max_bounds)))
        self._loss = None
        self._total_obj_fn_time = 0.0
        self._ncalls = 0

    def __call__(self, state):
        time1 = time.process_time()
        loss, gradient = self.loss_fn(state, self.measurements)
        self._total_obj_fn_time += time.process_time() - time1
        self._loss = loss
        self._ncalls += 1
        return loss, gradient

    # @classmethod
    # def Flexible(cls, measurements, gradient_fn, set_state_fn, project_gradient_to_state_fn,
    #                 min_bounds=None, max_bounds=None):
    #
    #     def loss_function(state, measurements):
    #         set_state_fn(state)
    #         loss, gradient, _ = gradient_fn(measurements)
    #         state_gradient = project_gradient_to_state_fn(state, gradient)
    #         return loss, state_gradient
    #
    #     return cls(measurements, loss_function, min_bounds=min_bounds,
    #                 max_bounds=max_bounds)

    #
    @classmethod
    def LevisApproxUncorrelatedL2(cls, measurements, solvers, forward_sensors, unknown_scatterers, set_state_fn,
                                  project_gradient_to_state, parallel_solve_kwargs={'n_jobs': 1, 'mpi_comm':None,
                                  'verbose':True, 'maxiter':100, 'init_solution':True},
                                  gradient_kwargs={'cost_function': 'L2', 'exact_single_scatter':True},
                                  uncertainty_kwargs={'add_noise': False},
                                  min_bounds=None, max_bounds=None):
        """
        Use the Levis approximation to the linearization of an least squares cost function.
        Only error covariances between Stokes components for the same pixel are supported.

        This is just a thin wrapper around at3d.gradient.levis_approx_uncorrelated_l2
        that shows how to define an `ObjectiveFunction` object.

        Parameters
        ----------
        cls : at3d.optimize.ObjectiveFunction
            This method will generate a special case of this class.
        measurements : at3d.containers.SensorsDict
            Contains the observations that constrain the optimization problem.
        solvers : at3d.containers.SolversDict
            Contains the solver.RTE objects that will be used to evaluate the
            data misfit term.
        forward_sensors : at3d.containers.SensorsDict
            The object that will store the evaluation of the Forward model
            (solver.RTE) at each iteration. This should be consistent with
            `measurements`
        unknown_scatterers : at3d.containers.UnknownScatterers
            The object which defines which unknowns to calculate derivatives
            for that are then used in the evaluation of the state gradient
            (see project_gradient_to_state). Also holds the method and data
            for evaluating the microphysical derivatives with respect to the
            unknowns.
        set_state_fn : callable
            A function which takes the state as input and updates the `solvers`
            so that they reflect the value of the state vector and
            `solvers`.parallel_solve is ready to be called. This function
            is typically defined at the script level.
        project_gradient_to_state : callable
            A function which evaluates the gradient for the abstract state vector
            using the gridded gradient output
            (see at3d.gradient.levis_approx_uncorrelated_l2) and the state
            vector.
        min_bounds : TODO

        max_bounds : TODO

        Returns
        -------
        An instance of ObjectiveFunction.

        Notes
        -----
        The callables `set_state_fn` and `project_gradient_to_state`
        must be defined consistently, this is almost certainly the easiest
        source of error in the setup of an optimization.
        """
        gradient_fun = at3d.gradient.LevisApproxGradientUncorrelated(
            measurements, solvers, forward_sensors, unknown_scatterers, parallel_solve_kwargs,
            gradient_kwargs, uncertainty_kwargs)

        def loss_function(state, measurements):

            set_state_fn(state)
            loss, gradient, _ = gradient_fun()
            # loss, gradient = at3d.gradient.levis_approx_uncorrelated_l2(measurements,
            #                                                              solvers,
            #                                                              forward_sensors,
            #                                                              unknown_scatterers,
            #                                                              n_jobs=n_jobs,
            #                                                              mpi_comm=mpi_comm,
            #                                                              verbose=verbose,
            #                                                              maxiter=maxiter,
            #                                                              init_solution=init_solution,
            #                                                              exact_single_scatter=exact_single_scatter)
            state_gradient = project_gradient_to_state(state, gradient)
            return loss, state_gradient

        return cls(measurements, loss_function, min_bounds=min_bounds, max_bounds=max_bounds)

    @property
    def loss(self):
        return self._loss

    @property
    def bounds(self):
        return self._bounds

class Optimizer:
    """
    Optmizer wraps the scipy optimization methods.
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
        #self._callback = None if callback_fn is None else self.callback
        self._iteration = None
        self._state = None

    def callback(self, state): #TODO check whether the callback function below should call state.
        """
        The callback function invokes the callbacks defined by the writer (if any).
        Additionally it keeps track of the iteration number.
        """
        self._iteration += 1
        if self._callback_fn[0] is not None:
            [function(optimizer=self) for function in self._callback_fn]

    def objective(self, state):
        """
        The callback function invokes the callbacks defined by the writer (if any).
        Additionally it keeps track of the iteration number.
        """
        self._state = state
        loss, gradient = self._objective_fn(state)
        if self._prior_fn is not None:
            p_loss = []
            p_gradient = []
            for prior in self._prior_fn:
                ploss, pgrad = prior(state, self._iteration)
                p_loss.append(ploss)
                p_gradient.append(pgrad)
            loss += np.stack(p_loss, axis=0).sum(axis=0)
            gradient += np.stack(p_gradient, axis=0).sum(axis=0)

        return loss, gradient

    def minimize(self, initial_state, iteration_step=0, **kwargs):
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
            'callback': self.callback
        }
        args.update(kwargs)
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