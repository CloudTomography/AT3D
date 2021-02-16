"""
This module contains objects which form and describe the inverse error
covariance matrix used in the inverse problem for weighting data misfits.
The choice of uncertainty is closely tied to the choice of cost function
(ie what misfits are being calculated.). See UPDATE_COSTFUNCTION in shdomsub4.f
for implemented cost functions and how theuncertainty values are used.
Uncertainties are currently independent per pixel but error correlations
between Stokes quantities are allowed for a given pixel.

New uncertainty objects should have a `calculate_uncertainties` and `add_noise`
method and inherit from the base class `Uncertainty`. If new cost functions are
implemented they should be added as valid cases to `Uncertainty`.
"""
import warnings
import numpy as np

import pyshdom.checks

class Uncertainty:
    """
    A base class for uncertainty models. A simple additive gaussian noise model
    which sets the template for these objects.

    Parameters
    ----------
    inverse_covariance : np.ndarray
        The inverse covariance matrix of errors between stokes components (I, Q, U, V).
        Should be square (ndim=2) and should have a shape consistent with
        `cost_function`.
    cost_function : str
        The two character string that picks the cost function. This is passed to
        solver.RTE.levisapprox_gradient() and used in UPDATE_COSTFUNCTION
        in src/shdomsub4.f. Valid values are hardcoded to be restricted here to
        ('L2', 'LL').

    Raises
    ------
    NotImplementedError
        If an invalid `cost_function` is supplied.
    ValueError
        If the shape of inverse_covariance is not consistent with the cost
        function.
    """
    def __init__(self, inverse_covariance, cost_function):


        self._valid_cost_functions = ('L2', 'LL')
        if not cost_function in self._valid_cost_functions:
            raise NotImplementedError("`cost_function` '{}' is not supported for "
                "this Uncertainty type. Valid values are '{}'".format(
                cost_function, self._valid_cost_functions)
                )
        self._cost_function = cost_function
        if self.cost_function == 'L2':
            self._num_uncertainty = 4
        elif self.cost_function == 'LL':
            self._num_uncertainty = 2

        if not (inverse_covariance.shape[0] == inverse_covariance.shape[1]) == self._num_uncertainty:
            raise ValueError(
                "`inverse_covariance` should be of shape=({0}, {0}) for `cost_function`"
                " '{1}''".format(
                    self._num_uncertainty, self.cost_function
                )
            )
        self._inverse_covariance = inverse_covariance
        try:
            self._covariance = np.linalg.inv(self.inverse_covariance)
        except np.linalg.LinAlgError:
            warnings.warn("`inverse_covariance` of errors in uninvertible."
                          "It cannot be used to generate noise and will "
                          "cause an exception to be raised if attempted.")
            self._covariance = None

    def calculate_uncertainties(self, sensor):
        """Replicates the inverse error-covariance matrix for all pixels and
        updates the sensor.

        Parameters
        ----------
        sensor : xr.Dataset
            A valid sensor containing Stokes components at specific locations
            and directions. See sensor.py for details.
        """
        pyshdom.checks.check_sensor(sensor)
        #NB Be aware that repeated dims cause errors so the second dim is set to
        #'num_uncertainty2' even though they are identical.
        #Issue 1378 on xarray.
        sensor['uncertainties'] = (['num_uncertainty', 'num_uncertainty2', 'npixels'],
                np.repeat(self.inverse_covariance[..., np.newaxis],
                    sensor.npixels.size, axis=-1)
                )

    def add_noise(self, sensor):
        """Calculates additive perturbations to the sensor's Stokes Components
        that are drawn from the error co-variance matrix.

        Parameters
        ----------
        sensor : xr.Dataset
            A valid sensor containing Stokes components at specific locations
            and directions. See sensor.py for details.

        Raises
        ------
        ValueError
            If there is no valid error-covariance matrix.
        KeyError
            If some of the required observables are not available.

        """
        if self._covariance is None:
            raise ValueError("Noise cannot be generated as `inverse_covariance`"
             " was not invertible.")

        perturbations = np.random.multivariate_normal(
            mean=np.zeros(self.num_uncertainty), cov=self._covariance, size=(sensor.sizes['npixels'])
            )
        pyshdom.checks.check_sensor(sensor)
        for i, has_stokes in enumerate(sensor.stokes.data):
            if has_stokes and (not sensor.stokes_index[i].data in sensor.data_vars):
                raise KeyError(
                    "Stokes component '{}' is not found in sensor even though it is an "
                    "observable. Noise perturbations for this observable cannot be generated.".format(
                    sensor.stokes_index[i].data
                    ))
            sensor[str(sensor.stokes_index[i].data)][:] += perturbations[i]

    @property
    def inverse_covariance(self):
        return self._inverse_covariance
    @property
    def cost_function(self):
        return self._cost_function
    @property
    def num_uncertainty(self):
        return self._num_uncertainty
    @property
    def _valid_cost_functions(self):
        return self._valid_cost_functions

class NullUncertainty(Uncertainty):
    """
    An equal weighting uncertainty model that is applied if none is supplied
    to ensure that a valid inverse error co-variance matrix is passed to
    solver.RTE.levisapprox_gradient.

    A diagonal matrix of the correct shape is generated to match the cost function
    that is going to be used.

    Parameters
    ----------
    cost_function : str
        The cost function that the gradient is set to use.
    """
    def __init__(self, cost_function):
        if cost_function == 'L2':
            inverse_covariance = np.diagonal(np.ones(4))
        elif cost_function == 'LL':
            inverse_covariance = np.diagonal(np.ones(2))
        super().__init__(inverse_covariance, cost_function)

    def add_noise(self, sensor):
        raise ValueError(
            "{} cannot be used to generate measurement noise. Please assign "
            "another uncertainty model.".format(type(self))
        )

class DoLPUncertainty(Uncertainty):
    """
    An unfinished Uncertainty class for the 'LL' cost function which uses
    DoLP as an observable rather than Q and U.
    """
    def __init__(self, dolp_abs_uncertainty, radiance_relative_uncertainty):
        inverse_covariance = np.diag(
            np.array([1.0/radiance_relative_uncertainty, 1.0/dolp_abs_uncertainty])
            )
        super().__init__(inverse_covariance, 'LL')

    def add_noise(self, sensor):
        raise NotImplementedError(
            "{} cannot yet be used to generate measurement noise. Please assign "
            "another uncertainty model or finish implementing this funciton.".format(type(self))
        )
        rad_perturbations = self._covariance[0,0]*np.random.normal(size=(sensor.sizes['npixels']))
        sensor['I'][:] = np.exp(np.log(sensor['I'].data) + rad_perturbations)

        dolp_perturbations = self._covariance[1,1]*np.random.normal(size=(sensor.sizes['npixels']))
        #need to assign perturbations to U/Q that produce the desired perturbations in DoLP
        #given the perturbations to radiance.
