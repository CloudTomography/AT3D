"""
TODO module docstring
"""
import warnings
import numpy as np

class Uncertainty:

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
        except np.linalg.LinAlgError as err:
            warnings.warn("`inverse_covariance` of errors in uninvertible."
                          "It cannot be used to generate noise and will "
                          "cause an exception to be raised if attempted.")
            self._covariance = None

    def calculate_uncertainties(self, sensor):
        #NB Be aware that repeated dims cause errors so the second dim is set to
        #'num_uncertainty2' even though they are identical.
        #Issue 1378 on xarray.
        sensor['uncertainties'] = (['num_uncertainty', 'num_uncertainty2', 'npixels'],
                np.repeat(self.inverse_covariance[..., np.newaxis],
                    sensor.npixels.size, axis=-1)
                )

    def add_noise(self, sensor):

        if self._covariance is None:
            raise ValueError("Noise cannot be generated as `inverse_covariance`"
             " was not invertible.")

        perturbations = np.random.multivariate_normal(
            mean=np.zeros(self.num_uncertainty), cov=self._covariance, size=(sensor.sizes['npixels'])
            )
        if 'I' in sensor.data_vars:
            sensor['I'][:] += perturbations[0]
        else:
            raise KeyError("Radiance 'I' is not found in sensor."
            " Cannot generate noise on measurements that do not exist.")
        if 'Q' in sensor.data_vars:
            sensor['Q'][:] += perturbations[1]
        if 'U' in sensor.data_vars:
            sensor['U'][:] += perturbations[2]

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
