"""
This module contains objects which form and describe the inverse error
covariance matrix used in the inverse problem for weighting data misfits.
The choice of uncertainty is closely tied to the choice of cost function
(ie what misfits are being calculated.). See UPDATE_COSTFUNCTION in shdomsub4.f
for implemented cost functions and how theuncertainty values are used.
Uncertainties are currently independent per pixel but error correlations
between Stokes quantities are allowed for a given pixel.

New uncertainty objects should have `_process_uncertainties` and `_process_noise`
methods which are called and inherit from the base class `Uncertainty`. If new cost functions are
implemented they should be added as valid cases to `Uncertainty`.
"""
import warnings
import numpy as np
import at3d.checks

import scipy.interpolate as si
import scipy.optimize as so
import numpy as np
import at3d.checks

class Uncertainty:
    """
    A base class for uncertainty and noise calculations that
    has the public methods for interfacing with valid at3d
    sensor dataset objects.

    Parameters
    ----------
    cost_function : str
        The type of cost function. Supported values are 'L2' for
        generalized least squares on up to 4 Stokes components
        and 'LL' which is a least squares in the logarithm of the
        Intensity and DOLP.
    """

    def __init__(self, cost_function):

        self._valid_cost_functions = ('L2', 'LL')
        if not cost_function in self._valid_cost_functions:
            raise NotImplementedError("`cost_function` '{}' is not supported for "
                "this Uncertainty type. Valid values are '{}'".format(
                cost_function, self._valid_cost_functions)
                )
        self._cost_function = cost_function
        if self._cost_function == 'L2':
            self._num_uncertainty = 4
        elif self._cost_function == 'LL':
            self._num_uncertainty = 2

    def calculate_uncertainties(self, sensor):
        """
        Updates the sensor in-place with the uncertainties calculated
        according to the uncertainty model.

        Parameters
        ----------
        sensor : xr.Dataset
            A valid sensor containing Stokes components at specific locations
            and directions. See sensor.py for details.
        """
        at3d.checks.check_sensor(sensor)
        big_uncertainties = self._process_uncertainties(sensor)
        #NB Be aware that repeated dims cause errors so the second dim is set to
        #'num_uncertainty2' even though they are identical.
        #Issue 1378 on xarray.
        sensor['uncertainties'] = (
            ['num_uncertainty', 'num_uncertainty2', 'npixels'], big_uncertainties
        )

    def add_noise(self, sensor, noise=True):
        """
        Calculates additive perturbations to the sensor's Stokes Components
        that are drawn according to the error model.

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
        perturbed_stokes = self._process_noise(sensor)
        at3d.checks.check_sensor(sensor)
        for i, has_stokes in enumerate(sensor.stokes.data):
            if has_stokes:
                if (not sensor.stokes_index.data[i] in sensor.data_vars):
                    raise KeyError(
                        "Stokes component '{}' is not found in sensor even though it is an "
                        "observable. Noise perturbations for this observable cannot be generated.".format(
                        sensor.stokes_index[i].data
                        ))
                sensor[str(sensor.stokes_index[i].data)][:] = perturbed_stokes[i]

    @property
    def cost_function(self):
        return self._cost_function
    @property
    def num_uncertainty(self):
        return self._num_uncertainty
    @property
    def valid_cost_functions(self):
        return self._valid_cost_functions

class NullUncertainty(Uncertainty):
    """
    A null uncertainty model for unweighted least squares.

    This uncertainty model includes a scaling factor that can be used to
    scale the size of the residuals and therefore the optimization problem.

    Parameters
    ----------
    cost_function : str
        The type of cost function. Supported values are 'L2' for
        generalized least squares on up to 4 Stokes components
        and 'LL' which is a least squares in the logarithm of the
        Intensity and DOLP.
    scaling_factor : float
        Multiplicatively increases the uncertainties to thereby reduce the size
        of the cost function and its gradient.
    """
    def __init__(self, cost_function, scaling_factor=1.0):
        super().__init__(cost_function)
        self._scaling_factor = scaling_factor

    def _process_uncertainties(self, sensor):

        big_uncertainties = self._scaling_factor*np.ones((self._num_uncertainty, self._num_uncertainty, sensor.npixels.size))
        return big_uncertainties

    def add_noise(self, sensor):
        raise ValueError(
            "{} cannot be used to generate measurement noise. Please assign "
            "another uncertainty model.".format(type(self))
        )


class RadiometricUncertainty(Uncertainty):
    """
    A simple radiometric noise model that includes a flat error variance
    and a variance proportional to the radiance.

    All noise is modeled as poisson.

    There are also calibration uncertainties that can be included.
    The total radiometric uncertainty is the sum in quadrature of the
    `absolute_calibration_uncertainty` and the
    `camera_to_camera_calibration_uncertainty`.

    Parameters
    ----------
    cost_function : str
        The type of cost function. Supported values are 'L2' for
        generalized least squares on up to 4 Stokes components
        and 'LL' which is a least squares in the logarithm of the
        Intensity and DOLP.
    sigma_floor : float
        The standard deviation of the noise when there is no signal.
    shot_parameter : float
        The coefficient that produces the error variance that is
        proportional to the signal.
    absolute_calibration_uncertainty : float
        An optional calibration uncertainty. The perturbation is globally
        generated, i.e. there is one perturbation per instance of this
        class. This is a fractional multiplication perturbation.
        Radiances are multiplied by a random number drawn from a normal
        distribution with a mean of 1.0 and a standard deviation of
        `absolute_calibration_uncertainty`.
    camera_to_camera_calibration_uncertainty : float
        The uncertainty in the radiometric scale from camera to camera.
        This perturbation is unique per `sensor`.
    """
    def __init__(self, cost_function, snr_func, sigma_floor,
                 absolute_calibration_uncertainty=0.0,
                 camera_to_camera_calibration_uncertainty=0.0,
                seed=None):

        super().__init__(cost_function)

        self._minimum_noise = sigma_floor
        self._snr_func = snr_func
        self._camera_to_camera_calibration_uncertainty = camera_to_camera_calibration_uncertainty

        self._absolute_calibration_uncertainty = absolute_calibration_uncertainty
        if not seed is None:
            np.random.seed(seed)
        self._absolute_calibration_perturbation = np.random.normal(
            loc=0.0, scale=absolute_calibration_uncertainty
        )

    def noise_curve(self, reflectance):
        """
        Generate the estimated radiometric noise for a given reflectance.

        Parameters
        ----------
        reflectance : np.ndarray
            The same units as the radiance / reflectance unit output from SHDOM.
            Assuming solarflux is unity this is (RADIANCE / SOLARMU*FSUN)

        Returns
        -------
        noise_out : np.ndarray
            The SNR of the noise at the specified reflectance values.
        """
        reflectance = np.atleast_1d(reflectance)
        SNR_out = self._snr_func(reflectance)

        return SNR_out

    def _process_noise(self, sensor, seed=None, camera_cal=True, noise=True, absolute_cal=True):
        """
        Perturbs the intensity measurement according to the
        radiometric noise model.

        Parameters
        ----------
        radiance : np.ndarray
            An array of radiance data as output by SHDOM to be perturbed.
        seed : int
            An optional seed for a random number generator for generating the
            noise perturbations.
        camera_cal : bool
            Whether to add perturbations due to camera-to-camera calibration uncertainty.
        noise : bool
            Whether to add perturbations to the sensor due to radiometric noise.
        absolute_cal : bool
            Whether to add perturbations due to absolute calibration uncertainty.

        Returns
        -------
        radiance : np.ndarray
            The same array modified in-place with error perturbations.
        """
        at3d.checks.check_sensor(sensor)
        radiance = sensor.I.data
        if seed is not None:
            np.random.seed(seed)
        if noise:
            signal_mean = radiance
            snr = self.noise_curve(radiance)
            noise_levels = np.maximum(radiance/snr, self._minimum_noise)
            noise_levels[np.where(snr == 0.0)] = self._minimum_noise
            signal_var = noise_levels**2

            poisson_rate = (signal_mean**2)/signal_var
            samples = np.random.poisson(poisson_rate)
            radiance = np.sqrt(samples*signal_var)

        if camera_cal:
            radiance *= np.random.normal(loc=1.0, scale=self._camera_to_camera_calibration_uncertainty)

        if absolute_cal:
            radiance *= (1.0+self._absolute_calibration_perturbation)

        big_noise = np.zeros((4, radiance.size))
        big_noise[0] = radiance
        return big_noise

    def _process_uncertainties(self, sensor, camera_cal=True, noise=True, absolute_cal=True):
        """
        Generates uncertainties for each radiance according to the
        radiometric error model.

        Parameters
        ----------
        radiance : np.ndarray
            An array of radiance data as output by SHDOM.
        seed : int
            An optional seed for a random number generator for generating the
            noise perturbations.
        camera_cal : bool
            Whether to add perturbations due to camera-to-camera calibration uncertainty.
        noise : bool
            Whether to add perturbations to the sensor due to radiometric noise.
        absolute_cal : bool
            Whether to add perturbations due to absolute calibration uncertainty.

        Returns
        -------
        big_uncertainties : np.ndarray
            The array of uncertainties. Only the intensity components are non-unity.
        """
        at3d.checks.check_sensor(sensor)
        radiance = sensor.I.data

        errors = []
        if noise:
            snr = self.noise_curve(radiance)
            noise_levels = np.maximum(radiance/snr, self._minimum_noise)
            noise_levels[np.where(snr==0.0)] = self._minimum_noise
            signal_var = noise_levels**2

            errors.append(signal_var)
        if camera_cal:
            errors.append((self._camera_to_camera_calibration_uncertainty*radiance)**2)
        if absolute_cal:
            errors.append((self._absolute_calibration_perturbation*radiance)**2)

        if errors:
            uncertainties = 1.0/sum(errors)
        else:
            uncertainties = np.ones(radiance.shape)
        big_uncertainties = np.ones((self._num_uncertainty, self._num_uncertainty, uncertainties.size))
        big_uncertainties[0, 0] = uncertainties

        return big_uncertainties


class TabulatedRadiometricUncertainty(RadiometricUncertainty):
    """
    A simple radiometric error model that is initialized
    from a Look Up Table (LUT) of prescribed SNR values at different
    radiance values.

    The method uses a flat extrapolation of the noise standard deviation
    below values in the LUT and fits a square variation of SNR
    with reflectance to extrapolate to values larger than in the LUT.

    There are also calibration uncertainties that can be included.
    The total radiometric uncertainty is the sum in quadrature of the
    `absolute_calibration_uncertainty` and the
    `camera_to_camera_calibration_uncertainty`.

    Parameters
    ----------
    cost_function : str
        The type of cost function. Supported values are 'L2' for
        generalized least squares on up to 4 Stokes components
        and 'LL' which is a least squares in the logarithm of the
        Intensity and DOLP.
    reflectance_values : np.ndarray
        The reflectance values at which the `SNR_values` are known.
    SNR_values : np.ndarray
        Defined as `reflectance_values` / standard_devation_of_reflectance
    absolute_calibration_uncertainty : float
        An optional calibration uncertainty. The perturbation is globally
        generated, i.e. there is one perturbation per instance of this
        class. This is a fractional multiplication perturbation.
        Radiances are multiplied by a random number drawn from a normal
        distribution with a mean of 1.0 and a standard deviation of
        `absolute_calibration_uncertainty`.
    camera_to_camera_calibration_uncertainty : float
        The uncertainty in the radiometric scale from camera to camera.
        This perturbation is unique per `sensor`.
    """
    def __init__(self, cost_function, reflectance_values, SNR_values,
                 absolute_calibration_uncertainty=0.0,
                 camera_to_camera_calibration_uncertainty=0.0,
                seed=None):

        noise = reflectance_values/SNR_values
        noise_spline = si.CubicSpline(reflectance_values, SNR_values, extrapolate=True)

        def func(x, a, b):
            return a + b*np.sqrt(x)

        popt, pcov = so.curve_fit(func, reflectance_values, SNR_values, p0=[1e-5, 10])

        def snr_func(x):
            reflectance = np.atleast_1d(x)
            snr_out = noise_spline(reflectance)
            condition = np.where(reflectance > reflectance_values.max())
            snr_out[condition] = func(reflectance[condition], *popt)

            condition2 = np.where(reflectance < reflectance_values.min())
            snr_out[condition2] = reflectance[condition2]/noise[0]
            return snr_out

        super().__init__(cost_function, snr_func, noise[0])

### Some existing instruments.
### At present all noise is intensity only, neglecting contributions do to polarization.

class ResearchScanningPolarimeter(RadiometricUncertainty):

    def __init__(self, cost_function):

        super().__init__(cost_function, 2e-5, 1e-7, camera_to_camera_calibration_uncertainty=0.0,
                        absolute_calibration_uncertainty=np.sqrt(0.015))

class TandemStereoCamera(TabulatedRadiometricUncertainty):
    """
    The SNR requirements for the camera proposed for the
    Tandem Stereo Camera for NASA's Atmospheric Observing System (AOS)
    mission.

    The suggested camera_to_camera_calibration_uncertainty and
    absolute_calibration_uncertainty come from NASA's MISR instrument.
    """
    def __init__(self, cost_function, camera_to_camera_calibration_uncertainty=0.01,
                 absolute_calibration_uncertainty=0.03):

        super().__init__(
            cost_function,
            np.array([0.01, 0.05, 0.1, 0.5, 1.0, 1.3])/np.pi,
            np.array([87.0, 201.0, 285.0, 639.0, 904.0, 1031.0]),
            absolute_calibration_uncertainty=absolute_calibration_uncertainty,
            camera_to_camera_calibration_uncertainty=camera_to_camera_calibration_uncertainty
        )
