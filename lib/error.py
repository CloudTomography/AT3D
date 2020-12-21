import numpy as np
import xarray as xr

class IndependentNoise:
    def __init__(self, noise_strength, relative=False):

        if np.atleast_1d(noise_strength).size == 1:
            self._error_diagonal = noise_strength
            self._max_nstokes = 1
        else:
            raise NotImplementedError
        self._relative = relative

    def add_noise(self, sensor, nstokes):
        if 'uncertainties' not in sensor:
            self.calculate_uncertainties(sensor, nstokes)
        sensor['I'] += np.random.normal(size=sensor.sizes['npixels'])*(sensor['uncertainties'][0,0].data**(-0.5))

    def calculate_uncertainties(self, sensor, nstokes):
        uncertainty_data = np.zeros((nstokes, nstokes,sensor.sizes['npixels']))
        if self._relative:
            uncertainty_data[0,0] = (sensor.I.data*self._error_diagonal)**2
        else:
            uncertainty_data[0,0] = (self._error_diagonal)**2
        uncertainty_data[np.where(uncertainty_data > 1e-9)] = 1.0/uncertainty_data[np.where(uncertainty_data > 1e-9)]
        uncertainty_data[np.where(uncertainty_data <= 1e-9)] = 1e9
        sensor['uncertainties'] = (['nstokes', 'nstokes2', 'npixels'], uncertainty_data)
