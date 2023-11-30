"""This module contains methods for defining the three different types
of sources used in SHDOM; `solar`, `thermal`, and `combined`, which
produce collimated flux from above, thermal emission or a combination of the two.
The source also defines the top boundary condition, which may be a
downwelling isotropic radiation field.
Note that only monochromatic sources are currently supported,
so spectral integration schemes must be implemented 'manually' at the script level.
Gradient calculations are only available for `solar` sources, not
`thermal` or `combined`.
"""
import xarray as xr
import numpy as np

def solar(wavelength, solarmu, solar_azimuth, solarflux=1.0, skyrad=0.0, volume_source=None):
    """Defines a solar source for use in solver.RTE (SHDOM).

    Parameters
    ----------
    wavelength : float
        The monochromatic (center) wavelength for use in the solar source in [microns]
    solarmu : float
        cosine of solar zenith angle. Points in the direction of propagation of
        the solar beam. Is required to be negative, if positive values are supplied
        they are set to negative.
    solar_azimuth : float
        The azimuthal angle [degrees] of the solar beam. Points in the direction of
        propagation of the solar beam.
    solarflux : float
        The flux of the solar beam on a horizontal surface.
    skyrad : float
        The RADIANCE of the downwelling isotropic radiance field at domain top.

    Returns
    -------
    source_dataset : xr.Dataset
        The dataset contains the flags, units and above parameters that define
        the source for solver.RTE (SHDOM).

    Raises
    ------
    ValueError
        If `solarmu` is outside of the range [-1.0, 0.0). Note that `solarmu`=0.0
        (purely horizontal solar beam) is not permitted.

    Notes
    -----
    Several of the numerical parameters used in SHDOM such as the splitting accuracy
    and spherical harmonics accuracy are absolute measures and so they should scale
    with the magnitude of the radiation field and therefore the prescribed source.
    To avoid having to retune the numerical parameters, it is advisible to set
    `solarflux`=1.0 for all purely solar problems.
    """
    solarmu = -1*np.abs(solarmu)
    if not (-1.0 <= solarmu) & (solarmu < 0.0):
        raise ValueError("solarmu must be in the range -1.0 <= solarmu < 0.0 not '{}'. "
                         "The SHDOM convention for solar direction is that it points"
                         "in the direction of the propagation of radiance.".format(solarmu))

    source_dataset = xr.Dataset(
        data_vars={
            'name': 'solar_source',
            'wavelength': wavelength,
            'solarflux': solarflux,
            'solarmu': solarmu,
            'solaraz': np.deg2rad(solar_azimuth),
            'srctype': 'S',
            'units': 'R',
            'wavenumber': [10000, 10001], #only used for CKD
            'skyrad': (['nstokes', 'nmu','nphi0max'], np.atleast_3d(skyrad)) #isotropic diffuse radiance from above
        }
    )
    if volume_source is None:
        volume_source = VolumeSource()
    if not isinstance(volume_source, VolumeSource):
        raise TypeError(
            "`volume_source` should inherit from {}".format(VolumeSource)
        )
    source_dataset['volume_source'] = volume_source
    return source_dataset

def thermal(wavelength, skyrad=0.0, units='radiance', volume_source=None):
    """Defines a thermal emission source for solver.RTE (SHDOM).

    Parameters
    ----------
    wavelength : float
        The monochromatic (center) wavelength for use in the solar source in [microns]
    skyrad : float
        The BRIGHTNESS TEMPERATURE of the downwelling isotropic radiance field at domain top.
    units : string
        Determines whether to set units to spectral radiance [W/m^2/steradian/micron] or
        to units of brightness temperature [K].

    Returns
    -------
    source_dataset : xr.Dataset
        The dataset contains the flags, units and above parameters that define
        the source for solver.RTE (SHDOM).
    Raises
    ------
    ValueError
        If the prescribed units are not either 'radiance' or 'brightness_temperature'

    Notes
    -----
    skyrad is now in units of BRIGHTNESS TEMPERATURE.
    When using a thermal source, it is not normalizable like the solar source.
    This means that care must be taken with the prescription of the accompanying
    numerical parameters such as the splitting accuracy.
    A good rule of thumb is to proportionally scale splitting accuracy and
    spherical harmonics accuracy based on the largest flux in the domain,
    which can be approximated as the flux emitted by the surface.
    """
    if units == 'radiance':
        units_flag = 'R'
    elif units == 'brightness_temperature':
        units_flag = 'T'
    else:
        raise ValueError("`units` should be either 'radiance' or 'brightness_temperature'.")

    source_dataset = xr.Dataset(
        data_vars={
            'name': 'solar_source',
            'wavelength': wavelength,
            'solarflux': 0.0,
            'solarmu': -0.5,
            'solaraz': 0.0,
            'srctype': 'T',
            'units': units_flag,
            'wavenumber': [10000, 10001], #only used for CKD
            'skyrad': skyrad #in thermal only this is brightness temperature of the isotropic diffuse radiance.
        }
    )
    if volume_source is None:
        volume_source = VolumeSource()
    if not isinstance(volume_source, VolumeSource):
        raise TypeError(
            "`volume_source` should inherit from {}".format(VolumeSource)
        )
    source_dataset['volume_source'] = volume_source
    return source_dataset

def combined(wavelength, solarmu, solar_azimuth, solarflux=1.0, skyrad=0.0, volume_source=None):
    """Defines a combined (both solar and thermal) source for use in solver.RTE (SHDOM).

    Parameters
    ----------
    wavelength : float
        The monochromatic (center) wavelength for use in the solar source in [microns]
    solarmu : float
        cosine of solar zenith angle. Points in the direction of propagation of
        the solar beam. Is required to be negative, if positive values are supplied
        they are set to negative.
    solar_azimuth : float
        The azimuthal angle of the solar beam. Points in the direction of
        propagation of the solar beam.
    solarflux : float
        The flux of the solar beam on a horizontal surface.
    skyrad : float
        The RADIANCE of the downwelling isotropic radiance field at domain top.

    Returns
    -------
    source_dataset : xr.Dataset
        The dataset contains the flags, units and above parameters that define
        the source for solver.RTE (SHDOM).

    Raises
    ------
    ValueError
        If `solarmu` is outside of the range [-1.0, 0.0). Note that `solarmu`=0.0
        (purely horizontal solar beam) is not permitted.

    Notes
    -----
    Units of output (including skyrad) are in radiance [W/m^2/steradian/micron].
    When using a combined source, the `solarflux` should be set to its non-normalized
    value so that its relative magnitude with the thermal emission is correct.
    In this case, the split accuracy and spherical harmonics accuracy should again be
    correspondingly scaled, these could be chosen as a proportionally scaling by
    the sum of `solarflux` and the flux emitted from the surface, for example.
    """
    solarmu = -1*np.abs(solarmu)
    if not (-1.0 <= solarmu) & (solarmu < 0.0):
        raise ValueError("solarmu must be in the range -1.0 <= solarmu < 0.0 not '{}'. "
                         "The SHDOM convention for solar direction is that it points"
                         "in the direction of the propagation of radiance.".format(solarmu))

    source_dataset = xr.Dataset(
        data_vars={
            'name': 'solar_source',
            'wavelength': wavelength,
            'solarflux': solarflux,
            'solarmu': solarmu,
            'solaraz': np.deg2rad(solar_azimuth),
            'srctype': 'B',
            'units': 'R',
            'wavenumber': [10000, 10001], #only used for CKD
            'skyrad': skyrad #In combined this is the same as solar (radiance)
        }
    )
    if volume_source is None:
        volume_source = VolumeSource()
    if not isinstance(volume_source, VolumeSource):
        raise TypeError(
            "`volume_source` should inherit from {}".format(VolumeSource)
        )
    source_dataset['volume_source'] = volume_source
    return source_dataset


class VolumeSource:
    """
    A generic class for defining (non-thermal) volume source.
    """
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def _check_inputs(self, *args):
        for arg in args:
            if not isinstance(arg, np.ndarray):
                raise TypeError(
                    "Coordinates should be passed as numpy arrays."
                )
            if not arg.ndim == 1:
                raise ValueError(
                    "Coordinates should be arrays of rank 1"
                )
        if not all([arg.size==args[0].size for arg in args]):
            raise ValueError(
                "Coordinates should all have the same size."
            )

    def _rad_func(self, x,y,z,mu,phi):
        return np.zeros(x.shape)

    def _check_output(self, rad):
        if np.any(~np.isfinite(rad)):
            raise ValueError(
                "Non-finite values in output of volume source model"
            )
    
    def __call__(self, x, y,z, mu, phi):
        self._check_inputs(x,y,z,mu,phi)
        rad_out = self._rad_func(x,y,z,mu,phi)
        self._check_output(rad_out)
        return rad_out

class IsotropicPointSource(VolumeSource):
    """
    An Isotropic Volume point source.

    Parameters
    ----------
    x0, y0, z0 : float
        The location of the point source.
    dr : float
        The radius of the point source. This parameter is provided so that there is some tolerance
        if the location does not perfectly match up with one of the radiative transfer grid's
        points.
    intensity : float
        The intensity of the radiance at every angle.
    """
    def __init__(self, x0=0.32,y0=0.32,z0=1.0,dr=0.05, intensity=1.0/(4*np.pi)):
        VolumeSource.__init__(self, x0=x0, y0=y0, z0=z0, dr=dr, intensity=intensity)

    def _rad_func(self, x, y, z, mu, phi):
        self._check_inputs(x,y,z,mu,phi)
        d = (x-self.x0)**2 + (y-self.y0)**2 + (z-self.z0)**2
        cond = np.where(np.sqrt(d)<= self.dr)
        rad_out = np.zeros((x.shape))
        rad_out[cond] = self.intensity
        return rad_out