import xarray as xr
import numpy as np
import warnings
import shdom
import pandas as pd
import core

@xr.register_dataarray_accessor("shdomSampling")
class shdomSampling(object):
    def __init__(self, data_array):
        self._obj = data_array

    def resample(self, grid, method='linear'):
        temp_data = np.zeros((len(grid['x']), len(grid['y']), len(grid['z'])))
        temp_grid_data = xr.DataArray(
            data=temp_data, coords=[grid['x'], grid['y'], grid['z']], dims=['x', 'y', 'z']
        )

        for key, coord in grid.items():
            if key in self._obj.coords:
                if coord.max() > self._obj.coords[key].max():
                    warnings.warn('{} coordinate exceeded the grid maximum. Data has been truncated.'.format(key))
                if coord.min() < self._obj.coords[key].max():
                    warnings.warn('{} coordinate exceeded the grid minimum. Data has been truncated.'.format(key))

        interped = self._obj.interp_like(temp_grid_data, method=method)
        resampled, temp = xr.broadcast(interped, temp_grid_data)

        return resampled


#make some data
xs = np.linspace(0,1.0,15)
ys = np.linspace(0,1.0,16)
zs = np.linspace(0,1.0,17)

zs2 = np.linspace(0,1.2,32)

data = np.random.normal(size=(len(xs),len(ys),len(zs)))
data_1d = np.arange(len(zs2))

#lwc = xr.DataArray(name='lwc',data=data,coords=[xs,ys,zs], dims=['x','y','z'])
#reff = xr.DataArray(name='reff',data=data_1d,coords=[zs2],dims=['z'])

# Air grid
z_a = np.arange(0, 20, 1.0)
z_c = np.arange(0.5, 1.5, 0.04)

# Cloud grid
x = np.linspace(0.0,1.0,15)
y = np.linspace(0.0,1.0,16)
z = shdom.combine_1d_grids(z_c, z_a)

# Atmosphere grid
grid = pd.Series(index=['x','y','z'], data=[x, y, z])

lwc = xr.DataArray(data=data,coords=[xs,ys,zs], dims=['x','y','z'], name='lwc')
reff = xr.DataArray(data=data_1d, coords=[zs2], dims=['z'], name='reff')

# Resampling inplace
lwc.shdomSampling.resample(grid)

# Return new object
reff = reff.shdomSampling.resample(grid, inplace=False)



### Mie


class Mie(object):
    """
    Mie monodisperse scattering for spherical particles.

    Parameters
    ----------
    particle_type: string
        Options are 'Water' or 'Aerosol'.
    wavelength_band: (float, float)
        (minimum, maximum) wavelength in microns.
        This defines the spectral band over which to integrate, if both are equal monochrome quantities are computed.
    minimum_effective_radius: float
        Minimum effective radius in microns. Used to compute minimum radius for integration.
    max_integration_radius: float
        Maximum radius in microns - cutoff for the size distribution integral
    wavelength_averaging: bool
        True - average scattering properties over the wavelength_band.
        False - scattering properties of the central wavelength.
    wavelength_resolution: float
        The distance between two wavelength samples in the band. Used only if wavelength_averaging is True.
    refractive_index: complex number or path to refractive index table (CSV)
    """

    def __init__(self, particle_type, wavelength_band, minimum_effective_radius=4.0, max_integration_radius=65.0,
                 wavelength_averaging=False, wavelength_resolution=0.001, refractive_index=None):

        # Set wavelength integration parameters
        self._wavelen1, self._wavelen2 = wavelength_band
        assert self._wavelen1 <= self._wavelen2, 'Minimum wavelength is smaller than maximum'

        avgflag = 'C'
        if self._wavelen1 == self._wavelen2:
            deltawave = -1
        elif wavelength_averaging:
            avgflag = 'A'
            deltawave = wavelength_resolution

        self._avgflag = avgflag
        self._deltawave = deltawave

        self._wavelencen = core.get_center_wavelen(
            wavelen1=self._wavelen1,
            wavelen2=self._wavelen2
        )

        # Set particle type properties
        if (particle_type == 'Water'):
            self._partype = 'W'
            self._rindex = core.get_refract_index(
                partype=self._partype,
                wavelen1=self._wavelen1,
                wavelen2=self._wavelen2
            )

        elif (particle_type == 'Aerosol'):
            self._partype = 'A'
            assert refractive_index is not None, "Refractive index is not specified. \
            This could be a complex number or string path to csv (a function of wavelength)"
            if isinstance(refractive_index, str):
                rindex_df = pd.read_csv('../ancillary_data/dust_volz_1972.ri', comment='#', sep=' ',
                                        names=['wavelength', 'n', 'k'], index_col='wavelength')
                refractive_index = xr.Dataset.from_dataframe(rindex_df).interp({'wavelength': self._wavelencen})
                self._rindex = np.complex(refractive_index['n'], - refractive_index['k'])
            else:
                self._rindex = refractive_index

        else:
            raise AttributeError('Particle type note implemented')

        # Set radius integration parameters
        self._nsize = core.get_nsize(
            sretab=minimum_effective_radius,
            maxradius=max_integration_radius,
            wavelen=self._wavelencen
        )

        self._radii = core.get_sizes(
            sretab=minimum_effective_radius,
            maxradius=max_integration_radius,
            wavelen=self._wavelencen,
            nsize=self._nsize
        )

        # Calculate the maximum size parameter and the max number of Legendre terms
        if self._avgflag == 'A':
            xmax = 2 * np.pi * max_integration_radius / self._wavelen1
        else:
            xmax = 2 * np.pi * max_integration_radius / self._wavelencen
        self._maxleg = int(np.round(2.0 * (xmax + 4.0 * xmax ** 0.3334 + 2.0)))

    def compute_table(self):
        """
        Compute monodisperse Mie scattering per radius.

        Notes
        -----
        This is a time consuming method.
        """
        extinct, scatter, nleg, legcoef, table_type = \
            core.compute_mie_all_sizes(
                nsize=self._nsize,
                maxleg=self._maxleg,
                wavelen1=self._wavelen1,
                wavelen2=self._wavelen2,
                deltawave=self._deltawave,
                wavelencen=self._wavelencen,
                radii=self._radii,
                rindex=self._rindex,
                avgflag=self._avgflag,
                partype=self._partype
            )

        table = xr.Dataset(
            data_vars={
                'extinction': (['radius'], extinct),
                'scatter': (['radius'], scatter),
                'nleg': (['radius'], nleg),
                'legendre': (['stokes_index', 'legendre_index', 'radius'], legcoef)
            },
            coords={'radius': self._radii},
            attrs={
                'Table type': table_type.decode(),
                'radius units': '[micron]'
            },
        )
        return table

mie = Mie(particle_type='Aerosol',
          refractive_index='../ancillary_data/dust_volz_1972.ri',
          wavelength_band=(0.8, 0.8), max_integration_radius=20)
table_aerosol = mie.compute_table()

mie = Mie(particle_type='Water',
          wavelength_band=(0.8, 0.8), max_integration_radius=20)
table_water = mie.compute_table()
