from unittest import TestCase
import os
import numpy as np
import xarray as xr
import pathlib

import at3d.aerosol

# This is a very limited test of the OPAC aerosol
# but internally there are a bunch of checks on things.
# This is basically a test that the tables are found etc.

class OPAC_mixture(TestCase):
    @classmethod
    def setUpClass(cls):

        aerosol_mixture = at3d.aerosol.OPACMixture()
        wavelengths = np.linspace(0.3, 0.8, 2)
        reference_aod = None

        for mixture_type in ['maritime_polluted']:

            humidity = xr.DataArray(
                data=np.ones(30)*0.5,
                dims='z',
                coords={
                    'z': np.linspace(0.0, 35.0, 30)
                    }
            )

            grid = at3d.grid.make_grid(0.05,2,0.05,2,z=humidity.z)
            humidity_on_grid = at3d.grid.resample_onto_grid(grid, humidity.to_dataset(name='humidity'))

            aerosol_scatterers = aerosol_mixture(atmosphere=humidity_on_grid ,
                                                 wavelengths=wavelengths, mixture_type=mixture_type,
                                                reference_aod=reference_aod, reference_wavelength=0.3)

            aaod = np.zeros(wavelengths.size)
            saod = np.zeros(wavelengths.size)

            for i, (key, scatterer) in enumerate(aerosol_scatterers.items()):
                absorb = (scatterer.extinction*(1.0 - scatterer.ssalb)).mean(['x', 'y'])
                scat = (scatterer.extinction*scatterer.ssalb).mean(['x', 'y'])
                aaod[i] = (0.5*(absorb[1:].data + absorb[:-1].data)*(scatterer.z[1:].data - scatterer.z[:-1].data)).sum()
                saod[i] = (0.5*(scat[1:].data + scat[:-1].data)*(scatterer.z[1:].data - scatterer.z[:-1].data)).sum()

            cls.aod = aaod[0] + saod[0]

    def test_aerosol_table_aod(self):
        self.assertTrue(np.allclose(0.16726762632365147, self.aod))
