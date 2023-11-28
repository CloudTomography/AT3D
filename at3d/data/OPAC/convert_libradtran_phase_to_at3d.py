# Convert opac table from libradtran to poly_table for use in at3d.
# This is necessary because SHDOM uses wigner function expansions of polarized
# phase functions rather than libradtran's legendre transforms.

# opac table to poly_table.

import at3d
import xarray as xr
import numpy as np
import os
import glob

import os
os.environ['at3d_DIR'] = '/Users/jesserl2/Documents/Code/pyshdom'
git_commit = at3d.util.github_version()
from datetime import datetime
date = datetime.now().date().strftime("%Y%m%d")

component_names={
            'inso': 'insoluble',
            'waso': 'water_soluble',
            'soot': 'soot',
            'ssam': 'sea_salt_accumulation_mode',
            'sscm': 'sea_salt_coarse_mode',
            'minm': 'mineral_nucleation_mode',
            'miam': 'mineral_accumulation_mode',
            'micm': 'mineral_coarse_mode',
            'mitr': 'mineral_transport',
            'suso': 'sulfate_soluble'
        }

directory='./aerosol'#'../../../data/OPAC/aerosol/'


for component_name, long_name in component_names.items():
    table = xr.load_dataset(os.path.join(directory,'optical_properties', component_name+ '.mie.cdf'))

    legcoef = np.zeros((6,table.pmom.nmommax.size,table.nhum.size, table.nlam.size))

    for i in range(table.nlam.size):
        for j in range(table.nhum.size):

            nrank = int(np.max(table.nmom[i,j])) - 1
            sample_points, weights = np.polynomial.legendre.leggauss(nrank)

            phase_components = []
            for k in range(4):

                mu = np.cos(np.deg2rad(table.theta[i,j,k].astype(np.float64)))
                p = table.phase[i,j,k].dropna('nthetamax').assign_coords({'nthetamax': mu.dropna('nthetamax')}).interp(
                    nthetamax=sample_points,
                    method='cubic'
                )
                phase_components.append(p)

            p1,p2,p3,p4 = phase_components
            coef = at3d.core.wigner_transform(
                p1=p1.data,
                p2=p2.data,
                p3=p3.data,
                p4=p4.data,
                mu=sample_points,
                nangle=p1.shape[0],
                wts=weights,
                nrank=nrank,
            )
            legcoef[:,:nrank+1,j,i] = coef

    legcoef[0,0,:,:] = 1.0

    indices = np.meshgrid(np.arange(table.nhum.size), np.arange(table.nlam.size), indexing='ij')
    table_index = np.ravel_multi_index(indices, dims=[table.nhum.size, table.nlam.size])

    poly_table = xr.Dataset(
        data_vars={
            'extinction': (['humidity','wavelength'], table.ext.data.T.astype(np.float32)),
            'ssalb': (['humidity','wavelength'], table.ssa.data.T.astype(np.float32)),
            'legcoef': (['stokes_index','legendre_index', 'humidity','wavelength'], legcoef.astype(np.float32))
        },
        coords={
            'wavelength': table.wavelen.data.astype(np.float32),
            'humidity': table.hum.data.astype(np.float32),
            'stokes_index': ['P11', 'P22', 'P33', 'P44', 'P12', 'P34'],
            #'table_index': (['humidity', 'wavelength'], table_index)
        },
        attrs={
            'file_info': 'OPAC aerosol optical properties extracted from libRadtran-2.0.3',
            'version': table.attrs['version'],
            'species_name': component_name,
            'long_name': long_name,
            'git_commit': git_commit,
            'creation_date': date,
        }
    )

    poly_table.to_netcdf(os.path.join(directory, 'optical_properties',component_name+ '.mie.nc'))
