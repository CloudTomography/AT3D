"""
TODO: description
"""

import numpy as np
import xarray as xr
from shdom import checks

@checks.dataset_checks(microphysics=(checks.check_positivity, 'density'))
def table_to_grid(microphysics, poly_table, exact_table=False, inverse_mode=False):
    """
    TODO
    Interpolates the poly table onto the microphysics grid.

    #TODO add checks that poly_table spans the microphysics. If only a single value
    is in the poly_table (from the size distribution) then nans fill optical properties
    silently as linear interpolation fails.
    """
    #CHECK FOR CORRECTLY NORMALIZED POLY_TABLE
    #check for missing microphysical coordinates.
    interp_names = set([name for name in poly_table.coords if name not in ('table_index', 'stokes_index')])
    microphysics_names = set([name for name in microphysics.variables.keys() if name not in 'density'])
    missing = interp_names - microphysics_names
    if len(list(missing)) > 0:
        raise KeyError("microphysics dataset is missing variables for interpolation of table onto grid.", *list(missing))

    interp_coords = {name:microphysics[name] for name in poly_table.coords if name not in ('table_index', 'stokes_index')}
    for interp_coord in interp_coords:
        if np.any(microphysics[interp_coord] <= poly_table[interp_coord].min()) or \
            np.any(microphysics[interp_coord] >= poly_table[interp_coord].max()):
            raise ValueError("Microphysical coordinate '{}' is not within the range of the mie table.".format(
                                                                                        interp_coord))

    interp_method = 'nearest' if exact_table else 'linear'
    ssalb = poly_table.ssalb.interp(interp_coords, method=interp_method)
    assert not np.any(np.isnan(ssalb.data)), 'Unexpected NaN in ssalb'
    extinction_efficiency = poly_table.extinction.interp(interp_coords, method=interp_method)

    if not inverse_mode:
        extinction = extinction_efficiency * microphysics.density
    else:
        extinction = extinction_efficiency
    assert not np.any(np.isnan(extinction.data)), 'Unexpected NaN in extinction'
    extinction.name = 'extinction'
    table_index = poly_table.coords['table_index'].interp(coords=interp_coords, method='nearest').round().astype(int)
    unique_table_indices, inverse = np.unique(table_index.data, return_inverse=True)
    subset_table_index = xr.DataArray(name=table_index.name,
                                        data = inverse.reshape(table_index.shape) + 1,
                                        dims=table_index.dims,
                                        coords=table_index.coords,
                                        )

    legendre_table_stack = poly_table['legcoef'].stack(table_index=interp_coords)
    subset_legcoef = legendre_table_stack.isel({'table_index':unique_table_indices})
    subset_legcoef = xr.DataArray(name='legcoef',
                                    data = subset_legcoef.data,
                                    dims=['stokes_index', 'legendre_index', 'table_index'],
                                    coords = {'stokes_index':subset_legcoef.coords['stokes_index']}
                                    )

    optical_properties = xr.merge([extinction, ssalb, subset_legcoef])
    optical_properties['density'] = microphysics.density
    table_coords = {'table_index': (['x','y','z'],subset_table_index.data)}

    optical_properties =optical_properties.assign_coords(table_coords)
    assert not np.any(np.isnan(optical_properties.table_index.data)), 'Unexpected NaN in table_index'
    assert not np.any(np.isnan(optical_properties.legcoef.data)), 'Unexpected NaN in legcoef'

    #inherit attributes.
    optical_properties = optical_properties.assign_attrs(poly_table.attrs)
    optical_properties = optical_properties.assign_attrs(microphysics.attrs)
    optical_properties['delx'] = microphysics.delx
    optical_properties['dely'] = microphysics.dely

    #TODO Check for whether these variables are 'in-place' modifications that
            #affect the earlier legendre_table/poly_table etc.
    return optical_properties


# def get_bounding_box(medium):
#     bounding_box = xr.DataArray(
#         data=[medium.x[0].data, medium.y[0].data, medium.z[0].data,
#               medium.x[-1].data, medium.y[-1].data, medium.z[-1].data],
#         coords={'bb_index': ['xmin', 'ymin', 'zmin', 'xmax', 'ymax', 'zmax']},
#         dims='bb_index',
#         attrs={'type': '3D bounding box'}
#     )
#     # TODO add units from the grid to the bounding box if exists
#     return bounding_box
