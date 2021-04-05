"""
This module contains functions to map microphysical variables to
optical properties. In the original SHDOM, these procedures would be performed
by the `propgen` program.
Currently only a Look-Up-Table approach is supported defined by `table_to_grid`
which uses linear interpolation and nearest interpolation for phase functions.
This method is also used to map partial derivatives of optical properties
from the Look-Up-Table onto the SHDOM grid.
"""
import numpy as np
import xarray as xr
import pyshdom.checks

def table_to_grid(microphysics, poly_table, exact_table=False, inverse_mode=False):
    """
    Calculates optical properties from microphysical properties using a Look-Up-Table.

    Optical properties are calculated on the 3D grid defined in `microphysics`
    using microphysical parameters defined in `microphysics` and a Look-Up-Table to
    map between microphysical variables and optical variables in `poly_table`
    using linear interpolation.

    Parameters
    ----------
    microphysics : xr.Dataset
        Should contain the microphysical parameters on a valid grid.
        See grid.py for details. The microphysical parameters should match
        those in `poly_table`.
    poly_table : xr.Dataset
        Should contain extinction efficiency, single scatter albedo and phase
        functions as a function microphysical parameters (e.g. effective radius).
    exact_table : bool
        Sets the interpolation method for the calculation of extinction and ssalb.
        linear if False, and nearest if True.
    inverse_mode : bool
        A flag for whether optical properties or their derivatives with respect
        to microphysical properties are being interpolated.
        The only difference is that extinction_efficiency passed instead of
        extinction.

    Returns
    -------
    optical_properties : xr.Dataset
        A dataset containing optical properties defined on an SHDOM grid
        ready to be used as input to solver.RTE.

    Notes
    -----
    The linear interpolation of extinction and single scatter albedo and
    nearest neighbor interpolation of phase functions used here differs
    from SHDOM's propgen program, which uses a different interpolation scheme
    and creates new phase functions if a specified accuracy is not met.
    This is not implemented as the Look-Up-Table approach implemented here
    is simple to apply to the inverse problem without the need to recompute
    any mie properties. To ensure high accuracy, the spacing of the table
    should be fine, see size_distribution.py.
    """
    pyshdom.checks.check_positivity(microphysics, 'density')
    pyshdom.checks.check_grid(microphysics)
    if not inverse_mode:
        pyshdom.checks.check_legendre(poly_table)

    interp_names = set([name for name in poly_table.coords
                        if name not in ('table_index', 'stokes_index')])
    microphysics_names = set([name for name in microphysics.variables.keys()
                              if name not in 'density'])
    missing = interp_names - microphysics_names
    if missing:
        raise KeyError(
            "microphysics dataset is missing variables "
            "for interpolation of table onto grid.", *list(missing)
            )
    interp_coords = {name:microphysics[name] for name in poly_table.coords
                     if name not in ('table_index', 'stokes_index')}
    for interp_coord in interp_coords:
        if np.any(microphysics[interp_coord] <= poly_table[interp_coord].min()) or \
            np.any(microphysics[interp_coord] >= poly_table[interp_coord].max()):
            raise ValueError(
                "Microphysical coordinate '{}' is not"
                " within the range of the mie table.".format(interp_coord)
                )
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
    table_index = poly_table.coords['table_index'].interp(
        coords=interp_coords, method='nearest'
        ).round().astype(int)
    unique_table_indices, inverse = np.unique(table_index.data, return_inverse=True)
    subset_table_index = xr.DataArray(name=table_index.name,
                                      data=inverse.reshape(table_index.shape) + 1,
                                      dims=table_index.dims,
                                      coords=table_index.coords,
                                      )

    legendre_table_stack = poly_table['legcoef'].stack(table_index=interp_coords)
    subset_legcoef = legendre_table_stack.isel({'table_index':unique_table_indices})
    subset_legcoef = xr.DataArray(name='legcoef',
                                  data=subset_legcoef.data,
                                  dims=['stokes_index', 'legendre_index', 'table_index'],
                                  coords={'stokes_index':subset_legcoef.coords['stokes_index']}
                                 )

    optical_properties = xr.merge([extinction, ssalb, subset_legcoef])
    optical_properties['density'] = microphysics.density
    table_coords = {'table_index': (['x', 'y', 'z'], subset_table_index.data)}

    optical_properties = optical_properties.assign_coords(table_coords)
    assert not np.any(np.isnan(optical_properties.table_index.data)), 'Unexpected NaN in table_index'
    assert not np.any(np.isnan(optical_properties.legcoef.data)), 'Unexpected NaN in legcoef'

    #inherit attributes.
    optical_properties = optical_properties.assign_attrs(poly_table.attrs)
    optical_properties = optical_properties.assign_attrs(microphysics.attrs)

    #transfer the grid variables. NOTE that delx, dely exist and be passed.
    # while nx/ny/nz are optional. delx/dely are checked for by check_grid.
    grid_variables = ('delx', 'dely', 'nx', 'ny', 'nz')
    for grid_variable in grid_variables:
        if grid_variable in microphysics.data_vars:
            optical_properties[grid_variable] = microphysics[grid_variable]

    return optical_properties
