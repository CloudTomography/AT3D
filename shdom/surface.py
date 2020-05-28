def fixed_lambertian_surface(albedo, ground_temperature=198.15):
    """
    TODO
    """
    return xr.Dataset(
        data_vars={
            'name': 'fixed_lambertian_surface',
            'sfctype':'FL',
            'gndalbedo':albedo,
            'gndtemp':ground_temperature,
            'maxfcpars':4,
            'nxsfc': 0,
            'nysfc': 0,
            'delxsfc': 0,
            'delysfc': 0,
            'nsfcpar': 1,
            'sfcparms': [],
            'sfcgridparms': []
        }
    )
