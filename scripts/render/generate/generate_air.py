"""TODO"""


import numpy as np
import shdom

def update_parser(parser):   
    parser.add_argument('--air_max_altitude',
                        default=20.0,
                        type=np.float32,
                        help='Maximum altitude for the air volume.')
    
    parser.add_argument('--air_num_grid_points',
                        default=20,
                        type=int,
                        help='Number of altitude grid points for the air volume.')
    
    return parser


def summer_mid_lat(args):
    """
    Interpolate temperatures and rayleigh scattering. 
    Temperature profile is taken from AFGL measurements of summer mid-lat.
    
    Parameters
    ----------
    args.wavelength: float
        Wavelength in [microns]
    args.air_max_altitude: float
        The maximum altitude for the air Medium
    args.air_num_grid_points: int
        Number of grid points for the air Medium.
        
    Returns
    -------
    TODO
    """
    temperatures = np.array([292.220, 292.040, 291.860, 291.680, 291.500, 291.320, 291.140, 290.960, 290.780, 
                             290.600, 290.420, 290.240, 290.060, 289.880, 289.700, 289.920, 290.140, 290.360, 
                             290.580, 290.800, 291.020, 291.240, 291.460, 291.680, 291.900])
   
    temp_grid = shdom.Grid(z=np.linspace(0.0, 20.0, len(temperatures)))
    temperature_profile = shdom.GridData(temp_grid, temperatures)
    rayleigh = shdom.Rayleigh(args.wavelength, temperature_profile)
    
    # Interpolate temperatures and rayleigh scattering up to max_altitude at grid resolution of
    altitudes = shdom.Grid(z=np.linspace(0.0, args.air_max_altitude, args.air_num_grid_points))
    extinction_a, albedo_a, phase_a = rayleigh.get_scattering_field(altitudes, phase_type='Grid')
    return extinction_a, albedo_a, phase_a