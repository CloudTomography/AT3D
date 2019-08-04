"""
Generate Mie tables 
-------------------

Generate mie tables, each table contains a Gamma distribution of drop size for multiple effective radii and variances.
Multiple tables can be computed by specifying wavelength as a list 
  -- For command line flags info: python scripts/generate_mie_tables -h

Example usage:
  python scripts/generate_mie_tables.py --wavelength 0.443 0.67 0.865
"""

import os, shdom, argparse
import numpy as np


def argument_parsing():
    """
    Handle all the argument parsing needed for this script.
    
    Returns
    -------
    args: arguments from argparse.ArgumentParser()
        The arguments requiered for this script.
    """
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--poly_dir', 
                        default='mie_tables/polydisperse/',
                        help='(default value: %(default)s) Path to an output directory where the Mie poly tables will be saved. \
                              If the folder doesnt exist, it will be created.')
    parser.add_argument('--mono_dir', 
                        default='mie_tables/monodisperse/',
                        help='(default value: %(default)s) Path to an output directory where the mie mono tables will be saved. \
                              If the folder doesnt exist, it will be created.')    
    parser.add_argument('--wavelength',
                        nargs='+',
                        type=np.float32,
                        help='Wavelengths [Micron] for which to compute the polarized mie table' \
                        'The output file name is formated as: Water_<wavelength[nm]>nm.scat<pol>')
    parser.add_argument('--num_reff',
                        default=50,
                        type=int,
                        help='(default value: %(default)s) Number of effective radii in the table')
    parser.add_argument('--start_reff',
                        default=4.0,
                        type=np.float32,
                        help='(default value: %(default)s) Starting effective radius [Micron]') 
    parser.add_argument('--end_reff',
                        default=25.0,
                        type=np.float32,
                        help='(default value: %(default)s) Ending effective radius [Micron]')
    parser.add_argument('--num_veff',
                        default=50,
                        type=int,
                        help='(default value: %(default)s) Number of effective variance in the table')
    parser.add_argument('--start_veff',
                        default=0.01,
                        type=np.float32,
                        help='(default value: %(default)s) Starting effective variance') 
    parser.add_argument('--end_veff',
                        default=0.2,
                        type=np.float32,
                        help='(default value: %(default)s) Ending effective variance')     
    parser.add_argument('--radius_cutoff',
                        default=65.0,
                        type=np.float32,
                        help='(default value: %(default)s) The cutoff radius for the pdf averaging [Micron]')    
    parser.add_argument('--polarized',
                        action='store_true',
                        help='Polarized/Unpolarized table.')    
    
    args = parser.parse_args()
    
    assert args.wavelength is not None, 'Requiered wavelength was not input. \n' \
           'For command line flags see: python generate_polarized_mie_tables.py -h'
    
    return args


def get_file_paths(wavelength, args):
    """TODO"""
    if not os.path.exists(args.mono_dir):
        os.makedirs(args.mono_path)    
    if not os.path.exists(args.poly_dir):
        os.makedirs(args.poly_path)   
    file_name = 'Water_{:3d}nm.scat'.format(shdom.int_round(wavelength))
    if args.polarized:
        file_name += 'pol'
    mono_path = os.path.join(args.mono_dir, file_name)
    poly_path = os.path.join(args.poly_dir, file_name)
    return mono_path, poly_path


def compute_mie_table(wavelength, args):
    """TODO"""
    mono_path, poly_path = get_file_paths(wavelength, args)

    mie_mono = shdom.MieMonodisperse()
    if os.path.exists(mono_path):
        mie_mono.read_table(mono_path)
    else:
        mie_mono.set_wavelength_integration(wavelength_band=(wavelength, wavelength))   
        mie_mono.set_radius_integration(minimum_effective_radius=args.start_reff, 
                                        max_integration_radius=args.radius_cutoff)
        mie_mono.compute_table()
        mie_mono.write_table(mono_path)
        
    size_dist = shdom.SizeDistribution(type='gamma')
    size_dist.set_parameters(reff=np.linspace(args.start_reff, args.end_reff, args.num_reff), 
                                     veff=np.linspace(args.start_veff, args.end_veff, args.num_veff))
    size_dist.compute_nd(radii=mie_mono.radii, particle_density=mie_mono.pardens)
    mie_poly = shdom.MiePolydisperse(mie_mono, size_dist)
    mie_poly.compute_table()
    mie_poly.write_table(poly_path)


if __name__ == "__main__":
    
    args = argument_parsing()
    for wavelength in args.wavelength:
        compute_mie_table(wavelength, args)
