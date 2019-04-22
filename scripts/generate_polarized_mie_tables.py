"""
Generate Polarized Mie tables 
-----------------------------

Generate mie tables, each table contains a Gamma distribution of drop size for multiple effective radii.
Multiple tables can be computed by specifying veff and wavelengths as a list 
  -- For command line flags info: python scripts/generate_polarized_mie_tables -h

Example usage:
  python scripts/generate_polarized_mie_tables.py --veff 0.01 0.02 0.05 0.1 --wavelength 0.443 0.67 0.865
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
    
    parser.add_argument('--output_dir', 
                        default='mie_tables',
                        help='(default value: %(default)s) Path to an output directory where the mie tables will be saved. \
                              If the folder doesnt exist, it will be created.')
    parser.add_argument('--veff',
                        nargs='+',
                        type=np.float32,
                        help='Effective droplet variance for which to compute the polarized mie table' \
                        'The output folder is formated as: <output_dir>/veff_<veff>')
    parser.add_argument('--wavelength',
                        nargs='+',
                        type=np.float32,
                        help='Wavelengths [Micron] for which to compute the polarized mie table' \
                        'The output file name is formated as: Water_<wavelength[nm]>nm_pol.scat')
    parser.add_argument('--num_reff',
                        default=50,
                        type=int,
                        help='(default value: %(default)s) Number of effective radii in the table')
    parser.add_argument('--start_reff',
                        default=8.0,
                        type=np.float32,
                        help='(default value: %(default)s) Starting effective radius [Micron]') 
    parser.add_argument('--end_reff',
                        default=20.0,
                        type=np.float32,
                        help='(default value: %(default)s) Ending effective radius [Micron]') 
    parser.add_argument('--radius_cutoff',
                        default=60.0,
                        type=np.float32,
                        help='(default value: %(default)s) The cutoff radius for the pdf averaging [Micron]')    
    parser.add_argument('--log_space',
                        action='store_true',
                        help='Logarithmic spacing for the effective radius')    
    
    args = parser.parse_args()
    
    assert args.veff is not None, 'Requiered effective variance was not input. \n' \
           'For command line flags see: python generate_polarized_mie_tables.py -h'
    
    assert args.wavelength is not None, 'Requiered wavelength was not input. \n' \
           'For command line flags see: python generate_polarized_mie_tables.py -h'
    
    return args


if __name__ == "__main__":
    
    args = argument_parsing()
    
    for veff in args.veff:
        alpha = 1.0/veff - 3.0
        directory = 'mie_tables/veff_{:.3f}'.format(veff)
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        for wavelen in args.wavelength:
            mie = shdom.MiePolarized()
            mie.set_parameters((wavelen, wavelen), 'Water', 'gamma', alpha)
            mie.compute_table(args.num_reff, args.start_reff, args.end_reff, args.radius_cutoff, args.log_space)
            output_path = os.path.join(directory, 'Water_{:d}nm_pol.scat'.format(int(np.round(wavelen*1000))))
            mie.write_table(output_path)
        