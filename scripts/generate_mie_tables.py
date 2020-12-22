"""
Generate Mie tables 
-------------------
Generate mie tables, each table contains a Gamma distribution of drop size for multiple effective radii and variances.
Multiple tables can be computed by specifying wavelength as a list 

For information about the command line flags see:
  python scripts/generate_mie_tables.py --help
  
For example usage see the README.md
"""

import os, pyshdom, argparse
import numpy as np


def argument_parsing():
    """
    Handle all the argument parsing needed for this script.
    
    Returns
    -------
    args: arguments from argparse.ArgumentParser()
        The arguments required for this script
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
    parser.add_argument('--bandwidth',
                        nargs='+',
                        type=np.float32,
                        help='Spectral bandwidth [nm]')
    parser.add_argument('--wavelength_resolution',
                        default=5,
                        type=np.float32,
                        help='(default value: %(default)s) Wavelength integration resolutions [nm]')
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
    parser.add_argument('--mie_base_name',
                        default='Water_<wavelength>nm.scat',
                        help='(default value: %(default)s) Mie table base file name. ' \
                        '<wavelength> will be replaced by the corresponding wavelengths.')      
    
    args = parser.parse_args()
    
    assert args.wavelength is not None, 'required wavelength was not input. \n' \
           'For command line flags see: python generate_polarized_mie_tables.py -h'
    
    return args


def get_file_paths(wavelength, args):
    """
    Retrieve the file paths according to the wavelength and base path argument.
    
    Parameters
    ----------
    wavelength: float
            Wavelength in microns
    args: arguments from argparse.ArgumentParser()
        Arguments required for this function
    
    Returns
    -------
    mono_path: str
        Path to the monodisperse table
    poly_path: str
        Path to the polydisperse table
    """
    if not os.path.exists(args.mono_dir):
        os.makedirs(args.mono_dir)    
    if not os.path.exists(args.poly_dir):
        os.makedirs(args.poly_dir)   
         
    file_name = args.mie_base_name.replace('<wavelength>', '{}'.format(pyshdom.int_round(wavelength)))
    if args.polarized:
        file_name += 'pol'
    mono_path = os.path.join(args.mono_dir, file_name)
    poly_path = os.path.join(args.poly_dir, file_name)
    return mono_path, poly_path


def compute_mie_table(wavelength, bandwidth, args):
    """
    Compute and save the monodisperse and polydisperse Mie tables.
    
    Parameters
    ----------
    wavelength: float
        Wavelength in microns
    bandwidth: float
        Wavelength bandwidth in nanometer
    args: arguments from argparse.ArgumentParser()
        Arguments required for this function
    """
    mono_path, poly_path = get_file_paths(wavelength, args)

    mie_mono = pyshdom.MieMonodisperse()
    if os.path.exists(mono_path):
        mie_mono.read_table(mono_path)
    else:
        if bandwidth is not None:
            wavelength_averaging = True
            wavelength_band = (wavelength - bandwidth/2000.0, wavelength + bandwidth/2000.0)
        else:
            wavelength_averaging = False
            wavelength_band = (wavelength, wavelength)
        mie_mono.set_wavelength_integration(wavelength_band=wavelength_band,
                                            wavelength_averaging=wavelength_averaging,
                                            wavelength_resolution=args.wavelength_resolution/1000.0)
        mie_mono.set_radius_integration(minimum_effective_radius=args.start_reff, 
                                        max_integration_radius=args.radius_cutoff)
        mie_mono.compute_table()
        mie_mono.write_table(mono_path)
        
    size_dist = pyshdom.SizeDistribution(type='gamma')
    size_dist.set_parameters(reff=np.linspace(args.start_reff, args.end_reff, args.num_reff),
                             veff=np.linspace(args.start_veff, args.end_veff, args.num_veff))
    size_dist.compute_nd(radii=mie_mono.radii, particle_density=mie_mono.pardens)
    mie_poly = pyshdom.MiePolydisperse(mie_mono, size_dist)
    mie_poly.compute_table()
    mie_poly.write_table(poly_path)


if __name__ == "__main__":
    args = argument_parsing()
    bandwidths = np.atleast_1d(args.bandwidth)
    if len(bandwidths) == 1:
        bandwidths = np.full(len(args.wavelength), bandwidths)
    else:
        assert len(bandwidths) == len(args.wavelength), \
            'Number of bandwidths should be either 1 or the same as wavelengths'
    for wavelength, bandwidth in zip(args.wavelength, bandwidths):
        compute_mie_table(wavelength, bandwidth, args)
