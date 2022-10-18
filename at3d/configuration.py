"""
Generates and saves/loads numerical parameters for SHDOM.
json is used as it is human readable rather than netCDF.

See ./shdom.txt for more details of the numerical parameters.
"""

import json
import xarray as xr

def make_config(config_file_name, x_boundary_condition='open',
                y_boundary_condition='open',
                num_mu_bins=16, num_phi_bins=32,
                split_accuracy=0.03, deltam=True,
                spherical_harmonics_accuracy=0.0,
                solution_accuracy=0.0001,
                acceleration_flag=True,
                max_total_mb=3000.0,
                adapt_grid_factor=5,
                num_sh_term_factor=1,
                cell_to_point_ratio=1.5,
                high_order_radiance=False,
                ip_flag=0,
                iterfixsh=30,
                tautol=0.1,
                angle_set=2,
                transcut=1e-5,
                transmin=1.0):

    shdom_parameters = make_config_data(
        x_boundary_condition=x_boundary_condition,
        y_boundary_condition=y_boundary_condition,
        num_mu_bins=num_mu_bins,
        num_phi_bins=num_phi_bins,
        split_accuracy=split_accuracy,
        deltam=deltam,
        spherical_harmonics_accuracy=spherical_harmonics_accuracy,
        solution_accuracy=solution_accuracy,
        acceleration_flag=acceleration_flag,
        max_total_mb=max_total_mb,
        adapt_grid_factor=adapt_grid_factor,
        num_sh_term_factor=num_sh_term_factor,
        cell_to_point_ratio=cell_to_point_ratio,
        high_order_radiance=high_order_radiance,
        ip_flag=ip_flag,
        iterfixsh=iterfixsh,
        tautol=tautol,
        angle_set=angle_set,
        transcut=transcut,
        transmin=transmin
    )
    json.dump(shdom_parameters, open(config_file_name +'.json', 'w'), indent=2)


def make_config_data(x_boundary_condition='open',
                y_boundary_condition='open',
                num_mu_bins=16, num_phi_bins=32,
                split_accuracy=0.03, deltam=True,
                spherical_harmonics_accuracy=0.0,
                solution_accuracy=0.0001,
                acceleration_flag=True,
                max_total_mb=3000.0,
                adapt_grid_factor=5,
                num_sh_term_factor=1,
                cell_to_point_ratio=1.5,
                high_order_radiance=False,
                ip_flag=0,
                iterfixsh=30,
                tautol=0.1,
                angle_set=2,
                transcut=1e-5,
                transmin=1.0
                ):
    """
    See default_config.json for description of parameters along with shdom.txt
    or read this function's code.
    """
    #if config_file_name == 'default_config':
        #raise ValueError("Config file cannot overwrite default")

    shdom_parameters = {
        'x_boundary_condition': {
            'default_value': x_boundary_condition,
            'description':'1. open - exiting radiance is lost. 2. periodic - exiting radiance returns from the opposite side.'
        },
        'y_boundary_condition': {
            'default_value': y_boundary_condition,
            'description':'1. open - exiting radiance is lost. 2. periodic - exiting radiance returns from the opposite side.'
        },
        'num_mu_bins': {
            'default_value': num_mu_bins,
            'description': '(NMU) number of discrete ordinates covering -1 < mu < 1.'
        },
        'num_phi_bins': {
            'default_value': num_phi_bins,
            'description': '(NPHI) number of discrete ordinates covering 0 < phi < 2pi'
        },
        'split_accuracy': {
            'default_value': split_accuracy,
            'description':
                ["(SPLITACC) cell splitting accuracy; grid cells that have the adaptive splitting criterion above this value are split.",
                 "This is an absolute measure, but cannot be easily associated with the resulting radiometric accuracy. Set to zero or negative for no adaptive cell splitting."]
         },
        'deltam': {
            'default_value': deltam,
            'description': '(DELTAM) True for delta-M scaling of medium and Nakajima and Tanaka method of computing radiances.'
        },
        'spherical_harmonics_accuracy': {
            'default_value': spherical_harmonics_accuracy,
            'description':
                ["(SHACC) adaptive spherical harmonic truncation accuracy; the spherical harmonic source function series is truncated after the terms are below this level.",
                 "Truncation can still happens if SHACC=0 (for 0 source terms).This is also an absolute measure, and is approximately the level of accuracy."]
        },
        'acceleration_flag': {
            'default_value':acceleration_flag,
            'description':'(ACCELFLAG) True to do the sequence acceleration. An acceleration extrapolation of the source function is done every other iteration.'
        },
        'solution_accuracy': {
            'default_value':solution_accuracy,
            'description': '(SOLACC) solution accuracy - tolerance for solution criterion.'
        },
        'max_total_mb': {
            'default_value': max_total_mb,
            'description': 'approximate maximum memory to use (MB for 4 byte reals)'
        },
        'adapt_grid_factor': {
            'default_value': adapt_grid_factor,
            'description': '(ADAPT_GRID_FACTOR) ratio of total grid points to base grid points'
        },
        'num_sh_term_factor': {
            'default_value': num_sh_term_factor,
            'description': '(NUM_SH_TERM_FACTOR) ratio of average number of spherical harmonic terms to total possible (NLM)'
        },
        'cell_to_point_ratio': {
            'default_value': cell_to_point_ratio,
            'description': '(CELL_TO_POINT_RATIO) ratio of number of grid cells to grid points'
        },
        'high_order_radiance': {
            'default_value': high_order_radiance,
            'description': 'True to keep the high order radiance field in memory for diagnostic purposes.'
        },
        'ip_flag':{
            'default_value': ip_flag,
            'description':
                ['(IPFLAG) Bit flags for independent pixel mode: 0 for 3D, 1 for independent (2D) scans in X, 2 for 2D scans in Y (X-Z planes)',
                 '3 for indepedent pixels (i.e. bit 0 for X and bit 1 for Y). Bit 2 of IPFLAG means do the direct beam in 3D,',
                 ' e.g. IPFLAG=7 means 3D direct beam but IP diffuse radiative transfer.']
        },
        'iterfixsh': {
            'default_value': iterfixsh,
            'description': 'The number of iterations after which the spherical harmonic truncation is fixed.'
        },
        'tautol': {
            'default_value': tautol,
            'description': 'The maximum optical length of each subgrid interval for integration of radiances. Only affects radiance accuracy.'
        },
        'angle_set': {
            'default_value': angle_set,
            'description': 'The angle set used for the discrete ordinates. Gaussian: 1, Reduced Gaussian: 2, Reduced Double Gaussian: 3.'
        },
        'transcut': {
            'default_value': transcut,
            'description': 'The minimum transmission to integrate to when calculating radiances.'
        },
        'transmin': {
            'default_value': transmin,
            'description': 'The minimum transmission to integrate in the short characteristics scheme. Values less than 1.0 turn 1-cell characteristics into potentially long characteristics.'
        }
    }
    return shdom_parameters

def get_config(config_file_name=None):
    """
    Load a json file of numerical parameters generated by at3d.configuration.make_config()
    as an xr.Dataset for use in solver.RTE

    See shdom.txt or configuration.py for description of numerical parameters.

    Parameters
    ----------
    config_file_name : string
        The name of the .json file to load numerical parameters from. See default_config.json
        for formatting.

    Returns
    -------
    config_dataset : xr.Dataset
        The numerical parameters.
    """
    if config_file_name is None:
        configuration = make_config_data()
        config_dataset = xr.Dataset(
            data_vars={key:attr['default_value'] for key, attr in configuration.items()}
        )
    else:
        with open(config_file_name, 'r') as f:
            configuration = json.load(f)

        config_dataset = xr.Dataset(
            data_vars={key:attr['default_value'] for key, attr in configuration.items()}
        )

    return config_dataset
