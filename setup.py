import setuptools
import os
import subprocess

NAME = "shdom"
EXTENSION_NAME = "pyshdom"
DESCRIPTION = "3D Radiative Transfer Inversion using the SHDOM forward algorithm"
LONG_DESCRIPTION ="""
The Spherical Harmonic Discrete Ordinate Method (SHDOM) for Atmospheric Radiative Transfer
Developed by Frank Evans (http://nit.colorado.edu/shdom.html).

The python wrapper for SHDOM was developed by Amit Aides.

The inversion algorithm for 3D ratiative transfer was developed and contributed to by:
Aviad Levis, Yoav Schechner, Amit Aides, Anthony B. Davis.
(https://www.aviadlevis.com/3d-remote-sensing).
"""

MAINTAINER = "Aviad Levis"
MAINTAINER_EMAIL = "aviad.levis@gmail.com"
URL = "https://github.com/aviadlevis/pyshdom"
LICENSE = "MIT"
VERSION = "3.0.0"

classifiers =  ['Development Status :: 3 - Alpha',
                'Programming Language :: Python',
                'License :: OSI Approved :: MIT License',
                'Intended Audience :: Science/Research',
                'Topic :: Scientific/Engineering',
                'Topic :: Scientific/Engineering :: Mathematics',
                'Operating System :: OS Independent']


#
# Set this to True for compiling the polarized 
# version of the SHDOM algorithm.
#
POLARIZED_SHDOM = True

#
# f2py stuff
#
F2PY_CMD = 'f2py'
F2PY_MODULE_NAME = 'core'
F2PY_SRC_PATH = 'src'
F2PY_SIGN_FILE = '{path}/core.pyf'.format(path=F2PY_SRC_PATH)

F2PY_SHDOM_FILES = ['fftpack.f', 'ocean_brdf.f', 'shdom_nompi.f', 'shdomsub5.f']

if POLARIZED_SHDOM:
    F2PY_SHDOM_FILES.extend(['polarized/shdom90.f90',
                             'polarized/shdomsub1.f',
                             'polarized/shdomsub2.f',
                             'polarized/shdomsub3.f',
                             'polarized/shdomsub4.f',
                             'polarized/make_mie_table.f90', 
                             'polarized/miewig.f',
                             'polarized/indexwatice.f'])
    
else:
    F2PY_SHDOM_FILES.extend(['unpolarized/shdom90.f90',
                             'unpolarized/shdomsub1.f',
                             'unpolarized/shdomsub2.f',
                             'unpolarized/shdomsub3.f',
                             'unpolarized/shdomsub4.f',
                             'unpolarized/make_mie_table.f90', 
                             'unpolarized/mieindsub.f'])

F2PY_SHDOM_FILES = [
            '{path}/{file_name}'.format(path=F2PY_SRC_PATH, file_name=file_name) for file_name in F2PY_SHDOM_FILES
]
  
F2PY_CORE_API = [
    'get_mie_table',
    'get_center_wavelen',
    'get_refract_index',
    'get_nsize',
    'get_sizes',
    'compute_mie_all_sizes',
    'make_multi_size_dist',
    'write_mono_table',
    'read_mono_table',
    'get_poly_table',
    'write_poly_table',
    'read_poly_table',
    'transform_leg_to_phase',
    'rayleigh_extinct',
    'rayleigh_phase_function',
    'start_mpi',
    'end_shdom_mpi',
    'check_input_parmeters',
    'new_grids',
    'init_cell_structure',
    'transfer_pa_to_grid',
    'init_solution',
    'solution_iterations',
    'make_direct',
    'make_direct_derivative',
    'render',
    'compute_top_radiances',
    'fixed_lambertian_boundary',
    'variable_lambertian_boundary',
    'gradient_normcorr',
    'gradient_l2',
    'space_carve',
    'precompute_phase_check'
]


#if POLARIZED_SHDOM:
    #F2PY_CORE_API.extend([])


def _run_command(cmd):
    proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out_file = proc.stdout.read()
    err_file = proc.stderr.read()                                                                     
    output = out_file + err_file
    proc.stdout.close()
    proc.stderr.close()
    # need this hack to get the exit status
    print('running ' + cmd)
    out_file = os.popen(cmd)
    if out_file.close():
        # close returns exit status of command.
        return ''
    else:
        # no errors, out_file.close() returns None.
        return output



def createSignatureFile():
    """Create the signature file for the f2py file."""
    signature_cmd = "{cmd} -h {sign} --overwrite-signature -m {module} {files} only: {api}".format(
        cmd=F2PY_CMD,
        sign=F2PY_SIGN_FILE,
        module=F2PY_MODULE_NAME,
        files=' '.join(F2PY_SHDOM_FILES),
        api=' '.join(F2PY_CORE_API)
    )
    _run_command(signature_cmd)
    

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(
        NAME,
        parent_package,
        top_path,
        version = VERSION,
        maintainer  = MAINTAINER,
        maintainer_email = MAINTAINER_EMAIL,
        description = DESCRIPTION,
        license = LICENSE,
        url = URL,
        long_description = LONG_DESCRIPTION
    )

    config.add_extension(
        name=F2PY_MODULE_NAME,
        sources=[F2PY_SIGN_FILE] + F2PY_SHDOM_FILES,
        f2py_options=[]
    )

    return config




def multiple_configurations(name):
    return lambda parent_package='', top_path=None: configuration(name, parent_package,top_path)

if __name__ == "__main__":

    from numpy.distutils.core import setup
    
    createSignatureFile()
    
    setup(
        configuration = configuration,
        packages = setuptools.find_packages(),
        include_package_data = True,
        platforms = ["any"],
        requires = ["numpy"],
        tests_require = ['nose',],
        test_suite = 'nose.collector',
        zip_safe = True,
        install_requires=[
            'tensorboardX'
        ],        
        classifiers = classifiers
    )