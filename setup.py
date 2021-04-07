import setuptools
import os
import subprocess

NAME = "pyshdom"
DESCRIPTION = "3D Radiative Transfer Inversion using the SHDOM forward algorithm"
LONG_DESCRIPTION ="""
Pyshdom performs 3D reconstruction of cloud microphysical properties from multi-angle, multi-spectral solar reflected
radiation using a non-linear optimization procedure. The core radiative transfer routines are sourced from the
Fortran SHDOM (Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer) code by Frank K. Evans [1].
The python package was created by Aviad Levis [2], Amit Aides (Technion - Israel Institute of Technology) and Jesse Loveridge (University of Illinois).
The inversion algorithms for can be found in the following papers:
 - `Levis, Aviad, et al. "Airborne Three-Dimensional Cloud Tomography." Proceedings of the IEEE International Conference on Computer Vision. 2015.`
 - `Levis, Aviad, et al. "Multiple-Scattering Microphysical Tomography." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2017.`

[1] http://nit.colorado.edu/shdom.html
[2] https://www.aviadlevis.com/3d-remote-sensing
"""

MAINTAINER = "Aviad Levis; Jesse Loveridge"
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

F2PY_SHDOM_FILES = ['fftpack.f', 'ocean_brdf.f', 'shdom_nompi.f', 'shdomsub5.f', 'surface.f','util.f90']

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
    'levisapprox_gradient',
    'space_carve',
    'precompute_phase_check',
    'precompute_phase_check_grad',
    'optical_depth',
    'prep_surface',
    'read_properties',
    'compute_netfluxdiv',
    'compute_sh',
    'min_optical_depth',
    'gradient_l2_old',
    'average_subpixel_rays',
    'pencil_beam_prop',
    'project',
    'util_integrate_rays',
    'util_locate_point',
    'util_get_interp_kernel2',
    'check_property_input',
    'nearest_binary',
    'cell_average',
    'prepare_diphaseind',
    'update_costfunction',
    'output_cell_split',
    'compute_radiance_grid',
    'compute_source_grid',
    'traverse_grid',
    'read_property_size',
    'adjoint_linear_interpolation',
    'get_shadow'
]

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
        package_path='lib',
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

if __name__ == "__main__":

    from numpy.distutils.core import setup
    createSignatureFile()

    setup(
        configuration = configuration,
        packages = setuptools.find_packages(),
        include_package_data = True,
        platforms = ["any"],
        requires = ["numpy", "scipy"],
        tests_require = ['nose',],
        test_suite = 'nose.collector',
        zip_safe = True,
        classifiers = classifiers
    )
