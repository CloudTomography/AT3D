import setuptools
import popen2
import os
import tempfile
import shlex
import numpy
import sys
import os


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
VERSION = "0.0.1"

classifiers =  ['Development Status :: 3 - Alpha',
                'Programming Language :: Python',
                'License :: OSI Approved :: MIT License',
                'Intended Audience :: Science/Research',
                'Topic :: Scientific/Engineering',
                'Topic :: Scientific/Engineering :: Mathematics',
                'Operating System :: OS Independent']

#
# Set this to True for compiling the parallel
# version of the SHDOM algorithm.
#
PARALLEL_SUPPORT = False

#
# f2py stuff
#
F2PY_CMD = 'f2py'
F2PY_MODULE_NAME = 'core'
F2PY_SRC_PATH = 'src'
F2PY_SIGN_FILE = '{path}/core.pyf'.format(path=F2PY_SRC_PATH)
F2PY_SHDOM_FILES = ['shdom90.f90',
                    'shdomsub1.f',
                    'shdomsub2.f',
                    'shdomsub3.f',
                    'shdomsub4.f',
                    'fftpack.f',
                    'ocean_brdf.f',
                    'make_mie_table.f90',
                    'mieindsub.f']
if PARALLEL_SUPPORT:
    F2PY_SHDOM_FILES += ['shdom_mpi.f'] 
else:
    F2PY_SHDOM_FILES += ['shdom_nompi.f']
F2PY_SHDOM_FILES = [
    
            '{path}/{file_name}'.format(path=F2PY_SRC_PATH, file_name=file_name) for file_name in F2PY_SHDOM_FILES
]

F2PY_CORE_API = [
    'get_mie_table',
    'get_center_wavelen',
    'write_mie_table',
    'read_mie_table',
    'rayleigh_extinct',
    'start_mpi',
    'end_shdom_mpi',
    'shdom_property_arrays',
    'check_input_parmeters',
    'new_grids',
    'init_cell_structure',
    'solve_rte',
    'render',
    'par_render',
    'precompute_phase',
    'ylmall',
    'compute_top_radiances',
    'fixed_lambertian_boundary',
    'variable_lambertian_boundary',
    'ext_gradient'
]


def uniq_arr(arr):
    """Remove repeated values from an array and return new array."""
    ret = []
    for i in arr:
        if i not in ret:
            ret.append(i)
    return ret


def _run_command(cmd):
    out_file, in_file, err_file = popen2.popen3(cmd)
    output = out_file.read() + err_file.read()
    out_file.close()
    in_file.close()
    err_file.close()
    # need this hack to get the exit status
    print 'running ' + cmd
    out_file = os.popen(cmd)
    if out_file.close():
        # close returns exit status of command.
        return ''
    else:
        # no errors, out_file.close() returns None.
        return output


def _get_mpi_cmd():
    """Returns the output of the command used to compile using
    mpicc."""
    # LAM/OPENMPI/MPICH2
    output = _run_command('mpicc -show')
    if output:
        return output

    # FIXME: If appears that MPICH actually needs this these hacks.

    # MPICH
    # works with MPICH version 1.2.1 (on Debian)
    output = _run_command('mpicc -compile_info -link_info')
    if output:
        return output

    # Old version of MPICH needs this hack.
    tmp_base = tempfile.mktemp()
    tmp_c = tmp_base + ".c"
    tmp_o = tmp_base + ".o"
    tmp_file = open(tmp_c, "w")
    tmp_file.write('#include "mpi.h"\nint main(){return 0;}\n')
    tmp_file.close()
    output = _run_command("mpicc -show;"\
                          "mpicc -echo -c %s -o %s"%(tmp_c, tmp_o))
    os.remove(tmp_c)
    if os.path.exists(tmp_o):
        os.remove(tmp_o)
    if output:
        return output
    else:
        return ''


def get_mpi_flags():
    output = _get_mpi_cmd()
    print output
    if not output:
        if sys.platform=='win32': # From Simon Frost
            #this didn't work on my machine (Vladimir Lazunin on April 7, 2009)
            #output = "gcc -L$MPICH_DIR\SDK.gcc\lib -lmpich -I$MPICH_DIR\SDK.gcc\include"

            #"MPICH_DIR" must be set manually in environment variables
            mpi_dir = os.getenv("MPICH_DIR")
            if mpi_dir == None:
                print 'MPICH_DIR environment variable must be set'
                exit()

            #for MPICH2
            sdk_prefix = mpi_dir
            lib_name = 'mpi'

            #for MPICH1
            if os.path.exists(sdk_prefix + '\\SDK'):
                sdk_prefix += '\\SDK'
                lib_name = 'mpich'
            output = 'gcc -L"%(sdk_prefix)s\lib" -l"%(lib_name)s" -I"%(sdk_prefix)s\include"' % {'sdk_prefix' : sdk_prefix, 'lib_name' : lib_name}
        else:
            output = 'cc -L/usr/opt/mpi -lmpi -lelan'


    # Now get the include, library dirs and the libs to link with.
    #flags = string.split(output)
    flags = shlex.split(output)
    flags = uniq_arr(flags) # Remove repeated values.
    inc_dirs = []
    lib_dirs = []
    libs = []
    def_macros = []
    undef_macros = []
    for f in flags:
        if f[:2] == '-I':
            inc_dirs.append(f[2:])
        elif f[:2] == '-L':
            lib_dirs.append(f[2:])
        elif f[:2] == '-l' and f[-1] != "'": # Patched by Michael McKerns July 2009
            libs.append(f[2:])
        elif f[:2] == '-U':
            undef_macros.append(f[2:])
        elif f[:2] == '-D':
            tmp = string.split(f[2:], '=')
            if len(tmp) == 1:
                def_macros.append((tmp[0], None))
            else:
                def_macros.append(tuple(tmp))
    return {'inc_dirs': inc_dirs, 'lib_dirs': lib_dirs, 'libs':libs,
            'def_macros': def_macros, 'undef_macros': undef_macros}


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


    if PARALLEL_SUPPORT:
        mpi_flags = get_mpi_flags()
        mpi_flags['inc_dirs'].append(numpy.get_include())
    
    
        # setting some extra compile flags for 64 bit architectures, utilizing
        # distutils.sysconfig to check which compiler to use
        if os.name == 'posix' and os.uname()[4] == 'x86_64':
            #Extra flags for 64 bit architectures
            extra_compile_args = ['-fPIC']
        else:
            extra_compile_args = None


        config.add_extension(
            name=F2PY_MODULE_NAME,
            sources=[F2PY_SIGN_FILE] + F2PY_SHDOM_FILES,
            include_dirs=mpi_flags['inc_dirs'],
            library_dirs=mpi_flags['lib_dirs'],
            libraries=mpi_flags['libs'],
            define_macros=mpi_flags['def_macros'],
            undef_macros=mpi_flags['undef_macros'],
            extra_compile_args=extra_compile_args,
            f2py_options=[]#'--debug-capi']            
        )
    else:
        config.add_extension(
            name=F2PY_MODULE_NAME,
            sources=[F2PY_SIGN_FILE] + F2PY_SHDOM_FILES,
            f2py_options=[]#s'--debug-capi']
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
        classifiers = classifiers
    )