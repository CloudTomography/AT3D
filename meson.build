project(
  'at3d',
  'c', 'fortran',
  version: '4.1.3',
  license: 'GPL-3.0-or-later',
  meson_version: '>=1.2.0',
  default_options: [
    'buildtype=release',
  ],
)

# --- Fortran compiler flags ---
# SHDOM's legacy F77 code passes arguments with mismatched types (e.g. REAL
# to INTEGER working arrays). GCC 10+ treats this as an error by default.
# The -fallow-argument-mismatch flag downgrades it back to a warning.
# It is harmless on F90 code that does not have mismatches.
fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc' and fc.version().version_compare('>=10')
  add_project_arguments('-fallow-argument-mismatch', language: 'fortran')
endif

# --- Python / NumPy / f2py setup ---
py = import('python').find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true,
).stdout().strip()

incdir_f2py = run_command(py,
  ['-c', 'import numpy.f2py; print(numpy.f2py.get_include())'],
  check: true,
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

# --- Build the at3d package ---
subdir('at3d')
