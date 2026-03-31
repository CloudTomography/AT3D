# --- Python source files ---
python_sources = [
  '__init__.py',
  'aerosol.py',
  'callback.py',
  'checks.py',
  'configuration.py',
  'containers.py',
  'exceptions.py',
  'gas_absorption.py',
  'generate.py',
  'gradient.py',
  'grid.py',
  'initialization.py',
  'initialization_workflow.py',
  'medium.py',
  'mie.py',
  'optimize.py',
  'parallel.py',
  'preprocessing.py',
  'rayleigh.py',
  'regularization.py',
  'sensor.py',
  'size_distribution.py',
  'solver.py',
  'source.py',
  'space_carve.py',
  'stereo.py',
  'surface.py',
  'transforms.py',
  'uncertainties.py',
  'util.py',
  'visualization.py',
]

py.install_sources(python_sources, subdir: 'at3d')

# --- Fortran source files (polarized SHDOM) ---
fortran_sources = files(
  # Base files
  '../src/fftpack.f',
  '../src/ocean_brdf.f',
  '../src/shdom_nompi.f',
  '../src/shdomsub5.f',
  '../src/surface.f',
  '../src/util.f90',
  '../src/microwave_gases.f90',
  '../src/ckdfu.f',
  # Polarized SHDOM files
  '../src/polarized/shdom90.f90',
  '../src/polarized/shdomsub1.f',
  '../src/polarized/shdomsub2.f',
  '../src/polarized/shdomsub3.f',
  '../src/polarized/shdomsub4.f',
  '../src/polarized/make_mie_table.f90',
  '../src/polarized/miewig.f',
  '../src/polarized/indexwatice.f',
)

# --- f2py wrapper generation ---
# Step A: Generate the .pyf signature file from Fortran sources.
# tools/generate_pyf.py is the single source of truth for the API list.
generate_pyf = custom_target(
  'core.pyf',
  input: fortran_sources,
  output: 'core.pyf',
  command: [
    py, meson.project_source_root() / 'tools' / 'generate_pyf.py',
    '--project-root', meson.project_source_root(),
    '--output-dir', '@OUTDIR@',
  ],
)

# Step B: Generate C wrapper and Fortran wrappers from the .pyf file.
core_wrappers = custom_target(
  'coremodule.c',
  input: generate_pyf,
  output: ['coremodule.c', 'core-f2pywrappers.f'],
  command: [py, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@'],
)

# --- Build the core extension module ---
py.extension_module(
  'core',
  [core_wrappers, fortran_sources, fortranobject_c],
  include_directories: inc_np,
  dependencies: py_dep,
  subdir: 'at3d',
  install: true,
)

# --- Package data ---
subdir('data')
