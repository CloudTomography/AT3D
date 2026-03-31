#!/usr/bin/env python
"""Generate the f2py signature file (core.pyf) for the AT3D Fortran extension.

This is the single source of truth for which Fortran subroutines are exposed
to Python via the ``at3d.core`` module. Meson calls this script at build time,
but it can also be run standalone for debugging:

    python tools/generate_pyf.py --output-dir src/

When adding a new Fortran subroutine to the Python interface:
  1. Add the subroutine name to F2PY_CORE_API below.
  2. If the subroutine lives in a new Fortran file, add that file to
     F2PY_SHDOM_FILES below.
  3. Rebuild: pip install .
"""
import argparse
import os
import subprocess
import sys

# ---------------------------------------------------------------------------
# Fortran source files (relative to the project root ``src/`` directory).
# Only the polarized SHDOM branch is compiled (unpolarized is unsupported).
# ---------------------------------------------------------------------------
F2PY_SRC_DIR = "src"

F2PY_SHDOM_FILES = [
    # Base files
    "fftpack.f",
    "ocean_brdf.f",
    "shdom_nompi.f",
    "shdomsub5.f",
    "surface.f",
    "util.f90",
    "microwave_gases.f90",
    "ckdfu.f",
    # Polarized SHDOM files
    "polarized/shdom90.f90",
    "polarized/shdomsub1.f",
    "polarized/shdomsub2.f",
    "polarized/shdomsub3.f",
    "polarized/shdomsub4.f",
    "polarized/make_mie_table.f90",
    "polarized/miewig.f",
    "polarized/indexwatice.f",
]

# ---------------------------------------------------------------------------
# Public API – Fortran subroutines/functions exposed to Python via f2py.
# Keep this list sorted by category for readability.
# ---------------------------------------------------------------------------
F2PY_CORE_API = [
    # Mie scattering
    "get_mie_table",
    "get_center_wavelen",
    "get_refract_index",
    "get_nsize",
    "get_sizes",
    "compute_mie_all_sizes",
    "make_multi_size_dist",
    "write_mono_table",
    "read_mono_table",
    "get_poly_table",
    "write_poly_table",
    "read_poly_table",
    "transform_leg_to_phase",
    # Rayleigh
    "rayleigh_extinct",
    "rayleigh_phase_function",
    # MPI / init
    "start_mpi",
    "end_shdom_mpi",
    "check_input_parmeters",
    # Grid
    "new_grids",
    "init_cell_structure",
    "transfer_pa_to_grid",
    # RTE solver
    "init_solution",
    "solution_iterations",
    "make_direct",
    "make_direct_derivative",
    # Radiances
    "render",
    "compute_top_radiances",
    "fixed_lambertian_boundary",
    "variable_lambertian_boundary",
    # Optimization / gradient
    "levisapprox_gradient",
    # Space carving / geometry
    "space_carve",
    # Phase functions
    "precompute_phase_check",
    "precompute_phase_check_grad",
    # Miscellaneous
    "optical_depth",
    "prep_surface",
    "read_properties",
    "compute_netfluxdiv",
    "compute_sh",
    "min_optical_depth",
    "gradient_l2_old",
    "average_subpixel_rays",
    "project",
    "util_integrate_rays",
    "util_locate_point",
    "util_get_interp_kernel2",
    "check_property_input",
    "nearest_binary",
    "cell_average",
    "update_costfunction",
    "output_cell_split",
    "compute_radiance_grid",
    "compute_source_grid",
    "traverse_grid",
    "read_property_size",
    "adjoint_linear_interpolation",
    "get_shadow",
    "transmission_integral",
    "quicksort_new",
    "construct_ptr",
    "ssort",
    "compute_dir_source",
    "eddrtf",
    "phase_function_mixing",
    "prepare_deriv_interps",
    "planck_function",
    "planck_derivative",
    "compute_gradient_oneproppoint",
    "compute_direct_beam_deriv",
    "extinction_derivative_point",
    "interpolate_point",
    "divide_cell",
    "grid_smoothing",
    "ylmall",
    "test_source",
    "wigner_transform",
    "sh_to_do_unpol",
    "surface_brdf",
    "sh_to_do",
    "make_angle_set",
    "surface_parm_interp",
    "do_to_sh",
    "do_to_sh_unpol",
    "make_sh_do_coef",
    "ross_li_thick_sparse",
    "calc_mw_gas_absorption",
    "integrate_thermal_source",
    "ckdfu",
]

MODULE_NAME = "core"


def generate_pyf(project_root, output_dir):
    """Run f2py to generate the .pyf signature file."""
    src_dir = os.path.join(project_root, F2PY_SRC_DIR)
    fortran_files = [
        os.path.join(src_dir, f) for f in F2PY_SHDOM_FILES
    ]
    pyf_path = os.path.join(output_dir, f"{MODULE_NAME}.pyf")

    cmd = [
        sys.executable, "-m", "numpy.f2py",
        "-h", pyf_path,
        "--overwrite-signature",
        "-m", MODULE_NAME,
    ] + fortran_files + [
        "only:",
    ] + F2PY_CORE_API

    print(f"Generating {pyf_path}")
    subprocess.check_call(cmd)
    return pyf_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory to write core.pyf into (default: <project_root>/src/)",
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Project root directory (default: auto-detect from script location)",
    )
    args = parser.parse_args()

    # Auto-detect project root: this script lives in <root>/tools/
    if args.project_root:
        project_root = os.path.abspath(args.project_root)
    else:
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir)
    else:
        output_dir = os.path.join(project_root, F2PY_SRC_DIR)

    os.makedirs(output_dir, exist_ok=True)
    generate_pyf(project_root, output_dir)


if __name__ == "__main__":
    main()
