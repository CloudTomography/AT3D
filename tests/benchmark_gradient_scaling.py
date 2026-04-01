"""Benchmark: gradient scaling with grid size.

Measures solve and gradient wall-clock time for the
SolarJacobianThinNoSurfaceExtinction configuration at multiple grid
resolutions and parallelism levels.  Results are appended to
``tests/benchmark_results.csv`` keyed by git commit hash.

Run with::

    pytest tests/benchmark_gradient_scaling.py -m benchmark -s

Skip during normal CI::

    pytest tests/ -m "not benchmark"
"""
import csv
import datetime
import os
import subprocess
import time
import warnings

import numpy as np
import pytest

import at3d
import at3d.parallel

from test_derivatives import cloud_solar

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------
CSV_PATH = os.path.join(os.path.dirname(__file__), "benchmark_results.csv")
CSV_COLUMNS = [
    "commit",
    "timestamp",
    "method",
    "resolution",
    "base_grid_pts",
    "actual_npts",
    "n_jobs",
    "solve_time_s",
    "gradient_time_s",
]


def _git_commit_hash():
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
    except Exception:
        return "unknown"


def append_to_csv(row: dict):
    write_header = not os.path.exists(CSV_PATH) or os.path.getsize(CSV_PATH) == 0
    with open(CSV_PATH, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


# ---------------------------------------------------------------------------
# Mie table (cached across parametrize calls)
# ---------------------------------------------------------------------------
_mie_cache = {}


def _get_mie_table():
    if "table" not in _mie_cache:
        table_path = os.path.join(os.path.dirname(__file__), "data", "mie_table_860nm.nc")
        if os.path.exists(table_path):
            import xarray as xr
            _mie_cache["table"] = xr.open_dataset(table_path)
        else:
            _mie_cache["table"] = at3d.mie.get_mono_table(
                "Water",
                (0.86, 0.86),
                max_integration_radius=65.0,
                minimum_effective_radius=0.1,
                relative_dir="./data",
                verbose=False,
            )
            _mie_cache["table"].to_netcdf(table_path)
    return _mie_cache["table"]


# ---------------------------------------------------------------------------
# Benchmark core
# ---------------------------------------------------------------------------
def run_benchmark(resolution, n_jobs, method="jacobian"):
    mie_mono_table = _get_mie_table()

    # ---- Setup (untimed): build solvers, sensors, solve RTE, get measurements ----
    solvers, Sensordict, cloud_poly_tables, _, rte_grid = cloud_solar(
        mie_mono_table,
        ext=0.1,
        veff=0.1,
        reff=10.0,
        ssalb=1.0,
        solarmu=1.0,
        surfacealb=0.0,
        ground_temperature=200.0,
        step=0.0,
        nmu=16,
        split=0.0,
        resolution=resolution,
        random=True,
        random_ssalb=True,
        perturb="extinct",
        deltam=True,
    )

    base_grid_pts = int(rte_grid.x.size * rte_grid.y.size * rte_grid.z.size)
    solver = solvers[0.86]
    indices_for_jacobian = np.where(solver.medium["cloud"].extinction.data > 0.0)

    # Uncertainty model + forward sensors (matches test_derivatives setup)
    Sensordict.add_uncertainty_model("MISR", at3d.uncertainties.NullUncertainty("L2"))
    for sensor in Sensordict["MISR"]["sensor_list"]:
        Sensordict["MISR"]["uncertainty_model"].calculate_uncertainties(sensor)
    forward_sensors = Sensordict.make_forward_sensors()

    # Unknown scatterers
    deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, "cloud", 0.86)
    unknown_scatterers = at3d.containers.UnknownScatterers(
        at3d.medium.UnknownScatterer(deriv_gen, "extinction")
    )

    # Build gradient object (needed for its levis_approximation_grad method)
    grad_kwargs = {
        "exact_single_scatter": True,
        "cost_function": "L2",
    }
    if method == "jacobian":
        grad_kwargs["indices_for_jacobian"] = indices_for_jacobian

    gradient_obj = at3d.gradient.LevisApproxGradientUncorrelated(
        Sensordict,
        solvers,
        forward_sensors,
        unknown_scatterers,
        parallel_solve_kwargs={
            "maxiter": 200,
            "n_jobs": n_jobs,
            "setup_grid": False,
            "verbose": False,
            "init_solution": False,
        },
        gradient_kwargs=grad_kwargs,
        uncertainty_kwargs={"add_noise": False},
    )

    # ---- Timed block 1: Solve (re-initialise so timing is meaningful) ----
    t0 = time.perf_counter()
    solvers.solve(
        maxiter=200,
        n_jobs=n_jobs,
        setup_grid=False,
        verbose=False,
        init_solution=True,
    )
    solve_time = time.perf_counter() - t0

    actual_npts = int(solvers[0.86]._npts)

    # ---- Timed block 2: Gradient (partial derivs + ray integration) ----
    t0 = time.perf_counter()
    solvers.calculate_microphysical_partial_derivatives(unknown_scatterers)
    solvers.calculate_direct_beam_derivative()
    rte_sensors, sensor_mapping = forward_sensors.sort_sensors(solvers, Sensordict)
    par_grad_kwargs = {
        "exact_single_scatter": True,
        "cost_function": "L2",
    }
    if method == "jacobian":
        par_grad_kwargs["indices_for_jacobian"] = indices_for_jacobian
    at3d.parallel.parallel_gradient(
        solvers,
        rte_sensors,
        sensor_mapping,
        forward_sensors,
        gradient_fun=gradient_obj.levis_approximation_grad,
        n_jobs=n_jobs,
        **par_grad_kwargs,
    )
    gradient_time = time.perf_counter() - t0

    return {
        "method": method,
        "resolution": resolution,
        "base_grid_pts": base_grid_pts,
        "actual_npts": actual_npts,
        "n_jobs": n_jobs,
        "solve_time_s": f"{solve_time:.3f}",
        "gradient_time_s": f"{gradient_time:.3f}",
    }


# ---------------------------------------------------------------------------
# Pytest entry points
# ---------------------------------------------------------------------------
@pytest.mark.benchmark
@pytest.mark.parametrize(
    "resolution,n_jobs,method",
    [
        (1, 1, "jacobian"),
        (1, 4, "jacobian"),
        (2, 1, "jacobian"),
        (2, 4, "jacobian"),
        (4, 1, "jacobian"),
        (4, 4, "jacobian"),
        (1, 1, "adjoint"),
        (1, 4, "adjoint"),
        (2, 1, "adjoint"),
        (2, 4, "adjoint"),
        (4, 1, "adjoint"),
        (4, 4, "adjoint"),
    ],
    ids=lambda val: str(val),
)
def test_gradient_scaling(resolution, n_jobs, method):
    result = run_benchmark(resolution, n_jobs, method=method)
    result["commit"] = _git_commit_hash()
    result["timestamp"] = datetime.datetime.now().isoformat()
    append_to_csv(result)
    print(
        f"\n  method={method}  resolution={resolution}  n_jobs={n_jobs}  "
        f"base_pts={result['base_grid_pts']}  actual_npts={result['actual_npts']}  "
        f"solve={result['solve_time_s']}s  gradient={result['gradient_time_s']}s"
    )
