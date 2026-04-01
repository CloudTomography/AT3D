"""Tests verifying the adjoint gradient path matches the Jacobian path.

Covers all derivative types (extinction, ssalb, legendre_0_1),
exact_single_scatter on/off, delta-M on/off, thermal source,
and L2/LL cost functions.
"""
import numpy as np
import sys
import os
import warnings

import pytest
import xarray as xr

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from test_derivatives import cloud_solar, cloud

import at3d

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Mie table caching
# ---------------------------------------------------------------------------
_mie_cache = {}


def _get_mie(wavelength_range, cache_name, cache_file):
    if cache_name not in _mie_cache:
        table_path = os.path.join(os.path.dirname(__file__), "data", cache_file)
        if os.path.exists(table_path):
            _mie_cache[cache_name] = xr.load_dataset(table_path)
        else:
            _mie_cache[cache_name] = at3d.mie.get_mono_table(
                "Water", wavelength_range,
                max_integration_radius=65.0,
                minimum_effective_radius=0.1,
                relative_dir="./data", verbose=False,
            )
            _mie_cache[cache_name].to_netcdf(table_path)
    return _mie_cache[cache_name]


def _get_solar_mie():
    return _get_mie((0.86, 0.86), "solar", "mie_table_860nm.nc")


def _get_thermal_mie():
    return _get_mie((11.0, 11.0), "thermal", "mie_table_11micron.nc")


# ---------------------------------------------------------------------------
# Helper: run both paths and compare
# ---------------------------------------------------------------------------
def _compare_adjoint_vs_jacobian(solvers, Sensordict, rte_grid,
                                 wavelength, unknown_variable,
                                 cost_function, exact_single_scatter,
                                 rtol=1e-4, atol=1e-10, label=""):
    """Run Jacobian and adjoint paths, compare cost and gradient."""
    solver = solvers[wavelength]
    idx = np.where(solver.medium['cloud'].extinction.data > 0.0)

    # Perturb measurements so cost/gradient are non-trivial
    for sensor in Sensordict['MISR']['sensor_list']:
        if cost_function == 'LL':
            # LL cost uses log(radiance); ensure all measurements are positive
            sensor['I'].values[:] = np.clip(sensor['I'].values, 1e-6, None)
        sensor['I'].values[:] *= 1.1

    Sensordict.add_uncertainty_model(
        'MISR', at3d.uncertainties.NullUncertainty(cost_function))
    for sensor in Sensordict['MISR']['sensor_list']:
        Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)

    deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', wavelength)
    unknown = at3d.containers.UnknownScatterers(
        at3d.medium.UnknownScatterer(deriv_gen, unknown_variable))

    psk = {'maxiter': 200, 'n_jobs': 1, 'setup_grid': False,
           'verbose': False, 'init_solution': False}

    # ---- Jacobian path ----
    forward_jac = Sensordict.make_forward_sensors()
    grad_jac_obj = at3d.gradient.LevisApproxGradientUncorrelated(
        Sensordict, solvers, forward_jac, unknown,
        parallel_solve_kwargs=psk,
        gradient_kwargs={'exact_single_scatter': exact_single_scatter,
                         'cost_function': cost_function,
                         'indices_for_jacobian': idx},
        uncertainty_kwargs={'add_noise': False})
    loss_jac, grad_jac, _ = grad_jac_obj()

    # ---- Adjoint path ----
    forward_adj = Sensordict.make_forward_sensors()
    grad_adj_obj = at3d.gradient.LevisApproxGradientUncorrelated(
        Sensordict, solvers, forward_adj, unknown,
        parallel_solve_kwargs=psk,
        gradient_kwargs={'exact_single_scatter': exact_single_scatter,
                         'cost_function': cost_function},
        uncertainty_kwargs={'add_noise': False})
    loss_adj, grad_adj, _ = grad_adj_obj()

    # ---- Assertions ----
    d1 = grad_jac['gradient'].data
    d2 = grad_adj['gradient'].data
    grad_norm = np.max(np.abs(d1)) + 1e-30
    maxdiff = np.max(np.abs(d1 - d2))
    reldiff = maxdiff / grad_norm

    cost_rd = abs(loss_jac - loss_adj) / (abs(loss_jac) + 1e-30)

    print(f"\n  [{label}]")
    print(f"  loss: jac={loss_jac:.4e}  adj={loss_adj:.4e}  reldiff={cost_rd:.2e}")
    print(f"  grad: norm={grad_norm:.4e}  maxdiff={maxdiff:.2e}  reldiff={reldiff:.2e}")

    assert grad_norm > 1e-20, \
        f"[{label}] Gradient is effectively zero (norm={grad_norm:.2e}); test is not meaningful"
    assert cost_rd < rtol, \
        f"[{label}] Cost mismatch: reldiff={cost_rd:.2e}"
    assert np.allclose(d1, d2, rtol=rtol, atol=atol), \
        f"[{label}] Gradient mismatch: maxdiff={maxdiff:.2e}  reldiff={reldiff:.2e}"


# ---------------------------------------------------------------------------
# Solar scene fixture
# ---------------------------------------------------------------------------
def _make_solar_scene(perturb, deltam):
    mie = _get_solar_mie()
    return cloud_solar(
        mie, ext=10.0, veff=0.1, reff=10.0, ssalb=1.0,
        solarmu=1.0, surfacealb=0.0, ground_temperature=200.0,
        step=0.0, nmu=16, split=0.0, resolution=1,
        random=True, random_ssalb=True, perturb=perturb, deltam=deltam,
    )


# ===========================================================================
# Solar / extinction derivative tests
# ===========================================================================
class TestAdjointSolarExtinction:

    def test_deltam_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("extinct", deltam=True)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "extinction",
            cost_function="L2", exact_single_scatter=True,
            label="solar/extinct/deltam/exact_ss")

    def test_deltam_no_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("extinct", deltam=True)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "extinction",
            cost_function="L2", exact_single_scatter=False,
            label="solar/extinct/deltam/no_exact_ss")

    def test_nodeltam_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("extinct", deltam=False)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "extinction",
            cost_function="L2", exact_single_scatter=True,
            label="solar/extinct/nodeltam/exact_ss")

    def test_nodeltam_no_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("extinct", deltam=False)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "extinction",
            cost_function="L2", exact_single_scatter=False,
            label="solar/extinct/nodeltam/no_exact_ss")

    def test_cost_LL(self):
        # LL cost uses log(radiance); use surface albedo to ensure
        # all pixels have positive radiance and measurements.
        mie = _get_solar_mie()
        solvers, sd, _, _, grid = cloud_solar(
            mie, ext=10.0, veff=0.1, reff=10.0, ssalb=1.0,
            solarmu=1.0, surfacealb=0.05, ground_temperature=200.0,
            step=0.0, nmu=16, split=0.0, resolution=1,
            random=True, random_ssalb=True, perturb="extinct", deltam=True,
        )
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "extinction",
            cost_function="LL", exact_single_scatter=True,
            label="solar/extinct/deltam/LL")


# ===========================================================================
# Solar / ssalb derivative tests
# ===========================================================================
class TestAdjointSolarAlbedo:

    def test_deltam_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("ssalb", deltam=True)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "ssalb",
            cost_function="L2", exact_single_scatter=True,
            label="solar/ssalb/deltam/exact_ss")

    def test_deltam_no_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("ssalb", deltam=True)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "ssalb",
            cost_function="L2", exact_single_scatter=False,
            label="solar/ssalb/deltam/no_exact_ss")

    def test_nodeltam_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("ssalb", deltam=False)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "ssalb",
            cost_function="L2", exact_single_scatter=True,
            label="solar/ssalb/nodeltam/exact_ss")


# ===========================================================================
# Solar / asymmetry (legendre_0_1) derivative tests
# ===========================================================================
class TestAdjointSolarAsymmetry:

    def test_deltam_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("g", deltam=True)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "legendre_0_1",
            cost_function="L2", exact_single_scatter=True,
            label="solar/g/deltam/exact_ss")

    def test_deltam_no_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("g", deltam=True)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "legendre_0_1",
            cost_function="L2", exact_single_scatter=False,
            label="solar/g/deltam/no_exact_ss")

    def test_nodeltam_exact_ss(self):
        solvers, sd, _, _, grid = _make_solar_scene("g", deltam=False)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 0.86, "legendre_0_1",
            cost_function="L2", exact_single_scatter=True,
            label="solar/g/nodeltam/exact_ss")


# ===========================================================================
# Thermal derivative tests
# ===========================================================================
class TestAdjointThermal:

    @staticmethod
    def _make_thermal_scene(ground_temperature):
        np.random.seed(1)
        mie = _get_thermal_mie()
        return cloud(
            mie, ext=10.0, veff=0.1, reff=10.0, ssalb=0.0,
            solarmu=1.0, surfacealb=0.0,
            ground_temperature=ground_temperature,
            step=0.0, nmu=2, split=0.5, resolution=1,
        )

    def test_no_surface(self):
        solvers, sd, _, _, grid = self._make_thermal_scene(ground_temperature=0.0)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 11.0, "extinction",
            cost_function="L2", exact_single_scatter=True,
            label="thermal/extinct/no_surface")

    def test_with_surface(self):
        solvers, sd, _, _, grid = self._make_thermal_scene(ground_temperature=200.0)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 11.0, "extinction",
            cost_function="L2", exact_single_scatter=True,
            label="thermal/extinct/with_surface")

    def test_no_exact_ss(self):
        solvers, sd, _, _, grid = self._make_thermal_scene(ground_temperature=0.0)
        _compare_adjoint_vs_jacobian(
            solvers, sd, grid, 11.0, "extinction",
            cost_function="L2", exact_single_scatter=False,
            label="thermal/extinct/no_exact_ss")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
