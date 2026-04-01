"""Quick test verifying adjoint gradient path matches Jacobian path."""
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from test_derivatives import cloud_solar
import at3d
import xarray as xr

def test_adjoint_vs_jacobian():
    """Gradient without Jacobian (adjoint path) must match gradient
    with Jacobian (old path) to within floating-point tolerance."""
    table_path = os.path.join(os.path.dirname(__file__), "data", "mie_table_860nm.nc")
    if os.path.exists(table_path):
        mie = xr.open_dataset(table_path)
    else:
        mie = at3d.mie.get_mono_table(
            "Water", (0.86, 0.86),
            max_integration_radius=65.0,
            minimum_effective_radius=0.1,
            relative_dir="./data", verbose=False,
        )
        mie.to_netcdf(table_path)

    # Create scene and get measurements
    solvers, Sensordict, cloud_poly, _, rte_grid = cloud_solar(
        mie, ext=10.0, veff=0.1, reff=10.0, ssalb=1.0,
        solarmu=1.0, surfacealb=0.0, ground_temperature=200.0,
        step=0.0, nmu=16, split=0.0, resolution=1,
        random=True, random_ssalb=True, perturb="extinct", deltam=True,
    )
    solver = solvers[0.86]
    idx = np.where(solver.medium['cloud'].extinction.data > 0.0)

    # Perturb the rendered 'I' values (used as measurements) to create a
    # non-trivial cost and gradient. This simulates a real inversion scenario.
    for sensor in Sensordict['MISR']['sensor_list']:
        sensor['I'].values[:] *= 1.1  # +10% bias in "measurements"

    Sensordict.add_uncertainty_model('MISR',
        at3d.uncertainties.NullUncertainty('L2'))
    for sensor in Sensordict['MISR']['sensor_list']:
        Sensordict['MISR']['uncertainty_model'].calculate_uncertainties(sensor)

    deriv_gen = at3d.medium.GridToOpticalProperties(rte_grid, 'cloud', 0.86)
    unknown = at3d.containers.UnknownScatterers(
        at3d.medium.UnknownScatterer(deriv_gen, 'extinction'))

    psk = {'maxiter': 200, 'n_jobs': 1, 'setup_grid': False,
           'verbose': False, 'init_solution': False}

    # ---- Jacobian path (MAKEJACOBIAN=True) ----
    forward_sensors_jac = Sensordict.make_forward_sensors()
    grad_jac_obj = at3d.gradient.LevisApproxGradientUncorrelated(
        Sensordict, solvers, forward_sensors_jac, unknown,
        parallel_solve_kwargs=psk,
        gradient_kwargs={'exact_single_scatter': True,
                         'cost_function': 'L2',
                         'indices_for_jacobian': idx},
        uncertainty_kwargs={'add_noise': False})

    loss_jac, grad_jac, jac = grad_jac_obj()
    print(f'Jacobian path: loss={loss_jac:.6e}')
    print(f'  gradient norm: {np.linalg.norm(grad_jac["gradient"].data):.6e}')

    # ---- Adjoint path (MAKEJACOBIAN=False) ----
    forward_sensors_adj = Sensordict.make_forward_sensors()
    grad_adj_obj = at3d.gradient.LevisApproxGradientUncorrelated(
        Sensordict, solvers, forward_sensors_adj, unknown,
        parallel_solve_kwargs=psk,
        gradient_kwargs={'exact_single_scatter': True,
                         'cost_function': 'L2'},
        uncertainty_kwargs={'add_noise': False})

    loss_adj, grad_adj, _ = grad_adj_obj()
    print(f'Adjoint path:  loss={loss_adj:.6e}')
    print(f'  gradient norm: {np.linalg.norm(grad_adj["gradient"].data):.6e}')

    # ---- Compare loss ----
    cost_diff = abs(loss_jac - loss_adj) / (abs(loss_jac) + 1e-30)
    print(f'Cost reldiff: {cost_diff:.2e}')
    assert cost_diff < 1e-4, f'Cost mismatch: reldiff={cost_diff}'

    # ---- Compare gradient datasets ----
    d1 = grad_jac['gradient'].data
    d2 = grad_adj['gradient'].data
    maxdiff = np.max(np.abs(d1 - d2))
    norm = np.max(np.abs(d1)) + 1e-30
    reldiff = maxdiff / norm
    print(f'gradient: maxdiff={maxdiff:.2e} reldiff={reldiff:.2e} norm={norm:.2e}')
    assert np.allclose(d1, d2, rtol=1e-4, atol=1e-10), \
        f'Gradient mismatch: maxdiff={maxdiff} reldiff={reldiff}'

    print('PASS: adjoint and Jacobian paths produce matching gradients')

if __name__ == '__main__':
    test_adjoint_vs_jacobian()
