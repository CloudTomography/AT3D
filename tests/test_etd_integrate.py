import numpy as np
from scipy.integrate import quad
import pytest

# Adjust this import based on your Meson build structure
# e.g., from at3d.core import fqetd
import at3d

# ==============================================================================
# FQETD UNIT TESTS
# ==============================================================================

def test_fqetd_branch_a_flat():
    """Test Branch A: Flat extinction gradient (|eta| < 1e-3)."""
    ext0, ext1 = 1.0, 1.0
    srcext0 = np.array([0.5], dtype=np.float64)
    srcext1 = np.array([0.5], dtype=np.float64)
    nq, so = 1, 0.1
    
    src, transcell, abscell = at3d.core.fqetd(ext0, ext1, srcext0, srcext1, nq, so)
    
    # Expected transmission = exp(-0.1)
    np.testing.assert_allclose(transcell, np.exp(-0.1), rtol=1e-7)
    assert np.all(np.isfinite(src))
    assert src[0] >= 0.0

def test_fqetd_branch_b_positive_gradient():
    """Test Branch B: Extinction increases along the ray (ext0 > ext1)."""
    ext0, ext1 = 2.0, 0.1
    srcext0 = np.array([1.0], dtype=np.float64)
    srcext1 = np.array([0.1], dtype=np.float64)
    nq, so = 1, 1.0
    
    src, transcell, abscell = at3d.core.fqetd(ext0, ext1, srcext0, srcext1, nq, so)
    
    assert np.all(np.isfinite(src))
    assert 0.0 <= transcell <= 1.0

def test_fqetd_branch_c_negative_gradient():
    """Test Branch C: Extinction decreases along the ray (ext0 < ext1)."""
    ext0, ext1 = 0.1, 2.0
    srcext0 = np.array([0.1], dtype=np.float64)
    srcext1 = np.array([1.0], dtype=np.float64)
    nq, so = 1, 1.0
    
    src, transcell, abscell = at3d.core.fqetd(ext0, ext1, srcext0, srcext1, nq, so)
    
    assert np.all(np.isfinite(src))
    assert 0.0 <= transcell <= 1.0

def test_fqetd_stokes_vector_rank1():
    """Test memory stability for multi-channel/Stokes arrays (NQ = 4)."""
    ext0, ext1 = 1.5, 1.2
    # I, Q, U, V components
    srcext0 = np.array([1.0, 0.2, -0.1, 0.0], dtype=np.float64)
    srcext1 = np.array([1.2, 0.1, 0.0, 0.05], dtype=np.float64)
    nq, so = 4, 0.5
    
    src, transcell, abscell = at3d.core.fqetd(ext0, ext1, srcext0, srcext1, nq, so)
    
    assert src.shape == (4,)
    assert np.all(np.isfinite(src))

def test_fqetd_zero_path_length():
    """
    Test the Division-by-Zero edge case. 
    If SO = 0.0, EXT_GRAD divides by zero. SHDOM grids can produce SO=0.
    """
    ext0, ext1 = 1.0, 1.0
    srcext0 = np.array([1.0], dtype=np.float64)
    srcext1 = np.array([1.0], dtype=np.float64)
    nq, so = 1, 0.0
    
    # Ideally, your wrapper or Fortran code should catch this.
    # If this segfaults or raises a FloatingPointError, you need to add an 
    # IF (SO .LE. 1.0D-8) RETURN at the top of the Fortran code.
    try:
        src, transcell, abscell = at3d.core.fqetd(ext0, ext1, srcext0, srcext1, nq, so)
        if np.isnan(src[0]):
            pytest.fail("SO=0.0 resulted in NaN. Needs a guard in Fortran.")
    except Exception as e:
        pytest.fail(f"SO=0.0 caused a crash: {e}")

def test_fqetd_optically_thick():
    """Test extreme optical depths to ensure exponentials do not overflow."""
    ext0, ext1 = 1000.0, 1000.0
    srcext0 = np.array([500.0], dtype=np.float64)
    srcext1 = np.array([500.0], dtype=np.float64)
    nq, so = 1, 10.0  # Tau = 10000
    
    src, transcell, abscell = at3d.core.fqetd(ext0, ext1, srcext0, srcext1, nq, so)
    
    # Transmission should safely underflow to exactly 0.0, not NaN
    assert transcell == 0.0
    assert np.all(np.isfinite(src))

def test_fqetd_negative_extinction():
    """Test non-physical negative extinction (can happen in iterative solvers)."""
    ext0, ext1 = -0.5, -0.2
    srcext0 = np.array([-0.1], dtype=np.float64)
    srcext1 = np.array([-0.05], dtype=np.float64)
    nq, so = 1, 0.5
    
    # We want to ensure it doesn't crash, even if the physics output is garbage
    src, transcell, abscell = at3d.core.fqetd(ext0, ext1, srcext0, srcext1, nq, so)
    assert np.all(np.isfinite(src))

import numpy as np
from scipy.integrate import quad

def python_reference_fqetd(ext0, ext1, srcext0, srcext1, so):
    # Optical depth and transmission
    ext_grad = (ext1 - ext0) / (2.0 * so)
    tau = 0.5 * (ext0 + ext1) * so
    transcell = np.exp(-tau)
    
    # INT_0: Integral of exp(-tau(s -> so))
    def integrand_int0(s):
        tau_s = quad(lambda sp: ext0 + 2.0 * ext_grad * sp, s, so)[0]
        return np.exp(-tau_s)
    
    int_0, _ = quad(integrand_int0, 0, so, epsabs=1e-13, epsrel=1e-13)
    
    # INT_1: Integral of s * exp(-tau(s -> so))
    def integrand_int1(s):
        tau_s = quad(lambda sp: ext0 + 2.0 * ext_grad * sp, s, so)[0]
        return s * np.exp(-tau_s)
        
    int_1, _ = quad(integrand_int1, 0, so, epsabs=1e-13, epsrel=1e-13)
    
    # Linear SRCEXT product assumption matching Fortran
    cqi = (srcext1 - srcext0) / so
    src_ref = srcext0 * int_0 + cqi * int_1
    
    return src_ref, transcell

@pytest.mark.skipif(at3d.core.fqetd is None, reason="Fortran module 'fqetd' not compiled or found.")
@pytest.mark.parametrize("ext0, ext1, so", [
    # --- BRANCH A: Small / Flat Extinction Gradients (|eta| < 1e-3) ---
    (1.0000, 1.0001, 0.05),  # Thin flat
    (2.5000, 2.5002, 0.20),  # Moderate flat
    
    # --- BRANCH B: Increasing Extinction (Ext1 > Ext0) ---
    (0.5000, 3.0000, 0.50),  # Moderate increase
    (0.1000, 8.0000, 0.30),  # Steep increase
    
    # --- BRANCH C: Decreasing Extinction (Ext1 < Ext0) ---
    (3.0000, 0.5000, 0.50),  # Moderate decrease
    (8.0000, 0.1000, 0.30),  # Steep decrease
    
    # --- EXTREME OPTICAL DEPTH SWEEP ---
    (0.0100, 0.0500, 1e-4),  # Optically ultra-thin (tests cancellation guard)
    (10.000, 20.000, 0.50),  # Optically thick
])
def test_fqetd_branches_against_quadrature(ext0, ext1, so):
    """
    Compares the Fast Padé ETD Fortran subroutine against 
    double-precision adaptive Gauss-Kronrod quadrature.
    """
    # Define sample source-extinction values at entry and exit
    srcext0 = 1.25
    srcext1 = 3.50
    
    # 1. Run Fortran FQETD subroutine
    nq = 1
    # Wrap scalars into 1-element arrays for f2py compatibility
    src_in_0 = np.array([srcext0], dtype=np.float64)
    src_in_1 = np.array([srcext1], dtype=np.float64)
    
    src_out, trans_out, abscell = at3d.core.fqetd(ext0, ext1, src_in_0, src_in_1, nq, so)
    
    # 2. Run High-Precision Python Reference
    src_ref, trans_ref = python_reference_fqetd(ext0, ext1, srcext0, srcext1, so)
    
    # 3. Assert agreement (allowing for the small approximation error of Padé/Taylor)
    # Fast Padé ETD has a theoretical max error around 1.5e-7 to 1e-8.
    np.testing.assert_allclose(trans_out, trans_ref, rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(src_out[0], src_ref, rtol=1e-5, atol=1e-6)


@pytest.mark.skipif(at3d.core.fqetd is None, reason="Fortran module 'fqetd' not compiled or found.")
def test_fqetd_multi_stokes_vector():
    """Validates that multi-channel / multi-Stokes vector arrays (NSTOKES > 1) 
    integrate correctly and independently against quadrature references."""
    ext0, ext1 = 1.5, 2.5
    so = 0.4
    nstokes = 3
    
    # Independent profiles for I, Q, U components
    srcext0 = np.array([1.0, 0.2, -0.1], dtype=np.float64)
    srcext1 = np.array([2.0, 0.4, -0.2], dtype=np.float64)
    
    src_out, trans_out, abscell = at3d.core.fqetd(ext0, ext1, srcext0, srcext1, nstokes, so)
    
    for i in range(nstokes):
        src_ref, trans_ref = python_reference_fqetd(ext0, ext1, srcext0[i], srcext1[i], so)
        np.testing.assert_allclose(trans_out, trans_ref, rtol=1e-6)
        np.testing.assert_allclose(src_out[i], src_ref, rtol=1e-5)