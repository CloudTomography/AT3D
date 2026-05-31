"""Tests for COMPUTE_GRAD_INNER_PRODUCT and helper subroutines.

Three-tier test suite:
  Phase 1: Unit tests for BUILD_TRILIN_MASS_MATRIX and EVAL_FIELD_AT_DIR
  Phase 2: Synthetic tests with small grids and analytic references
  Phase 4: Edge cases / regression
"""
import numpy as np
import pytest

import at3d.core


# ============================================================================
# Helpers
# ============================================================================

def _nlm(ml, mm):
    """Number of SH modes for given ML, MM."""
    n = 0
    for l in range(ml + 1):
        n += 2 * min(l, mm) + 1
    return n


def _lofj(ml, mm):
    """Build degree-L array indexed by sequential SH index J."""
    nlm = _nlm(ml, mm)
    lof = np.zeros(nlm, dtype=np.int32)
    j = 0
    for l in range(ml + 1):
        for m in range(-min(l, mm), min(l, mm) + 1):
            lof[j] = l
            j += 1
    return lof


def _sign_adj(ml, mm):
    """(-1)^L for each SH index."""
    lof = _lofj(ml, mm)
    return (-1.0) ** lof


def _ylmall(transpose, mu, phi, ml, mm, nstleg):
    """Call YLMALL with proper f2py convention. Returns yr(nstleg, nlm)."""
    nlm = _nlm(ml, mm)
    yr = np.zeros((nstleg, nlm), dtype=np.float32, order='F')
    yr = at3d.core.ylmall(int(transpose), float(mu), float(phi), ml, mm, yr)
    return yr


def _make_single_cell_grid(dx=1.0, dy=1.0, dz=1.0):
    """Create a single-cell grid with 8 corner points.

    Returns gridptr(8,1), treeptr(2,1), gridpos(3,8), cellflags(1).
    Node ordering: IOCT = 1 + BITX + 2*BITY + 4*BITZ.
    """
    gridptr = np.zeros((8, 1), dtype=np.int32, order='F')
    for i in range(8):
        gridptr[i, 0] = i + 1  # 1-based

    treeptr = np.zeros((2, 1), dtype=np.int32, order='F')
    treeptr[0, 0] = 0  # parent (unused)
    treeptr[1, 0] = 0  # child pointer = 0 → leaf cell

    gridpos = np.zeros((3, 8), dtype=np.float32, order='F')
    for i in range(8):
        ix = (i) % 2
        iy = (i // 2) % 2
        iz = (i) // 4
        gridpos[0, i] = ix * dx
        gridpos[1, i] = iy * dy
        gridpos[2, i] = iz * dz

    cellflags = np.zeros(1, dtype=np.int16)
    return gridptr, treeptr, gridpos, cellflags


def _make_uniform_sh_ptrs(npts, nlm):
    """Create uniform SH pointer arrays (same truncation at every point).

    Returns ptrs(npts+1) with ptrs[i] = i*nlm (0-based offsets).
    """
    ptrs = np.zeros(npts + 1, dtype=np.int32)
    for i in range(npts + 1):
        ptrs[i] = i * nlm
    return ptrs


def _build_mass_matrix(dx, dy, dz):
    """Call BUILD_TRILIN_MASS_MATRIX. Returns m(8,8)."""
    return at3d.core.build_trilin_mass_matrix(float(dx), float(dy), float(dz))


def _eval_field(coeffs, ylmdir, istokes):
    """Call EVAL_FIELD_AT_DIR. Returns scalar."""
    return at3d.core.eval_field_at_dir(
        np.asfortranarray(coeffs.astype(np.float32)),
        np.asfortranarray(ylmdir.astype(np.float32)),
        int(istokes),
    )


# ============================================================================
# Phase 1: Unit tests for helpers
# ============================================================================

class TestBuildTrilinMassMatrix:
    """Tests for BUILD_TRILIN_MASS_MATRIX."""

    def test_unit_cube_sum(self):
        """Sum of all entries equals cell volume (partition of unity)."""
        m = _build_mass_matrix(1.0, 1.0, 1.0)
        assert np.isclose(m.sum(), 1.0, rtol=1e-6)

    def test_unit_cube_symmetry(self):
        """Mass matrix is symmetric."""
        m = _build_mass_matrix(1.0, 1.0, 1.0)
        assert np.allclose(m, m.T, atol=1e-7)

    def test_unit_cube_row_sums(self):
        """Each row sums to 1/8 (integral of one hat function over cell)."""
        m = _build_mass_matrix(1.0, 1.0, 1.0)
        row_sums = m.sum(axis=1)
        assert np.allclose(row_sums, 1.0 / 8.0, rtol=1e-6)

    def test_scaled_cube_sum(self):
        """Volume = DX*DY*DZ for non-unit cube."""
        m = _build_mass_matrix(2.0, 3.0, 5.0)
        assert np.isclose(m.sum(), 30.0, rtol=1e-6)

    def test_unit_cube_exact_entries(self):
        """Known exact entries for the unit cube mass matrix."""
        m = _build_mass_matrix(1.0, 1.0, 1.0)
        # Node 1 (0,0,0) to itself: (2*2*2)/216
        assert np.isclose(m[0, 0], 8.0 / 216.0, rtol=1e-6)
        # Node 1 (0,0,0) to node 2 (1,0,0): share y,z → (1*2*2)/216
        assert np.isclose(m[0, 1], 4.0 / 216.0, rtol=1e-6)
        # Node 1 (0,0,0) to node 4 (1,1,0): share z → (1*1*2)/216
        assert np.isclose(m[0, 3], 2.0 / 216.0, rtol=1e-6)
        # Node 1 (0,0,0) to node 8 (1,1,1): all differ → (1*1*1)/216
        assert np.isclose(m[0, 7], 1.0 / 216.0, rtol=1e-6)


class TestEvalFieldAtDir:
    """Tests for EVAL_FIELD_AT_DIR."""

    def test_stokes_i_isotropic(self):
        """Isotropic field (L=0 only): result = coeff * Y00."""
        ml, mm = 2, 2
        nlm = _nlm(ml, mm)
        nstokes, nstleg = 1, 1
        coeff_val = 3.5

        coeffs = np.zeros((nstokes, nlm), dtype=np.float32, order='F')
        coeffs[0, 0] = coeff_val

        mu, phi = 0.6, 1.2
        ylmdir = _ylmall(False, mu, phi, ml, mm, nstleg)

        result = _eval_field(coeffs, ylmdir, 1)
        expected = coeff_val * ylmdir[0, 0]
        assert np.isclose(result, expected, rtol=1e-6)

    def test_stokes_i_dipole(self):
        """L=1, M=0 mode only: result depends on mu."""
        ml, mm = 2, 2
        nlm = _nlm(ml, mm)
        nstokes, nstleg = 1, 1

        # L=1, M=0 is J=3 in SHDOM ordering (L=0:M=0→J=1; L=1:M=-1→J=2,M=0→J=3,M=+1→J=4)
        # 0-indexed: index 2
        coeffs = np.zeros((nstokes, nlm), dtype=np.float32, order='F')
        coeffs[0, 2] = 2.0

        mu, phi = 0.7, 0.5
        ylmdir = _ylmall(False, mu, phi, ml, mm, nstleg)

        result = _eval_field(coeffs, ylmdir, 1)
        expected = 2.0 * ylmdir[0, 2]
        assert np.isclose(result, expected, rtol=1e-6)

    def test_stokes_q_coupling(self):
        """Stokes Q couples NSTLEG components 2 and 5."""
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        nstokes, nstleg = 3, 6

        rng = np.random.default_rng(42)
        coeffs = rng.standard_normal((nstokes, nlm)).astype(np.float32)
        coeffs = np.asfortranarray(coeffs)

        mu, phi = 0.5, 0.3
        ylmdir = _ylmall(False, mu, phi, ml, mm, nstleg)

        result = _eval_field(coeffs, ylmdir, 2)

        # Python reference: Q = sum_j coeffs[1,j]*ylm[1,j] + coeffs[2,j]*ylm[4,j]
        # (0-indexed: NSTLEG components 2→index1, 5→index4)
        q_ref = np.sum(coeffs[1, :] * ylmdir[1, :] +
                       coeffs[2, :] * ylmdir[4, :])
        assert np.isclose(result, q_ref, rtol=1e-5)

    def test_stokes_u_coupling(self):
        """Stokes U couples NSTLEG components 6 and 3."""
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        nstokes, nstleg = 3, 6

        rng = np.random.default_rng(123)
        coeffs = rng.standard_normal((nstokes, nlm)).astype(np.float32)
        coeffs = np.asfortranarray(coeffs)

        mu, phi = -0.3, 2.1
        ylmdir = _ylmall(False, mu, phi, ml, mm, nstleg)

        result = _eval_field(coeffs, ylmdir, 3)

        # Python reference: U = sum_j coeffs[1,j]*ylm[5,j] + coeffs[2,j]*ylm[2,j]
        # (0-indexed: NSTLEG components 6→index5, 3→index2)
        u_ref = np.sum(coeffs[1, :] * ylmdir[5, :] +
                       coeffs[2, :] * ylmdir[2, :])
        assert np.isclose(result, u_ref, rtol=1e-5)

    def test_stokes_v(self):
        """Stokes V uses NSTLEG component 4."""
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        nstokes, nstleg = 4, 6

        rng = np.random.default_rng(7)
        coeffs = rng.standard_normal((nstokes, nlm)).astype(np.float32)
        coeffs = np.asfortranarray(coeffs)

        mu, phi = 0.8, 4.0
        ylmdir = _ylmall(False, mu, phi, ml, mm, nstleg)

        result = _eval_field(coeffs, ylmdir, 4)

        # Python reference: V = sum_j coeffs[3,j]*ylm[3,j]
        # (0-indexed: NSTLEG component 4→index3)
        v_ref = np.sum(coeffs[3, :] * ylmdir[3, :])
        assert np.isclose(result, v_ref, rtol=1e-5)


# ============================================================================
# Phase 2: Synthetic tests with small grids
# ============================================================================

class TestInnerProductSynthetic:
    """Synthetic tests for COMPUTE_GRAD_INNER_PRODUCT on small grids."""

    def _call_inner_product(self, nstokes, nstleg, ml, mm, nlm, npts, ncells,
                            gridptr, treeptr, gridpos, cellflags,
                            rshptr, shptr, radiance, source, dirflux,
                            rshptr_pf, shptr_pf, radiance_pf, source_pf,
                            dirflux_pf, solarmu, solaraz, ylmsun,
                            det_stokes, cellwise_det, mu_pf_cell, phi_pf_cell):
        """Wrapper that calls the Fortran subroutine and returns POINT_SENS."""
        point_sens, ierr, errmsg = at3d.core.compute_grad_inner_product(
            nx=1, ny=1, nz=1,
            ml=ml, mm=mm,
            gridptr=gridptr, treeptr=treeptr,
            gridpos=gridpos, cellflags=cellflags,
            rshptr=rshptr, shptr=shptr,
            radiance=radiance, source=source, dirflux=dirflux,
            rshptr_pf=rshptr_pf, shptr_pf=shptr_pf,
            radiance_pf=radiance_pf, source_pf=source_pf,
            dirflux_pf=dirflux_pf,
            solarmu=float(solarmu), solaraz=float(solaraz), ylmsun=ylmsun,
            det_stokes=int(det_stokes), cellwise_det=int(cellwise_det),
            mu_pf_cell=mu_pf_cell, phi_pf_cell=phi_pf_cell,
        )
        assert ierr == 0, f"Fortran error: {errmsg}"
        return point_sens

    def test_uniform_isotropic_no_solar(self):
        """Uniform L=0 field, no solar → analytic result.

        With only L=0 mode populated and DIRFLUX=0:
          Term A: P_AA(k) = a * b (all corners same)
          A_IC(i) = a*b * sum_k M(i,k) = a*b * vol/8
          Term B: P_JA(i,k) = c * b
          B_IC(i) = -c*b * vol/8
          POINT_SENS = (a*b - c*b) * vol/8 = b*(a-c)*vol/8
        """
        nstokes, nstleg = 1, 1
        ml, mm = 0, 0
        nlm = 1
        npts, ncells = 8, 1
        dx, dy, dz = 2.0, 3.0, 4.0

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)

        a, b, c = 1.5, 2.0, 0.5  # rad, rad_pf, source values

        radiance = np.full((nstokes, npts * nlm), a, dtype=np.float32, order='F')
        source = np.full((nstokes, npts * nlm), c, dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)

        radiance_pf = np.full((nstokes, npts * nlm), b, dtype=np.float32, order='F')
        source_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        dirflux_pf = np.zeros(npts, dtype=np.float32)

        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        solarmu = 0.5
        solaraz = 0.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        point_sens = self._call_inner_product(
            nstokes, nstleg, ml, mm, nlm, npts, ncells,
            gridptr, treeptr, gridpos, cellflags,
            rshptr, shptr, radiance, source, dirflux,
            rshptr_pf, shptr_pf, radiance_pf, source_pf,
            dirflux_pf, solarmu, solaraz, ylmsun,
            1, 0, mu_pf_cell, phi_pf_cell,
        )

        vol = dx * dy * dz
        expected = b * (a - c) * vol / 8.0
        assert np.allclose(point_sens, expected, rtol=1e-5), \
            f"Expected {expected}, got {point_sens}"

    def test_parseval_random_nstokes1(self):
        """Random SH coefficients, NSTOKES=1: Term A matches numpy Parseval."""
        nstokes, nstleg = 1, 1
        ml, mm = 4, 4
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.0, 1.0, 1.0

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(99)
        radiance = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance = np.asfortranarray(radiance)
        radiance_pf = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance_pf = np.asfortranarray(radiance_pf)

        # Zero out source and dirflux to isolate Term A
        source = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)

        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        solarmu = 0.5
        solaraz = 0.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        point_sens = self._call_inner_product(
            nstokes, nstleg, ml, mm, nlm, npts, ncells,
            gridptr, treeptr, gridpos, cellflags,
            rshptr, shptr, radiance, source, dirflux,
            rshptr_pf, shptr_pf, radiance_pf, source_pf,
            dirflux_pf, solarmu, solaraz, ylmsun,
            1, 0, mu_pf_cell, phi_pf_cell,
        )

        # Compute reference: P_AA(k) = sum_j sign(j) * rad(1,j) * rad_pf(1,j)
        signs = _sign_adj(ml, mm).astype(np.float32)
        m = _build_mass_matrix(dx, dy, dz)

        p_aa = np.zeros(8, dtype=np.float32)
        for k in range(8):
            rad_k = radiance[0, rshptr[k]:rshptr[k] + nlm]
            rad_pf_k = radiance_pf[0, rshptr_pf[k]:rshptr_pf[k] + nlm]
            p_aa[k] = np.sum(signs * rad_k * rad_pf_k)

        expected = m @ p_aa
        assert np.allclose(point_sens, expected, rtol=1e-5), \
            f"Max error: {np.max(np.abs(point_sens - expected))}"

    def test_term_c_only(self):
        """Isolate Term C: adjoint Stokes-I at solar direction times DIRFLUX."""
        nstokes, nstleg = 1, 1
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.5, 2.0, 2.5

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)

        # Kill A, B: set forward radiance=0 and source=0
        radiance = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        source = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')

        # Kill D, E: set DIRFLUX_PF = 0
        dirflux_pf = np.zeros(npts, dtype=np.float32)

        # Set DIRFLUX nonzero
        rng = np.random.default_rng(55)
        dirflux = np.abs(rng.standard_normal(npts).astype(np.float32))

        # Set RADIANCE_PF to random
        radiance_pf = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance_pf = np.asfortranarray(radiance_pf)
        source_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')

        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        solarmu = 0.7
        solaraz = 1.5
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf_cell = np.array([0.3], dtype=np.float32)
        phi_pf_cell = np.array([0.8], dtype=np.float32)

        point_sens = self._call_inner_product(
            nstokes, nstleg, ml, mm, nlm, npts, ncells,
            gridptr, treeptr, gridpos, cellflags,
            rshptr, shptr, radiance, source, dirflux,
            rshptr_pf, shptr_pf, radiance_pf, source_pf,
            dirflux_pf, solarmu, solaraz, ylmsun,
            1, 0, mu_pf_cell, phi_pf_cell,
        )

        # Reference: C_IC(i) = sum_k M(i,k) * IADJ_SOL(k) * DIRFLUX(k)
        # IADJ_SOL(k) = sum_j SIGN(j) * YLMSUN(0,j) * RAD_PF(0, rshptr_pf[k]+j)
        signs = _sign_adj(ml, mm).astype(np.float32)
        iadj_sol = np.zeros(8, dtype=np.float32)
        for k in range(8):
            rad_pf_k = radiance_pf[0, rshptr_pf[k]:rshptr_pf[k] + nlm]
            iadj_sol[k] = np.sum(signs * ylmsun[0, :] * rad_pf_k)

        m = _build_mass_matrix(dx, dy, dz)
        expected = m @ (iadj_sol * dirflux)
        assert np.allclose(point_sens, expected, rtol=1e-5), \
            f"Max error: {np.max(np.abs(point_sens - expected))}"

    def test_term_d_only(self):
        """Isolate Term D: forward radiance at detector direction times DIRFLUX_PF."""
        nstokes, nstleg = 1, 1
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.0, 1.0, 1.0

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(77)

        # Kill A, B: RADIANCE_PF = 0
        radiance_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        # Kill C: DIRFLUX = 0
        dirflux = np.zeros(npts, dtype=np.float32)
        # Kill E: SOURCE = 0
        source = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')

        # Set RADIANCE to random and DIRFLUX_PF nonzero
        radiance = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance = np.asfortranarray(radiance)
        dirflux_pf = np.abs(rng.standard_normal(npts).astype(np.float32))

        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        solarmu = 0.5
        solaraz = 0.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        # Detector direction: pseudo-forward at (mu=0.4, phi=1.0)
        # Actual detector: mu_det=-0.4, phi_det=1.0+pi
        mu_pf = np.float32(0.4)
        phi_pf = np.float32(1.0)
        mu_pf_cell = np.array([mu_pf], dtype=np.float32)
        phi_pf_cell = np.array([phi_pf], dtype=np.float32)

        point_sens = self._call_inner_product(
            nstokes, nstleg, ml, mm, nlm, npts, ncells,
            gridptr, treeptr, gridpos, cellflags,
            rshptr, shptr, radiance, source, dirflux,
            rshptr_pf, shptr_pf, radiance_pf, source_pf,
            dirflux_pf, solarmu, solaraz, ylmsun,
            1, 0, mu_pf_cell, phi_pf_cell,
        )

        # Reference: D_IC(i) = sum_k M(i,k) * DIRFLUX_PF(k) * IFWD_DET(k)
        # Detector direction: -mu_pf, phi_pf + pi
        pi_val = np.float32(np.pi)
        ylm_det = _ylmall(False, -float(mu_pf), float(phi_pf) + float(pi_val),
                          ml, mm, nstleg)

        ifwd_det = np.zeros(8, dtype=np.float32)
        for k in range(8):
            rad_k = radiance[0, rshptr[k]:rshptr[k] + nlm]
            ifwd_det[k] = np.sum(rad_k * ylm_det[0, :])

        m = _build_mass_matrix(dx, dy, dz)
        expected = m @ (dirflux_pf * ifwd_det)
        assert np.allclose(point_sens, expected, rtol=1e-5), \
            f"Max error: {np.max(np.abs(point_sens - expected))}"

    def test_sign_flip_l1(self):
        """L=1 mode should contribute with negative sign in Parseval."""
        nstokes, nstleg = 1, 1
        ml, mm = 1, 1
        nlm = _nlm(ml, mm)  # = 4 (L=0:1, L=1:-1,0,+1)
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(1.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)
        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        # Only populate L=1, M=0 (J=3 in 1-based, index 2 in 0-based)
        # Ordering: L=0:M=0→J=1; L=1:M=-1→J=2, M=0→J=3, M=+1→J=4
        radiance = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        radiance_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        for k in range(npts):
            radiance[0, rshptr[k] + 2] = 1.0     # L=1, M=0
            radiance_pf[0, rshptr_pf[k] + 2] = 1.0

        source = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)

        solarmu = 0.5
        solaraz = 0.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        point_sens = self._call_inner_product(
            nstokes, nstleg, ml, mm, nlm, npts, ncells,
            gridptr, treeptr, gridpos, cellflags,
            rshptr, shptr, radiance, source, dirflux,
            rshptr_pf, shptr_pf, radiance_pf, source_pf,
            dirflux_pf, solarmu, solaraz, ylmsun,
            1, 0, mu_pf_cell, phi_pf_cell,
        )

        # L=1 sign = -1, so Term A = (-1)*1*1 = -1 at each corner
        # P_AA(k) = -1.0
        # A_IC(i) = sum_k M(i,k)*(-1) = -vol/8 = -1/8
        expected = -1.0 / 8.0
        assert np.allclose(point_sens, expected, rtol=1e-5)

    def test_stokes_v_sign(self):
        """NSTOKES=4: Stokes V contributes with negative sign in Parseval."""
        nstokes, nstleg = 4, 6
        ml, mm = 0, 0
        nlm = 1
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(1.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)
        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        # Only Stokes V (s=4, index 3) populated at L=0
        radiance = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        radiance_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        for k in range(npts):
            radiance[3, rshptr[k]] = 1.0  # V component
            radiance_pf[3, rshptr_pf[k]] = 1.0

        source = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)

        solarmu = 0.5
        solaraz = 0.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        point_sens = self._call_inner_product(
            nstokes, nstleg, ml, mm, nlm, npts, ncells,
            gridptr, treeptr, gridpos, cellflags,
            rshptr, shptr, radiance, source, dirflux,
            rshptr_pf, shptr_pf, radiance_pf, source_pf,
            dirflux_pf, solarmu, solaraz, ylmsun,
            1, 0, mu_pf_cell, phi_pf_cell,
        )

        # L=0 → sign=+1, but s=4 → extra -1
        # P_AA(k) = -1 * 1 * (+1) * 1 = -1
        # A_IC(i) = -1.0 * vol/8 = -1/8
        expected = -1.0 / 8.0
        assert np.allclose(point_sens, expected, rtol=1e-5)

        # Verify that Stokes I gives +1/8 with same setup
        radiance[:] = 0.0
        radiance_pf[:] = 0.0
        for k in range(npts):
            radiance[0, rshptr[k]] = 1.0
            radiance_pf[0, rshptr_pf[k]] = 1.0

        point_sens_i = self._call_inner_product(
            nstokes, nstleg, ml, mm, nlm, npts, ncells,
            gridptr, treeptr, gridpos, cellflags,
            rshptr, shptr, radiance, source, dirflux,
            rshptr_pf, shptr_pf, radiance_pf, source_pf,
            dirflux_pf, solarmu, solaraz, ylmsun,
            1, 0, mu_pf_cell, phi_pf_cell,
        )
        expected_i = +1.0 / 8.0
        assert np.allclose(point_sens_i, expected_i, rtol=1e-5)

    def test_all_five_terms_combined(self):
        """All terms active: verify against full Python reference."""
        nstokes, nstleg = 1, 1
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.0, 2.0, 3.0

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)
        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        rng = np.random.default_rng(2026)
        radiance = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance = np.asfortranarray(radiance)
        source = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        source = np.asfortranarray(source)
        radiance_pf = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance_pf = np.asfortranarray(radiance_pf)
        source_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')

        dirflux = np.abs(rng.standard_normal(npts).astype(np.float32))
        dirflux_pf = np.abs(rng.standard_normal(npts).astype(np.float32))

        solarmu = 0.6
        solaraz = 2.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf = np.float32(0.3)
        phi_pf = np.float32(0.7)
        mu_pf_cell = np.array([mu_pf], dtype=np.float32)
        phi_pf_cell = np.array([phi_pf], dtype=np.float32)

        point_sens = self._call_inner_product(
            nstokes, nstleg, ml, mm, nlm, npts, ncells,
            gridptr, treeptr, gridpos, cellflags,
            rshptr, shptr, radiance, source, dirflux,
            rshptr_pf, shptr_pf, radiance_pf, source_pf,
            dirflux_pf, solarmu, solaraz, ylmsun,
            1, 0, mu_pf_cell, phi_pf_cell,
        )

        # Full Python reference
        signs = _sign_adj(ml, mm).astype(np.float32)
        m = _build_mass_matrix(dx, dy, dz)

        pi_val = float(np.pi)
        ylm_det = _ylmall(False, -float(mu_pf), float(phi_pf) + pi_val,
                          ml, mm, nstleg)

        expected = np.zeros(8, dtype=np.float32)
        for i in range(8):
            # Term A: sum_k M(i,k) * P_AA(k)
            p_aa = np.zeros(8, dtype=np.float32)
            for k in range(8):
                rad_k = radiance[0, rshptr[k]:rshptr[k] + nlm]
                rpf_k = radiance_pf[0, rshptr_pf[k]:rshptr_pf[k] + nlm]
                p_aa[k] = np.sum(signs * rad_k * rpf_k)
            a_i = np.sum(m[i, :] * p_aa)

            # Term B: - sum_k M(i,k) * P_JA(i,k)
            src_i = source[0, shptr[i]:shptr[i] + nlm]
            b_i = 0.0
            for k in range(8):
                rpf_k = radiance_pf[0, rshptr_pf[k]:rshptr_pf[k] + nlm]
                p_ja = np.sum(signs * src_i * rpf_k)
                b_i += m[i, k] * p_ja
            b_i = -b_i

            # Term C: sum_k M(i,k) * IADJ_SOL(k) * DIRFLUX(k)
            iadj_sol = np.zeros(8, dtype=np.float32)
            for k in range(8):
                rpf_k = radiance_pf[0, rshptr_pf[k]:rshptr_pf[k] + nlm]
                iadj_sol[k] = np.sum(signs * ylmsun[0, :] * rpf_k)
            c_i = np.sum(m[i, :] * iadj_sol * dirflux)

            # Term D: sum_k M(i,k) * DIRFLUX_PF(k) * IFWD_DET(k)
            ifwd_det = np.zeros(8, dtype=np.float32)
            for k in range(8):
                rad_k = radiance[0, rshptr[k]:rshptr[k] + nlm]
                ifwd_det[k] = np.sum(rad_k * ylm_det[0, :])
            d_i = np.sum(m[i, :] * dirflux_pf * ifwd_det)

            # Term E: -sum_k M(i,k) * DIRFLUX_PF(k) * JFWD_DET(k)
            jfwd_det = np.zeros(8, dtype=np.float32)
            for k in range(8):
                src_k = source[0, shptr[k]:shptr[k] + nlm]
                jfwd_det[k] = np.sum(src_k * ylm_det[0, :])
            e_i = -np.sum(m[i, :] * dirflux_pf * jfwd_det)

            expected[i] = a_i + b_i + c_i + d_i + e_i

        assert np.allclose(point_sens, expected, rtol=1e-4), \
            f"Max error: {np.max(np.abs(point_sens - expected))}"


# ============================================================================
# Phase 4: Edge cases
# ============================================================================

class TestEdgeCases:
    """Edge case and regression tests."""

    def test_zero_volume_cell(self):
        """Degenerate cell (DX=0) contributes zero, no NaN."""
        nstokes, nstleg = 1, 1
        ml, mm = 2, 2
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(0.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)
        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        rng = np.random.default_rng(111)
        radiance = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance = np.asfortranarray(radiance)
        radiance_pf = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance_pf = np.asfortranarray(radiance_pf)
        source = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)

        solarmu = 0.5
        solaraz = 0.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        point_sens, ierr, errmsg = at3d.core.compute_grad_inner_product(
            nx=1, ny=1, nz=1,
            ml=ml, mm=mm,
            gridptr=gridptr, treeptr=treeptr,
            gridpos=gridpos, cellflags=cellflags,
            rshptr=rshptr, shptr=shptr,
            radiance=radiance, source=source, dirflux=dirflux,
            rshptr_pf=rshptr_pf, shptr_pf=shptr_pf,
            radiance_pf=radiance_pf, source_pf=source_pf,
            dirflux_pf=dirflux_pf,
            solarmu=float(solarmu), solaraz=float(solaraz), ylmsun=ylmsun,
            det_stokes=1, cellwise_det=0,
            mu_pf_cell=mu_pf_cell, phi_pf_cell=phi_pf_cell,
        )

        assert ierr == 0
        assert np.all(np.isfinite(point_sens))
        assert np.allclose(point_sens, 0.0)

    def test_all_zero_fields(self):
        """All-zero inputs → all-zero output."""
        nstokes, nstleg = 1, 1
        ml, mm = 2, 2
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(1.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)
        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        z = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        zv = np.zeros(npts, dtype=np.float32)

        solarmu = 0.5
        solaraz = 0.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        point_sens, ierr, errmsg = at3d.core.compute_grad_inner_product(
            nx=1, ny=1, nz=1,
            ml=ml, mm=mm,
            gridptr=gridptr, treeptr=treeptr,
            gridpos=gridpos, cellflags=cellflags,
            rshptr=rshptr, shptr=shptr,
            radiance=z.copy(), source=z.copy(), dirflux=zv.copy(),
            rshptr_pf=rshptr_pf, shptr_pf=shptr_pf,
            radiance_pf=z.copy(), source_pf=z.copy(),
            dirflux_pf=zv.copy(),
            solarmu=float(solarmu), solaraz=float(solaraz), ylmsun=ylmsun,
            det_stokes=1, cellwise_det=0,
            mu_pf_cell=mu_pf_cell, phi_pf_cell=phi_pf_cell,
        )

        assert ierr == 0
        assert np.allclose(point_sens, 0.0)

    def test_non_leaf_cells_skipped(self):
        """Non-leaf cells (TREEPTR(2,IC)!=0) are skipped."""
        nstokes, nstleg = 1, 1
        ml, mm = 2, 2
        nlm = _nlm(ml, mm)
        npts = 8
        ncells = 2  # 2 cells but both point to same 8 nodes

        # Two cells, first is a parent (treeptr(2,1)=1), second is leaf
        gridptr = np.zeros((8, ncells), dtype=np.int32, order='F')
        for i in range(8):
            gridptr[i, 0] = i + 1
            gridptr[i, 1] = i + 1
        treeptr = np.zeros((2, ncells), dtype=np.int32, order='F')
        treeptr[1, 0] = 1  # Non-leaf (has children)
        treeptr[1, 1] = 0  # Leaf
        cellflags = np.zeros(ncells, dtype=np.int16)

        gridpos = np.zeros((3, npts), dtype=np.float32, order='F')
        for i in range(8):
            gridpos[0, i] = ((i) % 2) * 1.0
            gridpos[1, i] = ((i // 2) % 2) * 1.0
            gridpos[2, i] = (i // 4) * 1.0

        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)
        rshptr_pf = rshptr.copy()
        shptr_pf = shptr.copy()

        rng = np.random.default_rng(333)
        radiance = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance = np.asfortranarray(radiance)
        radiance_pf = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance_pf = np.asfortranarray(radiance_pf)
        source = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)

        solarmu = 0.5
        solaraz = 0.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf_cell = np.array([0.5, 0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0, 0.0], dtype=np.float32)

        # With 2 cells (1 parent, 1 leaf): result should be same as
        # if only the leaf cell were present
        point_sens_2cells, ierr, _ = at3d.core.compute_grad_inner_product(
            nx=1, ny=1, nz=1,
            ml=ml, mm=mm,
            gridptr=gridptr, treeptr=treeptr,
            gridpos=gridpos, cellflags=cellflags,
            rshptr=rshptr, shptr=shptr,
            radiance=radiance, source=source, dirflux=dirflux,
            rshptr_pf=rshptr_pf, shptr_pf=shptr_pf,
            radiance_pf=radiance_pf, source_pf=source_pf,
            dirflux_pf=dirflux_pf,
            solarmu=float(solarmu), solaraz=float(solaraz), ylmsun=ylmsun,
            det_stokes=1, cellwise_det=0,
            mu_pf_cell=mu_pf_cell, phi_pf_cell=phi_pf_cell,
        )
        assert ierr == 0

        # Single cell (only the leaf)
        point_sens_1cell, ierr, _ = at3d.core.compute_grad_inner_product(
            nx=1, ny=1, nz=1,
            ml=ml, mm=mm,
            gridptr=np.asfortranarray(gridptr[:, 1:2].copy()),
            treeptr=np.asfortranarray(treeptr[:, 1:2].copy()),
            gridpos=gridpos, cellflags=cellflags[:1].copy(),
            rshptr=rshptr, shptr=shptr,
            radiance=radiance, source=source, dirflux=dirflux,
            rshptr_pf=rshptr_pf, shptr_pf=shptr_pf,
            radiance_pf=radiance_pf, source_pf=source_pf,
            dirflux_pf=dirflux_pf,
            solarmu=float(solarmu), solaraz=float(solaraz), ylmsun=ylmsun,
            det_stokes=1, cellwise_det=0,
            mu_pf_cell=mu_pf_cell[:1].copy(),
            phi_pf_cell=phi_pf_cell[:1].copy(),
        )
        assert ierr == 0

        assert np.allclose(point_sens_2cells, point_sens_1cell, rtol=1e-6)

    def test_mismatched_sh_truncation(self):
        """Different SH truncation for forward and pseudo-forward: no crash."""
        nstokes, nstleg = 1, 1
        ml, mm = 4, 4
        nlm = _nlm(ml, mm)
        nlm_short = _nlm(2, 2)  # Shorter truncation for PF
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(1.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)
        shptr = _make_uniform_sh_ptrs(npts, nlm)
        rshptr_pf = _make_uniform_sh_ptrs(npts, nlm_short)
        shptr_pf = _make_uniform_sh_ptrs(npts, nlm_short)

        rng = np.random.default_rng(444)
        radiance = rng.standard_normal((nstokes, npts * nlm)).astype(np.float32)
        radiance = np.asfortranarray(radiance)
        source = np.zeros((nstokes, npts * nlm), dtype=np.float32, order='F')
        radiance_pf = rng.standard_normal((nstokes, npts * nlm_short)).astype(np.float32)
        radiance_pf = np.asfortranarray(radiance_pf)
        source_pf = np.zeros((nstokes, npts * nlm_short), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)

        solarmu = 0.5
        solaraz = 0.0
        ylmsun = _ylmall(True, solarmu, solaraz, ml, mm, nstleg)

        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        point_sens, ierr, errmsg = at3d.core.compute_grad_inner_product(
            nx=1, ny=1, nz=1,
            ml=ml, mm=mm,
            gridptr=gridptr, treeptr=treeptr,
            gridpos=gridpos, cellflags=cellflags,
            rshptr=rshptr, shptr=shptr,
            radiance=radiance, source=source, dirflux=dirflux,
            rshptr_pf=rshptr_pf, shptr_pf=shptr_pf,
            radiance_pf=radiance_pf, source_pf=source_pf,
            dirflux_pf=dirflux_pf,
            solarmu=float(solarmu), solaraz=float(solaraz), ylmsun=ylmsun,
            det_stokes=1, cellwise_det=0,
            mu_pf_cell=mu_pf_cell, phi_pf_cell=phi_pf_cell,
        )

        assert ierr == 0
        assert np.all(np.isfinite(point_sens))

        # Verify it uses MIN(NR_FWD, NR_PF) = nlm_short modes
        signs = _sign_adj(ml, mm).astype(np.float32)
        m = _build_mass_matrix(1.0, 1.0, 1.0)

        p_aa = np.zeros(8, dtype=np.float32)
        for k in range(8):
            rad_k = radiance[0, rshptr[k]:rshptr[k] + nlm_short]
            rpf_k = radiance_pf[0, rshptr_pf[k]:rshptr_pf[k] + nlm_short]
            p_aa[k] = np.sum(signs[:nlm_short] * rad_k * rpf_k)

        expected = m @ p_aa
        assert np.allclose(point_sens, expected, rtol=1e-5)


# ============================================================================
# NSTOKES=1 comprehensive verification
# ============================================================================

class TestScalarNSTOKES1:
    """Comprehensive tests for COMPUTE_GRAD_INNER_PRODUCT with NSTOKES=1.

    Each term is tested in isolation and in combination with a full Python
    reference implementation. All tests use NSTOKES=1, NSTLEG=1.
    """

    def _call(self, ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
              cellflags, rshptr, shptr, radiance, source, dirflux,
              rshptr_pf, shptr_pf, radiance_pf, source_pf, dirflux_pf,
              ylmsun, mu_pf_cell, phi_pf_cell):
        """Shortcut for calling the Fortran with NSTOKES=1, NSTLEG=1."""
        point_sens, ierr, errmsg = at3d.core.compute_grad_inner_product(
            nx=1, ny=1, nz=1, ml=ml, mm=mm,
            gridptr=gridptr, treeptr=treeptr,
            gridpos=gridpos, cellflags=cellflags,
            rshptr=rshptr, shptr=shptr,
            radiance=radiance, source=source, dirflux=dirflux,
            rshptr_pf=rshptr_pf, shptr_pf=shptr_pf,
            radiance_pf=radiance_pf, source_pf=source_pf,
            dirflux_pf=dirflux_pf,
            solarmu=0.5, solaraz=0.0, ylmsun=ylmsun,
            det_stokes=1, cellwise_det=0,
            mu_pf_cell=mu_pf_cell, phi_pf_cell=phi_pf_cell,
        )
        assert ierr == 0, f"Fortran error: {errmsg}"
        return point_sens

    def _python_reference(self, ml, mm, nlm, npts, rshptr, shptr, rshptr_pf,
                          radiance, source, dirflux, radiance_pf, dirflux_pf,
                          ylmsun, mu_pf, phi_pf, dx, dy, dz):
        """Full Python reference for NSTOKES=1 single-cell inner product."""
        signs = _sign_adj(ml, mm).astype(np.float32)
        m = _build_mass_matrix(dx, dy, dz)
        pi_val = float(np.pi)
        ylm_det = _ylmall(False, -mu_pf, phi_pf + pi_val, ml, mm, 1)

        expected = np.zeros(npts, dtype=np.float64)
        for i in range(npts):
            # Term A
            p_aa = np.zeros(npts, dtype=np.float64)
            for k in range(npts):
                nr_fwd = rshptr[k + 1] - rshptr[k]
                nr_pf = rshptr_pf[k + 1] - rshptr_pf[k]
                nmin = min(nr_fwd, nr_pf)
                rad_k = radiance[0, rshptr[k]:rshptr[k] + nmin]
                rpf_k = radiance_pf[0, rshptr_pf[k]:rshptr_pf[k] + nmin]
                p_aa[k] = np.sum(signs[:nmin] * rad_k * rpf_k)
            a_i = np.sum(m[i, :] * p_aa)

            # Term B
            ns_fwd = shptr[i + 1] - shptr[i]
            src_i = source[0, shptr[i]:shptr[i] + ns_fwd]
            b_i = 0.0
            for k in range(npts):
                nr_pf = rshptr_pf[k + 1] - rshptr_pf[k]
                nmin = min(ns_fwd, nr_pf)
                rpf_k = radiance_pf[0, rshptr_pf[k]:rshptr_pf[k] + nmin]
                p_ja = np.sum(signs[:nmin] * src_i[:nmin] * rpf_k)
                b_i += m[i, k] * p_ja
            b_i = -b_i

            # Term C
            c_i = 0.0
            for k in range(npts):
                nr_pf = rshptr_pf[k + 1] - rshptr_pf[k]
                rpf_k = radiance_pf[0, rshptr_pf[k]:rshptr_pf[k] + nr_pf]
                iadj = np.sum(signs[:nr_pf] * ylmsun[0, :nr_pf] * rpf_k)
                c_i += m[i, k] * iadj * dirflux[k]

            # Term D
            d_i = 0.0
            for k in range(npts):
                nr_fwd = rshptr[k + 1] - rshptr[k]
                rad_k = radiance[0, rshptr[k]:rshptr[k] + nr_fwd]
                ifwd = np.sum(rad_k * ylm_det[0, :nr_fwd])
                d_i += m[i, k] * dirflux_pf[k] * ifwd

            # Term E
            e_i = 0.0
            for k in range(npts):
                ns_k = shptr[k + 1] - shptr[k]
                src_k = source[0, shptr[k]:shptr[k] + ns_k]
                jfwd = np.sum(src_k * ylm_det[0, :ns_k])
                e_i += m[i, k] * dirflux_pf[k] * jfwd
            e_i = -e_i

            expected[i] = a_i + b_i + c_i + d_i + e_i

        return expected.astype(np.float32)

    # ------------------------------------------------------------------
    # Individual term isolation
    # ------------------------------------------------------------------

    def test_term_a_only_uniform(self):
        """Term A alone with uniform radiance fields: analytic check."""
        ml, mm = 2, 2
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(1.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        # Uniform L=0 only: rad=2.0, rad_pf=3.0
        radiance = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        radiance_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        for k in range(npts):
            radiance[0, rshptr[k]] = 2.0
            radiance_pf[0, rshptr[k]] = 3.0

        source = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                        cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                        rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                        dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        # P_AA(k) = 2*3*sign(L=0)=6 for all k. A_IC(i)=6*vol/8 = 6/8
        assert np.allclose(ps, 6.0 / 8.0, rtol=1e-5)

    def test_term_a_only_random(self):
        """Term A with random multi-mode radiance vs Python reference."""
        ml, mm = 5, 5
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 2.0, 1.5, 0.8

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(101)
        radiance = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        radiance_pf = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))

        source = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                        cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                        rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                        dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        expected = self._python_reference(
            ml, mm, nlm, npts, rshptr, rshptr.copy(), rshptr.copy(),
            radiance, source, dirflux, radiance_pf, dirflux_pf,
            ylmsun, 0.5, 0.0, dx, dy, dz)
        assert np.allclose(ps, expected, rtol=1e-5)

    def test_term_b_only_random(self):
        """Term B alone: source cross-product with adjoint radiance."""
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.0, 1.0, 1.0

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(202)
        radiance = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        source = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        radiance_pf = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))

        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                        cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                        rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                        dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        # Reference: B_IC(i) = -sum_k M(i,k) * sum_j sign(j)*src(i,j)*rpf(k,j)
        signs = _sign_adj(ml, mm).astype(np.float32)
        m = _build_mass_matrix(dx, dy, dz)
        expected = np.zeros(npts, dtype=np.float32)
        for i in range(npts):
            src_i = source[0, rshptr[i]:rshptr[i] + nlm]
            b_i = 0.0
            for k in range(npts):
                rpf_k = radiance_pf[0, rshptr[k]:rshptr[k] + nlm]
                b_i += m[i, k] * np.sum(signs * src_i * rpf_k)
            expected[i] = -b_i

        assert np.allclose(ps, expected, rtol=1e-5)

    def test_term_b_uniform_analytic(self):
        """Term B with uniform fields: B = -source*rad_pf * vol/8."""
        ml, mm = 0, 0
        nlm = 1
        npts, ncells = 8, 1
        dx, dy, dz = 3.0, 2.0, 1.0

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        s_val, b_val = 4.0, 2.5
        radiance = np.zeros((1, npts), dtype=np.float32, order='F')
        source = np.full((1, npts), s_val, dtype=np.float32, order='F')
        radiance_pf = np.full((1, npts), b_val, dtype=np.float32, order='F')
        source_pf = np.zeros((1, npts), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                        cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                        rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                        dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        expected = -s_val * b_val * (dx * dy * dz) / 8.0
        assert np.allclose(ps, expected, rtol=1e-5)

    def test_term_e_only_random(self):
        """Term E alone: forward source evaluated at detector direction."""
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.0, 1.0, 1.0

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(303)
        # Kill A,B: radiance_pf=0; Kill C: dirflux=0; Kill D: radiance=0
        radiance = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        radiance_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)

        source = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        dirflux_pf = np.abs(rng.standard_normal(npts).astype(np.float32))
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')

        mu_pf, phi_pf = 0.6, 1.3
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([mu_pf], dtype=np.float32)
        phi_pf_cell = np.array([phi_pf], dtype=np.float32)

        ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                        cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                        rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                        dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        # Reference: E_IC(i) = -sum_k M(i,k) * dirflux_pf(k) * JFWD_DET(k)
        pi_val = float(np.pi)
        ylm_det = _ylmall(False, -mu_pf, phi_pf + pi_val, ml, mm, 1)
        m = _build_mass_matrix(dx, dy, dz)

        jfwd_det = np.zeros(npts, dtype=np.float32)
        for k in range(npts):
            src_k = source[0, rshptr[k]:rshptr[k] + nlm]
            jfwd_det[k] = np.sum(src_k * ylm_det[0, :nlm])

        expected = -(m @ (dirflux_pf * jfwd_det))
        assert np.allclose(ps, expected, rtol=1e-5)

    # ------------------------------------------------------------------
    # Combination tests
    # ------------------------------------------------------------------

    def test_terms_a_plus_b(self):
        """Terms A+B combined with Python reference."""
        ml, mm = 4, 4
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.0, 1.0, 1.0

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(404)
        radiance = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        source = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        radiance_pf = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))

        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                        cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                        rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                        dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        expected = self._python_reference(
            ml, mm, nlm, npts, rshptr, rshptr.copy(), rshptr.copy(),
            radiance, source, dirflux, radiance_pf, dirflux_pf,
            ylmsun, 0.5, 0.0, dx, dy, dz)
        assert np.allclose(ps, expected, rtol=1e-5)

    def test_terms_c_d_e_combined(self):
        """Terms C+D+E combined (direct-beam terms only)."""
        ml, mm = 4, 4
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.5, 2.5, 0.5

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(505)
        radiance_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        radiance = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        source = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        dirflux = np.abs(rng.standard_normal(npts).astype(np.float32))
        dirflux_pf = np.abs(rng.standard_normal(npts).astype(np.float32))

        mu_pf, phi_pf = 0.4, 2.1
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([mu_pf], dtype=np.float32)
        phi_pf_cell = np.array([phi_pf], dtype=np.float32)

        ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                        cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                        rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                        dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        expected = self._python_reference(
            ml, mm, nlm, npts, rshptr, rshptr.copy(), rshptr.copy(),
            radiance, source, dirflux, radiance_pf, dirflux_pf,
            ylmsun, mu_pf, phi_pf, dx, dy, dz)
        assert np.allclose(ps, expected, rtol=1e-4)

    def test_all_terms_high_order(self):
        """All five terms with ML=8: full reference check."""
        ml, mm = 8, 8
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.0, 1.0, 1.0

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(606)
        radiance = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        source = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        radiance_pf = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.abs(rng.standard_normal(npts).astype(np.float32))
        dirflux_pf = np.abs(rng.standard_normal(npts).astype(np.float32))

        mu_pf, phi_pf = 0.35, 0.9
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([mu_pf], dtype=np.float32)
        phi_pf_cell = np.array([phi_pf], dtype=np.float32)

        ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                        cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                        rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                        dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        expected = self._python_reference(
            ml, mm, nlm, npts, rshptr, rshptr.copy(), rshptr.copy(),
            radiance, source, dirflux, radiance_pf, dirflux_pf,
            ylmsun, mu_pf, phi_pf, dx, dy, dz)
        assert np.allclose(ps, expected, rtol=1e-4)

    # ------------------------------------------------------------------
    # Property / structural tests
    # ------------------------------------------------------------------

    def test_sign_alternation_m0_only(self):
        """Verify (-1)^L alternates correctly for L=0..4 with M=0 only."""
        ml, mm = 4, 0  # 5 modes (L=0,1,2,3,4 each with M=0)
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(1.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        # All coefficients = 1 → P_AA(k) = sum_L (-1)^L = +1-1+1-1+1 = 1
        radiance = np.ones((1, npts * nlm), dtype=np.float32, order='F')
        radiance_pf = np.ones((1, npts * nlm), dtype=np.float32, order='F')
        source = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                        cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                        rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                        dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        assert np.allclose(ps, 1.0 / 8.0, rtol=1e-5)

    def test_linearity_in_radiance(self):
        """POINT_SENS (Term A) scales linearly with forward radiance."""
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(1.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(707)
        radiance = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        radiance_pf = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        source = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.zeros(npts, dtype=np.float32)
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        ps1 = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                         cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                         rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                         dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        factor = np.float32(2.5)
        ps2 = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                         cellflags, rshptr, rshptr.copy(),
                         np.asfortranarray(radiance * factor), source, dirflux,
                         rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                         dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        assert np.allclose(ps2, ps1 * factor, rtol=1e-5)

    def test_linearity_in_radiance_pf(self):
        """Terms A,B,C scale linearly with pseudo-forward radiance."""
        ml, mm = 3, 3
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(1.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(808)
        radiance = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        source = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        radiance_pf = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.abs(rng.standard_normal(npts).astype(np.float32))
        dirflux_pf = np.zeros(npts, dtype=np.float32)  # Kill D,E
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([0.4], dtype=np.float32)
        phi_pf_cell = np.array([1.0], dtype=np.float32)

        ps1 = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                         cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                         rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                         dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        factor = np.float32(-3.0)
        ps2 = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                         cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                         rshptr.copy(), rshptr.copy(),
                         np.asfortranarray(radiance_pf * factor), source_pf,
                         dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        assert np.allclose(ps2, ps1 * factor, rtol=1e-5)

    def test_detector_direction_dependence(self):
        """Different detector directions give different Term D results."""
        ml, mm = 4, 4
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(1.0, 1.0, 1.0)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(909)
        radiance = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        radiance_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        source = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.zeros(npts, dtype=np.float32)
        dirflux_pf = np.ones(npts, dtype=np.float32)
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)

        ps1 = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                         cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                         rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                         dirflux_pf, ylmsun,
                         np.array([0.3], dtype=np.float32),
                         np.array([0.5], dtype=np.float32))

        ps2 = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                         cellflags, rshptr, rshptr.copy(), radiance, source, dirflux,
                         rshptr.copy(), rshptr.copy(), radiance_pf, source_pf,
                         dirflux_pf, ylmsun,
                         np.array([0.8], dtype=np.float32),
                         np.array([2.5], dtype=np.float32))

        assert not np.allclose(ps1, ps2, rtol=1e-3)

        # Both match Python reference
        exp1 = self._python_reference(ml, mm, nlm, npts, rshptr, rshptr.copy(),
                                      rshptr.copy(), radiance, source, dirflux,
                                      radiance_pf, dirflux_pf, ylmsun,
                                      0.3, 0.5, 1.0, 1.0, 1.0)
        exp2 = self._python_reference(ml, mm, nlm, npts, rshptr, rshptr.copy(),
                                      rshptr.copy(), radiance, source, dirflux,
                                      radiance_pf, dirflux_pf, ylmsun,
                                      0.8, 2.5, 1.0, 1.0, 1.0)
        assert np.allclose(ps1, exp1, rtol=1e-5)
        assert np.allclose(ps2, exp2, rtol=1e-5)

    def test_volume_scaling(self):
        """POINT_SENS scales with cell volume for uniform fields."""
        ml, mm = 0, 0
        nlm = 1
        npts, ncells = 8, 1
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([0.5], dtype=np.float32)
        phi_pf_cell = np.array([0.0], dtype=np.float32)

        results = []
        for (dx, dy, dz) in [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (2.0, 3.0, 4.0)]:
            gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
            rshptr = _make_uniform_sh_ptrs(npts, nlm)
            radiance = np.full((1, npts), 1.0, dtype=np.float32, order='F')
            radiance_pf = np.full((1, npts), 1.0, dtype=np.float32, order='F')
            source = np.zeros((1, npts), dtype=np.float32, order='F')
            source_pf = np.zeros((1, npts), dtype=np.float32, order='F')
            dirflux = np.zeros(npts, dtype=np.float32)
            dirflux_pf = np.zeros(npts, dtype=np.float32)

            ps = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                            cellflags, rshptr, rshptr.copy(), radiance, source,
                            dirflux, rshptr.copy(), rshptr.copy(), radiance_pf,
                            source_pf, dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)
            results.append(ps[0])

        vols = [1.0, 2.0, 24.0]
        for r, v in zip(results, vols):
            assert np.isclose(r, v / 8.0, rtol=1e-5)

    def test_superposition(self):
        """Sum of isolated terms equals the combined result."""
        ml, mm = 4, 4
        nlm = _nlm(ml, mm)
        npts, ncells = 8, 1
        dx, dy, dz = 1.0, 2.0, 1.5

        gridptr, treeptr, gridpos, cellflags = _make_single_cell_grid(dx, dy, dz)
        rshptr = _make_uniform_sh_ptrs(npts, nlm)

        rng = np.random.default_rng(1010)
        radiance = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        source = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        radiance_pf = np.asfortranarray(
            rng.standard_normal((1, npts * nlm)).astype(np.float32))
        source_pf = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        dirflux = np.abs(rng.standard_normal(npts).astype(np.float32))
        dirflux_pf = np.abs(rng.standard_normal(npts).astype(np.float32))

        mu_pf, phi_pf = 0.4, 1.2
        ylmsun = _ylmall(True, 0.5, 0.0, ml, mm, 1)
        mu_pf_cell = np.array([mu_pf], dtype=np.float32)
        phi_pf_cell = np.array([phi_pf], dtype=np.float32)

        z = np.zeros((1, npts * nlm), dtype=np.float32, order='F')
        zv = np.zeros(npts, dtype=np.float32)

        # Full
        ps_full = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                             cellflags, rshptr, rshptr.copy(), radiance, source,
                             dirflux, rshptr.copy(), rshptr.copy(), radiance_pf,
                             source_pf, dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        # A only: source=0, dirflux=0, dirflux_pf=0
        ps_a = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                          cellflags, rshptr, rshptr.copy(), radiance, z.copy(),
                          zv.copy(), rshptr.copy(), rshptr.copy(), radiance_pf,
                          source_pf, zv.copy(), ylmsun, mu_pf_cell, phi_pf_cell)

        # B only: radiance=0, dirflux=0, dirflux_pf=0
        ps_b = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                          cellflags, rshptr, rshptr.copy(), z.copy(), source,
                          zv.copy(), rshptr.copy(), rshptr.copy(), radiance_pf,
                          source_pf, zv.copy(), ylmsun, mu_pf_cell, phi_pf_cell)

        # C only: radiance=0, source=0, dirflux_pf=0
        ps_c = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                          cellflags, rshptr, rshptr.copy(), z.copy(), z.copy(),
                          dirflux, rshptr.copy(), rshptr.copy(), radiance_pf,
                          source_pf, zv.copy(), ylmsun, mu_pf_cell, phi_pf_cell)

        # D+E only: radiance_pf=0, dirflux=0
        ps_de = self._call(ml, mm, nlm, npts, ncells, gridptr, treeptr, gridpos,
                           cellflags, rshptr, rshptr.copy(), radiance, source,
                           zv.copy(), rshptr.copy(), rshptr.copy(), z.copy(),
                           source_pf, dirflux_pf, ylmsun, mu_pf_cell, phi_pf_cell)

        assert np.allclose(ps_full, ps_a + ps_b + ps_c + ps_de, rtol=1e-5)
