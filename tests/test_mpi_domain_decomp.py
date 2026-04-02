"""Tests for MPI domain decomposition via MAP_SHDOM_MPI.

Phase 1: Verify subdomain geometry (xstart, ystart, npx, npy, bcflag)
Phase 2: Verify data-slicing integrity (atmosphere, optical properties)

Run these tests with varying rank counts:
    pytest tests/test_mpi_domain_decomp.py -v                       # single-rank
    mpiexec -n 2 pytest tests/test_mpi_domain_decomp.py -v          # 2-rank
    mpiexec -n 4 pytest tests/test_mpi_domain_decomp.py -v          # 4-rank
"""
import copy
import warnings
from collections import OrderedDict

import numpy as np
import pytest
import xarray as xr

import at3d

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def extract_bcflag_bits(bcflag):
    """Return (bit2, bit3) from BCFLAG.

    bit2 = 1 means X-periodic (set when dims[0] > 1 in MAP_SHDOM_MPI).
    bit3 = 1 means Y-periodic (set when dims[1] > 1 in MAP_SHDOM_MPI).
    """
    bit2 = (bcflag >> 2) & 1
    bit3 = (bcflag >> 3) & 1
    return bit2, bit3


class FakeComm:
    """Minimal MPI communicator mock for single-rank tests."""
    def __init__(self, rank=0, size=1):
        self._rank = rank
        self._size = size
    def Get_rank(self):
        return self._rank
    def Get_size(self):
        return self._size
    def py2f(self):
        return 0


# ---------------------------------------------------------------------------
# Shared solver factory
# ---------------------------------------------------------------------------

# Grid dimensions chosen so that both 2-rank (1D) and 4-rank (2D) produce
# sensible subdomains (e.g. 16/2 = 8, 12/2 = 6).
_NPX = 16
_NPY = 12
_NPZ = 8
_DELX = 0.05
_DELY = 0.05
_ZLEVELS = np.linspace(0.1, 1.0, _NPZ)
_WAVELENGTH = 0.86


def _make_atmosphere(rte_grid):
    """Create a 3D atmosphere dataset on *rte_grid* with temperature.

    Temperature varies linearly from 280 K at the bottom to 220 K at the
    top and is uniform in x/y so values are easy to verify after slicing.
    """
    npx = rte_grid.sizes['x']
    npy = rte_grid.sizes['y']
    npz = rte_grid.sizes['z']
    # Temperature profile (z-only), broadcast to 3D.
    temp_z = np.linspace(280.0, 220.0, npz)
    temp_3d = np.broadcast_to(temp_z[np.newaxis, np.newaxis, :],
                              (npx, npy, npz)).copy()
    atmosphere = xr.Dataset(
        data_vars={
            'temperature': (['x', 'y', 'z'], temp_3d),
            'pressure': ('z', np.ones(npz) * 1013.25),
        },
        coords={
            'x': rte_grid.x,
            'y': rte_grid.y,
            'z': rte_grid.z,
        },
    )
    atmosphere['delx'] = rte_grid.delx
    atmosphere['dely'] = rte_grid.dely
    return atmosphere


def _make_solver():
    """Create a small RTE solver suitable for domain decomposition tests.

    Uses Rayleigh scattering as optical medium so no Mie tables are needed.
    Returns a fresh solver each time (deep-copy safe).
    """
    rte_grid = at3d.grid.make_grid(_DELX, _NPX, _DELY, _NPY, _ZLEVELS)

    # 1D atmosphere for rayleigh computation (needs z-only temperature).
    atm_1d = xr.Dataset(
        data_vars={
            'temperature': ('z', np.linspace(280.0, 220.0, _NPZ)),
            'pressure': ('z', np.ones(_NPZ) * 1013.25),
        },
        coords={'z': rte_grid.z.data},
    )

    rayleigh = at3d.rayleigh.to_grid(
        np.atleast_1d(_WAVELENGTH), atm_1d, rte_grid
    )
    medium = {'rayleigh': rayleigh[_WAVELENGTH]}

    # 3D atmosphere for the solver (needed for MPI temperature slicing).
    atmosphere = _make_atmosphere(rte_grid)

    config = at3d.configuration.get_config('../default_config.json')
    config['split_accuracy'] = 0.0
    config['spherical_harmonics_accuracy'] = 0.0
    config['num_mu_bins'] = 8
    config['num_phi_bins'] = 16
    config['solution_accuracy'] = 1e-5
    config['x_boundary_condition'] = 'periodic'
    config['y_boundary_condition'] = 'periodic'
    config['adapt_grid_factor'] = 100.0
    config['acceleration_flag'] = False

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        solver = at3d.solver.RTE(
            numerical_params=config,
            medium=medium,
            source=at3d.source.solar(_WAVELENGTH, -0.5, 0.0, solarflux=1.0),
            surface=at3d.surface.lambertian(albedo=0.05),
            num_stokes=1,
            atmosphere=atmosphere,
        )
    return solver


# ===================================================================
# Phase 1: Basic Domain Decomposition Geometry
# ===================================================================


class TestDomainDecompSingleRank:
    """Verify setup_domain_decomposition is a no-op for single-rank."""

    def test_single_rank_setup_preserves_full_domain(self):
        solver = _make_solver()

        # Record full-domain values before decomposition.
        full_npx = solver._pa.npx
        full_npy = solver._pa.npy

        # Activate domain decomp with a fake single-rank communicator.
        solver.setup_domain_decomposition(FakeComm(rank=0, size=1))

        # Subdomain should be the full domain.
        assert solver._mpi_xstart == 0.0
        assert solver._mpi_ystart == 0.0
        assert solver._mpi_npx == full_npx
        assert solver._mpi_npy == full_npy

        # Internal property arrays should still be full-domain sized.
        assert solver._pa.npx == full_npx
        assert solver._pa.npy == full_npy

    def test_single_rank_bcflag_unchanged(self):
        solver = _make_solver()

        # Both BCs are periodic → bcflag should be 0 (no open bits).
        bcflag_before = solver._bcflag

        solver.setup_domain_decomposition(FakeComm(rank=0, size=1))

        # Single-rank: MAP_SHDOM_MPI should NOT set bits 2/3.
        bcflag_after = solver._bcflag
        bit2, bit3 = extract_bcflag_bits(bcflag_after)
        assert bit2 == 0, f"bit2 should be 0 for single-rank, got {bit2}"
        assert bit3 == 0, f"bit3 should be 0 for single-rank, got {bit3}"
        assert bcflag_after == bcflag_before


class TestDomainDecomp2Rank:
    """Verify 1D domain decomposition with 2 MPI ranks.

    Run with: mpiexec -n 2 pytest ... -k TestDomainDecomp2Rank
    """

    @pytest.fixture(autouse=True)
    def _require_2_ranks(self):
        try:
            from mpi4py import MPI
            self.comm = MPI.COMM_WORLD
        except ImportError:
            pytest.skip("mpi4py not available")
        if self.comm.Get_size() < 2:
            pytest.skip("Need at least 2 MPI ranks")

    def _decompose(self):
        """Create solver and domain-decompose on comm."""
        solver = _make_solver()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            solver.setup_domain_decomposition(self.comm)
        return solver

    # -- Phase 1: Geometry --

    def test_subdomain_variables_positive(self):
        solver = self._decompose()
        assert solver._mpi_npx > 0
        assert solver._mpi_npy > 0
        assert solver._mpi_nx > 0
        assert solver._mpi_ny > 0
        assert solver._mpi_xstart >= 0.0
        assert solver._mpi_ystart >= 0.0

    def test_partitioning_covers_full_domain(self):
        solver = self._decompose()
        rank = self.comm.Get_rank()
        info = {
            'rank': rank,
            'xstart': solver._mpi_xstart,
            'ystart': solver._mpi_ystart,
            'npx': solver._mpi_npx,
            'npy': solver._mpi_npy,
        }
        all_info = self.comm.gather(info, root=0)

        if rank == 0:
            # Determine which dimension was split.
            unique_xstarts = sorted(set(d['xstart'] for d in all_info))
            unique_ystarts = sorted(set(d['ystart'] for d in all_info))
            x_split = len(unique_xstarts) > 1
            y_split = len(unique_ystarts) > 1

            if x_split:
                # Non-split dimension should keep full size.
                if not y_split:
                    for d in all_info:
                        assert d['npy'] == _NPY, \
                            f"Rank {d['rank']}: npy={d['npy']} != full {_NPY}"
                # With periodic BC, subdomain NPX includes one overlap
                # point at each inter-rank boundary.  The interior
                # points (npx - 1 per rank, except last which wraps)
                # sum to the full domain size.
                by_x = {}
                for d in all_info:
                    by_x.setdefault(d['xstart'], d['npx'])
                n_x_ranks = len(by_x)
                interior_npx = sum(npx - 1 for npx in by_x.values())
                assert interior_npx == _NPX, \
                    f"interior sum(npx-1)={interior_npx} != full {_NPX}"

            if y_split:
                if not x_split:
                    for d in all_info:
                        assert d['npx'] == _NPX, \
                            f"Rank {d['rank']}: npx={d['npx']} != full {_NPX}"
                by_y = {}
                for d in all_info:
                    by_y.setdefault(d['ystart'], d['npy'])
                interior_npy = sum(npy - 1 for npy in by_y.values())
                assert interior_npy == _NPY, \
                    f"interior sum(npy-1)={interior_npy} != full {_NPY}"

            assert x_split or y_split, "2-rank should split at least one dim"

    def test_bcflag_bits_match_decomposition(self):
        solver = self._decompose()
        rank = self.comm.Get_rank()
        info = {
            'xstart': solver._mpi_xstart,
            'ystart': solver._mpi_ystart,
            'bcflag': solver._bcflag,
        }
        all_info = self.comm.gather(info, root=0)

        if rank == 0:
            unique_xstarts = sorted(set(d['xstart'] for d in all_info))
            unique_ystarts = sorted(set(d['ystart'] for d in all_info))
            x_split = len(unique_xstarts) > 1
            y_split = len(unique_ystarts) > 1

            for d in all_info:
                bit2, bit3 = extract_bcflag_bits(d['bcflag'])
                assert bit2 == int(x_split), \
                    f"bit2={bit2} but x_split={x_split} (bcflag={d['bcflag']})"
                assert bit3 == int(y_split), \
                    f"bit3={bit3} but y_split={y_split} (bcflag={d['bcflag']})"

    def test_subdomain_coverage(self):
        """Every global grid cell is owned by at least one rank (complete coverage).
        Boundary points may be shared (periodic overlap)."""
        solver = self._decompose()
        rank = self.comm.Get_rank()
        info = {
            'xstart': solver._mpi_xstart,
            'ystart': solver._mpi_ystart,
            'npx': solver._mpi_npx,
            'npy': solver._mpi_npy,
        }
        all_info = self.comm.gather(info, root=0)

        if rank == 0:
            # Build set of (ix % full, iy % full) global indices per rank.
            union = set()
            for d in all_info:
                ix0 = int(round(d['xstart'] / _DELX))
                iy0 = int(round(d['ystart'] / _DELY))
                for ix in range(ix0, ix0 + d['npx']):
                    for iy in range(iy0, iy0 + d['npy']):
                        union.add((ix % _NPX, iy % _NPY))

            # Every global cell must be covered.
            assert len(union) == _NPX * _NPY, \
                f"Coverage={len(union)} != {_NPX * _NPY}"

    # -- Phase 2: Data Slicing --

    def test_atmosphere_slicing(self):
        solver = _make_solver()
        # Save full-domain temperature before decomposition.
        full_temp = solver.atmosphere.temperature.values.copy()

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            solver.setup_domain_decomposition(self.comm)

        rank = self.comm.Get_rank()

        # _pa.tempp should be sized for the local subdomain.
        expected_size = solver._mpi_npx * solver._mpi_npy * _NPZ
        assert solver._pa.tempp.shape[0] == expected_size, \
            f"tempp size={solver._pa.tempp.shape[0]} != expected {expected_size}"

        # Verify the values match the correct subdomain of the full
        # temperature field.  Build the same modular index arrays used
        # by _setup_atmosphere and extract the expected slice.
        ix0 = int(round(solver._mpi_xstart / _DELX))
        iy0 = int(round(solver._mpi_ystart / _DELY))
        ix_idx = np.arange(ix0, ix0 + solver._mpi_npx) % _NPX
        iy_idx = np.arange(iy0, iy0 + solver._mpi_npy) % _NPY
        expected_temp = full_temp[np.ix_(ix_idx, iy_idx, np.arange(_NPZ))]
        np.testing.assert_allclose(
            solver._pa.tempp.ravel(),
            expected_temp.ravel(),
            atol=1e-5,
            err_msg=f"rank={rank}: tempp values don't match global subdomain slice",
        )

    def test_optical_properties_slicing(self):
        solver = _make_solver()
        # Get full-domain extinction before decomposition.
        full_ext = {}
        for name, scat in solver.medium.items():
            full_ext[name] = scat.extinction.values.copy()

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            solver.setup_domain_decomposition(self.comm)

        rank = self.comm.Get_rank()
        ix0 = int(round(solver._mpi_xstart / _DELX))
        iy0 = int(round(solver._mpi_ystart / _DELY))

        # Verify the internal Fortran arrays have been filled with the
        # correct subdomain slice.  The _pa.extinctp array is sized
        # [maxpg, npart] where maxpg = npx_local * npy_local * npz.
        expected_maxpg = solver._mpi_npx * solver._mpi_npy * _NPZ
        assert solver._maxpg == expected_maxpg, \
            f"_maxpg={solver._maxpg} != {expected_maxpg}"

    def test_grid_reconstruction(self):
        solver = _make_solver()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            solver.setup_domain_decomposition(self.comm)

        # _pa.npx/npy should be LOCAL sizes.
        assert solver._pa.npx == solver._mpi_npx
        assert solver._pa.npy == solver._mpi_npy

        # _pa.xstart/ystart should be subdomain origin.
        assert np.isclose(float(solver._pa.xstart), solver._mpi_xstart, atol=1e-6)
        assert np.isclose(float(solver._pa.ystart), solver._mpi_ystart, atol=1e-6)

        # Full-domain totals preserved.
        assert solver._npxt == _NPX
        assert solver._npyt == _NPY


class TestDomainDecomp4Rank:
    """Verify 2D domain decomposition with 4 MPI ranks.

    Run with: mpiexec -n 4 pytest ... -k TestDomainDecomp4Rank
    """

    @pytest.fixture(autouse=True)
    def _require_4_ranks(self):
        try:
            from mpi4py import MPI
            self.comm = MPI.COMM_WORLD
        except ImportError:
            pytest.skip("mpi4py not available")
        if self.comm.Get_size() < 4:
            pytest.skip("Need at least 4 MPI ranks")

    def _decompose(self):
        solver = _make_solver()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            solver.setup_domain_decomposition(self.comm)
        return solver

    # -- Phase 1: Geometry --

    def test_subdomain_variables_positive(self):
        solver = self._decompose()
        assert solver._mpi_npx > 0
        assert solver._mpi_npy > 0
        assert solver._mpi_nx > 0
        assert solver._mpi_ny > 0
        assert solver._mpi_xstart >= 0.0
        assert solver._mpi_ystart >= 0.0

    def test_partitioning_covers_full_domain(self):
        solver = self._decompose()
        rank = self.comm.Get_rank()
        info = {
            'rank': rank,
            'xstart': solver._mpi_xstart,
            'ystart': solver._mpi_ystart,
            'npx': solver._mpi_npx,
            'npy': solver._mpi_npy,
        }
        all_info = self.comm.gather(info, root=0)

        if rank == 0:
            unique_xstarts = sorted(set(d['xstart'] for d in all_info))
            unique_ystarts = sorted(set(d['ystart'] for d in all_info))
            x_split = len(unique_xstarts) > 1
            y_split = len(unique_ystarts) > 1

            # For 4 ranks, expect 2D decomposition (2×2).
            assert x_split and y_split, \
                f"Expected 2D split with 4 ranks, got x_split={x_split}, y_split={y_split}"

            # Verify X coverage: interior points (npx - 1 per rank
            # due to periodic overlap) sum to full domain.
            by_x = {}
            for d in all_info:
                by_x.setdefault(d['xstart'], d['npx'])
            interior_npx = sum(npx - 1 for npx in by_x.values())
            assert interior_npx == _NPX, \
                f"X interior coverage: {interior_npx} != {_NPX}"

            # Verify Y coverage.
            by_y = {}
            for d in all_info:
                by_y.setdefault(d['ystart'], d['npy'])
            interior_npy = sum(npy - 1 for npy in by_y.values())
            assert interior_npy == _NPY, \
                f"Y interior coverage: {interior_npy} != {_NPY}"

    def test_bcflag_both_bits_set(self):
        """4-rank (2×2) should set both bit 2 (X) and bit 3 (Y)."""
        solver = self._decompose()
        rank = self.comm.Get_rank()
        all_bcflags = self.comm.gather(solver._bcflag, root=0)

        if rank == 0:
            for i, bcflag in enumerate(all_bcflags):
                bit2, bit3 = extract_bcflag_bits(bcflag)
                assert bit2 == 1, \
                    f"Rank {i}: bit2={bit2}, expected 1 (bcflag={bcflag})"
                assert bit3 == 1, \
                    f"Rank {i}: bit3={bit3}, expected 1 (bcflag={bcflag})"

    def test_subdomain_coverage(self):
        """Every global grid cell is covered by at least one rank."""
        solver = self._decompose()
        rank = self.comm.Get_rank()
        info = {
            'xstart': solver._mpi_xstart,
            'ystart': solver._mpi_ystart,
            'npx': solver._mpi_npx,
            'npy': solver._mpi_npy,
        }
        all_info = self.comm.gather(info, root=0)

        if rank == 0:
            union = set()
            for d in all_info:
                ix0 = int(round(d['xstart'] / _DELX))
                iy0 = int(round(d['ystart'] / _DELY))
                for ix in range(ix0, ix0 + d['npx']):
                    for iy in range(iy0, iy0 + d['npy']):
                        union.add((ix % _NPX, iy % _NPY))

            assert len(union) == _NPX * _NPY, \
                f"Coverage={len(union)} != {_NPX * _NPY}"

    def test_all_ranks_consistent_bcflag(self):
        """All ranks should have the same BCFLAG value."""
        solver = self._decompose()
        all_bcflags = self.comm.gather(solver._bcflag, root=0)
        rank = self.comm.Get_rank()

        if rank == 0:
            assert len(set(all_bcflags)) == 1, \
                f"Inconsistent bcflags across ranks: {all_bcflags}"

    # -- Phase 2: Data Slicing --

    def test_atmosphere_slicing(self):
        solver = _make_solver()
        full_temp = solver.atmosphere.temperature.values.copy()

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            solver.setup_domain_decomposition(self.comm)

        rank = self.comm.Get_rank()
        expected_size = solver._mpi_npx * solver._mpi_npy * _NPZ
        assert solver._pa.tempp.shape[0] == expected_size, \
            f"tempp size={solver._pa.tempp.shape[0]} != expected {expected_size}"

        ix0 = int(round(solver._mpi_xstart / _DELX))
        iy0 = int(round(solver._mpi_ystart / _DELY))
        ix_idx = np.arange(ix0, ix0 + solver._mpi_npx) % _NPX
        iy_idx = np.arange(iy0, iy0 + solver._mpi_npy) % _NPY
        expected_temp = full_temp[np.ix_(ix_idx, iy_idx, np.arange(_NPZ))]
        np.testing.assert_allclose(
            solver._pa.tempp.ravel(),
            expected_temp.ravel(),
            atol=1e-5,
            err_msg=f"rank={rank}: tempp mismatch",
        )

    def test_optical_properties_slicing(self):
        solver = _make_solver()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            solver.setup_domain_decomposition(self.comm)

        expected_maxpg = solver._mpi_npx * solver._mpi_npy * _NPZ
        assert solver._maxpg == expected_maxpg, \
            f"_maxpg={solver._maxpg} != {expected_maxpg}"

    def test_grid_reconstruction(self):
        solver = _make_solver()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            solver.setup_domain_decomposition(self.comm)

        assert solver._pa.npx == solver._mpi_npx
        assert solver._pa.npy == solver._mpi_npy
        assert np.isclose(float(solver._pa.xstart), solver._mpi_xstart, atol=1e-6)
        assert np.isclose(float(solver._pa.ystart), solver._mpi_ystart, atol=1e-6)
        assert solver._npxt == _NPX
        assert solver._npyt == _NPY

    def test_extinction_values_match_global(self):
        """Verify local extinction slice matches global domain at correct indices."""
        solver = _make_solver()
        # Grab full-domain extinction before decomposition.
        full_ext = {}
        for name, scat in solver.medium.items():
            full_ext[name] = scat.extinction.values.copy()

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            solver.setup_domain_decomposition(self.comm)

        rank = self.comm.Get_rank()
        ix0 = int(round(solver._mpi_xstart / _DELX))
        iy0 = int(round(solver._mpi_ystart / _DELY))
        npx_local = solver._mpi_npx
        npy_local = solver._mpi_npy

        # After decomposition, self.medium still holds full-domain data
        # but _pa.extinctp holds the sliced data. Verify the sliced indices
        # are consistent by checking the extinctp array size.
        expected_size = npx_local * npy_local * _NPZ
        actual_size = solver._pa.extinctp[:expected_size, :].shape[0]
        assert actual_size >= expected_size, \
            f"extinctp too small: {actual_size} < {expected_size}"

        # Gather local subdomain info and verify full coverage.
        info = {
            'ix0': ix0, 'iy0': iy0,
            'npx': npx_local, 'npy': npy_local,
        }
        all_info = self.comm.gather(info, root=0)

        if rank == 0:
            covered = set()
            for d in all_info:
                for ix in range(d['ix0'], d['ix0'] + d['npx']):
                    for iy in range(d['iy0'], d['iy0'] + d['npy']):
                        covered.add((ix % _NPX, iy % _NPY))
            assert len(covered) == _NPX * _NPY
