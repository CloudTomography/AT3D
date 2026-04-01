"""Tests verifying that MPI dependencies are installed and that the
Fortran MPI routines can be called through f2py.

These tests are skipped automatically when ``mpi4py`` is not installed.
They work both as a single process (``pytest tests/test_mpi.py``) and
under ``mpiexec`` (``mpiexec -n 2 python -m pytest tests/test_mpi.py``).
"""
import pytest

pytest.importorskip("mpi4py", reason="mpi4py not installed")
from mpi4py import MPI


class TestMPI4PyImport:
    """Verify that mpi4py itself is functional."""

    def test_comm_world_exists(self):
        comm = MPI.COMM_WORLD
        assert comm is not None

    def test_comm_world_size(self):
        size = MPI.COMM_WORLD.Get_size()
        assert size >= 1

    def test_comm_world_rank(self):
        rank = MPI.COMM_WORLD.Get_rank()
        assert 0 <= rank < MPI.COMM_WORLD.Get_size()

    def test_py2f_returns_int(self):
        handle = MPI.COMM_WORLD.py2f()
        assert isinstance(handle, int)


class TestCoreStartMPI:
    """Verify that ``at3d.core.start_mpi`` can be called via f2py."""

    def test_start_mpi_returns_masterproc(self):
        import at3d.core
        comm_handle = MPI.COMM_WORLD.py2f()
        # start_mpi returns (masterproc, comm) because COMM is intent(in,out)
        masterproc, _ = at3d.core.start_mpi(comm_handle)
        rank = MPI.COMM_WORLD.Get_rank()
        if rank == 0:
            assert masterproc is True or masterproc == 1
        else:
            assert masterproc is False or masterproc == 0

    def test_start_mpi_idempotent(self):
        """Calling start_mpi multiple times should not crash."""
        import at3d.core
        comm_handle = MPI.COMM_WORLD.py2f()
        m1, _ = at3d.core.start_mpi(comm_handle)
        m2, _ = at3d.core.start_mpi(comm_handle)
        assert m1 == m2
