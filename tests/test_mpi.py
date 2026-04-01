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


class TestAssignSolversToRanks:
    """Unit tests for :func:`at3d.parallel.assign_solvers_to_ranks`."""

    def test_more_ranks_than_solvers(self):
        """Extra ranks should be idle (empty list)."""
        from at3d.parallel import assign_solvers_to_ranks

        class FakeComm:
            def __init__(self, rank, size):
                self._rank = rank
                self._size = size
            def Get_rank(self):
                return self._rank
            def Get_size(self):
                return self._size

        assert assign_solvers_to_ranks(2, FakeComm(0, 4)) == [0]
        assert assign_solvers_to_ranks(2, FakeComm(1, 4)) == [1]
        assert assign_solvers_to_ranks(2, FakeComm(2, 4)) == []
        assert assign_solvers_to_ranks(2, FakeComm(3, 4)) == []

    def test_more_solvers_than_ranks(self):
        """Solvers should be distributed round-robin."""
        from at3d.parallel import assign_solvers_to_ranks

        class FakeComm:
            def __init__(self, rank, size):
                self._rank = rank
                self._size = size
            def Get_rank(self):
                return self._rank
            def Get_size(self):
                return self._size

        assert assign_solvers_to_ranks(5, FakeComm(0, 2)) == [0, 2, 4]
        assert assign_solvers_to_ranks(5, FakeComm(1, 2)) == [1, 3]

    def test_equal_ranks_and_solvers(self):
        """One solver per rank."""
        from at3d.parallel import assign_solvers_to_ranks

        class FakeComm:
            def __init__(self, rank, size):
                self._rank = rank
                self._size = size
            def Get_rank(self):
                return self._rank
            def Get_size(self):
                return self._size

        for r in range(3):
            assert assign_solvers_to_ranks(3, FakeComm(r, 3)) == [r]

    def test_single_rank(self):
        """Single rank should get all solvers."""
        from at3d.parallel import assign_solvers_to_ranks

        class FakeComm:
            def Get_rank(self):
                return 0
            def Get_size(self):
                return 1

        assert assign_solvers_to_ranks(4, FakeComm()) == [0, 1, 2, 3]

    def test_zero_solvers(self):
        """No solvers means no work for any rank."""
        from at3d.parallel import assign_solvers_to_ranks

        class FakeComm:
            def Get_rank(self):
                return 0
            def Get_size(self):
                return 4

        assert assign_solvers_to_ranks(0, FakeComm()) == []


class TestSplitCommForSolvers:
    """Test :func:`at3d.parallel.split_comm_for_solvers` with real MPI."""

    def test_split_returns_valid_comm(self):
        from at3d.parallel import split_comm_for_solvers
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        n_solvers = max(size, 1)
        solver_comm, my_indices = split_comm_for_solvers(n_solvers, comm)
        assert len(my_indices) >= 1
        assert solver_comm != MPI.COMM_NULL
        # Sub-communicator should have size 1 (one rank per solver).
        assert solver_comm.Get_size() == 1
        solver_comm.Free()

    def test_idle_rank_gets_comm_null(self):
        """When more ranks than solvers, extra ranks get COMM_NULL."""
        from at3d.parallel import split_comm_for_solvers
        comm = MPI.COMM_WORLD
        if comm.Get_size() < 2:
            pytest.skip("Need at least 2 ranks")
        # Only 1 solver — all ranks except rank 0 should be idle.
        solver_comm, my_indices = split_comm_for_solvers(1, comm)
        if comm.Get_rank() == 0:
            assert len(my_indices) == 1
            assert solver_comm != MPI.COMM_NULL
            solver_comm.Free()
        else:
            assert len(my_indices) == 0
            assert solver_comm == MPI.COMM_NULL

    def test_start_mpi_with_sub_comm(self):
        """Fortran start_mpi should accept the sub-communicator handle."""
        import at3d.core
        from at3d.parallel import split_comm_for_solvers
        comm = MPI.COMM_WORLD
        solver_comm, my_indices = split_comm_for_solvers(
            comm.Get_size(), comm)
        if my_indices:
            masterproc, _ = at3d.core.start_mpi(solver_comm.py2f())
            # Each sub-comm has 1 rank, so all are masterproc.
            assert masterproc is True or masterproc == 1
            solver_comm.Free()
