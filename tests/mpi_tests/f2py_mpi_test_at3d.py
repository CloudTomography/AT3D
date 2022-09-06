from mpi4py import MPI

comm = MPI.COMM_WORLD
fcomm = comm.py2f()

import at3d
#masterproc, fcomm = at3d.core.start_mpi(fcomm)

print('I am ', comm.Get_rank(), ' of ', comm.Get_size())

at3d.core.sayhello(fcomm)
