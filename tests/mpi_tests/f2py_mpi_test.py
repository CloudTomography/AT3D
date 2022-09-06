"""
compile

f2py --mpiexec=mpifort -c helloworld.f90 -m helloworld

and then

mpirun -n 2 python f2py_mpi_test.py

as a basic test of mpi installation.
"""

from mpi4py import MPI

import helloworld

comm = MPI.COMM_WORLD

fcomm = comm.py2f()


print('I am ', comm.Get_rank(), ' of ', comm.Get_size())


helloworld.sayhello(fcomm)
