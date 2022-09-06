! file: helloworld.f90
subroutine sayhello(comm)
  use mpi
  implicit none
  integer :: comm, rank, size, ierr, message_Item
!f2py intent(in,out) :: comm
  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)
  print *, 'Hello, World! I am process ',rank,' of ',size,'.'

  IF (rank == 0) THEN
      message_Item = 42
      call MPI_SEND(message_Item, 1, MPI_INT, 1, 1, comm, ierr)
      print *, "Sending message containing: ", message_Item
  ELSE IF (rank == 1) THEN
    call MPI_RECV(message_Item, 1, MPI_INT, 0, 1, comm, MPI_STATUS_IGNORE, ierr)
      print *, "Received message containing: ", message_Item
  END IF

end subroutine sayhello
