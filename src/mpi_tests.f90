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

! SUBROUTINE START_MPI (MASTERPROC, COMM)
!  ! Initializes the MPI system, and gets the number of processors in use
!  ! and the current processor.  If this is processor 0 then MASTERPROC
!  ! is returned true.  Also starts the logging files "shdom???.log"
!  ! for the standard output from the non-master processors.
!  !USE SHDOM_MPI_DATA
!   USE MPI
!   IMPLICIT NONE
!   INTEGER :: COMM
! !f2py intent(in,out) :: comm
!   LOGICAL, INTENT(OUT) :: MASTERPROC
! !f2py intent(out) :: MASTERPROC
!   INTEGER :: ierr
!
!   ! call MPI_INIT(ierr)
!   ! IF (ierr .ne. MPI_SUCCESS) THEN
!   !   WRITE (6,*) 'Error starting MPI version of SHDOM. Terminating.'
!   !   STOP
!   ! ENDIF
!
!   call MPI_COMM_SIZE (COMM, numproc, ierr)
!   call MPI_COMM_RANK (COMM, myproc, ierr)
!
!   call MPI_Errhandler_set (COMM, MPI_ERRORS_RETURN, ierr)
! !  call MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL, ierr)
!
!   StartTime = MPI_WTIME()
!
!   IF (myproc == 0) THEN
!     MASTERPROC = .TRUE.
!   ELSE
!     MASTERPROC = .FALSE.
!   ENDIF
! END SUBROUTINE START_MPI
