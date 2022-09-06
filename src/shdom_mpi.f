CMODULE SHDOM_MPI_DATA
C  INTEGER :: numproc, myproc, comm2d, iprocneigh(4), IXYPRP(2,2), IXYOUT(2,2)
C  LOGICAL, PARAMETER :: LoadBalance=.FALSE.
C  DOUBLE PRECISION :: StartTime
C  INTEGER, ALLOCATABLE :: BXYPTR(:,:,:), IBXYPTR(:,:,:)
C  REAL,    ALLOCATABLE :: BXYPOS(:,:,:)
CEND MODULE

SUBROUTINE START_MPI (MASTERPROC, COMM)
 ! Initializes the MPI system, and gets the number of processors in use
 ! and the current processor.  If this is processor 0 then MASTERPROC
 ! is returned true.  Also starts the logging files "shdom???.log"
 ! for the standard output from the non-master processors.
  USE SHDOM_MPI_DATA
  USE MPI
  IMPLICIT NONE
  INTEGER :: COMM
!f2py intent(in,out) :: comm
  LOGICAL, INTENT(OUT) :: MASTERPROC
!f2py intent(out) :: MASTERPROC
  INTEGER :: ierr

  ! call MPI_INIT(ierr)
  ! IF (ierr .ne. MPI_SUCCESS) THEN
  !   WRITE (6,*) 'Error starting MPI version of SHDOM. Terminating.'
  !   STOP
  ! ENDIF

  call MPI_COMM_SIZE (COMM, numproc, ierr)
  call MPI_COMM_RANK (COMM, myproc, ierr)

  call MPI_Errhandler_set (COMM, MPI_ERRORS_RETURN, ierr)
!  call MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL, ierr)

  StartTime = MPI_WTIME()

  IF (myproc == 0) THEN
    MASTERPROC = .TRUE.
  ELSE
    MASTERPROC = .FALSE.
  ENDIF
END SUBROUTINE START_MPI



CSUBROUTINE MAP_SHDOM_MPI (BCFLAG, NPX, NPY, NX, NY, NZ, DELX, DELY, PROPFILE, &
C                          XSTART, YSTART, RUNNAME, COMM)
! Sets up the mapping of processors for the SHDOM domain (creating the
! 2D communicator, comm2d).  Determines the subdomain property grid
! (NPX,NPY) and internal base grid (NX,NY) sizes, and outputs the
! starting position of this subdomain (XSTART,YSTART).  Bits 2 and 3
! of BCFLAG are set to indicate multiple processors in the X and Y
! directions, respectively.  Open boundary conditions (bits 0 and 1 in
! BCFLAG) are overridden by the multiple processor flags (i.e. multiple
! processor requires a periodic full domain).
!   The neighboring processors are stored in iprocneigh(:), with the
! following indexing scheme: 1=-X, 2=+X, 3-=-Y, 4=+Y.   The IXYPRP array
! (module global variable) is set with the starting and ending property
! grid indices of this subdomain (IXYPRP(:,1) is for X, and IXYPRP(:,2)
! is for Y).  The IXYOUT array is set with the starting index and number
! of grid points of the internal base grid for this subdomain (for the
! output arrays).
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(INOUT) :: BCFLAG, NPX, NPY, NX, NY, NZ, COMM
f2py intent(in, out) :: BCFLAG, NPX, NPY, NX, NY, NZ, COMM
C  REAL,    INTENT(IN) :: DELX, DELY
f2py intent(in) :: DELX, DELY
C  CHARACTER(LEN=80), INTENT(IN) :: PROPFILE, RUNNAME
f2py intent(in) :: PROPFILE, RUNNAME
C  REAL,    INTENT(OUT) :: XSTART, YSTART
f2py intent(out) :: XSTART, YSTART
C  INTEGER :: NPXT, NPYT, NXT, NYT
C  INTEGER :: i, ierr, dims(2), coords(2)
C  LOGICAL :: periodic(2) = (/ .TRUE., .TRUE. /), reorder=.FALSE.
C  REAL    :: XEND, YEND
C  CHARACTER(LEN=3) :: procstr
C  CHARACTER(LEN=72) :: logfile
C
C  IF (numproc == 1) THEN
C    XSTART = 0.0 ; YSTART = 0.0
C    RETURN
C  ENDIF
C
C  IF (myproc > 0) THEN
C    WRITE (procstr,'(I3.3)') myproc
C    logfile = TRIM(RUNNAME)//procstr//'.log'
C    OPEN (UNIT=6,FILE=logfile,STATUS='UNKNOWN')
C  ENDIF
C
C  NPXT=NPX ; NPYT=NPY ; NXT=NX ; NYT=NY
C
 ! Partition 2D Cartesian topology:
C  dims(:) = 0
C  call MPI_DIMS_CREATE (numproc, 2, dims, ierr)
C
  ! Create 2D Cartesian mapping of processors:
C  call MPI_CART_CREATE (COMM, 2, dims, periodic, reorder, &
C                        comm2d, ierr)
C
C  IF (dims(1) > 1) BCFLAG=IBSET(BCFLAG,2)
C  IF (dims(2) > 1) BCFLAG=IBSET(BCFLAG,3)
  ! Multi-processor BC takes precedence over open BCs
C  IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3))  BCFLAG=IAND(BCFLAG,12)
C
  ! Find where in the 2D grid of processors this is is at:
C  call MPI_CART_COORDS (comm2d, myproc, 2, coords, ierr)
C
   ! Make the optimal (or default) property grid divisions (in IXYPRP)
C  CALL OPTIMIZE_PROCESSOR_DOMAINS (dims, coords, NPXT, NPYT, PROPFILE, RUNNAME)
C
  ! Make the new property and base grid sizes and the IXYOUT arrays
C  IF (dims(1) > 1) THEN
C    NPX = IXYPRP(2,1)-IXYPRP(1,1)+1
C    NX = NINT(NXT*FLOAT(NPX-1)/NPXT)+1
C    IF (ABS(NPXT/FLOAT(NXT) - (NPX-1)/FLOAT(NX-1)) > 0.0001) THEN
C      CALL ABORT_SHDOM_MPI ('MAP_SHDOM_MPI: Even X base grid cannot be preserved for this NX and number of processors')
C    ENDIF
C    IXYOUT(1,1) = NINT(NXT*FLOAT(IXYPRP(1,1)-1)/NPXT)+1
C    IXYOUT(2,1) = NINT(NXT*FLOAT(NPX-1)/NPXT)
C    XSTART = (IXYPRP(1,1)-1)*DELX
C    XEND = (IXYPRP(2,1)-1)*DELX
C  ELSE
C    IXYOUT(1,1) = 1
C    IXYOUT(2,1) = NXT
C    XSTART = 0.0
C    XEND = NPXT*DELX
C  ENDIF
C
C  IF (dims(2) > 1) THEN
C    NPY = IXYPRP(2,2)-IXYPRP(1,2)+1
C    NY = NINT(NYT*FLOAT(NPY-1)/NPYT)+1
C    IF (ABS(NPYT/FLOAT(NYT) - (NPY-1)/FLOAT(NY-1)) > 0.0001) THEN
C      CALL ABORT_SHDOM_MPI ('MAP_SHDOM_MPI: Even Y base grid cannot be preserved for this NY and number of processors')
C    ENDIF
C    IXYOUT(1,2) = NINT(NYT*FLOAT(IXYPRP(1,2)-1)/NPYT)+1
C    IXYOUT(2,2) = NINT(NYT*FLOAT(NPY-1)/NPYT)
C    YSTART = (IXYPRP(1,2)-1)*DELY
C    YEND = (IXYPRP(2,2)-1)*DELY
C  ELSE
C    IXYOUT(1,2) = 1
C    IXYOUT(2,2) = NYT
C    YSTART = 0.0
C    YEND = NPYT*DELY
C  ENDIF
C
  ! Get the processors (rank numbers) of the neighbors
C  call MPI_CART_SHIFT (comm2d, 0, 1, iprocneigh(1), iprocneigh(2), ierr)
C  call MPI_CART_SHIFT (comm2d, 1, 1, iprocneigh(3), iprocneigh(4), ierr)
C
C  WRITE (6,'(A,2(1X,I2),A,I3)') 'Number of processors in X and Y:', dims(1:2),&
C                             '   This processor: ',myproc
C  WRITE (6,'(A,4(1X,I3))') 'Neighbors:',iprocneigh(1:4)
C  WRITE (6,'(2(A,2F7.3))')  'Processor X range: ',XSTART,XEND, '   Y range: ',YSTART,YEND
C  WRITE (6,*)
CEND SUBROUTINE MAP_SHDOM_MPI
C
CSUBROUTINE OPTIMIZE_PROCESSOR_DOMAINS (ProcDims, ProcCoord, NPXT, NPYT, &
C                                       PROPFILE, RUNNAME)
! Calculates the location of this processor (subdomain) in the property
! grid.  The starting and ending property grid point indices for X and Y
! are returned in IXYPRP.  Each processor sends its location in the
! 2D array of subdomains to the master process and the master process
! reads the RUNNAME//"_load_balance.inp" file (if there is one) and does
! the optimization.  If the "shdom_proc_times.inp" file of relative
! processor times exists then the optimization uses this information on
! heterogeneous processors to try to give more work to faster processors.
!
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(IN) :: ProcDims(2), ProcCoord(2), NPXT, NPYT
C  CHARACTER(LEN=80), INTENT(IN) :: PROPFILE, RUNNAME
C  INTEGER :: NPXTI, NPYTI
C  INTEGER :: IX, IY, JX, JY, Nxd, Nyd, N, I, ierr
C  INTEGER :: BndLines(ProcDims(1)+ProcDims(2)-2)
C  INTEGER :: IDX(ProcDims(1)+1), IDY(ProcDims(2)+1)
C  LOGICAL :: DoLoadBalancing
C  INTEGER, PARAMETER :: NCYCLE=30    ! number of annealing cooling cycles
C  INTEGER, PARAMETER :: NiterFac=100 ! number of iterations per length of optimization vector (BndLines)
C  REAL,    PARAMETER :: Temp0=0.1    ! initial temperature (units of change in objective function "energy")
C  REAL,    PARAMETER :: Tempfac=4.0/NCYCLE   ! exponential cooling rate
C  INTEGER :: CYCLE, ITER, NITER
C  INTEGER :: TrialLines(ProcDims(1)+ProcDims(2)-2), CurLines(ProcDims(1)+ProcDims(2)-2)
C  INTEGER :: iseed=-3
C  REAL    :: Temp, DelLines, u, ran
C  REAL    :: LOAD_OBJ_FUNC, Energy, NewEnergy, LowestEnergy, DeltaE, Boltz
C  INTEGER, ALLOCATABLE :: ProcCoordAll(:,:), Ngrid(:,:)
C  REAL,    ALLOCATABLE :: wtproc(:)
C  CHARACTER(LEN=80) :: PROPFILEI
C
C
  ! Initialize the subdomain boundary lines to evenly spaced
C  Nxd=ProcDims(1) ; Nyd=ProcDims(2) ; N=Nxd+Nyd-2
C  BndLines(1:Nxd-1) = (/ (NINT(I*FLOAT(NPXT)/Nxd)+1, I=1,Nxd-1) /)
C  BndLines(Nxd:N)   = (/ (NINT(I*FLOAT(NPYT)/Nyd)+1, I=1,Nyd-1) /)
C
C
C  IF (LoadBalance) THEN
C    DoLoadBalancing = .false.
C    IF (myproc == 0) THEN
      ! See if we can read the load balancing file and header matches
C      OPEN (UNIT=10, FILE=TRIM(RUNNAME)//'_load_balance.inp', STATUS='OLD', IOSTAT=ierr)
C      IF (ierr == 0) THEN
C        READ (10,*)
C        READ (10,'(A)') PROPFILEI
C        READ (10,*) NPXTI, NPYTI
C        IF (NPXTI /= NPXT .OR. NPYTI /= NPYT .OR. PROPFILEI /= PROPFILE) THEN
C          WRITE (6,'(A,A,A)') 'Load balancing file ',TRIM(RUNNAME),&
C               '_load_balance.inp is not compatible with property file:'
C          WRITE (6,*) PROPFILE
C          WRITE (6,*)
C        ELSE
C          DoLoadBalancing = .true.
C        ENDIF
C        CLOSE (10)
C      ENDIF
C    ENDIF
C
C
    ! Have the master process collect the locations of all
C    ALLOCATE (ProcCoordAll(2,0:numproc-1))
C    CALL MPI_GATHER (ProcCoord,2,MPI_INTEGER, ProcCoordAll,2,MPI_INTEGER, 0, comm2d, ierr)
C
C    IF (myproc==0 .AND. DoLoadBalancing) THEN
C       ! Read the load balancing file to get array of number of grid points
C      ALLOCATE (Ngrid(NPXT+1,NPYT+1))
C      OPEN (UNIT=10, FILE=TRIM(RUNNAME)//'_load_balance.inp', STATUS='OLD')
C      READ (10,*)
C      READ (10,*)
C      READ (10,*)
C      DO IX = 1, NPXT
C        DO IY = 1, NPYT
C          READ (10,*) JX, JY, Ngrid(IX,IY)
C        ENDDO
C      ENDDO
C      CLOSE (10)
C      Ngrid(NPXT+1,:) = Ngrid(1,:)
C      Ngrid(:,NPYT+1) = Ngrid(:,1)
C      Ngrid(NPXT+1,NPYT+1) = Ngrid(1,1)
C
C       ! Read the relative processor times if the "shdom_proc_times.inp" file exists
C      ALLOCATE (wtproc(numproc))
C      OPEN (UNIT=11, FILE='shdom_proc_times.inp', STATUS='OLD', IOSTAT=ierr)
C      IF (ierr == 0) THEN
C        DO i = 1, numproc
C          READ (11,*) wtproc(i)
C        ENDDO
C        CLOSE (11)
C        wtproc(:) = wtproc(:)*numproc/SUM(wtproc(:))
C      ELSE
C        wtproc(:) = 1.0
C      ENDIF
C
C       ! Do the simulated annealing to minimize the objective function
C      u = ran(iseed)
C      NITER = NiterFac*N
C      Energy = LOAD_OBJ_FUNC (BndLines, Nxd, Nyd, ProcCoordAll, wtproc, NPXT, NPYT, Ngrid)
C      LowestEnergy = Energy
C      CurLines(:) = BndLines(:)
C      ! print *, 'Annealing cycles: Cycle Temp Energy  BndLines'
C      ! print '(I3,2(1X,F9.6),10(1X,I3))', 0,Temp0,Energy, BndLines(:)
C
C       ! Loop through the cooling cycles
C      DO CYCLE = 1, NCYCLE
C        Temp = Temp0 * EXP(-Tempfac*(CYCLE-1))      ! lower the temp
C         ! Adjust the maximum distance a boundary line is perturbed with the cycle number
C        u = FLOAT(CYCLE-1)/(NCYCLE-1)
C        DelLines = (1-u)*0.5*MAX(NPXT/Nxd,NPYT/Nyd) + u*1.0
C
C         ! Do NITER iterations per temp cycle
C        DO ITER = 1, NITER
C           ! Perturb the boundary lines
C          TrialLines(:)=CurLines(:)
C          DO WHILE (ALL(TrialLines(:)==CurLines(:)))
C            DO I = 1, N
C              TrialLines(I) = NINT(CurLines(I) + (2*ran(iseed)-1.0)*DelLines)
C            ENDDO
C          ENDDO
C           ! Calculate the energy of these new boundary lines
C          NewEnergy = LOAD_OBJ_FUNC (TrialLines, Nxd, Nyd, ProcCoordAll, wtproc, NPXT, NPYT, Ngrid)
C          DeltaE = NewEnergy - Energy
C          IF (DeltaE < 0.0) THEN   ! If lower energy accept change
C            CurLines(:) = TrialLines(:)
C            Energy = NewEnergy
C          ELSE                      !  otherwise sometimes accept it anyway
C            Boltz = EXP(-DeltaE/Temp)
C            IF (ran(iseed) < Boltz) THEN
C              CurLines(:) = TrialLines(:)
C              Energy = NewEnergy
C            ENDIF
C          ENDIF
C          IF (Energy < LowestEnergy) THEN
C            BndLines(:) = CurLines(:)
C            LowestEnergy = Energy
C          ENDIF
C        ENDDO
C        ! print '(I3,2(1X,F9.6),10(1X,I3))', CYCLE,Temp,LowestEnergy, BndLines(:)
C      ENDDO  ! end of temperature cycles
C
C      DEALLOCATE (Ngrid)
C    ENDIF
C
C    CALL MPI_BCAST (BndLines, N, MPI_INTEGER, 0, comm2d, ierr)
C    DEALLOCATE (ProcCoordAll)
C    IF (myproc==0 .AND. DoLoadBalancing) THEN
C      WRITE (6,'(A)') 'Optimum load balanced boundary lines in the property grid:'
C      WRITE (6,'(A,16I4)') '  X lines:',BndLines(1:Nxd-1)
C      WRITE (6,'(A,16I4)') '  Y lines:',BndLines(Nxd:N)
C    ENDIF
C  ENDIF
C
C   ! Make the subdomain position array for this processor from BndLines
C  IDX(1)=1 ; IDX(Nxd+1)=NPXT+1 ; IF (Nxd==1) IDX(Nxd+1)=NPXT
C  IDY(1)=1 ; IDY(Nyd+1)=NPYT+1 ; IF (Nyd==1) IDY(Nyd+1)=NPYT
C  IDX(2:Nxd)=BndLines(1:Nxd-1) ; IDY(2:Nyd)=BndLines(Nxd:N)
C  IXYPRP(1:2,1)=IDX(ProcCoord(1)+1:ProcCoord(1)+2)
C  IXYPRP(1:2,2)=IDY(ProcCoord(2)+1:ProcCoord(2)+2)
CEND SUBROUTINE OPTIMIZE_PROCESSOR_DOMAINS
C
C
C
CREAL FUNCTION LOAD_OBJ_FUNC (BndLines, Nxd, Nyd, ProcCoordAll, wtproc, &
C                             NPXT, NPYT, Ngrid)
C ! Calculates the objective function to minimize in doing the processor
C ! load balancing optimization.  The boundary line indices in X and Y to
C ! the property grid, which is the optimization control vector is input
C ! in BndLines.  The objective function is the normalized sum of the
C ! squares of the difference between the number of grid points (Ngrid)
C ! of each domain and the mean number, weighted by the relative processor
C ! times.
C  INTEGER, INTENT(IN) :: Nxd, Nyd, BndLines(Nxd+Nyd-2)
C  INTEGER, INTENT(IN) :: ProcCoordAll(2,0:Nxd*Nyd-1)
C  REAL,    INTENT(IN) :: wtproc(0:Nxd*Nyd)
C  INTEGER, INTENT(IN) :: NPXT, NPYT, Ngrid(NPXT+1,NPYT+1)
C  INTEGER :: N, I, J, IX, IY, IDX(Nxd+1), IDY(Nyd+1)
C  REAL    :: OBJA
C  REAL, SAVE :: sumN=-1.0, a
C
C  IF (sumN < 0) THEN
C    sumN = SUM(Ngrid(:,:)) ;  a = FLOAT(Nxd*Nyd)/sumN
C  ENDIF
C  N=Nxd+Nyd-2
C  IDX(1)=1 ; IDX(Nxd+1)=NPXT+1 ; IF (Nxd==1) IDX(Nxd+1)=NPXT
C  IDY(1)=1 ; IDY(Nyd+1)=NPYT+1 ; IF (Nyd==1) IDY(Nyd+1)=NPYT
C
C  OBJA = 0.0
C  DO J = 0, Nxd*Nyd-1
C    IDX(2:Nxd)=BndLines(1:Nxd-1) ; IDY(2:Nyd)=BndLines(Nxd:N)
C    IX = ProcCoordAll(1,J)+1 ; IY = ProcCoordAll(2,J)+1
C    OBJA = OBJA + (1.0/Nxd*Nyd) * &
C         (a*wtproc(J)*SUM(Ngrid(IDX(IX):IDX(IX+1),IDY(IY):IDY(IY+1))) - 1.0)**2
C  ENDDO
C  LOAD_OBJ_FUNC = OBJA
CEND FUNCTION LOAD_OBJ_FUNC
C
C
C
C
CREAL FUNCTION ran(idum)
C ! Returns uniform random deviate between 0.0 and 1.0 (exclusive of endpoints).
C ! Call ran with idum negative to initialize; thereafter, do not alter idum.
C  IMPLICIT NONE
C  INTEGER, PARAMETER :: K4B=SELECTED_INT_KIND(9)
C  INTEGER(K4B), INTENT(INOUT) :: idum
C  INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
C  REAL, SAVE :: am
C  INTEGER(K4B), SAVE :: ix=-1, iy=-1, k
C
C  if (idum <= 0 .or. iy < 0) then
C    am = nearest(1.0,-1.0)/IM
C    iy = ior(ieor(888889999,abs(idum)),1)
C    ix = ieor(777755555,abs(idum))
C    idum = abs(idum)+1
C  endif
C  ix = ieor(ix,ishft(ix,13))
C  ix = ieor(ix,ishft(ix,-17))
C  ix = ieor(ix,ishft(ix,5))
C  k = iy/IQ
C  iy = IA*(iy-k*IQ)-IR*k
C  if (iy < 0) iy=iy+IM
C  ran = am*ior(iand(IM,ieor(ix,iy)),1)
CEND FUNCTION ran





 SUBROUTINE BROADCAST_USER_INPUT (RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
                        INSAVEFILE, OUTSAVEFILE,                      &
                        NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,       &
                        IPFLAG, KDIST, DELTAM, GRIDTYPE,              &
                        SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                        GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,   &
                        ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,  &
                        MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS,   &
                        MAX_TOTAL_MB, SPLITTING_FACTOR,               &
                        NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
  ! Broadcasts all the user input parameters (except unneeded file names)
  ! from the master processor to the slave processors.
   USE SHDOM_MPI_DATA
   USE MPI
   IMPLICIT NONE
   INTEGER, INTENT(INOUT) :: NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
   INTEGER, INTENT(IN)    :: MAXOUT, MAXPAR
   INTEGER, INTENT(INOUT) :: MAXITER, NUMOUT
   LOGICAL, INTENT(INOUT) :: KDIST, DELTAM, ACCELFLAG
   REAL,    INTENT(INOUT) :: SOLARFLUX, SOLARMU, SOLARAZ
   REAL,    INTENT(INOUT) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
   REAL,    INTENT(INOUT) :: SOLACC, SPLITACC, SHACC, OUTPARMS(MAXPAR,MAXOUT)
   REAL,    INTENT(INOUT) :: MAX_TOTAL_MB, SPLITTING_FACTOR
   REAL,    INTENT(INOUT) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
   CHARACTER(LEN=1), INTENT(INOUT) :: SRCTYPE, GRIDTYPE, UNITS, OUTTYPES(*)
   CHARACTER(LEN=80), INTENT(INOUT) :: RUNNAME, PROPFILE, SFCFILE, CKDFILE
   CHARACTER(LEN=80), INTENT(INOUT) :: INSAVEFILE, OUTSAVEFILE

   INTEGER :: ierr, intbuf(10)
   LOGICAL :: logbuf(3)
   REAL    :: realbuf(16)
   CHARACTER(LEN=1) :: char1buf(3)
   CHARACTER(LEN=80) :: char80buf(6)
   CHARACTER(LEN=3) :: procstr

   IF (numproc <= 1) RETURN

   intbuf(:) = (/ NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG, MAXITER, NUMOUT /)
   CALL MPI_BCAST (intbuf, 10, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   NSTOKES=intbuf(1) ; NX=intbuf(2) ; NY=intbuf(3) ; NZ=intbuf(4)
   NMU=intbuf(5) ; NPHI=intbuf(6) ; BCFLAG=intbuf(7) ; IPFLAG=intbuf(8)
   MAXITER=intbuf(9) ; NUMOUT=intbuf(10)

   logbuf(:) = (/ KDIST, DELTAM, ACCELFLAG /)
   CALL MPI_BCAST (logbuf, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
   KDIST=logbuf(1) ; DELTAM=logbuf(2) ; ACCELFLAG=logbuf(3)

   realbuf(:) = (/ SOLARFLUX, SOLARMU, SOLARAZ, GNDTEMP, GNDALBEDO, SKYRAD, &
              WAVENO(1:2), WAVELEN, SOLACC, SPLITACC, SHACC, &
       MAX_TOTAL_MB, SPLITTING_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO /)
   CALL MPI_BCAST (realbuf, 16, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   SOLARFLUX=realbuf(1) ; SOLARMU=realbuf(2) ; SOLARAZ=realbuf(3)
   GNDTEMP=realbuf(4) ; GNDALBEDO=realbuf(5) ; SKYRAD=realbuf(6)
   WAVENO(1:2)=realbuf(7:8) ; WAVELEN=realbuf(9)
   SOLACC=realbuf(10) ; SPLITACC=realbuf(11) ; SHACC=realbuf(12)
   MAX_TOTAL_MB=realbuf(13) ; SPLITTING_FACTOR=realbuf(14)
   NUM_SH_TERM_FACTOR=realbuf(15) ; CELL_TO_POINT_RATIO=realbuf(16)

   char1buf(:) = (/ SRCTYPE, GRIDTYPE, UNITS /)
   CALL MPI_BCAST (char1buf, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   SRCTYPE=char1buf(1) ; GRIDTYPE=char1buf(2) ; UNITS=char1buf(3)

   char80buf(:) = (/ PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, RUNNAME /)
   CALL MPI_BCAST (char80buf, 6*80, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   PROPFILE=char80buf(1) ; SFCFILE=char80buf(2) ; CKDFILE=char80buf(3)
   INSAVEFILE=char80buf(4) ; OUTSAVEFILE=char80buf(5) ; RUNNAME=char80buf(6)

   CALL MPI_BCAST (OUTPARMS, MAXPAR*NUMOUT, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   CALL MPI_BCAST (OUTTYPES, NUMOUT*1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

   WRITE (procstr,'(I3.3)') myproc
   IF (INSAVEFILE(1:4) .NE. 'NONE') INSAVEFILE = TRIM(INSAVEFILE)//procstr
   IF (OUTSAVEFILE(1:4) .NE. 'NONE') OUTSAVEFILE = TRIM(OUTSAVEFILE)//procstr
 END SUBROUTINE BROADCAST_USER_INPUT
C
C
C
CSUBROUTINE BROADCAST_PROPERTY_SIZE (NPX, NPY, NPZ, DELX, DELY, &
C                                    NUMPHASE, MAXLEG, MAXPGL)
C ! Broadcasts the property file size information from the master processor
C ! to the slave processors.
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(INOUT) :: NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
f2py intent(in, out) :: NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
C  REAL,    INTENT(INOUT) :: DELX, DELY
f2py intent(in, out) :: DELX, DELY
C
C  INTEGER :: ierr, intbuf(6)
C  REAL    :: realbuf(2)
C
C  IF (numproc <= 1) RETURN
C
C  intbuf(:) = (/ NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL /)
C  CALL MPI_BCAST (intbuf, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
C  NPX=intbuf(1) ; NPY=intbuf(2) ; NPZ=intbuf(3)
C  NUMPHASE=intbuf(4) ; MAXLEG=intbuf(5) ; MAXPGL=intbuf(6)
C
C  realbuf(:) = (/ DELX, DELY /)
C  CALL MPI_BCAST (realbuf, 2, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
C  DELX=realbuf(1) ; DELY=realbuf(2)
CEND SUBROUTINE BROADCAST_PROPERTY_SIZE


C
CSUBROUTINE SCATTER_PROPERTIES (NPXT, NPYT, NPX, NPY, NPZ, &
C                          NSTLEG, NLEG, NUMPHASE, &
C                          ZLEVELS, MAXASYM, &
C                          TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT, &
C                          TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP, &
C                          NPART)
C ! Sends portions of the optical properties on the property grid to each
C ! of the processors.  Deals with the duplicating the grid points on the
C ! boundaries between processors.
C
C ! Modified to use NPART.
C
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(IN) :: NPXT, NPYT, NPX, NPY, NPZ, NSTLEG, NUMPHASE, NPART
f2py intent(in) :: NPXT, NPYT, NPX, NPY, NPZ, NSTLEG, NUMPHASE, NPART
C  INTEGER, INTENT(INOUT) :: NLEG
f2py intent(in, out) :: NLEG
C  REAL,    INTENT(INOUT) :: ZLEVELS(NPZ), MAXASYM
f2py intent(in, out) :: ZLEVELS, MAXASYM
C  REAL,    INTENT(IN) :: TEMPPT(NPZ,NPYT,NPXT,NPART), EXTINCTPT(NPZ,NPYT,NPXT,NPART)
f2py intent(in) :: TEMPPT, EXTINCTPT
C  REAL,    INTENT(IN) :: ALBEDOPT(NPZ,NPYT,NPXT,NPART), LEGENPT(*)
f2py intent(in) :: ALBEDOPT, LEGENPT
C  INTEGER*2, INTENT(IN) :: IPHASEPT(NPZ,NPYT,NPXT,NPART)
f2py intent(in) :: IPHASEPT
C  REAL,    INTENT(OUT) :: TEMPP(NPZ,NPY,NPX,NPART), EXTINCTP(NPZ,NPY,NPX,NPART)
f2py intent(out) :: TEMPP, EXTINCTP
C  REAL,    INTENT(OUT) :: ALBEDOP(NPZ,NPY,NPX,NPART), LEGENP(*)
f2py intent(out) :: ALBEDOP, LEGENP
C  INTEGER*2, INTENT(OUT) :: IPHASEP(NPZ,NPY,NPX,NPART)
f2py intent(out) :: IPHASEP
C
C  INTEGER :: SX, EX, NX, SY, EY, NY, N, MAXN, MaxNX, MaxNY
C  INTEGER :: iproc, irq, ierr
C  INTEGER, ALLOCATABLE :: IXYALL(:,:,:), requests(:), status(:,:)
C  REAL,    ALLOCATABLE :: realbuf(:), rbuf(:,:,:,:)
C  REAL,    ALLOCATABLE :: tempbuf(:,:), extbuf(:,:), albbuf(:,:)
C  INTEGER*2, ALLOCATABLE :: ibuf(:,:,:,:), iphbuf(:,:)
C
C   ! Broadcast the small stuff that is going to all processors
C  ALLOCATE (realbuf(NPZ+1))
C  realbuf(1:NPZ+1) = (/ ZLEVELS(1:NPZ), MAXASYM /)
C  CALL MPI_BCAST (realbuf, NPZ+1, MPI_REAL, 0, comm2d, ierr)
C  ZLEVELS(1:NPZ)=realbuf(1:NPZ) ; MAXASYM=realbuf(NPZ+1)
C
C   ! If we have a tabulated phase function, then broadcast to all processors
C  IF (NUMPHASE > 0) THEN
C    IF (myproc == 0) LEGENP(1:NUMPHASE*(NLEG+1)*NSTLEG) = LEGENPT(1:NUMPHASE*(NLEG+1)*NSTLEG)
C    CALL MPI_BCAST (NLEG, 1, MPI_INTEGER, 0, comm2d, ierr)
C    CALL MPI_BCAST (LEGENP, NSTLEG*(NLEG+1)*NUMPHASE, MPI_REAL, 0, comm2d, ierr)
C  ENDIF
C
C   ! Get the IXYPRP arrays from all processors to tell where they are in property grid
C  ALLOCATE (IXYALL(2,2,0:numproc-1))
C  CALL MPI_GATHER (IXYPRP,4,MPI_INTEGER, IXYALL,4,MPI_INTEGER, 0, comm2d, ierr)
C
C  ALLOCATE (requests(4*numproc+4), status(MPI_STATUS_SIZE,4*numproc+4))
C
C  irq = 0
C  IF (myproc == 0) THEN
C     ! Allocate the send buffers for each processor
C    MAXN=1 ; MaxNX=1 ; MaxNY=1
C    DO iproc = 0, numproc-1
C      NX=IXYALL(2,1,iproc)-IXYALL(1,1,iproc)+1 ; MaxNX=MAX(MaxNX,NX)
C      NY=IXYALL(2,2,iproc)-IXYALL(1,2,iproc)+1 ; MaxNY=MAX(MaxNY,NY)
C      MAXN = MAX(MAXN,NPZ*NY*NX)
C    ENDDO
C    ALLOCATE (rbuf(NPZ,MaxNY,MaxNX,NPART), ibuf(NPZ,MaxNY,MaxNX,NPART))
C    ALLOCATE (tempbuf(MAXN*NPART,0:numproc-1), extbuf(MAXN*NPART,0:numproc-1))
C    ALLOCATE (albbuf(MAXN*NPART,0:numproc-1), iphbuf(MAXN*NPART,0:numproc-1))
C     ! The master processor reformats the four property arrays (to deal with
C     ! the wrap around at the full domain edge) and sends the arrays with MPI.
C    DO iproc = 0, numproc-1
C      SX=IXYALL(1,1,iproc) ; EX=IXYALL(2,1,iproc)
C      SY=IXYALL(1,2,iproc) ; EY=IXYALL(2,2,iproc)
C      NX=EX-SX+1 ; NY=EY-SY+1 ; N=NX*NY*NPZ*NPART
C      rbuf(:,1:NY-1,1:NX-1,:) = TEMPPT(:,SY:EY-1,SX:EX-1,:)
C      rbuf(:,1:NY-1,NX,:) = TEMPPT(:,SY:EY-1,MOD(EX-1,NPXT)+1,:)
C      rbuf(:,NY,1:NX-1,:) = TEMPPT(:,MOD(EY-1,NPYT)+1,SX:EX-1,:)
C      rbuf(:,NY,NX,:) = TEMPPT(:,MOD(EY-1,NPYT)+1,MOD(EX-1,NPXT)+1,:)
C      tempbuf(1:N,iproc) = RESHAPE(rbuf(1:NPZ,1:NY,1:NX,:), (/ N /))
C      CALL MPI_ISEND (tempbuf(1:N,iproc), N, MPI_REAL, iproc, 1, comm2d, requests(irq+1), ierr)
C
C      rbuf(:,1:NY-1,1:NX-1,:) = EXTINCTPT(:,SY:EY-1,SX:EX-1,:)
C      rbuf(:,1:NY-1,NX,:) = EXTINCTPT(:,SY:EY-1,MOD(EX-1,NPXT)+1,:)
C      rbuf(:,NY,1:NX-1,:) = EXTINCTPT(:,MOD(EY-1,NPYT)+1,SX:EX-1,:)
C      rbuf(:,NY,NX,:) = EXTINCTPT(:,MOD(EY-1,NPYT)+1,MOD(EX-1,NPXT)+1,:)
C      extbuf(1:N,iproc) = RESHAPE(rbuf(1:NPZ,1:NY,1:NX,:), (/ N /))
C      CALL MPI_ISEND (extbuf(1:N,iproc), N*NPART, MPI_REAL, iproc, 2, comm2d, requests(irq+2), ierr)
C
C      rbuf(:,1:NY-1,1:NX-1,:) = ALBEDOPT(:,SY:EY-1,SX:EX-1,:)
C      rbuf(:,1:NY-1,NX,:) = ALBEDOPT(:,SY:EY-1,MOD(EX-1,NPXT)+1,:)
C      rbuf(:,NY,1:NX-1,:) = ALBEDOPT(:,MOD(EY-1,NPYT)+1,SX:EX-1,:)
C      rbuf(:,NY,NX,:) = ALBEDOPT(:,MOD(EY-1,NPYT)+1,MOD(EX-1,NPXT)+1,:)
C      albbuf(1:N,iproc) = RESHAPE(rbuf(1:NPZ,1:NY,1:NX,:), (/ N /))
C      CALL MPI_ISEND (albbuf(1:N,iproc), N, MPI_REAL, iproc, 3, comm2d, requests(irq+3), ierr)
C
C      ibuf(:,1:NY-1,1:NX-1,:) = IPHASEPT(:,SY:EY-1,SX:EX-1,:)
C      ibuf(:,1:NY-1,NX,:) = IPHASEPT(:,SY:EY-1,MOD(EX-1,NPXT)+1,:)
C      ibuf(:,NY,1:NX-1,:) = IPHASEPT(:,MOD(EY-1,NPYT)+1,SX:EX-1,:)
C      ibuf(:,NY,NX,:) = IPHASEPT(:,MOD(EY-1,NPYT)+1,MOD(EX-1,NPXT)+1,:)
C      iphbuf(1:N,iproc) = RESHAPE(ibuf(1:NPZ,1:NY,1:NX,:), (/ N /))
C      CALL MPI_ISEND (iphbuf(1:N,iproc), N, MPI_INTEGER2, iproc, 4, comm2d, requests(irq+4), ierr)
C      irq = irq + 4
C    ENDDO
C  ENDIF
C
C  N=NPX*NPY*NPZ*NPART
C  CALL MPI_IRECV (TEMPP, N, MPI_REAL, 0, 1, comm2d, requests(irq+1), ierr)
C  CALL MPI_IRECV (EXTINCTP, N, MPI_REAL, 0, 2, comm2d, requests(irq+2), ierr)
C  CALL MPI_IRECV (ALBEDOP, N, MPI_REAL, 0, 3, comm2d, requests(irq+3), ierr)
C  CALL MPI_IRECV (IPHASEP, N, MPI_INTEGER2, 0, 4, comm2d, requests(irq+4), ierr)
C  irq = irq + 4
C  CALL MPI_WAITALL (irq, requests, status, ierr)
C
C  IF (myproc == 0) DEALLOCATE (tempbuf, extbuf, albbuf, iphbuf, rbuf, ibuf)
C  DEALLOCATE (requests, status, IXYALL, realbuf)
CEND SUBROUTINE SCATTER_PROPERTIES




 SUBROUTINE BROADCAST_SURFACE_SIZE (MAXSFCPTS, MAXSFCPARS)
  ! Broadcasts the surface array sizes from the master process to the others.
   USE SHDOM_MPI_DATA
   USE MPI
   IMPLICIT NONE
   INTEGER, INTENT(INOUT) :: MAXSFCPTS, MAXSFCPARS
 !f2py intent(in, out) :: MAXSFCPTS, MAXSFCPARS
   INTEGER :: ierr, intbuf(2)

   IF (numproc <= 1) RETURN

   intbuf(1:2) = (/ MAXSFCPTS, MAXSFCPARS/)
   CALL MPI_BCAST (intbuf, 2, MPI_INTEGER, 0, comm2d, ierr)
   MAXSFCPTS=intbuf(1) ; MAXSFCPARS=intbuf(2)
 END SUBROUTINE BROADCAST_SURFACE_SIZE


 SUBROUTINE BROADCAST_SURFACE_PARMS (SFCTYPE, NXSFC, NYSFC, NSFCPAR, &
                                 DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)
  ! Broadcasts the surface file parameters from the master process to the others.
   USE SHDOM_MPI_DATA
   USE MPI
   IMPLICIT NONE
   CHARACTER(LEN=2), INTENT(INOUT) :: SFCTYPE
 !f2py intent(in, out) :: SFCTYPE
   INTEGER, INTENT(INOUT) :: NXSFC, NYSFC, NSFCPAR
 !f2py intent(in,out) :: NXSFC, NYSFC, NSFCPAR
   REAL,    INTENT(INOUT) :: DELXSFC, DELYSFC, SFCPARMS(*), GNDTEMP, GNDALBEDO
 !f2py intent(in, out) :: DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO
   INTEGER :: ierr, intbuf(3), n
   REAL    :: realbuf(4)


   IF (numproc <= 1) RETURN

   CALL MPI_BCAST (SFCTYPE, 2, MPI_CHARACTER, 0, comm2d, ierr)

   intbuf(1:3) = (/ NXSFC, NYSFC, NSFCPAR /)
   CALL MPI_BCAST (intbuf, 3, MPI_INTEGER, 0, comm2d, ierr)
   NXSFC=intbuf(1) ; NYSFC=intbuf(2) ; NSFCPAR=intbuf(3)

   realbuf(1:4) = (/ DELXSFC, DELYSFC, GNDTEMP, GNDALBEDO /)
   CALL MPI_BCAST (realbuf, 4, MPI_REAL, 0, comm2d, ierr)
   DELXSFC=realbuf(1) ; DELYSFC=realbuf(2)
   GNDTEMP=realbuf(3) ; GNDALBEDO=realbuf(4)

   n = NSFCPAR*(NXSFC+1)*(NYSFC+1)
   CALL MPI_BCAST (SFCPARMS, n, MPI_REAL, 0, comm2d, ierr)
 END SUBROUTINE BROADCAST_SURFACE_PARMS




 SUBROUTINE BROADCAST_KDIST_SIZE (NG, NZCKD)
  ! Broadcasts the k-distribution array sizes from the master process to the others.
   USE SHDOM_MPI_DATA
   USE MPI
   IMPLICIT NONE
   INTEGER, INTENT(INOUT) :: NG, NZCKD
 !f2py intent(in, out) :: NG, NZCKD
   INTEGER :: ierr, intbuf(2)

   IF (numproc <= 1) RETURN

   intbuf(1:2) = (/ NG, NZCKD/)
   CALL MPI_BCAST (intbuf, 2, MPI_INTEGER, 0, comm2d, ierr)
   NG=intbuf(1) ; NZCKD=intbuf(2)
 END SUBROUTINE BROADCAST_KDIST_SIZE


 SUBROUTINE BROADCAST_KDIST_PARMS (SOLFLUX, NG, DELG, NZCKD, ZCKD, KABS)
  ! Broadcasts the k-distribution file parameters from the master process
  ! to the other processes.
   USE SHDOM_MPI_DATA
   USE MPI
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: NG, NZCKD
 !f2py intent(in) :: NG, NZKD
   REAL,    INTENT(INOUT) :: SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD,NG)
 !f2py intent(in, out) :: SOLDFLUX, DELG, ZCKD, KABS
   INTEGER :: ierr
   REAL, ALLOCATABLE :: realbuf(:)

   IF (numproc <= 1) RETURN

   ALLOCATE (realbuf(NG+NZCKD+1))
   realbuf(1) = SOLFLUX ; realbuf(2:NG+1) = DELG(:)
   realbuf(NG+2:NG+NZCKD+1) = ZCKD(:)
   CALL MPI_BCAST (realbuf, NG+NZCKD+1, MPI_REAL, 0, comm2d, ierr)
   SOLFLUX = realbuf(1) ; DELG(:) = realbuf(2:NG+1)
   ZCKD(:) = realbuf(NG+2:NG+NZCKD+1)

   CALL MPI_BCAST (KABS, NZCKD*NG, MPI_REAL, 0, comm2d, ierr)
 END SUBROUTINE BROADCAST_KDIST_PARMS



 SUBROUTINE READ_BROADCAST_MEM_PARAMS (MAX_TOTAL_MB, ADAPT_GRID_FACTOR, &
                           NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, RUNNAME)
  ! The master process reads the RUNNAME//"_mem_params.inp" file (if there
  ! is one) and broadcasts the memory parameters to all the other processors.
   USE SHDOM_MPI_DATA
   USE MPI
   IMPLICIT NONE
   REAL, INTENT(INOUT) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
 !f2py intent(in, out) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
   REAL, INTENT(INOUT) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
 !f2py intent(in, out) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
   CHARACTER(LEN=*), INTENT(IN) :: RUNNAME
 !f2py intent(in) :: RUNNAME
   INTEGER :: Nproc, ierr, i, j
   REAL, ALLOCATABLE :: MemParam(:,:)

   ALLOCATE (MemParam(4,0:numproc-1))
   MemParam(1,:) = MAX_TOTAL_MB       ; MemParam(2,:) = ADAPT_GRID_FACTOR
   MemParam(3,:) = NUM_SH_TERM_FACTOR ; MemParam(4,:) = CELL_TO_POINT_RATIO

   IF (myproc == 0) THEN
      ! See if we can read the memory parameter file and number of processors matches
     OPEN (UNIT=11, FILE=TRIM(RUNNAME)//'_mem_params.inp', STATUS='OLD', IOSTAT=ierr)
     IF (ierr == 0) THEN
       READ (11,*)
       READ (11,*) Nproc
       IF (Nproc /= numproc) THEN
           WRITE (6,'(A,A,A)') 'Memory parameter file ',TRIM(RUNNAME),&
                '_mem_params.inp is for a different number of processors.'
       ELSE
         DO i = 0, Nproc-1
           READ (11,*) MemParam(:,i)
         ENDDO
       ENDIF
       CLOSE (11)
     ENDIF
   ENDIF

   CALL MPI_BCAST (MemParam, 4*numproc, MPI_REAL, 0, comm2d, ierr)
   MAX_TOTAL_MB = MemParam(1,myproc)
   ADAPT_GRID_FACTOR = MemParam(2,myproc)
   NUM_SH_TERM_FACTOR = MemParam(3,myproc)
   CELL_TO_POINT_RATIO = MemParam(4,myproc)
   DEALLOCATE (MemParam)
   WRITE (6,'(F7.2,1X,F7.3,2(1X,F5.3),A)') MAX_TOTAL_MB, ADAPT_GRID_FACTOR, &
      NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, '  Memory parameters used'
 END SUBROUTINE READ_BROADCAST_MEM_PARAMS



CREAL FUNCTION SUM_CPU_TIME (cpuTimes)
C ! Sums the CPU time used across all processors
C  USE SHDOM_MPI_DATA
C  USE MPI
C  real cpuTimes
C  integer ierr
C
C  IF (numproc>1) THEN
C    CALL MPI_REDUCE(cpuTimes, SUM_CPU_TIME, 1, MPI_REAL, MPI_SUM, 0, comm2d, ierr)
C  ELSE
C    SUM_CPU_TIME = cpuTimes
C  ENDIF
CEND FUNCTION SUM_CPU_TIME
C
C
C
CSUBROUTINE GATHER_OUTPUT (FLUX_OUT, FLUXDIV_OUT, SH_OUT, &
C                          IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, &
C                          SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS, &
C                          NPXT, NPYT, DELX, DELY, XALLGRID, YALLGRID, &
C                          ALLFLUXES, ALLFLUXDIV, ALLSHTERMS, &
C                          NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)
C ! Gathers the base grid output data in SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV,
C ! and SUMSHTERMS from all the processors and assembles in the master
C ! processor output arrays ALLFLUXES, ALLFLUXDIV, ALLSHTERMS.  Gathers the
C ! IXYOUT array from each processor to determine where each subdomain is.
C ! Also gathers the full base grid position arrays (XALLGRID/YALLGRID) and
C ! sums NCELLS,NPTS,NSH over all processors.
C  USE SHDOM_MPI_DATA
C  USE MPI
C  LOGICAL, INTENT(IN) :: FLUX_OUT, FLUXDIV_OUT, SH_OUT
C  INTEGER, INTENT(IN) :: IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, NPXT, NPYT
C  REAL,    INTENT(IN) :: DELX, DELY
C  REAL,    INTENT(IN) :: SUMFLUXES(2,NBPTS), SUMDIRFLUX(NBPTS)
C  REAL,    INTENT(IN) :: SUMFLUXDIV(NBPTS), SUMSHTERMS(NSHOUT,NBPTS)
C  REAL,    INTENT(OUT) :: XALLGRID(NXT), YALLGRID(NYT)
C  REAL,    INTENT(OUT) :: ALLFLUXES(3,NZ,NYT,NXT), ALLFLUXDIV(NZ,NYT,NXT)
C  REAL,    INTENT(OUT) :: ALLSHTERMS(NSHOUT,NZ,NYT,NXT)
C  INTEGER, INTENT(IN) :: NCELLS, NPTS, NSH
C  INTEGER, INTENT(OUT) :: NCELLSTOT, NPTSTOT, NSHTOT
C  INTEGER :: IX, IY, NX, NY, NX1, NY1, MAXN, I, send3int(3), recv3int(3)
C  INTEGER :: iproc, ierr, irq
C  INTEGER, ALLOCATABLE :: IXYALL(:,:,:), requests(:), status(:,:)
C  REAL,    ALLOCATABLE :: sendbuf(:,:), recvbuf(:,:,:), buf(:,:,:,:)
C
C  ALLOCATE (IXYALL(2,2,0:numproc-1))
C  CALL MPI_GATHER (IXYOUT,4,MPI_INTEGER, IXYALL,4,MPI_INTEGER, 0, comm2d, ierr)
C
C  ALLOCATE (requests(numproc), status(MPI_STATUS_SIZE,numproc))
C
C  IF (FLUX_OUT) THEN
C    irq=0
C    IF (myproc > 0) THEN
C      ALLOCATE (sendbuf(3,NBPTS))
C      sendbuf(1:2,:) = SUMFLUXES(:,:) ; sendbuf(3,:) = SUMDIRFLUX(:)
C      irq = irq + 1
C      CALL MPI_ISEND (sendbuf, 3*NBPTS, MPI_REAL, 0, 1, comm2d, &
C                      requests(irq), ierr)
C    ELSE
C      MAXN = (MAXVAL(IXYALL(2,1,:))+1)*(MAXVAL(IXYALL(2,2,:))+1)*NZ
C      ALLOCATE (recvbuf(3,MAXN,1:numproc-1))
C      DO iproc = 1, numproc-1
C        irq = irq + 1
C        CALL MPI_IRECV (recvbuf(:,:,iproc), 3*MAXN, MPI_REAL, iproc, 1, &
C                        comm2d, requests(irq), ierr)
C      ENDDO
C    ENDIF
C    CALL MPI_WAITALL (irq, requests, status, ierr)
C
C    IF (myproc == 0) THEN
C      IX = IXYALL(1,1,0) ; NX = IXYALL(2,1,0)
C      IY = IXYALL(1,2,0) ; NY = IXYALL(2,2,0)
C      NX1=NX+1 ; NY1=NY+1
C      IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
C      IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
C      ALLOCATE (buf(3,NZ,NY1,NX1))
C      buf(1:2,:,:,:) = RESHAPE(SUMFLUXES(1:2,:),(/ 2,NZ,NY1,NX1 /) )
C      buf(3,:,:,:) = RESHAPE(SUMDIRFLUX(:),(/ NZ,NY1,NX1 /) )
C      ALLFLUXES(1:3,1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(:,1:NZ,1:NY,1:NX)
C      DEALLOCATE (buf)
C      DO iproc = 1, numproc-1
C        IX = IXYALL(1,1,iproc) ; NX = IXYALL(2,1,iproc)
C        IY = IXYALL(1,2,iproc) ; NY = IXYALL(2,2,iproc)
C        NX1=NX+1 ; NY1=NY+1
C        IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
C        IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
C        ALLOCATE (buf(3,NZ,NY1,NX1))
C        buf(:,:,:,:) = RESHAPE(recvbuf(:,1:NZ*NY1*NX1,iproc),(/ 3,NZ,NY1,NX1 /) )
C        ALLFLUXES(1:3,1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(:,1:NZ,1:NY,1:NX)
C        DEALLOCATE (buf)
C      ENDDO
C      DEALLOCATE (recvbuf)
C    ELSE
C      DEALLOCATE (sendbuf)
C    ENDIF
C  ENDIF
C
C
C  IF (FLUXDIV_OUT) THEN
C    irq=0
C    IF (myproc > 0) THEN
C      irq = irq + 1
C      CALL MPI_ISEND (SUMFLUXDIV, NBPTS, MPI_REAL, 0, 2, comm2d, &
C                      requests(irq), ierr)
C    ELSE
C      MAXN = (MAXVAL(IXYALL(2,1,:))+1)*(MAXVAL(IXYALL(2,2,:))+1)*NZ
C      ALLOCATE (recvbuf(1,MAXN,1:numproc-1))
C      DO iproc = 1, numproc-1
C        irq = irq + 1
C        CALL MPI_IRECV (recvbuf(:,:,iproc), MAXN, MPI_REAL, iproc, 2, &
C                        comm2d, requests(irq), ierr)
C      ENDDO
C    ENDIF
C    CALL MPI_WAITALL (irq, requests, status, ierr)
C
C    IF (myproc == 0) THEN
C      IX = IXYALL(1,1,0) ; NX = IXYALL(2,1,0)
C      IY = IXYALL(1,2,0) ; NY = IXYALL(2,2,0)
C      NX1=NX+1 ; NY1=NY+1
C      IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
C      IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
C      ALLOCATE (buf(1,NZ,NY1,NX1))
C      buf(1,:,:,:) = RESHAPE(SUMFLUXDIV(:),(/ NZ,NY1,NX1 /) )
C      ALLFLUXDIV(1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(1,1:NZ,1:NY,1:NX)
C      DEALLOCATE (buf)
C      DO iproc = 1, numproc-1
C        IX = IXYALL(1,1,iproc) ; NX = IXYALL(2,1,iproc)
C        IY = IXYALL(1,2,iproc) ; NY = IXYALL(2,2,iproc)
C        NX1=NX+1 ; NY1=NY+1
C        IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
C        IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
C        ALLOCATE (buf(1,NZ,NY1,NX1))
C        buf(:,:,:,:) = RESHAPE(recvbuf(:,1:NZ*NY1*NX1,iproc),(/ 1,NZ,NY1,NX1 /) )
C        ALLFLUXDIV(1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(1,1:NZ,1:NY,1:NX)
C        DEALLOCATE (buf)
C      ENDDO
C      DEALLOCATE (recvbuf)
C    ENDIF
C  ENDIF
C
C
C  IF (SH_OUT) THEN
C    irq=0
C    IF (myproc > 0) THEN
C      irq = irq + 1
C      CALL MPI_ISEND (SUMSHTERMS, NSHOUT*NBPTS, MPI_REAL, 0, 3, comm2d, &
C                      requests(irq), ierr)
C    ELSE
C      MAXN = (MAXVAL(IXYALL(2,1,:))+1)*(MAXVAL(IXYALL(2,2,:))+1)*NZ
C      ALLOCATE (recvbuf(NSHOUT,MAXN,1:numproc-1))
C      DO iproc = 1, numproc-1
C        irq = irq + 1
C        CALL MPI_IRECV (recvbuf(:,:,iproc), NSHOUT*MAXN, MPI_REAL, iproc, 3, &
C                        comm2d, requests(irq), ierr)
C      ENDDO
C    ENDIF
C    CALL MPI_WAITALL (irq, requests, status, ierr)
C
C    IF (myproc == 0) THEN
C      IX = IXYALL(1,1,0) ; NX = IXYALL(2,1,0)
C      IY = IXYALL(1,2,0) ; NY = IXYALL(2,2,0)
C      NX1=NX+1 ; NY1=NY+1
C      IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
C      IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
C      ALLOCATE (buf(NSHOUT,NZ,NY1,NX1))
C      buf(:,:,:,:) = RESHAPE(SUMSHTERMS(:,:),(/ NSHOUT,NZ,NY1,NX1 /) )
C      ALLSHTERMS(:,1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(:,1:NZ,1:NY,1:NX)
C      DEALLOCATE (buf)
C      DO iproc = 1, numproc-1
C        IX = IXYALL(1,1,iproc) ; NX = IXYALL(2,1,iproc)
C        IY = IXYALL(1,2,iproc) ; NY = IXYALL(2,2,iproc)
C        NX1=NX+1 ; NY1=NY+1
C        IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
C        IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
C        ALLOCATE (buf(NSHOUT,NZ,NY1,NX1))
C        buf(:,:,:,:) = RESHAPE(recvbuf(:,1:NZ*NY1*NX1,iproc),(/ NSHOUT,NZ,NY1,NX1 /) )
C        ALLSHTERMS(:,1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(:,1:NZ,1:NY,1:NX)
C        DEALLOCATE (buf)
C      ENDDO
C      DEALLOCATE (recvbuf)
C    ENDIF
C  ENDIF
C
C
C   ! Make the X and Y grids for the full domain
C  IF (myproc == 0) THEN
C    XALLGRID(1:NXT) = (/ (NPXT*DELX*(I-1)/NXT, I=1, NXT) /)
C    YALLGRID(1:NYT) = (/ (NPYT*DELY*(I-1)/NYT, I=1, NYT) /)
C  ENDIF
C
C   ! Get the sum of the NCELLS, NPTS, and NSH over the processors
C  send3int = (/ NCELLS, NPTS, NSH /)
C  CALL MPI_REDUCE (send3int, recv3int, 3, MPI_INTEGER, MPI_SUM, 0, comm2d, ierr)
C  NCELLSTOT=recv3int(1) ; NPTSTOT=recv3int(2) ; NSHTOT=recv3int(3)
C
C  DEALLOCATE (IXYALL, requests, status)
CEND SUBROUTINE GATHER_OUTPUT
C
C
CSUBROUTINE TOTAL_ALBEDO_MAX (ALBMAX)
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  REAL,    INTENT(INOUT) :: ALBMAX
f2py intent(in, out) :: ALBMAX
C  INTEGER :: ierr
C  REAL    :: sendrbuf
C
C  IF (numproc > 1) THEN
C     ! Call MPI_ALLREDUCE to find the max ALBMAX over all processors
C    sendrbuf = ALBMAX
C    call MPI_ALLREDUCE (sendrbuf, ALBMAX, 1, MPI_REAL, MPI_MAX, comm2d, ierr)
C    if (ierr /= MPI_SUCCESS) CALL ABORT_SHDOM_MPI ('TOTAL_ALBEDO_MAX: MPI_ALLREDUCE error')
C  ENDIF
CEND SUBROUTINE TOTAL_ALBEDO_MAX
C
C
CSUBROUTINE UNIFY_SPLITTING (DOSPLIT, STARTSPLITACC)
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  LOGICAL, INTENT(INOUT) :: DOSPLIT
f2py intent(in, out) :: DOSPLIT
C  REAL,    INTENT(INOUT) :: STARTSPLITACC
f2py intent(in, out) :: STARTSPLITACC
C  INTEGER :: ierr(2)
C  LOGICAL :: sendlbuf
C  REAL    :: sendrbuf
C
C  IF (numproc > 1) THEN
C     ! Call MPI_ALLREDUCE to
C    sendrbuf = STARTSPLITACC
C    call MPI_ALLREDUCE (sendrbuf, STARTSPLITACC, 1, MPI_REAL, MPI_MAX, comm2d, ierr(1))
C    sendlbuf = DOSPLIT
C    call MPI_ALLREDUCE (sendlbuf, DOSPLIT, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr(2))
C    if (ANY(ierr(:) /= MPI_SUCCESS)) CALL ABORT_SHDOM_MPI ('UNIFY_SPLITTING: MPI_ALLREDUCE error')
C  ENDIF
CEND SUBROUTINE UNIFY_SPLITTING
C
C
CSUBROUTINE TOTAL_SPLITCRIT_MAX (SPLITCRIT)
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  REAL,    INTENT(INOUT) :: SPLITCRIT
f2py intent(in, out) :: SPLITCRIT
C  INTEGER :: ierr
C  REAL    :: sendrbuf
C
C  IF (numproc > 1) THEN
C     ! Call MPI_ALLREDUCE to find the max SPLITCRIT over all processors
C    sendrbuf = SPLITCRIT
C    call MPI_ALLREDUCE (sendrbuf, SPLITCRIT, 1, MPI_REAL, MPI_MAX, comm2d, ierr)
C    if (ierr /= MPI_SUCCESS) CALL ABORT_SHDOM_MPI ('TOTAL_SPLITCRIT_MAX: MPI_ALLREDUCE error')
C  ENDIF
CEND SUBROUTINE TOTAL_SPLITCRIT_MAX
C
C
C
C
CSUBROUTINE MAKE_DIRECT_PAR (SPT, NPTS, BCFLAG, IPFLAG, &
C                            DELTAM, ML, NSTLEG, NLEG, &
C                            SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS, &
C                            NX, XGRID, NY, YGRID,  DIRFLUX)
C ! Makes the direct beam solar flux for the NPTS gridpoints started with
C ! index SPT having X/Y/Z positions in GRIDPOS.
C ! DIRFLUX is set to F*exp(-tau_sun).
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(IN) :: SPT, NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, NX, NY
C  LOGICAL, INTENT(IN) :: DELTAM
C  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, SOLARAZ
C  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS), XGRID(NX+1), YGRID(NY+1)
C  REAL,    INTENT(OUT) :: DIRFLUX(NPTS)
C
C  LOGICAL :: VALIDBEAM
C  REAL    :: UNIFZLEV, UZLEV, DIRPATH
C  INTEGER :: MAXNBEAM, Lupx, Lupy, Ldownx, Ldowny, Nbeam
C  INTEGER :: SIDE, NINVALID(0:4)
C  INTEGER :: IP, I, J, L
C  INTEGER :: ierr, irq, requests(12), status(MPI_STATUS_SIZE,12), iproc
C  INTEGER :: intsendbuf(2), intrecvbuf(2)
C  LOGICAL :: VALIDRAD, BTEST, AllDone
C  REAL    :: XO, YO, ZO, X, Y, Z, XMAX, YMAX, EPS, DIRJUNK
C  INTEGER, ALLOCATABLE :: Ndone(:)
C  INTEGER, ALLOCATABLE :: BndDirIntInfo(:,:,:)
C  REAL, ALLOCATABLE    :: BndDirRealInfo(:,:,:)
C  INTEGER, ALLOCATABLE :: DirPtrOut(:,:), DirPtrIn(:)
C  REAL, ALLOCATABLE    :: DirAllOut(:,:), DirAllIn(:)
C
C
C  MAXNBEAM = 2*NPTS
C  ALLOCATE (BndDirRealInfo(4,MAXNBEAM,4), BndDirIntInfo(2,MAXNBEAM,4))
C  ALLOCATE (DirAllOut(MAXNBEAM,0:numproc-1), DirPtrOut(MAXNBEAM,0:numproc-1))
C  ALLOCATE (DirAllIn(MAXNBEAM), DirPtrIn(MAXNBEAM), Ndone(0:numproc-1))
C
C  EPS = 1.0E-5*MAX(XGRID(2)-XGRID(1),YGRID(2)-YGRID(1))
C
C  IF (SPT == 1) THEN
C     ! Initialize the direct beam calculation
C    CALL DIRECT_BEAM_PROP (1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG, &
C             DELTAM, ML, NSTLEG, NLEG, SOLARFLUX,SOLARMU,SOLARAZ, DIRFLUX(1), &
C             UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
C     ! Find the max lowest uniform Z level over all processors and then
C     !   tell DIRECT_BEAM_PROP about it
C    UZLEV = UNIFZLEV
C    CALL MPI_ALLREDUCE (UZLEV, UNIFZLEV, 1, MPI_REAL, MPI_MAX, comm2d, ierr)
C    CALL DIRECT_BEAM_PROP (2, 0.0, 0.0, 0.0, BCFLAG, IPFLAG, &
C             DELTAM,ML, NSTLEG, NLEG, SOLARFLUX,SOLARMU,SOLARAZ, DIRFLUX(1), &
C             UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
C  ENDIF
C
C
C  Lupx=0 ; Lupy=0 ; Ldownx=0 ; Ldowny=0
C  IF (BTEST(BCFLAG,2)) THEN
C    IF (COS(SOLARAZ) .GT. 1.0E-5) THEN
C      Ldownx = 2 ; Lupx = 1
C    ELSE
C      Ldownx = 1 ; Lupx = 2
C    ENDIF
C    XMAX = XGRID(NX)
C  ELSE
C    XMAX = XGRID(NX+1)
C  ENDIF
C  IF (BTEST(BCFLAG,3)) THEN
C    IF (SIN(SOLARAZ) .GT. 1.0E-5) THEN
C      Ldowny = 4 ; Lupy = 3
C    ELSE
C      Ldowny = 3 ; Lupy = 4
C    ENDIF
C    YMAX = YGRID(NY)
C  ELSE
C    YMAX = YGRID(NY+1)
C  ENDIF
C
C   ! Do the initial direct beam tracing starting at the base grid points
C  NINVALID(:)=0
C  DO IP = SPT, NPTS
C    DIRPATH = 0.0
C    CALL DIRECT_BEAM_PROP (0, GRIDPOS(1,IP), GRIDPOS(2,IP), GRIDPOS(3,IP), &
C                 BCFLAG, IPFLAG, DELTAM, ML, NSTLEG, NLEG, &
C                 SOLARFLUX,SOLARMU,SOLARAZ,  DIRFLUX(IP), &
C                 UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
C     ! If we don't have a good direct beam then save the information
C     !   needed to pass to two upstream processors
C    IF (.NOT. VALIDBEAM) THEN
C      IF (SIDE /= Lupx .AND. SIDE /= Lupy) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Bad SIDE 1')
C      NINVALID(SIDE) = NINVALID(SIDE) + 1
C      IF (NINVALID(SIDE) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: NINVALID(SIDE) exceeds MAXNBEAM 1')
C      BndDirRealInfo(:,NINVALID(SIDE),SIDE) = (/ XO, YO, ZO, DIRPATH /)
C      BndDirIntInfo(:,NINVALID(SIDE),SIDE) = (/ myproc, IP /)
C    ENDIF
C  ENDDO
C
C
C   ! Cycle over the passing of info between processors and doing more
C   !  partial direct beam integrations until all rays have been finished
C  Ndone(:) = 0
C  AllDone = .FALSE.
C  DO WHILE (.NOT. AllDone)
C     ! Do the send and receives of the number of incomplete rays and
C     !   the real and integer info needed about the radiance integrations
C    irq=0
C    IF (BTEST(BCFLAG,2)) THEN
C      CALL MPI_ISEND (NINVALID(Lupx), 1, MPI_INTEGER, &
C                      iprocneigh(Lupx), 1, comm2d, requests(irq+1), ierr)
C      CALL MPI_ISEND (BndDirRealInfo(:,:,Lupx), 4*NINVALID(Lupx), MPI_REAL, &
C                      iprocneigh(Lupx), 2, comm2d, requests(irq+2), ierr)
C      CALL MPI_ISEND (BndDirIntInfo(:,:,Lupx), 2*NINVALID(Lupx), MPI_INTEGER, &
C                      iprocneigh(Lupx), 3, comm2d, requests(irq+3), ierr)
C      irq=irq+3
C    ENDIF
C    IF (BTEST(BCFLAG,3)) THEN
C      CALL MPI_ISEND (NINVALID(Lupy), 1, MPI_INTEGER, &
C                      iprocneigh(Lupy), 4, comm2d, requests(irq+1), ierr)
C      CALL MPI_ISEND (BndDirRealInfo(:,:,Lupy), 4*NINVALID(Lupy), MPI_REAL, &
C                      iprocneigh(Lupy), 5, comm2d, requests(irq+2), ierr)
C      CALL MPI_ISEND (BndDirIntInfo(:,:,Lupy), 2*NINVALID(Lupy), MPI_INTEGER, &
C                      iprocneigh(Lupy), 6, comm2d, requests(irq+3), ierr)
C      irq=irq+3
C    ENDIF
C    IF (BTEST(BCFLAG,2)) THEN
C      CALL MPI_IRECV (NINVALID(Ldownx), 1, MPI_INTEGER, &
C                      iprocneigh(Ldownx), 1, comm2d, requests(irq+1), ierr)
C      CALL MPI_IRECV (BndDirRealInfo(:,:,Ldownx), 4*MAXNBEAM, MPI_REAL, &
C                      iprocneigh(Ldownx), 2, comm2d, requests(irq+2), ierr)
C      CALL MPI_IRECV (BndDirIntInfo(:,:,Ldownx), 2*MAXNBEAM, MPI_INTEGER, &
C                      iprocneigh(Ldownx), 3, comm2d, requests(irq+3), ierr)
C      irq=irq+3
C    ENDIF
C    IF (BTEST(BCFLAG,3)) THEN
C      CALL MPI_IRECV (NINVALID(Ldowny), 1, MPI_INTEGER, &
C                      iprocneigh(Ldowny), 4, comm2d, requests(irq+1), ierr)
C      CALL MPI_IRECV (BndDirRealInfo(:,:,Ldowny), 4*MAXNBEAM, MPI_REAL, &
C                      iprocneigh(Ldowny), 5, comm2d, requests(irq+2), ierr)
C      CALL MPI_IRECV (BndDirIntInfo(:,:,Ldowny), 2*MAXNBEAM, MPI_INTEGER, &
C                      iprocneigh(Ldowny), 6, comm2d, requests(irq+3), ierr)
C      irq=irq+3
C    ENDIF
C    CALL MPI_WAITALL (irq, requests, status, ierr)
C
C    NINVALID(Lupx)=0 ; NINVALID(Lupy)=0
C    IF (BTEST(BCFLAG,2)) THEN
C       ! Continue the backwards ray integrations from the X boundary
C      DO J = 1, NINVALID(Ldownx)
C        IF (Ldownx == 1) BndDirRealInfo(1,J,Ldownx) = XGRID(1)
C        IF (Ldownx == 2) BndDirRealInfo(1,J,Ldownx) = XGRID(NX)
C        X=BndDirRealInfo(1,J,Ldownx) ; Y=BndDirRealInfo(2,J,Ldownx)
C        Z=BndDirRealInfo(3,J,Ldownx) ; DIRPATH=BndDirRealInfo(4,J,Ldownx)
C        IF (Y < YGRID(1)-EPS) Y = YMAX ; IF (Y > YMAX+EPS) Y = YGRID(1)
C        CALL DIRECT_BEAM_PROP (0, X, Y, Z, BCFLAG, IPFLAG, &
C                               DELTAM, ML, NSTLEG, NLEG, &
C                               SOLARFLUX,SOLARMU,SOLARAZ,  DIRJUNK, &
C                               UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
C         ! If we got a valid direct beam then store it in DirAllOut otherwise
C         ! save the information to continue passing to neighboring processors
C        IF (VALIDBEAM) THEN
C          iproc = BndDirIntInfo(1,J,Ldownx)
C          Ndone(iproc) = Ndone(iproc) + 1
C          IF (Ndone(iproc) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Ndone(iproc) exceeds MAXNBEAM x')
C          DirAllOut(Ndone(iproc),iproc) = DIRPATH
C          DirPtrOut(Ndone(iproc),iproc) = BndDirIntInfo(2,J,Ldownx)
C        ELSE
C          IF (SIDE /= Lupx .AND. SIDE /= Lupy) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Bad SIDE 2x')
C          NINVALID(SIDE) = NINVALID(SIDE) + 1
C          IF (NINVALID(SIDE) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: NINVALID(SIDE) exceeds MAXNBEAM 2x')
C          BndDirRealInfo(:,NINVALID(SIDE),SIDE) = (/ XO, YO, ZO, DIRPATH /)
C          BndDirIntInfo(:,NINVALID(SIDE),SIDE) = BndDirIntInfo(:,J,Ldownx)
C        ENDIF
C      ENDDO
C    ENDIF
C
C    IF (BTEST(BCFLAG,3)) THEN
C       ! Continue the backwards ray integrations for the Y boundary
C      DO J = 1, NINVALID(Ldowny)
C        IF (Ldowny == 3) BndDirRealInfo(2,J,Ldowny) = YGRID(1)
C        IF (Ldowny == 4) BndDirRealInfo(2,J,Ldowny) = YGRID(NY)
C        X=BndDirRealInfo(1,J,Ldowny) ; Y=BndDirRealInfo(2,J,Ldowny)
C        Z=BndDirRealInfo(3,J,Ldowny) ; DIRPATH=BndDirRealInfo(4,J,Ldowny)
C        IF (X < XGRID(1)-EPS) X = XMAX ; IF (X > XMAX+EPS) X = XGRID(1)
C        CALL DIRECT_BEAM_PROP (0, X, Y, Z, BCFLAG, IPFLAG, &
C                               DELTAM, ML, NSTLEG, NLEG, &
C                               SOLARFLUX,SOLARMU,SOLARAZ,  DIRJUNK, &
C                               UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
C         ! If we got a valid direct beam then store it in DirAllOut otherwise
C         ! save the information to continue passing to neighboring processors
C        IF (VALIDBEAM) THEN
C          iproc = BndDirIntInfo(1,J,Ldowny)
C          Ndone(iproc) = Ndone(iproc) + 1
C          IF (Ndone(iproc) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Ndone(iproc) exceeds MAXNBEAM y')
C          DirAllOut(Ndone(iproc),iproc) = DIRPATH
C          DirPtrOut(Ndone(iproc),iproc) = BndDirIntInfo(2,J,Ldowny)
C        ELSE
C          IF (SIDE /= Lupx .AND. SIDE /= Lupy) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Bad SIDE 2y')
C          NINVALID(SIDE) = NINVALID(SIDE) + 1
C          IF (NINVALID(SIDE) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: NINVALID(SIDE) exceeds MAXNBEAM 2x')
C          BndDirRealInfo(:,NINVALID(SIDE),SIDE) = (/ XO, YO, ZO, DIRPATH /)
C          BndDirIntInfo(:,NINVALID(SIDE),SIDE) = BndDirIntInfo(:,J,Ldowny)
C        ENDIF
C      ENDDO
C    ENDIF
C
C     ! See if there are any more invalid radiance rays on all the processors
C    intsendbuf = (/ NINVALID(Lupx), NINVALID(Lupy) /)
C    CALL MPI_ALLREDUCE (intsendbuf, intrecvbuf, 2, MPI_INTEGER, MPI_MAX, &
C                        comm2d, ierr)
C    AllDone = ALL(intrecvbuf(1:2) == 0)
C    IF (MAXVAL(intrecvbuf(1:2)) > MAXNBEAM .OR. MAXVAL(Ndone(:)) > MAXNBEAM) THEN
C      CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: MAXNBEAM exceeded 2')
C    ENDIF
C  ENDDO
C
C
C   ! Exchange the finished direct beam paths with all of the other processors
C  DO iproc = 0, numproc-1
C    IF (iproc /= myproc) THEN
C      CALL MPI_ISEND (Ndone(iproc), 1, MPI_INTEGER, &
C                      iproc, 7, comm2d, requests(1), ierr)
C      CALL MPI_ISEND (DirAllOut(:,iproc), Ndone(iproc), MPI_REAL, &
C                      iproc, 8, comm2d, requests(2), ierr)
C      CALL MPI_ISEND (DirPtrOut(:,iproc), Ndone(iproc), MPI_INTEGER, &
C                      iproc, 9, comm2d, requests(3), ierr)
C      CALL MPI_IRECV (Nbeam, 1, MPI_INTEGER, &
C                      iproc, 7, comm2d, requests(4), ierr)
C      CALL MPI_IRECV (DirAllIn, MAXNBEAM, MPI_REAL, &
C                      iproc, 8, comm2d, requests(5), ierr)
C      CALL MPI_IRECV (DirPtrIn, MAXNBEAM, MPI_INTEGER, &
C                      iproc, 9, comm2d, requests(6), ierr)
C      CALL MPI_WAITALL (6, requests, status, ierr)
C       ! Put the direct beam flux in the correct grid point
C      DO I = 1, Nbeam
C        IP = DirPtrIn(I)
C        DIRFLUX(IP) = SOLARFLUX*EXP(-DirAllIn(I))
C      ENDDO
C    ELSE
C       ! If we have the direct paths for this processor then avoid MPI and
C       !   put the direct beam flux in the grid point
C      DO I = 1, Ndone(myproc)
C        IP = DirPtrOut(I,myproc)
C        DIRFLUX(IP) = SOLARFLUX*EXP(-DirAllOut(I,myproc))
C      ENDDO
C    ENDIF
C  ENDDO
C
C  DEALLOCATE (DirAllOut, DirPtrOut, DirAllIn, DirPtrIn, Ndone)
C  DEALLOCATE (BndDirRealInfo, BndDirIntInfo)
CEND SUBROUTINE MAKE_DIRECT_PAR
C
C
C
CSUBROUTINE FIND_BOUNDARY_POINTS (BCFLAG, IPFLAG, NPTS, SWEEPORD, GRIDPTR, &
C                                 GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)
C ! Makes the pointers to and positions of the horizontal boundary grid
C ! points of the processors neighboring this one.
C ! For each of the potentially eight octants there are
C ! two sets of boundary pointers; 1 is for the current subdomain,
C ! and 2 is for the neighboring subdomain.
C ! Fills in the BXYPTR(:,IBC,JOCT) and BXYPOS(:,IBC,JOCT) arrays
C ! (module variables), where IBC is a boundary point index and JOCT is
C ! a sort of octant number. BXYPTR holds (iproc,IPT,SIDE) for the boundary
C ! grid points from grid point IPT in processor iproc that are on the SIDE
C ! boundary of this subdomain (SIDE: 1=-X, 2=+X, 3=-Y, 4=+Y).  BXYPOS holds
C ! the gridpoint positions (X,Y,Z) for these boundary points.
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(IN) :: BCFLAG, IPFLAG, NPTS, SWEEPORD(NPTS,*), GRIDPTR(8,*)
C  INTEGER, INTENT(IN) :: NX, NY, NZ
C  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS), XGRID(NX), YGRID(NY), ZGRID(NZ)
C
C  INTEGER :: tag, ierr, irq, requests(24), status(MPI_STATUS_SIZE,24)
C  INTEGER :: MAXNBXY, NOCT, JOCT, IORDER, IPCELL, ICORNER, IPT, SIDE
C  INTEGER :: L, LI, N(4), N1(8,4), N2(8,4), NS, NBX, NBY
C  INTEGER :: MAXBC, IBC, BITZ, SZ, EZ, DZ, LX, LY, IX, IY, IZ, J
C  INTEGER :: INOUT(4) = (/ 2, 1, 4, 3 /)
C  INTEGER, SAVE :: OLDNPTS=0
C  LOGICAL :: SAMEGRID
C  REAL    :: XEDGE(8), YEDGE(8), XYB(4)
C  INTEGER, ALLOCATABLE :: BNDPTR(:,:,:,:), INDX(:)
C  REAL,    ALLOCATABLE :: BNDPOS(:,:,:,:,:), ZARRAY(:)
C
C  IF (.NOT. (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3))) RETURN
C  IF (BTEST(IPFLAG,0) .AND. BTEST(IPFLAG,1)) RETURN
C  CALL  MPI_ALLREDUCE (NPTS==OLDNPTS, SAMEGRID, 1, MPI_LOGICAL, MPI_LAND, comm2d, ierr)
C  OLDNPTS = NPTS
C  IF (SAMEGRID) RETURN
C
C
   IOCT=1+BITX+2*BITY+4*BITZ ; JOCT=JOCTORDER(IOCT) = (/ 1,3,5,7,2,4,6,8 /)
      CX<0: BITX=1, JOCT=3,7,4,8      CY<0: BITY=1, JOCT=5,7,6,8
C  XEDGE(:) = (/ XGRID(NX),XGRID(NX),XGRID(1),XGRID(1),XGRID(NX),XGRID(NX),XGRID(1),XGRID(1) /)
C  YEDGE(:) = (/ YGRID(NY),YGRID(NY),YGRID(NY),YGRID(NY),YGRID(1),YGRID(1),YGRID(1),YGRID(1) /)
C  XYB(1) = XGRID(1) ; XYB(2) = XGRID(NX)
C  XYB(3) = YGRID(1) ; XYB(4) = YGRID(NY)
C
C  IF (BTEST(IPFLAG,1) .AND. BTEST(IPFLAG,0)) THEN
C    NOCT = 2
C  ELSE IF (BTEST(IPFLAG,1)) THEN
C    NOCT = 4
C  ELSE
C    NOCT = 8
C  ENDIF
C     ! Allocate grid point pointer and position arrays (second from last
C     !   index is 1 for local, 2 for neighbor; last index is boundary L)
C  IF (ALLOCATED(IBXYPTR)) THEN
C    MAXNBXY = NINT(3.0*MAXVAL(IBXYPTR(2,:,:)))
C  ELSE
C    MAXNBXY = 3.0*(NZ-2)*MAX(NX,NY)
C  ENDIF
C  ALLOCATE (BNDPTR(NOCT,MAXNBXY,2,4), BNDPOS(3,NOCT,MAXNBXY,2,4))
C
C  N1(:,:) = 0
C  DO JOCT = 1, NOCT
C    BITZ = 1-IBITS(JOCT,0,1)
C    N(:) = 0
C    DO IORDER = 1, NPTS
C      IPCELL = ISHFT(SWEEPORD(IORDER,JOCT),-3)
C      ICORNER = IBITS(SWEEPORD(IORDER,JOCT),0,3)+1
C      IPT = GRIDPTR(ICORNER,IPCELL)
C
C      IF ( (GRIDPOS(3,IPT) < ZGRID(NZ) .AND. BITZ==0) .OR. &
C           (GRIDPOS(3,IPT) > ZGRID(1)  .AND. BITZ==1) ) THEN
C        IF (BTEST(BCFLAG,2) .AND. GRIDPOS(1,IPT)==XEDGE(JOCT) &
C            .AND. .NOT. BTEST(IPFLAG,0)) THEN
C          SIDE = 2-IBITS(JOCT-1,1,1)
C          N(SIDE) = N(SIDE) + 1
C          IF (N(SIDE)>MAXNBXY) then
C            CALL ABORT_SHDOM_MPI ('FIND_BOUNDARY_POINTS: MAXNBXY exceeded by N(SIDE) 1x')
C          ENDIF
C          BNDPTR(JOCT,N(SIDE),1,SIDE) = IPT
C          BNDPOS(:,JOCT,N(SIDE),1,SIDE) = GRIDPOS(:,IPT)
C        ENDIF
C        IF (BTEST(BCFLAG,3) .AND. GRIDPOS(2,IPT)==YEDGE(JOCT) &
C            .AND. .NOT. BTEST(IPFLAG,1)) THEN
C          SIDE = 4-IBITS(JOCT-1,2,1)
C          N(SIDE) = N(SIDE) + 1
C          IF (N(SIDE)>MAXNBXY) then
C            CALL ABORT_SHDOM_MPI ('FIND_BOUNDARY_POINTS: MAXNBXY exceeded by N(SIDE) 1y')
C          ENDIF
C          BNDPTR(JOCT,N(SIDE),1,SIDE) = IPT
C          BNDPOS(:,JOCT,N(SIDE),1,SIDE) = GRIDPOS(:,IPT)
C        ENDIF
C      ENDIF
C    ENDDO
C    N1(JOCT,:) = N(:)
C  ENDDO
C
C
C  NS = MAXVAL(N1(:,:))
C  IF (NS > MAXNBXY) CALL ABORT_SHDOM_MPI ('FIND_BOUNDARY_POINTS: MAXNBXY might be exceeded')
C  N2(:,:) = 0
C  irq = 0
C  DO L = 1, 4
C    IF ((BTEST(BCFLAG,2) .AND. (L.EQ.1 .OR. L.EQ.2)) .OR. &
C        (BTEST(BCFLAG,3) .AND. (L.EQ.3 .OR. L.EQ.4))) THEN
C       ! Exchange the lists of boundary points with neighbors
C      LI = INOUT(L)
C      tag = L
C      CALL MPI_ISEND (N1(:,L), 8, MPI_INTEGER, &
C                      iprocneigh(L), tag, comm2d, requests(irq+1), ierr)
C      CALL MPI_IRECV (N2(:,LI), 8, MPI_INTEGER, &
C                      iprocneigh(LI), tag, comm2d, requests(irq+2), ierr)
C      tag = 1000+L
C      CALL MPI_ISEND (BNDPTR(:,:,1,L), NOCT*NS, MPI_INTEGER, &
C                      iprocneigh(L), tag, comm2d, requests(irq+3), ierr)
C      CALL MPI_IRECV (BNDPTR(:,:,2,LI), NOCT*MAXNBXY, MPI_INTEGER, &
C                      iprocneigh(LI), tag, comm2d, requests(irq+4), ierr)
C      tag = 2000+L
C      CALL MPI_ISEND (BNDPOS(:,:,:,1,L), 3*NOCT*NS, MPI_REAL, &
C                      iprocneigh(L), tag, comm2d, requests(irq+5), ierr)
C      CALL MPI_IRECV (BNDPOS(:,:,:,2,LI), 3*NOCT*MAXNBXY, MPI_REAL, &
C                      iprocneigh(LI), tag, comm2d, requests(irq+6), ierr)
C      irq = irq + 6
C    ENDIF
C  ENDDO
C   ! Wait for the transfers among all the neighbors to resolve
C  CALL MPI_WAITALL (irq, requests, status, ierr)
C  NS = MAXVAL(N2(:,:))
C  IF (NS > MAXNBXY) CALL ABORT_SHDOM_MPI ('FIND_BOUNDARY_POINTS: MAXNBXY was exceeded')
C
C
C  IF (ALLOCATED(BXYPTR)) DEALLOCATE (BXYPTR, BXYPOS, IBXYPTR)
C  MAXBC = MAXVAL( (/ N2(:,1)+N2(:,3), N2(:,1)+N2(:,4), N2(:,2)+N2(:,3), N2(:,2)+N2(:,4) /) )
C  ALLOCATE (BXYPTR(3,MAXBC,NOCT), BXYPOS(3,MAXBC,NOCT), IBXYPTR(2,NZ,NOCT))
C
C    ! Put the boundary points from the neighboring subdomains in the
C    ! BXYPTR and BXYPOS arrays.  Merge in X and Y boundary points one Z slab
C    ! at a time, and set pointers (IBXYPTR) to start and end of each slab.
C    !  DO IBC = IBXYPTR(1,IZ,JOCT), IBXYPTR(2,IZ,JOCT)
C    !    X=BXYPOS(1,IBC,JOCT) ; Y=BXYPOS(2,IBC,JOCT) ; Z=BXYPOS(3,IBC,JOCT)
C    !    iproc=BXYPTR(1,IBC,JOCT) ; IPT=BXYPTR(2,IBC,JOCT) ; SIDE=BXYPTR(2,IBC,JOCT)
C
C  DO JOCT = 1, NOCT
C    BITZ = 1-IBITS(JOCT,0,1)
C    DZ = 2*BITZ-1
C    LX = 1+IBITS(JOCT-1,1,1)
C    LY = 3+IBITS(JOCT-1,2,1)
C     ! Combine the indices to the X and Y boundary points and sort on Z
C    NBX=0;  IF (BTEST(BCFLAG,2)) NBX=N2(JOCT,LX)
C    NBY=0;  IF (BTEST(BCFLAG,3)) NBY=N2(JOCT,LY)
C    ALLOCATE (ZARRAY(NBX+NBY), INDX(NBX+NBY))
C    ZARRAY(1:NBX) = BNDPOS(3,JOCT,1:NBX,2,LX)
C    INDX(1:NBX) = (/ (J, J=1,NBX) /)
C    ZARRAY(NBX+1:NBX+NBY) = BNDPOS(3,JOCT,1:NBY,2,LY)
C    INDX(NBX+1:NBX+NBY) = (/ (-J, J=1,NBY) /)
C    CALL SSORT (ZARRAY, INDX, NBX+NBY, 2*DZ)
C
C     ! For each Z slab use Z sorted index to transfer the X or Y boundary points
C    IBXYPTR(:,:,JOCT) = 0
C    SZ = (NZ-1)*(1-BITZ)+1+DZ
C    EZ = (NZ-1)*(BITZ)+1
C    IBC=0 ; J=1
C    DO IZ = SZ, EZ, DZ
C      IBXYPTR(1,IZ,JOCT) = IBC+1
C       ! Store the boundary point info for this Z slab
C      DO WHILE (J <= NBX+NBY .AND. &
C               ( (ZARRAY(J)>=ZGRID(IZ) .AND. BITZ==0) .OR. &
C                 (ZARRAY(J)<=ZGRID(IZ) .AND. BITZ==1) ))
C        IBC = IBC + 1
C        IF (INDX(J) > 0) THEN
C          IX = INDX(J)
C          BXYPTR(:,IBC,JOCT) = (/ iprocneigh(LX), BNDPTR(JOCT,IX,2,LX), LX /)
C          BXYPOS(:,IBC,JOCT) = BNDPOS(:,JOCT,IX,2,LX)
C          BXYPOS(INT((LX+1)/2),IBC,JOCT) = XYB(LX)
C        ELSE
C          IY = ABS(INDX(J))
C          BXYPTR(:,IBC,JOCT) = (/ iprocneigh(LY), BNDPTR(JOCT,IY,2,LY), LY /)
C          BXYPOS(:,IBC,JOCT) = BNDPOS(:,JOCT,IY,2,LY)
C          BXYPOS(INT((LY+1)/2),IBC,JOCT) = XYB(LY)
C        ENDIF
C        J = J + 1
C      ENDDO
C      IBXYPTR(2,IZ,JOCT) = IBC
C    ENDDO
C    DEALLOCATE (ZARRAY, INDX)
C  ENDDO
C
C  DEALLOCATE (BNDPTR, BNDPOS)
CEND SUBROUTINE FIND_BOUNDARY_POINTS
C
C
C
C
C
CSUBROUTINE CALC_BOUNDARY_RADIANCES (BCFLAG, IPFLAG, JOCT, IZ, &
C                                   NX, NY, NZ, XGRID, YGRID, ZGRID, &
C                                   NA, NPTS, NCELLS, GRIDPTR, &
C                                   NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS, &
C                                   MU, PHI, EXTINCT, NSTOKES, SOURCE, &
C                                   KANG, GRIDRAD)
C ! Calculates the discrete ordinate boundary radiances for the KANG discrete
C ! ordinate (going in direction MU/PHI) for the grid points in the IZ slab.
C ! FIND_BOUNDARY_POINTS has already been called to set up the BXYPTR and
C ! BXYPOS arrays for the grid points on the boundaries, and the IBXYPTR
C ! array which has the starting and ending points for the IZ slab.  The
C ! grid points are done in the sweeping order, and so the octant (JOCT)
C ! is used to index into the IBXYPTR, BXYPTR, and BXYPOS arrays.  The
C ! grid points specified in BXYPTR/BXYPOS are on the boundaries of this
C ! processor, and need to be traced back along the discrete ordinate
C ! first within this subdomain.  The radiative transfer integrations
C ! with the extinction (EXTINCT) and source function (SOURCE) are done
C ! with the INTEGRATE_RAY routine.  If the top or bottom boundary is
C ! hit when tracing back then the radiance is saved, otherwise the
C ! required information (X,Y,Z, transmission, partial radiance, original
C ! processor, and starting grid point) are saved to be sent to the
C ! appropriate neighboring processor.  These data are received from the
C ! neighboring processors and the radiative transfer integrations are
C ! continued, looping until all the rays have been finished.  Finally,
C ! the finished radiance and original grid points are exchanged with
C ! all the other processors.
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(IN) :: BCFLAG, IPFLAG, JOCT, IZ
C  INTEGER, INTENT(IN) :: NSTOKES, NX, NY, NZ, NA, NPTS, NCELLS, KANG
C  REAL,    INTENT(IN) :: XGRID(NX), YGRID(NY), ZGRID(NZ)
C  INTEGER, INTENT(IN) :: GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS)
C  INTEGER, INTENT(IN) :: TREEPTR(2,NCELLS)
C  INTEGER*2, INTENT(IN) :: CELLFLAGS(NCELLS)
C  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS)
C  REAL,    INTENT(IN) :: MU, PHI
C  REAL,    INTENT(IN) :: EXTINCT(NPTS), SOURCE(NSTOKES,NA,NPTS)
C  REAL,    INTENT(INOUT) :: GRIDRAD(NSTOKES,NPTS)
C
C  INTEGER :: MAXNBND, Nbcpnts, Lupx, Lupy, Ldownx, Ldowny
C  INTEGER :: SIDE, NINVALID(0:4)
C  INTEGER :: IBC, IPT, I, J, L
C  INTEGER :: ierr(12), irq, requests(12), status(MPI_STATUS_SIZE,12), iproc
C  INTEGER :: intsendbuf(2), intrecvbuf(2)
C  LOGICAL :: VALIDRAD, BTEST, AllDone
C  REAL    :: X, Y, Z, RADIANCE(NSTOKES), TRANSMIT, xo,yo,zo
C  INTEGER, ALLOCATABLE :: Ndone(:)
C  INTEGER, ALLOCATABLE :: BndRadIntInfo(:,:,:)
C  REAL, ALLOCATABLE    :: BndRadRealInfo(:,:,:)
C  INTEGER, ALLOCATABLE :: RadPtrOut(:,:), RadPtrIn(:)
C  REAL, ALLOCATABLE    :: RadAllOut(:,:,:), RadAllIn(:,:)
C
C
C  MAXNBND = numproc*MAXVAL(IBXYPTR(2,:,:)-IBXYPTR(1,:,:)+1)
C  ALLOCATE (BndRadRealInfo(4+NSTOKES,MAXNBND,4), BndRadIntInfo(2,MAXNBND,4))
C  ALLOCATE (RadAllOut(NSTOKES,MAXNBND,0:numproc-1), RadPtrOut(MAXNBND,0:numproc-1))
C  ALLOCATE (RadAllIn(NSTOKES,MAXNBND), RadPtrIn(MAXNBND), Ndone(0:numproc-1))
C
C  Lupx=0 ; Lupy=0 ; Ldownx=0 ; Ldowny=0
C  IF (BTEST(BCFLAG,2)) THEN
C    IF (COS(PHI) .GT. 1.0E-5) THEN
C      Ldownx = 2 ; Lupx = 1
C    ELSE
C      Ldownx = 1 ; Lupx = 2
C    ENDIF
C  ENDIF
C  IF (BTEST(BCFLAG,3)) THEN
C    IF (SIN(PHI) .GT. 1.0E-5) THEN
C      Ldowny = 4 ; Lupy = 3
C    ELSE
C      Ldowny = 3 ; Lupy = 4
C    ENDIF
C  ENDIF
C
C
C   ! Do the initial ray tracing
C  Ndone(:) = 0
C  NINVALID(:)=0
C   ! Loop over the downstream boundary points from the neighbor subdomains
C  DO IBC = IBXYPTR(1,IZ,JOCT), IBXYPTR(2,IZ,JOCT)
C     ! Call INTEGRATE_RAY for this
C    X=BXYPOS(1,IBC,JOCT) ; Y=BXYPOS(2,IBC,JOCT) ; Z=BXYPOS(3,IBC,JOCT)
C    iproc=BXYPTR(1,IBC,JOCT) ; IPT=BXYPTR(2,IBC,JOCT) ; SIDE=BXYPTR(3,IBC,JOCT)
C    RADIANCE(:) = 0.0 ; TRANSMIT = 1.0
C    xo=x; yo=y; zo=z
C    CALL INTEGRATE_RAY (BCFLAG, IPFLAG, NX, NY, NZ, NA, &
C                        XGRID, YGRID, ZGRID, NPTS, NCELLS, &
C                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                        GRIDPOS, EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD, &
C                        MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
C     ! If we got a radiance then store it otherwise save
C     !  the information needed to pass to other processors
C    IF (VALIDRAD) THEN
C      Ndone(iproc) = Ndone(iproc) + 1
C      RadAllOut(:,Ndone(iproc),iproc) = RADIANCE(:)
C      RadPtrOut(Ndone(iproc),iproc) = IPT
C    ELSE
C      IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
C        WRITE (6,*) 'Bad SIDE 1:',myproc,mu,phi,SIDE,x,y,z,xo,yo,zo
C        CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES')
C      ENDIF
C      NINVALID(SIDE) = NINVALID(SIDE) + 1
C      BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ X, Y, Z, TRANSMIT, RADIANCE(:) /)
C      BndRadIntInfo(:,NINVALID(SIDE),SIDE) = (/ iproc, IPT /)
C    ENDIF
C  ENDDO
C
C
C   ! Cycle over the passing of info between processors and doing more
C   !  partial radiance integrations until all rays have been finished
C  AllDone = .FALSE.
C  DO WHILE (.NOT. AllDone)
C     ! Do the send and receives of the number of incomplete rays and
C     !   the real and integer info needed about the radiance integrations
C    irq=0
C    IF (BTEST(BCFLAG,2)) THEN
C      CALL MPI_ISEND (NINVALID(Lupx), 1, MPI_INTEGER, &
C                      iprocneigh(Lupx), 1, comm2d, requests(irq+1), ierr(irq+1))
C      CALL MPI_ISEND (BndRadRealInfo(:,:,Lupx), (4+NSTOKES)*NINVALID(Lupx), MPI_REAL, &
C                      iprocneigh(Lupx), 2, comm2d, requests(irq+2), ierr(irq+2))
C      CALL MPI_ISEND (BndRadIntInfo(:,:,Lupx), 2*NINVALID(Lupx), MPI_INTEGER, &
C                      iprocneigh(Lupx), 3, comm2d, requests(irq+3), ierr(irq+3))
C      irq=irq+3
C    ENDIF
C    IF (BTEST(BCFLAG,3)) THEN
C      CALL MPI_ISEND (NINVALID(Lupy), 1, MPI_INTEGER, &
C                      iprocneigh(Lupy), 4, comm2d, requests(irq+1), ierr(irq+1))
C      CALL MPI_ISEND (BndRadRealInfo(:,:,Lupy), (4+NSTOKES)*NINVALID(Lupy), MPI_REAL, &
C                      iprocneigh(Lupy), 5, comm2d, requests(irq+2), ierr(irq+2))
C      CALL MPI_ISEND (BndRadIntInfo(:,:,Lupy), 2*NINVALID(Lupy), MPI_INTEGER, &
C                      iprocneigh(Lupy), 6, comm2d, requests(irq+3), ierr(irq+3))
C      irq=irq+3
C    ENDIF
C    IF (BTEST(BCFLAG,2)) THEN
C      CALL MPI_IRECV (NINVALID(Ldownx), 1, MPI_INTEGER, &
C                      iprocneigh(Ldownx), 1, comm2d, requests(irq+1), ierr(irq+1))
C      CALL MPI_IRECV (BndRadRealInfo(:,:,Ldownx), (4+NSTOKES)*MAXNBND, MPI_REAL, &
C                      iprocneigh(Ldownx), 2, comm2d, requests(irq+2), ierr(irq+2))
C      CALL MPI_IRECV (BndRadIntInfo(:,:,Ldownx), 2*MAXNBND, MPI_INTEGER, &
C                      iprocneigh(Ldownx), 3, comm2d, requests(irq+3), ierr(irq+3))
C      irq=irq+3
C    ENDIF
C    IF (BTEST(BCFLAG,3)) THEN
C      CALL MPI_IRECV (NINVALID(Ldowny), 1, MPI_INTEGER, &
C                      iprocneigh(Ldowny), 4, comm2d, requests(irq+1), ierr(irq+1))
C      CALL MPI_IRECV (BndRadRealInfo(:,:,Ldowny), (4+NSTOKES)*MAXNBND, MPI_REAL, &
C                      iprocneigh(Ldowny), 5, comm2d, requests(irq+2), ierr(irq+2))
C      CALL MPI_IRECV (BndRadIntInfo(:,:,Ldowny), 2*MAXNBND, MPI_INTEGER, &
C                      iprocneigh(Ldowny), 6, comm2d, requests(irq+3), ierr(irq+3))
C      irq=irq+3
C    ENDIF
C    IF (ANY(ierr(1:irq) /= MPI_SUCCESS)) THEN
C      WRITE (6,*) 'MPI error: ',ierr(1:irq)
C      CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES: MPI error 1')
C    ENDIF
C    CALL MPI_WAITALL (irq, requests, status, ierr(1))
C    IF (ierr(1) /= MPI_SUCCESS) THEN
C      WRITE (6,*) 'MPI_WAITALL error: ',ierr(1)
C      CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES: MPI error 2')
C    ENDIF
C
C    NINVALID(Lupx)=0 ; NINVALID(Lupy)=0
C    IF (BTEST(BCFLAG,2)) THEN
C       ! Continue the backwards ray integrations for the X boundary
C      DO J = 1, NINVALID(Ldownx)
C        IF (Ldownx == 1) BndRadRealInfo(1,J,Ldownx) = XGRID(1)
C        IF (Ldownx == 2) BndRadRealInfo(1,J,Ldownx) = XGRID(NX)
C        X=BndRadRealInfo(1,J,Ldownx) ; Y=BndRadRealInfo(2,J,Ldownx)
C        Z=BndRadRealInfo(3,J,Ldownx)
C        SIDE = Ldownx ; TRANSMIT=BndRadRealInfo(4,J,Ldownx)
C        RADIANCE(:)=BndRadRealInfo(5:4+NSTOKES,J,Ldownx)
C        CALL INTEGRATE_RAY (BCFLAG, IPFLAG, NX, NY, NZ, NA, &
C                        XGRID, YGRID, ZGRID, NPTS, NCELLS, &
C                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                        GRIDPOS, EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD, &
C                        MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
C
C         ! If we got a radiance then store it in RadAllOut (if for another
C         !   processor) otherwise save the information needed to continue
C         !   passing to neighboring processors
C        IF (VALIDRAD) THEN
C          iproc = BndRadIntInfo(1,J,Ldownx)
C          Ndone(iproc) = Ndone(iproc) + 1
C          RadAllOut(:,Ndone(iproc),iproc) = RADIANCE(:)
C          RadPtrOut(Ndone(iproc),iproc) = BndRadIntInfo(2,J,Ldownx)
C        ELSE
C          IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
C            WRITE (6,*) 'Bad SIDE 2x:',myproc,mu,phi,SIDE,x,y,z
C            CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES')
C          ENDIF
C          NINVALID(SIDE) = NINVALID(SIDE) + 1
C          BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ X, Y, Z, TRANSMIT, RADIANCE(:) /)
C          BndRadIntInfo(:,NINVALID(SIDE),SIDE) = BndRadIntInfo(:,J,Ldownx)
C        ENDIF
C      ENDDO
C    ENDIF
C
C    IF (BTEST(BCFLAG,3)) THEN
C       ! Continue the backwards ray integrations for the Y boundary
C      DO J = 1, NINVALID(Ldowny)
C        IF (Ldowny == 3) BndRadRealInfo(2,J,Ldowny) = YGRID(1)
C        IF (Ldowny == 4) BndRadRealInfo(2,J,Ldowny) = YGRID(NY)
C        X=BndRadRealInfo(1,J,Ldowny) ; Y=BndRadRealInfo(2,J,Ldowny)
C        Z=BndRadRealInfo(3,J,Ldowny)
C        SIDE = Ldowny ; TRANSMIT=BndRadRealInfo(4,J,Ldowny)
C        RADIANCE(:)=BndRadRealInfo(5:4+NSTOKES,J,Ldowny)
C        CALL INTEGRATE_RAY (BCFLAG, IPFLAG, NX, NY, NZ, NA, &
C                        XGRID, YGRID, ZGRID, NPTS, NCELLS, &
C                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                        GRIDPOS, EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD, &
C                        MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
C         ! If we got a radiance then store it in RadAllOut (if for another
C         !   processor) otherwise save the information needed to continue
C         !   passing to neighboring processors
C        IF (VALIDRAD) THEN
C          iproc = BndRadIntInfo(1,J,Ldowny)
C          Ndone(iproc) = Ndone(iproc) + 1
C          RadAllOut(:,Ndone(iproc),iproc) = RADIANCE(:)
C          RadPtrOut(Ndone(iproc),iproc) = BndRadIntInfo(2,J,Ldowny)
C        ELSE
C          IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
C            WRITE (6,*) 'Bad SIDE 2y:',myproc,mu,phi,SIDE,x,y,z
C            CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES')
C          ENDIF
C          NINVALID(SIDE) = NINVALID(SIDE) + 1
C          BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ X, Y, Z, TRANSMIT, RADIANCE(:) /)
C          BndRadIntInfo(:,NINVALID(SIDE),SIDE) = BndRadIntInfo(:,J,Ldowny)
C        ENDIF
C      ENDDO
C    ENDIF
C
C     ! See if there are any more invalid radiance rays on all the processors
C    intsendbuf = (/ NINVALID(Lupx), NINVALID(Lupy) /)
C    CALL MPI_ALLREDUCE (intsendbuf, intrecvbuf, 2, MPI_INTEGER, MPI_MAX, &
C                        comm2d, ierr(1))
C    AllDone = ALL(intrecvbuf(1:2) == 0)
C  ENDDO
C
C   ! Exchange the finished boundary radiances with all of the other processors
C  DO iproc = 0, numproc-1
C    IF (iproc /= myproc) THEN
C      CALL MPI_ISEND (Ndone(iproc), 1, MPI_INTEGER, &
C                      iproc, 7, comm2d, requests(1), ierr(1))
C      CALL MPI_ISEND (RadAllOut(:,:,iproc), NSTOKES*Ndone(iproc), MPI_REAL, &
C                      iproc, 8, comm2d, requests(2), ierr(2))
C      CALL MPI_ISEND (RadPtrOut(:,iproc), Ndone(iproc), MPI_INTEGER, &
C                      iproc, 9, comm2d, requests(3), ierr(3))
C      CALL MPI_IRECV (Nbcpnts, 1, MPI_INTEGER, &
C                      iproc, 7, comm2d, requests(4), ierr(4))
C      CALL MPI_IRECV (RadAllIn, NSTOKES*MAXNBND, MPI_REAL, &
C                      iproc, 8, comm2d, requests(5), ierr(5))
C      CALL MPI_IRECV (RadPtrIn, MAXNBND, MPI_INTEGER, &
C                      iproc, 9, comm2d, requests(6), ierr(6))
C      CALL MPI_WAITALL (6, requests, status, ierr(7))
C      IF (ANY(ierr(1:7) /= MPI_SUCCESS)) THEN
C        WRITE (6,*) 'MPI error: ',ierr(1:7)
C        CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES: MPI error 3')
C      ENDIF
C      DO I = 1, Nbcpnts
C        IPT = RadPtrIn(I)
C        GRIDRAD(:,IPT) = RadAllIn(:,I)
C      ENDDO
C    ELSE
C      DO I = 1, Ndone(myproc)
C        IPT = RadPtrOut(I,myproc)
C        GRIDRAD(:,IPT) = RadAllOut(:,I,myproc)
C      ENDDO
C    ENDIF
C  ENDDO
C
C  DEALLOCATE (RadAllOut, RadPtrOut, RadAllIn, RadPtrIn, Ndone)
C  DEALLOCATE (BndRadRealInfo, BndRadIntInfo)
CEND SUBROUTINE CALC_BOUNDARY_RADIANCES
C
C
C
C
CSUBROUTINE INTEGRATE_RAY (BCFLAG, IPFLAG, NX, NY, NZ, NA, &
C                          XGRID, YGRID, ZGRID, NPTS, NCELLS,&
C                          GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                          GRIDPOS, EXTINCT, NSTOKES, SOURCE, KANG, &
C                          GRIDRAD, MU, PHI, &
C                          X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
C ! Integrates the radiative transfer equation tracing a ray starting
C ! at X,Y,Z on boundary (face) SIDE to get the radiance in direction MU,PHI.
C ! The source function and extinction are integrated backwards from the
C ! starting point.  As each cell is crossed the integration of the source
C ! function across the cell is done and added to the radiance.  If there
C ! are known radiances at the end of the ray path in GRIDRAD, then the
C ! radiance there is interpolated from the four grid points on the face
C ! and VALIDRAD is returned true. Otherwise VALIDRAD is returned false.
C ! TRANSMIT is returned with the transmission from the ray starting point
C ! to the ending point, and RADIANCE has either the total Stokes radiance or
C ! the accumulated radiance until the boundary (TRANSMIT and RADIANCE must be
C ! properly initialized on input).  The ending location X,Y,Z and the
C ! subdomain SIDE (1=-X,2=+X,3=-Y,4=+Y) are also returned.
C  USE SHDOM_MPI_DATA
C  IMPLICIT NONE
C  INTEGER, INTENT(IN) :: BCFLAG, IPFLAG, NX, NY, NZ, NA, NPTS, NCELLS
C  INTEGER, INTENT(IN) :: NSTOKES, KANG
C  REAL,    INTENT(IN) :: XGRID(NX), YGRID(NY), ZGRID(NZ)
C  INTEGER, INTENT(IN) :: GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS)
C  INTEGER, INTENT(IN) :: TREEPTR(2,NCELLS)
C  INTEGER*2, INTENT(IN) ::  CELLFLAGS(NCELLS)
C  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS)
C  REAL,    INTENT(IN) :: EXTINCT(NPTS)
C  REAL,    INTENT(IN) :: SOURCE(NSTOKES,NA,NPTS), GRIDRAD(NSTOKES,NPTS)
C  REAL,    INTENT(IN) :: MU, PHI
C  REAL,    INTENT(INOUT) :: X, Y, Z
C  INTEGER, INTENT(INOUT) :: SIDE
C  REAL,    INTENT(INOUT) :: TRANSMIT, RADIANCE(NSTOKES)
C  LOGICAL, INTENT(OUT) :: VALIDRAD
C
C  INTEGER :: BITX, BITY, BITZ, IOCT, NXC, NYC, BTEST
C  INTEGER :: IX, IY, IZ, IL, IM, IU, IPTR, DIR, ICELL, INEXTCELL
C  INTEGER :: I1, I2, I3, I4, IOPP, IFACE, JFACE, KFACE, IC
C  INTEGER :: GRIDFACE(4,6), OPPFACE(6)
C  LOGICAL :: IPINX, IPINY
C  DOUBLE PRECISION PI, CX, CY, CZ, CXINV, CYINV, CZINV
C  DOUBLE PRECISION XE, YE, ZE
C  DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
C  DOUBLE PRECISION U, V, F1, F2, F3, F4
C  DOUBLE PRECISION EXT, EXT0, EXT1
C  DOUBLE PRECISION SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
C  DOUBLE PRECISION EXT0P, SRCEXT0P(NSTOKES)
C  DOUBLE PRECISION TAU, TRANSCELL, ABSCELL, TRANS, RAD(NSTOKES), RAD0(NSTOKES)
C  DATA GRIDFACE/1,3,5,7, 2,4,6,8,  1,2,5,6, 3,4,7,8, 1,2,3,4, 5,6,7,8/
C  DATA OPPFACE/2,1,4,3,6,5/
C
C
C  EPS = 1.0E-3*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
C
C   ! Make the ray direction (opposite to the discrete ordinate direction)
C  PI = ACOS(-1.0)
C  CX = SQRT(1.0-MU**2)*COS(PHI+PI)
C  CY = SQRT(1.0-MU**2)*SIN(PHI+PI)
C  CZ = -MU
C  IF (ABS(CX) .GT. 1.0E-5) THEN
C    CXINV = 1.0D0/CX
C  ELSE
C    CX = 0.0
C    CXINV = 1.0E6
C  ENDIF
C  IF (ABS(CY) .GT. 1.0E-5) THEN
C    CYINV = 1.0D0/CY
C  ELSE
C    CY = 0.0
C    CYINV = 1.0E6
C  ENDIF
C  CZINV = 1.0D0/CZ
C
C   ! Setup the sweeping direction and the gridpoint location
C   !   BITc is 0 for positive direction trace back ray, 1 for negative.
C  IF (CX .LT. 0.0) THEN
C    BITX = 1
C  ELSE
C    BITX = 0
C  ENDIF
C  IF (CY .LT. 0.0) THEN
C    BITY = 1
C  ELSE
C    BITY = 0
C  ENDIF
C  IF (CZ .LT. -1.0E-3) THEN
C    BITZ = 1
C  ELSE IF (CZ .GT. 1.0E-3) THEN
C    BITZ = 0
C  ELSE
C    CALL ABORT_SHDOM_MPI ('INTEGRATE_RAY: Bad MU')
C  ENDIF
C   ! IOCT is the octant the ray come from or the discrete ordinate goes to.
C  IOCT = 1 + BITX + 2*BITY + 4*BITZ
C
C   ! Find the base grid cell that the starting point is in; first in Z
C  IL=0
C  IU=NZ
C  DO WHILE (IU-IL .GT. 1)
C    IM = (IU+IL)/2
C    IF (Z .GE. ZGRID(IM)) THEN
C      IL = IM
C    ELSE
C      IU=IM
C    ENDIF
C  ENDDO
C  IZ = MAX(IL,1)
C
C   ! Then get the base grid cell in X and Y assuming evenly spaced grid
C  NXC = NX
C  IF (BTEST(BCFLAG,2) .AND. .NOT. BTEST(IPFLAG,0)) NXC=NX-1
C  NYC = NY
C  IF (BTEST(BCFLAG,3) .AND. .NOT. BTEST(IPFLAG,1)) NYC=NY-1
C  IF (BTEST(IPFLAG,0)) THEN
C    IX = 1+NINT((NX-1)*(X-XGRID(1))/(XGRID(NX)-XGRID(1)))
C  ELSE
C    IX = 1+INT((NX-1)*(X-XGRID(1))/(XGRID(NX)-XGRID(1)))
C  ENDIF
C  IF (BTEST(IPFLAG,1)) THEN
C    IY = 1+NINT((NY-1)*(Y-YGRID(1))/(YGRID(NY)-YGRID(1)))
C  ELSE
C    IY = 1+INT((NY-1)*(Y-YGRID(1))/(YGRID(NY)-YGRID(1)))
C  ENDIF
C  if (ix < 1 .or. iy < 1 .or. ix > NX+1 .or. iy > NY+1) then
C    WRITE (6,*) 'Bad ix,iy',x,y,z,ix,iy,iz,xgrid(1),xgrid(nx),ygrid(1),ygrid(ny),myproc
C    WRITE (6,*) 'Bad ix,iy',side,transmit
C    CALL ABORT_SHDOM_MPI ('INTEGRATE_RAY:')
C  endif
C  IX = MIN(NXC,IX)
C  IY = MIN(NYC,IY)
C
C
C   ! Make the base grid cell pointer
C  ICELL = IZ + (NZ-1)*(IY-1) + (NZ-1)*NYC*(IX-1)
C
C   ! Trace down the tree: point to the positive child cell, get the grid
C   ! point that is on the negative (X,Y,Z) side, use the flags to find
C   ! which dimension (X,Y,Z) the parent cell was split in, then compare
C   ! the grid point with the test point (X,Y,Z) to find which child cell
C   ! (pos or neg) the test point is in.
C  DO WHILE (TREEPTR(2,ICELL) .GT. 0)
C    DIR = IBITS(INT(CELLFLAGS(ICELL)),2,2)
C    IC = TREEPTR(2,ICELL) + 1
C    IPTR = GRIDPTR(1,IC)
C    IF (DIR .EQ. 1) THEN
C      IF (X .LT. GRIDPOS(1,IPTR))  IC = IC - 1
C    ELSE IF (DIR .EQ. 2) THEN
C      IF (Y .LT. GRIDPOS(2,IPTR))  IC = IC - 1
C    ELSE IF (DIR .EQ. 3) THEN
C      IF (Z .LT. GRIDPOS(3,IPTR))  IC = IC - 1
C    ENDIF
C    ICELL = IC
C  ENDDO
C
C  IFACE = SIDE
C   ! Get the four gridpoints for this face
C  I1 = GRIDPTR(GRIDFACE(1,IFACE),ICELL)
C  I2 = GRIDPTR(GRIDFACE(2,IFACE),ICELL)
C  I3 = GRIDPTR(GRIDFACE(3,IFACE),ICELL)
C  I4 = GRIDPTR(GRIDFACE(4,IFACE),ICELL)
C  IPINX = BTEST(INT(CELLFLAGS(ICELL)),0)
C  IPINY = BTEST(INT(CELLFLAGS(ICELL)),1)
C   ! Compute the face interpolation factors
C  IF (IFACE .LE. 2) THEN
C    U = (Z-GRIDPOS(3,I1))/(GRIDPOS(3,I3)-GRIDPOS(3,I1))
C    IF (IPINY) THEN
C      V = 0.5
C    ELSE
C      V = (Y-GRIDPOS(2,I1))/(GRIDPOS(2,I2)-GRIDPOS(2,I1))
C    ENDIF
C  ELSE
C    U = (Z-GRIDPOS(3,I1))/(GRIDPOS(3,I3)-GRIDPOS(3,I1))
C    IF (IPINX) THEN
C      V = 0.5
C    ELSE
C      V = (X-GRIDPOS(1,I1))/(GRIDPOS(1,I2)-GRIDPOS(1,I1))
C    ENDIF
C  ENDIF
C
C   ! Bilinearly interpolate the starting extinction and source function
C  F1 = (1-U)*(1-V) ; F2 = (1-U)*V ; F3 = U*(1-V) ; F4 = U*V
C  EXT1 = F1*EXTINCT(I1) + F2*EXTINCT(I2) + F3*EXTINCT(I3) + F4*EXTINCT(I4)
C   ! Correctly interpolate source using extinction*source
C  SRCEXT1(:) = F1*SOURCE(:,KANG,I1)*EXTINCT(I1) + F2*SOURCE(:,KANG,I2)*EXTINCT(I2) &
C             + F3*SOURCE(:,KANG,I3)*EXTINCT(I3) + F4*SOURCE(:,KANG,I4)*EXTINCT(I4)
C  EXT1 = MAX(0.0D0,EXT1)
C  SRCEXT1(1) = MAX(0.0D0,SRCEXT1(1))
C
C   ! Loop until finding a face with known radiances or reaching the boundary
C  XE = X  ; YE = Y  ; ZE = Z
C  TRANS = TRANSMIT
C  RAD(:) = RADIANCE(:)
C  VALIDRAD = .FALSE.
C  DO WHILE (.NOT. VALIDRAD .AND. ICELL .GT. 0)
C    IPINX = BTEST(INT(CELLFLAGS(ICELL)),0)
C    IPINY = BTEST(INT(CELLFLAGS(ICELL)),1)
C
C     ! Find boundaries of the current cell
C     !   Find the three possible intersection planes (X,Y,Z)
C     !   from the coordinates of the opposite corner grid point
C    IOPP = GRIDPTR(9-IOCT,ICELL)
C     ! Get the distances to the 3 planes and select the closest
C    IF (IPINX) THEN
C      SOX = 1.0E20
C    ELSE
C      SOX = (GRIDPOS(1,IOPP)-XE)*CXINV
C    ENDIF
C    IF (IPINY) THEN
C      SOY = 1.0E20
C    ELSE
C      SOY = (GRIDPOS(2,IOPP)-YE)*CYINV
C    ENDIF
C    SOZ = (GRIDPOS(3,IOPP)-ZE)*CZINV
C    SO = MIN(SOX,SOY,SOZ)
C    IF (SO .LT. -EPS) THEN
C      WRITE (6,*) 'INTEGRATE_RAY: SO<0  ',MU,PHI,X,Y,Z,SOX,SOY,SOZ,SO,ICELL,TRANS
C      CALL ABORT_SHDOM_MPI ('INTEGRATE_RAY: SO<0')
C    ENDIF
C     ! Compute the coordinates of the cell exitting location
C    IF (.NOT. IPINX) XE = XE + SO*CX
C    IF (.NOT. IPINY) YE = YE + SO*CY
C    ZE = ZE + SO*CZ
C
C     ! Get the intersection face number (i.e. neighptr index)
C    IF (SOX .LE. SOZ .AND. SOX .LE. SOY) THEN
C      IFACE = 2-BITX
C      JFACE = 1
C    ELSE IF (SOY .LE. SOZ) THEN
C      IFACE = 4-BITY
C      JFACE = 2
C    ELSE
C      IFACE = 6-BITZ
C      JFACE = 3
C    ENDIF
C     ! Get the next cell to go to
C    INEXTCELL = NEIGHPTR(IFACE,ICELL)
C    IF (INEXTCELL .LT. 0) THEN
C      CALL NEXT_CELL (XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS, &
C                      GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
C    ENDIF
C
C     ! Get the face grid pointers
C    IF (NEIGHPTR(IFACE,ICELL) .GE. 0) THEN
C       ! If going to same or larger face then use previous face
C      KFACE = IFACE
C      IC = ICELL
C    ELSE
C       ! If going to smaller face then use next face (more accurate)
C      KFACE = OPPFACE(IFACE)
C      IC = INEXTCELL
C    ENDIF
C    I1 = GRIDPTR(GRIDFACE(1,KFACE),IC)
C    I2 = GRIDPTR(GRIDFACE(2,KFACE),IC)
C    I3 = GRIDPTR(GRIDFACE(3,KFACE),IC)
C    I4 = GRIDPTR(GRIDFACE(4,KFACE),IC)
C     ! Compute the face interpolation factors
C    IF (JFACE .EQ. 1) THEN
C      U = (ZE-GRIDPOS(3,I1))/(GRIDPOS(3,I3)-GRIDPOS(3,I1))
C      IF (IPINY) THEN
C        V = 0.5
C      ELSE
C        V = (YE-GRIDPOS(2,I1))/(GRIDPOS(2,I2)-GRIDPOS(2,I1))
C      ENDIF
C    ELSE IF (JFACE .EQ. 2) THEN
C      U = (ZE-GRIDPOS(3,I1))/(GRIDPOS(3,I3)-GRIDPOS(3,I1))
C      IF (IPINX) THEN
C        V = 0.5
C      ELSE
C        V = (XE-GRIDPOS(1,I1))/(GRIDPOS(1,I2)-GRIDPOS(1,I1))
C      ENDIF
C    ELSE
C      IF (IPINY) THEN
C        U = 0.5
C      ELSE
C        U = (YE-GRIDPOS(2,I1))/(GRIDPOS(2,I3)-GRIDPOS(2,I1))
C      ENDIF
C      IF (IPINX) THEN
C        V = 0.5
C      ELSE
C        V = (XE-GRIDPOS(1,I1))/(GRIDPOS(1,I2)-GRIDPOS(1,I1))
C      ENDIF
C    ENDIF
C    if (u < -1.0E-4 .or. v < -1.0E-4 .or. u > 1.0001 .or. v > 1.0001) then
C      WRITE (6,*) 'u,v<0 or u,v>1: ',mu,phi,x,y,z,xe,ye,ze,u,v,ipinx,ipiny
C      CALL ABORT_SHDOM_MPI ('INTEGRATE_RAY:')
C    endif
C
C     ! Get the location coordinate (does the boundary wrapping)
C    IF (INEXTCELL .GT. 0) THEN
C      IF (JFACE .EQ. 1) THEN
C        XE = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
C      ELSE IF (JFACE .EQ. 2) THEN
C        YE = GRIDPOS(2,GRIDPTR(IOCT,INEXTCELL))
C      ELSE
C        ZE = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
C      ENDIF
C    ENDIF
C
C     ! Interpolate extinction and source function at face intersection
C    F1 = (1-U)*(1-V) ; F2 = (1-U)*V ; F3 = U*(1-V) ; F4 = U*V
C    EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2) + F3*EXTINCT(I3) + F4*EXTINCT(I4)
C     ! Correctly interpolate source using extinction*source
C    SRCEXT0(:) = F1*SOURCE(:,KANG,I1)*EXTINCT(I1) + F2*SOURCE(:,KANG,I2)*EXTINCT(I2) &
C               + F3*SOURCE(:,KANG,I3)*EXTINCT(I3) + F4*SOURCE(:,KANG,I4)*EXTINCT(I4)
C    EXT0 = MAX(0.0D0,EXT0)
C    SRCEXT0(1) = MAX(0.0D0,SRCEXT0(1))
C     ! Compute the cell radiance: integration of the source function
C    EXT = 0.5*(EXT0+EXT1)
C    TAU=EXT*SO
C    IF (TAU .GE. 0.5) THEN
C      TRANSCELL = EXP(-TAU)
C      ABSCELL = 1.0 - TRANSCELL
C    ELSE
C      ABSCELL = TAU*(1.0-0.5*TAU*(1.0-0.33333333333*TAU*(1-0.25*TAU)))
C      TRANSCELL = 1.0 - ABSCELL
C    ENDIF
C    IF (TAU .LE. 2.0) THEN
C      IF (EXT .EQ. 0.0) THEN
C        SRC(:) = 0.0
C      ELSE
C         ! Linear extinction, linear source*extinction, to first order
C        SRC(:) = ( 0.5*(SRCEXT0(:)+SRCEXT1(:)) &
C                + 0.08333333333*(EXT0*SRCEXT1(:)-EXT1*SRCEXT0(:))*SO )/EXT
C      ENDIF
C    ELSE
C       ! Combined first order expansion and constant extinction formula
C      EXT0P = EXT0
C      SRCEXT0P(:) = SRCEXT0(:)
C      IF (TAU .GT. 4.0) THEN
C        EXT0P = EXT1 + (EXT0-EXT1)*4.0/TAU
C        IF (EXT0 .GT. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
C      ENDIF
C      SRC(:) = 1.0/(EXT0P+EXT1) *( SRCEXT0P(:)+SRCEXT1(:) &
C             + (EXT0P*SRCEXT1(:)-EXT1*SRCEXT0P(:))*2.0/(EXT0P+EXT1) &
C                 *(1-2/TAU+2*TRANSCELL/ABSCELL) )
C    ENDIF
C    SRC(1) = MAX(SRC(1),0.0D0)
C
C     ! Add in the cell radiance and update the transmission.
C     !  If this is a boundary or the transmission is below the
C     !  cutoff and we have a valid radiance then set the flag
C     !  to stop the tracing and add in the interpolated face radiance.
C    RAD(:) = RAD(:) + TRANS*SRC(:)*ABSCELL
C    TRANS = TRANS*TRANSCELL
C    if (RAD(1) < -1.0E-5) then
C      WRITE (6,'(3(1X,F5.3),8(1X,E12.5))'), &
C        xe,ye,ze,ext0,ext1,so,tau,trans,src,abscell,rad(:)
C    endif
C    IF (GRIDRAD(1,I1).GE.0.0 .AND. GRIDRAD(1,I2).GE.0.0 &
C        .AND. GRIDRAD(1,I3).GE.0.0 .AND. GRIDRAD(1,I4).GE.0.0) THEN
C      VALIDRAD = .TRUE.
C      RAD0(:) = F1*GRIDRAD(:,I1) + F2*GRIDRAD(:,I2) &
C              + F3*GRIDRAD(:,I3) + F4*GRIDRAD(:,I4)
C      RAD(:) = RAD(:) + TRANS*RAD0(:)
C    ELSE IF (TRANS .LT. 1.0E-5) THEN
C      VALIDRAD = .TRUE.
C    ELSE
C      EXT1 = EXT0
C      SRCEXT1(:) = SRCEXT0(:)
C      ICELL = INEXTCELL
C    ENDIF
C  ENDDO
C  X = XE ; Y = YE ; Z = ZE
C  RADIANCE(:) = RAD(:)
C  TRANSMIT = TRANS
C  SIDE = IFACE
CEND SUBROUTINE INTEGRATE_RAY
C
C
C
C
C
CSUBROUTINE COMPUTE_RADIANCE_PAR (NSTOKES, NX, NY, NZ, NPTS, NCELLS, &
C                  ML, MM, NLM, NSTLEG, NLEG, NUMPHASE, &
C                  NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
C                  BCFLAG, XDOMAIN, YDOMAIN, IPFLAG, &
C                  SRCTYPE, DELTAM, SOLARMU, SOLARAZ, &
C                  SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
C                  MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
C                  GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS, &
C                  XGRID, YGRID, ZGRID, GRIDPOS, &
C                  GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                  EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, &
C                  SHPTR, SOURCE, SOURCE1, GRIDRAD, &
C                  OUTPARMS,  NRAD, RADOUT)
C ! Computes the radiances for the locations and directions specified
C ! in OUTPARMS, putting the results in RADOUT and increasing NRAD by
C ! the number of radiances output.  Does the radiative transfer integrations
C ! across multiple processors.
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(IN) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
C  INTEGER, INTENT(IN) :: ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
C  INTEGER, INTENT(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
C  INTEGER, INTENT(INOUT) :: NRAD
C  INTEGER, INTENT(IN) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
C  INTEGER, INTENT(IN) :: GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
C  INTEGER, INTENT(IN) :: SHPTR(NPTS+1), BCPTR(MAXNBC,2)
C  INTEGER*2, INTENT(IN) :: CELLFLAGS(NCELLS), IPHASE(NPTS)
C  LOGICAL, INTENT(IN) :: DELTAM
C  REAL, INTENT(IN) :: SOLARMU, SOLARAZ
C  REAL, INTENT(IN) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
C  REAL, INTENT(IN) :: MU(NMU), PHI(NMU,NPHI0MAX), WTDO(NMU,NPHI0MAX)
C  REAL, INTENT(IN) :: XDOMAIN, YDOMAIN, XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)
C  REAL, INTENT(IN) :: GRIDPOS(3,NPTS)
C  REAL, INTENT(IN) :: SFCGRIDPARMS(*), BCRAD(NSTOKES,*)
C  REAL, INTENT(IN) :: EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG,0:NLEG,NPTS)
C  REAL, INTENT(IN) :: DIRFLUX(NPTS), FLUXES(2,NPTS), SOURCE(NSTOKES,*)
C  REAL, INTENT(INOUT) :: SOURCE1(NSTOKES,NPTS), GRIDRAD(NSTOKES,NPTS)
C  REAL, INTENT(IN) :: OUTPARMS(*)
C  REAL, INTENT(OUT) :: RADOUT(NSTOKES,*)
C  CHARACTER, INTENT(IN) ::  SRCTYPE*1, SFCTYPE*2, UNITS*1
C
C  INTEGER :: IBC, JX, JY, K
C  INTEGER :: NANGOUT, NXOUT, NYOUT
C  LOGICAL :: LAMBERTIAN, VALIDRAD
C  DOUBLE PRECISION :: X0,Y0,Z0, XE,YE,ZE, X,Y,Z, TRANSMIT, RADIANCE(NSTOKES)
C  REAL    :: MUOUT, PHIOUT, PHID
C  REAL    :: STARTX, STARTY, EPS=1.0E-5
C
C  INTEGER :: MAXNRAD, NDOWN, Lupx, Lupy, Ldownx, Ldowny
C  INTEGER :: SIDE, NINVALID(0:4)
C  INTEGER :: I, J, L, Ndone, intsendbuf(2), intrecvbuf(2)
C  INTEGER :: ierr, irq, iproc
C  LOGICAL :: AllDone
C  INTEGER, ALLOCATABLE :: requests(:), status(:,:), NradRecv(:)
C  INTEGER, ALLOCATABLE :: IndexSend(:), IndexRecv(:,:), BndRadIntInfo(:,:)
C  REAL, ALLOCATABLE    :: RadSend(:,:), RadRecv(:,:,:), BndRadRealInfo(:,:,:)
C
C
C   ! Make the isotropic radiances for the top boundary
C  CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN, &
C                              UNITS, NTOPPTS, NSTOKES, BCRAD(1,1))
C   ! Make the bottom boundary radiances for the Lambertian surfaces.
C   ! Compute the upwelling bottom radiances using the downwelling fluxes.
C  IF (SFCTYPE .EQ. 'FL') THEN
C    CALL FIXED_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2), &
C                  DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO, &
C                  WAVENO, WAVELEN, UNITS,  NSTOKES, BCRAD(1,1+NTOPPTS))
C  ELSE IF (SFCTYPE .EQ. 'VL') THEN
C    CALL VARIABLE_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2), &
C                    DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS, &
C                    NSTOKES, BCRAD(1,1+NTOPPTS))
C  ENDIF
C
C
C   ! Setup the regularly spaced radiance locations
C  STARTX = 0.0
C  IF (OUTPARMS(2) .LE. 0.0) THEN
C    NXOUT = 1
C  ELSE
C    NXOUT = MAX(1,NINT(XDOMAIN/OUTPARMS(2)))
C  ENDIF
C  STARTY = 0.0
C  IF (OUTPARMS(3) .LE. 0.0) THEN
C    NYOUT = 1
C  ELSE
C    NYOUT = MAX(1,NINT(YDOMAIN/OUTPARMS(3)))
C  ENDIF
C  Z0 = MIN( MAX(OUTPARMS(1),ZGRID(1)), ZGRID(NZ))
C  NANGOUT = NINT(OUTPARMS(6))
C
C
C  LAMBERTIAN = SFCTYPE(2:2) .EQ. 'L'
C
C
C  MAXNRAD = NINT(4.0*(NXOUT*NYOUT)/numproc)
  MAXNRAD = NXOUT*NYOUT
C  ALLOCATE (RadSend(NSTOKES,MAXNRAD), IndexSend(MAXNRAD))
C  IF (myproc==0) ALLOCATE (NradRecv(numproc-1), &
C              RadRecv(NSTOKES,MAXNRAD,numproc-1), IndexRecv(MAXNRAD,numproc-1))
C  ALLOCATE (BndRadRealInfo(4+NSTOKES,MAXNRAD,4), BndRadIntInfo(MAXNRAD,4))
C  ALLOCATE (requests(MAX(12,3*numproc)), status(MPI_STATUS_SIZE,MAX(12,3*numproc)))
C
C
C   ! Loop over the radiance directions
C  DO K = 1, NANGOUT
C    MUOUT = OUTPARMS(2*K+5)
C    PHID = OUTPARMS(2*K+6)
C    PHIOUT = PHID*ACOS(-1.0)/180.0
C    IF (MUOUT .EQ. 0.0 .OR. ABS(MUOUT) .GT. 1.0) THEN
C      WRITE (6,*) 'COMPUTE_RADIANCE: Bad mu for radiance',MUOUT
C    ELSE
C
C       ! Compute the source function throughout grid for this angle
C      CALL COMPUTE_ONE_SOURCE (NSTOKES, ML, MM, NLM, NSTLEG, &
C                NLEG, NUMPHASE, NPTS, DELTAM, MUOUT, PHIOUT, &
C                SRCTYPE, SOLARMU, SOLARAZ, ALBEDO, LEGEN, IPHASE, &
C                DIRFLUX, SHPTR, SOURCE,  SOURCE1)
C
C       ! Set the radiance field to -1, so we can determine valid radiances
C       !   Also set the source array to the extinction times the source
C       !   function, so the grid interpolation is correct.
C      DO I = 1, NPTS
C        GRIDRAD(1,I) = -1.0
C        SOURCE1(:,I) = SOURCE1(:,I)*EXTINCT(I)
C      ENDDO
C       ! Get boundary radiances: either top or bottom
C       ! Isotropic top boundary or Lambertian bottom boundary can use
C       !   the previously computed boundary radiances in BCRAD,
C       !   otherwise, compute the radiance for this angle
C       !   by integrating over the stored downwelling radiances.
C      IF (MUOUT .LT. 0.0) THEN
C        DO IBC = 1, NTOPPTS
C          I = BCPTR(IBC,1)
C          GRIDRAD(:,I) = BCRAD(:,IBC)
C        ENDDO
C      ELSE
C        IF (.NOT. LAMBERTIAN) THEN
C          CALL VARIABLE_BRDF_SURFACE (NBOTPTS,1,NBOTPTS,BCPTR(1,2), &
C                  NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MUOUT, PHIOUT, &
C                  SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX, &
C                  SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES, BCRAD(1,1+NTOPPTS))
C        ENDIF
C        DO IBC = 1, NBOTPTS
C          I = BCPTR(IBC,2)
C          GRIDRAD(:,I) = BCRAD(:,NTOPPTS+IBC)
C        ENDDO
C      ENDIF
C
C       ! Set up the upstream/downstream boundaries for this direction
C      Lupx=0 ; Lupy=0 ; Ldownx=0 ; Ldowny=0
C      IF (BTEST(BCFLAG,2)) THEN
C        IF (COS(PHIOUT) .GT. 1.0E-5) THEN
C          Ldownx = 2 ; Lupx = 1
C        ELSE
C          Ldownx = 1 ; Lupx = 2
C        ENDIF
C      ENDIF
C      IF (BTEST(BCFLAG,3)) THEN
C        IF (SIN(PHIOUT) .GT. 1.0E-5) THEN
C          Ldowny = 4 ; Lupy = 3
C        ELSE
C          Ldowny = 3 ; Lupy = 4
C        ENDIF
C      ENDIF
C
C      NINVALID(:)=0
C      Ndone = 0
C
C       ! Loop over all the radiance starting locations and integrate
C       !   backwards for those in this subdomain.
C      Y0 = STARTY
C      DO JY = 1, NYOUT
C        X0 = STARTX
C        DO JX = 1, NXOUT
C          NRAD = NRAD + 1
C          IF (X0 .GE. XGRID(1)-EPS .AND. X0 .LE. XGRID(NX)+EPS .AND. &
C              Y0 .GE. YGRID(1)-EPS .AND. Y0 .LE. YGRID(NY)+EPS) THEN
C            TRANSMIT = 1.0D0 ; RADIANCE(:) = 0.0D0
C            CALL INTEGRATE_SOURCE (BCFLAG, IPFLAG, NSTOKES, &
C                             NX, NY, NZ, NPTS, NCELLS, &
C                             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                             XGRID, YGRID, ZGRID, GRIDPOS, &
C                             MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1, &
C                             X0, Y0, Z0, XE,YE,ZE, SIDE, &
C                             TRANSMIT, RADIANCE, VALIDRAD)
C
C             ! If we got a radiance then store it in RadSend otherwise save
C             !  the information needed to pass to other processors
C            IF (VALIDRAD) THEN
C              Ndone = Ndone + 1
C              IF (Ndone>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 1')
C              RadSend(:,Ndone) = RADIANCE(:)
C              IndexSend(Ndone) = NRAD
C            ELSE
C              IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
C                WRITE (6,*) 'Bad SIDE 1:',myproc,muout,phiout,x0,y0,z0,SIDE,xe,ye,ze
C                CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR')
C              ENDIF
C              NINVALID(SIDE) = NINVALID(SIDE) + 1
C              IF (NINVALID(SIDE)>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 1')
C              BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ XE, YE, ZE, TRANSMIT, RADIANCE(:) /)
C              BndRadIntInfo(NINVALID(SIDE),SIDE) = NRAD
C            ENDIF
C          ENDIF
C          X0 = X0 + OUTPARMS(2)
C        ENDDO
C        Y0 = Y0 + OUTPARMS(3)
C      ENDDO
C
C       ! Cycle over the passing of info between processors and doing more
C       !  partial radiance integrations until all rays have been finished
C      AllDone = .FALSE.
C      DO WHILE (.NOT. AllDone)
C         ! Do the send and receives of the number of incomplete rays and
C         !   the real and integer info needed about the radiance integrations
C        irq=0
C        IF (BTEST(BCFLAG,2)) THEN
C          CALL MPI_ISEND (NINVALID(Lupx), 1, MPI_INTEGER, &
C                          iprocneigh(Lupx), 6*K+1, comm2d, requests(irq+1), ierr)
C          CALL MPI_ISEND (BndRadRealInfo(:,:,Lupx), (4+NSTOKES)*NINVALID(Lupx), MPI_DOUBLE_PRECISION, &
C                          iprocneigh(Lupx), 6*K+2, comm2d, requests(irq+2), ierr)
C          CALL MPI_ISEND (BndRadIntInfo(:,Lupx), NINVALID(Lupx), MPI_INTEGER, &
C                          iprocneigh(Lupx), 6*K+3, comm2d, requests(irq+3), ierr)
C          irq=irq+3
C        ENDIF
C        IF (BTEST(BCFLAG,3)) THEN
C          CALL MPI_ISEND (NINVALID(Lupy), 1, MPI_INTEGER, &
C                          iprocneigh(Lupy), 6*K+4, comm2d, requests(irq+1), ierr)
C          CALL MPI_ISEND (BndRadRealInfo(:,:,Lupy), (4+NSTOKES)*NINVALID(Lupy), MPI_DOUBLE_PRECISION, &
C                          iprocneigh(Lupy), 6*K+5, comm2d, requests(irq+2), ierr)
C          CALL MPI_ISEND (BndRadIntInfo(:,Lupy), NINVALID(Lupy), MPI_INTEGER, &
C                          iprocneigh(Lupy), 6*K+6, comm2d, requests(irq+3), ierr)
C          irq=irq+3
C        ENDIF
C        IF (BTEST(BCFLAG,2)) THEN
C          CALL MPI_IRECV (NINVALID(Ldownx), 1, MPI_INTEGER, &
C                          iprocneigh(Ldownx), 6*K+1, comm2d, requests(irq+1), ierr)
C          CALL MPI_IRECV (BndRadRealInfo(:,:,Ldownx), (4+NSTOKES)*MAXNRAD, MPI_DOUBLE_PRECISION, &
C                          iprocneigh(Ldownx), 6*K+2, comm2d, requests(irq+2), ierr)
C          CALL MPI_IRECV (BndRadIntInfo(:,Ldownx), MAXNRAD, MPI_INTEGER, &
C                          iprocneigh(Ldownx), 6*K+3, comm2d, requests(irq+3), ierr)
C          irq=irq+3
C        ENDIF
C        IF (BTEST(BCFLAG,3)) THEN
C          CALL MPI_IRECV (NINVALID(Ldowny), 1, MPI_INTEGER, &
C                          iprocneigh(Ldowny), 6*K+4, comm2d, requests(irq+1), ierr)
C          CALL MPI_IRECV (BndRadRealInfo(:,:,Ldowny), (4+NSTOKES)*MAXNRAD, MPI_DOUBLE_PRECISION, &
C                          iprocneigh(Ldowny), 6*K+5, comm2d, requests(irq+2), ierr)
C          CALL MPI_IRECV (BndRadIntInfo(:,Ldowny), MAXNRAD, MPI_INTEGER, &
C                          iprocneigh(Ldowny), 6*K+6, comm2d, requests(irq+3), ierr)
C          irq=irq+3
C        ENDIF
C        CALL MPI_WAITALL (irq, requests, status, ierr)
C
C
C        NINVALID(Lupx)=0 ; NINVALID(Lupy)=0
C        IF (BTEST(BCFLAG,2)) THEN
C           ! Continue the backwards ray integrations for the X boundary
C          DO J = 1, NINVALID(Ldownx)
C            IF (Ldownx == 1) BndRadRealInfo(1,J,Ldownx) = XGRID(1)
C            IF (Ldownx == 2) BndRadRealInfo(1,J,Ldownx) = XGRID(NX)
C            X=BndRadRealInfo(1,J,Ldownx) ; Y=BndRadRealInfo(2,J,Ldownx)
C            Z=BndRadRealInfo(3,J,Ldownx) ; TRANSMIT=BndRadRealInfo(4,J,Ldownx)
C            RADIANCE(:)=BndRadRealInfo(5:4+NSTOKES,J,Ldownx)
C            CALL INTEGRATE_SOURCE (BCFLAG, IPFLAG, NSTOKES, &
C                             NX, NY, NZ, NPTS, NCELLS, &
C                             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                             XGRID, YGRID, ZGRID, GRIDPOS, &
C                             MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1, &
C                             X, Y, Z, XE,YE,ZE, SIDE, &
C                             TRANSMIT, RADIANCE, VALIDRAD)
C             ! If we got a radiance then store it in RadSend, otherwise
C             ! save information needed to pass to neighboring processors
C            IF (VALIDRAD) THEN
C              Ndone = Ndone + 1
C              IF (Ndone>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 2')
C              RadSend(:,Ndone) = RADIANCE(:)
C              IndexSend(Ndone) = BndRadIntInfo(J,Ldownx)
C            ELSE
C              IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
C                WRITE (6,*) 'Bad SIDE 2x',myproc,muout,phiout,x,y,z,SIDE,xe,ye,ze
C                CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR')
C              ENDIF
C              NINVALID(SIDE) = NINVALID(SIDE) + 1
C              IF (NINVALID(SIDE)>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 2')
C              BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ XE, YE, ZE, TRANSMIT, RADIANCE(:) /)
C              BndRadIntInfo(NINVALID(SIDE),SIDE) = BndRadIntInfo(J,Ldownx)
C            ENDIF
C          ENDDO
C        ENDIF
C
C        IF (BTEST(BCFLAG,3)) THEN
C           ! Continue the backwards ray integrations for the Y boundary
C          DO J = 1, NINVALID(Ldowny)
C            IF (Ldowny == 3) BndRadRealInfo(2,J,Ldowny) = YGRID(1)
C            IF (Ldowny == 4) BndRadRealInfo(2,J,Ldowny) = YGRID(NY)
C            X=BndRadRealInfo(1,J,Ldowny) ; Y=BndRadRealInfo(2,J,Ldowny)
C            Z=BndRadRealInfo(3,J,Ldowny) ; TRANSMIT=BndRadRealInfo(4,J,Ldowny)
C            RADIANCE(:)=BndRadRealInfo(5:4+NSTOKES,J,Ldowny)
C            CALL INTEGRATE_SOURCE (BCFLAG, IPFLAG, NSTOKES, &
C                             NX, NY, NZ, NPTS, NCELLS, &
C                             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                             XGRID, YGRID, ZGRID, GRIDPOS, &
C                             MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1, &
C                             X, Y, Z, XE,YE,ZE, SIDE, &
C                             TRANSMIT, RADIANCE, VALIDRAD)
C             ! If we got a radiance then store it in RadSend, otherwise
C             ! save information needed to pass to neighboring processors
C            IF (VALIDRAD) THEN
C              Ndone = Ndone + 1
C              IF (Ndone>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 3')
C              RadSend(:,Ndone) = RADIANCE(:)
C              IndexSend(Ndone) = BndRadIntInfo(J,Ldowny)
C            ELSE
C              IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
C                WRITE (6,*) 'Bad SIDE 2y',myproc,muout,phiout,x,y,z,SIDE,xe,ye,ze
C                CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR')
C              ENDIF
C              NINVALID(SIDE) = NINVALID(SIDE) + 1
C              IF (NINVALID(SIDE)>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 3')
C              BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ XE, YE, ZE, TRANSMIT, RADIANCE(:) /)
C              BndRadIntInfo(NINVALID(SIDE),SIDE) = BndRadIntInfo(J,Ldowny)
C            ENDIF
C          ENDDO
C        ENDIF
C
C         ! See if there are any more invalid radiance rays on all the processors
C        intsendbuf = (/ NINVALID(Lupx), NINVALID(Lupy) /)
C        CALL MPI_ALLREDUCE (intsendbuf, intrecvbuf, 2, MPI_INTEGER, MPI_MAX, &
C                            comm2d, ierr)
C        AllDone = ALL(intrecvbuf(1:2) == 0)
C      ENDDO
C
C      IF (myproc > 0) THEN
C         ! Send the finished radiances to the master processor
C        CALL MPI_ISEND (Ndone, 1, MPI_INTEGER, 0, 6*K+7, comm2d, requests(1), ierr)
C        CALL MPI_ISEND (RadSend, NSTOKES*Ndone, MPI_REAL, 0, 6*K+8, comm2d, requests(2), ierr)
C        CALL MPI_ISEND (IndexSend, Ndone, MPI_INTEGER, 0, 6*K+9, comm2d, requests(3), ierr)
C        CALL MPI_WAITALL (3, requests, status, ierr)
C      ELSE
C         ! Put the radiances finished by the master processor in the list
C        DO i = 1, Ndone
C          RADOUT(:,IndexSend(i)) = RadSend(:,i)
C        ENDDO
C         ! Receive the radiances from all the other processors
C        irq=0
C        DO iproc = 1, numproc-1
C          CALL MPI_IRECV (NradRecv(iproc), 1, MPI_INTEGER, iproc, &
C                          6*K+7, comm2d, requests(irq+1), ierr)
C          CALL MPI_IRECV (RadRecv(:,:,iproc), NSTOKES*MAXNRAD, MPI_REAL, iproc, &
C                          6*K+8, comm2d, requests(irq+2), ierr)
C          CALL MPI_IRECV (IndexRecv(:,iproc), MAXNRAD, MPI_INTEGER, iproc, &
C                          6*K+9, comm2d, requests(irq+3), ierr)
C          irq = irq + 3
C        ENDDO
C        CALL MPI_WAITALL (irq, requests, status, ierr)
C        DO iproc = 1, numproc-1
C          DO i = 1, NradRecv(iproc)
C            RADOUT(:,IndexRecv(i,iproc)) = RadRecv(:,i,iproc)
C          ENDDO
C        ENDDO
C      ENDIF
C
C    ENDIF
C  ENDDO
C
C
C  DEALLOCATE (requests, status, RadSend, IndexSend)
C  IF (myproc==0) DEALLOCATE (NradRecv, RadRecv, IndexRecv)
C  DEALLOCATE (BndRadRealInfo, BndRadIntInfo)
CEND SUBROUTINE COMPUTE_RADIANCE_PAR
C
C
C
C
C
CSUBROUTINE VISUALIZE_RADIANCE_PAR (NSTOKES, NX, NY, NZ, NPTS, NCELLS, &
C                  ML, MM, NLM, NSTLEG, NLEG, NUMPHASE, &
C                  NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
C                  BCFLAG, XDOMAIN, YDOMAIN, IPFLAG, &
C                  SRCTYPE, DELTAM, SOLARMU, SOLARAZ, &
C                  SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
C                  MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
C                  GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS, &
C                  XGRID, YGRID, ZGRID, GRIDPOS, &
C                  GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                  EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, &
C                  SHPTR, SOURCE,  OUTPARMS,  IVIS, VISOUT)
C ! Computes Stokes radiances (output in VISOUT) for the visualization
C ! modes: 1) camera mode, and 2) cross track scanning.
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(IN) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
C  INTEGER, INTENT(IN) :: ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
C  INTEGER, INTENT(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
C  INTEGER, INTENT(INOUT) :: IVIS
C  INTEGER, INTENT(IN) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
C  INTEGER, INTENT(IN) :: GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
C  INTEGER, INTENT(IN) :: SHPTR(NPTS+1), BCPTR(MAXNBC,2)
C  INTEGER*2, INTENT(IN) :: CELLFLAGS(NCELLS), IPHASE(NPTS)
C  LOGICAL, INTENT(IN) :: DELTAM
C  REAL, INTENT(IN) :: SOLARMU, SOLARAZ
C  REAL, INTENT(IN) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
C  REAL, INTENT(IN) :: MU(NMU), PHI(NMU,NPHI0MAX), WTDO(NMU,NPHI0MAX)
C  REAL, INTENT(IN) :: XDOMAIN, YDOMAIN, XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)
C  REAL, INTENT(IN) :: GRIDPOS(3,NPTS)
C  REAL, INTENT(IN) :: SFCGRIDPARMS(*), BCRAD(NSTOKES,*)
C  REAL, INTENT(IN) :: EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG,0:NLEG,*)
C  REAL, INTENT(IN) :: DIRFLUX(NPTS), FLUXES(2,NPTS), SOURCE(NSTOKES,*)
C  REAL, INTENT(IN) :: OUTPARMS(*)
C  REAL, INTENT(OUT) :: VISOUT(NSTOKES,*)
C  CHARACTER, INTENT(IN) :: SRCTYPE*1, SFCTYPE*2, UNITS*1
C
C  INTEGER :: NSCATANGLE, NSTPHASE, I, J, L, N
C  INTEGER :: NL, NS, LINE, SAMP
C  LOGICAL :: CAMERA_MODE, VALIDRAD
C  REAL, PARAMETER  :: EPS=1.0E-5
C  DOUBLE PRECISION :: XE,YE,ZE, TRANSMIT, VISRAD(NSTOKES), RADIANCE(NSTOKES)
C  DOUBLE PRECISION :: X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, X, Y, Z
C  DOUBLE PRECISION :: THETA0, THETA1, PHIR, PHI0, COSSCAT, MU2, PHI2
C  DOUBLE PRECISION :: COSTH, SINTH, SINTH0, SINTH1
C  DOUBLE PRECISION :: U, V, UP, VP, COSDELPHI, PI, ROTANG, DEGRAD, R
C  DOUBLE PRECISION :: DIST, D, RX, RY, RZ, SCANANG
C  INTEGER, PARAMETER :: MAXSCATANG=721
C  REAL, ALLOCATABLE :: YLMSUN(:,:), PHASETAB(:,:,:)
C  REAL              :: MEAN, STD1, STD2
C  REAL, ALLOCATABLE :: AOLP(:)
C
C   ! Domain SIDEs: 1=-X, 2=+X, 3=-Y, 4=+Y
C  INTEGER :: MAXNRAD
C  INTEGER, PARAMETER :: RECVSIDE(4) = (/2, 1, 4, 3 /)
C  INTEGER :: SIDE, IBND, NINVALID(4), NINCOMING(4)
C  INTEGER :: Ndone, intrecvbuf(4)
C  INTEGER :: ierr, irq, iproc, K
C  LOGICAL :: AllDone
C  INTEGER, ALLOCATABLE :: requests(:), status(:,:), NradRecv(:)
C  INTEGER, ALLOCATABLE :: IndexSend(:), IndexRecv(:,:)
C  INTEGER, ALLOCATABLE :: BndSendRadIntInfo(:,:), BndRecvRadIntInfo(:,:)
C  REAL, ALLOCATABLE    :: RadSend(:,:), RadRecv(:,:,:)
C  REAL, ALLOCATABLE    :: BndSendRadRealInfo(:,:,:), BndRecvRadRealInfo(:,:,:)
C
C
C  ALLOCATE (YLMSUN(NSTLEG,NLM))
C
C  IF (SRCTYPE .NE. 'T') THEN
C    CALL YLMALL (.TRUE., SOLARMU, SOLARAZ, ML, MM, NSTLEG, YLMSUN)
C    IF (DELTAM .AND. NUMPHASE .GT. 0) THEN
C      NSCATANGLE = MAX(36,MIN(MAXSCATANG,2*NLEG))
C      NSTPHASE = MIN(NSTLEG,2)
C      ALLOCATE (PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE))
C      CALL PRECOMPUTE_PHASE (NSCATANGLE, NUMPHASE, NSTPHASE, &
C                             NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)
C    ENDIF
C  ENDIF
C
C
C   ! Make the isotropic radiances for the top boundary
C  CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN, &
C                              UNITS, NTOPPTS, NSTOKES, BCRAD(1,1))
C   ! Make the bottom boundary radiances for the Lambertian surfaces.
C   ! Compute the upwelling bottom radiances using the downwelling fluxes.
C  IF (SFCTYPE .EQ. 'FL') THEN
C    CALL FIXED_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2), &
C                   DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO, &
C                   WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1,1+NTOPPTS))
C  ELSE IF (SFCTYPE .EQ. 'VL') THEN
C    CALL VARIABLE_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2), &
C                   DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS, &
C                   NSTOKES, BCRAD(1,1+NTOPPTS))
C  ENDIF
C
C  CAMERA_MODE = NINT(OUTPARMS(1)) .EQ. 1
C  IF (CAMERA_MODE) THEN
C    NL = NINT(OUTPARMS(10))
C    NS = NINT(OUTPARMS(11))
C  ELSE
C    NL = 1 + INT( SQRT((OUTPARMS(4)-OUTPARMS(7))**2 &
C                   +(OUTPARMS(5)-OUTPARMS(8))**2 &
C                   +(OUTPARMS(6)-OUTPARMS(9))**2) /OUTPARMS(10) )
C    NS = 1 + INT(ABS(OUTPARMS(12)-OUTPARMS(11))/OUTPARMS(13))
C  ENDIF
C  PI = ACOS(-1.0D0)
C  DEGRAD = PI/180.
C
C   ! Allocate sending, receiving, etc. arrays for MPI
C  MAXNRAD = NL*NS
C  ALLOCATE (RadSend(NSTOKES,MAXNRAD), IndexSend(MAXNRAD))
C  IF (myproc==0) ALLOCATE (NradRecv(numproc-1), &
C              RadRecv(NSTOKES,MAXNRAD,numproc-1), IndexRecv(MAXNRAD,numproc-1))
C  ALLOCATE (BndSendRadRealInfo(6+NSTOKES,MAXNRAD,4), BndSendRadIntInfo(MAXNRAD,4))
C  ALLOCATE (BndRecvRadRealInfo(6+NSTOKES,MAXNRAD,4), BndRecvRadIntInfo(MAXNRAD,4))
C  ALLOCATE (requests(MAX(24,3*numproc)), status(MPI_STATUS_SIZE,MAX(24,3*numproc)))
C
C  NINVALID(:)=0
C  Ndone = 0
C
C   ! Loop over pixels in image
C  DO LINE = 1, NL
C  DO SAMP = 1, NS
C    IVIS = IVIS + 1
C
C    IF (CAMERA_MODE) THEN
C       ! Camera mode:
C       !   1, bytes, scale, X,Y,Z, theta, phi, rotang, NL, NS, delline, delsamp
C       ! Use spherical trig to find the pixel direction (MURAY,PHIRAY)
C       ! from the camera center (THETA0,PHI0) and the relative pixel angles (U,V).
C      UP = (SAMP-NS/2-1)*OUTPARMS(13)*DEGRAD
C      VP = (LINE-NL/2-1)*OUTPARMS(12)*DEGRAD
C      ROTANG = OUTPARMS(9)*DEGRAD
C      IF (ROTANG .EQ. 0.0) THEN
C        U = UP
C        V = VP
C      ELSE
C        U = COS(ROTANG)*UP - SIN(ROTANG)*VP
C        V = SIN(ROTANG)*UP + COS(ROTANG)*VP
C      ENDIF
C      THETA0 = DEGRAD*OUTPARMS(7)
C      PHI0 = DEGRAD*OUTPARMS(8)
C      THETA1 = THETA0 + V
C      IF (V .EQ. 0.0) THEN
C        COSTH = COS(U)*COS(THETA0)
C      ELSE
C        COSTH = COS(U)*(SIN(THETA1)*COS(V)-SIN(THETA0))/SIN(V)
C        COSTH = MIN(+1.0D0,MAX(-1.0D0,COSTH))
C      ENDIF
C      SINTH = SQRT(1-COSTH**2)
C      SINTH0 = SIN(THETA0)
C      SINTH1 = SIN(THETA1)
C      IF (ABS(SINTH) .LT. 1.0E-6) THEN
C        PHIR = 0.0
C      ELSE
C        IF (ABS(SINTH0).LT.1E-6 .AND. ABS(SINTH1).LE.1E-6) THEN
C          COSDELPHI = 0.0D0
C        ELSE IF (ABS(SINTH1) .GT. 1.0E-6) THEN
C          COSDELPHI = (COS(U)-COSTH*COS(THETA1))/(SINTH1*SINTH)
C        ELSE IF (ABS(SINTH0) .GT. 1.0E-6) THEN
C          COSDELPHI = (COS(U)*COS(V)-COSTH*COS(THETA0))/(SINTH0*SINTH)
C        ENDIF
C        COSDELPHI = MIN(+1.0D0,MAX(-1.0D0,COSDELPHI))
C        IF (U .GE. 0.0) THEN
C          PHIR = PHI0 - ACOS(COSDELPHI)
C        ELSE
C          PHIR = PHI0 + ACOS(COSDELPHI)
C        ENDIF
C      ENDIF
C
C      X0 = OUTPARMS(4)
C      Y0 = OUTPARMS(5)
C      Z0 = OUTPARMS(6)
C
C    ELSE
C
C       ! Cross track scanning in the vertical plane:
C       !  2, bytes, scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1, scan2, delscan
C       ! Start and end of scan angle (scan1, scan2) are +/- relative
C       ! to nadir, and positive is on right side.
C      X1 = OUTPARMS(4)
C      Y1 = OUTPARMS(5)
C      Z1 = OUTPARMS(6)
C      X2 = OUTPARMS(7)
C      Y2 = OUTPARMS(8)
C      Z2 = OUTPARMS(9)
C      DIST = SQRT( (X2-X1)**2 + (Y2-Y1)**2 + (Z2-Z1)**2)
C      RX = (X2-X1)/DIST
C      RY = (Y2-Y1)/DIST
C      RZ = (Z2-Z1)/DIST
C      D = (NL-LINE)*OUTPARMS(10)
C      X0 = X1 + D*RX
C      Y0 = Y1 + D*RY
C      Z0 = Z1 + D*RZ
C      SCANANG = DEGRAD*(OUTPARMS(11) + (SAMP-1)*OUTPARMS(13))
C      COSTH = COS(PI-ABS(SCANANG))
C      IF (SCANANG .GT. 0.0) THEN
C        PHIR = ATAN2(-RX,RY)
C      ELSE
C        PHIR = ATAN2(RX,-RY)
C      ENDIF
C    ENDIF
C
C     ! Extrapolate ray to domain top if above
C    IF (Z0 .GT. ZGRID(NZ)) THEN
C      IF (COSTH .GE. 0.0) THEN
C        Ndone = Ndone + 1
C        RadSend(:,Ndone) = 0.0
C        IndexSend(Ndone) = IVIS
C        CYCLE
C      ENDIF
C      R = (ZGRID(NZ) - Z0)/COSTH
C      X0 = X0 + R*SQRT(1-COSTH**2)*COS(PHIR)
C      Y0 = Y0 + R*SQRT(1-COSTH**2)*SIN(PHIR)
C      Z0 = ZGRID(NZ)
C    ELSE IF (Z0 .LT. ZGRID(1)) THEN
C      WRITE (6,*) 'VISUALIZE_RADIANCE: Level below domain'
C      STOP
C    ENDIF
C
C     ! Relocate the starting point in the overall periodic domain
C     ! (we know we have a periodic domain if doing MPI)
C    IF (X0 >= 0.0) THEN
C      X0 = MOD(X0,DBLE(XDOMAIN))
C    ELSE
C      X0 = MOD(X0+XDOMAIN*INT(1.0-X0/XDOMAIN),DBLE(XDOMAIN))
C    ENDIF
C    IF (Y0 >= 0.0) THEN
C      Y0 = MOD(Y0,DBLE(YDOMAIN))
C    ELSE
C      Y0 = MOD(Y0+YDOMAIN*INT(1.0-Y0/YDOMAIN),DBLE(YDOMAIN))
C    ENDIF
C
C
C     ! Start the radiance calculation for those rays starting in this subdomain
C    IF (X0 .GE. XGRID(1) .AND. X0 .LE. XGRID(NX) .AND. &
C        Y0 .GE. YGRID(1) .AND. Y0 .LE. YGRID(NY)) THEN
C
C       ! MU2,PHI2 is radiance travel direction (opposite the camera viewing direction)
C      MU2 = -COSTH
C      PHI2 = PHIR + PI
C
C       ! Integrate the extinction and source function along this ray
C       ! to calculate the Stokes radiance vector for this pixel
C      TRANSMIT = 1.0D0 ; RADIANCE(:) = 0.0D0
C      CALL INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES, NSTLEG, &
C                         NSTPHASE, NSCATANGLE, PHASETAB, &
C                         NX, NY, NZ, NPTS, NCELLS, &
C                         GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                         XGRID, YGRID, ZGRID, GRIDPOS, &
C                         ML, MM, NLM, NLEG, NUMPHASE, &
C                         NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
C                         DELTAM, SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, &
C                         EXTINCT, ALBEDO, LEGEN, IPHASE, &
C                         DIRFLUX, SHPTR, SOURCE, YLMSUN, &
C                         MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
C                         SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
C                         MU2, PHI2, X0,Y0,Z0, &
C                         XE,YE,ZE, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
C
C             ! If we got a radiance then store it in RadSend otherwise save
C             !  the information needed to pass to other processors
C      IF (VALIDRAD) THEN
C        Ndone = Ndone + 1
C        IF (Ndone>MAXNRAD) CALL ABORT_SHDOM_MPI ('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 1')
C        RadSend(:,Ndone) = RADIANCE(:)
C        IndexSend(Ndone) = IVIS
C      ELSE
C        IF (SIDE < 1 .OR. SIDE > 4) CALL ABORT_SHDOM_MPI ('VISUALIZE_RADIANCE_PAR: invalid SIDE 1')
C        NINVALID(SIDE) = NINVALID(SIDE) + 1
C        IF (NINVALID(SIDE)>MAXNRAD) CALL ABORT_SHDOM_MPI ('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 1')
C        BndSendRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ XE, YE, ZE, MU2, PHI2, TRANSMIT, RADIANCE(:) /)
C        BndSendRadIntInfo(NINVALID(SIDE),SIDE) = IVIS
C      ENDIF
C    ENDIF
C  ENDDO
C  ENDDO
C
C   ! Cycle over the passing of info between processors and doing more
C   !  partial radiance integrations until all rays have been finished
C  AllDone = .FALSE.
C  DO WHILE (.NOT. AllDone)
C    irq=0
C     ! Do the sends for all four sides of the subdomain.  The number of
C     ! incomplete rays and the real and integer info needed about the
C     ! radiance integrations is sent.
C    DO SIDE = 1, 4
C      IF ( (BTEST(BCFLAG,2) .AND. (SIDE==1 .OR. SIDE==2)) .OR. &
C           (BTEST(BCFLAG,3) .AND. (SIDE==3 .OR. SIDE==4)) ) THEN
C        K = 3*SIDE
C        CALL MPI_ISEND (NINVALID(SIDE), 1, MPI_INTEGER, &
C                      iprocneigh(SIDE), K+1, comm2d, requests(irq+1), ierr)
C        CALL MPI_ISEND (BndSendRadRealInfo(:,:,SIDE), (6+NSTOKES)*NINVALID(SIDE), MPI_DOUBLE_PRECISION, &
C                      iprocneigh(SIDE), K+2, comm2d, requests(irq+2), ierr)
C        CALL MPI_ISEND (BndSendRadIntInfo(:,SIDE), NINVALID(SIDE), MPI_INTEGER, &
C                      iprocneigh(SIDE), K+3, comm2d, requests(irq+3), ierr)
C        irq=irq+3
C      ENDIF
C    ENDDO
C     ! Do the receives for all four sides of the subdomain.
C    DO SIDE = 1, 4
C      IF ( (BTEST(BCFLAG,2) .AND. (SIDE==1 .OR. SIDE==2)) .OR. &
C           (BTEST(BCFLAG,3) .AND. (SIDE==3 .OR. SIDE==4)) ) THEN
C        K = 3*RECVSIDE(SIDE)
C        CALL MPI_IRECV (NINCOMING(SIDE), 1, MPI_INTEGER, &
C                      iprocneigh(SIDE), K+1, comm2d, requests(irq+1), ierr)
C        CALL MPI_IRECV (BndRecvRadRealInfo(:,:,SIDE), (6+NSTOKES)*MAXNRAD, MPI_DOUBLE_PRECISION, &
C                      iprocneigh(SIDE), K+2, comm2d, requests(irq+2), ierr)
C        CALL MPI_IRECV (BndRecvRadIntInfo(:,SIDE), MAXNRAD, MPI_INTEGER, &
C                      iprocneigh(SIDE), K+3, comm2d, requests(irq+3), ierr)
C        irq=irq+3
C      ENDIF
C    ENDDO
C    CALL MPI_WAITALL (irq, requests, status, ierr)
C
C     ! Continue the backwards ray integrations for the four boundaries
C    NINVALID(:)=0
C    DO IBND = 1, 4
C      DO J = 1, NINCOMING(IBND)
C        IF (IBND == 1) BndRecvRadRealInfo(1,J,IBND) = XGRID(1)
C        IF (IBND == 2) BndRecvRadRealInfo(1,J,IBND) = XGRID(NX)
C        IF (IBND == 3) BndRecvRadRealInfo(2,J,IBND) = YGRID(1)
C        IF (IBND == 4) BndRecvRadRealInfo(2,J,IBND) = YGRID(NY)
C        X=BndRecvRadRealInfo(1,J,IBND)
C        Y=BndRecvRadRealInfo(2,J,IBND)
C        Z=BndRecvRadRealInfo(3,J,IBND)
C        MU2 = BndRecvRadRealInfo(4,J,IBND)
C        PHI2 = BndRecvRadRealInfo(5,J,IBND)
C        TRANSMIT=BndRecvRadRealInfo(6,J,IBND)
C        RADIANCE(:)=BndRecvRadRealInfo(7:6+NSTOKES,J,IBND)
C        CALL INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES, NSTLEG, &
C                         NSTPHASE, NSCATANGLE, PHASETAB, &
C                         NX, NY, NZ, NPTS, NCELLS, &
C                         GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
C                         XGRID, YGRID, ZGRID, GRIDPOS, &
C                         ML, MM, NLM, NLEG, NUMPHASE, &
C                         NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
C                         DELTAM, SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, &
C                         EXTINCT, ALBEDO, LEGEN, IPHASE, &
C                         DIRFLUX, SHPTR, SOURCE, YLMSUN, &
C                         MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
C                         SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
C                         MU2, PHI2, X,Y,Z, &
C                         XE,YE,ZE, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
C
C           ! If we got a valid radiance then store it in RadSend, otherwise
C           ! save information needed to pass to neighboring processors
C        IF (VALIDRAD) THEN
C          Ndone = Ndone + 1
C          IF (Ndone>MAXNRAD) CALL ABORT_SHDOM_MPI ('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 2')
C          RadSend(:,Ndone) = RADIANCE(:)
C          IndexSend(Ndone) = BndRecvRadIntInfo(J,IBND)
C        ELSE
C          IF (SIDE < 1 .OR. SIDE > 4) CALL ABORT_SHDOM_MPI ('VISUALIZE_RADIANCE_PAR: invalid SIDE 2')
C          NINVALID(SIDE) = NINVALID(SIDE) + 1
C          IF (NINVALID(SIDE)>MAXNRAD) CALL ABORT_SHDOM_MPI ('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 2')
C          BndSendRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ XE, YE, ZE, MU2, PHI2, TRANSMIT, RADIANCE(:) /)
C          BndSendRadIntInfo(NINVALID(SIDE),SIDE) = BndRecvRadIntInfo(J,IBND)
C        ENDIF
C      ENDDO
C    ENDDO
C
C     ! See if there are any more invalid radiance rays on all the processors
C    CALL MPI_ALLREDUCE (NINVALID, intrecvbuf, 4, MPI_INTEGER, MPI_MAX, &
C                        comm2d, ierr)
C    AllDone = ALL(intrecvbuf(1:4) == 0)
C  ENDDO
  print *,'Done with ray tracing.', Ndone
C
C  IF (myproc > 0) THEN
C     ! Send the finished radiances to the master processor
C    CALL MPI_ISEND (Ndone, 1, MPI_INTEGER, 0, 21, comm2d, requests(1), ierr)
C    CALL MPI_ISEND (RadSend, NSTOKES*Ndone, MPI_REAL, 0, 22, comm2d, requests(2), ierr)
C    CALL MPI_ISEND (IndexSend, Ndone, MPI_INTEGER, 0, 23, comm2d, requests(3), ierr)
C    CALL MPI_WAITALL (3, requests, status, ierr)
C  ELSE
C     ! Put the radiances finished by the master processor in the list
C    DO i = 1, Ndone
C      VISOUT(:,IndexSend(i)) = RadSend(:,i)
C    ENDDO
C     ! Receive the radiances from all the other processors
C    irq=0
C    DO iproc = 1, numproc-1
C      CALL MPI_IRECV (NradRecv(iproc), 1, MPI_INTEGER, iproc, &
C                      21, comm2d, requests(irq+1), ierr)
C      CALL MPI_IRECV (RadRecv(:,:,iproc), NSTOKES*MAXNRAD, MPI_REAL, iproc, &
C                      22, comm2d, requests(irq+2), ierr)
C      CALL MPI_IRECV (IndexRecv(:,iproc), MAXNRAD, MPI_INTEGER, iproc, &
C                      23, comm2d, requests(irq+3), ierr)
C      irq = irq + 3
C    ENDDO
C    CALL MPI_WAITALL (irq, requests, status, ierr)
C    DO iproc = 1, numproc-1
C      DO i = 1, NradRecv(iproc)
C        VISOUT(:,IndexRecv(i,iproc)) = RadRecv(:,i,iproc)
C      ENDDO
C    ENDDO
C  ENDIF
C
C  DEALLOCATE (requests, status, RadSend, IndexSend)
C  IF (myproc==0) DEALLOCATE (NradRecv, RadRecv, IndexRecv)
C  DEALLOCATE (BndSendRadRealInfo, BndSendRadIntInfo)
C  DEALLOCATE (BndRecvRadRealInfo, BndRecvRadIntInfo)
C  DEALLOCATE (YLMSUN)
C  IF (ALLOCATED(PHASETAB))  DEALLOCATE (PHASETAB)
C
C   ! On the master processor convert VISOUT(2:,*) from Stokes parameters
C   ! to DOLP, AOLP, DOCP
C  IF (myproc==0) THEN
C    IF (NSTOKES .GT. 1) THEN
C      N = NL*NS
C      DO I = IVIS-N+1, IVIS
C        VISRAD(:) = VISOUT(:,I)
C        IF (VISRAD(1) > 0.0) THEN
C         IF (NSTOKES .GT. 1) THEN
C           ! Output degree (0 to 1) and angle (-180 to 180) of linear polarization
C          VISOUT(2,I) = SQRT(VISRAD(2)**2+VISRAD(3)**2)/VISRAD(1)
C          VISOUT(3,I) = (180/PI)*0.5*ATAN2(VISRAD(3),VISRAD(2))
C         ENDIF
C         IF (NSTOKES .EQ. 4) THEN
C            ! Output degree of circular polarization (-1 to 1)
C           VISOUT(4,I) = VISRAD(4)/VISRAD(1)
C         ENDIF
C        ENDIF
C      ENDDO
C       ! Choose the best range for the angle of linear polarization (-90 to 90 or 0 to 180)
C      ALLOCATE (AOLP(N))
C      AOLP(:) = VISOUT(3,IVIS-N+1:IVIS)
C      MEAN = SUM(AOLP(:))/N
C      STD1 = SQRT(SUM((AOLP(:)-MEAN)**2)/N)
C      WHERE (AOLP(:) < 0.0)
C        AOLP(:) = AOLP(:)+180.0
C      END WHERE
C      MEAN = SUM(AOLP(:))/N
C      STD2 = SQRT(SUM((AOLP(:)-MEAN)**2)/N)
C      IF (STD2 < STD1) THEN
C        VISOUT(3,IVIS-N+1:IVIS) = AOLP(:)
C      ENDIF
C      DEALLOCATE (AOLP)
C    ENDIF
C  ENDIF
CEND SUBROUTINE VISUALIZE_RADIANCE_PAR
C
C
C
C
C
CSUBROUTINE CALC_ACCEL_SOLCRIT (DOACCEL, DELJDOT, DELJOLD, DELJNEW, JNORM, &
C                               ACCELPAR, SOLCRIT)
C ! Does an MPI_ALLREDUCE to sum the delta source function vector dot
C ! products over all the processors, and from these calculates the
C ! full-domain-wide acceleration parameter and solution criterion.
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  LOGICAL, INTENT(IN) :: DOACCEL
C  REAL,    INTENT(INOUT) :: DELJDOT, DELJOLD, DELJNEW, JNORM
C  REAL,    INTENT(OUT) :: ACCELPAR, SOLCRIT
C  INTEGER :: ierr
C  REAL    ::  R, THETA, sendbuf(4), recvbuf(4)
C  REAL, SAVE :: A=0.0
C
C  IF (numproc > 1) THEN
C     ! Call MPI_ALLREDUCE to sum up the four dot products
C    sendbuf(:) = (/ DELJDOT, DELJOLD, DELJNEW, JNORM /)
C    call MPI_ALLREDUCE (sendbuf, recvbuf, 4, MPI_REAL, MPI_SUM, comm2d, ierr)
C    if (ierr .ne. MPI_SUCCESS) CALL ABORT_SHDOM_MPI ('CALC_ACCEL_SOLCRIT MPI_ALLREDUCE error')
C    DELJDOT=recvbuf(1) ; DELJOLD=recvbuf(2)
C    DELJNEW=recvbuf(3) ; JNORM=recvbuf(4)
C  ENDIF
C
C   ! Accelerate if desired, didn't last time, and things are converging.
C  IF (DOACCEL .AND. A .EQ. 0.0 .AND. DELJNEW .LT. DELJOLD) THEN
C    ! Compute the acceleration extrapolation factor and apply it.
C    R = SQRT(DELJNEW/DELJOLD)
C    THETA = ACOS(DELJDOT/SQRT(DELJOLD*DELJNEW))
C    A = (1 - R*COS(THETA) + R**(1+0.5*3.14159/THETA)) &
C              /(1 + R**2  - 2*R*COS(THETA))  - 1.0
C    A = MIN(10.0,MAX(0.0,A))
C    ! WRITE (6,'(1X,A,3(1X,F7.3))') '! Acceleration: ', A,R,THETA
C  ELSE
C    A = 0.0
C  ENDIF
C  ACCELPAR = A
C
C  IF (JNORM .GT. 0.0) THEN
C    SOLCRIT = SQRT(DELJNEW/JNORM)
C  ELSE IF (DELJNEW .EQ. 0.0) THEN
C    SOLCRIT = 0.0
C  ENDIF
CEND SUBROUTINE CALC_ACCEL_SOLCRIT
C
C
C
C
CSUBROUTINE END_SHDOM_MPI (NPTS, GRIDPOS, NPX, NPY, XSTART, YSTART, DELX, DELY,&
C                          NPXT, NPYT, PROPFILE, RUNNAME)
C ! Shuts down MPI at the end, but mainly finds the number of gridpoints
C ! around each property grid point column, sends that array to the master
C ! process, and the master process write out the RUNNAME//"_load_balance.out"
C ! file.
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  INTEGER, INTENT(IN) :: NPTS, NPX, NPY, NPXT, NPYT
C  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS), XSTART, YSTART, DELX, DELY
C  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE, RUNNAME
C  INTEGER :: I, IX, IY, MAXN, SX, SY, EX, EY, NX, NY
C  INTEGER :: iproc, ierr, irq
C  INTEGER, ALLOCATABLE :: NUMGRID(:,:), NUMGRIDALL(:,:), NORM(:,:), IXYALL(:,:,:)
C  INTEGER, ALLOCATABLE :: recvbuf(:,:), buf(:,:), requests(:), status(:,:)
C  REAL :: ElapsedTime
C
C  IF (LoadBalance) THEN
C     ! Calculate the number of grid points around each property grid column
C    ALLOCATE (NUMGRID(NPX,NPY))
C    NUMGRID(:,:) = 0
C    DO I = 1, NPTS
C      IX = NINT((GRIDPOS(1,I)-XSTART)/DELX)+1
C      IY = NINT((GRIDPOS(2,I)-YSTART)/DELY)+1
C      IF (IX <= NPX .AND. IY <= NPY) THEN
C        NUMGRID(IX,IY) = NUMGRID(IX,IY) + 1
C      ENDIF
C    ENDDO
C
C
C    IF (numproc > 1) THEN
C       ! Have the master process collect the property grid locations of all
C      ALLOCATE (IXYALL(2,2,0:numproc-1))
C      CALL MPI_GATHER (IXYPRP,4,MPI_INTEGER, IXYALL,4,MPI_INTEGER, 0, comm2d, ierr)
C
C       ! Send the NUMGRIDs to the master processor
C      ALLOCATE (requests(numproc), status(MPI_STATUS_SIZE,numproc))
C      irq=0
C      IF (myproc > 0) THEN
C        irq = irq + 1
C        CALL MPI_ISEND (NUMGRID, NPX*NPY, MPI_INTEGER, 0, 777, comm2d, &
C                        requests(irq), ierr)
C      ELSE
C        MAXN = (MAXVAL(IXYALL(2,1,:)-IXYALL(1,1,:))+1)*(MAXVAL(IXYALL(2,2,:)-IXYALL(1,2,:))+1)
C        ALLOCATE (recvbuf(MAXN,1:numproc-1))
C        DO iproc = 1, numproc-1
C          irq = irq + 1
C          CALL MPI_IRECV (recvbuf(:,iproc), MAXN, MPI_INTEGER, iproc, 777, &
C                          comm2d, requests(irq), ierr)
C        ENDDO
C      ENDIF
C      CALL MPI_WAITALL (irq, requests, status, ierr)
C
C       ! Put the NUMGRIDS into the full domain array on the master processor.
C       !   Use a normalization counter (NORM) to deal with grid points that
C       !   get contribution from two processors.
C      IF (myproc == 0) THEN
C        ALLOCATE (NUMGRIDALL(NPXT,NPYT), NORM(NPXT,NPYT))
C        NUMGRIDALL(:,:) = 0  ;  NORM(:,:) = 0
C        SX=IXYALL(1,1,0) ; EX=IXYALL(2,1,0) ; NX=EX-SX+1
C        SY=IXYALL(1,2,0) ; EY=IXYALL(2,2,0) ; NY=EY-SY+1
C        NUMGRIDALL(SX:EX,SY:EY) = NUMGRID(:,:)
C        NORM(SX:EX,SY:EY) = RESHAPE( (/ 1 /), (/ NX,NY /), (/ 1 /))
C        DO iproc = 1, numproc-1
C          SX=IXYALL(1,1,iproc) ; EX=IXYALL(2,1,iproc) ; NX=EX-SX+1
C          SY=IXYALL(1,2,iproc) ; EY=IXYALL(2,2,iproc) ; NY=EY-SY+1
C          ALLOCATE (buf(NX,NY))
C          buf(:,:) = RESHAPE(recvbuf(1:NX*NY,iproc),(/ NX,NY /) )
C          NUMGRIDALL(SX:EX-1,SY:EY-1) = NUMGRIDALL(SX:EX-1,SY:EY-1) + buf(1:NX-1,1:NY-1)
C          NORM(SX:EX-1,SY:EY-1) = NORM(SX:EX-1,SY:EY-1) + RESHAPE( (/ 1 /), (/ NX-1,NY-1 /), (/ 1 /))
C          IX=MOD(EX-1,NPXT)+1 ; IY=MOD(EY-1,NPYT)+1
C          NUMGRIDALL(IX,SY:EY-1) = NUMGRIDALL(IX,SY:EY-1) + buf(NX,1:NY-1)
C          NORM(IX,SY:EY-1) = NORM(IX,SY:EY-1) + RESHAPE( (/ 1 /), (/ NY-1 /), (/ 1 /))
C          NUMGRIDALL(SX:EX-1,IY) = NUMGRIDALL(SX:EX-1,IY) + buf(1:NX-1,NY)
C          NORM(SX:EX-1,IY) = NORM(SX:EX-1,IY) + RESHAPE( (/ 1 /), (/ NX-1 /), (/ 1 /))
C          NUMGRIDALL(IX,IY) = NUMGRIDALL(IX,IY) + buf(NX,NY)
C          NORM(IX,IY) = NORM(IX,IY) + 1
C          DEALLOCATE (buf)
C        ENDDO
C        DEALLOCATE (recvbuf)
C        NUMGRIDALL(:,:) = NUMGRIDALL(:,:)/NORM(:,:)
C        DEALLOCATE (NORM)
C      ENDIF
C      DEALLOCATE (requests, status, IXYALL)
C
C    ELSE
C      ALLOCATE (NUMGRIDALL(NPXT,NPYT))
C      NUMGRIDALL(1:NPX,1:NPY) = NUMGRID(:,:)
C    ENDIF
C    DEALLOCATE (NUMGRID)
C
C    IF (myproc == 0) THEN
C        ! Write out the array to the load balancing file
C      OPEN (UNIT=10, FILE=TRIM(RUNNAME)//'_load_balance.out', STATUS='UNKNOWN')
C      WRITE (10,'(A)') '! SHDOM grid point load balancing file for property file:'
C      WRITE (10,'(A)') PROPFILE
C      WRITE (10,'(2(1X,I4),A)') NPXT, NPYT, ' ! property grid size (Nx,Ny)'
C      DO IX = 1, NPXT
C        DO IY = 1, NPYT
C          WRITE (10,'(I4,1X,I4,1X,I5)') IX, IY, NUMGRIDALL(IX,IY)
C        ENDDO
C      ENDDO
C      CLOSE (10)
C      DEALLOCATE (NUMGRIDALL)
C    ENDIF
C  ENDIF
C
C  ElapsedTime = MPI_WTIME() - StartTime
C  WRITE (6,'(A,F9.2)') 'Total elapsed time (sec): ', ElapsedTime
C  WRITE (6,*) 'Finished with SHDOM run: ',TRIM(RUNNAME)
C
C  ! Shut down MPI
C  IF (ALLOCATED(BXYPTR)) DEALLOCATE (BXYPTR,IBXYPTR,BXYPOS)
C  CALL MPI_FINALIZE (ierr)
CEND SUBROUTINE END_SHDOM_MPI
C
C
CSUBROUTINE ABORT_SHDOM_MPI (ERRSTR)
C  USE SHDOM_MPI_DATA
C  USE MPI
C  IMPLICIT NONE
C  CHARACTER(LEN=*), INTENT(IN) :: ERRSTR
C  INTEGER :: ierr, errorcode=7
C
C  WRITE (6,*) ERRSTR
C  CALL MPI_ABORT (comm2d,errorcode,ierr)
CEND SUBROUTINE ABORT_SHDOM_MPI
