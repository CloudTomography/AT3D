 ! All of the MPI routines for the parallel version of SHDOM are
 ! in this file.  If you do not have MPI, compile with shdom_nompi.f90
 ! instead.

MODULE SHDOM_MPI_DATA
  INTEGER :: numproc, myproc, comm2d, iprocneigh(4), IXYPRP(2,2), IXYOUT(2,2)
  LOGICAL, PARAMETER :: LoadBalance=.FALSE.
  DOUBLE PRECISION :: StartTime
  INTEGER, ALLOCATABLE :: BXYPTR(:,:,:), IBXYPTR(:,:,:)
  REAL,    ALLOCATABLE :: BXYPOS(:,:,:)
END MODULE



SUBROUTINE START_MPI (MASTERPROC)
 ! Initializes the MPI system, and gets the number of processors in use
 ! and the current processor.  If this is processor 0 then MASTERPROC
 ! is returned true.  Also starts the logging files "shdom???.log" 
 ! for the standard output from the non-master processors.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  LOGICAL, INTENT(OUT) :: MASTERPROC
!f2py intent(out) MASTERPROC

  INTEGER :: ierr
  include 'mpif.h'

  call MPI_INIT(ierr)
  IF (ierr .ne. MPI_SUCCESS) THEN
    WRITE (6,*) 'Error starting MPI version of SHDOM. Terminating.'
    STOP
  ENDIF

  call MPI_COMM_SIZE (MPI_COMM_WORLD, numproc, ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, myproc, ierr)

  call MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
!  call MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL, ierr)

  StartTime = MPI_WTIME()

  IF (myproc == 0) THEN
    MASTERPROC = .TRUE.
  ELSE
    MASTERPROC = .FALSE.
  ENDIF
END SUBROUTINE START_MPI



SUBROUTINE MAP_SHDOM_MPI (BCFLAG, NPX, NPY, NX, NY, NZ, DELX, DELY, PROPFILE, &
                          XSTART, YSTART, RUNNAME)
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
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: BCFLAG, NPX, NPY, NX, NY, NZ
!f2py intent(in, out) :: BCFLAG, NPX, NPY, NX, NY, NZ
  REAL,    INTENT(IN) :: DELX, DELY
  CHARACTER(LEN=64), INTENT(IN) :: PROPFILE, RUNNAME
  REAL,    INTENT(OUT) :: XSTART, YSTART
  INTEGER :: NPXT, NPYT, NXT, NYT
  INTEGER :: i, ierr, dims(2), coords(2)
  LOGICAL :: periodic(2) = (/ .TRUE., .TRUE. /), reorder=.FALSE.
  REAL    :: XEND, YEND
  CHARACTER(LEN=3) :: procstr
  CHARACTER(LEN=72) :: logfile
  include 'mpif.h'

  IF (numproc == 1) THEN
    XSTART = 0.0 ; YSTART = 0.0
    RETURN
  ENDIF

  IF (myproc > 0) THEN
    WRITE (procstr,'(I3.3)') myproc
    logfile = TRIM(RUNNAME)//procstr//'.log'
    OPEN (UNIT=6,FILE=logfile,STATUS='UNKNOWN')
  ENDIF

  NPXT=NPX ; NPYT=NPY ; NXT=NX ; NYT=NY

  ! Partition 2D Cartesian topology:
  dims(:) = 0
  call MPI_DIMS_CREATE (numproc, 2, dims, ierr)

   ! Create 2D Cartesian mapping of processors:
  call MPI_CART_CREATE (MPI_COMM_WORLD, 2, dims, periodic, reorder, &
                        comm2d, ierr)

  IF (dims(1) > 1) BCFLAG=IBSET(BCFLAG,2)
  IF (dims(2) > 1) BCFLAG=IBSET(BCFLAG,3)
   ! Multi-processor BC takes precedence over open BCs
  IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3))  BCFLAG=IAND(BCFLAG,12)

   ! Find where in the 2D grid of processors this is is at:
  call MPI_CART_COORDS (comm2d, myproc, 2, coords, ierr)

   ! Make the optimal (or default) property grid divisions (in IXYPRP)
  CALL OPTIMIZE_PROCESSOR_DOMAINS (dims, coords, NPXT, NPYT, PROPFILE, RUNNAME)

   ! Make the new property and base grid sizes and the IXYOUT arrays
  IF (dims(1) > 1) THEN
    NPX = IXYPRP(2,1)-IXYPRP(1,1)+1
    NX = NINT(NXT*FLOAT(NPX-1)/NPXT)+1
    IF (ABS(NPXT/FLOAT(NXT) - (NPX-1)/FLOAT(NX-1)) > 0.0001) THEN
      CALL ABORT_SHDOM_MPI ('MAP_SHDOM_MPI: Even X base grid cannot be preserved for this NX and number of processors')
    ENDIF
    IXYOUT(1,1) = NINT(NXT*FLOAT(IXYPRP(1,1)-1)/NPXT)+1
    IXYOUT(2,1) = NINT(NXT*FLOAT(NPX-1)/NPXT)
    XSTART = (IXYPRP(1,1)-1)*DELX
    XEND = (IXYPRP(2,1)-1)*DELX
  ELSE
    IXYOUT(1,1) = 1
    IXYOUT(2,1) = NXT
    XSTART = 0.0
    XEND = NPXT*DELX
  ENDIF

  IF (dims(2) > 1) THEN
    NPY = IXYPRP(2,2)-IXYPRP(1,2)+1
    NY = NINT(NYT*FLOAT(NPY-1)/NPYT)+1
    IF (ABS(NPYT/FLOAT(NYT) - (NPY-1)/FLOAT(NY-1)) > 0.0001) THEN
      CALL ABORT_SHDOM_MPI ('MAP_SHDOM_MPI: Even Y base grid cannot be preserved for this NY and number of processors')
    ENDIF
    IXYOUT(1,2) = NINT(NYT*FLOAT(IXYPRP(1,2)-1)/NPYT)+1
    IXYOUT(2,2) = NINT(NYT*FLOAT(NPY-1)/NPYT)
    YSTART = (IXYPRP(1,2)-1)*DELY
    YEND = (IXYPRP(2,2)-1)*DELY
  ELSE
    IXYOUT(1,2) = 1
    IXYOUT(2,2) = NYT
    YSTART = 0.0  
    YEND = NPYT*DELY
  ENDIF

   ! Get the processors (rank numbers) of the neighbors
  call MPI_CART_SHIFT (comm2d, 0, 1, iprocneigh(1), iprocneigh(2), ierr)
  call MPI_CART_SHIFT (comm2d, 1, 1, iprocneigh(3), iprocneigh(4), ierr)

  WRITE (6,'(A,2(1X,I2),A,I3)') 'Number of processors in X and Y:', dims(1:2),&
                             '   This processor: ',myproc
  WRITE (6,'(A,4(1X,I3))') 'Neighbors:',iprocneigh(1:4)
  WRITE (6,'(2(A,2F7.3))')  'Processor X range: ',XSTART,XEND, '   Y range: ',YSTART,YEND
  WRITE (6,*)
END SUBROUTINE MAP_SHDOM_MPI



SUBROUTINE OPTIMIZE_PROCESSOR_DOMAINS (ProcDims, ProcCoord, NPXT, NPYT, &
                                       PROPFILE, RUNNAME)
 ! Calculates the location of this processor (subdomain) in the property
 ! grid.  The starting and ending property grid point indices for X and Y
 ! are returned in IXYPRP.  Each processor sends its location in the 
 ! 2D array of subdomains to the master process and the master process 
 ! reads the RUNNAME//"_load_balance.inp" file (if there is one) and does 
 ! the optimization.  If the "shdom_proc_times.inp" file of relative
 ! processor times exists then the optimization uses this information on
 ! heterogeneous processors to try to give more work to faster processors.
 ! 
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ProcDims(2), ProcCoord(2), NPXT, NPYT
  CHARACTER(LEN=64), INTENT(IN) :: PROPFILE, RUNNAME
  INTEGER :: NPXTI, NPYTI
  INTEGER :: IX, IY, JX, JY, Nxd, Nyd, N, I, ierr
  INTEGER :: BndLines(ProcDims(1)+ProcDims(2)-2)
  INTEGER :: IDX(ProcDims(1)+1), IDY(ProcDims(2)+1)
  LOGICAL :: DoLoadBalancing
  INTEGER, PARAMETER :: NCYCLE=30    ! number of annealing cooling cycles
  INTEGER, PARAMETER :: NiterFac=100 ! number of iterations per length of optimization vector (BndLines)
  REAL,    PARAMETER :: Temp0=0.1    ! initial temperature (units of change in objective function "energy")
  REAL,    PARAMETER :: Tempfac=4.0/NCYCLE   ! exponential cooling rate
  INTEGER :: CYCLE, ITER, NITER
  INTEGER :: TrialLines(ProcDims(1)+ProcDims(2)-2), CurLines(ProcDims(1)+ProcDims(2)-2)
  INTEGER :: iseed=-3
  REAL    :: Temp, DelLines, u, ran
  REAL    :: LOAD_OBJ_FUNC, Energy, NewEnergy, LowestEnergy, DeltaE, Boltz
  INTEGER, ALLOCATABLE :: ProcCoordAll(:,:), Ngrid(:,:)
  REAL,    ALLOCATABLE :: wtproc(:)
  CHARACTER(LEN=64) :: PROPFILEI
  include 'mpif.h'

   ! Initialize the subdomain boundary lines to evenly spaced
  Nxd=ProcDims(1) ; Nyd=ProcDims(2) ; N=Nxd+Nyd-2
  BndLines(1:Nxd-1) = (/ (NINT(I*FLOAT(NPXT)/Nxd)+1, I=1,Nxd-1) /)
  BndLines(Nxd:N)   = (/ (NINT(I*FLOAT(NPYT)/Nyd)+1, I=1,Nyd-1) /)


  IF (LoadBalance) THEN
    DoLoadBalancing = .false.
    IF (myproc == 0) THEN
       ! See if we can read the load balancing file and header matches
      OPEN (UNIT=10, FILE=TRIM(RUNNAME)//'_load_balance.inp', STATUS='OLD', IOSTAT=ierr)
      IF (ierr == 0) THEN
        READ (10,*) 
        READ (10,'(A)') PROPFILEI
        READ (10,*) NPXTI, NPYTI
        IF (NPXTI /= NPXT .OR. NPYTI /= NPYT .OR. PROPFILEI /= PROPFILE) THEN
          WRITE (6,'(A,A,A)') 'Load balancing file ',TRIM(RUNNAME),&
               '_load_balance.inp is not compatible with property file:'
          WRITE (6,*) PROPFILE
          WRITE (6,*)
        ELSE
          DoLoadBalancing = .true.
        ENDIF
        CLOSE (10)
      ENDIF
    ENDIF


     ! Have the master process collect the locations of all 
    ALLOCATE (ProcCoordAll(2,0:numproc-1))
    CALL MPI_GATHER (ProcCoord,2,MPI_INTEGER, ProcCoordAll,2,MPI_INTEGER, 0, comm2d, ierr)

    IF (myproc==0 .AND. DoLoadBalancing) THEN
       ! Read the load balancing file to get array of number of grid points
      ALLOCATE (Ngrid(NPXT+1,NPYT+1))
      OPEN (UNIT=10, FILE=TRIM(RUNNAME)//'_load_balance.inp', STATUS='OLD')
      READ (10,*) 
      READ (10,*)
      READ (10,*) 
      DO IX = 1, NPXT
        DO IY = 1, NPYT
          READ (10,*) JX, JY, Ngrid(IX,IY)
        ENDDO
      ENDDO
      CLOSE (10)
      Ngrid(NPXT+1,:) = Ngrid(1,:)
      Ngrid(:,NPYT+1) = Ngrid(:,1)
      Ngrid(NPXT+1,NPYT+1) = Ngrid(1,1)

       ! Read the relative processor times if the "shdom_proc_times.inp" file exists
      ALLOCATE (wtproc(numproc))
      OPEN (UNIT=11, FILE='shdom_proc_times.inp', STATUS='OLD', IOSTAT=ierr)
      IF (ierr == 0) THEN
        DO i = 1, numproc
          READ (11,*) wtproc(i)
        ENDDO
        CLOSE (11)
        wtproc(:) = wtproc(:)*numproc/SUM(wtproc(:))
      ELSE
        wtproc(:) = 1.0
      ENDIF

       ! Do the simulated annealing to minimize the objective function
      u = ran(iseed)
      NITER = NiterFac*N
      Energy = LOAD_OBJ_FUNC (BndLines, Nxd, Nyd, ProcCoordAll, wtproc, NPXT, NPYT, Ngrid)
      LowestEnergy = Energy
      CurLines(:) = BndLines(:)
      ! print *, 'Annealing cycles: Cycle Temp Energy  BndLines'
      ! print '(I3,2(1X,F9.6),10(1X,I3))', 0,Temp0,Energy, BndLines(:)

       ! Loop through the cooling cycles
      DO CYCLE = 1, NCYCLE
        Temp = Temp0 * EXP(-Tempfac*(CYCLE-1))      ! lower the temp
         ! Adjust the maximum distance a boundary line is perturbed with the cycle number
        u = FLOAT(CYCLE-1)/(NCYCLE-1)
        DelLines = (1-u)*0.5*MAX(NPXT/Nxd,NPYT/Nyd) + u*1.0

         ! Do NITER iterations per temp cycle
        DO ITER = 1, NITER
           ! Perturb the boundary lines
          TrialLines(:)=CurLines(:)
          DO WHILE (ALL(TrialLines(:)==CurLines(:)))
            DO I = 1, N
              TrialLines(I) = NINT(CurLines(I) + (2*ran(iseed)-1.0)*DelLines)
            ENDDO
          ENDDO
           ! Calculate the energy of these new boundary lines
          NewEnergy = LOAD_OBJ_FUNC (TrialLines, Nxd, Nyd, ProcCoordAll, wtproc, NPXT, NPYT, Ngrid)
          DeltaE = NewEnergy - Energy
          IF (DeltaE < 0.0) THEN   ! If lower energy accept change
            CurLines(:) = TrialLines(:)
            Energy = NewEnergy
          ELSE                      !  otherwise sometimes accept it anyway
            Boltz = EXP(-DeltaE/Temp)
            IF (ran(iseed) < Boltz) THEN
              CurLines(:) = TrialLines(:)
              Energy = NewEnergy
            ENDIF
          ENDIF
          IF (Energy < LowestEnergy) THEN
            BndLines(:) = CurLines(:)
            LowestEnergy = Energy
          ENDIF
        ENDDO
        ! print '(I3,2(1X,F9.6),10(1X,I3))', CYCLE,Temp,LowestEnergy, BndLines(:)
      ENDDO  ! end of temperature cycles

      DEALLOCATE (Ngrid)
    ENDIF

    CALL MPI_BCAST (BndLines, N, MPI_INTEGER, 0, comm2d, ierr)
    DEALLOCATE (ProcCoordAll)
    IF (myproc==0 .AND. DoLoadBalancing) THEN
      WRITE (6,'(A)') 'Optimum load balanced boundary lines in the property grid:'
      WRITE (6,'(A,16I4)') '  X lines:',BndLines(1:Nxd-1)
      WRITE (6,'(A,16I4)') '  Y lines:',BndLines(Nxd:N)
    ENDIF
  ENDIF

   ! Make the subdomain position array for this processor from BndLines
  IDX(1)=1 ; IDX(Nxd+1)=NPXT+1 ; IF (Nxd==1) IDX(Nxd+1)=NPXT
  IDY(1)=1 ; IDY(Nyd+1)=NPYT+1 ; IF (Nyd==1) IDY(Nyd+1)=NPYT
  IDX(2:Nxd)=BndLines(1:Nxd-1) ; IDY(2:Nyd)=BndLines(Nxd:N)
  IXYPRP(1:2,1)=IDX(ProcCoord(1)+1:ProcCoord(1)+2)
  IXYPRP(1:2,2)=IDY(ProcCoord(2)+1:ProcCoord(2)+2)
END SUBROUTINE OPTIMIZE_PROCESSOR_DOMAINS



REAL FUNCTION LOAD_OBJ_FUNC (BndLines, Nxd, Nyd, ProcCoordAll, wtproc, &
                             NPXT, NPYT, Ngrid)
 ! Calculates the objective function to minimize in doing the processor
 ! load balancing optimization.  The boundary line indices in X and Y to
 ! the property grid, which is the optimization control vector is input
 ! in BndLines.  The objective function is the normalized sum of the 
 ! squares of the difference between the number of grid points (Ngrid)
 ! of each domain and the mean number, weighted by the relative processor
 ! times.
  INTEGER, INTENT(IN) :: Nxd, Nyd, BndLines(Nxd+Nyd-2)
  INTEGER, INTENT(IN) :: ProcCoordAll(2,0:Nxd*Nyd-1)
  REAL,    INTENT(IN) :: wtproc(0:Nxd*Nyd)
  INTEGER, INTENT(IN) :: NPXT, NPYT, Ngrid(NPXT+1,NPYT+1)
  INTEGER :: N, I, J, IX, IY, IDX(Nxd+1), IDY(Nyd+1)
  REAL    :: OBJA
  REAL, SAVE :: sumN=-1.0, a

  IF (sumN < 0) THEN
    sumN = SUM(Ngrid(:,:)) ;  a = FLOAT(Nxd*Nyd)/sumN
  ENDIF
  N=Nxd+Nyd-2
  IDX(1)=1 ; IDX(Nxd+1)=NPXT+1 ; IF (Nxd==1) IDX(Nxd+1)=NPXT
  IDY(1)=1 ; IDY(Nyd+1)=NPYT+1 ; IF (Nyd==1) IDY(Nyd+1)=NPYT

  OBJA = 0.0
  DO J = 0, Nxd*Nyd-1
    IDX(2:Nxd)=BndLines(1:Nxd-1) ; IDY(2:Nyd)=BndLines(Nxd:N)
    IX = ProcCoordAll(1,J)+1 ; IY = ProcCoordAll(2,J)+1
    OBJA = OBJA + (1.0/Nxd*Nyd) * &
         (a*wtproc(J)*SUM(Ngrid(IDX(IX):IDX(IX+1),IDY(IY):IDY(IY+1))) - 1.0)**2
  ENDDO
  LOAD_OBJ_FUNC = OBJA
END FUNCTION LOAD_OBJ_FUNC




REAL FUNCTION ran(idum)
 ! Returns uniform random deviate between 0.0 and 1.0 (exclusive of endpoints).
 ! Call ran with idum negative to initialize; thereafter, do not alter idum.
  IMPLICIT NONE
  INTEGER, PARAMETER :: K4B=SELECTED_INT_KIND(9)
  INTEGER(K4B), INTENT(INOUT) :: idum
  INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
  REAL, SAVE :: am
  INTEGER(K4B), SAVE :: ix=-1, iy=-1, k
  
  if (idum <= 0 .or. iy < 0) then
    am = nearest(1.0,-1.0)/IM
    iy = ior(ieor(888889999,abs(idum)),1)
    ix = ieor(777755555,abs(idum))
    idum = abs(idum)+1
  endif
  ix = ieor(ix,ishft(ix,13))
  ix = ieor(ix,ishft(ix,-17))
  ix = ieor(ix,ishft(ix,5))
  k = iy/IQ
  iy = IA*(iy-k*IQ)-IR*k
  if (iy < 0) iy=iy+IM
  ran = am*ior(iand(IM,ieor(ix,iy)),1)
END FUNCTION ran




SUBROUTINE BROADCAST_USER_INPUT (RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
                       INSAVEFILE, OUTSAVEFILE,                      &
                       NX, NY, NZ, NMU, NPHI, BCFLAG,                &
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
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
  INTEGER, INTENT(IN)    :: MAXOUT, MAXPAR
  INTEGER, INTENT(INOUT) :: MAXITER, NUMOUT
  LOGICAL, INTENT(INOUT) :: KDIST, DELTAM, ACCELFLAG
  REAL,    INTENT(INOUT) :: SOLARFLUX, SOLARMU, SOLARAZ
  REAL,    INTENT(INOUT) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
  REAL,    INTENT(INOUT) :: SOLACC, SPLITACC, SHACC, OUTPARMS(MAXPAR,MAXOUT)
  REAL,    INTENT(INOUT) :: MAX_TOTAL_MB, SPLITTING_FACTOR
  REAL,    INTENT(INOUT) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
  CHARACTER(LEN=1), INTENT(INOUT) :: SRCTYPE, GRIDTYPE, UNITS, OUTTYPES(*)
  CHARACTER(LEN=64), INTENT(INOUT) :: RUNNAME, PROPFILE, SFCFILE, CKDFILE
  CHARACTER(LEN=64), INTENT(INOUT) :: INSAVEFILE, OUTSAVEFILE
  include 'mpif.h'
  INTEGER :: ierr, intbuf(9)
  LOGICAL :: logbuf(3)
  REAL    :: realbuf(16)
  CHARACTER(LEN=1) :: char1buf(3)
  CHARACTER(LEN=64) :: char64buf(6)
  CHARACTER(LEN=3) :: procstr

  IF (numproc <= 1) RETURN

  intbuf(:) = (/ NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG, MAXITER, NUMOUT /)
  CALL MPI_BCAST (intbuf, 9, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  NX=intbuf(1) ; NY=intbuf(2) ; NZ=intbuf(3) ; NMU=intbuf(4) ; NPHI=intbuf(5)
  BCFLAG=intbuf(6) ; IPFLAG=intbuf(7) ; MAXITER=intbuf(8) ; NUMOUT=intbuf(9)

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

  char64buf(:) = (/ PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, RUNNAME /)
  CALL MPI_BCAST (char64buf, 6*64, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  PROPFILE=char64buf(1) ; SFCFILE=char64buf(2) ; CKDFILE=char64buf(3)
  INSAVEFILE=char64buf(4) ; OUTSAVEFILE=char64buf(5) ; RUNNAME=char64buf(6)

  CALL MPI_BCAST (OUTPARMS, MAXPAR*NUMOUT, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST (OUTTYPES, NUMOUT*1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  WRITE (procstr,'(I3.3)') myproc
  IF (INSAVEFILE(1:4) .NE. 'NONE') INSAVEFILE = TRIM(INSAVEFILE)//procstr
  IF (OUTSAVEFILE(1:4) .NE. 'NONE') OUTSAVEFILE = TRIM(OUTSAVEFILE)//procstr
END SUBROUTINE BROADCAST_USER_INPUT



SUBROUTINE BROADCAST_PROPERTY_SIZE (NPX, NPY, NPZ, DELX, DELY, &
                                    NUMPHASE, MAXLEG, MAXPGL)
 ! Broadcasts the property file size information from the master processor
 ! to the slave processors.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
!f2py intent(in, out) :: NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
  REAL,    INTENT(INOUT) :: DELX, DELY
!f2py intent(in, out) :: DELX, DELY
  include 'mpif.h'
  INTEGER :: ierr, intbuf(6)
  REAL    :: realbuf(2)

  IF (numproc <= 1) RETURN

  intbuf(:) = (/ NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL /)
  CALL MPI_BCAST (intbuf, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  NPX=intbuf(1) ; NPY=intbuf(2) ; NPZ=intbuf(3)
  NUMPHASE=intbuf(4) ; MAXLEG=intbuf(5) ; MAXPGL=intbuf(6)

  realbuf(:) = (/ DELX, DELY /)
  CALL MPI_BCAST (realbuf, 2, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  DELX=realbuf(1) ; DELY=realbuf(2)
END SUBROUTINE BROADCAST_PROPERTY_SIZE



SUBROUTINE SCATTER_PROPERTIES (NPXT, NPYT, NPX, NPY, NPZ, NLEG, NUMPHASE, &
                          ZLEVELS, MAXASYM, &
                          TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT, &
                          TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP)
 ! Sends portions of the optical properties on the property grid to each
 ! of the processors.  Deals with the duplicating the grid points on the
 ! boundaries between processors.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NPXT, NPYT, NPX, NPY, NPZ, NUMPHASE
  INTEGER, INTENT(INOUT) :: NLEG
!f2py intent(in, out) :: NLEG
  REAL,    INTENT(INOUT) :: ZLEVELS(NPZ), MAXASYM
!f2py intent(in, out) :: ZLEVELS, MAXASYM
  REAL,    INTENT(IN) :: TEMPPT(NPZ,NPYT,NPXT), EXTINCTPT(NPZ,NPYT,NPXT)
  REAL,    INTENT(IN) :: ALBEDOPT(NPZ,NPYT,NPXT), LEGENPT(*)
  INTEGER*2, INTENT(IN) :: IPHASEPT(NPZ,NPYT,NPXT)
  REAL,    INTENT(OUT) :: TEMPP(NPZ,NPY,NPX), EXTINCTP(NPZ,NPY,NPX) 
  REAL,    INTENT(OUT) :: ALBEDOP(NPZ,NPY,NPX), LEGENP(*)
  INTEGER*2, INTENT(OUT) :: IPHASEP(NPZ,NPY,NPX)
  include 'mpif.h'
  INTEGER :: SX, EX, NX, SY, EY, NY, N, MAXN, MaxNX, MaxNY
  INTEGER :: iproc, irq, ierr
  INTEGER, ALLOCATABLE :: IXYALL(:,:,:), requests(:), status(:,:)
  REAL,    ALLOCATABLE :: realbuf(:), rbuf(:,:,:)
  REAL,    ALLOCATABLE :: tempbuf(:,:), extbuf(:,:), albbuf(:,:)
  INTEGER*2, ALLOCATABLE :: ibuf(:,:,:), iphbuf(:,:)

   ! Broadcast the small stuff that is going to all processors
  ALLOCATE (realbuf(NPZ+1))
  realbuf(1:NPZ+1) = (/ ZLEVELS(1:NPZ), MAXASYM /)
  CALL MPI_BCAST (realbuf, NPZ+1, MPI_REAL, 0, comm2d, ierr)
  ZLEVELS(1:NPZ)=realbuf(1:NPZ) ; MAXASYM=realbuf(NPZ+1)

   ! If we have a tabulated phase function, then broadcast to all processors
  IF (NUMPHASE > 0) THEN
    IF (myproc == 0) LEGENP(1:NUMPHASE*NLEG) = LEGENPT(1:NUMPHASE*NLEG)
    CALL MPI_BCAST (NLEG, 1, MPI_INTEGER, 0, comm2d, ierr)
    CALL MPI_BCAST (LEGENP, NLEG*NUMPHASE, MPI_REAL, 0, comm2d, ierr)
  ENDIF

   ! Get the IXYPRP arrays from all processors to tell where they are in property grid
  ALLOCATE (IXYALL(2,2,0:numproc-1))
  CALL MPI_GATHER (IXYPRP,4,MPI_INTEGER, IXYALL,4,MPI_INTEGER, 0, comm2d, ierr)

  ALLOCATE (requests(4*numproc+4), status(MPI_STATUS_SIZE,4*numproc+4))

  irq = 0
  IF (myproc == 0) THEN
     ! Allocate the send buffers for each processor
    MAXN=1 ; MaxNX=1 ; MaxNY=1
    DO iproc = 0, numproc-1
      NX=IXYALL(2,1,iproc)-IXYALL(1,1,iproc)+1 ; MaxNX=MAX(MaxNX,NX)
      NY=IXYALL(2,2,iproc)-IXYALL(1,2,iproc)+1 ; MaxNY=MAX(MaxNY,NY)
      MAXN = MAX(MAXN,NPZ*NY*NX)
    ENDDO
    ALLOCATE (rbuf(NPZ,MaxNY,MaxNX), ibuf(NPZ,MaxNY,MaxNX))
    ALLOCATE (tempbuf(MAXN,0:numproc-1), extbuf(MAXN,0:numproc-1))
    ALLOCATE (albbuf(MAXN,0:numproc-1), iphbuf(MAXN,0:numproc-1))
     ! The master processor reformats the four property arrays (to deal with
     ! the wrap around at the full domain edge) and sends the arrays with MPI.
    DO iproc = 0, numproc-1
      SX=IXYALL(1,1,iproc) ; EX=IXYALL(2,1,iproc)
      SY=IXYALL(1,2,iproc) ; EY=IXYALL(2,2,iproc)
      NX=EX-SX+1 ; NY=EY-SY+1 ; N=NX*NY*NPZ
      rbuf(:,1:NY-1,1:NX-1) = TEMPPT(:,SY:EY-1,SX:EX-1)
      rbuf(:,1:NY-1,NX) = TEMPPT(:,SY:EY-1,MOD(EX-1,NPXT)+1)
      rbuf(:,NY,1:NX-1) = TEMPPT(:,MOD(EY-1,NPYT)+1,SX:EX-1)
      rbuf(:,NY,NX) = TEMPPT(:,MOD(EY-1,NPYT)+1,MOD(EX-1,NPXT)+1)
      tempbuf(1:N,iproc) = RESHAPE(rbuf(1:NPZ,1:NY,1:NX), (/ N /))
      CALL MPI_ISEND (tempbuf(1:N,iproc), N, MPI_REAL, iproc, 1, comm2d, requests(irq+1), ierr)

      rbuf(:,1:NY-1,1:NX-1) = EXTINCTPT(:,SY:EY-1,SX:EX-1)
      rbuf(:,1:NY-1,NX) = EXTINCTPT(:,SY:EY-1,MOD(EX-1,NPXT)+1)
      rbuf(:,NY,1:NX-1) = EXTINCTPT(:,MOD(EY-1,NPYT)+1,SX:EX-1)
      rbuf(:,NY,NX) = EXTINCTPT(:,MOD(EY-1,NPYT)+1,MOD(EX-1,NPXT)+1)
      extbuf(1:N,iproc) = RESHAPE(rbuf(1:NPZ,1:NY,1:NX), (/ N /))
      CALL MPI_ISEND (extbuf(1:N,iproc), N, MPI_REAL, iproc, 2, comm2d, requests(irq+2), ierr)

      rbuf(:,1:NY-1,1:NX-1) = ALBEDOPT(:,SY:EY-1,SX:EX-1)
      rbuf(:,1:NY-1,NX) = ALBEDOPT(:,SY:EY-1,MOD(EX-1,NPXT)+1)
      rbuf(:,NY,1:NX-1) = ALBEDOPT(:,MOD(EY-1,NPYT)+1,SX:EX-1)
      rbuf(:,NY,NX) = ALBEDOPT(:,MOD(EY-1,NPYT)+1,MOD(EX-1,NPXT)+1)
      albbuf(1:N,iproc) = RESHAPE(rbuf(1:NPZ,1:NY,1:NX), (/ N /))
      CALL MPI_ISEND (albbuf(1:N,iproc), N, MPI_REAL, iproc, 3, comm2d, requests(irq+3), ierr)

      ibuf(:,1:NY-1,1:NX-1) = IPHASEPT(:,SY:EY-1,SX:EX-1)
      ibuf(:,1:NY-1,NX) = IPHASEPT(:,SY:EY-1,MOD(EX-1,NPXT)+1)
      ibuf(:,NY,1:NX-1) = IPHASEPT(:,MOD(EY-1,NPYT)+1,SX:EX-1)
      ibuf(:,NY,NX) = IPHASEPT(:,MOD(EY-1,NPYT)+1,MOD(EX-1,NPXT)+1)
      iphbuf(1:N,iproc) = RESHAPE(ibuf(1:NPZ,1:NY,1:NX), (/ N /))
      CALL MPI_ISEND (iphbuf(1:N,iproc), N, MPI_INTEGER2, iproc, 4, comm2d, requests(irq+4), ierr)
      irq = irq + 4
    ENDDO
  ENDIF

  N=NPX*NPY*NPZ
  CALL MPI_IRECV (TEMPP, N, MPI_REAL, 0, 1, comm2d, requests(irq+1), ierr)
  CALL MPI_IRECV (EXTINCTP, N, MPI_REAL, 0, 2, comm2d, requests(irq+2), ierr)
  CALL MPI_IRECV (ALBEDOP, N, MPI_REAL, 0, 3, comm2d, requests(irq+3), ierr)
  CALL MPI_IRECV (IPHASEP, N, MPI_INTEGER2, 0, 4, comm2d, requests(irq+4), ierr)
  irq = irq + 4
  CALL MPI_WAITALL (irq, requests, status, ierr)

  IF (myproc == 0) DEALLOCATE (tempbuf, extbuf, albbuf, iphbuf, rbuf, ibuf)
  DEALLOCATE (requests, status, IXYALL, realbuf)
END SUBROUTINE SCATTER_PROPERTIES




SUBROUTINE BROADCAST_SURFACE_SIZE (MAXSFCPTS, MAXSFCPARS)
 ! Broadcasts the surface array sizes from the master process to the others.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: MAXSFCPTS, MAXSFCPARS
  INTEGER :: ierr, intbuf(2)
  include 'mpif.h'

  IF (numproc <= 1) RETURN

  intbuf(1:2) = (/ MAXSFCPTS, MAXSFCPARS/)
  CALL MPI_BCAST (intbuf, 2, MPI_INTEGER, 0, comm2d, ierr)
  MAXSFCPTS=intbuf(1) ; MAXSFCPARS=intbuf(2)
END SUBROUTINE BROADCAST_SURFACE_SIZE


SUBROUTINE BROADCAST_SURFACE_PARMS (SFCTYPE, NXSFC, NYSFC, NSFCPAR, &
                                DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)
 ! Broadcasts the surface file parameters from the master process to the others.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  CHARACTER(LEN=2), INTENT(INOUT) :: SFCTYPE
!f2py intent(in, out) :: SFCTYPE
  INTEGER, INTENT(INOUT) :: NXSFC, NYSFC, NSFCPAR
!f2py intent(in, out) :: NXSFC, NYSFC, NSFCPAR
  REAL,    INTENT(INOUT) :: DELXSFC, DELYSFC, SFCPARMS(*), GNDTEMP, GNDALBEDO
!f2py intent(in, out) :: DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO
  INTEGER :: ierr, intbuf(3), n
  REAL    :: realbuf(4)
  include 'mpif.h'

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
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: NG, NZCKD
  INTEGER :: ierr, intbuf(2)
  include 'mpif.h'

  IF (numproc <= 1) RETURN

  intbuf(1:2) = (/ NG, NZCKD/)
  CALL MPI_BCAST (intbuf, 2, MPI_INTEGER, 0, comm2d, ierr)
  NG=intbuf(1) ; NZCKD=intbuf(2)
END SUBROUTINE BROADCAST_KDIST_SIZE


SUBROUTINE BROADCAST_KDIST_PARMS (SOLFLUX, NG, DELG, NZCKD, ZCKD, KABS)
 ! Broadcasts the k-distribution file parameters from the master process 
 ! to the other processes.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NG, NZCKD
  REAL,    INTENT(INOUT) :: SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD,NG)
  INTEGER :: ierr
  REAL, ALLOCATABLE :: realbuf(:)
  include 'mpif.h'

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
  IMPLICIT NONE
  REAL, INTENT(INOUT) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
!f2py intent(in, out) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
  REAL, INTENT(INOUT) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
!f2py intent(in, out) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
  CHARACTER(LEN=*), INTENT(IN) :: RUNNAME
!f2py intent(in) :: RUNNAME
  INTEGER :: Nproc, ierr, i, j
  REAL, ALLOCATABLE :: MemParam(:,:)
  include 'mpif.h'

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



REAL FUNCTION SUM_CPU_TIME (cpuTimes)
 ! Sums the CPU time used across all processors
  USE SHDOM_MPI_DATA
  real cpuTimes
  include 'mpif.h'
  integer ierr

  IF (numproc>1) THEN
    CALL MPI_REDUCE(cpuTimes, SUM_CPU_TIME, 1, MPI_REAL, MPI_SUM, 0, comm2d, ierr)
  ELSE
    SUM_CPU_TIME = cpuTimes
  ENDIF
END FUNCTION SUM_CPU_TIME



SUBROUTINE GATHER_OUTPUT (FLUX_OUT, FLUXDIV_OUT, SH_OUT, &
                          IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, &
                          SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS, &
                          NPXT, NPYT, DELX, DELY, XALLGRID, YALLGRID, &
                          ALLFLUXES, ALLFLUXDIV, ALLSHTERMS, &
                          NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)
 ! Gathers the base grid output data in SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, 
 ! and SUMSHTERMS from all the processors and assembles in the master 
 ! processor output arrays ALLFLUXES, ALLFLUXDIV, ALLSHTERMS.  Gathers the
 ! IXYOUT array from each processor to determine where each subdomain is.
 ! Also gathers the full base grid position arrays (XALLGRID/YALLGRID) and
 ! sums NCELLS,NPTS,NSH over all processors.
  USE SHDOM_MPI_DATA
  LOGICAL, INTENT(IN) :: FLUX_OUT, FLUXDIV_OUT, SH_OUT
  INTEGER, INTENT(IN) :: IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, NPXT, NPYT
  REAL,    INTENT(IN) :: DELX, DELY
  REAL,    INTENT(IN) :: SUMFLUXES(2,NBPTS), SUMDIRFLUX(NBPTS)
  REAL,    INTENT(IN) :: SUMFLUXDIV(NBPTS), SUMSHTERMS(NSHOUT,NBPTS)
  REAL,    INTENT(OUT) :: XALLGRID(NXT), YALLGRID(NYT)
  REAL,    INTENT(OUT) :: ALLFLUXES(3,NZ,NYT,NXT), ALLFLUXDIV(NZ,NYT,NXT)
  REAL,    INTENT(OUT) :: ALLSHTERMS(NSHOUT,NZ,NYT,NXT)
  INTEGER, INTENT(IN) :: NCELLS, NPTS, NSH
  INTEGER, INTENT(OUT) :: NCELLSTOT, NPTSTOT, NSHTOT
  include 'mpif.h'
  INTEGER :: IX, IY, NX, NY, NX1, NY1, MAXN, I, send3int(3), recv3int(3)
  INTEGER :: iproc, ierr, irq
  INTEGER, ALLOCATABLE :: IXYALL(:,:,:), requests(:), status(:,:)
  REAL,    ALLOCATABLE :: sendbuf(:,:), recvbuf(:,:,:), buf(:,:,:,:)

  ALLOCATE (IXYALL(2,2,0:numproc-1))
  CALL MPI_GATHER (IXYOUT,4,MPI_INTEGER, IXYALL,4,MPI_INTEGER, 0, comm2d, ierr)

  ALLOCATE (requests(numproc), status(MPI_STATUS_SIZE,numproc))
  
  IF (FLUX_OUT) THEN
    irq=0
    IF (myproc > 0) THEN
      ALLOCATE (sendbuf(3,NBPTS))
      sendbuf(1:2,:) = SUMFLUXES(:,:) ; sendbuf(3,:) = SUMDIRFLUX(:)
      irq = irq + 1
      CALL MPI_ISEND (sendbuf, 3*NBPTS, MPI_REAL, 0, 1, comm2d, &
                      requests(irq), ierr)
    ELSE
      MAXN = (MAXVAL(IXYALL(2,1,:))+1)*(MAXVAL(IXYALL(2,2,:))+1)*NZ
      ALLOCATE (recvbuf(3,MAXN,1:numproc-1))
      DO iproc = 1, numproc-1
        irq = irq + 1
        CALL MPI_IRECV (recvbuf(:,:,iproc), 3*MAXN, MPI_REAL, iproc, 1, &
                        comm2d, requests(irq), ierr)
      ENDDO
    ENDIF
    CALL MPI_WAITALL (irq, requests, status, ierr)

    IF (myproc == 0) THEN
      IX = IXYALL(1,1,0) ; NX = IXYALL(2,1,0)
      IY = IXYALL(1,2,0) ; NY = IXYALL(2,2,0)
      NX1=NX+1 ; NY1=NY+1
      IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
      IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
      ALLOCATE (buf(3,NZ,NY1,NX1))
      buf(1:2,:,:,:) = RESHAPE(SUMFLUXES(1:2,:),(/ 2,NZ,NY1,NX1 /) )
      buf(3,:,:,:) = RESHAPE(SUMDIRFLUX(:),(/ NZ,NY1,NX1 /) )
      ALLFLUXES(1:3,1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(:,1:NZ,1:NY,1:NX)
      DEALLOCATE (buf)
      DO iproc = 1, numproc-1
        IX = IXYALL(1,1,iproc) ; NX = IXYALL(2,1,iproc)
        IY = IXYALL(1,2,iproc) ; NY = IXYALL(2,2,iproc)
        NX1=NX+1 ; NY1=NY+1
        IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
        IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
        ALLOCATE (buf(3,NZ,NY1,NX1))
        buf(:,:,:,:) = RESHAPE(recvbuf(:,1:NZ*NY1*NX1,iproc),(/ 3,NZ,NY1,NX1 /) )
        ALLFLUXES(1:3,1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(:,1:NZ,1:NY,1:NX)
        DEALLOCATE (buf)
      ENDDO
      DEALLOCATE (recvbuf)
    ELSE
      DEALLOCATE (sendbuf)
    ENDIF
  ENDIF


  IF (FLUXDIV_OUT) THEN
    irq=0
    IF (myproc > 0) THEN
      irq = irq + 1
      CALL MPI_ISEND (SUMFLUXDIV, NBPTS, MPI_REAL, 0, 2, comm2d, &
                      requests(irq), ierr)
    ELSE
      MAXN = (MAXVAL(IXYALL(2,1,:))+1)*(MAXVAL(IXYALL(2,2,:))+1)*NZ
      ALLOCATE (recvbuf(1,MAXN,1:numproc-1))
      DO iproc = 1, numproc-1
        irq = irq + 1
        CALL MPI_IRECV (recvbuf(:,:,iproc), MAXN, MPI_REAL, iproc, 2, &
                        comm2d, requests(irq), ierr)
      ENDDO
    ENDIF
    CALL MPI_WAITALL (irq, requests, status, ierr)

    IF (myproc == 0) THEN
      IX = IXYALL(1,1,0) ; NX = IXYALL(2,1,0)
      IY = IXYALL(1,2,0) ; NY = IXYALL(2,2,0)
      NX1=NX+1 ; NY1=NY+1
      IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
      IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
      ALLOCATE (buf(1,NZ,NY1,NX1))
      buf(1,:,:,:) = RESHAPE(SUMFLUXDIV(:),(/ NZ,NY1,NX1 /) )
      ALLFLUXDIV(1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(1,1:NZ,1:NY,1:NX)
      DEALLOCATE (buf)
      DO iproc = 1, numproc-1
        IX = IXYALL(1,1,iproc) ; NX = IXYALL(2,1,iproc)
        IY = IXYALL(1,2,iproc) ; NY = IXYALL(2,2,iproc)
        NX1=NX+1 ; NY1=NY+1
        IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
        IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
        ALLOCATE (buf(1,NZ,NY1,NX1))
        buf(:,:,:,:) = RESHAPE(recvbuf(:,1:NZ*NY1*NX1,iproc),(/ 1,NZ,NY1,NX1 /) )
        ALLFLUXDIV(1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(1,1:NZ,1:NY,1:NX)
        DEALLOCATE (buf)
      ENDDO
      DEALLOCATE (recvbuf)
    ENDIF
  ENDIF


  IF (SH_OUT) THEN
    irq=0
    IF (myproc > 0) THEN
      irq = irq + 1
      CALL MPI_ISEND (SUMSHTERMS, NSHOUT*NBPTS, MPI_REAL, 0, 3, comm2d, &
                      requests(irq), ierr)
    ELSE
      MAXN = (MAXVAL(IXYALL(2,1,:))+1)*(MAXVAL(IXYALL(2,2,:))+1)*NZ
      ALLOCATE (recvbuf(NSHOUT,MAXN,1:numproc-1))
      DO iproc = 1, numproc-1
        irq = irq + 1
        CALL MPI_IRECV (recvbuf(:,:,iproc), NSHOUT*MAXN, MPI_REAL, iproc, 3, &
                        comm2d, requests(irq), ierr)
      ENDDO
    ENDIF
    CALL MPI_WAITALL (irq, requests, status, ierr)

    IF (myproc == 0) THEN
      IX = IXYALL(1,1,0) ; NX = IXYALL(2,1,0)
      IY = IXYALL(1,2,0) ; NY = IXYALL(2,2,0)
      NX1=NX+1 ; NY1=NY+1
      IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
      IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
      ALLOCATE (buf(NSHOUT,NZ,NY1,NX1))
      buf(:,:,:,:) = RESHAPE(SUMSHTERMS(:,:),(/ NSHOUT,NZ,NY1,NX1 /) )
      ALLSHTERMS(:,1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(:,1:NZ,1:NY,1:NX)
      DEALLOCATE (buf)
      DO iproc = 1, numproc-1
        IX = IXYALL(1,1,iproc) ; NX = IXYALL(2,1,iproc)
        IY = IXYALL(1,2,iproc) ; NY = IXYALL(2,2,iproc)
        NX1=NX+1 ; NY1=NY+1
        IF (BTEST(IPFLAG,0) .AND. .NOT.BTEST(BCFLAG,2)) NX1=NX
        IF (BTEST(IPFLAG,1) .AND. .NOT.BTEST(BCFLAG,3)) NY1=NY
        ALLOCATE (buf(NSHOUT,NZ,NY1,NX1))
        buf(:,:,:,:) = RESHAPE(recvbuf(:,1:NZ*NY1*NX1,iproc),(/ NSHOUT,NZ,NY1,NX1 /) )
        ALLSHTERMS(:,1:NZ,IY:IY+NY-1,IX:IX+NX-1) = buf(:,1:NZ,1:NY,1:NX)
        DEALLOCATE (buf)
      ENDDO
      DEALLOCATE (recvbuf)
    ENDIF
  ENDIF


   ! Make the X and Y grids for the full domain
  IF (myproc == 0) THEN
    XALLGRID(1:NXT) = (/ (NPXT*DELX*(I-1)/NXT, I=1, NXT) /)
    YALLGRID(1:NYT) = (/ (NPYT*DELY*(I-1)/NYT, I=1, NYT) /)
  ENDIF

   ! Get the sum of the NCELLS, NPTS, and NSH over the processors
  send3int = (/ NCELLS, NPTS, NSH /)
  CALL MPI_REDUCE (send3int, recv3int, 3, MPI_INTEGER, MPI_SUM, 0, comm2d, ierr)
  NCELLSTOT=recv3int(1) ; NPTSTOT=recv3int(2) ; NSHTOT=recv3int(3)

  DEALLOCATE (IXYALL, requests, status)
END SUBROUTINE GATHER_OUTPUT


SUBROUTINE TOTAL_ALBEDO_MAX (ALBMAX)
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  REAL,    INTENT(INOUT) :: ALBMAX
  INTEGER :: ierr
  REAL    :: sendrbuf
  include 'mpif.h'

  IF (numproc > 1) THEN
     ! Call MPI_ALLREDUCE to find the max ALBMAX over all processors 
    sendrbuf = ALBMAX
    call MPI_ALLREDUCE (sendrbuf, ALBMAX, 1, MPI_REAL, MPI_MAX, comm2d, ierr)
    if (ierr /= MPI_SUCCESS) CALL ABORT_SHDOM_MPI ('TOTAL_ALBEDO_MAX: MPI_ALLREDUCE error')
  ENDIF
END SUBROUTINE TOTAL_ALBEDO_MAX


SUBROUTINE UNIFY_SPLITTING (DOSPLIT, STARTSPLITACC)
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  LOGICAL, INTENT(INOUT) :: DOSPLIT
  REAL,    INTENT(INOUT) :: STARTSPLITACC
  INTEGER :: ierr(2)
  LOGICAL :: sendlbuf
  REAL    :: sendrbuf
  include 'mpif.h'

  IF (numproc > 1) THEN
     ! Call MPI_ALLREDUCE to 
    sendrbuf = STARTSPLITACC
    call MPI_ALLREDUCE (sendrbuf, STARTSPLITACC, 1, MPI_REAL, MPI_MAX, comm2d, ierr(1))
    sendlbuf = DOSPLIT
    call MPI_ALLREDUCE (sendlbuf, DOSPLIT, 1, MPI_LOGICAL, MPI_LOR, comm2d, ierr(2))
    if (ANY(ierr(:) /= MPI_SUCCESS)) CALL ABORT_SHDOM_MPI ('UNIFY_SPLITTING: MPI_ALLREDUCE error')
  ENDIF
END SUBROUTINE UNIFY_SPLITTING


SUBROUTINE TOTAL_SPLITCRIT_MAX (SPLITCRIT)
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  REAL,    INTENT(INOUT) :: SPLITCRIT
  INTEGER :: ierr
  REAL    :: sendrbuf
  include 'mpif.h'

  IF (numproc > 1) THEN
     ! Call MPI_ALLREDUCE to find the max SPLITCRIT over all processors 
    sendrbuf = SPLITCRIT
    call MPI_ALLREDUCE (sendrbuf, SPLITCRIT, 1, MPI_REAL, MPI_MAX, comm2d, ierr)
    if (ierr /= MPI_SUCCESS) CALL ABORT_SHDOM_MPI ('TOTAL_SPLITCRIT_MAX: MPI_ALLREDUCE error')
  ENDIF
END SUBROUTINE TOTAL_SPLITCRIT_MAX




SUBROUTINE MAKE_DIRECT_PAR (SPT, NPTS, BCFLAG, IPFLAG, DELTAM, ML, NLEG, &
                            SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS, &
                            NX, XGRID, NY, YGRID,  DIRFLUX)
 ! Makes the direct beam solar flux for the NPTS gridpoints started with 
 ! index SPT having X/Y/Z positions in GRIDPOS.
 ! DIRFLUX is set to F*exp(-tau_sun).
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: SPT, NPTS, BCFLAG, IPFLAG, ML, NLEG, NX, NY
  LOGICAL, INTENT(IN) :: DELTAM
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, SOLARAZ
  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS), XGRID(NX+1), YGRID(NY+1)
  REAL,    INTENT(OUT) :: DIRFLUX(NPTS)
  include 'mpif.h'
  LOGICAL :: VALIDBEAM
  REAL    :: UNIFZLEV, UZLEV, DIRPATH
  INTEGER :: MAXNBEAM, Lupx, Lupy, Ldownx, Ldowny, Nbeam
  INTEGER :: SIDE, NINVALID(0:4)
  INTEGER :: IP, I, J, L
  INTEGER :: ierr, irq, requests(12), status(MPI_STATUS_SIZE,12), iproc
  INTEGER :: intsendbuf(2), intrecvbuf(2)
  LOGICAL :: VALIDRAD, BTEST, AllDone
  REAL    :: XO, YO, ZO, X, Y, Z, XMAX, YMAX, EPS, DIRJUNK
  INTEGER, ALLOCATABLE :: Ndone(:)
  INTEGER, ALLOCATABLE :: BndDirIntInfo(:,:,:)
  REAL, ALLOCATABLE    :: BndDirRealInfo(:,:,:)
  INTEGER, ALLOCATABLE :: DirPtrOut(:,:), DirPtrIn(:)
  REAL, ALLOCATABLE    :: DirAllOut(:,:), DirAllIn(:)


  MAXNBEAM = 2*NPTS
  ALLOCATE (BndDirRealInfo(4,MAXNBEAM,4), BndDirIntInfo(2,MAXNBEAM,4))
  ALLOCATE (DirAllOut(MAXNBEAM,0:numproc-1), DirPtrOut(MAXNBEAM,0:numproc-1))
  ALLOCATE (DirAllIn(MAXNBEAM), DirPtrIn(MAXNBEAM), Ndone(0:numproc-1))

  EPS = 1.0E-5*MAX(XGRID(2)-XGRID(1),YGRID(2)-YGRID(1))

  IF (SPT == 1) THEN
     ! Initialize the direct beam calculation
    CALL DIRECT_BEAM_PROP (1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG, &
             DELTAM,ML,NLEG, SOLARFLUX,SOLARMU,SOLARAZ, DIRFLUX(1), &
             UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
     ! Find the max lowest uniform Z level over all processors and then
     !   tell DIRECT_BEAM_PROP about it
    UZLEV = UNIFZLEV
    CALL MPI_ALLREDUCE (UZLEV, UNIFZLEV, 1, MPI_REAL, MPI_MAX, comm2d, ierr)
    CALL DIRECT_BEAM_PROP (2, 0.0, 0.0, 0.0, BCFLAG, IPFLAG, &
             DELTAM,ML,NLEG, SOLARFLUX,SOLARMU,SOLARAZ, DIRFLUX(1), &
             UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
  ENDIF


  Lupx=0 ; Lupy=0 ; Ldownx=0 ; Ldowny=0
  IF (BTEST(BCFLAG,2)) THEN
    IF (COS(SOLARAZ) .GT. 1.0E-5) THEN
      Ldownx = 2 ; Lupx = 1
    ELSE
      Ldownx = 1 ; Lupx = 2
    ENDIF
    XMAX = XGRID(NX)
  ELSE
    XMAX = XGRID(NX+1)
  ENDIF
  IF (BTEST(BCFLAG,3)) THEN
    IF (SIN(SOLARAZ) .GT. 1.0E-5) THEN
      Ldowny = 4 ; Lupy = 3
    ELSE
      Ldowny = 3 ; Lupy = 4
    ENDIF
    YMAX = YGRID(NY)
  ELSE
    YMAX = YGRID(NY+1)
  ENDIF

   ! Do the initial direct beam tracing starting at the base grid points
  NINVALID(:)=0
  DO IP = SPT, NPTS
    DIRPATH = 0.0
    CALL DIRECT_BEAM_PROP (0, GRIDPOS(1,IP), GRIDPOS(2,IP), GRIDPOS(3,IP), &
                 BCFLAG, IPFLAG, DELTAM,ML,NLEG, &
                 SOLARFLUX,SOLARMU,SOLARAZ,  DIRFLUX(IP), &
                 UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
     ! If we don't have a good direct beam then save the information 
     !   needed to pass to two upstream processors
    IF (.NOT. VALIDBEAM) THEN
      IF (SIDE /= Lupx .AND. SIDE /= Lupy) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Bad SIDE 1')
      NINVALID(SIDE) = NINVALID(SIDE) + 1
      IF (NINVALID(SIDE) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: NINVALID(SIDE) exceeds MAXNBEAM 1')
      BndDirRealInfo(:,NINVALID(SIDE),SIDE) = (/ XO, YO, ZO, DIRPATH /)
      BndDirIntInfo(:,NINVALID(SIDE),SIDE) = (/ myproc, IP /)
    ENDIF
  ENDDO


   ! Cycle over the passing of info between processors and doing more
   !  partial direct beam integrations until all rays have been finished
  Ndone(:) = 0
  AllDone = .FALSE.
  DO WHILE (.NOT. AllDone)
     ! Do the send and receives of the number of incomplete rays and
     !   the real and integer info needed about the radiance integrations
    irq=0
    IF (BTEST(BCFLAG,2)) THEN
      CALL MPI_ISEND (NINVALID(Lupx), 1, MPI_INTEGER, &
                      iprocneigh(Lupx), 1, comm2d, requests(irq+1), ierr)
      CALL MPI_ISEND (BndDirRealInfo(:,:,Lupx), 4*NINVALID(Lupx), MPI_REAL, &
                      iprocneigh(Lupx), 2, comm2d, requests(irq+2), ierr)
      CALL MPI_ISEND (BndDirIntInfo(:,:,Lupx), 2*NINVALID(Lupx), MPI_INTEGER, &
                      iprocneigh(Lupx), 3, comm2d, requests(irq+3), ierr)
      irq=irq+3
    ENDIF
    IF (BTEST(BCFLAG,3)) THEN
      CALL MPI_ISEND (NINVALID(Lupy), 1, MPI_INTEGER, &
                      iprocneigh(Lupy), 4, comm2d, requests(irq+1), ierr)
      CALL MPI_ISEND (BndDirRealInfo(:,:,Lupy), 4*NINVALID(Lupy), MPI_REAL, &
                      iprocneigh(Lupy), 5, comm2d, requests(irq+2), ierr)
      CALL MPI_ISEND (BndDirIntInfo(:,:,Lupy), 2*NINVALID(Lupy), MPI_INTEGER, &
                      iprocneigh(Lupy), 6, comm2d, requests(irq+3), ierr)
      irq=irq+3
    ENDIF
    IF (BTEST(BCFLAG,2)) THEN
      CALL MPI_IRECV (NINVALID(Ldownx), 1, MPI_INTEGER, &
                      iprocneigh(Ldownx), 1, comm2d, requests(irq+1), ierr)
      CALL MPI_IRECV (BndDirRealInfo(:,:,Ldownx), 4*MAXNBEAM, MPI_REAL, &
                      iprocneigh(Ldownx), 2, comm2d, requests(irq+2), ierr)
      CALL MPI_IRECV (BndDirIntInfo(:,:,Ldownx), 2*MAXNBEAM, MPI_INTEGER, &
                      iprocneigh(Ldownx), 3, comm2d, requests(irq+3), ierr)
      irq=irq+3
    ENDIF
    IF (BTEST(BCFLAG,3)) THEN
      CALL MPI_IRECV (NINVALID(Ldowny), 1, MPI_INTEGER, &
                      iprocneigh(Ldowny), 4, comm2d, requests(irq+1), ierr)
      CALL MPI_IRECV (BndDirRealInfo(:,:,Ldowny), 4*MAXNBEAM, MPI_REAL, &
                      iprocneigh(Ldowny), 5, comm2d, requests(irq+2), ierr)
      CALL MPI_IRECV (BndDirIntInfo(:,:,Ldowny), 2*MAXNBEAM, MPI_INTEGER, &
                      iprocneigh(Ldowny), 6, comm2d, requests(irq+3), ierr)
      irq=irq+3
    ENDIF
    CALL MPI_WAITALL (irq, requests, status, ierr)

    NINVALID(Lupx)=0 ; NINVALID(Lupy)=0
    IF (BTEST(BCFLAG,2)) THEN
       ! Continue the backwards ray integrations from the X boundary
      DO J = 1, NINVALID(Ldownx)
        IF (Ldownx == 1) BndDirRealInfo(1,J,Ldownx) = XGRID(1)
        IF (Ldownx == 2) BndDirRealInfo(1,J,Ldownx) = XGRID(NX) 
        X=BndDirRealInfo(1,J,Ldownx) ; Y=BndDirRealInfo(2,J,Ldownx)
        Z=BndDirRealInfo(3,J,Ldownx) ; DIRPATH=BndDirRealInfo(4,J,Ldownx)
        IF (Y < YGRID(1)-EPS) Y = YMAX ; IF (Y > YMAX+EPS) Y = YGRID(1)
        CALL DIRECT_BEAM_PROP (0, X, Y, Z, BCFLAG, IPFLAG, DELTAM,ML,NLEG, &
                               SOLARFLUX,SOLARMU,SOLARAZ,  DIRJUNK, &
                               UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
         ! If we got a valid direct beam then store it in DirAllOut otherwise
         ! save the information to continue passing to neighboring processors
        IF (VALIDBEAM) THEN
          iproc = BndDirIntInfo(1,J,Ldownx)
          Ndone(iproc) = Ndone(iproc) + 1
          IF (Ndone(iproc) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Ndone(iproc) exceeds MAXNBEAM x')
          DirAllOut(Ndone(iproc),iproc) = DIRPATH
          DirPtrOut(Ndone(iproc),iproc) = BndDirIntInfo(2,J,Ldownx)
        ELSE
          IF (SIDE /= Lupx .AND. SIDE /= Lupy) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Bad SIDE 2x')
          NINVALID(SIDE) = NINVALID(SIDE) + 1
          IF (NINVALID(SIDE) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: NINVALID(SIDE) exceeds MAXNBEAM 2x')
          BndDirRealInfo(:,NINVALID(SIDE),SIDE) = (/ XO, YO, ZO, DIRPATH /)
          BndDirIntInfo(:,NINVALID(SIDE),SIDE) = BndDirIntInfo(:,J,Ldownx)
        ENDIF
      ENDDO
    ENDIF

    IF (BTEST(BCFLAG,3)) THEN
       ! Continue the backwards ray integrations for the Y boundary
      DO J = 1, NINVALID(Ldowny)
        IF (Ldowny == 3) BndDirRealInfo(2,J,Ldowny) = YGRID(1)
        IF (Ldowny == 4) BndDirRealInfo(2,J,Ldowny) = YGRID(NY) 
        X=BndDirRealInfo(1,J,Ldowny) ; Y=BndDirRealInfo(2,J,Ldowny)
        Z=BndDirRealInfo(3,J,Ldowny) ; DIRPATH=BndDirRealInfo(4,J,Ldowny)
        IF (X < XGRID(1)-EPS) X = XMAX ; IF (X > XMAX+EPS) X = XGRID(1)
        CALL DIRECT_BEAM_PROP (0, X, Y, Z, BCFLAG, IPFLAG, DELTAM,ML,NLEG, &
                               SOLARFLUX,SOLARMU,SOLARAZ,  DIRJUNK, &
                               UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
         ! If we got a valid direct beam then store it in DirAllOut otherwise
         ! save the information to continue passing to neighboring processors
        IF (VALIDBEAM) THEN
          iproc = BndDirIntInfo(1,J,Ldowny)
          Ndone(iproc) = Ndone(iproc) + 1
          IF (Ndone(iproc) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Ndone(iproc) exceeds MAXNBEAM y')
          DirAllOut(Ndone(iproc),iproc) = DIRPATH
          DirPtrOut(Ndone(iproc),iproc) = BndDirIntInfo(2,J,Ldowny)
        ELSE
          IF (SIDE /= Lupx .AND. SIDE /= Lupy) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: Bad SIDE 2y')
          NINVALID(SIDE) = NINVALID(SIDE) + 1
          IF (NINVALID(SIDE) > MAXNBEAM) CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: NINVALID(SIDE) exceeds MAXNBEAM 2x')
          BndDirRealInfo(:,NINVALID(SIDE),SIDE) = (/ XO, YO, ZO, DIRPATH /)
          BndDirIntInfo(:,NINVALID(SIDE),SIDE) = BndDirIntInfo(:,J,Ldowny)
        ENDIF
      ENDDO
    ENDIF

     ! See if there are any more invalid radiance rays on all the processors
    intsendbuf = (/ NINVALID(Lupx), NINVALID(Lupy) /)
    CALL MPI_ALLREDUCE (intsendbuf, intrecvbuf, 2, MPI_INTEGER, MPI_MAX, &
                        comm2d, ierr)
    AllDone = ALL(intrecvbuf(1:2) == 0)
    IF (MAXVAL(intrecvbuf(1:2)) > MAXNBEAM .OR. MAXVAL(Ndone(:)) > MAXNBEAM) THEN
      CALL ABORT_SHDOM_MPI ('MAKE_DIRECT_PAR: MAXNBEAM exceeded 2')
    ENDIF
  ENDDO


   ! Exchange the finished direct beam paths with all of the other processors
  DO iproc = 0, numproc-1
    IF (iproc /= myproc) THEN
      CALL MPI_ISEND (Ndone(iproc), 1, MPI_INTEGER, &
                      iproc, 7, comm2d, requests(1), ierr)
      CALL MPI_ISEND (DirAllOut(:,iproc), Ndone(iproc), MPI_REAL, &
                      iproc, 8, comm2d, requests(2), ierr)
      CALL MPI_ISEND (DirPtrOut(:,iproc), Ndone(iproc), MPI_INTEGER, &
                      iproc, 9, comm2d, requests(3), ierr)
      CALL MPI_IRECV (Nbeam, 1, MPI_INTEGER, &
                      iproc, 7, comm2d, requests(4), ierr)
      CALL MPI_IRECV (DirAllIn, MAXNBEAM, MPI_REAL, &
                      iproc, 8, comm2d, requests(5), ierr)
      CALL MPI_IRECV (DirPtrIn, MAXNBEAM, MPI_INTEGER, &
                      iproc, 9, comm2d, requests(6), ierr)
      CALL MPI_WAITALL (6, requests, status, ierr)
       ! Put the direct beam flux in the correct grid point
      DO I = 1, Nbeam
        IP = DirPtrIn(I)
        DIRFLUX(IP) = SOLARFLUX*EXP(-DirAllIn(I))
      ENDDO
    ELSE
       ! If we have the direct paths for this processor than avoid MPI and
       !   put the direct beam flux in the grid point
      DO I = 1, Ndone(myproc)
        IP = DirPtrOut(I,myproc)
        DIRFLUX(IP) = SOLARFLUX*EXP(-DirAllOut(I,myproc))
      ENDDO
    ENDIF
  ENDDO

  DEALLOCATE (DirAllOut, DirPtrOut, DirAllIn, DirPtrIn, Ndone)
  DEALLOCATE (BndDirRealInfo, BndDirIntInfo)
END SUBROUTINE MAKE_DIRECT_PAR



SUBROUTINE FIND_BOUNDARY_POINTS (BCFLAG, IPFLAG, NPTS, SWEEPORD, GRIDPTR, &
                                 GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)
 ! Makes the pointers to and positions of the horizontal boundary grid 
 ! points of the processors neighboring this one. 
 ! For each of the potentially eight octants there are
 ! two sets of boundary pointers; 1 is for the current subdomain,
 ! and 2 is for the neighboring subdomain.
 ! Fills in the BXYPTR(:,IBC,JOCT) and BXYPOS(:,IBC,JOCT) arrays 
 ! (module variables), where IBC is a boundary point index and JOCT is
 ! a sort of octant number. BXYPTR holds (iproc,IPT,SIDE) for the boundary 
 ! grid points from grid point IPT in processor iproc that are on the SIDE
 ! boundary of this subdomain (SIDE: 1=-X, 2=+X, 3=-Y, 4=+Y).  BXYPOS holds 
 ! the gridpoint positions (X,Y,Z) for these boundary points.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BCFLAG, IPFLAG, NPTS, SWEEPORD(NPTS,*), GRIDPTR(8,*)
  INTEGER, INTENT(IN) :: NX, NY, NZ
  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS), XGRID(NX), YGRID(NY), ZGRID(NZ)
  include 'mpif.h'
  INTEGER :: tag, ierr, irq, requests(24), status(MPI_STATUS_SIZE,24)
  INTEGER :: MAXNBXY, NOCT, JOCT, IORDER, IPCELL, ICORNER, IPT, SIDE
  INTEGER :: L, LI, N(4), N1(8,4), N2(8,4), NS, NBX, NBY
  INTEGER :: MAXBC, IBC, BITZ, SZ, EZ, DZ, LX, LY, IX, IY, IZ, J
  INTEGER :: INOUT(4) = (/ 2, 1, 4, 3 /)
  INTEGER, SAVE :: OLDNPTS=0
  LOGICAL :: SAMEGRID
  REAL    :: XEDGE(8), YEDGE(8), XYB(4)
  INTEGER, ALLOCATABLE :: BNDPTR(:,:,:,:), INDX(:)
  REAL,    ALLOCATABLE :: BNDPOS(:,:,:,:,:), ZARRAY(:)

  IF (.NOT. (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3))) RETURN
  IF (BTEST(IPFLAG,0) .AND. BTEST(IPFLAG,1)) RETURN
  CALL  MPI_ALLREDUCE (NPTS==OLDNPTS, SAMEGRID, 1, MPI_LOGICAL, MPI_LAND, comm2d, ierr)
  OLDNPTS = NPTS
  IF (SAMEGRID) RETURN


!   IOCT=1+BITX+2*BITY+4*BITZ ; JOCT=JOCTORDER(IOCT) = (/ 1,3,5,7,2,4,6,8 /)
!      CX<0: BITX=1, JOCT=3,7,4,8      CY<0: BITY=1, JOCT=5,7,6,8
  XEDGE(:) = (/ XGRID(NX),XGRID(NX),XGRID(1),XGRID(1),XGRID(NX),XGRID(NX),XGRID(1),XGRID(1) /)
  YEDGE(:) = (/ YGRID(NY),YGRID(NY),YGRID(NY),YGRID(NY),YGRID(1),YGRID(1),YGRID(1),YGRID(1) /)
  XYB(1) = XGRID(1) ; XYB(2) = XGRID(NX)
  XYB(3) = YGRID(1) ; XYB(4) = YGRID(NY)

  IF (BTEST(IPFLAG,1) .AND. BTEST(IPFLAG,0)) THEN
    NOCT = 2
  ELSE IF (BTEST(IPFLAG,1)) THEN
    NOCT = 4
  ELSE
    NOCT = 8
  ENDIF
     ! Allocate grid point pointer and position arrays (second from last
     !   index is 1 for local, 2 for neighbor; last index is boundary L)
  IF (ALLOCATED(IBXYPTR)) THEN
    MAXNBXY = NINT(3.0*MAXVAL(IBXYPTR(2,:,:)))
  ELSE
    MAXNBXY = 3.0*(NZ-2)*MAX(NX,NY)
  ENDIF
  ALLOCATE (BNDPTR(NOCT,MAXNBXY,2,4), BNDPOS(3,NOCT,MAXNBXY,2,4))

  N1(:,:) = 0  
  DO JOCT = 1, NOCT
    BITZ = 1-IBITS(JOCT,0,1)
    N(:) = 0
    DO IORDER = 1, NPTS
      IPCELL = ISHFT(SWEEPORD(IORDER,JOCT),-3)
      ICORNER = IBITS(SWEEPORD(IORDER,JOCT),0,3)+1
      IPT = GRIDPTR(ICORNER,IPCELL)

      IF ( (GRIDPOS(3,IPT) < ZGRID(NZ) .AND. BITZ==0) .OR. &
           (GRIDPOS(3,IPT) > ZGRID(1)  .AND. BITZ==1) ) THEN
        IF (BTEST(BCFLAG,2) .AND. GRIDPOS(1,IPT)==XEDGE(JOCT) &
            .AND. .NOT. BTEST(IPFLAG,0)) THEN
          SIDE = 2-IBITS(JOCT-1,1,1)
          N(SIDE) = N(SIDE) + 1
          IF (N(SIDE)>MAXNBXY) then
            CALL ABORT_SHDOM_MPI ('FIND_BOUNDARY_POINTS: MAXNBXY exceeded by N(SIDE) 1x')
          ENDIF
          BNDPTR(JOCT,N(SIDE),1,SIDE) = IPT
          BNDPOS(:,JOCT,N(SIDE),1,SIDE) = GRIDPOS(:,IPT)
        ENDIF
        IF (BTEST(BCFLAG,3) .AND. GRIDPOS(2,IPT)==YEDGE(JOCT) &
            .AND. .NOT. BTEST(IPFLAG,1)) THEN
          SIDE = 4-IBITS(JOCT-1,2,1)
          N(SIDE) = N(SIDE) + 1
          IF (N(SIDE)>MAXNBXY) then
            CALL ABORT_SHDOM_MPI ('FIND_BOUNDARY_POINTS: MAXNBXY exceeded by N(SIDE) 1y')
          ENDIF
          BNDPTR(JOCT,N(SIDE),1,SIDE) = IPT
          BNDPOS(:,JOCT,N(SIDE),1,SIDE) = GRIDPOS(:,IPT)
        ENDIF
      ENDIF
    ENDDO
    N1(JOCT,:) = N(:)
  ENDDO


  NS = MAXVAL(N1(:,:))
  IF (NS > MAXNBXY) CALL ABORT_SHDOM_MPI ('FIND_BOUNDARY_POINTS: MAXNBXY might be exceeded')
  N2(:,:) = 0
  irq = 0
  DO L = 1, 4
    IF ((BTEST(BCFLAG,2) .AND. (L.EQ.1 .OR. L.EQ.2)) .OR. &
        (BTEST(BCFLAG,3) .AND. (L.EQ.3 .OR. L.EQ.4))) THEN
       ! Exchange the lists of boundary points with neighbors
      LI = INOUT(L)
      tag = L
      CALL MPI_ISEND (N1(:,L), 8, MPI_INTEGER, &
                      iprocneigh(L), tag, comm2d, requests(irq+1), ierr)
      CALL MPI_IRECV (N2(:,LI), 8, MPI_INTEGER, &
                      iprocneigh(LI), tag, comm2d, requests(irq+2), ierr)
      tag = 1000+L
      CALL MPI_ISEND (BNDPTR(:,:,1,L), NOCT*NS, MPI_INTEGER, &
                      iprocneigh(L), tag, comm2d, requests(irq+3), ierr)
      CALL MPI_IRECV (BNDPTR(:,:,2,LI), NOCT*MAXNBXY, MPI_INTEGER, &
                      iprocneigh(LI), tag, comm2d, requests(irq+4), ierr)
      tag = 2000+L
      CALL MPI_ISEND (BNDPOS(:,:,:,1,L), 3*NOCT*NS, MPI_REAL, &
                      iprocneigh(L), tag, comm2d, requests(irq+5), ierr)
      CALL MPI_IRECV (BNDPOS(:,:,:,2,LI), 3*NOCT*MAXNBXY, MPI_REAL, &
                      iprocneigh(LI), tag, comm2d, requests(irq+6), ierr)
      irq = irq + 6
    ENDIF
  ENDDO
   ! Wait for the transfers among all the neighbors to resolve
  CALL MPI_WAITALL (irq, requests, status, ierr)
  NS = MAXVAL(N2(:,:))
  IF (NS > MAXNBXY) CALL ABORT_SHDOM_MPI ('FIND_BOUNDARY_POINTS: MAXNBXY was exceeded')


  IF (ALLOCATED(BXYPTR)) DEALLOCATE (BXYPTR, BXYPOS, IBXYPTR)
  MAXBC = MAXVAL( (/ N2(:,1)+N2(:,3), N2(:,1)+N2(:,4), N2(:,2)+N2(:,3), N2(:,2)+N2(:,4) /) )
  ALLOCATE (BXYPTR(3,MAXBC,NOCT), BXYPOS(3,MAXBC,NOCT), IBXYPTR(2,NZ,NOCT))

    ! Put the boundary points from the neighboring subdomains in the
    ! BXYPTR and BXYPOS arrays.  Merge in X and Y boundary points one Z slab
    ! at a time, and set pointers (IBXYPTR) to start and end of each slab.
    !  DO IBC = IBXYPTR(1,IZ,JOCT), IBXYPTR(2,IZ,JOCT)
    !    X=BXYPOS(1,IBC,JOCT) ; Y=BXYPOS(2,IBC,JOCT) ; Z=BXYPOS(3,IBC,JOCT)
    !    iproc=BXYPTR(1,IBC,JOCT) ; IPT=BXYPTR(2,IBC,JOCT) ; SIDE=BXYPTR(2,IBC,JOCT)

  DO JOCT = 1, NOCT
    BITZ = 1-IBITS(JOCT,0,1)
    DZ = 2*BITZ-1
    LX = 1+IBITS(JOCT-1,1,1)
    LY = 3+IBITS(JOCT-1,2,1)
     ! Combine the indices to the X and Y boundary points and sort on Z
    NBX=0;  IF (BTEST(BCFLAG,2)) NBX=N2(JOCT,LX)
    NBY=0;  IF (BTEST(BCFLAG,3)) NBY=N2(JOCT,LY)
    ALLOCATE (ZARRAY(NBX+NBY), INDX(NBX+NBY))
    ZARRAY(1:NBX) = BNDPOS(3,JOCT,1:NBX,2,LX)
    INDX(1:NBX) = (/ (J, J=1,NBX) /)
    ZARRAY(NBX+1:NBX+NBY) = BNDPOS(3,JOCT,1:NBY,2,LY)
    INDX(NBX+1:NBX+NBY) = (/ (-J, J=1,NBY) /)
    CALL SSORT (ZARRAY, INDX, NBX+NBY, 2*DZ)

     ! For each Z slab use Z sorted index to transfer the X or Y boundary points
    IBXYPTR(:,:,JOCT) = 0
    SZ = (NZ-1)*(1-BITZ)+1+DZ
    EZ = (NZ-1)*(BITZ)+1
    IBC=0 ; J=1
    DO IZ = SZ, EZ, DZ
      IBXYPTR(1,IZ,JOCT) = IBC+1
       ! Store the boundary point info for this Z slab
      DO WHILE (J <= NBX+NBY .AND. &
               ( (ZARRAY(J)>=ZGRID(IZ) .AND. BITZ==0) .OR. &
                 (ZARRAY(J)<=ZGRID(IZ) .AND. BITZ==1) ))
        IBC = IBC + 1
        IF (INDX(J) > 0) THEN
          IX = INDX(J)
          BXYPTR(:,IBC,JOCT) = (/ iprocneigh(LX), BNDPTR(JOCT,IX,2,LX), LX /)
          BXYPOS(:,IBC,JOCT) = BNDPOS(:,JOCT,IX,2,LX)
          BXYPOS(INT((LX+1)/2),IBC,JOCT) = XYB(LX)
        ELSE
          IY = ABS(INDX(J))
          BXYPTR(:,IBC,JOCT) = (/ iprocneigh(LY), BNDPTR(JOCT,IY,2,LY), LY /)
          BXYPOS(:,IBC,JOCT) = BNDPOS(:,JOCT,IY,2,LY)
          BXYPOS(INT((LY+1)/2),IBC,JOCT) = XYB(LY)
        ENDIF
        J = J + 1
      ENDDO
      IBXYPTR(2,IZ,JOCT) = IBC
    ENDDO
    DEALLOCATE (ZARRAY, INDX)
  ENDDO

  DEALLOCATE (BNDPTR, BNDPOS)
END SUBROUTINE FIND_BOUNDARY_POINTS





SUBROUTINE CALC_BOUNDARY_RADIANCES (BCFLAG, IPFLAG, JOCT, IZ, &
                                   NX, NY, NZ, XGRID, YGRID, ZGRID, &
                                   NA, NPTS, NCELLS, GRIDPTR, &
                                   NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS, &
                                   MU, PHI, EXTINCT, SOURCE, &
                                   KANG, GRIDRAD)
 ! Calculates the discrete ordinate boundary radiances for the KANG discrete
 ! ordinate (going in direction MU/PHI) for the grid points in the IZ slab.
 ! FIND_BOUNDARY_POINTS has already been called to set up the BXYPTR and 
 ! BXYPOS arrays for the grid points on the boundaries, and the IBXYPTR
 ! array which has the starting and ending points for the IZ slab.  The
 ! grid points are done in the sweeping order, and so the octant (JOCT)
 ! is used to index into the IBXYPTR, BXYPTR, and BXYPOS arrays.  The
 ! grid points specified in BXYPTR/BXYPOS are on the boundaries of this
 ! processor, and need to be traced back along the discrete ordinate 
 ! first within this subdomain.  The radiative transfer integrations
 ! with the extinction (EXTINCT) and source function (SOURCE) are done
 ! with the INTEGRATE_RAY routine.  If the top or bottom boundary is
 ! hit when tracing back then the radiance is saved, otherwise the 
 ! required information (X,Y,Z, transmission, partial radiance, original
 ! processor, and starting grid point) are saved to be sent to the 
 ! appropriate neighboring processor.  These data are received from the
 ! neighboring processors and the radiative transfer integrations are
 ! continued, looping until all the rays have been finished.  Finally,
 ! the finished radiance and original grid points are exchanged with
 ! all the other processors.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BCFLAG, IPFLAG, JOCT, IZ
  INTEGER, INTENT(IN) :: NX, NY, NZ, NA, NPTS, NCELLS, KANG
  REAL,    INTENT(IN) :: XGRID(NX), YGRID(NY), ZGRID(NZ)
  INTEGER, INTENT(IN) :: GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS)
  INTEGER, INTENT(IN) :: TREEPTR(2,NCELLS)
  INTEGER*2, INTENT(IN) :: CELLFLAGS(NCELLS)
  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS)
  REAL,    INTENT(IN) :: MU, PHI
  REAL,    INTENT(IN) :: EXTINCT(NPTS), SOURCE(NA,NPTS)
  REAL,    INTENT(INOUT) :: GRIDRAD(NPTS)
  include 'mpif.h'
  INTEGER :: MAXNBND, Nbcpnts, Lupx, Lupy, Ldownx, Ldowny
  INTEGER :: SIDE, NINVALID(0:4)
  INTEGER :: IBC, IPT, I, J, L
  INTEGER :: ierr(12), irq, requests(12), status(MPI_STATUS_SIZE,12), iproc
  INTEGER :: intsendbuf(2), intrecvbuf(2)
  LOGICAL :: VALIDRAD, BTEST, AllDone
  REAL    :: X, Y, Z, RADIANCE, TRANSMIT, xo,yo,zo
  INTEGER, ALLOCATABLE :: Ndone(:)
  INTEGER, ALLOCATABLE :: BndRadIntInfo(:,:,:)
  REAL, ALLOCATABLE    :: BndRadRealInfo(:,:,:)
  INTEGER, ALLOCATABLE :: RadPtrOut(:,:), RadPtrIn(:)
  REAL, ALLOCATABLE    :: RadAllOut(:,:), RadAllIn(:)


  MAXNBND = numproc*MAXVAL(IBXYPTR(2,:,:)-IBXYPTR(1,:,:)+1)
  ALLOCATE (BndRadRealInfo(5,MAXNBND,4), BndRadIntInfo(2,MAXNBND,4))
  ALLOCATE (RadAllOut(MAXNBND,0:numproc-1), RadPtrOut(MAXNBND,0:numproc-1))
  ALLOCATE (RadAllIn(MAXNBND), RadPtrIn(MAXNBND), Ndone(0:numproc-1))

  Lupx=0 ; Lupy=0 ; Ldownx=0 ; Ldowny=0
  IF (BTEST(BCFLAG,2)) THEN
    IF (COS(PHI) .GT. 1.0E-5) THEN
      Ldownx = 2 ; Lupx = 1
    ELSE
      Ldownx = 1 ; Lupx = 2
    ENDIF
  ENDIF
  IF (BTEST(BCFLAG,3)) THEN
    IF (SIN(PHI) .GT. 1.0E-5) THEN
      Ldowny = 4 ; Lupy = 3
    ELSE
      Ldowny = 3 ; Lupy = 4
    ENDIF
  ENDIF


   ! Do the initial ray tracing 
  Ndone(:) = 0
  NINVALID(:)=0
   ! Loop over the downstream boundary points from the neighbor subdomains
  DO IBC = IBXYPTR(1,IZ,JOCT), IBXYPTR(2,IZ,JOCT)
     ! Call INTEGRATE_RAY for this 
    X=BXYPOS(1,IBC,JOCT) ; Y=BXYPOS(2,IBC,JOCT) ; Z=BXYPOS(3,IBC,JOCT)
    iproc=BXYPTR(1,IBC,JOCT) ; IPT=BXYPTR(2,IBC,JOCT) ; SIDE=BXYPTR(3,IBC,JOCT)
    RADIANCE = 0.0 ; TRANSMIT = 1.0
    xo=x; yo=y; zo=z
    CALL INTEGRATE_RAY (BCFLAG, IPFLAG, NX, NY, NZ, NA, &
                        XGRID, YGRID, ZGRID, NPTS, NCELLS, &
                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                        GRIDPOS, EXTINCT, SOURCE, KANG, GRIDRAD, &
                        MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
     ! If we got a radiance then store it otherwise save
     !  the information needed to pass to other processors
    IF (VALIDRAD) THEN
      Ndone(iproc) = Ndone(iproc) + 1
      RadAllOut(Ndone(iproc),iproc) = RADIANCE
      RadPtrOut(Ndone(iproc),iproc) = IPT
    ELSE
      IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
        WRITE (6,*) 'Bad SIDE 1:',myproc,mu,phi,SIDE,x,y,z,xo,yo,zo
        CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES')
      ENDIF
      NINVALID(SIDE) = NINVALID(SIDE) + 1
      BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ X, Y, Z, TRANSMIT, RADIANCE /)
      BndRadIntInfo(:,NINVALID(SIDE),SIDE) = (/ iproc, IPT /)
    ENDIF
  ENDDO


   ! Cycle over the passing of info between processors and doing more
   !  partial radiance integrations until all rays have been finished
  AllDone = .FALSE.
  DO WHILE (.NOT. AllDone)
     ! Do the send and receives of the number of incomplete rays and
     !   the real and integer info needed about the radiance integrations
    irq=0
    IF (BTEST(BCFLAG,2)) THEN
      CALL MPI_ISEND (NINVALID(Lupx), 1, MPI_INTEGER, &
                      iprocneigh(Lupx), 1, comm2d, requests(irq+1), ierr(irq+1))
      CALL MPI_ISEND (BndRadRealInfo(:,:,Lupx), 5*NINVALID(Lupx), MPI_REAL, &
                      iprocneigh(Lupx), 2, comm2d, requests(irq+2), ierr(irq+2))
      CALL MPI_ISEND (BndRadIntInfo(:,:,Lupx), 2*NINVALID(Lupx), MPI_INTEGER, &
                      iprocneigh(Lupx), 3, comm2d, requests(irq+3), ierr(irq+3))
      irq=irq+3
    ENDIF
    IF (BTEST(BCFLAG,3)) THEN
      CALL MPI_ISEND (NINVALID(Lupy), 1, MPI_INTEGER, &
                      iprocneigh(Lupy), 4, comm2d, requests(irq+1), ierr(irq+1))
      CALL MPI_ISEND (BndRadRealInfo(:,:,Lupy), 5*NINVALID(Lupy), MPI_REAL, &
                      iprocneigh(Lupy), 5, comm2d, requests(irq+2), ierr(irq+2))
      CALL MPI_ISEND (BndRadIntInfo(:,:,Lupy), 2*NINVALID(Lupy), MPI_INTEGER, &
                      iprocneigh(Lupy), 6, comm2d, requests(irq+3), ierr(irq+3))
      irq=irq+3
    ENDIF
    IF (BTEST(BCFLAG,2)) THEN
      CALL MPI_IRECV (NINVALID(Ldownx), 1, MPI_INTEGER, &
                      iprocneigh(Ldownx), 1, comm2d, requests(irq+1), ierr(irq+1))
      CALL MPI_IRECV (BndRadRealInfo(:,:,Ldownx), 5*MAXNBND, MPI_REAL, &
                      iprocneigh(Ldownx), 2, comm2d, requests(irq+2), ierr(irq+2))
      CALL MPI_IRECV (BndRadIntInfo(:,:,Ldownx), 2*MAXNBND, MPI_INTEGER, &
                      iprocneigh(Ldownx), 3, comm2d, requests(irq+3), ierr(irq+3))
      irq=irq+3
    ENDIF
    IF (BTEST(BCFLAG,3)) THEN
      CALL MPI_IRECV (NINVALID(Ldowny), 1, MPI_INTEGER, &
                      iprocneigh(Ldowny), 4, comm2d, requests(irq+1), ierr(irq+1))
      CALL MPI_IRECV (BndRadRealInfo(:,:,Ldowny), 5*MAXNBND, MPI_REAL, &
                      iprocneigh(Ldowny), 5, comm2d, requests(irq+2), ierr(irq+2))
      CALL MPI_IRECV (BndRadIntInfo(:,:,Ldowny), 2*MAXNBND, MPI_INTEGER, &
                      iprocneigh(Ldowny), 6, comm2d, requests(irq+3), ierr(irq+3))
      irq=irq+3
    ENDIF
    IF (ANY(ierr(:) /= MPI_SUCCESS)) THEN
      WRITE (6,*) 'MPI error: ',ierr(1:12)
      CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES: MPI error 1')
    ENDIF
    CALL MPI_WAITALL (irq, requests, status, ierr(1))
    IF (ierr(1) /= MPI_SUCCESS) THEN
      WRITE (6,*) 'MPI_WAITALL error: ',ierr(1)
      CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES: MPI error 2')
    ENDIF

    NINVALID(Lupx)=0 ; NINVALID(Lupy)=0
    IF (BTEST(BCFLAG,2)) THEN
       ! Continue the backwards ray integrations for the X boundary
      DO J = 1, NINVALID(Ldownx)
        IF (Ldownx == 1) BndRadRealInfo(1,J,Ldownx) = XGRID(1)
        IF (Ldownx == 2) BndRadRealInfo(1,J,Ldownx) = XGRID(NX) 
        X=BndRadRealInfo(1,J,Ldownx) ; Y=BndRadRealInfo(2,J,Ldownx)
        Z=BndRadRealInfo(3,J,Ldownx)  
        SIDE = Ldownx ; TRANSMIT=BndRadRealInfo(4,J,Ldownx)
        RADIANCE=BndRadRealInfo(5,J,Ldownx)
        CALL INTEGRATE_RAY (BCFLAG, IPFLAG, NX, NY, NZ, NA, &
                        XGRID, YGRID, ZGRID, NPTS, NCELLS, &
                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                        GRIDPOS, EXTINCT, SOURCE, KANG, GRIDRAD, &
                        MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

         ! If we got a radiance then store it in RadAllOut (if for another 
         !   processor) otherwise save the information needed to continue 
         !   passing to neighboring processors
        IF (VALIDRAD) THEN
          iproc = BndRadIntInfo(1,J,Ldownx)
          Ndone(iproc) = Ndone(iproc) + 1
          RadAllOut(Ndone(iproc),iproc) = RADIANCE
          RadPtrOut(Ndone(iproc),iproc) = BndRadIntInfo(2,J,Ldownx)
        ELSE
          IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
            WRITE (6,*) 'Bad SIDE 2x:',myproc,mu,phi,SIDE,x,y,z
            CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES')
          ENDIF
          NINVALID(SIDE) = NINVALID(SIDE) + 1
          BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ X, Y, Z, TRANSMIT, RADIANCE /)
          BndRadIntInfo(:,NINVALID(SIDE),SIDE) = BndRadIntInfo(:,J,Ldownx)
        ENDIF
      ENDDO
    ENDIF

    IF (BTEST(BCFLAG,3)) THEN
       ! Continue the backwards ray integrations for the Y boundary
      DO J = 1, NINVALID(Ldowny)
        IF (Ldowny == 3) BndRadRealInfo(2,J,Ldowny) = YGRID(1)
        IF (Ldowny == 4) BndRadRealInfo(2,J,Ldowny) = YGRID(NY) 
        X=BndRadRealInfo(1,J,Ldowny) ; Y=BndRadRealInfo(2,J,Ldowny)
        Z=BndRadRealInfo(3,J,Ldowny)  
        SIDE = Ldowny ; TRANSMIT=BndRadRealInfo(4,J,Ldowny)
        RADIANCE=BndRadRealInfo(5,J,Ldowny)
        CALL INTEGRATE_RAY (BCFLAG, IPFLAG, NX, NY, NZ, NA, &
                        XGRID, YGRID, ZGRID, NPTS, NCELLS, &
                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                        GRIDPOS, EXTINCT, SOURCE, KANG, GRIDRAD, &
                        MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
         ! If we got a radiance then store it in RadAllOut (if for another 
         !   processor) otherwise save the information needed to continue 
         !   passing to neighboring processors
        IF (VALIDRAD) THEN
          iproc = BndRadIntInfo(1,J,Ldowny)
          Ndone(iproc) = Ndone(iproc) + 1
          RadAllOut(Ndone(iproc),iproc) = RADIANCE
          RadPtrOut(Ndone(iproc),iproc) = BndRadIntInfo(2,J,Ldowny)
        ELSE
          IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
            WRITE (6,*) 'Bad SIDE 2y:',myproc,mu,phi,SIDE,x,y,z
            CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES')
          ENDIF
          NINVALID(SIDE) = NINVALID(SIDE) + 1
          BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ X, Y, Z, TRANSMIT, RADIANCE /)
          BndRadIntInfo(:,NINVALID(SIDE),SIDE) = BndRadIntInfo(:,J,Ldowny)
        ENDIF
      ENDDO
    ENDIF

     ! See if there are any more invalid radiance rays on all the processors
    intsendbuf = (/ NINVALID(Lupx), NINVALID(Lupy) /)
    CALL MPI_ALLREDUCE (intsendbuf, intrecvbuf, 2, MPI_INTEGER, MPI_MAX, &
                        comm2d, ierr(1))
    AllDone = ALL(intrecvbuf(1:2) == 0)
  ENDDO

   ! Exchange the finished boundary radiances with all of the other processors
  DO iproc = 0, numproc-1
    IF (iproc /= myproc) THEN
      CALL MPI_ISEND (Ndone(iproc), 1, MPI_INTEGER, &
                      iproc, 7, comm2d, requests(1), ierr(1))
      CALL MPI_ISEND (RadAllOut(:,iproc), Ndone(iproc), MPI_REAL, &
                      iproc, 8, comm2d, requests(2), ierr(2))
      CALL MPI_ISEND (RadPtrOut(:,iproc), Ndone(iproc), MPI_INTEGER, &
                      iproc, 9, comm2d, requests(3), ierr(3))
      CALL MPI_IRECV (Nbcpnts, 1, MPI_INTEGER, &
                      iproc, 7, comm2d, requests(4), ierr(4))
      CALL MPI_IRECV (RadAllIn, MAXNBND, MPI_REAL, &
                      iproc, 8, comm2d, requests(5), ierr(5))
      CALL MPI_IRECV (RadPtrIn, MAXNBND, MPI_INTEGER, &
                      iproc, 9, comm2d, requests(6), ierr(6))
      CALL MPI_WAITALL (6, requests, status, ierr(7))
      IF (ANY(ierr(1:7) /= MPI_SUCCESS)) THEN
        WRITE (6,*) 'MPI error: ',ierr(1:7)
        CALL ABORT_SHDOM_MPI ('CALC_BOUNDARY_RADIANCES: MPI error 3')
      ENDIF
      DO I = 1, Nbcpnts
        IPT = RadPtrIn(I)
        GRIDRAD(IPT) = RadAllIn(I)
      ENDDO
    ELSE
      DO I = 1, Ndone(myproc)
        IPT = RadPtrOut(I,myproc)
        GRIDRAD(IPT) = RadAllOut(I,myproc)
      ENDDO
    ENDIF
  ENDDO

  DEALLOCATE (RadAllOut, RadPtrOut, RadAllIn, RadPtrIn, Ndone)
  DEALLOCATE (BndRadRealInfo, BndRadIntInfo)
END SUBROUTINE CALC_BOUNDARY_RADIANCES




SUBROUTINE INTEGRATE_RAY (BCFLAG, IPFLAG, NX, NY, NZ, NA, &
                          XGRID, YGRID, ZGRID, NPTS, NCELLS,&
                          GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                          GRIDPOS, EXTINCT, SOURCE, KANG, GRIDRAD, MU, PHI, &
                          X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
 ! Integrates the radiative transfer equation tracing a ray starting 
 ! at X,Y,Z on boundary (face) SIDE to get the radiance in direction MU,PHI.
 ! The source function and extinction are integrated backwards from the 
 ! starting point.  As each cell is crossed the integration of the source 
 ! function across the cell is done and added to the radiance.  If there 
 ! are known radiances at the end of the ray path in GRIDRAD, then the 
 ! radiance there is interpolated from the four grid points on the face
 ! and VALIDRAD is returned true. Otherwise VALIDRAD is returned false.
 ! TRANSMIT is returned with the transmission from the ray starting point
 ! to the ending point, and RADIANCE has either the total radiance or the
 ! accumulated radiance until the boundary (TRANSMIT and RADIANCE must be
 ! properly initialized on input).  The ending location X,Y,Z and the
 ! subdomain SIDE (1=-X,2=+X,3=-Y,4=+Y) are also returned.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BCFLAG, IPFLAG, NX, NY, NZ, NA, NPTS, NCELLS, KANG
  REAL,    INTENT(IN) :: XGRID(NX), YGRID(NY), ZGRID(NZ)
  INTEGER, INTENT(IN) :: GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS)
  INTEGER, INTENT(IN) :: TREEPTR(2,NCELLS)
  INTEGER*2, INTENT(IN) ::  CELLFLAGS(NCELLS)
  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS)
  REAL,    INTENT(IN) :: EXTINCT(NPTS), SOURCE(NA,NPTS), GRIDRAD(NPTS)
  REAL,    INTENT(IN) :: MU, PHI
  REAL,    INTENT(INOUT) :: X, Y, Z
  INTEGER, INTENT(INOUT) :: SIDE
  REAL,    INTENT(INOUT) :: TRANSMIT, RADIANCE
  LOGICAL, INTENT(OUT) :: VALIDRAD

  INTEGER :: BITX, BITY, BITZ, IOCT, NXC, NYC, BTEST
  INTEGER :: IX, IY, IZ, IL, IM, IU, IPTR, DIR, ICELL, INEXTCELL
  INTEGER :: I1, I2, I3, I4, IOPP, IFACE, JFACE, KFACE, IC
  INTEGER :: GRIDFACE(4,6), OPPFACE(6)
  LOGICAL :: IPINX, IPINY
  DOUBLE PRECISION PI, CX, CY, CZ, CXINV, CYINV, CZINV
  DOUBLE PRECISION XE, YE, ZE
  DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
  DOUBLE PRECISION U, V, F1, F2, F3, F4
  DOUBLE PRECISION EXT, EXT0, EXT1, SRC, SRCEXT0, SRCEXT1
  DOUBLE PRECISION EXT0P, SRCEXT0P
  DOUBLE PRECISION TAU, TRANSCELL, ABSCELL, RAD, RAD0, TRANS
  DATA GRIDFACE/1,3,5,7, 2,4,6,8,  1,2,5,6, 3,4,7,8, 1,2,3,4, 5,6,7,8/
  DATA OPPFACE/2,1,4,3,6,5/


  EPS = 1.0E-3*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))

   ! Make the ray direction (opposite to the discrete ordinate direction)
  PI = ACOS(-1.0)
  CX = SQRT(1.0-MU**2)*COS(PHI+PI)
  CY = SQRT(1.0-MU**2)*SIN(PHI+PI)
  CZ = -MU
  IF (ABS(CX) .GT. 1.0E-5) THEN
    CXINV = 1.0D0/CX
  ELSE
    CX = 0.0
    CXINV = 1.0E6
  ENDIF
  IF (ABS(CY) .GT. 1.0E-5) THEN
    CYINV = 1.0D0/CY
  ELSE
    CY = 0.0
    CYINV = 1.0E6
  ENDIF
  CZINV = 1.0D0/CZ

   ! Setup the sweeping direction and the gridpoint location
   !   BITc is 0 for positive direction trace back ray, 1 for negative.
  IF (CX .LT. 0.0) THEN
    BITX = 1
  ELSE
    BITX = 0
  ENDIF
  IF (CY .LT. 0.0) THEN
    BITY = 1
  ELSE
    BITY = 0
  ENDIF
  IF (CZ .LT. -1.0E-3) THEN
    BITZ = 1
  ELSE IF (CZ .GT. 1.0E-3) THEN
    BITZ = 0
  ELSE 
    CALL ABORT_SHDOM_MPI ('INTEGRATE_RAY: Bad MU')
  ENDIF
   ! IOCT is the octant the ray come from or the discrete ordinate goes to.
  IOCT = 1 + BITX + 2*BITY + 4*BITZ

   ! Find the base grid cell that the starting point is in; first in Z
  IL=0
  IU=NZ
  DO WHILE (IU-IL .GT. 1)
    IM = (IU+IL)/2
    IF (Z .GE. ZGRID(IM)) THEN
      IL = IM
    ELSE
      IU=IM
    ENDIF
  ENDDO
  IZ = MAX(IL,1)

   ! Then get the base grid cell in X and Y assuming evenly spaced grid
  NXC = NX
  IF (BTEST(BCFLAG,2) .AND. .NOT. BTEST(IPFLAG,0)) NXC=NX-1
  NYC = NY
  IF (BTEST(BCFLAG,3) .AND. .NOT. BTEST(IPFLAG,1)) NYC=NY-1
  IF (BTEST(IPFLAG,0)) THEN
    IX = 1+NINT((NX-1)*(X-XGRID(1))/(XGRID(NX)-XGRID(1)))
  ELSE
    IX = 1+INT((NX-1)*(X-XGRID(1))/(XGRID(NX)-XGRID(1)))
  ENDIF
  IF (BTEST(IPFLAG,1)) THEN
    IY = 1+NINT((NY-1)*(Y-YGRID(1))/(YGRID(NY)-YGRID(1)))
  ELSE
    IY = 1+INT((NY-1)*(Y-YGRID(1))/(YGRID(NY)-YGRID(1)))
  ENDIF
  if (ix < 1 .or. iy < 1 .or. ix > NX+1 .or. iy > NY+1) then
    WRITE (6,*) 'Bad ix,iy',x,y,z,ix,iy,iz,xgrid(1),xgrid(nx),ygrid(1),ygrid(ny),myproc
    WRITE (6,*) 'Bad ix,iy',side,transmit
    CALL ABORT_SHDOM_MPI ('INTEGRATE_RAY:')
  endif
  IX = MIN(NXC,IX)
  IY = MIN(NYC,IY)


   ! Make the base grid cell pointer
  ICELL = IZ + (NZ-1)*(IY-1) + (NZ-1)*NYC*(IX-1)

   ! Trace down the tree: point to the positive child cell, get the grid
   ! point that is on the negative (X,Y,Z) side, use the flags to find 
   ! which dimension (X,Y,Z) the parent cell was split in, then compare 
   ! the grid point with the test point (X,Y,Z) to find which child cell 
   ! (pos or neg) the test point is in.
  DO WHILE (TREEPTR(2,ICELL) .GT. 0)
    DIR = IBITS(INT(CELLFLAGS(ICELL)),2,2)
    IC = TREEPTR(2,ICELL) + 1
    IPTR = GRIDPTR(1,IC)
    IF (DIR .EQ. 1) THEN
      IF (X .LT. GRIDPOS(1,IPTR))  IC = IC - 1
    ELSE IF (DIR .EQ. 2) THEN
      IF (Y .LT. GRIDPOS(2,IPTR))  IC = IC - 1
    ELSE IF (DIR .EQ. 3) THEN
      IF (Z .LT. GRIDPOS(3,IPTR))  IC = IC - 1
    ENDIF
    ICELL = IC
  ENDDO

  IFACE = SIDE
   ! Get the four gridpoints for this face
  I1 = GRIDPTR(GRIDFACE(1,IFACE),ICELL)
  I2 = GRIDPTR(GRIDFACE(2,IFACE),ICELL)
  I3 = GRIDPTR(GRIDFACE(3,IFACE),ICELL)
  I4 = GRIDPTR(GRIDFACE(4,IFACE),ICELL)
  IPINX = BTEST(INT(CELLFLAGS(ICELL)),0)
  IPINY = BTEST(INT(CELLFLAGS(ICELL)),1)
   ! Compute the face interpolation factors
  IF (IFACE .LE. 2) THEN
    U = (Z-GRIDPOS(3,I1))/(GRIDPOS(3,I3)-GRIDPOS(3,I1))
    IF (IPINY) THEN
      V = 0.5
    ELSE
      V = (Y-GRIDPOS(2,I1))/(GRIDPOS(2,I2)-GRIDPOS(2,I1))
    ENDIF
  ELSE
    U = (Z-GRIDPOS(3,I1))/(GRIDPOS(3,I3)-GRIDPOS(3,I1))
    IF (IPINX) THEN
      V = 0.5
    ELSE
      V = (X-GRIDPOS(1,I1))/(GRIDPOS(1,I2)-GRIDPOS(1,I1))
    ENDIF
  ENDIF

   ! Bilinearly interpolate the starting extinction and source function
  F1 = (1-U)*(1-V) ; F2 = (1-U)*V ; F3 = U*(1-V) ; F4 = U*V
  EXT1 = F1*EXTINCT(I1) + F2*EXTINCT(I2) + F3*EXTINCT(I3) + F4*EXTINCT(I4)
   ! Correctly interpolate source using extinction*source
  SRCEXT1 = F1*SOURCE(KANG,I1)*EXTINCT(I1) + F2*SOURCE(KANG,I2)*EXTINCT(I2) &
          + F3*SOURCE(KANG,I3)*EXTINCT(I3) + F4*SOURCE(KANG,I4)*EXTINCT(I4)
  EXT1 = MAX(0.0D0,EXT1)
  SRCEXT1 = MAX(0.0D0,SRCEXT1)

   ! Loop until finding a face with known radiances or reaching the boundary
  XE = X  ; YE = Y  ; ZE = Z
  TRANS = TRANSMIT
  RAD = RADIANCE
  VALIDRAD = .FALSE.
  DO WHILE (.NOT. VALIDRAD .AND. ICELL .GT. 0)
    IPINX = BTEST(INT(CELLFLAGS(ICELL)),0)
    IPINY = BTEST(INT(CELLFLAGS(ICELL)),1)

     ! Find boundaries of the current cell
     !   Find the three possible intersection planes (X,Y,Z)
     !   from the coordinates of the opposite corner grid point
    IOPP = GRIDPTR(9-IOCT,ICELL)
     ! Get the distances to the 3 planes and select the closest
    IF (IPINX) THEN
      SOX = 1.0E20
    ELSE
      SOX = (GRIDPOS(1,IOPP)-XE)*CXINV
    ENDIF
    IF (IPINY) THEN
      SOY = 1.0E20
    ELSE
      SOY = (GRIDPOS(2,IOPP)-YE)*CYINV
    ENDIF
    SOZ = (GRIDPOS(3,IOPP)-ZE)*CZINV
    SO = MIN(SOX,SOY,SOZ)
    IF (SO .LT. -EPS) THEN
      WRITE (6,*) 'INTEGRATE_RAY: SO<0  ',MU,PHI,X,Y,Z,SOX,SOY,SOZ,SO,ICELL,TRANS
      CALL ABORT_SHDOM_MPI ('INTEGRATE_RAY: SO<0')
    ENDIF
     ! Compute the coordinates of the cell exitting location 
    IF (.NOT. IPINX) XE = XE + SO*CX
    IF (.NOT. IPINY) YE = YE + SO*CY
    ZE = ZE + SO*CZ

     ! Get the intersection face number (i.e. neighptr index)
    IF (SOX .LE. SOZ .AND. SOX .LE. SOY) THEN
      IFACE = 2-BITX
      JFACE = 1
    ELSE IF (SOY .LE. SOZ) THEN
      IFACE = 4-BITY
      JFACE = 2
    ELSE
      IFACE = 6-BITZ
      JFACE = 3
    ENDIF
     ! Get the next cell to go to
    INEXTCELL = NEIGHPTR(IFACE,ICELL)
    IF (INEXTCELL .LT. 0) THEN
      CALL NEXT_CELL (XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS, &
                      GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
    ENDIF

     ! Get the face grid pointers
    IF (NEIGHPTR(IFACE,ICELL) .GE. 0) THEN
       ! If going to same or larger face then use previous face
      KFACE = IFACE
      IC = ICELL
    ELSE
       ! If going to smaller face then use next face (more accurate)
      KFACE = OPPFACE(IFACE)
      IC = INEXTCELL
    ENDIF
    I1 = GRIDPTR(GRIDFACE(1,KFACE),IC)
    I2 = GRIDPTR(GRIDFACE(2,KFACE),IC)
    I3 = GRIDPTR(GRIDFACE(3,KFACE),IC)
    I4 = GRIDPTR(GRIDFACE(4,KFACE),IC)
     ! Compute the face interpolation factors
    IF (JFACE .EQ. 1) THEN
      U = (ZE-GRIDPOS(3,I1))/(GRIDPOS(3,I3)-GRIDPOS(3,I1))
      IF (IPINY) THEN
        V = 0.5
      ELSE
        V = (YE-GRIDPOS(2,I1))/(GRIDPOS(2,I2)-GRIDPOS(2,I1))
      ENDIF
    ELSE IF (JFACE .EQ. 2) THEN
      U = (ZE-GRIDPOS(3,I1))/(GRIDPOS(3,I3)-GRIDPOS(3,I1))
      IF (IPINX) THEN
        V = 0.5
      ELSE
        V = (XE-GRIDPOS(1,I1))/(GRIDPOS(1,I2)-GRIDPOS(1,I1))
      ENDIF
    ELSE
      IF (IPINY) THEN
        U = 0.5
      ELSE
        U = (YE-GRIDPOS(2,I1))/(GRIDPOS(2,I3)-GRIDPOS(2,I1))
      ENDIF
      IF (IPINX) THEN
        V = 0.5
      ELSE
        V = (XE-GRIDPOS(1,I1))/(GRIDPOS(1,I2)-GRIDPOS(1,I1))
      ENDIF
    ENDIF
    if (u < -1.0E-4 .or. v < -1.0E-4 .or. u > 1.0001 .or. v > 1.0001) then
      WRITE (6,*) 'u,v<0 or u,v>1: ',mu,phi,x,y,z,xe,ye,ze,u,v,ipinx,ipiny
      CALL ABORT_SHDOM_MPI ('INTEGRATE_RAY:')
    endif

     ! Get the location coordinate (does the boundary wrapping)
    IF (INEXTCELL .GT. 0) THEN
      IF (JFACE .EQ. 1) THEN
        XE = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
      ELSE IF (JFACE .EQ. 2) THEN
        YE = GRIDPOS(2,GRIDPTR(IOCT,INEXTCELL))
      ELSE
        ZE = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
      ENDIF
    ENDIF

     ! Interpolate extinction and source function at face intersection
    F1 = (1-U)*(1-V) ; F2 = (1-U)*V ; F3 = U*(1-V) ; F4 = U*V
    EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2) + F3*EXTINCT(I3) + F4*EXTINCT(I4)
     ! Correctly interpolate source using extinction*source
    SRCEXT0 = F1*SOURCE(KANG,I1)*EXTINCT(I1) + F2*SOURCE(KANG,I2)*EXTINCT(I2) &
            + F3*SOURCE(KANG,I3)*EXTINCT(I3) + F4*SOURCE(KANG,I4)*EXTINCT(I4)
    EXT0 = MAX(0.0D0,EXT0)
    SRCEXT0 = MAX(0.0D0,SRCEXT0)
     ! Compute the cell radiance: integration of the source function
    EXT = 0.5*(EXT0+EXT1)
    TAU=EXT*SO
    IF (TAU .GE. 0.5) THEN
      TRANSCELL = EXP(-TAU)
      ABSCELL = 1.0 - TRANSCELL
    ELSE
      ABSCELL = TAU*(1.0-0.5*TAU*(1.0-0.33333333333*TAU*(1-0.25*TAU)))
      TRANSCELL = 1.0 - ABSCELL
    ENDIF
    IF (TAU .LE. 2.0) THEN
      IF (EXT .EQ. 0.0) THEN
        SRC = 0.0
      ELSE
         ! Linear extinction, linear source*extinction, to first order
        SRC = ( 0.5*(SRCEXT0+SRCEXT1) &
                + 0.08333333333*(EXT0*SRCEXT1-EXT1*SRCEXT0)*SO )/EXT
      ENDIF
    ELSE
       ! Combined first order expansion and constant extinction formula
      EXT0P = EXT0
      SRCEXT0P = SRCEXT0
      IF (TAU .GT. 4.0) THEN 
        EXT0P = EXT1 + (EXT0-EXT1)*4.0/TAU
        IF (EXT0 .GT. 0.0) SRCEXT0P = SRCEXT0*EXT0P/EXT0
      ENDIF
      SRC = 1.0/(EXT0P+EXT1) *( SRCEXT0P+SRCEXT1 &
             + (EXT0P*SRCEXT1-EXT1*SRCEXT0P)*2.0/(EXT0P+EXT1) &
                 *(1-2/TAU+2*TRANSCELL/ABSCELL) )
    ENDIF
    SRC = MAX(SRC,0.0D0)

     ! Add in the cell radiance and update the transmission.
     !  If this is a boundary or the transmission is below the
     !  cutoff and we have a valid radiance then set the flag
     !  to stop the tracing and add in the interpolated face radiance.
    RAD = RAD + TRANS*SRC*ABSCELL
    TRANS = TRANS*TRANSCELL
    if (RAD < -1.0E-5) then
      WRITE (6,'(3(1X,F5.3),8(1X,E12.5))'), &
        xe,ye,ze,ext0,ext1,so,tau,trans,src,abscell,rad
    endif
    IF (GRIDRAD(I1).GE.0.0 .AND. GRIDRAD(I2).GE.0.0 &
        .AND. GRIDRAD(I3).GE.0.0 .AND. GRIDRAD(I4).GE.0.0) THEN
      VALIDRAD = .TRUE.
      RAD0 = F1*GRIDRAD(I1) + F2*GRIDRAD(I2) + F3*GRIDRAD(I3) + F4*GRIDRAD(I4)
      RAD = RAD + TRANS*RAD0
    ELSE IF (TRANS .LT. 1.0E-5) THEN
      VALIDRAD = .TRUE.
    ELSE
      EXT1 = EXT0
      SRCEXT1 = SRCEXT0
      ICELL = INEXTCELL
    ENDIF
  ENDDO
  X = XE ; Y = YE ; Z = ZE
  RADIANCE = RAD
  TRANSMIT = TRANS
  SIDE = IFACE
END SUBROUTINE INTEGRATE_RAY





SUBROUTINE COMPUTE_RADIANCE_PAR (NX, NY, NZ, NPTS, NCELLS, &
                  ML, MM, NCS, NLEG, NUMPHASE, &
                  NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                  BCFLAG, XDOMAIN, YDOMAIN, IPFLAG, &
                  SRCTYPE, DELTAM, SOLARMU, SOLARAZ, &
                  SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
                  MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
                  GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS, &
                  XGRID, YGRID, ZGRID, GRIDPOS, &
                  GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                  EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, &
                  SHPTR, SOURCE, SOURCE1, GRIDRAD, &
                  OUTPARMS,  NRAD, RADOUT)
 ! Computes the radiances for the locations and directions specified 
 ! in OUTPARMS, putting the results in RADOUT and increasing NRAD by
 ! the number of radiances output.  Does the radiative transfer integrations
 ! across multiple processors.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
  INTEGER, INTENT(IN) :: ML, MM, NCS, NLEG, NUMPHASE
  INTEGER, INTENT(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
  INTEGER, INTENT(INOUT) :: NRAD
!f2py intent(in, out) :: NRAD
  INTEGER, INTENT(IN) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
  INTEGER, INTENT(IN) :: GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
  INTEGER, INTENT(IN) :: SHPTR(NPTS+1), BCPTR(MAXNBC,2)
  INTEGER*2, INTENT(IN) :: CELLFLAGS(NCELLS), IPHASE(NPTS)
  LOGICAL, INTENT(IN) :: DELTAM
  REAL, INTENT(IN) :: SOLARMU, SOLARAZ
  REAL, INTENT(IN) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
  REAL, INTENT(IN) :: MU(NMU), PHI(NMU,NPHI0MAX), WTDO(NMU,NPHI0MAX)
  REAL, INTENT(IN) :: XDOMAIN, YDOMAIN, XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)
  REAL, INTENT(IN) :: GRIDPOS(3,NPTS)
  REAL, INTENT(IN) :: SFCGRIDPARMS(*), BCRAD(*)
  REAL, INTENT(IN) :: EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(0:NLEG,NPTS)
  REAL, INTENT(IN) :: DIRFLUX(NPTS), FLUXES(2,NPTS), SOURCE(*)
  REAL, INTENT(INOUT) :: SOURCE1(NPTS), GRIDRAD(NPTS)
!f2py intent(in, out) :: SOURCE1, GRIDRAD
  REAL, INTENT(IN) :: OUTPARMS(*)
  REAL, INTENT(OUT) :: RADOUT(*)
!f2py intent(in, out) :: RADOUT
  CHARACTER, INTENT(IN) ::  SRCTYPE*1, SFCTYPE*2, UNITS*1

  INCLUDE "mpif.h" 
  INTEGER :: IBC, JX, JY, K
  INTEGER :: NANGOUT, NXOUT, NYOUT
  LOGICAL :: LAMBERTIAN, VALIDRAD
  DOUBLE PRECISION :: X0,Y0,Z0, XE,YE,ZE, X,Y,Z, TRANSMIT, RADIANCE
  REAL    :: MUOUT, PHIOUT, PHID
  REAL    :: STARTX, STARTY, EPS=1.0E-5

  INTEGER :: MAXNRAD, NDOWN, Lupx, Lupy, Ldownx, Ldowny
  INTEGER :: SIDE, NINVALID(0:4)
  INTEGER :: I, J, L, Ndone, intsendbuf(2), intrecvbuf(2)
  INTEGER :: ierr, irq, iproc
  LOGICAL :: AllDone
  INTEGER, ALLOCATABLE :: requests(:), status(:,:), NradRecv(:)
  INTEGER, ALLOCATABLE :: IndexSend(:), IndexRecv(:,:), BndRadIntInfo(:,:)
  REAL, ALLOCATABLE    :: RadSend(:), RadRecv(:,:), BndRadRealInfo(:,:,:)


   ! Make the isotropic radiances for the top boundary
  CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN, &
                              UNITS, NTOPPTS, BCRAD(1))
   ! Make the bottom boundary radiances for the Lambertian surfaces.  
   ! Compute the upwelling bottom radiances using the downwelling fluxes.
  IF (SFCTYPE .EQ. 'FL') THEN
    CALL FIXED_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2), &
                  DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO, &
                  WAVENO, WAVELEN, UNITS,   BCRAD(1+NTOPPTS))
  ELSE IF (SFCTYPE .EQ. 'VL') THEN
    CALL VARIABLE_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2), &
                    DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS, &
                    BCRAD(1+NTOPPTS))
  ENDIF


   ! Setup the regularly spaced radiance locations
  STARTX = 0.0
  IF (OUTPARMS(2) .LE. 0.0) THEN
    NXOUT = 1
  ELSE
    NXOUT = MAX(1,NINT(XDOMAIN/OUTPARMS(2)))
  ENDIF
  STARTY = 0.0
  IF (OUTPARMS(3) .LE. 0.0) THEN
    NYOUT = 1
  ELSE
    NYOUT = MAX(1,NINT(YDOMAIN/OUTPARMS(3)))
  ENDIF
  Z0 = MIN( MAX(OUTPARMS(1),ZGRID(1)), ZGRID(NZ))
  NANGOUT = NINT(OUTPARMS(6))


  LAMBERTIAN = SFCTYPE(2:2) .EQ. 'L'


  MAXNRAD = NINT(4.0*(NXOUT*NYOUT)/numproc)
!  MAXNRAD = NXOUT*NYOUT
  ALLOCATE (RadSend(MAXNRAD), IndexSend(MAXNRAD))
  IF (myproc==0) ALLOCATE (NradRecv(numproc-1), &
                    RadRecv(MAXNRAD,numproc-1), IndexRecv(MAXNRAD,numproc-1))
  ALLOCATE (BndRadRealInfo(5,MAXNRAD,4), BndRadIntInfo(MAXNRAD,4))
  ALLOCATE (requests(MAX(12,3*numproc)), status(MPI_STATUS_SIZE,MAX(12,3*numproc)))


   ! Loop over the radiance directions
  DO K = 1, NANGOUT
    MUOUT = OUTPARMS(2*K+5)
    PHID = OUTPARMS(2*K+6)
    PHIOUT = PHID*ACOS(-1.0)/180.0
    IF (MUOUT .EQ. 0.0 .OR. ABS(MUOUT) .GT. 1.0) THEN
      WRITE (6,*) 'COMPUTE_RADIANCE: Bad mu for radiance',MUOUT
    ELSE

       ! Compute the source function throughout grid for this angle
      CALL COMPUTE_ONE_SOURCE (ML, MM, NCS, NLEG, NUMPHASE, &
                NPTS, DELTAM, MUOUT, PHIOUT, &
                SRCTYPE, SOLARMU, SOLARAZ, ALBEDO, LEGEN, IPHASE, &
                DIRFLUX, SHPTR, SOURCE,  SOURCE1)

       ! Set the radiance field to -1, so we can determine valid radiances
       !   Also set the source array to the extinction times the source
       !   function, so the grid interpolation is correct.
      DO I = 1, NPTS
        GRIDRAD(I) = -1.0
        SOURCE1(I) = SOURCE1(I)*EXTINCT(I)
      ENDDO
       ! Get boundary radiances: either top or bottom
       ! Isotropic top boundary or Lambertian bottom boundary can use
       !   the previously computed boundary radiances in BCRAD,
       !   otherwise, compute the radiance for this angle
       !   by integrating over the stored downwelling radiances.
      IF (MUOUT .LT. 0.0) THEN
        DO IBC = 1, NTOPPTS
          I = BCPTR(IBC,1)
          GRIDRAD(I) = BCRAD(IBC)
        ENDDO
      ELSE
        IF (.NOT. LAMBERTIAN) THEN
          CALL VARIABLE_BRDF_SURFACE (NBOTPTS,1,NBOTPTS,BCPTR(1,2), &
                  NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MUOUT, PHIOUT, &
                  SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX, &
                  SFCTYPE, NSFCPAR, SFCGRIDPARMS, BCRAD(1+NTOPPTS))
        ENDIF
        DO IBC = 1, NBOTPTS
          I = BCPTR(IBC,2)
          GRIDRAD(I) = BCRAD(NTOPPTS+IBC)
        ENDDO
      ENDIF

       ! Set up the upstream/downstream boundaries for this direction
      Lupx=0 ; Lupy=0 ; Ldownx=0 ; Ldowny=0
      IF (BTEST(BCFLAG,2)) THEN
        IF (COS(PHIOUT) .GT. 1.0E-5) THEN
          Ldownx = 2 ; Lupx = 1
        ELSE
          Ldownx = 1 ; Lupx = 2
        ENDIF
      ENDIF
      IF (BTEST(BCFLAG,3)) THEN
        IF (SIN(PHIOUT) .GT. 1.0E-5) THEN
          Ldowny = 4 ; Lupy = 3
        ELSE
          Ldowny = 3 ; Lupy = 4
        ENDIF
      ENDIF

      NINVALID(:)=0
      Ndone = 0

       ! Loop over all the radiance starting locations and integrate 
       !   backwards for those in this subdomain.
      Y0 = STARTY
      DO JY = 1, NYOUT
        X0 = STARTX
        DO JX = 1, NXOUT
          NRAD = NRAD + 1
          IF (X0 .GE. XGRID(1)-EPS .AND. X0 .LE. XGRID(NX)+EPS .AND. &
              Y0 .GE. YGRID(1)-EPS .AND. Y0 .LE. YGRID(NY)+EPS) THEN
            TRANSMIT = 1.0D0 ; RADIANCE = 0.0D0
            CALL INTEGRATE_SOURCE (BCFLAG, IPFLAG, &
                             NX, NY, NZ, NPTS, NCELLS, &
                             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                             XGRID, YGRID, ZGRID, GRIDPOS, &
                             MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1, &
                             X0, Y0, Z0, XE,YE,ZE, SIDE, &
                             TRANSMIT, RADIANCE, VALIDRAD)

             ! If we got a radiance then store it in RadSend otherwise save
             !  the information needed to pass to other processors
            IF (VALIDRAD) THEN
              Ndone = Ndone + 1
              IF (Ndone>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 1')
              RadSend(Ndone) = RADIANCE
              IndexSend(Ndone) = NRAD
            ELSE
              IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
                WRITE (6,*) 'Bad SIDE 1:',myproc,muout,phiout,x0,y0,z0,SIDE,xe,ye,ze
                CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR')
              ENDIF
              NINVALID(SIDE) = NINVALID(SIDE) + 1
              IF (NINVALID(SIDE)>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 1')
              BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ XE, YE, ZE, TRANSMIT, RADIANCE /)
              BndRadIntInfo(NINVALID(SIDE),SIDE) = NRAD
            ENDIF
          ENDIF
          X0 = X0 + OUTPARMS(2)
        ENDDO
        Y0 = Y0 + OUTPARMS(3)
      ENDDO

       ! Cycle over the passing of info between processors and doing more
       !  partial radiance integrations until all rays have been finished
      AllDone = .FALSE.
      DO WHILE (.NOT. AllDone)
         ! Do the send and receives of the number of incomplete rays and
         !   the real and integer info needed about the radiance integrations
        irq=0
        IF (BTEST(BCFLAG,2)) THEN
          CALL MPI_ISEND (NINVALID(Lupx), 1, MPI_INTEGER, &
                          iprocneigh(Lupx), 6*K+1, comm2d, requests(irq+1), ierr)
          CALL MPI_ISEND (BndRadRealInfo(:,:,Lupx), 5*NINVALID(Lupx), MPI_DOUBLE_PRECISION, &
                          iprocneigh(Lupx), 6*K+2, comm2d, requests(irq+2), ierr)
          CALL MPI_ISEND (BndRadIntInfo(:,Lupx), NINVALID(Lupx), MPI_INTEGER, &
                          iprocneigh(Lupx), 6*K+3, comm2d, requests(irq+3), ierr)
          irq=irq+3
        ENDIF
        IF (BTEST(BCFLAG,3)) THEN
          CALL MPI_ISEND (NINVALID(Lupy), 1, MPI_INTEGER, &
                          iprocneigh(Lupy), 6*K+4, comm2d, requests(irq+1), ierr)
          CALL MPI_ISEND (BndRadRealInfo(:,:,Lupy), 5*NINVALID(Lupy), MPI_DOUBLE_PRECISION, &
                          iprocneigh(Lupy), 6*K+5, comm2d, requests(irq+2), ierr)
          CALL MPI_ISEND (BndRadIntInfo(:,Lupy), NINVALID(Lupy), MPI_INTEGER, &
                          iprocneigh(Lupy), 6*K+6, comm2d, requests(irq+3), ierr)
          irq=irq+3
        ENDIF
        IF (BTEST(BCFLAG,2)) THEN
          CALL MPI_IRECV (NINVALID(Ldownx), 1, MPI_INTEGER, &
                          iprocneigh(Ldownx), 6*K+1, comm2d, requests(irq+1), ierr)
          CALL MPI_IRECV (BndRadRealInfo(:,:,Ldownx), 5*MAXNRAD, MPI_DOUBLE_PRECISION, &
                          iprocneigh(Ldownx), 6*K+2, comm2d, requests(irq+2), ierr)
          CALL MPI_IRECV (BndRadIntInfo(:,Ldownx), MAXNRAD, MPI_INTEGER, &
                          iprocneigh(Ldownx), 6*K+3, comm2d, requests(irq+3), ierr)
          irq=irq+3
        ENDIF
        IF (BTEST(BCFLAG,3)) THEN
          CALL MPI_IRECV (NINVALID(Ldowny), 1, MPI_INTEGER, &
                          iprocneigh(Ldowny), 6*K+4, comm2d, requests(irq+1), ierr)
          CALL MPI_IRECV (BndRadRealInfo(:,:,Ldowny), 5*MAXNRAD, MPI_DOUBLE_PRECISION, &
                          iprocneigh(Ldowny), 6*K+5, comm2d, requests(irq+2), ierr)
          CALL MPI_IRECV (BndRadIntInfo(:,Ldowny), MAXNRAD, MPI_INTEGER, &
                          iprocneigh(Ldowny), 6*K+6, comm2d, requests(irq+3), ierr)
          irq=irq+3
        ENDIF
        CALL MPI_WAITALL (irq, requests, status, ierr)


        NINVALID(Lupx)=0 ; NINVALID(Lupy)=0
        IF (BTEST(BCFLAG,2)) THEN
           ! Continue the backwards ray integrations for the X boundary
          DO J = 1, NINVALID(Ldownx)
            IF (Ldownx == 1) BndRadRealInfo(1,J,Ldownx) = XGRID(1)
            IF (Ldownx == 2) BndRadRealInfo(1,J,Ldownx) = XGRID(NX) 
            X=BndRadRealInfo(1,J,Ldownx) ; Y=BndRadRealInfo(2,J,Ldownx)
            Z=BndRadRealInfo(3,J,Ldownx) ; TRANSMIT=BndRadRealInfo(4,J,Ldownx)
            RADIANCE=BndRadRealInfo(5,J,Ldownx)
            CALL INTEGRATE_SOURCE (BCFLAG, IPFLAG, &
                             NX, NY, NZ, NPTS, NCELLS, &
                             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                             XGRID, YGRID, ZGRID, GRIDPOS, &
                             MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1, &
                             X, Y, Z, XE,YE,ZE, SIDE, &
                             TRANSMIT, RADIANCE, VALIDRAD)
             ! If we got a radiance then store it in RadSend, otherwise 
             ! save information needed to pass to neighboring processors
            IF (VALIDRAD) THEN
              Ndone = Ndone + 1
              IF (Ndone>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 2')
              RadSend(Ndone) = RADIANCE
              IndexSend(Ndone) = BndRadIntInfo(J,Ldownx)
            ELSE
              IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
                WRITE (6,*) 'Bad SIDE 2x',myproc,muout,phiout,x,y,z,SIDE,xe,ye,ze
                CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR')
              ENDIF
              NINVALID(SIDE) = NINVALID(SIDE) + 1
              IF (NINVALID(SIDE)>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 2')
              BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ XE, YE, ZE, TRANSMIT, RADIANCE /)
              BndRadIntInfo(NINVALID(SIDE),SIDE) = BndRadIntInfo(J,Ldownx)
            ENDIF
          ENDDO
        ENDIF

        IF (BTEST(BCFLAG,3)) THEN
           ! Continue the backwards ray integrations for the Y boundary
          DO J = 1, NINVALID(Ldowny)
            IF (Ldowny == 3) BndRadRealInfo(2,J,Ldowny) = YGRID(1)
            IF (Ldowny == 4) BndRadRealInfo(2,J,Ldowny) = YGRID(NY) 
            X=BndRadRealInfo(1,J,Ldowny) ; Y=BndRadRealInfo(2,J,Ldowny)
            Z=BndRadRealInfo(3,J,Ldowny) ; TRANSMIT=BndRadRealInfo(4,J,Ldowny)
            RADIANCE=BndRadRealInfo(5,J,Ldowny)
            CALL INTEGRATE_SOURCE (BCFLAG, IPFLAG, &
                             NX, NY, NZ, NPTS, NCELLS, &
                             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
                             XGRID, YGRID, ZGRID, GRIDPOS, &
                             MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1, &
                             X, Y, Z, XE,YE,ZE, SIDE, &
                             TRANSMIT, RADIANCE, VALIDRAD)
             ! If we got a radiance then store it in RadSend, otherwise 
             ! save information needed to pass to neighboring processors
            IF (VALIDRAD) THEN
              Ndone = Ndone + 1
              IF (Ndone>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 3')
              RadSend(Ndone) = RADIANCE
              IndexSend(Ndone) = BndRadIntInfo(J,Ldowny)
            ELSE
              IF (SIDE /= Lupx .AND. SIDE /= Lupy) THEN
                WRITE (6,*) 'Bad SIDE 2y',myproc,muout,phiout,x,y,z,SIDE,xe,ye,ze
                CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR')
              ENDIF
              NINVALID(SIDE) = NINVALID(SIDE) + 1
              IF (NINVALID(SIDE)>MAXNRAD) CALL ABORT_SHDOM_MPI ('COMPUTE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 3')
              BndRadRealInfo(:,NINVALID(SIDE),SIDE) = (/ XE, YE, ZE, TRANSMIT, RADIANCE /)
              BndRadIntInfo(NINVALID(SIDE),SIDE) = BndRadIntInfo(J,Ldowny)
            ENDIF
          ENDDO
        ENDIF

         ! See if there are any more invalid radiance rays on all the processors
        intsendbuf = (/ NINVALID(Lupx), NINVALID(Lupy) /)
        CALL MPI_ALLREDUCE (intsendbuf, intrecvbuf, 2, MPI_INTEGER, MPI_MAX, &
                            comm2d, ierr)
        AllDone = ALL(intrecvbuf(1:2) == 0)
      ENDDO

      IF (myproc > 0) THEN
         ! Send the finished radiances to the master processor
        CALL MPI_ISEND (Ndone, 1, MPI_INTEGER, 0, 6*K+7, comm2d, requests(1), ierr)
        CALL MPI_ISEND (RadSend, Ndone, MPI_REAL, 0, 6*K+8, comm2d, requests(2), ierr)
        CALL MPI_ISEND (IndexSend, Ndone, MPI_INTEGER, 0, 6*K+9, comm2d, requests(3), ierr)
        CALL MPI_WAITALL (3, requests, status, ierr)
      ELSE
         ! Put the radiances finished by the master processor in the list
        DO i = 1, Ndone
          RADOUT(IndexSend(i)) = RadSend(i)
        ENDDO
         ! Receive the radiances from all the other processors
        irq=0
        DO iproc = 1, numproc-1
          CALL MPI_IRECV (NradRecv(iproc), 1, MPI_INTEGER, iproc, &
                          6*K+7, comm2d, requests(irq+1), ierr)
          CALL MPI_IRECV (RadRecv(:,iproc), MAXNRAD, MPI_REAL, iproc, &
                          6*K+8, comm2d, requests(irq+2), ierr)
          CALL MPI_IRECV (IndexRecv(:,iproc), MAXNRAD, MPI_INTEGER, iproc, &
                          6*K+9, comm2d, requests(irq+3), ierr)
          irq = irq + 3
        ENDDO
        CALL MPI_WAITALL (irq, requests, status, ierr)
        DO iproc = 1, numproc-1
          DO i = 1, NradRecv(iproc)
            RADOUT(IndexRecv(i,iproc)) = RadRecv(i,iproc)
          ENDDO
        ENDDO
      ENDIF

    ENDIF
  ENDDO


  DEALLOCATE (requests, status, RadSend, IndexSend)
  IF (myproc==0) DEALLOCATE (NradRecv, RadRecv, IndexRecv)
  DEALLOCATE (BndRadRealInfo, BndRadIntInfo)
END SUBROUTINE COMPUTE_RADIANCE_PAR





SUBROUTINE CALC_ACCEL_SOLCRIT (DOACCEL, DELJDOT, DELJOLD, DELJNEW, JNORM, &
                               ACCELPAR, SOLCRIT)
 ! Does an MPI_ALLREDUCE to sum the delta source function vector dot 
 ! products over all the processors, and from these calculates the
 ! full-domain-wide acceleration parameter and solution criterion.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: DOACCEL
  REAL,    INTENT(INOUT) :: DELJDOT, DELJOLD, DELJNEW, JNORM
  REAL,    INTENT(OUT) :: ACCELPAR, SOLCRIT
  INTEGER :: ierr
  REAL    ::  R, THETA, sendbuf(4), recvbuf(4)
  REAL, SAVE :: A=0.0
  include 'mpif.h'

  IF (numproc > 1) THEN
     ! Call MPI_ALLREDUCE to sum up the four dot products
    sendbuf(:) = (/ DELJDOT, DELJOLD, DELJNEW, JNORM /)
    call MPI_ALLREDUCE (sendbuf, recvbuf, 4, MPI_REAL, MPI_SUM, comm2d, ierr)
    if (ierr .ne. MPI_SUCCESS) CALL ABORT_SHDOM_MPI ('CALC_ACCEL_SOLCRIT MPI_ALLREDUCE error')
    DELJDOT=recvbuf(1) ; DELJOLD=recvbuf(2)
    DELJNEW=recvbuf(3) ; JNORM=recvbuf(4)
  ENDIF

   ! Accelerate if desired, didn't last time, and things are converging. 
  IF (DOACCEL .AND. A .EQ. 0.0 .AND. DELJNEW .LT. DELJOLD) THEN
    ! Compute the acceleration extrapolation factor and apply it.
    R = SQRT(DELJNEW/DELJOLD)
    THETA = ACOS(DELJDOT/SQRT(DELJOLD*DELJNEW))
    A = (1 - R*COS(THETA) + R**(1+0.5*3.14159/THETA)) &
              /(1 + R**2  - 2*R*COS(THETA))  - 1.0
    A = MIN(10.0,MAX(0.0,A))
    ! WRITE (6,'(1X,A,3(1X,F7.3))') '! Acceleration: ', A,R,THETA
  ELSE
    A = 0.0
  ENDIF
  ACCELPAR = A

  IF (JNORM .GT. 0.0) THEN
    SOLCRIT = SQRT(DELJNEW/JNORM)
  ELSE IF (DELJNEW .EQ. 0.0) THEN
    SOLCRIT = 0.0
  ENDIF
END SUBROUTINE CALC_ACCEL_SOLCRIT




SUBROUTINE END_SHDOM_MPI (NPTS, GRIDPOS, NPX, NPY, XSTART, YSTART, DELX, DELY,&
                          NPXT, NPYT, PROPFILE, RUNNAME)
 ! Shuts down MPI at the end, but mainly finds the number of gridpoints
 ! around each property grid point column, sends that array to the master
 ! process, and the master process write out the RUNNAME//"_load_balance.out"
 ! file.
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NPTS, NPX, NPY, NPXT, NPYT
  REAL,    INTENT(IN) :: GRIDPOS(3,NPTS), XSTART, YSTART, DELX, DELY
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE, RUNNAME
  INTEGER :: I, IX, IY, MAXN, SX, SY, EX, EY, NX, NY
  INTEGER :: iproc, ierr, irq
  INTEGER, ALLOCATABLE :: NUMGRID(:,:), NUMGRIDALL(:,:), NORM(:,:), IXYALL(:,:,:)
  INTEGER, ALLOCATABLE :: recvbuf(:,:), buf(:,:), requests(:), status(:,:)
  REAL :: ElapsedTime
  include 'mpif.h'

  IF (LoadBalance) THEN
     ! Calculate the number of grid points around each property grid column
    ALLOCATE (NUMGRID(NPX,NPY))
    NUMGRID(:,:) = 0
    DO I = 1, NPTS
      IX = NINT((GRIDPOS(1,I)-XSTART)/DELX)+1
      IY = NINT((GRIDPOS(2,I)-YSTART)/DELY)+1
      IF (IX <= NPX .AND. IY <= NPY) THEN
        NUMGRID(IX,IY) = NUMGRID(IX,IY) + 1
      ENDIF
    ENDDO
 

    IF (numproc > 1) THEN
       ! Have the master process collect the property grid locations of all 
      ALLOCATE (IXYALL(2,2,0:numproc-1))
      CALL MPI_GATHER (IXYPRP,4,MPI_INTEGER, IXYALL,4,MPI_INTEGER, 0, comm2d, ierr)

       ! Send the NUMGRIDs to the master processor
      ALLOCATE (requests(numproc), status(MPI_STATUS_SIZE,numproc))
      irq=0
      IF (myproc > 0) THEN
        irq = irq + 1
        CALL MPI_ISEND (NUMGRID, NPX*NPY, MPI_INTEGER, 0, 777, comm2d, &
                        requests(irq), ierr)
      ELSE
        MAXN = (MAXVAL(IXYALL(2,1,:)-IXYALL(1,1,:))+1)*(MAXVAL(IXYALL(2,2,:)-IXYALL(1,2,:))+1)
        ALLOCATE (recvbuf(MAXN,1:numproc-1))
        DO iproc = 1, numproc-1
          irq = irq + 1
          CALL MPI_IRECV (recvbuf(:,iproc), MAXN, MPI_INTEGER, iproc, 777, &
                          comm2d, requests(irq), ierr)
        ENDDO
      ENDIF
      CALL MPI_WAITALL (irq, requests, status, ierr)

       ! Put the NUMGRIDS into the full domain array on the master processor.
       !   Use a normalization counter (NORM) to deal with grid points that 
       !   get contribution from two processors.
      IF (myproc == 0) THEN
        ALLOCATE (NUMGRIDALL(NPXT,NPYT), NORM(NPXT,NPYT))
        NUMGRIDALL(:,:) = 0  ;  NORM(:,:) = 0
        SX=IXYALL(1,1,0) ; EX=IXYALL(2,1,0) ; NX=EX-SX+1
        SY=IXYALL(1,2,0) ; EY=IXYALL(2,2,0) ; NY=EY-SY+1
        NUMGRIDALL(SX:EX,SY:EY) = NUMGRID(:,:)
        NORM(SX:EX,SY:EY) = RESHAPE( (/ 1 /), (/ NX,NY /), (/ 1 /))
        DO iproc = 1, numproc-1
          SX=IXYALL(1,1,iproc) ; EX=IXYALL(2,1,iproc) ; NX=EX-SX+1
          SY=IXYALL(1,2,iproc) ; EY=IXYALL(2,2,iproc) ; NY=EY-SY+1
          ALLOCATE (buf(NX,NY))
          buf(:,:) = RESHAPE(recvbuf(1:NX*NY,iproc),(/ NX,NY /) )
          NUMGRIDALL(SX:EX-1,SY:EY-1) = NUMGRIDALL(SX:EX-1,SY:EY-1) + buf(1:NX-1,1:NY-1)
          NORM(SX:EX-1,SY:EY-1) = NORM(SX:EX-1,SY:EY-1) + RESHAPE( (/ 1 /), (/ NX-1,NY-1 /), (/ 1 /))
          IX=MOD(EX-1,NPXT)+1 ; IY=MOD(EY-1,NPYT)+1
          NUMGRIDALL(IX,SY:EY-1) = NUMGRIDALL(IX,SY:EY-1) + buf(NX,1:NY-1)
          NORM(IX,SY:EY-1) = NORM(IX,SY:EY-1) + RESHAPE( (/ 1 /), (/ NY-1 /), (/ 1 /))
          NUMGRIDALL(SX:EX-1,IY) = NUMGRIDALL(SX:EX-1,IY) + buf(1:NX-1,NY)
          NORM(SX:EX-1,IY) = NORM(SX:EX-1,IY) + RESHAPE( (/ 1 /), (/ NX-1 /), (/ 1 /))
          NUMGRIDALL(IX,IY) = NUMGRIDALL(IX,IY) + buf(NX,NY)
          NORM(IX,IY) = NORM(IX,IY) + 1
          DEALLOCATE (buf)
        ENDDO
        DEALLOCATE (recvbuf)
        NUMGRIDALL(:,:) = NUMGRIDALL(:,:)/NORM(:,:)
        DEALLOCATE (NORM)
      ENDIF
      DEALLOCATE (requests, status, IXYALL)

    ELSE
      ALLOCATE (NUMGRIDALL(NPXT,NPYT))
      NUMGRIDALL(1:NPX,1:NPY) = NUMGRID(:,:)
    ENDIF
    DEALLOCATE (NUMGRID)

    IF (myproc == 0) THEN
        ! Write out the array to the load balancing file
      OPEN (UNIT=10, FILE=TRIM(RUNNAME)//'_load_balance.out', STATUS='UNKNOWN')
      WRITE (10,'(A)') '! SHDOM grid point load balancing file for property file:'
      WRITE (10,'(A)') PROPFILE
      WRITE (10,'(2(1X,I4),A)') NPXT, NPYT, ' ! property grid size (Nx,Ny)'
      DO IX = 1, NPXT
        DO IY = 1, NPYT
          WRITE (10,'(I4,1X,I4,1X,I5)') IX, IY, NUMGRIDALL(IX,IY)
        ENDDO
      ENDDO
      CLOSE (10)
      DEALLOCATE (NUMGRIDALL)
    ENDIF
  ENDIF

  ElapsedTime = MPI_WTIME() - StartTime
  WRITE (6,'(A,F9.2)') 'Total elapsed time (sec): ', ElapsedTime
  WRITE (6,*) 'Finished with SHDOM run: ',TRIM(RUNNAME)

  ! Shut down MPI
  IF (ALLOCATED(BXYPTR)) DEALLOCATE (BXYPTR,IBXYPTR,BXYPOS)
  CALL MPI_FINALIZE (ierr)
END SUBROUTINE END_SHDOM_MPI



SUBROUTINE ABORT_SHDOM_MPI (ERRSTR)
  USE SHDOM_MPI_DATA
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: ERRSTR
  INTEGER :: ierr, errorcode=7

  WRITE (6,*) ERRSTR
  CALL MPI_ABORT (comm2d,errorcode,ierr)
END SUBROUTINE ABORT_SHDOM_MPI


