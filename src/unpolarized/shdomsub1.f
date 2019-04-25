
        SUBROUTINE SOLVE_RTE (NSTOKES, NX, NY, NZ, NX1, NY1, NANG,
     .               ML, MM, NCS, NLM, NMU, NPHI, NLEG, NSTLEG, 
     .               NUMPHASE, NPHI0, MU, PHI, WTDO,
     .               MAXIV, MAXIC, MAXIG, MAXIDO, INRADFLAG, 
     .               BCFLAG, IPFLAG, DELTAM, SRCTYPE, HIGHORDERRAD,
     .               SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .               SFCTYPE, GNDTEMP, GNDALBEDO, 
     .               NXSFC, NYSFC, DELXSFC, DELYSFC, 
     .               NSFCPAR, SFCPARMS, SFCGRIDPARMS, 
     .               UNITS, WAVENO, WAVELEN, 
     .               ACCELFLAG, SOLACC, MAXITER, SOLCRIT, ITER,
     .               SPLITACC, SHACC, XGRID,YGRID,ZGRID,
     .               TEMP, PLANCK, EXTINCT, ALBEDO, LEGEN, IPHASE,
     .               MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .               NPTS, GRIDPOS, 
     .               NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .               RSHPTR, SHPTR, OSHPTR, WORK, WORK1, WORK2,
     .               SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX, 
     .               YLMSUN, VERBOSE, 
     .               NPX, NPY, NPZ, DELX, DELY,
     .               XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .               ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .               ZCKD, GASABS, OLDNPTS)
Cf2py threadsafe
C       Performs the SHDOM solution procedure.
C      Output is returned in SOURCE, RADIANCE, FLUXES, DIRFLUX.
C       The adaptive grid stucture (NCELLS,GRIDPTR,NEIGHPTR,TREEPTR,
C       CELLFLAGS,NPTS,GRIDPOS) is input and modified for new grid
C       cells and points on output.
      IMPLICIT NONE
      INTEGER NSTOKES, NX, NY, NZ, NX1, NY1, NXSFC, NYSFC, NSFCPAR
Cf2py intent(in) :: NSTOKES, NX, NY, NX1, NY1, NZ, NXSFC, NYSFC, NSFCPAR
      INTEGER ML, MM, NCS, NLM, NMU, NPHI, NANG, NLEG, NSTLEG, NUMPHASE
Cf2py intent(in) :: ML, MM, NCS, NLM, NMU, NPHI, NLEG, NSTLEG, NUMPHASE
Cf2py intent(out) :: NANG
      INTEGER NPHI0(NMU), MAXITER, ITER, BCFLAG, IPFLAG
Cf2py intent(in) :: MAXITER, BCFLAG, IPFLAG
Cf2py intent(out) :: NPHI0, ITER
      INTEGER MAXIV, MAXIC, MAXIG, MAXIDO
Cf2py intent(in) :: MAXIV, MAXIC, MAXIG, MAXIDO
      INTEGER MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR(MAXNBC,2)
Cf2py intent(out) :: NTOPPTS, NBOTPTS
Cf2py intent(in, out) :: BCPTR
Cf2py intent(in) :: MAXBCRAD
      INTEGER NPTS, NCELLS
Cf2py intent(in, out) :: NPTS, NCELLS
      INTEGER RSHPTR(*), SHPTR(*), OSHPTR(*)
Cf2py intent(in, out) :: RSHPTR, SHPTR, OSHPTR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in, out) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER*2 CELLFLAGS(*)
      INTEGER  IPHASE(*)
Cf2py intent(in, out) :: CELLFLAGS, IPHASE
      INTEGER WORK1(*)
Cf2py intent(in) :: WORK1
      LOGICAL INRADFLAG, DELTAM, ACCELFLAG, HIGHORDERRAD
Cf2py intent(in) :: INRADFLAG, DELTAM, ACCELFLAG, HIGHORDERRAD
      REAL    SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD
Cf2py intent(in) :: SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD
      REAL    GNDTEMP, GNDALBEDO
Cf2py intent(in) :: GNDTEMP, GNDALBEDO
      REAL    DELXSFC, DELYSFC, SFCPARMS(*), SFCGRIDPARMS(*)
Cf2py intent(in) :: DELXSFC, DELYSFC, SFCPARMS
Cf2py intent(in, out) :: SFCGRIDPARMS
      REAL    WAVENO(2), WAVELEN
Cf2py intent(in) :: WAVENO, WAVELEN
      REAL    SOLACC, SOLCRIT, SPLITACC, SHACC
Cf2py intent(in) :: SOLACC, SPLITACC, SHACC
Cf2py intent(in, out) :: SOLCRIT
      REAL    MU(*), WTDO(*), PHI(*)
Cf2py intent(in, out) :: MU, WTDO, PHI
      REAL    TEMP(*), PLANCK(*)
Cf2py intent(in, out) :: TEMP, PLANCK
      REAL    EXTINCT(MAXIG), ALBEDO(MAXIG), LEGEN(0:NLEG,NUMPHASE)
Cf2py intent(out) :: EXTINCT, ALBEDO, LEGEN
      REAL    XGRID(*), YGRID(*), ZGRID(*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID
      REAL    GRIDPOS(3,*)
Cf2py intent(in, out) :: GRIDPOS
      REAL    BCRAD(*)
Cf2py intent(in, out) :: BCRAD
      REAL    FLUXES(2,*), DIRFLUX(*)
Cf2py intent(in, out) :: FLUXES,  FLUXES1, DIRFLUX
      REAL    SOURCE(*), DELSOURCE(*), RADIANCE(*)
Cf2py intent(in, out) :: SOURCE, DELSOURCE, RADIANCE
      REAL    WORK(*), WORK2(*)
Cf2py intent(in) :: WORK, WORK2
      CHARACTER SRCTYPE*1, UNITS*1, SFCTYPE*2
Cf2py intent(in) :: SRCTYPE, UNITS, SFCTYPE
      LOGICAL VERBOSE
Cf2py intent(in) :: VERBOSE
      INTEGER MAXNLM, MAXNMU, MAXNPHI

      REAL    YLMSUN(*)
Cf2py intent(in, out) :: YLMSUN

      INTEGER NPX, NPY, NPZ
Cf2py intent(in) :: NPX, NPY, NPZ
      REAL DELX, DELY, XSTART, YSTART
Cf2py intent(in) :: DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)
Cf2py intent(in) :: ZLEVELS
      REAL TEMPP(*), EXTINCTP(*), ALBEDOP(*)
Cf2py intent(in) :: TEMPP, EXTINCTP, ALBEDOP
      REAL LEGENP(*), EXTDIRP(*)
Cf2py intent(in) :: LEGENP, EXTDIRP
      INTEGER IPHASEP(*)
Cf2py intent(in) :: IPHASEP
      INTEGER NZCKD
Cf2py intent(in) :: NZCKD
      REAL ZCKD(*), GASABS(*)
Cf2py intent(in) :: ZCKD, GASABS
      INTEGER OLDNPTS 
Cf2py intent(in,out) :: OLDNPTS

      REAL A
      INTEGER SP, STACK(50)
      DOUBLE PRECISION EXTMIN, SCATMIN
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      INTEGER IPDIRECT, DI, DJ, DK
      DOUBLE PRECISION EPSS, EPSZ, XDOMAIN, YDOMAIN
      DOUBLE PRECISION UNIFORMZLEV, DELXD,DELYD
      INTEGER IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ
      
      INTEGER NPHI0MAX, I, ORDINATESET, NBPTS, NBCELLS
      INTEGER MAXNANGBND, NANGBND(4,2)
      LOGICAL FIXSH, SPLITTESTING, DOSPLIT, OUTOFMEM, LAMBERTIAN

      LOGICAL, ALLOCATABLE :: FFTFLAG(:)
      REAL    ALBMAX
      REAL    STARTADAPTSOL, ENDADAPTSOL, ADAPTRANGE, SPLITCRIT
      REAL    STARTSPLITACC, CURSPLITACC, AVGSOLCRIT, BETA
      REAL    DELJDOT, DELJOLD, DELJNEW, JNORM, ACCELPAR
      
      REAL, ALLOCATABLE :: CMU1(:,:), CMU2(:,:)
      REAL, ALLOCATABLE :: CPHI1(:,:,:), CPHI2(:,:,:)
      REAL, ALLOCATABLE :: WTMU(:), WPHISAVE(:)
      REAL, ALLOCATABLE :: SFC_BRDF_DO(:,:,:,:,:,:)
C       Allocate local arrays
C         NSTLEG=1 for NSTOKES=1, otherwise NSTLEG=6
      ALLOCATE (CMU1(NLM,NMU), CMU2(NLM,NMU))
      ALLOCATE (CPHI1(-16:16,32,NMU), CPHI2(32,-16:16,NMU))
      ALLOCATE (FFTFLAG(NMU), WTMU(NMU), WPHISAVE(NMU*(3*NPHI+15)))
      
      
C       Set up some things before solution loop
C           Transfer the medium properties to the internal grid and add gas abs
      CALL INTERP_GRID (NPTS, NLEG, GRIDPOS,
     .                  TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .                  NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .                  XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .                  ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .                  ZCKD, GASABS, EXTMIN, SCATMIN)
 
C         If Delta-M then scale the extinction, albedo, and Legendre terms.
C         Put the Planck blackbody source in PLANCK.
      CALL PREPARE_PROP (ML, MM, NLEG, NPTS, DELTAM, NUMPHASE, SRCTYPE,
     .                   UNITS, WAVENO, WAVELEN, ALBMAX, 
     .                   EXTINCT, ALBEDO, LEGEN, TEMP, PLANCK, IPHASE)
C       Get the maximum single scattering albedo over all processors
      CALL TOTAL_ALBEDO_MAX (ALBMAX)
C           Compute the solar transmission in DIRFLUX. 
C           Precompute Ylm's for solar direction
      IF (SRCTYPE .NE. 'T') THEN
        IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3)) THEN
          CALL MAKE_DIRECT_PAR (1,NPTS, BCFLAG, IPFLAG, DELTAM,ML,NLEG,
     .                          SOLARFLUX, SOLARMU, SOLARAZ,  GRIDPOS, 
     .                          NX, XGRID, NY, YGRID, DIRFLUX)
        ELSE
          CALL MAKE_DIRECT (NPTS, BCFLAG, IPFLAG, DELTAM, ML, NLEG,
     .             SOLARFLUX, SOLARMU, SOLARAZ,  GRIDPOS, DIRFLUX, 
     .             NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .             XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .             ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .             ZCKD, GASABS, CX, CY, CZ, CXINV, CYINV,
     .             CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .             XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV)
        ENDIF
        CALL YLMALL (SOLARMU, SOLARAZ, ML, MM, NCS, YLMSUN)
      ENDIF
C           Make the discrete ordinates (angles) 
C             (set 2 is reduced gaussian, 3 is reduced double gauss)
      ORDINATESET = 2
      CALL MAKE_ANGLE_SET (NMU, NPHI, NCS, ORDINATESET,  NPHI0MAX,
     .         NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

C           Make the Ylm transform coefficients
      CALL MAKE_SH_DO_COEF (ML, MM, NLM, NMU, NPHI0, NCS, 
     .         NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG,
     .         CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)

      FIXSH = .FALSE.
      IF (.NOT. INRADFLAG) THEN
        OLDNPTS = 0
C           Initialize the radiance on the base grid using Eddington 
C             two-stream plane-parallel
        CALL INIT_RADIANCE (NX1*NY1, NZ, NCS, NLEG, RSHPTR, ZGRID,
     .         EXTINCT, ALBEDO, LEGEN, TEMP, NUMPHASE, IPHASE,
     .         SRCTYPE, SOLARFLUX, SOLARMU, GNDALBEDO, GNDTEMP, SKYRAD, 
     .         UNITS, WAVENO, WAVELEN, RADIANCE)
C           Interpolate the radiance on the non-base grid points from the 
C             base grid points
        NBPTS = NX1*NY1*NZ
        NBCELLS = (NZ-1)*(NX+IBITS(BCFLAG,0,1)-IBITS(BCFLAG,2,1))
     .                  *(NY+IBITS(BCFLAG,1,1)-IBITS(BCFLAG,3,1))
        CALL INTERP_RADIANCE (NBPTS, NPTS, RSHPTR, RADIANCE,
     .             NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
C           Initialize the source function from the radiance field
        
        CALL COMPUTE_SOURCE (ML,MM, NCS, NLM, NLEG, NUMPHASE, NPTS,
     .           FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO,
     .           LEGEN, IPHASE, PLANCK, DIRFLUX, SHACC, MAXIV,
     .           RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE, 
     .           .TRUE.,ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)
      ENDIF
      IF (ACCELFLAG) THEN
        OSHPTR(1:NPTS+1) = SHPTR(1:NPTS+1)
        DELSOURCE(1:OSHPTR(NPTS+1)) = 0.0
      ENDIF

      LAMBERTIAN = SFCTYPE(2:2) .EQ. 'L'
C         Make the pointers to the grid points on the top and bottom boundaries
      CALL BOUNDARY_PNTS (NPTS, NANG, LAMBERTIAN, MAXNBC, MAXBCRAD, 
     .                    ZGRID(1), ZGRID(NZ), GRIDPOS,
     .                    NTOPPTS, NBOTPTS, BCPTR)
      IF (SFCTYPE(1:1) .EQ. 'V') THEN
C             If this is a variable surface then bilinearly interpolate 
C               the surface temperature and reflection parameters to the
C               bottom grid points. 
        CALL SURFACE_PARM_INTERP (NBOTPTS, BCPTR(1,2), GRIDPOS,
     .             SRCTYPE, WAVENO, WAVELEN, UNITS,
     .             NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,
     .             SFCGRIDPARMS)
        IF (SFCTYPE(2:2) .EQ. 'O') THEN
          CALL CHECK_OCEAN_BRDF_INTEG (NXSFC, NYSFC, NSFCPAR, SFCPARMS,
     .                                 WAVELEN, SOLARMU, SOLARAZ,
     .                                 NMU, NPHI0, NPHI0MAX,
     .                                 MU, PHI, WTDO)
        ENDIF
      ENDIF


      ITER = 0
      IF (MAXITER .LE. ITER) RETURN
C         Starting values for the adaptive cell splitting controlling method.
C         The splitting accuracy used at a particular iteration is CURSPLITACC,
C         and it is started out high an lowered to the desired SPLITACC.
C         The cell splitting starts when the solution criterion (SOLCRIT)
C         goes below STARTADAPTSOL and ends when the solution criterion
C         goes below ENDADAPTSOL, though the current splitting accuracy
C         is designed to get to the final SPLITACC over a smaller range
C         (ADAPTRANGE).  The currect splitting accuracy is computed from
C         a geometric average of the solution criterion (AVGSOLCRIT).
C         The SPLIT_GRID routine is called and the splitting criterion returned
C         as long as SPLITTESTING is true.
      SPLITTESTING = .TRUE.
      OUTOFMEM = .FALSE.
      ENDADAPTSOL = 0.001
      STARTADAPTSOL = 0.1
      ADAPTRANGE = STARTADAPTSOL/(3.0*ENDADAPTSOL)
      CURSPLITACC = SPLITACC*ADAPTRANGE
      STARTSPLITACC = CURSPLITACC
      AVGSOLCRIT = SOLCRIT
      SPLITCRIT = 0.0
      A = 0.0
      
      IF (VERBOSE) THEN
        WRITE (6,*) '! Iter Log(Sol)  SplitCrit  Npoints  Nsh(avg)'
      ENDIF
 
C         Main solution loop
C           Iterate until the iterations exceed the limit, or the desired
C           solution criterion is reached and the current splitting accuracy
C           is down to the desired level.
      DO WHILE (ITER .LT. MAXITER .AND. (SOLCRIT .GT. SOLACC
     .        .OR. (SPLITCRIT .GT. SPLITACC .AND.
     .              CURSPLITACC .GT. SPLITACC .AND. .NOT. OUTOFMEM)) )
        ITER = ITER + 1

C           Don't do the adaptive grid cell splitting stuff unless wanted
        IF (SPLITACC .GT. 0.0) THEN
C             If in the right solution range then try to split some cells, and
C               interpolate the medium and source functions to the new points.
C             The current splitting accuracy is started high and then lowered.
          AVGSOLCRIT = SQRT(AVGSOLCRIT*SOLCRIT)
          DOSPLIT = SOLCRIT .LE. STARTADAPTSOL .AND. 
     .        (SOLCRIT .GT. ENDADAPTSOL .OR. CURSPLITACC .GT. SPLITACC)
     .         .AND. .NOT. OUTOFMEM
C           Make sure all processors are going to split cells if any want to
          CALL UNIFY_SPLITTING (DOSPLIT,STARTSPLITACC)
          BETA = LOG(STARTSPLITACC/SPLITACC)/LOG(ADAPTRANGE)
          CURSPLITACC = MIN(CURSPLITACC, MAX(SPLITACC,
     .                 SPLITACC*(AVGSOLCRIT/(3.0*ENDADAPTSOL))**BETA ))
          IF (SOLCRIT .LE. ENDADAPTSOL) CURSPLITACC = SPLITACC
          IF (SPLITTESTING) THEN
            CALL SPLIT_GRID (MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,
     .             DOSPLIT,OUTOFMEM, CURSPLITACC, SPLITCRIT,
     .             NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, 
     .             GRIDPOS, NX, XGRID, NY, YGRID,
     .             ML, MM, NLM, NCS, NLEG, NUMPHASE, 
     .             DELTAM, BCFLAG, IPFLAG, ACCELFLAG,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .             RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, 
     .             TEMP, PLANCK, DIRFLUX, SRCTYPE, 
     .             SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN, 
     .             UNITS, WAVENO, WAVELEN,  WORK, WORK2, 
     .             NPX, NPY, NPZ, DELX, DELY,
     .             XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .             ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .             ZCKD, GASABS, EXTMIN, SCATMIN, CX, CY, CZ, CXINV, 
     .             CYINV, CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .             XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV)
            IF (SOLCRIT .GT. STARTADAPTSOL)  STARTSPLITACC = SPLITCRIT
          ENDIF
C           Find the maximum splitting criterion over all processors
          CALL TOTAL_SPLITCRIT_MAX (SPLITCRIT)
          IF (SOLCRIT .LE. ENDADAPTSOL)  SPLITTESTING = .FALSE.
        ENDIF

C            Find the SH truncation to use for the radiance (RSHPTR)
        CALL RADIANCE_TRUNCATION (HIGHORDERRAD, 
     .           ML, MM, NCS, NLEG, NUMPHASE, 
     .           NPTS, ALBEDO, LEGEN, IPHASE, SHPTR, RADIANCE,
     .           MAXIV+MAXIG, FIXSH, SHACC, RSHPTR)

C           Integrate the source function along discrete ordinates 
C             to compute the radiance field.  The input source function 
C             and output radiance are in spherical harmonic expansions,
C             and the SH/DO transforms are done during the loop over 
C             discrete ordinates.
        CALL PATH_INTEGRATION (NX, NY, NZ, ML, MM, NLM, 
     .           NMU, NPHI0MAX, NPHI0, NANG, NCS, BCFLAG, IPFLAG,  
     .           NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS, 
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, 
     .           SRCTYPE, SOLARMU, SOLARAZ, SKYRAD, WAVENO, WAVELEN, 
     .           UNITS, SFCTYPE, GNDTEMP, GNDALBEDO, 
     .           NXSFC, NYSFC, DELXSFC, DELYSFC, 
     .           NSFCPAR, SFCPARMS, SFCGRIDPARMS, 
     .           MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .           MU, PHI, WTDO,
     .           CMU1, CPHI1, CMU2, CPHI2, WPHISAVE, FFTFLAG,
     .           DIRFLUX, FLUXES, EXTINCT, 
     .           SHPTR, SOURCE, RSHPTR, RADIANCE, 
     .           WORK, WORK1, WORK2, OLDNPTS, SP, STACK, IX, IY, 
     .           IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ)

C            Compute the source function from the radiance field,
C              do the adaptive spherical harmonics truncation, compute 
C              the solution criterion, and dot products for acceleration.
        IF (SOLCRIT .LT. ENDADAPTSOL .OR. ITER .GT. 30)  FIXSH = .TRUE.
        CALL COMPUTE_SOURCE (ML,MM, NCS, NLM, NLEG, NUMPHASE, NPTS,
     .           FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO,
     .           LEGEN, IPHASE, PLANCK, DIRFLUX, SHACC, MAXIV,
     .           RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE, 
     .           .FALSE.,ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)

C           Calculate the acceleration parameter and solution criterion
C            from all the processors
        CALL CALC_ACCEL_SOLCRIT (ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, 
     .                           JNORM, ACCELPAR, SOLCRIT, A)

C           Accelerate the convergence of the source function vector
        CALL ACCELERATE_SOLUTION (ACCELPAR, NPTS, SHPTR, OSHPTR, 
     .                            SOURCE, DELSOURCE)

C           If it is a not scattering medium then do only one iteration
        IF (ALBMAX .LT. SOLACC)  SOLCRIT = SOLACC
        
        IF (VERBOSE) THEN
C           Print out the log solution criterion, number of points, 
C             and average SH truncation.
          WRITE (6,'(2X,I4,F8.3,1X,E10.3,1X,I8,1X,F8.2,1X,F6.3)') 
     .          ITER, LOG10(MAX(SOLCRIT,1.0E-20)), SPLITCRIT,
     .          NPTS, FLOAT(SHPTR(NPTS+1))/NPTS, 
     .          FLOAT(SHPTR(NPTS+1))/(NPTS*NLM)
        ENDIF
      ENDDO
      
      IF (VERBOSE) THEN
        WRITE (6,'(1X,A,I6,A,F9.6)') '! Iterations: ', ITER,
     .                   '     Final Criterion: ', SOLCRIT
      ENDIF

C          Comment in for GLE graphics output of cell structure (see routine)
c      CALL VISUALIZE_CELLS ('cellfinal.gle', YGRID(1), 1, IPFLAG,
c     .           NX, XGRID(NX+1)-XGRID(1), ZGRID(NZ)-ZGRID(1),
c     .           NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, 
c     .           CELLFLAGS, GRIDPOS)

C          Comment in for output of cell splitting criterion for every cell
c      CALL OUTPUT_CELL_SPLIT ('cellsplit.dat',GRIDPTR,
c     .                        GRIDPOS, EXTINCT, SHPTR, SOURCE, NCELLS)

      RETURN
      END
 

      SUBROUTINE COMPUTE_SOURCE (ML,MM, NCS, NLM, NLEG, NUMPHASE,NPTS,
     .             FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO, 
     .             LEGEN, IPHASE, PLANCK, DIRFLUX, SHACC, MAXIV,
     .             RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,
     .             FIRST,ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)
C       Computes the source function (SOURCE) in spherical harmonic space 
C     for all the grid points.  The thermal source and/or solar 
C     pseudo-source (in PLANCK or DIRFLUX) is added to the scattering source 
C     (computed from LEGEN and the spherical harmonic expansion in RADIANCE).
C     Also computes the adaptive spherical harmonic truncation for
C     the source function, unless FIXSH is true.  The SH truncation is set
C     so that the higher terms of the source function are smaller in 
C     absolute value than SHACC. The truncation is to the nearest l,
C     so that a triangular (for M=L) truncation is preserved.  If there
C     is little or no scattering then the SH truncation can be set to
C     zero terms.  The source function SH truncation is output in terms
C     of the pointer SHPTR (note that RSHPTR refers to the RADIANCE array
C     and OSHPTR to the DELSOURCE array).
C       To save memory this routine also computes some things for the
C     series acceleration: the dot products of the difference in source 
C     function between successive iterations and the new difference source 
C     function vector (DELSOURCE).  To avoid using another array and 
C     because of the adaptive truncation pointer stuff, this requires 
C     two extra calculations of the source function.
C       Computes the sum of squares of the new source function for calculation
C     of the solution criterion (SOLCRIT, the RMS difference in succesive 
C     source function fields normalized by the RMS of the field).
      IMPLICIT NONE

      INTEGER NPTS, ML, MM, NCS, NLM, NLEG, MAXIV
Cf2py intent(in) :: NPTS, ML, MM, NCS, NLM, NLEG, MAXIV
      INTEGER NUMPHASE
      INTEGER RSHPTR(*), SHPTR(*), OSHPTR(*)
Cf2py intent(in) :: RSHPTR, SHPTR, OSHPTR
      INTEGER IPHASE(*)
Cf2py intent(in) :: IPHASE
      LOGICAL FIXSH, FIRST, ACCELFLAG
Cf2py intent(in) :: FIXSH, FIRST, ACCELFLAG
      REAL    SHACC, YLMSUN(*), SOLARMU
Cf2py intent(in) :: SHACC, YLMSUN, SOLARMU
      REAL    DELJDOT, DELJOLD, DELJNEW, JNORM
Cf2py intent(out) :: DELJDOT, DELJOLD, DELJNEW, JNORM
      REAL    ALBEDO(*), LEGEN(0:NLEG,*)
Cf2py intent(in) :: ALBEDO, LEGEN
      REAL    PLANCK(*), DIRFLUX(*)
Cf2py intent(in) :: PLANCK, DIRFLUX
      REAL    RADIANCE(*), SOURCE(*), DELSOURCE(*)
Cf2py intent(in) :: RADIANCE 
Cf2py intent(in,out) :: SOURCE, DELSOURCE
      CHARACTER SRCTYPE*1
Cf2py intent(in) :: SRCTYPE

      INTEGER IS, ISO, IR, I, J, JS, K, L, LS, M, MS, ME, NS, NR
      INTEGER, ALLOCATABLE :: LOFJ(:)
      REAL, ALLOCATABLE :: SOURCET(:)
      REAL    SRCMIN, C, SECMU0, D

      ALLOCATE (SOURCET(NLM), LOFJ(NLM))
      SOURCET = 0.0
      SRCMIN = SHACC
C         Make source function from scattering integral and real sources:
C           LEGEN has chil/(2l+1)
C             Solar source: DIRFLUX has F0*exp(-tau)
C             Thermal source: PLANCK has (1-omega)*B(T)
      SECMU0 = 1.0/ABS(SOLARMU)
      C = SQRT(4.0*ACOS(-1.0))
C         Make the l index as a function of SH term (J)
      J = 0
      DO L = 0, ML
        ME = MIN(L,MM)
        MS = (1-NCS)*ME
        DO M = MS, ME
          J = J + 1
          LOFJ(J) = L
        ENDDO
      ENDDO

C         Compute the dot products of the successive source differences
      IF (.NOT. FIRST) THEN
        DELJDOT = 0.0
        DELJOLD = 0.0
        DELJNEW = 0.0
        JNORM = 0.0
        DO I = 1, NPTS
C             Compute the temporary source function for this point
          IF (NUMPHASE .GT. 0) THEN
            K = IPHASE(I)
          ELSE
            K = I
          ENDIF
          IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
            D = DIRFLUX(I)*SECMU0
            DO J = 1, NLM
              SOURCET(J) = ALBEDO(I)*LEGEN(LOFJ(J),K)*D*YLMSUN(J)
            ENDDO
            IF (SRCTYPE .EQ. 'B') THEN
              SOURCET(1) = SOURCET(1) + C*PLANCK(I)
            ENDIF
          ELSE
            SOURCET(1) = C*PLANCK(I)
            DO J = 2, NLM
              SOURCET(J) = 0.0
            ENDDO
          ENDIF
          IR = RSHPTR(I)
          NR = RSHPTR(I+1)-IR
          DO J = 1, NR
            SOURCET(J) = SOURCET(J) 
     .                 + ALBEDO(I)*LEGEN(LOFJ(J),K)*RADIANCE(IR+J)
          ENDDO
C             Compute the dot product of the new and old source difference,
C               and the solution criterion (don't use new points for Jnorm).
          IF (ACCELFLAG) THEN
            IS = SHPTR(I)
            ISO = OSHPTR(I)
            NS = MIN(SHPTR(I+1)-IS,OSHPTR(I+1)-ISO)
            DO J = 1, NS
              DELJDOT=DELJDOT+(SOURCET(J)-SOURCE(IS+J))*DELSOURCE(ISO+J)
              DELJOLD = DELJOLD + DELSOURCE(ISO+J)**2
              DELJNEW = DELJNEW + (SOURCET(J)-SOURCE(IS+J))**2
              JNORM = JNORM + SOURCE(IS+J)**2
            ENDDO
          ELSE
            IS = SHPTR(I)
            NS = SHPTR(I+1)-IS
            DO J = 1, NS
              DELJNEW = DELJNEW + (SOURCET(J)-SOURCE(IS+J))**2
              JNORM = JNORM + SOURCE(IS+J)**2
            ENDDO
          ENDIF
        ENDDO
      ENDIF

C         Compute and store the new source function difference vector
      IF (.NOT. FIRST .AND. ACCELFLAG) THEN
        DO I = 1, NPTS
C             Compute the temporary source function for this point
          IF (NUMPHASE .GT. 0) THEN
            K = IPHASE(I)
          ELSE
            K = I
          ENDIF
          IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
            D = DIRFLUX(I)*SECMU0
            DO J = 1, NLM
              SOURCET(J) = ALBEDO(I)*LEGEN(LOFJ(J),K)*D*YLMSUN(J)
            ENDDO
            IF (SRCTYPE .EQ. 'B') THEN
              SOURCET(1) = SOURCET(1) + C*PLANCK(I)
            ENDIF
          ELSE
            SOURCET(1) = C*PLANCK(I)
            DO J = 2, NLM
              SOURCET(J) = 0.0
            ENDDO
          ENDIF
          IR = RSHPTR(I)
          NR = RSHPTR(I+1)-IR
          DO J = 1, NR
            SOURCET(J) = SOURCET(J) 
     .                 + ALBEDO(I)*LEGEN(LOFJ(J),K)*RADIANCE(IR+J)
          ENDDO
C             Make the source difference vector
          IS = SHPTR(I)
          NS = SHPTR(I+1)-IS
          OSHPTR(I) = IS
          DO J = 1, NS
            DELSOURCE(IS+J) = SOURCET(J) - SOURCE(IS+J)
          ENDDO
        ENDDO 
        OSHPTR(NPTS+1) = SHPTR(NPTS+1)
      ENDIF


C         Compute the source function in the temporary array.
C           Start with the thermal and/or solar source, then add
C           in the scattering from the radiance.
      IS = 0
      DO I = 1, NPTS
        IF (NUMPHASE .GT. 0) THEN
          K = IPHASE(I)
        ELSE
          K = I
        ENDIF
        IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
          D = DIRFLUX(I)*SECMU0
          DO J = 1, NLM
            SOURCET(J) = ALBEDO(I)*LEGEN(LOFJ(J),K)*D*YLMSUN(J)
          ENDDO
          IF (SRCTYPE .EQ. 'B') THEN
            SOURCET(1) = SOURCET(1) + C*PLANCK(I)
          ENDIF
        ELSE
          SOURCET(1) = C*PLANCK(I)
          DO J = 2, NLM
            SOURCET(J) = 0.0
          ENDDO
        ENDIF
        IR = RSHPTR(I)
        NR = RSHPTR(I+1)-IR
        DO J = 1, NR
          SOURCET(J) = SOURCET(J) 
     .                 + ALBEDO(I)*LEGEN(LOFJ(J),K)*RADIANCE(IR+J)
        ENDDO

        IF (FIXSH) THEN
          NS = SHPTR(I+1)-IS
        ELSE
C           Then do the adaptive spherical harmonics:
C             Find how many source function terms to use, and how many
C             terms wanted for the radiance.
          IF (SRCTYPE .EQ. 'S') THEN
            JS = 0
          ELSE
            JS = 1
          ENDIF
          DO J = 1, NLM
            IF (ABS(SOURCET(J)) .GT. SRCMIN)  JS = J
          ENDDO
C             Use all the m's for the last l for the source function.
          IF (JS .EQ. 0) THEN
            NS = 0
          ELSE 
            LS = LOFJ(JS)
            IF (LS .LE. MM) THEN
              NS = NCS*(LS*(LS+1))/2 + LS+1
            ELSE
              NS = (NCS*MM+1)*LS-(MM*(2+NCS*(MM-1)))/2 + MM+1
            ENDIF
          ENDIF
          SHPTR(I) = IS
        ENDIF
        IF (IS+NS .GT. MAXIV) THEN
          WRITE (6,*) 'COMPUTE_SOURCE: MAXIV exceeded.'
          WRITE (6,*) 'Out of memory for more spherical harmonic terms.'
          STOP
        ENDIF
C           Transfer the new source function
        DO J = 1, NS
          SOURCE(IS+J) = SOURCET(J)
        ENDDO
        IS = IS + NS        
      ENDDO 
      SHPTR(NPTS+1) = IS
      DEALLOCATE (SOURCET, LOFJ)
      RETURN
      END
 



      SUBROUTINE RADIANCE_TRUNCATION (HIGHORDERRAD, 
     .             ML, MM, NCS, NLEG, NUMPHASE, 
     .             NPTS, ALBEDO, LEGEN, IPHASE, SHPTR, RADIANCE,
     .             MAXIR, FIXSH, SHACC, RSHPTR)
C       Computes the next radiance spherical harmonic series truncation 
C     to use based on the current source function truncation and the 
C     scattering properties.  The radiance truncation is the smaller
C     of that based on a minimum SH term radiance and that needed for
C     the computation of the source function.  The truncation is always
C     made for a complete set in the L index.  The radiance SH truncation
C     is set to a few higher L index than the source function, so that
C     the source function truncation may grow (since initialized with
C     a low order truncation).  The minimum radiance truncation is 
C     L=1 so that the mean radiance and net flux terms may be computed,
C     unless HIGHORDERRAD is on, in which case all terms are kept.
      IMPLICIT NONE
      INTEGER NPTS, ML, MM, NCS, NLEG, MAXIR
Cf2py intent(in) :: NPTS, ML, MM, NCS, NLEG, MAXIR
      INTEGER NUMPHASE
      INTEGER SHPTR(*), RSHPTR(*)
Cf2py intent(in) ::SHPTR
Cf2py intent(in,out) :: RSHPTR
      INTEGER IPHASE(*)
Cf2py intent(in) :: IPHASE
      LOGICAL HIGHORDERRAD, FIXSH
Cf2py intent(in) ::  HIGHORDERRAD, FIXSH
      REAL    ALBEDO(*), LEGEN(0:NLEG,*), RADIANCE(*)
Cf2py intent(in) :: ALBEDO, LEGEN, RADIANCE
      REAL    SHACC
Cf2py intent(in) :: SHACC
      INTEGER NLM
      INTEGER IRO, IR, NRO, NR, I, J, K, L, LR, LS, M, MS, ME, NS
      INTEGER, ALLOCATABLE :: LOFJ(:)
      LOGICAL NOTEND
      REAL    RADMIN

      IF (FIXSH)  GOTO 190

      RADMIN = SHACC
      NOTEND = .TRUE.
C         Make the l index as a function of SH term (J)
      NLM = (2*MM+1)*(ML+1) - MM*(MM+1)
      ALLOCATE (LOFJ(NLM))
      J = 0
      DO L = 0, ML
        ME = MIN(L,MM)
        MS = (1-NCS)*ME
        DO M = MS, ME
          J = J + 1
          LOFJ(J) = L
        ENDDO
      ENDDO


      IRO = 0
      IR = 0
      RSHPTR(1) = 0
      DO I = 1, NPTS
C           Find how many terms wanted for the radiance based on the
C             radiance threshold (RADMIN) (use the mean radiance as an 
C             upper limit to the SH term radiance).
C             We don't have radiances for any new points, so don't use them
        NRO = RSHPTR(I+1)-IRO
        IF (NRO .EQ. 0)  NOTEND = .FALSE.
        IF (NOTEND) THEN 
          IF (NUMPHASE .GT. 0) THEN
            K = IPHASE(I)
          ELSE
            K = I
          ENDIF
          LR = 1
          DO L = 1, ML
            IF (ALBEDO(I)*LEGEN(L,K)*RADIANCE(IRO+1) .GT. RADMIN)
     .        LR = L
          ENDDO
        ELSE
          LR = ML
        ENDIF
        IRO = IRO + NRO
          
C           The radiance L truncation is the smaller of the radiance cutoff 
C             from above and a few more than the source truncation, but
C             it is at least 1 (so low order terms are always available).
C             For higher order radiance option, use all terms.
C             But if out of memory, use the source function truncation.
        NS = MAX(1,SHPTR(I+1)-SHPTR(I))
        LS = LOFJ(NS)
        LR = MIN(LR, LS+ML/8+2)
        IF (HIGHORDERRAD)  LR = ML
        IF (LR .LE. MM) THEN
          NR = NCS*(LR*(LR+1))/2 + LR+1
        ELSE
          NR = (NCS*MM+1)*LR-(MM*(2+NCS*(MM-1)))/2 + MM+1
        ENDIF
        IR = IR + NR
        RSHPTR(I+1) = IR
        IF (IR .GT. MAXIR) THEN
          WRITE (6,*) 
     .    'RADIANCE_TRUNCATION: Out of memory for more radiance terms.'
          GOTO 190
        ENDIF
      ENDDO 
      RSHPTR(NPTS+2) = IR
      DEALLOCATE (LOFJ)
      RETURN

C         If out of memory or the SH truncation is fixed then 
C           just set the radiance truncation to that of the source.
190   CONTINUE
      RSHPTR(1) = 0
      IR = 0
      DO I = 1, NPTS
        NR = MAX(NCS+2,SHPTR(I+1)-SHPTR(I))
        IF (HIGHORDERRAD) THEN
          IF (ML .LE. MM) THEN
            NR = NCS*(ML*(ML+1))/2 + ML+1
          ELSE
            NR = (NCS*MM+1)*ML-(MM*(2+NCS*(MM-1)))/2 + MM+1
          ENDIF
        ENDIF
        IR = IR + NR
        IF (IR .GT. MAXIR) THEN
          WRITE (6,*) 'RADIANCE_TRUNCATION: Really out of memory ',
     .      'for more radiance terms.'
          WRITE (6,*) 'Increase MAXIV.'
          STOP
        ENDIF
        RSHPTR(I+1) = IR
      ENDDO
      RSHPTR(NPTS+2) = IR
      RETURN
      END






      SUBROUTINE ACCELERATE_SOLUTION (ACCELPAR, NPTS, SHPTR, OSHPTR, 
     .                                SOURCE, DELSOURCE)
C       Accelerates the successive order of scattering series using the
C     input acceleration parameter (ACCELPAR).
      IMPLICIT NONE
      INTEGER NPTS, SHPTR(*), OSHPTR(*)
      REAL    ACCELPAR, SOURCE(*), DELSOURCE(*)
      INTEGER I, J, IS, ISD, NS, NSD, NSC

      IF (ACCELPAR .GT. 0.0) THEN
        DO I = 1, NPTS
          IS = SHPTR(I)
          NS = SHPTR(I+1)-IS
          ISD = OSHPTR(I)
          NSD = OSHPTR(I+1)-ISD
          NSC  = MIN(NS,NSD)
          DO J = 1, NSC
            SOURCE(IS+J) = SOURCE(IS+J) + ACCELPAR*DELSOURCE(ISD+J)
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END




      SUBROUTINE PATH_INTEGRATION (NX, NY, NZ, ML, MM, NLM, 
     .             NMU, NPHI0MAX, NPHI0, NANG, NCS, BCFLAG, IPFLAG,  
     .             NPTS, NCELLS,  XGRID, YGRID, ZGRID, GRIDPOS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, 
     .             SRCTYPE, SOLARMU, SOLARAZ, SKYRAD, WAVENO, WAVELEN, 
     .             UNITS, SFCTYPE, GNDTEMP, GNDALBEDO, 
     .             NXSFC, NYSFC, DELXSFC, DELYSFC, 
     .             NSFCPAR, SFCPARMS, SFCGRIDPARMS, 
     .             MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .             MU, PHI, WTDO,  
     .             CMU1, CPHI1, CMU2, CPHI2, WSAVE, FFTFLAG,
     .             DIRFLUX, FLUXES, EXTINCT, 
     .             SHPTR, SOURCE, RSHPTR, RADIANCE, 
     .             WORK, SWEEPORD, GRIDRAD, OLDNPTS, SP, STACK, 
     .       IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ)
C       Performs the path integrations through the medium specified by
C     the extinction (EXTINCT) and source function (SOURCE) in
C     spherical harmonic space.  The source function is transformed to
C     discrete ordinates one zenith angle at a time.  The integrations 
C     are done from the top boundary (with no reflection) down to the 
C     bottom boundary for the downwelling angles, and then up from the 
C     bottom to the top after computing the surface reflection/emission.  
C     After the path integrations produce the discrete angle radiance
C     it is transformed and added to output radiance spherical harmonic
C     expansion.  The hemispheric flux FLUXES is also computed.
C     The discrete ordinates in MU, PHI must be a complete set with the
C     downwelling (mu<0) angles first.
C       The surface reflection is handled differently for Lambertian
C     surfaces and more general surfaces specified by bidirectional 
C     reflection distribution functions (BRDF).  For the general BRDF,
C     the bottom boundary radiance must be computed for each upwelling 
C     ordinate, while for the Lambertian surface the boundary radiances 
c     can be computed just once, since they are isotropic.  For Lambertian 
C     cases the BCRAD array holds only the isotropic radiance for the 
C     boundaries at the top (downwelling) and bottom (upwelling) boundaries.
C     For the general BRDF, the BCRAD array in addition holds the 
C     downwelling radiances at all discrete ordinates for the bottom 
C     boundary points.  
      IMPLICIT NONE
      INTEGER NX, NY, NZ, ML, MM, NLM
Cf2py intent(in) :: NX, NY, NZ, ML, MM, NLM
      INTEGER NMU, NPHI0MAX, NPHI0(*), NANG, NCS
Cf2py intent(in) :: NMU, NPHI0MAX, NPHI0, NANG, NCS
      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) :: BCFLAG, IPFLAG
      INTEGER NPTS, NCELLS, MAXNBC, MAXBCRAD, NTOPPTS,NBOTPTS
Cf2py intent(in) :: NPTS, NCELLS, MAXNBC, MAXBCRAD
Cf2py intent(in,out) :: NTOPPTS,NBOTPTS
      INTEGER NXSFC, NYSFC, NSFCPAR
Cf2py intent(in,out) :: NXSFC, NYSFC, NSFCPAR
      INTEGER SHPTR(*), RSHPTR(*),  BCPTR(MAXNBC,2)
Cf2py intent(in,out) :: SHPTR, RSHPTR, BCPTR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in,out) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in,out) :: CELLFLAGS
      INTEGER SWEEPORD(NPTS,*)
Cf2py intent(in) :: SWEEPORD
      LOGICAL FFTFLAG(*)
Cf2py intent(in) :: FFTFLAG
      REAL    SOLARMU, SOLARAZ, SKYRAD
Cf2py intent(in) :: SOLARMU, SOLARAZ, SKYRAD
      REAL    GNDTEMP, GNDALBEDO
Cf2py intent(in) :: GNDTEMP, GNDALBEDO
      REAL    DELXSFC, DELYSFC
Cf2py intent(in) :: DELXSFC, DELYSFC
      REAL    SFCPARMS(NSFCPAR,*), SFCGRIDPARMS(NSFCPAR,*)
Cf2py intent(in,out) :: SFCPARMS, SFCGRIDPARMS
      REAL    WAVENO(2), WAVELEN
Cf2py intent(in) :: WAVENO, WAVELEN
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      REAL    MU(*), PHI(NMU,*), WTDO(NMU,*)
Cf2py intent(in) ::  MU, PHI, WTDO
      REAL    CMU1(*), CPHI1(*), CMU2(*), CPHI2(*), WSAVE(*)
Cf2py intent(in) :: CMU1, CPHI1, CMU2, CPHI2, WSAVE
      REAL    BCRAD(*)
Cf2py intent(in,out) :: BCRAD
      REAL    DIRFLUX(*), FLUXES(2,*), EXTINCT(*)
Cf2py intent(in,out) :: DIRFLUX, FLUXES
Cf2py intent(in) :: EXTINCT
      REAL    GRIDRAD(*), WORK(NPHI0MAX,NPTS), RADIANCE(*)
Cf2py intent(in) ::  GRIDRAD, WORK
Cf2py intent(in,out) :: RADIANCE
      REAL    SOURCE(*)
Cf2py intent(in,out) :: SOURCE
      CHARACTER SRCTYPE*1, UNITS*1, SFCTYPE*2
Cf2py intent(in) :: SRCTYPE, UNITS, SFCTYPE
      INTEGER OLDNPTS
Cf2py intent(in) :: OLDNPTS
      INTEGER I, K, NR, IPHI, IMU, IBC, IANG, IUPDOWN
      integer joct, ipt, ip, iorder, ipcell, icorner
      LOGICAL LAMBERTIAN, BTEST, theflag
      REAL    A, TRANSMIN
      INTEGER SP, STACK(50)
      INTEGER IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ
      
C         Set the minimum transmission for the cell tracing (1 for single cell)
      TRANSMIN = 1.00

C         Make the new grid cell/point sweeping order if need to
      IF (NPTS .NE. OLDNPTS) THEN
        CALL SWEEPING_ORDER (NX, NY, NZ, NPTS, NCELLS, 
     .                       GRIDPTR, TREEPTR, CELLFLAGS, 
     .                       BCFLAG, IPFLAG,  SWEEPORD,  GRIDRAD,
     .                       SP, STACK, IX, IY, IZ, SIX, SIY, SIZ,
     .                       EIX, EIY, EIZ, DIX, DIY, DIZ)
      ENDIF

      CALL FIND_BOUNDARY_POINTS (BCFLAG,IPFLAG, NPTS, SWEEPORD,GRIDPTR,
     .                       GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)


      LAMBERTIAN = SFCTYPE(2:2) .EQ. 'L'
      IF (NPTS .NE. OLDNPTS) THEN
C         Make the pointers to the grid points on the top and bottom boundaries
        CALL BOUNDARY_PNTS (NPTS, NANG, LAMBERTIAN, MAXNBC, MAXBCRAD, 
     .                      ZGRID(1), ZGRID(NZ), GRIDPOS,
     .                      NTOPPTS, NBOTPTS, BCPTR)
        IF (SFCTYPE(1:1) .EQ. 'V') THEN
C             If this is a variable surface then bilinearly interpolate 
C               the surface temperature and reflection parameters to the
C               bottom grid points. 
          CALL SURFACE_PARM_INTERP (NBOTPTS, BCPTR(1,2), GRIDPOS,
     .             SRCTYPE, WAVENO, WAVELEN, UNITS,
     .             NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,
     .             SFCGRIDPARMS)
        ENDIF
      ENDIF
C         Make the isotropic radiances for the top boundary
      CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN, 
     .                            UNITS, NTOPPTS, BCRAD(1))

      IF (SRCTYPE .NE. 'T') THEN
        CALL CHECK_SOLAR_SPECULAR_REFLECTION (SFCTYPE, 
     .               NMU, NPHI0MAX, NPHI0, MU, PHI, SOLARMU, SOLARAZ)
      ENDIF


C         Zero the flux arrays and set temporary radiance to -1 for not valid
      DO I = 1, NPTS
        FLUXES(1,I) = 0.0
        FLUXES(2,I) = 0.0
        GRIDRAD(I) = -1.0
      ENDDO
C         Zero the output spherical harmonic radiance
      NR = RSHPTR(NPTS+1)
      DO I = 1, NR
        RADIANCE(I) = 0.0
      ENDDO


C          Loop over the zenith angles (downwelling angles must be first)
      IANG = 1
      DO IMU = 1, NMU
        IF (IMU .EQ. NMU/2+1) THEN
C           For the first upwelling angle, make the bottom boundary radiances
C             for the Lambertian surfaces.  Compute the upwelling bottom 
C             radiances using the downwelling fluxes.
          IF (SFCTYPE .EQ. 'FL') THEN
            CALL FIXED_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2),
     .           DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO, 
     .           WAVENO, WAVELEN, UNITS, BCRAD(1+NTOPPTS))
          ELSE IF (SFCTYPE .EQ. 'VL') THEN
            CALL VARIABLE_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2),
     .           DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,
     .           BCRAD(1+NTOPPTS))
          ENDIF
        ENDIF
C           Transform the source to discrete ordinates for this zenith angle
        CALL SH_TO_DO (NPTS, ML,MM,NLM, NMU,NPHI0MAX,NPHI0(IMU), NCS,
     .         IMU, FFTFLAG, CMU1, CPHI1, WSAVE, SHPTR, SOURCE, WORK)

C           Enforce the discrete ordinate source function to be non-negative
        DO I = 1, NPTS
          DO IPHI = 1, NPHI0(IMU)
            WORK(IPHI,I) = MAX(0.0,WORK(IPHI,I))
          ENDDO
        ENDDO

C           Loop over the azimuthal angles
        DO IPHI = 1, NPHI0(IMU)
          IF (MU(IMU) .LT. 0.0) THEN
C               For downward ordinates, initialize the top boundary radiances
            IUPDOWN = 1
            DO IBC = 1, NTOPPTS
              I = BCPTR(IBC,1)
              GRIDRAD(I) = BCRAD(IBC)
            ENDDO
          ELSE
C               Initialize the bottom boundary radiances
            IUPDOWN = 2
            IF (.NOT. LAMBERTIAN) THEN
C                 If a Lambertian surface, use the radiances computed above,
C                   otherwise, compute the radiance for this angle
C                   by integrating over the stored downwelling radiances.
              CALL VARIABLE_BRDF_SURFACE (NBOTPTS,1,NBOTPTS,BCPTR(1,2),
     .               NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .               MU(IMU), PHI(IMU,IPHI), 
     .               SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX, 
     .               SFCTYPE, NSFCPAR, SFCGRIDPARMS, BCRAD(1+NTOPPTS))
            ENDIF
            DO IBC = 1, NBOTPTS
              I = BCPTR(IBC,2)
              GRIDRAD(I) = BCRAD(IBC+NTOPPTS)
            ENDDO
          ENDIF

C             Do backward path integration from each grid point not done yet
          IF (BTEST(IPFLAG,1) .AND. BTEST(IPFLAG,0)) THEN
            CALL BACK_INT_GRID1D (NX, NY, NZ, NPHI0MAX, NPTS, NCELLS, 
     .               GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .               GRIDPOS,SWEEPORD,MU(IMU),PHI(IMU,IPHI), TRANSMIN,
     .               EXTINCT, WORK, IPHI, GRIDRAD)
          ELSE IF (BTEST(IPFLAG,1) .AND. 
     .            .NOT. (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3))) THEN
            CALL BACK_INT_GRID2D (NX, NY, NZ, NPHI0MAX, NPTS, NCELLS, 
     .               GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, 
     .               GRIDPOS,SWEEPORD,MU(IMU),PHI(IMU,IPHI), TRANSMIN,
     .               EXTINCT, WORK, IPHI, GRIDRAD)
          ELSE
            CALL BACK_INT_GRID3D (NX, NY, NZ, NPHI0MAX, NPTS, NCELLS, 
     .               GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, 
     .               GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD, 
     .               MU(IMU),PHI(IMU,IPHI), TRANSMIN, BCFLAG, IPFLAG,
     .               EXTINCT, WORK, IPHI, GRIDRAD)
          ENDIF


C             Transfer the radiance for the (IMU,IPHI) discrete ordinate from 
C               GRIDRAD to the work array.
          DO I = 1, NPTS
            IF (GRIDRAD(I) .LT. -0.0001) THEN
              WRITE (6,*) 'PATH_INTEGRATION: Grid point has no value',
     .           I,MU(IMU),PHI(IMU,IPHI),GRIDRAD(I),XGRID(1),YGRID(1),
     .           GRIDPOS(1,I),GRIDPOS(2,I),GRIDPOS(3,I)
c              STOP
              GRIDRAD(I)=0.0
            ELSE
              WORK(IPHI,I) = GRIDRAD(I)
              GRIDRAD(I) = -1.0
            ENDIF
          ENDDO


C             Do the sum for the hemispheric fluxes for this angle
          A  = ABS(MU(IMU))*WTDO(IMU,IPHI)
          DO I = 1, NPTS
            FLUXES(IUPDOWN,I) = FLUXES(IUPDOWN,I) + A*WORK(IPHI,I)
          ENDDO

C             For a general BRDF surface save the downwelling radiance
C               at the bottom boundary points.
          IF (.NOT. LAMBERTIAN .AND. IMU .LE. NMU/2) THEN
            K = NTOPPTS + NBOTPTS*IANG
            DO IBC = 1, NBOTPTS
              I = BCPTR(IBC,2)
              BCRAD(IBC+K) = WORK(IPHI,I)
            ENDDO
          ENDIF

          IANG = IANG + 1
        ENDDO

C           Transform the radiance to spherical harmonics for this zenith angle
        CALL DO_TO_SH (NPTS, ML,MM, NLM, NMU,NPHI0MAX,NPHI0(IMU), NCS,
     .         IMU, FFTFLAG, CMU2,CPHI2, WSAVE, RSHPTR, WORK, RADIANCE)

      ENDDO

      OLDNPTS = NPTS

      RETURN
      END
 




      SUBROUTINE BOUNDARY_PNTS (NPTS,NANG, LAMBERTIAN, MAXNBC,MAXBCRAD,
     .                          ZBOT, ZTOP, GRIDPOS,
     .                          NTOPPTS, NBOTPTS, BCPTR)
C       Returns the BCPTR pointer array for the NTOPPTS top boundary
C     points (BCPTR(*,1)) and the NBOTPTS bottom boundary points (BCPTR(*,2)).
C     All the grid points are looked at to find the boundary ones.  
      IMPLICIT NONE
      INTEGER NPTS, NANG, MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS
      INTEGER BCPTR(MAXNBC,2)
      LOGICAL LAMBERTIAN
      REAL    ZBOT, ZTOP, GRIDPOS(3,NPTS)
      INTEGER I, IT, IB, NA

C         Loop over all points finding the top and bottom ones
      IF (LAMBERTIAN) THEN
        NA = 1
      ELSE
        NA = NANG/2+1
      ENDIF
      IT = 0
      IB = 0
      DO I = 1, NPTS
        IF (GRIDPOS(3,I) .GE. ZTOP) THEN
          IT = IT + 1
          IF (IT .GT. MAXNBC .OR. IT .GT. MAXBCRAD) THEN
            WRITE (6,*) 'BOUNDARY_PNTS: MAXNBC exceeded - top.'
            STOP
          ENDIF
          BCPTR(IT,1) = I
        ENDIF
        IF (GRIDPOS(3,I) .LE. ZBOT) THEN
          IB = IB + 1
          IF (IB .GT. MAXNBC .OR. IT+IB*NA .GT. MAXBCRAD) THEN
            WRITE (6,*) 'BOUNDARY_PNTS: MAXNBC exceeded - bottom.'
            STOP 
          ENDIF
          BCPTR(IB,2) = I
        ENDIF
      ENDDO
      NTOPPTS = IT
      NBOTPTS = IB
      RETURN
      END





      SUBROUTINE SURFACE_PARM_INTERP (NBOTPTS, BCPTR, GRIDPOS,
     .               SRCTYPE, WAVENO, WAVELEN, UNITS,
     .               NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,
     .               SFCGRIDPARMS)
C       Bilinearly interpolates the surface properties from the regular
C     surface grid to the boundary points of the adaptive grid.
C     The locations of the NBOTPTS grid points on the bottom boundary 
C     (pointed to by BCPTR) are in the GRIDPOS array.  The input regular
C     surface info is specified by the number of grid lines (NXSFC,NYSFC),
C     the grid spacing (DELXSFC,DELYSFC), the number of parameters at
C     each grid point (NSFCPAR), and the surface parameters (SFCPARMS).
C     The output interpolated surface parameters at the NBOTPTS boundary
C     points are in SFCGRIDPARMS.  The first parameter is input as
C     temperature, but converted to Planck function for output.
      IMPLICIT NONE
      INTEGER NBOTPTS, BCPTR(*),  NXSFC, NYSFC, NSFCPAR
      REAL    WAVENO(2), WAVELEN
      REAL    DELXSFC, DELYSFC
      REAL    GRIDPOS(3,*), SFCPARMS(NSFCPAR,NXSFC+1,NYSFC+1)
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS)
      CHARACTER  SRCTYPE*1, UNITS*1
      INTEGER I, IBC, IX, IY, J
      REAL    RX, RY, U, V, TEMP, PLANCK

C         Loop over all bottom points
      DO IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        RX = GRIDPOS(1,I)/DELXSFC
        RY = GRIDPOS(2,I)/DELYSFC
        IX = MAX(1,MIN(NXSFC,INT(RX)+1))
        IY = MAX(1,MIN(NYSFC,INT(RY)+1))
        U = MAX(0.0,MIN(1.0,RX-(IX-1)))
        V = MAX(0.0,MIN(1.0,RY-(IY-1)))
        DO J = 1, NSFCPAR
          SFCGRIDPARMS(J,IBC) = (1-U)*(1-V)*SFCPARMS(J,IX,IY)
     .        + (1-U)*V*SFCPARMS(J,IX,IY+1)
     .        + U*(1-V)*SFCPARMS(J,IX+1,IY) + U*V*SFCPARMS(J,IX+1,IY+1)
        ENDDO
        IF (SRCTYPE .EQ. 'S') THEN
          PLANCK = 0.0
        ELSE
          TEMP = SFCGRIDPARMS(1,IBC)
          CALL PLANCK_FUNCTION (TEMP, UNITS, WAVENO, WAVELEN, PLANCK)
        ENDIF
        SFCGRIDPARMS(1,IBC) = PLANCK
      ENDDO
      RETURN
      END



      SUBROUTINE COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD,WAVENO,WAVELEN,
     .                                  UNITS, NTOPPTS, BCRAD)
C       Returns the beginning of the BCRAD array with the downwelling 
C     radiance for the NTOPPTS top boundary points.  Currently all
C     points have the same isotropic radiance.
      IMPLICIT NONE
      INTEGER NTOPPTS
Cf2py intent(in) :: NTOPPTS
      REAL    SKYRAD, WAVENO(2), WAVELEN
Cf2py intent(in) :: SKYRAD, WAVENO, WAVELEN
      REAL    BCRAD(*)
Cf2py intent(in, out) :: BCRAD
      CHARACTER  SRCTYPE*1, UNITS*1
Cf2py intent(in) :: SRCTYPE, UNITS
      INTEGER IBC
      REAL    SKYRAD2

C           At top, boundary radiance is any isotropic sky radiation
      IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
        SKYRAD2 = SKYRAD
      ELSE
        CALL PLANCK_FUNCTION (SKYRAD, UNITS, WAVENO, WAVELEN, SKYRAD2)
      ENDIF
C         Loop over all points assigning the uniform radiance
      DO IBC = 1, NTOPPTS
        BCRAD(IBC) = SKYRAD2
      ENDDO
      RETURN
      END
 



      SUBROUTINE FIXED_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR,
     .             DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO, 
     .             WAVENO, WAVELEN, UNITS,   BCRAD)
C       Returns the BCRAD array with the upwelling radiance for a 
C     Lambertian surface for the bottom points (pointed to by BCPTR).  
C     If there is reflection the downwelling hemispheric flux at the 
C     bottom points is gotten (from FLUXES). There may also be thermal 
C     emission (depending on temperature and emissivity) or reflected 
C     collimated solar radiation (DIRFLUX).
      IMPLICIT NONE
      INTEGER NBOTPTS, BCPTR(*)
Cf2py intent(in) :: NBOTPTS, BCPTR
      REAL    DIRFLUX(*), FLUXES(2,*)
Cf2py intent(in) :: DIRFLUX, FLUXES
      REAL    GNDTEMP, GNDALBEDO, WAVENO(2), WAVELEN
Cf2py intent(in) ::  GNDTEMP, GNDALBEDO, WAVENO, WAVELEN
      REAL    BCRAD(*)
Cf2py intent(in, out) :: BCRAD
      CHARACTER  SRCTYPE*1, UNITS*1
Cf2py intent(in) :: SRCTYPE, UNITS
      INTEGER I, IBC
      REAL    ALB, GNDRAD

      ALB = GNDALBEDO/ACOS(-1.0)
      IF (SRCTYPE .EQ. 'T' .OR. SRCTYPE .EQ. 'B') THEN
        CALL PLANCK_FUNCTION (GNDTEMP, UNITS, WAVENO, WAVELEN, GNDRAD)
        GNDRAD = GNDRAD*(1.0-GNDALBEDO)
      ENDIF

C         Loop over all bottom points finding the upwelling radiance
      DO IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        IF (SRCTYPE .EQ. 'S') THEN
          BCRAD(IBC) = ALB*(DIRFLUX(I) + FLUXES(1,I))
        ELSE IF (SRCTYPE .EQ. 'T') THEN
          BCRAD(IBC) = GNDRAD + ALB*FLUXES(1,I)
        ELSE IF (SRCTYPE .EQ. 'B') THEN
          BCRAD(IBC) = ALB*(DIRFLUX(I) + FLUXES(1,I)) + GNDRAD
        ENDIF
        
      ENDDO
      RETURN
      END




      SUBROUTINE VARIABLE_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR,
     .             DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,
     .             BCRAD)
C       Returns the BCRAD array with the upwelling radiance for a 
C     variable Lambertian surface for the bottom points (pointed to by BCPTR).  
C     SFCGRIDPARMS contains the interpolated surface properties for the
C     bottom boundary points; the first is Planck function and the second is
C     albedo.  The reflection is computed from the downwelling hemispheric 
C     flux (in FLUXES).  There may also be thermal emission (depending on
C     temperature and emissivity) or reflected collimated solar radiation 
C     (from DIRFLUX).
      IMPLICIT NONE
      INTEGER NBOTPTS, BCPTR(*), NSFCPAR
Cf2py intent(in) :: NBOTPTS, BCPTR, NSFCPAR
      REAL    DIRFLUX(*), FLUXES(2,*)
Cf2py intent(in) :: DIRFLUX, FLUXES
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS)
Cf2py intent(in) :: SFCGRIDPARMS
      REAL    BCRAD(*)
Cf2py intent(in, out) :: BCRAD
      CHARACTER  SRCTYPE*1
Cf2py intent(in) :: SRCTYPE
      INTEGER I, IBC
      REAL    OPI, ALB, GNDRAD

      OPI = 1.0/ACOS(-1.0)

C         Loop over all bottom points
      DO IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        ALB = SFCGRIDPARMS(2,IBC)
        IF (SRCTYPE .EQ. 'S') THEN
          BCRAD(IBC) = OPI*ALB*(DIRFLUX(I) + FLUXES(1,I))
        ELSE
          GNDRAD = SFCGRIDPARMS(1,IBC)*(1-ALB)
          IF (SRCTYPE .EQ. 'T') THEN
            BCRAD(IBC) = GNDRAD + OPI*ALB*FLUXES(1,I)
          ELSE IF (SRCTYPE .EQ. 'B') THEN
            BCRAD(IBC) = OPI*ALB*(DIRFLUX(I) + FLUXES(1,I)) + GNDRAD
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END




      SUBROUTINE VARIABLE_BRDF_SURFACE (NBOTPTS, IBEG, IEND, BCPTR, 
     .               NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,
     .               SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,
     .               SFCTYPE, NSFCPAR, SFCGRIDPARMS, BCRAD)
C       Computes the upwelling radiance of the bottom boundary points
C     for one outgoing direction using the specified bidirectional 
C     reflectance distribution function.  The upwelling radiance
C     includes the reflection of the incident radiance, the thermal 
C     emission (emissivity*Planck function), and the reflected direct
C     solar flux (if applicable).  The upwelling radiance is the 
C     integral over all incident directions of the BRDF times the 
C     downwelling radiance, so a discrete sum is done and the integral 
C     weights (WTDO) are included. The general BRDF function is called 
C     to compute the reflectance for the particular BRDF type (SFCTYPE) 
C     with parameters (SFCGRIDPARMS) for the incident and outgoing 
C     directions.  The emissivity is computed implicitly from the integral
C     of the BRDF.  The outgoing direction is specified with (MU2,PHI2),
C     and the BRDF is computed for all incident directions (loop over 
C     JMU,JPHI).  The incident downwelling radiances are input in 
C     BCRAD(*,2...) and the outgoing upwelling radiances are output 
C     in BCRAD(*,1).  There is a special case for a specularly reflecting
C     surface, because the discrete ordinate integration cannot resolve 
C     the delta function in the BRDF.
      IMPLICIT NONE
      INTEGER NBOTPTS, IBEG, IEND, BCPTR(NBOTPTS)
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    MU2, PHI2, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(NBOTPTS,*)
      CHARACTER SRCTYPE*1, SFCTYPE*2
      INTEGER JMU, JPHI, IBC, I, JANG
      LOGICAL SPECULAR_SURFACE
      REAL    OPI, REFLECT, W, SUM0, SUM1, RADSPEC

      OPI = 1.0/ACOS(-1.0)

      DO IBC = IBEG, IEND

C         Initialize the upwelling boundary radiances to zero or to 
C           the reflected direct solar flux.
        IF (SRCTYPE .EQ. 'T') THEN
          BCRAD(IBC,1) = 0.0
        ELSE
          I = BCPTR(IBC)
          CALL SURFACE_BRDF (SFCTYPE(2:2), SFCGRIDPARMS(2,IBC),WAVELEN,
     .                       MU2, PHI2, SOLARMU,SOLARAZ,  REFLECT)
          BCRAD(IBC,1) = OPI*REFLECT*DIRFLUX(I)
        ENDIF

C         A specularly reflecting surface is a special case: instead
C           of the integral over incident radiances, the radiance from
C           the incident specular direction is interpolated.
        IF (SPECULAR_SURFACE(SFCTYPE(2:2))) THEN
C             Get the reflection coefficient for this outgoing direction
          CALL SURFACE_BRDF (SFCTYPE(2:2), SFCGRIDPARMS(2,IBC),WAVELEN,
     .                       MU2, PHI2, -MU2, PHI2, REFLECT)
C             Interpolate by weighting with the inverse cubed angular 
C               distance between discrete ordinate and incident directions.
          JANG = 1
          SUM0 = 0.0
          SUM1 = 0.0
          DO JMU = 1, NMU/2
            DO JPHI = 1, NPHI0(JMU)
              W = ACOS(-MU(JMU)*MU2+SQRT((1-MU(JMU)**2)*(1-MU2**2))
     .                 *COS(PHI2-PHI(JMU,JPHI)))
              W = 1.0/MAX(W**3,1.0E-8)
              SUM0 = SUM0 + W
              SUM1 = SUM1 + W*BCRAD(IBC,JANG+1)
              JANG = JANG + 1
            ENDDO
          ENDDO
          RADSPEC = SUM1/SUM0
          BCRAD(IBC,1) = BCRAD(IBC,1) + REFLECT*RADSPEC
     .                   + (1-REFLECT)*SFCGRIDPARMS(1,IBC)

        ELSE
C           For a nonspecular surface integrate over the incident discrete
C             ordinate directions (JMU,JPHI)
          JANG = 1
          DO JMU = 1, NMU/2
            DO JPHI = 1, NPHI0(JMU)
              CALL SURFACE_BRDF (SFCTYPE(2:2), SFCGRIDPARMS(2,IBC),
     .               WAVELEN, MU2,PHI2, MU(JMU),PHI(JMU,JPHI), REFLECT)
              W = OPI*ABS(MU(JMU))*WTDO(JMU,JPHI)
              BCRAD(IBC,1) = BCRAD(IBC,1) + W*REFLECT*BCRAD(IBC,JANG+1) 
     .                     + W*(1-REFLECT)*SFCGRIDPARMS(1,IBC)
              JANG = JANG + 1
            ENDDO
          ENDDO
        ENDIF

      ENDDO
      RETURN
      END


      SUBROUTINE CHECK_SOLAR_SPECULAR_REFLECTION (SFCTYPE, 
     .               NMU, NPHI0MAX, NPHI0, MU, PHI, SOLARMU, SOLARAZ)
C       For solar cases with specularly reflecting surface, check that
C     the reflected solar direction corresponds to a discrete ordinate 
C     so that the reflected solar beam is not lost.
      IMPLICIT NONE
      INTEGER NMU, NPHI0MAX, NPHI0(NMU)
      REAL    MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    SOLARMU, SOLARAZ
      CHARACTER SFCTYPE*2
      INTEGER JMU, JPHI
      LOGICAL SPECOK, SPECULAR_SURFACE

      IF (SPECULAR_SURFACE(SFCTYPE(2:2))) THEN
        SPECOK = .FALSE.
        DO JMU = NMU/2+1, NMU
          DO JPHI = 1, NPHI0(JMU)
            IF (ABS(-SOLARMU-MU(JMU)) .LT. 0.01 .AND. 
     .          ABS(SOLARAZ-PHI(JMU,JPHI)) .LT. 0.01)  SPECOK=.TRUE.
          ENDDO
        ENDDO
        IF (.NOT. SPECOK) THEN
          WRITE (6,*) 'Specularly reflecting surface must have solar'
          WRITE (6,*) 'reflected direction along a discrete ordinate.'
          STOP
        ENDIF
      ENDIF
      RETURN
      END





      SUBROUTINE SH_TO_DO (NPTS, ML,MM,NLM, NMU,NPHI0MAX,NPHI0, NCS, 
     .                     IMU, FFTFLAG, CMU1, CPHI1, WSAVE, 
     .                     SHPTR, INDATA, OUTDATA)
C       Transforms the input data from spherical harmonic space to discrete 
C     ordinate space for one zenith angle (IMU).  Successive transforms are 
C     done in zenith and then azimuth angle.  The Fourier transform in 
C     azimuth is done using a matrix multiply DFT for smaller number of 
C     angles (NPHI0*(NCS*MM+1)<160) and with an FFT for larger number of 
C     angles.  The FFT is used for any number of angles above the limit, 
C     but is much more efficient for NPHI0 whose factors are small primes
C     (powers of 2 are best).  The spherical harmonic truncation of the
C     input array is adaptive (with pointer SHPTR); the last term for 
C     each grid point must be the maximum m for the particular l term.
C     The number of floating point operations for the DFT is 2*Nmu*Nphi0*Nm.
C     The number of floating point operations for the zenith angle transform
C     is 2*Nmu*Nlm.
      IMPLICIT NONE
      INTEGER NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, NCS
      INTEGER IMU, SHPTR(NPTS+1)
      LOGICAL FFTFLAG(NMU)
      REAL    INDATA(*), OUTDATA(NPHI0MAX,NPTS)
      REAL    CMU1(NLM,NMU), CPHI1(-16:16,32,NMU)
      REAL    WSAVE(3*NPHI0MAX+15,NMU)
      INTEGER I, J, K, L, M, N, IPHI, MS, ME, MR, IS, NSH
      REAL    SUM1

      INTEGER, ALLOCATABLE :: MOFJ(:)
      REAL, ALLOCATABLE :: SUM(:)

      ALLOCATE (MOFJ(NLM), SUM(-MM:MM))

       
C         Make the M for each J array
      J = 0
      DO L = 0, ML
        DO M = (1-NCS)*MIN(L,MM), MIN(L,MM)
          J = J + 1
          MOFJ(J) = M
        ENDDO
      ENDDO

C         Do the transform for each grid point
      DO I = 1, NPTS
C             Determine the current SH truncation to use
        IS = SHPTR(I)
        NSH = SHPTR(I+1)-IS
C             Special case for isotropic field
        IF (NSH .EQ. 0) THEN
          DO K = 1, NPHI0
            OUTDATA(K,I) = 0.0
          ENDDO
        ELSE IF (NSH .EQ. 1) THEN
          SUM1 = CMU1(1,1)*INDATA(IS+1)
          DO K = 1, NPHI0
            OUTDATA(K,I) = SUM1
          ENDDO
        ELSE
C             Figure the max Fourier azimuthal mode we can do for this Nphi0
          IF (NCS .EQ. 1) THEN
            ME = MAX(0,MIN(NPHI0-2,MOFJ(NSH)))
            MS = 0
          ELSE
            ME = MAX(0,MIN((NPHI0/2)-1,MOFJ(NSH)))
            MS = -ME
          ENDIF
C             First do Legendre transform by summing over l for each m.
          DO M = MS, ME
            SUM(M) = 0.0
          ENDDO
          DO J = 1, NSH
            M = MOFJ(J)
            SUM(M) = SUM(M) + CMU1(J,IMU)*INDATA(IS+J)
          ENDDO
C             Then do the Fourier transform from m to phi for each mu
          IF (FFTFLAG(IMU)) THEN
            IF (NCS .EQ. 1) THEN
C               For cosine modes only use a cosine FFT
              OUTDATA(1,I) = SUM(0)
              DO M = 1, ME
                OUTDATA(M+1,I) = 0.5*SUM(M)
              ENDDO
              DO N = ME+2, NPHI0
                OUTDATA(N,I) = 0.0
              ENDDO
              CALL COST (NPHI0,OUTDATA(1,I),WSAVE(1,IMU))
            ELSE
C                 For sines and cosines use a real FFT, and rearrange terms
              OUTDATA(1,I) = SUM(0)
              MR = 2
              DO M = 1, ME
                OUTDATA(MR,I) = 0.5*SUM(M)
                OUTDATA(MR+1,I) = -0.5*SUM(-M)
                MR = MR + 2
              ENDDO
              DO N = MR, NPHI0
                OUTDATA(N,I) = 0.0
              ENDDO
              CALL RFFTB (NPHI0,OUTDATA(1,I),WSAVE(1,IMU))
            ENDIF
          ELSE
C               Else do a slow DFT
            K = 0
            DO IPHI = 1, NPHI0
              SUM1 = 0.0
              DO M = MS, ME
                SUM1 = SUM1 + CPHI1(M,IPHI,IMU)*SUM(M)
              ENDDO
              K = K + 1
              OUTDATA(K,I) = SUM1
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END




      SUBROUTINE DO_TO_SH (NPTS, ML,MM,NLM, NMU,NPHI0MAX,NPHI0, NCS, 
     .             IMU, FFTFLAG, CMU2, CPHI2, WSAVE, 
     .             RSHPTR, INDATA, OUTDATA)
C       Transforms the input field from discrete ordinate space to spherical 
C     harmonic space for one zenith angle (IMU) by keeping a running
C     sum in OUTDATA between successive calls.  (See SH_TO_DO for more).
      IMPLICIT NONE
      INTEGER NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, NCS
      INTEGER IMU, RSHPTR(NPTS+1)
      LOGICAL FFTFLAG(NMU)
      REAL    INDATA(NPHI0MAX,NPTS), OUTDATA(*)
      REAL    CMU2(NMU,NLM), CPHI2(32,-16:16,NMU)
      REAL    WSAVE(3*NPHI0MAX+15,NMU)
      INTEGER I, IPHI, J, K, L, M, MS, ME, MR, IS, NSH
      REAL    SUM1, DELPHI

      INTEGER, ALLOCATABLE :: MOFJ(:)
      REAL, ALLOCATABLE :: SUM(:)
      ALLOCATE (MOFJ(NLM), SUM(-MM:MM))
      IF (NCS .EQ. 1) THEN
        DELPHI = ACOS(-1.0)/MAX(1,NPHI0-1)
      ELSE
        DELPHI = 2.0*ACOS(-1.0)/NPHI0
      ENDIF
C         Make the M for each J array
      J = 0
      DO L = 0, ML
        DO M = (1-NCS)*MIN(L,MM), MIN(L,MM)
          J = J + 1
          MOFJ(J) = M
        ENDDO
      ENDDO

C         Do the transform for each grid point
      IS = 0
      DO I = 1, NPTS
C             Determine the current SH truncation to use
        NSH = RSHPTR(I+1)-RSHPTR(I)
        ME = MOFJ(NSH)
        DO M = (1-NCS)*ME, ME
          SUM(M) = 0.0
        ENDDO
C           Figure the max Fourier azimuthal mode we can do for this Nphi0
        IF (NCS .EQ. 1) THEN
          ME = MAX(0,MIN(NPHI0-2,MOFJ(NSH)))
          MS = 0
        ELSE
          ME = MAX(0,MIN((NPHI0/2)-1,MOFJ(NSH)))
          MS = -ME
        ENDIF
C           First do Fourier transform from phi to m for each mu
        IF (FFTFLAG(IMU)) THEN
          IF (NCS .EQ. 1) THEN
C               For cosine modes only use a cosine FFT
            CALL COST (NPHI0,INDATA(1,I),WSAVE(1,IMU))
            DO M = 0, ME
              SUM(M) = INDATA(M+1,I)*DELPHI
            ENDDO
          ELSE
C                 For sines and cosines use a real FFT, and rearrange terms
            CALL RFFTF (NPHI0,INDATA(1,I),WSAVE(1,IMU))
            SUM(0) = INDATA(1,I)*DELPHI
            MR = 2
            DO M = 1, ME
              SUM(M) = INDATA(MR,I)*DELPHI
              SUM(-M) = -INDATA(MR+1,I)*DELPHI
              MR = MR + 2
            ENDDO
          ENDIF
        ELSE
C               Else do a slow DFT
          DO M = MS, ME
            SUM1 = 0.0
            K = 0
            DO IPHI = 1, NPHI0
              K = K + 1
              SUM1 = SUM1 + CPHI2(IPHI,M,IMU)*INDATA(K,I)
            ENDDO
            SUM(M) = SUM1
          ENDDO
        ENDIF
C             Then do Legendre transform by adding to output for each l, m
        DO J = 1, NSH
          M = MOFJ(J)
          OUTDATA(IS+J) = OUTDATA(IS+J) + CMU2(IMU,J)*SUM(M)
        ENDDO
        IS = IS + NSH
      ENDDO
      RETURN
      END






      SUBROUTINE SWEEPING_ORDER (NX, NY, NZ, NPTS, NCELLS, 
     .                          GRIDPTR, TREEPTR, CELLFLAGS, 
     .                          BCFLAG, IPFLAG,  SWEEPORD, GRIDRAD,
     .       SP,STACK,IX,IY,IZ,SIX,SIY,SIZ,EIX,EIY,EIZ,DIX,DIY,DIZ)
C       Make the array of grid point sweeping order for each discrete
C     ordinate octant (SWEEPORD).  Rather than pointing to the grid points
C     directly, the elements of SWEEPORD point to the cell (bits 3-30) and 
C     the corner of the cell (bits 0-2).  The first index is the grid
C     point visiting order (1-NPTS), and the second index refers to the
C     octant, but is not the octant number.  The octant index is translated
C     to the octant number with IOCTORDER - this complication arises because
C     of the desire to use just 2 or 4 octants in 1D or 2D.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS, BCFLAG, IPFLAG
      INTEGER GRIDPTR(8,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SWEEPORD(NPTS,*)
      INTEGER*2 CELLFLAGS(NCELLS)
      REAL    GRIDRAD(NPTS)
      INTEGER NOCT, IOCT, JOCT, NXC, NYC
      INTEGER IPCELL, IPT, IORDER
      INTEGER SP, STACK(50)
      INTEGER INDEXCORN, ICORNER, CORNDOG(8,8), IOCTORDER(8)
      LOGICAL SWEEP_NEXT_CELL, BTEST
      INTEGER IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ
      DATA CORNDOG/8,7,6,5,4,3,2,1, 7,8,5,6,3,4,1,2, 6,5,8,7,2,1,4,3,
     .             5,6,7,8,1,2,3,4, 4,3,2,1,8,7,6,5, 3,4,1,2,7,8,5,6,
     .             2,1,4,3,6,5,8,7, 1,2,3,4,5,6,7,8/
      DATA IOCTORDER/1,5,2,6,3,7,4,8/

C         IOCT is the octant the ray come from or the discrete ordinate goes to.
C      IOCT = 1 + BITX + 2*BITY + 4*BITZ

      IF (BTEST(IPFLAG,1) .AND. BTEST(IPFLAG,0)) THEN
        NOCT = 2
      ELSE IF (BTEST(IPFLAG,1)) THEN
        NOCT = 4
      ELSE
        NOCT = 8
      ENDIF
      NXC = NX
      IF (BTEST(BCFLAG,0)) NXC = NX+1
      IF (BTEST(BCFLAG,2) .AND. .NOT. BTEST(IPFLAG,0)) NXC = NX-1
      NYC = NY
      IF (BTEST(BCFLAG,1)) NYC = NY+1
      IF (BTEST(BCFLAG,3) .AND. .NOT. BTEST(IPFLAG,1)) NYC = NY-1

C         Loop over the number of octants needed; get the octant number (1-8)
      DO JOCT = 1, NOCT
        IOCT = IOCTORDER(JOCT)

C           Set each grid point flag to indicate unvisited.
        DO IPT = 1, NPTS
          GRIDRAD(IPT) = -1.0
        ENDDO

C           Sweep through all the grid points
        IORDER = 1
        IPCELL = 0
        INDEXCORN = 8
100     CONTINUE
C         INDEXCORN is a logical ordering index to the 8 corners in a cell;
C           it is not the actual corner. ICORNER is the actual corner and 
C           depends on the ray direction (IOCT) and INDEXCORN.
C           Scan the 8 points in a cell according to the ray direction:
C           Do the ones furthest "upstream" first; loop until we find
C           a point that is not already done.
            IF (INDEXCORN .EQ. 8) THEN
C               If done with the current cell, then go to next
              INDEXCORN = 1
              IF (.NOT. SWEEP_NEXT_CELL (BCFLAG, NXC, NYC, NZ, IOCT, 
     .           IPCELL, TREEPTR, CELLFLAGS, SP, STACK, IX, IY, IZ,
     .       SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ))  GOTO 900
            ELSE
              INDEXCORN = INDEXCORN + 1
            ENDIF
C               Compute the corresponding point
C             ICORNER = IAND(NOT(IEOR(INDEXCORN-1,IOCT-1)),7) + 1
            ICORNER = CORNDOG(INDEXCORN,IOCT)
            IPT = GRIDPTR(ICORNER,IPCELL)
          IF (GRIDRAD(IPT) .GE. 0.0)  GOTO 100

          GRIDRAD(IPT) = 1.0
C             Fill in the next entry with the grid cell and corner of this point
          SWEEPORD(IORDER,JOCT) = IOR(ISHFT(IPCELL,3),ICORNER-1)
          IORDER = IORDER + 1
        GOTO 100
900     CONTINUE
      ENDDO

      RETURN
      END






      SUBROUTINE BACK_INT_GRID3D (NX,NY,NZ, NA, NPTS, NCELLS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD, 
     .             MU, PHI, TRANSMIN, BCFLAG, IPFLAG,
     .             EXTINCT, SOURCE, KANG, GRIDRAD)
C       Sweeps through the spatial grid computing radiances for the
C     discrete ordinate direction (MU,PHI) by integrating the source
C     function and extinction backward from each point.  The sweeping
C     order through the grid points is in SWEEPORD. If a point is
C     already done (nonnegative in GRIDRAD) then it is skipped.
C     For each grid point a ray opposite to the ordinate direction is traced
C     back to a face containing known (valid) radiances.  As each cell is
C     crossed (usually only one) the integration of the source function
C     across the cell is done and added to the radiance.  The radiance
C     at the end of the ray path is interpolated from the known grid point
C     radiances on the face.  The resulting radiances are put in GRIDRAD.
      INTEGER NX, NY, NZ, NA, NPTS, NCELLS, KANG, BCFLAG, IPFLAG
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SWEEPORD(NPTS,*)
      INTEGER*2 CELLFLAGS(NCELLS)
      REAL    GRIDPOS(3,NPTS), XGRID(*), YGRID(*), ZGRID(NZ)
      REAL    MU, PHI, TRANSMIN
      REAL    EXTINCT(NPTS), SOURCE(NA,NPTS), GRIDRAD(NPTS)
      INTEGER BITX, BITY, BITZ, IOCT, JOCT
      INTEGER IPCELL, ICELL, INEXTCELL, IPT, IFACE
      INTEGER I1, I2, I3, I4, IOPP, JFACE, KFACE, IC, IZ
      INTEGER GRIDFACE(4,6), OPPFACE(6)
      INTEGER IORDER, ICORNER, JOCTORDER3(8), JOCTORDER2(8)
      LOGICAL VALIDRAD, VALIDFACE, IPINX, IPINY, BTEST
      DOUBLE PRECISION PI, CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, path
      DOUBLE PRECISION U, V, F1, F2, F3, F4
      DOUBLE PRECISION EXT, EXT0, EXT1, SRC, SRCEXT0, SRCEXT1
      DOUBLE PRECISION EXT0P, SRCEXT0P
      DOUBLE PRECISION TAU, TRANSCELL, ABSCELL, RAD, RAD0, TRANSMIT
      DATA GRIDFACE/1,3,5,7, 2,4,6,8,  1,2,5,6, 3,4,7,8,
     .              1,2,3,4, 5,6,7,8/, OPPFACE/2,1,4,3,6,5/
      DATA JOCTORDER3/1,3,5,7,2,4,6,8/, JOCTORDER2/1,3,1,3,2,4,2,4/


      EPS = 1.0E-3*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))

C         Make the ray direction (opposite to the discrete ordinate direction)
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

C         Setup the sweeping direction and the gridpoint location
C         BITc is 0 for positive direction trace back ray, 1 for negative.
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
        STOP 'BACK_INT_GRID: Bad MU'
      ENDIF
C         IOCT is the octant the ray come from or the discrete ordinate goes to.
      IOCT = 1 + BITX + 2*BITY + 4*BITZ
      IF (BTEST(IPFLAG,1)) THEN
        JOCT = JOCTORDER2(IOCT)
      ELSE
        JOCT = JOCTORDER3(IOCT)
      ENDIF

      IZ = (NZ-1)*(1-BITZ)+1
C         Sweep through all the grid points
      DO IORDER = 1, NPTS
        IPCELL = ISHFT(SWEEPORD(IORDER,JOCT),-3)
        ICORNER = IBITS(SWEEPORD(IORDER,JOCT),0,3)+1
        IPT = GRIDPTR(ICORNER,IPCELL)
        
        IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3)) THEN
C           For multiple processors get the subdomain boundary radiances
C           for each Z slab
          IF ((GRIDPOS(3,IPT) .LT. ZGRID(IZ) .AND. BITZ.EQ.0) .OR.
     .        (GRIDPOS(3,IPT) .GT. ZGRID(IZ) .AND. BITZ.EQ.1)) THEN
            IZ = IZ + 2*BITZ-1
            CALL CALC_BOUNDARY_RADIANCES (BCFLAG, IPFLAG, JOCT, IZ,
     .                        NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                        NA, NPTS, NCELLS, GRIDPTR, 
     .                        NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS,
     .                        MU, PHI, EXTINCT, SOURCE, 
     .                        KANG, GRIDRAD)
          ENDIF
        ENDIF
        IF (GRIDRAD(IPT) .GE. 0.0)  GOTO 500

C           Start off at the grid point
        ICELL = IPCELL
        RAD = 0.0
        TRANSMIT = 1.0
        EXT1 = EXTINCT(IPT)
        SRCEXT1 = EXT1*SOURCE(KANG,IPT)
        XE = GRIDPOS(1,IPT)
        YE = GRIDPOS(2,IPT)
        ZE = GRIDPOS(3,IPT)

C           Loop until finding a face with known radiances
        VALIDRAD = .FALSE.
        DO WHILE (.NOT. VALIDRAD)
C             Make sure current cell is valid
          IF (ICELL .LE. 0) THEN
            WRITE (6,*) 'BACK_INT_GRID: ICELL=0  ', 
     .           MU, PHI, IPT, INEXTCELL,xe,ye,ze,
     .           gridpos(1,ipt),gridpos(2,ipt),gridpos(3,ipt)
            STOP
          ENDIF
          IPINX = BTEST(INT(CELLFLAGS(ICELL)),0)
          IPINY = BTEST(INT(CELLFLAGS(ICELL)),1)

C             Find boundaries of the current cell
C             Find the three possible intersection planes (X,Y,Z)
C               from the coordinates of the opposite corner grid point
          IOPP = GRIDPTR(9-IOCT,ICELL)
C             Get the distances to the 3 planes and select the closest
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
            WRITE (6,*) 'BACK_INT_GRID: SO<0  ', 
     .          MU,PHI,XE,YE,ZE,SOX,SOY,SOZ,SO,IPT,ICELL
            STOP
          ENDIF
C              Compute the coordinates of the cell exitting location 
          XE = XE + SO*CX
          YE = YE + SO*CY
          ZE = ZE + SO*CZ

C             Get the intersection face number (i.e. neighptr index)
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
C             Get the next cell to go to
          INEXTCELL = NEIGHPTR(IFACE,ICELL)
          IF (INEXTCELL .LT. 0) THEN
            CALL NEXT_CELL (XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
          ENDIF

C             Get the face grid pointers
          IF (NEIGHPTR(IFACE,ICELL) .GE. 0) THEN
C               If going to same or larger face then use previous face
            KFACE = IFACE
            IC = ICELL
          ELSE
C               If going to smaller face then use next face (more accurate)
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
          ENDIF
          I1 = GRIDPTR(GRIDFACE(1,KFACE),IC)
          I2 = GRIDPTR(GRIDFACE(2,KFACE),IC)
          I3 = GRIDPTR(GRIDFACE(3,KFACE),IC)
          I4 = GRIDPTR(GRIDFACE(4,KFACE),IC)
C             Compute the face interpolation factors
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
          if (u .lt. -0.0001 .or. v. lt. -0.0001 .or. 
     .        u .gt.  1.0001 .or. v .gt.  1.0001) then
            print *, 'BACK_INT_GRID3D: u,v<0 or u,v>1: ',
     .        mu,phi,xe,ye,ze,u,v,iface,jface,kface,
     .        gridpos(1,i1),gridpos(2,i1),gridpos(3,i1)
          endif

C             Get the location coordinate (does the boundary wrapping)
          IF (INEXTCELL .GT. 0) THEN
            IF (JFACE .EQ. 1) THEN
              XE = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
            ELSE IF (JFACE .EQ. 2) THEN
              YE = GRIDPOS(2,GRIDPTR(IOCT,INEXTCELL))
            ELSE
              ZE = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
            ENDIF
          ENDIF

C             Interpolate extinction and source function at face intersection
          F1 = (1-U)*(1-V)
          F2 = (1-U)*V  
          F3 = U*(1-V)
          F4 = U*V
          EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2)
     .         + F3*EXTINCT(I3) + F4*EXTINCT(I4)
C             Correctly interpolate source using extinction*source
          SRCEXT0 = (F1*SOURCE(KANG,I1)*EXTINCT(I1)
     .             + F2*SOURCE(KANG,I2)*EXTINCT(I2)
     .             + F3*SOURCE(KANG,I3)*EXTINCT(I3)
     .             + F4*SOURCE(KANG,I4)*EXTINCT(I4))
C             Compute the cell radiance: integration of the source function
          EXT = 0.5*(EXT0+EXT1)
          TAU=EXT*SO
          IF (TAU .GE. 0.5) THEN
            TRANSCELL = EXP(-TAU)
            ABSCELL = 1.0 - TRANSCELL
          ELSE
            ABSCELL = TAU*(1.0-0.5*TAU*
     .                 (1.0-0.33333333333*TAU*(1-0.25*TAU)))
            TRANSCELL = 1.0 - ABSCELL
          ENDIF
          IF (TAU .LE. 2.0) THEN
            IF (EXT .EQ. 0.0) THEN
              SRC = 0.0
            ELSE
C                 Linear extinction, linear source*extinction, to first order
              SRC = ( 0.5*(SRCEXT0+SRCEXT1) 
     .             + 0.08333333333*(EXT0*SRCEXT1-EXT1*SRCEXT0)*SO )/EXT
            ENDIF
          ELSE
C               Combined first order expansion and constant extinction formula
c            SRC = 0.5/EXT *( SRCEXT0+SRCEXT1
c     .             + (EXT0*SRCEXT1-EXT1*SRCEXT0)*SO
c     .                 *(1-2/TAU+2*TRANSCELL/ABSCELL)/TAU )
            EXT0P = EXT0
            SRCEXT0P = SRCEXT0
            IF (TAU .GT. 4.0) THEN 
              EXT0P = EXT1 + (EXT0-EXT1)*4.0/TAU
              IF (EXT0 .GT. 0.0) SRCEXT0P = SRCEXT0*EXT0P/EXT0
            ENDIF
            SRC = 1.0/(EXT0P+EXT1) *( SRCEXT0P+SRCEXT1
     .             + (EXT0P*SRCEXT1-EXT1*SRCEXT0P)*2.0/(EXT0P+EXT1)
     .                 *(1-2/TAU+2*TRANSCELL/ABSCELL) )
          ENDIF
          SRC = MAX(SRC,0.0D0)

C             Add in the cell radiance and update the transmission.
C             If this is a boundary or the transmission is below the
C             cutoff and we have a valid radiance then set the flag
C             to stop the tracing and add in the interpolated face radiance.
          RAD = RAD + TRANSMIT*SRC*ABSCELL
          TRANSMIT = TRANSMIT*TRANSCELL
          if (RAD < -1.0E-5) then
            print '(A,3(1X,F6.4),8(1X,E12.5))',
     .         'BACK_INT_GRID3D: RAD<0 ',
     .         xe,ye,ze,ext0,ext1,so,tau,transmit,src,abscell,rad
            CALL ABORT_SHDOM_MPI ('The end')
          endif
          VALIDFACE = (GRIDRAD(I1).GE.-0.1 .AND. GRIDRAD(I2).GE.-0.1
     .           .AND. GRIDRAD(I3).GE.-0.1 .AND. GRIDRAD(I4).GE.-0.1) 
          IF (INEXTCELL .LE. 0 .OR.  
     .         (TRANSMIT .LE. TRANSMIN .AND. VALIDFACE)) THEN
            IF (VALIDFACE) THEN
              VALIDRAD = .TRUE.
              RAD0 = F1*GRIDRAD(I1) + F2*GRIDRAD(I2)
     .             + F3*GRIDRAD(I3) + F4*GRIDRAD(I4)
              RAD = RAD + TRANSMIT*RAD0
            ELSE
              print *, 'BACK_INT_GRID3D: INEXTCELL=0: ', xe,ye,ze,icell
              print *, iz,joct,mu,phi
              print '(i7,4(1x,f7.4))', i1,GRIDRAD(I1),
     .                  gridpos(1,i1),gridpos(2,i1),gridpos(3,i1)
              print '(i7,4(1x,f7.4))', i2,GRIDRAD(I2),
     .                  gridpos(1,i2),gridpos(2,i2),gridpos(3,i2)
              print '(i7,4(1x,f7.4))', i3,GRIDRAD(I3),
     .                  gridpos(1,i3),gridpos(2,i3),gridpos(3,i3)
              print '(i7,4(1x,f7.4))', i4,GRIDRAD(I4),
     .                  gridpos(1,i4),gridpos(2,i4),gridpos(3,i4)
              print '(i7,8x,3(1x,f7.4))', ipt,
     .                  gridpos(1,ipt),gridpos(2,ipt),gridpos(3,ipt)
              STOP
            ENDIF
          ELSE
            EXT1 = EXT0
            SRCEXT1 = SRCEXT0
            ICELL = INEXTCELL
          ENDIF
        ENDDO

        GRIDRAD(IPT) = RAD
500     CONTINUE
      ENDDO

      RETURN
      END





      SUBROUTINE BACK_INT_GRID2D (NX,NY,NZ, NA, NPTS, NCELLS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             GRIDPOS, SWEEPORD, MU, PHI, TRANSMIN,
     .             EXTINCT, SOURCE, KANG, GRIDRAD)
C       Sweeps through the spatial grid computing radiances for the
C     discrete ordinate direction (MU,PHI) by integrating the source
C     function and extinction backward from each point.  If a point is
C     already done (nonnegative in GRIDRAD) then it is skipped.
C     For each grid point a ray opposite to the ordinate direction is traced
C     back to a face containing known (valid) radiances.  As each cell is
C     crossed (usually only one) the integration of the source function
C     across the cell is done and added to the radiance.  The radiance
C     ad the end of the ray path is interpolated from the known grid point
C     radiances on the face.  The resulting radiances are put in GRIDRAD.

      INTEGER NX, NY, NZ, NA, NPTS, NCELLS, KANG
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SWEEPORD(NPTS,*)
      INTEGER*2 CELLFLAGS(NCELLS)
      REAL    GRIDPOS(3,NPTS)
      REAL    MU, PHI, TRANSMIN
      REAL    EXTINCT(NPTS), SOURCE(NA,NPTS), GRIDRAD(NPTS)
      INTEGER BITX, BITZ, IOCT, JOCT
      INTEGER IPCELL, ICELL, INEXTCELL, IPT, IFACE
      INTEGER I1, I2, IOPP, JFACE, KFACE, IC
      INTEGER GRIDFACE(2,6), OPPFACE(6)
      INTEGER IORDER, ICORNER, JOCTORDER(8)
      LOGICAL VALIDRAD, VALIDFACE, IPINX, BTEST
      DOUBLE PRECISION PI, CX, CZ, CXINV, CZINV
      DOUBLE PRECISION XE, YE, ZE
      DOUBLE PRECISION SO, SOX, SOZ, EPS
      DOUBLE PRECISION U, F1, F2
      DOUBLE PRECISION EXT, EXT0, EXT1, SRC, SRCEXT0, SRCEXT1
      DOUBLE PRECISION EXT0P, SRCEXT0P
      DOUBLE PRECISION TAU, TRANSCELL, ABSCELL, RAD, RAD0, TRANSMIT
      DATA GRIDFACE/1,5, 2,6,  0,0, 0,0, 1,2, 5,6/
      DATA OPPFACE/2,1,4,3,6,5/
      DATA JOCTORDER/1,3,1,3,2,4,2,4/


      EPS = 1.0E-3*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))

C         Make the ray direction (opposite to the discrete ordinate direction)
      PI = ACOS(-1.0)
      CX = SQRT(1.0-MU**2)*COS(PHI+PI)
      CZ = -MU
      IF (ABS(CX) .GT. 1.0E-5) THEN
        CXINV = 1.0D0/CX
      ELSE
        CX = 0.0
        CXINV = 1.0E6
      ENDIF
      CZINV = 1.0D0/CZ

C         Setup the sweeping direction and the gridpoint location
C         BITc is 0 for positive direction trace back ray, 1 for negative.
      IF (CX .LT. 0.0) THEN
        BITX = 1
      ELSE
        BITX = 0
      ENDIF
      IF (CZ .LT. -1.0E-3) THEN
        BITZ = 1
      ELSE IF (CZ .GT. 1.0E-3) THEN
        BITZ = 0
      ELSE 
        STOP 'BACK_INT_GRID: Bad MU'
      ENDIF
C         IOCT is the octant the ray come from or the discrete ordinate goes to.
      IOCT = 1 + BITX + 4*BITZ
      JOCT = JOCTORDER(IOCT)



C         Sweep through all the grid points
      DO IORDER = 1, NPTS
        IPCELL = ISHFT(SWEEPORD(IORDER,JOCT),-3)
        ICORNER = IBITS(SWEEPORD(IORDER,JOCT),0,3)+1
        IPT = GRIDPTR(ICORNER,IPCELL)
        IF (GRIDRAD(IPT) .GE. 0.0)  GOTO 500

C           Start off at the grid point
        ICELL = IPCELL
        RAD = 0.0
        TRANSMIT = 1.0
        EXT1 = EXTINCT(IPT)
        SRCEXT1 = EXT1*SOURCE(KANG,IPT)
        XE = GRIDPOS(1,IPT)
        YE = GRIDPOS(2,IPT)
        ZE = GRIDPOS(3,IPT)


C           Loop until finding a face with known radiances
        VALIDRAD = .FALSE.
        DO WHILE (.NOT. VALIDRAD)
C             Make sure current cell is valid
          IF (ICELL .LE. 0) THEN
            WRITE (6,*) 'BACK_INT_GRID: ICELL=0  ', 
     .           MU, PHI, IPT, INEXTCELL
            STOP
          ENDIF
          IPINX = BTEST(INT(CELLFLAGS(ICELL)),0)

C             Find boundaries of the current cell
C             Find the three possible intersection planes (X,Y,Z)
C               from the coordinates of the opposite corner grid point
          IOPP = GRIDPTR(9-IOCT,ICELL)
C             Get the distances to the 3 planes and select the closest
          IF (IPINX) THEN
            SOX = 1.0E20
          ELSE
            SOX = (GRIDPOS(1,IOPP)-XE)*CXINV
          ENDIF
          SOZ = (GRIDPOS(3,IOPP)-ZE)*CZINV
          SO = MIN(SOX,SOZ)
          IF (SO .LT. -EPS) THEN
            WRITE (6,*) 'BACK_INT_GRID2D: SO<0  ', 
     .          MU,PHI,XE,ZE,SOX,SOZ,SO,IPT,ICELL
            STOP
          ENDIF
C              Compute the coordinates of the cell exitting location 
          XE = XE + SO*CX
          ZE = ZE + SO*CZ

C             Get the intersection face number (i.e. neighptr index)
          IF (SOX .LE. SOZ) THEN
            IFACE = 2-BITX
            JFACE = 1
          ELSE
            IFACE = 6-BITZ
            JFACE = 3
          ENDIF
C             Get the next cell to go to
          INEXTCELL = NEIGHPTR(IFACE,ICELL)
          IF (INEXTCELL .LT. 0) THEN
            CALL NEXT_CELL (XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
          ENDIF

C             Get the face grid pointers
          IF (NEIGHPTR(IFACE,ICELL) .GE. 0) THEN
C               If going to same or larger face then use previous face
            KFACE = IFACE
            IC = ICELL
          ELSE
C               If going to smaller face then use next face (more accurate)
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
          ENDIF
          I1 = GRIDPTR(GRIDFACE(1,KFACE),IC)
          I2 = GRIDPTR(GRIDFACE(2,KFACE),IC)
C             Compute the face interpolation factors
          IF (JFACE .EQ. 1) THEN
            U = (ZE-GRIDPOS(3,I1))/(GRIDPOS(3,I2)-GRIDPOS(3,I1))
          ELSE
            IF (IPINX) THEN
              U = 0.5
            ELSE
              U = (XE-GRIDPOS(1,I1))/(GRIDPOS(1,I2)-GRIDPOS(1,I1))
            ENDIF
          ENDIF

C             Get the location coordinate (does the boundary wrapping)
          IF (INEXTCELL .GT. 0) THEN
            IF (JFACE .EQ. 1) THEN
              XE = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
            ELSE
              ZE = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
            ENDIF
          ENDIF

C             Interpolate extinction and source function at face intersection
          F1 = 1-U
          F2 = U
          EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2)
C             Correctly interpolate source using extinction*source
          SRCEXT0 = (F1*SOURCE(KANG,I1)*EXTINCT(I1)
     .             + F2*SOURCE(KANG,I2)*EXTINCT(I2))
C             Compute the cell radiance: integration of the source function
          EXT = 0.5*(EXT0+EXT1)
          TAU=EXT*SO
          IF (TAU .GE. 0.5) THEN
            TRANSCELL = EXP(-TAU)
            ABSCELL = 1.0 - TRANSCELL
          ELSE
            ABSCELL = TAU*(1.0-0.5*TAU*
     .                 (1.0-0.33333333333*TAU*(1-0.25*TAU)))
            TRANSCELL = 1.0 - ABSCELL
          ENDIF
          IF (TAU .LE. 2.0) THEN
            IF (EXT .EQ. 0.0) THEN
              SRC = 0.0
            ELSE
C                 Linear extinction, linear source*extinction, to first order
              SRC = ( 0.5*(SRCEXT0+SRCEXT1) 
     .             + 0.08333333333*(EXT0*SRCEXT1-EXT1*SRCEXT0)*SO )/EXT
            ENDIF
          ELSE
C               Combined first order expansion and constant extinction formula
c            SRC = 0.5/EXT *( SRCEXT0+SRCEXT1
c     .             + (EXT0*SRCEXT1-EXT1*SRCEXT0)*SO
c     .                 *(1-2/TAU+2*TRANSCELL/ABSCELL)/TAU )
            EXT0P = EXT0
            SRCEXT0P = SRCEXT0
            IF (TAU .GT. 4.0) THEN 
              EXT0P = EXT1 + (EXT0-EXT1)*4.0/TAU
              IF (EXT0 .GT. 0.0) SRCEXT0P = SRCEXT0*EXT0P/EXT0
            ENDIF
            SRC = 1.0/(EXT0P+EXT1) *( SRCEXT0P+SRCEXT1
     .             + (EXT0P*SRCEXT1-EXT1*SRCEXT0P)*2.0/(EXT0P+EXT1)
     .                 *(1-2/TAU+2*TRANSCELL/ABSCELL) )
          ENDIF
          SRC = MAX(SRC,0.0D0)

C             Add in the cell radiance and update the transmission.
C             If this is a boundary or the transmission is below the
C             cutoff and we have a valid radiance then set the flag
C             to stop the tracing and add in the interpolated face radiance.
          RAD = RAD + TRANSMIT*SRC*ABSCELL
          TRANSMIT = TRANSMIT*TRANSCELL
          VALIDFACE = (GRIDRAD(I1) .GE.-0.1 .AND. GRIDRAD(I2).GE.-0.1)
          IF (INEXTCELL .LE. 0 .OR. 
     .         (TRANSMIT .LE. TRANSMIN .AND. VALIDFACE)) THEN
            IF (VALIDFACE) THEN
              VALIDRAD = .TRUE.
              RAD0 = F1*GRIDRAD(I1) + F2*GRIDRAD(I2)
              RAD = RAD + TRANSMIT*RAD0
            ELSE
              print *, 'BACK_INT_GRID2D: INEXTCELL=0: ', xe,ye,ze,icell
              print *, iz,joct,mu,phi
              print '(i7,4(1x,f7.4))', i1,GRIDRAD(I1),
     .                  gridpos(1,i1),gridpos(2,i1),gridpos(3,i1)
              print '(i7,4(1x,f7.4))', i2,GRIDRAD(I2),
     .                  gridpos(1,i2),gridpos(2,i2),gridpos(3,i2)
              print '(i7,8x,3(1x,f7.4))', ipt,
     .                  gridpos(1,ipt),gridpos(2,ipt),gridpos(3,ipt)
              STOP
            ENDIF
          ELSE
            EXT1 = EXT0
            SRCEXT1 = SRCEXT0
            ICELL = INEXTCELL
          ENDIF
        ENDDO

        GRIDRAD(IPT) = RAD
500     CONTINUE
      ENDDO 
      RETURN
      END






      SUBROUTINE BACK_INT_GRID1D (NX,NY,NZ, NA, NPTS, NCELLS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             GRIDPOS, SWEEPORD, MU, PHI, TRANSMIN,
     .             EXTINCT, SOURCE, KANG, GRIDRAD)
C       Sweeps through the spatial grid computing radiances for the
C     discrete ordinate direction (MU,PHI) by integrating the source
C     function and extinction backward from each point.  If a point is
C     already done (nonnegative in GRIDRAD) then it is skipped.
C     For each grid point a ray opposite to the ordinate direction is traced
C     back to a face containing known (valid) radiances.  As each cell is
C     crossed (usually only one) the integration of the source function
C     across the cell is done and added to the radiance.  The radiance
C     ad the end of the ray path is interpolated from the known grid point
C     radiances on the face.  The resulting radiances are put in GRIDRAD.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NA, NPTS, NCELLS, KANG
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SWEEPORD(NPTS,*)
      INTEGER*2 CELLFLAGS(NCELLS)
      REAL    GRIDPOS(3,NPTS)
      REAL    MU, PHI, TRANSMIN
      REAL    EXTINCT(NPTS), SOURCE(NA,NPTS), GRIDRAD(NPTS)
      INTEGER BITZ, IOCT, JOCT
      INTEGER IPCELL, ICELL, INEXTCELL, IPT, IFACE, I1, IOPP
      INTEGER GRIDFACE(6)
      INTEGER IORDER, ICORNER, JOCTORDER(8)
      LOGICAL VALIDRAD
      DOUBLE PRECISION PI, CZ, CZINV
      DOUBLE PRECISION XE, YE, ZE
      DOUBLE PRECISION SO, EPS
      DOUBLE PRECISION EXT, EXT0, EXT1, SRC, SRCEXT0, SRCEXT1
      DOUBLE PRECISION EXT0P, SRCEXT0P
      DOUBLE PRECISION TAU, TRANSCELL, ABSCELL, RAD, RAD0, TRANSMIT
      DATA GRIDFACE/0, 0, 0, 0, 1, 5/
      DATA JOCTORDER/1,1,1,1,2,2,2,2/


      EPS = 1.0E-3*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))

C         Make the ray direction (opposite to the discrete ordinate direction)
      PI = ACOS(-1.0)
      CZ = -MU
      CZINV = 1.0D0/CZ

C         Setup the sweeping direction and the gridpoint location
C         BITc is 0 for positive direction trace back ray, 1 for negative.
      IF (CZ .LT. -1.0E-3) THEN
        BITZ = 1
      ELSE IF (CZ .GT. 1.0E-3) THEN
        BITZ = 0
      ELSE 
        STOP 'BACK_INT_GRID3D: Bad MU'
      ENDIF
C         IOCT is the octant the ray come from or the discrete ordinate goes to.
      IOCT = 1 + 4*BITZ
      JOCT = JOCTORDER(IOCT)


C         Sweep through all the grid points
      DO IORDER = 1, NPTS
        IPCELL = ISHFT(SWEEPORD(IORDER,JOCT),-3)
        ICORNER = IBITS(SWEEPORD(IORDER,JOCT),0,3)+1
        IPT = GRIDPTR(ICORNER,IPCELL)
        IF (GRIDRAD(IPT) .GE. 0.0)  GOTO 500

C           Start off at the grid point
        ICELL = IPCELL
        RAD = 0.0
        TRANSMIT = 1.0
        EXT1 = EXTINCT(IPT)
        SRCEXT1 = EXT1*SOURCE(KANG,IPT)
        XE = GRIDPOS(1,IPT)
        YE = GRIDPOS(2,IPT)
        ZE = GRIDPOS(3,IPT)

C           Loop until finding a face with known radiances
        VALIDRAD = .FALSE.
        DO WHILE (.NOT. VALIDRAD)
C             Make sure current cell is valid
          IF (ICELL .LE. 0) THEN
            WRITE (6,*) 'BACK_INT_GRID: ICELL=0  ', 
     .           MU, PHI, IPT, INEXTCELL
            STOP
          ENDIF
C             Find boundaries of the current cell
C             Find the three possible intersection planes (X,Y,Z)
C               from the coordinates of the opposite corner grid point
          IOPP = GRIDPTR(9-IOCT,ICELL)
C             Get the distances to the 3 planes and select the closest
          SO = (GRIDPOS(3,IOPP)-ZE)*CZINV
          IF (SO .LT. -EPS) THEN
            WRITE (6,*) 'BACK_INT_GRID1D: SO<0  ', 
     .          MU,PHI,ZE,SO,IPT,ICELL
            STOP
          ENDIF
C              Compute the coordinates of the cell exitting location 
          ZE = ZE + SO*CZ

C             Get the intersection face number (i.e. neighptr index)
          IFACE = 6-BITZ
C             Get the next cell to go to
          INEXTCELL = NEIGHPTR(IFACE,ICELL)

C             Get the face grid pointers
          I1 = GRIDPTR(GRIDFACE(IFACE),ICELL)
C             Get the location coordinate (does the boundary wrapping)
          IF (INEXTCELL .GT. 0) THEN
            ZE = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
          ENDIF

C             Get extinction and source function at face intersection
          EXT0 = EXTINCT(I1)
          SRCEXT0 = EXTINCT(I1)*SOURCE(KANG,I1)
C             Compute the cell radiance: integration of the source function
          EXT = 0.5*(EXT0+EXT1)
          TAU=EXT*SO
          IF (TAU .GE. 0.5) THEN
            TRANSCELL = EXP(-TAU)
            ABSCELL = 1.0 - TRANSCELL
          ELSE
            ABSCELL = TAU*(1.0-0.5*TAU*
     .                 (1.0-0.33333333333*TAU*(1-0.25*TAU)))
            TRANSCELL = 1.0 - ABSCELL
          ENDIF
          IF (TAU .LE. 2.0) THEN
            IF (EXT .EQ. 0.0) THEN
              SRC = 0.0
            ELSE
C                 Linear extinction, linear source*extinction, to first order
              SRC = ( 0.5*(SRCEXT0+SRCEXT1) 
     .             + 0.08333333333*(EXT0*SRCEXT1-EXT1*SRCEXT0)*SO )/EXT
            ENDIF
          ELSE
C               Combined first order expansion and constant extinction formula
c            SRC = 0.5/EXT *( SRCEXT0+SRCEXT1
c     .             + (EXT0*SRCEXT1-EXT1*SRCEXT0)*SO
c     .                 *(1-2/TAU+2*TRANSCELL/ABSCELL)/TAU )
            EXT0P = EXT0
            SRCEXT0P = SRCEXT0
            IF (TAU .GT. 4.0) THEN 
              EXT0P = EXT1 + (EXT0-EXT1)*4.0/TAU
              IF (EXT0 .GT. 0.0) SRCEXT0P = SRCEXT0*EXT0P/EXT0
            ENDIF
            SRC = 1.0/(EXT0P+EXT1) *( SRCEXT0P+SRCEXT1
     .             + (EXT0P*SRCEXT1-EXT1*SRCEXT0P)*2.0/(EXT0P+EXT1)
     .                 *(1-2/TAU+2*TRANSCELL/ABSCELL) )
          ENDIF
          SRC = MAX(SRC,0.0D0)

C             See if there are valid radiances for the intersection face:
C             If so, set flag, interpolate radiance, and add to 
C               radiance from source function in cell
C             If not, only add in cell radiance and prepare for next cell
          IF (GRIDRAD(I1) .GE. -0.1) THEN
            VALIDRAD = .TRUE.
            RAD0 = GRIDRAD(I1)
            RAD = RAD + TRANSMIT*(RAD0*TRANSCELL + SRC*ABSCELL)
          ELSE
            RAD = RAD + TRANSMIT*SRC*ABSCELL
            TRANSMIT = TRANSMIT*TRANSCELL
            EXT1 = EXT0
            SRCEXT1 = SRCEXT0
            ICELL = INEXTCELL
          ENDIF
        ENDDO

        GRIDRAD(IPT) = RAD
500     CONTINUE
      ENDDO
      RETURN
      END


      SUBROUTINE NEXT_CELL (XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,
     .                GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXT)
C       Gets the next cell (INEXT) to go to that is adjacent to this face.

C     IFACE (1-6) has which face (-X,+X,-Y,+Y,-Z,+Z) of the current cell
C     (ICELL) we are at, and JFACE (1-3) is the direction (X,Y,Z).  
C     The current location (XE,YE,ZE) on the face is used to find the cell
C     if there are several adjacent to the current cell.
      IMPLICIT NONE
      INTEGER IFACE, JFACE, ICELL,  INEXT
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
      INTEGER*2 CELLFLAGS(*)
      REAL    GRIDPOS(3,*)
      DOUBLE PRECISION XE, YE, ZE
      INTEGER IC, IC1, DIR

C         Get the pointer to the neighbor cell
      INEXT = NEIGHPTR(IFACE,ICELL)
      IF (INEXT .LT. 0) THEN
C           If neighbor pointer is negative then we are going into a cell
C           that has subcells, and the absolute value of the pointer is the
C           cell to start the tree search.
        IC = ABS(INEXT)
C           Go down the tree until there is no children
        DO WHILE (TREEPTR(2,IC) .GT. 0)
C             Get the direction cell is split
          DIR = IBITS(INT(CELLFLAGS(IC)),2,2)
          IF (DIR .EQ. 0)  STOP 'NEXT_CELL: No split direction'
C             If split is in direction of face then choose the near child cell
C               which is the negative side cell if this is a positive face.
          IC1 = TREEPTR(2,IC)
          IF (DIR .EQ. JFACE) THEN
            IC = IC1 + 1 - MOD(IFACE-1,2)
          ELSE
C             If split is not in direction of face then find which of two cells
C               has a face that contains the current location.
C               Compare the appropriate coordinate (X,Y,Z) of the grid point
C                 on the split line with the current location.
            IC = IC1
            IF (DIR .EQ. 1) THEN
              IF (XE .GT. GRIDPOS(1,GRIDPTR(8,IC1)))  IC = IC + 1
            ELSE IF (DIR .EQ. 2) THEN
              IF (YE .GT. GRIDPOS(2,GRIDPTR(8,IC1)))  IC = IC + 1
            ELSE
              IF (ZE .GT. GRIDPOS(3,GRIDPTR(8,IC1)))  IC = IC + 1
            ENDIF
          ENDIF
        ENDDO
        INEXT = IC
      ENDIF

      RETURN
      END






      LOGICAL FUNCTION SWEEP_BASE_CELL (BCFLAG, NXC,NYC,NZ,IOCT, ICELL, 
     .         IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ)
C       Returns the next base grid cell.  The octant the ray direction
C     is coming from (IOCT) determines the order of sweeping.  ICELL is 0
C     for the first call to set things up.  The function returns false if 
C     we are done sweeping the grid.  Do loops over IZ, IY, and IX are 
C     simulated.  For open boundary condition do the endpoint cells 
C     (in X,Y) first, so the boundary grid points will be done with
C     independent pixel radiative transfer.
      IMPLICIT NONE
      INTEGER BCFLAG, NXC, NYC, NZ, IOCT, ICELL
      LOGICAL BTEST
      INTEGER IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ

      SWEEP_BASE_CELL = .TRUE.
C         If first time here then set up the loops
      IF (ICELL .EQ. 0) THEN
        IF (BTEST(IOCT-1,0)) THEN
          DIX = +1
          IF (BTEST(BCFLAG,0)) THEN
            SIX = NXC
            EIX = MAX(NXC-1,1)
          ELSE
            SIX = 1
            EIX = NXC
          ENDIF
        ELSE
          DIX = -1
          IF (BTEST(BCFLAG,0)) THEN
            SIX = 1
            EIX = MIN(2,NXC)
          ELSE
            SIX = NXC
            EIX = 1
          ENDIF
        ENDIF
        IF (BTEST(IOCT-1,1)) THEN
          DIY = +1
          IF (BTEST(BCFLAG,1)) THEN
            SIY = NYC
            EIY = MAX(NYC-1,1)
          ELSE
            SIY = 1
            EIY = NYC
          ENDIF
        ELSE
          DIY = -1
          IF (BTEST(BCFLAG,1)) THEN
            SIY = 1
            EIY = MIN(2,NYC)
          ELSE
            SIY = NYC
            EIY = 1
          ENDIF
        ENDIF
        IF (BTEST(IOCT-1,2)) THEN
          DIZ = +1
          SIZ = 1
          EIZ = NZ-1
        ELSE
          DIZ = -1
          SIZ = NZ-1
          EIZ = 1
        ENDIF
        IX = SIX
        IY = SIY
        IZ = SIZ
      ELSE
C         Sweep through grid, first in X, in Y, then in Z
        IF (IX .EQ. EIX) THEN
          IF (IY .EQ. EIY) THEN
            IF (IZ .EQ. EIZ) THEN
              SWEEP_BASE_CELL = .FALSE.
              RETURN
            ELSE
              IZ = IZ + DIZ
            ENDIF
            IY = SIY
          ELSE
            IY = MOD(IY+DIY+NYC-1,NYC) + 1
          ENDIF
          IX = SIX
        ELSE
          IX = MOD(IX+DIX+NXC-1,NXC) + 1
        ENDIF
      ENDIF

C         Compute the base grid cell 
      ICELL = IZ + (NZ-1)*(IY-1) + (NZ-1)*NYC*(IX-1)
      RETURN
      END




      LOGICAL FUNCTION SWEEP_NEXT_CELL (BCFLAG, NXC, NYC, NZ, IOCT,
     .                                  ICELL, TREEPTR, CELLFLAGS, 
     .  SP,STACK,IX,IY,IZ,SIX,SIY,SIZ,EIX,EIY,EIZ,DIX,DIY,DIZ)
C       Returns the next grid cell given the last cell (ICELL).  The 
C     octant the ray direction is coming from (IOCT) determines the
C     order of sweeping.  ICELL is 0 for the first call to set things up.
C     The function returns false if we are done sweeping the grid.
C     Within a base cell the tree structure is traced using a stack until 
C     we get to a new end cell.
      IMPLICIT NONE
      INTEGER BCFLAG, NXC, NYC, NZ, IOCT, ICELL, TREEPTR(2,*)
      INTEGER*2 CELLFLAGS(*)
      INTEGER IC, IDIR, SP, STACK(50)
      LOGICAL DONE, SWEEP_BASE_CELL, BTEST
      INTEGER IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ
      
      IC = ICELL

C         Special case for first time here
      IF (ICELL .EQ. 0) THEN      
        IC = 0
        SWEEP_NEXT_CELL = SWEEP_BASE_CELL (BCFLAG,NXC,NYC,NZ,IOCT, IC, 
     .         IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ)
        SP = 0
      ENDIF
C       The scanning direction is the direction (pos or neg) for each 
C         of X,Y,Z that the ray is going in (indicated by IOCT).

C         Go through the tree checking the cells:
C           The stack holds the second child cells, which are the 
C             next cell to go to if at an end cell.
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
C         If this is an end cell then we are either at the desired cell
C           or if this is the old cell then we need to pop the stack
        IF (TREEPTR(2,IC) .EQ. 0) THEN
          IF (IC .NE. ICELL) THEN
            DONE = .TRUE.
          ELSE
C             If no more to pop then we are at the base cell, so go to 
C               next base cell if there is one.
            IF (SP .EQ. 0) THEN
              SWEEP_NEXT_CELL = 
     .            SWEEP_BASE_CELL (BCFLAG, NXC, NYC, NZ, IOCT, IC, 
     .         IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ)
              IF (.NOT. SWEEP_NEXT_CELL)  RETURN
            ELSE
C               Otherwise, pop the stack to get the next cell to go to.
              IC = STACK(SP)
              SP = SP - 1
            ENDIF
          ENDIF
        ELSE
C           Not an end cell, so this cell has children.
C             Depending on the ordinate octant, decide which child to 
C             do first and put the second child on the stack.
          SP = SP + 1
          IF (SP .GT. 50) STOP 'SWEEP_NEXT_CELL: Stack exceeded'
          IDIR = IBITS(INT(CELLFLAGS(IC)),2,2)-1
          IF (BTEST(IOCT-1,IDIR)) THEN
            STACK(SP) = TREEPTR(2,IC)+1
            IC = TREEPTR(2,IC)
          ELSE
            STACK(SP) = TREEPTR(2,IC)
            IC = TREEPTR(2,IC)+1
          ENDIF
        ENDIF
      ENDDO

      ICELL = IC
      SWEEP_NEXT_CELL = .TRUE.
      RETURN
      END






      SUBROUTINE SPLIT_GRID (MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,
     .             DOSPLIT,OUTOFMEM, CURSPLITACC, SPLITCRIT,
     .             NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, 
     .             GRIDPOS, NX, XGRID, NY, YGRID,
     .             ML, MM, NLM, NCS, NLEG, NUMPHASE, 
     .             DELTAM, BCFLAG, IPFLAG, ACCELFLAG,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .             RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, 
     .             TEMP, PLANCK, DIRFLUX, SRCTYPE, 
     .             SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN, 
     .             UNITS, WAVENO, WAVELEN,  ADAPTCRIT, ADAPTIND, 
     .             NPX, NPY, NPZ, DELX, DELY,
     .             XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .             ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .             ZCKD, GASABS, EXTMIN, SCATMIN, CX, CY, CZ, CXINV, 
     .             CYINV, CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .             XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV)
C       Splits the cells that have a cell dividing criterion greater than
C     CURSPLITACC if DOSPLIT is true.  The current splitting criterion
C     achieved is returned in SPLITCRIT.  If we are at the end of the
C     grid point, grid cell, or spherical harmonic arrays (sizes given
C     by MAXIG, MAXIC, MAXIV, MAXIDO) then the OUTOFMEM flag is returned true.
      IMPLICIT NONE
      INTEGER MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX
Cf2py intent(in) :: MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX
      INTEGER NPTS, NCELLS, BCFLAG, IPFLAG
Cf2py intent(in,out) :: NPTS, NCELLS, BCFLAG, IPFLAG
      INTEGER NX, NY, ML, MM, NLM, NCS, NLEG
      INTEGER NUMPHASE
Cf2py intent(in) :: NX, NY, ML, MM, NLM, NCS, NLEG, NUMPHASE
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in,out) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER RSHPTR(*), SHPTR(*), OSHPTR(*), ADAPTIND(*)
Cf2py intent(in,out) :: RSHPTR, SHPTR, OSHPTR
Cf2py intent(in,out) :: ADAPTIND
      INTEGER*2 CELLFLAGS(*)
      INTEGER IPHASE(*)
Cf2py intent(in,out) :: CELLFLAGS, IPHASE
      LOGICAL DELTAM, ACCELFLAG, DOSPLIT, OUTOFMEM
Cf2py intent(in) :: DELTAM, ACCELFLAG, DOSPLIT 
Cf2py intent(in,out) :: OUTOFMEM
      REAL    CURSPLITACC, SPLITCRIT
Cf2py intent(in,out) :: CURSPLITACC, SPLITCRIT
      REAL    SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN(*)
Cf2py intent(in,out) :: SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN
      REAL    WAVENO(2), WAVELEN
Cf2py intent(in) :: WAVENO, WAVELEN
      REAL    GRIDPOS(3,*), XGRID(*), YGRID(*), ADAPTCRIT(*)
Cf2py intent(in,out) :: GRIDPOS, XGRID, YGRID
Cf2py intent(in) :: ADAPTCRIT
      REAL    EXTINCT(*), ALBEDO(*), LEGEN(0:NLEG,*)
Cf2py intent(in,out) :: EXTINCT, ALBEDO, LEGEN
      REAL    RADIANCE(*), SOURCE(*)
Cf2py intent(in,out) :: RADIANCE, SOURCE
      REAL    TEMP(*), PLANCK(*), DIRFLUX(*)
Cf2py intent(in,out) :: TEMP, PLANCK, DIRFLUX
      CHARACTER SRCTYPE*1, UNITS*1
Cf2py intent(in) :: SRCTYPE, UNITS
      INTEGER ICELL, ICELL1, IDIR, I, N, NEWPOINTS(3,4), OLDNPTS
      INTEGER MAXCELLS, MAXPTS, MAXWORK, MAXSH, MAXRAD
      LOGICAL OUTOFMEM0
      REAL    ADAPT(3), FRAC
      
      INTEGER NPX, NPY, NPZ
      REAL DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)
      REAL TEMPP(*), EXTINCTP(*), ALBEDOP(*)
      REAL LEGENP(*), EXTDIRP(*)
      INTEGER IPHASEP(*)
      INTEGER NZCKD
      REAL ZCKD(*), GASABS(*)
      DOUBLE PRECISION EXTMIN, SCATMIN
      
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      INTEGER IPDIRECT, DI, DJ, DK
      DOUBLE PRECISION EPSS, EPSZ, XDOMAIN, YDOMAIN
      DOUBLE PRECISION UNIFORMZLEV, DELXD,DELYD
      
      IF (DOSPLIT) THEN
        OLDNPTS = NPTS
C          Do the cells in batchs, so the newly split cells can be tested
C            for further splitting
        ICELL1 = 1
        DO WHILE (ICELL1 .LE. NCELLS)
C           For the current batch of cells, compute the adaptive grid
C             splitting criterion for 3 directions for each cell.
          N = 0
          DO ICELL = ICELL1, NCELLS
            IF (TREEPTR(2,ICELL) .EQ. 0) THEN
C               If this is a end node cell then compute adaptive criterion
              CALL CELL_SPLIT_TEST (GRIDPTR, GRIDPOS, EXTINCT,
     .                              SHPTR, SOURCE, ICELL, 
     .                              ADAPT, ADAPTCRIT(N+1), IDIR)
              ADAPTIND(N+1) = 4*ICELL + IDIR
              N = N + 1
            ENDIF
          ENDDO
          ICELL1 = NCELLS + 1

C             Sort the adaptive criterion in decreasing order
          CALL SSORT (ADAPTCRIT, ADAPTIND, N, -2)

C             Figure how many cells we can split, leaving 10% room 
C               for cells split during grid smoothing
          FRAC = 0.03
          MAXCELLS = MAXIC - FRAC*(MAXIC-NCELLS) - 2
          MAXPTS = MAXIG - FRAC*(MAXIG-NPTS) - 4
          MAXWORK = MAXIDO - FRAC*(MAXIDO-NPHI0MAX*NPTS) - 4*NPHI0MAX
          MAXSH = MAXIV - FRAC*(MAXIV-SHPTR(NPTS+1)) - 4*NLM
          MAXRAD = MAXIV+MAXIG- FRAC*(MAXIV+MAXIG-RSHPTR(NPTS+1))-4*NLM

C             Go down the list from highest adaptive cell criterion on down
          OUTOFMEM0 = OUTOFMEM
          I = 1
          DO WHILE (ADAPTCRIT(I) .GT. CURSPLITACC .AND. .NOT. OUTOFMEM0
     .              .AND. I .LE. N)
            IF (NCELLS .GT. MAXCELLS .OR. NPTS .GT. MAXPTS .OR.
     .          NPTS*NPHI0MAX .GT. MAXWORK .OR.
     .          SHPTR(NPTS+1) .GT. MAXSH .OR. 
     .          RSHPTR(NPTS+1) .GT. MAXRAD) THEN
              OUTOFMEM0 = .TRUE.
            ELSE
              ICELL = ADAPTIND(I)/4
              IDIR = IAND(ADAPTIND(I),3)
C                 Create the new grid points and make pointers
              CALL DIVIDE_CELL (NPTS, NCELLS, GRIDPTR, NEIGHPTR,
     .               TREEPTR,CELLFLAGS,GRIDPOS,ICELL,IDIR, NEWPOINTS)
C                 Interpolate the medium properties, radiance, and source
              CALL INTERPOLATE_POINT (NEWPOINTS, NPTS,
     .               ML,MM,NCS,NLM, NLEG, NUMPHASE, 
     .               DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS, 
     .               SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN, 
     .               UNITS, WAVENO, WAVELEN, 
     .               EXTINCT,ALBEDO,LEGEN,IPHASE, TEMP,PLANCK, DIRFLUX,
     .               RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, 
     .               NPX, NPY, NPZ, DELX, DELY,
     .               XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .               ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .               ZCKD, GASABS, EXTMIN, SCATMIN, CX, CY, CZ, CXINV, 
     .               CYINV, CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .               XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV)
            ENDIF
            I = I + 1
          ENDDO

          MAXCELLS = MAXIC  - 2
          MAXPTS = MAXIG - 4
          MAXWORK = MAXIDO - 4*NPHI0MAX
          MAXSH = MAXIV - 4*NLM
          MAXRAD = MAXIV+MAXIG - 4*NLM

C             Split some more cells to make the adaptive grid "smoother"
          DO WHILE (.NOT. OUTOFMEM .AND. I .LE. N)
            IF (NCELLS .GT. MAXCELLS .OR. NPTS .GT. MAXPTS .OR.
     .          NPTS*NPHI0MAX .GT. MAXWORK .OR.
     .          SHPTR(NPTS+1) .GT. MAXSH .OR. 
     .          RSHPTR(NPTS+1) .GT. MAXRAD) THEN
              OUTOFMEM = .TRUE.
            ELSE
C                 See if cell should be split to make better grid:
C                   IDIR=0 if no, IDIR=direction (1,2,3) otherwise
              ICELL = ADAPTIND(I)/4
              CALL GRID_SMOOTH_TEST (GRIDPTR, GRIDPOS, NEIGHPTR,
     .                       TREEPTR, CELLFLAGS, SHPTR, ICELL,  IDIR)
C                 Set IDIR=0 to stop adaptive grid smoothing
c              IDIR = 0
              IF (IDIR .GT. 0) THEN
C                   Create the new grid points and make pointers
                CALL DIVIDE_CELL (NPTS, NCELLS, GRIDPTR, NEIGHPTR, 
     .               TREEPTR,CELLFLAGS,GRIDPOS,ICELL,IDIR, NEWPOINTS)
C                   Interpolate the medium properties, radiance, and source
                CALL INTERPOLATE_POINT (NEWPOINTS, NPTS,
     .               ML,MM,NCS, NLM, NLEG, NUMPHASE, 
     .               DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS, 
     .               SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN, 
     .               UNITS, WAVENO, WAVELEN, 
     .               EXTINCT,ALBEDO,LEGEN,IPHASE, TEMP,PLANCK, DIRFLUX,
     .               RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, 
     .               NPX, NPY, NPZ, DELX, DELY,
     .               XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .               ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .               ZCKD, GASABS, EXTMIN, SCATMIN, CX, CY, CZ, CXINV, 
     .               CYINV, CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .               XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV)
              ENDIF
            ENDIF
            I = I + 1
          ENDDO

        ENDDO
        IF (OUTOFMEM0) THEN
          OUTOFMEM = .TRUE.
          WRITE (6,'(2A,I8,A,I8)') 
     .     'Further cell splitting requires larger array allocations:',
     .     '  Ncells=',NCELLS, '  Npnts=',NPTS
        ENDIF

C        If in parallel mode, calculate all the direct beam fluxes at once 
        IF (SRCTYPE .NE. 'T' .AND. 
     .      (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3))) THEN
          CALL MAKE_DIRECT_PAR (OLDNPTS+1, NPTS, BCFLAG, 
     .                            IPFLAG, DELTAM, ML, NLEG,
     .                            SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,
     .                            NX, XGRID, NY, YGRID,  DIRFLUX)
        ENDIF
      ENDIF


C         Find the max adaptive cell criterion after the cell divisions
      SPLITCRIT = 0.0
      DO ICELL = 1, NCELLS
        IF (TREEPTR(2,ICELL) .EQ. 0) THEN
          CALL CELL_SPLIT_TEST (GRIDPTR, GRIDPOS, EXTINCT,
     .                          SHPTR, SOURCE, ICELL, 
     .                          ADAPT, ADAPTCRIT(1), IDIR)
          SPLITCRIT = MAX(SPLITCRIT, ADAPTCRIT(1))
        ENDIF
      ENDDO
      
      RETURN
      END




      SUBROUTINE INTERPOLATE_POINT (NEWPOINTS, NPTS,
     .             ML, MM, NCS, NLM, NLEG, NUMPHASE, 
     .             DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS, 
     .             SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN, 
     .             UNITS, WAVENO, WAVELEN, 
     .             EXTINCT, ALBEDO, LEGEN,IPHASE, TEMP,PLANCK, 
     .             DIRFLUX, RSHPTR, RADIANCE, SHPTR, SOURCE,  
     .             OSHPTR, NPX, NPY, NPZ, DELX, DELY,
     .             XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .             ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .             ZCKD, GASABS, EXTMIN, SCATMIN, CX, CY, CZ, CXINV, 
     .             CYINV, CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .             XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV)
C       Interpolates the medium properties, radiance, and source function for
C     new grid points (specified in NEWPOINTS).  The medium properties are 
C     interpolated from the property grid, since interpolating from the 
C     internal grid may not faithfully represent the true medium.  
C     The delta-M scaling is done to the new grid point properties if needed.
C     The direct beam flux and the Planck function for the new points is 
C     recomputed (rather than interpolated).  The radiance is averaged
C     from the two parent grid points, and used to compute the source
C     function for the new points.
      IMPLICIT NONE
      INTEGER NEWPOINTS(3,4), NPTS
      INTEGER ML, MM, NCS, NLM, NLEG, BCFLAG, IPFLAG
      INTEGER NUMPHASE
      INTEGER RSHPTR(*), SHPTR(*), OSHPTR(*)
      INTEGER IPHASE(NPTS)
      LOGICAL DELTAM, ACCELFLAG
      REAL    SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN(*)
      REAL    WAVENO(2), WAVELEN
      REAL    EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(0:NLEG,NPTS)
      REAL    TEMP(NPTS), PLANCK(NPTS), DIRFLUX(NPTS)
      REAL    RADIANCE(*), SOURCE(*)
      REAL    GRIDPOS(3,*)
      CHARACTER SRCTYPE*1, UNITS*1
      INTEGER I, IP1, IP2, IP, J, K, L, M, ME, MS
      INTEGER IS, NS, NS1, NS2
      INTEGER IR, IR1, IR2, NR, NR1, NR2, SIDE
      LOGICAL VALIDBEAM
      REAL    X, Y, Z, F, BB, C, D, RAD1, RAD2, SECMU0
      REAL    UNIFZLEV, XO, YO, ZO, DIRPATH

      INTEGER, ALLOCATABLE :: LOFJ(:)
      
      INTEGER NPX, NPY, NPZ
      REAL DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)
      REAL TEMPP(*), EXTINCTP(*), ALBEDOP(*)
      REAL LEGENP(*), EXTDIRP(*)
      INTEGER IPHASEP(*)
      INTEGER NZCKD
      REAL ZCKD(*), GASABS(*)
      DOUBLE PRECISION  EXTMIN, SCATMIN
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      INTEGER IPDIRECT, DI, DJ, DK
      DOUBLE PRECISION EPSS, EPSZ, XDOMAIN, YDOMAIN
      DOUBLE PRECISION UNIFORMZLEV, DELXD,DELYD
      
      ALLOCATE (LOFJ(NLM))
C         Make the l index as a function of SH term (J)
      J = 0
      DO L = 0, ML
        ME = MIN(L,MM)
        MS = (1-NCS)*ME
        DO M = MS, ME
          J = J + 1
          LOFJ(J) = L
        ENDDO
      ENDDO
      C = SQRT(4.0*ACOS(-1.0))
      IF (SRCTYPE .NE. 'T') THEN
        SECMU0 = 1.0/ABS(SOLARMU)
      ENDIF


C         Loop over the (potentially) four new grid points
      DO I = 1, 4
        IF (NEWPOINTS(3,I) .GT. 0) THEN
C             Get the parent and new grid pointers
          IP1 = NEWPOINTS(1,I)
          IP2 = NEWPOINTS(2,I)
          IP  = NEWPOINTS(3,I)

C             Interpolate the medium properties from the property grid
          X = GRIDPOS(1,IP)
          Y = GRIDPOS(2,IP)
          Z = GRIDPOS(3,IP)
          CALL TRILIN_INTERP_PROP (X, Y, Z, .FALSE., NLEG, 
     .                             TEMP(IP), EXTINCT(IP), ALBEDO(IP), 
     .                             LEGEN(0,IP), IPHASE(IP), 
     .                      NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .                      XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .                      ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .                      ZCKD, GASABS, EXTMIN, SCATMIN)
C             Do the Delta-M scaling of extinction and albedo for this point
          IF (DELTAM) THEN
            IF (NUMPHASE .GT. 0) THEN
              F = LEGEN(ML+1,IPHASE(IP))
            ELSE
              F = LEGEN(ML+1,IP)
              DO L = 0, ML
                LEGEN(L,IP) = (LEGEN(L,IP) - F)/(1-F)
              ENDDO
            ENDIF
            EXTINCT(IP) = (1.0-ALBEDO(IP)*F)*EXTINCT(IP)
            ALBEDO(IP) = (1.0-F)*ALBEDO(IP)/(1.0-ALBEDO(IP)*F)
          ENDIF
C             Compute the new Planck source function (if needed)
          IF (SRCTYPE .NE. 'S') THEN
            CALL PLANCK_FUNCTION (TEMP(IP), UNITS,WAVENO,WAVELEN,BB)
            PLANCK(IP) = (1.0-ALBEDO(IP))*BB
          ENDIF
C             Compute the new direct beam flux (if needed)
          IF (SRCTYPE .NE. 'T') THEN
            IF (BTEST(BCFLAG,2) .OR. BTEST(BCFLAG,3)) THEN
C               If parallel mode then calculate a temporary approximation 
C                 to the direct flux (for source function calculation below)
              DIRFLUX(IP)=EXP(0.5*(LOG(DIRFLUX(IP1))+LOG(DIRFLUX(IP2))))
            ELSE
C               Otherwise, calculate the exact direct beam from property grid
              DIRPATH = 0.0
              CALL DIRECT_BEAM_PROP (0, X, Y, Z, BCFLAG, IPFLAG, 
     .                  DELTAM,ML,NLEG, SOLARFLUX,SOLARMU,SOLARAZ,
     .                  DIRFLUX(IP),
     .                  UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM, 
     .                  NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .                  XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .                  ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .                  ZCKD, GASABS, CX, CY, CZ, CXINV, CYINV,
     .                  CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .                  XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV)
            ENDIF
          ENDIF

C             Interpolate the radiance
          IR1 = RSHPTR(IP1)
          IR2 = RSHPTR(IP2)
          NR1 = RSHPTR(IP1+1)-RSHPTR(IP1)
          NR2 = RSHPTR(IP2+1)-RSHPTR(IP2)
          NR = MAX(NR1,NR2)
          IR = RSHPTR(IP)
          RSHPTR(IP+1) = IR + NR
          DO J = 1, NR
            IF (J .LE. NR1) THEN
              RAD1 = RADIANCE(IR1+J)
            ELSE
              RAD1 = 0.0
            ENDIF
            IF (J .LE. NR2) THEN
              RAD2 = RADIANCE(IR2+J)
            ELSE
              RAD2 = 0.0
            ENDIF
            RADIANCE(IR+J) = 0.5*(RAD1+RAD2)
          ENDDO

C             Compute the source function for the new points
          NS1 = SHPTR(IP1+1)-SHPTR(IP1)
          NS2 = SHPTR(IP2+1)-SHPTR(IP2)
          NS = MAX(NS1,NS2)
          IS = SHPTR(IP)
          SHPTR(IP+1) = IS + NS
          IF (ACCELFLAG)  OSHPTR(IP+1) = OSHPTR(IP)
          IF (NUMPHASE .GT. 0) THEN
            K = IPHASE(IP)
          ELSE
            K = IP
          ENDIF
          IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
            D = DIRFLUX(IP)*SECMU0
            DO J = 1, NS
              SOURCE(IS+J) = ALBEDO(IP)*LEGEN(LOFJ(J),K)*D*YLMSUN(J)
            ENDDO
            IF (SRCTYPE .EQ. 'B') THEN
              SOURCE(IS+1) = SOURCE(IS+1) + C*PLANCK(IP)
            ENDIF
          ELSE
            SOURCE(IS+1) = C*PLANCK(IP)
            DO J = 2, NS
              SOURCE(IS+J) = 0.0
            ENDDO
          ENDIF
          IR = RSHPTR(IP)
          NR = MIN(NS,RSHPTR(IP+1)-RSHPTR(IP))
          DO J = 1, NR
            SOURCE(IS+J) = SOURCE(IS+J) 
     .                 + ALBEDO(IP)*LEGEN(LOFJ(J),K)*RADIANCE(IR+J)
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END





      SUBROUTINE DIVIDE_CELL (NPTS, NCELLS, GRIDPTR, NEIGHPTR,
     .             TREEPTR, CELLFLAGS, GRIDPOS,  ICELL, IDIR, NEWPOINTS)
C       Divides a single cell (ICELL) in the IDIR direction (1=Z,2=Y,3=X).
C     Makes the pointers to the two new cells (TREEPTR), set the flags
C     (CELLFLAGS), creates the new grid points if they don't exist already,
C     and updates the neighbor pointers (NEIGHPTR) for the new cells and
C     their neighboring cells.
      IMPLICIT NONE
      INTEGER NPTS, NCELLS, ICELL, IDIR, NEWPOINTS(3,4)
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
      INTEGER*2 CELLFLAGS(*)
      REAL    GRIDPOS(3,*)
      INTEGER NEWCELL, IFACE, I
      LOGICAL BTEST

C         Only allow to divide an end cell (i.e. one not yet divided)
      IF (TREEPTR(2,ICELL) .NE. 0) THEN
        STOP 'DIVIDE_CELL: Cannot divide already split cell.'
      ENDIF

C         Do the tree pointers: parent points to negative child, 
C           children point back to parent.
      NEWCELL = NCELLS+1
      NCELLS = NCELLS + 2
      TREEPTR(2,ICELL) = NEWCELL
      TREEPTR(1,NEWCELL) = ICELL
      TREEPTR(2,NEWCELL) = 0
      TREEPTR(1,NEWCELL+1) = ICELL
      TREEPTR(2,NEWCELL+1) = 0
C         Do the flags: put the splitting direction in the parent and 
C           inherit the boundary flags for the children.
      CELLFLAGS(ICELL) = IOR(INT(CELLFLAGS(ICELL)),ISHFT(IDIR,2))
      CELLFLAGS(NEWCELL) = 0
      CELLFLAGS(NEWCELL+1) = 0
      DO I = 0, 1
        IF (BTEST(INT(CELLFLAGS(ICELL)),I)) THEN
          CELLFLAGS(NEWCELL) = IBSET(INT(CELLFLAGS(NEWCELL)),I)
          CELLFLAGS(NEWCELL+1) = IBSET(INT(CELLFLAGS(NEWCELL+1)),I)
        ENDIF
      ENDDO

C         Do the grid points:
      CALL NEW_GRID_POINTS (IDIR, ICELL, NEWCELL, NPTS, 
     .                      GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                      GRIDPOS,  NEWPOINTS)

C         Do the neighbor pointers:
      DO IFACE = 1, 6
C           Special case if the neighbor is the cell itself
        IF (NEIGHPTR(IFACE,ICELL) .EQ. ICELL) THEN 
          NEIGHPTR(IFACE,NEWCELL) = NEWCELL
          NEIGHPTR(IFACE,NEWCELL+1) = NEWCELL+1
C           For the faces in the split direction, the children either
C             point to the other child or inherit the parent's neighbor.
C             Also update the neighbor's neighbor pointer.
        ELSE IF (IFACE .EQ. 2*IDIR) THEN
          NEIGHPTR(IFACE,NEWCELL) = NEWCELL+1
          NEIGHPTR(IFACE,NEWCELL+1) = NEIGHPTR(IFACE,ICELL)
        ELSE IF (IFACE .EQ. 2*IDIR-1) THEN
          NEIGHPTR(IFACE,NEWCELL+1) = NEWCELL
          NEIGHPTR(IFACE,NEWCELL) = NEIGHPTR(IFACE,ICELL)
        ELSE
C           For the other four faces:
C             The children inherit the parent's neighbor pointers
          NEIGHPTR(IFACE,NEWCELL) = NEIGHPTR(IFACE,ICELL)
          NEIGHPTR(IFACE,NEWCELL+1) = NEIGHPTR(IFACE,ICELL)
        ENDIF
C         Then see if we can match faces further down the neighbor's tree;
C           also update the neighboring cells pointers
        CALL MATCH_NEIGHBOR_FACE (IFACE, NEWCELL, NEIGHPTR, 
     .              CELLFLAGS, TREEPTR, GRIDPTR, GRIDPOS)
        CALL MATCH_NEIGHBOR_FACE (IFACE, NEWCELL+1, NEIGHPTR, 
     .              CELLFLAGS, TREEPTR, GRIDPTR, GRIDPOS)
      ENDDO

      RETURN
      END





      SUBROUTINE MATCH_NEIGHBOR_FACE (IFACE, IC, NEIGHPTR, CELLFLAGS,
     .                                TREEPTR, GRIDPTR, GRIDPOS)
C       Resolves the neighbor pointers for a new cell IC (from a subdivision)
C     for the IFACE face.  Based on the size of this new face and
C     the corresponding faces of the neighboring cells (down its tree),
C     the new cell's NEIGHPTR and those of the neighboring cells are updated.
      IMPLICIT NONE
      INTEGER IFACE, IC, NEIGHPTR(6,*), TREEPTR(2,*), GRIDPTR(8,*)
      INTEGER*2 CELLFLAGS(*)
      REAL    GRIDPOS(3,*)
      INTEGER IN, INN, JFACE, IC1, IC8, IN1, IN8, DIR, DIR1, DIR2
      LOGICAL DONE
      REAL    POS(3)
      INTEGER OPPFACE(6)
      DATA    OPPFACE/2,1,4,3,6,5/

      IN = ABS(NEIGHPTR(IFACE,IC))
C         If no neighbor then return
      IF (IN .EQ. 0) RETURN
      JFACE = (IFACE+1)/2
C         Get the center of the new face
      IC1 = GRIDPTR(1,IC)
      IC8 = GRIDPTR(8,IC)
      DIR1 = MOD(JFACE-1+1,3)+1
      DIR2 = MOD(JFACE-1+2,3)+1
      POS(DIR1) = (GRIDPOS(DIR1,IC1)+GRIDPOS(DIR1,IC8))/2
      POS(DIR2) = (GRIDPOS(DIR2,IC1)+GRIDPOS(DIR2,IC8))/2
C         Go down the neighbors cell tree 
      DONE = .FALSE.
      DO WHILE ((.NOT. DONE) .AND. (IC.NE.IN))
C           Compare the coordinates of this neighbor cell to see
C             if the new face is contained within the neighbor face.
        IN1 = GRIDPTR(1,IN)
        IN8 = GRIDPTR(8,IN)
        IF (GRIDPOS(DIR1,IC1) .GE. GRIDPOS(DIR1,IN1) .AND. 
     .      GRIDPOS(DIR1,IC8) .LE. GRIDPOS(DIR1,IN8) .AND. 
     .      GRIDPOS(DIR2,IC1) .GE. GRIDPOS(DIR2,IN1) .AND. 
     .      GRIDPOS(DIR2,IC8) .LE. GRIDPOS(DIR2,IN8)) THEN
C             If this is the bottom of the tree then we have a face match
          IF (TREEPTR(2,IN) .EQ. 0) THEN
            NEIGHPTR(IFACE,IC) = IN
          ELSE 
            NEIGHPTR(IFACE,IC) = -IN
          ENDIF
        ENDIF
C           Compare the size of the neighbor face with the new face.
C             Update the neighbor's neighbor pointer to the new face
C             if the neighbor face is smaller, otherwise make the
C             pointer negative, because it no longer points to an end node.
        IF (GRIDPOS(DIR1,IN1) .GE. GRIDPOS(DIR1,IC1) .AND. 
     .      GRIDPOS(DIR1,IN8) .LE. GRIDPOS(DIR1,IC8) .AND. 
     .      GRIDPOS(DIR2,IN1) .GE. GRIDPOS(DIR2,IC1) .AND. 
     .      GRIDPOS(DIR2,IN8) .LE. GRIDPOS(DIR2,IC8)) THEN
            CALL INHERIT_NEIGHBOR (IN, OPPFACE(IFACE), IC,
     .                             NEIGHPTR, TREEPTR, CELLFLAGS)
        ELSE
          NEIGHPTR(OPPFACE(IFACE),IN)=-ABS(NEIGHPTR(OPPFACE(IFACE),IN))
        ENDIF

C         Go down the neighbors cell tree until there is no children, 
C           finding the cell that contains the center of the new face.
        IF (TREEPTR(2,IN) .EQ. 0) THEN
          DONE = .TRUE.
        ELSE
          DIR = IBITS(INT(CELLFLAGS(IN)),2,2)
          INN = TREEPTR(2,IN)
          IF (DIR .EQ. JFACE) THEN
            IN = INN + 1 - MOD(IFACE-1,2)
          ELSE
            IF (POS(DIR) .GT. GRIDPOS(DIR,GRIDPTR(8,INN))) THEN
              IN = INN + 1
            ELSE
              IN = INN
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END



      SUBROUTINE INHERIT_NEIGHBOR (ICELL, IFACE, IN,
     .                             NEIGHPTR, TREEPTR, CELLFLAGS)
C       Updates the ICELL cell's new neighbor cell (IN) for all 
C     children cells that border the face (IFACE).
      IMPLICIT NONE
      INTEGER ICELL, IFACE, IN, NEIGHPTR(6,*), TREEPTR(2,*)
      INTEGER*2 CELLFLAGS(*)
      INTEGER JFACE, IC, DIR, SP, STACK(50)
      LOGICAL DONE

C         Go through the tree checking the cells:
C           The stack holds the second child cells, which are the 
C             next cell to go to if at an end cell.
      JFACE = (IFACE+1)/2
      IC = ICELL
      SP = 0
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
        NEIGHPTR(IFACE,IC) = IN
C           If this is an end cell then pop the stack
        IF (TREEPTR(2,IC) .EQ. 0) THEN
C             If no more to pop then we are done 
          IF (SP .EQ. 0) THEN
            DONE = .TRUE.
          ELSE
C             Otherwise, pop the stack to get the next cell to go to
            IC = STACK(SP)
            SP = SP - 1
          ENDIF
        ELSE
C             Not and end cell, so this cell has children.
C             If split is in same direction as face then do only close child
          DIR = IBITS(INT(CELLFLAGS(IC)),2,2)
          IF (DIR .EQ. JFACE) THEN
            IC = TREEPTR(2,IC) + MOD(IFACE-1,2)
          ELSE
C             Otherwise need to go to both children so put the second child 
C               on the stack and set the pointer to go to the first child next.
            SP = SP + 1
            IF (SP .GT. 50) STOP 'INHERIT_NEIGHBOR: Stack exceeded'
            STACK(SP) = TREEPTR(2,IC)+1
            IC = TREEPTR(2,IC)
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END







      SUBROUTINE NEW_GRID_POINTS (IDIR, ICELL, NEWCELL, NPTS, 
     .                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                        GRIDPOS,  NEWPOINTS)
C       Sets up the pointers to the grid points for the two new cells.
C     Creates the four new points if they don't already exist, changing
C     NPTS and GRIDPOS accordingly.  Returns the NEWPOINTS(3,4) arrays
C     that gives information about the new points: NEWPOINTS(3,) is the
C     point address of the new points (or zero), and NEWPOINTS(1,) and
C     NEWPOINTS(2,) are the "parent" grid points these should be
C     interpolated from.
      IMPLICIT NONE
      INTEGER IDIR, ICELL, NEWCELL, NPTS,  NEWPOINTS(3,4)
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
      INTEGER*2 CELLFLAGS(*)
      REAL    GRIDPOS(3,*)
      INTEGER I, K, I1, I2, IP1, IP2, IPT, IPMATCH
      INTEGER IFACE1, IFACE2, ICELL2
      INTEGER GRIDCORNER(2,4,3), FACEGRID(2,4,3)
c     integer GRIDFACE(4,6)
      REAL    XP, YP, ZP, EPS
      DATA GRIDCORNER/1,2, 3,4, 5,6, 7,8,  1,3, 2,4, 5,7, 6,8,
     .                1,5, 2,6, 3,7, 4,8/
      DATA FACEGRID/3,5, 4,5, 3,6, 4,6,  1,5, 2,5, 1,6, 2,6,
     .              1,3, 2,3, 1,4, 2,4/

C     In general, when subdividing a cell four new points may be made. 
C     Each of these points belongs to two faces of the original cell.  
C     Use the MATCH_GRID_POINT routine to try to find matches in cells 
C     across the two faces that contain the point; also cross the face of 
C     one of these to get the third bounding cell.  If any of the appropriate
C     gridpoints of these three cells has the same coordinates then the point 
C     already exists, otherwise create the point.

C         Temporarily set new cells' grid pointer to the parent cell's ones
      DO K = 1, 8
        GRIDPTR(K,NEWCELL) = GRIDPTR(K,ICELL)
        GRIDPTR(K,NEWCELL+1) = GRIDPTR(K,ICELL)
      ENDDO

C         Loop over the four edges with potentially new points 
      DO I = 1, 4
C           Get the grid corners and coordinates for this edge and split dir
        I1 = GRIDCORNER(1,I,IDIR)
        I2 = GRIDCORNER(2,I,IDIR)
        IP1 = GRIDPTR(I1,ICELL)
        IP2 = GRIDPTR(I2,ICELL)
        XP = (GRIDPOS(1,IP1) + GRIDPOS(1,IP2))/2
        YP = (GRIDPOS(2,IP1) + GRIDPOS(2,IP2))/2
        ZP = (GRIDPOS(3,IP1) + GRIDPOS(3,IP2))/2

C          Look for gridpoints matching this new point in the three cells
C            that border this edge
        IFACE1 = FACEGRID(1,I,IDIR)
        CALL MATCH_GRID_POINT (XP, YP, ZP, ICELL, IFACE1, GRIDPOS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  IPMATCH)
        IF (IPMATCH .EQ. 0) THEN
          IFACE2 = FACEGRID(2,I,IDIR)
          CALL MATCH_GRID_POINT (XP,YP,ZP, ICELL, IFACE2, GRIDPOS,
     .               GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  IPMATCH)
        ENDIF
        IF (IPMATCH .EQ. 0) THEN
          ICELL2 = ABS(NEIGHPTR(IFACE1,ICELL))
          IF (ICELL2 .GT. 0) THEN
            CALL MATCH_GRID_POINT (XP,YP,ZP, ICELL2, IFACE2, GRIDPOS,
     .               GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  IPMATCH)
          ENDIF
        ENDIF

C           Set some new cells' grid pointers to the new or matching one
        IF (IPMATCH .EQ. 0) THEN
          NPTS = NPTS + 1
          GRIDPTR(I2,NEWCELL) = NPTS
          GRIDPTR(I1,NEWCELL+1) = NPTS
          GRIDPOS(1,NPTS) = XP
          GRIDPOS(2,NPTS) = YP
          GRIDPOS(3,NPTS) = ZP
          NEWPOINTS(1,I) = IP1
          NEWPOINTS(2,I) = IP2
          NEWPOINTS(3,I) = NPTS
        ELSE
          GRIDPTR(I2,NEWCELL) = IPMATCH
          GRIDPTR(I1,NEWCELL+1) = IPMATCH
          NEWPOINTS(3,I) = 0
        ENDIF
      ENDDO

      RETURN
      END



      SUBROUTINE MATCH_GRID_POINT (XP, YP, ZP, ICELL, IFACE,
     .                GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  
     .                IPMATCH)
C       Determines whether the neighboring cells to ICELL adjacent to 
C     face IFACE (1-6: -X,+X,-Y,+Y,-Z,+Z) have gridpoints that match
C     the point at XP,YP,ZP on face IFACE of cell ICELL.
C     If the face of ICELL has several neighboring cells, then the tree 
C     is traced down to find those cells that might contain XP,YP,ZP.
C     IPMATCH is returned 0 if there is no match, otherwise it has the 
C     number of the matching existing gridpoint.
      IMPLICIT NONE
      REAL    XP, YP, ZP
      INTEGER ICELL, IFACE
      REAL    GRIDPOS(3,*)
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
      INTEGER*2 CELLFLAGS(*)
      INTEGER IPMATCH
      INTEGER MAXSTACK
      PARAMETER (MAXSTACK=64)
      INTEGER CELLSTACK(MAXSTACK), SP
      INTEGER INEXT, IDIR, DIR, IC, IC1, KFACE, IPT, I
      INTEGER OPPFACE(6), GRIDFACE(4,6)
      DATA    OPPFACE/2,1,4,3,6,5/, 
     .   GRIDFACE/1,3,5,7, 2,4,6,8, 1,2,5,6, 3,4,7,8, 1,2,3,4, 5,6,7,8/

      IPMATCH = 0
      IDIR = (IFACE+1)/2
      KFACE=OPPFACE(IFACE)

C       Get the pointer to the neighbor cell
      IC = ABS(NEIGHPTR(IFACE,ICELL))
C       Don't process if the neighbor cell is a boundary
      IF (IC .EQ. 0) RETURN

C      Go down all relevant branches of the tree until there are no children
      SP = 0
      DO WHILE (.TRUE.)

C         If at an end cell then check for a gridpoint match at the four 
C           corners of the cell's face
        DO WHILE (TREEPTR(2,IC) .EQ. 0)
          DO I = 1, 4
            IPT = GRIDPTR(GRIDFACE(I,KFACE),IC)
            IF (XP .EQ. GRIDPOS(1,IPT) .AND. YP .EQ. GRIDPOS(2,IPT)
     .          .AND. ZP .EQ. GRIDPOS(3,IPT)) THEN
              IPMATCH = IPT
              RETURN
            ENDIF
          ENDDO
C           Return if there are no more branches on the stack, otherwise pop
          IF (SP .EQ. 0) RETURN
          IC = CELLSTACK(SP)
          SP = SP - 1
        ENDDO

C         Get the direction current cell is split
        DIR = IBITS(INT(CELLFLAGS(IC)),2,2)
        IF (DIR .EQ. 0)  STOP 'MATCH_GRID_POINT: No split direction'
C         If split is in direction of face then choose the near child cell
C           which is the negative side cell if this is a positive face.
        IC1 = TREEPTR(2,IC)
        IF (DIR .EQ. IDIR) THEN
          IC = IC1 + 1 - MOD(IFACE-1,2)
        ELSE
C          If split is not in direction of face then find which of two cells
C          has a face that contains the current location by comparing the 
C          appropriate coordinate (X,Y,Z) of the grid point on the split line 
C          with the current location.  If the split line matches a coordinate
C          then we have to trace down both cells so put one on the stack.
          IC = IC1
          IF (DIR .EQ. 1) THEN
            IF (XP .EQ. GRIDPOS(1,GRIDPTR(8,IC1))) THEN
              SP = SP + 1
              CELLSTACK(SP) = IC+1
            ELSE IF (XP .GT. GRIDPOS(1,GRIDPTR(8,IC1))) THEN
              IC = IC + 1
            ENDIF
          ELSE IF (DIR .EQ. 2) THEN
            IF (YP .EQ. GRIDPOS(2,GRIDPTR(8,IC1))) THEN
              SP = SP + 1
              CELLSTACK(SP) = IC+1
            ELSE IF (YP .GT. GRIDPOS(2,GRIDPTR(8,IC1))) THEN
              IC = IC + 1
            ENDIF
          ELSE
            IF (ZP .EQ. GRIDPOS(3,GRIDPTR(8,IC1))) THEN
              SP = SP + 1
              CELLSTACK(SP) = IC+1
            ELSE IF (ZP .GT. GRIDPOS(3,GRIDPTR(8,IC1))) THEN
              IC = IC + 1
            ENDIF
          ENDIF
          IF (SP .GE. MAXSTACK) 
     .      STOP 'MATCH_GRID_POINT: MAXSTACK exceeded'
        ENDIF
      ENDDO

      END





      SUBROUTINE CELL_SPLIT_TEST (GRIDPTR, GRIDPOS, EXTINCT, 
     .                            SHPTR, SOURCE, ICELL, 
     .                            ADAPTCRIT, MAXADAPT, IDIR)
C       Computes the adaptive cell splitting criterion for the X, Y, and Z
C     directions for cell (ICELL).  The three criterion are returned in 
C     ADAPTCRIT, the max criterion is returned in MAXADAPT, and the 
C     direction of max criterion in IDIR. 
C     The adaptive cell criterion is computed by averaging the division 
C     criterion over the four pairs of grid points that cross the potential 
C     split direction. The cell division criterion is based on the 
C     integration of the source function across the cell (including the 
C     optical depth/transmission factor) that gives the delta radiance.  
C     The criterion is actually the part of the integral that depends on
C     the difference in the extinction times source function across a cell.
C     The rms over all angles of the extinction/source function difference 
C     is computed in spherical harmonic space by summing the square of 
C     the difference over all terms.
      IMPLICIT NONE
      INTEGER ICELL, IDIR
      INTEGER SHPTR(*), GRIDPTR(8,*)
      REAL    GRIDPOS(3,*),  EXTINCT(*), SOURCE(*)
      REAL    ADAPTCRIT(3), MAXADAPT
      INTEGER ID, IE, IP1, IP2, IS1, IS2, NS1, NS2, NS
      INTEGER J, NUM
      INTEGER GRIDCORNER(2,4,3)
      REAL    SPLIT, SUM, JAY, TAU, EXT
      REAL    C0
      DATA GRIDCORNER/1,2, 3,4, 5,6, 7,8,  1,3, 2,4, 5,7, 6,8,
     .                1,5, 2,6, 3,7, 4,8/
      DATA     C0/0.282095/


      MAXADAPT = -1.0
C         Loop over the three possible directions to do the split (X,Y,Z)
      DO ID = 1, 3
        SUM = 0.0
        NUM = 0
C           Loop over four edges that cross the ID direction
        DO IE = 1, 4
C             Get the grid points
          IP1 = GRIDPTR(GRIDCORNER(1,IE,ID),ICELL)
          IP2 = GRIDPTR(GRIDCORNER(2,IE,ID),ICELL)
          IF (IP1 .NE. IP2) THEN
            NUM = NUM + 1
            IS1 = SHPTR(IP1)
            IS2 = SHPTR(IP2)
            NS1 = SHPTR(IP1+1)-IS1
            NS2 = SHPTR(IP2+1)-IS2
            JAY = 0.0
            NS = MIN(NS1,NS2)
C               Compute the source function rms difference integrated 
C               over all angles, which is sum of squares of spherical
C               harmonic term differences.
            DO J = 1, NS
              JAY = JAY + (EXTINCT(IP2)*SOURCE(IS2+J)
     .                     - EXTINCT(IP1)*SOURCE(IS1+J))**2
            ENDDO
            DO J = NS+1, NS1
              JAY = JAY + (EXTINCT(IP1)*SOURCE(IS1+J))**2
            ENDDO
            DO J = NS+1, NS2
              JAY = JAY + (EXTINCT(IP2)*SOURCE(IS2+J))**2
            ENDDO
            EXT = 0.5*(EXTINCT(IP1)+EXTINCT(IP2))
            IF (EXT .GT. 0.0) THEN
              JAY = C0*SQRT(JAY)/EXT
            ELSE
              JAY = 0.0
            ENDIF
            TAU = ABS(EXT*(GRIDPOS(ID,IP2)-GRIDPOS(ID,IP1)))
            SPLIT = ABS(JAY)*(1-EXP(-TAU))
            SUM = SUM + SPLIT
          ENDIF
        ENDDO

        IF (NUM .GT. 0) THEN
          ADAPTCRIT(ID) = SUM/NUM
        ELSE
          ADAPTCRIT(ID) = 0.0
        ENDIF
        IF (ADAPTCRIT(ID) .GT. MAXADAPT) THEN
          MAXADAPT = ADAPTCRIT(ID)
          IDIR = ID
        ENDIF
      ENDDO 

      RETURN
      END




      SUBROUTINE GRID_SMOOTH_TEST (GRIDPTR, GRIDPOS, NEIGHPTR,
     .                        TREEPTR, CELLFLAGS, SHPTR, ICELL,  IDIR)
C       Tests the grid cell ICELL for splitting in order to improve the
C     adaptive grid smoothness.  The direction to split the cell is returned
C     in IDIR (1-X, 2-Y, 3-Z, 0-don't).
C     The criterion for splitting the cell is 
C       1) cell has finer grid cells on either side
C       2) a neighboring cell has grid spacing more than 3 times finer
      IMPLICIT NONE
      INTEGER ICELL, IDIR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*), SHPTR(*)
      INTEGER*2 CELLFLAGS(*)
      REAL    GRIDPOS(3,*)
      INTEGER ID, IFACE, IP, IP1, IP2, IN, I, J, ISUM, IBITS
      INTEGER INCELL(2), DIRSPLIT(2)
      INTEGER FACECORNER(4,6), EDGECORNER(2,3)
      DOUBLE PRECISION XE, YE, ZE
      REAL    RATIO, SIZERATIO, CURGRIDSIZE(3), GRIDSIZE(2), GRIDSIZEINV
      DATA FACECORNER/1,3,5,7, 2,4,6,8, 1,2,5,6, 3,4,7,8,
     .                1,2,3,4, 5,6,7,8/
      DATA EDGECORNER/1,2, 1,3, 1,5/


      IDIR = 0
      SIZERATIO = 1.0
C         Loop over the three possible directions to do the split (X,Y,Z)
      DO ID = 1, 3
C           Get the grid points
        IP1 = GRIDPTR(EDGECORNER(1,ID),ICELL)
        IP2 = GRIDPTR(EDGECORNER(2,ID),ICELL)
        CURGRIDSIZE(ID) = 1.0E20
        IF (IP1 .NE. IP2) THEN
          CURGRIDSIZE(ID) = ABS(GRIDPOS(ID,IP1)-GRIDPOS(ID,IP2)) 
          GRIDSIZEINV = 1.0/CURGRIDSIZE(ID)
C             Loop over two opposite faces
          DO J = 1, 2
            DIRSPLIT(J) = 0
            IFACE = 2*(ID-1) + J
C               Get location of center of cell's face
            XE = 0.0
            YE = 0.0
            ZE = 0.0
            DO I = 1, 4
              IP = GRIDPTR(FACECORNER(I,IFACE),ICELL)
              XE = XE + GRIDPOS(1,IP)*0.25
              YE = YE + GRIDPOS(2,IP)*0.25
              ZE = ZE + GRIDPOS(3,IP)*0.25
            ENDDO
C               Get the neighboring cell and its size in this direction
            CALL NEXT_CELL (XE, YE, ZE, IFACE, ID, ICELL, GRIDPOS,
     .                 GRIDPTR,NEIGHPTR,TREEPTR,CELLFLAGS,  INCELL(J))
C               If the neighbor does not exist or is a base cell or
C                 is an independent pixel cell in X or Y then set
C                 its grid size the same as current cell (so don't split).
            IF (INCELL(J) .EQ. 0) THEN
              GRIDSIZE(J) = CURGRIDSIZE(ID) 
            ELSE
              GRIDSIZE(J) = ABS( GRIDPOS(ID,GRIDPTR(1,INCELL(J)))
     .                         - GRIDPOS(ID,GRIDPTR(8,INCELL(J))) )

              IF (TREEPTR(1,INCELL(J)) .EQ. 0) THEN
                GRIDSIZE(J) = CURGRIDSIZE(ID) 
              ENDIF
              IF (ID .LE. 2 .AND. BTEST(CELLFLAGS(INCELL(J)),ID-1)) THEN
                GRIDSIZE(J) = CURGRIDSIZE(ID) 
              ENDIF
            ENDIF
C               If there are multiple neighboring cells then get 
C                 direction the neighbor cell was split
            IN = NEIGHPTR(IFACE,ICELL)
            IF (IN .LT. 0) THEN
              DIRSPLIT(J) = IBITS(INT(CELLFLAGS(ABS(IN))),2,2)
            ENDIF
          ENDDO

C             Split cell this way if cells on opposite side are smaller in
C               this direction OR opposite neighbor cells are split in
C               same other direction.
          IF (GRIDSIZE(1)*GRIDSIZEINV .LT. 0.75 .AND. 
     .        GRIDSIZE(2)*GRIDSIZEINV .LT. 0.75) THEN
            IDIR = ID
          ENDIF
          IF (DIRSPLIT(1) .GT. 0 .AND. DIRSPLIT(1) .NE. ID .AND.
     .        DIRSPLIT(1) .EQ. DIRSPLIT(2)) THEN
            IDIR = DIRSPLIT(1)
          ENDIF

C             Split cell this way if ratio of neighboring grid cells' size
C               to this cells size (in ID direction) is less than 0.4.
          RATIO = MIN(GRIDSIZE(1)*GRIDSIZEINV, GRIDSIZE(2)*GRIDSIZEINV)
          IF (RATIO .LT. 0.4 .AND. RATIO .LT. SIZERATIO) THEN
C               Check if there is any source function for this cell, before
C                 deciding to split it (no point in dividing clear sky).
            ISUM = 0
            DO I = 1, 8
              IP = GRIDPTR(I,ICELL)
              ISUM = ISUM + SHPTR(IP+1)-SHPTR(IP)
            ENDDO
            IF (ISUM .GT. 0) THEN
              IDIR = ID
              SIZERATIO = RATIO
            ENDIF
          ENDIF

        ENDIF        
      ENDDO

      RETURN
      END


