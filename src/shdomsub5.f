C     Joint routines for both Polarized and Unpolarized sources
C     This contains adaptations of subroutines written by Frank Evans
C     for use in SHDOM.
C     https://nit.coloradolinux.com/~evans/shdom.html
C     See shdom.txt for documentation.

C     These subroutines have been written for use in pyshdom by
C     Aviad Levis, Technion Institute of Technology, 2019 and
C     later modified by Jesse Loveridge, University of Illinois at Urbana-Champaign,
C     2020-2021.

C     Most of these subroutines use the same core ray tracing.
C     They just differ in their output shape (and values.)
C     To ease maintainence of this code it would be better to merge
C     all methods.

      SUBROUTINE MIN_OPTICAL_DEPTH(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS, EXTINCT,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX, PATH, IPHASE,
     .             ALBEDO, LEGEN, DELTAMPATH, DELTAM,
     .             NPART, NSTLEG, NLEG,ML,PATHS_SIZE,
     .             INTERPMETHOD, PHASEMAX, PHASEINTERPWT,
     .             MAXNMICRO
     .             )
C     Integrates the extinction field (EXTINCT) along a ray and stores the minimum
C     optical path that intersects a cell adjacent to each gridpoint.
C     Stores only the minimum across all the set of rays.
C     -JRLoveridge 2021/02/22.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS, PATHS_SIZE
Cf2py intent(in) :: NX, NY, NZ, NPTS, NCELLS, PATHS_SIZE
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) :: BCFLAG, IPFLAG
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in) :: CELLFLAGS
      DOUBLE PRECISION    CAMX(*), CAMY(*), CAMZ(*), CAMMU(*)
      DOUBLE PRECISION    CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, NLEG, NPART, NSTLEG, MAXNMICRO
Cf2py intent(in) :: NPIX,NLEG,NPART,NSTLEG, MAXNMICRO
      REAL  PATH(NPTS,PATHS_SIZE)
Cf2py intent(out):: PATH
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
Cf2py intent(in) :: EXTINCT, ALBEDO
      REAL  LEGEN(NSTLEG,0:NLEG,*)
Cf2py intent(in) :: LEGEN
      INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART), ML
      REAL PHASEINTERPWT(8*MAXNMICRO,NPTS,NPART), PHASEMAX
Cf2py intent(in) :: IPHASE, ML, PHASEINTERPWT,PHASEMAX
      CHARACTER INTERPMETHOD*2
Cf2py intent(in) :: INTERPMETHOD
      LOGICAL DELTAMPATH, DELTAM
Cf2py intent(in) :: DELTAMPATH, DELTAM
      INTEGER N, K
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, R, PI

      PI = ACOS(-1.0D0)
      PATH = 999.0
C         Loop over pixels in image
      DO N = 1, NPIX
        X0 = CAMX(N)
        Y0 = CAMY(N)
        Z0 = CAMZ(N)
        MU2 = CAMMU(N)
        PHI2 = CAMPHI(N)
        MURAY = -MU2
        PHIRAY = PHI2 - PI

C             Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (MURAY .GE. 0.0) THEN
            GOTO 900
          ENDIF
          R = (ZGRID(NZ) - Z0)/MURAY
          X0 = X0 + R*SQRT(1-MURAY**2)*COS(PHIRAY)
          Y0 = Y0 + R*SQRT(1-MURAY**2)*SIN(PHIRAY)
          Z0 = ZGRID(NZ)
        ELSE IF (Z0 .LT. ZGRID(1)) THEN
          WRITE (6,*) 'SPACE_CARVING: Level', Z0, 'below domain',
     .		                      ZGRID(1)
          STOP
        ENDIF

        CALL MIN_OPTICAL_DEPTH_1RAY(NX, NY, NZ, NPTS, NCELLS, EXTINCT,
     .                          GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                          BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                          GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                X0, Y0, Z0, PATH, IPHASE,
     .                      ALBEDO, LEGEN, DELTAMPATH, DELTAM,
     .                      NPART, NSTLEG, NLEG,ML,N, PATHS_SIZE,
     .                      INTERPMETHOD, PHASEMAX, PHASEINTERPWT,
     .                      MAXNMICRO)
900     CONTINUE

      ENDDO

      RETURN
      END

      SUBROUTINE MIN_OPTICAL_DEPTH_1RAY(NX, NY, NZ, NPTS,
     .                      NCELLS, EXTINCT,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			             X0, Y0, Z0, PATH, IPHASE,
     .                      ALBEDO, LEGEN, DELTAMPATH, DELTAM,
     .                      NPART, NSTLEG, NLEG,ML,N, PATHS_SIZE,
     .                      INTERPMETHOD, PHASEMAX, PHASEINTERPWT,
     .                      MAXNMICRO)

C     Integrates the extinction field (EXTINCT) along a ray and stores the minimum
C     optical path that intersects a cell adjacent to each gridpoint.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS, PATHS_SIZE
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(*), N, MAXNMICRO
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      DOUBLE PRECISION    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0
      REAL    TAU, EXT1, EXTN
      REAL    PATH(NPTS, PATHS_SIZE)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, K
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), BCFLAG, IPFLAG
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    XM, YM
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, F(8)

      INTEGER NLEG,NPART, NSTLEG
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL  LEGEN(NSTLEG,0:NLEG,*), PHASEMAX
      REAL  PHASEINTERPWT(8*MAXNMICRO,NPTS,NPART)
      CHARACTER INTERPMETHOD*2
      INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART), ML
      LOGICAL DELTAMPATH, DELTAM

      REAL FM, UNSCALED_ALBEDO
      INTEGER J, IPART, Q

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/

      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)
      TAU = 0.0

C         Make the ray direction (opposite to the discrete ordinate direction)
      CX = SQRT(1.0-MURAY**2)*COS(PHIRAY)
      CY = SQRT(1.0-MURAY**2)*SIN(PHIRAY)
      CZ = MURAY
      IF (ABS(CX) .GT. 1.0E-6) THEN
        CXINV = 1.0D0/CX
      ELSE
        CX = 0.0
        CXINV = 1.0E6
      ENDIF
      IF (ABS(CY) .GT. 1.0E-6) THEN
        CYINV = 1.0D0/CY
      ELSE
        CY = 0.0
        CYINV = 1.0E6
      ENDIF
      IF (ABS(CZ) .GT. 1.0E-6) THEN
        CZINV = 1.0D0/CZ
      ELSE
        CZ = 0.0
        CZINV = 1.0E6
      ENDIF

C         Setup for the ray path direction
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
      IF (CZ .LT. 0.0) THEN
        BITZ = 1
      ELSE
        BITZ = 0
      ENDIF
      IOCT = 1 + BITX + 2*BITY + 4*BITZ
      XM = 0.5*(XGRID(1)+XGRID(NX))
      YM = 0.5*(YGRID(1)+YGRID(NY))

C         Start at the desired point, getting the extinction and source there
      XE = X0
      YE = Y0
      ZE = Z0
      CALL LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, XE, YE, ZE, ICELL)

      IFACE = 0
      NGRID = 0

C         Loop until reach a Z boundary or transmission is very small
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
C           Make sure current cell is valid
        IF (ICELL .LE. 0) THEN
          WRITE (6,*)'INTEGRATE_1RAY: ICELL=',ICELL,
     .                MURAY,PHIRAY,XE,YE,ZE
          STOP
        ENDIF
        NGRID = NGRID + 1

C         Interpolate the source and extinction to the current point
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XE, YE, ZE, F)

        EXT1 = 0.0
        DO IPART=1,NPART
          DO J=1,8
            IF (DELTAM .AND. .NOT. DELTAMPATH)THEN
              IF (PHASEINTERPWT(1,GRIDPTR(J,ICELL),IPART)
     .            .GE. PHASEMAX) THEN
                FM = LEGEN(1,ML+1,IPHASE(1,GRIDPTR(J,ICELL),IPART))
              ELSE
                FM = 0.0
                DO Q=1,8*MAXNMICRO
                  FM = FM +
     .              LEGEN(1,ML+1,IPHASE(Q,GRIDPTR(J,ICELL),IPART))
                ENDDO
              ENDIF
              UNSCALED_ALBEDO = ALBEDO(GRIDPTR(J,ICELL),IPART)/
     .                  (FM*(ALBEDO(GRIDPTR(J,ICELL),IPART)-1.0)+ 1.0)
              EXT1 = EXT1 + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)/
     .                   (1.0 - UNSCALED_ALBEDO*FM)
            ELSE
              EXT1 = EXT1 + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ENDIF
          ENDDO
       ENDDO

C           This cell is independent pixel if IP mode or open boundary
C             conditions and ray is leaving domain (i.e. not entering)
        IPINX = BTEST(INT(CELLFLAGS(ICELL)),0) .AND.
     .          .NOT. ( BTEST(BCFLAG,0) .AND.
     .          ((CX.GT.0.AND.XE.LT.XM) .OR. (CX.LT.0.AND.XE.GT.XM)) )
        IPINY = BTEST(INT(CELLFLAGS(ICELL)),1) .AND.
     .          .NOT. ( BTEST(BCFLAG,1) .AND.
     .          ((CY.GT.0.AND.YE.LT.YM) .OR. (CY.LT.0.AND.YE.GT.YM)) )

C           Find boundaries of the current cell
C           Find the three possible intersection planes (X,Y,Z)
C             from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9-IOCT,ICELL)
C           Get the distances to the 3 planes and select the closest
C             (always need to deal with the cell that is wrapped)
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
          WRITE (6,*) 'INTEGRATE_1RAY: SO<0  ',
     .      MURAY,PHIRAY,XE,YE,ZE,SO,ICELL
          STOP
        ENDIF
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ
C           Find the optical path across the grid cell and figure how
C             many subgrid intervals to use
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XN, YN, ZN, F)

        EXTN = 0.0
        DO IPART=1,NPART
          DO J=1,8
            IF (DELTAM .AND. .NOT. DELTAMPATH)THEN
              IF (PHASEINTERPWT(1,GRIDPTR(J,ICELL),IPART)
     .            .GE. PHASEMAX) THEN
                FM = LEGEN(1,ML+1,IPHASE(1,GRIDPTR(J,ICELL),IPART))
              ELSE
                FM = 0.0
                DO Q=1,8*MAXNMICRO
                  FM = FM +
     .              LEGEN(1,ML+1,IPHASE(Q,GRIDPTR(J,ICELL),IPART))
                ENDDO
              ENDIF
              UNSCALED_ALBEDO = ALBEDO(GRIDPTR(J,ICELL),IPART)/
     .                  (FM*(ALBEDO(GRIDPTR(J,ICELL),IPART)-1.0)+ 1.0)
              EXTN = EXTN + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)/
     .                   (1.0 - UNSCALED_ALBEDO*FM)
            ELSE
              EXTN = EXTN + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ENDIF
          ENDDO
       ENDDO

        TAU = TAU + SO*0.5*(EXT1+EXTN)
C        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
C     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        OUTOFDOMAIN = .FALSE.
        IF (.NOT.OUTOFDOMAIN) THEN
          IF (PATHS_SIZE .EQ. 1) THEN
            DO K=1,8
              PATH(GRIDPTR(K,ICELL),1) =
     .          MIN(TAU,PATH(GRIDPTR(K,ICELL),1))
            ENDDO
          ELSE
            DO K=1,8
              PATH(GRIDPTR(K,ICELL),N) =
     .            MIN(TAU,PATH(GRIDPTR(K,ICELL),N))
            ENDDO
          ENDIF
        ENDIF

C               Get the intersection face number (i.e. neighptr index)
        IF (SOX .LE. SOZ .AND. SOX .LE. SOY) THEN
          IFACE = 2-BITX
          JFACE = 1
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),0).AND.BTEST(BCFLAG,0)
        ELSE IF (SOY .LE. SOZ) THEN
          IFACE = 4-BITY
          JFACE = 2
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),1).AND.BTEST(BCFLAG,1)
        ELSE
          IFACE = 6-BITZ
          JFACE = 3
          OPENBCFACE=.FALSE.
        ENDIF
C            Get the next cell to go to
        INEXTCELL = NEIGHPTR(IFACE,ICELL)
        IF (INEXTCELL .LT. 0) THEN
          CALL NEXT_CELL (XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
        ENDIF
C             If going to same or larger face then use previous face
        IF (NEIGHPTR(IFACE,ICELL) .GE. 0 .AND. .NOT.OPENBCFACE) THEN
          KFACE = IFACE
          IC = ICELL
        ELSE
C             If going to smaller face then use next face (more accurate)
          KFACE = OPPFACE(IFACE)
          IC = INEXTCELL
          IFACE = 0
        ENDIF
C           Get the location coordinate
        IF (INEXTCELL .GT. 0) THEN
          IF (JFACE .EQ. 1) THEN
            XN = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
          ELSE IF (JFACE .EQ. 2) THEN
            YN = GRIDPOS(2,GRIDPTR(IOCT,INEXTCELL))
          ELSE
            ZN = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
          ENDIF
        ENDIF

C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
        IF ((INEXTCELL .EQ. 0).OR.(NGRID.GT.MAXCELLSCROSS)) THEN
          DONE = .TRUE.
        ELSE
          XE = XN
          YE = YN
          ZE = ZN
          ICELL = INEXTCELL
        ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE OPTICAL_DEPTH(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,EXTINCT,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX, PATH, IPHASE,
     .             ALBEDO, LEGEN, DELTAMPATH, DELTAM,
     .             NPART, NSTLEG, NLEG, ML,
     .             INTERPMETHOD, PHASEMAX, PHASEINTERPWT, MAXNMICRO)
C    Integrates the extinction field along the line of sight,
C    this may be delta-M scaled extinction or not depending on
C    DELTAMPATH (and DELTAM). - JRLoveridge 2021/02/22
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
Cf2py intent(in) :: NX, NY, NZ, NPTS, NCELLS
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) :: BCFLAG, IPFLAG
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in) :: CELLFLAGS
      DOUBLE PRECISION CAMX(*), CAMY(*), CAMZ(*), CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, NLEG, NSTLEG, NPART, MAXNMICRO
Cf2py intent(in) :: NPIX, NLEG, NSTLEG, NPART, MAXNMICRO
      REAL  PATH(NPIX)
Cf2py intent(out):: PATH
      INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART), ML
      REAL PHASEINTERPWT(8*MAXNMICRO,NPTS,NPART), PHASEMAX
      CHARACTER INTERPMETHOD*2
Cf2py intent(in) :: IPHASE, ML
      REAL EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
Cf2py intent(in) :: EXTINCT, ALBEDO
      REAL LEGEN(NSTLEG,0:NLEG,*)
Cf2py intent(in) :: LEGEN
      LOGICAL DELTAMPATH, DELTAM
Cf2py intent(in) :: DELTAMPATH, DELTAM

      REAL TAU
      INTEGER N, K
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, R, PI

      PI = ACOS(-1.0D0)
      PATH = 0.0
C         Loop over pixels in image
      DO N = 1, NPIX
        X0 = CAMX(N)
        Y0 = CAMY(N)
        Z0 = CAMZ(N)
        MU2 = CAMMU(N)
        PHI2 = CAMPHI(N)
        MURAY = -MU2
        PHIRAY = PHI2 - PI

C             Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (MURAY .GE. 0.0) THEN
            GOTO 900
          ENDIF
          R = (ZGRID(NZ) - Z0)/MURAY
          X0 = X0 + R*SQRT(1-MURAY**2)*COS(PHIRAY)
          Y0 = Y0 + R*SQRT(1-MURAY**2)*SIN(PHIRAY)
          Z0 = ZGRID(NZ)
        ELSE IF (Z0 .LT. ZGRID(1)) THEN
          WRITE (6,*) 'OPTICAL_DEPTH: Level', Z0, 'below domain',
     .		                      ZGRID(1)
          STOP
        ENDIF
        TAU = 0.0
        CALL OPTICAL_DEPTH_1RAY(NX, NY, NZ, NPTS, NCELLS, EXTINCT,
     .                          GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                          BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                          GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                X0, Y0, Z0, TAU, ALBEDO,LEGEN,IPHASE,
     .                      NSTLEG, NLEG,NPART, DELTAMPATH,DELTAM,
     .                      ML,INTERPMETHOD,PHASEINTERPWT,PHASEMAX,
     .                      MAXNMICRO)
        PATH(N) = TAU
900     CONTINUE

      ENDDO

      RETURN
      END

      SUBROUTINE OPTICAL_DEPTH_1RAY(NX, NY, NZ, NPTS, NCELLS, EXTINCT,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			             X0, Y0, Z0, TAU, ALBEDO, LEGEN, IPHASE,
     .                   NSTLEG, NLEG, NPART, DELTAMPATH, DELTAM,
     .                   ML,INTERPMETHOD,PHASEINTERPWT,PHASEMAX,
     .                   MAXNMICRO)

C       Integrates the source function through the extinction field
C     (EXTINCT) backward in the direction (MURAY,PHIRAY) to find the
C     outgoing radiance (RAD) at the point X0,Y0,Z0.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0
      INTEGER NPART, NSTLEG, NLEG, MAXNMICRO
      INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART)
      REAL    PHASEINTERPWT(8*MAXNMICRO,NPTS,NPART), PHASEMAX
      CHARACTER INTERPMETHOD*2
      REAL EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL LEGEN(NSTLEG,0:NLEG,*)
      LOGICAL DELTAMPATH, DELTAM
      REAL    TAU, EXT1, EXTN
C      REAL    PATH(*)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, K
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), BCFLAG, IPFLAG
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    XM, YM
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, F(8)

      REAL FM, UNSCALED_ALBEDO
      INTEGER J, IPART, ML, Q

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/


      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)
C      TAU = 0.0

C         Make the ray direction (opposite to the discrete ordinate direction)
      CX = SQRT(1.0-MURAY**2)*COS(PHIRAY)
      CY = SQRT(1.0-MURAY**2)*SIN(PHIRAY)
      CZ = MURAY
      IF (ABS(CX) .GT. 1.0E-6) THEN
        CXINV = 1.0D0/CX
      ELSE
        CX = 0.0
        CXINV = 1.0E6
      ENDIF
      IF (ABS(CY) .GT. 1.0E-6) THEN
        CYINV = 1.0D0/CY
      ELSE
        CY = 0.0
        CYINV = 1.0E6
      ENDIF
      IF (ABS(CZ) .GT. 1.0E-6) THEN
        CZINV = 1.0D0/CZ
      ELSE
        CZ = 0.0
        CZINV = 1.0E6
      ENDIF

C         Setup for the ray path direction
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
      IF (CZ .LT. 0.0) THEN
        BITZ = 1
      ELSE
        BITZ = 0
      ENDIF
      IOCT = 1 + BITX + 2*BITY + 4*BITZ
      XM = 0.5*(XGRID(1)+XGRID(NX))
      YM = 0.5*(YGRID(1)+YGRID(NY))

C         Start at the desired point, getting the extinction and source there
      XE = X0
      YE = Y0
      ZE = Z0
      CALL LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, XE, YE, ZE, ICELL)

      IFACE = 0
      NGRID = 0

C         Loop until reach a Z boundary or transmission is very small
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
C           Make sure current cell is valid
        IF (ICELL .LE. 0) THEN
          WRITE (6,*)'INTEGRATE_1RAY: ICELL=',ICELL,
     .                MURAY,PHIRAY,XE,YE,ZE
          STOP
        ENDIF
        NGRID = NGRID + 1

C         Interpolate the source and extinction to the current point
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XE, YE, ZE, F)
        EXT1 = 0.0
        DO IPART=1,NPART
          DO J=1,8
            IF (DELTAM .AND. .NOT. DELTAMPATH)THEN
              IF (PHASEINTERPWT(1,GRIDPTR(J,ICELL),IPART)
     .            .GE. PHASEMAX) THEN
                FM = LEGEN(1,ML+1,IPHASE(1,GRIDPTR(J,ICELL),IPART))
              ELSE
                FM = 0.0
                DO Q=1,8*MAXNMICRO
                  FM = FM +
     .              PHASEINTERPWT(Q,GRIDPTR(J,ICELL),IPART)*
     .              LEGEN(1,ML+1,IPHASE(Q,GRIDPTR(J,ICELL),IPART))
                ENDDO
              ENDIF
              UNSCALED_ALBEDO = ALBEDO(GRIDPTR(J,ICELL),IPART)/
     .                  (FM*(ALBEDO(GRIDPTR(J,ICELL),IPART)-1.0)+ 1.0)
              EXT1 = EXT1 + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)/
     .                   (1.0 - UNSCALED_ALBEDO*FM)
            ELSE
              EXT1 = EXT1 + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ENDIF
          ENDDO
       ENDDO
C           This cell is independent pixel if IP mode or open boundary
C             conditions and ray is leaving domain (i.e. not entering)
        IPINX = BTEST(INT(CELLFLAGS(ICELL)),0) .AND.
     .          .NOT. ( BTEST(BCFLAG,0) .AND.
     .          ((CX.GT.0.AND.XE.LT.XM) .OR. (CX.LT.0.AND.XE.GT.XM)) )
        IPINY = BTEST(INT(CELLFLAGS(ICELL)),1) .AND.
     .          .NOT. ( BTEST(BCFLAG,1) .AND.
     .          ((CY.GT.0.AND.YE.LT.YM) .OR. (CY.LT.0.AND.YE.GT.YM)) )

C           Find boundaries of the current cell
C           Find the three possible intersection planes (X,Y,Z)
C             from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9-IOCT,ICELL)
C           Get the distances to the 3 planes and select the closest
C             (always need to deal with the cell that is wrapped)
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
          WRITE (6,*) 'INTEGRATE_1RAY: SO<0  ',
     .      MURAY,PHIRAY,XE,YE,ZE,SO,ICELL
          STOP
        ENDIF
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ
C           Find the optical path across the grid cell and figure how
C             many subgrid intervals to use
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XN, YN, ZN, F)
        EXTN = 0.0
        DO IPART=1,NPART
          DO J=1,8
            IF (DELTAM .AND. .NOT. DELTAMPATH)THEN
              IF (PHASEINTERPWT(1,GRIDPTR(J,ICELL),IPART)
     .            .GE. PHASEMAX) THEN
                FM = LEGEN(1,ML+1,IPHASE(1,GRIDPTR(J,ICELL),IPART))
              ELSE
                FM = 0.0
                DO Q=1,8*MAXNMICRO
                  FM = FM +
     .              PHASEINTERPWT(Q,GRIDPTR(J,ICELL),IPART)*
     .              LEGEN(1,ML+1,IPHASE(Q,GRIDPTR(J,ICELL),IPART))
                ENDDO
              ENDIF
              UNSCALED_ALBEDO = ALBEDO(GRIDPTR(J,ICELL),IPART)/
     .                  (FM*(ALBEDO(GRIDPTR(J,ICELL),IPART)-1.0)+ 1.0)
              EXTN = EXTN + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)/
     .                   (1.0 - UNSCALED_ALBEDO*FM)
            ELSE
              EXTN = EXTN + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ENDIF
          ENDDO
       ENDDO
        TAU = TAU + SO*0.5*(EXT1+EXTN)
        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
C        IF (.NOT.OUTOFDOMAIN) THEN
C
C          DO K=1,8
C            PATH(GRIDPTR(K,ICELL)) = MIN(TAU,PATH(GRIDPTR(K,ICELL)))
C          ENDDO
C        ENDIF

C               Get the intersection face number (i.e. neighptr index)
        IF (SOX .LE. SOZ .AND. SOX .LE. SOY) THEN
          IFACE = 2-BITX
          JFACE = 1
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),0).AND.BTEST(BCFLAG,0)
        ELSE IF (SOY .LE. SOZ) THEN
          IFACE = 4-BITY
          JFACE = 2
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),1).AND.BTEST(BCFLAG,1)
        ELSE
          IFACE = 6-BITZ
          JFACE = 3
          OPENBCFACE=.FALSE.
        ENDIF
C            Get the next cell to go to
        INEXTCELL = NEIGHPTR(IFACE,ICELL)
        IF (INEXTCELL .LT. 0) THEN
          CALL NEXT_CELL (XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
        ENDIF
C             If going to same or larger face then use previous face
        IF (NEIGHPTR(IFACE,ICELL) .GE. 0 .AND. .NOT.OPENBCFACE) THEN
          KFACE = IFACE
          IC = ICELL
        ELSE
C             If going to smaller face then use next face (more accurate)
          KFACE = OPPFACE(IFACE)
          IC = INEXTCELL
          IFACE = 0
        ENDIF
C           Get the location coordinate
        IF (INEXTCELL .GT. 0) THEN
          IF (JFACE .EQ. 1) THEN
            XN = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
          ELSE IF (JFACE .EQ. 2) THEN
            YN = GRIDPOS(2,GRIDPTR(IOCT,INEXTCELL))
          ELSE
            ZN = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
          ENDIF
        ENDIF

C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
        IF ((INEXTCELL .EQ. 0)) THEN
          DONE = .TRUE.
        ELSE
          XE = XN
          YE = YN
          ZE = ZN
          ICELL = INEXTCELL
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE SPACE_CARVE(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX, VOLUME, FLAGS, WEIGHTS,
     .              LINEAR, COUNTS, DISTANCE_LIMITS)
C    Performs a space carving algorithm by counting intersections between
C    cells and rays. Can also backproject a set of weights. A nearest neighbor
C    interpolation kernel is assumed. The linear one does not operate correctly
C    and needs debugging. - JRLoveridge 2021/02/22
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
Cf2py intent(in) :: NX, NY, NZ, NPTS, NCELLS
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) :: BCFLAG, IPFLAG
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in) :: CELLFLAGS
      DOUBLE PRECISION   CAMX(*), CAMY(*), CAMZ(*), CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX
      INTEGER FLAGS(*)
Cf2py intent(in) :: FLAG
      DOUBLE PRECISION WEIGHTS(*)
Cf2py intent(in) :: WEIGHTS
      DOUBLE PRECISION DISTANCE_LIMITS(*)
Cf2py intent(in) :: DISTANCE_LIMITS
      LOGICAL LINEAR
Cf2py intent(in) LINEAR
      DOUBLE PRECISION  VOLUME(NPTS)
      INTEGER COUNTS(2,NPTS)
Cf2py intent(out):: VOLUME, COUNTS

      INTEGER N, K, FLAG
      DOUBLE PRECISION    MURAY, PHIRAY, MU2, PHI2, WEIGHT
      DOUBLE PRECISION X0, Y0, Z0, R, PI, DISTANCE_LIMIT

      PI = ACOS(-1.0D0)
      VOLUME = 0.0D0
      COUNTS = 0
C         Loop over pixels in image
      DO N = 1, NPIX
        X0 = CAMX(N)
        Y0 = CAMY(N)
        Z0 = CAMZ(N)
        MU2 = CAMMU(N)
        PHI2 = CAMPHI(N)
        FLAG = FLAGS(N)
        WEIGHT = WEIGHTS(N)
        DISTANCE_LIMIT = DISTANCE_LIMITS(N)
        MURAY = -MU2
        PHIRAY = PHI2 - PI
C             Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (MURAY .GE. 0.0) THEN
            GOTO 900
          ENDIF
          R = (ZGRID(NZ) - Z0)/MURAY
          X0 = X0 + R*SQRT(1-MURAY**2)*COS(PHIRAY)
          Y0 = Y0 + R*SQRT(1-MURAY**2)*SIN(PHIRAY)
          Z0 = ZGRID(NZ)
        ELSE IF (Z0 .LT. ZGRID(1)) THEN
          WRITE (6,*) 'SPACE_CARVING: Level', Z0, 'below domain',
     .		                      ZGRID(1)
          STOP
        ENDIF

        DISTANCE_LIMIT = DISTANCE_LIMIT - R

        CALL SPACE_CARVE_1RAY(NX, NY, NZ, NPTS, NCELLS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			             X0, Y0, Z0, VOLUME, FLAG, WEIGHT, LINEAR,
     .                   COUNTS, DISTANCE_LIMIT)
900     CONTINUE

      ENDDO

      RETURN
      END


      SUBROUTINE SPACE_CARVE_1RAY(NX, NY, NZ, NPTS, NCELLS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			             X0, Y0, Z0, VOLUME, FLAG, WEIGHT, LINEAR,
     .                   COUNTS, DISTANCE_LIMIT)

C     Performs the space carving by traversing the grid.
C     LINEAR mode is untested.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2, WEIGHT
      DOUBLE PRECISION DISTANCE_LIMIT
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION  VOLUME(NPTS), PATH_SUM
      LOGICAL LINEAR
      INTEGER FLAG, COUNTS(2,NPTS)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, K
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), BCFLAG, IPFLAG
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    XM, YM
      DOUBLE PRECISION XMID,YMID,ZMID, F1(8), F2(8)
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
      DOUBLE PRECISION DISTANCE
      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/

      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)

      DISTANCE = 0.0D0

C         Make the ray direction (opposite to the discrete ordinate direction)
      CX = SQRT(1.0-MURAY**2)*COS(PHIRAY)
      CY = SQRT(1.0-MURAY**2)*SIN(PHIRAY)
      CZ = MURAY
      IF (ABS(CX) .GT. 1.0E-6) THEN
        CXINV = 1.0D0/CX
      ELSE
        CX = 0.0
        CXINV = 1.0E6
      ENDIF
      IF (ABS(CY) .GT. 1.0E-6) THEN
        CYINV = 1.0D0/CY
      ELSE
        CY = 0.0
        CYINV = 1.0E6
      ENDIF
      IF (ABS(CZ) .GT. 1.0E-6) THEN
        CZINV = 1.0D0/CZ
      ELSE
        CZ = 0.0
        CZINV = 1.0E6
      ENDIF

C         Setup for the ray path direction
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
      IF (CZ .LT. 0.0) THEN
        BITZ = 1
      ELSE
        BITZ = 0
      ENDIF
      IOCT = 1 + BITX + 2*BITY + 4*BITZ
      XM = 0.5*(XGRID(1)+XGRID(NX))
      YM = 0.5*(YGRID(1)+YGRID(NY))

C         Start at the desired point, getting the extinction and source there
      XE = X0
      YE = Y0
      ZE = Z0
      CALL LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, XE, YE, ZE, ICELL)

      IFACE = 0
      NGRID = 0

C         Loop until reach a Z boundary or transmission is very small
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
C           Make sure current cell is valid
        IF (ICELL .LE. 0) THEN
          WRITE (6,*)'INTEGRATE_1RAY: ICELL=',ICELL,
     .                MURAY,PHIRAY,XE,YE,ZE
          STOP
        ENDIF
        NGRID = NGRID + 1

C         Interpolate the source and extinction to the current point
C        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XE, YE, ZE, F)

C           This cell is independent pixel if IP mode or open boundary
C             conditions and ray is leaving domain (i.e. not entering)
        IPINX = BTEST(INT(CELLFLAGS(ICELL)),0) .AND.
     .          .NOT. ( BTEST(BCFLAG,0) .AND.
     .          ((CX.GT.0.AND.XE.LT.XM) .OR. (CX.LT.0.AND.XE.GT.XM)) )
        IPINY = BTEST(INT(CELLFLAGS(ICELL)),1) .AND.
     .          .NOT. ( BTEST(BCFLAG,1) .AND.
     .          ((CY.GT.0.AND.YE.LT.YM) .OR. (CY.LT.0.AND.YE.GT.YM)) )

C           Find boundaries of the current cell
C           Find the three possible intersection planes (X,Y,Z)
C             from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9-IOCT,ICELL)
C           Get the distances to the 3 planes and select the closest
C             (always need to deal with the cell that is wrapped)
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
          WRITE (6,*) 'INTEGRATE_1RAY: SO<0  ',
     .      MURAY,PHIRAY,XE,YE,ZE,SO,ICELL
          STOP
        ENDIF
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ
        XMID = XE + SO*CX/2.0
        YMID = YE + SO*CY/2.0
        ZMID = ZE + SO*CZ/2.0
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS,
     .    XE,YE,ZE, F1)
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS,
     .    XN,YN,ZN, F2)

C	If this is not a boundary cell (currently assuming that the bc conditions are open)
        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        OUTOFDOMAIN = .FALSE.
        IF (.NOT. OUTOFDOMAIN) THEN
            DO K=1,8
              IF ((FLAG .EQ. 1) .AND.
     .          (DISTANCE .GE. DISTANCE_LIMIT)) THEN
                COUNTS(1,GRIDPTR(K,ICELL)) =
     .              COUNTS(1,GRIDPTR(K,ICELL)) + 1
              ELSE
!              ELSEIF (FLAG .EQ. 0) THEN
                 COUNTS(2,GRIDPTR(K,ICELL))=
     .                COUNTS(2,GRIDPTR(K,ICELL)) + 1
              ENDIF
              IF (LINEAR) THEN
                VOLUME(GRIDPTR(K,ICELL)) =
     .              VOLUME(GRIDPTR(K,ICELL)) +
     .              WEIGHT*(F1(K)+F2(K))*0.5*SO
              ELSE
                VOLUME(GRIDPTR(K,ICELL)) =
     .              VOLUME(GRIDPTR(K,ICELL)) + WEIGHT
              ENDIF
            ENDDO
        ENDIF

C               Get the intersection face number (i.e. neighptr index)
        IF (SOX .LE. SOZ .AND. SOX .LE. SOY) THEN
          IFACE = 2-BITX
          JFACE = 1
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),0).AND.BTEST(BCFLAG,0)
        ELSE IF (SOY .LE. SOZ) THEN
          IFACE = 4-BITY
          JFACE = 2
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),1).AND.BTEST(BCFLAG,1)
        ELSE
          IFACE = 6-BITZ
          JFACE = 3
          OPENBCFACE=.FALSE.
        ENDIF
C            Get the next cell to go to
        INEXTCELL = NEIGHPTR(IFACE,ICELL)
        IF (INEXTCELL .LT. 0) THEN
          CALL NEXT_CELL (XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
        ENDIF
C             If going to same or larger face then use previous face
        IF (NEIGHPTR(IFACE,ICELL) .GE. 0 .AND. .NOT.OPENBCFACE) THEN
          KFACE = IFACE
          IC = ICELL
        ELSE
C             If going to smaller face then use next face (more accurate)
          KFACE = OPPFACE(IFACE)
          IC = INEXTCELL
          IFACE = 0
        ENDIF
C           Get the location coordinate
        IF (INEXTCELL .GT. 0) THEN
          IF (JFACE .EQ. 1) THEN
            XN = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
          ELSE IF (JFACE .EQ. 2) THEN
            YN = GRIDPOS(2,GRIDPTR(IOCT,INEXTCELL))
          ELSE
            ZN = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
          ENDIF
        ENDIF

        DISTANCE = DISTANCE + SO

C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
C       .OR.(NGRID.GT.MAXCELLSCROSS)
        IF ((INEXTCELL .EQ. 0)) THEN
          DONE = .TRUE.
        ELSE
          XE = XN
          YE = YN
          ZE = ZN
          ICELL = INEXTCELL
        ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE PROJECT(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,WEIGHTS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX,PATH,
     .             PIXEL_INDICES, RAY_WEIGHTS, NRAYS)
C     Performs ray integration (PATH) of a set of WEIGHTS assuming
C     a linear interpolation kernel.
C     - JRLoveridge 2021/02/22
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
Cf2py intent(in) :: NX, NY, NZ, NPTS, NCELLS
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) :: BCFLAG, IPFLAG
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in) :: CELLFLAGS
      DOUBLE PRECISION  CAMX(NRAYS), CAMY(NRAYS), CAMZ(NRAYS)
      DOUBLE PRECISION  CAMMU(NRAYS), CAMPHI(NRAYS)
Cf2py intent(in) :: CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, PIXEL_INDICES(NRAYS), NRAYS
Cf2py intent(in) NPIX, PIXEL_INDICES, NRAYS
      DOUBLE PRECISION  PATH(NPIX)
      DOUBLE PRECISION RAY_WEIGHTS(NRAYS)
Cf2py intent(in) :: RAY_WEIGHTS
Cf2py intent(out):: PATH
      DOUBLE PRECISION WEIGHTS(NPTS)
Cf2py intent(in) :: WEIGHTS

      DOUBLE PRECISION TAU
      DOUBLE PRECISION RAY_WEIGHT
      INTEGER N, K,Q
      DOUBLE PRECISION    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, R, PI

      PI = ACOS(-1.0D0)
      PATH = 0.0
C      COUNTER = 1
C         Loop over pixels in image
      DO Q=1, NRAYS
        N = PIXEL_INDICES(Q) + 1
        X0 = CAMX(Q)
        Y0 = CAMY(Q)
        Z0 = CAMZ(Q)
        MU2 = CAMMU(Q)
        PHI2 = CAMPHI(Q)
        RAY_WEIGHT = RAY_WEIGHTS(Q)

        MURAY = -MU2
        PHIRAY = PHI2 - PI

C             Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (MURAY .GE. 0.0) THEN
            GOTO 900
          ENDIF
          R = (ZGRID(NZ) - Z0)/MURAY
          X0 = X0 + R*SQRT(1-MURAY**2)*COS(PHIRAY)
          Y0 = Y0 + R*SQRT(1-MURAY**2)*SIN(PHIRAY)
          Z0 = ZGRID(NZ)
        ELSE IF (Z0 .LT. ZGRID(1)) THEN
          WRITE (6,*) 'OPTICAL_DEPTH: Level', Z0, 'below domain',
     .		                      ZGRID(1)
          STOP
        ENDIF
        TAU = 0.0D0
        CALL PROJECT_1RAY(NX, NY, NZ, NPTS, NCELLS, WEIGHTS,
     .                      GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                      BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                      GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                X0, Y0, Z0, TAU ,N,
     .                      RAY_WEIGHT, NPIX, .FALSE.)
        PATH(N) = PATH(N) + TAU
900     CONTINUE

      ENDDO

      RETURN
      END

      SUBROUTINE PROJECT_1RAY(NX, NY, NZ, NPTS, NCELLS, WEIGHTS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                 X0, Y0, Z0, TAU, N,
     .                       RAY_WEIGHT, NPIX, FIRSTMODE)

C       Integrates the source function through the extinction field
C     (EXTINCT) backward in the direction (MURAY,PHIRAY) to find the
C     outgoing radiance (RAD) at the point X0,Y0,Z0.
      IMPLICIT NONE
      LOGICAL FIRSTMODE
      INTEGER NX, NY, NZ, NPTS, NCELLS
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      DOUBLE PRECISION    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION    WEIGHTS(NPTS)
      DOUBLE PRECISION RAY_WEIGHT
      INTEGER NPIX, COUNTER, N
      DOUBLE PRECISION    TAU, EXT1, EXTN
C      REAL    PATH(*)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, K
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), BCFLAG, IPFLAG
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    XM, YM
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, F1(8), F2(8)

      REAL FM
      INTEGER J

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/

      COUNTER = 1
      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)
      TAU = 0.0

C         Make the ray direction (opposite to the discrete ordinate direction)
      CX = SQRT(1.0-MURAY**2)*COS(PHIRAY)
      CY = SQRT(1.0-MURAY**2)*SIN(PHIRAY)
      CZ = MURAY
      IF (ABS(CX) .GT. 1.0E-6) THEN
        CXINV = 1.0D0/CX
      ELSE
        CX = 0.0
        CXINV = 1.0E6
      ENDIF
      IF (ABS(CY) .GT. 1.0E-6) THEN
        CYINV = 1.0D0/CY
      ELSE
        CY = 0.0
        CYINV = 1.0E6
      ENDIF
      IF (ABS(CZ) .GT. 1.0E-6) THEN
        CZINV = 1.0D0/CZ
      ELSE
        CZ = 0.0
        CZINV = 1.0E6
      ENDIF

C         Setup for the ray path direction
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
      IF (CZ .LT. 0.0) THEN
        BITZ = 1
      ELSE
        BITZ = 0
      ENDIF
      IOCT = 1 + BITX + 2*BITY + 4*BITZ
      XM = 0.5*(XGRID(1)+XGRID(NX))
      YM = 0.5*(YGRID(1)+YGRID(NY))

C         Start at the desired point, getting the extinction and source there
      XE = X0
      YE = Y0
      ZE = Z0
      CALL LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, XE, YE, ZE, ICELL)

      IFACE = 0
      NGRID = 0

C         Loop until reach a Z boundary or transmission is very small
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
C           Make sure current cell is valid
        IF (ICELL .LE. 0) THEN
          WRITE (6,*)'INTEGRATE_1RAY: ICELL=',ICELL,
     .                MURAY,PHIRAY,XE,YE,ZE
          STOP
        ENDIF
        NGRID = NGRID + 1

C         Interpolate the source and extinction to the current point
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XE, YE, ZE, F1)
        EXT1 = 0.0
        DO J=1,8
          EXT1 = EXT1 + WEIGHTS(GRIDPTR(J,ICELL))*F1(J)
        ENDDO
C           This cell is independent pixel if IP mode or open boundary
C             conditions and ray is leaving domain (i.e. not entering)
        IPINX = BTEST(INT(CELLFLAGS(ICELL)),0) .AND.
     .          .NOT. ( BTEST(BCFLAG,0) .AND.
     .          ((CX.GT.0.AND.XE.LT.XM) .OR. (CX.LT.0.AND.XE.GT.XM)) )
        IPINY = BTEST(INT(CELLFLAGS(ICELL)),1) .AND.
     .          .NOT. ( BTEST(BCFLAG,1) .AND.
     .          ((CY.GT.0.AND.YE.LT.YM) .OR. (CY.LT.0.AND.YE.GT.YM)) )

C           Find boundaries of the current cell
C           Find the three possible intersection planes (X,Y,Z)
C             from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9-IOCT,ICELL)
C           Get the distances to the 3 planes and select the closest
C             (always need to deal with the cell that is wrapped)
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
          WRITE (6,*) 'INTEGRATE_1RAY: SO<0  ',
     .      MURAY,PHIRAY,XE,YE,ZE,SO,ICELL
          STOP
        ENDIF
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ
C           Find the optical path across the grid cell and figure how
C             many subgrid intervals to use
C         Interpolate the source and extinction to the current point
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XN, YN, ZN, F2)
        EXTN = 0.0
        DO J=1,8
          EXTN = EXTN+ WEIGHTS(GRIDPTR(J,ICELL))*F2(J)
        ENDDO
        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        OUTOFDOMAIN=.FALSE.
        IF (.NOT.OUTOFDOMAIN) THEN
          IF (FIRSTMODE) THEN
            TAU = MAXVAL(WEIGHTS(GRIDPTR(:,ICELL)))*RAY_WEIGHT
            IF (TAU .GT. 1e-5 .AND. NGRID .GE. 1) DONE = .TRUE.
          ELSE
            TAU = TAU + SO*0.5*(EXT1+EXTN)*RAY_WEIGHT
          ENDIF
        ENDIF

C               Get the intersection face number (i.e. neighptr index)
        IF (SOX .LE. SOZ .AND. SOX .LE. SOY) THEN
          IFACE = 2-BITX
          JFACE = 1
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),0).AND.BTEST(BCFLAG,0)
        ELSE IF (SOY .LE. SOZ) THEN
          IFACE = 4-BITY
          JFACE = 2
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),1).AND.BTEST(BCFLAG,1)
        ELSE
          IFACE = 6-BITZ
          JFACE = 3
          OPENBCFACE=.FALSE.
        ENDIF
C            Get the next cell to go to
        INEXTCELL = NEIGHPTR(IFACE,ICELL)
        IF (INEXTCELL .LT. 0) THEN
          CALL NEXT_CELL (XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
        ENDIF
C             If going to same or larger face then use previous face
        IF (NEIGHPTR(IFACE,ICELL) .GE. 0 .AND. .NOT.OPENBCFACE) THEN
          KFACE = IFACE
          IC = ICELL
        ELSE
C             If going to smaller face then use next face (more accurate)
          KFACE = OPPFACE(IFACE)
          IC = INEXTCELL
          IFACE = 0
        ENDIF
C           Get the location coordinate
        IF (INEXTCELL .GT. 0) THEN
          IF (JFACE .EQ. 1) THEN
            XN = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
          ELSE IF (JFACE .EQ. 2) THEN
            YN = GRIDPOS(2,GRIDPTR(IOCT,INEXTCELL))
          ELSE
            ZN = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
          ENDIF
        ENDIF

C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
        IF ((INEXTCELL .EQ. 0).OR.(NGRID.GT.MAXCELLSCROSS)) THEN
          DONE = .TRUE.
        ELSE
          XE = XN
          YE = YN
          ZE = ZN
          ICELL = INEXTCELL
        ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE GET_SHADOW(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,WEIGHTS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX,PATH,
     .             PIXEL_INDICES, RAY_WEIGHTS, NRAYS
     .             )

      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
Cf2py intent(in) :: NX, NY, NZ, NPTS, NCELLS
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) :: BCFLAG, IPFLAG
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in) :: CELLFLAGS
      DOUBLE PRECISION SOLARMU, SOLARAZ
Cf2py intent(in) :: SOLARMU, SOLARAZ
      DOUBLE PRECISION  CAMX(NRAYS), CAMY(NRAYS), CAMZ(NRAYS)
      DOUBLE PRECISION  CAMMU(NRAYS), CAMPHI(NRAYS)
Cf2py intent(in) :: CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, PIXEL_INDICES(NRAYS), NRAYS
Cf2py intent(in) NPIX, PIXEL_INDICES, NRAYS
      DOUBLE PRECISION  PATH(NPIX)
      DOUBLE PRECISION RAY_WEIGHTS(NRAYS)
Cf2py intent(in) :: RAY_WEIGHTS
Cf2py intent(out):: PATH
      DOUBLE PRECISION WEIGHTS(NPTS)
Cf2py intent(in) :: WEIGHTS

      DOUBLE PRECISION TAU
      DOUBLE PRECISION RAY_WEIGHT
      INTEGER N, K,Q
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, R, PI

      PI = ACOS(-1.0D0)

      DO Q=1, NRAYS
        N = PIXEL_INDICES(Q) + 1
        X0 = CAMX(Q)
        Y0 = CAMY(Q)
        Z0 = CAMZ(Q)
        MU2 = CAMMU(Q)
        PHI2 = CAMPHI(Q)
        RAY_WEIGHT = RAY_WEIGHTS(Q)

        MURAY = -MU2
        PHIRAY = PHI2 - PI

C             Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (MURAY .GE. 0.0) THEN
            GOTO 900
          ENDIF
          R = (ZGRID(NZ) - Z0)/MURAY
          X0 = X0 + R*SQRT(1-MURAY**2)*COS(PHIRAY)
          Y0 = Y0 + R*SQRT(1-MURAY**2)*SIN(PHIRAY)
          Z0 = ZGRID(NZ)
        ELSE IF (Z0 .LT. ZGRID(1)) THEN
          WRITE (6,*) 'OPTICAL_DEPTH: Level', Z0, 'below domain',
     .		                      ZGRID(1)
          STOP
        ENDIF
        TAU = 0.0D0
        CALL PROJECT_1RAY(NX, NY, NZ, NPTS, NCELLS, WEIGHTS,
     .                      GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                      BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                      GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                X0, Y0, Z0, TAU ,N,
     .                      RAY_WEIGHT, NPIX, .TRUE.)
        PATH(N) = PATH(N) + TAU
900   CONTINUE
      ENDDO
      RETURN
      END


      SUBROUTINE RAYLEIGH_EXTINCT (NZT, ZLEVELS,TEMP, RAYSFCPRES,
     .                             RAYLCOEF, EXTRAYL)
C       Computes the molecular Rayleigh extinction profile EXTRAYL [/km]
C     from the temperature profile TEMP [K] at ZLEVELS [km].  Assumes
C     a linear lapse rate between levels to compute the pressure at
C     each level.  The Rayleigh extinction is proportional to air
C     density, with the coefficient RAYLCOEF in [K/(mb km)].
      IMPLICIT NONE
      INTEGER NZT
Cf2py intent(in) :: NZT
      REAL    ZLEVELS(NZT), TEMP(NZT), RAYSFCPRES, RAYLCOEF
      REAL    EXTRAYL(NZT)
Cf2py intent(in) :: ZLEVELS, TEMP, RAYSFCPRES, RAYLCOEF
Cf2py intent(out) :: EXTRAYL
      INTEGER I
      REAL    PRES, LAPSE, TS, DZ

C           Find surface pressure by integrating hydrostatic relation
C           for a dry atmosphere up to surface height.
      PRES = RAYSFCPRES
      TS = TEMP(1)
      LAPSE = 6.5*0.001
      PRES = PRES*(TS/(TS+LAPSE*ZLEVELS(1)*1000.))**(9.8/(287.*LAPSE))

C         Use layer mean temperature to compute fractional pressure change.
      DO I = 1, NZT-1
        EXTRAYL(I) = RAYLCOEF*PRES/TEMP(I)
        DZ = 1000.*(ZLEVELS(I+1)-ZLEVELS(I))
        LAPSE = (TEMP(I)-TEMP(I+1))/DZ
        IF (ABS(LAPSE) .GT. 0.00001) THEN
          PRES = PRES*(TEMP(I+1)/TEMP(I))**(9.8/(287.*LAPSE))
        ELSE
          PRES = PRES*EXP(-9.8*DZ/(287.*TEMP(I)))
        ENDIF
      ENDDO
      EXTRAYL(NZT) = RAYLCOEF*PRES/TEMP(NZT)

      RETURN
      END

      SUBROUTINE GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, X, Y, Z, F)
C     Compute trilinear interpolation kernel F(8)
      IMPLICIT NONE
      INTEGER ICELL, BCELL, IPT1, IPT2
      INTEGER GRIDPTR(8,*)
      REAL GRIDPOS(3,*)
      DOUBLE PRECISION U, V, W, F(8), DELX, DELY, DELZ
      DOUBLE PRECISION X, Y, Z,INVDELX,INVDELY,INVDELZ

      IPT1 = GRIDPTR(1,ICELL)
      IPT2 = GRIDPTR(8,ICELL)
      DELX = GRIDPOS(1,IPT2) - GRIDPOS(1,IPT1)
      IF (DELX .LE. 0.0) THEN
        INVDELX = 1.0
      ELSE
        INVDELX = 1.0D0/DELX
      ENDIF
      DELY = GRIDPOS(2,IPT2) - GRIDPOS(2,IPT1)
      IF (DELY .LE. 0.0) THEN
        INVDELY = 1.0
      ELSE
        INVDELY = 1.0D0/DELY
      ENDIF
      DELZ = GRIDPOS(3,IPT2) - GRIDPOS(3,IPT1)
      INVDELZ = 1.0D0/DELZ
      U = (X-GRIDPOS(1,IPT1))*INVDELX
      V = (Y-GRIDPOS(2,IPT1))*INVDELY
      W = (Z-GRIDPOS(3,IPT1))*INVDELZ

      F(1) = (1-W)*(1-V)*(1-U)
      F(2) = (1-W)*(1-V)*U
      F(3) = (1-W)*V*(1-U)
      F(4) = (1-W)*V*U
      F(5) = W*(1-V)*(1-U)
      F(6) = W*(1-V)*U
      F(7) = W*V*(1-U)
      F(8) = W*V*U

      RETURN
      END


      SUBROUTINE GET_BASE_GRID_CELL(BCELL, ICELL, TREEPTR)
      IMPLICIT NONE
      INTEGER BCELL, ICELL, TREEPTR(2,*)

      BCELL = ICELL
      DO WHILE (TREEPTR(1, BCELL) .GT. 0)
        BCELL = TREEPTR(1, BCELL)
      ENDDO

      RETURN
      END



      SUBROUTINE MAKE_DIRECT_DERIVATIVE (NPTS, BCFLAG, NPX, NPY, NPZ,
     .		         DELX, DELY, XSTART, YSTART, GRIDPOS, ZLEVELS,
     .		         IPDIRECT, DI, DJ, DK, CX, CY, CZ, CXINV,
     .		         CYINV, CZINV, EPSS, EPSZ, XDOMAIN, YDOMAIN,
     .		         UNIFORMZLEV, DELXD, DELYD, DPATH, DPTR,
     .             IERR, ERRMSG, LONGEST_PATH_PTS)
C     Calculates the sensitvity of the direct beam to the extinction
C     along the path from each grid point to the sun.
C     Actually calls DIRECT_BEAM_AND_PATHS_PROP to do all the hard work.
      IMPLICIT NONE
      INTEGER NPTS, BCFLAG
Cf2py intent(in) :: NPTS, BCFLAG
      INTEGER NPX, NPY, NPZ
Cf2py intent(in) :: NPX, NPY, NPZ
      REAL DELX, DELY, XSTART, YSTART, GRIDPOS(3,*), ZLEVELS(*)
Cf2py intent(in) :: DELX, DELY, XSTART, YSTART, GRIDPOS, ZLEVELS
      INTEGER IPDIRECT, DI, DJ, DK
Cf2py intent(in) :: IPDIRECT, DI, DJ, DK
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
Cf2py intent(in) ::   CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION EPSS, EPSZ, XDOMAIN, YDOMAIN
Cf2py intent(in) :: EPSS, EPSZ, XDOMAIN, YDOMAIN
      DOUBLE PRECISION UNIFORMZLEV, DELXD, DELYD
Cf2py intent(in) :: UNIFORMZLEV, DELXD, DELYD
      INTEGER LONGEST_PATH_PTS
Cf2py intent(in) :: LONGEST_PATH_PTS
      REAL DPATH(LONGEST_PATH_PTS,NPTS)
      INTEGER DPTR(LONGEST_PATH_PTS,NPTS)
Cf2py intent(out) :: DPATH, DPTR
      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG

      INTEGER SIDE, IP
      LOGICAL VALIDBEAM
      REAL    UNIFZLEV, XO, YO, ZO, DIR, DIRPATH
      IERR = 0
      DPTR = 0
      DPATH = 0.0

      DO IP = 1, NPTS
        DIRPATH = 0.0
        CALL DIRECT_BEAM_AND_PATHS_PROP(GRIDPOS(1,IP), GRIDPOS(2,IP),
     .            GRIDPOS(3,IP), BCFLAG, XO, YO, ZO, SIDE, NPX, NPY,
     .            NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS, CX, CY, CZ,
     .            CXINV, CYINV, CZINV, DI, DJ, DK, IPDIRECT, DELXD,
     .            DELYD, XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV,
     .            DPATH(:,IP), DPTR(:,IP), IERR, ERRMSG,
     .            LONGEST_PATH_PTS)
        IF (IERR .NE. 0) RETURN
      ENDDO
      RETURN
      END


      SUBROUTINE DIRECT_BEAM_AND_PATHS_PROP (XI, YI, ZI, BCFLAG,
     .           XO, YO, ZO, SIDE, NPX, NPY, NPZ, DELX, DELY,
     .           XSTART, YSTART, ZLEVELS, CX, CY, CZ, CXINV, CYINV,
     .           CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD, XDOMAIN,
     .           YDOMAIN, EPSS, EPSZ, UNIFORMZLEV, DPATH, DPTR,
     .           IERR, ERRMSG, LONGEST_PATH_PTS)
      IMPLICIT NONE

      REAL    DPATH(LONGEST_PATH_PTS)
      INTEGER DPTR(LONGEST_PATH_PTS)
      INTEGER BCFLAG, SIDE
      REAL    XO, YO, ZO, XI, YI, ZI
      INTEGER IERR, LONGEST_PATH_PTS
      CHARACTER ERRMSG*600

      INTEGER IX, IY, IZ, JZ, IL, IM, IU
      INTEGER I, J, K, IP, JP, I1, I2, I3, I4
      INTEGER IPDIRECT, DI, DJ, DK
      LOGICAL CONSTX, CONSTY, HITBOUNDARY, BTEST, OUTOFDOMAIN
      DOUBLE PRECISION UNIFORMZLEV
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION EPSS, EPSZ, XDOMAIN, YDOMAIN
      DOUBLE PRECISION XOFFS, YOFFS, DELXD, DELYD
      DOUBLE PRECISION X, Y, Z, XE, YE, ZE, XP, YP, ZP
      DOUBLE PRECISION X0, X1, Y0, Y1, Z0, Z1, SO, SOX, SOY, SOZ
      DOUBLE PRECISION U0, V0, W0, U1, V1, W1, AX, AY, AZ
      DOUBLE PRECISION U0M, V0M, W0M, U1M, V1M, W1M, DU, DV, DW
      DOUBLE PRECISION B1,B2,B3,B4,B5,B6,B7,B8,C1,C2,C3,C4,C5,C6,C7,C8
      DOUBLE PRECISION A1, A2, A3, A4, A5, A6, A7, A8
      DOUBLE PRECISION UV,UMV,UVM,UMVM,UW,UMW,UWM,UMWM,VW,VMW,VWM,VMWM
      DOUBLE PRECISION VWU,VWUM,UWV,UWVM,UVW,UVWM

      INTEGER NPX, NPY, NPZ, IDP
      REAL DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)

      OUTOFDOMAIN = .FALSE.

C         Here for computing the direct beam path for one starting point.
      Z = ZI
      X = XI - XSTART
      Y = YI - YSTART

C         Find the grid location
      IL=0
      IU=NPZ
      DO WHILE (IU-IL .GT. 1)
        IM = (IU+IL)/2
        IF (Z .GE. ZLEVELS(IM)) THEN
          IL = IM
        ELSE
          IU=IM
        ENDIF
      ENDDO
      K = MAX(IL,1)


      I = INT(X/DELXD) + 1
      IF (I .GT. NPX .AND. ABS(X-XDOMAIN) .LT. 0.001*DELXD) I = NPX
      IF (I .LT. 1 .OR. I .GT. NPX) THEN
        IERR = 1
        WRITE (ERRMSG,*) 'DIRECT_BEAM_AND_PATHS_PROP: Beyond X domain',
     .    I,XI,YI,ZI
        RETURN
      ENDIF
      J = INT(Y/DELYD) + 1
      IF (J .GT. NPY .AND. ABS(Y-YDOMAIN) .LT. 0.001*DELYD) J = NPY
      IF (J .LT. 1 .OR. J .GT. NPY) THEN
        IERR = 1
        WRITE (6,*) 'DIRECT_BEAM_AND_PATHS_PROP: Beyond Y domain',
     .    J,XI,YI,ZI
        RETURN
      ENDIF
      XE = X
      YE = Y
      ZE = Z
      XP = XE
      YP = YE
      ZP = ZE
      CONSTX = BTEST(IPDIRECT,0)
      CONSTY = BTEST(IPDIRECT,1)
      IF (CX .EQ. 0.0)  CONSTX = .TRUE.
      IF (CY .EQ. 0.0)  CONSTY = .TRUE.
      IF (BTEST(BCFLAG,0) .AND. (ABS(X) .LT. 0.01*DELXD
     .   .OR. ABS(X-(NPX-1)*DELXD) .LT. 0.01*DELXD)) THEN
        CONSTX = .TRUE.
        OUTOFDOMAIN = .TRUE.
      ENDIF
      IF (BTEST(BCFLAG,1) .AND. (ABS(Y) .LT. 0.01*DELYD
     .   .OR. ABS(Y-(NPY-1)*DELYD) .LT. 0.01*DELYD)) THEN
        CONSTY = .TRUE.
        OUTOFDOMAIN = .TRUE.
      ENDIF

C     If have multiple subdomains (processors) and ray is going outwards
C      from a boundary then set the HITBOUNDARY flag
      HITBOUNDARY = .FALSE.
      IF (BTEST(BCFLAG,2)) THEN
        IF (CX .GT. 0.0 .AND. ABS(X-XDOMAIN) .LT. 0.001*DELXD) THEN
          SIDE = 2
          HITBOUNDARY = .TRUE.
        ELSE IF (CX .LT. 0.0 .AND. ABS(X) .LT.  0.001*DELXD) THEN
          SIDE = 1
          HITBOUNDARY = .TRUE.
        ENDIF
      ENDIF
      IF (BTEST(BCFLAG,3)) THEN
        IF (CY .GT. 0.0 .AND. ABS(Y-YDOMAIN) .LT. 0.001*DELYD) THEN
          SIDE = 4
          HITBOUNDARY = .TRUE.
        ENDIF
        IF (CY .LT. 0.0 .AND. ABS(Y) .LT.  0.001*DELYD) THEN
          SIDE = 3
          HITBOUNDARY = .TRUE.
        ENDIF
      ENDIF
C           Grid cell loop begin
      IDP = 0
      DO WHILE (.NOT. HITBOUNDARY .AND. ABS(ZE-ZLEVELS(NPZ)) .GT. EPSZ)
C        PRINT *, I,J,K, XE, YE,ZE
        IP = I + 1
        IF (I .EQ. NPX) THEN
          IF (BTEST(BCFLAG,0) .OR. BTEST(BCFLAG,2)) THEN
            IP=NPX
          ELSE
            IP = 1
          ENDIF
        ENDIF
        JP = J + 1
        IF (J .EQ. NPY) THEN
          IF (BTEST(BCFLAG,1) .OR. BTEST(BCFLAG,3)) THEN
            JP=NPY
          ELSE
            JP = 1
          ENDIF
        ENDIF
        X0 = DELXD*(I-1)
        X1 = X0 + DELXD
        Y0 = DELYD*(J-1)
        Y1 = Y0 + DELYD
        Z0 = ZLEVELS(K)
        Z1 = ZLEVELS(K+1)

        IF (I .LT. 1 .OR. I .GT. NPX .OR.
     .      J .LT. 1 .OR. J .GT. NPY .OR.
     .      K .LT. 1 .OR. K .GE. NPZ) THEN
          WRITE (ERRMSG,'(A,3I4)') 'DIRECT_BEAM_AND_PATHS_PROP: ',
     .     'beyond grid!', I, J, K
          WRITE(ERRMSG,*) NPX, NPY, NPZ
          WRITE (ERRMSG,'(1(2X,3F9.5))') X, Y, Z
          RETURN
        ENDIF
C           Get the eight corner extinction values
        I1 = K + NPZ*(J-1)  + NPZ*NPY*(I-1)
        I2 = K + NPZ*(J-1)  + NPZ*NPY*(IP-1)
        I3 = K + NPZ*(JP-1) + NPZ*NPY*(I-1)
        I4 = K + NPZ*(JP-1) + NPZ*NPY*(IP-1)

C           Compute the distance to the next grid plane in  X, Y, and Z
C             If in horizontal uniform region or doing IP then fix X and/or Y.
C        IF (ZE .GE. UNIFORMZLEV) THEN
C          CONSTX = .TRUE.
C          CONSTY = .TRUE.
C        ENDIF
        IF (CONSTX) THEN
          SOX = 1.0E30
        ELSE IF (CX .GT. 0.0) THEN
          SOX = (X1-XE)*CXINV
          XP = X1
        ELSE
          SOX = (X0-XE)*CXINV
          XP = X0
        ENDIF
        IF (CONSTY) THEN
          SOY = 1.0E30
        ELSE IF (CY .GT. 0.0) THEN
          SOY = (Y1-YE)*CYINV
          YP = Y1
        ELSE
          SOY = (Y0-YE)*CYINV
          YP = Y0
        ENDIF
        IF (CZ .GT. 0.0) THEN
          SOZ = (Z1-ZE)*CZINV
          ZP = Z1
        ELSE IF (CZ .LT. 0.0) THEN
          SOZ = (Z0-ZE)*CZINV
          ZP = Z0
        ELSE
          SOZ = 1.0E30
        ENDIF

C           The shortest distance is the plane we stop at:
C             get the exitting location and increment the cell
        XOFFS = 0.0
        YOFFS = 0.0
        IF (SOZ .LE. SOX .AND. SOZ .LE. SOY) THEN
          SO = SOZ
          IF (.NOT. CONSTX) XP = XE + SO*CX
          IF (.NOT. CONSTY) YP = YE + SO*CY
          K = K + DK
        ELSE IF (SOX .LE. SOY) THEN
          SO = SOX
          IF (.NOT. CONSTY) YP = YE + SO*CY
          ZP = ZE + SO*CZ
          I = I + DI
C             If have reached a horizontal boundary then either wrap around
C               (periodic) or go into IP mode (open boundaries).
          IF (I .EQ. 0) THEN
            OUTOFDOMAIN = .TRUE.
            IF (BTEST(BCFLAG,0)) THEN
              I = 1
              CONSTX = .TRUE.
            ELSE IF (BTEST(BCFLAG,2)) THEN
              HITBOUNDARY = .TRUE.
              SIDE = 1
            ELSE
              I = NPX
              XOFFS = XDOMAIN
            ENDIF
          ELSE IF (I .GE. NPX .AND. BTEST(BCFLAG,2)) THEN
              HITBOUNDARY = .TRUE.
              SIDE = 2
          ELSE IF (I .EQ. NPX+1) THEN
            OUTOFDOMAIN = .TRUE.
            IF (BTEST(BCFLAG,0)) THEN
              I = NPX
              CONSTX = .TRUE.
            ELSE
              I = 1
              XOFFS = -XDOMAIN
            ENDIF
          ENDIF
        ELSE
          SO = SOY
          IF (.NOT. CONSTX) XP = XE + SO*CX
          ZP = ZE + SO*CZ
          J = J + DJ
          IF (J .EQ. 0) THEN
            OUTOFDOMAIN = .TRUE.
            IF (BTEST(BCFLAG,1)) THEN
              J = 1
              CONSTY = .TRUE.
            ELSE IF (BTEST(BCFLAG,3)) THEN
              HITBOUNDARY = .TRUE.
              SIDE = 3
            ELSE
              J = NPY
              YOFFS = YDOMAIN
            ENDIF
          ELSE IF (J .GE. NPY .AND. BTEST(BCFLAG,3)) THEN
              HITBOUNDARY = .TRUE.
              SIDE = 4
          ELSE IF (J .EQ. NPY+1) THEN
            OUTOFDOMAIN = .TRUE.
            IF (BTEST(BCFLAG,1)) THEN
              J = NPY
              CONSTY = .TRUE.
            ELSE
              J = 1
              YOFFS = -YDOMAIN
            ENDIF
          ENDIF
        ENDIF
        IF (SO .LT. -EPSS) THEN
          IERR = 1
          WRITE (ERRMSG,*) 'DIRECT_BEAM_AND_PATHS_PROP: SO<0',
     .      X,Y,Z,
     .         XE,YE,ZE, XP,YP,ZP, CX,CY,CZ, SOX,SOY,SOZ
          RETURN
        ENDIF

        SO = MAX(SO,0.0D0)
C        PRINT *, CONSTX, CONSTY, OUTOFDOMAIN
C        PRINT *, SO
C        PRINT *, I,J,K,IP,JP,ZP, ZE

C           Make the starting and ending interpolation factors
        AX = 1.0D0/(X1-X0)
        AY = 1.0D0/(Y1-Y0)
        AZ = 1.0D0/(Z1-Z0)
        U0 = (XE-X0)*AX
        V0 = (YE-Y0)*AY
        W0 = (ZE-Z0)*AZ
        U1 = (XP-X0)*AX
        V1 = (YP-Y0)*AY
        W1 = (ZP-Z0)*AZ
C           Compute the cubic polynomial extinction coefficients
        U0M = 1.0-U0
        V0M = 1.0-V0
        W0M = 1.0-W0
        U1M = 1.0-U1
        V1M = 1.0-V1
        W1M = 1.0-W1
        DU = U1-U0
        DV = V1-V0
        DW = W1-W0
        UV   = U0 *V0
        UMV  = U0M*V0
        UVM  = U0 *V0M
        UMVM = U0M*V0M
        UW   = U0 *W0
        UMW  = U0M*W0
        UWM  = U0 *W0M
        UMWM = U0M*W0M
        VW   = V0 *W0
        VMW  = V0M*W0
        VWM  = V0 *W0M
        VMWM = V0M*W0M

C       Unlike DIRECT_BEAM_PROP calculate A1 - A8
C       as VW is reused in the C1 - C8 calculations
C       leading to errors.
        A1 = U0M*VMWM
        A2 = U0*VMWM
        A3 = U0M*VWM
        A4 = U0*VWM
        A5 = U0M*VMW
        A6 = U0*VMW
        A7 = U0M*VW
        A8 = U0*VW

        B1 = -DU*VMWM - DV*UMWM - DW*UMVM
        B2 =  DU*VMWM - DV*UWM  - DW*UVM
        B3 = -DU*VWM  + DV*UMWM - DW*UMV
        B4 =  DU*VWM  + DV*UWM  - DW*UV
        B5 = -DU*VMW  - DV*UMW  + DW*UMVM
        B6 =  DU*VMW  - DV*UW   + DW*UVM
        B7 = -DU*VW   + DV*UMW  + DW*UMV
        B8 =  DU*VW   + DV*UW   + DW*UV

        VW = DV*DW
        VWU  = VW*U0
        VWUM = VW*U0M
        UW = DU*DW
        UWV  = UW*V0
        UWVM = UW*V0M
        UV = DU*DV
        UVW  = UV*W0
        UVWM = UV*W0M
        C1 = + VWUM + UWVM + UVWM
        C2 = + VWU  - UWVM - UVWM
        C3 = - VWUM + UWV  - UVWM
        C4 = - VWU  - UWV  + UVWM
        C5 = - VWUM - UWVM + UVW
        C6 = - VWU  + UWVM - UVW
        C7 = + VWUM - UWV  - UVW
        C8 = + VWU  + UWV  + UVW


C       Compute the inner derivatives (minus the path lengths)
C .NOT. OUTOFDOMAIN .AND.
C        IF (K .LE. NPZ+1) THEN
          IF (IDP+8 .GT. LONGEST_PATH_PTS) THEN
            IERR = 1
            WRITE(ERRMSG,*) 'DIRECT_BEAM_AND_PATHS_PROP: ',
     .      'Max number of property points to pass',
     .       'exceeded: IDP=', LONGEST_PATH_PTS
            RETURN
          ENDIF
        DPATH(IDP+1) = SO*(A1 + 0.5D0*B1 +
     .                 0.3333333333333333D0*C1 - 0.25D0*DU*DV*DW)
        DPATH(IDP+2) = SO*(A2 + 0.5D0*B2 +
     .                 0.3333333333333333D0*C2 + 0.25D0*DU*DV*DW)
        DPATH(IDP+3) = SO*(A3 + 0.5D0*B3 +
     .                 0.3333333333333333D0*C3 + 0.25D0*DU*DV*DW)
        DPATH(IDP+4) = SO*(A4 + 0.5D0*B4 +
     .                 0.3333333333333333D0*C4 - 0.25D0*DU*DV*DW)
        DPATH(IDP+5) = SO*(A5 + 0.5D0*B5 +
     .                 0.3333333333333333D0*C5 + 0.25D0*DU*DV*DW)
        DPATH(IDP+6) = SO*(A6 + 0.5D0*B6 +
     .                 0.3333333333333333D0*C6 - 0.25D0*DU*DV*DW)
        DPATH(IDP+7) = SO*(A7 + 0.5D0*B7 +
     .                 0.3333333333333333D0*C7 - 0.25D0*DU*DV*DW)
        DPATH(IDP+8) = SO*(A8 + 0.5D0*B8 +
     .                 0.3333333333333333D0*C8 + 0.25D0*DU*DV*DW)

        DPTR(IDP+1) = I1
        DPTR(IDP+2) = I2
        DPTR(IDP+3) = I3
        DPTR(IDP+4) = I4
        DPTR(IDP+5) = I1+1
        DPTR(IDP+6) = I2+1
        DPTR(IDP+7) = I3+1
        DPTR(IDP+8) = I4+1

        IDP = IDP + 8
C        ENDIF

        XE = XP + XOFFS
        YE = YP + YOFFS
        ZE = ZP

      ENDDO

      RETURN
      END

C    Store the source in sparse format.
C
C    Group input rays by angle, each with its magnitude.
C    Get a GIANT array of size MAXCELLSCROSS*8*NANGLES
C    For each angle:
C     Fill the GIANT arrays with the adjoint beams and pointers.
C
C    SORT the array by gridpoint, taking the values
C    and ANGLES along with it.
C    Form a pointer array for the start index of each grid point in
C    the data and angle arrays.
C    Set a flag for whether to precompute the YLMDIR for all of the angles.
C    If so:
C      Do the precomputing.
C    MODIFY CALC_SOURCE_PNT/CALC_SOURCE_PNT_UNPOL
C    to accomodate the new source.

C    Remove negativity check on radiance/SRCEXT during solution process
C    or make it optional.
C
C    Add a new SRCTYPE option and revise all of the previous conditions
C    because they assume only S,T,B.

C      SUBROUTINE MAKE_ADJOINT (CAMX, CAMY, CAMZ, CAMMU, CAMPHI,
C     .    CAMMAGNITUDE, ADJSOURCE, TOTAL_EXT, NX, NY, NZ, NCELLS,
C     .    NPTS, CELLFLAGS, XGRID, YGRID, ZGRID, GRIDPOS,
C     .    GRIDPTR, NEIGHPTR, TREEPTR, BCFLAG, IPFLAG,
C     .    NVIS, MAXCELLSCROSS, IERR, ERRMSG, ADJSOURCEPTR,
C     .    TAUTOL, ANGLEPTR)

C     Angle sorting must have been done in python.
C     Loop through rays
C       Update the sparse/big array with the non-zero entries of ADJTEMP
C       reset ADJTEMP when the angle changes.
C
C     SORT the sparse/big arrays to the gridpoint first sparse format.
C      IMPLICIT NONE
C      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) BCFLAG, IPFLAG
C      INTEGER NX, NY, NZ, NPTS, NCELLS, NVIS
Cf2py intent(in) NX, NY, NZ, NPTS, NCELLS, NVIS
C      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
Cf2py intent(in) GRIDPTR, NEIGHPTR, TREEPTR
C      INTEGER*2 CELLFLAGS(NCELLS)
Cf2py intent(in) CELLFLAGS
C      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
Cf2py intent(in) XGRID, YGRID, ZGRID, GRIDPOS
C      REAL CAMX(NVIS), CAMY(NVIS), CAMZ(NVIS), TOTAL_EXT(NPTS)
C      DOUBLE PRECISION CAMMU(NVIS), CAMPHI(NVIS), TAUTOL
Cf2py intent(in) :: TAUTOL, TOTAL_EXT
C      REAL CAMMAGNITUDE(NVIS)
Cf2py intent(in) CAMX, CAMY, CAMZ, CAMMU, CAMPHI, CAMMAGNITUDE
C      INTEGER IERR, MAXCELLSCROSS
Cf2py intent(in) MAXCELLSCROSS
Cf2py intent(out) IERR
C      CHARACTER ERRMSG*600
Cf2py intent(out) :: ERRMSG
C      REAL ADJSOURCE(8*MAXCELLSCROSS*NVIS)
C      INTEGER ANGLEPTRS(8*MAXCELLSCROSS*NVIS)
C      INTEGER PTSPTRS(8*MAXCELLSCROSS*NVIS)
C      INTEGER ADJSOURCEPTR(NPTS+1)
Cf2py intent(out) ADJSOURCE, ANGLEPTRS, ADJSOURCEPTR
C      INTEGER I, IVIS, NSOURCE, NSOURCETEMP, ANGLEPTR
Cf2py intent(out) ANGLEPTR
C      DOUBLE PRECISION MU, MUOLD, PHI, PHIOLD, MURAY,PHIRAY
C      DOUBLE PRECISION X0, Y0, Z0, PI, TRANSMIT, XE,YE,ZE,R
C      INTEGER SIDE, PTR, PTROLD, COUNTER
C      REAL ADJTEMP(NPTS), MAGNITUDE
C
C
C      PI = ACOS(-1.0D0)
C      IERR = 0
C      ADJTEMP = 0.0
C      NSOURCE = 1
C      ANGLEPTR = 1
C      MUOLD = CAMMU(1)
C      PHIOLD = CAMPHI(1)
C
C      DO IVIS=1,NVIS
C        MU = CAMMU(IVIS)
C        PHI = CAMPHI(IVIS)
C       Refresh the temporary adjoint source storage and
C       update the full storage if the angle has changed.
C        IF (ABS(MU - MUOLD) .GE. 1e-6 .OR. ABS(PHI - PHIOLD)
C     .    .GE. 1e-6) THEN
C
C          NSOURCETEMP = 1
C          DO I=1,NPTS
C            IF (ABS(ADJTEMP(I)) .GE. 1e-6) THEN
C             Update the coordinate list format sparse representation
C             of the adjoint source.
C              ADJSOURCE(NSOURCE + NSOURCETEMP) = ADJTEMP(I)
C              ANGLEPTRS(NSOURCE + NSOURCETEMP) = ANGLEPTR
C              PTSPTRS(NSOURCE + NSOURCETEMP) = I
C              NSOURCETEMP = NSOURCETEMP + 1
C            ENDIF
C          ENDDO
C          NSOURCE = NSOURCE + NSOURCETEMP
C          ANGLEPTR = ANGLEPTR + 1
C          ADJTEMP = 0.0
C        ENDIF
C        MUOLD = MU
C        PHIOLD = PHI
C
C        X0 = CAMX(IVIS)
C        Y0 = CAMY(IVIS)
C        Z0 = CAMZ(IVIS)
C        MAGNITUDE = CAMMAGNITUDE(IVIS)
C        MURAY = -MU
C        PHIRAY = PHI - PI

C             Extrapolate ray to domain top if above
C        IF (Z0 .GT. ZGRID(NZ)) THEN
C          IF (MURAY .GE. 0.0) THEN
C            GOTO 900
C          ENDIF
C          R = (ZGRID(NZ) - Z0)/MURAY
C          X0 = X0 + R*SQRT(1-MURAY**2)*COS(PHIRAY)
C          Y0 = Y0 + R*SQRT(1-MURAY**2)*SIN(PHIRAY)
C          Z0 = ZGRID(NZ)
C        ELSE IF (Z0 .LT. ZGRID(1)) THEN
C          IERR = 1
C          WRITE (ERRMSG,*) 'MAKE_ADJOINT: Level below domain'
C          RETURN
C        ENDIF

C       Trace a ray through the domain and update the temporary
C       adjoint source array (ADJTEMP).
C        TRANSMIT = 1.0D0
C        CALL PENCIL_BEAM_PROP (X0, Y0, Z0, BCFLAG, IPFLAG,
C     .    MAGNITUDE, MU, PHI, ADJTEMP, TOTAL_EXT, NX, NY,
C     .    NZ, NCELLS, NPTS, CELLFLAGS, XGRID, YGRID, ZGRID,
C     .    GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,
C     .    SIDE, XE,YE,ZE,TRANSMIT,TAUTOL, IERR, ERRMSG)
C        IF (IERR .NE. 0) RETURN
C
C900     CONTINUE
C      ENDDO

C     Add the last temporary array.
C      NSOURCETEMP = 1
C      DO I=1,NPTS
C        IF (ABS(ADJTEMP(I)) .GE. 1e-6) THEN
C             Update the coordinate list format sparse representation
C             of the adjoint source.
C          ADJSOURCE(NSOURCE + NSOURCETEMP) = ADJTEMP(I)
C          ANGLEPTRS(NSOURCE + NSOURCETEMP) = ANGLEPTR
C          PTSPTRS(NSOURCE + NSOURCETEMP) = I
C          NSOURCETEMP = NSOURCETEMP + 1
C        ENDIF
C      ENDDO
C      NSOURCE = NSOURCE + NSOURCETEMP
C      ANGLEPTR = ANGLEPTR + 1
C      ADJTEMP = 0.0

C     Sort the ADJSOURCE and ANGLEPTRS, PTSPTRS inplace so that
C     we are in a compressed sparse column format so that
C     values are grouped by point for efficient use in
C     the CALC_SOURCE_PNT/CALC_SOURCE_PNT_UNPOL.
C     ie PTSPTRS are sorted carrying along values in ANGLEPTRS
C     and ADJSOURCE. '2' means sort in increasing order carrying
C     other variables along.
C      CALL QUICKSORT_NEW (PTSPTRS,ADJSOURCE,ANGLEPTRS, NSOURCE, 2)
C
C     Construct the ADJSOURCEPTR from the sorted PTSPTRS
C     which holds the starting point of each gridpoint's
C     adjoint source, just like RSHPTR or SHPTR does for the
C     adaptive radiance/source spherical harmonic expansions.
C      CALL CONSTRUCT_PTR(ADJSOURCEPTR, PTSPTRS, NPTS,
C     .                   8*MAXCELLSCROSS*NVIS)
C
C      RETURN
C      END
C
C      SUBROUTINE CONSTRUCT_PTR(ADJSOURCEPTR, PTSPTRS,
C     .                         NPTS, SOURCESIZE)
C
C      IMPLICIT NONE
C      INTEGER SOURCESIZE, NPTS
C      INTEGER ADJSOURCEPTR(NPTS+1)
C      INTEGER PTSPTRS(SOURCESIZE)
Cf2py intent(in) NPTS
Cf2py intent(out) ADJSOURCEPTR
Cf2py intent(in) PTSPTRS
C
C      INTEGER COUNTER, I, PTR
C      ADJSOURCEPTR(1) = 0
C      COUNTER = 1
C      PTR = PTSPTRS(1)
C      DO I=1,NPTS
C        DO WHILE (PTR .EQ. I)
C          COUNTER = COUNTER + 1
C          PTR = PTSPTRS(COUNTER)
C        ENDDO
C        ADJSOURCEPTR(I+1) = COUNTER - 1
C      ENDDO
C
C      RETURN
C      END


C      SUBROUTINE SENSOR_CHARACTERISTIC(X,Y,Z,MU,PHI,
C     .    XS, YS, ZS, MAGNITUDE, INSUPPORT, CHARACTERISTIC_TYPE,
C     .    CHARACTERISTIC_PARAMS, NCHARPARAMS)
C
C      IMPLICIT NONE
C      INTEGER NCHARPARAMS
C      CHARACTER CHARACTERISTIC_TYPE*1
C      DOUBLE PRECISION CHARACTERISTIC_PARAMS(NCHARPARMS)
C      LOGICAL INSUPPORT
C      DOUBLE PRECISION MAGNITUDE
C      DOUBLE PRECISION X,Y,Z,MU,PHI,XS,YS,ZS
C
C      DOUBLE PRECISION CX,CY,CZ
C
C      INSUPPORT = .FALSE.
C      MAGNITUDE = 0.0
C
C test for support
C find closest approach of central ray to
C the reference position X,Y,Z (gridpoint position)
C      CX = SQRT(1.0D0-MU**2)*COS(PHI)
C      CY = SQRT(1.0D0-MU**2)*SIN(PHI)
C      CZ = MU
C      DX = X-XS
C      DY = Y-YS
C      DZ = Z-ZS
C      INNER = DX*CX + DY*CY + DZ*CZ
C      DISTANCE = SQRT((DX - CX*INNER)**2 +
C     .                (DY - CY*INNER)**2 +
C     .                (DZ - CZ*INNER)**2)
C
C      IF (DISTANCE .LE. CHARACTERISTIC_PARAMS(2)) THEN
C       INSUPPORT = .TRUE.
C      ELSE
C        RETURN
C      ENDIF
C
C     A box function characteristic.
C      IF (CHARACTERISTIC_TYPE .EQ. 'B') THEN
C        MAGNITUDE = CHARACTERISTIC_PARAMS(1)
C     Gaussian characteristic
C      ELSEIF (CHARACTERISTIC_TYPE .EQ. 'G') THEN
C        MAGNITUDE = CHARACTERISTIC_PARAMS(1)*
C     .    EXP(-(DISTANCE**2)/
C     .    (2*CHARACTERISTIC_PARAMS(3)**2))/
C     .    (SQRT(2*PI)*CHARACTERISTIC_PARAMS(3))
C     Add new types here.
C      ENDIF
C     MAGNITUDE SO FAR IS CALCULATED ALONG THE BEAM,
C     However, the SHDOM solar convention is that this is
C     the projection on a horizontal plane so we have to
C     multiply by cosine of zenith.
C      MAGNITUDE = MAGNITUDE*ABS(MU)
C
C      RETURN
C      END
C
C
C
C
C      SUBROUTINE MAKE_ADJOINT_SOURCE()
C
C        IMPLICIT NONE
C
C
C        DOUBLE PRECISION X,Y,Z
C       DEFINE EVERYTHING FOR DIRECT_BEAM_PROP HERE.
C      Call the INIT for DIRECT_BEAM_PROP HERE.
C        CALL DIRECT_BEAM_PROP (1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,
C     .            DELTAM, ML, NSTLEG, NLEGP,
C     .            SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),
C     .            UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM,
C     .            NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
C     .            XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
C     .            ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
C     .            ZCKD, GASABS, CX, CY, CZ, CXINV, CYINV,
C     .            CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
C     .            XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV,
C     .	          NPART, MAXPG, PHASEWTP, MAXNMICRO)
C
C        DIRFLUX = 0.0
C        DO IP = 1, NPTS
C          X = GRIDPOS(1,IP)
C          Y = GRIDPOS(2,IP)
C          Z = GRIDPOS(3,IP)
C          CALL SENSOR_CHARACTERISTIC(X,Y,Z, MU,PHI,
C     .      XS,YS,ZS, MAGNITUDE, INSUPPORT, CHARACTERISTIC_TYPE,
C     .      CHARACTERISTIC_PARAMS, NCHARPARAMS)
C
C          IF (INSUPPORT) THEN
C            DIRPATH = 0.0
C         We have replaced the solar angles with the
C         equivalent for the sensor.
C            CALL DIRECT_BEAM_PROP (0, X,Y,Z, BCFLAG, IPFLAG, DELTAM,
C     .        ML, NSTLEG, NLEGP, MAGNITUDE,MU,PHI,
C     .        DIRFLUX(IP),
C     .        UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM,
C     .        NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
C     .        XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
C     .        ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
C     .        ZCKD, GASABS, CX, CY, CZ, CXINV, CYINV,
C     .        CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
C     .        XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV,
C     .         NPART, MAXPG, PHASEWTP, MAXNMICRO)
C          ENDIF
C        ENDDO
C        RETURN
C      END
C
C
C      SUBROUTINE COMPUTE_RADIANCE_GRID(NSTOKES, RADIANCE, RSHPTR,
C     .           NPTS, MU, PHI, ML, MM, NSTLEG, NLM, GRIDRAD)
C
C      IMPLICIT NONE
C      INTEGER NSTOKES, NPTS, ML, MM, NSTLEG, NLM
Cf2py intent(in) :: NSTOKES, NPTS, ML, MM, NSTLEG, NLM
C      INTEGER RSHPTR(*)
Cf2py intent(in) :: RSHPTR
C      REAL MU, PHI, RADIANCE(NSTOKES, *)
Cf2py intent(in) :: MU, PHI, RADIANCE
C      REAL GRIDRAD(NSTOKES,NPTS)
Cf2py intent(out) :: GRIDRAD
C
C      INTEGER I, RS, RE,J
C      REAL YLMDIR(NSTLEG,NLM)
C      GRIDRAD = 0.0
C
C      CALL YLMALL (.FALSE.,MU,PHI,ML,MM,NSTLEG, YLMDIR)
C      DO I=1,NPTS
C        RS = RSHPTR(I)
C        RE = RSHPTR(I+1)-RS
C        DO J=1,RE
C          GRIDRAD(1,I) = GRIDRAD(1,I) +
C     .      RADIANCE(1,RS+J)*YLMDIR(1,J)
C        ENDDO
C      ENDDO
C      RETURN
C      END

      SUBROUTINE COMPUTE_DIR_SOURCE(NSTOKES, SOURCE, SHPTR,
     .    TOTAL_EXT, NPTS, MU, PHI, ML, MM, NSTLEG, NLM,
     .    SOURCE1, PHASETAB, NUMPHASE, EXTINCT, ALBEDO,
     .    LEGEN, IPHASE, DIRFLUX, NPART, INTERPMETHOD,
     .    PHASEINTERPWT, PHASEMAX, DELTAM, NLEG, SRCTYPE,
     .    SOLARMU, SOLARAZ,NSCATANGLE, NSTPHASE,
     .    IERR, ERRMSG)

      IMPLICIT NONE
      INTEGER NSTOKES, NSTLEG, NPTS, ML, MM, NLM, NLEG
Cf2py intent(in) :: NSTOKES, NSTLEG, NPTS, ML, MM, NLM, NLEG
      INTEGER NUMPHASE
Cf2py intent(in) :: NUMPHASE
      INTEGER NSCATANGLE, NSTPHASE
      INTEGER IPHASE(8,NPTS,NPART), NPART
Cf2py intent(in) :: IPHASE, NPART, NSCATANGLE, NSTPHASE
      REAL    PHASEINTERPWT(8,NPTS,NPART), PHASEMAX
Cf2py intent(in) :: PHASEINTERPWT, PHASEMAX
      CHARACTER INTERPMETHOD*2
Cf2py intent(in) :: INTERPMETHOD
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
Cf2py intent(in) :: EXTINCT, ALBEDO
      REAL    LEGEN(NSTLEG,0:NLEG,NUMPHASE), TOTAL_EXT(NPTS)
Cf2py intent(in) :: LEGEN, TOTAL_EXT
      REAL    DIRFLUX(NPTS), SOURCE(NSTOKES,*)
Cf2py intent(in) :: DIRFLUX, SOURCE
      REAL    SOURCE1(NSTOKES, NPTS)
Cf2py intent(out) :: SOURCE1
      REAL    PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
Cf2py intent(in) :: PHASETAB
      CHARACTER SRCTYPE*1
Cf2py intent(in) :: SRCTYPE
      REAL    MU, PHI, SOLARMU, SOLARAZ
Cf2py intent(in) :: MU, PHI, SOLARMU, SOLARAZ
      INTEGER IERR, SHPTR(*)
Cf2py intent(out) :: IERR
Cf2py intent(in) :: SHPTR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: ERRMSG

      REAL    SINGSCAT(NSTOKES,NUMPHASE)
      DOUBLE PRECISION SECMU0, COSSCAT, F
      INTEGER I, IS, NS, J, IPA, Q, MS, ME, L
      INTEGER JT, M
      REAL EXT, LEGENT(NSTLEG,0:NLEG), DA, A1,B1
      DOUBLE PRECISION W, PI
      REAL YLMDIR(NSTLEG,NLM), YLMSUN(NSTLEG,NLM)

      SECMU0 = 1.0D0/ABS(SOLARMU)
      PI = ACOS(-1.0D0)

      IERR = 0
      IF (ABS(MU) .GT. 1.0) THEN
        IERR = 1
        WRITE(ERRMSG, *) 'COMPUTE_DIR_SOURCE: ABS(MU)>1'
      ENDIF

C         Precompute Ylm's for output direction and solar direction
      CALL YLMALL (.FALSE., MU, PHI, ML, MM, NSTLEG, YLMDIR)


      IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN

        CALL YLMALL (.TRUE., SOLARMU, SOLARAZ, ML, MM, NSTLEG,
     .                YLMSUN)
        COSSCAT = SOLARMU*MU + SQRT((1.0D0-SOLARMU**2)*(1.0D0-MU**2))
     .                  *COS(SOLARAZ-PHI)
        COSSCAT = MAX(MIN(1.0D0, COSSCAT), -1.0D0)

        F = (NSCATANGLE-1)*(ACOS(COSSCAT)/PI) + 1
        J = MIN(NSCATANGLE-1,INT(F))
        F = F - J
        DO I = 1, NUMPHASE
          SINGSCAT(1:NSTPHASE,I)
     .               = (1-F)*PHASETAB(:,I,J) + F*PHASETAB(:,I,J+1)
          IF (NSTOKES .GT. 1) THEN
            CALL ROTATE_POL_PLANE (NSTOKES, COSSCAT, SOLARMU,
     .                   MU, SOLARAZ-PHI, SINGSCAT(:,I))
          ENDIF
        ENDDO
      ENDIF

C         Loop over all the grid points, computing the source function
C           at the viewing angle from the spherical harmonic source function.
      DO I = 1, NPTS
        EXT = TOTAL_EXT(I)
        IS = SHPTR(I)
        NS = SHPTR(I+1)-IS
C         Sum over the real generalized spherical harmonic series
C         of the source function
        SOURCE1(:,I) = 0.0
        DO J = 1, NS
          SOURCE1(1,I) = SOURCE1(1,I) + SOURCE(1,IS+J)*YLMDIR(1,J)
        ENDDO
        IF (NSTOKES .GT. 1) THEN
          DO J = 1, NS
            SOURCE1(2,I) = SOURCE1(2,I) + SOURCE(2,IS+J)*YLMDIR(2,J)
     .                                  + SOURCE(3,IS+J)*YLMDIR(5,J)
            SOURCE1(3,I) = SOURCE1(3,I) + SOURCE(2,IS+J)*YLMDIR(6,J)
     .                                  + SOURCE(3,IS+J)*YLMDIR(3,J)
          ENDDO
        ENDIF
        IF (NSTOKES .EQ. 4) THEN
          DO J = 1, NS
            SOURCE1(4,I) = SOURCE1(4,I) + SOURCE(4,IS+J)*YLMDIR(4,J)
          ENDDO
        ENDIF

C             Special case for solar source and Delta-M
        IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
  	    DO IPA = 1, NPART

  	      IF (EXT.EQ.0.0) THEN
  		        W = 1.0
  	      ELSE
  		        W = EXTINCT(I,IPA)/EXT
  	      ENDIF
  	      IF (W.EQ.0.0) CYCLE

          IF (INTERPMETHOD(2:2) .EQ. 'O' ) THEN
            LEGENT = LEGEN(:,:,IPHASE(1,I,IPA))
          ELSEIF (INTERPMETHOD(2:2) .EQ. 'N') THEN
            IF (PHASEINTERPWT(1,I,IPA) .GE. PHASEMAX) THEN
              LEGENT = LEGEN(:,:,IPHASE(1,I,IPA))
            ELSE
              LEGENT = 0.0
              DO Q=1,8
                IF (PHASEINTERPWT(Q,I,IPA) .LE. 1e-5) CYCLE
                LEGENT = LEGENT + LEGEN(:,:,IPHASE(Q,I,IPA))*
     .          PHASEINTERPWT(Q,I,IPA)
              ENDDO
            ENDIF
            F = LEGENT(1,ML+1)
            LEGENT = LEGENT/(1-F)
          ENDIF

C               First subtract off the truncated single scattering
  	      DA = ALBEDO(I,IPA)*DIRFLUX(I)*SECMU0*W
  	      J = 1

  	      DO L = 0, ML
  	        ME = MIN(L,MM)
  	        MS = -ME
  	        A1 = DA*LEGENT(1,L)
  	        B1 = DA*LEGENT(5,L)
  	        IF (J .LE. NS) THEN
  	          JT = J
  	          DO M = MS, ME
  	            SOURCE1(1,I) =SOURCE1(1,I)-A1*YLMDIR(1,J)*YLMSUN(1,J)
  	            J = J + 1
              ENDDO
  	          IF (NSTOKES .GT. 1) THEN
  	            J = JT
                DO M = MS, ME
  	              SOURCE1(2,I)=SOURCE1(2,I)-B1*YLMDIR(2,J)*YLMSUN(1,J)
  	              SOURCE1(3,I)=SOURCE1(3,I)-B1*YLMDIR(6,J)*YLMSUN(1,J)
  	              J = J + 1
               ENDDO
  	         ENDIF
  	       ENDIF
  	      ENDDO

C               Then add in the single scattering contribution for the
C               original unscaled phase function.
  	      IF (NUMPHASE .GT. 0) THEN
            IF (PHASEINTERPWT(1,I,IPA) .GE. PHASEMAX) THEN
              SOURCE1(:,I) = SOURCE1(:,I) +
     .          DA*SINGSCAT(:,IPHASE(1,I,IPA))
            ELSE
              DO Q=1,8
                IF (PHASEINTERPWT(Q,I,IPA) .LE. 1e-5) CYCLE
                SOURCE1(:,I) = SOURCE1(:,I) +
     .          DA*SINGSCAT(:,IPHASE(Q,I,IPA))*PHASEINTERPWT(Q,I,IPA)
              ENDDO
            ENDIF
  	      ELSE
            WRITE(ERRMSG,*) 'COMPUTE_DIR_SOURCE:',
     .      ' NUMPHASE=0 is not supported.'
            IERR = 1
            RETURN
  	      ENDIF
  	    ENDDO
        ENDIF
      ENDDO

      RETURN
      END



      SUBROUTINE COMPUTE_SOURCE_GRID(NSTOKES, SOURCE, SHPTR,
     .           NPTS, MU, PHI, ML, MM, NSTLEG, NLM, GRIDSOURCE)

      IMPLICIT NONE
      INTEGER NSTOKES, NPTS, ML, MM, NSTLEG, NLM
!f2py intent(in) :: NSTOKES, NPTS, ML, MM, NSTLEG, NLM
      INTEGER SHPTR(*)
!f2py intent(in) :: SHPTR
      REAL MU, PHI, SOURCE(NSTOKES, *)
!f2py intent(in) :: MU, PHI, SOURCE
      REAL GRIDSOURCE(NSTOKES,NPTS)
!f2py intent(out) :: GRIDSOURCE

      INTEGER I, RS, RE,J
      REAL YLMDIR(NSTLEG,NLM)
      GRIDSOURCE = 0.0

      CALL YLMALL (.FALSE.,MU,PHI,ML,MM,NSTLEG, YLMDIR)
      DO I=1,NPTS
        RS = SHPTR(I)
        RE = SHPTR(I+1)-RS
        DO J=1,RE
          GRIDSOURCE(1,I) = GRIDSOURCE(1,I) +
     .      SOURCE(1,RS+J)*YLMDIR(1,J)
        ENDDO
      ENDDO
      RETURN
      END

C      SUBROUTINE INNER_PRODUCT(NPTS,NCELLS,GRIDPTR, TREEPTR,
C     .    NEIGHPTR, CELLFLAGS, BCFLAG, IPFLAG, XGRID, YGRID,
C     .    ZGRID, GRIDPOS, FOR_RADIANCE, PSEUDO_RADIANCE,
C     .    FOR_RSHPTR, PSEUDO_RSHPTR, ML, MM, NLM, NSTLEG,
C     .    NSTOKES, OUTPUT)
CC    Loop over all of the deepest layer of cells and perform
CC    the inner product for each one. This will contribute to
CC    the derivative on the property grid.
C      IMPLICIT NONE
C      RETURN
C      END

      SUBROUTINE TRAVERSE_GRID(NX, NY, NZ, NPTS, NCELLS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                 X0, Y0, Z0, XE,YE,ZE, SIDE)
C    Finds the ending point of a ray integration so that we can
C    trace a ray from the end point to the starting point.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS,SIDE
Cf2py intent(in) :: NX,NY,NZ,NPTS,NCELLS
Cf2py intent(out) :: SIDE
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in) :: CELLFLAGS
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      REAL    MURAY, PHIRAY, MU2, PHI2
Cf2py intent(in) :: MURAY, PHIRAY, MU2,PHI2
      DOUBLE PRECISION X0, Y0, Z0
Cf2py intent(in) :: X0, Y0, Z0
      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) :: BCFLAG, IPFLAG

      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, K
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8)
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      INTEGER DONETHIS(8)
      REAL    XM, YM
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
Cf2py intent(out) :: XE,YE,ZE
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, F1(8), F2(8)

      REAL FM
      INTEGER J

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/
      DATA DONEFACE/0,0,0,0,0,0,0,0, 0,1,0,3,0,5,0,7, 2,0,4,0,6,0,8,0,
     .              0,0,1,2,0,0,5,6, 3,4,0,0,5,6,0,0,
     .              0,0,0,0,1,2,3,4, 5,6,7,8,0,0,0,0/

      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)

C         Make the ray direction (opposite to the discrete ordinate direction)
      CX = SQRT(1.0-MURAY**2)*COS(PHIRAY)
      CY = SQRT(1.0-MURAY**2)*SIN(PHIRAY)
      CZ = MURAY
      IF (ABS(CX) .GT. 1.0E-6) THEN
        CXINV = 1.0D0/CX
      ELSE
        CX = 0.0
        CXINV = 1.0E6
      ENDIF
      IF (ABS(CY) .GT. 1.0E-6) THEN
        CYINV = 1.0D0/CY
      ELSE
        CY = 0.0
        CYINV = 1.0E6
      ENDIF
      IF (ABS(CZ) .GT. 1.0E-6) THEN
        CZINV = 1.0D0/CZ
      ELSE
        CZ = 0.0
        CZINV = 1.0E6
      ENDIF

C         Setup for the ray path direction
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
      IF (CZ .LT. 0.0) THEN
        BITZ = 1
      ELSE
        BITZ = 0
      ENDIF
      IOCT = 1 + BITX + 2*BITY + 4*BITZ
      XM = 0.5*(XGRID(1)+XGRID(NX))
      YM = 0.5*(YGRID(1)+YGRID(NY))

C         Start at the desired point, getting the extinction and source there
      XE = X0
      YE = Y0
      ZE = Z0
      CALL LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, XE, YE, ZE, ICELL)

      IFACE = 0
      NGRID = 0

C         Loop until reach a Z boundary or transmission is very small
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
C           Make sure current cell is valid
        IF (ICELL .LE. 0) THEN
          WRITE (6,*)'TRAVERSE_GRID: ICELL=',ICELL,
     .                MURAY,PHIRAY,XE,YE,ZE
          STOP
        ENDIF
        NGRID = NGRID + 1

C           This cell is independent pixel if IP mode or open boundary
C             conditions and ray is leaving domain (i.e. not entering)
        IPINX = BTEST(INT(CELLFLAGS(ICELL)),0) .AND.
     .          .NOT. ( BTEST(BCFLAG,0) .AND.
     .          ((CX.GT.0.AND.XE.LT.XM) .OR. (CX.LT.0.AND.XE.GT.XM)) )
        IPINY = BTEST(INT(CELLFLAGS(ICELL)),1) .AND.
     .          .NOT. ( BTEST(BCFLAG,1) .AND.
     .          ((CY.GT.0.AND.YE.LT.YM) .OR. (CY.LT.0.AND.YE.GT.YM)) )

C           Find boundaries of the current cell
C           Find the three possible intersection planes (X,Y,Z)
C             from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9-IOCT,ICELL)
C           Get the distances to the 3 planes and select the closest
C             (always need to deal with the cell that is wrapped)
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
          WRITE (6,*) 'TRAVERSE_GRID: SO<0  ',
     .      MURAY,PHIRAY,XE,YE,ZE,SO,ICELL
          STOP
        ENDIF
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ

C               Get the intersection face number (i.e. neighptr index)
        IF (SOX .LE. SOZ .AND. SOX .LE. SOY) THEN
          IFACE = 2-BITX
          JFACE = 1
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),0).AND.BTEST(BCFLAG,0)
        ELSE IF (SOY .LE. SOZ) THEN
          IFACE = 4-BITY
          JFACE = 2
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),1).AND.BTEST(BCFLAG,1)
        ELSE
          IFACE = 6-BITZ
          JFACE = 3
          OPENBCFACE=.FALSE.
        ENDIF
C            Get the next cell to go to
        INEXTCELL = NEIGHPTR(IFACE,ICELL)
        IF (INEXTCELL .LT. 0) THEN
          CALL NEXT_CELL (XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
        ENDIF
C             If going to same or larger face then use previous face
        IF (NEIGHPTR(IFACE,ICELL) .GE. 0 .AND. .NOT.OPENBCFACE) THEN
          KFACE = IFACE
          IC = ICELL
        ELSE
C             If going to smaller face then use next face (more accurate)
          KFACE = OPPFACE(IFACE)
          IC = INEXTCELL
          IFACE = 0
        ENDIF
C           Get the location coordinate
        IF (INEXTCELL .GT. 0) THEN
          IF (JFACE .EQ. 1) THEN
            XN = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
          ELSE IF (JFACE .EQ. 2) THEN
            YN = GRIDPOS(2,GRIDPTR(IOCT,INEXTCELL))
          ELSE
            ZN = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
          ENDIF
        ENDIF

C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
        IF ((INEXTCELL .EQ. 0).OR.(NGRID.GT.MAXCELLSCROSS)) THEN
          DONE = .TRUE.
        ELSE
          XE = XN
          YE = YN
          ZE = ZN
          ICELL = INEXTCELL
        ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE PHASEFUNC_FROM_WIGNER(NANGLES, NLEG, LEGCOEF,
     .                                 PHASE, NSTPHASE, IERR,
     .    ERRMSG, NEGCHECK, NSTLEG)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: NEGCHECK
      INTEGER, INTENT(IN) :: NANGLES, NLEG, NSTLEG, NSTPHASE
      REAL,    INTENT(IN) :: LEGCOEF(NSTLEG,0:NLEG)
      REAL,    INTENT(OUT) :: PHASE(NSTPHASE, NANGLES)
      CHARACTER ERRMSG*600
      INTEGER IERR

      INTEGER J, L
      DOUBLE PRECISION  PI, OFOURPI, COSSCAT, FCT, F, X, A1, B1
      DOUBLE PRECISION, ALLOCATABLE :: DMM1(:), DMM2(:)
      ALLOCATE (DMM1(0:NLEG), DMM2(0:NLEG))

      PI = ACOS(-1.0D0)
      OFOURPI = 1.0/(4.0*PI)
      IERR = 0



      DO J = 1, NANGLES
        COSSCAT = COS(PI*DFLOAT(J-1)/(NANGLES-1))
        X = DBLE(COSSCAT)
        CALL WIGNERFCT (X, NLEG, 0, 0, DMM1)
        CALL WIGNERFCT (X, NLEG, 2, 0, DMM2)

C           Sum the first Wigner function series, which is actually a Legendre series
        A1 = 0.0D0
        DO L = 0, NLEG
          FCT = 1.0D0
!2.0D0*L + 1.0D0
          A1  = A1 + FCT*DBLE(LEGCOEF(1,L))*DMM1(L)
        ENDDO
        IF (NEGCHECK .AND. A1 .LE. 0.0) THEN
          IERR=1
          WRITE (ERRMSG,*) 'PHASEFUNC_FROM_WIGNER: ',
     .          'negative phase ',
     .         'function'
          RETURN
        ENDIF
        PHASE(1,J) = SNGL(A1*OFOURPI)

        B1 = 0.0D0
        DO L = 0, NLEG
          FCT = 1.0D0
!2.0D0*L + 1.0D0
          B1  = B1 - FCT*DBLE(LEGCOEF(5,L))*DMM2(L)
        ENDDO
        PHASE(2,J) = SNGL(B1*OFOURPI)

      ENDDO
      DEALLOCATE (DMM1, DMM2)
      RETURN
      END

      SUBROUTINE PHASEFUNC_FROM_LEGENDRE (NANGLES, NCOEF, LEGCOEF,
     .                                    PHASE)
       ! Calculates the phase function at NANGLES equally spaced angles
       ! from the Legendres series expansion in LEGCOEF.  The scalar (I,I)
       ! element of the phase matrix is expanded in Wigner d^l_{0,0} functions,
       ! which are the same as Legendre polynomials in \cos\Theta.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NANGLES, NCOEF
        REAL,    INTENT(IN) :: LEGCOEF(0:NCOEF)
        REAL,    INTENT(OUT) :: PHASE(NANGLES)
        INTEGER :: J, L
        REAL(8) :: RD, MU, SUM, PL, PL1, PL2

        RD = ACOS(-1.0D0)/180.D0
        DO J = 1, NANGLES
          MU = COS(RD*J*180.D0/NANGLES)
          SUM = 0.0
            ! Use upward recurrence to find Legendre polynomials
          PL1 = 1.0
          PL = 1.0
          DO L = 0, NCOEF
            IF (L .GT. 0) PL = (2*L-1)*MU*PL1/L-(L-1)*PL2/L
            SUM = SUM + LEGCOEF(L)*PL
            PL2 = PL1
            PL1 = PL
          ENDDO
          PHASE(J) = SUM
        ENDDO
      END SUBROUTINE PHASEFUNC_FROM_LEGENDRE


      SUBROUTINE TRANSMISSION_INTEGRAL(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,TOTAL_EXT,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX, PATH, TRANSCUT,
     .             TAUTOL, FIELD)
C    Integrates a field weighted by a EXTINCT*TRANSMISSION kernel.
C    This is equivalent to the radiance integration if the field
C    were the SOURCE evaluated in the appropriate direction.
C    This subroutine is used to propagate uncertainty
C    in the SOURCE to uncertainty in the radiances.
C    - JRLoveridge 2021/09/14
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
Cf2py intent(in) :: NX, NY, NZ, NPTS, NCELLS
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) :: BCFLAG, IPFLAG
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in) :: CELLFLAGS
      DOUBLE PRECISION CAMX(*), CAMY(*), CAMZ(*), CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX
      REAL  PATH(NPIX)
Cf2py intent(out):: PATH
      REAL TOTAL_EXT(NPTS), FIELD(NPTS)
Cf2py intent(in) :: TOTAL_EXT, FIELD
      DOUBLE PRECISION TRANSCUT, TAUTOL
Cf2py intent(in) :: TRANSCUT, TAUTOL

      REAL TAU
      INTEGER N, K
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, R, PI

      PI = ACOS(-1.0D0)
      PATH = 0.0
C         Loop over pixels in image
      DO N = 1, NPIX
        X0 = CAMX(N)
        Y0 = CAMY(N)
        Z0 = CAMZ(N)
        MU2 = CAMMU(N)
        PHI2 = CAMPHI(N)
        MURAY = -MU2
        PHIRAY = PHI2 - PI

C             Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (MURAY .GE. 0.0) THEN
            GOTO 900
          ENDIF
          R = (ZGRID(NZ) - Z0)/MURAY
          X0 = X0 + R*SQRT(1-MURAY**2)*COS(PHIRAY)
          Y0 = Y0 + R*SQRT(1-MURAY**2)*SIN(PHIRAY)
          Z0 = ZGRID(NZ)
        ELSE IF (Z0 .LT. ZGRID(1)) THEN
          WRITE (6,*) 'OPTICAL_DEPTH: Level', Z0, 'below domain',
     .		                      ZGRID(1)
          STOP
        ENDIF
        TAU = 0.0
        CALL TRANSMISSION_INTEGRAL_1RAY(NX, NY, NZ, NPTS, NCELLS,
     .                          GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                          BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                          GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                X0, Y0, Z0, TAU, TAUTOL, TRANSCUT,
     .                      TOTAL_EXT, FIELD)
        PATH(N) = TAU
900     CONTINUE

      ENDDO

      RETURN
      END

      SUBROUTINE TRANSMISSION_INTEGRAL_1RAY(NX, NY, NZ, NPTS, NCELLS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                 X0, Y0, Z0, INTEGRAL, TAUTOL, TRANSCUT,
     .                       TOTAL_EXT, FIELD)

      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0
      REAL TOTAL_EXT(NPTS), FIELD(NPTS)
      REAL    INTEGRAL, EXT0, EXT1, EXTN, FIELDEXT0, FIELDEXT1
C      REAL    PATH(*)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, K, SIDE
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), BCFLAG, IPFLAG
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    XM, YM
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, F(8)
      DOUBLE PRECISION TAUTOL, TAUGRID, S, DELS, TRANSCUT
      DOUBLE PRECISION EXT, TAU, TRANSCELL, ABSCELL, SRC
      DOUBLE PRECISION TRANSMIT

      INTEGER J

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/
      TRANSMIT = 1.0D0

      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)
C      TAU = 0.0

C         Make the ray direction (opposite to the discrete ordinate direction)
      CX = SQRT(1.0-MURAY**2)*COS(PHIRAY)
      CY = SQRT(1.0-MURAY**2)*SIN(PHIRAY)
      CZ = MURAY
      IF (ABS(CX) .GT. 1.0E-6) THEN
        CXINV = 1.0D0/CX
      ELSE
        CX = 0.0
        CXINV = 1.0E6
      ENDIF
      IF (ABS(CY) .GT. 1.0E-6) THEN
        CYINV = 1.0D0/CY
      ELSE
        CY = 0.0
        CYINV = 1.0E6
      ENDIF
      IF (ABS(CZ) .GT. 1.0E-6) THEN
        CZINV = 1.0D0/CZ
      ELSE
        CZ = 0.0
        CZINV = 1.0E6
      ENDIF

C         Setup for the ray path direction
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
      IF (CZ .LT. 0.0) THEN
        BITZ = 1
      ELSE
        BITZ = 0
      ENDIF
      IOCT = 1 + BITX + 2*BITY + 4*BITZ
      XM = 0.5*(XGRID(1)+XGRID(NX))
      YM = 0.5*(YGRID(1)+YGRID(NY))

C         Start at the desired point, getting the extinction and source there
      XE = X0
      YE = Y0
      ZE = Z0
      CALL LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, XE, YE, ZE, ICELL)

      IFACE = 0
      NGRID = 0

C         Loop until reach a Z boundary or transmission is very small
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
C           Make sure current cell is valid
        IF (ICELL .LE. 0) THEN
          WRITE (6,*)'INTEGRATE_1RAY: ICELL=',ICELL,
     .                MURAY,PHIRAY,XE,YE,ZE
          STOP
        ENDIF
        NGRID = NGRID + 1

C         Interpolate the source and extinction to the current point
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XE, YE, ZE, F)
        EXT1 = 0.0
        FIELDEXT1 = 0.0
        DO J=1,8
          EXT1 = EXT1 + TOTAL_EXT(J)*F(J)
          FIELDEXT1 = FIELDEXT1 + TOTAL_EXT(J)*F(J)*FIELD(J)
        ENDDO

C           This cell is independent pixel if IP mode or open boundary
C             conditions and ray is leaving domain (i.e. not entering)
        IPINX = BTEST(INT(CELLFLAGS(ICELL)),0) .AND.
     .          .NOT. ( BTEST(BCFLAG,0) .AND.
     .          ((CX.GT.0.AND.XE.LT.XM) .OR. (CX.LT.0.AND.XE.GT.XM)) )
        IPINY = BTEST(INT(CELLFLAGS(ICELL)),1) .AND.
     .          .NOT. ( BTEST(BCFLAG,1) .AND.
     .          ((CY.GT.0.AND.YE.LT.YM) .OR. (CY.LT.0.AND.YE.GT.YM)) )

C           Find boundaries of the current cell
C           Find the three possible intersection planes (X,Y,Z)
C             from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9-IOCT,ICELL)
C           Get the distances to the 3 planes and select the closest
C             (always need to deal with the cell that is wrapped)
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
          WRITE (6,*) 'INTEGRATE_1RAY: SO<0  ',
     .      MURAY,PHIRAY,XE,YE,ZE,SO,ICELL
          STOP
        ENDIF
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ
C           Find the optical path across the grid cell and figure how
C             many subgrid intervals to use
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XN, YN, ZN, F)
        EXTN = 0.0
        DO J=1,8
          EXTN = EXTN + F(J)*TOTAL_EXT(J)
        ENDDO

        TAUGRID = SO*0.5*(EXT1+EXTN)
        NTAU = MAX(1,1+INT(TAUGRID/TAUTOL))
        DELS = SO/NTAU

C           Loop over the subgrid cells
        DO IT = 1, NTAU
          S = IT*DELS
          XI = XE + S*CX
          YI = YE + S*CY
          ZI = ZE + S*CZ

          EXT0 = 0.0
          FIELDEXT0 = 0.0
          DO J=1,8
            EXT0 = EXT0 + F(J)*TOTAL_EXT(J)
            FIELDEXT0 = FIELDEXT0 + F(J)*TOTAL_EXT(J)*FIELD(J)
          ENDDO

C            Compute the subgrid radiance: integration of the source function
          EXT = 0.5*(EXT0+EXT1)
          IF (EXT .NE. 0.0) THEN
            TAU=EXT*DELS
            ABSCELL = TAU*(1.0-0.5*TAU*(1.0-0.33333333333*TAU))
            TRANSCELL = 1.0 - ABSCELL
C                 Linear extinction, linear source*extinction, to second order
            SRC = ( 0.5*(FIELDEXT0+FIELDEXT1)
     .           + 0.08333333333*(EXT0*FIELDEXT1-
     .                EXT1*FIELDEXT0)*DELS
     .                *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
          ENDIF

          INTEGRAL = INTEGRAL + TRANSMIT*SRC*ABSCELL
          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          FIELDEXT1 = FIELDEXT0
C                End of sub grid cell loop
        ENDDO

        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))


C               Get the intersection face number (i.e. neighptr index)
        IF (SOX .LE. SOZ .AND. SOX .LE. SOY) THEN
          IFACE = 2-BITX
          JFACE = 1
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),0).AND.BTEST(BCFLAG,0)
        ELSE IF (SOY .LE. SOZ) THEN
          IFACE = 4-BITY
          JFACE = 2
          OPENBCFACE=BTEST(INT(CELLFLAGS(ICELL)),1).AND.BTEST(BCFLAG,1)
        ELSE
          IFACE = 6-BITZ
          JFACE = 3
          OPENBCFACE=.FALSE.
        ENDIF
C            Get the next cell to go to
        INEXTCELL = NEIGHPTR(IFACE,ICELL)
        IF (INEXTCELL .LT. 0) THEN
          CALL NEXT_CELL (XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,  INEXTCELL)
        ENDIF
C             If going to same or larger face then use previous face
        IF (NEIGHPTR(IFACE,ICELL) .GE. 0 .AND. .NOT.OPENBCFACE) THEN
          KFACE = IFACE
          IC = ICELL
        ELSE
C             If going to smaller face then use next face (more accurate)
          KFACE = OPPFACE(IFACE)
          IC = INEXTCELL
          IFACE = 0
        ENDIF
C           Get the location coordinate
        IF (INEXTCELL .GT. 0) THEN
          IF (JFACE .EQ. 1) THEN
            XN = GRIDPOS(1,GRIDPTR(IOCT,INEXTCELL))
          ELSE IF (JFACE .EQ. 2) THEN
            YN = GRIDPOS(2,GRIDPTR(IOCT,INEXTCELL))
          ELSE
            ZN = GRIDPOS(3,GRIDPTR(IOCT,INEXTCELL))
          ENDIF
        ENDIF

C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
        IF ((INEXTCELL .EQ. 0).OR.(TRANSMIT .LT. TRANSCUT)) THEN
          DONE = .TRUE.
        ELSE
          XE = XN
          YE = YN
          ZE = ZN
          ICELL = INEXTCELL
        ENDIF
      ENDDO
      SIDE = IFACE

      RETURN
      END
