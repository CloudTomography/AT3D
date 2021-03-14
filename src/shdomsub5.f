C     Joint routines for both Polarized and Unpolarized sources
C     This contains adaptations of subroutines written by Frank Evans
C     for use in SHDOM.
C     https://nit.coloradolinux.com/~evans/shdom.html
C     See shdom.txt for documentation.

C     These subroutines have been written for use in pyshdom by
C     Aviad Levis, Technion Institute of Technology, 2019 and
C     later modified by Jesse Loveridge, University of Illinois at Urbana-Champaign,
C     2020-2021.

      SUBROUTINE MIN_OPTICAL_DEPTH(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS, EXTINCT,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX, PATH, IPHASE,
     .             ALBEDO, LEGEN, DELTAMPATH, DELTAM,
     .             NPART, NSTLEG, NLEG,ML,PATHS_SIZE)
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
      REAL    CAMX(*), CAMY(*), CAMZ(*), CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, NLEG, NPART, NSTLEG
Cf2py intent(in) :: NPIX,NLEG,NPART,NSTLEG
      REAL  PATH(NPTS,PATHS_SIZE)
Cf2py intent(out):: PATH
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
Cf2py intent(in) :: EXTINCT, ALBEDO
      REAL  LEGEN(NSTLEG,0:NLEG,*)
Cf2py intent(in) :: LEGEN
      INTEGER IPHASE(NPTS,NPART), ML
Cf2py intent(in) :: IPHASE, ML
      LOGICAL DELTAMPATH, DELTAM
Cf2py intent(in) :: DELTAMPATH, DELTAM
      INTEGER N, K
      REAL    MURAY, PHIRAY, MU2, PHI2
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
     .                      NPART, NSTLEG, NLEG,ML,N, PATHS_SIZE)
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
     .                      NPART, NSTLEG, NLEG,ML,N, PATHS_SIZE)

C     Integrates the extinction field (EXTINCT) along a ray and stores the minimum
C     optical path that intersects a cell adjacent to each gridpoint.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS, PATHS_SIZE
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(*), N
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    MURAY, PHIRAY, MU2, PHI2
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
      REAL  LEGEN(NSTLEG,0:NLEG,*)
      INTEGER IPHASE(NPTS,NPART), ML
      LOGICAL DELTAMPATH, DELTAM

      REAL FM, UNSCALED_ALBEDO
      INTEGER J, IPART

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
            IF (DELTAM .AND. .NOT. DELTAMPATH) THEN
              FM = LEGEN(1,ML+1,IPHASE(GRIDPTR(J,ICELL),IPART))
              UNSCALED_ALBEDO = ALBEDO(GRIDPTR(J,ICELL),IPART)/
     .                  (FM*ALBEDO(GRIDPTR(J,ICELL),IPART)+ 1.0)
              EXT1 = EXT1 + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)/
     .                   (1.0 - UNSCALED_ALBEDO*FM)
            ELSE IF (DELTAMPATH .AND. .NOT. DELTAM) THEN
              FM = LEGEN(1,ML+1,IPHASE(GRIDPTR(J,ICELL),IPART))
              EXT1 = EXT1 + F(J)*(1.0-ALBEDO(GRIDPTR(J,ICELL),IPART)
     .                  *FM)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ELSE
              EXT1 = EXT1 + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ENDIF
          ENDDO
       ENDDO

C        EXT1 = F(1)*TOTAL_EXT(GRIDPTR(1,ICELL)) +
C     .         F(2)*TOTAL_EXT(GRIDPTR(2,ICELL)) +
C     .         F(3)*TOTAL_EXT(GRIDPTR(3,ICELL)) +
C     .         F(4)*TOTAL_EXT(GRIDPTR(4,ICELL)) +
C     .         F(5)*TOTAL_EXT(GRIDPTR(5,ICELL)) +
C     .         F(6)*TOTAL_EXT(GRIDPTR(6,ICELL)) +
C     .         F(7)*TOTAL_EXT(GRIDPTR(7,ICELL)) +
C     .         F(8)*TOTAL_EXT(GRIDPTR(8,ICELL))

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
              FM = LEGEN(1,ML+1,IPHASE(GRIDPTR(J,ICELL),IPART))
              UNSCALED_ALBEDO = ALBEDO(GRIDPTR(J,ICELL),IPART)/
     .                  (FM*ALBEDO(GRIDPTR(J,ICELL),IPART)+ 1.0)
              EXTN = EXTN + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)/
     .                   (1.0 - UNSCALED_ALBEDO*FM)
            ELSE IF (DELTAMPATH .AND. .NOT. DELTAM) THEN
              FM = LEGEN(1,ML+1,IPHASE(GRIDPTR(J,ICELL),IPART))
              EXTN = EXTN + F(J)*(1.0-ALBEDO(GRIDPTR(J,ICELL),IPART)
     .                  *FM)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ELSE
              EXTN = EXTN + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ENDIF
          ENDDO
       ENDDO

C        EXTN = F(1)*TOTAL_EXT(GRIDPTR(1,ICELL)) +
C     .         F(2)*TOTAL_EXT(GRIDPTR(2,ICELL)) +
C     .         F(3)*TOTAL_EXT(GRIDPTR(3,ICELL)) +
C     .         F(4)*TOTAL_EXT(GRIDPTR(4,ICELL)) +
C     .         F(5)*TOTAL_EXT(GRIDPTR(5,ICELL)) +
C     .         F(6)*TOTAL_EXT(GRIDPTR(6,ICELL)) +
C     .         F(7)*TOTAL_EXT(GRIDPTR(7,ICELL)) +
C     .         F(8)*TOTAL_EXT(GRIDPTR(8,ICELL))
        TAU = TAU + SO*0.5*(EXT1+EXTN)
        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
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
     .             NPART, NSTLEG, NLEG, ML)
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
      REAL    CAMX(*), CAMY(*), CAMZ(*), CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, NLEG, NSTLEG, NPART
Cf2py intent(in) :: NPIX, NLEG, NSTLEG, NPART
      REAL  PATH(NPIX)
Cf2py intent(out):: PATH
      INTEGER IPHASE(NPTS,NPART), ML
Cf2py intent(in) :: IPHASE, ML
      REAL EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
Cf2py intent(in) :: EXTINCT, ALBEDO
      REAL LEGEN(NSTLEG,0:NLEG,*)
Cf2py intent(in) :: LEGEN
      LOGICAL DELTAMPATH, DELTAM
Cf2py intent(in) :: DELTAMPATH, DELTAM

      REAL TAU
      INTEGER N, K
      REAL    MURAY, PHIRAY, MU2, PHI2
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
     .                      ML)
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
     .                   ML)

C       Integrates the source function through the extinction field
C     (EXTINCT) backward in the direction (MURAY,PHIRAY) to find the
C     outgoing radiance (RAD) at the point X0,Y0,Z0.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0
      INTEGER NPART, NSTLEG, NLEG
      INTEGER IPHASE(NPTS,NPART)
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
      INTEGER J, IPART, ML

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
            IF (DELTAM .AND. .NOT. DELTAMPATH) THEN
              FM = LEGEN(1,ML+1,IPHASE(GRIDPTR(J,ICELL),IPART))
              UNSCALED_ALBEDO = ALBEDO(GRIDPTR(J,ICELL),IPART)/
     .                  (FM*ALBEDO(GRIDPTR(J,ICELL),IPART)+ 1.0)
              EXT1 = EXT1 + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)/
     .                   (1.0 - UNSCALED_ALBEDO*FM)
            ELSE IF (DELTAMPATH .AND. .NOT. DELTAM) THEN
              FM = LEGEN(1,ML+1,IPHASE(GRIDPTR(J,ICELL),IPART))
              EXT1 = EXT1 + F(J)*(1.0-ALBEDO(GRIDPTR(J,ICELL),IPART)
     .                  *FM)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ELSE
              EXT1 = EXT1 + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ENDIF
          ENDDO
       ENDDO

C        EXT1 = F(1)*TOTAL_EXT(GRIDPTR(1,ICELL)) +
C     .         F(2)*TOTAL_EXT(GRIDPTR(2,ICELL)) +
C     .         F(3)*TOTAL_EXT(GRIDPTR(3,ICELL)) +
C     .         F(4)*TOTAL_EXT(GRIDPTR(4,ICELL)) +
C     .         F(5)*TOTAL_EXT(GRIDPTR(5,ICELL)) +
C     .         F(6)*TOTAL_EXT(GRIDPTR(6,ICELL)) +
C     .         F(7)*TOTAL_EXT(GRIDPTR(7,ICELL)) +
C     .         F(8)*TOTAL_EXT(GRIDPTR(8,ICELL))

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
              FM = LEGEN(1,ML+1,IPHASE(GRIDPTR(J,ICELL),IPART))
              UNSCALED_ALBEDO = ALBEDO(GRIDPTR(J,ICELL),IPART)/
     .                  (FM*ALBEDO(GRIDPTR(J,ICELL),IPART)+ 1.0)
              EXTN = EXTN + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)/
     .                   (1.0 - UNSCALED_ALBEDO*FM)
            ELSE IF (DELTAMPATH .AND. .NOT. DELTAM) THEN
              FM = LEGEN(1,ML+1,IPHASE(GRIDPTR(J,ICELL),IPART))
              EXTN = EXTN + F(J)*(1.0-ALBEDO(GRIDPTR(J,ICELL),IPART)
     .                  *FM)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ELSE
              EXTN = EXTN + F(J)*EXTINCT(GRIDPTR(J,ICELL),IPART)
            ENDIF
          ENDDO
       ENDDO

C        EXTN = F(1)*TOTAL_EXT(GRIDPTR(1,ICELL)) +
C     .         F(2)*TOTAL_EXT(GRIDPTR(2,ICELL)) +
C     .         F(3)*TOTAL_EXT(GRIDPTR(3,ICELL)) +
C     .         F(4)*TOTAL_EXT(GRIDPTR(4,ICELL)) +
C     .         F(5)*TOTAL_EXT(GRIDPTR(5,ICELL)) +
C     .         F(6)*TOTAL_EXT(GRIDPTR(6,ICELL)) +
C     .         F(7)*TOTAL_EXT(GRIDPTR(7,ICELL)) +
C     .         F(8)*TOTAL_EXT(GRIDPTR(8,ICELL))
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


      SUBROUTINE SPACE_CARVE(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX, VOLUME, FLAGS, WEIGHTS,
     .              LINEAR)
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
      REAL    CAMX(*), CAMY(*), CAMZ(*), CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX
      INTEGER FLAGS(*)
Cf2py intent(in) :: FLAG
      REAL WEIGHTS(*)
Cf2py intent(in) :: WEIGHTS
      LOGICAL LINEAR
Cf2py intent(in) LINEAR
      REAL  VOLUME(3,NPTS)
Cf2py intent(out):: VOLUME

      INTEGER N, K, FLAG
      REAL    MURAY, PHIRAY, MU2, PHI2, WEIGHT
      DOUBLE PRECISION X0, Y0, Z0, R, PI

      PI = ACOS(-1.0D0)
      VOLUME = 0.0
C         Loop over pixels in image
      DO N = 1, NPIX
        X0 = CAMX(N)
        Y0 = CAMY(N)
        Z0 = CAMZ(N)
        MU2 = CAMMU(N)
        PHI2 = CAMPHI(N)
        FLAG = FLAGS(N)
        WEIGHT = WEIGHTS(N)
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

        CALL SPACE_CARVE_1RAY(NX, NY, NZ, NPTS, NCELLS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			             X0, Y0, Z0, VOLUME, FLAG, WEIGHT, LINEAR)
900     CONTINUE

      ENDDO

      RETURN
      END


      SUBROUTINE SPACE_CARVE_1RAY(NX, NY, NZ, NPTS, NCELLS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			             X0, Y0, Z0, VOLUME, FLAG, WEIGHT, LINEAR)

C     Performs the space carving by traversing the grid.
C     LINEAR mode is untested.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    MURAY, PHIRAY, MU2, PHI2, WEIGHT
      DOUBLE PRECISION X0, Y0, Z0
      REAL  VOLUME(3,NPTS), TEMP_VOLUME(NPTS), PATH_SUM
      LOGICAL LINEAR
      INTEGER FLAG
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, K
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), BCFLAG, IPFLAG
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    XM, YM
      DOUBLE PRECISION XMID,YMID,ZMID, F(8)
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/

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
     .    XMID, YMID, ZMID, F)

C	If this is not a boundary cell (currently assuming that the bc conditions are open)
C	    Check if the zlevel is within the zparticle levels
C		Compute the gradient field in direction  (MU2,PHI2)

        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        IF (.NOT.OUTOFDOMAIN) THEN
            DO K=1,8
              IF (FLAG .EQ. 1) THEN
                VOLUME(1,GRIDPTR(K,ICELL))=
     .              VOLUME(1,GRIDPTR(K,ICELL)) + 1
              ELSEIF (FLAG .EQ. 0) THEN
                 VOLUME(2,GRIDPTR(K,ICELL))=
     .                VOLUME(2,GRIDPTR(K,ICELL)) + 1
              ENDIF
              IF (LINEAR) THEN
                VOLUME(3,GRIDPTR(K,ICELL)) =
     .              VOLUME(3,GRIDPTR(K,ICELL)) + WEIGHT/(F(K)*SO)
              ELSE
                VOLUME(3,GRIDPTR(K,ICELL)) =
     .              VOLUME(3,GRIDPTR(K,ICELL)) + WEIGHT
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

      SUBROUTINE PROJECT(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,WEIGHTS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX,PATH, RETURN_MATRIX,
     .             MATRIX, MATRIX_PTRS, MATRIX_SIZE,
     .             PIXEL_INDICES, RAY_WEIGHTS, NRAYS)
C     Performs ray integration (PATH) of a set of WEIGHTS assuming
C     a linear interpolation kernel. Is also supposed to output
C     the matrix of describing this integration but it is inconsistent and
C     should not be used, likely due to boundary condition issues.
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
      REAL  CAMX(NRAYS), CAMY(NRAYS), CAMZ(NRAYS)
      REAL  CAMMU(NRAYS), CAMPHI(NRAYS)
Cf2py intent(in) :: CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, PIXEL_INDICES(NRAYS), NRAYS
Cf2py intent(in) NPIX, PIXEL_INDICES, NRAYS
      REAL  PATH(NPIX)
      DOUBLE PRECISION RAY_WEIGHTS(NRAYS)
Cf2py intent(in) :: RAY_WEIGHTS
Cf2py intent(out):: PATH
      REAL WEIGHTS(NPTS)
Cf2py intent(in) :: WEIGHTS
      LOGICAL RETURN_MATRIX
Cf2py intent(in) :: RETURN_MATRIX
      REAL MATRIX(MATRIX_SIZE, NPIX)
Cf2py intent(in,out) :: MATRIX
      INTEGER MATRIX_PTRS(MATRIX_SIZE, NPIX), MATRIX_SIZE
Cf2py intent(in,out) :: MATRIX_PTRS

      REAL TAU
      DOUBLE PRECISION RAY_WEIGHT
      INTEGER N, K,Q
      REAL    MURAY, PHIRAY, MU2, PHI2
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
        TAU = 0.0
        CALL PROJECT_1RAY(NX, NY, NZ, NPTS, NCELLS, WEIGHTS,
     .                      GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                      BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                      GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                X0, Y0, Z0, TAU, RETURN_MATRIX,
     .                      MATRIX, MATRIX_PTRS,MATRIX_SIZE,N,
     .                      RAY_WEIGHT, NPIX)
        PATH(N) = PATH(N) + TAU
900     CONTINUE

      ENDDO

      RETURN
      END

      SUBROUTINE PROJECT_1RAY(NX, NY, NZ, NPTS, NCELLS, WEIGHTS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			                 X0, Y0, Z0, TAU, RETURN_MATRIX,
     .                       MATRIX, MATRIX_PTRS,MATRIX_SIZE,N,
     .                       RAY_WEIGHT, NPIX)

C       Integrates the source function through the extinction field
C     (EXTINCT) backward in the direction (MURAY,PHIRAY) to find the
C     outgoing radiance (RAD) at the point X0,Y0,Z0.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0
      REAL    WEIGHTS(NPTS)
      DOUBLE PRECISION RAY_WEIGHT
      LOGICAL RETURN_MATRIX
      REAL    MATRIX(MATRIX_SIZE, NPIX)
      INTEGER MATRIX_SIZE, MATRIX_PTRS(MATRIX_SIZE, NPIX),N
      INTEGER NPIX, COUNTER
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
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XE, YE, ZE, F2)
        EXTN = 0.0
        DO J=1,8
          EXTN = EXTN+ WEIGHTS(GRIDPTR(J,ICELL))*F2(J)
        ENDDO

        TAU = TAU + SO*0.5*(EXT1+EXTN)*RAY_WEIGHT

        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        IF (.NOT. OUTOFDOMAIN) THEN
          IF (RETURN_MATRIX) THEN
            DO J=1,8
              MATRIX_PTRS(COUNTER, N) = GRIDPTR(J,ICELL)
              MATRIX(COUNTER, N) = MATRIX(COUNTER,N)
     .                + RAY_WEIGHT*SO*0.5*(F1(J) + F2(J))
              COUNTER = COUNTER + 1
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
     .		         UNIFORMZLEV, DELXD, DELYD, DPATH, DPTR)
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
      REAL DPATH(8*(NPX+NPY+NPZ),NPTS)
      INTEGER DPTR(8*(NPX+NPY+NPZ),NPTS)
Cf2py intent(out) :: DPATH, DPTR

      INTEGER SIDE, IP
      LOGICAL VALIDBEAM
      REAL    UNIFZLEV, XO, YO, ZO, DIR, DIRPATH

      DPTR = 0
      DPATH = 0.0

      DO IP = 1, NPTS
        DIRPATH = 0.0
        CALL DIRECT_BEAM_AND_PATHS_PROP(GRIDPOS(1,IP), GRIDPOS(2,IP),
     .            GRIDPOS(3,IP), BCFLAG, XO, YO, ZO, SIDE, NPX, NPY,
     .            NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS, CX, CY, CZ,
     .            CXINV, CYINV, CZINV, DI, DJ, DK, IPDIRECT, DELXD,
     .            DELYD, XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV,
     .            DPATH(:,IP), DPTR(:,IP))
      ENDDO
      RETURN
      END


      SUBROUTINE DIRECT_BEAM_AND_PATHS_PROP (XI, YI, ZI, BCFLAG,
     .           XO, YO, ZO, SIDE, NPX, NPY, NPZ, DELX, DELY,
     .           XSTART, YSTART, ZLEVELS, CX, CY, CZ, CXINV, CYINV,
     .           CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD, XDOMAIN,
     .           YDOMAIN, EPSS, EPSZ, UNIFORMZLEV, DPATH, DPTR)
      IMPLICIT NONE

      REAL    DPATH(8*(NPX+NPY+NPZ))
      INTEGER DPTR(8*(NPX+NPY+NPZ))
      INTEGER BCFLAG, SIDE
      REAL    XO, YO, ZO, XI, YI, ZI

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
        WRITE (6,*) 'DIRECT_BEAM_PROP: Beyond X domain',I,XI,YI,ZI
        STOP
      ENDIF
      J = INT(Y/DELYD) + 1
      IF (J .GT. NPY .AND. ABS(Y-YDOMAIN) .LT. 0.001*DELYD) J = NPY
      IF (J .LT. 1 .OR. J .GT. NPY) THEN
        WRITE (6,*) 'DIRECT_BEAM_PROP: Beyond Y domain',J,XI,YI,ZI
        STOP
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
          WRITE (6,'(A,3I4)') 'DIRECT_BEAM_PROP: beyond grid!', I, J, K
          WRITE(*,*) NPX, NPY, NPZ
          WRITE (6,'(1(2X,3F9.5))') X, Y, Z
          STOP
        ENDIF
C           Get the eight corner extinction values
        I1 = K + NPZ*(J-1)  + NPZ*NPY*(I-1)
        I2 = K + NPZ*(J-1)  + NPZ*NPY*(IP-1)
        I3 = K + NPZ*(JP-1) + NPZ*NPY*(I-1)
        I4 = K + NPZ*(JP-1) + NPZ*NPY*(IP-1)

C           Compute the distance to the next grid plane in  X, Y, and Z
C             If in horizontal uniform region or doing IP then fix X and/or Y.
        IF (ZE .GE. UNIFORMZLEV) THEN
          CONSTX = .TRUE.
          CONSTY = .TRUE.
        ENDIF
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
          WRITE (6,*) 'DIRECT_BEAM_PROP: SO<0', X,Y,Z,
     .         XE,YE,ZE, XP,YP,ZP, CX,CY,CZ, SOX,SOY,SOZ
          STOP
        ENDIF
        SO = MAX(SO,0.0D0)
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
        IF (.NOT. OUTOFDOMAIN .AND. K .LT. NPZ) THEN
          IF (IDP+8 .GT. 8*(NPX+NPY+NPZ)) THEN
            WRITE(*,*) 'ERROR DIRECT_BEAM_AND_PATHS_PROP: ',
     .                 'Max IDP exceeded: ', 8*(NPX+NPY+NPZ)
            STOP
          ENDIF
          DPATH(IDP+1) = MAX(0.0, SO*(U0M*VMWM + 0.5D0*B1 +
     .                 0.3333333333333333D0*C1 - 0.25D0*DU*DV*DW))
          DPATH(IDP+2) = MAX(0.0, SO*(U0*VMWM + 0.5D0*B2 +
     .                 0.3333333333333333D0*C2 + 0.25D0*DU*DV*DW))
          DPATH(IDP+3) = MAX(0.0, SO*(U0M*VWM + 0.5D0*B3 +
     .                 0.3333333333333333D0*C3 + 0.25D0*DU*DV*DW))
          DPATH(IDP+4) = MAX(0.0, SO*(U0*VWM + 0.5D0*B4 +
     .                 0.3333333333333333D0*C4 - 0.25D0*DU*DV*DW))
          DPATH(IDP+5) = MAX(0.0, SO*(U0M*VMW + 0.5D0*B5 +
     .                 0.3333333333333333D0*C5 + 0.25D0*DU*DV*DW))
          DPATH(IDP+6) = MAX(0.0, SO*(U0*VMW + 0.5D0*B6 +
     .                 0.3333333333333333D0*C6 - 0.25D0*DU*DV*DW))
          DPATH(IDP+7) = MAX(0.0, SO*(U0M*VW + 0.5D0*B7 +
     .                 0.3333333333333333D0*C7 - 0.25D0*DU*DV*DW))
          DPATH(IDP+8) = MAX(0.0, SO*(U0*VW + 0.5D0*B8 +
     .                 0.3333333333333333D0*C8 + 0.25D0*DU*DV*DW))

          DPTR(IDP+1) = I1
          DPTR(IDP+2) = I2
          DPTR(IDP+3) = I3
          DPTR(IDP+4) = I4
          DPTR(IDP+5) = I1+1
          DPTR(IDP+6) = I2+1
          DPTR(IDP+7) = I3+1
          DPTR(IDP+8) = I4+1

          IDP = IDP + 8
        ENDIF

        XE = XP + XOFFS
        YE = YP + YOFFS
        ZE = ZP

      ENDDO

      RETURN
      END

      SUBROUTINE PENCIL_BEAM_PROP (X0, Y0, Z0, BCFLAG, IPFLAG,
     .                   SOLARFLUX, SOLARMU, SOLARAZ,
     .                   DIRFLUX, TOTAL_EXT,
     .                   NX, NY, NZ, NCELLS, NPTS, CELLFLAGS,
     .                   XGRID, YGRID, ZGRID, GRIDPOS,GRIDPTR,
     .                   NEIGHPTR, TREEPTR)

C    Calculates the source function for a pencil-beam source.
C    Completely untested.
      IMPLICIT NONE
      INTEGER BCFLAG, IPFLAG
Cf2py intent(in) BCFLAG, IPFLAG
      INTEGER NX, NY, NZ, NPTS, NCELLS, SIDE
Cf2py intent(in) NX,NY,NPTS,NCELLS
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
Cf2py intent(in) GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER*2 CELLFLAGS(NCELLS)
Cf2py intent(in) CELLFLAGS
      LOGICAL VALIDRAD
      REAL    SOLARMU, SOLARAZ, SOLARFLUX
Cf2py intent(in) SOLARMU, SOLARAZ, SOLARFLUX
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
Cf2py intent(in) XGRID, YGRID, ZGRID, GRIDPOS
      REAL    DIRFLUX(NPTS), TOTAL_EXT(NPTS)
Cf2py intent(in) TOTAL_EXT
Cf2py intent(out) DIRFLUX
      DOUBLE PRECISION MU2, PHI2, X0,Y0,Z0, XE,YE,ZE
Cf2py intent(in) X0, Y0,Z0
      DOUBLE PRECISION TRANSMIT

      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
Cf2py intent(out) ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, J, L
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), DONETHIS(8)
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    EXT0, EXT1, EXTN
      REAL    XM,YM
      DOUBLE PRECISION PI, CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
      DOUBLE PRECISION TAUTOL, TAUGRID, S, DELS, TRANSCUT
      DOUBLE PRECISION EXT, TAU, TRANSCELL, ABSCELL
      DOUBLE PRECISION U,V,W, DELX,DELY,DELZ, INVDELX,INVDELY,INVDELZ
      DOUBLE PRECISION F(8), WEIGHTS(8)

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/
      DATA DONEFACE/0,0,0,0,0,0,0,0, 0,1,0,3,0,5,0,7, 2,0,4,0,6,0,8,0,
     .              0,0,1,2,0,0,5,6, 3,4,0,0,5,6,0,0,
     .              0,0,0,0,1,2,3,4, 5,6,7,8,0,0,0,0/

      TRANSMIT = 1.0D0
C         TRANSCUT is the transmission to stop the integration at
      TRANSCUT = 5.0E-5
C         TAUTOL is the maximum optical path for the subgrid intervals
      TAUTOL = 0.2
      DIRFLUX = 0.0

      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)
      PI = ACOS(-1.0D0)
      MU2 = SOLARMU
C         Make the ray direction (opposite to the outgoing direction)
      CX = SQRT(1.0D0-MU2**2)*COS(SOLARAZ-PI)
      CY = SQRT(1.0D0-MU2**2)*SIN(SOLARAZ-PI)
      CZ = -MU2
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
      PRINT *,  'x,y,z', X0,Y0,Z0
C         Start at the desired point, getting the extinction and source there
      XE = X0
      YE = Y0
      ZE = Z0
!      ZE = MAX(MIN(Z0,DBLE(ZGRID(NZ))),DBLE(ZGRID(1)))
      CALL LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, XE, YE, ZE,  ICELL)
c      print '(A,5(1X,F8.5),1X,I6)',
c     .     'INTEGRATE_1RAY:',X0,Y0,Z0,MU2,PHI2,ICELL
      IFACE = 0
      NGRID = 0
C         Loop until reach a Z boundary or transmission is very small
      VALIDRAD = .FALSE.
      DO WHILE (.NOT. VALIDRAD .AND. ICELL .GT. 0)
C           Make sure current cell is valid
        IF (ICELL .LE. 0) THEN
          WRITE (6,*)'INTEGRATE_1RAY: ICELL=',ICELL,
     .                MU2,PHI2,XE,YE,ZE
          STOP
        ENDIF
        NGRID = NGRID + 1

        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XE, YE, ZE, F)
        EXT1 = F(1)*TOTAL_EXT(GRIDPTR(1,ICELL)) +
     .         F(2)*TOTAL_EXT(GRIDPTR(2,ICELL)) +
     .         F(3)*TOTAL_EXT(GRIDPTR(3,ICELL)) +
     .         F(4)*TOTAL_EXT(GRIDPTR(4,ICELL)) +
     .         F(5)*TOTAL_EXT(GRIDPTR(5,ICELL)) +
     .         F(6)*TOTAL_EXT(GRIDPTR(6,ICELL)) +
     .         F(7)*TOTAL_EXT(GRIDPTR(7,ICELL)) +
     .         F(8)*TOTAL_EXT(GRIDPTR(8,ICELL))
        PRINT *, EXT1
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
     .      MU2,PHI2,XE,YE,ZE,SO,ICELL
          STOP
        ENDIF
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ
C           Find the optical path across the grid cell and figure how
C             many subgrid intervals to use

        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XN, YN, ZN, F)
        EXTN = F(1)*TOTAL_EXT(GRIDPTR(1,ICELL)) +
     .         F(2)*TOTAL_EXT(GRIDPTR(2,ICELL)) +
     .         F(3)*TOTAL_EXT(GRIDPTR(3,ICELL)) +
     .         F(4)*TOTAL_EXT(GRIDPTR(4,ICELL)) +
     .         F(5)*TOTAL_EXT(GRIDPTR(5,ICELL)) +
     .         F(6)*TOTAL_EXT(GRIDPTR(6,ICELL)) +
     .         F(7)*TOTAL_EXT(GRIDPTR(7,ICELL)) +
     .         F(8)*TOTAL_EXT(GRIDPTR(8,ICELL))

        TAUGRID = SO*0.5*(EXT1+EXTN)
        NTAU = MAX(1,1+INT(TAUGRID/TAUTOL))
        DELS = SO/NTAU

C           Loop over the subgrid cells
        DO IT = 1, NTAU
          PRINT *, TRANSMIT
          S = IT*DELS
          XI = XE + S*CX
          YI = YE + S*CY
          ZI = ZE + S*CZ

          CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XI, YI, ZI, F)
          EXT0 = F(1)*TOTAL_EXT(GRIDPTR(1,ICELL)) +
     .         F(2)*TOTAL_EXT(GRIDPTR(2,ICELL)) +
     .         F(3)*TOTAL_EXT(GRIDPTR(3,ICELL)) +
     .         F(4)*TOTAL_EXT(GRIDPTR(4,ICELL)) +
     .         F(5)*TOTAL_EXT(GRIDPTR(5,ICELL)) +
     .         F(6)*TOTAL_EXT(GRIDPTR(6,ICELL)) +
     .         F(7)*TOTAL_EXT(GRIDPTR(7,ICELL)) +
     .         F(8)*TOTAL_EXT(GRIDPTR(8,ICELL))

C            Compute the subgrid radiance: integration of the source function
          EXT = 0.5*(EXT0+EXT1)
          IF (EXT .NE. 0.0) THEN
            TAU=EXT*DELS
            ABSCELL = TAU*(1.0-0.5*TAU*(1.0-0.33333333333*TAU))
            TRANSCELL = 1.0 - ABSCELL

          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
          ENDIF
          DO I=1,8
            IF (EXT0 > 1e-5) THEN
              WEIGHTS(I) = F(I)*TOTAL_EXT(GRIDPTR(1,ICELL))/EXT0
            ELSE
              WEIGHTS(I) = F(I)*TOTAL_EXT(GRIDPTR(1,ICELL))*1e5
            ENDIF
            DIRFLUX(GRIDPTR(I,ICELL)) = DIRFLUX(GRIDPTR(I,ICELL)) +
     .            SOLARFLUX*TRANSMIT*WEIGHTS(I)
          ENDDO
          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
C                End of sub grid cell loop
        ENDDO
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
        PRINT *, 'inextcell',INEXTCELL
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
        IF (TRANSMIT .LT. TRANSCUT .OR. NGRID.GT.MAXCELLSCROSS) THEN
          VALIDRAD = .TRUE.
        ELSE IF (INEXTCELL .EQ. 0 .AND. IFACE .GE. 5) THEN
          VALIDRAD = .TRUE.
        ELSE
          ICELL = INEXTCELL
        ENDIF
        XE = XN
        YE = YN
        ZE = ZN
      ENDDO

      SIDE = IFACE

      RETURN
      END

      SUBROUTINE COMPUTE_RADIANCE_GRID(NSTOKES, RADIANCE, RSHPTR,
     .           NPTS, MU, PHI, ML, MM, NSTLEG, NLM, GRIDRAD)

      IMPLICIT NONE
      INTEGER NSTOKES, NPTS, ML, MM, NSTLEG, NLM
!f2py intent(in) :: NSTOKES, NPTS, ML, MM, NSTLEG, NLM
      INTEGER RSHPTR(*)
!f2py intent(in) :: RSHPTR
      REAL MU, PHI, RADIANCE(NSTOKES, *)
!f2py intent(in) :: MU, PHI, RADIANCE
      REAL GRIDRAD(NSTOKES,NPTS)
!f2py intent(out) :: GRIDRAD

      INTEGER I, RS, RE,J
      REAL YLMDIR(NSTLEG,NLM)
      GRIDRAD = 0.0

      CALL YLMALL (.FALSE.,MU,PHI,ML,MM,NSTLEG, YLMDIR)
      DO I=1,NPTS
        RS = RSHPTR(I)
        RE = RSHPTR(I+1)-RS
        DO J=1,RE
          GRIDRAD(1,I) = GRIDRAD(1,I) +
     .      RADIANCE(1,RS+J)*YLMDIR(1,J)
        ENDDO
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
