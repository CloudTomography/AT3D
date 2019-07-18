C     This file containts subroutines that were modified from their original purpose 
C     The original subroutines were written by Frank Evans for the Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer.
C     The modified subroutines were written by Aviad Levis, Technion Institute of Technology, 2019

      SUBROUTINE RENDER (NSTOKES, NX, NY, NZ, 
     .             NPTS, NCELLS, ML, MM, NCS, NLM, NSTLEG, NLEG,  
     .             NUMPHASE, NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, 
     .             SOLARAZ, SFCTYPE, NSFCPAR, SFCGRIDPARMS, 
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, 
     .             UNITS, XGRID, YGRID, ZGRID, GRIDPOS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, 
     .             FLUXES, SHPTR, SOURCE, CAMX, CAMY, CAMZ, CAMMU,  
     .             CAMPHI, NPIX, NPART, TOTAL_EXT, VISOUT,
     .             NSCATANGLE, YLMSUN, PHASETAB)
Cf2py threadsafe
      IMPLICIT NONE
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
      INTEGER ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE
Cf2py intent(in) :: ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE
      INTEGER NMU, NPHI0MAX, NPHI0(*)
Cf2py intent(in) :: NMU, NPHI0MAX, NPHI0
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
Cf2py intent(in) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER SHPTR(*), BCPTR(MAXNBC,2)
Cf2py intent(in) :: SHPTR, BCPTR
      INTEGER*2 CELLFLAGS(*)
      INTEGER IPHASE(NPTS,NPART)
Cf2py intent(in) :: CELLFLAGS, IPHASE
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    SOLARMU, SOLARAZ
Cf2py intent(in) :: SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
Cf2py intent(in) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN
      REAL    MU(*), PHI(NMU,*), WTDO(NMU,*)
Cf2py intent(in) :: MU, PHI, WTDO
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      REAL    SFCGRIDPARMS(*), BCRAD(*)
Cf2py intent(in) :: SFCGRIDPARMS, BCRAD
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    TOTAL_EXT(NPTS), LEGEN(0:NLEG,*)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN, TOTAL_EXT
      REAL    DIRFLUX(*), FLUXES(2,*), SOURCE(*)
Cf2py intent(in) :: DIRFLUX, FLUXES, SOURCE
      DOUBLE PRECISION  CAMX(*), CAMY(*), CAMZ(*)
      REAL CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, NPART
Cf2py intent(in) :: NPIX, NPART
      REAL VISOUT(NPIX)
Cf2py intent(out) :: VISOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS
      INTEGER NSCATANGLE
Cf2py intent(in) ::  NSCATANGLE
      REAL YLMSUN(NLM), PHASETAB(NUMPHASE,NSCATANGLE)
Cf2py intent(in) :: YLMSUN, PHASETAB
 
      
      INTEGER I, J, L, N
      REAL    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2
      DOUBLE PRECISION COSSCAT
      DOUBLE PRECISION U, R, PI 
      REAL, ALLOCATABLE :: YLMDIR(:), SINGSCAT(:)
      DOUBLE PRECISION, ALLOCATABLE :: SUNDIRLEG(:)

      ALLOCATE (YLMDIR(NLM))
      ALLOCATE (SUNDIRLEG(0:NLEG), SINGSCAT(NUMPHASE))
       
C         Make the isotropic radiances for the top boundary
      CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN, 
     .                            UNITS, NTOPPTS, BCRAD(1))
C         Make the bottom boundary radiances for the Lambertian surfaces.  
C          Compute the upwelling bottom radiances using the downwelling fluxes.
      IF (SFCTYPE .EQ. 'FL') THEN
          CALL FIXED_LAMBERTIAN_BOUNDARY(NBOTPTS, BCPTR(1,2), DIRFLUX,
     .                                   FLUXES, SRCTYPE, GNDTEMP,
     .                                   GNDALBEDO, WAVENO, WAVELEN, 
     .                                   UNITS, BCRAD(1+NTOPPTS))
      ELSE IF (SFCTYPE .EQ. 'VL') THEN
          CALL VARIABLE_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2),
     .                                       DIRFLUX, FLUXES, SRCTYPE, 
     .                                       NSFCPAR, SFCGRIDPARMS,
     .                                       BCRAD(1+NTOPPTS))
      ENDIF

      PI = ACOS(-1.0D0)
C         Loop over pixels in image
      DO N = 1, NPIX
        X0 = CAMX(N)
        Y0 = CAMY(N)
        Z0 = CAMZ(N)
        MU2 = CAMMU(N)
        PHI2 = CAMPHI(N)
        MURAY = -MU2
        PHIRAY = PHI2 - PI
  
C           Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (MURAY .GE. 0.0) THEN
              VISOUT(N) = 0.0
              GOTO 900
          ENDIF
          R = (ZGRID(NZ) - Z0)/MURAY
          X0 = X0 + R*SQRT(1-MURAY**2)*COS(PHIRAY)
          Y0 = Y0 + R*SQRT(1-MURAY**2)*SIN(PHIRAY)
          Z0 = ZGRID(NZ)
        ELSE IF (Z0 .LT. ZGRID(1)) THEN
          WRITE (6,*) 'VISUALIZE_RADIANCE: Level below domain'
          STOP
        ENDIF
  
        CALL YLMALL (MU2, PHI2, ML, MM, NCS, YLMDIR)
  
        IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
            COSSCAT = SOLARMU*MU2 + SQRT((1.0-SOLARMU**2)*(1.0-MU2**2))
     .               *COS(SOLARAZ-PHI2)
            IF (NUMPHASE .GT. 0) THEN
                U = (NSCATANGLE-1)*(ACOS(COSSCAT)/PI) + 1
                J = MIN(NSCATANGLE-1,INT(U))
                U = U - J
                DO I = 1, NUMPHASE
                    SINGSCAT(I)=(1-U)*PHASETAB(I,J)+U*PHASETAB(I,J+1)
                ENDDO
            ELSE
                CALL LEGENDRE_ALL (COSSCAT, NLEG, SUNDIRLEG)
                DO L = 0, NLEG
                    SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L+1)/(4*PI)
                ENDDO
            ENDIF
        ENDIF

        CALL INTEGRATE_1RAY (BCFLAG, IPFLAG, 
     .                 NX, NY, NZ, NPTS, NCELLS, 
     .                 GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                 XGRID, YGRID, ZGRID, GRIDPOS,
     .                 ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .                 NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                 DELTAM, SRCTYPE, WAVELEN, SOLARMU,SOLARAZ,
     .                 EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .                 SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG,
     .                 SINGSCAT, MAXNBC, NTOPPTS, NBOTPTS,
     .                 BCPTR, BCRAD, SFCTYPE, NSFCPAR, 
     .                 SFCGRIDPARMS,NPART, MURAY, PHIRAY, 
     .                 MU2, PHI2, X0, Y0, Z0, 
     .                 TOTAL_EXT, VISOUT(N))
C       WRITE(*,*) N, VISOUT(N)
  900   CONTINUE
      ENDDO
  
      RETURN
      END
      
      
      SUBROUTINE RAYLEIGH_PHASE_FUNCTION (WAVELEN, RAYLEGCOEF, 
     .                                    TABLE_TYPE)
      IMPLICIT NONE
      REAL      WAVELEN
Cf2py intent(in) :: WAVELEN
      INTEGER   NCOEF
      PARAMETER (NCOEF=2)
      REAL      RAYLEGCOEF(0:NCOEF)
Cf2py intent(out) :: RAYLEGCOEF
      CHARACTER(LEN=6) :: TABLE_TYPE
Cf2py intent(out) ::  TABLE_TYPE
    
      TABLE_TYPE = 'SCALAR'
      RAYLEGCOEF = (/1.0, 0.0, 0.5/)
      
      RETURN
      END
      
      SUBROUTINE BASE_GRID_PROJECTION (NBCELLS, NCELLS, NBPTS, 
     .                 GRIDPOS, GRIDPTR, TREEPTR, ARRAY, BGARRAY)

      IMPLICIT NONE
      INTEGER NBCELLS, NCELLS, NBPTS
Cf2py intent(in) :: NBCELLS, NCELLS, NBPTS
      REAL    GRIDPOS(3,*)
Cf2py intent(in) :: GRIDPOS
      INTEGER GRIDPTR(8,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, TREEPTR
      DOUBLE PRECISION  ARRAY(*)
Cf2py intent(in) :: ARRAY
      DOUBLE PRECISION BGARRAY(NBPTS)
Cf2py intent(out) :: BGARRAY
      
      DOUBLE PRECISION U, V, W, F(8), DELX, DELY, DELZ
      DOUBLE PRECISION INVDELX, INVDELY, INVDELZ
      INTEGER ICELL, BCELL, I, GPOINT, J, BGPOINT
      
      BGARRAY = ARRAY(1:NBPTS)
      DO ICELL = NBCELLS+1, NCELLS
C       Find base cell for icell
        BCELL = ICELL
        
        DO WHILE (TREEPTR(1, BCELL) .GT. 0)
          BCELL = TREEPTR(1, BCELL)
        ENDDO
        
        DELX = GRIDPOS(1,GRIDPTR(8,BCELL))-GRIDPOS(1,GRIDPTR(1,BCELL))
        DELY = GRIDPOS(2,GRIDPTR(8,BCELL))-GRIDPOS(2,GRIDPTR(1,BCELL))
        DELZ = GRIDPOS(3,GRIDPTR(8,BCELL))-GRIDPOS(3,GRIDPTR(1,BCELL))
      
C       loop over gridpoints belonging to icell and trilin interpolate
        DO I = 1,8
          GPOINT = GRIDPTR(I, ICELL)
          IF ((ARRAY(GPOINT) .NE. 0) .AND. (GPOINT .GT. NBPTS)) THEN
            U = 0.0
            V = 0.0
            W = 0.0
            IF (DELX .GT. 0.0) THEN 
              U = (GRIDPOS(1,GRIDPTR(8, BCELL))-GRIDPOS(1,GPOINT))/DELX
            ENDIF
            IF (DELY .GT. 0.0) THEN 
              V = (GRIDPOS(2,GRIDPTR(8, BCELL))-GRIDPOS(2,GPOINT))/DELY
            ENDIF
            IF (DELY .GT. 0.0) THEN 
              W = (GRIDPOS(3,GRIDPTR(8, BCELL))-GRIDPOS(3,GPOINT))/DELZ
            ENDIF
            F(8) = (1-U) * (1-V) * (1-W)
            F(7) =    U  * (1-V) * (1-W)
            F(6) = (1-U) *    V  * (1-W)
            F(5) =    U  *    V  * (1-W)
            F(4) = (1-U) * (1-V) * W
            F(3) =    U  * (1-V) * W
            F(2) = (1-U) *    V  * W
            F(1) =    U  *    V  * W
            
            DO J = 1,8 
               BGPOINT = GRIDPTR(J, BCELL)
               BGARRAY(BGPOINT) = BGARRAY(BGPOINT)+F(J)*ARRAY(GPOINT)
               ARRAY(GPOINT) = 0
            ENDDO
          ENDIF
        ENDDO
      ENDDO
  
      RETURN
      END
      
      
      
      SUBROUTINE SPACE_CARVE(NX, NY, NZ, NPTS, NCELLS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             BCFLAG, IPFLAG, CAMX, CAMY, CAMZ,
     .             CAMMU, CAMPHI, NPIX, VOLUME)

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
      INTEGER  VOLUME(NPTS)
Cf2py intent(out):: VOLUME

      INTEGER N, K
      REAL    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, R, PI

      PI = ACOS(-1.0D0)
      VOLUME = (/ (0 , K = 1, NPTS) /)
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

        CALL SPACE_CARVE_1RAY(NX, NY, NZ, NPTS, NCELLS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			     X0, Y0, Z0, VOLUME)
900     CONTINUE

      ENDDO

      RETURN
      END
      
      
       SUBROUTINE SPACE_CARVE_1RAY(NX, NY, NZ, NPTS, NCELLS,
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       BCFLAG, IPFLAG, XGRID, YGRID, ZGRID,
     .                       GRIDPOS, MURAY, PHIRAY, MU2, PHI2,
     .			     X0, Y0, Z0, VOLUME)

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
      INTEGER  VOLUME(NPTS)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, K
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), BCFLAG, IPFLAG
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    XM, YM, F0(8), F1(8)
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
      DOUBLE PRECISION U,V,W, DELX,DELY,DELZ, INVDELX,INVDELY,INVDELZ
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
        U = (XE-GRIDPOS(1,IPT1))*INVDELX
        V = (YE-GRIDPOS(2,IPT1))*INVDELY
        W = (ZE-GRIDPOS(3,IPT1))*INVDELZ

        F1(1) = (1-W)*(1-V)*(1-U)
        F1(2) = (1-W)*(1-V)*U
        F1(3) = (1-W)*V*(1-U)
        F1(4) = (1-W)*V*U
        F1(5) = W*(1-V)*(1-U)
        F1(6) = W*(1-V)*U
        F1(7) = W*V*(1-U)
        F1(8) = W*V*U


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

C	If this is not a boundary cell (currently assuming that the bc conditions are open)
C	    Check if the zlevel is within the zparticle levels
C		Compute the gradient field in direction  (MU2,PHI2)

        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        IF (.NOT.OUTOFDOMAIN) THEN
            DO K=1,8
                VOLUME(GRIDPTR(K,ICELL))=1
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
      
      
      
      
      SUBROUTINE GRADIENT(NSTOKES, NX, NY, NZ, NPTS, NBPTS, NCELLS,
     .           NBCELLS, ML, MM, NCS, NLM, NSTLEG, NLEG, NUMPHASE, 
     .           NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .           BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .           SFCTYPE, NSFCPAR, SFCGRIDPARMS, MAXNBC, NTOPPTS,
     .           NBOTPTS, BCPTR, BCRAD, GNDTEMP, GNDALBEDO, SKYRAD,
     .           WAVENO, WAVELEN, UNITS, XGRID, YGRID, ZGRID, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, EXTINCT,  
     .           ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, SHPTR, 
     .           SOURCE, CAMX, CAMY, CAMZ, CAMMU, CAMPHI, NPIX, 
     .           GRADOUT, COST,  MEASUREMENTS, RSHPTR, VISOUT,   
     .           NPART, TOTAL_EXT,RADIANCE, NUMDER, PARTDER, DEXT, 
     .           DALB, DIPHASE, DLEG, NSCATANGLE, YLMSUN, PHASETAB, 
     .           DPHASETAB, DNUMPHASE) 
Cf2py threadsafe
      IMPLICIT NONE
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS
      INTEGER NBPTS, NCELLS, NBCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NBPTS, NCELLS, NBCELLS
      INTEGER ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
Cf2py intent(in) :: ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE
      INTEGER  DNUMPHASE
Cf2py intent(in) :: DNUMPHASE
      INTEGER NMU, NPHI0MAX, NPHI0(*)
Cf2py intent(in) :: NMU, NPHI0MAX, NPHI0
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
Cf2py intent(in) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER SHPTR(*), BCPTR(MAXNBC,2), RSHPTR(*)
Cf2py intent(in) :: SHPTR, BCPTR, RSHPTR
      INTEGER*2 CELLFLAGS(*)
      INTEGER IPHASE(NPTS,NPART)
Cf2py intent(in) :: CELLFLAGS, IPHASE
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    SOLARMU, SOLARAZ
Cf2py intent(in) :: SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
Cf2py intent(in) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN
      REAL    MU(*), PHI(NMU,*), WTDO(NMU,*)
Cf2py intent(in) :: MU, PHI, WTDO
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      REAL    SFCGRIDPARMS(*), BCRAD(*)
Cf2py intent(in) :: SFCGRIDPARMS, BCRAD
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    TOTAL_EXT(NPTS), LEGEN(0:NLEG,*)
Cf2py intent(in) :: EXTINCT, ALBEDO, TOTAL_EXT, LEGEN
      REAL    DIRFLUX(*), FLUXES(2,*), SOURCE(*), RADIANCE(*)
Cf2py intent(in) :: DIRFLUX, FLUXES, SOURCE, RADIANCE
      REAL    CAMX(*), CAMY(*), CAMZ(*), CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX 
      REAL   MEASUREMENTS(*), DLEG(0:NLEG,DNUMPHASE)
      REAL   DEXT(NBPTS,NUMDER), DALB(NBPTS,NUMDER)
      INTEGER DIPHASE(NBPTS,NUMDER)
Cf2py intent(in) :: MEASUREMENTS, DEXT ,DALB, DIPHASE, DLEG
      REAL VISOUT(NPIX)
      DOUBLE PRECISION  GRADOUT(NBPTS,NUMDER), COST
Cf2py intent(out) :: GRADOUT, COST, VISOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS
      INTEGER NUMDER, PARTDER(NUMDER)
Cf2py intent(in) :: NUMDER, PARTDER
      INTEGER NSCATANGLE
      REAL    YLMSUN(NLM), PHASETAB(NUMPHASE, NSCATANGLE)
      REAL    DPHASETAB(DNUMPHASE, NSCATANGLE)
Cf2py intent(in) :: NSCATANGLE, YLMSUN, PHASETAB, DPHASETAB

      INTEGER I, J, L, K
      INTEGER N, M, ME, MS, ND
      DOUBLE PRECISION PIXEL_ERROR, RAYGRAD(NBPTS,NUMDER)
      REAL    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2
      DOUBLE PRECISION THETA0, THETA1, PHIR, PHI0, COSSCAT
      DOUBLE PRECISION U, R, PI
      REAL, ALLOCATABLE :: YLMDIR(:)
      INTEGER, ALLOCATABLE :: LOFJ(:)
      REAL, ALLOCATABLE :: SINGSCAT(:), DSINGSCAT(:)
      DOUBLE PRECISION, ALLOCATABLE :: SUNDIRLEG(:)
      INTEGER MAXSCATANG
      PARAMETER (MAXSCATANG=721)
      
      ALLOCATE (YLMDIR(NLM), LOFJ(NLM), SUNDIRLEG(0:NLEG))
      ALLOCATE (SINGSCAT(NUMPHASE), DSINGSCAT(NUMPHASE))
      
      GRADOUT = 0.0
      
      J = 0
      DO L = 0, ML
        ME = MIN(L,MM)
        MS = (1-NCS)*ME
        DO M = MS, ME
          J = J + 1
          LOFJ(J) = L
        ENDDO
      ENDDO

      
C         Make the isotropic radiances for the top boundary
      CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN, 
     .                            UNITS, NTOPPTS, BCRAD(1))
C         Make the bottom boundary radiances for the Lambertian surfaces.  
C          Compute the upwelling bottom radiances using the downwelling fluxes.
      IF (SFCTYPE .EQ. 'FL') THEN
        CALL FIXED_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2),
     .             DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO, 
     .             WAVENO, WAVELEN, UNITS, BCRAD(1+NTOPPTS))
      ELSE IF (SFCTYPE .EQ. 'VL') THEN
        CALL VARIABLE_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2),
     .               DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,
     .               BCRAD(1+NTOPPTS))
      ENDIF

      PI = ACOS(-1.0D0)
C         Loop over pixels in image
      COST = 0
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
            VISOUT(N) = 0.0
            GOTO 900
          ENDIF
          R = (ZGRID(NZ) - Z0)/MURAY
          X0 = X0 + R*SQRT(1-MURAY**2)*COS(PHIRAY)
          Y0 = Y0 + R*SQRT(1-MURAY**2)*SIN(PHIRAY)
          Z0 = ZGRID(NZ)
        ELSE IF (Z0 .LT. ZGRID(1)) THEN
          WRITE (6,*) 'VISUALIZE_RADIANCE: Level below domain'
          STOP
        ENDIF

        CALL YLMALL (MU2, PHI2, ML, MM, NCS, YLMDIR)

        IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
          COSSCAT = SOLARMU*MU2 + SQRT((1.0-SOLARMU**2)*(1.0-MU2**2))
     .                  *COS(SOLARAZ-PHI2)
          IF (NUMPHASE .GT. 0) THEN
            U = (NSCATANGLE-1)*(ACOS(COSSCAT)/PI) + 1
            J = MIN(NSCATANGLE-1,INT(U))
            U = U - J
            DO I = 1, NUMPHASE
              SINGSCAT(I) = (1-U)*PHASETAB(I,J) + U*PHASETAB(I,J+1)
            ENDDO
            DO I = 1, DNUMPHASE
              DSINGSCAT(I) = (1-U)*DPHASETAB(I,J)+U*DPHASETAB(I,J+1)
            ENDDO
          ELSE
            CALL LEGENDRE_ALL (COSSCAT, NLEG, SUNDIRLEG)
            DO L = 0, NLEG
              SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L+1)/(4*PI)
            ENDDO
          ENDIF
        ENDIF
        
        RAYGRAD = 0.0
        CALL GRAD_INTEGRATE_1RAY (BCFLAG, IPFLAG, 
     .               NX, NY, NZ, NPTS, NCELLS, 
     .               GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .               XGRID, YGRID, ZGRID, GRIDPOS,
     .               ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .               NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .               DELTAM, SRCTYPE, WAVELEN, SOLARMU,SOLARAZ,
     .               EXTINCT(:NPTS,:), ALBEDO(:NPTS,:), LEGEN,
     .               IPHASE(:NPTS,:), DIRFLUX, SHPTR, SOURCE, 
     .               YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .               MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .               SFCTYPE, NSFCPAR, SFCGRIDPARMS,NPART,
     .               MURAY, PHIRAY, MU2, PHI2, X0, Y0, Z0, 
     .               VISOUT(N), RAYGRAD, RSHPTR, TOTAL_EXT, 
     .               RADIANCE, LOFJ, PARTDER, NUMDER,DSINGSCAT, 
     .               DEXT, DALB, DIPHASE, DLEG, NBPTS, DNUMPHASE)
900     CONTINUE
        
        PIXEL_ERROR = VISOUT(N) - MEASUREMENTS(N)

        GRADOUT = GRADOUT + PIXEL_ERROR*RAYGRAD
        COST = COST + 0.5*PIXEL_ERROR**2
C        WRITE(*,*) N, VISOUT(N), MEASUREMENTS(N), COST
      ENDDO
   
      RETURN
      END
      
      
      
      SUBROUTINE GRAD_INTEGRATE_1RAY (BCFLAG, IPFLAG, 
     .              NX, NY, NZ, NPTS, NCELLS, 
     .              GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .              XGRID, YGRID, ZGRID, GRIDPOS,
     .              ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .              NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .              DELTAM, SRCTYPE, WAVELEN,SOLARMU,SOLARAZ,
     .              EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .              DIRFLUX, SHPTR, SOURCE, YLMDIR, 
     .              YLMSUN, SUNDIRLEG, SINGSCAT, MAXNBC, 
     .              NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .              SFCTYPE, NSFCPAR, SFCGRIDPARMS,NPART,
     .              MURAY,PHIRAY, MU2, PHI2, X0, Y0, Z0, 
     .              RADOUT, RAYGRAD, RSHPTR, TOTAL_EXT,  
     .              RADIANCE, LOFJ, PARTDER, NUMDER, DSINGSCAT,
     .              DEXT, DALB, DIPHASE, DLEG, NBPTS, DNUMPHASE)

C       Integrates the source function through the extinction field 
C     (EXTINCT) backward in the direction (MURAY,PHIRAY) to find the 
C     outgoing radiance (RAD) at the point X0,Y0,Z0.
      IMPLICIT NONE
      INTEGER BCFLAG, IPFLAG, NX, NY, NZ, NPTS, NCELLS
      INTEGER ML, MM, NCS, NLM, NLEG, NUMPHASE, NBPTS
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, KK, GRIDPOINT
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SHPTR(NPTS+1), RSHPTR(NPTS+1)
      INTEGER*2 CELLFLAGS(NCELLS)
      INTEGER IPHASE(NPTS,NPART), NPART, DNUMPHASE
      INTEGER BCPTR(MAXNBC,2), LOFJ(NLM)
      LOGICAL DELTAM
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    WAVELEN, SOLARMU, SOLARAZ
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    DIRFLUX(NPTS), TOTAL_EXT(NPTS), LEGEN(0:NLEG,NPTS) 
      DOUBLE PRECISION RAYGRAD(NBPTS,NUMDER)
      REAL    SOURCE(*), RADIANCE(*), DSINGSCAT(DNUMPHASE)
      REAL    YLMDIR(NLM), YLMSUN(NLM), SINGSCAT(NUMPHASE)
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
      REAL    MURAY, PHIRAY, MU2, PHI2, RADOUT
      DOUBLE PRECISION X0, Y0, Z0
      CHARACTER SRCTYPE*1, SFCTYPE*2
      REAL      DLEG(0:NLEG,DNUMPHASE), DEXT(NBPTS,NUMDER)
      REAL      DALB(NBPTS,NUMDER)
      INTEGER   PARTDER(NUMDER), NUMDER, DIPHASE(NBPTS,NUMDER)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, BCELL
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID, K
      INTEGER OPPFACE(6), OLDIPTS(8), DONETHIS(8)
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    OEXTINCT8(8), OSRCEXT8(8), EXTINCT8(8), SRCEXT8(8)
      REAL    EXT0, EXT1, EXTN, SRCEXT0, SRCEXT1, RADBND
      REAL    XM,YM, GRAD8(8,NUMDER), OGRAD8(8,NUMDER)
      REAL    GRAD0(8,NUMDER), GRAD1(8,NUMDER), SRCGRAD(8,NUMDER)
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, F(8), FB(8)
      DOUBLE PRECISION TAUTOL, TAUGRID, S, DELS, TRANSCUT, ABSCELL1
      DOUBLE PRECISION EXT, SRC, TAU, TRANSCELL,ABSCELL, TRANSMIT, RAD
      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/
      DATA DONEFACE/0,0,0,0,0,0,0,0, 0,1,0,3,0,5,0,7, 2,0,4,0,6,0,8,0,
     .              0,0,1,2,0,0,5,6, 3,4,0,0,5,6,0,0,
     .              0,0,0,0,1,2,3,4, 5,6,7,8,0,0,0,0/


C         TRANSCUT is the transmission to stop the integration at
      TRANSCUT = 5.0E-5
C         TAUTOL is the maximum optical path for the subgrid intervals
      TAUTOL = 0.2

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
      RAD = 0.0D0
      TRANSMIT = 1.0D0

C         Loop until reach a Z boundary or transmission is very small
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
      
        CALL GET_BASE_GRID_CELL(BCELL, ICELL, TREEPTR) 
        
C           Make sure current cell is valid
        IF (ICELL .LE. 0) THEN
          WRITE (6,*)'INTEGRATE_1RAY: ICELL=',ICELL,
     .                MURAY,PHIRAY,XE,YE,ZE
          STOP
        ENDIF
        NGRID = NGRID + 1

C           Decide which of the eight grid points we need the source function
        DO I = 1, 8
          DONETHIS(I) = DONEFACE(I,IFACE+1)
          IF (NX .EQ. 1 .AND. ONEX(I) .LT. 0) DONETHIS(I) = ONEX(I)
          IF (NY .EQ. 1 .AND. ONEY(I) .LT. 0) DONETHIS(I) = ONEY(I)
          OEXTINCT8(I) = EXTINCT8(I)
          OSRCEXT8(I) = SRCEXT8(I)
          OGRAD8(I,:) = GRAD8(I,:)
        ENDDO
     
C         Compute the source function times extinction in direction (MU2,PHI2)
C         In addition compute the gradient field in direction (MU2, PHI2)
        CALL COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR, ML,
     .            MM, NCS, NLM, NLEG, NUMPHASE, NPTS, DELTAM, 
     .            SRCTYPE, SOLARMU, EXTINCT, ALBEDO, LEGEN, 
     .            IPHASE, DIRFLUX, SHPTR, RSHPTR, SOURCE, 
     .            YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT, DONETHIS,  
     .            OLDIPTS, OEXTINCT8, OSRCEXT8, EXTINCT8, 
     .            SRCEXT8, TOTAL_EXT, NPART, RADIANCE, OGRAD8,
     .            GRAD8, LOFJ, CELLFLAGS, PARTDER, NUMDER, DNUMPHASE,
     .            DEXT, DALB, DIPHASE, DLEG, NBPTS, BCELL, DSINGSCAT)
     
C         Interpolate the source and extinction to the current point
        CALL GET_INTERP_KERNEL(ICELL,GRIDPTR,GRIDPOS,XE,YE,ZE,F)
        SRCEXT1 = F(1)*SRCEXT8(1) + F(2)*SRCEXT8(2) +
     .            F(3)*SRCEXT8(3) + F(4)*SRCEXT8(4) +
     .            F(5)*SRCEXT8(5) + F(6)*SRCEXT8(6) +
     .            F(7)*SRCEXT8(7) + F(8)*SRCEXT8(8)
        EXT1 = F(1)*EXTINCT8(1) + F(2)*EXTINCT8(2) +
     .         F(3)*EXTINCT8(3) + F(4)*EXTINCT8(4) +
     .         F(5)*EXTINCT8(5) + F(6)*EXTINCT8(6) +
     .         F(7)*EXTINCT8(7) + F(8)*EXTINCT8(8)

        CALL GET_INTERP_KERNEL(BCELL,GRIDPTR,GRIDPOS,XE,YE,ZE,FB)
        DO KK=1,8
          GRAD1(KK,:) = FB(KK)*GRAD8(KK,:)
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
        CALL GET_INTERP_KERNEL(ICELL,GRIDPTR,GRIDPOS,XN,YN,ZN,F)
        EXTN = F(1)*EXTINCT8(1) + F(2)*EXTINCT8(2) +
     .         F(3)*EXTINCT8(3) + F(4)*EXTINCT8(4) +
     .         F(5)*EXTINCT8(5) + F(6)*EXTINCT8(6) +
     .         F(7)*EXTINCT8(7) + F(8)*EXTINCT8(8)
     
        TAUGRID = SO*0.5*(EXT1+EXTN)
        NTAU = MAX(1,1+INT(TAUGRID/TAUTOL))
        DELS = SO/NTAU 

C           Loop over the subgrid cells
        DO IT = 1, NTAU
          S = IT*DELS
          IF (IPINX) THEN
            XI = XE
          ELSE
            XI = XE + S*CX
          ENDIF
          IF (IPINY) THEN
            YI = YE
          ELSE
            YI = YE + S*CY
          ENDIF
          ZI = ZE + S*CZ
          
C            Interpolate extinction and source function along path
          CALL GET_INTERP_KERNEL(ICELL,GRIDPTR,GRIDPOS,XI,YI,ZI,F)
          SRCEXT0 = F(1)*SRCEXT8(1) + F(2)*SRCEXT8(2) +
     .              F(3)*SRCEXT8(3) + F(4)*SRCEXT8(4) +
     .              F(5)*SRCEXT8(5) + F(6)*SRCEXT8(6) +
     .              F(7)*SRCEXT8(7) + F(8)*SRCEXT8(8)
          IF (IT .NE. NTAU) THEN
            EXT0 = F(1)*EXTINCT8(1) + F(2)*EXTINCT8(2) +
     .             F(3)*EXTINCT8(3) + F(4)*EXTINCT8(4) +
     .             F(5)*EXTINCT8(5) + F(6)*EXTINCT8(6) +
     .             F(7)*EXTINCT8(7) + F(8)*EXTINCT8(8)
          ELSE
            EXT0 = EXTN
          ENDIF
          SRCEXT0 = MAX(0.0,SRCEXT0)
          
          CALL GET_INTERP_KERNEL(BCELL,GRIDPTR,GRIDPOS,XI,YI,ZI,FB)
          DO KK=1,8
            GRAD0(KK,:) = FB(KK)*GRAD8(KK,:)
          ENDDO
        
C            Compute the subgrid radiance: integration of the source function
          EXT = 0.5*(EXT0+EXT1)
          IF (EXT .NE. 0.0) THEN
            TAU=EXT*DELS
            ABSCELL = TAU*(1.0-0.5*TAU*(1.0-0.33333333333*TAU))
            TRANSCELL = 1.0 - ABSCELL
C                 Linear extinction, linear source*extinction, to second order
            SRC = ( 0.5*(SRCEXT0+SRCEXT1) 
     .                + 0.08333333333*(EXT0*SRCEXT1-EXT1*SRCEXT0)*DELS
     .                    *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
          
            SRCGRAD = ( 0.5*(GRAD0+GRAD1) 
     .        + 0.08333333333*(EXT0*GRAD1-EXT1*GRAD0)*DELS
     .         *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT

          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
            SRC = 0.0
            SRCGRAD = 0.0
          ENDIF
            
          RAD = RAD + TRANSMIT*SRC*ABSCELL
          DO KK = 1, 8
            GRIDPOINT = GRIDPTR(KK, BCELL)
            RAYGRAD(GRIDPOINT,:) = RAYGRAD(GRIDPOINT,:) +
     .        TRANSMIT*SRCGRAD(KK,:)*ABSCELL
          ENDDO
          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          SRCEXT1 = SRCEXT0
          GRAD1 = GRAD0

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
          DONE = .TRUE.
        ELSE IF (INEXTCELL .EQ. 0) THEN
          DONE = .TRUE.
          CALL FIND_BOUNDARY_RADIANCE (XN, YN, MU2, PHI2, 
     .                      IC, KFACE, GRIDPTR, GRIDPOS,
     .                      MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                      NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                      SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, DIRFLUX, 
     .                      SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                      RADBND)
          
          RAD = RAD + TRANSMIT*RADBND
        ELSE
          XE = XN
          YE = YN
          ZE = ZN
          ICELL = INEXTCELL
        ENDIF

      ENDDO

      RADOUT = RAD
      RETURN
      END
      
      SUBROUTINE COMPUTE_SOURCE_GRAD_1CELL (ICELL,  
     .            GRIDPTR, ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .            NPTS, DELTAM, SRCTYPE, SOLARMU, EXTINCT,
     .            ALBEDO, LEGEN, IPHASE, DIRFLUX, SHPTR, 
     .            RSHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, 
     .            SINGSCAT, DONETHIS, OLDIPTS, OEXTINCT8,
     .            OSRCEXT8, EXTINCT8, SRCEXT8, TOTAL_EXT, NPART,
     .            RADIANCE, OGRAD8, GRAD8, LOFJ, CELLFLAGS, PARTDER,
     .            NUMDER, DNUMPHASE, DEXT, DALB, DIPHASE, DLEG, 
     .            NBPTS, BCELL, DSINGSCAT)
C       Computes the source function times extinction for gridpoints 
C     belonging to cell ICELL in the direction (MU,PHI).  The results 
C     are returned in SRCEXT8 and EXTINCT8.
C     The spherical harmonic source function series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
      IMPLICIT NONE
      INTEGER ICELL, NPTS, ML, MM, NCS, NLM, NLEG, BCELL
      INTEGER NUMPHASE, NPART, NUMDER, NBPTS
Cf2py intent(in) :: ICELL, NPTS, ML, MM, NCS, NLM, NLEG, NUMPHASE
      INTEGER GRIDPTR(8,*), SHPTR(*), RSHPTR(*)
Cf2py intent(in) :: GRIDPTR, SHPTR      
      INTEGER DONETHIS(8), OLDIPTS(8)
Cf2py intent(in) :: DONETHIS, OLDIPTS      
      INTEGER IPHASE(NPTS,NPART), DNUMPHASE
Cf2py intent(in) :: IPHASE    
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    SOLARMU
Cf2py intent(in) :: SOLARMU
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(0:NLEG,*),  TOTAL_EXT(NPTS)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN
      REAL    DIRFLUX(*), SOURCE(*), RADIANCE(*)
Cf2py intent(in) :: DIRFLUX, SOURCE
      REAL    YLMDIR(*), YLMSUN(*), SINGSCAT(*)
      REAL    DSINGSCAT(DNUMPHASE)
Cf2py intent(in) :: YLMDIR, YLMSUN, SINGSCAT
      REAL    OEXTINCT8(8), OSRCEXT8(8)
Cf2py intent(in) :: OEXTINCT8, OSRCEXT8
      REAL    EXTINCT8(8), SRCEXT8(8), EXT
Cf2py intent(out) :: EXTINCT8, SRCEXT8
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
Cf2py intent(in) :: SUNDIRLEG
      CHARACTER SRCTYPE*1
Cf2py intent(in) :: SRCTYPE
      INTEGER*2 CELLFLAGS(*)
      LOGICAL OUTOFDOMAIN
      REAL    GRAD8(8,NUMDER), OGRAD8(8,NUMDER)
      REAL    TRUNC_SINGSCAT, FULL_SINGSCAT
      REAL    DLEG(0:NLEG,DNUMPHASE), DEXT(NBPTS,NUMDER)
      REAL    DALB(NBPTS,NUMDER), DEXTM, DALBM, DLEGM
      INTEGER DIPHASE(NBPTS,NUMDER), KK
      INTEGER IP, J, L, M, MS, ME, K, IS, NS, N, I, IB
      INTEGER IPA, LOFJ(NLM), PARTDER(NUMDER), RNS, RIS, IDR
      DOUBLE PRECISION DA, F, A, SECMU0, W
      
      GRAD8 = 0.0
      SECMU0 = 1.0D0/ABS(SOLARMU)

C         Loop over the grid points, computing the source function 
C           at the viewing angle from the spherical harmonic source function.
      DO N = 1, 8
        IP = GRIDPTR(N,ICELL)
        IB = GRIDPTR(N,BCELL)
        
        I = DONETHIS(N)
        IF (I .GT. 0 .AND. IP .EQ. OLDIPTS(N)) THEN 
          EXTINCT8(N) = OEXTINCT8(I)
          SRCEXT8(N) = OSRCEXT8(I)
          GRAD8(N,:) = OGRAD8(I,:)
        ELSE IF (I .LT. 0) THEN
          EXTINCT8(N) = EXTINCT8(ABS(I))
          SRCEXT8(N) = SRCEXT8(ABS(I))
          GRAD8(N,:) = GRAD8(ABS(I),:)
        ELSE
        
          EXT = TOTAL_EXT(IP)
          OLDIPTS(N) = IP
          IS = SHPTR(IP)
          NS = SHPTR(IP+1)-IS

          OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                   BTEST(INT(CELLFLAGS(ICELL)),1))
           
          IF (.NOT.OUTOFDOMAIN) THEN 
            RIS = RSHPTR(IP)
            RNS = RSHPTR(IP+1)-RIS
            DO IDR = 1, NUMDER
            
              IPA = PARTDER(IDR)
              IF (NUMPHASE .GT. 0) THEN
                K = IPHASE(IP,IPA)
              ELSE
                K = IP
              ENDIF
              
              DA = DIRFLUX(IP)*SECMU0
              FULL_SINGSCAT = 0.0
              
              IF (DELTAM) THEN
                F = LEGEN(ML+1,K)
              ELSE 
                F = 0.0
              ENDIF 
              
              DEXTM = (1.0-ALBEDO(IP,IPA)*F) * DEXT(IB,IDR)
              DALBM = (1.0-F)*DALB(IB,IDR)/(1.0-ALBEDO(IP,IPA)*F)
              
              DO J = 1, RNS
                L = LOFJ(J)
                DLEGM = DLEG(L,DIPHASE(IB,IDR))/(1-F)
                
                GRAD8(N,IDR) = GRAD8(N,IDR) + 
     .            RADIANCE(RIS+J)*YLMDIR(J)*(
     .              DEXTM*(ALBEDO(IP,IPA)*LEGEN(L,K)-1.0) +
     .              EXTINCT(IP,IPA)*DALBM*LEGEN(L,K) + 
     .              EXTINCT(IP,IPA)*ALBEDO(IP,IDR)*DLEGM)
     
                IF (NUMPHASE .EQ. 0) THEN
                  IF (L .LE. ML) THEN
                    A = LEGEN(L,K) + F/(1-F)
                  ELSE
                    A = LEGEN(L,K)/(1-F)
                  ENDIF
                  FULL_SINGSCAT = FULL_SINGSCAT + DA*(
     .               DEXTM*ALBEDO(IP,IPA)*A*SUNDIRLEG(L) + 
     .               EXTINCT(IP,IPA)*DALBM*A*SUNDIRLEG(L) +
     .               EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM*SUNDIRLEG(L))
                ENDIF
              ENDDO

C             Add in the single scattering contribution for the original unscaled phase function.  
C             For DELTA M and L<=ML this requires reconstructing the original phase function values (LEGEN) by unscaling.  
C             Need to put the inverse of the tau scaling in the source function because extinction is still scaled.
              IF (NUMPHASE .GT. 0) THEN
                FULL_SINGSCAT = DA*(
     .               DEXTM*ALBEDO(IP,IPA)*SINGSCAT(K) + 
     .               EXTINCT(IP,IPA)*DALBM*SINGSCAT(K) +
     .               EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*
     .               DSINGSCAT(DIPHASE(IB,IDR)))
              ENDIF
              GRAD8(N,IDR) = GRAD8(N,IDR) + FULL_SINGSCAT
            ENDDO
          ENDIF
          
C             Sum over the spherical harmonic series of the source function
          SRCEXT8(N) = 0.0
          DO J = 1, NS
            SRCEXT8(N) = SRCEXT8(N) + SOURCE(IS+J)*YLMDIR(J)
          ENDDO

          
C             Special case for solar source and Delta-M
          IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
            
            DO IPA = 1, NPART
              IF (EXT.EQ.0.0) THEN
                W = 1.0D0
              ELSE
                W = EXTINCT(IP,IPA)/EXT
              ENDIF
              IF (W.EQ.0.0) CYCLE
            
              IF (NUMPHASE .GT. 0) THEN
                K = IPHASE(IP,IPA)
              ELSE
                K = IP
              ENDIF

C               First subtract off the truncated single scattering 
              DA = ALBEDO(IP,IPA)*DIRFLUX(IP)*SECMU0*W
              J = 1
              DO L = 0, ML
                ME = MIN(L,MM)
                MS = (1-NCS)*ME
                A = DA*LEGEN(L,K)
                IF (J .LE. NS) THEN
                  DO M = MS, ME
                    SRCEXT8(N) = SRCEXT8(N) - A*YLMDIR(J)*YLMSUN(J)
                    J = J + 1
                  ENDDO
                ENDIF
              ENDDO
C               Then add in the single scattering contribution for the
C               original unscaled phase function.  For L<=ML this requires
C               reconstructing the original phase function values (LEGEN) by
C               unscaling.  Need to put the inverse of the tau scaling in the
C               source function because extinction is still scaled.
              IF (NUMPHASE .GT. 0) THEN
                SRCEXT8(N) = SRCEXT8(N) + DA*SINGSCAT(K)
              ELSE
                F = LEGEN(ML+1,K)
                DO L = 0, NLEG
                  IF (L .LE. ML) THEN
                    A = DA*(LEGEN(L,K) + F/(1-F))
                  ELSE
                    A = DA*LEGEN(L,K)/(1-F)
                  ENDIF
                  SRCEXT8(N) = SRCEXT8(N) + A*SUNDIRLEG(L)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
          SRCEXT8(N) = SRCEXT8(N)*EXT
          EXTINCT8(N) = EXT
        ENDIF
      ENDDO
 
      RETURN
      END
       
      
      SUBROUTINE GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, 
     .                             X, Y, Z, F)
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
      
      
      SUBROUTINE DIRECT_BEAM_GRAD (XI, YI, ZI, BCFLAG, IPFLAG, 
     .                SOLARFLUX, SOLARMU, SOLARAZ, XO, YO, ZO,
     .                DIRPATH, SIDE, VALIDBEAM, NPX, NPY, NPZ,
     .                DELX, DELY, XSTART, YSTART, ZLEVELS,
     .                EXTDIRP, CX, CY, CZ, CXINV, CYINV, CZINV,
     .                DI, DJ, DK, IPDIRECT, DELXD, DELYD,XDOMAIN,
     .                YDOMAIN, EPSS, EPSZ, UNIFORMZLEV, DIRGRAD)
      IMPLICIT NONE
      INTEGER BCFLAG, IPFLAG, SIDE 
      LOGICAL VALIDBEAM
      REAL    XI, YI, ZI, SOLARFLUX, SOLARMU, SOLARAZ
      REAL    DIRFLUX, XO, YO, ZO, DIRPATH, DIRGRAD
      INTEGER IX, IY, IZ, JZ, IL, IM, IU
      INTEGER I, J, K, IP, JP, I1, I2, I3, I4, IPA
      INTEGER IPDIRECT, DI, DJ, DK
      LOGICAL CONSTX, CONSTY, HITBOUNDARY, BTEST
      DOUBLE PRECISION EXTINCT, ALBEDO, F
      DOUBLE PRECISION SUNMU, SUNAZ, PATH
      DOUBLE PRECISION UNIFORMZLEV
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION EPSS, EPSZ, XDOMAIN, YDOMAIN
      DOUBLE PRECISION XOFFS, YOFFS, DELXD, DELYD
      DOUBLE PRECISION X, Y, Z, XE, YE, ZE, XP, YP, ZP
      DOUBLE PRECISION X0, X1, Y0, Y1, Z0, Z1, SO, SOX, SOY, SOZ
      DOUBLE PRECISION U0, V0, W0, U1, V1, W1, AX, AY, AZ
      DOUBLE PRECISION U0M, V0M, W0M, U1M, V1M, W1M, DU, DV, DW
      DOUBLE PRECISION E1,E2,E3,E4,E5,E6,E7,E8,  A, B, C, D
      DOUBLE PRECISION B1,B2,B3,B4,B5,B6,B7,B8,C1,C2,C3,C4,C5,C6,C7,C8 
      DOUBLE PRECISION UV,UMV,UVM,UMVM,UW,UMW,UWM,UMWM,VW,VMW,VWM,VMWM
      DOUBLE PRECISION VWU,VWUM,UWV,UWVM,UVW,UVWM
      
      INTEGER NPX, NPY, NPZ
      REAL DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)
      REAL EXTDIRP(*)


!         Here for computing the direct beam path for one starting point.
      Z = ZI
      X = XI - XSTART
      Y = YI - YSTART

!         Find the grid location 
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
     .  .OR. ABS(X-(NPX-1)*DELXD) .LT. 0.01*DELXD))  CONSTX = .TRUE.
      IF (BTEST(BCFLAG,1) .AND. (ABS(Y) .LT. 0.01*DELYD  
     .  .OR. ABS(Y-(NPY-1)*DELYD) .LT. 0.01*DELYD))  CONSTY = .TRUE.

       ! If have multiple subdomains (processors) and ray is going outwards 
       !   from a boundary then set the HITBOUNDARY flag
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

!           Grid cell loop begin
      PATH = DIRPATH
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
!           Get the eight corner extinction values
        I1 = K + NPZ*(J-1)  + NPZ*NPY*(I-1)
        I2 = K + NPZ*(J-1)  + NPZ*NPY*(IP-1)
        I3 = K + NPZ*(JP-1) + NPZ*NPY*(I-1)
        I4 = K + NPZ*(JP-1) + NPZ*NPY*(IP-1)
        E1 = EXTDIRP(I1)
        E2 = EXTDIRP(I2)
        E3 = EXTDIRP(I3)
        E4 = EXTDIRP(I4)
        E5 = EXTDIRP(I1+1)
        E6 = EXTDIRP(I2+1)
        E7 = EXTDIRP(I3+1)
        E8 = EXTDIRP(I4+1)

!           Compute the distance to the next grid plane in  X, Y, and Z
!             If in horizontal uniform region or doing IP then fix X and/or Y.
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
        
!           The shortest distance is the plane we stop at:
!             get the exitting location and increment the cell
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
!             If have reached a horizontal boundary then either wrap around
!               (periodic) or go into IP mode (open boundaries).
          IF (I .EQ. 0) THEN
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
     .        XE,YE,ZE, XP,YP,ZP, CX,CY,CZ, SOX,SOY,SOZ
          STOP
        ENDIF
        SO = MAX(SO,0.0D0)
!           Make the starting and ending interpolation factors
        AX = 1.0D0/(X1-X0)
        AY = 1.0D0/(Y1-Y0)
        AZ = 1.0D0/(Z1-Z0)
        U0 = (XE-X0)*AX
        V0 = (YE-Y0)*AY
        W0 = (ZE-Z0)*AZ
        U1 = (XP-X0)*AX
        V1 = (YP-Y0)*AY
        W1 = (ZP-Z0)*AZ
!           Compute the cubic polynomial extinction coefficients
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
        A =  (E1*U0M + E2*U0)*VMWM 
     .     + (E3*U0M + E4*U0)*VWM 
     .     + (E5*U0M + E6*U0)*VMW 
     .     + (E7*U0M + E8*U0)*VW
        B1 = -DU*VMWM - DV*UMWM - DW*UMVM
        B2 =  DU*VMWM - DV*UWM  - DW*UVM
        B3 = -DU*VWM  + DV*UMWM - DW*UMV
        B4 =  DU*VWM  + DV*UWM  - DW*UV
        B5 = -DU*VMW  - DV*UMW  + DW*UMVM
        B6 =  DU*VMW  - DV*UW   + DW*UVM
        B7 = -DU*VW   + DV*UMW  + DW*UMV
        B8 =  DU*VW   + DV*UW   + DW*UV
        B = B1*E1 +B2*E2 +B3*E3 +B4*E4 +B5*E5 +B6*E6 +B7*E7 +B8*E8 
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
        C = C1*E1 +C2*E2 +C3*E3 +C4*E4 +C5*E5 +C6*E6 +C7*E7 +C8*E8 
        D = DU*DV*DW*(E2+E3+E5+E8-E1-E4-E6-E7)
!            Compute the path through the cell: integration of extinction
        PATH = PATH + SO*(A +0.5D0*B +0.3333333333333333D0*C +0.25D0*D)

        XE = XP + XOFFS
        YE = YP + YOFFS
        ZE = ZP
      ENDDO
      DIRFLUX = SOLARFLUX*EXP(-PATH)
      DIRPATH = PATH
      XO=XE+XSTART ; YO=YE+YSTART ; ZO=ZE
      VALIDBEAM = ABS(ZE-ZLEVELS(NPZ)) .LT. EPSZ
      IF (VALIDBEAM) SIDE=6
      RETURN
      END
      
      
      SUBROUTINE PRECOMPUTE_PHASE_CHECK(NSCATANGLE, NUMPHASE, ML, NLM,
     .                          NLEG, LEGEN, PHASETAB, DELTAM)
C       Precomputes the phase function as a function of scattering angle
C     for all the tabulated phase functions.
      IMPLICIT NONE
      INTEGER NSCATANGLE, ML, NLM, NLEG, NUMPHASE
Cf2py intent(in) :: NSCATANGLE, ML, NLM, NLEG, NUMPHASE
      REAL    LEGEN(0:NLEG,NUMPHASE)
Cf2py intent(in) :: LEGEN
      REAL    PHASETAB(NUMPHASE,NSCATANGLE)
Cf2py intent(out) :: PHASETAB
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      INTEGER I, J, L, MAXLEG
      DOUBLE PRECISION PI, SUM, F, A, COSSCAT, LEGSCAT(0:NLEG)
      
      PI = ACOS(-1.0D0)
      DO J = 1, NSCATANGLE
        COSSCAT = COS(PI*DFLOAT(J-1)/(NSCATANGLE-1))
C       Compute the Legendre polynomials for the scattering angle 
C       for the untruncated solar single scattering computation.
        CALL LEGENDRE_ALL (COSSCAT, NLEG, LEGSCAT)
        DO L = 0, NLEG
          LEGSCAT(L) = LEGSCAT(L)*(2*L+1)/(4*PI)
        ENDDO
        DO I = 1, NUMPHASE
          SUM = 0.0D0  
          F = LEGEN(ML+1,I)
          DO L = 0, NLEG
            IF (L .LE. ML .AND. DELTAM) THEN
              A = (LEGEN(L,I) + F/(1-F))
            ELSE IF (DELTAM) THEN
              A = LEGEN(L,I)/(1-F)
            ELSE 
              A = LEGEN(L,I)
            ENDIF
            SUM = SUM + A*LEGSCAT(L)
          ENDDO
          IF (SUM .LT. 0.0) THEN
            WRITE (6,*) 'PRECOMPUTE_PHASE: negative source',
     .          ' function for tabulated phase function: ',
     .          I, J, SUM
            STOP
          ENDIF
          PHASETAB(I,J) = SNGL(SUM)
        ENDDO
      ENDDO
       

      RETURN
      END