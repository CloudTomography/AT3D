!  SHDOM: Spherical harmonic discrete ordinate radiative transfer method.
!      See shdom.txt for documentation.
!      Fortran 90 version of the main program for using allocatable arrays.

      SUBROUTINE TRILIN_INTERP_PROP (X, Y, Z, INIT, NSTLEG, NLEG, &
                     TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE, &
                     NPX, NPY, NPZ, NUMPHASE, DELX, DELY, &
                     XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP, &
                     ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD, &
                     ZCKD, GASABS, EXTMIN, SCATMIN)
!      Trilinearly interpolates the quantities on the input property
!     grid at the single point (X,Y,Z) to get the output TEMP,EXTINCT,
!     ALBEDO, and LEGEN or IPHASE.  Interpolation is done on the 
!     volume coefficients.  Also adds in the separate gaseous absorption.
!     Divides the phase function Legendre coefficients LEGEN by 2*l+1.
!     The phase function pointer IPHASE (for tabulated phase functions) 
!     is that of the maximum weighted scattering property grid point.
!     If INIT=.TRUE. then transfers the tabulated phase functions.
      IMPLICIT NONE
      INTEGER NSTLEG, NLEG
      INTEGER IPHASE
      LOGICAL INIT
      REAL    X, Y, Z,  TEMP, EXTINCT, ALBEDO, LEGEN(NSTLEG,0:NLEG,*)
      INTEGER IX, IXP, IY, IYP, IZ, L, IL, IM, IU, J
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I
      DOUBLE PRECISION U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
      DOUBLE PRECISION SCAT1,SCAT2,SCAT3,SCAT4,SCAT5,SCAT6,SCAT7,SCAT8
      DOUBLE PRECISION SCATTER, MAXSCAT, KG, EXTMIN, SCATMIN
      
      INTEGER NPX, NPY, NPZ
      INTEGER NUMPHASE
      REAL DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)
      REAL TEMPP(*), EXTINCTP(*), ALBEDOP(*)
      REAL LEGENP(*), EXTDIRP(*)
      INTEGER IPHASEP(*)
      INTEGER NZCKD
      REAL ZCKD(*), GASABS(*)
      
      IF (INIT) THEN
!         If there are tabulated phase functions, then transfer them
        DO I = 1, NUMPHASE
          DO L = 0, NLEG
            DO J = 1, NSTLEG
              LEGEN(J,L,I) = LEGENP(J+NSTLEG*(L+(NLEG+1)*(I-1)))/(2*L+1)
            ENDDO
          ENDDO
        ENDDO
        EXTMIN = 1.0E-5 / ( (ZLEVELS(NPZ)-ZLEVELS(1))/NPZ)
        SCATMIN = 0.1*EXTMIN
        RETURN
      ENDIF

!         Find the grid location and compute the interpolation factors
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
      IZ = MAX(IL,1)
      W = DBLE(Z - ZLEVELS(IZ))/(ZLEVELS(IZ+1) - ZLEVELS(IZ))
      W = MAX( MIN( W, 1.0D0), 0.0D0)
      IX = INT((X-XSTART)/DELX) + 1
      IF (ABS(X-XSTART-NPX*DELX) .LT. 0.01*DELX) IX = NPX
      IF (IX .LT. 1 .OR. IX .GT. NPX) THEN
        WRITE (6,*) 'TRILIN: Beyond X domain',IX,NPX,X,XSTART
        STOP
      ENDIF
      IXP = MOD(IX,NPX) + 1
      U = DBLE(X-XSTART-DELX*(IX-1))/DELX
      U = MAX( MIN( U, 1.0D0), 0.0D0)
      IF (U .LT. 1.0D-5) U = 0.0D0
      IF (U .GT. 1.0D0-1.0D-5) U = 1.0D0
      IY = INT((Y-YSTART)/DELY) + 1
      IF (ABS(Y-YSTART-NPY*DELY) .LT. 0.01*DELY) IY = NPY
      IF (IY .LT. 1 .OR. IY .GT. NPY) THEN
        WRITE (6,*) 'TRILIN: Beyond Y domain',IY,NPY,Y,YSTART
        STOP
      ENDIF
      IYP = MOD(IY,NPY) + 1
      V = DBLE(Y-YSTART-DELY*(IY-1))/DELY
      V = MAX( MIN( V, 1.0D0), 0.0D0)
      IF (V .LT. 1.0D-5) V = 0.0D0
      IF (V .GT. 1.0D0-1.0D-5) V = 1.0D0

      F1 = (1-U)*(1-V)*(1-W)
      F2 =    U *(1-V)*(1-W)
      F3 = (1-U)*   V *(1-W)
      F4 =    U *   V *(1-W)
      F5 = (1-U)*(1-V)*   W
      F6 =    U *(1-V)*   W
      F7 = (1-U)*   V *   W
      F8 =    U *   V *   W
      I1 = IZ + NPZ*(IY-1) + NPZ*NPY*(IX-1)
      I2 = IZ + NPZ*(IY-1) + NPZ*NPY*(IXP-1)
      I3 = IZ + NPZ*(IYP-1) + NPZ*NPY*(IX-1)
      I4 = IZ + NPZ*(IYP-1) + NPZ*NPY*(IXP-1)
      I5 = I1+1
      I6 = I2+1
      I7 = I3+1
      I8 = I4+1

!         Trilinearly interpolate the temperature, extinction, scattering
      TEMP = F1*TEMPP(I1) + F2*TEMPP(I2) + F3*TEMPP(I3) + F4*TEMPP(I4) &
          + F5*TEMPP(I5) + F6*TEMPP(I6) + F7*TEMPP(I7) + F8*TEMPP(I8)
      EXTINCT = F1*EXTINCTP(I1) + F2*EXTINCTP(I2) + F3*EXTINCTP(I3) &
              + F4*EXTINCTP(I4) + F5*EXTINCTP(I5) + F6*EXTINCTP(I6) &
              + F7*EXTINCTP(I7) + F8*EXTINCTP(I8)
      SCAT1 = F1*EXTINCTP(I1)*ALBEDOP(I1)
      SCAT2 = F2*EXTINCTP(I2)*ALBEDOP(I2)
      SCAT3 = F3*EXTINCTP(I3)*ALBEDOP(I3)
      SCAT4 = F4*EXTINCTP(I4)*ALBEDOP(I4)
      SCAT5 = F5*EXTINCTP(I5)*ALBEDOP(I5)
      SCAT6 = F6*EXTINCTP(I6)*ALBEDOP(I6)
      SCAT7 = F7*EXTINCTP(I7)*ALBEDOP(I7)
      SCAT8 = F8*EXTINCTP(I8)*ALBEDOP(I8)
      SCATTER = SCAT1+SCAT2+SCAT3+SCAT4+SCAT5+SCAT6+SCAT7+SCAT8
      IF (EXTINCT .GT. EXTMIN) THEN
        ALBEDO = SCATTER/EXTINCT
      ELSE
        ALBEDO = SCATTER/EXTMIN
      ENDIF
            
!         For tabulated phase functions pick the one we are on top of
!         or the one with the most scattering weight.
      IF (NUMPHASE .GT. 0) THEN
        MAXSCAT = -1.0
        IF (SCAT1 .GT. MAXSCAT .OR. ABS(F1-1) .LT. 0.001) THEN
          MAXSCAT = SCAT1
          IPHASE = IPHASEP(I1)
        ENDIF
        IF (SCAT2 .GT. MAXSCAT .OR. ABS(F2-1) .LT. 0.001) THEN
          MAXSCAT = SCAT2
          IPHASE = IPHASEP(I2)
        ENDIF
        IF (SCAT3 .GT. MAXSCAT .OR. ABS(F3-1) .LT. 0.001) THEN
          MAXSCAT = SCAT3
          IPHASE = IPHASEP(I3)
        ENDIF
        IF (SCAT4 .GT. MAXSCAT .OR. ABS(F4-1) .LT. 0.001) THEN
          MAXSCAT = SCAT4
          IPHASE = IPHASEP(I4)
        ENDIF
        IF (SCAT5 .GT. MAXSCAT .OR. ABS(F5-1) .LT. 0.001) THEN
          MAXSCAT = SCAT5
          IPHASE = IPHASEP(I5)
        ENDIF
        IF (SCAT6 .GT. MAXSCAT .OR. ABS(F6-1) .LT. 0.001) THEN
          MAXSCAT = SCAT6
          IPHASE = IPHASEP(I6)
        ENDIF
        IF (SCAT7 .GT. MAXSCAT .OR. ABS(F7-1) .LT. 0.001) THEN
          MAXSCAT = SCAT7
          IPHASE = IPHASEP(I7)
        ENDIF
        IF (SCAT8 .GT. MAXSCAT .OR. ABS(F8-1) .LT. 0.001) THEN
          MAXSCAT = SCAT8
          IPHASE = IPHASEP(I8)
        ENDIF
      ELSE
!        Standard property file format: unpolarized
        LEGEN(1,0,1) = 1.0
        DO L = 1, NLEG
          LEGEN(1,L,1) = SCAT1*LEGENP(L+(NLEG+1)*(I1-1)) &
                    + SCAT2*LEGENP(L+(NLEG+1)*(I2-1)) &
                    + SCAT3*LEGENP(L+(NLEG+1)*(I3-1)) &
                    + SCAT4*LEGENP(L+(NLEG+1)*(I4-1)) &
                    + SCAT5*LEGENP(L+(NLEG+1)*(I5-1)) &
                    + SCAT6*LEGENP(L+(NLEG+1)*(I6-1)) &
                    + SCAT7*LEGENP(L+(NLEG+1)*(I7-1)) &
                    + SCAT8*LEGENP(L+(NLEG+1)*(I8-1))
          IF (SCATTER .GT. SCATMIN) THEN
            LEGEN(1,L,1) = LEGEN(1,L,1)/SCATTER
          ELSE
            LEGEN(1,L,1) = LEGEN(1,L,1)/SCATMIN
          ENDIF
          LEGEN(1,L,1) = LEGEN(1,L,1)/(2*L+1)
        ENDDO
      ENDIF

!         Add in the gaseous absorption to extinction and albedo
      IF (NZCKD .GT. 0) THEN
        IL = 1
        IU = NZCKD
        DO WHILE (IU-IL .GT. 1)
          IM=(IU+IL)/2
          IF (Z .LE. ZCKD(IM)) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
        I = MIN(MAX(IL,1),NZCKD-1)
        F = (Z-ZCKD(I))/(ZCKD(I+1)-ZCKD(I))
        F = MIN( MAX(F,0.0D0), 1.0D0)
        KG = (1.0-F)*GASABS(I) + F*GASABS(I+1)
        IF (EXTINCT+KG .GT. 0.0) THEN
          ALBEDO = ALBEDO*EXTINCT /(EXTINCT + KG)
        ELSE
          ALBEDO = 0.0
        ENDIF
        EXTINCT = EXTINCT + KG
      ENDIF
      RETURN
      END
 
 


 
      SUBROUTINE DIRECT_BEAM_PROP (INIT, XI, YI, ZI, BCFLAG, IPFLAG, &
                     DELTAM, ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, & 
                     DIRFLUX,  UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM, &
                     NPX, NPY, NPZ, NUMPHASE, DELX, DELY, &
                     XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP, &
                     ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD, &
                     ZCKD, GASABS, CX, CY, CZ, CXINV, CYINV, &
                     CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD, &
                     XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV, NPART,&
                     NPTS)
!       Computes the direct beam flux at point (XI,YI,ZI) by integrating
!     the extinction through the property grid.  If called with 
!     INIT=1 then the property grid extinction array, solar direction
!     terms, and lowest uniform level are computed and saved.  If called
!     with INIT=2 then the input lowest uniform level (UNIFZLEV) is stored.
!     Call with INIT=0 to do the path integration and return the direct
!     beam flux on the horizontal (DIRFLUX) at XI,YI,ZI.
!     The DELTAM flag is used to determine whether to delta-M scale the 
!     extinction for the direct beam computation.  The extinction includes 
!     the gaseous absorption.  If the IPFLAG has bit 2 set then the direct
!     beam tracing is done in 3D, otherwise the lower two bits determine
!     the type of tracing: 0 for 3D tracing, 1 for XZ only tracing, 
!     2 for YZ only tracing, and 3 for Z only tracing.  If BCFLAG bits 0 
!     or 1 are set then have open boundary conditions in X and/or Y.  In 
!     this case when the ray tracing back to the sun reaches the boundary 
!     then independent pixel mode is entered so that only Z grid 
!     intersections occur.
!       For use with multiple processors (bit 2 or 3 set in BCFLAG), 
!     the VALIDBEAM flag is returned true if the ray to the sun made it 
!     to the top of the domain before hitting the subdomain side.  If the
!     flag is false SIDE is returned with the boundary hit (1=-X, 2=+X,
!     3=-Y, 4=+Y).  XE,YE,ZE returns the location of the exitting ray,
!     and path is the optical path from XI,YI,ZI to the sun.
      IMPLICIT NONE
      INTEGER INIT, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, SIDE, NPTS
      LOGICAL DELTAM, VALIDBEAM
      REAL    XI, YI, ZI, SOLARFLUX, SOLARMU, SOLARAZ
      REAL    DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH

      INTEGER IX, IY, IZ, JZ, IL, IM, IU
      INTEGER I, J, K, L, IPH, IP, JP, I1, I2, I3, I4, IPA
      INTEGER IPDIRECT, DI, DJ, DK
      LOGICAL CONSTX, CONSTY, HITBOUNDARY, BTEST
      DOUBLE PRECISION EXTINCT, ALBEDO, F
      DOUBLE PRECISION SUNMU, SUNAZ, PATH
      DOUBLE PRECISION EXTBEAMCUT, UNIFORMZLEV
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
      REAL GASEXT(NPZ), EXTMIN(NPZ), EXTMAX(NPZ)
      
      INTEGER NPX, NPY, NPZ
      INTEGER NUMPHASE, NPART
      REAL DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)
      REAL TEMPP(*), EXTINCTP(NPTS,NPART), ALBEDOP(NPTS,NPART)
      REAL LEGENP(*), EXTDIRP(*)
      INTEGER IPHASEP(NPTS,NPART)
      INTEGER NZCKD
      REAL ZCKD(*), GASABS(*)


      IF (INIT .EQ. 9) THEN
        RETURN
      ENDIF

      IF (INIT .EQ. 1) THEN
!           Get the gaseous extinction at the property grid levels
        DO IZ = 1, NPZ
          IF (NZCKD .GT. 0) THEN
            IL = 1
            IU = NZCKD
            DO WHILE (IU-IL .GT. 1)
              IM=(IU+IL)/2
              IF (ZLEVELS(IZ) .LE. ZCKD(IM)) THEN
                IL = IM
              ELSE
                IU = IM
              ENDIF
            ENDDO
            I = MIN(MAX(IL,1),NZCKD-1)
            W0 = (ZLEVELS(IZ)-ZCKD(I))/(ZCKD(I+1)-ZCKD(I))
            W0 = MIN( MAX(W0,0.0D0), 1.0D0)
            GASEXT(IZ) = (1.0-W0)*GASABS(I) + W0*GASABS(I+1)
          ELSE
            GASEXT(IZ) = 0.0
          ENDIF
        ENDDO
!           First make the property grid extinction field, on which all
!           the direct beam paths will be computed.
        IP = 0
        DO IX = 1, NPX
          DO IY = 1, NPY
            DO IZ = 1, NPZ
              IP = IP + 1
              EXTDIRP(IP) = 0.0
              DO IPA = 1, NPART
                EXTINCT = EXTINCTP(IP,IPA)
                ALBEDO = ALBEDOP(IP,IPA)

!                 Add in the gaseous absorption to extinction
                IF (GASEXT(IZ) .GT. 0.0) THEN
                  ALBEDO = ALBEDO*EXTINCT/(EXTINCT + GASEXT(IZ))
                  EXTINCT = EXTINCT + GASEXT(IZ)
                ENDIF
!                 Do the Delta-M scaling if needed
                IF (DELTAM) THEN
                  IF (NUMPHASE .GT. 0) THEN
                    IPH = IPHASEP(IP,IPA)
                  ELSE
                    IPH = IP
                  ENDIF
                  L = ML+1
                  F = LEGENP(1+NSTLEG*(L+(NLEG+1)*(IPH-1)))/(2*L+1)
                  EXTINCT = (1.0-ALBEDO*F)*EXTINCT
                ENDIF
                EXTDIRP(IP) = EXTDIRP(IP) + EXTINCT
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!           Bit 2 of IPFLAG means do the direct beam in 3D
        IPDIRECT = IPFLAG
        IF (BTEST(IPFLAG,2)) IPDIRECT = 0

!           Make the ray direction (opposite to the solar beam)
!           SOLAR?? is direction of beam, SUN?? is direction to sun (opposite)
        SUNMU = -SOLARMU
        SUNAZ = SOLARAZ + ACOS(-1.0)
        CX = SQRT(1.0-SUNMU**2)*COS(SUNAZ)
        CY = SQRT(1.0-SUNMU**2)*SIN(SUNAZ)
        CZ = ABS(SUNMU)
        IF (ABS(CX) .GT. 1.0E-6) THEN
          CXINV = 1.0D0/CX
        ELSE
          CX = 0.0
          CXINV = 1.0E20
        ENDIF
        IF (ABS(CY) .GT. 1.0E-6) THEN
          CYINV = 1.0D0/CY
        ELSE
          CY = 0.0
          CYINV = 1.0E20
        ENDIF
        IF (ABS(CZ) .GT. 1.0E-6) THEN
          CZINV = 1.0D0/CZ
        ELSE
          CZ = 0.0
          CZINV = 1.0E20
        ENDIF
        DI = NINT(SIGN(1.0D0,CX))
        DJ = NINT(SIGN(1.0D0,CY))
        DK = NINT(SIGN(1.0D0,CZ))
        EPSZ = 1.0E-6*(ZLEVELS(NPZ)-ZLEVELS(1))
        EPSS = 1.0E-3*(ZLEVELS(NPZ)-ZLEVELS(1))/NPZ
        IF (.NOT. BTEST(IPDIRECT,0))  EPSS = MAX(EPSS,1.0D-4*DELX)
        IF (.NOT. BTEST(IPDIRECT,1))  EPSS = MAX(EPSS,1.0D-4*DELY)
        DELXD = DBLE(DELX)
        DELYD = DBLE(DELY)
        EPSS = MAX(0.001*DELXD,0.001*DELYD,EPSS)
        XDOMAIN = DELXD*NPX
        IF (BTEST(BCFLAG,2)) XDOMAIN = DELXD*(NPX-1)
        YDOMAIN = DELYD*NPY
        IF (BTEST(BCFLAG,3)) YDOMAIN = DELYD*(NPY-1)

!           Find the Z level above which the medium is plane-parallel
        EXTBEAMCUT = 1.0E-4
        DO IZ = 1, NPZ
          EXTMIN(IZ) = 1.0E20
          EXTMAX(IZ) = 0.0
        ENDDO
        IP = 0
        DO IX = 1, NPX
          DO IY = 1, NPY
            DO IZ = 1, NPZ
              IP = IP + 1
              EXTMIN(IZ) = MIN(SUM(EXTINCTP(IP,:)), EXTMIN(IZ))
              EXTMAX(IZ) = MAX(SUM(EXTINCTP(IP,:)), EXTMAX(IZ))
            ENDDO
          ENDDO
        ENDDO
        JZ = 0
        DO IZ = 1, NPZ
          IF (EXTMAX(IZ)-EXTMIN(IZ) .GT. EXTBEAMCUT)  JZ = IZ
        ENDDO
        JZ = MIN(NPZ,JZ+1)
        UNIFORMZLEV = ZLEVELS(JZ)
        UNIFZLEV = UNIFORMZLEV
        RETURN
!           Done with initialization
      ENDIF

      IF (INIT .EQ. 2) THEN
        UNIFORMZLEV = UNIFZLEV
        RETURN
      ENDIF



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
      IF (BTEST(BCFLAG,0) .AND. (ABS(X) .LT. 0.01*DELXD  &
         .OR. ABS(X-(NPX-1)*DELXD) .LT. 0.01*DELXD))  CONSTX = .TRUE.
      IF (BTEST(BCFLAG,1) .AND. (ABS(Y) .LT. 0.01*DELYD  &
         .OR. ABS(Y-(NPY-1)*DELYD) .LT. 0.01*DELYD))  CONSTY = .TRUE.

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
        IF (I .LT. 1 .OR. I .GT. NPX .OR. &
            J .LT. 1 .OR. J .GT. NPY .OR. &
            K .LT. 1 .OR. K .GE. NPZ) THEN
          WRITE (6,'(A,3I4)') 'DIRECT_BEAM_PROP: beyond grid!', I, J, K
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
          WRITE (6,*) 'DIRECT_BEAM_PROP: SO<0', X,Y,Z, &
              XE,YE,ZE, XP,YP,ZP, CX,CY,CZ, SOX,SOY,SOZ
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
        A =  (E1*U0M + E2*U0)*VMWM &
           + (E3*U0M + E4*U0)*VWM &
           + (E5*U0M + E6*U0)*VMW &
           + (E7*U0M + E8*U0)*VW
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
      
