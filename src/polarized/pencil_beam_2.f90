      SUBROUTINE DIRECT_BEAM_PROP (INIT, XI, YI, ZI, BCFLAG, IPFLAG, &
                     DELTAM, ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, &
                     DIRFLUX,  UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM, &
                     NPX, NPY, NPZ, NUMPHASE, DELX, DELY, &
                     XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP, &
                     ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD, &
                     ZCKD, GASABS, CX, CY, CZ, CXINV, CYINV, &
                     CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD, &
                     XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV, NPART,&
                     NBPTS)

       IMPLICIT NONE
       INTEGER INIT, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, SIDE, NBPTS
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
       REAL TEMPP(*), EXTINCTP(NBPTS,NPART), ALBEDOP(NBPTS,NPART)
       REAL LEGENP(*), EXTDIRP(*)
       INTEGER IPHASEP(NBPTS,NPART)
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
