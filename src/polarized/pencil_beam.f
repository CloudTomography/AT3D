      SUBROUTINE PENCIL_BEAM_PROP (X0, Y0, Z0, BCFLAG, IPFLAG,
     .                   SOLARFLUX, SOLARMU, SOLARAZ,
     .                   DIRFLUX, TOTAL_EXT,
     .                   NX, NY, NZ, NCELLS, NPTS, CELLFLAGS,
     .                   XGRID, YGRID, ZGRID, GRIDPOS,GRIDPTR,
     .                   NEIGHPTR, TREEPTR)

C    Calculates the source function for a pencil-beam source.
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
      PRINT *, 'x,y,z', X0,Y0,Z0
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

      SUBROUTINE LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, X0, Y0, Z0,  ICELL)
C       Locates the grid cell in the tree structure containing the
C     specified point (X0,Y0,Z0), and returns the cell pointer ICELL.
C     First the base grid cell is found (using NX,NY,NZ, XGRID,YGRID,ZGRID),
C     and the tree is traced down to find the smallest cell containing
C     the point.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NCELLS, ICELL, BCFLAG, IPFLAG
      INTEGER GRIDPTR(8,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(NCELLS)
      DOUBLE PRECISION X0, Y0, Z0
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)
      REAL    GRIDPOS(3,*)
      INTEGER IL, IU, IM, IX, IY, IZ, NXC, NYC, IC, IPTR, DIR, IBITS
      LOGICAL BTEST
      DOUBLE PRECISION XDOMAIN, YDOMAIN

C           Move point to within horizontal domain if periodic
      IF (.NOT. (BTEST(BCFLAG,0) .OR. BTEST(BCFLAG,2))) THEN
        XDOMAIN = XGRID(NX+1)-XGRID(1)
        IF (X0 .LT. XGRID(1)) THEN
          X0 = X0 - XDOMAIN*(INT((X0-XGRID(1))/XDOMAIN)-1)
        ELSE IF (X0 .GT. XGRID(NX+1)) THEN
          X0 = X0 - XDOMAIN*INT((X0-XGRID(1))/XDOMAIN)
        ENDIF
      ENDIF
      IF (.NOT. (BTEST(BCFLAG,1) .OR. BTEST(BCFLAG,3))) THEN
        YDOMAIN = YGRID(NY+1)-YGRID(1)
        IF (Y0 .LT. YGRID(1)) THEN
          Y0 = Y0 - YDOMAIN*(INT((Y0-YGRID(1))/YDOMAIN)-1)
        ELSE IF (Y0 .GT. YGRID(NY+1)) THEN
          Y0 = Y0 - YDOMAIN*INT((Y0-YGRID(1))/YDOMAIN)
        ENDIF
      ENDIF

C         Find the base grid cell (base grid does not have to be even spaced,
C           so do binary search in X, Y, and Z)
C           Special case for independent pixel: get closest grid column.
      IL=0
      IF (BTEST(IPFLAG,0)) THEN
        IU=NX
        DO WHILE (IU-IL .GT. 1)
          IM = (IU+IL)/2
          IF (X0 .GE. 0.5*(XGRID(IM)+XGRID(IM+1))) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
        IL = IL + 1
      ELSE
        IU=NX+1
        DO WHILE (IU-IL .GT. 1)
          IM = (IU+IL)/2
          IF (X0 .GE. XGRID(IM)) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
      ENDIF
      IX = MAX(IL,1)

      IL=0
      IF (BTEST(IPFLAG,1)) THEN
        IU=NY
        DO WHILE (IU-IL .GT. 1)
          IM = (IU+IL)/2
          IF (Y0 .GE. 0.5*(YGRID(IM)+YGRID(IM+1))) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
        IL = IL + 1
      ELSE
        IU=NY+1
        DO WHILE (IU-IL .GT. 1)
          IM = (IU+IL)/2
          IF (Y0 .GE. YGRID(IM)) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
      ENDIF
      IY = MAX(IL,1)

      IL=0
      IU=NZ
      DO WHILE (IU-IL .GT. 1)
        IM = (IU+IL)/2
        IF (Z0 .GE. ZGRID(IM)) THEN
          IL = IM
        ELSE
          IU=IM
        ENDIF
      ENDDO
      IZ = MAX(IL,1)

C         Special case for open boundary conditions: if outside the 1 to NX,NY
C           domain then use the IP grid cells at either end
      NXC = NX
      IF (BTEST(BCFLAG,0)) THEN
        NXC = NX+1
        IF (X0 .LT. XGRID(1)) THEN
C          X0 = XGRID(1)
          IX = 1
        ELSE IF (X0 .GT. XGRID(NX)) THEN
C          X0 = XGRID(NX)
          IX = NX + 1
        ELSE
          IX = IX + 1
        ENDIF
      ENDIF
      IF (BTEST(BCFLAG,2) .AND. .NOT. BTEST(IPFLAG,0)) THEN
        NXC=NX-1
        IX = MIN(IX,NXC)
      ENDIF
      NYC = NY
      IF (BTEST(BCFLAG,1)) THEN
        NYC = NY+1
        IF (Y0 .LT. YGRID(1)) THEN
C          Y0 = YGRID(1)
          IY = 1
        ELSE IF (Y0 .GT. YGRID(NY)) THEN
C          Y0 = YGRID(NY)
          IY = NY + 1
        ELSE
          IY = IY + 1
        ENDIF
      ENDIF
      IF (BTEST(BCFLAG,3) .AND. .NOT. BTEST(IPFLAG,1)) THEN
        NYC=NY-1
        IY = MIN(IY,NYC)
      ENDIF

C         Get the base grid cell pointer
      ICELL = IZ + (NZ-1)*(IY-1) + (NZ-1)*NYC*(IX-1)

C         Trace down the tree: point to the positive child cell, get the
C           grid point that is on the negative (X,Y,Z) side, use the
C           flags to find which dimension (X,Y,Z) the parent cell was
C           split in, then compare the grid point with the test point
C           (X0,Y0,Z0) to find which child cell (pos or neg) the test
C           point is in.
      DO WHILE (TREEPTR(2,ICELL) .GT. 0)
        DIR = IBITS(INT(CELLFLAGS(ICELL)),2,2)
        IC = TREEPTR(2,ICELL) + 1
        IPTR = GRIDPTR(1,IC)
        IF (DIR .EQ. 1) THEN
          IF (X0 .LT. GRIDPOS(1,IPTR))  IC = IC - 1
        ELSE IF (DIR .EQ. 2) THEN
          IF (Y0 .LT. GRIDPOS(2,IPTR))  IC = IC - 1
        ELSE IF (DIR .EQ. 3) THEN
          IF (Z0 .LT. GRIDPOS(3,IPTR))  IC = IC - 1
        ENDIF
        ICELL = IC
      ENDDO

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
