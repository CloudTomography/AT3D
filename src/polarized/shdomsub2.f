
      SUBROUTINE NEW_GRIDS (BCFLAG, GRIDTYPE, NPX,NPY,NPZ, NX, NY, NZ,
     .                      XSTART, YSTART, DELXP, DELYP, ZLEVELS, 
     .                      XGRID, YGRID, ZGRID)
C       Makes the XGRID, YGRID, and ZGRID arrays from the input Z levels
C     and X and Y spacing from the property file according to the GRIDTYPE.
      IMPLICIT NONE
      INTEGER BCFLAG, NPX, NPY, NPZ, NX, NY, NZ
Cf2py intent(in) :: BCFLAG, NPX, NPY, NPZ, NX, NY, NZ
      REAL    XSTART, YSTART, DELXP, DELYP, ZLEVELS(NPZ)
Cf2py intent(in) :: XSTART, YSTART, DELXP, DELYP, ZLEVELS
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)
Cf2py intent(out) :: XGRID, YGRID, ZGRID
      CHARACTER GRIDTYPE*1
Cf2py intent(in) :: GRIDTYPE
      INTEGER IX, IY, IZ

C         Make the X, Y, and Z grids, depending on type:
C           E for evenly spaced, P from property file, F from zgrid.inp
C           Only Even allowed horizontally
      IF (DELXP .LE. 0.0) DELXP = 1.0
      IF (BTEST(BCFLAG,2)) THEN
        DO IX = 1, NX
          XGRID(IX) = XSTART + (IX-1)*(DELXP*(NPX-1))/FLOAT(NX-1)
        ENDDO
      ELSE
        DO IX = 1, NX+1
          XGRID(IX) = XSTART + (IX-1)*(DELXP*NPX)/FLOAT(NX)
        ENDDO
      ENDIF

      IF (DELYP .LE. 0.0) DELYP = 1.0
      IF (BTEST(BCFLAG,3)) THEN
        DO IY = 1, NY
          YGRID(IY) = YSTART + (IY-1)*(DELYP*(NPY-1))/FLOAT(NY-1)
        ENDDO
      ELSE
        DO IY = 1, NY+1
          YGRID(IY) = YSTART + (IY-1)*(DELYP*NPY)/FLOAT(NY)
        ENDDO
      ENDIF

      IF (GRIDTYPE .EQ. 'P') THEN
        IF (NZ .NE. NPZ) THEN
          WRITE (6,*) 
     .    'NEW_GRIDS: For internal Z grid same as property file grid,'
          WRITE (6,*) 'must have same number of grid points (NZ=NPZ).'
          STOP
        ENDIF
        DO IZ = 1, NZ
          ZGRID(IZ) = ZLEVELS(IZ)
        ENDDO
      ELSE IF (GRIDTYPE .EQ. 'E') THEN
        DO IZ = 1, NZ
          ZGRID(IZ) = ZLEVELS(1)
     .             + (IZ-1)*(ZLEVELS(NPZ)-ZLEVELS(1))/FLOAT(NZ-1)
        ENDDO
      ELSE IF (GRIDTYPE .EQ. 'F') THEN
        OPEN (UNIT=3, FILE='zgrid.inp', STATUS='OLD')
        DO IZ = 1, NZ
          READ (3,*) ZGRID(IZ)
        ENDDO
        CLOSE (3)
      ELSE
        STOP 'NEW_GRIDS: Illegal grid type in Z'
      ENDIF

      RETURN
      END





      SUBROUTINE INIT_CELL_STRUCTURE (BCFLAG, IPFLAG, 
     .                     NX, NY, NZ, NX1, NY1, NPTS, NCELLS,  
     .                     XGRID, YGRID, ZGRID, GRIDPOS,
     .                     GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
C     Make the initial grid cell data structure for the base grid.
C     For periodic boundary conditions there are NX*NY*(NZ-1) grid cells 
C     and (NX+1)*(NY+1)*NZ grid points.  The IX=NX+1 column is the same 
C     as the IX=1 column, because of the periodicity, but the grid points 
C     are repeated.  The last grid cell on the "right" has the last column 
C     of grid points on the "left" and the repeated "virtual" set of grid 
C     points on the "right".  For open boundary conditions in X and Y 
C     there are (NX+1)*(NY+1)*(NZ-1) grid cells and NX*NY*NZ grid points.
C     Open boundary conditions are done with the lowest "left most" and
C     "right most" X and Y columns being independent pixels with
C     their own grid cells (zero width), therefore one extra grid cell
C     is required.   Multiple processor boundary conditions have 
C     (NX-1)*(NY-1)*(NZ-1) grid cells and NX*NY*NZ grid points.
C     Independent pixel mode has NX*NY*(NZ-1) grid cells and NX*NY*NZ 
C     grid points
C       Variable       Description
C      IPFLAG         Independent pixel flag - 0 for none, 1 for X, 2 for Y,
C                       3 for X&Y.  
C      BCFLAG         Boundary condition flag - 0 for periodic, 1 for open in X,
C                       2 for open in Y, 3 for open in X&Y.  4 for multiple 
C                       processor in X, 8 for MP in Y, 12 for MP in both.
C      NPTS           total number of grid points
C      NCELLS         total number of cells
C      GRIDPOS(3,IP)  x,y,z coordinate values of grid point locations
C      GRIDPTR(8,IC)  pointer to 8 grid points belonging to a cell.
C                       Order of pointers is based on binary digits zyx:
C                       IOCT = 1 + BITX + 2*BITY + 4*BITZ, where 
C                       BITc is 0 for negative side and 1 for positive side.
C      NEIGHPTR(6,IC) pointer to 6 neighboring cells of a cell.
C                       Order of pointers are -X,+X,-Y,+Y,-Z,+Z (1-6).
C                       Is 0 if there is no neighbor cell (boundary).
C                       Is positive if there is just one cell across face
C                       of current cell (unique neighbor).  Is negative
C                       (-index) if there is more than one cell across 
C                       the face of current cell.  In this case the
C                       cell pointed to is the parent cell that has a 
C                       has same or larger face adjoining current cell.
C      TREEPTR(2,IC)  pointer to parent and children cells in tree structure.
C                       First is pointer to parent cell, 0 for a base cell.
C                       Second is pointer to first of the two children cells.
C                       This pointer is 0 if there are no children cells.
C                       The children cells are in successive order, with
C                       the negative side first.  See CELLFLAGS for direction
C                       of division of split.
C      CELLFLAGS(IC)  a collection of binary flags for the cells:
C                       bit  For bit set
C                        0  Do independent pixel for X
C                        1  Do independent pixel for Y
C                       2,3 Direction cell is split (0=no split, 1=X, 2=Y, 3=Z)
C        IP is a grid point index, IC is a cell index.
      IMPLICIT NONE
      INTEGER IPFLAG, BCFLAG, NX, NY, NZ, NX1, NY1, NPTS, NCELLS
Cf2py intent(in) :: IPFLAG, BCFLAG, NX, NY, NZ, NX1, NY1
Cf2py intent(out) :: NPTS, NCELLS
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in, out) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in, out) :: CELLFLAGS
      REAL    XGRID(NX1), YGRID(NY1), ZGRID(NZ), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID
Cf2py intent(in, out) :: GRIDPOS
      INTEGER IX, IY, IZ, IX0, IX1, IY0, IY1, I, NXC, NYC
      LOGICAL BTEST

C         First, make the grid point positions
      I = 0
      DO IX = 1, NX1
        DO IY = 1, NY1
          DO IZ = 1, NZ
            I = I + 1
            GRIDPOS(1,I) = XGRID(IX)
            GRIDPOS(2,I) = YGRID(IY)
            GRIDPOS(3,I) = ZGRID(IZ)
          ENDDO
        ENDDO
      ENDDO
      NPTS = I

      NXC = NX
      IF (BTEST(BCFLAG,0)) NXC = NX + 1
      IF (BTEST(BCFLAG,2) .AND. .NOT. BTEST(IPFLAG,0)) NXC = NX - 1
      NYC = NY
      IF (BTEST(BCFLAG,1)) NYC = NY + 1
      IF (BTEST(BCFLAG,3) .AND. .NOT. BTEST(IPFLAG,1)) NYC = NY - 1

C         Then setup the cell structures
      I = 0
      DO IX = 1, NXC
        IF (BTEST(IPFLAG,0)) THEN
          IX0 = IX
          IX1 = IX
        ELSE IF (BTEST(BCFLAG,0)) THEN
          IX0 = MAX(1,IX-1)
          IX1 = MIN(NX,IX)
        ELSE
          IX0 = IX
          IX1 = IX+1
        ENDIF
        DO IY = 1, NYC
          IF (BTEST(IPFLAG,1)) THEN
            IY0 = IY
            IY1 = IY
          ELSE IF (BTEST(BCFLAG,1)) THEN
            IY0 = MAX(1,IY-1)
            IY1 = MIN(NY,IY)
          ELSE
            IY0 = IY
            IY1 = IY+1
          ENDIF
          DO IZ = 1, NZ-1
            I = I + 1
            GRIDPTR(1,I) = IZ   + NZ*(IY0-1) + NZ*NY1*(IX0-1)
            GRIDPTR(2,I) = IZ   + NZ*(IY0-1) + NZ*NY1*(IX1-1)
            GRIDPTR(3,I) = IZ   + NZ*(IY1-1) + NZ*NY1*(IX0-1)
            GRIDPTR(4,I) = IZ   + NZ*(IY1-1) + NZ*NY1*(IX1-1)
            GRIDPTR(5,I) = IZ+1 + NZ*(IY0-1) + NZ*NY1*(IX0-1)
            GRIDPTR(6,I) = IZ+1 + NZ*(IY0-1) + NZ*NY1*(IX1-1)
            GRIDPTR(7,I) = IZ+1 + NZ*(IY1-1) + NZ*NY1*(IX0-1)
            GRIDPTR(8,I) = IZ+1 + NZ*(IY1-1) + NZ*NY1*(IX1-1)
            CELLFLAGS(I) = 0
C               Initialize the neighbor pointer; to self for IP mode or
C                 if open boundary conditions and at either end
            IF (BTEST(IPFLAG,0)) THEN
              NEIGHPTR(1,I) = I
              NEIGHPTR(2,I) = I
              CELLFLAGS(I) = IBSET(INT(CELLFLAGS(I)),0)
            ELSE IF (BTEST(BCFLAG,0)) THEN
              IF (IX .EQ. 1) THEN
                NEIGHPTR(1,I) = I
                NEIGHPTR(2,I) = I + (NZ-1)*NYC
                CELLFLAGS(I) = IBSET(INT(CELLFLAGS(I)),0)

              ELSE IF (IX .EQ. NXC) THEN
                NEIGHPTR(1,I) = I - (NZ-1)*NYC
                NEIGHPTR(2,I) = I
                CELLFLAGS(I) = IBSET(INT(CELLFLAGS(I)),0)
              ELSE
                NEIGHPTR(1,I) = I - (NZ-1)*NYC
                NEIGHPTR(2,I) = I + (NZ-1)*NYC
              ENDIF
            ELSE IF (BTEST(BCFLAG,2)) THEN
              IF (IX .EQ. 1) THEN
                NEIGHPTR(1,I) = 0
                NEIGHPTR(2,I) = I + (NZ-1)*NYC
              ELSE IF (IX .EQ. NXC) THEN
                NEIGHPTR(1,I) = I - (NZ-1)*NYC
                NEIGHPTR(2,I) = 0
              ELSE
                NEIGHPTR(1,I) = I - (NZ-1)*NYC
                NEIGHPTR(2,I) = I + (NZ-1)*NYC
              ENDIF
            ELSE
              IF (IX .EQ. 1) THEN
                NEIGHPTR(1,I) = I + (NXC-1)*(NZ-1)*NYC
              ELSE 
                NEIGHPTR(1,I) = I - (NZ-1)*NYC
              ENDIF
              IF (IX .EQ. NX) THEN
                NEIGHPTR(2,I) = I - (NXC-1)*(NZ-1)*NYC
              ELSE
                NEIGHPTR(2,I) = I + (NZ-1)*NYC
              ENDIF
            ENDIF
            IF (BTEST(IPFLAG,1)) THEN
              NEIGHPTR(3,I) = I
              NEIGHPTR(4,I) = I
              CELLFLAGS(I) = IBSET(INT(CELLFLAGS(I)),1)
            ELSE IF (BTEST(BCFLAG,1)) THEN
              IF (IY .EQ. 1) THEN
                NEIGHPTR(3,I) = I
                NEIGHPTR(4,I) = I + (NZ-1)
                CELLFLAGS(I) = IBSET(INT(CELLFLAGS(I)),1)
              ELSE IF (IY .EQ. NYC) THEN
                NEIGHPTR(3,I) = I - (NZ-1)
                NEIGHPTR(4,I) = I
                CELLFLAGS(I) = IBSET(INT(CELLFLAGS(I)),1)
              ELSE
                NEIGHPTR(3,I) = I - (NZ-1)
                NEIGHPTR(4,I) = I + (NZ-1)
              ENDIF
            ELSE IF (BTEST(BCFLAG,3)) THEN
              IF (IY .EQ. 1) THEN
                NEIGHPTR(3,I) = 0
                NEIGHPTR(4,I) = I + (NZ-1)
              ELSE IF (IY .EQ. NYC) THEN
                NEIGHPTR(3,I) = I - (NZ-1)
                NEIGHPTR(4,I) = 0
              ELSE
                NEIGHPTR(3,I) = I - (NZ-1)
                NEIGHPTR(4,I) = I + (NZ-1)
              ENDIF
            ELSE
              IF (IY .EQ. 1) THEN
                NEIGHPTR(3,I) = I + (NYC-1)*(NZ-1)
              ELSE 
                NEIGHPTR(3,I) = I - (NZ-1)
              ENDIF
              IF (IY .EQ. NY) THEN
                NEIGHPTR(4,I) = I - (NYC-1)*(NZ-1)
              ELSE
                NEIGHPTR(4,I) = I + (NZ-1)
              ENDIF
            ENDIF
            IF (IZ .EQ. 1) THEN
              NEIGHPTR(5,I) = 0
            ELSE 
              NEIGHPTR(5,I) = I - 1
            ENDIF
            IF (IZ .EQ. NZ-1) THEN
              NEIGHPTR(6,I) = 0
            ELSE
              NEIGHPTR(6,I) = I + 1
            ENDIF
C               Base grid cells have no parents or children
            TREEPTR(1,I) = 0
            TREEPTR(2,I) = 0
          ENDDO
        ENDDO
      ENDDO
      NCELLS = I
      RETURN
      END




 
      SUBROUTINE INTERP_GRID (NPTS, NSTLEG, NLEG, GRIDPOS,
     .               TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)
C       Calls TRILIN_INTERP_PROP to interpolate the input arrays from 
C     the property grid to each internal grid point. 
      IMPLICIT NONE
      INTEGER NPTS, NSTLEG, NLEG
      INTEGER IPHASE(*)
      REAL    GRIDPOS(3,NPTS)
      REAL    TEMP(*), EXTINCT(*), ALBEDO(*), LEGEN(NSTLEG,0:NLEG,*)
      INTEGER IP

C         Initialize: transfer the tabulated phase functions
      CALL TRILIN_INTERP_PROP (0.0, 0.0, 0.0, .TRUE., NSTLEG, NLEG, 
     .                         TEMP, EXTINCT, ALBEDO, 
     .                         LEGEN(1,0,1), IPHASE)

C         Trilinearly interpolate from the property grid to the adaptive grid
      DO IP = 1, NPTS
        CALL TRILIN_INTERP_PROP 
     .          (GRIDPOS(1,IP), GRIDPOS(2,IP), GRIDPOS(3,IP), 
     .           .FALSE., NSTLEG, NLEG, TEMP(IP), EXTINCT(IP),
     .            ALBEDO(IP), LEGEN(1,0,IP), IPHASE(IP))
      ENDDO 
      RETURN
      END
 




      SUBROUTINE MAKE_DIRECT (NPTS, BCFLAG, IPFLAG, DELTAM, 
     .                ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, 
     .                GRIDPOS, DIRFLUX)
C       Makes the direct beam solar flux for the internal base grid.
C     DIRFLUX is set to F*exp(-tau_sun).
C     Actually calls DIRECT_BEAM_PROP to do all the hard work.
      IMPLICIT NONE
      INTEGER NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG
      LOGICAL DELTAM
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
      REAL    GRIDPOS(3,NPTS), DIRFLUX(NPTS)
      INTEGER SIDE, IP
      LOGICAL VALIDBEAM
      REAL    UNIFZLEV, XO, YO, ZO, DIR, DIRPATH

      CALL DIRECT_BEAM_PROP (1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,
     .            DELTAM, ML, NSTLEG, NLEG, 
     .            SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),
     .            UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)

      DO IP = 1, NPTS
        DIRPATH = 0.0
        CALL DIRECT_BEAM_PROP 
     .           (0, GRIDPOS(1,IP), GRIDPOS(2,IP), GRIDPOS(3,IP),
     .            BCFLAG, IPFLAG, DELTAM, ML, NSTLEG, NLEG,
     .            SOLARFLUX, SOLARMU, SOLARAZ,   DIRFLUX(IP),
     .            UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
      ENDDO
      RETURN
      END
 

 
 
  
      SUBROUTINE PREPARE_PROP (ML, MM, NSTLEG, NLEG, NPTS, DELTAM, 
     .              NUMPHASE, SRCTYPE, UNITS, WAVENO, WAVELEN, ALBMAX, 
     .              EXTINCT, ALBEDO, LEGEN, TEMP, PLANCK, IPHASE)
C       Prepares the grid arrays for the iterative solution process.
C       If doing Delta-M scaling then the extinction, albedo, and Legendre
C       terms are scaled first; only the 0 to ML LEGEN terms are scaled.
C       Outputs PLANCK with (1-omega)*B(T) for thermal source, where B(T) is
C       the Planck function (meaning depends on UNITS).  
C       TEMP array is unchanged.
      IMPLICIT NONE
      INTEGER ML, MM, NSTLEG, NLEG, NPTS, NUMPHASE
      INTEGER IPHASE(NPTS)
      LOGICAL DELTAM
      REAL  WAVENO(2), WAVELEN, ALBMAX
      REAL  EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG,0:NLEG,*)
      REAL  TEMP(NPTS), PLANCK(NPTS)
      CHARACTER*1 SRCTYPE, UNITS
      INTEGER I, IPH, L
      REAL    F, BB

      ALBMAX = 0.0 

      IF (DELTAM .AND. NUMPHASE .GT. 0) THEN
        DO IPH = 1, NUMPHASE
          F = LEGEN(1,ML+1,IPH)
          DO L = 0, ML
C            Scale the diagonal and off-diagonal phase matrix elements differently
            LEGEN(1,L,IPH) = (LEGEN(1,L,IPH) - F)/(1-F)
            IF (NSTLEG .GT. 1) THEN
              LEGEN(2,L,IPH) = (LEGEN(2,L,IPH) - F)/(1-F)
              LEGEN(3,L,IPH) = (LEGEN(3,L,IPH) - F)/(1-F)
              LEGEN(4,L,IPH) = (LEGEN(4,L,IPH) - F)/(1-F)
              LEGEN(5,L,IPH) = LEGEN(5,L,IPH)/(1-F)
              LEGEN(6,L,IPH) = LEGEN(6,L,IPH)/(1-F)
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      DO I = 1, NPTS
        IF (DELTAM) THEN
          IF (NUMPHASE .GT. 0) THEN
            F = LEGEN(1,ML+1,IPHASE(I))
          ELSE
C             NUMPHASE=0 is standard property format, which is unpolarized
            F = LEGEN(1,ML+1,I)
            DO L = 0, ML
              LEGEN(1,L,I) = (LEGEN(1,L,I) - F)/(1-F)
            ENDDO
          ENDIF
          EXTINCT(I) = (1.0-ALBEDO(I)*F)*EXTINCT(I)
          ALBEDO(I) = (1.0-F)*ALBEDO(I)/(1.0-ALBEDO(I)*F)
        ENDIF
        ALBMAX = MAX(ALBMAX,ALBEDO(I))
        IF (SRCTYPE .NE. 'S') THEN
          CALL PLANCK_FUNCTION (TEMP(I), UNITS,WAVENO,WAVELEN,BB)
          PLANCK(I) = (1.0-ALBEDO(I))*BB
        ENDIF
      ENDDO
      RETURN
      END





      SUBROUTINE INIT_RADIANCE (NSTOKES, NXY, NZ, NSTLEG, NLEG, 
     .             RSHPTR, ZGRID,
     .             EXTINCT, ALBEDO, LEGEN, TEMP, NUMPHASE, IPHASE,
     .             SRCTYPE, SOLARFLUX, SOLARMU, GNDALBEDO, GNDTEMP, 
     .             SKYRAD, UNITS, WAVENO, WAVELEN,  RADIANCE)
C       Initializes radiance field by solving plane-parallel two-stream.
C     Solves the L=1 M=0 SH system by transforming the pentadiagonal
C     system to tridiagonal and calling solver.
C     Does the initialization for the NXY columns of the base grid.
C     Only the I Stokes parameter is initialized, the rest are zeroed.
      IMPLICIT NONE
      INTEGER NSTOKES, NXY, NZ, NSTLEG, NLEG, NUMPHASE, RSHPTR(*)
      INTEGER IPHASE(NZ,NXY)
      REAL    GNDALBEDO, GNDTEMP, SKYRAD, SOLARFLUX, SOLARMU
      REAL    WAVENO(2), WAVELEN
      REAL    ZGRID(NZ)
      REAL    EXTINCT(NZ,NXY), ALBEDO(NZ,NXY), LEGEN(NSTLEG,0:NLEG,*)
      REAL    TEMP(NZ,NXY)
      REAL    RADIANCE(NSTOKES,*)
      CHARACTER SRCTYPE*1, UNITS*1
      INTEGER NLAYER, I, IZ, IR, J, K, L
      LOGICAL DELTAM
      REAL    PI, C0, C1
      REAL    EXT0, EXT1, SCAT0, SCAT1, G0, G1, GNDEMIS
      REAL, ALLOCATABLE :: OPTDEPTHS(:), ALBEDOS(:), ASYMMETRIES(:)
      REAL, ALLOCATABLE :: TEMPS(:), FLUXES(:,:)
      CHARACTER SRCT*1

 
      ALLOCATE (OPTDEPTHS(NZ), ALBEDOS(NZ), ASYMMETRIES(NZ))
      ALLOCATE (TEMPS(NZ), FLUXES(3,NZ))

      NLAYER = NZ-1
      DELTAM=.FALSE.
      SRCT = SRCTYPE
      IF (SRCT .EQ. 'B') SRCT='T'
      PI = ACOS(-1.0)
      C0 = SQRT(1.0/PI)
      C1 = SQRT(3.0/(4*PI))

C         Loop over all the columns in the base grid
      IR = 0 
      RSHPTR(1) = 0
      DO I = 1, NXY
C           Make layer properties for the Eddington routine
        DO IZ = 1, NLAYER
          L = NZ-IZ
          EXT0 = EXTINCT(IZ,I)
          EXT1 = EXTINCT(IZ+1,I)
          SCAT0 = ALBEDO(IZ,I)*EXT0
          SCAT1 = ALBEDO(IZ+1,I)*EXT1
          OPTDEPTHS(L) = (ZGRID(IZ+1)-ZGRID(IZ))* (EXT0+EXT1)/2
          IF (EXT0+EXT1 .GT. 0.0) THEN
            ALBEDOS(L) = (SCAT0+SCAT1)/(EXT0+EXT1)
          ELSE
            ALBEDOS(L) = 0.0
          ENDIF
          IF (NUMPHASE .GT. 0) THEN
            J = IPHASE(IZ,I)
            G0 = LEGEN(1,1,J)
            J = IPHASE(IZ+1,I)
            G1 = LEGEN(1,1,J)
          ELSE
            J = NZ*(I-1)+IZ
            G0 = LEGEN(1,1,J)
            G1 = LEGEN(1,1,J+1)
          ENDIF
          IF (SCAT0+SCAT1 .GT. 0.0) THEN
            ASYMMETRIES(L) = (SCAT0*G0+SCAT1*G1)/(SCAT0+SCAT1)
          ELSE
            ASYMMETRIES(L) = 0.0
          ENDIF
          TEMPS(L+1) = TEMP(IZ,I)
        ENDDO
        TEMPS(1) = TEMP(NZ,I)
        GNDEMIS = 1.0-GNDALBEDO
C           Call the Eddington flux routine
        CALL EDDRTF (NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES, 
     .              TEMPS, DELTAM, SRCTYPE, SOLARFLUX, SOLARMU, 
     .              GNDTEMP, GNDEMIS, SKYRAD, UNITS, WAVENO,WAVELEN,
     .              FLUXES)
C           Convert fluxes to first two moments of spherical harmonics
        DO IZ = 1, NZ
          L = NZ+1-IZ
          RADIANCE(:,IR+1) = 0.0
          RADIANCE(:,IR+2) = 0.0
          RADIANCE(:,IR+3) = 0.0
          RADIANCE(:,IR+4) = 0.0
          RADIANCE(1,IR+1) = C0*(FLUXES(1,L)+FLUXES(2,L))
          RADIANCE(1,IR+3) = C1*(FLUXES(1,L)-FLUXES(2,L))
          IR = IR + 4
          RSHPTR(IZ+NZ*(I-1)+1) = IR
        ENDDO
      ENDDO
      DEALLOCATE (OPTDEPTHS, ALBEDOS, ASYMMETRIES, TEMPS, FLUXES)
      RETURN
      END


 
 
      SUBROUTINE EDDRTF (NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES, 
     .              TEMPS, DELTAM, SRCTYPE, SOLARFLUX, SOLARMU, 
     .              GNDTEMP, GNDEMIS, SKYRAD, UNITS, WAVENO,WAVELEN,
     .              FLUXES)
C
C       EDDRTF computes the layer interface fluxes for a plane-parallel
C     atmosphere with either solar or thermal source of radiation using 
C     the Eddington approximation.  The medium is specified by a number 
C     of homogeneous layers.  For a thermal source the Planck function is
C     linear with optical depth, while for a solar source it is exponential.
C     The temperatures, optical depth, single scattering albedo, and 
C     asymmetry parameter are specified for each layer.  The boundary
C     conditions such as the solar flux, and reflection and/or emission from 
C     ground surface are also specified. Delta Eddington scaling may be 
C     used.  The diffuse Eddington fluxes and the solar direct flux at 
C     each level are returned.
C       The model works by calculating the reflection, transmission, and
C     source terms for each layer from the input properties.  A
C     tri-diagonal matrix solver is then used to compute the diffuse fluxes 
C     at each layer interface from the applied boundary conditions.
C
C     Parameters:
C       Input:
C     NLAYER         integer      Number of homogenous layers
C                                  (layers are specified from the top down)
C     OPTDEPTHS      real array   Optical thickness of layers
C     ALBEDOS        real array   Single scattering albedos
C     ASYMMETRIES    real array   Asymmetry parameters
C     TEMPS          real array   Temperatures (K) at layer interfaces
C                                  (e.g. TEMPS(1) is at top of top layer, 
C                                   TEMPS(2) is at bottom of top layer).
C     DELTAM         logical      True for delta-Eddington scaling
C     SRCTYPE        character    'S' for solar source, 'T' for thermal source
C     SOLARFLUX      real         Incident solar flux on horizonal plane
C     SOLARMU        real         Cosine of the solar zenith angle
C     GNDTEMP        real         Ground temperature (Kelvin)
C     GNDEMIS        real         Ground emissivity (1-albedo)
C     SKYRAD         real         Radiance (for solar) or brightness 
C                                  temperature (for thermal) of isotropic 
C                                  incident radiation from above
C     UNITS          character    Units of flux for thermal source
C                                 'T' for brightness temperature, 
C                                 'R' for W/(m^2 micron), 
C                                 'B' for W/(m^2)  for band integration.
C     WAVENO(2)      real         Wavenumber range (cm^-1) (UNITS='B')
C     WAVELEN        real         Wavelength (micron) (UNITS='R')
C
C       Output:
C     FLUXES         real         Eddington fluxes at layer interfaces.
C                                   FLUXES(1,L) is upwelling diffuse, 
C                                   FLUXES(2,L) is downwelling diffuse,
C                                   FLUXES(3,L) is downwelling direct,
C                                   L=1 is top, L=NUML+1 is bottom
      IMPLICIT NONE
      INTEGER   NLAYER
      LOGICAL   DELTAM
      REAL      TEMPS(*)
      REAL      OPTDEPTHS(*), ALBEDOS(*), ASYMMETRIES(*)
      REAL      GNDTEMP, GNDEMIS, SKYRAD, SOLARFLUX, SOLARMU
      REAL      WAVENO(2), WAVELEN
      CHARACTER*1 SRCTYPE, UNITS
      REAL      FLUXES(3,*)
 
      INTEGER   N, L, I
      DOUBLE PRECISION DELTAU, G, OMEGA, F
      DOUBLE PRECISION LAMBDA, R, T, D, CP, CM, A, B, X1, X2
      DOUBLE PRECISION REFLECT, TRANS, SOURCEP, SOURCEM
      DOUBLE PRECISION RADP1P, RADP1M, RADP2P, RADP2M
      DOUBLE PRECISION PI, MU0, SKYFLUX,GNDFLUX, PLANCK1,PLANCK2,C,TAU
      DOUBLE PRECISION EXLP, EXLM, V, DS, B1, B2, SOLPP, SOLPM
      DOUBLE PRECISION, ALLOCATABLE :: LOWER(:),UPPER(:),DIAG(:),RHS(:)
      REAL      BBRAD
      PARAMETER (PI=3.1415926535)


      N = 2*NLAYER+2
      ALLOCATE (LOWER(N), UPPER(N), DIAG(N), RHS(N))
 
C               Compute the reflection, transmission, and source
C               coefficients for each layer for the diffuse Eddington
C               two stream problem.
      IF (SRCTYPE .EQ. 'T') THEN
        CALL PLANCK_FUNCTION (TEMPS(1),UNITS,WAVENO,WAVELEN,BBRAD)
        PLANCK1 = PI*BBRAD
      ENDIF
      MU0 = ABS(SOLARMU)
      TAU = 0.0
      I = 2
      DO L = 1, NLAYER
        DELTAU = OPTDEPTHS(L)
        IF (DELTAU .LT. 0.0) STOP 'EDDRTF: TAU<0'
C            Special case for zero optical depth
        IF (DELTAU .EQ. 0.0) THEN
          TRANS = 1.0
          REFLECT = 0.0
          SOURCEP = 0.0
          SOURCEM = 0.0
        ELSE
          OMEGA = ALBEDOS(L)
          G = ASYMMETRIES(L)
          IF (DELTAM) THEN
            F = G**2
            DELTAU = (1-OMEGA*F)*DELTAU
            OMEGA = (1-F)*OMEGA/(1-OMEGA*F)
            G = (G-F)/(1-F)
          ENDIF
          R = ( 1.0 - OMEGA*(4.0-3.0*G) )/4.0
          T = ( 7.0 - OMEGA*(4.0+3.0*G) )/4.0
          LAMBDA = SQRT( 3.0*(1.0-OMEGA)*(1.0-OMEGA*G) )
C              Special case for conservative scattering (lambda=0)
          IF (LAMBDA .EQ. 0.0) THEN
            D = 1.0/(1.0+T*DELTAU)
            TRANS = D
            REFLECT = -R*DELTAU*D
          ELSE
            X1 = -R
            X2 = LAMBDA + T
            EXLP = DEXP(MIN(LAMBDA*DELTAU,75.D0))
            EXLM = 1.0/EXLP
            TRANS = 2.*LAMBDA/(X2*EXLP + (LAMBDA-T)*EXLM)
            REFLECT = X1*(EXLP - EXLM) *TRANS /(2.*LAMBDA)
            D = 1.0/(X2**2 *EXLP - X1**2 *EXLM)
          ENDIF

          IF (SRCTYPE .EQ. 'T') THEN
C               Calculate thermal source terms
            CALL PLANCK_FUNCTION(TEMPS(L+1),UNITS,WAVENO,WAVELEN,BBRAD)
            PLANCK2 = PI*BBRAD
            V = 2.0*(PLANCK2-PLANCK1)/(3.0*(1.-OMEGA*G)*DELTAU)
            RADP1P = -V + PLANCK1
            RADP2M =  V + PLANCK2
            RADP2P = -V + PLANCK2
            RADP1M =  V + PLANCK1
            IF (LAMBDA .EQ. 0.0) THEN
              A =  (R*DELTAU*RADP1P - RADP2M) *D
              B = -(R*RADP1P + T*RADP2M) *D
              SOURCEP = (B - T*(A+B*DELTAU))/R + RADP2P
              SOURCEM = A + RADP1M
            ELSE
              CP  =  (X1*EXLM*RADP1P - X2*RADP2M) *D
              CM = (-X2*EXLP*RADP1P + X1*RADP2M) *D
              SOURCEP = X1*CP*EXLP + X2*CM*EXLM + RADP2P
              SOURCEM = X2*CP + X1*CM + RADP1M
            ENDIF
            PLANCK1 = PLANCK2
            FLUXES(3,L) = 0.0
          ELSE
C               Calculate solar source terms
            FLUXES(3,L) = SOLARFLUX*EXP(-TAU/MU0)
            DS = 1.0/(LAMBDA**2-1.0/MU0**2)
            B1 = 0.5*OMEGA*(SOLARFLUX/MU0)*EXP(-TAU/MU0) *DS
            B2 = 0.5*OMEGA*(SOLARFLUX/MU0)*EXP(-(TAU+DELTAU)/MU0) *DS
            SOLPP =  1.0 + 1.5*G*MU0
            SOLPM = -1.0 + 1.5*G*MU0
            RADP1P = ( (T+1.0/MU0)*SOLPP + R*SOLPM )*B1
            RADP2M = ((-T+1.0/MU0)*SOLPM - R*SOLPP )*B2
            RADP2P = ( (T+1.0/MU0)*SOLPP + R*SOLPM )*B2
            RADP1M = ((-T+1.0/MU0)*SOLPM - R*SOLPP )*B1
            IF (LAMBDA .EQ. 0.0) THEN
              A =  (R*DELTAU*RADP1P - RADP2M) *D
              B = -(R*RADP1P + T*RADP2M) *D
              SOURCEP = (B - T*(A+B*DELTAU))/R + RADP2P
              SOURCEM = A + RADP1M
            ELSE
              CP  =  (X1*EXLM*RADP1P - X2*RADP2M) *D
              CM = (-X2*EXLP*RADP1P + X1*RADP2M) *D
              SOURCEP = X1*CP*EXLP + X2*CM*EXLM + RADP2P
              SOURCEM = X2*CP + X1*CM + RADP1M
            ENDIF
            TAU = TAU + DELTAU
          ENDIF
        ENDIF
        DIAG(I) = -REFLECT
        DIAG(I+1) = -REFLECT
        LOWER(I) = 1.0
        LOWER(I+1) = -TRANS
        UPPER(I) = -TRANS
        UPPER(I+1) = 1.0
        RHS(I) = SOURCEM
        RHS(I+1) = SOURCEP
        I = I + 2
      ENDDO

C           Set up boundary radiances
      IF (SRCTYPE .EQ. 'S') THEN
        FLUXES(3,NLAYER+1) = SOLARFLUX*EXP(-TAU/MU0)
        GNDFLUX = (1.0-GNDEMIS)*SOLARFLUX*EXP(-TAU/MU0)
        SKYFLUX = PI*SKYRAD
      ELSE
        FLUXES(3,NLAYER+1) = 0.0
        CALL PLANCK_FUNCTION (GNDTEMP,UNITS,WAVENO,WAVELEN,BBRAD)
        GNDFLUX = PI*BBRAD*GNDEMIS
        CALL PLANCK_FUNCTION (SKYRAD,UNITS,WAVENO,WAVELEN,BBRAD)
        SKYFLUX = PI*BBRAD
      ENDIF
C           Setup for and call the tri-diagonal matrix solver
      RHS(1) = SKYFLUX
      DIAG(1) = 0.0
      UPPER(1) = 1.0
      DIAG(N) = -(1.0-GNDEMIS)
      LOWER(N) = 1.0
      RHS(N) = GNDFLUX
      CALL TRIDIAG (N, LOWER, DIAG, UPPER, RHS)

C           Put the fluxes in the output array
      IF (UNITS .EQ. 'T') THEN
        C = 1.0/PI
      ELSE
        C = 1.0
      ENDIF
      I = 1
      DO L = 1, NLAYER+1 
        FLUXES(1,L) = C*RHS(I)
        FLUXES(2,L) = C*RHS(I+1)
        I = I + 2
      ENDDO
 
      DEALLOCATE (LOWER, UPPER, DIAG, RHS)
      RETURN
      END
 


      SUBROUTINE TRIDIAG (N, LOWER, DIAG, UPPER, RHS)
C       Computes the solution to a tridiagonal system. 
C       N is order of the matrix.  LOWER(2..N) is the subdiagonal,
C       DIAG(1..N) is the diagonal, and UPPER(1..N-1) is the 
C       superdiagonal.  On input RHS is the right hand side, while
C       on output it is the solution vector.  Everything is destroyed.
C       Hacked from Linpack DGTSL.
      IMPLICIT NONE
      INTEGER N 
      DOUBLE PRECISION LOWER(*), DIAG(*), UPPER(*), RHS(*)
      INTEGER K, KB
      DOUBLE PRECISION T

      IF (N .EQ. 1) THEN
        IF (DIAG(1) .EQ. 0.0) GOTO 990
        RHS(1) = RHS(1)/DIAG(1)
      ENDIF
      LOWER(1) = DIAG(1)
      DIAG(1) = UPPER(1)
      UPPER(1) = 0.0
      UPPER(N) = 0.0
      DO K = 1, N-1
C              Interchange this and next row to the get the largest pivot.
        IF (ABS(LOWER(K+1)) .GE. ABS(LOWER(K))) THEN
          T = LOWER(K+1)
          LOWER(K+1) = LOWER(K)
          LOWER(K) = T
          T = DIAG(K+1)
          DIAG(K+1) = DIAG(K)
          DIAG(K) = T
          T = UPPER(K+1)
          UPPER(K+1) = UPPER(K)
          UPPER(K) = T
          T = RHS(K+1)
          RHS(K+1) = RHS(K)
          RHS(K) = T
        ENDIF
        IF (LOWER(K) .EQ. 0.0) GOTO 990
        T = -LOWER(K+1)/LOWER(K)
        LOWER(K+1) = DIAG(K+1) + T*DIAG(K)
        DIAG(K+1) = UPPER(K+1) + T*UPPER(K)
        UPPER(K+1) = 0.0
        RHS(K+1) = RHS(K+1) + T*RHS(K)
      ENDDO
      IF (LOWER(N) .EQ. 0.0) GOTO 990

C           Back substitute
      RHS(N) = RHS(N)/LOWER(N)
      RHS(N-1) = (RHS(N-1) - DIAG(N-1)*RHS(N))/LOWER(N-1)
      DO KB = 1, N-2
        K = N - 2 - KB + 1
        RHS(K) = (RHS(K) -DIAG(K)*RHS(K+1) -UPPER(K)*RHS(K+2))/LOWER(K)
      ENDDO
      RETURN

990   CONTINUE
        STOP 'Singular matrix in TRIDIAG'
      END
 


 
 

      SUBROUTINE MAKE_ANGLE_SET (NMU, NPHI, ITYPE, NPHI0MAX,
     .             NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)
C       Make the set of angles for the discrete space representation.
C     The number of mu's (cosine zenith angles) and maximum number of
C     phi's (azimuth angles) is input.  The actual number of azimuthal 
C     angles for each mu is output (NPHI0).  There are three types of
C     discrete ordinate sets: ITYPE=1 is a gaussian grid, 2 is a reduced
C     gaussian grid, and 3 is a reduced double gaussian set.
C     The output is the NMU mu values, the NPHI0 phi values for each mu,
C     and the integration weight for each ordinate. The first NMU/2 mu 
C     angles are the downwelling (mu<0) angles.  Also output are the
C     maximum number of azimuthal angles (NPHI0MAX), the flags for doing 
C     an azimuthal FFT for each mu and the total number of angles.
      IMPLICIT NONE
      INTEGER NMU, NPHI, ITYPE,  NPHI0MAX, NPHI0(NMU), NANG
      LOGICAL FFTFLAG(NMU)
      REAL    MU(NMU), PHI(NMU,*), WTMU(NMU), WTDO(NMU,*)
      INTEGER MAXNPHI, J, K, MM
      REAL    DELPHI
      PARAMETER (MAXNPHI=256)
      INTEGER GOODNFFT(MAXNPHI)
      DATA GOODNFFT/1,2,3,4,5,6,7,8,   9,10,11,12,13,14,16,16,
     .      18,18,20,20,24,24,24,24,  25,27,27,30,30,30,32,32,
     .      36,36,36,36,40,40,40,40,  45,45,45,45,45,48,48,48,
     .      50,50,54,54,54,54,60,60,  60,60,60,60,64,64,64,64,
     .      72,72,72,72,72,72,72,72,  80,80,80,80,80,80,80,80,
     .      81,90,90,90,90,90,90,90,  90,90,96,96,96,96,96,96,
     .          100,100,100,100,108,108,108,108, 
     .          108,108,108,108,120,120,120,120,
     .          120,120,120,120,120,120,120,120,
     .          128,128,128,128,128,128,128,128,
     .          144,144,144,144,144,144,144,144,
     .          144,144,144,144,144,144,144,144,
     .          160,160,160,160,160,160,160,160,
     .          160,160,160,160,160,160,160,160,
     .          180,180,180,180,180,180,180,180,180,180,
     .          180,180,180,180,180,180,180,180,180,180,
     .          192,192,192,192,192,192,192,192,192,192,192,192,
     .          200,200,200,200,200,200,200,200,
     .          216,216,216,216,216,216,216,216,
     .          216,216,216,216,216,216,216,216,
     .          240,240,240,240,240,240,240,240,240,240,240,240,
     .          240,240,240,240,240,240,240,240,240,240,240,240,
     .          256,256,256,256,256,256,256,256,
     .          256,256,256,256,256,256,256,256/

      MM = MAX(0,INT(NPHI/2)-1)
      NANG = 0
      IF (ITYPE .GE. 1 .OR. ITYPE .LE. 3) THEN
        IF (ITYPE .LE. 2) THEN
          CALL GAUSQUADS (NMU, MU, WTMU)
        ELSE IF (ITYPE .EQ. 3) THEN
          CALL DGAUSQUADS (NMU, MU, WTMU)
        ENDIF
        DO J = 1, NMU
          NPHI0MAX = NPHI
          IF (ITYPE .EQ. 1) THEN
            NPHI0(J) = NPHI0MAX
          ELSE
            IF (NPHI0MAX .GT. MAXNPHI) 
     .          STOP 'MAKE_ANGLE_SET: exceeded GOODNFFT'
            NPHI0(J) = INT(0.9+NPHI0MAX*SQRT(1-MU(J)**2))
            NPHI0(J) = MIN(NPHI0MAX,GOODNFFT(NPHI0(J)))
          ENDIF
          DELPHI = 2.0*ACOS(-1.0)/NPHI0(J)
          DO K = 1, NPHI0(J)
            PHI(J,K) = (K-1)*DELPHI
            WTDO(J,K) = DELPHI*WTMU(J)
          ENDDO
          NANG = NANG + NPHI0(J)
          FFTFLAG(J) = (NPHI0(J) .GT. 14) .OR. (MM .GT. 15)
        ENDDO
        
      ELSE
        STOP 'MAKE_ANGLE_SET: invalid discrete ordinate type'
      ENDIF
      RETURN
      END




      SUBROUTINE MAKE_SH_DO_COEF (NSTLEG, ML, MM, NLM, NMU, NPHI0,
     .             NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG,
     .             CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)
C       Makes the transformation coefficients for the real generalized
C     spherical harmonic transform.  The main coefficients are output in 
C     four arrays: 1 is for the SH_TO_DO forward transform, 2 is for the
C     DO_TO_SH back transform, which contains the discrete angle integration 
C     weights.  The coefficients are divided up into the mu dependent set
C      CMUn (function of l, m, mu_j), and the phi dependent set CPHIn 
C     (function of m, phi_k) for each mu_j.
C     The FFTPACK phase coefficients for the FFT in azimuth are also 
C     output in WPHISAVE.
      IMPLICIT NONE
      INTEGER NSTLEG, ML, MM, NLM, NMU, NPHI0(NMU), NPHI0MAX
      LOGICAL FFTFLAG(NMU)
      REAL    MU(NMU), PHI(NMU,*), WTMU(NMU), WTDO(NMU,*)
      REAL    CMU1(NSTLEG,NLM,NMU), CMU2(NSTLEG,NMU,NLM)
      REAL    CPHI1(-16:16,32,NMU), CPHI2(32,-16:16,NMU)
      REAL    WPHISAVE(3*NPHI0MAX+15,NMU)
      INTEGER I, J, K, M, Q
      REAL    X, W
      REAL, ALLOCATABLE :: PRC(:,:)

C        Make the mu part of the real generalized spherical harmonics 
      ALLOCATE (PRC(6,NLM))
      DO I = 1, NMU
        X = MU(I)
        PRC(:,:) = 0.0
        CALL PLMALL (.FALSE., X, ML, MM, PRC)
        DO J = 1, NLM
          DO Q = 1, NSTLEG
            CMU1(Q,J,I) = PRC(Q,J)
          END DO
        END DO
        PRC(:,:) = 0.0
        CALL PLMALL (.TRUE., X, ML, MM, PRC)
        DO J = 1, NLM
          DO Q = 1, NSTLEG
            CMU2(Q,I,J) = PRC(Q,J)*WTMU(I)
          END DO
        END DO
      END DO
      DEALLOCATE (PRC)

C         Make the to and from Fourier coefficients for each mu
      DO I = 1, NMU
        IF (.NOT. FFTFLAG(I)) THEN
C             If not doing an FFT for this mu then make the DFT coefficient
          W = 1.0/WTMU(I)
          DO K = 1, NPHI0(I)
            IF (NPHI0(I) .GT. 32 .OR. MM .GT. 16)  STOP 
     .        'MAKE_SH_DO_COEF: Fourier coefficients array exceeded'
            CPHI1(0,K,I) = 1.0
            CPHI2(K,0,I) = WTDO(I,K)*W
            DO M = 1, MM
              CPHI1(M,K,I)  = COS(M*PHI(I,K))
              CPHI1(-M,K,I) = SIN(M*PHI(I,K))
              CPHI2(K,M,I)  = CPHI1(M,K,I)*WTDO(I,K)*W
              CPHI2(K,-M,I) = CPHI1(-M,K,I)*WTDO(I,K)*W
            ENDDO
          ENDDO
        ELSE
C             Precompute the phase factors for the FFTs
          CALL RFFTI (NPHI0(I),WPHISAVE(1,I))
        ENDIF
      ENDDO
      RETURN
      END
 

 


      SUBROUTINE SURFACE_BRDF (SFCTYPE, REFPARMS, WAVELEN,
     .                         MU2, PHI2, MU1, PHI1, NSTOKES, REFLECT)
C       Returns the reflection matrix for the general bidirectional
C     reflection distribution function of the specified type (SFCTYPE).
C     The incident direction is (MU1,PHI1), and the outgoing direction
C     is (MU2,PHI2) (MU is cosine of zenith angle, and PHI is the azimuthal
C     angle in radians).  The incident directions have mu<0 (downward), 
C     while the outgoing directions have mu>0 (upward). The reflection 
C     function is normalized so that for a Lambertian surface (uniform 
C     reflection) the returned value (REFLECT) is simply the albedo.
C       This routine calls the desired BRDF function, passing the 
C     appropriate parameters.  More BRDF surface types may be added easily 
C     by putting in the appropriate function calls.
C            Type        Parameters
C       L  Lambertian    albedo
C       W  WaveFresnel   Real, Imaginary index of refraction, wind speed (m/s)
C       D  Diner et al   a, k, b, zeta, sigma
C       O  Ocean         Wind speed (m/s), Pigmentation (mg/m^3)
C       R  RPV-original  rho0, k, Theta
      IMPLICIT NONE
      INTEGER NSTOKES
      REAL    REFPARMS(*), WAVELEN, MU1, PHI1, MU2, PHI2, REFLECT(4,4)
      CHARACTER  SFCTYPE*1
      INTEGER K1, K2
      REAL   PI
      REAL   RPV_REFLECTION

      PI = ACOS(-1.0)
      IF (SFCTYPE .EQ. 'L' .OR. SFCTYPE .EQ. 'l') THEN
C         L or l: Lambertian surface BRDF is constant.
C           (for testing, as done more efficiently by Lambertian routines)
        REFLECT(1:NSTOKES,1:NSTOKES) = 0.0
        REFLECT(1,1) = REFPARMS(1)
      ELSE IF (SFCTYPE .EQ. 'W') THEN
C         W: Fresnel reflection with Gaussian distribution of slopes
        CALL WAVE_FRESNEL_REFLECTION (REFPARMS(1), REFPARMS(2), 
     .            REFPARMS(3), MU1, MU2, PHI1, PHI2, NSTOKES, REFLECT)
      ELSE IF (SFCTYPE .EQ. 'D') THEN
C         D: Diner et al., depolarizing modified RPV + Fresnel reflection 
C            from microfacets with uniform or Gaussian distribution of slopes
        CALL DINER_REFLECTION (REFPARMS(1), REFPARMS(2), REFPARMS(3),
     .                         REFPARMS(4), REFPARMS(5),
     .                         -MU1, MU2, PHI2-PHI1, NSTOKES, REFLECT)
      ELSE IF (SFCTYPE .EQ. 'R') THEN
C         R: Rahman, Pinty, and Verstraete
        IF (NSTOKES .GT. 1) THEN
          PRINT *, 'RPV BRDF is only for unpolarized case'
          STOP
        ENDIF
        REFLECT(1,1) = RPV_REFLECTION (REFPARMS(1),REFPARMS(2),
     .                           REFPARMS(3), -MU1, MU2, PHI1-PHI2-PI)
      ELSE IF (SFCTYPE .EQ. 'O') THEN
C         O: Ocean BRDF from 6S modified by Norm Loeb
        IF (NSTOKES .GT. 1) THEN
          PRINT *, 'Ocean BRDF from 6S is only for unpolarized case'
          STOP
        ENDIF
        CALL ocean_brdf_sw (REFPARMS(1), -1., REFPARMS(2), WAVELEN,
     .                      -MU1, MU2, PHI1-PHI2, PHI1, REFLECT(1,1))
      ELSE
        STOP 'SURFACE_BRDF: Unknown BRDF type'
      ENDIF

      RETURN
      END


  
   

      SUBROUTINE WAVE_FRESNEL_REFLECTION (MRE, MIM, WINDSPEED, 
     .                         MUI, MUR, PHII, PHIR, NSTOKES, REFLECT)
C       Computes the polarized Fresnel reflection from a dielectric
C     interface with a Gaussian distribution of slopes, including shadowing.
C     The complex index of refraction of the dielectric is specified with 
C     MRE,MIM (MIM is negative for an absorbing medium) and the wind speed 
C     is in m/s.  The incident direction is (MUI,PHII), and the outgoing 
C     direction is (MUR,PHIR) (MU is cosine of zenith angle with MU>0 and 
C     PHI is the azimuthal angle in radians). 
C       The code is courtesy of Michael Mishchenko, with some modifications.
C     The mean square surface slope (s^2) or SIGMA2 in Mishchenko and Travis
C     is obtained from the Cox and Munk formula in the paper.
C       For all formulas and definitions, see the paper:
C     M. I. Mishchenko and L. D. Travis, 1997: Satellite retrieval
C     of aerosol properties over the ocean using polarization as well as
C     intensity of reflected sunlight.  J. Geophys. Res. 102, 16989-17013.
C
      IMPLICIT NONE
      INTEGER NSTOKES
      REAL    MRE, MIM, WINDSPEED, MUI, MUR, PHII, PHIR, REFLECT(4,4)
      INTEGER I, J
      REAL*8  SIGMA2, DMUI, DMUR, DCOSI, DSINI, DCOSR, DSINR, DSI, DSR
      REAL*8  VI1, VI2, VI3, VR1, VR2, VR3, UNIT1, UNIT2, UNIT3
      REAL*8  FACT1, FACTOR, XI1
      REAL*8  TI1, TI2, TI3, TR1, TR2, TR3, PI1, PI2, PI3, PR1,PR2,PR3
      REAL*8  PIKR, PRKI, TIKR, TRKI, E1, E2, E3, E4
      REAL*8  VP1, VP2, VP3, DMOD, RDZ2, RDZ4, DCOEFF, DEX
      REAL*8  AF, AF11, AF12, AF21, AF22
      REAL*8  P, S1, S2, S3, XI, XXI, DCOT, T1, T2
      REAL*8  SHADOWI, SHADOWR, SHADOW, DERFC
      COMPLEX*16 CN1, CN2, CXI2, C1, C2, CRPER, CRPAR
      COMPLEX*16 CF11, CF12, CF21, CF22
      COMPLEX*16 CI, C21, C22, CTTTP, CTTPT, CTTPP, CTPPT, CTPPP, CPTPP

      CN1 = (1.0D0,0.0D0)
      CN2 = CMPLX(MRE,MIM)
      SIGMA2 = MAX(0.0005D0,0.0015D0+0.00256D0*WINDSPEED)
      DMUI = DABS(DBLE(MUI))
      DMUR = DBLE(MUR)
 
C       Cartesian components of the unit vectors of the incident and
C       scattered beams
      IF (DABS(DMUI-1E0) .LT. 1D-10) DMUI=0.999999999999D0
      IF (DABS(DMUR-1E0) .LT. 1D-10) DMUR=0.999999999999D0
      DCOSI=DCOS(DBLE(PHII))
      DSINI=DSIN(DBLE(PHII))
      DCOSR=DCOS(DBLE(PHIR))
      DSINR=DSIN(DBLE(PHIR))
      DSI=DSQRT(1D0-DMUI*DMUI)
      DSR=DSQRT(1D0-DMUR*DMUR)
      VI1=DSI*DCOSI
      VI2=DSI*DSINI
      VI3=-DMUI
      VR1=DSR*DCOSR
      VR2=DSR*DSINR
      VR3=DMUR
C       Local surface normal for specular reflection
      UNIT1=VI1-VR1
      UNIT2=VI2-VR2
      UNIT3=VI3-VR3
      FACT1=UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR=SQRT(1D0/FACT1)
 
C      Fresnel reflection coefficients
      XI1=FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2=CDSQRT(1D0 - (1D0-XI1*XI1)*CN1*CN1/(CN2*CN2))
      C1=CN1*XI1
      C2=CN2*CXI2
      CRPER=(C1-C2)/(C1+C2)
      C1=CN2*XI1
      C2=CN1*CXI2
      CRPAR=(C1-C2)/(C1+C2)
 
C      Calculation of the amplitude scattering matrix
      TI1=-DMUI*DCOSI
      TI2=-DMUI*DSINI
      TI3=-DSI
      TR1=DMUR*DCOSR
      TR2=DMUR*DSINR
      TR3=-DSR
      PI1=-DSINI
      PI2=DCOSI
      PI3=0D0
      PR1=-DSINR
      PR2=DCOSR
      PR3=0D0
      PIKR=PI1*VR1+PI2*VR2+PI3*VR3
      PRKI=PR1*VI1+PR2*VI2+PR3*VI3
      TIKR=TI1*VR1+TI2*VR2+TI3*VR3
      TRKI=TR1*VI1+TR2*VI2+TR3*VI3
      E1=PIKR*PRKI
      E2=TIKR*TRKI
      E3=TIKR*PRKI
      E4=PIKR*TRKI
      CF11=E1*CRPER+E2*CRPAR
      CF12=-E3*CRPER+E4*CRPAR
      CF21=-E4*CRPER+E3*CRPAR
      CF22=E2*CRPER+E1*CRPAR
 
C      Calculation of the Stokes reflection matrix
      VP1=VI2*VR3-VI3*VR2
      VP2=VI3*VR1-VI1*VR3
      VP3=VI1*VR2-VI2*VR1
      DMOD=(VP1*VP1+VP2*VP2+VP3*VP3)**2
      RDZ2=UNIT3*UNIT3
      RDZ4=RDZ2*RDZ2
      DEX=DEXP(-(UNIT1*UNIT1 + UNIT2*UNIT2)/(2*SIGMA2*RDZ2))
      DCOEFF=FACT1**2 *DEX/(4*DMUI*DMUR*DMOD*RDZ4*2*SIGMA2)
      AF=0.5D0*DCOEFF
      AF11=CDABS(CF11)**2
      AF12=CDABS(CF12)**2
      AF21=CDABS(CF21)**2
      AF22=CDABS(CF22)**2
 
      REFLECT(1,1)=(AF11+AF12+AF21+AF22)*AF
      IF (NSTOKES .GE. 2) THEN 
        REFLECT(1,2)=(AF11-AF12+AF21-AF22)*AF
        REFLECT(2,1)=(AF11-AF22+AF12-AF21)*AF
        REFLECT(2,2)=(AF11-AF12-AF21+AF22)*AF
      ENDIF

      CI=(0.0D0,-1.0D0)
      C21=CONJG(CF21)
      C22=CONJG(CF22)
      CTTTP=CF11*CONJG(CF12)
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22
 
      IF (NSTOKES .GE. 3) THEN 
        REFLECT(1,3)=    (-CTTTP-CPTPP)*DCOEFF
        REFLECT(2,3)=    (-CTTTP+CPTPP)*DCOEFF
        REFLECT(3,1)=    (-CTTPT-CTPPP)*DCOEFF
        REFLECT(3,2)=    (-CTTPT+CTPPP)*DCOEFF
        REFLECT(3,3)=    ( CTTPP+CTPPT)*DCOEFF
      ENDIF  
      IF (NSTOKES .GE. 4) THEN 
        REFLECT(1,4)=-CI*( CTTTP+CPTPP)*DCOEFF
        REFLECT(2,4)=-CI*( CTTTP-CPTPP)*DCOEFF
        REFLECT(3,4)= CI*( CTTPP-CTPPT)*DCOEFF
        REFLECT(4,1)= CI*( CTTPT+CTPPP)*DCOEFF
        REFLECT(4,2)= CI*( CTTPT-CTPPP)*DCOEFF
        REFLECT(4,3)=-CI*( CTTPP+CTPPT)*DCOEFF
        REFLECT(4,4)=    ( CTTPP-CTPPT)*DCOEFF
      ENDIF

C      Shadowing
      P=DACOS(-1.0D0)
      S1=DSQRT(2*SIGMA2/P)
      S3=1D0/(DSQRT(2*SIGMA2))
      S2=S3*S3
      DCOT=DMUI/SQRT(1D0-DMUI**2)
      T1=DEXP(-S2*DCOT**2)
      T2=DERFC(DCOT*S3)
      SHADOWI=0.5D0*(S1*T1/DCOT-T2)
      DCOT=DMUR/DSQRT(1D0-DMUR**2)
      T1=DEXP(-S2*DCOT**2)
      T2=DERFC(DCOT*S3)
      SHADOWR=0.5*(S1*T1/DCOT-T2)
      SHADOW=1D0/(1D0+SHADOWI+SHADOWR)
      REFLECT(1:NSTOKES,1:NSTOKES) = REFLECT(1:NSTOKES,1:NSTOKES)*SHADOW
      RETURN
      END   



      SUBROUTINE DINER_REFLECTION (A, K, B, ZETA, SIGMA, 
     .                         MU1, MU2, PHI, NSTOKES, REFLECT)
C       Computes polarized reflection from the Diner et al. (2012) model
C     which two terms: 1) completely depolarizing volumetric scattering with 
C     the modified RPV model, and 2) Fresnel reflection from randomly 
C     oriented microfacets with either a uniform or Gaussian distribution
C     of slopes.  The three parameters of the modified-RPV model are A, K,
C     and B.  The weight given to the Fresnel-microfacet term is ZETA (if
C     ZETA=0 then only the REFLECT(1,1) element is non-zero with the
C     modified-RPV model).  SIGMA is the standard deviation of microfacet
C     slopes for the Gaussian orientation distribution, but if SIGMA<=0
C     the uniform orientation distribution is used.  The complex index of 
C     refraction of the Fresnel microfacets is specified in this subroutine.
C     The Diner et al. (2012) model leaves out the hotspot factor of the
C     modified RPV model (and the original RPV), but the hotspot factor is
C     included here.
C       The cosine of the incident zenith angle is MU1 (MU1>0), the 
C     cosine of the outgoing zenith angle is MU2 (MU2>0), and the relative 
C     azimuthal angle is PHI=PHI2-PHI1 (radians).
C       For all formulas and definitions, see the paper:
C     Diner, D. J., F. Xu, J. V. Martonchik, B. E. Rheingans, S. Geier,
C     V. M. Jovanovic, A. Davis, R. A. Chipman, S. C. McClain, 2012:
C     Exploration of a polarized surface bidirectional reflectance model 
C     using the ground-based multiangle spectropolarimetric imager.
C     Atmosphere 2012, 3, 591-619; doi:10.3390/atmos3040591.
C
      IMPLICIT NONE
      INTEGER NSTOKES
      REAL    A, K, B, ZETA, SIGMA, MU1, MU2, PHI
      REAL    REFLECT(4,4)
      INTEGER K1, K2
      REAL    SINTH1, SINTH2, COSPHI, COSSCATANG, TAN1,TAN2, CAPG, HOT
      REAL    GAMMA, COSGAMMA, F11, F12, F33, F34
      REAL    COSBETA, H, SINPHI, ALPHA1, ALPHA2
      REAL    COS2ALPHA1, SIN2ALPHA1, COS2ALPHA2, SIN2ALPHA2
      COMPLEX SFCINDEX, EPSILON, D, RP, RS

C      Another free parameter of the Fresnel microfacet model is the 
C        index of refraction specified here
      SFCINDEX = CMPLX(1.5,0.0)

      REFLECT(1:NSTOKES,1:NSTOKES) = 0.0

C      The modified RPV model (Martonchik et al., 1998, IEEE-TGARS, 36, 1266)
C      is implemented here.  Equation 6 in Diner et al. leaves out the hotspot factor.
      SINTH1 = SQRT(1.0-MU1**2)
      SINTH2 = SQRT(1.0-MU2**2)
      COSPHI = COS(PHI)
      COSSCATANG = -MU1*MU2 + SINTH1*SINTH2*COSPHI
      COSSCATANG = MIN(1.0,MAX(-1.0,COSSCATANG))
      TAN1 = SINTH1/MU1
      TAN2 = SINTH2/MU2
      CAPG = SQRT(ABS(TAN1**2 + TAN2**2 + 2*TAN1*TAN2*COSPHI))
      HOT = 1 + (1-A)/(1+CAPG)
      REFLECT(1,1) = A*((MU1+MU2)*MU1*MU2)**(K-1) *EXP(B*COSSCATANG)
      REFLECT(1,1) = REFLECT(1,1)*HOT
      IF (ZETA .LT. 0.0) RETURN

C      Calculate incident/reflectance angle relative to microfacet (eq 16a)
      GAMMA = 0.5*ACOS(-COSSCATANG)
      COSGAMMA = COS(GAMMA)

C      Calculate the Fresnel reflection coefficients (equivalent to eq 17)
C        [D = n*cos(gamma')]
      EPSILON = SFCINDEX**2
      D = SQRT(EPSILON - 1.0 + COSGAMMA**2)
      RP =(EPSILON*COSGAMMA - D) / (EPSILON*COSGAMMA + D)
      RS = (COSGAMMA - D) / (COSGAMMA + D)

C      Calculate the Mueller matrix elements for Fresnel reflection (eq 18)
      F11 = 0.5*(CABS(RP)**2 + CABS(RS)**2)
      F12 = 0.5*(CABS(RP)**2 - CABS(RS)**2)
      F33 = REAL(RP*CONJG(RS))
      F34 = -IMAG(RP*CONJG(RS))

C      Calculate facet tilt angle, beta (eq 16b)
      COSBETA = 0.5*(MU1+MU2)/COSGAMMA

C      Calculate the factor before the M*P*M matrix including the distribution
C       of microfacet orientations (Gaussian for SIGMA>0) (eq 13 & 14)
C       Note: Diner's BRDF is multiplied by pi here
      IF (SIGMA .GT. 0.0) THEN
        H = ZETA *EXP(-0.5*(1/COSBETA**2-1)/SIGMA**2) 
     .        /(8*SIGMA**2*MU2*MU1*COSBETA**4)
      ELSE
        H = ZETA/(8*MU2*MU1*COSBETA)
      ENDIF

      REFLECT(1,1) = REFLECT(1,1) + H*F11

      IF (NSTOKES .GE. 2) THEN
C        Calculate the polarization rotation angles (eq 20 & 21)
        SINPHI = SIN(PHI)
        ALPHA1 = ATAN(SINTH2*SINPHI/(MU2*SINTH1+SINTH2*MU1*COSPHI))
        ALPHA2 = ATAN(SINTH1*SINPHI/(SINTH2*MU1+MU2*SINTH1*COSPHI))
        COS2ALPHA1 = COS(2*ALPHA1)
        SIN2ALPHA1 = SIN(2*ALPHA1)
        COS2ALPHA2 = COS(2*ALPHA2)
        SIN2ALPHA2 = SIN(2*ALPHA2)

C        Calculate the rotated Mueller matrix (M*P*M) for Fresnel reflection
C         and multiply by the H factor to get the final BRDF
        REFLECT(1,2) = H*F12*COS2ALPHA1
        REFLECT(2,1) = H*F12*COS2ALPHA2
        REFLECT(2,2) = H*(F11*COS2ALPHA1*COS2ALPHA2 
     .                  + F33*SIN2ALPHA1*SIN2ALPHA2)
      ENDIF
      IF (NSTOKES .GE. 3) THEN
        REFLECT(1,3) = -H*F12*SIN2ALPHA1
        REFLECT(2,3) = H*(-F11*SIN2ALPHA1*COS2ALPHA2
     .                   + F33*COS2ALPHA1*SIN2ALPHA2)
        REFLECT(3,1) = -H*F12*SIN2ALPHA2
        REFLECT(3,2) = H*(-F11*COS2ALPHA1*SIN2ALPHA2
     .                   + F33*SIN2ALPHA1*COS2ALPHA2)
        REFLECT(3,3) =  H*(F11*SIN2ALPHA1*SIN2ALPHA2
     .                   + F33*COS2ALPHA1*COS2ALPHA2)
      ENDIF
      IF (NSTOKES .GE. 4) THEN
        REFLECT(1,4) = 0.0
        REFLECT(2,4) = H*F34*SIN2ALPHA2
        REFLECT(3,4) = H*F34*COS2ALPHA2
        REFLECT(4,1) = 0.0
        REFLECT(4,2) = -H*F34*SIN2ALPHA1
        REFLECT(4,3) = -H*F34*COS2ALPHA1
        REFLECT(4,4) = H*F33
      ENDIF
      RETURN
      END




      REAL FUNCTION RPV_REFLECTION (RHO0, K, THETA, MU1, MU2, PHI)
C       Computes the Rahman, Pinty, Verstraete BRDF.  The incident
C     and outgoing cosine zenith angles are MU1 and MU2, respectively,
C     and the relative azimuthal angle is PHI.  In this case the incident
C     direction is where the radiation is coming from (i.e. opposite of the
C     discrete ordinate), so MU1>0 and the hot spot is MU2=MU1 and PHI=0.
C     The reference is:
C       Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere 
C       Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable 
C       With NOAA Advanced Very High Resolution Radiometer Data,
C       J. Geophys. Res., 98, 20791-20801.
C
C     Implementation modified by Alexei Lyapustin to prevent negative
C     values in SQRT function in calculating CAPG, and to prevent
C     unphysical growth of RPV BRDF at small mu or mu0.

      REAL RHO0, K, THETA, MU1, MU2, PHI, mu_min, x1, x2
      REAL M, F, H, COSPHI, SIN1, SIN2, COSG, TAN1, TAN2, CAPG

      mu_min = 0.03
      x1 = MU1
      if(x1.LT.mu_min) x1 = mu_min
      x2 = MU2
      if(x2.LT.mu_min) x2 = mu_min

      M = (x1 * x2 * (x1 + x2)) ** (K-1)
      COSPHI = COS(PHI)
      SIN1 = SQRT(1.0-x1**2)
      SIN2 = SQRT(1.0-x2**2)
      COSG = x1*x2 + SIN1*SIN2*COSPHI
      F = (1-THETA**2) / (1 + 2*THETA*COSG + THETA**2)**1.5
      TAN1 = SIN1/x1
      TAN2 = SIN2/x2
      CAPG = SQRT( ABS( TAN1**2 + TAN2**2 - 2*TAN1*TAN2*COSPHI ))
      H = 1 + (1-RHO0)/(1+CAPG)
      RPV_REFLECTION = RHO0 * M * F * H

      RETURN
      END   







      SUBROUTINE SUM_OUTPUT (IFLAG, WT, N, DATA, SUMDATA) 
C         Sums an output array over the k-distribution. 
C     If IFLAG=0 then the output array is set to the input times the weight.
      IMPLICIT NONE
      INTEGER IFLAG, N,  I
      REAL    WT, DATA(N), SUMDATA(N)
      IF (IFLAG .EQ. 0) THEN
        DO I = 1, N
          SUMDATA(I) = WT*DATA(I)
        ENDDO
      ELSE
        DO I = 1, N
          SUMDATA(I) = SUMDATA(I) + WT*DATA(I)
        ENDDO
      ENDIF
      RETURN
      END




      SUBROUTINE COMPUTE_NETFLUXDIV (NSTOKES, NPTS, RSHPTR, SRCTYPE, 
     .             SOLARMU, EXTINCT, ALBEDO, PLANCK, DIRFLUX,
     .             RADIANCE,  NETFLUXDIV)
C       Computes the net flux divergence at every grid point.
      IMPLICIT NONE
      INTEGER NSTOKES, NPTS, RSHPTR(NPTS+1)
      REAL    PLANCK(NPTS), DIRFLUX(NPTS), SOLARMU
      REAL    EXTINCT(NPTS), ALBEDO(NPTS)
      REAL    RADIANCE(NSTOKES,*), NETFLUXDIV(NPTS)
      CHARACTER SRCTYPE*1
      INTEGER I, IR
      REAL    PI, C, SECMU, HR

      PI = ACOS(-1.0)
      C = SQRT(4.0*PI)
      SECMU = 1.0/ABS(SOLARMU)
      DO I = 1, NPTS
        IR = RSHPTR(I)
        HR = C*EXTINCT(I)*(1-ALBEDO(I))*RADIANCE(1,IR+1)
        IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
          HR = HR + EXTINCT(I)*(1.0-ALBEDO(I))*DIRFLUX(I)*SECMU
        ENDIF
        IF (SRCTYPE .EQ. 'T' .OR. SRCTYPE .EQ. 'B') THEN
          HR = HR - 4*PI*EXTINCT(I)*PLANCK(I)
        ENDIF
        NETFLUXDIV(I) = -HR
      ENDDO
      RETURN
      END




      SUBROUTINE COMPUTE_SH (NSHOUT, NSTOKES, NPTS, SRCTYPE, 
     .                       SOLARMU, SOLARAZ, DIRFLUX,
     .                       RSHPTR, ML, MM, RADIANCE, SHTERMS)
C       Computes the quantities for the spherical harmonic output.
C     At each grid point the mean radiance and net flux (Fx, Fy, Fz)
C     are computed.  The output includes the diffuse and direct solar
C     components.  If NSHOUT=5 then the HIGHORDERRAD flag is set and the 
C     rms of the higher order terms in the radiance series is computed.
      IMPLICIT NONE
      INTEGER NSTOKES, NPTS, RSHPTR(NPTS+1), ML, MM, NSHOUT
      REAL    SOLARMU, SOLARAZ, DIRFLUX(NPTS)
      REAL    RADIANCE(NSTOKES,*), SHTERMS(NSHOUT,NPTS)
      CHARACTER SRCTYPE*1
      INTEGER I, IR, J, NLM
      REAL    PI, A, C, C0, IMEAN, FX, FY, FZ, HORMS
      REAL    SECSOL, SINSOL, SOLX, SOLY, SOLM, F0
      REAL, ALLOCATABLE :: YLMSUN(:)
 
      NLM = (2*MM+1)*(ML+1) - MM*(MM+1)
      ALLOCATE (YLMSUN(NLM))

C             Spherical Harmonic output: individual SH terms
      PI = ACOS(-1.0)
      A = SQRT(2.0*PI/3.0)
      C = SQRT(4.0*PI/3.0)
      C0 = 1.0/SQRT(4.0*PI)
      SECSOL = 1.0/ABS(SOLARMU)
      SINSOL = SQRT(1.0-SOLARMU**2)
      SOLX = SINSOL*COS(SOLARAZ)*SECSOL
      SOLY = SINSOL*SIN(SOLARAZ)*SECSOL
      SOLM = SECSOL/(4.0*PI)
      F0 = 0.0
      IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
        CALL YLMALL (.TRUE., SOLARMU, SOLARAZ, ML, MM, 1, YLMSUN)
      ENDIF
      DO I = 1, NPTS
        IR = RSHPTR(I)
        IMEAN = C0*RADIANCE(1,IR+1)
        FZ = C*RADIANCE(1,IR+3)
        IF (MM .GT. 0) THEN
c          FX = -C*RADIANCE(1,IR+4)   ! old spherical harmonics basis
          FX = -A*(RADIANCE(1,IR+4)+RADIANCE(1,IR+2))   ! real generalized spherical harmonics basis
        ELSE
          FX = 0.0
        ENDIF
        IF (MM .GT. 0) THEN
c          FY = -C*RADIANCE(1,IR+2)
          FY = -A*(RADIANCE(1,IR+2)-RADIANCE(1,IR+4))   ! real generalized spherical harmonics basis
        ELSE
          FY = 0.0
        ENDIF
        IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
          IMEAN = IMEAN + DIRFLUX(I)*SOLM
          FX = FX + DIRFLUX(I)*SOLX
          FY = FY + DIRFLUX(I)*SOLY
          FZ = FZ - DIRFLUX(I)
        ENDIF
        SHTERMS(1,I) = IMEAN
        SHTERMS(2,I) = FX
        SHTERMS(3,I) = FY
        SHTERMS(4,I) = FZ
        IF (NSHOUT .EQ. 5) THEN
          HORMS = 0.0
          IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
            F0 = SECSOL*DIRFLUX(I)
          ENDIF
          DO J = 5, RSHPTR(I+1)-IR
            HORMS = HORMS + (F0*YLMSUN(J)+RADIANCE(1,IR+J))**2
          ENDDO
          SHTERMS(5,I) = SQRT(HORMS)
        ENDIF
      ENDDO
      DEALLOCATE (YLMSUN)
      RETURN
      END





      SUBROUTINE INTERP_OUTPUT (OLDNPTS, NPTS, 
     .                  FLUX_OUT, FLUXDIV_OUT, SH_OUT,
     .                  FLUXES, DIRFLUX, FLUXDIV, SHTERMS,
     .                  NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
C       Interpolates the output data to the new grid points (those >OLDNPTS)
C     from the old grid points.  The output data is in FLUXES, DIRFLUX,
C     FLUXDIV, and SHTERMS.  OLDNPTS is set to zero to do no interpolation.
C     Each new point is located in a cell, and the two parent grid points
C     are found to interpolate the data from.
      IMPLICIT NONE
      INTEGER OLDNPTS, NPTS, NBCELLS, NCELLS
      INTEGER GRIDPTR(8,NCELLS), TREEPTR(2,NCELLS)
      LOGICAL FLUX_OUT, FLUXDIV_OUT, SH_OUT
      REAL    FLUXES(2,NPTS), DIRFLUX(NPTS), FLUXDIV(NPTS)
      REAL    SHTERMS(5,NPTS)
      REAL    GRIDPOS(3,*)
      INTEGER IP, IC, ICELL, IPARENT, IDIR, I, J, DIR1, DIR2, IP1, IP2
      INTEGER GRIDCORNER(2,4,3)
      DATA GRIDCORNER/1,2, 3,4, 5,6, 7,8,  1,3, 2,4, 5,7, 6,8,
     .                1,5, 2,6, 3,7, 4,8/

      IF (OLDNPTS .LE. 0)  RETURN

C         Go through all the non-base grid cells
      DO ICELL = NBCELLS+1, NCELLS
        DO IC = 1, 8
          IP = GRIDPTR(IC,ICELL)
          IF (IP .GT. OLDNPTS) THEN
            IPARENT = TREEPTR(1,ICELL)
C             Go through all 12 edges of parent cell to find which one has point
            DO IDIR = 1, 3
              DIR1 = MOD(IDIR-1+1,3)+1
              DIR2 = MOD(IDIR-1+2,3)+1
              DO I = 1, 4
                IP1 = GRIDPTR(GRIDCORNER(1,I,IDIR),IPARENT)
                IP2 = GRIDPTR(GRIDCORNER(2,I,IDIR),IPARENT)
                IF (GRIDPOS(DIR1,IP1) .EQ. GRIDPOS(DIR1,IP) .AND.
     .              GRIDPOS(DIR2,IP1) .EQ. GRIDPOS(DIR2,IP))  GOTO 120
              ENDDO            
            ENDDO
            STOP 'INTERP_OUTPUT: point not on an edge'
120         CONTINUE
C            Interpolate the data from the two parent points at the
C              ends of the edge; assume the cell was split in half.
            IF (FLUX_OUT) THEN
              FLUXES(1,IP) = 0.5*(FLUXES(1,IP1)+FLUXES(1,IP2))
              FLUXES(2,IP) = 0.5*(FLUXES(2,IP1)+FLUXES(2,IP2))
              DIRFLUX(IP) = 0.5*(DIRFLUX(IP1)+DIRFLUX(IP2))
            ENDIF
            IF (FLUXDIV_OUT) THEN
              FLUXDIV(IP) = 0.5*(FLUXDIV(IP1)+FLUXDIV(IP2))
            ENDIF
            IF (SH_OUT) THEN
              DO J = 1, 5
                SHTERMS(J,IP) = 0.5*(SHTERMS(J,IP1)+SHTERMS(J,IP2))
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END




      SUBROUTINE INTERP_RADIANCE (NSTOKES, OLDNPTS, NPTS, RSHPTR, 
     .          RADIANCE, NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
C       Interpolates the spherical harmonic radiance array to the new
C     grid points (>OLDNPTS) from the old grid points.
      IMPLICIT NONE
      INTEGER NSTOKES, OLDNPTS, NPTS, NBCELLS, NCELLS
      INTEGER RSHPTR(NPTS+1), GRIDPTR(8,NCELLS), TREEPTR(2,NCELLS)
      REAL    RADIANCE(NSTOKES,*), GRIDPOS(3,*)
      INTEGER IP, IC, ICELL, IPARENT, IDIR, I, J, K, DIR1,DIR2, IP1,IP2
      INTEGER IR, IR1, IR2, NR, NR1, NR2
      INTEGER GRIDCORNER(2,4,3)
      REAL    RAD1(NSTOKES), RAD2(NSTOKES)
      DATA GRIDCORNER/1,2, 3,4, 5,6, 7,8,  1,3, 2,4, 5,7, 6,8,
     .                1,5, 2,6, 3,7, 4,8/

      IF (OLDNPTS .LE. 0)  RETURN

C         Go through all the non-base grid cells
      DO ICELL = NBCELLS+1, NCELLS
        DO IC = 1, 8
          IP = GRIDPTR(IC,ICELL)
          IF (IP .GT. OLDNPTS) THEN
            IPARENT = TREEPTR(1,ICELL)
C             Go through all 12 edges of parent cell to find which one has point
            DO IDIR = 1, 3
              DIR1 = MOD(IDIR-1+1,3)+1
              DIR2 = MOD(IDIR-1+2,3)+1
              DO I = 1, 4
                IP1 = GRIDPTR(GRIDCORNER(1,I,IDIR),IPARENT)
                IP2 = GRIDPTR(GRIDCORNER(2,I,IDIR),IPARENT)
                IF (GRIDPOS(DIR1,IP1) .EQ. GRIDPOS(DIR1,IP) .AND.
     .              GRIDPOS(DIR2,IP1) .EQ. GRIDPOS(DIR2,IP))  GOTO 120
              ENDDO            
            ENDDO
            STOP 'INTERP_OUTPUT: point not on an edge'
120         CONTINUE
C            Interpolate the data from the two parent points at the
C              ends of the edge; assume the cell was split in half.
            IR1 = RSHPTR(IP1)
            IR2 = RSHPTR(IP2)
            NR1 = RSHPTR(IP1+1)-RSHPTR(IP1)
            NR2 = RSHPTR(IP2+1)-RSHPTR(IP2)
            NR = MAX(NR1,NR2)
            IR = RSHPTR(IP)
            RSHPTR(IP+1) = IR + NR
            DO J = 1, NR
              IF (J .LE. NR1) THEN
                RAD1(:) = RADIANCE(:,IR1+J)
              ELSE
                RAD1(:) = 0.0
              ENDIF
              IF (J .LE. NR2) THEN
                RAD2(:) = RADIANCE(:,IR2+J)
              ELSE
                RAD2(:) = 0.0
              ENDIF
              RADIANCE(:,IR+J) = 0.5*(RAD1(:)+RAD2(:))
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END







      SUBROUTINE VISUALIZE_RADIANCE (NSTOKES, NX, NY, NZ, NPTS, NCELLS,
     .             ML, MM, NLM, NSTLEG, NLEG, NUMPHASE, 
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, 
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,
     .             XGRID, YGRID, ZGRID, GRIDPOS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,
     .             SHPTR, SOURCE,  OUTPARMS,  IVIS, VISOUT)
C       Computes Stokes radiances (output in VISOUT) for the visualization 
C      modes: 1) camera mode, and 2) cross track scanning.
      IMPLICIT NONE
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
      INTEGER ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), IVIS
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SHPTR(NPTS+1), BCPTR(MAXNBC,2)
      INTEGER*2 CELLFLAGS(NCELLS), IPHASE(NPTS)
      LOGICAL DELTAM
      REAL    SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
      REAL    MU(NMU), PHI(NMU,NPHI0MAX), WTDO(NMU,NPHI0MAX)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    SFCGRIDPARMS(*), BCRAD(NSTOKES,*)
      REAL    EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG,0:NLEG,*)
      REAL    DIRFLUX(NPTS), FLUXES(2,NPTS), SOURCE(NSTOKES,*)
      REAL    OUTPARMS(*), VISOUT(NSTOKES,*)
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
 
      INTEGER NSCATANGLE, NSTPHASE, I, J, L, N, SIDE
      INTEGER NL, NS, LINE, SAMP
      LOGICAL CAMERA_MODE, VALIDRAD
      REAL    MURAY, PHIRAY
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT, VISRAD(NSTOKES)
      DOUBLE PRECISION X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2
      DOUBLE PRECISION THETA0, THETA1, PHIR, PHI0
      DOUBLE PRECISION COSTH, SINTH, SINTH0, SINTH1, MU2, PHI2
      DOUBLE PRECISION U, V, UP, VP, COSDELPHI, PI, ROTANG, DEGRAD, R
      DOUBLE PRECISION DIST, D, RX, RY, RZ, SCANANG
      INTEGER MAXSCATANG
      PARAMETER (MAXSCATANG=721)
      REAL, ALLOCATABLE :: YLMSUN(:,:), PHASETAB(:,:,:)
      REAL              :: MEAN, STD1, STD2
      REAL, ALLOCATABLE :: AOLP(:)


      ALLOCATE (YLMSUN(NSTLEG,NLM))

      IF (SRCTYPE .NE. 'T') THEN
        CALL YLMALL (.TRUE., SOLARMU, SOLARAZ, ML, MM, NSTLEG, YLMSUN)
        IF (DELTAM .AND. NUMPHASE .GT. 0) THEN
          NSCATANGLE = MAX(36,MIN(MAXSCATANG,2*NLEG))
          NSTPHASE = MIN(NSTLEG,2)
          ALLOCATE (PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE))
          CALL PRECOMPUTE_PHASE (NSCATANGLE, NUMPHASE, NSTPHASE, 
     .                    NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)
        ENDIF
      ENDIF


C         Make the isotropic radiances for the top boundary
      CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN, 
     .                            UNITS, NTOPPTS, NSTOKES, BCRAD(1,1))
C         Make the bottom boundary radiances for the Lambertian surfaces.  
C          Compute the upwelling bottom radiances using the downwelling fluxes.
      IF (SFCTYPE .EQ. 'FL') THEN
        CALL FIXED_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2),
     .             DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO, 
     .             WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1,1+NTOPPTS))
      ELSE IF (SFCTYPE .EQ. 'VL') THEN
        CALL VARIABLE_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2),
     .               DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,
     .               NSTOKES, BCRAD(1,1+NTOPPTS))
      ENDIF


      CAMERA_MODE = NINT(OUTPARMS(1)) .EQ. 1
      IF (CAMERA_MODE) THEN
        NL = NINT(OUTPARMS(10))
        NS = NINT(OUTPARMS(11))
      ELSE
        NL = 1 + INT( SQRT((OUTPARMS(4)-OUTPARMS(7))**2
     .               +(OUTPARMS(5)-OUTPARMS(8))**2
     .               +(OUTPARMS(6)-OUTPARMS(9))**2) /OUTPARMS(10) )
        NS = 1 + INT(ABS(OUTPARMS(12)-OUTPARMS(11))/OUTPARMS(13))
      ENDIF
      PI = ACOS(-1.0D0)
      DEGRAD = PI/180.


C         Loop over pixels in image
      DO LINE = 1, NL
      DO SAMP = 1, NS
        IVIS = IVIS + 1

        IF (CAMERA_MODE) THEN
C         Camera mode:
C          1, bytes, scale, X,Y,Z, theta, phi, rotang, NL, NS, delline, delsamp
C
C           Use spherical trig to find the pixel direction (MURAY,PHIRAY)
C             from the camera center (THETA0,PHI0) and the relative pixel
C             angles (U,V).
          UP = (SAMP-NS/2-1)*OUTPARMS(13)*DEGRAD
          VP = (LINE-NL/2-1)*OUTPARMS(12)*DEGRAD
          ROTANG = OUTPARMS(9)*DEGRAD
          IF (ROTANG .EQ. 0.0) THEN
            U = UP
            V = VP
          ELSE
            U = COS(ROTANG)*UP - SIN(ROTANG)*VP
            V = SIN(ROTANG)*UP + COS(ROTANG)*VP
          ENDIF
          THETA0 = DEGRAD*OUTPARMS(7)
          PHI0 = DEGRAD*OUTPARMS(8)
          THETA1 = THETA0 + V
          IF (V .EQ. 0.0) THEN
            COSTH = COS(U)*COS(THETA0)
          ELSE
            COSTH = COS(U)*(SIN(THETA1)*COS(V)-SIN(THETA0))/SIN(V)
            COSTH = MIN(+1.0D0,MAX(-1.0D0,COSTH))
          ENDIF
          SINTH = SQRT(1-COSTH**2)
          SINTH0 = SIN(THETA0)
          SINTH1 = SIN(THETA1)
          IF (ABS(SINTH) .LT. 1.0E-6) THEN
            PHIR = 0.0
          ELSE 
            IF (ABS(SINTH0).LT.1E-6 .AND. ABS(SINTH1).LE.1E-6) THEN
              COSDELPHI = 0.0D0
            ELSE IF (ABS(SINTH1) .GT. 1.0E-6) THEN
              COSDELPHI = (COS(U)-COSTH*COS(THETA1))/(SINTH1*SINTH)
            ELSE IF (ABS(SINTH0) .GT. 1.0E-6) THEN
              COSDELPHI = (COS(U)*COS(V)-COSTH*COS(THETA0))
     .                     /(SINTH0*SINTH)
            ENDIF
            COSDELPHI = MIN(+1.0D0,MAX(-1.0D0,COSDELPHI))
            IF (U .GE. 0.0) THEN
              PHIR = PHI0 - ACOS(COSDELPHI)
            ELSE
              PHIR = PHI0 + ACOS(COSDELPHI)
            ENDIF
          ENDIF

          X0 = OUTPARMS(4)
          Y0 = OUTPARMS(5)
          Z0 = OUTPARMS(6)

        ELSE

C         Cross track scanning in the vertical plane:
C          2, bytes, scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1, scan2, delscan
C            Start and end of scan angle (scan1, scan2) are +/- relative
C            to nadir, and positive is on right side.
          X1 = OUTPARMS(4)
          Y1 = OUTPARMS(5)
          Z1 = OUTPARMS(6)
          X2 = OUTPARMS(7)
          Y2 = OUTPARMS(8)
          Z2 = OUTPARMS(9)
          DIST = SQRT( (X2-X1)**2 + (Y2-Y1)**2 + (Z2-Z1)**2)
          RX = (X2-X1)/DIST
          RY = (Y2-Y1)/DIST
          RZ = (Z2-Z1)/DIST
c         D = (LINE-1)*OUTPARMS(10)
          D = (NL-LINE)*OUTPARMS(10)
          X0 = X1 + D*RX
          Y0 = Y1 + D*RY
          Z0 = Z1 + D*RZ
          SCANANG = DEGRAD*(OUTPARMS(11) + (SAMP-1)*OUTPARMS(13))
          COSTH = COS(PI-ABS(SCANANG))
          IF (SCANANG .GT. 0.0) THEN
            PHIR = ATAN2(-RX,RY)
          ELSE
            PHIR = ATAN2(RX,-RY)
          ENDIF
        ENDIF


C             Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (COSTH .GE. 0.0) THEN
            VISRAD(:) = 0.0
            GOTO 900
          ENDIF
          R = (ZGRID(NZ) - Z0)/COSTH
          X0 = X0 + R*SQRT(1-COSTH**2)*COS(PHIR)
          Y0 = Y0 + R*SQRT(1-COSTH**2)*SIN(PHIR)
          Z0 = ZGRID(NZ)
        ELSE IF (Z0 .LT. ZGRID(1)) THEN
          WRITE (6,*) 'VISUALIZE_RADIANCE: Level below domain'
          STOP
        ENDIF

C         MURAY,PHRAY is camera pixel viewing direction; 
C         MU2,PHI2 is radiance travel direction (the opposite)
        MURAY = SNGL(COSTH)
        PHIRAY = SNGL(PHIR)
        MU2 = -COSTH
        PHI2 = PHIR + PI

C         Integrate the extinction and source function along this ray
C         to calculate the Stokes radiance vector for this pixel
        TRANSMIT = 1.0D0 ; VISRAD(:) = 0.0D0
        CALL INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES, NSTLEG, 
     .                       NSTPHASE, NSCATANGLE, PHASETAB,
     .                       NX, NY, NZ, NPTS, NCELLS, 
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       XGRID, YGRID, ZGRID, GRIDPOS,
     .                       ML, MM, NLM, NLEG, NUMPHASE,
     .                       NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                       DELTAM, SRCTYPE, WAVELEN, SOLARMU,SOLARAZ,
     .                       EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .                       DIRFLUX, SHPTR, SOURCE, YLMSUN, 
     .                       MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                       SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                       MU2, PHI2, X0,Y0,Z0, 
     .                     XE,YE,ZE, SIDE, TRANSMIT, VISRAD, VALIDRAD)
900     CONTINUE
c        WRITE (6,'(1X,2F8.4,1X,2F11.7,4(1X,F11.6))') 
c     .         X0,Y0,MURAY,PHIRAY,VISRAD(:)
        VISOUT(1,IVIS) = VISRAD(1)
        IF (VISRAD(1) .GT. 0.0) THEN
         IF (NSTOKES .GT. 1) THEN
C           Output degree (0 to 1) and angle (-180 to 180) of linear polarization
           VISOUT(2,IVIS) = SQRT(VISRAD(2)**2+VISRAD(3)**2)/VISRAD(1)
           VISOUT(3,IVIS) = (180/PI)*0.5*ATAN2(VISRAD(3),VISRAD(2))
         ENDIF
         IF (NSTOKES .EQ. 4) THEN
C           Output degree of circular polarization (-1 to 1)
           VISOUT(4,IVIS) = VISRAD(4)/VISRAD(1)
         ENDIF    
        ELSE
          VISOUT(2:,IVIS) = 0.0
        ENDIF
      ENDDO
      ENDDO

      DEALLOCATE (YLMSUN)
      IF (ALLOCATED(PHASETAB))  DEALLOCATE (PHASETAB)
      IF (NSTOKES .GT. 1) THEN
C        Choose the best range for the angle of linear polarization (-90 to 90 or 0 to 180)
        N = NL*NS
        ALLOCATE (AOLP(N))
        AOLP(:) = VISOUT(3,IVIS-N+1:IVIS)
        MEAN = SUM(AOLP(:))/N
        STD1 = SQRT(SUM((AOLP(:)-MEAN)**2)/N)
        WHERE (AOLP(:) < 0.0)
          AOLP(:) = AOLP(:)+180.0
        END WHERE
        MEAN = SUM(AOLP(:))/N
        STD2 = SQRT(SUM((AOLP(:)-MEAN)**2)/N)
        IF (STD2 < STD1) THEN
          VISOUT(3,IVIS-N+1:IVIS) = AOLP(:)
        ENDIF
        DEALLOCATE (AOLP)
      ENDIF
      RETURN
      END




      SUBROUTINE INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES, NSTLEG, 
     .                        NSTPHASE, NSCATANGLE, PHASETAB,
     .                        NX, NY, NZ, NPTS, NCELLS, 
     .                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                        XGRID, YGRID, ZGRID, GRIDPOS,
     .                        ML, MM, NLM, NLEG, NUMPHASE,
     .                        NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                        DELTAM, SRCTYPE, WAVELEN,SOLARMU,SOLARAZ,
     .                        EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .                        DIRFLUX, SHPTR, SOURCE, YLMSUN, 
     .                        MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                        SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                        MU2, PHI2, X0,Y0,Z0, 
     .                   XE,YE,ZE, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
C       Integrates the source function through the extinction field 
C     (EXTINCT) backward from the outgoing direction (MU2,PHI2) to find the 
C     radiance (RADIANCE) at the point X0,Y0,Z0.
C     The transmission and radiance of the ray so far (TRANSMIT, RADIANCE)
C     are input and returned after the integration along with the exitting
C     ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
C     5=-Z,6=+Z).
      IMPLICIT NONE
      INTEGER BCFLAG, IPFLAG, NSTOKES, NSTLEG, NSTPHASE, NSCATANGLE
      INTEGER NX, NY, NZ, NPTS, NCELLS, SIDE
      INTEGER ML, MM, NLM, NLEG, NUMPHASE
      INTEGER MAXNBC, NTOPPTS, NBOTPTS
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SHPTR(NPTS+1)
      INTEGER*2 CELLFLAGS(NCELLS)
      INTEGER IPHASE(NPTS)
      INTEGER BCPTR(MAXNBC,2)
      LOGICAL DELTAM, VALIDRAD
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    WAVELEN, SOLARMU, SOLARAZ
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG,0:NLEG,NPTS)
      REAL    DIRFLUX(NPTS), SOURCE(NSTOKES,*)
      REAL    YLMSUN(NSTLEG,NLM)
      REAL    PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
      DOUBLE PRECISION MU2, PHI2, X0,Y0,Z0, XE,YE,ZE
      DOUBLE PRECISION TRANSMIT, RADIANCE(NSTOKES)
      CHARACTER SRCTYPE*1, SFCTYPE*2

      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, J, L
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), DONETHIS(8)
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    OEXTINCT8(8), OSRCEXT8(NSTOKES,8)
      REAL    EXTINCT8(8), SRCEXT8(NSTOKES,8)
      REAL    EXT0, EXT1, EXTN
      REAL    SRCEXT0(NSTOKES), SRCEXT1(NSTOKES), RADBND(NSTOKES)
      REAL    XM,YM, F
      DOUBLE PRECISION PI, CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
      DOUBLE PRECISION TAUTOL, TAUGRID, S, DELS, TRANSCUT
      DOUBLE PRECISION EXT, TAU, TRANSCELL, ABSCELL
      DOUBLE PRECISION SRC(NSTOKES)
      DOUBLE PRECISION U,V,W, DELX,DELY,DELZ, INVDELX,INVDELY,INVDELZ
      DOUBLE PRECISION COSSCAT
      REAL, ALLOCATABLE :: YLMDIR(:,:), SINGSCAT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: SUNDIRLEG(:)

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
      PI = ACOS(-1.0D0)

C       Calculate the generalized spherical harmonics for this direction
      ALLOCATE (YLMDIR(NSTLEG,NLM))
      CALL YLMALL (.FALSE.,SNGL(MU2),SNGL(PHI2),ML,MM,NSTLEG, YLMDIR)

      IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
C          Get the solar single scattering Stokes vector for the outgoing
C          direction my interpolating in scattering angle in the PHASETAB 
C          table and then rotating the Q/U polarization to the outgoing plane.
        COSSCAT = SOLARMU*MU2 + SQRT((1.0-SOLARMU**2)*(1.0-MU2**2))
     .                  *COS(SOLARAZ-PHI2)
        IF (NUMPHASE .GT. 0) THEN
          ALLOCATE (SINGSCAT(NSTOKES,NUMPHASE))
          F = (NSCATANGLE-1)*(ACOS(COSSCAT)/PI) + 1
          J = MIN(NSCATANGLE-1,INT(F))
          F = F - J
          DO I = 1, NUMPHASE
            SINGSCAT(1:NSTPHASE,I) 
     .               = (1-F)*PHASETAB(:,I,J) + F*PHASETAB(:,I,J+1)
            IF (NSTOKES .GT. 1) THEN
              CALL ROTATE_POL_PLANE (NSTOKES, COSSCAT, SOLARMU,
     .                   SNGL(MU2), SOLARAZ-SNGL(PHI2), SINGSCAT(:,I))
            ENDIF
          ENDDO
        ELSE
          ALLOCATE (SUNDIRLEG(0:NLEG))
          CALL LEGENDRE_ALL (COSSCAT, NLEG, SUNDIRLEG)
          DO L = 0, NLEG
            SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L+1)/(4*PI)
          ENDDO
        ENDIF
      ENDIF

C         Make the ray direction (opposite to the outgoing direction)
      CX = SQRT(1.0D0-MU2**2)*COS(PHI2-PI)
      CY = SQRT(1.0D0-MU2**2)*SIN(PHI2-PI)
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

C           Decide which of the eight grid points we need the source function
        DO I = 1, 8
          DONETHIS(I) = DONEFACE(I,IFACE+1)
          IF (NX .EQ. 1 .AND. ONEX(I) .LT. 0) DONETHIS(I) = ONEX(I)
          IF (NY .EQ. 1 .AND. ONEY(I) .LT. 0) DONETHIS(I) = ONEY(I)
          OEXTINCT8(I) = EXTINCT8(I)
          OSRCEXT8(:,I) = SRCEXT8(:,I)
        ENDDO
C         Compute the source function times extinction in direction (MU2,PHI2)
        IF (NSTOKES .EQ. 1) THEN
          CALL COMPUTE_SOURCE_1CELL_UNPOL (ICELL, GRIDPTR, 
     .             ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU, 
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, 
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8, 
     .             EXTINCT8, SRCEXT8)
        ELSE
          CALL COMPUTE_SOURCE_1CELL (ICELL, GRIDPTR, 
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU, 
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, 
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8, 
     .             EXTINCT8, SRCEXT8)
        ENDIF

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
        SRCEXT1(:)=(1-W)*((1-V)*((1-U)*SRCEXT8(:,1) + U*SRCEXT8(:,2))
     .                     + V*((1-U)*SRCEXT8(:,3) + U*SRCEXT8(:,4))) 
     .              + W*((1-V)*((1-U)*SRCEXT8(:,5) + U*SRCEXT8(:,6))
     .                     + V*((1-U)*SRCEXT8(:,7) + U*SRCEXT8(:,8)))
        SRCEXT1(1) = MAX(0.0,SRCEXT1(1))
        EXT1 = (1-W)*((1-V)*((1-U)*EXTINCT8(1) + U*EXTINCT8(2))
     .                  + V*((1-U)*EXTINCT8(3) + U*EXTINCT8(4))) 
     .           + W*((1-V)*((1-U)*EXTINCT8(5) + U*EXTINCT8(6))
     .                  + V*((1-U)*EXTINCT8(7) + U*EXTINCT8(8)))

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
        U = (XN-GRIDPOS(1,IPT1))*INVDELX
        V = (YN-GRIDPOS(2,IPT1))*INVDELY
        W = (ZN-GRIDPOS(3,IPT1))*INVDELZ
        EXTN = (1-W)*((1-V)*((1-U)*EXTINCT8(1) + U*EXTINCT8(2))
     .                  + V*((1-U)*EXTINCT8(3) + U*EXTINCT8(4))) 
     .           + W*((1-V)*((1-U)*EXTINCT8(5) + U*EXTINCT8(6))
     .                  + V*((1-U)*EXTINCT8(7) + U*EXTINCT8(8)))
        TAUGRID = SO*0.5*(EXT1+EXTN)
        NTAU = MAX(1,1+INT(TAUGRID/TAUTOL))
        DELS = SO/NTAU 

C           Loop over the subgrid cells
        DO IT = 1, NTAU
          S = IT*DELS
          XI = XE + S*CX
          YI = YE + S*CY
          ZI = ZE + S*CZ
C            Interpolate extinction and source function along path
          U = (XI-GRIDPOS(1,IPT1))*INVDELX
          V = (YI-GRIDPOS(2,IPT1))*INVDELY
          W = (ZI-GRIDPOS(3,IPT1))*INVDELZ
          EXT0 = (1-W)*((1-V)*((1-U)*EXTINCT8(1) + U*EXTINCT8(2))
     .                    + V*((1-U)*EXTINCT8(3) + U*EXTINCT8(4))) 
     .             + W*((1-V)*((1-U)*EXTINCT8(5) + U*EXTINCT8(6))
     .                    + V*((1-U)*EXTINCT8(7) + U*EXTINCT8(8)))
          SRCEXT0(:)=(1-W)*((1-V)*((1-U)*SRCEXT8(:,1) + U*SRCEXT8(:,2))
     .                    + V*((1-U)*SRCEXT8(:,3) + U*SRCEXT8(:,4))) 
     .             + W*((1-V)*((1-U)*SRCEXT8(:,5) + U*SRCEXT8(:,6))
     .                    + V*((1-U)*SRCEXT8(:,7) + U*SRCEXT8(:,8)))
          SRCEXT0(1) = MAX(0.0,SRCEXT0(1))

C            Compute the subgrid radiance: integration of the source function
          EXT = 0.5*(EXT0+EXT1)
          IF (EXT .NE. 0.0) THEN
            TAU=EXT*DELS
            ABSCELL = TAU*(1.0-0.5*TAU*(1.0-0.33333333333*TAU))
            TRANSCELL = 1.0 - ABSCELL
C                 Linear extinction, linear source*extinction, to second order
            SRC(:) = ( 0.5*(SRCEXT0(:)+SRCEXT1(:)) 
     .           + 0.08333333333*(EXT0*SRCEXT1(:)-EXT1*SRCEXT0(:))*DELS
     .                *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
            SRC(:) = 0.0
          ENDIF

          RADIANCE(:) = RADIANCE(:) + TRANSMIT*SRC(:)*ABSCELL
          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          SRCEXT1(:) = SRCEXT0(:)
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
          VALIDRAD = .TRUE.
        ELSE IF (INEXTCELL .EQ. 0 .AND. IFACE .GE. 5) THEN
          VALIDRAD = .TRUE.
          CALL FIND_BOUNDARY_RADIANCE (NSTOKES, XN, YN, 
     .                      SNGL(MU2), SNGL(PHI2), 
     .                      IC, KFACE, GRIDPTR, GRIDPOS,
     .                      MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                      NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                      SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, DIRFLUX, 
     .                      SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                      RADBND)
          RADIANCE(:) = RADIANCE(:) + TRANSMIT*RADBND(:)
        ELSE
          ICELL = INEXTCELL
        ENDIF
        XE = XN
        YE = YN
        ZE = ZN

      ENDDO

      SIDE = IFACE

      IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
        IF (NUMPHASE .GT. 0) THEN
          DEALLOCATE (SINGSCAT)
        ELSE
          DEALLOCATE (SUNDIRLEG)
        ENDIF
      ENDIF
      DEALLOCATE (YLMDIR)
      RETURN
      END




      SUBROUTINE FIND_BOUNDARY_RADIANCE (NSTOKES, XB, YB, MU2, PHI2, 
     .                      ICELL, KFACE, GRIDPTR, GRIDPOS,
     .                      MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                      NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                      SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, DIRFLUX, 
     .                      SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                      RADBND)
C       Returns the interpolated Stokes radiance at the boundary (RADBND).
C     Inputs are the boundary location (XB,YB), ray direction away from
C     boundary (MU2,PHI2), cell number (ICELL) and face (KFACE) at
C     the boundary point.
      IMPLICIT NONE
      INTEGER NSTOKES, ICELL, KFACE, MAXNBC, NTOPPTS, NBOTPTS
      INTEGER GRIDPTR(8,*), BCPTR(MAXNBC,2)
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
      REAL    MU2, PHI2, RADBND(NSTOKES)
      DOUBLE PRECISION XB, YB
      REAL    GRIDPOS(3,*)
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(NSTOKES,*)
      CHARACTER SRCTYPE*1, SFCTYPE*2

      INTEGER IL, IM, IU, IP, IBC, J
      LOGICAL LAMBERTIAN
      REAL    X(4), Y(4), RAD(NSTOKES,4), U, V
      INTEGER GRIDFACE(4,6)
      DATA    GRIDFACE/1,3,5,7, 2,4,6,8,  1,2,5,6, 3,4,7,8,
     .                 1,2,3,4, 5,6,7,8/

      LAMBERTIAN = SFCTYPE(2:2) .EQ. 'L'

C       Loop over the four gridpoints of the cell face on the boundary
      DO J = 1, 4
        IP = GRIDPTR(GRIDFACE(J,KFACE),ICELL)
        X(J) = GRIDPOS(1,IP)
        Y(J) = GRIDPOS(2,IP)

        IF (MU2 .LT. 0.0) THEN
C           Do a binary search to locate the top boundary point
          IL = 1
          IU = NTOPPTS
          DO WHILE (IU-IL .GT. 1)
            IM = (IU+IL)/2
            IF (IP .GE. BCPTR(IM,1)) THEN
              IL = IM
            ELSE
              IU = IM
            ENDIF
          ENDDO
          IBC = IL
          IF (BCPTR(IBC,1) .NE. IP)  IBC=IU
          IF (BCPTR(IBC,1) .NE. IP)  
     .      STOP 'FIND_BOUNDARY_RADIANCE: Not at boundary'
          RAD(:,J) = BCRAD(:,IBC)

        ELSE

C           Do a binary search to locate the bottom boundary point
          IL = 1
          IU = NBOTPTS
          DO WHILE (IU-IL .GT. 1)
            IM = (IU+IL)/2
            IF (IP .GE. BCPTR(IM,2)) THEN
              IL = IM
            ELSE
              IU = IM
            ENDIF
          ENDDO
          IBC = IL
          IF (BCPTR(IBC,2) .NE. IP)  IBC=IU
          IF (BCPTR(IBC,2) .NE. IP)
     .      STOP 'FIND_BOUNDARY_RADIANCE: Not at boundary'

          IF (.NOT. LAMBERTIAN) THEN
            CALL VARIABLE_BRDF_SURFACE (NBOTPTS,IBC,IBC, BCPTR(1,2), 
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,
     .             SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX, 
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES,
     .             BCRAD(:,1+NTOPPTS))
          ENDIF
          RAD(:,J) = BCRAD(:,NTOPPTS+IBC)
        ENDIF
      ENDDO
      IF (X(2)-X(1) .GT. 0.0) THEN
        U = (XB-X(1))/(X(2)-X(1))
      ELSE
        U = 0.0
      ENDIF
      IF (Y(3)-Y(1) .GT. 0.0) THEN
        V = (YB-Y(1))/(Y(3)-Y(1))
      ELSE
        V = 0.0
      ENDIF
      RADBND(:) = (1-U)*(1-V)*RAD(:,1) + U*(1-V)*RAD(:,2)
     .              + (1-U)*V*RAD(:,3) +     U*V*RAD(:,4)
      RETURN
      END



 
      SUBROUTINE COMPUTE_SOURCE_1CELL (ICELL, GRIDPTR, 
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU, 
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, 
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8, 
     .             EXTINCT8, SRCEXT8)
C       Computes the source function times extinction for gridpoints 
C     belonging to cell ICELL in the direction (MU,PHI).  The results 
C     are returned in SRCEXT8 and EXTINCT8.
C     The spherical harmonic source function series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
      IMPLICIT NONE
      INTEGER ICELL, NSTOKES, NSTLEG, NPTS, ML, MM, NLM, NLEG, NUMPHASE
      INTEGER GRIDPTR(8,*), SHPTR(NPTS+1)
      INTEGER DONETHIS(8), OLDIPTS(8)
      INTEGER*2 IPHASE(NPTS)
      LOGICAL DELTAM
      REAL    SOLARMU
      REAL    EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG,0:NLEG,NPTS)
      REAL    DIRFLUX(NPTS), SOURCE(NSTOKES,*)
      REAL    YLMDIR(NSTLEG,NLM), YLMSUN(NSTLEG,NLM)
      REAL    SINGSCAT(NSTOKES,NUMPHASE)
      REAL    OEXTINCT8(8), OSRCEXT8(NSTOKES,8)
      REAL    EXTINCT8(8), SRCEXT8(NSTOKES,8)
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
      CHARACTER SRCTYPE*1

      INTEGER IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I
      REAL    SECMU0, F, DA, A, A1, B1
 
      SECMU0 = 1.0D0/ABS(SOLARMU)

C         Loop over the grid points, computing the source function 
C           at the viewing angle from the spherical harmonics source function.
      DO N = 1, 8
        IP = GRIDPTR(N,ICELL)
        I = DONETHIS(N)
        IF (I .GT. 0 .AND. IP .EQ. OLDIPTS(N)) THEN 
          EXTINCT8(N) = OEXTINCT8(I)
          SRCEXT8(:,N) = OSRCEXT8(:,I)
        ELSE IF (I .LT. 0) THEN
          EXTINCT8(N) = EXTINCT8(ABS(I))
          SRCEXT8(:,N) = SRCEXT8(:,ABS(I))
        ELSE

          OLDIPTS(N) = IP
          IS = SHPTR(IP)
          NS = SHPTR(IP+1)-IS
C             Sum over the real generalized spherical harmonic series 
C             of the source function
          SRCEXT8(:,N) = 0.0
          DO J = 1, NS
            SRCEXT8(1,N) = SRCEXT8(1,N) + SOURCE(1,IS+J)*YLMDIR(1,J)
          ENDDO
          IF (NSTOKES .GT. 1) THEN
            DO J = 1, NS
              SRCEXT8(2,N) = SRCEXT8(2,N) + SOURCE(2,IS+J)*YLMDIR(2,J)
     .                                    + SOURCE(3,IS+J)*YLMDIR(5,J)
              SRCEXT8(3,N) = SRCEXT8(3,N) + SOURCE(2,IS+J)*YLMDIR(6,J)
     .                                    + SOURCE(3,IS+J)*YLMDIR(3,J)
            ENDDO
          ENDIF
          IF (NSTOKES .EQ. 4) THEN
            DO J = 1, NS
              SRCEXT8(4,N) = SRCEXT8(4,N) + SOURCE(4,IS+J)*YLMDIR(4,J)
            ENDDO
          ENDIF

C             Special case for solar source and Delta-M
          IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
            IF (NUMPHASE .GT. 0) THEN
              IPH = IPHASE(IP)
            ELSE
              IPH = IP
            ENDIF
C               First subtract off the truncated single scattering 
            DA = ALBEDO(IP)*DIRFLUX(IP)*SECMU0
            J = 1
            DO L = 0, ML
              ME = MIN(L,MM)
              MS = -ME
              A1 = DA*LEGEN(1,L,IPH)
              B1 = DA*LEGEN(5,L,IPH)
              IF (J .LE. NS) THEN
                JT = J
                DO M = MS, ME
                  SRCEXT8(1,N) =SRCEXT8(1,N)-A1*YLMDIR(1,J)*YLMSUN(1,J)
                  J = J + 1
                ENDDO
                IF (NSTOKES .GT. 1) THEN
                  J = JT
                  DO M = MS, ME
                    SRCEXT8(2,N)=SRCEXT8(2,N)-B1*YLMDIR(2,J)*YLMSUN(1,J)
                    SRCEXT8(3,N)=SRCEXT8(3,N)-B1*YLMDIR(6,J)*YLMSUN(1,J)
                    J = J + 1
                  ENDDO
                ENDIF  
              ENDIF    
            ENDDO
C               Then add in the single scattering contribution for the
C               original unscaled phase function.  
            IF (NUMPHASE .GT. 0) THEN
              SRCEXT8(:,N) = SRCEXT8(:,N) + DA*SINGSCAT(:,IPH)
            ELSE IF (NSTOKES .EQ. 1) THEN
              F = LEGEN(1,ML+1,IPH)
              DO L = 0, NLEG
                IF (L .LE. ML) THEN
                  A = DA*(LEGEN(1,L,IPH) + F/(1-F))
                ELSE
                  A = DA*LEGEN(1,L,IPH)/(1-F)
                ENDIF
                SRCEXT8(1,N) = SRCEXT8(1,N) + A*SUNDIRLEG(L)
              ENDDO
            ENDIF
          ENDIF

          SRCEXT8(:,N) = SRCEXT8(:,N)*EXTINCT(IP)
          EXTINCT8(N) = EXTINCT(IP)
        ENDIF
      ENDDO
 
      RETURN
      END
 


      SUBROUTINE COMPUTE_SOURCE_1CELL_UNPOL (ICELL, GRIDPTR, 
     .             ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU, 
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, 
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8, 
     .             EXTINCT8, SRCEXT8)
C       Computes the source function times extinction for gridpoints 
C     belonging to cell ICELL in the direction (MU,PHI).  The results 
C     are returned in SRCEXT8 and EXTINCT8.
C     The spherical harmonic source function series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
      IMPLICIT NONE
      INTEGER ICELL, NPTS, ML, MM, NLM, NLEG, NUMPHASE
      INTEGER GRIDPTR(8,*), SHPTR(NPTS+1)
      INTEGER DONETHIS(8), OLDIPTS(8)
      INTEGER*2 IPHASE(NPTS)
      LOGICAL DELTAM
      REAL    SOLARMU
      REAL    EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(0:NLEG,NPTS)
      REAL    DIRFLUX(NPTS), SOURCE(*)
      REAL    YLMDIR(NLM), YLMSUN(NLM)
      REAL    SINGSCAT(NUMPHASE)
      REAL    OEXTINCT8(8), OSRCEXT8(8)
      REAL    EXTINCT8(8), SRCEXT8(8)
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
      CHARACTER SRCTYPE*1

      INTEGER IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I
      REAL    SECMU0, F, DA, A, A1, B1
 
      SECMU0 = 1.0D0/ABS(SOLARMU)

C         Loop over the grid points, computing the source function 
C           at the viewing angle from the spherical harmonics source function.
      DO N = 1, 8
        IP = GRIDPTR(N,ICELL)
        I = DONETHIS(N)
        IF (I .GT. 0 .AND. IP .EQ. OLDIPTS(N)) THEN 
          EXTINCT8(N) = OEXTINCT8(I)
          SRCEXT8(N) = OSRCEXT8(I)
        ELSE IF (I .LT. 0) THEN
          EXTINCT8(N) = EXTINCT8(ABS(I))
          SRCEXT8(N) = SRCEXT8(ABS(I))
        ELSE

          OLDIPTS(N) = IP
          IS = SHPTR(IP)
          NS = SHPTR(IP+1)-IS
C             Sum over the real generalized spherical harmonic series 
C             of the source function
          SRCEXT8(N) = 0.0
          DO J = 1, NS
            SRCEXT8(N) = SRCEXT8(N) + SOURCE(IS+J)*YLMDIR(J)
          ENDDO

C             Special case for solar source and Delta-M
          IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
            IF (NUMPHASE .GT. 0) THEN
              IPH = IPHASE(IP)
            ELSE
              IPH = IP
            ENDIF
C               First subtract off the truncated single scattering 
            DA = ALBEDO(IP)*DIRFLUX(IP)*SECMU0
            J = 1
            DO L = 0, ML
              ME = MIN(L,MM)
              MS = -ME
              A1 = DA*LEGEN(L,IPH)
              IF (J .LE. NS) THEN
                JT = J
                DO M = MS, ME
                  SRCEXT8(N) =SRCEXT8(N)-A1*YLMDIR(J)*YLMSUN(J)
                  J = J + 1
                ENDDO
              ENDIF    
            ENDDO
C               Then add in the single scattering contribution for the
C               original unscaled phase function.  
            IF (NUMPHASE .GT. 0) THEN
              SRCEXT8(N) = SRCEXT8(N) + DA*SINGSCAT(IPH)
            ELSE
              F = LEGEN(ML+1,IPH)
              DO L = 0, NLEG
                IF (L .LE. ML) THEN
                  A = DA*(LEGEN(L,IPH) + F/(1-F))
                ELSE
                  A = DA*LEGEN(L,IPH)/(1-F)
                ENDIF
                SRCEXT8(N) = SRCEXT8(N) + A*SUNDIRLEG(L)
              ENDDO
            ENDIF
          ENDIF

          SRCEXT8(N) = SRCEXT8(N)*EXTINCT(IP)
          EXTINCT8(N) = EXTINCT(IP)
        ENDIF
      ENDDO
 
      RETURN
      END
 


      SUBROUTINE PRECOMPUTE_PHASE (NSCATANGLE, NUMPHASE, NSTPHASE, 
     .                      NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)
C       Precomputes the phase matrix elements I-I and I-Q as a function
C     of scattering angle for solar direct scattering for all the 
C     tabulated phase functions. Output is in PHASETAB.
C     NSTPHASE=1 for NSTOKES=1 and NSTPHASE=2 for NSTOKES>1.  
      IMPLICIT NONE
      INTEGER NSCATANGLE, NUMPHASE, NSTPHASE
      INTEGER NSTOKES, NSTLEG, ML, NLEG
      REAL    LEGEN(NSTLEG,0:NLEG,NUMPHASE)
      REAL    PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
      INTEGER IPH, I, J, L
      DOUBLE PRECISION  PI, OFOURPI, COSSCAT, FCT, F, X, A1, B1
      DOUBLE PRECISION, ALLOCATABLE :: UNSCLEGEN(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: DMM1(:), DMM2(:)

      ALLOCATE (UNSCLEGEN(2,0:NLEG), DMM1(0:NLEG), DMM2(0:NLEG)) 
      PI = ACOS(-1.0D0)
      OFOURPI = 1.0/(4.0*PI)

      DO J = 1, NSCATANGLE
        COSSCAT = COS(PI*DFLOAT(J-1)/(NSCATANGLE-1))
        X = DBLE(COSSCAT)
        CALL WIGNERFCT (X, NLEG, 0, 0, DMM1)
        IF (NSTOKES .GT. 1) THEN
          CALL WIGNERFCT (X, NLEG, 2, 0, DMM2)
        ENDIF

C          Unscaled the needed Legendre coefficients, but divide by 1-f
C          because extinction is still scaled (for TMS method)
        DO IPH = 1, NUMPHASE
          F = LEGEN(1,ML+1,IPH)
          DO L = 0, NLEG
            IF (L .LE. ML) THEN
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH) + F/(1-F)
            ELSE
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)/(1-F)
            ENDIF
            IF (NSTLEG .GT. 1) THEN
              IF (L .LE. ML) THEN
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)
              ELSE
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)/(1-F)
              ENDIF
            ENDIF
          ENDDO

C           Sum the first Wigner function series, which is actually a Legendre series
          A1 = 0.0D0
          DO L = 0, NLEG
            FCT = 2.0D0*L + 1.0D0
            A1  = A1 + FCT*DBLE(UNSCLEGEN(1,L))*DMM1(L)
          ENDDO
          IF (A1 .LE. 0.0) THEN
            WRITE (6,*) 'PRECOMPUTE_PHASE: negative phase function',
     .          ' for tabulated phase function: ',IPH, J, A1
            STOP
          ENDIF
          PHASETAB(1,IPH,J) = SNGL(A1*OFOURPI)

          IF (NSTOKES .GT. 1) THEN
C           If doing polarization, sum the second Wigner function series
            B1 = 0.0D0
            DO L = 0, NLEG
              FCT = 2.0D0*L + 1.0D0
              B1  = B1 - FCT*DBLE(UNSCLEGEN(2,L))*DMM2(L)
            ENDDO
            PHASETAB(2,IPH,J) = SNGL(B1*OFOURPI)
          ENDIF
        ENDDO
      ENDDO

      DEALLOCATE (UNSCLEGEN, DMM1, DMM2) 
      RETURN
      END

 

      SUBROUTINE ROTATE_POL_PLANE (NSTOKES, COSSCAT,
     .                          SOLARMU, MU, DELPHI, SCATVECT)
C      Rotates the plane of polarization from the scattering plane
C      to the outgoing meridional plane for the solar single
C      scattered Stokes vector in SINGSCAT.
      IMPLICIT NONE
      INTEGER NSTOKES
      DOUBLE PRECISION COSSCAT
      REAL    SOLARMU, MU, DELPHI
      REAL    SCATVECT(NSTOKES)
      REAL    B1
      DOUBLE PRECISION SIN_SCAT, SIN_THETA1, SIN_THETA2
      DOUBLE PRECISION SINPHI, COSPHI, SIN2, COS2, SIN22, COS22

      IF (NSTOKES .GT. 1) THEN
        B1 = SCATVECT(2)
        SIN_SCAT = DSQRT(MAX(0.0D0,1.D0-COSSCAT**2))
        SIN_THETA1 = DSQRT(1.D0-SOLARMU**2)
        SIN_THETA2 = DSQRT(1.D0-MU**2)
        SINPHI = DSIN(DBLE(DELPHI))   
        COSPHI = DCOS(DBLE(DELPHI))   
        IF (SIN_SCAT .EQ. 0.0) THEN   
          SIN2 = 0.0D0
          COS2 = -1.0D0
        ELSE
          SIN2 = SIN_THETA1*SINPHI /SIN_SCAT
          COS2 = (SIN_THETA2*SOLARMU - SIN_THETA1*MU*COSPHI)/SIN_SCAT
        ENDIF
        SIN22 = 2.0D0*SIN2*COS2
        COS22 = 1.0D0 - 2.0D0*SIN2**2
        SCATVECT(2) = B1*COS22
        SCATVECT(3) = B1*SIN22
      ENDIF
      IF (NSTOKES .EQ. 4) THEN
        SCATVECT(4) = 0.0
      ENDIF
      RETURN
      END







      SUBROUTINE COMPUTE_RADIANCE (NSTOKES, NX, NY, NZ, NPTS, NCELLS,
     .             ML, MM, NLM, NSTLEG, NLEG, NUMPHASE, 
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, 
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,
     .             XGRID, YGRID, ZGRID, GRIDPOS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,
     .             SHPTR, SOURCE, SOURCE1, GRIDRAD,
     .             OUTPARMS,  IRAD, RADOUT)
C       Computes the radiances for the specified locations and directions.
      IMPLICIT NONE
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
      INTEGER ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), IRAD
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SHPTR(NPTS+1), BCPTR(MAXNBC,2)
      INTEGER*2 CELLFLAGS(NCELLS)
      INTEGER IPHASE(NPTS)
      LOGICAL DELTAM
      REAL    SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
      REAL    MU(NMU), PHI(NMU,NPHI0MAX), WTDO(NMU,NPHI0MAX)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    SFCGRIDPARMS(*), BCRAD(NSTOKES,*)
      REAL    EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG,0:NLEG,NPTS)
      REAL    DIRFLUX(NPTS), FLUXES(2,NPTS)
      REAL    SOURCE1(NSTOKES,NPTS), GRIDRAD(NSTOKES,NPTS)
      REAL    SOURCE(NSTOKES,*)
      REAL    OUTPARMS(*), RADOUT(NSTOKES,*)
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
 
      INTEGER I, IBC, IANG, JX, JY, K, SIDE
      INTEGER NANGOUT, NXOUT, NYOUT
      LOGICAL LAMBERTIAN, VALIDRAD
      DOUBLE PRECISION X0,Y0,Z0, XE,YE,ZE, TRANSMIT, RADIANCE(4)
      REAL    MUOUT, PHIOUT, PHID
      REAL    XDOMAIN, YDOMAIN, STARTX, STARTY


C         Make the isotropic radiances for the top boundary
      CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN, 
     .                            UNITS, NTOPPTS, NSTOKES, BCRAD(1,1))
C         Make the bottom boundary radiances for the Lambertian surfaces.  
C          Compute the upwelling bottom radiances using the downwelling fluxes.
      IF (SFCTYPE .EQ. 'FL') THEN
        CALL FIXED_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2),
     .             DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO, 
     .             WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1,1+NTOPPTS))
      ELSE IF (SFCTYPE .EQ. 'VL') THEN
        CALL VARIABLE_LAMBERTIAN_BOUNDARY (NBOTPTS, BCPTR(1,2),
     .               DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,
     .               NSTOKES, BCRAD(1,1+NTOPPTS))
      ENDIF


C         Setup the regularly spaced radiance locations
      XDOMAIN = XGRID(NX+1-IBITS(BCFLAG,0,1))-XGRID(1) + 2*OUTPARMS(4)
      STARTX = -OUTPARMS(4)
      IF (OUTPARMS(2) .EQ. 0.0) THEN
        NXOUT = 1
      ELSE
        NXOUT = MAX(1,NINT(XDOMAIN/OUTPARMS(2)))
      ENDIF
      YDOMAIN = YGRID(NY+1-IBITS(BCFLAG,1,1))-YGRID(1) + 2*OUTPARMS(5)
      STARTY = -OUTPARMS(5)
      IF (OUTPARMS(3) .EQ. 0.0) THEN
        NYOUT = 1
      ELSE
        NYOUT = MAX(1,NINT(YDOMAIN/OUTPARMS(3)))
      ENDIF
      Z0 = MIN( MAX(OUTPARMS(1),ZGRID(1)), ZGRID(NZ))
      NANGOUT = NINT(OUTPARMS(6))


      LAMBERTIAN = SFCTYPE(2:2) .EQ. 'L'

C         Loop over the radiance directions
      DO IANG = 1, NANGOUT
        MUOUT = OUTPARMS(2*IANG+5)
        PHID = OUTPARMS(2*IANG+6)
        PHIOUT = PHID*ACOS(-1.0)/180.0
        IF (MUOUT .EQ. 0.0 .OR. ABS(MUOUT) .GT. 1.0) THEN
          WRITE (6,*) 'COMPUTE_RADIANCE: Bad mu for radiance',MUOUT
        ELSE          

C             Compute the source function throughout grid for this angle
          CALL COMPUTE_ONE_SOURCE (NSTOKES, ML, MM, NLM, NSTLEG,NLEG,
     .                 NUMPHASE, NPTS, DELTAM, MUOUT, PHIOUT, SRCTYPE,
     .                 SOLARMU, SOLARAZ, ALBEDO, LEGEN,
     .                 IPHASE, DIRFLUX, SHPTR, SOURCE,  SOURCE1)

C             Set the radiance field to -1, so we can determine valid radiances
C               Also set the source array to the extinction times the source
C               function, so the grid interpolation is correct.
          DO I = 1, NPTS
            GRIDRAD(1,I) = -1.0
            SOURCE1(:,I) = SOURCE1(:,I)*EXTINCT(I)
          ENDDO
C             Get boundary radiances: either top or bottom
C             Isotropic top boundary or Lambertian bottom boundary can use
C               the previously computed boundary radiances in BCRAD,
C               otherwise, compute the radiance for this angle
C               by integrating over the stored downwelling radiances.
          IF (MUOUT .LT. 0.0) THEN
            DO IBC = 1, NTOPPTS
              I = BCPTR(IBC,1)
              GRIDRAD(:,I) = BCRAD(:,IBC)
            ENDDO
          ELSE
            IF (.NOT. LAMBERTIAN) THEN
              CALL VARIABLE_BRDF_SURFACE (NBOTPTS,1,NBOTPTS,BCPTR(1,2),
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MUOUT, PHIOUT,
     .             SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX, 
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES, 
     .             BCRAD(1,1+NTOPPTS))
            ENDIF
            DO IBC = 1, NBOTPTS
              I = BCPTR(IBC,2)
              GRIDRAD(:,I) = BCRAD(:,NTOPPTS+IBC)
            ENDDO
          ENDIF

C             Integrate backward from the location to get the radiance
          Y0 = STARTY
          DO JY = 1, NYOUT
            X0 = STARTX
            DO JX = 1, NXOUT
              IRAD = IRAD + 1
              TRANSMIT = 1.0D0 ; RADIANCE(:) = 0.0D0
              CALL INTEGRATE_SOURCE (BCFLAG, IPFLAG, NSTOKES,
     .                        NX, NY, NZ, NPTS, NCELLS, 
     .                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                        XGRID, YGRID, ZGRID, GRIDPOS,
     .                        MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1,
     .                        X0, Y0, Z0, XE,YE,ZE, SIDE,
     .                        TRANSMIT, RADIANCE, VALIDRAD)
              RADOUT(:,IRAD) = RADIANCE(1:NSTOKES)
              X0 = X0 + OUTPARMS(2)
            ENDDO
            Y0 = Y0 + OUTPARMS(3)
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END


 
 
      SUBROUTINE COMPUTE_ONE_SOURCE (NSTOKES, 
     .             ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,
     .             NPTS, DELTAM, MU, PHI, SRCTYPE, SOLARMU, SOLARAZ,
     .             ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .             SHPTR, SOURCE,  SOURCE1)
C       Computes the source function (SOURCE1) in the direction (MU,PHI)
C     for the whole domain (NX,NZ).  The spherical harmonic source function
C     series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
      IMPLICIT NONE
      INTEGER NSTOKES, NPTS, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
      INTEGER SHPTR(NPTS+1)
      INTEGER*2 IPHASE(NPTS)
      LOGICAL DELTAM
      REAL    SOLARMU, SOLARAZ, MU, PHI
      REAL    ALBEDO(NPTS), LEGEN(NSTLEG,0:NLEG,NPTS), DIRFLUX(NPTS)
      REAL    SOURCE(NSTOKES,*), SOURCE1(NSTOKES,NPTS)
      CHARACTER SRCTYPE*1
      INTEGER I, J, JT, K, L, M, N, IPH, IS, NS, ME, MS
      REAL    OFOURPI, DA, F, A, A1, B1, SECMU0, COSSCAT
      REAL, ALLOCATABLE :: YLMDIR(:,:), YLMSUN(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: SUNDIRLEG(:)
      REAL, ALLOCATABLE :: UNSCLEGEN(:,:), SINGSCAT(:,:)
      LOGICAL TMS

      ALLOCATE (YLMDIR(NSTLEG,NLM), YLMSUN(NSTLEG,NLM))
      ALLOCATE (SUNDIRLEG(0:NLEG))
      ALLOCATE (UNSCLEGEN(2,0:NLEG), SINGSCAT(NSTOKES,NUMPHASE))
      IF (ABS(MU) .GT. 1.0) STOP 'COMPUTE_ONE_SOURCE: Bad mu'

      TMS = DELTAM
c      TMS = .FALSE.

C         Precompute Ylm's for output direction and solar direction
      CALL YLMALL (.FALSE., MU, PHI, ML, MM, NSTLEG, YLMDIR)
      IF (SRCTYPE .NE. 'T') THEN
        OFOURPI = 1.0/(4.0*ACOS(-1.0))
        CALL YLMALL (.TRUE., SOLARMU, SOLARAZ, ML, MM, NSTLEG,
     .               YLMSUN)
        SECMU0 = 1.0/ABS(SOLARMU)
        COSSCAT = SOLARMU*MU + SQRT((1.0D0-SOLARMU**2)*(1.0D0-MU**2))
     .              *COS(SOLARAZ-PHI)
        IF (NUMPHASE .EQ. 0 .AND. TMS) THEN
C           If non-tabulate phase functions (unpolarized) then
C           compute the Legendre polynomials for the scattering angle 
C           for the untruncated solar single scattering computation.
          CALL LEGENDRE_ALL (DBLE(COSSCAT), NLEG, SUNDIRLEG)
          DO L = 0, NLEG
            SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L+1)*OFOURPI
          ENDDO
        ENDIF
      ENDIF

C         If solar source and tabulated phase functions then precompute
C         the single scattered solar Stokes vector to outgoing direction
      IF (SRCTYPE .NE. 'T' .AND. TMS .AND. NUMPHASE .GT. 0) THEN
        DO IPH = 1, NUMPHASE
C           Unscaled the needed Legendre coefficients, but divide by 1-f 
C           because extinction is still scaled (for TMS method)
          F = LEGEN(1,ML+1,IPH)
          DO L = 0, NLEG
            IF (L .LE. ML) THEN
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH) + F/(1-F)
            ELSE
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)/(1-F)
            ENDIF
            IF (NSTLEG .GT. 1) THEN
              IF (L .LE. ML) THEN
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)
              ELSE
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)/(1-F)
              ENDIF
            ENDIF
          ENDDO
C           Sum the Wigner series (for A1 and B1), rotate the plane of 
C           polarization to the outgoing direction, and return the Stokes vector
          CALL SOLAR_SINGSCAT_VECT (NLEG, UNSCLEGEN, COSSCAT, 
     .                              SOLARMU, MU, SOLARAZ-PHI, NSTOKES,
     .                              SINGSCAT(1,IPH))
          IF (SINGSCAT(1,IPH) .LE. 0.0) THEN
            WRITE (6,*) 'COMPUTE_ONE_SOURCE: negative source function',
     .      ' for tabulated phase function:',IPH,MU,PHI,SINGSCAT(1,IPH)
            STOP
          ENDIF
        ENDDO
      ENDIF


C         Loop over all the grid points, computing the source function 
C           at the viewing angle from the spherical harmonic source function.
      DO I = 1, NPTS
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

C           Special case for solar source and Delta-M
        IF (SRCTYPE .NE. 'T' .AND. TMS) THEN
          IF (NUMPHASE .GT. 0) THEN
            IPH = IPHASE(I)
          ELSE
            IPH = I
          ENDIF
C             First subtract off the truncated single scattering 
          DA = ALBEDO(I)*DIRFLUX(I)*SECMU0
          J = 1
          DO L = 0, ML
            ME = MIN(L,MM)
            MS = -ME
            A1 = DA*LEGEN(1,L,IPH)
            B1 = DA*LEGEN(5,L,IPH)
            IF (J .LE. NS) THEN
              JT = J
              DO M = MS, ME
                SOURCE1(1,I) = SOURCE1(1,I) -A1*YLMDIR(1,J)*YLMSUN(1,J)
                J = J + 1
              ENDDO
              IF (NSTOKES .GT. 1) THEN
                J = JT
                DO M = MS, ME
                  SOURCE1(2,I)=SOURCE1(2,I) -B1*YLMDIR(2,J)*YLMSUN(1,J)
                  SOURCE1(3,I)=SOURCE1(3,I) -B1*YLMDIR(6,J)*YLMSUN(1,J)
                  J = J + 1
                ENDDO
              ENDIF
            ENDIF  
          ENDDO    

C             Then add in the single scattering contribution for the
C             original unscaled phase matrix.
          IF (NUMPHASE .GT. 0) THEN
            SOURCE1(:,I) = SOURCE1(:,I) + DA*SINGSCAT(:,IPH)
          ELSE IF (NSTOKES .EQ. 1) THEN
C             If we have a phase function for each grid point (which
C             should be the unpolarized case) then unscale the Legendre 
C             coefficients for L<=ML and sum over the Legendre series
C             in cosine of the scattering angle.
            F = LEGEN(1,ML+1,IPH)
            DO L = 0, NLEG
              IF (L .LE. ML) THEN
                A = DA*(LEGEN(1,L,IPH) + F/(1-F))
              ELSE
                A = DA*LEGEN(1,L,IPH)/(1-F)
              ENDIF
              SOURCE1(1,I) = SOURCE1(1,I) + A*SUNDIRLEG(L)
            ENDDO
          ENDIF
        ENDIF
      ENDDO
 
      DEALLOCATE (YLMDIR, YLMSUN, SUNDIRLEG, UNSCLEGEN, SINGSCAT)
      RETURN
      END
 


      SUBROUTINE SOLAR_SINGSCAT_VECT (NLEG, SINGLEGEN, COSSCAT, 
     .                               SOLARMU, MU, DELPHI, NSTOKES,
     .                               SCATVECT)
C      Calculates the Stokes vector for a single scattering of unpolarized 
C      sunlight from the solar direction (SOLARMU) to the outgoing 
C      direction (MU) and with cosine scattering angle (COSSCAT).  
C      The A1 and B1 elements of the phase matrix are calculated from 
C      the Wigner function expansion coefficients (SINGLEGEN). 
C      The plane of polarization is rotated from the scattering plane to 
C      the outgoing meridional plane.  Also divides the output vector by 4*pi.
      IMPLICIT NONE
      INTEGER NLEG, NSTOKES
      REAL    SINGLEGEN(2,0:NLEG), COSSCAT, SOLARMU, MU, DELPHI
      REAL    SCATVECT(NSTOKES)
      INTEGER I, J, K, L
      DOUBLE PRECISION  OFOURPI, X, FCT, A1, B1
      DOUBLE PRECISION, ALLOCATABLE :: DMM1(:)
      DOUBLE PRECISION SIN_SCAT, SIN_THETA1, SIN_THETA2
      DOUBLE PRECISION SINPHI, COSPHI, SIN2, COS2, SIN22, COS22
      PARAMETER (OFOURPI = 0.07957747151)


      ALLOCATE (DMM1(0:NLEG))

C       Sum the first Wigner function series, which is actually a Legendre series
      X = DBLE(COSSCAT)
      CALL WIGNERFCT (X, NLEG, 0, 0, DMM1)
      A1 = 0.0D0
      DO L = 0, NLEG
        FCT = 2.0D0*L + 1.0D0
        A1  = A1 + FCT*DBLE(SINGLEGEN(1,L))*DMM1(L)
      ENDDO
      SCATVECT(1) = A1*OFOURPI

      IF (NSTOKES .GT. 1) THEN
C        If doing polarization, sum the second Wigner function series
        CALL WIGNERFCT (X, NLEG, 2, 0, DMM1)
        B1 = 0.0D0
        DO L = 0, NLEG
          FCT = 2.0D0*L + 1.0D0
          B1  = B1 - FCT*DBLE(SINGLEGEN(2,L))*DMM1(L)
        ENDDO

C        Rotate the plane of polarization from the scattering plane 
C        to the outgoing meridional plane.
        SIN_SCAT = DSQRT(MAX(0.0D0,1.D0-COSSCAT**2))
        SIN_THETA1 = DSQRT(1.D0-SOLARMU**2)
        SIN_THETA2 = DSQRT(1.D0-MU**2)
        SINPHI = DSIN(DBLE(DELPHI))
        COSPHI = DCOS(DBLE(DELPHI))
        IF (SIN_SCAT .EQ. 0.0) THEN
          SIN2 = 0.0D0
          COS2 = -1.0D0
        ELSE
          SIN2 = SIN_THETA1*SINPHI /SIN_SCAT
          COS2 = (SIN_THETA2*SOLARMU - SIN_THETA1*MU*COSPHI)/SIN_SCAT
        ENDIF
        SIN22 = 2.0D0*SIN2*COS2
        COS22 = 1.0D0 - 2.0D0*SIN2**2
        SCATVECT(2) = B1*COS22*OFOURPI
        SCATVECT(3) = B1*SIN22*OFOURPI
      ENDIF
      IF (NSTOKES .EQ. 4) THEN
        SCATVECT(4) = 0.0
      ENDIF
      DEALLOCATE (DMM1)
      RETURN
      END




      SUBROUTINE INTEGRATE_SOURCE (BCFLAG, IPFLAG, NSTOKES,
     .                        NX, NY, NZ, NPTS, NCELLS, 
     .                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                        XGRID, YGRID, ZGRID, GRIDPOS,
     .                        MU, PHI, GRIDRAD, EXTINCT, SOURCE, 
     .                        X0, Y0, Z0,  XE, YE, ZE, SIDE, 
     .                        TRANSMIT, RADIANCE, VALIDRAD)
C       Integrates the source function (SOURCE) through the extinction
C     field (EXTINCT) backward in the direction opposite to (MU,PHI) 
C     to find the outgoing radiance (RADIANCE) at the point X0,Y0,Z0.
C     The radiances at the opposite boundary are input in GRIDRAD.
C     The transmission and radiance of the ray so far (TRANSMIT, RADIANCE)
C     are input and returned after the integration along with the exitting
C     ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
C     5=-Z,6=+Z).
      IMPLICIT NONE
      INTEGER BCFLAG, IPFLAG, NSTOKES, NX, NY, NZ, NPTS, NCELLS, SIDE
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER*2 CELLFLAGS(NCELLS)
      LOGICAL VALIDRAD
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    MU, PHI
      DOUBLE PRECISION X0,Y0,Z0, XE,YE,ZE, TRANSMIT, RADIANCE(NSTOKES)
      REAL    EXTINCT(NPTS)
      REAL    SOURCE(NSTOKES,NPTS), GRIDRAD(NSTOKES,NPTS)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER I1, I2, I3, I4, IOPP, NTAU, IT, K
      LOGICAL IPINX, IPINY, OPENBCFACE, BTEST
      INTEGER GRIDFACE(4,6), OPPFACE(6), JFACE, KFACE, IC
      REAL    EXT0, EXT1, EXTN, XM, YM
      REAL    SRCEXT0(NSTOKES), SRCEXT1(NSTOKES), RAD0(NSTOKES)
      DOUBLE PRECISION PI, CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
      DOUBLE PRECISION TAUTOL, TAUGRID, S, DELS
      DOUBLE PRECISION EXT, TAU, TRANSCELL,ABSCELL
      DOUBLE PRECISION SRC(NSTOKES), RAD(NSTOKES)
      DATA GRIDFACE/1,3,5,7, 2,4,6,8,  1,2,5,6, 3,4,7,8,
     .              1,2,3,4, 5,6,7,8/, OPPFACE/2,1,4,3,6,5/

C         TAUTOL is the maximum optical path for the subgrid intervals
      TAUTOL = 0.1
      EPS = 1.0E-4*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))

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
      IF (CZ .LT. -1.0E-3) THEN
        BITZ = 1
      ELSE IF (CZ .GT. 1.0E-3) THEN
        BITZ = 0
      ELSE 
        STOP 'INTEGRATE_SOURCE: Bad MU'
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
     .                  BCFLAG, IPFLAG, XE, YE, ZE,  ICELL)
c      print '(A,5(1X,F8.5),1X,I6)', 
c     .     'INTEGRATE_SOURCE:',X0,Y0,Z0,MU,PHI,ICELL
      CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                        XE, YE, ZE, 1, EXTINCT, 1, 1, EXT1)
      CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                XE, YE, ZE, NSTOKES, SOURCE, 1, NSTOKES, SRCEXT1)
      SRCEXT1(1) = MAX(0.0,SRCEXT1(1))

C         Loop until finding a face with known radiances
      VALIDRAD = .FALSE.
      DO WHILE (.NOT. VALIDRAD .AND. ICELL .GT. 0)
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
        IF (IPINX .OR. CX .EQ. 0.0) THEN
          SOX = 1.0E20
        ELSE
          SOX = (GRIDPOS(1,IOPP)-XE)*CXINV
        ENDIF
        IF (IPINY .OR. CY .EQ. 0.0) THEN
          SOY = 1.0E20
        ELSE
          SOY = (GRIDPOS(2,IOPP)-YE)*CYINV
        ENDIF
        SOZ = (GRIDPOS(3,IOPP)-ZE)*CZINV
        SO = MIN(SOX,SOY,SOZ)
        IF (SO .LT. -EPS) THEN
          WRITE (6,*) 'INTEGRATE_SOURCE: SO<0  ', 
     .      MU,PHI,X0,Y0,Z0,XE,YE,ZE,SO,SOX,SOY,SOZ,XGRID(1),YGRID(1),
     .      GRIDPOS(1,IOPP),GRIDPOS(2,IOPP),GRIDPOS(3,IOPP)
          STOP
        ENDIF
        IF (IPINX) THEN
          XN = XE
        ELSE
          XN = XE + SO*CX
        ENDIF
        IF (IPINY) THEN
          YN = YE
        ELSE
          YN = YE + SO*CY
        ENDIF
        ZN = ZE + SO*CZ

C           Find the optical path across the grid cell and figure how
C             many subgrid intervals to use
        CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                          XN, YN, ZN, 1, EXTINCT, 1, 1, EXTN)
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
          IF (IT .EQ. NTAU) THEN
            EXT0 = EXTN
          ELSE
            CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                          XI, YI, ZI, 1, EXTINCT, 1, 1, EXT0)
          ENDIF
          CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .               XI, YI, ZI, NSTOKES, SOURCE, 1, NSTOKES, SRCEXT0)
          SRCEXT0(1) = MAX(0.0,SRCEXT0(1))

C            Compute the subgrid radiance: integration of the source function
          EXT = 0.5*(EXT0+EXT1)
          IF (EXT .NE. 0.0) THEN
            TAU=EXT*DELS
            ABSCELL = TAU*(1.0-0.5*TAU*(1.0-0.33333333333*TAU))
            TRANSCELL = 1.0 - ABSCELL
C                 Linear extinction, linear source*extinction, to second order
            SRC(:) = ( 0.5*(SRCEXT0(:)+SRCEXT1(:)) + 0.08333333333*
     .                   (EXT0*SRCEXT1(:)-EXT1*SRCEXT0(:))*DELS
     .                    *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
            SRC(:) = 0.0
          ENDIF

          RADIANCE(:) = RADIANCE(:) + TRANSMIT*SRC(:)*ABSCELL
          SRCEXT1(:) = SRCEXT0(:)
          EXT1 = EXT0
          TRANSMIT = TRANSMIT*TRANSCELL
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
C           Get the four grid points on the intersection face for
C             telling if at the boundary
C             If going to same or larger face then use previous face
        IF (NEIGHPTR(IFACE,ICELL) .GE. 0 .AND. .NOT.OPENBCFACE) THEN
          KFACE = IFACE
          IC = ICELL
        ELSE
C             If going to smaller face then use next face (more accurate)
          KFACE = OPPFACE(IFACE)
          IC = INEXTCELL
        ENDIF
        I1 = GRIDPTR(GRIDFACE(1,KFACE),IC)
        I2 = GRIDPTR(GRIDFACE(2,KFACE),IC)
        I3 = GRIDPTR(GRIDFACE(3,KFACE),IC)
        I4 = GRIDPTR(GRIDFACE(4,KFACE),IC)
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

C           See if there are valid radiances for the intersection face:
C           If so, set flag, interpolate radiance, and add in 
C           If not, only prepare for next cell
        IF (TRANSMIT .LT. 1.0E-5) THEN
          VALIDRAD = .TRUE.
        ELSE IF (GRIDRAD(1,I1).GE.0.0 .AND. GRIDRAD(1,I2).GE.0.0 .AND.
     .      GRIDRAD(1,I3) .GE. 0.0 .AND. GRIDRAD(1,I4) .GE. 0.0) THEN
          VALIDRAD = .TRUE.
          CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, XI, YI, ZI,
     .                            NSTOKES, GRIDRAD, 1, NSTOKES, RAD0)
          RADIANCE(:) = RADIANCE(:) + TRANSMIT*RAD0(:)
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




      SUBROUTINE INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                            X, Y, Z, ND, FIELD, ID1, NV, VALUES)
C       Trilinearly interpolates the field (FIELD) given a grid cell
C     pointer (ICELL), to compute the field values (VALUES) at the desired
C     point (X,Y,Z). ND is the leading index of the FIELD (1 for a 
C     scalar field), ID1 is the first element of the vector to interpolate,
C     and NV elements are output.
      IMPLICIT NONE
      INTEGER ICELL, GRIDPTR(8,*), ND, ID1, NV
      DOUBLE PRECISION X, Y, Z
      REAL    GRIDPOS(3,*), FIELD(ND,*), VALUES(NV)
      INTEGER IPT1, IPT2, ID
      DOUBLE PRECISION U, V, W, UM, VM, WM, DELX, DELY, DELZ

      IPT1 = GRIDPTR(1,ICELL)
      IPT2 = GRIDPTR(8,ICELL)
      DELX = GRIDPOS(1,IPT2) - GRIDPOS(1,IPT1)
      IF (DELX .LE. 0.0)  DELX = 1.0
      DELY = GRIDPOS(2,IPT2) - GRIDPOS(2,IPT1)
      IF (DELY .LE. 0.0)  DELY = 1.0
      DELZ = GRIDPOS(3,IPT2) - GRIDPOS(3,IPT1)
      U = (X-GRIDPOS(1,IPT1))/DELX
      V = (Y-GRIDPOS(2,IPT1))/DELY
      W = (Z-GRIDPOS(3,IPT1))/DELZ
      UM = 1-U
      VM = 1-V
      WM = 1-W

      DO ID = ID1, ID1+NV-1      
        VALUES(ID-ID1+1) = WM*(VM*(UM*FIELD(ID,GRIDPTR(1,ICELL))
     .               + U*FIELD(ID,GRIDPTR(2,ICELL)))
     .           + V*(UM*FIELD(ID,GRIDPTR(3,ICELL))
     .               + U*FIELD(ID,GRIDPTR(4,ICELL)))) 
     .       + W*(VM*(UM*FIELD(ID,GRIDPTR(5,ICELL))
     .               + U*FIELD(ID,GRIDPTR(6,ICELL)))
     .           + V*(UM*FIELD(ID,GRIDPTR(7,ICELL))
     .               + U*FIELD(ID,GRIDPTR(8,ICELL))))
      ENDDO
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





      SUBROUTINE LEGENDRE_ALL (COSSCAT, NLEG, P)
C       This subroutine computes a set of Legendre polynomials for
C     a particular scattering angle COSSCAT.  NLEG is the maximum term.
C     The Legendre functions evaluated at COSSCAT are returned in 
C     P, starting at l=0 and ending with l=NLEG  (NLEG+1 terms).
      IMPLICIT NONE
      INTEGER NLEG
      DOUBLE PRECISION COSSCAT, P(0:NLEG)
      INTEGER L
      DOUBLE PRECISION X, PL, PL1, PL2

      X = DBLE(COSSCAT)
      IF (X*X .GT. 1.) STOP 'LEGENDRE_ALL: |COSSCAT| larger than 1'
C         Use the stable upward recursion on l, starting from P_0
      PL2 = 1.0D0
      P(0) = PL2
      IF (NLEG .GT. 1) THEN
        PL1 = X
        P(1) = X
      ENDIF
      DO L = 2, NLEG
        PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
        P(L) = PL
        PL2 = PL1
        PL1 = PL
      ENDDO
      RETURN
      END




      SUBROUTINE YLMALL (TRANSPOSE, MU, PHI, ML, MM, NSTLEG, YR)
C       This subroutine computes a set of normalized real generalized
C     spherical harmonic functions (YR) for a particular direction mu,phi.
C     (Doicu et al. 2013; http://dx.doi.org/10.1016/j.jqsrt.2012.12.009)
C     ML is the maximum meridional mode, MM is the maximum azimuthal mode.
C     NSTLEG is the number of matrix elements returned, either 1 for the
C     1,1 term or 6 for the 6 unique elements of the 4x4 matrix.
C     J = 2*(L*(L+1))/2 + M+1  for L<=MM
C     J = (2*MM+1)*L-MM*(2+2*(MM-1))/2 + M+1  for L>MM
      IMPLICIT NONE
      INTEGER ML, MM, NSTLEG
Cf2py intent(in) :: ML, MM, NSTLEG
      REAL    MU, PHI
Cf2py intent(in) :: MU, PHI
      LOGICAL TRANSPOSE
Cf2py intent(in) :: TRANSPOSE
      REAL    YR(NSTLEG,*)
Cf2py intent(in,out) :: YR
      INTEGER           J, M, MABS, L
      DOUBLE PRECISION  X, PI, FCT, P1, P2, P3, COSM, SINM, SIGN
      DOUBLE PRECISION, ALLOCATABLE :: DM0(:), DM2P(:), DM2M(:)


      IF (NSTLEG .EQ. 1) THEN
        CALL YLMALL_UNPOL (MU, PHI, ML, MM, YR)
        RETURN
      ENDIF

      ALLOCATE (DM0(0:ML), DM2P(0:ML), DM2M(0:ML))
      X = DBLE(MU)

      PI  = DACOS(- 1.D0)
      FCT = 1.D0/DSQRT(2.D0*PI)

C       sign accounting for transpose matrices
      IF (.NOT. TRANSPOSE) THEN
        SIGN =  1.D0
      ELSE
        SIGN = -1.D0
      END IF

C     ..................................................................
C                                M = 0
      M = 0
      CALL WIGNERFCT02P2M_NORMALIZED (X, ML, M,  DM0, DM2P, DM2M)
      DO L = 0, ML
        IF (L .LE. MM) THEN
          J = L*(L+1) + M + 1
        ELSE
          J = (2*MM+1)*L - MM**2 + M + 1
        ENDIF
        P1 = FCT*DM0(L)
        P2 = -0.5D0*FCT*(DM2P(L) + DM2M(L))
        P3 = -0.5D0*FCT*(DM2P(L) - DM2M(L))

C       --- Unm = Pnm*cos(m*phi) ---
        YR(1,J) =  P1
        IF (NSTLEG .EQ. 6) THEN
          YR(2,J) = P2
          YR(3,J) = P2
          YR(4,J) = P1
          YR(5,J) = P3
          YR(6,J) = P3
        ENDIF
      END DO

C     ..................................................................
C                                   M =/ 0
      DO MABS = 1, MM
        CALL WIGNERFCT02P2M_NORMALIZED (X, ML, MABS,DM0, DM2P, DM2M)
        COSM = COS(MABS*PHI)
        SINM = SIN(MABS*PHI)
        DO L = MABS, ML
C         ..............................................................
C                                     M > 0
          M = MABS
          IF (L .LE. MM) THEN
            J = L * (L+1) + M + 1
          ELSE
            J = (2*MM+1) * L - MM**2 + M + 1
          ENDIF
          P1 = FCT*DM0(L)
          P2 = -0.5D0*FCT*(DM2P(L) + DM2M(L))
          P3 = -0.5D0*FCT*(DM2P(L) - DM2M(L))
C
C         --- Unm = Pnm*cos(m*phi) - D*Pnm*sin(m*phi) ---
          YR(1,J) = P1*COSM - P1*SINM
          IF (NSTLEG .EQ. 6) THEN
            YR(2,J) = P2*COSM - P2*SINM
            YR(3,J) = P2*COSM + P2*SINM
            YR(4,J) = P1*COSM + P1*SINM
            YR(5,J) = P3*COSM - SIGN*P3*SINM
            YR(6,J) = P3*COSM + SIGN*P3*SINM
          ENDIF
C         ..............................................................
C                                      M < 0
          M = -MABS
          IF (L .LE. MM) THEN
            J = L*(L+1) + M + 1
          ELSE
            J = (2*MM + 1)*L - MM**2 + M + 1
          ENDIF

C         --- Vnm = Pnm*sin(m*phi) + D*Pnm*cos(m*phi) ---
          YR(1,J) = P1*SINM + P1*COSM
          IF (NSTLEG .EQ. 6) THEN
            YR(2,J) = P2*SINM + P2*COSM
            YR(3,J) = P2*SINM - P2*COSM
            YR(4,J) = P1*SINM - P1*COSM
            YR(5,J) = P3*SINM + SIGN*P3*COSM
            YR(6,J) = P3*SINM - SIGN*P3*COSM
          ENDIF
        END DO
      END DO
      DEALLOCATE (DM0, DM2P, DM2M)
      RETURN
      END


      SUBROUTINE WIGNERFCT02P2M_NORMALIZED (MIU, NRANK, M, 
     .                                      DM0, DM2P, DM2M)
C       -----------------------------------------------------------------
C       |  The routine computes the DMM1N vector coefficients for       |
C       |               M >= 0; M1 = 0, 2, -2; N = 0,...,NRANK, and     |
C       |              -1 <= MIU <= 1                                   |
C       |  For M1 = -2, we use the symmetry relation                    |
C       |               DM-M1(MIU) = (-1)**(N+M) * DMM1(-MIU)           |
C       -----------------------------------------------------------------
        IMPLICIT NONE
        INTEGER           NRANK, M
        DOUBLE PRECISION  MIU
        DOUBLE PRECISION  DM0(0:NRANK), DM2P(0:NRANK), DM2M(0:NRANK)
        INTEGER           N0, N
        DOUBLE PRECISION  XP, XM, FACT1, FACT2, FACTP, FACTM, DMM1_N0
C
        XP =   MIU !for M1 = 0 and M1 = 2
        XM = - MIU !for M1 = - 2
C       ................................................................
C                                DM0 for M1 = 0
C       ................................................................
        DM0 = 0.D0
        IF( M .LE. NRANK ) THEN
          IF ( M .EQ. 0 ) THEN
            DM0(0) = 1.D0
            DM0(1) = XP
            DO N = 1, NRANK - 1
              FACT1 = DBLE( 2 * N + 1 ) * XP / DBLE( N + 1 )
              FACT2 = DBLE( N ) / DBLE( N + 1 )
              DM0(N+1) = FACT1 * DM0(N) - FACT2 * DM0(N-1)
            END DO
          ELSE
            DM0(M) = DMM1_N0( XP, M, 0 )
            DO N = M, NRANK - 1
              FACT1 = DBLE( N * (N + 1) ) * XP
              FACT1 = FACT1 / DSQRT( DBLE( (N + 1)**2 - M**2   ) )
              FACT1 = FACT1 / DBLE( N + 1 )
              FACT1 = FACT1 * DBLE( 2 * N + 1 ) / DBLE( N )
C
              FACT2 = DSQRT( DBLE( N**2 - M**2   ) ) * DBLE( N )
              FACT2 = FACT2 / DSQRT( DBLE( (N + 1)**2 - M**2   ) )
              FACT2 = FACT2 / DBLE( N + 1 )
              FACT2 = FACT2 * DBLE( N + 1 ) / DBLE( N )
C
              DM0(N+1) = FACT1 * DM0(N) - FACT2 * DM0(N-1)
            END DO
          END IF
          DO N = 0, NRANK
            DM0(N) = DSQRT( N + 0.5D0 ) * DM0(N)
          END DO
        END IF
C       ................................................................
C                    DM2P and DM2M for M1 = 2 and M1 = -2
C       ................................................................
        N0   = MAX( M, 2 )
        DM2P = 0.D0
        DM2M = 0.D0
        IF( M .LE. NRANK .AND. NRANK .GE. 2 ) THEN
          DM2P(N0) = DMM1_N0( XP, M, 2 )
          DM2M(N0) = DMM1_N0( XM, M, 2 )
          DO N = N0, NRANK - 1
            FACTP = DBLE( N * (N + 1) ) * XP - DBLE( 2 * M  )
            FACTM = DBLE( N * (N + 1) ) * XM - DBLE( 2 * M  )
C
            FACT1 = 1.D0  / DSQRT( DBLE( (N + 1)**2 - M**2   ) )
            FACT1 = FACT1 / DSQRT( DBLE( (N + 1)**2 - 4 ) )
            FACT1 = FACT1 * DBLE( 2 * N + 1 ) / DBLE( N )
C
            FACT2 = DSQRT( DBLE( N**2 - M**2   ) )
     &            * DSQRT( DBLE( N**2 - 4 ) )
            FACT2 = FACT2 / DSQRT( DBLE( (N + 1)**2 - M**2   ) )
            FACT2 = FACT2 / DSQRT( DBLE( (N + 1)**2 - 4 ) )
            FACT2 = FACT2 * DBLE( N + 1 ) / DBLE( N )
C
            DM2P(N+1) = FACTP * FACT1 * DM2P(N) - FACT2 * DM2P(N-1)
            DM2M(N+1) = FACTM * FACT1 * DM2M(N) - FACT2 * DM2M(N-1)
          END DO
          DO N = 0, NRANK
            DM2M(N) = ( - 1.D0 )**( N + M ) * DM2M(N)
C
            DM2P(N) = DSQRT( N + 0.5D0 ) * DM2P(N)
            DM2M(N) = DSQRT( N + 0.5D0 ) * DM2M(N)
          END DO
        END IF
        RETURN
        END



        FUNCTION DMM1_N0( X, M, M1 ) RESULT( DMM1N0 )
C       ---------------------------------------------------------------
C       | THE ROUTINE COMPUTES THE WIGNER FUNCTIONS FOR               |
C       |     M >= 0, M1 >= 0, AND N0 = MAX(M,M1)                     |
C       ---------------------------------------------------------------
        IMPLICIT NONE
        INTEGER           M, M1
        DOUBLE PRECISION  X, DMM1N0
C
        INTEGER           P, MAXM, MINM
        DOUBLE PRECISION  PROD, FACT, CMM1
C
        IF ( M .EQ. M1 ) THEN
          FACT   = ( ( 1.D0 + X ) / 2.D0 )**M
          DMM1N0 = FACT
        ELSE
          IF ( M1 .GT. M ) THEN
            CMM1 = 1.D0
          ELSE
            CMM1 = ( - 1.D0 )**( M - M1 )
          END IF
          MAXM = MAX(M,M1)
          MINM = MIN(M,M1)
          PROD = 1.D0
          DO P = 1, MAXM - MINM
            FACT = DSQRT( DBLE( M + M1 + P ) / DBLE( P ) )
            PROD = PROD * FACT
          END DO
          FACT   = DSQRT( ( 1.D0 - X ) / 2.D0 )
          DMM1N0 = CMM1 * PROD * FACT**( MAXM - MINM )
          FACT   = DSQRT( ( 1.D0 + X ) / 2.D0 )
          DMM1N0 = DMM1N0 * FACT**( MAXM + MINM )
        END IF
        RETURN
        END



      SUBROUTINE YLMALL_UNPOL (MU, PHI, ML, MM, YR)
C       This subroutine computes a set of the unpolarized (I-I) element of
C     normalized real generalized spherical harmonic functions (YR) for 
C     a particular direction mu,phi.  (Doicu et al. 2013; 
C     http://dx.doi.org/10.1016/j.jqsrt.2012.12.009)
C     ML is the maximum meridional mode, MM is the maximum azimuthal mode.
C       J = 2*(L*(L+1))/2 + M+1  for L<=MM
C       J = (2*MM+1)*L-MM*(2+2*(MM-1))/2 + M+1  for L>MM
      IMPLICIT NONE
      INTEGER ML, MM
      REAL    MU, PHI
      REAL    YR(*)
      INTEGER           J, M, L
      DOUBLE PRECISION  X, PI, FCT, COSM, SINM
      DOUBLE PRECISION, ALLOCATABLE :: DM0(:)

      ALLOCATE (DM0(0:ML))
      X = DBLE(MU)

      PI  = DACOS(-1.D0)
      FCT = 1.D0/DSQRT(2.D0*PI)

      DO M = 0, MM
        CALL WIGNERFCT_DM0 (X, ML, M, DM0)
        DM0(:) = FCT*DM0(:)
        IF (M .GT. 0) THEN
          COSM = COS(M*PHI)
          SINM = SIN(M*PHI)
        ELSE
          COSM = 1.0D0
          SINM = 0.0D0
        ENDIF
        DO L = M, ML
          IF (L .LE. MM) THEN
            J = L*(L+1) + M+1
          ELSE
            J = (2*MM+1)*L - MM**2 + M+1
          ENDIF

C           M>=0:  Unm = Pnm*cos(m*phi) - D*Pnm*sin(m*phi)
          YR(J) = (COSM - SINM)*DM0(L)

          J = J - 2*M
C            M<0:  Vnm = Pnm*sin(m*phi) + D*Pnm*cos(m*phi)
          YR(J) = (COSM + SINM)*DM0(L)
        END DO
      END DO
      DEALLOCATE (DM0)
      RETURN
      END


      SUBROUTINE WIGNERFCT_DM0 (X, NRANK, M, DM0)
C       -----------------------------------------------------------------
C       |  The routine computes the DMM1N vector coefficients for       |
C       |               M >= 0; M1 = 0; N = 0,...,NRANK, and     |
C       |              -1 <= X <= 1                                   |
C       -----------------------------------------------------------------
        IMPLICIT NONE
        INTEGER           NRANK, M
        DOUBLE PRECISION  X
        DOUBLE PRECISION  DM0(0:NRANK)
        INTEGER           N
        DOUBLE PRECISION  DM_M10_N0

C                   DM0 for M1 = 0
        DM0(:) = 0.D0
        IF (M .LE. NRANK) THEN
          IF (M .EQ. 0) THEN
            DM0(0) = 1.0D0
            DM0(1) = X
            DO N = 1, NRANK-1
              DM0(N+1) = ( (2*N+1)*X*DM0(N) - N*DM0(N-1) )/(N+1)
            ENDDO
          ELSE
            DM0(M) = DM_M10_N0 (X, M)
            DO N = M, NRANK-1
              DM0(N+1) = ( (2*N+1)*X*DM0(N) 
     .                            -  DSQRT(DBLE(N**2-M**2))*DM0(N-1) )
     .                        /DSQRT(DBLE((N+1)**2-M**2)) 
            ENDDO
          END IF
          DO N = 0, NRANK
            DM0(N) = DSQRT(N+0.5D0)*DM0(N)
          ENDDO
        ENDIF
        RETURN
        END


        FUNCTION DM_M10_N0 (X, M)  RESULT(DMM0N0)
C        The routine computes the Wigner functions for M >= 0, M1 = 0, and
C        N0 = M          
        IMPLICIT NONE
        INTEGER           M
        DOUBLE PRECISION  X, DMM0N0
        INTEGER           P
        DOUBLE PRECISION  CM, PROD

        IF ( M .EQ. 0 ) THEN
          DMM0N0 = 1.D0
        ELSE
          CM = (-1.D0)**M
          PROD = 1.D0
          DO P = 1, M
            PROD = PROD * DSQRT((M+P)/DBLE(P))
          ENDDO
          DMM0N0 = CM*PROD* ( 0.5D0*DSQRT((1.D0-X)*(1.D0+X)) )**M
        ENDIF
        RETURN
        END



        SUBROUTINE WIGNERFCT( X, NRANK, M, M1, DMM1 )
C       -----------------------------------------------------------------
C       |  The routine computes the DMM1N vector coefficients for       |
C       |               M >= 0, M1 >= 0,  N = 0,...,NRANK, and          |
C       |               -1 < X= COS(BETA) < 1                           |
C       -----------------------------------------------------------------
        IMPLICIT NONE
        INTEGER           NRANK, M, M1
        DOUBLE PRECISION  X, DMM1( 0:NRANK )
C
        INTEGER           N0, N
        DOUBLE PRECISION  FACT1, FACT2, DMM1_N0
C
        N0   = MAX( M, M1 )
        DMM1 = 0.D0
        IF ( N0 .EQ. 0 ) THEN
          DMM1(0) = 1.D0
          DMM1(1) = X
          DO N = 1, NRANK - 1
            FACT1 = DBLE( 2 * N + 1 ) * X / DBLE( N + 1 )
            FACT2 = DBLE( N ) / DBLE( N + 1 )
            DMM1(N+1) = FACT1 * DMM1(N) - FACT2 * DMM1(N-1)
          END DO
        ELSE
          DMM1(N0) = DMM1_N0( X, M, M1 )
          DO N = N0, NRANK - 1
            FACT1 = DBLE( N * (N + 1) ) * X - DBLE( M * M1  )
            FACT1 = FACT1 / DSQRT( DBLE( (N + 1)**2 - M**2   ) )
            FACT1 = FACT1 / DSQRT( DBLE( (N + 1)**2 - M1**2 ) )
            FACT1 = FACT1 * DBLE( 2 * N + 1 ) / DBLE( N )
C
            FACT2 = DSQRT( DBLE( N**2 - M**2   ) ) 
     &             * DSQRT( DBLE( N**2 - M1**2 ) )
            FACT2 = FACT2 / DSQRT( DBLE( (N + 1)**2 - M**2   ) )
            FACT2 = FACT2 / DSQRT( DBLE( (N + 1)**2 - M1**2 ) )
            FACT2 = FACT2 * DBLE( N + 1 ) / DBLE( N )
C
            DMM1(N+1) = FACT1 * DMM1(N) - FACT2 * DMM1(N-1)
          END DO
        END IF
        RETURN
        END




      SUBROUTINE PLMALL (TRANSPOSE, MU, ML, MM, PRC)
      IMPLICIT NONE
      INTEGER ML, MM
      REAL    MU
      LOGICAL TRANSPOSE
      REAL    PRC(6,*)
C     --- local variables ---
      INTEGER           J, M, MABS, L
      DOUBLE PRECISION  X, PI, FCT, P1, P2, P3, SIGN
      DOUBLE PRECISION, ALLOCATABLE :: DM0( : ), DM2P( : ),
     &                  DM2M( : )

C     --- allocate ---
      ALLOCATE( DM0( 0:ML ), DM2P( 0:ML ), DM2M( 0:ML ) )

      X = DBLE( MU )

C     --- normalization constant ---
      PI  = DACOS( - 1.D0 )
      FCT = 1.D0 / DSQRT( 2.D0 * PI )

C     --- sign accounting for transpose matrices ---
      IF ( .NOT. TRANSPOSE ) THEN
        SIGN =  1.D0
      ELSE
        SIGN = -1.D0
      END IF
C     ..................................................................
C                                M = 0
C     ..................................................................
      M = 0
      CALL WIGNERFCT02P2M_NORMALIZED
     I   ( X, ML, M,
     O     DM0, DM2P, DM2M )
      DO L = 0, ML
        IF ( L .LE. MM ) THEN
          J = L * ( L + 1 ) + M + 1
        ELSE
          J = ( 2*MM + 1 ) * L - MM**2 + M + 1
        ENDIF
        P1 =   FCT * DM0(L)
        P2 = - 0.5D0 * FCT * ( DM2P(L) + DM2M(L) )
        P3 = - 0.5D0 * FCT * ( DM2P(L) - DM2M(L) )
C
C       --- PRCnm = Pnm ---
        PRC(1,J) =  P1
        PRC(2,J) =  P2
        PRC(3,J) =  P2
        PRC(4,J) =  P1
        PRC(5,J) =  P3
        PRC(6,J) =  P3
      END DO
C     ..................................................................
C                                   M =/ 0
C     ..................................................................
      DO MABS = 1, MM
        CALL WIGNERFCT02P2M_NORMALIZED
     I     ( X, ML, MABS,
     O       DM0, DM2P, DM2M )
        DO L = MABS, ML
C         ..............................................................
C                                     M > 0
C         ..............................................................
          M = MABS
          IF ( L .LE. MM ) THEN
            J = L * ( L + 1 ) + M + 1
          ELSE
            J = ( 2*MM + 1 ) * L - MM**2 + M + 1
          ENDIF
          P1 =   FCT * DM0(L)
          P2 = - 0.5D0 * FCT * ( DM2P(L) + DM2M(L) )
          P3 = - 0.5D0 * FCT * ( DM2P(L) - DM2M(L) )
C
C         --- PRCnm = Pnm ---
          PRC(1,J) =  P1
          PRC(2,J) =  P2
          PRC(3,J) =  P2
          PRC(4,J) =  P1
          PRC(5,J) =  P3
          PRC(6,J) =  P3
C         ..............................................................
C                                      M < 0
C         ..............................................................
          M = - MABS
          IF ( L .LE. MM ) THEN
            J = L * ( L + 1 ) + M + 1
          ELSE
            J = ( 2*MM + 1 ) * L - MM**2 + M + 1
          ENDIF
C
C         --- PRCnm = D*Pnm ---
          PRC(1,J) =  P1
          PRC(2,J) =  P2
          PRC(3,J) = -P2
          PRC(4,J) = -P1
          PRC(5,J) =  SIGN * P3
          PRC(6,J) = -SIGN * P3
        END DO
      END DO
      DEALLOCATE( DM0, DM2P, DM2M )
      RETURN
      END


 
 
      SUBROUTINE PLANCK_FUNCTION (TEMP, UNITS, WAVENO, WAVELEN, PLANCK)
C        Calculates the Planck blackbody radiance. If UNITS='T' then
C     using brightness temperature units and the temperature is simply
C     returned. If UNITS='B' then doing a band integration and the
C     Planck blackbody radiance in [Watts /(meter^2 ster)] over a 
C     wavenumber range [cm^-1] is returned. Otherwise, the Planck 
C     blackbody radiance in [Watts /(meter^2 ster micron)] for a 
C     temperature in [Kelvins] at a wavelength in [microns] is returned.
      IMPLICIT NONE
      REAL  TEMP, WAVENO(2), WAVELEN, PLANCK
      CHARACTER*1  UNITS
      DOUBLE PRECISION X1, X2, F
 
      IF (UNITS .EQ. 'T') THEN
        PLANCK = TEMP
      ELSE IF (UNITS .EQ. 'B') THEN
        IF (TEMP .GT. 0.0) THEN
          X1 = 1.4388D0*WAVENO(1)/TEMP
          X2 = 1.4388D0*WAVENO(2)/TEMP
          CALL INTEGRATE_PLANCK (X1, X2, F)
          PLANCK = 1.1911D-8*(TEMP/1.4388D0)**4 *F
        ELSE
          PLANCK = 0.0
        ENDIF
      ELSE
        IF (TEMP .GT. 0.0) THEN
          PLANCK = 1.1911E8 / WAVELEN**5
     .            / (EXP(1.4388E4/(WAVELEN*TEMP)) - 1)
        ELSE
          PLANCK = 0.0
        ENDIF
      ENDIF
      RETURN
      END
 


      SUBROUTINE INTEGRATE_PLANCK (X1, X2, F)
C       Returns integral of x**3/(exp(x)-1) from x1 to x2.
C       Accurate to better than 1 part in 10**9 for integral from
C       0 to infinity (pi**4/15).
      IMPLICIT NONE
      DOUBLE PRECISION X1, X2, F
      DOUBLE PRECISION INTLOW, INTHIGH
      DOUBLE PRECISION C

      C = 1.0D0
      IF (X1 .LT. C .AND. X2 .LT. C) THEN
        F = INTLOW(X2) - INTLOW(X1)
      ELSE IF (X1 .LT. C .AND. X2 .GE. C) THEN
        F = INTLOW(C) - INTLOW(X1) + INTHIGH(C) - INTHIGH(X2)
      ELSE IF (X1 .GE. C .AND. X2 .GE. C) THEN
        F = INTHIGH(X1) - INTHIGH(X2)
      ELSE
        STOP 'X1 and X2 out of order'
      ENDIF
      RETURN
      END


      DOUBLE PRECISION FUNCTION INTLOW (X)
C       Integral of x**3/(exp(x)-1) from 0 to x.
C       Accurate for x less than about 1.
C       Uses Taylor series expansion around x=0.
      IMPLICIT NONE
      DOUBLE PRECISION X
      INTEGER N
      DOUBLE PRECISION SUM, F, A(29)
      DATA    A/0.0D0,
     .          0.0D0,  0.333333333333333333333D0,
     .       -0.125D0,  0.016666666666666666667D0,
     .          0.0D0, -0.0001984126984126984127D0,
     .          0.0D0,  0.36743092298647854203D-5,
     .          0.0D0, -0.75156325156325156325D-7,
     .          0.0D0,  0.16059043836821614599D-8,
     .          0.0D0, -0.35227934257916621232D-10,
     .          0.0D0,  0.78720803121674581370D-12,
     .          0.0D0, -0.17840422612224120352D-13,
     .          0.0D0,  0.40886009791799259829D-15,
     .          0.0D0, -0.94559508632959211871D-17,
     .          0.0D0,  0.22036011313440918061D-18,
     .          0.0D0, -0.51683202540046382743D-20,
     .          0.0D0,  0.12188644964239543006D-20/

      SUM = A(4)*X**4
      F = X
      DO 100 N = 3, 29, 2
          F = X*X*F
          SUM = SUM + A(N)*F
100   CONTINUE
      INTLOW = SUM
      RETURN
      END


      DOUBLE PRECISION FUNCTION INTHIGH (X)
C       Integral of x**3/(exp(x)-1) from x to infinity.
C       Accurate for x greater than about 1.
      IMPLICIT NONE
      DOUBLE PRECISION X
      INTEGER N
      DOUBLE PRECISION  SUM

      SUM = 0.0D0
      DO 100 N = 1, 15
          SUM = SUM + EXP(-N*X)
     .           *(X**3/N + 3*X**2/N**2 + 6*X/N**3 + 6.0D0/N**4)
100   CONTINUE
      INTHIGH = SUM
      RETURN
      END





      SUBROUTINE GAUSQUADS (N, XA, WT)
C         Generates the abscissas (X) and weights (W) for an N point
C       Gauss-Legendre quadrature.  The XA are returned in this order: 
C       -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK, 
C       K=N/2, N must be even).
      IMPLICIT NONE
      INTEGER  N
      REAL     XA(*), WT(*)
      INTEGER  K, I, J, L
      DOUBLE PRECISION  X, XP, PL, PL1, PL2, DPL, TINY
      PARAMETER (TINY=3.0D-7)

      K = N/2
      IF (2*K .NE. N) STOP 'GAUSQUADS: N must be even'
      DO J = 1, K
        X = COS(3.141592654*(J-.25)/(N+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
          ENDDO
          DPL = N*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        XA(J)     = -X
        XA(J+K) = X
        WT(J  )   = 2.0D0/((1.0D0-X*X)*DPL*DPL)
        WT(J+K) = WT(J)
      ENDDO
      RETURN
      END




      SUBROUTINE DGAUSQUADS (N, XA, WT)
C         Generates the abscissas (X) and weights (W) for an N point
C       Double-Gauss-Legendre quadrature.  The XA are returned in this order: 
C       -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK, 
C       K=N/2, N must be even).
      IMPLICIT NONE
      INTEGER  N
      REAL     XA(*), WT(*)
      INTEGER  K, I, J, L, N2
      DOUBLE PRECISION  X, XP, PL, PL1, PL2, DPL, TINY
      PARAMETER (TINY=3.0D-7)

      N2 = N/2
      IF (2*N2 .NE. N) STOP 'DGAUSQUADS: N must be even'
      K = (N2+1)/2
      DO J = 1, K
        X = COS(3.141592654*(J-.25)/(N2+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO L = 2, N2
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
          ENDDO
          DPL = N2*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        XA(J+N2) = (1+X)/2
        XA(N2+1-J+N2) = (1-X)/2
        WT(J+N2) = 1/((1-X*X)*DPL*DPL)
        WT(N2+1-J+N2) = WT(J+N2)
        XA(J) = -XA(J+N2)
        XA(N2+1-J) = -XA(N2+1-J+N2)
        WT(J) = WT(J+N2)
        WT(N2+1-J) = WT(J+N2)
      ENDDO
      RETURN
      END




      SUBROUTINE SSORT (X, Y, N, KFLAG)
C***BEGIN PROLOGUE  SSORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   SSORT sorts array X and optionally makes the same interchanges in
C   array Y.  The array X may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      X - array of values to be sorted   (usually abscissas)
C      Y - array to be (optionally) carried along
C      N - number of values in array X to be sorted
C      KFLAG - control parameter
C            =  2  means sort X in increasing order and carry Y along.
C            =  1  means sort X in increasing order (ignoring Y)
C            = -1  means sort X in decreasing order (ignoring Y)
C            = -2  means sort X in decreasing order and carry Y along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***END PROLOGUE  SSORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
c      REAL X(*), Y(*)
      REAL X(*)
      INTEGER Y(*)
C     .. Local Scalars ..
c      REAL R, T, TT, TTY, TY
      REAL R, T, TT
      INTEGER TY, TTY
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(51), IU(51)
C     .. External Subroutines ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***First executable statement  SSORT
      NN = N
      IF (NN .LT. 1) THEN
         STOP 'The number of values to be sorted is not positive.'
      ENDIF

      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
        STOP 'The sort control parameter, K, is not 2, 1, -1, or -2.'
      ENDIF

C     Alter array X to get decreasing order if needed
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF

      IF (KK .EQ. 2) GO TO 100

C     Sort X only
      M = 1
      I = 1
      J = NN
      R = 0.375E0

   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF

   30 K = I

C     Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)

C     If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J

C     If last element of array is less than than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)

C        If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF

C     Find an element in the second half of the array which is smaller
C     than T
   40 L = L-1
      IF (X(L) .GT. T) GO TO 40

C     Find an element in the first half of the array which is greater
C     than T
   50 K = K+1
      IF (X(K) .LT. T) GO TO 50

C     Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF

C     Save upper and lower subscripts of the array yet to be sorted
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70

C     Begin again on another portion of the unsorted array
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1

   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = X(I+1)
      IF (X(I) .LE. T) GO TO 80
      K = I

   90 X(K+1) = X(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 90
      X(K+1) = T
      GO TO 80

C     Sort X and carry Y along
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0

  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
  120 K = I

C     Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)

C     If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)
      ENDIF
      L = J

C     If last element of array is less than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)

C        If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
         ENDIF
      ENDIF

C     Find an element in the second half of the array which is smaller
C     than T
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130

C     Find an element in the first half of the array which is greater
C     than T
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140

C     Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY
         GO TO 130
      ENDIF

C     Save upper and lower subscripts of the array yet to be sorted
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160

C     Begin again on another portion of the unsorted array
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1

  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)
      TY = Y(I+1)
      IF (X(I) .LE. T) GO TO 170
      K = I

  180 X(K+1) = X(K)
      Y(K+1) = Y(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y(K+1) = TY
      GO TO 170

C     Clean up
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
      RETURN
      END



