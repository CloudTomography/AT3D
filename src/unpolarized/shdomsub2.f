
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


      SUBROUTINE INIT_CELL_STRUCTURE (BCFLAG, IPFLAG, NX, NY,
     .                     NZ, NX1, NY1, NPTS, NCELLS, XGRID,  
     .                     YGRID, ZGRID, GRIDPOS, MAXIC, MAXIG, 
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
      INTEGER MAXIC, MAXIG
Cf2py intent(in) :: IPFLAG, BCFLAG, NX, NY, NZ, NX1, NY1, MAXIC, MAXIG
Cf2py intent(out) :: NPTS, NCELLS
      INTEGER GRIDPTR(8,MAXIC), NEIGHPTR(6,MAXIC), TREEPTR(2,MAXIC)
Cf2py intent(out) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER*2 CELLFLAGS(MAXIC)
Cf2py intent(out) :: CELLFLAGS
      REAL    XGRID(NX1), YGRID(NY1), ZGRID(NZ), GRIDPOS(3,MAXIG)
Cf2py intent(in) :: XGRID, YGRID, ZGRID
Cf2py intent(out) :: GRIDPOS
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


      SUBROUTINE INTERP_GRID (NPTS, NLEG, GRIDPOS,
     .               TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .               NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .               XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .               ALBEDOP, LEGENP, IPHASEP, NZCKD,
     .               ZCKD, GASABS, EXTMIN, SCATMIN,NPART, 
     .		     TOTAL_EXT, NBPTS)
C       Calls TRILIN_INTERP_PROP to interpolate the input arrays from 
C     the property grid to each internal grid point. 
      IMPLICIT NONE
      INTEGER NPTS, NLEG, NPART, NBPTS
      INTEGER IPHASE(NPTS,NPART)  
      REAL    GRIDPOS(3,NPTS), TOTAL_EXT(NPTS)
      REAL    TEMP(*), EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(0:NLEG,*)
      INTEGER IP, IPA
      
      INTEGER NPX, NPY, NPZ
      INTEGER NUMPHASE
      REAL DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)
      REAL TEMPP(*), EXTINCTP(NBPTS,NPART), ALBEDOP(NBPTS,NPART)
      REAL LEGENP(*)
      INTEGER IPHASEP(NBPTS,NPART)
      INTEGER NZCKD
      REAL ZCKD(*), GASABS(*)
      DOUBLE PRECISION EXTMIN, SCATMIN
      
C         Initialize: transfer the tabulated phase functions
      CALL TRILIN_INTERP_PROP (0.0, 0.0, 0.0, .TRUE., NLEG, 
     .                         TEMP, EXTINCT, ALBEDO, 
     .                         LEGEN(0,1), IPHASE, 
     .                         NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .                         XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .                         ALBEDOP, LEGENP, IPHASEP, NZCKD,
     .                         ZCKD, GASABS, EXTMIN, SCATMIN)

C         Trilinearly interpolate from the property grid to the adaptive grid
      TOTAL_EXT(:NPTS) = 0.0
      DO IPA = 1, NPART
	DO IP = 1, NPTS
	  CALL TRILIN_INTERP_PROP 
     .          (GRIDPOS(1,IP), GRIDPOS(2,IP), GRIDPOS(3,IP), 
     .           .FALSE., NLEG, TEMP(IP), EXTINCT(IP,IPA), 
     .            ALBEDO(IP,IPA), LEGEN(0,IP), IPHASE(IP,IPA), 
     .            NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .            XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP(:,IPA),
     .            ALBEDOP(:,IPA),LEGENP, IPHASEP(:,IPA),
     .            NZCKD, ZCKD, GASABS, EXTMIN, SCATMIN)
	  TOTAL_EXT(IP) = TOTAL_EXT(IP) + EXTINCT(IP,IPA)
	ENDDO
      ENDDO
       
      RETURN
      END
 




      SUBROUTINE MAKE_DIRECT (NPTS, BCFLAG, IPFLAG, DELTAM, 
     .                ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, 
     .                SOLARAZ, GRIDPOS, DIRFLUX, 
     .                NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .                XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .                ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .                ZCKD, GASABS, CX, CY, CZ, CXINV, CYINV,
     .                CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .                XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV,
     .		      NPART, NBPTS)
C       Makes the direct beam solar flux for the internal base grid.
C     DIRFLUX is set to F*exp(-tau_sun).
C     Actually calls DIRECT_BEAM_PROP to do all the hard work.
      IMPLICIT NONE
      INTEGER NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, NBPTS
Cf2py intent(in) :: NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, NBPTS
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS(3,*)
Cf2py intent(in) :: SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS
      INTEGER NPX, NPY, NPZ, NPART, NUMPHASE
Cf2py intent(in) :: NPX, NPY, NPZ, NPART, NUMPHASE
      REAL DELX, DELY, XSTART, YSTART
Cf2py intent(in) :: DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*), TEMPP(*)
      REAL EXTINCTP(NBPTS,NPART), ALBEDOP(NBPTS,NPART)
Cf2py intent(in) :: ZLEVELS, TEMPP, EXTINCTP, ALBEDOP
      REAL LEGENP(*), ZCKD(*), GASABS(*)
      INTEGER IPHASEP(NBPTS,NPART), NZCKD
Cf2py intent(in) :: LEGENP, ZCKD, GASABS, IPHASEP, NZCKD
    
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
Cf2py intent(out) :: CX, CY, CZ, CXINV, CYINV, CZINV
      INTEGER IPDIRECT, DI, DJ, DK
Cf2py intent(out) :: IPDIRECT, DI, DJ, DK
      DOUBLE PRECISION EPSS, EPSZ, XDOMAIN, YDOMAIN
Cf2py intent(out) :: EPSS, EPSZ, XDOMAIN, YDOMAIN
      DOUBLE PRECISION UNIFORMZLEV, DELXD,DELYD
Cf2py intent(out) :: UNIFORMZLEV, DELXD, DELYD
      REAL    DIRFLUX(*), EXTDIRP(*)
Cf2py intent(in, out) :: DIRFLUX, EXTDIRP

      INTEGER SIDE, IP
      LOGICAL VALIDBEAM
      REAL    UNIFZLEV, XO, YO, ZO, DIR, DIRPATH


      CALL DIRECT_BEAM_PROP (1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,
     .         DELTAM,ML,NLEG, SOLARFLUX,SOLARMU,SOLARAZ, DIRFLUX(1),
     .         UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM, 
     .         NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .         XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .         ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .         ZCKD, GASABS, CX, CY, CZ, CXINV, CYINV,
     .         CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .         XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV, 
     .	       NPART, NBPTS)

      DO IP = 1, NPTS
        DIRPATH = 0.0
        CALL DIRECT_BEAM_PROP 
     .           (0, GRIDPOS(1,IP), GRIDPOS(2,IP), GRIDPOS(3,IP),
     .            BCFLAG, IPFLAG, DELTAM,ML,NLEG,
     .            SOLARFLUX,SOLARMU,SOLARAZ, DIRFLUX(IP),
     .            UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM, 
     .            NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .            XSTART, YSTART, ZLEVELS, TEMPP, EXTINCTP,
     .            ALBEDOP, LEGENP, EXTDIRP, IPHASEP, NZCKD,
     .            ZCKD, GASABS, CX, CY, CZ, CXINV, CYINV,
     .            CZINV, DI, DJ, DK, IPDIRECT, DELXD, DELYD,
     .            XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV, 
     .		  NPART, NBPTS)
      ENDDO
      RETURN
      END
 

 
 
  
      SUBROUTINE PREPARE_PROP (ML, MM, NLEG, NPTS, DELTAM, NUMPHASE, 
     .               SRCTYPE, UNITS, WAVENO, WAVELEN, ALBMAX, 
     .               EXTINCT, ALBEDO, LEGEN, TEMP, PLANCK, IPHASE,
     .		     NPART, TOTAL_EXT)
C       Prepares the grid arrays for the iterative solution process.
C       If doing Delta-M scaling then the extinction, albedo, and Legendre
C       terms are scaled first; only the 0 to ML LEGEN terms are scaled.
C       Outputs PLANCK with (1-omega)*B(T) for thermal source, where B(T) is
C       the Planck function (meaning depends on UNITS).  
C       TEMP array is unchanged.
      IMPLICIT NONE
      INTEGER ML, MM, NLEG, NPTS, NUMPHASE
      INTEGER IPHASE(NPTS,NPART), NPART
      LOGICAL DELTAM
      REAL  WAVENO(2), WAVELEN, ALBMAX, TOTAL_EXT(NPTS)
      REAL  EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL  LEGEN(0:NLEG,*), TEMP(NPTS), PLANCK(NPTS,NPART)
      CHARACTER*1 SRCTYPE, UNITS
      INTEGER I, K, L, IPA
      REAL    F, BB

      ALBMAX = 0.0 

      IF (DELTAM .AND. NUMPHASE .GT. 0) THEN
        DO K = 1, NUMPHASE
          F = LEGEN(ML+1,K)
          DO L = 0, ML
            LEGEN(L,K) = (LEGEN(L,K) - F)/(1-F)
          ENDDO
        ENDDO
      ENDIF

      IF (DELTAM) TOTAL_EXT(:NPTS) = 0.0
      
      DO IPA = 1, NPART
	DO I = 1, NPTS
	  IF (DELTAM) THEN
	    IF (NUMPHASE .GT. 0) THEN
	      F = LEGEN(ML+1,IPHASE(I,IPA))
	    ELSE
	      F = LEGEN(ML+1,I)
	      DO L = 0, ML
		LEGEN(L,I) = (LEGEN(L,I) - F)/(1-F)
	      ENDDO
	    ENDIF
	    EXTINCT(I,IPA) = (1.0-ALBEDO(I,IPA)*F)*EXTINCT(I,IPA)
	    ALBEDO(I,IPA) = (1.0-F)*ALBEDO(I,IPA)/
     .			    (1.0-ALBEDO(I,IPA)*F)
            TOTAL_EXT(I) = TOTAL_EXT(I) + EXTINCT(I,IPA)
	  ENDIF
	  ALBMAX = MAX(ALBMAX, ALBEDO(I,IPA))
	  IF (SRCTYPE .NE. 'S') THEN
	    CALL PLANCK_FUNCTION (TEMP(I), UNITS,WAVENO,WAVELEN,BB)
	    PLANCK(I,IPA) = (1.0-ALBEDO(I,IPA))*BB
	  ENDIF
	ENDDO
      ENDDO
      RETURN
      END





      SUBROUTINE INIT_RADIANCE (NXY, NZ, NCS, NLEG, RSHPTR, ZGRID,
     .             EXTINCT, ALBEDO, LEGEN, TEMP, NUMPHASE, IPHASE,
     .             SRCTYPE, SOLARFLUX, SOLARMU, GNDALBEDO, GNDTEMP, 
     .             SKYRAD, UNITS, WAVENO, WAVELEN, RADIANCE, NPART,
     .		   TOTAL_EXT)
C       Initializes radiance field by solving plane-parallel two-stream.
C     Solves the L=1 M=0 SH system by transforming the pentadiagonal
C     system to tridiagonal and calling solver.
C     Does the initialization for the NXY columns of the base grid.
      IMPLICIT NONE
      INTEGER NXY, NZ, NCS, NLEG, NUMPHASE, RSHPTR(*)
      INTEGER IPHASE(NZ,NXY,NPART), NPART
      REAL    GNDALBEDO, GNDTEMP, SKYRAD, SOLARFLUX, SOLARMU
      REAL    WAVENO(2), WAVELEN
      REAL    ZGRID(NZ), TOTAL_EXT(NZ,NXY)
      REAL    EXTINCT(NZ,NXY,NPART), ALBEDO(NZ,NXY,NPART)
      REAL    LEGEN(0:NLEG,*)
      REAL    TEMP(NZ,NXY)
      REAL    RADIANCE(*)
      CHARACTER SRCTYPE*1, UNITS*1
      INTEGER NLAYER, I, IZ, IR, J, L

      LOGICAL DELTAM
      REAL    PI, C0, C1
      REAL    EXT0, EXT1, SCAT0, SCAT1, G0, G1, GNDEMIS
      REAL, ALLOCATABLE :: OPTDEPTHS(:), ALBEDOS(:)
      REAL, ALLOCATABLE :: TEMPS(:), FLUXES(:,:), ASYMMETRIES(:)
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
          EXT0 = TOTAL_EXT(IZ,I)
          EXT1 = TOTAL_EXT(IZ+1,I)
          SCAT0 = SUM(ALBEDO(IZ,I,:)*EXTINCT(IZ,I,:))
          SCAT1 = SUM(ALBEDO(IZ+1,I,:)*EXTINCT(IZ+1,I,:)) 
          OPTDEPTHS(L) = (ZGRID(IZ+1)-ZGRID(IZ))* (EXT0+EXT1)/2
          IF (EXT0+EXT1 .GT. 0.0) THEN
            ALBEDOS(L) = (SCAT0+SCAT1)/(EXT0+EXT1)
          ELSE
            ALBEDOS(L) = 0.0
          ENDIF
          IF (NUMPHASE .GT. 0) THEN
            G0 = SUM(ALBEDO(IZ,I,:)*EXTINCT(IZ,I,:)*
     .			LEGEN(1,IPHASE(IZ,I,:)))
	    G1 = SUM(ALBEDO(IZ+1,I,:)*EXTINCT(IZ+1,I,:)*
     .			LEGEN(1,IPHASE(IZ+1,I,:)))
          ELSE
            J = NZ*(I-1)+IZ
            G0 = LEGEN(1,J)
            G1 = LEGEN(1,J+1)
          ENDIF
          IF (SCAT0+SCAT1 .GT. 0.0) THEN
            ASYMMETRIES(L) = (G0 + G1)/(SCAT0+SCAT1)
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
          RADIANCE(IR+1) = C0*(FLUXES(1,L)+FLUXES(2,L))
          RADIANCE(IR+2) = 0.0
          RADIANCE(IR+1+NCS) = C1*(FLUXES(1,L)-FLUXES(2,L))
          RADIANCE(IR+2+NCS) = 0.0
          IR = IR + 2+NCS
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
      REAL      WAVENO, WAVELEN
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
 


 
 

      SUBROUTINE MAKE_ANGLE_SET (NMU, NPHI, NCS, ITYPE, NPHI0MAX,
     .             NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)
C       Make the set of angles for the discrete space representation.
C     The number of mu's (cosine zenith angles) and maximum number of
C     phi's (azimuth angles) is input.  The actual number of azimuthal 
C     angles for each mu is output (NPHI0).  There are three types of
C     discrete ordinate sets: ITYPE=1 is a gaussian grid, 2 is a reduced
C     gaussian grid, and 3 is a reduced double gaussian set.
C     If NCS=1 then only about half the azimuthal angles (from 0 to pi) 
C     are used because the radiance is even in phi (cosine terms).  
C     The output is the NMU mu values, the NPHI0 phi values for each mu,
C     and the integration weight for each ordinate. The first NMU/2 mu 
C     angles are the downwelling (mu<0) angles.  Also output are the
C     maximum number of azimuthal angles (NPHI0MAX), the flags for doing 
C     an azimuthal FFT for each mu and the total number of angles.
      IMPLICIT NONE
      INTEGER NMU, NPHI, NCS, ITYPE,  NPHI0MAX, NPHI0(NMU), NANG
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
          IF (NCS .EQ. 1) THEN
            NPHI0MAX = INT((NPHI+2)/2)
            IF (ITYPE .EQ. 1) THEN
              NPHI0(J) = NPHI0MAX
            ELSE
C               For the reduced gaussian set, make the smaller NPHI0 values
C               still be good for the FFT (but don't let NPHI0MAX be exceeded).
              IF (NPHI0MAX .GT. MAXNPHI) 
     .          STOP 'MAKE_ANGLE_SET: exceeded GOODNFFT'
              NPHI0(J) = INT(0.9+1+(NPHI0MAX-1)*SQRT(1-MU(J)**2))
              IF (NPHI0(J) .GT. 1)  
     .          NPHI0(J) = MIN(NPHI0MAX,GOODNFFT(NPHI0(J)-1)+1)
            ENDIF
C               Compute the azimuth angles and weights
            DELPHI = ACOS(-1.0)/MAX(1,NPHI0(J)-1)
            DO K = 1, NPHI0(J)
              PHI(J,K) = (K-1)*DELPHI
              IF ((K.EQ.1 .OR. K.EQ.NPHI0(J)) .AND. NPHI0(J).NE.1) THEN
                WTDO(J,K) = DELPHI*WTMU(J)
              ELSE
                WTDO(J,K) = 2.0*DELPHI*WTMU(J)
              ENDIF
            ENDDO
          ELSE
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
          ENDIF
          NANG = NANG + NPHI0(J)
          FFTFLAG(J) = (NPHI0(J) .GT. 14) .OR. (MM .GT. 15)
        ENDDO
        
      ELSE
        STOP 'MAKE_ANGLE_SET: invalid discrete ordinate type'
      ENDIF
      RETURN
      END



      SUBROUTINE MAKE_SH_DO_COEF (ML, MM, NLM, NMU, NPHI0, NCS, 
     .             NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG,
     .             CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)
C       Makes the transformation coefficients for the spherical harmonic
C     transform.  The main coefficients are output in four arrays: 1 is for
C     the SH_TO_DO forward transform, 2 is for the DO_TO_SH back transform,
C     which contains the discrete angle integration weights.
C     The coefficients are divided up into the mu dependent set CMUn 
C     (function of l, m, mu_j), and the phi dependent set CPHIn 
C     (function of m, phi_k) for each mu_j.
C     The FFTPACK phase coefficients for the FFT in azimuth are also 
C     output in WPHISAVE.
      IMPLICIT NONE
      INTEGER ML, MM, NLM, NMU, NPHI0(NMU), NPHI0MAX, NCS
      LOGICAL FFTFLAG(NMU)
      REAL    MU(NMU), PHI(NMU,*), WTMU(NMU), WTDO(NMU,*)
      REAL    CMU1(NLM,NMU), CMU2(NMU,NLM)
      REAL    CPHI1(-16:16,32,NMU), CPHI2(32,-16:16,NMU)
      REAL    WPHISAVE(3*NPHI0MAX+15,NMU)
      INTEGER I, J, K, M
      REAL    W

C         Make the to and from associate Legendre coefficients       
      DO I = 1, NMU
        CALL YLMALL (MU(I), 0.0, ML, MM, -NCS, CMU1(1,I))
      ENDDO
      DO I = 1, NMU
        DO J = 1, NLM
          CMU2(I,J) = CMU1(J,I)*WTMU(I)
        ENDDO
      ENDDO


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
          IF (NCS .EQ. 1) THEN
            CALL COSTI (NPHI0(I),WPHISAVE(1,I))
          ELSE
            CALL RFFTI (NPHI0(I),WPHISAVE(1,I))
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END
 
 
 


      SUBROUTINE SURFACE_BRDF (SFCTYPE, REFPARMS, WAVELEN,
     .                         MU2, PHI2, MU1, PHI1, REFLECT)
C       Returns the reflection coefficient for the general bidirectional
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
C       F  Fresnel       Real, Imaginary part of index of refraction
C       R  Rahman et al  rho0, k, Theta
C       O  Ocean         Wind Speed (m/s), Pigmentation (mg/m^3)
      IMPLICIT NONE
      REAL   REFPARMS(*), WAVELEN, MU1, PHI1, MU2, PHI2, REFLECT
      CHARACTER  SFCTYPE*1
      REAL   PI
      REAL   FRESNEL_REFLECTION, RPV_REFLECTION

      PI = ACOS(-1.0)
      IF (SFCTYPE .EQ. 'L' .OR. SFCTYPE .EQ. 'l') THEN
C         L or l: Lambertian surface BRDF is constant.
C           (for testing, as done more efficiently by Lambertian routines)
        REFLECT = REFPARMS(1)
      ELSE IF (SFCTYPE .EQ. 'F') THEN
C         F: Fresnel reflection for a dielectric interface
        REFLECT = FRESNEL_REFLECTION (REFPARMS(1), REFPARMS(2), 
     .                                MU1, MU2, PHI1-PHI2)
      ELSE IF (SFCTYPE .EQ. 'R') THEN
C         R: Rahman, Pinty, and Verstraete
        REFLECT = RPV_REFLECTION (REFPARMS(1),REFPARMS(2),REFPARMS(3),
     .                            -MU1, MU2, PHI1-PHI2-PI)
      ELSE IF (SFCTYPE .EQ. 'O') THEN
C         O: Ocean BRDF from 6S modified by Norm Loeb
        CALL ocean_brdf_sw (REFPARMS(1), -1., REFPARMS(2), WAVELEN,
     .                      -MU1, MU2, PHI1-PHI2, PHI1, REFLECT)
      ELSE
        STOP 'SURFACE_BRDF: Unknown BRDF type'
      ENDIF

      RETURN
      END



      LOGICAL FUNCTION SPECULAR_SURFACE (SFCTYPE)
C       Returns true for specularly reflecting surfaces, and false otherwise.
C     Needed because specular surfaces are a special case, because the
C     discrete ordinate integration cannot resolve the delta function
C     reflectance function.  For now, only the Fresnel surface is specular.
      CHARACTER SFCTYPE*1
      SPECULAR_SURFACE = SFCTYPE .EQ. 'F'
      RETURN
      END



      REAL FUNCTION FRESNEL_REFLECTION (MRE, MIM, MU1, MU2, PHI)
C       Computes the unpolarized Fresnel reflection from a dielectric
C     interface.  The complex index of refraction of the dielectric
C     is specified with MRE,MIM.  The reflection is zero if the MU and
C     PHI of the incident and outgoing directions are not within 0.01.
C     The reflection coefficient is the average of the horizontal and
C     vertical polarizations.
      REAL    MRE, MIM, MU1, MU2, PHI
      REAL    MU, REF
      COMPLEX EPSILON, D, RH, RV

      MU = ABS(MU1)
      IF (ABS(MU-MU2) .LT. 0.01 .AND. ABS(PHI) .LT. 0.01) THEN
        EPSILON = CMPLX(MRE,MIM)**2
        D = SQRT(EPSILON - 1.0 + MU**2)
        RH = (MU - D) / (MU + D)
        RV =(EPSILON*MU - D) / (EPSILON*MU + D)
        REF = (CABS(RH)**2 + CABS(RV)**2 )/2.0
      ELSE 
        REF = 0.0
      ENDIF
      FRESNEL_REFLECTION = REF
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






      SUBROUTINE CHECK_OCEAN_BRDF_INTEG (NXSFC,NYSFC,NSFCPAR,SFCPARMS,
     .                                   WAVELEN, SOLARMU, SOLARAZ,
     .                                   NMU, NPHI0, NPHI0MAX,
     .                                   MU, PHI, WTDO)
C       Checks that there is enough resolution in the discrete ordinates
C      to resolve the specular peak in the ocean BRDF.  The width of the
C      peak depends mainly on the windspeed, so the minimum windspeed in
C      the variable surface is used for the test.  The ocean BRDF is 
C      integrated over the outgoing directions for the incident solar
C      direction to find the albedo.  The ocean albedo calculated with
C      the existing discrete ordinates must agree within MAXALBEDODIF with 
C      that calculated from a high resolution ordinate set (NMU2 mu's).
      IMPLICIT NONE
      INTEGER NXSFC, NYSFC, NSFCPAR
      INTEGER NMU, NPHI0MAX, NPHI0(NMU)
      REAL    SFCPARMS(NSFCPAR,NXSFC+1,NYSFC+1)
      REAL    WAVELEN, SOLARMU, SOLARAZ
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      INTEGER IX, IY, JMU, JPHI, NMU2
      PARAMETER (NMU2=64)
      REAL    WINDSPD, PIGMENT, PI, OPI, REFLECT, W, ALBEDO1, ALBEDO2
      REAL    MU2(NMU2), WT2(NMU2), DELPHI, PHI2, MAXALBEDODIF
      PARAMETER (MAXALBEDODIF=0.002)

      IF (WAVELEN .LT. 0.25 .OR. WAVELEN .GT. 4.00) THEN
        WRITE (6,'(A,A)') 'CHECK_OCEAN_BRDF_INTEG: Ocean BRDF requires',
     .                    ' wavelength between 0.25 and 4.0 microns.'
        STOP
      ENDIF

C       Find the minimum wind speed and pigment content in ocean surface
      WINDSPD = 1000.0
      PIGMENT = 1.0E6
      DO IY = 1, NYSFC
        DO IX = 1, NXSFC
          WINDSPD = MIN(WINDSPD,SFCPARMS(2,IX,IY))
          PIGMENT = MIN(PIGMENT,SFCPARMS(3,IX,IY))
        ENDDO
      ENDDO
      IF (WINDSPD .LT. 0.5) THEN
        WRITE (6,'(A,A)') 'CHECK_OCEAN_BRDF_INTEG: Ocean BRDF requires',
     .                    ' wind speed no less than 0.5 m/s.'
        STOP
      ENDIF

C       Integrate the ocean BRDF for incident Sun over outgoing directions
C         to get the albedo using the original discrete ordinate set.
      PI = ACOS(-1.0)
      OPI = 1.0/PI
      ALBEDO1 = 0.0
      DO JMU = NMU/2+1, NMU
        DO JPHI = 1, NPHI0(JMU)
          CALL ocean_brdf_sw (WINDSPD, -1., PIGMENT, WAVELEN, 
     .                        -SOLARMU, MU(JMU), 
     .                        SOLARAZ-PHI(JMU,JPHI), SOLARAZ, REFLECT)
          W = OPI*ABS(MU(JMU))*WTDO(JMU,JPHI)
          ALBEDO1 = ALBEDO1  + W*REFLECT
        ENDDO
      ENDDO

C       Integrate the ocean BRDF to get the albedo for the high resolution
C         discrete ordinate set.
      CALL DGAUSQUADS (NMU2, MU2, WT2)
      DELPHI = 2*PI/(2*NMU2)
      ALBEDO2 = 0.0
      DO JMU = NMU2/2+1, NMU2
        DO JPHI = 1, 2*NMU2
          PHI2 = DELPHI*(JPHI-1)
          CALL ocean_brdf_sw (WINDSPD, -1., PIGMENT, WAVELEN, 
     .                        -SOLARMU, MU2(JMU), 
     .                        SOLARAZ-PHI2, SOLARAZ, REFLECT)
          W = OPI*ABS(MU2(JMU))*WT2(JMU)*DELPHI
          ALBEDO2 = ALBEDO2  + W*REFLECT
        ENDDO
      ENDDO

      IF (ABS(ALBEDO1-ALBEDO2) .GT. MAXALBEDODIF) THEN
        WRITE (6,*)
        WRITE (6,'(A,A)') 'CHECK_OCEAN_BRDF_INTEG: Warning - ',
     .                ' insufficient discrete ordinate resolution to'
        WRITE (6,'(A,A,2F6.3)') 'resolve peak in ocean BRDF.',
     .        '  Reference and approx albedos:', ALBEDO2, ALBEDO1
        WRITE (6,*)
      ENDIF
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




      SUBROUTINE COMPUTE_NETFLUXDIV (NPTS, RSHPTR, SRCTYPE, 
     .             SOLARMU, EXTINCT, ALBEDO, PLANCK, DIRFLUX,
     .             RADIANCE,  NETFLUXDIV)
C       Computes the net flux divergence at every grid point.
      IMPLICIT NONE
      INTEGER NPTS, RSHPTR(NPTS+1)
      REAL    PLANCK(NPTS), DIRFLUX(NPTS), SOLARMU
      REAL    EXTINCT(NPTS), ALBEDO(NPTS)
      REAL    RADIANCE(*), NETFLUXDIV(NPTS)
      CHARACTER SRCTYPE*1
      INTEGER I, IR
      REAL    PI, C, SECMU, HR

      PI = ACOS(-1.0)
      C = SQRT(4.0*PI)
      SECMU = 1.0/ABS(SOLARMU)
      DO I = 1, NPTS
        IR = RSHPTR(I)
        HR = C*EXTINCT(I)*(1-ALBEDO(I))*RADIANCE(IR+1)
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




      SUBROUTINE COMPUTE_SH (NSHOUT, NPTS, SRCTYPE, 
     .                       SOLARMU, SOLARAZ, DIRFLUX,
     .                       RSHPTR, ML,MM,NCS, RADIANCE, SHTERMS)
C       Computes the quantities for the spherical harmonic output.
C     At each grid point the mean radiance and net flux (Fx, Fy, Fz)
C     are computed.  The output includes the diffuse and direct solar
C     components.  If NSHOUT=5 then the HIGHORDERRAD flag is set and the 
C     rms of the higher order terms in the radiance series is computed.
      IMPLICIT NONE
      INTEGER NPTS, RSHPTR(NPTS+1), ML, MM, NCS, NSHOUT
      REAL    SOLARMU, SOLARAZ, DIRFLUX(NPTS)
      REAL    RADIANCE(*), SHTERMS(NSHOUT,NPTS)
      CHARACTER SRCTYPE*1
      INTEGER I, IR, J, NLM
      REAL    PI, C, C0, IMEAN, FX, FY, FZ, HORMS
      REAL    SECSOL, SINSOL, SOLX, SOLY, SOLM, F0   
      REAL, ALLOCATABLE :: YLMSUN(:)
 
      NLM = (2*MM+1)*(ML+1) - MM*(MM+1)
      ALLOCATE (YLMSUN(NLM))

C             Spherical Harmonic output: individual SH terms
      PI = ACOS(-1.0)
      C = SQRT(4.0*PI/3.0)
      C0 = 1.0/SQRT(4.0*PI)
      SECSOL = 1.0/ABS(SOLARMU)
      SINSOL = SQRT(1.0-SOLARMU**2)
      SOLX = SINSOL*COS(SOLARAZ)*SECSOL
      SOLY = SINSOL*SIN(SOLARAZ)*SECSOL
      SOLM = SECSOL/(4.0*PI)
      F0 = 0.0
      IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
        CALL YLMALL (SOLARMU, SOLARAZ, ML, MM, NCS, YLMSUN)
      ENDIF
      DO I = 1, NPTS
        IR = RSHPTR(I)
        IMEAN = C0*RADIANCE(IR+1)
        FZ = C*RADIANCE(IR+1+NCS)
        IF (MM .GT. 0) THEN
          FX = -C*RADIANCE(IR+2+NCS)
        ELSE
          FX = 0.0
        ENDIF
        IF (MM .GT. 0 .AND. NCS .NE. 1) THEN
          FY = -C*RADIANCE(IR+2)
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
          DO J = NCS+2+1, RSHPTR(I+1)-IR
            HORMS = HORMS + (F0*YLMSUN(J)+RADIANCE(IR+J))**2
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




      SUBROUTINE INTERP_RADIANCE (OLDNPTS, NPTS, RSHPTR, RADIANCE,
     .               NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
C       Interpolates the spherical harmonic radiance array to the new
C     grid points (>OLDNPTS) from the old grid points.
      IMPLICIT NONE
      INTEGER OLDNPTS, NPTS, NBCELLS, NCELLS
      INTEGER RSHPTR(NPTS+1), GRIDPTR(8,NCELLS), TREEPTR(2,NCELLS)
      REAL    RADIANCE(*), GRIDPOS(3,*)
      INTEGER IP, IC, ICELL, IPARENT, IDIR, I, J, DIR1, DIR2, IP1, IP2
      INTEGER IR, IR1, IR2, NR, NR1, NR2
      INTEGER GRIDCORNER(2,4,3)
      REAL    RAD1, RAD2
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
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END



      SUBROUTINE INTEGRATE_1RAY (BCFLAG, IPFLAG, 
     .                        NX, NY, NZ, NPTS, NCELLS, 
     .                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                        XGRID, YGRID, ZGRID, GRIDPOS,
     .                        ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .                        NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                        DELTAM, SRCTYPE, WAVELEN,SOLARMU,SOLARAZ,
     .                        EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .                        DIRFLUX, SHPTR, SOURCE, 
     .                        YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .                        MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                        SFCTYPE, NSFCPAR, SFCGRIDPARMS, NPART, 
     .                        MURAY,PHIRAY, MU2,PHI2, X0,Y0,Z0,
     .			              TOTAL_EXT, RADOUT)

C       Integrates the source function through the extinction field 
C     (EXTINCT) backward in the direction (MURAY,PHIRAY) to find the 
C     outgoing radiance (RAD) at the point X0,Y0,Z0.

      IMPLICIT NONE
      INTEGER BCFLAG, IPFLAG, NX, NY, NZ, NPTS, NCELLS
      INTEGER ML, MM, NCS, NLM, NLEG
      INTEGER NUMPHASE, NPART
      INTEGER MAXNBC, NTOPPTS, NBOTPTS
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SHPTR(NPTS+1)
      INTEGER*2 CELLFLAGS(NCELLS)
      INTEGER IPHASE(NPTS,NPART)
      INTEGER BCPTR(MAXNBC,2)
      LOGICAL DELTAM
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    WAVELEN, SOLARMU, SOLARAZ
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    DIRFLUX(NPTS), LEGEN(0:NLEG,NPTS)
      REAL    SOURCE(*), TOTAL_EXT(NPTS)
      REAL    YLMDIR(NLM), YLMSUN(NLM), SINGSCAT(NUMPHASE)
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
      REAL    MURAY, PHIRAY, MU2, PHI2, RADOUT
      DOUBLE PRECISION X0, Y0, Z0
      CHARACTER SRCTYPE*1, SFCTYPE*2
      

      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), DONETHIS(8)
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    OEXTINCT8(8), OSRCEXT8(8), EXTINCT8(8), SRCEXT8(8)
      REAL    EXT0, EXT1, EXTN, SRCEXT0, SRCEXT1, RADBND
      REAL    XM,YM
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
      DOUBLE PRECISION TAUTOL, TAUGRID, S, DELS, TRANSCUT
      DOUBLE PRECISION EXT, SRC, TAU, TRANSCELL,ABSCELL, TRANSMIT, RAD
      DOUBLE PRECISION U,V,W, DELX,DELY,DELZ, INVDELX,INVDELY,INVDELZ
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
        ENDDO

C         Compute the source function times extinction in direction (MU2,PHI2)
        CALL COMPUTE_SOURCE_1CELL (ICELL, GRIDPTR, 
     .        ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .        NPTS, DELTAM, SRCTYPE, SOLARMU, EXTINCT,
     .        ALBEDO, LEGEN, IPHASE,
     .        DIRFLUX, SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, 
     .        SINGSCAT,DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8, 
     .        EXTINCT8, SRCEXT8, TOTAL_EXT, NPART)

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
        SRCEXT1 = (1-W)*((1-V)*((1-U)*SRCEXT8(1) + U*SRCEXT8(2))
     .                     + V*((1-U)*SRCEXT8(3) + U*SRCEXT8(4))) 
     .              + W*((1-V)*((1-U)*SRCEXT8(5) + U*SRCEXT8(6))
     .                     + V*((1-U)*SRCEXT8(7) + U*SRCEXT8(8)))

        SRCEXT1 = MAX(0.0,SRCEXT1)
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
     .      MURAY,PHIRAY,XE,YE,ZE,SO,ICELL,IOPP,SOX,SOY,SOZ, 
     .      GRIDPOS(1,IOPP),GRIDPOS(2,IOPP),GRIDPOS(3,IOPP)
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
          SRCEXT0 = (1-W)*((1-V)*((1-U)*SRCEXT8(1) + U*SRCEXT8(2))
     .                    + V*((1-U)*SRCEXT8(3) + U*SRCEXT8(4))) 
     .             + W*((1-V)*((1-U)*SRCEXT8(5) + U*SRCEXT8(6))
     .                    + V*((1-U)*SRCEXT8(7) + U*SRCEXT8(8)))
          SRCEXT0 = MAX(0.0,SRCEXT0)
	  
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
          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
            SRC = 0.0
          ENDIF

          RAD = RAD + TRANSMIT*SRC*ABSCELL
          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          SRCEXT1 = SRCEXT0
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




      SUBROUTINE FIND_BOUNDARY_RADIANCE (XB, YB, MU2, PHI2, 
     .                      ICELL, KFACE, GRIDPTR, GRIDPOS,
     .                      MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                      NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                      SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, DIRFLUX, 
     .                      SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                      RADBND)
C       Returns the interpolated radiance at the boundary (RADBND).
C     Inputs are the boundary location (XB,YB), ray direction away from
C     boundary (MU2,PHI2), cell number (ICELL) and face (KFACE) at
C     the boundary point.
      IMPLICIT NONE
      INTEGER ICELL, KFACE, MAXNBC, NTOPPTS, NBOTPTS
      INTEGER GRIDPTR(8,*), BCPTR(MAXNBC,2)
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
      REAL    MU2, PHI2, RADBND
      DOUBLE PRECISION XB, YB
      REAL    GRIDPOS(3,*)
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(*)
      CHARACTER SRCTYPE*1, SFCTYPE*2

      INTEGER IL, IM, IU, IP, IBC, J
      LOGICAL LAMBERTIAN
      REAL    X(4), Y(4), RAD(4), U, V
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
          RAD(J) = BCRAD(IBC)

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
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, BCRAD(1+NTOPPTS))
          ENDIF
          RAD(J) = BCRAD(NTOPPTS+IBC)
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
      RADBND = (1-U)*(1-V)*RAD(1) + U*(1-V)*RAD(2)
     .         + (1-U)*V*RAD(3) + U*V*RAD(4)

      RETURN
      END



 
      SUBROUTINE COMPUTE_SOURCE_1CELL (ICELL, GRIDPTR, 
     .            ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .            NPTS, DELTAM, SRCTYPE, SOLARMU, EXTINCT,
     .            ALBEDO, LEGEN, IPHASE, DIRFLUX, SHPTR, 
     .            SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .            DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8, 
     .            EXTINCT8, SRCEXT8, TOTAL_EXT, NPART)
C       Computes the source function times extinction for gridpoints 
C     belonging to cell ICELL in the direction (MU,PHI).  The results 
C     are returned in SRCEXT8 and EXTINCT8.
C     The spherical harmonic source function series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
      IMPLICIT NONE
      INTEGER ICELL, NPTS, ML, MM, NCS, NLM, NLEG
      INTEGER NUMPHASE, NPART
Cf2py intent(in) :: ICELL, NPTS, ML, MM, NCS, NLM, NLEG, NUMPHASE
      INTEGER GRIDPTR(8,*), SHPTR(*)
Cf2py intent(in) :: GRIDPTR, SHPTR      
      INTEGER DONETHIS(8), OLDIPTS(8)
Cf2py intent(in) :: DONETHIS, OLDIPTS      
      INTEGER IPHASE(NPTS,NPART)
Cf2py intent(in) :: IPHASE    
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    SOLARMU
Cf2py intent(in) :: SOLARMU
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(0:NLEG,*),  TOTAL_EXT(NPTS)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN
      REAL    DIRFLUX(*), SOURCE(*)
Cf2py intent(in) :: DIRFLUX, SOURCE
      REAL    YLMDIR(*), YLMSUN(*), SINGSCAT(*)
Cf2py intent(in) :: YLMDIR, YLMSUN, SINGSCAT
      REAL    OEXTINCT8(8), OSRCEXT8(8)
Cf2py intent(in) :: OEXTINCT8, OSRCEXT8
      REAL    EXTINCT8(8), SRCEXT8(8), EXT
Cf2py intent(out) :: EXTINCT8, SRCEXT8
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
Cf2py intent(in) :: SUNDIRLEG
      CHARACTER SRCTYPE*1
Cf2py intent(in) :: SRCTYPE
      INTEGER IP, J, L, M, MS, ME, K, IS, NS, N, I
      INTEGER IPA
      DOUBLE PRECISION DA, F, A, SECMU0, W
 
      SECMU0 = 1.0D0/ABS(SOLARMU)

C         Loop over the grid points, computing the source function 
C           at the viewing angle from the spherical harmonic source function.
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
	
	      EXT = TOTAL_EXT(IP)
          OLDIPTS(N) = IP
          IS = SHPTR(IP)
          NS = SHPTR(IP+1)-IS
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
 


 
      SUBROUTINE PRECOMPUTE_PHASE (NSCATANGLE, NUMPHASE, NSTPHASE, 
     .                      NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)
C       Precomputes the phase function as a function of scattering angle
C     for all the tabulated phase functions.
      IMPLICIT NONE
      INTEGER NSCATANGLE, ML, NLEG
      INTEGER NUMPHASE, NSTOKES, NSTLEG, NSTPHASE
Cf2py intent(in) :: NSCATANGLE, NUMPHASE, ML, NLEG,  NSTOKES, NSTLEG, NSTPHASE
      REAL    LEGEN(0:NLEG,NUMPHASE)
Cf2py intent(in) :: LEGEN
      REAL    PHASETAB(NUMPHASE,NSCATANGLE)
Cf2py intent(out) :: PHASETAB

      INTEGER I, J, L, MAXLEG
      DOUBLE PRECISION PI, SUM, F, A, COSSCAT, LEGSCAT(0:NLEG)

 
      PI = ACOS(-1.0D0)
      DO J = 1, NSCATANGLE
        COSSCAT = COS(PI*DFLOAT(J-1)/(NSCATANGLE-1))
C           Compute the Legendre polynomials for the scattering angle 
C             for the untruncated solar single scattering computation.
        CALL LEGENDRE_ALL (COSSCAT, NLEG, LEGSCAT)
        DO L = 0, NLEG
          LEGSCAT(L) = LEGSCAT(L)*(2*L+1)/(4*PI)
        ENDDO
        DO I = 1, NUMPHASE
          SUM = 0.0D0
          F = LEGEN(ML+1,I)
          DO L = 0, NLEG
            IF (L .LE. ML) THEN
              A = (LEGEN(L,I) + F/(1-F))
            ELSE
              A = LEGEN(L,I)/(1-F)
            ENDIF
            SUM = SUM + A*LEGSCAT(L)
          ENDDO
          IF (SUM .LE. 0.0) THEN
            WRITE (6,*) 'PRECOMPUTE_PHASE: negative source function',
     .          ' for tabulated phase function: ',I, J, SUM
            STOP
          ENDIF
          PHASETAB(I,J) = SNGL(SUM)
        ENDDO
      ENDDO

      RETURN
      END
 

      SUBROUTINE COMPUTE_RADIANCE (NX, NY, NZ, NPTS, NCELLS,
     .             ML, MM, NCS, NLEG, NUMPHASE, 
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, 
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,
     .             XGRID, YGRID, ZGRID, GRIDPOS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,
     .             SHPTR, SOURCE, SOURCE1, GRIDRAD,
     .             OUTPARMS, NRAD, NANGOUT, RADOUT)
C       Computes the radiances for the specified locations and directions.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS, NANGOUT
Cf2py intent(in) :: NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS, NANGOUT
      INTEGER ML, MM, NCS, NLEG
      INTEGER NUMPHASE
Cf2py intent(in) :: ML, MM, NCS, NLEG, NUMPHASE
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NRAD
Cf2py intent(in) :: NMU, NPHI0MAX, NPHI0
Cf2py intent(in, out) :: NRAD
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
Cf2py intent(in) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER SHPTR(NPTS+1), BCPTR(MAXNBC,2)
Cf2py intent(in) :: SHPTR, BCPTR      
      INTEGER*2 CELLFLAGS(NCELLS)
      INTEGER IPHASE(NPTS)
Cf2py intent(in) :: CELLFLAGS, IPHASE
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    SOLARMU, SOLARAZ
Cf2py intent(in) :: SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
Cf2py intent(in) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN
      REAL    MU(NMU), PHI(NMU,NPHI0MAX), WTDO(NMU,NPHI0MAX)
Cf2py intent(in) :: MU, PHI, WTDO
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      REAL    SFCGRIDPARMS(*), BCRAD(*)
Cf2py intent(in) :: SFCGRIDPARMS, BCRAD
      REAL    EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(0:NLEG,*)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN
      REAL    DIRFLUX(NPTS), FLUXES(2,NPTS)
Cf2py intent(in) :: DIRFLUX, FLUXES
      REAL    SOURCE1(NPTS), GRIDRAD(NPTS), SOURCE(*)
Cf2py intent(in) :: GRIDRAD, SOURCE
Cf2py intent(in) :: SOURCE1
      REAL    OUTPARMS(*)
Cf2py intent(in) :: OUTPARMS
      REAL    RADOUT(NANGOUT,*)
Cf2py intent(in, out) :: RADOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS

      INTEGER I, IBC, JX, JY, K, SIDE
      INTEGER NXOUT, NYOUT
      LOGICAL LAMBERTIAN, VALIDRAD
      DOUBLE PRECISION X0,Y0,Z0, XE,YE,ZE, TRANSMIT, RADIANCE
      REAL    MUOUT, PHIOUT, PHID
      REAL    XDOMAIN, YDOMAIN, STARTX, STARTY


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


      LAMBERTIAN = SFCTYPE(2:2) .EQ. 'L'

C         Loop over the radiance directions
      DO K = 1, NANGOUT
        MUOUT = OUTPARMS(2*K+5)
        PHID = OUTPARMS(2*K+6)
        PHIOUT = PHID*ACOS(-1.0)/180.0
        IF (MUOUT .EQ. 0.0 .OR. ABS(MUOUT) .GT. 1.0) THEN
          WRITE (6,*) 'COMPUTE_RADIANCE: Bad mu for radiance',MUOUT
        ELSE          

C             Compute the source function throughout grid for this angle
          CALL COMPUTE_ONE_SOURCE (ML, MM, NCS, NLEG, NUMPHASE,
     .           NPTS, DELTAM, MUOUT, PHIOUT, 
     .           SRCTYPE, SOLARMU, SOLARAZ, ALBEDO, LEGEN, IPHASE, 
     .           DIRFLUX, SHPTR, SOURCE, SOURCE1)

C             Set the radiance field to -1, so we can determine valid radiances
C               Also set the source array to the extinction times the source
C               function, so the grid interpolation is correct.
          DO I = 1, NPTS
            GRIDRAD(I) = -1.0
            SOURCE1(I) = SOURCE1(I)*EXTINCT(I)
          ENDDO
C             Get boundary radiances: either top or bottom
C             Isotropic top boundary or Lambertian bottom boundary can use
C               the previously computed boundary radiances in BCRAD,
C               otherwise, compute the radiance for this angle
C               by integrating over the stored downwelling radiances.
          IF (MUOUT .LT. 0.0) THEN
            DO IBC = 1, NTOPPTS
              I = BCPTR(IBC,1)
              GRIDRAD(I) = BCRAD(IBC)
            ENDDO
          ELSE
            IF (.NOT. LAMBERTIAN) THEN
              CALL VARIABLE_BRDF_SURFACE (NBOTPTS,1,NBOTPTS,BCPTR(1,2),
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MUOUT, PHIOUT,
     .             SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX, 
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, BCRAD(1+NTOPPTS))
            ENDIF
            DO IBC = 1, NBOTPTS
              I = BCPTR(IBC,2)
              GRIDRAD(I) = BCRAD(NTOPPTS+IBC)
            ENDDO
          ENDIF
	  NRAD=0
C             Integrate backward from the location to get the radiance
          Y0 = STARTY
          DO JY = 1, NYOUT
            X0 = STARTX
            DO JX = 1, NXOUT
              NRAD = NRAD + 1
              TRANSMIT = 1.0D0 ; RADIANCE = 0.0D0
              CALL INTEGRATE_SOURCE (BCFLAG, IPFLAG, 
     .                        NX, NY, NZ, NPTS, NCELLS, 
     .                        GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                        XGRID, YGRID, ZGRID, GRIDPOS,
     .                        MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1,
     .                        X0, Y0, Z0, XE,YE,ZE, SIDE,
     .                        TRANSMIT, RADIANCE, VALIDRAD)
              RADOUT(K,NRAD) = RADIANCE
              X0 = X0 + OUTPARMS(2)
            ENDDO
            Y0 = Y0 + OUTPARMS(3)
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END

 
      SUBROUTINE COMPUTE_ONE_SOURCE (ML, MM, NLM, NCS, NLEG, NUMPHASE,
     .             NPTS, DELTAM, MU, PHI, SRCTYPE, SOLARMU, SOLARAZ,
     .             ALBEDO, LEGEN, IPHASE, DIRFLUX, SHPTR,
     .             SOURCE, SOURCE1)
C       Computes the source function (SOURCE1) in the direction (MU,PHI)
C     for the whole domain (NX,NZ).  The spherical harmonic source function
C     series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
      IMPLICIT NONE
      INTEGER NPTS, ML, MM, NLM, NCS, NLEG, SHPTR(*)
      INTEGER NUMPHASE
Cf2py intent(in) :: NPTS,ML,MM,NLM,NCS,NLEG,NUMPHASE,SHPTR   
      INTEGER IPHASE(*)
Cf2py intent(in) :: IPHASE
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    SOLARMU, SOLARAZ, MU, PHI
Cf2py intent(in) :: SOLARMU, SOLARAZ, MU, PHI
      REAL    ALBEDO(*), LEGEN(0:NLEG,*), DIRFLUX(*)
Cf2py intent(in) :: ALBEDO, LEGEN, DIRFLUX
      REAL    SOURCE(*)
Cf2py intent(in) :: SOURCE
      REAL    SOURCE1(*)
Cf2py intent(in, out) :: SOURCE1    
      CHARACTER SRCTYPE*1
Cf2py intent(in) :: SRCTYPE

      INTEGER I, J, L, M, MS, ME, K, IS, NS
      REAL, ALLOCATABLE :: YLMDIR(:), YLMSUN(:)
      DOUBLE PRECISION, ALLOCATABLE :: SUNDIRLEG(:)
      REAL, ALLOCATABLE :: SINGSCAT(:)
      DOUBLE PRECISION DA, F, A, SECMU0, COSSCAT
 
      ALLOCATE (YLMDIR(NLM), YLMSUN(NLM))
      ALLOCATE (SUNDIRLEG(0:NLEG), SINGSCAT(NUMPHASE))

      IF (ABS(MU).GT.1.0) STOP 'COMPUTE_ONE_SOURCE: Bad mu'
 

C         Precompute Ylm's for output direction and solar direction
      CALL YLMALL (MU, PHI, ML, MM, NCS, YLMDIR)
      IF (SRCTYPE .NE. 'T') THEN
        CALL YLMALL (SOLARMU, SOLARAZ, ML, MM, NCS, YLMSUN)
        SECMU0 = 1.0/ABS(SOLARMU)
        IF (DELTAM) THEN
C           Compute the Legendre polynomials for the scattering angle 
C           for the untruncated solar single scattering computation.
          COSSCAT = SOLARMU*MU + SQRT((1.0D0-SOLARMU**2)*(1.0D0-MU**2))
     .              *COS(SOLARAZ-PHI)
          CALL LEGENDRE_ALL (COSSCAT, NLEG, SUNDIRLEG)
          A = 1.0D0/(4.0D0*ACOS(-1.0D0))
          DO L = 0, NLEG
            SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L+1)*A
          ENDDO
        ENDIF
      ENDIF
C           If tabulated phase functions then precompute the 
C           single scattering contribution except for
C           the direct beam transmission weighting times albedo.
      IF (SRCTYPE .NE. 'T' .AND. DELTAM .AND. NUMPHASE .GT. 0) THEN
        DO I = 1, NUMPHASE
          SINGSCAT(I) = 0.0
          F = LEGEN(ML+1,I)
          DO L = 0, NLEG
            IF (L .LE. ML) THEN
              A = (LEGEN(L,I) + F/(1-F))
            ELSE
              A = LEGEN(L,I)/(1-F)
            ENDIF
            SINGSCAT(I) = SINGSCAT(I) + A*SUNDIRLEG(L)
          ENDDO
          IF (SINGSCAT(I) .LE. 0.0) THEN
            WRITE (6,*) 'COMPUTE_ONE_SOURCE: negative source function',
     .        'for tabulated phase function: ',I
            STOP
          ENDIF
        ENDDO
      ENDIF


C         Loop over all the grid points, computing the source function 
C           at the viewing angle from the spherical harmonic source function.
      DO I = 1, NPTS
        IS = SHPTR(I)
        NS = SHPTR(I+1)-IS
C           Sum over the spherical harmonic series of the source function
        SOURCE1(I) = 0.0
        DO J = 1, NS
           SOURCE1(I) = SOURCE1(I) + SOURCE(IS+J)*YLMDIR(J)
        ENDDO

C           Special case for solar source and Delta-M
        IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
          IF (NUMPHASE .GT. 0) THEN
            K = IPHASE(I)
          ELSE
            K = I
          ENDIF
C             First subtract off the truncated single scattering 
          DA = ALBEDO(I)*DIRFLUX(I)*SECMU0
          J = 1
          DO L = 0, ML
            ME = MIN(L,MM)
            MS = (1-NCS)*ME
            A = DA*LEGEN(L,K)
            IF (J .LE. NS) THEN
              DO M = MS, ME
                SOURCE1(I) = SOURCE1(I) - A*YLMDIR(J)*YLMSUN(J)
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
            SOURCE1(I) = SOURCE1(I) + DA*SINGSCAT(K)
          ELSE
            F = LEGEN(ML+1,K)
            DO L = 0, NLEG
              IF (L .LE. ML) THEN
                A = DA*(LEGEN(L,K) + F/(1-F))
              ELSE
                A = DA*LEGEN(L,K)/(1-F)
              ENDIF
              SOURCE1(I) = SOURCE1(I) + A*SUNDIRLEG(L)
            ENDDO
          ENDIF
        ENDIF
      ENDDO
 
      DEALLOCATE (YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT)
      RETURN
      END
 



      SUBROUTINE INTEGRATE_SOURCE (BCFLAG, IPFLAG, 
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
C     The tranmsission and radiance of the ray so far (TRANSMIT, RADIANCE)
C     are input and returned after the integration along with the exitting
C     ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
C     5=-Z,6=+Z).
      IMPLICIT NONE
      INTEGER BCFLAG, IPFLAG, NX, NY, NZ, NPTS, NCELLS, SIDE
Cf2py intent(in) :: BCFLAG, IPFLAG, NX, NY, NZ, NPTS, NCELLS, SIDE      
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR  
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in) :: CELLFLAGS     
      LOGICAL VALIDRAD
Cf2py intent(in) :: VALIDRAD     
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      REAL    MU, PHI
Cf2py intent(in) :: MU, PHI
      DOUBLE PRECISION X0, Y0, Z0, XE,YE,ZE, TRANSMIT
Cf2py intent(in) :: X0, Y0, Z0, XE,YE,ZE, TRANSMIT
      DOUBLE PRECISION RADIANCE
Cf2py intent(in) :: RADIANCE
      REAL    EXTINCT(*), SOURCE(*), GRIDRAD(*)
Cf2py intent(in) :: EXTINCT, SOURCE, GRIDRAD

      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER I1, I2, I3, I4, IOPP, NTAU, IT
      LOGICAL IPINX, IPINY, OPENBCFACE, BTEST
      INTEGER GRIDFACE(4,6), OPPFACE(6), JFACE, KFACE, IC
      REAL    EXT0, EXT1, EXTN, SRCEXT0, SRCEXT1, RAD0, XM,YM
      DOUBLE PRECISION PI, CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS
      DOUBLE PRECISION TAUTOL, TAUGRID, S, DELS
      DOUBLE PRECISION EXT, SRC, TAU, TRANSCELL,ABSCELL, RAD
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
      CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                        XE, YE, ZE, 1, 1, EXTINCT, EXT1)
      CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                        XE, YE, ZE, 1, 1, SOURCE, SRCEXT1)
      SRCEXT1 = MAX(0.0,SRCEXT1)

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
     .                          XN, YN, ZN, 1, 1, EXTINCT, EXTN)
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
     .                              XI, YI, ZI, 1, 1, EXTINCT, EXT0)
          ENDIF
          CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                            XI, YI, ZI, 1, 1, SOURCE, SRCEXT0)
          SRCEXT0 = MAX(0.0,SRCEXT0)

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
          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
            SRC = 0.0
          ENDIF

          RADIANCE = RADIANCE + TRANSMIT*SRC*ABSCELL
          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          SRCEXT1 = SRCEXT0
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
        ELSE IF (GRIDRAD(I1) .GE. 0.0 .AND. GRIDRAD(I2) .GE. 0.0 .AND.
     .      GRIDRAD(I3) .GE. 0.0 .AND. GRIDRAD(I4) .GE. 0.0) THEN
          VALIDRAD = .TRUE.
          CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                            XI, YI, ZI, 1, 1, GRIDRAD, RAD0)
          RADIANCE = RADIANCE + TRANSMIT*RAD0
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
     .                              X, Y, Z, ND, ID, FIELD, VALUE)
C       Trilinearly interpolates the field (FIELD) given a grid cell
C     pointer (ICELL), to compute the field value (VALUE) at the desired
C     point (X,Y,Z). ND is the leading index of the FIELD (1 for a 
C     scalar field), and the ID is the element to interpolate (if vector).
      IMPLICIT NONE
      INTEGER ICELL, GRIDPTR(8,*), ND, ID
Cf2py intent(in) :: ICELL, GRIDPTR, ND, ID
      DOUBLE PRECISION X, Y, Z
Cf2py intent(in) :: X, Y, Z     
      REAL    GRIDPOS(3,*), FIELD(ND,*)
Cf2py intent(in) :: GRIDPOS, FIELD   
      REAL VALUE 
Cf2py intent(out) :: VALUE
      INTEGER IPT1, IPT2
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
      VALUE = WM*(VM*(UM*FIELD(ID,GRIDPTR(1,ICELL))
     .               + U*FIELD(ID,GRIDPTR(2,ICELL)))
     .           + V*(UM*FIELD(ID,GRIDPTR(3,ICELL))
     .               + U*FIELD(ID,GRIDPTR(4,ICELL)))) 
     .       + W*(VM*(UM*FIELD(ID,GRIDPTR(5,ICELL))
     .               + U*FIELD(ID,GRIDPTR(6,ICELL)))
     .           + V*(UM*FIELD(ID,GRIDPTR(7,ICELL))
     .               + U*FIELD(ID,GRIDPTR(8,ICELL))))
      RETURN
      END


 

      SUBROUTINE LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID, 
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, X0, Y0, Z0, ICELL)
C       Locates the grid cell in the tree structure containing the
C     specified point (X0,Y0,Z0), and returns the cell pointer ICELL.
C     First the base grid cell is found (using NX,NY,NZ, XGRID,YGRID,ZGRID),
C     and the tree is traced down to find the smallest cell containing 
C     the point.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NCELLS, BCFLAG, IPFLAG
Cf2py intent(in) :: NX, NY, NZ, NCELLS, BCFLAG, IPFLAG
      INTEGER ICELL
Cf2py intent(out) :: ICELL
      INTEGER GRIDPTR(8,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, TREEPTR      
      INTEGER*2 CELLFLAGS(*)
Cf2py intent(in) :: CELLFLAGS      
      DOUBLE PRECISION X0, Y0, Z0
Cf2py intent(in) :: X0, Y0, Z0     
      REAL    XGRID(*), YGRID(*), ZGRID(*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID  
      REAL    GRIDPOS(3,*)
Cf2py intent(in) :: GRIDPOS      
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
Cf2py intent(in) :: NLEG       
      DOUBLE PRECISION COSSCAT 
Cf2py intent(in) :: COSSCAT
      DOUBLE PRECISION P(0:NLEG)
Cf2py intent(in,out) :: P
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




      SUBROUTINE YLMALL (MU, PHI, ML, MM, NCS, P)
C       This subroutine computes a set of normalized spherical harmonic 
C     functions, P(J), for a particular direction mu,phi. 
C     ML is the maximum meridional mode, MM is the maximum azimuthal mode,
C     and NCS is the azimuthal mode flag (|NCS|=1 for cosine only, |NCS|=2 for 
C     sines and cosines).  Returns normalized associated Legendre functions 
C     only if NCS<0. The set is returned for triangular truncation: 
C     J = NCS*(L*(L+1))/2 + M+1  for L<=MM
C     J = (NCS*MM+1)*L-MM*(2+NCS*(MM-1))/2 + M+1  for L>MM
      IMPLICIT NONE
      INTEGER ML, MM, NCS
Cf2py intent(in) :: ML, MM, NCS
      REAL    MU, PHI
Cf2py intent(in) :: MU, PHI
      REAL    P(*)
Cf2py intent(in,out) :: P
      INTEGER J, L, M, C
      DOUBLE PRECISION X, Y, A, PMM, PL, PL1, PL2, PHI8

      C = ABS(NCS)
      IF (C .NE. 1 .AND. C .NE. 2)  STOP 'YLMALL: bad NCS'
      IF (MM .GT. ML)  STOP 'YLMALL: MM greater than LM'
      IF (MU*MU .GT. 1.) THEN
	    WRITE(*,*) MU
	    STOP 'YLMALL: |MU| larger than 1'
      ENDIF
      X = DBLE(MU)
      Y = SQRT(1.0D0-X*X)
C         Use the stable upward recursion on l, starting from P^m_m
C         Put in the spherical harmonic normalization as it goes
      PMM = 1.0D0/SQRT(4.0D0*ACOS(-1.0D0))
      DO M = 0, MM
        IF (M .GT. 0)  PMM = -PMM*Y*SQRT((2*M+1.0D0)/(2.0D0*M))
        J = C*(M*(M+1))/2 + M+1
        P(J) = PMM
        PL2 = PMM
        IF (M .LT. ML) THEN
          IF (M+1.LE.MM) J=C*((M+1)*(M+2))/2 +M+1
          IF (M+1.GT.MM) J=(C*MM+1)*(M+1)-(MM*(2+C*(MM-1)))/2+M+1
          PL1 = SQRT(2*M+3.0D0)*X*PMM
          P(J) = PL1
        ENDIF
        DO L = M+1, ML-1
          IF (L+1.LE.MM) J=C*((L+1)*(L+2))/2 +M+1
          IF (L+1.GT.MM) J=(C*MM+1)*(L+1)-(MM*(2+C*(MM-1)))/2+M+1
          A = 1.0D0/((L+M+1.D0)*(L-M+1.D0))
          PL = SQRT((2*L+1)*A*(2*L+3)) *X*PL1
     .             - SQRT((2*L+3)*A*(L+M)*(L-M)/(2*L-1.)) *PL2
          P(J) = PL
          PL2 = PL1
          PL1 = PL
        ENDDO
        IF (M .EQ. 0) PMM = PMM*SQRT(2.0D0)
      ENDDO
C         If there are M<0 terms then fill them in
      IF (C .EQ. 2) THEN
        DO L = 0, ML
          DO M = 1, MIN(L,MM)
            IF (L .LE. MM) J = L*(L+1) +M+1
            IF (L .GT. MM) J = MM*(2*L-MM) +L+M+1
            P(J-2*M) = P(J)
          ENDDO
        ENDDO
      ENDIF
C         Put in the azimuthal dependence
      PHI8 = DBLE(PHI)
      IF (NCS .GT. 0) THEN
        DO M = (1-NCS)*MM, MM
          IF (M .LT. 0) THEN
            A = SIN(-M*PHI8)
          ELSE IF (M .GT. 0) THEN
            A = COS(M*PHI8)
          ELSE
            A = 1.0D0
          ENDIF
          DO L = ABS(M), ML
            IF (L .LE. MM) THEN
              J = C*(L*(L+1))/2 +M+1
            ELSE
              J = (C*MM+1)*L-(MM*(2+C*(MM-1)))/2 + M+1
            ENDIF
            P(J) = A*P(J)
          ENDDO
        ENDDO
      ENDIF
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






      SUBROUTINE DSORT (X, Y, N, KFLAG)
C***BEGIN PROLOGUE  DSORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   DSORT sorts array X and optionally makes the same interchanges in
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
      REAL*8 X(*)
      INTEGER Y(*)
C     .. Local Scalars ..
c      REAL R, T, TT, TTY, TY
      REAL*8 R, T, TT
      INTEGER TY, TTY
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(51), IU(51)
C     .. External Subroutines ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***First executable statement  DSORT
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