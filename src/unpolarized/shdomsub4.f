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
     .             NSCATANGLE, YLMSUN, PHASETAB, NSTPHASE)
Cf2py threadsafe
      IMPLICIT NONE
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
      INTEGER ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
Cf2py intent(in) :: ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
      INTEGER NMU, NPHI0MAX, NPHI0(*), NSTPHASE
Cf2py intent(in) :: NMU, NPHI0MAX, NPHI0, NSTPHASE
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
      REAL CAMX(*), CAMY(*), CAMZ(*)
      DOUBLE PRECISION CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX
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
      DOUBLE PRECISION X0, Y0, Z0
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
     .                 SFCGRIDPARMS, NPART, MURAY, PHIRAY, 
     .                 MU2, PHI2, X0, Y0, Z0, 
     .                 TOTAL_EXT, VISOUT(N))
C       WRITE(*,*) N, VISOUT(N)
  900   CONTINUE
      ENDDO
      DEALLOCATE(YLMDIR, SUNDIRLEG, SINGSCAT)
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
      
      SUBROUTINE GRADIENT_NORMCORR(NSTOKES, NX, NY, NZ, NPTS, NBPTS,
     .           NCELLS,NBCELLS, ML, MM, NCS, NLM, NSTLEG, NLEG,
     .           NUMPHASE, NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .           BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .           SFCTYPE, NSFCPAR, SFCGRIDPARMS, MAXNBC, NTOPPTS,
     .           NBOTPTS, BCPTR, BCRAD, GNDTEMP, GNDALBEDO, SKYRAD,
     .           WAVENO, WAVELEN, UNITS, XGRID, YGRID, ZGRID, GRIDPOS,
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, EXTINCT,
     .           ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, SHPTR, SOURCE,
     .           CAMX, CAMY, CAMZ, CAMMU, CAMPHI, NPIX, GRAD1, GRAD2,
     .           NORM1, NORM2, COST, MEASUREMENTS, RSHPTR, STOKESOUT,
     .           NPART, TOTAL_EXT, RADIANCE, NUMDER, PARTDER, DEXT,
     .           DALB, DIPHASE, DLEG, NSCATANGLE, YLMSUN, PHASETAB,
     .           NSTPHASE, DPHASETAB, DNUMPHASE, SOLARFLUX, NPX, NPY,
     .           NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS, EXTDIRP,
     .           UNIFORMZLEV, DPATH, DPTR, EXACT_SINGLE_SCATTER,
     .           WEIGHTS)
Cf2py threadsafe
      IMPLICIT NONE
      LOGICAL EXACT_SINGLE_SCATTER
Cf2py intent(in) ::  EXACT_SINGLE_SCATTER
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS
      INTEGER NBPTS, NCELLS, NBCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS, NBPTS, NBCELLS
      INTEGER NPX, NPY, NPZ
      REAL    DELX, DELY, XSTART, YSTART
      REAL    ZLEVELS(*)
      REAL    EXTDIRP(*)
      DOUBLE PRECISION UNIFORMZLEV
Cf2py intent(in) :: DELX, DELY, XSTART, YSTART NPX, NPY, NPZ, ZLEVELS, EXTDIRP, UNIFORMZLEV
      INTEGER ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
Cf2py intent(in) :: ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
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
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
Cf2py intent(in) :: SOLARFLUX, SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
Cf2py intent(in) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN
      REAL    MU(*), PHI(NMU,*), WTDO(NMU,*)
Cf2py intent(in) :: MU, PHI, WTDO
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      REAL    SFCGRIDPARMS(*), BCRAD(NSTOKES,*)
Cf2py intent(in) :: SFCGRIDPARMS, BCRAD
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    TOTAL_EXT(NPTS), LEGEN(NSTLEG,0:NLEG,*)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN, TOTAL_EXT
      REAL    DIRFLUX(*), FLUXES(2,*)
      REAL    SOURCE(NSTOKES, *), RADIANCE(NSTOKES,*)
Cf2py intent(in) :: DIRFLUX, FLUXES, SOURCE, RADIANCE
      REAL CAMX(*), CAMY(*), CAMZ(*)
      DOUBLE PRECISION CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX
      REAL   MEASUREMENTS(NSTOKES,*), DLEG(NSTLEG,0:NLEG,DNUMPHASE)
      REAL   DEXT(NBPTS,NUMDER), DALB(NBPTS,NUMDER)
      INTEGER DIPHASE(NBPTS,NUMDER)
Cf2py intent(in) :: MEASUREMENTS, DEXT ,DALB, DIPHASE, DLEG

      REAL  STOKESOUT(NSTOKES,NPIX), COST(NSTOKES)
      REAL  NORM1(NSTOKES), NORM2(NSTOKES)
      REAL  GRAD1(NSTOKES,NBPTS,NUMDER)
      REAL  GRAD2(NSTOKES,NBPTS,NUMDER)

Cf2py intent(out) :: GRAD1, GRAD2, NORM1, NORM2, COST, STOKESOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS
      INTEGER NUMDER, PARTDER(NUMDER)
Cf2py intent(in) :: NUMDER, PARTDER
      INTEGER NSCATANGLE, NSTPHASE
      REAL YLMSUN(NSTLEG,NLM), PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
      REAL DPHASETAB(NSTPHASE,DNUMPHASE,NSCATANGLE)
Cf2py intent(in) :: YLMSUN, PHASETAB
Cf2py intent(in) :: NSCATANGLE, YLMSUN, PHASETAB, DPHASETAB, NSTPHASE
      REAL DPATH(8*(NPX+NPY+NPZ),*), WEIGHTS(NSTOKES)
      INTEGER DPTR(8*(NPX+NPY+NPZ),*)
Cf2py intent(in) :: DPATH, DPTR, WEIGHTS

      REAL PIXEL_ERROR(NSTOKES)
      DOUBLE PRECISION RAYGRAD(NSTOKES,NBPTS,NUMDER), VISRAD(NSTOKES)
      INTEGER I, J, L, SIDE, IVIS
      LOGICAL VALIDRAD
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION U, R, PI
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT
      INTEGER M, ME, MS, NS

      INTEGER, ALLOCATABLE :: LOFJ(:)

      WRITE(*,*) 'NotImplementedError'
      STOP

      RETURN
      END

      
      SUBROUTINE GRADIENT_L2(NSTOKES, NX, NY, NZ, NPTS, NBPTS, NCELLS,
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
     .           NSTPHASE, DPHASETAB, DNUMPHASE, SOLARFLUX, NPX, NPY,
     .           NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS, EXTDIRP, 
     .           UNIFORMZLEV, DPATH, DPTR, EXACT_SINGLE_SCATTER,
     .           WEIGHTS)
Cf2py threadsafe
      IMPLICIT NONE
      LOGICAL EXACT_SINGLE_SCATTER
Cf2py intent(in) ::  EXACT_SINGLE_SCATTER
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS
      INTEGER NBPTS, NCELLS, NBCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NBPTS, NCELLS, NBCELLS
      INTEGER NPX, NPY, NPZ
      REAL    DELX, DELY, XSTART, YSTART
      REAL    ZLEVELS(*)
      REAL    EXTDIRP(*)
      DOUBLE PRECISION UNIFORMZLEV
Cf2py intent(in) :: DELX, DELY, XSTART, YSTART NPX, NPY, NPZ, ZLEVELS, EXTDIRP, UNIFORMZLEV
      INTEGER ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
Cf2py intent(in) :: ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
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
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
Cf2py intent(in) :: SOLARFLUX, SOLARMU, SOLARAZ
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
      REAL    CAMX(*), CAMY(*), CAMZ(*)
      DOUBLE PRECISION CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX 
      REAL   MEASUREMENTS(*), DLEG(0:NLEG,DNUMPHASE)
      REAL   DEXT(NBPTS,NUMDER), DALB(NBPTS,NUMDER)
      INTEGER DIPHASE(NBPTS,NUMDER)
Cf2py intent(in) :: MEASUREMENTS, DEXT ,DALB, DIPHASE, DLEG
      REAL VISOUT(NPIX), PIXEL_ERROR
      DOUBLE PRECISION  GRADOUT(NBPTS,NUMDER), COST
Cf2py intent(out) :: GRADOUT, COST, VISOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS
      INTEGER NUMDER, PARTDER(NUMDER)
Cf2py intent(in) :: NUMDER, PARTDER
      INTEGER NSCATANGLE, NSTPHASE
      REAL    YLMSUN(NLM), PHASETAB(NUMPHASE,NSCATANGLE)
      REAL    DPHASETAB(DNUMPHASE,NSCATANGLE)
Cf2py intent(in) :: NSCATANGLE, YLMSUN, PHASETAB, DPHASETAB, NSTPHASE
      REAL DPATH(8*(NPX+NPY+NPZ),*), WEIGHTS(NSTOKES)
      INTEGER DPTR(8*(NPX+NPY+NPZ),*)
Cf2py intent(in) :: DPATH, DPTR, WEIGHTS

      INTEGER I, J, L, K
      INTEGER N, M, ME, MS, ND
      DOUBLE PRECISION  RAYGRAD(NBPTS,NUMDER)
      REAL    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, COSSCAT
      DOUBLE PRECISION U, R, PI
      REAL, ALLOCATABLE :: YLMDIR(:)
      INTEGER, ALLOCATABLE :: LOFJ(:)
      REAL, ALLOCATABLE :: SINGSCAT(:), DSINGSCAT(:)
      DOUBLE PRECISION, ALLOCATABLE :: SUNDIRLEG(:)

      ALLOCATE (YLMDIR(NLM), LOFJ(NLM), SUNDIRLEG(0:NLEG))
      ALLOCATE (SINGSCAT(NUMPHASE), DSINGSCAT(DNUMPHASE))
      
      GRADOUT = 0.0D0
      
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
      COST = 0.0D0
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

        IF (SRCTYPE .NE. 'T') THEN
          COSSCAT = SOLARMU*MU2 + SQRT((1.0-SOLARMU**2)*(1.0-MU2**2))
     .                  *COS(SOLARAZ-PHI2)
          U = (NSCATANGLE-1)*(ACOS(COSSCAT)/PI) + 1
          J = MIN(NSCATANGLE-1,INT(U))
          U = U - J
          DO I = 1, NUMPHASE
            SINGSCAT(I) = (1-U)*PHASETAB(I,J) + U*PHASETAB(I,J+1)
          ENDDO
          DO I = 1, DNUMPHASE
            DSINGSCAT(I) = (1-U)*DPHASETAB(I,J)+U*DPHASETAB(I,J+1)
          ENDDO
          
          IF (NUMPHASE .GT. 0) THEN
            CALL LEGENDRE_ALL (COSSCAT, NLEG, SUNDIRLEG)
            DO L = 0, NLEG
              SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L+1)/(4*PI)
            ENDDO
          ENDIF
        ENDIF
        
        RAYGRAD = 0.0D0; VISOUT(N) = 0.0D0;
        CALL GRAD_INTEGRATE_1RAY (BCFLAG, IPFLAG, 
     .             NX, NY, NZ, NPTS, NCELLS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .             DELTAM, SRCTYPE, WAVELEN, SOLARMU,SOLARAZ,
     .             EXTINCT, ALBEDO, LEGEN,
     .             IPHASE, DIRFLUX, SHPTR,
     .             SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS,NPART,
     .             MURAY, PHIRAY, MU2, PHI2, X0, Y0, Z0, 
     .             VISOUT(N), RAYGRAD, RSHPTR, TOTAL_EXT, 
     .             RADIANCE, LOFJ, PARTDER, NUMDER, DSINGSCAT,
     .             DEXT, DALB, DIPHASE, DLEG, NBPTS, DNUMPHASE,
     .             SOLARFLUX, NPX, NPY, NPZ, DELX, DELY, XSTART,
     .             YSTART, ZLEVELS, EXTDIRP, UNIFORMZLEV,
     .             DPATH, DPTR, EXACT_SINGLE_SCATTER)
900     CONTINUE

        PIXEL_ERROR = VISOUT(N) - MEASUREMENTS(N)
        GRADOUT = GRADOUT + PIXEL_ERROR*RAYGRAD
        COST = COST + 0.5*PIXEL_ERROR**2
C        WRITE(*,*) N, VISOUT(N), MEASUREMENTS(N), COST
      ENDDO
      DEALLOCATE(YLMDIR, LOFJ, SUNDIRLEG, SINGSCAT, DSINGSCAT)
      RETURN
      END
      
      
      
      SUBROUTINE GRAD_INTEGRATE_1RAY (BCFLAG, IPFLAG, 
     .             NX, NY, NZ, NPTS, NCELLS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .             DELTAM, SRCTYPE, WAVELEN,SOLARMU,SOLARAZ,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .             DIRFLUX, SHPTR, SOURCE, YLMDIR, 
     .             YLMSUN, SUNDIRLEG, SINGSCAT, MAXNBC, 
     .             NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS,NPART,
     .             MURAY,PHIRAY, MU2, PHI2, X0, Y0, Z0, 
     .             RADOUT, RAYGRAD, RSHPTR, TOTAL_EXT,  
     .             RADIANCE, LOFJ, PARTDER, NUMDER, DSINGSCAT,
     .             DEXT, DALB, DIPHASE, DLEG, NBPTS, DNUMPHASE,
     .             SOLARFLUX, NPX, NPY, NPZ, DELX, DELY, XSTART,
     .             YSTART, ZLEVELS, EXTDIRP, UNIFORMZLEV,
     .             DPATH, DPTR, EXACT_SINGLE_SCATTER)

C       Integrates the source function through the extinction field 
C     (EXTINCT) backward in the direction (MURAY,PHIRAY) to find the 
C     outgoing radiance (RAD) at the point X0,Y0,Z0.
      IMPLICIT NONE
      LOGICAL EXACT_SINGLE_SCATTER
      INTEGER NPX, NPY, NPZ
      REAL    DELX, DELY, XSTART, YSTART, SOLARFLUX
      REAL    ZLEVELS(*)
      REAL    EXTDIRP(*)
      DOUBLE PRECISION UNIFORMZLEV
      INTEGER BCFLAG, IPFLAG, NX, NY, NZ, NPTS, NCELLS
      INTEGER ML, MM, NCS, NLM, NLEG, NUMPHASE, NBPTS
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, KK, II, GRIDPOINT
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR, SSP, IDR, IPA
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SHPTR(NPTS+1), RSHPTR(NPTS+1)
      INTEGER*2 CELLFLAGS(NCELLS)
      INTEGER IPHASE(NPTS,NPART), NPART, DNUMPHASE
      INTEGER BCPTR(MAXNBC,2), LOFJ(NLM), KI
      LOGICAL DELTAM
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    WAVELEN, SOLARMU, SOLARAZ
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    DIRFLUX(NPTS), TOTAL_EXT(NPTS), LEGEN(0:NLEG,NPTS) 
      DOUBLE PRECISION RAYGRAD(NBPTS,NUMDER), DIRGRAD(NBPTS,NUMDER)
      REAL    SOURCE(*), RADIANCE(*), DSINGSCAT(DNUMPHASE)
      REAL    YLMDIR(NLM), YLMSUN(NLM), SINGSCAT(NUMPHASE)
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
      REAL    MURAY, PHIRAY, MU2, PHI2, RADOUT, SRCSINGSCAT(8)
      REAL    SINGSCAT0(8), SINGSCAT1(8)
      DOUBLE PRECISION X0, Y0, Z0
      CHARACTER SRCTYPE*1, SFCTYPE*2
      REAL      DLEG(0:NLEG,DNUMPHASE), DEXT(NBPTS,NUMDER)
      REAL      DALB(NBPTS,NUMDER), SINGSCAT8(8), OSINGSCAT8(8)
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
      REAL DPATH(8*(NPX+NPY+NPZ),*), DEXTM
      INTEGER DPTR(8*(NPX+NPY+NPZ),*)
      
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
          OSINGSCAT8(I) = SINGSCAT8(I)
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
     .            DEXT, DALB, DIPHASE, DLEG, NBPTS, BCELL, DSINGSCAT,
     .            SINGSCAT8, OSINGSCAT8)

C         Interpolate the source and extinction to the current point
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XE, YE, ZE, F)
        SRCEXT1 = F(1)*SRCEXT8(1) + F(2)*SRCEXT8(2) +
     .            F(3)*SRCEXT8(3) + F(4)*SRCEXT8(4) +
     .            F(5)*SRCEXT8(5) + F(6)*SRCEXT8(6) +
     .            F(7)*SRCEXT8(7) + F(8)*SRCEXT8(8)
        SRCEXT1 = MAX(0.0,SRCEXT1)
        
        EXT1 = F(1)*EXTINCT8(1) + F(2)*EXTINCT8(2) +
     .         F(3)*EXTINCT8(3) + F(4)*EXTINCT8(4) +
     .         F(5)*EXTINCT8(5) + F(6)*EXTINCT8(6) +
     .         F(7)*EXTINCT8(7) + F(8)*EXTINCT8(8)

        
        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        IF (.NOT. OUTOFDOMAIN) THEN
          CALL GET_INTERP_KERNEL(BCELL,GRIDPTR,GRIDPOS,XE,YE,ZE,FB)
          DO KK=1,8
            GRAD1(KK,:) = FB(KK)*GRAD8(KK,:)
            SINGSCAT1(KK) = FB(KK)*SINGSCAT8(KK)
            SINGSCAT1(KK) = MAX(0.0,SINGSCAT1(KK))
          ENDDO
         ELSE
            GRAD1 = 0.0
         ENDIF
         
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
          CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XI, YI, ZI, F)
          SRCEXT0 = F(1)*SRCEXT8(1) + F(2)*SRCEXT8(2) +
     .              F(3)*SRCEXT8(3) + F(4)*SRCEXT8(4) +
     .              F(5)*SRCEXT8(5) + F(6)*SRCEXT8(6) +
     .              F(7)*SRCEXT8(7) + F(8)*SRCEXT8(8)
            EXT0 = F(1)*EXTINCT8(1) + F(2)*EXTINCT8(2) +
     .             F(3)*EXTINCT8(3) + F(4)*EXTINCT8(4) +
     .             F(5)*EXTINCT8(5) + F(6)*EXTINCT8(6) +
     .             F(7)*EXTINCT8(7) + F(8)*EXTINCT8(8)
          SRCEXT0 = MAX(0.0,SRCEXT0)

          IF (.NOT. OUTOFDOMAIN) THEN
            CALL GET_INTERP_KERNEL(BCELL,GRIDPTR,GRIDPOS,XI,YI,ZI,FB)
            DO KK=1,8
              GRAD0(KK,:) = FB(KK)*GRAD8(KK,:)
              SINGSCAT0(KK) = FB(KK)*SINGSCAT8(KK)
              SINGSCAT0(KK) = MAX(0.0, SINGSCAT0(KK))
            ENDDO
          ELSE
            GRAD0 = 0.0
          ENDIF

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

            IF (.NOT. OUTOFDOMAIN) THEN
              SRCGRAD = ( 0.5*(GRAD0+GRAD1)
     .          + 0.08333333333*(EXT0*GRAD1-EXT1*GRAD0)*DELS
     .           *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT

              IF (EXACT_SINGLE_SCATTER) THEN
                SRCSINGSCAT = ( 0.5*(SINGSCAT0+SINGSCAT1)
     .            + 0.08333333333*(EXT0*SINGSCAT1-EXT1*SINGSCAT0)*DELS
     .                    *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
              ENDIF
            ENDIF 

          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
            SRC = 0.0
            SRCGRAD = 0.0
            SRCSINGSCAT = 0.0
          ENDIF

          RAD = RAD + TRANSMIT*SRC*ABSCELL
          IF (.NOT. OUTOFDOMAIN) THEN
            DO KK = 1, 8
              GRIDPOINT = GRIDPTR(KK,BCELL)
              RAYGRAD(GRIDPOINT,:) = RAYGRAD(GRIDPOINT,:) +
     .                     TRANSMIT*SRCGRAD(KK,:)*ABSCELL
     
C             Add gradient component due to the direct solar beam propogation
              IF (EXACT_SINGLE_SCATTER) THEN
                II = 1
                GRIDPOINT = GRIDPTR(KK,BCELL)
                DO WHILE (DPTR(II,GRIDPOINT) .GT. 0)
                  SSP = DPTR(II,GRIDPOINT)
                  DO IDR = 1, NUMDER
                    IPA = PARTDER(IDR)

                    IF (DELTAM) THEN
                      IF (NUMPHASE .GT. 0) THEN
                        DEXTM = (1.0-ALBEDO(SSP,IPA)*
     .                     LEGEN(ML+1,IPHASE(SSP,IPA)))*DEXT(SSP,IDR)
                      ELSE
                        DEXTM = (1.0-ALBEDO(SSP,IPA)*LEGEN(ML+1,SSP))
     .                           *DEXT(SSP,IDR)
                      ENDIF
                    ELSE
                      DEXTM = DEXT(SSP,IDR)
                    ENDIF
                    RAYGRAD(SSP,IDR) = RAYGRAD(SSP,IDR) -
     .                  DPATH(II,GRIDPOINT)*DEXTM*ABSCELL*TRANSMIT*
     .                  SRCSINGSCAT(KK)
                  ENDDO
                  II = II + 1
                ENDDO
              ENDIF
            ENDDO
          ENDIF 
          
          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          SRCEXT1 = SRCEXT0
          GRAD1 = GRAD0
          SINGSCAT1 = SINGSCAT0
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
     .           GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
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
     .            NBPTS, BCELL, DSINGSCAT, SINGSCAT8, OSINGSCAT8)
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
      REAL    DSINGSCAT(*)
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
      REAL    FULL_SINGSCAT, SINGSCAT8(8), OSINGSCAT8(8)
      REAL    DLEG(0:NLEG,DNUMPHASE), DEXT(NBPTS,NUMDER)
      REAL    DALB(NBPTS,NUMDER), DEXTM, DALBM, DLEGM, W
      INTEGER DIPHASE(NBPTS,NUMDER), KK
      INTEGER IP, J, L, M, MS, ME, K, IS, NS, N, I, IB
      INTEGER IPA, LOFJ(NLM), PARTDER(NUMDER), RNS, RIS, IDR
      DOUBLE PRECISION DA, F, A, SECMU0
      
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
          SINGSCAT8(N) = OSINGSCAT8(I)
          GRAD8(N,:) = OGRAD8(I,:)
        ELSE IF (I .LT. 0) THEN
          EXTINCT8(N) = EXTINCT8(ABS(I))
          SRCEXT8(N) = SRCEXT8(ABS(I))
          SINGSCAT8(N) = SINGSCAT8(ABS(I))
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
     .              EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM)

                IF (NUMPHASE .EQ. 0) THEN
                  IF (L .LE. ML) THEN
                    A = LEGEN(L,K) + F/(1-F)
                  ELSE
                    A = LEGEN(L,K)/(1-F)
                  ENDIF
                  FULL_SINGSCAT = FULL_SINGSCAT + DA*SUNDIRLEG(L)*(
     .               DEXTM*ALBEDO(IP,IPA)*A + EXTINCT(IP,IPA)*DALBM*A
     .               + EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM)
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
          SINGSCAT8(N) = 0.0
          DO J = 1, NS
            SRCEXT8(N) = SRCEXT8(N) + SOURCE(IS+J)*YLMDIR(J)
          ENDDO

          DO IPA = 1, NPART
          
            IF (EXT.EQ.0.0) THEN
              W = 1.0
            ELSE
              W = EXTINCT(IP,IPA)/EXT
            ENDIF
            IF (W.EQ.0.0) CYCLE
          
            IF (NUMPHASE .GT. 0) THEN
              K = IPHASE(IP,IPA)
            ELSE
              K = IP
            ENDIF

            DA = ALBEDO(IP,IPA)*DIRFLUX(IP)*SECMU0*W
            IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
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
            ENDIF
            
            IF (NUMPHASE .GT. 0) THEN
              SINGSCAT8(N) = SINGSCAT8(N) + DA*SINGSCAT(K)
              IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
                SRCEXT8(N) = SRCEXT8(N) + DA*SINGSCAT(K)
              ENDIF
            ELSE
              F = LEGEN(ML+1,K)
              DO L = 0, NLEG
                IF (L .LE. ML) THEN
                  A = DA*(LEGEN(L,K) + F/(1-F))
                ELSE
                  A = DA*LEGEN(L,K)/(1-F)
                ENDIF
                SINGSCAT8(N) = SINGSCAT8(N) + A*SUNDIRLEG(L)
                IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
                  SRCEXT8(N) = SRCEXT8(N) + A*SUNDIRLEG(L)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          SINGSCAT8(N) = SINGSCAT8(N)*EXT
          SRCEXT8(N) = SRCEXT8(N)*EXT
          EXTINCT8(N) = EXT
        ENDIF
      ENDDO
 
      RETURN
      END

      
      SUBROUTINE PRECOMPUTE_PHASE_CHECK(NSCATANGLE, NUMPHASE, NSTPHASE,
     .                              NSTOKES, ML, NLM, NSTLEG, NLEG,
     .                              LEGEN, PHASETAB, DELTAM, NEGCHECK)
C       Precomputes the phase function as a function of scattering angle
C     for all the tabulated phase functions.
      IMPLICIT NONE
      INTEGER NSCATANGLE, NSTPHASE, NSTOKES, NSTLEG
      INTEGER ML, NLM, NLEG, NUMPHASE
Cf2py intent(in) :: NSCATANGLE, NSTPHASE, NSTOKES, NSTLEG, ML, NLM, NLEG, NUMPHASE
      REAL    LEGEN(0:NLEG,NUMPHASE)
Cf2py intent(in) :: LEGEN
      REAL    PHASETAB(NUMPHASE,NSCATANGLE)
Cf2py intent(out) :: PHASETAB
      LOGICAL DELTAM, NEGCHECK
Cf2py intent(in) :: DELTAM, NEGCHECK
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
          IF (NEGCHECK .AND. SUM .LT. 0.0) THEN
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
      