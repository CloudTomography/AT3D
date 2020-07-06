C     This file containts subroutines that were modified from their original purpose
C     The original subroutines were written by Frank Evans for the Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer.
C     The modified subroutines were written by Aviad Levis, Technion Institute of Technology, 2019

      SUBROUTINE RENDER (NSTOKES, NX, NY, NZ,
     .                   NPTS, NCELLS, ML, MM, NCS, NLM, NSTLEG, NLEG,
     .                   NUMPHASE, NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .                   BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU,
     .                   SOLARAZ, SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                   MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .                   GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN,
     .                   UNITS, XGRID, YGRID, ZGRID, GRIDPOS,
     .                   GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                   EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .                   FLUXES, SHPTR, SOURCE, CAMX, CAMY, CAMZ, CAMMU,
     .                   CAMPHI, NPIX, NPART, TOTAL_EXT, STOKES,
     .                   NSCATANGLE, YLMSUN, PHASETAB, NSTPHASE)
Cf2py threadsafe
      IMPLICIT NONE
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
      INTEGER ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
Cf2py intent(in) :: ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
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
      REAL    SFCGRIDPARMS(*), BCRAD(NSTOKES, *)
Cf2py intent(in) :: SFCGRIDPARMS, BCRAD
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    TOTAL_EXT(NPTS), LEGEN(NSTLEG,0:NLEG,*)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN, TOTAL_EXT

      REAL    DIRFLUX(*), FLUXES(2,*), SOURCE(NSTOKES, *)
Cf2py intent(in) :: DIRFLUX, FLUXES, SOURCE
      REAL CAMX(*), CAMY(*), CAMZ(*)
      DOUBLE PRECISION CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, NSTPHASE, NSCATANGLE
Cf2py intent(in) :: NPIX, NSTPHASE, NSCATANGLE
      REAL   STOKES(NSTOKES, NPIX)
Cf2py intent(out) :: STOKES
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS
      REAL YLMSUN(NSTLEG, NLM), PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
Cf2py intent(in) :: YLMSUN, PHASETAB

      INTEGER I, J, L, SIDE
      INTEGER IVIS
      LOGICAL VALIDRAD
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION U, R, PI
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT, VISRAD(NSTOKES)

      REAL  MEAN, STD1, STD2

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

      PI = ACOS(-1.0D0)
C         Loop over pixels in image
      DO IVIS = 1, NPIX
        X0 = CAMX(IVIS)
        Y0 = CAMY(IVIS)
        Z0 = CAMZ(IVIS)
        MU2 = CAMMU(IVIS)
        PHI2 = CAMPHI(IVIS)
        MURAY = -MU2
        PHIRAY = PHI2 - PI

C             Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (MURAY .GE. 0.0) THEN
            VISRAD(:) = 0.0
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
     .                       XE,YE,ZE, SIDE, TRANSMIT, VISRAD, VALIDRAD,
     .   	                 TOTAL_EXT, NPART)
900   CONTINUE

C      WRITE (6,'(1X,2F8.4,1X,2F11.7,4(1X,F11.6))')
C     .         X0,Y0,Z0,MURAY,PHIRAY,VISRAD(:)

        STOKES(:, IVIS) = VISRAD(:)
C        IF (VISRAD(1) .GT. 0.0) THEN
C         IF (NSTOKES .GT. 1) THEN
C           Output degree (0 to 1) and angle (-180 to 180) of linear polarization
C           DOLP(IVIS) = SQRT(VISRAD(2)**2+VISRAD(3)**2)/VISRAD(1)
C           AOLP(IVIS) = (180/PI)*0.5*ATAN2(VISRAD(3),VISRAD(2))
C         ENDIF
C         IF (NSTOKES .EQ. 4) THEN
C           Output degree of circular polarization (-1 to 1)
C           DOCP(IVIS) = VISRAD(4)/VISRAD(1)
C         ENDIF
C        ELSE
C          DOLP(IVIS) = 0.0
C          DOCP(IVIS) = 0.0
C          AOLP(IVIS) = 0.0
C        ENDIF
      ENDDO

C      IF (NSTOKES .GT. 1) THEN
C        Choose the best range for the angle of linear polarization (-90 to 90 or 0 to 180)
C        MEAN = SUM(AOLP(:))/NPIX
C        STD1 = SQRT(SUM((AOLP(:)-MEAN)**2)/NPIX)
C        WHERE (AOLP(:) < 0.0)
C          AOLP1(:) = AOLP(:)+180.0
C        END WHERE
C        MEAN = SUM(AOLP1(:))/NPIX
C        STD2 = SQRT(SUM((AOLP1(:)-MEAN)**2)/NPIX)
C        IF (STD2 < STD1) THEN
C          AOLP = AOLP1
C        ENDIF
C      ENDIF

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
      ALLOCATE (LOFJ(NLM))

      J = 0
      DO L = 0, ML
        ME = MIN(L,MM)
        MS = -ME
        DO M = MS, ME
          J = J + 1
          LOFJ(J) = L
        ENDDO
      ENDDO

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

      PI = ACOS(-1.0D0)
C         Loop over pixels in image
      COST = 0.0; NORM1 = 0.0; NORM2 = 0.0; GRAD1 = 0.0; GRAD2 = 0.0
      DO IVIS = 1, NPIX
        X0 = CAMX(IVIS)
        Y0 = CAMY(IVIS)
        Z0 = CAMZ(IVIS)
        MU2 = CAMMU(IVIS)
        PHI2 = CAMPHI(IVIS)
        MURAY = -MU2
        PHIRAY = PHI2 - PI

C             Extrapolate ray to domain top if above
        IF (Z0 .GT. ZGRID(NZ)) THEN
          IF (MURAY .GE. 0.0) THEN
            VISRAD(:) = 0.0
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

C         Integrate the extinction and source function along this ray
C         to calculate the Stokes radiance vector for this pixel
        TRANSMIT = 1.0D0 ; VISRAD = 0.0D0; RAYGRAD = 0.0D0
        CALL GRAD_INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES, NSTLEG,
     .             NSTPHASE, NSCATANGLE, PHASETAB,
     .             NX, NY, NZ, NPTS, NCELLS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             ML, MM, NLM, NLEG, NUMPHASE,
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             DELTAM, SRCTYPE, WAVELEN, SOLARMU,SOLARAZ,
     .             EXTINCT, ALBEDO, LEGEN,
     .             IPHASE, DIRFLUX, SHPTR, SOURCE,
     .             YLMSUN, MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, MU2, PHI2,
     .             X0, Y0, Z0, XE, YE, ZE, SIDE, TRANSMIT, VISRAD,
     .   	       VALIDRAD, TOTAL_EXT, NPART, RAYGRAD, RSHPTR,
     .             RADIANCE, LOFJ, PARTDER, NUMDER, DEXT, DALB,
     .             DIPHASE, DLEG, NBPTS, DNUMPHASE, SOLARFLUX, NPX,
     .             NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER)
900     CONTINUE
        STOKESOUT(:,IVIS) = VISRAD(:)
        DO NS = 1, NSTOKES
          NORM1(NS) = NORM1(NS) + STOKESOUT(NS,IVIS)**2
          NORM2(NS) = NORM2(NS) + MEASUREMENTS(NS,IVIS)**2
          GRAD1(NS,:,:) = GRAD1(NS,:,:)  + STOKESOUT(NS,IVIS) *
     .                                                 RAYGRAD(NS,:,:)
          GRAD2(NS,:,:)  = GRAD2(NS,:,:)  + MEASUREMENTS(NS,IVIS) *
     .                                                 RAYGRAD(NS,:,:)
          COST(NS) = COST(NS) + MEASUREMENTS(NS,IVIS) * VISRAD(NS)
        ENDDO

      ENDDO

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
     .           GRADOUT, COST, MEASUREMENTS, RSHPTR, STOKESOUT,
     .           NPART, TOTAL_EXT, RADIANCE, NUMDER, PARTDER, DEXT,
     .           DALB, DIPHASE, DLEG, NSCATANGLE, YLMSUN, PHASETAB,
     .           NSTPHASE, DPHASETAB, DNUMPHASE, SOLARFLUX, NPX, NPY,
     .           NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS, EXTDIRP,
     .           UNIFORMZLEV, DPATH, DPTR, EXACT_SINGLE_SCATTER,
     .           UNCERTAINTIES, JACOBIAN,MAKEJACOBIAN,
     .           JACOBIANPTR, COUNTER, RAYS_PER_PIXEL, RAY_WEIGHTS,
     .           STOKES_WEIGHTS)
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
      REAL   STOKESOUT(NSTOKES,NPIX)
Cf2py intent(out) :: STOKESOUT
      DOUBLE PRECISION  GRADOUT(NBPTS,NUMDER), COST
Cf2py intent(out) :: GRADOUT, COST, STOKESOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS
      INTEGER NUMDER, PARTDER(NUMDER)
Cf2py intent(in) :: NUMDER, PARTDER
      INTEGER NSCATANGLE, NSTPHASE
      REAL YLMSUN(NSTLEG,NLM), PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
      REAL DPHASETAB(NSTPHASE,DNUMPHASE,NSCATANGLE)
Cf2py intent(in) :: YLMSUN, PHASETAB
Cf2py intent(in) :: NSCATANGLE, YLMSUN, PHASETAB, DPHASETAB, NSTPHASE
      REAL DPATH(8*(NPX+NPY+NPZ),*), WEIGHT
      REAL UNCERTAINTIES(NSTOKES,NSTOKES,*)
      INTEGER DPTR(8*(NPX+NPY+NPZ),*)
Cf2py intent(in) :: DPATH, DPTR, WEIGHTS, UNCERTAINTIES

      REAL JACOBIAN(NSTOKES,NUMDER,*)
Cf2py intent(in,out) :: JACOBIAN
      LOGICAL MAKEJACOBIAN
Cf2py intent(in) :: MAKEJACOBIAN
      INTEGER JI, JACOBIANPTR(2,*)
Cf2py intent(in,out) JACOBIANPTR
      INTEGER COUNTER
Cf2py intent(out) COUNTER
      INTEGER RAYS_PER_PIXEL(*)
Cf2py intent(in) RAYS_PER_PIXEL
      REAL  RAY_WEIGHTS(*), STOKES_WEIGHTS(NSTOKES, *)
Cf2py intent(in) RAY_WEIGHTS, STOKES_WEIGHTS

      DOUBLE PRECISION SYNTHETIC_MEASUREMENT(NSTOKES)
      DOUBLE PRECISION PIXEL_ERROR
      DOUBLE PRECISION RAYGRAD(NSTOKES,NBPTS,NUMDER), VISRAD(NSTOKES)
      DOUBLE PRECISION RAYGRAD_PIXEL(NSTOKES,NBPTS,NUMDER)
      INTEGER I, J, L, SIDE, IVIS
      LOGICAL VALIDRAD
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION U, R, PI
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT
      INTEGER M, ME, MS, NS, NS1
      INTEGER RAY_COUNTER, I2

      INTEGER, ALLOCATABLE :: LOFJ(:)
      ALLOCATE (LOFJ(NLM))

      GRADOUT = 0.0D0
      COUNTER = 1

      J = 0
      DO L = 0, ML
        ME = MIN(L,MM)
        MS = -ME
        DO M = MS, ME
          J = J + 1
          LOFJ(J) = L
        ENDDO
      ENDDO

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

      PI = ACOS(-1.0D0)
C         Loop over pixels in image
      COST = 0.0D0
      DO I = 1, NPIX
        SYNTHETIC_MEASUREMENT = 0.0D0
        RAYGRAD_PIXEL = 0.0D0
        DO I2=1,RAYS_PER_PIXEL(I)
          IVIS = RAY_COUNTER + I2
          X0 = CAMX(IVIS)
          Y0 = CAMY(IVIS)
          Z0 = CAMZ(IVIS)
          MU2 = CAMMU(IVIS)
          PHI2 = CAMPHI(IVIS)
          MURAY = -MU2
          PHIRAY = PHI2 - PI

C             Extrapolate ray to domain top if above
          IF (Z0 .GT. ZGRID(NZ)) THEN
            IF (MURAY .GE. 0.0) THEN
              VISRAD(:) = 0.0
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

C         Integrate the extinction and source function along this ray
C         to calculate the Stokes radiance vector for this pixel
          TRANSMIT = 1.0D0 ; VISRAD = 0.0D0; RAYGRAD = 0.0D0
          CALL GRAD_INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES, NSTLEG,
     .             NSTPHASE, NSCATANGLE, PHASETAB,
     .             NX, NY, NZ, NPTS, NCELLS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             ML, MM, NLM, NLEG, NUMPHASE,
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             DELTAM, SRCTYPE, WAVELEN, SOLARMU,SOLARAZ,
     .             EXTINCT, ALBEDO, LEGEN,
     .             IPHASE, DIRFLUX, SHPTR, SOURCE,
     .             YLMSUN, MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, MU2, PHI2,
     .             X0, Y0, Z0, XE, YE, ZE, SIDE, TRANSMIT, VISRAD,
     .   	       VALIDRAD, TOTAL_EXT, NPART, RAYGRAD, RSHPTR,
     .             RADIANCE, LOFJ, PARTDER, NUMDER, DEXT, DALB,
     .             DIPHASE, DLEG, NBPTS, DNUMPHASE, SOLARFLUX, NPX,
     .             NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER)
  900     CONTINUE

          DO NS=1,NSTOKES
            STOKESOUT(NS,IVIS) = VISRAD(NS)
            SYNTHETIC_MEASUREMENT(NS) = SYNTHETIC_MEASUREMENT(NS) +
     .          VISRAD(NS)*RAY_WEIGHTS(IVIS)*
     .          STOKES_WEIGHTS(NS,IVIS)
            RAYGRAD_PIXEL(NS,:,:) = RAYGRAD_PIXEL(NS,:,:) +
     .          RAYGRAD(NS,:,:)*RAY_WEIGHTS(IVIS)*
     .          STOKES_WEIGHTS(NS,IVIS)
          ENDDO
        ENDDO

        DO NS = 1, NSTOKES
          PIXEL_ERROR = SYNTHETIC_MEASUREMENT(NS) -
     .                    MEASUREMENTS(NS, I)

          DO NS1 = 1, NSTOKES
            WEIGHT = UNCERTAINTIES(NS,NS1,IVIS)
            GRADOUT = GRADOUT + WEIGHT*PIXEL_ERROR*RAYGRAD_PIXEL(NS,:,:)
            COST = COST + 0.5 * WEIGHT*PIXEL_ERROR**2
          ENDDO
        ENDDO
        RAY_COUNTER = RAY_COUNTER + RAYS_PER_PIXEL(I)

        IF (MAKEJACOBIAN .EQV. .TRUE.) THEN
          DO JI = 1, NBPTS
            IF (ANY(ABS(RAYGRAD_PIXEL(1,JI,:))>0.0)) THEN
              JACOBIANPTR(1,COUNTER) = IVIS
              JACOBIANPTR(2,COUNTER) = JI
              JACOBIAN(:,:,COUNTER) = JACOBIAN(:,:,COUNTER) +
     .        RAYGRAD_PIXEL(:,JI,:)
              COUNTER = COUNTER + 1
            ENDIF
          ENDDO
        ENDIF

      ENDDO
      RETURN
      END


      SUBROUTINE GRAD_INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES, NSTLEG,
     .             NSTPHASE, NSCATANGLE, PHASETAB, NX, NY, NZ,
     .             NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR,
     .             CELLFLAGS, XGRID, YGRID, ZGRID, GRIDPOS,
     .             ML, MM, NLM, NLEG, NUMPHASE, NMU, NPHI0MAX,
     .             NPHI0, MU, PHI, WTDO, DELTAM, SRCTYPE,
     .             WAVELEN, SOLARMU, SOLARAZ, EXTINCT, ALBEDO,
     .             LEGEN, IPHASE, DIRFLUX, SHPTR, SOURCE, YLMSUN,
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, MU2, PHI2,
     .             X0, Y0, Z0, XE, YE, ZE, SIDE, TRANSMIT,
     .		       RADOUT, VALIDRAD, TOTAL_EXT, NPART, RAYGRAD,
     .		       RSHPTR, RADIANCE, LOFJ, PARTDER, NUMDER, DEXT,
     .             DALB, DIPHASE, DLEG, NBPTS, DNUMPHASE, SOLARFLUX,
     .             NPX, NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER)

C       Integrates the source function through the extinction field
C     (EXTINCT) backward from the outgoing direction (MU2,PHI2) to find the
C     radiance (RADOUT) at the point X0,Y0,Z0.
C     The transmission and radiance of the ray so far (TRANSMIT, RADOUT)
C     are input and returned after the integration along with the exitting
C     ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
C     5=-Z,6=+Z).
      IMPLICIT NONE
      LOGICAL EXACT_SINGLE_SCATTER
      INTEGER NPX, NPY, NPZ, NBPTS, BCELL
      REAL    DELX, DELY, XSTART, YSTART, SOLARFLUX
      REAL    ZLEVELS(*)
      REAL    EXTDIRP(*)
      DOUBLE PRECISION UNIFORMZLEV
      INTEGER BCFLAG, IPFLAG, NSTOKES, NSTLEG, NSTPHASE, NSCATANGLE
      INTEGER NX, NY, NZ, NPTS, NCELLS, SIDE
      INTEGER ML, MM, NLM, NLEG, NUMPHASE, NPART
      INTEGER MAXNBC, NTOPPTS, NBOTPTS
      INTEGER NMU, NPHI0MAX, NPHI0(*), NSFCPAR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
      INTEGER SHPTR(*), RSHPTR(*)
      INTEGER*2 CELLFLAGS(*)
      INTEGER IPHASE(NPTS,NPART), DNUMPHASE
      INTEGER BCPTR(MAXNBC,2), LOFJ(NLM)
      LOGICAL DELTAM, VALIDRAD, OUTOFDOMAIN
      REAL    WTDO(NMU,*), MU(*), PHI(NMU,*)
      REAL    WAVELEN, SOLARMU, SOLARAZ
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(*)
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(NSTLEG,0:NLEG,*)
      DOUBLE PRECISION RAYGRAD(NSTOKES,NBPTS,NUMDER)
      REAL    DIRFLUX(*), SOURCE(NSTOKES,*), RADIANCE(NSTOKES,*)
      REAL    TOTAL_EXT(*), YLMSUN(NSTLEG,*)
      REAL    PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
      REAL    DPHASETAB(NSTPHASE,DNUMPHASE,NSCATANGLE)
      DOUBLE PRECISION MU2, PHI2, X0,Y0,Z0, XE,YE,ZE
      DOUBLE PRECISION TRANSMIT, RADOUT(NSTOKES)
      CHARACTER SRCTYPE*1, SFCTYPE*2
      REAL      DLEG(NSTLEG,0:NLEG,*), DEXT(NBPTS,NUMDER)
      REAL      DALB(NBPTS,NUMDER)
      REAL      UNSCALED_ALBEDO
      REAL      SINGSCAT8(NSTOKES,8), OSINGSCAT8(NSTOKES,8)
      INTEGER   PARTDER(NUMDER), NUMDER, DIPHASE(NBPTS,NUMDER)
      REAL      SINGSCAT0(NSTOKES,8), SINGSCAT1(NSTOKES,8)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE, SSP
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, J, L, II, IDR, IPA, KK
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
      INTEGER OPPFACE(6), OLDIPTS(8), DONETHIS(8)
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8), GRIDPOINT
      REAL    OEXTINCT8(8), OSRCEXT8(NSTOKES,8)
      REAL    EXTINCT8(8), SRCEXT8(NSTOKES,8)
      REAL    EXT0, EXT1, EXTN
      REAL    SRCEXT0(NSTOKES), SRCEXT1(NSTOKES), RADBND(NSTOKES)
      REAL    XM, YM, F
      REAL    GRAD8(NSTOKES,8,NUMDER), OGRAD8(NSTOKES,8,NUMDER)
      REAL    GRAD0(NSTOKES,8,NUMDER), GRAD1(NSTOKES,8,NUMDER)
      REAL    SRCGRAD(NSTOKES,8,NUMDER), SRCSINGSCAT(NSTOKES,8)
      REAL    DPATH(8*(NPX+NPY+NPZ),*), DEXTM
      INTEGER DPTR(8*(NPX+NPY+NPZ),*)
      DOUBLE PRECISION PI, CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, FC(8), FB(8)
      DOUBLE PRECISION TAUTOL, TAUGRID, S, DELS, TRANSCUT
      DOUBLE PRECISION EXT, TAU, TRANSCELL, ABSCELL
      DOUBLE PRECISION SRC(NSTOKES)
      DOUBLE PRECISION COSSCAT
      REAL YLMDIR(NSTLEG,NLM)
      REAL, ALLOCATABLE ::  SINGSCAT(:,:), DSINGSCAT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: SUNDIRLEG(:)

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/
      DATA DONEFACE/0,0,0,0,0,0,0,0, 0,1,0,3,0,5,0,7, 2,0,4,0,6,0,8,0,
     .              0,0,1,2,0,0,5,6, 3,4,0,0,5,6,0,0,
     .              0,0,0,0,1,2,3,4, 5,6,7,8,0,0,0,0/


C         TRANSCUT is the transmission to stop the integration at
      TRANSCUT = 0.0D0
C5.0E-5
C         TAUTOL is the maximum optical path for the subgrid intervals
      TAUTOL = 0.2

      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)
      PI = ACOS(-1.0D0)
C       Calculate the generalized spherical harmonics for this direction
      CALL YLMALL (.FALSE.,SNGL(MU2),SNGL(PHI2),ML,MM,NSTLEG, YLMDIR)

      IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
C          Get the solar single scattering Stokes vector for the outgoing
C          direction my interpolating in scattering angle in the PHASETAB
C          table and then rotating the Q/U polarization to the outgoing plane.
        COSSCAT = SOLARMU*MU2 + SQRT((1.0-SOLARMU**2)*(1.0-MU2**2))
     .                  *COS(SOLARAZ-PHI2)
        COSSCAT = MAX(MIN(1.0D0, COSSCAT), -1.0D0)
        IF (NUMPHASE .GT. 0) THEN
          ALLOCATE (SINGSCAT(NSTOKES,NUMPHASE))
          ALLOCATE (DSINGSCAT(NSTOKES,DNUMPHASE))
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
          DO I = 1, DNUMPHASE
            DSINGSCAT(1:NSTPHASE,I)
     .               = (1-F)*DPHASETAB(:,I,J) + F*DPHASETAB(:,I,J+1)
            IF (NSTOKES .GT. 1) THEN
              CALL ROTATE_POL_PLANE (NSTOKES, COSSCAT, SOLARMU,
     .                   SNGL(MU2), SOLARAZ-SNGL(PHI2), DSINGSCAT(:,I))
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

        CALL GET_BASE_GRID_CELL(BCELL, ICELL, TREEPTR)

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
          OSINGSCAT8(:,I) = SINGSCAT8(:,I)
          OGRAD8(:,I,:) = GRAD8(:,I,:)
        ENDDO
C         Compute the source function times extinction in direction (MU2,PHI2)


        IF (NSTOKES .EQ. 1) THEN
          CALL COMPUTE_SOURCE_GRAD_1CELL_UNPOL (ICELL, GRIDPTR,
     .            ML, MM, NLM, NLEG, NUMPHASE,
     .            NPTS, DELTAM, SRCTYPE, SOLARMU,
     .            EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .            SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .            DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,
     .            EXTINCT8, SRCEXT8, TOTAL_EXT, NPART, RSHPTR, RADIANCE,
     .            OGRAD8, GRAD8, LOFJ, CELLFLAGS, PARTDER, NUMDER,
     .            DNUMPHASE, DEXT, DALB, DIPHASE, DLEG, NBPTS, BCELL,
     .            DSINGSCAT, SINGSCAT8, OSINGSCAT8)
        ELSE
          CALL COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR,
     .            NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .            NPTS, DELTAM, SRCTYPE, SOLARMU,
     .            EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .            SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .            DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,
     .            EXTINCT8, SRCEXT8, TOTAL_EXT, NPART, RSHPTR, RADIANCE,
     .            OGRAD8, GRAD8, LOFJ, CELLFLAGS, PARTDER, NUMDER,
     .            DNUMPHASE, DEXT, DALB, DIPHASE, DLEG, NBPTS, BCELL,
     .            DSINGSCAT, SINGSCAT8, OSINGSCAT8)
        ENDIF

C         Interpolate the source and extinction to the current point
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XE, YE, ZE, FC)
        SRCEXT1 = FC(1)*SRCEXT8(:,1) + FC(2)*SRCEXT8(:,2) +
     .            FC(3)*SRCEXT8(:,3) + FC(4)*SRCEXT8(:,4) +
     .            FC(5)*SRCEXT8(:,5) + FC(6)*SRCEXT8(:,6) +
     .            FC(7)*SRCEXT8(:,7) + FC(8)*SRCEXT8(:,8)
        SRCEXT1(1) = MAX(0.0, SRCEXT1(1))

        EXT1 = FC(1)*EXTINCT8(1) + FC(2)*EXTINCT8(2) +
     .         FC(3)*EXTINCT8(3) + FC(4)*EXTINCT8(4) +
     .         FC(5)*EXTINCT8(5) + FC(6)*EXTINCT8(6) +
     .         FC(7)*EXTINCT8(7) + FC(8)*EXTINCT8(8)

        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        IF (.NOT. OUTOFDOMAIN) THEN
          CALL GET_INTERP_KERNEL(BCELL,GRIDPTR,GRIDPOS,XE,YE,ZE,FB)
          DO KK=1,8
            GRAD1(:,KK,:) = FB(KK)*GRAD8(:,KK,:)
            SINGSCAT1(:,KK) = FB(KK)*SINGSCAT8(:,KK)
            SINGSCAT1(1,KK) = MAX(0.0, SINGSCAT1(1,KK))
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
     .      MU2,PHI2,XE,YE,ZE,SO,ICELL
          STOP
        ENDIF
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ
C           Find the optical path across the grid cell and figure how
C             many subgrid intervals to use
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XN, YN, ZN, FC)
        EXTN = FC(1)*EXTINCT8(1) + FC(2)*EXTINCT8(2) +
     .         FC(3)*EXTINCT8(3) + FC(4)*EXTINCT8(4) +
     .         FC(5)*EXTINCT8(5) + FC(6)*EXTINCT8(6) +
     .         FC(7)*EXTINCT8(7) + FC(8)*EXTINCT8(8)

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
          CALL GET_INTERP_KERNEL(ICELL,GRIDPTR,GRIDPOS,XI,YI,ZI,FC)
          SRCEXT0(:) = FC(1)*SRCEXT8(:,1) + FC(2)*SRCEXT8(:,2) +
     .                 FC(3)*SRCEXT8(:,3) + FC(4)*SRCEXT8(:,4) +
     .                 FC(5)*SRCEXT8(:,5) + FC(6)*SRCEXT8(:,6) +
     .                 FC(7)*SRCEXT8(:,7) + FC(8)*SRCEXT8(:,8)


          IF (IT .NE. NTAU) THEN
            EXT0 = FC(1)*EXTINCT8(1) + FC(2)*EXTINCT8(2) +
     .             FC(3)*EXTINCT8(3) + FC(4)*EXTINCT8(4) +
     .             FC(5)*EXTINCT8(5) + FC(6)*EXTINCT8(6) +
     .             FC(7)*EXTINCT8(7) + FC(8)*EXTINCT8(8)
          ELSE
            EXT0 = EXTN
          ENDIF
          SRCEXT0(1) = MAX(0.0, SRCEXT0(1))

          IF (.NOT. OUTOFDOMAIN) THEN
            CALL GET_INTERP_KERNEL(BCELL,GRIDPTR,GRIDPOS,XI,YI,ZI,FB)
            DO KK=1, 8
              GRAD0(:,KK,:) = FB(KK)*GRAD8(:,KK,:)
              SINGSCAT0(:,KK) = FB(KK)*SINGSCAT8(:,KK)
              SINGSCAT0(1,KK) = MAX(0.0, SINGSCAT0(1,KK))
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
            SRC(:) = ( 0.5*(SRCEXT0(:)+SRCEXT1(:))
     .           + 0.08333333333*(EXT0*SRCEXT1(:)-EXT1*SRCEXT0(:))*DELS
     .                *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT

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
            SRC(:) = 0.0
            SRCGRAD = 0.0
            SRCSINGSCAT = 0.0
          ENDIF

          RADOUT(:) = RADOUT(:) + TRANSMIT*SRC(:)*ABSCELL
          IF (.NOT. OUTOFDOMAIN) THEN
            DO KK = 1, 8
              GRIDPOINT = GRIDPTR(KK,BCELL)
              RAYGRAD(:,GRIDPOINT,:) = RAYGRAD(:,GRIDPOINT,:) +
     .                     TRANSMIT*SRCGRAD(:,KK,:)*ABSCELL
C             Add gradient component due to the direct solar beam propogation
              IF (EXACT_SINGLE_SCATTER) THEN
                II = 1
                GRIDPOINT = GRIDPTR(KK,BCELL)
                DO WHILE (DPTR(II,GRIDPOINT).GT.0)
                  SSP = DPTR(II,GRIDPOINT)
                  DO IDR = 1, NUMDER
                    IPA = PARTDER(IDR)

                    IF (DELTAM) THEN
                      IF (NUMPHASE .GT. 0) THEN
                        UNSCALED_ALBEDO = ALBEDO(SSP,IPA)/
     .                    (LEGEN(1,ML+1,IPHASE(SSP,IPA))*
     .                    (ALBEDO(SSP,IPA)- 1.0)+ 1.0)
                        DEXTM = (1.0-UNSCALED_ALBEDO*
     .                    LEGEN(1,ML+1,IPHASE(SSP,IPA)))*DEXT(SSP,IDR)
                      ELSE
                        UNSCALED_ALBEDO = ALBEDO(SSP,IPA)/
     .          (LEGEN(1,ML+1,SSP)*(ALBEDO(SSP,IPA)- 1.0)+ 1.0)
                        DEXTM = (1.0-ALBEDO(SSP,IPA)*LEGEN(1,ML+1,SSP))
     .                           *DEXT(SSP,IDR)
                      ENDIF
                    ELSE
                      DEXTM = DEXT(SSP,IDR)
                    ENDIF
                    RAYGRAD(:,SSP,IDR) = RAYGRAD(:,SSP,IDR) -
     .                  DPATH(II,GRIDPOINT)*DEXTM*ABSCELL*TRANSMIT*
     .                  SRCSINGSCAT(:,KK)
                  ENDDO
                  II = II + 1
                ENDDO
              ENDIF
            ENDDO
          ENDIF

          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          SRCEXT1(:) = SRCEXT0(:)
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
          RADOUT(:) = RADOUT(:) + TRANSMIT*RADBND(:)
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
          DEALLOCATE (SINGSCAT, DSINGSCAT)
        ELSE
          DEALLOCATE (SUNDIRLEG)
        ENDIF
      ENDIF
      RETURN
      END



      SUBROUTINE COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR,
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,
     .             EXTINCT8, SRCEXT8, TOTAL_EXT, NPART, RSHPTR,
     .             RADIANCE, OGRAD8, GRAD8, LOFJ, CELLFLAGS, PARTDER,
     .             NUMDER, DNUMPHASE, DEXT, DALB, DIPHASE, DLEG, NBPTS,
     .             BCELL, DSINGSCAT, SINGSCAT8, OSINGSCAT8)
C       Computes the source function times extinction for gridpoints
C     belonging to cell ICELL in the direction (MU,PHI).  The results
C     are returned in SRCEXT8 and EXTINCT8.
C     The spherical harmonic source function series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
      IMPLICIT NONE
      INTEGER ICELL, NSTOKES, NSTLEG, NPTS, ML, MM, NLM, NLEG, NUMPHASE
      INTEGER GRIDPTR(8,*), SHPTR(*), RSHPTR(*)
      INTEGER DONETHIS(8), OLDIPTS(8), BCELL, DNUMPHASE
      INTEGER IPHASE(NPTS,NPART), NPART, NUMDER, NBPTS
      LOGICAL DELTAM
      REAL    SOLARMU
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(NSTLEG,0:NLEG,*), TOTAL_EXT(*)
      REAL    DIRFLUX(*), SOURCE(NSTOKES,*)
      REAL    RADIANCE(NSTOKES,*)
      REAL    YLMDIR(NSTLEG,*), YLMSUN(NSTLEG,*)
      REAL    SINGSCAT(NSTOKES,NUMPHASE)
      REAL    OEXTINCT8(8), OSRCEXT8(NSTOKES,8)
      REAL    EXTINCT8(8), SRCEXT8(NSTOKES,8), W
      DOUBLE PRECISION SUNDIRLEG(0:NLEG), UNSCALED_ALBEDO
      CHARACTER SRCTYPE*1
      INTEGER*2 CELLFLAGS(*)
      LOGICAL OUTOFDOMAIN
      REAL    GRAD8(NSTOKES,8,NUMDER), OGRAD8(NSTOKES,8,NUMDER)
      REAL    SINGSCAT8(NSTOKES,8), OSINGSCAT8(NSTOKES,8)
      REAL    FULL_SINGSCAT(NSTOKES)
      REAL    DLEG(NSTLEG,0:NLEG,DNUMPHASE), DEXT(NBPTS,NUMDER)
      REAL    DALB(NBPTS,NUMDER), DEXTM, DALBM, DLEGM(NSTLEG)
      INTEGER DIPHASE(NBPTS,NUMDER), KK
      REAL    DSINGSCAT(NSTOKES,DNUMPHASE)
      INTEGER LOFJ(*), PARTDER(NUMDER), RNS, RIS, IDR, IB
      INTEGER IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I,IPA,K
      REAL    SECMU0, F, DA, A, A1, B1, EXT

      GRAD8 = 0.0
      SECMU0 = 1.0D0/ABS(SOLARMU)

C         Loop over the grid points, computing the source function
C           at the viewing angle from the spherical harmonics source function.
      DO N = 1, 8
        IP = GRIDPTR(N,ICELL)
        IB = GRIDPTR(N,BCELL)
        I = DONETHIS(N)
        IF (I .GT. 0 .AND. IP .EQ. OLDIPTS(N)) THEN
          EXTINCT8(N) = OEXTINCT8(I)
          SRCEXT8(:,N) = OSRCEXT8(:,I)
          SINGSCAT8(:,N) = OSINGSCAT8(:,I)
          GRAD8(:,N,:) = OGRAD8(:,I,:)
        ELSE IF (I .LT. 0) THEN
          EXTINCT8(N) = EXTINCT8(ABS(I))
          SRCEXT8(:,N) = SRCEXT8(:,ABS(I))
          SINGSCAT8(:,N) = SINGSCAT8(:,ABS(I))
          GRAD8(:,N,:) = GRAD8(:,ABS(I),:)
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
                F = LEGEN(1,ML+1,K)
              ELSE
                F = 0.0
              ENDIF

              UNSCALED_ALBEDO = ALBEDO(IP,IPA)/
     .          (F*(ALBEDO(IP,IPA)- 1.0)+ 1.0)
              DEXTM = (1.0-UNSCALED_ALBEDO*F) * DEXT(IB,IDR)
              DALBM = (1.0-F)*DALB(IB,IDR)/
     .          ((1.0-UNSCALED_ALBEDO*F)**2)

C             Sum over the real generalized spherical harmonic series
C             of the source function
              DO J = 1, RNS
                L = LOFJ(J)
                DLEGM(:) = DLEG(:,L,DIPHASE(IB,IDR))/(1-F)

                GRAD8(1,N,IDR) = GRAD8(1,N,IDR) +
     .            RADIANCE(1,RIS+J)*YLMDIR(1,J)*(
     .              DEXTM*(ALBEDO(IP,IPA)*LEGEN(1,L,K)-1.0) +
     .              EXTINCT(IP,IPA)*DALBM*LEGEN(1,L,K) +
     .              EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(1))

                IF (NSTOKES .GT. 1) THEN
                  GRAD8(1,N,IDR) = GRAD8(1,N,IDR) +
     .              RADIANCE(2,RIS+J)*YLMDIR(1,J)*(
     .              DEXTM*ALBEDO(IP,IPA)*LEGEN(5,L,K) +
     .              EXTINCT(IP,IPA)*DALBM*LEGEN(5,L,K) +
     .              EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(5))

                  IF (J .GE. 5) THEN
                    GRAD8(2,N,IDR) = GRAD8(2,N,IDR) +
     .                RADIANCE(1,RIS+J)*YLMDIR(2,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(5,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(5,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(5))
                    GRAD8(2,N,IDR) = GRAD8(2,N,IDR) +
     .                RADIANCE(2,RIS+J)*YLMDIR(2,J)*(
     .                DEXTM*(ALBEDO(IP,IPA)*LEGEN(2,L,K)-1.0) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(2,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(2))
                    GRAD8(2,N,IDR) = GRAD8(2,N,IDR) +
     .                RADIANCE(3,RIS+J)*YLMDIR(5,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(3,L,K)  +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(3,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(3))

                    GRAD8(3,N,IDR) = GRAD8(3,N,IDR) +
     .                RADIANCE(1,RIS+J)*YLMDIR(6,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(5,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(5,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(5))
                    GRAD8(3,N,IDR) = GRAD8(3,N,IDR) +
     .                RADIANCE(2,RIS+J)*YLMDIR(6,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(2,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(2,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(2))
                    GRAD8(3,N,IDR) = GRAD8(3,N,IDR) +
     .                RADIANCE(3,RIS+J)*YLMDIR(3,J)*(
     .                DEXTM*(ALBEDO(IP,IPA)*LEGEN(3,L,K)-1.0) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(3,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(3))
                  ENDIF
                ENDIF
                IF (NSTOKES .EQ. 4) THEN
                    GRAD8(2,N,IDR) = GRAD8(2,N,IDR) +
     .                RADIANCE(4,RIS+J)*YLMDIR(5,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(6))

                    GRAD8(3,N,IDR) = GRAD8(3,N,IDR) +
     .                RADIANCE(4,RIS+J)*YLMDIR(3,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(6))

                    GRAD8(4,N,IDR) = GRAD8(4,N,IDR) -
     .                RADIANCE(3,RIS+J)*YLMDIR(4,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(6))

                    GRAD8(4,N,IDR) = GRAD8(4,N,IDR) +
     .                RADIANCE(4,RIS+J)*YLMDIR(4,J)*(
     .                DEXTM*(ALBEDO(IP,IPA)*LEGEN(4,L,K)-1.0) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(4,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(4))
                ENDIF
              ENDDO
C             Add in the single scattering contribution for the original unscaled phase function.
C             For DELTA M and L<=ML this requires reconstructing the original phase function values (LEGEN) by unscaling.
C             Need to put the inverse of the tau scaling in the source function because extinction is still scaled.
              IF (NUMPHASE .GT. 0) THEN
                FULL_SINGSCAT(:) = DA*(
     .               DEXTM*ALBEDO(IP,IPA)*SINGSCAT(:,K) +
     .               EXTINCT(IP,IPA)*DALBM*SINGSCAT(:,K) +
     .               EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*
     .               DSINGSCAT(:,DIPHASE(IB,IDR))/(1.0 - F))
              ELSE
                WRITE(*,*) 'NUMPHASE=', NUMPHASE, ' NOT SUPPORTED'
                STOP
              ENDIF
              GRAD8(:,N,IDR) = GRAD8(:,N,IDR) + FULL_SINGSCAT(:)
            ENDDO
          ENDIF

C             Sum over the real generalized spherical harmonic series
C             of the source function
          SRCEXT8(:,N) = 0.0
          SINGSCAT8(:,N) = 0.0
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

C
          DO IPA = 1, NPART

            IF (EXT.EQ.0.0) THEN
              W = 1.0
            ELSE
              W = EXTINCT(IP,IPA)/EXT
            ENDIF
            IF (W.EQ.0.0) CYCLE

            IF (NUMPHASE .GT. 0) THEN
              IPH = IPHASE(IP,IPA)
            ELSE
              IPH = IP
            ENDIF

C             First subtract off the truncated single scattering
            IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
              DA = ALBEDO(IP,IPA)*DIRFLUX(IP)*SECMU0*W
              J = 1

              DO L = 0, ML
                ME = MIN(L,MM)
                MS = -ME
                A1 = DA*LEGEN(1,L,IPH)
                B1 = DA*LEGEN(5,L,IPH)
                IF (J .LE. NS) THEN
                  JT = J
                  DO M = MS, ME
                    SRCEXT8(1,N) =SRCEXT8(1,N) -
     .                                       A1*YLMDIR(1,J)*YLMSUN(1,J)
                    J = J + 1
                  ENDDO
                  IF (NSTOKES .GT. 1) THEN
                    J = JT
                    DO M = MS, ME
                      SRCEXT8(2,N)=SRCEXT8(2,N) -
     .                                       B1*YLMDIR(2,J)*YLMSUN(1,J)
                      SRCEXT8(3,N)=SRCEXT8(3,N) -
     .                                       B1*YLMDIR(6,J)*YLMSUN(1,J)
                      J = J + 1
                    ENDDO
                  ENDIF
                ENDIF
              ENDDO
            ENDIF

C             Then add in the single scattering contribution for the
C             original unscaled phase function.
            IF (NUMPHASE .GT. 0) THEN
              SINGSCAT8(:,N) = SINGSCAT8(:,N) + DA*SINGSCAT(:,IPH)
              IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
                SRCEXT8(:,N) = SRCEXT8(:,N) + DA*SINGSCAT(:,IPH)
              ENDIF
            ELSE IF (NSTOKES .EQ. 1) THEN
              F = LEGEN(1,ML+1,IPH)
              DO L = 0, NLEG
                IF (L .LE. ML) THEN
                  A = DA*(LEGEN(1,L,IPH) + F/(1-F))
                ELSE
                  A = DA*LEGEN(1,L,IPH)/(1-F)
                ENDIF
                SINGSCAT8(1,N) = SINGSCAT8(1,N) + A*SUNDIRLEG(L)
                IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
                  SRCEXT8(1,N) = SRCEXT8(1,N) + A*SUNDIRLEG(L)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          SINGSCAT8(:,N) = SINGSCAT8(:,N)*EXT
          SRCEXT8(:,N) = SRCEXT8(:,N)*EXT
          EXTINCT8(N) = EXT
        ENDIF
      ENDDO

      RETURN
      END



      SUBROUTINE COMPUTE_SOURCE_GRAD_1CELL_UNPOL (ICELL, GRIDPTR,
     .             ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,
     .             EXTINCT8, SRCEXT8, TOTAL_EXT, NPART, RSHPTR,
     .             RADIANCE, OGRAD8, GRAD8, LOFJ, CELLFLAGS, PARTDER,
     .             NUMDER, DNUMPHASE, DEXT, DALB, DIPHASE, DLEG, NBPTS,
     .             BCELL, DSINGSCAT, SINGSCAT8, OSINGSCAT8)
C       Computes the source function times extinction for gridpoints
C     belonging to cell ICELL in the direction (MU,PHI).  The results
C     are returned in SRCEXT8 and EXTINCT8.
C     The spherical harmonic source function series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
      IMPLICIT NONE
      INTEGER ICELL, NPTS, ML, MM, NLM, NLEG, NUMPHASE, BCELL
      INTEGER GRIDPTR(8,*), SHPTR(NPTS+1), NBPTS, RSHPTR(*)
      INTEGER DONETHIS(8), OLDIPTS(8), NUMDER
      INTEGER IPHASE(NPTS,NPART), NPART, DNUMPHASE
      LOGICAL DELTAM
      REAL    SOLARMU
      INTEGER*2 CELLFLAGS(*)
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(0:NLEG,*),  TOTAL_EXT(NPTS)
      REAL    DIRFLUX(*), SOURCE(*), RADIANCE(*)
      REAL    YLMDIR(*), YLMSUN(*)
      REAL    SINGSCAT(NUMPHASE)
      REAL    OEXTINCT8(8), OSRCEXT8(8)
      REAL    EXTINCT8(8), SRCEXT8(8)
      REAL    GRAD8(8,NUMDER), OGRAD8(8,NUMDER)
      REAL    FULL_SINGSCAT, SINGSCAT8(8), OSINGSCAT8(8)
      REAL    DLEG(0:NLEG,DNUMPHASE), DEXT(NBPTS,NUMDER)
      REAL    DALB(NBPTS,NUMDER), DEXTM, DALBM, DLEGM
      LOGICAL OUTOFDOMAIN
      REAL    DSINGSCAT(DNUMPHASE), W
      INTEGER DIPHASE(NBPTS,NUMDER), KK
      DOUBLE PRECISION SUNDIRLEG(0:NLEG), UNSCALED_ALBEDO
      CHARACTER SRCTYPE*1
      INTEGER LOFJ(NLM), PARTDER(NUMDER), RNS, RIS, IDR, IB
      INTEGER IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, IPA, I, K
      REAL    SECMU0, F, DA, A, A1, B1, EXT

      SECMU0 = 1.0D0/ABS(SOLARMU)
      GRAD8 = 0.0

C         Loop over the grid points, computing the source function
C           at the viewing angle from the spherical harmonics source function.
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
              UNSCALED_ALBEDO = ALBEDO(IP,IPA)/
     .          (F*(ALBEDO(IP,IPA)- 1.0)+ 1.0)
               DEXTM = (1.0-UNSCALED_ALBEDO*F) * DEXT(IB,IDR)
               DALBM = (1.0-F)*DALB(IB,IDR)/
     .          ((1.0-UNSCALED_ALBEDO*F)**2)

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
     .               DEXTM*ALBEDO(IP,IPA)*A + EXTINCT(IP,IPA)*DALBM*A +
     .               EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM)
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
     .               DSINGSCAT(DIPHASE(IB,IDR))/(1.0-F))
              ENDIF
              GRAD8(N,IDR) = GRAD8(N,IDR) + FULL_SINGSCAT
            ENDDO
          ENDIF

C             Sum over the real generalized spherical harmonic series
C             of the source function
          SRCEXT8(N) = 0.0
          DO J = 1, NS
            SRCEXT8(N) = SRCEXT8(N) + SOURCE(IS+J)*YLMDIR(J)
          ENDDO

C             Special case for solar source and Delta-M
          IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
            DO IPA = 1, NPART
              IF (EXT.EQ.0.0) THEN
                W = 1.0
              ELSE
                W = EXTINCT(IP,IPA)/EXT
              ENDIF
              IF (W.EQ.0.0) CYCLE
              IF (NUMPHASE .GT. 0) THEN
                IPH = IPHASE(IP,IPA)
              ELSE
                IPH = IP
              ENDIF
C               First subtract off the truncated single scattering
            DA = ALBEDO(IP,IPA)*DIRFLUX(IP)*SECMU0*W
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
            ENDDO
          ENDIF

          SRCEXT8(N) = SRCEXT8(N)*EXT
          EXTINCT8(N) = EXT
        ENDIF
      ENDDO

      RETURN
      END



      SUBROUTINE RAYLEIGH_PHASE_FUNCTION (WAVELEN, RAYLEGCOEF,
     .                                    TABLE_TYPE)
C     Returns the Wigner d-function coefficients for either the Rayleigh
C     scalar phase function or polarized phase matrix for molecular
C     scattering by air.  Includes the wavelength depolarization factor.
C     From Mishchenko's book "Multiple scattering of light by particles:
C     Radiative Transfer and Coherent Backscattering", Cambridge, 2006.
C     Thanks to Adrian Doicu.   WAVELEN is the wavelength in microns.
      IMPLICIT NONE
      REAL      WAVELEN
Cf2py intent(in) :: WAVELEN
      INTEGER   NCOEF
      PARAMETER (NCOEF=2)
      REAL      RAYLEGCOEF(6, 0:NCOEF)
Cf2py intent(out) :: RAYLEGCOEF
      CHARACTER(LEN=6) :: TABLE_TYPE
Cf2py intent(out) ::  TABLE_TYPE
      DOUBLE PRECISION AKING, BKING, CKING
      PARAMETER (AKING=1.0469541D0,
     .           BKING=3.2503153D-04, CKING=3.8622851D-05)
      REAL :: FKING, DEPOL, DELTA, DELTAP

      TABLE_TYPE = 'VECTOR'
      RAYLEGCOEF(1:6,:) = 0.0
      FKING = AKING + BKING/WAVELEN**2 + CKING/WAVELEN**4
      DEPOL = 6.D0*(FKING-1.D0) / (3.D0 + 7.D0*FKING)

      DELTA = (1.0 - DEPOL) / (1.0 + 0.5*DEPOL)
      DELTAP = (1.0 - 2.0*DEPOL) / (1.0 - DEPOL)
      RAYLEGCOEF(1,0) = 1.0  ; RAYLEGCOEF(1,2) = 0.5*DELTA
      RAYLEGCOEF(2,2) = 3.0*DELTA
      RAYLEGCOEF(4,1) = 1.5*DELTAP*DELTA
      RAYLEGCOEF(5,2) = SQRT(1.5)*DELTA

      RETURN
      END


      SUBROUTINE PRECOMPUTE_PHASE_CHECK(NSCATANGLE, NUMPHASE, NSTPHASE,
     .                              NSTOKES, ML, NLM, NSTLEG, NLEG,
     .                              LEGEN, PHASETAB, DELTAM, NEGCHECK)
C       Precomputes the phase matrix elements I-I and I-Q as a function
C     of scattering angle for solar direct scattering for all the
C     tabulated phase functions. Output is in PHASETAB.
C     NSTPHASE=1 for NSTOKES=1 and NSTPHASE=2 for NSTOKES>1.
      IMPLICIT NONE
      INTEGER NSCATANGLE, NSTPHASE, NSTOKES, NSTLEG
      INTEGER ML, NLM, NLEG, NUMPHASE
Cf2py intent(in) :: NSCATANGLE, NSTPHASE, NSTOKES, NSTLEG, ML, NLM, NLEG, NUMPHASE
      REAL    LEGEN(NSTLEG,0:NLEG,NUMPHASE)
Cf2py intent(in) :: LEGEN
      REAL    PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
Cf2py intent(out) :: PHASETAB
      LOGICAL DELTAM, NEGCHECK
Cf2py intent(in) :: DELTAM, NEGCHECK
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
            IF (L .LE. ML .AND. DELTAM) THEN
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH) + F/(1-F)
            ELSE IF (DELTAM) THEN
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)/(1-F)
            ELSE
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)
            ENDIF
            IF (NSTLEG .GT. 1) THEN
              IF (L .LE. ML .AND. DELTAM) THEN
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)
              ELSE IF (DELTAM) THEN
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)/(1-F)
              ELSE
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)
              ENDIF
            ENDIF
          ENDDO

C           Sum the first Wigner function series, which is actually a Legendre series
          A1 = 0.0D0
          DO L = 0, NLEG
            FCT = 2.0D0*L + 1.0D0
            A1  = A1 + FCT*DBLE(UNSCLEGEN(1,L))*DMM1(L)
          ENDDO
          IF (NEGCHECK .AND. A1 .LE. 0.0) THEN
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


      SUBROUTINE PRECOMPUTE_PHASE_CHECK_GRAD(NSCATANGLE, DNUMPHASE,
     .                              NSTPHASE,
     .                              NSTOKES, ML, NLM, NSTLEG, NLEG,
     .                              DLEG, DPHASETAB, DELTAM, NEGCHECK)
C       Precomputes the phase matrix elements I-I and I-Q as a function
C     of scattering angle for solar direct scattering for all the
C     tabulated phase functions. Output is in PHASETAB.
C     NSTPHASE=1 for NSTOKES=1 and NSTPHASE=2 for NSTOKES>1.
      IMPLICIT NONE
      INTEGER NSCATANGLE, NSTPHASE, NSTOKES, NSTLEG
      INTEGER ML, NLM, NLEG, DNUMPHASE
Cf2py intent(in) :: NSCATANGLE, NSTPHASE, NSTOKES, NSTLEG, ML, NLM, NLEG, DNUMPHASE
      REAL    DLEG(NSTLEG,0:NLEG,DNUMPHASE)
Cf2py intent(in) :: DLEG
      REAL    DPHASETAB(NSTPHASE,DNUMPHASE,NSCATANGLE)
Cf2py intent(out) :: DPHASETAB
      LOGICAL DELTAM, NEGCHECK
Cf2py intent(in) :: DELTAM, NEGCHECK
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
        DO IPH = 1, DNUMPHASE
C          F = LEGEN(1,ML+1,IPH)
          DO L = 0, NLEG
            IF (L .LE. ML .AND. DELTAM) THEN
              UNSCLEGEN(1,L) = DLEG(1,L,IPH)
            ELSE IF (DELTAM) THEN
              UNSCLEGEN(1,L) = DLEG(1,L,IPH)
            ELSE
              UNSCLEGEN(1,L) = DLEG(1,L,IPH)
            ENDIF
            IF (NSTLEG .GT. 1) THEN
              IF (L .LE. ML .AND. DELTAM) THEN
                UNSCLEGEN(2,L) = DLEG(5,L,IPH)
              ELSE IF (DELTAM) THEN
                UNSCLEGEN(2,L) = DLEG(5,L,IPH)
              ELSE
                UNSCLEGEN(2,L) = DLEG(5,L,IPH)
              ENDIF
            ENDIF
          ENDDO

C           Sum the first Wigner function series, which is actually a Legendre series
          A1 = 0.0D0
          DO L = 0, NLEG
            FCT = 2.0D0*L + 1.0D0
            A1  = A1 + FCT*DBLE(UNSCLEGEN(1,L))*DMM1(L)
          ENDDO
          IF (NEGCHECK .AND. A1 .LE. 0.0) THEN
            WRITE (6,*) 'PRECOMPUTE_PHASE: negative phase function',
     .          ' for tabulated phase function: ',IPH, J, A1
            STOP
          ENDIF
          DPHASETAB(1,IPH,J) = SNGL(A1*OFOURPI)

          IF (NSTOKES .GT. 1) THEN
C           If doing polarization, sum the second Wigner function series
            B1 = 0.0D0
            DO L = 0, NLEG
              FCT = 2.0D0*L + 1.0D0
              B1  = B1 - FCT*DBLE(UNSCLEGEN(2,L))*DMM2(L)
            ENDDO
            DPHASETAB(2,IPH,J) = SNGL(B1*OFOURPI)
          ENDIF
        ENDDO
      ENDDO

      DEALLOCATE (UNSCLEGEN, DMM1, DMM2)
      RETURN
      END
