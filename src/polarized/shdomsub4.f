C     This file contains fortran subroutines adapted from subroutines written
C     Frank Evans for SHDOM.
C     https://nit.coloradolinux.com/~evans/shdom.html
C     See shdom.txt for documentation of SHDOM variables.
C     These subroutines have been written for use in pyshdom by
C     Aviad Levis, Technion Institute of Technology, 2019 and
C     later modified by Jesse Loveridge, University of Illinois at Urbana-Champaign,
C     2020-2021.
C     This file contains the subroutines that are used to evaluate cost functions
C     gradients (using the Levis approximation) for solutions to retrieval problems.
C     -JRLoveridge 2021/02/22

      SUBROUTINE UPDATE_COSTFUNCTION (STOKESOUT, RAYGRAD_PIXEL,
     .             GRADOUT,COST, UNCERTAINTIES, COSTFUNC, NSTOKES,
     .             MAXPG, NUMDER, NCOST, NGRAD, MEASUREMENT,
     .             NUNCERTAINTY, IERR, ERRMSG)
C    Updates the cost function (COST) and its gradient (GRADOUT)
C    with the contribution from one pixel using the forward modeled StokesVector
C    (STOKESOUT), and the Frechet derivative for that pixel (RAYGRAD_PIXEL)
C    and the inverse error co-variance for this pixel which should have been
C    prepared in uncertainties.py to match the choice of observables/cost function
C    (COSTFUNC).
C    To add a new cost function, add a new block here for a new value
C    of COSTFUNC, add a test to verify it works as intended to tests/test_derivatives.py,
C    and update gradient.py and uncertainties.py to accept the new gradient flag.
      CHARACTER*2 COSTFUNC
Cf2py intent(in) :: COSTFUNC
      INTEGER MAXPG, NUMDER, NSTOKES, NGRAD, NCOST
      INTEGER NUNCERTAINTIY
      DOUBLE PRECISION COST(NCOST)
      DOUBLE PRECISION GRADOUT(MAXPG,NUMDER,NGRAD)
Cf2py intent(in,out) :: COST, GRADOUT
      DOUBLE PRECISION RAYGRAD_PIXEL(NSTOKES,MAXPG,NUMDER)
Cf2py intent(in) :: RAYGRAD_PIXEL
      DOUBLE PRECISION STOKESOUT(NSTOKES), MEASUREMENT(NSTOKES)
Cf2py intent(in) :: STOKESOUT, MEASUREMENT
      DOUBLE PRECISION UNCERTAINTIES(NUNCERTAINTY,NUNCERTAINTY)
Cf2py intetn(in) :: UNCERTAINTIES
      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG

      DOUBLE PRECISION PIXEL_ERROR
      INTEGER I, J
      IERR = 0
      IF (COSTFUNC .EQ. 'L2') THEN
C       A weighted least squares on each of the stokes components.
        DO I=1,NSTOKES
          PIXEL_ERROR = STOKESOUT(I) - MEASUREMENT(I)
          DO J=1,NSTOKES
            COST(1) = COST(1) + 0.5D0*
     .        UNCERTAINTIES(I,J) * PIXEL_ERROR**2
            GRADOUT(:,:,1) = GRADOUT(:,:,1) + UNCERTAINTIES(I,J)*
     .        PIXEL_ERROR*RAYGRAD_PIXEL(I,:,:)
           ENDDO
         ENDDO

      ELSEIF (COSTFUNC .EQ. 'LL') THEN
C       Logs of radiance and degree of linear polarization are observables
C       which is standard in atmospheric remote sensing.
C       e.g. Dubovik et al. 2011 https://doi.org/10.5194/amt-4-975-2011.
        RADERROR = LOG(STOKESOUT(1)) - LOG(MEASUREMENT(1))

        COST(1) = COST(1) + 0.5D0 *(RADERROR**2 *
     .      UNCERTAINTIES(1,1))
        GRADOUT(:,:,1) = GRADOUT(:,:,1) + RADERROR*
     .      UNCERTAINTIES(1,1)*RAYGRAD_PIXEL(1,:,:)/STOKESOUT(1)

        IF (NSTOKES .GT. 1) THEN
          DOLP1 = SQRT(STOKESOUT(2)**2 + STOKESOUT(3)**2) /
     .      STOKESOUT(1)
          DOLP2 = SQRT(MEASUREMENT(2)**2 + MEASUREMENT(3)**2) /
     .      MEASUREMENT(1)
          DOLPERR = LOG(DOLP1) - LOG(DOLP2)

          COST(1) = COST(1) + 0.5D0 *(DOLPERR**2 *
     .      UNCERTAINTIES(2,2))
          GRADOUT(:,:,1) = GRADOUT(:,:,1) + DOLPERR*
     .      UNCERTAINTIES(2,2)*(STOKESOUT(2)*
     .      RAYGRAD_PIXEL(2,:,:) +
     .      STOKESOUT(3)*RAYGRAD_PIXEL(3,:,:))/
     .      (STOKESOUT(2)**2 + STOKESOUT(3)**2)
        ENDIF
      ELSE
        IERR = 1
        WRITE(ERRMSG,*) 'UPDATE_COSTFUNCTION: Bad cost function',
     .    COSTFUNC
        RETURN
      ENDIF
      RETURN
      END

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
     .                   NSCATANGLE, YLMSUN, PHASETAB, NSTPHASE,
     .                  IERR, ERRMSG, INTERPMETHOD, PHASEINTERPWT,
     .                   PHASEMAX)
C    Calculates the Stokes Vector at the given directions (CAMMU, CAMPHI)
C    and positions CAMX,CAMY,CAMZ by integrating the source function.

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
      INTEGER IPHASE(8,NPTS,NPART)
Cf2py intent(in) :: CELLFLAGS, IPHASE
      REAL PHASEINTERPWT(8,NPTS,NPART), PHASEMAX
Cf2py intent(in) :: PHASEINTERPWT, PHASEMAX
      CHARACTER INTERPMETHOD*2
Cf2py intent(in) :: INTERPMETHOD
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
      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG
      INTEGER I, J, L, SIDE
      INTEGER IVIS
      LOGICAL VALIDRAD
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION U, R, PI
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT, VISRAD(NSTOKES)

      REAL  MEAN, STD1, STD2
      IERR = 0
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
          IERR = 1
          WRITE (ERRMSG,*) 'RENDER: Level below domain'
          RETURN
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
     .   	                 TOTAL_EXT, NPART, IERR,ERRMSG, INTERPMETHOD,
     .                     PHASEINTERPWT, PHASEMAX)
      IF (IERR .NE. 0) RETURN
900   CONTINUE


        STOKES(:, IVIS) = VISRAD(:)
      ENDDO

      RETURN
      END

      SUBROUTINE LEVISAPPROX_GRADIENT(NSTOKES, NX, NY, NZ, NPTS,
     .           NCELLS,MAXPG,
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
     .           JACOBIANPTR, NUM_JACOBIAN_PTS, RAYS_PER_PIXEL,
     .           RAY_WEIGHTS, STOKES_WEIGHTS, DIPHASEIND,
     .           COSTFUNC, NCOST, NGRAD, NUNCERTAINTY, PLANCK,
     .           TAUTOL, NODIFFUSE, IERR, ERRMSG)
C    Calculates the cost function and its gradient using the Levis approximation
C    to the Frechet derivatives of the radiative transfer equation.
C    Calculates the Stokes Vector at the given directions (CAMMU, CAMPHI)
C    and positions CAMX,CAMY,CAMZ by integrating the source function for evaluation
C    of the cost function while simultaneously calculating the gradient.
C    MEASUREMENTS and STOKESOUT are used to evaluate the cost function.
C    First loops through pixels and rays in each pixel and calculates ray contributions
C    to pixel quantities weighted by RAY_WEIGHTS.
C    DEXT, DALB, DIPHASE, DLEG, DPHASETAB provide partial derivatives of
C    (possibly delta-M scaled) optical properties with respect to the unknowns.
C    for evaluation of the Frechet derivatives.
C    UNCERTAINTIES holds the inverse error-covariance matrix or other weighting matrix
C    for evaluation of the cost function.
C    Will also output specific Frechet derivative values if MAKEJACOBIAN is TRUE

Cf2py threadsafe
      IMPLICIT NONE
      CHARACTER COSTFUNC*2
Cf2py intent(in) :: COSTFUNC
      INTEGER NCOST, NGRAD
Cf2py intent(in) :: NCOST, NGRAD
      LOGICAL EXACT_SINGLE_SCATTER
Cf2py intent(in) ::  EXACT_SINGLE_SCATTER
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS
      INTEGER MAXPG, NCELLS, NBCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG,
Cf2py intent(in) :: NPTS, NCELLS, MAXPG, NBCELLS
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
      REAL    PLANCK(NPTS,NPART)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN, TOTAL_EXT
Cf2py intent(in) :: PLANCK
      REAL    DIRFLUX(*), FLUXES(2,*)
      REAL    SOURCE(NSTOKES, *), RADIANCE(NSTOKES,*)
Cf2py intent(in) :: DIRFLUX, FLUXES, SOURCE, RADIANCE
      REAL CAMX(*), CAMY(*), CAMZ(*)
      DOUBLE PRECISION CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX
      REAL   MEASUREMENTS(NSTOKES,*), DLEG(NSTLEG,0:NLEG,DNUMPHASE)
      REAL   DEXT(MAXPG,NUMDER), DALB(MAXPG,NUMDER)
      INTEGER DIPHASE(MAXPG,NUMDER)
Cf2py intent(in) :: MEASUREMENTS, DEXT ,DALB, DIPHASE, DLEG
      REAL  STOKESOUT(NSTOKES,NPIX)
Cf2py intent(out) :: STOKESOUT
      DOUBLE PRECISION  GRADOUT(MAXPG,NUMDER,NGRAD), COST(NCOST)
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
      REAL DPATH(8*(NPX+NPY+NPZ),*)
      INTEGER :: NUNCERTAINTY
Cf2py intent(in) :: NUNCERTAINTY
      DOUBLE PRECISION UNCERTAINTIES(NUNCERTAINTY,NUNCERTAINTY,*)
      INTEGER DPTR(8*(NPX+NPY+NPZ),*)
Cf2py intent(in) :: DPATH, DPTR, UNCERTAINTIES
      INTEGER DIPHASEIND(NPTS,NUMDER)
Cf2py intent(in) :: DIPHASEIND
      REAL JACOBIAN(NSTOKES,NUMDER,NUM_JACOBIAN_PTS,*)
Cf2py intent(in,out) :: JACOBIAN
      LOGICAL MAKEJACOBIAN
Cf2py intent(in) :: MAKEJACOBIAN
      INTEGER JI,JJ,JK, NUM_JACOBIAN_PTS, JACOBIANPTR(NUM_JACOBIAN_PTS)
Cf2py intent(in) :: NUM_JACOBIAN_PTS, JACOBIANPTR
      INTEGER RAYS_PER_PIXEL(*)
Cf2py intent(in) :: RAYS_PER_PIXEL
      DOUBLE PRECISION   RAY_WEIGHTS(*), STOKES_WEIGHTS(NSTOKES, *)
Cf2py intent(in) :: RAY_WEIGHTS, STOKES_WEIGHTS
      DOUBLE PRECISION TAUTOL
Cf2py intent(in) :: TAUTOL
      LOGICAL NODIFFUSE
Cf2py intent(in) :: NODIFFUSE
      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG


      DOUBLE PRECISION WEIGHT
      DOUBLE PRECISION PIXEL_ERROR
      DOUBLE PRECISION RAYGRAD(NSTOKES,MAXPG,NUMDER), VISRAD(NSTOKES)
      DOUBLE PRECISION RAYGRAD_PIXEL(NSTOKES,MAXPG,NUMDER)
      INTEGER IPIX, J, L, SIDE, IRAY
      LOGICAL VALIDRAD
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION U, R, PI
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT
      INTEGER M, ME, MS, NS, NS1
      INTEGER I2

      INTEGER, ALLOCATABLE :: LOFJ(:)
      ALLOCATE (LOFJ(NLM))
      IERR = 0
      GRADOUT = 0.0D0
      STOKESOUT = 0.0D0

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
      IRAY = 0
      DO IPIX = 1, NPIX
        RAYGRAD_PIXEL = 0.0D0
        DO I2=1 ,RAYS_PER_PIXEL(IPIX)
          IRAY = IRAY + 1
          X0 = CAMX(IRAY)
          Y0 = CAMY(IRAY)
          Z0 = CAMZ(IRAY)
          MU2 = CAMMU(IRAY)
          PHI2 = CAMPHI(IRAY)
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
            IERR = 1
            WRITE (ERRMSG,*) 'LEVISAPPROX_GRADIENT: Level below domain'
            RETURN
          ENDIF
C         Integrate the extinction and source function along this ray
C         to calculate the Stokes radiance vector for this pixel.
C         Simultaneously calculate the approximate Frechet derivatives
C         while traversing the SHDOM grid.
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
     .             DIPHASE, DLEG, MAXPG, DNUMPHASE, SOLARFLUX, NPX,
     .             NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER, DIPHASEIND, PLANCK,
     .             TAUTOL, NODIFFUSE,
     .             GNDALBEDO, IERR, ERRMSG)
          IF (IERR .NE. 0) RETURN
  900     CONTINUE
          DO NS=1,NSTOKES
            STOKESOUT(NS,IPIX) = STOKESOUT(NS,IPIX) + VISRAD(NS)*
     .        RAY_WEIGHTS(IRAY)*STOKES_WEIGHTS(NS,IPIX)
            RAYGRAD_PIXEL(NS,:,:) = RAYGRAD_PIXEL(NS,:,:) +
     .        RAYGRAD(NS,:,:)*RAY_WEIGHTS(IRAY)*STOKES_WEIGHTS(NS,IPIX)
          ENDDO
        ENDDO
        CALL UPDATE_COSTFUNCTION(DBLE(STOKESOUT(:,IPIX)), RAYGRAD_PIXEL,
     .             GRADOUT, COST, UNCERTAINTIES(:,:,IPIX), COSTFUNC,
     .             NSTOKES, MAXPG, NUMDER, NCOST, NGRAD,
     .             DBLE(MEASUREMENTS(:,IPIX)), NUNCERTAINTY, IERR,
     .              ERRMSG)
        IF (IERR .NE. 0) RETURN

        IF (MAKEJACOBIAN .EQV. .TRUE.) THEN
          DO JI = 1,NUM_JACOBIAN_PTS
            JACOBIAN(:,:,JI,IPIX) =
     .      RAYGRAD_PIXEL(:,JACOBIANPTR(JI),:)
          END DO
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
     .             DALB, DIPHASE, DLEG, MAXPG, DNUMPHASE, SOLARFLUX,
     .             NPX, NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER, DIPHASEIND, PLANCK,
     .              TAUTOL, NODIFFUSE,
     .             GNDALBEDO, IERR, ERRMSG)
C       Integrates the source function through the extinction field
C     (EXTINCT) backward from the outgoing direction (MU2,PHI2) to find the
C     radiance (RADOUT) at the point X0,Y0,Z0.
C     The transmission and radiance of the ray so far (TRANSMIT, RADOUT)
C     are input and returned after the integration along with the exitting
C     ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
C     5=-Z,6=+Z).
C     Updates RAYGRAD with the approximate Frechet derivatives calculated using
C     the partial derivatives DEXT, DALB, DIPHASE, DLEG, DPHASETAB.

      IMPLICIT NONE
      LOGICAL EXACT_SINGLE_SCATTER
      INTEGER NPX, NPY, NPZ, MAXPG, BCELL
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
      REAL    LEGEN(NSTLEG,0:NLEG,*), PLANCK(NPTS,NPART)
      DOUBLE PRECISION RAYGRAD(NSTOKES,MAXPG,NUMDER)
      REAL    DIRFLUX(*), SOURCE(NSTOKES,*), RADIANCE(NSTOKES,*)
      REAL    TOTAL_EXT(*), YLMSUN(NSTLEG,*)
      REAL    PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
      REAL    DPHASETAB(NSTPHASE,DNUMPHASE,NSCATANGLE)
      DOUBLE PRECISION MU2, PHI2, X0,Y0,Z0, XE,YE,ZE
      DOUBLE PRECISION TRANSMIT, RADOUT(NSTOKES)
      CHARACTER SRCTYPE*1, SFCTYPE*2
      REAL      DLEG(NSTLEG,0:NLEG,*), DEXT(MAXPG,NUMDER)
      REAL      DALB(MAXPG,NUMDER)
      REAL      UNSCALED_ALBEDO
      REAL      SINGSCAT8(NSTOKES,8), OSINGSCAT8(NSTOKES,8)
      INTEGER   PARTDER(NUMDER), NUMDER, DIPHASE(MAXPG,NUMDER)
      INTEGER   DIPHASEIND(NPTS,NUMDER)
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
      REAL    GRAD8(NSTOKES,8,8,NUMDER), OGRAD8(NSTOKES,8,8,NUMDER)
      DOUBLE PRECISION BASEADAPTINTERP(8,8)
      DOUBLE PRECISION OBASEADAPTINTERP(8,8)
      REAL GRAD0(NSTOKES,8,NUMDER), GRAD1(NSTOKES,8,NUMDER)
      REAL    SRCGRAD(NSTOKES,8,NUMDER), SRCSINGSCAT(NSTOKES,8)
      REAL    DPATH(8*(NPX+NPY+NPZ),*), DEXTM, SECMU0
      INTEGER DPTR(8*(NPX+NPY+NPZ),*), N
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
      LOGICAL :: NODIFFUSE
      DOUBLE PRECISION :: DIRRAD(NSTOKES,4), BOUNDINTERP(4)
      INTEGER :: BOUNDPTS(4), IB,IP
      REAL :: GNDALBEDO
      INTEGER IERR
      CHARACTER ERRMSG*600

      INTEGER, ALLOCATABLE :: PASSEDPOINTS(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDRAD(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDABSCELL(:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDTRANSMIT(:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDINTERP0(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDINTERP1(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDDELS(:)
      INTEGER MAXSUBGRIDINTS, NPASSED, K, MAXBYOPT
      REAL MAXTAU, DZMAX
      DOUBLE PRECISION  RADGRAD(NSTOKES), RAD0(NSTOKES)
      DOUBLE PRECISION  RAD1(NSTOKES)

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/
      DATA DONEFACE/0,0,0,0,0,0,0,0, 0,1,0,3,0,5,0,7, 2,0,4,0,6,0,8,0,
     .              0,0,1,2,0,0,5,6, 3,4,0,0,5,6,0,0,
     .              0,0,0,0,1,2,3,4, 5,6,7,8,0,0,0,0/

C         TRANSCUT is the transmission to stop the integration at
      TRANSCUT = 5.0E-5
C0.0D0
C         TAUTOL is the maximum optical path for the subgrid intervals
C      it is now a numerical parameter.
C      TAUTOL = 0.2

      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)

C     Set up the arrays for storing information about passed points so that
C     the radiance at each subgrid interval can be accumulated and
C     used to evaluate the gradient.

C     Estimate the maximum number of subgrid intervals that we will pass.
C     MAXCELLSCROSS already restricts the number of cells. The maximum
C     optical path then sets the maximum number of subgrid intervals
C     per cell.
      DZMAX = -1.0
      DO I=1,NZ-1
        IF (ZGRID(I+1) - ZGRID(I) .GE. DZMAX) THEN
          DZMAX = ZGRID(I+1) - ZGRID(I)
        ENDIF
      ENDDO
      MAXTAU = MAXVAL(EXTINCT)*MAX(DELX, DELY, DZMAX)
      MAXBYOPT = INT(MAXCELLSCROSS*MAXTAU/TAUTOL)
      MAXSUBGRIDINTS = INT(MAX(MAXCELLSCROSS,MAXBYOPT)) + 1
      ALLOCATE (PASSEDPOINTS(8,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDINTERP0(8,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDINTERP1(8,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDDELS(MAXSUBGRIDINTS))
      ALLOCATE (PASSEDRAD(NSTOKES,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDABSCELL(MAXSUBGRIDINTS))
      ALLOCATE (PASSEDTRANSMIT(MAXSUBGRIDINTS))

C     initialize counter for subgrid intervals.
      NPASSED = 1
C     initialize radiance at passed subgrid intervals as zero
      PASSEDRAD = 0.0
C     initialize all other passed points to horrific numbers so its obvious if there
C     is an error.
      PASSEDPOINTS = 0
      PASSEDINTERP0 = 99999.0
      PASSEDINTERP1 = 99999.0
      PASSEDDELS = 99999.0
      PASSEDABSCELL = -99999.0
      PASSEDTRANSMIT = 99999.0

      PI = ACOS(-1.0D0)
C       Calculate the generalized spherical harmonics for this direction
      CALL YLMALL (.FALSE.,SNGL(MU2),SNGL(PHI2),ML,MM,NSTLEG, YLMDIR)
      SECMU0 = 1.0D0/ABS(SOLARMU)
      IF (SRCTYPE .NE. 'T') THEN
C       This is modified so that its done even for no delta-M so the
C       'exact single scatter' derivative can still be calculated.
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
          IERR = 1
          WRITE (ERRMSG,*)'GRAD_INTEGRATE_1RAY: ICELL=',ICELL,
     .                MU2,PHI2,XE,YE,ZE
          RETURN
        ENDIF
        NGRID = NGRID + 1

C           Decide which of the eight grid points we need the source function
C           These are also the points we need to recalculate gradient
C           terms for.
        DO I = 1, 8
          DONETHIS(I) = DONEFACE(I,IFACE+1)
          IF (NX .EQ. 1 .AND. ONEX(I) .LT. 0) DONETHIS(I) = ONEX(I)
          IF (NY .EQ. 1 .AND. ONEY(I) .LT. 0) DONETHIS(I) = ONEY(I)
          OEXTINCT8(I) = EXTINCT8(I)
          OSRCEXT8(:,I) = SRCEXT8(:,I)
          OSINGSCAT8(:,I) = SINGSCAT8(:,I)
          OGRAD8(:,:,I,:) = GRAD8(:,:,I,:)
          OBASEADAPTINTERP(:,I) = BASEADAPTINTERP(:,I)
        ENDDO

C         Compute the source function times extinction in direction (MU2,PHI2)
C     Compute the quantities at each base and adaptive grid combination
C     that will be used for calculating sensitivities with respect to the
C     scattering kernel/emission/solar source.
        CALL COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR,
     .            NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .            NPTS, DELTAM, SRCTYPE, SOLARMU,
     .            EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .            SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .            DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,
     .            EXTINCT8, SRCEXT8, TOTAL_EXT, NPART, RSHPTR, RADIANCE,
     .            LOFJ, CELLFLAGS, PARTDER, NUMDER,
     .            DNUMPHASE, DEXT, DALB, DIPHASE, DLEG, MAXPG, BCELL,
     .            DSINGSCAT, SINGSCAT8, OSINGSCAT8, DIPHASEIND,
     .            PLANCK, NODIFFUSE,
     .            GRAD8, OGRAD8, BASEADAPTINTERP,
     .            OBASEADAPTINTERP, GRIDPOS)

C         Interpolate the source and extinction to the current point

C       Note that this recalculation of SRCEXT1 when moving between cells
C       rather than inheritance is the cause of discrepancy in radiance
C       calculation between
C       INTEGRATE_SOURCE (original SHDOM) and INTEGRATE_1RAY/GRAD_INTEGRATE_1RAY
C       which is used in pyshdom and only for the 'CAMERA_MODE' calculations
C       in original SHDOM.
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

C       OUTOFDOMAIN tests whether we are in a cell which has
C       an independent column boundary condition in either X (first bit)
C       or Y. Gradient contributions are not calculated for these components
C       This is now removed as gradients do appear to be calculated
C       without larger errors at boundary points or in 1D mode.
C        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
C     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        OUTOFDOMAIN = .FALSE.
        IF (.NOT. OUTOFDOMAIN) THEN
C         Base grid cell interpolation for radiance gradient contribution.
          CALL GET_INTERP_KERNEL(BCELL,GRIDPTR,GRIDPOS,XE,YE,ZE,FB)
          PASSEDINTERP1(:,NPASSED) = FB(:)
C         Single scatter for EXACT_SINGLE_SCATTER contribution to gradient.
          DO N=1,8
              SINGSCAT1(:,N) = FC(N)*SINGSCAT8(:,N)
              SINGSCAT1(1,N) = MAX(0.0, SINGSCAT1(1,N))
          ENDDO
C         Sum 'SOURCE' gradient contributions from each adaptive grid cell.
C         to get the total for each base grid cell.
          GRAD1 = 0.0
          DO IDR=1,NUMDER
            IPA = PARTDER(IDR)
            DO KK=1,8
              DO N=1,8
                IF (EXTINCT(GRIDPTR(N,ICELL),IPA) .NE. 0.0) THEN
                GRAD1(:,KK,IDR) = GRAD1(:,KK,IDR) + FC(N)*
     .            GRAD8(:,KK,N,IDR)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ELSE
          PASSEDINTERP1(:,NPASSED) = 0.0
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
          IERR = 1
          WRITE (ERRMSG,*) 'GRAD_INTEGRATE_1RAY: SO<0  ',
     .      MU2,PHI2,XE,YE,ZE,SO,ICELL
          RETURN
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

C         save the interpolation kernel and pointers for each subgrid
C         interval so that the radiance gradient can be calculated.
          PASSEDDELS(NPASSED) = DELS
          IF (.NOT. OUTOFDOMAIN) THEN
            PASSEDINTERP1(:,NPASSED) = FB(:)
          ELSE
            PASSEDINTERP1(:,NPASSED) = 0.0
          ENDIF
          PASSEDPOINTS(:,NPASSED) = GRIDPTR(:,BCELL)
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
C           Compute the base grid cell interpolation factors which are
C           used in the radiance part of the gradient.
            CALL GET_INTERP_KERNEL(BCELL,GRIDPTR,GRIDPOS,XI,YI,ZI,FB)
            PASSEDINTERP0(:,NPASSED) = FB(:)
C           Single scatter is used for EXACT_SINGLE_SCATTER portions of gradient.
            DO N=1,8
                SINGSCAT0(:,N) = FC(N)*SINGSCAT8(:,N)
                SINGSCAT0(1,N) = MAX(0.0, SINGSCAT0(1,N))
            ENDDO
C         Sum over the adaptive points to give the gradient contribution
C         to each base point from the 'SOURCE'.
            GRAD0= 0.0
            DO IDR=1,NUMDER
              IPA = PARTDER(IDR)
              DO KK=1,8
                DO N=1,8
                IF (EXTINCT(GRIDPTR(N,ICELL),IPA) .NE. 0.0) THEN
                  GRAD0(:,KK,IDR) = GRAD0(:,KK,IDR) + FC(N)*
     .            GRAD8(:,KK,N,IDR)
                ENDIF
                ENDDO
              ENDDO
            ENDDO
          ELSE
            PASSEDINTERP0(:,NPASSED) = 0.0
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
C             Inner product between scattering kernel/emission/singlescatter (SOURCE)
C             related gradient terms and the 'adjoint direct beam' which is
C             just a weighting by transmission to sensor.
C             The neglect of other terms in the integral here is the Levis
C             approximation to the Frechet derivative of SHDOM.
C             The portion of the inner product from integration of radiance
C             multiplied by extinction is calculated at the end of the
C             subroutine.
              SRCGRAD = ( 0.5*(GRAD0+GRAD1)
     .          + 0.08333333333*(EXT0*GRAD1-EXT1*GRAD0)*DELS
     .           *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
               IF (EXACT_SINGLE_SCATTER .AND.
     .           SRCTYPE .NE. 'T') THEN
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
C         Add the radiance increment to the gridpoints we have passed.
C         so that component of the gradient can be calculated.
C         (This is done finally at the end of this subroutine).
          PASSEDABSCELL(NPASSED) = ABSCELL
          PASSEDTRANSMIT(NPASSED) = TRANSMIT
          DO KK=1,NPASSED
            PASSEDRAD(:,KK) = PASSEDRAD(:,KK) +
     .        TRANSMIT*SRC(:)*ABSCELL/PASSEDTRANSMIT(KK)
          ENDDO

          IF (.NOT. OUTOFDOMAIN) THEN
            DO KK = 1, 8
              GRIDPOINT = GRIDPTR(KK,BCELL)
              RAYGRAD(:,GRIDPOINT,:) = RAYGRAD(:,GRIDPOINT,:) +
     .           TRANSMIT*SRCGRAD(:,KK,:)*ABSCELL

C             Add gradient component due to the direct solar beam.
C             Sensitivity of single scattered radiation to
C             extinction along the path between the gridpoint and the sun.
              IF (EXACT_SINGLE_SCATTER .AND.
     .          SRCTYPE .NE. 'T') THEN
                II = 1
                GRIDPOINT = GRIDPTR(KK,BCELL)
                DO IDR=1,NUMDER
                  DO WHILE (DPTR(II,GRIDPOINT).GT.0)
                    SSP = DPTR(II,GRIDPOINT)
                    RAYGRAD(:,SSP,IDR) = RAYGRAD(:,SSP,IDR) -
     .                 DPATH(II,GRIDPOINT)*DEXT(SSP,IDR)*
     .                 ABSCELL*TRANSMIT*
     .                 SRCSINGSCAT(:,KK)*DIRFLUX(SSP)*SECMU0
                    II = II + 1
                    ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF
          NPASSED = NPASSED + 1
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
C
C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
        IF (TRANSMIT .LT. TRANSCUT .OR. NGRID.GT.MAXCELLSCROSS) THEN
          VALIDRAD = .TRUE.

        ELSE IF (INEXTCELL .EQ. 0 .AND. IFACE .GE. 5) THEN
          VALIDRAD = .TRUE.
          CALL FIND_BOUNDARY_RADIANCE_GRAD (NSTOKES, XN, YN,
     .                      SNGL(MU2), SNGL(PHI2),
     .                      IC, KFACE, GRIDPTR, GRIDPOS,
     .                      MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .                      NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .                      SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, DIRFLUX,
     .                      SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                      RADBND, BOUNDPTS, BOUNDINTERP, DIRRAD,
     .                      GNDALBEDO)
          RADOUT(:) = RADOUT(:) + TRANSMIT*RADBND(:)
          PASSEDTRANSMIT(NPASSED) = TRANSMIT
          PASSEDABSCELL(NPASSED) = 0.0
          DO KK=1,NPASSED
            PASSEDRAD(:,KK) = PASSEDRAD(:,KK) +
     .        TRANSMIT*RADBND(:)/PASSEDTRANSMIT(KK)
          ENDDO
C         Calculate the sensitivity of the surface reflection to
C         extinction along the path between the surface and the sun.
C         This does not include sensitivity of surface reflection to
C         downwelling diffuse radiation, only the reflection of the direct
C         solar beam.
          IF (EXACT_SINGLE_SCATTER .AND.
     .          SRCTYPE .NE. 'T') THEN
            DO KK=1,4
              II = 1
              GRIDPOINT = BOUNDPTS(KK)
              DO IDR=1,NUMDER
                DO WHILE (DPTR(II,GRIDPOINT).GT.0)
                  SSP = DPTR(II,GRIDPOINT)
                  RAYGRAD(:,SSP,IDR) = RAYGRAD(:,SSP,IDR) -
     .                 DPATH(II,GRIDPOINT)*DEXT(SSP,IDR)*
     .                 TRANSMIT*DIRRAD(:,KK)*BOUNDINTERP(KK)
                  II = II + 1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
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
C     Add in the gradient components due to the radiance that is now
C     fully calculated using the saved properties from each
C     subgrid integration interval.
      DO IDR=1,NUMDER
        IPA = PARTDER(IDR)
        DO KK=1,NPASSED-1
          EXT0 = 0.0
          EXT1 = 0.0
          DELS = PASSEDDELS(KK)
          DO K=1,8
            EXT0 = EXT0 + TOTAL_EXT(PASSEDPOINTS(K,KK))*
     .        PASSEDINTERP0(K,KK)
            EXT1 = EXT1 + TOTAL_EXT(PASSEDPOINTS(K,KK))*
     .        PASSEDINTERP1(K,KK)
          ENDDO
          EXT = 0.5*(EXT0+EXT1)
          IF (EXT .NE. 0.0) THEN
            DO K=1,8
              RAD0(:) = -1*PASSEDRAD(:,KK+1)*
     .        DEXT(PASSEDPOINTS(K,KK),IDR)*
     .        PASSEDINTERP0(K,KK)
              RAD1(:) = -1*PASSEDRAD(:,KK)*
     .        DEXT(PASSEDPOINTS(K,KK),IDR)*
     .        PASSEDINTERP1(K,KK)
              RADGRAD(:) = ( 0.5*(RAD0+RAD1)
     .           + 0.08333333333*(EXT0*RAD1-EXT1*RAD0)*DELS
     .                *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
              GRIDPOINT = PASSEDPOINTS(K,KK)

              RAYGRAD(:,GRIDPOINT,IDR) = RAYGRAD(:,GRIDPOINT,IDR)+
     .          RADGRAD(:)*PASSEDTRANSMIT(KK)*PASSEDABSCELL(KK)
            ENDDO
          ENDIF

        ENDDO
      ENDDO
      DEALLOCATE (PASSEDPOINTS, PASSEDRAD,PASSEDINTERP0,
     .            PASSEDINTERP1, PASSEDDELS, PASSEDABSCELL,
     .            PASSEDTRANSMIT)
      RETURN
      END

      SUBROUTINE COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR,
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,
     .             EXTINCT8, SRCEXT8, TOTAL_EXT, NPART, RSHPTR,
     .             RADIANCE, LOFJ, CELLFLAGS, PARTDER,
     .             NUMDER, DNUMPHASE, DEXT, DALB, DIPHASE, DLEG, MAXPG,
     .             BCELL, DSINGSCAT, SINGSCAT8, OSINGSCAT8,DIPHASEIND,
     .             PLANCK, NODIFFUSE,
     .             GRAD8,OGRAD8, BASEADAPTINTERP,
     .             OBASEADAPTINTERP, GRIDPOS)
C       Computes the source function times extinction for gridpoints
C     belonging to cell ICELL in the direction (MU,PHI).  The results
C     are returned in SRCEXT8 and EXTINCT8.
C     The spherical harmonic source function series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
C     Evaluate the 'source of the linearized RTE' for each combination of
C     of base and adaptive grid point.
C     This is unapproximated (apart from practicalities of discretization).
      IMPLICIT NONE
      INTEGER ICELL, NSTOKES, NSTLEG, NPTS, ML, MM, NLM, NLEG, NUMPHASE
      INTEGER GRIDPTR(8,*), SHPTR(*), RSHPTR(*)
      INTEGER DONETHIS(8), OLDIPTS(8), BCELL, DNUMPHASE
      INTEGER IPHASE(NPTS,NPART), NPART, NUMDER, MAXPG
      INTEGER DIPHASEIND(NPTS,NUMDER)
      LOGICAL DELTAM
      REAL    SOLARMU
      REAL    GRIDPOS(3,*)
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(NSTLEG,0:NLEG,*), TOTAL_EXT(*)
      REAL    PLANCK(NPTS,NPART)
      REAL    DIRFLUX(*), SOURCE(NSTOKES,*)
      REAL    RADIANCE(NSTOKES,*)
      REAL    LONGRADIANCE(NSTOKES, NPTS)
      CHARACTER USELONGRAD
      REAL    YLMDIR(NSTLEG,*), YLMSUN(NSTLEG,*)
      REAL    SINGSCAT(NSTOKES,NUMPHASE)
      REAL    OEXTINCT8(8), OSRCEXT8(NSTOKES,8)
      REAL    EXTINCT8(8), SRCEXT8(NSTOKES,8), W
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
      CHARACTER SRCTYPE*1
      INTEGER*2 CELLFLAGS(*)
      LOGICAL OUTOFDOMAIN, NODIFFUSE
      REAL    SINGSCAT8(NSTOKES,8), OSINGSCAT8(NSTOKES,8)
      REAL    DLEG(NSTLEG,0:NLEG,DNUMPHASE), DEXT(MAXPG,NUMDER)
      REAL    DALB(MAXPG,NUMDER)
      INTEGER DIPHASE(MAXPG,NUMDER), KK
      REAL    DSINGSCAT(NSTOKES,DNUMPHASE)
      INTEGER LOFJ(*), PARTDER(NUMDER), RNS, RIS, IDR, IB,NB
      INTEGER IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I,IPA,K
      REAL    SECMU0, F, DA, A, A1, B1, EXT
      INTEGER IDP
      REAL GRAD8(NSTOKES,8,8,NUMDER), SOURCET(NSTOKES)
      REAL OGRAD8(NSTOKES,8,8,NUMDER)
      REAL PHAGRAD(NSTOKES)
      DOUBLE PRECISION BASEADAPTINTERP(8,8)
      DOUBLE PRECISION OBASEADAPTINTERP(8,8)
      DOUBLE PRECISION X0,Y0,Z0
      SECMU0 = 1.0D0/ABS(SOLARMU)

C         Loop over the grid points, computing the source function
C           at the viewing angle from the spherical harmonics source function.
C     Also calculate the contribution to the gradient at each pair of adaptive
C     and base grid points.
      DO N = 1, 8
        IP = GRIDPTR(N,ICELL)
        I = DONETHIS(N)
C     We reuse calculations from gridpoints that were part of the previous
C     cell.
        IF (I .GT. 0 .AND. IP .EQ. OLDIPTS(N)) THEN
          EXTINCT8(N) = OEXTINCT8(I)
          SRCEXT8(:,N) = OSRCEXT8(:,I)
          SINGSCAT8(:,N) = OSINGSCAT8(:,I)
          GRAD8(:,:,N,:) = OGRAD8(:,:,I,:)
          BASEADAPTINTERP(:,N) = OBASEADAPTINTERP(:,I)
        ELSE IF (I .LT. 0) THEN
          EXTINCT8(N) = EXTINCT8(ABS(I))
          SRCEXT8(:,N) = SRCEXT8(:,ABS(I))
          SINGSCAT8(:,N) = SINGSCAT8(:,ABS(I))
          GRAD8(:,:,N,:) = GRAD8(:,:,ABS(I),:)
          BASEADAPTINTERP(:,N) = BASEADAPTINTERP(:,ABS(I))
        ELSE
C         Calculate new points. First the gradient.
          EXT = TOTAL_EXT(IP)
          OLDIPTS(N) = IP
          IS = SHPTR(IP)
          NS = SHPTR(IP+1)-IS
          GRAD8(:,:,N,:) = 0.0

C         See comment in above subroutine (GRAD_INTEGRATE_1RAY)
C         for removal of this flag.
C          OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
C     .                   BTEST(INT(CELLFLAGS(ICELL)),1))
          OUTOFDOMAIN = .FALSE.
          BASEADAPTINTERP(:,N) = 0.0D0

          IF (.NOT.OUTOFDOMAIN) THEN
C           do interpolation from base to adaptive grid points
C           if we are in a non base grid cell.
            IF (BCELL .NE. ICELL) THEN
              X0 = GRIDPOS(1,IP)
              Y0 = GRIDPOS(2,IP)
              Z0 = GRIDPOS(3,IP)
              CALL GET_INTERP_KERNEL(BCELL,GRIDPTR,GRIDPOS,
     .          X0,Y0,Z0,BASEADAPTINTERP(:,N))
            ENDIF
C          overwrite any interpolation factors for grid points
C          that are also base grid points just to doubly  make sure
C          that they are 1.0 and rest are zero.
C          This is also the default behaviour when we are in a
C          base grid cell.
            DO NB=1,8
             IB = GRIDPTR(NB,BCELL)
             IF (IB .EQ. IP) THEN
               BASEADAPTINTERP(:,N) = 0.0
               BASEADAPTINTERP(N,N) = 1.0
             ENDIF
            ENDDO

            RIS = RSHPTR(IP)
            RNS = RSHPTR(IP+1)-RIS
            DO IDR = 1, NUMDER

              IPA = PARTDER(IDR)
              IF (NUMPHASE .GT. 0) THEN
                K = IPHASE(IP,IPA)
              ELSE
                K = IP
              ENDIF
C             Compute the radiance/phase convolutions for each derivative
C             species and adaptive point which is used for the extinction
C             and albedo derivatives. This doesn't need to be redone
C             for each base grid point as well.
              SOURCET = 0.0
C             For diagnostic purposes we might want to neglect the diffuse
C             contribution, for example, if wanting to confirm that the
C             single scatter derivative is exact. Normally NODIFFUSE = .FALSE..
              IF (.NOT. NODIFFUSE) THEN
                DO J = 1, RNS
                  SOURCET(1) = SOURCET(1) + RADIANCE(1,RIS+J)*
     .              YLMDIR(1,J)*LEGEN(1,LOFJ(J),K)
                ENDDO
                IF (NSTOKES .GT. 1) THEN
                  DO J = 1, RNS
                    SOURCET(1) = SOURCET(1) + LEGEN(5,LOFJ(J),K)
     .              *RADIANCE(2,RIS+J)*YLMDIR(1,J)
                  ENDDO
                  DO J = 5, RNS
                    SOURCET(2) = SOURCET(2) +LEGEN(5,LOFJ(J),K)*
     .                       RADIANCE(1,J)*YLMDIR(2,J)
     .                 + LEGEN(2,LOFJ(J),K)*RADIANCE(2,J)*YLMDIR(2,J)
     .                 + LEGEN(3,LOFJ(J),K)*RADIANCE(3,RIS+J)*
     .                       YLMDIR(5,J)
                    SOURCET(3) = SOURCET(3) + LEGEN(5,LOFJ(J),K)*
     .                      RADIANCE(1,RIS+J)*YLMDIR(6,J)
     .                      + LEGEN(2,LOFJ(J),K)*RADIANCE(2,RIS+J)
     .                        *YLMDIR(6,J)
     .                      + LEGEN(3,LOFJ(J),K)*RADIANCE(3,RIS+J)*
     .                         YLMDIR(3,J)
                  ENDDO
                ENDIF
                IF (NSTOKES .EQ. 4) THEN
                  DO J = 1, RNS
                    SOURCET(2) = SOURCET(2) + LEGEN(6,LOFJ(J),K)*
     .                  RADIANCE(4,RIS+J)*YLMDIR(5,J)
                    SOURCET(3) =SOURCET(3) + LEGEN(6,LOFJ(J),K)*
     .                  RADIANCE(4,RIS+J)*YLMDIR(3,J)
                    SOURCET(4) =SOURCET(4) - LEGEN(6,LOFJ(J),K)*
     .                    RADIANCE(3,RIS+J)*YLMDIR(4,J)
     .                +LEGEN(4,LOFJ(J),K)*RADIANCE(4,RIS+J)*
     .                    YLMDIR(4,J)
                  ENDDO
                ENDIF
              ENDIF
C             Add the single scatter contribution to the source.
              IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
                IF (NUMPHASE .GT. 0) THEN
                  SOURCET(:) = SOURCET(:) + DIRFLUX(IP)*SECMU0
     .              *SINGSCAT(:,K)
                ELSE
                  WRITE(*,*) 'NUMPHASE=', NUMPHASE, ' NOT SUPPORTED'
                  STOP
                ENDIF
              ELSEIF (SRCTYPE .NE. 'T') THEN
                DA = ALBEDO(IP,IPA)*DIRFLUX(IP)*SECMU0*W
                J = 1
                DO L = 0, ML
                  ME = MIN(L,MM)
                  MS = -ME
                  A1 = DA*LEGEN(1,L,K)
                  B1 = DA*LEGEN(5,L,K)
                  IF (J .LE. NS) THEN
                    JT = J
                    DO M = MS, ME
                      SOURCET(1) =SOURCET(1) +
     .                  A1*YLMDIR(1,J)*YLMSUN(1,J)
                      J = J + 1
                    ENDDO
                    IF (NSTOKES .GT. 1) THEN
                      J = JT
                      DO M = MS, ME
                        SOURCET(2)=SOURCET(2) +
     .                    B1*YLMDIR(2,J)*YLMSUN(1,J)
                        SOURCET(3)=SOURCET(3) +
     .                    B1*YLMDIR(6,J)*YLMSUN(1,J)
                        J = J + 1
                      ENDDO
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF

              SOURCET(1) = MAX(0.0, SOURCET(1))
              DO NB=1,8
                IB = GRIDPTR(NB,BCELL)
C               Calculate the derivatives of extinction and albedo for each
C               combination of base and adaptive grid point and scatter species.
                GRAD8(:,NB,N,IDR) = BASEADAPTINTERP(NB,N)*SOURCET(:)*(
     .            ALBEDO(IB,IPA)*DEXT(IB,IDR) +
     .           EXTINCT(IB,IPA)*DALB(IB,IDR))

C              Add the contributions for thermal emission to the
C              extinction and albedo derivatives.
                IF (SRCTYPE .NE. 'S' .AND. (1.0 - ALBEDO(IP,IPA) .GE.
     .              1e-9) ) THEN
C                 Make sure no divide by zero error. We could also
C                 recalculate Planck radiance here from temperature.
                  GRAD8(1,NB,N,IDR) = GRAD8(1,NB,N,IDR) +
     .                (BASEADAPTINTERP(NB,N)*PLANCK(IP,IPA)
     .                /(1.0 - ALBEDO(IP,IPA))) *
     .                (DEXT(IB,IDR)*(1.0 - ALBEDO(IB,IPA)) -
     .                EXTINCT(IB,IDR)*DALB(IB,IDR))
                ENDIF
C               Check if we need a phase derivative for this pair of adaptive and
C               base grid cell. The interpolation rule for phase functions is non-smooth;
C               each adaptive point simply inherits the phase function of the
C               base point with the maximum scattering coefficient from the 8
C               possible candidates.
C               If the base point is the inheriting point then we add it.
C               The equivalent to BASEADAPTINTERP is just 1.0 for the inheriting
C               point and zero elsewhere.
                PHAGRAD = 0.0
                IF (DIPHASEIND(IP,IDR) .EQ. IB) THEN
                  IDP = DIPHASE(IB,IDR)
                  IF (.NOT. NODIFFUSE) THEN
                    DO J = 1, RNS
                      GRAD8(1,NB,N,IDR) = GRAD8(1,NB,N,IDR) +
     .                 RADIANCE(1,RIS+J)*YLMDIR(1,J)*DLEG(1,LOFJ(J),IDP)
                    ENDDO
                    IF (NSTOKES .GT. 1) THEN
                      DO J = 1, RNS
                        GRAD8(1,NB,N,IDR) =GRAD8(1,NB,N,IDR)
     .                     +DLEG(5,LOFJ(J),IDP)*
     .                       RADIANCE(2,RIS+J)*YLMDIR(1,J)
                      ENDDO
                      DO J = 5, RNS
                        GRAD8(2,NB,N,IDR) = GRAD8(2,NB,N,IDR)
     .                     +DLEG(5,LOFJ(J),IDP)*
     .                       RADIANCE(1,J)*YLMDIR(2,J)
     .                 + DLEG(2,LOFJ(J),IDP)*RADIANCE(2,J)*YLMDIR(2,J)
     .                 + DLEG(3,LOFJ(J),IDP)*RADIANCE(3,RIS+J)*
     .                       YLMDIR(5,J)
                        GRAD8(3,NB,N,IDR)=GRAD8(3,NB,N,IDR) +
     .                        DLEG(5,LOFJ(J),IDP)*
     .                      RADIANCE(1,RIS+J)*YLMDIR(6,J)
     .                      + DLEG(2,LOFJ(J),IDP)*RADIANCE(2,RIS+J)
     .                        *YLMDIR(6,J)
     .                      + DLEG(3,LOFJ(J),IDP)*RADIANCE(3,RIS+J)*
     .                         YLMDIR(3,J)
                      ENDDO
                    ENDIF
                    IF (NSTOKES .EQ. 4) THEN
                      DO J = 1, RNS
                        GRAD8(2,NB,N,IDR) =GRAD8(2,NB,N,IDR)
     .                     +DLEG(6,LOFJ(J),IDP)*
     .                    RADIANCE(4,RIS+J)*YLMDIR(5,J)
                        GRAD8(3,NB,N,IDR) =GRAD8(3,NB,N,IDR) +
     .                      DLEG(6,LOFJ(J),IDP)*
     .                        RADIANCE(4,RIS+J)*YLMDIR(3,J)
                        GRAD8(4,NB,N,IDR) =GRAD8(4,NB,N,IDR) -
     .                      DLEG(6,LOFJ(J),IDP)*
     .                        RADIANCE(3,RIS+J)*YLMDIR(4,J)
     .                + DLEG(4,LOFJ(J),IDP)*RADIANCE(4,RIS+J)*
     .                    YLMDIR(4,J)
                      ENDDO
                    ENDIF
                  ENDIF
                  GRAD8(:,NB,N,IDR) = GRAD8(:,NB,N,IDR) +
     .              PHAGRAD*ALBEDO(IP,IPA)*EXTINCT(IP,IPA)
C             Add the single scatter contribution to the phase derivative.
                  IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
                    IF (NUMPHASE .GT. 0) THEN
                      GRAD8(:,NB,N,IDR) = GRAD8(:,NB,N,IDR) +
     .                  DIRFLUX(IP)*SECMU0*EXTINCT(IP,IPA)*
     .                  ALBEDO(IP,IPA)*DSINGSCAT(:,IDP)
                    ELSE
                      WRITE(*,*) 'NUMPHASE=', NUMPHASE, ' NOT SUPPORTED'
                      STOP
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF
C         Now we perform the calculations for the forward model, calculating thw
C         SRCEXT product. This involves some recalculations so its not that
C         efficient. This could be optimized further.

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
              IF (SRCTYPE .NE. 'T') THEN
                DA = ALBEDO(IP,IPA)*DIRFLUX(IP)*SECMU0*W
                SINGSCAT8(:,N) = SINGSCAT8(:,N) + DA*SINGSCAT(:,IPH)
                IF (DELTAM) THEN
                  SRCEXT8(:,N) = SRCEXT8(:,N) + DA*SINGSCAT(:,IPH)
                ENDIF
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
C         The SRCEXT8 is used to calculate the radiance, which is itself
C         used in the gradient. If we only want the single scatter (NODIFFUSE) then
C         the source should only include the single scatter so we set that
C         here. This doesn't save computation time at all though.
          IF (NODIFFUSE) THEN
            SRCEXT8(:,N) = SINGSCAT8(:,N)*EXT
          ELSE
            SRCEXT8(:,N) = SRCEXT8(:,N)*EXT
          ENDIF
          EXTINCT8(N) = EXT
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE COMPUTE_SOURCE_GRAD_1CELL_OLD (ICELL, GRIDPTR,
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,
     .             EXTINCT8, SRCEXT8, TOTAL_EXT, NPART, RSHPTR,
     .             RADIANCE, OGRAD8, GRAD8, LOFJ, CELLFLAGS, PARTDER,
     .             NUMDER, DNUMPHASE, DEXT, DALB, DIPHASE, DLEG, NBPTS,
     .             BCELL, DSINGSCAT, SINGSCAT8, OSINGSCAT8,DIPHASEIND,
     .             PLANCK, LONGRADIANCE, USELONGRAD, NODIFFUSE)
C       Computes the source function times extinction for gridpoints
C     belonging to cell ICELL in the direction (MU,PHI).  The results
C     are returned in SRCEXT8 and EXTINCT8.
C     The spherical harmonic source function series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
C     Evaluate the 'source of the linearized RTE' at each of the grid points.
C     This is unapproximated (apart from practicalities of discretization).
      IMPLICIT NONE
      INTEGER ICELL, NSTOKES, NSTLEG, NPTS, ML, MM, NLM, NLEG, NUMPHASE
      INTEGER GRIDPTR(8,*), SHPTR(*), RSHPTR(*)
      INTEGER DONETHIS(8), OLDIPTS(8), BCELL, DNUMPHASE
      INTEGER IPHASE(NPTS,NPART), NPART, NUMDER, NBPTS
      INTEGER DIPHASEIND(NPTS,NUMDER)
      LOGICAL DELTAM
      REAL    SOLARMU
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(NSTLEG,0:NLEG,*), TOTAL_EXT(*)
      REAL    PLANCK(NPTS,NPART)
      REAL    DIRFLUX(*), SOURCE(NSTOKES,*)
      REAL    RADIANCE(NSTOKES,*)
      REAL    LONGRADIANCE(NSTOKES, NPTS)
      CHARACTER USELONGRAD
      REAL    YLMDIR(NSTLEG,*), YLMSUN(NSTLEG,*)
      REAL    SINGSCAT(NSTOKES,NUMPHASE)
      REAL    OEXTINCT8(8), OSRCEXT8(NSTOKES,8)
      REAL    EXTINCT8(8), SRCEXT8(NSTOKES,8), W
      DOUBLE PRECISION SUNDIRLEG(0:NLEG), UNSCALED_ALBEDO
      CHARACTER SRCTYPE*1
      INTEGER*2 CELLFLAGS(*)
      LOGICAL OUTOFDOMAIN, NODIFFUSE
      REAL    GRAD8(NSTOKES,8,8,NUMDER), OGRAD8(NSTOKES,8,8,NUMDER)
      REAL    SINGSCAT8(NSTOKES,8), OSINGSCAT8(NSTOKES,8)
      REAL    FULL_SINGSCAT(NSTOKES)
      REAL    DLEG(NSTLEG,0:NLEG,DNUMPHASE), DEXT(NBPTS,NUMDER)
      REAL    DALB(NBPTS,NUMDER), DEXTM, DALBM, DLEGM(NSTLEG)
      INTEGER DIPHASE(NBPTS,NUMDER), KK
      REAL    DSINGSCAT(NSTOKES,DNUMPHASE)
      INTEGER LOFJ(*), PARTDER(NUMDER), RNS, RIS, IDR, IB,NB
      INTEGER IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I,IPA,K
      REAL    SECMU0, F, DA, A, A1, B1, EXT
      REAL    C
      PARAMETER (C=3.544907703)

      GRAD8 = 0.0
      SECMU0 = 1.0D0/ABS(SOLARMU)

C         Loop over the grid points, computing the source function
C           at the viewing angle from the spherical harmonics source function.
C     Also calculate the contribution to the gradient at each pair of adaptive
C     and base grid points.
      DO N = 1, 8

        IP = GRIDPTR(N,ICELL)
        I = DONETHIS(N)
        IF (I .GT. 0 .AND. IP .EQ. OLDIPTS(N)) THEN
          EXTINCT8(N) = OEXTINCT8(I)
          SRCEXT8(:,N) = OSRCEXT8(:,I)
          SINGSCAT8(:,N) = OSINGSCAT8(:,I)
          GRAD8(:,:,N,:) = OGRAD8(:,:,I,:)
        ELSE IF (I .LT. 0) THEN
          EXTINCT8(N) = EXTINCT8(ABS(I))
          SRCEXT8(:,N) = SRCEXT8(:,ABS(I))
          SINGSCAT8(:,N) = SINGSCAT8(:,ABS(I))
C         I<0 occurs when nx/ny=1 and because GRAD8
C         is zeroed at the start of this subroutine old values will not
C         be passed correctly and the gradients will likely be incorrect.
C         -JRLoveridge 2021/02/10
          GRAD8(:,:,N,:) = GRAD8(:,:,ABS(I),:)
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

              IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
                DA = DIRFLUX(IP)*SECMU0
                F = LEGEN(1,ML+1,K)
              ELSE
                F=0.0
              ENDIF

              DO NB=1,8
                FULL_SINGSCAT = 0.0
                IB = GRIDPTR(NB,BCELL)
                DEXTM = DEXT(IB,IDR)*EXTINCT(IP,IPA)
                DALBM = DALB(IB,IDR)*EXTINCT(IB,IPA)

C             Sum over the real generalized spherical harmonic series
C             and calculate the source for the linearized
C             RTE (not the Source function). The contribution from each base grid point and
C             (adpative) grid point pair is calculated. - JRLoveridge (2021/02/10)
                IF (.NOT. NODIFFUSE) THEN
                DO J = 1, RNS
                  L = LOFJ(J)
C               The phase derivative is only non-zero for the
C               base grid point which supplied the phase function to
C               the adaptive point. This is due to the non-smooth
C               'max-scatter' interpolation rule currently used for
C               phase functions (TRILIN_INTERP_PROP).
C               To avoid this non-smoothness the interpolation rule
C               should be updated. - JRLoveridge 2021/02/09
                  IF (DIPHASEIND(IP,IDR) .EQ. IB) THEN
                    DLEGM(:) = DLEG(:,L,DIPHASE(IB,IDR))
                  ELSE
                    DLEGM = 0.0
                  ENDIF

                  GRAD8(1,NB,N,IDR) = GRAD8(1,NB,N,IDR) +
     .            RADIANCE(1,RIS+J)*YLMDIR(1,J)*(
     .              DEXTM*(ALBEDO(IP,IPA)*LEGEN(1,L,K)) +
     .              EXTINCT(IP,IPA)*DALBM*LEGEN(1,L,K) +
     .              EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(1))

                  IF (NSTOKES .GT. 1) THEN
                    GRAD8(1,NB,N,IDR) = GRAD8(1,NB,N,IDR) +
     .              RADIANCE(2,RIS+J)*YLMDIR(1,J)*(
     .              DEXTM*ALBEDO(IP,IPA)*LEGEN(5,L,K) +
     .              EXTINCT(IP,IPA)*DALBM*LEGEN(5,L,K) +
     .              EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(5))

                    IF (J .GE. 5) THEN
                      GRAD8(2,NB,N,IDR) = GRAD8(2,NB,N,IDR) +
     .                RADIANCE(1,RIS+J)*YLMDIR(2,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(5,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(5,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(5))
                      GRAD8(2,NB,N,IDR) = GRAD8(2,NB,N,IDR) +
     .                RADIANCE(2,RIS+J)*YLMDIR(2,J)*(
     .                DEXTM*(ALBEDO(IP,IPA)*LEGEN(2,L,K)) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(2,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(2))
                      GRAD8(2,NB,N,IDR) = GRAD8(2,NB,N,IDR) +
     .                RADIANCE(3,RIS+J)*YLMDIR(5,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(3,L,K)  +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(3,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(3))

                      GRAD8(3,NB,N,IDR) = GRAD8(3,NB,N,IDR) +
     .                RADIANCE(1,RIS+J)*YLMDIR(6,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(5,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(5,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(5))
                      GRAD8(3,NB,N,IDR) = GRAD8(3,NB,N,IDR) +
     .                RADIANCE(2,RIS+J)*YLMDIR(6,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(2,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(2,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(2))
                      GRAD8(3,NB,N,IDR) = GRAD8(3,NB,N,IDR) +
     .                RADIANCE(3,RIS+J)*YLMDIR(3,J)*(
     .                DEXTM*(ALBEDO(IP,IPA)*LEGEN(3,L,K)) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(3,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(3))
                    ENDIF
                  ENDIF
                  IF (NSTOKES .EQ. 4) THEN
                      GRAD8(2,NB,N,IDR) = GRAD8(2,NB,N,IDR) +
     .                RADIANCE(4,RIS+J)*YLMDIR(5,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(6))

                      GRAD8(3,NB,N,IDR) = GRAD8(3,NB,N,IDR) +
     .                RADIANCE(4,RIS+J)*YLMDIR(3,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(6))

                      GRAD8(4,NB,N,IDR) = GRAD8(4,NB,N,IDR) -
     .                RADIANCE(3,RIS+J)*YLMDIR(4,J)*(
     .                DEXTM*ALBEDO(IP,IPA)*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(6,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(6))

                      GRAD8(4,NB,N,IDR) = GRAD8(4,NB,N,IDR) +
     .                RADIANCE(4,RIS+J)*YLMDIR(4,J)*(
     .                DEXTM*(ALBEDO(IP,IPA)*LEGEN(4,L,K)) +
     .                EXTINCT(IP,IPA)*DALBM*LEGEN(4,L,K) +
     .                EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*DLEGM(4))
                  ENDIF
                ENDDO
                ENDIF
C             Add in the single scattering contribution for the original unscaled phase function.
C             For DELTA M and L<=ML this requires reconstructing the original phase function values
C              Note that this is part of the SOURCE derivative.
                IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
                  IF (NUMPHASE .GT. 0) THEN
C                   The phase function derivative is only non-zero for the
C                   base grid point which supplied the phase function to the (adaptive)
C                   grid point.
                    FULL_SINGSCAT(:) = DA*SINGSCAT(:,K)*(
     .                EXTINCT(IP,IPA)*DALBM + DEXTM*
     .                ALBEDO(IP,IPA))
                    IF (DIPHASEIND(IP,IDR) .EQ. IB) THEN
                      FULL_SINGSCAT(:) = FULL_SINGSCAT(:) +
     .                DA*EXTINCT(IP,IPA)*ALBEDO(IP,IPA)*
     .                DSINGSCAT(:,DIPHASE(IB,IDR))
                    ENDIF
                  ELSE
                    WRITE(*,*) 'NUMPHASE=', NUMPHASE, ' NOT SUPPORTED'
                    STOP
                  ENDIF
                ENDIF
                IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
                  GRAD8(:,NB,N,IDR) = GRAD8(:,NB,N,IDR)
     .                + FULL_SINGSCAT(:)
                ENDIF
                IF (USELONGRAD .EQ. 'T') THEN
                  GRAD8(:,NB,N,IDR) = GRAD8(:,NB,N,IDR) -
     .              DEXTM*LONGRADIANCE(:,IP)
                ENDIF
                IF (USELONGRAD .EQ. 'F') THEN
                  DO J=1,RNS
                    GRAD8(:,NB,N,IDR) = GRAD8(:,NB,N,IDR) -
     .                DEXTM*RADIANCE(:,RIS+J)*YLMDIR(:,J)
                  ENDDO
                ENDIF
                IF (SRCTYPE .NE. 'S') THEN
                  GRAD8(:,NB,N,IDR) = GRAD8(:,NB,N,IDR) +
     .              DEXTM*PLANCK(IP, IPA)
                  IF (1.0 - ALBEDO(IP,IPA) .GE. 1e-9) THEN
                    GRAD8(:,NB,N,IDR) = GRAD8(:,NB,N,IDR) -
     .              EXTINCT(IP,IPA)*DALBM*PLANCK(IP, IPA)/
     .              (1.0 - ALBEDO(IP, IPA))
                  ENDIF
                ENDIF
              ENDDO

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
              IF (SRCTYPE .NE. 'T') THEN
                DA = ALBEDO(IP,IPA)*DIRFLUX(IP)*SECMU0*W
                SINGSCAT8(:,N) = SINGSCAT8(:,N) + DA*SINGSCAT(:,IPH)
                IF (DELTAM) THEN
                  SRCEXT8(:,N) = SRCEXT8(:,N) + DA*SINGSCAT(:,IPH)
                ENDIF
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
          IF (NODIFFUSE) THEN
            SRCEXT8(:,N) = SINGSCAT8(:,N)*EXT
          ELSE
            SRCEXT8(:,N) = SRCEXT8(:,N)*EXT
          ENDIF
          EXTINCT8(N) = EXT
        ENDIF
      ENDDO

      RETURN
      END


C     This subroutine has not been updated with the correct
C     microphysical derivative interpolation. It is no longer supported.
C     -JRLoveridge 2021/02/10
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


      SUBROUTINE FIND_BOUNDARY_RADIANCE_GRAD (NSTOKES, XB, YB, MU2,
     .                      PHI2,
     .                      ICELL, KFACE, GRIDPTR, GRIDPOS,
     .                      MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .                      NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .                      SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, DIRFLUX,
     .                      SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                      RADBND, BOUNDPTS, BOUNDINTERP,DIRRAD,
     .                      GNDALBEDO)
C       Returns the interpolated Stokes radiance at the boundary (RADBND).
C     Inputs are the boundary location (XB,YB), ray direction away from
C     boundary (MU2,PHI2), cell number (ICELL) and face (KFACE) at
C     the boundary point.
      IMPLICIT NONE
      INTEGER NSTOKES, ICELL, KFACE, MAXNBC, NTOPPTS, NBOTPTS
      INTEGER GRIDPTR(8,*), BCPTR(MAXNBC,2)
      INTEGER NMU, NPHI0MAX, NPHI0(*), NSFCPAR
      REAL    MU2, PHI2
      INTEGER BOUNDPTS(4)
      DOUBLE PRECISION BOUNDINTERP(4), DIRRAD(NSTOKES,4)
      DOUBLE PRECISION XB, YB
      REAL    GRIDPOS(3,*), RADBND(NSTOKES)
      REAL    WTDO(NMU,*), MU(NMU), PHI(NMU,*)
      REAL    WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
      REAL    SFCGRIDPARMS(NSFCPAR,*), BCRAD(NSTOKES,*)
      CHARACTER SRCTYPE*1, SFCTYPE*2
      REAL GNDALBEDO

      DOUBLE PRECISION U, V
      INTEGER IL, IM, IU, IP, IBC, J
      LOGICAL LAMBERTIAN
      REAL    X(4), Y(4), RAD(NSTOKES,4), OPI
      REAL    REFLECT(4,4)
      INTEGER GRIDFACE(4,6)
      DATA    GRIDFACE/1,3,5,7, 2,4,6,8,  1,2,5,6, 3,4,7,8,
     .                 1,2,3,4, 5,6,7,8/

      LAMBERTIAN = SFCTYPE(2:2) .EQ. 'L'
      OPI = 1.0/ACOS(-1.0)
C       Loop over the four gridpoints of the cell face on the boundary
      DIRRAD = 0.0
      DO J = 1, 4
        IP = GRIDPTR(GRIDFACE(J,KFACE),ICELL)
        BOUNDPTS(J) = IP
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

          IF (LAMBERTIAN) THEN
            IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
              IF (SFCTYPE(1:1) .EQ. 'V') THEN
                DIRRAD(1,J) = OPI*SFCGRIDPARMS(2,IBC)*DIRFLUX(IP)
              ELSEIF (SFCTYPE(1:1) .EQ. 'F') THEN
                DIRRAD(1,J) = OPI*GNDALBEDO*DIRFLUX(IP)
              ENDIF
            ENDIF
          ELSE
            CALL VARIABLE_BRDF_SURFACE (NBOTPTS,IBC,IBC, BCPTR(1,2),
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,
     .             SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES,
     .             BCRAD(:,1+NTOPPTS))

            CALL SURFACE_BRDF (SFCTYPE(2:2), SFCGRIDPARMS(2,IBC),
     .                    WAVELEN, MU2, PHI2, SOLARMU,SOLARAZ,
     .                      NSTOKES, REFLECT)
            DIRRAD(:,J) = OPI*REFLECT(1:NSTOKES,1)*DIRFLUX(IP)

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
      BOUNDINTERP(1) = (1-U)*(1-V)
      BOUNDINTERP(2) = U*(1-V)
      BOUNDINTERP(3) = (1-U)*V
      BOUNDINTERP(4) = U*V
      RADBND(:) = (1-U)*(1-V)*RAD(:,1) + U*(1-V)*RAD(:,2)
     .              + (1-U)*V*RAD(:,3) +     U*V*RAD(:,4)
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
     .                              LEGEN, PHASETAB, DELTAM, NEGCHECK,
     .                              ERRMSG,IERR)
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
      CHARACTER ERRMSG*600
Cf2py intent(out) :: ERRMSG
      INTEGER IERR
Cf2py intent(out) :: IERR
      INTEGER IPH, I, J, L
      DOUBLE PRECISION  PI, OFOURPI, COSSCAT, FCT, F, X, A1, B1
      DOUBLE PRECISION, ALLOCATABLE :: UNSCLEGEN(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: DMM1(:), DMM2(:)
      IERR = 0
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
C       We are changed to using an un-scaled input so actually we just
C       need to divide by 1-f now because extinction is scaled
C       for the TMS method. -JRLoveridge 2021/04/05
        DO IPH = 1, NUMPHASE
          F = LEGEN(1,ML+1,IPH)
          DO L = 0, NLEG
            IF (L .LE. ML .AND. DELTAM) THEN
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)/(1-F)
            ELSE IF (DELTAM) THEN
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)/(1-F)
            ELSE
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)
            ENDIF
            IF (NSTLEG .GT. 1) THEN
              IF (L .LE. ML .AND. DELTAM) THEN
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)/(1-F)
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
            IERR = 1
            WRITE (ERRMSG,*) 'PRECOMPUTE_PHASE_CHECK: negative phase ',
     .        'function',
     .          ' for tabulated phase function: ','IPH',IPH,
     .        'J', J, 'A1', A1
            RETURN
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
     .                              DLEG, DPHASETAB, DELTAM, NEGCHECK,
     .                            IERR, ERRMSG)
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
      CHARACTER ERRMSG*600
Cf2py intent(out) :: ERRMSG
      INTEGER IERR
Cf2py intent(out) :: IERR
      INTEGER IPH, I, J, L
      DOUBLE PRECISION  PI, OFOURPI, COSSCAT, FCT, F, X, A1, B1
      DOUBLE PRECISION, ALLOCATABLE :: UNSCLEGEN(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: DMM1(:), DMM2(:)
      IERR = 0
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

C       Inputs are actually pre-scaled by (1-F) for DELTAM
C       so no need for anything fancy here.
        DO IPH = 1, DNUMPHASE
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
            IERR=1
            WRITE (ERRMSG,*) 'PRECOMPUTE_PHASE_CHECK_GRAD: ',
     .          'negative phase ',
     .         'function for tabulated phase function: ',
     .          'IPH',IPH,'J', J,'A1', A1
            RETURN
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

      SUBROUTINE GRADIENT_L2_OLD(NSTOKES, NX, NY, NZ, NPTS, NBPTS,
     .           NCELLS,
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
     .           WEIGHTS, UNCERTAINTIES, JACOBIAN,MAKEJACOBIAN,
     .           JACOBIANPTR, COUNTER)
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
      REAL DPATH(8*(NPX+NPY+NPZ),*), WEIGHTS(NSTOKES), WEIGHT
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

      REAL PIXEL_ERROR
      DOUBLE PRECISION RAYGRAD(NSTOKES,NBPTS,NUMDER), VISRAD(NSTOKES)
      INTEGER I, J, L, SIDE, IVIS
      LOGICAL VALIDRAD
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION U, R, PI
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT
      INTEGER M, ME, MS, NS, NS1

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

        DO NS = 1, NSTOKES
          STOKESOUT(NS,IVIS) = VISRAD(NS)
          PIXEL_ERROR = STOKESOUT(NS,IVIS) - MEASUREMENTS(NS,IVIS)
C          IF (MAKEJACOBIAN .EQV. .TRUE.) THEN
C              DO JI = 1, NBPTS
C                IF (ANY(ABS(RAYGRAD(1,JI,:))>0.0)) THEN
C                  JACOBIANPTR(1,COUNTER) = IVIS
C                  JACOBIANPTR(2,COUNTER) = JI
C                  JACOBIAN(:,COUNTER) = JACOBIAN(:,COUNTER) +
C     .            WEIGHTS(NS)*RAYGRAD(NS,JI,:)
C                  IF (NS .EQ. 1) THEN
C                    COUNTER = COUNTER + 1
C                  ENDIF
C                ENDIF
C              ENDDO
C          ENDIF
          DO NS1 = 1, NSTOKES
            WEIGHT = UNCERTAINTIES(NS,NS1,IVIS)*WEIGHTS(NS)
            GRADOUT = GRADOUT + WEIGHT*PIXEL_ERROR*RAYGRAD(NS,:,:)
            COST = COST + 0.5 * WEIGHT*PIXEL_ERROR**2
          ENDDO
        ENDDO

        IF (MAKEJACOBIAN .EQV. .TRUE.) THEN
          DO JI = 1, NBPTS
            IF (ANY(ABS(RAYGRAD(1,JI,:))>0.0)) THEN
              JACOBIANPTR(1,COUNTER) = IVIS
              JACOBIANPTR(2,COUNTER) = JI
              JACOBIAN(:,:,COUNTER) = JACOBIAN(:,:,COUNTER) +
     .        RAYGRAD(:,JI,:)
              COUNTER = COUNTER + 1
            ENDIF
          ENDDO
        ENDIF

      ENDDO
      RETURN
      END

      SUBROUTINE PREPARE_DIPHASEIND(GRIDPOS, DIPHASEIND,
     .                     NPX, NPY, NPZ, NPTS, MAXPG, NUMPHASE, DELX,
     .                      DELY, XSTART, YSTART, ZLEVELS, EXTINCTP,
     .                      ALBEDOP, NPART, NUMDER)
CPrepares the DIPHASEIND for use in the gradient calculation.
CThis holds the pointer to the property grid point from which the
C(adaptive) grid point's phase function comes from.
        IMPLICIT NONE
        INTEGER NPTS, NUMDER, NPART, MAXPG
Cf2py intent(in) :: NPTS, NUMDER, NPART, MAXPG
        INTEGER NPX, NPY, NPZ, NUMPHASE
Cf2py intent(in) :: NPX, NPY, NPZ, NUMPHASE
        REAL GRIDPOS(3, NPTS)
Cf2py intent(in) :: GRIDPOS
        INTEGER DIPHASEIND(NPTS,NUMDER)
Cf2py intent(out) :: DIPHASEIND
        REAL ALBEDOP(MAXPG, NPART), EXTINCTP(MAXPG, NPART)
Cf2py intent(in) :: ALBEDOP, EXTINCTP
        REAL DELX, DELY, XSTART, YSTART
Cf2py intent(in) :: DELX, DELY, XSTART, YSTART
        REAL ZLEVELS(*)
Cf2py intent(in) :: ZLEVELS

        INTEGER IP, IPA

        DO IPA=1,NUMDER
          DO IP=1,NPTS
            CALL GET_PHASE_INTERP_INDS(GRIDPOS(1,IP), GRIDPOS(2,IP),
     .         GRIDPOS(3,IP), NPX, NPY, NPZ, NUMPHASE, DELX,
     .         DELY, XSTART, YSTART, ZLEVELS, EXTINCTP(:,IPA),
     .         ALBEDOP(:,IPA), DIPHASEIND(IP,IPA))
          ENDDO
        ENDDO
        RETURN
      END

      SUBROUTINE GET_PHASE_INTERP_INDS (X, Y, Z,
     .                 NPX, NPY, NPZ, NUMPHASE, DELX, DELY,
     .                 XSTART, YSTART, ZLEVELS, EXTINCTP,
     .                 ALBEDOP, DIPHASEIND)
C    Adapated from TRILIN_INTERP_PROP - JRLoveridge 2021/02/10
      IMPLICIT NONE
      INTEGER DIPHASE
      REAL    X, Y, Z,  DEXT, DALB, DSCAT
      INTEGER IX, IXP, IY, IYP, IZ, IL, IM, IU, J
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I
      DOUBLE PRECISION U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
      DOUBLE PRECISION SCAT1,SCAT2,SCAT3,SCAT4,SCAT5,SCAT6,SCAT7,SCAT8
      DOUBLE PRECISION SCATTER, MAXSCAT

      INTEGER NPX, NPY, NPZ
      INTEGER NUMPHASE, DIPHASEIND
      REAL DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)
      REAL EXTINCTP(*), ALBEDOP(*)
      REAL  EXTINCT, ALBEDO

C         Find the grid location and compute the interpolation factors
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
        WRITE (6,*) 'INTERP_DERIVS: Beyond X domain',IX,NPX,X,XSTART
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
        WRITE (6,*) 'INTERP_DERIVS: Beyond Y domain',IY,NPY,Y,YSTART
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

      SCAT1 = F1*EXTINCTP(I1)*ALBEDOP(I1)
      SCAT2 = F2*EXTINCTP(I2)*ALBEDOP(I2)
      SCAT3 = F3*EXTINCTP(I3)*ALBEDOP(I3)
      SCAT4 = F4*EXTINCTP(I4)*ALBEDOP(I4)
      SCAT5 = F5*EXTINCTP(I5)*ALBEDOP(I5)
      SCAT6 = F6*EXTINCTP(I6)*ALBEDOP(I6)
      SCAT7 = F7*EXTINCTP(I7)*ALBEDOP(I7)
      SCAT8 = F8*EXTINCTP(I8)*ALBEDOP(I8)
      SCATTER = SCAT1+SCAT2+SCAT3+SCAT4+SCAT5+SCAT6+SCAT7+SCAT8

C         For tabulated phase functions pick the one we are on top of
C         or the one with the most scattering weight.

      IF (NUMPHASE .GT. 0) THEN
        MAXSCAT = -1.0
        IF (SCAT1 .GT. MAXSCAT .OR. ABS(F1-1) .LT. 0.001) THEN
          MAXSCAT = SCAT1
          DIPHASEIND = I1
        ENDIF
        IF (SCAT2 .GT. MAXSCAT .OR. ABS(F2-1) .LT. 0.001) THEN
          MAXSCAT = SCAT2
          DIPHASEIND = I2
        ENDIF
        IF (SCAT3 .GT. MAXSCAT .OR. ABS(F3-1) .LT. 0.001) THEN
          MAXSCAT = SCAT3
          DIPHASEIND = I3
        ENDIF
        IF (SCAT4 .GT. MAXSCAT .OR. ABS(F4-1) .LT. 0.001) THEN
          MAXSCAT = SCAT4
          DIPHASEIND = I4
        ENDIF
        IF (SCAT5 .GT. MAXSCAT .OR. ABS(F5-1) .LT. 0.001) THEN
          MAXSCAT = SCAT5
          DIPHASEIND = I5
        ENDIF
        IF (SCAT6 .GT. MAXSCAT .OR. ABS(F6-1) .LT. 0.001) THEN
          MAXSCAT = SCAT6
          DIPHASEIND = I6
        ENDIF
        IF (SCAT7 .GT. MAXSCAT .OR. ABS(F7-1) .LT. 0.001) THEN
          MAXSCAT = SCAT7
          DIPHASEIND = I7
        ENDIF
        IF (SCAT8 .GT. MAXSCAT .OR. ABS(F8-1) .LT. 0.001) THEN
          MAXSCAT = SCAT8
          DIPHASEIND = I8
        ENDIF
      ELSE
        WRITE (6,*) 'PHASE_INTERP_INDS: Gridded phase',
     .                 'is not supported in pyshdom.'
        STOP
      ENDIF
      RETURN
      END
