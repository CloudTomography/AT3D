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
     .                   PHASEMAX, MAXNMICRO)
C    Calculates the Stokes Vector at the given directions (CAMMU, CAMPHI)
C    and positions CAMX,CAMY,CAMZ by integrating the source function.

Cf2py threadsafe
      IMPLICIT NONE
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
      INTEGER ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
Cf2py intent(in) :: ML, MM, NCS, NSTLEG, NLM, NLEG, NUMPHASE, NPART
      INTEGER NMU, NPHI0MAX, NPHI0(*), MAXNMICRO
Cf2py intent(in) :: NMU, NPHI0MAX, NPHI0, MAXNMICRO
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
Cf2py intent(in) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER SHPTR(*), BCPTR(MAXNBC,2)
Cf2py intent(in) :: SHPTR, BCPTR
      INTEGER*2 CELLFLAGS(*)
      INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART)
Cf2py intent(in) :: CELLFLAGS, IPHASE
      REAL PHASEINTERPWT(8*MAXNMICRO,NPTS,NPART), PHASEMAX
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
     .                     PHASEINTERPWT, PHASEMAX,MAXNMICRO)
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
     .           DALB, DLEG, NSCATANGLE, YLMSUN, PHASETAB,
     .           NSTPHASE, DPHASETAB, DNUMPHASE, SOLARFLUX, NPX, NPY,
     .           NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS, EXTDIRP,
     .           UNIFORMZLEV, DPATH, DPTR, EXACT_SINGLE_SCATTER,
     .           UNCERTAINTIES, JACOBIAN,MAKEJACOBIAN,
     .           JACOBIANPTR, NUM_JACOBIAN_PTS, RAYS_PER_PIXEL,
     .           RAY_WEIGHTS, STOKES_WEIGHTS,
     .           COSTFUNC, NCOST, NGRAD, NUNCERTAINTY,
     .           TAUTOL, NODIFFUSE, IERR, ERRMSG, INTERPMETHOD,
     .           MAXNMICRO, DIPHASEP, IPHASEP, PHASEWTP, DPHASEWTP,
     .           PHASEINTERPWT, ALBEDOP, EXTINCTP, OPTINTERPWT,
     .           INTERPPTR, EXTMIN, SCATMIN, DOEXACT, TEMP, PHASEMAX,
     .           DTEMP)
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
      INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART)
Cf2py intent(in) :: CELLFLAGS, IPHASE
      REAL PHASEINTERPWT(8*MAXNMICRO,NPTS,NPART), PHASEMAX
Cf2py intent(in) ::PHASEINTERPWT, PHASEMAX
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
      INTEGER  NPIX, MAXNMICRO
Cf2py intent(in) :: NPIX, MAXNMICRO
      REAL   MEASUREMENTS(NSTOKES,*), DLEG(NSTLEG,0:NLEG,DNUMPHASE)
      REAL   DEXT(MAXPG,NUMDER), DALB(MAXPG,NUMDER)
Cf2py intent(in) :: MEASUREMENTS, DEXT ,DALB,  DLEG
      INTEGER DIPHASEP(MAXNMICRO,MAXPG,NUMDER)
      INTEGER IPHASEP(MAXNMICRO,MAXPG,NPART)
      REAL    PHASEWTP(MAXNMICRO,MAXPG,NPART)
      REAL    DPHASEWTP(MAXNMICRO,MAXPG,NUMDER)
      REAL    DTEMP(MAXPG,NUMDER)
Cf2py intent(in) :: DIPHASEP, IPHASEP, PHASEWTP, DPHASEWTP,DTEMP
      REAL    OPTINTERPWT(8,NPTS)
      INTEGER INTERPPTR(8,NPTS)
Cf2py intent(in) :: OPTINTERPWT, INTERPPTR
      REAL ALBEDOP(MAXPG,NPART), EXTINCTP(MAXPG,NPART)
Cf2py intent(in) :: ALBEDOP, EXTINCTP
      REAL TEMP(NPTS)
Cf2py intent(in) :: TEMP
      INTEGER DOEXACT(NUMDER)
Cf2py intent(in) :: DOEXACT
      DOUBLE PRECISION EXTMIN, SCATMIN
Cf2py intent(in) :: EXTMIN, SCATMIN

      REAL  STOKESOUT(NSTOKES,NPIX)
Cf2py intent(out) :: STOKESOUT
      DOUBLE PRECISION  GRADOUT(MAXPG,NUMDER,NGRAD), COST(NCOST)
Cf2py intent(out) :: GRADOUT, COST, STOKESOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1, INTERPMETHOD*2
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS, INTERPMETHOD
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
     .		       VISRAD, VALIDRAD, TOTAL_EXT, NPART, RAYGRAD,
     .		       RSHPTR, RADIANCE, LOFJ, PARTDER, NUMDER, DEXT,
     .             DALB, DLEG, MAXPG, DNUMPHASE, SOLARFLUX,
     .             NPX, NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER, TAUTOL, NODIFFUSE,
     .             GNDALBEDO, IERR, ERRMSG, TEMP, PHASEMAX,
     .             PHASEINTERPWT, OPTINTERPWT, INTERPPTR,
     .             MAXNMICRO, EXTINCTP, ALBEDOP, PHASEWTP, DPHASEWTP,
     .             IPHASEP, DIPHASEP, WAVENO, UNITS, EXTMIN, SCATMIN,
     .             DOEXACT, INTERPMETHOD, DTEMP)
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
     .             DALB, DLEG, MAXPG, DNUMPHASE, SOLARFLUX,
     .             NPX, NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER, TAUTOL, NODIFFUSE,
     .             GNDALBEDO, IERR, ERRMSG, TEMP, PHASEMAX,
     .             PHASEINTERPWT, OPTINTERPWT, INTERPPTR,
     .             MAXNMICRO, EXTINCTP, ALBEDOP, PHASEWTP, DPHASEWTP,
     .             IPHASEP, DIPHASEP, WAVENO, UNITS, EXTMIN, SCATMIN,
     .             DOEXACT, INTERPMETHOD, DTEMP)
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
      INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART), DNUMPHASE
      REAL    PHASEINTERPWT(8*MAXNMICRO,NPTS,NPART)
      REAL    TEMP(NPTS), PHASEMAX
      INTEGER BCPTR(MAXNBC,2), LOFJ(NLM)
      LOGICAL DELTAM, VALIDRAD, OUTOFDOMAIN
      REAL    WTDO(NMU,*), MU(*), PHI(NMU,*)
      REAL    WAVELEN, SOLARMU, SOLARAZ
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(*)
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(NSTLEG,0:NLEG,*)
      DOUBLE PRECISION RAYGRAD(NSTOKES,MAXPG,NUMDER)
      REAL    DIRFLUX(*), SOURCE(NSTOKES,*), RADIANCE(NSTOKES,*)
      REAL    TOTAL_EXT(*), YLMSUN(NSTLEG,*)
      REAL    PHASETAB(NSTPHASE,NUMPHASE,NSCATANGLE)
      REAL    DPHASETAB(NSTPHASE,DNUMPHASE,NSCATANGLE)
      DOUBLE PRECISION MU2, PHI2, X0,Y0,Z0, XE,YE,ZE
      DOUBLE PRECISION TRANSMIT, RADOUT(NSTOKES)
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1, INTERPMETHOD*2
      REAL      DLEG(NSTLEG,0:NLEG,*), DEXT(MAXPG,NUMDER)
      REAL      DALB(MAXPG,NUMDER), DTEMP(MAXPG,NUMDER)

      REAL      OPTINTERPWT(8,NPTS)
      INTEGER   INTERPPTR(8,NPTS), MAXNMICRO
      REAL      EXTINCTP(MAXPG,NPART), ALBEDOP(MAXPG,NPART)
      REAL      PHASEWTP(MAXNMICRO,MAXPG,NPART)
      REAL      DPHASEWTP(MAXNMICRO,MAXPG,NUMDER)
      INTEGER   IPHASEP(MAXNMICRO,MAXPG,NPART)
      INTEGER   DIPHASEP(MAXNMICRO,MAXPG,NUMDER)
      REAL      WAVENO(2)
      DOUBLE PRECISION EXTMIN, SCATMIN
      INTEGER   DOEXACT(NUMDER)

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
      REAL GRAD0(NSTOKES,8,8,NUMDER), GRAD1(NSTOKES,8,8,NUMDER)
      REAL    SRCGRAD(NSTOKES,8,8,NUMDER), SRCSINGSCAT(NSTOKES,8)
      REAL    DPATH(8*(NPX+NPY+NPZ),*), DEXTM, SECMU0
      INTEGER DPTR(8*(NPX+NPY+NPZ),*), N
      DOUBLE PRECISION PI, CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, FC(8), FB(8),FCN(8)
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
        ENDDO
C         Compute the source function times extinction in direction (MU2,PHI2)
C     Compute the quantities at each base and adaptive grid combination
C     that will be used for calculating sensitivities with respect to the
C     scattering kernel/emission/solar source.
        CALL COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR,
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, PHASEINTERPWT,
     .             DIRFLUX,
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,
     .             EXTINCT8, SRCEXT8, TOTAL_EXT, NPART, RSHPTR,
     .             RADIANCE, LOFJ, CELLFLAGS, PARTDER,
     .             NUMDER, DNUMPHASE, DEXT, DALB, DLEG, MAXPG,
     .             DSINGSCAT, GRAD8,OGRAD8, OPTINTERPWT, INTERPPTR,
     .             ALBEDOP, EXTINCTP, TEMP, PHASEWTP, IPHASEP,
     .             DIPHASEP, DPHASEWTP, SINGSCAT8, OSINGSCAT8,
     .             MAXNMICRO, DOEXACT, EXTMIN, SCATMIN,
     .             UNITS, WAVELEN, WAVENO, PHASEMAX, INTERPMETHOD,
     .             NODIFFUSE, DTEMP)

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
C       or Y. This is hardcoded to False now that we include the boundary
C       points for gradient calculations.
        OUTOFDOMAIN = .FALSE.
        IF (.NOT. OUTOFDOMAIN) THEN
          PASSEDINTERP1(:,NPASSED) = FC(:)
C         Single scatter for EXACT_SINGLE_SCATTER contribution to gradient.
          DO N=1,8
              SINGSCAT1(:,N) = FC(N)*SINGSCAT8(:,N)
              SINGSCAT1(1,N) = MAX(0.0, SINGSCAT1(1,N))
          ENDDO
C         Sum 'SOURCE' gradient contributions from each adaptive grid cell.
C         to get the total for each base grid cell.
          GRAD1 = 0.0
          DO N=1,8
            GRAD1(:,:,N,:) = FC(N)*GRAD8(:,:,N,:)
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
        CALL GET_INTERP_KERNEL(ICELL, GRIDPTR, GRIDPOS, XN, YN, ZN, FCN)
        EXTN = FCN(1)*EXTINCT8(1) + FCN(2)*EXTINCT8(2) +
     .         FCN(3)*EXTINCT8(3) + FCN(4)*EXTINCT8(4) +
     .         FCN(5)*EXTINCT8(5) + FCN(6)*EXTINCT8(6) +
     .         FCN(7)*EXTINCT8(7) + FCN(8)*EXTINCT8(8)

        TAUGRID = SO*0.5*(EXT1+EXTN)
        NTAU = MAX(1,1+INT(TAUGRID/TAUTOL))
        DELS = SO/NTAU
C           Loop over the subgrid cells
        DO IT = 1, NTAU

C         save the interpolation kernel and pointers for each subgrid
C         interval so that the radiance gradient can be calculated.
          PASSEDDELS(NPASSED) = DELS
          IF (.NOT. OUTOFDOMAIN) THEN
            PASSEDINTERP1(:,NPASSED) = FC(:)
          ELSE
            PASSEDINTERP1(:,NPASSED) = 0.0
          ENDIF
          PASSEDPOINTS(:,NPASSED) = GRIDPTR(:,ICELL)
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
            PASSEDINTERP0(:,NPASSED) = FC(:)
C           Single scatter is used for EXACT_SINGLE_SCATTER portions of gradient.
            DO N=1,8
                SINGSCAT0(:,N) = FC(N)*SINGSCAT8(:,N)
                SINGSCAT0(1,N) = MAX(0.0, SINGSCAT0(1,N))
            ENDDO

            GRAD0= 0.0
            DO N=1,8
              GRAD0(:,:,N,:) = FC(N)*GRAD8(:,:,N,:)
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
            DO KK=1,8
              IP = GRIDPTR(KK,ICELL)
              DO K=1,8
                IB = INTERPPTR(K,IP)
                RAYGRAD(:,IB,:) = RAYGRAD(:,IB,:) +
     .            TRANSMIT*SRCGRAD(:,K,KK,:)*ABSCELL
              ENDDO
              IF (EXACT_SINGLE_SCATTER .AND. (SRCTYPE .EQ. 'S'
     .            .OR. SRCTYPE .EQ. 'B')) THEN
C             Add gradient component due to the direct solar beam.
C             Sensitivity of single scattered radiation to
C             extinction along the path between the gridpoint and the sun.
                CALL COMPUTE_DIRECT_BEAM_DERIV(DPATH(:,IP),
     .            DPTR(:,IP),
     .            NUMDER, PARTDER, NPART,DEXT, TRANSMIT,
     .            ABSCELL, SRCSINGSCAT, NSTOKES, NPX,NPY,
     .            NPZ,MAXPG,ALBEDOP,PHASEWTP,LEGEN,NSTLEG,
     .            NUMPHASE,MAXNMICRO,IPHASEP,RAYGRAD,
     .            ML, NLEG, DELTAM, DNUMPHASE, DIPHASEP,
     .            DALB, EXTINCTP, DPHASEWTP, DLEG, DOEXACT)
              ENDIF
            ENDDO
          ENDIF
          NPASSED = NPASSED + 1
          IF (NPASSED .GT. MAXSUBGRIDINTS) THEN
            IERR = 1
            WRITE(ERRMSG, *) "GRAD_INTEGRATE_1RAY: The maximum ",
     .        "number of ",
     .        " subgrid intervals for calculation of the",
     .        " radiance along the",
     .        " ray path has been exceeded.", "NPASSED=", NPASSED,
     .        "MAXSUBGRDINTS=", MAXSUBGRIDINTS
            RETURN
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
C
C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
C        .OR. NGRID.GT.MAXCELLSCROSS Exceeding MAXCELLSCROSS is no
C        longer a valid reason for stopping ray integration as this
C        may compromise physical accuracy.
        IF (TRANSMIT .LT. TRANSCUT) THEN
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
          IF (EXACT_SINGLE_SCATTER .AND. (SRCTYPE .EQ. 'S'
     .            .OR. SRCTYPE .EQ. 'B')) THEN
C             Add gradient component due to the direct solar beam.
C             Sensitivity of single scattered radiation to
C             extinction along the path between the gridpoint and the sun.
            DO KK=1,4
              IP = BOUNDPTS(KK)
              CALL COMPUTE_DIRECT_BEAM_DERIV(DPATH(:,IP),
     .          DPTR(:,IP),
     .          NUMDER, PARTDER, NPART,DEXT, TRANSMIT,
     .          1.0D0, SNGL(BOUNDINTERP(KK)*DIRRAD(:,KK)),
     .          NSTOKES, NPX,NPY,
     .          NPZ,MAXPG,ALBEDOP,PHASEWTP,LEGEN,NSTLEG,
     .          NUMPHASE,MAXNMICRO,IPHASEP,RAYGRAD,
     .          ML, NLEG, DELTAM, DNUMPHASE, DIPHASEP,
     .          DALB, EXTINCTP, DPHASEWTP, DLEG, DOEXACT)
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
      CALL COMPUTE_RADIANCE_DERIVATIVE(NUMDER, PARTDER,
     .  PASSEDDELS,
     .  PASSEDPOINTS, PASSEDINTERP0, PASSEDINTERP1,
     .  TOTAL_EXT, PASSEDRAD, RAYGRAD, PASSEDTRANSMIT,
     .  PASSEDABSCELL, DEXT, DALB, LEGEN, DLEG, DIPHASEP,
     .  MAXNMICRO,
     .  IPHASEP, PHASEWTP, DPHASEWTP, INTERPMETHOD, ALBEDOP,
     .  EXTINCTP, DOEXACT, NPTS, NSTLEG, NLEG, DNUMPHASE,
     .  NUMPHASE, OPTINTERPWT, INTERPPTR, MAXPG, NPART,
     .  NSTOKES, NPASSED, DELTAM, ML, MAXSUBGRIDINTS)

      DEALLOCATE (PASSEDPOINTS, PASSEDRAD,PASSEDINTERP0,
     .            PASSEDINTERP1, PASSEDDELS, PASSEDABSCELL,
     .            PASSEDTRANSMIT)
      RETURN
      END

      SUBROUTINE COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR,
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, PHASEINTERPWT,
     .             DIRFLUX,
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,
     .             EXTINCT8, SRCEXT8, TOTAL_EXT, NPART, RSHPTR,
     .             RADIANCE, LOFJ, CELLFLAGS, PARTDER,
     .             NUMDER, DNUMPHASE, DEXT, DALB, DLEG, MAXPG,
     .             DSINGSCAT, GRAD8,OGRAD8, OPTINTERPWT, INTERPPTR,
     .             ALBEDOP, EXTINCTP, TEMP, PHASEWTP, IPHASEP,
     .             DIPHASEP, DPHASEWTP, SINGSCAT8, OSINGSCAT8,
     .             MAXNMICRO, DOEXACT, EXTMIN, SCATMIN,
     .             UNITS, WAVELEN, WAVENO, PHASEMAX, INTERPMETHOD,
     .             NODIFFUSE, DTEMP)
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
      INTEGER MAXNMICRO
      INTEGER GRIDPTR(8,*), SHPTR(*), RSHPTR(*)
      INTEGER DONETHIS(8), OLDIPTS(8), DNUMPHASE
      INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART), NPART, NUMDER, MAXPG
      REAL    PHASEINTERPWT(8*MAXNMICRO,NPTS,NPART)
      LOGICAL DELTAM
      REAL    SOLARMU, PHASEMAX
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    LEGEN(NSTLEG,0:NLEG,*), TOTAL_EXT(*)
      REAL    DIRFLUX(*), SOURCE(NSTOKES,*)
      REAL    RADIANCE(NSTOKES,*)
      REAL    YLMDIR(NSTLEG,*), YLMSUN(NSTLEG,*)
      REAL    SINGSCAT(NSTOKES,NUMPHASE)
      REAL    OEXTINCT8(8), OSRCEXT8(NSTOKES,8)
      REAL    EXTINCT8(8), SRCEXT8(NSTOKES,8), W
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
      CHARACTER SRCTYPE*1, INTERPMETHOD*2, UNITS*1
      INTEGER*2 CELLFLAGS(*)
      LOGICAL OUTOFDOMAIN, NODIFFUSE
      REAL    SINGSCAT8(NSTOKES,8), OSINGSCAT8(NSTOKES,8)
      REAL    DLEG(NSTLEG,0:NLEG,DNUMPHASE), DEXT(MAXPG,NUMDER)
      REAL    DALB(MAXPG,NUMDER)
      REAL    PHASEWTP(MAXNMICRO,MAXPG,NPART)
      REAL    DPHASEWTP(MAXNMICRO,MAXPG,NUMDER)
      REAL    DTEMP(MAXPG,NUMDER)
      REAL    EXTINCTP(MAXPG,NPART), ALBEDOP(MAXPG, NPART)
      REAL    TEMP(NPTS)
      REAL    OPTINTERPWT(8,NPTS)
      INTEGER INTERPPTR(8,NPTS)
      INTEGER IPHASEP(MAXNMICRO,MAXPG,NPART)
      INTEGER DIPHASEP(MAXNMICRO,MAXPG,NUMDER), KK
      REAL    DSINGSCAT(NSTOKES,DNUMPHASE)
      INTEGER LOFJ(*), PARTDER(NUMDER), RNS, RIS, IDR, IB,NB
      INTEGER IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I,IPA,K
      REAL    SECMU0, F, DA, A, A1, B1, EXT
      INTEGER IDP
      REAL GRAD8(NSTOKES,8,8,NUMDER), SOURCET(NSTOKES)
      REAL OGRAD8(NSTOKES,8,8,NUMDER)
      INTEGER DOEXACT(IDR)
      REAL WAVELEN, WAVENO(2)

      REAL XI
      REAL LEGENT(NSTLEG,0:NLEG)
      REAL SINGSCATJ(NSTOKES)
      REAL SCATTERJ, EXTINCTJ, ALBEDOJ
      DOUBLE PRECISION EXTMIN, SCATMIN
      REAL SPATIAL_WEIGHT
      REAL DIVIDE
      INTEGER Q

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
        ELSE IF (I .LT. 0) THEN
          EXTINCT8(N) = EXTINCT8(ABS(I))
          SRCEXT8(:,N) = SRCEXT8(:,ABS(I))
          SINGSCAT8(:,N) = SINGSCAT8(:,ABS(I))
          GRAD8(:,:,N,:) = GRAD8(:,:,ABS(I),:)
        ELSE
C         Calculate data for new points. First the gradient.
          EXT = TOTAL_EXT(IP)
          OLDIPTS(N) = IP
          IS = SHPTR(IP)
          NS = SHPTR(IP+1)-IS
          GRAD8(:,:,N,:) = 0.0

          RIS = RSHPTR(IP)
          RNS = RSHPTR(IP+1)-RIS

          DO IDR = 1, NUMDER

            IPA = PARTDER(IDR)
C             Compute the radiance/phase convolutions for each derivative
C             species and adaptive point which is used for the extinction
C             and albedo derivatives. This doesn't need to be redone
C             for each base grid point as well.
            SOURCET = 0.0
            LEGENT = 0.0

C           Compute values at radiative transfer grid points used in the
C           gradient computation at each property grid point.
            SCATTERJ = 0.0
            EXTINCTJ = 0.0
            SINGSCATJ = 0.0
            F=0.0

            DO NB=1,8
              IB = INTERPPTR(NB,IP)
              XI = OPTINTERPWT(NB,IP)
              SPATIAL_WEIGHT = XI*ALBEDOP(IB,IPA)*EXTINCTP(IB,IPA)
              SCATTERJ = SCATTERJ + SPATIAL_WEIGHT
              EXTINCTJ = EXTINCTJ + XI*EXTINCTP(IB,IPA)
              IF (SPATIAL_WEIGHT .LE. 1e-6) CYCLE

              DO Q=1,MAXNMICRO
                IF (PHASEWTP(Q,IB,IPA) .LE. 1e-6) CYCLE
                LEGENT = LEGENT +
     .              SPATIAL_WEIGHT*PHASEWTP(Q,IB,IPA)*
     .              LEGEN(:,:,IPHASEP(Q,IB,IPA))
                IF (DELTAM) THEN
                  SINGSCATJ = SINGSCATJ + SPATIAL_WEIGHT*
     .               PHASEWTP(Q,IB,IPA)*
     .               SINGSCAT(:,IPHASEP(Q,IB,IPA))
                ENDIF
              ENDDO
            ENDDO


C           DIVIDE is used in albedo gradient below.
            IF (EXTINCTJ .GT. EXTMIN) THEN
              ALBEDOJ = SCATTERJ/EXTINCTJ
            ELSE
              ALBEDOJ = SCATTERJ/EXTMIN
            ENDIF
            IF (SCATTERJ .GT. SCATMIN) THEN
              LEGENT = LEGENT/SCATTERJ
              SINGSCATJ = SINGSCATJ/SCATTERJ
            ELSE
              LEGENT = LEGENT/SCATMIN
              SINGSCATJ = SINGSCATJ/SCATMIN
            ENDIF

C           Do delta-M scaling.
            IF (DELTAM) THEN
              F = LEGENT(1,ML+1)
              IF (INTERPMETHOD(2:2) .EQ. 'N') THEN
                LEGENT(:,0:ML) = LEGENT(:,0:ML)/(1-F)
              ENDIF
            ENDIF
            DIVIDE = 1.0D0/(1 - F*ALBEDOJ)
C           Calculate the phase/radiance component of source function.
C           for extinction/albedo gradients.
C           (this is just for one SPECIES so we don't reuse it.
C           This would be a wasted calculation in the special case of one
C           SPECIES (NPART=1) so that SOURCE already contains what we need.)

            IF (SCATTERJ .GT. SCATMIN) THEN
C             don't bother with the computation if the scatter coefficient is
C             gonna be zero anyway.
              CALL COMPUTE_SOURCE_DIRECTION(LEGENT, SOURCET,
     .          NLM, LOFJ, RIS, RNS, NSTOKES, RADIANCE, YLMDIR,
     .          NSTLEG, DELTAM, SRCTYPE, ML, MM, NS, YLMSUN,
     .          DIRFLUX(IP), SECMU0, NLEG)
              IF (DELTAM .AND. (SRCTYPE .EQ. 'S' .OR.
     .            SRCTYPE .EQ. 'B')) THEN
                SOURCET = SOURCET + DIRFLUX(IP)*SINGSCATJ*SECMU0/(1-F)
              ENDIF
              SOURCET(1) = MAX(0.0, SOURCET(1))
            ENDIF

C           Undo delta-M scaling so we have the unscaled legendre/wigner coefficients
C           on the RTE grid for phase function gradients.
            IF (DELTAM) THEN
              LEGENT(:,0:ML) = LEGENT(:,0:ML)*(1-F)
              LEGENT(1,0:ML) = LEGENT(1,0:ML) + F
              IF (NSTLEG .GT. 1) THEN
                LEGENT(2:4,0:ML) = LEGENT(2:4,0:ML) + F
              ENDIF
            ENDIF

C           Loop again through each property grid point, this time calculating the
C           derivatives using what we prepared. Gradient is updated in GRAD8
C           for the specified property grid point and species.
            DO NB = 1,8
              IB = INTERPPTR(NB,IP)
              XI = OPTINTERPWT(NB,IP)
              CALL COMPUTE_GRADIENT_ONEPROPPOINT(XI, ALBEDOP(IB,IPA),
     .          EXTINCTP(IB,IPA), DEXT(IB,IDR), DALB(IB,IDR), MAXNMICRO,
     .          DLEG, DIPHASEP(:,IB,IDR), PHASEWTP(:,IB,IPA),
     .          IPHASEP(:,IB,IPA), DPHASEWTP(:,IB,IDR), ALBEDO(IP,IPA),
     .          EXTINCT(IP,IPA), LEGEN, NUMPHASE, DNUMPHASE, NLEG,
     .          NSTLEG, RADIANCE, NSTOKES, DELTAM, DOEXACT(IDR),
     .          NLM, LOFJ, RIS,
     .          RNS, ML, MM, NS, YLMSUN, DIRFLUX(IP), SECMU0, YLMDIR,
     .          SRCTYPE, SINGSCAT, DSINGSCAT, DTEMP(IB,IDR),
     .          TEMP(IP), INTERPMETHOD, GRAD8(:,NB,N,IDR), F,
     .          ALBEDOJ, DIVIDE, SINGSCATJ, SOURCET, LEGENT, UNITS,
     .          WAVELEN, WAVENO)
            ENDDO
          ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         End of gradient calculations.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C         Now we perform the calculations for the forward model, calculating thw
C         SRCEXT product. This makes use of the already calculated SOURCE
C         and only has to consider DELTAM correction.

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

C         Special case for solar source and delta-M
          IF ((SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B')
     .      .AND. DELTAM) THEN

           DO IPA = 1, NPART
             IF (EXT.EQ.0.0) THEN
                 W = 1.0
             ELSE
                 W = EXTINCT(IP,IPA)/EXT
             ENDIF
             IF (W.EQ.0.0) CYCLE

             IF (INTERPMETHOD(2:2) .EQ. 'O' ) THEN
               LEGENT = LEGEN(:,:,IPHASE(1,IP,IPA))
               F = LEGENT(1,ML+1)
             ELSEIF (INTERPMETHOD(2:2) .EQ. 'N') THEN
               IF (PHASEINTERPWT(1,IP,IPA) .GE. PHASEMAX) THEN
                 LEGENT = LEGEN(:,:,IPHASE(1,IP,IPA))
               ELSE
                 LEGENT = 0.0
                 DO Q=1,8*MAXNMICRO
                   IF (PHASEINTERPWT(Q,IP,IPA) .LE. 1e-5) CYCLE
                   LEGENT = LEGENT + LEGEN(:,:,IPHASE(Q,IP,IPA))*
     .              PHASEINTERPWT(Q,IP,IPA)
                 ENDDO
               ENDIF
               F = LEGENT(1,ML+1)
               LEGENT(:,0:ML) = LEGENT(:,0:ML)/(1-F)
             ENDIF
C               First subtract off the truncated single scattering
             DA = ALBEDO(IP,IPA)*DIRFLUX(IP)*SECMU0*W
             J = 1

             DO L = 0, ML
               ME = MIN(L,MM)
               MS = -ME
               A1 = DA*LEGENT(1,L)
               B1 = DA*LEGENT(5,L)
               IF (J .LE. NS) THEN
                 JT = J
                 DO M = MS, ME
                   SRCEXT8(1,N) =SRCEXT8(1,N)
     .               -A1*YLMDIR(1,J)*YLMSUN(1,J)
                   J = J + 1
                 ENDDO
                 IF (NSTOKES .GT. 1) THEN
                   J = JT
                   DO M = MS, ME
                     SRCEXT8(2,N)=SRCEXT8(2,N)
     .                 -B1*YLMDIR(2,J)*YLMSUN(1,J)
                     SRCEXT8(3,N)=SRCEXT8(3,N)
     .                 -B1*YLMDIR(6,J)*YLMSUN(1,J)
                     J = J + 1
                  ENDDO
                ENDIF
              ENDIF
             ENDDO

C               Then add in the single scattering contribution for the
C               original unscaled phase function.
             IF (NUMPHASE .GT. 0) THEN
               IF (PHASEINTERPWT(1,IP,IPA) .GE. PHASEMAX) THEN
                 SINGSCAT8(:,N) = SINGSCAT8(:,N) +
     .          DA*SINGSCAT(:,IPHASE(1,IP,IPA))/(1-F)
               ELSE
                 DO Q=1,8*MAXNMICRO
                   IF (PHASEINTERPWT(Q,IP,IPA) .LE. 1e-5) CYCLE
                   SINGSCAT8(:,N) = SINGSCAT8(:,N) +
     .              DA*SINGSCAT(:,IPHASE(Q,IP,IPA))*
     .              PHASEINTERPWT(Q,IP,IPA)/(1-F)
                 ENDDO
               ENDIF
             ELSE
               WRITE(6,*) 'NUMPHASE=0 is not supported.'
               STOP
             ENDIF
           ENDDO
           SRCEXT8(:,N) = SRCEXT8(:,N) + SINGSCAT8(:,N)
          ENDIF
          SINGSCAT8(:,N) = SINGSCAT8(:,N)*EXT
C         The SRCEXT8 is used to calculate the radiance, which is itself
C         used in the gradient. If we only want the single scatter (NODIFFUSE) then
C         the source should only include the single scatter so we set that
C         here. This doesn't save computation time at all though.
          IF (NODIFFUSE) THEN
            SRCEXT8(:,N) = SINGSCAT8(:,N)
          ELSE
            SRCEXT8(:,N) = SRCEXT8(:,N)*EXT
          ENDIF
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
C       for the TMS method. -JRLoveridge 2021/04/05'

C       Why are we unscaling LEGEN for each scat angle????
C       Isn't this a waste of time?
        DO IPH = 1, NUMPHASE
CF = LEGEN(1,ML+1,IPH)/(2*(ML+1)+1)
C          IF (J .EQ. 1) THEN
C            PRINT *, 'HELLO',IPH, F, LEGEN(1,ML+1,IPH)
C          ENDIF
          DO L = 0, NLEG
            IF (L .LE. ML .AND. DELTAM) THEN
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)/(2*L+1)
C((2*L+1)*(1-F))
            ELSE IF (DELTAM) THEN
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)/(2*L+1)
            ELSE
              UNSCLEGEN(1,L) = LEGEN(1,L,IPH)/(2*L+1)
            ENDIF
C            IF (J .EQ. 1) THEN
C              PRINT *, UNSCLEGEN(1,L), LEGEN(1,L,IPH)
C            ENDIF
            IF (NSTLEG .GT. 1) THEN
              IF (L .LE. ML .AND. DELTAM) THEN
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)/(2*L+1)
              ELSE IF (DELTAM) THEN
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)/(2*L+1)
              ELSE
                UNSCLEGEN(2,L) = LEGEN(5,L,IPH)/(2*L+1)
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


      SUBROUTINE COMPUTE_GRADIENT_ONEPROPPOINT(XI, ALBEDOP, EXTINCTP,
     .  DEXT, DALB, MAXNMICRO, DLEG, DIPHASEP, PHASEWTP,
     .  IPHASEP, DPHASEWTP, ALBEDO, EXTINCT,
     .  LEGEN, NUMPHASE, DNUMPHASE, NLEG, NSTLEG, RADIANCE,
     .  NSTOKES, DELTAM, DOEXACT, NLM, LOFJ, RIS, RNS,
     .  ML, MM, NS, YLMSUN, DIRFLUX, SECMU0, YLMDIR,
     .  SRCTYPE, SINGSCAT, DSINGSCAT, DTEMP, TEMP,
     . INTERPMETHOD, GRADTEMP, F, ALBEDOJ, DIVIDE, SINGSCATJ,
     . SOURCET, LEGENT, UNITS, WAVELEN, WAVENO)
C      Calculates the gradient of the source/extinction product at a given RTE grid point
C      for a given species with respect to unknowns on a given property grid point
C      using partial derivatives of optical properties on property grid with
C      respect to unknowns (DEXT, DALB, DLEG, DIPHASEP, DPHASEWTP, DSINGSCAT)
C      in combination with optical properties on property grid
C      (ALBEDOP, EXTINCTP, IPAHSEP, PHASEWTP, LEGEN, TEMPP).
C      XI is the interpolation weight between the property grid point and the
C      RTE grid point.
C      (SOURCET, LEGENT, SINGSCATJ, ALBEDOJ, F, PLANCK, DIVIDE) are RTE grid point values.
C      ALBEDO and EXTINCT are delta-M scaled RTE grid point values.
C      The gradient contribution is updated in GRADTEMP so GRADTEMP needs to be
C      initialized outside of this subroutine.
C      The complexity of the gradient calculations comes from the fact that
C      delta-M scaling is done on the RTE grid, which is distinct from the property grid.
C      Note that at zero albedo, the derivative with respect to albedo is not
C      currently well defined.
C      This is because the phase function (and delta-M scaling) are undefined as well
C      and both change simultaneously when albedo becomes non-zero.
C      The gradient works fine for very very small albedo (e.g. 1e-7).
C      but care should be taken when transitioning between exactly non-scattering
C      and partially non-scattering.
      IMPLICIT NONE

      INTEGER MAXNMICRO, NUMPHASE, DNUMPHASE, NLEG, NSTLEG
Cf2py intent(in) :: MAXNMICRO, NUMPHASE, DNUMPHASE, NLEG, NSTLEG
      INTEGER ML, MM, NS, NLM, DOEXACT, NSTOKES
Cf2py intent(in) :: ML, MM, NS, NLM, DOEXACT, NSTOKES
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL DIRFLUX, SECMU0, YLMDIR(NSTLEG, NLM), YLMSUN(NSTLEG, NLM)
Cf2py intent(in) :: DIRFLUX, SECMU0, YLMDIR, YLMSUN
      REAL XI, ALBEDOP, EXTINCTP, DEXT, DALB, DTEMP, TEMP
Cf2py intent(in) :: XI, ALBEDOP, EXTINCTP, DEXT, DALB, DTEMP, TEMP
      REAL ALBEDO, EXTINCT
Cf2py intent(in) :: ALBEDO, EXTINCT
      INTEGER DIPHASEP(MAXNMICRO), IPHASEP(MAXNMICRO)
Cf2py intent(in) :: DIPHASEP, IPHASEP
      REAL PHASEWTP(MAXNMICRO), DPHASEWTP(MAXNMICRO)
Cf2py intent(in) :: PHASEWTP, DPHASEWTP
      REAL SINGSCAT(NSTOKES, NUMPHASE), DSINGSCAT(NSTOKES, DNUMPHASE)
Cf2py intent(in) :: SINGSCAT, DSINGSCAT
      REAL RADIANCE(NSTOKES, *)
Cf2py intent(in) :: RADIANCE
      INTEGER RIS, RNS, LOFJ(NLM)
Cf2py intent(in) :: RIS, RNS, LOFJ
      CHARACTER SRCTYPE*1, INTERPMETHOD*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, INTERPMETHOD, UNITS
      REAL LEGEN(NSTLEG,0:NLEG, NUMPHASE)
Cf2py intent(in) :: LEGEN
      REAL DLEG(NSTLEG,0:NLEG, DNUMPHASE)
Cf2py intent(in) :: DLEG
      REAL LEGENT(NSTLEG, 0:NLEG), SOURCET(NSTOKES)
Cf2py intent(in) :: LEGENT, SOURCET
      REAL SINGSCATJ(NSTOKES), F, ALBEDOJ, DIVIDE
Cf2py intent(in) :: SINGSCATJ, F, ALBEDOJ, DIVIDE
      REAL WAVELEN, WAVENO(2)
Cf2py intent(in) :: WAVELEN, WAVENO
      REAL GRADTEMP(NSTOKES)
Cf2py intent(in,out) :: GRADTEMP

      INTEGER Q
      REAL DSINGSCATP(NSTOKES), FA, FP, DSOURCE(NSTOKES)
      REAL SINGSCATP(NSTOKES), DPLANCK, PLANCK
      REAL DFJ, EXTINCT_GRAD, ALBEDO_GRAD, FTEMP
      REAL, ALLOCATABLE :: DLEGP(:,:), UNSCALED_LEGEN(:,:), DLEGT(:,:)
      REAL, ALLOCATABLE :: LEGENP(:,:)
      ALLOCATE(DLEGP(NSTLEG, 0:NLEG), UNSCALED_LEGEN(NSTLEG,0:NLEG))
      ALLOCATE(DLEGT(NSTLEG, 0:NLEG), LEGENP(NSTLEG, 0:NLEG))

C     Prepare phase function derivative on property grid (DLEGP),
C     derivative of deltam scaling on property grid (FA)
C     deltam scaling on property grid (FP)
C     and single scatter derivative on property grid (if DELTAM) (DSINGSCATP)
      DLEGP = 0.0
      LEGENP = 0.0
      FA = 0.0
      FP = 0.0
      DSINGSCATP = 0.0

C     Loop over each of the phase functions whose weighted sum form the
C     phase function at the property grid.
      DO Q=1,MAXNMICRO
C       gradient of phase function at this point due to
C       individual gradients in phase functions. This is either
C       one sum over all appropriate phase function gradients (DLEG)
C       if a single exact phase function is supplied for each property
C       grid point for this variable.
C       Otherwise the phase function at this point is represented by
C       linear interpolation for this variable. In this case we calculate
C       derivative using a weighted sum over phase functions with the
C       weights held in DPHASEWTP.
        IF (DOEXACT .EQ. 1) THEN
          IF (DELTAM) THEN
            FA = FA + PHASEWTP(Q)*
     .        DLEG(1,ML+1,DIPHASEP(Q))
            DSINGSCATP = DSINGSCATP + PHASEWTP(Q)*
     .         DSINGSCAT(:,DIPHASEP(Q))
          ENDIF
          DLEGP = DLEGP + PHASEWTP(Q)*
     .        DLEG(:,:,DIPHASEP(Q))
        ENDIF

C       Prepare the unscaled LEGEN for this property grid point.
        UNSCALED_LEGEN = LEGEN(:,:,IPHASEP(Q))
        IF (DELTAM) THEN
          FTEMP = UNSCALED_LEGEN(1,ML+1)
          IF (INTERPMETHOD(2:2) .EQ. 'O') THEN
            UNSCALED_LEGEN(1,0:ML) =
     .        UNSCALED_LEGEN(1,0:ML)*(1-FTEMP)
          ENDIF
          UNSCALED_LEGEN(1,0:ML) =
     .      UNSCALED_LEGEN(1,0:ML) + FTEMP
          IF (NSTLEG .GT. 1) THEN
            UNSCALED_LEGEN(2:4,0:ML) =
     .        UNSCALED_LEGEN(2:4,0:ML) + FTEMP
          ENDIF
        ENDIF
        LEGENP = LEGENP +
     .    PHASEWTP(Q)*UNSCALED_LEGEN
        IF (DOEXACT .EQ. 0) THEN
          DLEGP = DLEGP +
     .        DPHASEWTP(Q)*UNSCALED_LEGEN

          IF (DELTAM) THEN
            FA = FA + DPHASEWTP(Q)*LEGEN(1,ML+1,IPHASEP(Q))
            DSINGSCATP = DSINGSCATP +
     .        DPHASEWTP(Q)*SINGSCAT(:,IPHASEP(Q))
          ENDIF
        ENDIF
      ENDDO

C     Calculate some components of extinct / albedo / phase
C     gradients that will be used again.

C      Calculate the delta-M scaling factor for this property grid point.
      IF (DELTAM) THEN
        FP = LEGENP(1,ML+1)
      ENDIF

C     Partial derivative of delta-M scaled extinction with respect
C     to the unknown.
      EXTINCT_GRAD = DEXT*(1-FP*ALBEDOP) -DALB*FP*EXTINCTP
     .            -EXTINCTP*ALBEDOP*FA

C     Partial derivative of delta-M scaled albedo MULTIPLIED
C     BY THE DELTAM SCALED EXTINCTION.
      ALBEDO_GRAD = DIVIDE*(
     .           DEXT*((1-F)*(ALBEDOP - ALBEDOJ)
     .            +(ALBEDOJ - 1)*ALBEDOP*(FP - F))
     .          +DALB*((1-F)*EXTINCTP
     .            +(ALBEDOJ - 1)*EXTINCTP*(FP - F))
     .          +FA*(ALBEDOJ - 1)*EXTINCTP*ALBEDOP
     .          )

C     Initialize the phase part of the gradient.
      DSOURCE = 0.0

C     Calculate the delta-M part of the phase gradient.
      DFJ = DEXT*(FP-F)*ALBEDOP + DALB*(FP-F)*EXTINCTP +
     .       FA*EXTINCTP*ALBEDOP

C     DLEGT is the phase derivative MULTIPLIED BY THE DELTAM SCALED
C     EXTINCTION/ALBEDO product.
      DLEGT = DEXT*(LEGENP-LEGENT)*ALBEDOP +
     .            DALB*(LEGENP-LEGENT)*EXTINCTP+
     .            DLEGP*EXTINCTP*ALBEDOP +
     .            (LEGENT-1)*DFJ/(1-F)

C     Weight each phase legendre/Wigner coefficient derivative
C     by the corresponding radiance harmonics.
C     Note that this may be slightly inconsistent as we use
C     all of the radiance harmonics in the derivative calculation,
C     which is typically only a few more terms than
C     in the SOURCE function but may be much more if the
C     source function truncation is active while HIGHORDERRAD is True.
C     (HIGHORDERRAD keeps the full set of NLM radiance harmonics.)
C     On the other hand, the pseudo-solar source
C     uses the exact same number of terms as SOURCE.
      CALL COMPUTE_SOURCE_DIRECTION(DLEGT, DSOURCE,
     .   NLM, LOFJ, RIS, RNS, NSTOKES, RADIANCE, YLMDIR,
     .   NSTLEG, DELTAM, SRCTYPE, ML, MM, NS, YLMSUN,
     .   DIRFLUX, SECMU0, NLEG)

C      Add the extinction/albedo/phase gradients.
      GRADTEMP(:) = GRADTEMP(:) +
     .  XI*(SOURCET(:)*(ALBEDO*EXTINCT_GRAD + ALBEDO_GRAD)
     .      + DSOURCE(:))

C     Now calculate the solar single scatter component of the gradient.
C     The non delta-M case is covered inside DSOURCE, see
C     CALCULATE_SOURCE_DIRECTION. Like DLEGT, this is already multiplied
C     by the deltam scaled albedo/extinct product.
      IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
        IF (DELTAM) THEN
          SINGSCATP = 0.0
          DO Q=1,MAXNMICRO
            SINGSCATP(:) = SINGSCATP(:)
     .        + PHASEWTP(Q)*SINGSCAT(:,IPHASEP(Q))
          ENDDO
          GRADTEMP(:) = GRADTEMP(:) +
     .      DIRFLUX*SECMU0*XI*(SINGSCATJ*DFJ/(1-F) +
     .          DSINGSCATP*EXTINCTP*ALBEDOP
     .        + DALB*(SINGSCATP - SINGSCATJ)*EXTINCTP
     .        + DEXT*(SINGSCATP - SINGSCATJ)*ALBEDOP)
        ENDIF
      ENDIF

C      Calculate thermal component of the gradient.
      IF (SRCTYPE .EQ. 'T' .OR. SRCTYPE .EQ. 'B') THEN
C        DPLANCK COULD BE PREPROCESSED.
         CALL PLANCK_DERIVATIVE(TEMP, UNITS, WAVENO,
     .                          WAVELEN, DPLANCK)
C        To capture the case when ALBEDO ~ 1, but we still
C        want a derivative we just recalculate PLANCK
C        without the (1-ALBEDO) factor, rather than
C        trying to divide, as it is otherwise undefined.
         CALL PLANCK_FUNCTION(TEMP, UNITS, WAVENO,
     .                          WAVELEN, PLANCK)
C       Catch NaNs that occur when PLANCK is very small.
         IF (PLANCK .LT. 1e-9) THEN
           DPLANCK = 0.0
         ENDIF
          GRADTEMP(1) = GRADTEMP(1) +
     .      XI*(EXTINCT*(1.0 - ALBEDO)*DPLANCK*DTEMP
     .          - PLANCK*ALBEDO_GRAD
     .         + PLANCK*(1.0-ALBEDO)*EXTINCT_GRAD)

      ENDIF

      DEALLOCATE(DLEGP, UNSCALED_LEGEN, DLEGT, LEGENP)
      RETURN
      END

      SUBROUTINE COMPUTE_RADIANCE_DERIVATIVE(NUMDER, PARTDER,
     .  PASSEDDELS, PASSEDPOINTS, PASSEDINTERP0, PASSEDINTERP1,
     .  TOTAL_EXT, PASSEDRAD, RAYGRAD, PASSEDTRANSMIT, PASSEDABSCELL,
     .  DEXT, DALB, LEGEN, DLEG, DIPHASEP, MAXNMICRO,
     .  IPHASEP, PHASEWTP, DPHASEWTP, INTERPMETHOD, ALBEDOP,
     .  EXTINCTP, DOEXACT, NPTS, NSTLEG, NLEG, DNUMPHASE,
     .  NUMPHASE, OPTINTERPWT, INTERPPTR, MAXPG, NPART,
     .  NSTOKES, NPASSED, DELTAM, ML, MAXSUBGRIDINTS)
      IMPLICIT NONE
C     Computes the inner product of the adjoint direct beam with
C     product of radiance and extinction (which is a part of the sensitivity
C     source.
C     Variables with PASSED are defined at each of the subgrid intervals
C     used in the ray integration. By the end of the ray integration
C     PASSEDRAD contains the radiance at the edge of each subgrid interval.
C     The same integration formula is used as in the source function integration
C     which calculates the adjoint direct beam assuming that extinction varies
C     linearly and assumes that the extinction/radiance product also varies linearly.
C     The assumption of linear radiance variation will induce error if the
C     subgrid intervals are not sufficiently small as it is not perfectly
C     consistent with the assumption of the source-extinction product varying
C     linearly over the characteristic. At least, it does not
C     appear to meaningfully cause error with tautol=0.2 or smaller.
      INTEGER NPTS, NSTLEG, NLEG, MAXNMICRO, MAXPG, NUMDER, NPART
      INTEGER PARTDER(NUMDER), NSTOKES, NUMPHASE, DNUMPHASE
      INTEGER NPASSED, ML,MAXSUBGRIDINTS
      LOGICAL DELTAM
      CHARACTER INTERPMETHOD*2
      REAL DEXT(MAXPG,NUMDER), DALB(MAXPG,NUMDER)
      REAL ALBEDOP(MAXPG,NPART), EXTINCTP(MAXPG,NPART)
      REAL OPTINTERPWT(8,NPTS)
      INTEGER INTERPPTR(8,NPTS)
      REAL LEGEN(NSTLEG,0:NLEG,NUMPHASE)
      REAL PHASEWTP(MAXNMICRO,MAXPG,NPART)
      REAL DLEG(NSTLEG,0:NLEG,DNUMPHASE)
      REAL DPHASEWTP(MAXNMICRO,MAXPG,NUMDER)
      INTEGER DIPHASEP(MAXNMICRO,MAXPG,NUMDER)
      INTEGER IPHASEP(MAXNMICRO,MAXPG,NPART)
      INTEGER DOEXACT(NUMDER)
      REAL TOTAL_EXT(NPTS)
      DOUBLE PRECISION RAYGRAD(NSTOKES,MAXPG,NUMDER)
      DOUBLE PRECISION PASSEDRAD(NSTOKES,MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDABSCELL(MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDTRANSMIT(MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDINTERP0(8,MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDINTERP1(8,MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDDELS(MAXSUBGRIDINTS)
      INTEGER PASSEDPOINTS(8,MAXSUBGRIDINTS)

      INTEGER IDR, IPA, IP, IB, K, NB, KK
      DOUBLE PRECISION EXT0, EXT1, DELS, EXT
      REAL XI, EXTGRAD, RADGRAD0(NSTOKES)
      REAL RADGRAD1(NSTOKES), RADGRAD(NSTOKES)

      DO IDR =1,NUMDER
        IPA=PARTDER(IDR)
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
            RADGRAD0 = 0.0
            RADGRAD1 = 0.0
            DO K=1,8
              IP = PASSEDPOINTS(K,KK)
              DO NB=1,8
                IB = INTERPPTR(NB,IP)
                XI = OPTINTERPWT(NB,IP)
                CALL EXTINCTION_DERIVATIVE_POINT(
     .            DEXT(IB,IDR), DOEXACT(IDR), XI, LEGEN,
     .            NLEG, NSTLEG, NUMPHASE, DNUMPHASE,
     .            PHASEWTP(:,IB,IPA),
     .            DELTAM, IPHASEP(:,IB,IPA),
     .            DIPHASEP(:,IB,IDR),EXTGRAD, DLEG,
     .            MAXNMICRO, DALB(IB,IDR), ALBEDOP(IB,IPA),
     .            EXTINCTP(IB,IPA), DPHASEWTP(:,IB,IDR), ML)

                RADGRAD0(:) = -1*PASSEDRAD(:,KK+1)*
     .                    EXTGRAD*PASSEDINTERP0(K,KK)
                RADGRAD1(:) = -1*PASSEDRAD(:,KK)*
     .                    EXTGRAD*PASSEDINTERP1(K,KK)
                RADGRAD = ( 0.5*(RADGRAD0+RADGRAD1)
     .           + 0.08333333333*(EXT0*RADGRAD1-EXT1*RADGRAD0)*DELS
     .                *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
                RAYGRAD(:,IB,IDR) = RAYGRAD(:,IB,IDR) +
     .           RADGRAD*PASSEDTRANSMIT(KK)*PASSEDABSCELL(KK)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END


      SUBROUTINE EXTINCTION_DERIVATIVE_POINT(DEXT, DOEXACT, XI,
     .      LEGEN, NLEG, NSTLEG, NUMPHASE, DNUMPHASE, PHASEWTP,
     .      DELTAM, IPHASEP, DIPHASEP, EXTGRAD, DLEG, MAXNMICRO,
     .      DALB, ALBEDOP, EXTINCTP, DPHASEWTP, ML)
C     Evaluate the derivative of the delta-M scaled extinction
C     with respect to the unscaled optical properties on the
C     property grid.
C     This is only used in the computation of the direct beam
C     derivatives (COMPUTE_DIRECT_BEAM_DERIV) though it is
C     identical to that used in COMPUTE_GRADIENT_ONEPROPPOINT.
C     The extinction gradient could be precomputed as it only requires
C     property grid variables and is typically recalculated several
C     times at each grid point for the direct beam derivatives.
      IMPLICIT NONE
      INTEGER NLEG, NSTLEG, NUMPHASE, DNUMPHASE, MAXNMICRO
Cf2py intent(in) :: NLEG, NSTLEG, NUMPHASE, DNUMPHASE, MAXNMICRO
      REAL LEGEN(NSTLEG,0:NLEG,NUMPHASE), DLEG(NSTLEG,0:NLEG,DNUMPHASE)
Cf2py intent(in) :: LEGEN, DLEG
      REAL PHASEWTP(MAXNMICRO), DPHASEWTP(MAXNMICRO)
Cf2py intent(in) :: PHASEWTP, DPHASEWTP
      INTEGER IPHASEP(MAXNMICRO), DIPHASEP(MAXNMICRO)
Cf2py intent(in) :: IPHASEP, DIPHASEP
      INTEGER DOEXACT, ML
Cf2py intent(in) :: DOEXACT, ML
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL DEXT, DALB, ALBEDOP, EXTINCTP, XI, EXTGRAD
Cf2py intent(in) :: DEXT, DALB, ALBEDOP, EXTINCTP, XI
Cf2py intent(out) :: EXTGRAD

      REAL FP, DFP
      INTEGER Q
C     FP is the delta-M scale factor on the property grid.
C     DFP is the derivative of the delta-M scale factor with
C     respect to the unknowns.
      FP = 0.0
      DFP = 0.0
      DO Q=1,MAXNMICRO
        IF (PHASEWTP(Q) .LE. 1E-6) CYCLE
        IF (DELTAM) THEN
          FP = FP + PHASEWTP(Q)*
     .                LEGEN(1,ML+1,IPHASEP(Q))
          IF (DOEXACT .EQ. 1) THEN
            DFP = DFP + PHASEWTP(Q)*
     .                  DLEG(1,ML+1,DIPHASEP(Q))
          ELSEIF (DOEXACT .EQ. 0) THEN
            DFP = DFP + DPHASEWTP(Q)*
     .                  LEGEN(1,ML+1,IPHASEP(Q))
          ENDIF
        ENDIF
      ENDDO
      EXTGRAD = XI*(DEXT*(1-FP*ALBEDOP) -
     .            DALB*FP*EXTINCTP -
     .            EXTINCTP*ALBEDOP*DFP)
      RETURN
      END


      SUBROUTINE COMPUTE_DIRECT_BEAM_DERIV(DPATH, DPTR,
     .     NUMDER, PARTDER, NPART, DEXT, TRANSMIT, ABSCELL,
     .     INPUTWEIGHT, NSTOKES, NPX,NPY,NPZ,MAXPG,
     .  ALBEDOP, PHASEWTP, LEGEN, NSTLEG, NUMPHASE,
     .  MAXNMICRO, IPHASEP, RAYGRAD, ML, NLEG, DELTAM,
     .  DNUMPHASE, DIPHASEP, DALB, EXTINCTP, DPHASEWTP,
     .  DLEG, DOEXACT)
      IMPLICIT NONE
C     Computes the derivatives of the direct beam with respect
C     to the unknowns. This ensures that the derivative is
C     exact in the single scatter case.
C     This makes use of precomputed weights (DPATH) and pointers
C     (DPTR) for the integration of the direct solar beam through
C     the property grid.
C     The derivative of the extinction with respect to the unknowns
C     at each point is evaluated and used to calculate the
C     derivative. The INPUTWEIGHT is either the surface reflected radiance
C     due to the direct solar beam or the single-scatter source
C     at the grid point.
      INTEGER NSTOKES, NUMDER, NPART, PARTDER(NUMDER)
Cf2py intent(in) :: NSTOKES, NUMDER, NPART, PARTDER
      INTEGER NPX, NPY, NPZ, MAXPG, MAXNMICRO
Cf2py intent(in) :: NPX, NPY, NPZ, MAXPG, MAXNMICRO
      INTEGER NSTLEG, NUMPHASE, NLEG, ML, DNUMPHASE
Cf2py intent(in) :: NSTLEG, NUMPHASE, NLEG, ML, DNUMPHASE
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      INTEGER DOEXACT(NUMDER)
Cf2py intent(in) :: DOEXACT
      REAL DPATH(8*(NPX+NPY+NPZ))
Cf2py intent(in) :: DPATH
      INTEGER DPTR(8*(NPX+NPY+NPZ))
Cf2py intent(in) :: DPTR
      REAL DEXT(MAXPG, NUMDER), DALB(MAXPG,NUMDER)
Cf2py intent(in) :: DEXT, DALB
      REAL ALBEDOP(MAXPG,NPART), EXTINCTP(MAXPG,NPART)
Cf2py intent(in) :: ALBEDOP, EXTINCTP
      REAL PHASEWTP(MAXNMICRO,MAXPG,NPART)
Cf2py intent(in) :: PHASEWTP
      INTEGER IPHASEP(MAXNMICRO,MAXPG,NPART)
Cf2py intent(in) :: IPHASEP
      REAL DPHASEWTP(MAXNMICRO,MAXPG,NUMDER)
Cf2py intent(in) :: DPHASEWTP
      INTEGER DIPHASEP(MAXNMICRO, MAXPG, NUMDER)
Cf2py intent(in) :: DIPHASEP
      REAL LEGEN(NSTLEG,0:NLEG,NUMPHASE)
Cf2py intent(in) :: LEGEN
      REAL DLEG(NSTLEG,0:NLEG,DNUMPHASE)
Cf2py intent(in) :: DLEG
      DOUBLE PRECISION TRANSMIT, ABSCELL
Cf2py intent(in) :: TRANSMIT, ABSCELL
      REAL INPUTWEIGHT(NSTOKES)
Cf2py intent(in) :: INPUTWEIGHT
      DOUBLE PRECISION RAYGRAD(NSTOKES,MAXPG,NUMDER)
Cf2py intent(in,out) :: RAYGRAD

      INTEGER IB, II, IDR, Q, IPA
      REAL FP, EXTGRAD

      II=1
      DO WHILE (DPTR(II) .GT. 0)
        IB = DPTR(II)
        DO IDR=1,NUMDER
          IPA = PARTDER(IDR)
          CALL EXTINCTION_DERIVATIVE_POINT(
     .     DEXT(IB,IDR), DOEXACT(IDR), 1.0, LEGEN,
     .     NLEG, NSTLEG, NUMPHASE, DNUMPHASE,
     .     PHASEWTP(:,IB,IPA),
     .     DELTAM, IPHASEP(:,IB,IPA),
     .     DIPHASEP(:,IB,IDR),EXTGRAD, DLEG,
     .     MAXNMICRO, DALB(IB,IDR), ALBEDOP(IB,IPA),
     .     EXTINCTP(IB,IPA), DPHASEWTP(:,IB,IDR), ML)

          RAYGRAD(:,IB,IDR) = RAYGRAD(:,IB,IDR) -
     .      EXTGRAD*DPATH(II)*ABSCELL*TRANSMIT*
     .      INPUTWEIGHT(:)
        ENDDO
        II = II + 1
      ENDDO
      RETURN
      END

      SUBROUTINE COMPUTE_SOURCE_DIRECTION(LEGEN, SOURCET,
     .  NLM, LOFJ, RIS, RNS, NSTOKES, RADIANCE, YLMDIR,
     .  NSTLEG, DELTAM, SRCTYPE, ML, MM, NS, YLMSUN,
     .  DIRFLUX, SECMU0, NLEG)
      IMPLICIT NONE
C     Compute the convolution of the RADIANCE with a set of
C     Legendre/Wigner coefficients (LEGEN) and evaluates it in a
C     particular direction (with weights in YLMDIR).
C     Outputs in SOURCET. Also includes the
C     single scatter if delta-M is not used.
C     This is used in the gradient computation where LEGEN
C     may either be expansion coefficients of either the
C     phase function OR its gradient with respect to the
C     unknowns.
      INTEGER NLM, RIS, RNS, NSTOKES, NSTLEG, ML, MM, NS
Cf2py intent(in) :: NLM, RIS, RNS, NSTOKES, NSTLEG, ML,MM,NS
      REAL YLMDIR(NSTLEG, NLM), LEGEN(NSTLEG, 0:NLEG)
Cf2py intent(in) :: YLMDIR, LEGEN
      INTEGER LOFJ(NLM), NLEG
Cf2py intent(in) :: LOFJ, NLEG
      REAL RADIANCE(NSTOKES,*), SOURCET(NSTOKES)
Cf2py intent(in) :: RADIANCE, SOURCET
      LOGICAL DELTAM
      CHARACTER SRCTYPE*1
      REAL YLMSUN(NSTLEG, NLM), DIRFLUX, SECMU0

      INTEGER L, ME, MS, M, JT, J
      REAL A1, B1, DA

       DO J = 1, RNS
         SOURCET(1) = SOURCET(1)
     .      + LEGEN(1,LOFJ(J))*RADIANCE(1,RIS+J)*YLMDIR(1,J)
       ENDDO

       IF (NSTOKES .GT. 1) THEN
         DO J = 1, RNS
           SOURCET(1) = SOURCET(1)
     .       + LEGEN(5,LOFJ(J))*RADIANCE(2,RIS+J)*YLMDIR(1,J)
         ENDDO
         DO J = 5, RNS
           SOURCET(2) = SOURCET(2)
     .       + LEGEN(5,LOFJ(J))*RADIANCE(1,RIS+J)*YLMDIR(2,J)
     .       + LEGEN(2,LOFJ(J))*RADIANCE(2,RIS+J)*YLMDIR(2,J)
     .       + LEGEN(3,LOFJ(J))*RADIANCE(3,RIS+J)*YLMDIR(5,J)
           SOURCET(3) = SOURCET(3)
     .      + LEGEN(5,LOFJ(J))*RADIANCE(1,RIS+J)*YLMDIR(6,J)
     .      + LEGEN(2,LOFJ(J))*RADIANCE(2,RIS+J)*YLMDIR(6,J)
     .      + LEGEN(3,LOFJ(J))*RADIANCE(3,RIS+J)*YLMDIR(3,J)
         ENDDO
       ENDIF
       IF (NSTOKES .EQ. 4) THEN
         DO J = 1, RNS
           SOURCET(2) = SOURCET(2)
     .      + LEGEN(6,LOFJ(J))*RADIANCE(4,RIS+J)*YLMDIR(5,J)
           SOURCET(3) = SOURCET(3)
     .      + LEGEN(6,LOFJ(J))*RADIANCE(4,RIS+J)*YLMDIR(3,J)
           SOURCET(4) = SOURCET(4)
     .        - LEGEN(6,LOFJ(J))*RADIANCE(3,RIS+J)*YLMDIR(4,J)
     .        + LEGEN(4,LOFJ(J))*RADIANCE(4,RIS+J)*YLMDIR(4,J)
         ENDDO
       ENDIF

       IF ((.NOT. DELTAM) .AND.
     .    (SRCTYPE .EQ. 'S') .OR. SRCTYPE .EQ. 'B') THEN
        J = 1
        DA = DIRFLUX*SECMU0
        DO L = 0, ML
          ME = MIN(L,MM)
          MS = -ME
          A1 = DA*LEGEN(1,L)
          B1 = DA*LEGEN(5,L)
          IF (J .LE. NS) THEN
            JT = J
            DO M = MS, ME
              SOURCET(1) = SOURCET(1) + A1*YLMDIR(1,J)*YLMSUN(1,J)
              J = J + 1
            ENDDO
            IF (NSTOKES .GT. 1) THEN
              J = JT
              DO M = MS, ME
                SOURCET(2)=SOURCET(2) + B1*YLMDIR(2,J)*YLMSUN(1,J)
                SOURCET(3)=SOURCET(3) + B1*YLMDIR(6,J)*YLMSUN(1,J)
                J = J + 1
              ENDDO
            ENDIF
          ENDIF
        ENDDO

       ENDIF

       RETURN
       END


      SUBROUTINE PREPARE_DERIV_INTERPS(GRIDPOS, NPTS,
     .    NPX, NPY, NPZ, MAXPG, DELX, DELY, XSTART, YSTART,
     .    ZLEVELS, OPTINTERPWT, INTERPPTR, IERR, ERRMSG)
C     This subroutine just precomputes the interpolation weights
C     (OPTINTERPWT) and pointers (INTERPPTR) from the property
C     grid onto the RTE grid in SHDOM analagously to in
C     TRILIN_INTERP_PROP. This subroutine could be expanded
C     to include all pre-computed quantities for the gradient
C     calculation.
        IMPLICIT NONE
        INTEGER NPTS, MAXPG
Cf2py intent(in) :: NPTS, MAXPG
        INTEGER NPX, NPY, NPZ
Cf2py intent(in) :: NPX, NPY, NPZ
        REAL GRIDPOS(3, NPTS)
Cf2py intent(in) :: GRIDPOS
        REAL OPTINTERPWT(8, NPTS)
Cf2py intent(out) :: OPTINTERPWT
        INTEGER INTERPPTR(8,NPTS)
Cf2py intent(out) :: INTERPPTR
        REAL DELX, DELY, XSTART, YSTART
Cf2py intent(in) :: DELX, DELY, XSTART, YSTART
        REAL ZLEVELS(*)
Cf2py intent(in) :: ZLEVELS
        INTEGER IERR
        CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG

        INTEGER IP

        IERR = 0
        DO IP=1,NPTS
          CALL GET_DERIV_INTERP(GRIDPOS(1,IP), GRIDPOS(2,IP),
     .       GRIDPOS(3,IP), NPX, NPY, NPZ, DELX,
     .       DELY, XSTART, YSTART, ZLEVELS,
     .       IERR, ERRMSG,
     .       OPTINTERPWT(:,IP), INTERPPTR(:,IP))
          IF (IERR .NE. 0) RETURN
        ENDDO
      RETURN
      END

      SUBROUTINE GET_DERIV_INTERP (X, Y, Z,
     .      NPX, NPY, NPZ, DELX, DELY,
     .      XSTART, YSTART, ZLEVELS,
     .      IERR, ERRMSG,
     .      OPTINTERPWT, INTERPPTR)
C     Modified from TRILIN_INTERP_PROP to only return the interpolation
C     weights (OPTINTERPWT) and pointers (INTERPPTR)from the property grid
C     onto the RTE grid for use in gradient calculations.
      IMPLICIT NONE
      INTEGER INTERPPTR(8)
      REAL OPTINTERPWT(8)
      REAL    X, Y, Z
      INTEGER IX, IXP, IY, IYP, IZ, L, IL, IM, IU, J
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I
      DOUBLE PRECISION U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
      INTEGER IERR
      CHARACTER ERRMSG*600

      INTEGER NPX, NPY, NPZ
      INTEGER NUMPHASE
      REAL DELX, DELY, XSTART, YSTART, INTERPTEMP
      REAL ZLEVELS(*)
C
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
        IERR = 1
        WRITE (ERRMSG,*) 'TRILIN: Beyond X domain',IX,NPX,X,XSTART
        RETURN
      ENDIF
      IXP = MOD(IX,NPX) + 1
      U = DBLE(X-XSTART-DELX*(IX-1))/DELX
      U = MAX( MIN( U, 1.0D0), 0.0D0)
      IF (U .LT. 1.0D-5) U = 0.0D0
      IF (U .GT. 1.0D0-1.0D-5) U = 1.0D0
      IY = INT((Y-YSTART)/DELY) + 1
      IF (ABS(Y-YSTART-NPY*DELY) .LT. 0.01*DELY) IY = NPY
      IF (IY .LT. 1 .OR. IY .GT. NPY) THEN
        IERR = 1
        WRITE (ERRMSG,*) 'TRILIN: Beyond Y domain',IY,NPY,Y,YSTART
        RETURN
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

      INTERPPTR(1) = I1
      INTERPPTR(2) = I2
      INTERPPTR(3) = I3
      INTERPPTR(4) = I4
      INTERPPTR(5) = I5
      INTERPPTR(6) = I6
      INTERPPTR(7) = I7
      INTERPPTR(8) = I8

      OPTINTERPWT(1) = F1
      OPTINTERPWT(2) = F2
      OPTINTERPWT(3) = F3
      OPTINTERPWT(4) = F4
      OPTINTERPWT(5) = F5
      OPTINTERPWT(6) = F6
      OPTINTERPWT(7) = F7
      OPTINTERPWT(8) = F8

C         Trilinearly interpolate the extinction, scattering
C     purely so we can get the phase interpolation weights.
C
C      EXTINCT = F1*EXTINCTP(I1) + F2*EXTINCTP(I2) + F3*EXTINCTP(I3)
C     .          + F4*EXTINCTP(I4) + F5*EXTINCTP(I5) + F6*EXTINCTP(I6)
C     .          + F7*EXTINCTP(I7) + F8*EXTINCTP(I8)
C      SCAT1 = F1*EXTINCTP(I1)*ALBEDOP(I1)
C      SCAT2 = F2*EXTINCTP(I2)*ALBEDOP(I2)
C      SCAT3 = F3*EXTINCTP(I3)*ALBEDOP(I3)
C      SCAT4 = F4*EXTINCTP(I4)*ALBEDOP(I4)
C      SCAT5 = F5*EXTINCTP(I5)*ALBEDOP(I5)
C      SCAT6 = F6*EXTINCTP(I6)*ALBEDOP(I6)
C      SCAT7 = F7*EXTINCTP(I7)*ALBEDOP(I7)
C      SCAT8 = F8*EXTINCTP(I8)*ALBEDOP(I8)
C      SCATTER = SCAT1+SCAT2+SCAT3+SCAT4+SCAT5+SCAT6+SCAT7+SCAT8
C
C      IF (EXTINCT .GT. EXTMIN) THEN
C        ALBEDO = SCATTER/EXTINCT
C      ELSE
C        ALBEDO = SCATTER/EXTMIN
C      ENDIF
C
C      IPHASE(1:MAXNMICRO) = IPHASEP(:,I1)
C      IPHASE(MAXNMICRO+1:2*MAXNMICRO) = IPHASEP(:,I2)
C      IPHASE(2*MAXNMICRO+1:3*MAXNMICRO) = IPHASEP(:,I3)
C      IPHASE(3*MAXNMICRO+1:4*MAXNMICRO) = IPHASEP(:,I4)
C      IPHASE(4*MAXNMICRO+1:5*MAXNMICRO) = IPHASEP(:,I5)
C      IPHASE(5*MAXNMICRO+1:6*MAXNMICRO) = IPHASEP(:,I6)
C      IPHASE(6*MAXNMICRO+1:7*MAXNMICRO) = IPHASEP(:,I7)
C      IPHASE(7*MAXNMICRO+1:) = IPHASEP(:,I8)
C
C      IF (INTERPMETHOD(2:2) .EQ. 'N') THEN
C        IF (SCATTER .GE. SCATMIN) THEN
C          PHASEINTERPWT(1:MAXNMICRO) = PHASEWTP(:,I1)*SCAT1/SCATTER
C          PHASEINTERPWT(MAXNMICRO+1:2*MAXNMICRO) = PHASEWTP(:,I2)
C     .      *SCAT2/SCATTER
C          PHASEINTERPWT(2*MAXNMICRO+1:3*MAXNMICRO) = PHASEWTP(:,I3)
C     .      *SCAT3*SCAT2/SCATTER
C          PHASEINTERPWT(3*MAXNMICRO+1:4*MAXNMICRO) = PHASEWTP(:,I4)
C     .      *SCAT4/SCATTER
C          PHASEINTERPWT(4*MAXNMICRO+1:5*MAXNMICRO) = PHASEWTP(:,I5)
C     .      *SCAT5/SCATTER
C          PHASEINTERPWT(5*MAXNMICRO+1:6*MAXNMICRO) = PHASEWTP(:,I6)
C     .      *SCAT6/SCATTER
C          PHASEINTERPWT(6*MAXNMICRO+1:7*MAXNMICRO) = PHASEWTP(:,I7)
C     .      *SCAT7/SCATTER
C          PHASEINTERPWT(7*MAXNMICRO+1:) = PHASEWTP(:,I8)
C     .      *SCAT8/SCATTER
C        ELSE
C          PHASEINTERPWT(:MAXNMICRO) = PHASEWTP(:,I1)*SCAT1/SCATMIN
C          PHASEINTERPWT(MAXNMICRO+1:2*MAXNMICRO) = PHASEWTP(:,I2)
C     .      *SCAT2/SCATMIN
C          PHASEINTERPWT(2*MAXNMICRO+1:3*MAXNMICRO) = PHASEWTP(:,I3)
C     .      *SCAT3/SCATMIN
C          PHASEINTERPWT(3*MAXNMICRO+1:4*MAXNMICRO) = PHASEWTP(:,I4)
C     .      *SCAT4/SCATMIN
C          PHASEINTERPWT(4*MAXNMICRO+1:5*MAXNMICRO) = PHASEWTP(:,I5)
C     .      *SCAT5/SCATMIN
C          PHASEINTERPWT(5*MAXNMICRO+1:6*MAXNMICRO) = PHASEWTP(:,I6)
C     .      *SCAT6/SCATMIN
C          PHASEINTERPWT(6*MAXNMICRO+1:7*MAXNMICRO) = PHASEWTP(:,I7)
C     .      *SCAT7/SCATMIN
C          PHASEINTERPWT(7*MAXNMICRO+1:) = PHASEWTP(:,I8)
C     .      *SCAT8/SCATMIN
C        ENDIF

C     Unlike TRILIN_INTERP_PROP, we don't consolidate weights here so we can
C     use these for derivatives more easily.

C      ELSEIF (INTERPMETHOD(2:2) .EQ. 'O') THEN
C        MAXSCAT = -1.0
C        PHASEINTERPWT(:) = 0.0
C        PHASEINTERPWT(1) = 1.0
C        IF (SCAT1 .GT. MAXSCAT .OR. ABS(F1-1) .LT. 0.001) THEN
C          CALL SSORT(PHASEWTP(:,I1), IPHASEP(:,I1), MAXNMICRO, -2)
C          MAXSCAT = SCAT1
C          IPHASE(1) = IPHASEP(1,I1)
C        ENDIF
C        IF (SCAT2 .GT. MAXSCAT .OR. ABS(F2-1) .LT. 0.001) THEN
C          CALL SSORT(PHASEWTP(:,I2), IPHASEP(:,I2), MAXNMICRO, -2)
C          MAXSCAT = SCAT2
C          IPHASE(1) = IPHASEP(1,I2)
C        ENDIF
C        IF (SCAT3 .GT. MAXSCAT .OR. ABS(F3-1) .LT. 0.001) THEN
C          CALL SSORT(PHASEWTP(:,I3), IPHASEP(:,I3), MAXNMICRO, -2)
C          MAXSCAT = SCAT3
C          IPHASE(1) = IPHASEP(1,I3)
C        ENDIF
C        IF (SCAT4 .GT. MAXSCAT .OR. ABS(F4-1) .LT. 0.001) THEN
C          CALL SSORT(PHASEWTP(:,I4), IPHASEP(:,I4), MAXNMICRO, -2)
C          MAXSCAT = SCAT4
C          IPHASE(1) = IPHASEP(1,I4)
C        ENDIF
C        IF (SCAT5 .GT. MAXSCAT .OR. ABS(F5-1) .LT. 0.001) THEN
C          CALL SSORT(PHASEWTP(:,I5), IPHASEP(:,I5), MAXNMICRO, -2)
C          MAXSCAT = SCAT5
C          IPHASE(1) = IPHASEP(1,I5)
C        ENDIF
C        IF (SCAT6 .GT. MAXSCAT .OR. ABS(F6-1) .LT. 0.001) THEN
C          CALL SSORT(PHASEWTP(:,I6), IPHASEP(:,I6), MAXNMICRO, -2)
C          MAXSCAT = SCAT6
C          IPHASE(1) = IPHASEP(1,I6)
C        ENDIF
C        IF (SCAT7 .GT. MAXSCAT .OR. ABS(F7-1) .LT. 0.001) THEN
C          CALL SSORT(PHASEWTP(:,I7), IPHASEP(:,I7), MAXNMICRO, -2)
C          MAXSCAT = SCAT7
C          IPHASE(1) = IPHASEP(1,I7)
C        ENDIF
C        IF (SCAT8 .GT. MAXSCAT .OR. ABS(F8-1) .LT. 0.001) THEN
C          CALL SSORT(PHASEWTP(:,I8), IPHASEP(:,I8), MAXNMICRO, -2)
C          MAXSCAT = SCAT8
C          IPHASE(1) = IPHASEP(1,I8)
C        ENDIF
C      ENDIF

      RETURN
      END



      SUBROUTINE PLANCK_DERIVATIVE (TEMP, UNITS, WAVENO, WAVELEN,
     .                              PLANCK)
C        Calculates the derivative of the Planck blackbody radiance with
C      respect to temperature. Derivative is returned in PLANCK.
C      -JRLoveridge 2021/04/21.
C     If UNITS='T' then
C     using brightness temperature units and the temperature is simply
C     returned. If UNITS='B' then doing a band integration and the
C     Planck blackbody radiance in [Watts /(meter^2 ster)] over a
C     wavenumber range [cm^-1] is returned. Otherwise, the Planck
C     blackbody radiance in [Watts /(meter^2 ster micron)] for a
C     temperature in [Kelvins] at a wavelength in [microns] is returned.
      IMPLICIT NONE
      REAL  TEMP, WAVENO(2), WAVELEN, PLANCK
      CHARACTER*1  UNITS
Cf2py intent(in) :: TEMP, WAVENO, WAVELEN, UNITS
Cf2py intent(out) :: PLANCK
      DOUBLE PRECISION X1, X2, F, P1, P2, T

      IF (UNITS .EQ. 'T') THEN
        PLANCK = 1
      ELSE IF (UNITS .EQ. 'B') THEN
        IF (TEMP .GT. 0.0) THEN
C         Use central difference for derivative.
          T = TEMP - 1.0D-4
          X1 = 1.4388D0*WAVENO(1)/T
          X2 = 1.4388D0*WAVENO(2)/T
          CALL INTEGRATE_PLANCK (X1, X2, F)
          P1= 1.1911D-8*(T/1.4388D0)**4 *F

          T = TEMP + 1.0D-4
          X1 = 1.4388D0*WAVENO(1)/T
          X2 = 1.4388D0*WAVENO(2)/T
          CALL INTEGRATE_PLANCK (X1, X2, F)
          P2= 1.1911D-8*(T/1.4388D0)**4 *F
          PLANCK = (P2 - P1) / (2.0D-4)
        ELSE
          PLANCK = 0.0
        ENDIF
      ELSE
        IF (TEMP .GT. 0.0) THEN
          PLANCK = (1.1911E8 *1.4388E4) / WAVELEN**6 *
     .      EXP(1.4388E4/(WAVELEN*TEMP)) /
     .      (TEMP*TEMP*(EXP(1.4388E4/(WAVELEN*TEMP)) - 1)**2)

        ELSE
          PLANCK = 0.0
        ENDIF
      ENDIF
      RETURN
      END
