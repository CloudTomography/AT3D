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
     .                   PHASEMAX, MAXNMICRO, TAUTOL,NOSURFACE,
     .                  CORRECTINTERPOLATE, TRANSCUT, SINGLESCATTER,
     .                  SFCGRIDRAD, NANG)
C    Calculates the Stokes Vector at the given directions (CAMMU, CAMPHI)
C    and positions CAMX,CAMY,CAMZ by integrating the source function.

Cf2py threadsafe
      IMPLICIT NONE
      INTEGER NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
      LOGICAL CORRECTINTERPOLATE, SINGLESCATTER, NOSURFACE
Cf2py intent(in) :: CORRECTINTERPOLATE, SINGLESCATTER, NOSURFACE
      DOUBLE PRECISION TRANSCUT
Cf2py intent(in) :: TRANSCUT
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
      REAL    GNDTEMP, GNDALBEDO
      REAL    SKYRAD(NSTOKES,NMU/2,NPHI0MAX)
      REAL    WAVENO(2), WAVELEN
Cf2py intent(in) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN
      REAL    MU(*), PHI(NMU,*), WTDO(NMU,*)
Cf2py intent(in) :: MU, PHI, WTDO
      REAL    XGRID(*), YGRID(*), ZGRID(*), GRIDPOS(3,*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID, GRIDPOS
      REAL    SFCGRIDPARMS(*), BCRAD(NSTOKES, *)
Cf2py intent(in) :: SFCGRIDPARMS
Cf2py intent(in, out) :: BCRAD
      REAL    EXTINCT(NPTS,NPART), ALBEDO(NPTS,NPART)
      REAL    TOTAL_EXT(NPTS), LEGEN(NSTLEG,0:NLEG,*)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN, TOTAL_EXT

      REAL    SFCGRIDRAD(NANG/2 + 1, *)
Cf2py intent(in) :: SFCGRIDRAD
      INTEGER NANG
Cf2py intent(in) :: NANG

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
      DOUBLE PRECISION TAUTOL
Cf2py intent(in) :: TAUTOL
      INTEGER I, J, L, SIDE
      INTEGER IVIS, ITOP
      LOGICAL VALIDRAD
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION U, R, PI
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT, VISRAD(NSTOKES)

      REAL TIME1, TIME2, TIME_SOURCE, TIME_RAY_TOTAL
      REAL  MEAN, STD1, STD2

      TIME_RAY_TOTAL = 0.0
      TIME_SOURCE = 0.0
      IERR = 0
C         Make the isotropic radiances for the top boundary
C      CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN,
C     .                            UNITS, NTOPPTS, NSTOKES, BCRAD(1,1))


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

        IF (MURAY .GT. 0.0) THEN
          CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD,WAVENO,WAVELEN,
     .                              UNITS, NTOPPTS, NSTOKES, BCRAD(1,1),
     .                              NPHI0MAX, NMU, 1, 1,
     .                              REAL(MU2), REAL(PHI2), MU, PHI, 
     .                              NPHI0,1)
        ELSE
          DO ITOP=1,NTOPPTS
            BCRAD(:,ITOP) = 0.0
          ENDDO
        ENDIF


C         Integrate the extinction and source function along this ray
C         to calculate the Stokes radiance vector for this pixel
        TRANSMIT = 1.0D0 ; VISRAD(:) = 0.0D0
C        CALL CPU_TIME(TIME1)
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
     .                     PHASEINTERPWT, PHASEMAX,MAXNMICRO, TAUTOL,
     .                     TIME_SOURCE, CORRECTINTERPOLATE,
     .                     TRANSCUT, SINGLESCATTER, NOSURFACE,
     .                     SFCGRIDRAD, NANG, UNITS, WAVENO)
      IF (IERR .NE. 0) RETURN
C       CALL CPU_TIME(TIME2)
C       TIME_RAY_TOTAL = TIME_RAY_TOTAL + TIME2-TIME1
900   CONTINUE


        STOKES(:, IVIS) = VISRAD(:)
      ENDDO
C      PRINT *, 'NRAYS', NPIX
C      PRINT *, 'TIME_RAY_TOTAL', TIME_RAY_TOTAL
C      PRINT *, 'TIME_SOURCE', TIME_SOURCE
      RETURN
      END

      SUBROUTINE LEVISAPPROX_GRADIENT(NSTOKES, NX, NY, NZ, NPTS,
     .           NCELLS,MAXPG,
     .           ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,
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
     .           TAUTOL, SINGLESCATTER, IERR, ERRMSG, INTERPMETHOD,
     .           MAXNMICRO, DIPHASEP, IPHASEP, PHASEWTP, DPHASEWTP,
     .           PHASEINTERPWT, ALBEDOP, EXTINCTP, OPTINTERPWT,
     .           INTERPPTR, EXTMIN, SCATMIN, DOEXACT, TEMP, PHASEMAX,
     .           DTEMP, DEXTM, DALBM, DFJ, MAXSUBGRIDINTS,
     .           TRANSCUT, LONGEST_PATH_PTS,DERIV_MAXNMICRO, 
     .            NSFCDER,NSFCPTS,SFCDER,DELXSFC,
     .             DELYSFC, NANG,
     .             SFCGRIDRAD)
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
      INTEGER MAXPG, NCELLS
Cf2py intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG,
Cf2py intent(in) :: NPTS, NCELLS, MAXPG
      INTEGER NPX, NPY, NPZ
      REAL    DELX, DELY, XSTART, YSTART
      REAL    ZLEVELS(*)
      REAL    EXTDIRP(*)
      DOUBLE PRECISION UNIFORMZLEV
Cf2py intent(in) :: DELX, DELY, XSTART, YSTART, NPX, NPY, NPZ, ZLEVELS, EXTDIRP, UNIFORMZLEV
      INTEGER ML, MM, NSTLEG, NLM, NLEG, NUMPHASE, NPART
Cf2py intent(in) :: ML, MM, NSTLEG, NLM, NLEG, NUMPHASE, NPART
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
      REAL    GNDTEMP, GNDALBEDO
      REAL    SKYRAD(NSTOKES,NMU/2,NPHI0MAX)
      REAL    WAVENO(2), WAVELEN
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
      INTEGER DERIV_MAXNMICRO
Cf2py intent(in) :: DERIV_MAXNMICRO
      INTEGER DIPHASEP(DERIV_MAXNMICRO,MAXPG,NUMDER)
      INTEGER IPHASEP(MAXNMICRO,MAXPG,NPART)
      REAL    PHASEWTP(MAXNMICRO,MAXPG,NPART)
      REAL    DPHASEWTP(DERIV_MAXNMICRO,MAXPG,NUMDER)
      REAL    DTEMP(MAXPG,NUMDER)
Cf2py intent(in) :: DIPHASEP, IPHASEP, PHASEWTP, DPHASEWTP,DTEMP
      REAL    OPTINTERPWT(8,NPTS)
      INTEGER INTERPPTR(8,NPTS)
Cf2py intent(in) :: OPTINTERPWT, INTERPPTR
      REAL DEXTM(MAXPG,NUMDER), DALBM(8,NPTS,NUMDER)
Cf2py intent(in) :: DEXTM, DALBM
      REAL DFJ(8,NPTS,NUMDER)
Cf2py intent(in) :: DFJ
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
      REAL DPATH(LONGEST_PATH_PTS,*)
      INTEGER :: NUNCERTAINTY
Cf2py intent(in) :: NUNCERTAINTY
      DOUBLE PRECISION UNCERTAINTIES(NUNCERTAINTY,NUNCERTAINTY,*)
      INTEGER DPTR(LONGEST_PATH_PTS,*)
Cf2py intent(in) :: DPATH, DPTR, UNCERTAINTIES
      REAL JACOBIAN(NSTOKES,NUMDER,NUM_JACOBIAN_PTS,*)
Cf2py intent(in,out) :: JACOBIAN
      LOGICAL MAKEJACOBIAN
Cf2py intent(in) :: MAKEJACOBIAN
      INTEGER JI,NUM_JACOBIAN_PTS, JACOBIANPTR(NUM_JACOBIAN_PTS)
Cf2py intent(in) :: NUM_JACOBIAN_PTS, JACOBIANPTR
      INTEGER RAYS_PER_PIXEL(*)
Cf2py intent(in) :: RAYS_PER_PIXEL
      DOUBLE PRECISION   RAY_WEIGHTS(*), STOKES_WEIGHTS(NSTOKES, *)
Cf2py intent(in) :: RAY_WEIGHTS, STOKES_WEIGHTS
      DOUBLE PRECISION TAUTOL
Cf2py intent(in) :: TAUTOL
      DOUBLE PRECISION TRANSCUT
Cf2py intent(in) :: TRANSCUT
      LOGICAL SINGLESCATTER
Cf2py intent(in) :: SINGLESCATTER
      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG
      INTEGER MAXSUBGRIDINTS, LONGEST_PATH_PTS
Cf2py intent(in) :: MAXSUBGRIDINTS, LONGEST_PATH_PTS
      INTEGER NSFCDER, NSFCPTS, NXSFC, NYSFC, NANG
Cf2py intent(in) :: NSFCDER, NSFCPTS, NXSFC, NYSFC, NANG
      INTEGER SFCDER(NSFCDER, NSFCPTS)
Cf2py intent(in) :: SFCDER
      REAL SFCGRIDRAD(NANG/2 + 1,*)
      REAL DELXSFC, DELYSFC
Cf2py intent(int) :: SFCGRIDRAD, DELXSFC, DELYSFC     


C     At the moment these surface variables do nothing.
      DOUBLE PRECISION SFCGRAD_RAY(NSTOKES,NSFCDER, NSFCPTS)
      DOUBLE PRECISION SFCSRCGRAD_RAY(NSTOKES,NANG/2 + 1, NSFCPTS)

C     Local variables
      DOUBLE PRECISION WEIGHT
      DOUBLE PRECISION VISRAD(NSTOKES)
      DOUBLE PRECISION RAYGRAD(NSTOKES,MAXPG,NUMDER)
      DOUBLE PRECISION RAYGRAD_PIXEL(NSTOKES,MAXPG,NUMDER)
      DOUBLE PRECISION, ALLOCATABLE :: ADJ_WEIGHTS(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: BEAM_WEIGHT(:)
      DOUBLE PRECISION WEIGHT_VEC(NSTOKES)
      INTEGER IPIX, J, L, SIDE, IRAY, IP
      LOGICAL VALIDRAD, VERBOSE
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION R, PI
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT
      INTEGER M, ME, MS, NS, NS1
      INTEGER I2, ITOP
      REAL TIME_SOURCE(3), TIME_DIRECT_POINT, TIME_DIRECT_SURFACE
      REAL TIME_RADIANCE, TIME_RAY_TOTAL, TIME1, TIME2
      REAL TIME_SUBGRID, TIME_ALLOCATE
      LOGICAL CORRECTINTERP_LOCAL, NOSURFACE_LOCAL
      INTEGER, ALLOCATABLE :: LOFJ(:)
      ALLOCATE (LOFJ(NLM))
      ALLOCATE (ADJ_WEIGHTS(NSTOKES, NPIX))
      ALLOCATE (BEAM_WEIGHT(NPTS))
      IERR = 0
      GRADOUT = 0.0D0
      STOKESOUT = 0.0D0
      TIME_SOURCE = 0.0
      TIME_DIRECT_POINT = 0.0
      TIME_DIRECT_SURFACE = 0.0
      TIME_RADIANCE = 0.0
      TIME_RAY_TOTAL = 0.0
      TIME_SUBGRID = 0.0
      TIME_ALLOCATE = 0.0
      CORRECTINTERP_LOCAL = .TRUE.
      NOSURFACE_LOCAL = .FALSE.

C      CALL CPU_TIME(TIME1)
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
C      CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD, WAVENO, WAVELEN,
C     .                            UNITS, NTOPPTS, NSTOKES, BCRAD(1,1))
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

      IF (MAKEJACOBIAN) THEN
C     ===== JACOBIAN PATH: Use original GRAD_INTEGRATE_1RAY =====
C     The Jacobian requires per-pixel RAYGRAD_PIXEL, which is
C     incompatible with the adjoint accumulation approach.
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
          VERBOSE = .FALSE.
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
          IF (MURAY .GT. 0.0) THEN
            CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD,WAVENO,WAVELEN,
     .                              UNITS, NTOPPTS, NSTOKES, BCRAD(1,1),
     .                              NPHI0MAX, NMU, 1, 1,
     .                              MU2, PHI2, MU, PHI, NPHI0,
     .                              1)
          ELSE
            DO ITOP=1,NTOPPTS
              BCRAD(:,ITOP) = 0.0
            ENDDO
          ENDIF
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
     .             VISRAD, VALIDRAD, TOTAL_EXT, NPART, RAYGRAD,
     .             RSHPTR, RADIANCE, LOFJ, PARTDER, NUMDER, DEXT,
     .             DALB, DLEG, MAXPG, DNUMPHASE, SOLARFLUX,
     .             NPX, NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER, TAUTOL, SINGLESCATTER,
     .             GNDALBEDO, IERR, ERRMSG, TEMP, PHASEMAX,
     .             PHASEINTERPWT, OPTINTERPWT, INTERPPTR,
     .             MAXNMICRO, EXTINCTP, ALBEDOP, PHASEWTP, DPHASEWTP,
     .             IPHASEP, DIPHASEP, WAVENO, UNITS, EXTMIN, SCATMIN,
     .             DOEXACT, INTERPMETHOD, DTEMP, DEXTM, DALBM,
     .             DFJ, MAXSUBGRIDINTS, TRANSCUT,LONGEST_PATH_PTS,
     .             TIME_SOURCE,DERIV_MAXNMICRO,
     .             TIME_DIRECT_POINT, TIME_DIRECT_SURFACE,
     .             TIME_RADIANCE, TIME_SUBGRID, TIME_ALLOCATE,
     .             VERBOSE, NSFCDER,NSFCPTS,SFCDER,DELXSFC,
     .             DELYSFC, SFCGRAD_RAY, SFCSRCGRAD_RAY, NANG,
     .             SFCGRIDRAD)
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
        DO JI = 1,NUM_JACOBIAN_PTS
          JACOBIAN(:,:,JI,IPIX) =
     .      RAYGRAD_PIXEL(:,JACOBIANPTR(JI),:)
        END DO
      ENDDO

      ELSE
C     ===== ADJOINT PATH: Three-phase gradient computation =====

C     --- Phase 1: Forward radiance pass ---
      IRAY = 0
      DO IPIX = 1, NPIX
        DO I2=1 ,RAYS_PER_PIXEL(IPIX)
          IRAY = IRAY + 1
          X0 = CAMX(IRAY)
          Y0 = CAMY(IRAY)
          Z0 = CAMZ(IRAY)
          MU2 = CAMMU(IRAY)
          PHI2 = CAMPHI(IRAY)
          MURAY = -MU2
          PHIRAY = PHI2 - PI
          IF (Z0 .GT. ZGRID(NZ)) THEN
            IF (MURAY .GE. 0.0) THEN
              VISRAD(:) = 0.0
              GOTO 800
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
          IF (MURAY .GT. 0.0) THEN
            CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD,WAVENO,WAVELEN,
     .                              UNITS, NTOPPTS, NSTOKES, BCRAD(1,1),
     .                              NPHI0MAX, NMU, 1, 1,
     .                              MU2, PHI2, MU, PHI, NPHI0,
     .                              1)
          ELSE
            DO ITOP=1,NTOPPTS
              BCRAD(:,ITOP) = 0.0
            ENDDO
          ENDIF
          TRANSMIT = 1.0D0 ; VISRAD = 0.0D0
          CALL INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES, NSTLEG,
     .             NSTPHASE, NSCATANGLE, PHASETAB, NX, NY, NZ,
     .             NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR,
     .             CELLFLAGS, XGRID, YGRID, ZGRID, GRIDPOS,
     .             ML, MM, NLM, NLEG, NUMPHASE, NMU, NPHI0MAX,
     .             NPHI0, MU, PHI, WTDO, DELTAM, SRCTYPE,
     .             WAVELEN, SOLARMU, SOLARAZ, EXTINCT, ALBEDO,
     .             LEGEN, IPHASE, DIRFLUX, SHPTR, SOURCE, YLMSUN,
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .             MU2, PHI2, X0, Y0, Z0,
     .             XE, YE, ZE, SIDE, TRANSMIT, VISRAD, VALIDRAD,
     .             TOTAL_EXT, NPART, IERR, ERRMSG, INTERPMETHOD,
     .             PHASEINTERPWT, PHASEMAX, MAXNMICRO, TAUTOL,
     .             TIME_RADIANCE, CORRECTINTERP_LOCAL,
     .             TRANSCUT, SINGLESCATTER, NOSURFACE_LOCAL,
     .             SFCGRIDRAD, NANG, SKYRAD, UNITS, WAVENO)
          IF (IERR .NE. 0) RETURN
  800     CONTINUE
          DO NS=1,NSTOKES
            STOKESOUT(NS,IPIX) = STOKESOUT(NS,IPIX) + VISRAD(NS)*
     .        RAY_WEIGHTS(IRAY)*STOKES_WEIGHTS(NS,IPIX)
          ENDDO
        ENDDO
      ENDDO

C     --- Phase 2: Adjoint weights + cost ---
      COST = 0.0D0
      ADJ_WEIGHTS = 0.0D0
      DO IPIX = 1, NPIX
        CALL COMPUTE_ADJOINT_WEIGHTS(DBLE(STOKESOUT(:,IPIX)),
     .    DBLE(MEASUREMENTS(:,IPIX)), UNCERTAINTIES(:,:,IPIX),
     .    COSTFUNC, NSTOKES, NUNCERTAINTY,
     .    ADJ_WEIGHTS(:,IPIX), COST, NCOST, IERR, ERRMSG)
        IF (IERR .NE. 0) RETURN
      ENDDO

C     --- Phase 3: Adjoint derivative pass ---
      GRADOUT = 0.0D0
      BEAM_WEIGHT = 0.0D0
      IRAY = 0
      DO IPIX = 1, NPIX
        DO I2=1 ,RAYS_PER_PIXEL(IPIX)
          IRAY = IRAY + 1
          X0 = CAMX(IRAY)
          Y0 = CAMY(IRAY)
          Z0 = CAMZ(IRAY)
          MU2 = CAMMU(IRAY)
          PHI2 = CAMPHI(IRAY)
          MURAY = -MU2
          PHIRAY = PHI2 - PI
          VERBOSE = .FALSE.
          IF (Z0 .GT. ZGRID(NZ)) THEN
            IF (MURAY .GE. 0.0) THEN
              GOTO 850
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
          IF (MURAY .GT. 0.0) THEN
            CALL COMPUTE_TOP_RADIANCES (SRCTYPE, SKYRAD,WAVENO,WAVELEN,
     .                              UNITS, NTOPPTS, NSTOKES, BCRAD(1,1),
     .                              NPHI0MAX, NMU, 1, 1,
     .                              MU2, PHI2, MU, PHI, NPHI0,
     .                              1)
          ELSE
            DO ITOP=1,NTOPPTS
              BCRAD(:,ITOP) = 0.0
            ENDDO
          ENDIF
C         Compute per-ray adjoint weight vector
          DO NS=1,NSTOKES
            WEIGHT_VEC(NS) = ADJ_WEIGHTS(NS,IPIX)*
     .        RAY_WEIGHTS(IRAY)*STOKES_WEIGHTS(NS,IPIX)
          ENDDO
          TRANSMIT = 1.0D0
          CALL ADJOINT_INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES,
     .             NSTLEG,
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
     .             VALIDRAD, TOTAL_EXT, NPART,
     .             WEIGHT_VEC, GRADOUT, NGRAD, BEAM_WEIGHT,
     .             RSHPTR, RADIANCE, LOFJ, PARTDER, NUMDER, DEXT,
     .             DALB, DLEG, MAXPG, DNUMPHASE, SOLARFLUX,
     .             NPX, NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER, TAUTOL, SINGLESCATTER,
     .             GNDALBEDO, IERR, ERRMSG, TEMP, PHASEMAX,
     .             PHASEINTERPWT, OPTINTERPWT, INTERPPTR,
     .             MAXNMICRO, EXTINCTP, ALBEDOP, PHASEWTP, DPHASEWTP,
     .             IPHASEP, DIPHASEP, WAVENO, UNITS, EXTMIN, SCATMIN,
     .             DOEXACT, INTERPMETHOD, DTEMP, DEXTM, DALBM,
     .             DFJ, MAXSUBGRIDINTS, TRANSCUT,LONGEST_PATH_PTS,
     .             TIME_SOURCE,DERIV_MAXNMICRO,
     .             TIME_DIRECT_POINT, TIME_DIRECT_SURFACE,
     .             TIME_RADIANCE, TIME_SUBGRID, TIME_ALLOCATE,
     .             VERBOSE, NSFCDER,NSFCPTS,SFCDER,DELXSFC,
     .             DELYSFC, SFCGRAD_RAY, SFCSRCGRAD_RAY, NANG,
     .             SFCGRIDRAD)
          IF (IERR .NE. 0) RETURN
  850     CONTINUE
        ENDDO
      ENDDO

C     --- Phase 4: Batched direct beam derivatives ---
      IF (EXACT_SINGLE_SCATTER .AND. (SRCTYPE .EQ. 'S' .OR.
     .    SRCTYPE .EQ. 'B')) THEN
        DO IP = 1, NPTS
          IF (BEAM_WEIGHT(IP) .NE. 0.0D0) THEN
            CALL COMPUTE_DIRECT_BEAM_DERIV_ADJOINT(DPATH(:,IP),
     .        DPTR(:,IP), NUMDER, DEXTM, BEAM_WEIGHT(IP),
     .        MAXPG, GRADOUT, NGRAD, LONGEST_PATH_PTS)
          ENDIF
        ENDDO
      ENDIF

      ENDIF
C     ===== END OF JACOBIAN/ADJOINT BRANCH =====

      DEALLOCATE (LOFJ, ADJ_WEIGHTS, BEAM_WEIGHT)
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
     .		         RADOUT, VALIDRAD, TOTAL_EXT, NPART, RAYGRAD,
     .		         RSHPTR, RADIANCE, LOFJ, PARTDER, NUMDER, DEXT,
     .             DALB, DLEG, MAXPG, DNUMPHASE, SOLARFLUX,
     .             NPX, NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER, TAUTOL, SINGLESCATTER,
     .             GNDALBEDO, IERR, ERRMSG, TEMP, PHASEMAX,
     .             PHASEINTERPWT, OPTINTERPWT, INTERPPTR,
     .             MAXNMICRO, EXTINCTP, ALBEDOP, PHASEWTP, DPHASEWTP,
     .             IPHASEP, DIPHASEP, WAVENO, UNITS, EXTMIN, SCATMIN,
     .             DOEXACT, INTERPMETHOD, DTEMP, DEXTM,DALBM,
     .             DFJ, MAXSUBGRIDINTS, TRANSCUT, LONGEST_PATH_PTS,
     .             TIME_SOURCE, DERIV_MAXNMICRO,
     .             TIME_DIRECT_POINT, TIME_DIRECT_SURFACE,
     .             TIME_RADIANCE, TIME_SUBGRID, TIME_ALLOCATE,
     .          VERBOSE, NSFCDER,NSFCPTS,SFCDER,DELXSFC,
     .             DELYSFC, SFCGRAD_RAY, SFCSRCGRAD_RAY, NANG,
     .             SFCGRIDRAD)
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
      LOGICAL EXACT_SINGLE_SCATTER, VERBOSE
      INTEGER NPX, NPY, NPZ, MAXPG, BCELL, NANG
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
      INTEGER   DERIV_MAXNMICRO
      REAL      EXTINCTP(MAXPG,NPART), ALBEDOP(MAXPG,NPART)
      REAL      PHASEWTP(MAXNMICRO,MAXPG,NPART)
      REAL      DPHASEWTP(DERIV_MAXNMICRO,MAXPG,NUMDER)
      INTEGER   IPHASEP(MAXNMICRO,MAXPG,NPART)
      INTEGER   DIPHASEP(DERIV_MAXNMICRO,MAXPG,NUMDER)
      REAL      WAVENO(2)
      DOUBLE PRECISION EXTMIN, SCATMIN
      INTEGER   DOEXACT(NUMDER)
      INTEGER   MAXSUBGRIDINTS

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
      INTEGER LONGEST_PATH_PTS
      REAL    DPATH(LONGEST_PATH_PTS,*), SECMU0
      INTEGER DPTR(LONGEST_PATH_PTS,*), N
      REAL DEXTM(MAXPG,NUMDER), DALBM(8,NPTS,NUMDER)
      REAL DFJ(8,NPTS,NUMDER)
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
      LOGICAL :: SINGLESCATTER
      DOUBLE PRECISION :: DIRRAD(NSTOKES,4), BOUNDINTERP(4)
      INTEGER :: BOUNDPTS(4), IB,IP
      REAL :: GNDALBEDO
      INTEGER IERR
      CHARACTER ERRMSG*600

      INTEGER NSFCDER, NSFCPTS, NXSFC, NYSFC
      INTEGER SFCDER(NSFCDER, NSFCPTS)
      DOUBLE PRECISION SFCGRAD_RAY(NSTOKES,NSFCDER, NSFCPTS)
      DOUBLE PRECISION SFCSRCGRAD_RAY(NSTOKES,NANG/2 + 1, NSFCPTS)
      REAL SFCGRIDRAD(NANG/2 + 1,*)
      REAL DELXSFC, DELYSFC      

      REAL TIME_SOURCE(3), TIME_DIRECT_POINT, TIME_DIRECT_SURFACE
      REAL TIME_RADIANCE, TIME_SUBGRID, TIME_ALLOCATE

      INTEGER, ALLOCATABLE :: PASSEDPOINTS(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDRAD(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDABSCELL(:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDTRANSMIT(:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDINTERP0(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDINTERP1(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDDELS(:)
      INTEGER NPASSED, K, MAXBYOPT
      REAL TIME1, TIME2
      LOGICAL DOSINGSCAT

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/
      DATA DONEFACE/0,0,0,0,0,0,0,0, 0,1,0,3,0,5,0,7, 2,0,4,0,6,0,8,0,
     .              0,0,1,2,0,0,5,6, 3,4,0,0,5,6,0,0,
     .              0,0,0,0,1,2,3,4, 5,6,7,8,0,0,0,0/

C         TRANSCUT is the transmission to stop the integration at
C      TRANSCUT = 5.0E-5
C0.0D0
C         TAUTOL is the maximum optical path for the subgrid intervals
C      it is now a numerical parameter.
C      TAUTOL = 0.2

      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)

      MAXSUBGRIDINTS = MAX(MAXCELLSCROSS, MAXSUBGRIDINTS)


      ALLOCATE (PASSEDPOINTS(8,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDINTERP0(8,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDINTERP1(8,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDDELS(MAXSUBGRIDINTS))
      ALLOCATE (PASSEDRAD(NSTOKES,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDABSCELL(MAXSUBGRIDINTS))
      ALLOCATE (PASSEDTRANSMIT(MAXSUBGRIDINTS))

C     initialize counter for subgrid intervals.
      NPASSED = 1

      PI = ACOS(-1.0D0)
C       Calculate the generalized spherical harmonics for this direction
      CALL YLMALL (.FALSE.,SNGL(MU2),SNGL(PHI2),ML,MM,NSTLEG, YLMDIR)
      SECMU0 = 1.0D0/ABS(SOLARMU)
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

      IFACE = 0
      NGRID = 0
C         Loop until reach a Z boundary or transmission is very small
      VALIDRAD = .FALSE.
C      CALL CPU_TIME(TIME2)
C      TIME_ALLOCATE = TIME_ALLOCATE + TIME2 - TIME1
      DO WHILE (.NOT. VALIDRAD .AND. ICELL .GT. 0)

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
C      CALL CPU_TIME(TIME1)

        CALL COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR,
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, PHASEINTERPWT,
     .             DIRFLUX,DERIV_MAXNMICRO,
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
     .             SINGLESCATTER, DTEMP, DEXTM, DALBM, DFJ, TIME_SOURCE)

C         CALL CPU_TIME(TIME2)
C         TIME_SOURCE = TIME_SOURCE + TIME2 - TIME1
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

        SRCSINGSCAT = 0.0

C           Loop over the subgrid cells
        DO IT = 1, NTAU

          DOSINGSCAT = .TRUE.
C         save the interpolation kernel and pointers for each subgrid
C         interval so that the radiance gradient can be calculated.
C          CALL CPU_TIME(TIME1)
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
C             multiplied by extinction derivative is calculated at the end of the
C             subroutine.
              SRCGRAD = ( 0.5*(GRAD0+GRAD1)
     .          + 0.08333333333*(EXT0*GRAD1-EXT1*GRAD0)*DELS
     .           *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT

               IF (EXACT_SINGLE_SCATTER .AND.
     .           SRCTYPE .NE. 'T') THEN
                 SRCSINGSCAT = SRCSINGSCAT + TRANSMIT*ABSCELL*
     .                    ( 0.5*(SINGSCAT0+SINGSCAT1)
     .            + 0.08333333333*(EXT0*SINGSCAT1-EXT1*SINGSCAT0)*DELS
     .                    *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
               ENDIF
            ENDIF

            RADOUT(:) = RADOUT(:) + TRANSMIT*SRC(:)*ABSCELL
C         The stuff below is only done when there is extinction. Strictly
C         speaking, it should only need to be performed when there is no
C         extinction from any of the particle species for which derivatives
C         are calculated for.

C         Store the radiance contribution for this subgrid interval.
C         The cumulative sum needed by COMPUTE_RADIANCE_DERIVATIVE
C         is built in a single O(N) backward pass after the ray loop.
            PASSEDABSCELL(NPASSED) = ABSCELL
            PASSEDTRANSMIT(NPASSED) = TRANSMIT
            PASSEDRAD(:,NPASSED) = TRANSMIT*SRC(:)*ABSCELL
C          CALL CPU_TIME(TIME2)
C          TIME_SUBGRID = TIME_SUBGRID + TIME2-TIME1

C          CALL CPU_TIME(TIME1)
            IF (.NOT. OUTOFDOMAIN) THEN
              DO KK=1,8
                IP = GRIDPTR(KK,ICELL)
                DO K=1,8
                  IB = INTERPPTR(K,IP)
                  RAYGRAD(:,IB,:) = RAYGRAD(:,IB,:) +
     .            TRANSMIT*SRCGRAD(:,K,KK,:)*ABSCELL
                ENDDO
C                IF (EXACT_SINGLE_SCATTER .AND. (SRCTYPE .EQ. 'S'
C     .            .OR. SRCTYPE .EQ. 'B')) THEN
CC             Add gradient component due to the direct solar beam.
CC             Sensitivity of single scattered radiation to
CC             extinction along the path between the gridpoint and the sun.
C                  CALL COMPUTE_DIRECT_BEAM_DERIV(DPATH(:,IP),
C     .            DPTR(:,IP),
C     .            NUMDER, DEXTM, TRANSMIT,
C     .            ABSCELL, SRCSINGSCAT, NSTOKES,
C     .            MAXPG,RAYGRAD, LONGEST_PATH_PTS)
C
C                ENDIF
              ENDDO
            ENDIF
C          CALL CPU_TIME(TIME2)
C          TIME_DIRECT_POINT= TIME_DIRECT_POINT +
C     .      TIME2 - TIME1

            NPASSED = NPASSED + 1
            IF (NPASSED .GT. MAXSUBGRIDINTS) THEN
              IERR = 1
              WRITE(ERRMSG, *) "GRAD_INTEGRATE_1RAY: The maximum ",
     .        "number of ",
     .        " subgrid intervals for calculation of the",
     .        " radiance along the",
     .        " ray path has been exceeded.", "NPASSED=", NPASSED,
     .        "MAXSUBGRIDINTS=", MAXSUBGRIDINTS
              RETURN
            ENDIF

          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
            SRC(:) = 0.0
            SRCGRAD = 0.0
            SRCSINGSCAT = 0.0

          ENDIF

          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          SRCEXT1(:) = SRCEXT0(:)
          GRAD1 = GRAD0
          SINGSCAT1 = SINGSCAT0
C                End of sub grid cell loop
        ENDDO

C       C             Add gradient component due to the direct solar beam.
C             Sensitivity of single scattered radiation to
C             extinction along the path between the gridpoint and the sun.
        IF (EXACT_SINGLE_SCATTER .AND. (SRCTYPE .EQ. 'S'
     .       .OR. SRCTYPE .EQ. 'B')) THEN
          DO KK =1,8
            IP = GRIDPTR(KK,ICELL)
                  CALL COMPUTE_DIRECT_BEAM_DERIV(DPATH(:,IP),
     .            DPTR(:,IP),
     .            NUMDER, DEXTM, 1.0d0,
     .            1.0d0, SRCSINGSCAT(:,KK), NSTOKES,
     .            MAXPG,RAYGRAD, LONGEST_PATH_PTS)

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
C
C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
C        .OR. NGRID.GT.MAXCELLSCROSS Exceeding MAXCELLSCROSS is no
C        longer a valid reason for stopping ray integration as this
C        may compromise physical accuracy.
        IF ((TRANSMIT .LT. TRANSCUT)) THEN
          VALIDRAD = .TRUE.
C         No boundary radiance; initialise the sentinel entry.
          PASSEDRAD(:,NPASSED) = 0.0D0
          PASSEDTRANSMIT(NPASSED) = 1.0D0

        ELSE IF (INEXTCELL .EQ. 0 .AND. IFACE .GE. 5) THEN
          VALIDRAD = .TRUE.
C          CALL CPU_TIME(TIME1)
          CALL FIND_BOUNDARY_RADIANCE_GRAD (NSTOKES, XN, YN,
     .                      SNGL(MU2), SNGL(PHI2),
     .                      IC, KFACE, GRIDPTR, GRIDPOS,
     .                      MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .                      NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .                      SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, DIRFLUX,
     .                      SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                      RADBND, SFCGRIDRAD, NANG,
     .                      UNITS, WAVENO, BOUNDPTS, BOUNDINTERP, 
     .                      DIRRAD,
     .                      GNDALBEDO,SFCGRAD_RAY,NSFCDER,SFCDER,
     .                      TRANSMIT, SFCSRCGRAD_RAY,
     .                      DELXSFC,DELYSFC,DTEMP,TEMP)
          RADOUT(:) = RADOUT(:) + TRANSMIT*RADBND(:)
C         These values aren't actually used. set to -1.0 for
C         debugging purposes.
          PASSEDTRANSMIT(NPASSED) = TRANSMIT
          PASSEDABSCELL(NPASSED) = -1.0
C         Store boundary radiance contribution in the last entry.
C         The backward cumsum after the ray loop will propagate it.
          PASSEDRAD(:,NPASSED) = TRANSMIT*RADBND(:)
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
     .            DPTR(:,IP),
     .            NUMDER, DEXTM, TRANSMIT,
     .            1.0D0, SNGL(BOUNDINTERP(KK)*DIRRAD(:,KK)),
     .            NSTOKES, MAXPG,RAYGRAD,
     .            LONGEST_PATH_PTS)
            ENDDO
          ENDIF
C          CALL CPU_TIME(TIME2)
C          TIME_DIRECT_SURFACE = TIME_DIRECT_SURFACE +
C     .      TIME2 - TIME1
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
C     Build the cumulative future-radiance array via backward sum.
C     After this loop PASSEDRAD(:,k) = (1/T_k) * sum_{n=k}^{N} F_n
C     where F_n = T_n * S_n * A_n, matching the original semantics.
      DO KK = NPASSED-1, 1, -1
        PASSEDRAD(:,KK) = PASSEDRAD(:,KK) + PASSEDRAD(:,KK+1)
      ENDDO
      DO KK = 1, NPASSED
        PASSEDRAD(:,KK) = PASSEDRAD(:,KK) / PASSEDTRANSMIT(KK)
      ENDDO

C     Add in the gradient components due to the radiance that is now
C     fully calculated using the saved properties from each
C     subgrid integration interval.
      CALL COMPUTE_RADIANCE_DERIVATIVE(NUMDER,
     .  PASSEDDELS, PASSEDPOINTS, PASSEDINTERP0,
     .  PASSEDINTERP1, TOTAL_EXT, PASSEDRAD, RAYGRAD,
     .  PASSEDTRANSMIT, PASSEDABSCELL, NPTS, OPTINTERPWT,
     .  INTERPPTR, MAXPG, NSTOKES, NPASSED, DELTAM,
     .  MAXSUBGRIDINTS, DEXTM, IERR, ERRMSG)
      IF (IERR .NE. 0) RETURN
C      CALL CPU_TIME(TIME2)
C      TIME_RADIANCE = TIME_RADIANCE + TIME2 - TIME1
      DEALLOCATE (PASSEDPOINTS, PASSEDRAD,PASSEDINTERP0,
     .            PASSEDINTERP1, PASSEDDELS, PASSEDABSCELL,
     .            PASSEDTRANSMIT)
      RETURN
      END

      SUBROUTINE COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR,
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, PHASEINTERPWT,
     .             DIRFLUX,DERIV_MAXNMICRO,
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
     .             SINGLESCATTER, DTEMP, DEXTM, DALBM, DFJ, TIME_SOURCE)
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
      REAL TIME_SOURCE(3)
      INTEGER ICELL, NSTOKES, NSTLEG, NPTS, ML, MM, NLM, NLEG, NUMPHASE
      INTEGER MAXNMICRO, DERIV_MAXNMICRO
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
      LOGICAL OUTOFDOMAIN, SINGLESCATTER
      REAL    SINGSCAT8(NSTOKES,8), OSINGSCAT8(NSTOKES,8)
      REAL    DLEG(NSTLEG,0:NLEG,DNUMPHASE), DEXT(MAXPG,NUMDER)
      REAL    DALB(MAXPG,NUMDER)
      REAL    PHASEWTP(MAXNMICRO,MAXPG,NPART)
      REAL    DPHASEWTP(DERIV_MAXNMICRO,MAXPG,NUMDER)
      REAL    DTEMP(MAXPG,NUMDER)
      REAL    EXTINCTP(MAXPG,NPART), ALBEDOP(MAXPG, NPART)
      REAL    DEXTM(MAXPG,NUMDER), DALBM(8,NPTS,NUMDER)
      REAL    DFJ(8,NPTS,NUMDER)
      REAL    TEMP(NPTS)
      REAL    OPTINTERPWT(8,NPTS)
      INTEGER INTERPPTR(8,NPTS)
      INTEGER IPHASEP(MAXNMICRO,MAXPG,NPART)
      INTEGER DIPHASEP(DERIV_MAXNMICRO,MAXPG,NUMDER), KK
      REAL    DSINGSCAT(NSTOKES,DNUMPHASE)
      INTEGER LOFJ(*), PARTDER(NUMDER), RNS, RIS, IDR, IB,NB
      INTEGER IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I,IPA,K
      REAL    SECMU0, F, DA, A, A1, B1, EXT
      INTEGER IDP
      REAL GRAD8(NSTOKES,8,8,NUMDER), SOURCET(NSTOKES)
      REAL OGRAD8(NSTOKES,8,8,NUMDER)
      INTEGER DOEXACT(NUMDER)
      REAL WAVELEN, WAVENO(2)

      REAL XI, TRUNCSINGSCAT(NSTOKES)
      REAL LEGENT(NSTLEG,0:NLEG)
      REAL SINGSCATJ(NSTOKES)
      REAL SCATTERJ
      DOUBLE PRECISION EXTMIN, SCATMIN
      REAL SPATIAL_WEIGHT
      INTEGER Q
      REAL TIME1, TIME2, TIME3, TIME4

C     Locals for inlined COMPUTE_GRADIENT_ONEPROPPOINT
      REAL DLEGP(NSTLEG, 0:NLEG), UNSCALED_LEGEN(NSTLEG, 0:NLEG)
      REAL DLEGT(NSTLEG, 0:NLEG), LEGENP(NSTLEG, 0:NLEG)
      REAL DSINGSCATP(NSTOKES), DSOURCE(NSTOKES)
      REAL SINGSCATP(NSTOKES), DPLANCK, PLANCK, FTEMP
      REAL LEG_DIFF(NSTLEG, 0:NLEG)
      INTEGER LAST_IPA

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
C         Calculate data for new points. First the forward radiance then
C         the gradient.
C          CALL CPU_TIME(TIME1)
          EXT = TOTAL_EXT(IP)
          OLDIPTS(N) = IP
          IS = SHPTR(IP)
          NS = SHPTR(IP+1)-IS
          GRAD8(:,:,N,:) = 0.0

          RIS = RSHPTR(IP)
          RNS = RSHPTR(IP+1)-RIS

C             Sum over the real generalized spherical harmonic series
C             of the source function
          SRCEXT8(:,N) = 0.0
          SINGSCAT8(:,N) = 0.0
          DO J = 1, NS
            SRCEXT8(1,N) = SRCEXT8(1,N) + SOURCE(1,IS+J)*YLMDIR(1,J)
          ENDDO
C
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
C         Special case for solar source and delta-M
          IF ((SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B')) THEN

           DO IPA = 1, NPART
             IF (EXT.EQ.0.0) THEN
                 W = 1.0
             ELSE
                 W = EXTINCT(IP,IPA)/EXT
             ENDIF
             IF (W.EQ.0.0) CYCLE

             IF (INTERPMETHOD(2:2) .EQ. 'O' ) THEN
               LEGENT = LEGEN(:,:,IPHASE(1,IP,IPA))
C               F = LEGENT(1,ML+1)
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
C               F = LEGENT(1,ML+1)
             ENDIF
             IF (DELTAM) THEN
               F = LEGENT(1,ML+1)
               LEGENT(:,0:ML) = LEGENT(:,0:ML)/(1-F)
             ENDIF
C               First subtract off the truncated single scattering
             DA = ALBEDO(IP,IPA)*DIRFLUX(IP)*SECMU0*W
             J = 1
             TRUNCSINGSCAT = 0.0
             DO L = 0, ML
               ME = MIN(L,MM)
               MS = -ME
               A1 = DA*LEGENT(1,L)
               B1 = DA*LEGENT(5,L)
               IF (J .LE. NS) THEN
                 JT = J
                 DO M = MS, ME
                   TRUNCSINGSCAT(1) = TRUNCSINGSCAT(1) +
C                   SRCEXT8(1,N) = SRCEXT8(1,N) -
     .              A1*YLMDIR(1,J)*YLMSUN(1,J)
                   J = J + 1
                 ENDDO
                 IF (NSTOKES .GT. 1) THEN
                   J = JT
                   DO M = MS, ME
                   TRUNCSINGSCAT(2) = TRUNCSINGSCAT(2) +
C                     SRCEXT8(2,N) = SRCEXT8(2,N) -
     .                 B1*YLMDIR(2,J)*YLMSUN(1,J)
                   TRUNCSINGSCAT(3) = TRUNCSINGSCAT(3) +
C                     SRCEXT8(3,N) = SRCEXT8(3,N) -
     .                 B1*YLMDIR(6,J)*YLMSUN(1,J)
                     J = J + 1
                  ENDDO
                ENDIF
              ENDIF
             ENDDO
C
C               Then add in the single scattering contribution for the
C               original unscaled phase function.
             IF (DELTAM) THEN
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
               SRCEXT8(:,N) = SRCEXT8(:,N) - TRUNCSINGSCAT(:)
             ELSE
               SINGSCAT8(:,N) = SINGSCAT8(:,N) + TRUNCSINGSCAT(:)
             ENDIF
           ENDDO
           IF (DELTAM) THEN
              SRCEXT8(:,N) = SRCEXT8(:,N) + SINGSCAT8(:,N)
           ENDIF
          ENDIF

C          CALL CPU_TIME(TIME2)
C          TIME_SOURCE(1) = TIME_SOURCE(1) + TIME2 - TIME1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         BEGIN gradient calculations.
C         Restructured: group by species, hoist invariants,
C         inline COMPUTE_GRADIENT_ONEPROPPOINT.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C         Precompute PLANCK/DPLANCK once per RTE grid point (thermal).
          IF (SRCTYPE .EQ. 'T' .OR. SRCTYPE .EQ. 'B') THEN
            CALL PLANCK_DERIVATIVE(TEMP(IP), UNITS, WAVENO,
     .                             WAVELEN, DPLANCK)
            CALL PLANCK_FUNCTION(TEMP(IP), UNITS, WAVENO,
     .                           WAVELEN, PLANCK)
            IF (PLANCK .LT. 1e-7) DPLANCK = 0.0
          ENDIF

          LAST_IPA = -1
          DO IDR = 1, NUMDER
            IPA = PARTDER(IDR)

C           === Per-species setup (skip if same IPA as previous IDR) ===
            IF (IPA .NE. LAST_IPA) THEN
              LAST_IPA = IPA

C             Compute SCATTERJ and SINGSCATJ on RTE grid for this species.
              SCATTERJ = 0.0
              SINGSCATJ = 0.0
              IF (DELTAM) THEN
                DO NB=1,8
                  IB = INTERPPTR(NB,IP)
                  XI = OPTINTERPWT(NB,IP)
                  SPATIAL_WEIGHT = XI*ALBEDOP(IB,IPA)*EXTINCTP(IB,IPA)
                  SCATTERJ = SCATTERJ + SPATIAL_WEIGHT
                  IF (SPATIAL_WEIGHT .LE. 1e-6) CYCLE
                  DO Q=1,MAXNMICRO
                    IF (PHASEWTP(Q,IB,IPA) .LE. 1e-6) CYCLE
                    SINGSCATJ = SINGSCATJ + SPATIAL_WEIGHT*
     .                PHASEWTP(Q,IB,IPA)*
     .                SINGSCAT(:,IPHASEP(Q,IB,IPA))
                  ENDDO
                ENDDO
                IF (SCATTERJ .GT. SCATMIN) THEN
                  SINGSCATJ = SINGSCATJ/SCATTERJ
                ELSE
                  SINGSCATJ = SINGSCATJ/SCATMIN
                ENDIF
              ENDIF

C             Compute SOURCET and LEGENT for this species on RTE grid.
              IF (NPART .EQ. 1) THEN
                IF (ALBEDO(IP,IPA) .GT. 1e-8) THEN
                  SOURCET = SRCEXT8(:,N)/ALBEDO(IP,IPA)
                ELSE
                  SOURCET = 0.0
                ENDIF
              ELSE
                SOURCET = 0.0
                LEGENT = 0.0
                F = 0.0
                DO NB=1,8
                  IB = INTERPPTR(NB,IP)
                  XI = OPTINTERPWT(NB,IP)
                  SPATIAL_WEIGHT = XI*ALBEDOP(IB,IPA)*EXTINCTP(IB,IPA)
                  IF (SPATIAL_WEIGHT .LE. 1e-6) CYCLE
                  DO Q=1,MAXNMICRO
                    IF (PHASEWTP(Q,IB,IPA) .LE. 1e-6) CYCLE
                    LEGENT = LEGENT +
     .                SPATIAL_WEIGHT*PHASEWTP(Q,IB,IPA)*
     .                LEGEN(:,:,IPHASEP(Q,IB,IPA))
                  ENDDO
                ENDDO
                IF (SCATTERJ .GT. SCATMIN) THEN
                  LEGENT = LEGENT/SCATTERJ
                ELSE
                  LEGENT = LEGENT/SCATMIN
                ENDIF
                IF (DELTAM) THEN
                  F = LEGENT(1,ML+1)
                  IF (INTERPMETHOD(2:2) .EQ. 'N') THEN
                    LEGENT(:,0:ML) = LEGENT(:,0:ML)/(1-F)
                  ENDIF
                ENDIF
                IF (SCATTERJ .GT. SCATMIN) THEN
                  CALL COMPUTE_SOURCE_DIRECTION(LEGENT, SOURCET,
     .              NLM, LOFJ, RIS, RNS, NSTOKES, RADIANCE, YLMDIR,
     .              NSTLEG, DELTAM, SRCTYPE, ML, MM, NS, YLMSUN,
     .              DIRFLUX(IP), SECMU0, NLEG)
                  IF (DELTAM .AND. (SRCTYPE .EQ. 'S' .OR.
     .                SRCTYPE .EQ. 'B')) THEN
                    SOURCET = SOURCET +
     .                DIRFLUX(IP)*SINGSCATJ*SECMU0/(1-F)
                  ENDIF
                ENDIF
              ENDIF
              SOURCET(1) = MAX(0.0, SOURCET(1))

C             Undo delta-M scaling to get unscaled Legendre coefficients.
              IF (DELTAM) THEN
                LEGENT(:,0:ML) = LEGENT(:,0:ML)*(1-F)
                LEGENT(1,0:ML) = LEGENT(1,0:ML) + F
                IF (NSTLEG .GT. 1) THEN
                  LEGENT(2:4,0:ML) = LEGENT(2:4,0:ML) + F
                ENDIF
              ENDIF
            ENDIF
C           === End per-species setup ===

C           Loop over property grid points. For each NB, compute
C           LEGENP and SINGSCATP (invariant across IDRs for same IPA),
C           then compute per-IDR gradient (DLEGP, DLEGT, DSOURCE).
            DO NB = 1, 8
              IB = INTERPPTR(NB,IP)
              XI = OPTINTERPWT(NB,IP)
              IF (XI .LT. 1e-7) CYCLE

C             --- Build LEGENP at this property grid point ---
              LEGENP = 0.0
              IF (DELTAM) THEN
                DO Q = 1, MAXNMICRO
                  UNSCALED_LEGEN = LEGEN(:,:,IPHASEP(Q,IB,IPA))
                  FTEMP = UNSCALED_LEGEN(1,ML+1)
                  IF (INTERPMETHOD(2:2) .EQ. 'O') THEN
                    UNSCALED_LEGEN(1,0:ML) =
     .                UNSCALED_LEGEN(1,0:ML)*(1-FTEMP)
                  ENDIF
                  UNSCALED_LEGEN(1,0:ML) =
     .              UNSCALED_LEGEN(1,0:ML) + FTEMP
                  IF (NSTLEG .GT. 1) THEN
                    UNSCALED_LEGEN(2:4,0:ML) =
     .                UNSCALED_LEGEN(2:4,0:ML) + FTEMP
                  ENDIF
                  LEGENP = LEGENP +
     .              PHASEWTP(Q,IB,IPA)*UNSCALED_LEGEN
                ENDDO
              ELSE
                DO Q = 1, MAXNMICRO
                  LEGENP = LEGENP +
     .              PHASEWTP(Q,IB,IPA)*LEGEN(:,:,IPHASEP(Q,IB,IPA))
                ENDDO
              ENDIF

C             Precompute LEG_DIFF = LEGENP - LEGENT (reused across IDRs).
              LEG_DIFF = LEGENP - LEGENT

C             --- Build SINGSCATP at this property point (solar+deltam) ---
              IF ((SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B')
     .            .AND. DELTAM) THEN
                SINGSCATP = 0.0
                DO Q = 1, MAXNMICRO
                  SINGSCATP(:) = SINGSCATP(:)
     .              + PHASEWTP(Q,IB,IPA)*
     .                SINGSCAT(:,IPHASEP(Q,IB,IPA))
                ENDDO
              ENDIF

C             --- Per-IDR gradient computation (inlined) ---
C             Build DLEGP and DSINGSCATP (IDR-specific).
              DLEGP = 0.0
              DSINGSCATP = 0.0
              DO Q = 1, DERIV_MAXNMICRO
                IF (DOEXACT(IDR) .EQ. 1) THEN
                  IF (DELTAM) THEN
                    DSINGSCATP = DSINGSCATP +
     .                PHASEWTP(Q,IB,IPA)*
     .                DSINGSCAT(:,DIPHASEP(Q,IB,IDR))
                  ENDIF
                  DLEGP = DLEGP +
     .              PHASEWTP(Q,IB,IPA)*
     .              DLEG(:,:,DIPHASEP(Q,IB,IDR))
                ELSE IF (DOEXACT(IDR) .EQ. 0) THEN
                  UNSCALED_LEGEN = LEGEN(:,:,IPHASEP(Q,IB,IPA))
                  IF (DELTAM) THEN
                    DSINGSCATP = DSINGSCATP +
     .                DPHASEWTP(Q,IB,IDR)*
     .                SINGSCAT(:,IPHASEP(Q,IB,IPA))
                    FTEMP = UNSCALED_LEGEN(1,ML+1)
                    IF (INTERPMETHOD(2:2) .EQ. 'O') THEN
                      UNSCALED_LEGEN(1,0:ML) =
     .                  UNSCALED_LEGEN(1,0:ML)*(1-FTEMP)
                    ENDIF
                    UNSCALED_LEGEN(1,0:ML) =
     .                UNSCALED_LEGEN(1,0:ML) + FTEMP
                    IF (NSTLEG .GT. 1) THEN
                      UNSCALED_LEGEN(2:4,0:ML) =
     .                  UNSCALED_LEGEN(2:4,0:ML) + FTEMP
                    ENDIF
                  ENDIF
                  DLEGP = DLEGP +
     .              DPHASEWTP(Q,IB,IDR)*UNSCALED_LEGEN
                ENDIF
              ENDDO

C             Compute DLEGT = phase derivative * extinction/albedo.
              DLEGT = DEXT(IB,IDR)*LEG_DIFF*ALBEDOP(IB,IPA) +
     .          DALB(IB,IDR)*LEG_DIFF*EXTINCTP(IB,IPA) +
     .          DLEGP*EXTINCTP(IB,IPA)*ALBEDOP(IB,IPA) +
     .          (LEGENT-1)*DFJ(NB,IP,IDR)

C             Convolve DLEGT with radiance harmonics.
              DSOURCE = 0.0
              CALL COMPUTE_SOURCE_DIRECTION(DLEGT, DSOURCE,
     .          NLM, LOFJ, RIS, RNS, NSTOKES, RADIANCE, YLMDIR,
     .          NSTLEG, DELTAM, SRCTYPE, ML, MM, NS, YLMSUN,
     .          DIRFLUX(IP), SECMU0, NLEG)

C             Add extinction/albedo/phase gradients.
              GRAD8(:,NB,N,IDR) = GRAD8(:,NB,N,IDR) +
     .          XI*(SOURCET(:)*(ALBEDO(IP,IPA)*DEXTM(IB,IDR)
     .            + DALBM(NB,IP,IDR)) + DSOURCE(:))

C             Solar single scatter component (delta-M only).
              IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
                IF (DELTAM) THEN
                  GRAD8(:,NB,N,IDR) = GRAD8(:,NB,N,IDR) +
     .              DIRFLUX(IP)*SECMU0*XI*(
     .                SINGSCATJ*DFJ(NB,IP,IDR)
     .              + DSINGSCATP*EXTINCTP(IB,IPA)*ALBEDOP(IB,IPA)
     .              + DALB(IB,IDR)*(SINGSCATP-SINGSCATJ)*
     .                EXTINCTP(IB,IPA)
     .              + DEXT(IB,IDR)*(SINGSCATP-SINGSCATJ)*
     .                ALBEDOP(IB,IPA))
                ENDIF
              ENDIF

C             Thermal component.
              IF (SRCTYPE .EQ. 'T' .OR. SRCTYPE .EQ. 'B') THEN
                GRAD8(1,NB,N,IDR) = GRAD8(1,NB,N,IDR) +
     .            XI*(EXTINCT(IP,IPA)*(1.0-ALBEDO(IP,IPA))*
     .              DPLANCK*DTEMP(IB,IDR)
     .            - PLANCK*DALBM(NB,IP,IDR)
     .            + PLANCK*(1.0-ALBEDO(IP,IPA))*DEXTM(IB,IDR))
              ENDIF

            ENDDO
          ENDDO
C          CALL CPU_TIME(TIME2)
C          TIME_SOURCE(2) = TIME_SOURCE(2) + TIME2 - TIME1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         End of gradient calculations.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Final part of the forward model where we multiply the source
C      by the extinction.
          SINGSCAT8(:,N) = SINGSCAT8(:,N)*EXT
C         The SRCEXT8 is used to calculate the radiance, which is itself
C         used in the gradient. If we only want the single scatter (SINGLESCATTER) then
C         the source should only include the single scatter so we set that
C         here. This doesn't save computation time at all though.
          IF (SINGLESCATTER) THEN
            SRCEXT8(:,N) = SINGSCAT8(:,N)
          ELSE
            SRCEXT8(:,N) = SRCEXT8(:,N)*EXT
          ENDIF
          EXTINCT8(N) = EXT
        ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE COMPUTE_TOP_RADIANCES_GRAD (SRCTYPE, SKYRAD,WAVENO,
     .                                  WAVELEN,
     .                                  UNITS, NTOPPTS, NSTOKES, BCRAD,
     .                                  NPHI0MAX, NMU, IMU, IPHI,
     .                                  MU, PHI, MUS, PHIS, NPHI0,
     .                                  INTERPOLATE_FLAG,RADGRAD)
C       Returns the beginning of the BCRAD array with the downwelling
C     radiance for the NTOPPTS top boundary points.  Currently all
C     points have the same unpolarized isotropic radiance.
      IMPLICIT NONE
      INTEGER NTOPPTS, NSTOKES
Cf2py intent(in) :: NTOPPTS, NSTOKES
      INTEGER NMU, IMU, NPHI0MAX, IPHI
Cf2py intent(in) :: NMI, IMU, NPHI0MAX, IPHI
      REAL MU, PHI, MUS(NMU), PHIS(NMU, NPHI0MAX)
Cf2py intent(in) :: MU, PHI, MUS, PHIS
      INTEGER NPHI0(NMU)
Cf2py intent(in) :: NPHI0
      INTEGER INTERPOLATE_FLAG
Cf2py intent(in) :: INTERPOLATE_FLAG
      REAL    SKYRAD(NSTOKES,NMU/2,NPHI0MAX), WAVENO(2), WAVELEN
Cf2py intent(in) :: SKYRAD, WAVENO, WAVELEN
      REAL    BCRAD(NSTOKES, *)
Cf2py intent(in, out) :: BCRAD
      CHARACTER  SRCTYPE*1, UNITS*1
Cf2py intent(in) :: SRCTYPE, UNITS
      REAL    RADGRAD(NSTOKES,*)
Cf2py intent(out) :: RADGRAD

      INTEGER IBC, K
      REAL    SKYRAD2(NSTOKES), SKYRAD3(NSTOKES)

      INTEGER I,J, IANG
      DOUBLE PRECISION POWER, DISTANCE, WEIGHTEDSUM(NSTOKES)
      DOUBLE PRECISION WEIGHTSUM, WEIGHT

C           At top, boundary radiance either directly evaluated
C       at discrete ordinate points or interpolated to a specific
C       MU, PHI. The latter is only done when evaluating a radiance
C       (not during the solution iterations) and only affects
C       upward looking (ground based) instruments.

      IF (INTERPOLATE_FLAG .EQ. 1) THEN
C     Inverse distance weighting interpolation (Cubic).
C     Distance is based on the scattering angle between the two angles.
        POWER = 3.0D0
        WEIGHTEDSUM = 0.0D0
        WEIGHTSUM = 0.0D0
        DO I = 1, NMU/2
          DO J=1, NPHI0(I)
            DISTANCE = ACOS(MU*MUS(I) +
     .          SQRT((1.0-MU**2)*(1.0-MUS(I)**2))*
     .          COS(PHI-PHIS(I,J)))
            IF (ABS(DISTANCE) .LT. 1E-6) THEN
              WEIGHT = 1.0D8
            ELSE
              WEIGHT = 1.0D0/(DISTANCE**POWER)
            ENDIF
            WEIGHTEDSUM(:) = WEIGHTEDSUM(:) + SKYRAD(:,I,J)*WEIGHT
            WEIGHTSUM = WEIGHTSUM + WEIGHT
          ENDDO
        ENDDO
        SKYRAD3 = WEIGHTEDSUM/WEIGHTSUM
      ELSE IF (INTERPOLATE_FLAG .EQ. 2) THEN
C       For surface.
        IANG = 1
        POWER = 3.0D0
        WEIGHTEDSUM = 0.0D0
        WEIGHTSUM = 0.0D0
        DO I = NMU/2 + 1, NMU
          DO J=1, NPHI0(I)
            DISTANCE = ACOS(MU*MUS(I) +
     .          SQRT((1.0-MU**2)*(1.0-MUS(I)**2))*
     .          COS(PHI-PHIS(I,J)))
            IF (ABS(DISTANCE) .LT. 1E-6) THEN
              WEIGHT = 1.0D8
            ELSE
              WEIGHT = 1.0D0/(DISTANCE**POWER)
            ENDIF
            WEIGHTEDSUM(:) = WEIGHTEDSUM(:) + 
     .              SKYRAD(:,I - NMU/2,J)*WEIGHT
            WEIGHTSUM = WEIGHTSUM + WEIGHT
            RADGRAD(1,IANG) = WEIGHT
            IANG = IANG + 1
          ENDDO
        ENDDO
        SKYRAD3 = WEIGHTEDSUM/WEIGHTSUM
        RADGRAD(:,:IANG) = RADGRAD(:,:IANG)/WEIGHTSUM
      ELSE
        SKYRAD3 = SKYRAD(:,IMU,IPHI)
      ENDIF
      IF (SRCTYPE .EQ. 'T') THEN
        SKYRAD2(1:) = 0.0
          CALL PLANCK_FUNCTION (SKYRAD3(1), UNITS,
     .                        WAVENO, WAVELEN, SKYRAD2(1))
      ELSE
        SKYRAD2 = SKYRAD3
      ENDIF

C         Loop over all points assigning the uniform radiance
      DO IBC = 1, NTOPPTS
        BCRAD(:,IBC) = SKYRAD2(:)
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
     .                      RADBND, SFCGRIDRAD, NANG,
     .                      UNITS, WAVENO, BOUNDPTS, BOUNDINTERP,DIRRAD,
     .                      GNDALBEDO,SFCGRAD_RAY,NSFCDER,SFCDER,
     .                      TRANSMIT, SFCSRCGRAD_RAY,
     .                      DELXSFC,DELYSFC, DTEMP, TEMP)
C       Returns the interpolated Stokes radiance at the boundary (RADBND).
C     Inputs are the boundary location (XB,YB), ray direction away from
C     boundary (MU2,PHI2), cell number (ICELL) and face (KFACE) at
C     the boundary point.
      IMPLICIT NONE
      INTEGER NSTOKES, ICELL, KFACE, MAXNBC, NTOPPTS, NBOTPTS
      INTEGER GRIDPTR(8,*), BCPTR(MAXNBC,2)
      INTEGER NMU, NPHI0MAX, NPHI0(*), NSFCPAR, NANG
      REAL    MU2, PHI2
      INTEGER BOUNDPTS(4), SFCDER(NSFCDER), NSFCDER
      INTEGER NXSFC, NYSFC
      DOUBLE PRECISION BOUNDINTERP(4), DIRRAD(NSTOKES,4)
      DOUBLE PRECISION XB, YB, TRANSMIT
      DOUBLE PRECISION    SFCGRAD_RAY(NSTOKES,NSFCDER,*)
      REAL    SFCGRAD(NSTOKES,NSFCDER,4)
      REAL    GRIDPOS(3,*), RADBND(NSTOKES)
      REAL    WTDO(NMU,*), MU(NMU), PHI(NMU,*)
      REAL    WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
      REAL    SFCGRIDPARMS(NSFCPAR,*), BCRAD(NSTOKES,*)
      REAL    SFCGRIDRAD(NANG/2 + 1, *), WAVENO(2)
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
      REAL GNDALBEDO, DTEMP(*), TEMP(*)
      DOUBLE PRECISION SFCSRCGRAD_RAY(NSTOKES,NANG/2 + 1, *)
      REAL    DELXSFC, DELYSFC

      REAL    SFCRAD_TEMP(NSTOKES,NMU/2,NPHI0MAX)
      DOUBLE PRECISION U, V
      INTEGER IL, IM, IU, IP, IBC, J,I,IANG,JA
      LOGICAL LAMBERTIAN
      REAL    X(4), Y(4), RAD(NSTOKES,4), OPI
      REAL    RADEMIS(NSTOKES,4)
      REAL RADEMISGRAD(NSTOKES,NANG/2 + 1,4)
      REAL    REFLECT(4,4), DPLANCK
      REAL RX, RY, SFCU, SFCV, SFCINTERP(4)
      INTEGER SFCPTR(4)
      INTEGER IX, IY, ISFC, Q

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
C           Need to get DTEMP for that property grid point
C           Need surface 
            CALL PLANCK_DERIVATIVE(SFCGRIDPARMS(1,IBC), UNITS, WAVENO,
     .                          WAVELEN, DPLANCK)
            CALL VARIABLE_BRDF_SURFACE_GRAD (NBOTPTS,IBC,IBC, 
     .             BCPTR(1,2),
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,
     .             SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES,
     .             BCRAD(:,1+NTOPPTS), DIRRAD(1,J), NSFCDER,
     .             SFCDER, SFCGRAD(:,:,J), DPLANCK)

            CALL SURFACE_BRDF (SFCTYPE(2:2), SFCGRIDPARMS(2,IBC),
     .                    WAVELEN, MU2, PHI2, SOLARMU,SOLARAZ,
     .                      NSTOKES, REFLECT)
            DIRRAD(:,J) = OPI*REFLECT(1:NSTOKES,1)*DIRFLUX(IP)
          ENDIF
C         Hack to use COMPUTE_TOP_RADIANCES to compute the surface
C         emission.
          IANG = 1
C          SFCRAD_TEMP = 0.0
          DO I = 1, NMU/2
            DO JA = 1, NPHI0(I)
              SFCRAD_TEMP(1,I,JA) = SFCGRIDRAD(IANG+1,IBC)
              IANG = IANG + 1
            ENDDO
          ENDDO
          CALL COMPUTE_TOP_RADIANCES_GRAD(SRCTYPE,SFCRAD_TEMP,WAVENO,
     .                  WAVELEN,
     .                  UNITS,1, NSTOKES, RADEMIS(1,J), NPHI0MAX,
     .                  NMU,-1,-1,MU2,PHI2,MU,PHI,NPHI0,2,
     .                  RADEMISGRAD)
          RAD(:,J) = RADEMIS(:,J) + BCRAD(:,NTOPPTS+IBC)

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

      DO J=1,4
        IP = GRIDPTR(GRIDFACE(J,KFACE),ICELL)
        BOUNDPTS(J) = IP
        X(J) = GRIDPOS(1,IP)
        Y(J) = GRIDPOS(2,IP)

        RX = X(J)/DELXSFC
        RY = Y(J)/DELYSFC
        IX = MAX(1,MIN(NXSFC,INT(RX)+1))
        IY = MAX(1,MIN(NYSFC,INT(RY)+1))

        SFCPTR(1) = IX   + (NXSFC+1)*IY
        SFCPTR(2) = IX   + (NXSFC+1)*(IY+1)
        SFCPTR(3) = IX+1 + (NXSFC+1)*IY
        SFCPTR(4) = IX+1 + (NXSFC+1)*(IY+1)

        SFCU = MAX(0.0,MIN(1.0,RX-(IX-1)))
        SFCV = MAX(0.0,MIN(1.0,RY-(IY-1)))

        SFCINTERP(1) = (1-SFCU)*(1-SFCV)
        SFCINTERP(2) = SFCU*(1-SFCV)
        SFCINTERP(3) = (1-SFCU)*SFCV
        SFCINTERP(4) = SFCU*SFCV
        DO Q=1,4
          ISFC =SFCPTR(Q)
          SFCGRAD_RAY(:,:,ISFC) = SFCGRAD_RAY(:,:,ISFC) +
     .          TRANSMIT*BOUNDINTERP(J)*SFCGRAD(:,:,J)*
     .          SFCINTERP(Q)
          SFCSRCGRAD_RAY(:,:,ISFC) = 
     .          SFCSRCGRAD_RAY(:,:,ISFC) + 
     .          TRANSMIT*BOUNDINTERP(J)*RADEMISGRAD(:,:,J)*
     .          SFCINTERP(Q)
        ENDDO
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


      SUBROUTINE COMPUTE_RADIANCE_DERIVATIVE(NUMDER,
     .  PASSEDDELS, PASSEDPOINTS, PASSEDINTERP0, PASSEDINTERP1,
     .  TOTAL_EXT, PASSEDRAD, RAYGRAD, PASSEDTRANSMIT, PASSEDABSCELL,
     .  NPTS, OPTINTERPWT, INTERPPTR, MAXPG,NSTOKES, NPASSED,
     .  DELTAM, MAXSUBGRIDINTS, DEXTM, IERR, ERRMSG)
C     .  DEXT, DALB, LEGEN, DLEG, DIPHASEP, MAXNMICRO,
C     .  IPHASEP, PHASEWTP, DPHASEWTP, INTERPMETHOD, ALBEDOP,
C     .  EXTINCTP, DOEXACT,
C     .  NPTS, OPTINTERPWT, INTERPPTR, MAXPG,NPART
C     .  ML,, NSTLEG, NLEG, DNUMPHASE,
C.  NUMPHASE,
C     .  NSTOKES, NPASSED, DELTAM, MAXSUBGRIDINTS,
C     .  DEXTM)
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
C      INTEGER NPTS, NSTLEG, NLEG, MAXNMICRO, MAXPG, NUMDER, NPART
      INTEGER NSTOKES, NPTS, MAXPG, NUMDER
C      INTEGER PARTDER(NUMDER), NSTOKES, NUMPHASE, DNUMPHASE,ML
      INTEGER NPASSED, MAXSUBGRIDINTS
      LOGICAL DELTAM
      INTEGER IERR
      CHARACTER INTERPMETHOD*2, ERRMSG*600
C      REAL DEXT(MAXPG,NUMDER), DALB(MAXPG,NUMDER)
C      REAL ALBEDOP(MAXPG,NPART), EXTINCTP(MAXPG,NPART)
      REAL OPTINTERPWT(8,NPTS)
      INTEGER INTERPPTR(8,NPTS)
      REAL DEXTM(MAXPG,NUMDER)
C      REAL LEGEN(NSTLEG,0:NLEG,NUMPHASE)
C      REAL PHASEWTP(MAXNMICRO,MAXPG,NPART)
C      REAL DLEG(NSTLEG,0:NLEG,DNUMPHASE)
C      REAL DPHASEWTP(MAXNMICRO,MAXPG,NUMDER)
C      INTEGER DIPHASEP(MAXNMICRO,MAXPG,NUMDER)
C      INTEGER IPHASEP(MAXNMICRO,MAXPG,NPART)
C      INTEGER DOEXACT(NUMDER)
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
            DO IDR=1,NUMDER
              DO NB=1,8
                IB = INTERPPTR(NB,IP)
                XI = OPTINTERPWT(NB,IP)
                EXTGRAD = DEXTM(IB,IDR)*XI
C                CALL EXTINCTION_DERIVATIVE_POINT(
C     .            DEXT(IB,IDR), DOEXACT(IDR), XI, LEGEN,
C     .            NLEG, NSTLEG, NUMPHASE, DNUMPHASE,
C     .            PHASEWTP(:,IB,IPA),
C     .            DELTAM, IPHASEP(:,IB,IPA),
C     .            DIPHASEP(:,IB,IDR),EXTGRAD, DLEG,
C     .            MAXNMICRO, DALB(IB,IDR), ALBEDOP(IB,IPA),
C     .            EXTINCTP(IB,IPA), DPHASEWTP(:,IB,IDR), ML)

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
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END


      SUBROUTINE EXTINCTION_DERIVATIVE_POINT(DEXT, DOEXACT, XI,
     .      LEGEN, NLEG, NSTLEG, NUMPHASE, DNUMPHASE, PHASEWTP,
     .      DELTAM, IPHASEP, DIPHASEP, EXTGRAD, DLEG, MAXNMICRO,
     .      DALB, ALBEDOP, EXTINCTP, DPHASEWTP, ML,
     .      DERIV_MAXNMICRO)
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
      REAL PHASEWTP(MAXNMICRO), DPHASEWTP(DERIV_MAXNMICRO)
Cf2py intent(in) :: PHASEWTP, DPHASEWTP
      INTEGER IPHASEP(MAXNMICRO), DIPHASEP(DERIV_MAXNMICRO)
Cf2py intent(in) :: IPHASEP, DIPHASEP
      INTEGER DOEXACT, ML, DERIV_MAXNMICRO
Cf2py intent(in) :: DOEXACT, ML, DERIV_MAXNMICRO
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
      IF (DELTAM) THEN
        DO Q=1,MAXNMICRO
          IF (PHASEWTP(Q) .LE. 1E-6) CYCLE
          FP = FP + PHASEWTP(Q)*
     .            LEGEN(1,ML+1,IPHASEP(Q))
        ENDDO     

        IF (DOEXACT .EQ. 0) THEN
          DO Q=1,DERIV_MAXNMICRO
            DFP = DFP + DPHASEWTP(Q)*
     .             LEGEN(1,ML+1,IPHASEP(Q))
          ENDDO
        ELSE IF (DOEXACT .EQ. 1) THEN
          DO Q=1,DERIV_MAXNMICRO
            IF (PHASEWTP(Q) .LE. 1E-6) CYCLE
            DFP = DFP + PHASEWTP(Q)*
     .        DLEG(1,ML+1,DIPHASEP(Q))
          ENDDO
        ENDIF
      ENDIF

    !   DO Q=1,MAXNMICRO
    !     IF (PHASEWTP(Q) .LE. 1E-6) CYCLE
    !     IF (DELTAM) THEN
    !       FP = FP + PHASEWTP(Q)*
    !  .                LEGEN(1,ML+1,IPHASEP(Q))
    !       IF (DOEXACT .EQ. 1) THEN
    !         DFP = DFP + PHASEWTP(Q)*
    !  .                  DLEG(1,ML+1,DIPHASEP(Q))
    !       ELSEIF (DOEXACT .EQ. 0) THEN
    !         DFP = DFP + DPHASEWTP(Q)*
    !  .                  LEGEN(1,ML+1,IPHASEP(Q))
    !       ENDIF
    !     ENDIF
    !   ENDDO
      EXTGRAD = XI*(DEXT*(1-FP*ALBEDOP) -
     .            DALB*FP*EXTINCTP -
     .            EXTINCTP*ALBEDOP*DFP)
      RETURN
      END


      SUBROUTINE COMPUTE_DIRECT_BEAM_DERIV(DPATH, DPTR,
     .     NUMDER,DEXTM, TRANSMIT, ABSCELL,
     .     INPUTWEIGHT, NSTOKES, MAXPG,
     .     RAYGRAD, LONGEST_PATH_PTS)
C     .  ALBEDOP, PHASEWTP, LEGEN, NSTLEG, NUMPHASE,
C     .  MAXNMICRO, IPHASEP, RAYGRAD, ML, NLEG, DELTAM,
C     .  DNUMPHASE, DIPHASEP, DALB, EXTINCTP, DPHASEWTP,
C     .  DLEG, DOEXACT)
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
      INTEGER NSTOKES, NUMDER
Cf2py intent(in) :: NSTOKES, NUMDER
      INTEGER MAXPG
Cf2py intent(in) :: MAXPG
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL DPATH(LONGEST_PATH_PTS)
Cf2py intent(in) :: DPATH
      INTEGER DPTR(LONGEST_PATH_PTS)
Cf2py intent(in) :: DPTR
      REAL DEXTM(MAXPG,NUMDER)
Cf2py intent(in) :: DEXTM
      DOUBLE PRECISION TRANSMIT, ABSCELL
Cf2py intent(in) :: TRANSMIT, ABSCELL
      REAL INPUTWEIGHT(NSTOKES)
Cf2py intent(in) :: INPUTWEIGHT
      DOUBLE PRECISION RAYGRAD(NSTOKES,MAXPG,NUMDER)
Cf2py intent(in,out) :: RAYGRAD
      INTEGER LONGEST_PATH_PTS
Cf2py intent(in) :: LONGEST_PATH_PTS

      INTEGER IB, II, IDR, Q, IPA
      REAL FP, EXTGRAD

      II=1
      DO WHILE (DPTR(II) .GT. 0)
        IB = DPTR(II)
        DO IDR=1,NUMDER

          RAYGRAD(:,IB,IDR) = RAYGRAD(:,IB,IDR) -
     .      DEXTM(IB,IDR)*DPATH(II)*ABSCELL*
     .      TRANSMIT*INPUTWEIGHT(:)
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

      INTEGER J


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
     .    (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B')) THEN
        DO J = 1, NLM
          SOURCET(1) = SOURCET(1)
     .      + DIRFLUX*SECMU0*LEGEN(1,LOFJ(J))*YLMSUN(1,J)*YLMDIR(1,J)
        ENDDO
        IF (NSTOKES .GT. 1) THEN
          DO J = 5, NLM
            SOURCET(2) = SOURCET(2)
     .      + DIRFLUX*SECMU0*LEGEN(5,LOFJ(J))*YLMSUN(1,J)*YLMDIR(2,J)
          ENDDO
        ENDIF

       ENDIF

       RETURN
       END


      SUBROUTINE PREPARE_DERIV_INTERPS(GRIDPOS, NPTS,
     .    NPX, NPY, NPZ, MAXPG, DELX, DELY, XSTART, YSTART,
     .    ZLEVELS, OPTINTERPWT, INTERPPTR, IERR, ERRMSG,
     .    LEGEN, NUMPHASE, NSTLEG, DLEG, DNUMPHASE,
     .    PHASEWTP, DPHASEWTP, IPHASEP, DIPHASEP,
     .    DALBM, DEXTM, DFJ, NLEG, MAXNMICRO,
     .    ALBEDOP, EXTINCTP, NPART, DALB, DEXT,
     .    NUMDER, PARTDER, DOEXACT, ML, DELTAM,
     .    ALBEDO, IPHASE, PHASEINTERPWT, PHASEMAX,
     .    INTERPMETHOD, DERIV_MAXNMICRO)
C     Precomputes the interpolation weights (OPTINTERPWT) and pointers
C     (INTERPPTR) from the property grid onto the RTE grid, and computes
C     delta-M scaled derivative arrays DEXTM, DALBM, DFJ.
C
C     Restructured into three phases:
C       Phase 1: Compute INTERPPTR/OPTINTERPWT once per RTE point.
C       Phase 2: Compute FP, DFP, DEXTM on property grid directly.
C       Phase 3: Compute DALBM, DFJ using precomputed FP/DFP + RTE-grid F.
      IMPLICIT NONE
      INTEGER NPTS, MAXPG, NUMDER, PARTDER(NUMDER)
Cf2py intent(in) :: NPTS, MAXPG, NUMDER, PARTDER
      INTEGER NPX, NPY, NPZ, ML, NLEG, NPART
Cf2py intent(in) :: NPX, NPY, NPZ, ML, NLEG, NPART
      REAL GRIDPOS(3, NPTS)
Cf2py intent(in) :: GRIDPOS
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL OPTINTERPWT(8, NPTS)
Cf2py intent(out) :: OPTINTERPWT
      INTEGER INTERPPTR(8,NPTS)
Cf2py intent(out) :: INTERPPTR
      REAL DALBM(8,NPTS,NUMDER)
Cf2py intent(out) :: DALBM
      REAL DEXTM(MAXPG,NUMDER)
Cf2py intent(out) :: DEXTM
      REAL DFJ(8, NPTS,NUMDER)
Cf2py intent(out) :: DFJ
      REAL DELX, DELY, XSTART, YSTART
Cf2py intent(in) :: DELX, DELY, XSTART, YSTART
      REAL ZLEVELS(*)
Cf2py intent(in) :: ZLEVELS
      INTEGER NUMPHASE, NSTLEG, DNUMPHASE, MAXNMICRO
Cf2py intent(in) :: NUMPHASE, NSTLEG, DNUMPHASE, MAXNMICRO
      INTEGER DOEXACT(NUMDER), DERIV_MAXNMICRO
Cf2py intent(in) :: DOEXACT, DERIV_MAXNMOCRO
      REAL LEGEN(NSTLEG,0:NLEG,NUMPHASE)
      REAL DLEG(NSTLEG,0:NLEG, DNUMPHASE)
Cf2py intent(in) :: LEGEN, DLEG
      REAL PHASEWTP(MAXNMICRO, MAXPG,NPART)
      REAL DPHASEWTP(DERIV_MAXNMICRO,MAXPG,NUMDER)
      INTEGER IPHASEP(MAXNMICRO,MAXPG,NPART)
      INTEGER DIPHASEP(DERIV_MAXNMICRO,MAXPG,NUMDER)
Cf2py intent(in) :: DPHASEWTP, DPHASEWTP, IPHASEP, DIPHASEP
      REAL ALBEDOP(MAXPG, NPART), EXTINCTP(MAXPG,NPART)
      REAL DALB(MAXPG,NUMDER), DEXT(MAXPG,NUMDER)
Cf2py intent(in) :: ALBEDOP, EXTINCTP, DALB, DEXT
      REAL ALBEDO(NPTS,NPART)
Cf2py intent(in) :: ALBEDO
      REAL PHASEINTERPWT(8*MAXNMICRO, NPTS,NPART)
      INTEGER IPHASE(8*MAXNMICRO,NPTS,NPART)
      REAL PHASEMAX
Cf2py intent(in) :: PHASEINTERPWT, IPHASE, PHASEMAX
      CHARACTER INTERPMETHOD*2
Cf2py intent(in) :: INTERPMETHOD
      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG

      INTEGER IP, IDR, Q, IPA, NB, IB
      REAL F, ALBEDOJ, DIVIDE, FP_VAL, DFP_VAL
      REAL, ALLOCATABLE :: FP_BUF(:,:), DFP_BUF(:,:)

      IERR = 0

C     === Phase 1: Compute INTERPPTR/OPTINTERPWT once per RTE point ===
      DO IP = 1, NPTS
        CALL COMPUTE_INTERP_WEIGHTS(GRIDPOS(1,IP), GRIDPOS(2,IP),
     .    GRIDPOS(3,IP), NPX, NPY, NPZ, DELX, DELY, XSTART,
     .    YSTART, ZLEVELS, INTERPPTR(:,IP), OPTINTERPWT(:,IP),
     .    IERR, ERRMSG)
        IF (IERR .NE. 0) RETURN
      ENDDO

C     === Phase 2: Compute FP, DFP, DEXTM on property grid directly ===
      ALLOCATE(FP_BUF(MAXPG, NUMDER), DFP_BUF(MAXPG, NUMDER))
      DO IDR = 1, NUMDER
        IPA = PARTDER(IDR)
        DO IB = 1, MAXPG
          FP_VAL = 0.0
          DFP_VAL = 0.0
          IF (DELTAM) THEN
            DO Q = 1, DERIV_MAXNMICRO
              FP_VAL = FP_VAL + PHASEWTP(Q,IB,IPA)*
     .          LEGEN(1,ML+1,IPHASEP(Q,IB,IPA))
              IF (DOEXACT(IDR) .EQ. 1) THEN
                DFP_VAL = DFP_VAL + PHASEWTP(Q,IB,IPA)*
     .            DLEG(1,ML+1,DIPHASEP(Q,IB,IDR))
              ELSEIF (DOEXACT(IDR) .EQ. 0) THEN
                DFP_VAL = DFP_VAL + DPHASEWTP(Q,IB,IDR)*
     .            LEGEN(1,ML+1,IPHASEP(Q,IB,IPA))
              ENDIF
            ENDDO
          ENDIF
          FP_BUF(IB,IDR) = FP_VAL
          DFP_BUF(IB,IDR) = DFP_VAL
          DEXTM(IB,IDR) = DEXT(IB,IDR)*(1 - FP_VAL*
     .      ALBEDOP(IB,IPA)) - DALB(IB,IDR)*FP_VAL*
     .      EXTINCTP(IB,IPA) - EXTINCTP(IB,IPA)*
     .      ALBEDOP(IB,IPA)*DFP_VAL
        ENDDO
      ENDDO

C     === Phase 3: Compute DALBM, DFJ using RTE-grid F + precomputed FP/DFP ===
      DO IDR = 1, NUMDER
        IPA = PARTDER(IDR)
        DO IP = 1, NPTS

C         Compute F and ALBEDOJ on the RTE grid
          F = 0.0
          IF (DELTAM) THEN
            IF (INTERPMETHOD(2:2) .EQ. 'O') THEN
              F = LEGEN(1,ML+1,IPHASE(1,IP,IPA))
            ELSEIF (INTERPMETHOD(2:2) .EQ. 'N') THEN
              IF (PHASEINTERPWT(1,IP,IPA) .GE. PHASEMAX) THEN
                F = LEGEN(1,ML+1,IPHASE(1,IP,IPA))
              ELSE
                DO Q = 1, 8*MAXNMICRO
                  IF (PHASEINTERPWT(Q,IP,IPA) .LT. 1e-7) CYCLE
                  F = F + PHASEINTERPWT(Q,IP,IPA)*
     .              LEGEN(1,ML+1,IPHASE(Q,IP,IPA))
                ENDDO
              ENDIF
            ENDIF
            ALBEDOJ = ALBEDO(IP,IPA)/(F*(ALBEDO(IP,IPA)-1) + 1)
          ELSE
            ALBEDOJ = ALBEDO(IP,IPA)
          ENDIF
          DIVIDE = 1.0/(1.0 - ALBEDOJ*F)

          DO NB = 1, 8
            IB = INTERPPTR(NB,IP)
            DALBM(NB,IP,IDR) = DIVIDE*(
     .        DEXT(IB,IDR)*((1-F)*(ALBEDOP(IB,IPA) - ALBEDOJ)
     .          + (ALBEDOJ-1)*ALBEDOP(IB,IPA)*
     .            (FP_BUF(IB,IDR)-F))
     .        + DALB(IB,IDR)*((1-F)*EXTINCTP(IB,IPA)
     .          + (ALBEDOJ-1)*EXTINCTP(IB,IPA)*
     .            (FP_BUF(IB,IDR)-F))
     .        + DFP_BUF(IB,IDR)*(ALBEDOJ-1)*
     .          EXTINCTP(IB,IPA)*ALBEDOP(IB,IPA))
            DFJ(NB,IP,IDR) = (DEXT(IB,IDR)*
     .        (FP_BUF(IB,IDR)-F)*ALBEDOP(IB,IPA)
     .        + DALB(IB,IDR)*(FP_BUF(IB,IDR)-F)*
     .          EXTINCTP(IB,IPA)
     .        + DFP_BUF(IB,IDR)*EXTINCTP(IB,IPA)*
     .          ALBEDOP(IB,IPA))/(1-F)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(FP_BUF, DFP_BUF)
      RETURN
      END


      SUBROUTINE COMPUTE_INTERP_WEIGHTS (X, Y, Z, NPX, NPY, NPZ,
     .  DELX, DELY, XSTART, YSTART, ZLEVELS,
     .  INTERPPTR, OPTINTERPWT, IERR, ERRMSG)
C     Computes trilinear interpolation weights and pointers from the
C     property grid onto the RTE grid for a single RTE grid point.
      IMPLICIT NONE
      REAL    X, Y, Z
      INTEGER NPX, NPY, NPZ
      REAL    DELX, DELY, XSTART, YSTART
      REAL    ZLEVELS(*)
      INTEGER INTERPPTR(8)
      REAL    OPTINTERPWT(8)
      INTEGER IERR
      CHARACTER ERRMSG*600

      INTEGER IX, IXP, IY, IYP, IZ, IL, IM, IU
      INTEGER I1, I2, I3, I4
      DOUBLE PRECISION U, V, W

C     Find the Z grid location via binary search
      IL = 0
      IU = NPZ
      DO WHILE (IU-IL .GT. 1)
        IM = (IU+IL)/2
        IF (Z .GE. ZLEVELS(IM)) THEN
          IL = IM
        ELSE
          IU = IM
        ENDIF
      ENDDO
      IZ = MAX(IL,1)
      W = DBLE(Z - ZLEVELS(IZ))/(ZLEVELS(IZ+1) - ZLEVELS(IZ))
      W = MAX( MIN( W, 1.0D0), 0.0D0)

C     Find X grid location
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

C     Find Y grid location
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

C     Compute trilinear interpolation weights
      OPTINTERPWT(1) = (1-U)*(1-V)*(1-W)
      OPTINTERPWT(2) =    U *(1-V)*(1-W)
      OPTINTERPWT(3) = (1-U)*   V *(1-W)
      OPTINTERPWT(4) =    U *   V *(1-W)
      OPTINTERPWT(5) = (1-U)*(1-V)*   W
      OPTINTERPWT(6) =    U *(1-V)*   W
      OPTINTERPWT(7) = (1-U)*   V *   W
      OPTINTERPWT(8) =    U *   V *   W

C     Compute interpolation pointers
      I1 = IZ + NPZ*(IY-1) + NPZ*NPY*(IX-1)
      I2 = IZ + NPZ*(IY-1) + NPZ*NPY*(IXP-1)
      I3 = IZ + NPZ*(IYP-1) + NPZ*NPY*(IX-1)
      I4 = IZ + NPZ*(IYP-1) + NPZ*NPY*(IXP-1)
      INTERPPTR(1) = I1
      INTERPPTR(2) = I2
      INTERPPTR(3) = I3
      INTERPPTR(4) = I4
      INTERPPTR(5) = I1+1
      INTERPPTR(6) = I2+1
      INTERPPTR(7) = I3+1
      INTERPPTR(8) = I4+1

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

      SUBROUTINE ADJOINT_INTEGRATE_1RAY (BCFLAG, IPFLAG, NSTOKES,
     .             NSTLEG,
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
     .             VALIDRAD, TOTAL_EXT, NPART,
     .             ADJ_WEIGHT, GRADOUT, NGRAD, BEAM_WEIGHT,
     .             RSHPTR, RADIANCE, LOFJ, PARTDER, NUMDER, DEXT,
     .             DALB, DLEG, MAXPG, DNUMPHASE, SOLARFLUX,
     .             NPX, NPY, NPZ, DELX, DELY, XSTART, YSTART, ZLEVELS,
     .             EXTDIRP, UNIFORMZLEV, DPHASETAB, DPATH, DPTR,
     .             EXACT_SINGLE_SCATTER, TAUTOL, SINGLESCATTER,
     .             GNDALBEDO, IERR, ERRMSG, TEMP, PHASEMAX,
     .             PHASEINTERPWT, OPTINTERPWT, INTERPPTR,
     .             MAXNMICRO, EXTINCTP, ALBEDOP, PHASEWTP, DPHASEWTP,
     .             IPHASEP, DIPHASEP, WAVENO, UNITS, EXTMIN, SCATMIN,
     .             DOEXACT, INTERPMETHOD, DTEMP, DEXTM,DALBM,
     .             DFJ, MAXSUBGRIDINTS, TRANSCUT, LONGEST_PATH_PTS,
     .             TIME_SOURCE, DERIV_MAXNMICRO,
     .             TIME_DIRECT_POINT, TIME_DIRECT_SURFACE,
     .             TIME_RADIANCE, TIME_SUBGRID, TIME_ALLOCATE,
     .          VERBOSE, NSFCDER,NSFCPTS,SFCDER,DELXSFC,
     .             DELYSFC, SFCGRAD_RAY, SFCSRCGRAD_RAY, NANG,
     .             SFCGRIDRAD)
C     Adjoint version of ADJOINT_INTEGRATE_1RAY.
C     Traverses the grid backward from (MU2,PHI2) computing the same
C     forward quantities (TRANSMIT, SRC, ABSCELL) but contracts the
C     NSTOKES dimension using ADJ_WEIGHT and accumulates directly into
C     GRADOUT(MAXPG,NUMDER,NGRAD) instead of per-ray RAYGRAD.
C     Also accumulates BEAM_WEIGHT(NPTS) for batched direct beam
C     derivatives (called once per grid point after all rays).

      IMPLICIT NONE
      LOGICAL EXACT_SINGLE_SCATTER, VERBOSE
      INTEGER NPX, NPY, NPZ, MAXPG, BCELL, NANG
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
      DOUBLE PRECISION ADJ_WEIGHT(NSTOKES)
      DOUBLE PRECISION GRADOUT(MAXPG,NUMDER,NGRAD)
      DOUBLE PRECISION BEAM_WEIGHT(NPTS)
      INTEGER NGRAD
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
      INTEGER   DERIV_MAXNMICRO
      REAL      EXTINCTP(MAXPG,NPART), ALBEDOP(MAXPG,NPART)
      REAL      PHASEWTP(MAXNMICRO,MAXPG,NPART)
      REAL      DPHASEWTP(DERIV_MAXNMICRO,MAXPG,NUMDER)
      INTEGER   IPHASEP(MAXNMICRO,MAXPG,NPART)
      INTEGER   DIPHASEP(DERIV_MAXNMICRO,MAXPG,NUMDER)
      REAL      WAVENO(2)
      DOUBLE PRECISION EXTMIN, SCATMIN
      INTEGER   DOEXACT(NUMDER)
      INTEGER   MAXSUBGRIDINTS

      REAL      SINGSCAT8(NSTOKES,8), OSINGSCAT8(NSTOKES,8)
      INTEGER   PARTDER(NUMDER), NUMDER, DIPHASE(MAXPG,NUMDER)
      INTEGER   DIPHASEIND(NPTS,NUMDER)
      REAL      SINGSCAT0(NSTOKES,8), SINGSCAT1(NSTOKES,8)
      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE, SSP
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2, J, L, II, IDR, IPA, KK
      INTEGER NS
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
      DOUBLE PRECISION CONTRIB
      INTEGER LONGEST_PATH_PTS
      REAL    DPATH(LONGEST_PATH_PTS,*), SECMU0
      INTEGER DPTR(LONGEST_PATH_PTS,*), N
      REAL DEXTM(MAXPG,NUMDER), DALBM(8,NPTS,NUMDER)
      REAL DFJ(8,NPTS,NUMDER)
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
      LOGICAL :: SINGLESCATTER
      DOUBLE PRECISION :: DIRRAD(NSTOKES,4), BOUNDINTERP(4)
      INTEGER :: BOUNDPTS(4), IB,IP
      REAL :: GNDALBEDO
      INTEGER IERR
      CHARACTER ERRMSG*600

      INTEGER NSFCDER, NSFCPTS, NXSFC, NYSFC
      INTEGER SFCDER(NSFCDER, NSFCPTS)
      DOUBLE PRECISION SFCGRAD_RAY(NSTOKES,NSFCDER, NSFCPTS)
      DOUBLE PRECISION SFCSRCGRAD_RAY(NSTOKES,NANG/2 + 1, NSFCPTS)
      REAL SFCGRIDRAD(NANG/2 + 1,*)
      REAL DELXSFC, DELYSFC      

      REAL TIME_SOURCE(3), TIME_DIRECT_POINT, TIME_DIRECT_SURFACE
      REAL TIME_RADIANCE, TIME_SUBGRID, TIME_ALLOCATE

      INTEGER, ALLOCATABLE :: PASSEDPOINTS(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDRAD(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDABSCELL(:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDTRANSMIT(:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDINTERP0(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDINTERP1(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: PASSEDDELS(:)
      INTEGER NPASSED, K, MAXBYOPT
      REAL TIME1, TIME2
      LOGICAL DOSINGSCAT

      DATA OPPFACE/2,1,4,3,6,5/
      DATA ONEY/0,0,-1,-2,0,0,-5,-6/, ONEX/0,-1,0,-3,0,-5,0,-7/
      DATA DONEFACE/0,0,0,0,0,0,0,0, 0,1,0,3,0,5,0,7, 2,0,4,0,6,0,8,0,
     .              0,0,1,2,0,0,5,6, 3,4,0,0,5,6,0,0,
     .              0,0,0,0,1,2,3,4, 5,6,7,8,0,0,0,0/

C         TRANSCUT is the transmission to stop the integration at
C      TRANSCUT = 5.0E-5
C0.0D0
C         TAUTOL is the maximum optical path for the subgrid intervals
C      it is now a numerical parameter.
C      TAUTOL = 0.2

      EPS = 1.0E-5*(GRIDPOS(3,GRIDPTR(8,1))-GRIDPOS(3,GRIDPTR(1,1)))
      MAXCELLSCROSS = 50*MAX(NX,NY,NZ)

      MAXSUBGRIDINTS = MAX(MAXCELLSCROSS, MAXSUBGRIDINTS)


      ALLOCATE (PASSEDPOINTS(8,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDINTERP0(8,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDINTERP1(8,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDDELS(MAXSUBGRIDINTS))
      ALLOCATE (PASSEDRAD(NSTOKES,MAXSUBGRIDINTS))
      ALLOCATE (PASSEDABSCELL(MAXSUBGRIDINTS))
      ALLOCATE (PASSEDTRANSMIT(MAXSUBGRIDINTS))

C     initialize counter for subgrid intervals.
      NPASSED = 1
      RADOUT = 0.0D0

      PI = ACOS(-1.0D0)
C       Calculate the generalized spherical harmonics for this direction
      CALL YLMALL (.FALSE.,SNGL(MU2),SNGL(PHI2),ML,MM,NSTLEG, YLMDIR)
      SECMU0 = 1.0D0/ABS(SOLARMU)
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

      IFACE = 0
      NGRID = 0
C         Loop until reach a Z boundary or transmission is very small
      VALIDRAD = .FALSE.
C      CALL CPU_TIME(TIME2)
C      TIME_ALLOCATE = TIME_ALLOCATE + TIME2 - TIME1
      DO WHILE (.NOT. VALIDRAD .AND. ICELL .GT. 0)

C           Make sure current cell is valid
        IF (ICELL .LE. 0) THEN
          IERR = 1
          WRITE (ERRMSG,*)'ADJOINT_INTEGRATE_1RAY: ICELL=',ICELL,
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
C      CALL CPU_TIME(TIME1)

        CALL COMPUTE_SOURCE_GRAD_1CELL (ICELL, GRIDPTR,
     .             NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, PHASEINTERPWT,
     .             DIRFLUX,DERIV_MAXNMICRO,
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
     .             SINGLESCATTER, DTEMP, DEXTM, DALBM, DFJ, TIME_SOURCE)

C         CALL CPU_TIME(TIME2)
C         TIME_SOURCE = TIME_SOURCE + TIME2 - TIME1
C         Interpolate the source and extinction to the current point
C       Note that this recalculation of SRCEXT1 when moving between cells
C       rather than inheritance is the cause of discrepancy in radiance
C       calculation between
C       INTEGRATE_SOURCE (original SHDOM) and INTEGRATE_1RAY/ADJOINT_INTEGRATE_1RAY
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
          WRITE (ERRMSG,*) 'ADJOINT_INTEGRATE_1RAY: SO<0  ',
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

        SRCSINGSCAT = 0.0

C           Loop over the subgrid cells
        DO IT = 1, NTAU

          DOSINGSCAT = .TRUE.
C         save the interpolation kernel and pointers for each subgrid
C         interval so that the radiance gradient can be calculated.
C          CALL CPU_TIME(TIME1)
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
C             multiplied by extinction derivative is calculated at the end of the
C             subroutine.
              SRCGRAD = ( 0.5*(GRAD0+GRAD1)
     .          + 0.08333333333*(EXT0*GRAD1-EXT1*GRAD0)*DELS
     .           *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT

               IF (EXACT_SINGLE_SCATTER .AND.
     .           SRCTYPE .NE. 'T') THEN
                 SRCSINGSCAT = SRCSINGSCAT + TRANSMIT*ABSCELL*
     .                    ( 0.5*(SINGSCAT0+SINGSCAT1)
     .            + 0.08333333333*(EXT0*SINGSCAT1-EXT1*SINGSCAT0)*DELS
     .                    *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
               ENDIF
            ENDIF

            RADOUT(:) = RADOUT(:) + TRANSMIT*SRC(:)*ABSCELL
C         The stuff below is only done when there is extinction. Strictly
C         speaking, it should only need to be performed when there is no
C         extinction from any of the particle species for which derivatives
C         are calculated for.

C         Store the radiance contribution for this subgrid interval.
C         The cumulative sum needed by COMPUTE_RADIANCE_DERIVATIVE
C         is built in a single O(N) backward pass after the ray loop.
            PASSEDABSCELL(NPASSED) = ABSCELL
            PASSEDTRANSMIT(NPASSED) = TRANSMIT
            PASSEDRAD(:,NPASSED) = TRANSMIT*SRC(:)*ABSCELL
C          CALL CPU_TIME(TIME2)
C          TIME_SUBGRID = TIME_SUBGRID + TIME2-TIME1

C          CALL CPU_TIME(TIME1)
            IF (.NOT. OUTOFDOMAIN) THEN
              DO KK=1,8
                IP = GRIDPTR(KK,ICELL)
                DO K=1,8
                  IB = INTERPPTR(K,IP)
                  DO IDR=1,NUMDER
                    CONTRIB = 0.0D0
                    DO NS=1,NSTOKES
                      CONTRIB = CONTRIB +
     .                  ADJ_WEIGHT(NS)*SRCGRAD(NS,K,KK,IDR)
                    ENDDO
                    GRADOUT(IB,IDR,1) = GRADOUT(IB,IDR,1) +
     .              TRANSMIT*CONTRIB*ABSCELL
                  ENDDO
                ENDDO
C                IF (EXACT_SINGLE_SCATTER .AND. (SRCTYPE .EQ. 'S'
C     .            .OR. SRCTYPE .EQ. 'B')) THEN
CC             Add gradient component due to the direct solar beam.
CC             Sensitivity of single scattered radiation to
CC             extinction along the path between the gridpoint and the sun.
C                  CALL COMPUTE_DIRECT_BEAM_DERIV(DPATH(:,IP),
C     .            DPTR(:,IP),
C     .            NUMDER, DEXTM, TRANSMIT,
C     .            ABSCELL, SRCSINGSCAT, NSTOKES,
C     .            MAXPG,RAYGRAD, LONGEST_PATH_PTS)
C
C                ENDIF
              ENDDO
            ENDIF
C          CALL CPU_TIME(TIME2)
C          TIME_DIRECT_POINT= TIME_DIRECT_POINT +
C     .      TIME2 - TIME1

            NPASSED = NPASSED + 1
            IF (NPASSED .GT. MAXSUBGRIDINTS) THEN
              IERR = 1
              WRITE(ERRMSG, *) "ADJOINT_INTEGRATE_1RAY: The maximum ",
     .        "number of ",
     .        " subgrid intervals for calculation of the",
     .        " radiance along the",
     .        " ray path has been exceeded.", "NPASSED=", NPASSED,
     .        "MAXSUBGRIDINTS=", MAXSUBGRIDINTS
              RETURN
            ENDIF

          ELSE
            ABSCELL = 0.0
            TRANSCELL = 1.0
            SRC(:) = 0.0
            SRCGRAD = 0.0
            SRCSINGSCAT = 0.0

          ENDIF

          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          SRCEXT1(:) = SRCEXT0(:)
          GRAD1 = GRAD0
          SINGSCAT1 = SINGSCAT0
C                End of sub grid cell loop
        ENDDO

C       C             Add gradient component due to the direct solar beam.
C             Sensitivity of single scattered radiation to
C             extinction along the path between the gridpoint and the sun.
        IF (EXACT_SINGLE_SCATTER .AND. (SRCTYPE .EQ. 'S'
     .       .OR. SRCTYPE .EQ. 'B')) THEN
          DO KK =1,8
            IP = GRIDPTR(KK,ICELL)
C           Accumulate adjoint-contracted single scatter weight
C           for batched COMPUTE_DIRECT_BEAM_DERIV_ADJOINT call.
            DO NS=1,NSTOKES
              BEAM_WEIGHT(IP) = BEAM_WEIGHT(IP) +
     .          ADJ_WEIGHT(NS)*SRCSINGSCAT(NS,KK)
            ENDDO
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
C
C           If the transmission is greater than zero and not at a
C             boundary then prepare for next cell
C        .OR. NGRID.GT.MAXCELLSCROSS Exceeding MAXCELLSCROSS is no
C        longer a valid reason for stopping ray integration as this
C        may compromise physical accuracy.
        IF ((TRANSMIT .LT. TRANSCUT)) THEN
          VALIDRAD = .TRUE.
C         No boundary radiance; initialise the sentinel entry.
          PASSEDRAD(:,NPASSED) = 0.0D0
          PASSEDTRANSMIT(NPASSED) = 1.0D0

        ELSE IF (INEXTCELL .EQ. 0 .AND. IFACE .GE. 5) THEN
          VALIDRAD = .TRUE.
C          CALL CPU_TIME(TIME1)
          CALL FIND_BOUNDARY_RADIANCE_GRAD (NSTOKES, XN, YN,
     .                      SNGL(MU2), SNGL(PHI2),
     .                      IC, KFACE, GRIDPTR, GRIDPOS,
     .                      MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .                      NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .                      SRCTYPE, WAVELEN, SOLARMU,SOLARAZ, DIRFLUX,
     .                      SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                      RADBND, SFCGRIDRAD, NANG,
     .                      UNITS, WAVENO, BOUNDPTS, BOUNDINTERP, 
     .                      DIRRAD,
     .                      GNDALBEDO,SFCGRAD_RAY,NSFCDER,SFCDER,
     .                      TRANSMIT, SFCSRCGRAD_RAY,
     .                      DELXSFC,DELYSFC,DTEMP,TEMP)
          RADOUT(:) = RADOUT(:) + TRANSMIT*RADBND(:)
C         These values aren't actually used. set to -1.0 for
C         debugging purposes.
          PASSEDTRANSMIT(NPASSED) = TRANSMIT
          PASSEDABSCELL(NPASSED) = -1.0
C         Store boundary radiance contribution in the last entry.
C         The backward cumsum after the ray loop will propagate it.
          PASSEDRAD(:,NPASSED) = TRANSMIT*RADBND(:)
C         Calculate the sensitivity of the surface reflection to
C         extinction along the path between the surface and the sun.
C         This does not include sensitivity of surface reflection to
C         downwelling diffuse radiation, only the reflection of the direct
C         solar beam.
          IF (EXACT_SINGLE_SCATTER .AND. (SRCTYPE .EQ. 'S'
     .            .OR. SRCTYPE .EQ. 'B')) THEN
C           Accumulate adjoint-contracted boundary direct beam weight.
            DO KK=1,4
              IP = BOUNDPTS(KK)
              DO NS=1,NSTOKES
                BEAM_WEIGHT(IP) = BEAM_WEIGHT(IP) +
     .            ADJ_WEIGHT(NS)*TRANSMIT*
     .            BOUNDINTERP(KK)*DIRRAD(NS,KK)
              ENDDO
            ENDDO
          ENDIF
C          CALL CPU_TIME(TIME2)
C          TIME_DIRECT_SURFACE = TIME_DIRECT_SURFACE +
C     .      TIME2 - TIME1
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
C     Build the cumulative future-radiance array via backward sum.
C     After this loop PASSEDRAD(:,k) = (1/T_k) * sum_{n=k}^{N} F_n
C     where F_n = T_n * S_n * A_n, matching the original semantics.
      DO KK = NPASSED-1, 1, -1
        PASSEDRAD(:,KK) = PASSEDRAD(:,KK) + PASSEDRAD(:,KK+1)
      ENDDO
      DO KK = 1, NPASSED
        PASSEDRAD(:,KK) = PASSEDRAD(:,KK) / PASSEDTRANSMIT(KK)
      ENDDO

C     Add in the gradient components due to the radiance that is now
C     fully calculated using the saved properties from each
C     subgrid integration interval.
      CALL COMPUTE_RADIANCE_DERIVATIVE_ADJOINT(NUMDER,
     .  PASSEDDELS, PASSEDPOINTS, PASSEDINTERP0,
     .  PASSEDINTERP1, TOTAL_EXT, PASSEDRAD, ADJ_WEIGHT,
     .  GRADOUT, NGRAD,
     .  PASSEDTRANSMIT, PASSEDABSCELL, NPTS, OPTINTERPWT,
     .  INTERPPTR, MAXPG, NSTOKES, NPASSED, DELTAM,
     .  MAXSUBGRIDINTS, DEXTM, IERR, ERRMSG)
      IF (IERR .NE. 0) RETURN
C      CALL CPU_TIME(TIME2)
C      TIME_RADIANCE = TIME_RADIANCE + TIME2 - TIME1
      DEALLOCATE (PASSEDPOINTS, PASSEDRAD,PASSEDINTERP0,
     .            PASSEDINTERP1, PASSEDDELS, PASSEDABSCELL,
     .            PASSEDTRANSMIT)
      RETURN
      END

      SUBROUTINE COMPUTE_ADJOINT_WEIGHTS(STOKESOUT, MEASUREMENT,
     .  UNCERTAINTIES, COSTFUNC, NSTOKES, NUNCERTAINTY,
     .  ADJ_WEIGHT, COST, NCOST, IERR, ERRMSG)
C     Computes adjoint weights for the gradient computation from the
C     cost function. These weights contract the NSTOKES dimension
C     early, enabling the adjoint gradient approach.
C     For L2: ADJ_WEIGHT(I) = sum_J U(I,J) * (R_I - M_I)
C     For LL: chain-rule through log transforms
      IMPLICIT NONE
      CHARACTER*2 COSTFUNC
      INTEGER NSTOKES, NUNCERTAINTY, NCOST
      DOUBLE PRECISION STOKESOUT(NSTOKES), MEASUREMENT(NSTOKES)
      DOUBLE PRECISION UNCERTAINTIES(NUNCERTAINTY,NUNCERTAINTY)
      DOUBLE PRECISION ADJ_WEIGHT(NSTOKES)
      DOUBLE PRECISION COST(NCOST)
      INTEGER IERR
      CHARACTER ERRMSG*600

      DOUBLE PRECISION PIXEL_ERROR, RADERROR, DOLP1, DOLP2, DOLPERR
      INTEGER I, J

      IERR = 0
      ADJ_WEIGHT = 0.0D0

      IF (COSTFUNC .EQ. 'L2') THEN
        DO I=1,NSTOKES
          PIXEL_ERROR = STOKESOUT(I) - MEASUREMENT(I)
          DO J=1,NSTOKES
            COST(1) = COST(1) + 0.5D0*
     .        UNCERTAINTIES(I,J) * PIXEL_ERROR**2
            ADJ_WEIGHT(I) = ADJ_WEIGHT(I) +
     .        UNCERTAINTIES(I,J)*PIXEL_ERROR
          ENDDO
        ENDDO

      ELSEIF (COSTFUNC .EQ. 'LL') THEN
        RADERROR = LOG(STOKESOUT(1)) - LOG(MEASUREMENT(1))
        COST(1) = COST(1) + 0.5D0*(RADERROR**2 *
     .      UNCERTAINTIES(1,1))
        ADJ_WEIGHT(1) = ADJ_WEIGHT(1) + RADERROR*
     .      UNCERTAINTIES(1,1)/STOKESOUT(1)

        IF (NSTOKES .GT. 1) THEN
          DOLP1 = SQRT(STOKESOUT(2)**2 + STOKESOUT(3)**2) /
     .      STOKESOUT(1)
          DOLP2 = SQRT(MEASUREMENT(2)**2 + MEASUREMENT(3)**2) /
     .      MEASUREMENT(1)
          DOLPERR = LOG(DOLP1) - LOG(DOLP2)

          COST(1) = COST(1) + 0.5D0*(DOLPERR**2 *
     .      UNCERTAINTIES(2,2))
          ADJ_WEIGHT(2) = ADJ_WEIGHT(2) + DOLPERR*
     .      UNCERTAINTIES(2,2)*STOKESOUT(2)/
     .      (STOKESOUT(2)**2 + STOKESOUT(3)**2)
          ADJ_WEIGHT(3) = ADJ_WEIGHT(3) + DOLPERR*
     .      UNCERTAINTIES(2,2)*STOKESOUT(3)/
     .      (STOKESOUT(2)**2 + STOKESOUT(3)**2)
        ENDIF
      ELSE
        IERR = 1
        WRITE(ERRMSG,*) 'COMPUTE_ADJOINT_WEIGHTS: Bad cost function',
     .    COSTFUNC
        RETURN
      ENDIF
      RETURN
      END


      SUBROUTINE COMPUTE_RADIANCE_DERIVATIVE_ADJOINT(NUMDER,
     .  PASSEDDELS, PASSEDPOINTS, PASSEDINTERP0, PASSEDINTERP1,
     .  TOTAL_EXT, PASSEDRAD, ADJ_WEIGHT, GRADOUT, NGRAD,
     .  PASSEDTRANSMIT, PASSEDABSCELL,
     .  NPTS, OPTINTERPWT, INTERPPTR, MAXPG, NSTOKES, NPASSED,
     .  DELTAM, MAXSUBGRIDINTS, DEXTM, IERR, ERRMSG)
C     Adjoint version of COMPUTE_RADIANCE_DERIVATIVE.
C     Contracts the NSTOKES dimension using ADJ_WEIGHT and
C     accumulates directly into GRADOUT instead of RAYGRAD.
      IMPLICIT NONE
      INTEGER NSTOKES, NPTS, MAXPG, NUMDER, NGRAD
      INTEGER NPASSED, MAXSUBGRIDINTS
      LOGICAL DELTAM
      INTEGER IERR
      CHARACTER ERRMSG*600
      REAL OPTINTERPWT(8,NPTS)
      INTEGER INTERPPTR(8,NPTS)
      REAL DEXTM(MAXPG,NUMDER)
      REAL TOTAL_EXT(NPTS)
      DOUBLE PRECISION ADJ_WEIGHT(NSTOKES)
      DOUBLE PRECISION GRADOUT(MAXPG,NUMDER,NGRAD)
      DOUBLE PRECISION PASSEDRAD(NSTOKES,MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDABSCELL(MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDTRANSMIT(MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDINTERP0(8,MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDINTERP1(8,MAXSUBGRIDINTS)
      DOUBLE PRECISION PASSEDDELS(MAXSUBGRIDINTS)
      INTEGER PASSEDPOINTS(8,MAXSUBGRIDINTS)

      INTEGER IDR, IP, IB, K, NB, KK, NS
      DOUBLE PRECISION EXT0, EXT1, DELS, EXT
      REAL XI, EXTGRAD
      DOUBLE PRECISION GBASE0, GBASE1, ADJ_GBASE

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
            IP = PASSEDPOINTS(K,KK)
C           Precompute adjoint-contracted base for this (K,KK).
C           Factor out EXTGRAD since it depends only on (NB,IDR).
            GBASE0 = 0.0D0
            GBASE1 = 0.0D0
            DO NS=1,NSTOKES
              GBASE0 = GBASE0 - ADJ_WEIGHT(NS)*
     .          PASSEDRAD(NS,KK+1)*PASSEDINTERP0(K,KK)
              GBASE1 = GBASE1 - ADJ_WEIGHT(NS)*
     .          PASSEDRAD(NS,KK)*PASSEDINTERP1(K,KK)
            ENDDO
            ADJ_GBASE = ( 0.5*(GBASE0+GBASE1)
     .       + 0.08333333333*(EXT0*GBASE1-EXT1*GBASE0)*DELS
     .            *(1.0 - 0.05*(EXT1-EXT0)*DELS) )/EXT
            ADJ_GBASE = ADJ_GBASE *
     .        PASSEDTRANSMIT(KK)*PASSEDABSCELL(KK)

            DO IDR=1,NUMDER
              DO NB=1,8
                IB = INTERPPTR(NB,IP)
                XI = OPTINTERPWT(NB,IP)
                EXTGRAD = DEXTM(IB,IDR)*XI
                GRADOUT(IB,IDR,1) = GRADOUT(IB,IDR,1) +
     .            EXTGRAD*ADJ_GBASE
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END


      SUBROUTINE COMPUTE_DIRECT_BEAM_DERIV_ADJOINT(DPATH, DPTR,
     .     NUMDER, DEXTM, SCALAR_WEIGHT,
     .     MAXPG, GRADOUT, NGRAD, LONGEST_PATH_PTS)
C     Adjoint version of COMPUTE_DIRECT_BEAM_DERIV.
C     Uses a scalar weight (already contracted over NSTOKES)
C     and accumulates into GRADOUT instead of RAYGRAD.
      IMPLICIT NONE
      INTEGER NUMDER, MAXPG, NGRAD, LONGEST_PATH_PTS
      REAL DPATH(LONGEST_PATH_PTS)
      INTEGER DPTR(LONGEST_PATH_PTS)
      REAL DEXTM(MAXPG,NUMDER)
      DOUBLE PRECISION SCALAR_WEIGHT
      DOUBLE PRECISION GRADOUT(MAXPG,NUMDER,NGRAD)

      INTEGER IB, II, IDR

      II=1
      DO WHILE (DPTR(II) .GT. 0)
        IB = DPTR(II)
        DO IDR=1,NUMDER
          GRADOUT(IB,IDR,1) = GRADOUT(IB,IDR,1) -
     .      DEXTM(IB,IDR)*DPATH(II)*SCALAR_WEIGHT
        ENDDO
        II = II + 1
      ENDDO
      RETURN
      END
