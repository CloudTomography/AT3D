C     TODO: Write description for this file 

      SUBROUTINE RENDER (NX, NY, NZ, NPTS, NCELLS,
     .                   ML, MM, NCS, NLM, NLEG, NUMPHASE, 
     .                   NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .                   BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, 
     .                   SOLARAZ, SFCTYPE, NSFCPAR, SFCGRIDPARMS, 
     .                   MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .                   GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, 
     .                   UNITS, XGRID, YGRID, ZGRID, GRIDPOS, 
     .                   GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                   EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, 
     .                   FLUXES, SHPTR, SOURCE, CAMX, CAMY, CAMZ, CAMMU, 
     .                   CAMPHI, NPIX, VISOUT)
      IMPLICIT NONE
      INTEGER NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
Cf2py intent(in) :: NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
      INTEGER ML, MM, NCS, NLM, NLEG, NUMPHASE
Cf2py intent(in) :: ML, MM, NCS, NLM, NLEG, NUMPHASE
      INTEGER NMU, NPHI0MAX, NPHI0(*)
Cf2py intent(in) :: NMU, NPHI0MAX, NPHI0
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
Cf2py intent(in) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER SHPTR(*), BCPTR(MAXNBC,2)
Cf2py intent(in) :: SHPTR, BCPTR
      INTEGER*2 CELLFLAGS(*)
      INTEGER IPHASE(*)
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
      REAL    EXTINCT(*), ALBEDO(*), LEGEN(0:NLEG,*)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN
      REAL    DIRFLUX(*), FLUXES(2,*), SOURCE(*)
Cf2py intent(in) :: DIRFLUX, FLUXES, SOURCE
      REAL    CAMX(*), CAMY(*), CAMZ(*), CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX 
      REAL   VISOUT(NPIX)
Cf2py intent(out) :: VISOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS

      INTEGER NSCATANGLE, I, J, L
      INTEGER N  
      REAL    MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2
      DOUBLE PRECISION THETA0, THETA1, PHIR, PHI0, COSSCAT
      DOUBLE PRECISION U, R, PI
      INTEGER MAXNLM, MAXLEG, MAXPHASE, MAXSCATANG
      PARAMETER (MAXNLM=16384,MAXLEG=10000) 
      PARAMETER (MAXPHASE=500000,MAXSCATANG=361)
      REAL    YLMDIR(MAXNLM), YLMSUN(MAXNLM)
      REAL   PHASETAB(MAXPHASE,MAXSCATANG), SINGSCAT(MAXPHASE)
      DOUBLE PRECISION SUNDIRLEG(0:MAXLEG)
  
  
      IF ((ML+1)**2-(2-NCS)*(ML*(ML+1))/2 .GT. MAXNLM)
     .    STOP 'VISUALIZE_RADIANCE: MAXNLM exceeded'
          IF (NLEG .GT. MAXLEG) STOP 'RENDER: MAXLEG exceeded'
          IF (NUMPHASE .GT. MAXPHASE) THEN
              WRITE(*,*) 'NUMPHASE=', NUMPHASE, 'MAXPHASE=', MAXPHASE
              STOP 'RENDER: MAXPHASE exceeded.'
          ENDIF
  
          IF (SRCTYPE .NE. 'T') THEN
              CALL YLMALL (SOLARMU, SOLARAZ, ML, MM, NCS, YLMSUN)
              IF (DELTAM .AND. NUMPHASE .GT. 0) THEN
                  NSCATANGLE = MAX(36,MIN(MAXSCATANG,2*NLEG))
                  CALL PRECOMPUTE_PHASE(MAXPHASE, NSCATANGLE, NUMPHASE, 
     .                                  ML, NLEG, LEGEN, PHASETAB)
              ENDIF
          ENDIF
  
  
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
  
      IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
          COSSCAT = SOLARMU*MU2 + SQRT((1.0-SOLARMU**2)*(1.0-MU2**2))
     .                 *COS(SOLARAZ-PHI2)
          IF (NUMPHASE .GT. 0) THEN
              U = (NSCATANGLE-1)*(ACOS(COSSCAT)/PI) + 1
              J = MIN(NSCATANGLE-1,INT(U))
              U = U - J
              DO I = 1, NUMPHASE
                  SINGSCAT(I) = (1-U)*PHASETAB(I,J) + U*PHASETAB(I,J+1)
              ENDDO
          ELSE
              CALL LEGENDRE_ALL (COSSCAT, NLEG, SUNDIRLEG)
              DO L = 0, NLEG
                  SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L+1)/(4*PI)
              ENDDO
          ENDIF
          ENDIF
          CALL INTEGRATE_1RAY (BCFLAG, IPFLAG, 
     .                         NX, NY, NZ, NPTS, NCELLS, 
     .                         GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                         XGRID, YGRID, ZGRID, GRIDPOS,
     .                         ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .                         NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                         DELTAM, SRCTYPE, WAVELEN, SOLARMU,
     .                         SOLARAZ, EXTINCT, ALBEDO, LEGEN,  
     .                         IPHASE, DIRFLUX, SHPTR, SOURCE, 
     .                         YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .                         MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                         SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                         MURAY, PHIRAY, MU2, PHI2, X0, Y0, Z0, 
     .                         VISOUT(N))
  900     CONTINUE
      ENDDO
  
      RETURN
      END
      
      
      
      SUBROUTINE RAYLEIGH_EXTINCT (NZT, ZLEVELS,TEMP, RAYLCOEF, EXTRAYL)
C       Computes the molecular Rayleigh extinction profile EXTRAYL [/km]
C     from the temperature profile TEMP [K] at ZLEVELS [km].  Assumes
C     a linear lapse rate between levels to compute the pressure at
C     each level.  The Rayleigh extinction is proportional to air
C     density, with the coefficient RAYLCOEF in [K/(mb km)].
      IMPLICIT NONE
      INTEGER NZT
Cf2py intent(in) :: NZT
      REAL    ZLEVELS(NZT), TEMP(NZT), RAYLCOEF, EXTRAYL(NZT)
Cf2py intent(in) :: ZLEVELS, TEMP, RAYLCOEF
Cf2py intent(out) :: EXTRAYL
      INTEGER I
      REAL    PRES, LAPSE, TS, DZ

C           Find surface pressure by integrating hydrostatic relation
C           for a dry atmosphere up to surface height.
      PRES = 1013.
      TS = TEMP(1)
      LAPSE = 6.5*0.001
      PRES = PRES*(TS/(TS+LAPSE*ZLEVELS(1)*1000.))**(9.8/(287.*LAPSE))

C         Use layer mean temperature to compute fractional pressure change.
      DO I = 1, NZT-1
        EXTRAYL(I) = RAYLCOEF*PRES/TEMP(I)
        DZ = 1000.*(ZLEVELS(I+1)-ZLEVELS(I))
        LAPSE = (TEMP(I)-TEMP(I+1))/DZ
        IF (ABS(LAPSE) .GT. 0.00001) THEN
          PRES = PRES*(TEMP(I+1)/TEMP(I))**(9.8/(287.*LAPSE))
        ELSE
          PRES = PRES*EXP(-9.8*DZ/(287.*TEMP(I)))
        ENDIF
      ENDDO  
      EXTRAYL(NZT) = RAYLCOEF*PRES/TEMP(NZT)
      
      RETURN
      END