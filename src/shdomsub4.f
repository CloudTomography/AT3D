C     This file containts subroutines that were modified from their original purpose 
C     The original subroutines were written by Frank Evans for the Spherical Harmonic Discrete Ordinate Method for 3D Atmospheric Radiative Transfer.
C     The modified subroutines were written by Aviad Levis, Technion Institute of Technology, 2019

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
      PARAMETER (MAXPHASE=50000,MAXSCATANG=361)
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
      
      
      SUBROUTINE EXT_GRADIENT (NX, NY, NZ, NPTS, NBPTS, NCELLS,
     .             NBCELLS, ML, MM, NCS, NLM, NLEG, NUMPHASE, 
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, 
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,
     .             XGRID, YGRID, ZGRID, GRIDPOS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,
     .             SHPTR, SOURCE, CAMX, CAMY, CAMZ, CAMMU, CAMPHI, 
     .             NPIX, GRADOUT, COST, MEASUREMENTS, RSHPTR, 
     .             RADIANCE)

      IMPLICIT NONE
      INTEGER NX, NY, NZ, BCFLAG, IPFLAG, NPTS
      INTEGER NBPTS, NCELLS, NBCELLS
Cf2py intent(in) :: NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NBPTS, NCELLS, NBCELLS
      INTEGER ML, MM, NCS, NLM, NLEG
      INTEGER NUMPHASE
Cf2py intent(in) :: ML, MM, NCS, NLM, NLEG, NUMPHASE
      INTEGER NMU, NPHI0MAX, NPHI0(*)
Cf2py intent(in) :: NMU, NPHI0MAX, NPHI0
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
Cf2py intent(in) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER SHPTR(*), BCPTR(MAXNBC,2), RSHPTR(*)
Cf2py intent(in) :: SHPTR, BCPTR, RSHPTR
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
      REAL    DIRFLUX(*), FLUXES(2,*), SOURCE(*), RADIANCE(*)
Cf2py intent(in) :: DIRFLUX, FLUXES, SOURCE, RADIANCE
      REAL    CAMX(*), CAMY(*), CAMZ(*), CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX 
      REAL   MEASUREMENTS(*)
Cf2py intent(in) :: MEASUREMENTS
      REAL VISOUT
      DOUBLE PRECISION  GRAD(NPTS), GRADOUT(NBPTS), COST
Cf2py intent(out) :: GRADOUT, COST
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS

      INTEGER NSCATANGLE, I, J, L, K
      INTEGER N  
      DOUBLE PRECISION PIXEL_ERROR, RAYGRAD(NPTS)
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

      GRAD = 0.0
      IF ((ML+1)**2-(2-NCS)*(ML*(ML+1))/2 .GT. MAXNLM)
     .    STOP 'VISUALIZE_RADIANCE: MAXNLM exceeded'
      IF (NLEG .GT. MAXLEG)  STOP 'COMPUTE_ONE_SOURCE: MAXLEG exceeded'
      IF (NUMPHASE .GT. MAXPHASE)  
     .    STOP 'COMPUTE_ONE_SOURCE: MAXPHASE exceeded'

      IF (SRCTYPE .NE. 'T') THEN
        CALL YLMALL (SOLARMU, SOLARAZ, ML, MM, NCS, YLMSUN)
        IF (DELTAM .AND. NUMPHASE .GT. 0) THEN
          NSCATANGLE = MAX(36,MIN(MAXSCATANG,2*NLEG))
          CALL PRECOMPUTE_PHASE (MAXPHASE, NSCATANGLE, NUMPHASE, ML,
     .                           NLEG, LEGEN, PHASETAB)
        ENDIF
      ENDIF

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
      COST = 0
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
            VISOUT = 0.0
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
     .                  *COS(SOLARAZ-PHI2)
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
        
        RAYGRAD = 0.0
        CALL EXTGRAD_INTEGRATE_1RAY (BCFLAG, IPFLAG, 
     .                       NX, NY, NZ, NPTS, NCELLS, 
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       XGRID, YGRID, ZGRID, GRIDPOS,
     .                       ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .                       NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                       DELTAM, SRCTYPE, WAVELEN, SOLARMU,SOLARAZ,
     .                       EXTINCT, ALBEDO, LEGEN, IPHASE, 
     .                       DIRFLUX, SHPTR, SOURCE, 
     .                       YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .                       MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                       SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                       MURAY, PHIRAY, MU2, PHI2, X0, Y0, Z0, 
     .                       VISOUT, RAYGRAD, RSHPTR, RADIANCE)
900     CONTINUE
        
        PIXEL_ERROR = (VISOUT - MEASUREMENTS(N))
        
C        WRITE(*,*) N, VISOUT, MEASUREMENTS(N), PIXEL_ERROR, 
C     .              MINVAL(RAYGRAD), MAXVAL(RAYGRAM)

        GRAD = GRAD + PIXEL_ERROR*RAYGRAD
        COST = COST + 0.5*PIXEL_ERROR**2
      ENDDO
    
      CALL BASE_GRID_PROJECTION(NBCELLS, NCELLS, NBPTS, 
     .              GRIDPOS, GRIDPTR, TREEPTR, GRAD, GRADOUT)
   
      RETURN
      END
      
      SUBROUTINE EXTGRAD_INTEGRATE_1RAY (BCFLAG, IPFLAG, 
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
     .                        SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                        MURAY,PHIRAY, MU2, PHI2, X0, Y0, Z0, 
     .                        RADOUT, RAYGRAD, RSHPTR, RADIANCE)

C       Integrates the source function through the extinction field 
C     (EXTINCT) backward in the direction (MURAY,PHIRAY) to find the 
C     outgoing radiance (RAD) at the point X0,Y0,Z0.
      IMPLICIT NONE
      INTEGER BCFLAG, IPFLAG, NX, NY, NZ, NPTS, NCELLS
      INTEGER ML, MM, NCS, NLM, NLEG, NUMPHASE
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, KK, GRIDPOINT
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
      INTEGER SHPTR(NPTS+1), RSHPTR(NPTS+1)
      INTEGER*2 CELLFLAGS(NCELLS)
      INTEGER IPHASE(NPTS)
      INTEGER BCPTR(MAXNBC,2)
      LOGICAL DELTAM
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    WAVELEN, SOLARMU, SOLARAZ
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(*)
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(0:NLEG,NPTS) 
      REAL    DIRFLUX(NPTS)
      DOUBLE PRECISION RAYGRAD(*)
      REAL    SOURCE(*), RADIANCE(*)
      REAL    YLMDIR(NLM), YLMSUN(NLM), SINGSCAT(NUMPHASE)
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
      REAL    MURAY, PHIRAY, MU2, PHI2, RADOUT
      DOUBLE PRECISION X0, Y0, Z0
      CHARACTER SRCTYPE*1, SFCTYPE*2

      INTEGER BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
      INTEGER IOPP, NTAU, IT, I, IPT1, IPT2
      LOGICAL DONE, IPINX, IPINY, OPENBCFACE, OUTOFDOMAIN
      INTEGER JFACE, KFACE, IC, MAXCELLSCROSS, NGRID, K
      INTEGER OPPFACE(6), OLDIPTS(8), DONETHIS(8)
      INTEGER DONEFACE(8,7), ONEX(8), ONEY(8)
      REAL    OEXTINCT8(8), OSRCEXT8(8), EXTINCT8(8), SRCEXT8(8)
      REAL    EXT0, EXT1, EXTN, SRCEXT0, SRCEXT1, RADBND
      REAL    XM,YM, GRADFIELD8(8), OGRADFIELD8(8)
      DOUBLE PRECISION CX, CY, CZ, CXINV, CYINV, CZINV
      DOUBLE PRECISION XE, YE, ZE, XN, YN, ZN, XI, YI, ZI
      DOUBLE PRECISION SO, SOX, SOY, SOZ, EPS, F0(8), F1(8)
      DOUBLE PRECISION TAUTOL, TAUGRID, S, DELS, TRANSCUT, ABSCELL1
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
     .                  BCFLAG, IPFLAG, XE, YE, ZE,  ICELL)
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
          OGRADFIELD8(I) = GRADFIELD8(I)
        ENDDO
     
C         Compute the source function times extinction in direction (MU2,PHI2)
        CALL COMPUTE_SOURCE_1CELL (ICELL, GRIDPTR, 
     .             ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU, 
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, 
     .             SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, 
     .             SINGSCAT, DONETHIS, OLDIPTS, OEXTINCT8,
     .             OSRCEXT8, EXTINCT8, SRCEXT8)
     
     
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

        F1(1) = (1-W)*(1-V)*(1-U)
        F1(2) = (1-W)*(1-V)*U
        F1(3) = (1-W)*V*(1-U)
        F1(4) = (1-W)*V*U
        F1(5) = W*(1-V)*(1-U)
        F1(6) = W*(1-V)*U
        F1(7) = W*V*(1-U)
        F1(8) = W*V*U

        SRCEXT1 = F1(1)*SRCEXT8(1) + F1(2)*SRCEXT8(2) +
     .            F1(3)*SRCEXT8(3) + F1(4)*SRCEXT8(4) +
     .            F1(5)*SRCEXT8(5) + F1(6)*SRCEXT8(6) +
     .            F1(7)*SRCEXT8(7) + F1(8)*SRCEXT8(8)

                
        EXT1 = F1(1)*EXTINCT8(1) + F1(2)*EXTINCT8(2) +
     .         F1(3)*EXTINCT8(3) + F1(4)*EXTINCT8(4) +
     .         F1(5)*EXTINCT8(5) + F1(6)*EXTINCT8(6) +
     .         F1(7)*EXTINCT8(7) + F1(8)*EXTINCT8(8)


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
        
C	If this is not a boundary cell (currently assuming that the bc conditions are open) 
C	    Check if the zlevel is within the zparticle levels
C		Compute the gradient field in direction  (MU2,PHI2)
        OUTOFDOMAIN = (BTEST(INT(CELLFLAGS(ICELL)),0).OR.
     .                 BTEST(INT(CELLFLAGS(ICELL)),1))
        IF (.NOT.OUTOFDOMAIN) THEN
           CALL EXT_GRADFIELD_1CELL (ICELL, GRIDPTR, 
     .             ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU, 
     .             ALBEDO, LEGEN, IPHASE, DIRFLUX, 
     .             YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .             DONETHIS, OLDIPTS, SOURCE, RADIANCE,
     .             SHPTR, RSHPTR, OGRADFIELD8, GRADFIELD8)
        ENDIF



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
          U = (XI-GRIDPOS(1,IPT1))*INVDELX
          V = (YI-GRIDPOS(2,IPT1))*INVDELY
          W = (ZI-GRIDPOS(3,IPT1))*INVDELZ 
          F0(1) = (1-W)*(1-V)*(1-U)
          F0(2) = (1-W)*(1-V)*U
          F0(3) = (1-W)*V*(1-U)
          F0(4) = (1-W)*V*U
          F0(5) = W*(1-V)*(1-U)
          F0(6) = W*(1-V)*U
          F0(7) = W*V*(1-U)
          F0(8) = W*V*U

          SRCEXT0 = F0(1)*SRCEXT8(1) + F0(2)*SRCEXT8(2) +
     .              F0(3)*SRCEXT8(3) + F0(4)*SRCEXT8(4) +
     .              F0(5)*SRCEXT8(5) + F0(6)*SRCEXT8(6) +
     .              F0(7)*SRCEXT8(7) + F0(8)*SRCEXT8(8)
          IF (IT .NE. NTAU) THEN
            EXT0 = F0(1)*EXTINCT8(1) + F0(2)*EXTINCT8(2) +
     .             F0(3)*EXTINCT8(3) + F0(4)*EXTINCT8(4) +
     .             F0(5)*EXTINCT8(5) + F0(6)*EXTINCT8(6) +
     .             F0(7)*EXTINCT8(7) + F0(8)*EXTINCT8(8)
          ELSE
            EXT0 = EXTN
          ENDIF
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

C                Update the gradient 
        IF (.NOT.OUTOFDOMAIN) THEN
            DO KK=1,8
              GRIDPOINT = GRIDPTR(KK,ICELL)
             
             IF (EXT .NE. 0.0) THEN
                 ABSCELL1=ABSCELL/EXT
             ELSE
                 ABSCELL1 = 1.0
             ENDIF
             
             RAYGRAD(GRIDPOINT) = RAYGRAD(GRIDPOINT) +
     .	     0.5*(F0(KK)+F1(KK))*DELS*TRANSMIT*ABSCELL1
     .	     *GRADFIELD8(KK)
            
            ENDDO
          ENDIF

          RAD = RAD + TRANSMIT*SRC*ABSCELL
          TRANSMIT = TRANSMIT*TRANSCELL
          EXT1 = EXT0
          SRCEXT1 = SRCEXT0
          IF (IT .LT. NTAU) THEN
            DO KK=1,8
                F1(KK) = F0(KK)
            ENDDO
          ENDIF
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
      
      
        SUBROUTINE EXT_GRADFIELD_1CELL (ICELL, GRIDPTR,
     .             ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .             NPTS, DELTAM, SRCTYPE, SOLARMU, ALBEDO,
     .             LEGEN, IPHASE, DIRFLUX, YLMDIR, YLMSUN,
     .             SUNDIRLEG, SINGSCAT, DONETHIS, OLDIPTS,
     .             SOURCE, RADIANCE, SHPTR, RSHPTR, 
     .             OGRADFIELD8, GRADFIELD8) 
C       Computes the source function times extinction for gridpoints
C     belonging to cell ICELL in the direction (MU,PHI).  The results
C     are returned in SRCEXT8 and EXTINCT8.
C     The spherical harmonic source function series is input in SOURCE.
C     For a solar source if delta-M then use Nakajima and Tanaka TMS
C     procedure, replacing delta-M single scattering with single scattering
C     for unscaled untruncated phase function.
      IMPLICIT NONE
      INTEGER ICELL, NPTS, ML, MM, NCS, NLM, NLEG
      INTEGER NUMPHASE
      INTEGER GRIDPTR(8,*), RSHPTR(NPTS+1), SHPTR(NPTS+1)
      INTEGER DONETHIS(8), OLDIPTS(8)
      INTEGER IPHASE(NPTS)
      LOGICAL DELTAM
      REAL    SOLARMU
      REAL    ALBEDO(NPTS), LEGEN(0:NLEG,NPTS)
      REAL    DIRFLUX(NPTS), RADIANCE(*), SOURCE(*)
      REAL    YLMDIR(NLM), YLMSUN(NLM), SINGSCAT(NUMPHASE)
      REAL    GRADFIELD8(8), OGRADFIELD8(8)
      DOUBLE PRECISION SUNDIRLEG(0:NLEG)
      CHARACTER SRCTYPE*1

      INTEGER IP, J, L, M, MS, ME, K, SIS, RIS, SNS, RNS, N, I
      DOUBLE PRECISION DA, F, A, SECMU0

      SECMU0 = 1.0D0/ABS(SOLARMU)

C         Loop over the grid points, computing the source function
C           at the viewing angle from the spherical harmonic source function.
      DO N = 1, 8
        IP = GRIDPTR(N,ICELL)
        I = DONETHIS(N)
        IF (I .GT. 0 .AND. IP .EQ. OLDIPTS(N)) THEN
          GRADFIELD8(N) = OGRADFIELD8(I)
        ELSE IF (I .LT. 0) THEN
          GRADFIELD8(N) = GRADFIELD8(ABS(I))
        ELSE

          OLDIPTS(N) = IP
          SIS = SHPTR(IP)
          SNS = SHPTR(IP+1)-SIS
          RIS = RSHPTR(IP)
          RNS = RSHPTR(IP+1)-RIS
C             Sum over the spherical harmonic series of the source function
          GRADFIELD8(N) = 0.0
          DO J = 1, MIN(RNS,SNS)
            GRADFIELD8(N) = GRADFIELD8(N) + 
     .			(SOURCE(SIS+J)-RADIANCE(RIS+J))*YLMDIR(J)
          ENDDO
          IF (SNS .GT. RNS) THEN
            DO J = MIN(RNS,SNS)+1, MAX(RNS,SNS)
              GRADFIELD8(N) = GRADFIELD8(N) + 
     .			      SOURCE(SIS+J)*YLMDIR(J)
            ENDDO
          ELSEIF (SNS .LT. RNS) THEN
            DO J = MIN(RNS,SNS)+1, MAX(RNS,SNS)
              GRADFIELD8(N) = GRADFIELD8(N) - 
     .			      RADIANCE(SIS+J)*YLMDIR(J)
            ENDDO
          ENDIF

C             Special case for solar source and Delta-M
          IF (SRCTYPE .NE. 'T' .AND. DELTAM) THEN
            IF (NUMPHASE .GT. 0) THEN
              K = IPHASE(IP)
            ELSE
              K = IP
            ENDIF
C               First subtract off the truncated single scattering
            DA = ALBEDO(IP)*DIRFLUX(IP)*SECMU0
            J = 1
            DO L = 0, ML
              ME = MIN(L,MM)
              MS = (1-NCS)*ME
              A = DA*LEGEN(L,K)
              IF (J .LE. SNS) THEN
                DO M = MS, ME
                  GRADFIELD8(N) = GRADFIELD8(N) - A*YLMDIR(J)*YLMSUN(J)
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
              GRADFIELD8(N) = GRADFIELD8(N) + DA*SINGSCAT(K)
            ELSE
              F = LEGEN(ML+1,K)
              DO L = 0, NLEG
                IF (L .LE. ML) THEN
                  A = DA*(LEGEN(L,K) + F/(1-F))
                ELSE
                  A = DA*LEGEN(L,K)/(1-F)
                ENDIF
                GRADFIELD8(N) = GRADFIELD8(N) + A*SUNDIRLEG(L)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE PAR_RENDER (NX, NY, NZ, NPTS, NCELLS,
     .             ML, MM, NCS, NLM, NLEG, NUMPHASE, 
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS, 
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,
     .             XGRID, YGRID, ZGRID, GRIDPOS, 
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,
     .             SHPTR, SOURCE, CAMX, CAMY, CAMZ, CAMMU, CAMPHI, 
     .             NPIX, VISOUT, PHASETAB, YLMSUN, NSCATANGLE)

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
      DOUBLE PRECISION  CAMX(*), CAMY(*), CAMZ(*)
      REAL CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX
Cf2py intent(in) :: NPIX 
      REAL   VISOUT(NPIX)
Cf2py intent(out) :: VISOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS
      INTEGER NSCATANGLE, I, J, L 
Cf2py intent(in) :: NSCATANGLE
      REAL    MURAY, PHIRAY 
      INTEGER MAXNLM, MAXLEG, MAXPHASE, MAXSCATANG
      PARAMETER (MAXNLM=16384,MAXLEG=10000)
      PARAMETER (MAXPHASE=500000, MAXSCATANG=361)
      REAL    YLMDIR(MAXNLM), YLMSUN(MAXNLM)
cf2py intent(in) :: YLMSUN
      REAL    PHASETAB(MAXPHASE,MAXSCATANG), SINGSCAT(MAXPHASE)
cf2py intent(in) :: PHASETAB
      DOUBLE PRECISION SUNDIRLEG(0:MAXLEG)
      
      DOUBLE PRECISION X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2
      DOUBLE PRECISION COSSCAT
      DOUBLE PRECISION U, R, PI
      REAL    MU2, PHI2
      INTEGER N  
  
  
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
     .                       NX, NY, NZ, NPTS, NCELLS, 
     .                       GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .                       XGRID, YGRID, ZGRID, GRIDPOS,
     .                       ML, MM, NCS, NLM, NLEG, NUMPHASE,
     .                       NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, 
     .                       DELTAM, SRCTYPE, WAVELEN, SOLARMU,
     .                       SOLARAZ, EXTINCT, ALBEDO, LEGEN,  
     .                       IPHASE, DIRFLUX, SHPTR, SOURCE, 
     .                       YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,
     .                       MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, 
     .                       SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .                       MURAY, PHIRAY, MU2, PHI2, X0, Y0, Z0, 
     .                       VISOUT(N))
  900   CONTINUE
      ENDDO
  
      RETURN
      END
      
      
      SUBROUTINE BASE_GRID_PROJECTION (NBCELLS, NCELLS, NBPTS, 
     .                 GRIDPOS, GRIDPTR, TREEPTR, ARRAY, BGARRAY)

      IMPLICIT NONE
      INTEGER NBCELLS, NCELLS, NBPTS
Cf2py intent(in) :: NBCELLS, NCELLS, NBPTS
      REAL    GRIDPOS(3,*)
Cf2py intent(in) :: GRIDPOS
      INTEGER GRIDPTR(8,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, TREEPTR
      DOUBLE PRECISION  ARRAY(*)
Cf2py intent(in) :: ARRAY
      DOUBLE PRECISION BGARRAY(NBPTS)
Cf2py intent(out) :: BGARRAY
      
      DOUBLE PRECISION U, V, W, F(8), DELX, DELY, DELZ
      DOUBLE PRECISION INVDELX, INVDELY, INVDELZ
      INTEGER ICELL, BCELL, I, GPOINT, J, BGPOINT
      
      BGARRAY = ARRAY(1:NBPTS)
      DO ICELL = NBCELLS+1, NCELLS
C       Find base cell for icell
        BCELL = ICELL
        
        DO WHILE (TREEPTR(1, BCELL) .GT. 0)
          BCELL = TREEPTR(1, BCELL)
        ENDDO
        
        DELX = GRIDPOS(1,GRIDPTR(8,BCELL))-GRIDPOS(1,GRIDPTR(1,BCELL))
        DELY = GRIDPOS(2,GRIDPTR(8,BCELL))-GRIDPOS(2,GRIDPTR(1,BCELL))
        DELZ = GRIDPOS(3,GRIDPTR(8,BCELL))-GRIDPOS(3,GRIDPTR(1,BCELL))
      
C       loop over gridpoints belonging to icell and trilin interpolate
        DO I = 1,8
          GPOINT = GRIDPTR(I, ICELL)
          IF ((ARRAY(GPOINT) .NE. 0) .AND. (GPOINT .GT. NBPTS)) THEN
            U = 0.0
            V = 0.0
            W = 0.0
            IF (DELX .GT. 0.0) THEN 
              U = (GRIDPOS(1,GRIDPTR(8, BCELL))-GRIDPOS(1,GPOINT))/DELX
            ENDIF
            IF (DELY .GT. 0.0) THEN 
              V = (GRIDPOS(2,GRIDPTR(8, BCELL))-GRIDPOS(2,GPOINT))/DELY
            ENDIF
            IF (DELY .GT. 0.0) THEN 
              W = (GRIDPOS(3,GRIDPTR(8, BCELL))-GRIDPOS(3,GPOINT))/DELZ
            ENDIF
            F(8) = (1-U) * (1-V) * (1-W)
            F(7) =    U  * (1-V) * (1-W)
            F(6) = (1-U) *    V  * (1-W)
            F(5) =    U  *    V  * (1-W)
            F(4) = (1-U) * (1-V) * W
            F(3) =    U  * (1-V) * W
            F(2) = (1-U) *    V  * W
            F(1) =    U  *    V  * W
            
            DO J = 1,8 
               BGPOINT = GRIDPTR(J, BCELL)
               BGARRAY(BGPOINT) = BGARRAY(BGPOINT)+F(J)*ARRAY(GPOINT)
               ARRAY(GPOINT) = 0
            ENDDO
          ENDIF
        ENDDO
      ENDDO
  
      RETURN
      END