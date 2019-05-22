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
     .                   CAMPHI, NPIX, NPART, TOTAL_EXT, STOKES)
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
      REAL    TOTAL_EXT(NPTS), LEGEN(NSTLEG, 0:NLEG,*)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN, TOTAL_EXT

      REAL    DIRFLUX(*), FLUXES(2,*), SOURCE(NSTOKES, *)
Cf2py intent(in) :: DIRFLUX, FLUXES, SOURCE
      REAL CAMX(*), CAMY(*), CAMZ(*)
      DOUBLE PRECISION CAMMU(*), CAMPHI(*)
Cf2py intent(in) ::  CAMX, CAMY, CAMZ, CAMMU, CAMPHI
      INTEGER  NPIX, NSTPHASE, NSCATANGLE
Cf2py intent(in) :: NPIX
      REAL   STOKES(NSTOKES, NPIX)
Cf2py intent(out) :: STOKES
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS

      INTEGER I, J, L, SIDE
      INTEGER IVIS
      LOGICAL VALIDRAD
      DOUBLE PRECISION MURAY, PHIRAY, MU2, PHI2
      DOUBLE PRECISION U, R, PI
      DOUBLE PRECISION X0, Y0, Z0
      DOUBLE PRECISION XE,YE,ZE, TRANSMIT, VISRAD(NSTOKES)
      INTEGER MAXSCATANG
      PARAMETER (MAXSCATANG=721)
      REAL, ALLOCATABLE :: YLMSUN(:,:), PHASETAB(:,:,:)
      REAL  MEAN, STD1, STD2
  
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
     .   	             TOTAL_EXT, NPART)
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
