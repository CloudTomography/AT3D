
      SUBROUTINE PREP_SURFACE (MAXSFCPTS, MAXSFCPARS,
     .                     SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC,
     .                     NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO,
     .                     GRID_COORDS, PARMS_IN, IERR, ERRMSG)
C    MODIFIED FROM READ_SURFACE (shdomsub2.f) by Jesse Loveridge. This
C    function now prepares the surface properties based on arrays passed
C    from python rather than read from file.

C       Reads the surface properties from the file into the surface arrays.
C     There are currently six surface types (Lambertian, Wave-Fresnel,
C     Ocean, RPV, and Diner polarized RPV), though it is simple to add more.
C     In all cases the surface properties are specified on a regular, evenly
C     spaced grid of size NXSFC by NYSFC and spacing DELXSFC by DELYSFC.
C     The properties are returned on a NXSFC+1 by NYSFC+1 grid that is
C     periodic (NX+1 same as 1, etc).  The surface properties at each
C     grid point are returned in SFCPARMS, and the first parameter must
C     be temperature.  The number of parameters for each grid point
C     (including temperature) are returned in NSFCPAR.
C     Also returned are the domain average ground albedo and temperature
C     (GNDALBEDO, GNDTEMP) for use with initialization (hence they can be
C     highly approximate).
C            Type        Parameters
C       L  Lambertian    albedo
C       W  Wave-Fresnel  Real, Imaginary index of refraction, wind speed (m/s)
C       O  Ocean         Wind speed (m/s), Pigmentation (mg/m^3)
C       D  Diner et al   a, k, b, zeta, sigma
C       R  RPV           rho0, k, Theta
C
      IMPLICIT NONE
      INTEGER MAXSFCPTS, MAXSFCPARS, NXSFC, NYSFC, NSFCPAR
Cf2py intent(in) :: MAXSFCPTS, MAXSFCPARS, NXSFC, NYSFC
Cf2py intent(out) :: NSFCPAR
      REAL    DELXSFC, DELYSFC
Cf2py intent(in) :: DELXSFC, DELYSFC
      REAL    SFCPARMS(MAXSFCPTS*MAXSFCPARS), GNDTEMP, GNDALBEDO
Cf2py intent(out) :: SFCPARMS, GNDTEMP, GNDALBEDO
      CHARACTER SFCTYPE*2
Cf2py intent(in) :: SFCTYPE
      INTEGER :: GRID_COORDS(2, NXSFC*NYSFC)
Cf2py intent(in) :: GRIDCOORDS
      REAL :: PARMS_IN(MAXSFCPARS, NXSFC*NYSFC)
Cf2py intent(in) ::PARMS_IN
      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG
      INTEGER COUNT
      INTEGER N, I, I0, IX, IY, J
      REAL    ALB, TEMP, MRE, MIM, RHO0, KR, THETA, WSPD, PCL
      REAL    A, B, ZETA, SIGMA, HB,BR,FISO,FGEO,FVOL

      IERR = 0
      GNDTEMP = 0.0
      GNDALBEDO = 0.0
      N = 0

C           Surface file type L is for variable Lambertian surface
      IF (SFCTYPE .EQ. 'VL' .OR. SFCTYPE .EQ. 'Vl') THEN
        NSFCPAR = 2
        IF (NSFCPAR .GT. MAXSFCPARS ) THEN
          IERR = 1
          WRITE (ERRMSG,*) 'PREP_SURFACE: MAXSFCPARS exceeded'
          RETURN
        ENDIF
        DO COUNT =1,NXSFC*NYSFC
          IX = GRID_COORDS(1,COUNT)
          IY = GRID_COORDS(2,COUNT)
          TEMP = PARMS_IN(1,COUNT)
          ALB = PARMS_IN(2,COUNT)

          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
            I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
            SFCPARMS(I+1) = TEMP
            SFCPARMS(I+2) = ALB
            IF (ALB .LT. 0.0 .OR. ALB .GT. 1.0) THEN
              IERR = 1
              WRITE(ERRMSG, *),
     .        'PREP_SURFACE: Illegal surface albedo'
              RETURN
            ENDIF
            IF (TEMP .LT. 150. .OR. TEMP .GT. 350.)
     .        WRITE (*,*) 'PREP_SURFACE: Warning -',
     .          ' surface temperature: ',IX,IY,TEMP
            GNDTEMP = GNDTEMP + TEMP
            GNDALBEDO = GNDALBEDO + ALB
            N = N + 1
          ENDIF
        ENDDO

C           Surface file type W is for Wave-Fresnel surface
      ELSE IF (SFCTYPE .EQ. 'VW') THEN
        NSFCPAR = 4
        IF (NSFCPAR .GT. MAXSFCPARS ) THEN
          IERR = 1
          WRITE (ERRMSG,*) 'PREP_SURFACE: MAXSFCPARS exceeded'
          RETURN
        ENDIF
        DO COUNT=1,NXSFC*NYSFC
          IX = GRID_COORDS(1,COUNT)
          IY = GRID_COORDS(2,COUNT)
          TEMP = PARMS_IN(1,COUNT)
          MRE = PARMS_IN(2,COUNT)
          MIM = PARMS_IN(3,COUNT)
          WSPD = PARMS_IN(4,COUNT)

          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
            I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
            SFCPARMS(I+1) = TEMP
            SFCPARMS(I+2) = MRE
            SFCPARMS(I+3) = MIM
            SFCPARMS(I+4) = WSPD
            GNDTEMP = GNDTEMP + TEMP
            N = N + 1
          ENDIF
        ENDDO

C           Surface file type O is for Ocean surface BRDF
C             Norm Loeb's modification of 6S ocean reflectance module.
      ELSE IF (SFCTYPE .EQ. 'VO') THEN
        NSFCPAR = 3
        IF (NSFCPAR .GT. MAXSFCPARS ) THEN
          IERR = 1
          WRITE (ERRMSG,*) 'READ_SURFACE: MAXSFCPARS exceeded'
          RETURN
        ENDIF
        DO COUNT=1,NXSFC*NYSFC
          IX = GRID_COORDS(1,COUNT)
          IY = GRID_COORDS(2,COUNT)
          TEMP = PARMS_IN(1,COUNT)
          WSPD = PARMS_IN(2,COUNT)
          PCL = PARMS_IN(3,COUNT)
          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
            I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
            SFCPARMS(I+1) = TEMP
            SFCPARMS(I+2) = WSPD
            SFCPARMS(I+3) = PCL
            GNDTEMP = GNDTEMP + TEMP
            N = N + 1
          ENDIF
        ENDDO

C           Surface file type R is for RPV surface:
C             (Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere
C              Reflectance (CSAR) Model. 2. Semiempirical Surface Model
C              Usable With NOAA Advanced Very High Resolution Radiometer
C              Data, J. Geophys. Res., 98, 20791-20801.)
      ELSE IF (SFCTYPE .EQ. 'VR') THEN
        NSFCPAR = 4
        IF (NSFCPAR .GT. MAXSFCPARS ) THEN
          IERR = 1
          WRITE (ERRMSG,*) 'PREP_SURFACE: MAXSFCPARS exceeded'
          RETURN
        ENDIF
        DO COUNT=1,NXSFC*NYSFC
          IX = GRID_COORDS(1,COUNT)
          IY = GRID_COORDS(2,COUNT)
          TEMP = PARMS_IN(1,COUNT)
          RHO0 = PARMS_IN(2,COUNT)
          KR = PARMS_IN(3,COUNT)
          THETA =PARMS_IN(4,COUNT)
          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
            I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
            SFCPARMS(I+1) = TEMP
            SFCPARMS(I+2) = RHO0
            SFCPARMS(I+3) = KR
            SFCPARMS(I+4) = THETA
            GNDTEMP = GNDTEMP + TEMP
            GNDALBEDO = GNDALBEDO + RHO0
            N = N + 1
          ENDIF
        ENDDO

C           Surface file type D is for Diner et al surface reflection:
C     (Diner, D. J., F. Xu, J. V. Martonchik, B. E. Rheingans, S. Geier,
C     V. M. Jovanovic, A. Davis, R. A. Chipman, S. C. McClain, 2012:
C     Exploration of a polarized surface bidirectional reflectance model
C     using the ground-based multiangle spectropolarimetric imager.
C     Atmosphere 2012, 3, 591-619; doi:10.3390/atmos3040591.)
      ELSE IF (SFCTYPE .EQ. 'VD') THEN
        NSFCPAR = 6
        IF (NSFCPAR .GT. MAXSFCPARS ) THEN
          IERR = 1
          WRITE (ERRMSG,*) 'PREP_SURFACE: MAXSFCPARS exceeded'
          RETURN
        ENDIF
        DO COUNT=1,NXSFC*NYSFC
          IX = GRID_COORDS(1,COUNT)
          IY = GRID_COORDS(2,COUNT)
          TEMP = PARMS_IN(1,COUNT)
          A = PARMS_IN(2,COUNT)
          KR = PARMS_IN(3,COUNT)
          B = PARMS_IN(4, COUNT)
          ZETA = PARMS_IN(5,COUNT)
          SIGMA = PARMS_IN(6,COUNT)
          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
            I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
            SFCPARMS(I+1) = TEMP
            SFCPARMS(I+2) = A
            SFCPARMS(I+3) = KR
            SFCPARMS(I+4) = B
            SFCPARMS(I+5) = ZETA
            SFCPARMS(I+6) = SIGMA
            GNDTEMP = GNDTEMP + TEMP
            GNDALBEDO = GNDALBEDO + A
            N = N + 1
          ENDIF
        ENDDO
      ELSE IF(SFCTYPE .EQ. 'VM') THEN
        NSFCPAR = 6
        IF (NSFCPAR .GT. MAXSFCPARS ) THEN
          IERR = 1
          WRITE (ERRMSG,*) 'PREP_SURFACE: MAXSFCPARS exceeded'
          RETURN
        ENDIF
        DO COUNT=1,NXSFC*NYSFC
          IX = GRID_COORDS(1,COUNT)
          IY = GRID_COORDS(2,COUNT)
          TEMP = PARMS_IN(1,COUNT)
          FISO = PARMS_IN(2,COUNT)
          FGEO = PARMS_IN(3,COUNT)
          FVOL = PARMS_IN(4, COUNT)
          HB = PARMS_IN(5,COUNT)
          BR = PARMS_IN(6,COUNT)
          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
            I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
            SFCPARMS(I+1) = TEMP
            SFCPARMS(I+2) = FISO
            SFCPARMS(I+3) = FGEO
            SFCPARMS(I+4) = FVOL
            SFCPARMS(I+5) = HB
            SFCPARMS(I+6) = BR
            GNDTEMP = GNDTEMP + TEMP
            GNDALBEDO = GNDALBEDO + FISO
            N = N + 1
          ENDIF
        ENDDO
      ELSE
        IERR = 1
        WRITE(ERRMSG,*) 'PREP_SURFACE: Unknown BRDF type'
        RETURN
      ENDIF


290   CONTINUE

C         Copy the edge points to the opposite side for periodic boundary
      DO IX = 1, NXSFC
        I0 = NSFCPAR*(IX-1 + (NXSFC+1)*(1-1))
        I  = NSFCPAR*(IX-1 + (NXSFC+1)*(NYSFC+1-1))
        DO J = 1, NSFCPAR
          SFCPARMS(I+J) = SFCPARMS(I0+J)
        ENDDO
      ENDDO
      DO IY = 1, NYSFC+1
        I0 = NSFCPAR*(        (NXSFC+1)*(IY-1))
        I  = NSFCPAR*(NXSFC + (NXSFC+1)*(IY-1))
        DO J = 1, NSFCPAR
          SFCPARMS(I+J) = SFCPARMS(I0+J)
        ENDDO
      ENDDO
      GNDALBEDO = GNDALBEDO/N
      GNDTEMP = GNDTEMP/N
      RETURN
      END

      SUBROUTINE VARIABLE_BRDF_SURFACE_GRAD (NBOTPTS, IBEG, IEND, BCPTR,
     .               NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,
     .               SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,
     .               SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES, BCRAD,
     .               SOLARREFLECT, SFCGRIDGRAD, NSFCDER, SFCDER,
     .               DPLANCK)
C       Computes the upwelling Stokes vector of the bottom boundary points
C     for one outgoing direction using the specified bidirectional
C     reflectance distribution function.  The upwelling Stokes vector
C     includes the reflection of the incident radiance, the thermal
C     emission (emissivity*Planck function), and the reflected direct
C     solar flux (if applicable).  The upwelling radiance vector is the
C     integral over all incident directions of the BRDF times the
C     downwelling radiance, so a discrete sum is done and the integral
C     weights (WTDO) are included. The general BRDF function is called
C     to compute the Stokes reflectance matrix for the particular BRDF
C     type (SFCTYPE) with parameters (SFCGRIDPARMS) for the incident and
C     outgoing directions.  The emissivity is computed implicitly from the
C     integral of the BRDF.  The outgoing direction is specified with
C     (MU2,PHI2), and the BRDF is computed for all incident directions
C     (loop over JMU,JPHI).  The incident downwelling radiance vectors are
C     input in BCRAD(*,*,2...) and the outgoing upwelling radiance vectors
C     are output in BCRAD(*,*,1).

C     Note that when this is called for the surface gradient IBEG=IEND
C     so there is no loop over NBOTPTS, for this reason SFCGRIDGRAD
C     and DPLANCK do not depend on NBOTPTS unlike the other variables.
      IMPLICIT NONE
      INTEGER NBOTPTS, IBEG, IEND, BCPTR(NBOTPTS)
      INTEGER NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR, NSTOKES
      REAL    WTDO(NMU,NPHI0MAX), MU(NMU), PHI(NMU,NPHI0MAX)
      REAL    MU2, PHI2, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
      REAL    SFCGRIDPARMS(NSFCPAR,NBOTPTS), BCRAD(NSTOKES,NBOTPTS,*)
      INTEGER NSFCDER, SFCDER(NSFCDER)
      REAL    SFCGRIDGRAD(NSTOKES,NSFCDER)
      REAL    DPLANCK
      REAL    SOLARREFLECT(NSTOKES)
      CHARACTER SRCTYPE*1, SFCTYPE*2

      INTEGER JMU, JPHI, IBC, I, JANG, K, K1, K2, Q
      REAL    OPI, REFLECT(4,4), W, SUM0, SUM1, RADSPEC(4)
      REAL    REFLECTGRAD(4,4,NSFCDER)

      OPI = 1.0/ACOS(-1.0)

      DO IBC = IBEG, IEND

C         Initialize the upwelling boundary radiances to zero or to
C           the reflection of direct unpolarized solar flux.
        BCRAD(:,IBC,1) = 0.0
        SOLARREFLECT(:) = 0.0
        SFCGRIDGRAD(:,:) = 0.0

        IF (SRCTYPE .NE. 'T') THEN
          I = BCPTR(IBC)
          CALL SURFACE_BRDF (SFCTYPE(2:2), SFCGRIDPARMS(2,IBC),WAVELEN,
     .                    MU2, PHI2, SOLARMU,SOLARAZ, NSTOKES, REFLECT)
          SOLARREFLECT(:) = OPI*REFLECT(1:NSTOKES,1)*DIRFLUX(I)
          BCRAD(:,IBC,1) = BCRAD(:,IBC,1)
     .                    + SOLARREFLECT(:)
          CALL SURFACE_BRDF_GRAD(SFCTYPE(2:2), SFCGRIDPARMS(2,IBC),
     .        WAVELEN,
     .        MU2, PHI2, SOLARMU, SOLARAZ, NSTOKES, REFLECTGRAD,
     .        NSFCDER, SFCDER)
          SFCGRIDGRAD(:,:) = OPI*REFLECTGRAD(1:NSTOKES,1,:)*DIRFLUX(I)
        ENDIF

C         Integrate over the incident discrete ordinate directions (JMU,JPHI)
        JANG = 1
        DO JMU = 1, NMU/2
          DO JPHI = 1, NPHI0(JMU)
            CALL SURFACE_BRDF (SFCTYPE(2:2), SFCGRIDPARMS(2,IBC),
     .                     WAVELEN, MU2,PHI2, MU(JMU),PHI(JMU,JPHI),
     .                     NSTOKES, REFLECT)
          CALL SURFACE_BRDF_GRAD(SFCTYPE(2:2), SFCGRIDPARMS(2,IBC),
     .          WAVELEN, MU2, PHI2, SOLARMU, SOLARAZ, NSTOKES,
     .          REFLECTGRAD, NSFCDER, SFCDER)
            W = OPI*ABS(MU(JMU))*WTDO(JMU,JPHI)
C             Add in the polarized reflection
            DO K1 = 1, NSTOKES
              BCRAD(:,IBC,1) = BCRAD(:,IBC,1)
     .                   + W*REFLECT(1:NSTOKES,K1)*BCRAD(K1,IBC,JANG+1)
              SFCGRIDGRAD(:,:) = SFCGRIDGRAD(:,:)
     .          + W*REFLECTGRAD(1:NSTOKES,K1,:)*BCRAD(K1,IBC,JANG+1)
            ENDDO
C             Add in the polarized thermal emission
            BCRAD(1,IBC,1) = BCRAD(1,IBC,1)
     .                     + W*(1-REFLECT(1,1))*SFCGRIDPARMS(1,IBC)
            SFCGRIDGRAD(1,:) = SFCGRIDGRAD(1,1)
     .          + W*(1-REFLECT(1,1))*DPLANCK
     .          - W*SFCGRIDPARMS(1,IBC)*REFLECTGRAD(1,1,:)
            BCRAD(2:,IBC,1) = BCRAD(2:,IBC,1)
     .                   - W*REFLECT(2:NSTOKES,1)*SFCGRIDPARMS(1,IBC)
            SFCGRIDGRAD(2:,:) = SFCGRIDGRAD(2:,:)
     .          - W*REFLECTGRAD(2:NSTOKES,1,:)*SFCGRIDPARMS(1,IBC)
            DO Q=1,NSFCDER
              SFCGRIDGRAD(2:,Q) = SFCGRIDGRAD(2:,Q)
     .          - W*REFLECT(2:NSTOKES,1)*DPLANCK
            ENDDO
            JANG = JANG + 1
          ENDDO
        ENDDO

      ENDDO
      RETURN
      END

      SUBROUTINE SURFACE_BRDF_GRAD (SFCTYPE, REFPARMS, WAVELEN,
     .                         MU2, PHI2, MU1, PHI1, NSTOKES,
     .                         REFLECTGRAD, NSFCDER, SFCDER)
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
      INTEGER NSTOKES, NSFCDER, SFCDER(NSFCDER)
      REAL    REFPARMS(*), WAVELEN, MU1, PHI1, MU2, PHI2
      REAL    REFLECTGRAD(4,4,NSFCDER)
      CHARACTER  SFCTYPE*1

      REAL   PI, KGEO,KVOL
      REAL   REFLECT

      INTEGER I, IPARAM

      PI = ACOS(-1.0)
      IF (SFCTYPE .EQ. 'L' .OR. SFCTYPE .EQ. 'l') THEN
C         L or l: Lambertian surface BRDF is constant.
C           (for testing, as done more efficiently by Lambertian routines)
        REFLECTGRAD(1:NSTOKES,1:NSTOKES,:) = 0.0
        DO I=1,NSFCDER
          IPARAM = SFCDER(I)
          REFLECTGRAD(1,1,I) = 1.0
        ENDDO
      ELSE IF (SFCTYPE .EQ. 'M') THEN
C       M:  RossThickLiSpare-R Surface BRDF kernel.
        CALL ROSS_THICK_LI_SPARSE(REFPARMS(1), REFPARMS(2), 
     .                              REFPARMS(3), 2.0,
     .                              1.0, -MU1,MU2,PHI1-PHI2,
     .                              REFLECT,KGEO,KVOL)
        REFLECTGRAD(1,1,1) = 1.0
        REFLECTGRAD(1,1,2) = KGEO
        REFLECTGRAD(1,1,3) = KVOL 
C     THIS IS WHERE YOU SHOULD ADD CALLS TO LINEARIZATIONS OF THE
C     DIFFERENT SURFACES IF THEY ARE AVAILABLE.
      ELSE
        PRINT *, 'SURFACE BRDF derivatives are not yet',
     .    ' supported for surface type', SFCTYPE
        STOP
      ENDIF




      RETURN
      END