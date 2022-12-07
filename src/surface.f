
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
      CHARACTER SFCFILE*80, SFCTYPE*2
Cf2py intent(in) :: SFCTYPE
      INTEGER :: GRID_COORDS(2, NXSFC*NYSFC)
Cf2py intent(in) :: GRIDCOORDS
      REAL :: PARMS_IN(MAXSFCPARS, NXSFC*NYSFC)
Cf2py intent(in) ::PARMS_IN
      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG
      INTEGER COUNT
      INTEGER N, I, I0, IX, IY, J, Q
      REAL    ALB, TEMP, MRE, MIM, RHO0, KR, THETA, WSPD, PCL
      REAL    A, B, ZETA, SIGMA

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

      ELSE IF (SFCTYPE .EQ. 'VP') THEN
        NSFCPAR = MAXSFCPARS
        DO COUNT=1,NXSFC*NYSFC
          IX = GRID_COORDS(1,COUNT)
          IY = GRID_COORDS(2,COUNT)
          TEMP = PARMS_IN(1,COUNT)
          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
              I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
              SFCPARMS(I+1) = TEMP
              DO Q=2,NSFCPAR
                SFCPARMS(I+Q) = PARMS_IN(Q,COUNT)
              ENDDO
              N = N + 1
          ENDIF
        ENDDO

      ELSE
        IERR = 1
        WRITE(ERRMSG,*), 'PREP_SURFACE: Unknown BRDF type'
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
