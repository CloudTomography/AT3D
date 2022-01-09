!     This file contains fortran subroutines that perform mie calculations
!     written by Frank Evans.
!     https://nit.coloradolinux.com/~evans/shdom.html
!     Many of these subroutines have been modified for use in pyshdom by
!     Aviad Levis, Technion Institute of Technology, 2019 and
!     Jesse Loveridge, University of Illinois at Urbana-Champaign, 2020-2021.
!     - JRLoveridge 2021/02/22


SUBROUTINE GET_MIE_TABLE (NRETAB, MAXLEG, WAVELEN1, WAVELEN2, WAVELENCEN, DELTAWAVE, &
                          PARDENS, SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS, RINDEX, PARTYPE, &
                          AVGFLAG, DISTFLAG, REFF, EXTINCT,SSALB,NLEG,LEGCOEF,LOGRE, &
                          VERBOSE, IERR, ERRMSG)
!
! Does Mie computations to create a scattering table as a function of
! effective radius for gamma, modified gamma, or lognormal size distributions
! of spherical particles.  The particles may be water or ice (in which case
! the program provides the index of refraction depending on wavelength) or
! "aerosols" (in which case the index of refraction is user specified).
! For water or ice particles the scattering properties may be averaged
! over the desired spectral range with Planck function weighting.
! The phase functions in the output scattering table are represented
! with Legendre series.  For polarized output, the six unique elements of
! the phase matrix are represented with Wigner d-function expansions
! (Doicu et al., 2013, JQSRT, http://dx.doi.org/10.1016/j.jqsrt.2012.12.009).
!
!  compile: pgf90 -fast -o make_mie_table  make_mie_table.f90
!                 indexwatice.f  miewig.f
!
!    Frank Evans    University of Colorado       May 2003
!  Modified for polarization                     March 2012
!  Modified for modified gamma distribution      September 2012
!  Modified for Wigner d-function coefficients   July 2012
  IMPLICIT NONE
  INTEGER :: NRETAB
  LOGICAL :: LOGRE, VERBOSE
!  f2py intent(in) :: NRETAB, LOGRE, VERBOSE
  REAL :: WAVELEN1, WAVELEN2, DELTAWAVE, PARDENS
!  f2py intent(in) :: WAVELEN1, WAVELEN2, DELTAWAVE, PARDENS
  REAL :: SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS
!  f2py intent(in) :: SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS
  COMPLEX :: RINDEX
  CHARACTER(LEN=1), INTENT(IN) :: PARTYPE, AVGFLAG, DISTFLAG
!  f2py intent(in) :: PARTYPE, AVGFLAG, DISTFLAG
  INTEGER :: NSIZE, I, J, L, NL
  INTEGER, INTENT(IN) :: MAXLEG
!  f2py intent(in) :: MAXLEG
  REAL :: SCATTER, WAVELENCEN
!  f2py intent(in) :: WAVELENCEN
  REAL, INTENT(OUT) :: REFF(NRETAB), EXTINCT(NRETAB), SSALB(NRETAB)
!  f2py intent(out) :: REFF, EXTINCT, SSALB
  REAL, INTENT(OUT) :: NLEG(NRETAB), LEGCOEF(6,0:MAXLEG,NRETAB)
!  f2py intent(out) :: NLEG, LEGCOEF'
  INTEGER :: IERR
  CHARACTER(LEN=600) :: ERRMSG
  INTEGER, ALLOCATABLE :: NLEG1(:)
  REAL, ALLOCATABLE :: RADII(:), ND(:)
  REAL, ALLOCATABLE :: QEXT(:), QSCA(:)
  REAL, ALLOCATABLE :: EXTINCT1(:), SCATTER1(:), LEGCOEF1(:, :,:)

  IF (PARTYPE == 'W') THEN
    PARDENS = 1.0
  ELSE IF (PARTYPE == 'I') THEN
    PARDENS = 0.916
  ENDIF

   ! Get the average index of refraction for water or ice
  IF (PARTYPE /= 'A') THEN
    CALL GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2, RINDEX)
  ENDIF

   ! Figure the number of radii there will be
  CALL GET_NSIZE (SRETAB, MAXRADIUS, WAVELENCEN, NSIZE)

   ! Allocate all the arrays here
  ALLOCATE (RADII(NSIZE), ND(NSIZE), NLEG1(NSIZE))
  ALLOCATE (EXTINCT1(NSIZE), SCATTER1(NSIZE), LEGCOEF1(6, 0:MAXLEG,NSIZE))


   ! Make up the discrete particle radii to use
  CALL GET_SIZES (SRETAB, MAXRADIUS, WAVELENCEN, NSIZE, RADII)

   ! Do the Mie computations for each radius, which may involve several
   !   Mie calculation over the wavelength integration

  CALL COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, PARTYPE, &
                              WAVELENCEN, RINDEX, NSIZE, RADII, MAXLEG, &
                              EXTINCT1, SCATTER1, NLEG1, LEGCOEF1, VERBOSE, &
                              IERR, ERRMSG)


  ! Loop over the number of output tabulated effective radii
  DO I = 1, NRETAB
    ! Set tabulated effective radius
    IF (NRETAB <= 1) THEN
      REFF(I) = SRETAB
    ELSE IF (LOGRE) THEN
      REFF(I) = EXP((LOG(ERETAB)-LOG(SRETAB))*FLOAT(I-1)/(NRETAB-1) +LOG(SRETAB))
    ELSE
      REFF(I) = (ERETAB-SRETAB)*FLOAT(I-1)/(NRETAB-1) + SRETAB
    ENDIF
    ! Calculate the discrete size number concentrations (ND), which vary
    ! according to a truncated gamma, modified gamma, or lognormal
    ! distribution that gives the desired effective radius (REFF) and LWC (1 g/m^3).
    CALL MAKE_SIZE_DIST (DISTFLAG, PARDENS, NSIZE, RADII, REFF(I), ALPHA, GAMMA, &
                         ND, IERR, ERRMSG)

    ! Sum the scattering properties over the discrete size distribution
    EXTINCT(I) = 0.0
    SCATTER = 0.0
    LEGCOEF(:,:,I) = 0.0
    NL = 1
    DO J = 1, NSIZE
      EXTINCT(I) = EXTINCT(I) + ND(J)*EXTINCT1(J)
      SCATTER = SCATTER + ND(J)*SCATTER1(J)
      NL = MAX(NL,NLEG1(J))
      LEGCOEF(:,0:NL,I) = LEGCOEF(:,0:NL,I) + ND(J)*LEGCOEF1(:,0:NL,J)
    ENDDO
    DO L = 0, NL
      LEGCOEF(:,L,I) = LEGCOEF(:,L,I)/SCATTER
      IF (LEGCOEF(1,L,I) .GT. 0.5E-5) NLEG(I) = L
    ENDDO
    IF (ABS(LEGCOEF(1,0,I)-1.0) > 0.0001) THEN
      PRINT *,'Phase function not normalized for Reff=',REFF,LEGCOEF(1,0,I)
      STOP
    ENDIF
    IF (EXTINCT(I) > 0.0) THEN
      SSALB(I) = SCATTER/EXTINCT(I)
    ENDIF
    EXTINCT(I) = 0.001*EXTINCT(I)

  ENDDO  ! end of effective radius loop

END SUBROUTINE GET_MIE_TABLE





SUBROUTINE USER_INPUT (POLTAB, WAVELEN1,WAVELEN2, PARTYPE, RINDEX, PARDENS, &
                       AVGFLAG, DELTAWAVE, DISTFLAG, ALPHA, GAMMA, &
                       NRETAB, SRETAB, ERETAB, LOGRE, MAXRADIUS,  MIETABFILE)
 ! Reads the input parameters from the standard input
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: NRETAB
  LOGICAL, INTENT(OUT) :: POLTAB, LOGRE
  REAL,    INTENT(OUT) :: WAVELEN1, WAVELEN2, DELTAWAVE, PARDENS
  REAL,    INTENT(OUT) :: SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS
  COMPLEX, INTENT(OUT) :: RINDEX
  CHARACTER(LEN=1), INTENT(OUT) :: PARTYPE, AVGFLAG, DISTFLAG
  CHARACTER(LEN=*), INTENT(OUT) :: MIETABFILE

  WRITE(*,*) 'Making Mie scattering tables for spherical particles'

  WRITE (*,*) 'Make polarized Mie table (T or F)'
  READ (*,*) POLTAB
    WRITE (*,'(L)') POLTAB

  WRITE (*,*) 'Wavelength range (micron)'
  READ (*,*) WAVELEN1, WAVELEN2
    WRITE (*,'(2(1X,F9.3))') WAVELEN1, WAVELEN2
  IF (WAVELEN1 > WAVELEN2) STOP 'USER_INPUT: wavelength1 must be <= wavelength2'

  WRITE(*,*) 'Water, Ice, or Aerosol spherical particles (W,I,A)'
  READ(*,'(A1)') PARTYPE
    WRITE(*,'(1X,A1)') PARTYPE

  IF (PARTYPE == 'W' .OR. PARTYPE == 'I') THEN
    WRITE(*,*) 'Average Mie properties over wavelength or use Planck weighted'
    WRITE(*,*) '  index of refraction at center wavelength (A or C)'
    READ(*,'(A1)') AVGFLAG
    WRITE (*,'(1X,A1)') AVGFLAG
    IF (AVGFLAG == 'A') THEN
      WRITE(*,*) 'Wavelength interval for averaging (micron)'
      READ(*,*) DELTAWAVE
      WRITE(*,'(1X,F6.3)') DELTAWAVE
    ELSE
      AVGFLAG = 'C'
      DELTAWAVE = -1.0
    ENDIF
  ELSE
    PARTYPE = 'A'
    WRITE (*,*) 'Aerosol complex index of refraction (negative imaginary part)'
    READ (*,*) RINDEX
    WRITE (*,*) RINDEX
    WRITE (*,*) 'Aerosol particle bulk density (g/cm^3)'
    READ (*,*) PARDENS
    WRITE (*,'(1X,F5.3)') PARDENS
    AVGFLAG = 'C'
  ENDIF

  WRITE (*,*) 'Particle size distribution type: G = Gamma, M = modified gamma, L = Lognormal'
  READ (*,*) DISTFLAG
  WRITE (*,*) DISTFLAG
  IF (DISTFLAG == 'L') THEN
    WRITE (*,*) 'Log normal size distribution log standard deviation'
    READ (*,*) ALPHA
    WRITE (*,'(1X,F6.3)') ALPHA
  ELSE IF (DISTFLAG == 'G') THEN
    WRITE(*,*) 'Gamma size distribution shape parameter (alpha)'
    READ (*,*) ALPHA
    WRITE (*,'(1X,F6.3)') ALPHA
  ELSE IF (DISTFLAG == 'M') THEN
    WRITE(*,*) 'Modified gamma size distribution shape parameters (alpha & gamma)'
    READ (*,*) ALPHA, GAMMA
    WRITE (*,'(2(1X,F6.3))') ALPHA, GAMMA
  ELSE
    WRITE (*,*) 'Unrecognized size distribution type'
    STOP
  ENDIF

  WRITE(*,*) 'Number, starting, and ending tabulated effective radius (micron)'
  READ(*,*) NRETAB, SRETAB, ERETAB
    WRITE (*,'(1X,I3,2(1X,F7.2))') NRETAB, SRETAB, ERETAB

  WRITE(*,*) 'Log-spaced effective radius (T or F) (F for evenly spaced)'
  READ(*,*) LOGRE
    WRITE (*,'(L)') LOGRE

  WRITE(*,*) 'Maximum particle radius in size distribution (micron)'
  READ(*,*) MAXRADIUS
    WRITE (*,'(1X,F7.2)') MAXRADIUS

  WRITE (*,*) 'Output Mie scattering table name'
  READ (*,*) MIETABFILE
    WRITE(*,'(1X,A70)') MIETABFILE
END SUBROUTINE USER_INPUT



SUBROUTINE GET_NSIZE (SRETAB, MAXRADIUS, WAVELEN, NSIZE)
 ! Calculates the number of radii for which the Mie computation will be run.
 ! The formula and spacing in size parameter can be changed to trade
 ! off size distribution integration accuracy vs. computer time.
  IMPLICIT NONE
  REAL,    INTENT(IN)  :: SRETAB, MAXRADIUS, WAVELEN
  INTEGER, INTENT(OUT) :: NSIZE
  REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD, DELRADMAX

  TWOPI = 2.0*ACOS(-1.0)
  RADMIN = 0.02*SRETAB
  RAD = RADMIN
  DELRADMAX = 0.002*MAXRADIUS
  NSIZE = 1
  DO WHILE (RAD < MAXRADIUS)
    X = TWOPI*RAD/WAVELEN
    DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
    DELRAD = MIN(DELX*WAVELEN/TWOPI,DELRADMAX)
    RAD = RAD + DELRAD
    NSIZE = NSIZE + 1
  ENDDO
END SUBROUTINE GET_NSIZE


SUBROUTINE GET_SIZES (SRETAB, MAXRADIUS, WAVELEN, NSIZE, RADII)
 ! Calculates the radii for which the Mie computation will be run and
 ! from which all the size distributions will be computed.
 ! The formula and spacing in size parameter can be changed to trade
 ! off size distribution integration accuracy vs. computer time.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: SRETAB, MAXRADIUS, WAVELEN
  REAL,    INTENT(OUT) :: RADII(NSIZE)
  INTEGER :: N
  REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD, DELRADMAX

  TWOPI = 2.0*ACOS(-1.0)
  RAD = 0.02*SRETAB
  DELRADMAX = 0.002*MAXRADIUS
  RADII(1) = RAD
  DO N = 2, NSIZE
    X = TWOPI*RAD/WAVELEN
    DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
    DELRAD = MIN(DELX*WAVELEN/TWOPI,DELRADMAX)
    RAD = RAD + DELRAD
    RADII(N) = RAD
  ENDDO
END SUBROUTINE GET_SIZES



SUBROUTINE GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
!  Returns the Planck weighted center wavelength averaged over the
! wavelength interval (WAVELEN1 < WAVELEN2 [microns]).  A solar
! blackbody temperature of 5800 K is used for the Planck weighting
! if the average wavelength is less than 3 microns, no Planck weighting
! is done for an average wavelength between 3 and 5 microns, and a
! blackbody temperature of 270 K is done for an average wavelength
! greater than 5 microns.
  IMPLICIT NONE
  REAL, INTENT(IN)  :: WAVELEN1, WAVELEN2
!  f2py intent(in) :: WAVELEN1, WAVELEN2
  REAL, INTENT(OUT) :: WAVELENCEN
!  f2py intent(out) :: WAVELENCEN
  REAL :: WAVECEN, DELWAVE, WAVE, BBTEMP, PLANCK, SUMP, SUMW

  IF (WAVELEN1 == WAVELEN2) THEN
    WAVELENCEN = WAVELEN1
  ELSE
    WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
    IF (WAVECEN < 3.0) THEN
      BBTEMP = 5800.0
    ELSE IF (WAVECEN > 5.0) THEN
      BBTEMP = 270.0
    ELSE
      BBTEMP = -1.0
      PLANCK = 1.0
    ENDIF
    DELWAVE = MIN(WAVECEN/100.,0.1*ABS(WAVELEN2-WAVELEN1))
    SUMP = 0.0
    SUMW = 0.0
    WAVE = WAVELEN1
    DO WHILE (WAVE <= WAVELEN2)
      IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
      SUMP = SUMP + PLANCK
      SUMW = SUMW + PLANCK*WAVE
      WAVE = WAVE + DELWAVE
    ENDDO
    WAVELENCEN = 0.001*NINT(1000*SUMW/SUMP)
  ENDIF
END SUBROUTINE GET_CENTER_WAVELEN




SUBROUTINE GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2, RINDEX)
 ! Returns the index of refraction for water or ice averaged over
 ! the wavelength interval (WAVELEN1 < WAVELEN2 [microns]).   The
 ! averaging is done at 0.05 micron intervals and is weighted by
 ! a Planck function at 5800 K temperature for central wavelengths
 ! less than 3 microns, a flat function between 3 and 5 microns, and
 ! 270 K Planck function beyond 5 microns.
 ! The index of refraction is at -30 C for ice and +10 C for water
 ! (the temperature dependence is important in the microwave).
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN) :: PARTYPE
  REAL, INTENT(IN) :: WAVELEN1, WAVELEN2
  COMPLEX, INTENT(OUT) :: RINDEX
  REAL :: WAVECEN, WAVECUT, DELWAVE, WAVE, BBTEMP, PLANCK
  REAL :: MRE, MIM, SUMP, SUMMR, SUMMI, A

  WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
  IF (WAVECEN < 3.0) THEN
    BBTEMP = 5800.0
  ELSE IF (WAVECEN > 5.0) THEN
    BBTEMP = 270.0
  ELSE
    BBTEMP = -1.0
    PLANCK = 1.0
  ENDIF
  DELWAVE = MIN(WAVECEN/100.,0.1*ABS(WAVELEN2-WAVELEN1))
  DELWAVE = MAX(DELWAVE,WAVECEN*1.0E-5)
  SUMP = 0.0
  SUMMR = 0.0
  SUMMI = 0.0
  WAVE = WAVELEN1
  DO WHILE (WAVE <= WAVELEN2)
    IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
    SUMP = SUMP + PLANCK
    IF (PARTYPE == 'I') THEN
      CALL REFICE (0, WAVE, 243.0, MRE, MIM, A, A)
    ELSE
      CALL REFWAT (0, WAVE, 283.0, MRE, MIM, A, A)
    ENDIF
    SUMMR = SUMMR + PLANCK*MRE
    SUMMI = SUMMI + PLANCK*MIM
    WAVE = WAVE + DELWAVE
  ENDDO
  MRE = SUMMR/SUMP
  MIM = SUMMI/SUMP
  RINDEX = CMPLX(MRE,-MIM)
END SUBROUTINE GET_REFRACT_INDEX




SUBROUTINE COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, &
                                PARTYPE, WAVELENCEN, RINDEX, NSIZE, RADII, &
                                MAXLEG, EXTINCT1, SCATTER1, NLEG1, LEGCOEF1,&
				                        VERBOSE, IERR, ERRMSG)
 ! Does a Mie computation for each particle radius in RADII and returns the
 ! optical properties in arrays EXTINCT1, SCATTER1, NLEG1, and LEGCOEF1.
 ! For AVGFLAG='C' the computation is done at a single wavelength (WAVELENCEN),
 ! using the input index of refraction (RINDEX).  For AVGFLAG='A' an
 ! integration of the Mie properties over wavelength is performed for
 ! each radius.  For each wavelength, with spacing DELTAWAVE, the water
 ! or ice (depending on PARTYPE) index of refraction is obtained and
 ! used in the Mie computation for that wavelength, and the Mie optical
 ! properties are averaged with Planck function weighting (blackbody
 ! temperature depends on wavelength).  The Wigner d coefficients are
 ! returned with the product of the phase function times the scattering
 ! coefficient.
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: VERBOSE
  INTEGER, INTENT(IN) :: NSIZE, MAXLEG
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE, WAVELENCEN
  REAL,    INTENT(IN) :: RADII(NSIZE)
  COMPLEX, INTENT(IN) :: RINDEX
  CHARACTER(LEN=1), INTENT(IN) :: AVGFLAG, PARTYPE
  INTEGER, INTENT(OUT) :: NLEG1(NSIZE)
  REAL,    INTENT(OUT) :: EXTINCT1(NSIZE), SCATTER1(NSIZE)
  REAL,    INTENT(OUT) :: LEGCOEF1(6,0:MAXLEG,NSIZE)
  INTEGER,  INTENT(OUT) :: IERR
  CHARACTER(LEN=600), INTENT(OUT) :: ERRMSG
!f2py intent(out) :: ERRMSG, IERR, EXTINCT1, SCATTER1, LEGCOEF1
!f2py intent(out) :: NLEG1
  CHARACTER(LEN=6) :: TABLE_TYPE
  INTEGER :: I, NL
  REAL    :: WAVECEN, WAVE, BBTEMP, PLANCK, SUMP, A
  REAL    :: MRE, MIM, EXT, SCAT, COEF(6,0:MAXLEG)
  COMPLEX :: REFIND
  IERR = 0
  IF (VERBOSE) THEN
    WRITE(*,*) 'Computing mie scattering for all sizes, &
      this may take a while...'
  ENDIF
  TABLE_TYPE = 'VECTOR'
  IF (AVGFLAG == 'C') THEN
     ! For using one central wavelength: just call Mie routine for each radius
    DO I = 1, NSIZE
      IF (VERBOSE) THEN
        WRITE(*,*) 'Computing mie for radius: ', RADII(I), ' microns'
      ENDIF
      CALL MIE_ONE (WAVELENCEN, RINDEX, RADII(I), MAXLEG, &
                    EXTINCT1(I), SCATTER1(I), NLEG1(I), LEGCOEF1(1,0,I),&
                     IERR, ERRMSG)
      IF (IERR .NE. 0) RETURN
    ENDDO

  ELSE
     ! For averaging over wavelength range:
    WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
    IF (WAVECEN < 3.0) THEN
      BBTEMP = 5800.0
    ELSE IF (WAVECEN > 5.0) THEN
      BBTEMP = 270.0
    ELSE
      BBTEMP = -1.0
      PLANCK = 1.0
    ENDIF
    EXTINCT1(:) = 0.0
    SCATTER1(:) = 0.0
    NLEG1(:) = 1
    LEGCOEF1(:,:,:) = 0.0
    SUMP = 0.0
    WAVE = WAVELEN1
    DO WHILE (WAVE <= WAVELEN2)   ! Loop over the wavelengths
      IF (VERBOSE) THEN
        WRITE(*,*) 'Computing mie for wavelength: ', WAVE, ' microns'
      ENDIF
      IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
      SUMP = SUMP + PLANCK
      IF (PARTYPE == 'I') THEN   ! Get the index of refraction of water or ice
        CALL REFICE (0, WAVE, 243.0, MRE, MIM, A, A)
      ELSE
        CALL REFWAT (0, WAVE, 283.0, MRE, MIM, A, A)
      ENDIF
      REFIND = CMPLX(MRE,-MIM)
      DO I = 1, NSIZE
        IF (VERBOSE) THEN
          WRITE(*,*) 'Computing mie for radius: ', RADII(I), &
                  ' microns [wavelength = ', WAVE, 'microns]'
        ENDIF
        CALL MIE_ONE (WAVE, REFIND, RADII(I), MAXLEG, EXT, SCAT, NL, COEF, &
                IERR, ERRMSG)
        IF (IERR .NE. 0) RETURN
        EXTINCT1(I) = EXTINCT1(I) + PLANCK*EXT
        SCATTER1(I) = SCATTER1(I) + PLANCK*SCAT
        NLEG1(I) = MAX(NLEG1(I),NL)
        LEGCOEF1(:,0:NL,I) = LEGCOEF1(:,0:NL,I) + PLANCK*COEF(:,0:NL)
      ENDDO
      WAVE = WAVE + DELTAWAVE
    ENDDO
    EXTINCT1(:) = EXTINCT1(:)/SUMP
    SCATTER1(:) = SCATTER1(:)/SUMP
    LEGCOEF1(:,:,:) = LEGCOEF1(:,:,:)/SUMP
  ENDIF
END SUBROUTINE COMPUTE_MIE_ALL_SIZES



SUBROUTINE MAKE_MULTI_SIZE_DIST(DISTFLAG, PARDENS, NSIZE, RADII, REFF, ALPHA,&
				GAMMA, ND, NDIST, IERR, ERRMSG)
 ! Calculates the number concentrations (ND in cm^-3) for the NSIZE
 ! discrete particle radii (micron) of a gamma or lognormal size distribution
 ! with an effective radius of REFF (micron), gamma shape parameter or
 ! lognormal standard deviation of ALPHA, and mass content of 1 g/m^3.
 ! Also handles modified gamma distributions specified by ALPHA and GAMMA.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NSIZE, NDIST
  REAL,    INTENT(IN)  :: RADII(NSIZE), REFF(NDIST), ALPHA(NDIST), GAMMA, PARDENS
  REAL,    INTENT(OUT) :: ND(NSIZE, NDIST)
  CHARACTER(LEN=600), INTENT(OUT) :: ERRMSG
  INTEGER, INTENT(OUT) :: IERR
  CHARACTER(LEN=1), INTENT(IN) :: DISTFLAG
  INTEGER :: K
  IERR = 0
  DO K = 1, NDIST
    CALL MAKE_SIZE_DIST(DISTFLAG, PARDENS, NSIZE, RADII, REFF(K), ALPHA(K),&
   			GAMMA, ND(1,K), IERR, ERRMSG)
    IF (IERR .NE. 0) RETURN
  ENDDO
END SUBROUTINE MAKE_MULTI_SIZE_DIST




SUBROUTINE MAKE_SIZE_DIST(DISTFLAG, PARDENS, NSIZE, RADII, REFF, ALPHA, GAMMA, ND,&
                          IERR, ERRMSG)
 ! Calculates the number concentrations (ND in cm^-3) for the NSIZE
 ! discrete particle radii (micron) of a gamma or lognormal size distribution
 ! with an effective radius of REFF (micron), gamma shape parameter or
 ! lognormal standard deviation of ALPHA, and mass content of 1 g/m^3.
 ! Also handles modified gamma distributions specified by ALPHA and GAMMA.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NSIZE
  REAL,    INTENT(IN)  :: RADII(NSIZE), REFF, ALPHA, GAMMA, PARDENS
  REAL,    INTENT(OUT) :: ND(NSIZE)
  CHARACTER(LEN=1), INTENT(IN) :: DISTFLAG
  CHARACTER(LEN=600) :: ERRMSG
  INTEGER :: IERR
  REAL, PARAMETER :: TOL=0.001  ! fractional tolerance in achieving Reff
  INTEGER :: I
  REAL    :: TRUERE, F, REHI, RELO, REMID

   ! See if the true effective radius is already close enough
  CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, GAMMA, REFF, &
                     NSIZE, RADII, ND, TRUERE)
  IF (ABS(TRUERE-REFF) < TOL*REFF) RETURN
  F = REFF/TRUERE
  IF (TRUERE < REFF) THEN
    ! Find Reff that gives true Reff above desired value
    RELO = REFF
    REHI = F*REFF
    I = 0
    TRUERE = REFF/F
    DO WHILE (TRUERE <= REFF .AND. I < 8)
      REHI = F*REHI
      I = I + 1
      CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, GAMMA, REHI, &
                         NSIZE, RADII, ND, TRUERE)
    ENDDO
    IF (TRUERE <= REFF) THEN
      IERR = 1
      WRITE(ERRMSG, *) 'MAKE_SIZE_DIST: effective radius cannot be achieved', &
          REFF,TRUERE,ALPHA
      RETURN
    ENDIF
  ELSE
    ! Find Reff that gives true Reff below desired value
    REHI = REFF
    RELO = F*REFF
    I = 0
    TRUERE = REFF/F
    DO WHILE (TRUERE >= REFF .AND. I < 8)
      RELO = F*RELO
      I = I + 1
      CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, GAMMA, RELO, &
                         NSIZE, RADII, ND, TRUERE)
    ENDDO
    IF (TRUERE >= REFF) THEN
      IERR = 1
      WRITE(ERRMSG, *) 'MAKE_SIZE_DIST: effective radius cannot be achieved', &
          REFF,TRUERE,ALPHA
      RETURN
      ! PRINT *, 'MAKE_SIZE_DIST: effective radius cannot be achieved',REFF,TRUERE,ALPHA
      ! STOP
    ENDIF
  ENDIF
  ! Do bisection to get correct effective radius
  DO WHILE (ABS(TRUERE-REFF) > TOL*REFF)
    REMID = 0.5*(RELO+REHI)
    CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, GAMMA, REMID, &
                       NSIZE, RADII, ND, TRUERE)
    IF (TRUERE < REFF) THEN
      RELO = REMID
    ELSE
      REHI = REMID
    ENDIF
  ENDDO
END SUBROUTINE MAKE_SIZE_DIST



SUBROUTINE DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, GAMMA, RE, NSIZE, RADII, &
                         ND, TRUERE)
 ! For the input effective radius (RE) [um], returns the number concentrations
 ! ND [cm^-3] and the calculated effective radius TRUERE [um] for a
 ! gamma, modified gamma, or lognormal size distribution with mass content
 ! of 1 g/m^3.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: PARDENS, ALPHA, GAMMA, RE, RADII(NSIZE)
  CHARACTER(LEN=1), INTENT(IN) :: DISTFLAG
  REAL,    INTENT(OUT) :: ND(NSIZE), TRUERE
  INTEGER :: J
  REAL :: GAMMLN
  DOUBLE PRECISION :: DENS, PI, A, B, LWC, R, DELR, SUM2, SUM3, DGAMMA
  DOUBLE PRECISION :: BDERIV, ADERIV

  PI = ACOS(-1.0)
  IF (DISTFLAG == 'G') THEN
    B = (ALPHA+3)/RE
    A = 1.E6/( (4*PI/3.)*PARDENS*B**(-ALPHA-4)*DGAMMA( DBLE(ALPHA)+4.0D0) )
  ELSE IF (DISTFLAG == 'M') THEN
    B = (EXP(GAMMLN((ALPHA+4.)/GAMMA)-GAMMLN((ALPHA+3.)/GAMMA)) /RE)**GAMMA
    A = 1.E6*GAMMA *B**((ALPHA+4)/GAMMA) &
        /( (4*PI/3.)*PARDENS *EXP(GAMMLN((ALPHA+4.)/GAMMA)) )
  ELSE IF (DISTFLAG == 'L') THEN
    B = RE*EXP(-2.5*ALPHA**2)
    A = 1.E6/( (4*PI/3.)*PARDENS *SQRT(2*PI)*ALPHA * B**3 *EXP(4.5*ALPHA**2) )
! A test of the derivative of the lognormal distribution with respect to Reff
  ELSE IF (DISTFLAG == 'D') THEN
    B = RE*EXP(-2.5*ALPHA**2)
    BDERIV = EXP(-2.5*ALPHA**2)
    A = 1.E6/( (4*PI/3.)*PARDENS *SQRT(2*PI)*ALPHA * B**3 *EXP(4.5*ALPHA**2) )
    ADERIV = -3*A/B * BDERIV
  ENDIF
  LWC = 0.0
  SUM2 = 0.0
  SUM3 = 0.0
  DO J = 1, NSIZE
    R = RADII(J)
    DELR = SQRT(RADII(J)*RADII(MIN(NSIZE,J+1))) &
         - SQRT(RADII(J)*RADII(MAX(1,J-1)))
    IF (DISTFLAG == 'G') THEN
      ND(J) = A* R**ALPHA *EXP(-B*R) *DELR
    ELSE IF (DISTFLAG == 'M') THEN
      ND(J) = A* R**ALPHA *EXP(-B*R**GAMMA) *DELR
    ELSE IF (DISTFLAG == 'L') THEN
      ND(J) = (A/R)*EXP(-0.5*(LOG(R/B))**2/ALPHA**2) *DELR
    ELSE IF (DISTFLAG == 'D') THEN
      ND(J) = DELR*( (ADERIV/R)*EXP(-0.5*(LOG(R/B))**2/ALPHA**2) + &
              (A*LOG(R/B)*EXP(-0.5*(LOG(R/B))**2/ALPHA**2)/(R*B*ALPHA**2))*BDERIV )
    ENDIF
    IF (DISTFLAG == 'D') THEN
      LWC = LWC + 1.0E-6*PARDENS*(A/R)*EXP(-0.5*(LOG(R/B))**2/ALPHA**2)*DELR*(4*PI/3)*R**3
    ELSE
      LWC = LWC + 1.0E-6*PARDENS*ND(J)*(4*PI/3)*R**3
    ENDIF
    SUM2 = SUM2 + ND(J)*R**2
    SUM3 = SUM3 + ND(J)*R**3
  ENDDO
  ND(:) = (1.0/LWC)*ND(:)
  TRUERE = SUM3/SUM2
END SUBROUTINE DO_SIZE_DIST




SUBROUTINE WRITE_POLY_TABLE (MIETABFILE, WAVELEN1, WAVELEN2, DELTAWAVE, &
                            PARTYPE, PARDENS, RINDEX, &
                            DISTFLAG, ALPHA, GAMMA, NRETAB,&
			    NVETAB, REFF, VEFF, EXTINCT, SSALB,&
			    NLEG, LEGCOEF, NDIST, MAXLEG)
 ! Writes the table of Mie scattering properties as a function of
 ! effective radius.  There are two types of tables output:
 ! unpolarized with Legendre (Wigner d_l_00) series of a single phase
 ! matrix element (P11) or polarized with six Wigner d-function series
 ! of the Mie phase matrix elements.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NRETAB, NVETAB, NDIST
  INTEGER, INTENT(IN) :: MAXLEG, NLEG(NDIST)
  COMPLEX, INTENT(IN) :: RINDEX
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE
  REAL,    INTENT(IN) :: PARDENS, ALPHA(NVETAB), GAMMA
  REAL,    INTENT(IN) :: REFF(NRETAB), VEFF(NVETAB), EXTINCT(NDIST)
  REAL,    INTENT(IN) :: LEGCOEF(6,0:MAXLEG,NDIST), SSALB(NDIST)
  CHARACTER(LEN=1), INTENT(IN) :: PARTYPE, DISTFLAG
  CHARACTER(LEN=*), INTENT(IN) :: MIETABFILE
  INTEGER :: I, J, M, N, JP, K, L, NL, NPHASEPOL

  OPEN (UNIT=3, FILE=MIETABFILE, STATUS='REPLACE')
  NPHASEPOL = 6
  WRITE (3,'(A)') '! Polarized Mie scattering table vs. effective radius (LWC=1 g/m^3)'
  IF (DELTAWAVE < 0.0) THEN
    WRITE (3,'(2(1X,F8.3),A)') WAVELEN1, WAVELEN2, '  wavelength range (micron)'
  ELSE
    WRITE (3,'(3(1X,F8.3),A)') WAVELEN1, WAVELEN2, DELTAWAVE, '  wavelength range and averaging step (micron)'
  ENDIF
  WRITE (3,'(1X,F5.3,2X,A1,A)') PARDENS, PARTYPE, '   particle density (g/cm^3) and type (Water, Ice, Aerosol)'
  WRITE (3,'(2(1X,E13.6),A)') RINDEX, '  particle index of refraction'
  IF (DISTFLAG == 'L') THEN
    WRITE (3,'(A)') 'lognormal log standard deviation'
  ELSE
    WRITE (3,'(A)') 'gamma size distribution shape parameter'
  ENDIF
  WRITE (3,'(1X,I3,2(1X,F8.3),A)') NRETAB,REFF(1),REFF(NRETAB),&
        '  number, starting, ending effective radius'
  WRITE (3,'(1X,I3,2(1X,F8.3),A)') NVETAB,VEFF(1),VEFF(NVETAB),&
        '  number, starting, ending effective variance'
    WRITE (3,'(1X,I3,2(1X,F8.3),A)') NVETAB,ALPHA(1),ALPHA(NVETAB),&
        '  number, starting, ending shape parameter alpha'
  WRITE (3,'(1X,I5,A)') MAXLEG, '  maximum rank allowed'
  DO N = 1, NVETAB
    DO I = 1, NRETAB
      M = (N-1)*NRETAB + I
      WRITE (3,'(1X,F8.4,1X,F8 .4,1X,E12.5,1X,F8.6,1X,I6,A)') &
        REFF(I), VEFF(N), EXTINCT(M), SSALB(M), NLEG(M),&
			'  Reff	 Veff  Ext  Alb  Nleg'
      DO K = 1, NPHASEPOL
        WRITE (3,'(I1,1X,201(1X,F10.5))') K, (LEGCOEF(K,L,M), L=0,MIN(NLEG(M),200))
      	DO J = 200, NLEG(M)-1, 200
          WRITE (3,'(2X,200(1X,F10.5))') (LEGCOEF(K,J+L,M),L=1,MIN(200,NLEG(M)-J))
      	ENDDO
      ENDDO
    ENDDO
  ENDDO
  CLOSE (3)
END SUBROUTINE WRITE_POLY_TABLE


SUBROUTINE READ_POLY_TABLE(MIETABFILE, NRETAB, NVETAB, NDIST, REFF, &
                           VEFF, EXTINCT, SSALB, NLEG, LEGCOEF, MAXLEG,&
			   TABLE_TYPE)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: MIETABFILE
  INTEGER, INTENT(IN) :: NRETAB, NVETAB, NDIST
  REAL, INTENT(OUT) :: REFF(NRETAB), VEFF(NVETAB)
  REAL, INTENT(OUT) :: EXTINCT(NDIST), SSALB(NDIST)
  INTEGER, INTENT(OUT) :: NLEG(NDIST)
  REAL, INTENT(OUT) :: LEGCOEF(6,0:MAXLEG,NDIST)
  INTEGER, INTENT(IN) :: MAXLEG
  CHARACTER(LEN=6), INTENT(OUT) :: TABLE_TYPE
  INTEGER :: I, J, L, M, N, D, K, NL, NPHASEPOL

  TABLE_TYPE = 'VECTOR'
  NPHASEPOL = 6
  OPEN (UNIT=1, FILE=MIETABFILE, STATUS='OLD')
    READ (1,*)
    READ (1,*)
    READ (1,*)
    READ (1,*)
    READ (1,*)
    READ (1,*)
    READ (1,*)
    READ (1,*)
    READ (1,*)
    DO N = 1, NVETAB
      DO I = 1, NRETAB
	M = (N-1)*NRETAB + I
      	READ (1,'(1X,F8.4,1X,F8.4,1X,E12.5,1X,F8.6,1X,I6,A)') &
          REFF(I), VEFF(N), EXTINCT(M), SSALB(M), NLEG(M)
      	DO K = 1, NPHASEPOL
      	  READ (1,'(I1,1X,201(1X,F10.5))'), D, (LEGCOEF(K,L,M), L=0,MIN(NLEG(M),200))
          DO J = 200, NLEG(M)-1, 200
            READ (1,'(2X,200(1X,F10.5))') (LEGCOEF(K,J+L,M),L=1,MIN(200,NLEG(M)-J))
	  ENDDO
        ENDDO
      ENDDO
    ENDDO
  CLOSE (1)
END SUBROUTINE READ_POLY_TABLE


SUBROUTINE TRANSFORM_LEG_TO_PHASE (MAXLEG, NPHASEPOL, &
				   PELEM, NLEG, LEGCOEF, &
				   NANGLE, ANGLE, PHASE)
 ! Transforms the phase matrix element (PELEM=1 to 6) from the Wigner
 ! d-function based coefficients of the scattering matrix (LEGCOEF) to
 ! a function of angle, PHASE(NANGLE) (at scattering angles in degrees
 ! in ANGLE(:)).  The order of the six elements in LEGCOEF is the four
 ! diagonal ones (alpha1, alpha2, alpha3, alpha4) followed by the IQ and
 ! UV coefficients (beta1, beta2).  The phase matrix elements indexed by
 ! PELEM are P11, P22, P33, P44, P12, P34.  If PELEM<0 then the ABS(PELEM)
 ! phase matrix element is normalized by the P11 element on output in PHASE.
 ! (Doicu et al., 2013, JQSRT, http://dx.doi.org/10.1016/j.jqsrt.2012.12.009).
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MAXLEG, NPHASEPOL, PELEM, NLEG, NANGLE
  REAL,    INTENT(IN) :: LEGCOEF(NPHASEPOL,0:MAXLEG), ANGLE(NANGLE)
  REAL,    INTENT(OUT) :: PHASE(NANGLE)
  INTEGER :: J, N
  DOUBLE PRECISION :: PI, X, XM, A1, A2, A3, A4, B1, B2, A2P3, A2M3
  DOUBLE PRECISION, ALLOCATABLE :: D00(:), D22P(:), D22M(:), D20(:)
  REAL, ALLOCATABLE :: NORM(:)

  IF (NPHASEPOL == 1 .AND. PELEM > 1) THEN
    PRINT '(A,I1,A)', 'Cannot do phase matrix element ',PELEM, &
       ' for unpolarized scattering table'
    STOP
  ENDIF

   ! Allocate arrays for Wigner functions D20, D22, D2-2 and D00
  ALLOCATE (D00(0:NLEG), D22P(0:NLEG), D22M(0:NLEG), D20(0:NLEG))
  ALLOCATE (NORM(NANGLE))

  PI = DACOS(-1.0D0)
  DO J = 1, NANGLE
    X = DCOS(ANGLE(J)*PI/180)
    XM = -X

     ! Calculate the desired scattering matrix element for this angle (X)
     ! from the Wigner d-function coefficients
    SELECT CASE (ABS(PELEM))
    CASE (1)
      CALL WIGNERFCT (X, NLEG, 0, 0, D00)
      A1 = 0.0D0
      DO N = 0, NLEG
	A1 = A1 + LEGCOEF(1,N) * D00(N)
      ENDDO
      PHASE(J) = A1

    CASE (2:3)
      CALL WIGNERFCT (X, NLEG, 2, 2, D22P)
      CALL WIGNERFCT (XM, NLEG, 2, 2, D22M) ! multiply by (-1)**N
      A2P3 = 0.0D0
      A2M3 = 0.0D0
      DO N = 2, NLEG
	A2P3 = A2P3 + (LEGCOEF(2,N) + LEGCOEF(3,N)) * D22P(N)
	A2M3 = A2M3 + (LEGCOEF(2,N) - LEGCOEF(3,N)) *(-1)**N * D22M(N)
      ENDDO
      A2 = 0.5D0 *(A2P3 + A2M3)
      A3 = 0.5D0 *(A2P3 - A2M3)
      IF (ABS(PELEM) == 2) THEN
	PHASE(J) = A2
      ELSE
	PHASE(J) = A3
      ENDIF

    CASE (4)
      CALL WIGNERFCT (X, NLEG, 0, 0, D00)
      A4 = 0.0D0
      DO N = 0, NLEG
	A4 = A4 + LEGCOEF(4,N) * D00(N)
      ENDDO
      PHASE(J) = A4

    CASE (5)
      CALL WIGNERFCT (X, NLEG, 2, 0, D20)
      B1 = 0.0D0
      DO N = 2, NLEG
	B1 = B1 - LEGCOEF(5,N) * D20(N)
      ENDDO
      PHASE(J) = B1

    CASE (6)
      CALL WIGNERFCT (X, NLEG, 2, 0, D20)
      B2 = 0.0D0
      DO N = 2, NLEG
	B2 = B2 - LEGCOEF(6,N) * D20(N)
      ENDDO
      PHASE(J) = B2

    CASE DEFAULT
      PRINT *, 'Illegal PELEM phase matrix element number'
      STOP
    END SELECT

     ! Calculate P11 function if normalization by it is desired
    IF (PELEM < 0) THEN
      CALL WIGNERFCT (X, NLEG, 0, 0, D00)
      A1 = 0.0D0
      DO N = 0, NLEG
	A1 = A1 + LEGCOEF(1,N) * D00(N)
      ENDDO
      NORM(J) = A1
    ENDIF
  ENDDO

   ! Do the P11 normalization if desired
  IF (PELEM < 0) THEN
    PHASE(:) = PHASE(:)/NORM(:)
  ENDIF
  DEALLOCATE (D00, D22P, D22M, D20, NORM)
END SUBROUTINE TRANSFORM_LEG_TO_PHASE



SUBROUTINE READ_MONO_TABLE(MIETABFILE, NRTAB, RADIUS, EXTINCT, &
		     SCATTER, NLEG, LEGCOEF, MAXLEG, TABLE_TYPE)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: MIETABFILE
  INTEGER, INTENT(IN) :: NRTAB
  REAL, INTENT(OUT) :: RADIUS(NRTAB), EXTINCT(NRTAB), SCATTER(NRTAB)
  INTEGER, INTENT(OUT) :: NLEG(NRTAB)
  REAL, INTENT(OUT) :: LEGCOEF(6,0:MAXLEG,NRTAB)
  CHARACTER(LEN=6), INTENT(OUT) :: TABLE_TYPE
  INTEGER, INTENT(IN) :: MAXLEG
  INTEGER :: I, J, L, M, K, NL, NPHASEPOL

  TABLE_TYPE = 'VECTOR'
  NPHASEPOL = 6
  OPEN (UNIT=1, FILE=MIETABFILE, STATUS='OLD')
    READ (1,*)
    READ (1,*)
    READ (1,*)
    READ (1,*)
    READ (1,*)
    READ (1,*)
    DO I = 1, NRTAB
      READ (1,'(1X,F8.4,1X,E12.5,1X,E12.5,1X,I6,A)') &
	  RADIUS(I), EXTINCT(I), SCATTER(I), NLEG(I)
      DO J = 1, NPHASEPOL
	READ (1,'(I1,1X,201(1X,F15.5))'), M, (LEGCOEF(J,L,I), L=0,MIN(NLEG(I),200))
	DO K = 200, NLEG(I)-1, 200
	  READ (1,'(2X,200(1X,F15.5))') (LEGCOEF(J,K+L,I),L=1,MIN(200,NLEG(I)-K))
	ENDDO
      ENDDO
    ENDDO
  CLOSE (1)
END SUBROUTINE READ_MONO_TABLE


SUBROUTINE WRITE_MONO_TABLE (MIETABFILE, WAVELEN1, WAVELEN2, DELTAWAVE, &
			    PARTYPE, PARDENS, RINDEX, NRTAB, &
			    RADII, EXTINCT, SCATTER, NLEG, MAXLEG, LEGCOEF)

 ! Writes the table of Mie monodisperse scattering properties as a function of
 ! radius.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NRTAB
  INTEGER, INTENT(IN) :: MAXLEG, NLEG(NRTAB)
  COMPLEX, INTENT(IN) :: RINDEX
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE
  REAL,    INTENT(IN) :: PARDENS
  REAL,    INTENT(IN) :: RADII(NRTAB), EXTINCT(NRTAB), SCATTER(NRTAB)
  REAL,    INTENT(IN) :: LEGCOEF(6,0:MAXLEG, NRTAB)
  CHARACTER(LEN=1), INTENT(IN) :: PARTYPE
  CHARACTER(LEN=*), INTENT(IN) :: MIETABFILE

  INTEGER :: I, J, JP, K, L, NL, NPHASEPOL
  OPEN (UNIT=3, FILE=MIETABFILE, STATUS='REPLACE')
  NPHASEPOL = 6
  WRITE (3,'(A)') '! Polarized Monodisperse Mie scattering table vs. radius'
  IF (DELTAWAVE < 0.0) THEN
    WRITE (3,'(2(1X,F8.3),A)') WAVELEN1, WAVELEN2, '  wavelength range (micron)'
  ELSE
    WRITE (3,'(3(1X,F8.3),A)') WAVELEN1, WAVELEN2, DELTAWAVE, '  wavelength range and averaging step (micron)'
  ENDIF
  WRITE (3,'(1X,F5.3,2X,A1,A)') PARDENS, PARTYPE, '   particle density (g/cm^3) and type (Water, Ice, Aerosol)'
  WRITE (3,'(2(1X,E13.6),A)') RINDEX, '  particle index of refraction'
  WRITE (3,'(1X,I5,2(1X,F8.3),A)') NRTAB, RADII(1), RADII(NRTAB), &
	'  number, starting, ending radius'
  WRITE (3,'(1X,I5,A)') MAXLEG, '  maximum rank allowed'
  DO I = 1, NRTAB
    WRITE (3,'(1X,F8.4,1X,E12.5,1X,E12.5,1X,I6,A)') &
        RADII(I), EXTINCT(I), SCATTER(I), NLEG(I), '  Radius  Ext  Scat  Nleg'
    DO J = 1, NPHASEPOL
    	WRITE (3,'(I1,1X,201(1X,F15.5))') J, (LEGCOEF(J,L,I), L=0,MIN(NLEG(I),200))
      DO K = 200, NLEG(I)-1, 200
        WRITE (3,'(2X,200(1X,F15.5))') (LEGCOEF(J,K+L,I),L=1,MIN(200,NLEG(I)-K))
      ENDDO
    ENDDO
  ENDDO
  CLOSE (3)
END SUBROUTINE WRITE_MONO_TABLE


SUBROUTINE GET_POLY_TABLE (ND, NDIST, NSIZE, MAXLEG, NLEG1, EXTINCT1,&
		      SCATTER1, LEGCOEF1, EXTINCT, SSALB, NLEG, LEGCOEF)
  IMPLICIT NONE
  REAL, INTENT(IN) :: ND(NSIZE, NDIST)
  INTEGER, INTENT(IN) :: NDIST, NSIZE
  INTEGER, INTENT(IN) :: MAXLEG
  INTEGER, INTENT(IN) :: NLEG1(NSIZE)
  REAL, INTENT(IN) :: EXTINCT1(NSIZE), SCATTER1(NSIZE)
  REAL, INTENT(IN) :: LEGCOEF1(6,0:MAXLEG, NSIZE)

  REAL, INTENT(OUT) :: EXTINCT(NDIST), SSALB(NDIST)
  REAL, INTENT(OUT) :: NLEG(NDIST), LEGCOEF(6,0:MAXLEG,NDIST)

  INTEGER :: I, J, L, NL
  REAL :: SCATTER

  ! Loop over the number of output tabulated effective radii
  DO I = 1, NDIST
    ! Sum the scattering properties over the discrete size distribution
    EXTINCT(I) = 0.0
    SCATTER = 0.0
    LEGCOEF(:,:,I) = 0.0
    NL = 1
    DO J = 1, NSIZE
      EXTINCT(I) = EXTINCT(I) + ND(J,I)*EXTINCT1(J)
      SCATTER = SCATTER + ND(J,I)*SCATTER1(J)
      NL = MAX(NL,NLEG1(J))
      LEGCOEF(:,0:NL,I) = LEGCOEF(:,0:NL,I) + ND(J,I)*LEGCOEF1(:,0:NL,J)
    ENDDO
    DO L = 0, NL
      LEGCOEF(:,L,I) = LEGCOEF(:,L,I)/SCATTER
      IF (LEGCOEF(1,L,I) .GT. 0.5E-5) NLEG(I) = L
    ENDDO
    IF (ABS(LEGCOEF(1,0,I)-1.0) > 0.0001) THEN
      PRINT *,'Phase function not normalized',I,LEGCOEF(1,0,I)
      STOP
    ENDIF
    IF (EXTINCT(I) > 0.0) THEN
      SSALB(I) = SCATTER/EXTINCT(I)
    ENDIF
    EXTINCT(I) = 0.001*EXTINCT(I)

  ENDDO  ! end of effective radius loop

END SUBROUTINE GET_POLY_TABLE
