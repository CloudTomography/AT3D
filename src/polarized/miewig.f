
      SUBROUTINE MIE_ONE (WAVELENGTH, MINDEX, RADIUS,
     .                   MAXRANK, EXTINCTION, SCATTER, NRANK, WIGCOEF)
C       Computes the Mie scattering properties for a single homogeneous
C     sphere of radius RADIUS.  The six phase matrix elements times the 
C     scattering coefficient is returned as Wigner d-function series 
C     coefficients.
      IMPLICIT NONE
      INTEGER     MAXRANK, NRANK
      REAL        WAVELENGTH, RADIUS
      REAL        EXTINCTION, SCATTER, WIGCOEF(6,0:*)
      COMPLEX     MINDEX
      INTEGER     MAXN
      PARAMETER   (MAXN=10000)
      INTEGER     NMIE, NQUAD, LASTNQUAD
      INTEGER     I, J, L
      DOUBLE PRECISION X, PI, QEXT, QSCAT, GEOMAREA
      DOUBLE PRECISION MU(MAXN), WTS(MAXN)
      DOUBLE PRECISION P1, P2, P3, P4
      DOUBLE PRECISION D00(0:MAXN), D20(0:MAXN)
      DOUBLE PRECISION D22P(0:MAXN), D22M(0:MAXN)
      DOUBLE PRECISION COEF(6,0:MAXN), F, A2, A3
      DOUBLE COMPLEX MSPHERE
      DOUBLE COMPLEX A(MAXN), B(MAXN)
      SAVE MU, WTS, LASTNQUAD
      DATA LASTNQUAD/-1/
      

      PI = ACOS(-1.0D0)      
      X = 2.0D0*PI*RADIUS/WAVELENGTH
      GEOMAREA = PI*RADIUS**2
      MSPHERE = MINDEX

C         Compute the An's and Bn's and the cross sections
      NMIE = 0
      CALL MIECALC (NMIE, X, MSPHERE, A, B)
      CALL MIECROSS (NMIE, X, A, B, QEXT, QSCAT)
      EXTINCTION = GEOMAREA*QEXT
      SCATTER = GEOMAREA*QSCAT

C         Get the Gauss-Legendre quadrature abscissas and weights
      NRANK = MIN(MAXRANK,2*NMIE)
      NQUAD  = (NRANK + 2*NMIE + 2)/2
      IF (NQUAD .LE. LASTNQUAD .AND. NQUAD.GE.NINT(0.8*LASTNQUAD)) THEN
        NQUAD = LASTNQUAD
      ELSE
        NQUAD = MIN(NINT(1.25*NQUAD),MAXRANK)
        IF (NQUAD .GT. MAXN) THEN
          PRINT *, 'MIE_ONE: MAXN exceeded by NQUAD'
          STOP
        ENDIF
        CALL GAUSQUAD (NQUAD, MU, WTS)
        LASTNQUAD = NQUAD
      ENDIF


C         Calculate the phase function at the quadrature angles and then
C         integrate over angle to get the Wigner d-function coefficients. 
      DO L = 0, NRANK
        DO J = 1, 6
          COEF(J,L) = 0.0D0
        ENDDO
      ENDDO
      DO I = 1, NQUAD
        CALL MIEANGLE (NMIE, A, B, MU(I), P1, P2, P3, P4)

C         Calculate the Wigner functions D20, D22, D2-2 and D00 for
C            l=0 to NRANK for MU(I)
        CALL WIGNERFCT (MU(I), NRANK, 0, 0, D00)
        CALL WIGNERFCT (MU(I), NRANK, 2, 2, D22P)
        CALL WIGNERFCT (-MU(I), NRANK, 2, 2, D22M)
        CALL WIGNERFCT (MU(I), NRANK, 2, 0, D20)
C         Calculate the series expansion coefficients
        DO L = 0, NRANK
          F = (-1.0D0)**L
          COEF(1,L) = COEF(1,L) + (L+0.5D0)*WTS(I)*P1*D00(L)
          COEF(2,L) = COEF(2,L) + (L+0.5D0)*WTS(I)*(P1+P3)*D22P(L)
          COEF(3,L) = COEF(3,L) + (L+0.5D0)*WTS(I)*(P1-P3)*F*D22M(L)
          COEF(4,L) = COEF(4,L) + (L+0.5D0)*WTS(I)*P3*D00(L)
          COEF(5,L) = COEF(5,L) - (L+0.5D0)*WTS(I)*P2*D20(L)
          COEF(6,L) = COEF(6,L) - (L+0.5D0)*WTS(I)*P4*D20(L)
        ENDDO
      ENDDO
      DO L = 0, NRANK
        A2 = 0.5D0*(COEF(2,L) + COEF(3,L))
        A3 = 0.5D0*(COEF(2,L) - COEF(3,L))
        COEF(2,L) = A2
        COEF(3,L) = A3
      ENDDO
      DO L = 0, NRANK
        DO J = 1, 6
          WIGCOEF(J,L) = SNGL( (WAVELENGTH**2/PI) *COEF(J,L) )
        ENDDO
      ENDDO
       
      RETURN
      END
 
 

 
 
      SUBROUTINE MIECALC (NTERMS, X, MN, A, B)
C        MIECALC calculates the complex Mie coefficients An and Bn
C      given the dimensionless size parameter X and the complex
C      index of refraction (Mre,Mim).  The number of terms calculated
C      is given by NTERMS unless NTERMS <= 0 or in which case the
C      appropriate number is calculated and returned in NTERMS.
      IMPLICIT NONE
      INTEGER   NTERMS
      DOUBLE PRECISION  X
      DOUBLE COMPLEX  MN, A(*), B(*)
      INTEGER     MAXTERMS,  NSTOP, N, NN
      PARAMETER   (MAXTERMS=10000)
      DOUBLE PRECISION  PSIN, PSIM, CHIN, CHIM, TMP
      DOUBLE PRECISION  DCOS, DSIN
      DOUBLE COMPLEX  M, Y, D(MAXTERMS+15), XIN, XIM, CTMP
      DOUBLE COMPLEX  DCMPLX
 
 
C           If NTERMS is not specified calculate it
      NSTOP = X + 4.0*X**0.3334 + 2
      IF (NTERMS .LE. 0)  NTERMS = NSTOP
      IF (NTERMS .GT. MAXTERMS) THEN
          WRITE (*,*)
     .       'Mie calculation requires more terms than available.'
          STOP
      ENDIF
 
C           Generate the Dn's by down recurrence  D = d(log(PSI(y)))/dy
      M = DCONJG(MN)
      Y = M*X
      NN = NTERMS + 15
      D(NN) = DCMPLX (0.0D0, 0.0D0)
      DO N = NN, 2, -1
          D(N-1) = N/Y - 1.0/ (D(N) + N/Y)
      ENDDO
 
C           Generate the PSIn's and XIn'S by upward recurrence
C             and calculate the An's and Bn's from them.
C             (PSIN = PSI(n), PSIM = PSI(n-1), same for CHI)
      PSIM = DCOS(X)
      PSIN = DSIN(X)
      CHIM = -DSIN(X)
      CHIN = DCOS(X)
      DO N = 1, NTERMS
          TMP = PSIN
          PSIN = (2*N-1)/X *PSIN - PSIM
          PSIM = TMP
          TMP = CHIN
          CHIN = (2*N-1)/X *CHIN - CHIM
          CHIM = TMP
          XIN = DCMPLX (PSIN, -CHIN)
          XIM = DCMPLX (PSIM, -CHIM)
          CTMP = D(N)/M + N/X
          A(N) = (CTMP*PSIN - PSIM) / (CTMP*XIN - XIM)
          CTMP = M*D(N) + N/X
          B(N) = (CTMP*PSIN - PSIM) / (CTMP*XIN - XIM)
      ENDDO
 
      RETURN
      END
 
 
 

      SUBROUTINE MIECROSS (NTERMS, X, A, B, QEXT, QSCAT)
C        MIECROSS calculates the extinction, scattering, and
C      backscatter efficiencies given the Mie coefficients An and Bn
C      and the size parameter X.
      IMPLICIT NONE
      INTEGER    NTERMS
      DOUBLE PRECISION  X, QEXT, QSCAT
      DOUBLE COMPLEX  A(*), B(*)
      INTEGER    N
      DOUBLE PRECISION  SUM1, SUM2, DREAL
 
      SUM1 = 0.0D0
      SUM2 = 0.0D0
      DO N = 1, NTERMS
          SUM1 = SUM1 + (2*N+1)*( DREAL(A(N)) + DREAL(B(N)) )
          SUM2 = SUM2 + (2*N+1)*( DREAL(A(N)*DCONJG(A(N)))
     .                          + DREAL(B(N)*DCONJG(B(N))) )
      ENDDO
      QEXT = 2.0D0/X**2 * SUM1
      QSCAT = 2.0D0/X**2 * SUM2

      RETURN
      END
 
 
 
 
      SUBROUTINE MIEANGLE (NTERMS, A, B, MU, P1, P2, P3, P4)
C        MIEANGLE calculates the scattering matrix elements (P1,P2,P3,P4)
C      for a particular value of MU (cos(theta)) from the Mie coefficients 
C      An's and Bn's.  The matrix elements are for the Stokes vector
C      (I,Q,U,V) and are calculated from the complex scattering 
C      amplitudes S1 and S2.
      IMPLICIT NONE
      INTEGER    NTERMS
      DOUBLE PRECISION  MU, P1, P2, P3, P4
      DOUBLE COMPLEX  A(*), B(*)
      INTEGER    N
      DOUBLE PRECISION TMP, PIN, PIM, TAUN, C
      DOUBLE COMPLEX  S1, S2
 
      S1 = DCMPLX(0.0,0.0)
      S2 = DCMPLX(0.0,0.0)
C               Sum up the series using the An's and Bn's
      PIN = 1.0
      PIM = 0.0
      DO N = 1, NTERMS
          TAUN = N*MU*PIN - (N+1)*PIM
C               Calculate the scattering functions at +mu and -mu
C                 using the PIn's and the TAUn's.
          C = (2*N+1) / DFLOAT(N*(N+1))
          S1 = S1 + C*( A(N)*PIN + B(N)*TAUN)
          S2 = S2 + C*( B(N)*PIN + A(N)*TAUN)
C               Calculate the angular function PIn by up recurrence
          TMP = PIN
          PIN = ( (2*N+1)*MU*PIN - (N+1)*PIM ) / N
          PIM = TMP
      ENDDO
C           Calculate the first Stokes parameter scattering matrix element
      P1 = 0.5*( CDABS(S2)**2 + CDABS(S1)**2 )
      P2 = 0.5*( CDABS(S2)**2 - CDABS(S1)**2 )
      P3 = DREAL( DCONJG(S1)*S2 )
      P4 = -DIMAG( DCONJG(S1)*S2 )
      RETURN
      END

 

      SUBROUTINE GAUSQUAD (N, XA, WT)
C        Generates the abscissas (X) and weights (W) for an N point
C      Gauss-Legendre quadrature.  
      IMPLICIT NONE
      INTEGER  N
      DOUBLE PRECISION  XA(*), WT(*)
      INTEGER  K, I, J, L
      REAL*8   X, XP, PL, PL1, PL2, DPL, TINY
      PARAMETER (TINY=3.0D-13)

      K = (N+1)/2
      DO 130 J = 1, K
        X = COS(3.141592654*(J-.25)/(N+.5))
        I = 0
100     CONTINUE
          PL1 = 1
          PL = X
          DO 120 L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
120       CONTINUE
          DPL = N*(X*PL-PL1)/(X*X-1)
          XP = X
          X = XP - PL/DPL
          I = I+1
        IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
        XA(J)     = -X
        XA(N-J+1) = X
        WT(J  )   = 2.0D0/((1.0D0-X*X)*DPL*DPL)
        WT(N-J+1) = WT(J)
130   CONTINUE

      RETURN
      END



