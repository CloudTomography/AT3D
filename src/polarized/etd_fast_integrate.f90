C ======================================================================
C     FAST QUADRATIC ETD CELL INTEGRATOR (FQETD)
C     Computes exact cell emission and transmission assuming linear
C     extinction and source gradients.
C
C     Arguments:
C       EXT0       : Extinction at upstream face (entry, s = 0)
C       EXT1       : Extinction at grid point (exit, s = SO)
C       SRCEXT0    : Source*Extinction array at upstream face (NSTOKES)
C       SRCEXT1    : Source*Extinction array at grid point (NSTOKES)
C       NSTOKES    : Number of Stokes components / channels
C       SO         : Path length through the cell segment
C       SRC        : Output integrated cell radiance array (NSTOKES)
C       TRANSCELL  : Output transmission across the cell
C ======================================================================
      SUBROUTINE FQETD(EXT0, EXT1, SRCEXT0, SRCEXT1, NSTOKES, SO,
     .                 SRC, TRANSCELL, ABSCELL)
      IMPLICIT NONE
      INTEGER NSTOKES
      DOUBLE PRECISION EXT0, EXT1, SO, TRANSCELL, ABSCELL
      DOUBLE PRECISION SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
      DOUBLE PRECISION SRC(NSTOKES)

Cf2py intent(in)  :: EXT0, EXT1, NSTOKES, SO
Cf2py intent(in)  :: SRCEXT0, SRCEXT1
Cf2py intent(out) :: SRC, TRANSCELL, ABSCELL
Cf2py depend(NSTOKES) :: SRCEXT0, SRCEXT1, SRC

      DOUBLE PRECISION EXT_GRAD, TAU, Z_SCALED, ETA_SCALED
      DOUBLE PRECISION ETA2, ETA3, ETA4, CQI
      DOUBLE PRECISION SQRT_GRAD, U0, U1, INT_0, INT_1
      DOUBLE PRECISION K_GRAD, SQRT_K, V0, V1, PI
      DOUBLE PRECISION PHI(10)
      INTEGER I
      
      DOUBLE PRECISION FAST_DAWSN, FAST_ERFCX
      EXTERNAL FAST_DAWSN, FAST_ERFCX

C     --- Edge Case Guard: Zero path length ---
      IF (SO .LE. 1.0D-12) THEN
          TRANSCELL = 1.0D0
          DO I = 1, NSTOKES
              SRC(I) = 0.0D0
          END DO
          RETURN
      END IF

      PI = 3.14159265358979323846D0
      EXT_GRAD = (EXT1 - EXT0) / (2.0D0 * SO)
      TAU = 0.5D0 * (EXT0 + EXT1) * SO
      IF (TAU .GE. 0.5D0) THEN
          TRANSCELL = DEXP(-TAU)
          ABSCELL = 1.0D0 - TRANSCELL
      ELSE
          ! Stable series expansion for small optical paths to prevent cancellation
          ABSCELL = TAU * (1.0D0 - 0.5D0 * TAU * 
     .              (1.0D0 - 0.333333333333333D0 * TAU * 
     .              (1.0D0 - 0.25D0 * TAU)))
          TRANSCELL = 1.0D0 - ABSCELL
      END IF

C     Scaled variables mapped to start of ray (EXT0)
      Z_SCALED = EXT0 * SO
      ETA_SCALED = EXT_GRAD * SO * SO

C     ==================================================================
C     Branch A: Small extinction gradient (|eta| < 1e-3)
C     ==================================================================
      IF (DABS(ETA_SCALED) .LT. 5.0D-2) THEN
          CALL COMPUTE_PHI_FAST(Z_SCALED, PHI)
          ETA2 = ETA_SCALED * ETA_SCALED
          ETA3 = ETA2 * ETA_SCALED
          ETA4 = ETA3 * ETA_SCALED

          INT_0 = SO * (PHI(1) - ETA_SCALED * PHI(3) 
     .          + 0.5D0 * ETA2 * PHI(5) - (ETA3 / 6.0D0) * PHI(7) 
     .          + (ETA4 / 24.0D0) * PHI(9))

          INT_1 = (SO * SO) * (PHI(2) - ETA_SCALED * PHI(4) 
     .          + 0.5D0 * ETA2 * PHI(6) - (ETA3 / 6.0D0) * PHI(8) 
     .          + (ETA4 / 24.0D0) * PHI(10))

C     ==================================================================
C     Branch B: Extinction decreases along ray (EXT1 < EXT0)
C     ==================================================================
      ELSE IF (EXT_GRAD .LT. 0.0D0) THEN
          K_GRAD = -EXT_GRAD
          SQRT_K = DSQRT(K_GRAD)
          U0 = -EXT0 / (2.0D0 * SQRT_K)
          U1 = SQRT_K * SO + U0

          INT_0 = (DSQRT(PI) / (2.0D0 * SQRT_K)) * 
     .         (TRANSCELL * ERFC_SCALED(U0) - ERFC_SCALED(U1))          
          INT_1 = (1.0D0 - TRANSCELL - EXT0 * INT_0) / 
     .            (2.0D0 * EXT_GRAD)

C     ==================================================================
C     Branch C: Extinction increases along ray (EXT1 > EXT0)
C     ==================================================================
      ELSE
          SQRT_GRAD = DSQRT(EXT_GRAD)
          V0 = EXT0 / (2.0D0 * SQRT_GRAD)
          V1 = SQRT_GRAD * SO + V0

          INT_0 = (1.0D0 / SQRT_GRAD) * (FAST_DAWSN(V1) 
     .          - TRANSCELL * FAST_DAWSN(V0))
          
          INT_1 = (1.0D0 - TRANSCELL - EXT0 * INT_0) / 
     .            (2.0D0 * EXT_GRAD)
      END IF

C     ==================================================================
C     Vectorized Source Integration over NSTOKES
C     ==================================================================
      DO I = 1, NSTOKES
          CQI = (SRCEXT1(I) - SRCEXT0(I)) / SO
          SRC(I) = SRCEXT0(I) * INT_0 + CQI * INT_1
      END DO

      RETURN
      END


C ======================================================================
C     LOOP-FREE DAWSONS INTEGRAL
C ======================================================================
      DOUBLE PRECISION FUNCTION FAST_DAWSN(X)
      IMPLICIT NONE
      DOUBLE PRECISION X, XABS, X2, NUM, DEN

      XABS = DABS(X)
      X2 = XABS * XABS

      IF (XABS .LT. 2.5D0) THEN
          NUM = 1.0D0 + X2 * (0.1051564023405313D0 + X2 * 
     .          (0.0345759132104932D0 + X2 * (0.0036034120395301D0 + 
     .          X2 * 0.0001095209304910D0)))
                  
          DEN = 1.0D0 + X2 * (0.7718230690071980D0 + X2 * 
     .          (0.2844781031940930D0 + X2 * (0.0531536031201940D0 + 
     .          X2 * 0.0042456010034020D0)))
                  
          FAST_DAWSN = XABS * (NUM / DEN)
      ELSE
          FAST_DAWSN = (1.0D0 + 0.5D0/X2 + 0.75D0/(X2*X2)) / 
     .                 (2.0D0 * XABS)
      END IF

      IF (X .LT. 0.0D0) FAST_DAWSN = -FAST_DAWSN

      RETURN
      END


C ======================================================================
C     LOOP-FREE SCALED ERROR FUNCTION
C ======================================================================
      DOUBLE PRECISION FUNCTION FAST_ERFCX(X)
      IMPLICIT NONE
      DOUBLE PRECISION X, T
      DOUBLE PRECISION P, A1, A2, A3, A4, A5
      PARAMETER (P  = 0.3275911D0)
      PARAMETER (A1 = 0.254829592D0)
      PARAMETER (A2 = -0.284496736D0)
      PARAMETER (A3 = 1.421413741D0)
      PARAMETER (A4 = -1.453152027D0)
      PARAMETER (A5 = 1.061405429D0)

      T = 1.0D0 / (1.0D0 + P * X)
      FAST_ERFCX = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))

      RETURN
      END


C ======================================================================
C     LOOP-FREE PHI FUNCTION EVALUATOR
C ======================================================================
      SUBROUTINE COMPUTE_PHI_FAST(Z_SCALED, PHI)
      IMPLICIT NONE
      DOUBLE PRECISION Z_SCALED, PHI(10)
      DOUBLE PRECISION EXP_Z, Z2, Z3, Z4, Z5, KREAL
      INTEGER K

      IF (Z_SCALED .LT. 0.5D0) THEN
          Z2 = Z_SCALED * Z_SCALED
          Z3 = Z2 * Z_SCALED
          Z4 = Z3 * Z_SCALED
          Z5 = Z4 * Z_SCALED
          
          DO K = 1, 10
              KREAL = DBLE(K)
              PHI(K) = (1.0D0 / KREAL)
     .               - (Z_SCALED / (KREAL + 1.0D0))
     .               + (Z2 / (2.0D0 * (KREAL + 2.0D0)))
     .               - (Z3 / (6.0D0 * (KREAL + 3.0D0)))
     .               + (Z4 / (24.0D0 * (KREAL + 4.0D0)))
     .               - (Z5 / (120.0D0 * (KREAL + 5.0D0)))
          END DO
      ELSE
          EXP_Z = DEXP(-Z_SCALED)
          PHI(1) = (1.0D0 - EXP_Z) / Z_SCALED
          DO K = 2, 10
              PHI(K) = (1.0D0 - DBLE(K - 1) * PHI(K - 1)) / Z_SCALED
          END DO
      END IF

      RETURN
      END
