C     Forward-adjoint inner product for SHDOM extinction derivatives.
C     Implements the five-term decomposition (Doicu & Efremenko 2019,
C     Eqs. 44-49) directly in the real generalized spherical harmonic
C     basis. The adjoint field is obtained from the pseudo-forward solve
C     by direction reversal ((-1)^L parity flip on SH coefficients and
C     Stokes-V sign negation for NSTOKES=4).
C
C     Angular integrals (Terms A, B) are computed via Parseval (SH
C     coefficient dot products). For the Stokes inner product, the simple
C     Parseval sum_J sum_s a(s,J)*b(s,J) is exact because the cross-
C     coupling between Stokes 2,3 through NSTLEG components 2,3,5,6
C     cancels when summing |Q|^2 + |U|^2.
C
C     Directional evaluations (Terms C, D, E) use the proper
C     reconstruction formula matching COMPUTE_ONE_SOURCE, coupling
C     NSTLEG components to NSTOKES Stokes parameters.
C
C     Spatial integrals use the trilinear mass matrix on leaf cells of
C     the adaptive grid (TREEPTR(2,IC)==0).
C
C     -JRLoveridge / derived from plan 2026

      SUBROUTINE COMPUTE_GRAD_INNER_PRODUCT (NSTOKES, NSTLEG,
     .           NX, NY, NZ, ML, MM, NLM, NPTS, NCELLS,
     .           GRIDPTR, TREEPTR, GRIDPOS, CELLFLAGS,
     .           RSHPTR, SHPTR, RADIANCE, SOURCE, DIRFLUX,
     .           RSHPTR_PF, RADIANCE_PF,
     .           DIRFLUX_PF, RADSOLAR_PF,
     .           SOLARMU, SOLARAZ, YLMSUN,
     .           DET_STOKES, CELLWISE_DET,
     .           MU_PF_CELL, PHI_PF_CELL,
     .           USE_LOGLINEAR_DIRFLUX,
     .           MAXPG,NUMDER, INTERPPTR,OPTINTERPWT, DEXTM,
     .           POINT_SENS_EXT, IERR, ERRMSG,
     .           POINT_SENS_RTGRID)
C     Computes the per-grid-point sensitivity POINT_SENS(1:NPTS) such
C     that for any extinction parameter x_p:
C       dI/dx_p = sum_g POINT_SENS(g) * d_sigma_ext(g)/dx_p
C
C     The sensitivity is the sum of five terms (A+B+C+D+E) accumulated
C     over leaf cells of the adaptive grid.
C
C     Inputs:
C       Forward solve arrays: RADIANCE, SOURCE, DIRFLUX, RSHPTR, SHPTR
C       Pseudo-forward solve arrays: RADIANCE_PF, SOURCE_PF, DIRFLUX_PF,
C                                    RSHPTR_PF, SHPTR_PF
C       Grid: GRIDPTR, TREEPTR, GRIDPOS, CELLFLAGS
C       Solar direction: SOLARMU, SOLARAZ, YLMSUN
C       Detector geometry: DET_STOKES (which Stokes component measured,
C         1=I,2=Q,3=U,4=V), CELLWISE_DET (T/F),
C         MU_PF_CELL/PHI_PF_CELL (pseudo-fwd beam direction per cell;
C         detector viewing direction = reversed from this)
C     Output:
C       POINT_SENS(NPTS) - per-grid-point extinction sensitivity
C
Cf2py threadsafe
      IMPLICIT NONE
C     --- Arguments ---
      INTEGER NSTOKES, NSTLEG, NX, NY, NZ
Cf2py intent(in) :: NSTOKES, NSTLEG, NX, NY, NZ
      INTEGER ML, MM, NLM, NPTS, NCELLS
Cf2py intent(in) :: ML, MM, NLM, NPTS, NCELLS
      INTEGER GRIDPTR(8,NCELLS), TREEPTR(2,NCELLS)
Cf2py intent(in) :: GRIDPTR, TREEPTR
      REAL    GRIDPOS(3,NPTS)
Cf2py intent(in) :: GRIDPOS
      INTEGER*2 CELLFLAGS(NCELLS)
Cf2py intent(in) :: CELLFLAGS
      INTEGER RSHPTR(NPTS+1), SHPTR(NPTS+1)
Cf2py intent(in) :: RSHPTR, SHPTR
      REAL    RADIANCE(NSTOKES,*), SOURCE(NSTOKES,*)
Cf2py intent(in) :: RADIANCE, SOURCE
      REAL    DIRFLUX(NPTS)
Cf2py intent(in) :: DIRFLUX
      INTEGER RSHPTR_PF(NPTS+1)
Cf2py intent(in) :: RSHPTR_PF
      REAL    RADIANCE_PF(NSTOKES,*), RADSOLAR_PF(NPTS)
Cf2py intent(in) :: RADIANCE_PF, RADSOLAR_PF
      REAL    DIRFLUX_PF(NPTS)
Cf2py intent(in) :: DIRFLUX_PF
      REAL    SOLARMU, SOLARAZ
Cf2py intent(in) :: SOLARMU, SOLARAZ
      REAL    YLMSUN(NSTLEG,NLM)
Cf2py intent(in) :: YLMSUN
      INTEGER DET_STOKES
Cf2py intent(in) :: DET_STOKES
      LOGICAL CELLWISE_DET
Cf2py intent(in) :: CELLWISE_DET
      REAL    MU_PF_CELL(NCELLS), PHI_PF_CELL(NCELLS)
Cf2py intent(in) :: MU_PF_CELL, PHI_PF_CELL
      LOGICAL USE_LOGLINEAR_DIRFLUX
Cf2py intent(in) :: USE_LOGLINEAR_DIRFLUX
      INTEGER INTERPPTR(8,NPTS)
Cf2py intent(in) :: INTERPPTR
      REAL    OPTINTERPWT(8,NPTS)
Cf2py intent(in) :: OPTINTERPWT
      INTEGER NUMDER
Cf2py intent(in) :: NUMDER
      INTEGER MAXPG
Cf2py intent(in) :: MAXPG
      REAL    DEXTM(MAXPG,NUMDER)
Cf2py intent(in) :: DEXTM
      REAL    POINT_SENS_EXT(MAXPG,NUMDER)
Cf2py intent(out) :: POINT_SENS_EXT
      REAL    POINT_SENS_RTGRID(NPTS)
Cf2py intent(out) :: POINT_SENS_RTGRID

      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG

C     --- Local variables ---
      INTEGER IC, I, K, J, L, S, IP(8), IDR, IB
      INTEGER NR_FWD, NR_PF, NS_FWD, NMIN
      REAL    DX, DY, DZ, VOL
      REAL    M3(8,8,8)
      REAL    MEXP(8,8), MEXP_PF(8,8)
      REAL    F0_CORNERS(8)
      REAL    P_AA(8,8), P_JA(8,8)
      REAL    IADJ_SOL(8), IFWD_DET(8), JFWD_DET(8)
      REAL    IADJ_SOL_FILT(8)
      REAL    A_IC(8), B_IC(8), C_IC(8), D_IC(8), E_IC(8)
      REAL    TA_IC(8), TB_IC(8)
      REAL    P_AA_FILT(8,8)
      REAL    ACC, ACC_FILT
      REAL    YLM_DET(NSTLEG,NLM)
      REAL    SIGN_ADJ(NLM)
      INTEGER LOFJ(NLM)
      REAL    PI_VAL
      PARAMETER (PI_VAL = 3.14159265358979)

      IERR = 0
      ERRMSG = ' '

C     --- Build LOFJ: sequential SH index -> degree L ---
C     Matches SHDOM ordering: L=0..ML, M=-min(L,MM)..+min(L,MM)
      J = 0
      DO L = 0, ML
        DO I = -MIN(L,MM), MIN(L,MM)
          J = J + 1
          LOFJ(J) = L
        ENDDO
      ENDDO

C     --- Build parity sign table: SIGN_ADJ(J) = (-1)^L(J) ---
      DO J = 1, NLM
        SIGN_ADJ(J) = 1.0 - 2.0*MOD(LOFJ(J), 2)
      ENDDO

C     --- Compute YLM at detector direction if global ---
C     Direction reversal: pseudo-forward propagates TOWARD medium from
C     detector, so detector viewing direction is reversed.
C     mu -> -mu, phi -> phi + pi
      IF (.NOT. CELLWISE_DET) THEN
        CALL YLMALL (.FALSE., -MU_PF_CELL(1),
     .               PHI_PF_CELL(1) + PI_VAL,
     .               ML, MM, NSTLEG, YLM_DET)
      ENDIF

C     --- Initialize output ---
      DO I = 1, MAXPG
        DO J = 1, NUMDER
          POINT_SENS_EXT(I,J) = 0.0
        ENDDO
      ENDDO

      DO I = 1, NPTS
        POINT_SENS_RTGRID(I) = 0.0
      ENDDO

C     ================================================================
C     Main loop over leaf cells
C     ================================================================
      DO IC = 1, NCELLS
C       Skip non-leaf cells (parents that have been split)
        IF (TREEPTR(2,IC) .NE. 0) CYCLE

C       Gather 8 corner grid point indices
        DO I = 1, 8
          IP(I) = GRIDPTR(I, IC)
        ENDDO

C       Compute cell dimensions from corner positions
C       GRIDPTR ordering: IOCT = 1 + BITX + 2*BITY + 4*BITZ
C       so point 1 is (xmin,ymin,zmin), point 8 is (xmax,ymax,zmax)
        DX = GRIDPOS(1,IP(2)) - GRIDPOS(1,IP(1))
        DY = GRIDPOS(2,IP(3)) - GRIDPOS(2,IP(1))
        DZ = GRIDPOS(3,IP(5)) - GRIDPOS(3,IP(1))
        VOL = DX * DY * DZ

C       Skip degenerate cells (zero-volume, e.g. open BC placeholders)
        IF (VOL .LE. 0.0) CYCLE

C       Build the 8x8x8 trilinear mass matrix M3(i,j,k) = integral(Li*Lj*Lk*dV)
        CALL BUILD_TRILIN_MASS_MATRIX_3(DX, DY, DZ, M3)

C       Compute YLM at detector direction (per-cell if cellwise)
        IF (CELLWISE_DET) THEN
          CALL YLMALL (.FALSE., -MU_PF_CELL(IC),
     .                 PHI_PF_CELL(IC) + PI_VAL,
     .                 ML, MM, NSTLEG, YLM_DET)
        ENDIF

C       ============================================================
C       Term A: diffuse-diffuse Parseval with cross-node products
C       P_AA(j,k) = integral_S2 I_d(j) * I_adj(k) dOmega
C         = sum_J sum_s RAD(s,J@j) * SIGN_ADJ(J) * RAD_PF(s,J@k)
C       For s=1,2,3: RSTOKES=+1; for s=4: RSTOKES=-1 (Stokes V flip)
C       A_i = -sum_j sum_k M3(i,j,k) * P_AA(j,k)
C       Also computes P_AA_FILT for albedo gradient (Term TA):
C       P_AA_FILT(j,k) = sum_J CHI_NORM(J,j)*RAD(J@j)*SIGN_ADJ(J)*RAD_PF(J@k)
C       ============================================================
        DO I = 1, 8
          DO K = 1, 8
            NR_FWD = RSHPTR(IP(I)+1) - RSHPTR(IP(I))
            NR_PF  = RSHPTR_PF(IP(K)+1) - RSHPTR_PF(IP(K))
            NMIN   = MIN(NR_FWD, NR_PF)
            ACC = 0.0
            IF (NSTOKES .LT. 4) THEN
              DO J = 1, NMIN
                DO S = 1, NSTOKES
                  ACC = ACC + RADIANCE(S, RSHPTR(IP(I))+J)
     .                      * SIGN_ADJ(J)
     .                      * RADIANCE_PF(S, RSHPTR_PF(IP(K))+J)
                ENDDO
              ENDDO
            ELSE
              DO J = 1, NMIN
                DO S = 1, 3
                  ACC = ACC + RADIANCE(S, RSHPTR(IP(I))+J)
     .                      * SIGN_ADJ(J)
     .                      * RADIANCE_PF(S, RSHPTR_PF(IP(K))+J)
                ENDDO
                ACC = ACC - RADIANCE(4, RSHPTR(IP(I))+J)
     .                    * SIGN_ADJ(J)
     .                    * RADIANCE_PF(4, RSHPTR_PF(IP(K))+J)
              ENDDO
            ENDIF
            P_AA(I,K) = ACC
          ENDDO
        ENDDO

C       A_ic(i) = -sum_j sum_k M3(i,j,k) * P_AA(j,k)
        DO I = 1, 8
          ACC = 0.0
          DO J = 1, 8
            DO K = 1, 8
              ACC = ACC + M3(I,J,K) * P_AA(J,K)
            ENDDO
          ENDDO
          A_IC(I) = -ACC
        ENDDO

C       ============================================================
C       Term B: cross-Parseval with cross-node products
C       P_JA(j,k) = integral_S2 J_fwd(j) * I_adj(k) dOmega
C       B_i = +sum_j sum_k M3(i,j,k) * P_JA(j,k)
C       ============================================================
C       TODO: Need to subtract filtered direct beam and evaluate
C             Exact derivatives using exact phase function and
C             solar transmission derivatives.
C             Put this option under delta-M.
        DO I = 1, 8
          DO K = 1, 8
            NS_FWD = SHPTR(IP(I)+1) - SHPTR(IP(I))
            NR_PF  = RSHPTR_PF(IP(K)+1) - RSHPTR_PF(IP(K))
            NMIN   = MIN(NS_FWD, NR_PF)
            ACC = 0.0

            IF (NSTOKES .LT. 4) THEN
              DO J = 1, NMIN
                DO S = 1, NSTOKES
                  ACC = ACC + SOURCE(S, SHPTR(IP(I))+J)
     .                      * SIGN_ADJ(J)
     .                      * RADIANCE_PF(S, RSHPTR_PF(IP(K))+J)
                ENDDO
              ENDDO
            ELSE
              DO J = 1, NMIN
                DO S = 1, 3
                  ACC = ACC + SOURCE(S, SHPTR(IP(I))+J)
     .                      * SIGN_ADJ(J)
     .                      * RADIANCE_PF(S, RSHPTR_PF(IP(K))+J)
                ENDDO
                ACC = ACC - SOURCE(4, SHPTR(IP(I))+J)
     .                    * SIGN_ADJ(J)
     .                    * RADIANCE_PF(4, RSHPTR_PF(IP(K))+J)
              ENDDO
            ENDIF
            P_JA(I,K) = ACC
          ENDDO
        ENDDO

C       B_ic(i) = +sum_j sum_k M3(i,j,k) * P_JA(j,k)
        DO I = 1, 8
          ACC = 0.0
          DO J = 1, 8
            DO K = 1, 8
              ACC = ACC + M3(I,J,K) * P_JA(J,K)
            ENDDO
          ENDDO
          B_IC(I) = ACC
        ENDDO

C       ============================================================
C       Terms C, D, E: directional SH evaluations
C       These use the proper reconstruction formula from
C       COMPUTE_ONE_SOURCE, coupling NSTLEG to Stokes.
C       Term C: adjoint diffuse, Stokes-I only (solar unpolarized)
C       Term D: forward diffuse, Stokes DET_STOKES component
C       Term E: forward source, Stokes DET_STOKES component
C       ============================================================
        DO K = 1, 8
C         ---- Term C: adjoint Stokes-I at solar direction ----
C         hat_I_adj_1(Omega_0) = sum_J SIGN_ADJ(J)*RAD_PF(1,J)
C                                      *YLMSUN(1,J)
C         Only Stokes-I needed; YLMSUN(1,J) unaffected by TRANSPOSE.
C         Also compute filtered version for albedo gradient (Term TB).

          IADJ_SOL(K) = -RADSOLAR_PF(IP(K))

C         ---- Term D: forward diffuse at detector direction ----
C         Evaluate Stokes component DET_STOKES of forward radiance
C         at the detector direction (YLMALL with TRANSPOSE=.FALSE.)
          NR_FWD = RSHPTR(IP(K)+1) - RSHPTR(IP(K))
          CALL EVAL_FIELD_AT_DIR(RADIANCE(1, RSHPTR(IP(K))+1),
     .         NR_FWD, NSTOKES, NSTLEG, NLM, YLM_DET,
     .         DET_STOKES, IFWD_DET(K))

C         ---- Term E: forward source at detector direction ----
          NS_FWD = SHPTR(IP(K)+1) - SHPTR(IP(K))
          CALL EVAL_FIELD_AT_DIR(SOURCE(1, SHPTR(IP(K))+1),
     .         NS_FWD, NSTOKES, NSTLEG, NLM, YLM_DET,
     .         DET_STOKES, JFWD_DET(K))
        ENDDO


C       Build exponentially-weighted mass matrices for direct beams.
C       Three modes: USE_LOGLINEAR_DIRFLUX selects beam-traced analytical
C       (exact for any tau) vs Gauss-based log-linear vs trilinear M3.
C       C_i = sum_j MEXP(i,j) * IADJ_SOL(j)
C       D_i = sum_j MEXP_PF(i,j) * IFWD_DET(j)
C       E_i = -sum_j MEXP_PF(i,j) * JFWD_DET(j)
        IF (USE_LOGLINEAR_DIRFLUX) THEN
C         Beam-traced analytical exponential mass matrix
          DO K = 1, 8
            F0_CORNERS(K) = DIRFLUX(IP(K)) / ABS(SOLARMU)
          ENDDO
          CALL BUILD_EXP_MASS_MATRIX(DX, DY, DZ, F0_CORNERS, MEXP)
C     .         SOLARMU, SOLARAZ, MEXP)
          DO K = 1, 8
            F0_CORNERS(K) = DIRFLUX_PF(IP(K)) / ABS(MU_PF_CELL(IC))
          ENDDO
          CALL BUILD_EXP_MASS_MATRIX(DX, DY, DZ, F0_CORNERS, MEXP_PF)
C     .         MU_PF_CELL(IC), PHI_PF_CELL(IC), MEXP_PF)

          DO I = 1, 8
            C_IC(I) = 0.0
            D_IC(I) = 0.0
            E_IC(I) = 0.0
            DO J = 1, 8
              C_IC(I) = C_IC(I) + MEXP(I,J) * IADJ_SOL(J)
              D_IC(I) = D_IC(I) + MEXP_PF(I,J) * IFWD_DET(J)
              E_IC(I) = E_IC(I) - MEXP_PF(I,J) * JFWD_DET(J)
            ENDDO
          ENDDO
        ELSE
C         Trilinear (M3-based) direct beam: sum_j sum_k M3(i,j,k)*f(k)*field(j)
          DO I = 1, 8
            C_IC(I) = 0.0
            D_IC(I) = 0.0
            E_IC(I) = 0.0
            DO J = 1, 8
              DO K = 1, 8
                ACC = M3(I,J,K) * DIRFLUX(IP(K)) / ABS(SOLARMU)
                C_IC(I) = C_IC(I) + ACC * IADJ_SOL(J)
                ACC = M3(I,J,K) * DIRFLUX_PF(IP(K))
     .                           / ABS(MU_PF_CELL(IC))
                D_IC(I) = D_IC(I) + ACC * IFWD_DET(J)
                E_IC(I) = E_IC(I) - ACC * JFWD_DET(J)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        DO I = 1, 8
          POINT_SENS_RTGRID(IP(I)) = POINT_SENS_RTGRID(IP(I)) 
     .      + A_IC(I) + B_IC(I) + C_IC(I)
        ENDDO

C       Scatter per-corner contributions to global point sensitivity
        DO IDR = 1, NUMDER
          DO IB = 1, 8
            DO I = 1, 8
              POINT_SENS_EXT(INTERPPTR(IB,IP(I)),IDR) = 
     .          POINT_SENS_EXT(INTERPPTR(IB,IP(I)),IDR)
     .        + OPTINTERPWT(IB,IP(I))*DEXTM(IB,IDR)*(
     .        + A_IC(I) + B_IC(I)
     .        + C_IC(I))
C .        + D_IC(I) + E_IC(I) )
            ENDDO
          ENDDO
        ENDDO

      ENDDO
C     End of leaf cell loop

      RETURN
      END


C     ================================================================
      SUBROUTINE BUILD_TRILIN_MASS_MATRIX_3 (DX, DY, DZ, M3)
C     Builds the 8x8x8 rank-3 trilinear mass matrix for a rectangular
C     cell.  M3(i,j,k) = integral over cell of L_i(r)*L_j(r)*L_k(r) dV
C     where L_i are trilinear basis functions.
C
C     This is the tensor product of three 1D cubic integrals:
C       integral_0^h phi_a*phi_b*phi_c dx
C     where phi_0(x)=1-x/h, phi_1(x)=x/h.  The values are:
C       (0,0,0) -> h/4,  (1,1,1) -> h/4
C       all other combinations -> h/12
C
C     Node ordering follows SHDOM convention:
C       IOCT = 1 + BITX + 2*BITY + 4*BITZ
C
      IMPLICIT NONE
      REAL DX, DY, DZ, M3(8,8,8)
Cf2py intent(in) :: DX, DY, DZ
Cf2py intent(out) :: M3
      REAL W3X(2,2,2), W3Y(2,2,2), W3Z(2,2,2)
      INTEGER I, J, K, IX, IY, IZ, JX, JY, JZ, KX, KY, KZ
      INTEGER A, B, C

C     Fill 1D cubic mass weights (unit interval values * 12)
C     w(a,b,c): all-same = 3, otherwise = 1
      DO A = 1, 2
        DO B = 1, 2
          DO C = 1, 2
            IF (A.EQ.B .AND. B.EQ.C) THEN
              W3X(A,B,C) = 3.0
              W3Y(A,B,C) = 3.0
              W3Z(A,B,C) = 3.0
            ELSE
              W3X(A,B,C) = 1.0
              W3Y(A,B,C) = 1.0
              W3Z(A,B,C) = 1.0
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C     M3(i,j,k) = (DX*DY*DZ)/(12^3) * Wx(ix,jx,kx)*Wy(iy,jy,ky)*Wz(iz,jz,kz)
      DO I = 1, 8
        IX = MOD(I-1, 2) + 1
        IY = MOD((I-1)/2, 2) + 1
        IZ = (I-1)/4 + 1
        DO J = 1, 8
          JX = MOD(J-1, 2) + 1
          JY = MOD((J-1)/2, 2) + 1
          JZ = (J-1)/4 + 1
          DO K = 1, 8
            KX = MOD(K-1, 2) + 1
            KY = MOD((K-1)/2, 2) + 1
            KZ = (K-1)/4 + 1
            M3(I,J,K) = (DX * DY * DZ) / 1728.0
     .                * W3X(IX,JX,KX) * W3Y(IY,JY,KY)
     .                * W3Z(IZ,JZ,KZ)
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END


C     ================================================================
      SUBROUTINE EVAL_FIELD_AT_DIR (COEFFS, NSH, NSTOKES, NSTLEG,
     .                              NLM, YLMDIR, ISTOKES, RESULT)
C     Evaluates one Stokes component of an SH-expanded field at a
C     direction defined by YLMDIR. Uses the coupling formula from
C     COMPUTE_ONE_SOURCE:
C       Stokes 1: sum_J COEFFS(1,J) * YLMDIR(1,J)
C       Stokes 2: sum_J COEFFS(2,J)*YLMDIR(2,J)+COEFFS(3,J)*YLMDIR(5,J)
C       Stokes 3: sum_J COEFFS(2,J)*YLMDIR(6,J)+COEFFS(3,J)*YLMDIR(3,J)
C       Stokes 4: sum_J COEFFS(4,J) * YLMDIR(4,J)
C
C     COEFFS(NSTOKES, NSH) - SH coefficients at a grid point
C     YLMDIR(NSTLEG, NLM)  - real generalized SH at target direction
C     ISTOKES              - which Stokes component to evaluate (1-4)
C     RESULT               - output scalar value
C
      IMPLICIT NONE
      INTEGER NSH, NSTOKES, NSTLEG, NLM, ISTOKES
Cf2py intent(in) :: COEFFS, YLMDIR, ISTOKES
      REAL    COEFFS(NSTOKES, NSH)
      REAL    YLMDIR(NSTLEG, NLM)
      REAL    RESULT
Cf2py intent(out) :: RESULT
      INTEGER J
      REAL    ACC

      ACC = 0.0

      IF (ISTOKES .EQ. 1) THEN
C       Stokes I: only NSTLEG component 1
        DO J = 1, NSH
          ACC = ACC + COEFFS(1,J) * YLMDIR(1,J)
        ENDDO

      ELSEIF (ISTOKES .EQ. 2) THEN
C       Stokes Q: NSTLEG components 2 and 5
        IF (NSTOKES .GE. 3) THEN
          DO J = 1, NSH
            ACC = ACC + COEFFS(2,J) * YLMDIR(2,J)
     .                + COEFFS(3,J) * YLMDIR(5,J)
          ENDDO
        ENDIF

      ELSEIF (ISTOKES .EQ. 3) THEN
C       Stokes U: NSTLEG components 6 and 3
        IF (NSTOKES .GE. 3) THEN
          DO J = 1, NSH
            ACC = ACC + COEFFS(2,J) * YLMDIR(6,J)
     .                + COEFFS(3,J) * YLMDIR(3,J)
          ENDDO
        ENDIF

      ELSEIF (ISTOKES .EQ. 4) THEN
C       Stokes V: only NSTLEG component 4
        IF (NSTOKES .GE. 4) THEN
          DO J = 1, NSH
            ACC = ACC + COEFFS(4,J) * YLMDIR(4,J)
          ENDDO
        ENDIF
      ENDIF

      RESULT = ACC
      RETURN
      END


C     ================================================================
      SUBROUTINE BUILD_EXP_MASS_MATRIX (DX, DY, DZ, F0, MEXP)
C     Builds the 8x8 exponentially-weighted mass matrix:
C       MEXP(i,j) = integral over cell of L_i(r)*L_j(r)*F0_loglin(r) dV
C     where F0_loglin(r) = exp(sum_k L_k(r)*ln(F0(k))) is the log-linear
C     interpolation of the direct flux. This avoids the Jensen inequality
C     bias from trilinear interpolation of exponentially-varying fields.
C
C     Evaluated via 2x2x2 Gauss-Legendre quadrature on the unit cube.
C     Gauss points: g = (1 +/- 1/sqrt(3))/2
C     Basis function values at Gauss points are precomputed constants.
C
C     When any F0(k) <= 0, that corner's log is clamped to -50
C     (exp(-50) ~ 2e-22, effectively zero contribution).
C
      IMPLICIT NONE
      REAL DX, DY, DZ, F0(8), MEXP(8,8)
Cf2py intent(in) :: DX, DY, DZ, F0
Cf2py intent(out) :: MEXP

C     Gauss-Legendre points on [0,1]: (1 +/- 1/sqrt(3))/2
      REAL    GP(2)
      REAL    VOL, LOGF0(8), CVAL, FQ
      REAL    PHI_X(2,2), PHI_Y(2,2), PHI_Z(2,2)
      REAL    BF(8)
      INTEGER I, J, K, QX, QY, QZ, IX, IY, IZ
      REAL    CLAMP
      PARAMETER (CLAMP = -50.0)

C     Gauss points on [0,1]
      GP(1) = 0.5 - 0.5/SQRT(3.0)
      GP(2) = 0.5 + 0.5/SQRT(3.0)

      VOL = DX * DY * DZ

C     Compute log of corner fluxes with clamping
      DO K = 1, 8
        IF (F0(K) .GT. 0.0) THEN
          LOGF0(K) = LOG(F0(K))
        ELSE
          LOGF0(K) = CLAMP
        ENDIF
      ENDDO

C     Precompute 1D basis function values at the 2 Gauss points
C     phi_0(g) = 1-g,  phi_1(g) = g
      DO I = 1, 2
        PHI_X(1,I) = 1.0 - GP(I)
        PHI_X(2,I) = GP(I)
        PHI_Y(1,I) = 1.0 - GP(I)
        PHI_Y(2,I) = GP(I)
        PHI_Z(1,I) = 1.0 - GP(I)
        PHI_Z(2,I) = GP(I)
      ENDDO

C     Initialize
      DO I = 1, 8
        DO J = 1, 8
          MEXP(I,J) = 0.0
        ENDDO
      ENDDO

C     2x2x2 Gauss quadrature
      DO QZ = 1, 2
        DO QY = 1, 2
          DO QX = 1, 2
C           Compute 8 basis function values at this Gauss point
C           Node K: ix=MOD(K-1,2)+1, iy=MOD((K-1)/2,2)+1, iz=(K-1)/4+1
            DO K = 1, 8
              IX = MOD(K-1, 2) + 1
              IY = MOD((K-1)/2, 2) + 1
              IZ = (K-1)/4 + 1
              BF(K) = PHI_X(IX,QX) * PHI_Y(IY,QY) * PHI_Z(IZ,QZ)
            ENDDO

C           Log-linear interpolation: F0(r_q) = exp(sum_k BF(k)*logF0(k))
            CVAL = 0.0
            DO K = 1, 8
              CVAL = CVAL + BF(K) * LOGF0(K)
            ENDDO
            FQ = EXP(CVAL)

C           Accumulate: MEXP(i,j) += w_q * BF(i) * BF(j) * FQ * VOL
C           Weight w_q = 1/8 for each of 8 Gauss points on unit cube
            DO I = 1, 8
              DO J = 1, 8
                MEXP(I,J) = MEXP(I,J)
     .                    + BF(I) * BF(J) * FQ
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C     Multiply by VOL * w_q = VOL/8
      DO I = 1, 8
        DO J = 1, 8
          MEXP(I,J) = MEXP(I,J) * VOL / 8.0
        ENDDO
      ENDDO

      RETURN
      END


C     ================================================================
      SUBROUTINE BUILD_EXP_MASS_MATRIX_BEAM (DX, DY, DZ, F0,
     .                                       BEAMMU, BEAMAZ, MEXP)
C     Builds the 8x8 exponentially-weighted mass matrix using analytical
C     integration along the beam direction with 2D Gauss quadrature on
C     the beam entry face. Exact for any cell optical depth.
C
C     The direct flux is modeled as F0(r) decaying exponentially along
C     the beam propagation direction. The effective tau along each ray
C     is computed from the corner F0 values (bilinear interpolation of
C     ln(F0) on entry and exit faces).
C
C     Method:
C       1. Determine entry face from beam direction
C       2. For each 2x2 Gauss point on entry face, trace ray across cell
C       3. Along each ray, L_i(t)*L_j(t) is degree-6 polynomial in t
C       4. Integrate L_i*L_j*exp(-tau*t) analytically using In(tau)
C       5. Sum over entry-face quadrature points
C
C     BEAMMU: cosine of beam zenith angle (negative = downward)
C     BEAMAZ: beam azimuth angle in radians
C
      IMPLICIT NONE
      REAL DX, DY, DZ, F0(8), MEXP(8,8)
      REAL BEAMMU, BEAMAZ
Cf2py intent(in) :: DX, DY, DZ, F0, BEAMMU, BEAMAZ
Cf2py intent(out) :: MEXP

C     Local variables
      REAL    SX, SY, SZ, SINTH
      REAL    GP(2), WGP
      REAL    ENTRY_U, ENTRY_V
      REAL    X0, Y0, Z0, SMAX
      REAL    TX, TY, TZ, TMIN
      REAL    LOGF0(8), LF_ENTRY, LF_EXIT
      REAL    TAU, ETAU, F0_ENTRY
      REAL    EIN(0:6)
      REAL    AX(8), BX(8), AY(8), BY(8), AZ(8), BZ(8)
      REAL    LI_A, LI_B, LJ_A, LJ_B
      REAL    POLY_I(0:3), POLY_J(0:3), POLY_IJ(0:6)
      REAL    CONTRIB
      INTEGER I, J, K, Q1, Q2, IX, IY, IZ
      INTEGER ENTRY_FACE
      REAL    CLAMP, PI_VAL
      PARAMETER (CLAMP = -50.0, PI_VAL = 3.14159265358979)
      REAL    FACE_AREA
      REAL    UQ, VQ
      REAL    XE, YE, ZE

C     Gauss points on [0,1]
      GP(1) = 0.5 - 0.5/SQRT(3.0)
      GP(2) = 0.5 + 0.5/SQRT(3.0)
      WGP = 0.25

C     Beam direction vector (unit vector, pointing in propagation dir)
      SINTH = SQRT(1.0 - BEAMMU*BEAMMU)
      SX = SINTH * COS(BEAMAZ)
      SY = SINTH * SIN(BEAMAZ)
      SZ = BEAMMU

C     Determine entry face: face with largest |s_d/d_d|
C     Entry face codes: 1=x0, 2=x1, 3=y0, 4=y1, 5=z0, 6=z1
C     Beam enters the face it hits FIRST traveling in +s direction
      TX = 1.0E30
      TY = 1.0E30
      TZ = 1.0E30
      IF (ABS(SX) .GT. 1.0E-10) THEN
        IF (SX .GT. 0.0) THEN
          TX = 0.0
        ELSE
          TX = DX / ABS(SX)
        ENDIF
        TX = DX / ABS(SX)
      ENDIF
      IF (ABS(SY) .GT. 1.0E-10) THEN
        IF (SY .GT. 0.0) THEN
          TY = 0.0
        ELSE
          TY = DY / ABS(SY)
        ENDIF
        TY = DY / ABS(SY)
      ENDIF
      IF (ABS(SZ) .GT. 1.0E-10) THEN
        IF (SZ .GT. 0.0) THEN
          TZ = 0.0
        ELSE
          TZ = DZ / ABS(SZ)
        ENDIF
        TZ = DZ / ABS(SZ)
      ENDIF

C     Entry face is the one the beam crosses in minimum time
C     For a downward beam (SZ<0), it enters through z=DZ (top face)
C     Determine by comparing traversal times for each axis
      IF (ABS(SX)/DX .GE. ABS(SY)/DY .AND.
     .    ABS(SX)/DX .GE. ABS(SZ)/DZ) THEN
        IF (SX .GT. 0.0) THEN
          ENTRY_FACE = 1
        ELSE
          ENTRY_FACE = 2
        ENDIF
      ELSEIF (ABS(SY)/DY .GE. ABS(SX)/DX .AND.
     .        ABS(SY)/DY .GE. ABS(SZ)/DZ) THEN
        IF (SY .GT. 0.0) THEN
          ENTRY_FACE = 3
        ELSE
          ENTRY_FACE = 4
        ENDIF
      ELSE
        IF (SZ .GT. 0.0) THEN
          ENTRY_FACE = 5
        ELSE
          ENTRY_FACE = 6
        ENDIF
      ENDIF

C     Compute log of corner fluxes
      DO K = 1, 8
        IF (F0(K) .GT. 0.0) THEN
          LOGF0(K) = LOG(F0(K))
        ELSE
          LOGF0(K) = CLAMP
        ENDIF
      ENDDO

C     Path length for a ray fully traversing the cell along beam dir
C     SMAX = min time to exit any face (from entry face)
      TMIN = 1.0E30
      IF (ABS(SX) .GT. 1.0E-10) TMIN = MIN(TMIN, DX/ABS(SX))
      IF (ABS(SY) .GT. 1.0E-10) TMIN = MIN(TMIN, DY/ABS(SY))
      IF (ABS(SZ) .GT. 1.0E-10) TMIN = MIN(TMIN, DZ/ABS(SZ))
      SMAX = TMIN

C     Initialize output
      DO I = 1, 8
        DO J = 1, 8
          MEXP(I,J) = 0.0
        ENDDO
      ENDDO

C     2x2 Gauss quadrature on entry face
      DO Q1 = 1, 2
        DO Q2 = 1, 2
          UQ = GP(Q1)
          VQ = GP(Q2)

C         Compute entry point (x0,y0,z0) and face area based on entry face
          IF (ENTRY_FACE .EQ. 1) THEN
C           x=0 face, u=y, v=z
            X0 = 0.0
            Y0 = UQ * DY
            Z0 = VQ * DZ
            FACE_AREA = DY * DZ
          ELSEIF (ENTRY_FACE .EQ. 2) THEN
C           x=DX face, u=y, v=z
            X0 = DX
            Y0 = UQ * DY
            Z0 = VQ * DZ
            FACE_AREA = DY * DZ
          ELSEIF (ENTRY_FACE .EQ. 3) THEN
C           y=0 face, u=x, v=z
            X0 = UQ * DX
            Y0 = 0.0
            Z0 = VQ * DZ
            FACE_AREA = DX * DZ
          ELSEIF (ENTRY_FACE .EQ. 4) THEN
C           y=DY face, u=x, v=z
            X0 = UQ * DX
            Y0 = DY
            Z0 = VQ * DZ
            FACE_AREA = DX * DZ
          ELSEIF (ENTRY_FACE .EQ. 5) THEN
C           z=0 face, u=x, v=y
            X0 = UQ * DX
            Y0 = VQ * DY
            Z0 = 0.0
            FACE_AREA = DX * DY
          ELSE
C           z=DZ face (top), u=x, v=y
            X0 = UQ * DX
            Y0 = VQ * DY
            Z0 = DZ
            FACE_AREA = DX * DY
          ENDIF

C         Compute actual path length for this ray (may exit through
C         a different face than the opposite one)
          TMIN = 1.0E30
          IF (ABS(SX) .GT. 1.0E-10) THEN
            IF (SX .GT. 0.0) THEN
              TMIN = MIN(TMIN, (DX - X0)/SX)
            ELSE
              TMIN = MIN(TMIN, -X0/SX)
            ENDIF
          ENDIF
          IF (ABS(SY) .GT. 1.0E-10) THEN
            IF (SY .GT. 0.0) THEN
              TMIN = MIN(TMIN, (DY - Y0)/SY)
            ELSE
              TMIN = MIN(TMIN, -Y0/SY)
            ENDIF
          ENDIF
          IF (ABS(SZ) .GT. 1.0E-10) THEN
            IF (SZ .GT. 0.0) THEN
              TMIN = MIN(TMIN, (DZ - Z0)/SZ)
            ELSE
              TMIN = MIN(TMIN, -Z0/SZ)
            ENDIF
          ENDIF
          SMAX = TMIN
          IF (SMAX .LE. 0.0) CYCLE

C         Compute F0 at entry point (bilinear on entry face from logF0)
C         and tau along this ray from corner values
          XE = X0 + SMAX*SX
          YE = Y0 + SMAX*SY
          ZE = Z0 + SMAX*SZ
C         Bilinear interp of logF0 at entry and exit points
          CALL TRILIN_INTERP_SCALAR(X0/DX, Y0/DY, Z0/DZ,
     .                              LOGF0, LF_ENTRY)
          CALL TRILIN_INTERP_SCALAR(XE/DX, YE/DY, ZE/DZ,
     .                              LOGF0, LF_EXIT)
          F0_ENTRY = EXP(LF_ENTRY)
          TAU = LF_ENTRY - LF_EXIT
          IF (TAU .LT. 0.0) TAU = 0.0

C         Compute In(tau) for n=0..6
          CALL COMPUTE_EIN(TAU, EIN)

C         For each basis function, compute linear coefficients along ray
C         L_k(r(t)) = phi_kx(x0+sx*t) * phi_ky(y0+sy*t) * phi_kz(z0+sz*t)
C         phi_0(x) = 1-x/dx, phi_1(x) = x/dx
C         Along ray: phi_kx(x0+sx*t) = a_kx + b_kx*(t/smax)
          DO K = 1, 8
            IX = MOD(K-1, 2)
            IY = MOD((K-1)/2, 2)
            IZ = (K-1)/4
C           x-component: ix=0 -> phi=1-x/dx, ix=1 -> phi=x/dx
            IF (IX .EQ. 0) THEN
              AX(K) = 1.0 - X0/DX
              BX(K) = -SX*SMAX/DX
            ELSE
              AX(K) = X0/DX
              BX(K) = SX*SMAX/DX
            ENDIF
C           y-component
            IF (IY .EQ. 0) THEN
              AY(K) = 1.0 - Y0/DY
              BY(K) = -SY*SMAX/DY
            ELSE
              AY(K) = Y0/DY
              BY(K) = SY*SMAX/DY
            ENDIF
C           z-component
            IF (IZ .EQ. 0) THEN
              AZ(K) = 1.0 - Z0/DZ
              BZ(K) = -SZ*SMAX/DZ
            ELSE
              AZ(K) = Z0/DZ
              BZ(K) = SZ*SMAX/DZ
            ENDIF
          ENDDO

C         For each (i,j) pair, compute polynomial L_i(u)*L_j(u) and
C         integrate against exp(-tau*u) using EIN
          DO I = 1, 8
C           L_i(u) = (ai_x + bi_x*u)(ai_y + bi_y*u)(ai_z + bi_z*u)
C           = cubic polynomial. Compute coefficients.
            CALL CUBIC_FROM_LINEAR3(AX(I), BX(I), AY(I), BY(I),
     .                              AZ(I), BZ(I), POLY_I)
            DO J = 1, 8
              CALL CUBIC_FROM_LINEAR3(AX(J), BX(J), AY(J), BY(J),
     .                                AZ(J), BZ(J), POLY_J)
C             Multiply two cubics -> degree 6
              CALL POLY_MULT_3_3(POLY_I, POLY_J, POLY_IJ)
C             Integrate: sum_n POLY_IJ(n) * EIN(n) * F0_entry * smax * facearea * wgp
              CONTRIB = 0.0
              DO K = 0, 6
                CONTRIB = CONTRIB + POLY_IJ(K) * EIN(K)
              ENDDO
              MEXP(I,J) = MEXP(I,J)
     .          + WGP * FACE_AREA * SMAX * F0_ENTRY * CONTRIB
            ENDDO
          ENDDO

        ENDDO
      ENDDO

      RETURN
      END


C     ================================================================
      SUBROUTINE COMPUTE_EIN (TAU, EIN)
C     Computes the exponential integrals In(tau) = int_0^1 u^n e^{-tau*u} du
C     for n = 0, 1, ..., 6.
C     Uses Taylor expansion for small tau, recurrence for moderate/large tau.
C
      IMPLICIT NONE
      REAL TAU, EIN(0:6)
Cf2py intent(in) :: TAU
Cf2py intent(out) :: EIN
      INTEGER N
      REAL ETAU, TAU2, TAU3, TAU4

      IF (TAU .LT. 1.0E-4) THEN
C       Taylor expansion: In ≈ 1/(n+1) - tau/(n+2) + tau^2/(2(n+3))
C                                - tau^3/(6(n+4)) + tau^4/(24(n+5))
        TAU2 = TAU*TAU
        TAU3 = TAU2*TAU
        TAU4 = TAU3*TAU
        DO N = 0, 6
          EIN(N) = 1.0/REAL(N+1)
     .           - TAU/REAL(N+2)
     .           + TAU2/(2.0*REAL(N+3))
     .           - TAU3/(6.0*REAL(N+4))
     .           + TAU4/(24.0*REAL(N+5))
        ENDDO
      ELSEIF (TAU .GT. 500.0) THEN
C       Asymptotic: In ≈ n!/tau^(n+1) for large tau
        EIN(0) = 1.0/TAU
        DO N = 1, 6
          EIN(N) = REAL(N)*EIN(N-1)/TAU
        ENDDO
      ELSE
C       Direct computation with recurrence
C       I0 = (1 - e^{-tau})/tau
C       In = (n*I_{n-1} - e^{-tau})/tau
        ETAU = EXP(-TAU)
        EIN(0) = (1.0 - ETAU)/TAU
        DO N = 1, 6
          EIN(N) = (REAL(N)*EIN(N-1) - ETAU)/TAU
        ENDDO
      ENDIF

      RETURN
      END


C     ================================================================
      SUBROUTINE CUBIC_FROM_LINEAR3 (A1, B1, A2, B2, A3, B3, POLY)
C     Computes the cubic polynomial coefficients of:
C       (a1 + b1*u) * (a2 + b2*u) * (a3 + b3*u)
C     Result: POLY(0) + POLY(1)*u + POLY(2)*u^2 + POLY(3)*u^3
C
      IMPLICIT NONE
      REAL A1, B1, A2, B2, A3, B3, POLY(0:3)
Cf2py intent(in) :: A1, B1, A2, B2, A3, B3
Cf2py intent(out) :: POLY
      REAL P12_0, P12_1, P12_2

C     First multiply (a1+b1*u)*(a2+b2*u) -> quadratic
      P12_0 = A1*A2
      P12_1 = A1*B2 + B1*A2
      P12_2 = B1*B2

C     Then multiply by (a3+b3*u) -> cubic
      POLY(0) = P12_0*A3
      POLY(1) = P12_0*B3 + P12_1*A3
      POLY(2) = P12_1*B3 + P12_2*A3
      POLY(3) = P12_2*B3

      RETURN
      END


C     ================================================================
      SUBROUTINE POLY_MULT_3_3 (P, Q, R)
C     Multiplies two cubic polynomials P(0:3) and Q(0:3) to produce
C     a degree-6 polynomial R(0:6).
C
      IMPLICIT NONE
      REAL P(0:3), Q(0:3), R(0:6)
Cf2py intent(in) :: P, Q
Cf2py intent(out) :: R
      INTEGER I, J

      DO I = 0, 6
        R(I) = 0.0
      ENDDO
      DO I = 0, 3
        DO J = 0, 3
          R(I+J) = R(I+J) + P(I)*Q(J)
        ENDDO
      ENDDO

      RETURN
      END


C     ================================================================
      SUBROUTINE TRILIN_INTERP_SCALAR (U, V, W, VALS, RESULT)
C     Trilinear interpolation of 8 corner values at normalized
C     coordinates (u,v,w) in [0,1]^3.
C     Corner ordering: IOCT = 1 + BITX + 2*BITY + 4*BITZ
C
      IMPLICIT NONE
      REAL U, V, W, VALS(8), RESULT
Cf2py intent(in) :: U, V, W, VALS
Cf2py intent(out) :: RESULT

      RESULT = (1.0-U)*(1.0-V)*(1.0-W)*VALS(1)
     .       + U*(1.0-V)*(1.0-W)*VALS(2)
     .       + (1.0-U)*V*(1.0-W)*VALS(3)
     .       + U*V*(1.0-W)*VALS(4)
     .       + (1.0-U)*(1.0-V)*W*VALS(5)
     .       + U*(1.0-V)*W*VALS(6)
     .       + (1.0-U)*V*W*VALS(7)
     .       + U*V*W*VALS(8)

      RETURN
      END


C     ================================================================
      SUBROUTINE BUILD_TRILIN_MASS_MATRIX (DX, DY, DZ, M)
C     Builds the 8x8 trilinear mass matrix for a rectangular cell.
C     M(i,k) = integral over cell of L_i(r) * L_k(r) dV
C     where L_i are trilinear basis functions on the unit cube scaled
C     by DX, DY, DZ.
C
C     The result is the tensor product of three 1D mass matrices:
C       1D mass matrix for interval h: (h/6)*[2 1; 1 2]
C     so the 3D entry for nodes (ix,iy,iz) and (jx,jy,jz) is:
C       (DX/6)(DY/6)(DZ/6) * w_x(ix,jx) * w_y(iy,jy) * w_z(iz,jz)
C     where w(0,0)=2, w(0,1)=w(1,0)=1, w(1,1)=2
C
C     Node ordering follows SHDOM convention:
C       IOCT = 1 + BITX + 2*BITY + 4*BITZ
C       node 1: (0,0,0), node 2: (1,0,0), node 3: (0,1,0), node 4: (1,1,0)
C       node 5: (0,0,1), node 6: (1,0,1), node 7: (0,1,1), node 8: (1,1,1)
C
      IMPLICIT NONE
      REAL DX, DY, DZ, M(8,8)
Cf2py intent(in) :: DX, DY, DZ
Cf2py intent(out) :: M
      REAL WX(2,2), WY(2,2), WZ(2,2)
      REAL FACTOR
      INTEGER I, K, IX, IY, IZ, JX, JY, JZ

C     1D FEM mass matrix weights (for unit interval, *6)
C     w(1,1)=2, w(1,2)=1, w(2,1)=1, w(2,2)=2
      WX(1,1) = 2.0
      WX(1,2) = 1.0
      WX(2,1) = 1.0
      WX(2,2) = 2.0
      WY(1,1) = 2.0
      WY(1,2) = 1.0
      WY(2,1) = 1.0
      WY(2,2) = 2.0
      WZ(1,1) = 2.0
      WZ(1,2) = 1.0
      WZ(2,1) = 1.0
      WZ(2,2) = 2.0

      FACTOR = (DX * DY * DZ) / 216.0

C     Build M using the tensor product structure
C     Node I maps to (IX,IY,IZ) via IOCT = 1+BITX+2*BITY+4*BITZ
      DO I = 1, 8
        IX = MOD(I-1, 2) + 1
        IY = MOD((I-1)/2, 2) + 1
        IZ = (I-1)/4 + 1
        DO K = 1, 8
          JX = MOD(K-1, 2) + 1
          JY = MOD((K-1)/2, 2) + 1
          JZ = (K-1)/4 + 1
          M(I,K) = FACTOR * WX(IX,JX) * WY(IY,JY) * WZ(IZ,JZ)
        ENDDO
      ENDDO

      RETURN
      END
