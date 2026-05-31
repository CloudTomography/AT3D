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
     .           RSHPTR_PF, SHPTR_PF, RADIANCE_PF, SOURCE_PF,
     .           DIRFLUX_PF,
     .           SOLARMU, SOLARAZ, YLMSUN,
     .           DET_STOKES, CELLWISE_DET,
     .           MU_PF_CELL, PHI_PF_CELL,
     .           POINT_SENS, IERR, ERRMSG)
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
      INTEGER RSHPTR_PF(NPTS+1), SHPTR_PF(NPTS+1)
Cf2py intent(in) :: RSHPTR_PF, SHPTR_PF
      REAL    RADIANCE_PF(NSTOKES,*), SOURCE_PF(NSTOKES,*)
Cf2py intent(in) :: RADIANCE_PF, SOURCE_PF
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
      REAL    POINT_SENS(NPTS)
Cf2py intent(out) :: POINT_SENS
      INTEGER IERR
      CHARACTER ERRMSG*600
Cf2py intent(out) :: IERR, ERRMSG

C     --- Local variables ---
      INTEGER IC, I, K, J, L, S, IP(8)
      INTEGER NR_FWD, NR_PF, NS_FWD, NMIN
      REAL    DX, DY, DZ, VOL
      REAL    M(8,8)
      REAL    P_AA(8), P_JA(8,8)
      REAL    IADJ_SOL(8), IFWD_DET(8), JFWD_DET(8)
      REAL    A_IC(8), B_IC(8), C_IC(8), D_IC(8), E_IC(8)
      REAL    SUMC, SUMD, SUME, ACC
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
      DO I = 1, NPTS
        POINT_SENS(I) = 0.0
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

C       Build the 8x8 trilinear mass matrix M(i,k) = integral(Li*Lk*dV)
        CALL BUILD_TRILIN_MASS_MATRIX(DX, DY, DZ, M)

C       Compute YLM at detector direction (per-cell if cellwise)
        IF (CELLWISE_DET) THEN
          CALL YLMALL (.FALSE., -MU_PF_CELL(IC),
     .                 PHI_PF_CELL(IC) + PI_VAL,
     .                 ML, MM, NSTLEG, YLM_DET)
        ENDIF

C       ============================================================
C       Term A: diffuse-diffuse Parseval (inline sign flip)
C       Parseval: integral_S2 sum_s I_d(s) * I_adj(s) dOmega
C         = sum_J sum_s RAD(s,J) * SIGN_ADJ(J) * RSTOKES(s) * RAD_PF(s,J)
C       For s=1,2,3: RSTOKES=+1; for s=4: RSTOKES=-1 (Stokes V flip)
C       ============================================================
        DO K = 1, 8
          NR_FWD = RSHPTR(IP(K)+1) - RSHPTR(IP(K))
          NR_PF  = RSHPTR_PF(IP(K)+1) - RSHPTR_PF(IP(K))
          NMIN   = MIN(NR_FWD, NR_PF)
          ACC = 0.0
          IF (NSTOKES .LT. 4) THEN
C           No Stokes V: simple dot product with sign flip
            DO J = 1, NMIN
              DO S = 1, NSTOKES
                ACC = ACC + RADIANCE(S, RSHPTR(IP(K))+J)
     .                    * SIGN_ADJ(J)
     .                    * RADIANCE_PF(S, RSHPTR_PF(IP(K))+J)
              ENDDO
            ENDDO
          ELSE
C           NSTOKES=4: s=4 gets extra -1 for Stokes V direction reversal
            DO J = 1, NMIN
              DO S = 1, 3
                ACC = ACC + RADIANCE(S, RSHPTR(IP(K))+J)
     .                    * SIGN_ADJ(J)
     .                    * RADIANCE_PF(S, RSHPTR_PF(IP(K))+J)
              ENDDO
              ACC = ACC - RADIANCE(4, RSHPTR(IP(K))+J)
     .                  * SIGN_ADJ(J)
     .                  * RADIANCE_PF(4, RSHPTR_PF(IP(K))+J)
            ENDDO
          ENDIF
          P_AA(K) = ACC
        ENDDO

C       A_ic(i) = sum_k M(i,k) * P_AA(k)
        DO I = 1, 8
          ACC = 0.0
          DO K = 1, 8
            ACC = ACC + M(I,K) * P_AA(K)
          ENDDO
          A_IC(I) = ACC
        ENDDO

C       ============================================================
C       Term B: cross-Parseval (forward SOURCE at node i,
C               adjoint radiance at corner k)
C       Same Parseval identity applies for the angular integral.
C       ============================================================
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
C         B_ic(i) = - sum_k M(i,k) * P_JA(i,k)
          ACC = 0.0
          DO K = 1, 8
            ACC = ACC + M(I,K) * P_JA(I,K)
          ENDDO
          B_IC(I) = -ACC
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
          NR_PF = RSHPTR_PF(IP(K)+1) - RSHPTR_PF(IP(K))
          ACC = 0.0
          DO J = 1, NR_PF
            ACC = ACC + SIGN_ADJ(J) * YLMSUN(1,J)
     .                * RADIANCE_PF(1, RSHPTR_PF(IP(K))+J)
          ENDDO
          IADJ_SOL(K) = ACC

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

C       Combine C, D, E with mass matrix and direct beams
        DO I = 1, 8
          SUMC = 0.0
          SUMD = 0.0
          SUME = 0.0
          DO K = 1, 8
            SUMC = SUMC + M(I,K) * IADJ_SOL(K) * DIRFLUX(IP(K))
            SUMD = SUMD + M(I,K) * DIRFLUX_PF(IP(K)) * IFWD_DET(K)
            SUME = SUME + M(I,K) * DIRFLUX_PF(IP(K)) * JFWD_DET(K)
          ENDDO
          C_IC(I) = SUMC
          D_IC(I) = SUMD
          E_IC(I) = -SUME
        ENDDO

C       Scatter per-corner contributions to global point sensitivity
        DO I = 1, 8
          POINT_SENS(IP(I)) = POINT_SENS(IP(I))
     .        + A_IC(I) + B_IC(I) + C_IC(I) + D_IC(I) + E_IC(I)
        ENDDO

      ENDDO
C     End of leaf cell loop

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
