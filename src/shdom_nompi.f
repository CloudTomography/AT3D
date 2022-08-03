C     This file contains fortran subroutines that perform the
C     original SHDOM solution written by Frank Evans.
C     https://nit.coloradolinux.com/~evans/shdom.html
C     See shdom.txt for documentation.
C     Many of these subroutines have been modified for use in pyshdom by
C     Aviad Levis, Technion Institute of Technology, 2019 and
C     Jesse Loveridge, University of Illinois at Urbana-Champaign, 2020-2021.
C     Subroutines only have f2py directives added.
C Dummy routines to allow parallel SHDOM to compile without MPI.

      SUBROUTINE START_MPI (MASTERPROC)
      IMPLICIT NONE
      LOGICAL MASTERPROC
Cf2py intent(out) MASTERPROC
      MASTERPROC = .TRUE.
      RETURN
      END



      SUBROUTINE MAP_SHDOM_MPI (BCFLAG, NPX, NPY, NX, NY, NZ, DELX,
     .                 DELY, PROPFILE, XSTART, YSTART, RUNNAME)
      IMPLICIT NONE
      INTEGER BCFLAG, NPX, NPY, NX, NY, NZ
Cf2py intent(in, out) :: BCFLAG, NPX, NPY, NX, NY, NZ
      REAL    DELX, DELY
Cf2py intent(in) :: DELX, DELY
      CHARACTER PROPFILE*64, RUNNAME*64
Cf2py intent(in) :: PROPFILE, RUNNAME
      REAL    XSTART, YSTART
Cf2py intent(out) :: XSTART, YSTART
      XSTART = 0.0
      YSTART = 0.0
      RETURN
      END




      SUBROUTINE BROADCAST_USER_INPUT (PROPFILE, SFCFILE, CKDFILE,
     .                 INSAVEFILE, OUTSAVEFILE,
     .                 NX, NY, NZ, NMU, NPHI, BCFLAG,
     .                 IPFLAG, KDIST, DELTAM, GRIDTYPE,
     .                 SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .                 GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,
     .                 ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,
     .                 MAXOUT,MAXPAR,NUMOUT,OUTTYPES,OUTPARMS,OUTFILES)
      INTEGER NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
      INTEGER MAXOUT, MAXPAR, MAXITER, NUMOUT
      LOGICAL KDIST, DELTAM, ACCELFLAG
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
      REAL    SOLACC, SPLITACC, SHACC, OUTPARMS(MAXPAR,MAXOUT)
      CHARACTER SRCTYPE*1, GRIDTYPE*1, UNITS*1, OUTTYPES(*)*1
      CHARACTER PROPFILE*64, SFCFILE*64, CKDFILE*64
      CHARACTER OUTFILES(*)*64, INSAVEFILE*64, OUTSAVEFILE*64

      RETURN
      END



      SUBROUTINE BROADCAST_PROPERTY_SIZE (NPX, NPY, NPZ, DELX, DELY,
     .                                    NUMPHASE, MAXLEG, MAXPGL)
      INTEGER NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
Cf2py intent(in, out) :: NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
      REAL    DELX, DELY
Cf2py intent(in, out) :: DELX, DELY
      RETURN
      END


      SUBROUTINE SCATTER_PROPERTIES (NPXT,NPYT, NPX,NPY,NPZ, NLEG,
     .                    NUMPHASE, ZLEVELS, MAXASYM, TEMPPT,
     .                    EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT,
     .                    TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP)
      INTEGER   NPXT, NPYT, NPX, NPY, NPZ, NUMPHASE, NLEG
Cf2py intent(in) :: NPXT, NPYT, NPX, NPY, NPZ, NUMPHASE
Cf2py intent(in, out) :: NLEG
      REAL      ZLEVELS(NPZ), MAXASYM
Cf2py intent(in, out) :: ZLEVELS, MAXASYM
      REAL      TEMPPT(NPZ,NPYT,NPXT), EXTINCTPT(NPZ,NPYT,NPXT)
Cf2py intent(in) :: TEMPPT, EXTINCTPT
      REAL      ALBEDOPT(NPZ,NPYT,NPXT), LEGENPT(*)
Cf2py intent(in) :: ALBEDOPT, LEGENPT
      INTEGER IPHASEPT(NPZ,NPYT,NPXT)
Cf2py intent(in) :: IPHASEPT
      REAL     TEMPP(NPZ,NPY,NPX), EXTINCTP(NPZ,NPY,NPX)
Cf2py intent(out) :: TEMPP, EXTINCTP
      REAL     ALBEDOP(NPZ,NPY,NPX), LEGENP(*)
Cf2py intent(out) :: ALBEDOP, LEGENP
      INTEGER IPHASEP(NPZ,NPY,*)
Cf2py intent(out) :: IPHASEP

      RETURN
      END


      SUBROUTINE READ_BROADCAST_MEM_PARAMS (MAX_TOTAL_MB,
     .                ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR,
     .                CELL_TO_POINT_RATIO, RUNNAME)
      REAL MAX_TOTAL_MB, ADAPT_GRID_FACTOR
Cf2py intent(in, out) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
      REAL NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
Cf2py intent(in, out) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
      CHARACTER(*) RUNNAME
Cf2py intent(in) :: RUNNAME

      RETURN
      END


      SUBROUTINE BROADCAST_SURFACE_SIZE (MAXSFCPTS, MAXSFCPARS)
      INTEGER MAXSFCPTS, MAXSFCPARS

      RETURN
      END


      SUBROUTINE BROADCAST_SURFACE_PARMS (SFCTYPE, NXSFC,NYSFC,NSFCPAR,
     .                  DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)
      CHARACTER SFCTYPE*2
Cf2py intent(in, out) :: SFCTYPE
      INTEGER   NXSFC, NYSFC, NSFCPAR
Cf2py intent(in, out) :: NXSFC, NYSFC, NSFCPAR
      REAL      DELXSFC, DELYSFC, SFCPARMS(*), GNDTEMP, GNDALBEDO
Cf2py intent(in, out) :: DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO

      RETURN
      END




      SUBROUTINE BROADCAST_KDIST_SIZE (NG, NZCKD)
      INTEGER NG, NZCKD

      RETURN
      END


      SUBROUTINE BROADCAST_KDIST_PARMS (SOLFLUX, NG, DELG,
     .                                  NZCKD, ZCKD, KABS)
      INTEGER NG, NZCKD
      REAL    SOLFLUX, DELG(NG), ZCKD(*), KABS(NZCKD,NG)

      RETURN
      END




      SUBROUTINE GATHER_OUTPUT (FLUX_OUT, FLUXDIV_OUT, SH_OUT,
     .                   IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT,
     .                   SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,
     .                   NPX, NPY, DELX, DELY, XALLGRID, YALLGRID,
     .                   ALLFLUXES, ALLFLUXDIV, ALLSHTERMS,
     .                   NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)
      LOGICAL FLUX_OUT, FLUXDIV_OUT, SH_OUT
Cf2py intent(in) :: FLUX_OUT, FLUXDIV_OUT, SH_OUT
      INTEGER IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, NPX, NPY
Cf2py intent(in) :: IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, NPX, NPY
      REAL    DELX, DELY
Cf2py intent(in) :: DELX, DELY
      REAL    SUMFLUXES(2,*), SUMDIRFLUX(*)
Cf2py intent(in) :: SUMFLUXES, SUMDIRFLUX
      REAL    SUMFLUXDIV(*), SUMSHTERMS(NSHOUT,*)
Cf2py intent(in) :: SUMFLUXDIV, SUMSHTERMS
      REAL    XALLGRID(*), YALLGRID(*)
Cf2py intent(out) :: XALLGRID, YALLGRID
      REAL    ALLFLUXES(3,NZ,NYT,NXT), ALLFLUXDIV(NZ,NYT,NXT)
Cf2py intent(out) :: ALLFLUXES, ALLFLUXDIV
      REAL    ALLSHTERMS(NSHOUT,NZ,NYT,NXT)
Cf2py intent(out) :: ALLSHTERMS
      INTEGER NCELLS, NPTS, NSH
Cf2py intent(in) :: NCELLS, NPTS, NSH
      INTEGER NCELLSTOT, NPTSTOT, NSHTOT
Cf2py intent(in) :: NCELLSTOT, NPTSTOT, NSHTOT

      RETURN
      END


      REAL FUNCTION SUM_CPU_TIME (cpuTimes)
        real cpuTimes
        SUM_CPU_TIME = cpuTimes
        RETURN
      END


      SUBROUTINE TOTAL_ALBEDO_MAX (ALBMAX)
      REAL ALBMAX
      RETURN
      END


      SUBROUTINE UNIFY_SPLITTING (DOSPLIT, STARTSPLITACC)
      LOGICAL DOSPLIT
      REAL    STARTSPLITACC
      RETURN
      END

      SUBROUTINE TOTAL_SPLITCRIT_MAX (SPLITCRIT)
      REAL SPLITCRIT
      RETURN
      END


      SUBROUTINE MAKE_DIRECT_PAR (SPT, NPTS, BCFLAG, IPFLAG, DELTAM,
     .                ML, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,
     .                NX, XGRID, NY, YGRID,  DIRFLUX)
      IMPLICIT NONE
      INTEGER SPT, NPTS, BCFLAG, IPFLAG, ML, NLEG, NX, NY
      LOGICAL DELTAM
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
      REAL    GRIDPOS(3,*), XGRID(*), YGRID(*)
      REAL    DIRFLUX(*)

      RETURN
      END




      SUBROUTINE FIND_BOUNDARY_POINTS (BCFLAG, IPFLAG, NPTS, SWEEPORD,
     .               GRIDPTR, GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)
      IMPLICIT NONE
      INTEGER BCFLAG, IPFLAG, NPTS, SWEEPORD(NPTS,*), GRIDPTR(8,*)
      INTEGER NX, NY, NZ
      REAL GRIDPOS(3,*), XGRID(*), YGRID(*), ZGRID(*)

      RETURN
      END



      SUBROUTINE CALC_BOUNDARY_RADIANCES (BCFLAG, IPFLAG, JOCT, IZ,
     .                           NX, NY, NZ, XGRID, YGRID, ZGRID,
     .                           NA, NPTS, NCELLS, GRIDPTR,
     .                           NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS,
     .                           MU, PHI, EXTINCT, SOURCE,
     .                           KANG, GRIDRAD)
      INTEGER BCFLAG, IPFLAG, JOCT, IZ
      INTEGER NX, NY, NZ, NA, NPTS, NCELLS, KANG
      REAL    XGRID(*), YGRID(*), ZGRID(*)
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
      INTEGER*2 CELLFLAGS(*)
      REAL    GRIDPOS(3,*), MU, PHI
      REAL    EXTINCT(*), SOURCE(NA,*), GRIDRAD(*)

      RETURN
      END



      SUBROUTINE COMPUTE_RADIANCE_PAR (NX, NY, NZ, NPTS, NCELLS,
     .             ML, MM, NCS, NLEG, NUMPHASE,
     .             NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,
     .             BCFLAG, XDOMAIN, YDOMAIN, IPFLAG,
     .             SRCTYPE, DELTAM, SOLARMU, SOLARAZ,
     .             SFCTYPE, NSFCPAR, SFCGRIDPARMS,
     .             MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,
     .             GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,
     .             SHPTR, SOURCE, SOURCE1, GRIDRAD,
     .             OUTPARMS,  NRAD, NANGOUT, RADOUT)
      IMPLICIT NONE
      INTEGER NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS, NANGOUT
Cf2py intent(in) :: NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS, NANGOUT
      INTEGER ML, MM, NCS, NLEG, NUMPHASE
Cf2py intent(in) :: ML, MM, NCS, NLEG, NUMPHASE
      INTEGER NMU, NPHI0MAX, NPHI0(*), NRAD
Cf2py intent(in) :: NMU, NPHI0MAX, NPHI0
Cf2py intent(in, out) :: NRAD
      INTEGER MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
Cf2py intent(in) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER SHPTR(*), BCPTR(MAXNBC,2)
Cf2py intent(in) :: SHPTR, BCPTR
      INTEGER*2 CELLFLAGS(NCELLS)
      INTEGER   IPHASE(NPTS)
Cf2py intent(in) :: CELLFLAGS, IPHASE
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    SOLARMU, SOLARAZ
Cf2py intent(in) :: SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
Cf2py intent(in) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN
      REAL    MU(NMU), PHI(NMU,NPHI0MAX), WTDO(NMU,NPHI0MAX)
Cf2py intent(in) :: MU, PHI, WTDO
      REAL    XDOMAIN, YDOMAIN, XGRID(*), YGRID(*), ZGRID(*)
Cf2py intent(in) :: XDOMAIN, YDOMAIN, XGRID, YGRID, ZGRID
      REAL    GRIDPOS(3,NPTS)
Cf2py intent(in) :: GRIDPOS
      REAL    SFCGRIDPARMS(*), BCRAD(*)
Cf2py intent(in) :: SFCGRIDPARMS, BCRAD
      REAL    EXTINCT(*), ALBEDO(*), LEGEN(0:NLEG,*)
Cf2py intent(in) :: EXTINCT, ALBEDO, LEGEN
      REAL    DIRFLUX(*), FLUXES(2,*), SOURCE(*)
Cf2py intent(in) :: DIRFLUX, FLUXES, SOURCE
      REAL    SOURCE1(*), GRIDRAD(*), OUTPARMS(*)
Cf2py intent(in, out) :: SOURCE1, GRIDRAD
Cf2py intent(in) :: OUTPARMS
      REAL    RADOUT(NANGOUT,*)
Cf2py intent(in, out) :: RADOUT
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1
Cf2py intent(in) :: SRCTYPE, SFCTYPE, UNITS

      RETURN
      END



      SUBROUTINE CALC_ACCEL_SOLCRIT (DOACCEL, DELJDOT, DELJOLD, DELJNEW,
     .                               JNORM, ACCELPAR, SOLCRIT, A)
C     Calculates the acceleration parameter and solution criterion from
C     the delta source function vector dot  products.
      IMPLICIT NONE
      LOGICAL DOACCEL
      REAL    DELJDOT, DELJOLD, DELJNEW, JNORM, ACCELPAR, SOLCRIT
      REAL    R, THETA, A

C       Accelerate if desired, didn't last time, and things are converging.
      IF (DOACCEL .AND. A .EQ. 0.0 .AND. DELJNEW .LT. DELJOLD) THEN
C       Compute the acceleration extrapolation factor and apply it.
        R = SQRT(DELJNEW/DELJOLD)
        THETA = ACOS(DELJDOT/SQRT(DELJOLD*DELJNEW))
        A = (1 - R*COS(THETA) + R**(1+0.5*3.14159/THETA))
     .         /(1 + R**2  - 2*R*COS(THETA))  - 1.0
        A = MIN(10.0,MAX(0.0,A))
C         WRITE (*,'(1X,A,3(1X,F7.3))') '! Acceleration: ', A,R,THETA
      ELSE
        A = 0.0
      ENDIF
      ACCELPAR = A

      IF (JNORM .GT. 0.0) THEN
        SOLCRIT = SQRT(DELJNEW/JNORM)
      ELSE IF (DELJNEW .EQ. 0.0) THEN
        SOLCRIT = 0.0
      ENDIF
      RETURN
      END




      SUBROUTINE END_SHDOM_MPI (NPTS, GRIDPOS, NPX,NPY, XSTART,YSTART,
     .                          DELX, DELY, NPXT, NPYT)
      IMPLICIT NONE
      INTEGER NPTS, NPX, NPY, NPXT, NPYT
Cf2py intent(in) :: NPTS, NPX, NPY, NPXT, NPYT
      REAL    GRIDPOS(3,*), XSTART, YSTART, DELX, DELY
Cf2py intent(in) :: GRIDPOS, XSTART, YSTART, DELX, DELY

      RETURN
      END


      SUBROUTINE ABORT_SHDOM_MPI (ERRSTR)
      CHARACTER(*) ERRSTR

      WRITE (6,*) ERRSTR
      stop
      END
