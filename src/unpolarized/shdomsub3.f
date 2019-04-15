      SUBROUTINE NAMELIST_INPUT 
     .               (dRUNNAME, dPROPFILE, dSFCFILE, dCKDFILE, 
     .                dINSAVEFILE, dOUTSAVEFILE,
     .                dNX, dNY, dNZ, dNMU, dNPHI, dBCFLAG, dIPFLAG, 
     .                dKDIST, dDELTAM, dGRIDTYPE, 
     .                dSRCTYPE, dSOLARFLUX,dSOLARMU,dSOLARAZ, dSKYRAD,
     .                dGNDTEMP, dGNDALBEDO, dUNITS, dWAVENO, dWAVELEN, 
     .                dACCELFLAG, dSOLACC, dMAXITER, dSPLITACC,dSHACC,
     .       dMAXOUT, dMAXPAR, dNUMOUT, dOUTTYPES, dOUTPARMS, dOUTFILES,
     .                dOutFileNC, dMAX_TOTAL_MB, dADAPT_GRID_FACTOR,
     .                dNUM_SH_TERM_FACTOR, dCELL_TO_POINT_RATIO)
C       Obtains the input parameters for the program from a namelist 
C     read from stdin.  See the overall program documentation for the 
C     list of input parameters.  The subroutine dummy variables can't
C     be in the namelist, hence the two sets of variables.
      INTEGER dNX, dNY, dNZ, dNMU, dNPHI, dBCFLAG, dIPFLAG, dMAXPAR
      INTEGER dMAXITER, dNUMOUT, dMAXOUT
      LOGICAL dKDIST, dDELTAM, dACCELFLAG
      REAL    dSOLARFLUX, dSOLARMU, dSOLARAZ
      REAL    dGNDTEMP, dGNDALBEDO, dSKYRAD
      REAL    dWAVENO(2), dWAVELEN, dSOLACC, dSPLITACC, dSHACC
      REAL    dOUTPARMS(dMAXPAR,dMAXOUT)
      REAL    dMAX_TOTAL_MB, dADAPT_GRID_FACTOR
      REAL    dNUM_SH_TERM_FACTOR, dCELL_TO_POINT_RATIO
      CHARACTER dSRCTYPE*1, dGRIDTYPE*1, dUNITS*1
      CHARACTER dOUTTYPES(*)*1, dOUTFILES(*)*64
      CHARACTER*64 dOutFileNC, dRUNNAME
      CHARACTER*64 dPROPFILE, dSFCFILE, dCKDFILE
      CHARACTER*64 dINSAVEFILE, dOUTSAVEFILE
      INTEGER I, J, N

      INTEGER MAXOUT, MAXPAR
      PARAMETER (MAXOUT=20, MAXPAR=400)
      INTEGER NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
      INTEGER MAXITER, NUMOUT
      LOGICAL DELTAM, ACCELFLAG
      REAL  SOLARFLUX, SOLARMU, SOLARAZ
      REAL  GNDTEMP, GNDALBEDO, SKYRAD
      REAL  WAVENO(2), WAVELEN, SOLACC, SPLITACC, SHACC
      REAL  OUTPARMS(MAXPAR,MAXOUT)
      CHARACTER SRCTYPE*1, GRIDTYPE*1, UNITS*1
      CHARACTER OUTTYPES(MAXOUT)*1, OUTFILES(MAXOUT)*64
      CHARACTER*64 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
      CHARACTER*64 OutFileNC, RUNNAME
      NAMELIST /SHDOMINPUT/ RUNNAME, PROPFILE, SFCFILE, CKDFILE, 
     .                 INSAVEFILE, OUTSAVEFILE,
     .                 NX, NY, NZ, NMU, NPHI,
     .                 BCFLAG, IPFLAG, DELTAM, GRIDTYPE, 
     .                 SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .                 GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN, 
     .                 ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,
     .                 NUMOUT, OUTTYPES, OUTPARMS, OUTFILES, OutFileNC,
     .                 MAX_TOTAL_MB, ADAPT_GRID_FACTOR, 
     .                 NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO

      CHARACTER*64 NMLFILE

      IF (dMAXOUT.NE.MAXOUT) STOP 'NAMELIST_INPUT: MAXOUT not the same.'
      IF (dMAXPAR.NE.MAXPAR) STOP 'NAMELIST_INPUT: MAXPAR not the same.'

      SOLARMU = -1.0     
      GNDALBEDO = 0.0 
      WAVELEN = 1.0
      WAVENO(1) = 10000.
      WAVENO(2) = 10001.

C         Read in the namelist
      WRITE (*,'(1X,A)') 'Namelist file name'
      READ (*,'(A)') NMLFILE
      OPEN (UNIT=5, FILE=NMLFILE, STATUS='OLD')
      READ (5,NML=SHDOMINPUT)
      CLOSE (5)

C         Make the assignments to the dummy variables
      dRUNNAME = RUNNAME
      dPROPFILE = PROPFILE
      dSFCFILE = SFCFILE
      dCKDFILE = CKDFILE
      dKDIST = CKDFILE(1:2) .NE. 'NO'
      dINSAVEFILE = INSAVEFILE
      dOUTSAVEFILE = OUTSAVEFILE
      dNX = NX
      dNY = NY
      dNZ = NZ
      dNMU = NMU
      dNPHI = NPHI
      dIPFLAG = IPFLAG
      dBCFLAG = BCFLAG
      dDELTAM = DELTAM
      dGRIDTYPE = GRIDTYPE

      dSRCTYPE = SRCTYPE
      IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
        dSOLARFLUX = SOLARFLUX
        dSOLARMU = -ABS(SOLARMU)
        dSOLARAZ = SOLARAZ*ACOS(-1.0)/180.0
      ELSE
        dSOLARFLUX = 0.0
        dSOLARMU = 1.0
        dSOLARAZ = 0.0
      ENDIF
      dSKYRAD = SKYRAD
      dGNDTEMP = GNDTEMP
      dGNDALBEDO = GNDALBEDO

      IF (dKDIST) THEN
        dUNITS = 'B'
        dWAVENO(1) = WAVENO(1)
        dWAVENO(2) = WAVENO(2)
        dWAVELEN = 1.0E4/(0.5*(WAVENO(1)+WAVENO(2)))
      ELSE 
        dUNITS = UNITS
        dWAVELEN = WAVELEN
        IF (SRCTYPE .NE. 'T') THEN
          dUNITS = 'R'
        ENDIF
      ENDIF

      dSPLITACC = SPLITACC
      dSHACC = SHACC
      dACCELFLAG = ACCELFLAG
      dSOLACC = SOLACC
      dMAXITER = MAXITER

      dNUMOUT = NUMOUT
      IF (NUMOUT .GT. MAXOUT) STOP 'NAMELIST_INPUT: MAXOUT exceeded'
      DO I = 1, NUMOUT
        dOUTTYPES(I) = OUTTYPES(I)
        IF (OUTTYPES(I)(1:1) .EQ. 'R') THEN
          N = NINT(OUTPARMS(6,I))
          DO J = 1, MIN(2*N+6,MAXPAR)
            dOUTPARMS(J,I) = OUTPARMS(J,I)
          ENDDO
        ELSE IF (OUTTYPES(I)(1:1) .EQ. 'V') THEN
          DO J = 1, 13
            dOUTPARMS(J,I) = OUTPARMS(J,I)
          ENDDO
        ELSE IF (OUTTYPES(I)(1:1) .EQ. 'F') THEN
          dOUTPARMS(1,I) = OUTPARMS(1,I)
          IF (NINT(OUTPARMS(1,I)) .EQ. 2) THEN
            DO J = 2, 4
              dOUTPARMS(J,I) = OUTPARMS(J,I)
            ENDDO
          ENDIF
        ELSE IF (OUTTYPES(I)(1:1) .EQ. 'H') THEN
          dOUTPARMS(1,I) = OUTPARMS(1,I)
        ELSE IF (OUTTYPES(I)(1:1) .EQ. 'S') THEN
          dOUTPARMS(1,I) = OUTPARMS(1,I)
        ELSE IF (OUTTYPES(I)(1:1) .EQ. 'J') THEN
          dOUTPARMS(1,I) = OUTPARMS(1,I)
          dOUTPARMS(2,I) = OUTPARMS(2,I)
          dOUTPARMS(3,I) = OUTPARMS(3,I)
        ELSE IF (OUTTYPES(I)(1:1) .EQ. 'M') THEN
          dOUTPARMS(1,I) = OUTPARMS(1,I)
        ENDIF
        dOUTFILES(I) = OUTFILES(I)
      ENDDO
      dOutFileNC = OutFileNC

      dMAX_TOTAL_MB = MAX_TOTAL_MB
      dADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR
      dNUM_SH_TERM_FACTOR = NUM_SH_TERM_FACTOR
      dCELL_TO_POINT_RATIO = CELL_TO_POINT_RATIO

      RETURN
      END
 


      SUBROUTINE CHECK_INPUT_PARMETERS (NMU, NPHI, DELTAM, MAXASYM,
     .               GRIDTYPE, SRCTYPE, UNITS, SOLARFLUX, SOLARMU, 
     .               WAVENO, WAVELEN, GNDTEMP, GNDALBEDO, 
     .               SPLITACC, SHACC)
C       Checks for validity of the input parameters. Gives warning for
C     parameters that might not be correct. 
C     Check the letter flags.
C     Delta-M should be used for solar transfer with forward peaked 
C       phase functions.
C     Check for reasonable ranges on SPLITACC and SHACC.
C     Check for reasonable range of GNDTEMP and valid range of GNDALBEDO.
C     Check for SOLARMU=0.
C     Check for NPHI<NMU.
      INTEGER NMU, NPHI
Cf2py intent(in) :: NMU, NPHI
      LOGICAL DELTAM
Cf2py intent(in) :: DELTAM
      REAL    MAXASYM, SOLARFLUX, SOLARMU, WAVENO(2), WAVELEN
Cf2py intent(in) :: MAXASYM, SOLARFLUX, SOLARMU, WAVENO(2), WAVELEN
      REAL    GNDTEMP, GNDALBEDO, SPLITACC, SHACC
Cf2py intent(in) :: GNDTEMP, GNDALBEDO, SPLITACC, SHACC
      CHARACTER*1 GRIDTYPE, SRCTYPE, UNITS
Cf2py intent(in) :: GRIDTYPE, SRCTYPE, UNITS
      REAL    FLUX

      IF (SRCTYPE .NE.'S' .AND. SRCTYPE .NE.'T' .AND. SRCTYPE .NE.'B')
     .    STOP 'CHECK_INPUT_PARMETERS: Illegal source type'
      IF (UNITS .NE. 'R' .AND. UNITS .NE. 'T' .AND. UNITS .NE. 'B')
     .    STOP 'CHECK_INPUT_PARMETERS: Illegal units'
      IF (GRIDTYPE.NE.'E' .AND. GRIDTYPE.NE.'P' .AND. GRIDTYPE.NE.'F')
     .    STOP 'CHECK_INPUT_PARMETERS: Illegal grid type'

      IF (SRCTYPE .NE.'T' .AND. SOLARMU .EQ. 0.0) THEN
        WRITE (*,*) 'CHECK_INPUT_PARMETERS: cosine solar zenith angle'
        WRITE (*,*) ' (SOLARMU) cannot be zero.'
        STOP
      ENDIF

      IF (.NOT.DELTAM .AND. SRCTYPE.NE.'T' .AND. MAXASYM.GT.0.5) THEN
        WRITE (*,*) 'CHECK_INPUT_PARMETERS: Delta-M should be used'
        WRITE (*,*) ' (DELTAM=.TRUE.) for solar radiative transfer'
        WRITE (*,*) ' with highly peaked phase functions.'
        WRITE (*,*)
      ENDIF

C         Check for SPLITACC and SHACC magnitudes
      IF (SRCTYPE .EQ. 'T' .OR. SRCTYPE .EQ. 'B') THEN
        CALL PLANCK_FUNCTION (GNDTEMP,UNITS,WAVENO,WAVELEN,FLUX)
        FLUX = 3.14*FLUX
      ELSE
        FLUX = 0.0
      ENDIF
      IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
        FLUX = FLUX + SOLARFLUX
      ENDIF
      IF (SPLITACC .GT. 0.0 .AND. 
     .    (SPLITACC .LT. 0.001*FLUX .OR. SPLITACC .GT. 0.1*FLUX)) THEN
        WRITE (*,*)'CHECK_INPUT_PARMETERS: Splitting accuracy parameter'
        WRITE (*,*) ' is not a fraction, but scales with fluxes'
        WRITE (*,*)
      ENDIF
      IF (SHACC .GT. 0.0 .AND. SHACC .GT. 0.03*FLUX) THEN
        WRITE (*,*)'CHECK_INPUT_PARMETERS: Spherical harmonic accuracy'
        WRITE (*,*) ' (SHACC) is not a fraction, but scales with fluxes'
        WRITE (*,*)
      ENDIF

      IF (SRCTYPE .NE. 'S' .AND. 
     .    (GNDTEMP .LT. 150. .OR. GNDTEMP .GT. 350.)) THEN
        WRITE (*,*)'CHECK_INPUT_PARMETERS: Warning - '
        WRITE (*,*) '  ground temperature out of Earth range.'
        WRITE (*,*)
      ENDIF
      IF (GNDALBEDO .LT. 0.0 .OR. GNDALBEDO .GT. 1.0) THEN
        WRITE (*,*)'CHECK_INPUT_PARMETERS: Ground albedo must be 0 to 1'
        STOP
      ENDIF

      IF (NPHI .LT. NMU) THEN
        WRITE (*,*)'CHECK_INPUT_PARMETERS: Usually want NPHI >= NMU'
        WRITE (*,*)
      ENDIF

      RETURN
      END

 
 


 
      SUBROUTINE READ_PROPERTY_SIZE (PROPFILE, NLEG, NPX, NPY, NPZ, 
     .                            NUMPHASE, MAXLEG, MAXPGL, DELX, DELY)
C       Reads parts of the property file to get the maximum array sizes
C     needed for allocatable arrays.  For extinction only and tabulated
C     phase function formats, only the header is read, while for the
C     standard format, the whole file must be read to determine MAXLEG.
      INTEGER NPX, NPY, NPZ
Cf2py intent(out) NPX, NPY, NPZ
      INTEGER NUMPHASE
      INTEGER NLEG, MAXLEG, MAXPGL
Cf2py intent(out) NUMPHASE, MAXLEG, MAXPGL
      REAL    DELX, DELY
Cf2py intent(out) DELX, DELY
      CHARACTER PROPFILE*64
Cf2py intent(in) PROPFILE
      INTEGER NUML, IX, IY, IZ, I, K, L
      REAL    ZLEVELS, TMP, EXT, ALB, CHI
      CHARACTER PROPTYPE*1
 
C          Open the file, figure out the type, and get the grid size
      OPEN (UNIT=1, FILE=PROPFILE, STATUS='OLD')
      READ (1,'(A1)') PROPTYPE
      IF (PROPTYPE .NE. 'E' .AND. PROPTYPE .NE. 'T')  REWIND (1)
      READ (1,*) NPX, NPY, NPZ
      READ (1,*) DELX, DELY, (ZLEVELS, K=1,NPZ)
 
C           Property file type E is for extinction only format
      IF (PROPTYPE .EQ. 'E') THEN
        NUMPHASE = 1
        READ (1,*) (TMP, IZ=1,NPZ)
        READ (1,*) ALB, NUML
        MAXLEG = MAX(NLEG,NUML)
        MAXPGL = MAXLEG

C           Property file type T is for tabulated phase function format
      ELSE IF (PROPTYPE .EQ. 'T') THEN
        MAXLEG = NLEG
        READ (1,*) NUMPHASE
        DO I = 1, NUMPHASE
          READ (1,*) NUML, (CHI, L=1,NUML)
          MAXLEG = MAX(NUML,MAXLEG)
        ENDDO
        MAXPGL = NUMPHASE*MAXLEG

      ELSE
C           Standard property file has everything variable
        MAXLEG = NLEG
        DO WHILE (.TRUE.)
          READ (1,*,END=190) IX, IY, IZ, TMP, EXT, ALB, NUML,
     .                         (CHI, L=1,NUML)
          MAXLEG = MAX(NUML,MAXLEG)
        ENDDO
190     CONTINUE
        MAXPGL = NPX*NPY*NPZ*MAXLEG
      ENDIF
      CLOSE (1)
      RETURN
      END



      SUBROUTINE READ_PROPERTIES (PROPFILE, NPX, NPY, NPZ, 
     .             MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM, 
     .             PROPTYPE, DELX, DELY, ZLEVELS, MAXASYM,
     .             TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)
C       Reads the medium properties from the file into the property arrays.
C     There are three formats: the standard format specifies the extinction,
C     scattering albedo, temperature, and Legendre phase function expansion
C     at every grid point; the extinction only format only specifies
C     extinction at every grid point, while the other quantities are given
C     in the header part; and the tabulated phase function format specifies
C     all properties except the phase function for each grid point, while
C     many phase functions are tabulated in the header.  If the format is 
C     extinction only or tabulated phase function then NUMPHASE is the 
C     number of phase functions, otherwise NUMPHASE=0.
C     If doing delta-M then NLEG is the minimum number of Legendre terms on
C     input and the actual maximum number of terms on output, except that
C     it may not exceed MAXLEG; otherwise NLEG is the number of Legendre
C     terms to be used (regardless of what is in property file).
C     See the overall program documentation for the three file formats.
      INTEGER NPX, NPY, NPZ, MAXLEG, NLEG, NUMPHASE
Cf2py intent(in, out) NPX, NPY, NPZ
Cf2py intent(in) MAXLEG
Cf2py intent(in, out) NLEG, NUMPHASE
      INTEGER MAXNZ, MAXPG, MAXPGL
Cf2py intent(in) MAXNZ, MAXPG, MAXPGL
      INTEGER IPHASEP(MAXPG)
Cf2py intent(out) IPHASEP
      LOGICAL DELTAM
Cf2py intent(in) DELTAM
      REAL    DELX, DELY, ZLEVELS(NPZ)
Cf2py intent(out) DELX, DELY
Cf2py intent(out) ZLEVELS
      REAL    TEMPP(MAXPG), EXTINCTP(MAXPG)
Cf2py intent(out) TEMPP, EXTINCTP
      REAL    ALBEDOP(MAXPG), LEGENP(MAXPGL), MAXASYM
Cf2py intent(out) ALBEDOP, LEGENP
Cf2py intent(out) MAXASYM
      CHARACTER PROPFILE*64
Cf2py intent(in) PROPFILE
      CHARACTER PROPTYPE*1
Cf2py intent(out) PROPTYPE
      INTEGER I, IX, IY, IZ, K, KL, L, IL, NUML, IPH, N
      REAL    TMP, EXT, ALB, CHI
 

C          Open the file, figure out the type, and get the grid size
      OPEN (UNIT=1, FILE=PROPFILE, STATUS='OLD')
      READ (1,'(A1)') PROPTYPE
      IF (PROPTYPE .NE. 'E' .AND. PROPTYPE .NE. 'T')  REWIND (1)
      READ (1,*) NPX, NPY, NPZ
      IF (NPZ .GT. MAXNZ)  STOP 'READ_PROPERTIES: MAXNZ exceeded'
      READ (1,*) DELX, DELY, (ZLEVELS(K), K=1,NPZ)
 
C           Property file type E is for extinction only format
      IF (PROPTYPE .EQ. 'E') THEN
        NUMPHASE = 1
        READ (1,*) (TMP, IZ=1,NPZ)
        READ (1,*) ALB, NUML
        IF (DELTAM) THEN
          NLEG = MAX(NLEG,NUML)
          NLEG = MIN(NLEG,MAXLEG)
        ENDIF
        IF (NPX*NPY*NPZ .GT. MAXPG)
     .    STOP 'READ_PROPERTIES: MAXPG exceeded'
        IF (NLEG .GT. MAXPGL)
     .    STOP 'READ_PROPERTIES: MAXPGL exceeded'
        REWIND (1)
        READ (1,*)
        READ (1,*) NPX, NPY, NPZ
        READ (1,*) DELX, DELY, (ZLEVELS(K), K=1,NPZ)
        READ (1,*) (TEMPP(IZ), IZ=1,NPZ)
        READ (1,*) ALB, NUML, (LEGENP(L), L=1,MIN(NUML,NLEG)), 
     .                        (CHI, L=NLEG+1,NUML)
        MAXASYM = LEGENP(1)/3.0
        DO L = NUML+1, NLEG
          LEGENP(L) = 0.0
        ENDDO
        IY = 1
        DO WHILE (.TRUE.)
          IF (NPY .EQ. 1) THEN
              READ (1,*,END=290) IX, IZ, EXT
          ELSE
              READ (1,*,END=290) IX, IY, IZ, EXT
          ENDIF
          IF (IX .GE. 1 .AND. IX .LE. NPX .AND.
     .        IY .GE. 1 .AND. IY .LE. NPY .AND.
     .        IZ .GE. 1 .AND. IZ .LE. NPZ) THEN
              K = IZ + NPZ*(IY-1) + NPZ*NPY*(IX-1)
              TEMPP(K) = TEMPP(IZ)
              EXTINCTP(K) = EXT
              ALBEDOP(K) = ALB
              IPHASEP(K) = 1
          ENDIF
        ENDDO


C           Property file type T is for tabulated phase function format
      ELSE IF (PROPTYPE .EQ. 'T') THEN
        READ (1,*) NUMPHASE
C             If delta-M then find the largest number of Legendre terms
        IF (DELTAM) THEN
          DO I = 1, NUMPHASE
            READ (1,*) NUML, (CHI, L=1,NUML)
            NLEG = MAX(NUML,NLEG)
          ENDDO
          IF (NLEG .GT. MAXLEG) THEN
            WRITE (*,'(A,I4,A,I4,A)')
     .        'Warning: truncating phase function from ',
     .         NLEG,' to ',MAXLEG,' Legendre terms.'
          ENDIF
          NLEG = MIN(NLEG,MAXLEG)
        ENDIF
        N = NPX*NPY*NPZ
        IF (N .GT. MAXPG)  STOP 'READ_PROPERTIES: MAXPG exceeded'
        IF (NUMPHASE*NLEG .GT. MAXPGL)
     .    STOP 'READ_PROPERTIES: MAXPGL exceeded'
        REWIND (1)
        READ (1,*)
        READ (1,*) NPX, NPY, NPZ
        READ (1,*) DELX, DELY, (ZLEVELS(K), K=1,NPZ)
        READ (1,*) NUMPHASE
        MAXASYM = -1.0
        DO I = 1, NUMPHASE
          IL = NLEG*(I-1)
          READ (1,*) NUML, (LEGENP(L+IL), L=1,MIN(NUML,NLEG)), 
     .                     (CHI, L=NLEG+1,NUML)
          MAXASYM = MAX(MAXASYM,LEGENP(IL+1)/3.0)
          DO L = NUML+1, NLEG
            LEGENP(L+IL) = 0.0
          ENDDO
        ENDDO

        DO WHILE (.TRUE.)
          READ (1,*,END=290) IX, IY, IZ, TMP, EXT, ALB, IPH
          IF (IX .GE. 1 .AND. IX .LE. NPX .AND.
     .        IY .GE. 1 .AND. IY .LE. NPY .AND.
     .        IZ .GE. 1 .AND. IZ .LE. NPZ) THEN
              K = IZ + NPZ*(IY-1) + NPZ*NPY*(IX-1)
              TEMPP(K) = TMP
              EXTINCTP(K) = EXT
              ALBEDOP(K) = ALB
              IPHASEP(K) = IPH
          ENDIF
        ENDDO

      ELSE

 
C           Standard property file has everything variable
C             If delta-M then find the largest number of Legendre terms
        NUMPHASE = 0
        IF (DELTAM) THEN
          DO WHILE (.TRUE.)
            READ (1,*,END=190) IX, IY, IZ, TMP, EXT, ALB, NUML,
     .                         (CHI, L=1,NUML)
            NLEG = MAX(NUML,NLEG)
          ENDDO
190       CONTINUE
          NLEG = MIN(NLEG,MAXLEG)
          REWIND (1)
          READ (1,*) NPX, NPY, NPZ
          READ (1,*) DELX, DELY, (ZLEVELS(K), K=1,NPZ)
        ENDIF
        N = NPX*NPY*NPZ
        IF (N .GT. MAXPG)  STOP 'READ_PROPERTIES: MAXPG exceeded'
        IF (N*NLEG .GT. MAXPGL)
     .    STOP 'READ_PROPERTIES: MAXPGL exceeded'
        MAXASYM = -1.0
        DO WHILE (.TRUE.)
          READ (1,*,END=290) IX, IY, IZ, TMP, EXT, ALB, NUML
          IF (IX .GE. 1 .AND. IX .LE. NPX .AND.
     .        IY .GE. 1 .AND. IY .LE. NPY .AND.
     .        IZ .GE. 1 .AND. IZ .LE. NPZ) THEN
              K = IZ + NPZ*(IY-1) + NPZ*NPY*(IX-1)
              KL = NLEG*(K-1)
              BACKSPACE(1)
              READ (1,*,END=290) IX, IY, IZ, TEMPP(K), EXTINCTP(K),
     .            ALBEDOP(K), NUML, (LEGENP(L+KL), L=1,MIN(NUML,NLEG)),
     .                              (CHI, L=NLEG+1,NUML)
              MAXASYM = MAX(MAXASYM,LEGENP(KL+1)/3.0)
              DO L = NUML+1, NLEG
                LEGENP(L+KL) = 0.0
              ENDDO
          ENDIF
        ENDDO

      ENDIF
 
290   CONTINUE
      CLOSE (1)
      RETURN
      END

 


      SUBROUTINE CHECK_PROPERTY_INPUT (NPX, NPY, NPZ, NLEG,  
     .             DELX, DELY, ZLEVELS, 
     .             TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)
C       Checks the input medium property values for correctness.
C       Warning if temperature below 150 K or above 350 K.
C       Stop if negative extinction, single scattering albedo below 0 or
C       above 1, asymmetry parameter below -1 or above 1, or IPHASE
C       out of range (of tabulated phase functions).  
C       Warning if Chi1=1, for people specifying Chi0.
C       Warning if optical depth across a grid cell exceeds 2.
      INTEGER NPX, NPY, NPZ, NLEG
      INTEGER NUMPHASE
Cf2py intent(in) :: NPX, NPY, NPZ, NLEG, NUMPHASE
      INTEGER IPHASEP(NPX*NPY*NPZ)
Cf2py intent(in) :: IPHASEP
      REAL    DELX, DELY, ZLEVELS(NPZ)
Cf2py intent(in) :: DELX, DELY, ZLEVELS
      REAL    TEMPP(NPX*NPY*NPZ), EXTINCTP(NPX*NPY*NPZ)
Cf2py intent(in) :: TEMPP, EXTINCTP
      REAL    ALBEDOP(NPX*NPY*NPZ), LEGENP(NLEG,NUMPHASE)
Cf2py intent(in) :: ALBEDOP, LEGENP
      INTEGER I, IZ, N, NTHICK
      LOGICAL LEGWARN
      REAL    TMIN, TMAX, TAU, MAXTAU, DZ

C         NPX and NPY must be 1 or greater, NPZ must be 2 or greater
      IF (NPX .LT. 1 .OR. NPY .LT. 1 .OR. NPZ .LT. 2) THEN
        WRITE (*,*) 'CHECK_PROPERTY_INPUT: Bad grid size - ', 
     .         NPX, NPY, NPZ
        STOP
      ENDIF
      IF (DELX .LE. 0.0 .OR. DELY .LE. 0.0) THEN
        WRITE (*,*) 'CHECK_PROPERTY_INPUT: Bad grid spacing - ', 
     .         'DELX,DELY: ', DELX, DELY
        STOP
      ENDIF
C         The vertical levels must be increasing
      DO IZ = 1, NPZ-1
        IF (ZLEVELS(IZ+1) .LE. ZLEVELS(IZ)) THEN
          WRITE (*,*) 'CHECK_PROPERTY_INPUT: ZLEVELS must increase : ',
     .         IZ, ZLEVELS(IZ), ZLEVELS(IZ+1)
          STOP
        ENDIF
      ENDDO
C         Loop over the property grid points, checking parameters
C           Warn about temperatures outside of Earth atmosphere temperatures
      N = NPX*NPY*NPZ
      TMIN = 1000.0
      TMAX = -1000.0
      DO I = 1, N
        TMIN = MIN(TMIN,TEMPP(I))
        TMAX = MAX(TMAX,TEMPP(I))
      ENDDO
      IF (TMIN .LT. 150.0 .OR. TMAX .GT. 350.0) THEN
        WRITE (*,*) 'CHECK_PROPERTY_INPUT: Warning - ',
     .   'temperature min and max (Kelvin): ', I, TEMPP(I), TMIN, TMAX
        WRITE (*,*)
      ENDIF
      DO I = 1, N
        IF (EXTINCTP(I) .LT. 0.0) THEN
          WRITE (*,*) 'CHECK_PROPERTY_INPUT: Negative extinction - ',
     .        EXTINCTP(I), ' at grid point ', I
          STOP
        ENDIF
        IF (ALBEDOP(I) .LT. 0.0 .OR. ALBEDOP(I) .GT. 1.0) THEN
          WRITE (*,*) 'CHECK_PROPERTY_INPUT: Single scattering ',
     .      'albedo out of range - ', ALBEDOP(I), ' at grid point ', I
          STOP
        ENDIF
      ENDDO

      IF (NUMPHASE .GT. 0)  N = NUMPHASE
      LEGWARN = .FALSE.
      DO I = 1, N
        IF (LEGENP(1,I) .LT. -3.0 .OR. LEGENP(1,I) .GT. 3.0) THEN
          WRITE (*,*) 'CHECK_PROPERTY_INPUT: Asymmetry parameter ',
     .      'out of range - ', LEGENP(1,I)/3.0, ' at grid point ', I
          STOP
        ENDIF
        IF (LEGENP(1,I) .EQ. 1.0)  LEGWARN = .TRUE.
      ENDDO
      IF (LEGWARN) THEN
        WRITE (*,*) 'CHECK_PROPERTY_INPUT: Warning - ',
     .        'an asymmetry parameter is 1/3 (Chi1=1); ',
     .        '  Chi0 should not be specified for Chi1.'
        WRITE (*,*)
      ENDIF

      IF (NUMPHASE .GT. 0) THEN
        N = NPX*NPY*NPZ
        DO I = 1, N
          IF (IPHASEP(I) .LT. 1 .OR. IPHASEP(I) .GT. NUMPHASE) THEN
            WRITE (*,*) 'CHECK_PROPERTY_INPUT: Phase function pointer',
     .      ' out of range - ', IPHASEP(I), ' at grid point ', I
            STOP
          ENDIF
        ENDDO
      ENDIF

      N = NPX*NPY*NPZ
      NTHICK = 0
      MAXTAU = 0.0
      DO IZ = 1, NPZ-1
        DZ = ZLEVELS(IZ+1) - ZLEVELS(IZ)
        DO I = 0, N-1, NPZ
          TAU = 0.5*(EXTINCTP(I+IZ)+EXTINCTP(I+IZ+1)) *DZ
          IF (TAU .GT. 2.0) THEN
            NTHICK = NTHICK + 1
            MAXTAU = MAX(MAXTAU,TAU)
          ENDIF
        ENDDO
      ENDDO
      IF (NTHICK .GT. 0) THEN
        WRITE (*,*) 'CHECK_PROPERTY_INPUT: Warning -'
        WRITE (*,'(A,A,I7)') '  Number of property grid cells',
     .     ' with optical depth greater than 2 : ', NTHICK
        WRITE (*,*) ' Max cell optical depth : ', MAXTAU
        WRITE (*,*)
      ENDIF

      RETURN
      END
      




 
      SUBROUTINE READ_SURFACE_SIZE (SFCFILE, MAXSFCPTS, MAXSFCPARS)
C       Gets the size of the surface property inputs (temperature is
C       always included).
      INTEGER MAXSFCPTS, MAXSFCPARS
Cf2py intent(out) :: MAXSFCPTS, MAXSFCPARS
      CHARACTER SFCFILE*64
Cf2py intent(in) :: SFCFILE
      INTEGER NXSFC, NYSFC
      CHARACTER SFCTYPE*2
 
C          Open the file, figure out the type, and get the grid size
      SFCTYPE(1:1)='V'
      OPEN (UNIT=1, FILE=SFCFILE, STATUS='OLD')
      READ (1,'(A1)') SFCTYPE(2:2)
      READ (1,*) NXSFC, NYSFC
      MAXSFCPTS = (NXSFC+1)*(NYSFC+1)

      IF (SFCTYPE .EQ. 'VL' .OR. SFCTYPE .EQ. 'Vl') THEN
C           Surface file type L is for variable Lambertian surface
        MAXSFCPARS = 2
      ELSE IF (SFCTYPE .EQ. 'VF') THEN
C           Surface file type F is for Fresnel surface
        MAXSFCPARS = 3
      ELSE IF (SFCTYPE .EQ. 'VR') THEN
C           Surface file type R is for RPV surface:
        MAXSFCPARS = 4
      ELSE IF (SFCTYPE .EQ. 'VO') THEN
C           Surface file type O is for Ocean BRDF surface:
        MAXSFCPARS = 3
      ELSE
        STOP 'READ_SURFACE_SIZE: Unknown BRDF type'
      ENDIF
      CLOSE (1)
      RETURN
      END



 
      SUBROUTINE READ_SURFACE (SFCFILE, MAXSFCPTS, MAXSFCPARS, 
     .                     SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC, 
     .                     NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO)
C       Reads the surface properties from the file into the surface arrays.
C     There are currently four surface types (Lambertian, Fresnel, 
C     Rahman et al, and Ocean), though it is simple to add more.  In all cases
C     the surface properties are specified on a regular, evenly spaced 
C     grid of size NXSFC by NYSFC and spacing DELXSFC by DELYSFC.  
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
C       F  Fresnel       Real, Imaginary part of index of refraction
C       R  RPV           rho0, k, Theta
C       O  Ocean         Wind Speed (m/s), Pigmentation (mg/m^3)
C
      INTEGER MAXSFCPTS, MAXSFCPARS, NXSFC, NYSFC, NSFCPAR
Cf2py intent(in) :: MAXSFCPTS, MAXSFCPARS
Cf2py intent(out) :: NXSFC, NYSFC, NSFCPAR
      REAL    DELXSFC, DELYSFC
Cf2py intent(out) :: DELXSFC, DELYSFC
      REAL    SFCPARMS(MAXSFCPARS*MAXSFCPTS), GNDTEMP, GNDALBEDO
Cf2py intent(out) :: SFCPARMS, GNDTEMP, GNDALBEDO
      CHARACTER SFCFILE*64, SFCTYPE*2
Cf2py intent(in) :: SFCFILE
Cf2py intent(out) :: SFCTYPE
      INTEGER N, I, I0, IX, IY, J
      REAL    ALB, TEMP, MRE, MIM, RHO0, KR, THETA, WS, PCL
 
C       Surface type is variable something
      SFCTYPE(1:1)='V'
C          Open the file, figure out the type, and get the grid size
      OPEN (UNIT=1, FILE=SFCFILE, STATUS='OLD')
      READ (1,'(A1)') SFCTYPE(2:2)
      READ (1,*) NXSFC, NYSFC, DELXSFC, DELYSFC
      IF ( (NXSFC+1)*(NYSFC+1) .GT. MAXSFCPTS ) THEN
        WRITE (*,*) 'READ_SURFACE: MAXSFCPTS exceeded'
        STOP
      ENDIF

      GNDTEMP = 0.0
      GNDALBEDO = 0.0
      N = 0

C           Surface file type L is for variable Lambertian surface
      IF (SFCTYPE .EQ. 'VL' .OR. SFCTYPE .EQ. 'Vl') THEN
        NSFCPAR = 2
        IF (NSFCPAR .GT. MAXSFCPARS ) THEN
          WRITE (*,*) 'READ_SURFACE: MAXSFCPARS exceeded'
          STOP
        ENDIF
        DO WHILE (.TRUE.)
          READ (1,*,END=290) IX, IY, TEMP, ALB
          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
            I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
            SFCPARMS(I+1) = TEMP
            SFCPARMS(I+2) = ALB
            IF (ALB .LT. 0.0 .OR. ALB .GT. 1.0)
     .        STOP 'READ_SURFACE: Illegal surface albedo'
            IF (TEMP .LT. 150. .OR. TEMP .GT. 350.)
     .        WRITE (*,*) 'READ_SURFACE: Warning -',
     .          ' surface temperature: ',IX,IY,TEMP
            GNDTEMP = GNDTEMP + TEMP
            GNDALBEDO = GNDALBEDO + ALB
            N = N + 1
          ENDIF
        ENDDO

C           Surface file type F is for Fresnel surface
      ELSE IF (SFCTYPE .EQ. 'VF') THEN
        NSFCPAR = 3
        IF (NSFCPAR .GT. MAXSFCPARS ) THEN
          WRITE (*,*) 'READ_SURFACE: MAXSFCPARS exceeded'
          STOP
        ENDIF
        DO WHILE (.TRUE.)
          READ (1,*,END=290) IX, IY, TEMP, MRE, MIM
          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
            I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
            SFCPARMS(I+1) = TEMP
            SFCPARMS(I+2) = MRE
            SFCPARMS(I+3) = MIM
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
          WRITE (*,*) 'READ_SURFACE: MAXSFCPARS exceeded'
          STOP
        ENDIF
        DO WHILE (.TRUE.)
          READ (1,*,END=290) IX, IY, TEMP, RHO0, KR, THETA
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

C           Surface file type O is for Ocean surface BRDF 
C             Norm Loeb's modification of 6S ocean reflectance module.
      ELSE IF (SFCTYPE .EQ. 'VO') THEN
        NSFCPAR = 3
        IF (NSFCPAR .GT. MAXSFCPARS ) THEN
          WRITE (*,*) 'READ_SURFACE: MAXSFCPARS exceeded'
          STOP
        ENDIF
        DO WHILE (.TRUE.)
          READ (1,*,END=290) IX, IY, TEMP, WS, PCL
          IF (IX .GE. 1 .AND. IX .LE. NXSFC .AND.
     .        IY .GE. 1 .AND. IY .LE. NYSFC) THEN
            I = NSFCPAR*(IX-1 + (NXSFC+1)*(IY-1))
            SFCPARMS(I+1) = TEMP
            SFCPARMS(I+2) = WS
            SFCPARMS(I+3) = PCL
            GNDTEMP = GNDTEMP + TEMP
            N = N + 1
          ENDIF
        ENDDO

      ELSE
        STOP 'READ_SURFACE: Unknown BRDF type'
      ENDIF


290   CONTINUE
      CLOSE (1)

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




      SUBROUTINE READ_CKD_SIZE (CKDFILE, WAVENO, NG, NZCKD)
C       Find the size of the k-distribution arrays for this wavenumber band.
      INTEGER NG, NZCKD
      REAL    WAVENO(2)
      CHARACTER CKDFILE*64
      INTEGER NB, IB, JB, KB, N
      REAL    WAVENUM1, WAVENUM2, SF

      OPEN (UNIT=1, FILE=CKDFILE, STATUS='OLD')
      READ (1,*)
      READ (1,*) NB
      READ (1,*)
C         Read band information until the right one is found
      JB = 0
      DO IB = 1, NB
        READ (1,*) KB, WAVENUM1, WAVENUM2, SF, N
        IF ( ABS(WAVENUM1-WAVENO(1)) .LT. 1.0 .AND. 
     .       ABS(WAVENUM2-WAVENO(2)) .LT. 1.0 ) THEN
          BACKSPACE (1)
          NG = N
          READ (1,*) JB, WAVENUM1, WAVENUM2, SF, N
        ENDIF
      ENDDO
      IF (JB .EQ. 0) THEN
        WRITE (*,*) ' Wavenumber range not found in CKD file : ',
     .              WAVENO(1), WAVENO(2), CKDFILE
        STOP
      ENDIF
      READ (1,*) NZCKD
      CLOSE (1)
      RETURN
      END



      SUBROUTINE READ_CKD (MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX, 
     .                     NG, DELG, NZCKD, ZCKD, KABS)
C       Reads the information appropriate for one band from a 
C     correlated k-distribution file.   The wavenumber range (cm^-1)
C     in the file must match the desired in WAVENO.  The solar flux
C     is returned in SOLFLUX, the number of "g"'s in NG, the weights
C     or delta g's in DELG, the number of levels in NZCKD, the
C     levels from the top down in ZCKD, and the absorption coefficients
C     for each g and level in KABS.  The array sizes MAXNG and MAXNZ
C     are checked.
      INTEGER MAXNG, MAXNZ, NG, NZCKD
      REAL    WAVENO(2), SOLFLUX, DELG(*), ZCKD(*), KABS(*)
      CHARACTER CKDFILE*64
      INTEGER NB, IB, JB, KB, IG, I, J, N
      REAL    WAVENUM1, WAVENUM2, SF
 
      OPEN (UNIT=1, FILE=CKDFILE, STATUS='OLD')
      READ (1,*)
      READ (1,*) NB
      READ (1,*)
C         Read band information until the right one is found
      JB = 0
      DO IB = 1, NB
        READ (1,*) KB, WAVENUM1, WAVENUM2, SF, N
        IF ( ABS(WAVENUM1-WAVENO(1)) .LT. 1.0 .AND. 
     .       ABS(WAVENUM2-WAVENO(2)) .LT. 1.0 ) THEN
          BACKSPACE (1)
          NG = N
          IF (NG .GT. MAXNG)  STOP 'READ_CKD: MAXNG exceeded'
          READ (1,*) JB, WAVENO(1), WAVENO(2), SOLFLUX,
     .               N, (DELG(IG), IG=1, NG)
        ENDIF
      ENDDO
      IF (JB .EQ. 0) THEN
        WRITE (*,*) ' Wavenumber range not found in CKD file : ',
     .              WAVENO(1), WAVENO(2), CKDFILE
        STOP
      ENDIF
C         Read the Z level heights
      READ (1,*) NZCKD
      IF (NZCKD .GT. MAXNZ)  STOP 'READ_CKD: MAXNZ exceeded'
      READ (1,*)
      READ (1,*)
      DO I = 1, NZCKD
        READ (1,*) ZCKD(I)
      ENDDO
C         Skip over the irrelevant bands and read the absorption coefficients
      READ (1,*)
      DO I = 1, (JB-1)*NZCKD
        READ (1,*)
      ENDDO
      DO I = 1, NZCKD
        READ (1,*) IB, J, (KABS(I+(IG-1)*NZCKD), IG = 1, NG)
        DO IG = 1, NG
          IF (KABS(I+(IG-1)*NZCKD) .LT. 0.0) THEN
            WRITE (*,'(1X,A,A,I3,A,I3,A,I2)')
     .        'Negative absorption in CKD file.  ',
     .        'Band: ',IB, '   Level: ',I, '   k: ',IG
          ENDIF
        ENDDO
      ENDDO
      CLOSE (1)
      RETURN
      END


 


      SUBROUTINE RESTORE_STATE (INSAVEFILE, NX,NY,NZ, ML,MM,NCS,NLM,
     .       INRADFLAG, NEWGRIDFLAG, XGRID, YGRID, ZGRID,
     .       NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR, 
     .       CELLFLAGS,  FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)
C       Reads the binary save file containing a previous state of
C     the SHDOM solution. This includes the cell tree structure,
C     the spherical harmonic expansion of the source function,
C     the low order terms of the radiance expansion, and the flux array.
C     If the file name is 'NONE' then the INRADFLAG is set to false, 
C     otherwise it is true.  The input file may have any spherical 
C     harmonic truncation smaller than the current (ML,MM,NLM).
C     The flag NEWGRIDFLAG is set to false if reading succeeds.
      INTEGER NX, NY, NZ, ML, MM, NCS, NLM, NPTS, NCELLS
      INTEGER SHPTR(*), RSHPTR(*)
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
      INTEGER*2 CELLFLAGS(*)
      LOGICAL INRADFLAG, NEWGRIDFLAG
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)
      REAL    GRIDPOS(3,*), FLUXES(2,*), SOURCE(*), RADIANCE(*)
      CHARACTER INSAVEFILE*64
      INTEGER NXI, NYI, NZI, MLI, MMI, NCSI, NLMI
      INTEGER I, J, NS, NR

      INRADFLAG = .FALSE.
      NEWGRIDFLAG = .TRUE. 
      IF (INSAVEFILE(1:4) .EQ. 'NONE')  RETURN

      OPEN (UNIT=3, FILE=INSAVEFILE, STATUS='OLD',
     .              FORM='UNFORMATTED', ERR=900)
      READ (3) NXI, NYI, NZI, MLI, MMI, NCSI, NLMI, NPTS, NCELLS
      IF (NXI .NE. NX .OR. NYI .NE. NY .OR. NZI .NE. NZ) THEN
        WRITE (6,*) 'RESTORE_STATE: wrong grid size:', NXI, NYI, NZI
        CLOSE (3)
        RETURN
      ENDIF
      IF (MLI .GT. ML .OR. MMI .GT. MM .OR. NCSI .NE. NCS) THEN
        WRITE (6,*) 'RESTORE_STATE: incompatible SH:', MLI, MMI, NCSI
        CLOSE (3)
        RETURN
      ENDIF
      READ (3) (XGRID(I), I=1, NX+1)
      READ (3) (YGRID(I), I=1, NY+1)
      READ (3) (ZGRID(I), I=1, NZ)
      READ (3) ((GRIDPTR(J,I), J=1,8), I=1,NCELLS)
      READ (3) ((NEIGHPTR(J,I), J=1,6), I=1,NCELLS)
      READ (3) ((TREEPTR(J,I), J=1,2), I=1,NCELLS)
      READ (3) (CELLFLAGS(I), I=1,NCELLS)

      READ (3) ((GRIDPOS(J,I), J=1,3), I=1,NPTS)
      READ (3) ((FLUXES(J,I), J=1,2), I=1,NPTS)

      READ (3) (SHPTR(I), I=1,NPTS+1)
      NS = SHPTR(NPTS+1)
      READ (3) (SOURCE(I), I=1,NS)

      NR = (NCS+2)*NPTS
      READ (3) (RADIANCE(I), I=1,NR)
      RSHPTR(1) = 0
      DO I = 1, NPTS
        RSHPTR(I+1) = RSHPTR(I) + NCS+2
      ENDDO

      CLOSE (3)
      INRADFLAG = .TRUE.
      NEWGRIDFLAG = .FALSE. 
      RETURN

900   CONTINUE
      WRITE (6,*) 'Input radiance vector file not found.'
      INRADFLAG = .FALSE.
      NEWGRIDFLAG = .TRUE.
      RETURN
      END




      SUBROUTINE SAVE_STATE (OUTSAVEFILE, NX,NY,NZ, ML,MM,NCS, NLM,
     .       WORK, XGRID, YGRID, ZGRID,
     .       NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR, 
     .       CELLFLAGS,  FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)
C       Write a binary save file containing the current state of
C     the SHDOM solution. This includes the cell tree structure,
C     the spherical harmonic expansion of the source function,
C     the low order terms of the radiance expansion, and the flux array.
C     If the file name is 'NONE' then no file is written.
      INTEGER NX, NY, NZ, ML, MM, NCS, NLM, NPTS, NCELLS
Cf2py intent(in) :: NX, NY, NZ, ML, MM, NCS, NLM, NCELLS
      INTEGER SHPTR(NPTS+1), RSHPTR(NPTS)
Cf2py intent(in) :: SHPTR, RSHPTR
      INTEGER GRIDPTR(8,NCELLS), NEIGHPTR(6,NCELLS), TREEPTR(2,NCELLS)
Cf2py intent(in) :: GRIDPTR, NEIGHPTR, TREEPTR
      INTEGER*2 CELLFLAGS(NCELLS)
Cf2py intent(in) :: CELLFLAGS
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ),  WORK(*)
Cf2py intent(in) :: XGRID, YGRID, ZGRID,  WORK
      REAL    GRIDPOS(3,NPTS), FLUXES(2,NPTS), SOURCE(*), RADIANCE(*)
Cf2py intent(in) :: GRIDPOS, FLUXES, SOURCE, RADIANCE
      CHARACTER OUTSAVEFILE*64
Cf2py intent(in) :: OUTSAVEFILE
      INTEGER I, J, K, IR, NS, NR
 
      IF (OUTSAVEFILE(1:4) .EQ. 'NONE')  RETURN
 
      OPEN (UNIT=3, FILE=OUTSAVEFILE, STATUS='UNKNOWN',
     .              FORM='UNFORMATTED')
      WRITE (3) NX, NY, NZ, ML, MM, NCS, NLM, NPTS, NCELLS
 
      WRITE (3) (XGRID(I), I=1, NX+1)
      WRITE (3) (YGRID(I), I=1, NY+1)
      WRITE (3) (ZGRID(I), I=1, NZ)
      WRITE (3) ((GRIDPTR(J,I), J=1,8), I=1,NCELLS)
      WRITE (3) ((NEIGHPTR(J,I), J=1,6), I=1,NCELLS)
      WRITE (3) ((TREEPTR(J,I), J=1,2), I=1,NCELLS)
      WRITE (3) (CELLFLAGS(I), I=1,NCELLS)

      WRITE (3) ((GRIDPOS(J,I), J=1,3), I=1,NPTS)
      WRITE (3) ((FLUXES(J,I), J=1,2), I=1,NPTS)

      WRITE (3) (SHPTR(I), I=1,NPTS+1)
      NS = SHPTR(NPTS+1)
      WRITE (3) (SOURCE(I), I=1,NS)
   
      K = 1
      DO I = 1, NPTS
        IR = RSHPTR(I)
        DO J = 1, NCS+2
          WORK(K) = RADIANCE(IR+J)
          K = K + 1
        ENDDO
      ENDDO
      NR = (NCS+2)*NPTS
      WRITE (3) (WORK(I), I=1,NR)

      CLOSE (3)
      RETURN
      END
 
 



      SUBROUTINE OUTPUT_CELL_SPLIT (OUTFILE, GRIDPTR,GRIDPOS, EXTINCT, 
     .                              SHPTR, SOURCE, NCELLS)
C       Outputs the final cell splitting criterion for each grid cell
C     and direction (X,Y,Z).  Calls CELL_SPLIT_TEST.
C     Mainly for debugging purposes.
      INTEGER NCELLS
      INTEGER SHPTR(*), GRIDPTR(8,*)
      REAL    GRIDPOS(3,*), EXTINCT(*), SOURCE(*)
      CHARACTER OUTFILE*72
      INTEGER IP1, IP2, ICELL, IDIR
      REAL    ADAPTCRIT(3), MAXADAPT, XC, YC, ZC, EXT

      OPEN (UNIT=9, FILE=OUTFILE, STATUS='UNKNOWN')
      WRITE (9,*) '!  ICELL   Xc     Yc     Zc  ',
     .            'Extinct SplitX  SplitY  SplitZ'

C         Loop over all the cells
      DO ICELL = 1, NCELLS
C           Get the adaptive grid criterion for the 3 directions
        CALL CELL_SPLIT_TEST (GRIDPTR, GRIDPOS, EXTINCT,
     .                        SHPTR, SOURCE, ICELL,
     .                        ADAPTCRIT, MAXADAPT, IDIR)
C           Compute the cell center location and mean extinction
        IP1 = GRIDPTR(1,ICELL)
        IP2 = GRIDPTR(8,ICELL)
        XC = 0.5*(GRIDPOS(1,IP1)+GRIDPOS(1,IP2))
        YC = 0.5*(GRIDPOS(2,IP1)+GRIDPOS(2,IP2))
        ZC = 0.5*(GRIDPOS(3,IP1)+GRIDPOS(3,IP2))
        EXT = (EXTINCT(GRIDPTR(1,ICELL)) + EXTINCT(GRIDPTR(2,ICELL))
     .       + EXTINCT(GRIDPTR(3,ICELL)) + EXTINCT(GRIDPTR(4,ICELL))
     .       + EXTINCT(GRIDPTR(5,ICELL)) + EXTINCT(GRIDPTR(6,ICELL))
     .       + EXTINCT(GRIDPTR(7,ICELL)) + EXTINCT(GRIDPTR(8,ICELL)))/8
        WRITE (9,'(1X,I5,3(1X,F6.3),1X,E11.4,3(1X,F7.5))') 
     .        ICELL, XC, YC, ZC, EXT,
     .        ADAPTCRIT(1), ADAPTCRIT(2), ADAPTCRIT(3) 
      ENDDO

      CLOSE (9)
      RETURN
      END




      SUBROUTINE VISUALIZE_CELLS (OUTFILE, YVAL, OUTTYPE, 
     .      IPFLAG, NX, XDOMAIN, ZDOMAIN, 
     .      NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS)
C       Visualizes the cell structure by outputting GLE code to draw
C     a diagram of the cells in the XZ plane.  Only the end cells are
C     output.  The outlines of the cells are drawn, and arrows between
C     the cell centers show the neighbor relationships.
C     OUTTYPE=0 for cell outlines only, OUTTYPE=1 for cells and arrows.
C     There is a special output for IP mode.
      INTEGER OUTTYPE, IPFLAG, NX, NCELLS
      INTEGER GRIDPTR(8,*), NEIGHPTR(6,*), TREEPTR(2,*)
      INTEGER*2 CELLFLAGS(*)
      REAL    YVAL, XDOMAIN, ZDOMAIN, GRIDPOS(3,*)
      CHARACTER OUTFILE*72
      INTEGER IC
      INTEGER J, IN
      LOGICAL BTEST
      REAL    XSIZE, ZSIZE, S, DELX
      REAL    X1,Z1, X2,Z2, X3,Z3, X4,Z4, XC,ZC
      REAL    DX,DZ, XN,ZN, XA,ZA

      S = MIN(24.0/XDOMAIN,18.0/ZDOMAIN)
      XSIZE = S*XDOMAIN
      ZSIZE = S*ZDOMAIN
      OPEN (UNIT=9, FILE=OUTFILE, STATUS='UNKNOWN')
      WRITE (9,'(A,2(1X,F6.2))') 'size ', XSIZE, ZSIZE
      WRITE (9,'(A)') 'set just cc'
      WRITE (9,'(A)') 'amove 0 0'
      WRITE (9,'(A)') 'begin origin'
      WRITE (9,'(A)') 'begin scale 1 1'
      WRITE (9,'(A,1X,F7.3)') '! Y value: ', YVAL
      WRITE (9,'(A)') 
      IF (BTEST(IPFLAG,0)) THEN
        DELX = 0.5*XDOMAIN/NX
        WRITE (9,'(A)') 'sub drawbox x1 y1 x3 y3 xc yc n'
        WRITE (9,'(A)') ' set color black'
        WRITE (9,'(A)') ' amove x1 y1'
        WRITE (9,'(A,2F7.3)') ' rline ', S*DELX, 0
        WRITE (9,'(A)') ' amove x3 y3'
        WRITE (9,'(A,2F7.3)') ' rline ', S*DELX, 0
        WRITE (9,'(A)') ' amove xc yc'
        WRITE (9,'(A)') ' set hei 0.3'
        WRITE (9,'(A)') ' write num$(n)'
        WRITE (9,'(A)') 'end sub'
      ELSE
        WRITE (9,'(A)') 'sub drawbox x1 y1 x2 y2 x3 y3 x4 y4  xc yc n'
        WRITE (9,'(A)') ' set color black'
        WRITE (9,'(A)') ' amove x1 y1'
        WRITE (9,'(A)') ' aline x2 y2'
        WRITE (9,'(A)') ' aline x3 y3'
        WRITE (9,'(A)') ' aline x4 y4'
        WRITE (9,'(A)') ' aline x1 y1'
        WRITE (9,'(A)') ' amove xc yc'
        WRITE (9,'(A)') ' set hei 0.3'
        WRITE (9,'(A)') ' write num$(n)'
        WRITE (9,'(A)') 'end sub'
      ENDIF
      WRITE (9,'(A)') 
      IF (OUTTYPE .EQ. 1) THEN
        WRITE (9,'(A)') 'sub drawarrow x1 y1 x2 y2 i'
        WRITE (9,'(A)') ' if (i<0) then'
        WRITE (9,'(A)') '   set color red'
        WRITE (9,'(A)') ' else'
        WRITE (9,'(A)') '   set color green'
        WRITE (9,'(A)') ' end if'
        WRITE (9,'(A)') ' amove x1 y1'
        WRITE (9,'(A)') ' set hei 0.2'
        WRITE (9,'(A)') ' aline x2 y2 arrow end'
        WRITE (9,'(A)') ' set color black'
        WRITE (9,'(A)') 'end sub'
        WRITE (9,'(A)') 
      ENDIF

      DO IC = 1, NCELLS
C           Output if it is an end cell and crosses the Y value
        IF (TREEPTR(2,IC) .EQ. 0 .AND.
     .      GRIDPOS(2,GRIDPTR(1,IC)) .LE. YVAL .AND. 
     .      GRIDPOS(2,GRIDPTR(3,IC)) .GE. YVAL) THEN 
C             Get the center of the cell
          XC = (GRIDPOS(1,GRIDPTR(1,IC))+GRIDPOS(1,GRIDPTR(8,IC)))/2
          ZC = (GRIDPOS(3,GRIDPTR(1,IC))+GRIDPOS(3,GRIDPTR(8,IC)))/2
C             Draw the outline in XZ (four lines and four grid points)
C               and the cell number at the center using the gle subroutine
          X1 = GRIDPOS(1,GRIDPTR(1,IC))
          Z1 = GRIDPOS(3,GRIDPTR(1,IC))
          X2 = GRIDPOS(1,GRIDPTR(2,IC))
          Z2 = GRIDPOS(3,GRIDPTR(2,IC))
          X3 = GRIDPOS(1,GRIDPTR(6,IC))
          Z3 = GRIDPOS(3,GRIDPTR(6,IC))
          X4 = GRIDPOS(1,GRIDPTR(5,IC))
          Z4 = GRIDPOS(3,GRIDPTR(5,IC))
          IF (BTEST(IPFLAG,0)) THEN
            WRITE (9,'(A,6(1X,F6.3),I8)') '@drawbox ', 
     .       S*X1,S*Z1, S*X3,S*Z3, S*XC,S*ZC, IC
          ELSE
            WRITE (9,'(A,10(1X,F6.3),I8)') '@drawbox ', 
     .       S*X1,S*Z1, S*X2,S*Z2, S*X3,S*Z3, S*X4,S*Z4, S*XC,S*ZC, IC
          ENDIF
          IF (OUTTYPE .EQ. 1) THEN
C             Draw the neighbor arrows
            DX = 0.1*(X2-X1)
            DZ = 0.1*(Z4-Z1)
            DO J = 1, 6
             IF (NEIGHPTR(J,IC) .NE. 0 .AND. J.NE.3 .AND. J.NE.4) THEN
              IN = ABS(NEIGHPTR(J,IC))   
              XN=(GRIDPOS(1,GRIDPTR(1,IN))+GRIDPOS(1,GRIDPTR(8,IN)))/2
              ZN=(GRIDPOS(3,GRIDPTR(1,IN))+GRIDPOS(3,GRIDPTR(8,IN)))/2
              XA = XC + 1.0*(XN-XC)
              ZA = ZC + 1.0*(ZN-ZC)
              IF (.NOT. (XN-XC .GT. 0.0 .AND. J .EQ. 1 .OR. 
     .                   XN-XC .LT. 0.0 .AND. J .EQ. 2) ) THEN
                WRITE (9,'(A,4(1X,F6.3),1X,I8)') '@drawarrow ', 
     .            S*XC, S*ZC, S*XA, S*ZA, NEIGHPTR(J,IC)
              ENDIF
             ENDIF
            ENDDO
          ENDIF

        ENDIF
      ENDDO
      WRITE (9,'(A)') 'end scale'
      WRITE (9,'(A)') 'end origin'

      CLOSE (9)
      RETURN
      END





      SUBROUTINE OUTPUT_RESULTS (NX, NY, NZ, NBPTS, NPTS, NCELLS,
     .             NSH, ML,MM,NLM, NLEG, NUMPHASE, NMU,NPHI,NANG, NG, 
     .             PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE,
     .             BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,
     .             SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .             SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, 
     .             SPLITACC, SHACC, SOLACC, MAXITER, TOTITER, 
     .             XGRID, YGRID, ZGRID, GRIDPOS,
     .             TREEPTR, GRIDPTR, CELLFLAGS,
     .             EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP,  
     .             FLUXES, DIRFLUX, FLUXDIV, 
     .             IRAD, RADOUT, NSHOUT, SHTERMS, SOURCE1OUT,
     .             OUTTYPE, OUTPARMS, OUTFILE)
C       Writes the desired type of output file from the output fields
C     with a header giving the input parameters.
C     There are five types (OUTTYPE) of output: 'R' - for radiance,  
C     'F' - for hemispheric flux, 'H' - for heating rate (net flux 
C     convergence), 'S' - for spherical harmonic terms, and 'M' for
C     medium properties.
C     See the overall program documentation for output parameters.

      INTEGER NX, NY, NZ, NBPTS, NPTS, NCELLS, NSH
      INTEGER ML, MM, NLM, NLEG, NMU, NPHI, NANG, NG
      INTEGER NUMPHASE
      INTEGER BCFLAG, IPFLAG
      INTEGER MAXITER, TOTITER, IRAD, NSHOUT
      INTEGER TREEPTR(2,NCELLS), GRIDPTR(8,NCELLS)
      INTEGER*2 CELLFLAGS(NCELLS)
      INTEGER IPHASE(NPTS)
      LOGICAL DELTAM
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD
      REAL    SOLACC, SPLITACC, SHACC
      REAL    WAVENO(2), WAVELEN
      REAL    XGRID(NX+1), YGRID(NY+1), ZGRID(NZ), GRIDPOS(3,NPTS)
      REAL    OUTPARMS(*)
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1, GRIDTYPE*1, OUTTYPE*1
      CHARACTER*64 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE
      CHARACTER*64 OUTFILE
      REAL    EXTINCT(NPTS),ALBEDO(NPTS),LEGEN(0:NLEG,NPTS), TEMP(NPTS)
      REAL    FLUXES(2,NPTS), DIRFLUX(NPTS), FLUXDIV(NPTS)
      REAL    RADOUT(*), SHTERMS(NSHOUT,NPTS), SOURCE1OUT(*)
 
      INTEGER IZ, I, K, N, NX2, NY2, JX, JY, ICELL
      INTEGER NANGOUT, NXOUT, NYOUT
      LOGICAL ALLGRID
      REAL    MUOUT, PHID, C, PI, SUM, SUM1, SUM2, SUM3
      REAL    XDOMAIN, YDOMAIN, STARTX, STARTY
      REAL    FUP, FDN, FDIR,   ASYM
      DOUBLE PRECISION X0, Y0, Z0
      CHARACTER*64 FORM
      CHARACTER*32 GRIDNAME, SOURCENAME, UNITSNAME, OUTNAME, SFCNAME
 

      GRIDNAME = 'EVEN (X,Y)  '
      IF (GRIDTYPE .EQ. 'P') THEN
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
      ELSE IF (GRIDTYPE .EQ. 'F') THEN
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
      ELSE
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
      ENDIF
 
      IF (SRCTYPE .EQ. 'S') THEN
          SOURCENAME = 'SOLAR'
      ELSE IF (SRCTYPE .EQ. 'T') THEN
          SOURCENAME = 'THERMAL'
      ELSE IF (SRCTYPE .EQ. 'B') THEN
          SOURCENAME = 'SOLAR/THERMAL'
      ENDIF

      IF (SFCTYPE .EQ. 'FL') THEN 
        SFCNAME = 'FIXED LAMBERTIAN'
      ELSE IF (SFCTYPE .EQ. 'VL') THEN 
        SFCNAME = 'VARIABLE LAMBERTIAN'
      ELSE IF (SFCTYPE .EQ. 'VF') THEN 
        SFCNAME = 'VARIABLE FRESNEL'
      ELSE IF (SFCTYPE .EQ. 'VR') THEN 
        SFCNAME = 'VARIABLE RPV'
      ELSE IF (SFCTYPE(1:1) .EQ. 'V') THEN 
        SFCNAME = 'VARIABLE OTHER'
      ENDIF
 
      IF (UNITS .EQ. 'T') THEN
          UNITSNAME = 'KELVIN'
      ELSE IF (UNITS .EQ. 'B') THEN
        IF (OUTTYPE .EQ. 'F') THEN
            UNITSNAME = 'WATTS/(M^2)'
        ELSE IF (OUTTYPE .EQ. 'H') THEN
            UNITSNAME = 'WATTS/(M^2 KM)'
        ELSE
            UNITSNAME = 'WATTS/(M^2 STER)'
        ENDIF
      ELSE
        IF (OUTTYPE .EQ. 'F' .OR. OUTTYPE .EQ. 'S') THEN
            UNITSNAME = 'WATTS/(M^2 MICRON)'
        ELSE IF (OUTTYPE .EQ. 'H') THEN
            UNITSNAME = 'WATTS/(M^2 MICRON KM)'
        ELSE
            UNITSNAME = 'WATTS/(M^2 MICRON STER)'
        ENDIF
      ENDIF
      IF (OUTTYPE .EQ. 'R')  OUTNAME = 'RADIANCE'
      IF (OUTTYPE .EQ. 'F')  OUTNAME = 'FLUX'
      IF (OUTTYPE .EQ. 'H')  OUTNAME = 'NET_FLUX_DIV'
      IF (OUTTYPE .EQ. 'S')  OUTNAME = 'SPHERICAL-HARMONIC'
      IF (OUTTYPE .EQ. 'J')  OUTNAME = 'SOURCE_FUNCTION'
      IF (OUTTYPE .EQ. 'M')  OUTNAME = 'MEDIUM-PROPERTIES'
 
 
      OPEN (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')
      WRITE (2,'(A,A)') '! Spherical Harmonic Discrete Ordinate',
     .                  ' Radiative Transfer Output'
      WRITE (2,'(2(A,I3),A,I5,2(A,I3),A,I5,A,I10)') 
     .    '!  L=',ML,'  M=',MM, '  NLM=',NLM, 
     .    '   NMU=',NMU, '  NPHI=',NPHI, '  NANG=', NANG, 
     .    '   NSH=', NSH
      WRITE (2,'(3(A,I4),2(A,I8))') '!  NX=',NX,'   NY=',NY,
     .    '   NZ=',NZ, '    NPTS=',NPTS, '   NCELLS=',NCELLS
      WRITE (2,'(A,A32)')     '!  PROPERTY_FILE=', PROPFILE
      WRITE (2,'(A,A32,A,I2)') '!  CORRELATED_K-DIST_FILE=', CKDFILE,
     .                        '   NUM_G=', NG
      WRITE (2,'(A,A32)')     '!  INPUT_SAVE_FILE=', INSAVEFILE
      WRITE (2,'(A,A32)')     '!  OUTPUT_SAVE_FILE=', OUTSAVEFILE
      IF (DELTAM) THEN
        WRITE (2,'(A,A14,A)')  '!  SOURCE_TYPE=', SOURCENAME,
     .                        '      DELTA-M METHOD'
      ELSE
        WRITE (2,'(A,A14)')    '!  SOURCE_TYPE=', SOURCENAME
      ENDIF
      WRITE (2,'(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,
     .                         '   INDEPENDENT_PIXEL=', IPFLAG
      WRITE (2,'(A,A20,A,I1)') '!  SURFACE_TYPE=', SFCNAME,
     .                         '   HORIZ_BOUNDARY_COND=', BCFLAG
      IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
        IF (SFCTYPE(1:1) .EQ. 'V') THEN
          WRITE (2,'(A,A32,A,E12.5)') '!  SURFACE_FILE=', SFCFILE,
     .                                 '  SKY_RAD=', SKYRAD
        ELSE
          IF (SRCTYPE .EQ. 'B') THEN
            WRITE (2,'(A,F9.7,A,F8.3,A,E12.5)') 
     .        '!  GROUND_ALBEDO=',GNDALBEDO, '  GROUND_TEMP=',GNDTEMP,
     .        '  SKY_RAD=',SKYRAD
          ELSE
            WRITE (2,'(A,F9.7,A,E12.5)') 
     .        '!  GROUND_ALBEDO=',GNDALBEDO, '  SKY_RAD=',SKYRAD
          ENDIF
        ENDIF
        WRITE (2,'(A,E13.6,A,F10.7,A,F8.3)')
     .         '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,
     .         '   SOLAR_AZ=', SOLARAZ*180.0/ACOS(-1.0)
      ELSE
        IF (SFCTYPE(1:1) .EQ. 'V') THEN
          WRITE (2,'(A,A32)') '!  SURFACE_FILE=', SFCFILE
        ELSE
          WRITE (2,'(A,F8.3,A,F9.7)') 
     .      '!  GROUND_TEMP=',GNDTEMP, '  GROUND_EMIS=',1.0-GNDALBEDO
        ENDIF
        WRITE (2,'(A,F8.3)') '!  SKY_TEMP=', SKYRAD
      ENDIF

      IF (UNITS .EQ. 'T') THEN
        WRITE (2,'(A,A24)')   '!  UNITS=', UNITSNAME
      ELSE IF (UNITS .EQ. 'B') THEN
        WRITE (2,'(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,
     .           '   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2)
      ELSE
        WRITE (2,'(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,
     .           '   WAVELENGTH=', WAVELEN
      ENDIF
      WRITE (2,'(2(A,E10.3))') '!  SPLITTING_ACCURACY=', SPLITACC,
     .                         '   SPHERICAL_HARMONIC_ACCURACY=',SHACC
      WRITE (2,'(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
      WRITE (2,'(2(A,I4))')    '!  MAXIMUM_ITERATIONS=', MAXITER,
     .                         '   NUMBER_ITERATIONS=', TOTITER
      WRITE (2,'(A,A20)')      '!  OUTPUT_TYPE=', OUTNAME
 
      PI = ACOS(-1.0) 


      IF (OUTTYPE .EQ. 'R') THEN
C             Radiance output
        XDOMAIN = XGRID(NX+1-IBITS(BCFLAG,0,1))-XGRID(1) +2*OUTPARMS(4)
        STARTX = -OUTPARMS(4)
        IF (OUTPARMS(2) .EQ. 0.0) THEN
          NXOUT = 1
        ELSE
          NXOUT = MAX(1,NINT(XDOMAIN/OUTPARMS(2)))
        ENDIF
        YDOMAIN = YGRID(NY+1-IBITS(BCFLAG,1,1))-YGRID(1) +2*OUTPARMS(5)
        STARTY = -OUTPARMS(5)
        IF (OUTPARMS(3) .EQ. 0.0) THEN
          NYOUT = 1
        ELSE
          NYOUT = MAX(1,NINT(YDOMAIN/OUTPARMS(3)))
        ENDIF
        Z0 = MIN( MAX(OUTPARMS(1),ZGRID(1)), ZGRID(NZ))
        NANGOUT = NINT(OUTPARMS(6))

        WRITE (2,'(A,F7.3,3(A,I4))') '!    RADIANCE AT Z=', Z0,
     .      '    NXO=',NXOUT, '    NYO=',NYOUT, '    NDIR=',NANGOUT
        IF (UNITS .EQ. 'T') THEN
          FORM = '(2(1X,F7.3),2X,F6.2)'
          WRITE (2,'(A)') '!   X       Y       RADIANCE'
        ELSE
          FORM = '(2(1X,F7.3),2X,E12.5)'
          WRITE (2,'(A)') '!   X       Y        RADIANCE'
        ENDIF

        DO K = 1, NANGOUT
          MUOUT = OUTPARMS(2*K+5)
          PHID = OUTPARMS(2*K+6)
          WRITE (2,'(A,1X,F8.5,1X,F6.2,2X,A)') 
     .          '! ', MUOUT, PHID, '<- (mu,phi)'
          Y0 = STARTY
          DO JY = 1, NYOUT
            X0 = STARTX
            DO JX = 1, NXOUT
              IRAD = IRAD + 1
              WRITE (2,FORM) X0, Y0, RADOUT(IRAD)
              X0 = X0 + OUTPARMS(2)
            ENDDO
            Y0 = Y0 + OUTPARMS(3)
          ENDDO
        ENDDO


      ELSE IF (OUTTYPE .EQ. 'F') THEN
C           Hemispheric flux output
C             Format 1 is flux at top and bottom of medium at grid points
        IF (NINT(OUTPARMS(1)) .LE. 1) THEN
          WRITE (2,'(A,F7.3)') '!    UPWELLING FLUX:     Z=',ZGRID(NZ)
          WRITE (2,'(A,F7.3)') '!    DOWNWELLING FLUX:   Z=',ZGRID(1)
          IF (UNITS .EQ. 'T') THEN
            FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
            WRITE (2,'(A)') '!   X       Y        UP     DOWN'
          ELSE
            IF (SRCTYPE .NE. 'T') THEN
              FORM = '(2(1X,F7.3),1X,3(2X,E12.5))'
              WRITE (2,'(A)')
     .     '!   X       Y           UP       DOWN_DIFFUSE   DOWN_DIRECT'
            ELSE
             FORM = '(2(1X,F7.3),1X,2(2X,E12.5))'
             WRITE (2,'(A)') '!   X       Y           UP           DOWN'
            ENDIF
          ENDIF
          C = 1.0
          IF (UNITS .EQ. 'T')  C = 1.0/PI
          DO I = 1, NBPTS, NZ
            IF (GRIDPOS(1,I) .LT. XGRID(NX+1) .AND.
     .          GRIDPOS(2,I) .LT. YGRID(NY+1)) THEN
              IF (SRCTYPE .NE. 'T') THEN
                WRITE (2,FORM) GRIDPOS(1,I), GRIDPOS(2,I),
     .            C*FLUXES(2,I+NZ-1), C*FLUXES(1,I), DIRFLUX(I)
              ELSE
                WRITE (2,FORM) GRIDPOS(1,I), GRIDPOS(2,I),
     .              C*FLUXES(2,I+NZ-1), C*FLUXES(1,I)
              ENDIF
            ENDIF
          ENDDO 

        ELSE IF (NINT(OUTPARMS(1)) .EQ. 2) THEN
C             Format 2 is flux at a given level at regular locations
          Z0 = MIN( MAX(OUTPARMS(2),ZGRID(1)), ZGRID(NZ))
          IF (OUTPARMS(3) .EQ. 0.0)  OUTPARMS(3) = 1.0
          NX2 = NX+1-IBITS(BCFLAG,0,1)
          NXOUT = MAX(1,NINT((XGRID(NX2)-XGRID(1))/OUTPARMS(3)))
          IF (OUTPARMS(4) .EQ. 0.0)  OUTPARMS(4) = 1.0
          NY2 = NY+1-IBITS(BCFLAG,1,1)
          NYOUT = MAX(1,NINT((YGRID(NY2)-YGRID(1))/OUTPARMS(4)))
          WRITE (2,'(A,F7.3,2(A,I4))')'!    UPWELLING FLUX:     Z=',Z0,
     .        '    NXO=',NXOUT, '    NYO=',NYOUT
          WRITE (2,'(A,F7.3)') '!    DOWNWELLING FLUX:   Z=',Z0
          IF (UNITS .EQ. 'T') THEN
            FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
            WRITE (2,'(A)') '!   X       Y        UP     DOWN'
          ELSE
            IF (SRCTYPE .NE. 'T') THEN
              FORM = '(2(1X,F7.3),1X,3(2X,E12.5))'
              WRITE (2,'(A)')
     .     '!   X       Y           UP       DOWN_DIFFUSE   DOWN_DIRECT'
            ELSE
             FORM = '(2(1X,F7.3),1X,2(2X,E12.5))'
             WRITE (2,'(A)')'!   X       Y           UP           DOWN'
            ENDIF
          ENDIF
          C = 1.0
          IF (UNITS .EQ. 'T')  C = 1.0/PI
          Y0 = YGRID(1)
          DO JY = 1, NYOUT
            X0 = XGRID(1)
            DO JX = 1, NXOUT
              CALL LOCATE_GRID_CELL (NX, NY, NZ, XGRID, YGRID, ZGRID, 
     .                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,
     .                  BCFLAG, IPFLAG, X0, Y0, Z0,  ICELL)
              CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                                X0, Y0, Z0, 2, 1, FLUXES, FDN)
              CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                                X0, Y0, Z0, 2, 2, FLUXES, FUP)
              CALL INTERPOLATE_FIELD (ICELL, GRIDPTR, GRIDPOS, 
     .                                X0, Y0, Z0, 1, 1, DIRFLUX, FDIR)
              IF (SRCTYPE .NE. 'T') THEN
                WRITE (2,FORM) X0, Y0, C*FUP, C*FDN, FDIR
              ELSE
                WRITE (2,FORM) X0, Y0, C*FUP, C*FDN
              ENDIF
              X0 = X0 + OUTPARMS(3)
            ENDDO 
            Y0 = Y0 + OUTPARMS(4)
          ENDDO 

        ELSE IF (NINT(OUTPARMS(1)) .EQ. 3) THEN
C             Format 3 is domain averaged vertical profile
          IF (UNITS .EQ. 'T') THEN
            FORM = '(1X,F7.3,1X,2(2X,F6.2))'
            WRITE (2,'(A)') '!   Z        UP     DOWN'
          ELSE
            IF (SRCTYPE .NE. 'T') THEN
              FORM = '(1X,F7.3,1X,3(2X,E12.5))'
              WRITE (2,'(A)')
     .            '!    Z          UP       DOWN_DIFFUSE   DOWN_DIRECT'
            ELSE
              FORM = '(1X,F7.3,1X,2(2X,E12.5))'
              WRITE (2,'(A)') '!    Z          UP           DOWN'
            ENDIF
          ENDIF
          C = 1.0
          IF (UNITS .EQ. 'T')  C = 1.0/PI
          DO IZ = 1, NZ
            SUM1 = 0.0
            SUM2 = 0.0
            SUM3 = 0.0
            N = 0
            DO I = IZ, NBPTS, NZ
              IF (GRIDPOS(1,I) .LT. XGRID(NX+1) .AND.
     .            GRIDPOS(2,I) .LT. YGRID(NY+1)) THEN
                SUM1 = SUM1 + FLUXES(1,I)
                SUM2 = SUM2 + FLUXES(2,I)
                SUM3 = SUM3 + DIRFLUX(I)
                N = N + 1
              ENDIF
            ENDDO
            IF (SRCTYPE .NE. 'T') THEN
              WRITE (2,FORM) ZGRID(IZ), C*SUM2/N, C*SUM1/N, SUM3/N
            ELSE
              WRITE (2,FORM) ZGRID(IZ), C*SUM2/N, C*SUM1/N
            ENDIF
          ENDDO

        ELSE IF (NINT(OUTPARMS(1)) .GE. 4) THEN
C             Format 4 is fluxes at every base grid point
C             Format 5 is fluxes at every grid point
          IF (NINT(OUTPARMS(1)) .EQ. 5) THEN
            N = NPTS
            ALLGRID = .TRUE.
          ELSE
            N = NBPTS
            ALLGRID = .FALSE.
          ENDIF
          IF (UNITS .EQ. 'T') THEN
            FORM = '(3(1X,F7.3),1X,2(2X,F6.2))'
            WRITE (2,'(A)') '!   X       Y       Z        UP     DOWN'
          ELSE
            IF (SRCTYPE .NE. 'T') THEN
              FORM = '(3(1X,F7.3),1X,3(2X,E12.5))'
              WRITE (2,'(A,A)') '!   X       Y       Z',
     .          '           UP       DOWN_DIFFUSE   DOWN_DIRECT'
            ELSE
              FORM = '(3(1X,F7.3),1X,2(2X,E12.5))'
              WRITE (2,'(A)') 
     .       '!   X       Y       Z           UP           DOWN'
            ENDIF
          ENDIF
          C = 1.0
          IF (UNITS .EQ. 'T')  C = 1.0/PI
          DO I = 1, N
            IF ((GRIDPOS(1,I) .LT. XGRID(NX+1) .AND.
     .           GRIDPOS(2,I) .LT. YGRID(NY+1)) .OR. ALLGRID) THEN
              IF (SRCTYPE .NE. 'T') THEN
                WRITE (2,FORM) GRIDPOS(1,I),GRIDPOS(2,I),GRIDPOS(3,I),
     .            C*FLUXES(2,I), C*FLUXES(1,I), DIRFLUX(I)
              ELSE
                WRITE (2,FORM) GRIDPOS(1,I),GRIDPOS(2,I),GRIDPOS(3,I),
     .            C*FLUXES(2,I), C*FLUXES(1,I)
              ENDIF
            ENDIF
          ENDDO
        ENDIF 


      ELSE IF (OUTTYPE .EQ. 'H') THEN
 
C             Heating output: net flux convergence 
        IF (NINT(OUTPARMS(1)) .EQ. 1) THEN
C             Format 1 is domain averaged vertical profile
          WRITE (2,'(A,A)') '!    Z       -DIV(Fnet)'
          DO IZ = 1, NZ
            SUM = 0.0
            N = 0
            DO I = IZ, NBPTS, NZ
              IF (GRIDPOS(1,I) .LT. XGRID(NX+1) .AND.
     .            GRIDPOS(2,I) .LT. YGRID(NY+1)) THEN
                SUM = SUM + FLUXDIV(I)
                N = N + 1
              ENDIF
            ENDDO
            WRITE (2,'(1X,F7.3,2X,E13.5)') ZGRID(IZ), -SUM/N
          ENDDO

        ELSE 
          IF (NINT(OUTPARMS(1)) .EQ. 3) THEN
            N = NPTS
            ALLGRID = .TRUE.
          ELSE
            N = NBPTS
            ALLGRID = .FALSE.
          ENDIF
C             Format 2 is flux convergence for every base grid point
C             Format 3 is flux convergence for every grid point
          WRITE (2,'(A)') '!    X       Y        Z       -DIV(Fnet)'
          DO I = 1, N
            IF ((GRIDPOS(1,I) .LT. XGRID(NX+1) .AND.
     .           GRIDPOS(2,I) .LT. YGRID(NY+1)) .OR. ALLGRID) THEN
              WRITE (2,'(3(1X,F7.3),2X,E13.5)') 
     .          GRIDPOS(1,I), GRIDPOS(2,I), GRIDPOS(3,I), -FLUXDIV(I)
            ENDIF
          ENDDO
        ENDIF
 
 
      ELSE IF (OUTTYPE .EQ. 'S') THEN
 
C          Spherical Harmonic output: 
C            Output mean intensity and net flux (x,y,z components).
C               and maybe normalized rms of higher order terms 
        IF (NINT(OUTPARMS(1)) .EQ. 2) THEN
          N = NPTS
          ALLGRID = .TRUE.
        ELSE
          N = NBPTS
          ALLGRID = .FALSE.
        ENDIF
        IF (NSHOUT .EQ. 5) THEN
          WRITE (2,'(A,A)') '!    X       Y        Z         Imean',
     .     '          Fx           Fy           Fz       HOrms'
          DO I = 1, N
            IF ((GRIDPOS(1,I) .LT. XGRID(NX+1) .AND.
     .         GRIDPOS(2,I) .LT. YGRID(NY+1)) .OR. ALLGRID) THEN
              WRITE (2,'(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')
     .          GRIDPOS(1,I), GRIDPOS(2,I), GRIDPOS(3,I), 
     .          SHTERMS(1,I), SHTERMS(2,I),SHTERMS(3,I),SHTERMS(4,I),
     .          SHTERMS(5,I)/MAX(1E-20,SHTERMS(1,I))
            ENDIF
          ENDDO
        ELSE
          WRITE (2,'(A,A)') '!    X       Y        Z         Imean',
     .     '          Fx           Fy           Fz '
          DO I = 1, N
            IF ((GRIDPOS(1,I) .LT. XGRID(NX+1) .AND.
     .         GRIDPOS(2,I) .LT. YGRID(NY+1)) .OR. ALLGRID) THEN
              WRITE (2,'(3(1X,F7.3),1X,4(1X,E12.5))')
     .          GRIDPOS(1,I), GRIDPOS(2,I), GRIDPOS(3,I), 
     .          SHTERMS(1,I), SHTERMS(2,I),SHTERMS(3,I),SHTERMS(4,I)
            ENDIF
          ENDDO
        ENDIF


      ELSE IF (OUTTYPE .EQ. 'J') THEN
 
C          J source function output: source function for a specified angle
C            at every grid point 
        IF (NINT(OUTPARMS(1)) .EQ. 2) THEN
          N = NPTS
          ALLGRID = .TRUE.
        ELSE
          N = NBPTS
          ALLGRID = .FALSE.
        ENDIF
        MUOUT = OUTPARMS(2)
        PHID = OUTPARMS(3)
        WRITE (2,'(A,1X,F8.5,1X,F6.2,2X,A)') 
     .          '! ', MUOUT, PHID, '<- (mu,phi)'
        WRITE (2,'(A,A)') '!    X       Y       Z    ',
     .                    ' Extinction   Source Function'
        DO I = 1, N
          IF ((GRIDPOS(1,I) .LT. XGRID(NX+1) .AND.
     .         GRIDPOS(2,I) .LT. YGRID(NY+1)) .OR. ALLGRID) THEN
            WRITE (2,'(3(1X,F7.3),2(1X,E12.5))')
     .           GRIDPOS(1,I), GRIDPOS(2,I), GRIDPOS(3,I), 
     .           EXTINCT(I), SOURCE1OUT(I)
          ENDIF
        ENDDO


      ELSE IF (OUTTYPE .EQ. 'M') THEN
 
C          Medium output: grid point medium properties -
C            Output extinction, single scattering albedo, 
C            asymmetry parameter, and temperature.
        IF (NINT(OUTPARMS(1)) .EQ. 2) THEN
          N = NPTS
          ALLGRID = .TRUE.
        ELSE
          N = NBPTS
          ALLGRID = .FALSE.
        ENDIF
        WRITE (2,'(A,A)') '!    X       Y       Z',
     .         '      Extinct  Albedo  Asymmetry Temp'
        DO I = 1, N
          IF (NUMPHASE .GT. 0) THEN
            ASYM = LEGEN(1,IPHASE(I))
          ELSE
            ASYM = LEGEN(1,I)
          ENDIF
          IF ((GRIDPOS(1,I) .LT. XGRID(NX+1) .AND.
     .         GRIDPOS(2,I) .LT. YGRID(NY+1)) .OR. ALLGRID) THEN
            WRITE (2,'(3(1X,F7.3),1X,F10.6,1X,F8.6,1X,F8.5,1X,F7.2)')
     .           GRIDPOS(1,I), GRIDPOS(2,I), GRIDPOS(3,I), 
     .           EXTINCT(I), ALBEDO(I), ASYM, TEMP(I)
          ENDIF
        ENDDO

      ENDIF
 
      CLOSE (2)
      RETURN
      END
 




      SUBROUTINE OUTPUT_RESULTS_PAR (NX, NY, NZ, NPTS, NCELLS,
     .             NSH, ML,MM,NLM, NMU,NPHI,NANG, NG, 
     .             PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE,
     .             BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,
     .             SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .             SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, 
     .             SPLITACC, SHACC, SOLACC, MAXITER, TOTITER, 
     .             XDOMAIN, YDOMAIN, XGRID, YGRID, ZGRID, 
     .             FLUXES, FLUXDIV, NSHOUT, SHTERMS, IRAD, RADOUT, 
     .             OUTTYPE, OUTPARMS, OUTFILE)
C       Writes the desired type of output file from the output fields
C     with a header giving the input parameters.
C     There are four supported types (OUTTYPE) of output: 'R' - for radiance,  
C     'F' - for hemispheric flux, 'H' - for heating rate (net flux 
C     convergence), 'S' - for spherical harmonic terms
C     See the overall program documentation for output parameters.

      INTEGER NX, NY, NZ, NPTS, NCELLS, NSH
      INTEGER ML, MM, NLM, NMU, NPHI, NANG, NG
      INTEGER BCFLAG, IPFLAG
      INTEGER MAXITER, TOTITER, IRAD, NSHOUT
      LOGICAL DELTAM
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD
      REAL    SOLACC, SPLITACC, SHACC
      REAL    WAVENO(2), WAVELEN
      REAL    XDOMAIN, YDOMAIN, XGRID(NX), YGRID(NY), ZGRID(NZ)
      REAL    OUTPARMS(*)
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1, GRIDTYPE*1, OUTTYPE*1
      CHARACTER*64 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE
      CHARACTER*64 OUTFILE
      REAL    FLUXES(3,NZ,NY,NX), FLUXDIV(NZ,NY,NX)
      REAL    RADOUT(*), SHTERMS(NSHOUT,NZ,NY,NX)
 
      INTEGER I, J, K, N, JX, JY
      INTEGER NANGOUT, NXOUT, NYOUT
      REAL    MUOUT, PHID, C, PI, SUM, SUM1, SUM2, SUM3
      REAL    STARTX, STARTY
      DOUBLE PRECISION X0, Y0, Z0
      CHARACTER*64 FORM
      CHARACTER*32 GRIDNAME, SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

      GRIDNAME = 'EVEN (X,Y)  '
      IF (GRIDTYPE .EQ. 'P') THEN
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
      ELSE IF (GRIDTYPE .EQ. 'F') THEN
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
      ELSE
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
      ENDIF
 
      IF (SRCTYPE .EQ. 'S') THEN
          SOURCENAME = 'SOLAR'
      ELSE IF (SRCTYPE .EQ. 'T') THEN
          SOURCENAME = 'THERMAL'
      ELSE IF (SRCTYPE .EQ. 'B') THEN
          SOURCENAME = 'SOLAR/THERMAL'
      ENDIF

      IF (SFCTYPE .EQ. 'FL') THEN 
        SFCNAME = 'FIXED LAMBERTIAN'
      ELSE IF (SFCTYPE .EQ. 'VL') THEN 
        SFCNAME = 'VARIABLE LAMBERTIAN'
      ELSE IF (SFCTYPE .EQ. 'VF') THEN 
        SFCNAME = 'VARIABLE FRESNEL'
      ELSE IF (SFCTYPE .EQ. 'VR') THEN 
        SFCNAME = 'VARIABLE RPV'
      ELSE IF (SFCTYPE(1:1) .EQ. 'V') THEN 
        SFCNAME = 'VARIABLE OTHER'
      ENDIF
 
      IF (UNITS .EQ. 'T') THEN
          UNITSNAME = 'KELVIN'
      ELSE IF (UNITS .EQ. 'B') THEN
        IF (OUTTYPE .EQ. 'F') THEN
            UNITSNAME = 'WATTS/(M^2)'
        ELSE IF (OUTTYPE .EQ. 'H') THEN
            UNITSNAME = 'WATTS/(M^2 KM)'
        ELSE
            UNITSNAME = 'WATTS/(M^2 STER)'
        ENDIF
      ELSE
        IF (OUTTYPE .EQ. 'F' .OR. OUTTYPE .EQ. 'S') THEN
            UNITSNAME = 'WATTS/(M^2 MICRON)'
        ELSE IF (OUTTYPE .EQ. 'H') THEN
            UNITSNAME = 'WATTS/(M^2 MICRON KM)'
        ELSE
            UNITSNAME = 'WATTS/(M^2 MICRON STER)'
        ENDIF
      ENDIF
      IF (OUTTYPE .EQ. 'R')  OUTNAME = 'RADIANCE'
      IF (OUTTYPE .EQ. 'F')  OUTNAME = 'FLUX'
      IF (OUTTYPE .EQ. 'H')  OUTNAME = 'NET_FLUX_DIV'
      IF (OUTTYPE .EQ. 'S')  OUTNAME = 'SPHERICAL-HARMONIC'
 
 
      OPEN (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')
      WRITE (2,'(A,A)') '! Spherical Harmonic Discrete Ordinate',
     .                  ' Radiative Transfer Output'
      WRITE (2,'(2(A,I3),A,I5,2(A,I3),A,I5,A,I10)') 
     .    '!  L=',ML,'  M=',MM, '  NLM=',NLM, 
     .    '   NMU=',NMU, '  NPHI=',NPHI, '  NANG=', NANG, 
     .    '   NSH=', NSH
      WRITE (2,'(3(A,I4),2(A,I8))') '!  NX=',NX,'   NY=',NY,
     .    '   NZ=',NZ, '    NPTS=',NPTS, '   NCELLS=',NCELLS
      WRITE (2,'(A,A32)')     '!  PROPERTY_FILE=', PROPFILE
      WRITE (2,'(A,A32,A,I2)') '!  CORRELATED_K-DIST_FILE=', CKDFILE,
     .                        '   NUM_G=', NG
      WRITE (2,'(A,A32)')     '!  INPUT_SAVE_FILE=', INSAVEFILE
      WRITE (2,'(A,A32)')     '!  OUTPUT_SAVE_FILE=', OUTSAVEFILE
      IF (DELTAM) THEN
        WRITE (2,'(A,A14,A)')  '!  SOURCE_TYPE=', SOURCENAME,
     .                        '      DELTA-M METHOD'
      ELSE
        WRITE (2,'(A,A14)')    '!  SOURCE_TYPE=', SOURCENAME
      ENDIF
      WRITE (2,'(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,
     .                         '   INDEPENDENT_PIXEL=', IPFLAG
      WRITE (2,'(A,A20,A,I2)') '!  SURFACE_TYPE=', SFCNAME,
     .                         '   HORIZ_BOUNDARY_COND=', BCFLAG
      IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
        IF (SFCTYPE(1:1) .EQ. 'V') THEN
          WRITE (2,'(A,A32,A,E12.5)') '!  SURFACE_FILE=', SFCFILE,
     .                                 '  SKY_RAD=', SKYRAD
        ELSE
          IF (SRCTYPE .EQ. 'B') THEN
            WRITE (2,'(A,F9.7,A,F8.3,A,E12.5)') 
     .        '!  GROUND_ALBEDO=',GNDALBEDO, '  GROUND_TEMP=',GNDTEMP,
     .        '  SKY_RAD=',SKYRAD
          ELSE
            WRITE (2,'(A,F9.7,A,E12.5)') 
     .        '!  GROUND_ALBEDO=',GNDALBEDO, '  SKY_RAD=',SKYRAD
          ENDIF
        ENDIF
        WRITE (2,'(A,E13.6,A,F10.7,A,F8.3)')
     .         '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,
     .         '   SOLAR_AZ=', SOLARAZ*180.0/ACOS(-1.0)
      ELSE
        IF (SFCTYPE(1:1) .EQ. 'V') THEN
          WRITE (2,'(A,A32)') '!  SURFACE_FILE=', SFCFILE
        ELSE
          WRITE (2,'(A,F8.3,A,F9.7)') 
     .      '!  GROUND_TEMP=',GNDTEMP, '  GROUND_EMIS=',1.0-GNDALBEDO
        ENDIF
        WRITE (2,'(A,F8.3)') '!  SKY_TEMP=', SKYRAD
      ENDIF

      IF (UNITS .EQ. 'T') THEN
        WRITE (2,'(A,A24)')   '!  UNITS=', UNITSNAME
      ELSE IF (UNITS .EQ. 'B') THEN
        WRITE (2,'(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,
     .           '   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2)
      ELSE
        WRITE (2,'(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,
     .           '   WAVELENGTH=', WAVELEN
      ENDIF
      WRITE (2,'(2(A,E10.3))') '!  SPLITTING_ACCURACY=', SPLITACC,
     .                         '   SPHERICAL_HARMONIC_ACCURACY=',SHACC
      WRITE (2,'(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
      WRITE (2,'(2(A,I4))')    '!  MAXIMUM_ITERATIONS=', MAXITER,
     .                         '   NUMBER_ITERATIONS=', TOTITER
      WRITE (2,'(A,A20)')      '!  OUTPUT_TYPE=', OUTNAME
 
      PI = ACOS(-1.0) 


      IF (OUTTYPE .EQ. 'R') THEN
C             Radiance output
        STARTX = 0.0
        IF (OUTPARMS(2) .LE. 0.0) THEN
          NXOUT = 1
        ELSE
          NXOUT = MAX(1,NINT(XDOMAIN/OUTPARMS(2)))
        ENDIF
        STARTY = 0.0
        IF (OUTPARMS(3) .LE. 0.0) THEN
          NYOUT = 1
        ELSE
          NYOUT = MAX(1,NINT(YDOMAIN/OUTPARMS(3)))
        ENDIF
        Z0 = MIN( MAX(OUTPARMS(1),ZGRID(1)), ZGRID(NZ))
        NANGOUT = NINT(OUTPARMS(6))

        WRITE (2,'(A,F7.3,3(A,I4))') '!    RADIANCE AT Z=', Z0,
     .      '    NXO=',NXOUT, '    NYO=',NYOUT, '    NDIR=',NANGOUT
        IF (UNITS .EQ. 'T') THEN
          FORM = '(2(1X,F7.3),2X,F6.2)'
          WRITE (2,'(A)') '!   X       Y      RADIANCE'
        ELSE
          FORM = '(2(1X,F7.3),2X,E12.5)'
          WRITE (2,'(A)') '!   X       Y        RADIANCE'
        ENDIF

        DO K = 1, NANGOUT
          MUOUT = OUTPARMS(2*K+5)
          PHID = OUTPARMS(2*K+6)
          WRITE (2,'(A,1X,F8.5,1X,F6.2,2X,A)') 
     .          '! ', MUOUT, PHID, '<- (mu,phi)'
          Y0 = STARTY
          DO JY = 1, NYOUT
            X0 = STARTX
            DO JX = 1, NXOUT
              IRAD = IRAD + 1
              WRITE (2,FORM) X0, Y0, RADOUT(IRAD)
              X0 = X0 + OUTPARMS(2)
            ENDDO
            Y0 = Y0 + OUTPARMS(3)
          ENDDO
        ENDDO


      ELSE IF (OUTTYPE .EQ. 'F') THEN
C           Hemispheric flux output
C             Format 1 is flux at top and bottom of medium at grid points
        IF (NINT(OUTPARMS(1)) .LE. 1) THEN
          WRITE (2,'(A,F7.3)') '!    UPWELLING FLUX:     Z=',ZGRID(NZ)
          WRITE (2,'(A,F7.3)') '!    DOWNWELLING FLUX:   Z=',ZGRID(1)
          IF (UNITS .EQ. 'T') THEN
            FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
            WRITE (2,'(A)') '!   X       Y        UP     DOWN'
          ELSE
            IF (SRCTYPE .NE. 'T') THEN
              FORM = '(2(1X,F7.3),1X,3(2X,E12.5))'
              WRITE (2,'(A)')
     .     '!   X       Y           UP       DOWN_DIFFUSE   DOWN_DIRECT'
            ELSE
             FORM = '(2(1X,F7.3),1X,2(2X,E12.5))'
             WRITE (2,'(A)') '!   X       Y           UP           DOWN'
            ENDIF
          ENDIF
          C = 1.0
          IF (UNITS .EQ. 'T')  C = 1.0/PI
          DO I = 1, NX
            DO J = 1, NY
              IF (SRCTYPE .NE. 'T') THEN
                WRITE (2,FORM) XGRID(I), YGRID(J),
     .            C*FLUXES(2,NZ,J,I), C*FLUXES(1,1,J,I), FLUXES(3,1,J,I)
              ELSE
                WRITE (2,FORM) XGRID(I), YGRID(J),
     .              C*FLUXES(2,NZ,J,I), C*FLUXES(1,1,J,I)
              ENDIF
            ENDDO 
          ENDDO 

        ELSE IF (NINT(OUTPARMS(1)) .EQ. 2) THEN
C             Format 2 is flux at a given level at regular locations
C             Not supported for multiple processors          

        ELSE IF (NINT(OUTPARMS(1)) .EQ. 3) THEN
C             Format 3 is domain averaged vertical profile
          IF (UNITS .EQ. 'T') THEN
            FORM = '(1X,F7.3,1X,2(2X,F6.2))'
            WRITE (2,'(A)') '!   Z        UP     DOWN'
          ELSE
            IF (SRCTYPE .NE. 'T') THEN
              FORM = '(1X,F7.3,1X,3(2X,E12.5))'
              WRITE (2,'(A)')
     .            '!    Z          UP       DOWN_DIFFUSE   DOWN_DIRECT'
            ELSE
              FORM = '(1X,F7.3,1X,2(2X,E12.5))'
              WRITE (2,'(A)') '!    Z          UP           DOWN'
            ENDIF
          ENDIF
          C = 1.0
          IF (UNITS .EQ. 'T')  C = 1.0/PI
          N = NX*NY
          DO K = 1, NZ
            SUM1 = 0.0
            SUM2 = 0.0
            SUM3 = 0.0
            DO I = 1, NX
              DO J = 1, NY
                SUM1 = SUM1 + FLUXES(1,K,J,I)
                SUM2 = SUM2 + FLUXES(2,K,J,I)
                SUM3 = SUM3 + FLUXES(3,K,J,I)
              ENDDO
            ENDDO
            IF (SRCTYPE .NE. 'T') THEN
              WRITE (2,FORM) ZGRID(K), C*SUM2/N, C*SUM1/N, SUM3/N
            ELSE
              WRITE (2,FORM) ZGRID(K), C*SUM2/N, C*SUM1/N
            ENDIF
          ENDDO

        ELSE
C             Format 4 is fluxes at every base grid point
C             Format 5 is fluxes at every grid point (not supported)
          IF (UNITS .EQ. 'T') THEN
            FORM = '(3(1X,F7.3),1X,2(2X,F6.2))'
            WRITE (2,'(A)') '!   X       Y       Z        UP     DOWN'
          ELSE
            IF (SRCTYPE .NE. 'T') THEN
              FORM = '(3(1X,F7.3),1X,3(2X,E12.5))'
              WRITE (2,'(A,A)') '!   X       Y       Z',
     .          '           UP       DOWN_DIFFUSE   DOWN_DIRECT'
            ELSE
              FORM = '(3(1X,F7.3),1X,2(2X,E12.5))'
              WRITE (2,'(A)') 
     .       '!   X       Y       Z           UP           DOWN'
            ENDIF
          ENDIF
          C = 1.0
          IF (UNITS .EQ. 'T')  C = 1.0/PI
          DO I = 1, NX
           DO J = 1, NY
            DO K = 1, NZ
              IF (SRCTYPE .NE. 'T') THEN
                WRITE (2,FORM) XGRID(I),YGRID(J),ZGRID(K),
     .            C*FLUXES(2,K,J,I), C*FLUXES(1,K,J,I), FLUXES(3,K,J,I)
              ELSE
                WRITE (2,FORM) XGRID(I),YGRID(J),ZGRID(K),
     .            C*FLUXES(2,K,J,I), C*FLUXES(1,K,J,I)
              ENDIF
            ENDDO
           ENDDO
          ENDDO
        ENDIF 


      ELSE IF (OUTTYPE .EQ. 'H') THEN
 
C             Heating output: net flux convergence 
        IF (NINT(OUTPARMS(1)) .EQ. 1) THEN
C             Format 1 is domain averaged vertical profile
          WRITE (2,'(A,A)') '!    Z       -DIV(Fnet)'
          DO K = 1, NZ
            SUM = 0.0
            DO I = 1, NX
              DO J = 1, NY
                SUM = SUM + FLUXDIV(K,J,I)
              ENDDO
            ENDDO
            WRITE (2,'(1X,F7.3,2X,E13.5)') ZGRID(K), -SUM/(NX*NY)
          ENDDO

        ELSE 
C             Format 2 is flux convergence for every base grid point
C             Format 3 is flux convergence for every grid point (not supported)
          WRITE (2,'(A)') '!    X       Y        Z       -DIV(Fnet)'
          DO I = 1, NX
           DO J = 1, NY
            DO K = 1, NZ
              WRITE (2,'(3(1X,F7.3),2X,E13.5)') 
     .             XGRID(I), YGRID(J), ZGRID(K), -FLUXDIV(K,J,I)
            ENDDO
           ENDDO
          ENDDO
        ENDIF
 
 
      ELSE IF (OUTTYPE .EQ. 'S') THEN
 
C          Spherical Harmonic output: 
C            Output mean intensity and net flux (x,y,z components).
C               and maybe normalized rms of higher order terms 
        IF (NSHOUT .EQ. 5) THEN
          WRITE (2,'(A,A)') '!    X       Y        Z         Imean',
     .     '          Fx           Fy           Fz       HOrms'
          DO I = 1, NX
           DO J = 1, NY
            DO K = 1, NZ
              WRITE (2,'(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')
     .          XGRID(I), YGRID(J), ZGRID(K), SHTERMS(1,K,J,I),
     .          SHTERMS(2,K,J,I),SHTERMS(3,K,J,I),SHTERMS(4,K,J,I),
     .          SHTERMS(5,K,J,I)/MAX(1E-20,SHTERMS(1,K,J,I))
            ENDDO
           ENDDO
          ENDDO
        ELSE
          WRITE (2,'(A,A)') '!    X       Y        Z         Imean',
     .     '          Fx           Fy           Fz '
          DO I = 1, NX
           DO J = 1, NY
            DO K = 1, NZ
              WRITE (2,'(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')
     .          XGRID(I), YGRID(J), ZGRID(K), SHTERMS(1,K,J,I),
     .          SHTERMS(2,K,J,I),SHTERMS(3,K,J,I),SHTERMS(4,K,J,I)
            ENDDO
           ENDDO
          ENDDO
        ENDIF

      ELSE IF (OUTTYPE .EQ. 'J') THEN
 
C          J source function output: source function for a specified angle
C            at every grid point (not supported)

      ELSE IF (OUTTYPE .EQ. 'M') THEN
 
C          Medium output: grid point medium properties -
C            Output extinction, single scattering albedo, 
C            asymmetry parameter, and temperature.  (not supported)

      ENDIF
 
      CLOSE (2)
      RETURN
      END
 






      SUBROUTINE OUTPUT_IMAGE (NX, NY, NZ, NPTS, NCELLS,
     .             NSH, ML,MM,NLM, NMU, NPHI, NANG, NG, 
     .             PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE,
     .             BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,
     .             SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,
     .             SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, 
     .             SPLITACC, SHACC, SOLACC, MAXITER, TOTITER, 
     .             IVIS, VISOUT, OUTPARMS, OUTFILE)
C       Writes the visualization image in PDS format: an ascii header 
C     giving the image parameters and the SHDOM input parameters,
C     followed by the byte or integer*2 image pixels.
      IMPLICIT NONE
      INTEGER NX, NY, NZ, NPTS, NCELLS, NSH
      INTEGER ML, MM, NLM, NMU, NPHI, NANG, NG
      INTEGER BCFLAG, IPFLAG
      INTEGER MAXITER, TOTITER, IVIS
      LOGICAL DELTAM
      REAL    SOLARFLUX, SOLARMU, SOLARAZ
      REAL    GNDTEMP, GNDALBEDO, SKYRAD
      REAL    SOLACC, SPLITACC, SHACC
      REAL    WAVENO(2), WAVELEN
      REAL    OUTPARMS(*)
      CHARACTER SRCTYPE*1, SFCTYPE*2, UNITS*1, GRIDTYPE*1
      CHARACTER*64 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE,OUTSAVEFILE
      CHARACTER*64 OUTFILE
      REAL    VISOUT(*)

      INTEGER MAXNS
      PARAMETER (MAXNS=2048) 
      INTEGER   MAXHEAD, NL, NS, LABRECS, RECLEN, NBYTES, FILESIZE
      INTEGER   DN, I, J, K
      LOGICAL   CAMERA_MODE
      INTEGER*1 BUFFER1(MAXNS)
      INTEGER*2 BUFFER2(MAXNS)
      CHARACTER*32 GRIDNAME, SOURCENAME, UNITSNAME, SFCNAME
      CHARACTER*1 CRLF(2)
      CHARACTER   BUF*2000
 
      MAXHEAD = 2000
      DO I = 1, MAXHEAD
        BUF(I:I) = ' '
      ENDDO

C         Find the image size
      CAMERA_MODE = NINT(OUTPARMS(1)) .EQ. 1
      IF (CAMERA_MODE) THEN
        NL = NINT(OUTPARMS(10))
        NS = NINT(OUTPARMS(11))
      ELSE IF (OUTPARMS(1) .EQ. 2) THEN
        NL = 1 + SQRT((OUTPARMS(4)-OUTPARMS(7))**2
     .               +(OUTPARMS(5)-OUTPARMS(8))**2
     .               +(OUTPARMS(6)-OUTPARMS(9))**2) /OUTPARMS(10)
        NS = 1 + ABS(OUTPARMS(12)-OUTPARMS(11))/OUTPARMS(13)
      ENDIF
      NBYTES = NINT(OUTPARMS(2))
      IF (NS .GT. MAXNS) STOP 'MAXNS exceeded in OUTPUT_IMAGE'

C       Open the PDF image
      RECLEN = NS*NBYTES
      OPEN (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN',
     .      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=RECLEN)

C       Make the PDS part of the header
      LABRECS = MAXHEAD/RECLEN + 1
      FILESIZE = (LABRECS+NL)*RECLEN - 20
      CRLF(1) = CHAR(13)
      CRLF(2) = CHAR(10)
      I = 1
      WRITE (BUF(I:I+79),'(A11,I9.9,A28,2A1)')  
     .  'NJPL1I00PDS',FILESIZE,'            = PDS_SFDU_LABEL', CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(A39,2A1)') 
     .  'FILE_TYPE                       = IMAGE', CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(A46,2A1)')
     .  'RECORD_TYPE                     = FIXED_LENGTH', CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),11)
     .  'RECORD_BYTES                    = ', RECLEN, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),11)
     .  'FILE_RECORDS                    = ', LABRECS+NL, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),11)
     .  'LABEL_RECORD_BYTES              = ', RECLEN, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),11)
     .  'LABEL_RECORDS                   = ', LABRECS, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),11)
     .  'IMAGE_RECORD_BYTES              = ', RECLEN, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),11)
     .  'IMAGE_RECORDS                   = ', NL, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),11)
     .  'IMAGE_LINES                     = ', NL, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),11)
     .  'LINE_SAMPLES                    = ', NS, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),11)
     .  'SAMPLE_BITS                     = ', NBYTES*8, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(A30,2A1)') 
     .  'SHDOM_VISUALIZATION_HEADER = "', CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
11    FORMAT (A34,I4,2A1)

C       Make the SHDOM part of the header
      GRIDNAME = 'EVEN (X,Y)  '
      IF (GRIDTYPE .EQ. 'P') THEN
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
      ELSE IF (GRIDTYPE .EQ. 'F') THEN
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
      ELSE
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
      ENDIF
 
      IF (SRCTYPE .EQ. 'S') THEN
          SOURCENAME = 'SOLAR'
      ELSE IF (SRCTYPE .EQ. 'T') THEN
          SOURCENAME = 'THERMAL'
      ELSE IF (SRCTYPE .EQ. 'B') THEN
          SOURCENAME = 'SOLAR/THERMAL'
      ENDIF

      IF (SFCTYPE .EQ. 'FL') THEN 
        SFCNAME = 'FIXED LAMBERTIAN'
      ELSE IF (SFCTYPE .EQ. 'VL') THEN 
        SFCNAME = 'VARIABLE LAMBERTIAN'
      ELSE IF (SFCTYPE .EQ. 'VF') THEN 
        SFCNAME = 'VARIABLE FRESNEL'
      ELSE IF (SFCTYPE .EQ. 'VR') THEN 
        SFCNAME = 'VARIABLE RPV'
      ELSE IF (SFCTYPE(1:1) .EQ. 'V') THEN 
        SFCNAME = 'VARIABLE OTHER'
      ENDIF
 
      IF (UNITS .EQ. 'T') THEN
        UNITSNAME = 'KELVIN'
      ELSE IF (UNITS .EQ. 'B') THEN
        UNITSNAME = 'WATTS/(M^2 STER)'
      ELSE
        UNITSNAME = 'WATTS/(M^2 MICRON STER)'
      ENDIF

      WRITE (BUF(I:I+79),'(2(A,I3),A,I5,2(A,I3),A,I5,A,I10,2A1)') 
     .    '!  L=',ML,'  M=',MM, '  NLM=',NLM, 
     .    '   NMU=',NMU, '  NPHI=',NPHI, '  NANG=', NANG, 
     .    '   NSH=', NSH, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(3(A,I4),2(A,I8),2A1)') 
     .  '!  NX=',NX,'   NY=',NY, '  NZ=',NZ, 
     .  '    NPTS=',NPTS, '   NCELLS=',NCELLS, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(A,A32,2A1)')
     .  '!  PROPERTY_FILE=', PROPFILE, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(A,A32,A,I2,2A1)')
     .  '!  CORRELATED_K-DIST_FILE=', CKDFILE, '   NUM_G=', NG, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(A,A32,2A1)')
     .  '!  INPUT_SAVE_FILE=', INSAVEFILE, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(A,A32,2A1)')
     .  '!  OUTPUT_SAVE_FILE=', OUTSAVEFILE, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      IF (DELTAM) THEN
        WRITE (BUF(I:I+79),'(A,A14,A,2A1)')
     .    '!  SOURCE_TYPE=', SOURCENAME, '    DELTA-M METHOD', CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
      ELSE
        WRITE (BUF(I:I+79),'(A,A14,2A1)')
     .    '!  SOURCE_TYPE=', SOURCENAME, CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
      ENDIF
      WRITE (BUF(I:I+79),'(A,A32,A,I1,2A1)') '!  GRID_TYPE=', GRIDNAME,
     .                         '   INDEPENDENT_PIXEL=', IPFLAG, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(A,A20,A,I1,2A1)') 
     . '!  SURFACE_TYPE=',SFCNAME,'   HORIZ_BOUNDARY_COND=',BCFLAG,CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
        IF (SFCTYPE(1:1) .EQ. 'V') THEN
          WRITE (BUF(I:I+79),'(A,A32,A,E12.5,2A1)')
     .        '!  SURFACE_FILE=', SFCFILE,
     .        '  SKY_RAD=', SKYRAD, CRLF
          I = I + INDEX(BUF(I:I+79),CHAR(10))
        ELSE
          IF (SRCTYPE .EQ. 'B') THEN
            WRITE (BUF(I:I+79),'(A,F9.7,A,F8.3,A,E12.5,2A1)') 
     .        '!  GROUND_ALBEDO=',GNDALBEDO, '  GROUND_TEMP=',GNDTEMP,
     .        '  SKY_RAD=',SKYRAD, CRLF
            I = I + INDEX(BUF(I:I+79),CHAR(10))
          ELSE
            WRITE (BUF(I:I+79),'(A,F9.7,A,E12.5,2A1)') 
     .        '!  GROUND_ALBEDO=',GNDALBEDO, '  SKY_RAD=',SKYRAD, CRLF
            I = I + INDEX(BUF(I:I+79),CHAR(10))
          ENDIF
        ENDIF
        WRITE (BUF(I:I+79),'(A,E13.6,A,F10.7,A,F8.3,2A1)')
     .         '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,
     .         '   SOLAR_AZ=', SOLARAZ*180.0/ACOS(-1.0), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
      ELSE
        IF (SFCTYPE(1:1) .EQ. 'V') THEN
          WRITE (BUF(I:I+79),'(A,A32,2A1)')
     .      '!  SURFACE_FILE=', SFCFILE, CRLF
          I = I + INDEX(BUF(I:I+79),CHAR(10))
        ELSE
          WRITE (BUF(I:I+79),'(A,F8.3,A,F9.7,2A1)') 
     .      '!  GROUND_TEMP=',GNDTEMP, 
     .      '  GROUND_EMIS=',1.0-GNDALBEDO, CRLF
          I = I + INDEX(BUF(I:I+79),CHAR(10))
        ENDIF
        WRITE (BUF(I:I+79),'(A,F8.3,2A1)') '!  SKY_TEMP=', SKYRAD,CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
      ENDIF

      IF (UNITS .EQ. 'T') THEN
        WRITE (BUF(I:I+79),'(A,A24,2A1)')   '!  UNITS=',UNITSNAME,CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
      ELSE IF (UNITS .EQ. 'B') THEN
        WRITE (BUF(I:I+79),'(A,A24,A,2(1X,F10.2),2A1)') 
     .    '!  UNITS=', UNITSNAME,
     .    '   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
      ELSE
        WRITE (BUF(I:I+79),'(A,A24,A,F10.2,2A1)') 
     .    '!  UNITS=',UNITSNAME, '   WAVELENGTH=', WAVELEN, CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
      ENDIF
      WRITE (BUF(I:I+79),'(2(A,E10.3),2A1)') 
     .  '!  SPLITTING_ACCURACY=', SPLITACC,
     .  '   SPHERICAL_HARMONIC_ACCURACY=',SHACC, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(A,E10.3,2A1)') 
     .  '!  SOLUTION_ACCURACY=', SOLACC, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+79),'(2(A,I4),2A1)') 
     .  '!  MAXIMUM_ITERATIONS=', MAXITER,
     .  '   NUMBER_ITERATIONS=', TOTITER, CRLF
      I = I + INDEX(BUF(I:I+79),CHAR(10))
      WRITE (BUF(I:I+7),'(A,2A1)') '!  ', CRLF
      I = I + INDEX(BUF(I:I+7),CHAR(10))

      IF (CAMERA_MODE) THEN
        WRITE (BUF(I:I+79),'(A,2A1)') '!  CAMERA MODE:', CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
C          1, bytes, scale, X,Y,Z, theta, phi, rotang, NL, NS, delline, delsamp
        WRITE (BUF(I:I+79),'(A,F9.2,2A1)') 
     .      '!  RADIANCE_SCALE=', OUTPARMS(3), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,3(1X,F7.3),2A1)') '!  CAMERA_LOCATION=',
     .      OUTPARMS(4),OUTPARMS(5),OUTPARMS(6),  CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,F10.5,2A1)')
     .      '!  CAMERA_ZENITH_ANGLE=', OUTPARMS(7), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,F10.5,2A1)')
     .      '!  CAMERA_AZIMUTH_ANGLE=', OUTPARMS(8), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,F10.5,2A1)')
     .      '!  CAMERA_ROTATION_ANGLE=', OUTPARMS(9), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,F8.5,2A1)') 
     .      '!  CAMERA_LINE_SPACING_ANGLE=', OUTPARMS(12), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,F8.5,2A1)') 
     .      '!  CAMERA_SAMPLE_SPACING_ANGLE=', OUTPARMS(13), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))

      ELSE
        WRITE (BUF(I:I+79),'(A,2A1)') 
     .      '!  CROSS TRACK SCANNING MODE:', CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
C          2, bytes, scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1, scan2, delscan
        WRITE (BUF(I:I+79),'(A,F9.2,2A1)') 
     .      '!  RADIANCE_SCALE=', OUTPARMS(3), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,3(1X,F7.3),2A1)') '!  START_LOCATION=',
     .      OUTPARMS(4),OUTPARMS(5),OUTPARMS(6), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,3(1X,F7.3),2A1)') '!  END_LOCATION=',
     .      OUTPARMS(7),OUTPARMS(8),OUTPARMS(9), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,F8.4,2A1)') 
     .      '!  SPACING_ON_TRACK=', OUTPARMS(10), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,F8.3,2A1)') 
     .      '!  START_SCAN_ANGLE=',OUTPARMS(11), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,F8.3,2A1)') 
     .      '!  END_SCAN_ANGLE=',OUTPARMS(12), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
        WRITE (BUF(I:I+79),'(A,F7.4,2A1)') 
     .      '!  SCAN_ANGLE_SPACING=',OUTPARMS(13), CRLF
        I = I + INDEX(BUF(I:I+79),CHAR(10))
      ENDIF

      WRITE (BUF(I:I+7),'(A1,2A1,A3,2A1)') '"', CRLF, 'END', CRLF

C       Write the header records
      DO I = 1, LABRECS
        J = RECLEN*(I-1)
        WRITE (2, REC=I) (BUF(K:K), K=J+1,MIN(J+RECLEN,MAXHEAD))
      ENDDO

C       Output the image, converting from real to byte or integer
      K = IVIS+1
      DO J = 1, NL
        IF (NBYTES .EQ. 1) THEN
          DO I = 1, NS
            DN = NINT(OUTPARMS(3)*VISOUT(K))
            DN = MAX(0,MIN(255,DN))
            IF (DN .GT. 127) THEN
              BUFFER1(I) = DN - 256
            ELSE
              BUFFER1(I) = DN
            ENDIF
            K = K + 1
          ENDDO
          WRITE (2,REC=LABRECS+J) (BUFFER1(I), I=1,NS)
        ELSE
          DO I = 1, NS
            BUFFER2(I) = NINT(OUTPARMS(3)*VISOUT(K))
            K = K + 1
          ENDDO
          WRITE (2,REC=LABRECS+J) (BUFFER2(I), I=1,NS)
        ENDIF
      ENDDO
      IVIS = IVIS + NL*NS

      CLOSE (2)
      RETURN
      END