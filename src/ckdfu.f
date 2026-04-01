!       PROGRAM CKDFU
! C       Prepares a correlated k-distribution file (for use with SHDOM)
! C     from Qiang Fu's correlated k-distribution code.
! C     Reads in a McClatchey atmosphere profile and computes the 
! C     absorption coefficient for each band, each "g" value (or k), and each
! C     level over the 18 bands of the shortwave and longwave regions.
! C       There are three sections of the output file.
! C     1. The band information: number of bands and for each band
! C        the wavenumber range (cm^-1), solar flux (W/m^2), 
! C        number of g's (or k's), and delta g values (weights).
! C     2. The atmosphere information: the concentration of CO2, CH4, N2O,
! C        and for each level (from the top down) the height (km), 
! C        pressure (mb), temperature (K), mixing ratio of water vapor and ozone.
! C     3. The absorption coefficients: each line contains the
! C        band number, level number (from top down), and k's for the g's.

!       INTEGER MAXNZ, MAXNB, MAXNG, NB
!       PARAMETER (MAXNZ=201, MAXNB=27, MAXNG=12, NB=15)
!       INTEGER NLEV, NZ
!       INTEGER I, IB, IG
!       REAL    KABS(MAXNZ,MAXNG,MAXNB),  DELG(MAXNG,MAXNB)
!       REAL    HEIGHT(MAXNZ), PRES(MAXNZ), TEMP(MAXNZ)
!       REAL    RHOAIR(MAXNZ), RHOH2O(MAXNZ), RHOO3(MAXNZ)
!       REAL    KGC(MAXNZ), KG(MAXNZ)
!       REAL    PP(MAXNZ), PT(MAXNZ), PH(MAXNZ), PO(MAXNZ),  Z(MAXNZ)
!       REAL    UMCO2, UMCH4, UMN2O
!       CHARACTER  ATMFILE*72, OUTFILE*72

!       INTEGER NG(27)
!       REAL    WAVENO(2,27), SOLARFLUX(27)
!       DATA NG/1,1,1,1,1,1,1,1,1,1,8,12,7,12,5,2,3,4,4,3,5,2,10,12,7,7,8/
!       DATA WAVENO/44500,57000,41000,44500,35000,41000,33500,35000,
!      .            31008,33500,27972,31008,22857,27972,20101,22857,
!      .            16807,20101,14500,16807,
!      .            7700,14500,5250,7700,4000,5250,
!      .            2850,4000,2500,2850,1900,2200,1700,1900,1400,1700,
!      .            1250,1400,1100,1250,980,1100,800,980,670,800,540,670,
!      .            400,540,280,400,0,280/
! C           1361 W/m2
!       DATA SOLARFLUX/0.81,0.93,6.25,6.09,15.27,33.41,109.81,120.54
!      .             ,179.06,151.55,489.46,152.52,50.36,28.07,5.45,12*0.0/

!       COMMON /ATMOS/ PP, PT, PH, PO, NZ
!       COMMON /UMCON/ UMCO2, UMCH4, UMN2O


!       UMCO2=330
!       UMCH4=1.6
!       UMN2O=0.28

!       WRITE (*,*) 'Atmospheric profile file '
!       READ (*,'(A)') ATMFILE

!       WRITE (*,*) 'Output k-distribution file '
!       READ (*,'(A)') OUTFILE

!       WRITE (*,'(A,A)') ' Concentration of CO2, CH4, and N2O',
!      .                  ' (ppmv) (try: 360, 1.7, 0.3)'
!       READ (*,*) UMCO2, UMCH4, UMN2O


! C           Read in the McClatchey atmospheric profile
!       CALL READ_ATMOS (ATMFILE, NLEV, HEIGHT, PRES, TEMP,
!      .                 RHOAIR, RHOH2O, RHOO3)
!       IF (NLEV .GT. MAXNZ) STOP 'MAXNZ exceeded'

! C           Interpolate the McClatchey atmosphere to the atmos profile arrays
!       DO I = 1, NLEV
!         Z(I) = HEIGHT(I)
!         PP(I) = PRES(I)
!         PT(I) = TEMP(I)
!         PH(I) = RHOH2O(I)/RHOAIR(I)
!         PO(I) = RHOO3(I)/RHOAIR(I)
!       ENDDO
!       NZ = NLEV


! C         Loop over bands
!       DO IB = 1, NB
! C           Get the water vapor continuum for this band
!         CALL GASCON(IB, KGC)
! C           Loop over the g's for this band
!         DO IG = 1, NG(IB)
! C             Get the layer optical depths and compute the absorption coef.
!           CALL GASES (IB, IG, DELG(IG,IB), KG)
!           DO I = 1, NZ
!             KABS(I,IG,IB) = KG(I) + KGC(I)
! 			IF ((I .eq. 1) .and. (IB .eq. 1)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 3.8336E+00
! 			ELSEIF ((I .eq. 1) .and. (IB .eq. 2)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 1.8233E+01
! 			ELSEIF ((I .eq. 1) .and. (IB .eq. 3)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 2.0454E+01
! 			ELSEIF ((I .eq. 1) .and. (IB .eq. 4)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 3.1557E+00
! 			ELSEIF ((I .eq. 1) .and. (IB .eq. 5)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 2.5597E-01
! 			ELSEIF ((I .eq. 1) .and. (IB .eq. 6)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 9.7944E-03
! 			ELSEIF ((I .eq. 1) .and. (IB .eq. 7)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 8.6257E-05
! 			ELSEIF ((I .eq. 1) .and. (IB .eq. 8)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 1.4610E-03
! 			ELSEIF ((I .eq. 1) .and. (IB .eq. 9)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 9.8763E-03
! 			ELSEIF ((I .eq. 1) .and. (IB .eq. 10)) THEN
! 			  KABS(I,IG,IB) = KABS(I,IG,IB) + 9.2101E-03
! 			ENDIF
!             IF (KABS(I,IG,IB) .LT. 0.0) THEN
!               WRITE (*,'(1X,A,A,E10.3)') 'Surely something is wrong,',
!      .               ' negative absorption: ',KABS(I,IG,IB)
!               WRITE (*,'(1X,A,I2,A,I3,A,I2)')
!      .           'Band: ',IB, '   Level:',I, '   k:',IG
!               STOP
!             ENDIF
!           ENDDO
!         ENDDO
!       ENDDO



! C         Output CKD file
!       OPEN (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')

!       WRITE (2,'(A)') '! Correlated k-distribution file: method of Fu'
!       WRITE (2,'(I3,A)') NB, '  ! number of bands'
!       WRITE (2,'(A,A)') '!Band  Wavenums    SolF  Ng  Delta g'
!       DO IB = 1, NB
!         WRITE (2,'(1X,I2,2(1X,F6.0),1X,F6.2,1X,I2,12(1X,F6.4))')
!      .         IB, WAVENO(1,IB), WAVENO(2,IB), SOLARFLUX(IB),
!      .         NG(IB), (DELG(IG,IB), IG=1, NG(IB))
!       ENDDO

!       WRITE (2,'(I3,A)') NZ, '  ! number of levels'
!       WRITE (2,'(1X,F5.1,2(1X,F5.3),A)')  UMCO2, UMCH4, UMN2O,  
!      .       '  ! concentration of CO2, CH4, N2O in ppmv'
!       WRITE (2,'(A)') '!  Z     Pres  Temp      qH2O        qO3'
!       DO I = 1, NZ
!         WRITE (2,'(1X,F7.3,1X,F6.1,1X,F5.1,2(1X,E11.4))')
!      .        Z(I), PP(I), PT(I), PH(I), PO(I)

!       ENDDO

!       WRITE (2,'(A)') '!IB IZ    Kabs (km^-1)'
!       DO IB = 1, NB
!         DO I = 1, NZ
!           WRITE (2,'(1X,I2,1X,I3,12(1X,E12.5))')
!      .        IB, I, (KABS(I,IG,IB), IG=1, NG(IB))
!         ENDDO
!       ENDDO
!       CLOSE (2)      

!       END


	SUBROUTINE CKDFU(NLEV,HEIGHT,PRES,TEMP,RHOH2O,RHOO3,
     .               UMCO2_IN, UMCH4_IN,UMN2O_IN,DELG,
     .				 KABS)


	  IMPLICIT NONE
      INTEGER MAXNZ, MAXNB, MAXNG, NB
      PARAMETER (MAXNZ=201,MAXNB=27, MAXNG=12, NB=15)
      INTEGER NLEV, NZ
      INTEGER I, IB, IG
      REAL    KABS(MAXNZ,MAXNG,MAXNB),  DELG(MAXNG,MAXNB)
!f2py intent(out) :: KABS, DELG
      REAL    HEIGHT(NLEV), PRES(NLEV), TEMP(NLEV)
!f2py intent(in) :: HEIGHT, PRES, TEMP
      REAL    RHOH2O(NLEV), RHOO3(NLEV)
!f2py intent(in) :: RHOH2O, RHOO3
      REAL    KGC(MAXNZ), KG(MAXNZ)
      REAL    PP(MAXNZ), PT(MAXNZ), PH(MAXNZ), PO(MAXNZ),  Z(MAXNZ)
      REAL    UMCO2_IN, UMCH4_IN, UMN2O_IN
!f2py intent(in) :: UMCO2_IN, UMCH4_IN, UMN2O_IN
      REAL    UMCO2, UMCH4, UMN2O
      INTEGER NG(MAXNB)
      REAL    WAVENO(2,MAXNB), SOLARFLUX(MAXNB)
!f2py intent(out) :: SOLARFLUX
      DATA NG/1,1,1,1,1,1,1,1,1,1,8,12,7,12,5,2,3,4,4,3,5,2,10,12,7,7,8/
      DATA WAVENO/44500,57000,41000,44500,35000,41000,33500,35000,
     .            31008,33500,27972,31008,22857,27972,20101,22857,
     .            16807,20101,14500,16807,
     .            7700,14500,5250,7700,4000,5250,
     .            2850,4000,2500,2850,1900,2200,1700,1900,1400,1700,
     .            1250,1400,1100,1250,980,1100,800,980,670,800,540,670,
     .            400,540,280,400,0,280/
C           1361 W/m2
      DATA SOLARFLUX/0.81,0.93,6.25,6.09,15.27,33.41,109.81,120.54
     .             ,179.06,151.55,489.46,152.52,50.36,28.07,5.45,12*0.0/

      COMMON /ATMOS/ PP, PT, PH, PO, NZ
      COMMON /UMCON/ UMCO2, UMCH4, UMN2O
	
	  UMCO2 = UMCO2_IN
	  UMCH4 = UMCH4_IN
	  UMN2O = UMN2O_IN

      DO I = 1, NLEV
        Z(I) = HEIGHT(I)
        PP(I) = PRES(I)
        PT(I) = TEMP(I)
        PH(I) = RHOH2O(I)!/RHOAIR(I)
        PO(I) = RHOO3(I)!/RHOAIR(I)
      ENDDO
      NZ = NLEV

C         Loop over bands
      DO IB = 1, NB
C           Get the water vapor continuum for this band
        CALL GASCON(IB, KGC)
C           Loop over the g's for this band
        DO IG = 1, NG(IB)
C             Get the layer optical depths and compute the absorption coef.
          CALL GASES (IB, IG, DELG(IG,IB), KG)
          DO I = 1, NZ
            KABS(I,IG,IB) = KG(I) + KGC(I)
			IF ((I .eq. 1) .and. (IB .eq. 1)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 3.8336E+00
			ELSEIF ((I .eq. 1) .and. (IB .eq. 2)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 1.8233E+01
			ELSEIF ((I .eq. 1) .and. (IB .eq. 3)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 2.0454E+01
			ELSEIF ((I .eq. 1) .and. (IB .eq. 4)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 3.1557E+00
			ELSEIF ((I .eq. 1) .and. (IB .eq. 5)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 2.5597E-01
			ELSEIF ((I .eq. 1) .and. (IB .eq. 6)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 9.7944E-03
			ELSEIF ((I .eq. 1) .and. (IB .eq. 7)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 8.6257E-05
			ELSEIF ((I .eq. 1) .and. (IB .eq. 8)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 1.4610E-03
			ELSEIF ((I .eq. 1) .and. (IB .eq. 9)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 9.8763E-03
			ELSEIF ((I .eq. 1) .and. (IB .eq. 10)) THEN
			  KABS(I,IG,IB) = KABS(I,IG,IB) + 9.2101E-03
			ENDIF
            IF (KABS(I,IG,IB) .LT. 0.0) THEN
              WRITE (*,'(1X,A,A,E10.3)') 'Surely something is wrong,',
     .               ' negative absorption: ',KABS(I,IG,IB)
              WRITE (*,'(1X,A,I2,A,I3,A,I2)')
     .           'Band: ',IB, '   Level:',I, '   k:',IG
              STOP
            ENDIF
          ENDDO
        ENDDO
      ENDDO

	END SUBROUTINE




C     Correlated k-distribution routines and tables for the method
C   described in Qiang Fu and K. N. Liou, 1992: On the correlated
C   k-distribution method for radiative transfer in nonhomogeneous 
C   atmospheres, J. Atmos. Sci., 49, 2139-2156.


C Band  Region (/cm)  N g's  Gases  /Overlap method
C   1  50000-14500    10       O3   
C   2  14500-7700      8       H2O  
C   3   7700-5250     12       H2O  
C   4   5250-4000      7       H2O  
C   5   4000-2850     12       H2O  
C   6   2850-2500      5       H2O  
C
C   7   2200-1900      2       H2O
C   8   1900-1700      3       H2O
C   9   1700-1400      4       H2O
C  10   1400-1250      4       H2O, CH4, N2O    1
C  11   1250-1100      3       H2O, CH4, N2O    1
C  12   1100- 980      5       H2O, O3          1
C  13    980- 800      2       H2O
C  14    800- 670     10       H2O, CO2         2
C  15    670- 540     12       H2O, CO2         2
C  16    540- 400      7       H2O 
C  17    400- 280      7       H2O
C  18    280- 000      8       H2O


C  Input and output are through common blocks.
C  Input:
C    Atmospheric profile specified in 
C     common /atmos/ pp(MAXNZ),pt(MAXNZ),ph(MAXNZ),po(MAXNZ), nz
C    Minor gas concentrations specified in /umcon/ umco2, umch4, umn2o.
C
C    integer ng(18)/10,8,12,7,12,5,2,3,4,4,3,5,2,10,12,7,7,8/
C  Call:
C    gases (ib, ig, hk, kg)
C       Input is band number (1-18) ib and g number ig (1-nk(ib)).
C       Output is g interval (weight) hk and gaseous 
C         absorption coefficient at levels kg(nz)
C    gascon(ib, kgc)  
C       Input is band number (1-18) ib. 
C       Output is continuum absorption coefficient at levels kgc(nz).


C   Common blocks:
C     common /atmos/ pp(MAXNZ),pt(MAXNZ),ph(MAXNZ),po(MAXNZ), nz
C       nz levels of pressure (mb), temperature (K), 
C       humidity (water vapor mixing ratio g/g), ozone (mixing ratio g/g)
C       in order of increasing pressure.  
C     /umcon/ umco2, umch4, umn2o   
C       Concentrations of CO2, CH4, N2O in ppmv.
C       "Standard" values are 330, 1.6, and 0.28.
C
C     /bandN/ delta g, cNG(nT,nP,ng)
C      Nth band (N=1-18), G=gas (o3,h2o,co2,n2o,ch4), 
C      nT = num temperatures, nP = num pressures, ng = num g's

C   Subroutines:
C    gascon(ib, kgc)  
C       Computes absorption coefficient for water vapor continuum for ib band
C    qopcon ( vv, kgc )    
C       Does the actual continuum calculation

C    gases ( ib, ig, hk )
C       Computes volume absorption coefficient (km^-1) from gases 
C       for ib band, for ig k,

C    qks ( coefks, fkg )
C       Computes the gaseous absorption coefficients in units 
C       of (cm-atm)**-1 for a given cumulative probability in nz levels. 
C       coefks(3,11) are the coefficients to calculate the absorption 
C       coefficient at the temperature T for the 11 pressures.
C       Interpolates to the level values of temperature and pressure.
C       This is called for solar bands (1-6).  (T0=245K)
C    qki ( coefki, fkg )      Similar to qks, but for infrared bands.
C    qkio3 ( coefki, fkg )    Similar to qkio3, but for ozone band (T0=250K)

C       The following routine compute the absorption coefficient (km^-1)
C       for each level for various gases from the absorption coefs fkg.
C    qoph2o ( fkg, kg )       For water vapor
C    qopo3i ( fkg, kg1 )      For ozone
C    qopch4 ( fkg, kg2 )      For methane
C    qopn2o ( fkg, kg3 )      For nitrous oxide
C    qophc ( fkg, kg )        For CO2, H2O overlap 
   




	subroutine gascon ( ib, kgc )
c *********************************************************************
c kgc(nz) are the volume absorption coefficients (km^-1) to water 
c vapor continuum absorption at nz levels for a given band ib. 
c We include continuum absorption in the 280 to 1250 cm**-1 region. 
c vv(11)-vv(17) are the central wavenumbers of each band in this region. 
c *********************************************************************
        real kgc(*)
        integer MAXNZ, nz, ib, i
        real    pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
        real vv(27)
	data vv / 19*0.0, 1175.0, 1040.0, 890.0, 735.0, 
     1	          605.0, 470.0, 340.0, 0.0 /
	if ( ib .gt. 19 .and. ib .lt. 27 ) then
	   call qopcon ( vv(ib), kgc )
	else
	   do 10 i = 1, nz
              kgc(i) = 0.0
10	   continue
	endif
	return
	end



	subroutine gases (ib, ig, hk, kg)
c *********************************************************************
c kg(nz) are the volume absorption coefficients (km^-1) due to nongray gaseous 
c absorption, at nz levels for a given band ib and cumulative probability ig. 
c *********************************************************************
        real hk, kg(nz)
        integer MAXNZ, nz, ib, ig,  i
        real    pp, pt, ph, po, umco2, umch4, umn2o
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
        common /umcon/ umco2, umch4, umn2o 
	real    fk, fkg(MAXNZ), fkga(MAXNZ), fkgb(MAXNZ), pq(MAXNZ)
	real    kg1(MAXNZ), kg2(MAXNZ), kg3(MAXNZ)
        real    hk101,hk102,hk103,hk104,hk105
        real    hk106,hk107,hk108,hk109,hk110
        real    hk2, hk3, hk4, hk5, hk6, hk7, hk8, hk9
        real    hk10, hk11, hk12, hk13, hk14, hk15, hk16, hk17, hk18
        real    fk101o3,fk102o3,fk103o3,fk104o3,fk105o3 
        real    fk106o3,fk107o3,fk108o3,fk109o3,fk110o3
        real    c2h2o, c3h2o, c4h2o, c5h2o, c6h2o
        real    c7h2o, c8h2o, c9h2o, c10h2o, c10ch4, c10n2o
        real    c11h2o, c11ch4, c11n2o, c12o3, c12h2o, c13h2o
        real    c14hca, c14hcb, c15hca, c15hcb
        real    c16h2o, c17h2o, c18h2o
	common /band101/ hk101(1), fk101o3(1)
	common /band102/ hk102(1), fk102o3(1)
	common /band103/ hk103(1), fk103o3(1)
	common /band104/ hk104(1), fk104o3(1)
	common /band105/ hk105(1), fk105o3(1)
	common /band106/ hk106(1), fk106o3(1)
	common /band107/ hk107(1), fk107o3(1)
	common /band108/ hk108(1), fk108o3(1)
	common /band109/ hk109(1), fk109o3(1)
	common /band110/ hk110(1), fk110o3(1)
	common /band2/ hk2(8), c2h2o(3,11,8)
	common /band3/ hk3(12), c3h2o(3,11,12)
	common /band4/ hk4(7), c4h2o(3,11,7)
	common /band5/ hk5(12), c5h2o(3,11,12)
	common /band6/ hk6(5), c6h2o(3,11,5)
	common /band7/ hk7(2), c7h2o(3,19,2)
	common /band8/ hk8(3), c8h2o(3,19,3)
	common /band9/ hk9(4), c9h2o(3,19,4)
	common /band10/ hk10(4),c10h2o(3,19,4),c10ch4(3,19),c10n2o(3,19)
	common /band11/ hk11(3),c11h2o(3,19,3),c11ch4(3,19),c11n2o(3,19)
	common /band12/ hk12(5), c12o3(3,19,5), c12h2o(3,19)
	common /band13/ hk13(2), c13h2o(3,19,2)
	common /band14/ hk14(10), c14hca(3,19,10), c14hcb(3,19,10)
	common /band15/ hk15(12), c15hca(3,19,12), c15hcb(3,19,12)
	common /band16/ hk16(7), c16h2o(3,19,7)
	common /band17/ hk17(7), c17h2o(3,19,7)
	common /band18/ hk18(8), c18h2o(3,19,8)

	goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) ib
	stop 'Band out of range'

c In this band ( 50000 - 14500 cm**-1 ), we have considered the nongray
c gaseous absorption of O3.
1	fk = fk101o3(ig)
	call qopo3s ( fk, kg )
	hk = hk101(ig)
	goto 30

2	fk = fk102o3(ig)
	call qopo3s ( fk, kg )
	hk = hk102(ig)
	goto 30

3	fk = fk103o3(ig)
	call qopo3s ( fk, kg )
	hk = hk103(ig)
	goto 30

4	fk = fk104o3(ig)
	call qopo3s ( fk, kg )
	hk = hk104(ig)
	goto 30

5	fk = fk105o3(ig)
	call qopo3s ( fk, kg )
	hk = hk105(ig)
	goto 30

6	fk = fk106o3(ig)
	call qopo3s ( fk, kg )
	hk = hk106(ig)
	goto 30

7	fk = fk107o3(ig)
	call qopo3s ( fk, kg )
	hk = hk107(ig)
	goto 30 

8	fk = fk108o3(ig)
	call qopo3s ( fk, kg )
	hk = hk108(ig)
	goto 30

9	fk = fk109o3(ig)
	call qopo3s ( fk, kg )
	hk = hk109(ig)
	goto 30

10	fk = fk110o3(ig)
	call qopo3s ( fk, kg )
	hk = hk101(ig)
	goto 30

c In this band ( 14500 - 7700 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.  
11	call qks ( c2h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk2(ig)
	goto 30

c In this band ( 7700 - 5250 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O. 
12	call qks ( c3h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk3(ig)
	goto 30

c In this band ( 5250 - 4000 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O. 
13	call qks ( c4h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk4(ig)
	goto 30

c In this band ( 4000 - 2850 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O. 
14	call qks ( c5h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk5(ig)
	goto 30

c In this band ( 2850 - 2500 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O. 
15	call qks ( c6h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk6(ig)
	goto 30

c In this band ( 2200 - 1900 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
16	call qki ( c7h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk7(ig)
	goto 30

c In this band ( 1900 - 1700 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
17	call qki ( c8h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk8(ig)
	goto 30

c In this band ( 1700 - 1400 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
18	call qki ( c9h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk9(ig)
	goto 30

c In this band ( 1400 - 1250 cm**-1 ), we have considered the overlapping
c absorption of H2O, CH4, and N2O by approach one of Fu(1991).
19	call qki ( c10h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg1 )
	call qki ( c10ch4, fkg )
	call qopch4 ( fkg, kg2 )
	call qki ( c10n2o, fkg )
	call qopn2o ( fkg, kg3 )
	do 205 i = 1, nz
	   kg(i) = kg1(i) + kg2(i)*umch4 + kg3(i)*umn2o
205	continue
	hk = hk10(ig)
	goto 30

c In this band ( 1250 - 1100 cm**-1 ), we have considered the overlapping
c absorption of H2O, CH4, and N2O by approach one of Fu(1991).
20	call qki ( c11h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg1 )
	call qki ( c11ch4, fkg )
	call qopch4 ( fkg, kg2 )
	call qki ( c11n2o, fkg )
	call qopn2o ( fkg, kg3 )
	do 215 i = 1, nz
	   kg(i) = kg1(i) + kg2(i)*umch4 + kg3(i)*umn2o
215	continue
	hk = hk11(ig)
	goto 30

c In this band ( 1100 - 980 cm**-1 ), we have considered the overlapping
c absorption of H2O and O3 by approach one of Fu(1991).
21	call qkio3 ( c12o3(1,1,ig), fkg )
	call qopo3i ( fkg, kg1 )
	call qki ( c12h2o, fkg )
	call qoph2o ( fkg, kg2 )
	do 225 i = 1, nz
	   kg(i) = kg1(i) + kg2(i)
225	continue
	hk = hk12(ig)
	goto 30

c In this band ( 980 - 800 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
22	call qki ( c13h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk13(ig)
	goto 30

c In this band ( 800 - 670 cm**-1), we have considered the overlapping
c absorption of H2O and CO2 by approach two of Fu(1991).
23      do 333 i = 1, nz
	   if ( pp(i) .ge. 63.1 ) then
             pq(i) = ph(i)
           else
             pq(i) = 0.0
       	   endif
333	continue
	call qki ( c14hca(1,1,ig), fkga )
	call qki ( c14hcb(1,1,ig), fkgb )
	do 343 i = 1, nz
	   fkg(i) = fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
343	continue
	call qophc ( fkg, kg)
	hk = hk14(ig)
	goto 30

c In this band ( 670 - 540 cm**-1), we have considered the overlapping
c absorption of H2O and CO2 by approach two of Fu(1991).
24      do 353 i = 1, nz
	   if ( pp(i) .ge. 63.1 ) then
             pq(i) = ph(i)
           else
             pq(i) = 0.0
       	   endif
353	continue
	call qki ( c15hca(1,1,ig), fkga )
	call qki ( c15hcb(1,1,ig), fkgb )
	do 363 i = 1, nz
	   fkg(i) = fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
363	continue
	call qophc ( fkg, kg)
	hk = hk15(ig)
	goto 30

c In this band ( 540 - 400 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
25	call qki ( c16h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk16(ig)
	goto 30

c In this band ( 400 - 280 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
27	call qki ( c17h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk17(ig)
	goto 30

c In this band ( 280 - 000 cm**-1 ), we have considered the nongray
c gaseous absorption of H2O.
28	call qki ( c18h2o(1,1,ig), fkg )
	call qoph2o ( fkg, kg )
	hk = hk18(ig)

30	continue
	return
	end



	subroutine qks ( coefks, fkg )
c *********************************************************************
c fkg(nz) are the gaseous absorption coefficients in units of (cm-atm)
c **-1 for a given cumulative probability in nz layers. coefks(3,11)
c are the coefficients to calculate the absorption coefficient at the
c temperature t for the 11 pressures by
c         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
c and the absorption coefficient at conditions other than those eleven
c pressures is interpolated linearly with pressure (Fu, 1991).
c *********************************************************************
	real  coefks(3,11),  fkg(*)
        integer MAXNZ, nz, i, i1
        real pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	real stanp(11), x1, x2, y1
      	data stanp / 10.0, 15.8, 25.1, 39.8, 63.1, 100.0,
     1	             158.0, 251.0, 398.0, 631.0, 1000.0 /
	i1 = 1
	do 5 i = 1, nz
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefks(1,1) + coefks(2,1) * ( pt(i) - 245.0 )
     1       + coefks(3,1) * ( pt(i) - 245.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(11) ) then
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefks(1,10) + coefks(2,10) * ( pt(i) - 245.0 )
     1	     + coefks(3,10) * y1 )
    	     x2 = exp ( coefks(1,11) + coefks(2,11) * ( pt(i) - 245.0 )
     1	     + coefks(3,11) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(11) - stanp(10) )
     1	     * ( pp(i) - stanp(10) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefks(1,i1-1) + coefks(2,i1-1) * (pt(i)-245.0)
     1	     + coefks(3,i1-1) * y1 )
	     x2 = exp ( coefks(1,i1) + coefks(2,i1) * ( pt(i) - 245.0 )
     1	     + coefks(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )
     1	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end


	subroutine qki ( coefki, fkg )
c *********************************************************************
c fkg(nz) are the gaseous absorption coefficients in units of (cm-atm)
c **-1 for a given cumulative probability in nz layers. coefki(3,19)
c are the coefficients to calculate the absorption coefficient at the
c temperature t for the 19 pressures by
c         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
c and the absorption coefficient at  conditions  other  than  those 19
c pressures is interpolated linearly with pressure (Fu, 1991).
c *********************************************************************
	real coefki(3,19), fkg(*)
        integer MAXNZ, nz, i, i1
        real pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	real stanp(19), x1, x2, y1
	data stanp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 
     1	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,
     1	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /
	i1 = 1
	do 5 i = 1, nz
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefki(1,1) + coefki(2,1) * ( pt(i) - 245.0 )
     1       + coefki(3,1) * ( pt(i) - 245.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(19) ) then
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefki(1,18) + coefki(2,18) * ( pt(i) - 245.0 )
     1	     + coefki(3,18) * y1 )
    	     x2 = exp ( coefki(1,19) + coefki(2,19) * ( pt(i) - 245.0 )
     1	     + coefki(3,19) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )
     1	     * ( pp(i) - stanp(18) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefki(1,i1-1) + coefki(2,i1-1) * (pt(i)-245.0)
     1	     + coefki(3,i1-1) * y1 )
	     x2 = exp ( coefki(1,i1) + coefki(2,i1) * ( pt(i) - 245.0 )
     1	     + coefki(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )
     1	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end


	subroutine qkio3 ( coefki, fkg )
c *********************************************************************
c fkg(nz) are the gaseous absorption coefficients in units of (cm-atm)
c **-1 for a given cumulative probability in nz layers. coefki(3,19)
c are the coefficients to calculate the absorption coefficient at the
c temperature t for the 19 pressures by
c         ln k = a + b * ( t - 250 ) + c * ( t - 250 ) ** 2
c and the absorption coefficient at  conditions  other  than  those 19
c pressures is interpolated linearly with pressure (Fu, 1991).
c *********************************************************************
	real coefki(3,19), fkg(*)
        integer MAXNZ, nz, i, i1
        real pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	real stanp(19), x1, x2, y1
	data stanp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51, 
     1	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,
     1	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /
	i1 = 1
	do 5 i = 1, nz
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefki(1,1) + coefki(2,1) * ( pt(i) - 250.0 )
     1       + coefki(3,1) * ( pt(i) - 250.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(19) ) then
	     y1 = ( pt(i) - 250.0 ) * ( pt(i) - 250.0 )
	     x1 = exp ( coefki(1,18) + coefki(2,18) * ( pt(i) - 250.0 )
     1	     + coefki(3,18) * y1 )
    	     x2 = exp ( coefki(1,19) + coefki(2,19) * ( pt(i) - 250.0 )
     1	     + coefki(3,19) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )
     1	     * ( pp(i) - stanp(18) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( pt(i) - 250.0 ) * ( pt(i) - 250.0 )
	     x1 = exp ( coefki(1,i1-1) + coefki(2,i1-1) * (pt(i)-250.0)
     1	     + coefki(3,i1-1) * y1 )
	     x2 = exp ( coefki(1,i1) + coefki(2,i1) * ( pt(i) - 250.0 )
     1	     + coefki(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )
     1	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end



	subroutine qopo3s ( fk, kg )
	real fk, kg(*)
        integer MAXNZ, nz, i
        real pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	do 10 i = 1, nz
 	   kg(i) = po(i) * pp(i)/pt(i) * fk * 16265.
10	continue
c         16265 = 2.241E4/(M*287.05) *1E4, where M = 48.00 for O3.
	return
	end

	subroutine qoph2o ( fkg, kg )
	real fkg(*), kg(*)
        integer MAXNZ, nz, i
        real pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	do 10 i = 1, nz
	   kg(i) = fkg(i) * ph(i) * pp(i)/pt(i) *43336.
10	continue
c         43336 = 2.241E4/(M*287.05) *1E4, where M = 18.015 for H2O.
	return
	end

	subroutine qopch4 ( fkg, kg )
	real fkg(*), kg(*)
        integer MAXNZ, nz, i
        real pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	do 10 i = 1, nz
	   kg(i) = fkg(i) * pp(i)/pt(i) * 0.026954
10	continue
c         0.026954 = 2.241E4/(28.964*287.05) *1.0E-6 *1E4
	return
	end

	subroutine qopn2o ( fkg, kg )
	real fkg(*), kg(*)
        integer MAXNZ, nz, i
        real pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	do 10 i = 1, nz
	   kg(i) = fkg(i) * pp(i)/pt(i) * 0.026954
10	continue
c         0.026954 = 2.241E4/(28.964*287.05) *1.0E-6 *1E4
	return
	end

	subroutine qopo3i ( fkg, kg )
	real fkg(*), kg(*)
        integer MAXNZ, nz, i
        real pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	do 10 i = 1, nz
	   kg(i) = fkg(i) * po(i) * pp(i)/pt(i) * 16265.
10	continue
c         16265 = 2.241E4/(M*287.05) *1E4, where M = 48.00 for O3.
	return
	end

	subroutine qophc ( fkg, kg )
	real fkg(*), kg(*)
        integer MAXNZ, nz, i
        real pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	do 10 i = 1, nz
	   kg(i) = fkg(i) * pp(i)/pt(i) *34.16
10	continue
c         34.16 = 9.80665/287.05 *1E3
c See page 86 of Fu (1991).
	return
	end


	subroutine qopcon ( vv, kg )
	real vv, kg(*)
        integer MAXNZ, nz, i
        real x, y, z, r, s, pe, w, ff, pp, pt, ph, po
        PARAMETER (MAXNZ=201)
	common /atmos/ pp(MAXNZ), pt(MAXNZ), ph(MAXNZ), po(MAXNZ), nz
	x = 4.18
	y = 5577.8
	z = 0.00787
	r = 0.002
	s = ( x + y * exp ( - z * vv ) ) / 1013.25
	do 5 i = 1, nz
 	   pe = pp(i) * ph(i) / ( 0.622 + 0.378 * ph(i) )
	   w = exp ( 1800.0 / pt(i) - 6.08108 )
	   ff = s * ( pe + r * pp(i) ) * w
	   kg(i) = ff * ph(i) * pp(i)/pt(i) * 34.84
5	continue
c         34.84 = 1E4/287.05
	return
	end







      SUBROUTINE READ_ATMOS (ATMFILE, NLEV, HEIGHT, PRES, TEMP,
     .                       RHOAIR, RHOH2O, RHOO3)
C       Reads in a McClatchey atmospheric profile.  The levels are
C       in order of decreasing height.
C        NLEV        number of levels (interfaces)
C        HEIGHT      height of levels (km)
C        PRES        pressure of levels (mb)
C        TEMP        temperature of levels (Kelvin)
C        RHOAIR      air density at levels (g/m^3)
C        RHOH2O      water vapor density at levels (g/m^3)
C        RHOO3       ozone density at levels (g/m^3)
      INTEGER NLEV
      REAL    HEIGHT(*), PRES(*), TEMP(*)
      REAL    RHOAIR(*), RHOH2O(*), RHOO3(*)
      CHARACTER*(*) ATMFILE
      INTEGER I, MAXLEV
      PARAMETER (MAXLEV=201)

      OPEN (UNIT=1, FILE=ATMFILE, STATUS='OLD')
      I = 1
100   CONTINUE
          READ (1,*,END=120) HEIGHT(I), PRES(I), TEMP(I),
     .             RHOAIR(I), RHOH2O(I), RHOO3(I)
          I = I + 1
      IF (I .LT. MAXLEV) GOTO 100 
120   CONTINUE
      NLEV = I - 1
      CLOSE (1)
      RETURN
      END


