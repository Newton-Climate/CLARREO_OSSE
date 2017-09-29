      SUBROUTINE DISORT(CORINT, DELTAM, ONLYFL,                         &
     &     MSFLAG,ACCUR, PRNT, HEADER,                                  &
     &     MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI,                      &
     &     PLANK, USRANG,                                               &
     &     USRTAU, FBEAM, UMU, UMU0, PHI0, DTAUC, FISOT,                &
     &     ALBEDO, BTEMP, IBCND, NCOEF, NLYR,                           &
     &     NPHI, NSTR, NTAU, NUMU, SSALB, PHI,                          &
     &     PMOM, UTAU, TEMIS, TTEMP, TEMPER, WVNMLO,                    &
     &     WVNMHI, RFLDIR, RFLDN,                                       &
     &     FLDIR, FLDN,                                                 &
     &     FLUP, UAVG, DFDT, UU, U0U, ALBMED, TRNMED,                   &
     &     WN0, S0CMS, T0CMS, FDNSRT, FDNTRT)
! **********************************************************************
!       PLANE-PARALLEL DISCRETE ORDINATES RADIATIVE TRANSFER PROGRAM
!             (SEE DISORT.DOC FOR COMPLETE DOCUMENTATION)
! **********************************************************************

!+---------------------------------------------------------------------+
!------------------    I/O VARIABLE SPECIFICATIONS     -----------------
!+---------------------------------------------------------------------+
!MOD  PARAMETER (MAXCLY=35, MAXULV=35)
!MOD  PARAMETER (MAXUMU=8, MAXCMU=4, MAXCOE=10, MAXPHI=6)
      INTEGER MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXCOE, MAXPHI
      PARAMETER (MAXCOE=16)
      CHARACTER  HEADER*127
      LOGICAL  CORINT, DELTAM, PLANK, ONLYFL, USRANG, USRTAU
      INTEGER  IBCND, NCOEF, NLYR, NUMU, NSTR, NPHI, NTAU
      REAL     ALBEDO, BTEMP, DTAUC(MAXCLY), FBEAM, FISOT,              &
     &     PHI(MAXPHI), PMOM(0:MAXCOE,1:MAXCLY),                        &
     &     SSALBM,SSALB(MAXCLY), TEMIS, TTEMP, WVNMLO, WVNMHI,          &
     &     UMU(MAXUMU), UMU0, UTAU(MAXULV)

      REAL     RFLDIR(MAXULV), RFLDN(MAXULV), FLUP(MAXULV),             &
     &     UAVG(MAXULV), DFDT(MAXULV),                                  &
     &     UU(MAXUMU, MAXULV, MAXPHI)
      REAL WN0, S0CMS(MAXUMU, MAXULV), T0CMS(MAXUMU,MAXULV)
      REAL BBFN,PLKAVG
      DOUBLE PRECISION DISBBF
      LOGICAL MSFLAG
      INTEGER LEV

      REAL     ACCUR, PHI0, TEMPER(0:MAXCLY),                           &
     &     ALBMED(MAXUMU), U0U(MAXUMU, MAXULV), TRNMED(MAXUMU)
      INTEGER MXNLYR

!+---------------------------------------------------------------------+
!      ROUTINES CALLED (IN ORDER):  SLFTST, ZEROAL, CHEKIN, SETDIS,
!                                   PRTINP, ALBTRN, LEPOLY, SURFAC,
!                                   SOLEIG, UPBEAM, UPISOT, TERPEV,
!                                   TERPSO, SETMTX, SOLVE0, FLUXES,
!                                   USRINT, PRAVIN, PRTINT
!+---------------------------------------------------------------------+

!  INDEX CONVENTIONS (FOR ALL DO-LOOPS AND ALL VARIABLE DESCRIPTIONS):

!     IU     :  FOR USER POLAR ANGLES

!  IQ,JQ,KQ  :  FOR COMPUTATIONAL POLAR ANGLES ('QUADRATURE ANGLES')

!   IQ/2     :  FOR HALF THE COMPUTATIONAL POLAR ANGLES (JUST THE ONES
!               IN EITHER 0-90 DEGREES, OR 90-180 DEGREES)

!     J      :  FOR USER AZIMUTHAL ANGLES

!     K,L    :  FOR LEGENDRE EXPANSION COEFFICIENTS OR, ALTERNATIVELY,
!               SUBSCRIPTS OF ASSOCIATED LEGENDRE POLYNOMIALS

!     LU     :  FOR USER LEVELS

!     LC     :  FOR COMPUTATIONAL LAYERS (EACH HAVING A DIFFERENT
!               SINGLE-SCATTER ALBEDO AND/OR PHASE FUNCTION)

!    LEV     :  FOR COMPUTATIONAL LEVELS

!    MAZIM   :  FOR AZIMUTHAL COMPONENTS IN FOURIER COSINE EXPANSION
!               OF INTENSITY AND PHASE FUNCTION

!+---------------------------------------------------------------------+
!               I N T E R N A L    V A R I A B L E S

!   AMB(IQ/2,IQ/2)    FIRST MATRIX FACTOR IN REDUCED EIGENVALUE PROBLEM
!                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN SOLEIG)

!   APB(IQ/2,IQ/2)    SECOND MATRIX FACTOR IN REDUCED EIGENVALUE PROBLEM
!                     OF EQS. SS(12), STWJ(8E)  (USED ONLY IN SOLEIG)

!   ARRAY(IQ,IQ)      SCRATCH MATRIX FOR SOLEIG, UPBEAM AND UPISOT
!                     (SEE EACH SUBROUTINE FOR DEFINITION)

!   B()               RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO
!                     SOLVE0,1;  RETURNS AS SOLUTION VECTOR
!                     VECTOR  L, THE CONSTANTS OF INTEGRATION

!   BDR(IQ/2,0:IQ/2)  BOTTOM-BOUNDARY BIDIRECTIONAL REFLECTIVITY FOR A
!                     GIVEN AZIMUTHAL COMPONENT.  FIRST INDEX ALWAYS
!                     REFERS TO A COMPUTATIONAL ANGLE.  SECOND INDEX:
!                     IF ZERO, REFERS TO INCIDENT BEAM ANGLE UMU0;
!                     IF NON-ZERO, REFERS TO A COMPUTATIONAL ANGLE.

!   BEM(IQ/2)         BOTTOM-BOUNDARY DIRECTIONAL EMISSIVITY AT COMPU-
!                     TATIONAL ANGLES.

!   BPLANK            INTENSITY EMITTED FROM BOTTOM BOUNDARY

!   CBAND()           MATRIX OF LEFT-HAND SIDE OF THE LINEAR SYSTEM
!                     EQ. SC(5), SCALED BY EQ. SC(12);  IN BANDED
!                     FORM REQUIRED BY LINPACK SOLUTION ROUTINES

!   CC(IQ,IQ)         C-SUB-IJ IN EQ. SS(5)

!   CMU(IQ)           COMPUTATIONAL POLAR ANGLES (GAUSSIAN)

!   CWT(IQ)           QUADRATURE WEIGHTS CORRESPONDING TO CMU

!   DELM0             KRONECKER DELTA, DELTA-SUB-M0, WHERE M=MAZIM
!                     IS THE NUMBER OF THE FOURIER COMPONENT IN THE
!                     AZIMUTH COSINE EXPANSION

!   EMU(IU)           BOTTOM-BOUNDARY DIRECTIONAL EMISSIVITY AT USER
!                     ANGLES.

!   EVAL(IQ)          TEMPORARY STORAGE FOR EIGENVALUES OF EQ. SS(12)

!   EVECC(IQ,IQ)      COMPLETE EIGENVECTORS OF SS(7) ON RETURN FROM
!                     SOLEIG; STORED PERMANENTLY IN  GC

!   EXPBEA(LC)        TRANSMISSION OF DIRECT BEAM IN DELTA-M OPTICAL
!                     DEPTH COORDINATES

!   FLYR(LC)          TRUNCATED FRACTION IN DELTA-M METHOD

!   GL(K,LC)          PHASE FUNCTION LEGENDRE POLYNOMIAL EXPANSION
!                     COEFFICIENTS, CALCULATED FROM PMOM BY
!                     INCLUDING SINGLE-SCATTERING ALBEDO, FACTOR
!                     2K+1, AND (IF DELTAM=TRUE) THE DELTA-M
!                     SCALING

!   GC(IQ,IQ,LC)      EIGENVECTORS AT POLAR QUADRATURE ANGLES,
!                     G  IN EQ. SC(1)

!   GU(IU,IQ,LC)      EIGENVECTORS INTERPOLATED TO USER POLAR ANGLES
!                     (G  IN EQS. SC(3) AND S1(8-9), I.E.
!                       G WITHOUT THE L FACTOR)

!   HLPR()            LEGENDRE COEFFICIENTS OF BOTTOM BIDIRECTIONAL
!                     REFLECTIVITY (AFTER INCLUSION OF 2K+1 FACTOR)

!   IPVT(LC*IQ)       INTEGER VECTOR OF PIVOT INDICES FOR LINPACK
!                     ROUTINES

!   KK(IQ,LC)         EIGENVALUES OF COEFF. MATRIX IN EQ. SS(7)

!   KCONV             COUNTER IN AZIMUTH CONVERGENCE TEST

!   LAYRU(LU)         COMPUTATIONAL LAYER IN WHICH USER OUTPUT LEVEL
!                     UTAU(LU) IS LOCATED

!   LL(IQ,LC)         CONSTANTS OF INTEGRATION L IN EQ. SC(1),
!                     OBTAINED BY SOLVING SCALED VERSION OF EQ. SC(5)

!   LYRCUT            TRUE, RADIATION IS ASSUMED ZERO BELOW LAYER
!                     NCUT BECAUSE OF ALMOST COMPLETE ABSORPTION

!   NAZ               NUMBER OF AZIMUTHAL COMPONENTS CONSIDERED

!   NCUT              COMPUTATIONAL LAYER NUMBER IN WHICH ABSORPTION
!                     OPTICAL DEPTH FIRST EXCEEDS ABSCUT

!   OPRIM(LC)         SINGLE SCATTERING ALBEDO AFTER DELTA-M SCALING

!   PASS1             TRUE ON FIRST ENTRY, FALSE THEREAFTER

!   PKAG(0:LC)        INTEGRATED PLANCK FUNCTION FOR INTERNAL EMISSION

!   PSI(IQ)           SUM JUST AFTER SQUARE BRACKET IN  EQ. SD(9)

!   RMU(IU,0:IQ)      BOTTOM-BOUNDARY BIDIRECTIONAL REFLECTIVITY FOR A
!                     GIVEN AZIMUTHAL COMPONENT.  FIRST INDEX ALWAYS
!                     REFERS TO A USER ANGLE.  SECOND INDEX:
!                     IF ZERO, REFERS TO INCIDENT BEAM ANGLE UMU0;
!                     IF NON-ZERO, REFERS TO A COMPUTATIONAL ANGLE.

!   TAUC(0:LC)        CUMULATIVE OPTICAL DEPTH (UN-DELTA-M-SCALED)

!   TAUCPR(0:LC)      CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED IF
!                     DELTAM=TRUE, OTHERWISE EQUAL TO TAUC)

!   TPLANK            INTENSITY EMITTED FROM TOP BOUNDARY

!   UUM(IU,LU)        COMPONENTS OF THE INTENSITY (U-SUPER-M) WHEN
!                     EXPANDED IN FOURIER COSINE SERIES IN AZIMUTH ANGLE

!   U0C(IQ,LU)        AZIMUTHALLY-AVERAGED INTENSITY

!   UTAUPR(LU)        OPTICAL DEPTHS OF USER OUTPUT LEVELS IN DELTA-M
!                     COORDINATES;  EQUAL TO  UTAU(LU) IF NO DELTA-M

!   WK()              SCRATCH ARRAY

!   XR0(LC)           X-SUB-ZERO IN EXPANSION OF THERMAL SOURCE FUNC-
!                     TION PRECEDING EQ. SS(14) (HAS NO MU-DEPENDENCE)

!   XR1(LC)           X-SUB-ONE IN EXPANSION OF THERMAL SOURCE FUNC-
!                     TION;  SEE  EQS. SS(14-16)

!   YLM0(L)           NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
!                     OF SUBSCRIPT L AT THE BEAM ANGLE (NOT SAVED
!                     AS FUNCTION OF SUPERSCIPT M)

!   YLMC(L,IQ)        NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
!                     OF SUBSCRIPT L AT THE COMPUTATIONAL ANGLES
!                     (NOT SAVED AS FUNCTION OF SUPERSCIPT M)

!   YLMU(L,IU)        NORMALIZED ASSOCIATED LEGENDRE POLYNOMIAL
!                     OF SUBSCRIPT L AT THE USER ANGLES
!                     (NOT SAVED AS FUNCTION OF SUPERSCIPT M)

!   Z()               SCRATCH ARRAY USED IN  SOLVE0,1  TO SOLVE A
!                     LINEAR SYSTEM FOR THE CONSTANTS OF INTEGRATION

!   Z0(IQ)            SOLUTION VECTORS Z-SUB-ZERO OF EQ. SS(16)

!   Z0U(IU,LC)        Z-SUB-ZERO IN EQ. SS(16) INTERPOLATED TO USER
!                     ANGLES FROM AN EQUATION DERIVED FROM SS(16)

!   Z1(IQ)            SOLUTION VECTORS Z-SUB-ONE  OF EQ. SS(16)

!   Z1U(IU,LC)        Z-SUB-ONE IN EQ. SS(16) INTERPOLATED TO USER
!                     ANGLES FROM AN EQUATION DERIVED FROM SS(16)

!   ZBEAM(IU,LC)      PARTICULAR SOLUTION FOR BEAM SOURCE

!   ZJ(IQ)            RIGHT-HAND SIDE VECTOR  X-SUB-ZERO IN
!                     EQ. SS(19), ALSO THE SOLUTION VECTOR
!                     Z-SUB-ZERO AFTER SOLVING THAT SYSTEM

!   ZZ(IQ,LC)         PERMANENT STORAGE FOR THE BEAM SOURCE VECTORS ZJ

!   ZPLK0(IQ,LC)      PERMANENT STORAGE FOR THE THERMAL SOURCE
!                     VECTORS  Z0  OBTAINED BY SOLVING  EQ. SS(16)

!   ZPLK1(IQ,LC)      PERMANENT STORAGE FOR THE THERMAL SOURCE
!                     VECTORS  Z1  OBTAINED BY SOLVING  EQ. SS(16)

!+---------------------------------------------------------------------+
!  LOCAL SYMBOLIC DIMENSIONS (HAVE BIG EFFECT ON STORAGE REQUIREMENTS):

!     MXCLY =MAX NO. OF COMPUTATIONAL LAYERS
!     MXULV =MAX NO. OF OUTPUT LEVELS
!     MXCMU =MAX NO. OF COMPUTATION POLAR ANGLES
!     MXUMU =MAX NO. OF OUTPUT POLAR ANGLES
!     MXPHI =MAX NO. OF OUTPUT AZIMUTHAL ANGLES
!MOD
!MOD  INTEGER MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI
!MOD  PARAMETER (MXCLY=35, MXULV=35, MXCMU=4, MXUMU=8,
!MOD $     MXPHI=6, MI=MXCMU/2, MI9M2=9*MI-2,
!MOD $     NNLYRI=MXCMU*MXCLY)
      INCLUDE 'PARAMS.h'

!     /SURFWV/
!       LAMBER  LOGICAL FLAG, .TRUE. FOR LAMBERTIAN SURFACE.
!       TPTEMP  TARGET-PIXEL SURFACE TEMPERATURES [K].
!       TPHDIR  TARGET-PIXEL HEMISPHERE DIRECTIONAL REFLECTANCE AT
!               VIEWING ANGLE.
!       TPBRDF  TARGET-PIXEL BIDIRECTIONAL REFLECTANCE DISTRIBUTION
!               FUNCTION AT VIEWING AND SUN ANGLE.
!       AATEMP  AREA-AVERAGED GROUND SURFACE TEMPERATURES [K].
!       AASALB  AREA-AVERAGED GROUND SURFACE ALBEDO.
!       AADREF  AREA-AVERAGED GROUND SURFACE DIRECTIONAL REFLECTIVITY
!               AT THE SOLAR ZENITH ANGLE.
!       EMU     GROUND DIRECTIONAL EMISSIVITY AT VIEWING ANGLE.
!       BEM     GROUND DIRECTIONAL EMISSIVITY AT QUADRATURE ANGLE.
!       RMU     GROUND BRDF AZIMUTH COMPONENTS AT VIEWING ANGLE
!               AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
!       BDR     GROUND BRDF AZIMUTH COMPONENTS AT QUADRATURE ANGLE
!               AND AT SUN (=0) OR QUADRATURE (>0) ANGLE.
      LOGICAL LAMBER
      REAL TPTEMP,TPHDIR,TPBRDF,AATEMP,AASALB,AADREF,EMU,BEM,RMU,BDR
      COMMON/SURFWV/LAMBER,TPTEMP,TPHDIR,TPBRDF,AATEMP,AASALB,AADREF,   &
     &  EMU(MXUMU),BEM(MI),RMU(1:MXUMU,0:MI,0:MAZ),BDR(1:MI,0:MI,0:MAZ)

      INTEGER MXCLY,MXULV,MXPHI,MI9M2,NNLYRI
      PARAMETER(MXCLY=LAYDIM,MXULV=LAYDIM+1,MXPHI=1,MI9M2=9*MI-2,       &
     &  NNLYRI=MXCMU*MXCLY)

!+---------------------------------------------------------------------+

      LOGICAL LYRCUT, PASS1, PRNT(7)
      INTEGER IPVT(NNLYRI), LAYRU(MXULV)
      REAL    AMB(MI,MI), APB(MI,MI),                                   &
     &     ARRAY(MXCMU,MXCMU), B(NNLYRI),                               &
     &     CBAND(MI9M2,NNLYRI), CC(MXCMU,MXCMU),                        &
     &     CMU(MXCMU), CWT(MXCMU), EVAL(MI),                            &
     &     EVECC(MXCMU,MXCMU), EXPBEA(0:MXCLY), FLYR(MXCLY),            &
     &     FLDN(MAXULV), FLDIR(MAXULV), GL(0:MXCMU,MXCLY),              &
     &     GC(MXCMU,MXCMU,MXCLY), GU(MXUMU,MXCMU,MXCLY),                &
     &     KK(MXCMU,MXCLY), LL(MXCMU,MXCLY),                            &
     &     OPRIM(MXCLY), PHASA(MXCLY), PHASE(MXCLY),                    &
     &     PHASM(MXCLY), PHIRAD(MXPHI), PKAG(0:MXCLY),                  &
     &     PSI(MXCMU), TAUC(0:MXCLY),                                   &
     &     TAUCPR(0:MXCLY), U0C(MXCMU,MXULV), UTAUPR(MXULV),            &
     &     UUM(MXUMU,MXULV), WK(MXCMU), XR0(MXCLY),                     &
     &     XR1(MXCLY), YLM0(0:MXCMU), YLMC(0:MXCMU,MXCMU),              &
     &     YLMU(0:MXCMU,MXUMU), Z(NNLYRI), Z0(MXCMU),                   &
     &     Z0U(MXUMU,MXCLY), Z1(MXCMU), Z1U(MXUMU,MXCLY),               &
     &     ZBEAM(MXUMU,MXCLY), ZJ(MXCMU), ZZ(MXCMU,MXCLY),              &
     &     ZPLK0(MXCMU,MXCLY), ZPLK1(MXCMU,MXCLY)
      REAL SAVEEXPBEA(0:MXCLY), SAVELL(MXCMU,MXCLY)
      REAL EPSIL, DUM, DELM0, SGN, DELTAT,AZERR, COSPHI, AZTERM
      REAL TPLANK, BPLANK
      REAL RATIO
      INTEGER II   !DRF ADDITION
      INTEGER NCUT, NN, KCONV, NAZ, MAZIM, L, IQ, IU, LC,               &
     &     NCOL, LU, J, JK,IK
      REAL DUMMY0(MXCMU,MXCLY), DUMMY1(MXCMU,MXCLY), DUMBEM(MI),        &
     &     USAVE(MXUMU),                                                &
     &     CBANDS(MI9M2,NNLYRI), CBANDT(MI9M2,NNLYRI),                  &
     &     BEAMMS(MXUMU,MXCLY),Z0UMS(MXUMU,MXCLY),Z1UMS(MXUMU,MXCLY),   &
     &     FDNSRT,FDNTRT

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     LOCAL ARRAY:
!       UMU0X    ARRAY OF DIMENSION 1 CONTAINING NEGATIVE OF UMU0
      REAL UMU0X(1)
!     DATA  PASS1 / .TRUE. /

      MXNLYR=MXCMU*MAXCLY
!                  INITIALIZE UNUSED U0U, ALBMED, TRNMED, HEADER, TEMPER
!MOD  PHI0 =0.
!MOD  ACCUR=0.
      PASS1=.FALSE.
      EPSIL=20.*RRIGHT
      IF (PASS1)  THEN
         PI=2. * ASIN(1.0)
!MOD     EPSIL=10.*R1MACH(4)
!MOD     RPD=PI / 180.0
!                                      INSERT INPUT VALUES FOR SELF-TEST
!                                   NOTE: SELF-TEST MUST NOT USE IBCND=1

      CALL  SLFTST(ACCUR, ALBEDO, BPLANK, BTEMP, DELTAM, DTAUC(1),      &
     &        FBEAM, FISOT, IBCND, CORINT, LAMBER, NLYR, PLANK,         &
     &        NPHI, NUMU, NSTR, NTAU, ONLYFL, PHI(1), PHI0,             &
     &        PKAG, PMOM, PRNT, SSALB(1), TEMIS, TEMPER,                &
     &        TPLANK, TTEMP, UMU(1), USRANG, USRTAU,                    &
     &        UTAU(1), UMU0, WVNMHI, WVNMLO, .FALSE.,                   &
     &        DUM, DUM, DUM, DUM)

      END IF

 10   CONTINUE

!-----------------------------------------------------------------------
!-------------------FOR MODTRAN (THE CODE BELOW)
!            **  INITIALIZE FOR MULTIPLE SCATTERING SOURCE FUNCTIONS USE
      IF (.NOT.PASS1) THEN
         DO IU = 1,NUMU
            USAVE(IU) = UMU(IU)
         ENDDO
         CALL ZEROIT(BEAMMS,MXUMU*MAXCLY)
         CALL ZEROIT(Z0UMS,MXUMU*MAXCLY)
         CALL ZEROIT(Z1UMS,MXUMU*MAXCLY)
      ENDIF
!-------------------FOR MODTRAN (THE CODE ABOVE)
!-----------------------------------------------------------------------

!                          CALCULATE CUMULATIVE OPTICAL DEPTH AND DITHER
!                             SINGLE-SCATTER ALBEDO TO IMPROVE NUMERICAL
!                              BEHAVIOR OF EIGENVALUE/VECTOR COMPUTATION

      CALL  ZEROIT(TAUC, MAXCLY+1)
      SSALBM=1.- 8*EPSIL
      DO LC=1,NLYR
          IF(SSALB(LC).GT.SSALBM)SSALB(LC)=SSALBM
          TAUC(LC)=TAUC(LC-1)+DTAUC(LC)
      ENDDO
!                                   CHECK INPUT DIMENSIONS AND VARIABLES
      CALL  CHEKIN(ACCUR, ALBEDO, BTEMP, DTAUC, FBEAM, FISOT,           &
     &     IBCND, LAMBER, MAXCLY, MAXCMU, MAXPHI,                       &
     &     MAXULV, MAXUMU, MXCLY, MXCMU, MXPHI, MXULV, MXUMU,           &
     &     NLYR, NPHI, NSTR, NTAU, NUMU, ONLYFL, PHI, PHI0,             &
     &     PLANK, PMOM, SSALB, TAUC, TEMIS, TEMPER, TTEMP,              &
     &     UMU, UMU0, USRANG, USRTAU, UTAU, WVNMHI, WVNMLO)
!                          ZERO SOME ARRAYS (NOT STRICTLY NECESSARY, BUT
!                      OTHERWISE UNUSED PARTS OF ARRAYS COLLECT GARBAGE)

      CALL  ZEROAL(MXCLY, EXPBEA(1), FLYR, OPRIM, PHASA, PHASE, PHASM,  &
     &     TAUCPR(1), XR0, XR1,                                         &
     &     MXCMU, CMU, CWT, PSI, WK, Z0, Z1, ZJ,                        &
     &     MXCMU+1, YLM0,                                               &
     &     MXCMU**2, ARRAY, CC, EVECC,                                  &
     &     (MXCMU+1)*MXCLY, GL,                                         &
     &     (MXCMU+1)*MXCMU, YLMC,                                       &
     &     (MXCMU+1)*MXUMU, YLMU,                                       &
     &     MXCMU*MXCLY, KK, LL, ZZ, ZPLK0, ZPLK1,                       &
     &     MXCMU**2*MXCLY, GC,                                          &
     &     MXULV, LAYRU, UTAUPR,                                        &
     &     MXUMU*MXCMU*MXCLY, GU,                                       &
     &     MXUMU*MXCLY, Z0U, Z1U, ZBEAM,                                &
     &     MI, EVAL,                                                    &
     &     MI**2, AMB, APB,                                             &
     &     NNLYRI, IPVT, Z)

!                                       PERFORM VARIOUS SETUP OPERATIONS

!!BEGIN DRF MODIFICATION

      IF (FBEAM.GT.0.0)  THEN
         NN = NSTR/2
         CALL  QGAUSN(NN, CMU, CWT)
         DO II=1,NN
            IF (ABS(UMU0-CMU(II))/UMU0 .LT. 1.E-4) THEN
               UMU0=UMU0*1.001  
               !WRITE(*,*) "CHANGING BEAM ANGLE (DRF)"
            ENDIF
         END DO
         !RE-ZERO CMU AND CWT
         DO II=1,MXCMU
            CMU(II) = 0.
            CWT(II) = 0.
         END DO
      ENDIF

!!END DRF MODIFICATION

      CALL  SETDIS(CMU, CWT, DELTAM, DTAUC, EXPBEA, FBEAM, FLYR,        &
     &     GL, IBCND, LAYRU, LYRCUT, MAXCOE,                            &
     &     MAXUMU, MXCMU, NCUT, NLYR, NTAU, NN, NSTR,                   &
     &     PLANK, NUMU, ONLYFL, OPRIM, PMOM, SSALB, TAUC,               &
     &     TAUCPR, UTAU, UTAUPR, UMU, UMU0, USRTAU, USRANG)

!                                                PRINT INPUT INFORMATION
      IF (PRNT(1))                                                      &
     &     CALL PRTINP(HEADER, NLYR, DTAUC, SSALB, PMOM, TEMPER,        &
     &     WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,                 &
     &     NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT,                  &
     &     LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,                         &
     &     DELTAM, PLANK, ONLYFL, ACCUR, FLYR, LYRCUT,                  &
     &     OPRIM, TAUC, TAUCPR, MAXCOE, PRNT(7))

!                             HANDLE SPECIAL CASE FOR GETTING ALBEDO AND
!                  TRANSMISSIVITY OF MEDIUM FOR MANY BEAM ANGLES AT ONCE

      IF (IBCND.EQ.1)  THEN
         CALL  ALBTRN(ALBEDO, AMB, APB, ARRAY, B, CBAND, CC,            &
     &        CMU, CWT, EVAL, EVECC, GL, GC, GU, IPVT, KK,              &
     &        LL, NLYR, NN, NSTR, NUMU, PRNT, TAUCPR, UMU,              &
     &        U0U, WK, YLMC, YLMU, Z, MI9M2, MAXULV,                    &
     &        MAXUMU, NNLYRI, WN0, ALBMED, TRNMED)
         RETURN
      END IF

!-------------------THE CODE BELOW IS FOR MODTRAN
!                                   ** CALCULATE PLANCK FUNCTIONS
      IF(.NOT.PLANK)THEN
         BPLANK = 0.0
         TPLANK = 0.0
         CALL ZEROIT(PKAG,MAXCLY+1)
!                      **  USE DIFFERENT PLANCK FUNCTIONS PLKAVG OR BBFN
      ELSEIF(PASS1)THEN
         TPLANK = TEMIS*PLKAVG(WVNMLO,WVNMHI,TTEMP)
         BPLANK = PLKAVG(WVNMLO,WVNMHI,BTEMP)
         DO LEV = 0,NLYR
            PKAG(LEV) = PLKAVG(WVNMLO,WVNMHI,TEMPER(LEV))
         ENDDO
      ELSE
         TPLANK = TEMIS*BBFN(TTEMP,WN0)
         BPLANK = BBFN(BTEMP,WN0)
         DO LEV = 0,NLYR
            PKAG(LEV) = REAL(DISBBF(DBLE(TEMPER(LEV)),DBLE(WN0)))
         ENDDO
      ENDIF
!------------------THE CODE ABOVE FOR MODTRAN

! ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
!           (EQ STWJ 5)

      KCONV=0
      NAZ=NSTR-1
!                                               AZIMUTH-INDEPENDENT CASE

      IF (FBEAM.EQ.0.0 .OR. (1.-UMU0).LT.1.E-5 .OR. ONLYFL .OR.         &
     &     (NUMU.EQ.1.AND.(1.-UMU(1)).LT.1.E-5))                        &
     &     NAZ=0
      CALL  ZEROIT(UU, MAXUMU*MAXULV*MAXPHI)
      IF(FBEAM.GT.0.)UMU0X(1)=-UMU0
      DO MAZIM=0, NAZ

         IF (MAZIM.EQ.0)  DELM0=1.0
         IF (MAZIM.GT.0)  DELM0=0.0

!                                     GET NORMALIZED ASSOCIATED LEGENDRE
!                             POLYNOMIALS FOR INCIDENT BEAM ANGLE COSINE

         IF(FBEAM.GT.0.)CALL LEPOLY(1,MAZIM,MXCMU,NSTR-1,UMU0X,YLM0)

!                         GET NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
!                         FOR COMPUTATIONAL AND USER POLAR ANGLE COSINES

         IF (.NOT.ONLYFL .AND. USRANG)                                  &
     &        CALL  LEPOLY(NUMU, MAZIM, MXCMU, NSTR-1, UMU, YLMU)
         CALL  LEPOLY(NN,   MAZIM, MXCMU, NSTR-1, CMU, YLMC)

!-------------------BEGIN MODTRAN
         IF(.NOT.PASS1)THEN
            DO IU = 1,NUMU
               UMU(IU) = USAVE(IU)
            ENDDO
            CALL LEPOLY(NUMU,MAZIM,MXCMU,NSTR-1,UMU,YLMU)
         ENDIF
!-------------------END MODTRAN

!                                EVALUATE NORMALIZED ASSOCIATED LEGENDRE
!                               POLYNOMIALS WITH NEGATIVE CMU FROM THOSE
!                              WITH POSITIVE CMU; DAVE/ARMSTRONG EQ.(15)
         SGN =- 1.0
         DO L=MAZIM, NSTR-1
            SGN=- SGN
            DO IQ=NN+1, NSTR
               YLMC(L,IQ)=SGN * YLMC(L,IQ-NN)
            ENDDO
         ENDDO

!            SPECIFY USERS BOTTOM REFLECTIVITY AND EMISSIVITY PROPERTIES

!mod     IF (.NOT.LYRCUT)
!mod $        CALL  SURFAC(ALBEDO, CMU, FBEAM, LAMBER, MAZIM,
!mod $        MI, MXUMU, NN, NUMU, ONLYFL, PI,
!mod $        UMU, UMU0, USRANG, BDR, EMU, BEM, RMU)

! ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============

         DO LC=1, NCUT

!                           SOLVE EIGENFUNCTION PROBLEM IN EQ. STWJ(8B);
!                                    RETURN EIGENVALUES AND EIGENVECTORS

            CALL  SOLEIG(AMB, APB, ARRAY, CMU, CWT, GL(0,LC),MI,MAZIM,  &
     &           MXCMU, NN, NSTR, WK, YLMC, CC, EVECC, EVAL,            &
     &           KK(1,LC), GC(1,1,LC))

!                                      CALCULATE PARTICULAR SOLUTIONS OF
!                                     EQ.SS(18) FOR INCIDENT BEAM SOURCE
            IF (FBEAM.GT.0.0)                                           &
     &           CALL  UPBEAM(ARRAY, CC, CMU, DELM0, FBEAM, GL(0,LC),   &
     &           IPVT, MAZIM, MXCMU, NN, NSTR, PI, UMU0, WK,            &
     &           YLM0, YLMC, ZJ, ZZ(1,LC))

!                                      CALCULATE PARTICULAR SOLUTIONS OF
!                                  EQ.SS(15) FOR THERMAL EMISSION SOURCE

            IF (PLANK .AND. MAZIM.EQ.0)  THEN
               DELTAT=TAUCPR(LC) - TAUCPR(LC-1)
               XR1(LC)=0.0
               IF (DELTAT.GT.0.0) XR1(LC)=(PKAG(LC) - PKAG(LC-1))       &
     &              / DELTAT
               XR0(LC)=PKAG(LC-1) - XR1(LC) * TAUCPR(LC-1)
               CALL UPISOT(ARRAY, CC, CMU, IPVT, MXCMU, NN, NSTR,       &
     &              OPRIM(LC), WK, XR0(LC), XR1(LC), Z0, Z1,            &
     &              ZPLK0(1,LC), ZPLK1(1,LC))
            END IF
!                                INTERPOLATE EIGENVECTORS TO USER ANGLES
            IF ((.NOT.ONLYFL .AND. USRANG).OR.MSFLAG)  THEN
               CALL  TERPEV(CWT, EVECC, GL(0,LC), GU(1,1,LC), MAZIM,    &
     &              MXCMU, MXUMU, NN, NSTR, NUMU, WK, YLMC, YLMU)

!                                INTERPOLATE SOURCE TERMS TO USER ANGLES

               CALL  TERPSO(CWT, DELM0, FBEAM, GL(0,LC), MAZIM, MXCMU,  &
     &              PLANK, NUMU, NSTR, OPRIM(LC), PI,                   &
     &              YLM0, YLMC, YLMU, PSI, XR0(LC), XR1(LC),            &
     &              Z0, ZJ, ZBEAM(1,LC), Z0U(1,LC), Z1U(1,LC),          &
     &              Z0UMS(1,LC),Z1UMS(1,LC),BEAMMS(1,LC))
            END IF
         ENDDO

! ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============

!                          SET COEFFICIENT MATRIX OF EQUATIONS COMBINING
!                                BOUNDARY AND LAYER INTERFACE CONDITIONS

         CALL  SETMTX(BDR(1,0,MAZIM), CBAND, CMU, CWT, DELM0, GC, KK,   &
     &        LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT, NNLYRI,     &
     &        NN, NSTR, TAUCPR, WK)

!-------------------FOR MODTRAN BELOW
         DO IK = 1,MXNLYR
            DO JK = 1,MI9M2
               CBANDS(JK,IK) = CBAND(JK,IK)
               CBANDT(JK,IK) = CBAND(JK,IK)
            ENDDO
         ENDDO

!-------------------FOR MODTRAN ABOVE

!                      SOLVE FOR CONSTANTS OF INTEGRATION IN HOMOGENEOUS
!                                 SOLUTION (GENERAL BOUNDARY CONDITIONS)

         CALL  SOLVE0(B, BDR(1,0,MAZIM), BEM, BPLANK, CBAND, CMU, CWT,  &
     &        EXPBEA, FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM, MI,&
     &        MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, MXNLYR, PI,           &
     &        TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1)

!                                     COMPUTE UPWARD AND DOWNWARD FLUXES
         IF (MAZIM.EQ.0)                                                &
     &        CALL FLUXES(CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,   &
     &        MAXULV, MXCMU, MXULV, NCUT, NN, NSTR, NTAU,               &
     &        PI, PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,              &
     &        XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP,                   &
     &        FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C)

!-------------------THE STUFF BELOW IS INSERTED FOR MODTRAN
!-----------------------------------------------------------------------
!            **  COMPUTE MULTIPLE SCATTERING SOURCE FUNCTIONS SEPARATELY

         IF (.NOT.(ONLYFL.AND.MSFLAG)) PASS1=.TRUE.
         IF((.NOT.PASS1)) THEN

            CALL ZEROIT(DUMMY0,MXNLYR)
!                            **  FIRST CALL FOR SOLAR MS SOURCE FUNCTION
            IF(FBEAM.GT.0.)THEN

               CALL ZEROIT(DUMMY1,MXNLYR)
               CALL ZEROIT(DUMBEM,MI)
               CALL SOLVE0(B,BDR(1,0,MAZIM),DUMBEM,0.0,CBANDS,CMU,CWT,  &
     &              EXPBEA,FBEAM,FISOT,IPVT,LAMBER,LL,LYRCUT,MAZIM,MI,  &
     &              MI9M2,MXCMU,NCOL,NCUT,NN,NSTR,MXNLYR,PI,0.0,        &
     &              TAUCPR,UMU0,Z,ZZ,DUMMY0,DUMMY1)
               CALL MSSOLR(CMU,CWT,FBEAM,GC,GU,KK,LAYRU,LL,LYRCUT,      &
     &              MAXUMU,MXCMU,MXUMU,NCUT,NN,NSTR,NTAU,NUMU,PI,       &
     &              TAUCPR,UMU0,UTAU,UTAUPR,BEAMMS,ZZ,FDNSRT,           &
     &              S0CMS)
            ENDIF
!                         **  SECOND CALL FOR THERMAL MS SOURCE FUNCTION

!           COMPUTE THERMAL MS SOURCE FUNCTION

            IF(PLANK)THEN
               CALL ZEROIT(DUMMY0, MXNLYR)

               CALL SAV1D0(EXPBEA,SAVEEXPBEA,MXCLY)
               CALL SAV1D1(LL,SAVELL,MXCMU*MXCLY)

               DO JK = 0,MXCLY
                  EXPBEA(JK) = 0.0
               ENDDO

               CALL SOLVE0(B,BDR(1,0,MAZIM),BEM,BPLANK,CBANDT,CMU,CWT,  &
     &              EXPBEA,FBEAM,FISOT,IPVT,LAMBER,LL,LYRCUT,MAZIM,MI,  &
     &              MI9M2,MXCMU,NCOL,NCUT,NN,NSTR,MXNLYR,PI,            &
     &              TPLANK,TAUCPR,UMU0,Z,DUMMY0,ZPLK0,ZPLK1)

               CALL MSTHML(CMU,CWT,GC,GU,KK,LAYRU,LL,LYRCUT,MAXUMU,     &
     &              MXCMU,MXUMU,NCUT,NN,NSTR,NTAU,NUMU,OPRIM,PI,        &
     &              TAUCPR,UTAUPR,XR0,XR1,Z0UMS,Z1UMS,ZPLK0,            &
     &              ZPLK1,FDNTRT,T0CMS)

               CALL SAV1D1(SAVELL,LL,MXCMU*MXCLY)
               CALL SAV1D0(SAVEEXPBEA,EXPBEA,MXCLY)

               CALL ZEROIT(B,MXNLYR)
               CALL ZEROIT(DUMMY0,MXNLYR)
             ENDIF
         ENDIF
         IF (.NOT.(ONLYFL.AND.MSFLAG)) PASS1=.FALSE.
!-----------------------------------------------------------------------
!-------------------THE ABOVE INSERTED FOR MODTRAN

         IF (ONLYFL)  THEN
            IF (MAXUMU.GE.NSTR)  THEN

!             SAVE AZIMUTHALLY AVERAGED INTENSITIES AT QUADRATURE ANGLES

               DO LU=1, NTAU
                  DO IQ=1, NSTR
                     U0U(IQ,LU)=U0C(IQ,LU)
                  ENDDO
               ENDDO

            ELSE
               CALL  ZEROIT(U0U, MAXUMU*MAXULV)
            END IF
            GO TO 300
         END IF

         CALL  ZEROIT(UUM, MXUMU*MAXULV)
         IF (USRANG)  THEN

!                  COMPUTE AZIMUTHAL INTENSITY COMPONENTS AT USER ANGLES

            CALL  USRINT(BPLANK, CMU, CWT, DELM0, EMU, EXPBEA, FBEAM,   &
     &           FISOT, GC, GU, KK, LAMBER, LAYRU, LL, LYRCUT,          &
     &           MAZIM, MXCMU, MXUMU, NCUT, NLYR, NN, NSTR,             &
     &           PLANK, NUMU, NTAU, PI, RMU(1,0,MAZIM), TAUCPR, TPLANK, &
     &           UMU, UMU0, UTAUPR, WK, ZBEAM, Z0U, Z1U, ZZ,            &
     &           ZPLK0, ZPLK1, UUM)
         ELSE
!            COMPUTE AZIMUTHAL INTENSITY COMPONENTS AT QUADRATURE ANGLES

            CALL  CMPINT(FBEAM, GC, KK, LAYRU, LL, LYRCUT, MAZIM, MXCMU,&
     &           MXUMU, NCUT, NN, NSTR, PLANK, NTAU, TAUCPR,            &
     &           UMU0, UTAUPR, ZZ, ZPLK0, ZPLK1, UUM)
         END IF

         IF (MAZIM.EQ.0)  THEN
!                                  SAVE AZIMUTHALLY AVERAGED INTENSITIES
            DO LU=1, NTAU
               DO IU=1, NUMU
                  U0U(IU,LU)=UUM(IU,LU)
                  DO J=1, NPHI
                     UU(IU,LU,J)=UUM(IU,LU)
                  ENDDO
               ENDDO
            ENDDO
!                  PRINT AZIMUTHALLY AVERAGED INTENSITIES AT USER ANGLES

            IF (PRNT(4))                                                &
     &           CALL  PRAVIN(UMU, NUMU, MAXUMU, UTAU, NTAU, U0U)
            IF (NAZ.GT.0)  THEN
               CALL  ZEROIT(PHIRAD, MXPHI)
               DO J=1, NPHI
                  PHIRAD(J)= (PHI(J) - PHI0) / DEG
               ENDDO
            END IF

         ELSE
!                               INCREMENT INTENSITY BY CURRENT AZIMUTHAL
!                            COMPONENT (FOURIER COSINE SERIES); EQ.SD(2)
            AZERR=0.0
            DO J=1, NPHI
               COSPHI=COS(MAZIM * PHIRAD(J))
               DO LU=1, NTAU
                  DO IU=1, NUMU
                     AZTERM=UUM(IU,LU) * COSPHI
                     UU(IU,LU,J)=UU(IU,LU,J) + AZTERM
                     AZERR=AMAX1(RATIO(ABS(AZTERM), ABS(UU(IU,LU,J))),  &
     &                    AZERR)
                  ENDDO
               ENDDO
            ENDDO
            IF (AZERR.LE.ACCUR)  KCONV=KCONV + 1
            IF (KCONV.GE.2)      GO TO 300
         END IF

      ENDDO

! ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============

 300  CONTINUE
!                                  NAKAJIMA/TANAKA INTENSITY CORRECTIONS

      IF (CORINT)                                                       &
     &     CALL  INTCOR(EPSIL, FBEAM, FLYR, LAYRU, LYRCUT, MAXCOE,      &
     &     MAXULV, MAXUMU, NCOEF, NCUT, NPHI, NSTR,                     &
     &     NTAU, NUMU, OPRIM, PHASA, PHASE, PHASM,                      &
     &     PHIRAD, PI, PMOM, SSALB, TAUC, TAUCPR,                       &
     &     UMU, UMU0, UTAU, UTAUPR, UU)

!                                                      PRINT INTENSITIES
      IF (PRNT(5) .AND. .NOT.ONLYFL)                                    &
     &     CALL  PRTINT(UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI,           &
     &     MAXULV, MAXUMU)

!        COMPARE TEST CASE RESULTS WITH CORRECT ANSWERS AND ABORT IF BAD

      IF (PASS1)  THEN
         CALL  SLFTST(ACCUR, ALBEDO, BPLANK, BTEMP, DELTAM, DTAUC(1),   &
     &        FBEAM, FISOT, IBCND, CORINT, LAMBER, NLYR, PLANK,         &
     &        NPHI, NUMU, NSTR, NTAU, ONLYFL, PHI(1), PHI0,             &
     &        PKAG, PMOM, PRNT, SSALB(1), TEMIS, TEMPER,                &
     &        TPLANK, TTEMP, UMU(1), USRANG, USRTAU,                    &
     &        UTAU(1), UMU0, WVNMHI, WVNMLO, .TRUE.,                    &
     &        FLUP(1), RFLDIR(1), RFLDN(1), UU(1,1,1))
         PASS1=.FALSE.
         GO TO 10
      END IF
      RETURN
      END

      SUBROUTINE SAV1D0(X,XSAVE,LENGTH)
      INTEGER I, LENGTH
      REAL X(0:LENGTH),XSAVE(0:LENGTH)

!     THE STARTING  INDEX IS 0; THIS MUST BE NOTED
!     SAVE WHAT IS IN X IN XSAVE

      DO I=0, LENGTH
         XSAVE(I)=X(I)
      ENDDO
      END

      SUBROUTINE SAV1D1(X,XSAVE,LENGTH)
      INTEGER I, LENGTH
      REAL X(1:LENGTH),XSAVE(1:LENGTH)

!     THE STARTING  INDEX IS 1; THIS MUST BE NOTED
!     SAVE WHAT IS IN X IN XSAVE

      DO I=1, LENGTH
         XSAVE(I)=X(I)
      ENDDO
      END
