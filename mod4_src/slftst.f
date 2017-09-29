      SUBROUTINE  SLFTST(ACCUR, ALBEDO, BPLANK, BTEMP, DELTAM, DTAUC,   &
     &                    FBEAM, FISOT, IBCND, CORINT, LAMBER, NLYR,    &
     &                    PLANK, NPHI, NUMU, NSTR, NTAU, ONLYFL, PHI,   &
     &                    PHI0, PKAG, PMOM, PRNT, SSALB, TEMIS, TEMPER, &
     &                    TPLANK, TTEMP, UMU, USRANG, USRTAU, UTAU,     &
     &                    UMU0, WVNMHI, WVNMLO, COMPAR,                 &
     &                    FLUP, RFLDIR, RFLDN, UU)

!       If  COMPAR=FALSE, save user input values that would otherwise
!       be destroyed and replace them with input values for self-test.
!       If  COMPAR=TRUE, compare self-test case results with correct
!       answers and restore user input values if test is passed.

!       (See file 'DisORT.Doc' for variable definitions.)

!                 I N T E R N A L    V A R I A B L E S:

!         ACC     Relative accuracy required for passing self-test
!         ERRORn  Relative errors in 'DisORT' output variables
!         OK      Logical variable for determining failure of self-test
!         All variables ending in 'S':  Temporary 'S'torage for input
!+---------------------------------------------------------------------+
      REAL     PMOM(0:*), TEMPER(0:*), PKAG(0:*)
      LOGICAL  COMPAR, CORINT, DELTAM, LAMBER, PLANK, OK, ONLYFL,       &
     &         PRNT(*), USRANG, USRTAU
      REAL     PMOMS(0:4), TEMPES (0:1), PKAGS(0:1)
      LOGICAL  CORINS, DELTAS, LAMBES, PLANKS, ONLYFS, PRNTS (7),       &
     &         USRANS, USRTAS, TSTBAD
      SAVE
      DATA     ACC / 1.E-4 /

      IF  (.NOT.COMPAR)  THEN
!                                                 Save user input values
         NLYRS =NLYR
         DTAUCS=DTAUC
         SSALBS=SSALB
         DO 10 N=0, 4
10         PMOMS(N)=PMOM(N)
         NSTRS =NSTR
         USRANS=USRANG
         NUMUS =NUMU
         UMUS  =UMU
         USRTAS=USRTAU
         NTAUS =NTAU
         UTAUS =UTAU
         NPHIS =NPHI
         PHIS  =PHI
         IBCNDS=IBCND
         CORINS=CORINT
         FBEAMS=FBEAM
         UMU0S =UMU0
         PHI0S =PHI0
         FISOTS=FISOT
         LAMBES=LAMBER
         ALBEDS=ALBEDO
         DELTAS=DELTAM
         ONLYFS=ONLYFL
         ACCURS=ACCUR
         PLANKS=PLANK
         WVNMLS=WVNMLO
         WVNMHS=WVNMHI
         BTEMPS=BTEMP
         TTEMPS=TTEMP
         TEMISS=TEMIS
         BPLNKS=BPLANK
         TPLNKS=TPLANK
         TEMPES(0)=TEMPER(0)
         TEMPES(1)=TEMPER(1)
         PKAGS(0)=PKAG(0)
         PKAGS(1)=PKAG(1)
         DO 20 I=1, 7
20         PRNTS(I)=PRNT(I)
!                                         Set input values for self-test
         NLYR =1
         DTAUC=1.0
         SSALB=0.9
!                                                         Haze-L moments
         PMOM(0)=1.0
         PMOM(1)=0.8042
         PMOM(2)=0.646094
         PMOM(3)=0.481851
         PMOM(4)=0.359056
         NSTR  =4
         USRANG=.TRUE.
         NUMU  =1
         UMU   =0.5
         USRTAU=.TRUE.
         NTAU  =1
         UTAU  =0.5
         NPHI  =1
         PHI   =90.0
         IBCND =0
         CORINT=.FALSE.
         FBEAM =3.1415927
         UMU0  =0.866
         PHI0  =0.0
         FISOT =1.0
         LAMBER=.TRUE.
         ALBEDO=0.7
         DELTAM=.TRUE.
         ONLYFL=.FALSE.
         ACCUR =1.E-4
         PLANK =.TRUE.
         WVNMLO=0.0
         WVNMHI=50000.
         BTEMP =300.0
         TTEMP =100.0
         TEMIS =0.8
         TEMPER(0)=210.0
         TEMPER(1)=200.0
         DO 30 I=1, 7
30         PRNT(I)=.FALSE.
         TPLANK=TEMIS * PLKAVG(WVNMLO, WVNMHI, TTEMP)
         BPLANK=        PLKAVG(WVNMLO, WVNMHI, BTEMP)
         PKAG(0)=PLKAVG(WVNMLO, WVNMHI, TEMPER(0))
         PKAG(1)=PLKAVG(WVNMLO, WVNMHI, TEMPER(1))

      ELSE

!        Compare test case results with correct answers and abort if bad

         OK=.TRUE.
         ERROR1=(UU  - 47.86005) / 47.86005
         ERROR2=(RFLDIR - 1.527286) / 1.527286
         ERROR3=(RFLDN - 28.37223) / 28.37223
         ERROR4=(FLUP   - 152.5853) / 152.5853
         IF(ABS(ERROR1).GT.ACC) OK=TSTBAD('UU',     ERROR1)
         IF(ABS(ERROR2).GT.ACC) OK=TSTBAD('RFLDIR', ERROR2)
         IF(ABS(ERROR3).GT.ACC) OK=TSTBAD('RFLDN',  ERROR3)
         IF(ABS(ERROR4).GT.ACC) OK=TSTBAD('FLUP',   ERROR4)

         IF(.NOT. OK)                                                   &
     &       CALL ERRMSG('DISORT--SELF-TEST FAILED', .TRUE.)

!                                              Restore user input values
         NLYR =NLYRS
         DTAUC=DTAUCS
         SSALB=SSALBS
         DO 40 N=0, 4
40         PMOM(N)=PMOMS(N)
         NSTR  =NSTRS
         USRANG=USRANS
         NUMU  =NUMUS
         UMU   =UMUS
         USRTAU=USRTAS
         NTAU  =NTAUS
         UTAU  =UTAUS
         NPHI  =NPHIS
         PHI   =PHIS
         IBCND =IBCNDS
         CORINT=CORINS
         FBEAM =FBEAMS
         UMU0  =UMU0S
         PHI0  =PHI0S
         FISOT =FISOTS
         LAMBER=LAMBES
         ALBEDO=ALBEDS
         DELTAM=DELTAS
         ONLYFL=ONLYFS
         ACCUR =ACCURS
         PLANK =PLANKS
         WVNMLO=WVNMLS
         WVNMHI=WVNMHS
         BTEMP =BTEMPS
         TTEMP =TTEMPS
         TEMIS =TEMISS
         BPLANK=BPLNKS
         TPLANK=TPLNKS
         TEMPER(0)=TEMPES(0)
         TEMPER(1)=TEMPES(1)
         PKAG(0)=PKAGS(0)
         PKAG(1)=PKAGS(1)
         DO 50 I=1, 7
50         PRNT(I)=PRNTS(I)
      END IF

      RETURN
      END
