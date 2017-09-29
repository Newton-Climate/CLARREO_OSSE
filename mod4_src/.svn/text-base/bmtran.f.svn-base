      REAL FUNCTION BMTRAN(DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2)

!     THIS FUNCTION RETURNS PATH TRANSMITTANCE
!     AVERAGED OVER A SPECTRAL INTERVAL.

!     OUTPUTS:
!       BMTRAN  =   PATH TRANSMITTANCE [DIMENSIONLESS]

!     INPUTS:
!                    _
!       DEPTH   =   >_  SU/D,     THE ABSORPTION COEFFICIENT (S/D)
!                                 AND COLUMN DENSITY (U) PRODUCT,
!                                 SUMMED OVER THE PATH [DIMENSIONLESS]

!       ODBAR   =   <1/D>,        THE RECIPROCAL OF THE LINE
!                                 SPACING, (D), PATH AVERAGED [CM]

!       ADBAR   =   <ALFDOP/D>,   THE DOPPLER LINE WIDTH (ALFDOP)
!                                 OVER THE AVERAGE LINE SPACING (D),
!                                 PATH AVERAGED [DIMENSIONLESS]

!       ACBAR   =   <ALFCOL/D>,   THE COLLISION (LORENTZ) LINE WIDTH
!                                 (ALFCOL) OVER THE AVERAGE LINE SPACING
!                                 (D), PATH AVERAGED [DIMENSIONLESS]

!                          2  2
!       ACBAR2  =   <ALFCOL /D >, THE COLLISION (LORENTZ) LINE WIDTH
!                                 (ALFCOL) SQUARED OVER THE AVERAGE
!                                 LINE SPACING (D) SQUARED, PATH
!                                 AVERAGED [DIMENSIONLESS] (OZONE ONLY)

!     ARGUMENTS:
      REAL DEPTH,ODBAR,ADBAR,ACBAR,ACBAR2

!     INCLUDE PARAMETERS:
!       ONE6TH   ONE-SIXTH.
      REAL ONE6TH
      PARAMETER(ONE6TH=1./6.)
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'BMHEAD.h'
!     LOCAL VARIABLES:
!       RLINES   EFFECTIVE NUMBER OF LINES IN THE INTERVAL.
!       ELINES   EFFECTIVE NUMBER OF LINES FROM LINE CENTER
!                TO THE NEAREST INTERVAL EDGE.
!       SUGCPI   PRODUCT OF LINE STRENGTH (S), COLUMN AMOUNT (U), AND
!                LORENTZ HALF-WIDTH (GAMC) ALL DIVIDED BY PI [CM-2].
!       GAMD2    DOPPLER HALF-WIDTH SQUARED DIVIDED BY (2 LN2) [CM-2].
!       GAMC2    LORENTZ (COLLISION) HALF-WIDTH SQUARED [CM-2].
      REAL WDW,STORE,DENOM,DOPTRM,WDRATD,WDRATL,WDV,                    &
     &  WDTAIL,WDFRAC,ONEML2,RLINES,ELINES,GAMC2,GAMD2,SUGCPI,          &
     &  F1,F2,F3,STRPI,RATIO,RHO,RHOP1,RHO2,RHO2P1,RATIO2

!     FUNCTIONS:
!       BMTAIL   THE ABSORPTION FROM THE TAIL OF A VOIGT LINE DIVIDED
!                BY THE DISTANCE BETWEEN LINE CENTER AND LINE TAIL.
      REAL BMTAIL

!     THE RODGERS AND WILLIAMS APPROXIMATION IS USED TO
!     CALCULATE THE TOTAL VOIGT EQUIVALENT WIDTH (WDV):

!                2            2            2            2          2
!       (WDV/WDW)  = (WDL/WDW)  + (WDD/WDW)  - (WDL/WDW)  (WDD/WDW)

!     WHERE WDW IS THE WEAK-LINE EQUIVALENT WIDTH,
!           WDL IS THE LORENTZ EQUIVALENT WIDTH, AND
!           WDD IS THE DOPPLER EQUIVALENT WIDTH.

!     THE APPROXIMATIONS FOR THE EQUIVALENT WIDTHS FOLLOW:
!                                       _
!          WDW   =   <SU>         =  [ >_ SU/D ] / <1/D>

!                             2        A + X (B + 4 X)
!       WDRATL   =   (WDL/WDW)   =  ---------------------
!                                   A + X [C + X (D + X)]

!                             2
!       WDRATD   =   (WDD/WDW)   =  LN (1 + DOPTRM) / DOPTRM

!     WHERE

!                    X  =  WDW/<ALFCOL>
!                                                   2
!                    D  =  PI (95 PI - 208) / (32 PI  - 124 PI + 96)

!                    B  =  4 D - 2 PI

!                    C  =  8 PI D - 3 B - 32 PI

!                    A  =  2 PI (C - B)

!                                                   2
!               DOPTRM  =  [LN(2)/2]  (WDW/<ALFDOP>)

!     THE HALF-WIDTHS ARE CALCULATED FROM THE FOLLOWING EXPRESSIONS:
!                                                                 _
!       <ALFCOL>/WDW = [<ALFCOL/D> / <1/D>] / WDW = <ALFCOL/D> / >_ SU/D
!                                                                 _
!       <ALFDOP>/WDW = [<ALFDOP/D> / <1/D>] / WDW = <ALFDOP/D> / >_ SU/D

!     FOR OZONE AN IMPROVED APPROXIMATION FOR LORENTZ EQUIVALENT WIDTH
!     IS REQUIRED BASED ON THE WORK OF R. M. GOODY, WITH REFINEMENTS
!     BY L. S. BERNSTEIN TO ELIMINATE DECREASES IN THE CURVE-OF-GROWTH

!     DATA:
!        PT5LN2  =  .5 * LN(2)
!        P5OLN2     .5 / LN(2)
      REAL PT5LN2,P5OLN2,PI,PADEA,PADEB,PADEC,PADED
      DATA PT5LN2,P5OLN2,PI/.34657359,.72134752,3.141592654/,           &
     &  PADEA,PADEB,PADEC,PADED/258.45674,44.756426,85.891094,12.759903/

!     EFFECTIVE NUMBER OF LINES:
      ELINES=EDGENR*ODBAR
      RLINES=IBNDWD*ODBAR

!     WEAK-LINE EQUIVALENT WIDTH [CM-1]:
      WDW=DEPTH/ODBAR

!     DENOMINATOR IN THE LORENTZ EQUIVALENT WIDTH EXPRESSION:
      STORE=DEPTH/ACBAR
      RATIO=1.
      IF(ACBAR2.NE.0.)THEN

!         REFINEMENT FOR OZONE ONLY
          RHO=ACBAR2/ACBAR**2
          RHOP1=RHO+1.
          RHO2=RHO**2
          RHO2P1=RHO2+1.
          DENOM=RHOP1*(RHO2+RHOP1)
          F1=RHO*(RHO-1)*RHO2P1/(2*DENOM)
          DENOM=DENOM*RHOP1**2
          F2=RHO2*RHO2P1/DENOM

!         THE F3 EQUATION HAS BEEN IMPROVED (THE CURVE-OF-GROWTH
!         WAS DECREASING FOR SOME DOWN LOOKING PATHS ON THE WING
!         OF THE 9.6 MICRON O3 BAND; THIS NO LONGER OCCURS).
          F3=.7200*RHO-.5314
          STRPI=STORE/PI
          RATIO=1.-STRPI*F1/(2.+STRPI*(F3+STRPI*F2))
          STORE=RATIO*STORE
      ENDIF
      DENOM=PADEA+STORE*(PADEC+STORE*(PADED+STORE))
      RATIO2=RATIO**2

!     CALCULATE TERM IN EXPRESSION FOR DOPPLER EQUIVALENT WIDTH.
      DOPTRM=PT5LN2*(DEPTH/ADBAR)**2

!     VOIGT TOTAL EQUIVALENT WIDTH (RODGERS-WILLIAMS APPROXIMATION):
      IF(DOPTRM.GT..01)THEN

!         STANDARD EXPRESSIONS
          WDRATD=LOG(1.+DOPTRM)/DOPTRM
          WDRATL=(PADEA+STORE*(PADEB+4.*STORE))/DENOM*RATIO2
          WDV=WDW*SQRT(WDRATD+WDRATL-WDRATD*WDRATL)
      ELSE

!         WDRATD IS NEAR ONE.  CALCULATE ONE MINUS WDRATL.
          ONEML2=(PADEA*(1.-RATIO)*(1.+RATIO)+STORE*                    &
     &      (PADEC-RATIO2*PADEB+STORE*(PADED-4.*RATIO2+STORE)))/DENOM
          IF(DOPTRM.GT..0001)THEN

!             IF DOPTRM IS SMALL, WDRATD IS NEAR ONE.  REPLACE
!             THE LOG AND THE SQRT WITH POWER SERIES EXPANSIONS.
              STORE=DOPTRM*(.25-DOPTRM*(ONE6TH-.125*DOPTRM))*ONEML2
              WDV=WDW*(1.-STORE*(1.+.5*STORE*(1.+STORE)))
          ELSE

!             IF DOPTRM IS VERY SMALL, TRUNCATED EXPANSIONS SUFFICE.
              WDV=WDW*(1.-.25*DOPTRM*ONEML2)
          ENDIF
      ENDIF

!     SUBTRACT LINE TAIL CONTRIBUTIONS ASSUMING HALF-WIDTHS
!     SMALL COMPARED TO LINE CENTER TO LINE TAIL DISTANCE:
      GAMC2=ACBAR**2
      GAMD2=P5OLN2*ADBAR**2
      SUGCPI=DEPTH*ACBAR/PI
      WDTAIL=EDGENR*BMTAIL(SUGCPI,GAMD2,GAMC2,ELINES**2)                &
     &  +EDGEFR*BMTAIL(SUGCPI,GAMD2,GAMC2,(RLINES-ELINES)**2)
      WDFRAC=(WDV-WDTAIL)/IBNDWD

!     CALCULATE TRANSMITTANCE USING POWER LAW EXPRESSION.
      IF(WDFRAC.GE.1.)THEN
          BMTRAN=0.
      ELSE
          BMTRAN=(1.-WDFRAC)**RLINES
      ENDIF

!     RETURN TO BMOD:
      RETURN
      END
      REAL FUNCTION BMTRN(DEPTH,ODBAR,ADBAR,ACBAR)

!     THIS FUNCTION RETURNS PATH TRANSMITTANCE
!     AVERAGED OVER A SPECTRAL INTERVAL.

!     OUTPUTS:
!       BMTRN   =   PATH TRANSMITTANCE [DIMENSIONLESS]

!     INPUTS:
!                    _
!       DEPTH   =   >_  SU/D,     THE ABSORPTION COEFFICIENT (S/D)
!                                 AND COLUMN DENSITY (U) PRODUCT,
!                                 SUMMED OVER THE PATH [DIMENSIONLESS]

!       ODBAR   =   <1/D>,        THE RECIPROCAL OF THE LINE
!                                 SPACING, (D), PATH AVERAGED [CM]

!       ADBAR   =   <ALFDOP/D>,   THE DOPPLER LINE WIDTH (ALFDOP)
!                                 OVER THE AVERAGE LINE SPACING (D),
!                                 PATH AVERAGED [DIMENSIONLESS]

!       ACBAR   =   <ALFCOL/D>,   THE COLLISION (LORENTZ) LINE WIDTH
!                                 (ALFCOL) OVER THE AVERAGE LINE SPACING
!                                 (D), PATH AVERAGED [DIMENSIONLESS]

!     ARGUMENTS:
      REAL DEPTH,ODBAR,ADBAR,ACBAR

!     INCLUDE PARAMETERS:
!       ONE6TH   ONE-SIXTH.
      REAL ONE6TH
      PARAMETER(ONE6TH=1./6.)
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'BMHEAD.h'
!     LOCAL VARIABLES:
!       RLINES   EFFECTIVE NUMBER OF LINES IN THE INTERVAL.
!       ELINES   EFFECTIVE NUMBER OF LINES FROM LINE CENTER
!                TO THE NEAREST INTERVAL EDGE.
!       SUGCPI   PRODUCT OF LINE STRENGTH (S), COLUMN AMOUNT (U), AND
!                LORENTZ HALF-WIDTH (GAMC) ALL DIVIDED BY PI [CM-2].
!       GAMD2    DOPPLER HALF-WIDTH SQUARED DIVIDED BY (2 LN2) [CM-2].
!       GAMC2    LORENTZ (COLLISION) HALF-WIDTH SQUARED [CM-2].
      REAL WDW,STORE,DENOM,DOPTRM,WDRATD,WDRATL,WDV,                    &
     &  WDTAIL,WDFRAC,ONEML2,RLINES,ELINES,GAMC2,GAMD2,SUGCPI

!     FUNCTIONS:
!       BMTAIL   THE ABSORPTION FROM THE TAIL OF A VOIGT LINE DIVIDED
!                BY THE DISTANCE BETWEEN LINE CENTER AND LINE TAIL.
!       VGTWID   VOIGT EQUIVALENT WIDTH IN A FINITE SPECTRAL
!                INTERVAL OVER THE INTERVAL WIDTH.
      REAL BMTAIL,VGTWID

!     THE RODGERS AND WILLIAMS APPROXIMATION IS USED TO
!     CALCULATE THE TOTAL VOIGT EQUIVALENT WIDTH (WDV):

!                2            2            2            2          2
!       (WDV/WDW)  = (WDL/WDW)  + (WDD/WDW)  - (WDL/WDW)  (WDD/WDW)

!     WHERE WDW IS THE WEAK-LINE EQUIVALENT WIDTH,
!           WDL IS THE LORENTZ EQUIVALENT WIDTH, AND
!           WDD IS THE DOPPLER EQUIVALENT WIDTH.

!     THE APPROXIMATIONS FOR THE EQUIVALENT WIDTHS FOLLOW:
!                                       _
!          WDW   =   <SU>         =  [ >_ SU/D ] / <1/D>

!                             2        A + X (B + 4 X)
!       WDRATL   =   (WDL/WDW)   =  ---------------------
!                                   A + X [C + X (D + X)]

!                             2
!       WDRATD   =   (WDD/WDW)   =  LN (1 + DOPTRM) / DOPTRM

!     WHERE

!                    X  =  WDW/<ALFCOL>
!                                                   2
!                    D  =  PI (95 PI - 208) / (32 PI  - 124 PI + 96)

!                    B  =  4 D - 2 PI

!                    C  =  8 PI D - 3 B - 32 PI

!                    A  =  2 PI (C - B)

!                                                   2
!               DOPTRM  =  [LN(2)/2]  (WDW/<ALFDOP>)

!     THE HALF-WIDTHS ARE CALCULATED FROM THE FOLLOWING EXPRESSIONS:
!                                                                 _
!       <ALFCOL>/WDW = [<ALFCOL/D> / <1/D>] / WDW = <ALFCOL/D> / >_ SU/D
!                                                                 _
!       <ALFDOP>/WDW = [<ALFDOP/D> / <1/D>] / WDW = <ALFDOP/D> / >_ SU/D

!     DATA:
!        PT5LN2  =  .5 * LN(2)
!        P5OLN2     .5 / LN(2)
!        PADECB  =  PADEC - PADEB
!        PADED4  =  PADED - 4
      REAL PT5LN2,P5OLN2,PI,PADEA,PADEB,PADEC,PADED,PADECB,PADED4
      DATA PT5LN2,P5OLN2,PI/.34657359,.72134752,3.141592654/,           &
     &  PADEA,PADEB,PADEC,PADED,PADECB,PADED4/                          &
     &  258.45674,44.756426,85.891094,12.759903,41.134668,8.759903/

!     EFFECTIVE NUMBER OF LINES:
      ELINES=EDGENR*ODBAR
      RLINES=IBNDWD*ODBAR

!     CHECK DOPPLER WIDTH:
      IF(ADBAR.GT..2*ELINES)THEN

!         DOPPLER WIDTH IS MORE THAN ONE-FIFTH OF DISTANCE FROM LINE
!         CENTER TO BIN EDGE.  THE LINE TAIL EXPANSION IS NOT ACCURATE
!         UNDER THESE CONDITIONS.  COMPUTE THE VOIGT FINITE BIN
!         EQUIVALENT WIDTH BY EXPLICITLY PERFORMING THE NUMERICALLY
!         INTEGRATION.  SINCE THE LINE WIDTH IS LARGE, THE INTEGRATION
!         CAN BE PERFORMED ON A REASONABLE COARSE SPECTRAL GRID.
          WDFRAC=VGTWID(DEPTH,ACBAR,ADBAR,RLINES,-ELINES)
      ELSE

!         WEAK-LINE EQUIVALENT WIDTH [CM-1]:
          WDW=DEPTH/ODBAR

!         DENOMINATOR IN THE LORENTZ EQUIVALENT WIDTH EXPRESSION:
          STORE=DEPTH/ACBAR
          DENOM=PADEA+STORE*(PADEC+STORE*(PADED+STORE))

!         CALCULATE TERM IN EXPRESSION FOR DOPPLER EQUIVALENT WIDTH.
          DOPTRM=PT5LN2*(DEPTH/ADBAR)**2

!         VOIGT TOTAL EQUIVALENT WIDTH (RODGERS-WILLIAMS APPROXIMATION):
          IF(DOPTRM.GT..01)THEN

!             STANDARD EXPRESSIONS
              WDRATD=LOG(1.+DOPTRM)/DOPTRM
              WDRATL=(PADEA+STORE*(PADEB+4.*STORE))/DENOM
              WDV=WDW*SQRT(WDRATD+WDRATL-WDRATD*WDRATL)
          ELSE

!             WDRATD IS NEAR ONE.  CALCULATE ONE MINUS WDRATL.
              ONEML2=STORE*(PADECB+STORE*(PADED4+STORE))/DENOM
              IF(DOPTRM.GT..0001)THEN

!                 IF DOPTRM IS SMALL, WDRATD IS NEAR ONE.  REPLACE
!                 THE LOG AND THE SQRT WITH POWER SERIES EXPANSIONS.
                  STORE=DOPTRM*(.25-DOPTRM*(ONE6TH-.125*DOPTRM))*ONEML2
                  WDV=WDW*(1.-STORE*(1.+.5*STORE*(1.+STORE)))
              ELSE

!                 IF DOPTRM IS VERY SMALL, TRUNCATED EXPANSIONS SUFFICE.
                  WDV=WDW*(1.-.25*DOPTRM*ONEML2)
              ENDIF
          ENDIF

!         SUBTRACT LINE TAIL CONTRIBUTIONS ASSUMING HALF-WIDTHS
!         SMALL COMPARED TO LINE CENTER TO LINE TAIL DISTANCE:
          GAMC2=ACBAR**2
          GAMD2=P5OLN2*ADBAR**2
          SUGCPI=DEPTH*ACBAR/PI
          WDTAIL=EDGENR*BMTAIL(SUGCPI,GAMD2,GAMC2,ELINES**2)            &
     &      +EDGEFR*BMTAIL(SUGCPI,GAMD2,GAMC2,(RLINES-ELINES)**2)
          WDFRAC=(WDV-WDTAIL)/IBNDWD
      ENDIF

!     CALCULATE TRANSMITTANCE USING POWER LAW EXPRESSION.
      IF(WDFRAC.GE.1.)THEN
          BMTRN=0.
      ELSE
          BMTRN=(1.-WDFRAC)**RLINES
      ENDIF

!     RETURN TO BMOD:
      RETURN
      END
      REAL FUNCTION BMTAIL(SUGCPI,GAMD2,GAMC2,DFREQ2)

!     BMTAIL COMPUTES THE ABSORPTION FROM THE TAIL OF A VOIGT LINE
!     DIVIDED BY THE DISTANCE BETWEEN LINE CENTER AND LINE TAIL.  THE
!     CALCULATIONS ASSUME THAT THE LORENTZ AND DOPPLER HALF-WIDTHS ARE
!     SMALL WHEN COMPARED TO THE LINE TAIL TO LINE CENTER DISTANCE.

!                       INF
!                   1   /   /              2  \
!          BMTAIL = --  |   | 1 - EXP[-T(N) ] |  D N
!                   NU  /   \                 /
!                       NU

!     WHERE
!                                                          2
!                                1/2  INF        - (LN 2) Z
!                 2    2 / LN 2 \     /         E            DZ
!             T(N)  = X  | ---- |     |     -----------------------   ,
!                        \  PI  /     /         2                 2
!                                    -INF   GAMC  + ( N - GAMD Z )

!              2   S U GAMC
!             X  = --------   ,
!                     PI

!             S = LINE STRENGTH   ,   AND

!             U = COLUMN DENSITY

!     FOR LORENTZ (GAMC) AND DOPPLER (GAMD) HALF-WIDTHS SMALL
!     COMPARED TO THE DISTANCE TO LINE TAIL (NU),

!                    2  /          2     2                       \
!               2   X   |      3 GD    GD  /      2         2 \  |
!          T(NU)  = --  |  1 + ----- + --- | 15 GD  - 4 GAMC  |  |   ,
!                    2  |         2      4 \                  /  |
!                   D   \        D      D                        /

!              2     2       2
!             D  = NU  + GAMC    ,

!               2       2
!             GD  = GAMD  /  (2 LN 2)   ,

!                     /                2   \
!          BMTAIL = - | 1 - EXP [ -T(N)  ] |
!                     \                    /

!                    T(NU)            /           2       2     2
!                2 X   /          2   |     / 3 GD  - GAMC  \  T
!              + ---   |   EXP (-T )  | 1 + | ------------- |  --
!                 NU   /              |     \       2       /   2
!                      0              \                        X

!                                4        2     2       4     4  \
!                         / 15 GD  - 10 GD  GAMC  - GAMC  \  T   |
!                       + | ----------------------------- |  --  |  D T
!                         \                8              /   4  |
!                                                            X   /

!     THE ERROR FUNCTION INTEGRALS ARE EVALUATED DIFFERENTLY
!     DEPENDING ON THE MAGNITUDE OF T(NU)**2 [VARIABLE "OPTDEP"].

!     INPUT ARGUMENTS (NOTE THAT THE ARGUMENTS CAN ALL BE SCALED BY
!     AN ARBITRARY POSITIVE CONSTANT WITHOUT CHANGING THE FUNCTION).
!       SUGCPI   PRODUCT OF LINE STRENGTH (S), COLUMN AMOUNT (U),
!                AND LORENTZ HALF-WIDTH (GAMC) ALL DIVIDED BY PI
!                TIMES THE LINE SPACING.
!       GAMD2    DOPPLER HALF-WIDTH SQUARED DIVIDED BY THE PRODUCT
!                OF (2 LN2) AND THE LINE SPACING SQUARED.
!       GAMC2    LORENTZ (COLLISION) HALF-WIDTH SQUARED DIVIDED
!                BY THE LINE SPACING SQUARED.
!       DFREQ2   DISTANCE BETWEEN LINE CENTER AND LINE TAIL
!                SQUARED DIVIDED BY THE LINE SPACING SQUARED.
      REAL SUGCPI,GAMD2,GAMC2,DFREQ2

!     PARAMETERS:
!       PCON     CONSTANT USED IN RATIONAL APPROXIMATION TO THE ERROR
!                FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A1       LINEAR COEFFICIENT USED IN RATIONAL APPROXIMATION TO
!                THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A2       QUADRATIC COEFFICIENT USED IN RATIONAL APPROXIMATION
!                TO THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A3       CUBIC COEFFICIENT USED IN RATIONAL APPROXIMATION TO
!                THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A4       QUARTIC COEFFICIENT USED IN RATIONAL APPROXIMATION TO
!                THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       A5       5TH-ORDER COEFFICIENT USED IN RATIONAL APPROXIMATION
!                TO THE ERROR FUNCTION [ABRAMOWITZ AND STEGUN, 7.1.26].
!       RTPI     SQUARE ROOT OF PI.
!       RTPI4    SQUARE ROOT OF PI DIVIDED BY 4.
!       RTPIC    SQUARE ROOT OF PI TIMES 3/32.
!       REC3     RECIPROCAL OF 3.
!       REC6     RECIPROCAL OF 6.
!       REC10    RECIPROCAL OF 10.
!       REC28    RECIPROCAL OF 28.
!       REC40    RECIPROCAL OF 40.
!       REC56    RECIPROCAL OF 56.
!       REC144   RECIPROCAL OF 144.
      REAL PCON,A1,A2,A3,A4,A5,RTPI,RTPI4,RTPIC,                        &
     &  REC3,REC6,REC10,REC28,REC40,REC56,REC144
      PARAMETER(PCON=.3275911,A1=.451673692,A2=-.504257335,             &
     &  A3=2.519390259,A4=-2.575644906,A5=1.881292140,                  &
     &  RTPI=A1+A2+A3+A4+A5,RTPI4=RTPI/4,RTPIC=3*RTPI/32,               &
     &  REC3=1./3,REC6=1./6,REC10=1./10,REC28=1./28,                    &
     &  REC40=1./40,REC56=1./56,REC144=1./144)

!     LOCAL VARIABLES:
!       COEF2    COEFFICIENT OF QUADRATIC TERM IN VOIGT LINE TAIL
!                SERIES EXPANSION [CM-2].
!       COEF4    COEFFICIENT OF QUARTIC TERM IN VOIGT LINE TAIL
!                SERIES EXPANSION [CM-4].
!       C4RAT    COEF4 DIVIDED BY SUGCPI [CM-2].
!       DENOM    LORENTZ HALF-WIDTH SQUARED PLUS THE SQUARE OF THE
!                DISTANCE BETWEEN THE LINE TAIL AND LINE CENTER [CM-2].
!       VOIGT    VOIGT LINE SHAPE EVALUATE AT DFREQ DIVIDED BY THE
!                LORENTZ HALF-WIDTH SQUARED [CM2].
!       OPTDEP   OPTICAL DEPTH AT FREQUENCY DISPLACEMENT DFREQ.
!       RTDEP    SQUARE ROOT OF THE OPTICAL DEPTH.
!       PHI      PRODUCT OF RTPI, THE EXPONENTIAL FUNCTION EVALUATED AT
!                OPTDEP, AND THE COMPLEMENTARY ERROR FUNCTION OF RTDEP.
!       RTRAT    SQUARE ROOT OF SUGCPI OVER DFREQ2.
!       RAPP     RATIO USED IN RATIONAL APPROXIMATION TO ERROR FUNCTION.
      REAL DENOM,VOIGT,OPTDEP,RTDEP,COEF2,COEF4,C4RAT,PHI,RTRAT,RAPP

!     DEFINE REQUIRED LOCAL VARIABLES:
      COEF2=3*GAMD2-GAMC2
      COEF4=5*(COEF2-GAMC2)*GAMD2-GAMC2**2
      DENOM=GAMC2+DFREQ2
      VOIGT=(1.+GAMD2*(3.+(5*COEF2+GAMC2)/DENOM)/DENOM)/DENOM
      OPTDEP=VOIGT*SUGCPI

!     BRANCH BASED ON MAGNITUDE OF THE OPTICAL DEPTH:
      IF(OPTDEP.LE..01)THEN

!         SMALL OPTICAL DEPTH:
          BMTAIL=2*SUGCPI*SQRT(VOIGT/DFREQ2)                            &
     &      *(1.-(REC3-REC10*OPTDEP)*OPTDEP                             &
     &       +VOIGT*(COEF2*(REC6-(REC10-OPTDEP*REC28)*OPTDEP)           &
     &              +VOIGT*COEF4*(REC40-(REC56-REC144*OPTDEP)*OPTDEP))) &
     &      -OPTDEP*(1.-(.5-REC6*OPTDEP)*OPTDEP)
      ELSE

!         LARGE OPTICAL DEPTH CONTRIBUTION:
          C4RAT=COEF4/SUGCPI
          RTRAT=SQRT(SUGCPI/DFREQ2)
          BMTAIL=RTRAT*(RTPI+(RTPI4*COEF2+RTPIC*C4RAT)/SUGCPI)-1.
          IF(OPTDEP.LT.10.24)THEN

!             INTERMEDIATE OPTICAL DEPTH:
              RTDEP=SQRT(OPTDEP)
              IF(RTDEP.LT.2.33)THEN

!                 RATIONAL APPROXIMATION TO ERROR FUNCTION:
                  RAPP=1./(1.+PCON*RTDEP)
                  PHI=RAPP*(A1+RAPP*(A2+RAPP*(A3+RAPP*(A4+RAPP*A5))))
              ELSE

!                 CONTINUED FRACTION REPRESENTATION FOR ERROR FUNCTION:
                  PHI=(2.+OPTDEP*(4.5+OPTDEP))                          &
     &              /(RTDEP*(3.75+OPTDEP*(5+OPTDEP)))
              ENDIF
              BMTAIL=BMTAIL+EXP(-OPTDEP)                                &
     &          *(1.-RTRAT*(PHI+(2*COEF2*(2*RTDEP+PHI)                  &
     &          +C4RAT*(RTDEP*(OPTDEP+1.5)+.75*PHI))/(8*SUGCPI)))
          ENDIF
      ENDIF
      RETURN
      END
      REAL FUNCTION VGTWID(DEPTH,ACBAR,ADBAR,RLINES,FMIN)

!     VGTWID RETURNS THE VOIGT EQUIVALENT WIDTH IN A FINITE
!     SPECTRAL INTERVAL OVER THE INTERVAL WIDTH.

!     ARGUMENTS:
!       DEPTH    PRODUCT OF LINE STRENGTH AND COLUMN AMOUNT OVER
!                LINE SPACING.
!       ACBAR    LORENTZ HALF-WIDTH OVER THE LINE SPACING.
!       ADBAR    HALF-WIDTH OVER THE LINE SPACING.
!       RLINES   EFFECTIVE NUMBER OF LINES IN SPECTRAL BIN (EQUAL
!                TO THE SPECTRAL BIN WIDTH OVER THE LINE SPACING).
!       FMIN     SPECTRAL BIN MINIMUM FREQUENCY OVER THE LINE SPACING.
        REAL DEPTH,ACBAR,ADBAR,RLINES,FMIN

!     LOCAL VARIABLES:
!       NPOINT   ONE LESS THAN THE NUMBER OF INTEGRATION POINTS.
!       IPOINT   INDEX FOR INTEGRATION POINTS.
!       COEF     COEFFICIENTS USED IN EQUIVALENT WIDTH CALCULATIONS:
!       X        REAL PART OF THE COMPLEX PROBABILITY FUNCTION
!                ARGUMENT [EQUAL TO SQRT(LN 2) TIMES THE FREQUENCY
!                DISPLACEMENT OVER THE DOPPLER HALF-WIDTH].
!       DELX     THE INTERGRATION STEP SIZE.
!       Y1P5     IMAGINARY PART OF THE COMPLEX PROBABILITY FUNCTION
!                ARGUMENT [EQUAL TO SQRT(LN 2) TIMES THE LORENTZ
!                HALF-WIDTH OVER THE DOPPLER HALF-WIDTH] PLUS 1.5.
!       Y1P5SQ   Y1P5 SQUARED.
      INTEGER NPOINT,IPOINT
      REAL COEF,X,DELX,Y1P5,Y1P5SQ

!     FUNCTIONS:
!       RCPF12   REAL PART OF THE COMPLEX PROBABILITY FUNCTION.
      REAL RCPF12

!     DATA:
!       RTPI     SQUARE ROOT OF PI.
!       RTLN2    SQUARE ROOT OF LN(2).
      REAL RTPI,RTLN2
      DATA RTPI,RTLN2/1.7724539,.83255461/

!     INITIALIZATIONS:
      NPOINT=4*(INT(RLINES/MAX(ACBAR,ADBAR))+1)
      COEF=RTLN2/ADBAR
      X=COEF*FMIN
      DELX=COEF*RLINES/NPOINT
      Y1P5=COEF*ACBAR+1.5
      Y1P5SQ=Y1P5**2
      COEF=-DEPTH*COEF/RTPI
      VGTWID=EXP(COEF*RCPF12(X,Y1P5,Y1P5SQ))/2

!     SET Y1P5 TO ZERO TO INDICATE THAT ITS VALUE
!     WAS PREVIOUSLY PASSED TO FUNCTION RCPF12.
      Y1P5=0.

!     LOOP OVER INTERMEDIATE INTEGRATION POINTS:
      DO 10 IPOINT=2,NPOINT
          X=X+DELX
          VGTWID=VGTWID+EXP(COEF*RCPF12(X,Y1P5,Y1P5SQ))
   10 CONTINUE

!     ADD FINAL INTEGRATION POINT
      X=X+DELX
      VGTWID=VGTWID+EXP(COEF*RCPF12(X,Y1P5,Y1P5SQ))/2
      VGTWID=1.-VGTWID/NPOINT
      RETURN
      END
      REAL FUNCTION RCPF12(X,Y1P5,Y1P5SQ)

!     RCPF12 RETURNS THE REAL PART OF THE COMPLEX PROBABILITY
!     FUNCTION [ W(Z)=EXP(-Z**2)*ERFC(-I*Z);   Z = X + IY ]
!     IN THE UPPER HALF-PLANE (I.E., FOR Y >= 0.).  THE MAXIMUM
!     RELATIVE ERROR OF RCPF12 IS <2.E-6.  THIS ROUTINE WAS
!     DEVELOPED BY J.HUMLICEK, JQSRT, VOL 21, P309 (1980)

!     ARGUMENTS:
!       X        REAL PART OF ARGUMENT.
!       Y1P5     IMAGINARY PART OF ARGUMENT PLUS 1.5.
!       Y1P5SQ   Y1P5 SQUARED.
      REAL X,Y1P5,Y1P5SQ

!     LOCAL VARIABLES:
!       XMH      X MINUS THE ROOT OF A HERMITE POLYNOMIAL.
!       XPH      X PLUS THE ROOT OF A HERMITE POLYNOMIAL.
!       CY1P5    COTANGENT DATA TIMES Y1P5.
      REAL XMH,XPH,CY1P5(6)
      SAVE CY1P5

!     DATA:
!       HERM12   ZEROES OF THE 12TH ORDER HERMITE POLYNOMIAL
!                (ABRAMOWITZ AND STEGUN, TABLE 25.10).

!                  WHRM12   9/4
!       SINDAT  =  ------  E     SIN( 3 HERM12 )  ,  WHERE WHRM12 IS
!                    PI

!                THE WEIGHT OF THE 12TH ORDER HERMITE POLYNOMIAL.
!       CTNDAT   COTANGENT DATA (DERIVED FROM COSDAT AND SINDAT).
      REAL HERM12(6),CTNDAT(6),SINDAT(6)
      SAVE HERM12,CTNDAT,SINDAT
      DATA HERM12/ 0.314240376,     0.947788391,     1.59768264,        &
     &             2.27950708,      3.02063703,      3.8897249     /,   &
     &     CTNDAT/ 0.72617082,     -3.25314144,     -0.08083430,        &
     &             1.61167866,     -2.63380063,     -0.798131224   /,   &
     &     SINDAT/ 1.393237,        0.231152406,    -0.155351466,       &
     &             6.21836624E-03,  9.19082986E-05, -6.27525958E-07/
!OLD
!OLD  DATA COSDAT/ 1.01172805,     -0.75197147,      1.2557727E-02,
!OLD 1             1.00220082E-02, -2.42068135E-04,  5.00848061E-07/,

!     Y1P5 IS INPUT AS ZERO TO INDICATE THAT ARRAY CY1P5
!     HAS ALREADY BEEN DEFINED:
      IF(Y1P5.GT.0.)THEN

!         LOOP OVER DATA:
          CY1P5(1)=CTNDAT(1)*Y1P5
          XMH=X-HERM12(1)
          XPH=X+HERM12(1)
          RCPF12=SINDAT(1)*((CY1P5(1)-XMH)/(XMH**2+Y1P5SQ)              &
     &                     +(CY1P5(1)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(2)=CTNDAT(2)*Y1P5
          XMH=X-HERM12(2)
          XPH=X+HERM12(2)
          RCPF12=RCPF12+SINDAT(2)*((CY1P5(2)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(2)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(3)=CTNDAT(3)*Y1P5
          XMH=X-HERM12(3)
          XPH=X+HERM12(3)
          RCPF12=RCPF12+SINDAT(3)*((CY1P5(3)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(3)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(4)=CTNDAT(4)*Y1P5
          XMH=X-HERM12(4)
          XPH=X+HERM12(4)
          RCPF12=RCPF12+SINDAT(4)*((CY1P5(4)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(4)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(5)=CTNDAT(5)*Y1P5
          XMH=X-HERM12(5)
          XPH=X+HERM12(5)
          RCPF12=RCPF12+SINDAT(5)*((CY1P5(5)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(5)+XPH)/(XPH**2+Y1P5SQ))
          CY1P5(6)=CTNDAT(6)*Y1P5
          XMH=X-HERM12(6)
          XPH=X+HERM12(6)
          RCPF12=RCPF12+SINDAT(6)*((CY1P5(6)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(6)+XPH)/(XPH**2+Y1P5SQ))
      ELSE

!         LOOP OVER DATA:
          XMH=X-HERM12(1)
          XPH=X+HERM12(1)
          RCPF12=SINDAT(1)*((CY1P5(1)-XMH)/(XMH**2+Y1P5SQ)              &
     &                     +(CY1P5(1)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(2)
          XPH=X+HERM12(2)
          RCPF12=RCPF12+SINDAT(2)*((CY1P5(2)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(2)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(3)
          XPH=X+HERM12(3)
          RCPF12=RCPF12+SINDAT(3)*((CY1P5(3)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(3)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(4)
          XPH=X+HERM12(4)
          RCPF12=RCPF12+SINDAT(4)*((CY1P5(4)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(4)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(5)
          XPH=X+HERM12(5)
          RCPF12=RCPF12+SINDAT(5)*((CY1P5(5)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(5)+XPH)/(XPH**2+Y1P5SQ))
          XMH=X-HERM12(6)
          XPH=X+HERM12(6)
          RCPF12=RCPF12+SINDAT(6)*((CY1P5(6)-XMH)/(XMH**2+Y1P5SQ)       &
     &                            +(CY1P5(6)+XPH)/(XPH**2+Y1P5SQ))
      ENDIF
      RETURN
      END
