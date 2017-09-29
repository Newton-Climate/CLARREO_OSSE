      SUBROUTINE BMLOAD

!     BMLOAD IS CALLED BY BMOD TO LOAD BAND MODEL DATA FOR A
!     SINGLE PARAMETER SET INTO THE MATRICES SD, OD AND ALF0.

!     ORIGINALLY, THE BAND MODEL PARAMETERS ARE STORE IN
!     COMPACT PARAMETER SETS.  THERE ARE TWO TYPES OF PARAMETER
!     SETS CORRESPONDING TO LINE CENTER AND LINE TAIL
!     CONTRIBUTIONS.  IF IMOL IS 12 OR LESS, THE PARAMETER SET
!     CONTAINS LINE CENTER MOLECULAR ABSORPTION INFORMATION:

!      IBIN = THE BIN NUMBER (IBIN*IBNDWD = BIN CENTER FREQUENCY)
!      IMOL =  1 FOR H2O
!           =  2 FOR CO2
!           =  3 FOR O3
!           =  4 FOR N2O
!           =  5 FOR CO
!           =  6 FOR CH4
!           =  7 FOR O2
!           =  8 FOR NO
!           =  9 FOR SO2
!           = 10 FOR NO2
!           = 11 FOR NH3
!           = 12 FOR HNO3

!        SDZ(IT)= AVERAGE MOLECULAR ABSORPTION COEFFICIENT PARAMETER
!                 FOR THE IT'TH TEMPERATURE (CM-1/AMAGAT)
!          IALF = LORENTZ LINE WIDTH PARAMETER AT STANDARD PRESSURE AND
!                 TEMPERATURE (1.E-04 CM-1/ATM)
!        ODZ(IT)= LINE DENSITY PARAMETER FOR THE TEMPERATURE IT (CM-1)

!     IF IMOL IS 12 OR MORE, THE PARAMETER SET CONTAINS CONTINUUM
!     ABSORPTION INFORMATION AND

!                  IBIN = THE BIN NUMBER
!             IMOL,IALF = 13 FOR  H2O TAIL CONTRIBUTION
!                       = 14 FOR  CO2 TAIL CONTRIBUTION
!                       = 15 FOR   O3 TAIL CONTRIBUTION
!                       = 16 FOR  N2O TAIL CONTRIBUTION
!                       = 17 FOR   CO TAIL CONTRIBUTION
!                       = 18 FOR  CH4 TAIL CONTRIBUTION
!                       = 19 FOR   O2 TAIL CONTRIBUTION
!                       = 20 FOR   NO TAIL CONTRIBUTION
!                       = 21 FOR  SO2 TAIL CONTRIBUTION
!                       = 22 FOR  NO2 TAIL CONTRIBUTION
!                       = 23 FOR  NH3 TAIL CONTRIBUTION
!                       = 24 FOR HNO3 TAIL CONTRIBUTION
!               SDZ(IT) = MOLECULAR LINE WING CONTINUUM ABSORPTION
!                         COEFFICIENT FOR SPECIES IMOL AT THE
!                         IT'TH TEMPERATURE (CM-1/AMAGAT)
!               ODZ(IT) = MOLECULAR LINE WING CONTINUUM ABSORPTION
!                         COEFFICIENT FOR SPECIES IALF AT THE
!                         IT'TH TEMPERATURE (CM-1/AMAGAT)
      INCLUDE 'PARAMS.h'

!     LIST COMMONS:
      INCLUDE 'BMHEAD.h'
      INCLUDE 'BMDAT.h'

!     LOCAL VARIABLES:
      INTEGER ITLSUB,IDIVSV,IM,IT,JM

!     LIST DATA:
      INTEGER IDIV
      SAVE IDIV
      DATA IDIV/1/

!     FILL THE SD, OD AND ALF0 MATRICES.
      IF(IMOL.LE.NMOL)THEN

!         LINE CENTER CONTRIBUTIONS
          IM=IMOL
          DO 20 IT=1,NTEMP
              SD(IT,IM,0)=SDZ(IT)
              OD(IT,IM)=ODZ(IT)
              DO 10 ITLSUB=1,NTLSUB
                  SD(IT,IM,ITLSUB)=0.
   10         CONTINUE
   20     CONTINUE
          ALF0(IM)=1.E-04*IALF
      ELSE

!         LINE TAIL CONTRIBUTIONS
          IM=IMOL-NMOL
          IF(IALF.EQ.0)THEN

!             ONE TAIL
              DO 30 IT=1,NTEMP
                  SD(IT,IM,IDIV)=SDZ(IT)
   30         CONTINUE
              IDIV=1
          ELSE

!             TWO TAILS
              JM=IALF-NMOL
              IDIVSV=IDIV
              IF(IDIV.EQ.NTLSUB)THEN
                  IDIV=1
              ELSE
                  IDIV=IDIV+1
              ENDIF
              DO 40 IT=1,NTEMP
                  SD(IT,IM,IDIVSV)=SDZ(IT)
                  SD(IT,JM,IDIV  )=ODZ(IT)
   40         CONTINUE
              IF(IDIV.EQ.NTLSUB)THEN
                  IDIV=1
              ELSE
                  IDIV=IDIV+1
              ENDIF
          ENDIF
      ENDIF
      RETURN
      END
