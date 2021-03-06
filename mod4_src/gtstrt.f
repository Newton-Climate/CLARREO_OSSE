      SUBROUTINE GTSTRT(IV,IREC,ITB,LSTREC)

!     FIND THE RECORD (IREC) WHERE FREQUENCY IV STARTS IN THE
!     BAND MODEL PARAMETER TAPE.  ITB IS THE UNIT NUMBER.

!     DECLARE INPUTS/OUTPUTS
      INTEGER IV,IREC,ITB,LSTREC

!     DECLARE LOCAL VARIABLES
      INTEGER IV2,IRCLO,IRCHI,IRCTST,IDR
      IRCLO=2

!     SET IRCHI TO THE MAXIMUM RECORD.
      IRCHI=LSTREC

!     START BISECTION LOOP.  FIND FREQUENCY AT THE MID-RECORD
!     BETWEEN IRCLO AND IRCHI.  RESET IRCLO AND IRCHI.
   10 CONTINUE
      IDR=IRCHI-IRCLO
      IF(IDR.GT.1)THEN
          IRCTST=IRCLO+IDR/2
          READ(ITB,REC=IRCTST,ERR=20)IV2
          IF(IV2.LT.IV)THEN
              IRCLO=IRCTST
              GOTO10
          ENDIF
   20     CONTINUE
          IRCHI=IRCTST
          GOTO10
      ENDIF
      READ(ITB,REC=IRCLO)IV2
      IF(IV2.EQ.IV)THEN
         IREC=IRCLO
      ELSE
         IREC=IRCHI
      ENDIF
      END
