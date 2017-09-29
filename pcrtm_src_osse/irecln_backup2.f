      INTEGER FUNCTION IRECLN(NVAR,LABEL)

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     FIND RECORD LENGTH FOR "NVAR" INTEGER*4 VARIABLES.  FIRST FIND
!     RECORD LENGTH FOR 2 INTEGERS AND THEN FOR 3 TO FIND ANY OFFSET.

!     INPUT ARGUMENTS:
!       NVAR     NUMBER OF REAL*4 OR INTEGER*4 VARIABLES
      INTEGER NVAR
      CHARACTER LABEL*6

!     BEGIN DRF MODIFICATION
!     COMMON BLOCKS
!     /SUNFLG/
!       LNFLRT   LENGTH OF I/O FILE ROOT NAME, 0 IF NO mod5root.in FILE.
!       LRDSUN   SOLAR DATA FLAG, TRUE IF IRRADIANCES IN COMMON BLOCK
!                /SOL01/ HAVE BEEN MODIFIED FROM THE BLOCK DATA.
      INTEGER LNFLRT
      COMMON/SUNFLG/LNFLRT

!     /SUNNAM/
!       SUNFIL   NAME OF DEFAULT FILE CONTAINING SOLAR IRRADIANCE DATA.
!       FLRT     ROOT NAME FOR ALL I/O FILES.
      CHARACTER FLRT*(NAMLEN-4)
      COMMON/SUNNAM/FLRT
!     END DRF MODIFICATION

!     LOCAL VARIABLES:
!       LUNIT    UNIT NUMBER FOR SCRATCH FILE.
!       IRECL2   TWO INTEGER RECORD LENGTH LOOP INDEX.
!       IRECL3   THREE INTEGER RECORD LENGTH LOOP INDEX.
!       ITEST0   TEST COMPARE VARIABLE FOR IREC.
!       ITEST1   TEST COMPARE VARIABLE FOR IREC+1.
!       ITEST2   TEST COMPARE VARIABLE FOR IREC+2.
!       IOS      RESULT OF IOSTAT TEST.
      INTEGER LUNIT,IRECL2,IRECL3,ITEST0,ITEST1,ITEST2,IOS

!     FUNCTIONS:
!       NUNIT    RETURNS AN UNUSED FILE UNIT NUMBER.
      EXTERNAL NUNIT
      INTEGER NUNIT

!     FETCH FILE UNIT NUMBER AND REMOVE OLD FILE:
      LUNIT=NUNIT()

!     LOOP OVER POSSIBLE RECORD LENGTHS FOR TWO INTEGERS:
      DO IRECL2=2,24

!        WRITE 2 INTEGERS TO FIRST 2 RECORDS OF SCRATCH FILE:
!         CALL OPNFL(LUNIT,IRECL2,'IRECLN_FILE'//LABEL,'UNKNOWN',        &
!     &     'UNFORMATTED','IRECLN')

         !BEGIN DRF MODIFICATION
        CALL OPNFL(LUNIT,IRECL2,'/tmp/'//FLRT(3:LNFLRT)//'.IRECLN_FILE',&
     &              'UNKNOWN','UNFORMATTED','IRECLN')
         !END DRF MODIFICATION

         WRITE(LUNIT,REC=1,ERR=20)IRECL2,IRECL2+1
         WRITE(LUNIT,REC=2,ERR=20)IRECL2,IRECL2+1
         CLOSE(LUNIT)

!        CHECK IF 2 INTEGERS CAN BE READ WITH CURRENT RECORD LENGTH.
         CALL OPNFL(LUNIT,IRECL2,'IRECLN_FILE'//LABEL,'OLD',            & 
     &              'UNFORMATTED','IRECLN')

         !BEGIN DRF MODIFICATION
        CALL OPNFL(LUNIT,IRECL2,'/tmp/'//FLRT(3:LNFLRT)//'.IRECLN_FILE',&
     &              'OLD','UNFORMATTED','IRECLN')
         !END DRF MODIFICATION

         READ(LUNIT,REC=2,ERR=20)ITEST0,ITEST1
         IF(ITEST0.EQ.IRECL2 .AND. ITEST1.EQ.IRECL2+1)THEN

!           SUCCESSFULLY READ:  DETERMINE RECORD LENGTH FOR 3 INTEGERS.
            DO IRECL3=IRECL2+1,IRECL2+4

!              WRITE 3 INTEGERS TO FIRST 2 RECORDS OF SCRATCH FILE:
               CLOSE(LUNIT,STATUS='DELETE',IOSTAT=IOS)
               CALL OPNFL(LUNIT,IRECL3,'IRECLN_FILE'//LABEL,'UNKNOWN',  &
     &           'UNFORMATTED','IRECLN')

!BEGIN DRF MODIFICATION
               CALL OPNFL(LUNIT,IRECL3,                                 &
     &                  '/tmp/'//FLRT(3:LNFLRT)//'.IRECLN_FILE',        &
     &                  'UNKNOWN','UNFORMATTED','IRECLN')
!END DRF MODIFICATION
               WRITE(LUNIT,REC=1,ERR=10)IRECL3,IRECL3+1,IRECL3+2
               WRITE(LUNIT,REC=2,ERR=10)IRECL3,IRECL3+1,IRECL3+2
               CLOSE(LUNIT)

!              CHECK IF 3 INTEGERS CAN BE READ.

!BEGIN DRF MODIFICATION
!               CALL OPNFL(LUNIT,IRECL3,                                 &
!     &                    '/tmp/'//FLRT(3:LNFLRT)//'.IRECLN_FILE',      &
!     &                     'OLD','UNFORMATTED','IRECLN')

!DRF               CALL OPNFL(LUNIT,IRECL3,'IRECLN_FILE','OLD',             &
!DRF     &           'UNFORMATTED','IRECLN')
!END DRF MODIFICATION


               READ(LUNIT,REC=2,ERR=10)ITEST0,ITEST1,ITEST2
               IF(ITEST0.EQ.IRECL3 .AND. ITEST1.EQ.IRECL3+1             &
     &                             .AND. ITEST2.EQ.IRECL3+2)THEN

!                 RETURN RECORD LENGTH AS OFFSET (3*IRECL2-2*IRECL3)
!                 PLUS NVAR TIMES VARIABLE LENGTH (IRECL3-IRECL2):
                  CLOSE(LUNIT,STATUS='DELETE')
                  IRECLN=NVAR*(IRECL3-IRECL2)+3*IRECL2-2*IRECL3
                  RETURN
               ENDIF
   10          CONTINUE
            ENDDO
            STOP 'IRECLN could not determine record length.'
         ENDIF
   20    CONTINUE
         CLOSE(LUNIT,STATUS='DELETE',IOSTAT=IOS)
      ENDDO
      STOP 'IRECLN could not determine record length.'
      END
