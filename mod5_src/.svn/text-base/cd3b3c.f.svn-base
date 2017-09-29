      SUBROUTINE CD3B3C(TP5_FLAG,TP5_FILE,TP5_LINENUMBER)

!     PROCESS CARD3B OR CARD3C INPUTS:
      IMPLICIT NONE

!     PARAMETERS:
      INCLUDE 'PARAMS.h'

!     ARGUMENTS:
!DRF    TP5_FLAG LOGICAL FLAG IF READING IN TP5 FILE
!DRF    TP5_FILE STRING CONTAINING TP5 FILE INFORMATION
!DRF    TP5_LINENUMBER (I/O) CURRENT LINENUMBER FOR TP5 READING
      INTEGER TP5_LINENUMBER
      LOGICAL TP5_FLAG
      CHARACTER TP5_FILE(1000)*110
      INTEGER NUM_REPEATS,DUMMY,NUM_REMAINDER

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'SOLS.h'

!     /USRDTA/
      INTEGER NAERF,NWLF
      REAL WLF,F
      COMMON/USRDTA/NAERF,NWLF,WLF(MWLF),F(MAERF,MANGLS,MWLF)

!     LOCAL VARIABLES:
!       IANG     LOOP INDEX FOR SCATTERING ANGLE GRID.
!       IWAV     LOOP INDEX FOR SPECTRAL WAVELENGTH GRID.
!       IAER     LOOP INDEX FOR AEROSOL TYPE.
      INTEGER IANG,IWAV,IAER

!     USER-DEFINED PHASE FUNCTIONS ARE SPECIFIED FOR ALL 4 AEROSOLS:
      NAERF=4
      IF(LJMASS)THEN
          CALL INITCD('CARD3B1')
          IF(NWLF.LE.0)THEN
              CALL INITCD('CARD3B2')
              NWLF=1
              WLF(1)=.55
          ELSE

!             "CARD3C" INCLUDES ALL THE CARD3CN DATA:
               CALL INITCD('CARD3C')
          ENDIF
          RETURN
      ENDIF

!     CARD 3B1 USER DEFINED PHASE FUNCTION:
      IF (TP5_FLAG) THEN
         READ(IRD,'(3(I5))')NANGLS,NWLF
      ELSE
         READ(TP5_FILE(TP5_LINENUMBER),'(3(I5))')NANGLS,NWLF
         TP5_LINENUMBER=TP5_LINENUMBER+1
      ENDIF
      WRITE(IPR,'(A,3(1X,I4))')' CARD 3B1*****',NANGLS,NWLF
      IF(NANGLS.GT.MANGLS)THEN
          WRITE(IPR,'(A,I6,A,/18X,A,I6,A)')' Error in CD3B3C: '//       &
     &      ' Input number of phase function angles (',NANGLS,          &
     &      ') exceeds the maximum',' allowed (',MANGLS,                &
     &      ').  Increase parameter MANGLS in routine "PARAMS.h".'
          STOP 'CD3B3C error:  Too many phase function angles.'
      ENDIF
      IF(NWLF.GT.MWLF)THEN
          WRITE(IPR,'(A,I6,A,/18X,A,I6,A)')' Error in CD3B3C: '//       &
     &      ' Input number of phase function spectral points (',NWLF,   &
     &      ') exceeds the maximum',' allowed (',MWLF,                  &
     &      ').  Increase parameter MWLF in routine "PARAMS.h".'
          STOP 'CD3B3C error:  Too many phase function spectral points.'
      ENDIF
      IF(NWLF.LE.0)THEN

!         CARD 3B2:
          IF (TP5_FLAG) THEN
          READ(IRD,'((5F10.0))')(ANGF(IANG),F(1,IANG,1),F(2,IANG,1),    &
     &      F(3,IANG,1),F(4,IANG,1),IANG=1,NANGLS)
          ELSE
            READ(TP5_FILE(TP5_LINENUMBER),                              &
     &         '((5F10.0))')(ANGF(IANG),F(1,IANG,1),F(2,IANG,1),        &
     &      F(3,IANG,1),F(4,IANG,1),IANG=1,NANGLS)
            TP5_LINENUMBER=TP5_LINENUMBER+1
          ENDIF
          WRITE(IPR,'(/A,0P,F10.5,1P,4E10.3,/(15X,0P,F10.5,1P,4E10.3))')&
     &      ' CARD 3B2***** ',(ANGF(IANG),F(1,IANG,1),F(2,IANG,1),      &
     &      F(3,IANG,1),F(4,IANG,1),IANG=1,NANGLS)
          NWLF=1
          WLF(1)=.55
      ELSE

!         CARD 3C1 (READ ANGLES IN DEGREES)
!         CARD 3C2 (READ WAVELENGTHS IN MICRONS)
!         CARD 3C3 (1ST AEROSOL)
!         CARD 3C3 (READ FOR 1ST ANGLE AND ALL WAVELENGTHS)
!         CARD 3C3 (READ FOR 2ND ANGLE AND ALL WAVELENGTHS)
!         ...
!         CARDS 3C4-3C6 (REPEAT FOR 2ND, 3RD & 4TH AEROSOLS)
!         ...
          IF (TP5_FLAG) THEN
             READ(IRD,'((8F10.0))')(ANGF(IANG),IANG=1,NANGLS)
          ELSE
             !DRF NANGLS MUST BE 50
             NUM_REPEATS = 6
             NUM_REMAINDER = NANGLS-8*NUM_REPEATS
             DO DUMMY=1,NUM_REPEATS
             READ(TP5_FILE(TP5_LINENUMBER),                             &
     &          '((8F10.0))')(ANGF(IANG),                               &
     &            IANG=1+(DUMMY-1)*8,DUMMY*8)
             TP5_LINENUMBER=TP5_LINENUMBER+1
             END DO
             READ(TP5_FILE(TP5_LINENUMBER),                             &
     &          '((2F10.0))')(ANGF(IANG),                               &
     &            IANG=49,50)
             TP5_LINENUMBER=TP5_LINENUMBER+1
          ENDIF
          WRITE(IPR,'(/A,/(8(1X,F9.5)))')                               &
     &      ' CARD 3C1***** ',(ANGF(IANG),IANG=1,NANGLS)
          IF (TP5_FLAG) THEN
             READ(IRD,'((8F10.0))')(WLF(IWAV),IWAV=1,NWLF)
          ELSE
             READ(TP5_FILE(TP5_LINENUMBER),                             &
     &           '((8F10.0))')(WLF(IWAV),IWAV=1,8)
             TP5_LINENUMBER=TP5_LINENUMBER+1
             READ(TP5_FILE(TP5_LINENUMBER),                             &
     &           '((7F10.0))')(WLF(IWAV),IWAV=9,NWLF)
             TP5_LINENUMBER=TP5_LINENUMBER+1
          ENDIF
          WRITE(IPR,'(/A,/(8(1X,F9.4)))')                               &
     &       ' CARD 3C2***** ',(WLF(IWAV),IWAV=1,NWLF)
          DO IAER=1,NAERF
              WRITE(IPR,'(/A,I1,A)')' CARDS 3C',IAER+2,'**** '
              DO IANG=1,NANGLS
                  IF (TP5_FLAG) THEN
                  READ(IRD,'((8F10.0))')(F(IAER,IANG,IWAV),IWAV=1,NWLF)
                  ELSE
                    READ(TP5_FILE(TP5_LINENUMBER),                      &
     &                 '((8F10.0))')(F(IAER,IANG,IWAV),IWAV=1,8)
                    TP5_LINENUMBER=TP5_LINENUMBER+1
                    READ(TP5_FILE(TP5_LINENUMBER),                      &
     &                 '((8F10.0))')(F(IAER,IANG,IWAV),IWAV=9,NWLF)
                    TP5_LINENUMBER=TP5_LINENUMBER+1
                  ENDIF
                  WRITE(IPR,'(1P,(8(1X,E9.3)))')                        &
     &              (F(IAER,IANG,IWAV),IWAV=1,NWLF)
              ENDDO
          ENDDO
      ENDIF

!     RETURN TO CD3A:
      RETURN
      END
