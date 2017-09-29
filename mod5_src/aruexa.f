      SUBROUTINE ARUEXA(TP5_FLAG,TP5_FILE,TP5_LINENUMBER)

!     READ IN USER DEFINED EXTINCTION, ABSORPTION COEFFICIENTS AND
!     ASYMMETRY PARAMETERS FOR THE USS OPTION
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

!     COMMONS:
      INCLUDE 'BASE.h'
      INCLUDE 'IFIL.h'
      LOGICAL LARUSS
      INTEGER NARSPC
      REAL VARSPC
      COMMON/USSPC/NARSPC(4),VARSPC(4,NWAVLN),LARUSS
      DOUBLE PRECISION ALTB
      INTEGER IREG,IREGC
      COMMON/CARD2D/IREG(4),ALTB(4),IREGC(4)

!     LOCAL VARIABLES
      CHARACTER AERNAM(18)*4
      INTEGER ISPC,ISPCM1,ISPC3,LIMSPC,IHC,IK,IOS

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     ANNOUNCE READING OF AEROSOL SPECTRAL DATA
      IF(.NOT.LJMASS)WRITE(IPR,'(/A)')                                  &
     &  ' USER-SUPPLIED SPECTRAL DATA (USS ENHANCEMENT)'

!     FOR THE USUAL AEROSOLS (NON-USS CASE)
!     THE 4 IREG'S ARE 1 OR 0 AND READ IF IHAZE=7 OR ICLD=11.
!     BUT FOR THE USS SCHEME, IREG(IK)=0 HAS THE SAME MEANING.
!     THAT IS, IREG(IK)=0 MEANS USE BUILTIN MODEL SPECTRAL DATA

!     BUT IREG(IK) = INTEGER MEANS READ USER-DEFINED SPECTRAL DATA.
!     AFTER THE DATA IS READ, IREG(IK) WILL BE RESET TO 1.

      IF(LJMASS)THEN
          CALL INITCD('CARD2D')
      ELSE
          IF (TP5_FLAG) THEN
            READ(IRD,'(4I5)')(IREG(IK),IK=1,4)
          ELSE
            READ(TP5_FILE(TP5_LINENUMBER),'(4I5)')(IREG(IK),IK=1,4)
            TP5_LINENUMBER=TP5_LINENUMBER+1
          ENDIF
          WRITE(IPR,'(A,4I5,A)')' CARD 2D *****',(IREG(IK),IK=1,4),     &
     &      ' (NUMBER OF USER-DEFINED SPECTRAL POINTS FOR AEROSOLS)'
      ENDIF

      DO IHC=1,4
          NARSPC(IHC)=IREG(IHC)

!         CHECK THAT NARSPC IS NOT TOO LARGE
          IF(NARSPC(IHC).GT.MXWVLN)THEN
              WRITE(IPR,'(/A,I3,A,/26X,A,I3,A)')                        &
     &          ' Error in routine ARUEXA:  Input number of aerosol'//  &
     &          ' spectral data points (NARSPC =',NARSPC(IHC),          &
     &          ') exceeds',' the maximum number (Parameter MXWVLN =',  &
     &          MXWVLN, ' in file PARAMS.h).'
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' ARUEXA: NUMBER OF SPECTRAL POINTS IS TOO LARGE.'
          ENDIF

          IF(NARSPC(IHC).LT.0)THEN
              IF(LJMASS)CALL WRTBUF(FATAL)
              STOP ' ARUEXA: NUMBER OF SPECTRAL DATA MUST BE POSITIVE.'
          ENDIF
          IF(NARSPC(IHC).NE.0)THEN
              IF(LJMASS)THEN
                  CALL INITCD('CARD2D1')
              ELSE
                  IF (TP5_FLAG) THEN
                    READ(IRD,'(F10.0,18A4)')AWCCON(IHC),AERNAM
                  ELSE
                    READ(TP5_FILE(TP5_LINENUMBER),                      &
     &                '(F10.0,18A4)')AWCCON(IHC),AERNAM
                    TP5_LINENUMBER=TP5_LINENUMBER+1
                  ENDIF
                  WRITE(IPR,'(/A,1P,E10.3,18A4)')' CARD 2D1 ****'       &
     &              //' EQUIVALENT WATER = ',AWCCON(IHC),AERNAM
                  WRITE(IPR,'(/A)')' CARD 2D2 ****'
                  ISPCM1=0
                  DO ISPC3=1,NARSPC(IHC),3
                      LIMSPC=MIN(NARSPC(IHC),ISPC3+2)
                      IF (TP5_FLAG) THEN
                      READ(IRD,'(3(F6.2,2F7.5,F6.4))',IOSTAT=IOS)       &
     &                  (VARSPC(IHC,ISPC),EXTC(IHC,ISPC),ABSC(IHC,ISPC),&
     &                  ASYM(IHC,ISPC),ISPC=ISPC3,LIMSPC)
                      ELSE
                      READ(TP5_FILE(TP5_LINENUMBER),                    &
     &                  '(3(F6.2,2F7.5,F6.4))',IOSTAT=IOS)              &
     &                  (VARSPC(IHC,ISPC),EXTC(IHC,ISPC),ABSC(IHC,ISPC),&
     &                  ASYM(IHC,ISPC),ISPC=ISPC3,LIMSPC)
                        TP5_LINENUMBER=TP5_LINENUMBER+1
                      ENDIF
                      IF(IOS.NE.0)THEN
                          WRITE(IPR,'(/A,I3,A)')'Error in ARUEXA:  '//  &
     &                    'Unable to read line',1+ISPC3/3,' of CARD 2D2'
                          STOP'Error in ARUEXA: Unable to read CARD 2D2'
                      ENDIF
                      WRITE(IPR,'(3(F6.2,1X,2(F7.5,1X),F6.4,1X))')      &
     &                  (VARSPC(IHC,ISPC),EXTC(IHC,ISPC),ABSC(IHC,ISPC),&
     &                  ASYM(IHC,ISPC),ISPC=ISPC3,LIMSPC)
                      DO ISPC=ISPC3,LIMSPC
                          IF(ISPC.EQ.1)THEN
                              IF(VARSPC(IHC,1).LT.0.)THEN
                                  WRITE(IPR,'(/A)')                     &
     &                              'Error in ARUEXA:  Negative entry'  &
     &                              //' in CARD 2D2 spectral grid.'
                                  STOP'Error in ARUEXA: Wavelength < 0.'
                              ENDIF
                          ELSEIF(VARSPC(IHC,ISPC)                       &
     &                       .LE.VARSPC(IHC,ISPCM1))THEN
                              WRITE(IPR,'(/A)')'Error in ARUEXA: '//    &
     &                          ' Decreasing CARD 2D2 spectral grid.'
                              STOP 'Decreasing CARD 2D2 spectral grid.'
                          ENDIF
                          IF(EXTC(IHC,ISPC).LT.0.)THEN
                              WRITE(IPR,'(/A)')'Error in ARUEXA:  Neg'//&
     &                          'ative CARD 2D2 extinction coefficient.'
                              STOP 'CARD 2D2 extinction coef < 0.'
                          ELSEIF(ABSC(IHC,ISPC).LT.0.)THEN
                              WRITE(IPR,'(/A)')'Error in ARUEXA:  Neg'//&
     &                          'ative CARD 2D2 absorption coefficient.'
                              STOP 'CARD 2D2 absorption coef < 0.'
                          ELSEIF(ABSC(IHC,ISPC).GT.EXTC(IHC,ISPC))THEN
                              WRITE(IPR,'(/A)')'Error in ARUEXA:  CARD' &
     &                          //' 2D2 absorption exceeds extinction.'
                              STOP 'CARD 2D2 absorption > extinction.'
                          ELSEIF(ABS(ASYM(IHC,ISPC)).GE.1.)THEN
                              WRITE(IPR,'(/A)')'Error in ARUEXA:  As'// &
     &                          'ymmetry factor magnitude exceeds one.'
                              STOP 'CARD 2D2 |asymmetry factor| > 1.'
                          ENDIF
                          ISPCM1=ISPC
                      ENDDO
                  ENDDO
              ENDIF

!JMASS        AERNAM NOT USED IN JMASS, OTHER VARIABLES AVAILABLE IN
!             COMMON BLOCK, NWAVLN DEFINED IN PARAMS.h

!             SET NON-ZERO IREG VALUES TO 1 TO TRIGGER USER-DEFINED
!             OPTION.  THIS IS UNNECESSARY BUT DOES NOT HURT
!             THE USER-DEFINED VALUE IS USED FOR ANY NON-ZER0 IREG(IHC)
!             AND IREG(IHC) IS NON-ZERO HERE ANYWAY
              IREG(IHC)=1
          ENDIF
      ENDDO
      RETURN
      END
