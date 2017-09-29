      SUBROUTINE WRTFLT(IEMSCT,IVXMIN,IVXMAX)

!     THIS ROUTINE PERFORMS THE FILTER RESPONSE SUMS FOR EACH CHANNEL.

!     ARGUMENTS:
!       IEMSCT   MODTRAN RADIATION TRANSPORT MODE FLAG
!                  0 FOR TRANSMITTANCE ONLY
!                  1 FOR THERMAL RADIANCE
!                  2 FOR THERMAL + SOLAR RADIANCE
!                  3 FOR TRANSMITTED SOLAR IRRADIANCE
!       IVXMIN   MINIMUM CALCULATION FREQUENCY [CM-1].
!       IVXMAX   MAXIMUM CALCULATION FREQUENCY [CM-1].
      INTEGER IEMSCT,IVXMIN,IVXMAX

!     PARAMETERS:
      INCLUDE 'CHANLS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     FUNCTIONS:
      INTEGER NUNIT

!     LOCAL VARIABLES:
!       WLUNIT   WAVELENGTH UNIT FOR SPECTRAL RADIANCE OR IRRADIANCE.
!       UNITS    UNITS USED FOR WIDTH AND SPECTRAL MINIMUM AND MAXIMUM.
!       ICHAN    CHANNEL INDEX.
!       LEXIST   FILE EXIST FLAG.
!       LOPEN    FILE OPEN FLAG.
      CHARACTER WLUNIT*12,UNITS*44
      INTEGER ICHAN
      LOGICAL LEXIST,LOPEN

!     DATA
!       ZERO     THE NUMBER 0.
      REAL ZERO
      DATA ZERO/0./

!     CHECK UNTFLG AND DEFINE UNITS:
      IF(UNTFLG.EQ.'M')THEN
          WLUNIT='(PER MICRON)'
          UNITS='    (CM-1)   (MICRONS)  (MICRONS)  (MICRONS)'
      ELSEIF(UNTFLG.EQ.'N')THEN
          WLUNIT='    (PER NM)'
          UNITS='    (CM-1)      (NM)       (NM)       (NM)  '
      ELSEIF(UNTFLG.EQ.'W')THEN
          WLUNIT='(PER MICRON)'
          UNITS='    (CM-1)   (MICRONS)    (CM-1)     (CM-1) '
      ELSE
          WRITE(IPR,'(/4A,/10X,A)')' WARNING: ',                        &
     &      ' No write to file "',CHNOUT(1:LENCHN),'" because UNTFLG',  &
     &      ' in "CHANLS.h" does not equal "M", "N" or "W".'
          RETURN
      ENDIF

!     OPEN CHANNEL OUTPUT FILE IF NOT ALREADY OPENED.
      INQUIRE(FILE=CHNOUT,EXIST=LEXIST,OPENED=LOPEN)
      IF(.NOT.LEXIST)THEN
          ICHNUN=NUNIT()
          OPEN(ICHNUN,FILE=CHNOUT,STATUS='NEW')
      ELSEIF(.NOT.LOPEN)THEN
          WRITE(IPR,'(/3A)')                                            &
     &      ' WARNING:  Deleting old "',CHNOUT(1:LENCHN),'" file.'
          ICHNUN=NUNIT()
          OPEN(ICHNUN,FILE=CHNOUT,STATUS='UNKNOWN')
          CLOSE(ICHNUN,STATUS='DELETE')
          OPEN(ICHNUN,FILE=CHNOUT,STATUS='NEW')
      ENDIF

!     HEADER AND RESULTS:
!       IF NPR=-1, PRINT ALL CHANNELS;
!       IF NPR= 0, ONLY PRINT CHANNELS WITH A POSITIVE VLCHAN; AND
!       IF NPR= 1, ONLY PRINT FULLY INTEGRATED CHANNELS.
      IF(IEMSCT.EQ.0)THEN
          IF(.NOT.LJMASS) WRITE(ICHNUN,'(/(2A))')                       &
     &      ' CHAN      AVERAGE        CHANNEL  ',                      &
     &      '      FULL CHANNEL      SPECTRAL   SPECTRAL   CHANNEL',    &
     &      '  NEL    EXTINCTION     EXTINCTION ',                      &
     &      '    EQUIVALENT WIDTH     MINIMUM    MAXIMUM   DESCRIPTION',&
     &      '  NO.   (1 - TRANS)       (CM-1)   ',UNITS,                &
     &      ' ----  -------------  -------------',                      &
     &      '  ---------  ---------  ---------  ---------  -----------'
          DO 10 ICHAN=1,NCHAN
              IF(VLCHAN(ICHAN).GT.0. .OR. NPR.LT.0)THEN
                  IF(NFRQLO(ICHAN).LT.IVXMIN .OR.                       &
     &              NFRQHI(ICHAN).GT.IVXMAX)THEN
                      IF(NPR.LE.0)                                      &
     &                  WRITE(ICHNUN,'(I5,2F15.7,4F11.4,2X,A)')         &
     &                  -ICHAN,ZERO,VLCHAN(ICHAN),                      &
     &                  WDFREQ(ICHAN),WDWAVE(ICHAN),SPECLO(ICHAN),      &
     &                  SPECHI(ICHAN),NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ELSE
                      WRITE(ICHNUN,'(I5,2F15.7,4F11.4,2X,A)')ICHAN,     &
     &                  VLCHAN(ICHAN)/WDFREQ(ICHAN),VLCHAN(ICHAN),      &
     &                  WDFREQ(ICHAN),WDWAVE(ICHAN),SPECLO(ICHAN),      &
     &                  SPECHI(ICHAN),NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ENDIF
              ENDIF
              VLCHAN(ICHAN)=0.
   10     CONTINUE
      ELSEIF(IEMSCT.EQ.3)THEN
          WRITE(ICHNUN,'(2(/2A),/(4A))')                                &
     &      ' CHAN    TRANSMITTED SPECTRAL SOLAR   TRANSMITTED ',       &
     &      '      FULL CHANNEL      SPECTRAL   SPECTRAL   CHANNEL',    &
     &      '  NEL    IRRADIANCE (W CM-2 / XXXX)   SOLAR IRRAD.',       &
     &      '    EQUIVALENT WIDTH     MINIMUM    MAXIMUM   DESCRIPTION',&
     &      '  NO.     (PER CM-1)   ',WLUNIT,'     (W CM-2)  ',UNITS,   &
     &      ' ----  -------------  -------------  -------------',       &
     &      '  ---------  ---------  ---------  ---------  -----------'
          DO 20 ICHAN=1,NCHAN
              IF(VLCHAN(ICHAN).GT.0. .OR. NPR.LT.0)THEN
                  IF(NFRQLO(ICHAN).LT.IVXMIN .OR.                       &
     &              NFRQHI(ICHAN).GT.IVXMAX)THEN
                      IF(NPR.LE.0)                                      &
     &                  WRITE(ICHNUN,'(I5,1P3E15.6,0P4F11.4,2X,A)')     &
     &                  -ICHAN,ZERO,ZERO,VLCHAN(ICHAN),                 &
     &                  WDFREQ(ICHAN),WDWAVE(ICHAN),SPECLO(ICHAN),      &
     &                  SPECHI(ICHAN),NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ELSE
                      WRITE(ICHNUN,'(I5,1P3E15.6,0P4F11.4,2X,A)')       &
     &                  ICHAN,VLCHAN(ICHAN)/WDFREQ(ICHAN),              &
     &                  VLCHAN(ICHAN)/WDWAVE(ICHAN),VLCHAN(ICHAN),      &
     &                  WDFREQ(ICHAN),WDWAVE(ICHAN),SPECLO(ICHAN),      &
     &                  SPECHI(ICHAN),NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ENDIF
              ENDIF
              VLCHAN(ICHAN)=0.
   20     CONTINUE
      ELSE
          IF(.NOT.LJMASS) WRITE(ICHNUN,'(2(/2A),/(4A))')                &
     &      ' CHAN        SPECTRAL  RADIANCE          CHANNEL  ',       &
     &      '      FULL CHANNEL      SPECTRAL   SPECTRAL   CHANNEL',    &
     &      '  NEL       (W SR-1 CM-2 / XXXX)        RADIANCE  ',       &
     &      '    EQUIVALENT WIDTH     MINIMUM    MAXIMUM   DESCRIPTION',&
     &      '  NO.     (PER CM-1)   ',WLUNIT,'  (W SR-1 CM-2)',UNITS,   &
     &      ' ----  -------------  -------------  -------------',       &
     &      '  ---------  ---------  ---------  ---------  -----------'
          DO 30 ICHAN=1,NCHAN
              IF(VLCHAN(ICHAN).GT.0. .OR. NPR.LT.0)THEN
                  IF(NFRQLO(ICHAN).LT.IVXMIN .OR.                       &
     &              NFRQHI(ICHAN).GT.IVXMAX)THEN
                      IF(NPR.LE.0)                                      &
     &                  WRITE(ICHNUN,'(I5,1P3E15.6,0P4F11.4,2X,A)')     &
     &                  -ICHAN,ZERO,ZERO,VLCHAN(ICHAN),                 &
     &                  WDFREQ(ICHAN),WDWAVE(ICHAN),SPECLO(ICHAN),      &
     &                  SPECHI(ICHAN),NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ELSE
                      WRITE(ICHNUN,'(I5,1P3E15.6,0P4F11.4,2X,A)')       &
     &                  ICHAN,VLCHAN(ICHAN)/WDFREQ(ICHAN),              &
     &                  VLCHAN(ICHAN)/WDWAVE(ICHAN),VLCHAN(ICHAN),      &
     &                  WDFREQ(ICHAN),WDWAVE(ICHAN),SPECLO(ICHAN),      &
     &                  SPECHI(ICHAN),NMCHAN(ICHAN)(1:LNCHAN(ICHAN))
                  ENDIF
              ENDIF
              VLCHAN(ICHAN)=0.
   30     CONTINUE
      ENDIF

!     TABLE COMPLETE:
      RETURN
      END
