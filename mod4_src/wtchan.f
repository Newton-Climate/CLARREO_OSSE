      LOGICAL FUNCTION WTCHAN(LSTBIN,SPECMN,WGTMN,SPECMX,WGTMX)

!     THIS FUNCTION DETERMINES CHANNEL CONTRIBUTIONS (WTCHN) TO
!     SPECTRAL BINS.  A VALUE OF FALSE IS RETURNED IF THE ROUTINE
!     IS UNSUCCESSFUL.  THE CONTRIBUTIONS ARE CALCULATED FREQUENCY
!     INTEGRALS FROM V0 TO V1 OVER A TABULATED WEIGHTING FUNCTIONS:

!                   / V1
!         WTCHN  =  |      W(V)  DV
!                   / V0

!     IF SPECMN AND SPECMX ARE IN WAVENUMBERS [CM-1],

!                          V - SPECMN
!        W(V) = WGTMN + --------------- (WGTMX - WGTMN)   .
!                       SPECMX - SPECMN

!     IF SPECMN AND SPECMX ARE IN WAVELENGTH UNITS,

!                       CONVRT/V - SPECMX
!        W(V) = WGTMX + ----------------- (WGTMX - WGTMN)   .
!                        SPECMN - SPECMX

!     WHERE CONVRT IS THE CONVERSION FACTOR BETWEEN THE WAVELENGTH
!     UNIT AND WAVENUMBERS.

!     ARGUMENTS:
!       LSTBIN   LAST WAVENUMBER SPECTRAL BIN FOR WHICH A CHANNEL
!                CONTRIBUTION WAS CALCULATED [=-1 FOR NEW CHANNELS].
!       SPECMN   SPECTRAL GRID MINIMUM VALUE [UNITS FROM UNTFLG].
!       WGTMN    FILTER RESPONSE WEIGHT AT SPECMN.
!       SPECMX   SPECTRAL GRID MAXIMUM VALUE [UNITS FROM UNTFLG].
!       WGTMX    FILTER RESPONSE WEIGHT AT SPECMX.
        INTEGER LSTBIN
        REAL SPECMN,WGTMN,SPECMX,WGTMX

!     PARAMETERS:
      INCLUDE 'PARAMS.h'
      INCLUDE 'CHANLS.h'

!     COMMONS:
      INCLUDE 'IFIL.h'
      INCLUDE 'BMHEAD.h'

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      INTEGER NBINMN,NBINMX,NBIN
      REAL WGT,FREQMN,FREQMX,DFREQ,FREQHI,DWGT,SLOPEW,SLOPEF,           &
     &  WGT0,WGT0D,SLOPE0,CONVRT,DWAVE,DRATIO,EXPAND,X

!     INTERNAL FUNCTION
      EXPAND(X)=.5+X*(.33333333+.25*X)

!     INITIALIZE WTCHAN TO .FALSE. AND CHECK INPUTS:
      WTCHAN=.FALSE.
      IF(SPECMN.GT.SPECMX)THEN
          WRITE(IPR,'(/1X,A)')                                          &
     &      'WARNING:  SPECMX is less than SPECMN in function WTCHAN'
          RETURN
      ELSEIF(WGTMN.LT.0. .OR. WGTMX.LT.0.)THEN
          WRITE(IPR,'(/1X,A)')                                          &
     &      'WARNING:  Negative weight input to function WTCHAN'
          RETURN
      ENDIF

!     BRANCH BASED ON UNTFLG:
      IF(UNTFLG.EQ.'W')THEN

!         FREQUENCY GRID [CM-1] - DETERMINE MINIMUM
!         AND MAXIMUM SPECTRAL WAVENUMBER BINS:
          NBINMN=INT(SPECMN/IBNDWD+.5)
          NBINMX=INT(SPECMX/IBNDWD+.5)
          IF(NBINMN.LT.MNBIN)THEN
              WRITE(IPR,'(/2A,I6,A,/10X,A,I6)')' WARNING: ',            &
     &          ' Parameter MNBIN (=',MNBIN,') must be decreased.',     &
     &          ' Current NBINMN in function WTCHAN equals',NBINMN
              RETURN
          ELSEIF(NBINMX.GT.MXBIN)THEN
              WRITE(IPR,'(/2A,I6,A,/10X,A,I6)')' WARNING: ',            &
     &          ' Parameter MXBIN (=',MXBIN,') must be increased.',     &
     &          ' Current NBINMX in function WTCHAN equals',NBINMX
              RETURN
          ENDIF

!         INCREMENT INTEGRATED WIDTH OF CHANNEL:
          WDFREQ(NCHAN)=WDFREQ(NCHAN)+.5*(WGTMX+WGTMN)*(SPECMX-SPECMN)
          SLOPEF=(WGTMX-WGTMN)/(SPECMX-SPECMN)
          IF(SPECMN.GT.0.)THEN
              WDWAVE(NCHAN)=WDWAVE(NCHAN)+10000.                        &
     &          *(WGTMN/SPECMN-WGTMX/SPECMX+SLOPEF*LOG(SPECMX/SPECMN))
          ELSE
              WDWAVE(NCHAN)=BIGNUM
          ENDIF

!         DETERMINE CONTRIBUTIONS TO EACH SPECTRAL BIN.
          IF(NBINMN.EQ.NBINMX)THEN

!             ONE BIN ONLY; CALCULATE WEIGHT.
              WGT=.5*(WGTMN+WGTMX)*(SPECMX-SPECMN)
              IF(NBINMN.EQ.LSTBIN)THEN

!                 ADD CURRENT WEIGHT TO OLD VALUE:
                  WTCHN(NUMCHN(NBINMN),NBINMN)=                         &
     &              WTCHN(NUMCHN(NBINMN),NBINMN)+WGT
              ELSEIF(NUMCHN(NBINMN).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINMN CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ELSE

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBINMN)=NUMCHN(NBINMN)+1
                  LSTCHN(NUMCHN(NBINMN),NBINMN)=NCHAN
                  WTCHN(NUMCHN(NBINMN),NBINMN)=WGT
              ENDIF
          ELSE

!             MULTIPLE BINS; DEFINE FIRST WEIGHT.
              DFREQ=IBNDWD*(NBINMN+.5)-SPECMN
              WGT=DFREQ*(WGTMN+.5*DFREQ*SLOPEF)
              IF(NBINMN.EQ.LSTBIN)THEN

!                 ADD CURRENT WEIGHT TO OLD VALUE:
                  WTCHN(NUMCHN(NBINMN),NBINMN)=                         &
     &              WTCHN(NUMCHN(NBINMN),NBINMN)+WGT
              ELSEIF(NUMCHN(NBINMN).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINMN CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ELSE

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBINMN)=NUMCHN(NBINMN)+1
                  LSTCHN(NUMCHN(NBINMN),NBINMN)=NCHAN
                  WTCHN(NUMCHN(NBINMN),NBINMN)=WGT
              ENDIF

!             LOOP OVER INTERMEDIATE BINS:
              SLOPE0=IBNDWD*SLOPEF
              WGT0=IBNDWD*WGTMN
              DO 10 NBIN=NBINMN+1,NBINMX-1
                  IF(NUMCHN(NBIN).GE.MXNCHN)THEN

!                     SPECTRAL BIN NBIN CONTRIBUTES TO TOO MANY CHANNELS
                      WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',      &
     &                  ' Parameter MXNCHN must be increased. ',        &
     &                  ' Single spectral',' bins contribute',          &
     &                  ' to more than',MXNCHN,' filter channels.'
                      RETURN
                  ENDIF

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBIN)=NUMCHN(NBIN)+1
                  LSTCHN(NUMCHN(NBIN),NBIN)=NCHAN
                  WTCHN(NUMCHN(NBIN),NBIN)                              &
     &              =WGT0+(IBNDWD*NBIN-SPECMN)*SLOPE0
   10         CONTINUE

!             LAST BIN:  INCREMENT COUNTER, LIST CHANNEL
!             AND DEFINE WEIGHT AND VARIABLE LSTBIN.
              IF(NUMCHN(NBINMX).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINMX CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ENDIF
              DFREQ=SPECMX-(NBINMX-.5)*IBNDWD
              NUMCHN(NBINMX)=NUMCHN(NBINMX)+1
              LSTCHN(NUMCHN(NBINMX),NBINMX)=NCHAN
              WTCHN(NUMCHN(NBINMX),NBINMX)=DFREQ*(WGTMX-.5*DFREQ*SLOPEF)
          ENDIF

!         SET LSTBIN TO THE MAXIMUM BIN FOR FREQUENCY FILTERS.
          LSTBIN=NBINMX
      ELSE

!         WAVELENGTH GRID - DETERMINE MINIMUM
!         AND MAXIMUM FREQUENCIES [CM-1]:
          IF(UNTFLG.EQ.'N')THEN

!             NANOMETERS:
              CONVRT=1.E7
          ELSEIF(UNTFLG.EQ.'M')THEN

!             MICRONS:
              CONVRT=1.E4
          ELSE
              WRITE(IPR,'(/2A,2(/11X,A),2A)')' WARNING:  ',             &
     &          'Input UNTFLG to routine WTCHAN must equal "W" for',    &
     &          'WAVENUMBERS, "M" for Microns, or "N" for Nanometers.', &
     &          'Currently, UNTFLG = "',UNTFLG,'"'
              RETURN
          ENDIF
          FREQMN=CONVRT/SPECMX
          FREQMX=CONVRT/SPECMN

!         DETERMINE MINIMUM AND MAXIMUM SPECTRAL WAVENUMBER BINS:
          NBINMN=INT(FREQMN/IBNDWD+.5)
          NBINMX=INT(FREQMX/IBNDWD+.5)
          IF(NBINMN.LT.MNBIN)THEN
              WRITE(IPR,'(/2A,I6,A,/10X,A,I6)')' WARNING: ',            &
     &          ' Parameter MNBIN (=',MNBIN,') must be decreased.',     &
     &          ' Current NBINMN in function WTCHAN equals',NBINMN
              RETURN
          ELSEIF(NBINMX.GT.MXBIN)THEN
              WRITE(IPR,'(/2A,I6,A,/10X,A,I6)')' WARNING: ',            &
     &          ' Parameter MXBIN (=',MXBIN,') must be increased.',     &
     &          ' Current NBINMX in function WTCHAN equals',NBINMX
              RETURN
          ENDIF

!         INCREMENT INTEGRATED WIDTH OF CHANNEL:
          WDWAVE(NCHAN)=WDWAVE(NCHAN)+.5*(WGTMX+WGTMN)*(SPECMX-SPECMN)

!         THE WEIGHTING FUNCTION SLOPE:
          DWAVE=SPECMX-SPECMN
          DWGT=WGTMX-WGTMN
          SLOPEW=DWGT/DWAVE
          SLOPEF=CONVRT*SLOPEW
          IF(NBINMN.EQ.NBINMX)THEN

!             ONE BIN ONLY; CALCULATE WEIGHT.
              DRATIO=DWAVE/SPECMX
              IF(DRATIO.GT..01)THEN
                  WGT=(FREQMX*WGTMN-FREQMN*WGTMX)-SLOPEF*LOG(1.-DRATIO)
              ELSE
                  WGT=DRATIO*(WGTMN*FREQMX+DWGT*FREQMN*EXPAND(DRATIO))
              ENDIF
              IF(NBINMN.EQ.LSTBIN)THEN

!                 ADD CURRENT WEIGHT TO OLD VALUE:
                  WTCHN(NUMCHN(NBINMN),NBINMN)=                         &
     &              WTCHN(NUMCHN(NBINMN),NBINMN)+WGT
              ELSEIF(NUMCHN(NBINMN).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINMN CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ELSE

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBINMN)=NUMCHN(NBINMN)+1
                  LSTCHN(NUMCHN(NBINMN),NBINMN)=NCHAN
                  WTCHN(NUMCHN(NBINMN),NBINMN)=WGT
              ENDIF
              WDFREQ(NCHAN)=WDFREQ(NCHAN)+WGT
          ELSE

!             MULTIPLE BINS; DEFINE FIRST WEIGHT.
              FREQHI=(NBINMN+.5)*IBNDWD
              DFREQ=FREQHI-FREQMN
              DRATIO=DFREQ/FREQMN
              WGT0=(WGTMN*SPECMX-WGTMX*SPECMN)/DWAVE
              IF(DRATIO.GT..01)THEN
                  WGT=DFREQ*WGT0+SLOPEF*LOG(1.+DRATIO)
              ELSE
                  WGT=DFREQ*(WGTMX-SPECMX*SLOPEW*DRATIO*EXPAND(-DRATIO))
              ENDIF
              IF(NUMCHN(NBINMN).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINMN CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ENDIF

!             INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
              NUMCHN(NBINMN)=NUMCHN(NBINMN)+1
              LSTCHN(NUMCHN(NBINMN),NBINMN)=NCHAN
              WTCHN(NUMCHN(NBINMN),NBINMN)=WGT
              WDFREQ(NCHAN)=WDFREQ(NCHAN)+WGT

!             LOOP OVER INTERMEDIATE BINS:
              WGT0D=IBNDWD*WGT0
              DO 20 NBIN=NBINMN+1,NBINMX-1
                  IF(NUMCHN(NBIN).GE.MXNCHN)THEN

!                     SPECTRAL BIN NBIN CONTRIBUTES TO TOO MANY CHANNELS
                      WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',      &
     &                  ' Parameter MXNCHN must be increased. ',        &
     &                  ' Single spectral',' bins contribute',          &
     &                  ' to more than',MXNCHN,' filter channels.'
                      RETURN
                  ENDIF

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBIN)=NUMCHN(NBIN)+1
                  LSTCHN(NUMCHN(NBIN),NBIN)=NCHAN
                  FREQHI=(NBIN+.5)*IBNDWD
                  DRATIO=IBNDWD/FREQHI
                  IF(DRATIO.GT..01)THEN
                      WGT=WGT0D-SLOPEF*LOG(1.-DRATIO)
                  ELSE
                      WGT=WGT0D+SLOPEF*DRATIO*(1.+DRATIO*EXPAND(DRATIO))
                  ENDIF
                  WTCHN(NUMCHN(NBIN),NBIN)=WGT
                  WDFREQ(NCHAN)=WDFREQ(NCHAN)+WGT
   20         CONTINUE

!             LAST BIN:  INCREMENT COUNTER, LIST CHANNEL
!             AND DEFINE WEIGHT AND VARIABLE LSTBIN.
              DFREQ=FREQMX-FREQHI
              DRATIO=DFREQ/FREQMX
              IF(DRATIO.GT..01)THEN
                  WGT=DFREQ*WGT0-SLOPEF*LOG(1.-DRATIO)
              ELSE
                  WGT=DFREQ*(WGTMN+SPECMN*SLOPEW*DRATIO*EXPAND(DRATIO))
              ENDIF
              IF(NBINMX.EQ.LSTBIN)THEN

!                 ADD CURRENT WEIGHT TO OLD VALUE:
                  WTCHN(NUMCHN(NBINMX),NBINMX)=                         &
     &              WTCHN(NUMCHN(NBINMX),NBINMX)+WGT
              ELSEIF(NUMCHN(NBINMX).GE.MXNCHN)THEN

!                 SPECTRAL BIN NBINMX CONTRIBUTES TO TOO MANY CHANNELS.
                  WRITE(IPR,'(/3A,/10X,2A,I4,A)')' WARNING: ',          &
     &              ' Parameter MXNCHN must be increased. ',            &
     &              ' Single spectral',' bins contribute',              &
     &              ' to more than',MXNCHN,' filter channels.'
                  RETURN
              ELSE

!                 INCREMENT COUNTER, LIST CHANNEL AND DEFINE WEIGHT.
                  NUMCHN(NBINMX)=NUMCHN(NBINMX)+1
                  LSTCHN(NUMCHN(NBINMX),NBINMX)=NCHAN
                  WTCHN(NUMCHN(NBINMX),NBINMX)=WGT
              ENDIF
              WDFREQ(NCHAN)=WDFREQ(NCHAN)+WGT
          ENDIF

!         SET LSTBIN TO THE MINIMUM BIN FOR WAVELENGTH FILTERS.
          LSTBIN=NBINMN
      ENDIF

!     CHANGE WTCHAN TO .TRUE. AND RETURN:
      WTCHAN=.TRUE.
      RETURN
      END
