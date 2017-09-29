      SUBROUTINE NOVAER(IHAZE,VIS,ALTNOV,RHNOV,PNOV,TNOV,DENNOV)
      INCLUDE 'PARAMS.h'
      INCLUDE 'ERROR.h'
      INTEGER I,I1,J,IHAZE,NNOV,NWLNOV
      REAL VIS,ALTNOV(MLNOV),RHNOV(MLNOV),TNOV(MLNOV),PNOV(MLNOV)

      REAL EXTNOV(MNOV,MXWVLN),ABSNOV(MNOV,MXWVLN),                     &
     &     ASMNOV(MNOV,MXWVLN),WLNOV(MXWVLN),DENNOV(MNOV,MLNOV),        &
     &     TMPALT(MNOV),TMPRH(MNOV),TMPT(MNOV),TMPP(MNOV)
      LOGICAL LNOVAM

      COMMON/COMNOV/LNOVAM,EXTNOV,ABSNOV,ASMNOV,WLNOV,NNOV,NWLNOV

!     THIS ROUTINE READS THE OUTPUT OF NOVAM AND
!     PUTS THE OUTPUT IN APPROPRIATE ARRAYS.

!     THERE IS A CHANGE FROM USUAL MODTRAN NOMENCLATURE AS FOLLOWS:
!     NOVAM HAS FOUR AEROSOLS.
!     NORMALLY, WE'D HAVE 4 PROFILES; EACH WITH N LB'S (LYR BOUNDARIES)
!     N .LE. LAYDIM WHICH IS SET IN PARAMS.h
!     NOVAM, WE WILL ASSUME, HAS NO MORE THAN 11 LB'S (10 LAYERS)
!     FOR A LOWER ESTIMATE OF THE NUMBER OF VARIABLES, LET N = 11.
!     WE WOULD HAVE 4 ARRAYS OF DIMENSION N OF AEROSOL CONCENTRATIONS.
!     EACH AEROSOL WILL HAVE EXT COEFF AT THE NWVLEN WAVELENGTHS.
!     THE NUMBER OF WAVELENGTH VALUES .LE. MXWVLN (SEE PARAMS.h)
!     FOR THE BUILTIN AEROSOLS, THIS VALUE USED TO BE 47,
!     OLD NWAVLN IN PARAMS.h.
!     FOR NOVAM, WE MAY MAKE NWLNOV SAME AS MXWVLN
!     AT THE MOMENT, IT IS STILL 47.

!     THE UNITLESS OD (OPTICAL DEPTH) = CONC X EXT X PATH LENGTH.
!     UNITLESS BECAUSE (#/L**3)(L**2)(L) IS UNITLESS.
!     THE NUMBER OF VARIABLES FOR PROFILE IS 4*N.
!     THE NUMBER FOR EXTINCTION USED IS 4*47.
!     THE TOTAL IS 4*(11+47)=232 IF N=11 AND NWVLEN=47.
!     NOW THAT NWVLEN IS 788, SO THIS NUMBER CAN BE 3196.

!     BUT NOVAM ALREADY COMPUTES CONC*EXTINCTION PRODUCTS IN 1/KM,
!     ALREADY ADDED UP FOR ALL FOUR AEROSOLS.
!     THERE ARE NWVLEN VALUES OF THIS AT EACH WAVELENGTH FOR EACH LB.
!     TOTAL NUMBER IS NWVLEN*N = 517 IF N=11, NWVLEN=47.
!     TOTAL NUMBER IS NWVLEN*N = 8668 IF N=11, NWVLEN=788.

!     NOTE THAT THE NO. OF AEROSOL LYRS IS REALLY LOT LESS THAN LAYDIM.
!     NOVAM DOES NOT GO ABOVE 2 KM (BOUNDARY LAYER).
!     WE WILL ASSUME THAT NOVAM HAS AT MOST 11 LB'S OR 10 LAYERS.

!     WE WANT NOVAM TO PARALLEL OTHER AEROSOLS IN MODTRAN.
!     MODTRAN METHOD:  PRODUCT OF CONC AND EXT IS IN 1/KM

!     SO WE WILL DO THIS:
!     THE I-TH AEROSOL OCCUPIES ONLY THE I-TH LAYER WHERE CONC IS 1
!     CONC IS ZERO ELSEWHERE
!     I-TH AEROSOL'S EXTINCTION COEFF IS EXT (1/KM) OF THE I-TH LAYER.
!     NOTE THAT PRODUCT OF CONC * EXT IS IMPORTANT AND MUST BE 1/KM.
!     I-TH LAYER'S COEFF = THE AVERAGE OF TOP AND BOTTOM NOVAM VALUE.

!     STILL THERE IS A SLIGHT PROBLEM IN CALCULATING ABSORBER AMOUNTS.
!     MODTRAN WON'T SET TO 0 THE I-TH SPECIE IMMEDIATELY BELOW ITH LYR.
!     SO EACH LB WILL REALLY BE A PAIR OF CLOSELY SPACED (5 M) LB'S.
!     A LAYER WILL BE MODELED BY FOUR LAYER BOUNDARIES AS SHOWN:

!     ----- LB (DENSITY IS 0, TOP LAYER PAIR)
!     ----- LB (DENSITY IS NONZERO, ACTUALLY 1)

!     ----- LB (DENSITY IS NONZERO, ACTUALLY 1, BOTTOM LAYER PAIR)
!     ----- LB (DENSITY IS 0)

!     SO FOR 10 AEROSOLS, WE NEED (10+1)*2=22 LAYER BOUNDARIES

!     IHAZE  = THE TYPE OF BOUNDARY LAYER AEROSOL (EXCLUDING NOVAM)

!     OUTPUT FROM NOVAER:

!     ALTNOV = HEIGHT IN KM
!     RHNOV  = RELATIVE HUMIDITY IN PERCENT
!     TNOV   = TEMPERATURE IN K
!     PNOV   = PRESSURE IN MB
!     ABSNOV = ABSORPTION (NWVLEN VALUES AT EACH LAYER, NOT LB)
!     EXTNOV = EXTINCTION (NWVLEN ...)
!     ASMNOV = ASYMMETRY  (NWVLWN ...)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

      REAL HALF,DELTA,FAC,EXPINT
      INTEGER IA, IB, IMIN1,ITOP,IBOT
      DATA DELTA/0.005/

!*****
!     WRITE CODE TO RUN THE NOVAM DRIVER
!*****

!     THE FOLLOWING CODE IS TO READ NOVAM OUTPUT (NOVAM.OUT)

      OPEN(FILE='NOVAM.OUT',UNIT=60,STATUS='OLD',ERR=99)
      IF(.NOT.LJMASS)                                                   &
     &     WRITE(IPR,'(/,A,A)')' *** NOVAM.OUT FILE IS OPENED ***',     &
     &       ' EXTRA NOVAM ALTITUDES WILL BE INSERTED ***'

!     THE FIRST SET IN NOVAM.OUT IS THE NWVLEN WAVELENGTHS IN MICRONS
      READ(60,*)NWLNOV
      READ(60,*)(WLNOV(I),I=1,NWLNOV)

!     NOW READ THE ALTITUDES IN M AND CONVERT TO KM; ALSO READ RH.
      READ(60,*)NNOV
      IF (NNOV*2.GT.MLNOV) THEN
          IF(LJMASS)CALL WRTBUF(FATAL)
          STOP 'NOVAER:  NOVAM OUTPUT HAS TOO MANY LAYERS'
      ENDIF
      READ(60,*)(TMPALT(I),I=1,NNOV)
      DO 10 I =1, NNOV

!        CONVERT TO METERS
         TMPALT(I)=TMPALT(I)/1000.0
 10   CONTINUE
      READ(60,*)(TMPT(I),I=1,NNOV)
      READ(60,*)(TMPP(I),I=1,NNOV)
      READ(60,*)(TMPRH(I),I=1,NNOV)

!     READ EXT, ABS AND ASYMMETRY PARAMS
      DO 110 J=1, NWLNOV
         READ(60,*)(EXTNOV(I,J),I=1,NNOV)
         DO 100 I=1, NNOV
            IF (EXTNOV(I,J).LT.0.0) THEN
               NWLNOV=J-1
               GO TO 115
            ENDIF
 100     CONTINUE
         READ(60,*)(ABSNOV(I,J),I=1,NNOV)
         READ(60,*)(ASMNOV(I,J),I=1,NNOV)
 110  CONTINUE
 115  CONTINUE
      CLOSE(60)

!     INITIALIZE AND DEFINE DENNOV AND DEFINE ALTS
!     ALSO DEFINE TEMP, P AND RH
!     SINCE DELTA IS SMALL, PROBABLY NO INTERPOLATION IS REQUIRED
!     BUT WE WILL DO IT ANYWAY
      HALF=DELTA/2
      DO 30 I=1, MNOV
         IMIN1=I-1
         IF (I .LE. NNOV) THEN
            J=IMIN1*2+1
            ALTNOV(J)=TMPALT(I)-HALF

!           INTERPOLATE FOR T, P AND RH
            IA=I
            IB=I-1
            IF (IA.EQ.1) IB=2
            FAC=(ALTNOV(J)-TMPALT(IA))/(TMPALT(IB)-TMPALT(IA))
            TNOV(J)=(TMPT(IB)-TMPT(IA))*FAC+TMPT(IA)
            RHNOV(J)=(TMPRH(IB)-TMPRH(IA))*FAC+TMPRH(I)
            PNOV(J)=EXPINT(TMPP(IA),TMPP(IB),FAC)

            J=J+1
            ALTNOV(J)=TMPALT(I)+HALF

!           INTERPOLATE FOR T, P AND RH
            IA=I
            IB=I+1
            IF (IA.EQ.NNOV) IB=I-1
            FAC=(ALTNOV(J)-TMPALT(IA))/(TMPALT(IB)-TMPALT(IA))
            TNOV(J)=(TMPT(IB)-TMPT(IA))*FAC+TMPT(IA)
            RHNOV(J)=(TMPRH(IB)-TMPRH(IA))*FAC+TMPRH(I)
            PNOV(J)=EXPINT(TMPP(IA),TMPP(IB),FAC)
         ENDIF
         IBOT=IMIN1*2+2
         ITOP=IBOT+1
         DO 50 J=1, MLNOV

!           WITH DENNOV, 1ST INDEX IS THE AEROSOL INDEX
!           J IS THE LAYER INDEX
            DENNOV(I,J)=0.0
            IF (J.EQ.ITOP.OR.J.EQ.IBOT)DENNOV(I,J)=1.0
 50      CONTINUE
 30   CONTINUE

!     THE FOLLOWING IS FROM VSA.F - DON'T THINK IT'S NEEDED
      IF(VIS.LE.0.0)THEN
!        DEFAULT FOR VISIBILITY DEPENDS ON THE VALUE OF IHAZE.
         IF(IHAZE.EQ.8)VIS = 0.2
         IF(IHAZE.EQ.9)VIS = 0.5
         IF(IHAZE.EQ.2 .OR. IHAZE.EQ.5)VIS = 5.0
         IF(IHAZE.EQ.1 .OR. IHAZE.EQ.4 .OR. IHAZE.EQ.7)VIS = 23.0
         IF(IHAZE.EQ.6)VIS = 50.0
!        IF(IHAZE.EQ.3)VIS= OR IHAZE = 10 VIS IS DETERMINED ELSEWHERE
      ENDIF

!     CONVERT NNOV TO NUMBER OF NOVAM AEROSOLS
      NNOV=NNOV-1

!     CONVERT TO LAYER VALUES
      DO 130 I=1, NNOV
         I1=I+1
         DO 150 J=1,NWLNOV
            EXTNOV(I,J)=(EXTNOV(I,J)+EXTNOV(I1,J))/2
            ABSNOV(I,J)=(ABSNOV(I,J)+ABSNOV(I1,J))/2
            ASMNOV(I,J)=(ASMNOV(I,J)+ASMNOV(I1,J))/2
 150     CONTINUE
 130  CONTINUE
      RETURN

!     ERROR HANDLING
 99   CONTINUE
      LNOVAM=.FALSE.
      WRITE(IPR,'(//,A,A)')'ALTHOUGH NOVAM AEROSOLS WERE SELECTED, ',   &
     &     '''NOVAM.OUT'' FILE DOES NOT EXIST.'
      WRITE(IPR,'(A,A,/,A)')'CREATE IT BY RUNNING NOVAM OFFLINE,',      &
     &     ' OR DON''T USE NOVAM BY ',                                  &
     &     'REPLACING ''N'' WITH A BLANK IN CARD 2 OF TAPE5.'
      IF(LJMASS)CALL WRTBUF(FATAL)
      STOP
      END
