      BLOCK DATA CRD5BD

!     COMMONS:

!     /JM5/
!       IRPT     REPEAT INPUT FLAG (0=NONE, 1=ALL, 3=GEOM, 4=SPEC).
!       IFAC     CURRENT COLUMN SCALING FACTOR INDEX.
!       NFACMN   NUMBER OF COLUMN SCALING FACTOR LESS THAN 1.
!       NFACMX   NUMBER OF COLUMN SCALING FACTOR GREATER THAN 1.
!       FACMC    CURRENT COLUMN SCALING FACTOR.
!       SCALMN   MINIMUM COLUMN SCALING FACTOR.
!       SCALMX   MAXIMUM COLUMN SCALING FACTOR.
      INTEGER IRPT,IFAC,NFACMN,NFACMX
      REAL FACMC
      DOUBLE PRECISION SCALMN,SCALMX
      COMMON/JM5/SCALMN,SCALMX,IRPT,IFAC,NFACMN,NFACMX,FACMC

!     /CJM5/
!       AMOD3D   FLAG INDICATING OUTPUT DATABASE FILE TYPE:
      CHARACTER AMOD3D*1
      COMMON/CJM5/AMOD3D

!     DATA:
      DATA IRPT,IFAC,NFACMN,NFACMX/4*0/,FACMC/1./,SCALMN,SCALMX/2*1.D0/
      DATA AMOD3D/' '/
      END
      SUBROUTINE CARD5(LNFLRT,FLRT)

!     ROUTINE TO READ IN AFTER PROCESS CARD5 INPUTS:

!     INPUT ARGUMENTS:
!       LNFLRT   LENGTH OF ROOT NAME FOR ALL I/O FILES.
!       FLRT     ROOT NAME FOR ALL I/O FILES.
      INTEGER LNFLRT
      CHARACTER FLRT*(*)

!     PARAMETERS:
!       ONE      THE NUMBER ONE IN DOUBLE PRECISION.
      DOUBLE PRECISION ONE
      PARAMETER(ONE=1.D0)

!     COMMONS:
      INCLUDE 'IFIL.h'

!     /JM5/
!       IRPT     REPEAT INPUT FLAG (0=NONE, 1=ALL, 3=GEOM, 4=SPEC).
!       IFAC     CURRENT COLUMN SCALING FACTOR INDEX.
!       NFACMN   NUMBER OF COLUMN SCALING FACTOR LESS THAN 1.
!       NFACMX   NUMBER OF COLUMN SCALING FACTOR GREATER THAN 1.
!       SCALMN   MINIMUM COLUMN SCALING FACTOR.
!       FACMC    CURRENT COLUMN SCALING FACTOR.
!       SCALMX   MAXIMUM COLUMN SCALING FACTOR.
      INTEGER IRPT,IFAC,NFACMN,NFACMX
      REAL FACMC
      DOUBLE PRECISION SCALMN,SCALMX
      COMMON/JM5/SCALMN,SCALMX,IRPT,IFAC,NFACMN,NFACMX,FACMC

!     /CJM5/
!       AMOD3D   FLAG INDICATING OUTPUT DATABASE FILE TYPE:
      CHARACTER AMOD3D*1
      COMMON/CJM5/AMOD3D

!     /CNTRL/
!       IKMAX    NUMBER OF PATH SEGMENTS ALONG LINE-OF-SIGHT.
!       ML       NUMBER OF ATMOSPHERIC PROFILE LEVELS.
!       MLFLX    NUMBER OF LEVELS FOR WHICH FLUX VALUES ARE WRITTEN.
!       ISSGEO   LINE-OF-SIGHT FLAG (0 = SENSOR PATH, 1 = SOLAR PATHS).
!       IMULT    MULTIPLE SCATTERING FLAG
!                  (0=NONE, 1=AT SENSOR, -1=AT FINAL OR TANGENT POINT).
      INTEGER IKMAX,ML,MLFLX,ISSGEO,IMULT
      COMMON/CNTRL/IKMAX,ML,MLFLX,ISSGEO,IMULT

!     FUNCTIONS:
!       IRECLN   RETURNS THE RECORD LENGTH FOR A DIRECT ACCESS FILES
!                CONTAINING A KNOWN NUMBER OF VARIABLES IN EACH RECORD.
      INTEGER IRECLN

!     LOCAL VARIABLES:
!       IOS      RESULT OF IOSTAT CHECK.
!       LRECLN   RECORD LENGTH OF DIRECT ACCESS FILES.
!       EXTEN    FILE EXTENSION NAME.
!       CHAR80   INPUT STRING.
!       FULNAM   FULL PATH NAME.
      INTEGER IOS,LRECLN
      CHARACTER EXTEN*8,CHAR80*80,FULNAM*88

!     SAVED VARIABLES:
!       SNUMER   NUMERATOR FACTOR USED TO DEFINE COLUMN SCALING.
!       SDENOM   DENOMINATOR FACTOR USED TO DEFINE COLUMN SCALING.
      DOUBLE PRECISION SNUMER,SDENOM
      SAVE SNUMER,SDENOM

!     CHECK FOR NEW READ:
      IF(IFAC.LT.NFACMX)THEN
          IFAC=IFAC+1
      ELSE
!DRF          READ(IRD,'(A80)')CHAR80
	  CHAR80 = '    0                                               '
          READ(CHAR80,'(I5,50X,A1,2(I3,F9.0))',IOSTAT=IOS)              &
     &      IRPT,AMOD3D,NFACMN,SCALMN,NFACMX,SCALMX
          IF(IOS.NE.0)THEN

!             OLD FORMAT:
              READ(CHAR80(1:5),'(I5)')IRPT
              AMOD3D=' '
              FACMC=1.
              IFAC=0
          ELSEIF(AMOD3D.EQ.'T' .OR. AMOD3D.EQ.'t')THEN

!             MOD3D DATABASE:
              AMOD3D='T'
              FACMC=1.
              NFACMX=0
              IFAC=0
          ELSEIF(AMOD3D.EQ.'C' .OR. AMOD3D.EQ.'c')THEN

!             MC CONTINUA DATABASE:
              AMOD3D='C'
              CLOSE(IDBOUT,STATUS='DELETE')
              LRECLN=IRECLN(10*IKMAX,LNFLRT,FLRT)
              IF(LNFLRT.LE.0)THEN
                  OPEN(IDBOUT,FILE='mc.bin',FORM='UNFORMATTED',         &
     &              ACCESS='DIRECT',RECL=LRECLN)
              ELSE
                  FULNAM=FLRT(1:LNFLRT)//'.mcb'
                  OPEN(IDBOUT,FILE=FULNAM,FORM='UNFORMATTED',           &
     &              ACCESS='DIRECT',RECL=LRECLN)
              ENDIF
              NFACMX=0
              IFAC=0
          ELSEIF(AMOD3D.EQ.'M' .OR. AMOD3D.EQ.'m')THEN

!             MC MOLECULAR TRANSMITTANCE DATABASE:
              AMOD3D='M'
              IFAC=-NFACMN
              IF(NFACMN.LE.0 .OR. SCALMN.GE.ONE .OR.                    &
     &           NFACMX.LE.0 .OR. SCALMX.LE.ONE)STOP 'Bad CARD5 input'
              SDENOM=NFACMN*NFACMX*(SCALMX-SCALMN)
              SNUMER=(SCALMX*(ONE-SCALMN)*NFACMX                        &
     &               -SCALMN*(SCALMX-ONE)*NFACMN)/SDENOM
              SDENOM=((ONE-SCALMN)*NFACMX-(SCALMX-ONE)*NFACMN)/SDENOM
              WRITE(IDBOUT,'(/(I13,A,/F13.4,A))')                       &
     &          NFACMN,' NUMBER OF SCALE FACTORS LESS THAN 1.',         &
     &          SCALMN,' MINIMUM SCALE FACTOR',                         &
     &          NFACMX,' NUMBER OF SCALE FACTORS GREATER THAN 1.',      &
     &          SCALMX,' MAXIMUM SCALE FACTOR'
          ELSE
              AMOD3D=' '
              FACMC=1.
              IFAC=0
          ENDIF
      ENDIF
      IF(AMOD3D.EQ.'M')THEN
          FACMC=REAL((ONE+IFAC*SNUMER)/(ONE+IFAC*SDENOM))
          REWIND(IPR)
          REWIND(IPU)
          CLOSE(IDBOUT)
          EXTEN(1:2)='.M'
          WRITE(EXTEN(3:8),'(I6.6)')INT(100*FACMC+.5)
          LRECLN=IRECLN(6*IKMAX+6,LNFLRT,FLRT)
          IF(LNFLRT.LE.0)THEN
              OPEN(IDBOUT,FILE=EXTEN(2:8),STATUS='UNKNOWN')
              CLOSE(IDBOUT,STATUS='DELETE')
              OPEN(IDBOUT,FILE=EXTEN(2:8),STATUS='NEW',ACCESS='DIRECT', &
     &          FORM='UNFORMATTED',RECL=LRECLN)
          ELSE
              FULNAM=FLRT(1:LNFLRT)//EXTEN
              OPEN(IDBOUT,FILE=FULNAM,STATUS='UNKNOWN')
              CLOSE(IDBOUT,STATUS='DELETE')
              OPEN(IDBOUT,FILE=FULNAM,STATUS='NEW',ACCESS='DIRECT',     &
     &          FORM='UNFORMATTED',RECL=LRECLN)
          ENDIF
      ENDIF
      RETURN
      END