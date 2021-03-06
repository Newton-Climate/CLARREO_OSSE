
!     /BMHEAD/  BAND MODEL PARAMETER FILE HEADER:
!       IBNDWD   SPECTRAL WIDTH OF BAND MODEL INTERVAL [CM-1].
!       MXFREQ   MAXIMUM FREQUENCY OF BAND MODEL DATA.
!       LSTREC   RECORD NUMBER OF LAST RECORD IN FILE.
!       NTLSUB   NUMBER OF LINE TAIL PARAMETERS PER BAND.
!       EDGENR   NEAR INTERVAL EDGE TO LINE CENTER DISTANCE [CM-1].
!       EDGEFR   FAR INTERVAL EDGE TO LINE CENTER DISTANCE [CM-1].
!       NTEMP    NUMBER OF TEMPERATURES IN BAND MODEL PARAMETER GRID.
!       TBAND    BAND MODEL PARAMETER TEMPERATURE GRID.
      INTEGER IBNDWD,MXFREQ,LSTREC,NTLSUB,NTEMP
      REAL EDGENR,EDGEFR,TBAND
      COMMON/BMHEAD/IBNDWD,MXFREQ,LSTREC,NTLSUB,EDGENR,EDGEFR,          &
     &  NTEMP,TBAND(MXTEMP)
