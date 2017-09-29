      SUBROUTINE FRN296(V1C,FH2O)
      REAL V1C,FH2O

!     LOADS FOREIGN CONTINUUM 296K
      INTEGER NPT
      REAL V1,DV,F296
      COMMON/FH2O/V1,DV,NPT,F296(2003)

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      SAVE /FH2O/
      EXTERNAL BFH2O
      CALL SINT(V1,V1C,DV,NPT,F296,FH2O)
      RETURN
      END