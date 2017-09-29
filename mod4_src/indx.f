      SUBROUTINE INDX(WAVL,TC,KEY,REIL,ZIMAG)
      REAL AB,DOP,R1,R2,REIL,TC,WAVL,ZIMAG
      INTEGER KEY
!*** END OF DECLARATIONS INSERTED BY SPAG
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * *
! * * WAVELENGTH IS IN CENTIMETERS.  TEMPERATURE IS IN DEG. C.      * *
! * *
! * * KEY IS SET TO 1 IN SUBROUTINE GAMFOG                          * *
! * *
! * * KEY IS SET TO 0 IN SUBROUTINE GAMFOG    FOR CIRRUS            * *
! * *
! * * REAL IS THE REAL PART OF THE REFRACTIVE INDEX.                * *
! * *
! * * ZIMAG IS THE IMAGINARY PART OF THE REFRACTIVE INDEX IT IS     * *
! * *
! * * RETURNED NEG. I.E.  M= REAL - I*ZIMAG  .                      * *
! * *
! * * A SERIES OF CHECKS ARE MADE AND WARNINGS GIVEN.               * *
! * *
! * * RAY APPLIED OPTICS VOL 11,NO.8,AUG 72, PG. 1836-1844          * *
! * *
! * * CORRECTIONS HAVE BEEN MADE TO RAYS ORIGINAL PAPER             * *
! * *
! * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD
      R1 = 0.0
      R2 = 0.0
      IF(WAVL.LT..0001)WRITE(IPR,'(//2A,/)')' Warning from INDX: ',     &
     &    ' Attempting to evaluate for a wavelength below 1 micron.'
      IF(TC.LT.-20.)WRITE(IPR,'(//2A,/)')' Warning from INDX: ',        &
     &    ' Attempting to evaluate for a temperature below -20 deg C.'
      CALL DEBYE(WAVL,TC,KEY,REIL,ZIMAG)

! * *  TABLE 3 WATER PG. 1840
      IF(WAVL.GT..034)THEN
         IF(WAVL.LE..1)THEN
            R2 = DOP(WAVL,1.83899,1639.,52340.4,10399.2,588.24,345005., &
     &           259913.,161.29,43319.7,27661.2)
            R2 = R2+R2*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)
            REIL = REIL*(WAVL-.034)/.066+R2*(.1-WAVL)/.066
         ENDIF
      ELSEIF(WAVL.GT..0006)THEN
         REIL = DOP(WAVL,1.83899,1639.,52340.4,10399.2,588.24,345005.,  &
     &          259913.,161.29,43319.7,27661.2)
         REIL = REIL+REIL*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)
         IF(WAVL.LE..0007)THEN
            R1 = DOP(WAVL,1.79907,3352.27,99.914E+04,15.1963E+04,1639., &
     &           50483.5,9246.27,588.24,84.4697E+04,10.7615E+05)
            R1 = R1+R1*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)
            REIL = R1*(.0007-WAVL)/.0001+REIL*(WAVL-.0006)/.0001
         ENDIF
      ELSE
         REIL = DOP(WAVL,1.79907,3352.27,99.914E+04,15.1963E+04,1639.,  &
     &          50483.5,9246.27,588.24,84.4697E+04,10.7615E+05)
         REIL = REIL+REIL*(TC-25.)*.0001*EXP((.000025*WAVL)**.25)
      ENDIF

! * *  TABLE 2 WATER PG. 1840
      IF(WAVL.LT..3)THEN
         IF(WAVL.GE..03)THEN
            ZIMAG = ZIMAG+AB(WAVL,.25,300.,.47,3.)                      &
     &              +AB(WAVL,.39,17.,.45,1.3)+AB(WAVL,.41,62.,.35,1.7)
         ELSEIF(WAVL.GE..0062)THEN
            ZIMAG = ZIMAG+AB(WAVL,.41,62.,.35,1.7)                      &
     &              +AB(WAVL,.39,17.,.45,1.3)+AB(WAVL,.25,300.,.4,2.)
         ELSEIF(WAVL.GE..0017)THEN
            ZIMAG = ZIMAG+AB(WAVL,.39,17.,.45,1.3)                      &
     &              +AB(WAVL,.41,62.,.22,1.8)+AB(WAVL,.25,300.,.4,2.)
         ELSEIF(WAVL.GE..00061)THEN
            ZIMAG = ZIMAG+AB(WAVL,.12,6.1,.042,.6)                      &
     &              +AB(WAVL,.39,17.,.165,2.4)+AB(WAVL,.41,62.,.22,1.8)
         ELSEIF(WAVL.GE..000495)THEN
            ZIMAG = ZIMAG+AB(WAVL,.01,4.95,.05,1.)                      &
     &              +AB(WAVL,.12,6.1,.009,2.)
         ELSEIF(WAVL.GE..000297)THEN
            ZIMAG = ZIMAG+AB(WAVL,.27,2.97,.04,2.)                      &
     &              +AB(WAVL,.01,4.95,.06,1.)
         ELSE
            ZIMAG = ZIMAG+AB(WAVL,.27,2.97,.025,2.)                     &
     &              +AB(WAVL,.01,4.95,.06,1.)
         ENDIF
      ENDIF
      RETURN
      END
