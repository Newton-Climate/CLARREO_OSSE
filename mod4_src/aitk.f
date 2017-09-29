      FUNCTION AITK(ARG,VAL,X,NDIM)

!     IBM SCIENTIFIC FUNCTION: AITKEN INTERPOLATION (ASSUMES NDIM>0).

!     ARGUMENTS:
      INTEGER NDIM
      REAL X,ARG(NDIM),VAL(NDIM)

!     START OF AITKEN-LOOP
      DO J=2,NDIM
          IEND=J-1
          DO I=1,IEND
              IF(ARG(I).EQ.ARG(J))THEN

!                 TWO IDENTICAL ARGUMENT VALUES IN VECTOR ARG
                  AITK=VAL(IEND)
                  RETURN
              ENDIF
              VAL(J)=VAL(J)+(VAL(I)-VAL(J))*(X-ARG(J))/(ARG(I)-ARG(J))
          ENDDO
      ENDDO
      AITK=VAL(NDIM)
      RETURN
      END
