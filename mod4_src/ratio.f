      REAL FUNCTION  RATIO(A, B)

!               Calculate ratio A/B with over- and under-flow protection

         IF (ABS(A).LT.1.0E-8 .AND. ABS(B).LT.1.0E-8)  THEN
            RATIO=1.0
         ELSE IF (B.EQ.0.0)  THEN
            RATIO=1.E+20
         ELSE
            RATIO=A / B
         END IF

      RETURN
      END
