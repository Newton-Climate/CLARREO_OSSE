      SUBROUTINE  ASYMTX(AA, EVEC, EVAL, M, IA, IEVEC, IER, WK)

!    =======  S I N G L E    P R E C I S I O N    V E R S I O N  ======

!       SOLVES EIGENFUNCTION PROBLEM FOR REAL ASYMMETRIC MATRIX
!       FOR WHICH IT IS KNOWN A PRIORI THAT THE EIGENVALUES ARE REAL.

!       THIS IS AN ADAPTATION OF A SUBROUTINE EIGRF IN THE IMSL
!       LIBRARY TO USE REAL INSTEAD OF COMPLEX ARITHMETIC, ACCOUNTING
!       FOR THE KNOWN FACT THAT THE EIGENVALUES AND EIGENVECTORS IN
!       THE DISCRETE ORDINATE SOLUTION ARE REAL.  OTHER CHANGES INCLUDE
!       PUTTING ALL THE CALLED SUBROUTINES IN-LINE, DELETING THE
!       PERFORMANCE INDEX CALCULATION, UPDATING MANY DO-LOOPS
!       TO FORTRAN77, AND IN CALCULATING THE MACHINE PRECISION
!       DRIGHT INSTEAD OF SPECIFYING IT IN A DATA STATEMENT.

!       EIGRF IS BASED PRIMARILY ON EISPACK ROUTINES.  THE MATRIX IS
!       FIRST BALANCED USING THE PARLETT-REINSCH ALGORITHM.  THEN
!       THE MARTIN-WILKINSON ALGORITHM IS APPLIED.

!       REFERENCES:
!          DONGARRA, J. AND C. MOLER, EISPACK -- A PACKAGE FOR SOLVING
!             MATRIX EIGENVALUE PROBLEMS, IN COWELL, ED., 1984:
!             SOURCES AND DEVELOPMENT OF MATHEMATICAL SOFTWARE,
!             PRENTICE-HALL, ENGLEWOOD CLIFFS, NJ
!         PARLETT AND REINSCH, 1969: BALANCING A MATRIX FOR CALCULATION
!             OF EIGENVALUES AND EIGENVECTORS, NUM. MATH. 13, 293-304
!         WILKINSON, J., 1965: THE ALGEBRAIC EIGENVALUE PROBLEM,
!             CLARENDON PRESS, OXFORD

!   I N P U T    V A R I A B L E S:

!        A    :  INPUT ASYMMETRIC MATRIX, DESTROYED AFTER SOLVED
!        M    :  ORDER OF  A
!       IA    :  FIRST DIMENSION OF  A
!    IEVEC    :  FIRST DIMENSION OF  EVEC

!   O U T P U T    V A R I A B L E S:

!       EVEC  :  (UNNORMALIZED) EIGENVECTORS OF  A
!                   ( COLUMN J CORRESPONDS TO EVAL(J) )

!       EVAL  :  (UNORDERED) EIGENVALUES OF  A ( DIMENSION AT LEAST M )

!       IER   :  IF .NE. 0, SIGNALS THAT EVAL(IER) FAILED TO CONVERGE;
!                   IN THAT CASE EIGENVALUES IER+1,IER+2,...,M  ARE
!                   CORRECT BUT EIGENVALUES 1,...,IER ARE SET TO ZERO.

!   S C R A T C H   V A R I A B L E S:

!       AAD   :  DOUBLE PRECISION STAND-IN FOR -A-
!       EVECD :  DOUBLE PRECISION STAND-IN FOR -EVEC-
!       EVALD :  DOUBLE PRECISION STAND-IN FOR -EVAL-
!       WKD   :  DOUBLE PRECISION STAND-IN FOR -WK-
!+---------------------------------------------------------------------+

      REAL              AA(IA,*),  WK(*),  EVAL(*), EVEC(IEVEC,*)
      LOGICAL           NOCONV, NOTLAS

!     COMMON/RCNSTN/
!       PI       THE CONSTANT PI.
!       DEG      NUMBER OF DEGREES IN ONE RADIAN.
!       BIGNUM   MAXIMUM SINGLE PRECISION NUMBER.
!       BIGEXP   MAXIMUM EXPONENTIAL ARGUMENT WITHOUT OVERFLOW.
!       RRIGHT   SMALLEST SINGLE PRECISION REAL ADDED TO 1 EXCEEDS 1.
      REAL PI,DEG,BIGNUM,BIGEXP,RRIGHT
      COMMON/RCNSTN/PI,DEG,BIGNUM,BIGEXP,RRIGHT
      DATA     C1/ 0.4375 /, C2/ 0.5 /, C3/ 0.75 /, C4/ 0.95 /,         &
     &         C5/ 16.0 /, C6/ 256.0 /, ZERO / 0.0 /, ONE / 1.0 /

      IER=0
      IF (M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M)                          &
     &     CALL ERRMSG('ASYMTX--BAD INPUT VARIABLE(S)', .TRUE.)

!                           ** HANDLE 1X1 AND 2X2 SPECIAL CASES
      IF (M.EQ.1)  THEN
         EVAL(1)=AA(1,1)
         EVEC(1,1)=ONE
         RETURN

      ELSE IF (M.EQ.2)  THEN
         DISCRI=(AA(1,1) - AA(2,2))**2 + 4. * AA(1,2) * AA(2,1)
         IF (DISCRI.LT.ZERO)                                            &
     &        CALL ERRMSG('ASYMTX--COMPLEX EVALS IN 2X2 CASE', .TRUE.)
         SGN=ONE
         IF (AA(1,1).LT.AA(2,2))  SGN=- ONE
         EVAL(1)=0.5 * (AA(1,1) + AA(2,2) + SGN*SQRT(DISCRI))
         EVAL(2)=0.5 * (AA(1,1) + AA(2,2) - SGN*SQRT(DISCRI))
         EVEC(1,1)=ONE
         EVEC(2,2)=ONE
         IF (AA(1,1).EQ.AA(2,2) .AND.                                   &
     &        (AA(2,1).EQ.ZERO.OR.AA(1,2).EQ.ZERO))  THEN
            RNORM=  ABS(AA(1,1)) + ABS(AA(1,2)) + ABS(AA(2,1))          &
     &              + ABS(AA(2,2))
            W=RRIGHT * RNORM
            EVEC(2,1)=  AA(2,1) / W
            EVEC(1,2)=- AA(1,2) / W
         ELSE
            EVEC(2,1)=AA(2,1) / (EVAL(1) - AA(2,2))
            EVEC(1,2)=AA(1,2) / (EVAL(2) - AA(1,1))
         ENDIF

         RETURN

      END IF
!                                        ** INITIALIZE OUTPUT VARIABLES
      IER=0
      DO 20 I=1, M
         EVAL(I)=ZERO
         DO 10 J=1, M
            EVEC(I,J)=ZERO
10       CONTINUE
         EVEC(I,I)=ONE
20    CONTINUE
!                  ** BALANCE THE INPUT MATRIX AND REDUCE ITS NORM BY
!                  ** DIAGONAL SIMILARITY TRANSFORMATION STORED IN WK;
!                  ** THEN SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                  ** AND PUSH THEM DOWN
      RNORM=ZERO
      L =1
      K =M

30    KKK=K
         DO 70  J=KKK, 1, -1
            ROW=ZERO
            DO 40 I=1, K
               IF (I.NE.J) ROW=ROW + ABS(AA(J,I))
40          CONTINUE
            IF (ROW.EQ.ZERO) THEN
               WK(K)=J
               IF (J.NE.K) THEN
                  DO 50 I=1, K
                     REPL  =AA(I,J)
                     AA(I,J)=AA(I,K)
                     AA(I,K)=REPL
50                CONTINUE
                  DO 60 I=L, M
                     REPL  =AA(J,I)
                     AA(J,I)=AA(K,I)
                     AA(K,I)=REPL
60                CONTINUE
               END IF
               K=K - 1
               GO TO 30
            END IF
70       CONTINUE
!                                     ** SEARCH FOR COLUMNS ISOLATING AN
!                                       ** EIGENVALUE AND PUSH THEM LEFT
80    LLL=L
         DO 120 J=LLL, K
            COL=ZERO
            DO 90 I=L, K
               IF (I.NE.J) COL=COL + ABS(AA(I,J))
90          CONTINUE
            IF (COL.EQ.ZERO) THEN
               WK(L)=J
               IF (J.NE.L) THEN
                  DO 100 I=1, K
                     REPL  =AA(I,J)
                     AA(I,J)=AA(I,L)
                     AA(I,L)=REPL
100               CONTINUE
                  DO 110 I=L, M
                     REPL  =AA(J,I)
                     AA(J,I)=AA(L,I)
                     AA(L,I)=REPL
110               CONTINUE
               END IF
               L=L + 1
               GO TO 80
            END IF
120      CONTINUE
!                           ** BALANCE THE SUBMATRIX IN ROWS L THROUGH K
      DO 130 I=L, K
         WK(I)=ONE
130   CONTINUE

140   NOCONV=.FALSE.
         DO 200 I=L, K
            COL=ZERO
            ROW=ZERO
            DO 150 J=L, K
               IF (J.NE.I) THEN
                  COL=COL + ABS(AA(J,I))
                  ROW=ROW + ABS(AA(I,J))
               END IF
150         CONTINUE
            F=ONE
            G=ROW / C5
            H=COL + ROW
160         IF (COL.LT.G) THEN
               F  =F * C5
               COL=COL * C6
               GO TO 160
            END IF
            G=ROW * C5
170         IF (COL.GE.G) THEN
               F  =F / C5
               COL=COL / C6
               GO TO 170
            END IF
!                                                         ** NOW BALANCE
            IF ((COL+ROW) / F .LT. C4 * H) THEN
               WK(I) =WK(I) * F
               NOCONV=.TRUE.
               DO 180 J=L, M
                  AA(I,J)=AA(I,J) / F
180            CONTINUE
               DO 190 J=1, K
                  AA(J,I)=AA(J,I) * F
190            CONTINUE
            END IF
200      CONTINUE

      IF (NOCONV) GO TO 140
!                                  ** IS -A- ALREADY IN HESSENBERG FORM?
      IF (K-1.LT.L+1) GO TO 350
!                                   ** TRANSFER -A- TO A HESSENBERG FORM
      DO 290 N=L+1, K-1
         H      =ZERO
         WK(N+M)=ZERO
         SCALE  =ZERO
!                                                        ** SCALE COLUMN
         DO 210 I=N, K
            SCALE=SCALE + ABS(AA(I,N-1))
210      CONTINUE
         IF (SCALE.NE.ZERO) THEN
            DO 220 I=K, N, -1
               WK(I+M)=AA(I,N-1) / SCALE
               H=H + WK(I+M) * WK(I+M)
220         CONTINUE
            G=- SIGN(SQRT(H),WK(N+M))
            H=H - WK(N+M) * G
            WK(N+M)=WK(N+M) - G
!                                                 ** FORM (I-(U*UT)/H)*A
            DO 250 J=N, M
               F=ZERO
               DO 230  I=K, N, -1
                  F=F + WK(I+M) * AA(I,J)
230            CONTINUE
               DO 240 I=N, K
                  AA(I,J)=AA(I,J) - WK(I+M) * F / H
240            CONTINUE
250         CONTINUE
!                                    ** FORM (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 280 I=1, K
               F=ZERO
               DO 260  J=K, N, -1
                  F=F + WK(J+M) * AA(I,J)
260            CONTINUE
               DO 270 J=N, K
                  AA(I,J)=AA(I,J) - WK(J+M) * F / H
270            CONTINUE
280         CONTINUE
            WK(N+M) =SCALE * WK(N+M)
            AA(N,N-1)=SCALE * G
         END IF
290   CONTINUE

      DO 340  N=K-2, L, -1
         N1=N + 1
         N2=N + 2
         F =AA(N+1,N)
         IF (F.NE.ZERO) THEN
!mod
!mod        Protect against F*WK(N+1+M) being 0.
!mod        F =F * WK(N+1+M)
            F1 = WK(N+1+M)
            DO 300 I=N+2, K
               WK(I+M)=AA(I,N)
300         CONTINUE
            IF (N+1.LE.K) THEN
               DO 330 J=1, M
                  G=ZERO
                  DO 310 I=N+1, K
                     G=G + WK(I+M) * EVEC(I,J)
310               CONTINUE
                  G=G / F
!mod
!mod              Incorporate F1 into the denominator now.
                  G=G / F1
                  DO 320 I=N+1, K
                     EVEC(I,J)=EVEC(I,J) + G * WK(I+M)
320               CONTINUE
330            CONTINUE
            END IF
         END IF
340   CONTINUE

350   CONTINUE
      N=1
      DO 370 I=1, M
         DO 360 J=N, M
            RNORM=RNORM + ABS(AA(I,J))
360      CONTINUE
         N=I
         IF (I.LT.L .OR. I.GT.K) EVAL(I)=AA(I,I)
370   CONTINUE
      N=K
      T=ZERO
!                                         ** SEARCH FOR NEXT EIGENVALUES
380   IF (N.LT.L) GO TO 530
      IN=0
      N1=N - 1
      N2=N - 2
!                          ** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
390   CONTINUE
      DO 400 I=L, N
         LB=N+L - I
         IF (LB.EQ.L) GO TO 410
         S=ABS(AA(LB-1,LB-1)) + ABS(AA(LB,LB))
         IF (S.EQ.ZERO) S=RNORM
         IF (ABS(AA(LB,LB-1)) .LE. RRIGHT * S) GO TO 410
400   CONTINUE
!                                                      ** ONE EVAL FOUND
410   X=AA(N,N)
      IF (LB.EQ.N ) THEN
         AA(N,N) =X + T
         EVAL(N)=AA(N,N)
         N=N1
         GO TO 380
      END IF
!                                                     ** TWO EVALS FOUND
      Y=AA(N1,N1)
      W=AA(N,N1) * AA(N1,N)
      IF (LB.EQ.N1) THEN
         P=(Y-X) * C2
         Q=P * P + W
         Z=SQRT(ABS(Q))
         AA(N,N)=X + T
         X=AA(N,N)
         AA(N1,N1)=Y + T
!                                                           ** REAL PAIR
         Z=P + SIGN(Z,P)
         EVAL(N1)=X + Z
         EVAL(N) =EVAL(N1)
         IF (Z.NE.ZERO) EVAL(N)=X - W / Z
         X=AA(N,N1)
!                                  ** EMPLOY SCALE FACTOR IN CASE
!                                  ** X AND Z ARE VERY SMALL
         R=SQRT(X * X + Z * Z)
         P=X / R
         Q=Z / R
!                                                    ** ROW MODIFICATION
         DO 420 J=N1, M
            Z=AA(N1,J)
            AA(N1,J)=Q * Z + P * AA(N,J)
            AA(N,J) =Q * AA(N,J) - P * Z
420      CONTINUE
!                                                 ** COLUMN MODIFICATION
         DO 430 I=1, N
            Z=AA(I,N1)
            AA(I,N1)=Q * Z + P * AA(I,N)
            AA(I,N) =Q * AA(I,N) - P * Z
430      CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
         DO 440 I=L, K
            Z=EVEC(I,N1)
            EVEC(I,N1)=Q * Z + P * EVEC(I,N)
            EVEC(I,N) =Q * EVEC(I,N) - P * Z
440      CONTINUE
         N=N2
         GO TO 380
      END IF
!                    ** NO CONVERGENCE AFTER 30 ITERATIONS; SET ERROR
!                    ** INDICATOR TO THE INDEX OF THE CURRENT EIGENVALUE

      IF (IN.EQ.30) THEN
         IER=128 + N
         RETURN
      END IF
!                                                          ** FORM SHIFT
      IF (IN.EQ.10 .OR. IN.EQ.20) THEN
         T=T + X
         DO 450 I=L, N
            AA(I,I)=AA(I,I) - X
450      CONTINUE
         S=ABS(AA(N,N1)) + ABS(AA(N1,N2))
         X=C3 * S
         Y=X
         W=-C1 * S * S
      END IF

      IN=IN + 1
!                ** LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS

      DO 460 J=LB, N2
         I=N2+LB - J
         Z=AA(I,I)
         R=X - Z
         S=Y - Z
         P=(R * S-W) / AA(I+1,I) + AA(I,I+1)
         Q=AA(I+1,I+1) - Z - R - S
         R=AA(I+2,I+1)
         S=ABS(P) + ABS(Q) + ABS(R)
         P=P / S
         Q=Q / S
         R=R / S
         IF (I.EQ.LB) GO TO 470
         UU=ABS(AA(I,I-1)) * (ABS(Q) + ABS(R))
         VV=ABS(P) * (ABS(AA(I-1,I-1)) + ABS(Z) + ABS(AA(I+1,I+1)))
         IF (UU .LE. RRIGHT*VV) GO TO 470
460   CONTINUE

470   CONTINUE
      AA(I+2,I)=ZERO
      DO 480 J=I+3, N
         AA(J,J-2)=ZERO
         AA(J,J-3)=ZERO
480   CONTINUE

!             ** DOUBLE QR STEP INVOLVING ROWS K TO N AND COLUMNS M TO N

      DO 520 KA=I, N1
         NOTLAS=KA.NE.N1
         IF (KA.EQ.I) THEN
            S=SIGN(SQRT(P*P + Q*Q + R*R),P)
            IF (LB.NE.I) AA(KA,KA-1)=- AA(KA,KA-1)
         ELSE
            P=AA(KA,KA-1)
            Q=AA(KA+1,KA-1)
            R=ZERO
            IF (NOTLAS) R=AA(KA+2,KA-1)
            X=ABS(P) + ABS(Q) + ABS(R)
            IF (X.EQ.ZERO) GO TO 520
            P=P / X
            Q=Q / X
            R=R / X
            S=SIGN(SQRT(P*P + Q*Q + R*R),P)
            AA(KA,KA-1)=- S * X
         END IF
         P=P + S
         X=P / S
         Y=Q / S
         Z=R / S
         Q=Q / P
         R=R / P
!                                                    ** ROW MODIFICATION
         DO 490 J=KA, M
            P=AA(KA,J) + Q * AA(KA+1,J)
            IF (NOTLAS) THEN
               P=P + R * AA(KA+2,J)
               AA(KA+2,J)=AA(KA+2,J) - P * Z
            END IF
            AA(KA+1,J)=AA(KA+1,J) - P * Y
            AA(KA,J)  =AA(KA,J)   - P * X
490      CONTINUE
!                                                 ** COLUMN MODIFICATION
         DO 500 II=1, MIN0(N,KA+3)
            P=X * AA(II,KA) + Y * AA(II,KA+1)
            IF (NOTLAS) THEN
               P=P + Z * AA(II,KA+2)
               AA(II,KA+2)=AA(II,KA+2) - P * R
            END IF
            AA(II,KA+1)=AA(II,KA+1) - P * Q
            AA(II,KA)  =AA(II,KA) - P
500      CONTINUE
!                                          ** ACCUMULATE TRANSFORMATIONS
         DO 510 II=L, K
            P=X * EVEC(II,KA) + Y * EVEC(II,KA+1)
            IF (NOTLAS) THEN
               P=P + Z * EVEC(II,KA+2)
               EVEC(II,KA+2)=EVEC(II,KA+2) - P * R
            END IF
            EVEC(II,KA+1)=EVEC(II,KA+1) - P * Q
            EVEC(II,KA)  =EVEC(II,KA) - P
510      CONTINUE
520   CONTINUE
      GO TO 390
!                     ** ALL EVALS FOUND, NOW BACKSUBSTITUTE REAL VECTOR
530   CONTINUE
      IF (RNORM.NE.ZERO) THEN
         DO 560  N=M, 1, -1
            N2=N
            AA(N,N)=ONE
            DO 550  I=N-1, 1, -1
               W=AA(I,I) - EVAL(N)
               IF (W.EQ.ZERO) W=RRIGHT * RNORM
               R=AA(I,N)
               DO 540 J=N2, N-1
                  R=R + AA(I,J) * AA(J,N)
540            CONTINUE
               AA(I,N)=-R / W
               N2=I
550         CONTINUE
560      CONTINUE
!                      ** END BACKSUBSTITUTION VECTORS OF ISOLATED EVALS

         DO 580 I=1, M
            IF (I.LT.L .OR. I.GT.K) THEN
               DO 570 J=I, M
                  EVEC(I,J)=AA(I,J)
570            CONTINUE
            END IF
580      CONTINUE
!                                   ** MULTIPLY BY TRANSFORMATION MATRIX
         IF (K.NE.0) THEN
            DO 610  J=M, L, -1
               DO 600 I=L, K
                  Z=ZERO
                  DO 590 N=L, MIN0(J,K)
                     Z=Z + EVEC(I,N) * AA(N,J)
590               CONTINUE
                  EVEC(I,J)=Z
600            CONTINUE
610         CONTINUE
         END IF

      END IF

      DO 620 I=L, K
         DO 620 J=1, M
            EVEC(I,J)=EVEC(I,J) * WK(I)
620   CONTINUE
!                           ** INTERCHANGE ROWS IF PERMUTATIONS OCCURRED
      DO 640  I=L-1, 1, -1
         J=INT(WK(I))
         IF (I.NE.J) THEN
            DO 630 N=1, M
               REPL     =EVEC(I,N)
               EVEC(I,N)=EVEC(J,N)
               EVEC(J,N)=REPL
630         CONTINUE
         END IF
640   CONTINUE

      DO 660 I=K+1, M
         J=INT(WK(I))
         IF (I.NE.J) THEN
            DO 650 N=1, M
               REPL     =EVEC(I,N)
               EVEC(I,N)=EVEC(J,N)
               EVEC(J,N)=REPL
650         CONTINUE
         END IF
660   CONTINUE

      RETURN
      END
