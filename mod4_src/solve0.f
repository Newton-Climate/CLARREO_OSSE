      SUBROUTINE  SOLVE0(B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,  &
     &                    FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT,       &
     &                    MAZIM, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR,&
     &                    NNLYRI, PI, TPLANK, TAUCPR, UMU0, Z, ZZ,      &
     &                    ZPLK0, ZPLK1)

!        CONSTRUCT RIGHT-HAND SIDE VECTOR -B- FOR GENERAL BOUNDARY
!        CONDITIONS STWJ(17) AND SOLVE SYSTEM OF EQUATIONS OBTAINED
!        FROM THE BOUNDARY CONDITIONS AND THE
!        CONTINUITY-OF-INTENSITY-AT-LAYER-INTERFACE EQUATIONS.
!        THERMAL EMISSION CONTRIBUTES ONLY IN AZIMUTHAL INDEPENDENCE.

!     ROUTINES CALLED:  SGBCO, SGBSL, ZEROIT

!     I N P U T      V A R I A B L E S:

!       BDR      :  SURFACE BIDIRECTIONAL REFLECTIVITY
!       BEM      :  SURFACE BIDIRECTIONAL EMISSIVITY
!       BPLANK   :  BOTTOM BOUNDARY THERMAL EMISSION
!       CBAND    :  LEFT-HAND SIDE MATRIX OF LINEAR SYSTEM EQ. SC(5),
!                   SCALED BY EQ. SC(12); IN BANDED FORM REQUIRED
!                   BY LINPACK SOLUTION ROUTINES
!       CMU      :  ABSCISSAE FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       CWT      :  WEIGHTS FOR GAUSS QUADRATURE OVER ANGLE COSINE
!       EXPBEA   :  TRANSMISSION OF INCIDENT BEAM, EXP(-TAUCPR/UMU0)
!       LYRCUT   :  LOGICAL FLAG FOR TRUNCATION OF COMPUT. LAYER
!       MAZIM    :  ORDER OF AZIMUTHAL COMPONENT
!       NCOL     :  COUNTS OF COLUMNS IN -CBAND-
!       NN       :  ORDER OF DOUBLE-GAUSS QUADRATURE (NSTR/2)
!       NCUT     :  TOTAL NUMBER OF COMPUTATIONAL LAYERS CONSIDERED
!       TPLANK   :  TOP BOUNDARY THERMAL EMISSION
!       TAUCPR   :  CUMULATIVE OPTICAL DEPTH (DELTA-M-SCALED)
!       ZZ       :  BEAM SOURCE VECTORS IN EQ. SS(19)
!       ZPLK0    :  THERMAL SOURCE VECTORS -Z0-, BY SOLVING EQ. SS(16)
!       ZPLK1    :  THERMAL SOURCE VECTORS -Z1-, BY SOLVING EQ. SS(16)
!       (REMAINDER ARE 'DISORT' INPUT VARIABLES)

!   O U T P U T     V A R I A B L E S:

!       B        :  RIGHT-HAND SIDE VECTOR OF EQ. SC(5) GOING INTO
!                   *SGBSL*; RETURNS AS SOLUTION VECTOR OF EQ.
!                   SC(12), CONSTANTS OF INTEGRATION WITHOUT
!                   EXPONENTIAL TERM
!      LL        :  PERMANENT STORAGE FOR -B-, BUT RE-ORDERED

!   I N T E R N A L    V A R I A B L E S:

!       IPVT     :  INTEGER VECTOR OF PIVOT INDICES
!       IT       :  POINTER FOR POSITION IN  -B-
!       NCD      :  NUMBER OF DIAGONALS BELOW OR ABOVE MAIN DIAGONAL
!       RCOND    :  INDICATOR OF SINGULARITY FOR -CBAND-
!       Z        :  SCRATCH ARRAY REQUIRED BY *SGBCO*
!+---------------------------------------------------------------------+

      LOGICAL  LAMBER, LYRCUT
      INTEGER  IPVT(*)
      REAL     B(*), BDR(MI,0:*), BEM(*), CBAND(MI9M2,NNLYRI),          &
     &         CMU(*), CWT(*), EXPBEA(0:*), LL(MXCMU,*),                &
     &         TAUCPR(0:*), Z(*), ZZ(MXCMU,*), ZPLK0(MXCMU,*),          &
     &         ZPLK1(MXCMU,*)

      CALL  ZEROIT(B, NNLYRI)
!                             ** CONSTRUCT -B-,  STWJ(20A,C) FOR
!                             ** PARALLEL BEAM + BOTTOM REFLECTION +
!                             ** THERMAL EMISSION AT TOP AND/OR BOTTOM

      IF (MAZIM.GT.0 .AND. FBEAM.GT.0.0)  THEN

!                                         ** AZIMUTH-DEPENDENT CASE
!                                         ** (NEVER CALLED IF FBEAM=0)
         IF (LYRCUT .OR. LAMBER) THEN

!               ** NO AZIMUTHAL-DEPENDENT INTENSITY FOR LAMBERT SURFACE;
!               ** NO INTENSITY COMPONENT FOR TRUNCATED BOTTOM LAYER

            DO 10  IQ=1, NN
!                                                     ** TOP BOUNDARY
               B(IQ)=- ZZ(NN+1-IQ,1)
!                                                  ** BOTTOM BOUNDARY
               B(NCOL-NN+IQ)=- ZZ(IQ+NN,NCUT) * EXPBEA(NCUT)
10          CONTINUE

         ELSE

            DO 20  IQ=1, NN
               B(IQ)=- ZZ(NN+1-IQ,1)

               SUM=0.
               DO 15  JQ=1, NN
                  SUM=SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)              &
     &                        * ZZ(NN+1-JQ,NCUT) * EXPBEA(NCUT)
15             CONTINUE
               B(NCOL-NN+IQ)=SUM
               IF (FBEAM.GT.0.0)                                        &
     &              B(NCOL-NN+IQ)=SUM + (BDR(IQ,0) * UMU0*FBEAM/PI      &
     &                                 - ZZ(IQ+NN,NCUT)) * EXPBEA(NCUT)
20          CONTINUE
         END IF
!                             ** CONTINUITY CONDITION FOR LAYER
!                             ** INTERFACES OF EQ. STWJ(20B)
         IT=NN
         DO 40  LC=1, NCUT-1
            DO 30  IQ=1, NSTR
               IT   =IT + 1
               B(IT)=(ZZ(IQ,LC+1) - ZZ(IQ,LC)) * EXPBEA(LC)
30          CONTINUE
40       CONTINUE

      ELSE
!                                   ** AZIMUTH-INDEPENDENT CASE
         IF (FBEAM.EQ.0.0)  THEN

            DO 50 IQ=1, NN
!                                      ** TOP BOUNDARY

               B(IQ)=- ZPLK0(NN+1-IQ,1) + FISOT + TPLANK
50          CONTINUE

            IF (LYRCUT) THEN
!                               ** NO INTENSITY COMPONENT FOR TRUNCATED
!                               ** BOTTOM LAYER
               DO 60 IQ=1, NN
!                                      ** BOTTOM BOUNDARY

                  B(NCOL-NN+IQ)=- ZPLK0(IQ+NN,NCUT)                     &
     &                            - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)
60             CONTINUE

            ELSE

               DO 80 IQ=1, NN

                  SUM=0.
                  DO 70 JQ=1, NN
                     SUM=SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)           &
     &                          * (ZPLK0(NN+1-JQ,NCUT)                  &
     &                            + ZPLK1(NN+1-JQ,NCUT) * TAUCPR(NCUT))
70                CONTINUE
                  B(NCOL-NN+IQ)=2.*SUM + BEM(IQ) * BPLANK               &
     &                            - ZPLK0(IQ+NN,NCUT)                   &
     &                            - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)
80             CONTINUE
            END IF
!                             ** CONTINUITY CONDITION FOR LAYER
!                             ** INTERFACES, STWJ(20B)
            IT=NN
            DO 100  LC=1, NCUT-1
               DO 90  IQ=1, NSTR
                  IT   =IT + 1
                  B(IT)=ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC) +                 &
     &                  (ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC)) * TAUCPR(LC)
90             CONTINUE
100         CONTINUE

         ELSE

            DO 150 IQ=1, NN
               B(IQ)=- ZZ(NN+1-IQ,1) - ZPLK0(NN+1-IQ,1) +FISOT +TPLANK
150         CONTINUE

            IF (LYRCUT) THEN

               DO 160 IQ=1, NN
                  B(NCOL-NN+IQ)=- ZZ(IQ+NN,NCUT) * EXPBEA(NCUT)         &
     &                            - ZPLK0(IQ+NN,NCUT)                   &
     &                            - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)
160            CONTINUE

            ELSE

               DO 180 IQ=1, NN

                  SUM=0.
                  DO 170 JQ=1, NN
                     SUM=SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)           &
     &                          * (ZZ(NN+1-JQ,NCUT) * EXPBEA(NCUT)      &
     &                            + ZPLK0(NN+1-JQ,NCUT)                 &
     &                            + ZPLK1(NN+1-JQ,NCUT) * TAUCPR(NCUT))
170               CONTINUE
                  B(NCOL-NN+IQ)=2.*SUM + (BDR(IQ,0) * UMU0*FBEAM/PI     &
     &                                 - ZZ(IQ+NN,NCUT)) * EXPBEA(NCUT) &
     &                            + BEM(IQ) * BPLANK                    &
     &                            - ZPLK0(IQ+NN,NCUT)                   &
     &                            - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)
180            CONTINUE
            END IF

            IT=NN
            DO 200  LC=1, NCUT-1
               DO 190  IQ=1, NSTR
                  IT   =IT + 1
                  B(IT)=(ZZ(IQ,LC+1) - ZZ(IQ,LC)) * EXPBEA(LC)          &
     &                    + ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC) +             &
     &                    (ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC)) * TAUCPR(LC)
190            CONTINUE
200         CONTINUE

         END IF

      END IF
!                     ** FIND L-U (LOWER/UPPER TRIANGULAR) DECOMPOSITION
!                     ** OF BAND MATRIX -CBAND- AND TEST IF IT IS NEARLY
!                     ** SINGULAR (NOTE: -CBAND- IS DESTROYED)
!                     ** (-CBAND- IS IN LINPACK PACKED FORMAT)
      RCOND=0.0
      NCD=3*NN - 1
      CALL  SGBCO(CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z)
      IF (1.0+RCOND .EQ. 1.0)  CALL  ERRMSG                             &
     &   ('SOLVE0--SGBCO SAYS MATRIX NEAR SINGULAR',.FALSE.)

!                   ** SOLVE LINEAR SYSTEM WITH COEFF MATRIX -CBAND-
!                   ** AND R.H. SIDE(S) -B- AFTER -CBAND- HAS BEEN L-U
!                   ** DECOMPOSED.  SOLUTION IS RETURNED IN -B-.

      CALL  SGBSL(CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0)

!                   ** ZERO -CBAND- (IT MAY CONTAIN 'FOREIGN'
!                   ** ELEMENTS UPON RETURNING FROM LINPACK);
!                   ** NECESSARY TO PREVENT ERRORS

      CALL  ZEROIT(CBAND, MI9M2*NNLYRI)

      DO 220  LC=1, NCUT
         IPNT=LC*NSTR - NN
         DO 220  IQ=1, NN
            LL(NN+1-IQ,LC)=B(IPNT+1-IQ)
            LL(IQ+NN,  LC)=B(IQ+IPNT)
220   CONTINUE

      RETURN
      END
