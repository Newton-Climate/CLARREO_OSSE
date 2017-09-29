      SUBROUTINE CKD(V,SH2OT0,SH2OT1,FH2O)

!     FROM LBLRTM MODULE "CONTNM.F" REVISION 5.12 (3/24/00)
!     USES CKD_2.4, I.E. VERSION 2.4 OF THE CLOUGH KNEIZYS DATA

!     INPUT ARGUMENTS:
!       V        SPECTRAL FREQUENCY [CM-1]
      REAL V

!     OUTPUT ARGUMENTS:
!       SH2OT0   SELF-BROADENED H2O CONTINUA DATA AT 298K.
!       SH2OT1   SELF-BROADENED H2O CONTINUA DATA AT 260K.
!       FH2O     FOREIGN-BROADENED H2O CONTINUA DATA.
      REAL SH2OT0,SH2OT1,FH2O

!     PARAMETERS:
      REAL   HWSQS1,BETAS1,FACS1,SNUM1,                                 &
     &  V0S2,HWSQS2,       FACS2,SNUM2,V0S3,HWSQS3,BETAS3,FACS3,SNUM3,  &
     &  V0F0,HWSQF0,BETAF0,FACF0,FNUM0,V0F1,HWSQF1,BETAF1,FACF1,FNUM1,  &
     &  V0F2,HWSQF2,BETAF2,FACF2,FNUM2,V0F3,HWSQF3,BETAF3,FACF3,FNUM3
      PARAMETER(   HWSQS1=100.*100.,BETAS1=0.0001,FACS1= .688,          &
     &  V0S2=1050.,HWSQS2=200.*200.,              FACS2=-.2333,         &
     &  V0S3=1310.,HWSQS3=120.*120.,BETAS3=5.E-06,FACS3=-.15,           &
     &  V0F0= 350.,HWSQF0=200.*200.,BETAF0=5.E-09,FACF0=-.7,            &
     &  V0F1= 630.,HWSQF1= 65.* 65.,BETAF1=2.E-08,FACF1= .75,           &
     &  V0F2=1130.,HWSQF2=330.*330.,BETAF2=8.E-11,FACF2=-.97,           &
     &  V0F3=1975.,HWSQF3=250.*250.,BETAF3=5.E-06,FACF3=-.65,           &
     &  SNUM1=FACS1*HWSQS1,SNUM2=FACS2*HWSQS2,SNUM3=FACS3*HWSQS3,       &
     &  FNUM0=FACF0*HWSQF0,FNUM1=FACF1*HWSQF1,FNUM2=FACF2*HWSQF2,       &
     &  FNUM3=FACF3*HWSQF3)

!     LOCAL VARIABLES:
!       VFAC     FREQUENCY COEFFICIENT TIMES UNIT CONVERSION [CM-1].
      REAL VFAC,V0D2S1,V0D2S3,SFAC,V0D2F0,V0D2F1,V0D2F2,V0D2F3

!     DATA:

!     SELF-CONTINUUM MODIFICATION FACTORS FROM 700-1200 CM-1
      REAL XFAC(0:47)
      SAVE XFAC
      DATA (XFAC(I),I=0,47)/                                            &
     &    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,              &
     &    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,              &
     &    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,              &
     &    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,              &
     &    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,              &
     &    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,              &
     &    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,              &
     &    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000/

!     NO H2O CONTINUA DATA BELOW 0.5 MICRONS
      IF(V.GT.20000.)THEN
          SH2OT0=0.
          SH2OT1=0.
          FH2O=0.
      ELSE

!         RADIATION FIELD FREQUENCY COEFFICIENT TIMES UNIT CONVERSION:
          VFAC=1.E-20*V

!         SELF-BROADENED:
          V0D2S1=V**2
          V0D2S3=(V-V0S3)**2

!         INCORPORATE SPECTRAL FREQUENCY FROM RADIATION FIELD TERM:
          SFAC=VFAC*(1+SNUM1/(HWSQS1+     V0D2S1*(1+BETAS1*V0D2S1)))    &
     &             *(1+SNUM2/(HWSQS2+(V-V0S2)**2                  ))    &
     &             *(1+SNUM3/(HWSQS3+     V0D2S3*(1+BETAS3*V0D2S3)))
          IF(V.GT.700. .AND.  V.LT.1170.)SFAC=XFAC(INT(V/10-69.99))*SFAC

!         TEMPERATURE DEPENDENT TERMS:
          CALL SLF296(V,SH2OT0)
          CALL SLF260(V,SH2OT1)
          SH2OT0=SFAC*SH2OT0
          SH2OT1=SFAC*SH2OT1

!         FOREIGN-BROADENED:
          CALL FRN296(V,FH2O)
          V0D2F0=(V-V0F0)**2
          V0D2F1=(V-V0F1)**2
          V0D2F2=(V-V0F2)**2
          V0D2F3=(V-V0F3)**2
          FH2O=VFAC*FH2O*(1+FNUM0/(HWSQF0+V0D2F0*(1+BETAF0*V0D2F0**2))) &
     &                  *(1+FNUM1/(HWSQF1+V0D2F1*(1+BETAF1*V0D2F1**2))) &
     &                  *(1+FNUM2/(HWSQF2+V0D2F2*(1+BETAF2*V0D2F2**2))) &
     &                  *(1+FNUM3/(HWSQF3+V0D2F3*(1+BETAF3*V0D2F3   )))
      ENDIF
      RETURN
      END
