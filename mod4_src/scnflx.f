      SUBROUTINE SCNFLX(NUM_LRES,SOLAR_FLUX,DIFFUSE_FLUX,               &
     &                NUM_LEVS,FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)

!     SCNFLX WRITES OUT SPECTRAL FLUX, DEGRADED WITH
!     A SCANNING FUNCTION, TO UNIT IFLUX.

!     LIST PARAMETERS:
      INCLUDE 'PARAMS.h'

!     COMMONS:
      INCLUDE 'BMHEAD.h'

!     /CARD4/
!       IV1      LOWEST FREQUENCY OUTPUT [CM-1].
!       IV2      HIGHEST FREQUENCY OUTPUT [CM-1].
!       IDV      PRINTOUT FREQUENCY STEP SIZE [CM-1].
!       IFWHM    TRIANGULAR SLIT FULL-WIDTH-HALF-MAXIMUM [CM-1].
!       IVX      CURRENT COMPUTATION FREQUENCY [CM-1].
!       IVOFF    OFFSET BETWEEN COMPUTATION AND OUTPUT FREQUENCIES,
!                REQUIRED FOR SLIT FUNCTION [CM-1].
!       IWRITE   COMPUTATION FREQUENCY OF NEXT WRITE [CM-1].
!       NSPCDT   NUMBER OF OUTPUT SPECTRAL DATA POINTS.
!       NWGT     NUMBER OF SPECTRAL BINS CONTRIBUTING TO SLIT FUNCTION.
!       WGT      NORMALIZED WEIGHTS USED TO DEFINE THE SLIT FUNCTION.
      INTEGER IV1,IV2,IDV,IFWHM,IVX,IVOFF,IWRITE,NSPCDT,NWGT
      REAL WGT
      COMMON/CARD4/IV1,IV2,IDV,IFWHM,IVX,IVOFF,IWRITE,NSPCDT,           &
     &  NWGT,WGT(NBINS)

      !DRF
      INTEGER NUM_LRES
      REAL SOLAR_FLUX(NUM_LRES),DIFFUSE_FLUX(NUM_LRES)
      INTEGER NUM_LEVS
      REAL FLX_UPDIFFUSE(NUM_LEVS,NUM_LRES)
      REAL FLX_DNDIFFUSE(NUM_LEVS,NUM_LRES) 
      REAL FLX_DNDIRECT(NUM_LEVS,NUM_LRES)
      !DRF

!     COMMON/SCAN/
!       V1       LOWER BOUND ON SPECTRAL RANGE [CHUNIT DEFINES UNIT].
!       V2       UPPER BOUND ON SPECTRAL RANGE [CHUNIT DEFINES UNIT].
!       DV       SPECTRAL STEP SIZE FOR OUTPUT [CHUNIT DEFINES UNIT].
!       FWHM     FULL-WIDTH-AT-HALF-MAXIMUM [CHUNIT DEFINES UNIT].
!       FWHMSQ   TRIANGULAR SLIT NORMALIZATION FACTOR
!                (EQUALS FWHM SQUARED) [CHUNIT DEFINES UNIT].
!       VOUT     CURRENT SPECTRAL OUTPUT [CHUNIT DEFINES UNIT].
      REAL V1,V2,DV,FWHM,FWHMSQ,VOUT
      COMMON/SCAN/V1,V2,DV,FWHM,FWHMSQ,VOUT

!     COMMON/CSCAN/
!       CHUNIT   UNIT FLAG ('W'=WAVENUMBERS;'M'=MICRONS;'N'=NANOMETERS).
!       RELABS   SPECTRAL RESOLUTION FLAG('A'=ABSOLUTE;'R'=RELATIVE[%]).
!       LNFEED   LINE FEED FLAG FOR .FLX FILE ('T' FOR 80 CHARACTER
!                  LINES, 'F' FOR LONG LINES, ' ' FOR NO .FLX FILE).
      CHARACTER*1 CHUNIT,RELABS,LNFEED
      COMMON/CSCAN/CHUNIT,RELABS,LNFEED

!     FUNCTIONS:
!       WTRGT    RIGHT SIDE WAVELENGTH TRIANGULAR SLIT INTEGRATION.
!       WTLFT    LEFT SIDE WAVELENGTH TRIANGULAR SLIT INTEGRATION.
      REAL WTRGT,WTLFT

!     LOCAL VARIABLES:
!       JFLUX    SPECTRAL INDEX FOR OUTPUT FLUX VALUES.
!       JOFF     SPECTRAL INDEX OFFSET (EQUALS NUMBER OF SLIT
!                FUNCTIONS TERMINATING IN CURRENT SPECTRAL BIN).
!       LAST     LOGICAL FLAG, TRUE IF CURRENT SLIT FUNCTION
!                TERMINATES IN CURRENT SPECTRAL BIN.
!       CONVRT   SPECTRAL CONVERSION FACTOR [WAVELENGTHS / CM].
!       VOUTCN   SLIT FUNCTION CENTRAL WAVELENGTH [MICRONS OR NM].
!       VOUTMN   SLIT FUNCTION MINIMUM WAVELENGTH [MICRONS OR NM].
!       VOUTMX   SLIT FUNCTION MAXIMUM WAVELENGTH [MICRONS OR NM].
!       WNCEN    SLIT FUNCTION CENTRAL FREQUENCY [CM-1].
!       WNMIN    SLIT FUNCTION MINIMUM FREQUENCY [CM-1].
!       WNMAX    SLIT FUNCTION MAXIMUM FREQUENCY [CM-1].
!       WNFWHM   FULL-WIDTH-AT-HALF-MAXIMUM [CM-1].
!       BINMIN   SPECTRAL BIN MINIMUM FREQUENCY [CM-1].
!       BINMAX   SPECTRAL BIN MAXIMUM FREQUENCY [CM-1].
!       CONFLX   FLUX SPECTRAL CONVERSION FACTOR [CM-1 / SPECTRAL UNIT].
!       DIFF     FREQUENCY DIFFERENCE OF INTEGRAL LIMITS [CM-1].
!       COEF     COEFFICIENT USED TO EVALUATE INTEGRAL [WAVELENGTH CM].
!       AVG      AVERAGE OF INTEGRAL LIMITS [CM-1].
!       WT       WEIGHT OF SPECTRAL BIN IN EVALUATING SLIT FUNCTION.
      INTEGER JFLUX,JOFF
      LOGICAL LAST
      REAL CONVRT,VOUTCN,VOUTMN,VOUTMX,WNCEN,WNMIN,WNMAX,WNFWHM,        &
     &  BINMIN,BINMAX,CONFLX,DIFF,COEF,AVG,WT

!     THERE ARE 10 POSSIBLE OVERLAP SCENARIOS BETWEEN THE TRIANGULAR
!     SLIT AND THE SPECTRAL BIN (THE ABSCISSA IS WAVELENGTH):

!                                  / \
!                                 / | \
!                                /  |  \
!                               /   |   \
!                              /    |    \
!                             /     |     \
!                            /      |      \
!                           /       |       \
!                          /        |        \
!                         /         |         \

!     CASE A:                                      |__|

!     CASE B:                            |____________|

!     CASE C:               |_________________________|

!     CASE D:       |_________________________________|

!     CASE E:                            |__|

!     CASE F:               |_______________|

!     CASE C:       |_______________________|

!     CASE E:               |__|

!     CASE I:       |__________|

!     CASE J:       |__|

!     SKIP WRITING FLUX FILE IF LNFEED IS BLANK:
      !WRITE(*,*) 'into scnflx'
      IF(LNFEED.EQ.' ')RETURN

!     BRANCH BASED ON UNITS:
      IF(CHUNIT.EQ.'W')THEN

!         WAVENUMBERS.

!         INITIALIZE OUTPUT:
          JFLUX=0
          JOFF=0
          LAST=.TRUE.
          WNCEN=VOUT
          BINMIN=IVX-.5*IBNDWD
          BINMAX=IVX+.5*IBNDWD
          CONFLX=1.
   10     CONTINUE
          IF(WNCEN.GT.V2)RETURN
          IF(RELABS.EQ.'R')THEN

!             FWHM IS A RELATIVE [%] RESOLUTION.
              WNFWHM=.01*FWHM*WNCEN
              FWHMSQ=WNFWHM**2
              WNMIN=(1.-.01*FWHM)*WNCEN
              WNMAX=(1.+.01*FWHM)*WNCEN
          ELSE

!             FWHM IS A ABSOLUTE.
              WNFWHM=FWHM
              WNMIN=WNCEN-FWHM
              WNMAX=WNCEN+FWHM
          ENDIF

!         BRANCH BASED ON SPECTRAL OVERLAP:
          IF(BINMAX.LE.WNMIN)THEN

!             CASE A:  BIN DOES NOT CONTRIBUTE TO SLIT.
              RETURN
          ELSEIF(BINMIN.LE.WNMIN)THEN

!             INTEGRATION BEGINS AT WNMIN:
              IF(BINMAX.LE.WNCEN)THEN

!                  CASE B:  WNMIN TO BINMAX INTEGRAL.
                   DIFF=BINMAX-WNMIN
                   AVG=.5*DIFF
                   WT=DIFF*AVG/FWHMSQ
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
              ELSEIF(BINMAX.LT.WNMAX)THEN

!                  CASE C:  WNCEN TO BINMAX INTEGRAL + 1/2.
                   DIFF=BINMAX-WNCEN
                   AVG=WNFWHM-.5*DIFF
                   WT=.5+DIFF*AVG/FWHMSQ
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
               ELSE

!                  CASE D:  WNMIN TO WNMAX INTEGRAL.
                   WT=1.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
                   JOFF=JFLUX
                   VOUT=VOUT+DV
               ENDIF
           ELSEIF(BINMIN.LT.WNCEN)THEN

!              INTEGRATION BEGINS AT BINMIN:
               IF(BINMAX.LE.WNCEN)THEN

!                  CASE E:  BINMIN TO BINMAX INTEGRAL - RIGHT SIDE.
                   DIFF=FLOAT(IBNDWD)
                   AVG=FLOAT(IVX)-WNMIN
                   WT=DIFF*AVG/FWHMSQ
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
               ELSEIF(BINMAX.LT.WNMAX)THEN

!                  CASE F:  BINMIN TO BINMAX INTEGRAL - BOTH SIDES.
                   DIFF=WNCEN-BINMIN
                   AVG=WNFWHM-.5*DIFF
                   WT=DIFF*AVG/FWHMSQ
                   DIFF=BINMAX-WNCEN
                   AVG=WNFWHM-.5*DIFF
                   WT=WT+DIFF*AVG/FWHMSQ
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
               ELSE

!                  CASE G:  BINMIN TO WNCEN INTEGRAL + 1/2.
                   DIFF=WNCEN-BINMIN
                   AVG=WNFWHM-.5*DIFF
                   WT=.5+DIFF*AVG/FWHMSQ
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
                   JOFF=JFLUX
                   VOUT=VOUT+DV
               ENDIF
           ELSEIF(BINMIN.LT.WNMAX)THEN

!              INTEGRATION BEGINS AT BINMIN:
               IF(BINMAX.LT.WNMAX)THEN

!                  CASE H:  BINMIN TO BINMAX INTEGRAL - LEFT SIDE.
                   DIFF=FLOAT(IBNDWD)
                   AVG=WNMAX-FLOAT(IVX)
                   WT=DIFF*AVG/FWHMSQ
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
               ELSE

!                  CASE I:  BINMIN TO WNMAX INTEGRAL
                   DIFF=WNMAX-BINMIN
                   AVG=.5*DIFF
                   WT=DIFF*AVG/FWHMSQ
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
                   JOFF=JFLUX
                   VOUT=VOUT+DV
               ENDIF
           ELSE

!             CASE J:  BIN DOES NOT CONTRIBUTE TO SLIT.
              RETURN
           ENDIF
           WNCEN=WNCEN+DV
           GOTO10
      ELSE

!         MICRONS OR NANOMETERS.

!         INITIALIZE OUTPUT:
          JFLUX=0
          JOFF=0
          LAST=.TRUE.
          VOUTCN=VOUT
          BINMIN=IVX-.5*IBNDWD
          BINMAX=IVX+.5*IBNDWD
   20     CONTINUE
          IF(VOUTCN.LT.V1)RETURN
          IF(CHUNIT.EQ.'M')THEN
              CONVRT=1.E4
          ELSE
              CONVRT=1.E7
          ENDIF
          IF(RELABS.EQ.'R')THEN

!             FWHM IS A RELATIVE [%] RESOLUTION.
              FWHMSQ=(.01*FWHM*VOUTCN)**2
              VOUTMN=(1.-.01*FWHM)*VOUTCN
              VOUTMX=(1.+.01*FWHM)*VOUTCN
          ELSE

!             FWHM IS A ABSOLUTE.
              VOUTMN=VOUTCN-FWHM
              VOUTMX=VOUTCN+FWHM
          ENDIF
          WNCEN=CONVRT/VOUTCN
          WNMIN=CONVRT/VOUTMX
          WNMAX=CONVRT/VOUTMN
          CONFLX=WNCEN/VOUTCN

!         BRANCH BASED ON SPECTRAL OVERLAP:
          IF(BINMAX.LE.WNMIN)THEN

!             CASE A:  BIN DOES NOT CONTRIBUTE TO SLIT.
              RETURN
          ELSEIF(BINMIN.LE.WNMIN)THEN

!             INTEGRATION BEGINS AT WNMIN:
              IF(BINMAX.LE.WNCEN)THEN

!                  CASE B:  WNMIN TO BINMAX INTEGRAL.
                   DIFF=BINMAX-WNMIN
                   COEF=VOUTMX/BINMAX
                   AVG=.5*(WNMIN+BINMAX)
                   WT=WTRGT(DIFF,COEF,AVG,VOUTMX,FWHMSQ)
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
              ELSEIF(BINMAX.LT.WNMAX)THEN

!                  CASE C:  WNCEN TO BINMAX INTEGRAL + 1/2.
                   DIFF=BINMAX-WNCEN
                   COEF=VOUTCN/BINMAX
                   AVG=.5*(WNCEN+BINMAX)
                   WT=.5+WTLFT(DIFF,COEF,AVG,VOUTMN,FWHMSQ)
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
               ELSE

!                  CASE D:  WNMIN TO WNMAX INTEGRAL.
                   WT=1.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
                   JOFF=JFLUX
                   VOUT=VOUT-DV
               ENDIF
           ELSEIF(BINMIN.LT.WNCEN)THEN

!              INTEGRATION BEGINS AT BINMIN:
               IF(BINMAX.LE.WNCEN)THEN

!                  CASE E:  BINMIN TO BINMAX INTEGRAL - RIGHT SIDE.
                   DIFF=FLOAT(IBNDWD)
                   COEF=CONVRT/(BINMIN*BINMAX)
                   AVG=FLOAT(IVX)
                   WT=WTRGT(DIFF,COEF,AVG,VOUTMX,FWHMSQ)
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
               ELSEIF(BINMAX.LT.WNMAX)THEN

!                  CASE F:  BINMIN TO BINMAX INTEGRAL - BOTH SIDES.
                   DIFF=WNCEN-BINMIN
                   COEF=VOUTCN/BINMIN
                   AVG=.5*(BINMIN+WNCEN)
                   WT=WTRGT(DIFF,COEF,AVG,VOUTMX,FWHMSQ)
                   DIFF=BINMAX-WNCEN
                   COEF=VOUTCN/BINMAX
                   AVG=.5*(WNCEN+BINMAX)
                   WT=WT+WTLFT(DIFF,COEF,AVG,VOUTMN,FWHMSQ)
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
               ELSE

!                  CASE G:  BINMIN TO WNCEN INTEGRAL + 1/2.
                   DIFF=WNCEN-BINMIN
                   COEF=VOUTCN/BINMIN
                   AVG=.5*(BINMIN+WNCEN)
                   WT=.5+WTRGT(DIFF,COEF,AVG,VOUTMX,FWHMSQ)
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
                   JOFF=JFLUX
                   VOUT=VOUT-DV
               ENDIF
           ELSEIF(BINMIN.LT.WNMAX)THEN

!              INTEGRATION BEGINS AT BINMIN:
               IF(BINMAX.LT.WNMAX)THEN

!                  CASE H:  BINMIN TO BINMAX INTEGRAL - LEFT SIDE.
                   DIFF=FLOAT(IBNDWD)
                   COEF=CONVRT/(BINMIN*BINMAX)
                   AVG=FLOAT(IVX)
                   WT=WTLFT(DIFF,COEF,AVG,VOUTMN,FWHMSQ)
                   LAST=.FALSE.
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
               ELSE

!                  CASE I:  BINMIN TO WNMAX INTEGRAL
                   DIFF=WNMAX-BINMIN
                   COEF=VOUTMN/BINMIN
                   AVG=.5*(BINMIN+WNMAX)
                   WT=WTLFT(DIFF,COEF,AVG,VOUTMN,FWHMSQ)
                   CALL WTSUM(JFLUX,JOFF,CONFLX,WT,LAST,NUM_LRES,       &
     &                       SOLAR_FLUX,DIFFUSE_FLUX,NUM_LEVS,          &
     &                       FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
                   JOFF=JFLUX
                   VOUT=VOUT-DV
               ENDIF
           ELSE

!             CASE J:  BIN DOES NOT CONTRIBUTE TO SLIT.
              RETURN
           ENDIF
           VOUTCN=VOUTCN-DV
           GOTO20
      ENDIF
      END
