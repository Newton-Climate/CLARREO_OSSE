#include <misc.h>
#include <params.h>
      SUBROUTINE SUN_CALC( YEAR, DAY, HOUR, LAT, LONG, &
                       AZ, ZEN, SOLDIA, SOLDST )

!     Calculates the local solar azimuth and elevation angles, and
!     the distance to and angle subtended by the Sun, at a specific 
!     location and time using approximate formulas in The Astronomical 
!     Almanac.  Accuracy of angles is 0.01 deg or better (the angular 
!     width of the Sun is about 0.5 deg, so 0.01 deg is more than
!     sufficient for most applications).

!     Unlike many GCM (and other) sun angle routines, this
!     one gives slightly different sun angles depending on
!     the year.  The difference is usually down in the 4th
!     significant digit but can slowly creep up to the 3rd
!     significant digit after several decades to a century.

!     A refraction correction appropriate for the "US Standard
!     Atmosphere" is added, so that the returned sun position is
!     the APPARENT one.  The correction is below 0.1 deg for solar
!     elevations above 9 deg.  To remove it, comment out the code
!     block where variable REFRAC is set, and replace it with
!     REFRAC = 0.0.

!   References:

!     Michalsky, J., 1988: The Astronomical Almanac's algorithm for
!        approximate solar position (1950-2050), Solar Energy 40,
!        227-235 (but the version of this program in the Appendix
!        contains errors and should not be used)

!     The Astronomical Almanac, U.S. Gov't Printing Office, Washington,
!        D.C. (published every year): the formulas used from the 1995
!        version are as follows:
!        p. A12: approximation to sunrise/set times
!        p. B61: solar elevation ("altitude") and azimuth
!        p. B62: refraction correction
!        p. C24: mean longitude, mean anomaly, ecliptic longitude,
!                obliquity of ecliptic, right ascension, declination,
!                Earth-Sun distance, angular diameter of Sun
!        p. L2:  Greenwich mean sidereal time (ignoring T^2, T^3 terms)


!c   Authors:  Dr. Joe Michalsky (joe@asrc.albany.edu)
!c             Dr. Lee Harrison (lee@asrc.albany.edu)
!c             Atmospheric Sciences Research Center
!c             State University of New York
!c             Albany, New York

!c   Modified by:  Dr. Warren Wiscombe (wiscombe@climate.gsfc.nasa.gov)
!c                 NASA Goddard Space Flight Center
!c                 Code 913
!c                 Greenbelt, MD 20771


!   WARNING:  Do not use this routine outside the year range
!             1950-2050.  The approximations become increasingly
!             worse, and the calculation of Julian date becomes
!             more involved.

!   Input:

!      YEAR     year (INTEGER; range 1950 to 2050)

!      DAY      day of year at LAT-LONG location (INTEGER; range 1-366)

!      HOUR     hour of DAY [GMT or UT] (REAL; range -13.0 to 36.0)
!               = (local hour) + (time zone number)
!                 + (Daylight Savings Time correction; -1 or 0)
!               where (local hour) range is 0 to 24,
!               (time zone number) range is -12 to +12, and
!               (Daylight Time correction) is -1 if on Daylight Time
!               (summer half of year), 0 otherwise;  
!               Example: 8:30 am Eastern Daylight Time would be
!
!                           HOUR = 8.5 + 5 - 1 = 12.5

!      LAT      latitude [degrees]
!               (REAL; range -90.0 to 90.0; north is positive)

!      LONG     longitude [degrees]
!               (REAL; range -180.0 to 180.0; east is positive)


!   Output:

!      AZ       solar azimuth angle (measured east from north,
!               0 to 360 degs)

!      ZEN      solar zenith angle [-90 to 90 degs]; 
!               solar elevation angle = 90 - ZEN

!      SOLDIA   solar diameter [degs]

!      SOLDST   distance to sun [Astronomical Units, AU]
!               (1 AU = mean Earth-sun distance = 1.49597871E+11 m
!                in IAU 1976 System of Astronomical Constants)


!   Local Variables:

!     DEC       Declination (radians)

!     ECLONG    Ecliptic longitude (radians)

!     GMST      Greenwich mean sidereal time (hours)

!     HA        Hour angle (radians, -pi to pi)

!     JD        Modified Julian date (number of days, including 
!               fractions thereof, from Julian year J2000.0);
!               actual Julian date is JD + 2451545.0

!     LMST      Local mean sidereal time (radians)

!     MNANOM    Mean anomaly (radians, normalized to 0 to 2*pi)

!     MNLONG    Mean longitude of Sun, corrected for aberration 
!               (deg; normalized to 0-360)

!     OBLQEC    Obliquity of the ecliptic (radians)

!     RA        Right ascension  (radians)

!     REFRAC    Refraction correction for US Standard Atmosphere (degs)

! --------------------------------------------------------------------
!   Uses double precision for safety and because Julian dates can
!   have a large number of digits in their full form (but in practice
!   this version seems to work fine in single precision if you only
!   need about 3 significant digits and aren't doing precise climate
!   change or solar tracking work).
! --------------------------------------------------------------------

!   Why does this routine require time input as Greenwich Mean Time 
!   (GMT; also called Universal Time, UT) rather than "local time"?
!   Because "local time" (e.g. Eastern Standard Time) can be off by
!   up to half an hour from the actual local time (called "local mean
!   solar time").  For society's convenience, "local time" is held 
!   constant across each of 24 time zones (each technically 15 longitude
!   degrees wide although some are distorted, again for societal 
!   convenience).  Local mean solar time varies continuously around a 
!   longitude belt;  it is not a step function with 24 steps.  
!   Thus it is far simpler to calculate local mean solar time from GMT,
!   by adding 4 min for each degree of longitude the location is
!   east of the Greenwich meridian or subtracting 4 min for each degree
!   west of it.  

! --------------------------------------------------------------------

!   TIME
!   
!   The measurement of time has become a complicated topic.  A few
!   basic facts are:
!   
!   (1) The Gregorian calendar was introduced in 1582 to replace 
!   Julian calendar; in it, every year divisible by four is a leap 
!   year just as in the Julian calendar except for centurial years
!   which must be exactly divisible by 400 to be leap years.  Thus 
!   2000 is a leap year, but not 1900 or 2100.

!   (2) The Julian day begins at Greenwich noon whereas the calendar 
!   day begins at the preceding midnight;  and Julian years begin on
!   "Jan 0" which is really Greenwich noon on Dec 31.  True Julian 
!   dates are a continous count of day numbers beginning with JD 0 on 
!   1 Jan 4713 B.C.  The term "Julian date" is widely misused and few
!   people understand it; it is best avoided unless you want to study
!   the Astronomical Almanac and learn to use it correctly.

!   (3) Universal Time (UT), the basis of civil timekeeping, is
!   defined by a formula relating UT to GMST (Greenwich mean sidereal
!   time).  UTC (Coordinated Universal Time) is the time scale 
!   distributed by most broadcast time services.  UT, UTC, and other
!   related time measures are within a few sec of each other and are
!   frequently used interchangeably.

!   (4) Beginning in 1984, the "standard epoch" of the astronomical
!   coordinate system is Jan 1, 2000, 12 hr TDB (Julian date 
!   2,451,545.0, denoted J2000.0).  The fact that this routine uses
!   1949 as a point of reference is merely for numerical convenience.
! --------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8


!     .. Scalar Arguments ..

      INTEGER   YEAR, DAY
      real(r8)      AZ, ZEN, HOUR, LAT, LONG, SOLDIA, SOLDST,EL
!     ..
!     .. Local Scalars ..

      LOGICAL   PASS1
      INTEGER   DELTA, LEAP
      real(r8)  DEC, DEN, ECLONG, GMST, HA, JD, LMST, &
                       MNANOM, MNLONG, NUM, OBLQEC, PI, RA, &
                       RPD, REFRAC, TIME, TWOPI
!     ..
!     .. Intrinsic Functions ..

      !INTRINSIC AINT, ASIN, ATAN, COS, MOD, SIN, TAN
!     ..
!     .. Data statements ..

      !real(r8) PI, TWOPI, RPD
      !SAVE     PASS1, PI, TWOPI, RPD
      !DATA     PASS1 /.True./
!     ..

      IF( YEAR.LT.1950 .OR. YEAR.GT.2050 ) &
         STOP 'SUNAE--bad input variable YEAR'
      IF( DAY.LT.1 .OR. DAY.GT.366 ) &
         STOP 'SUNAE--bad input variable DAY'
      IF( HOUR.LT.-13.0 .OR. HOUR.GT.36.0 ) &
         STOP 'SUNAE--bad input variable HOUR'
      IF( LAT.LT.-90.0 .OR. LAT.GT.90.0 ) &
         STOP 'SUNAE--bad input variable LAT'
      IF( LONG.LT.-180.0 .OR. LONG.GT.180.0 ) &
         STOP 'SUNAE--bad input variable LONG'

      !IF(PASS1) THEN
         PI     = 2.*ASIN( 1.0 )
         TWOPI  = 2.*PI
         RPD    = PI / 180.
      !   PASS1 = .False.
      !ENDIF

!                    ** current Julian date (actually add 2,400,000 
!                    ** for true JD);  LEAP = leap days since 1949;
!                    ** 32916.5 is midnite 0 jan 1949 minus 2.4e6

      DELTA  = YEAR - 1949
      LEAP   = DELTA / 4
      JD     = 32916.5 + (DELTA*365 + LEAP + DAY) + HOUR / 24.

!                    ** last yr of century not leap yr unless divisible
!                    ** by 400 (not executed for the allowed YEAR range,
!                    ** but left in so our successors can adapt this for 
!                    ** the following 100 years)

      IF( MOD( YEAR, 100 ).EQ.0 .AND. &
         MOD( YEAR, 400 ).NE.0 ) JD = JD - 1.

!                     ** ecliptic coordinates
!                     ** 51545.0 + 2.4e6 = noon 1 jan 2000

      TIME  = JD - 51545.0

!                    ** force mean longitude between 0 and 360 degs

      MNLONG = 280.460 + 0.9856474*TIME
      MNLONG = MOD( MNLONG, 360. )
      IF( MNLONG.LT.0. ) MNLONG = MNLONG + 360.

!                    ** mean anomaly in radians between 0 and 2*pi

      MNANOM = 357.528 + 0.9856003*TIME
      MNANOM = MOD( MNANOM, 360. )
      IF( MNANOM.LT.0.) MNANOM = MNANOM + 360.

      MNANOM = MNANOM*RPD

!                    ** ecliptic longitude and obliquity 
!                    ** of ecliptic in radians

      ECLONG = MNLONG + 1.915*SIN( MNANOM ) + 0.020*SIN( 2.*MNANOM )
      ECLONG = MOD( ECLONG, 360. )
      IF( ECLONG.LT.0. ) ECLONG = ECLONG + 360.

      OBLQEC = 23.439 - 0.0000004*TIME
      ECLONG = ECLONG*RPD
      OBLQEC = OBLQEC*RPD

!                    ** right ascension

      NUM  = COS( OBLQEC )*SIN( ECLONG )
      DEN  = COS( ECLONG )
      RA   = ATAN( NUM / DEN )

!                    ** Force right ascension between 0 and 2*pi

      IF( DEN.LT.0.0 ) THEN
         RA  = RA + PI
      ELSE IF( NUM.LT.0.0 ) THEN
         RA  = RA + TWOPI
      END IF

!                    ** declination

      DEC  = ASIN( SIN( OBLQEC )*SIN( ECLONG ) )

!                    ** Greenwich mean sidereal time in hours

      GMST = 6.697375 + 0.0657098242*TIME + HOUR

!                    ** Hour not changed to sidereal time since 
!                    ** 'time' includes the fractional day

      GMST  = MOD( GMST, 24. )
      IF( GMST.LT.0. ) GMST   = GMST + 24.

!                    ** local mean sidereal time in radians

      LMST  = GMST + LONG / 15.
      LMST  = MOD( LMST, 24. )
      IF( LMST.LT.0. ) LMST   = LMST + 24.

      LMST   = LMST*15.*RPD

!                    ** hour angle in radians between -pi and pi

      HA  = LMST - RA

      IF( HA.LT.- PI ) HA  = HA + TWOPI
      IF( HA.GT.PI )   HA  = HA - TWOPI

!                    ** solar azimuth and elevation

      EL  = ASIN( SIN( DEC )*SIN( LAT*RPD ) + &
                 COS( DEC )*COS( LAT*RPD )*COS( HA ) )

      AZ  = ASIN( - COS( DEC )*SIN( HA ) / COS( EL ) )

!                    ** Put azimuth between 0 and 2*pi radians

      IF( SIN( DEC ) - SIN( EL )*SIN( LAT*RPD ).GE.0. ) THEN

         IF( SIN(AZ).LT.0.) AZ  = AZ + TWOPI

      ELSE

         AZ  = PI - AZ

      END IF

!                     ** Convert elevation and azimuth to degrees
      EL  = EL / RPD
      AZ  = AZ / RPD

!  ======== Refraction correction for U.S. Standard Atmos. ==========
!      (assumes elevation in degs) (3.51823=1013.25 mb/288 K)

      IF( EL.GE.19.225 ) THEN

         REFRAC = 0.00452*3.51823 / TAN( EL*RPD )

      ELSE IF( EL.GT.-0.766 .AND. EL.LT.19.225 ) THEN

         REFRAC = 3.51823 * ( 0.1594 + EL*(0.0196 + 0.00002*EL) ) / &
                 ( 1. + EL*(0.505 + 0.0845*EL) )

      ELSE IF( EL.LE.-0.766 ) THEN

         REFRAC = 0.0

      END IF

      EL  = EL + REFRAC
      ZEN = 90. - EL
! ===================================================================

!                   ** distance to sun in A.U. & diameter in degs

      SOLDST = 1.00014 - 0.01671*COS(MNANOM) - 0.00014*COS( 2.*MNANOM )
      SOLDIA = 0.5332 / SOLDST

      IF( EL.LT.-90.0 .OR. EL.GT.90.0 ) &
         STOP 'SUNAE--output argument EL out of range'
      IF( AZ.LT.0.0 .OR. AZ.GT.360.0 ) &
         STOP 'SUNAE--output argument AZ out of range'

      RETURN

      END


