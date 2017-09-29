#include <misc.h>
      subroutine calendr(nstep   ,dtime   ,mdbase  ,msbase  ,mbdate  ,
     &                   mbsec   ,mdcur   ,mscur   ,mcdate  ,mcsec   ,
     &                   calday  )
!-----------------------------------------------------------------------
!
! Compute current date and day information for history tape header.
! Compute current julian day (including fraction) for sun angle 
! computations and time interpolation of boundary dataset information.
! One year is defined as 365 days.
!
! Computational notes: 
!
! 86400 is the number of seconds in 1 day.
!
! Dividing an integer by 10**n has the effect of right-shifting the           
! decimal digits n positions (ex: 861231/100 = 008612).
!
! mod(integer,10**n) has the effect of extracting the rightmost     
! n decimal digits of the integer (ex: mod(861231,10000) = 1231).
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Boville, August 1992
! Reviewed:          J. Kiehl, April 1996
! Reviewed:          B. Boville, April 1996
! Modified the type of the arguments to the mod function  -CAA-6/22/09
!   (There may be an issue here with the typecasting) 
!
!-----------------------------------------------------------------------
!
! $Id: calendr.F,v 1.2 1997/09/16 17:33:38 mvertens Exp $
! $Author: mvertens $
!
      use shr_kind_mod, only: r8 => shr_kind_r8  
      implicit none
!
! Input arguments
!
      integer nstep        ! current time step (0, 1, ...)
      real(r8) dtime       ! length of time step (seconds)
      integer mdbase       ! base day of run (e.g., 0)
      integer msbase       ! base seconds of base day (e.g., 0)
      real(r8) msbase_real ! base seconds of base day (e.g., 0) DRF
      integer mbdate       ! base date (yyyymmdd format) of run (e.g., 000901)
      integer mbsec        ! base seconds of base date (e.g., 0)
      real(r8) mbsec_real  ! base seconds of base date (e.g., 0) DRF
!
! Output arguments
!
      integer mdcur        ! current day (0, 1, ...)
      real(r8) mdcur_real  ! current day (0, 1, ...) DRF
      integer mscur        ! current seconds of current day (0, ..., 86400)
      real(r8) mscur_real  ! current seconds of current day (0, ..., 86400) DRF
      integer mcdate       ! current date (yyyymmdd format) (e.g., 021105)
      real(r8) mcdate_real ! current date (yyyymmdd format) (e.g., 021105) DRF
      integer mcsec        ! current seconds of current date (0, ..., 86400)
      real(r8) mcsec_real  ! current seconds of current date (0, ..., 86400) DRF
      real(r8) calday      ! current julian day including fraction
!
!---------------------------Local workspace-----------------------------
!
      integer nsecs        ! run time in seconds  
      real(r8) nsecs_real  ! run time in seconds DRF
      integer ndays        ! ndays + nyears = number of day changes since start
      integer nyears       ! ndays + nyears = number of day changes since start
      integer mcyear       ! current year of current date (0-99)
      integer mcmnth       ! current month of current date (1-12)
      integer mcday        ! current day of current date (1-31)
      integer jday         ! Julian day (1-365)
      integer mbmnth       ! base month (1-12) validity check
      integer mbday        ! base day of base date (1-31) validity check
      integer ndm(12)      ! number of days per month
      integer jdcon(12)    ! convert month index to julian day
!
      save ndm,jdcon
      data ndm/31,28,31,30,31,30,31,31,30,31,30,31/
      data jdcon/0,31,59,90,120,151,181,212,243,273,304,334/
!
!-----------------------------------------------------------------------
!
! Check validity of input data
!
      mbmnth = mod(mbdate,10000)/100
      mbday = mod(mbdate,100)
      if (mbmnth.lt.1 .or. mbmnth.gt.12) then
        write(6,*)' CALENDR: Invalid base month input:',mbmnth
        call endrun
      end if
      if (mbday.lt.1 .or. mbday.gt.ndm(mbmnth)) then
        write(6,*)' CALENDR: Invalid base day of base date input:',mbday
        call endrun
      end if
      if (msbase.lt.0 .or. msbase.ge.86400) then
        write(6,*)' CALENDR: Invalid base seconds(msbase):',msbase
        call endrun
      end if
      if (mbsec.lt.0 .or. mbsec.ge.86400) then
        write(6,*)' CALENDR: Invalid base seconds(mbsec):',mbsec
        call endrun
      end if
!
      !write(*,*) "dtime = ",nint(dtime)
      nsecs = nstep*nint(dtime)
      nsecs_real = dble(nstep)*dble(dtime)
      !write(*,*) "nsecs = ",nsecs
      !write(*,*) "nsecs = ",nsecs_real

!
! First: current day, seconds
!
      mdcur_real = dble(mdbase) + (nsecs_real+dble(msbase))/86400.
      mdcur = nint(mdcur_real)

      !write(*,*) "mdcur = ",mdcur
      !write(*,*) "mdcur = ",mdcur_real
      
      !DRF mdcur = mdbase + (nsecs+msbase)/86400
      mscur_real = mod((nsecs_real+dble(msbase)),86400.)
      !mscur = mod((nsecs+msbase),86400_r8)
      mscur = nint(mscur_real)

      !write(*,*) "mscur = ",mscur
      !write(*,*) "mscur = ",mscur_real
      
!
! Next: current date, seconds.
! Uncommenting the next line will have the effect of modifying  
! nsecs to include base day and base seconds.  This is the way             
! the current date was computed in CCM1.
!
!      nsecs = nsecs + mdbase*86400 + msbase

      !DRF ndays = (nsecs+mbsec)/86400
      ndays = nint((nsecs_real+dble(mbsec))/86400.)-15
      nyears = ndays/365
      ndays = mod(ndays,365)
      

      !write(*,*) "nyears = ",nyears
      !write(*,*) "ndays_two = ",ndays

!
! Current seconds of current date
!
      mcsec = mod(nsecs+mbsec,86400_r8)
!
! Initialize current year, month, day.
!
      mcyear = mbdate/10000 + nyears
      !write(*,*) "mbdate = ",mbdate
      !write(*,*) "mcyear = ",mcyear
      mcmnth = mod(mbdate,10000)/100
      mcday = mod(mbdate,100) + ndays
!
! Now loop through months, converting yyyy, mm, and ddd to yyyymmdd.
! ex: 791235 becomes 800104. 190001370 becomes 19010105.
!
   10 if (mcday.gt.ndm(mcmnth)) then
        mcday = mcday - ndm(mcmnth)
        mcmnth = mcmnth + 1
        if (mcmnth.eq.13) then             ! add a year
          mcyear = mcyear + 1
          mcmnth = 1
        end if
        go to 10
      end if
      mcdate = mcyear*10000 + mcmnth*100 + mcday
!
! Convert current month, day, seconds to Julian day + fraction
!
      jday = jdcon(mcmnth) + mcday
      calday = float(jday) + float(mcsec)/86400.

      return
      end
 
