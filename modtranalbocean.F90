MODULE ZJIN_OCEAN_WGTLUT
  implicit none
  real*4 alball(24, 16*15*7*5), rflxall(24, 16*15 )
  logical :: lzjlut = .false.

  !=================================================================
  contains
  
  subroutine zjin_wgtlut_init
    integer i,j
    real*4 alball_local, rflxall_local

    open(83,file='ocnalbtab24bnd.bin',form='Unformatted',access='Direct',recl=829440)
    !!829440 = 24*4*16*15*7*5 + 24*4*16*15

    read ( 83,rec=1)alball,rflxall
    close(83)
    do i=1,24
       do j=1,8400
          CALL native_4byte_real( alball(i,j), alball_local)
          alball(i,j) = alball_local
       end do
    end do
    do i=1,24
       do j=1,240
          CALL native_4byte_real( rflxall(i,j), rflxall_local)
          rflxall(i,j) = rflxall_local
       end do
    end do
  end subroutine zjin_wgtlut_init

  !==============================================================
  subroutine zjin_wgtlut( tau_in ,u0 , wind_in ,chl_in , fuspecalb,fu_wvl )
     integer ,parameter :: nb=24, nt=16, ns=15, nw=7, nc=5

     real*4 ,intent(IN) :: tau_in  !cloud/aerosol optical depth
     real*4, intent(IN) :: u0      !cosine solar zenith angle
     real*4, intent(IN) :: wind_in !wind speed (m/s)
     real*4, intent(IN) :: chl_in  !chlorophyll concentration in mg/m^3
     real*4, intent(OUT):: fuspecalb(15) !fu-liou spectral albedo output
     real*4, intent(OUT):: fu_wvl(15) !fu-liou spectral albedo output

     integer ifu,ib,ib1,ib2,itt,iss,iww,icc,irr,irec
     integer it,is,ic,iw
     real*4 albedo,wtb,twtb,wtt,wtf,wls,wle
     real*4 tau,wind,chl
     integer i
     real*4,dimension(nb)::  alb,rflx, albb,wtflx

     real*4  taunode(16), u0node(15), windnode(7), chlnode(5)
     real*4  wlnode(25)
     real*4,dimension(2)::  rt, rs, rw, rc

     real*4  fuwl(16) ! FULIOU SW BAND Bounds
     data fuwl/0.1754,0.2247,0.2439,0.2857,0.2985,&
	   0.3225,0.3575,0.4375,0.4975,0.5950,&
	   0.6897,1.2987,1.9048,2.5000,3.5088,4.000/

     data taunode /0.0, 0.05, 0.1, 0.16, 0.24, 0.35, 0.5, 0.7, 0.99, &
                      1.3, 1.8, 2.5, 5.0, 9.0, 15.0, 25.0 /
     data u0node /0.05, 0.09, 0.15, 0.21, 0.27, 0.33, 0.39, 0.45,&
                     0.52, 0.60, 0.68, 0.76, 0.84, 0.92, 1.0 /
     data windnode /0., 3., 6., 9., 12., 15., 18. /
     data chlnode /0.0, 0.1, 0.5, 2.0, 12.0/
     data wlnode /0.25, 0.30, 0.33, 0.36,  0.40,  0.44, 0.48, 0.52, &
                 0.57, 0.64, 0.69, 0.752, 0.780, 0.867, 1.0, 1.096, &
                 1.19, 1.276, 1.534, 1.645, 2.128, 2.381, 2.907,  3.425, 4.0/

!    now find the albedo corresponding to the 4 parameters above:

      if( .not. lzjlut ) call zjin_wgtlut_init

!-------

!!if(tau.lt.0. .or. (u0.lt.0. .or. u0.gt.1.) .or. wind.lt.0. &
!!.or. chl.lt.0.)stop 'Err: input parameters wrong!'
	 tau  = tau_in
	 wind = wind_in
	 chl  = chl_in
	if(tau .gt. taunode(nt)) tau=taunode(nt)
	if(wind .gt. windnode(nw))wind=windnode(nw)
	if(chl .gt. 45.)chl=45.

	call locate(taunode,nt,tau,it)
	call locate(u0node,ns,u0,is)
	call locate(windnode,nw,wind,iw)
	call locate(chlnode,nc,chl,ic)

        rt(2) = (tau-taunode(it))/(taunode(it+1)-taunode(it)) ; rt(1) = 1.0-rt(2)
        rs(2) = (u0-u0node(is))/(u0node(is+1)-u0node(is)) ;  rs(1) = 1.0-rs(2)
        rw(2) = (wind-windnode(iw))/(windnode(iw+1)-windnode(iw)) ; rw(1) = 1.0-rw(2)
        rc(2) = (chl-chlnode(ic))/(chlnode(ic+1)-chlnode(ic)) ;rc(1) = 1.0-rc(2)

        do i=1,15
           fu_wvl(i) = fuwl(i)
        end do
! ** get alb(ib1:ib2) by 4 dimensional linear interpolaton **
	
        albb(1:24) = 0.0
        wtflx(1:24) = 0.0

	do  itt=it,it+1
	do  iss=is,is+1
	 do  iww=iw,iw+1
	 do  icc=ic,ic+1
	  irec = (itt-1)*15*7*5 + (iss-1)*7*5 + (iww-1)*5 + icc

	  alb(1:24) = alball(1:24,irec)

	  wtt = rt(itt-it+1)*rs(iss-is+1)*rw(iww-iw+1)*rc(icc-ic+1)

!	  do ib=ib1,ib2
	  do ib=1,24
	    albb(ib) = albb(ib) + wtt*alb(ib)
	  enddo
		

	enddo
	enddo

! *** get 24 band down flux weights ***
	 irr = (itt-1)*15 + iss

	  rflx(1:24) = rflxall(1:24,irr)

         wtf = rt(itt-it+1)*rs(iss-is+1)

	 do ib=1,24
	   wtflx(ib) = wtflx(ib) + wtf*rflx(ib)
	 enddo

	enddo
	enddo
!--------------------------------------------------

FUBANDS : do ifu = 1,15 

wls = fuwl(ifu)
wle = fuwl(ifu+1)

!!if(wls .gt. wle)stop 'Err: Start wavelength should be smaller.'

	if(wls .lt. wlnode(1))wls=wlnode(1)
	if(wle .gt. wlnode(25))wle=wlnode(25)
	call locate(wlnode,25,wls,ib1)
	call locate(wlnode,25,wle,ib2)


!** get albedo in the specified band by weighted sum **
	twtb = 0.0
	albedo = 0.0
	do ib=ib1,ib2
	  if(ib .eq. ib1)then
	    wtb = (wlnode(ib1+1)-wls)/(wlnode(ib1+1)-wlnode(ib1))
	  else if (ib .eq. ib2)then
	    wtb = (wle-wlnode(ib2))/(wlnode(ib2+1)-wlnode(ib2))
	  else
	    wtb = 1.0
	  endif
     	  albedo = albedo + wtb*wtflx(ib)*albb(ib)
	  twtb = twtb + wtb*wtflx(ib)
	enddo
	if (twtb.eq.0. .and. ib1.eq.ib2)then
	   albedo = albb(ib1)
	else
	   albedo = albedo/twtb
	endif


fuspecalb(ifu) = albedo
enddo FUBANDS

return

end subroutine zjin_wgtlut

!=======================================================================
      subroutine locate(xx,n,x,j)
!
! purpose:  given an array xx of length n, and given a value X, returns
!           a value J such that X is between xx(j) and xx(j+1). xx must
!           be monotonic, either increasing of decreasing. this function
!           returns j=1 or j=n-1 if x is out of range.
!
! input:
!   xx      monitonic table
!   n       size of xx
!   x       single floating point value perhaps within the range of xx
!
      integer j,n
      real*4 x,xx(n)
      integer jl,jm,ju

      if(x.eq.xx(1)) then
        j=1
        return
      endif
      if(x.eq.xx(n)) then
        j=n-1
        return
      endif
      jl=1
      ju=n
10    if(ju-jl.gt.1) then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      end subroutine locate
!=======================================================================



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           FILE: SUBR_native_4byte_real.f90
!     SUBPROGRAM: native_4byte_real
!
!         AUTHOR: David Stepaniak, NCAR/CGD/CAS
! DATE INITIATED: 29 April 2003 
!  LAST MODIFIED: 29 April 2003
!
!       SYNOPSIS: Converts a 32 bit, 4 byte, REAL from big Endian to
!                 little Endian, or conversely from little Endian to big
!                 Endian.
!
!    DESCRIPTION: This subprogram allows one to convert a 32 bit, 4 byte,
!                 REAL data element that was generated with, say, a big
!                 Endian processor (e.g. Sun/sparc, SGI/R10000, etc.) to its
!                 equivalent little Endian representation for use on little
!                 Endian processors (e.g. PC/Pentium running Linux). The
!                 converse, little Endian to big Endian, also holds.
!                 This conversion is accomplished by writing the 32 bits of
!                 the REAL data element into a generic 32 bit INTEGER space
!                 with the TRANSFER intrinsic, reordering the 4 bytes with
!                 the MVBITS intrinsic, and writing the reordered bytes into
!                 a new 32 bit REAL data element, again with the TRANSFER
!                 intrinsic. The following schematic illustrates the
!                 reordering process
!
!
!                  --------    --------    --------    --------
!                 |    D   |  |    C   |  |    B   |  |    A   |  4 Bytes
!                  --------    --------    --------    --------
!                                                             |
!                                                              -> 1 bit
!                                       ||
!                                     MVBITS
!                                       ||
!                                       \/
!
!                  --------    --------    --------    --------
!                 |    A   |  |    B   |  |    C   |  |    D   |  4 Bytes
!                  --------    --------    --------    --------
!                         |           |           |           |
!                         24          16          8           0   <- bit
!                                                                 position
!
!          INPUT: realIn,  a single 32 bit, 4 byte REAL data element.
!         OUTPUT: realOut, a single 32 bit, 4 byte REAL data element, with
!                 reverse byte order to that of realIn.
!    RESTRICTION: It is assumed that the default REAL data element is
!                 32 bits / 4 bytes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE native_4byte_real( realIn, realOut )
  IMPLICIT NONE
  REAL*4, INTENT(IN)                              :: realIn
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element
  REAL*4, INTENT(OUT)                             :: realOut
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element, with
                                                   ! reverse byte order to
                                                   ! that of realIn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Local variables (generic 32 bit INTEGER spaces):

  INTEGER                                       :: i_element
  INTEGER                                       :: i_element_br
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transfer 32 bits of realIn to generic 32 bit INTEGER space:
  i_element = TRANSFER( realIn, 0 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reverse order of 4 bytes in 32 bit INTEGER space:
  CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
  CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
  CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
  CALL MVBITS( i_element,  0, 8, i_element_br, 24 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transfer reversed order bytes to 32 bit REAL space (realOut):
  realOut = TRANSFER( i_element_br, 0.0 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ZJIN_OCEAN_WGTLUT

