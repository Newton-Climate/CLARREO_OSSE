!DRF!
!DRF!! Interpolate ozone volume mixing ratio to model levels
!DRF!! !
!DRF !!   if (build_ozone) then
!DRF !!      call radozn(lchnk   ,ncol    ,pmid    ,o3vmr   )
!DRF !!   endif
!DRF !!   call outfld('O3VMR   ',o3vmr ,pcols, lchnk)

!DRF!!
!!DRF!!
!!DRF!! Set chunk dependent radiation input
!!DRF!!
 !DRF!!  call radinp(lchnk   ,ncol    ,                                &
!DRF!!              pmid    ,pint    ,o3vmr   , pbr     ,&
!DRF!!#ifndef OFFLINE
 !DRF!!                             pnm     ,eccf    ,o3mmr   )
!DRF!!#else
 !DRF!!              pnm     ,eccf    ,o3mmr   , calday)
!DRF!!
 !DRF!!  if (sol_ann_mean) then
 !DRF!!     eccf = 1.0_r8
 !DRF!!  endif
!DRF!!#endif

!DRF      call aqsat(state%t, state%pmid, esat, qsat, pcols, &
!DRF                 ncol, pver, 1, pver)

      ! calculate relative humidity
!DRF      rh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / qsat(1:ncol,1:pver) * &
!DRF         ((1.0 - epsilo) * qsat(1:ncol,1:pver) + epsilo) / &
!DRF         ((1.0 - epsilo) * state%q(1:ncol,1:pver,1) + epsilo)


 !DRF write(*,*) "SW:lchnk = ", lchnk
 !DRF write(*,*) "SW:ncol = ", ncol
 !DRF write(*,*) "SW:pnm = ", pnm 
 !DRF write(*,*) "SW:pbr = ", pbr 
 !DRF write(*,*) "SW:qm1(1,1,1)= ", qm1(1,1,1)
 !DRF write(*,*) "SW:rh = ", rh 
 !DRF write(*,*) "SW:o3mmr = ", o3mmr 
 !DRF write(*,*) "SW:aerosol(1,1,1)= ", aerosol(1,1,1)
 !DRF write(*,*) "SW:cld = ", cld 
 !DRF write(*,*) "SW:cicewp = ", cicewp 
 !DRF write(*,*) "SW:cliqwp = ", cliqwp 
 !DRF write(*,*) "SW:rel = ", rel 
 !DRF write(*,*) "SW:rei=  ", rei 
 !DRF write(*,*) "SW:eccf= ", eccf
 !DRF write(*,*) "SW:coszrs(1)= ", coszrs(1)
 !DRF write(*,*) "SW:scon= ", scon
 !DRF write(*,*) "SW:solin(1)= ", solin(1)
 !DRF write(*,*) "SW:asdir= ", asdir
 !DRF write(*,*) "SW:adif= ", asdif
 !DRF write(*,*) "SW:aldir= ", aldir
 !DRF write(*,*) "SW:aldif= ", aldif
 !DRF write(*,*) "SW:nmxrgnrf= ", nmxrgnrf
 !DRF write(*,*) "SW:pmxrgnrf= ", pmxrgnrf
 !DRF write(*,*) "SW:qrs= ", qrs
 !DRF write(*,*) "SW:fsnt= ", fsnt
 !DRF write(*,*) "SW:fsntc= ", fsntc
 !DRF write(*,*) "SW:fsntoa= ", fsntoa
 !DRF write(*,*) "SW:fsntoac= ", fsntoac
 !DRF write(*,*) "SW:fsnirt= ", fsnirt
 !DRF write(*,*) "SW:fsnrtc= ", fsnrtc
 !DRF write(*,*) "SW:fsnirtsq= ", fsnirtsq
 !DRF write(*,*) "SW:fsns= ", fsns
 !DRF write(*,*) "SW:sols= ", sols
 !DRF write(*,*) "SW:soll= ",soll
 !DRF write(*,*) "SW:solsd= ",solsd
 !DRF write(*,*) "SW:solld= ",solld
 !DRF write(*,*) "SW:frc_day= ",frc_day
 !DRF write(*,*) "SW:aertau= ",aertau
 !DRF write(*,*) "SW:aerssa= ",aerssa
 !DRF write(*,*) "SW:aerasm= ",aerasm
 !DRF write(*,*) "SW:aerfwd= ",aerfwd
 !DRF write(*,*) "SW:fsn= ",fsn
 !DRF write(*,*) "SW:fsnc= ",fsnc
!DRF write(*,*) "SW:no_o2_abs= ",no_o2_abs
!write(*,*) "calling radcswmx 1st time"



!DRF         call cnst_get_ind('N2O'  , in2o)
!DRF         call cnst_get_ind('CH4'  , ich4)
!DRF         call cnst_get_ind('CFC11', if11)
!DRF         call cnst_get_ind('CFC12', if12)
!DRF         call t_startf("radclwmx")
! DRF write(*,*) "in radctl, o3vmr(1,1) = ",o3vmr(1,1)


!DRF variables
!write(*,*) "lchnk = ", lchnk
!write(*,*) "ncol = ", ncol
!write(*,*) "lwupcgs= ",lwupcgs 
!write(*,*) "t= ", t
!write(*,*) "qm1(1,1,1)= ",qm1(1,1,1)
!write(*,*) "o3vmr= ", o3vmr
!write(*,*) "pbr_lw= ", pbr
!write(*,*) "pnm= ", pnm
!write(*,*) "pmln= ", pmln
!write(*,*) "piln= ", piln
!write(*,*) "qm1(1,1,in2o)= ", qm1(1,1,in2o)
!write(*,*) "qm1(1,1,ich4)= ", qm1(1,1,ich4)
!write(*,*) "qm1(1,1,if11)= ", qm1(1,1,if11)
!write(*,*) "qm1(1,1,if12)= ", qm1(1,1,if12)
!write(*,*) "cld= ", cld
!write(*,*) "emis= ", emis
!write(*,*) "pmxrgn= ", pmxrgn
!write(*,*) "nmxrgn= ", nmxrgn
!write(*,*) "qrl= ", qrl
!write(*,*) "flns= ", flns
!write(*,*) "flnt= ", flnt
!write(*,*) "flnsc= ", flnsc
!write(*,*) "flntc= ", flntc
!write(*,*) "flwds= ", flwds
!write(*,*) "flut= ", flut
!write(*,*) "flutc= ", flutc
!write(*,*) "aerosl(:,:,idxVOLC)= ", aerosol(:,:,idxVOLC)
!write(*,*) "fln= ", fln
!write(*,*) "flnc= ", flnc
