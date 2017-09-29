      SUBROUTINE FNAMES(CLRT,PLTOUT,PLTSCN,SPFLUX,DBOUT,                &
     &  M3DGEN,M3DCNT,M3DDAT,M3DMOL,LNFLRT,FLRT,FILE_NAME,IN_LENGTH)

!     THIS ROUTINE ASSIGNS FILE NAMES TO THE ABOVE CHARACTER VARIABLES,

!     ARGUMENTS:
!       M3DGEN   NAME OF GENERAL INFORMATION FILE FOR MOD3D.
!       M3DCNT   NAME OF CONTINUA DATA FILE FOR MOD3D.
!       M3DDAT   NAME OF MOLECULAR EXTINCTION FILE SORTED BY ATMOSPHERE.
!       M3DMOL   NAME OF MOLECULAR EXTINCTION DATA FILE FOR MOD3D.
      CHARACTER*(*) CLRT,PLTOUT,PLTSCN,SPFLUX,DBOUT,                    &
     &  M3DGEN,M3DCNT,M3DDAT,M3DMOL
      CHARACTER*(120) FILE_NAME
      INTEGER IN_LENGTH

!     INCLUDE 'CHANLS.h' FOR DECLARATION OF CHNOUT:
      INCLUDE 'PARAMS.h'
      INCLUDE 'CHANLS.h'
      INCLUDE 'ERROR.h'

!     COMMONS:
      INCLUDE 'IFIL.h'

!     DECLARE BLOCK DATA ROUTINES EXTERNAL:
      EXTERNAL DEVCBD

!     LOCAL VARIABLES:
      CHARACTER FLRT*(NAMLEN-4)
      INTEGER LUNIT,NUNIT,LNFLRT,LENSTR
      LOGICAL LROOT,LXIST7
      LOGICAL DEBUG_FLAG

!     DATA:
!       CFRMT    CHARACTER STRING INPUT FORMAT.
      CHARACTER*9 CFRMT
      SAVE CFRMT
      DATA CFRMT/'((A    ))'/

!     FOR THE JMASS API OPTION, JUST OPEN 'tape6' FOR
!     WARNING MESSAGES AND RETURN.
	DEBUG_FLAG = .FALSE.
      IF(LJMASS)THEN
          CALL OPNFL(IPR,'tape6','UNKNOWN')
          M3DGEN='mod3d.g3d'
          M3DCNT='mod3d.c3d'
          M3DDAT='mod3d.d3d'
          M3DMOL='mod3d.m3d'
          RETURN
      ENDIF
      LNFLRT=0
      LUNIT=NUNIT()
      INQUIRE(FILE='modroot.in',EXIST=LROOT)
      IF(LROOT)THEN
          CALL OPNFL(LUNIT,'modroot.in','OLD')
      ELSE
          INQUIRE(FILE='MODROOT.IN',EXIST=LROOT)
          IF (LROOT) CALL OPNFL(LUNIT,'MODROOT.IN','OLD')
      ENDIF
      IF(LROOT)THEN
          WRITE(CFRMT(4:7),'(I4.4)')LEN(FLRT)
          READ(LUNIT,CFRMT)FLRT
          CLOSE(LUNIT)
          LNFLRT=LENSTR(FLRT)
          !DRF
          LNFLRT = IN_LENGTH 
          FLRT = FILE_NAME
          !WRITE(*,*) 'filename = ',FILE_NAME
      ENDIF
      IF(.NOT.LROOT.OR.LNFLRT.EQ.0)THEN
          CALL OPNFL(IRD,'tape5','OLD')
          CALL OPNFL(IPR,'tape6','UNKNOWN')
          CALL OPNFL(IPU,'tape7','UNKNOWN')
          CALL OPNFL(IPR1,'tape8','UNKNOWN')
          CALL OPNFL(IDBOUT,'mc.dat','UNKNOWN')
          PLTOUT='pltout'
          CLRT='clrates'
          CALL OPNFL(IPUSCN,'tape7.scn','UNKNOWN')
          PLTSCN='pltout.scn'
          SPFLUX='specflux'
          CHNOUT='channels.out'
          LENCHN=12
          M3DGEN='mod3d.g3d'
          M3DCNT='mod3d.c3d'
          M3DDAT='mod3d.d3d'
          M3DMOL='mod3d.m3d'
          INQUIRE(FILE='tape7.scr',EXIST=LXIST7)
          IF(.NOT.LXIST7)THEN
              CALL OPNFL(IPUSCR,'tape7.scr','NEW')
          ELSE
              CALL OPNFL(IPUSCR,'tape7.scr','OLD')
          ENDIF
      ELSE
	IF (DEBUG_FLAG) THEN
          CALL OPNFL(IRD,FLRT(1:LNFLRT)//'.tp5','OLD')
          CALL OPNFL(IPR,FLRT(1:LNFLRT)//'.tp6','UNKNOWN')
          CALL OPNFL(IPU,FLRT(1:LNFLRT)//'.tp7','UNKNOWN')
          CALL OPNFL(IPR1,FLRT(1:LNFLRT)//'.tp8','UNKNOWN')
!SSI      CALL OPNFL(IDBIN,FLRT(1:LNFLRT)//'.dbi','UNKNOWN')
          CALL OPNFL(IDBOUT,FLRT(1:LNFLRT)//'.mc','UNKNOWN')
          DBOUT=FLRT(1:LNFLRT)//'.dbo'
          PLTOUT=FLRT(1:LNFLRT)//'.plt'
          CLRT=FLRT(1:LNFLRT)//'.clr'
          CALL OPNFL(IPUSCN,FLRT(1:LNFLRT)//'.7sc','UNKNOWN')
          PLTSCN=FLRT(1:LNFLRT)//'.psc'
          SPFLUX=FLRT(1:LNFLRT)//'.flx'
          CHNOUT=FLRT(1:LNFLRT)//'.chn'
          LENCHN=LNFLRT+4
          M3DGEN=FLRT(1:LNFLRT)//'.g3d'
          M3DCNT=FLRT(1:LNFLRT)//'.c3d'
          M3DDAT=FLRT(1:LNFLRT)//'.d3d'
          M3DMOL=FLRT(1:LNFLRT)//'.m3d'
          INQUIRE(FILE=FLRT(1:LNFLRT)//'.7sr',EXIST=LXIST7)
          IF(.NOT.LXIST7)THEN
              CALL OPNFL(IPUSCR,FLRT(1:LNFLRT)//'.7sr','NEW')
          ELSE
              CALL OPNFL(IPUSCR,FLRT(1:LNFLRT)//'.7sr','OLD')
          ENDIF
	ELSE
          !CALL OPNFL(IRD,FLRT(1:LNFLRT)//'.tp5','UNKNOWN')
          !CALL OPNFL(IPR,FLRT(1:LNFLRT)//'.tp6','NEW')
          
          CALL OPNFL(IPR,'/dev/null','OLD')
          CALL OPNFL(IPU,'/dev/null','OLD')
          CALL OPNFL(IPR1,'/dev/null','OLD')
          CALL OPNFL(IDBOUT,'/dev/null','OLD')
          DBOUT=FLRT(1:LNFLRT)//'.dbo'
          PLTOUT=FLRT(1:LNFLRT)//'.plt'
          CLRT=FLRT(1:LNFLRT)//'.clr'
          !CALL OPNFL(IPUSCN,FLRT(1:LNFLRT)//'.7sc','UNKNOWN')
          !DRF CALL OPNFL(IPUSCN,FLRT(1:LNFLRT)//'.7sc','NEW')
          CALL OPNFL(IPUSCN,'/dev/null','OLD')
          PLTSCN=FLRT(1:LNFLRT)//'.psc'
          SPFLUX=FLRT(1:LNFLRT)//'.flx'
          CHNOUT=FLRT(1:LNFLRT)//'.chn'
          LENCHN=LNFLRT+4
          M3DGEN=FLRT(1:LNFLRT)//'.g3d'
          M3DCNT=FLRT(1:LNFLRT)//'.c3d'
          M3DDAT=FLRT(1:LNFLRT)//'.d3d'
          M3DMOL=FLRT(1:LNFLRT)//'.m3d'
          INQUIRE(FILE='/tmp/'//FLRT(1:LNFLRT)//'.7sr',EXIST=LXIST7)
          IF(.NOT.LXIST7)THEN
              CALL OPNFL(IPUSCR,'/tmp/'//FLRT(1:LNFLRT)//'.7sr','NEW')
          ELSE
              CALL OPNFL(IPUSCR,'/tmp/'//FLRT(1:LNFLRT)//'.7sr','OLD')
          ENDIF
	ENDIF
      ENDIF
      RETURN
      END
