      SUBROUTINE INIT_MODPCRTM(FILE_NAME,FLENGTH,TP5_FLAG,TP5_FILE,     &
     &   NUM_LEVS                                                       &
     &                  ,WVLNGTH_LRES,RADIANCE_LRES,WV_HRES,            &
     &                  RADIANCE_HRES,NUM_LRES,NUM_HRES,NUM_AC_WVL,     &
     &                  SOLAR_FLUX,DIFFUSE_FLUX,                        &
     &                  BB_UPDIFFUSE,BB_DNDIFFUSE,BB_DNDIRECT,          &
     &                  FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)

            use mpi

      implicit none
      logical PASS1
      common /GTBRDFC/PASS1

      integer ierr,num_procs,id

!BEGIN DRF MODIFICATIONS
!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER NUM_LEVS
      INTEGER NUM_LRES,NUM_HRES,NUM_AC_WVL
      LOGICAL TP5_FLAG
      CHARACTER TP5_FILE(1000)*110
      CHARACTER DUMMY*80
      INTEGER TP5_LINENUMBER
      CHARACTER*120 FILE_NAME
      INTEGER FLENGTH
      REAL WVLNGTH_LRES(NUM_LRES), RADIANCE_LRES(NUM_LRES)
      REAL WV_HRES(NUM_HRES), RADIANCE_HRES(NUM_HRES)
      REAL SOLAR_FLUX(NUM_LRES),DIFFUSE_FLUX(NUM_LRES)
      REAL BB_UPDIFFUSE(NUM_LEVS)
      REAL BB_DNDIFFUSE(NUM_LEVS) 
      REAL BB_DNDIRECT(NUM_LEVS)
      REAL FLX_UPDIFFUSE(NUM_LEVS,NUM_LRES)
      REAL FLX_DNDIFFUSE(NUM_LEVS,NUM_LRES) 
      REAL FLX_DNDIRECT(NUM_LEVS,NUM_LRES)
      INTEGER CLD_INDEX_LOW   !FOR CTHIK CALCULATIONS
      INTEGER CLD_INDEX_HIGH   !FOR CTHIK CALCULATIONS
      INTEGER CLD_INDEX            !FOR CTHIK CALCULATIONS
!!!!!!!!!!!!!!!!!!!!!!!!
      !END DRF MODIFICATIONS

!      call mpi_init(ierr)
         call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
!         call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
        pass1 = .true.

      CALL DRIVER(ID,FILE_NAME,FLENGTH,TP5_FLAG,TP5_FILE,NUM_LEVS &
     &                  ,WVLNGTH_LRES,RADIANCE_LRES,WV_HRES,            &
     &                  RADIANCE_HRES,NUM_LRES,NUM_HRES,NUM_AC_WVL,     &
     &                  SOLAR_FLUX,DIFFUSE_FLUX,                        &
     &                  BB_UPDIFFUSE,BB_DNDIFFUSE,BB_DNDIRECT,          &
     &                  FLX_UPDIFFUSE,FLX_DNDIFFUSE,FLX_DNDIRECT)
!      call mpi_finalize(ierr)

      return
      end
