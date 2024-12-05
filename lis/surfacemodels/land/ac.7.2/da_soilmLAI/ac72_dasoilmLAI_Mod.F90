!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module ac72_dasoilmLAI_Mod
!BOP
!
! !MODULE: ac72_dasoilmLAI_Mod
!
! !DESCRIPTION:
!  
! !REVISION HISTORY:

! Sujay Kumar; Initial Code
! 9 Sep 2016: Mahdi Navari; Modified for ac72 !
! 18 Jun 2021: Michel Bechtold: SM and LAI updating with S1 backscatter w/ WCM
! !USES:        
  use ESMF
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_logMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ac72_dasoilmLAI_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: ac72_dasm_struc
!EOP

 type, public :: dasm_dec
     real,    allocatable       :: model_xrange(:,:,:)
     real,    allocatable       :: model_cdf(:,:,:)
     real,    allocatable       :: model_mu(:)

     integer                :: nbins
     integer                :: ntimes
     integer                :: scal

  end type dasm_dec
  
  type(dasm_dec), allocatable :: ac72_dasm_struc(:)

contains
!BOP
! 
! !ROUTINE: ac72_dasoilmLAI_init
! \label{ac72_dasoilmLAI_init}
! 
! !INTERFACE:
  subroutine ac72_dasoilmLAI_init(k)
! !USES:
! !DESCRIPTION:        
!
!EOP
    

    implicit none
    integer                :: k
    integer                :: n 
    character*100          :: modelcdffile(LIS_rc%nnest)
    integer                :: status
    integer                :: ngrid

    if(.not.allocated(ac72_dasm_struc)) then 
       allocate(ac72_dasm_struc(LIS_rc%nnest))
    endif
    
!TBD: SVK
#if 0 
    if(LIS_rc%dascaloption(k).eq."Linear scaling") then 
       call ESMF_ConfigFindLabel(LIS_config,"ac72 soil moisture CDF file:",&
            rc=status)
       do n=1,LIS_rc%nnest
          call ESMF_ConfigGetAttribute(LIS_config,modelcdffile(n),rc=status)
          call LIS_verify(status, 'ac72 soil moisture CDF file: not defined')
       enddo
       
       do n=1,LIS_rc%nnest
       
!Hardcoded for now.
          ac72_dasm_struc(n)%nbins = 100
          
          call LIS_getCDFattributes(modelcdffile(n),&
               ac72_dasm_struc(n)%ntimes, ngrid)
          
          allocate(ac72_dasm_struc(n)%model_xrange(&
               LIS_rc%ngrid(n), ac72_dasm_struc(n)%ntimes, &
               ac72_dasm_struc(n)%nbins))
          allocate(ac72_dasm_struc(n)%model_cdf(&
               LIS_rc%ngrid(n), ac72_dasm_struc(n)%ntimes, &
               ac72_dasm_struc(n)%nbins))
          
          call LIS_readCDFdata(n,&
               ac72_dasm_struc(n)%nbins, &
               ac72_dasm_struc(n)%ntimes, &
               ngrid, &
               modelcdffile(n), &
               "SoilMoist",&
               ac72_dasm_struc(n)%model_xrange,&
               ac72_dasm_struc(n)%model_cdf)
       enddo
    endif
#endif

  end subroutine ac72_dasoilmLAI_init
end module ac72_dasoilmLAI_Mod
