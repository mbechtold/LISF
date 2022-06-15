!-----------------------BEGIN---------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module WCMRTM_multi_Mod
!BOP
!
! !MODULE: WCMRTM_multi_Mod
!
! !DESCRIPTION:
!    This module provides the routines to control the execution of 
!    the Water Cloud model (Ulaby et al., 1990). The routine allows to compute
!    both the backscatter in VV and the backscatter in VH. It is written
!    considering a static incidence angle of 37°. This can be improved in cas
!    a local incidence angle for each observation is available. Table of
!    calibrated parameters (different for the two polarizations are added also
!    in the configuration file). Each pixel of the study area needs to have
!    assigned A,B,C,D parameters for each polarization, associated with
!    longitude and latitude of the same pixel (pixel where the soil moisture is
!    equa to nan need to be removed)
!
! !HISTORY:
! 28 Aug 2020: Sara Modanesi
! 26 Mar 2021 Sara Modanesi: Added specifications for forward states
! !USES:        


#if (defined RTMS)

  use ESMF
  use LIS_coreMod
  use LIS_RTMMod
  use LIS_logMod

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: WCMRTM_multi_initialize
  public :: WCMRTM_multi_f2t
  public :: WCMRTM_multi_run
  public :: WCMRTM_multi_output
  public :: WCMRTM_multi_geometry 
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: wcm_struc
!EOP
  type, public ::  wcm_type_dec 
   
     character(len=256) :: AA_VV_037_D_tbl_name !! from NoahMP36_lsmMod.F90
     character(len=256) :: BB_VV_037_D_tbl_name
     character(len=256) :: CC_VV_037_D_tbl_name
     character(len=256) :: DD_VV_037_D_tbl_name

     character(len=256) :: AA_VH_037_D_tbl_name
     character(len=256) :: BB_VH_037_D_tbl_name
     character(len=256) :: CC_VH_037_D_tbl_name
     character(len=256) :: DD_VH_037_D_tbl_name

     real, allocatable :: AA_VV_037_D(:)
     real, allocatable :: BB_VV_037_D(:)
     real, allocatable :: CC_VV_037_D(:)
     real, allocatable :: DD_VV_037_D(:)

     real, allocatable :: AA_VH_037_D(:)
     real, allocatable :: BB_VH_037_D(:)
     real, allocatable :: CC_VH_037_D(:)
     real, allocatable :: DD_VH_037_D(:)

     real, allocatable :: lone(:)
     real, allocatable :: late(:)
     !-------output------------!   
     real, allocatable :: Sig0VV_037_D(:)
     real, allocatable :: Sig0VH_037_D(:)

     ! -------------------- 088_A ----
     character(len=256) :: AA_VV_088_A_tbl_name !! from NoahMP36_lsmMod.F90
     character(len=256) :: BB_VV_088_A_tbl_name
     character(len=256) :: CC_VV_088_A_tbl_name
     character(len=256) :: DD_VV_088_A_tbl_name

     character(len=256) :: AA_VH_088_A_tbl_name
     character(len=256) :: BB_VH_088_A_tbl_name
     character(len=256) :: CC_VH_088_A_tbl_name
     character(len=256) :: DD_VH_088_A_tbl_name

     real, allocatable :: AA_VV_088_A(:)
     real, allocatable :: BB_VV_088_A(:)
     real, allocatable :: CC_VV_088_A(:)
     real, allocatable :: DD_VV_088_A(:)

     real, allocatable :: AA_VH_088_A(:)
     real, allocatable :: BB_VH_088_A(:)
     real, allocatable :: CC_VH_088_A(:)
     real, allocatable :: DD_VH_088_A(:)

     !-------output------------!   
     real, allocatable :: Sig0VV_088_A(:)
     real, allocatable :: Sig0VH_088_A(:)

     ! -------------------- 110_D ----
     character(len=256) :: AA_VV_110_D_tbl_name !! from NoahMP36_lsmMod.F90
     character(len=256) :: BB_VV_110_D_tbl_name
     character(len=256) :: CC_VV_110_D_tbl_name
     character(len=256) :: DD_VV_110_D_tbl_name

     character(len=256) :: AA_VH_110_D_tbl_name
     character(len=256) :: BB_VH_110_D_tbl_name
     character(len=256) :: CC_VH_110_D_tbl_name
     character(len=256) :: DD_VH_110_D_tbl_name

     real, allocatable :: AA_VV_110_D(:)
     real, allocatable :: BB_VV_110_D(:)
     real, allocatable :: CC_VV_110_D(:)
     real, allocatable :: DD_VV_110_D(:)

     real, allocatable :: AA_VH_110_D(:)
     real, allocatable :: BB_VH_110_D(:)
     real, allocatable :: CC_VH_110_D(:)
     real, allocatable :: DD_VH_110_D(:)

     !-------output------------!   
     real, allocatable :: Sig0VV_110_D(:)
     real, allocatable :: Sig0VH_110_D(:)

     ! -------------------- 139_D ----
     character(len=256) :: AA_VV_139_D_tbl_name !! from NoahMP36_lsmMod.F90
     character(len=256) :: BB_VV_139_D_tbl_name
     character(len=256) :: CC_VV_139_D_tbl_name
     character(len=256) :: DD_VV_139_D_tbl_name

     character(len=256) :: AA_VH_139_D_tbl_name
     character(len=256) :: BB_VH_139_D_tbl_name
     character(len=256) :: CC_VH_139_D_tbl_name
     character(len=256) :: DD_VH_139_D_tbl_name

     real, allocatable :: AA_VV_139_D(:)
     real, allocatable :: BB_VV_139_D(:)
     real, allocatable :: CC_VV_139_D(:)
     real, allocatable :: DD_VV_139_D(:)

     real, allocatable :: AA_VH_139_D(:)
     real, allocatable :: BB_VH_139_D(:)
     real, allocatable :: CC_VH_139_D(:)
     real, allocatable :: DD_VH_139_D(:)

     !-------output------------!   
     real, allocatable :: Sig0VV_139_D(:)
     real, allocatable :: Sig0VH_139_D(:)

     ! -------------------- 161_A ----
     character(len=256) :: AA_VV_161_A_tbl_name !! from NoahMP36_lsmMod.F90
     character(len=256) :: BB_VV_161_A_tbl_name
     character(len=256) :: CC_VV_161_A_tbl_name
     character(len=256) :: DD_VV_161_A_tbl_name

     character(len=256) :: AA_VH_161_A_tbl_name
     character(len=256) :: BB_VH_161_A_tbl_name
     character(len=256) :: CC_VH_161_A_tbl_name
     character(len=256) :: DD_VH_161_A_tbl_name

     real, allocatable :: AA_VV_161_A(:)
     real, allocatable :: BB_VV_161_A(:)
     real, allocatable :: CC_VV_161_A(:)
     real, allocatable :: DD_VV_161_A(:)

     real, allocatable :: AA_VH_161_A(:)
     real, allocatable :: BB_VH_161_A(:)
     real, allocatable :: CC_VH_161_A(:)
     real, allocatable :: DD_VH_161_A(:)

     !-------output------------!   
     real, allocatable :: Sig0VV_161_A(:)
     real, allocatable :: Sig0VH_161_A(:)

  end type wcm_type_dec

  type(wcm_type_dec), allocatable :: wcm_struc(:) 

  SAVE

contains
!BOP
! 
! !ROUTINE: WCMRTM_multi_initialize
! \label{WCMRTM_multi_initialize}
! 
! !INTERFACE:
   subroutine WCMRTM_multi_initialize()
! !USES:

! !DESCRIPTION:        
!
!  This routine creates the datatypes and allocates memory for noahMP3.6-specific
!  variables. It also invokes the routine to read the runtime specific options
!  for noahMP3.6 from the configuration file. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readWCMRTM_multicrd](\ref{readWCMRTM_multicrd}) \newline
!
!EOP
   implicit none
    
   integer :: rc
   integer :: n,t
   integer :: ierr
   integer , parameter :: OPEN_OK = 0
   character*128 :: message

!allocate memory for nest
   allocate(wcm_struc(LIS_rc%nnest))

   do n=1,LIS_rc%nnest
!allocate memory for all tile in current nest


      allocate(wcm_struc(n)%AA_VV_037_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VV_037_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VV_037_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VV_037_D(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%AA_VH_037_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VH_037_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VH_037_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VH_037_D(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%lone(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%late(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%Sig0VV_037_D(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      allocate(wcm_struc(n)%Sig0VH_037_D(LIS_rc%npatch(n,LIS_rc%lsm_index)))

      ! ---------------------- 088_A
      allocate(wcm_struc(n)%AA_VV_088_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VV_088_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VV_088_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VV_088_A(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%AA_VH_088_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VH_088_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VH_088_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VH_088_A(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%Sig0VV_088_A(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      allocate(wcm_struc(n)%Sig0VH_088_A(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      ! ----------------------
      ! ---------------------- 110_D
      allocate(wcm_struc(n)%AA_VV_110_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VV_110_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VV_110_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VV_110_D(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%AA_VH_110_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VH_110_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VH_110_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VH_110_D(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%Sig0VV_110_D(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      allocate(wcm_struc(n)%Sig0VH_110_D(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      ! ----------------------
      ! ---------------------- 139_D
      allocate(wcm_struc(n)%AA_VV_139_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VV_139_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VV_139_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VV_139_D(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%AA_VH_139_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VH_139_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VH_139_D(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VH_139_D(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%Sig0VV_139_D(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      allocate(wcm_struc(n)%Sig0VH_139_D(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      ! ----------------------
      ! ---------------------- 161_A
      allocate(wcm_struc(n)%AA_VV_161_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VV_161_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VV_161_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VV_161_A(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%AA_VH_161_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%BB_VH_161_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%CC_VH_161_A(LIS_rc%glbngrid(n)))
      allocate(wcm_struc(n)%DD_VH_161_A(LIS_rc%glbngrid(n)))

      allocate(wcm_struc(n)%Sig0VV_161_A(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      allocate(wcm_struc(n)%Sig0VH_161_A(LIS_rc%npatch(n,LIS_rc%lsm_index)))
      ! ----------------------

      call add_sfc_fields(n,LIS_sfcState(n), "Soil Moisture Content")
      call add_sfc_fields(n,LIS_sfcState(n), "Leaf Area Index")

   enddo 

!----------------A,B,C and D parameter tables for VV_037_D pol---------------------!
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VV_037_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest 
      call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VV_037_D_tbl_name, rc=rc)
      call LIS_verify(rc, "WCMRTM_multi AA_VV_037_D parameter table: not defined")
   enddo
   

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VV_037_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest   
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VV_037_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VV_037_D parameter table: not defined")
   enddo
    
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VV_037_D parameter table:",rc = rc)    
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VV_037_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VV_037_D parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VV_037_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VV_037_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VV_037_D parameter table: not defined")
   enddo



!----------------A,B,C and D parameter tables for VH_037_D pol---------------------!
   
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VH_037_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VH_037_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi AA_VH_037_D parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VH_037_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VH_037_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VH_037_D parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VH_037_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VH_037_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VH_037_D parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VH_037_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VH_037_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VH_037_D parameter table: not defined")
   enddo
!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VV_037_D pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VV_037_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VV_037_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VV_037_D(t),wcm_struc(n)%lone(t),&
          wcm_struc(n)%late(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VV_037_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VV_037_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VV_037_D_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VV_037_D(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VV_037_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VV_037_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VV_037_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VV_037_D(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VV_037_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VV_037_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VV_037_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VV_037_D(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo


!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VH_037_D pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VH_037_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VH_037_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VH_037_D(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VH_037_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VH_037_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VH_037_D_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VH_037_D(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VH_037_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VH_037_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VH_037_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VH_037_D(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VH_037_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VH_037_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VH_037_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VH_037_D(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo

   do n=1,LIS_rc%nnest !added fields to State 26032021        
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VV_037_D")
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VH_037_D")
   enddo


!---------------------------------------------------------------------------------------------
!----------------A,B,C and D parameter tables for VV_088_A pol---------------------!
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VV_088_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest 
      call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VV_088_A_tbl_name, rc=rc)
      call LIS_verify(rc, "WCMRTM_multi AA_VV_088_A parameter table: not defined")
   enddo
   

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VV_088_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest   
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VV_088_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VV_088_A parameter table: not defined")
   enddo
    
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VV_088_A parameter table:",rc = rc)    
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VV_088_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VV_088_A parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VV_088_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VV_088_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VV_088_A parameter table: not defined")
   enddo



!----------------A,B,C and D parameter tables for VH_088_A pol---------------------!
   
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VH_088_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VH_088_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi AA_VH_088_A parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VH_088_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VH_088_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VH_088_A parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VH_088_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VH_088_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VH_088_A parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VH_088_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VH_088_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VH_088_A parameter table: not defined")
   enddo
!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VV_088_A pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VV_088_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VV_088_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VV_088_A(t),wcm_struc(n)%lone(t),&
          wcm_struc(n)%late(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VV_088_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VV_088_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VV_088_A_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VV_088_A(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VV_088_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VV_088_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VV_088_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VV_088_A(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VV_088_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VV_088_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VV_088_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VV_088_A(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo


!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VH_088_A pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VH_088_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VH_088_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VH_088_A(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VH_088_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VH_088_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VH_088_A_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VH_088_A(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VH_088_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VH_088_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VH_088_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VH_088_A(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VH_088_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VH_088_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VH_088_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VH_088_A(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo

   do n=1,LIS_rc%nnest !added fields to State 26032021        
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VV_088_A")
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VH_088_A")
   enddo

!---------------------------------------------------------------------------------------------
!----------------A,B,C and D parameter tables for VV_110_D pol---------------------!
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VV_110_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest 
      call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VV_110_D_tbl_name, rc=rc)
      call LIS_verify(rc, "WCMRTM_multi AA_VV_110_D parameter table: not defined")
   enddo
   

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VV_110_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest   
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VV_110_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VV_110_D parameter table: not defined")
   enddo
    
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VV_110_D parameter table:",rc = rc)    
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VV_110_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VV_110_D parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VV_110_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VV_110_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VV_110_D parameter table: not defined")
   enddo



!----------------A,B,C and D parameter tables for VH_110_D pol---------------------!
   
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VH_110_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VH_110_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi AA_VH_110_D parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VH_110_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VH_110_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VH_110_D parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VH_110_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VH_110_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VH_110_D parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VH_110_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VH_110_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VH_110_D parameter table: not defined")
   enddo
!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VV_110_D pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VV_110_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VV_110_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VV_110_D(t),wcm_struc(n)%lone(t),&
          wcm_struc(n)%late(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VV_110_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VV_110_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VV_110_D_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VV_110_D(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VV_110_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VV_110_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VV_110_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VV_110_D(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VV_110_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VV_110_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VV_110_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VV_110_D(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo


!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VH_110_D pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VH_110_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VH_110_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VH_110_D(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VH_110_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VH_110_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VH_110_D_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VH_110_D(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VH_110_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VH_110_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VH_110_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VH_110_D(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VH_110_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VH_110_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VH_110_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VH_110_D(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo

   do n=1,LIS_rc%nnest !added fields to State 26032021        
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VV_110_D")
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VH_110_D")
   enddo

!---------------------------------------------------------------------------------------------
!----------------A,B,C and D parameter tables for VV_139_D pol---------------------!
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VV_139_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest 
      call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VV_139_D_tbl_name, rc=rc)
      call LIS_verify(rc, "WCMRTM_multi AA_VV_139_D parameter table: not defined")
   enddo
   

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VV_139_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest   
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VV_139_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VV_139_D parameter table: not defined")
   enddo
    
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VV_139_D parameter table:",rc = rc)    
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VV_139_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VV_139_D parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VV_139_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VV_139_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VV_139_D parameter table: not defined")
   enddo



!----------------A,B,C and D parameter tables for VH_139_D pol---------------------!
   
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VH_139_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VH_139_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi AA_VH_139_D parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VH_139_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VH_139_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VH_139_D parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VH_139_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VH_139_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VH_139_D parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VH_139_D parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VH_139_D_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VH_139_D parameter table: not defined")
   enddo
!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VV_139_D pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VV_139_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VV_139_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VV_139_D(t),wcm_struc(n)%lone(t),&
          wcm_struc(n)%late(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VV_139_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VV_139_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VV_139_D_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VV_139_D(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VV_139_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VV_139_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VV_139_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VV_139_D(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VV_139_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VV_139_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VV_139_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VV_139_D(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo


!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VH_139_D pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VH_139_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VH_139_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VH_139_D(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VH_139_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VH_139_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VH_139_D_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VH_139_D(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VH_139_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VH_139_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VH_139_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VH_139_D(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VH_139_D pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VH_139_D_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VH_139_D_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VH_139_D(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo

   do n=1,LIS_rc%nnest !added fields to State 26032021        
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VV_139_D")
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VH_139_D")
   enddo

!---------------------------------------------------------------------------------------------
!----------------A,B,C and D parameter tables for VV_161_A pol---------------------!
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VV_161_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest 
      call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VV_161_A_tbl_name, rc=rc)
      call LIS_verify(rc, "WCMRTM_multi AA_VV_161_A parameter table: not defined")
   enddo
   

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VV_161_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest   
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VV_161_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VV_161_A parameter table: not defined")
   enddo
    
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VV_161_A parameter table:",rc = rc)    
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VV_161_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VV_161_A parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VV_161_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VV_161_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VV_161_A parameter table: not defined")
   enddo



!----------------A,B,C and D parameter tables for VH_161_A pol---------------------!
   
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi AA_VH_161_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%AA_VH_161_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi AA_VH_161_A parameter table: not defined")
   enddo

   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi BB_VH_161_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%BB_VH_161_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi BB_VH_161_A parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi CC_VH_161_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%CC_VH_161_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi CC_VH_161_A parameter table: not defined")
   enddo
   call ESMF_ConfigFindLabel(LIS_config, "WCMRTM_multi DD_VH_161_A parameter table:",rc = rc)
   do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,wcm_struc(n)%DD_VH_161_A_tbl_name, rc=rc)
       call LIS_verify(rc, "WCMRTM_multi DD_VH_161_A parameter table: not defined")
   enddo
!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VV_161_A pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VV_161_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VV_161_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VV_161_A(t),wcm_struc(n)%lone(t),&
          wcm_struc(n)%late(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VV_161_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VV_161_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VV_161_A_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VV_161_A(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VV_161_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VV_161_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VV_161_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VV_161_A(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VV_161_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VV_161_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VV_161_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VV_161_A(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo


!-------------------------------------------------------------------------
!!!!READ IN A PARAMETER TABLE for VH_161_A pol
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%AA_VH_161_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening AA_VH_161_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%AA_VH_161_A(t)
       enddo 
       CLOSE (19)
   enddo
!----------------------------------------------------------------------------
!-----READ IN B PARAMETER TABLE for VH_161_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%BB_VH_161_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening BB_VH_161_A_PARM.TBL'
         CALL wrf_error_fatal ( message ) 
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*)wcm_struc(n)%BB_VH_161_A(t)
       enddo

   ! READ (19,*) B_val
       CLOSE (19)
   enddo
!-----------------------------------------------------------------------------
!-----READ IN C PARAMETER TABLE for VH_161_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%CC_VH_161_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening CC_VH_161_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)  
           READ (19,*) wcm_struc(n)%CC_VH_161_A(t)
       enddo

!       READ (19,*) C_val
       CLOSE (19)
   enddo
!-------------------------------------------------------------------------------
!-----READ IN D PARAMETER TABLE for VH_161_A pol 
   do n=1, LIS_rc%nnest
       OPEN(19, FILE=trim(wcm_struc(n)%DD_VH_161_A_tbl_name),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr) 
       IF(ierr .NE. OPEN_OK ) THEN
         WRITE(message,FMT='(A)') &
         'failure opening DD_VH_161_A_PARM.TBL'
         CALL wrf_error_fatal ( message )
       END IF
       do t=1,LIS_rc%glbngrid(n)
           READ (19,*) wcm_struc(n)%DD_VH_161_A(t)
       enddo

!        READ (19,*) D_val
       CLOSE (19)
   enddo

   do n=1,LIS_rc%nnest !added fields to State 26032021        
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VV_161_A")
       call add_fields_toState(n,LIS_forwardState(n),"WCM_Sig0VH_161_A")
   enddo

   end subroutine WCMRTM_multi_initialize

   subroutine add_fields_toState(n, inState,varname) !added subroutine add-fields_toState 26032021

    use LIS_logMod,   only : LIS_verify
    use LIS_coreMod,  only : LIS_vecTile

    implicit none

    integer            :: n
    type(ESMF_State)   :: inState
    character(len=*)   :: varname

    type(ESMF_Field)     :: varField
    type(ESMF_ArraySpec) :: arrspec
    integer              :: status
    real :: sum
    call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
    call LIS_verify(status)

    varField = ESMF_FieldCreate(arrayspec=arrSpec, &
         grid=LIS_vecTile(n), name=trim(varname), &
         rc=status)
    call LIS_verify(status, 'Error in field_create of '//trim(varname))

    call ESMF_StateAdd(inState, (/varField/), rc=status)
    call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

   end subroutine add_fields_toState
!!--------------------------------------------------------------------------------

   subroutine add_sfc_fields(n, sfcState,varname)

   implicit none 

   integer            :: n 
   type(ESMF_State)   :: sfcState
   character(len=*)   :: varname

   type(ESMF_Field)     :: varField
   type(ESMF_ArraySpec) :: arrspec
   integer              :: status
   real :: sum
   call ESMF_ArraySpecSet(arrspec,rank=1,typekind=ESMF_TYPEKIND_R4,&
         rc=status)
   call LIS_verify(status)

   varField = ESMF_FieldCreate(arrayspec=arrSpec, & 
         grid=LIS_vecTile(n), name=trim(varname), &
         rc=status)
   call LIS_verify(status, 'Error in field_create of '//trim(varname))
    
   call ESMF_StateAdd(sfcState, (/varField/), rc=status)
   call LIS_verify(status, 'Error in StateAdd of '//trim(varname))

   end subroutine add_sfc_fields


   subroutine WCMRTM_multi_f2t(n)

   implicit none

   integer, intent(in)    :: n 

   end subroutine WCMRTM_multi_f2t


   subroutine WCMRTM_multi_geometry(n)
   implicit none
   integer, intent(in)    :: n

   end subroutine WCMRTM_multi_geometry 
  !Do nothing for now
   subroutine WCMRTM_multi_run(n)
   use LIS_histDataMod
! !USES: 
   implicit none

   integer, intent(in) :: n

   integer             :: t,p
   integer             :: status
   integer             :: col,row
   real                :: A_VV_037_D_cal,B_VV_037_D_cal,C_VV_037_D_cal,D_VV_037_D_cal,A_VH_037_D_cal,&
                          B_VH_037_D_cal, C_VH_037_D_cal, D_VH_037_D_cal,lon,lat,lon1,lat1
   real, pointer       :: sm(:), lai(:)
   real                :: sigmabare_VV_037_D,sigmabare_VH_037_D,s0VV_037_D_s_db, s0VH_037_D_s_dB, &
                        sigmacan_VV_037_D, sigmacan_VH_037_D,sigmasoil_VV_037_D,sigmasoil_VH_037_D,&
                        tt_VV_037_D, tt_VH_037_D
! ----------------------------- 088_A ----
   real                :: A_VV_088_A_cal,B_VV_088_A_cal,C_VV_088_A_cal,D_VV_088_A_cal,A_VH_088_A_cal,&
                          B_VH_088_A_cal, C_VH_088_A_cal, D_VH_088_A_cal
   real                :: sigmabare_VV_088_A,sigmabare_VH_088_A,s0VV_088_A_s_db, s0VH_088_A_s_dB, &
                        sigmacan_VV_088_A, sigmacan_VH_088_A,sigmasoil_VV_088_A,sigmasoil_VH_088_A,&
                        tt_VV_088_A, tt_VH_088_A
! ----------------------------- 110_D ----
   real                :: A_VV_110_D_cal,B_VV_110_D_cal,C_VV_110_D_cal,D_VV_110_D_cal,A_VH_110_D_cal,&
                          B_VH_110_D_cal, C_VH_110_D_cal, D_VH_110_D_cal
   real                :: sigmabare_VV_110_D,sigmabare_VH_110_D,s0VV_110_D_s_db, s0VH_110_D_s_dB, &
                        sigmacan_VV_110_D, sigmacan_VH_110_D,sigmasoil_VV_110_D,sigmasoil_VH_110_D,&
                        tt_VV_110_D, tt_VH_110_D
! ----------------------------- 139_D ----
   real                :: A_VV_139_D_cal,B_VV_139_D_cal,C_VV_139_D_cal,D_VV_139_D_cal,A_VH_139_D_cal,&
                          B_VH_139_D_cal, C_VH_139_D_cal, D_VH_139_D_cal
   real                :: sigmabare_VV_139_D,sigmabare_VH_139_D,s0VV_139_D_s_db, s0VH_139_D_s_dB, &
                        sigmacan_VV_139_D, sigmacan_VH_139_D,sigmasoil_VV_139_D,sigmasoil_VH_139_D,&
                        tt_VV_139_D, tt_VH_139_D
! ----------------------------- 161_A ----
   real                :: A_VV_161_A_cal,B_VV_161_A_cal,C_VV_161_A_cal,D_VV_161_A_cal,A_VH_161_A_cal,&
                          B_VH_161_A_cal, C_VH_161_A_cal, D_VH_161_A_cal
   real                :: sigmabare_VV_161_A,sigmabare_VH_161_A,s0VV_161_A_s_db, s0VH_161_A_s_dB, &
                        sigmacan_VV_161_A, sigmacan_VH_161_A,sigmasoil_VV_161_A,sigmasoil_VH_161_A,&
                        tt_VV_161_A, tt_VH_161_A

   real                :: theta, ctheta
   real, pointer       :: sig0val(:) !added for forward states 26032021

   ! theta = 0.6458 !incidence angle in radians (37 deg)
   theta = 0. !incidence angle in radians (i.e., rad(37deg for backscatter or rad(0deg) for gamma0)
   ctheta = cos(theta)


!   map surface properties to SFC    
   call getsfcvar(LIS_sfcState(n), "Soil Moisture Content",&
         sm)
   call getsfcvar(LIS_sfcState(n), "Leaf Area Index", &
         lai)

!---------------------------------------------
! Tile loop 
!--------------------------------------------
   do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
       col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
       lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
       lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
       do p=1,LIS_rc%glbngrid(n)
          lon1= wcm_struc(n)%lone(p)
          lat1= wcm_struc(n)%late(p)
          if (lon1 .eq. lon .and. lat1 .eq. lat) then
             A_VV_037_D_cal=wcm_struc(n)%AA_VV_037_D(p)
             B_VV_037_D_cal=wcm_struc(n)%BB_VV_037_D(p)
             C_VV_037_D_cal=wcm_struc(n)%CC_VV_037_D(p)
             D_VV_037_D_cal=wcm_struc(n)%DD_VV_037_D(p)        
             A_VH_037_D_cal=wcm_struc(n)%AA_VH_037_D(p)
             B_VH_037_D_cal=wcm_struc(n)%BB_VH_037_D(p)
             C_VH_037_D_cal=wcm_struc(n)%CC_VH_037_D(p)
             D_VH_037_D_cal=wcm_struc(n)%DD_VH_037_D(p)
          endif
       enddo
       
      if(.not.isNaN(sm(t)).and. sm(t).ne.LIS_rc%udef) then
         !bare soil backscatter in db
          s0VV_037_D_s_db=C_VV_037_D_cal+D_VV_037_D_cal*sm(t)
          s0VH_037_D_s_db=C_VH_037_D_cal+D_VH_037_D_cal*sm(t)
         !bare soil backscatter in linear units      
          sigmabare_VV_037_D=10.**(s0VV_037_D_s_db/10.)
          sigmabare_VH_037_D=10.**(s0VH_037_D_s_db/10.)
         !attenuation
          tt_VV_037_D=exp(-2.*B_VV_037_D_cal*lai(t)/ctheta)
          tt_VH_037_D=exp(-2.*B_VH_037_D_cal*lai(t)/ctheta)
         !attenuated soil backscatter
          sigmasoil_VV_037_D=tt_VV_037_D*sigmabare_VV_037_D
          sigmasoil_VH_037_D=tt_VH_037_D*sigmabare_VH_037_D
         !vegetation backscatter in linear units
          sigmacan_VV_037_D=(1.-tt_VV_037_D)*ctheta*(A_VV_037_D_cal*lai(t))
          sigmacan_VH_037_D=(1.-tt_VH_037_D)*ctheta*(A_VH_037_D_cal*lai(t))
         !total backscatter
          wcm_struc(n)%Sig0VV_037_D(t)=10.*log10(sigmacan_VV_037_D+sigmasoil_VV_037_D)
          wcm_struc(n)%Sig0VH_037_D(t)=10.*log10(sigmacan_VH_037_D+sigmasoil_VH_037_D)
       else
          wcm_struc(n)%Sig0VV_037_D(t)=LIS_rc%udef
          wcm_struc(n)%Sig0VH_037_D(t)=LIS_rc%udef

       endif

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VV,value=       &
          wcm_struc(n)%Sig0VV_037_D(t),             &
          vlevel=1, unit="dB",direction="-")

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VH,value=       &
          wcm_struc(n)%Sig0VH_037_D(t),             &
          vlevel=1, unit="dB",direction="-")
   enddo

   call getsfcvar(LIS_forwardState(n), "WCM_Sig0VV_037_D", sig0val) !added for forward states 26032021
   sig0val = wcm_struc(n)%Sig0VV_037_D

   call getsfcvar(LIS_forwardState(n),"WCM_Sig0VH_037_D", sig0val)
   sig0val = wcm_struc(n)%Sig0VH_037_D

! ----------------------------------- 088_A
   do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
       col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
       lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
       lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
       do p=1,LIS_rc%glbngrid(n)
          lon1= wcm_struc(n)%lone(p)
          lat1= wcm_struc(n)%late(p)
          if (lon1 .eq. lon .and. lat1 .eq. lat) then
             A_VV_088_A_cal=wcm_struc(n)%AA_VV_088_A(p)
             B_VV_088_A_cal=wcm_struc(n)%BB_VV_088_A(p)
             C_VV_088_A_cal=wcm_struc(n)%CC_VV_088_A(p)
             D_VV_088_A_cal=wcm_struc(n)%DD_VV_088_A(p)        
             A_VH_088_A_cal=wcm_struc(n)%AA_VH_088_A(p)
             B_VH_088_A_cal=wcm_struc(n)%BB_VH_088_A(p)
             C_VH_088_A_cal=wcm_struc(n)%CC_VH_088_A(p)
             D_VH_088_A_cal=wcm_struc(n)%DD_VH_088_A(p)
          endif
       enddo
       
      if(.not.isNaN(sm(t)).and. sm(t).ne.LIS_rc%udef) then
         !bare soil backscatter in db
          s0VV_088_A_s_db=C_VV_088_A_cal+D_VV_088_A_cal*sm(t)
          s0VH_088_A_s_db=C_VH_088_A_cal+D_VH_088_A_cal*sm(t)
         !bare soil backscatter in linear units      
          sigmabare_VV_088_A=10.**(s0VV_088_A_s_db/10.)
          sigmabare_VH_088_A=10.**(s0VH_088_A_s_db/10.)
         !attenuation
          tt_VV_088_A=exp(-2.*B_VV_088_A_cal*lai(t)/ctheta)
          tt_VH_088_A=exp(-2.*B_VH_088_A_cal*lai(t)/ctheta)
         !attenuated soil backscatter
          sigmasoil_VV_088_A=tt_VV_088_A*sigmabare_VV_088_A
          sigmasoil_VH_088_A=tt_VH_088_A*sigmabare_VH_088_A
         !vegetation backscatter in linear units
          sigmacan_VV_088_A=(1.-tt_VV_088_A)*ctheta*(A_VV_088_A_cal*lai(t))
          sigmacan_VH_088_A=(1.-tt_VH_088_A)*ctheta*(A_VH_088_A_cal*lai(t))
         !total backscatter
          wcm_struc(n)%Sig0VV_088_A(t)=10.*log10(sigmacan_VV_088_A+sigmasoil_VV_088_A)
          wcm_struc(n)%Sig0VH_088_A(t)=10.*log10(sigmacan_VH_088_A+sigmasoil_VH_088_A)
       else
          wcm_struc(n)%Sig0VV_088_A(t)=LIS_rc%udef
          wcm_struc(n)%Sig0VH_088_A(t)=LIS_rc%udef

       endif

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VV,value=       &
          wcm_struc(n)%Sig0VV_088_A(t),             &
          vlevel=1, unit="dB",direction="-")

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VH,value=       &
          wcm_struc(n)%Sig0VH_088_A(t),             &
          vlevel=1, unit="dB",direction="-")
   enddo

   call getsfcvar(LIS_forwardState(n), "WCM_Sig0VV_088_A", sig0val) !added for forward states 26032021
   sig0val = wcm_struc(n)%Sig0VV_088_A

   call getsfcvar(LIS_forwardState(n),"WCM_Sig0VH_088_A", sig0val)
   sig0val = wcm_struc(n)%Sig0VH_088_A

! ----------------------------------- 110_D
   do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
       col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
       lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
       lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
       do p=1,LIS_rc%glbngrid(n)
          lon1= wcm_struc(n)%lone(p)
          lat1= wcm_struc(n)%late(p)
          if (lon1 .eq. lon .and. lat1 .eq. lat) then
             A_VV_110_D_cal=wcm_struc(n)%AA_VV_110_D(p)
             B_VV_110_D_cal=wcm_struc(n)%BB_VV_110_D(p)
             C_VV_110_D_cal=wcm_struc(n)%CC_VV_110_D(p)
             D_VV_110_D_cal=wcm_struc(n)%DD_VV_110_D(p)        
             A_VH_110_D_cal=wcm_struc(n)%AA_VH_110_D(p)
             B_VH_110_D_cal=wcm_struc(n)%BB_VH_110_D(p)
             C_VH_110_D_cal=wcm_struc(n)%CC_VH_110_D(p)
             D_VH_110_D_cal=wcm_struc(n)%DD_VH_110_D(p)
          endif
       enddo
       
      if(.not.isNaN(sm(t)).and. sm(t).ne.LIS_rc%udef) then
         !bare soil backscatter in db
          s0VV_110_D_s_db=C_VV_110_D_cal+D_VV_110_D_cal*sm(t)
          s0VH_110_D_s_db=C_VH_110_D_cal+D_VH_110_D_cal*sm(t)
         !bare soil backscatter in linear units      
          sigmabare_VV_110_D=10.**(s0VV_110_D_s_db/10.)
          sigmabare_VH_110_D=10.**(s0VH_110_D_s_db/10.)
         !attenuation
          tt_VV_110_D=exp(-2.*B_VV_110_D_cal*lai(t)/ctheta)
          tt_VH_110_D=exp(-2.*B_VH_110_D_cal*lai(t)/ctheta)
         !attenuated soil backscatter
          sigmasoil_VV_110_D=tt_VV_110_D*sigmabare_VV_110_D
          sigmasoil_VH_110_D=tt_VH_110_D*sigmabare_VH_110_D
         !vegetation backscatter in linear units
          sigmacan_VV_110_D=(1.-tt_VV_110_D)*ctheta*(A_VV_110_D_cal*lai(t))
          sigmacan_VH_110_D=(1.-tt_VH_110_D)*ctheta*(A_VH_110_D_cal*lai(t))
         !total backscatter
          wcm_struc(n)%Sig0VV_110_D(t)=10.*log10(sigmacan_VV_110_D+sigmasoil_VV_110_D)
          wcm_struc(n)%Sig0VH_110_D(t)=10.*log10(sigmacan_VH_110_D+sigmasoil_VH_110_D)
       else
          wcm_struc(n)%Sig0VV_110_D(t)=LIS_rc%udef
          wcm_struc(n)%Sig0VH_110_D(t)=LIS_rc%udef

       endif

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VV,value=       &
          wcm_struc(n)%Sig0VV_110_D(t),             &
          vlevel=1, unit="dB",direction="-")

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VH,value=       &
          wcm_struc(n)%Sig0VH_110_D(t),             &
          vlevel=1, unit="dB",direction="-")
   enddo

   call getsfcvar(LIS_forwardState(n), "WCM_Sig0VV_110_D", sig0val) !added for forward states 26032021
   sig0val = wcm_struc(n)%Sig0VV_110_D

   call getsfcvar(LIS_forwardState(n),"WCM_Sig0VH_110_D", sig0val)
   sig0val = wcm_struc(n)%Sig0VH_110_D

! ----------------------------------- 139_D
   do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
       col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
       lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
       lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
       do p=1,LIS_rc%glbngrid(n)
          lon1= wcm_struc(n)%lone(p)
          lat1= wcm_struc(n)%late(p)
          if (lon1 .eq. lon .and. lat1 .eq. lat) then
             A_VV_139_D_cal=wcm_struc(n)%AA_VV_139_D(p)
             B_VV_139_D_cal=wcm_struc(n)%BB_VV_139_D(p)
             C_VV_139_D_cal=wcm_struc(n)%CC_VV_139_D(p)
             D_VV_139_D_cal=wcm_struc(n)%DD_VV_139_D(p)        
             A_VH_139_D_cal=wcm_struc(n)%AA_VH_139_D(p)
             B_VH_139_D_cal=wcm_struc(n)%BB_VH_139_D(p)
             C_VH_139_D_cal=wcm_struc(n)%CC_VH_139_D(p)
             D_VH_139_D_cal=wcm_struc(n)%DD_VH_139_D(p)
          endif
       enddo
       
      if(.not.isNaN(sm(t)).and. sm(t).ne.LIS_rc%udef) then
         !bare soil backscatter in db
          s0VV_139_D_s_db=C_VV_139_D_cal+D_VV_139_D_cal*sm(t)
          s0VH_139_D_s_db=C_VH_139_D_cal+D_VH_139_D_cal*sm(t)
         !bare soil backscatter in linear units      
          sigmabare_VV_139_D=10.**(s0VV_139_D_s_db/10.)
          sigmabare_VH_139_D=10.**(s0VH_139_D_s_db/10.)
         !attenuation
          tt_VV_139_D=exp(-2.*B_VV_139_D_cal*lai(t)/ctheta)
          tt_VH_139_D=exp(-2.*B_VH_139_D_cal*lai(t)/ctheta)
         !attenuated soil backscatter
          sigmasoil_VV_139_D=tt_VV_139_D*sigmabare_VV_139_D
          sigmasoil_VH_139_D=tt_VH_139_D*sigmabare_VH_139_D
         !vegetation backscatter in linear units
          sigmacan_VV_139_D=(1.-tt_VV_139_D)*ctheta*(A_VV_139_D_cal*lai(t))
          sigmacan_VH_139_D=(1.-tt_VH_139_D)*ctheta*(A_VH_139_D_cal*lai(t))
         !total backscatter
          wcm_struc(n)%Sig0VV_139_D(t)=10.*log10(sigmacan_VV_139_D+sigmasoil_VV_139_D)
          wcm_struc(n)%Sig0VH_139_D(t)=10.*log10(sigmacan_VH_139_D+sigmasoil_VH_139_D)
       else
          wcm_struc(n)%Sig0VV_139_D(t)=LIS_rc%udef
          wcm_struc(n)%Sig0VH_139_D(t)=LIS_rc%udef

       endif

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VV,value=       &
          wcm_struc(n)%Sig0VV_139_D(t),             &
          vlevel=1, unit="dB",direction="-")

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VH,value=       &
          wcm_struc(n)%Sig0VH_139_D(t),             &
          vlevel=1, unit="dB",direction="-")
   enddo

   call getsfcvar(LIS_forwardState(n), "WCM_Sig0VV_139_D", sig0val) !added for forward states 26032021
   sig0val = wcm_struc(n)%Sig0VV_139_D

   call getsfcvar(LIS_forwardState(n),"WCM_Sig0VH_139_D", sig0val)
   sig0val = wcm_struc(n)%Sig0VH_139_D

! ----------------------------------- 161_A
   do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
       row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
       col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
       lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
       lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
       do p=1,LIS_rc%glbngrid(n)
          lon1= wcm_struc(n)%lone(p)
          lat1= wcm_struc(n)%late(p)
          if (lon1 .eq. lon .and. lat1 .eq. lat) then
             A_VV_161_A_cal=wcm_struc(n)%AA_VV_161_A(p)
             B_VV_161_A_cal=wcm_struc(n)%BB_VV_161_A(p)
             C_VV_161_A_cal=wcm_struc(n)%CC_VV_161_A(p)
             D_VV_161_A_cal=wcm_struc(n)%DD_VV_161_A(p)        
             A_VH_161_A_cal=wcm_struc(n)%AA_VH_161_A(p)
             B_VH_161_A_cal=wcm_struc(n)%BB_VH_161_A(p)
             C_VH_161_A_cal=wcm_struc(n)%CC_VH_161_A(p)
             D_VH_161_A_cal=wcm_struc(n)%DD_VH_161_A(p)
          endif
       enddo
       
      if(.not.isNaN(sm(t)).and. sm(t).ne.LIS_rc%udef) then
         !bare soil backscatter in db
          s0VV_161_A_s_db=C_VV_161_A_cal+D_VV_161_A_cal*sm(t)
          s0VH_161_A_s_db=C_VH_161_A_cal+D_VH_161_A_cal*sm(t)
         !bare soil backscatter in linear units      
          sigmabare_VV_161_A=10.**(s0VV_161_A_s_db/10.)
          sigmabare_VH_161_A=10.**(s0VH_161_A_s_db/10.)
         !attenuation
          tt_VV_161_A=exp(-2.*B_VV_161_A_cal*lai(t)/ctheta)
          tt_VH_161_A=exp(-2.*B_VH_161_A_cal*lai(t)/ctheta)
         !attenuated soil backscatter
          sigmasoil_VV_161_A=tt_VV_161_A*sigmabare_VV_161_A
          sigmasoil_VH_161_A=tt_VH_161_A*sigmabare_VH_161_A
         !vegetation backscatter in linear units
          sigmacan_VV_161_A=(1.-tt_VV_161_A)*ctheta*(A_VV_161_A_cal*lai(t))
          sigmacan_VH_161_A=(1.-tt_VH_161_A)*ctheta*(A_VH_161_A_cal*lai(t))
         !total backscatter
          wcm_struc(n)%Sig0VV_161_A(t)=10.*log10(sigmacan_VV_161_A+sigmasoil_VV_161_A)
          wcm_struc(n)%Sig0VH_161_A(t)=10.*log10(sigmacan_VH_161_A+sigmasoil_VH_161_A)
       else
          wcm_struc(n)%Sig0VV_161_A(t)=LIS_rc%udef
          wcm_struc(n)%Sig0VH_161_A(t)=LIS_rc%udef

       endif

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VV,value=       &
          wcm_struc(n)%Sig0VV_161_A(t),             &
          vlevel=1, unit="dB",direction="-")

       call LIS_diagnoseRTMOutputVar(n, t,LIS_MOC_RTM_Sig0VH,value=       &
          wcm_struc(n)%Sig0VH_161_A(t),             &
          vlevel=1, unit="dB",direction="-")
   enddo

   call getsfcvar(LIS_forwardState(n), "WCM_Sig0VV_161_A", sig0val) !added for forward states 26032021
   sig0val = wcm_struc(n)%Sig0VV_161_A

   call getsfcvar(LIS_forwardState(n),"WCM_Sig0VH_161_A", sig0val)
   sig0val = wcm_struc(n)%Sig0VH_161_A

   end subroutine WCMRTM_multi_run


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   subroutine getsfcvar(sfcState, varname, var)
! !USES: 
    
   implicit none
    
   type(ESMF_State)      :: sfcState
   type(ESMF_Field)      :: varField
   character(len=*)      :: varname
   real, pointer         :: var(:)
   integer               :: status

   call ESMF_StateGet(sfcState, trim(varname), varField, rc=status)
   call LIS_verify(status, 'Error in StateGet: CMEM3_handlerMod '//trim(varname))
   call ESMF_FieldGet(varField, localDE=0,farrayPtr=var, rc=status)
   call LIS_verify(status, 'Error in FieldGet: CMEM3_handlerMod '//trim(varname))

   end subroutine getsfcvar

!!!!BOP
!!!! !ROUTINE: WCMRTM_multi_output
!!!! \label{WCMRTM_multi_output}
!!!!
!!!! !INTERFACE: 
   subroutine WCMRTM_multi_output(n)
   integer, intent(in) :: n 
   end subroutine WCMRTM_multi_output
#endif
end module WCMRTM_multi_Mod



