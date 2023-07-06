!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: ac70_setvegvars
! \label{ac70_setvegvars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac70_setvegvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use ac70_lsmMod
  use AC_VEG_PARAMETERS_70
  use LIS_logMod, only : LIS_logunit, LIS_verify
  use ac_kinds, only: dp

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
! 
! !DESCRIPTION:
! 
!  This routine assigns the progognostic variables to noah's
!  model space. The state vector consists of veg
! 
!EOP

  type(ESMF_Field)       :: AC70BIOMASSField

  integer                :: t
  integer                :: status
  real, pointer          :: AC70BIOMASS(:)
 
  ! BIOMASS
  call ESMF_StateGet(LSM_State,"AC70 BIOMASS",AC70BIOMASSField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: AC70BIOMASS failed in ac70_setsoilmLAI")
  call ESMF_FieldGet(AC70BIOMASSField,localDE=0,farrayPtr=AC70BIOMASS,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: AC70BIOMASS failed in ac70_setsoilmLAI")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     AC70_struc(n)%ac70(t)%SumWaBal%Biomass = AC70BIOMASS(t)
  enddo
  
end subroutine ac70_setvegvars


