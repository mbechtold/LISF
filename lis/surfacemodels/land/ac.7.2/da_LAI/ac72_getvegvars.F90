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
! !ROUTINE: ac72_getvegvars
! \label{ac72_getvegvars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac72_getvegvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use ac72_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the related state prognostic variables for
!  veg data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)       :: AC72BIOMASSField

  integer                :: t
  integer                :: status
  real, pointer          :: AC72BIOMASS(:)
 
  call ESMF_StateGet(LSM_State,"AC72 BIOMASS",AC72BIOMASSField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for AC72BIOMASS in ac72_getsoilmLAI')
  call ESMF_FieldGet(AC72BIOMASSField,localDE=0,farrayPtr=AC72BIOMASS,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for AC72BIOMASS in ac72_getsoilmLAI')


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     AC72BIOMASS(t) = AC72_struc(n)%ac72(t)%SumWaBal%Biomass
  enddo

end subroutine ac72_getvegvars

