!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: Ac70_sfc2wcm
!  \label{Ac70_sfc2wcm}
!
! !REVISION HISTORY:
!  4 Sep 2020: Sara Modanesi; Initial Specification
!  23 Nov 2022: Michel Bechtold; Soil moisture aggregation needed for Aquacrop physics
! !INTERFACE:
subroutine ac70_sfc2wcm(n, sfcState)
! !USES:      
  use ESMF
  use LIS_coreMod
  use LIS_logMod,    only : LIS_verify
  use LIS_constantsMod,  only : LIS_CONST_RHOFW

  use Ac70_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  type(ESMF_State)    :: sfcState

! FUNCTIONS

! 
! !DESCRIPTION: 
! This subroutine assigns the ac70 specific surface variables
! to WCM. 
!
!EOP
  type(ESMF_Field)    :: smcField, laiField
  real, pointer       :: soil_moisture_content(:), leaf_area_index(:)
  real                :: w1, w2, w3
  integer             :: t,status

  call ESMF_StateGet(sfcState,"Soil Moisture Content",smcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=soil_moisture_content,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Leaf Area Index",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=leaf_area_index, rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     w2 = 0.15
     w3 = 0.1
     w1 = 1 - w2 - w3
     soil_moisture_content(t) = w1*AC70_struc(n)%ac70(t)%smc(1) + &
                                w2*AC70_struc(n)%ac70(t)%smc(2) + &
                                w3*AC70_struc(n)%ac70(t)%smc(3)
     Leaf_Area_Index(t) = AC70_struc(n)%ac70(t)%WCMV1V2
  enddo

end subroutine ac70_sfc2wcm
