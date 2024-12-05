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
!BOP
! !ROUTINE: ac72_sfc2cmem3
! \label{ac72_sfc2cmem3}
!
! !REVISION HISTORY:
!  37 Mar 2009: Sujay Kumar; Initial Code
!  2  May 2011: Yudong Tian; Modified for Noah 3.2 
!  22 Jun 2018: Peter Shellito; Modified for ac72. 
!
! !INTERFACE:
subroutine ac72_sfc2cmem3(n, sfcState)
! !USES:      
  use ESMF
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_verify
  use LIS_historyMod, only : LIS_patch2tile
  use LIS_constantsMod,  only : LIS_CONST_RHOFW

  use ac72_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  type(ESMF_State)    :: sfcState
  
! FUNCTIONS

! 
! !DESCRIPTION: 
! This subroutine assigns the ac72 specific surface variables
! to CMEM3.
!
!EOP

! Initialize variables
  integer             :: t, status
! These vars are of type ESMF_Field
! The next 3 lines contain the fields but without relsmcField
!  type(ESMF_Field)    :: windField, lcField, ltField, ltempField, stempField, &
!       smcField, stcField, cmcField, vegFracField, &
!       scField, snowhField, snodField, laiField
  type(ESMF_Field)    :: windField, lcField, ltField, ltempField, &
       smcField, relsmcField, cmcField, vegFracField, &
       scField, laiField
! These vars are pointers
  real, pointer       :: wind_speed(:), land_coverage(:), land_temperature(:), &
       soil_moisture_content(:), soil_temperature(:), canopy_water_content(:), &
       vegetation_fraction(:), land_type(:), &
       leaf_area_index(:), &
       relative_soil_moisture_content(:)
! the line above doesn't have relative_soil_moisture_content
  real, allocatable       :: lsm_var(:)
  real, allocatable       :: land_cov(:)

! And here we assign/determine the actual bounds for the array
  allocate(lsm_var(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(land_cov(LIS_rc%npatch(n,LIS_rc%lsm_index)))

! I don't have to do anything with these lines. 
  call ESMF_StateGet(sfcState,"Wind Speed",windField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(windField,localDE=0,farrayPtr=wind_speed,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Land Coverage",lcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(lcField,localDE=0,farrayPtr=land_coverage,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Land Type",ltField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(ltField,localDE=0,farrayPtr=land_type,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Land Temperature",ltempField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(ltempField,localDE=0,farrayPtr=land_temperature,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Soil Moisture Content",smcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=soil_moisture_content,&
       rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Canopy Water Content",cmcField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(cmcField,localDE=0,farrayPtr=canopy_water_content,&
       rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Vegetation Fraction",vegfracField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(vegfracField,localDE=0,farrayPtr=vegetation_fraction,&
       rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Leaf Area Index",laiField,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=leaf_area_index,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(sfcState,"Relative Soil Moisture Content",relsmcField,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(relsmcField,localDE=0,&
       farrayPtr=relative_soil_moisture_content,rc=status)
  call LIS_verify(status)

! In this section, we loop through the "patches" to accumulate the variable for each patch, then the routine LIS_patch2tile accumulates that patch over each tile and saves it as the variable named in the second-to-last parameter there.
! The LIS_patch2tile is found in LIS_historyMod.F90 in /core
! public :: LIS_patch2tile ! convert a patchspace variable to tilespace
! subroutine LIS_patch2tile(n,m,tvar,pvar):
! n is the index of the domain or nest
! m is unknown, probably another index
! tvar is the variable dimensioned in the tile space
! pvar is the variable dimensioned in the patch space

  !REWRITE THIS ROUTINE AS AN ESMF_STATE + FIELD STYLE. 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ! obviously much work needs to be done here. 
     !gross type of surface determined by coverage
     ! assume only land points
     ! Surface type independent data
     lsm_var(t) = sqrt(AC72_struc(n)%ac72(t)%wind_e*& 
          AC72_struc(n)%ac72(t)%wind_e+& 
          AC72_struc(n)%ac72(t)%wind_n*&
          AC72_struc(n)%ac72(t)%wind_n)
  enddo
  !put in wind direction using math on u and v; deg. E. of N.
  call LIS_patch2tile(n,LIS_rc%lsm_index,wind_speed,lsm_var)
     ! Snow surface type data
     ! Snow depth are m in CMEM3 and in noah
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     land_cov(t) = 1.0
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,land_coverage,land_cov)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t)     = AC72_struc(n)%ac72(t)%vegetype 
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,land_type, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = AC72_struc(n)%ac72(t)%trad
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,land_temperature, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = AC72_struc(n)%ac72(t)%trad
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,soil_moisture_content, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = AC72_struc(n)%ac72(t)%fveg
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,vegetation_fraction, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ! --------------------------------------------------------
     ! Canopy water content (mm): only includes canopy interception (cmc). 
     ! VWC will be calculated in CMEM module with lai  
     ! --------------------------------------------------------
     lsm_var(t) = AC72_struc(n)%ac72(t)%canliq
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,canopy_water_content, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     lsm_var(t) = AC72_struc(n)%ac72(t)%lai
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,leaf_area_index, lsm_var)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           lsm_var(t) = &  
                (AC72_struc(n)%ac72(t)%smc(1) - &
                AC72_struc(n)%ac72(t)%smcwlt) / &
                (AC72_struc(n)%ac72(t)%smcmax - &
                AC72_struc(n)%ac72(t)%smcwlt)
           if ( lsm_var(t) > 1.0 ) then
              lsm_var(t)  = 1.0
           endif
           if ( lsm_var(t)  < 0.01 ) then
              lsm_var(t)  = 0.01
           endif
  enddo
  call LIS_patch2tile(n,LIS_rc%lsm_index,relative_soil_moisture_content, &
       lsm_var)

  deallocate(lsm_var)
  deallocate(land_cov)

end subroutine ac72_sfc2cmem3
