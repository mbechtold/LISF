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
! !ROUTINE: ac72_getsoilmLAI
! \label{ac72_getsoilmLAI}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 1 Aug 2016: Mahdi Navari; Modified for ac72 
!   To do: makes it general for x layers (currently hard coded for 4 layers)
! 18 Jun 2021: Michel Bechtold: SM and Biomass updating with S1 backscatter w/ WCM
! !INTERFACE:
subroutine ac72_getsoilmLAI(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use ac72_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture and Biomass related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sm1Field
  real, pointer          :: soilm1(:)
  type(ESMF_Field)       :: sm2Field
  real, pointer          :: soilm2(:)
  type(ESMF_Field)       :: sm3Field
  real, pointer          :: soilm3(:)
  type(ESMF_Field)       :: AC72BIOMASSField
  real, pointer          :: AC72BIOMASS(:)
  integer                :: t
  integer                :: status
  character*100          :: lsm_state_objs(2)


  !Layer 1
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in ac72_getsoilmLAI')

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm1 in ac72_getsoilmLAI')

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = AC72_struc(n)%ac72(t)%smc(1)
  enddo
  
  !Layer 2
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm2 in ac72_getsoilmLAI')

  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm2 in ac72_getsoilmLAI')

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm2(t) = AC72_struc(n)%ac72(t)%smc(2)
  enddo

  !Layer 3
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm3 in ac72_getsoilmLAI')

  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm3 in ac72_getsoilmLAI')

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm3(t) = AC72_struc(n)%ac72(t)%smc(3)
  enddo

  call ESMF_StateGet(LSM_State,"AC72 BIOMASS",AC72BIOMASSField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for AC72BIOMASS in ac72_getsoilmLAI')
  call ESMF_FieldGet(AC72BIOMASSField,localDE=0,farrayPtr=AC72BIOMASS,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for AC72BIOMASS in ac72_getsoilmLAI')


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     AC72BIOMASS(t) = AC72_struc(n)%ac72(t)%SumWaBal%Biomass
  enddo

end subroutine ac72_getsoilmLAI

