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
! !ROUTINE: ac70_updatesoilmLAI
!  \label{ac70_updatesoilmLAI}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 9 Sep 2016: Mahdi Navari; Modified for ac70 
!   To do: makes it general for x layers (currently hard coded for 4 layers)
! 18 Jun 2021: Michel Bechtold: SM and LAI updating with S1 backscatter w/ WCM
!
! !INTERFACE:
subroutine ac70_updatesoilmLAI(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use ac70_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP

  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm1IncrField
  real, pointer          :: soilm1(:)
  real, pointer          :: soilmIncr1(:)

  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm2IncrField
  real, pointer          :: soilm2(:)
  real, pointer          :: soilmIncr2(:)

  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm3IncrField
  real, pointer          :: soilm3(:)
  real, pointer          :: soilmIncr3(:)

  type(ESMF_Field)       :: AC70BIOMASSField
  type(ESMF_Field)       :: AC70BIOMASSIncrField
  real, pointer          :: AC70BIOMASS(:)
  real, pointer          :: AC70BIOMASSincr(:)
  integer                :: t,i,m,gid
  integer                :: status

  real                   :: AC70BIOMASStmp,AC70BIOMASSmax,AC70BIOMASSmin

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: AC70BIOMASSmean(LIS_rc%ngrid(n))
  integer                :: nAC70BIOMASSmean(LIS_rc%ngrid(n))

  ! Layer 1
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in ac70_updatesoilmLAI")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in ac70_updatesoilmLAI")

  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 1",sm1IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in ac70_updatesoilmLAI")

  call ESMF_FieldGet(sm1IncrField,localDE=0,farrayPtr=soilmIncr1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in ac70_updatesoilmLAI")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = soilm1(t) + soilmIncr1(t)
  enddo

  ! Layer 2
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 2 failed in ac70_updatesoilmLAI")

  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in ac70_updatesoilmLAI")

  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 2",sm2IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 2 failed in ac70_updatesoilmLAI")

  call ESMF_FieldGet(sm2IncrField,localDE=0,farrayPtr=soilmIncr2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in ac70_updatesoilmLAI")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm2(t) = soilm2(t) + soilmIncr2(t)
  enddo

  ! Layer 3
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 3 failed in ac70_updatesoilmLAI")

  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in ac70_updatesoilmLAI")

  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 3",sm3IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 3 failed in ac70_updatesoilmLAI")

  call ESMF_FieldGet(sm3IncrField,localDE=0,farrayPtr=soilmIncr3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in ac70_updatesoilmLAI")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm3(t) = soilm3(t) + soilmIncr3(t)
  enddo

  ! BIOMASS
  call ESMF_StateGet(LSM_State,"AC70 BIOMASS",AC70BIOMASSField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: LSM_State, failed in ac70_updatesoilmLAI")
  call ESMF_FieldGet(AC70BIOMASSField,localDE=0,farrayPtr=AC70BIOMASS,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: AC70BIOMASSField failed in ac70_updatesoilmLAI")
 
  call ESMF_StateGet(LSM_Incr_State,"AC70 BIOMASS",AC70BIOMASSIncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: LSM_Incr_State AC70BIOMASS failed in ac70_updatesoilmLAI")
  call ESMF_FieldGet(AC70BIOMASSIncrField,localDE=0,farrayPtr=AC70BIOMASSincr,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: AC70BIOMASSIncrField failed in ac70_updatesoilmLAI")

  call ESMF_AttributeGet(AC70BIOMASSField,"Max Value",AC70BIOMASSmax,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: AC70BIOMASSField Max Value failed in ac70_updatesoilmLAI")
  call ESMF_AttributeGet(AC70BIOMASSField,"Min Value",AC70BIOMASSmin,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: AC70BIOMASSField Min Value failed in ac70_updatesoilmLAI")


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = soilm1(t) + soilmIncr1(t)
  enddo


  update_flag    = .true.
  perc_violation = 0.0
  AC70BIOMASSmean       = 0.0
  nAC70BIOMASSmean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     AC70BIOMASStmp =  AC70BIOMASS(t) + AC70BIOMASSincr(t)


     if(AC70BIOMASStmp.lt.AC70BIOMASSmin.or.AC70BIOMASStmp.gt.AC70BIOMASSmax) then
        update_flag(gid) = .false.
        perc_violation(gid) = perc_violation(gid) +1
     endif

  enddo

  do gid=1,LIS_rc%ngrid(n)
     perc_violation(gid) = perc_violation(gid)/LIS_rc%nensem(n)
  enddo

! For ensembles that are unphysical, compute the
! ensemble average after excluding them. This
! is done only if the majority of the ensemble
! members are good (>60%)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)
     if(.not.update_flag(gid)) then
        if(perc_violation(gid).lt.0.8) then
           if((AC70BIOMASS(t)+AC70BIOMASSincr(t).gt.AC70BIOMASSmin).and.&
                (AC70BIOMASS(t)+AC70BIOMASSincr(t).lt.AC70BIOMASSmax)) then 
              AC70BIOMASSmean(gid) = AC70BIOMASSmean(gid) + &
                   AC70BIOMASS(t) + AC70BIOMASSincr(t)
              nAC70BIOMASSmean(gid) = nAC70BIOMASSmean(gid) + 1
           endif
        endif
     endif
  enddo

 do gid=1,LIS_rc%ngrid(n)
     if(nAC70BIOMASSmean(gid).gt.0) then
        AC70BIOMASSmean(gid) = AC70BIOMASSmean(gid)/nAC70BIOMASSmean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     AC70BIOMASStmp =  AC70BIOMASS(t) + AC70BIOMASSincr(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        AC70BIOMASS(t) = AC70BIOMASStmp
     elseif(perc_violation(gid).lt.0.8) then
        if(AC70BIOMASStmp.lt.AC70BIOMASSmin.or.AC70BIOMASStmp.gt.AC70BIOMASSmax) then
           AC70BIOMASS(t) = AC70BIOMASSmean(gid)
        else
           AC70BIOMASS(t) = AC70BIOMASS(t) + AC70BIOMASSincr(t)
        endif
     endif
  enddo

end subroutine ac70_updatesoilmLAI

