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
! !ROUTINE: ac72_updatevegvars
! \label{ac72_updatevegvars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac72_updatevegvars(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use ac72_lsmMod
  use LIS_logMod,   only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Returns the snow related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \item[LSM\_Incr\_State] ESMF State container for LSM state increments \newline
!  \end{description}
!
!EOP

  type(ESMF_Field)       :: AC72BIOMASSField
  type(ESMF_Field)       :: AC72BIOMASSIncrField

  integer                :: t,gid
  integer                :: status
  real, pointer          :: AC72BIOMASS(:)
  real, pointer          :: AC72BIOMASSincr(:)
  real                   :: AC72BIOMASStmp,AC72BIOMASSmax,AC72BIOMASSmin

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: AC72BIOMASSmean(LIS_rc%ngrid(n))
  integer                :: nAC72BIOMASSmean(LIS_rc%ngrid(n))


  ! BIOMASS
  call ESMF_StateGet(LSM_State,"AC72 BIOMASS",AC72BIOMASSField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: LSM_State, failed in ac72_updatesoilmLAI")
  call ESMF_FieldGet(AC72BIOMASSField,localDE=0,farrayPtr=AC72BIOMASS,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: AC72BIOMASSField failed in ac72_updatesoilmLAI")
 
  call ESMF_StateGet(LSM_Incr_State,"AC72 BIOMASS",AC72BIOMASSIncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: LSM_Incr_State AC72BIOMASS failed in ac72_updatesoilmLAI")
  call ESMF_FieldGet(AC72BIOMASSIncrField,localDE=0,farrayPtr=AC72BIOMASSincr,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: AC72BIOMASSIncrField failed in ac72_updatesoilmLAI")

  call ESMF_AttributeGet(AC72BIOMASSField,"Max Value",AC72BIOMASSmax,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: AC72BIOMASSField Max Value failed in ac72_updatesoilmLAI")
  call ESMF_AttributeGet(AC72BIOMASSField,"Min Value",AC72BIOMASSmin,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: AC72BIOMASSField Min Value failed in ac72_updatesoilmLAI")



  update_flag    = .true.
  perc_violation = 0.0
  AC72BIOMASSmean       = 0.0
  nAC72BIOMASSmean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     AC72BIOMASStmp =  AC72BIOMASS(t) + AC72BIOMASSincr(t)


     if(AC72BIOMASStmp.lt.AC72BIOMASSmin.or.AC72BIOMASStmp.gt.AC72BIOMASSmax) then
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
           if((AC72BIOMASS(t)+AC72BIOMASSincr(t).gt.AC72BIOMASSmin).and.&
                (AC72BIOMASS(t)+AC72BIOMASSincr(t).lt.AC72BIOMASSmax)) then 
              AC72BIOMASSmean(gid) = AC72BIOMASSmean(gid) + &
                   AC72BIOMASS(t) + AC72BIOMASSincr(t)
              nAC72BIOMASSmean(gid) = nAC72BIOMASSmean(gid) + 1
           endif
        endif
     endif
  enddo

 do gid=1,LIS_rc%ngrid(n)
     if(nAC72BIOMASSmean(gid).gt.0) then
        AC72BIOMASSmean(gid) = AC72BIOMASSmean(gid)/nAC72BIOMASSmean(gid)
     endif
  enddo


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     AC72BIOMASStmp =  AC72BIOMASS(t) + AC72BIOMASSincr(t)

! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not
! update.

     if(update_flag(gid)) then
        AC72BIOMASS(t) = AC72BIOMASStmp
     elseif(perc_violation(gid).lt.0.8) then
        if(AC72BIOMASStmp.lt.AC72BIOMASSmin.or.AC72BIOMASStmp.gt.AC72BIOMASSmax) then
           AC72BIOMASS(t) = AC72BIOMASSmean(gid)
        else
           AC72BIOMASS(t) = AC72BIOMASS(t) + AC72BIOMASSincr(t)
        endif
     endif
  enddo


end subroutine ac72_updatevegvars

