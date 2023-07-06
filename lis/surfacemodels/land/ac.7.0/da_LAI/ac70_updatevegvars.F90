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
! !ROUTINE: ac70_updatevegvars
! \label{ac70_updatevegvars}
!
! !REVISION HISTORY:
! 22 Dec 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine ac70_updatevegvars(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use ac70_lsmMod
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

  type(ESMF_Field)       :: AC70BIOMASSField
  type(ESMF_Field)       :: AC70BIOMASSIncrField

  integer                :: t,gid
  integer                :: status
  real, pointer          :: AC70BIOMASS(:)
  real, pointer          :: AC70BIOMASSincr(:)
  real                   :: AC70BIOMASStmp,AC70BIOMASSmax,AC70BIOMASSmin

  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: AC70BIOMASSmean(LIS_rc%ngrid(n))
  integer                :: nAC70BIOMASSmean(LIS_rc%ngrid(n))


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


end subroutine ac70_updatevegvars

