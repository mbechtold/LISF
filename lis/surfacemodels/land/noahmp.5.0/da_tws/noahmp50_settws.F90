!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp50_settws
!  \label{noahmp50_settws}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!
! !INTERFACE:
subroutine noahmp50_settws(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_constantsMod
  use LIS_logMod
  use NoahMP50_lsmMod
  use MicroTopoCorrectionMod, only : CalcWaterTableFromSoil_m, CalcEquilibriumProfile_m

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture and groundwater prognostic variables
!  to NoahMP's model space.
! 
!EOP
  integer, parameter     :: PEAT_SOILTYPE = 17
  real, parameter        :: MIN_GWS_THRESHOLD = 0.00
  real, parameter        :: MAX_GWS_THRESHOLD = 7000.0
  real, parameter        :: MAX_WA = 7000.0
  real, parameter        :: ZSOIL = 2 !mm
  real, parameter        :: ROUS = 0.2 ! specific yield
  !Bailing changed this to be WLTSMC
!  real, parameter        :: MIN_THRESHOLD = 0.02 
  real                   :: MIN_THRESHOLD 
  real                   :: MAX_THRESHOLD
  real                   :: sm_threshold
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  type(ESMF_Field)       :: gwField
  type(ESMF_Field)       :: sweField

  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: gws(:)
  real, pointer          :: swe(:)

  real                   :: delta
  logical                :: diffCheck(LIS_rc%ngrid(n))
  logical                :: ensCheck(LIS_rc%ngrid(n))
  logical                :: largeSM(LIS_rc%ngrid(n))
  real                   :: snodens(LIS_rc%npatch(n,LIS_rc%lsm_index))
  integer                :: i, c,r,t,m,gid,j
  integer                :: SOILTYP           ! soil type index [-]
  real                   :: W_soil, z_col_bot, wtd_bg, wtd_an
  real                   :: thetas_l, he_l, b_l
  real                   :: bgWsoil(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real                   :: ice
  real                   :: sm_eq_arr(NoahMP50_struc(n)%nsoil)
  real, parameter        :: WTD_EQUIL_THRESHOLD = 0.5  ! [m] Richards/equil split (peat physics)
  real                   :: dsneqv,dsnowh,snowh_new
  real                   :: TWS1, TWS2, TWSd,delta1
  real                   :: remainder
  integer                :: status
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: rc1,rc2,rc3,rc4,rc5
  
  
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in noahmp50_settws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 2 failed in noahmp50_settws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 3 failed in noahmp50_settws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 4 failed in noahmp50_settws")
  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp50_settws")
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: SWE failed in noahmp50_settws")


  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp50_settws")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in noahmp50_settws")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in noahmp50_settws")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in noahmp50_settws")
  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Groundwater Storage failed in noahmp50_settws")
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: SWE failed in noahmp50_settws")


  ensCheck = .true.
  diffCheck = .false.
  largeSM  = .false.

  ! Peatland: snapshot background column soil water [m] (pre-analysis) so the
  ! water-table update at the end carries only the applied DA increment.
  if(NoahMP50_struc(n)%peat_opt.eq.1) then
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        bgWsoil(t) = 0.0
        if(NoahMP50_struc(n)%noahmp50(t)%soiltype.eq.PEAT_SOILTYPE) then
           do j=1,NoahMP50_struc(n)%nsoil
              bgWsoil(t) = bgWsoil(t) + NoahMP50_struc(n)%noahmp50(t)%smc(j)*&
                   NoahMP50_struc(n)%sldpth(j)
           enddo
        endif
     enddo
  endif
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     
     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     i = LIS_domain(n)%gindex(c,r)

     SOILTYP = NoahMP50_struc(n)%noahmp50(t)%soiltype
     MAX_THRESHOLD = NoahMP50_struc(n)%noahmp50(t)%param%SMCMAX(1)  !SMCMAX_TABLE(SOILTYP)

     !locations with large soil moisture values are ice points.
     !we turn off the increments in such locations.
     !For peat, SM routinely exceeds 0.50 (high porosity), so only check against SMCMAX.
     if(SOILTYP.eq.PEAT_SOILTYPE) then
        if(NoahMP50_struc(n)%noahmp50(t)%smc(1).gt.MAX_THRESHOLD) &
             largeSM(i) = .true.
     else
        if(NoahMP50_struc(n)%noahmp50(t)%smc(1).gt.MAX_THRESHOLD.or.&
             NoahMP50_struc(n)%noahmp50(t)%smc(1).gt.0.50) &
             largeSM(i) = .true.
     endif

     if(NoahMP50_struc(n)%noahmp50(t)%snowh.gt.0) then
        snodens(t) = NoahMP50_struc(n)%noahmp50(t)%sneqv/&
             NoahMP50_struc(n)%noahmp50(t)%snowh
     else
        snodens(t) = 0.0
     endif

  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     
     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     i = LIS_domain(n)%gindex(c,r)
     if(largeSM(i)) then 
        soilm1(t) = NoahMP50_struc(n)%noahmp50(t)%smc(1)
        soilm2(t) = NoahMP50_struc(n)%noahmp50(t)%smc(2)
        soilm3(t) = NoahMP50_struc(n)%noahmp50(t)%smc(3)
        soilm4(t) = NoahMP50_struc(n)%noahmp50(t)%smc(4)
        gws(t)    = NoahMP50_struc(n)%noahmp50(t)%wa
        swe(t)    = NoahMP50_struc(n)%noahmp50(t)%sneqv
     endif
  enddo
  
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row
     i = LIS_domain(n)%gindex(c,r)
     
     SOILTYP = NoahMP50_struc(n)%noahmp50(t)%soiltype        
     MAX_THRESHOLD = NoahMP50_struc(n)%noahmp50(t)%param%SMCMAX(1)  !SMCMAX_TABLE(SOILTYP) 
     MIN_THRESHOLD = 0.02 !SMCWLT_TABLE(SOILTYP)
     sm_threshold  = MAX_THRESHOLD

     if((soilm1(t).lt.MIN_THRESHOLD.or.&
          soilm1(t).gt.MAX_THRESHOLD).or.&
          (soilm2(t).lt.MIN_THRESHOLD.or.&
          soilm2(t).gt.MAX_THRESHOLD).or.&
          (soilm3(t).lt.MIN_THRESHOLD.or.&
          soilm3(t).gt.MAX_THRESHOLD).or.&
          (soilm4(t).lt.MIN_THRESHOLD.or.&
          soilm4(t).gt.MAX_THRESHOLD)) then
        ensCheck(i) = .false.
     endif
     ! The groundwater (wa) bound is not applied to peat tiles: wa is not a
     ! peat prognostic and is not updated for peat, so it must not flag a
     ! valid (<= SMCMAX) saturated peat column for ensemble reordering.
     if(SOILTYP.ne.PEAT_SOILTYPE) then
        if(gws(t).lt.MIN_GWS_THRESHOLD.or.gws(t).gt.MAX_GWS_THRESHOLD) then
           ensCheck(i) = .false.
        endif
     endif
     if((soilm1(t).ne.soilm1(i*LIS_rc%nensem(n))).and.&
          (soilm2(t).ne.soilm2(i*LIS_rc%nensem(n))).and.&
          (soilm3(t).ne.soilm3(i*LIS_rc%nensem(n))).and.&
          (soilm4(t).ne.soilm4(i*LIS_rc%nensem(n))).and.&
          (gws(t).ne.gws(i*LIS_rc%nensem(n)))) then 
        diffCheck(i) = .true.
     endif
  enddo

  do i=1,LIS_rc%ngrid(n)
     rc1 = .true.
     rc2 = .true.
     rc3 = .true.
     rc4 = .true.
     rc5 = .true. 
     if(.not.ensCheck(i).and.diffCheck(i).and.(.not.largeSM(i))) then
        call noahmp50_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             soilm1((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD,rc1)
        call noahmp50_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             soilm2((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD,rc2)
        call noahmp50_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             soilm3((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD,rc3)
        call noahmp50_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             soilm4((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_THRESHOLD, MAX_THRESHOLD,rc4)
        call noahmp50_tws_reorderEnsForOutliers(i,&
             LIS_rc%nensem(n),&
             gws((i-1)*LIS_rc%nensem(n)+1:i*LIS_rc%nensem(n)),&
             MIN_GWS_THRESHOLD, MAX_GWS_THRESHOLD,rc5)
     endif
     if(.not.rc1.or.&
          .not.rc2.or.&
          .not.rc3.or.&
          .not.rc4.or.&
          .not.rc5) then

        do m=1,LIS_rc%nensem(n)
           t = (i-1)*LIS_rc%nensem(n)+m
           soilm1(t) = NoahMP50_struc(n)%noahmp50(t)%smc(1)
           soilm2(t) = NoahMP50_struc(n)%noahmp50(t)%smc(2)
           soilm3(t) = NoahMP50_struc(n)%noahmp50(t)%smc(3)
           soilm4(t) = NoahMP50_struc(n)%noahmp50(t)%smc(4)
           gws(t)    = NoahMP50_struc(n)%noahmp50(t)%wa
           swe(t)    = NoahMP50_struc(n)%noahmp50(t)%sneqv
        enddo
     endif
  enddo
        
  ! Per-layer increment clipping (replaces the former all-or-nothing,
  ! per-grid-cell update_flag rejection).  Each soil layer absorbs the
  ! physically-admissible portion of its analysis increment independently,
  ! bounded by [MIN_THRESHOLD, SMCMAX] on the liquid water content and
  ! without ever being pushed further outside those bounds than it already
  ! is.  This lets unsaturated layers be updated even when other layers in
  ! the same column (e.g. deep, fully-saturated peat) are already at the upper
  ! bound, instead of discarding the whole-cell update.  For peat tiles the
  ! resulting soil-moisture change is propagated to the water table depth (and
  ! hence the diagnostic surface water storage) below via increment->zwt
  ! routing, since smc is re-equilibrated from zwt at the next physics step.
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     MAX_THRESHOLD = NoahMP50_struc(n)%noahmp50(t)%param%SMCMAX(1)
     MIN_THRESHOLD = 0.02
     sm_threshold  = MAX_THRESHOLD

     call noahmp50_tws_applyClippedIncr(soilm1(t),&
          NoahMP50_struc(n)%noahmp50(t)%smc(1),&
          NoahMP50_struc(n)%noahmp50(t)%sh2o(1),MIN_THRESHOLD,sm_threshold)
     call noahmp50_tws_applyClippedIncr(soilm2(t),&
          NoahMP50_struc(n)%noahmp50(t)%smc(2),&
          NoahMP50_struc(n)%noahmp50(t)%sh2o(2),MIN_THRESHOLD,sm_threshold)
     call noahmp50_tws_applyClippedIncr(soilm3(t),&
          NoahMP50_struc(n)%noahmp50(t)%smc(3),&
          NoahMP50_struc(n)%noahmp50(t)%sh2o(3),MIN_THRESHOLD,sm_threshold)
     call noahmp50_tws_applyClippedIncr(soilm4(t),&
          NoahMP50_struc(n)%noahmp50(t)%smc(4),&
          NoahMP50_struc(n)%noahmp50(t)%sh2o(4),MIN_THRESHOLD,sm_threshold)

     ! Push the part of the soil-moisture increment that the (clipped) soil
     ! layers could not absorb onto SWE, mirroring the previous fallback but
     ! using the actual unabsorbed remainder.  After clipping, the unabsorbed
     ! volumetric increment of a layer is (soilm - smc); convert to mm and
     ! sum.  Only applied where there is substantial snow and the remainder
     ! is small relative to it (snow-free tiles, e.g. peat, are zeroed in the
     ! snow-update loop below, so the remainder is simply lost there).
     remainder = ((soilm1(t)-NoahMP50_struc(n)%noahmp50(t)%smc(1))*&
          NoahMP50_struc(n)%sldpth(1)+&
          (soilm2(t)-NoahMP50_struc(n)%noahmp50(t)%smc(2))*&
          NoahMP50_struc(n)%sldpth(2)+&
          (soilm3(t)-NoahMP50_struc(n)%noahmp50(t)%smc(3))*&
          NoahMP50_struc(n)%sldpth(3)+&
          (soilm4(t)-NoahMP50_struc(n)%noahmp50(t)%smc(4))*&
          NoahMP50_struc(n)%sldpth(4))*LIS_CONST_RHOFW

     if(NoahMP50_struc(n)%noahmp50(t)%sneqv.gt.5.0.and.&
          swe(t)+remainder.gt.0.0) then
        if(abs(remainder)/NoahMP50_struc(n)%noahmp50(t)%sneqv.lt.0.10) then
           swe(t) = swe(t) + remainder
        endif
     endif

     if(NoahMP50_struc(n)%noahmp50(t)%soiltype.ne.PEAT_SOILTYPE) then
        NoahMP50_struc(n)%noahmp50(t)%wa = gws(t)
     endif

  enddo

!  if(LIS_localPet.eq.387) then
!     gid = LIS_domain(n)%gindex(&
!          LIS_surface(n,LIS_rc%lsm_index)%tile(16068)%col,&
!          LIS_surface(n,LIS_rc%lsm_index)%tile(16068)%row)
!     print*, 'tw2 ',NoahMP50_struc(n)%noahmp50(16068)%smc,&
!          NoahMP50_struc(n)%noahmp50(16068)%sh2o,&
!          NoahMP50_struc(n)%noahmp50(16068)%sneqv,&
!          update_flag(gid)
!  endif

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(snodens(t).eq.0) then
        swe(t) = 0.0
     endif
     dsneqv =  swe(t) - NoahMP50_struc(n)%noahmp50(t)%sneqv

     snowh_new = 0
     if(snodens(t).gt.0) then
        snowh_new = swe(t)/snodens(t)
     endif

     dsnowh = snowh_new - NoahMP50_struc(n)%noahmp50(t)%sneqv

     call noahmp50_snow_update(n, t, dsneqv, dsnowh)
  enddo

!-----------------------------------------------------------------------------------------
! Peatland: route the APPLIED soil-moisture increment to the water table depth (zwt).
! In the PEATCLSM/microtopography scheme zwt is the prognostic memory variable and smc is
! reset to the hydrostatic-equilibrium profile every physics step, so an increment that is
! not carried into zwt is erased.  We move zwt only by the equilibrium response to the
! change in column soil water actually applied this analysis (anWsoil - bgWsoil), as a
! difference of equilibrium inversions: a zero increment leaves zwt untouched (an OL run is
! unaffected) and any systematic offset between the stored profile and the storage operator
! cancels.  Surface (ponded) storage follows from zwt diagnostically.
!-----------------------------------------------------------------------------------------
  if(NoahMP50_struc(n)%peat_opt.eq.1) then
     z_col_bot = sum(NoahMP50_struc(n)%sldpth(1:NoahMP50_struc(n)%nsoil))
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        if(NoahMP50_struc(n)%noahmp50(t)%soiltype.eq.PEAT_SOILTYPE) then
           W_soil = 0.0
           do j=1,NoahMP50_struc(n)%nsoil
              W_soil = W_soil + NoahMP50_struc(n)%noahmp50(t)%smc(j)*&
                   NoahMP50_struc(n)%sldpth(j)
           enddo
           if(abs(W_soil-bgWsoil(t)).gt.0.0) then
              ! copy Campbell params (kind_noahmp) into default-real locals for the call
              thetas_l = NoahMP50_struc(n)%noahmp50(t)%param%SMCMAX(1)
              he_l     = abs(NoahMP50_struc(n)%noahmp50(t)%param%PSISAT(1))
              b_l      = NoahMP50_struc(n)%noahmp50(t)%param%BEXP(1)
              call CalcWaterTableFromSoil_m(bgWsoil(t), thetas_l, he_l, b_l,&
                   z_col_bot, NoahMP50_struc(n)%noahmp50(t)%zwt, wtd_bg)
              call CalcWaterTableFromSoil_m(W_soil, thetas_l, he_l, b_l,&
                   z_col_bot, NoahMP50_struc(n)%noahmp50(t)%zwt, wtd_an)
              NoahMP50_struc(n)%noahmp50(t)%zwt = &
                   NoahMP50_struc(n)%noahmp50(t)%zwt + (wtd_an - wtd_bg)

              ! In the wet/equilibrium regime (water table shallower than the
              ! equilibrium threshold) the peat physics assumes hydrostatic
              ! equilibrium, so reset the profile here to the equilibrium for the
              ! new zwt -- keeping smc and zwt mutually consistent and removing
              ! the soil<->surface double-representation that otherwise lets them
              ! diverge.  In the deeper Richards regime smc is a genuine
              ! prognostic and is left as is.
              if(NoahMP50_struc(n)%noahmp50(t)%zwt.lt.WTD_EQUIL_THRESHOLD) then
                 call CalcEquilibriumProfile_m(NoahMP50_struc(n)%noahmp50(t)%zwt,&
                      thetas_l, he_l, b_l, z_col_bot,&
                      NoahMP50_struc(n)%nsoil, NoahMP50_struc(n)%sldpth, sm_eq_arr)
                 do j=1,NoahMP50_struc(n)%nsoil
                    ice = NoahMP50_struc(n)%noahmp50(t)%smc(j) - &
                          NoahMP50_struc(n)%noahmp50(t)%sh2o(j)
                    NoahMP50_struc(n)%noahmp50(t)%sh2o(j) = sm_eq_arr(j)
                    NoahMP50_struc(n)%noahmp50(t)%smc(j)  = sm_eq_arr(j) + ice
                 enddo
              endif
           endif
        endif
     enddo
  endif

!  write(101,fmt='(I4.4, 1x, I2.2, 1x, I2.2, 1x, I2.2, 1x, I2.2,1x,10E14.6)') &
!       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr,LIS_rc%mn,&
!       NoahMP50_struc(n)%noahmp50(991:1000)%sneqv

end subroutine noahmp50_settws


subroutine noahmp50_tws_reorderEnsForOutliers(i,nensem, statevec, &
     minvalue,maxvalue, status)
  
  implicit none
  integer              :: i
  integer              :: nensem
  real                 :: statevec(nensem)
  real                 :: minvalue,maxvalue
  logical              :: status
  
  real                 :: minvT, maxvT, minvG, maxvG
  integer              :: k
  real                 :: spread_total, spread_good, spread_ratio
  
  !Ensemble spread (total and with 'good' ensemble members
  minvT = 1E10
  maxvT = -1E10
  minvG = 1E10
  maxvG = -1E10
  status = .true. 
  
  do k=1,nensem

     if(statevec(k).lt.minvT) then 
        minvT = statevec(k)
     endif
     if(statevec(k).gt.maxvT) then 
        maxvT = statevec(k)
     endif

     if(statevec(k).gt.minvalue.and.statevec(k).lt.maxvalue) then 
        if(statevec(k).lt.minvG) then 
           minvG = statevec(k)
        endif
        if(statevec(k).gt.maxvG) then 
           maxvG = statevec(k)
        endif
     endif
  enddo
  
  if(minvG.eq.1E10.and.maxvG.eq.-1E10) then
     !all members are unphysical.

     statevec = minvalue
     status = .false.
     
  else
     spread_total = (maxvT - minvT)
     spread_good  = (maxvG - minvG)
     
     spread_ratio = spread_good/spread_total
     
     !rescale the ensemble 
     
     do k=1,nensem-1
        statevec(k) = statevec(nensem) + &
             (statevec(k) - statevec(nensem))*spread_ratio 
     enddo
  endif

end subroutine noahmp50_tws_reorderEnsForOutliers


!BOP
! !ROUTINE: noahmp50_tws_applyClippedIncr
!  \label{noahmp50_tws_applyClippedIncr}
!
! !INTERFACE:
subroutine noahmp50_tws_applyClippedIncr(soilm, smc, sh2o, lo, hi)
! !DESCRIPTION:
!  Apply the admissible portion of one soil layer's analysis increment.
!  The proposed increment (soilm - smc) is added to both the total (smc)
!  and liquid (sh2o) soil moisture, with the liquid water clipped so that
!  it is never driven further outside the admissible range [lo,hi] than it
!  already is.  A layer already above hi (e.g. near-saturated peat) thus
!  accepts only drying increments, and a layer already below lo only
!  wetting increments, while in-range layers are clipped to [lo,hi].
!
!   soilm : analysis (proposed) total soil moisture for the layer [in]
!   smc   : model total soil moisture for the layer [inout]
!   sh2o  : model liquid soil moisture for the layer [inout]
!   lo,hi : admissible liquid-water bounds (MIN_THRESHOLD, SMCMAX-0.02)
!EOP
  implicit none
  real :: soilm
  real :: smc
  real :: sh2o
  real :: lo
  real :: hi
  real :: delta, cur, new_sh2o, applied

  delta    = soilm - smc
  cur      = sh2o
  new_sh2o = cur + delta
  if(new_sh2o .gt. max(hi,cur)) new_sh2o = max(hi,cur)
  if(new_sh2o .lt. min(lo,cur)) new_sh2o = min(lo,cur)
  applied  = new_sh2o - cur
  smc  = smc  + applied
  sh2o = sh2o + applied

end subroutine noahmp50_tws_applyClippedIncr
