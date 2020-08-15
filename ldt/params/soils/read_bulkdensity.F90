!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_STATSGOFAO_texture
! \label{read_STATSGOFAO_texture}
!
! !REVISION HISTORY:
!  05 March  2019: Maertens Michiel ; Initial Specification
!
! !INTERFACE:
subroutine read_bulkdensity(n,array)

! !USES:
  use LDT_coreMod,    only : LDT_rc
  use LDT_logMod,     only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,  only : readLISdata, LDT_checkDomainExtents 
  use LDT_paramTileInputMod, only: param_index_fgrdcalc

  implicit none
! !ARGUMENTS: 
  integer,intent(in)    :: n
 real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
!
! !DESCRIPTION:
!  This subroutine retrieves HWSD-STASGO pisat  data and reprojects
!  it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved soil texture
!  \end{description}
!EOP
  integer :: c, r, t 
  integer :: water_class
  integer :: ftn
  real    :: temp(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
 

 ! real    :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  logical :: file_exists
! ___________________________________________________________________
  
    

  inquire(file=trim(LDT_rc%bulkdensityfile(n)), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "bulkdensity  map ",trim(LDT_rc%weltingpointfile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  write(unit=LDT_logunit,fmt=*) "[INFO] Reading bulkdensity file: ",&
                                trim(LDT_rc%bulkdensityfile(n))
!- Open file:
   ftn = LDT_getNextUnitNumber()
   open( ftn,file=LDT_rc%bulkdensityfile(n),form='unformatted',status='old',&
        access='direct',recl=4 )

      call readLISdata(n, ftn, LDT_rc%soils_proj, LDT_rc%soils_gridtransform(n), &
                       LDT_rc%soil_gridDesc(n,:), 1, temp )  ! 1 indicates 2D layer



   array(:,:,:) = temp(:,:,:)

! ---
   call LDT_releaseUnitNumber(ftn)
   write(unit=LDT_logunit,fmt=*) "[INFO] Done reading bukdensity file."

 end subroutine read_bulkdensity

