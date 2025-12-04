module WaterStateBulkType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Water state variables
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar, only : nlevgrnd, nlevsno
  use clm_varcon, only : ispval, nan => spval
  use decompMod, only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  !PUBLIC DATA TYPES:

  type, public :: waterstatebulk_type

    real(r8), pointer :: h2osoi_liq_col   (:,:)  ! col liquid water (kg H2O/m2) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_ice_col   (:,:)  ! col ice lens (kg H2O/m2) (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: h2osoi_vol_col   (:,:)  ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
    real(r8), pointer :: h2osfc_col       (:)    ! col surface water (mm H2O)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type waterstatebulk_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(waterstatebulk_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(waterstatebulk_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    integer :: begc, endc   ! Column indices
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc

    allocate (this%h2osoi_liq_col   (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq_col   (:,:) = nan
    allocate (this%h2osoi_ice_col   (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice_col   (:,:) = nan
    allocate (this%h2osoi_vol_col   (begc:endc,1:nlevgrnd))          ; this%h2osoi_vol_col   (:,:) = nan
    allocate (this%h2osfc_col       (begc:endc))                     ; this%h2osfc_col       (:)   = nan

  end subroutine InitAllocate

end module WaterStateBulkType
