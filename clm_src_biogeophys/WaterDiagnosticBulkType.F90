module WaterDiagnosticBulkType

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

  type, public :: waterdiagnosticbulk_type

    real(r8), pointer :: q_ref2m_patch(:)      ! patch 2 m height surface specific humidity (kg/kg)
    real(r8), pointer :: frac_sno_eff_col(:)   ! col fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: bw_col(:,:)           ! col partial density of water in the snow pack (ice + liquid) [kg/m3]

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type waterdiagnosticbulk_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(waterdiagnosticbulk_type) :: this
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
    class(waterdiagnosticbulk_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    integer :: begc, endc   ! Column indices
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc

    allocate (this%q_ref2m_patch    (begp:endp))                     ; this%q_ref2m_patch    (:)   = nan
    allocate (this%frac_sno_eff_col (begc:endc))                     ; this%frac_sno_eff_col (:)   = nan
    allocate (this%bw_col           (begc:endc,-nlevsno+1:0))        ; this%bw_col           (:,:) = nan

  end subroutine InitAllocate

end module WaterDiagnosticBulkType
