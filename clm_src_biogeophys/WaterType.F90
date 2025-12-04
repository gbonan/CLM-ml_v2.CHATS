module WaterType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Water state variables
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varcon, only : ispval, nan => spval
  use decompMod, only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  !PUBLIC DATA TYPES:

  type, public :: water_type

    real(r8), pointer :: h2osno_col(:)    ! col snow water (mm H2O)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type water_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(water_type) :: this
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
    class(water_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    integer :: begc, endc   ! Column indices
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc

    allocate (this%h2osno_col (begc:endc)) ; this%h2osno_col(:) = nan

  end subroutine InitAllocate

end module WaterType
