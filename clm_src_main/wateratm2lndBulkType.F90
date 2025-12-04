module wateratm2lndBulkType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! atm -> land variables
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar, only : numrad
  use decompMod, only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  !PUBLIC DATA TYPES:

  type, public :: wateratm2lndbulk_type

    ! atm -> land: downscaled to column
    real(r8), pointer :: forc_q_downscaled_col    (:)    ! Atmospheric specific humidity (kg/kg)
    real(r8), pointer :: forc_rain_downscaled_col (:)    ! Rainfall rate (mm/s)
    real(r8), pointer :: forc_snow_downscaled_col (:)    ! Snowfall rate (mm/s)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type wateratm2lndbulk_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(wateratm2lndbulk_type) :: this
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
    class(wateratm2lndbulk_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! Initial value
    integer  :: begg, endg      ! Grid cell indices
    integer  :: begc, endc      ! Column indices
    !---------------------------------------------------------------------

    begg = bounds%begg ; endg = bounds%endg
    begc = bounds%begc ; endc = bounds%endc

    allocate (this%forc_q_downscaled_col    (begc:endc))        ; this%forc_q_downscaled_col    (:)   = ival
    allocate (this%forc_rain_downscaled_col (begc:endc))        ; this%forc_rain_downscaled_col (:)   = ival
    allocate (this%forc_snow_downscaled_col (begc:endc))        ; this%forc_snow_downscaled_col (:)   = ival

  end subroutine InitAllocate

end module wateratm2lndBulkType
