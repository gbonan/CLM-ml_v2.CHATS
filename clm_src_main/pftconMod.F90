module pftconMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation (PFT) parameters and method to
  ! read and initialize them
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varctl, only : iulog
  use clm_varpar, only : mxpft, numrad, ivis, inir
  use MLclm_varctl, only : pftcon_val
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! PFT parameters
  !
  type, public :: pftcon_type

    ! CLM pft parameters
    real(r8), allocatable :: dleaf            (:)   ! Characteristic leaf dimension (m)
    real(r8), allocatable :: c3psn            (:)   ! Photosynthetic pathway: 0. = C4, 1. = C3
    real(r8), allocatable :: xl               (:)   ! Leaf/stem orientation index
    real(r8), allocatable :: rhol             (:,:) ! Leaf reflectance: 1=vis, 2=nir
    real(r8), allocatable :: rhos             (:,:) ! Stem reflectance: 1=vis, 2=nir
    real(r8), allocatable :: taul             (:,:) ! Leaf transmittance: 1=vis, 2=nir
    real(r8), allocatable :: taus             (:,:) ! Stem transmittance: 1=vis, 2=nir
    real(r8), allocatable :: rootprof_beta    (:)   ! Jackson1996 rooting distribution parameter (-)
    real(r8), allocatable :: slatop           (:)   ! Specific leaf area at top of canopy (m2/gC)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, private :: InitRead

  end type pftcon_type

  type(pftcon_type), public :: pftcon
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this)

    class(pftcon_type) :: this

    call this%InitAllocate()
    call this%InitRead()

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate (this)
    !
    ! !DESCRIPTION:
    ! Allocate memory for pft data structure
    !
    ! !ARGUMENTS:
    class(pftcon_type) :: this
    !---------------------------------------------------------------------

    ! CLM pft parameters
    allocate (this%dleaf            (0:mxpft))
    allocate (this%c3psn            (0:mxpft))
    allocate (this%xl               (0:mxpft))
    allocate (this%rhol             (0:mxpft,numrad))
    allocate (this%rhos             (0:mxpft,numrad))
    allocate (this%taul             (0:mxpft,numrad))
    allocate (this%taus             (0:mxpft,numrad))
    allocate (this%rootprof_beta    (0:mxpft))
    allocate (this%slatop           (0:mxpft))

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitRead (this)
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) parameters
    !
    ! !ARGUMENTS:
    class(pftcon_type) :: this
    !---------------------------------------------------------------------

    ! PFTs
    ! -------------------------------
    !    0 => not_vegetated
    !    1 => needleleaf_evergreen_temperate_tree
    !    2 => needleleaf_evergreen_boreal_tree
    !    3 => needleleaf_deciduous_boreal_tree
    !    4 => broadleaf_evergreen_tropical_tree
    !    5 => broadleaf_evergreen_temperate_tree
    !    6 => broadleaf_deciduous_tropical_tree
    !    7 => broadleaf_deciduous_temperate_tree
    !    8 => broadleaf_deciduous_boreal_tree
    !    9 => broadleaf_evergreen_shrub
    !   10 => broadleaf_deciduous_temperate_shrub
    !   11 => broadleaf_deciduous_boreal_shrub
    !   12 => c3_arctic_grass
    !   13 => c3_non-arctic_grass
    !   14 => c4_grass
    !   15 => c3_crop
    !   16 => c3_irrigated
    !   17 => temperate_corn
    !   18 => irrigated_temperate_corn
    !   19 => spring_wheat
    !   20 => irrigated_spring_wheat
    !   21 => winter_wheat
    !   22 => irrigated_winter_wheat
    !   23 => temperate_soybean
    !   24 => irrigated_temperate_soybean
    !   25 => barley
    !   26 => irrigated_barley
    !   27 => winter_barley
    !   28 => irrigated_winter_barley
    !   29 => rye
    !   30 => irrigated_rye
    !   31 => winter_rye
    !   32 => irrigated_winter_rye
    !   33 => cassava
    !   34 => irrigated_cassava
    !   35 => citrus
    !   36 => irrigated_citrus
    !   37 => cocoa
    !   38 => irrigated_cocoa
    !   39 => coffee
    !   40 => irrigated_coffee
    !   41 => cotton
    !   42 => irrigated_cotton
    !   43 => datepalm
    !   44 => irrigated_datepalm
    !   45 => foddergrass
    !   46 => irrigated_foddergrass
    !   47 => grapes
    !   48 => irrigated_grapes
    !   49 => groundnuts
    !   50 => irrigated_groundnuts
    !   51 => millet
    !   52 => irrigated_millet
    !   53 => oilpalm
    !   54 => irrigated_oilpalm
    !   55 => potatoes
    !   56 => irrigated_potatoes
    !   57 => pulses
    !   58 => irrigated_pulses
    !   59 => rapeseed
    !   60 => irrigated_rapeseed
    !   61 => rice
    !   62 => irrigated_rice
    !   63 => sorghum
    !   64 => irrigated_sorghum
    !   65 => sugarbeet
    !   66 => irrigated_sugarbeet
    !   67 => sugarcane
    !   68 => irrigated_sugarcane
    !   69 => sunflower
    !   70 => irrigated_sunflower
    !   71 => miscanthus
    !   72 => irrigated_miscanthus
    !   73 => switchgrass
    !   74 => irrigated_switchgrass
    !   75 => tropical_corn
    !   76 => irrigated_tropical_corn
    !   77 => tropical_soybean
    !   78 => irrigated_tropical_soybean
    ! -------------------------------

    ! Leaf dimension

    this%dleaf(:) = -999._r8
    this%dleaf(1:16) = 0.04_r8

    ! Photosynthetic pathway: 1. = C3 plant and 0. = C4 plant

    this%c3psn(:) = -999._r8
    this%c3psn( 1:13) = 1._r8
    this%c3psn(14:14) = 0._r8
    this%c3psn(15:16) = 1._r8

    ! Leaf angle

    this%xl(:) = -999._r8
    this%xl( 1: 3) = 0.01_r8
    this%xl( 4: 5) = 0.10_r8
    this%xl( 6: 6) = 0.01_r8
    this%xl( 7: 8) = 0.25_r8
    this%xl( 9: 9) = 0.01_r8
    this%xl(10:11) = 0.25_r8
    this%xl(12:16) = -0.30_r8

    ! Leaf reflectance: visible and near-infrared

    this%rhol(:,:) = -999._r8

    this%rhol( 1: 3,ivis) = 0.07_r8
    this%rhol( 4: 8,ivis) = 0.10_r8
    this%rhol( 9: 9,ivis) = 0.07_r8
    this%rhol(10:11,ivis) = 0.10_r8
    this%rhol(12:16,ivis) = 0.11_r8

    this%rhol( 1: 3,inir) = 0.35_r8
    this%rhol( 4: 8,inir) = 0.45_r8
    this%rhol( 9: 9,inir) = 0.35_r8
    this%rhol(10:11,inir) = 0.45_r8
    this%rhol(12:16,inir) = 0.35_r8

    ! Stem reflectance: visible and near-infrared

    this%rhos(:,:) = -999._r8

    this%rhos( 1:11,ivis) = 0.16_r8
    this%rhos(12:16,ivis) = 0.31_r8

    this%rhos( 1:11,inir) = 0.39_r8
    this%rhos(12:16,inir) = 0.53_r8

    ! Leaf transmittance: visible and near-infrared

    this%taul(:,:) = -999._r8

    this%taul(1:16,ivis) = 0.05_r8

    this%taul( 1: 3,inir) = 0.10_r8
    this%taul( 4: 8,inir) = 0.25_r8
    this%taul( 9: 9,inir) = 0.10_r8
    this%taul(10:11,inir) = 0.25_r8
    this%taul(12:16,inir) = 0.34_r8

    ! Stem transmittance: visible and near-infrared

    this%taus(:,:) = -999._r8

    this%taus( 1:11,ivis) = 0.001_r8
    this%taus(12:16,ivis) = 0.12_r8

    this%taus( 1:11,inir) = 0.001_r8
    this%taus(12:16,inir) = 0.25_r8

    ! Jackson1996 rooting distribution parameters (-)

    this%rootprof_beta(:) = -999._r8
    this%rootprof_beta( 1: 1) = 0.976_r8
    this%rootprof_beta( 2: 3) = 0.943_r8
    this%rootprof_beta( 4: 4) = 0.993_r8
    this%rootprof_beta( 5: 5) = 0.966_r8
    this%rootprof_beta( 6: 6) = 0.993_r8
    this%rootprof_beta( 7: 7) = 0.966_r8
    this%rootprof_beta( 8: 8) = 0.943_r8
    this%rootprof_beta( 9:10) = 0.964_r8
    this%rootprof_beta(11:12) = 0.914_r8
    this%rootprof_beta(13:16) = 0.943_r8

    ! Specific leaf area at top of canopy, projected area basis (m2/gC)

    this%slatop(:) = -999._r8
    this%slatop( 1) = 0.010_r8
    this%slatop( 2) = 0.008_r8
    this%slatop( 3) = 0.024_r8
    this%slatop( 4) = 0.012_r8
    this%slatop( 5) = 0.012_r8
    this%slatop( 6) = 0.030_r8
    this%slatop( 7) = 0.030_r8
    this%slatop( 8) = 0.030_r8
    this%slatop( 9) = 0.012_r8
    this%slatop(10) = 0.030_r8
    this%slatop(11) = 0.030_r8
    this%slatop(12) = 0.030_r8
    this%slatop(13) = 0.030_r8
    this%slatop(14) = 0.030_r8
    this%slatop(15) = 0.030_r8
    this%slatop(16) = 0.030_r8

    ! Tower site adjustments

    if (pftcon_val == 1) then
       write (iulog,*) 'pftconMod ... using non-default values'

       ! Adjust optical properties for Majasalmi & Bright (2019)

!      this%xl(7) = 0.59_r8                  ! CHATS: BDT temperate
       this%xl(7) = 0.53_r8                  ! CHATS: walnut

!      this%rhol(7,ivis) = 0.08_r8           ! CHATS: BDT temperate
       this%rhol(7,ivis) = 0.06_r8           ! CHATS: walnut

!      this%rhol(7,inir) = 0.42_r8           ! CHATS: BDT temperate
       this%rhol(7,inir) = 0.42_r8           ! CHATS: walnut

       this%rhos(7,ivis) = 0.21_r8           ! CHATS: deciduous bark
       this%rhos(7,inir) = 0.49_r8           ! CHATS: deciduous bark

!      this%taul(7,ivis) = 0.06_r8           ! CHATS: BDT temperate
       this%taul(7,ivis) = 0.04_r8           ! CHATS: walnut

!      this%taul(7,inir) = 0.43_r8           ! CHATS: BDT temperate
       this%taul(7,inir) = 0.43_r8           ! CHATS: walnut
    end if

  end subroutine InitRead

end module pftconMod
