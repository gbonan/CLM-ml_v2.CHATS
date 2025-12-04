module MLpftconMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing vegetation (PFT) parameters unique to the multilayer
  ! canopy and method to read and initialize them. Additional PFT parameters
  ! are provided by CLM in pftcon.
  !
  ! !USES:
  use spmdMod, only : masterproc
  use clm_varctl, only : iulog
  use clm_varpar, only : mxpft
  use MLclm_varctl, only : pftcon_val
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! PFT parameters
  !
  type, public :: MLpftcon_type

    real(r8), allocatable :: vcmaxpft         (:)   ! Maximum carboxylation rate at 25C (umol/m2/s)
    real(r8), allocatable :: gplant_SPA       (:)   ! Stem (xylem-to-leaf) hydraulic conductance (mmol H2O/m2 leaf area/s/Mpa)
    real(r8), allocatable :: capac_SPA        (:)   ! Plant capacitance (mmol H2O/m2 leaf area/MPa)
    real(r8), allocatable :: iota_SPA         (:)   ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    real(r8), allocatable :: root_radius_SPA  (:)   ! Fine root radius (m)
    real(r8), allocatable :: root_density_SPA (:)   ! Fine root density (g biomass / m3 root)
    real(r8), allocatable :: root_resist_SPA  (:)   ! Hydraulic resistivity of root tissue (MPa.s.g/mmol H2O)
    real(r8), allocatable :: gsmin_SPA        (:)   ! Minimum stomatal conductance (mol H2O/m2/s)
    real(r8), allocatable :: g0_BB            (:)   ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
    real(r8), allocatable :: g1_BB            (:)   ! Ball-Berry slope of conductance-photosynthesis relationship
    real(r8), allocatable :: g0_MED           (:)   ! Medlyn minimum leaf conductance (mol H2O/m2/s)
    real(r8), allocatable :: g1_MED           (:)   ! Medlyn slope of conductance-photosynthesis relationship
    real(r8), allocatable :: psi50_gs         (:)   ! Leaf water potential at which 50% of stomatal conductance is lost (MPa)
    real(r8), allocatable :: shape_gs         (:)   ! Shape parameter for stomatal conductance in relation to leaf water potential (-)
    real(r8), allocatable :: emleaf           (:)   ! Leaf emissivity (-)
    real(r8), allocatable :: clump_fac        (:)   ! Foliage clumping index (-)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, private :: InitRead

  end type MLpftcon_type

  type(MLpftcon_type), public :: MLpftcon
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this)

    class(MLpftcon_type) :: this

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
    class(MLpftcon_type) :: this
    !---------------------------------------------------------------------

    allocate (this%vcmaxpft         (0:mxpft))
    allocate (this%gplant_SPA       (0:mxpft))
    allocate (this%capac_SPA        (0:mxpft))
    allocate (this%iota_SPA         (0:mxpft))
    allocate (this%root_radius_SPA  (0:mxpft))
    allocate (this%root_density_SPA (0:mxpft))
    allocate (this%root_resist_SPA  (0:mxpft))
    allocate (this%gsmin_SPA        (0:mxpft))
    allocate (this%g0_BB            (0:mxpft))
    allocate (this%g1_BB            (0:mxpft))
    allocate (this%g0_MED           (0:mxpft))
    allocate (this%g1_MED           (0:mxpft))
    allocate (this%psi50_gs         (0:mxpft))
    allocate (this%shape_gs         (0:mxpft))
    allocate (this%emleaf           (0:mxpft))
    allocate (this%clump_fac        (0:mxpft))

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitRead (this)
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) parameters
    !
    ! !ARGUMENTS:
    class(MLpftcon_type) :: this
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

    if (masterproc) then
       write (iulog,*) 'Attempting to initialize MLpftcon .....'
    end if

    ! vcmax

    this%vcmaxpft(:) = -999._r8
    this%vcmaxpft( 1) = 62.5_r8
    this%vcmaxpft( 2) = 62.5_r8
    this%vcmaxpft( 3) = 39.1_r8
    this%vcmaxpft( 4) = 41.0_r8
    this%vcmaxpft( 5) = 61.4_r8
    this%vcmaxpft( 6) = 41.0_r8
    this%vcmaxpft( 7) = 57.7_r8
    this%vcmaxpft( 8) = 57.7_r8
    this%vcmaxpft( 9) = 61.7_r8
    this%vcmaxpft(10) = 54.0_r8
    this%vcmaxpft(11) = 54.0_r8
    this%vcmaxpft(12) = 78.2_r8
    this%vcmaxpft(13) = 78.2_r8
    this%vcmaxpft(14) = 51.6_r8
    this%vcmaxpft(15) = 100.7_r8
    this%vcmaxpft(16) = 100.7_r8

    ! Plant hydraulics

    this%gplant_SPA(:) = -999._r8
    this%gplant_SPA(1:16) = 4._r8

    this%capac_SPA(:) = -999._r8
    this%capac_SPA( 1:11) = 2500._r8
    this%capac_SPA(12:16) = 500._r8

    ! Stomatal optimization

    this%iota_SPA(:) = -999._r8
    this%iota_SPA(1: 1) = 750._r8
    this%iota_SPA(2: 3) = 1500._r8
    this%iota_SPA(4: 4) = 500._r8
    this%iota_SPA(5:16) = 750._r8

    ! Root hydraulics

    this%root_radius_SPA(:) = -999._r8
    this%root_density_SPA(:) = -999._r8
    this%root_resist_SPA(:) = -999._r8

    this%root_radius_SPA(1:16) = 0.29e-03_r8
    this%root_density_SPA(1:16) = 0.31e06_r8
    this%root_resist_SPA(1:16) = 25._r8

    ! Minimum stomatal conductance

    this%gsmin_SPA(:)= -999._r8
    this%gsmin_SPA(1:16)= 0.002_r8

    ! Ball-Berry stomatal conductance parameters

    this%g0_BB(:)= -999._r8
    this%g1_BB(:)= -999._r8

    this%g0_BB( 1:13) = 0.01_r8
    this%g0_BB(14:14) = 0.04_r8
    this%g0_BB(15:16) = 0.01_r8

    this%g1_BB( 1:13) = 9._r8
    this%g1_BB(14:14) = 4._r8
    this%g1_BB(15:16) = 9._r8

    ! Medlyn stomatal conductance parameters

    this%g0_MED(:)= -999._r8
    this%g1_MED(:)= -999._r8

    this%g0_MED(1:16) = 0.0001_r8

    this%g1_MED( 1) = 2.35_r8
    this%g1_MED( 2) = 2.35_r8
    this%g1_MED( 3) = 2.35_r8
    this%g1_MED( 4) = 4.12_r8
    this%g1_MED( 5) = 4.12_r8
    this%g1_MED( 6) = 4.45_r8
    this%g1_MED( 7) = 4.45_r8
    this%g1_MED( 8) = 4.45_r8
    this%g1_MED( 9) = 4.70_r8
    this%g1_MED(10) = 4.70_r8
    this%g1_MED(11) = 4.70_r8
    this%g1_MED(12) = 2.22_r8
    this%g1_MED(13) = 5.25_r8
    this%g1_MED(14) = 1.62_r8
    this%g1_MED(15) = 5.79_r8
    this%g1_MED(16) = 5.79_r8

    ! Leaf water potential at which 50% of stomatal conductance is lost

    this%psi50_gs(:) = -999._r8
    this%psi50_gs(1:16) = -2.3_r8

    ! Shape parameter for stomatal conductance in relation to leaf water potential

    this%shape_gs(:) = -999._r8
    this%shape_gs(1:16) = 40._r8

    ! Leaf emissivity

    this%emleaf(:) = -999._r8
    this%emleaf(1:16) = 0.98_r8

    ! Foliage clumping index

    this%clump_fac(:) = -999._r8
    this%clump_fac(1:16) = 1._r8

    ! Tower site adjustments

    if (pftcon_val == 1) then
       write (iulog,*) 'MLpftcon ... using non-default values'

       this%vcmaxpft(7) = 125._r8            ! CHATS: Rosati et al. (2006)
       this%gplant_SPA(7) = 7._r8            ! CHATS: Tyree et al. (1993)
       this%iota_SPA(7) = 375._r8            ! CHATS: Rosati et al. (2006)
       this%root_resist_SPA(7) = 14._r8      ! CHATS: Tyree et al. (1994)

!      this%psi50_gs(7) = -0.45_r8           ! CHATS: Rosati et al. (2006)
!      this%shape_gs(7) = 2.60_r8            ! CHATS: Rosati et al. (2006)
!      this%psi50_gs(7) = -1.21_r8           ! CHATS: Cochard et al. (2002)
!      this%shape_gs(7) = 6.05_r8            ! CHATS: Cochard et al. (2002)
       this%psi50_gs(7) = -1.60_r8           ! CHATS: SPA-walnut
       this%shape_gs(7) = 40.0_r8            ! CHATS: SPA-walnut
    end if

    if (masterproc) then
       write (iulog,*) 'Successfuly initialized MLpftcon'
    end if

  end subroutine InitRead

end module MLpftconMod
