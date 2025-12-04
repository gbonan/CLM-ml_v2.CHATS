module MLLeafBoundaryLayerMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf boundary layer conductance
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LeafBoundaryLayer
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LeafBoundaryLayer (num_filter, filter, il, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Leaf boundary layer conductance
    ! See Bonan (2019) Climate Change and Terrestrial Ecosystem Modeling (Chapter 10)
    !
    ! !USES:
    use clm_varcon, only : tfrz, grav
    use PatchType, only : patch
    use pftconMod, only : pftcon
    use MLclm_varcon, only : visc0, dh0, dv0, dc0, gb_factor, gbh_min
    use MLclm_varctl, only : gb_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter        ! Number of patches in filter
    integer, intent(in) :: filter(:)         ! Patch filter
    integer, intent(in) :: il                ! Sunlit (1) or shaded (2) leaf index
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                           ! Filter index
    integer  :: p                            ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                           ! Aboveground layer index
    real(r8) :: visc                         ! Kinematic viscosity (m2/s)
    real(r8) :: dh                           ! Molecular diffusivity, heat (m2/s)
    real(r8) :: dv                           ! Molecular diffusivity, H2O (m2/s)
    real(r8) :: dc                           ! Molecular diffusivity, CO2 (m2/s)
    real(r8) :: fac                          ! Correction factor for temperature and pressure
    real(r8) :: nu                           ! Nusselt number (dimensionless)
    real(r8) :: pr                           ! Prandtl number (dimensionless)
    real(r8) :: re                           ! Reynolds number (dimensionless)
    real(r8) :: gr                           ! Grashof number (dimensionless)
    real(r8) :: gbh_lam, gbv_lam, gbc_lam    ! Forced convection - laminar: conductances (mol/m2/s)
    real(r8) :: gbh_turb, gbv_turb, gbc_turb ! Forced convection - turbulent: conductances (mol/m2/s)
    real(r8) :: gbh_free, gbv_free, gbc_free ! Free convection: conductances (mol/m2/s)
    !---------------------------------------------------------------------

    associate ( &
                                                   ! *** Input ***
    dleaf     => pftcon%dleaf                 , &  ! CLM: Leaf dimension (m)
    tref      => mlcanopy_inst%tref_forcing   , &  ! Air temperature at reference height (K)
    pref      => mlcanopy_inst%pref_forcing   , &  ! Air pressure at reference height (Pa)
    rhomol    => mlcanopy_inst%rhomol_forcing , &  ! Molar density at reference height (mol/m3)
    ncan      => mlcanopy_inst%ncan_canopy    , &  ! Number of aboveground layers
    dpai      => mlcanopy_inst%dpai_profile   , &  ! Canopy layer plant area index (m2/m2)
    wind      => mlcanopy_inst%wind_profile   , &  ! Canopy layer wind speed (m/s)
    tair      => mlcanopy_inst%tair_profile   , &  ! Canopy layer air temperature (K)
    tleaf     => mlcanopy_inst%tleaf_leaf     , &  ! Leaf temperature (K)
                                                   ! *** Output ***
    gbh       => mlcanopy_inst%gbh_leaf       , &  ! Leaf boundary layer conductance: heat (mol/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv_leaf       , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    gbc       => mlcanopy_inst%gbc_leaf         &  ! Leaf boundary layer conductance: CO2 (mol CO2/m2 leaf/s)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Adjust diffusivity for temperature and pressure

       fac = 101325._r8 / pref(p) * (tref(p) / tfrz)**1.81_r8
       visc = visc0 * fac
       dh = dh0 * fac
       dv = dv0 * fac
       dc = dc0 * fac

       do ic = 1, ncan(p)

          if (dpai(p,ic) > 0._r8) then

             ! Reynolds number, Prandtl number, and Grashof number

             re = wind(p,ic) * dleaf(patch%itype(p)) / visc
             pr  = visc / dh
             gr = grav * dleaf(patch%itype(p))**3 * max(tleaf(p,ic,il)-tair(p,ic), 0._r8) / (tair(p,ic) * visc * visc)

             ! Nusselt number depends on convection regime

             ! a. Forced convection
             ! Note use of minimum conductance applied to gbh (needed for low wind speed)

             ! (i) Laminar flow

                nu = gb_factor * 0.66_r8 * pr**0.33_r8 * re**0.5_r8
                gbh_lam = (dh * nu / dleaf(patch%itype(p))) * rhomol(p) ; gbh_lam = max(gbh_lam, gbh_min)
                gbv_lam = gbh_lam * (dv / dh)**0.67_r8
                gbc_lam = gbh_lam * (dc / dh)**0.67_r8

             ! (ii) Turbulent flow

                nu = gb_factor * 0.036_r8 * pr**0.33_r8 * re**0.8_r8
                gbh_turb = (dh * nu / dleaf(patch%itype(p))) * rhomol(p) ; gbh_turb = max(gbh_turb, gbh_min)
                gbv_turb = gbh_turb * (dv / dh)**0.67_r8
                gbc_turb = gbh_turb * (dc / dh)**0.67_r8

             ! b. Free convection

                nu = 0.54_r8 * pr**0.25_r8 * gr**0.25_r8
                gbh_free = (dh * nu / dleaf(patch%itype(p))) * rhomol(p)
                gbv_free = gbh_free * (dv / dh)**0.75_r8
                gbc_free = gbh_free * (dc / dh)**0.75_r8

             ! Choose flow regimes to use

             select case (gb_type)
             case (1)

                ! Use only laminar flow

                gbh(p,ic,il) = gbh_lam
                gbv(p,ic,il) = gbv_lam
                gbc(p,ic,il) = gbc_lam

             case (2)

                ! Use laminar and turbulent flow

                gbh(p,ic,il) = max(gbh_lam, gbh_turb)
                gbv(p,ic,il) = max(gbv_lam, gbv_turb)
                gbc(p,ic,il) = max(gbc_lam, gbc_turb)

             case (3)

                ! Both forced and free convection occur together

                gbh(p,ic,il) = max(gbh_lam, gbh_turb) + gbh_free
                gbv(p,ic,il) = max(gbv_lam, gbv_turb) + gbv_free
                gbc(p,ic,il) = max(gbc_lam, gbc_turb) + gbc_free

             case default

                call endrun (msg=' ERROR: LeafBoundaryLayer: gb_type not valid')

             end select

          else

             gbh(p,ic,il) = 0._r8
             gbv(p,ic,il) = 0._r8
             gbc(p,ic,il) = 0._r8

          end if

       end do
    end do

    end associate
  end subroutine LeafBoundaryLayer

end module MLLeafBoundaryLayerMod
