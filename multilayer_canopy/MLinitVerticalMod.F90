module MLinitVerticalMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initialize multilayer canopy vertical structure and profiles
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog
  use decompMod, only : bounds_type
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initVerticalStructure  ! Define canopy layer vertical structure
  public :: initVerticalProfiles   ! Initialize vertical profiles and canopy states
  public :: getPADparameters       ! Initialize parameters for plant area density profile
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initVerticalStructure (bounds, num_filter, filter, &
  canopystate_inst, frictionvel_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Define canopy layer vertical structure
    !
    ! !USES:
    use PatchType, only : patch
    use MLclm_varctl, only : dz_tall, dz_short, dz_param, nlayer_above, nlayer_within, dpai_min
    use MLclm_varpar, only : nlevmlcan
    use MLMathToolsMod, only : beta_distribution_cdf
    use CanopyStateType, only : canopystate_type
    use FrictionVelocityMod, only : frictionvel_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_filter            ! Number of patches in filter
    integer, intent(in) :: filter(:)             ! Patch filter

    type(canopystate_type), intent(in) :: canopystate_inst
    type(frictionvel_type), intent(in) :: frictionvel_inst
    type(mlcanopy_type), intent(out) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                               ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                               ! Aboveground layer index
    integer  :: iflag                            ! Error flag
    real(r8) :: ztop_to_zref                     ! Atmospheric reference height - canopy height (m)
    real(r8) :: dz_within                        ! Height increment within canopy (m)
    real(r8) :: dz_above                         ! Height increment above canopy (m)
    integer  :: nabove                           ! Number of above-canopy layers
    real(r8) :: zrel                             ! Height relative to canopy top
    real(r8) :: beta_cdf_bot, beta_cdf_top       ! Beta distribution: cumulative distribution at bottom and top heights for layer
    real(r8) :: lai_err, sai_err, pai_err        ! Sum of leaf, stem, or plant area index over all layers (m2/m2)
    real(r8) :: lai_miss, sai_miss               ! Missing leaf (stem) area after imposing dpai_min constraint on layers (m2/m2)

    real(r8) :: dlai(bounds%begp:bounds%endp,1:nlevmlcan)  ! Layer leaf area index (m2/m2)
    real(r8) :: dsai(bounds%begp:bounds%endp,1:nlevmlcan)  ! Layer stem area index (m2/m2)

    real(r8) :: unit_lai = 1.0_r8                          ! Unit leaf area index of canopy (m2/m2)
    real(r8) :: unit_sai = 1.0_r8                          ! Unit stem area index of canopy (m2/m2)
    !---------------------------------------------------------------------

    associate ( &
                                                          ! *** Input ***
    forc_hgt_u  => frictionvel_inst%forc_hgt_u_patch , &  ! CLM: Atmospheric reference height (m)
    htop        => canopystate_inst%htop_patch       , &  ! CLM: Canopy height (m)
    pbeta_lai   => mlcanopy_inst%pbeta_lai_canopy    , &  ! Parameters for the leaf area density 2-parameter beta distribution (-)
    pbeta_sai   => mlcanopy_inst%pbeta_sai_canopy    , &  ! Parameters for the stem area density 2-parameter beta distribution (-)
                                                          ! *** Output ***
    zref        => mlcanopy_inst%zref_forcing        , &  ! Atmospheric reference height (m)
    ztop        => mlcanopy_inst%ztop_canopy         , &  ! Canopy foliage top height (m)
    zbot        => mlcanopy_inst%zbot_canopy         , &  ! Canopy foliage bottom height (m)
    ncan        => mlcanopy_inst%ncan_canopy         , &  ! Number of aboveground layers
    ntop        => mlcanopy_inst%ntop_canopy         , &  ! Index for top leaf layer
    nbot        => mlcanopy_inst%nbot_canopy         , &  ! Index for bottom leaf layer
    dlai_frac   => mlcanopy_inst%dlai_frac_profile   , &  ! Canopy layer leaf area index (fraction of canopy total)
    dsai_frac   => mlcanopy_inst%dsai_frac_profile   , &  ! Canopy layer stem area index (fraction of canopy total)
    zs          => mlcanopy_inst%zs_profile          , &  ! Canopy layer height for scalar concentration and source (m)
    zw          => mlcanopy_inst%zw_profile          , &  ! Canopy height at interface between two adjacent layers (m)
    dz          => mlcanopy_inst%dz_profile            &  ! Canopy layer thickness (m)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Atmospheric forcing height

       zref(p) = forc_hgt_u(p)

       ! Top height of canopy foliage

       ztop(p) = htop(p)

       ! Distance from canopy to zref

       ztop_to_zref = zref(p) - ztop(p)

       ! Determine the number of layers and their thickness using a
       ! user-specified number of layers or a user-specified thickness

       if (nlayer_within > 0 .and. nlayer_above > 0) then

          ! Use specifed number of layers within (nlayer_within) and above (nlayer_above) the canopy

          ntop(p) = nlayer_within
          dz_within = ztop(p) / float(ntop(p))
          nabove = nlayer_above
          ncan(p) = ntop(p) + nabove
          dz_above = ztop_to_zref / float(nabove)

       else if (nlayer_within == 0 .or. nlayer_above == 0) then

          ! Determine number of within-canopy layers by setting height increment
          ! based on tall or short canopies as specified by dz_param. Define ntop
          ! as the number of within-canopy layers and adjust dz_within if needed.

          if (ztop(p) > dz_param) then
             dz_within = dz_tall
          else
             dz_within = dz_short
          end if

          ntop(p) = nint(ztop(p) / dz_within)
          dz_within = ztop(p) / float(ntop(p))

          ! Now determine the number of above-canopy layers and their thickness

          dz_above = dz_within
          nabove = nint(ztop_to_zref / dz_above)

          ! Uncomment these two lines to run 1L at all sites with same above canopy layers as ML
          ! ntop(p) = 1
          ! dz_within = ztop(p) / float(ntop(p))

          ncan(p) = ntop(p) + nabove
          dz_above = ztop_to_zref / float(nabove)

       else
          call endrun (msg=' ERROR: initVerticalStructure: invalid canopy specification')
       end if

       if (ncan(p) > nlevmlcan) then
          call endrun (msg=' ERROR: initVerticalStructure: ncan > nlevmlcan')
       end if

       ! Calculate heights at layer interfaces (zw). They are defined for
       ! ic = 0 (ground) to ic = ncan (top above-canopy layer). For layer ic,
       ! zw(p,ic-1) is the bottom height and zw(p,ic) is the top height.
       ! First calculate within-canopy heights (for ic = 0 to ntop).

       ic = ntop(p)
       zw(p,ic) = ztop(p)
       do ic = ntop(p)-1, 0, -1
          zw(p,ic) = zw(p,ic+1) - dz_within
       end do

       if (abs(zw(p,0)) > 1.e-10_r8) then
          call endrun (msg=' ERROR: initVerticalStructure: zw(p,0) improperly defined')
       end if

       zw(p,0) = max(zw(p,0), 0._r8)   ! Set small negative value to zero
       zw(p,0) = min(zw(p,0), 0._r8)   ! Set small positive value to zero

       ! Now calculate the above-canopy heights (for ic = ntop+1 to ncan)

       ic = ncan(p)
       zw(p,ic) = zref(p)
       do ic = ncan(p)-1, ntop(p)+1, -1
          zw(p,ic) = zw(p,ic+1) - dz_above
       end do

       ! Thickness of each layer

       do ic = 1, ncan(p)
          dz(p,ic) = zw(p,ic) - zw(p,ic-1)
       end do

       ! Determine heights of the scalar concentration and scalar source
       ! (zs). These are physically centered in the layer.

       do ic = 1, ncan(p)
          zs(p,ic) = 0.5_r8 * (zw(p,ic) + zw(p,ic-1))
       end do

       ! Calculate leaf area index at each height from the leaf area density profile.
       ! Use the cumulative distribution function evaluated at the the bottom and top
       ! heights for the layer. The vertical profile of leaf area density is specified
       ! using the 2-parameter beta distribution.
       !
       ! See Bonan et al. (2018) Geosci. Model Dev., 11, 1467-1496, doi:10.5194/gmd-11-1467-2018, eq. (28)

       do ic = 1, ntop(p)

          ! Lower height at bottom of layer

          zrel = min(zw(p,ic-1)/ztop(p), 1._r8)
          beta_cdf_bot = beta_distribution_cdf(pbeta_lai(p,1),pbeta_lai(p,2),zrel)

          ! Upper height at top of layer

          zrel = min(zw(p,ic)/ztop(p), 1._r8)
          beta_cdf_top = beta_distribution_cdf(pbeta_lai(p,1),pbeta_lai(p,2),zrel)

          ! Leaf area index (m2/m2)

          dlai(p,ic) = (beta_cdf_top - beta_cdf_bot) * unit_lai

       end do

       ! Repeat these calculations for stem area index. This allows for a different
       ! profile of stem area index.

       do ic = 1, ntop(p)
          zrel = min(zw(p,ic-1)/ztop(p), 1._r8)
          beta_cdf_bot = beta_distribution_cdf(pbeta_sai(p,1),pbeta_sai(p,2),zrel)
          zrel = min(zw(p,ic)/ztop(p), 1._r8)
          beta_cdf_top = beta_distribution_cdf(pbeta_sai(p,1),pbeta_sai(p,2),zrel)
          dsai(p,ic) = (beta_cdf_top - beta_cdf_bot) * unit_sai
       end do

       ! Check to make sure sum of dlai+dsai matches canopy lai+sai

       pai_err = sum(dlai(p,1:ntop(p))) + sum(dsai(p,1:ntop(p)))
       if (abs(pai_err - (unit_lai+unit_sai)) > 1.e-06_r8) then
          call endrun (msg=' ERROR: initVerticalStructure: plant area profile does not sum to canopy total')
       end if

       ! Now determine the number of layers with plant area. Set the layers with
       ! plant area < dpai_min to zero and sum the "missing" leaf and stem area
       ! so that this can be distributed back across the layers.

       lai_miss = 0._r8 ; sai_miss = 0._r8
       do ic = 1, ntop(p)
          if ((dlai(p,ic)+dsai(p,ic)) < dpai_min) then
             lai_miss = lai_miss + dlai(p,ic)
             sai_miss = sai_miss + dsai(p,ic)
             dlai(p,ic) = 0._r8
             dsai(p,ic) = 0._r8
          end if
       end do

       ! Distribute the missing leaf area across layers in proportion to the leaf area profile

       if (lai_miss > 0._r8) then
          lai_err = sum(dlai(p,1:ntop(p)))
          do ic = 1, ntop(p)
             dlai(p,ic) = dlai(p,ic) + lai_miss * (dlai(p,ic) / lai_err)
          end do
       end if

       ! Now do the same for stem area

       if (sai_miss > 0._r8) then
          sai_err = sum(dsai(p,1:ntop(p)))
          do ic = 1, ntop(p)
             dsai(p,ic) = dsai(p,ic) + sai_miss * (dsai(p,ic) / sai_err)
          end do
       end if

       ! Find the lowest leaf/stem layer. nbot is the first layer above the
       ! ground with leaves or stems and is used to accommodate open layers
       ! beneath the bottom of the canopy. It is used mostly in radiative
       ! transfer, where the ground is referenced as ic=0 but the first
       ! leaf/stem layer may not be ic=1.

       nbot(p) = 0
       do ic = ntop(p), 1, -1
          if ((dlai(p,ic)+dsai(p,ic)) > 0._r8) nbot(p) = ic
       end do

       if (nbot(p) == 0) then
          call endrun (msg=' ERROR: initVerticalStructure: nbot not defined')
       end if

       ! Bottom height of canopy is at the bottom of layer nbot (equal to top of nbot-1)

       zbot(p) = zw(p,nbot(p)-1)

       ! Error check

       lai_err = sum(dlai(p,1:ntop(p)))
       if (abs(lai_err - unit_lai) > 1.e-06_r8) then
          call endrun (msg=' ERROR: initVerticalStructure: leaf area profile does not sum to canopy total after redistribution')
       end if

       sai_err = sum(dsai(p,1:ntop(p)))
       if (abs(sai_err - unit_sai) > 1.e-06_r8) then
          call endrun (msg=' ERROR: initVerticalStructure: stem area profile does not sum to canopy total after redistribution')
       end if

       ! Zero out above-canopy layers

       do ic = ntop(p)+1, ncan(p)
          dlai(p,ic) = 0._r8
          dsai(p,ic) = 0._r8
       end do

       ! Normalize profiles

       do ic = 1, ncan(p)
          dlai_frac(p,ic) = dlai(p,ic) / unit_lai
          dsai_frac(p,ic) = dsai(p,ic) / unit_sai
       end do

       ! Check to make sure all canopy layers have dpai > 0

       iflag = 0
       do ic = nbot(p), ntop(p)
          if ((dlai_frac(p,ic) + dsai_frac(p,ic)) <= 0._r8) iflag = 1
       end do
       if (iflag == 1) then
          call endrun (msg=' ERROR: initVerticalStructure: canopy layer has zero plant area index')
       end if

    end do

    end associate
  end subroutine initVerticalStructure

  !-----------------------------------------------------------------------
  subroutine initVerticalProfiles (num_filter, filter, &
  atm2lnd_inst, wateratm2lndbulk_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Initialize vertical profiles and canopy states
    !
    ! !USES:
    use PatchType, only : patch
    use MLclm_varcon, only : mmh2o, mmdry
    use MLclm_varpar, only : isun, isha
    use atm2lndType, only : atm2lnd_type
    use Wateratm2lndBulkType, only : wateratm2lndbulk_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter            ! Number of patches in filter
    integer, intent(in) :: filter(:)             ! Patch filter

    type(atm2lnd_type), intent(in) :: atm2lnd_inst
    type(wateratm2lndbulk_type), intent(in) :: wateratm2lndbulk_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                               ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                ! Column index for CLM g/l/c/p hierarchy
    integer  :: g                                ! Grid cell index for CLM g/l/c/p hierarchy
    integer  :: ic                               ! Aboveground layer index
    !---------------------------------------------------------------------

    associate ( &
                                                                   ! *** Input ***
    forc_u      => atm2lnd_inst%forc_u_grc                    , &  ! CLM: Atmospheric wind speed in east direction (m/s)
    forc_v      => atm2lnd_inst%forc_v_grc                    , &  ! CLM: Atmospheric wind speed in north direction (m/s)
    forc_pco2   => atm2lnd_inst%forc_pco2_grc                 , &  ! CLM: Atmospheric CO2 partial pressure (Pa)
    forc_t      => atm2lnd_inst%forc_t_downscaled_col         , &  ! CLM: Atmospheric temperature (K)
    forc_q      => wateratm2lndbulk_inst%forc_q_downscaled_col, &  ! CLM: Atmospheric specific humidity (kg/kg)
    forc_pbot   => atm2lnd_inst%forc_pbot_downscaled_col      , &  ! CLM: Atmospheric pressure (Pa)
    ncan        => mlcanopy_inst%ncan_canopy                  , &  ! Number of aboveground layers
                                                                   ! *** Output ***
    tg          => mlcanopy_inst%tg_soil                      , &  ! Soil surface temperature (K)
    wind        => mlcanopy_inst%wind_profile                 , &  ! Canopy layer wind speed (m/s)
    tair        => mlcanopy_inst%tair_profile                 , &  ! Canopy layer air temperature (K)
    eair        => mlcanopy_inst%eair_profile                 , &  ! Canopy layer vapor pressure (Pa)
    cair        => mlcanopy_inst%cair_profile                 , &  ! Canopy layer atmospheric CO2 (umol/mol)
    h2ocan      => mlcanopy_inst%h2ocan_profile               , &  ! Canopy layer intercepted water (kg H2O/m2)
    lwp         => mlcanopy_inst%lwp_leaf                     , &  ! Leaf water potential (MPa)
    tleaf       => mlcanopy_inst%tleaf_leaf                     &  ! Leaf temperature (K)
    )

    ! Initialization of vertical profiles and canopy states

    do fp = 1, num_filter
       p = filter(fp)
       c = patch%column(p)
       g = patch%gridcell(p)

       do ic = 1, ncan(p)
          wind(p,ic) = sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g))
          tair(p,ic) = forc_t(c)
          eair(p,ic) = forc_q(c) * forc_pbot(c) / (mmh2o/mmdry + (1._r8-mmh2o/mmdry) * forc_q(c))
          cair(p,ic) = forc_pco2(g) / forc_pbot(c) * 1.e06_r8

          tleaf(p,ic,isun) = forc_t(c) ; tleaf(p,ic,isha) = forc_t(c)
          lwp(p,ic,isun) = -0.1_r8 ; lwp(p,ic,isha) = -0.1_r8
          h2ocan(p,ic) = 0._r8
       end do

       tg(p) = forc_t(c)
    end do

    end associate
  end subroutine initVerticalProfiles

  !-----------------------------------------------------------------------
  subroutine getPADparameters (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Get parameters that define plant area density vertical profile. These
    ! should be read in from the CLM surface dataset but for convenience are
    ! set here based on PFT. The code is only executed if the parameters have
    ! not been previously set to user specified values.
    !
    ! !USES:
    use clm_varpar, only : mxpft
    use PatchType, only : patch
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter            ! Number of patches in filter
    integer, intent(in) :: filter(:)             ! Patch filter
    type(mlcanopy_type), intent(out) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                               ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    real(r8) :: pbeta_lai_pft(0:mxpft,2)         ! PFT parameters for the leaf area density beta distribution (-)
    real(r8) :: pbeta_sai_pft(0:mxpft,2)         ! PFT parameters for the stem area density beta distribution (-)
    !---------------------------------------------------------------------

    associate ( &
    pbeta_lai   => mlcanopy_inst%pbeta_lai_canopy , &  ! Parameters for the leaf area density 2-parameter beta distribution (-)
    pbeta_sai   => mlcanopy_inst%pbeta_sai_canopy   &  ! Parameters for the stem area density 2-parameter beta distribution (-)
    )

    ! Parameters for the leaf/stem area density beta distribution.
    ! NOTE: pbeta(1) = pbeta(2) = 1 gives a uniform distribution.

    ! Values for the CLM PFTs

    pbeta_lai_pft(:,:) = -999._r8
    pbeta_sai_pft(:,:) = -999._r8

    ! 1 => needleleaf_evergreen_temperate_tree
    pbeta_lai_pft(1,1) = 11.5_r8; pbeta_lai_pft(1,2) = 3.5_r8

    ! 2 => needleleaf_evergreen_boreal_tree
    ! 3 => needleleaf_deciduous_boreal_tree
    pbeta_lai_pft(2:3,1) = 3.5_r8; pbeta_lai_pft(2:3,2) = 2.0_r8

    ! 4 => broadleaf_evergreen_tropical_tree
    ! 5 => broadleaf_evergreen_temperate_tree
    pbeta_lai_pft(4:5,1) = 3.5_r8; pbeta_lai_pft(4:5,2) = 2.0_r8

    ! 6 => broadleaf_deciduous_tropical_tree
    ! 7 => broadleaf_deciduous_temperate_tree
    ! 8 => broadleaf_deciduous_boreal_tree
    pbeta_lai_pft(6:8,1) = 3.5_r8; pbeta_lai_pft(6:8,2) = 2.0_r8

    !  9 => broadleaf_evergreen_shrub
    ! 10 => broadleaf_deciduous_temperate_shrub
    ! 11 => broadleaf_deciduous_boreal_shrub
    pbeta_lai_pft(9:11,1) = 3.5_r8; pbeta_lai_pft(9:11,2) = 2.0_r8

    ! 12 => c3_arctic_grass
    ! 13 => c3_non-arctic_grass
    ! 14 => c4_grass
    ! 15 => c3_crop
    ! 16 => c3_irrigated
    pbeta_lai_pft(12:16,1) = 2.5_r8; pbeta_lai_pft(12:16,2) = 2.5_r8

    ! Stem area index
    pbeta_sai_pft(:,:) = pbeta_lai_pft(:,:)

    ! Assign values to each patch based on PFT
    !
    ! Note: pbeta_lai and pbeta_sai are assigned as -spval on initialization. If they
    ! have beeen set to valid values > 0, then those assigned values are used and not
    ! the PFT values. 

    do fp = 1, num_filter
       p = filter(fp)
       if (pbeta_lai(p,1) < 0._r8 .or. pbeta_lai(p,2) < 0._r8 .or. pbeta_sai(p,1) < 0._r8 .or. pbeta_sai(p,2) < 0._r8) then
          pbeta_lai(p,1) = pbeta_lai_pft(patch%itype(p),1)
          pbeta_lai(p,2) = pbeta_lai_pft(patch%itype(p),2)
          pbeta_sai(p,1) = pbeta_sai_pft(patch%itype(p),1)
          pbeta_sai(p,2) = pbeta_sai_pft(patch%itype(p),2)
       end if
    end do

    end associate
  end subroutine getPADparameters

end module MLinitVerticalMod
