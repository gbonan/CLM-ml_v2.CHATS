module MLCanopyFluxesMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate  multilayer canopy fluxes
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use abortutils             , only : endrun
  use clm_varctl             , only : iulog
  use decompMod              , only : bounds_type
  use atm2lndType            , only : atm2lnd_type
  use CanopyStateType        , only : canopystate_type
  use ColumnType             , only : col
  use EnergyFluxType         , only : energyflux_type
  use FrictionVelocityMod    , only : frictionvel_type
  use GridcellType           , only : grc
  use PatchType              , only : patch
  use SoilStateType          , only : soilstate_type
  use SolarAbsorbedType      , only : solarabs_type
  use SurfaceAlbedoType      , only : surfalb_type
  use TemperatureType        , only : temperature_type
  use Wateratm2lndBulkType   , only : wateratm2lndbulk_type
  use WaterDiagnosticBulkType, only : waterdiagnosticbulk_type
  use WaterFluxBulkType      , only : waterfluxbulk_type
  use WaterStateBulkType     , only : waterstatebulk_type
  use MLCanopyFluxesType     , only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PRIVATE TYPES:
  integer, parameter :: nvar1d = 23     ! Number of single-level fluxes to accumulate over multilayer canopy timesteps
  integer, parameter :: nvar2d = 14     ! Number of multi-level profile fluxes to accumulate over multilayer canopy timesteps
  integer, parameter :: nvar3d = 12     ! Number of multi-level leaf fluxes to accumulate over multilayer canopy timesteps
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: MLCanopyFluxes              ! Compute canopy fluxes
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: GetCLMVar                  ! Copy CLM variables to multilayer canopy variables
  private :: MLTimeStepFluxIntegration  ! Integrate fluxes over multilayer canopy timesteps to the CLM timestep
  private :: CanopyFluxesDiagnostics    ! Sum leaf and soil fluxes and other canopy diagnostics
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine MLCanopyFluxes (bounds, num_exposedvegp, filter_exposedvegp, &
  atm2lnd_inst, canopystate_inst, soilstate_inst, temperature_inst, waterstatebulk_inst, &
  waterfluxbulk_inst, energyflux_inst, frictionvel_inst, surfalb_inst, solarabs_inst, &
  mlcanopy_inst, wateratm2lndbulk_inst, waterdiagnosticbulk_inst)
    !
    ! !DESCRIPTION:
    ! Compute fluxes for sunlit and shaded leaves at each level
    ! and for soil surface
    !
    ! !USES:
    use clm_time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_varcon, only : grav, spval
    use clm_varpar, only : ivis, inir
    use spmdMod, only : masterproc
    use MLclm_varcon, only : mmh2o, rgas
    use MLclm_varctl, only : mlcan_to_clm, dtime_ml, ml_vert_init, runge_kutta_type, nrk, met_type
    use MLclm_varpar, only : isun, isha, nlevmlcan, nleaf
    use MLCanopyNitrogenProfileMod, only : CanopyNitrogenProfile
    use MLCanopyTurbulenceMod, only : CanopyTurbulence
    use MLFluxProfileSolutionMod, only : FluxProfileSolution
    use MLCanopyWaterMod, only : CanopyWettedFraction, CanopyInterception, CanopyEvaporation
    use MLinitVerticalMod, only : initVerticalProfiles, initVerticalStructure, getPADparameters
    use MLLeafBoundaryLayerMod, only : LeafBoundaryLayer
    use MLLeafHeatCapacityMod, only : LeafHeatCapacity
    use MLLeafPhotosynthesisMod, only : LeafPhotosynthesis
    use MLLongwaveRadiationMod, only : LongwaveRadiation
    use MLPlantHydraulicsMod, only: SoilResistance, PlantResistance, LeafWaterPotential
    use MLSolarRadiationMod, only: SolarRadiation
    use MLWaterVaporMod, only : LatVap
    use MLGetAtmForcingMod, only : GetAtmForcing
    use MLRungeKuttaMod, only : RungeKuttaIni, RungeKuttaUpdate
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_exposedvegp           ! Number of non-snow-covered veg patches in CLM patch filter
    integer, intent(in) :: filter_exposedvegp(:)     ! CLM patch filter for non-snow-covered vegetation

    type(atm2lnd_type)            , intent(in)    :: atm2lnd_inst
    type(canopystate_type)        , intent(inout) :: canopystate_inst
    type(soilstate_type)          , intent(inout) :: soilstate_inst
    type(temperature_type)        , intent(inout) :: temperature_inst
    type(waterstatebulk_type)     , intent(inout) :: waterstatebulk_inst
    type(waterfluxbulk_type)      , intent(inout) :: waterfluxbulk_inst
    type(energyflux_type)         , intent(inout) :: energyflux_inst
    type(frictionvel_type)        , intent(inout) :: frictionvel_inst
    type(surfalb_type)            , intent(inout) :: surfalb_inst
    type(solarabs_type)           , intent(inout) :: solarabs_inst
    type(mlcanopy_type)           , intent(inout) :: mlcanopy_inst
    type(wateratm2lndbulk_type)   , intent(in)    :: wateratm2lndbulk_inst
    type(waterdiagnosticbulk_type), intent(inout) :: waterdiagnosticbulk_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: num_mlcan                               ! Number of vegetated patches for multilayer canopy
    integer  :: filter_mlcan(bounds%endp-bounds%begp+1) ! Patch filter for multilayer canopy
    integer  :: fp                                      ! Filter index
    integer  :: p                                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                       ! Column index for CLM g/l/c/p hierarchy
    integer  :: g                                       ! Gridcell index for CLM g/l/c/p hierarchy
    integer  :: ic                                      ! Aboveground layer index
    integer  :: nstep                                   ! Current CLM timestep number
    integer  :: num_ml_steps                            ! Number of multilayer canopy timesteps within a CLM timestep
    integer  :: nstep_ml                                ! Current multilayer timestep number
    real(r8) :: dtime_clm                               ! CLM timestep (s)
    real(r8) :: curr_calday_end                         ! Current calendar day at end of CLM timestep (equals 1.000 on 0Z January 1 of current year)
    real(r8) :: curr_calday_beg                         ! Current calendar day at beginning of CLM timestep
    real(r8) :: calday_interp_bef                       ! Calendar day for previous atmospheric forcing (used in interpolation to ML timestep)
    real(r8) :: calday_interp_cur                       ! Calendar day for current atmospheric forcing (used in interpolation to ML timestep)
    real(r8) :: calday_interp_next                      ! Calendar day for next atmospheric forcing (used in interpolation to ML timestep)
    real(r8) :: calday_interp_ml                        ! Calendar day for current ML timestep
    real(r8) :: totpai                                  ! Canopy lai+sai for error check (m2/m2)

    ! These are used to accumulate flux variables over the multilayer canopy timesteps to the CLM timestep.
    ! The last dimension is the number of variables

    real(r8) :: flux_accumulator(bounds%begp:bounds%endp,nvar1d)                          ! Single-level fluxes
    real(r8) :: flux_accumulator_profile(bounds%begp:bounds%endp,1:nlevmlcan+1,nvar2d)    ! Multi-level profile fluxes
    real(r8) :: flux_accumulator_leaf(bounds%begp:bounds%endp,1:nlevmlcan,1:nleaf,nvar3d) ! Multi-level leaf fluxes

    integer  :: irk                                     ! Runge-Kutta step
    integer  :: nrk_steps                               ! Number of Runge-Kutta steps
    real(r8) :: ark(nrk,nrk), brk(nrk), crk(nrk)        ! Runge-Kutta parameters
    !---------------------------------------------------------------------

    associate ( &
                                                                         ! *** CLM variables ***
    elai           => canopystate_inst%elai_patch                   , &  ! Leaf area index of canopy (m2/m2)
    esai           => canopystate_inst%esai_patch                   , &  ! Stem area index of canopy (m2/m2)
    smp_l          => soilstate_inst%smp_l_col                      , &  ! Soil layer matric potential (mm)
    t_soisno       => temperature_inst%t_soisno_col                 , &  ! Soil temperature (K)
    eflx_lh_tot    => energyflux_inst%eflx_lh_tot_patch             , &  ! OUTPUT TO CLM: patch total latent heat flux (W/m2)
    eflx_sh_tot    => energyflux_inst%eflx_sh_tot_patch             , &  ! OUTPUT TO CLM: patch total sensible heat flux (W/m2)
    eflx_lwrad_out => energyflux_inst%eflx_lwrad_out_patch          , &  ! OUTPUT TO CLM: patch emitted infrared (longwave) radiation (W/m2)
    taux           => energyflux_inst%taux_patch                    , &  ! OUTPUT TO CLM: patch wind (shear) stress: e-w (kg/m/s2)
    tauy           => energyflux_inst%tauy_patch                    , &  ! OUTPUT TO CLM: patch wind (shear) stress: n-s (kg/m/s2)
    fv             => frictionvel_inst%fv_patch                     , &  ! OUTPUT TO CLM: patch friction velocity (m/s)
    u10_clm        => frictionvel_inst%u10_clm_patch                , &  ! OUTPUT TO CLM: patch 10-m wind (m/s)
    fsa            => solarabs_inst%fsa_patch                       , &  ! OUTPUT TO CLM: patch solar radiation absorbed (total) (W/m2)
    albd           => surfalb_inst%albd_patch                       , &  ! OUTPUT TO CLM: patch surface albedo (direct)
    albi           => surfalb_inst%albi_patch                       , &  ! OUTPUT TO CLM: patch surface albedo (diffuse)
    t_ref2m        => temperature_inst%t_ref2m_patch                , &  ! OUTPUT TO CLM: patch 2 m height surface air temperature (K)
    qflx_evap_tot  => waterfluxbulk_inst%qflx_evap_tot_patch        , &  ! OUTPUT TO CLM: patch total evapotranspiration flux (kg H2O/m2/s)
    q_ref2m        => waterdiagnosticbulk_inst%q_ref2m_patch        , &  ! OUTPUT TO CLM: patch 2 m height surface specific humidity (kg/kg)

                                                                  ! *** Multilayer canopy variables ***
    zref           => mlcanopy_inst%zref_forcing             , &  ! Atmospheric reference height (m)
    tref_bef       => mlcanopy_inst%tref_bef_forcing         , &  ! Air temperature at reference height (K) [previous CLM timestep]
    tref_cur       => mlcanopy_inst%tref_cur_forcing         , &  ! Air temperature at reference height (K) [current CLM timestep]
    qref_bef       => mlcanopy_inst%qref_bef_forcing         , &  ! Specific humidity at reference height (kg/kg) [previous CLM timestep]
    qref_cur       => mlcanopy_inst%qref_cur_forcing         , &  ! Specific humidity at reference height (kg/kg) [current CLM timestep]
    uref_bef       => mlcanopy_inst%uref_bef_forcing         , &  ! Wind speed at reference height (m/s) [previous CLM timestep]
    uref_cur       => mlcanopy_inst%uref_cur_forcing         , &  ! Wind speed at reference height (m/s) [current CLM timestep]
    pref_bef       => mlcanopy_inst%pref_bef_forcing         , &  ! Air pressure at reference height (Pa) [previous CLM timestep]
    pref_cur       => mlcanopy_inst%pref_cur_forcing         , &  ! Air pressure at reference height (Pa) [current CLM timestep]
    co2ref_bef     => mlcanopy_inst%co2ref_bef_forcing       , &  ! Atmospheric CO2 at reference height (umol/mol) [previous CLM timestep]
    co2ref_cur     => mlcanopy_inst%co2ref_cur_forcing       , &  ! Atmospheric CO2 at reference height (umol/mol) [current CLM timestep]
    swskyb_bef     => mlcanopy_inst%swskyb_bef_forcing       , &  ! Atmospheric direct beam solar radiation (W/m2) [previous CLM timestep]
    swskyb_cur     => mlcanopy_inst%swskyb_cur_forcing       , &  ! Atmospheric direct beam solar radiation (W/m2) [current CLM timestep]
    swskyd_bef     => mlcanopy_inst%swskyd_bef_forcing       , &  ! Atmospheric diffuse solar radiation (W/m2) [previous CLM timestep]
    swskyd_cur     => mlcanopy_inst%swskyd_cur_forcing       , &  ! Atmospheric diffuse solar radiation (W/m2) [current CLM timestep]
    lwsky_bef      => mlcanopy_inst%lwsky_bef_forcing        , &  ! Atmospheric longwave radiation (W/m2) [previous CLM timestep]
    lwsky_cur      => mlcanopy_inst%lwsky_cur_forcing        , &  ! Atmospheric longwave radiation (W/m2) [current CLM timestep]
    ncan           => mlcanopy_inst%ncan_canopy              , &  ! Number of aboveground layers
    lai            => mlcanopy_inst%lai_canopy               , &  ! Leaf area index of canopy (m2/m2)
    sai            => mlcanopy_inst%sai_canopy               , &  ! Stem area index of canopy (m2/m2)
    swveg          => mlcanopy_inst%swveg_canopy             , &  ! Absorbed solar radiation: vegetation (W/m2)
    lwup           => mlcanopy_inst%lwup_canopy              , &  ! Upward longwave radiation above canopy (W/m2)
    shflx          => mlcanopy_inst%shflx_canopy             , &  ! Total sensible heat flux, including soil (W/m2)
    lhflx          => mlcanopy_inst%lhflx_canopy             , &  ! Total latent heat flux, including soil (W/m2)
    etflx          => mlcanopy_inst%etflx_canopy             , &  ! Total water vapor flux, including soil (mol H2O/m2/s)
    ustar          => mlcanopy_inst%ustar_canopy             , &  ! Friction velocity (m/s)
    swsoi          => mlcanopy_inst%swsoi_soil               , &  ! Absorbed solar radiation: ground (W/m2)
    lwsoi          => mlcanopy_inst%lwsoi_soil               , &  ! Absorbed longwave radiation: ground (W/m2)
    rnsoi          => mlcanopy_inst%rnsoi_soil               , &  ! Net radiation: ground (W/m2)
    tg             => mlcanopy_inst%tg_soil                  , &  ! Soil surface temperature (K)
    tg_bef         => mlcanopy_inst%tg_bef_soil              , &  ! Soil surface temperature for previous timestep (K)
    rhg            => mlcanopy_inst%rhg_soil                 , &  ! Relative humidity of airspace at soil surface (fraction)
    dlai           => mlcanopy_inst%dlai_profile             , &  ! Canopy layer leaf area index (m2/m2)
    dsai           => mlcanopy_inst%dsai_profile             , &  ! Canopy layer stem area index (m2/m2)
    dpai           => mlcanopy_inst%dpai_profile             , &  ! Canopy layer plant area index (m2/m2)
    dlai_frac      => mlcanopy_inst%dlai_frac_profile        , &  ! Canopy layer leaf area index (fraction of canopy total)
    dsai_frac      => mlcanopy_inst%dsai_frac_profile        , &  ! Canopy layer stem area index (fraction of canopy total)
    fracsun        => mlcanopy_inst%fracsun_profile          , &  ! Canopy layer sunlit fraction (-)
    tair           => mlcanopy_inst%tair_profile             , &  ! Canopy layer air temperature (K)
    tair_bef       => mlcanopy_inst%tair_bef_profile         , &  ! Canopy layer air temperature for previous timestep (K)
    eair           => mlcanopy_inst%eair_profile             , &  ! Canopy layer vapor pressure (Pa)
    eair_bef       => mlcanopy_inst%eair_bef_profile         , &  ! Canopy layer vapor pressure for previous timestep (Pa)
    cair           => mlcanopy_inst%cair_profile             , &  ! Canopy layer atmospheric CO2 (umol/mol)
    cair_bef       => mlcanopy_inst%cair_bef_profile         , &  ! Canopy layer atmospheric CO2 for previous timestep (umol/mol)
    h2ocan         => mlcanopy_inst%h2ocan_profile           , &  ! Canopy layer intercepted water (kg H2O/m2)
    h2ocan_bef     => mlcanopy_inst%h2ocan_bef_profile       , &  ! Canopy layer intercepted water for previous timestep (kg H2O/m2)
    swleaf         => mlcanopy_inst%swleaf_leaf              , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    lwleaf         => mlcanopy_inst%lwleaf_leaf              , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf         => mlcanopy_inst%rnleaf_leaf              , &  ! Leaf net radiation (W/m2 leaf)
    tleaf          => mlcanopy_inst%tleaf_leaf               , &  ! Leaf temperature (K)
    tleaf_bef      => mlcanopy_inst%tleaf_bef_leaf           , &  ! Leaf temperature for previous timestep (K)
    tleaf_hist     => mlcanopy_inst%tleaf_hist_leaf          , &  ! Leaf temperature (not sun/shade average) for history files (K)
    lwp            => mlcanopy_inst%lwp_leaf                 , &  ! Leaf water potential (MPa)
    lwp_bef        => mlcanopy_inst%lwp_bef_leaf             , &  ! Leaf water potential for previous timestep (MPa)
    lwp_hist       => mlcanopy_inst%lwp_hist_leaf              &  ! Leaf water potential (not sun/shade average) for history files (MPa)
    )

    ! Get current step counter (nstep) and step size (dtime_clm) from CLM

    nstep = get_nstep()
    dtime_clm = get_step_size()

    ! Get current calendar day at end of CLM timestep

    curr_calday_end = get_curr_calday(offset=0)

    ! Calendar day at beginning of CLM timestep

    curr_calday_beg = get_curr_calday(offset=-int(dtime_clm))

    ! Set calendar days used in interpolation of atmospheric forcing
    ! from CLM timestep to ML timestep
    !
    ! ======================
    ! 2-point interpolation: 
    ! ======================
    ! Uses atmospheric forcing at current CLM timestep and previous CLM timestep.
    ! The corresponding calendar day is at the end of the timestep (as in CLM).
    ! This option is needed to couple with CAM or to couple with CLM.
    !
    ! ======================
    ! 3-point interpolation: 
    ! ======================
    ! Uses atmospheric forcing at current, previous, and next CLM timesteps.
    ! The corresponding calendar day is at the middle of the timestep (as in CHATS).
    ! This option is valid for tower forcing that is centered in the timestep.

    select case (met_type)
    case (0)
       ! No interpolation
       calday_interp_cur = 0._r8
       calday_interp_bef = 0._r8
       calday_interp_next = 0._r8
    case (2)
       ! 2-point interpolation
       call endrun (msg=' ERROR: met_type not valid')
       calday_interp_cur = curr_calday_end
       calday_interp_bef = calday_interp_cur - dtime_clm / 86400._r8
       calday_interp_next = 0._r8
    case (3)
       ! 3-point interpolation
       calday_interp_cur = 0.5_r8 * (curr_calday_end + curr_calday_beg)
       calday_interp_bef = calday_interp_cur - dtime_clm / 86400._r8
       calday_interp_next = calday_interp_cur + dtime_clm / 86400._r8
    end select

    ! Number of multilayer canopy timesteps within a CLM timestep

    num_ml_steps = int(dtime_clm / dtime_ml)

    ! Build filter of patches to process with multilayer canopy.
    ! Use this filter to exclude some CLM patches if so desired
    ! (e.g., exclude short vegetation).

    num_mlcan = 0
    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       g = patch%gridcell(fp)
       num_mlcan = num_mlcan + 1
       filter_mlcan(num_mlcan) = p
    end do

    ! Initialize multilayer canopy model. The call for initialization is
    ! triggered by any zref = spval, which only happens on the first CLM
    ! timestep.

    ml_vert_init = 0
    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       if (zref(p) == spval) ml_vert_init = 1
    end do

    if (ml_vert_init == 1) then

       ! Read parameters for leaf area density and stem area density vertical profiles.
       ! These should come from the CLM surface dataset, but are set here instead.

       if (masterproc) then
          write (iulog,*) 'Attempting to initialize multilayer canopy pad parameters .....'
       end if

       call getPADparameters (num_mlcan, filter_mlcan, mlcanopy_inst)

       if (masterproc) then
          write (iulog,*) 'Successfuly initialized multilayer canopy pad parameters'
       end if

       ! Initialize multilayer canopy vertical structure and profiles. This is only
       ! done once (on first CLM timestep) because the forcing height (and therefore
       ! the volume of air in the surface layer) changes between CLM timesteps.

       if (masterproc) then
          write (iulog,*) 'Attempting to initialize multilayer canopy vertical structure .....'
       end if

       call initVerticalStructure (bounds, num_mlcan, filter_mlcan, &
       canopystate_inst, frictionvel_inst, mlcanopy_inst)

       call initVerticalProfiles (num_mlcan, filter_mlcan, &
       atm2lnd_inst, wateratm2lndbulk_inst, mlcanopy_inst)

       if (masterproc) then
          write (iulog,*) 'Successfully initialized multilayer canopy vertical structure'
       end if

       ! Initialize Runge-Kutta parameters

       call RungeKuttaIni (ark, brk, crk)

    end if

    ! Copy CLM variables to multilayer canopy variables

     call GetCLMVar (nstep, dtime_clm, num_mlcan, filter_mlcan, atm2lnd_inst, & 
     soilstate_inst, temperature_inst, surfalb_inst, wateratm2lndbulk_inst, &
     mlcanopy_inst)

    ! Initialize atmospheric forcing for previous CLM timestep. This is
    ! needed for interpolation to the multilayer canopy timestep.

    if (ml_vert_init == 1) then
       do fp = 1, num_mlcan
          p = filter_mlcan(fp)
          uref_bef(p) = uref_cur(p)
          tref_bef(p) = tref_cur(p)
          qref_bef(p) = qref_cur(p)
          pref_bef(p) = pref_cur(p)
          co2ref_bef(p) = co2ref_cur(p)
          swskyb_bef(p,ivis) = swskyb_cur(p,ivis)
          swskyb_bef(p,inir) = swskyb_cur(p,inir)
          swskyd_bef(p,ivis) = swskyd_cur(p,ivis)
          swskyd_bef(p,inir) = swskyd_cur(p,inir)
          lwsky_bef(p) = lwsky_cur(p)
       end do
    end if

    ! Update leaf and stem area profile for current values

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)

       ! Get canopy values from CLM

       lai(p) = elai(p)
       sai(p) = esai(p)

       ! Scale vertical profiles for current lai and sai. This assumes that
       ! the plant area density profile does not change with time, only the 
       ! total lai and sai change

       do ic = 1, ncan(p)
          dlai(p,ic) = dlai_frac(p,ic) * lai(p)
          dsai(p,ic) = dsai_frac(p,ic) * sai(p)

          ! Reset values to minimum: Needed for CLM5/6 US-NR1 testing
!         if (dlai(p,ic) > 0._r8) dlai(p,ic) = max(dlai(p,ic), 0.01_r8)
!         if (dsai(p,ic) > 0._r8) dsai(p,ic) = max(dsai(p,ic), 0.01_r8)

          dpai(p,ic) = dlai(p,ic) + dsai(p,ic)
       end do

       ! Comment out this error check if using minimum dlai and dsai

       totpai = sum(dpai(p,1:ncan(p)))
       if (abs(totpai - (lai(p)+sai(p))) > 1.e-06_r8) then
          call endrun (msg=' ERROR: MLCanopyFluxes: plant area index not updated correctly')
       end if

    end do

    ! Terms for plant hydraulics: soil resistance and leaf-specific conductance

    call SoilResistance (num_mlcan, filter_mlcan, &
    soilstate_inst, waterstatebulk_inst, mlcanopy_inst)

    call PlantResistance (num_mlcan, filter_mlcan, mlcanopy_inst)

    ! Leaf heat capacity

    call LeafHeatCapacity (num_mlcan, filter_mlcan, mlcanopy_inst)

    ! Relative humidity in soil airspace using CLM smp_l and t_soisno

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       c = patch%column(p)
       rhg(p) = exp(grav * mmh2o * smp_l(c,1)*1.e-03_r8 / (rgas * t_soisno(c,1)))
    end do

    ! Calculate fluxes (timestep = dtime_ml) over the CLM timestep (dtime_clm)

    do nstep_ml = 1, num_ml_steps

       ! Calendar day for multilayer canopy timestep

       select case (met_type)
       case (0, 2)
          ! Calendar is at end of time interval
          calday_interp_ml = curr_calday_beg + float(nstep_ml) * (dtime_ml / 86400._r8)
       case (3)
          ! Calendar is centered in time interval
          calday_interp_ml = curr_calday_beg + (float(nstep_ml) - 0.5_r8) * (dtime_ml / 86400._r8)
       end select

       ! Save values for previous multilayer canopy timestep

       do fp = 1, num_mlcan
          p = filter_mlcan(fp)
          tg_bef(p) = tg(p)
          do ic = 1, ncan(p)
             tair_bef(p,ic) = tair(p,ic)
             eair_bef(p,ic) = eair(p,ic)
             cair_bef(p,ic) = cair(p,ic)
             h2ocan_bef(p,ic) = h2ocan(p,ic)
             tleaf_bef(p,ic,isun) = tleaf(p,ic,isun)
             tleaf_bef(p,ic,isha) = tleaf(p,ic,isha)
             lwp_bef(p,ic,isun) = lwp(p,ic,isun)
             lwp_bef(p,ic,isha) = lwp(p,ic,isha)
          end do
       end do

       ! Atmospheric forcing for current mutilayer canopy timestep

       call GetAtmForcing (calday_interp_bef, calday_interp_cur, calday_interp_next, &
       calday_interp_ml, num_mlcan, filter_mlcan, mlcanopy_inst)

       ! Solar radiation transfer through the canopy

       call SolarRadiation (bounds, num_mlcan, filter_mlcan, mlcanopy_inst)

       ! Canopy profile of photosynthetic capacity (uses fracsun)

       call CanopyNitrogenProfile (num_mlcan, filter_mlcan, mlcanopy_inst)

       ! Problem: Fluxes depend on state variables but state variables depend on fluxes.
       ! Solution: Use Runge-Kutta methods to solve flux-profile equations and dependencies
       ! with intermediate estimates for state variables (tair, eair, tleaf, tg, lwp, h2ocan).

       ! case 10: Euler (1st-order) with no intermediate estimates for state variables.
       ! Advances through a full timestep using state values from previous timestep.

       ! case 20s, 30s, and 40s: 2nd-, 3rd-, and 4th-order Runge-Kutta methods.
       ! Use n = 2, 3, or 4 steps to obtain n intermediate estimates for state variables.
       ! Use one more step with weighted estimates for state variables to advance through
       ! a full timestep. Runge-Kutta methods are defined in RungeKuttaIni.

       select case (runge_kutta_type)
       case (10)
          nrk_steps = 0
       case (20 :)
          nrk_steps = runge_kutta_type / 10
       end select

       do irk = 1, nrk_steps+1

          ! Canopy wetted fraction
          !
          ! Dependency: Wetted (fwet) and dry green (fdry) fractions of canopy depend
          ! on intercepted water (h2ocan), which changes over the timestep

          call CanopyWettedFraction (num_mlcan, filter_mlcan, mlcanopy_inst)

          ! Longwave radiation transfer through the canopy
          !
          ! Dependency: Longwave flux (lwleaf, lwsoi) depends on leaf (tleaf)
          ! and ground (tg) temperatures, which depend on longwave flux

          call LongwaveRadiation (bounds, num_mlcan, filter_mlcan, mlcanopy_inst)

          ! Net radiation at each layer and at ground

          do fp = 1, num_mlcan
             p = filter_mlcan(fp)
             do ic = 1, ncan(p)
                rnleaf(p,ic,isun) = swleaf(p,ic,isun,ivis) + swleaf(p,ic,isun,inir) + lwleaf(p,ic,isun)
                rnleaf(p,ic,isha) = swleaf(p,ic,isha,ivis) + swleaf(p,ic,isha,inir) + lwleaf(p,ic,isha)
             end do
             rnsoi(p) = swsoi(p,ivis) + swsoi(p,inir) + lwsoi(p)
          end do

          ! Canopy turbulence, wind speed, and aerodynamic conductances
          !
          ! Dependency: 
          !
          ! 1. Obukhov length (obu) depends on canopy temperature (tair(ntop) = taf)
          ! and canopy specific humidity (eair(ntop) -> qaf)
          !
          ! 2. Turbulence quantities depend on Obukhov length:
          !    ustar - Friction velocity
          !    beta  - Value of u* / u at canopy top
          !    PrSc  - Prandtl (Schmidt) number at canopy top
          !    zdisp - Displacement height (m)
          !    wind  - Wind speed
          !    gac   - Aerodynamic conductances

          call CanopyTurbulence (nstep_ml, num_mlcan, filter_mlcan, mlcanopy_inst)

          ! Leaf boundary layer conductance
          !
          ! Dependency: Conductances depend on wind speed (wind), which changes
          ! over the timestep based on Obukhov length (i.e., on fluxes). Also leaf
          ! temperature (tleaf) and air temperature (tair) are needed for Grashof number

          call LeafBoundaryLayer (num_mlcan, filter_mlcan, isun, mlcanopy_inst)
          call LeafBoundaryLayer (num_mlcan, filter_mlcan, isha, mlcanopy_inst)

          ! Photosynthesis and stomatal conductance
          !
          ! Dependency: Stomatal conductance (gs) depends on leaf temperature
          ! (tleaf), VPD (eair, tleaf), and leaf water potential (lwp), which
          ! depend on fluxes

          call LeafPhotosynthesis (num_mlcan, filter_mlcan, isun, mlcanopy_inst)
          call LeafPhotosynthesis (num_mlcan, filter_mlcan, isha, mlcanopy_inst)

          ! Source/sink fluxes for leaves and soil, and concentration profiles.
          ! Solves for tair, eair, tleaf, and tg and associated fluxes
          !
          ! Dependency: Fluxes and concentration profiles depend on net radiation,
          ! leaf boundary layer conductance, leaf stomatal conductance, wetted and
          ! dry green fraction, and aerodynamic conductance calculated in prior
          ! subroutine calls
          !
          ! Note: See associated MLSoilTemperatureMod, which uses the calculated soil
          ! heat flux to update soil temperature. This is not compatible with CLM's
          ! soil temperature.

          call FluxProfileSolution (num_mlcan, filter_mlcan, mlcanopy_inst)

          ! Update leaf water potential for the current transpiration rate
          !
          ! Dependency: Leaf water potential (lwp) depends on transpiration (trleaf),
          ! which depends on leaf water potential (through stomatal conductance)

          call LeafWaterPotential (num_mlcan, filter_mlcan, isun, mlcanopy_inst)
          call LeafWaterPotential (num_mlcan, filter_mlcan, isha, mlcanopy_inst)

          ! Canopy intercepted water
          !
          ! Dependency: Canopy water (h2ocan) depends on canopy evaporation
          ! (evleaf), which depends on wetted fraction of canopy. Dew depends
          ! on both evaporation (evleaf) and transpiraion (trleaf).

          call CanopyInterception (num_mlcan, filter_mlcan, mlcanopy_inst)
          call CanopyEvaporation (num_mlcan, filter_mlcan, mlcanopy_inst)

          ! Update states for next Runge-Kutta step (but not for Euler)

          if (nrk_steps > 0 .and. irk <= nrk_steps) then
             call RungeKuttaUpdate (irk, ark, brk, crk, num_mlcan, filter_mlcan, mlcanopy_inst)
          end if

       end do

       ! Write 5-min output (same format as flux.out file written in CLMml_driver)

       go to 100
       write (36,'(f12.7,17f10.3)') calday_interp_ml, &
       (mlcanopy_inst%swbeam_profile(p,mlcanopy_inst%ntop_canopy(p),ivis) + &          ! net radiation (W/m2)
       mlcanopy_inst%swbeam_profile(p,mlcanopy_inst%ntop_canopy(p),inir) + &
       mlcanopy_inst%swdwn_profile(p,mlcanopy_inst%ntop_canopy(p),ivis) + &
       mlcanopy_inst%swdwn_profile(p,mlcanopy_inst%ntop_canopy(p),inir) + &
       mlcanopy_inst%lwdwn_profile(p,mlcanopy_inst%ntop_canopy(p))) - &
       (mlcanopy_inst%swupw_profile(p,mlcanopy_inst%ntop_canopy(p),ivis) + &
       mlcanopy_inst%swupw_profile(p,mlcanopy_inst%ntop_canopy(p),inir) + &
       mlcanopy_inst%lwupw_profile(p,mlcanopy_inst%ntop_canopy(p))), &
       sum(mlcanopy_inst%stair_profile(p,1:mlcanopy_inst%ncan_canopy(p))), &           ! canopy air storage flux (W/m2)
       mlcanopy_inst%shair_profile(p,mlcanopy_inst%ncan_canopy(p)), &                  ! sensible heat flux (W/m2)
       mlcanopy_inst%etair_profile(p,mlcanopy_inst%ncan_canopy(p))*LatVap(mlcanopy_inst%tref_forcing(p)), &  ! latent heat flux (W/m2)
       sum(mlcanopy_inst%agross_leaf(p,1:mlcanopy_inst%ntop_canopy(p),isun)*fracsun(p,1:mlcanopy_inst%ntop_canopy(p))* &  ! gpp
       dpai(p,1:mlcanopy_inst%ntop_canopy(p))*mlcanopy_inst%fdry_profile(p,1:mlcanopy_inst%ntop_canopy(p))/ &
       (1._r8-mlcanopy_inst%fwet_profile(p,1:mlcanopy_inst%ntop_canopy(p))))+ &
       sum(mlcanopy_inst%agross_leaf(p,1:mlcanopy_inst%ntop_canopy(p),isha)*(1._r8-fracsun(p,1:mlcanopy_inst%ntop_canopy(p)))* &
       dpai(p,1:mlcanopy_inst%ntop_canopy(p))*mlcanopy_inst%fdry_profile(p,1:mlcanopy_inst%ntop_canopy(p))/ &
       (1._r8-mlcanopy_inst%fwet_profile(p,1:mlcanopy_inst%ntop_canopy(p)))), &
       mlcanopy_inst%ustar_canopy(p), &                                                ! friction velocity (m/s)
       (mlcanopy_inst%swupw_profile(p,mlcanopy_inst%ntop_canopy(p),ivis) + &           ! reflected solar radiation (W/m2)
       mlcanopy_inst%swupw_profile(p,mlcanopy_inst%ntop_canopy(p),inir)), &
       lwup(p), &                                                                      ! upward longwave radiation (W/m2)
       mlcanopy_inst%tair_profile(p,mlcanopy_inst%ntop_canopy(p)), &                   ! air temperature at canopy top (K)
       mlcanopy_inst%gsoi_soil(p), mlcanopy_inst%rnsoi_soil(p), mlcanopy_inst%shsoi_soil(p), mlcanopy_inst%lhsoi_soil(p), &
       (sum(mlcanopy_inst%trleaf_leaf(p,1:mlcanopy_inst%ntop_canopy(p),isun)*&          ! transpiration (W/m2)
       fracsun(p,1:mlcanopy_inst%ntop_canopy(p))*dpai(p,1:mlcanopy_inst%ntop_canopy(p))) + &
       sum(mlcanopy_inst%trleaf_leaf(p,1:mlcanopy_inst%ntop_canopy(p),isha)*&
       (1._r8-fracsun(p,1:mlcanopy_inst%ntop_canopy(p)))*dpai(p,1:mlcanopy_inst%ntop_canopy(p))))*LatVap(mlcanopy_inst%tref_forcing(p)), &
       (sum(mlcanopy_inst%evleaf_leaf(p,1:mlcanopy_inst%ntop_canopy(p),isun)*&          ! canopy evaporation (W/m2)
       fracsun(p,1:mlcanopy_inst%ntop_canopy(p))*dpai(p,1:mlcanopy_inst%ntop_canopy(p))) + &
       sum(mlcanopy_inst%evleaf_leaf(p,1:mlcanopy_inst%ntop_canopy(p),isha)*&
       (1._r8-fracsun(p,1:mlcanopy_inst%ntop_canopy(p)))*dpai(p,1:mlcanopy_inst%ntop_canopy(p))))*LatVap(mlcanopy_inst%tref_forcing(p)), &
       mlcanopy_inst%beta_canopy(p), &                                                 ! beta (-)
       sum(mlcanopy_inst%stleaf_leaf(p,1:mlcanopy_inst%ntop_canopy(p),isun)* &         ! vegetation storage flux (W/m2)
       fracsun(p,1:mlcanopy_inst%ntop_canopy(p))*dpai(p,1:mlcanopy_inst%ntop_canopy(p))) + &
       sum(mlcanopy_inst%stleaf_leaf(p,1:mlcanopy_inst%ntop_canopy(p),isha)* &
       (1._r8-fracsun(p,1:mlcanopy_inst%ntop_canopy(p)))*dpai(p,1:mlcanopy_inst%ntop_canopy(p)))
100    continue

       ! Fluxes need to be accumulated to the CLM timestep over all multilayer canopy
       ! timesteps. Other variables are instantaneous for the final multilayer timestep.

       call MLTimeStepFluxIntegration (nstep_ml, num_ml_steps, num_mlcan, filter_mlcan, &
       flux_accumulator, flux_accumulator_profile, flux_accumulator_leaf, mlcanopy_inst)

    end do    ! End multilayer canopy time-stepping

    ! Save current CLM atmospheric forcing to use in ML timestep interpolation at next CLM timestep

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       uref_bef(p) = uref_cur(p)
       tref_bef(p) = tref_cur(p)
       qref_bef(p) = qref_cur(p)
       pref_bef(p) = pref_cur(p)
       co2ref_bef(p) = co2ref_cur(p)
       swskyb_bef(p,ivis) = swskyb_cur(p,ivis)
       swskyb_bef(p,inir) = swskyb_cur(p,inir)
       swskyd_bef(p,ivis) = swskyd_cur(p,ivis)
       swskyd_bef(p,inir) = swskyd_cur(p,inir)
       lwsky_bef(p) = lwsky_cur(p)
    end do

    ! Sum leaf and soil fluxes and other canopy diagnostics

    call CanopyFluxesDiagnostics (num_mlcan, filter_mlcan, mlcanopy_inst)

    ! Leaf temperature and leaf water potential are prognostic variables
    ! for sunlit and shaded leaves. But sun/shade fractions change over
    ! time. Merge temperature and leaf water potential for sunlit and
    ! shaded leaves to layer-average value, which is used at next CLM timestep.
    ! Use this method (rather than retaining sun/shade states) because some
    ! shade leaf becomes sun leaf (and vice versa) between CLM timesteps as
    ! fracsun changes.

    do fp = 1, num_mlcan
       p = filter_mlcan(fp)
       do ic = 1, ncan(p)

          ! First save sun/shade leaves for model output

          tleaf_hist(p,ic,isun) = tleaf(p,ic,isun)
          tleaf_hist(p,ic,isha) = tleaf(p,ic,isha)
          lwp_hist(p,ic,isun) = lwp(p,ic,isun)
          lwp_hist(p,ic,isha) = lwp(p,ic,isha)

          ! Now merge sun/shade leaves

          if (dpai(p,ic) > 0._r8) then
             tleaf(p,ic,isun) = tleaf(p,ic,isun) * fracsun(p,ic) + tleaf(p,ic,isha) * (1._r8 - fracsun(p,ic))
             tleaf(p,ic,isha) = tleaf(p,ic,isun)
             lwp(p,ic,isun) = lwp(p,ic,isun) * fracsun(p,ic) + lwp(p,ic,isha) * (1._r8 - fracsun(p,ic))
             lwp(p,ic,isha) = lwp(p,ic,isun)
          end if

       end do
    end do

    ! Copy multilayer canopy variables to CLM variables. These are
    ! passed from CLM to CAM. 

    if (mlcan_to_clm == 1) then
       do fp = 1, num_mlcan
          p = filter_mlcan(fp)
          albd(p,ivis) = 0._r8 ; albd(p,inir) = 0._r8
          albi(p,ivis) = 0._r8 ; albi(p,inir) = 0._r8
          taux(p) = 0._r8
          tauy(p) = 0._r8
          eflx_lh_tot(p) = lhflx(p)
          eflx_sh_tot(p) = shflx(p)
          eflx_lwrad_out(p) = lwup(p)
          qflx_evap_tot(p) = etflx(p) * mmh2o
          fv(p) = ustar(p)
          u10_clm(p) = 0._r8
          t_ref2m(p) = 0._r8
          q_ref2m(p) = 0._r8
          fsa(p) = swveg(p,ivis) + swveg(p,inir) + swsoi(p,ivis) + swsoi(p,inir)
       end do
    end if

    end associate
  end subroutine MLCanopyFluxes

  !-----------------------------------------------------------------------
  subroutine GetCLMVar (nstep, dtime_clm, num_filter, filter, atm2lnd_inst, & 
  soilstate_inst, temperature_inst, surfalb_inst, wateratm2lndbulk_inst, &
  mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Copy CLM variables to multilayer canopy variables
    !
    ! !USES:
    use clm_varpar, only : ivis, inir
    use clm_varcon, only : pi => rpi
    use clm_time_manager, only : get_curr_calday
    use clm_varorb, only : eccen, obliqr, lambm0, mvelpp
    use shr_orb_mod, only : shr_orb_decl, shr_orb_cosz
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: nstep         ! Current CLM timestep number
    real(r8), intent(in) :: dtime_clm     ! CLM timestep (s)
    integer, intent(in)  :: num_filter    ! Number of patches in filter
    integer, intent(in)  :: filter(:)     ! Patch filter

    type(atm2lnd_type)         , intent(in)    :: atm2lnd_inst
    type(soilstate_type)       , intent(in)    :: soilstate_inst
    type(temperature_type)     , intent(in)    :: temperature_inst
    type(surfalb_type)         , intent(in)    :: surfalb_inst
    type(wateratm2lndbulk_type), intent(in)    :: wateratm2lndbulk_inst
    type(mlcanopy_type)        , intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                       ! Column index for CLM g/l/c/p hierarchy
    integer  :: g                       ! Gridcell index for CLM g/l/c/p hierarchy
    real(r8) :: lat, lon                ! Latitude and longitude (radians)
    real(r8) :: coszen                  ! Cosine solar zenith angle
    real(r8) :: caldaym1                ! Calendar day for zenith angle (1.000 on 0Z January 1 of current year)
    real(r8) :: declinm1                ! Solar declination angle for zenith angle (radians)
    real(r8) :: eccf                    ! Earth orbit eccentricity factor
    !---------------------------------------------------------------------

    associate ( &
                                                                         ! *** Input ***
    forc_u         => atm2lnd_inst%forc_u_grc                       , &  ! CLM: Atmospheric wind speed in east direction (m/s)
    forc_v         => atm2lnd_inst%forc_v_grc                       , &  ! CLM: Atmospheric wind speed in north direction (m/s)
    forc_pco2      => atm2lnd_inst%forc_pco2_grc                    , &  ! CLM: Atmospheric CO2 partial pressure (Pa)
    forc_po2       => atm2lnd_inst%forc_po2_grc                     , &  ! CLM: Atmospheric O2 partial pressure (Pa)
    forc_solad_col => atm2lnd_inst%forc_solad_downscaled_col        , &  ! CLM: Atmospheric direct beam radiation (W/m2)
    forc_solai     => atm2lnd_inst%forc_solai_grc                   , &  ! CLM: Atmospheric diffuse radiation (W/m2)
    forc_t         => atm2lnd_inst%forc_t_downscaled_col            , &  ! CLM: Atmospheric temperature (K)
    forc_pbot      => atm2lnd_inst%forc_pbot_downscaled_col         , &  ! CLM: Atmospheric pressure (Pa)
    forc_lwrad     => atm2lnd_inst%forc_lwrad_downscaled_col        , &  ! CLM: Atmospheric longwave radiation (W/m2)
    forc_q         => wateratm2lndbulk_inst%forc_q_downscaled_col   , &  ! CLM: Atmospheric specific humidity (kg/kg)
    forc_rain      => wateratm2lndbulk_inst%forc_rain_downscaled_col, &  ! CLM: Rainfall rate (mm/s)
    forc_snow      => wateratm2lndbulk_inst%forc_snow_downscaled_col, &  ! CLM: Snowfall rate (mm/s)
    albgrd         => surfalb_inst%albgrd_col                       , &  ! CLM: Direct beam albedo of ground (soil)
    albgri         => surfalb_inst%albgri_col                       , &  ! CLM: Diffuse albedo of ground (soil)
    soilresis      => soilstate_inst%soilresis_col                  , &  ! CLM: Soil evaporative resistance (s/m)
    thk            => soilstate_inst%thk_col                        , &  ! CLM: Soil layer thermal conductivity (W/m/K)
    t_a10_patch    => temperature_inst%t_a10_patch                  , &  ! CLM: 10-day running mean of the 2-m temperature (K)
    t_soisno       => temperature_inst%t_soisno_col                 , &  ! CLM: Soil temperature (K)
    snl            => col%snl                                       , &  ! CLM: Number of snow layers
    z              => col%z                                         , &  ! CLM: Soil layer depth (m)
    zi             => col%zi                                        , &  ! CLM: Soil layer depth at layer interface (m)
                                                                  ! *** Output ***
    tref_cur       => mlcanopy_inst%tref_cur_forcing         , &  ! Air temperature at reference height (K) [current CLM timestep]
    qref_cur       => mlcanopy_inst%qref_cur_forcing         , &  ! Specific humidity at reference height (kg/kg) [current CLM timestep]
    uref_cur       => mlcanopy_inst%uref_cur_forcing         , &  ! Wind speed at reference height (m/s) [current CLM timestep]
    pref_cur       => mlcanopy_inst%pref_cur_forcing         , &  ! Air pressure at reference height (Pa) [current CLM timestep]
    co2ref_cur     => mlcanopy_inst%co2ref_cur_forcing       , &  ! Atmospheric CO2 at reference height (umol/mol) [current CLM timestep]
    o2ref          => mlcanopy_inst%o2ref_forcing            , &  ! Atmospheric O2 at reference height (mmol/mol)
    solar_zen      => mlcanopy_inst%solar_zen_forcing        , &  ! Solar zenith angle (radians)
    swskyb_cur     => mlcanopy_inst%swskyb_cur_forcing       , &  ! Atmospheric direct beam solar radiation (W/m2) [current CLM timestep]
    swskyd_cur     => mlcanopy_inst%swskyd_cur_forcing       , &  ! Atmospheric diffuse solar radiation (W/m2) [current CLM timestep]
    lwsky_cur      => mlcanopy_inst%lwsky_cur_forcing        , &  ! Atmospheric longwave radiation (W/m2) [current CLM timestep]
    qflx_rain      => mlcanopy_inst%qflx_rain_forcing        , &  ! Rainfall (mm H2O/s = kg H2O/m2/s)
    qflx_snow      => mlcanopy_inst%qflx_snow_forcing        , &  ! Snowfall (mm H2O/s = kg H2O/m2/s)
    tacclim        => mlcanopy_inst%tacclim_forcing          , &  ! Average air temperature for acclimation (K)
    albsoib        => mlcanopy_inst%albsoib_soil             , &  ! Direct beam albedo of ground (-)
    albsoid        => mlcanopy_inst%albsoid_soil             , &  ! Diffuse albedo of ground (-)
    soilres        => mlcanopy_inst%soilres_soil             , &  ! Soil evaporative resistance (s/m)
    soil_t         => mlcanopy_inst%soil_t_soil              , &  ! Temperature of first snow/soil layer (K)
    soil_dz        => mlcanopy_inst%soil_dz_soil             , &  ! Depth to temperature of first snow/soil layer (m)
    soil_tk        => mlcanopy_inst%soil_tk_soil               &  ! Thermal conductivity of first snow/soil layer (W/m/K)
    )

    ! Copy CLM variables to multilayer canopy variables. Note the
    ! distinction between grid cell (g), column (c), and patch (p)
    ! variables. All multilayer canopy variables are for patches.

    do fp = 1, num_filter
       p = filter(fp)
       c = patch%column(p)
       g = patch%gridcell(p)

       ! Atmospheric forcing: CLM grid cell (g) variables -> patch (p) variables

       uref_cur(p) = sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g))
       swskyd_cur(p,ivis) = forc_solai(g,ivis)
       swskyd_cur(p,inir) = forc_solai(g,inir)

       ! Atmospheric forcing: CLM column (c) variables -> patch (p) variables

       tref_cur(p) = forc_t(c)
       qref_cur(p) = forc_q(c)
       pref_cur(p) = forc_pbot(c)
       lwsky_cur(p) = forc_lwrad(c)
       qflx_rain(p) = forc_rain(c)
       qflx_snow(p) = forc_snow(c)
       swskyb_cur(p,ivis) = forc_solad_col(c,ivis)
       swskyb_cur(p,inir) = forc_solad_col(c,inir)

       ! CO2 and O2: CLM grid cell (g) -> patch (p). Note that the units
       ! conversion requires pbot.

       co2ref_cur(p) = forc_pco2(g) / forc_pbot(c) * 1.e06_r8  ! Pa -> umol/mol
       o2ref(p) = forc_po2(g) / forc_pbot(c) * 1.e03_r8        ! Pa -> mmol/mol

       ! Miscellaneous

       tacclim(p) = t_a10_patch(p)

       ! Ground (soil) albedos: CLM column (c) -> patch (p)

       albsoib(p,ivis) = albgrd(c,ivis) ; albsoib(p,inir) = albgrd(c,inir)
       albsoid(p,ivis) = albgri(c,ivis) ; albsoid(p,inir) = albgri(c,inir)

       ! Soil evaporative resistance: CLM column (c) -> patch (p)

       soilres(p) = soilresis(c)

       ! Properties of first soil layer needed for soil heat flux:
       ! CLM column (c) -> patch (p)

       !!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!
       ! These variables are required to calculate the heat flux   !
       ! into the soil (gsoi_soil) as part of the canopy and soil  !
       ! energy balance. When run in stand-alone mode uncoupled    !
       ! from CLM, soil temperature (t_soisno) is updated with     !
       ! gsoi_soil as the boundary condition using:                !
       ! MLSoilTemperatureMod. When coupled to CLM, the CLM soil   !
       ! temperature (SoilTemperatureMod) DOES NOT use gsoi_soil   !
       ! and therefore t_soisno is calculated INDEPENDENT of the   !
       ! multilayer canopy soil energy balance.                    !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       soil_t(p) = t_soisno(c,snl(c)+1)          ! Temperature of first snow/soil layer (K)
       soil_dz(p) = (z(c,snl(c)+1)-zi(c,snl(c))) ! Depth to temperature of first snow/soil layer (m)
       soil_tk(p) = thk(c,snl(c)+1)              ! Thermal conductivity of first snow/soil layer (W/m/K)

    end do

    ! Solar zenith angle. Need to subtract one timestep (-dtime_clm) because
    ! zenith angle is calculated for the beginning of the timestep. So use
    ! calendar day at beginning of the timestep (nstep-1).

    caldaym1 = get_curr_calday(offset=-int(dtime_clm))
    call shr_orb_decl (caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf)

    do fp = 1, num_filter
       p = filter(fp)
       c = patch%column(p)
       g = patch%gridcell(p)

       lat = grc%latdeg(g) * pi / 180._r8
       lon = grc%londeg(g) * pi / 180._r8
       coszen = shr_orb_cosz (caldaym1, lat, lon, declinm1)
       solar_zen(p) = acos(max(0.01_r8,coszen))
    end do

    end associate
  end subroutine GetCLMVar

  !-----------------------------------------------------------------------
  subroutine MLTimeStepFluxIntegration (nstep_ml, num_ml_steps, num_filter, filter, &
  flux_accumulator, flux_accumulator_profile, flux_accumulator_leaf, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Integrate fluxes over multilayer canopy timesteps to the CLM timestep
    !
    ! !USES:
    use clm_varpar, only : ivis, inir
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nstep_ml                             ! Current multilayer timestep number
    integer, intent(in) :: num_ml_steps                         ! Number of multilayer canopy timesteps within a CLM timestep
    integer, intent(in) :: num_filter                           ! Number of patches in filter
    integer, intent(in) :: filter(:)                            ! Patch filter
    real(r8), intent(inout) :: flux_accumulator(:,:)            ! Single-level flux accumulator variable
    real(r8), intent(inout) :: flux_accumulator_profile(:,:,:)  ! Multi-level profile flux accumulator variable
    real(r8), intent(inout) :: flux_accumulator_leaf(:,:,:,:)   ! Multi-level leaf flux accumulator variable
    type(mlcanopy_type), intent(in) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                                              ! Filter index
    integer  :: p                                               ! Patch index for CLM g/l/c/p hierarchy
    integer  :: i,j,k                                           ! Variable index
    !---------------------------------------------------------------------

    associate ( &
    swskyb      => mlcanopy_inst%swskyb_forcing       , &  ! Atmospheric direct beam solar radiation (W/m2) [interpolated to ML timestep]
    swskyd      => mlcanopy_inst%swskyd_forcing       , &  ! Atmospheric diffuse solar radiation (W/m2) [interpolated to ML timestep]
    lwsky       => mlcanopy_inst%lwsky_forcing        , &  ! Atmospheric longwave radiation (W/m2) [interpolated to ML timestep]
    ncan        => mlcanopy_inst%ncan_canopy          , &  ! Number of aboveground layers
    ustar       => mlcanopy_inst%ustar_canopy         , &  ! Friction velocity (m/s)
    beta        => mlcanopy_inst%beta_canopy          , &  ! Value of u* / u at canopy top (-)
    obu         => mlcanopy_inst%obu_canopy           , &  ! Obukhov length (m)
    z0m         => mlcanopy_inst%z0m_canopy           , &  ! Roughness length for momentum (m)
    zdisp       => mlcanopy_inst%zdisp_canopy         , &  ! Displacement height (m)
    lwup        => mlcanopy_inst%lwup_canopy          , &  ! Upward longwave radiation above canopy (W/m2)
    qflx_intr   => mlcanopy_inst%qflx_intr_canopy     , &  ! Intercepted precipitation (kg H2O/m2/s)
    qflx_tflrain => mlcanopy_inst%qflx_tflrain_canopy , &  ! Total rain throughfall onto ground (kg H2O/m2/s)
    qflx_tflsnow => mlcanopy_inst%qflx_tflsnow_canopy , &  ! Total snow throughfall onto ground (kg H2O/m2/s)
    swsoi       => mlcanopy_inst%swsoi_soil           , &  ! Absorbed solar radiation: ground (W/m2)
    lwsoi       => mlcanopy_inst%lwsoi_soil           , &  ! Absorbed longwave radiation: ground (W/m2)
    rnsoi       => mlcanopy_inst%rnsoi_soil           , &  ! Net radiation: ground (W/m2)
    shsoi       => mlcanopy_inst%shsoi_soil           , &  ! Sensible heat flux: ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi_soil           , &  ! Latent heat flux: ground (W/m2)
    etsoi       => mlcanopy_inst%etsoi_soil           , &  ! Water vapor flux: ground (mol H2O/m2/s)
    gsoi        => mlcanopy_inst%gsoi_soil            , &  ! Soil heat flux (W/m2)
    gac0        => mlcanopy_inst%gac0_soil            , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    shair       => mlcanopy_inst%shair_profile        , &  ! Canopy layer air sensible heat flux (W/m2)
    etair       => mlcanopy_inst%etair_profile        , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair       => mlcanopy_inst%stair_profile        , &  ! Canopy layer air storage heat flux (W/m2)
    mflx        => mlcanopy_inst%mflx_profile         , &  ! Canopy layer momentum flux (m2/s2)
    kc_eddy     => mlcanopy_inst%kc_eddy_profile      , &  ! Canopy layer eddy diffusivity from Harman and Finnigan (m2/s)
    gac         => mlcanopy_inst%gac_profile          , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    swupw       => mlcanopy_inst%swupw_profile        , &  ! Upward diffuse solar flux above canopy layer (W/m2)
    swdwn       => mlcanopy_inst%swdwn_profile        , &  ! Downward diffuse solar flux above canopy layer (W/m2)
    swbeam      => mlcanopy_inst%swbeam_profile       , &  ! Direct beam solar flux above canopy layer (W/m2)
    lwupw       => mlcanopy_inst%lwupw_profile        , &  ! Upward longwave flux above canopy layer (W/m2)
    lwdwn       => mlcanopy_inst%lwdwn_profile        , &  ! Downward longwave flux above canopy layer (W/m2)
    swleaf      => mlcanopy_inst%swleaf_leaf          , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    lwleaf      => mlcanopy_inst%lwleaf_leaf          , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf      => mlcanopy_inst%rnleaf_leaf          , &  ! Leaf net radiation (W/m2 leaf)
    shleaf      => mlcanopy_inst%shleaf_leaf          , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf      => mlcanopy_inst%lhleaf_leaf          , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf      => mlcanopy_inst%trleaf_leaf          , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf      => mlcanopy_inst%evleaf_leaf          , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    stleaf      => mlcanopy_inst%stleaf_leaf          , &  ! Leaf storage heat flux (W/m2 leaf)
    anet        => mlcanopy_inst%anet_leaf            , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    agross      => mlcanopy_inst%agross_leaf          , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    gs          => mlcanopy_inst%gs_leaf                &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Initialize flux variables that are summed over ML timesteps

       if (nstep_ml == 1) then
          flux_accumulator(p,:) = 0._r8
          flux_accumulator_profile(p,:,:) = 0._r8
          flux_accumulator_leaf(p,:,:,:) = 0._r8
       end if

       ! Accumulate fluxes over ML timesteps for temporal averaging
       ! NOTE: State variables are not averaged over ML timesteps

       i = 0
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + ustar(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + beta(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + obu(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + z0m(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + zdisp(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lwup(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + swsoi(p,ivis)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + swsoi(p,inir)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lwsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + rnsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + shsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lhsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + etsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + gsoi(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + gac0(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + qflx_intr(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + qflx_tflrain(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + qflx_tflsnow(p)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + swskyb(p,ivis)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + swskyb(p,inir)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + swskyd(p,ivis)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + swskyd(p,inir)
       i = i + 1; flux_accumulator(p,i) = flux_accumulator(p,i) + lwsky(p)

       j = 0
       j = j + 1; flux_accumulator_profile(p,1:ncan(p),j) = flux_accumulator_profile(p,1:ncan(p),j) + shair(p,1:ncan(p))
       j = j + 1; flux_accumulator_profile(p,1:ncan(p),j) = flux_accumulator_profile(p,1:ncan(p),j) + etair(p,1:ncan(p))
       j = j + 1; flux_accumulator_profile(p,1:ncan(p),j) = flux_accumulator_profile(p,1:ncan(p),j) + stair(p,1:ncan(p))
       j = j + 1; flux_accumulator_profile(p,1:ncan(p),j) = flux_accumulator_profile(p,1:ncan(p),j) + mflx(p,1:ncan(p))
       j = j + 1; flux_accumulator_profile(p,1:ncan(p),j) = flux_accumulator_profile(p,1:ncan(p),j) + kc_eddy(p,1:ncan(p))
       j = j + 1; flux_accumulator_profile(p,1:ncan(p),j) = flux_accumulator_profile(p,1:ncan(p),j) + gac(p,1:ncan(p))

       ! Note special indexing for radiative fluxes, which are dimensioned 0 to ncan
       j = j + 1; flux_accumulator_profile(p,1:ncan(p)+1,j) = flux_accumulator_profile(p,1:ncan(p)+1,j) + swupw(p,0:ncan(p),ivis)
       j = j + 1; flux_accumulator_profile(p,1:ncan(p)+1,j) = flux_accumulator_profile(p,1:ncan(p)+1,j) + swupw(p,0:ncan(p),inir)
       j = j + 1; flux_accumulator_profile(p,1:ncan(p)+1,j) = flux_accumulator_profile(p,1:ncan(p)+1,j) + swdwn(p,0:ncan(p),ivis)
       j = j + 1; flux_accumulator_profile(p,1:ncan(p)+1,j) = flux_accumulator_profile(p,1:ncan(p)+1,j) + swdwn(p,0:ncan(p),inir)
       j = j + 1; flux_accumulator_profile(p,1:ncan(p)+1,j) = flux_accumulator_profile(p,1:ncan(p)+1,j) + swbeam(p,0:ncan(p),ivis)
       j = j + 1; flux_accumulator_profile(p,1:ncan(p)+1,j) = flux_accumulator_profile(p,1:ncan(p)+1,j) + swbeam(p,0:ncan(p),inir)
       j = j + 1; flux_accumulator_profile(p,1:ncan(p)+1,j) = flux_accumulator_profile(p,1:ncan(p)+1,j) + lwupw(p,0:ncan(p))
       j = j + 1; flux_accumulator_profile(p,1:ncan(p)+1,j) = flux_accumulator_profile(p,1:ncan(p)+1,j) + lwdwn(p,0:ncan(p))

       k = 0
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + swleaf(p,:,:,ivis)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + swleaf(p,:,:,inir)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + lwleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + rnleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + shleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + lhleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + trleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + evleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + stleaf(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + anet(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + agross(p,:,:)
       k = k + 1; flux_accumulator_leaf(p,:,:,k) = flux_accumulator_leaf(p,:,:,k) + gs(p,:,:)

       if (i > nvar1d .or. j > nvar2d .or. k > nvar3d) then
          call endrun (msg=' ERROR: MLTimeStepFluxIntegration: nvar error')
       end if

       ! Average fluxes over ML timesteps

       if (nstep_ml == num_ml_steps) then

          ! Time averaging

          flux_accumulator(p,:) = flux_accumulator(p,:) / float(num_ml_steps)
          flux_accumulator_profile(p,:,:) = flux_accumulator_profile(p,:,:) / float(num_ml_steps)
          flux_accumulator_leaf(p,:,:,:) = flux_accumulator_leaf(p,:,:,:) / float(num_ml_steps)

          ! Map fluxes to variables: variables must be in the same order as above

          i = 0
          i = i + 1; ustar(p) = flux_accumulator(p,i)
          i = i + 1; beta(p) = flux_accumulator(p,i)
          i = i + 1; obu(p) = flux_accumulator(p,i)
          i = i + 1; z0m(p) = flux_accumulator(p,i)
          i = i + 1; zdisp(p) = flux_accumulator(p,i)
          i = i + 1; lwup(p) = flux_accumulator(p,i)
          i = i + 1; swsoi(p,ivis) = flux_accumulator(p,i)
          i = i + 1; swsoi(p,inir) = flux_accumulator(p,i)
          i = i + 1; lwsoi(p) = flux_accumulator(p,i)
          i = i + 1; rnsoi(p) = flux_accumulator(p,i)
          i = i + 1; shsoi(p) = flux_accumulator(p,i)
          i = i + 1; lhsoi(p) = flux_accumulator(p,i)
          i = i + 1; etsoi(p) = flux_accumulator(p,i)
          i = i + 1; gsoi(p) = flux_accumulator(p,i)
          i = i + 1; gac0(p) = flux_accumulator(p,i)
          i = i + 1; qflx_intr(p) = flux_accumulator(p,i)
          i = i + 1; qflx_tflrain(p) = flux_accumulator(p,i)
          i = i + 1; qflx_tflsnow(p) = flux_accumulator(p,i)
          i = i + 1; swskyb(p,ivis) = flux_accumulator(p,i)
          i = i + 1; swskyb(p,inir) = flux_accumulator(p,i)
          i = i + 1; swskyd(p,ivis) = flux_accumulator(p,i)
          i = i + 1; swskyd(p,inir) = flux_accumulator(p,i)
          i = i + 1; lwsky(p) = flux_accumulator(p,i)

          j = 0
          j = j + 1; shair(p,1:ncan(p)) = flux_accumulator_profile(p,1:ncan(p),j)
          j = j + 1; etair(p,1:ncan(p)) = flux_accumulator_profile(p,1:ncan(p),j)
          j = j + 1; stair(p,1:ncan(p)) = flux_accumulator_profile(p,1:ncan(p),j)
          j = j + 1; mflx(p,1:ncan(p)) = flux_accumulator_profile(p,1:ncan(p),j)
          j = j + 1; kc_eddy(p,1:ncan(p)) = flux_accumulator_profile(p,1:ncan(p),j)
          j = j + 1; gac(p,1:ncan(p)) = flux_accumulator_profile(p,1:ncan(p),j)

          ! Note special indexing for radiative fluxes, which are dimensioned 0 to ncan
          j = j + 1; swupw(p,0:ncan(p),ivis) = flux_accumulator_profile(p,1:ncan(p)+1,j)
          j = j + 1; swupw(p,0:ncan(p),inir) = flux_accumulator_profile(p,1:ncan(p)+1,j)
          j = j + 1; swdwn(p,0:ncan(p),ivis) = flux_accumulator_profile(p,1:ncan(p)+1,j)
          j = j + 1; swdwn(p,0:ncan(p),inir) = flux_accumulator_profile(p,1:ncan(p)+1,j)
          j = j + 1; swbeam(p,0:ncan(p),ivis) = flux_accumulator_profile(p,1:ncan(p)+1,j)
          j = j + 1; swbeam(p,0:ncan(p),inir) = flux_accumulator_profile(p,1:ncan(p)+1,j)
          j = j + 1; lwupw(p,0:ncan(p)) = flux_accumulator_profile(p,1:ncan(p)+1,j)
          j = j + 1; lwdwn(p,0:ncan(p)) = flux_accumulator_profile(p,1:ncan(p)+1,j)

          k = 0
          k = k + 1; swleaf(p,:,:,ivis) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; swleaf(p,:,:,inir) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; lwleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; rnleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; shleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; lhleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; trleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; evleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; stleaf(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; anet(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; agross(p,:,:) = flux_accumulator_leaf(p,:,:,k)
          k = k + 1; gs(p,:,:) = flux_accumulator_leaf(p,:,:,k)

          if (i > nvar1d .or. j > nvar2d .or. k > nvar3d) then
             call endrun (msg=' ERROR: MLTimeStepFluxIntegration: nvar error')
          end if

       end if

    end do

    end associate
  end subroutine MLTimeStepFluxIntegration

  !-----------------------------------------------------------------------
  subroutine CanopyFluxesDiagnostics (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Sum leaf and soil fluxes to get canopy fluxes and calculate
    ! other canopy diagnostics
    !
    ! !USES:
    use clm_varpar, only : numrad, ivis, inir
    use MLclm_varctl, only : flux_profile_type
    use MLclm_varpar, only : isun, isha
    use MLWaterVaporMod, only : LatVap
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                      ! Aboveground layer index
    integer  :: ib                      ! Waveband index
    real(r8) :: err                     ! Energy imbalance (W/m2)
    real(r8) :: radin                   ! Incoming radiation (W/m2)
    real(r8) :: radout                  ! Outgoing radiation (W/m2)
    real(r8) :: avail                   ! Available energy (W/m2)
    real(r8) :: flux                    ! Turbulent fluxes + storage (W/m2)
    real(r8) :: fracgreen               ! Green fraction of plant area index: lai/(lai+sai)
    real(r8) :: minlwp                  ! Minimum leaf water potential for canopy water stress diagnostic (MPa)
    !---------------------------------------------------------------------

    associate ( &
                                                         ! *** Input ***
    tref        => mlcanopy_inst%tref_forcing       , &  ! Air temperature at reference height (K)
    swskyb      => mlcanopy_inst%swskyb_forcing     , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd_forcing     , &  ! Atmospheric diffuse solar radiation (W/m2)
    lwsky       => mlcanopy_inst%lwsky_forcing      , &  ! Atmospheric longwave radiation (W/m2)
    ncan        => mlcanopy_inst%ncan_canopy        , &  ! Number of aboveground layers
    ntop        => mlcanopy_inst%ntop_canopy        , &  ! Index for top leaf layer
    lai         => mlcanopy_inst%lai_canopy         , &  ! Leaf area index of canopy (m2/m2)
    sai         => mlcanopy_inst%sai_canopy         , &  ! Stem area index of canopy (m2/m2)
!   swveg       => mlcanopy_inst%swveg_canopy       , &  ! Absorbed solar radiation: vegetation (W/m2)
!   albcan      => mlcanopy_inst%albcan_canopy      , &  ! Albedo above canopy (-)
    lwup        => mlcanopy_inst%lwup_canopy        , &  ! Upward longwave radiation above canopy (W/m2)
    shsoi       => mlcanopy_inst%shsoi_soil         , &  ! Sensible heat flux: ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi_soil         , &  ! Latent heat flux: ground (W/m2)
    gsoi        => mlcanopy_inst%gsoi_soil          , &  ! Soil heat flux (W/m2)
    swsoi       => mlcanopy_inst%swsoi_soil         , &  ! Absorbed solar radiation: ground (W/m2)
    lwsoi       => mlcanopy_inst%lwsoi_soil         , &  ! Absorbed longwave radiation: ground (W/m2)
    etsoi       => mlcanopy_inst%etsoi_soil         , &  ! Water vapor flux: ground (mol H2O/m2/s)
    dpai        => mlcanopy_inst%dpai_profile       , &  ! Canopy layer plant area index (m2/m2)
    fwet        => mlcanopy_inst%fwet_profile       , &  ! Canopy layer fraction of plant area index that is wet
    fdry        => mlcanopy_inst%fdry_profile       , &  ! Canopy layer fraction of plant area index that is green and dry
    tair        => mlcanopy_inst%tair_profile       , &  ! Canopy layer air temperature (K)
    wind        => mlcanopy_inst%wind_profile       , &  ! Canopy layer wind speed (m/s)
    shair       => mlcanopy_inst%shair_profile      , &  ! Canopy layer air sensible heat flux (W/m2)
    etair       => mlcanopy_inst%etair_profile      , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair       => mlcanopy_inst%stair_profile      , &  ! Canopy layer air storage heat flux (W/m2)
    lwupw       => mlcanopy_inst%lwupw_profile      , &  ! Upward longwave flux above canopy layer (W/m2)
    lwdwn       => mlcanopy_inst%lwdwn_profile      , &  ! Downward longwave flux above canopy layer (W/m2)
    swupw       => mlcanopy_inst%swupw_profile      , &  ! Upward diffuse solar flux above canopy layer (W/m2)
    swdwn       => mlcanopy_inst%swdwn_profile      , &  ! Downward diffuse solar flux above canopy layer (W/m2)
    swbeam      => mlcanopy_inst%swbeam_profile     , &  ! Direct beam solar flux above canopy layer (W/m2)
    fracsun     => mlcanopy_inst%fracsun_profile    , &  ! Canopy layer sunlit fraction (-)
    vcmax25_profile => mlcanopy_inst%vcmax25_profile, &  ! Canopy layer leaf maximum carboxylation rate at 25C (umol/m2/s)
    lwleaf      => mlcanopy_inst%lwleaf_leaf        , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf      => mlcanopy_inst%rnleaf_leaf        , &  ! Leaf net radiation (W/m2 leaf)
    stleaf      => mlcanopy_inst%stleaf_leaf        , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf      => mlcanopy_inst%shleaf_leaf        , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf      => mlcanopy_inst%lhleaf_leaf        , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf      => mlcanopy_inst%trleaf_leaf        , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf      => mlcanopy_inst%evleaf_leaf        , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    swleaf      => mlcanopy_inst%swleaf_leaf        , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    agross      => mlcanopy_inst%agross_leaf        , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    apar        => mlcanopy_inst%apar_leaf          , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    anet        => mlcanopy_inst%anet_leaf          , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    gs          => mlcanopy_inst%gs_leaf            , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    tleaf       => mlcanopy_inst%tleaf_leaf         , &  ! Leaf temperature (K)
    lwp         => mlcanopy_inst%lwp_leaf           , &  ! Leaf water potential (MPa)
    vcmax25_leaf=> mlcanopy_inst%vcmax25_leaf       , &  ! Leaf maximum carboxylation rate at 25C (umol/m2/s)
                                                         ! *** Output ***
    rnet        => mlcanopy_inst%rnet_canopy        , &  ! Total net radiation, including soil (W/m2)
    stflx_air   => mlcanopy_inst%stflx_air_canopy   , &  ! Canopy air storage heat flux (W/m2)
    stflx_veg   => mlcanopy_inst%stflx_veg_canopy   , &  ! Canopy biomass storage heat flux (W/m2)
    shflx       => mlcanopy_inst%shflx_canopy       , &  ! Total sensible heat flux, including soil (W/m2)
    lhflx       => mlcanopy_inst%lhflx_canopy       , &  ! Total latent heat flux, including soil (W/m2)
    etflx       => mlcanopy_inst%etflx_canopy       , &  ! Total water vapor flux, including soil (mol H2O/m2/s)
    albcan      => mlcanopy_inst%albcan_canopy      , &  ! Albedo above canopy (-)
    swveg       => mlcanopy_inst%swveg_canopy       , &  ! Absorbed solar radiation: vegetation (W/m2)
    swvegsun    => mlcanopy_inst%swvegsun_canopy    , &  ! Absorbed solar radiation: sunlit canopy (W/m2)
    swvegsha    => mlcanopy_inst%swvegsha_canopy    , &  ! Absorbed solar radiation: shaded canopy (W/m2)
    lwveg       => mlcanopy_inst%lwveg_canopy       , &  ! Absorbed longwave radiation: vegetation (W/m2)
    lwvegsun    => mlcanopy_inst%lwvegsun_canopy    , &  ! Absorbed longwave radiation: sunlit canopy (W/m2)
    lwvegsha    => mlcanopy_inst%lwvegsha_canopy    , &  ! Absorbed longwave radiation: shaded canopy (W/m2)
    shveg       => mlcanopy_inst%shveg_canopy       , &  ! Sensible heat flux: vegetation (W/m2)
    shvegsun    => mlcanopy_inst%shvegsun_canopy    , &  ! Sensible heat flux: sunlit canopy (W/m2)
    shvegsha    => mlcanopy_inst%shvegsha_canopy    , &  ! Sensible heat flux: shaded canopy (W/m2)
    lhveg       => mlcanopy_inst%lhveg_canopy       , &  ! Latent heat flux: vegetation (W/m2)
    lhvegsun    => mlcanopy_inst%lhvegsun_canopy    , &  ! Latent heat flux: sunlit canopy (W/m2)
    lhvegsha    => mlcanopy_inst%lhvegsha_canopy    , &  ! Latent heat flux: shaded canopy (W/m2)
    etveg       => mlcanopy_inst%etveg_canopy       , &  ! Water vapor flux: vegetation (mol H2O/m2/s)
    etvegsun    => mlcanopy_inst%etvegsun_canopy    , &  ! Water vapor flux: sunlit canopy (mol H2O/m2/s)
    etvegsha    => mlcanopy_inst%etvegsha_canopy    , &  ! Water vapor flux: shaded canopy (mol H2O/m2/s)
    trveg       => mlcanopy_inst%trveg_canopy       , &  ! Water vapor flux: transpiration (mol H2O/m2/s)
    evveg       => mlcanopy_inst%evveg_canopy       , &  ! Water vapor flux: canopy evaporation (mol H2O/m2/s)
    gppveg      => mlcanopy_inst%gppveg_canopy      , &  ! Gross primary production: vegetation (umol CO2/m2/s)
    gppvegsun   => mlcanopy_inst%gppvegsun_canopy   , &  ! Gross primary production: sunlit canopy (umol CO2/m2/s)
    gppvegsha   => mlcanopy_inst%gppvegsha_canopy   , &  ! Gross primary production: shaded canopy (umol CO2/m2/s)
    vcmax25veg  => mlcanopy_inst%vcmax25veg_canopy  , &  ! Vcmax at 25C: total canopy (umol/m2/s)
    vcmax25sun  => mlcanopy_inst%vcmax25sun_canopy  , &  ! Vcmax at 25C: sunlit canopy (umol/m2/s)
    vcmax25sha  => mlcanopy_inst%vcmax25sha_canopy  , &  ! Vcmax at 25C: shaded canopy (umol/m2/s)
    gsveg       => mlcanopy_inst%gsveg_canopy       , &  ! Stomatal conductance: canopy (mol H2O/m2/s)
    gsvegsun    => mlcanopy_inst%gsvegsun_canopy    , &  ! Stomatal conductance: sunlit canopy (mol H2O/m2/s)
    gsvegsha    => mlcanopy_inst%gsvegsha_canopy    , &  ! Stomatal conductance: shaded canopy (mol H2O/m2/s)
    windveg     => mlcanopy_inst%windveg_canopy     , &  ! Wind speed: canopy (m/s)
    windvegsun  => mlcanopy_inst%windvegsun_canopy  , &  ! Wind speed: sunlit canopy (m/s)
    windvegsha  => mlcanopy_inst%windvegsha_canopy  , &  ! Wind speed: shaded canopy (m/s)
    tlveg       => mlcanopy_inst%tlveg_canopy       , &  ! Leaf temperature: canopy (K)
    tlvegsun    => mlcanopy_inst%tlvegsun_canopy    , &  ! Leaf temperature: sunlit canopy (K)
    tlvegsha    => mlcanopy_inst%tlvegsha_canopy    , &  ! Leaf temperature: shaded canopy (K)
    taveg       => mlcanopy_inst%taveg_canopy       , &  ! Air temperature: canopy (K)
    tavegsun    => mlcanopy_inst%tavegsun_canopy    , &  ! Air temperature: sunlit canopy (K)
    tavegsha    => mlcanopy_inst%tavegsha_canopy    , &  ! Air temperature: shaded canopy (K)
    laisun      => mlcanopy_inst%laisun_canopy      , &  ! Canopy plant area index (lai+sai): sunlit canopy (m2/m2) 
    laisha      => mlcanopy_inst%laisha_canopy      , &  ! Canopy plant area index (lai+sai): shaded canopy (m2/m2) 
    fracminlwp  => mlcanopy_inst%fracminlwp_canopy  , &  ! Fraction of canopy that is water-stressed
    swsrc       => mlcanopy_inst%swsrc_profile      , &  ! Canopy layer source/sink flux: absorbed solar radiation (W/m2)
    lwsrc       => mlcanopy_inst%lwsrc_profile      , &  ! Canopy layer source/sink flux: absorbed longwave radiation (W/m2)
    rnsrc       => mlcanopy_inst%rnsrc_profile      , &  ! Canopy layer source/sink flux: net radiation (W/m2)
    stsrc       => mlcanopy_inst%stsrc_profile      , &  ! Canopy layer source/sink flux: storage heat flux (W/m2)
    shsrc       => mlcanopy_inst%shsrc_profile      , &  ! Canopy layer source/sink flux: sensible heat (W/m2)
    lhsrc       => mlcanopy_inst%lhsrc_profile      , &  ! Canopy layer source/sink flux: latent heat (W/m2)
    etsrc       => mlcanopy_inst%etsrc_profile      , &  ! Canopy layer source/sink flux: water vapor (mol H2O/m2/s)
    trsrc       => mlcanopy_inst%trsrc_profile      , &  ! Canopy layer source/sink flux: transpiration water vapor (mol H2O/m2/s)
    evsrc       => mlcanopy_inst%evsrc_profile      , &  ! Canopy layer source/sink flux: evaporation water vapor (mol H2O/m2/s)
    fco2src     => mlcanopy_inst%fco2src_profile    , &  ! Canopy layer source/sink flux: CO2 (umol CO2/m2/s)
    swleaf_mean => mlcanopy_inst%swleaf_mean_profile, &  ! Canopy layer weighted mean: leaf absorbed solar radiation (W/m2 leaf)
    lwleaf_mean => mlcanopy_inst%lwleaf_mean_profile, &  ! Canopy layer weighted mean: leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf_mean => mlcanopy_inst%rnleaf_mean_profile, &  ! Canopy layer weighted mean: leaf net radiation (W/m2 leaf)
    stleaf_mean => mlcanopy_inst%stleaf_mean_profile, &  ! Canopy layer weighted mean: leaf storage heat flux (W/m2 leaf)
    shleaf_mean => mlcanopy_inst%shleaf_mean_profile, &  ! Canopy layer weighted mean: leaf sensible heat flux (W/m2 leaf)
    lhleaf_mean => mlcanopy_inst%lhleaf_mean_profile, &  ! Canopy layer weighted mean: leaf latent heat flux (W/m2 leaf)
    etleaf_mean => mlcanopy_inst%etleaf_mean_profile, &  ! Canopy layer weighted mean: leaf water vapor flux (mol H2O/m2 leaf/s)
    trleaf_mean => mlcanopy_inst%trleaf_mean_profile, &  ! Canopy layer weighted mean: leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf_mean => mlcanopy_inst%evleaf_mean_profile, &  ! Canopy layer weighted mean: leaf evaporation flux (mol H2O/m2 leaf/s)
    fco2_mean   => mlcanopy_inst%fco2_mean_profile  , &  ! Canopy layer weighted mean: leaf net photosynthesis (umol CO2/m2 leaf/s)
    apar_mean   => mlcanopy_inst%apar_mean_profile  , &  ! Canopy layer weighted mean: leaf absorbed PAR (umol photon/m2 leaf/s)
    gs_mean     => mlcanopy_inst%gs_mean_profile    , &  ! Canopy layer weighted mean: leaf stomatal conductance (mol H2O/m2 leaf/s)
    tleaf_mean  => mlcanopy_inst%tleaf_mean_profile , &  ! Canopy layer weighted mean: leaf temperature (K)
    lwp_mean    => mlcanopy_inst%lwp_mean_profile     &  ! Canopy layer weighted mean: leaf water potential (MPa)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Leaf flux profiles

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then

             ! Leaf fluxes/states (per unit leaf area)

             lwleaf_mean(p,ic) = lwleaf(p,ic,isun)*fracsun(p,ic) + lwleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             swleaf_mean(p,ic,ivis) = swleaf(p,ic,isun,ivis)*fracsun(p,ic) + swleaf(p,ic,isha,ivis)*(1._r8 - fracsun(p,ic))
             swleaf_mean(p,ic,inir) = swleaf(p,ic,isun,inir)*fracsun(p,ic) + swleaf(p,ic,isha,inir)*(1._r8 - fracsun(p,ic))
             rnleaf_mean(p,ic) = rnleaf(p,ic,isun)*fracsun(p,ic) + rnleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             stleaf_mean(p,ic) = stleaf(p,ic,isun)*fracsun(p,ic) + stleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             shleaf_mean(p,ic) = shleaf(p,ic,isun)*fracsun(p,ic) + shleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             lhleaf_mean(p,ic) = lhleaf(p,ic,isun)*fracsun(p,ic) + lhleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             etleaf_mean(p,ic) = (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) &
                               + (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * (1._r8 - fracsun(p,ic))
             trleaf_mean(p,ic) = trleaf(p,ic,isun)*fracsun(p,ic) + trleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             evleaf_mean(p,ic) = evleaf(p,ic,isun)*fracsun(p,ic) + evleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             fco2_mean(p,ic) = anet(p,ic,isun)*fracsun(p,ic) + anet(p,ic,isha)*(1._r8 - fracsun(p,ic))

             apar_mean(p,ic) = apar(p,ic,isun)*fracsun(p,ic) + apar(p,ic,isha)*(1._r8 - fracsun(p,ic))
             gs_mean(p,ic) = gs(p,ic,isun)*fracsun(p,ic) + gs(p,ic,isha)*(1._r8 - fracsun(p,ic))
             tleaf_mean(p,ic) = tleaf(p,ic,isun)*fracsun(p,ic) + tleaf(p,ic,isha)*(1._r8 - fracsun(p,ic))
             lwp_mean(p,ic) = lwp(p,ic,isun)*fracsun(p,ic) + lwp(p,ic,isha)*(1._r8 - fracsun(p,ic))

             ! Source fluxes (per unit ground area)

             lwsrc(p,ic) = lwleaf_mean(p,ic) * dpai(p,ic)
             swsrc(p,ic,ivis) = swleaf_mean(p,ic,ivis) * dpai(p,ic)
             swsrc(p,ic,inir) = swleaf_mean(p,ic,inir) * dpai(p,ic)
             rnsrc(p,ic) = rnleaf_mean(p,ic) * dpai(p,ic)
             stsrc(p,ic) = stleaf_mean(p,ic) * dpai(p,ic)
             shsrc(p,ic) = shleaf_mean(p,ic) * dpai(p,ic)
             lhsrc(p,ic) = lhleaf_mean(p,ic) * dpai(p,ic)
             etsrc(p,ic) = etleaf_mean(p,ic) * dpai(p,ic)
             trsrc(p,ic) = trleaf_mean(p,ic) * dpai(p,ic)
             evsrc(p,ic) = evleaf_mean(p,ic) * dpai(p,ic)
             fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
             fco2src(p,ic) = (anet(p,ic,isun)*fracsun(p,ic) + anet(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic) * fracgreen

          else

             lwleaf_mean(p,ic) = 0._r8
             swleaf_mean(p,ic,ivis) = 0._r8
             swleaf_mean(p,ic,inir) = 0._r8
             rnleaf_mean(p,ic) = 0._r8
             stleaf_mean(p,ic) = 0._r8
             shleaf_mean(p,ic) = 0._r8
             lhleaf_mean(p,ic) = 0._r8
             etleaf_mean(p,ic) = 0._r8
             trleaf_mean(p,ic) = 0._r8
             evleaf_mean(p,ic) = 0._r8
             fco2_mean(p,ic) = 0._r8

             apar_mean(p,ic) = 0._r8
             gs_mean(p,ic) = 0._r8
             tleaf_mean(p,ic) = 0._r8
             lwp_mean(p,ic) = 0._r8

             lwsrc(p,ic) = 0._r8
             swsrc(p,ic,ivis) = 0._r8
             swsrc(p,ic,inir) = 0._r8
             rnsrc(p,ic) = 0._r8
             stsrc(p,ic) = 0._r8
             shsrc(p,ic) = 0._r8
             lhsrc(p,ic) = 0._r8
             etsrc(p,ic) = 0._r8
             trsrc(p,ic) = 0._r8
             evsrc(p,ic) = 0._r8
             fco2src(p,ic) = 0._r8
          end if
       end do

       ! Sum leaf fluxes over canopy

       swveg(p,ivis) = 0._r8
       swveg(p,inir) = 0._r8
       lwveg(p) = 0._r8
       stflx_veg(p) = 0._r8
       shveg(p) = 0._r8
       lhveg(p) = 0._r8
       etveg(p) = 0._r8
       trveg(p) = 0._r8
       evveg(p) = 0._r8
       gppveg(p) = 0._r8
       vcmax25veg(p) = 0._r8
       gsveg(p) = 0._r8

       do ic = 1, ncan(p)
          swveg(p,ivis) = swveg(p,ivis) + swsrc(p,ic,ivis)
          swveg(p,inir) = swveg(p,inir) + swsrc(p,ic,inir)
          lwveg(p) = lwveg(p) + lwsrc(p,ic)
          stflx_veg(p) = stflx_veg(p) + stsrc(p,ic)
          shveg(p) = shveg(p) + shsrc(p,ic)
          lhveg(p) = lhveg(p) + lhsrc(p,ic)
          etveg(p) = etveg(p) + etsrc(p,ic)
          trveg(p) = trveg(p) + trsrc(p,ic)
          evveg(p) = evveg(p) + evsrc(p,ic)
          if (dpai(p,ic) > 0._r8) then
             fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
             gppveg(p) = gppveg(p) &
                       + (agross(p,ic,isun)*fracsun(p,ic) + agross(p,ic,isha)*(1._r8 - fracsun(p,ic))) * dpai(p,ic) * fracgreen
             gsveg(p) = gsveg(p) + gs_mean(p,ic) * dpai(p,ic)
          end if
          vcmax25veg(p) = vcmax25veg(p) + vcmax25_profile(p,ic) * dpai(p,ic)
       end do

       ! Check vegetation energy balance for conservation

       err = swveg(p,ivis) + swveg(p,inir) + lwveg(p) - shveg(p) - lhveg(p) - stflx_veg(p)
       if (abs(err) >= 1.e-03_r8) then
          call endrun (msg=' ERROR: CanopyFluxesDiagnostics: energy conservation error (1)')
       end if

       ! Albedo

       do ib = 1, numrad
          radin = swskyb(p,ib) + swskyd(p,ib)
          if (radin > 0._r8) then
             albcan(p,ib) = swupw(p,ntop(p),ib) / radin
          else
             albcan(p,ib) = 0._r8
          end if
       end do

       ! Turbulent fluxes

       select case (flux_profile_type)
       case (0, -1)
          ! Sum of vegetation and soil fluxes
          shflx(p) = shveg(p) + shsoi(p)
          etflx(p) = etveg(p) + etsoi(p)
          lhflx(p) = lhveg(p) + lhsoi(p)
       case (1)
          ! Turbulent fluxes are at the top layer
          shflx(p) = shair(p,ncan(p))
          etflx(p) = etair(p,ncan(p))
          lhflx(p) = etair(p,ncan(p)) * LatVap(tref(p))
       case default
          call endrun (msg=' ERROR: CanopyFluxesDiagnostics: turbulence type not valid')
       end select

       ! Air heat storage flux

       stflx_air(p) = 0._r8
       do ic = 1, ncan(p)
          stflx_air(p) = stflx_air(p) + stair(p,ic)
       end do

       ! Overall energy balance check:
       ! radiation in - radiation out - soil heat = available energy = turbulent flux + canopy storage flux

       rnet(p) = swveg(p,ivis) + swveg(p,inir) + swsoi(p,ivis) + swsoi(p,inir) + lwveg(p) + lwsoi(p)
       radin = swskyb(p,ivis) + swskyd(p,ivis) + swskyb(p,inir) + swskyd(p,inir) + lwsky(p)
       radout = albcan(p,ivis)*(swskyb(p,ivis)+swskyd(p,ivis)) + albcan(p,inir)*(swskyb(p,inir)+swskyd(p,inir)) + lwup(p)

       err = rnet(p) - (radin - radout)
       if (abs(err) > 0.001_r8) then
          call endrun (msg=' ERROR: CanopyFluxesDiagnostics: energy conservation error (2)')
       end if

       avail = radin - radout - gsoi(p)
       flux = shflx(p) + lhflx(p) + stflx_air(p) + stflx_veg(p)
       err = avail - flux
       if (abs(err) > 0.01_r8) then
          call endrun (msg=' ERROR: CanopyFluxesDiagnostics: energy conservation error (3)')
       end if

       ! Check radiative fluxes at top of canopy

       ic = ntop(p)
       radin = swbeam(p,ic,ivis) + swbeam(p,ic,inir) + swdwn(p,ic,ivis) + swdwn(p,ic,inir) + lwdwn(p,ic)
       radout = swupw(p,ic,ivis) + swupw(p,ic,inir) + lwupw(p,ic)
       err = (radin - radout) - rnet(p)
       if (abs(err) > 0.001_r8) then
          call endrun (msg=' ERROR: CanopyFluxesDiagnostics: energy conservation error (4)')
       end if

       ! Sunlit and shaded canopy fluxes

       laisun(p) = 0._r8 ; laisha(p) = 0._r8
       swvegsun(p,1:numrad) = 0._r8 ; swvegsha(p,1:numrad) = 0._r8
       lwvegsun(p) = 0._r8 ; lwvegsha(p) = 0._r8
       shvegsun(p) = 0._r8 ; shvegsha(p) = 0._r8
       lhvegsun(p) = 0._r8 ; lhvegsha(p) = 0._r8
       etvegsun(p) = 0._r8 ; etvegsha(p) = 0._r8
       gppvegsun(p) = 0._r8 ; gppvegsha(p) = 0._r8
       vcmax25sun(p) = 0._r8 ; vcmax25sha(p) = 0._r8
       gsvegsun(p) = 0._r8 ; gsvegsha(p) = 0._r8

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             laisun(p) = laisun(p) + fracsun(p,ic) * dpai(p,ic)
             laisha(p) = laisha(p) + (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             swvegsun(p,ivis) = swvegsun(p,ivis) + swleaf(p,ic,isun,ivis) * fracsun(p,ic) * dpai(p,ic)
             swvegsun(p,inir) = swvegsun(p,inir) + swleaf(p,ic,isun,inir) * fracsun(p,ic) * dpai(p,ic)
             swvegsha(p,ivis) = swvegsha(p,ivis) + swleaf(p,ic,isha,ivis) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             swvegsha(p,inir) = swvegsha(p,inir) + swleaf(p,ic,isha,inir) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             lwvegsun(p) = lwvegsun(p) + lwleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             lwvegsha(p) = lwvegsha(p) + lwleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             shvegsun(p) = shvegsun(p) + shleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             shvegsha(p) = shvegsha(p) + shleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             lhvegsun(p) = lhvegsun(p) + lhleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             lhvegsha(p) = lhvegsha(p) + lhleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             etvegsun(p) = etvegsun(p) + (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic)
             etvegsha(p) = etvegsha(p) + (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             fracgreen = fdry(p,ic) / (1._r8 - fwet(p,ic))
             gppvegsun(p) = gppvegsun(p) + agross(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) * fracgreen
             gppvegsha(p) = gppvegsha(p) + agross(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) * fracgreen
             vcmax25sun(p) = vcmax25sun(p) + vcmax25_leaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             vcmax25sha(p) = vcmax25sha(p) + vcmax25_leaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             gsvegsun(p) = gsvegsun(p) + gs(p,ic,isun) * fracsun(p,ic) * dpai(p,ic)
             gsvegsha(p) = gsvegsha(p) + gs(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
          end if
       end do

       ! Sunlit and shaded canopy temperatures and wind speed are weighted for sun/shade leaf area

       windveg(p) = 0._r8 ; windvegsun(p) = 0._r8 ; windvegsha(p) = 0._r8
       tlveg(p) = 0._r8 ; tlvegsun(p) = 0._r8 ; tlvegsha(p) = 0._r8
       taveg(p) = 0._r8 ; tavegsun(p) = 0._r8 ; tavegsha(p) = 0._r8
       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             windveg(p) = windveg(p) + wind(p,ic) * dpai(p,ic) / (laisun(p) + laisha(p))
             windvegsun(p) = windvegsun(p) + wind(p,ic) * fracsun(p,ic) * dpai(p,ic) / laisun(p)
             windvegsha(p) = windvegsha(p) + wind(p,ic) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) / laisha(p)

             tlveg(p) = tlveg(p) + tleaf_mean(p,ic) * dpai(p,ic) / (laisun(p) + laisha(p))
             tlvegsun(p) = tlvegsun(p) + tleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) / laisun(p)
             tlvegsha(p) = tlvegsha(p) + tleaf(p,ic,isha) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) / laisha(p)

             taveg(p) = taveg(p) + tair(p,ic) * dpai(p,ic) / (laisun(p) + laisha(p))
             tavegsun(p) = tavegsun(p) + tair(p,ic) * fracsun(p,ic) * dpai(p,ic) / laisun(p)
             tavegsha(p) = tavegsha(p) + tair(p,ic) * (1._r8 - fracsun(p,ic)) * dpai(p,ic) / laisha(p)
          end if
       end do

       ! Diagnose fraction of the canopy that is water stressed

       minlwp = -2._r8
       fracminlwp(p) = 0._r8

       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) then
             if (lwp_mean(p,ic) <= minlwp) then
                fracminlwp(p) = fracminlwp(p) + dpai(p,ic)
             end if
          end if
       end do

       if ((lai(p) + sai(p)) > 0._r8) then
          fracminlwp(p) = fracminlwp(p) / (lai(p) + sai(p))
       end if

    end do

    end associate
  end subroutine CanopyFluxesDiagnostics

end module MLCanopyFluxesMod
