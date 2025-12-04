module MLFluxProfileSolutionMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Source/sink fluxes for leaves & soil and concentration profiles
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
  public :: FluxProfileSolution           ! Main routine for source/sink fluxes and concentration profiles
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: ImplicitFluxProfileSolution  ! Implicit solution for source/sink fluxes and concentration profiles
  private :: WellMixed                    ! Canopy scalar profiles equal reference height values
  private :: ErrorCheck01                 ! Implicit solution error checks
  private :: ErrorCheck02                 ! Conservation error checks
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine FluxProfileSolution (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Source/sink fluxes for leaves & soil and concentration profiles
    !
    ! !USES:
    use MLclm_varctl, only : flux_profile_type
    use MLCanopyFluxesType, only : mlcanopy_type
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
    !---------------------------------------------------------------------

    associate ( &
                                                     ! *** Input ***
    co2ref    => mlcanopy_inst%co2ref_forcing   , &  ! Atmospheric CO2 at reference height (umol/mol)
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
                                                     ! *** Output ***
    cair      => mlcanopy_inst%cair_profile       &  ! Canopy layer atmospheric CO2 (umol/mol)
    )

    select case (flux_profile_type)
    case (0, -1)

       ! Use the well-mixed assumption, or read profile data from dataset

       call WellMixed (num_filter, filter, mlcanopy_inst)

       ! This is legacy code and is not supported. It is included here in case it is
       ! ever needed again.

       call endrun (msg=' ERROR: WellMixed option not valid')

    case (1)

       ! Source/sink fluxes and concentration profiles using implicit solution

       do fp = 1, num_filter
          p = filter(fp)
    
          call ImplicitFluxProfileSolution (p, mlcanopy_inst)
 
          ! No profile for CO2
 
          do ic = 1, ncan(p)
             cair(p,ic) = co2ref(p)
          end do

       end do

    case default

       call endrun (msg=' ERROR: FluxProfileSolution: flux_profile_type not valid')

    end select

    end associate
  end subroutine FluxProfileSolution

  !-----------------------------------------------------------------------
  subroutine ImplicitFluxProfileSolution (p, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute source/sink fluxes for leaves & soil and concentration
    ! profiles. This uses an implicit solution for temperature and vapor
    ! pressure. The boundary conditions are the above-canopy scalar values at
    ! the reference height (T, q) and the temperature of the first soil layer.
    ! Vegetation and ground temperatures and fluxes are calculated as part of
    ! the implicit solution.
    !
    ! !USES:
    use MLclm_varctl, only : dtime_ml
    use MLclm_varpar, only : isun, isha, nlevmlcan, nleaf
    use MLWaterVaporMod, only : SatVap, LatVap
    use MLMathToolsMod, only: tridiag_2eq
    use MLLeafFluxesMod, only : LeafFluxes
    use MLSoilFluxesMod, only : SoilFluxes
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p                  ! Patch index for CLM g/l/c/p hierarchy
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                            ! Aboveground layer index
    integer  :: il                            ! Sunlit (1) or shaded (2) leaf index
    real(r8) :: dtime                         ! Multilayer canopy timestep (s)
    real(r8) :: lambda                        ! Latent heat of vaporization (J/mol)
    real(r8) :: esat                          ! Saturation vapor pressure (Pa)
    real(r8) :: desat                         ! Temperature derivative of saturation vapor pressure (Pa/K)
    real(r8) :: qsat                          ! Saturation vapor pressure (mol/mol)
    real(r8) :: den                           ! Intermediate calculation
    real(r8) :: pai                           ! Plant area index of layer (m2/m2)
    real(r8) :: gac_ic_minus_one              ! Special case for ic-1 = 0: use soil conductance not ic-1
    real(r8) :: qsat0                         ! Saturation vapor pressure at ground (mol/mol)
    real(r8) :: dqsat0                        ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: gsw                           ! Soil conductance for water vapor (mol H2O/m2/s)
    real(r8) :: gs0                           ! Total soil-to-air conductance for water vapor (mol H2O/m2/s)
    real(r8) :: alpha0                        ! Coefficient for ground temperature (dimensionless)
    real(r8) :: beta0                         ! Coefficient for ground temperature (K mol/mol)
    real(r8) :: delta0                        ! Coefficient for ground temperature (K)
    real(r8) :: c01                           ! Soil heat flux term (W/m2)
    real(r8) :: c02                           ! Soil heat flux term (W/m2/K)
    real(r8) :: sh0                           ! Implicit solution: Ground sensible heat flux (W/m2)
    real(r8) :: et0                           ! Implicit solution: Ground evaporation flux (mol H2O/m2/s)
    real(r8) :: g0                            ! Implicit solution: Soil heat flux (W/m2)
    real(r8) :: t0                            ! Implicit solution: Soil surface temperature (K)
    real(r8) :: e0                            ! Implicit solution: Soil surface vapor pressure (Pa)

    real(r8) :: rho_dz_over_dt(nlevmlcan)     ! Intermediate calculation for canopy air storage
    real(r8) :: gleaf_sh(nlevmlcan,nleaf)     ! Leaf conductance for sensible heat (mol/m2/s)
    real(r8) :: gleaf_et(nlevmlcan,nleaf)     ! Leaf conductance for water vapor (mol/m2/s)
    real(r8) :: heatcap(nlevmlcan,nleaf)      ! Heat capacity of leaves (J/m2/K)
    real(r8) :: avail_energy(nlevmlcan,nleaf) ! Available energy for leaf (W/m2)
    real(r8) :: dqsat(nlevmlcan,nleaf)        ! Temperature derivative of saturation vapor pressure (mol/mol/K)
    real(r8) :: qsat_term(nlevmlcan,nleaf)    ! Intermediate calculation for saturation vapor pressure (mol/mol)
    real(r8) :: alpha(nlevmlcan,nleaf)        ! Coefficient for leaf temperature (dimensionless)
    real(r8) :: beta(nlevmlcan,nleaf)         ! Coefficient for leaf temperature (K mol/mol)
    real(r8) :: delta(nlevmlcan,nleaf)        ! Coefficient for leaf temperature (K)

    real(r8) :: a1(nlevmlcan)                 ! Coefficient for canopy air temperature
    real(r8) :: b11(nlevmlcan)                ! Coefficient for canopy air temperature
    real(r8) :: b12(nlevmlcan)                ! Coefficient for canopy air temperature
    real(r8) :: c1(nlevmlcan)                 ! Coefficient for canopy air temperature
    real(r8) :: d1(nlevmlcan)                 ! Coefficient for canopy air temperature

    real(r8) :: a2(nlevmlcan)                 ! Coefficient for canopy air water vapor mole fraction
    real(r8) :: b21(nlevmlcan)                ! Coefficient for canopy air water vapor mole fraction
    real(r8) :: b22(nlevmlcan)                ! Coefficient for canopy air water vapor mole fraction
    real(r8) :: c2(nlevmlcan)                 ! Coefficient for canopy air water vapor mole fraction
    real(r8) :: d2(nlevmlcan)                 ! Coefficient for canopy air water vapor mole fraction

    real(r8) :: tleaf_implic(nlevmlcan,nleaf) ! Implicit solution: Leaf temperature (K)
    real(r8) :: storage_sh(nlevmlcan)         ! Implicit solution: Heat storage flux in air (W/m2)
    real(r8) :: storage_et(nlevmlcan)         ! Implicit solution: Water vapor storage flux in air (mol H2O/m2/s)
    real(r8) :: stveg(nlevmlcan)              ! Implicit solution: Canopy layer leaf storage heat flux (W/m2)
    real(r8) :: shsrc(nlevmlcan)              ! Implicit solution: Canopy layer leaf sensible heat flux (W/m2)
    real(r8) :: etsrc(nlevmlcan)              ! Implicit solution: Canopy layer leaf water vapor flux (mol H2O/m2/s)
    !---------------------------------------------------------------------

    associate ( &
                                                     ! *** Input ***
    tref      => mlcanopy_inst%tref_forcing     , &  ! Air temperature at reference height (K)
    thref     => mlcanopy_inst%thref_forcing    , &  ! Atmospheric potential temperature at reference height (K)
    eref      => mlcanopy_inst%eref_forcing     , &  ! Vapor pressure at reference height (Pa)
    pref      => mlcanopy_inst%pref_forcing     , &  ! Air pressure at reference height (Pa)
    rhomol    => mlcanopy_inst%rhomol_forcing   , &  ! Molar density at reference height (mol/m3)
    cpair     => mlcanopy_inst%cpair_forcing    , &  ! Specific heat of air (constant pressure) at reference height (J/mol/K)
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    gac0      => mlcanopy_inst%gac0_soil        , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    tg_bef    => mlcanopy_inst%tg_bef_soil      , &  ! Soil surface temperature for previous timestep (K)
    rhg       => mlcanopy_inst%rhg_soil         , &  ! Relative humidity of airspace at soil surface (fraction)
    rnsoi     => mlcanopy_inst%rnsoi_soil       , &  ! Net radiation: ground (W/m2)
    soilres   => mlcanopy_inst%soilres_soil     , &  ! Soil evaporative resistance (s/m)
    soil_t    => mlcanopy_inst%soil_t_soil      , &  ! Temperature of first snow/soil layer (K)
    soil_dz   => mlcanopy_inst%soil_dz_soil     , &  ! Depth to temperature of first snow/soil layer (m)
    soil_tk   => mlcanopy_inst%soil_tk_soil     , &  ! Thermal conductivity of first snow/soil layer (W/m/K)
    dz        => mlcanopy_inst%dz_profile       , &  ! Canopy layer thickness (m)
    dpai      => mlcanopy_inst%dpai_profile     , &  ! Canopy layer plant area index (m2/m2)
    fwet      => mlcanopy_inst%fwet_profile     , &  ! Canopy layer fraction of plant area index that is wet
    fdry      => mlcanopy_inst%fdry_profile     , &  ! Canopy layer fraction of plant area index that is green and dry
    fracsun   => mlcanopy_inst%fracsun_profile  , &  ! Canopy layer sunlit fraction (-)
    cpleaf    => mlcanopy_inst%cpleaf_profile   , &  ! Canopy layer leaf heat capacity (J/m2 leaf/K)
    gac       => mlcanopy_inst%gac_profile      , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    tair_bef  => mlcanopy_inst%tair_bef_profile , &  ! Canopy layer air temperature for previous timestep (K)
    eair_bef  => mlcanopy_inst%eair_bef_profile , &  ! Canopy layer vapor pressure for previous timestep (Pa)
    gbh       => mlcanopy_inst%gbh_leaf         , &  ! Leaf boundary layer conductance: heat (mol/m2 leaf/s)
    gbv       => mlcanopy_inst%gbv_leaf         , &  ! Leaf boundary layer conductance: H2O (mol H2O/m2 leaf/s)
    gs        => mlcanopy_inst%gs_leaf          , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    rnleaf    => mlcanopy_inst%rnleaf_leaf      , &  ! Leaf net radiation (W/m2 leaf)
    tleaf_bef => mlcanopy_inst%tleaf_bef_leaf   , &  ! Leaf temperature for previous timestep (K)
                                                     ! *** Output ***
    tair      => mlcanopy_inst%tair_profile     , &  ! Canopy layer air temperature (K)
    eair      => mlcanopy_inst%eair_profile     , &  ! Canopy layer vapor pressure (Pa)
    shair     => mlcanopy_inst%shair_profile    , &  ! Canopy layer air sensible heat flux (W/m2)
    etair     => mlcanopy_inst%etair_profile    , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair     => mlcanopy_inst%stair_profile    , &  ! Canopy layer air storage heat flux (W/m2)
                                                     ! *** Output from LeafFluxes
    tleaf     => mlcanopy_inst%tleaf_leaf       , &  ! Leaf temperature (K)
    stleaf    => mlcanopy_inst%stleaf_leaf      , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf    => mlcanopy_inst%shleaf_leaf      , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf    => mlcanopy_inst%lhleaf_leaf      , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf    => mlcanopy_inst%trleaf_leaf      , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf    => mlcanopy_inst%evleaf_leaf      , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
                                                     ! *** Output from SoilFluxes
    shsoi     => mlcanopy_inst%shsoi_soil       , &  ! Sensible heat flux: ground (W/m2)
    lhsoi     => mlcanopy_inst%lhsoi_soil       , &  ! Latent heat flux: ground (W/m2)
    etsoi     => mlcanopy_inst%etsoi_soil       , &  ! Water vapor flux: ground (mol H2O/m2/s)
    gsoi      => mlcanopy_inst%gsoi_soil        , &  ! Soil heat flux (W/m2)
    eg        => mlcanopy_inst%eg_soil          , &  ! Soil surface vapor pressure (Pa)
    tg        => mlcanopy_inst%tg_soil            &  ! Soil surface temperature (K)
    )

    ! Timestep (s)

    dtime = dtime_ml

    ! Latent heat of vaporization

    lambda = LatVap(tref(p))

    ! Terms for ground temperature, which is calculated from the energy balance:
    !
    ! Rn0 - H0 - lambda*E0 - G = 0
    !
    ! and is rewritten in relation to T and q as:
    !
    ! T0 = alpha0*T(1) + beta0*q(1) + delta0
    !
    ! by substituting the flux equations for H0, E0, and G
    !
    ! See Bonan et al. (2018) Geosci. Model Dev., 11, 1467-1496, doi:10.5194/gmd-11-1467-2018, eqs. (14)-(15), (A6)-(A9)

    call SatVap (tg_bef(p), esat, desat)               ! Vapor pressure (Pa) at ground temperature
    qsat0 = esat / pref(p) ; dqsat0 = desat / pref(p)  ! Pa -> mol/mol

    gsw = (1._r8 / soilres(p)) * rhomol(p)             ! Soil conductance for water vapor: s/m -> mol H2O/m2/s
    gs0 = gac0(p) * gsw / (gac0(p) + gsw)              ! Total soil-to-air conductance, including aerodynamic term

    c02 = soil_tk(p) / soil_dz(p)                      ! Soil heat flux term (W/m2/K)
    c01 = -c02 * soil_t(p)                             ! Soil heat flux term (W/m2)

    den = cpair(p) * gac0(p) + lambda * rhg(p) * gs0 * dqsat0 + c02
    alpha0 = cpair(p) * gac0(p) / den
    beta0 = lambda * gs0 / den
    delta0 = (rnsoi(p) - lambda * rhg(p) * gs0 * (qsat0 - dqsat0 * tg_bef(p)) - c01) / den

    ! Similarly, re-arrange the leaf energy balance to calculate leaf
    ! temperature in relation to T and q:
    !
    ! Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
    ! Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)
    !
    ! See Bonan et al. (2018), eqs. (10)-(13), (A1)-(A5)

    do ic = 1, ncan(p)

       ! Calculate terms for sunlit and shaded leaves

       if (dpai(p,ic) > 0._r8) then

          do il = 1, nleaf

             ! Leaf conductances

             gleaf_sh(ic,il) = 2._r8 * gbh(p,ic,il)
             gleaf_et(ic,il) = gs(p,ic,il)*gbv(p,ic,il)/(gs(p,ic,il)+gbv(p,ic,il)) * fdry(p,ic) + gbv(p,ic,il) * fwet(p,ic)

             ! Heat capacity of leaves

             heatcap(ic,il) = cpleaf(p,ic)

             ! Available energy: net radiation

             avail_energy(ic,il) = rnleaf(p,ic,il)

             ! Saturation vapor pressure and derivative for leaf temperature at time n: Pa -> mol/mol

             call SatVap (tleaf_bef(p,ic,il), esat, desat)
             qsat = esat / pref(p) ; dqsat(ic,il) = desat / pref(p)

             ! Term for linearized vapor pressure at leaf temperature:
             ! qsat(tleaf) = qsat(tleaf_bef) + dqsat * (tleaf - tleaf_bef)
             ! Here, qsat_term contains the terms with tleaf_bef

             qsat_term(ic,il) = qsat - dqsat(ic,il) * tleaf_bef(p,ic,il)

             ! alpha, beta, delta coefficients for leaf temperature

             den = heatcap(ic,il) / dtime + gleaf_sh(ic,il) * cpair(p) + gleaf_et(ic,il) * lambda * dqsat(ic,il)
             alpha(ic,il) = gleaf_sh(ic,il) * cpair(p) / den
             beta(ic,il) = gleaf_et(ic,il) * lambda / den
             delta(ic,il) = avail_energy(ic,il) / den &
                          - lambda * gleaf_et(ic,il) * qsat_term(ic,il) / den &
                          + heatcap(ic,il) / dtime * tleaf_bef(p,ic,il) / den

             ! Now scale flux terms for leaf area so that fluxes are for the canopy layer

             if (il == isun) then
                pai = dpai(p,ic) * fracsun(p,ic)
             else if (il == isha) then
                pai = dpai(p,ic) * (1._r8 - fracsun(p,ic))
             end if

             gleaf_sh(ic,il) = gleaf_sh(ic,il) * pai
             gleaf_et(ic,il) = gleaf_et(ic,il) * pai
             heatcap(ic,il) = heatcap(ic,il) * pai
             avail_energy(ic,il) = avail_energy(ic,il) * pai

          end do

       else

          ! Zero out terms

          do il = 1, nleaf
             gleaf_sh(ic,il) = 0._r8
             gleaf_et(ic,il) = 0._r8
             heatcap(ic,il) = 0._r8
             avail_energy(ic,il) = 0._r8
             dqsat(ic,il) = 0._r8
             qsat_term(ic,il) = 0._r8
             alpha(ic,il) = 0._r8
             beta(ic,il) = 0._r8
             delta(ic,il) = 0._r8
          end do

       end if

    end do

    ! The system of equations for air temperature (K) and water vapor (mol/mol)
    ! at each layer is:
    !
    ! a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
    ! a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
    !
    ! These equations are obtained by substituting the equations:
    !
    ! Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
    ! Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)
    !
    ! and:
    !
    ! T0 = alpha0*T(1) + beta0*q(1) + delta0
    !
    ! into the one-dimensional scalar conservation equations for T and q.
    !
    ! Calculate the:
    ! a1, b11, b12, c1, d1 coefficients for air temperature
    ! a2, b21, b22, c2, d2 coefficients for water vapor mole fraction
    !
    ! See Bonan et al. (2018), eqs. (16)-(17) and (S10)-(S31)

    do ic = 1, ncan(p)

       ! Storage term

       rho_dz_over_dt(ic) = rhomol(p) * dz(p,ic) / dtime

       ! a1,b11,b12,c1,d1 coefficients for air temperature

       if (ic == 1) then
          gac_ic_minus_one = gac0(p)
       else
          gac_ic_minus_one = gac(p,ic-1)
       end if

       a1(ic) = -gac_ic_minus_one
       b11(ic) = rho_dz_over_dt(ic) + gac_ic_minus_one + gac(p,ic) &
               + gleaf_sh(ic,isun) * (1._r8 - alpha(ic,isun)) + gleaf_sh(ic,isha) * (1._r8 - alpha(ic,isha))
       b12(ic) = -gleaf_sh(ic,isun) * beta(ic,isun) - gleaf_sh(ic,isha) * beta(ic,isha)
       c1(ic) = -gac(p,ic)
       d1(ic) = rho_dz_over_dt(ic) * tair_bef(p,ic) + gleaf_sh(ic,isun) * delta(ic,isun) + gleaf_sh(ic,isha) * delta(ic,isha)

       ! Special case for top layer

       if (ic == ncan(p)) then
          c1(ic) = 0._r8
          d1(ic) = d1(ic) + gac(p,ic) * thref(p)
       end if

       ! Special case for first canopy layer (i.e., immediately above the ground)

       if (ic == 1) then
          a1(ic) = 0._r8
          b11(ic) = b11(ic) - gac0(p) * alpha0
          b12(ic) = b12(ic) - gac0(p) * beta0
          d1(ic) = d1(ic) + gac0(p) * delta0
       end if

       ! a2,b21,b22,c2,d2 coefficients for water vapor mole fraction

       if (ic == 1) then
          gac_ic_minus_one = gs0
       else
          gac_ic_minus_one = gac(p,ic-1)
       end if

       a2(ic) = -gac_ic_minus_one
       b21(ic) = -gleaf_et(ic,isun) * dqsat(ic,isun) * alpha(ic,isun) - gleaf_et(ic,isha) * dqsat(ic,isha) * alpha(ic,isha)
       b22(ic) = rho_dz_over_dt(ic) + gac_ic_minus_one + gac(p,ic) &
               + gleaf_et(ic,isun) * (1._r8 - dqsat(ic,isun) * beta(ic,isun)) &
               + gleaf_et(ic,isha) * (1._r8 - dqsat(ic,isha) * beta(ic,isha))
       c2(ic) = -gac(p,ic)
       d2(ic) = rho_dz_over_dt(ic) * (eair_bef(p,ic) / pref(p)) &
              + gleaf_et(ic,isun) * (dqsat(ic,isun) * delta(ic,isun) + qsat_term(ic,isun)) &
              + gleaf_et(ic,isha) * (dqsat(ic,isha) * delta(ic,isha) + qsat_term(ic,isha)) 

       ! Special case for top layer

       if (ic == ncan(p)) then
          c2(ic) = 0._r8
          d2(ic) = d2(ic) + gac(p,ic) * (eref(p) / pref(p))
       end if

       ! Special case for first canopy layer (i.e., immediately above the ground)

       if (ic == 1) then
          a2(ic) = 0._r8
          b21(ic) = b21(ic) - gs0 * rhg(p) * dqsat0 * alpha0
          b22(ic) = b22(ic) - gs0 * rhg(p) * dqsat0 * beta0
          d2(ic) = d2(ic) + gs0 * rhg(p) * (qsat0 + dqsat0 * (delta0 - tg_bef(p)))
       end if

    end do

    ! Solve for air temperature and water vapor (mol/mol):
    !
    ! a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
    ! a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
    !
    ! Note that as used here eair = mol/mol

    call tridiag_2eq (a1, b11, b12, c1, d1, a2, b21, b22, c2, d2, tair(p,:), eair(p,:), ncan(p))

    ! Soil surface temperature (K) and vapor pressure (mol/mol)

    t0 = alpha0 * tair(p,1) + beta0 * eair(p,1) + delta0
    e0 = rhg(p) * (qsat0 + dqsat0 * (t0 - tg_bef(p)))

    ! Leaf temperature

    do ic = 1, ncan(p)
       tleaf_implic(ic,isun) = alpha(ic,isun)*tair(p,ic) + beta(ic,isun)*eair(p,ic) + delta(ic,isun)
       tleaf_implic(ic,isha) = alpha(ic,isha)*tair(p,ic) + beta(ic,isha)*eair(p,ic) + delta(ic,isha)
    end do

    ! Convert water vapor from mol/mol to Pa

    do ic = 1, ncan(p)
       eair(p,ic) = eair(p,ic) * pref(p)
    end do
    e0 = e0 * pref(p)

    ! Use LeafFluxes to calculate leaf fluxes (per unit leaf area) for the
    ! current air temperature and vapor pressure profiles in the canopy. The
    ! flux and leaf temperature calculations there are the same as here, so
    ! the answers are the same in both routines. Compare answers to make sure
    ! the implicit flux-profile solution is correct.

    do ic = 1, ncan(p)
       call LeafFluxes (p, ic, isun, mlcanopy_inst)
       call LeafFluxes (p, ic, isha, mlcanopy_inst)
    end do

    ! Use SoilFluxes to calculate soil fluxes. The flux and temperature
    ! calculations there are the same as here, so the answers are the
    ! same in both routines. Compare answers to make sure the
    ! implicit flux-profile solution is correct.

    call SoilFluxes (p, mlcanopy_inst)

    ! Vertical sensible heat and water vapor fluxes between layers

    do ic = 1, ncan(p)-1
       shair(p,ic) = -cpair(p) * (tair(p,ic+1) - tair(p,ic)) * gac(p,ic)
       etair(p,ic) = -(eair(p,ic+1) - eair(p,ic)) / pref(p) * gac(p,ic)
    end do
    ic = ncan(p)
    shair(p,ic) = -cpair(p) * (thref(p) - tair(p,ic)) * gac(p,ic)
    etair(p,ic) = -(eref(p) - eair(p,ic)) / pref(p) * gac(p,ic)

    ! Canopy air storage flux (W/m2) and its component terms

    do ic = 1, ncan(p)
       storage_sh(ic) = cpair(p) * (tair(p,ic) - tair_bef(p,ic)) * rho_dz_over_dt(ic)
       storage_et(ic) = (eair(p,ic) - eair_bef(p,ic)) / pref(p) * rho_dz_over_dt(ic)
       stair(p,ic) = storage_sh(ic) + storage_et(ic) * lambda
    end do

    !-----------------------------------------------------------------------
    ! Error checks: compare solution here with LeafFluxes and SoilFluxes
    ! and also check energy balance
    !-----------------------------------------------------------------------

    ! Source fluxes for each layer as calculated from implicit solution. Remember
    ! that here the fluxes for sunlit/shaded leaves are multiplied by their leaf area.

    do ic = 1, ncan(p)
       shsrc(ic) = 0._r8; etsrc(ic) = 0._r8; stveg(ic) = 0._r8

       if (dpai(p,ic) > 0._r8) then
          do il = 1, nleaf
             shsrc(ic) = shsrc(ic) + cpair(p) * (tleaf_implic(ic,il) - tair(p,ic)) * gleaf_sh(ic,il)
             call SatVap (tleaf_bef(p,ic,il), esat, desat)
             etsrc(ic) = etsrc(ic) + (esat + desat * (tleaf_implic(ic,il) - tleaf_bef(p,ic,il)) - eair(p,ic)) / pref(p) &
                       * gleaf_et(ic,il)
             stveg(ic) = stveg(ic) + heatcap(ic,il) * (tleaf_implic(ic,il) - tleaf_bef(p,ic,il)) / dtime
          end do
       end if
    end do

    ! Soil fluxes for each layer as calculated from implicit solution

    sh0 = -cpair(p) * (tair(p,1) - t0) * gac0(p)
    et0 = -(eair(p,1) - e0) / pref(p) * gs0
    g0 = -soil_tk(p) / soil_dz(p) * soil_t(p) + soil_tk(p) / soil_dz(p) * t0

    ! Compare implicit solution withe LeafFluxes and SoilFluxes

    call ErrorCheck01 (p, lambda, shsrc, etsrc, stveg, tleaf_implic, &
    sh0, et0, g0, t0, e0, mlcanopy_inst)

    ! Conservation error checks

    call ErrorCheck02 (p, lambda, avail_energy, shsrc, etsrc, stveg, storage_sh, storage_et, &
    sh0, et0, g0, mlcanopy_inst)

    end associate
  end subroutine ImplicitFluxProfileSolution

  !-----------------------------------------------------------------------
  subroutine ErrorCheck01 (p, lambda, shsrc, etsrc, stveg, tleaf_implic, &
  sh0, et0, g0, t0, e0, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Implicit solution error checks
    !
    ! !USES:
    use MLclm_varpar, only : isun, isha, nlevmlcan, nleaf
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p                  ! Patch index for CLM g/l/c/p hierarchy
    real(r8) :: lambda                        ! Latent heat of vaporization (J/mol)
    real(r8) :: shsrc(nlevmlcan)              ! Imlicit solution: Canopy layer leaf sensible heat flux (W/m2)
    real(r8) :: etsrc(nlevmlcan)              ! Imlicit solution: Canopy layer leaf water vapor flux (mol H2O/m2/s)
    real(r8) :: stveg(nlevmlcan)              ! Imlicit solution: Canopy layer leaf storage heat flux (W/m2)
    real(r8) :: tleaf_implic(nlevmlcan,nleaf) ! Imlicit solution: Leaf temperature (K)
    real(r8) :: sh0                           ! Imlicit solution: Ground sensible heat flux (W/m2)
    real(r8) :: et0                           ! Imlicit solution: Ground evaporation flux (mol H2O/m2/s)
    real(r8) :: g0                            ! Imlicit solution: Soil heat flux (W/m2)
    real(r8) :: t0                            ! Imlicit solution: Soil surface temperature (K)
    real(r8) :: e0                            ! Imlicit solution: Soil surface vapor pressure (Pa)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                            ! Aboveground layer index
    real(r8) :: shsrc_leaf                    ! LeafFluxes: Canopy layer leaf sensible heat flux (W/m2)
    real(r8) :: etsrc_leaf                    ! LeafFluxes: Canopy layer leaf water vapor flux (mol H2O/m2/s)
    real(r8) :: stveg_leaf                    ! LeafFluxes: Canopy layer leaf storage heat flux (W/m2)
    !---------------------------------------------------------------------

    associate ( &
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    dpai      => mlcanopy_inst%dpai_profile     , &  ! Canopy layer plant area index (m2/m2)
    fracsun   => mlcanopy_inst%fracsun_profile  , &  ! Canopy layer sunlit fraction (-)
    tleaf     => mlcanopy_inst%tleaf_leaf       , &  ! LeafFluxes: Leaf temperature (K)
    stleaf    => mlcanopy_inst%stleaf_leaf      , &  ! LeafFluxes: Leaf storage heat flux (W/m2 leaf)
    shleaf    => mlcanopy_inst%shleaf_leaf      , &  ! LeafFluxes: Leaf sensible heat flux (W/m2 leaf)
    trleaf    => mlcanopy_inst%trleaf_leaf      , &  ! LeafFluxes: Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf    => mlcanopy_inst%evleaf_leaf      , &  ! LeafFluxes: Leaf evaporation flux (mol H2O/m2 leaf/s)
    shsoi     => mlcanopy_inst%shsoi_soil       , &  ! SoilFluxes: Sensible heat flux: ground (W/m2)
    etsoi     => mlcanopy_inst%etsoi_soil       , &  ! SoilFluxes: Water vapor flux: ground (mol H2O/m2/s)
    gsoi      => mlcanopy_inst%gsoi_soil        , &  ! SoilFluxes: Soil heat flux (W/m2)
    eg        => mlcanopy_inst%eg_soil          , &  ! SoilFluxes: Soil surface vapor pressure (Pa)
    tg        => mlcanopy_inst%tg_soil            &  ! SoilFluxes: Soil surface temperature (K)
    )

    ! Compare leaf fluxes

    do ic = 1, ncan(p)
       if (dpai(p,ic) > 0._r8) then

          ! Layer fluxes as calculated by LeafFluxes

          shsrc_leaf = (shleaf(p,ic,isun) * fracsun(p,ic) + shleaf(p,ic,isha) * (1._r8 - fracsun(p,ic))) * dpai(p,ic)
          etsrc_leaf = (trleaf(p,ic,isun) + evleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic) &
                     + (trleaf(p,ic,isha) + evleaf(p,ic,isha)) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
          stveg_leaf = (stleaf(p,ic,isun) * fracsun(p,ic) + stleaf(p,ic,isha) * (1._r8 - fracsun(p,ic))) * dpai(p,ic)

          ! Error checks

          if (abs(shsrc(ic)-shsrc_leaf) > 0.001_r8) then
             call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Leaf sensible heat flux error')
          end if

          if (abs(lambda*(etsrc(ic)-etsrc_leaf)) > 0.001_r8) then
             call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Leaf latent heat flux error')
          end if

          if (abs(stveg(ic)-stveg_leaf) > 0.001_r8) then
             call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Leaf heat storage error')
          end if

          if (abs(tleaf(p,ic,isun)-tleaf_implic(ic,isun)) > 1.e-06_r8) then
             call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Leaf temperature error (sunlit)')
          end if

          if (abs(tleaf(p,ic,isha)-tleaf_implic(ic,isha)) > 1.e-06_r8) then
             call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Leaf temperature error (shaded)')
          end if
       end if
    end do

    ! Compare soil fluxes

    if (abs(shsoi(p)-sh0) > 0.001_r8) then
       call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Soil sensible heat flux error')
    end if

    if (abs(lambda*(etsoi(p)-et0)) > 0.001_r8) then
       call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Soil latent heat flux error')
    end if

    if (abs(gsoi(p)-g0) > 0.001_r8) then
       call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Soil heat flux error')
    end if

    if (abs(tg(p)-t0) > 1.e-06_r8) then
       call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Soil surface temperature error')
    end if

    if (abs(eg(p)-e0) > 1.e-06_r8) then
       call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Soil surface vapor pressure error')
    end if

    end associate
  end subroutine ErrorCheck01

  !-----------------------------------------------------------------------
  subroutine ErrorCheck02 (p, lambda, avail_energy, shsrc, etsrc, stveg, storage_sh, storage_et, &
  sh0, et0, g0, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Conservation error checks
    !
    ! !USES:
    use MLclm_varpar, only : isun, isha, nlevmlcan, nleaf
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p                  ! Patch index for CLM g/l/c/p hierarchy
    real(r8) :: lambda                        ! Latent heat of vaporization (J/mol)
    real(r8) :: avail_energy(nlevmlcan,nleaf) ! Available energy for leaf (W/m2)
    real(r8) :: shsrc(nlevmlcan)              ! Canopy layer leaf sensible heat flux (W/m2)
    real(r8) :: etsrc(nlevmlcan)              ! Canopy layer leaf water vapor flux (mol H2O/m2/s)
    real(r8) :: stveg(nlevmlcan)              ! Canopy layer leaf storage heat flux (W/m2)
    real(r8) :: storage_sh(nlevmlcan)         ! Heat storage flux in air (W/m2)
    real(r8) :: storage_et(nlevmlcan)         ! Water vapor storage flux in air (mol H2O/m2/s)
    real(r8) :: sh0                           ! Ground sensible heat flux (W/m2)
    real(r8) :: et0                           ! Ground evaporation flux (mol H2O/m2/s)
    real(r8) :: g0                            ! Soil heat flux (W/m2)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                            ! Aboveground layer index
    real(r8) :: err                           ! Energy imbalance (W/m2)
    real(r8) :: sum_src                       ! Sum of source flux over all layers
    real(r8) :: sum_storage                   ! Sum of storage flux over all layers
    !---------------------------------------------------------------------

    associate ( &
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    ntop      => mlcanopy_inst%ntop_canopy      , &  ! Index for top leaf layer
    rnsoi     => mlcanopy_inst%rnsoi_soil       , &  ! Net radiation: ground (W/m2)
    shair     => mlcanopy_inst%shair_profile    , &  ! Canopy layer air sensible heat flux (W/m2)
    etair     => mlcanopy_inst%etair_profile      &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    )

    ! Vegetation flux energy balance

    do ic = 1, ncan(p)
       err = avail_energy(ic,isun) + avail_energy(ic,isha) - shsrc(ic) - lambda * etsrc(ic) - stveg(ic)
       if (abs(err) > 0.001_r8) then
          call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Leaf energy balance error')
       end if
    end do

    ! Flux conservation at each layer

    do ic = 1, ncan(p)
       if (ic == 1) then
          err = storage_sh(ic) - (sh0 + shsrc(ic) - shair(p,ic))
       else
          err = storage_sh(ic) - (shair(p,ic-1) + shsrc(ic) - shair(p,ic))
       end if
       if (abs(err) > 0.001_r8) then
          call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Sensible heat layer conservation error')
       end if

       if (ic == 1) then
          err = storage_et(ic) - (et0 + etsrc(ic) - etair(p,ic))
       else
          err = storage_et(ic) - (etair(p,ic-1) + etsrc(ic) - etair(p,ic))
       end if
       err = err * lambda
       if (abs(err) > 0.001_r8) then
          call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Latent heat layer conservation error')
       end if
    end do

    ! Flux conservation for canopy sensible heat and latent heat. This is to
    ! check canopy conservation equation (so the sum is to ntop not ncan).

    sum_src = 0._r8 ; sum_storage = 0._r8
    do ic = 1, ntop(p)
       sum_src = sum_src + shsrc(ic)
       sum_storage = sum_storage + storage_sh(ic)
    end do

    err = (sh0 + sum_src - sum_storage) - shair(p,ntop(p))
    if (abs(err) > 0.001_r8) then
       call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Sensible heat canopy conservation error')
    end if

    sum_src = 0._r8 ; sum_storage = 0._r8
    do ic = 1, ntop(p)
       sum_src = sum_src + etsrc(ic)
       sum_storage = sum_storage + storage_et(ic)
    end do

    err = (et0 + sum_src - sum_storage) - etair(p,ntop(p))  ! mol H2O/m2/s
    err = err * lambda                                      ! W/m2
    if (abs(err) > 0.001_r8) then
       call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Latent heat canopy conservation error')
    end if

    ! Ground energy balance conservation

    err = rnsoi(p) - sh0 - lambda * et0 - g0
    if (abs(err) > 0.001_r8) then
       call endrun (msg=' ERROR: ImplicitFluxProfileSolution: Ground temperature energy balance error')
    end if

    end associate
  end subroutine ErrorCheck02

  !-----------------------------------------------------------------------
  subroutine WellMixed (num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Canopy scalar profiles equal reference height values
    ! (well-mixed assumption) or are read in from dataset
    !
    ! !USES:
    use MLclm_varctl, only : flux_profile_type
    use MLclm_varpar, only : isun, isha
    use MLLeafFluxesMod, only : LeafFluxes
    use MLSoilFluxesMod, only : SoilFluxes
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer :: fp                       ! Filter index
    integer :: p                        ! Patch index for CLM g/l/c/p hierarchy
    integer :: ic                       ! Aboveground layer index
    !---------------------------------------------------------------------

    associate ( &
                                                     ! *** Input ***
    uref      => mlcanopy_inst%uref_forcing     , &  ! Wind speed at reference height (m/s)
    tref      => mlcanopy_inst%tref_forcing     , &  ! Air temperature at reference height (K)
    eref      => mlcanopy_inst%eref_forcing     , &  ! Vapor pressure at reference height (Pa)
    co2ref    => mlcanopy_inst%co2ref_forcing   , &  ! Atmospheric CO2 at reference height (umol/mol)
    qref      => mlcanopy_inst%qref_forcing     , &  ! Specific humidity at reference height (kg/kg)
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    wind_data => mlcanopy_inst%wind_data_profile, &  ! Canopy layer wind speed FROM DATASET (m/s)
    tair_data => mlcanopy_inst%tair_data_profile, &  ! Canopy layer air temperature FROM DATASET (K)
    eair_data => mlcanopy_inst%eair_data_profile, &  ! Canopy layer vapor pressure FROM DATASET (Pa)
                                                     ! *** Output ***
    uaf       => mlcanopy_inst%uaf_canopy       , &  ! Wind speed at canopy top (m/s)
    taf       => mlcanopy_inst%taf_canopy       , &  ! Air temperature at canopy top for ObuFunc (K)
    qaf       => mlcanopy_inst%qaf_canopy       , &  ! Specific humidity at canopy top for ObuFunc (kg/kg)
    gac0      => mlcanopy_inst%gac0_soil        , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    wind      => mlcanopy_inst%wind_profile     , &  ! Canopy layer wind speed (m/s)
    tair      => mlcanopy_inst%tair_profile     , &  ! Canopy layer air temperature (K)
    eair      => mlcanopy_inst%eair_profile     , &  ! Canopy layer vapor pressure (Pa)
    cair      => mlcanopy_inst%cair_profile     , &  ! Canopy layer atmospheric CO2 (umol/mol)
    gac       => mlcanopy_inst%gac_profile      , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    shair     => mlcanopy_inst%shair_profile    , &  ! Canopy layer air sensible heat flux (W/m2)
    etair     => mlcanopy_inst%etair_profile    , &  ! Canopy layer air water vapor flux (mol H2O/m2/s)
    stair     => mlcanopy_inst%stair_profile    , &  ! Canopy layer air storage heat flux (W/m2)
    mflx      => mlcanopy_inst%mflx_profile     , &  ! Canopy layer momentum flux (m2/s2)
    ! From LeafFluxes
    tleaf     => mlcanopy_inst%tleaf_leaf       , &  ! Leaf temperature (K)
    stleaf    => mlcanopy_inst%stleaf_leaf      , &  ! Leaf storage heat flux (W/m2 leaf)
    shleaf    => mlcanopy_inst%shleaf_leaf      , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf    => mlcanopy_inst%lhleaf_leaf      , &  ! Leaf latent heat flux (W/m2 leaf)
    trleaf    => mlcanopy_inst%trleaf_leaf      , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf    => mlcanopy_inst%evleaf_leaf      , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    ! From SoilFluxes
    shsoi     => mlcanopy_inst%shsoi_soil       , &  ! Sensible heat flux: ground (W/m2)
    lhsoi     => mlcanopy_inst%lhsoi_soil       , &  ! Latent heat flux: ground (W/m2)
    etsoi     => mlcanopy_inst%etsoi_soil       , &  ! Water vapor flux: ground (mol H2O/m2/s)
    gsoi      => mlcanopy_inst%gsoi_soil        , &  ! Soil heat flux (W/m2)
    eg        => mlcanopy_inst%eg_soil          , &  ! Soil surface vapor pressure (Pa)
    tg        => mlcanopy_inst%tg_soil            &  ! Soil surface temperature (K)
    )

    do fp = 1, num_filter
       p = filter(fp)
       do ic = 1, ncan(p)

          ! Scalar profiles

          cair(p,ic) = co2ref(p)

          select case (flux_profile_type)

          ! Use the well-mixed assumption

          case (0)
             wind(p,ic) = uref(p)
             tair(p,ic) = tref(p)
             eair(p,ic) = eref(p)

          ! Use profiles from dataset

          case (-1)
             wind(p,ic) = wind_data(p,ic)
             tair(p,ic) = tair_data(p,ic)
             eair(p,ic) = eair_data(p,ic)
             ! Set each profile individually to WMA if desired
             ! wind(p,ic) = uref(p)
             ! tair(p,ic) = tref(p)
             ! eair(p,ic) = eref(p)

          end select

          ! Verticcal fluxes

          shair(p,ic) = 0._r8
          etair(p,ic) = 0._r8
          stair(p,ic) = 0._r8
          mflx(p,ic) = 0._r8

          ! Calculate leaf fluxes (per unit leaf area)

          call LeafFluxes (p, ic, isun, mlcanopy_inst)
          call LeafFluxes (p, ic, isha, mlcanopy_inst)

       end do

       ! Calculate soil fluxes, but need gac0. Use a large resistance
       ! to approximate bare ground and so that soil fluxes are small.

       gac0(p) = (1._r8 / 100._r8) * 42.3_r8

       call SoilFluxes (p, mlcanopy_inst)

       ! Only needed for output files. These cannot be zero
       ! for analysis package to work.

       uaf(p) = uref(p)
       taf(p) = tref(p)
       qaf(p) = qref(p)
       do ic = 1, ncan(p)
          gac(p,ic) = (1._r8 / 10._r8) * 42.3_r8  ! small non-zero resistance
       end do

    end do

    end associate
  end subroutine WellMixed

end module MLFluxProfileSolutionMod
