module MLCanopyTurbulenceMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Canopy turbulence parameterization
  !
  ! !USES:
  use abortutils, only : endrun
  use clm_varctl, only : iulog, rslfile
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyTurbulence        ! Main routine for canopy turbulence parameterization
  public :: LookupPsihatINI         ! Initialize the RSL psihat look-up tables
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: HF2008                 ! Harman & Finnigan (2008) roughness sublayer theory
  private :: GetObu                 ! Solve for the Obukhov length using subroutine ObuFunc
  private :: ObuFunc                ! Subroutine to solve for the Obukhov length
  private :: GetBeta                ! beta = u* / u at canopy top
  private :: GetPrSc                ! Prandlt number (Pr) and Schmidt number (Sc) at canopy top
  private :: GetPsiRSL              ! RSL-modified stability functions
  private :: phim_monin_obukhov     ! Monin-Obukhov phi stability function for momentum
  private :: phic_monin_obukhov     ! Monin-Obukhov phi stability function for scalars
  private :: psim_monin_obukhov     ! Monin-Obukhov psi stability function for momentum
  private :: psic_monin_obukhov     ! Monin-Obukhov psi stability function for scalars
  private :: LookupPsihat           ! Determines the RSL function psihat as provided through a look-up table
  private :: RoughnessLength        ! Roughness length for momentum
  private :: WindProfile            ! Wind speed profile
  private :: AerodynamicConductance ! Aerodynamic conductances
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CanopyTurbulence (nstep_ml, num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Canopy turbulence parameterization
    !
    ! !USES:
    use MLclm_varctl, only : turb_type
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nstep_ml     ! Current multilayer timestep number
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    select case (turb_type)
    case (1)

       ! Use Harman & Finnigan (2008) roughness sublayer theory

       call HF2008 (nstep_ml, num_filter, filter, mlcanopy_inst)

    case default

       call endrun (msg=' ERROR: CanopyTurbulence: turb_type not valid')

    end select

  end subroutine CanopyTurbulence

  !-----------------------------------------------------------------------
  subroutine HF2008 (nstep_ml, num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Canopy turbulence, wind speed, and aerodynamic conductances using the
    ! Harman and Finnigan (2008) roughness sublayer (RSL) parameterization 
    !
    ! See Bonan et al. (2018) Geosci. Model Dev., 11, 1467-1496, doi:10.5194/gmd-11-1467-2018
    ! with revisions by Bonan et al. (2025)
    !
    ! !USES:
    use clm_time_manager, only : get_nstep
    use MLclm_varcon, only : mmh2o, mmdry, cd, eta_max
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nstep_ml     ! Current multilayer timestep number
    integer, intent(in) :: num_filter   ! Number of patches in filter
    integer, intent(in) :: filter(:)    ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                      ! Aboveground layer index
    real(r8) :: lm                      ! Length scale for momentum (m)
    real(r8) :: lm_over_beta            ! lm / beta
    real(r8) :: eta                     ! eta parameter for wind profile
    integer  :: nstep                   ! Current CLM timestep number
    !---------------------------------------------------------------------

    associate ( &
                                                     ! *** Input ***
    pref      => mlcanopy_inst%pref_forcing     , &  ! Air pressure at reference height (Pa)
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    ntop      => mlcanopy_inst%ntop_canopy      , &  ! Index for top leaf layer
    ztop      => mlcanopy_inst%ztop_canopy      , &  ! Canopy foliage top height (m)
    lai       => mlcanopy_inst%lai_canopy       , &  ! Leaf area index of canopy (m2/m2)
    sai       => mlcanopy_inst%sai_canopy       , &  ! Stem area index of canopy (m2/m2)
    zw        => mlcanopy_inst%zw_profile       , &  ! Canopy height at interface between two adjacent layers (m)
    tair      => mlcanopy_inst%tair_profile     , &  ! Canopy layer air temperature (K)
    eair      => mlcanopy_inst%eair_profile     , &  ! Canopy layer vapor pressure (Pa)
                                                     ! *** Output ***
    ! From HF2008
    Lc        => mlcanopy_inst%Lc_canopy        , &  ! Canopy density length scale (m)
    taf       => mlcanopy_inst%taf_canopy       , &  ! Air temperature at canopy top for ObuFunc (K)
    qaf       => mlcanopy_inst%qaf_canopy       , &  ! Specific humidity at canopy top for ObuFunc (kg/kg)
    mflx      => mlcanopy_inst%mflx_profile     , &  ! Canopy layer momentum flux (m2/s2)
    ! From GetObu
    zdisp     => mlcanopy_inst%zdisp_canopy     , &  ! Displacement height (m)
    beta      => mlcanopy_inst%beta_canopy      , &  ! Value of u* / u at canopy top (-)
    PrSc      => mlcanopy_inst%PrSc_canopy      , &  ! Prandtl (Schmidt) number at canopy top (-)
    ustar     => mlcanopy_inst%ustar_canopy     , &  ! Friction velocity (m/s)
    gac_to_hc => mlcanopy_inst%gac_to_hc_canopy , &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    obu       => mlcanopy_inst%obu_canopy       , &  ! Obukhov length (m)
    ! From RoughnessLength
    z0m       => mlcanopy_inst%z0m_canopy       , &  ! Roughness length for momentum (m)
    ! From WindProfile
    uaf       => mlcanopy_inst%uaf_canopy       , &  ! Wind speed at canopy top (m/s)
    wind      => mlcanopy_inst%wind_profile     , &  ! Canopy layer wind speed (m/s)
    ! From AerodynamicConductance
    gac0      => mlcanopy_inst%gac0_soil        , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    gac       => mlcanopy_inst%gac_profile      , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    kc_eddy   => mlcanopy_inst%kc_eddy_profile    &  ! Canopy layer eddy diffusivity from Harman and Finnigan (m2/s)
    )

    ! Get current step counter (nstep)

    nstep = get_nstep()

    do fp = 1, num_filter
       p = filter(fp)

       ! Canopy density length scale

       Lc(p) = ztop(p) / (cd * (lai(p) + sai(p)))

       ! Temperature and vapor pressure at top of canopy (used in ObuFunc)

       taf(p) = tair(p,ntop(p))
       qaf(p) = mmh2o/mmdry * eair(p,ntop(p)) / (pref(p) - (1._r8-mmh2o/mmdry) * eair(p,ntop(p)))

       ! Calculate Obukhov length

       call GetObu (p, mlcanopy_inst)

       ! Calculate roughness length for momentum

       call RoughnessLength (p, mlcanopy_inst)

       ! The roughness sublayer parameterization uses the expression lm / beta
       ! to calculate wind speed and conductances. Use a constrained value for
       ! lm / beta based on maximum value for eta. eta is the parameter in the
       ! equivalent expression for wind speed: u(zs)=u(ztop)*exp(-eta(1-zs/ztop))
       !
       ! See Bonan (2019), Climate Change and Terrestial Ecosystem Modeling, eqs. (16.8), (16.14)
       !
       ! Note: The current code issues a warning rather than set a maximum value

       lm = 2._r8 * beta(p)**3 * Lc(p)
!      lm_over_beta = lm / beta(p)

       eta = beta(p) / lm * ztop(p)
       if (eta >= eta_max) then
          write (iulog,*) ' Warning: HF2008: lm/beta error'
          write (iulog,*) ' nstep = ', nstep, ' nstep_ml = ', nstep_ml
          write (iulog,*) ' eta = ', eta
       endif 
!      eta = min (eta, eta_max)
       lm_over_beta = ztop(p) / eta

       ! Wind speed profile

       call WindProfile (p, lm_over_beta, mlcanopy_inst)

       ! Momentum flux profile: profile is defined at zw (similar to vertical fluxes)

       do ic = 1, ncan(p)
          if (zw(p,ic) > ztop(p)) then
             ! above canopy
             mflx(p,ic) = -ustar(p)**2
          else
             ! within canopy
             mflx(p,ic) = -(ustar(p)**2) * exp(2._r8*(zw(p,ic) - ztop(p)) / lm_over_beta)
          end if
       end do

       ! Aerodynamic conductances

       call AerodynamicConductance (p, lm_over_beta, mlcanopy_inst)

    end do

    end associate
  end subroutine HF2008

  !-----------------------------------------------------------------------
  subroutine GetObu (p, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Solve for the Obukhov length using subroutine ObuFunc
    !
    ! !USES:
    use MLMathToolsMod, only : hybrid, bisection
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: p               ! Patch index for CLM g/l/c/p hierarchy
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                          ! Aboveground layer index
    integer  :: il                          ! Sunlit (1) or shaded (2) leaf index
    real(r8) :: obu0, obu1                  ! Initial estimates for Obukhov length (m)
    real(r8) :: tol                         ! Accuracy tolerance for Obukhov length (m)
    real(r8) :: dummy                       ! Dummy argument
    !---------------------------------------------------------------------

    associate (obu => mlcanopy_inst%obu_canopy) ! Obukhov length (m)

    ! These are not used, but are needed to pass into the hybrid root solver

    ic = 0 ; il = 0

    ! Calculate Obukhov length (obu) using the subroutine ObuFunc to iterate until
    ! the change in obu is < tol. Do not use the returned value (dummy), and
    ! instead use the value used to calculate ustar

    obu0 = 100._r8          ! Initial estimate for Obukhov length (m)
    obu1 = -100._r8         ! Initial estimate for Obukhov length (m)
    tol = 0.1_r8            ! Accuracy tolerance for Obukhov length (m)

    dummy = hybrid ('GetObu', p, ic, il, mlcanopy_inst, ObuFunc, obu0, obu1, tol)

    ! Code for bisection

!   obu0 = 10000000._r8; obu1 = -10000000._r8
!   dummy = bisection ('GetObu', p, ic, il, mlcanopy_inst, ObuFunc, obu0, obu1, tol)

    ! Set Obukhov length to a near-neutral value

!   obu(p) = -1000._r8
!   call ObuFunc (p, ic, il, mlcanopy_inst, obu(p), dummy)

    end associate
  end subroutine GetObu

  !-----------------------------------------------------------------------
  subroutine ObuFunc (p, ic, il, mlcanopy_inst, obu_val, obu_dif)
    !
    ! !DESCRIPTION:
    ! Solve for the Obukhov length. For the current estimate of the Obukhov length
    ! (obu_val), calculate u*, T*, and q* and then the new length (obu). The subroutine
    ! returns the change in Obukhov length (obu_dif), which equals zero when the
    ! Obukhov length does not change value between iterations.
    !
    ! !USES:
    use clm_varcon, only : grav, vkc
    use MLclm_varctl, only : sparse_canopy_type
    use MLclm_varcon, only : beta_neutral_max, cr, z0mg, LcL_min, LcL_max, aH12
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: p               ! Patch index for CLM g/l/c/p hierarchy
    integer, intent(in)  :: ic              ! Aboveground layer index
    integer, intent(in)  :: il              ! Sunlit (1) or shaded (2) leaf index
    real(r8), intent(in) :: obu_val         ! Input value for Obukhov length (m)
    real(r8), intent(out) :: obu_dif        ! Difference in Obukhov length (m)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: obu_cur                     ! Current value for Obukhov length
    real(r8) :: obu_new                     ! New value for Obukhov length
    real(r8) :: c1                          ! Parameter to calculate beta_neutral
    real(r8) :: beta_neutral                ! Neutral value for beta = u*/u(h)
    real(r8) :: hc_minus_d                  ! Displacement height below canopy top (ztop - zdisp)
    real(r8) :: psim, psic                  ! psi functions for momentum and scalars
    real(r8) :: dum1, dum2                  ! dummy variables (not used)
    real(r8) :: zlog                        ! log((zref-zdisp)/(ztop-zdisp))
    real(r8) :: tstar                       ! Temperature scale (K)
    real(r8) :: qstar                       ! Water vapor scale (kg/kg)
    real(r8) :: tvstar                      ! Virtual potential temperature scale (K)
    real(r8) :: obu_min_stable              ! Minimum value (stable) for Obukhov length (obu)
    real(r8) :: obu_max_unstable            ! Maximum value (unstable) for Obukhov length (obu)
    real(r8) :: LcL                         ! Canopy density scale (Lc) / Obukhov length (obu)
    real(r8) :: beta_HF                     ! RSL: value of u* / u at canopy top (-)
    real(r8) :: beta_norsl                  ! No RSL: value of u* / u at canopy top (-)
    !---------------------------------------------------------------------

    associate ( &
                                                    ! *** Input ***
    zref      => mlcanopy_inst%zref_forcing    , &  ! Atmospheric reference height (m)
    uref      => mlcanopy_inst%uref_forcing    , &  ! Wind speed at reference height (m/s)
    thref     => mlcanopy_inst%thref_forcing   , &  ! Atmospheric potential temperature at reference height (K)
    thvref    => mlcanopy_inst%thvref_forcing  , &  ! Atmospheric virtual potential temperature at reference height (K)
    qref      => mlcanopy_inst%qref_forcing    , &  ! Specific humidity at reference height (kg/kg)
    rhomol    => mlcanopy_inst%rhomol_forcing  , &  ! Molar density at reference height (mol/m3)
    ztop      => mlcanopy_inst%ztop_canopy     , &  ! Canopy foliage top height (m)
    lai       => mlcanopy_inst%lai_canopy      , &  ! Leaf area index of canopy (m2/m2)
    sai       => mlcanopy_inst%sai_canopy      , &  ! Stem area index of canopy (m2/m2)
    Lc        => mlcanopy_inst%Lc_canopy       , &  ! Canopy density length scale (m)
    taf       => mlcanopy_inst%taf_canopy      , &  ! Air temperature at canopy top for ObuFunc (K)
    qaf       => mlcanopy_inst%qaf_canopy      , &  ! Specific humidity at canopy top for ObuFunc (kg/kg)
                                                    ! *** Output ***
    zdisp     => mlcanopy_inst%zdisp_canopy    , &  ! Displacement height (m)
    beta      => mlcanopy_inst%beta_canopy     , &  ! Value of u* / u at canopy top (-)
    PrSc      => mlcanopy_inst%PrSc_canopy     , &  ! Prandtl (Schmidt) number at canopy top (-)
    ustar     => mlcanopy_inst%ustar_canopy    , &  ! Friction velocity (m/s)
    gac_to_hc => mlcanopy_inst%gac_to_hc_canopy, &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    obu       => mlcanopy_inst%obu_canopy        &  ! Obukhov length (m)
    )

    ! Use this current value of Obukhov length

    obu_cur = obu_val

    ! Limit Obukhov length so that values approaching zero are excluded. Criteria uses LcL

    obu_min_stable = Lc(p) / LcL_max
    obu_max_unstable = Lc(p) / LcL_min
    if (obu_cur >= 0._r8) then
       obu_cur = max(obu_cur,obu_min_stable)
    else
       obu_cur = min(obu_cur,obu_max_unstable)
    end if
    LcL = Lc(p) / obu_cur

    ! Determine neutral value of beta
    ! See Bonan et al. (2018), eqs. (A31)-(A32)

    c1 = (vkc / log((ztop(p) + z0mg)/z0mg))**2
    beta_neutral = min(sqrt(c1 + cr*(lai(p)+sai(p))), beta_neutral_max)

    ! Calculate beta = u* / u(h) for current Obukhov length.
    ! Use the parameterization of Harman (2012, eq. 13), which
    ! adds in the value of beta with no RSL influence.

    call GetBeta (beta_neutral, LcL, beta_HF)
    call GetBeta (vkc/2._r8, LcL, beta_norsl)
    if (LcL > aH12(2)) then
       beta(p) = beta_HF
    else
       beta(p) = beta_norsl + (beta_HF - beta_norsl) / (1 + aH12(1) * abs(LcL - aH12(2))**aH12(3))
    end if

    ! Displacement height, and then adjust for canopy sparseness
    ! See Bonan et al. (2018), eqs. (A28), (A33)

    hc_minus_d = beta(p)**2 * Lc(p)
    select case (sparse_canopy_type)
    case (1)
       hc_minus_d = hc_minus_d * (1._r8 - exp(-0.25_r8*(lai(p)+sai(p))/beta(p)**2))
    end select
    hc_minus_d = min(ztop(p), hc_minus_d)
    zdisp(p) = ztop(p) - hc_minus_d

    if ((zref(p) - zdisp(p)) < 0._r8) then
       call endrun (msg=' ERROR: ObuFunc: zdisp height > zref')
    end if

    ! Calculate Prandlt number (Pr) and Schmidt number (Sc) at canopy
    ! top for current Obukhov length. Only need Pr because Sc = Pr.

    call GetPrSc (beta_neutral, beta_neutral_max, LcL, PrSc(p))

    ! Calculate the stability functions psi for momentum and scalars. The
    ! returned function values (psim, psic) are the Monin-Obukhov psi functions
    ! and additionally include the roughness sublayer psi_hat functions.
    ! These are evaluated at the reference height and at the canopy height.

    call GetPsiRSL (zref(p), ztop(p), zdisp(p), obu_cur, beta(p), PrSc(p), psim, psic, dum1, dum2)

    ! Friction velocity
    ! See Bonan et al. (2018), eq. (19)

    zlog = log((zref(p)-zdisp(p)) / (ztop(p)-zdisp(p)))
    ustar(p) = uref(p) * vkc / (zlog + psim)

    ! Temperature scale
    ! See Bonan et al. (2018), eq. (20)

    tstar = (thref(p) - taf(p)) * vkc / (zlog + psic)

    ! Water vapor scale - use units of specific humidity (kg/kg)

    qstar = (qref(p) - qaf(p)) * vkc / (zlog + psic)

    ! Aerodynamic conductance to canopy height

    gac_to_hc(p) = rhomol(p) * vkc * ustar(p) / (zlog + psic)

    ! Save value for obu used to calculate ustar

    obu(p) = obu_cur

    ! Calculate new Obukhov length (m)
    ! See Bonan et al. (2018), eqs. (A29), (A30)

    tvstar = tstar + 0.61_r8 * thref(p) * qstar
    obu_new = ustar(p)**2 * thvref(p) / (vkc * grav * tvstar)

    ! Change in Obukhov length (m)

    obu_dif = obu_new - obu_val

    end associate
  end subroutine ObuFunc

  !-----------------------------------------------------------------------
  subroutine GetBeta (beta_neutral, LcL, beta)
    !
    ! !DESCRIPTION:
    ! Calculate beta = u* / u(h) for current Obukhov length
    ! See Bonan et al. (2018), eqs. (A22)-(A24)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: beta_neutral    ! Neutral value for beta = u*/u(h)
    real(r8), intent(in)  :: LcL             ! Canopy density scale (Lc) / Obukhov length (obu)
    real(r8), intent(out) :: beta            ! Value of u*/u(h) at canopy top
    !
    ! !LOCAL VARIABLES:
    real(r8) :: aa, bb, cc, dd, qq, rr       ! Terms for quadratic or cubic solution
    real(r8) :: y, fy, err                   ! Used for error check
    !---------------------------------------------------------------------

    if (LcL <= 0._r8) then

       ! Unstable case: quadratic equation for beta^2 at LcL

       aa = 1._r8
       bb = 16._r8 * LcL * beta_neutral**4     
       cc = -beta_neutral**4
       beta = sqrt((-bb + sqrt(bb**2 - 4._r8 * aa * cc)) / (2._r8 * aa))

    else

       ! Stable case: cubic equation for beta at LcL

       aa = 5._r8 * LcL
       bb = 0._r8
       cc = 1._r8
       dd = -beta_neutral
       qq = (2._r8*bb**3 - 9._r8*aa*bb*cc + 27._r8*(aa**2)*dd)**2 - 4._r8*(bb**2 - 3._r8*aa*cc)**3
       qq = sqrt(qq)
       rr = 0.5_r8 * (qq + 2._r8*bb**3 - 9._r8*aa*bb*cc + 27._r8*(aa**2)*dd)
       rr = rr**(1._r8/3._r8)
       beta = -(bb+rr)/(3._r8*aa) - (bb**2 - 3._r8*aa*cc)/(3._r8*aa*rr)    

    end if

    ! Error check

    y = LcL * beta**2 
    fy = phim_monin_obukhov(y)
    err = beta * fy - beta_neutral
    if (abs(err) > 1.e-06_r8) then
       call endrun (msg=' ERROR: GetBeta: beta error')
    end if

  end subroutine GetBeta

  !-----------------------------------------------------------------------
  subroutine GetPrSc (beta_neutral, beta_neutral_max, LcL, PrSc)
    !
    ! !DESCRIPTION:
    ! Calculate Prandlt number (Pr) and Schmidt number (Sc) at canopy
    ! top for current Obukhov length
    ! See Bonan et al. (2018), eqs. (A25), (A34)
    !
    ! !USES:
    use MLclm_varctl, only : sparse_canopy_type
    use MLclm_varcon, only : Pr0, Pr1, Pr2
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: beta_neutral     ! Neutral value for beta = u*/u(h)
    real(r8), intent(in)  :: beta_neutral_max ! Maximum value for beta in neutral conditions
    real(r8), intent(in)  :: LcL              ! Canopy density scale (Lc) / Obukhov length (obu)
    real(r8), intent(out) :: PrSc             ! Prandtl (Schmidt) number at canopy top
    !---------------------------------------------------------------------
    
    PrSc = Pr0 + Pr1 * tanh(Pr2*LcL)

    ! Adjust for canopy sparseness

    select case (sparse_canopy_type)
    case (1)
       PrSc = (1._r8 - beta_neutral/beta_neutral_max) * 1._r8 + (beta_neutral/beta_neutral_max) * PrSc
    end select

  end subroutine GetPrSc

  !-----------------------------------------------------------------------
  subroutine GetPsiRSL (za, hc, disp, obu, beta, PrSc, psim, psic, psim2, psim_hat2)
    !
    ! !DESCRIPTION:
    ! Calculate the stability functions psi for momentum and scalars. The
    ! returned function values (psim, psic) are the Monin-Obukhov psi functions
    ! and additionally include the roughness sublayer psihat functions.
    ! These are evaluated between the height za and the canopy height hc.
    ! See Bonan et al. (2018), appendix A2
    !
    ! !USES:
    use clm_varcon, only : vkc
    use MLclm_varcon, only : c2, dtLgridM, zdtgridM, psigridM, dtLgridH, zdtgridH, psigridH
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: za        ! Atmospheric height (m)
    real(r8), intent(in)  :: hc        ! Canopy height (m)
    real(r8), intent(in)  :: disp      ! Displacement height (m)
    real(r8), intent(in)  :: obu       ! Obukhov length (m)
    real(r8), intent(in)  :: beta      ! Value of u* / u at canopy top
    real(r8), intent(in)  :: PrSc      ! Prandtl (Schmidt) number at canopy top
    real(r8), intent(out) :: psim      ! psi function for momentum including RSL influence
    real(r8), intent(out) :: psic      ! psi function for scalars  including RSL influence
    real(r8), intent(out) :: psim2     ! Monin-Obukhov psi function for momentum evaluated at hc
    real(r8), intent(out) :: psim_hat2 ! RSL psihat function for momentum evaluated at hc
    !
    ! !LOCAL VARIABLES:
    real(r8) :: dt                     ! Displacement height below canopy top (dt = ztop - zdisp)
    real(r8) :: phim                   ! Monin-Obukhov phi function for momentum at canopy top
    real(r8) :: phic                   ! Monin-Obukhov phi function for scalars at canopy top
    real(r8) :: c1                     ! RSL magnitude multiplier
    real(r8) :: psim1, psic1           ! Monin-Obukhov psi functions (momentum, scalars) evaluated at za
    real(r8) :: psic2                  ! Monin-Obukhov psi function (scalars) evaluated at hc
    real(r8) :: psim_hat1, psic_hat1   ! RSL psihat functions (momentum, scalars) evaluated at za
    real(r8) :: psic_hat2              ! RSL psihat function (scalars) evaluated at hc
    !---------------------------------------------------------------------

    ! In the RSL theory, c1 and c2 represent the scaled magnitude and
    ! height over which the RSL theory modifies MOST via the psihat functions:
    !
    !
    !     z
    !     ^
    !     |                                    . ___
    !     |                                    .  ^
    !     |                                   .   |
    !     |                                  .    |
    !     |                                 .     |
    !     |                               .       |
    !     |                             .        c2
    !     |                          .            |
    !     |                       .               |
    !     |                  .                    |
    !     |           .                           |
    !     |   .                                 _\/_
    !     -------------------------------------------> u
    !         |<-------------- c1 --------------->|

    ! See Bonan et al. (2018), eq. (A16)

    phim = phim_monin_obukhov((hc-disp)/obu)
    c1 = (1._r8 - vkc / (2._r8 * beta * phim)) * exp(0.5_r8*c2)

    ! Evaluate the roughness sublayer psihat function for momentum at
    ! the height za and at the canopy height hc. Values for psihat are obtained
    ! from a look-up table. Here, heights are above the canopy for compatibility
    ! with the supplied look-up table. These heights are also scaled by dt = hc-d
    ! so that the look-up table uses (za-hc)/dt and (hc-hc)/dt. Also the term
    ! dt in the integration of psihat is scaled by the Obukhov length L (dt/obu).
    ! This means that the returned psihat value needs to be scaled (multiplied) by
    ! c1 before it fully represents psihat as it appears in the RSL equations.
    ! See Bonan et al. (2018), eq. (A21)

    dt = hc - disp
    call LookupPsihat ((za-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM, psim_hat1)
    call LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridM, dtLgridM, psigridM, psim_hat2)
    psim_hat1 = psim_hat1 * c1
    psim_hat2 = psim_hat2 * c1

    ! Evaluate the Monin-Obukhov psi function for momentum at the height za
    ! and at the canopy height hc

    psim1 = psim_monin_obukhov((za-disp)/obu)
    psim2 = psim_monin_obukhov((hc-disp)/obu)

    ! psi function for momentum including RSL influence
    ! See Bonan et al. (2018), eq. (19)

    psim = -psim1 + psim2 + psim_hat1 - psim_hat2 + vkc / beta

    ! Now do the same for scalars
    ! Corresponding equations are (A19), (A21), and (20)

    phic = phic_monin_obukhov((hc-disp)/obu)
    c1 = (1._r8 - PrSc*vkc / (2._r8 * beta * phic)) * exp(0.5_r8*c2)

    call LookupPsihat ((za-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH, psic_hat1)
    call LookupPsihat ((hc-hc)/dt, dt/obu, zdtgridH, dtLgridH, psigridH, psic_hat2)
    psic_hat1 = psic_hat1 * c1
    psic_hat2 = psic_hat2 * c1

    psic1 = psic_monin_obukhov((za-disp)/obu)
    psic2 = psic_monin_obukhov((hc-disp)/obu)

    psic = -psic1 + psic2 + psic_hat1 - psic_hat2

  end subroutine GetPsiRSL

  !-----------------------------------------------------------------------
  function phim_monin_obukhov (zeta) result(phi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov phi stability function for momentum
    ! See Bonan et al. (2018), eq. (A10)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: phi               ! phi for momentum
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then               ! --- unstable
       phi = 1._r8 / sqrt(sqrt(1._r8 - 16._r8 * zeta))
    else                                 ! --- stable
       phi = 1._r8 + 5._r8 * zeta
    end if

  end function phim_monin_obukhov

  !-----------------------------------------------------------------------
  function phic_monin_obukhov (zeta) result(phi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov phi stability function for scalars
    ! See Bonan et al. (2018), eq. (A11)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: phi               ! phi for scalars
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then               ! --- unstable
       phi = 1._r8 / sqrt(1._r8 - 16._r8 * zeta)
    else                                 ! --- stable
       phi = 1._r8 + 5._r8 * zeta
    end if

  end function phic_monin_obukhov

  !-----------------------------------------------------------------------
  function psim_monin_obukhov (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov psi stability function for momentum
    ! See Bonan et al. (2018), eq. (A12)
    !
    ! !USES:
    use clm_varcon, only : pi => rpi
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x                 ! (1 - 16*zeta)**1/4
    real(r8) :: psi               ! psi for momentum
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then               ! --- unstable
       x = sqrt(sqrt(1._r8 - 16._r8 * zeta))
       psi = 2._r8 * log((1._r8+x)/2._r8) + log((1._r8+x*x)/2._r8) - 2._r8*atan(x) + pi * 0.5_r8
    else                                 ! --- stable
       psi = -5._r8 * zeta
    end if

  end function psim_monin_obukhov

  !-----------------------------------------------------------------------
  function psic_monin_obukhov (zeta) result(psi)
    !
    ! !DESCRIPTION:
    ! Monin-Obukhov psi stability function for scalars
    ! See Bonan et al. (2018), eq. (A13)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zeta  ! Monin-Obukhov stability parameter
    !
    ! !LOCAL VARIABLES:
    real(r8) :: x                 ! (1 - 16*zeta)**1/4
    real(r8) :: psi               ! psi for scalars
    !---------------------------------------------------------------------

    if (zeta < 0._r8) then               ! --- unstable
       x = sqrt(sqrt(1._r8 - 16._r8 * zeta))
       psi = 2._r8 * log((1._r8+x*x)/2._r8)
    else                                 ! --- stable
       psi = -5._r8 * zeta
    end if

  end function psic_monin_obukhov

  !-----------------------------------------------------------------------
  subroutine LookupPsihat (zdt, dtL, zdtgrid, dtLgrid, psigrid, psihat)
    !
    ! !DESCRIPTION:
    ! Evaluate the RSL psihat function as provided through a look-up table for
    ! input values of zdt = (z-h)/(h-d) and dtL = (h-d)/L. Linearly interpolate
    ! between values supplied on the look-up table grid defined by zdtgrid, dtLgrid,
    ! and psigrid.
    !
    ! NOTE: The psihat presented in Harman and Finnigan (2007,2008) and Harman (2012)
    ! has been re-written in non-dimensional form such that it now appears as:
    !
    ! psihat(z) = c1 * A[(z-h)/(h-d),(h-d)/L]
    !
    ! This routine gets the value of A from a look-up table. The returned psihat
    ! value needs to be scaled (multiplied) by c1 before it fully represents
    ! psihat as it appears in the RSL equations.
    !
    ! NOTE: dt = beta^2*Lc = h-d
    !
    ! See Bonan et al. (2018), eq. (A21)
    !
    ! !USES:
    use MLclm_varcon, only : nZ, nL
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: zdt            ! Height above canopy (z-h) normalized by displacment height below canopy (h-d)
    real(r8), intent(in) :: dtL            ! Displacement height below canopy (h-d) / Obukhov length
    real(r8), intent(in) :: zdtgrid(nZ,1)  ! Grid of zdt on which psihat is given
    real(r8), intent(in) :: dtLgrid(1,nL)  ! Grid of dtL on which psihat is given
    real(r8), intent(in) :: psigrid(nZ,nL) ! Grid of psihat values
    real(r8), intent(out):: psihat         ! Returned value for psihat
    !
    !LOCAL VARIABLES
    integer  :: ii, jj                     ! Looping indices
    integer  :: L1, L2, Z1, Z2             ! Grid indices for psihat sought
    real(r8) :: wL1, wL2, wZ1, wZ2         ! Weights for averaging
    !---------------------------------------------------------------------

    ! Find indices and weights for dtL values which bracket the specified dtL

    L1 = 0 ; L2 = 0
    if (dtL <= dtLgrid(1,1)) then
       L1 = 1
       L2 = 1
       wL1 = 0.5_r8
       wL2 = 0.5_r8
    else if (dtL > dtLgrid(1,nL)) then
       L1 = nL
       L2 = nL
       wL1 = 0.5_r8
       wL2 = 0.5_r8
    else
       do jj = 1, nL-1
          if ((dtL <= dtLgrid(1,jj+1)) .and. (dtL > dtLgrid(1,jj))) then
             L1 = jj
             L2 = jj + 1
             wL1 = (dtLgrid(1,L2) - dtL) / (dtLgrid(1,L2) - dtLgrid(1,L1))
             wL2 = 1._r8 - wL1
          end if
       end do
    end if

    if (L1 == 0 .or. L2 == 0) then
       call endrun (msg=' ERROR: LookupPsihat error: indices L1 and L2 not found')
    end if

    ! Find indices and weights for zdt values which bracket the specified zdt

    Z1 = 0 ; Z2 = 0
    if (zdt > zdtgrid(1,1)) then
       Z1 = 1
       Z2 = 1
       wZ1 = 0.5_r8
       wZ2 = 0.5_r8
    else if (zdt < zdtgrid(nZ,1)) then
       Z1 = nZ
       Z2 = nZ
       wZ1 = 0.5_r8
       wZ2 = 0.5_r8
    else
       do ii = 1, nZ-1
          if ((zdt >= zdtgrid(ii+1,1)) .and. (zdt < zdtgrid(ii,1))) then
             Z1 = ii
             Z2 = ii + 1
             wZ1 = (zdt - zdtgrid(ii+1,1)) / (zdtgrid(ii,1) - zdtgrid(ii+1,1))
             wZ2 = 1._r8 - wZ1
          end if
       end do
    end if

    if (Z1 == 0 .or. Z2 == 0) then
       call endrun (msg=' ERROR: LookupPsihat error: indices Z1 and Z2 not found')
    end if

    ! Calculate psihat as a weighted average of the values of psihat on the grid

    psihat = wZ1 * wL1 * psigrid(Z1,L1) + wZ2 * wL1 * psigrid(Z2,L1) &
           + wZ1 * wL2 * psigrid(Z1,L2) + wZ2 * wL2 * psigrid(Z2,L2)

  end subroutine LookupPsihat

  !-----------------------------------------------------------------------
  subroutine RoughnessLength (p, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Calculate roughness length for momentum
    !
    ! !USES:
    use clm_varcon, only : vkc
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p                ! Patch index for CLM g/l/c/p hierarchy
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: hc_minus_d                  ! Displacement height below canopy top (ztop - zdisp)
    real(r8) :: psim, psic                  ! psi functions for momentum and scalars
    real(r8) :: psim_hc                     ! Monin-Obukhov psi function for momentum evaluated at hc
    real(r8) :: psim_hat_hc                 ! RSL psihat function for momentum evaluated at hc
    real(r8) :: exp1, exp2                  ! Exponential terms
    real(r8) :: aval, bval                  ! Bracketed estimates for roughness length (m)
    real(r8) :: cval                        ! Refined roughness length (m)
    real(r8) :: z0m_aval                    ! z0m evaluated with aval
    real(r8) :: z0m_bval                    ! z0m evaluated with bval
    real(r8) :: z0m_cval                    ! z0m evaluated with cval
    real(r8) :: psim_z0m_aval               ! Monin-Obukhov psi function for momentum evaluated at z0m = aval
    real(r8) :: psim_z0m_bval               ! Monin-Obukhov psi function for momentum evaluated at z0m = bval
    real(r8) :: psim_z0m_cval               ! Monin-Obukhov psi function for momentum evaluated at z0m = cval
    real(r8) :: fa, fb, fc                  ! Error (m)
    real(r8) :: err                         ! Error tolerance (m)
    integer  :: n                           ! Iteration counter
    integer  :: nmax                        ! Maximum number of iterations
    !---------------------------------------------------------------------

    associate ( &
                                                    ! *** Input ***
    zref      => mlcanopy_inst%zref_forcing    , &  ! Atmospheric reference height (m)
    ztop      => mlcanopy_inst%ztop_canopy     , &  ! Canopy foliage top height (m)
    zdisp     => mlcanopy_inst%zdisp_canopy    , &  ! Displacement height (m)
    obu       => mlcanopy_inst%obu_canopy      , &  ! Obukhov length (m)
    beta      => mlcanopy_inst%beta_canopy     , &  ! Value of u* / u at canopy top (-)
    PrSc      => mlcanopy_inst%PrSc_canopy     , &  ! Prandtl (Schmidt) number at canopy top (-)
                                                    ! *** Output ***
    z0m       => mlcanopy_inst%z0m_canopy        &  ! Roughness length for momentum (m)
    )

    ! Get psi and psi_hat terms at canopy top

    call GetPsiRSL (zref(p), ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim, psic, psim_hc, psim_hat_hc)

    ! Other common terms

    hc_minus_d = ztop(p) - zdisp(p)
    exp1 = exp(-vkc / beta(p))
    exp2 = exp(psim_hat_hc)

    ! Use bisection to find z0m, which lies between aval and bval, and refine the
    ! estimate until the difference is less than err or after nmax itertions

    aval = ztop(p)
    bval = 0._r8
    err = 0.001_r8
    nmax = 20

    psim_z0m_aval = psim_monin_obukhov(aval/obu(p))
    z0m_aval = hc_minus_d * exp1 * exp(-psim_hc + psim_z0m_aval) * exp2
    fa = z0m_aval - aval

    psim_z0m_bval = psim_monin_obukhov(bval/obu(p))
    z0m_bval = hc_minus_d * exp1 * exp(-psim_hc + psim_z0m_bval) * exp2
    fb = z0m_bval - bval

    if (fa * fb > 0._r8) then
       call endrun (msg=' ERROR: RoughnessLength: bisection error - f(a) and f(b) do not have opposite signs')
    end if

    n = 1
    do while (abs(bval-aval) > err .and. n <= nmax)
       cval = (aval + bval) / 2._r8
       psim_z0m_cval = psim_monin_obukhov(cval/obu(p))
       z0m_cval = hc_minus_d * exp1 * exp(-psim_hc + psim_z0m_cval) * exp2
       fc = z0m_cval - cval
       if (fa * fc < 0._r8) then
          bval = cval; fb = fc
       else
          aval = cval; fa = fc
       end if
       n = n + 1
    end do

    if (n > nmax) then
       call endrun (msg=' ERROR: RoughnessLength: maximum iteration exceeded')
    end if 

    z0m(p) = cval

    end associate
  end subroutine RoughnessLength

  !-----------------------------------------------------------------------
  subroutine WindProfile (p, lm_over_beta, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Wind speed profile above and within canopy
    !
    ! !USES:
    use clm_varcon, only : vkc
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p              ! Patch index for CLM g/l/c/p hierarchy
    real(r8), intent(in) :: lm_over_beta  ! lm / beta
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                        ! Aboveground layer index
    real(r8) :: psim                      ! psi function for momentum
    real(r8) :: psic                      ! psi function for scalars
    real(r8) :: dum1, dum2                ! dummy variables (not used)
    real(r8) :: zlog_m                    ! log height
    !---------------------------------------------------------------------

    associate ( &
                                                 ! *** Input ***
    ncan      => mlcanopy_inst%ncan_canopy  , &  ! Number of aboveground layers
    ntop      => mlcanopy_inst%ntop_canopy  , &  ! Index for top leaf layer
    zs        => mlcanopy_inst%zs_profile   , &  ! Canopy layer height for scalar concentration and source (m)
    ztop      => mlcanopy_inst%ztop_canopy  , &  ! Canopy foliage top height (m)
    zdisp     => mlcanopy_inst%zdisp_canopy , &  ! Displacement height (m)
    obu       => mlcanopy_inst%obu_canopy   , &  ! Obukhov length (m)
    beta      => mlcanopy_inst%beta_canopy  , &  ! Value of u* / u at canopy top (-)
    PrSc      => mlcanopy_inst%PrSc_canopy  , &  ! Prandtl (Schmidt) number at canopy top (-)
    ustar     => mlcanopy_inst%ustar_canopy , &  ! Friction velocity (m/s)
                                                 ! *** Output ***
    uaf      => mlcanopy_inst%uaf_canopy    , &  ! Wind speed at canopy top (m/s)
    wind     => mlcanopy_inst%wind_profile    &  ! Canopy layer wind speed (m/s)
    )

    ! Above-canopy wind profile: wind speed is defined at zs
    ! See Bonan et al. (2018) Geosci. Model Dev., 11, 1467-1496, doi:10.5194/gmd-11-1467-2018, eq. (19)

    do ic = ntop(p)+1, ncan(p)
       call GetPsiRSL (zs(p,ic), ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim, psic, dum1, dum2)
       zlog_m = log((zs(p,ic)-zdisp(p)) / (ztop(p)-zdisp(p)))
       wind(p,ic) = ustar(p) / vkc * (zlog_m + psim)
    end do

    ! Wind speed at top of canopy
    ! See Bonan et al. (2018), eq. (21)

    uaf(p) = ustar(p) / beta(p)

    ! Within-canopy wind profile: wind speed is defined at zs

    do ic = 1, ntop(p)
       wind(p,ic) = uaf(p) * exp((zs(p,ic) - ztop(p)) / lm_over_beta)
    end do

    end associate
  end subroutine WindProfile

  !-----------------------------------------------------------------------
  subroutine AerodynamicConductance (p, lm_over_beta, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Aerodynamic conductances above and within the canopy.
    ! Conductances are defined between zs(i) and zs(i+1).
    !
    ! !USES:
    use clm_varcon, only : vkc
    use MLclm_varctl, only : HF_extension_type
    use MLclm_varcon, only : z0mg, ra_max
    use MLCanopyFluxesType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: p              ! Patch index for CLM g/l/c/p hierarchy
    real(r8), intent(in) :: lm_over_beta  ! lm / beta
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                        ! Aboveground layer index
    real(r8) :: psim1, psim2              ! psi function for momentum
    real(r8) :: psic, psic1, psic2        ! psi function for scalars
    real(r8) :: dum1, dum2                ! dummy variables (not used)
    real(r8) :: zlog_m, zlog_c            ! log height
    real(r8) :: zu, zl                    ! Upper and lower heights for within canopy resistances (m)
    real(r8) :: res                       ! Resistance (s/m)
    real(r8) :: ustar_g                   ! Friction velocity at ground (m/s)
    real(r8) :: z0cg                      ! Roughness length of ground for scalars (m)
    real(r8) :: sumres                    ! Sum of aerodynamic resistances above canopy
    real(r8) :: gac_above_foliage         ! Aerodynamic conductance for top/bottom foliage layer
    real(r8) :: gac_below_foliage         ! Aerodynamic conductance for top/bottom foliage layer
    !---------------------------------------------------------------------

    associate ( &
                                                     ! *** Input ***
    zref      => mlcanopy_inst%zref_forcing     , &  ! Atmospheric reference height (m)
    rhomol    => mlcanopy_inst%rhomol_forcing   , &  ! Molar density at reference height (mol/m3)
    ncan      => mlcanopy_inst%ncan_canopy      , &  ! Number of aboveground layers
    ntop      => mlcanopy_inst%ntop_canopy      , &  ! Index for top leaf layer
    ztop      => mlcanopy_inst%ztop_canopy      , &  ! Canopy foliage top height (m)
    zdisp     => mlcanopy_inst%zdisp_canopy     , &  ! Displacement height (m)
    obu       => mlcanopy_inst%obu_canopy       , &  ! Obukhov length (m)
    beta      => mlcanopy_inst%beta_canopy      , &  ! Value of u* / u at canopy top (-)
    PrSc      => mlcanopy_inst%PrSc_canopy      , &  ! Prandtl (Schmidt) number at canopy top (-)
    ustar     => mlcanopy_inst%ustar_canopy     , &  ! Friction velocity (m/s)
    gac_to_hc => mlcanopy_inst%gac_to_hc_canopy , &  ! Aerodynamic conductance for a scalar above canopy (mol/m2/s)
    zs        => mlcanopy_inst%zs_profile       , &  ! Canopy layer height for scalar concentration and source (m)
    wind      => mlcanopy_inst%wind_profile     , &  ! Canopy layer wind speed (m/s)
                                                     ! *** Output ***
    gac0      => mlcanopy_inst%gac0_soil        , &  ! Aerodynamic conductance for soil fluxes (mol/m2/s)
    gac       => mlcanopy_inst%gac_profile      , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    kc_eddy   => mlcanopy_inst%kc_eddy_profile    &  ! Canopy layer eddy diffusivity from Harman and Finnigan (m2/s)
    )

    ! -------------------------------------
    ! Above-canopy aerodynamic conductances 
    ! conductance: mol/m3 * m/s = mol/m2/s
    ! -------------------------------------

    ! See Bonan et al. (2018) Geosci. Model Dev., 11, 1467-1496, doi:10.5194/gmd-11-1467-2018, eq. (24)

    do ic = ntop(p)+1, ncan(p)-1
       call GetPsiRSL (zs(p,ic),   ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim1, psic1, dum1, dum2)
       call GetPsiRSL (zs(p,ic+1), ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim2, psic2, dum1, dum2)
       ! equivalent to: -psi_c_z2 + psi_c_z1 + psi_c_rsl_z2 - psi_c_rsl_z1
       psic = psic2 - psic1
       zlog_c = log((zs(p,ic+1)-zdisp(p)) / (zs(p,ic)-zdisp(p)))
       gac(p,ic) = rhomol(p) * vkc * ustar(p) / (zlog_c + psic)
    end do

    ! Special case for the top layer to the reference height

    ic = ncan(p)
    call GetPsiRSL (zs(p,ic), ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim1, psic1, dum1, dum2)
    call GetPsiRSL (zref(p),  ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim2, psic2, dum1, dum2)
    psic = psic2 - psic1
    zlog_c = log((zref(p)-zdisp(p)) / (zs(p,ic)-zdisp(p)))
    gac(p,ic) = rhomol(p) * vkc * ustar(p) / (zlog_c + psic)

    ! The top foliage layer includes terms for within and above the canopy. Here,
    ! calculate the conductance from top of foliage at height ztop to the layer
    ! immediately above it at height zs(ntop+1)

    ic = ntop(p)
    call GetPsiRSL (ztop(p),    ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim1, psic1, dum1, dum2)
    call GetPsiRSL (zs(p,ic+1), ztop(p), zdisp(p), obu(p), beta(p), PrSc(p), psim2, psic2, dum1, dum2)
    psic = psic2 - psic1
    zlog_c = log((zs(p,ic+1)-zdisp(p)) / (ztop(p)-zdisp(p)))
    gac_above_foliage = rhomol(p) * vkc * ustar(p) / (zlog_c + psic)

    ! Make sure above-canopy aerodynamic resistances sum to 1/gac_to_hc

    sumres = 1._r8 / gac_above_foliage
    do ic = ntop(p)+1, ncan(p)
       sumres = sumres + 1._r8 / gac(p,ic)
    end do

    if (abs(1._r8/sumres - gac_to_hc(p)) > 1.e-06_r8) then
       call endrun (msg=' ERROR: AerodynamicConductance: above-canopy aerodynamic conductance error')
    end if

    ! --------------------------------------
    ! Within-canopy aerodynamic conductances
    ! res = resistance: s/m
    ! gac = conductance: (mol/m3) / (s/m) = mol/m2/s
    ! See Bonan et al. (2018), eq. (26)
    ! --------------------------------------

    do ic = 1, ntop(p)-1
       zl = zs(p,ic) - ztop(p)
       zu = zs(p,ic+1) - ztop(p)
       res = PrSc(p) / (beta(p) * ustar(p)) * (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
       gac(p,ic) = rhomol(p) / res
    end do

    ! Special case for top foliage layer: conductance from zs to ztop ...

    ic = ntop(p)
    zl = zs(p,ic) - ztop(p)
    zu = ztop(p) - ztop(p)
    res = PrSc(p) / (beta(p) * ustar(p)) * (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
    gac_below_foliage = rhomol(p) / res

    ! ... and now include additional conductance from top of foliage to next layer above

    gac(p,ic) = 1._r8 / (1._r8 / gac_below_foliage + 1._r8 / gac_above_foliage)

    ! Aerodynamic conductance at ground

    z0cg = 0.1_r8 * z0mg
    if (z0mg > zs(p,1) .or. z0cg > zs(p,1)) then
       call endrun (msg=' ERROR: AerodynamicConductance: soil roughness error')
    end if

    select case (HF_extension_type)

    case (1)

       ! Extend HF exponential profile to ground (taken as z0cg)
       ! as in Bonan et al. (2021)

       zl = z0cg - ztop(p)
       zu = zs(p,1) - ztop(p)
       res = PrSc(p) / (beta(p) * ustar(p)) * (exp(-zl/lm_over_beta) - exp(-zu/lm_over_beta))
       gac0(p) = rhomol(p) / res

    case (2)

       ! Use log profile to ground (taken as z0mg, where wind speed is zero)
       ! Note the minimum limit for wind speed at the ground to prevent small gac0
       ! See Bonan et al. (2018), eq. (27)

       zlog_m = log(zs(p,1)/z0mg)
       ustar_g = max(wind(p,1),0.1_r8) * vkc / zlog_m
       gac0(p) = rhomol(p) * vkc * ustar_g / zlog_m

!      zlog_m = log(zs(p,1)/z0mg)                     !!! CLMml v0 CODE !!!
!      ustar_g = wind(p,1) * vkc / zlog_m             !!! CLMml v0 CODE !!!
!      ustar_g = max(ustar_g, 0.01_r8)                !!! CLMml v0 CODE !!!
!      z0cg = 0.1_r8 * z0mg                           !!! CLMml v0 CODE !!!
!      zlog_c = log(zs(p,1)/z0cg)                     !!! CLMml v0 CODE !!!
!      gac0(p) = rhomol(p) * vkc * ustar_g / zlog_c   !!! CLMml v0 CODE !!!

    end select

    ! Limit resistances to < 500 s/m

    res = min (rhomol(p)/gac0(p), ra_max)
    gac0(p) = rhomol(p) / res

    do ic = 1, ncan(p)
       res = min (rhomol(p)/gac(p,ic), ra_max)
       gac(p,ic) = rhomol(p) / res
    end do

    ! Eddy diffusivity diagnosed from conductance

    do ic = 1, ncan(p)
       if (ic == ncan(p)) then
          kc_eddy(p,ic) = gac(p,ic) / rhomol(p) * (zref(p) - zs(p,ic))
       else
          kc_eddy(p,ic) = gac(p,ic) / rhomol(p) * (zs(p,ic+1) - zs(p,ic))
       end if
    end do

    end associate
  end subroutine AerodynamicConductance

  !-----------------------------------------------------------------------
  subroutine LookupPsihatINI
    !
    ! !DESCRIPTION:
    ! Initialize the look-up tables needed to calculate the RSL psihat functions.
    ! Remember that in a netcdf file the dimensions appear in the opposite order:
    ! netcdf: psigridM_nc(nL,nZ) -> Fortran: psigridM(nZ,nL)
    ! netcdf: psigridH_nc(nL,nZ) -> Fortran: psigridH(nZ,nL)
    !
    ! Notes: zdt = (z-h)/(h-d) and dtL = (h-d)/L
    !
    ! !USES:
    use fileutils, only : getfil
    use ncdio_pio, only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t
    use ncdio_pio, only : ncd_inqdid, ncd_inqdlen
    use spmdMod, only : masterproc
    use MLclm_varcon, only : nZ, nL, dtLgridM, zdtgridM, psigridM, dtLgridH, zdtgridH, psigridH
    !
    ! !ARGUMENTS:
    implicit none
    !
    !LOCAL VARIABLES
    character(len=256) :: locfn        ! Local file name
    type(file_desc_t) :: ncid          ! pio netCDF file id
    integer :: dimid                   ! netCDF dimension id
    logical :: readv                   ! read variable in or not

    real(r8) :: zdtgridM_nc(nZ)        ! netcdf data: Grid of zdt on which psihat is given for momentum
    real(r8) :: dtLgridM_nc(nL)        ! netcdf data: Grid of dtL on which psihat is given for momentum
    real(r8) :: psigridM_nc(nL,nZ)     ! netcdf data: Grid of psihat values for momentum
    real(r8) :: zdtgridH_nc(nZ)        ! netcdf data: Grid of zdt on which psihat is given for heat
    real(r8) :: dtLgridH_nc(nL)        ! netcdf data: Grid of dtL on which psihat is given for heat
    real(r8) :: psigridH_nc(nL,nZ)     ! netcdf data: Grid of psihat values for heat
    integer :: nZ_nc, nL_nc            ! netcdf data: dimensions
    integer :: ii, jj                  ! Looping indices
    !---------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read RSL look-up table .....'
    end if

    ! Get netcdf file

    call getfil (rslfile, locfn, 0)

    ! Open netcdf file

    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Check dimensions

    call ncd_inqdid (ncid, 'nZ', dimid)
    call ncd_inqdlen (ncid, dimid, nZ_nc)

    if (nZ_nc /= nZ) then
       call endrun (msg=' ERROR: LookupPsihatINI: nZ does not equal expected value')
    end if

    call ncd_inqdid (ncid, 'nL', dimid)
    call ncd_inqdlen (ncid, dimid, nL_nc)

    if (nL_nc /= nL) then
       call endrun (msg=' ERROR: LookupPsihatINI: nL does not equal expected value')
    end if

    ! Read variables

    call ncd_io('dtLgridM', dtLgridM_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading dtLgridM')

    call ncd_io('zdtgridM', zdtgridM_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading zdtgridM')

    call ncd_io('psigridM', psigridM_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading psigridM')

    call ncd_io('dtLgridH', dtLgridH_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading dtLgridH')

    call ncd_io('zdtgridH', zdtgridH_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading zdtgridH')

    call ncd_io('psigridH', psigridH_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading psigridH')

    ! Close netcdf file

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*) 'Successfully read RSL look-up table'
    end if

    ! Copy netcdf variables

    do jj = 1, nL
       dtLgridM(1,jj) = dtLgridM_nc(jj)
       dtLgridH(1,jj) = dtLgridH_nc(jj)
    end do

    do ii = 1, nZ
       zdtgridM(ii,1) = zdtgridM_nc(ii)
       zdtgridH(ii,1) = zdtgridH_nc(ii)
    end do

    do ii = 1, nZ
       do jj = 1, nL
          psigridM(ii,jj) = psigridM_nc(jj,ii)
          psigridH(ii,jj) = psigridH_nc(jj,ii)
       end do
    end do

    return

  end subroutine LookupPsihatINI

end module MLCanopyTurbulenceMod
