module MLGetAtmForcingMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Atmospheric forcing for current mutilayer canopy timestep
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use abortutils        , only : endrun
  use clm_varctl        , only : iulog
  use MLCanopyFluxesType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: GetAtmForcing          ! Atmospheric forcing for current mutilayer canopy timestep
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: TimeInterpolation2    ! 2-point linear interpolation of atmospheric forcing at CLM timestep to ML timestep
  private :: TimeInterpolation3    ! 3-point linear interpolation of atmospheric forcing at CLM timestep to ML timestep
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  function TimeInterpolation2 (x0, x1, t0, t1, tx) result (x)
    !
    ! !DESCRIPTION:
    ! 2-point linear interpolation of atmospheric forcing at CLM timestep
    ! to ML timestep. Forcing is linear from x0 to x1 with x1 attained at t1.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: t0           ! Calendar day for preceeding CLM timestep
    real(r8), intent(in)  :: t1           ! Calendar day for current CLM timestep
    real(r8), intent(in)  :: tx           ! Calendar day for ML interpolation
    real(r8), intent(in)  :: x0           ! Value for preceeding CLM timestep (at time t0)
    real(r8), intent(in)  :: x1           ! Value for current CLM timestep (at time t1)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: b0                        ! Intercept term for linear interpolation
    real(r8) :: b1                        ! Slope term for linear interpolation
    real(r8) :: x                         ! Interpolated value at time tx
    !---------------------------------------------------------------------

    b1 = (x1 - x0) / (t1 - t0)
    b0 = x1 - b1 * t1
    x = b0 + b1 * tx

  end function TimeInterpolation2

  !-----------------------------------------------------------------------
  function TimeInterpolation3 (x0, x1, x2, t0, t1, t2, tx) result (x)
    !
    ! !DESCRIPTION:
    ! 3-point linear interpolation of atmospheric forcing at CLM timestep to ML timestep
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: t0           ! Calendar day for preceeding CLM timestep
    real(r8), intent(in)  :: t1           ! Calendar day for current CLM timestep
    real(r8), intent(in)  :: t2           ! Calendar day for next CLM timestep
    real(r8), intent(in)  :: tx           ! Calendar day for ML interpolation
    real(r8), intent(in)  :: x0           ! Value for preceeding CLM timestep (at time t0)
    real(r8), intent(in)  :: x1           ! Value for current CLM timestep (at time t1)
    real(r8), intent(in)  :: x2           ! Value for next CLM timestep (at time t2)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: b0                        ! Intercept term for linear interpolation
    real(r8) :: b1                        ! Slope term for linear interpolation
    real(r8) :: x                         ! Interpolated value at time tx
    !---------------------------------------------------------------------

    if (tx < t1) then
       b1 = (x1 - x0) / (t1 - t0)
       b0 = x0 - b1 * t0
    else if (tx > t1) then
       b1 = (x2 - x1) / (t2 - t1)
       b0 = x1 - b1 * t1
    end if
    x = b0 + b1 * tx

  end function TimeInterpolation3

  !-----------------------------------------------------------------------
  subroutine GetAtmForcing (time_bef, time_cur, time_next, time_ml, num_filter, &
  filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Atmospheric forcing for current mutilayer canopy timestep
    !
    ! !USES:
    use clm_varpar, only : ivis, inir
    use MLclm_varctl, only : met_type
    use MLclm_varcon, only : mmh2o, mmdry, cpd, cpw, rgas, lapse_rate, wind_forc_min
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: time_bef      ! Calendar day for preceeding CLM timestep
    real(r8), intent(in) :: time_cur      ! Calendar day for current CLM timestep
    real(r8), intent(in) :: time_next     ! Calendar day for next CLM timestep
    real(r8), intent(in) :: time_ml       ! Calendar day for ML interpolation
    integer, intent(in)  :: num_filter    ! Number of patches in filter
    integer, intent(in)  :: filter(:)     ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                        ! Filter index
    integer  :: p                         ! Patch index for CLM g/l/c/p hierarchy
    !---------------------------------------------------------------------

    associate ( &
                                                               ! *** Input ***
    tref_bef       => mlcanopy_inst%tref_bef_forcing      , &  ! Air temperature at reference height (K) [previous CLM timestep]
    tref_cur       => mlcanopy_inst%tref_cur_forcing      , &  ! Air temperature at reference height (K) [current CLM timestep]
    tref_next      => mlcanopy_inst%tref_next_forcing     , &  ! Air temperature at reference height (K) [next CLM timestep]
    qref_bef       => mlcanopy_inst%qref_bef_forcing      , &  ! Specific humidity at reference height (kg/kg) [previous CLM timestep]
    qref_cur       => mlcanopy_inst%qref_cur_forcing      , &  ! Specific humidity at reference height (kg/kg) [current CLM timestep]
    qref_next      => mlcanopy_inst%qref_next_forcing     , &  ! Specific humidity at reference height (kg/kg) [next CLM timestep]
    uref_bef       => mlcanopy_inst%uref_bef_forcing      , &  ! Wind speed at reference height (m/s) [previous CLM timestep]
    uref_cur       => mlcanopy_inst%uref_cur_forcing      , &  ! Wind speed at reference height (m/s) [current CLM timestep]
    uref_next      => mlcanopy_inst%uref_next_forcing     , &  ! Wind speed at reference height (m/s) [next CLM timestep]
    pref_bef       => mlcanopy_inst%pref_bef_forcing      , &  ! Air pressure at reference height (Pa) [previous CLM timestep]
    pref_cur       => mlcanopy_inst%pref_cur_forcing      , &  ! Air pressure at reference height (Pa) [current CLM timestep]
    pref_next      => mlcanopy_inst%pref_next_forcing     , &  ! Air pressure at reference height (Pa) [next CLM timestep]
    co2ref_bef     => mlcanopy_inst%co2ref_bef_forcing    , &  ! Atmospheric CO2 at reference height (umol/mol) [previous CLM timestep]
    co2ref_cur     => mlcanopy_inst%co2ref_cur_forcing    , &  ! Atmospheric CO2 at reference height (umol/mol) [current CLM timestep]
    co2ref_next    => mlcanopy_inst%co2ref_next_forcing   , &  ! Atmospheric CO2 at reference height (umol/mol) [next CLM timestep]
    swskyb_bef     => mlcanopy_inst%swskyb_bef_forcing    , &  ! Atmospheric direct beam solar radiation (W/m2) [previous CLM timestep]
    swskyb_cur     => mlcanopy_inst%swskyb_cur_forcing    , &  ! Atmospheric direct beam solar radiation (W/m2) [current CLM timestep]
    swskyb_next    => mlcanopy_inst%swskyb_next_forcing   , &  ! Atmospheric direct beam solar radiation (W/m2) [next CLM timestep]
    swskyd_bef     => mlcanopy_inst%swskyd_bef_forcing    , &  ! Atmospheric diffuse solar radiation (W/m2) [previous CLM timestep]
    swskyd_cur     => mlcanopy_inst%swskyd_cur_forcing    , &  ! Atmospheric diffuse solar radiation (W/m2) [current CLM timestep]
    swskyd_next    => mlcanopy_inst%swskyd_next_forcing   , &  ! Atmospheric diffuse solar radiation (W/m2) [next CLM timestep]
    lwsky_bef      => mlcanopy_inst%lwsky_bef_forcing     , &  ! Atmospheric longwave radiation (W/m2) [previous CLM timestep]
    lwsky_cur      => mlcanopy_inst%lwsky_cur_forcing     , &  ! Atmospheric longwave radiation (W/m2) [current CLM timestep]
    lwsky_next     => mlcanopy_inst%lwsky_next_forcing    , &  ! Atmospheric longwave radiation (W/m2) [next CLM timestep]
    zref           => mlcanopy_inst%zref_forcing          , &  ! Atmospheric reference height (m)
                                                               ! *** Output ***
    tref           => mlcanopy_inst%tref_forcing          , &  ! Interpolated: Air temperature at reference height (K)
    qref           => mlcanopy_inst%qref_forcing          , &  ! Interpolated: Specific humidity at reference height (kg/kg)
    uref           => mlcanopy_inst%uref_forcing          , &  ! Interpolated: Wind speed at reference height (m/s)
    pref           => mlcanopy_inst%pref_forcing          , &  ! Interpolated: Air pressure at reference height (Pa)
    co2ref         => mlcanopy_inst%co2ref_forcing        , &  ! Interpolated: Atmospheric CO2 at reference height (umol/mol)
    swskyb         => mlcanopy_inst%swskyb_forcing        , &  ! Interpolated: Atmospheric direct beam solar radiation (W/m2)
    swskyd         => mlcanopy_inst%swskyd_forcing        , &  ! Interpolated: Atmospheric diffuse solar radiation (W/m2)
    lwsky          => mlcanopy_inst%lwsky_forcing         , &  ! Interpolated: Atmospheric longwave radiation (W/m2)
    thref          => mlcanopy_inst%thref_forcing         , &  ! Derived: Atmospheric potential temperature at reference height (K)
    thvref         => mlcanopy_inst%thvref_forcing        , &  ! Derived: Atmospheric virtual potential temperature at reference height (K)
    eref           => mlcanopy_inst%eref_forcing          , &  ! Derived: Vapor pressure at reference height (Pa)
    rhoair         => mlcanopy_inst%rhoair_forcing        , &  ! Derived: Air density at reference height (kg/m3)
    rhomol         => mlcanopy_inst%rhomol_forcing        , &  ! Derived: Molar density at reference height (mol/m3)
    mmair          => mlcanopy_inst%mmair_forcing         , &  ! Derived: Molecular mass of air at reference height (kg/mol)
    cpair          => mlcanopy_inst%cpair_forcing           &  ! Derived: Specific heat of air (constant pressure) at reference height (J/mol/K)
    )

    do fp = 1, num_filter
       p = filter(fp)

       ! Atmospheric forcing

       select case (met_type)
       case (0)
          ! No interpolation
          uref(p) = uref_cur(p)
          tref(p) = tref_cur(p)
          qref(p) = qref_cur(p)
          pref(p) = pref_cur(p)
          co2ref(p) = co2ref_cur(p)
          swskyb(p,ivis) = swskyb_cur(p,ivis)
          swskyd(p,ivis) = swskyd_cur(p,ivis)
          swskyb(p,inir) = swskyb_cur(p,inir)
          swskyd(p,inir) = swskyd_cur(p,inir)
          lwsky(p) = lwsky_cur(p)
       case (2)
          ! Interpolation: 2-point linear interpolation from "bef" to "cur"
          call endrun (msg=' ERROR: met_type not valid')
          uref(p) = TimeInterpolation2 (uref_bef(p), uref_cur(p), time_bef, time_cur, time_ml)
          tref(p) = TimeInterpolation2 (tref_bef(p), tref_cur(p), time_bef, time_cur, time_ml)
          qref(p) = TimeInterpolation2 (qref_bef(p), qref_cur(p), time_bef, time_cur, time_ml)
          pref(p) = TimeInterpolation2 (pref_bef(p), pref_cur(p), time_bef, time_cur, time_ml)
          co2ref(p) = TimeInterpolation2 (co2ref_bef(p), co2ref_cur(p), time_bef, time_cur, time_ml)
          swskyb(p,ivis) = TimeInterpolation2 (swskyb_bef(p,ivis), swskyb_cur(p,ivis), time_bef, time_cur, time_ml)
          swskyd(p,ivis) = TimeInterpolation2 (swskyd_bef(p,ivis), swskyd_cur(p,ivis), time_bef, time_cur, time_ml)
          swskyb(p,inir) = TimeInterpolation2 (swskyb_bef(p,inir), swskyb_cur(p,inir), time_bef, time_cur, time_ml)
          swskyd(p,inir) = TimeInterpolation2 (swskyd_bef(p,inir), swskyd_cur(p,inir), time_bef, time_cur, time_ml)
          lwsky(p) = TimeInterpolation2 (lwsky_bef(p), lwsky_cur(p), time_bef, time_cur, time_ml)
       case (3)
          ! Interpolation: 3-point linear interpolation using "bef", "cur", and "next"
          uref(p) = TimeInterpolation3 (uref_bef(p), uref_cur(p), uref_next(p), time_bef, time_cur, time_next, time_ml)
          tref(p) = TimeInterpolation3 (tref_bef(p), tref_cur(p), tref_next(p), time_bef, time_cur, time_next, time_ml)
          qref(p) = TimeInterpolation3 (qref_bef(p), qref_cur(p), qref_next(p), time_bef, time_cur, time_next, time_ml)
          pref(p) = TimeInterpolation3 (pref_bef(p), pref_cur(p), pref_next(p), time_bef, time_cur, time_next, time_ml)
          co2ref(p) = TimeInterpolation3 (co2ref_bef(p), co2ref_cur(p), co2ref_next(p), time_bef, time_cur, time_next, time_ml)
          swskyb(p,ivis) = TimeInterpolation3 (swskyb_bef(p,ivis), swskyb_cur(p,ivis), swskyb_next(p,ivis), time_bef, time_cur, time_next, time_ml)
          swskyd(p,ivis) = TimeInterpolation3 (swskyd_bef(p,ivis), swskyd_cur(p,ivis), swskyd_next(p,ivis), time_bef, time_cur, time_next, time_ml)
          swskyb(p,inir) = TimeInterpolation3 (swskyb_bef(p,inir), swskyb_cur(p,inir), swskyb_next(p,inir), time_bef, time_cur, time_next, time_ml)
          swskyd(p,inir) = TimeInterpolation3 (swskyd_bef(p,inir), swskyd_cur(p,inir), swskyd_next(p,inir), time_bef, time_cur, time_next, time_ml)
          lwsky(p) = TimeInterpolation3 (lwsky_bef(p), lwsky_cur(p), lwsky_next(p), time_bef, time_cur, time_next, time_ml)
       end select

       ! Minimum wind speed restriction

       uref(p) = max (wind_forc_min, uref(p))

       ! Derived atmospheric variables

       eref(p) = qref(p) * pref(p) / (mmh2o / mmdry + (1._r8 - mmh2o / mmdry) * qref(p))
       rhomol(p) = pref(p) / (rgas * tref(p))
       rhoair(p) = rhomol(p) * mmdry * (1._r8 - (1._r8 - mmh2o/mmdry) * eref(p) / pref(p))
       mmair(p) = rhoair(p) / rhomol(p)
       cpair(p) = cpd * (1._r8 + (cpw/cpd - 1._r8) * qref(p)) * mmair(p)
       thref(p) = tref(p) + lapse_rate * zref(p)
       thvref(p) = thref(p) * (1._r8 + 0.61_r8 * qref(p))

    end do

    end associate
  end subroutine GetAtmForcing

end module MLGetAtmForcingMod
