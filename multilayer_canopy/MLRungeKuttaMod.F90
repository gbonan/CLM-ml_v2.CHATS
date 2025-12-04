module MLRungeKuttaMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Runge-Kutta method for state updates
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
  public :: RungeKuttaIni      ! Initialize Butcher tableu for Runge-Kutta methods
  public :: RungeKuttaUpdate   ! Runge-Kutta state update
  !-----------------------------------------------------------------------

  contains

  !-----------------------------------------------------------------------
  subroutine RungeKuttaUpdate (irk, a, b, c, num_filter, filter, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Runge-Kutta state update: Update states for next Runge-Kutta step
    !
    ! !USES:
    use MLclm_varpar, only : isun, isha
    use MLclm_varctl, only : nrk
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: irk         ! Current Runge-Kutta step
    real(r8), intent(in) :: a(nrk,nrk)  ! Runge-Kutta parameters
    real(r8), intent(in) :: b(nrk)      ! Runge-Kutta parameters
    real(r8), intent(in) :: c(nrk)      ! Runge-Kutta parameters (not used)
    integer, intent(in)  :: num_filter  ! Number of patches in filter
    integer, intent(in)  :: filter(:)   ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                      ! Filter index
    integer  :: p                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                      ! Aboveground layer index
    integer  :: j                       ! Runge-Kutta step index
    !---------------------------------------------------------------------

    associate ( &
                                                         ! *** Input ***
    ncan        => mlcanopy_inst%ncan_canopy        , &  ! Number of aboveground layers
    tair_bef    => mlcanopy_inst%tair_bef_profile   , &  ! Canopy layer air temperature for previous timestep (K)
    dtair       => mlcanopy_inst%dtair_profile      , &  ! Change in canopy layer air temperature over Runge-Kutta step (K)
    eair_bef    => mlcanopy_inst%eair_bef_profile   , &  ! Canopy layer vapor pressure for previous timestep (Pa)
    deair       => mlcanopy_inst%deair_profile      , &  ! Change in canopy layer vapor pressure over Runge-Kutta step (Pa)
    h2ocan_bef  => mlcanopy_inst%h2ocan_bef_profile , &  ! Canopy layer intercepted water for previous timestep (kg H2O/m2)
    dh2ocan     => mlcanopy_inst%dh2ocan_profile    , &  ! Change in canopy layer intercepted water over Runge-Kutta step (kg H2O/m2)
    tleaf_bef   => mlcanopy_inst%tleaf_bef_leaf     , &  ! Leaf temperature for previous timestep (K)
    dtleaf      => mlcanopy_inst%dtleaf_leaf        , &  ! Change in leaf temperature over Runge-Kutta step (K)
    lwp_bef     => mlcanopy_inst%lwp_bef_leaf       , &  ! Leaf water potential for previous timestep (MPa)
    dlwp        => mlcanopy_inst%dlwp_leaf          , &  ! Change in leaf water potential over Runge-Kutta step (MPa)
    tg_bef      => mlcanopy_inst%tg_bef_soil        , &  ! Soil surface temperature for previous timestep (K)
    dtg         => mlcanopy_inst%dtg_soil           , &  ! Change in soil surface temperature over Runge-Kutta step (K)
                                                         ! *** Input/Output ***
    tair        => mlcanopy_inst%tair_profile       , &  ! Canopy layer air temperature (K)
    eair        => mlcanopy_inst%eair_profile       , &  ! Canopy layer vapor pressure (Pa)
    h2ocan      => mlcanopy_inst%h2ocan_profile     , &  ! Canopy layer intercepted water (kg H2O/m2)
    tleaf       => mlcanopy_inst%tleaf_leaf         , &  ! Leaf temperature (K) 
    lwp         => mlcanopy_inst%lwp_leaf           , &  ! Leaf water potential (MPa)
    tg          => mlcanopy_inst%tg_soil              &  ! Soil surface temperature (K)
    )

    ! General algorithm for p-order Runge-Kutta to solve the differential
    ! equation dy/dt = f(t,y) with h the step size:
    !
    ! k(1) = h * f[t_n, y_n]                                                : slope (dy) at the beginning of the interval, using  y
    ! k(2) = h * f[t_n + c(2) h, y_n + a(21) k(1)]                          : slope (dy) at intermediate point in the interval, using y and k(1)
    ! k(3) = h * f[t_n + c(3) h, y_n + a(31) k(1) + a(32) k(2)]             : slope (dy) at intermediate point in the interval, using y and k(1) and k2)
    !                         ...
    ! k(p) = h * f[t_n + c(p) h, y_n + a(p,1) k(1) + ... + a(p,p-1) k(p-1)] : slope (dy) at intermediate point in the interval, using y and k(1), ..., k(p-1)
    !
    ! y_{n+1} = y_n + [b(1) k(1) + b(2) k(2) + b(3) k(3) + ... + b(p) k(p)] : y is updated using weighted estimate of slopes

    do fp = 1, num_filter
       p = filter(fp)

       ! tair, eair, tleaf, lwp, and h2ocan

       do ic = 1, ncan(p)

          ! Change in state variables at Runge-Kutta step irk

          dtair(p,ic,irk) = tair(p,ic) - tair_bef(p,ic)
          deair(p,ic,irk) = eair(p,ic) - eair_bef(p,ic)
          dh2ocan(p,ic,irk) = h2ocan(p,ic) - h2ocan_bef(p,ic)
          dtleaf(p,ic,isun,irk) = tleaf(p,ic,isun) - tleaf_bef(p,ic,isun)
          dtleaf(p,ic,isha,irk) = tleaf(p,ic,isha) - tleaf_bef(p,ic,isha)
          dlwp(p,ic,isun,irk) = lwp(p,ic,isun) - lwp_bef(p,ic,isun)
          dlwp(p,ic,isha,irk) = lwp(p,ic,isha) - lwp_bef(p,ic,isha)

          ! Update state variables

          if (irk < nrk) then

             ! Intermediate state update for Runge-Kutta step irk+1 using derivatives from steps 1, ..., irk

             tair(p,ic) = tair_bef(p,ic)
             eair(p,ic) = eair_bef(p,ic)
             h2ocan(p,ic) = h2ocan_bef(p,ic)
             tleaf(p,ic,isun) = tleaf_bef(p,ic,isun)
             tleaf(p,ic,isha) = tleaf_bef(p,ic,isha)
             lwp(p,ic,isun) = lwp_bef(p,ic,isun)
             lwp(p,ic,isha) = lwp_bef(p,ic,isha)

             do j = 1, irk
                tair(p,ic) = tair(p,ic) + a(irk+1,j) * dtair(p,ic,j)
                eair(p,ic) = eair(p,ic) + a(irk+1,j) * deair(p,ic,j)
                h2ocan(p,ic) = h2ocan(p,ic) + a(irk+1,j) * dh2ocan(p,ic,j)
                tleaf(p,ic,isun) = tleaf(p,ic,isun) + a(irk+1,j) * dtleaf(p,ic,isun,j)
                tleaf(p,ic,isha) = tleaf(p,ic,isha) + a(irk+1,j) * dtleaf(p,ic,isha,j)
                lwp(p,ic,isun) = lwp(p,ic,isun) + a(irk+1,j) * dlwp(p,ic,isun,j)
                lwp(p,ic,isha) = lwp(p,ic,isha) + a(irk+1,j) * dlwp(p,ic,isha,j)
             end do

          else if (irk == nrk) then

             ! Weighted state values used in final subroutine calls

             tair(p,ic) = tair_bef(p,ic)
             eair(p,ic) = eair_bef(p,ic)
             h2ocan(p,ic) = h2ocan_bef(p,ic)
             tleaf(p,ic,isun) = tleaf_bef(p,ic,isun)
             tleaf(p,ic,isha) = tleaf_bef(p,ic,isha)
             lwp(p,ic,isun) = lwp_bef(p,ic,isun)
             lwp(p,ic,isha) = lwp_bef(p,ic,isha)

             do j = 1, nrk
                tair(p,ic) = tair(p,ic) + b(j) * dtair(p,ic,j)
                eair(p,ic) = eair(p,ic) + b(j) * deair(p,ic,j)
                h2ocan(p,ic) = h2ocan(p,ic) + b(j) * dh2ocan(p,ic,j)
                tleaf(p,ic,isun) = tleaf(p,ic,isun) + b(j) * dtleaf(p,ic,isun,j)
                tleaf(p,ic,isha) = tleaf(p,ic,isha) + b(j) * dtleaf(p,ic,isha,j)
                lwp(p,ic,isun) = lwp(p,ic,isun) + b(j) * dlwp(p,ic,isun,j)
                lwp(p,ic,isha) = lwp(p,ic,isha) + b(j) * dlwp(p,ic,isha,j)
             end do

          end if

       end do

       ! Now for tg

       dtg(p,irk) = tg(p) - tg_bef(p)

       if (irk < nrk) then

          tg(p) = tg_bef(p)
          do j = 1, irk
             tg(p) = tg(p) + a(irk+1,j) * dtg(p,j)
          end do

       else if (irk == nrk) then

          tg(p) = tg_bef(p)
          do j = 1, nrk
             tg(p) = tg(p) + b(j) * dtg(p,j)
          end do

       end if

    end do

    end associate
  end subroutine RungeKuttaUpdate

  !-----------------------------------------------------------------------
  subroutine RungeKuttaIni (a, b, c)
    !
    ! !DESCRIPTION:
    ! Initialize Butcher tableu for Runge-Kutta methods
    !
    ! !USES:
    use clm_varcon, only : spval
    use MLclm_varctl, only : runge_kutta_type, nrk
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(out) :: a(nrk,nrk), b(nrk), c(nrk) ! Runge-Kutta parameters
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    ! runge_kutta_type defines the Runge-Kutta method

    ! References:
    ! Kutta, W. (1901). Beitrag zur näherungsweisen Integration totaler Differentialgleichungen.
    ! Z. Math. Phys. (Zeitschrift für Mathematik und Physik), 46, 435-453.
    !
    ! Butcher, J. C. (1996). A history of Runge-Kutta methods. Applied Numerical
    ! Mathematics, 20, 247-260.
    !
    ! Ralston, A. (1962). Runge-Kutta methods with minimum error bounds.
    ! Mathematics of Computation, 16, 431-437.

    ! General algorithm for p-order Runge-Kutta to solve the differential
    ! equation dy/dt = f(t,y) with h the step size:
    !
    ! k(1) = h * f[t_n, y_n]
    ! k(2) = h * f[t_n + c(2) h, y_n + a(21) k(1)]
    ! k(3) = h * f[t_n + c(3) h, y_n + a(31) k(1) + a(32) k(2)]
    !                         ...
    ! k(p) = h * f[t_n + c(p) h, y_n + a(p,1) k(1) + ... + a(p,p-1) k(p-1)]
    !
    ! y_{n+1} = y_n + [b(1) k(1) + b(2) k(2) + b(3) k(3) + ... + b(p) k(p)]

    select case (runge_kutta_type)
    case (21)
       ! 2nd-order (trapezoidal rule; Butcher 1996)
       a(1,1) = spval
       a(1,2) = spval
       a(2,1) = 1._r8
       a(2,2) = spval

       b(1) = 1._r8 / 2._r8
       b(2) = 1._r8 / 2._r8

       c(1) = 0._r8
       c(2) = 1._r8
    case (22)
       ! 2nd-order (midpoint rule; Butcher 1996)
       a(1,1) = spval
       a(1,2) = spval
       a(2,1) = 1._r8 / 2._r8
       a(2,2) = spval

       b(1) = 0._r8
       b(2) = 1._r8

       c(1) = 0._r8
       c(2) = 1._r8 / 2._r8
    case (23)
       ! 2nd-order (Ralston 1962)
       a(1,1) = spval
       a(1,2) = spval
       a(2,1) = 2._r8 / 3._r8
       a(2,2) = spval

       b(1) = 1._r8 / 4._r8
       b(2) = 3._r8 / 4._r8

       c(1) = 0._r8
       c(2) = 2._r8 / 3._r8
    case (31)
       ! 3rd-order (Huen's method; Butcher 1996)
       a(1,1) = spval
       a(1,2) = spval
       a(1,3) = spval
       a(2,1) = 1._r8 / 3._r8
       a(2,2) = spval
       a(2,3) = spval
       a(3,1) = 0._r8
       a(3,2) = 2._r8 / 3._r8
       a(3,3) = spval

       b(1) = 1._r8 / 4._r8
       b(2) = 0._r8
       b(3) = 3._r8 / 4._r8

       c(1) = 0._r8
       c(2) = 1._r8 / 3._r8
       c(3) = 2._r8 / 3._r8
    case (32)
       ! 3rd-order (Ralston 1962)
       a(1,1) = spval
       a(1,2) = spval
       a(1,3) = spval
       a(2,1) = 1._r8 / 2._r8
       a(2,2) = spval
       a(2,3) = spval
       a(3,1) = 0._r8
       a(3,2) = 3._r8 / 4._r8
       a(3,3) = spval

       b(1) = 2._r8 / 9._r8
       b(2) = 3._r8 / 9._r8
       b(3) = 4._r8 / 9._r8

       c(1) = 0._r8
       c(2) = 1._r8 / 2._r8
       c(3) = 3._r8 / 4._r8
    case (33)
       ! 3rd-order (Kutta's method; Kutta 1901)
       a(1,1) = spval
       a(1,2) = spval
       a(1,3) = spval
       a(2,1) = 1._r8 / 2._r8
       a(2,2) = spval
       a(2,3) = spval
       a(3,1) = -1._r8
       a(3,2) = 2._r8
       a(3,3) = spval

       b(1) = 1._r8 / 6._r8
       b(2) = 4._r8 / 6._r8
       b(3) = 1._r8 / 6._r8

       c(1) = 0._r8
       c(2) = 1._r8 / 2._r8
       c(3) = 1._r8
    case (41)
       ! 4th-order (Kutta's method; Butcher 1996 & Kutta 1901)
       a(1,1) = spval
       a(1,2) = spval
       a(1,3) = spval
       a(1,4) = spval
       a(2,1) = 1._r8 / 2._r8
       a(2,2) = spval
       a(2,3) = spval
       a(2,4) = spval
       a(3,1) = 0._r8
       a(3,2) = 1._r8 / 2._r8
       a(3,3) = spval
       a(3,4) = spval
       a(4,1) = 0._r8
       a(4,2) = 0._r8
       a(4,3) = 1._r8
       a(4,4) = spval

       b(1) = 1._r8 / 6._r8
       b(2) = 2._r8 / 6._r8
       b(3) = 2._r8 / 6._r8
       b(4) = 1._r8 / 6._r8

       c(1) = 0._r8
       c(2) = 1._r8 / 2._r8
       c(3) = 1._r8 / 2._r8
       c(4) = 1._r8
    end select

  end subroutine RungeKuttaIni

end module MLRungeKuttaMod
