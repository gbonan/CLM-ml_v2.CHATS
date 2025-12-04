module clm_driver

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Main CLM model driver to calculate fluxes
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  use ColumnType, only : col
  use decompMod, only : bounds_type
  use clm_instMod
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_drv
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_drv (bounds, time_indx, fin1, fin2)
    !
    ! !DESCRIPTION:
    ! Main CLM model driver to calculate fluxes
    !
    ! !USES:
    use clm_varpar, only : nlevgrnd, nlevsno
    use clmDataMod, only : clmData
    use filterMod, only : filter, setExposedvegpFilter
    use SurfaceAlbedoMod, only : SoilAlbedo
    use SurfaceResistanceMod, only : calc_soilevap_resis
    use MLSoilTemperatureMod, only : SoilTemperature, SoilThermProp
    use SoilWaterMovementMod, only : SoilWater
    use MLCanopyFluxesMod, only : MLCanopyFluxes  !!! CLMml !!!
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds     ! CLM bounds
    integer, intent(in) :: time_indx            ! Time index from reference date (0Z January 1 of current year, when calday = 1.000)
    character(len=256) :: fin1, fin2            ! File name
    !
    ! !LOCAL VARIABLES:
    integer  :: f                                ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                ! Column index for CLM g/l/c/p hierarchy

    real(r8) :: cv (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! CLM: soil heat capacity (J/m2/K)
    real(r8) :: tk (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! CLM: soil thermal conductivity at layer interface (W/m/K)
    real(r8) :: tk_h2osfc(bounds%begc:bounds%endc)                 ! CLM: thermal conductivity of h2osfc (W/m/K)
    !---------------------------------------------------------------------

    associate ( &
    snl            => col%snl                                  , &  ! Number of snow layers
    frac_veg_nosno => canopystate_inst%frac_veg_nosno_patch    , &  ! Fraction of vegetation not covered by snow (0 or 1)
    frac_sno_eff   => waterdiagnosticbulk_inst%frac_sno_eff_col, &  ! Effective fraction of ground covered by snow (0 to 1)
    h2osno         => water_inst%h2osno_col                    , &  ! Total snow water (kg H2O/m2)
    h2osfc         => waterstatebulk_inst%h2osfc_col             &  ! Surface water (kg H2O/m2)
    )

    ! Read CLM data for current time slice

    call clmData (fin1, fin2, time_indx, bounds%begp, bounds%endp, bounds%begc, bounds%endc, &
    soilstate_inst, waterstatebulk_inst, canopystate_inst, surfalb_inst)

    ! Set CLM frac_veg_nosno and its filter (filter%exposedvegp)

    do p = bounds%begp, bounds%endp
       frac_veg_nosno(p) = 1
    end do
    call setExposedvegpFilter (filter, frac_veg_nosno)

    ! Calculate CLM soil albedo

    call SoilAlbedo (bounds, filter%num_nourbanc, filter%nourbanc, waterstatebulk_inst, surfalb_inst)

    ! Calculate CLM moisture stress/resistance for soil evaporation

    call calc_soilevap_resis (bounds, filter%num_nolakec, filter%nolakec, &
    soilstate_inst, waterstatebulk_inst, temperature_inst)

    ! Zero out snow and surface water

    do f = 1, filter%num_nolakec
       c = filter%nolakec(f)
       snl(c) = 0
       frac_sno_eff(c) = 0._r8
       h2osno(c) = 0._r8
       h2osfc(c) = 0._r8
    end do

    ! Calculate CLM thermal conductivity and heat capacity. Only need
    ! thk(c,snl(c)+1), which is the thermal conductivity of the first
    ! snow/soil layer.

    call SoilThermProp (bounds, filter%num_nolakec, filter%nolakec, tk(bounds%begc:bounds%endc,:), &
    cv(bounds%begc:bounds%endc,:), tk_h2osfc(bounds%begc:bounds%endc), temperature_inst, &
    waterdiagnosticbulk_inst, waterstatebulk_inst, water_inst, soilstate_inst)

    ! CLM hydraulic conductivity and soil matric potential

    call SoilWater (bounds, filter%num_hydrologyc, filter%hydrologyc, &
    soilstate_inst, waterstatebulk_inst)

    ! CLMml: Multilayer canopy and soil fluxes

    call MLCanopyFluxes (bounds, filter%num_exposedvegp, filter%exposedvegp, &
    atm2lnd_inst, canopystate_inst, soilstate_inst, temperature_inst, &
    waterstatebulk_inst, waterfluxbulk_inst, &
    energyflux_inst, frictionvel_inst, surfalb_inst, solarabs_inst, &
    mlcanopy_inst, wateratm2lndbulk_inst, waterdiagnosticbulk_inst)

    ! Update CLM soil temperatures

    !!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This soil temperature routine is specific to the multilayer !
    ! canopy and is not the same is in CLM                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call SoilTemperature (bounds, filter%num_nolakec, filter%nolakec, &
    soilstate_inst, temperature_inst, waterdiagnosticbulk_inst, &
    waterstatebulk_inst, water_inst, mlcanopy_inst)

    end associate
  end subroutine clm_drv

end module clm_driver
