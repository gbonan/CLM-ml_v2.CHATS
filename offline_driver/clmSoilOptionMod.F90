module clmSoilOptionMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! The offline code uncoupled to CLM reads a file with soil moisture on
  ! CLM's vertical grid. Snow/soil layers differ for CLM4.5 and CLM5. clm_phys
  ! is only needed for offline simulations uncoupled to CLM and is used for
  ! backward compatibility with older simulations that used CLM4.5 soils.
  ! Option to also adjust the soil moisture values. Values are set in namelist.

  implicit none

  character(len=6) :: clm_phys  ! CLM snow/soil layers. Options: 'CLM4_5' or 'CLM5_0'
  integer :: nlev_soil_adjust   ! Number of soil layers to apply soil moisture adjustment (use 0 to turn off)

end module clmSoilOptionMod
