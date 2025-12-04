module clm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing run control variables
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none

  integer, save :: iulog = 6                                              ! "stdout" log file unit number
  character(len=256), public :: rslfile = '../rsl_lookup_tables/psihat.nc'  ! RSL psihat look-up tables

end module clm_varctl
