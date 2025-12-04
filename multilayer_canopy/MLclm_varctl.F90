module MLclm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing multilayer canopy model run control variables
  !
  ! !USES:
  use clm_varcon, only : ispval
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none

  ! Stomatal conductance and photosynthesis

  integer  :: gs_type = 2                 ! Stomatal conductance: Medlyn (0), Ball-Berry (1), or WUE optimization (2)
  integer  :: gspot_type = 1              ! Stomatal conductance: use potential conductance (0) or water-stressed conductance (1)
  integer  :: gs_solver = 2               ! Stomatal conductance: numerical solver for WUE optimization uses (1) Brent or (2) bisection
  integer  :: colim_type = 1              ! Photosynthesis: minimum rate (0) or co-limited rate (1)
  integer  :: acclim_type = 1             ! Photosynthesis: temperature acclimation off (0) or on (1)
  real(r8) :: kn_val = -999._r8           ! Canopy nitrogen profile: for a user-specified Kn, set kn_val to desired value > 0

  ! Canopy turbulence

  integer  :: turb_type = 1               ! Turbulence parameterization: H&F roughness sublayer (1)
  integer  :: sparse_canopy_type = 1      ! H&F roughness sublayer: sparse canopy off (0) or on (1)
  integer  :: HF_extension_type = 2       ! Aerodynamic conductance at ground: extend H&F to ground (1) or use log profile (2)
  integer  :: flux_profile_type = 1       ! Flux-profile solution: dataset (-1) or well-mixed assumption (0) or implicit (1)
  integer  :: gb_type = 3                 ! Boundary layer conductance: laminar only (1) or also turbulent (2) and free convection (3)

  ! Radiative transfer

  integer  :: light_type = 2              ! Solar radiative transfer: Norman (1) or two-stream approximation (2)
  integer  :: leaf_optics_type = 0        ! Leaf optical properties: Constant with height (0) or vary with height (1)
  integer  :: longwave_type = 1           ! Longwave radiative transfer: Norman (1)

  ! Multilayer canopy timestepping (must be evenly divisible into CLM timestep: num_ml_steps = dtime_clm / dtime_ml)

  real(r8) :: dtime_ml = 300._r8          ! Multilayer canopy timestep (s)

  ! Runge-Kutta method

  integer, parameter :: runge_kutta_type = 41       ! Euler (10)
                                                    ! 2nd-order (21:trapezoidal, 22:midpoint; 23:Ralston)
                                                    ! 3rd-order (31:Huen, 32:Ralston, 33:Kutta)
                                                    ! 4th-order (41:Kutta)
  integer, parameter :: nrk = (runge_kutta_type/10) ! Number of Runge-Kutta steps

  ! Use dz_tall or dz_short to determine the number of canopy layers ...

  real(r8) :: dz_tall = 0.5_r8            ! Height increment for tall canopies > dz_param (m)
  real(r8) :: dz_short = 0.1_r8           ! Height increment for short canopies <= dz_param (m)
  real(r8) :: dz_param = 2._r8            ! Height above which a canopy is tall (m)

  ! ... or directly specify the number of layers (used if these are > 0)

  integer  :: nlayer_above = 0            ! The number of above-canopy layers
  integer  :: nlayer_within = 0           ! The number of within-canopy layers

  ! Miscellaneous

  integer  :: mlcan_to_clm = 0            ! Pass multilayer canopy fluxes to CLM for use by CAM: no (0), yes (1)
  integer  :: ml_vert_init = ispval       ! Flag used to initialize multilayer canopy vertical structure and profiles

  ! The next variables can also be set in the clmML_inparm namelist (when uncoupled to CLM).
  ! These are the default settings. Tower-specific settings are found in the tower namelist file.

  integer  :: met_type = 0                ! Meteorological forcing for multilayer canopy timestep: 
                                          ! 0 = no interpolation (uses standard CLM calendar)
                                          ! 2 = 2-point interpolation (not supported)
                                          ! 3 = 3-point interpolation (time is centered in timestep as in CHATS)
  real(r8) :: dpai_min = 0.01_r8          ! Minimum plant area index (normalized) to be considered a vegetation layer (m2/m2)
  integer  :: pftcon_val = 0              ! PFT parameters: use default values (0) or override for CHATS (1)

end module MLclm_varctl
