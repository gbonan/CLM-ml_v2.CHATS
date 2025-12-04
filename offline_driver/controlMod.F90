module controlMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initialize namelist run control variables
  !
  ! !USES:
  use abortutils,   only : endrun
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: control
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine control (ntim, clm_start_ymd, clm_start_tod, fin_tower, fin_clm, fin_soil_adjust, dirout)
    !
    ! !DESCRIPTION:
    ! Initialize run control variables from namelist
    !
    ! !USES:
    use clm_time_manager, only : start_date_ymd, start_date_tod, dtstep
    use clm_varctl, only : iulog
    use clmSoilOptionMod, only : clm_phys, nlev_soil_adjust
    use TowerDataMod, only : ntower, tower_id, tower_num, tower_time
    use MLclm_varctl, only : met_type, dpai_min, pftcon_val
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(out) :: ntim                       ! Number of time steps to process
    integer, intent(out) :: clm_start_ymd              ! CLM file start date (yyyymmdd format)
    integer, intent(out) :: clm_start_tod              ! CLM file start time-of-day (seconds past 0Z UTC)
    character(len=256), intent(out) :: fin_tower       ! Tower meteorology file name
    character(len=256), intent(out) :: fin_clm         ! CLM file name
    character(len=256), intent(out) :: fin_soil_adjust ! Soil moisture adjustment factor file name
    character(len=256), intent(out) :: dirout          ! Model output file directory path
    !
    ! !LOCAL VARIABLES:
    character(len=6) :: tower_name                ! Flux tower site to process
    character(len=6) :: stop_option               ! Character flag to specify run length
    integer :: start_ymd                          ! Run start date in yyyymmdd format
    integer :: start_tod                          ! Time-of-day (UTC) of the start date (seconds past 0Z; 0 to 86400)
    integer :: stop_n                             ! Length of simulation

    integer :: i                                  ! Index
    integer :: steps_per_day                      ! Number of time steps per day

    namelist /clmML_inparm/ tower_name, start_ymd, start_tod, stop_option, &
    stop_n, fin_tower, fin_clm, clm_start_ymd, clm_start_tod, clm_phys, &
    fin_soil_adjust, nlev_soil_adjust, dirout, met_type, dpai_min, pftcon_val
    !---------------------------------------------------------------------

    ! Default namelist variables

    tower_name = ' '      ! Flux tower site to process
    start_ymd = 0         ! Run start date in yyyymmdd format
    start_tod = 0         ! Time-of-day (UTC) of the start date (seconds past 0Z; 0 to 86400)
    stop_option = ' '     ! Sets the run length as days ('ndays') or timesteps ('nsteps')
    stop_n = 0            ! Sets the length of the run (days or timesteps depending on stop_option)
    clm_start_ymd = 0     ! CLM file start date (yyyymmdd format)
    clm_start_tod = 0     ! CLM file start time-of-day (seconds past 0Z UTC)
    fin_tower = ' '       ! Tower meteorology file name
    fin_clm = ' '         ! CLM file name
    clm_phys = ' '        ! CLM snow/soil layers. Options: 'CLM4_5' or 'CLM5_0'
    fin_soil_adjust = ' ' ! Soil moisture adjustment factor file name
    nlev_soil_adjust = 0  ! Number of soil layers to apply soil moisture adjustment (use 0 to turn off)
    dirout = ' '          ! Model output file directory path

    ! These three variables are set to default values in MLclm_varctl.F90. Use the namelist to set
    ! them to tower-specific values.
    !
    ! met_type = 0          ! Meteorological forcing for multilayer canopy timestep:
                            ! 0 = no interpolation (uses standard CLM calendar)
                            ! 2 = 2-point interpolation (not supported)
                            ! 3 = 3-point interpolation (time is centered in timestep as in CHATS)
    ! dpai_min = 0.01D0     ! Minimum plant area index (normalized) to be considered a vegetation layer (m2/m2)
    ! pftcon_val = 0        ! PFT parameters: use default values (0) or override for CHATS (1)

    ! Read namelist file

    write(iulog,*) 'Attempting to read namelist file .....'
    read (5, clmML_inparm)
    write(iulog,*) 'Successfully read namelist file'

    ! Set calendar variables

    start_date_ymd = start_ymd
    start_date_tod = start_tod

    ! Match tower site to correct index for TowerDataMod arrays

    tower_num = 0
    do i = 1, ntower
       if (tower_name == tower_id(i)) then
          tower_num = i
          exit
       else
          cycle
       end if
    end do

    if (tower_num == 0) then
       write (iulog,*) 'control error: tower site = ',tower_name, ' not found'
       call endrun()
    end if

    ! Time step of forcing data (in seconds). This varies among tower sites.

    dtstep = tower_time(tower_num) * 60

    ! Set length of simulation

    if (stop_option == 'nsteps') then
       ntim = stop_n                       ! Number of time steps to execute
    else if (stop_option == 'ndays') then
       steps_per_day = 86400 / dtstep      ! Number of time steps per day
       ntim = steps_per_day * stop_n       ! Number of time steps to execute
    end if

  end subroutine control

end module controlMod
