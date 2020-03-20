!> Implements a constant increment nudging approach to data assimilation
module da_const_inc

  use MOM_diag_remap,   only : diag_remap_ctrl,
  use MOM_grids,        only : ocean_grid_type
  use MOM_time_manager, only : time_type
  use MOM_verticalGrid, only : verticalGrid_type

  implicit none; private
#include <MOM_memory.h>

  type, public :: da_const_inc_type ; private
    integer :: dt_da        !< Length of time between prior and posterior update
    integer :: dt_inc       !< Length of time over which to apply an increment
    integer :: inc_accum    !< Accumulated time since the start of increment applciation
    integer :: nk_obs       !< Number of levels in the observations
    logical :: do_temp      !< If true, send/receive increment for temperature
    logical :: do_salt      !< If true, send/receive increment for salinity
    real, dimension(:,:,:), allocatable :: temp_prior !< Temperature from the model remapped onto the obs grid
    real, dimension(:,:,:), allocatable :: salt_prior !< Salinity from the model remapped onto the obs grid
    real, dimension(:,:,:), allocatable :: temp_obs   !< Temperature from the observations
    real, dimension(:,:,:), allocatable :: salt_obs   !< Salinity from the observations
    real, dimension(:,:,:), allocatable :: temp_inc   !< Novelty used based on the data assimiliation
    real, dimension(:,:,:), allocatable :: salt_inc   !< Novelty used based on the data assimiliation
    real, dimension(:,:,:), allocatable :: obs_mask   !< 3d mask for the obs grid

    type(diag_remap_ctrl) :: remap_CS !< Control structure to hold diagnostic grids
  end type da_const_inc_type

  public da_const_inc_init, send_prior, receive_increment

contains

  subroutine da_const_inc_init(CS, GV, US, param_file, T, S)
    type(da_const_inc_type), intent(in   ) :: CS !< Control structure for data assimilation
    type(verticalGrid_type), intent(in   ) :: GV !< ocean vertical grid structure
    type(unit_scale_type),   intent(in   ) :: US !< A dimensional unit scaling type
    type(param_file_type),   intent(in   ) :: param_file  !< Structure indicating parameter file to parse
    real, dimension(NIMEM_, NJMEM_, NKMEM_), target, intent(in   ) :: T !< Prognostic temperature
    real, dimension(NIMEM_, NJMEM_, NKMEM_), target, intent(in   ) :: S !< Progonostic salinity

    character(len=40) :: mdl = "MOM_DA_CONST_INC"
    character(len=240), allocatable :: diag_coords(:)

    call get_param(param_file, mdl, "DT_DA", CS%dt_da,       &
        "Length of time (s) between DA updates", default=21600)
    call get_param(param_file, mdl, "DA_TEMP", CS%do_temp,           &
        "Data assimilate temperature", default=.true.)
    call get_param(param_file, mdl, "DA_SALT", CS%do_salt,           &
        "Data assimilate salinity", default=.true.)
    call get_param(param_file, mdl, 'DA_COORDS', diag_coords, &
                 'A list of string tuples associating a coordinate definition for '//&
                 'data assimilation diagnostics. Each string \n'//&
                 'is of the form "MODULE_SUFFIX PARAMETER_SUFFIX COORDINATE_NAME".', &
                 default='z Z ZSTAR')

    ! Initialize the regridding
    call diag_remap_init(CS%remap_CS, diag_coords)
    call initialize_regridding(CS%remap_cs%regrid_cs, GV, US, GV%max_depth, param_file, mdl, &
             trim(CS%remap_cs%vertical_coord_name), "DA_COORDS", trim(CS%remap_cs%diag_coord_name))
    call set_regrid_params(CS%remap_cs%regrid_cs, min_thickness=0., integrate_downward_for_e=.false.)
    call initialize_remapping(remap_cs%remap_cs, 'PPM_IH4', boundary_extrapolation=.false., &
                               answers_2018=remap_cs%answers_2018)
    CS%remap_cs%nz = get_regrid_size(CS%remap_cs%regrid_cs)
    CS%nk_obs = CS%remap_cs%nz
    allocate(CS%obs_mask(SZI_(G),SZJ_(G),CS%nk_obs))

    call diag_remap_calc_hmask(CS%remap_cs, G, CS%obs_mask)
    call diag_remap_update(CS%remap_cs, G, GV, US, h, T, S, equn_of_state)

    ! Allocate all arrays used in data assimilation
    if (CS%do_temp) then
      allocate( CS%temp_prior(SZI_(G),SZJ_(G),CS%nk_obs) )
      allocate( CS%temp_obs(SZI_(G),SZJ_(G),CS%nk_obs) )
      allocate( CS%temp_inc(SZI_(G),SZJ_(G),G%ke) )
      CS%temp_prior(:,:,:) = 0.
      CS%temp_obs(:,:,:) = 0.
      CS%temp_inc(:,:,:) = 0.
    endif
    if (CS%do_salt) then
      allocate( CS%salt_prior(SZI_(G),SZJ_(G),CS%nk_obs) )
      allocate( CS%salt_obs(SZI_(G),SZJ_(G),CS%nk_obs) )
      allocate( CS%salt_inc(SZI_(G),SZJ_(G),G%ke) )
      CS%salt_prior(:,:,:) = 0.
      CS%salt_obs(:,:,:) = 0.
      CS%salt_inc(:,:,:) = 0.
    endif

   end subroutine da_const_inc_init

   subroutine calculate_increment(CS, G, h, T, S)
     type(da_const_inc_type), intent(in   ) :: CS !< Control structure for data assimilation
     type(ocean_grid_type),   intent(in   ) :: G  !< Ocean grid structure
     real, dimension(:,:,:),  intent(in)    :: h  !< The model's prognostic layer thicknesses
     real, dimension(:,:,:),  intent(in)    :: T  !< The model's prognostic layer temperature
     real, dimension(:,:,:),  intent(in)    :: S  !< The model's prognostic layer salinity

     integer :: isc, iec, jsc, jec
     real :: Idt
     isc = G%isc; iec = G%iec
     jsc = G%jsc; jec = G%jec

     Idt = 1./CS%dt_da

     if (CS%do_temp) then
       call get_array(CS%ssc_data_client, "obs_temp", CS%temp_obs)
       ! Remap to the model grid
       do i=isc,iec ; do j=jsc,jec; do k=1,G%ke
         CS%temp_inc(i,j,k) = Idt*(CS%temp_inc(i,j,k) - T(i,j,k))
       enddo; enddo; enddo;
     endif
     if (CS%do_salt) then
       call get_array(CS%ssc_data_client, "obs_salt", CS%salt_obs)
       ! Remap to the model grid
       do i=isc,iec ; do j=jsc,jec; do k=1,G%ke
         CS%salt_inc(i,j,k) = Idt*(CS%salt_inc(i,j,k) - S(i,j,k))
       enddo; enddo; enddo;
     endif
   end subroutine calculate_increment

   subroutine apply_increment(CS, G, dt, T, S)
     type(da_const_inc_type), intent(in   ) :: CS !< Control structure for data assimilation
     type(ocean_grid_type),   intent(in   ) :: G  !< Ocean grid structure
     real,                    intent(in   ) :: dt !< Thermodynamic timestep
     real, dimension(:,:,:),  intent(inout) :: T  !< The model's prognostic layer temperature
     real, dimension(:,:,:),  intent(inout) :: S  !< The model's prognostic layer salinity
     
     integer :: isc, iec, jsc, jec
     isc = G%isc; iec = G%iec
     jsc = G%jsc; jec = G%jec

     CS%inc_accum = mod(CS%inc_accum + dt,CS%dt_da)
     if (CS%inc_accum < CS%dt_inc) t hen
       if (CS%do_temp) then
         do i=isc,iec; do j=jsc,jec; do k=1,G%ke
           T(i,j,k) = T(i,j,k) + CS%temp_inc*dt
         enddo; enddo; enddo
       endif
       if (CS%do_salt) then
         do i=isc,iec; do j=jsc,jec; do k=1,G%ke
           S(i,j,k) = S(i,j,k) + CS%salt_inc*dt
         enddo; enddo; enddo
       endif
     endif

   end subroutine apply_increment

   subroutine send_prior(CS, G, h, T, S)
     type(da_const_inc_type), intent(in   ) :: CS !< Control structure for data assimilation
     type(ocean_grid_type),   intent(in   ) :: G  !< Ocean grid structure
     type(verticalGrid_type), intent(in   ) :: GV !< ocean vertical grid structure
     real, dimension(:,:,:),  intent(in)    :: h  !< The model's prognostic layer thicknesses
     real, dimension(:,:,:),  intent(in)    :: T  !< The model's prognostic layer temperature
     real, dimension(:,:,:),  intent(in)    :: S  !< The model's prognostic layer salinity

     if (CS%do_temp) then
       call diag_remap_do_remap(CS%remap_cs, G, GV, h, .false., .false., CS%obs_mask, 0., T, CS%temp_prior)
       call put_array(CS%ssc_data_client, "prior_temp", CS%salt_temp)
     endif
     if (CS%do_salt) then
       call diag_remap_do_remap(CS%remap_cs, G, GV, h, .false., .false., CS%obs_mask, 0., S, CS%salt_prior)
       call put_array(CS%ssc_data_client, "prior_salt", CS%salt_prior)
     endif

   end subroutine send_prior

end module da_constant_increment
