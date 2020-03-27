!> Implements a constant increment nudging approach to data assimilation
module da_const_inc

  use iso_c_binding,      only : c_ptr
  use client_fortran_api, only : init_ssc_client, put_array, get_array
  use mpp_mod,            only : mpp_npes, mpp_get_current_pelist
  use MOM_coms,           only : PE_here
  use MOM_diag_remap,     only : diag_remap_ctrl
  use MOM_error_handler,  only : is_root_pe
  use MOM_EOS,            only : EOS_type
  use MOM_grid,           only : ocean_grid_type
  use MOM_file_parser,    only : get_param, param_file_type
  use MOM_remapping,      only : initialize_remapping
  use MOM_regridding,     only : set_regrid_params, get_regrid_size
  use MOM_time_manager,   only : time_type, get_date
  use MOM_unit_scaling,   only : unit_scale_type, unit_scaling_init
  use MOM_verticalGrid,   only : verticalGrid_type

  implicit none; private
#include <MOM_memory.h>

  type, public :: da_const_inc_type ; private

    type(c_ptr) :: ssc_data_client !< Pointer to data client
    real :: dt_da        !< Length of time between prior and posterior update
    real :: da_accum     !< Accumulated time since the last update
    real :: dt_inc       !< Length of time over which to apply an increment
    real :: inc_accum    !< Accumulated time since the start of increment applciation
    integer :: nk_obs       !< Number of levels in the observations
    logical :: do_temp      !< If true, send/receive increment for temperature
    logical :: do_salt      !< If true, send/receive increment for salinity
    logical :: use_ssc      !< If true, use the data client to get the observations
    real, dimension(:,:,:), allocatable :: temp_prior !< Temperature from the model remapped onto the obs grid
    real, dimension(:,:,:), allocatable :: salt_prior !< Salinity from the model remapped onto the obs grid
    real, dimension(:,:,:), allocatable :: temp_inc   !< Temperature increment from the data assimiliation (on obs grid)
    real, dimension(:,:,:), allocatable :: salt_inc   !< Salinity increment from the data assimiliation (on obs grid)
    real, dimension(:,:,:), allocatable :: obs_mask   !< 3d mask for the obs grid
    character(len=7)  :: id !< Identifier used as a key prefix
    character(len=50) :: sent_prior_key = "sent-prior" !< Signals whether the model has sent the prior
    character(len=50) :: sent_inc_key   = "sent-inc"   !< Signals whether python client has sent the increment
    character(len=50) :: temp_prior_key = "temp-prior" !< Name of key to set the prior for temperature
    character(len=50) :: temp_inc_key   = "temp-inc"   !< Name of key to receive the increment for temperature
    character(len=50) :: salt_prior_key = "salt-prior" !< Name of key to set the prior for salt
    character(len=50) :: salt_inc_key   = "salt-inc"   !< Name of key to receive the increment for salt

    type(diag_remap_ctrl) :: diag_remap !< Control structure to hold diagnostic grids
    type(EOS_type), pointer    :: eqn_of_state !< A pointer to the equation of state
  end type da_const_inc_type

  public da_const_inc_init, send_prior_recv_increment, apply_increment

contains

  subroutine da_const_inc_init(CS, G, GV, US, Time, param_file, h, T, S, eqn_of_state)
    type(da_const_inc_type), intent(inout) :: CS !< Control structure for data assimilation
    type(ocean_grid_type),   intent(in   ) :: G  !< Ocean grid structure
    type(verticalGrid_type), intent(in   ) :: GV !< ocean vertical grid structure
    type(unit_scale_type),   intent(in   ) :: US !< A dimensional unit scaling type
    type(time_type),         intent(in   ) :: Time !< Time at initialization
    type(param_file_type),   intent(in   ) :: param_file  !< Structure indicating parameter file to parse
    real, dimension(NIMEM_, NJMEM_, NKMEM_), intent(in   ) :: h !< Prognostic thicknesses
    real, dimension(NIMEM_, NJMEM_, NKMEM_), intent(in   ) :: T !< Prognostic temperature
    real, dimension(NIMEM_, NJMEM_, NKMEM_), intent(in   ) :: S !< Progonostic salinity
    type(EOS_type), target                 :: eqn_of_state

    character(len=40) :: mdl = trim("MOM_DA_CONST_INC")
    character(len=240), allocatable :: diag_coords(:)
    integer, dimension(:), allocatable :: pe_list
    integer, dimension(6) :: array_meta
    integer :: ni, nj, nk

    call get_param(param_file, mdl, "DT_DA", CS%dt_da,       &
        "Length of time (s) between DA updates", default=21600.)
    call get_param(param_file, mdl, "DA_TEMP", CS%do_temp,           &
        "Data assimilate temperature", default=.true.)
    call get_param(param_file, mdl, "DA_SALT", CS%do_salt,           &
        "Data assimilate salinity", default=.true.)
    call get_param(param_file, mdl, 'DA_COORDS', diag_coords, &
                 'A list of string tuples associating a coordinate definition for '//&
                 'data assimilation diagnostics. Each string \n'//&
                 'is of the form "MODULE_SUFFIX PARAMETER_SUFFIX COORDINATE_NAME".', &
                 default='z Z ZSTAR')
    call get_param(param_file, mdl, "DA_USE_SSC", CS%use_ssc, &
                  "If true, get the data for assimilation from a data client", default=.false.)
    if (CS%use_ssc) CS%ssc_data_client = init_ssc_client()
    CS%eqn_of_state => eqn_of_state
    ! Initialize the regridding
    call diag_remap_init(CS%diag_remap, diag_coords)
    call initialize_regridding(CS%diag_remap%regrid_cs, GV, US, GV%max_depth, param_file, mdl, &
             trim(CS%diag_remap%vertical_coord_name), "DA_COORDS", trim(CS%diag_remap%diag_coord_name))
    call set_regrid_params(CS%diag_remap%regrid_cs, min_thickness=0., integrate_downward_for_e=.false.)
    call initialize_remapping(CS%diag_remap%remap_cs, 'PPM_IH4', boundary_extrapolation=.false. )
    CS%diag_remap%nz = get_regrid_size(CS%diag_remap%regrid_cs)
    CS%nk_obs = CS%diag_remap%nz
    allocate(CS%obs_mask(SZI_(G),SZJ_(G),CS%nk_obs))

    call diag_remap_calc_hmask(CS%diag_remap, G, CS%obs_mask)
    call diag_remap_update(CS%diag_remap, G, GV, US, h, T, S, CS%eqn_of_state)

    ! Allocate all arrays used in data assimilation
    if (CS%do_temp) then
      allocate( CS%temp_prior(SZI_(G),SZJ_(G),CS%nk_obs) )
      CS%temp_prior(:,:,:) = 0.
      CS%temp_inc(:,:,:) = 0.
    endif
    if (CS%do_salt) then
      allocate( CS%salt_prior(SZI_(G),SZJ_(G),CS%nk_obs) )
      allocate( CS%salt_inc(SZI_(G),SZJ_(G),G%ke) )
      CS%salt_prior(:,:,:) = 0.
      CS%salt_inc(:,:,:) = 0.
    endif

    write(CS%id,"(I6.6,A)") PE_here(), "_"
    CS%sent_inc_key   = trim(CS%id)//CS%sent_inc_key
    CS%sent_prior_key = trim(CS%id)//CS%sent_prior_key
    CS%temp_inc_key   = trim(CS%id)//CS%temp_inc_key
    CS%temp_prior_key = trim(CS%id)//CS%temp_prior_key
    CS%salt_inc_key   = trim(CS%id)//CS%salt_inc_key
    CS%salt_prior_key = trim(CS%id)//CS%salt_prior_key

    ! Send all necessary metadata
    ni = G%ied-G%isd+1; nj = G%jed-G%jsd+1; nk = G%ke
    array_meta = [G%isd_global, G%jsd_global, 1, ni, nj, nk]
    call send_array(CS%ssc_data_client, trim(CS%id)//"array-meta", array_meta) 
    if (is_root_pe()) then
      allocate(pe_list(mpp_npes()))
      call mpp_get_current_pelist(pe_list)
      call send_array(CS%ssc_data_client, "rank-ids", pe_list)
      deallocate(pe_list)
      call send_time_to_client(CS,Time,"initial-time")
    endif
   end subroutine da_const_inc_init

   subroutine apply_increment(CS, G, GV, US, dt, h, T, S)
     type(da_const_inc_type), intent(inout) :: CS !< Control structure for data assimilation
     type(ocean_grid_type),   intent(in   ) :: G  !< Ocean grid structure
     type(verticalGrid_type), intent(in   ) :: GV !< ocean vertical grid structure
     type(unit_scale_type),   intent(in   ) :: US !< A dimensional unit scaling type
     real,                    intent(in   ) :: dt !< Thermodynamic timestep
     real, dimension(:,:,:),  intent(inout) :: h  !< The model's prognostic layer temperature
     real, dimension(:,:,:),  intent(inout) :: T  !< The model's prognostic layer temperature
     real, dimension(:,:,:),  intent(inout) :: S  !< The model's prognostic layer salinity

     real :: wt
     integer :: isc, iec, jsc, jec, i, j
     real, dimension(:,:,:), pointer :: h_obs
     real, dimension(SZK_(G)) :: inc_model

     isc = G%isc; iec = G%iec
     jsc = G%jsc; jec = G%jec
     wt = dt/CS%dt_da

     CS%inc_accum = mod(CS%inc_accum + dt,CS%dt_da)
     if (CS%inc_accum < CS%dt_inc) then
       call diag_remap_update(CS%diag_remap, G, GV, US, h, T, S, CS%eqn_of_state)
       if (CS%do_temp) then
         ! Remap increments onto the current model grid
         do i=isc,iec ; do j=jsc,jec
           call remapping_core_h(CS%diag_remap%remap_cs, CS%nk_obs, CS%diag_remap%h, CS%temp_inc, G%ke, h, inc_model)
           T(i,j,:) = T(i,j,:)+wt*inc_model(:)
         enddo; enddo
       endif
       if (CS%do_salt) then
         ! Remap increments onto the current model grid
         do i=isc,iec ; do j=jsc,jec
           call remapping_core_h(CS%diag_remap%remap_cs, CS%nk_obs, CS%diag_remap%h, CS%salt_inc, G%ke, h, inc_model)
           S(i,j,:) = S(i,j,:)+wt*inc_model(:)
         enddo; enddo
       endif
     endif

   end subroutine apply_increment

   !> Send the model state (the prior) and prepare to receive the increments for the posterior integration
   subroutine send_prior_recv_increment(CS, Time, G, GV, h, T, S)
     type(da_const_inc_type), intent(inout) :: CS !< Control structure for data assimilation
     type(time_type),         intent(in   ) :: Time !< Timestamp of the prior to be sent
     type(ocean_grid_type),   intent(in   ) :: G  !< Ocean grid structure
     type(verticalGrid_type), intent(in   ) :: GV !< ocean vertical grid structure
     real, dimension(:,:,:),  intent(in   ) :: h  !< The model's prognostic layer thicknesses
     real, dimension(:,:,:),  intent(in   ) :: T  !< The model's prognostic layer temperature
     real, dimension(:,:,:),  intent(in   ) :: S  !< The model's prognostic layer salinity

     integer :: isc, iec, jsc, jec, i, j, k
     real :: Idt

     Idt = 1./CS%dt_da
     isc = G%isc; iec = G%iec
     jsc = G%jsc; jec = G%jec

     if (is_root_pe()) then
       call send_time_to_client(CS,Time,"simulation-time")
       call set_scalar(CS%ssc_data_client, "sent-time", 1)
     endif
     ! Remap the model onto the observational grid 
     if (CS%do_temp) call diag_remap_do_remap(CS%diag_remap, G, GV, h, .false., .false., CS%obs_mask, 0., T, CS%temp_prior)
     if (CS%do_salt) call diag_remap_do_remap(CS%diag_remap, G, GV, h, .false., .false., CS%obs_mask, 0., S, CS%salt_prior)

     ! Check to make sure that the python client is ready to receive the prior (i.e. that it sent the increment from
     ! the last timestamp
     call poll_key_and_check_value(CS%ssc_data_client, CS%sent_inc_key, 0, 100, -1)
     if (CS%do_temp) call put_array(CS%ssc_data_client, CS%temp_prior_key, CS%temp_prior)
     if (CS%do_salt) call put_array(CS%ssc_data_client, CS%salt_prior_key, CS%salt_prior)
     ! Signal that prior has been sent and wait until the increment was sent
     call put_scalar(CS%ssc_data_client, CS%sent_prior_key, 1)
     call poll_key_and_check_value(CS%ssc_data_client, CS%sent_inc_key, 1, 100, -1)
     if (CS%do_temp) call get_array(CS%ssc_data_client, CS%temp_inc_key, CS%temp_inc)
     if (CS%do_salt) call get_array(CS%ssc_data_client, CS%salt_inc_key, CS%salt_inc)
     ! Reset sent prior_key for next loop 
     call put_scalar(CS%ssc_data_client, CS%sent_prior_key, 0)

     call pass_var(temp_inc,G%domain)
     call pass_var(salt_inc,G%domain)
     if (is_root_pe()) call set_scalar(CS%ssc_data_client, "sent-time", 0)

   end subroutine send_prior_recv_increment

   subroutine send_time_to_client( CS, Time, key ) 
     type(da_const_inc_type), intent(in) :: CS   !< Control structure for data assimilation
     type(time_type),         intent(in) :: Time !< Time to send
     character(len=*),        intent(in) :: key  !< Key used for the database

     integer :: yr, mon, day, hr, min, sec
     integer*8, dimension(6) :: date_vector
     call get_date(Time, yr, mon, day, hr, min, sec)
     date_vector = [ yr, mon, day, hr, min, sec ]
     call put_array(CS%ssc_data_client,key,date_vector)

   end subroutine send_time_to_client

end module da_const_inc
