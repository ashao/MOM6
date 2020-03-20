!> Implements a constant increment nudging approach to data assimilation
module da_const_inc

  use iso_c_binding,      only : c_ptr
  use client_fortran_api, only : init_ssc_client, put_array, get_array
  use MOM_diag_remap,     only : diag_remap_ctrl
  use MOM_grids,          only : ocean_grid_type
  use MOM_file_parser,    only : get_param, param_file_type
  use MOM_time_manager,   only : time_type
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

    type(diag_remap_ctrl) :: diag_remap !< Control structure to hold diagnostic grids
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
    call get_param(param_file, mdl, "DA_USE_SSC", CS%use_ssc, &
                  "If true, get the data for assimilation from a data client", default=.false.)
    if (CS%usc_ssc) CS%ssc_data_client = init_ssc_client
    ! Initialize the regridding
    call diag_remap_init(CS%diag_remap, diag_coords)
    call initialize_regridding(CS%diag_remap%remap_cs%regrid_cs, GV, US, GV%max_depth, param_file, mdl, &
             trim(CS%diag_remap%remap_cs%vertical_coord_name), "DA_COORDS", trim(CS%diag_remap%remap_cs%diag_coord_name))
    call set_regrid_params(CS%diag_remap%remap_cs%regrid_cs, min_thickness=0., integrate_downward_for_e=.false.)
    call initialize_remapping(CS%diag_remap%remap_cs, 'PPM_IH4', boundary_extrapolation=.false. )
    CS%diag_remap%remap_cs%nz = get_regrid_size(CS%diag_remap%remap_cs%regrid_cs)
    CS%nk_obs = CS%diag_remap%remap_cs%nz
    allocate(CS%obs_mask(SZI_(G),SZJ_(G),CS%nk_obs))

    call diag_remap_calc_hmask(CS%diag_remap, G, CS%obs_mask)
    call diag_remap_update(CS%diag_remap, G, GV, US, h, T, S, equn_of_state)

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

   end subroutine da_const_inc_init

   subroutine apply_increment(CS, G, dt, h, T, S)
     type(da_const_inc_type), intent(in   ) :: CS !< Control structure for data assimilation
     type(ocean_grid_type),   intent(in   ) :: G  !< Ocean grid structure
     real,                    intent(in   ) :: dt !< Thermodynamic timestep
     real, dimension(:,:,:),  intent(inout) :: h  !< The model's prognostic layer temperature
     real, dimension(:,:,:),  intent(inout) :: T  !< The model's prognostic layer temperature
     real, dimension(:,:,:),  intent(inout) :: S  !< The model's prognostic layer salinity

     integer :: isc, iec, jsc, jec
     real, dimension(:), pointer :: h_obs
     rela, dimension(SZK_(G)) :: inc_model

     h_obs => CS%diag_remap%h

     isc = G%isc; iec = G%iec
     jsc = G%jsc; jec = G%jec

     CS%inc_accum = mod(CS%inc_accum + dt,CS%dt_da)

     if (CS%inc_accum < CS%dt_inc) then
       if (CS%do_temp) then
         ! Remap increments onto the current model grid
         do i=isc,iec ; do j=jsc,jec
           call remapping_core_h(CS%diag_remap%remap_cs, CS%nk_obs, h_obs, CS%temp_inc, G%ke, h, inc_model)
           T(i,j,:) = T(i,j,:)+dt*inc_model(:)
         enddo; enddo
       endif
       if (CS%do_salt) then
         ! Remap increments onto the current model grid
         do i=isc,iec ; do j=jsc,jec
           call remapping_core_h(CS%diag_remap%remap_cs, CS%nk_obs, h_obs, CS%salt_inc, G%ke, h, inc_model)
           S(i,j,:) = S(i,j,:)+dt*inc_model(:)
         enddo; enddo
       endif
     endif

   end subroutine apply_increment

   !> Send the model state (the prior) and prepare to receive the 'posterior' increments
   subroutine send_prior_recv_increment(CS, G, GV, dt, h, T, S)
     type(da_const_inc_type), intent(in) :: CS !< Control structure for data assimilation
     type(ocean_grid_type),   intent(in) :: G  !< Ocean grid structure
     type(verticalGrid_type), intent(in) :: GV !< ocean vertical grid structure
     real                                :: dt !< Timestep over which this process is applied
     real, dimension(:,:,:),  intent(in) :: h  !< The model's prognostic layer thicknesses
     real, dimension(:,:,:),  intent(in) :: T  !< The model's prognostic layer temperature
     real, dimension(:,:,:),  intent(in) :: S  !< The model's prognostic layer salinity

     integer :: isc, iec, jsc, jec
     real :: Idt

     Idt = 1./dt
     isc = G%isc; iec = G%iec
     jsc = G%jsc; jec = G%jec


     CS%da_accum = MOD(CS%da_accum+dt,CS%dt_da)
     if (CS%da_accum == 0.) then
       if (CS%do_temp) then
         ! Set status key for send prior to not done
         ! Set status key for receive increment to not ready
         call diag_remap_do_remap(CS%remap_cs, G, GV, h, .false., .false., CS%obs_mask, 0., T, CS%temp_prior)
         call put_array(CS%ssc_data_client, "temp_prior", CS%temp_prior)
         ! Set status key of send to done
         call get_array(CS%ssc_data_client, "temp_inc", CS%temp_inc)
         do i=isc,iec; do j=jsc,jec; do k=1,G%ke
           CS%temp_inc(i,j,k) = CS%temp_inc(i,j,k)*Idt
         enddo; enddo; enddo
       endif
       if (CS%do_salt) then
         ! Set status key for send prior to not done
         ! Set status key for receive increment to not ready
         call diag_remap_do_remap(CS%remap_cs, G, GV, h, .false., .false., CS%obs_mask, 0., S, CS%salt_prior)
         call put_array(CS%ssc_data_client, "salt_prior", CS%salt_prior)
         ! Set status key of send to done
         call get_array(CS%ssc_data_client, "salt_inc", CS%salt_inc)
         do i=isc,iec; do j=jsc,jec; do k=1,G%ke
           CS%salt_inc(i,j,k) = CS%salt_inc(i,j,k)*Idt
         enddo; enddo; enddo
       endif
     endif

   end subroutine update_prior_posterior

end module da_constant_increment
