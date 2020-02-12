!> Contains the routines necessary to stream data into and out of the model to
!! a compatible orchestrator
module MOM_data_client
  use MOM_diag_mediator,        only : diag_ctrl, register_diag_field, post_data
  use MOM_error_handler,        only : MOM_error, FATAL
  use MOM_time_manager,         only : time_type, time_type_to_real
  use MOM_grid,                 only : ocean_grid_type
  use MOM_verticalGrid,         only : verticalGrid_type, get_thickness_units
  use MOM_unit_scaling,         only : unit_scale_type
  
  ! This file is part of MOM6. See LICENSE.md for the license.
  
  implicit none ; private
#include <MOM_memory.h>
  
  !> Contains pointers to arrays in the model that might be potentially streamed
  !! into or out of the model
  type, public :: data_client_type ; private
  
    real, dimension(:,:,:), pointer :: h   !< Layer thicknesses
    real, dimension(:,:,:), pointer :: T   !< Temperature
    real, dimension(:,:,:), pointer :: S   !< Salinity
    real, dimension(:,:,:), pointer :: uh  !< Mass (thickness) transport u-flux
    real, dimension(:,:,:), pointer :: vh  !< Mass (thickness) transport v-flux
    real, dimension(:,:)  , pointer :: SSH !< SSH
  
    integer :: id_streamin_h    = -1 !< ID for streaming in h 
    integer :: id_streamin_T    = -1 !< ID for streaming in T
    integer :: id_streamin_S    = -1 !< ID for streaming in S
    integer :: id_streamin_uh   = -1 !< ID for streaming in uh
    integer :: id_streamin_vh   = -1 !< ID for streaming in vh
    integer :: id_streamin_SSH  = -1 !< ID for streaming in SSH
    integer :: id_streamout_h   = -1 !< ID for streaming out h    
    integer :: id_streamout_T   = -1 !< ID for streaming out T    
    integer :: id_streamout_S   = -1 !< ID for streaming out S    
    integer :: id_streamout_uh  = -1 !< ID for streaming out uh   
    integer :: id_streamout_vh  = -1 !< ID for streaming out vh   
    integer :: id_streamout_SSH = -1 !< ID for streaming out SSH  
  
  end type data_client_type

  public :: MOM_data_client_init, streamout_data, streamin_data

  contains

  subroutine MOM_data_client_init(CS, GV, diag, Time, US, h, T, S, uh, vh, SSH)
    type(data_client_type), pointer,  intent(inout) :: CS    !< data client control structure
    type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
    type(diag_ctrl),         intent(inout) :: diag  !< regulates diagnostic output
    type(time_type),         intent(in)    :: Time  !< current model time
    type(unit_scale_type),   pointer       :: US    !< Pointer to a structure containing
                                                    ! various unit conversion factors
    real, dimension(NIMEM_,NJMEM_,NKMEM_), target :: h !< Layer thicknesses 
    real, dimension(NIMEM_,NJMEM_,NKMEM_), target :: T !< Layer temperatures
    real, dimension(NIMEM_,NJMEM_,NKMEM_), target :: S !< Layer saliniites
    real, dimension(NIMEM_,NJMEM_,NKMEM_), target :: uh !< Zonal mass transport (x) 
    real, dimension(NIMEM_,NJMEM_,NKMEM_), target :: vh !< Meridional mass transport
    real, dimension(NIMEM_,NJMEM_),        target :: SSH !< Sea surface height

    real :: H_convert
    character(len=48) :: thickness_units

    thickness_units = get_thickness_units(GV)
    if (GV%Boussinesq) then
      H_convert = GV%H_to_m
    else
      H_convert = GV%H_to_kg_m2
    endif


    if (associated(CS)) then
      call MOM_error(FATAL,"Data client CS is already initialized")
    endif

    allocate(CS)

    CS%id_streamin_h = register_diag_field('ocean_model', 'streamin_h', diag%axesTL, Time, &
      'Layer thicknesses passed into MOM6 from the orchestrator', 'm', &
      v_extensive=.true., conversion=H_convert)
    CS%id_streamin_T = register_diag_field('ocean_model', 'streamin_T', diag%axesTL, Time, &
        'Temperature passed into MOM6 from the orchestrator', 'C')
    CS%id_streamin_S = register_diag_field('ocean_model', 'streamin_S', diag%axesTL, Time, &
        'Salinity passed into MOM6 from the orchestrator', 'psu')
    CS%id_streamin_uh = register_diag_field('ocean_model', 'streamin_uh', diag%axesCuL, Time, &
        'Ocean Mass X Transport passed into MOM6', 'kg s-1', v_extensive=.true.)
    CS%id_streamin_vh = register_diag_field('ocean_model', 'streamin_vh', diag%axesCvL, Time, &
        'Ocean Mass Y Transport passed into MOM6', 'kg s-1', v_extensive=.true.)
    CS%id_streamin_ssh = register_diag_field('ocean_model', 'streamin_SSH', diag%axesT1, &
        Time, 'Instantaneous Sea Surface Height passed into MOM6 from the orchestrator', 'm')
    CS%id_streamout_h = register_diag_field('ocean_model', 'streamout_h', diag%axesTL, Time, &
      'Layer thicknesses passed out of MOM6 to the orchestrator', 'm', &
      v_extensive=.true., conversion=H_convert)
    CS%id_streamout_T = register_diag_field('ocean_model', 'streamout_T', diag%axesTL, Time, &
        'Temperature passed out of MOM6 to the orchestrator', 'C')
    CS%id_streamout_S = register_diag_field('ocean_model', 'streamout_S', diag%axesTL, Time, &
        'Salinity passed out of MOM6 to the orchestrator', 'psu')
    CS%id_streamout_uh = register_diag_field('ocean_model', 'streamout_uh', diag%axesCuL, Time, &
        'Ocean Mass X Transport passed out of MOM6', 'kg s-1', v_extensive=.true.)
    CS%id_streamout_vh = register_diag_field('ocean_model', 'streamout_vh', diag%axesCvL, Time, &
        'Ocean Mass Y Transport passed out of MOM6', 'kg s-1', v_extensive=.true.)
    CS%id_streamout_ssh = register_diag_field('ocean_model', 'streamout_SSH', diag%axesT1, &
        Time, 'Instantaneous Sea Surface Height passed out of MOM6 to the orchestrator', 'm')

    if ( (CS%id_streamin_h > 0) .or. (CS%id_streamout_h > 0) ) then
      CS%h => h 
    endif
    if ( (CS%id_streamin_T > 0) .or. (CS%id_streamout_T > 0) ) then
      CS%T => T 
    endif
    if ( (CS%id_streamin_S > 0) .or. (CS%id_streamout_S > 0) ) then
      CS%S => S 
    endif
    if ( (CS%id_streamin_uh > 0) .or. (CS%id_streamout_uh > 0) ) then
      CS%uh => uh 
    endif
    if ( (CS%id_streamin_vh > 0) .or. (CS%id_streamout_vh > 0) ) then
      CS%vh => vh 
    endif
    if ( (CS%id_streamin_ssh > 0) .or. (CS%id_streamout_ssh > 0) ) then
      CS%ssh => ssh 
    endif

  end subroutine MOM_data_client_init

  subroutine streamout_data(CS, Time, G, GV)
    type(data_client_type),    intent(inout) :: CS    !< data client control structure
    type(time_type),         intent(in)    :: Time  !< current model time
    type(ocean_grid_type),   intent(in)    :: G     !< ocean grid structure
    type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure

    integer :: isc, iec, jsc, jec, nk
    real :: time_real

    isc = G%isc; iec = G%iec
    jsc = G%jsc; jec = G%jec
    nk = GV%ke

    time_real = time_type_to_real(Time)
    if (CS%id_streamout_h > 0)   call send_array_3d(CS%h,   CS%id_streamout_h,   time_real, isc, iec, jsc, jec, nk) 
    if (CS%id_streamout_T > 0)   call send_array_3d(CS%T,   CS%id_streamout_T,   time_real, isc, iec, jsc, jec, nk)
    if (CS%id_streamout_S > 0)   call send_array_3d(CS%S,   CS%id_streamout_S,   time_real, isc, iec, jsc, jec, nk)
    !if (CS%id_streamout_uh > 0)  call send_array_3d(CS%uh,  CS%id_streamout_uh,  time_real, isc, iec, jsc, jec, nk)
    !if (CS%id_streamout_vh > 0)  call send_array_3d(CS%vh,  CS%id_streamout_vh,  time_real, isc, iec, jsc, jec, nk)
    if (CS%id_streamout_ssh > 0) call put_2d_array_double("ssh",CS%ssh, isc, iec, jsc, jec, .true.)

  end subroutine streamout_data
  
  subroutine streamin_data(CS, Time, G, GV)
    type(data_client_type),    intent(inout) :: CS    !< data client control structure
    type(time_type),         intent(in)    :: Time  !< current model time
    type(ocean_grid_type),   intent(in)    :: G     !< ocean grid structure
    type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure

   integer :: isc, iec, jsc, jec, nk
    real :: time_real

    isc = G%isc; iec = G%iec
    jsc = G%jsc; jec = G%jec
    nk = GV%ke

    time_real = time_type_to_real(Time)
    !if (CS%id_streamin_h)   call recv_array_3d(CS%h,   CS%id_streamin_h,   time_real, isc, iec, jsc, jec, nk) 
    !if (CS%id_streamin_T)   call recv_array_3d(CS%T,   CS%id_streamin_T,   time_real, isc, iec, jsc, jec, nk)
    !if (CS%id_streamin_S)   call recv_array_3d(CS%S,   CS%id_streamin_S,   time_real, isc, iec, jsc, jec, nk)
    !if (CS%id_streamin_uh)  call recv_array_3d(CS%uh,  CS%id_streamin_uh,  time_real, isc, iec, jsc, jec, nk)
    !if (CS%id_streamin_vh)  call recv_array_3d(CS%vh,  CS%id_streamin_vh,  time_real, isc, iec, jsc, jec, nk)
    !if (CS%id_streamin_ssh) call recv_array_2d(CS%ssh, CS%id_streamin_ssh, time_real, isc, iec, jsc, jec)

  end subroutine streamin_data

end module MOM_data_client
