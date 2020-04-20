!> Contains the routines necessary to stream data into and out of the model to
!! a compatible orchestrator
module MOM_data_client
  use MOM_domains,        only : pass_var
  use MOM_diag_mediator,        only : diag_ctrl, register_diag_field, post_data
  use MOM_error_handler,        only : MOM_error, FATAL
  use MOM_time_manager,         only : time_type, time_type_to_real
  use MOM_grid,                 only : ocean_grid_type
  use MOM_verticalGrid,         only : verticalGrid_type, get_thickness_units
  use MOM_unit_scaling,         only : unit_scale_type
  use iso_c_binding,            only : c_ptr
  use client_fortran_api,       only : init_ssc_client, put_array, get_array
  use MOM_coms,                 only : PE_here
  implicit none ; private

  ! This file is part of MOM6. See LICENSE.md for the license.

#include <MOM_memory.h>

  !> Contains pointers to arrays in the model that might be potentially streamed
  !! into or out of the model
  type, public :: data_client_type ; private
    type(c_ptr)                     :: ssc_data_client !< Pointer to data client
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

    integer :: accumulated_time

    real :: timing_current !< Store the timing for the current send
    real :: timing_total   !< Store the total time spent sending/receiving

    integer :: num_elements !< Total number of numbers in the current send
    integer :: timing_unit  !< Unit to write output of timing information

  end type data_client_type

  public :: MOM_data_client_init, streamout_data, streamin_data

  contains
  subroutine MOM_data_client_init(CS, G, GV, diag, Time, US, h, T, S, uh, vh, SSH)
    type(data_client_type), pointer,  intent(inout) :: CS    !< data client control structure
    type(ocean_grid_type), intent(in)      :: G     !< ocean horizontal vertical grid structure
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
    character(len=1024) :: key, fname

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
!      CS%h => h
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

    CS%ssc_data_client = init_ssc_client()
    CS%accumulated_time = 0
    CS%timing_total = 0.

    ! Send metadata needed to reconstruct global array
    write(key,("(I6.6,A)")) PE_here(),"_rank-meta"
    call put_array(CS%ssc_data_client, key, REAL([G%isd_global, G%idg_offset, G%jsd_global, G%jdg_offset]))

    ! Open file to be used for writing timing
    write(fname,("(A,I6.6)")) "ssc-timing_",PE_here()
 !   open(newunit=CS%timing_unit, file=trim(fname), status='REPLACE')
 !   write(CS%timing_unit,*) "MODEL_TIMESTAMP TOTAL_TIME_ELAPSED TIME_ELAPSED NUMBER_OF_ELEMENTS_SENT"
  end subroutine MOM_data_client_init

  subroutine streamout_data(CS, Time, G, GV, dt, h)
    type(data_client_type),    intent(inout) :: CS    !< data client control structure
    type(time_type),         intent(in)    :: Time  !< current model time
    type(ocean_grid_type),   intent(in)    :: G     !< ocean grid structure
    type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
    real,                    intent(in)    :: dt    !< Timestep
    real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: h
    real, dimension(SZI_(G), SZJ_(G), SZK_(G)) :: temp_t
    real :: time_start, time_end

    character(len=1024) :: key_prefix

    integer :: isc, iec, jsc, jec, nk, i, j, k
    real :: time_real

    isc = G%isc; iec = G%iec
    jsc = G%jsc; jec = G%jec
    nk = GV%ke
    time_real = time_type_to_real(Time)
    write(key_prefix,("(I6.6,A,F16.2,A)")) PE_here(),"_",time_real,"_"

    call stripspaces(key_prefix)
    if ( mod(CS%accumulated_time,86400) == 0) then
      ! Timing related quantities
      CS%num_elements = 0
      CS%timing_current = 0.

      call pass_var(h,G%Domain)
      ! Check to see whether we should actually send data
      if (CS%id_streamout_h > 0)   call put_array(CS%ssc_data_client, trim(key_prefix)//"h"  , h)
      if (CS%id_streamout_T > 0)   call put_array(CS%ssc_data_client, trim(key_prefix)//"T"  , CS%T)
      if (CS%id_streamout_S > 0)   call put_array(CS%ssc_data_client, trim(key_prefix)//"S"  , CS%S)
      if (CS%id_streamout_uh > 0)  call put_array(CS%ssc_data_client, trim(key_prefix)//"uh" , CS%uh)
      if (CS%id_streamout_vh > 0)  call put_array(CS%ssc_data_client, trim(key_prefix)//"vh" , CS%vh)
      if (CS%id_streamout_ssh > 0) call put_array(CS%ssc_data_client, trim(key_prefix)//"ssh", CS%ssh)
      CS%accumulated_time = 0
      CS%timing_total = CS%timing_total + CS%timing_current
  !    write(CS%timing_unit,*) NINT(time_real/86400.), CS%timing_total, CS%timing_current, CS%num_elements
    endif
    CS%accumulated_time = CS%accumulated_time + dt*2

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

  subroutine StripSpaces(string)
  character(len=*) :: string
  integer :: stringLen
  integer :: last, actual

  stringLen = len (string)
  last = 1
  actual = 1

  do while (actual < stringLen)
      if (string(last:last) == ' ') then
          actual = actual + 1
          string(last:last) = string(actual:actual)
          string(actual:actual) = ' '
      else
          last = last + 1
          if (actual < last) &
              actual = last
      endif
  end do

  end subroutine

  !> Wrapper for SSC fortran client put command
  subroutine put_array_here( CS, key, array )
    type(data_client_type)      :: CS      !< Data client type
    character(len=*)            :: key     !< The key used in the database
    real(kind=8), dimension(..) :: array   !< Data to be sent

    real :: time_start, time_end


    call cpu_time(time_start)
    call put_array(CS%ssc_data_client, trim(key), array)
    call cpu_time(time_end)

    CS%num_elements = CS%num_elements + size(array)
    CS%timing_current = CS%timing_current + (time_end - time_start)
  end subroutine put_array_here

end module MOM_data_client
