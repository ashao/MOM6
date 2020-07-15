module MOM_smartsim_connector

use client_fortran_api, only : put_array, init_ssc_client
use MOM_coms,           only : pe_here
use MOM_domains,        only : pass_var
use MOM_grid,           only : ocean_grid_type
use MOM_time_manager,   only : time_type, get_date
use iso_c_binding,      only : c_ptr
implicit none; private

#include <MOM_memory.h>

integer, parameter :: PE_STRING_LENGTH = 5

type, public :: smartsim_type; private
  type(c_ptr)                     :: smartsim_client    !< Pointer to an initialized Smartsim Client
  character(len=PE_STRING_LENGTH) :: unique_id          !< A string used to identify this PE
  logical                         :: send_h = .true.    !< If true, send layer thicknesses (h)
  real                            :: averaging_interval = 86400. !< How frequently to send the data in seconds
  real                            :: accumulated_time   !< Time since the last time data was sent to the database
  real ALLOCABLE_, dimension(NIMEM_,NJMEM_,NKMEM_) :: h_avg
end type smartsim_type

public :: initialize_smartsim_connector
public :: send_variables_to_database

contains

!> Initalize the SmartSim client and send any static metadata
subroutine initialize_smartsim_connector(CS, G)
  type(smartsim_type), pointer,  intent(inout) :: CS !< Control structure for the smartsim connector
  type(ocean_grid_type),         intent(in)    :: G  !< Contains grid metrics

  integer, dimension(2) :: global_start_indices
  global_start_indices(1) =  G%isd_global
  global_start_indices(2) =  G%jsd_global

  allocate(CS)
  CS%smartsim_client = init_ssc_client()
  write(CS%unique_id,'(I5.5)') pe_here()
  call put_array(CS%smartsim_client, CS%unique_id//"_rank-meta", global_start_indices)
  CS%accumulated_time = 0.

  ALLOC_(CS%h_avg(G%isd:G%ied,G%jsd:G%jed,G%ke)); CS%h_avg(:,:,:) = 0.

end subroutine initialize_smartsim_connector

!> Send the current state of temperature and salinity to the database
subroutine send_variables_to_database(CS, G, time, dt, h)
  type(smartsim_type),      intent(inout) :: CS   !< Control structure for the smartsim connectori
  type(ocean_grid_type),    intent(in   ) :: G    !< Contains grid metrics
  type(time_type),          intent(in   ) :: time !< Current time of the model
  real,                     intent(in   ) :: dt   !< Length of the timestep
  real, dimension(:,:,:),   intent(in   ) :: h    !< 3D array of layer thicknesses

  character(len=10) :: datestring
  character(len=17) :: prefix
  real :: wt
  integer :: month, day, year, minute, hour, second

  ! Update the accumulated time
  CS%accumulated_time = CS%accumulated_time + dt

  ! Update the arrays containing averages
  wt = dt/CS%averaging_interval
  if (CS%send_h) CS%h_avg(:,:,:) = CS%h_avg(:,:,:) + h(:,:,:)*wt

  ! Check whether it's time to send data to the database
  if (CS%accumulated_time == CS%averaging_interval) then

    ! Construct a key prefix based on the time and the pe rank
    ! Write the time as as string (YYYY-MM-DD)
    call get_date(Time, year, month, day, hour, minute, second)
    write(datestring,'(I4.4,A,I2.2,A,I2.2)') year, '-', month, '-', day
    write(prefix,'(A,A,A,A)') CS%unique_id, '_', datestring, '_'

    ! Send any of the requested fields
    if (CS%send_h) then
      call pass_var(CS%h_avg,G%Domain)
      call put_array(CS%smartsim_client, prefix//'h', CS%h_avg)
      CS%h_avg(:,:,:) = 0.
    endif
    CS%accumulated_time = 0
  endif
end subroutine send_variables_to_database

end module mom_smartsim_connector
