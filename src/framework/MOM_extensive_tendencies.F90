!> Provides buffers that can dynamically grow as needed. These are primarily intended for the
!! diagnostics which need to store intermediate or partial states of state variables
module MOM_extensive_tendencies
! This file is part of MOM6. See LICENSE.md for the license.

use iso_c_binding,     only : c_ptr, c_loc
use MOM_diag_buffers,  only : diag_buffer_3d
use MOM_diag_mediator, only : diag_ctrl, post_data
use MOM_grid,          only : ocean_grid_type
use MOM_verticalGrid,  only : verticalGrid_type

implicit none ; private

!> A base type to use for extensive tendency remapping
type, abstract :: extensive_tendency_remapper_base ; private
  type(diag_buffer_3d) :: h_ante !< Stores the thickness fields before a process [H ~> m or kg m-2]

  contains

  procedure, public :: prepare_tendency_calculation !< Perform tasks needed prior to calculating tendencies
  procedure, public :: finalize_tendency_calculation !< Perform tasks need to finalize the calculation of tendencies
  procedure :: hash_string !< Hash a string to an integer number
  procedure(a_calculate_tendency), deferred : calculate_tendency
end type extensive_tendency_remapper_base

!> Implements the "original" algorithm for remapping tendencies where all tendencies are remapped to the
!! vertical grid prior to the process
type, extends(extensive_tendency_remapper_base) ::extensive_tendency_remapper_original ; private

  contains

  procedure :: calculate_tendency => calculate_tendency_original

end type extensive_tendency_remapper_original

contains

!> Calculate and remap a tendency using the new (correct) algorithm. Compute

!! before the process occurs. Some logic in here controls how the tendency is actually calculated.
subroutine post_tendency_new(this, G, GV, diag_cs, tag, coeff, h_ante, h_post, field_ante, &
                                  field_post, diag_field_id)
  class(extensive_tendency_remapper_original), intent(in) :: this !< The tendency remapper
  type(ocean_grid_type),                       intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in) :: GV !< The ocean's vertical grid structure.
  type(diag_ctrl),                             intent(in) :: diag_cs !< Diagnostic control structure
  character(len=*), intent(in) :: tag   !< The tag used to identify the process causing a change
  real,                                        intent(in) :: coeff !< A coefficient used to multiply the
                                                                   !! subtracted fields, usually Idt but could
                                                                   !! be something like a heat capacity [arbitrary]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target, &
                    intent(in) :: h_ante     !< The vertical grid before the process [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target, &
                    intent(in) :: h_post     !< The vertical grid after the process [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                    intent(in) :: field_ante !< The extensive quantity before the process occurred [arbitrary]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                    intent(in) :: field_post !< The extensive quantity after the process occurred [arbitrary]
  integer,          intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                !! previous call to register_diag_field.


  integer :: i, j, k, is, ie, js, je, nz, ax_idx
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: extensive_ante, extensive_post, dz_ante, dz_post
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV),diag_CS%max_ke) :: remapped_ante, remapped_post, tendency
  logical :: compute_dz
  logical :: is_h_tendency

  is = G%isc; ie = G%iec; js = G%jsc; je = G%jec; nz = GV%ke

  !! Calculate the extensive field before and after the process
  ! Check if the thickness tendency is being calculated by seeing if the h and field arrays are the same
  is_h_tendency = (c_loc(h_ante) == c_loc(field_ante)) .and. (c_loc(h_post) == c_loc(field_post))

  if(is_h_tendency) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      extensive_ante(i,j,k) = coeff*h_ante(i,j,k)
      extensive_post(i,j,k) = coeff*h_post(i,j,k)
    enddo; enddo; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      extensive_ante(i,j,k) = (coeff*h_ante(i,j,k))*field_ante(i,j,k)
      extensive_post(i,j,k) = (coeff*h_post(i,j,k))*field_post(i,j,k)
    enddo; enddo; enddo
  endif

  !! Post the results
  diag => diag_cs%diags(diag_field_id)
  ! Check to see if we are remapping to a z-based coordinate
  do while (associated(diag)) then
    if(.not. diag%axes%is_native) then
      call thickness_to_dz(h_ante, diag_cs%tv, dz_ante, G, GV, diag_cs%US, halo_size=1)
      call thickness_to_dz(h_post, diag_cs%tv, dz_post, G, GV, diag_cs%US, halo_size=1)
      exit
    endif
    diag => diag%next
  enddo

  diag => diag_cs%diags(diag_field_id)

  do while (associated(diag))
    ! No reintegration needed for native grid and can just be directly posted
    if (diag%axes%is_native) then
      nz = GV%ke
      tendency(:,:,1:nz) = (extensive_post(:,:,:) - extensive_ante(:,:,:))
      call post_data_3d_low(diag, tendency(:,:,1:nz), diag_cs, .false., diag%axes%mask3d)
    else
      ax_idx = diag%axes%vertical_coordinate_number
      nz = diag%axes%nz
      if (diag_cs%diag_remap_cs(ax_idx)%Z_based_coord) then
        call vertically_reintegrate_diag_field(                                    &
                diag_cs%diag_remap_cs(ax_idx), diag_cs%G, &
                dz_ante, diag_cs%diag_remap_cs(ax_idx)%h_ante, & ! TODO: Update this with diag_buffer for ante-grid
                .false., .false., diag%axes%mask3d, extensive_ante, remapped_ante(:,:,1:nz))
        call vertically_reintegrate_diag_field(                                    &
                diag_cs%diag_remap_cs(ax_idx), diag_cs%G, &
                dz_post, diag_cs%diag_remap_cs(ax_idx)%h_post, & ! TODO: Update this with diag_buffer for post-grid
                .false., .false., diag%axes%mask3d, extensive_post, remapped_post(:,:,1:nz))
      else
        call vertically_reintegrate_diag_field(                                    &
                diag_cs%diag_remap_cs(ax_idx), diag_cs%G, &
                h_ante, diag_cs%diag_remap_cs(ax_idx)%h_ante, & ! TODO: Update this with diag_buffer for ante-grid
                .false., .false., diag%axes%mask3d, extensive_ante, remapped_ante(:,:,1:nz))
        call vertically_reintegrate_diag_field(                                    &
                diag_cs%diag_remap_cs(ax_idx), diag_cs%G, &
                h_post, diag_cs%diag_remap_cs(ax_idx)%h_post, & ! TODO: Update this with diag_buffer for post-grid
                .false., .false., diag%axes%mask3d, extensive_post, remapped_post(:,:,1:nz))
      endif
      tendency(:,:,1:nz) = (remapped_post(:,:,1:nz) - remapped_ante(:,:,1:nz))
      if (associated(diag%axes%mask3d)) then
        call post_data_3d_low(diag, tendency_remapped(:,:,1:nz), diag_cs, is_static, &
                              mask=diag%axes%mask3d)
      else
        call post_data_3d_low(diag, tendency_remapped(:,:,1:nz), diag_cs, is_static)
      endif
    endif
    diag => diag%next
  enddo
end subroutine post_tendency_new

!> Calculate and remap a tendency using the original algorithm. Only remap to the vertical grid
!! before the process occurs. Some logic in here controls how the tendency is actually calculated.
subroutine post_tendency_original(this, G, GV, diag_cs, tag, coeff, h_ante, h_post, field_ante, &
                                  field_post, diag_field_id)
  class(extensive_tendency_remapper_original), intent(in) :: this !< The tendency remapper
  type(ocean_grid_type),                       intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                     intent(in) :: GV !< The ocean's vertical grid structure.
  type(diag_ctrl),                             intent(in) :: diag_cs !< Diagnostic control structure
  character(len=*), intent(in) :: tag   !< The tag used to identify the process causing a change
  real,                                        intent(in) :: coeff !< A coefficient used to multiply the
                                                                   !! subtracted fields, usually Idt but could
                                                                   !! be something like a heat capacity [arbitrary]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target, &
                    intent(in) :: h_ante     !< The vertical grid before the process [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), target, &
                    intent(in) :: h_post     !< The vertical grid after the process [H ~> m or kg m-2] or [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                    intent(in) :: field_ante !< The extensive quantity before the process occurred [arbitrary]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                    intent(in) :: field_post !< The extensive quantity after the process occurred [arbitrary]
  integer,          intent(in) :: diag_field_id !< The id for an output variable returned by a
                                                !! previous call to register_diag_field.


  integer :: i, j, k, is, ie, js, je, nz, ax_idx
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: tendency_native, dz_ante
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV),diag_CS%max_ke) :: remapped_field
  logical :: compute_dz
  logical :: h_unchanged, is_h_tendency

  is = G%isc; ie = G%iec; js = G%jsc; je = G%jec; nz = GV%ke

  !! Calculate the tendency on the native grid
  ! Check if h does not change by comparing whether the ante/post h field had the same input arguments
  h_unchanged = c_loc(h_ante) == c_loc(h_post)
  ! Check if the thickness tendency is being calculated by seeing if the h and field arrays are the same
  is_h_tendency = (c_loc(h_ante) == c_loc(field_ante)) .and. (c_loc(h_post) == c_loc(field_post))

  if(h_unchanged) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      tendency_native(i,j,k) = (coeff*h_ante(i,j,k))*(field_post(i,j,k)-field_ante(i,j,k))
    enddo; enddo; enddo
  elseif(is_h_tendency) then
    do k=1,nz ; do j=js,je ; do i=is,ie
      tendency_native(i,j,k) = coeff*(h_post(i,j,k)-h_ante(i,j,k))
    enddo; enddo; enddo
  else
    do k=1,nz ; do j=js,je ; do i=is,ie
      tendency_native(i,j,k) = coeff*(field_post(i,j,k)*h_post(i,j,k)-field_ante(i,j,k)*h_ante(i,j,k))
    enddo; enddo; enddo
  endif

  !! Post the results
  diag => diag_cs%diags(diag_field_id)
  ! Check to see if we are remapping to a z-based coordinate
  do while (associated(diag)) then
    if(.not. diag%axes%is_native) then
      call thickness_to_dz(h_ante, diag_cs%tv, dz_ante, G, GV, diag_cs%US, halo_size=1)
      exit
    endif
    diag => diag%next
  enddo

  diag => diag_cs%diags(diag_field_id)
  do while (associated(diag))
    ! No reintegration needed for native grid and can just be directly posted
    if (diag%axes%is_native) then
      call post_data_3d_low(diag, tendency_native, diag_cs, .false., diag%axes%mask3d)
    else
      ax_idx = diag%axes%vertical_coordinate_number
      nz = diag%axes%nz
      if (diag_cs%diag_remap_cs(ax_idx)%Z_based_coord) then
        call vertically_reintegrate_diag_field(                                    &
                diag_cs%diag_remap_cs(ax_idx), diag_cs%G, &
                dz_ante, diag_cs%diag_remap_cs(ax_idx)%h_extensive, &
                .false., .false., diag%axes%mask3d, tendency_native, tendency_remapped(:,:,1:nz))
      else
        call vertically_reintegrate_diag_field(                                    &
                diag_cs%diag_remap_cs(ax_idx), diag_cs%G, &
                h_ante, diag_cs%diag_remap_cs(ax_idx)%h_extensive, &
                .false., .false., diag%axes%mask3d, tendency_native, tendency_remapped(:,:,1:nz))
      endif
      if (associated(diag%axes%mask3d)) then
        call post_data_3d_low(diag, tendency_remapped(:,:,1:nz), diag_cs, is_static, &
                              mask=diag%axes%mask3d)
      else
        call post_data_3d_low(diag, tendency_remapped(:,:,1:nz), diag_cs, is_static)
      endif
    endif
    diag => diag%next
  enddo
end subroutine post_tendency_original

!> Prepare for calculating a series of extensive tendencies due to a particular process
subroutine prepare_tendency_calculation(this, G, GV, h, tag)
  class(extensive_tendency_remapper_base), intent(inout) :: this !< The tendency remapper
  type(ocean_grid_type),                   intent(in   ) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                 intent(in   ) :: GV   !< The ocean's vertical grid structure.
  character(len=*), intent(in) :: tag !< The tag used to identify the process causing a change
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                    intent(in) :: h   !< Model vertical grid (thicknesses) [H ~> m or kg m-2] or [Z ~> m]
                                      !! used to define the vertical extents of cells
  integer :: hash, slot

  hash = this%hash_string(tag)
  slot = this%h_ante%check_capacity_by_id(hash)
  call this%h_ante(h, tag)
end subroutine prepare_tendency_calculation

!> Finalize the calculation of extensive tendencies due to a particular process
subroutine finalize_tendency_calculation(this, tag)
  class(extensive_tendency_remapper_base), intent(inout) :: this !< The tendency remapper
  character(len=*), intent(in) :: tag !< The tag used to identify the process causing a change

  integer :: hash
  hash = this%hash_string(tag)
  call this%h_ante%mark_available(hash)
end subroutine finalize_tendency_calculation

!> Hash a string following Fowler–Noll–Vo hash function
pure function hash_string(this, string) result(hash)
  class(extensive_tendency_remapper_base), intent(in) :: this !< The tendency remapper
  character(len=*), intent(in) :: string !< The string to be hashed
  integer :: hash
  integer :: i

  ! The following are valid for a 32-bit hash
  integer, parameter :: offset = 2166136261
  integer, parameter :: prime = 16777619

  hash = offset_basis

  do i = 1,len_trim(str)
      hash = ieor(hash, ichar(str(i:i)))
      hash = hash * FNV_prime
  end do

end function hash_string

end module MOM_extensive_tendencies

