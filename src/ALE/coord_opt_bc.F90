!> Regrid columns with a grid spacing that is optimally defined to resolve vertical baroclinic modes
module coord_opt_bc

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, NOTE, FATAL
use MOM_remapping,     only : remapping_CS, remapping_core_h
use MOM_verticalGrid,  only : verticalGrid_type
use MOM_EOS,           only : EOS_type, calculate_density_derivs
use MOM_EOS,           only : EOS_manual_init, EOS_LINEAR
use regrid_interp,     only : interp_CS_type, build_and_interpolate_grid, DEGREE_MAX

implicit none ; private

integer, parameter, public :: OPT_BC_CHEBYSHEV = 1 !< A faux-enum for selecting a sampling method based on Chebyshev polynomials
integer, parameter, public :: OPT_BC_COSINE = 2    !< A faux-enum for selecting a sampling method based on a cosine function

!> Control structure containing required parameters for the opt_bc coordinate
type, public :: opt_bc_CS ; private

  !> Number of layers
  integer :: nk

  !> Minimum thickness allowed for layers, often in [H ~> m or kg m-2]
  real :: min_thickness = 0.

  ! TODO: Add min_N2 get_param call
  !> Minimum (positive) buoyancy frequency (s^-2, this is from Jeffrey Early's code)
  real :: min_N2 = 1.e-10
  real :: max_N2 = 1.e-4

  ! TODO: Add sample_method to get_param call
  !> Which basis to use (cosine or chebyshev)
  integer :: sample_method = OPT_BC_COSINE

  !> Interpolation control structure
  type(interp_CS_type) :: interp_CS

  real :: PI = 4.0*atan(1.0)

  ! Stewart Z* grid parameters
  logical :: hybridize_stewart !< If True, merge the Gauss-Lobatto grid with z* grid as described in
                               !! Stewart et al. [Ocean Modelling, 2017]
  real :: stewart_min_dz       !< The minimum layer thickness used to calculate the Stewart grid
  real :: stewart_max_dz       !< The maximum layer thickness used to calculate the Stewart grid
  real :: stewart_S_h          !< The shape parameter for the vertical tanh function
  real :: stewart_H_max        !< The maximum estimated boundary layer depth of the ocean
  integer :: stewart_nk        !< The number of layers that are in the calculated Stewart grid
  real, dimension(:), allocatable :: stewart_z_interface

end type opt_bc_CS

public init_coord_opt_bc, set_opt_bc_params, build_opt_bc_column
public adjust_opt_bc_surface
public coord_opt_bc_unit_tests
public end_coord_opt_bc
public create_stewart_boundary_layer_grid, merge_opt_bc_stewart, initialize_stewart_grid

contains

!> Initialise a opt_bc_CS with pointers to parameters
subroutine init_coord_opt_bc(CS, nk)
  type(opt_bc_CS),         pointer    :: CS !< Unassociated pointer to hold the control structure
  integer,              intent(in) :: nk !< Number of layers in the grid


  if (associated(CS)) call MOM_error(FATAL, "init_coord_opt_bc: CS already associated!")
  allocate(CS)

  CS%nk = nk


end subroutine init_coord_opt_bc

!> This subroutine deallocates memory in the control structure for the coord_opt_bc module
subroutine end_coord_opt_bc(CS)
  type(opt_bc_CS), pointer :: CS !< Coordinate control structure

  ! nothing to do
  if (.not. associated(CS)) return
  deallocate(CS)
end subroutine end_coord_opt_bc

subroutine initialize_stewart_grid( CS )
  type(opt_bc_CS),         pointer    :: CS !< Unassociated pointer to hold the control structure
  real :: total_depth
  real :: new_dz
  integer :: k

  if (CS%hybridize_stewart) then
    ! Calculate the Stewart grid twice. Once to get the number of layers and the second to actually
    ! store the interface heights
    total_depth = 0.
    CS%stewart_nk = 0

    do while ( total_depth < CS%stewart_H_max )
     new_dz = CS%stewart_max_dz*TANH(CS%PI*total_depth/(CS%stewart_S_h*CS%stewart_H_max)) + CS%stewart_min_dz
     total_depth = total_depth + new_dz
     CS%stewart_nk = CS%stewart_nk + 1
    end do

    if (CS%stewart_nk > CS%nk) call MOM_error(FATAL, &
      "The number of layers needed to accommodate the Stewart grid is larger\n"//&
      "than the number of vertical levels. Increase NK or modify parameters of\n"//&
      "the Stewart grid")

    ! Allocate the required size of array
    allocate(CS%stewart_z_interface(CS%stewart_nk+1))
    CS%stewart_z_interface(1) = 0.

    ! Calculate the Stewart grid again and now store the interface heights
    total_depth = 0.
    do k = 1,CS%stewart_nk
     new_dz = CS%stewart_max_dz*TANH(CS%PI*total_depth/(CS%stewart_S_h*CS%stewart_H_max)) + CS%stewart_min_dz
     total_depth = total_depth + new_dz
     CS%stewart_z_interface(k+1) = total_depth
    end do
  endif

end subroutine initialize_stewart_grid

!> This subroutine can be used to set the parameters for the coord_opt_bc module
subroutine set_opt_bc_params(CS, min_thickness, min_N2, max_N2, sample_method, &
  hybridize_stewart, stewart_min_dz, stewart_max_dz, stewart_S_h, stewart_H_max )
  type(opt_bc_CS),      pointer    :: CS !< Coordinate control structure
  real,    optional, intent(in) :: min_thickness !< Minimum allowed thickness [H ~> m or kg m-2]
  real,    optional, intent(in) :: min_N2 !< Minimum N2 allowed
  real,    optional, intent(in) :: max_N2 !< Maximum N2 allowed
  integer, optional, intent(in) :: sample_method !< The way to sample the 's' coordinate
  logical, optional, intent(in) :: hybridize_stewart !< If True, merge the Gauss-Lobatto grid with z* grid as described in
                                   !! Stewart et al. [Ocean Modelling, 2017]
  real,    optional, intent(in) :: stewart_min_dz !< The minimum layer thickness used to calculate the Stewart grid
  real,    optional, intent(in) :: stewart_max_dz !< The maximum layer thickness used to calculate the Stewart grid
  real,    optional, intent(in) :: stewart_S_h    !< The shape parameter for the vertical tanh function
  real,    optional, intent(in) :: stewart_H_max  !< The maximum estimated boundary layer depth of the ocean

  if (.not. associated(CS)) call MOM_error(FATAL, "set_opt_bc_params: CS not associated")

  if (present(min_thickness)) CS%min_thickness = min_thickness
  if (present(min_N2)) CS%min_N2 = min_N2
  if (present(max_N2)) CS%max_N2 = max_N2
  if (present(sample_method)) CS%sample_method = sample_method
  if( present(hybridize_stewart))  CS%hybridize_stewart = hybridize_stewart
  if( present(stewart_min_dz))     CS%stewart_min_dz = stewart_min_dz
  if( present(stewart_max_dz))     CS%stewart_max_dz = stewart_max_dz
  if( present(stewart_S_h   ))     CS%stewart_S_h    = stewart_S_h
  if( present(stewart_H_max ))     CS%stewart_H_max  = stewart_H_max

end subroutine set_opt_bc_params

!> Set the interface heights based on the previously calculated Stewart grid that intersects the diagnosed boundary layer depth)
subroutine create_stewart_boundary_layer_grid( CS, GV, boundary_layer_depth, nk_boundary_layer )
  type(opt_bc_CS),         intent(in)    :: CS !< coord_opt_bc control structure
  type(verticalGrid_type), intent(in)    :: GV !< Vertical grid structure
  real,                    intent(in)    :: boundary_layer_depth !< The depth of the boundary layer
  integer,                 intent(  out) :: nk_boundary_layer

  real :: total_depth
  integer :: k

  if (CS%hybridize_stewart .and. allocated(CS%stewart_z_interface)) then
    nk_boundary_layer = 0
    do k=1,CS%stewart_nk+1
      if (CS%stewart_z_interface(k) < boundary_layer_depth) then
        nk_boundary_layer = nk_boundary_layer + 1
      else
        exit
      endif
    enddo
  endif

end subroutine create_stewart_boundary_layer_grid

!> Build a opt_bc coordinate column
!!
!! 1. Density profiles are calculated on the source grid.
!! 2. Positions of target densities (for interfaces) are found by interpolation.
subroutine build_opt_bc_column(CS, GV, nz, nk_boundary_layer, h, T, S, eta_orig, z_interface, EOS)
  type(opt_bc_CS),        intent(in)    :: CS !< coord_opt_bc control structure
  type(verticalGrid_type), intent(in) :: GV !< Vertical grid structure
  integer,             intent(in)    :: nz !< Number of levels on source grid (i.e. length of  h, T, S)
  integer,             intent(in)    :: nk_boundary_layer !< Number of layers in the boundary layer
  real, dimension(nz), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(nz), intent(in)    :: T  !< Temperature for source column [degC]
  real, dimension(nz), intent(in)    :: S  !< Salinity for source column [ppt]
  real, dimension(nz+1), intent(in)  :: eta_orig !< Absolute positions of interfaces
  real, dimension(CS%nk-nk_boundary_layer+1), &
                       intent(inout) :: z_interface !< Absolute positions of interfaces
  type(EOS_type),      pointer       :: EOS !< Control structure for equation of state

  ! Local variables
  integer :: k, kidx
  integer :: k_thin
  real, dimension(nz+1) :: pres     ! Pressures used to calculate density [R L2 T-2 ~> Pa]
  real, dimension(nz+1) :: T_int, S_int! Derivatives of density
  real, dimension(nz+1) :: drho_dT, drho_dS ! Derivatives of density
  real, dimension(nz+1) :: N2, N
  real, dimension(nz+1) :: s_unscaled
  real, dimension(CS%nk-nk_boundary_layer+1) :: s_scaled
  real, dimension(nz+1) :: I_dz_int
  real, dimension(CS%nk-nk_boundary_layer) :: h_new ! New thicknesses [H ~> m or kg m-2]
  real :: I_nk


  real :: z0_top, eta ! Thicknesses or heights [Z ~> m] or [H ~> m or kg m-2]
  real :: range_s

  ! For now, only do the calculation of the Chebyshev grid over all layers
  ! Revisit this calculation later to see if we need to be a bit more careful
  ! in how to deal with vanished layers

  ! Calculate pressure, temperature, and salinity at layer interfaces.
  ! For now, simply average the layer quantites for T/S, later we could
  ! support higher order interpolations use the other ALE code
  pres(1) = 0.
  T_int(1) = T(1)
  S_int(1) = S(1)
  do k=2,nz
    pres(k) = pres(k-1) + GV%H_to_Pa*h(k)
    T_int(k) = 0.5*(T(k-1)+T(k))
    S_int(k) = 0.5*(S(k-1)+S(k))
  enddo
  T_int(nz+1) = T(nz)
  S_int(nz+1) = S(nz)
  pres(nz+1) = pres(nz)+h(nz)*GV%H_to_Pa

  ! Original code has a limiter for a minimum N^2, but not a maximum
  ! value which might be needed in the limit of thin, but extant layers.

  call calculate_density_derivs(T_int, S_int, pres, drho_dT, drho_dS, 2, nz, EOS)

  !   Set up I_dz_int as the inverse of the distance between
  ! adjacent layer centers.
  I_dz_int(1) = 2.0 / (h(1)*GV%H_to_Z)
  do K=2,nz
    I_dz_int(K) = 2.0 / (GV%H_to_Z*(h(k-1) + h(k)))
  enddo
  I_dz_int(nz+1) = 2.0 / (GV%H_to_Z*h(nz))

  N(1) = SQRT(CS%min_N2) ; N(nz+1) = SQRT(CS%min_N2)
  do K=2,nz
    N(k) = SQRT(MAX((CS%min_N2), (GV%g_Earth/GV%Rho0)*(-I_dz_int(K) * (drho_dT(K) * (T(k-1)-T(k)) + drho_dS(K) * (S(k-1)-S(k))))))
    N(k) = MIN(N(k), SQRT(CS%max_N2))
  enddo

  ! Calculate the stretched coordinate from the bottom up s.t. s = 0 corresponds to the bottom
  s_unscaled(nz+1) = 0.
  do K=nz,1,-1
    s_unscaled(k) = s_unscaled(k+1) + (0.5*h(k))*(N(k)+N(k+1))
  enddo

  range_s = MAXVAL(s_unscaled(1:nz))
  I_nk = 1./(CS%nk)

  ! Gauss-Lobatto grid
  select case (CS%sample_method)
    case(OPT_BC_COSINE)
      kidx = 0
      do k=CS%nk-nk_boundary_layer,0,-1
        kidx = kidx + 1
        s_scaled(kidx) = range_s*(k*I_nk)
      enddo
    case(OPT_BC_CHEBYSHEV)
      do K=0,CS%nk-nk_boundary_layer
        s_scaled(k+1) = (0.5*range_s)*(cos( CS%pi*(K*I_nk)) + 1)
      enddo
  end select

  ! Interpolate from eta on the scaled grid back to eta on the s-grid
  z_interface(:) = interp_linear( s_unscaled, eta_orig, s_scaled )
  z_interface(CS%nk-nk_boundary_layer+1) = eta_orig(nz+1)
  z_interface(1) = eta_orig(1)

end subroutine build_opt_bc_column

!> Merge the z*-like stewart grids with the previously computed Gauss-Lobatto grid
subroutine merge_opt_bc_stewart( CS, GV, z_bottom, z_gl, nk_st, z_interface_out)
  type(opt_bc_CS),          intent(in) :: CS !< coord_opt_bc control structure
  type(verticalGrid_type),  intent(in) :: GV !< Vertical grid structure
  real,                     intent(in) :: z_bottom !< The depth of the column
  real, dimension(CS%nk+1), intent(in) :: z_gl !< The interfaces from the Gauss-Lobatto grid
  integer,                  intent(in) :: nk_st !< The number of layers of the Stewart grid within the boundary layer
  real, dimension(GV%ke+1), intent(  out) :: z_interface_out !< The merged Stewart and Gauss-Lobatto grids

  integer :: k_gl, k_st, k

  k_gl = 2
  k_st = 2

  if (CS%hybridize_stewart) then
    z_interface_out(1) = 0.
    do k=2,GV%ke+1
      if (k_st <= nk_st+1) then
        if ( min(CS%stewart_z_interface(k_st),z_bottom) <= z_gl(k_gl)) then
          z_interface_out(k) = min(CS%stewart_z_interface(k_st), z_bottom)
          k_st = k_st + 1
        else
          z_interface_out(k) = z_gl(k_gl)
          k_gl = k_gl + 1
        endif
      else
        z_interface_out(k) = z_gl(k_gl)
        k_gl = k_gl + 1
      endif
    enddo
  endif
  ! Manual unit test
  !  st    gl
  !  [0    0]
  !  [1    2]
  !  [3    3]
  !  [4    4]
  !        5
  !        6
  !
  !k_st = 2 k_gl = 2
  !z_interface(2) = 1
  !k_st = 3 k_gl = 2
  !z_interface(3) = 2
  !k_st = 3 k_gl = 3
  !z_interface(4) = 3
  !k_st = 3 k_gl = 4
  !z_interface(5) = 3
  !k_st = 4 k_gl = 4
  !z_interface(6) = 4
  !k_st = 4 k_gl = 5
  !z_interface(7) = 4
  !k_st = 5 k_gl = 5
  !z_interface(8) = 5
  !k_st = 5 k_gl = 6
  !z_interface(9) = 6

end subroutine merge_opt_bc_stewart

!> Take care of the cases where the interpolated interfaces are less than the desired minimum thickness
!! by linearly spacing interfaces at the top and bottom.
subroutine adjust_opt_bc_surface( CS, nk, h )
  type(opt_bc_CS),        intent(in)    :: CS !< coord_opt_bc control structure
  integer,                intent(in)    :: nk !< Number of levels on source grid (i.e. length of  h, T, S)
  real, dimension(nk),    intent(inout) :: h  !< The thicknesses to be adjusted

  integer :: k, k_thin
  real :: htot

  htot = 0.
  k_thin = 1

  do k = 1,CS%nk
    htot = htot + h(k)
    k_thin = k
    if ( htot >= k*CS%min_thickness ) then
      exit
    endif
  enddo

  do k=1,k_thin
    h(k) = htot/k_thin
  enddo

end subroutine adjust_opt_bc_surface

!> Perform simple tests of the opt_bc/cosine scaled coordinates
function coord_opt_bc_unit_tests( verbose ) result( failed )
  logical, intent(in) :: verbose !< If true, print information for each unit test
  logical :: failed !< If true, one of the unit tests has failed

  type(EOS_type), pointer :: EOS
  type(opt_bc_CS) :: CS
  type(verticalGrid_type) :: GV

  integer, parameter :: nk = 4
  real, dimension(nk) :: temp, salt, h
  real, dimension(nk+1) :: eta_orig, eta_new, eta_expected

  ! Initialize basic variables in this unit test
  failed = .false.
  CS%nk = nk
  allocate(EOS)
  call EOS_manual_init( EOS, form_of_EOS = EOS_linear, dRho_dT = -1., dRho_dS = 0., Compressible = .false. )
  GV%ke = nk
  GV%H_to_Pa = 10000. ! Assumes that H is in meters with 1m = 1dbar
  GV%H_to_Z = 1.
  GV%g_Earth = 1.
  GV%Rho0 = 1.

  ! Test 1: N2 = 1, cosine basis
  CS%sample_method =  OPT_BC_COSINE
  temp(:) = [4., 3., 2., 1.]
  salt(:) = [0., 0., 0., 0.]
  h(:)    = [1., 1., 1., 1.]
  eta_orig(:) = [0., 1., 2., 3., 4.]
  eta_expected(:) = [0., 1.25, 2., 2.75, 4.]

  call build_opt_bc_column( CS, GV, nk, 0, h, temp, salt, eta_orig, eta_new, EOS)
  failed = any( ABS(eta_orig - eta_expected) > 5.e-3 )
  if (failed) call MOM_error( NOTE, "FAILED: coord_opt_bc N2 = 1, cosine basis" )

  ! Test 2: N2 = 1, chebyshev basis
  CS%sample_method = OPT_BC_CHEBYSHEV
  eta_expected = [0., 0.8786738, 2., 3.121326, 4.]
  call build_opt_bc_column( CS, GV, nk, 0, h, temp, salt, eta_orig, eta_new, EOS)
  failed = any( ABS(eta_orig - eta_expected) > 5.e-3 )
  if (failed) call MOM_error( NOTE, "FAILED: coord_opt_bc N2 = 1, Chebyshev basis" )

end function coord_opt_bc_unit_tests

!> Interpolate a vector of points
function interp_linear( xin, fin, xout ) result( fout )
  real, dimension(:) :: xin !< Coordinates of the original function
  real, dimension(:) :: fin !< The values of the original function
  real, dimension(:) :: xout !< The points at which to interpolate
  real, dimension(size(xout)) :: fout !< The interpolated function

  integer :: nin, nout
  integer :: kout, kin, kin1, kin2
  logical :: xin_increases
  real :: wt

  nin = size(xin)
  nout = size(xout)

  xin_increases = xin(1) < xin(nin)

  do kout = 1,nout

    ! If out of bounds, then assume constant extrapolation
    if (xin_increases) then
      if (xout(kout) >= xin(nin)) then
        fout(kout) = fin(nin)
        cycle
      elseif (xout(kout) <= xin(1)) then
        fout(kout) = fin(1)
        cycle
      endif
    else
      if (xout(kout) <= xin(nin)) then
        fout(kout) = fin(nin)
        cycle
      elseif (xout(kout) >= xin(1)) then
        fout(kout) = fin(1)
        cycle
      endif
    endif
    ! Otherwise we can find the interval to interpolate between
    call find_interval( xin, xout(kout), kin1, kin2 )
    wt = ( xout(kout) - xin(kin1) )/(xin(kin2) - xin(kin1))
    fout(kout) = (1.-wt)*fin(kin1) + wt*fin(kin2)

  enddo

end function interp_linear

!> Find the indices k1 and k2 of xin in which xout lies
subroutine find_interval( xin, xout, k1, k2 )
  real, dimension(:), intent(in   ) :: xin  !< The values of the input coordinate
  real,               intent(in   ) :: xout !< The coordinate to find the interval
  integer,            intent(  out) :: k1   !< The index of the left bound
  integer,            intent(  out) :: k2   !< The index of the right bound

  integer :: nin, k

  logical :: x_increases

  nin = size(xin)

  x_increases = xin(1) < xin(nin)

  ! First capture the out of bound cases
  if ( xout > MAXVAL(xin) ) then
    if (x_increases) then
      k1 = nin
      k2 = nin + 1
    else
      k1 = 0
      k2 = 1
    endif
  elseif ( xout < MINVAL(xin) ) then
    if (x_increases) then
      k1 = 0
      k2 = 1
    else
      k1 = nin
      k2 = nin + 1
    endif
  else ! Find the first interval that xout exists
    do k=1,nin-1
      if ( ((xout >= xin(k)) .and. (xout < xin(k+1))) .or. &
           ((xout <= xin(k)) .and. (xout > xin(k+1))) ) then
            k1 = k
            k2 = k+1
            exit
      endif
    enddo
  endif

end subroutine find_interval


end module coord_opt_bc
