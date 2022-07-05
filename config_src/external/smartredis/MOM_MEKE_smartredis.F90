!> Contains routines that contain dummy routines for the smart
module MOM_MEKE_smartredis

use MOM_diag_mediator,     only : diag_ctrl
use MOM_error_handler,     only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_grid,              only : ocean_grid_type
use MOM_file_parser,       only : param_file_type
use MOM_smartredis,        only : smartredis_CS_type
use MOM_unit_scaling,      only : unit_scale_type
use MOM_variables,         only : thermo_var_ptrs
use MOM_verticalGrid,      only : verticalGrid_type

implicit none; private

#include <MOM_memory.h>

public meke_smartredis_init, infer_meke

type, public :: meke_smartredis_CS_type; private

end type meke_smartredis_CS_type

contains

!> Initializer for the SmartRedis MEKE module that uses ML to predict eddy kinetic energy
subroutine meke_smartredis_init(diag, G, US, param_file, smartredis_CS, CS)
  type(diag_ctrl), target, intent(inout) :: diag       !< Diagnostics structure.
  type(ocean_grid_type),         intent(inout) :: G          !< The ocean's grid structure.
  type(unit_scale_type),         intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),         intent(in)    :: param_file !< Parameter file parser structure.
  type(smartredis_CS_type),      intent(in)    :: smartredis_CS !< SmartRedis client
  type(meke_smartredis_CS_type), intent(inout) :: CS         !< Control structure for this module

  call MOM_error(FATAL,"meke_smartredis_init was compiled using the dummy module. Recompile"//&
                       "with source code from https://github.com/CrayLabs/MOM6-smartredis")
end subroutine meke_smartredis_init

!> Use the SmartRedis client to call a machine learning to predict eddy kinetic energy
subroutine infer_meke(G, GV, MEKE, u, v, tv, h, dt, CS)
  type(ocean_grid_type),                     intent(inout) :: G  !< Ocean grid
  type(verticalGrid_type),                   intent(in)    :: GV !< Ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G)), intent(  out) :: MEKE !< Vertically averaged eddy kinetic energy [L2 T-2 ~> m2 s-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(G)), intent(inout) :: u  !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G),SZJB_(G),SZK_(G)), intent(inout) :: v  !< Meridional velocity [L T-1 ~> m s-1]
  type(thermo_var_ptrs),                     intent(in)    :: tv !< Type containing thermodynamic variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  real,                                      intent(in)    :: dt !< Model(baroclinic) time-step [T ~> s].
  type(meke_smartredis_CS_type),             intent(in)    :: CS !< Control structure for inferring MEKE using SmartRedis

  call MOM_error(FATAL,"infer_meke was compiled using the dummy module. Recompile"//&
                       "with source code from https://github.com/CrayLabs/MOM6-smartredis")

end subroutine infer_meke

end module MOM_MEKE_smartredis