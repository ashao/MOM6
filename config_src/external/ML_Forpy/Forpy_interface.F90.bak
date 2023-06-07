module Forpy_interface

use MOM_error_handler,         only : MOM_error, WARNING
use MOM_grid,                  only : ocean_grid_type

implicit none; private

public :: forpy_run_python_init,forpy_run_python,forpy_run_python_finalize

!> Control structure for Python interface
type, public :: python_interface ; private
end type

contains

!> Initialize Forpy with specify Python script and directory
subroutine forpy_run_python_init(CS,python_dir,python_file)
    character(len=*),         intent(in)  :: python_dir   !< The directory in which python scripts are found
    character(len=*),         intent(in)  :: python_file  !< The name of the Python script to read
    type(python_interface),   intent(inout) :: CS         !< Python interface object
    call MOM_error(WARNING,"Forpy was compiled using the dummy module. If this was &
                            a mistake, please follow the instructions in: &
                            MOM6/config_src/external/ML_Forpy/README.md")

end subroutine forpy_run_python_init

!> !> Send variables to a python script and output the results
subroutine forpy_run_python(in1, out1, CS, TopLayer, G)
    type(python_interface),        intent(in)  :: CS     !< Python interface object
    type(ocean_grid_type),         intent(in)    :: G     !< The ocean's grid structure.
  ! Local Variables
    logical, intent(in) :: TopLayer             !< If true, only top layer is used.
    real, dimension(:,:,:,:), &
                                    intent(in) :: in1     ! input variables.
    real, dimension(:,:,:,:), &
                                    intent(inout) :: out1      ! output variables.

end subroutine forpy_run_python 

!> Finalize Forpy
subroutine forpy_run_python_finalize(CS)
    type(python_interface), intent(inout) :: CS !< Python interface object

end subroutine forpy_run_python_finalize

end module Forpy_interface