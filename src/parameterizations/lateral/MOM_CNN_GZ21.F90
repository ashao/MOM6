module MOM_CNN_GZ21

use MOM_grid,                  only : ocean_grid_type
use MOM_verticalGrid,          only : verticalGrid_type
use MOM_domains,               only : clone_MOM_domain,MOM_domain_type
use MOM_domains,               only : pass_vector,pass_var,BGRID_NE
use MOM_domains,               only : create_group_pass,do_group_pass,group_pass_type
use MOM_domains,               only : To_North, To_East
use MOM_coms,                  only : reproducing_sum
use MOM_diag_mediator,         only : post_data,register_diag_field
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_unit_scaling,          only : unit_scale_type
use MOM_file_parser,           only : get_param,param_file_type
use MOM_string_functions,      only : lowercase
use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
!use Forpy_interface,           only : forpy_run_python, python_interface
use SmartSim_interface,        only : smartsim_run_python, smartsim_python_interface


implicit none; private

#include <MOM_memory.h>
#ifdef STATIC_MEMORY_
#  ifndef BTHALO_
#    define BTHALO_ 0
#  endif
#  define WHALOI_ MAX(BTHALO_-NIHALO_,0)
#  define WHALOJ_ MAX(BTHALO_-NJHALO_,0)
#  define NIMEMW_   1-WHALOI_:NIMEM_+WHALOI_
#  define NJMEMW_   1-WHALOJ_:NJMEM_+WHALOJ_
#  define NIMEMBW_  -WHALOI_:NIMEM_+WHALOI_
#  define NJMEMBW_  -WHALOJ_:NJMEM_+WHALOJ_
#  define SZIW_(G)  NIMEMW_
#  define SZJW_(G)  NJMEMW_
#  define SZIBW_(G) NIMEMBW_
#  define SZJBW_(G) NJMEMBW_
#else
#  define NIMEMW_   :
#  define NJMEMW_   :
#  define NIMEMBW_  :
#  define NJMEMBW_  :
#  define SZIW_(G)  G%isdw:G%iedw
#  define SZJW_(G)  G%jsdw:G%jedw
#  define SZIBW_(G) G%isdw-1:G%iedw
#  define SZJBW_(G) G%jsdw-1:G%jedw
#endif

public :: CNN_init,CNN_inference

!> Control structure for CNN
type, public :: CNN_CS ; private
  type(MOM_domain_type), pointer :: CNN_Domain => NULL()  !< Domain for inputs/outputs for CNN
  integer :: isdw !< The lower i-memory limit for the wide halo arrays.
  integer :: iedw !< The upper i-memory limit for the wide halo arrays.
  integer :: jsdw !< The lower j-memory limit for the wide halo arrays.
  integer :: jedw !< The upper j-memory limit for the wide halo arrays.
  integer :: CNN_halo_size  !< Halo size at each side of subdomains

  logical :: CNN_BT    !< If true, momentum forcing from CNN is barotropic.

  type(diag_ctrl), pointer :: diag => NULL() !< A type that regulates diagnostics output
  !>@{ Diagnostic handles
  integer :: id_CNNu = -1, id_CNNv = -1, id_KE_CNN = -1
  integer :: id_Sxmean = -1, id_Symean = -1
  integer :: id_Sxstd = -1, id_Systd = -1
  !>@}

  ! Clock ids
  integer :: id_cnn_pre
  integer :: id_cnn_inference   !< Clock id to time initialization of the client
  integer :: id_cnn_post
  integer :: id_cnn_post1
  integer :: id_cnn_post2
  integer :: id_cnn_post3
  integer :: id_cnn_post4
end type CNN_CS

contains

!> Prepare CNN input variables with wide halos
subroutine CNN_init(Time,G,GV,US,param_file,diag,CS)
  type(time_type),               intent(in)    :: Time       !< The current model time.
  type(ocean_grid_type),         intent(in)    :: G     !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)    :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),         intent(in)    :: US         !< A dimensional unit scaling type
  type(param_file_type),         intent(in)    :: param_file !< Parameter file parser structure.
  type(diag_ctrl), target,       intent(inout) :: diag  !< Diagnostics structure.
  type(CNN_CS),                  intent(inout) :: CS    !< Control structure for CNN
  ! Local Variables
  integer :: wd_halos(2) ! Varies with CNN
  character(len=40)  :: mdl = "MOM_CNN"  ! module name

  ! Register fields for output from this module.
  CS%diag => diag

  CS%id_CNNu = register_diag_field('ocean_model', 'CNNu', diag%axesCuL, Time, &
      'Zonal Acceleration from CNN model', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_CNNv = register_diag_field('ocean_model', 'CNNv', diag%axesCvL, Time, &
      'Meridional Acceleration from CNN model', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_KE_CNN = register_diag_field('ocean_model', 'KE_CNN', diag%axesTL, Time, &
      'Kinetic Energy Source from Horizontal Viscosity of CNN model', &
      'm3 s-3', conversion=GV%H_to_m*(US%L_T_to_m_s**2)*US%s_to_T)
  CS%id_Sxmean = register_diag_field('ocean_model', 'Sxmean', diag%axesTL, Time, &
      'Zonal Acceleration from CNN model mean part', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_Symean = register_diag_field('ocean_model', 'Symean', diag%axesTL, Time, &
      'Meridional Acceleration from CNN model mean part', 'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_Sxstd = register_diag_field('ocean_model', 'Sxstd', diag%axesTL, Time, &
      'Zonal Acceleration from CNN model standard deviation part', &
      'm s-2', conversion=US%L_T2_to_m_s2)
  CS%id_Systd = register_diag_field('ocean_model', 'Systd', diag%axesTL, Time, &
      'Meridional Acceleration from CNN model standard deviation part', &
      'm s-2', conversion=US%L_T2_to_m_s2)
      
  call get_param(param_file, mdl, "CNN_BT", CS%CNN_BT, &
      "If true, momentum forcing from CNN is barotropic, otherwise baroclinic (default).", &
      default=.false.)
  call get_param(param_file, mdl, "CNN_HALO_SIZE", CS%CNN_halo_size, &
      "Halo size at each side of subdomains, depends on CNN architecture.", & 
      units="nondim", default=10)

  wd_halos(1) = CS%CNN_halo_size
  wd_halos(2) = CS%CNN_halo_size
  call clone_MOM_domain(G%Domain, CS%CNN_Domain, min_halo=wd_halos, symmetric=.true.)
  CS%isdw = G%isc-wd_halos(1) ; CS%iedw = G%iec+wd_halos(1)
  CS%jsdw = G%jsc-wd_halos(2) ; CS%jedw = G%jec+wd_halos(2)

  ! Set various clock ids
  CS%id_cnn_pre         = cpu_clock_id('(CNN before inference)', grain=CLOCK_ROUTINE)
  CS%id_cnn_inference   = cpu_clock_id('(CNN total inference)', grain=CLOCK_ROUTINE)
  CS%id_cnn_post        = cpu_clock_id('(CNN after inference)', grain=CLOCK_ROUTINE)
  CS%id_cnn_post1       = cpu_clock_id('(CNN after inference 1)', grain=CLOCK_ROUTINE)
  CS%id_cnn_post2       = cpu_clock_id('(CNN after inference 2)', grain=CLOCK_ROUTINE)
  CS%id_cnn_post3       = cpu_clock_id('(CNN after inference 3)', grain=CLOCK_ROUTINE)
  CS%id_cnn_post4       = cpu_clock_id('(CNN after inference 4)', grain=CLOCK_ROUTINE)

end subroutine CNN_init

!> Manage input and output of CNN model
subroutine CNN_inference(u, v, h, diffu, diffv, G, GV, FP_CS, SS_CS, CNN, python_bridge_lib)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(smartsim_python_interface),        intent(in)  :: FP_CS  !< Forpy Python interface object
  type(smartsim_python_interface),intent(in)  :: SS_CS  !< SmartSim Python interface object
  type(CNN_CS),                  intent(in)  :: CNN    !< Control structure for CNN
  character(len=*),              intent(in)  :: python_bridge_lib !< The library used for language bridging
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in) :: h       !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(inout) :: diffu  !< Zonal acceleration due to convergence of
                                                         !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(inout) :: diffv  !< Meridional acceleration due to convergence
                                                         !! of along-coordinate stress tensor [L T-2 ~> m s-2].
  ! Local Variables
  real, dimension(SZIW_(CNN),SZJW_(CNN),SZK_(GV)) &
                                    :: WH_u     ! The zonal velocity with a wide halo [L T-1 ~> m s-1].
  real, dimension(SZIW_(CNN),SZJW_(CNN),SZK_(GV)) &
                                    :: WH_v     ! The meridional velocity with a wide halo [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: Sx,  &     ! CNN output Sx
                                               Sy,  &     ! CNN output Sy
                                               Sxmean,  & ! CNN output Sxmean
                                               Symean,  & ! CNN output Symean
                                               Sxstd,   & ! CNN output Sxstd
                                               Systd      ! CNN output Systd
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)):: fx     ! CNN output Sx at cell faces
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)):: fy     ! CNN output Sy at cell faces
  ! real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)):: fxmean     ! CNN output Sxmean at cell faces
  ! real, dimension(SZI_(G),SZJB_(G),SZK_(GV)):: fymean     ! CNN output Symean at cell faces
  ! real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)):: fxstd     ! CNN output Sxstd at cell faces
  ! real, dimension(SZI_(G),SZJB_(G),SZK_(GV)):: fystd     ! CNN output Systd at cell faces
  real, dimension(2,SZIW_(CNN),SZJW_(CNN),SZK_(GV)) :: WH_uv     ! CNN input
  real, dimension(6,SZI_(G),SZJ_(G),SZK_(GV)) :: Sxy! CNN output
  type(group_pass_type) :: pass_CNN

  integer :: i, j, k
  integer :: is, ie, js, je, nz, nztemp
  integer :: isdw, iedw, jsdw, jedw

  call cpu_clock_begin(CNN%id_cnn_pre)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isdw = CNN%isdw; iedw = CNN%iedw; jsdw = CNN%jsdw; jedw = CNN%jedw
  if (CNN%CNN_BT) then
    nztemp = 1
  else
    nztemp = nz
  endif

  WH_u = 0.0; WH_v = 0.0;
  do k=1,nztemp
    do j=js,je ; do i=is,ie
      WH_u(i,j,k) = 0.5*( u(I-1,j,k) + u(I,j,k) ) ! Copy the computational section from u into cell center
    enddo ; enddo
    do j=js,je ; do i=is,ie
      WH_v(i,j,k) = 0.5*( v(i,J-1,k) + v(i,J,k) ) ! Copy the computational section from v into cell center
    enddo ; enddo
  enddo

  ! Update the wide halos of WH_u WH_v
  call pass_var(WH_u, CNN%CNN_Domain)
  call pass_var(WH_v, CNN%CNN_Domain)

  ! Combine arrays for CNN input
  WH_uv = 0.0
  do k=1,nztemp
    do j=jsdw,jedw ; do i=isdw,iedw 
      WH_uv(1,i,j,k) = WH_u(i,j,k)
      WH_uv(2,i,j,k) = WH_v(i,j,k)
    enddo ; enddo 
  enddo

  ! Run Python script for CNN inference
  Sxy = 0.0
  call cpu_clock_end(CNN%id_cnn_pre)
  call cpu_clock_begin(CNN%id_cnn_inference)
  select case (lowercase(python_bridge_lib))
  case("forpy")
 !   call forpy_run_python(WH_uv, Sxy, FP_CS, CNN%CNN_BT, G)
  case("smartsim")
    call smartsim_run_python(WH_uv, Sxy, SS_CS, CNN%CNN_BT, CNN%CNN_halo_size)
  end select
  call cpu_clock_end(CNN%id_cnn_inference)

  !Extract data from CNN output
  call cpu_clock_begin(CNN%id_cnn_post)
  call cpu_clock_begin(CNN%id_cnn_post1)
  Sx=0.0; Sy=0.0; Sxmean=0.0; Symean=0.0; Sxstd=0.0; Systd=0.0;
  do k=1,nztemp
    do j=js,je ; do i=is,ie 
      Sx(i,j,k) = Sxy(1,i,j,k)
      Sy(i,j,k) = Sxy(2,i,j,k)
      Sxmean(i,j,k) = Sxy(3,i,j,k)
      Symean(i,j,k) = Sxy(4,i,j,k)
      Sxstd(i,j,k) = Sxy(5,i,j,k)
      Systd(i,j,k) = Sxy(6,i,j,k)
    enddo ; enddo 
  enddo
  call cpu_clock_end(CNN%id_cnn_post1)

  ! Update the halos of Sx Sy
  call cpu_clock_begin(CNN%id_cnn_post2)
  ! call pass_var(Sx, G%Domain)
  ! call pass_var(Sy, G%Domain)
  ! call pass_var(Sxmean, G%Domain)
  ! call pass_var(Symean, G%Domain)
  ! call pass_var(Sxstd, G%Domain)
  ! call pass_var(Systd, G%Domain) 
  call create_group_pass(pass_CNN,Sx,G%Domain)
  call create_group_pass(pass_CNN,Sy,G%Domain)
  call create_group_pass(pass_CNN,Sxmean,G%Domain)
  call create_group_pass(pass_CNN,Symean,G%Domain)
  call create_group_pass(pass_CNN,Sxstd,G%Domain)
  call create_group_pass(pass_CNN,Systd,G%Domain)
  call do_group_pass(pass_CNN,G%Domain)
  call cpu_clock_end(CNN%id_cnn_post2)
 
  call cpu_clock_begin(CNN%id_cnn_post3)
  fx = 0.0; fy = 0.0; 
  do k=1,nz
    do j=js,je ; do I=is-1,ie
      if (CNN%CNN_BT) then
        fx(I,j,k) = 0.5*(Sx(i,j,1) + Sx(i+1,j,1))
      else
        fx(I,j,k) = 0.5*(Sx(i,j,k) + Sx(i+1,j,k))
      endif
      diffu(I,j,k) = diffu(I,j,k) + fx(I,j,k) ! Update diffu with Sx
    enddo ; enddo
    do J=js-1,je ; do i=is,ie
      if (CNN%CNN_BT) then
        fy(i,J,k) = 0.5*(Sy(i,j,1) + Sy(i,j+1,1))
      else
        fy(i,J,k) = 0.5*(Sy(i,j,k) + Sy(i,j+1,k))
      endif
      diffv(i,J,k) = diffv(i,J,k) + fy(i,J,k) ! Update diffv with Sy
    enddo ; enddo
  enddo
  call cpu_clock_end(CNN%id_cnn_post3)

  call cpu_clock_begin(CNN%id_cnn_post4)
  if (CNN%id_CNNu>0)   call post_data(CNN%id_CNNu, fx, CNN%diag)
  if (CNN%id_CNNv>0)   call post_data(CNN%id_CNNv, fy, CNN%diag)
  if (CNN%id_Sxmean>0)   call post_data(CNN%id_Sxmean, Sxmean, CNN%diag)
  if (CNN%id_Symean>0)   call post_data(CNN%id_Symean, Symean, CNN%diag)
  if (CNN%id_Sxstd>0)   call post_data(CNN%id_Sxstd, Sxstd, CNN%diag)
  if (CNN%id_Systd>0)   call post_data(CNN%id_Systd, Systd, CNN%diag)
  if (CNN%id_KE_CNN>0) call compute_energy_source(u, v, h, fx, fy, G, GV, CNN)
  call cpu_clock_end(CNN%id_cnn_post4)

  call cpu_clock_end(CNN%id_cnn_post)

end subroutine CNN_inference

! This is copy-paste from MOM_diagnostics.F90, specifically 1125 line
subroutine compute_energy_source(u, v, h, fx, fy, G, GV, CS)
  type(ocean_grid_type),         intent(in)  :: G      !< The ocean's grid structure.
  type(verticalGrid_type),       intent(in)  :: GV     !< The ocean's vertical grid structure.
  type(CNN_CS),                  intent(in)  :: CS    !< Control structure for CNN

  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in)  :: u      !< The zonal velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in)  :: v      !< The meridional velocity [L T-1 ~> m s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in) :: h       !< Layer thicknesses [H ~> m or kg m-2].
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                                 intent(in) :: fx     !< Zonal acceleration due to convergence of
                                                      !! along-coordinate stress tensor [L T-2 ~> m s-2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                                 intent(in) :: fy     !< Meridional acceleration due to convergence
                                                      !! of along-coordinate stress tensor [L T-2 ~> m s-2]
  
  real :: KE_term(SZI_(G),SZJ_(G),SZK_(GV)) ! A term in the kinetic energy budget
                                 ! [H L2 T-3 ~> m3 s-3 or W m-2]
  real :: tmp(SZI_(G),SZJ_(G),SZK_(GV)) ! temporary array for integration
  real :: KE_u(SZIB_(G),SZJ_(G)) ! The area integral of a KE term in a layer at u-points
                                 ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]
  real :: KE_v(SZI_(G),SZJB_(G)) ! The area integral of a KE term in a layer at v-points
                                 ! [H L4 T-3 ~> m5 s-3 or kg m2 s-3]

  type(group_pass_type) :: pass_KE_uv !< A handle used for group halo passes

  integer :: is, ie, js, je, Isq, Ieq, Jsq, Jeq, nz
  integer :: i, j, k

  real :: uh !< Transport through zonal faces = u*h*dy,
             !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: vh !< Transport through meridional faces = v*h*dx,
             !! [H L2 T-1 ~> m3 s-1 or kg s-1].
  real :: global_integral !< Global integral of the energy effect of CNN model [W]

  is  = G%isc  ; ie  = G%iec  ; js  = G%jsc  ; je  = G%jec ; nz = GV%ke
  Isq = G%IscB ; Ieq = G%IecB ; Jsq = G%JscB ; Jeq = G%JecB
  
  call create_group_pass(pass_KE_uv, KE_u, KE_v, G%Domain, To_North+To_East)
  
  KE_term(:,:,:) = 0.
  tmp(:,:,:) = 0.
  ! Calculate the KE source from Zanna-Bolton2020 [H L2 T-3 ~> m3 s-3].
  do k=1,nz
    KE_u(:,:) = 0.
    KE_v(:,:) = 0.
    do j=js,je ; do I=Isq,Ieq
      uh = u(I,j,k) * 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i+1,j)*h(i+1,j,k)) * &
        G%dyCu(I,j)
      KE_u(I,j) = uh * G%dxCu(I,j) * fx(I,j,k)
    enddo ; enddo
    do J=Jsq,Jeq ; do i=is,ie
      vh = v(i,J,k) * 0.5 * (G%mask2dT(i,j)*h(i,j,k) + G%mask2dT(i,j+1)*h(i,j+1,k)) * &
        G%dxCv(i,J)
      KE_v(i,J) = vh * G%dyCv(i,J) * fy(i,J,k)
    enddo ; enddo
    call do_group_pass(pass_KE_uv, G%domain)
    do j=js,je ; do i=is,ie
      KE_term(i,j,k) = 0.5 * G%IareaT(i,j) &
          * (KE_u(I,j) + KE_u(I-1,j) + KE_v(i,J) + KE_v(i,J-1))
      ! copy-paste from MOM_spatial_means.F90, line 42
      tmp(i,j,k) = KE_term(i,j,k) * G%areaT(i,j) * G%mask2dT(i,j)
    enddo ; enddo
  enddo

  global_integral = reproducing_sum(tmp)

  call post_data(CS%id_KE_CNN, KE_term, CS%diag)

end subroutine compute_energy_source
  
end module MOM_CNN_GZ21
