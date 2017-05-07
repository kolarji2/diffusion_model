!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module defines used types
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module types
		!constants
		integer, parameter:: dp=kind(0.d0)
		integer, parameter:: charlen=128
		integer, parameter:: longcharlen=256
		!model data
		type physical_prop
			character(len=charlen) :: nameid
			real(dp) :: D_gas ! difusivity in gas
			real(dp) :: D_pol ! difusivity in polymer
			real(dp) :: Henry_const ! Equilibrium constant		
		end type
		type model_data
			character(len=charlen), allocatable :: foam_files(:)
			character(len=charlen) :: output_file
			integer :: mainaxis !index of axis, which parallel to diff grad
			real(dp) :: c0 ! concentration i=1
			real(dp) :: c1 ! concentration i=max_i
			real(dp) :: h ! size of voxel		
			real(dp) :: p ! pressure Pa
		end type
		
		type agmg_config
			integer :: ijob,iprint,nrest,iter
			real(dp) :: tol				
		end type
		
end module types

