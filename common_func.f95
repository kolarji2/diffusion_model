!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains common functions for the program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module common_func
use types, only: dp,physical_prop
implicit none
private
public ijk_to_l, wrap_ind, get_concentration, get_diffusivity
contains
	integer function ijk_to_l(i,j,k,arr_dim) result(l)
	!get position in row vector from 3D array indices i,j,k
	!prefered ordering:
	! (i,j,k) = (1,1,1), (1,1,2), ... (1,1,k_max) ... (1,j_max,k_max) ..
	! .. (i_max,j_max,k_max)
	integer, intent(in) :: i,j,k,arr_dim
	l=(i-1)*arr_dim**2+(j-1)*arr_dim+k	
	end function
	
	integer function wrap_ind(ind,arr_dim) result(wind)
		integer,intent(in) :: ind,arr_dim
		wind=modulo((ind-1),arr_dim)+1
	end function
	
	real(dp) function get_concentration(conc,phase,physical_properties) result(phaseconc)
		real(dp), intent(in) :: conc
		integer, intent(in) ::  phase
		type(physical_prop), intent(in) :: physical_properties
		!concentration is set to polymer phase
		if (phase==1) then
			phaseconc=conc/physical_properties%Henry_const
		else
			phaseconc=conc
		end if
	end function
	
	real(dp) function get_diffusivity(phase1,phase2,physical_properties) result(D)
		integer, intent(in) ::  phase1,phase2
		type(physical_prop), intent(in) :: physical_properties
		real(dp) :: D1,D2
			D1=get_diffusivity_(phase1,physical_properties)
			D2=get_diffusivity_(phase2,physical_properties)
			D=2.0*D1/(D1+D2)*D2		
	end function
	
	real(dp) function get_diffusivity_(phase,physical_properties) result(D)
		integer, intent(in) ::  phase
		type(physical_prop), intent(in) :: physical_properties
		if (phase==1) then
			D=physical_properties%D_pol/physical_properties%Henry_const
		else
			D=physical_properties%D_gas
		end if
	end function	
end module
