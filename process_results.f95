!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module can process solution from AGMG method (concentration at each voxel)
! and computes effective diffusivity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module process_results
use types, only: dp,physical_prop,model_data
use common_func, only: ijk_to_l, wrap_ind, get_concentration, get_diffusivity
implicit none
private
public compute_effective_diffusivity
contains
subroutine compute_effective_diffusivity(conc,foam_struc,foam_dim,model,physical_properties,deff,vrc)
	real(dp), allocatable, intent(in) :: conc(:)
	integer, allocatable, intent(in) :: foam_struc(:,:,:)
	integer, intent(in) :: foam_dim
	type(model_data), intent(in) :: model
	type(physical_prop), intent(in) :: physical_properties
	real(dp), intent(out) :: deff,vrc
	real(dp), allocatable :: Jiaveg(:)
	integer i,j,k,l1,l2,phase1,phase2
	real(dp) :: D12,Jiijk,c1,c2,deffi
	allocate(Jiaveg(foam_dim-1))
	Jiaveg=0
	deff=0 !average
	vrc=0 !variance
	do i=1,foam_dim-1 ! foam_dim-1
	do j=1,foam_dim
	do k=1,foam_dim
		l1=ijk_to_l(i,j,k,foam_dim)
		l2=ijk_to_l(i+1,j,k,foam_dim)
		phase1=foam_struc(k,j,i)
		phase2=foam_struc(k,j,i+1)
		D12=get_diffusivity(phase1,phase2,physical_properties)
		Jiijk=-D12/model%h*(conc(l2)-conc(l1))
		Jiaveg(i)=Jiaveg(i)+Jiijk
	end do
	end do
	Jiaveg(i)=Jiaveg(i)/foam_dim**2
	deffi=-Jiaveg(i)*(foam_dim-1)*model%h/(model%c1-model%c0)
	deff=deff+deffi
	vrc=vrc+deffi**2
	end do
	!compute vector J of molar flux intensity in i,j,k
	deff=deff/(foam_dim-1)
	vrc=sqrt(abs(vrc/(foam_dim-1)-deff**2))/deff*100
	!average vector Ji = sum at layer
	!deff =1
	
end subroutine

end module
