!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module generate matrix representing implicit euler method for 
! calculating steady state diffusion in foam with fixed boundaries 
! on index i and periodic boundaries on indices j,k. 
!	Ac=b
! Resulting sparse matrix is coded in the form a,ja,ia,b suitable for
! AGMG method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fvm_iem_matrix
use types, only: dp,physical_prop
use common_func, only: ijk_to_l, wrap_ind, get_concentration, get_diffusivity
implicit none
private
public generate_matrix
contains

subroutine generate_matrix(c0,c1,h,foam_struc,foam_dim,physical_properties,a,ja,ia,b,x)
	real(dp), intent(in) :: c0,c1,h
	integer, allocatable, intent(in) :: foam_struc(:,:,:)
	integer, intent(in) :: foam_dim
	type(physical_prop), intent(in) :: physical_properties
	real(dp), allocatable, intent(out) :: a(:),b(:),x(:)
	integer, allocatable, intent(out) :: ja(:)
	integer, allocatable, intent(out) :: ia(:)
	integer :: ijkdash(3)
	integer :: phase1,phase2
	integer i,j,k,ntot_voxels,l,m,ii,dii
	integer ntot_aelem,n_aelem_fix,n_aelem_cont,curr_aelem,curr_row
	real(dp) :: D_gas,D_pol,Henry_const,Dsum,D12,c_linaprox
	D_gas=physical_properties%D_gas
	D_pol=physical_properties%D_pol
	Henry_const=physical_properties%Henry_const
	ntot_voxels=foam_dim**3
	! 2*foam_dim**2 matrix elemets for fixed conditions
	! 7*(foam_dim-2)*foam_dim**2 matrix elemets for other voxels
	n_aelem_fix=2*foam_dim**2
	n_aelem_cont=7*(foam_dim-2)*foam_dim**2
	ntot_aelem=n_aelem_fix+n_aelem_cont
	if (.not. allocated(a)) allocate(a(ntot_aelem)) ! allocate number of sparse matrix elements
	if (.not. allocated(ja)) allocate(ja(ntot_aelem)) ! allocate number of sparse matrix elements
	if (.not. allocated(ia)) allocate(ia(ntot_voxels+1)) !allocate number of sparse matrix row+1
	if (.not. allocated(b)) allocate(b(ntot_voxels)) !allocate number of sparse matrix row
	if (.not. allocated(x)) allocate(x(ntot_voxels))
	!assign fixed boundary to index i=0 i=foam_dim
	curr_aelem=1
	curr_row=1
	ia(curr_row)=curr_aelem
	do i=1,foam_dim
		c_linaprox=c0+(i-1)*(c1-c0)/(foam_dim-1)
	do j=1,foam_dim
	do k=1,foam_dim
		l=ijk_to_l(i,j,k,foam_dim) ! mapping to vector
		phase1=foam_struc(k,j,i) !structure saved as slices [:,j,i]
		x(l)=c_linaprox !linear guess
		if (i==1 .OR. i==foam_dim) then
			!assign fixed boundary condition to voxel (i,j,k) and elem b_l
			a(curr_aelem)=1.0_dp
			ja(curr_aelem)=l
			if (i==1) then
				b(l)=c0
			else
				b(l)=c1
			end if
			curr_aelem=curr_aelem+1
		else
			!assign balance condition to voxel (i,j,k) and elem b_l
			b(l)=0		
			Dsum=0
			!print *, "-------------------------------------------------------------"
			do ii=1,3
				do dii=-1,1,2
					ijkdash=[i,j,k]
					ijkdash(ii)=wrap_ind(ijkdash(ii)+dii,foam_dim) !neighbouring voxel
					m=ijk_to_l(ijkdash(1),ijkdash(2),ijkdash(3),foam_dim)
					phase2=foam_struc(ijkdash(3),ijkdash(2),ijkdash(1))
					D12=get_diffusivity(phase1,phase2,physical_properties)
					Dsum=Dsum+D12
					a(curr_aelem)=-D12/h**2		
					!print *, m-l,	D12/h**2	
					ja(curr_aelem)=m
					curr_aelem=curr_aelem+1
				end do			
			end do
			!for voxel i,j,k
			a(curr_aelem)=Dsum/h**2
			!print *, l, -Dsum
			ja(curr_aelem)=l
			curr_aelem=curr_aelem+1
		end if
		curr_row=curr_row+1
		ia(curr_row)=curr_aelem
	end do
	end do
	end do
end subroutine

	
end module
