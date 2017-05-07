!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module can generate simple structures useful for model validation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module structure_generator
use types, only: dp,physical_prop,model_data
use common_func, only: ijk_to_l, wrap_ind, get_concentration, get_diffusivity
implicit none
private
public gen_walls,gen_cubic
contains
subroutine gen_walls(dspace,dwall,nrepeat,vtk_data,data_dim)
	integer, intent(in) :: dspace,dwall, nrepeat
	integer, allocatable, intent(out) :: vtk_data(:,:,:)
	integer, intent(out) :: data_dim
	integer :: i,j
	data_dim=(dspace+dwall)*nrepeat
	if (allocated(vtk_data)) deallocate(vtk_data)
	allocate(vtk_data(data_dim,data_dim,data_dim))
	do i=1,nrepeat
		do j=1,dspace
			vtk_data(:,:,(i-1)*(dspace+dwall)+j)=0
		end do
		
		do j=1,dwall
			vtk_data(:,:,(i-1)*(dspace+dwall)+dspace+j)=1
		end do		
	end do	
end subroutine

subroutine gen_cubic(dspace,dwall,nrepeat,vtk_data,data_dim)
	integer, intent(in) :: dspace, dwall, nrepeat
	integer, allocatable, intent(out) :: vtk_data(:,:,:)
	integer, intent(out) :: data_dim
	integer :: i,j,k
	data_dim=(dspace+dwall)*nrepeat
	if (allocated(vtk_data)) deallocate(vtk_data)
	allocate(vtk_data(data_dim,data_dim,data_dim))
	vtk_data=0
	do i=1,nrepeat
		do j=1,dwall
			vtk_data(:,:,(i-1)*(dspace+dwall)+dspace+j)=1
		end do
		do j=1,dwall
			vtk_data(:,(i-1)*(dspace+dwall)+dspace+j,:)=1
		end do	
		do j=1,dwall
			vtk_data((i-1)*(dspace+dwall)+dspace+j,:,:)=1
		end do	
	end do	
	
end subroutine


end module
