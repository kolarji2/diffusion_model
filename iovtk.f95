!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module handle loading and saving of foam structure 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module iovtk
use types, only: dp
use common_func, only: wrap_ind,ijk_to_l
implicit none
private
public load_vtk_txt, write_vtk_txt,write_vtkdata_to_vtk,rotate_box,write_solution_to_file
contains

subroutine load_vtk_txt(file_name,vtk_data,data_dim)
	!note that vtk_data is returned in array as vtk_data(k,j,i)
	!because of faster acces to slices vtk_data(:,:,i)
	CHARACTER(len=128), intent(in) :: file_name
	integer, allocatable, intent(out) :: vtk_data(:,:,:)
	integer, intent(out) :: data_dim
	integer :: fid
	integer :: i,j,k
	integer :: voxel_val,max_dim
	character(len=1000) :: s
	open(newunit=fid, file=file_name, status="old")
	max_dim=1000
	!read header
	read(fid, *) data_dim
	if (allocated(vtk_data)) deallocate(vtk_data)
	if (data_dim<max_dim) then
		allocate(vtk_data(data_dim,data_dim,data_dim))
	else
		return
	end if
	read(fid, *)
	do i=1,data_dim
		do j=1,data_dim
			read(fid, *) s
			do k=1,data_dim
				if (s(k:k)=='1') then
					vtk_data(k,j,i)=1
				else
					vtk_data(k,j,i)=0
				end if
			end do
		end do
	read(fid, *)
	end do	
	close(fid)
end subroutine

subroutine write_vtk_txt(file_name,vtk_data,data_dim)
	CHARACTER(len=*), intent(in) :: file_name
	integer, allocatable, intent(in) :: vtk_data(:,:,:)
	integer, intent(in) :: data_dim
	integer :: fid
	integer :: i,j,k
	integer :: voxel_val
	open(newunit=fid, file=file_name, status="replace")
	!write header
	write(fid,"(I3)") data_dim
	write(fid,"(A1)") '#'
	do i=1,data_dim
		do j=1,data_dim
			write(fid,"(1000I1)") vtk_data(:,j,i)
		end do
	write(fid,"(A1)") '#'
	end do	
	close(fid)
end subroutine

subroutine write_solution_to_file(file_name,conc,data_dim,c0,c1)
	CHARACTER(len=*), intent(in) :: file_name
	real(dp), allocatable, intent(in) :: conc(:)
	integer, intent(in) :: data_dim
	real(dp), intent(in) :: c0,c1
	integer :: fid
	integer :: i,l
	integer :: voxel_val
	open(newunit=fid, file=file_name, status="replace")
	!write header
	write(fid,"(A26)") "# vtk DataFile Version 4.1"
	write(fid,"(A10)") "vtk output"
	write(fid,"(A5)") "ASCII"
	write(fid,"(A25)") "DATASET STRUCTURED_POINTS"
	write(fid,"(A11,3(I3A1))") "DIMENSIONS ", data_dim," ",data_dim," ",data_dim
	write(fid,"(A13)") "SPACING 1 1 1"
	write(fid,"(A12)") "ORIGIN 0 0 0"
	write(fid,"(A11,I20)") "POINT_DATA ", data_dim**3
	write(fid,"(A26)") "COLOR_SCALARS voxel_data 1"
	i=0
	do l=1,(data_dim**3-1),2
		if (l==1 .and. modulo(data_dim**3,2)/=0) then
			write(fid,"(F10.8A1F10.8A1F10.8)") gc(conc(l),c0,c1)," ", &
						gc(conc(l+1),c0,c1)," ",gc(conc(l+2),c0,c1)
			i=1
		else
			write(fid,"(F10.8A1F10.8)") gc(conc(l+1),c0,c1)," ", gc(conc(l+i+1),c0,c1)
		end if
	end do
	
	close(fid)
end subroutine

subroutine write_vtkdata_to_vtk(file_name,vtk_data,data_dim)
	CHARACTER(len=*), intent(in) :: file_name
	integer, allocatable, intent(in) :: vtk_data(:,:,:)
	integer, intent(in) :: data_dim
	integer :: i,j,k,l
	real(dp), allocatable :: conc(:)
	allocate(conc(data_dim**3))
	do i=1,data_dim
		do j=1,data_dim
			do k=1,data_dim
				l=ijk_to_l(i,j,k,data_dim)
				conc(l)=1.0*vtk_data(k,j,i)
			end do
		end do
	end do	
	call write_solution_to_file(file_name,conc,data_dim,0.0_dp,1.0_dp)
end subroutine


real(dp) function gc(conc,c0,c1) result(res)
		real(dp), intent(in) ::  conc,c0,c1
		res=(conc-c0)/(c1-c0)
	end function

subroutine rotate_box(vtk_data,data_dim,ind)
	integer, allocatable, intent(inout) :: vtk_data(:,:,:)
	integer, intent(in) :: data_dim
	integer, intent(in) :: ind
	integer, allocatable :: rotated_vtk_data(:,:,:)
	integer :: i,j,k
	integer :: ijkdash(3)
	if (ind==1) return
	print *, "Rotating box"
	allocate(rotated_vtk_data(data_dim,data_dim,data_dim))
	do i=1,data_dim
		do j=1,data_dim			
			do k=1,data_dim
					ijkdash(wrap_ind(ind,3))=i
					ijkdash(wrap_ind(ind+1,3))=j
					ijkdash(wrap_ind(ind+2,3))=k			
					rotated_vtk_data(k,j,i)=vtk_data(ijkdash(3),ijkdash(2),ijkdash(1))	
			end do
		end do
	end do
	vtk_data=rotated_vtk_data
end subroutine
end module
