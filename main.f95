!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program can compute steady state differential equation on 3D foam
! structure. Output is concentration in each voxel and effective diffusivity
! Call other modules from here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
	use types, only: dp,charlen,physical_prop,model_data,agmg_config
	use error_handler, only: error
	use common_func, only: wrap_ind,ijk_to_l
	use structure_generator, only: gen_walls, gen_cubic
	use iovtk, only: load_vtk_txt, write_vtk_txt, rotate_box, write_solution_to_file, write_vtkdata_to_vtk
	use load_configuration, only: load_from_json
	use fvm_iem_matrix, only: generate_matrix
	use process_results, only: compute_effective_diffusivity
	implicit none
	CHARACTER(len=charlen) :: json_file,file_name_out
	integer,allocatable :: vtk_data(:,:,:)
	integer :: data_dim,outstream
	real(dp) :: c0,c1,h
	!coded matrix and vector b for AGMG method
	real(dp), allocatable :: a(:),b(:),x(:)
	integer, allocatable :: ja(:)
	integer, allocatable :: ia(:)
	integer :: n,iter,iprint,ijob,nrest
	integer :: i,ifile, dspace,dwall,nrepeat
    real (dp) :: tol,deff,vrcdeff
	type(physical_prop), allocatable :: physical_properties(:)
	type(model_data) :: model
	type(agmg_config) :: agmg_conf
	!load agruments
	if (iargc()/=1) call error("Specify one argument: input json configuration file!")
    CALL getarg(1, json_file)
	
	!	load configuration
	!	load model data concentration, voxel size
	!	load solver settings
	!	load physical properties for each gas component
	call load_from_json(json_file,model,agmg_conf,physical_properties)
	open(newunit=outstream,file=model%output_file, status="replace")
	
	do ifile=1,size(model%foam_files)
		!	load foam structure and rotate for diffusion along set mainaxis 
		!file_name_out="teststruc.txt"
		!call write_vtk_txt(file_name_out,vtk_data,data_dim)
		call load_vtk_txt(model%foam_files(ifile),vtk_data,data_dim)
		call rotate_box(vtk_data,data_dim,model%mainaxis)
		!generate matrix
		c0=model%c0
		c1=model%c1
		h=model%h
		write(*,*) "File:", model%foam_files(ifile)
		write(outstream,*) "#", model%foam_files(ifile)
		write(outstream,"(A18,A2,A16,A2,A6)") "D_eff",", ","rel vrc[%]",", ","nameid"
		write(*,"(A18,A2,A16,A2,A6)") "D_eff",", ","rel vrc[%]",", ","nameid"
		do i=1,size(physical_properties)
			!	prepare matrix
			call generate_matrix(c0,c1,h,vtk_data,data_dim,physical_properties(i),a,ja,ia,b,x)
			!	agmg solver
			iter=agmg_conf%iter		
			tol=agmg_conf%tol
			iprint=agmg_conf%iprint
			ijob=agmg_conf%ijob
			nrest=agmg_conf%nrest
			n=data_dim**3
			if (ijob/=0 .and. ijob/=10) then
				call dagmg(N,a,ja,ia,b,x,1,iprint,nrest,iter,tol)
			end if
			call dagmg(N,a,ja,ia,b,x,ijob,iprint,nrest,iter,tol)
			call write_solution_to_file("conc_file.vtk",x,data_dim,c0,c1)
			call write_vtkdata_to_vtk("struc_file.vtk",vtk_data,data_dim)
			!	process results
			!	compute effective difffusivity
			call compute_effective_diffusivity(x,vtk_data,data_dim,model,physical_properties(i),deff,vrcdeff)
			!
			write(*,"(ES18.9,A2,F16.10,A2,A16)") deff,", ",vrcdeff,", ",physical_properties(i)%nameid
			write(outstream,"(ES18.9,A2,F16.10,A2,A20)") deff,", ",vrcdeff,", ",physical_properties(i)%nameid
			
			!
			!      uncomment the following lines to write solution on disk for checking
			!
			!       open(outstream,file='sol.out',form='formatted')
			!       write(10,'(e22.15)') x(1:n)
			!       close(10)
			!process results
			
		end do
	end do
	close(outstream)
end 
