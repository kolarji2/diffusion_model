!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program run parametric study testing various resolution
! and cell-wall ratio on space-wall and cubic grid.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_simple
	use types, only: dp,charlen,physical_prop,model_data,agmg_config
	use error_handler, only: error
	use common_func, only: wrap_ind,ijk_to_l
	use structure_generator, only: gen_walls, gen_cubic
	use iovtk, only: load_vtk_txt, write_vtk_txt, rotate_box
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
	integer :: dwall_arr(5), dspace_arr(5)
	integer :: n,iter,iprint,ijob,nrest
	integer :: m,i,j,k,dspace,dwall,nrepeat
	logical :: dowall
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
	!	load foam structure and rotate for diffusion along set mainaxis 
	
	!file_name_out="teststruc.txt"
	!call write_vtk_txt(file_name_out,vtk_data,data_dim)
	!call load_vtk_txt(model%foam_file,vtk_data,data_dim)
	!call rotate_box(vtk_data,data_dim,model%mainaxis)
	!generate matrix
	c0=model%c0
	c1=model%c1
	h=model%h
	!dwall_arr=[1,2,3,4,5]
	!dspace_arr=[5,10,20,50,95]
	dwall_arr=[1,2,4,8,16]
	dspace_arr=[5,10,20,40,80]
	
	open(newunit=outstream,file=model%output_file, status="replace")
	dowall=.true.
	do m=1,2
	do j=1,size(dspace_arr)
	do k=1,size(dwall_arr)
		dspace=dspace_arr(j)
		dwall=dwall_arr(k)
		nrepeat=min(6*(dspace+dwall),233)/(dspace+dwall)
		if (dowall) then			
			call gen_walls(dspace,dwall,nrepeat,vtk_data,data_dim)
			write(outstream,*) "GEN WALLS dspace: ",dspace," dwall: ",dwall," res: ",data_dim,"n:",nrepeat
			write(*,*) "GEN WALLS dspace: ", dspace, " dwall: ", dwall, " res: ", data_dim,"n:",nrepeat
		else
			call gen_cubic(dspace,dwall,nrepeat,vtk_data,data_dim)
			write(outstream,*) "GEN CUBIC dspace: ",dspace," dwall: ",dwall," res: ", data_dim,"n:",nrepeat
			write(*,*) "GEN CUBIC dspace: ", dspace, " dwall: ", dwall, " res: ", data_dim,"n:",nrepeat
		end if
		write(outstream,"(A18,A2,A16,A2,A6)") "D_eff",", ","rel vrc[%]",", ","nameid"
		!write(*,"(A18,A2,A16,A2,A6)") "D_eff",", ","rel vrc[%]",", ","nameid"
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
			
			!	process results
			!	compute effective difffusivity
			call compute_effective_diffusivity(x,vtk_data,data_dim,model,physical_properties(i),deff,vrcdeff)
			!
			!rite(*,"(ES18.9,A2,F16.10,A2,A16)") deff,", ",vrcdeff,", ",physical_properties(i)%nameid
			write(outstream,"(ES18.9,A2,F16.10,A2,A20)") deff,", ",vrcdeff,", ",physical_properties(i)%nameid
			!
			!      uncomment the following lines to write solution on disk for checking
			!
			!       open(outstream,file='sol.out',form='formatted')
			!       write(10,'(e22.15)') x(1:n)
			!       close(10)
			!process results
			
		end do
		call flush(outstream)
	end do
	end do
	dowall=.false.
	end do
	close(outstream)
end 
