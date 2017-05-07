!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module can load configuration from json file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module load_configuration
use types, only: dp,charlen,physical_prop,model_data,agmg_config
use json_module
use error_handler, only: error
implicit none
private
public load_from_json
contains

subroutine load_from_json(file_name,model,agmg_conf,physical_properties)
	character(len=*), intent(in) :: file_name
	type(model_data), intent(out) :: model
	type(agmg_config), intent(out) :: agmg_conf
	type(physical_prop),allocatable, intent(out) :: physical_properties(:)
	type(json_file) :: json
	type(json_value), pointer :: prop_pointer
	logical :: found
	integer :: propid,maxpropid,i
	integer, allocatable :: ilen_foam_files(:)
	character(len=charlen) :: propertyid,propertyidpath
	character(len=:), allocatable :: foam_files(:),nameid,output_file
	! initialize the class
    call json%initialize()
    !load file
	call json%load_file(filename =file_name)
	!load module info
	call json%get('model.foam-files', foam_files, ilen_foam_files,found)
		if ( .not. found ) call error("model foam-files not found in json file")
	call json%get('model.output-file', output_file, found)
		if ( .not. found ) call error("model output-file not found in json file")
	call json%get('model.mainaxis', model%mainaxis, found)
        if ( .not. found ) call error("model mainaxis not found in json file")
	call json%get('model.c0', model%c0, found)
        if ( .not. found ) call error("model c0 not found in json file")
    call json%get('model.c1', model%c1, found)
        if ( .not. found ) call error("model c1 not found in json file")
    call json%get('model.h', model%h, found)
        if ( .not. found ) call error("model h not found in json file")
    call json%get('model.p', model%p, found)
        if ( .not. found ) call error("model p not found in json file")    
    allocate(model%foam_files(size(foam_files)))
    do i=1,size(foam_files)
		write(*,*) foam_files(i)
		model%foam_files(i)=foam_files(i)
    end do    
	model%output_file=output_file
    
    !load configuration for AGMG solver
    call json%get('agmg-solver.ijob', agmg_conf%ijob, found)
        if ( .not. found ) call error("AGMG ijob not found in json file")
    call json%get('agmg-solver.iprint', agmg_conf%iprint, found)
        if ( .not. found ) call error("AGMG iprint not found in json file")
    call json%get('agmg-solver.nrest', agmg_conf%nrest, found)
        if ( .not. found ) call error("AGMG nrest not found in json file")
    call json%get('agmg-solver.iter', agmg_conf%iter, found)
        if ( .not. found ) call error("AGMG iter not found in json file")
    call json%get('agmg-solver.tol', agmg_conf%tol, found)
        if ( .not. found ) call error("AGMG tol not found in json file")
    !load physical_properties
    maxpropid=0
    !get maxpropid
    do propid=1,100
		write(propertyid,"(A20I1A)") 'physical-properties(',propid,')'
		!property list id i exist
		call json%get(propertyid, prop_pointer, found)
			if ( .not. found ) exit
		maxpropid=maxpropid+1
	end do
	if (maxpropid==0) call error("No physical properties in json file!") 
    allocate(physical_properties(maxpropid))
    do propid=1,maxpropid
		write(propertyid,"(A20I1A)") 'physical-properties(',propid,')'
		!load physical-properties
		!property nameid
			write(propertyidpath,"(A22A7)") propertyid,".nameid"
			call json%get(propertyidpath,nameid,found)
				if ( .not. found ) call error("Physical property nameid not found in json file")
			physical_properties(propid)%nameid=nameid
		!property dgas
			write(propertyidpath,"(A22A5)") propertyid,".dgas"
			call json%get(propertyidpath,physical_properties(propid)%D_gas,found)
				if ( .not. found ) call error("Physical property dgas not found in json file")
		!property dpol
			write(propertyidpath,"(A22A5)") propertyid,".dpol"
			call json%get(propertyidpath,physical_properties(propid)%D_pol,found)
				if ( .not. found ) call error("Physical property dpol not found in json file")
		!property H
			write(propertyidpath,"(A22A2)") propertyid,".H"
			call json%get(propertyidpath,physical_properties(propid)%Henry_const,found)
				if ( .not. found ) call error("Physical property H not found in json file")		
	end do
	call json%destroy()
end subroutine
end module
