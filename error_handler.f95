!author: Jiri Kolar
!email: jiri1kolar@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module is responsible for error handling and reporting info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module error_handler
use types, only: dp, longcharlen
implicit none
private
public error, warning
contains

subroutine error(msg)
	character(len=*), intent(in) :: msg 
	write(*,*) "Error: ",msg
	call Exit(0)
end subroutine

subroutine warning(msg)
	character(len=*), intent(in) :: msg 
	write(*,*) "Warning: ",msg
end subroutine

end module
