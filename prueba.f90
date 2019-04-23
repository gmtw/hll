program pruebas
	implicit none
	integer :: i 
	integer, parameter ::nx=2
	real :: u(1,1:nx+1)

	do i=0,nx
		u(1,i+1)=i
	end do
	print*, u(1,0)

    
end program pruebas