module check_angles_mod

use matrix_operations

implicit none

public	:: check_angles


contains

subroutine check_angles(cell,a,b,c,angles)
implicit none
real(kind=8), dimension(:,:), intent(in)		:: cell
real(kind=8), dimension(:), intent(in)			:: a, b, c
logical, intent(out)							:: angles
real(kind=8), dimension(3)						:: a2, b2, c2, norm
real(kind=8)									:: theta, alpha
real(kind=8), parameter							:: pi = acos(-1.0d0)

a2 = matmul(a,cell)
b2 = matmul(b,cell)
c2 = matmul(c,cell)

angles = .false.

if (norm2(a2) > 0.1d0 .and. norm2(b2) > 0.1d0 .and. norm2(c2) > 0.1d0) then
	theta = acos(dot_product(a2,b2)/(norm2(a2)*norm2(b2)))*180.0d0/pi
	!print*, "theta", theta
	norm = cross_product(a2,b2)
	norm = norm/norm2(norm)

	alpha = asin(abs(dot_product(c2,norm))/(norm2(c2)*norm2(norm)))*180.0d0/pi 
	!print*, "alpha", alpha
	if (theta >= 60.0d0 .and. theta <= 120.0d0 .and.&
		alpha >= 60.0d0 .and. alpha <= 120.0d0) then

		angles = .true.
	endif

endif



end subroutine


end module
