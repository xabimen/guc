module distf1_mod

implicit none

public		:: distf1

contains

subroutine distf1(n,n_atoms,coor,coorref)
!	n -> number of different elements
!	n_atoms -> size(0:n)
implicit none
integer, intent(in)								        :: n
integer, dimension(0:),intent(in)				        :: n_atoms
real(kind=8), dimension(:,:), intent(in)		        :: coor
real(kind=8), dimension(:,:), allocatable, intent(out)  :: coorref
integer											        :: i, j
real(kind=8), dimension(3)						        :: ref

allocate(coorref(sum(n_atoms),3))


do i = 1, n
	ref = coor(n_atoms(i-1)+1,:)

	do j = sum(n_atoms(1:i-1)) + 1, sum(n_atoms(1:i))
		coorref(j,:) = modulo(coor(j,:) - ref,1.0d0)
	enddo
enddo

end subroutine

end module
