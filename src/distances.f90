module distances_mod

implicit none

public	:: distances

contains

subroutine distances(mat,cell,vec)
real(kind=8), dimension(:,:), intent(in)				:: mat, cell
real(kind=8), dimension(:), allocatable, intent(out)	:: vec
integer													:: i
allocate(vec(size(mat,1)))

do i = 1, size(vec)
	vec(i) = norm2(matmul(mat(i,:),cell))
enddo
end subroutine

end module
