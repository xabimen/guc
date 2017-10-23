module intersection_mod

implicit none

public  :: intersection

contains

subroutine intersection(n_atoms8,coor8,tol,inter_vec)
integer, dimension(0:), intent(in)						:: n_atoms8
real(kind=8), dimension(:,:), intent(in)				:: coor8
real(kind=8), intent(in)								:: tol
real(kind=8), dimension(:,:), allocatable, intent(out)	:: inter_vec
real(kind=8), dimension(size(coor8,1),3)					:: aux
real(kind=8)											:: val
integer													:: i, j, k, kont
!logical, dimension(2:size(n_atoms8)-1)					:: found
logical													:: found

aux = 999.0
kont = 1

do i = 1, n_atoms8(1)

	found = .false.
	do j = 2, size(n_atoms8) - 1

		do k = sum(n_atoms8(1:j-1)) + 1, sum(n_atoms8(1:j))
			if (all(coor8(k,:) > coor8(i,:)-tol .and. &
				coor8(k,:) < coor8(i,:)+tol)) then
				!found(j) = .true.
				found = .true.
				exit
			endif

		enddo

	enddo

	if (found) then
		aux(kont,:) = coor8(i,:)

		!Do honekin 0.99=0, 1.99=1 eta negatiboak ere
		do j = 1, 3
			val = abs(aux(kont,j)) - nint(abs(aux(kont,j)))
			if (val > -tol .and. val < 0.0d0) then
				aux(kont,j) = aux(kont,j)! - abs(aux(kont,j))/aux(kont,j)
			endif
		enddo

		kont = kont + 1
	endif
enddo

kont = kont - 1


allocate(inter_vec(kont,3))

do i = 1, kont
	inter_vec(i,:) = aux(i,:)
enddo


end subroutine

end module
