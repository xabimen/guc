module check_cell_mod

use matrix_operations
use supercell_mod

implicit none

public	:: check_cell

contains

subroutine check_cell(newcell,cell,coor,n_atoms,unit_cell,coor_barne,n_atoms_barne)
real(kind=8), dimension(:,:), intent(in)			:: cell, coor, newcell
real(kind=8)										:: tol
integer, dimension(0:), intent(in)					:: n_atoms
logical, intent(out)								:: unit_cell
real(kind=8), dimension(:,:), allocatable, intent(out)	:: coor_barne
integer, dimension(:), allocatable, intent(out)	:: n_atoms_barne
real(kind=8), dimension(:,:), allocatable			:: coor8
integer, dimension(:), allocatable					:: n_atoms8
real(kind=8), dimension(3,3)						:: newcell_inv
real(kind=8), dimension(sum(n_atoms),3)				:: coor2
integer												:: i, j, k
real(kind=8), dimension(3)							:: r_frac
logical												:: errepikatuta
allocate(n_atoms_barne(0:size(n_atoms)-1),n_atoms8(0:size(n_atoms)-1))

tol = 0.1d0

n_atoms_barne(0) = 0

!Jarri kordenatuak lehenengoaren mempe eta denak positiboak
do i = 1, sum(n_atoms)
	coor2(i,:) = modulo(coor(i,:) - coor(1,:),1.0d0)
enddo
!Sortu supergelaxka
call supercell(n_atoms,coor2,n_atoms8,coor8)

!Pasatu denak kartestarretara

do i = 1, sum(n_atoms8)
	coor8(i,:) = matmul(coor8(i,:),cell)
enddo

call inverse(newcell,newcell_inv,3)

call barruan(n_atoms8,coor8,newcell_inv,tol,n_atoms_barne,coor_barne)

!print*, "barruan"
!print*, n_atoms_barne
!do i = 1, sum(n_atoms_barne)
!    print*, coor_barne(i,:)
!enddo
unit_cell = .true.
outter: do i = 1, size(n_atoms) - 1
	do j = sum(n_atoms(1:i-1)) + 1 , sum(n_atoms(1:i))
		r_frac = matmul(matmul(coor2(j,:),cell),newcell_inv)
		r_frac = r_frac - nint(r_frac + tol - 0.5d0)

		errepikatuta = .false.
		do k = sum(n_atoms_barne(1:i-1)) + 1, sum(n_atoms_barne(1:i))
			if (all(abs(r_frac - coor_barne(k,:)) < tol))then
				errepikatuta = .true.
			endif
		enddo
		if (errepikatuta .eqv. .false.) then
			unit_cell = .false.
			exit outter
		endif

	enddo
enddo outter

end subroutine

subroutine barruan(n_atoms,coor,cell_inv,tol,n_atoms_barne,coor_barne)
integer, dimension(0:), intent(in)						:: n_atoms
real(kind=8), dimension(:,:), intent(in)				:: coor, cell_inv
real(kind=8), intent(in)								:: tol
integer, dimension(0:), intent(out)						:: n_atoms_barne
real(kind=8), dimension(:,:), allocatable, intent(out)	:: coor_barne
integer													:: k, i, j, l
real(kind=8), dimension(3)								:: r_frac
logical													:: errepikatuta
real(kind=8), dimension(1000,3)							:: barne_aux

n_atoms_barne = 0
k = 1

do i = 1, size(n_atoms) - 1
	do j = sum(n_atoms(1:i-1)) + 1, sum(n_atoms(1:i))
		
		r_frac = matmul(coor(j,:),cell_inv)

		if (all(r_frac > 0.0d0 - tol) .and. &
			all(r_frac < 1.0d0 + tol)) then
		
			r_frac = r_frac - nint(r_frac + tol  - 0.5d0)

			errepikatuta = .false.
			do l = sum(n_atoms_barne(1:i-1)) + 1 , k-1
				if (all(abs(r_frac - barne_aux(l,:)) < tol)) then
					errepikatuta = .true.
					exit
				endif
			enddo

			if (errepikatuta .eqv. .false.) then
				barne_aux(k,:) = r_frac
				k = k + 1
			endif
		endif

	enddo

	n_atoms_barne(i) = k - 1 - sum(n_atoms_barne(1:i-1))
enddo

allocate(coor_barne(sum(n_atoms_barne),3))

do i = 1, sum(n_atoms_barne)
	coor_barne(i,:) = barne_aux(i,:)
enddo







end subroutine
end module
