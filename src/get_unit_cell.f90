program get_unit_cell

use read_data
use distf1_mod
use supercell_mod
use intersection_mod
use distances_mod
use check_angles_mod
use check_cell_mod
use output

implicit none

!Deklarazioa
!******************************************************
integer                                         :: n_elements, i, i1, i2, i3, j
integer, dimension(1)							:: k
integer, dimension(100000,3)						:: loc
integer, dimension(:), allocatable              :: n_atoms, n_atoms8, n_atoms_barne
real(kind=8), dimension(3,3)                    :: cell, newcell
real(kind=8), dimension(:,:), allocatable       :: coor, coorref, coor8, coorref8, inter_vec, coor_barne
real(kind=8), dimension(:), allocatable			:: dist
real(kind=8), dimension(100000)					:: fit
real(kind=8), dimension(3)						:: a, b, c, a2, b2, c2, norm
real(kind=8)                                    :: tol, theta, alpha
character(len=4), dimension(:), allocatable     :: names
character(len=20)                               :: file_name
logical											:: angles, unit_cell
real(kind=8), parameter						:: pi = acos(-1.0d0)
character(len=10), dimension(2)				:: arg
!******************************************************


call get_command_argument(1,file_name)
call get_command_argument(2,arg(1))
call get_command_argument(3,arg(2))
read(arg(1),fmt=*) n_elements
read(arg(2), fmt=*) tol
print*, "number of elemets", n_elements
print*, "tolerance", tol
allocate(n_atoms(0:n_elements),n_atoms8(0:n_elements), names(n_elements))
!n_atoms bektorea 0-tik definituta dago holan gero iteratzeko
!orduan errezago izan dadin.

call read_data_file(file_name,n_elements,cell,names,n_atoms,coor)

!do i = 1, sum(n_atoms)
!    print*, coor(i,:)
!enddo

!print*,

call distf1(n_elements,n_atoms,coor,coorref)

!do i = 1, sum(n_atoms)
!    print*, coorref(i,:)
!enddo

call supercell(n_atoms,coor,n_atoms8,coor8)
call supercell(n_atoms,coorref,n_atoms8,coorref8)
!print*, "COORREF8"
!do i = 1, sum(n_atoms8)
!    print*, i, coorref8(i,:)
!enddo

call intersection(n_atoms8,coorref8,tol,inter_vec)

call distances(inter_vec,cell,dist)

!print*, "INTER_VEC"
!do i = 1, size(inter_vec,1)
!	print*, i, inter_vec(i,:), dist(i)
!enddo

!GET POSSIBLE CELLS AND THEIR FITNESS
fit = 10000.0
j = 1
do i1 = 2, size(inter_vec,1)
	do i2 = 2, size(inter_vec,1)
		do i3 = 2, size(inter_vec,1)
			a = inter_vec(i1,:)
			b = inter_vec(i2,:)
			c = inter_vec(i3,:)
			if (dist(i1) < maxval(norm2(cell,1))/2.0d0) then
				call check_angles(cell,a,b,c,angles)
				if (angles) then 
					loc(j,:) = (/i1,i2,i3/)
					fit(j) = norm2(a)*norm2(b)*norm2(c)
					j = j + 1
				endif
			endif

		enddo
	enddo
enddo


do i = 1, 100
	k = minloc(fit)
	i1 = loc(k(1),1)
	i2 = loc(k(1),2)
	i3 = loc(k(1),3)
	fit(k(1)) = 10000.0

	newcell(1,:) = matmul(inter_vec(i1,:),cell)
	newcell(2,:) = matmul(inter_vec(i2,:),cell)
	newcell(3,:) = matmul(inter_vec(i3,:),cell)
	a2 = newcell(1,:)
	b2 = newcell(2,:)
	c2 = newcell(3,:)
	!print*, i1, i2, i3
	print*, "Possible unit cell"
	print*, "a", inter_vec(i1,:),"norm-> ", norm2(newcell(1,:))
	print*, "b", inter_vec(i2,:), "norm-> ", norm2(newcell(2,:))
	print*, "c", inter_vec(i3,:), "norm-> ", norm2(newcell(3,:))

	!print*, "theta", theta
	!norm = cross_product(a2,b2)
	!norm = norm/norm2(norm)

	!alpha = asin(abs(dot_product(c2,norm))/(norm2(c2)*norm2(norm)))*180.0d0/pi 
	!print*, "alpha", alpha

	call check_cell(newcell,cell,coor,n_atoms,unit_cell,coor_barne,n_atoms_barne)

	if (unit_cell) then
		print*, "Translational symmetry found"
		print*, "Output in file 'output.vasp'"
		call write_output_vasp(newcell,coor_barne,n_atoms_barne,names)
		exit
	else
		print*, "Translational symmetry not found"
	endif

	print*, "************************************"
	if (i==100) then
		print*,"ERROR: Maximum iterations reached"
	endif
enddo

end program
