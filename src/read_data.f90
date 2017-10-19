module read_data


implicit none


public              :: read_data_file


contains

subroutine read_data_file(file_name,n,cell,names,n_atoms,coor)
!   n -> number of differents elements
!   cell -> cell parameters
!   names -> names of atoms
!   n_atoms -> number of atoms of each type
!   coor -> position of the atoms in fractional coordinates
implicit none
character(len=20), intent(in)                           :: file_name
integer, intent(in)                                     :: n
real(kind=8), dimension(:,:), intent(out)               :: cell
character(len=4), dimension(:),intent(out)              :: names
integer, dimension(0:), intent(out)                      :: n_atoms
real(kind=8), dimension(:,:), allocatable, intent(out)  :: coor
integer                                                 :: i

n_atoms = 0

!Irakurri sarrerako fitxategia
!***************************************************
open(unit=123, file=file_name, status='old', action='read')
read(unit=123,fmt=*)
read(unit=123,fmt=*)
read(unit=123,fmt=*) cell(1,:)
read(unit=123,fmt=*) cell(2,:)
read(unit=123,fmt=*) cell(3,:)
read(unit=123,fmt=*) names
read(unit=123,fmt=*) n_atoms(1:)
read(unit=123,fmt=*)
print*, names
print*, n_atoms(1:)

allocate(coor(sum(n_atoms),3))

do i = 1, sum(n_atoms)
    read(unit=123,fmt=*) coor(i,:)
enddo

close(unit=123)

end subroutine
end module
