module output

implicit none

public  :: write_output_vasp


contains

subroutine  write_output_vasp(cell,coor,n_atoms,names)
implicit none
real(kind=8), dimension(:,:), intent(in)        :: cell, coor
integer, dimension(0:), intent(in)               :: n_atoms
character(len=4), dimension(:), intent(in)      :: names
integer                                         :: i
open(unit=123, file='output.vasp', status='replace', action='write')

write(unit=123,fmt='(a)') "Output"
write(unit=123,fmt='(a)') "1.0"
write(unit=123,fmt='(3f20.10)') cell(1,:)
write(unit=123,fmt='(3f20.10)') cell(2,:)
write(unit=123,fmt='(3f20.10)') cell(3,:)
write(unit=123,fmt='(10a)') names
write(unit=123,fmt='(100i5)') n_atoms(1:)
write(unit=123,fmt='(a)') "Direct"

do i = 1, sum(n_atoms)
    write(unit=123,fmt='(3f20.10)') coor(i,:)
enddo

close(unit=123)

end subroutine


end module
