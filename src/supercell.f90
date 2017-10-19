module supercell_mod

implicit none

public :: supercell

contains

subroutine supercell(n_atoms,coor,n_atoms8,coor8)
integer, dimension(0:), intent(in)              :: n_atoms
real(kind=8), dimension(:,:), intent(in)        :: coor
integer, dimension(0:)                          :: n_atoms8
real(kind=8), dimension(:,:), allocatable       :: coor8
integer                                         :: i, j, k

allocate(coor8(sum(n_atoms)*8,3))

n_atoms8 = n_atoms*8

k = 1
do i = 1, size(n_atoms) - 1
    do j = sum(n_atoms(1:i-1)) + 1 , sum(n_atoms(1:i))
        coor8(k,:)      = coor(j,:)
        coor8(k+1,:)    = coor(j,:) - (/1.0d0,0.0d0,0.0d0/)
        coor8(k+2,:)    = coor(j,:) - (/0.0d0,1.0d0,0.0d0/)
        coor8(k+3,:)    = coor(j,:) - (/0.0d0,0.0d0,1.0d0/)
        coor8(k+4,:)    = coor(j,:) - (/1.0d0,1.0d0,0.0d0/)
        coor8(k+5,:)    = coor(j,:) - (/1.0d0,0.0d0,1.0d0/)
        coor8(k+6,:)    = coor(j,:) - (/0.0d0,1.0d0,1.0d0/)
        coor8(k+7,:)    = coor(j,:) - (/1.0d0,1.0d0,1.0d0/)

        k = k + 8
    enddo

enddo


end subroutine

end module
