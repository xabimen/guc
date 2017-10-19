!This is file : inverse.f90

module matrix_operations

implicit none

public            :: inverse, cross_product, M33DET

contains

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
integer, intent(in)                         :: n
real(kind=8), dimension(:,:), intent(in)    :: a
real(kind=8), dimension(:,:), intent(out)   :: c
real(kind=8), dimension(n,n)                :: L, U, a2
real(kind=8), dimension(n)                  :: b, d, x
real(kind=8)                                :: coeff
integer                                     :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0
a2 = a
! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a2(i,k)/a2(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a2(i,j) = a2(i,j)-coeff*a2(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a2(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

function cross_product(a,b) result(ans)
real(kind=8), dimension(:), intent(in)      :: a, b
real(kind=8), dimension(size(a))            :: ans

ans(1) = a(2)*b(3) - b(2)*a(3)
ans(2) = a(3)*b(1) - b(3)*a(1)
ans(3) = a(1)*b(2) - b(1)*a(2)

end function

function projection(r1,r2) result(ans)
real(kind=8), dimension(:), intent(in)      :: r1, r2
real(kind=8)                                :: ans

ans = dot_product(r1,r2)/norm2(r2)

end function


FUNCTION M33DET (A) RESULT (DET)

      IMPLICIT NONE

      real(kind=8), DIMENSION(3,3), INTENT(IN)  :: A

      real(kind=8) :: DET


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      RETURN

END FUNCTION M33DET


end module
