! Start of Main Program
! Plotting Characteristic Function (Polynomial)
! Example 4.2.

    program main

! Declaration of variables
    IMPLICIT NONE

    INTEGER(4), PARAMETER :: max_iter = 3000
    INTEGER(4) :: i
    REAL(8) :: a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3,x,f,f1,f2
    REAL(8), PARAMETER :: xl = -1.0d0, dx = 1.0d-3
    

! Open output file
    OPEN(unit=10, file="Example4_2.out", status = "unknown")

! Link Coefficients
    a11 = 0.0d0; a12 = -2.0d0/5.0d0; a13 = -2.0d0/5.0d0; b1 = 9.0d0
    a21 = 1.0d0/3.0d0; a22 = 0.0d0; a23 = 1.0d0/2.0d0; b2 = -1.0d0
    a31 = -7.0d0; a32 = -2.0d0; a33 = 0.0d0; b3 = 10.0d0

! Construct polynomial

    DO i = 1,max_iter
      x = xl + (i-1)*dx
      f = (a11-x)*((a22-x)*(a33-x)-a23*a32) - a12*(a21*(a33-x)-a23*a31) + a13*(a21*a32-(a22-x)*a31)
      f1 = (x-1.638)*(x*x+1.638*x+1.013)
      WRITE(*,*) x, f,f1
      WRITE(10,*) x, f
    ENDDO

    CLOSE(unit=10)

    END Program main
