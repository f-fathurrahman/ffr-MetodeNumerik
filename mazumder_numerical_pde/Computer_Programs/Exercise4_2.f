! Start of Main Program
! 3 by 3 equation solution using Gauss-Seidel
! Exercise 4.2

    program main

! Declaration of variables
    IMPLICIT NONE

    INTEGER(4), PARAMETER :: max_iter = 100
    INTEGER(4) :: i,iter
    REAL(8) :: tol = 1.0d-6
    REAL(8) :: a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3
    REAL(8) :: xn,yn,zn,xo,yo,zo,res,sumres
    

! Open output file
    OPEN(unit=10, file="Exercise4_2.rsl", status = "unknown")

! Link Coefficients
    a11 = 5.0d0; a12 = -1.0d0; a13 = -3.0d0; b1 = 1.0d0
    a21 = -2.0d0; a22 = 6.0d0; a23 = -2.0d0; b2 = 2.0d0
    a31 = -1.0d0; a32 = 0.0d0; a33 = 2.0d0; b3 = 1.0d0


! Start of iteration loop

    xn = 0.0d0; yn = 0.0d0; zn = 0.0d0

    DO iter = 1,max_iter

      xn = (b1-a12*yn-a13*zn)/a11
      yn = (b2-a21*xn-a23*zn)/a22
      zn = (b3-a31*xn-a32*yn)/a33

      sumres = (a11*xn+a12*yn+a13*zn-b1)**2 + &
               (a21*xn+a22*yn+a23*zn-b2)**2 + &
               (a31*xn+a32*yn+a33*zn-b3)**2

      res = SQRT(sumres)

      WRITE(*,10) iter, res, xn,yn,zn
      WRITE(10,10) iter, res, xn,yn,zn

     IF(res < 1.0D-6) EXIT

    ENDDO  !iterations
     

    10 format(I4,4(1x,F13.6))

    CLOSE(unit=10)

    END Program main
