! Start of Main Program
! This is a probram to solve the transcendental equation x*TAN(x) = 1
! Exercise 8.5

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: iter_max=1000
    INTEGER(4) :: iter,n
    REAL(8), PARAMETER :: one=1.0D0,half=0.5d0, &
                        & zero= 0.0d0,tolerance=1.0D-10,tenth=0.1D0
    REAL(8) :: x,dx,func,deriv,pi
    

! Opening output file

    OPEN(unit=10, file = "Exercise8_5.out", status = "unknown")
    WRITE(10,*) "Roots of x*TAN(x) = 1"

    pi = 4.0d0*datan(one)

    DO n = 1,10 ! Calculating the first 10 roots of the equation

      x = (n-1)*pi+tenth*pi    ! Initial guess for each root 

! Start of Iteration loop for Newton-Raphson Method
      DO iter = 1,iter_max
  
        func = x*sin(x) - cos(x)        ! Function
        deriv = 2.0d0*sin(x) + x*cos(x) ! Derivative of function
        dx = -func/deriv                  ! Change in z
        x = x + dx                        ! Update of z
        print *, n,iter
        IF(ABS(func) < tolerance) EXIT    ! Checking for tolerance
      ENDDO
      IF(iter > iter_max) print *, "Solution did not converge for n = ", n

!   Printing final solution to file roots.dat
      WRITE(10,*) n, x, (n-1)*pi  
      WRITE(*,*) n, x, (n-1)*pi

    ENDDO

    CLOSE(unit = 10)

    END Program main
