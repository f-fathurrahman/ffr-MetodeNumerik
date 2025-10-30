! Start of Main Program
! This is a program for integrating a function in [1,2] using the
! Trapezoidal method.
! Exercise 8.3 (a)

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 5
    INTEGER(4) :: i,j,k
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0,tiny=1.0d-10
    REAL(8), DIMENSION(n) :: x,f
    REAL(8) :: dx,itrap,iexact,sumf,up,down
    
! Open output file
    OPEN(unit=10, file="Exercise8_3_trap.out", status = "unknown")

! Input data
    up = two; down = tiny
    dx = (up-down)/(n-1)
    DO i = 1,n
      x(i) = down+(i-1)*dx
      f(i) = exp(x(i))/((exp(x(i))-one)**2)
    ENDDO

! Compute Integral
    
    sumf = zero
    DO i = 1,n-1
      sumf = sumf + dx*(f(i+1)+f(i))
    ENDDO

    itrap = sumf/two

    iexact = ((exp(down)+one)/(exp(down)-one)-(exp(up)+one)/(exp(up)-one))/two

    WRITE(10,*) iexact,itrap,100.0d0*(iexact-itrap)/iexact
    WRITE(*,*) iexact,itrap,100.0d0*(iexact-itrap)/iexact

 10 format(3(F14.8,1x))

    CLOSE(unit=10)

    END Program main
