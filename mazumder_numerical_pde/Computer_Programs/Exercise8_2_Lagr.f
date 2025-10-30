! Start of Main Program
! This is a program for finding the Lagrange polynomial
! fit to the function f(x) = 1/x.
! Exercise 8.2 (b)

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 10, m = 1001
    INTEGER(4) :: i,j,k
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), DIMENSION(n) :: x,y,l
    REAL(8) :: xx,dx,yexact,ypoly,sump,prod
    
! Open output file
    OPEN(unit=10, file="Exercise8_2_Lagr.out", status = "unknown")

! Input data
    DO i = 1,n
      x(i) = one*i
      y(i) = one/x(i)
    ENDDO

! Print results
    dx = (10.5d0-0.5d0)/(m-1)
    DO k = 1,m
      xx = 0.5d0 + (k-1)*dx
      yexact = one/xx
      sump = zero
      DO i = 1,n
        prod = one
        DO j = 1,n
          IF(j == i)CYCLE
          prod = prod*(xx-x(j))/(x(i)-x(j))
        ENDDO
        l(i) = prod
        sump = sump + y(i)*l(i)
      ENDDO
      ypoly = sump
      WRITE(10,10) xx,yexact,ypoly,100.0d0*(yexact-ypoly)/yexact
      WRITE(*,10) xx,yexact,ypoly,100.0d0*(yexact-ypoly)/yexact
    ENDDO

 10 format(F13.4,1x,2(F14.8,1x),F14.6)

    CLOSE(unit=10)

    END Program main
