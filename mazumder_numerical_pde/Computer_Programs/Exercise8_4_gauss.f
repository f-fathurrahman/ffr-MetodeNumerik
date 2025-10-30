! Start of Main Program
! This is a program for integrating a function in [1,2] using
! Gaussian quadrature
! Exercise 8.4 (d)

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 5
    INTEGER(4) :: i,j,k
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), DIMENSION(n) :: x,w
    REAL(8) :: igauss,iexact,sumf,c1,c2,func,zi,up,down
    
! Open output file
    OPEN(unit=10, file="Exercise8_3_gauss.out", status = "unknown")

! Input data
    x(1) = 0.0d0; w(1) = 0.568888888888889
    x(2) = 0.538469310105683d0; w(2) = 0.478628670499366
    x(3) = 0.906179845938664d0; w(3) = 0.236926885056189
    x(4) = -0.538469310105683d0; w(4) = 0.478628670499366
    x(5) = -0.906179845938664d0; w(5) = 0.236926885056189
   
! Compute Integral
    up = two; down = zero
    sumf = zero
    c1 = (up-down)/two
    c2 = (up+down)/two
    DO i = 1,n
      zi = c1*x(i)+c2
      func = exp(zi)/((exp(zi)-one)**2)
      sumf = sumf + w(i)*func
    ENDDO

    igauss = c1*sumf

    iexact = ((exp(down)+one)/(exp(down)-one)-(exp(up)+one)/(exp(up)-one))/two

    WRITE(10,10) iexact,igauss,100.0d0*(iexact-igauss)/iexact
    WRITE(*,10) iexact,igauss,100.0d0*(iexact-igauss)/iexact

 10 format(3(F14.8,1x))

    CLOSE(unit=10)

    END Program main
