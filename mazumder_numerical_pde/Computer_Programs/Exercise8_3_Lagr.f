! Start of Main Program
! This is a program for integrating a function in [1,2] using the
! Lagarnge polynomial formula.
! Exercise 8.3 (b)

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 5
    INTEGER(4) :: i,j,k
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), DIMENSION(n) :: x,f
    REAL(8) :: dx,iLagr,iexact,sumf
    REAL(8) :: a,b,c,d,num,den,w(n),up,down
    
! Open output file
    OPEN(unit=10, file="Exercise8_3_Lagr.out", status = "unknown")

! Input data
    up = two; down = one
    dx = (up-down)/(n-1)
    DO i = 1,n
      x(i) = down+(i-1)*dx
      f(i) = exp(x(i))/((exp(x(i))-one)**2)
    ENDDO

! Compute weights
    DO i = 1,n
      
      IF(i == 1) THEN
        a = x(2); b = x(3); c = x(4); d = x(5)
      ELSEIF(i == 2) THEN
        a = x(1); b = x(3); c = x(4); d = x(5)
      ELSEIF(i == 3) THEN
        a = x(1); b = x(2); c = x(4); d = x(5)
      ELSEIF(i == 4) THEN
        a = x(1); b = x(2); c = x(3); d = x(5)
      ELSEIF(i == 5) THEN
        a = x(1); b = x(2); c = x(3); d = x(4)
      ENDIF
    
      num = (up**5)/5.0d0 + a*b*c*d*up -(a+b+c+d)*(up**4)/4.0d0 + &
          (c*d+a*b+(a+b)*(c+d))*(up**3)/3.0d0 - &
          ((a+b)*c*d+(c+d)*a*b)*(up**2)/two - &
          ((down**5)/5.0d0 + a*b*c*d*down -(a+b+c+d)*(down**4)/4.0d0 + &
          (c*d+a*b+(a+b)*(c+d))*(down**3)/3.0d0 - &
          ((a+b)*c*d+(c+d)*a*b)*(down**2)/two)
      den = (x(i)-a)*(x(i)-b)*(x(i)-c)*(x(i)-d)
      w(i) = num/den
    ENDDO

! Compute Integral
    
    sumf = zero
    DO i = 1,n
      sumf = sumf + f(i)*w(i)
    ENDDO

    iLagr = sumf

    iexact = ((exp(one)+one)/(exp(one)-one)-(exp(two)+one)/(exp(two)-one))/two

    WRITE(10,10) iexact,iLagr,100.0d0*(iexact-iLagr)/iexact
    WRITE(*,10) iexact,iLagr,100.0d0*(iexact-iLagr)/iexact

 10 format(3(F14.8,1x))

    CLOSE(unit=10)

    END Program main
