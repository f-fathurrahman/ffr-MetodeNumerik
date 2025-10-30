! Start of Main Program
! This is a program for integrating a function in [1,2] using 
! spline interpolation.
! Exercise 8.3 (c)

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 5
    INTEGER(4) :: i,j,k
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), DIMENSION(n) :: x,f,eta,dd,rr,uu,del,aa,cc,yy
    REAL(8) :: dx,ispline,iexact,sumf
    REAL(8) :: a,b,c,d,num,den,w(n),up,down
    
! Open output file
    OPEN(unit=10, file="Exercise8_3_spline.out", status = "unknown")

! Input data
    up = two; down = one
    dx = (up-down)/(n-1)
    DO i = 1,n
      x(i) = down+(i-1)*dx
      eta(i) = x(i)
      f(i) = exp(x(i))/((exp(x(i))-one)**2)
    ENDDO

! Compute spline coeeficients
    DO i = 1,n-1
      del(i) = eta(i+1)-eta(i)
      yy(i) = f(i+1)-f(i)
    ENDDO
    
    dd(:)=zero; aa(:)=zero;cc(:)=zero;rr(:)=zero;uu(:)=zero
    dd(1) = one
    DO i = 2,n-1
     dd(i) = two*(del(i)+del(i-1))
     cc(i) = del(i)
     aa(i-1) = del(i-1)
     rr(i) = 6.0d0*(yy(i)/del(i)-yy(i-1)/del(i-1))
    ENDDO
    dd(n) = one

    CALL TRI(n,aa,dd,cc,rr,uu)

! Compute Integral
    
    sumf = zero
    DO i = 1,n-1
      sumf = sumf + (f(i+1)+f(i))*dx/two
    ENDDO
    DO i = 1,n-1
      sumf = sumf - (uu(i+1)+uu(i))*(dx**3)/24.0d0
    ENDDO

    ispline = sumf

    iexact = ((exp(one)+one)/(exp(one)-one)-(exp(two)+one)/(exp(two)-one))/two

    WRITE(10,10) iexact,ispline,100.0d0*(iexact-ispline)/iexact
    WRITE(*,10) iexact,ispline,100.0d0*(iexact-ispline)/iexact

 10 format(3(F14.8,1x))

    CLOSE(unit=10)

    END Program main

!--------------------------------------------------------------------------

! Tridiagonal solver (uses Thomas Algorithm)

    SUBROUTINE TRI(n,a,d,c,b,x)

    INTEGER(4), INTENT(in) :: n
    INTEGER(4) :: i
    REAL(8) :: xmult
    REAL(8), DIMENSION(n) :: a,b,c,d,x

    DO i = 2,n
      xmult = a(i-1)/d(i-1)
      d(i) = d(i) - xmult*c(i-1)
      b(i) = b(i) - xmult*b(i-1)
    ENDDO

    x(n) = b(n)/d(n)
    DO i = n-1,1,-1
      x(i) = (b(i)-c(i)*x(i+1))/d(i)
    ENDDO

    END SUBROUTINE TRI

