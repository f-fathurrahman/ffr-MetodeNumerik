! Start of Main Program
! This is a program for finding the natural cubic spline
! fit to the function f(x) = 1/x.
! Exercise 8.1 (c)

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 10, m = 901
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), DIMENSION(n) :: x,y
    REAL(8), DIMENSION(n) :: dd,aa,cc,rr,uu,del,yy,eta
    REAL(8) :: sx,dx,xx,yexact
    
! Open output file
    OPEN(unit=10, file="Exercise8_1_spline.out", status = "unknown")

! Input data
    DO i = 1,n
      x(i) = one*i
      y(i) = one/x(i)
      eta(i) = x(i)
    ENDDO

    DO i = 1,n-1
      del(i) = eta(i+1)-eta(i)
      yy(i) = y(i+1)-y(i)
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
       
! Print results
    dx = (x(n)-x(1))/(m-1)
    DO j = 1,m
      xx = x(1) + (j-1)*dx
      yexact = one/xx
      DO i = 1,n-1
        IF( xx < eta(i+1) .and. xx >= eta(i)) THEN
          sx = uu(i)*((eta(i+1)-xx)**3)/(6*del(i)) + uu(i+1)*((xx-eta(i))**3)/(6*del(i)) + &
                   (y(i+1)/del(i)-uu(i+1)*del(i)/6)*(xx-eta(i)) + &
                   (y(i)/del(i)-uu(i)*del(i)/6)*(eta(i+1)-xx)
          WRITE(10,10) xx,yexact,sx,100*(yexact-sx)/yexact
          WRITE(*,10) xx,yexact,sx,100*(yexact-sx)/yexact
        ENDIF
      ENDDO
    ENDDO

 10 format(F13.4,1x,2(F14.8,1x),F14.6)

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

