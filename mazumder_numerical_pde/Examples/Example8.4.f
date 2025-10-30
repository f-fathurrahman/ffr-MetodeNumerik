! Start of Main Program
! Example 8.4

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n = 5, m = 401
    INTEGER(4) :: i,j,iter
    REAL(8), PARAMETER :: one = 1.0D0,two = 2.0d0, &
                             & zero = 0.0d0,half = 0.5d0
    REAL(8), DIMENSION(n) :: x,y
    REAL(8), DIMENSION(n) :: dd,aa,cc,rr,uu,del,yy,eta
    REAL(8) :: sx,dx,xx,sp,sumi
    
! Open output file
    OPEN(unit=10, file="Example8_4.out", status = "unknown")

    x(1) = -2; x(2) = -1; x(3) = 0; x(4) = 1; x(5) = 2
    DO i = 1,n
      y(i) = exp(x(i))*x(i)
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
      DO i = 1,n-1
        IF( xx <= eta(i+1) .and. xx >= eta(i)) THEN
          sx = uu(i)*((eta(i+1)-xx)**3)/(6*del(i)) + uu(i+1)*((xx-eta(i))**3)/(6*del(i)) + &
                   (y(i+1)/del(i)-uu(i+1)*del(i)/6)*(xx-eta(i)) + &
                   (y(i)/del(i)-uu(i)*del(i)/6)*(eta(i+1)-xx)
          WRITE(10,*) xx,sx,exp(xx)*xx
          WRITE(*,*) xx,sx,exp(xx)*xx
        ENDIF
      ENDDO
    ENDDO

    print *, uu

! Compute Integral
    sumi = zero
    DO i = 1,n-1
      sumi = sumi - (del(i)**3)*(uu(i+1)+uu(i))/24.0d0 + del(i)*(y(i+1)+y(i))/2.0d0
    ENDDO
    print *, sumi

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

