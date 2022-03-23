SUBROUTINE newton_interp(N, x, y, xi, yint, ea)
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: x(0:N)
  REAL(8) :: y(0:N)
  REAL(8) :: ea(0:N-1)
  REAL(8) :: yint(0:N)
  REAL(8) :: xi
  !
  REAL(8) :: fdd(0:N,0:N)
  INTEGER :: i, j, order
  REAL(8) :: yint2, xterm

  DO i = 0,n
    fdd(i,0) = y(i)
  END DO 

  DO j = 1,n
    DO i = 0, n-j
      fdd(i,j) = ( fdd(i+1,j-1) - fdd(i,j-1) ) / ( x(i+j) - x(i) )
    END DO 
  END DO 

  xterm = 1.d0
  yint(0) = fdd(0,0)
  DO order = 1,N
    xterm = xterm * ( xi - x(order-1) )
    yint2 = yint(order-1) + fdd(0,order)*xterm
    ea(order-1) = yint2 - yint(order-1)
    yint(order) = yint2 
  END DO 

  RETURN 
END  SUBROUTINE 

