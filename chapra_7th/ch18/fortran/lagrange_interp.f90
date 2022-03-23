SUBROUTINE lagrange_interp(N, x, y, xx, res)
  IMPLICIT NONE 
  INTEGER :: N 
  REAL(8) :: x(0:N), y(0:N)
  REAL(8) :: xx, res
  !
  REAL(8) :: ss, pp
  INTEGER :: i, j

  ss = 0.d0
  DO i = 0,N
    pp = y(i)
    DO j = 0,N
      IF( i /= j ) THEN 
        pp = pp*( xx - x(j) ) / ( x(i) - x(j) )
      END IF 
    enddo
    ss = ss + pp
  END DO 
  res = ss
  RETURN 
END SUBROUTINE 

