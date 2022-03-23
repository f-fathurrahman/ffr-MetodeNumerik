SUBROUTINE interp_nat_cubic_spline( N, x, y, d2x, xu, yu, dy, d2y )
  !
  IMPLICIT NONE 
  !
  INTEGER :: N
  REAL(8) :: x(0:N), y(0:N)
  REAL(8) :: d2x(0:N)
  REAL(8) :: xu
  REAL(8) :: yu
  REAL(8) :: dy, d2y
  !
  INTEGER :: i
  LOGICAL :: flag, is_in_interval
  REAL(8) :: c1, c2, c3, c4, t1, t2, t3, t4

  flag = .false.
  i = 1
  DO WHILE(.true.)
    !
    is_in_interval = (xu >= x(i-1)) .and. (xu <= x(i))
    !
    IF( is_in_interval ) THEN 
      !
      !write(*,*) 'interval = ', i
      !
      c1 = d2x(i-1)/6.d0/( x(i) - x(i-1) )
      c2 = d2x(i)/6.d0/( x(i) - x(i-1) )
      c3 = y(i-1)/(x(i) - x(i-1)) - d2x(i-1)*(x(i) - x(i-1))/6.d0
      c4 = y(i)/(x(i) - x(i-1)) - d2x(i)*(x(i) - x(i-1))/6.d0
      !
      t1 = c1*( x(i) - xu )**3
      t2 = c2*( xu - x(i-1) )**3
      t3 = c3*( x(i) - xu )
      t4 = c4*( xu - x(i-1) )
      !
      yu = t1 + t2 + t3 + t4
      !
      t1 = -3.d0*c1*( x(i) - xu )**2
      t2 = 3.d0*c2*( xu - x(i-1) )**2
      t3 = -c3
      t4 = c4
      !
      dy = t1 + t2 + t3 + t4
      !
      t1 = 6.d0*c1*( x(i) - xu )
      t2 = 6.d0*c2*( xu - x(i-1) )
      !
      d2y = t1 + t2
      !
      flag = .true.
      !
    ELSE 
      ! 
      i = i + 1
    END IF 
    !
    IF( i == N + 1 .or. flag ) THEN 
      EXIT ! break out of the loop
    END IF 
  END DO 

  IF( .not. flag ) THEN 
    WRITE(*,*) "ERROR: xu is outside range of spline"
    RETURN
  END IF

  RETURN 
END SUBROUTINE 


SUBROUTINE decomp_trid(N, e, f, g)
  IMPLICIT  NONE 
  INTEGER :: N
  REAL(8) :: e(N), f(N), g(N)
  INTEGER :: k
  DO k = 2,N
    e(k) = e(k)/f(k-1)
    f(k) = f(k) - e(k)*g(k-1)
  END DO 
  RETURN 
END SUBROUTINE 

! should be called after calling decomp_trid
SUBROUTINE subs_trid(N, e, f, g, r, x)
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: e(N), f(N), g(N), r(N)
  REAL(8) :: x(N) ! output
  INTEGER :: k

  ! Forward subs
  DO k = 2,N
    r(k) = r(k) - e(k)*r(k-1)
  END DO 

  ! back subs
  x(N) = r(N)/f(N)
  DO k = N-1,1,-1
    x(k) = ( r(k) - g(k)*x(k+1) ) / f(k)
  END DO 
  RETURN 
END SUBROUTINE 

SUBROUTINE gen_trid_matrix(N, x, y, e, f, g, r)
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: x(0:N), y(0:N)
  REAL(8) :: e(N-1), f(N-1), g(N-1), r(N-1)
  INTEGER :: i

  f(1) = 2.d0*( x(2) - x(0) )
  g(1) = x(2) - x(1)
  r(1) = 6.d0/( x(2) - x(1) ) * ( y(2) - y(1) )
  r(1) = r(1) + 6.d0/( x(1) - x(0) ) * ( y(0) - y(1) )
  !
  DO i = 2,N-2
    e(i) = x(i) - x(i-1)
    f(i) = 2.d0*( x(i+1) - x(i-1) )
    g(i) = x(i+1) - x(i)
    r(i) = 6.d0/( x(i+1) - x(i) ) * ( y(i+1) - y(i) )
    r(i) = r(i) + 6.d0/( x(i) - x(i-1) ) * ( y(i-1) - y(i) )
  END DO 
  !
  e(n-1) = x(n-1) - x(n-2)
  f(n-1) = 2.d0*( x(n) - x(n-2) )
  r(n-1) = 6.d0/( x(n) - x(n-1) ) * ( y(n) - y(n-1) )
  r(n-1) = r(n-1) + 6.d0/( x(n-1) - x(n-2) ) * ( y(n-2) - y(n-1) )

  RETURN 
END SUBROUTINE 

