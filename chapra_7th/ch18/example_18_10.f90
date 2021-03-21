PROGRAM example
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: x(0:3), y(0:3)
  REAL(8), ALLOCATABLE :: e(:), f(:), g(:), r(:)
  REAL(8), ALLOCATABLE :: d2x(:)
  REAL(8) :: xu, yu, dy, d2y
  INTEGER :: NptsPlot
  REAL(8) :: xi, xf, dx
  INTEGER :: i

  x = (/ 3.d0, 4.5d0, 7.d0, 9.d0 /)
  y = (/ 2.5d0, 1.d0, 2.5d0, 0.5d0 /)

  N = 3
  ALLOCATE( e(N-1), f(N-1), g(N-1), r(N-1) )
  ALLOCATE( d2x(0:N) )
  ! Natural spline
  d2x(0) = 0.d0
  d2x(N) = 0.d0

  CALL gen_trid_matrix( N, x, y, e, f, g, r )

  !write(*,*) 'e = ', e
  !write(*,*) 'f = ', f
  !write(*,*) 'g = ', g
  !write(*,*) 'r = ', r

  CALL decomp_trid(N-1, e, f, g)
  CALL subs_trid( N-1, e, f, g, r, d2x(1:N-1) ) ! only solve for N-1 unknowns

  !write(*,*) 'd2x = ', d2x ! last element should be zero

  !xu = 7.1d0
  !call interp_nat_cubic_spline( N, x, y, d2x, xu, yu, dy, d2y )
  !write(*,*) 'yu = ', yu

  xi = x(0) ! initial point in the whole interval
  xf = x(N) ! final point in the whole interval
  NptsPlot = 50
  dx = (xf - xi)/(NptsPlot-1)
  DO i = 1,NptsPlot
    xu = 3.d0 + (i-1)*dx
    CALL interp_nat_cubic_spline( N, x, y, d2x, xu, yu, dy, d2y )
    WRITE(*,'(1x,2F18.10)') xu, yu
  END DO

  DEALLOCATE( e, f, g, r )
  DEALLOCATE( d2x )

END PROGRAM 

