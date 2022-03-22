program example
  implicit none
  integer :: N_input
  integer :: N
  real(8), allocatable :: x(:), y(:)
  real(8), allocatable :: e(:), f(:), g(:), r(:)
  real(8), allocatable :: d2x(:)
  real(8) :: xu, yu, dy, d2y, yu_true
  integer :: NptsPlot
  real(8) :: xi, xf, dx
  integer :: i
  real(8), parameter :: pi=4.d0*atan(1.d0)
  real(8), parameter :: L = 1.d0
  
  N_input = 6
  allocate(x(0:N_input), y(0:N_input))
  dx = L/dble(N_input)
  write(*,*) 'dx = ', dx
  !x = (/ 0.d0, 0.3d0*L, 0.8d0*L, L /)
  do i = 0,N_input
    x(i) = i*dx
    y(i) = cos(2.d0*pi*x(i)/L)
  enddo

  N = size(x) - 1
  allocate( e(N-1), f(N-1), g(N-1), r(N-1) )

  write(*,*) 'Data: '
  do i = 0,N_input
    write(*,*) x(i), y(i)
  enddo

  call gen_trid_matrix( N, x, y, e, f, g, r )

  write(*,*) 'After gen_trid_matrix:'
  write(*,*) 'e = ', e
  write(*,*) 'f = ', f
  write(*,*) 'g = ', g
  write(*,*) 'r = ', r

  call decomp_trid(N-1, e, f, g)
  write(*,*) 'After decomp_trid'
  write(*,*) 'e = ', e
  write(*,*) 'f = ', f
  write(*,*) 'g = ', g
  write(*,*) 'r = ', r

  allocate( d2x(0:N) )
  ! Natural spline
  d2x(0) = 0.d0
  d2x(N) = 0.d0
  call subs_trid( N-1, e, f, g, r, d2x(1:N-1) ) ! only solve for N-1 unknowns

  write(*,*) 'd2x = ', d2x ! last element should be zero
  !stop 'ffr 51'

  xu = 0.51d0
  call interp_nat_cubic_spline( N, x, y, d2x, xu, yu, dy, d2y )
  write(*,*) 'yu  = ', yu
  write(*,*) 'dy  = ', dy
  write(*,*) 'd2y = ', d2y
  stop 'ffr 60'

  xi = x(0) ! initial point in the whole interval
  xf = x(N) ! final point in the whole interval
  NptsPlot = 50
  dx = (xf - xi)/(NptsPlot-1)
  do i = 1,NptsPlot
    xu = xi + (i-1)*dx
    call interp_nat_cubic_spline( N, x, y, d2x, xu, yu, dy, d2y )
    yu_true = cos(2.d0*pi*xu/L)
    write(*,'(1x,3F18.10)') xu, yu, yu_true
  enddo

  deallocate( x, y )
  deallocate( e, f, g, r )
  deallocate( d2x )

end program