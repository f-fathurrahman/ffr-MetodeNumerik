subroutine interp_nat_cubic_spline( N, x, y, d2x, xu, yu, dy, d2y )
  !
  implicit none 
  !
  integer :: N
  real(8) :: x(0:N), y(0:N)
  real(8) :: d2x(0:N)
  real(8) :: xu
  real(8) :: yu
  real(8) :: dy, d2y
  !
  integer :: i
  logical :: flag, is_in_interval
  real(8) :: c1, c2, c3, c4, t1, t2, t3, t4

  !write(*,*) 'size d2x = ', size(d2x)
  !write(*,*) 'd2x = ', d2x

  flag = .false.
  i = 1
  do while(.true.)
    !
    is_in_interval = (xu >= x(i-1)) .and. (xu <= x(i))
    !
    if( is_in_interval ) then 
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
    else 
      ! 
      i = i + 1
    endif 
    !
    if( i == N + 1 .or. flag ) then 
      exit ! break out of the loop
    endif
  enddo 

  if( .not. flag ) then
    write(*,*) "ERROR: xu is outside range of spline"
    return
  endif

  return
end subroutine 


subroutine decomp_trid(N, e, f, g)
  implicit none
  integer :: N
  real(8) :: e(N), f(N), g(N)
  integer :: k
  do k = 2,N
    e(k) = e(k)/f(k-1)
    f(k) = f(k) - e(k)*g(k-1)
  enddo
  return
end subroutine

! should be called after calling decomp_trid
subroutine subs_trid(N, e, f, g, r, x)
  implicit none
  integer :: N
  real(8) :: e(N), f(N), g(N), r(N)
  real(8) :: x(N) ! output
  integer :: k

  ! Forward subs
  do k = 2,N
    r(k) = r(k) - e(k)*r(k-1)
  enddo

  ! back subs
  x(N) = r(N)/f(N)
  do k = N-1,1,-1
    x(k) = ( r(k) - g(k)*x(k+1) ) / f(k)
  enddo
  return
end subroutine

subroutine gen_trid_matrix(N, x, y, e, f, g, r)
  implicit none
  integer :: N
  real(8) :: x(0:N), y(0:N)
  real(8) :: e(N-1), f(N-1), g(N-1), r(N-1)
  integer :: i

  f(1) = 2.d0*( x(2) - x(0) )
  g(1) = x(2) - x(1)
  r(1) = 6.d0/( x(2) - x(1) ) * ( y(2) - y(1) )
  r(1) = r(1) + 6.d0/( x(1) - x(0) ) * ( y(0) - y(1) )
  !
  do i = 2,N-2
    e(i) = x(i) - x(i-1)
    f(i) = 2.d0*( x(i+1) - x(i-1) )
    g(i) = x(i+1) - x(i)
    r(i) = 6.d0/( x(i+1) - x(i) ) * ( y(i+1) - y(i) )
    r(i) = r(i) + 6.d0/( x(i) - x(i-1) ) * ( y(i-1) - y(i) )
  enddo
  !
  e(n-1) = x(n-1) - x(n-2)
  f(n-1) = 2.d0*( x(n) - x(n-2) )
  r(n-1) = 6.d0/( x(n) - x(n-1) ) * ( y(n) - y(n-1) )
  r(n-1) = r(n-1) + 6.d0/( x(n-1) - x(n-2) ) * ( y(n-2) - y(n-1) )


  return
end subroutine