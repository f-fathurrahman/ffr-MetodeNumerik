subroutine newton_interp(N, x, y, xi, yint, ea)
  implicit none
  integer :: N
  real(8) :: x(0:N)
  real(8) :: y(0:N)
  real(8) :: ea(0:N-1)
  real(8) :: yint(0:N)
  real(8) :: xi
  !
  real(8) :: fdd(0:N,0:N)
  integer :: i, j, order
  real(8) :: yint2, xterm

  do i = 0,n
    fdd(i,0) = y(i)
  enddo

  do j = 1,n
    do i = 0, n-j
      fdd(i,j) = ( fdd(i+1,j-1) - fdd(i,j-1) ) / ( x(i+j) - x(i) )
    enddo
  enddo

  xterm = 1.d0
  yint(0) = fdd(0,0)
  do order = 1,N
    xterm = xterm * ( xi - x(order-1) )
    yint2 = yint(order-1) + fdd(0,order)*xterm
    ea(order-1) = yint2 - yint(order-1)
    yint(order) = yint2 
  enddo

  return
end subroutine