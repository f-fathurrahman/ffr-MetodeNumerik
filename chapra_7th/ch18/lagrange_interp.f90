subroutine lagrange_interp(N, x, y, xx, res)
  implicit none
  integer :: N 
  real(8) :: x(0:N), y(0:N)
  real(8) :: xx, res
  !
  real(8) :: ss, pp
  integer :: i, j

  ss = 0.d0
  do i = 0,N
    pp = y(i)
    do j = 0,N
      if( i /= j ) then
        pp = pp*( xx - x(j) ) / ( x(i) - x(j) )
      endif
    enddo
    ss = ss + pp
  enddo
  res = ss
  return
end subroutine 