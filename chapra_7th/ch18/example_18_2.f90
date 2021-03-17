program example
  implicit none
  integer :: N
  real(8) :: x(3), y(3)
  real(8) :: xi
  real(8), allocatable :: ea(:)
  real(8), allocatable :: yint(:)
  !
  real(8), allocatable :: fdd(:,:)

  x = (/ 1.d0, 4.d0, 6.d0 /)
  y = (/ 0.d0, 1.386294d0, 1.791759d0 /)

  N = 2
  xi = 2.d0
  
  allocate( ea(0:N-1) )
  allocate( yint(0:N) )

  call newton_interp(N, x, y, xi, yint, ea)
  
  write(*,*) 'yint = ', yint
  write(*,*) 'ea   = ', ea

  deallocate( ea )
  deallocate( yint )

end program 