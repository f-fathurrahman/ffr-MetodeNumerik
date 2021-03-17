program example
  implicit none
  integer :: N
  real(8) :: x(8), y(8)
  real(8) :: xi
  real(8), allocatable :: ea(:)
  real(8), allocatable :: yint(:)
  !
  real(8), allocatable :: fdd(:,:)
  integer, allocatable :: idx_choose(:)

  x = (/ 1.d0, 4.d0, 6.d0 , 5.d0, 3.d0, 1.5d0, 2.5d0, 3.5d0 /)
  y = (/ 0.d0, 1.386294d0, 1.791759d0, 1.609438d0, 1.0986123d0, &
         0.4054641d0, 0.9162907d0, 1.2527630d0 /)

  N = 3
  allocate( idx_choose(N) )
  idx_choose = (/ 1, 6, 7, 5 /)
  xi = 2.d0
  
  allocate( ea(0:N-1) )
  allocate( yint(0:N) )

  write(*,*) 'using the following data:'
  write(*,*) 'x = ', x(idx_choose)
  write(*,*) 'y = ', y(idx_choose)
  write(*,*)
  call newton_interp(N, x(idx_choose), y(idx_choose), xi, yint, ea)
  
  write(*,*)
  write(*,*) 'yint = ', yint
  write(*,*) 'ea   = ', ea

  write(*,*) abs( yint(N) - log(xi) )

  deallocate(idx_choose)
  deallocate( ea )
  deallocate( yint )

end program 