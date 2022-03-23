PROGRAM example_18_3
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: x(4), y(4)
  REAL(8) :: xi
  REAL(8), ALLOCATABLE :: ea(:)
  REAL(8), ALLOCATABLE :: yint(:)
  !
  REAL(8), ALLOCATABLE :: fdd(:,:)

  x = (/ 1.d0, 4.d0, 6.d0 , 5.d0 /)
  y = (/ 0.d0, 1.386294d0, 1.791759d0, 1.609438d0 /)

  N = 3
  xi = 2.d0
  
  ALLOCATE( ea(0:N-1) )
  ALLOCATE( yint(0:N) )

  CALL newton_interp(N, x, y, xi, yint, ea)
  
  WRITE(*,*) 'yint = ', yint
  WRITE(*,*) 'ea   = ', ea

  WRITE(*,*) abs( yint(N) - log(xi) )

  DEALLOCATE( ea )
  DEALLOCATE( yint )

END PROGRAM 

