PROGRAM example_18_5
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: x(8), y(8)
  REAL(8) :: xi
  REAL(8), ALLOCATABLE :: ea(:)
  REAL(8), ALLOCATABLE :: yint(:)
  !
  REAL(8), ALLOCATABLE :: fdd(:,:)
  INTEGER, ALLOCATABLE :: idx_choose(:)

  x = (/ 1.d0, 4.d0, 6.d0 , 5.d0, 3.d0, 1.5d0, 2.5d0, 3.5d0 /)
  y = (/ 0.d0, 1.386294d0, 1.791759d0, 1.609438d0, 1.0986123d0, &
         0.4054641d0, 0.9162907d0, 1.2527630d0 /)

  N = 4
  ALLOCATE( idx_choose(N) )
  idx_choose = (/ 1, 6, 7, 5 /)
  WRITE(*,*) idx_choose(4)
  WRITE(*,*) 'idx_choose = ', idx_choose
  WRITE(*,*) 'size(idx_choose) = ', size(idx_choose)
  xi = 2.d0
  
  ALLOCATE( ea(0:N-1) )
  ALLOCATE( yint(0:N) )

  WRITE(*,*) 'using the following data:'
  WRITE(*,*) 'x = ', x(idx_choose)
  WRITE(*,*) 'y = ', y(idx_choose)
  WRITE(*,*)
  CALL newton_interp(N, x(idx_choose), y(idx_choose), xi, yint, ea)
  
  WRITE(*,*)
  WRITE(*,*) 'yint = ', yint
  WRITE(*,*) 'ea   = ', ea

  WRITE(*,*) abs( yint(N) - log(xi) )

  DEALLOCATE(idx_choose)
  DEALLOCATE( ea )
  DEALLOCATE( yint )

END PROGRAM  

