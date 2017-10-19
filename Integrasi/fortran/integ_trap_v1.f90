SUBROUTINE integ_trap( N, F, a, b, integF )
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: F(N)
  REAL(8) :: h, a, b
  REAL(8) :: integF
  INTEGER :: i

  h = (b-a)/(N-1)
  integF = F(1) + F(N)
  DO i = 2,N-1
    integF = integF + 2.d0*F(i)
  ENDDO 
  integF = 0.5d0*h*integF
END SUBROUTINE 


FUNCTION test_func(x) RESULT(f)
  IMPLICIT NONE 
  REAL(8) :: x, f
  f = sqrt(9.81d0*68.1d0/0.25d0)*tanh(sqrt(9.81d0*0.25d0/68.1d0)*x)
END FUNCTION 

PROGRAM test
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: a, b, dx, x, integF
  REAL(8), ALLOCATABLE :: F(:)
  INTEGER :: i
  ! function
  REAL(8) :: test_func

  a = 0.d0
  b = 3.d0
  N = 21
  dx = (b-a)/(N-1)
  ALLOCATE( F(N) )
  DO i = 1,N
    x = a + (i-1)*dx
    F(i) = test_func(x)
  ENDDO 

  CALL integ_trap( N, F, a, b, integF )
  WRITE(*,'(1x,I8,F18.10)') N, integF

  DEALLOCATE( F )
  
END PROGRAM 
