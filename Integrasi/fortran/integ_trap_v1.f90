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


SUBROUTINE test()
  IMPLICIT NONE 
  INTEGER :: N
  REAL(8) :: a, b, dx, x, integF
  REAL(8), ALLOCATABLE :: F(:)
  INTEGER :: i
  ! function
  REAL(8) :: test_func

  a = 0.d0
  b = 3.d0
  N = 21  ! definition of N is quite different from the one used in m_integ_trap
          ! to get similar result with m_integ_trap with N=20, we add 1 point (to match step size)
  dx = (b-a)/(N-1)
  ALLOCATE( F(N) )
  DO i = 1,N
    x = a + (i-1)*dx
    F(i) = test_func(x)
  ENDDO 

  CALL integ_trap( N, F, a, b, integF )
  WRITE(*,'(1x,I8,F18.10)') N, integF

  DEALLOCATE( F )
  
END SUBROUTINE 



FUNCTION num(x) RESULT(f)
  IMPLICIT NONE 
  REAL(8) :: x, f
  f = x**2 / ( exp(x) - 1.d0 )
END FUNCTION 


FUNCTION denum(x) RESULT(f)
  IMPLICIT NONE 
  REAL(8) :: x, f
  f = x**3 / ( exp(x) - 1.d0 )
END FUNCTION 


SUBROUTINE eval_u_x_g( x_g, integ_f )
  IMPLICIT NONE 
  ! arguments
  REAL(8) :: x_g
  REAL(8) :: integ_f
  ! local
  REAL(8) :: A, B
  REAL(8) :: x1, dx1, x2, dx2
  ! functions
  REAL(8) :: num, denum
  ! arrays
  REAL(8), ALLOCATABLE :: f_num(:), f_denum(:)
  !
  INTEGER :: N, i
  REAL(8) :: SMALL, LARGE

  SMALL = epsilon(0.d0)
  LARGE = 1.d4

  N = 30*int(LARGE)

  dx1 = (LARGE-x_g)/(N-1)
  dx2 = (LARGE-SMALL)/(N-1)

  ALLOCATE( f_num(N) )
  ALLOCATE( f_denum(N) )

  DO i = 1,N
    x1 = x_g + (i-1)*dx1
    f_num(i) = num(x1)
    !
    x2 = SMALL + (i-1)*dx2
    f_denum(i) = denum(x2)
  ENDDO 

  CALL integ_trap( N, f_num, x_g, LARGE, A )
  CALL integ_trap( N, f_denum, SMALL, LARGE, B ) 

  integ_f = x_g*A/B

END SUBROUTINE 


PROGRAM main
  IMPLICIT NONE 
  INTEGER, PARAMETER :: Npoints = 100
  INTEGER :: i
  REAL(8) :: x_start, x_end, dx, x, u_x_g

  x_start = epsilon(0.d0)  ! small number, but not zero
  x_end = 6.d0
  dx = (x_end - x_start)/(Npoints-1)

  DO i = 1,Npoints
    x = x_start + (i-1)*dx
    CALL eval_u_x_g( x, u_x_g )
    WRITE(*,*) x, u_x_g
  ENDDO 

END PROGRAM 
