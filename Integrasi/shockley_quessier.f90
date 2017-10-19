
! a module to wrap integ_trap function
MODULE m_integ_trap

CONTAINS 

FUNCTION integ_trap( F, a, b, N ) RESULT(integF)
  IMPLICIT NONE 
  ! arguments
  INTERFACE
    FUNCTION F(x) RESULT(s)
      REAL(8) :: s, x
    END FUNCTION 
  END INTERFACE 
  REAL(8) :: a, b
  INTEGER :: N
  ! local
  INTEGER :: i
  REAL(8) :: s, x, h
  REAL(8) :: integF
  
  h = (b-a)/N
  s = F(a)
  x = a
  DO i = 1, N-1
    x = x + h
    s = s + 2.d0*F(x)
  ENDDO 
  s = s + F(b)
  integF = 0.5d0*h*s
END FUNCTION 

END MODULE 


!
! a module to wrap various functions to be integrated
!
MODULE m_funcs

CONTAINS 

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


FUNCTION test_func(x) RESULT(f)
  IMPLICIT NONE 
  REAL(8) :: x, f
  f = sqrt(9.81d0*68.1d0/0.25d0)*tanh(sqrt(9.81d0*0.25d0/68.1d0)*x)
END FUNCTION 

END MODULE 



FUNCTION u_x_g( x_g ) RESULT(f)
  USE m_integ_trap
  USE m_funcs
  IMPLICIT NONE 
  REAL(8) :: x_g, f
  REAL(8) :: SMALL, LARGE
  REAL(8) :: A, B
  !
  SMALL = epsilon(0.d0) ! relative small but not zero
  LARGE = 1.d4
  A = integ_trap( num, x_g, LARGE, 30*int(LARGE) )
  B = integ_trap( denum, SMALL, LARGE, 30*int(LARGE) )
  f = x_g*A/B
END FUNCTION 

! a subroutine to test whether integ_trap is working correctly
SUBROUTINE test_integ
  USE m_integ_trap
  USE m_funcs
  IMPLICIT NONE 
  REAL(8) :: s
  INTEGER :: i, N
  
  DO i = 1,30
    N = i*20
    s = integ_trap( test_func, 0.d0, 3.d0, N )
    WRITE(*,*) i, s
  ENDDO 
END SUBROUTINE 


PROGRAM shockley_quessier
  USE m_integ_trap
  USE m_funcs
  IMPLICIT NONE 
  ! function
  REAL(8) :: u_x_g
  INTEGER, PARAMETER :: Npoints = 100
  INTEGER :: i
  REAL(8) :: x_start, x_end, dx, x
  REAL(8) :: SMALL, LARGE

  x_start = epsilon(0.d0)  ! small number, but not zero
  x_end = 6.d0
  dx = (x_end - x_start)/(Npoints-1)

  DO i = 1,Npoints
    x = x_start + (i-1)*dx
    WRITE(*,*) x, u_x_g(x)
  ENDDO 

!  SMALL = epsilon(1.d0) ! relative small but not zero
!  LARGE = 1.d4
!  WRITE(*,*) integ_trap( denum, SMALL, LARGE, 30*int(LARGE) )
END PROGRAM 

