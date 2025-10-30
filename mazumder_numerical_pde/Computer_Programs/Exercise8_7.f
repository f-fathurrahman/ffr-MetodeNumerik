! Start of Main Program
! This is a program to solve a nonlinear 1D diffusion equation subject
! to two Dirichlet Boundary conditions using finite difference
! with direct Newton iterations
! Exercise 8.7

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 41 , max_iter = 100
    INTEGER(4) :: iter,n_eq,i,j
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,phi_left=0.0D0,phi_right=1.0D0, &
                             & length=1.0D0, tiny=1.0D-40, tol=1.0D-10, &
                             & con = 4.0d0
    REAL(8), DIMENSION(n_nodes) :: f,x,phi,rhs,dphi
    REAL(8), DIMENSION(n_nodes,n_nodes) :: jac
    REAL(8) :: dx,dx2,err_i,sumres,res,fde
    
! Open output file
    OPEN(unit=10, file="Exercise8_7.out", status = "unknown")
    OPEN(unit=20, file="Exercise8_7.rsl", status = "unknown")

    dx = length/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO

! Start of Newton Iteration
    phi(:) = zero   ! Initial Guess for phi

    DO iter = 1,max_iter

! Calculate function
      f(1) = phi(1) - phi_left
      DO i = 2,n_nodes-1
        f(i) = -two*phi(i)/dx2 + phi(i+1)/dx2 + phi(i-1)/dx2 - &
             & exp(con*phi(i))
      ENDDO
      f(n_nodes) = phi(n_nodes) - phi_right
      rhs(:) = -f(:) ! Right hand side of equations

! Calculate Jacobian (Coefficient Matrix)
      jac(:,:) = zero
      jac(1,1) = one
      DO i = 2,n_nodes-1
        DO j = 1, n_nodes
          IF(i == j) THEN
            jac(i,j) = -two/dx2 - con*exp(con*phi(i))
          ELSEIF(j == i+1 .or. j == i-1) THEN
            jac(i,j) = one/dx2
          ENDIF
        ENDDO
      ENDDO
      jac(n_nodes,n_nodes) = one

! Solve system of equations using Gaussian Elimination
      CALL gauss_elim(n_nodes,jac,rhs,dphi)

! Update Solution
      phi(:) = phi(:) + dphi(:)

! Calculate L2NORM
      sumres = zero
      DO i = 2,n_nodes-1
        fde = (phi(i+1)-two*phi(i)+phi(i-1))/dx2-exp(con*phi(i))
        sumres = sumres + fde*fde
      ENDDO
      res = sqrt(MAX(zero,sumres))

      WRITE(*,*) iter, res
      WRITE(20,*) iter, res

      IF(res < tol) EXIT    ! Convergence check

    ENDDO ! iteration loop

! Print results
    DO i = 1,n_nodes
      WRITE(*,*) x(i),phi(i)
      WRITE(10,*) x(i),phi(i)
    ENDDO

    CLOSE(unit=10)
    CLOSE(unit=20)

    END Program main

!-----------------------------------------------------------------------
  SUBROUTINE gauss_elim(n,a,b,x)
!-----------------------------------------------------------------------
!
! Purpose: Solves a set of linear algebraic equations.
!          Using Gaussian Elimination
!
!          n = number of equations (input)
!          a = nxn coefficient matrix (input)
!          b = right hand vector (input)
!          x = solution vector (output
!
!-----------------------------------------------------------------------

   IMPLICIT NONE
   INTEGER(4) :: i,j,k
   INTEGER(4) , INTENT(IN) :: n
   REAL(8) :: sum1,xmult
   REAL(8) :: b(n),a(n,n)
   REAL(8), INTENT(OUT) :: x(n)
!-----------------------------------------------------------------------

   DO k = 1,n-1
     DO i = k+1,n
       xmult = a(i,k)/a(k,k)
       DO j=k+1,n
         a(i,j) = a(i,j)-xmult*a(k,j)
       ENDDO
       a(i,k) = xmult
       b(i) = b(i) - xmult*b(k)
     ENDDO
   ENDDO
   x(n) = b(n)/a(n,n)
   DO i = n-1,1,-1
      sum1 = b(i)
      DO j = i+1,n
         sum1 = sum1 - a(i,j)*x(j)
      ENDDO
      x(i) = sum1/a(i,i)
   ENDDO

   END SUBROUTINE gauss_elim
