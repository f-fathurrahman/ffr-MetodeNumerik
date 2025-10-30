! Start of Main Program
! This is a program to solve a second order linear elliptic PDE
! using the finite-difference method and Gaussian Elimination
! Exercise 2.6

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes=41, m_nodes=41, neq = n_nodes*m_nodes
    INTEGER(4) :: i,j,k
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,phi_left=0.0D0,phi_right=1.0D0, &
                             & length=1.0D0, tiny=1.0D-40, tol=1.0D-6, &
                             & con = 2.0d0, width=1.0d0
    REAL(8), DIMENSION(n_nodes) :: x
    REAL(8), DIMENSION(m_nodes) :: y
    REAL(8), DIMENSION(n_nodes,m_nodes) :: phia,phi
    REAL(8), DIMENSION(neq) :: phi1D,rhs
    REAL(8), DIMENSION(neq,neq) :: acoef
    REAL(8) :: dx,dx2,dy,dy2,err,sij,diag
    
! Open output file
    OPEN(unit=10, file="Exercise2_6.out", status = "unknown")

    dx = length/(n_nodes-1)
    dy = width/(m_nodes-1)
    dx2 = dx*dx
    dy2 = dy*dy
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
    ENDDO
    DO j = 1,m_nodes
      y(j) = (j-1)*dy
    ENDDO


! Boundary Conditions
    DO i = 1,n_nodes
      phi(i,1) = 100.0d0*x(i)+500.0d0*exp(-50.0*(one-x(i))**2) !bottom
      phi(i,m_nodes) = 500.0d0*exp(-50.0*((one-x(i))**2+one))  !top
    ENDDO
    DO j = 1,m_nodes
      phi(1,j) = 500.0d0*exp(-50.0*(one+y(j)**2)) !left
      phi(n_nodes,j) = 100*(one-y(j))+500.0d0*exp(-50.0*y(j)**2)
    ENDDO

! Set up linear system
    acoef(:,:) = zero; rhs(:) = zero

! Interior nodes
    DO i = 2,n_nodes-1
      DO j = 2,m_nodes-1
         k = (j-1)*n_nodes+i
         acoef(k,k) = two/dx2+two/dy2
         acoef(k,k+1) = -one/dx2
         acoef(k,k-1) = -one/dx2
         acoef(k,k+n_nodes) = -one/dy2
         acoef(k,k-n_nodes) = -one/dy2
         sij = 5.0D4*(100.0d0*((one-x(i))**2+y(j)**2)-two)* &
                exp(-50.0d0*((one-x(i))**2+y(j)**2))
         rhs(k) = -sij
      ENDDO
    ENDDO
! Left Boundary
    i = 1
    DO j = 1,m_nodes
       k = (j-1)*n_nodes+i
       acoef(k,k) = one
       rhs(k) = phi(i,j)
    ENDDO
! Right Boundary
    i = n_nodes
    DO j = 1,m_nodes
       k = (j-1)*n_nodes+i
       acoef(k,k) = one
       rhs(k) = phi(i,j)
    ENDDO
! Bottom Boundary
    j = 1
    DO i = 1,n_nodes
       k = (j-1)*n_nodes+i
       acoef(k,k) = one
       rhs(k) = phi(i,j)
    ENDDO
! Top Boundary
    j = m_nodes
    DO i = 1,n_nodes
       k = (j-1)*n_nodes+i
       acoef(k,k) = one
       rhs(k) = phi(i,j)
    ENDDO
    
! Solve linear system
    CALL gauss_elim(neq,acoef,rhs,phi1D)
      
! Print results
    WRITE(10,*) 'VARIABLES = "X", "Y", "PHIAN", "PHI", "Error"'
    WRITE(10,*) 'ZONE I=', n_nodes, ', J=', m_nodes, ', DATAPACKING=POINT'
    DO j = 1,m_nodes
      DO i = 1,n_nodes
        k = (j-1)*n_nodes+i
        phia(i,j) = 500.0d0*exp(-50.0*((one-x(i))**2+y(j)**2)) + &
                    100.0d0*x(i)*(one-y(j))
        WRITE(10,10) x(i),y(j),phia(i,j),phi1D(k),phia(i,j)-phi1D(k)
      ENDDO
    ENDDO

 10 format(5(1x,F14.6))

    CLOSE(unit=10)

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
