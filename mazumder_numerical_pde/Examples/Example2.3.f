! Start of Main Program
! Example 2.3

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 9
    INTEGER(4) :: iter,n_eq,i,approach
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0, &
                             & zero= 0.0d0,phi_right=1.0D0, &
                             & length=1.0D0,phi_left = 0.0d0
    REAL(8), DIMENSION(n_nodes) :: b,x,phi,phi_an,phi1,phi2,s
    REAL(8), DIMENSION(n_nodes,n_nodes) :: a
    REAL(8) :: dx,dx2,err_1,err_2,c1,c2
    
! Open output file
    OPEN(unit=10, file="Example2_3.out", status = "unknown")

    dx = length/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
      s(i) = exp(x(i))
    ENDDO

    approach = 1

 11 CONTINUE

    a(:,:) = zero; b(:) = zero

! Set up coefficient matrix and source vector

! Boundary conditions

    IF(approach == 1)THEN  ! Second order

      a(1,1) = one
      b(1) = phi_left
      DO i = 2,n_nodes-1
        b(i) = s(i)
        a(i,i-1) = one/dx2
        a(i,i+1) = one/dx2
        a(i,i) = -two/dx2
      ENDDO
      a(n_nodes,n_nodes) = one
      b(n_nodes) = phi_right

    ELSE  ! Fourth order

      a(1,1) = one
      b(1) = phi_left
      a(2,2) = -two/dx2
      a(2,1) = one/dx2
      a(2,3) = one/dx2
      b(2) = s(2)
      DO i = 3,n_nodes-2
        b(i) = s(i)
        a(i,i-1) = one/dx2 + one/(3.0d0*dx2)
        a(i,i+1) = one/dx2 + one/(3.0d0*dx2)
        a(i,i) = -two/dx2 - 0.5d0/dx2
        a(i,i+2) = -one/(12.0d0*dx2)
        a(i,i-2) = -one/(12.0d0*dx2)
      ENDDO
      a(n_nodes-1,n_nodes-1) = -two/dx2
      a(n_nodes-1,n_nodes-2) = one/dx2
      a(n_nodes-1,n_nodes) = one/dx2
      b(n_nodes-1) = s(n_nodes-1)
      a(n_nodes,n_nodes) = one
      b(n_nodes) = phi_right

    ENDIF


! Solve system of equations using the Gaussian Elimination

    CALL gauss_elim(n_nodes,a,b,phi)

    IF(approach == 1)THEN
      phi1(:) = phi(:)
      approach = 2
      GO TO 11
    ELSE
      phi2(:) = phi(:)
    ENDIF

! Print results
    DO i = 1,n_nodes
      c2 = phi_left - one
      c1 = phi_right - c2 - exp(one)
      phi_an(i) = exp(x(i)) + c1*x(i) + c2
      err_1 = 100.0d0*(phi_an(i)-phi1(i))/phi_right
      err_2 = 100.0d0*(phi_an(i)-phi2(i))/phi_right
      WRITE(*,10) x(i),phi_an(i),phi1(i),err_1,phi2(i),err_2
      WRITE(10,10) x(i),phi_an(i),phi1(i),err_1,phi2(i),err_2
    ENDDO

 10 format(6(F13.6))

    CLOSE(unit=10)

    END Program main

!--------------------------------------------------------------------------

! Tridiagonal solver (uses Thomas Algorithm)

    SUBROUTINE TRI(n,a,d,c,b,x)

    INTEGER(4), INTENT(in) :: n
    INTEGER(4) :: i
    REAL(8) :: xmult
    REAL(8), DIMENSION(n) :: a,b,c,d,x

    DO i = 2,n
      xmult = a(i-1)/d(i-1)
      d(i) = d(i) - xmult*c(i-1)
      b(i) = b(i) - xmult*b(i-1)
    ENDDO

    x(n) = b(n)/d(n)
    DO i = n-1,1,-1
      x(i) = (b(i)-c(i)*x(i+1))/d(i)
    ENDDO

    END SUBROUTINE TRI

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

