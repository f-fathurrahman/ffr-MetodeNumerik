! Start of Main Program
! Example 6.1
! Finite Difference Method

    program main

! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 9
    INTEGER(4) :: i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,zero= 0.0d0
    REAL(8), PARAMETER :: phi_left=0.0d0, phi_right=0.0D0
    REAL(8), DIMENSION(n_nodes) :: a,d,c,b,s,phi,x,phian
    REAL(8) :: dx,dx2,length,J_l,J_r,c1,c2,sphi,imb
    REAL(8) :: J_an_l,J_an_r,sphi_an,imb_an
    
! Open output file
    OPEN(unit=10, file="Example6_1_fdm.out", status = "unknown")

    length = one
    dx = length/(n_nodes-1)
    dx2 = dx*dx
    DO i = 1,n_nodes
      x(i) = (i-1)*dx
      s(i) = cos(x(i))
    ENDDO

! Set up coefficient matrix and source vector

    a(:) = zero; b(:) = zero; c(:) = zero; d(:) = zero; phi(:) = zero

! Boundary conditions

    d(1) = one
    b(1) = phi_left
    DO i = 2,n_nodes-1
      b(i) = -s(i)
      d(i) = two/dx2
      c(i) = -one/dx2
      a(i-1) = -one/dx2
    ENDDO
    d(n_nodes) = one
    b(n_nodes) = phi_right

! Solve system of equations using the TDMA

    CALL TRI(n_nodes,a,d,c,b,phi)

! Print results
! Phi
    c2 = one
    c1 = (-c2+cos(length))/length
    DO i = 1,n_nodes
      phian(i) = -cos(x(i)) + c1*x(i) + c2
      WRITE(*,10) x(i),phian(i),phi(i),phian(i)-phi(i)
      WRITE(10,10) x(i),phian(i),phi(i),phian(i)-phi(i)
    ENDDO
! Fluxes
    J_l = -(4.0d0*phi(2)-phi(3)-3.0d0*phi(1))/(two*dx)
    J_r = (4.0d0*phi(n_nodes-1)-phi(n_nodes-2)-3.0d0*phi(n_nodes))/(two*dx)
    sphi = s(1)*dx/two
    DO i = 2,n_nodes-1
      sphi = sphi + s(i)*dx
    ENDDO
    sphi = sphi + s(n_nodes)*dx/two
    sphi = -sphi

    imb = J_l-J_r+sphi  ! Net in

    WRITE(*,*) "  Left          ","Right         ","Source         ","Imbalance   "
    WRITE(10,*)"  Left          ","Right         ","Source         ","Imbalance   "
    WRITE(*,10) J_l,J_r,sphi,imb
    WRITE(10,10) J_l,J_r,sphi,imb

    J_an_l = -(c1)
    J_an_r = -(sin(length)+c1)
    sphi_an = -sin(length)
    imb_an = J_an_l-J_an_r+sphi_an
    WRITE(*,*) "Analytical Solution"
    WRITE(*,10) J_an_l, J_an_r, sphi_an,imb_an
    WRITE(10,*) "Analytical Solution"
    WRITE(10,10) J_an_l, J_an_r, sphi_an,imb_an

 10 format(4(F14.6))

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


