! Start of Main Program
! This is a program to solve the one-dimensional linear 2nd order
! ODE using an unequal mesh and the FDM
! Exercise 2.5

    program main


! Declaration of variables
    IMPLICIT NONE
    INTEGER(4), PARAMETER :: n_nodes = 81
    INTEGER(4) :: iter,n_eq,i
    REAL(8), PARAMETER :: one=1.0D0,two=2.0d0,pi = 3.1415927d0, &
                             & zero= 0.0d0,phi_left=0.0D0,phi_right=1.0D0, &
                             & length=1.0D0,st=1.02d0
    REAL(8), DIMENSION(n_nodes) :: a,b,c,d,s,x,phi,dx
    REAL(8) :: phi_an,dx0,dx2,st_fac,dx12,q_left,q_right
    
! Open output file
    OPEN(unit=10, file="Exercise2_5.out", status = "unknown")

    n_eq = n_nodes
    st_fac = (st - st**n_nodes)/(one-st)
    dx0 = length/st_fac
    DO i = 1,n_nodes-1
      dx(i) = dx0*(st**i)
    ENDDO
    x(1) = zero
    DO i = 2,n_nodes-1
      x(i) = x(i-1)+dx(i-1)
    ENDDO
    x(n_nodes) = length

! Set up coefficient matrix and source vector
    DO i = 2,n_nodes-1
      b(i) = (2.0d0*x(i)-one)*(dx(i)+dx(i-1))*dx(i-1)*dx(i)/2.0d0
      a(i-1) = dx(i)
      c(i) = dx(i-1)
      d(i) = -(dx(i)+dx(i-1))
    ENDDO

! Boundary conditions
    d(1) = one
    c(1) = zero
    b(1) = phi_left
    d(n_nodes) = one
    a(n_nodes-1) = zero
    c(n_nodes) = zero
    b(n_nodes) = phi_right

! Solve system of equations using the Thomas Algorithm
    CALL TRI(n_eq,a,d,c,b,phi)

! Print results
    DO i = 1,n_nodes
      phi_an = (x(i)**3)/3.0d0 - (x(i)**2)/2.0d0 + 7.0d0*x(i)/6.0d0
      WRITE(10,10) x(i),phi(i),phi_an,100.0d0*(phi_an-phi(i))/(phi_an+1.0D-30)
    ENDDO
    dx12 = dx(1)+dx(2)
    q_left = (dx12*dx12*phi(2)-dx(1)*dx(1)*phi(3)- &
             (dx12*dx12-dx(1)*dx(1))*phi(1))/(dx12*dx12*dx(1)-dx12*dx(1)*dx(1))
    dx12 = dx(n_nodes-1)+dx(n_nodes-2)
    q_right = -(dx12*dx12*phi(n_nodes-1)-dx(n_nodes-1)*dx(n_nodes-1)* &
              phi(n_nodes-2)-(dx12*dx12-dx(n_nodes-1)*dx(n_nodes-1))* &
              phi(n_nodes))/(dx12*dx12*dx(n_nodes-1)- &
              dx12*dx(n_nodes-1)*dx(n_nodes-1))
    print *, q_left,q_right,q_left-q_right

    CLOSE(unit=10)

 10 format(4(1x,F14.6))

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
